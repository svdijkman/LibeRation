struct ObjectiveTape {
  CppAD::ADFun<double> fun;
  std::vector<std::string> domain_names;
  std::vector<std::string> dynamic_columns;
  std::vector<int> dynamic_observed_rows;
  std::vector<int> structural_dvid;
  std::vector<double> dynamic_values;
  int n_rows = 0;
  libertad::SparseHessianCache hessian_cache;
};

class TapePathChange : public std::runtime_error {
 public:
  explicit TapePathChange(const std::string& context,
                          std::vector<double> point = std::vector<double>())
      : std::runtime_error("CppAD tape path changed in " + context +
                           "; automatic retaping is required."),
        point_(std::move(point)) {}

  const std::vector<double>& point() const { return point_; }

 private:
  std::vector<double> point_;
};

void require_unchanged_path(CppAD::ADFun<double>& fun,
                            const std::string& context) {
  if (fun.compare_change_number() != 0U) throw TapePathChange(context);
}
struct PredictionTape {
  CppAD::ADFun<double> fun;
  std::vector<std::string> domain_names;
  std::vector<std::string> dynamic_columns;
  std::vector<double> dynamic_values;
  int n_rows = 0;
  std::string propagation_kernel;
  std::size_t operation_count = 0;
  std::size_t variable_count = 0;
  std::string derivative_strategy = "not-evaluated";
  std::size_t jacobian_nonzeros = 0;
};

bool structural_data_input(const std::string& name) {
  static const std::unordered_map<std::string, bool> structural = {
    {"ID", true}, {"TIME", true}, {"AMT", true}, {"RATE", true},
    {"II", true}, {"ADDL", true}, {"EVID", true}, {"CMT", true},
    {"SS", true}, {"MIXNUM", true}, {"DVID", true}, {"DV", true},
    {"MDV", true}, {"LLOQ", true}, {"BLQ", true}, {"CENS", true},
    {".ID_INDEX", true}, {".OCC_INDEX", true}
  };
  return structural.find(name) != structural.end();
}

bool data_backed_model_input(const std::string& name) {
  if (structural_data_input(name) || name == "F" || name == "T" ||
      name == "MIXNUM" || starts_with(name, "THETA_") ||
      starts_with(name, "ETA_") || starts_with(name, "SIGMA_") ||
      starts_with(name, "ERR_") || starts_with(name, "A_")) {
    return false;
  }
  return true;
}

std::vector<std::string> prediction_dynamic_columns(
    const ModelEngine& engine, const Rcpp::DataFrame& data) {
  std::vector<std::string> columns;
  std::unordered_map<std::string, bool> seen;
  auto append = [&](const std::vector<std::string>& inputs) {
    for (const std::string& name : inputs) {
      if (!data_backed_model_input(name) ||
          !data.containsElementNamed(name.c_str()) || seen[name]) continue;
      for (int row = 0; row < data.nrows(); ++row) {
        if (!std::isfinite(data_value(data, name, row))) {
          throw std::domain_error("Dynamic model input '" + name +
                                  "' contains a non-finite value.");
        }
      }
      seen[name] = true;
      columns.push_back(name);
    }
  };
  append(engine.pred->input_names);
  if (engine.des) append(engine.des->input_names);
  if (engine.post_pred) append(engine.post_pred->input_names);
  return columns;
}

std::vector<double> prediction_dynamic_values(
    const std::vector<std::string>& columns, const Rcpp::DataFrame& data,
    int expected_rows = -1) {
  if (expected_rows >= 0 && data.nrows() != expected_rows) {
    throw std::invalid_argument("Dynamic prediction data has a different row count.");
  }
  std::vector<double> values;
  values.reserve(columns.size() * static_cast<std::size_t>(data.nrows()));
  for (const std::string& name : columns) {
    if (!data.containsElementNamed(name.c_str())) {
      throw std::invalid_argument("Dynamic prediction data is missing column '" + name + "'.");
    }
    for (int row = 0; row < data.nrows(); ++row) {
      const double value = data_value(data, name, row);
      if (!std::isfinite(value)) {
        throw std::domain_error("Dynamic prediction input '" + name +
                                "' contains a non-finite value.");
      }
      values.push_back(value);
    }
  }
  return values;
}

std::vector<int> fo_observed_rows(const Rcpp::DataFrame& data) {
  Rcpp::NumericVector dv = data["DV"];
  Rcpp::NumericVector evid = data["EVID"];
  Rcpp::NumericVector mdv = data["MDV"];
  std::vector<int> observed;
  for (int row = 0; row < data.nrows(); ++row) {
    if (evid[row] == 0.0 && mdv[row] == 0.0 && std::isfinite(dv[row])) {
      observed.push_back(row);
    }
  }
  return observed;
}

std::vector<int> fo_dvid_values(const Rcpp::DataFrame& data) {
  std::vector<int> result(static_cast<std::size_t>(data.nrows()), 1);
  if (!data.containsElementNamed("DVID")) return result;
  Rcpp::NumericVector dvid = data["DVID"];
  for (int row = 0; row < data.nrows(); ++row) {
    result[static_cast<std::size_t>(row)] =
      std::max(1, static_cast<int>(dvid[row]));
  }
  return result;
}

std::vector<double> fo_dynamic_values(const ObjectiveTape& tape,
                                      const Rcpp::DataFrame& data) {
  if (tape.n_rows != data.nrows()) {
    throw std::invalid_argument("A shared FO tape received a different number of rows.");
  }
  if (fo_observed_rows(data) != tape.dynamic_observed_rows) {
    throw std::invalid_argument("A shared FO tape received a different observation pattern.");
  }
  if (fo_dvid_values(data) != tape.structural_dvid) {
    throw std::invalid_argument("A shared FO tape received a different DVID pattern.");
  }
  std::vector<double> values = prediction_dynamic_values(
    tape.dynamic_columns, data, tape.n_rows);
  Rcpp::NumericVector dv = data["DV"];
  values.reserve(values.size() + tape.dynamic_observed_rows.size());
  for (int row : tape.dynamic_observed_rows) {
    const double value = dv[row];
    if (!std::isfinite(value)) {
      throw std::domain_error("A shared FO tape received a non-finite observation.");
    }
    values.push_back(value);
  }
  return values;
}

void set_fo_dynamic(ObjectiveTape& tape, const Rcpp::DataFrame& data) {
  const std::vector<double> values = fo_dynamic_values(tape, data);
  if (values.size() != tape.fun.size_dyn_ind()) {
    throw std::logic_error("Shared FO dynamic data do not match the recorded tape.");
  }
  if (!values.empty()) tape.fun.new_dynamic(values);
  tape.dynamic_values = values;
}

std::vector<double> flatten_parameters(const Rcpp::NumericVector& theta,
                                       const Rcpp::NumericMatrix& eta,
                                       const Rcpp::NumericVector& sigma) {
  std::vector<double> result;
  result.reserve(theta.size() + eta.size() + sigma.size());
  for (double value : theta) result.push_back(value);
  for (int row = 0; row < eta.nrow(); ++row) {
    for (int column = 0; column < eta.ncol(); ++column) result.push_back(eta(row, column));
  }
  for (double value : sigma) result.push_back(value);
  return result;
}

std::vector<std::string> parameter_names(int n_theta, int n_subjects,
                                         int n_eta, int n_sigma) {
  std::vector<std::string> names;
  for (int i = 0; i < n_theta; ++i) names.push_back("THETA_" + std::to_string(i + 1));
  for (int subject = 0; subject < n_subjects; ++subject) {
    for (int i = 0; i < n_eta; ++i) {
      names.push_back("ETA_" + std::to_string(subject + 1) + "_" + std::to_string(i + 1));
    }
  }
  for (int i = 0; i < n_sigma; ++i) names.push_back("SIGMA_" + std::to_string(i + 1));
  return names;
}

std::unique_ptr<PredictionTape> record_prediction_tape(
    const ModelEngine& engine, const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta, const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma) {
  const int minimum_eta_columns = required_eta_columns(engine, data);
  const int between_eta = engine.n_eta - engine.iov;
  if (theta.size() != engine.n_theta || eta.ncol() < minimum_eta_columns ||
      (!engine.re_enabled && engine.iov > 0 &&
       (eta.ncol() - between_eta) % engine.iov != 0)) {
    throw std::invalid_argument("Prediction tape parameter dimensions are inconsistent with the model.");
  }
  std::vector<double> point = flatten_parameters(theta, eta, sigma);
  std::vector<CppAD::AD<double>> independent(point.begin(), point.end());
  const std::vector<std::string> dynamic_columns =
    prediction_dynamic_columns(engine, data);
  const std::vector<double> dynamic_values =
    prediction_dynamic_values(dynamic_columns, data);
  std::vector<CppAD::AD<double>> dynamic(dynamic_values.begin(), dynamic_values.end());
  if (dynamic.empty()) CppAD::Independent(independent);
  else CppAD::Independent(independent, dynamic);
  DynamicDataT<CppAD::AD<double>> dynamic_data;
  dynamic_data.n_rows = data.nrows();
  dynamic_data.values = dynamic;
  for (std::size_t column = 0; column < dynamic_columns.size(); ++column) {
    dynamic_data.column_positions[dynamic_columns[column]] = column;
  }
  std::size_t cursor = 0;
  std::vector<CppAD::AD<double>> theta_ad(static_cast<std::size_t>(theta.size()));
  for (auto& value : theta_ad) value = independent[cursor++];
  std::vector<CppAD::AD<double>> eta_ad(static_cast<std::size_t>(eta.size()));
  for (auto& value : eta_ad) value = independent[cursor++];
  std::vector<CppAD::AD<double>> sigma_ad(static_cast<std::size_t>(sigma.size()));
  for (auto& value : sigma_ad) value = independent[cursor++];
  std::vector<CppAD::AD<double>> predictions = simulate_analytical_t(
    engine, data, theta_ad, eta_ad, sigma_ad, std::vector<int>(),
    dynamic.empty() ? nullptr : &dynamic_data);
  auto tape = std::make_unique<PredictionTape>();
  tape->fun.Dependent(independent, predictions);
  tape->fun.optimize();
  tape->operation_count = tape->fun.size_op();
  tape->variable_count = tape->fun.size_var();
  tape->domain_names = parameter_names(theta.size(), eta.nrow(), eta.ncol(), sigma.size());
  tape->dynamic_columns = dynamic_columns;
  tape->dynamic_values = dynamic_values;
  tape->n_rows = data.nrows();
  tape->propagation_kernel = propagation_kernel_name(engine);
  return tape;
}
std::vector<double> prediction_point(PredictionTape& tape,
                                     const Rcpp::NumericVector& point) {
  if (point.size() != static_cast<R_xlen_t>(tape.domain_names.size())) {
    throw std::invalid_argument("Prediction tape point has the wrong length.");
  }
  return Rcpp::as<std::vector<double>>(point);
}

template <class Scalar>
Scalar scalar_floor_t(const Scalar& value, double floor) {
  return libertad::choose_gt(value, Scalar(floor), value, Scalar(floor));
}

inline double scalar_erf_t(double value) { return std::erf(value); }
inline CppAD::AD<double> scalar_erf_t(const CppAD::AD<double>& value) {
  return CppAD::erf(value);
}

template <class Scalar>
Scalar normal_cdf_t(const Scalar& value) {
  const Scalar probability = Scalar(0.5) *
    (Scalar(1.0) + scalar_erf_t(value * Scalar(0.7071067811865475244)));
  return scalar_floor_t(probability, 1e-300);
}

template <class Scalar>
Scalar residual_variance_t(const ModelEngine& engine,
                           const Scalar& prediction,
                           const std::vector<Scalar>& sigma,
                           int dvid) {
  int per_response = engine.error_type == "combined" ? 2 :
    (engine.error_type == "power" ? 2 : 1);
  int offset = std::max(0, dvid - 1) * per_response;
  if (offset + per_response > static_cast<int>(sigma.size())) offset = 0;
  if (sigma.empty()) {
    throw std::domain_error("Likelihood evaluation requires residual SIGMA parameters.");
  }
  const Scalar s1 = sigma[static_cast<std::size_t>(offset)];
  const auto sigma_variance = [&](const Scalar& value) {
    return engine.sigma_parameterization == "variance" ? value : value * value;
  };
  Scalar variance;
  if (engine.error_type == "additive" || engine.error_type == "exponential") {
    variance = sigma_variance(s1);
  } else if (engine.error_type == "proportional") {
    variance = sigma_variance(s1) * prediction * prediction;
  } else if (engine.error_type == "power") {
    if (sigma.size() <= static_cast<std::size_t>(offset + 1)) {
      throw std::domain_error("Power residual error requires two SIGMA parameters.");
    }
    const Scalar magnitude = scalar_floor_t(libertad::scalar_abs(prediction), 1e-12);
    variance = sigma_variance(s1) * libertad::scalar_pow(
      magnitude, Scalar(2.0) * sigma[static_cast<std::size_t>(offset + 1)]);
  } else {
    if (sigma.size() <= static_cast<std::size_t>(offset + 1)) {
      throw std::domain_error("Combined residual error requires two SIGMA parameters.");
    }
    const Scalar s2 = sigma[static_cast<std::size_t>(offset + 1)];
    variance = sigma_variance(s1) * prediction * prediction + sigma_variance(s2);
  }
  return scalar_floor_t(variance, 1e-16);
}

template <class Scalar>
MatrixT<Scalar> omega_matrix_t(const ModelEngine& engine,
                               const std::vector<Scalar>& omega) {
  if (omega.size() != engine.omega_rows.size()) {
    throw std::invalid_argument("OMEGA parameter vector has the wrong length.");
  }
  MatrixT<Scalar> covariance = MatrixT<Scalar>::Zero(engine.n_eta, engine.n_eta);
  for (std::size_t i = 0; i < omega.size(); ++i) {
    const int row = engine.omega_rows[i];
    const int column = engine.omega_cols[i];
    if (row < 0 || column < 0 || row >= engine.n_eta || column >= engine.n_eta) {
      throw std::logic_error("OMEGA covariance index is outside the ETA dimension.");
    }
    covariance(row, column) = omega[i];
    covariance(column, row) = omega[i];
  }
  return covariance;
}

template <class Scalar>
MatrixT<Scalar> expanded_omega_t(const ModelEngine& engine,
                                 const Rcpp::DataFrame& data,
                                 const MatrixT<Scalar>& base,
                                 int expanded_dimension) {
  if (engine.re_enabled) {
    MatrixT<Scalar> output = MatrixT<Scalar>::Zero(
      expanded_dimension, expanded_dimension);
    int offset = 0;
    for (std::size_t block_index = 0; block_index < engine.re_blocks.size(); ++block_index) {
      const std::vector<int>& indices = engine.re_blocks[block_index];
      const std::string total_name = ".RE_TOTAL_" + std::to_string(block_index + 1U);
      const int units = static_cast<int>(data_value(data, total_name, 0));
      for (int unit = 0; unit < units; ++unit) {
        for (std::size_t row = 0; row < indices.size(); ++row) {
          for (std::size_t column = 0; column < indices.size(); ++column) {
            output(offset + unit * static_cast<int>(indices.size()) + static_cast<int>(row),
                   offset + unit * static_cast<int>(indices.size()) + static_cast<int>(column)) =
              base(indices[row], indices[column]);
          }
        }
      }
      offset += units * static_cast<int>(indices.size());
    }
    if (offset != expanded_dimension) {
      throw std::invalid_argument("Expanded random-effect covariance has the wrong dimension.");
    }
    return output;
  }
  if (engine.iov == 0) return base;
  const int between = engine.n_eta - engine.iov;
  if (expanded_dimension < between ||
      (expanded_dimension - between) % engine.iov != 0) {
    throw std::invalid_argument("Expanded IOV covariance has an invalid layout.");
  }
  const int occasions = (expanded_dimension - between) / engine.iov;
  MatrixT<Scalar> output = MatrixT<Scalar>::Zero(
    expanded_dimension, expanded_dimension);
  if (between) output.topLeftCorner(between, between) =
    base.topLeftCorner(between, between);
  for (int occasion = 0; occasion < occasions; ++occasion) {
    const int target = between + occasion * engine.iov;
    output.block(target, target, engine.iov, engine.iov) =
      base.bottomRightCorner(engine.iov, engine.iov);
  }
  return output;
}

template <class Scalar>
Scalar omega_subject_prior_t(const MatrixT<Scalar>& covariance,
                             const VectorT<Scalar>& eta) {
  const Eigen::Index n = covariance.rows();
  MatrixT<Scalar> lower = MatrixT<Scalar>::Zero(n, n);
  Scalar logdet = Scalar(0.0);
  for (Eigen::Index row = 0; row < n; ++row) {
    for (Eigen::Index column = 0; column <= row; ++column) {
      Scalar value = covariance(row, column);
      for (Eigen::Index k = 0; k < column; ++k) {
        value -= lower(row, k) * lower(column, k);
      }
      if (row == column) {
        if (!(scalar_value(value) > 1e-14)) {
          throw std::domain_error("OMEGA matrix is not positive definite at the recording point.");
        }
        lower(row, column) = libertad::scalar_sqrt(value);
        logdet += Scalar(2.0) * libertad::scalar_log(lower(row, column));
      } else {
        lower(row, column) = value / lower(column, column);
      }
    }
  }
  VectorT<Scalar> standardized(n);
  for (Eigen::Index row = 0; row < n; ++row) {
    Scalar value = eta[row];
    for (Eigen::Index column = 0; column < row; ++column) {
      value -= lower(row, column) * standardized[column];
    }
    standardized[row] = value / lower(row, row);
  }
  Scalar quadratic = Scalar(0.0);
  for (Eigen::Index i = 0; i < n; ++i) quadratic += standardized[i] * standardized[i];
  return logdet + quadratic;
}

template <class Scalar>
std::vector<Scalar> evaluate_error_outputs_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data, int row, int subject,
    const std::vector<Scalar>& theta, const std::vector<Scalar>& eta,
    int eta_columns, const std::vector<Scalar>& sigma, int mixture_number,
    const Scalar& prediction, bool first, double previous_dv,
    double previous_time, const std::vector<std::size_t>& selected,
    const char* context,
    const ParametersT<Scalar>* input_overrides = nullptr) {
  if (!engine.error || selected.empty()) {
    throw std::logic_error(std::string("Compiled ") + context + " outputs are missing.");
  }
  const double current_time = row_optional(data, "TIME", row, 0.0);
  const double current_dv = row_optional(data, "DV", row, NA_REAL);
  ParametersT<Scalar> parameters = evaluate_parameters_t(
    engine, data, row, subject, theta, eta, eta_columns, sigma,
    mixture_number);
  std::vector<Scalar> inputs(engine.error->input_names.size(), Scalar(0.0));
  for (std::size_t index = 0; index < engine.error->input_names.size(); ++index) {
    const std::string& name = engine.error->input_names[index];
    int parameter_index = indexed_name(name, "THETA_");
    if (parameter_index >= 0) {
      if (parameter_index >= static_cast<int>(theta.size())) {
        throw std::out_of_range("THETA index exceeds values in user likelihood.");
      }
      inputs[index] = theta[static_cast<std::size_t>(parameter_index)];
      continue;
    }
    parameter_index = indexed_name(name, "ETA_");
    if (parameter_index >= 0) {
      const int column = eta_column(
        engine, data, row, parameter_index, eta_columns);
      const std::size_t position =
        static_cast<std::size_t>(subject * eta_columns + column);
      if (position >= eta.size()) {
        throw std::out_of_range("ETA index exceeds values in user likelihood.");
      }
      inputs[index] = eta[position];
      continue;
    }
    parameter_index = indexed_name(name, "SIGMA_");
    if (parameter_index >= 0) {
      if (parameter_index >= static_cast<int>(sigma.size())) {
        throw std::out_of_range("SIGMA index exceeds values in user likelihood.");
      }
      inputs[index] = sigma[static_cast<std::size_t>(parameter_index)];
      continue;
    }
    if (name == "F" || name == "PRED" || name == "IPRED") {
      inputs[index] = prediction;
      continue;
    }
    if (name == "DV") {
      inputs[index] = Scalar(current_dv);
      continue;
    }
    if (name == "PREV_DV") {
      inputs[index] = Scalar(previous_dv);
      continue;
    }
    if (name == "PREV_TIME") {
      inputs[index] = Scalar(previous_time);
      continue;
    }
    if (name == "DT") {
      inputs[index] = Scalar(first ? 0.0 : current_time - previous_time);
      continue;
    }
    if (name == "FIRST") {
      inputs[index] = Scalar(first ? 1.0 : 0.0);
      continue;
    }
    if (name == "MIXNUM") {
      inputs[index] = Scalar(mixture_number);
      continue;
    }
    if (input_overrides != nullptr) {
      const auto overridden = input_overrides->find(name);
      if (overridden != input_overrides->end()) {
        inputs[index] = overridden->second;
        continue;
      }
    }
    const auto assigned = parameters.find(name);
    if (assigned != parameters.end()) {
      inputs[index] = assigned->second;
      continue;
    }
    inputs[index] = Scalar(data_value(data, name, row));
    if (!std::isfinite(scalar_value(inputs[index]))) {
      throw std::domain_error("User likelihood input '" + name +
                              "' is non-finite at row " +
                              std::to_string(row + 1) + ".");
    }
  }
  const std::vector<Scalar> output = engine.error->eval_outputs(inputs, selected);
  for (const Scalar& value : output) {
    if (!std::isfinite(scalar_value(value))) {
      throw std::domain_error(std::string(context) +
                              " returned a non-finite value at row " +
                              std::to_string(row + 1) + ".");
    }
  }
  return output;
}

template <class Scalar>
Scalar user_likelihood_nll_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data, int row, int subject,
    const std::vector<Scalar>& theta, const std::vector<Scalar>& eta,
    int eta_columns, const std::vector<Scalar>& sigma, int mixture_number,
    const Scalar& prediction, bool first, double previous_dv,
    double previous_time) {
  if (engine.likelihood_output.size() != 1U) {
    throw std::logic_error("Compiled user likelihood is missing.");
  }
  const std::vector<Scalar> output = evaluate_error_outputs_t(
    engine, data, row, subject, theta, eta, eta_columns, sigma,
    mixture_number, prediction, first, previous_dv, previous_time,
    engine.likelihood_output, "User likelihood");
  if (engine.likelihood_scale == "log") return Scalar(-2.0) * output[0];
  if (!(scalar_value(output[0]) > 0.0)) {
    throw std::domain_error("LIK must be positive at row " +
                            std::to_string(row + 1) + ".");
  }
  return Scalar(-2.0) * libertad::scalar_log(
    scalar_floor_t(output[0], 1e-300));
}

template <class Scalar>
Scalar log_sum_exp_t(const std::vector<Scalar>& values) {
  if (values.empty()) throw std::invalid_argument("log-sum-exp requires values.");
  Scalar maximum = values.front();
  for (std::size_t index = 1; index < values.size(); ++index) {
    maximum = libertad::choose_gt(values[index], maximum, values[index], maximum);
  }
  Scalar total = Scalar(0.0);
  for (const Scalar& value : values) total += libertad::scalar_exp(value - maximum);
  return maximum + libertad::scalar_log(total);
}

template <class Scalar>
std::vector<Scalar> hmm_normalize_weights_t(
    const std::vector<Scalar>& values, const std::string& scale,
    const char* label, int row) {
  if (values.empty()) throw std::invalid_argument("Hidden Markov weights are empty.");
  std::vector<Scalar> result(values.size());
  if (scale == "log") {
    const Scalar normalizer = log_sum_exp_t(values);
    for (std::size_t index = 0; index < values.size(); ++index) {
      result[index] = libertad::scalar_exp(values[index] - normalizer);
    }
    return result;
  }
  Scalar total = Scalar(0.0);
  for (const Scalar& value : values) {
    if (scalar_value(value) < 0.0) {
      throw std::domain_error(std::string("HMM ") + label +
                              " weights must be non-negative at row " +
                              std::to_string(row + 1) + ".");
    }
    total += value;
  }
  if (!(scalar_value(total) > 0.0)) {
    throw std::domain_error(std::string("HMM ") + label +
                            " weights sum to zero at row " +
                            std::to_string(row + 1) + ".");
  }
  for (std::size_t index = 0; index < values.size(); ++index) {
    result[index] = values[index] / total;
  }
  return result;
}

template <class Scalar>
struct HmmRowComponents {
  std::vector<Scalar> initial;
  std::vector<std::vector<Scalar>> transition;
  std::vector<Scalar> emission;
  std::vector<Scalar> log_emission;
};

template <class Scalar>
HmmRowComponents<Scalar> hmm_row_components_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data, int row, int subject,
    const std::vector<Scalar>& theta, const std::vector<Scalar>& eta,
    int eta_columns, const std::vector<Scalar>& sigma, int mixture_number,
    const Scalar& prediction, bool first, double previous_dv,
    double previous_time) {
  const int states = engine.hmm_states;
  const std::size_t state_count = static_cast<std::size_t>(states);
  const std::size_t transition_count = engine.hmm_continuous ?
    state_count * (state_count - 1U) : state_count * state_count;
  const std::size_t expected = state_count + transition_count + state_count;
  if (!engine.hmm_enabled || states < 2 || engine.hmm_outputs.size() != expected) {
    throw std::logic_error("Compiled hidden Markov likelihood is inconsistent.");
  }
  const std::vector<Scalar> output = evaluate_error_outputs_t(
    engine, data, row, subject, theta, eta, eta_columns, sigma,
    mixture_number, prediction, first, previous_dv, previous_time,
    engine.hmm_outputs, "Hidden Markov likelihood");
  HmmRowComponents<Scalar> result;
  if (first) {
    result.initial = hmm_normalize_weights_t(
      std::vector<Scalar>(output.begin(), output.begin() + states),
      engine.hmm_initial_scale, "initial", row);
  } else {
    const std::size_t transition_offset = state_count;
    result.transition.resize(state_count);
    if (engine.hmm_continuous) {
      MatrixT<Scalar> generator = MatrixT<Scalar>::Zero(states, states);
      std::size_t rate_index = transition_offset;
      for (int from = 0; from < states; ++from) {
        for (int to = 0; to < states; ++to) {
          if (from == to) continue;
          Scalar rate = output[rate_index++];
          if (engine.hmm_rate_scale == "log") rate = libertad::scalar_exp(rate);
          if (scalar_value(rate) < 0.0) {
            throw std::domain_error(
              "Continuous-time HMM rates must be non-negative at row " +
              std::to_string(row + 1) + ".");
          }
          generator(from, to) = rate;
          generator(from, from) -= rate;
        }
      }
      const double current_time = row_optional(data, "TIME", row, 0.0);
      const double dt = current_time - previous_time;
      if (!std::isfinite(dt) || dt < 0.0) {
        throw std::domain_error(
          "Continuous-time HMM observations must be ordered by non-decreasing TIME within sequence.");
      }
      const MatrixT<Scalar> transition = matrix_exp_pade(
        MatrixT<Scalar>(generator * Scalar(dt)));
      for (int from = 0; from < states; ++from) {
        result.transition[static_cast<std::size_t>(from)].resize(state_count);
        for (int to = 0; to < states; ++to) {
          result.transition[static_cast<std::size_t>(from)]
                           [static_cast<std::size_t>(to)] = transition(from, to);
        }
      }
    } else {
      for (int from = 0; from < states; ++from) {
        const auto begin = output.begin() + static_cast<std::ptrdiff_t>(
          transition_offset + static_cast<std::size_t>(from * states));
        result.transition[static_cast<std::size_t>(from)] =
          hmm_normalize_weights_t(
            std::vector<Scalar>(begin, begin + states),
            engine.hmm_transition_scale, "transition", row);
      }
    }
  }
  const std::size_t emission_offset = state_count + transition_count;
  result.emission.resize(state_count);
  result.log_emission.resize(state_count);
  for (std::size_t state = 0; state < state_count; ++state) {
    const Scalar value = output[emission_offset + state];
    if (engine.hmm_emission_scale == "log") {
      result.log_emission[state] = value;
      result.emission[state] = Scalar(0.0);
    } else {
      if (scalar_value(value) < 0.0) {
        throw std::domain_error(
          "HMM emission likelihoods must be non-negative at row " +
          std::to_string(row + 1) + ".");
      }
      result.emission[state] = value;
      result.log_emission[state] = libertad::scalar_log(
        scalar_floor_t(value, 1e-300));
    }
  }
  return result;
}

template <class Scalar>
Scalar hmm_row_nll_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data, int row, int subject,
    const std::vector<Scalar>& theta, const std::vector<Scalar>& eta,
    int eta_columns, const std::vector<Scalar>& sigma, int mixture_number,
    const Scalar& prediction, bool first, double previous_dv,
    double previous_time, std::vector<Scalar>& filtered,
    HmmRowComponents<Scalar>* components_output = nullptr) {
  const int states = engine.hmm_states;
  const std::size_t state_count = static_cast<std::size_t>(states);
  HmmRowComponents<Scalar> components = hmm_row_components_t(
    engine, data, row, subject, theta, eta, eta_columns, sigma,
    mixture_number, prediction, first, previous_dv, previous_time);
  std::vector<Scalar> prior(state_count);
  if (first) {
    prior = components.initial;
  } else {
    if (filtered.size() != state_count) {
      throw std::logic_error("Hidden Markov filter state has the wrong dimension.");
    }
    std::fill(prior.begin(), prior.end(), Scalar(0.0));
    for (int from = 0; from < states; ++from) {
      for (int to = 0; to < states; ++to) {
        prior[static_cast<std::size_t>(to)] +=
          filtered[static_cast<std::size_t>(from)] *
          components.transition[static_cast<std::size_t>(from)]
                               [static_cast<std::size_t>(to)];
      }
    }
  }
  if (components_output != nullptr) *components_output = components;
  filtered.assign(state_count, Scalar(0.0));
  if (engine.hmm_emission_scale == "log") {
    std::vector<Scalar> log_weight(state_count);
    for (std::size_t state = 0; state < state_count; ++state) {
      log_weight[state] = libertad::scalar_log(
        scalar_floor_t(prior[state], 1e-300)) + components.log_emission[state];
    }
    const Scalar log_likelihood = log_sum_exp_t(log_weight);
    for (std::size_t state = 0; state < state_count; ++state) {
      filtered[state] = libertad::scalar_exp(log_weight[state] - log_likelihood);
    }
    return Scalar(-2.0) * log_likelihood;
  }
  Scalar likelihood = Scalar(0.0);
  for (std::size_t state = 0; state < state_count; ++state) {
    filtered[state] = prior[state] * components.emission[state];
    likelihood += filtered[state];
  }
  if (!(scalar_value(likelihood) > 0.0)) {
    throw std::domain_error("HMM observation likelihood is zero at row " +
                            std::to_string(row + 1) + ".");
  }
  for (Scalar& value : filtered) value /= likelihood;
  return Scalar(-2.0) * libertad::scalar_log(
    scalar_floor_t(likelihood, 1e-300));
}

template <class Scalar>
Scalar positive_definite_gaussian_nll_t(
    const MatrixT<Scalar>& covariance, const VectorT<Scalar>& residual,
    const std::string& context);

inline int residual_group_endpoint(const ResidualGroupSpec& group, double dvid) {
  for (std::size_t index = 0; index < group.dvid.size(); ++index) {
    if (group.dvid[index] == dvid) return static_cast<int>(index);
  }
  return -1;
}

inline int residual_group_for_dvid(const ModelEngine& engine, double dvid) {
  for (std::size_t index = 0; index < engine.residual_groups.size(); ++index) {
    if (residual_group_endpoint(engine.residual_groups[index], dvid) >= 0) {
      return static_cast<int>(index);
    }
  }
  return -1;
}

template <class Scalar>
Scalar residual_group_correlation_t(
    const ResidualGroupSpec& group, int row, int column,
    const std::vector<Scalar>& theta, const std::vector<Scalar>& sigma) {
  const int dimension = static_cast<int>(group.dvid.size());
  if (row < 0 || column < 0 || row >= dimension || column >= dimension) {
    throw std::out_of_range("Residual-group endpoint index is outside its correlation matrix.");
  }
  const std::size_t position = static_cast<std::size_t>(row * dimension + column);
  const std::string& source = group.source[position];
  if (source == "fixed") return Scalar(group.value[position]);
  const std::vector<Scalar>& parameters = source == "theta" ? theta : sigma;
  const int index = group.index[position];
  if (index < 0 || index >= static_cast<int>(parameters.size())) {
    throw std::out_of_range("Residual correlation parameter index is outside its parameter vector.");
  }
  Scalar correlation = parameters[static_cast<std::size_t>(index)];
  if (group.transform == "tanh") correlation = libertad::scalar_tanh(correlation);
  if (!(std::abs(scalar_value(correlation)) < 1.0)) {
    throw std::domain_error("Cross-endpoint residual correlations must be strictly between -1 and 1.");
  }
  return correlation;
}

template <class Scalar>
std::vector<Scalar> residual_grouped_subject_nll_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data,
    const std::vector<Scalar>& prediction, const std::vector<Scalar>& theta,
    const std::vector<Scalar>& sigma,
    const std::vector<Scalar>* variance_prediction = nullptr) {
  Rcpp::NumericVector dv = data["DV"];
  Rcpp::IntegerVector evid = data["EVID"];
  Rcpp::IntegerVector mdv = data["MDV"];
  Rcpp::IntegerVector subjects = data[".ID_INDEX"];
  const bool has_dvid = data.containsElementNamed("DVID");
  int n_subjects = 0;
  for (int value : subjects) n_subjects = std::max(n_subjects, value);
  std::vector<Scalar> result(static_cast<std::size_t>(n_subjects), Scalar(0.0));
  std::vector<bool> processed(static_cast<std::size_t>(data.nrows()), false);
  auto observed = [&](int row) {
    return evid[row] == 0 && mdv[row] == 0 && std::isfinite(dv[row]);
  };
  for (int row = 0; row < data.nrows(); ++row) {
    if (!observed(row) || processed[static_cast<std::size_t>(row)]) continue;
    const int subject = subjects[row] - 1;
    const double dvid = has_dvid ? row_optional(data, "DVID", row, 1.0) : 1.0;
    const int group_index = residual_group_for_dvid(engine, dvid);
    std::vector<int> rows{row};
    if (group_index >= 0) {
      rows.clear();
      const double time = row_optional(data, "TIME", row, 0.0);
      std::vector<bool> endpoint_seen(
        engine.residual_groups[static_cast<std::size_t>(group_index)].dvid.size(), false);
      for (int candidate = 0; candidate < data.nrows(); ++candidate) {
        if (!observed(candidate) || subjects[candidate] - 1 != subject ||
            row_optional(data, "TIME", candidate, 0.0) != time) continue;
        const double candidate_dvid = has_dvid ?
          row_optional(data, "DVID", candidate, 1.0) : 1.0;
        const int endpoint = residual_group_endpoint(
          engine.residual_groups[static_cast<std::size_t>(group_index)], candidate_dvid);
        if (endpoint < 0) continue;
        if (endpoint_seen[static_cast<std::size_t>(endpoint)]) {
          throw std::domain_error(
            "A correlated residual group has duplicate DVID observations at one subject/time.");
        }
        endpoint_seen[static_cast<std::size_t>(endpoint)] = true;
        rows.push_back(candidate);
      }
    }
    const Eigen::Index dimension = static_cast<Eigen::Index>(rows.size());
    VectorT<Scalar> residual(dimension);
    VectorT<Scalar> variance(dimension);
    std::vector<int> endpoint(rows.size(), -1);
    for (Eigen::Index position = 0; position < dimension; ++position) {
      const int selected_row = rows[static_cast<std::size_t>(position)];
      const Scalar f = prediction[static_cast<std::size_t>(selected_row)];
      const Scalar scale_prediction = variance_prediction == nullptr ? f :
        variance_prediction->at(static_cast<std::size_t>(selected_row));
      const double selected_dvid = has_dvid ?
        row_optional(data, "DVID", selected_row, 1.0) : 1.0;
      variance[position] = residual_variance_t(
        engine, scale_prediction, sigma, std::max(1, static_cast<int>(selected_dvid)));
      if (engine.error_type == "exponential") {
        if (!(dv[selected_row] > 0.0)) {
          throw std::domain_error("Exponential residual likelihood requires positive DV.");
        }
        residual[position] = Scalar(std::log(dv[selected_row])) -
          libertad::scalar_log(scalar_floor_t(f, 1e-300));
      } else {
        residual[position] = Scalar(dv[selected_row]) - f;
      }
      if (group_index >= 0) {
        endpoint[static_cast<std::size_t>(position)] = residual_group_endpoint(
          engine.residual_groups[static_cast<std::size_t>(group_index)], selected_dvid);
      }
      processed[static_cast<std::size_t>(selected_row)] = true;
    }
    MatrixT<Scalar> covariance(dimension, dimension);
    for (Eigen::Index first = 0; first < dimension; ++first) {
      for (Eigen::Index second = 0; second < dimension; ++second) {
        Scalar correlation = first == second ? Scalar(1.0) : Scalar(0.0);
        if (group_index >= 0 && first != second) {
          correlation = residual_group_correlation_t(
            engine.residual_groups[static_cast<std::size_t>(group_index)],
            endpoint[static_cast<std::size_t>(first)],
            endpoint[static_cast<std::size_t>(second)], theta, sigma);
        }
        covariance(first, second) = correlation *
          libertad::scalar_sqrt(variance[first] * variance[second]);
      }
    }
    result[static_cast<std::size_t>(subject)] += positive_definite_gaussian_nll_t(
      covariance, residual, "Cross-endpoint residual covariance");
  }
  return result;
}

template <class Scalar>
struct KalmanFilterState {
  VectorT<Scalar> mean;
  MatrixT<Scalar> covariance;
};

template <class Scalar>
struct KalmanRowComponents {
  VectorT<Scalar> predicted_mean;
  MatrixT<Scalar> predicted_covariance;
  VectorT<Scalar> filtered_mean;
  MatrixT<Scalar> filtered_covariance;
  MatrixT<Scalar> transition;
  MatrixT<Scalar> smoother_cross_covariance;
  VectorT<Scalar> observation;
  Scalar observation_variance = Scalar(0.0);
  Scalar innovation = Scalar(0.0);
  Scalar innovation_variance = Scalar(0.0);
  std::vector<VectorT<Scalar>> particle_states;
  std::vector<Scalar> particle_weights;
  std::vector<int> particle_ancestors;
  std::vector<int> particle_regimes;
};

template <class Scalar>
MatrixT<Scalar> state_space_cholesky_t(const MatrixT<Scalar>& covariance,
                                       const char* context) {
  const Eigen::Index dimension = covariance.rows();
  MatrixT<Scalar> lower = MatrixT<Scalar>::Zero(dimension, dimension);
  for (Eigen::Index row = 0; row < dimension; ++row) {
    for (Eigen::Index column = 0; column <= row; ++column) {
      Scalar value = Scalar(0.5) *
        (covariance(row, column) + covariance(column, row));
      for (Eigen::Index k = 0; k < column; ++k) {
        value -= lower(row, k) * lower(column, k);
      }
      if (row == column) {
        if (!(scalar_value(value) > -1e-10)) {
          throw std::domain_error(std::string(context) + " is not positive semidefinite.");
        }
        lower(row, column) = libertad::scalar_sqrt(scalar_floor_t(value, 1e-12));
      } else {
        lower(row, column) = value / lower(column, column);
      }
    }
  }
  return lower;
}

inline std::uint64_t state_space_hash(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31U);
}

inline double state_space_uniform(int seed, int row, int particle,
                                  int dimension, int stream = 0) {
  std::uint64_t value = static_cast<std::uint64_t>(static_cast<std::uint32_t>(seed));
  value ^= static_cast<std::uint64_t>(row + 1) * 0x9e3779b97f4a7c15ULL;
  value ^= static_cast<std::uint64_t>(particle + 1) * 0xbf58476d1ce4e5b9ULL;
  value ^= static_cast<std::uint64_t>(dimension + 1) * 0x94d049bb133111ebULL;
  value ^= static_cast<std::uint64_t>(stream + 1) * 0xd6e8feb86659fd93ULL;
  const std::uint64_t hashed = state_space_hash(value);
  return (static_cast<double>((hashed >> 11U) + 1ULL)) /
    9007199254740993.0;
}

inline double state_space_normal(int seed, int row, int particle,
                                 int dimension, int stream = 0) {
  const double first = state_space_uniform(seed, row, particle, 2 * dimension, stream);
  const double second = state_space_uniform(seed, row, particle, 2 * dimension + 1, stream);
  return std::sqrt(-2.0 * std::log(std::max(first, 1e-15))) *
    std::cos(2.0 * 3.14159265358979323846 * second);
}

template <class Scalar>
struct NonlinearStateOutputs {
  VectorT<Scalar> initial_mean;
  MatrixT<Scalar> initial_covariance;
  VectorT<Scalar> transition;
  MatrixT<Scalar> process_covariance;
  Scalar observation = Scalar(0.0);
  Scalar observation_variance = Scalar(0.0);
};

template <class Scalar>
struct ParticleFilterState {
  std::vector<VectorT<Scalar>> particles;
  std::vector<Scalar> weights;
  std::vector<int> regimes;
};

template <class Scalar>
NonlinearStateOutputs<Scalar> nonlinear_state_raw_outputs_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data, int row, int subject,
    const std::vector<Scalar>& theta, const std::vector<Scalar>& eta,
    int eta_columns, const std::vector<Scalar>& sigma, int mixture_number,
    const Scalar& prediction, bool first, double previous_dv,
    double previous_time, const VectorT<Scalar>& state) {
  ParametersT<Scalar> overrides;
  for (int index = 0; index < engine.kalman_states; ++index) {
    overrides[engine.kalman_state_inputs[static_cast<std::size_t>(index)]] = state[index];
  }
  const std::vector<Scalar> output = evaluate_error_outputs_t(
    engine, data, row, subject, theta, eta, eta_columns, sigma,
    mixture_number, prediction, first, previous_dv, previous_time,
    engine.kalman_outputs, "Nonlinear state-space likelihood", &overrides);
  const int states = engine.kalman_states;
  const std::size_t state_count = static_cast<std::size_t>(states);
  const std::size_t expected = 2U * state_count * state_count +
    2U * state_count + 2U;
  if (output.size() != expected) {
    throw std::logic_error("Compiled nonlinear state-space outputs are inconsistent.");
  }
  std::size_t cursor = 0U;
  NonlinearStateOutputs<Scalar> result;
  result.initial_mean.resize(states);
  result.transition.resize(states);
  for (int index = 0; index < states; ++index) result.initial_mean[index] = output[cursor++];
  auto read_matrix = [&]() {
    MatrixT<Scalar> matrix(states, states);
    for (int matrix_row = 0; matrix_row < states; ++matrix_row) {
      for (int matrix_column = 0; matrix_column < states; ++matrix_column) {
        matrix(matrix_row, matrix_column) = output[cursor++];
      }
    }
    return matrix;
  };
  result.initial_covariance = read_matrix();
  for (int index = 0; index < states; ++index) result.transition[index] = output[cursor++];
  result.process_covariance = read_matrix();
  result.observation = output[cursor++];
  result.observation_variance = output[cursor++];
  return result;
}

template <class Scalar>
NonlinearStateOutputs<Scalar> nonlinear_state_outputs_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data, int row, int subject,
    const std::vector<Scalar>& theta, const std::vector<Scalar>& eta,
    int eta_columns, const std::vector<Scalar>& sigma, int mixture_number,
    const Scalar& prediction, bool first, double previous_dv,
    double previous_time, const VectorT<Scalar>& state) {
  NonlinearStateOutputs<Scalar> result = nonlinear_state_raw_outputs_t(
    engine, data, row, subject, theta, eta, eta_columns, sigma,
    mixture_number, prediction, first, previous_dv, previous_time, state);
  if (engine.kalman_dynamics != "sde" || first) return result;
  const double current_time = row_optional(data, "TIME", row, previous_time);
  const double interval = std::max(0.0, current_time - previous_time);
  const double step = interval / static_cast<double>(engine.kalman_sde_substeps);
  VectorT<Scalar> current = state;
  MatrixT<Scalar> covariance = MatrixT<Scalar>::Zero(
    engine.kalman_states, engine.kalman_states);
  if (step > 0.0) {
    for (int substep = 0; substep < engine.kalman_sde_substeps; ++substep) {
      const NonlinearStateOutputs<Scalar> local = nonlinear_state_raw_outputs_t(
        engine, data, row, subject, theta, eta, eta_columns, sigma,
        mixture_number, prediction, false, previous_dv,
        previous_time + substep * step, current);
      MatrixT<Scalar> drift_jacobian = MatrixT<Scalar>::Zero(
        engine.kalman_states, engine.kalman_states);
      for (int column = 0; column < engine.kalman_states; ++column) {
        const double delta = engine.kalman_jacobian_step *
          std::max(1.0, std::abs(scalar_value(current[column])));
        VectorT<Scalar> plus = current;
        VectorT<Scalar> minus = current;
        plus[column] += Scalar(delta);
        minus[column] -= Scalar(delta);
        const NonlinearStateOutputs<Scalar> upper = nonlinear_state_raw_outputs_t(
          engine, data, row, subject, theta, eta, eta_columns, sigma,
          mixture_number, prediction, false, previous_dv,
          previous_time + substep * step, plus);
        const NonlinearStateOutputs<Scalar> lower = nonlinear_state_raw_outputs_t(
          engine, data, row, subject, theta, eta, eta_columns, sigma,
          mixture_number, prediction, false, previous_dv,
          previous_time + substep * step, minus);
        drift_jacobian.col(column) =
          (upper.transition - lower.transition) / Scalar(2.0 * delta);
      }
      const MatrixT<Scalar> evolution =
        MatrixT<Scalar>::Identity(engine.kalman_states, engine.kalman_states) +
        Scalar(step) * drift_jacobian;
      covariance = evolution * covariance * evolution.transpose() +
        Scalar(step) * local.process_covariance *
        local.process_covariance.transpose();
      current += Scalar(step) * local.transition;
    }
  }
  result.transition = current;
  result.process_covariance = Scalar(0.5) *
    (covariance + covariance.transpose());
  return result;
}

template <class Scalar>
struct SwitchingStateOutputs {
  std::vector<Scalar> initial;
  MatrixT<Scalar> regime_transition;
  std::vector<NonlinearStateOutputs<Scalar>> regime;
};

template <class Scalar>
void normalize_switching_weights(std::vector<Scalar>& value,
                                 const std::string& scale,
                                 const char* context) {
  Scalar total = Scalar(0.0);
  if (scale == "log") {
    double maximum = -std::numeric_limits<double>::infinity();
    for (const Scalar& item : value) maximum = std::max(maximum, scalar_value(item));
    for (Scalar& item : value) {
      item = libertad::scalar_exp(item - Scalar(maximum));
      total += item;
    }
  } else {
    for (const Scalar& item : value) {
      if (scalar_value(item) < 0.0) {
        throw std::domain_error(std::string(context) + " contains a negative weight.");
      }
      total += item;
    }
  }
  if (!(scalar_value(total) > 0.0)) {
    throw std::domain_error(std::string(context) + " has zero total weight.");
  }
  for (Scalar& item : value) item /= total;
}

template <class Scalar>
SwitchingStateOutputs<Scalar> switching_state_raw_outputs_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data, int row, int subject,
    const std::vector<Scalar>& theta, const std::vector<Scalar>& eta,
    int eta_columns, const std::vector<Scalar>& sigma, int mixture_number,
    const Scalar& prediction, bool first, double previous_dv,
    double previous_time, const VectorT<Scalar>& state) {
  ParametersT<Scalar> overrides;
  for (int index = 0; index < engine.kalman_states; ++index) {
    overrides[engine.kalman_state_inputs[static_cast<std::size_t>(index)]] = state[index];
  }
  const std::vector<Scalar> output = evaluate_error_outputs_t(
    engine, data, row, subject, theta, eta, eta_columns, sigma,
    mixture_number, prediction, first, previous_dv, previous_time,
    engine.switching_outputs, "Switching state-space likelihood", &overrides);
  const int regimes = engine.switching_regimes;
  const int states = engine.kalman_states;
  const std::size_t expected = static_cast<std::size_t>(regimes + regimes * regimes +
    regimes * (states + states * states + 2));
  if (output.size() != expected) {
    throw std::logic_error("Compiled switching state-space outputs are inconsistent.");
  }
  std::size_t cursor = 0U;
  SwitchingStateOutputs<Scalar> result;
  result.initial.resize(static_cast<std::size_t>(regimes));
  for (int regime = 0; regime < regimes; ++regime) result.initial[regime] = output[cursor++];
  normalize_switching_weights(
    result.initial, engine.switching_initial_scale, "Initial regime distribution");
  result.regime_transition.resize(regimes, regimes);
  for (int from = 0; from < regimes; ++from) {
    std::vector<Scalar> local(static_cast<std::size_t>(regimes));
    for (int to = 0; to < regimes; ++to) local[to] = output[cursor++];
    normalize_switching_weights(
      local, engine.switching_transition_scale, "Regime transition row");
    for (int to = 0; to < regimes; ++to) result.regime_transition(from, to) = local[to];
  }
  result.regime.resize(static_cast<std::size_t>(regimes));
  for (int regime = 0; regime < regimes; ++regime) {
    NonlinearStateOutputs<Scalar>& local = result.regime[static_cast<std::size_t>(regime)];
    local.transition.resize(states);
    for (int index = 0; index < states; ++index) local.transition[index] = output[cursor++];
    local.process_covariance.resize(states, states);
    for (int matrix_row = 0; matrix_row < states; ++matrix_row) {
      for (int matrix_column = 0; matrix_column < states; ++matrix_column) {
        local.process_covariance(matrix_row, matrix_column) = output[cursor++];
      }
    }
    local.observation = output[cursor++];
    local.observation_variance = output[cursor++];
  }
  return result;
}

template <class Scalar>
int sample_switching_regime(const std::vector<Scalar>& probability,
                            double uniform) {
  Scalar cumulative = probability.front();
  int selected = 0;
  while (selected + 1 < static_cast<int>(probability.size()) &&
         path_lt(cumulative, Scalar(uniform))) {
    ++selected;
    cumulative += probability[static_cast<std::size_t>(selected)];
  }
  return selected;
}

template <class Scalar>
Scalar particle_row_nll_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data, int row, int subject,
    const std::vector<Scalar>& theta, const std::vector<Scalar>& eta,
    int eta_columns, const std::vector<Scalar>& sigma, int mixture_number,
    const Scalar& prediction, bool first, double previous_dv,
    double previous_time, ParticleFilterState<Scalar>& filter,
    KalmanRowComponents<Scalar>* components_output = nullptr) {
  const int states = engine.kalman_states;
  const int particles = engine.kalman_particles;
  const VectorT<Scalar> zero = VectorT<Scalar>::Zero(states);
  const NonlinearStateOutputs<Scalar> base = nonlinear_state_outputs_t(
    engine, data, row, subject, theta, eta, eta_columns, sigma,
    mixture_number, prediction, first, previous_dv, previous_time,
    first ? zero : filter.particles.front());
  MatrixT<Scalar> root = MatrixT<Scalar>::Zero(states, states);
  if (first || !engine.switching_enabled) {
    const MatrixT<Scalar> covariance = first ?
      base.initial_covariance : base.process_covariance;
    root = state_space_cholesky_t(
      covariance, first ? "Particle initial covariance" : "Particle process covariance");
  }
  std::vector<VectorT<Scalar>> propagated(static_cast<std::size_t>(particles));
  std::vector<int> propagated_regime(static_cast<std::size_t>(particles), 0);
  std::vector<Scalar> regime_importance(
    static_cast<std::size_t>(particles), Scalar(1.0));
  std::vector<Scalar> prior_weight(static_cast<std::size_t>(particles),
                                   Scalar(1.0 / static_cast<double>(particles)));
  if (!first) {
    if (filter.particles.size() != static_cast<std::size_t>(particles) ||
        filter.weights.size() != static_cast<std::size_t>(particles) ||
        (engine.switching_enabled &&
         filter.regimes.size() != static_cast<std::size_t>(particles))) {
      throw std::logic_error("Particle filter state dimension changed within a sequence.");
    }
    prior_weight = filter.weights;
  }
  for (int particle = 0; particle < particles; ++particle) {
    VectorT<Scalar> normal(states);
    for (int state = 0; state < states; ++state) {
      normal[state] = Scalar(state_space_normal(
        engine.kalman_seed + 7919 * subject, row, particle, state, first ? 1 : 2));
    }
    if (engine.switching_enabled) {
      const VectorT<Scalar>& current = first ? zero :
        filter.particles[static_cast<std::size_t>(particle)];
      const SwitchingStateOutputs<Scalar> switching = switching_state_raw_outputs_t(
        engine, data, row, subject, theta, eta, eta_columns, sigma,
        mixture_number, prediction, first, previous_dv, previous_time, current);
      std::vector<Scalar> regime_probability;
      if (first) {
        regime_probability = switching.initial;
      } else {
        regime_probability.resize(static_cast<std::size_t>(engine.switching_regimes));
        const int previous_regime = filter.regimes[static_cast<std::size_t>(particle)];
        for (int regime = 0; regime < engine.switching_regimes; ++regime) {
          regime_probability[static_cast<std::size_t>(regime)] =
            switching.regime_transition(previous_regime, regime);
        }
      }
      const int target_regime = particle % engine.switching_regimes;
      propagated_regime[static_cast<std::size_t>(particle)] = target_regime;
      const int target_count = particles / engine.switching_regimes +
        (target_regime < particles % engine.switching_regimes ? 1 : 0);
      const double proposal_probability =
        static_cast<double>(target_count) / static_cast<double>(particles);
      regime_importance[static_cast<std::size_t>(particle)] =
        regime_probability[static_cast<std::size_t>(target_regime)] /
        Scalar(proposal_probability);
    }
    if (first) {
      propagated[static_cast<std::size_t>(particle)] = base.initial_mean + root * normal;
    } else if (engine.kalman_dynamics == "sde") {
      VectorT<Scalar> current = filter.particles[static_cast<std::size_t>(particle)];
      const double current_time = row_optional(data, "TIME", row, previous_time);
      const double interval = std::max(0.0, current_time - previous_time);
      const double step = interval / static_cast<double>(engine.kalman_sde_substeps);
      for (int substep = 0; substep < engine.kalman_sde_substeps; ++substep) {
        const NonlinearStateOutputs<Scalar> local = engine.switching_enabled ?
          switching_state_raw_outputs_t(
            engine, data, row, subject, theta, eta, eta_columns, sigma,
            mixture_number, prediction, false, previous_dv,
            previous_time + substep * step, current).regime[
              static_cast<std::size_t>(propagated_regime[static_cast<std::size_t>(particle)])] :
          nonlinear_state_raw_outputs_t(
            engine, data, row, subject, theta, eta, eta_columns, sigma,
            mixture_number, prediction, false, previous_dv,
            previous_time + substep * step, current);
        VectorT<Scalar> increment(states);
        for (int state = 0; state < states; ++state) {
          increment[state] = Scalar(std::sqrt(step) * state_space_normal(
            engine.kalman_seed + 7919 * subject, row, particle,
            substep * states + state, 6));
        }
        VectorT<Scalar> next = current + Scalar(step) * local.transition +
          local.process_covariance * increment;
        if (engine.kalman_sde_method == "milstein") {
          for (int first_state = 0; first_state < states; ++first_state) {
            for (int second_state = 0; second_state < states; ++second_state) {
              if (first_state != second_state &&
                  std::abs(scalar_value(local.process_covariance(first_state, second_state))) > 1e-12) {
                throw std::domain_error(
                  "Milstein propagation currently requires diagonal diffusion.");
              }
            }
            const double derivative_step = engine.kalman_jacobian_step *
              std::max(1.0, std::abs(scalar_value(current[first_state])));
            VectorT<Scalar> plus = current;
            VectorT<Scalar> minus = current;
            plus[first_state] += Scalar(derivative_step);
            minus[first_state] -= Scalar(derivative_step);
            const int regime = propagated_regime[static_cast<std::size_t>(particle)];
            const Scalar upper = engine.switching_enabled ?
              switching_state_raw_outputs_t(
                engine, data, row, subject, theta, eta, eta_columns, sigma,
                mixture_number, prediction, false, previous_dv,
                previous_time + substep * step, plus).regime[
                  static_cast<std::size_t>(regime)].process_covariance(first_state, first_state) :
              nonlinear_state_raw_outputs_t(
                engine, data, row, subject, theta, eta, eta_columns, sigma,
                mixture_number, prediction, false, previous_dv,
                previous_time + substep * step, plus).process_covariance(first_state, first_state);
            const Scalar lower = engine.switching_enabled ?
              switching_state_raw_outputs_t(
                engine, data, row, subject, theta, eta, eta_columns, sigma,
                mixture_number, prediction, false, previous_dv,
                previous_time + substep * step, minus).regime[
                  static_cast<std::size_t>(regime)].process_covariance(first_state, first_state) :
              nonlinear_state_raw_outputs_t(
                engine, data, row, subject, theta, eta, eta_columns, sigma,
                mixture_number, prediction, false, previous_dv,
                previous_time + substep * step, minus).process_covariance(first_state, first_state);
            const Scalar derivative = (upper - lower) / Scalar(2.0 * derivative_step);
            const Scalar diffusion = local.process_covariance(first_state, first_state);
            next[first_state] += Scalar(0.5) * diffusion * derivative *
              (increment[first_state] * increment[first_state] - Scalar(step));
          }
        }
        current = next;
      }
      propagated[static_cast<std::size_t>(particle)] = current;
    } else {
      const VectorT<Scalar>& current = filter.particles[static_cast<std::size_t>(particle)];
      const NonlinearStateOutputs<Scalar> point = engine.switching_enabled ?
        switching_state_raw_outputs_t(
          engine, data, row, subject, theta, eta, eta_columns, sigma,
          mixture_number, prediction, first, previous_dv, previous_time, current).regime[
            static_cast<std::size_t>(propagated_regime[static_cast<std::size_t>(particle)])] :
        nonlinear_state_outputs_t(
          engine, data, row, subject, theta, eta, eta_columns, sigma,
          mixture_number, prediction, first, previous_dv, previous_time, current);
      const MatrixT<Scalar> local_root = engine.switching_enabled ?
        state_space_cholesky_t(point.process_covariance,
                               "Switching particle process covariance") : root;
      propagated[static_cast<std::size_t>(particle)] =
        point.transition + local_root * normal;
    }
  }
  KalmanRowComponents<Scalar> components;
  std::vector<Scalar> predictive_weight = prior_weight;
  if (engine.switching_enabled) {
    Scalar predictive_total = Scalar(0.0);
    for (int particle = 0; particle < particles; ++particle) {
      predictive_weight[static_cast<std::size_t>(particle)] *=
        regime_importance[static_cast<std::size_t>(particle)];
      predictive_total += predictive_weight[static_cast<std::size_t>(particle)];
    }
    if (!(scalar_value(predictive_total) > 0.0)) {
      throw std::domain_error("Switching-state predictive regime mass is zero.");
    }
    for (Scalar& weight : predictive_weight) weight /= predictive_total;
  }
  components.predicted_mean = VectorT<Scalar>::Zero(states);
  for (int particle = 0; particle < particles; ++particle) {
    components.predicted_mean += predictive_weight[static_cast<std::size_t>(particle)] *
      propagated[static_cast<std::size_t>(particle)];
  }
  components.predicted_covariance = MatrixT<Scalar>::Zero(states, states);
  for (int particle = 0; particle < particles; ++particle) {
    const VectorT<Scalar> centered = propagated[static_cast<std::size_t>(particle)] -
      components.predicted_mean;
    components.predicted_covariance += predictive_weight[static_cast<std::size_t>(particle)] *
      centered * centered.transpose();
  }
  std::vector<Scalar> log_weight(static_cast<std::size_t>(particles));
  std::vector<Scalar> observation_value(static_cast<std::size_t>(particles));
  double maximum = -std::numeric_limits<double>::infinity();
  const Scalar baseline = engine.kalman_prediction_baseline ? prediction : Scalar(0.0);
  const Scalar observed = Scalar(row_optional(data, "DV", row, NA_REAL));
  Scalar predicted_observation = Scalar(0.0);
  for (int particle = 0; particle < particles; ++particle) {
    const NonlinearStateOutputs<Scalar> point = engine.switching_enabled ?
      switching_state_raw_outputs_t(
        engine, data, row, subject, theta, eta, eta_columns, sigma,
        mixture_number, prediction, first, previous_dv, previous_time,
        propagated[static_cast<std::size_t>(particle)]).regime[
          static_cast<std::size_t>(propagated_regime[static_cast<std::size_t>(particle)])] :
      nonlinear_state_outputs_t(
        engine, data, row, subject, theta, eta, eta_columns, sigma,
        mixture_number, prediction, first, previous_dv, previous_time,
        propagated[static_cast<std::size_t>(particle)]);
    if (!(scalar_value(point.observation_variance) > 1e-14)) {
      throw std::domain_error("Particle observation variance must be positive.");
    }
    observation_value[static_cast<std::size_t>(particle)] = baseline + point.observation;
    predicted_observation += predictive_weight[static_cast<std::size_t>(particle)] *
      observation_value[static_cast<std::size_t>(particle)];
    const Scalar residual = observed - observation_value[static_cast<std::size_t>(particle)];
    log_weight[static_cast<std::size_t>(particle)] =
      libertad::scalar_log(scalar_floor_t(
        prior_weight[static_cast<std::size_t>(particle)] *
          regime_importance[static_cast<std::size_t>(particle)], 1e-300)) -
      Scalar(0.5) * (libertad::scalar_log(point.observation_variance) +
                     residual * residual / point.observation_variance);
    maximum = std::max(maximum, scalar_value(log_weight[static_cast<std::size_t>(particle)]));
  }
  Scalar normalizer = Scalar(0.0);
  for (const Scalar& value : log_weight) normalizer += libertad::scalar_exp(value - Scalar(maximum));
  const Scalar log_likelihood = Scalar(maximum) + libertad::scalar_log(normalizer);
  filter.weights.resize(static_cast<std::size_t>(particles));
  for (int particle = 0; particle < particles; ++particle) {
    filter.weights[static_cast<std::size_t>(particle)] =
      libertad::scalar_exp(log_weight[static_cast<std::size_t>(particle)] - log_likelihood);
  }
  components.innovation = observed - predicted_observation;
  components.innovation_variance = Scalar(0.0);
  for (int particle = 0; particle < particles; ++particle) {
    const Scalar centered = observation_value[static_cast<std::size_t>(particle)] -
      predicted_observation;
    components.innovation_variance += predictive_weight[static_cast<std::size_t>(particle)] *
      centered * centered;
  }
  components.transition = MatrixT<Scalar>::Identity(states, states);
  components.smoother_cross_covariance = MatrixT<Scalar>::Zero(states, states);
  VectorT<Scalar> filtered_mean = VectorT<Scalar>::Zero(states);
  for (int particle = 0; particle < particles; ++particle) {
    filtered_mean += filter.weights[static_cast<std::size_t>(particle)] *
      propagated[static_cast<std::size_t>(particle)];
  }
  MatrixT<Scalar> filtered_covariance = MatrixT<Scalar>::Zero(states, states);
  for (int particle = 0; particle < particles; ++particle) {
    const VectorT<Scalar> centered = propagated[static_cast<std::size_t>(particle)] - filtered_mean;
    filtered_covariance += filter.weights[static_cast<std::size_t>(particle)] *
      centered * centered.transpose();
  }
  Scalar squared_weight = Scalar(0.0);
  for (const Scalar& weight : filter.weights) squared_weight += weight * weight;
  std::vector<int> ancestors(static_cast<std::size_t>(particles));
  std::iota(ancestors.begin(), ancestors.end(), 0);
  // Keep both the resampling decision and systematic ancestry comparisons on
  // the CppAD comparison tape. A changed ancestry then raises TapePathChange
  // on replay and the population engine records the correct particle graph at
  // the new parameter point instead of silently reusing stale indices.
  const Scalar resampling_boundary = Scalar(
    1.0 / (engine.kalman_ess_threshold * static_cast<double>(particles)));
  if (squared_weight > resampling_boundary) {
    const double offset = state_space_uniform(
      engine.kalman_seed + 104729 * subject, row, 0, 0, 9) /
      static_cast<double>(particles);
    std::vector<VectorT<Scalar>> resampled(static_cast<std::size_t>(particles));
    std::vector<int> resampled_regime(static_cast<std::size_t>(particles));
    int selected = 0;
    Scalar cumulative = filter.weights[0];
    for (int particle = 0; particle < particles; ++particle) {
      const double target = offset + static_cast<double>(particle) /
        static_cast<double>(particles);
      while (selected + 1 < particles && cumulative < Scalar(target)) {
        ++selected;
        cumulative += filter.weights[static_cast<std::size_t>(selected)];
      }
      resampled[static_cast<std::size_t>(particle)] =
        propagated[static_cast<std::size_t>(selected)];
      ancestors[static_cast<std::size_t>(particle)] = selected;
      resampled_regime[static_cast<std::size_t>(particle)] =
        propagated_regime[static_cast<std::size_t>(selected)];
    }
    filter.particles = std::move(resampled);
    if (engine.switching_enabled) filter.regimes = std::move(resampled_regime);
    std::fill(filter.weights.begin(), filter.weights.end(),
              Scalar(1.0 / static_cast<double>(particles)));
  } else {
    filter.particles = std::move(propagated);
    if (engine.switching_enabled) filter.regimes = std::move(propagated_regime);
  }
  if (components_output != nullptr) {
    components.observation = VectorT<Scalar>::Zero(states);
    components.filtered_mean = filtered_mean;
    components.filtered_covariance = filtered_covariance;
    components.particle_states = filter.particles;
    components.particle_weights = filter.weights;
    components.particle_ancestors = ancestors;
    components.particle_regimes = filter.regimes;
    *components_output = components;
  }
  return Scalar(-2.0) * log_likelihood;
}

template <class Scalar>
Scalar nonlinear_kalman_row_nll_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data, int row, int subject,
    const std::vector<Scalar>& theta, const std::vector<Scalar>& eta,
    int eta_columns, const std::vector<Scalar>& sigma, int mixture_number,
    const Scalar& prediction, bool first, double previous_dv,
    double previous_time, KalmanFilterState<Scalar>& filtered,
    KalmanRowComponents<Scalar>* components_output) {
  const int states = engine.kalman_states;
  const VectorT<Scalar> zero = VectorT<Scalar>::Zero(states);
  NonlinearStateOutputs<Scalar> base = nonlinear_state_outputs_t(
    engine, data, row, subject, theta, eta, eta_columns, sigma,
    mixture_number, prediction, first, previous_dv, previous_time,
    first ? zero : filtered.mean);
  KalmanRowComponents<Scalar> components;
  components.transition = MatrixT<Scalar>::Identity(states, states);
  components.smoother_cross_covariance = MatrixT<Scalar>::Zero(states, states);
  if (first) {
    components.predicted_mean = base.initial_mean;
    components.predicted_covariance = Scalar(0.5) *
      (base.initial_covariance + base.initial_covariance.transpose());
  } else if (engine.kalman_filter_type == "ekf") {
    const VectorT<Scalar> previous_mean = filtered.mean;
    const MatrixT<Scalar> previous_covariance = filtered.covariance;
    components.predicted_mean = base.transition;
    MatrixT<Scalar> jacobian(states, states);
    for (int column = 0; column < states; ++column) {
      const double step = engine.kalman_jacobian_step *
        std::max(1.0, std::abs(scalar_value(previous_mean[column])));
      VectorT<Scalar> plus = previous_mean;
      VectorT<Scalar> minus = previous_mean;
      plus[column] += Scalar(step);
      minus[column] -= Scalar(step);
      const VectorT<Scalar> upper = nonlinear_state_outputs_t(
        engine, data, row, subject, theta, eta, eta_columns, sigma,
        mixture_number, prediction, first, previous_dv, previous_time, plus).transition;
      const VectorT<Scalar> lower = nonlinear_state_outputs_t(
        engine, data, row, subject, theta, eta, eta_columns, sigma,
        mixture_number, prediction, first, previous_dv, previous_time, minus).transition;
      jacobian.col(column) = (upper - lower) / Scalar(2.0 * step);
    }
    components.transition = jacobian;
    components.smoother_cross_covariance = previous_covariance * jacobian.transpose();
    components.predicted_covariance = jacobian * previous_covariance *
      jacobian.transpose() + base.process_covariance;
  } else if (engine.kalman_filter_type == "ukf") {
    const VectorT<Scalar> previous_mean = filtered.mean;
    const MatrixT<Scalar> previous_covariance = filtered.covariance;
    const double lambda = engine.kalman_ukf_alpha * engine.kalman_ukf_alpha *
      (static_cast<double>(states) + engine.kalman_ukf_kappa) -
      static_cast<double>(states);
    const double scale = static_cast<double>(states) + lambda;
    const int points = 2 * states + 1;
    std::vector<double> mean_weight(static_cast<std::size_t>(points),
                                    1.0 / (2.0 * scale));
    std::vector<double> covariance_weight = mean_weight;
    mean_weight[0] = lambda / scale;
    covariance_weight[0] = mean_weight[0] +
      (1.0 - engine.kalman_ukf_alpha * engine.kalman_ukf_alpha +
       engine.kalman_ukf_beta);
    const MatrixT<Scalar> root = state_space_cholesky_t(
      previous_covariance, "UKF filtered covariance") * Scalar(std::sqrt(scale));
    std::vector<VectorT<Scalar>> sigma_points(static_cast<std::size_t>(points));
    std::vector<VectorT<Scalar>> propagated(static_cast<std::size_t>(points));
    sigma_points[0] = previous_mean;
    for (int column = 0; column < states; ++column) {
      sigma_points[static_cast<std::size_t>(1 + column)] = previous_mean + root.col(column);
      sigma_points[static_cast<std::size_t>(1 + states + column)] = previous_mean - root.col(column);
    }
    components.predicted_mean = VectorT<Scalar>::Zero(states);
    for (int point = 0; point < points; ++point) {
      propagated[static_cast<std::size_t>(point)] = nonlinear_state_outputs_t(
        engine, data, row, subject, theta, eta, eta_columns, sigma,
        mixture_number, prediction, first, previous_dv, previous_time,
        sigma_points[static_cast<std::size_t>(point)]).transition;
      components.predicted_mean += Scalar(mean_weight[static_cast<std::size_t>(point)]) *
        propagated[static_cast<std::size_t>(point)];
    }
    components.predicted_covariance = base.process_covariance;
    for (int point = 0; point < points; ++point) {
      const VectorT<Scalar> before = sigma_points[static_cast<std::size_t>(point)] - previous_mean;
      const VectorT<Scalar> after = propagated[static_cast<std::size_t>(point)] -
        components.predicted_mean;
      components.predicted_covariance +=
        Scalar(covariance_weight[static_cast<std::size_t>(point)]) *
        after * after.transpose();
      components.smoother_cross_covariance +=
        Scalar(covariance_weight[static_cast<std::size_t>(point)]) *
        before * after.transpose();
    }
  } else {
    throw std::logic_error("Particle filtering uses the dedicated particle state path.");
  }
  components.predicted_covariance = Scalar(0.5) *
    (components.predicted_covariance + components.predicted_covariance.transpose());

  VectorT<Scalar> observation_gradient(states);
  Scalar observation_mean = Scalar(0.0);
  if (engine.kalman_filter_type == "ukf") {
    const double lambda = engine.kalman_ukf_alpha * engine.kalman_ukf_alpha *
      (static_cast<double>(states) + engine.kalman_ukf_kappa) -
      static_cast<double>(states);
    const double scale = static_cast<double>(states) + lambda;
    const int points = 2 * states + 1;
    std::vector<double> mean_weight(static_cast<std::size_t>(points),
                                    1.0 / (2.0 * scale));
    std::vector<double> covariance_weight = mean_weight;
    mean_weight[0] = lambda / scale;
    covariance_weight[0] = mean_weight[0] +
      (1.0 - engine.kalman_ukf_alpha * engine.kalman_ukf_alpha +
       engine.kalman_ukf_beta);
    const MatrixT<Scalar> root = state_space_cholesky_t(
      components.predicted_covariance, "UKF predicted covariance") *
      Scalar(std::sqrt(scale));
    std::vector<VectorT<Scalar>> points_state(static_cast<std::size_t>(points));
    std::vector<Scalar> values(static_cast<std::size_t>(points));
    points_state[0] = components.predicted_mean;
    for (int column = 0; column < states; ++column) {
      points_state[static_cast<std::size_t>(1 + column)] =
        components.predicted_mean + root.col(column);
      points_state[static_cast<std::size_t>(1 + states + column)] =
        components.predicted_mean - root.col(column);
    }
    for (int point = 0; point < points; ++point) {
      values[static_cast<std::size_t>(point)] = nonlinear_state_outputs_t(
        engine, data, row, subject, theta, eta, eta_columns, sigma,
        mixture_number, prediction, first, previous_dv, previous_time,
        points_state[static_cast<std::size_t>(point)]).observation;
      observation_mean += Scalar(mean_weight[static_cast<std::size_t>(point)]) *
        values[static_cast<std::size_t>(point)];
    }
    components.innovation_variance = base.observation_variance;
    VectorT<Scalar> cross = VectorT<Scalar>::Zero(states);
    for (int point = 0; point < points; ++point) {
      const Scalar centered = values[static_cast<std::size_t>(point)] - observation_mean;
      components.innovation_variance +=
        Scalar(covariance_weight[static_cast<std::size_t>(point)]) * centered * centered;
      cross += Scalar(covariance_weight[static_cast<std::size_t>(point)]) *
        (points_state[static_cast<std::size_t>(point)] - components.predicted_mean) * centered;
    }
    observation_gradient = cross / components.innovation_variance;
  } else {
    const NonlinearStateOutputs<Scalar> at_mean = nonlinear_state_outputs_t(
      engine, data, row, subject, theta, eta, eta_columns, sigma,
      mixture_number, prediction, first, previous_dv, previous_time,
      components.predicted_mean);
    observation_mean = at_mean.observation;
    for (int column = 0; column < states; ++column) {
      const double step = engine.kalman_jacobian_step *
        std::max(1.0, std::abs(scalar_value(components.predicted_mean[column])));
      VectorT<Scalar> plus = components.predicted_mean;
      VectorT<Scalar> minus = components.predicted_mean;
      plus[column] += Scalar(step);
      minus[column] -= Scalar(step);
      const Scalar upper = nonlinear_state_outputs_t(
        engine, data, row, subject, theta, eta, eta_columns, sigma,
        mixture_number, prediction, first, previous_dv, previous_time, plus).observation;
      const Scalar lower = nonlinear_state_outputs_t(
        engine, data, row, subject, theta, eta, eta_columns, sigma,
        mixture_number, prediction, first, previous_dv, previous_time, minus).observation;
      observation_gradient[column] = (upper - lower) / Scalar(2.0 * step);
    }
    components.innovation_variance =
      (observation_gradient.transpose() * components.predicted_covariance *
       observation_gradient)[0] + at_mean.observation_variance;
    observation_gradient = components.predicted_covariance *
      observation_gradient / components.innovation_variance;
  }
  const Scalar baseline = engine.kalman_prediction_baseline ? prediction : Scalar(0.0);
  components.innovation = Scalar(row_optional(data, "DV", row, NA_REAL)) -
    baseline - observation_mean;
  if (!(scalar_value(components.innovation_variance) > 1e-14)) {
    throw std::domain_error("Nonlinear state-space innovation variance must be positive.");
  }
  const VectorT<Scalar> gain = observation_gradient;
  filtered.mean = components.predicted_mean + gain * components.innovation;
  filtered.covariance = components.predicted_covariance -
    gain * components.innovation_variance * gain.transpose();
  filtered.covariance = Scalar(0.5) *
    (filtered.covariance + filtered.covariance.transpose());
  components.filtered_mean = filtered.mean;
  components.filtered_covariance = filtered.covariance;
  if (components_output != nullptr) *components_output = components;
  return libertad::scalar_log(components.innovation_variance) +
    components.innovation * components.innovation /
      components.innovation_variance;
}

template <class Scalar>
Scalar kalman_row_nll_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data, int row, int subject,
    const std::vector<Scalar>& theta, const std::vector<Scalar>& eta,
    int eta_columns, const std::vector<Scalar>& sigma, int mixture_number,
    const Scalar& prediction, bool first, double previous_dv,
    double previous_time, KalmanFilterState<Scalar>& filtered,
    KalmanRowComponents<Scalar>* components_output = nullptr) {
  if (engine.kalman_filter_type != "linear") {
    return nonlinear_kalman_row_nll_t(
      engine, data, row, subject, theta, eta, eta_columns, sigma,
      mixture_number, prediction, first, previous_dv, previous_time,
      filtered, components_output);
  }
  const int states = engine.kalman_states;
  const std::size_t state_count = static_cast<std::size_t>(states);
  const std::size_t matrix_count = state_count * state_count;
  const std::size_t expected = 2U * state_count + 3U * matrix_count + 1U;
  if (!engine.kalman_enabled || states < 1 || engine.kalman_outputs.size() != expected) {
    throw std::logic_error("Compiled Kalman state-space likelihood is inconsistent.");
  }
  const std::vector<Scalar> output = evaluate_error_outputs_t(
    engine, data, row, subject, theta, eta, eta_columns, sigma,
    mixture_number, prediction, first, previous_dv, previous_time,
    engine.kalman_outputs, "Kalman state-space likelihood");
  std::size_t cursor = 0U;
  VectorT<Scalar> initial_mean(states);
  for (int state = 0; state < states; ++state) initial_mean[state] = output[cursor++];
  auto read_matrix = [&]() {
    MatrixT<Scalar> matrix(states, states);
    for (int matrix_row = 0; matrix_row < states; ++matrix_row) {
      for (int matrix_column = 0; matrix_column < states; ++matrix_column) {
        matrix(matrix_row, matrix_column) = output[cursor++];
      }
    }
    return matrix;
  };
  const MatrixT<Scalar> initial_covariance = read_matrix();
  const MatrixT<Scalar> transition = read_matrix();
  const MatrixT<Scalar> process_covariance = read_matrix();
  VectorT<Scalar> observation(states);
  for (int state = 0; state < states; ++state) observation[state] = output[cursor++];
  const Scalar observation_variance = output[cursor++];
  if (!(scalar_value(observation_variance) >= 0.0)) {
    throw std::domain_error(
      "Kalman observation variance must be non-negative at row " +
      std::to_string(row + 1) + ".");
  }
  KalmanRowComponents<Scalar> components;
  components.transition = transition;
  components.smoother_cross_covariance = MatrixT<Scalar>::Zero(states, states);
  components.observation = observation;
  components.observation_variance = observation_variance;
  if (first) {
    components.predicted_mean = initial_mean;
    components.predicted_covariance = Scalar(0.5) *
      (initial_covariance + initial_covariance.transpose());
  } else {
    if (filtered.mean.size() != states || filtered.covariance.rows() != states ||
        filtered.covariance.cols() != states) {
      throw std::logic_error("Kalman filter state has the wrong dimension.");
    }
    components.smoother_cross_covariance = filtered.covariance * transition.transpose();
    components.predicted_mean = transition * filtered.mean;
    components.predicted_covariance = transition * filtered.covariance *
      transition.transpose() + process_covariance;
    components.predicted_covariance = Scalar(0.5) *
      (components.predicted_covariance + components.predicted_covariance.transpose());
  }
  const Scalar baseline = engine.kalman_prediction_baseline ? prediction : Scalar(0.0);
  const Scalar observation_mean = baseline +
    (observation.transpose() * components.predicted_mean)[0];
  components.innovation = Scalar(row_optional(data, "DV", row, NA_REAL)) -
    observation_mean;
  components.innovation_variance =
    (observation.transpose() * components.predicted_covariance * observation)[0] +
    observation_variance;
  if (!(scalar_value(components.innovation_variance) > 1e-14)) {
    throw std::domain_error(
      "Kalman innovation variance is not positive at row " +
      std::to_string(row + 1) + ".");
  }
  const VectorT<Scalar> gain =
    components.predicted_covariance * observation / components.innovation_variance;
  filtered.mean = components.predicted_mean + gain * components.innovation;
  const MatrixT<Scalar> identity = MatrixT<Scalar>::Identity(states, states);
  const MatrixT<Scalar> update = identity - gain * observation.transpose();
  filtered.covariance = update * components.predicted_covariance * update.transpose() +
    gain * observation_variance * gain.transpose();
  filtered.covariance = Scalar(0.5) *
    (filtered.covariance + filtered.covariance.transpose());
  components.filtered_mean = filtered.mean;
  components.filtered_covariance = filtered.covariance;
  if (components_output != nullptr) *components_output = components;
  return libertad::scalar_log(components.innovation_variance) +
    components.innovation * components.innovation /
      components.innovation_variance;
}

template <class Scalar>
Scalar ar1_rho_t(const ModelEngine& engine,
                 const std::vector<Scalar>& theta,
                 const std::vector<Scalar>& sigma) {
  if (engine.ar1_parameter_source == "fixed") return Scalar(engine.ar1_rho);
  const std::vector<Scalar>& parameters =
    engine.ar1_parameter_source == "theta" ? theta : sigma;
  if (engine.ar1_parameter_index < 0 ||
      engine.ar1_parameter_index >= static_cast<int>(parameters.size())) {
    throw std::out_of_range("Estimated AR(1) parameter index is outside the parameter vector.");
  }
  Scalar rho = parameters[static_cast<std::size_t>(engine.ar1_parameter_index)];
  if (engine.ar1_transform == "tanh") rho = libertad::scalar_tanh(rho);
  if (!(std::abs(scalar_value(rho)) < 1.0)) {
    throw std::domain_error("Estimated AR(1) correlation must be strictly between -1 and 1.");
  }
  return rho;
}

template <class Scalar>
std::vector<Scalar> residual_subject_nll_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data,
    const std::vector<Scalar>& prediction, const std::vector<Scalar>& theta,
    const std::vector<Scalar>& eta, const std::vector<Scalar>& sigma,
    const std::vector<int>& mixture_assignment = std::vector<int>(),
    const std::vector<Scalar>* variance_prediction = nullptr) {
  if (!engine.residual_groups.empty()) {
    if (engine.error_type == "likelihood" || engine.sigma_correlation != "independent" ||
        engine.blq_method != "none") {
      throw std::logic_error("Correlated residual groups were combined with an incompatible likelihood option.");
    }
    return residual_grouped_subject_nll_t(
      engine, data, prediction, theta, sigma, variance_prediction);
  }
  Rcpp::NumericVector dv = data["DV"];
  Rcpp::IntegerVector evid = data["EVID"];
  Rcpp::IntegerVector mdv = data["MDV"];
  Rcpp::IntegerVector subjects = data[".ID_INDEX"];
  int n_subjects = 0;
  for (int value : subjects) n_subjects = std::max(n_subjects, value);
  std::vector<Scalar> result(static_cast<std::size_t>(n_subjects), Scalar(0.0));
  const bool has_dvid = data.containsElementNamed("DVID");
  const bool has_lloq = data.containsElementNamed("LLOQ");
  const bool has_blq = data.containsElementNamed("BLQ");
  const bool has_cens = data.containsElementNamed("CENS");
  int previous_subject = -1;
  std::unordered_map<int, Scalar> previous_standardized_residual;
  std::unordered_map<int, double> previous_outcome;
  std::unordered_map<int, double> previous_outcome_time;
  std::unordered_map<int, std::vector<Scalar>> hmm_filtered;
  std::unordered_map<int, KalmanFilterState<Scalar>> kalman_filtered;
  std::unordered_map<int, ParticleFilterState<Scalar>> particle_filtered;
  const int eta_columns = n_subjects > 0 ?
    static_cast<int>(eta.size() / static_cast<std::size_t>(n_subjects)) : 0;
  const Scalar ar1_rho = ar1_rho_t(engine, theta, sigma);

  for (int row = 0; row < data.nrows(); ++row) {
    const int subject = subjects[row] - 1;
    const int dvid = has_dvid ?
      static_cast<int>(row_optional(data, "DVID", row, 1.0)) : 1;
    if (subject != previous_subject) {
      previous_outcome.clear();
      previous_outcome_time.clear();
      hmm_filtered.clear();
      kalman_filtered.clear();
      particle_filtered.clear();
      previous_standardized_residual.clear();
    }
    if (evid[row] == 0 && mdv[row] == 0 && std::isfinite(dv[row])) {
      const Scalar f = prediction[static_cast<std::size_t>(row)];
      if (engine.error_type == "likelihood") {
        const int sequence =
          (engine.hmm_enabled && !engine.hmm_by_dvid) ||
          (engine.kalman_enabled && !engine.kalman_by_dvid) ? 1 : dvid;
        const auto previous = previous_outcome.find(sequence);
        const bool first = previous == previous_outcome.end();
        const double previous_value = first ? dv[row] : previous->second;
        const double previous_time = first ? row_optional(data, "TIME", row, 0.0) :
          previous_outcome_time.at(sequence);
        const int mixture_number = mixture_assignment.empty() ?
          static_cast<int>(row_optional(data, "MIXNUM", row, 1.0)) :
          mixture_assignment.at(static_cast<std::size_t>(subject));
        if (engine.hmm_enabled) {
          result[static_cast<std::size_t>(subject)] += hmm_row_nll_t(
            engine, data, row, subject, theta, eta, eta_columns, sigma,
            mixture_number, f, first, previous_value, previous_time,
            hmm_filtered[sequence]);
        } else if (engine.kalman_enabled) {
          if (engine.kalman_filter_type == "particle") {
            result[static_cast<std::size_t>(subject)] += particle_row_nll_t(
              engine, data, row, subject, theta, eta, eta_columns, sigma,
              mixture_number, f, first, previous_value, previous_time,
              particle_filtered[sequence]);
          } else {
            result[static_cast<std::size_t>(subject)] += kalman_row_nll_t(
              engine, data, row, subject, theta, eta, eta_columns, sigma,
              mixture_number, f, first, previous_value, previous_time,
              kalman_filtered[sequence]);
          }
        } else {
          result[static_cast<std::size_t>(subject)] += user_likelihood_nll_t(
            engine, data, row, subject, theta, eta, eta_columns, sigma,
            mixture_number, f, first, previous_value, previous_time);
        }
        previous_standardized_residual.erase(dvid);
      } else {
      const Scalar scale_prediction = variance_prediction == nullptr ? f :
        variance_prediction->at(static_cast<std::size_t>(row));
      const Scalar variance = residual_variance_t(
        engine, scale_prediction, sigma, dvid);
      const Scalar sd = libertad::scalar_sqrt(variance);
      double limit = has_lloq ? row_optional(data, "LLOQ", row, engine.lloq) : engine.lloq;
      bool censored = false;
      if (engine.blq_method != "none" && std::isfinite(limit)) {
        if (has_blq) censored = row_optional(data, "BLQ", row, 0.0) == 1.0;
        if (has_cens) censored = censored || row_optional(data, "CENS", row, 0.0) == 1.0;
        if (!has_blq && !has_cens) censored = dv[row] < limit;
      }
      if (censored) {
        Scalar z;
        if (engine.error_type == "exponential") {
          z = (Scalar(std::log(std::max(limit, 1e-300))) -
            libertad::scalar_log(scalar_floor_t(f, 1e-300))) / sd;
        } else {
          z = (Scalar(limit) - f) / sd;
        }
        Scalar probability = normal_cdf_t(z);
        if (engine.blq_method == "m4" && engine.error_type != "exponential") {
          const Scalar below_zero = normal_cdf_t((-f) / sd);
          probability = scalar_floor_t(
            (probability - below_zero) /
              scalar_floor_t(Scalar(1.0) - below_zero, 1e-300),
            1e-300);
        }
        result[static_cast<std::size_t>(subject)] -=
          Scalar(2.0) * libertad::scalar_log(probability);
        previous_standardized_residual.erase(dvid);
      } else {
        Scalar residual = Scalar(dv[row]) - f;
        if (engine.error_type == "exponential") {
          if (!(dv[row] > 0.0)) {
            throw std::domain_error("Exponential residual likelihood requires positive DV.");
          }
          residual = Scalar(std::log(dv[row])) -
            libertad::scalar_log(scalar_floor_t(f, 1e-300));
        }
        const auto previous_residual = previous_standardized_residual.find(dvid);
        if (engine.sigma_correlation == "ar1" &&
            previous_residual != previous_standardized_residual.end()) {
          const Scalar standardized = residual / sd;
          const Scalar innovation = standardized -
            ar1_rho * previous_residual->second;
          const Scalar innovation_variance =
            Scalar(1.0) - ar1_rho * ar1_rho;
          result[static_cast<std::size_t>(subject)] +=
            libertad::scalar_log(variance) +
            libertad::scalar_log(innovation_variance) +
            innovation * innovation / innovation_variance;
        } else {
          result[static_cast<std::size_t>(subject)] +=
            libertad::scalar_log(variance) + residual * residual / variance;
        }
        previous_standardized_residual[dvid] = residual / sd;
      }
      }
    }
    if (evid[row] == 0 && std::isfinite(dv[row]) &&
        ((!engine.hmm_enabled && !engine.kalman_enabled) || mdv[row] == 0)) {
      const int sequence =
        (engine.hmm_enabled && !engine.hmm_by_dvid) ||
        (engine.kalman_enabled && !engine.kalman_by_dvid) ? 1 : dvid;
      previous_outcome[sequence] = dv[row];
      previous_outcome_time[sequence] = row_optional(data, "TIME", row, 0.0);
    }
    previous_subject = subject;
  }
  return result;
}

template <class Scalar>
Scalar population_joint_nll_t(const ModelEngine& engine,
                              const Rcpp::DataFrame& data,
                              const std::vector<Scalar>& theta,
                              const std::vector<Scalar>& eta,
                              const std::vector<Scalar>& sigma,
                              const std::vector<Scalar>& omega,
                              bool interaction = true) {
  Rcpp::IntegerVector subjects = data[".ID_INDEX"];
  int n_subjects = 0;
  for (int value : subjects) n_subjects = std::max(n_subjects, value);
  Scalar total = Scalar(0.0);
  if (engine.mixture_probabilities.empty()) {
    const std::vector<Scalar> prediction = simulate_analytical_t(
      engine, data, theta, eta, sigma);
    std::vector<Scalar> variance_prediction;
    if (!interaction && engine.error_type != "likelihood") {
      const std::vector<Scalar> zero_eta(eta.size(), Scalar(0.0));
      variance_prediction = simulate_analytical_t(
        engine, data, theta, zero_eta, sigma);
    }
    const std::vector<Scalar> residual = residual_subject_nll_t(
      engine, data, prediction, theta, eta, sigma, std::vector<int>(),
      interaction || engine.error_type == "likelihood" ? nullptr : &variance_prediction);
    for (const Scalar& value : residual) total += value;
  } else {
    std::vector<std::vector<Scalar>> component_nll;
    component_nll.reserve(engine.mixture_probabilities.size());
    for (std::size_t component = 0; component < engine.mixture_probabilities.size(); ++component) {
      std::vector<int> assignment(static_cast<std::size_t>(n_subjects),
                                  static_cast<int>(component + 1));
      const std::vector<Scalar> prediction = simulate_analytical_t(
        engine, data, theta, eta, sigma, assignment);
      std::vector<Scalar> variance_prediction;
      if (!interaction && engine.error_type != "likelihood") {
        const std::vector<Scalar> zero_eta(eta.size(), Scalar(0.0));
        variance_prediction = simulate_analytical_t(
          engine, data, theta, zero_eta, sigma, assignment);
      }
      component_nll.push_back(residual_subject_nll_t(
        engine, data, prediction, theta, eta, sigma, assignment,
        interaction || engine.error_type == "likelihood" ? nullptr : &variance_prediction));
    }
    for (int subject = 0; subject < n_subjects; ++subject) {
      std::vector<Scalar> log_component(engine.mixture_probabilities.size());
      std::size_t maximum_index = 0;
      for (std::size_t component = 0; component < log_component.size(); ++component) {
        log_component[component] = Scalar(std::log(engine.mixture_probabilities[component])) -
          Scalar(0.5) * component_nll[component][static_cast<std::size_t>(subject)];
        if (scalar_value(log_component[component]) > scalar_value(log_component[maximum_index])) {
          maximum_index = component;
        }
      }
      const Scalar maximum = log_component[maximum_index];
      Scalar sum = Scalar(0.0);
      for (const Scalar& value : log_component) {
        sum += libertad::scalar_exp(value - maximum);
      }
      total -= Scalar(2.0) * (maximum + libertad::scalar_log(sum));
    }
  }

  if (engine.n_eta > 0) {
    const MatrixT<Scalar> covariance = omega_matrix_t(engine, omega);
    int n_subjects = 0;
    for (int value : subjects) n_subjects = std::max(n_subjects, value);
    if (n_subjects < 1 || eta.size() % static_cast<std::size_t>(n_subjects) != 0U) {
      throw std::invalid_argument("ETA vector cannot be divided into subject rows for its prior.");
    }
    const int eta_columns = static_cast<int>(eta.size() / static_cast<std::size_t>(n_subjects));
    for (int subject = 0; subject < n_subjects; ++subject) {
      if (engine.re_enabled) {
        int offset = 0;
        for (std::size_t block_index = 0; block_index < engine.re_blocks.size(); ++block_index) {
          const std::vector<int>& indices = engine.re_blocks[block_index];
          const std::string total_name = ".RE_TOTAL_" + std::to_string(block_index + 1U);
          if (!data.containsElementNamed(total_name.c_str())) {
            throw std::invalid_argument("Random-effect unit totals are missing from compiled data.");
          }
          const int units = static_cast<int>(data_value(data, total_name, 0));
          MatrixT<Scalar> block_covariance(indices.size(), indices.size());
          for (std::size_t row = 0; row < indices.size(); ++row) {
            for (std::size_t column = 0; column < indices.size(); ++column) {
              block_covariance(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) =
                covariance(indices[row], indices[column]);
            }
          }
          for (int unit = 0; unit < units; ++unit) {
            VectorT<Scalar> effect(static_cast<Eigen::Index>(indices.size()));
            for (std::size_t index = 0; index < indices.size(); ++index) {
              const int column = offset + unit * static_cast<int>(indices.size()) +
                static_cast<int>(index);
              if (column >= eta_columns) {
                throw std::invalid_argument("Expanded random-effect columns exceed the ETA matrix.");
              }
              effect[static_cast<Eigen::Index>(index)] =
                eta[static_cast<std::size_t>(subject * eta_columns + column)];
            }
            total += omega_subject_prior_t(block_covariance, effect);
          }
          offset += units * static_cast<int>(indices.size());
        }
        if (offset != eta_columns) {
          throw std::invalid_argument("Expanded random-effect columns do not match the design.");
        }
        continue;
      }
      if (engine.iov <= 0) {
        if (eta_columns != engine.n_eta) {
          throw std::invalid_argument("ETA matrix columns do not match OMEGA dimension.");
        }
        VectorT<Scalar> effect(engine.n_eta);
        for (int index = 0; index < engine.n_eta; ++index) {
          effect[index] = eta[static_cast<std::size_t>(subject * eta_columns + index)];
        }
        total += omega_subject_prior_t(covariance, effect);
        continue;
      }
      const int between = engine.n_eta - engine.iov;
      if (eta_columns < between || (eta_columns - between) % engine.iov != 0) {
        throw std::invalid_argument("IOV ETA columns do not match the occasion layout.");
      }
      if (between > 0) {
        MatrixT<Scalar> between_covariance = covariance.topLeftCorner(between, between);
        VectorT<Scalar> between_effect(between);
        for (int index = 0; index < between; ++index) {
          between_effect[index] = eta[static_cast<std::size_t>(subject * eta_columns + index)];
        }
        total += omega_subject_prior_t(between_covariance, between_effect);
      }
      const MatrixT<Scalar> iov_covariance = covariance.bottomRightCorner(engine.iov, engine.iov);
      const int occasions = (eta_columns - between) / engine.iov;
      for (int occasion = 0; occasion < occasions; ++occasion) {
        VectorT<Scalar> occasion_effect(engine.iov);
        for (int index = 0; index < engine.iov; ++index) {
          const int column = between + occasion * engine.iov + index;
          occasion_effect[index] = eta[static_cast<std::size_t>(subject * eta_columns + column)];
        }
        total += omega_subject_prior_t(iov_covariance, occasion_effect);
      }
    }
  }
  return total;
}

template <class Scalar>
Scalar positive_definite_gaussian_nll_t(
    const MatrixT<Scalar>& covariance, const VectorT<Scalar>& residual,
    const std::string& context) {
  const Eigen::Index dimension = covariance.rows();
  if (covariance.cols() != dimension || residual.size() != dimension) {
    throw std::invalid_argument(context + " dimensions are inconsistent.");
  }
  if (!dimension) return Scalar(0.0);
  MatrixT<Scalar> lower = MatrixT<Scalar>::Zero(dimension, dimension);
  Scalar logdet = Scalar(0.0);
  for (Eigen::Index row = 0; row < dimension; ++row) {
    for (Eigen::Index column = 0; column <= row; ++column) {
      Scalar value = covariance(row, column);
      for (Eigen::Index inner = 0; inner < column; ++inner) {
        value -= lower(row, inner) * lower(column, inner);
      }
      if (row == column) {
        if (!(scalar_value(value) > 1e-14)) {
          throw std::domain_error(context + " is not positive definite at the recording point.");
        }
        lower(row, column) = CppAD::sqrt(value);
        logdet += Scalar(2.0) * CppAD::log(lower(row, column));
      } else {
        lower(row, column) = value / lower(column, column);
      }
    }
  }
  VectorT<Scalar> forward(dimension);
  for (Eigen::Index row = 0; row < dimension; ++row) {
    Scalar value = residual[row];
    for (Eigen::Index column = 0; column < row; ++column) {
      value -= lower(row, column) * forward[column];
    }
    forward[row] = value / lower(row, row);
  }
  VectorT<Scalar> solution(dimension);
  for (Eigen::Index offset = 0; offset < dimension; ++offset) {
    const Eigen::Index row = dimension - offset - 1;
    Scalar value = forward[row];
    for (Eigen::Index column = row + 1; column < dimension; ++column) {
      value -= lower(column, row) * solution[column];
    }
    solution[row] = value / lower(row, row);
  }
  Scalar quadratic = Scalar(0.0);
  for (Eigen::Index row = 0; row < dimension; ++row) {
    quadratic += residual[row] * solution[row];
  }
  return logdet + quadratic;
}

template <class Scalar>
Scalar positive_definite_logdet_t(
    const MatrixT<Scalar>& covariance, const std::string& context) {
  const Eigen::Index dimension = covariance.rows();
  if (covariance.cols() != dimension) {
    throw std::invalid_argument(context + " must be square.");
  }
  MatrixT<Scalar> lower = MatrixT<Scalar>::Zero(dimension, dimension);
  Scalar logdet = Scalar(0.0);
  for (Eigen::Index row = 0; row < dimension; ++row) {
    for (Eigen::Index column = 0; column <= row; ++column) {
      Scalar value = covariance(row, column);
      for (Eigen::Index inner = 0; inner < column; ++inner) {
        value -= lower(row, inner) * lower(column, inner);
      }
      if (row == column) {
        if (!(scalar_value(value) > 1e-14)) {
          throw std::domain_error(context + " is not positive definite at the recording point.");
        }
        lower(row, column) = CppAD::sqrt(value);
        logdet += Scalar(2.0) * CppAD::log(lower(row, column));
      } else {
        lower(row, column) = value / lower(column, column);
      }
    }
  }
  return logdet;
}

std::unique_ptr<ObjectiveTape> record_fo_tape(
    const ModelEngine& engine, PredictionTape& prediction_tape,
    const Rcpp::DataFrame& data, const Rcpp::NumericVector& theta,
    const Rcpp::NumericVector& sigma, const Rcpp::NumericVector& omega) {
  const int n_theta = theta.size();
  const int n_sigma = sigma.size();
  const int n_omega = omega.size();
  const int n_eta = static_cast<int>(prediction_tape.domain_names.size()) -
    n_theta - n_sigma;
  if (n_eta < 0 || n_omega != static_cast<int>(engine.omega_rows.size())) {
    throw std::invalid_argument("FO tape parameter dimensions are inconsistent with the model.");
  }
  Rcpp::NumericVector dv = data["DV"];
  Rcpp::NumericVector dvid = data.containsElementNamed("DVID") ?
    Rcpp::NumericVector(data["DVID"]) : Rcpp::NumericVector(data.nrows(), 1.0);
  const std::vector<int> observed = fo_observed_rows(data);
  std::vector<double> dynamic_values = prediction_tape.dynamic_values;
  dynamic_values.reserve(dynamic_values.size() + observed.size());
  for (int row : observed) dynamic_values.push_back(dv[row]);

  std::vector<double> point;
  point.reserve(static_cast<std::size_t>(n_theta + n_sigma + n_omega));
  for (double value : theta) point.push_back(value);
  for (double value : sigma) point.push_back(value);
  for (double value : omega) point.push_back(value);
  std::vector<CppAD::AD<double>> independent(point.begin(), point.end());
  std::vector<CppAD::AD<double>> dynamic(
    dynamic_values.begin(), dynamic_values.end());
  if (dynamic.empty()) CppAD::Independent(independent);
  else CppAD::Independent(independent, dynamic);
  std::vector<CppAD::AD<double>> theta_ad(
    independent.begin(), independent.begin() + n_theta);
  std::vector<CppAD::AD<double>> sigma_ad(
    independent.begin() + n_theta, independent.begin() + n_theta + n_sigma);
  std::vector<CppAD::AD<double>> omega_ad(
    independent.begin() + n_theta + n_sigma, independent.end());

  std::vector<CppAD::AD<double>> prediction_point;
  prediction_point.reserve(prediction_tape.domain_names.size());
  prediction_point.insert(prediction_point.end(), theta_ad.begin(), theta_ad.end());
  prediction_point.insert(
    prediction_point.end(), static_cast<std::size_t>(n_eta), CppAD::AD<double>(0.0));
  prediction_point.insert(prediction_point.end(), sigma_ad.begin(), sigma_ad.end());
  auto prediction_ad = prediction_tape.fun.base2ad();
  if (!prediction_tape.dynamic_values.empty()) {
    std::vector<CppAD::AD<double>> prediction_dynamic(
      dynamic.begin(),
      dynamic.begin() + static_cast<std::ptrdiff_t>(prediction_tape.dynamic_values.size()));
    prediction_ad.new_dynamic(prediction_dynamic);
  }
  std::ostringstream messages;
  const std::vector<CppAD::AD<double>> predictions =
    prediction_ad.Forward(0, prediction_point, messages);
  MatrixT<CppAD::AD<double>> eta_jacobian(
    static_cast<Eigen::Index>(predictions.size()), n_eta);
  std::vector<CppAD::AD<double>> direction(
    prediction_tape.domain_names.size(), CppAD::AD<double>(0.0));
  for (int eta = 0; eta < n_eta; ++eta) {
    direction[static_cast<std::size_t>(n_theta + eta)] = CppAD::AD<double>(1.0);
    const std::vector<CppAD::AD<double>> derivative =
      prediction_ad.Forward(1, direction, messages);
    direction[static_cast<std::size_t>(n_theta + eta)] = CppAD::AD<double>(0.0);
    for (std::size_t row = 0; row < derivative.size(); ++row) {
      eta_jacobian(static_cast<Eigen::Index>(row), eta) = derivative[row];
    }
  }

  const Eigen::Index n_observed = static_cast<Eigen::Index>(observed.size());
  VectorT<CppAD::AD<double>> residual(n_observed);
  VectorT<CppAD::AD<double>> variance(n_observed);
  MatrixT<CppAD::AD<double>> jacobian(n_observed, n_eta);
  for (Eigen::Index index = 0; index < n_observed; ++index) {
    const int row = observed[static_cast<std::size_t>(index)];
    const CppAD::AD<double> prediction = predictions[static_cast<std::size_t>(row)];
    const CppAD::AD<double> observation = dynamic[
      prediction_tape.dynamic_values.size() + static_cast<std::size_t>(index)];
    variance[index] = residual_variance_t(
      engine, prediction, sigma_ad, std::max(1, static_cast<int>(dvid[row])));
    if (engine.error_type == "exponential") {
      if (!(dv[row] > 0.0) || !(scalar_value(prediction) > 0.0)) {
        throw std::domain_error("FO exponential likelihood requires positive DV and predictions.");
      }
      residual[index] = CppAD::log(observation) - CppAD::log(prediction);
      for (int eta = 0; eta < n_eta; ++eta) {
        jacobian(index, eta) = eta_jacobian(row, eta) / prediction;
      }
    } else {
      residual[index] = observation - prediction;
      for (int eta = 0; eta < n_eta; ++eta) {
        jacobian(index, eta) = eta_jacobian(row, eta);
      }
    }
  }

  MatrixT<CppAD::AD<double>> base_omega =
    MatrixT<CppAD::AD<double>>::Zero(engine.n_eta, engine.n_eta);
  for (int index = 0; index < n_omega; ++index) {
    const int row = engine.omega_rows[static_cast<std::size_t>(index)];
    const int column = engine.omega_cols[static_cast<std::size_t>(index)];
    base_omega(row, column) = omega_ad[static_cast<std::size_t>(index)];
    base_omega(column, row) = omega_ad[static_cast<std::size_t>(index)];
  }
  MatrixT<CppAD::AD<double>> effect_omega = expanded_omega_t(
    engine, data, base_omega, n_eta);
  if (effect_omega.rows() != n_eta) {
    throw std::invalid_argument("FO random-effect covariance has the wrong dimension.");
  }

  MatrixT<CppAD::AD<double>> residual_covariance(n_observed, n_observed);
  for (Eigen::Index row = 0; row < n_observed; ++row) {
    for (Eigen::Index column = 0; column < n_observed; ++column) {
      CppAD::AD<double> correlation = row == column ?
        CppAD::AD<double>(1.0) : CppAD::AD<double>(0.0);
      if (engine.sigma_correlation == "ar1" && dvid[observed[row]] == dvid[observed[column]]) {
        const CppAD::AD<double> rho = ar1_rho_t(engine, theta_ad, sigma_ad);
        const Eigen::Index first = std::min(row, column);
        const Eigen::Index last = std::max(row, column);
        int lag = 0;
        for (Eigen::Index position = first + 1; position <= last; ++position) {
          if (dvid[observed[position]] == dvid[observed[row]]) ++lag;
        }
        correlation = CppAD::pow(rho, lag);
      }
      if (!engine.residual_groups.empty() && row != column &&
          row_optional(data, "TIME", observed[row], 0.0) ==
            row_optional(data, "TIME", observed[column], 0.0)) {
        const int group_index = residual_group_for_dvid(engine, dvid[observed[row]]);
        if (group_index >= 0 && residual_group_for_dvid(engine, dvid[observed[column]]) == group_index) {
          const ResidualGroupSpec& group =
            engine.residual_groups[static_cast<std::size_t>(group_index)];
          correlation = residual_group_correlation_t(
            group, residual_group_endpoint(group, dvid[observed[row]]),
            residual_group_endpoint(group, dvid[observed[column]]),
            theta_ad, sigma_ad);
        }
      }
      residual_covariance(row, column) = correlation *
        CppAD::sqrt(variance[row] * variance[column]);
    }
  }
  MatrixT<CppAD::AD<double>> marginal = residual_covariance +
    jacobian * effect_omega * jacobian.transpose();
  std::vector<CppAD::AD<double>> dependent(1);
  dependent[0] = positive_definite_gaussian_nll_t(
    marginal, residual, "FO marginal covariance");
  auto tape = std::make_unique<ObjectiveTape>();
  tape->fun.Dependent(independent, dependent);
  tape->fun.optimize();
  for (int index = 0; index < n_theta; ++index) {
    tape->domain_names.push_back("THETA_" + std::to_string(index + 1));
  }
  for (int index = 0; index < n_sigma; ++index) {
    tape->domain_names.push_back("SIGMA_" + std::to_string(index + 1));
  }
  for (int index = 0; index < n_omega; ++index) {
    tape->domain_names.push_back("OMEGA_" + std::to_string(index + 1));
  }
  tape->dynamic_columns = prediction_tape.dynamic_columns;
  tape->dynamic_observed_rows = observed;
  tape->structural_dvid = fo_dvid_values(data);
  tape->dynamic_values = dynamic_values;
  tape->n_rows = data.nrows();
  return tape;
}
std::unique_ptr<ObjectiveTape> record_curvature_tape(
    const ModelEngine& engine, PredictionTape& prediction_tape,
    ObjectiveTape& objective_tape, const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta, const Rcpp::NumericVector& eta,
    const Rcpp::NumericVector& sigma, const Rcpp::NumericVector& omega,
    const std::string& approximation) {
  if (approximation != "foce" && approximation != "focei" &&
      approximation != "laplace") {
    throw std::invalid_argument("Unknown conditional-curvature approximation.");
  }
  const int n_theta = theta.size();
  const int n_eta = eta.size();
  const int n_sigma = sigma.size();
  const int n_omega = omega.size();
  std::vector<double> point;
  point.reserve(static_cast<std::size_t>(n_theta + n_eta + n_sigma + n_omega));
  for (double value : theta) point.push_back(value);
  for (double value : eta) point.push_back(value);
  for (double value : sigma) point.push_back(value);
  for (double value : omega) point.push_back(value);
  if (point.size() != objective_tape.domain_names.size() ||
      prediction_tape.domain_names.size() !=
        static_cast<std::size_t>(n_theta + n_eta + n_sigma)) {
    throw std::invalid_argument("Curvature tape parameter dimensions are inconsistent.");
  }
  std::vector<CppAD::AD<double>> independent(point.begin(), point.end());
  CppAD::Independent(independent);
  std::ostringstream messages;
  MatrixT<CppAD::AD<double>> curvature(n_eta, n_eta);

  if (approximation == "laplace") {
    auto objective_ad = objective_tape.fun.base2ad();
    objective_ad.Forward(0, independent, messages);
    std::vector<CppAD::AD<double>> direction(
      independent.size(), CppAD::AD<double>(0.0));
    const std::vector<CppAD::AD<double>> weight(1, CppAD::AD<double>(1.0));
    for (int column = 0; column < n_eta; ++column) {
      const std::size_t position = static_cast<std::size_t>(n_theta + column);
      direction[position] = CppAD::AD<double>(1.0);
      objective_ad.Forward(1, direction, messages);
      direction[position] = CppAD::AD<double>(0.0);
      const std::vector<CppAD::AD<double>> reverse =
        objective_ad.Reverse(2, weight);
      for (int row = 0; row < n_eta; ++row) {
        const std::size_t row_position = static_cast<std::size_t>(n_theta + row);
        curvature(row, column) = reverse[row_position * 2U + 1U];
      }
    }
    curvature = CppAD::AD<double>(0.5) *
      MatrixT<CppAD::AD<double>>(curvature + curvature.transpose());
  } else {
    std::vector<CppAD::AD<double>> prediction_point(
      independent.begin(), independent.begin() + n_theta + n_eta + n_sigma);
    auto prediction_ad = prediction_tape.fun.base2ad();
    const std::vector<CppAD::AD<double>> prediction =
      prediction_ad.Forward(0, prediction_point, messages);
    MatrixT<CppAD::AD<double>> eta_jacobian(
      static_cast<Eigen::Index>(prediction.size()), n_eta);
    std::vector<CppAD::AD<double>> direction(
      prediction_point.size(), CppAD::AD<double>(0.0));
    for (int column = 0; column < n_eta; ++column) {
      direction[static_cast<std::size_t>(n_theta + column)] = CppAD::AD<double>(1.0);
      const std::vector<CppAD::AD<double>> derivative =
        prediction_ad.Forward(1, direction, messages);
      direction[static_cast<std::size_t>(n_theta + column)] = CppAD::AD<double>(0.0);
      for (std::size_t row = 0; row < derivative.size(); ++row) {
        eta_jacobian(static_cast<Eigen::Index>(row), column) = derivative[row];
      }
    }
    std::vector<CppAD::AD<double>> scale_prediction = prediction;
    if (approximation == "foce") {
      std::vector<CppAD::AD<double>> zero_eta_point = prediction_point;
      for (int column = 0; column < n_eta; ++column) {
        zero_eta_point[static_cast<std::size_t>(n_theta + column)] =
          CppAD::AD<double>(0.0);
      }
      scale_prediction = prediction_ad.Forward(0, zero_eta_point, messages);
    }
    std::vector<CppAD::AD<double>> sigma_ad(
      independent.begin() + n_theta + n_eta,
      independent.begin() + n_theta + n_eta + n_sigma);
    Rcpp::NumericVector dv = data["DV"];
    Rcpp::NumericVector evid = data["EVID"];
    Rcpp::NumericVector mdv = data["MDV"];
    Rcpp::NumericVector dvid = data.containsElementNamed("DVID") ?
      Rcpp::NumericVector(data["DVID"]) : Rcpp::NumericVector(data.nrows(), 1.0);
    curvature.setZero();
    for (int row = 0; row < data.nrows(); ++row) {
      if (evid[row] != 0.0 || mdv[row] != 0.0 || !std::isfinite(dv[row])) continue;
      const CppAD::AD<double> variance = residual_variance_t(
        engine, scale_prediction[static_cast<std::size_t>(row)], sigma_ad,
        std::max(1, static_cast<int>(dvid[row])));
      for (int first = 0; first < n_eta; ++first) {
        for (int second = 0; second < n_eta; ++second) {
          curvature(first, second) += CppAD::AD<double>(2.0) *
            eta_jacobian(row, first) * eta_jacobian(row, second) / variance;
        }
      }
    }

    MatrixT<CppAD::AD<double>> base_omega =
      MatrixT<CppAD::AD<double>>::Zero(engine.n_eta, engine.n_eta);
    const std::size_t omega_offset = static_cast<std::size_t>(n_theta + n_eta + n_sigma);
    for (int index = 0; index < n_omega; ++index) {
      const int row = engine.omega_rows[static_cast<std::size_t>(index)];
      const int column = engine.omega_cols[static_cast<std::size_t>(index)];
      const CppAD::AD<double> value = independent[omega_offset + index];
      base_omega(row, column) = value;
      base_omega(column, row) = value;
    }
    MatrixT<CppAD::AD<double>> effect_omega = expanded_omega_t(
      engine, data, base_omega, n_eta);
    MatrixT<CppAD::AD<double>> identity =
      MatrixT<CppAD::AD<double>>::Identity(n_eta, n_eta);
    const MatrixT<CppAD::AD<double>> omega_inverse = solve_linear(
      effect_omega, identity, "Conditional OMEGA curvature");
    curvature += CppAD::AD<double>(2.0) * omega_inverse;
  }

  std::vector<CppAD::AD<double>> dependent(1);
  dependent[0] = positive_definite_logdet_t(
    curvature, "Conditional curvature determinant");
  auto tape = std::make_unique<ObjectiveTape>();
  tape->fun.Dependent(independent, dependent);
  tape->fun.optimize();
  tape->domain_names = objective_tape.domain_names;
  return tape;
}
std::unique_ptr<ObjectiveTape> record_objective_tape(
    const ModelEngine& engine, const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta, const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma, const Rcpp::NumericVector& omega,
    bool interaction) {
  std::vector<double> point = flatten_parameters(theta, eta, sigma);
  for (double value : omega) point.push_back(value);
  std::vector<CppAD::AD<double>> independent(point.begin(), point.end());
  CppAD::Independent(independent);
  std::size_t cursor = 0;
  std::vector<CppAD::AD<double>> theta_ad(static_cast<std::size_t>(theta.size()));
  for (auto& value : theta_ad) value = independent[cursor++];
  std::vector<CppAD::AD<double>> eta_ad(static_cast<std::size_t>(eta.size()));
  for (auto& value : eta_ad) value = independent[cursor++];
  std::vector<CppAD::AD<double>> sigma_ad(static_cast<std::size_t>(sigma.size()));
  for (auto& value : sigma_ad) value = independent[cursor++];
  std::vector<CppAD::AD<double>> omega_ad(static_cast<std::size_t>(omega.size()));
  for (auto& value : omega_ad) value = independent[cursor++];
  std::vector<CppAD::AD<double>> dependent(1);
  dependent[0] = population_joint_nll_t(
    engine, data, theta_ad, eta_ad, sigma_ad, omega_ad, interaction);
  auto tape = std::make_unique<ObjectiveTape>();
  tape->fun.Dependent(independent, dependent);
  tape->fun.optimize();
  tape->domain_names = parameter_names(theta.size(), eta.nrow(), eta.ncol(), sigma.size());
  for (int i = 0; i < omega.size(); ++i) {
    tape->domain_names.push_back("OMEGA_" + std::to_string(i + 1));
  }
  return tape;
}
// End of likelihood implementation.
