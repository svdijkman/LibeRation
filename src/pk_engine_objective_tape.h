// Internal implementation fragment included inside namespace liberation by
// pk_engine.cpp. It intentionally has no include guard or namespace declaration.
// Contains objective/prediction tape types and shared tape plumbing.

struct ObjectiveTape {
  struct EtaOptimizerState {
    std::vector<double> dynamic_values;
    Matrix inverse_hessian;
  };

  CppAD::ADFun<double> fun;
  std::vector<std::string> domain_names;
  std::vector<std::string> dynamic_columns;
  std::vector<int> dynamic_observed_rows;
  std::vector<int> structural_dvid;
  std::vector<double> dynamic_values;
  int n_rows = 0;
  bool fo_low_rank = false;
  bool fo_low_rank_fallback = false;
  std::string fo_low_rank_reason;
  double fo_low_rank_relative_difference = 0.0;
  libertad::SparseHessianCache hessian_cache;
  std::unordered_map<std::size_t, std::vector<EtaOptimizerState>>
    eta_optimizer_states;
  std::size_t eta_optimizer_state_hits = 0U;
  std::size_t eta_optimizer_state_updates = 0U;
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
    const ModelEngine& engine, const EventDataView& data) {
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
    const std::vector<std::string>& columns, const EventDataView& data,
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

std::vector<int> fo_observed_rows(const EventDataView& data) {
  auto dv = data.values("DV");
  auto evid = data.values("EVID");
  auto mdv = data.values("MDV");
  std::vector<int> observed;
  for (int row = 0; row < data.nrows(); ++row) {
    if (evid[row] == 0.0 && mdv[row] == 0.0 && std::isfinite(dv[row])) {
      observed.push_back(row);
    }
  }
  return observed;
}

std::vector<int> fo_dvid_values(const EventDataView& data) {
  std::vector<int> result(static_cast<std::size_t>(data.nrows()), 1);
  if (!data.containsElementNamed("DVID")) return result;
  auto dvid = data.values("DVID");
  for (int row = 0; row < data.nrows(); ++row) {
    result[static_cast<std::size_t>(row)] =
      std::max(1, static_cast<int>(dvid[row]));
  }
  return result;
}

std::vector<double> fo_dynamic_values(const ObjectiveTape& tape,
                                      const EventDataView& data) {
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
  auto dv = data.values("DV");
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

template <class Tape>
void set_tape_dynamic_values(Tape& tape, const std::vector<double>& values,
                             const std::string& context);

void set_fo_dynamic(ObjectiveTape& tape, const EventDataView& data) {
  const std::vector<double> values = fo_dynamic_values(tape, data);
  set_tape_dynamic_values(tape, values, "Shared FO tape");
}

template <class Tape>
void set_tape_dynamic_values(Tape& tape, const std::vector<double>& values,
                             const std::string& context) {
  if (values.size() != tape.fun.size_dyn_ind()) {
    throw std::invalid_argument(
      context + " dynamic-data length does not match the recorded tape.");
  }
  // new_dynamic() invalidates CppAD's current Taylor state.  Dedicated
  // subject tapes commonly receive the same covariates/observations for many
  // consecutive objective evaluations, so replaying an identical dynamic
  // vector only adds work and prevents later Forward/Reverse reuse.
  if (tape.dynamic_values == values) return;
  if (!values.empty()) tape.fun.new_dynamic(values);
  tape.dynamic_values = values;
}

std::vector<double> shared_objective_dynamic_values(
    const ObjectiveTape& tape, const EventDataView& data) {
  if (tape.n_rows != data.nrows()) {
    throw std::invalid_argument(
      "A shared conditional-objective tape received a different number of rows.");
  }
  if (fo_observed_rows(data) != tape.dynamic_observed_rows) {
    throw std::invalid_argument(
      "A shared conditional-objective tape received a different observation pattern.");
  }
  if (fo_dvid_values(data) != tape.structural_dvid) {
    throw std::invalid_argument(
      "A shared conditional-objective tape received a different DVID pattern.");
  }
  std::vector<double> values = prediction_dynamic_values(
    tape.dynamic_columns, data, tape.n_rows);
  auto dv = data.values("DV");
  values.reserve(values.size() + tape.dynamic_observed_rows.size());
  for (int row : tape.dynamic_observed_rows) {
    const double value = dv[row];
    if (!std::isfinite(value)) {
      throw std::domain_error(
        "A shared conditional-objective tape received a non-finite observation.");
    }
    values.push_back(value);
  }
  return values;
}

void set_shared_objective_dynamic(ObjectiveTape& tape,
                                  const EventDataView& data) {
  if (tape.dynamic_values.empty() && tape.dynamic_columns.empty() &&
      tape.dynamic_observed_rows.empty()) return;
  const std::vector<double> values = shared_objective_dynamic_values(tape, data);
  set_tape_dynamic_values(
    tape, values, "Shared conditional-objective tape");
}

inline std::vector<double> objective_dynamic_values(
    const ObjectiveTape& tape, SEXP input) {
  if (tape.fun.size_dyn_ind() == 0U && tape.dynamic_columns.empty() &&
      tape.dynamic_observed_rows.empty()) {
    return std::vector<double>();
  }
  if (Rf_isNumeric(input) && !Rf_inherits(input, "data.frame")) {
    return Rcpp::as<std::vector<double>>(input);
  }
  return shared_objective_dynamic_values(tape, event_data_view(input));
}

inline void set_objective_dynamic_input(ObjectiveTape& tape, SEXP input) {
  set_tape_dynamic_values(
    tape, objective_dynamic_values(tape, input), "Objective tape");
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
    const ModelEngine& engine, const EventDataView& data,
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
