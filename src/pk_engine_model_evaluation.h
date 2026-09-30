// Internal implementation fragment included inside namespace liberation by
// pk_engine.cpp. It intentionally has no include guard or namespace declaration.
// Contains ETA mapping, parameter evaluation, and algebraic model evaluation.

int eta_column(const ModelEngine& engine, const EventDataView& data,
               int row, int eta_index, int eta_columns) {
  if (eta_index < 0 || eta_index >= engine.n_eta) {
    throw std::out_of_range("ETA index exceeds the model ETA definitions.");
  }
  if (engine.re_enabled) {
    const std::string column_name = ".ETA_COLUMN_" + std::to_string(eta_index + 1);
    if (!data.containsElementNamed(column_name.c_str())) {
      throw std::invalid_argument("General random-effect execution requires compiled ETA mapping columns.");
    }
    const int column = static_cast<int>(data_value(data, column_name, row)) - 1;
    if (column < 0 || column >= eta_columns) {
      throw std::out_of_range("Mapped random-effect ETA index exceeds the supplied ETA matrix.");
    }
    return column;
  }
  if (engine.iov <= 0 || eta_index < engine.n_eta - engine.iov) return eta_index;
  if (!data.containsElementNamed(".OCC_INDEX")) {
    throw std::invalid_argument("IOV execution requires compiled .OCC_INDEX data.");
  }
  const int between = engine.n_eta - engine.iov;
  const int occasion = static_cast<int>(data_value(data, ".OCC_INDEX", row)) - 1;
  const int column = between + occasion * engine.iov + (eta_index - between);
  if (occasion < 0 || column < 0 || column >= eta_columns) {
    throw std::out_of_range("Occasion-specific ETA index exceeds the supplied ETA matrix.");
  }
  return column;
}

int required_eta_columns(const ModelEngine& engine,
                         const EventDataView& data) {
  if (engine.re_enabled) {
    int maximum = 0;
    for (int eta = 1; eta <= engine.n_eta; ++eta) {
      const std::string name = ".ETA_COLUMN_" + std::to_string(eta);
      if (!data.containsElementNamed(name.c_str())) {
        throw std::invalid_argument("Compiled random-effect ETA mapping is missing.");
      }
      for (int row = 0; row < data.nrows(); ++row) {
        maximum = std::max(maximum,
          static_cast<int>(data_value(data, name, row)));
      }
    }
    return maximum;
  }
  if (engine.iov <= 0) return engine.n_eta;
  if (!data.containsElementNamed(".OCC_INDEX")) {
    throw std::invalid_argument("IOV execution requires compiled .OCC_INDEX data.");
  }
  int count = 0;
  for (int row = 0; row < data.nrows(); ++row) {
    count = std::max(count,
      static_cast<int>(data_value(data, ".OCC_INDEX", row)));
  }
  return engine.n_eta - engine.iov + count * engine.iov;
}

Parameters evaluate_parameters(const ModelEngine& engine,
                               const EventDataView& data,
                               int row, int subject,
                               const Rcpp::NumericVector& theta,
                               const Rcpp::NumericMatrix& eta,
                               const Rcpp::NumericVector& sigma) {
  std::vector<double> inputs(engine.pred->input_names.size(), 0.0);
  for (std::size_t i = 0; i < engine.pred->input_names.size(); ++i) {
    const CompiledInputBinding& binding = engine.pred_inputs[i];
    const std::string& name = binding.name;
    int index = binding.theta;
    if (index >= 0) {
      if (index >= theta.size()) throw std::out_of_range("THETA index exceeds supplied values.");
      inputs[i] = theta[index];
      continue;
    }
    index = binding.eta;
    if (index >= 0) {
      inputs[i] = eta(subject, eta_column(engine, data, row, index, eta.ncol()));
      continue;
    }
    index = binding.sigma;
    if (index >= 0) {
      if (index >= sigma.size()) throw std::out_of_range("SIGMA index exceeds supplied values.");
      inputs[i] = sigma[index];
      continue;
    }
    if (binding.zero_error && name != "F") {
      inputs[i] = 0.0;
      continue;
    }
    if (name == "F") {
      inputs[i] = 0.0;
      continue;
    }
    if (name == "MIXNUM") {
      inputs[i] = data.containsElementNamed("MIXNUM") ?
        data_value(data, "MIXNUM", row) : 1.0;
      continue;
    }
    inputs[i] = data_value(data, name, row);
    if (!std::isfinite(inputs[i])) {
      throw std::domain_error("PRED input '" + name + "' is non-finite at row " +
                              std::to_string(row + 1) + ".");
    }
  }
  std::vector<double> output = engine.pred->eval_outputs(inputs, engine.all_outputs);
  Parameters parameters;
  for (std::size_t i = 0; i < output.size(); ++i) {
    parameters[engine.pred->output_names[i]] = output[i];
  }
  return parameters;
}

double evaluate_post_prediction(
    const ModelEngine& engine, const EventDataView& data,
    int row, int subject, double time, const Vector& state,
    const Rcpp::NumericVector& theta, const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma, double advan_prediction,
    Parameters& parameters) {
  if (!engine.post_pred) return advan_prediction;
  std::vector<double> inputs(engine.post_pred->input_names.size(), 0.0);
  for (std::size_t i = 0; i < engine.post_pred->input_names.size(); ++i) {
    const CompiledInputBinding& binding = engine.post_pred_inputs[i];
    const std::string& name = binding.name;
    int index = binding.theta;
    if (index >= 0) {
      if (index >= theta.size()) throw std::out_of_range("THETA index exceeds values.");
      inputs[i] = theta[index];
      continue;
    }
    index = binding.eta;
    if (index >= 0) {
      inputs[i] = eta(subject, eta_column(engine, data, row, index, eta.ncol()));
      continue;
    }
    index = binding.sigma;
    if (index >= 0) {
      if (index >= sigma.size()) throw std::out_of_range("SIGMA index exceeds values.");
      inputs[i] = sigma[index];
      continue;
    }
    index = binding.state;
    if (index >= 0) {
      if (index >= state.size()) throw std::out_of_range("$PRED A() index exceeds state dimension.");
      inputs[i] = state[index];
      continue;
    }
    if (name == "F_ADVAN") { inputs[i] = advan_prediction; continue; }
    if (name == "T" || name == "TIME") { inputs[i] = time; continue; }
    if (name == "MIXNUM") {
      inputs[i] = data.containsElementNamed("MIXNUM") ?
        data_value(data, "MIXNUM", row) : 1.0;
      continue;
    }
    const auto assigned = parameters.find(name);
    if (assigned != parameters.end()) {
      inputs[i] = assigned->second;
      continue;
    }
    inputs[i] = data_value(data, name, row);
    if (!std::isfinite(inputs[i])) {
      throw std::domain_error(
        "Post-ADVAN $PRED input '" + name + "' is non-finite at row " +
        std::to_string(row + 1) + ".");
    }
  }
  const std::vector<double> output =
    engine.post_pred->eval_outputs(inputs, engine.post_all_outputs);
  for (std::size_t i = 0; i < output.size(); ++i) {
    parameters[engine.post_pred->output_names[i]] = output[i];
  }
  const auto prediction = parameters.find("F");
  if (prediction == parameters.end() || !std::isfinite(prediction->second)) {
    throw std::domain_error("Post-ADVAN $PRED did not produce a finite F.");
  }
  return prediction->second;
}

Vector evaluate_algebraic_residuals(
    const ModelEngine& engine, const EventDataView& data,
    int row, int subject, double time, const Vector& state,
    const Parameters& parameters, const Rcpp::NumericVector& theta,
    const Rcpp::NumericMatrix& eta, const Rcpp::NumericVector& sigma,
    const Vector& algebraic) {
  if (!engine.alg) throw std::logic_error("DAE algebraic residual program is missing.");
  std::vector<double> inputs(engine.alg->input_names.size(), 0.0);
  for (std::size_t i = 0; i < engine.alg->input_names.size(); ++i) {
    const CompiledInputBinding& binding = engine.alg_inputs[i];
    const std::string& name = binding.name;
    int index = binding.state;
    if (index >= 0) {
      if (index >= state.size()) throw std::out_of_range("A() index exceeds the DAE state dimension.");
      inputs[i] = state[index];
      continue;
    }
    if (name == "T") { inputs[i] = time; continue; }
    if (binding.algebraic >= 0) {
      inputs[i] = algebraic[binding.algebraic];
      continue;
    }
    auto parameter = parameters.find(name);
    if (parameter != parameters.end()) { inputs[i] = parameter->second; continue; }
    index = binding.theta;
    if (index >= 0) { inputs[i] = theta[index]; continue; }
    index = binding.eta;
    if (index >= 0) {
      inputs[i] = eta(subject, eta_column(engine, data, row, index, eta.ncol()));
      continue;
    }
    index = binding.sigma;
    if (index >= 0) { inputs[i] = sigma[index]; continue; }
    if (binding.zero_error) continue;
    inputs[i] = data_value(data, name, row);
  }
  const std::vector<double> output = engine.alg->eval_outputs(inputs, engine.algebraic_outputs);
  Vector residual(static_cast<Eigen::Index>(output.size()));
  for (std::size_t i = 0; i < output.size(); ++i) residual[static_cast<Eigen::Index>(i)] = output[i];
  return residual;
}

Vector solve_algebraic(
    const ModelEngine& engine, const EventDataView& data,
    int row, int subject, double time, const Vector& state,
    const Parameters& parameters, const Rcpp::NumericVector& theta,
    const Rcpp::NumericMatrix& eta, const Rcpp::NumericVector& sigma) {
  Vector value(static_cast<Eigen::Index>(engine.dae_initial.size()));
  for (std::size_t i = 0; i < engine.dae_initial.size(); ++i) value[static_cast<Eigen::Index>(i)] = engine.dae_initial[i];
  Vector residual;
  for (int iteration = 0; iteration < engine.dae_maxit; ++iteration) {
    residual = evaluate_algebraic_residuals(
      engine, data, row, subject, time, state, parameters, theta, eta, sigma, value);
    if (!residual.allFinite()) throw std::domain_error("DAE residual is non-finite.");
    if (residual.cwiseAbs().maxCoeff() <= engine.dae_tolerance) return value;
    Matrix jacobian = Matrix::Zero(value.size(), value.size());
    for (Eigen::Index column = 0; column < value.size(); ++column) {
      const double delta = engine.dae_jacobian_step * std::max(1.0, std::abs(value[column]));
      Vector plus = value; plus[column] += delta;
      Vector minus = value; minus[column] -= delta;
      const Vector upper = evaluate_algebraic_residuals(
        engine, data, row, subject, time, state, parameters, theta, eta, sigma, plus);
      const Vector lower = evaluate_algebraic_residuals(
        engine, data, row, subject, time, state, parameters, theta, eta, sigma, minus);
      jacobian.col(column) = (upper - lower) / (2.0 * delta);
      if (!engine.dae_sparsity.empty()) {
        for (Eigen::Index row_index = 0; row_index < value.size(); ++row_index) {
          if (!engine.dae_sparsity[static_cast<std::size_t>(row_index * value.size() + column)]) {
            jacobian(row_index, column) = 0.0;
          }
        }
      }
    }
    Vector update = Vector::Zero(value.size());
    for (std::size_t block = 0; block < engine.dae_block_rows.size(); ++block) {
      const auto& rows = engine.dae_block_rows[block];
      const auto& columns = engine.dae_block_columns[block];
      Matrix local(rows.size(), columns.size());
      Vector rhs(rows.size());
      for (std::size_t local_row = 0; local_row < rows.size(); ++local_row) {
        rhs[static_cast<Eigen::Index>(local_row)] = -residual[rows[local_row]];
        for (std::size_t local_column = 0; local_column < columns.size(); ++local_column) {
          local(static_cast<Eigen::Index>(local_row),
                static_cast<Eigen::Index>(local_column)) =
            jacobian(rows[local_row], columns[local_column]);
        }
      }
      Eigen::FullPivLU<Matrix> lu(local);
      if (!lu.isInvertible()) throw std::runtime_error("DAE Newton Jacobian block is singular.");
      const Vector local_update = lu.solve(rhs);
      for (std::size_t local_column = 0; local_column < columns.size(); ++local_column) {
        update[columns[local_column]] = local_update[static_cast<Eigen::Index>(local_column)];
      }
    }
    if (!update.allFinite()) throw std::runtime_error("DAE Newton update is non-finite.");
    value += update;
    if (update.cwiseAbs().maxCoeff() <= engine.dae_tolerance) {
      residual = evaluate_algebraic_residuals(
        engine, data, row, subject, time, state, parameters, theta, eta, sigma, value);
      if (residual.cwiseAbs().maxCoeff() <= 10.0 * engine.dae_tolerance) return value;
    }
  }
  throw std::runtime_error("DAE algebraic Newton solve did not converge.");
}
