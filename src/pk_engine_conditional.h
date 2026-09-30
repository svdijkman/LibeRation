// Internal implementation fragment included inside namespace liberation by
// pk_engine.cpp. It intentionally has no include guard or namespace declaration.
// Contains FOCE, FOCEI, and Laplace conditional objective/curvature recorders.

std::unique_ptr<ObjectiveTape> record_curvature_tape(
    const ModelEngine& engine, PredictionTape& prediction_tape,
    ObjectiveTape& objective_tape, const EventDataView& data,
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
    auto dv = data.values("DV");
    auto evid = data.values("EVID");
    auto mdv = data.values("MDV");
    auto dvid = data.values("DVID", 1.0);
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
    const ModelEngine& engine, const EventDataView& data,
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

// Record one analytical Gaussian conditional objective with observations and
// model covariates as CppAD dynamic parameters.  Subjects with the same event
// and observation design can therefore share the operation sequence while
// retaining their own data.  Eligibility is deliberately narrow: complex
// likelihood, mixture, IOV, and hierarchical random-effect paths continue to
// use the general per-subject recorder above.
std::unique_ptr<ObjectiveTape> record_shared_fo_objective_tape(
    const ModelEngine& engine, PredictionTape& prediction_tape,
    const EventDataView& data, const Rcpp::NumericVector& theta,
    const Rcpp::NumericMatrix& eta, const Rcpp::NumericVector& sigma,
    const Rcpp::NumericVector& omega) {
  if (engine.error_type == "likelihood" || !engine.residual_groups.empty() ||
      engine.sigma_correlation != "independent" ||
      engine.blq_method != "none" || !engine.mixture_probabilities.empty() ||
      engine.iov != 0 || engine.re_enabled || eta.nrow() != 1 ||
      eta.ncol() != engine.n_eta) {
    throw std::invalid_argument(
      "This model is not eligible for a structurally shared FO conditional tape.");
  }
  if (prediction_tape.domain_names.size() !=
      static_cast<std::size_t>(theta.size() + eta.size() + sigma.size())) {
    throw std::invalid_argument(
      "Shared FO prediction and conditional-objective dimensions differ.");
  }

  std::vector<double> point = flatten_parameters(theta, eta, sigma);
  for (double value : omega) point.push_back(value);
  std::vector<CppAD::AD<double>> independent(point.begin(), point.end());
  const std::vector<int> observed = fo_observed_rows(data);
  std::vector<double> dynamic_values = prediction_dynamic_values(
    prediction_tape.dynamic_columns, data, data.nrows());
  auto dv = data.values("DV");
  dynamic_values.reserve(dynamic_values.size() + observed.size());
  for (int row : observed) dynamic_values.push_back(dv[row]);
  std::vector<CppAD::AD<double>> dynamic(
    dynamic_values.begin(), dynamic_values.end());
  if (dynamic.empty()) CppAD::Independent(independent);
  else CppAD::Independent(independent, dynamic);

  std::size_t cursor = 0U;
  std::vector<CppAD::AD<double>> theta_ad(
    static_cast<std::size_t>(theta.size()));
  for (auto& value : theta_ad) value = independent[cursor++];
  std::vector<CppAD::AD<double>> eta_ad(
    static_cast<std::size_t>(eta.size()));
  for (auto& value : eta_ad) value = independent[cursor++];
  std::vector<CppAD::AD<double>> sigma_ad(
    static_cast<std::size_t>(sigma.size()));
  for (auto& value : sigma_ad) value = independent[cursor++];
  std::vector<CppAD::AD<double>> omega_ad(
    static_cast<std::size_t>(omega.size()));
  for (auto& value : omega_ad) value = independent[cursor++];

  auto prediction_ad = prediction_tape.fun.base2ad();
  const std::size_t prediction_dynamic =
    prediction_tape.dynamic_values.size();
  if (prediction_dynamic) {
    std::vector<CppAD::AD<double>> values(
      dynamic.begin(), dynamic.begin() +
        static_cast<std::ptrdiff_t>(prediction_dynamic));
    prediction_ad.new_dynamic(values);
  }
  std::vector<CppAD::AD<double>> prediction_point;
  prediction_point.reserve(
    static_cast<std::size_t>(theta.size() + eta.size() + sigma.size()));
  prediction_point.insert(
    prediction_point.end(), theta_ad.begin(), theta_ad.end());
  prediction_point.insert(
    prediction_point.end(), eta_ad.begin(), eta_ad.end());
  prediction_point.insert(
    prediction_point.end(), sigma_ad.begin(), sigma_ad.end());
  std::ostringstream messages;
  const std::vector<CppAD::AD<double>> prediction =
    prediction_ad.Forward(0, prediction_point, messages);

  const std::vector<int> dvid = fo_dvid_values(data);
  CppAD::AD<double> objective = CppAD::AD<double>(0.0);
  std::size_t observation = prediction_dynamic;
  for (int row : observed) {
    const CppAD::AD<double> outcome = dynamic[observation++];
    const CppAD::AD<double> fitted = prediction[static_cast<std::size_t>(row)];
    const CppAD::AD<double> variance = residual_variance_t(
      engine, fitted, sigma_ad, dvid[static_cast<std::size_t>(row)]);
    CppAD::AD<double> residual = outcome - fitted;
    if (engine.error_type == "exponential") {
      residual = CppAD::log(outcome) -
        CppAD::log(scalar_floor_t(fitted, 1e-300));
    }
    objective += CppAD::log(variance) + residual * residual / variance;
  }
  if (engine.n_eta > 0) {
    const MatrixT<CppAD::AD<double>> covariance =
      omega_matrix_t(engine, omega_ad);
    VectorT<CppAD::AD<double>> effect(engine.n_eta);
    for (int index = 0; index < engine.n_eta; ++index) {
      effect[index] = eta_ad[static_cast<std::size_t>(index)];
    }
    objective += omega_subject_prior_t(covariance, effect);
  }

  std::vector<CppAD::AD<double>> dependent(1U, objective);
  auto tape = std::make_unique<ObjectiveTape>();
  tape->fun.Dependent(independent, dependent);
  tape->fun.optimize();
  tape->domain_names = parameter_names(
    theta.size(), eta.nrow(), eta.ncol(), sigma.size());
  for (int index = 0; index < omega.size(); ++index) {
    tape->domain_names.push_back("OMEGA_" + std::to_string(index + 1));
  }
  tape->dynamic_columns = prediction_tape.dynamic_columns;
  tape->dynamic_observed_rows = observed;
  tape->structural_dvid = dvid;
  tape->dynamic_values = dynamic_values;
  tape->n_rows = data.nrows();
  return tape;
}
// End of likelihood implementation.
