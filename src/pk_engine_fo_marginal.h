// Internal implementation fragment included inside namespace liberation by
// pk_engine.cpp. It intentionally has no include guard or namespace declaration.
// Contains FO dense and low-rank Woodbury marginal likelihood recorders.

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

// Evaluate R + G OMEGA G' without forming or factoring the observation-sized
// marginal covariance.  With independent residual errors R is diagonal and
// the determinant lemma/Woodbury identity reduce the factorisation to the
// random-effect dimension:
//
//   U = G chol(OMEGA), M = I + U' R^-1 U
//   log|R + U U'| = log|R| + log|M|
//   e'(R + U U')^-1 e = e'R^-1e - b'M^-1b.
//
// The calculation is algebraically identical to the dense FO likelihood but
// has a different floating-point evaluation order.  Selection is guarded by
// covariance conditioning and an explicit comparison with the dense result at
// tape creation, so it is safe to use under either numerical policy.
template <class Scalar>
Scalar fo_low_rank_gaussian_nll_t(
    const VectorT<Scalar>& variance, const MatrixT<Scalar>& jacobian,
    const MatrixT<Scalar>& omega, const VectorT<Scalar>& residual,
    const std::string& context) {
  const Eigen::Index observations = residual.size();
  const Eigen::Index effects = omega.rows();
  if (variance.size() != observations || jacobian.rows() != observations ||
      jacobian.cols() != effects || omega.cols() != effects) {
    throw std::invalid_argument(context + " dimensions are inconsistent.");
  }
  if (!observations) return Scalar(0.0);

  MatrixT<Scalar> omega_lower = MatrixT<Scalar>::Zero(effects, effects);
  for (Eigen::Index row = 0; row < effects; ++row) {
    for (Eigen::Index column = 0; column <= row; ++column) {
      Scalar value = omega(row, column);
      for (Eigen::Index inner = 0; inner < column; ++inner) {
        value -= omega_lower(row, inner) * omega_lower(column, inner);
      }
      if (row == column) {
        if (!(scalar_value(value) > 1e-14)) {
          throw std::domain_error(
            context + " random-effect covariance is not positive definite "
            "at the recording point.");
        }
        omega_lower(row, column) = CppAD::sqrt(value);
      } else {
        omega_lower(row, column) = value / omega_lower(column, column);
      }
    }
  }

  const MatrixT<Scalar> low_rank = jacobian * omega_lower;
  MatrixT<Scalar> information = MatrixT<Scalar>::Identity(effects, effects);
  VectorT<Scalar> right = VectorT<Scalar>::Zero(effects);
  Scalar residual_logdet = Scalar(0.0);
  Scalar base_quadratic = Scalar(0.0);
  for (Eigen::Index observation = 0; observation < observations; ++observation) {
    if (!(scalar_value(variance[observation]) > 1e-14)) {
      throw std::domain_error(
        context + " residual variance is not positive at the recording point.");
    }
    const Scalar inverse_variance = Scalar(1.0) / variance[observation];
    residual_logdet += CppAD::log(variance[observation]);
    base_quadratic += residual[observation] * residual[observation] *
      inverse_variance;
    for (Eigen::Index row = 0; row < effects; ++row) {
      right[row] += low_rank(observation, row) * residual[observation] *
        inverse_variance;
      for (Eigen::Index column = 0; column <= row; ++column) {
        information(row, column) += low_rank(observation, row) *
          low_rank(observation, column) * inverse_variance;
      }
    }
  }
  for (Eigen::Index row = 0; row < effects; ++row) {
    for (Eigen::Index column = 0; column < row; ++column) {
      information(column, row) = information(row, column);
    }
  }

  MatrixT<Scalar> lower = MatrixT<Scalar>::Zero(effects, effects);
  Scalar information_logdet = Scalar(0.0);
  for (Eigen::Index row = 0; row < effects; ++row) {
    for (Eigen::Index column = 0; column <= row; ++column) {
      Scalar value = information(row, column);
      for (Eigen::Index inner = 0; inner < column; ++inner) {
        value -= lower(row, inner) * lower(column, inner);
      }
      if (row == column) {
        if (!(scalar_value(value) > 1e-14)) {
          throw std::domain_error(
            context + " low-rank information matrix is not positive definite "
            "at the recording point.");
        }
        lower(row, column) = CppAD::sqrt(value);
        information_logdet += Scalar(2.0) * CppAD::log(lower(row, column));
      } else {
        lower(row, column) = value / lower(column, column);
      }
    }
  }
  VectorT<Scalar> forward(effects);
  for (Eigen::Index row = 0; row < effects; ++row) {
    Scalar value = right[row];
    for (Eigen::Index column = 0; column < row; ++column) {
      value -= lower(row, column) * forward[column];
    }
    forward[row] = value / lower(row, row);
  }
  Scalar correction = Scalar(0.0);
  for (Eigen::Index row = 0; row < effects; ++row) {
    correction += forward[row] * forward[row];
  }
  return residual_logdet + information_logdet + base_quadratic - correction;
}

template <class Scalar>
bool fo_low_rank_conditioned(
    const VectorT<Scalar>& variance, const MatrixT<Scalar>& jacobian,
    const MatrixT<Scalar>& omega, double tolerance, std::string& reason) {
  const Eigen::Index observations = variance.size();
  const Eigen::Index effects = omega.rows();
  if (omega.cols() != effects || jacobian.rows() != observations ||
      jacobian.cols() != effects) {
    reason = "inconsistent low-rank dimensions";
    return false;
  }
  Matrix omega_value(effects, effects);
  Matrix jacobian_value(observations, effects);
  Vector variance_value(observations);
  for (Eigen::Index row = 0; row < observations; ++row) {
    variance_value[row] = scalar_value(variance[row]);
    if (!(variance_value[row] > 1e-14) ||
        !std::isfinite(variance_value[row])) {
      reason = "non-positive or non-finite residual variance";
      return false;
    }
    for (Eigen::Index column = 0; column < effects; ++column) {
      jacobian_value(row, column) = scalar_value(jacobian(row, column));
      if (!std::isfinite(jacobian_value(row, column))) {
        reason = "non-finite ETA Jacobian";
        return false;
      }
    }
  }
  for (Eigen::Index row = 0; row < effects; ++row) {
    for (Eigen::Index column = 0; column < effects; ++column) {
      omega_value(row, column) = scalar_value(omega(row, column));
      if (!std::isfinite(omega_value(row, column))) {
        reason = "non-finite random-effect covariance";
        return false;
      }
    }
  }
  const auto omega_eigen = libertad::detail::self_adjoint_eigen(
    omega_value, false);
  if (omega_eigen.info != Eigen::Success || !omega_eigen.values.size()) {
    reason = "random-effect covariance eigendecomposition failed";
    return false;
  }
  const double omega_max = omega_eigen.values.maxCoeff();
  const double omega_min = omega_eigen.values.minCoeff();
  if (!(omega_max > 0.0) || !(omega_min > std::max(1e-14, tolerance * omega_max))) {
    reason = "random-effect covariance is singular or ill-conditioned";
    return false;
  }
  Eigen::LLT<Matrix> omega_llt(omega_value);
  if (omega_llt.info() != Eigen::Success) {
    reason = "random-effect covariance Cholesky factorisation failed";
    return false;
  }
  const Matrix low_rank = jacobian_value * Matrix(omega_llt.matrixL());
  Matrix information = Matrix::Identity(effects, effects);
  for (Eigen::Index row = 0; row < observations; ++row) {
    information.noalias() += low_rank.row(row).transpose() *
      low_rank.row(row) / variance_value[row];
  }
  const auto information_eigen = libertad::detail::self_adjoint_eigen(
    information, false);
  if (information_eigen.info != Eigen::Success ||
      !information_eigen.values.size()) {
    reason = "low-rank information eigendecomposition failed";
    return false;
  }
  const double information_max = information_eigen.values.maxCoeff();
  const double information_min = information_eigen.values.minCoeff();
  if (!(information_min > std::max(1e-14, tolerance * information_max))) {
    reason = "low-rank information matrix is ill-conditioned";
    return false;
  }
  return true;
}

std::unique_ptr<ObjectiveTape> record_fo_tape(
    const ModelEngine& engine, PredictionTape& prediction_tape,
    const EventDataView& data, const Rcpp::NumericVector& theta,
    const Rcpp::NumericVector& sigma, const Rcpp::NumericVector& omega,
    bool low_rank = false, double low_rank_tolerance = 1e-9,
    double low_rank_condition_tolerance = 1e-12) {
  const int n_theta = theta.size();
  const int n_sigma = sigma.size();
  const int n_omega = omega.size();
  const int n_eta = static_cast<int>(prediction_tape.domain_names.size()) -
    n_theta - n_sigma;
  if (n_eta < 0 || n_omega != static_cast<int>(engine.omega_rows.size())) {
    throw std::invalid_argument("FO tape parameter dimensions are inconsistent with the model.");
  }
  auto dv = data.values("DV");
  auto dvid = data.values("DVID", 1.0);
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

  if (!(low_rank_tolerance >= 0.0) || !std::isfinite(low_rank_tolerance) ||
      !(low_rank_condition_tolerance > 0.0) ||
      !std::isfinite(low_rank_condition_tolerance)) {
    throw std::invalid_argument("FO low-rank tolerances are invalid.");
  }
  auto dense_nll = [&]() {
    MatrixT<CppAD::AD<double>> residual_covariance(n_observed, n_observed);
    for (Eigen::Index row = 0; row < n_observed; ++row) {
      for (Eigen::Index column = 0; column < n_observed; ++column) {
        CppAD::AD<double> correlation = row == column ?
          CppAD::AD<double>(1.0) : CppAD::AD<double>(0.0);
        if (engine.sigma_correlation == "ar1" &&
            dvid[observed[row]] == dvid[observed[column]]) {
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
          const int group_index = residual_group_for_dvid(
            engine, dvid[observed[row]]);
          if (group_index >= 0 && residual_group_for_dvid(
                engine, dvid[observed[column]]) == group_index) {
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
    return positive_definite_gaussian_nll_t(
      marginal, residual, "FO marginal covariance");
  };

  std::vector<CppAD::AD<double>> dependent(1);
  const bool candidate = low_rank && n_eta > 0 && n_eta < n_observed &&
    engine.sigma_correlation == "independent" &&
    engine.residual_groups.empty();
  bool use_low_rank = false;
  bool low_rank_fallback = false;
  std::string low_rank_reason;
  double low_rank_relative_difference = 0.0;
  if (candidate && fo_low_rank_conditioned(
        variance, jacobian, effect_omega, low_rank_condition_tolerance,
        low_rank_reason)) {
    try {
      const CppAD::AD<double> low_rank_nll = fo_low_rank_gaussian_nll_t(
        variance, jacobian, effect_omega, residual,
        "FO marginal covariance");
      const CppAD::AD<double> dense = dense_nll();
      const double low_rank_value = scalar_value(low_rank_nll);
      const double dense_value = scalar_value(dense);
      low_rank_relative_difference = std::abs(low_rank_value - dense_value) /
        std::max(1.0, std::abs(dense_value));
      if (std::isfinite(low_rank_relative_difference) &&
          low_rank_relative_difference <= low_rank_tolerance) {
        dependent[0] = low_rank_nll;
        use_low_rank = true;
        low_rank_reason = "conditioned and dense-equivalent";
      } else {
        dependent[0] = dense;
        low_rank_fallback = true;
        low_rank_reason = "dense-equivalence tolerance exceeded";
      }
    } catch (const std::exception& error) {
      dependent[0] = dense_nll();
      low_rank_fallback = true;
      low_rank_reason = std::string("low-rank construction failed: ") +
        error.what();
    }
  } else {
    dependent[0] = dense_nll();
    low_rank_fallback = candidate;
    if (!candidate && low_rank) {
      low_rank_reason = "model structure requires dense covariance";
    } else if (!low_rank) {
      low_rank_reason = "low-rank route disabled";
    }
  }
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
  tape->fo_low_rank = use_low_rank;
  tape->fo_low_rank_fallback = low_rank_fallback;
  tape->fo_low_rank_reason = low_rank_reason;
  tape->fo_low_rank_relative_difference = low_rank_relative_difference;
  return tape;
}
