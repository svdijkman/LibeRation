// Internal implementation fragment included inside namespace liberation by
// pk_engine.cpp. It intentionally has no include guard or namespace declaration.
// Contains the native adaptive and fixed Gaussian-quadrature coordinator.

class NativeGqCoordinator {
 public:
  NativeGqCoordinator(
      StochasticEtaCollection* context, SEXP retained_context,
      const Rcpp::List& map_config, const Rcpp::NumericMatrix& nodes,
      const Rcpp::NumericVector& log_measure,
      const Rcpp::NumericVector& measure_sign, bool adaptive,
      int eta_maxit, double tolerance)
      : context_(context), retained_context_(retained_context), map_(map_config),
        adaptive_(adaptive), eta_maxit_(eta_maxit), tolerance_(tolerance) {
    if (!context_ || eta_maxit_ < 1 || !std::isfinite(tolerance_) ||
        tolerance_ <= 0.0 || nodes.nrow() < 1 ||
        nodes.ncol() != context_->eta_dimension() ||
        log_measure.size() != nodes.nrow() ||
        measure_sign.size() != nodes.nrow()) {
      throw std::invalid_argument("Native GQ coordinator inputs are inconsistent.");
    }
    nodes_.resize(nodes.nrow(), nodes.ncol());
    for (int row = 0; row < nodes.nrow(); ++row) {
      for (int column = 0; column < nodes.ncol(); ++column) {
        nodes_(row, column) = nodes(row, column);
      }
    }
    log_measure_ = Eigen::Map<const Vector>(
      log_measure.begin(), static_cast<Eigen::Index>(log_measure.size()));
    measure_sign_ = Eigen::Map<const Vector>(
      measure_sign.begin(), static_cast<Eigen::Index>(measure_sign.size()));
    if (!nodes_.allFinite() || !log_measure_.allFinite() ||
        !measure_sign_.allFinite()) {
      throw std::invalid_argument("Native GQ nodes and weights must be finite.");
    }
    R_PreserveObject(retained_context_);
  }

  ~NativeGqCoordinator() {
    if (retained_context_ != R_NilValue) R_ReleaseObject(retained_context_);
  }

  NativeGqCoordinator(const NativeGqCoordinator&) = delete;
  NativeGqCoordinator& operator=(const NativeGqCoordinator&) = delete;

  Rcpp::List evaluate(const Rcpp::NumericVector& encoded, bool gradient) {
    Vector point = Eigen::Map<const Vector>(
      encoded.begin(), static_cast<Eigen::Index>(encoded.size()));
    const NativeGqEvaluation& current = evaluate_native(point, gradient);
    return evaluation_list(current, point, gradient);
  }

  Rcpp::List optimize(int maxit, int trace, bool exact_refinement) {
    if (maxit < 1) {
      throw std::invalid_argument("Native GQ optimization requires maxit >= 1.");
    }
    const Rcpp::NumericVector start = Rcpp::wrap(map_.start());
    const Rcpp::NumericVector lower = Rcpp::wrap(map_.lower());
    const Rcpp::NumericVector upper = Rcpp::wrap(map_.upper());
    const auto score_value = [this](const Vector& point) {
      return evaluate_native(point, true).value;
    };
    const auto score_gradient = [this](const Vector& point) {
      return evaluate_native(point, true).native_gradient;
    };
    Rcpp::List score;
    if (!map_.dimension()) {
      const NativeGqEvaluation& fixed = evaluate_native(Vector(0), true);
      score = fixed_optimizer(fixed.value);
    } else {
      score = native_optimizer_core(
        score_value, score_gradient, start, lower, upper,
        maxit, tolerance_, trace);
    }
    score["backend"] = "native-bfgs-gq-score";
    score["coordinator"] = "persistent-native-cpp-gq";
    score["objective_backend"] =
      "persistent-cpp-gaussian-quadrature-objective";

    Rcpp::List selected = Rcpp::clone(score);
    bool refinement_accepted = false;
    if (exact_refinement && map_.dimension()) {
      const Rcpp::NumericVector refinement_start = score["par"];
      const auto exact_value = [this](const Vector& point) {
        return evaluate_native(point, false).value;
      };
      const auto finite_gradient = [this](const Vector& point) {
        return finite_difference(point);
      };
      Rcpp::List refined = native_optimizer_core(
        exact_value, finite_gradient, refinement_start, lower, upper,
        maxit, tolerance_, trace);
      const double score_value_at_end = Rcpp::as<double>(score["value"]);
      const double refined_value = Rcpp::as<double>(refined["value"]);
      if (std::isfinite(refined_value) && refined_value <=
          score_value_at_end + tolerance_ *
            std::max(1.0, std::abs(score_value_at_end))) {
        refined["score_search"] = Rcpp::List::create(
          Rcpp::Named("value") = score["value"],
          Rcpp::Named("convergence") = score["convergence"],
          Rcpp::Named("objective_evaluations") = score["objective_evaluations"],
          Rcpp::Named("gradient_evaluations") = score["gradient_evaluations"]);
        refined["backend"] =
          "native-bfgs-gq-score+native-bfgs-exact-grid-refinement";
        refined["coordinator"] = "persistent-native-cpp-gq";
        refined["objective_backend"] =
          "persistent-cpp-gaussian-quadrature-objective";
        selected = refined;
        refinement_accepted = true;
      }
    }
    const Rcpp::NumericVector selected_par = selected["par"];
    const Vector final_point = Eigen::Map<const Vector>(
      selected_par.begin(), static_cast<Eigen::Index>(selected_par.size()));
    const NativeGqEvaluation& final = evaluate_native(final_point, false);
    selected["population_objective"] = telemetry();
    selected["gradient_fallbacks"] = 0;
    selected["gradient_fallback_evaluations"] =
      static_cast<double>(finite_difference_evaluations_);
    return Rcpp::List::create(
      Rcpp::Named("optimizer") = selected,
      Rcpp::Named("modes") = libertad::eigen_matrix_to_r(final.modes),
      Rcpp::Named("effective_quadrature_points") = final.effective_points,
      Rcpp::Named("quadrature_cancellation_ratio") = final.cancellation_ratio,
      Rcpp::Named("exact_finite_grid_refinement") = refinement_accepted,
      Rcpp::Named("telemetry") = telemetry());
  }

  Rcpp::List telemetry() const {
    return Rcpp::List::create(
      Rcpp::Named("backend") = "persistent-native-cpp-gq-coordinator",
      Rcpp::Named("adaptive") = adaptive_,
      Rcpp::Named("value_requests") = static_cast<double>(value_requests_),
      Rcpp::Named("gradient_requests") = static_cast<double>(gradient_requests_),
      Rcpp::Named("parameter_evaluations") =
        static_cast<double>(parameter_evaluations_),
      Rcpp::Named("cache_hits") = static_cast<double>(cache_hits_),
      Rcpp::Named("proposal_refreshes") =
        static_cast<double>(proposal_refreshes_),
      Rcpp::Named("quadrature_node_evaluations") =
        static_cast<double>(node_evaluations_),
      Rcpp::Named("finite_difference_evaluations") =
        static_cast<double>(finite_difference_evaluations_),
      Rcpp::Named("stochastic_context") = context_->telemetry());
  }

 private:
  StochasticEtaCollection* context_ = nullptr;
  SEXP retained_context_ = R_NilValue;
  StochasticBayesMap map_;
  Matrix nodes_;
  Vector log_measure_, measure_sign_;
  bool adaptive_ = true;
  int eta_maxit_ = 100;
  double tolerance_ = 1e-7;
  bool cache_valid_ = false;
  bool cache_gradient_valid_ = false;
  Vector cache_key_;
  NativeGqEvaluation cache_;
  long long value_requests_ = 0;
  long long gradient_requests_ = 0;
  long long parameter_evaluations_ = 0;
  long long cache_hits_ = 0;
  long long proposal_refreshes_ = 0;
  long long node_evaluations_ = 0;
  long long finite_difference_evaluations_ = 0;

  bool same_key(const Vector& point) const {
    return cache_valid_ && point.size() == cache_key_.size() &&
      (point.array() == cache_key_.array()).all();
  }

  static double regularize(Matrix& matrix, const std::string& context) {
    matrix = 0.5 * (matrix + matrix.transpose()).eval();
    if (!matrix.allFinite()) throw std::domain_error(context + " is not finite.");
    const auto eigen = libertad::detail::self_adjoint_eigen(matrix, false);
    if (eigen.info != Eigen::Success || !eigen.values.allFinite()) {
      throw std::runtime_error(context + " decomposition failed.");
    }
    const double largest = std::max(eigen.values.cwiseAbs().maxCoeff(), 1.0);
    const double jitter = std::max(
      0.0, largest * 1e-9 - eigen.values.minCoeff());
    if (jitter > largest * 1e-2) {
      throw std::domain_error(context + " is not sufficiently positive definite.");
    }
    matrix.diagonal().array() += jitter;
    return jitter;
  }

  NativeGqEvaluation evaluate_uncached(const Vector& point, bool gradient) {
    if (!map_.in_bounds(point)) {
      NativeGqEvaluation invalid;
      invalid.native_gradient = Vector::Zero(point.size());
      return invalid;
    }
    const StochasticBayesParameters parameters = map_.decode(point);
    const int subjects = context_->subjects();
    const int n_eta = context_->eta_dimension();
    Matrix modes = Matrix::Zero(subjects, n_eta);
    std::vector<Matrix> roots(static_cast<std::size_t>(subjects));
    if (adaptive_) {
      Rcpp::NumericMatrix starts(subjects, n_eta);
      const Rcpp::List proposal = context_->laplace_proposal(
        Rcpp::wrap(parameters.theta), starts, Rcpp::wrap(parameters.sigma),
        Rcpp::wrap(parameters.omega), eta_maxit_, tolerance_);
      const Rcpp::NumericMatrix proposal_modes = proposal["modes"];
      modes = Eigen::Map<const Matrix>(
        proposal_modes.begin(), proposal_modes.nrow(), proposal_modes.ncol());
      const Rcpp::List root_list = proposal["roots"];
      if (root_list.size() != subjects) {
        throw std::runtime_error("Native GQ proposal count changed.");
      }
      for (int subject = 0; subject < subjects; ++subject) {
        const Rcpp::NumericMatrix current_root = root_list[subject];
        roots[static_cast<std::size_t>(subject)] = Eigen::Map<const Matrix>(
          current_root.begin(), current_root.nrow(), current_root.ncol());
      }
    } else {
      Matrix covariance = map_.omega_covariance(parameters);
      if (covariance.rows() != n_eta || covariance.cols() != n_eta) {
        throw std::invalid_argument(
          "Fixed native GQ currently requires an unexpanded OMEGA dimension.");
      }
      regularize(covariance, "fixed native GQ OMEGA covariance");
      Eigen::LLT<Matrix> factor(covariance);
      if (factor.info() != Eigen::Success) {
        throw std::runtime_error("Fixed native GQ OMEGA factorization failed.");
      }
      const Matrix root = Matrix(factor.matrixL());
      std::fill(roots.begin(), roots.end(), root);
    }
    ++proposal_refreshes_;
    NativeGqEvaluation result = context_->quadrature(
      parameters, nodes_, log_measure_, measure_sign_, modes, roots, gradient);
    node_evaluations_ += result.node_evaluations;
    Vector prior_gradient;
    const double prior = map_.prior_nll(
      parameters, gradient ? &prior_gradient : nullptr);
    if (!result.valid || !std::isfinite(prior) || prior >= 1e99) {
      result.value = 1e100;
      result.native_gradient = Vector::Zero(point.size());
      result.valid = false;
      return result;
    }
    result.value += prior;
    if (gradient) {
      Vector population = Vector::Zero(
        static_cast<Eigen::Index>(parameters.theta.size() +
          parameters.sigma.size() + parameters.omega.size()));
      Eigen::Index cursor = 0;
      for (std::size_t index = 0; index < parameters.theta.size(); ++index) {
        population[cursor++] = result.native_gradient[
          static_cast<Eigen::Index>(index)];
      }
      const Eigen::Index sigma_domain = static_cast<Eigen::Index>(
        parameters.theta.size() + static_cast<std::size_t>(n_eta));
      for (std::size_t index = 0; index < parameters.sigma.size(); ++index) {
        population[cursor++] = result.native_gradient[
          sigma_domain + static_cast<Eigen::Index>(index)];
      }
      const Eigen::Index omega_domain = sigma_domain +
        static_cast<Eigen::Index>(parameters.sigma.size());
      for (std::size_t index = 0; index < parameters.omega.size(); ++index) {
        population[cursor++] = result.native_gradient[
          omega_domain + static_cast<Eigen::Index>(index)];
      }
      population += prior_gradient;
      result.native_gradient = map_.outer_gradient(
        point, parameters, population);
    } else {
      result.native_gradient.resize(0);
    }
    return result;
  }

  const NativeGqEvaluation& evaluate_native(
      const Vector& point, bool gradient) {
    if (gradient) ++gradient_requests_; else ++value_requests_;
    if (same_key(point) && (!gradient || cache_gradient_valid_)) {
      ++cache_hits_;
      return cache_;
    }
    cache_ = evaluate_uncached(point, gradient);
    cache_key_ = point;
    cache_valid_ = true;
    cache_gradient_valid_ = gradient;
    ++parameter_evaluations_;
    return cache_;
  }

  Vector finite_difference(const Vector& point) {
    const double baseline = evaluate_native(point, false).value;
    Vector result(point.size());
    const std::vector<double>& lower = map_.lower();
    const std::vector<double>& upper = map_.upper();
    const auto valid = [](double value) {
      return std::isfinite(value) && value < 1e99;
    };
    for (Eigen::Index index = 0; index < point.size(); ++index) {
      const double step = 1e-5 * std::max(std::abs(point[index]), 1.0);
      Vector low = point;
      Vector high = point;
      low[index] = std::max(
        lower[static_cast<std::size_t>(index)], point[index] - step);
      high[index] = std::min(
        upper[static_cast<std::size_t>(index)], point[index] + step);
      const double low_value = low[index] < point[index] ?
        evaluate_native(low, false).value : baseline;
      if (low[index] < point[index]) ++finite_difference_evaluations_;
      const double high_value = high[index] > point[index] ?
        evaluate_native(high, false).value : baseline;
      if (high[index] > point[index]) ++finite_difference_evaluations_;
      if (low[index] < point[index] && high[index] > point[index] &&
          valid(low_value) && valid(high_value)) {
        result[index] = (high_value - low_value) /
          (high[index] - low[index]);
      } else if (high[index] > point[index] && valid(baseline) &&
                 valid(high_value)) {
        result[index] = (high_value - baseline) / (high[index] - point[index]);
      } else if (low[index] < point[index] && valid(baseline) &&
                 valid(low_value)) {
        result[index] = (baseline - low_value) / (point[index] - low[index]);
      } else {
        throw std::runtime_error("Native GQ finite-grid derivative is not finite.");
      }
    }
    return result;
  }

  Rcpp::List evaluation_list(
      const NativeGqEvaluation& value, const Vector& point,
      bool gradient) const {
    Rcpp::NumericVector gradient_value;
    SEXP gradient_output = R_NilValue;
    if (gradient) {
      gradient_value = libertad::eigen_vector_to_r(value.native_gradient);
      gradient_output = gradient_value;
    }
    return Rcpp::List::create(
      Rcpp::Named("value") = value.value,
      Rcpp::Named("gradient") = gradient_output,
      Rcpp::Named("par") = libertad::eigen_vector_to_r(point),
      Rcpp::Named("modes") = libertad::eigen_matrix_to_r(value.modes),
      Rcpp::Named("effective_quadrature_points") = value.effective_points,
      Rcpp::Named("quadrature_cancellation_ratio") = value.cancellation_ratio,
      Rcpp::Named("valid") = value.valid);
  }

  static Rcpp::List fixed_optimizer(double value) {
    return Rcpp::List::create(
      Rcpp::Named("par") = Rcpp::NumericVector(),
      Rcpp::Named("value") = value,
      Rcpp::Named("convergence") = 0,
      Rcpp::Named("message") = R_NilValue,
      Rcpp::Named("counts") = Rcpp::IntegerVector::create(
        Rcpp::Named("function") = 1, Rcpp::Named("gradient") = 0),
      Rcpp::Named("iterations") = 0,
      Rcpp::Named("objective_evaluations") = 1,
      Rcpp::Named("gradient_evaluations") = 0,
      Rcpp::Named("telemetry") = R_NilValue);
  }
};

// Persistent evaluator used by the compatibility SAEM M-step.  It retains
// tape pointers, fixed ETAs, and reusable point buffers across R optim()
// callbacks while leaving parameter decoding, priors, summation order, RNG,
// and the optimizer itself on the R side.
