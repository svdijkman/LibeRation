// Internal implementation fragment included inside namespace liberation by
// pk_engine.cpp. It intentionally has no include guard or namespace declaration.
// Contains the persistent stochastic ETA context shared by native BAYES, f-SAEM, and Gaussian quadrature.

class StochasticEtaCollection {
 public:
  StochasticEtaCollection(
      SEXP engine_pointer, const Rcpp::List& tape_pointers,
      const Rcpp::List& subject_data, int n_theta, int n_eta, int n_sigma,
      int n_omega, bool use_ode, const Rcpp::NumericVector& initial_theta,
      const Rcpp::NumericVector& initial_sigma,
      const Rcpp::NumericVector& initial_omega, double guard_radius,
      bool fused_values, int native_threads)
      : n_theta_(n_theta), n_eta_(n_eta), n_sigma_(n_sigma),
        n_omega_(n_omega), use_ode_(use_ode), guard_radius_(guard_radius),
        fused_requested_(fused_values),
        native_threads_(std::max(1, native_threads)) {
    if (tape_pointers.size() < 1 ||
        subject_data.size() != tape_pointers.size() ||
        n_theta < 0 || n_eta < 0 || n_sigma < 0 || n_omega < 0 ||
        initial_theta.size() != n_theta || initial_sigma.size() != n_sigma ||
        initial_omega.size() != n_omega || !std::isfinite(guard_radius_) ||
        guard_radius_ <= 0.0 || native_threads < 1) {
      throw std::invalid_argument(
        "Persistent stochastic subject inputs are inconsistent.");
    }
    Rcpp::XPtr<ModelEngine> engine(engine_pointer);
    engine_ = engine.get();
    retained_engine_ = engine_pointer;
    retained_subject_data_ = subject_data;
    R_PreserveObject(retained_engine_);
    R_PreserveObject(retained_subject_data_);
    domain_ = static_cast<std::size_t>(
      n_theta_ + n_eta_ + n_sigma_ + n_omega_);
    tapes_.reserve(static_cast<std::size_t>(tape_pointers.size()));
    subject_data_.reserve(static_cast<std::size_t>(subject_data.size()));
    points_.resize(static_cast<std::size_t>(tape_pointers.size()));
    anchors_.resize(static_cast<std::size_t>(tape_pointers.size()));
    owned_tapes_.resize(static_cast<std::size_t>(tape_pointers.size()));
    owned_subject_data_.resize(static_cast<std::size_t>(tape_pointers.size()));
    fused_subject_.assign(static_cast<std::size_t>(tape_pointers.size()), false);
    eta_positions_.reserve(static_cast<std::size_t>(n_eta_));
    for (int effect = 0; effect < n_eta_; ++effect) {
      eta_positions_.push_back(static_cast<std::size_t>(n_theta_ + effect));
    }
    StochasticBayesParameters initial_parameters;
    initial_parameters.theta = Rcpp::as<std::vector<double>>(initial_theta);
    initial_parameters.sigma = Rcpp::as<std::vector<double>>(initial_sigma);
    initial_parameters.omega = Rcpp::as<std::vector<double>>(initial_omega);
    const bool fused_model = fused_requested_ && engine_ &&
      use_specialized_advan(*engine_) && !use_ode_;
    for (int subject = 0; subject < tape_pointers.size(); ++subject) {
      subject_data_.push_back(subject_data[subject]);
      Rcpp::XPtr<ObjectiveTape> tape(tape_pointers[subject]);
      if (tape->domain_names.size() != domain_) {
        throw std::invalid_argument(
          "A persistent stochastic tape has an inconsistent domain.");
      }
      set_objective_dynamic_input(*tape, subject_data[subject]);
      tapes_.push_back(tape.get());
      points_[static_cast<std::size_t>(subject)].assign(domain_, 0.0);
      Vector eta = Vector::Zero(n_eta_);
      anchors_[static_cast<std::size_t>(subject)] =
        native_point(initial_parameters, eta);
      if (fused_model) {
        try {
          const EventDataView source = event_data_view(subject_data[subject]);
          owned_subject_data_[static_cast<std::size_t>(subject)] =
            std::make_shared<const OwnedEventTable>(
              source.parent(), source.start(), source.nrows());
          std::ostringstream messages;
          const auto agrees = [&](const StochasticBayesParameters& parameters,
                                  const Vector& current_eta) {
            const double direct = direct_value_raw(
              static_cast<std::size_t>(subject), parameters, current_eta);
            const std::vector<double> point = native_point(parameters, current_eta);
            const std::vector<double> recorded = tape->fun.Forward(
              0, point, messages);
            require_unchanged_path(tape->fun, "fused ADVAN eligibility check");
            const double reference = recorded.empty() ?
              std::numeric_limits<double>::infinity() : recorded[0];
            const double scale = 1.0 +
              std::max(std::abs(direct), std::abs(reference));
            return std::isfinite(direct) && std::isfinite(reference) &&
              std::abs(direct - reference) <= 5e-12 * scale;
          };
          bool eligible = agrees(initial_parameters, eta);
          if (eligible && n_eta_ > 0) {
            Vector positive = eta;
            Vector negative = eta;
            for (int effect = 0; effect < n_eta_; ++effect) {
              positive[effect] = 0.05;
              negative[effect] = -0.05;
            }
            eligible = agrees(initial_parameters, positive) &&
              agrees(initial_parameters, negative);
          }
          if (eligible) {
            StochasticBayesParameters shifted = initial_parameters;
            for (double& value : shifted.theta) {
              value += 0.01 * std::max(std::abs(value), 1.0);
            }
            for (double& value : shifted.sigma) value *= 1.01;
            for (double& value : shifted.omega) value *= 1.01;
            eligible = agrees(shifted, eta);
          }
          if (eligible) {
            fused_subject_[static_cast<std::size_t>(subject)] = true;
          } else {
            ++fused_guard_failures_;
          }
        } catch (const std::exception& error) {
          ++fused_guard_failures_;
          if (fused_fallback_reason_.empty()) fused_fallback_reason_ = error.what();
          owned_subject_data_[static_cast<std::size_t>(subject)].reset();
        }
      }
    }
    fused_enabled_ = !fused_subject_.empty() &&
      std::all_of(fused_subject_.begin(), fused_subject_.end(),
                  [](bool value) { return value; });
    if (fused_requested_ && !fused_enabled_ && fused_fallback_reason_.empty()) {
      fused_fallback_reason_ = fused_model ?
        "direct and recorded objectives did not agree at the eligibility point" :
        "model is not an eligible specialised analytical ADVAN";
    }
    if (!fused_enabled_) native_threads_ = 1;
    if (fused_enabled_ && native_threads_ > 1) {
      subject_pool_ = std::make_unique<NativeSubjectPool>(
        std::min(native_threads_, static_cast<int>(tapes_.size())));
    }
  }

  ~StochasticEtaCollection() {
    // Join native workers before releasing the engine and owned event tables
    // they may read during a dispatch.
    subject_pool_.reset();
    if (retained_engine_ != R_NilValue) R_ReleaseObject(retained_engine_);
    if (retained_subject_data_ != R_NilValue) {
      R_ReleaseObject(retained_subject_data_);
    }
  }

  StochasticEtaCollection(const StochasticEtaCollection&) = delete;
  StochasticEtaCollection& operator=(const StochasticEtaCollection&) = delete;

  int subjects() const { return static_cast<int>(tapes_.size()); }
  int eta_dimension() const { return n_eta_; }

  Rcpp::NumericVector evaluate(
      const Rcpp::NumericVector& theta, const Rcpp::NumericMatrix& eta,
      const Rcpp::NumericVector& sigma, const Rcpp::NumericVector& omega) {
    validate_parameters(theta, eta, sigma, omega);
    StochasticBayesParameters parameters;
    parameters.theta = Rcpp::as<std::vector<double>>(theta);
    parameters.sigma = Rcpp::as<std::vector<double>>(sigma);
    parameters.omega = Rcpp::as<std::vector<double>>(omega);
    Rcpp::NumericVector values(static_cast<R_xlen_t>(tapes_.size()));
    if (fused_enabled_) {
      Matrix eta_native(eta.nrow(), eta.ncol());
      for (int subject = 0; subject < eta.nrow(); ++subject) {
        for (int effect = 0; effect < eta.ncol(); ++effect) {
          eta_native(subject, effect) = eta(subject, effect);
        }
      }
      std::vector<double> native_values(tapes_.size());
      parallel_subjects(tapes_.size(), [&](std::size_t subject) {
        native_values[subject] = direct_value_raw(
          subject, parameters,
          eta_native.row(static_cast<Eigen::Index>(subject)).transpose());
      });
      for (std::size_t subject = 0; subject < native_values.size(); ++subject) {
        values[static_cast<R_xlen_t>(subject)] = native_values[subject];
      }
      evaluations_ += static_cast<long long>(native_values.size());
      fused_evaluations_ += static_cast<long long>(native_values.size());
      return values;
    }
    for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
      fill_point(subject, theta, eta, sigma, omega);
      Vector subject_eta(n_eta_);
      for (int effect = 0; effect < n_eta_; ++effect) {
        subject_eta[effect] = eta(static_cast<int>(subject), effect);
      }
      const double value = guarded_value(
        subject, parameters, subject_eta, points_[subject],
        "persistent stochastic objective");
      values[static_cast<R_xlen_t>(subject)] = value;
      if ((subject + 1U) % 256U == 0U) Rcpp::checkUserInterrupt();
    }
    return values;
  }

  Rcpp::List laplace_proposal(
      const Rcpp::NumericVector& theta,
      const Rcpp::NumericMatrix& starts,
      const Rcpp::NumericVector& sigma,
      const Rcpp::NumericVector& omega,
      int maxit, double tolerance) {
    validate_parameters(theta, starts, sigma, omega);
    if (maxit < 1 || !std::isfinite(tolerance) || tolerance <= 0.0) {
      throw std::invalid_argument(
        "Persistent f-SAEM mode controls are invalid.");
    }
    Rcpp::NumericMatrix modes(
      static_cast<int>(tapes_.size()), n_eta_);
    Rcpp::NumericVector values(static_cast<R_xlen_t>(tapes_.size()));
    Rcpp::List roots(static_cast<int>(tapes_.size()));
    Rcpp::List precisions(static_cast<int>(tapes_.size()));
    Rcpp::NumericVector jitters(static_cast<R_xlen_t>(tapes_.size()));
    StochasticBayesParameters native_parameters;
    native_parameters.theta = Rcpp::as<std::vector<double>>(theta);
    native_parameters.sigma = Rcpp::as<std::vector<double>>(sigma);
    native_parameters.omega = Rcpp::as<std::vector<double>>(omega);
    int total_iterations = 0;
    int total_evaluations = 0;
    int total_gradient_evaluations = 0;
    int restarts = 0;
    int relative_convergence = 0;
    int newton_convergence = 0;
    for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
      fill_point(subject, theta, starts, sigma, omega);
      Rcpp::NumericVector start(n_eta_);
      for (int effect = 0; effect < n_eta_; ++effect) {
        start[effect] = starts(static_cast<int>(subject), effect);
      }
      Rcpp::List mode;
      Matrix accepted_hessian;
      bool has_hessian = false;
      bool converged = false;
      for (int restart = 0; restart < 8 && !converged; ++restart) {
        // Any curvature calculated below belongs to the candidate produced by
        // this mode search only.  A restarted search must never reuse it.
        has_hessian = false;
        try {
          mode = objective_eta_mode(
            *tapes_[subject], points_[subject], eta_positions_, start,
            maxit, tolerance, false);
        } catch (const TapePathChange& change) {
          Vector anchor(n_eta_);
          for (int effect = 0; effect < n_eta_; ++effect) {
            const std::size_t position = eta_positions_[
              static_cast<std::size_t>(effect)];
            anchor[effect] = change.point().size() > position ?
              change.point()[position] : start[effect];
            start[effect] = anchor[effect];
          }
          record_subject_tape(subject, native_parameters, anchor, true);
          fill_point(subject, theta, starts, sigma, omega);
          ++restarts;
          continue;
        }
        total_iterations += Rcpp::as<int>(mode["iterations"]);
        total_evaluations += Rcpp::as<int>(mode["evaluations"]);
        total_gradient_evaluations +=
          Rcpp::as<int>(mode["gradient_evaluations"]);
        const Rcpp::NumericVector par = mode["par"];
        const Rcpp::NumericVector gradient = mode["gradient"];
        const double value = Rcpp::as<double>(mode["value"]);
        std::vector<double> mode_point = points_[subject];
        Vector eta_eigen(n_eta_), gradient_eigen(n_eta_);
        double gradient_norm = 0.0;
        for (int effect = 0; effect < n_eta_; ++effect) {
          eta_eigen[effect] = par[effect];
          gradient_eigen[effect] = gradient[effect];
          gradient_norm = std::max(
            gradient_norm, std::abs(gradient_eigen[effect]));
          mode_point[eta_positions_[static_cast<std::size_t>(effect)]] =
            eta_eigen[effect];
        }
        converged = Rcpp::as<int>(mode["convergence"]) == 0;
        if (!converged && std::isfinite(value) &&
            std::isfinite(gradient_norm) &&
            gradient_norm <= std::max(
              10.0 * tolerance, tolerance * (1.0 + std::abs(value)))) {
          converged = true;
          ++relative_convergence;
        }
        if (!converged && std::isfinite(value) &&
            gradient_eigen.allFinite()) {
          try {
            accepted_hessian = objective_eta_hessian(
              *tapes_[subject], mode_point, eta_positions_);
            const double jitter = regularize_curvature(
              accepted_hessian, "f-SAEM conditional curvature");
            const Vector displacement =
              accepted_hessian.ldlt().solve(gradient_eigen);
            has_hessian = true;
            if (displacement.allFinite() &&
                displacement.lpNorm<Eigen::Infinity>() <=
                  std::sqrt(tolerance) *
                  (1.0 + eta_eigen.lpNorm<Eigen::Infinity>())) {
              converged = true;
              ++newton_convergence;
              jitters[static_cast<R_xlen_t>(subject)] = jitter;
            }
          } catch (const std::exception&) {
            has_hessian = false;
          }
        }
        if (!converged) {
          start = Rcpp::clone(par);
          ++restarts;
        }
      }
      if (!converged) {
        throw std::runtime_error(
          "Persistent f-SAEM conditional mode failed for subject " +
          std::to_string(subject + 1U) + ".");
      }
      const Rcpp::NumericVector par = mode["par"];
      std::vector<double> mode_point = points_[subject];
      for (int effect = 0; effect < n_eta_; ++effect) {
        modes(static_cast<int>(subject), effect) = par[effect];
        mode_point[eta_positions_[static_cast<std::size_t>(effect)]] =
          par[effect];
      }
      if (!has_hessian) {
        try {
          accepted_hessian = objective_eta_hessian(
            *tapes_[subject], mode_point, eta_positions_);
        } catch (const TapePathChange&) {
          Vector anchor(n_eta_);
          for (int effect = 0; effect < n_eta_; ++effect) {
            anchor[effect] = par[effect];
          }
          record_subject_tape(subject, native_parameters, anchor, true);
          accepted_hessian = objective_eta_hessian(
            *tapes_[subject], mode_point, eta_positions_);
        }
        jitters[static_cast<R_xlen_t>(subject)] = regularize_curvature(
          accepted_hessian, "f-SAEM conditional curvature");
      }
      Eigen::LLT<Matrix> precision_factor(accepted_hessian);
      if (precision_factor.info() != Eigen::Success) {
        throw std::runtime_error(
          "f-SAEM conditional precision factorization failed.");
      }
      Matrix covariance = 2.0 * precision_factor.solve(
        Matrix::Identity(n_eta_, n_eta_));
      covariance = 0.5 * (covariance + covariance.transpose()).eval();
      Eigen::LLT<Matrix> covariance_factor(covariance);
      if (covariance_factor.info() != Eigen::Success) {
        throw std::runtime_error(
          "f-SAEM proposal covariance factorization failed.");
      }
      roots[static_cast<int>(subject)] = libertad::eigen_matrix_to_r(
        Matrix(covariance_factor.matrixL()));
      precisions[static_cast<int>(subject)] = libertad::eigen_matrix_to_r(
        0.5 * accepted_hessian);
      values[static_cast<R_xlen_t>(subject)] =
        Rcpp::as<double>(mode["value"]);
      if ((subject + 1U) % 32U == 0U) Rcpp::checkUserInterrupt();
    }
    laplace_refreshes_ += 1;
    laplace_mode_evaluations_ += total_evaluations;
    return Rcpp::List::create(
      Rcpp::Named("modes") = modes,
      Rcpp::Named("values") = values,
      Rcpp::Named("roots") = roots,
      Rcpp::Named("precisions") = precisions,
      Rcpp::Named("jitter") = jitters,
      Rcpp::Named("mode_iterations") = total_iterations,
      Rcpp::Named("mode_evaluations") = total_evaluations,
      Rcpp::Named("gradient_evaluations") = total_gradient_evaluations,
      Rcpp::Named("restarts") = restarts,
      Rcpp::Named("relative_convergence") = relative_convergence,
      Rcpp::Named("newton_convergence") = newton_convergence,
      Rcpp::Named("backend") =
        "persistent-cpp-laplace-proposal");
  }

  // Evaluate a complete fixed or adaptive Gaussian-quadrature grid.  Proposal
  // modes and roots are supplied by the native coordinator; the signed
  // log-sum-exp reduction and normalized score derivative are kept here next
  // to the retained subject tapes so no per-node R objects or callbacks are
  // needed during outer optimization.
  NativeGqEvaluation quadrature(
      const StochasticBayesParameters& parameters, const Matrix& nodes,
      const Vector& log_measure, const Vector& measure_sign,
      const Matrix& modes, const std::vector<Matrix>& roots,
      bool gradient) {
    const Eigen::Index draws = nodes.rows();
    if (draws < 1 || nodes.cols() != n_eta_ ||
        log_measure.size() != draws || measure_sign.size() != draws ||
        modes.rows() != static_cast<Eigen::Index>(tapes_.size()) ||
        modes.cols() != n_eta_ || roots.size() != tapes_.size()) {
      throw std::invalid_argument("Native GQ proposal dimensions are inconsistent.");
    }
    for (const Matrix& root : roots) {
      if (root.rows() != n_eta_ || root.cols() != n_eta_ ||
          !root.allFinite()) {
        throw std::invalid_argument("A native GQ proposal root is invalid.");
      }
    }
    NativeGqEvaluation result;
    result.value = 0.0;
    result.native_gradient = Vector::Zero(static_cast<Eigen::Index>(domain_));
    result.modes = modes;
    result.effective_points.assign(tapes_.size(), 0.0);
    result.cancellation_ratio.assign(tapes_.size(), 0.0);
    result.valid = true;
    const double log_two_pi = std::log(2.0 * std::acos(-1.0));
    const std::vector<double> reverse_weight(1U, 1.0);

    for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
      const Matrix& root = roots[subject];
      double logdet = 0.0;
      for (int effect = 0; effect < n_eta_; ++effect) {
        const double diagonal = root(effect, effect);
        if (!(diagonal > 0.0) || !std::isfinite(diagonal)) {
          throw std::domain_error("A native GQ proposal factor is singular.");
        }
        logdet += 2.0 * std::log(diagonal);
      }
      std::vector<double> log_integrand(static_cast<std::size_t>(draws));
      Matrix sample_gradient;
      if (gradient) {
        sample_gradient = Matrix::Zero(draws, static_cast<Eigen::Index>(domain_));
      }
      double maximum = -std::numeric_limits<double>::infinity();
      for (Eigen::Index draw = 0; draw < draws; ++draw) {
        const Vector eta = modes.row(static_cast<Eigen::Index>(subject)).transpose() +
          root * nodes.row(draw).transpose();
        std::vector<double>& point = points_[subject];
        fill_point_native(point, parameters, eta);
        if (material_movement(subject, parameters, eta)) {
          record_subject_tape(subject, parameters, eta, true);
        }
        double objective = std::numeric_limits<double>::infinity();
        Vector derivative;
        bool evaluated = false;
        for (int attempt = 0; attempt < 4 && !evaluated; ++attempt) {
          try {
            std::ostringstream messages;
            const std::vector<double> value = tapes_[subject]->fun.Forward(
              0, point, messages);
            require_unchanged_path(tapes_[subject]->fun, "native GQ objective");
            objective = value.empty() ?
              std::numeric_limits<double>::infinity() : value[0];
            if (gradient) {
              const std::vector<double> current = tapes_[subject]->fun.Reverse(
                1, reverse_weight);
              require_unchanged_path(tapes_[subject]->fun, "native GQ score");
              if (current.size() != domain_) {
                throw std::runtime_error(
                  "A native GQ tape returned an invalid gradient length.");
              }
              derivative = Eigen::Map<const Vector>(
                current.data(), static_cast<Eigen::Index>(current.size()));
            }
            evaluated = true;
          } catch (const TapePathChange&) {
            record_subject_tape(subject, parameters, eta, true);
          }
        }
        if (!evaluated) {
          throw std::runtime_error(
            "A GQ objective remained structurally unstable after retaping.");
        }
        ++evaluations_;
        ++recorded_evaluations_;
        ++result.node_evaluations;
        const double log_proposal = -0.5 * (
          static_cast<double>(n_eta_) * log_two_pi + logdet +
          nodes.row(draw).squaredNorm());
        const double integrand = -0.5 * objective - log_proposal +
          log_measure[draw];
        log_integrand[static_cast<std::size_t>(draw)] = integrand;
        if (std::isfinite(integrand) && std::isfinite(measure_sign[draw]) &&
            measure_sign[draw] != 0.0) {
          maximum = std::max(maximum, integrand);
        }
        if (gradient && derivative.size() == static_cast<Eigen::Index>(domain_)) {
          sample_gradient.row(draw) = derivative.transpose();
        }
      }

      double signed_total = 0.0;
      double absolute_total = 0.0;
      std::vector<double> scaled(static_cast<std::size_t>(draws), 0.0);
      if (std::isfinite(maximum)) {
        for (Eigen::Index draw = 0; draw < draws; ++draw) {
          const double integrand = log_integrand[static_cast<std::size_t>(draw)];
          if (!std::isfinite(integrand) ||
              !std::isfinite(measure_sign[draw]) ||
              measure_sign[draw] == 0.0) continue;
          const double current = measure_sign[draw] *
            std::exp(integrand - maximum);
          scaled[static_cast<std::size_t>(draw)] = current;
          signed_total += current;
          absolute_total += std::abs(current);
        }
      }
      const bool valid = std::isfinite(signed_total) &&
        std::isfinite(absolute_total) && signed_total >
          std::numeric_limits<double>::epsilon() *
            std::max(1.0, absolute_total);
      if (!valid) {
        result.value = 1e100;
        result.native_gradient.setZero();
        result.valid = false;
        return result;
      }
      result.value += -2.0 * (maximum + std::log(signed_total));
      double squared_absolute_weights = 0.0;
      for (Eigen::Index draw = 0; draw < draws; ++draw) {
        const double absolute_weight =
          std::abs(scaled[static_cast<std::size_t>(draw)]) / absolute_total;
        squared_absolute_weights += absolute_weight * absolute_weight;
        if (gradient) {
          const double normalized =
            scaled[static_cast<std::size_t>(draw)] / signed_total;
          result.native_gradient.noalias() +=
            normalized * sample_gradient.row(draw).transpose();
        }
      }
      result.effective_points[subject] =
        squared_absolute_weights > 0.0 ? 1.0 / squared_absolute_weights : 0.0;
      result.cancellation_ratio[subject] = signed_total / absolute_total;
      if ((subject + 1U) % 32U == 0U) Rcpp::checkUserInterrupt();
    }
    ++calls_;
    return result;
  }

  Rcpp::List random_walk(
      const Rcpp::NumericVector& theta, const Rcpp::NumericMatrix& eta_input,
      const Rcpp::NumericVector& sigma, const Rcpp::NumericVector& omega,
      const Rcpp::List& proposal_roots, const Rcpp::NumericMatrix& normals,
      const Rcpp::NumericVector& log_uniforms, int mcmc_steps,
      double step_scale,
      Rcpp::Nullable<Rcpp::NumericVector> current_values_input = R_NilValue) {
    validate_sampler(
      theta, eta_input, sigma, omega, proposal_roots, normals,
      log_uniforms, mcmc_steps);
    if (!std::isfinite(step_scale) || step_scale <= 0.0) {
      throw std::invalid_argument(
        "Persistent stochastic random-walk scale must be positive.");
    }
    return metropolis(
      theta, eta_input, sigma, omega, proposal_roots, R_NilValue,
      R_NilValue, normals, log_uniforms, mcmc_steps, step_scale,
      current_values_input, false, R_PosInf, R_NilValue);
  }

  Rcpp::List laplace_independence(
      const Rcpp::NumericVector& theta, const Rcpp::NumericMatrix& eta_input,
      const Rcpp::NumericVector& sigma, const Rcpp::NumericVector& omega,
      const Rcpp::NumericMatrix& proposal_modes,
      const Rcpp::List& proposal_roots,
      const Rcpp::List& proposal_precisions,
      const Rcpp::NumericMatrix& normals,
      const Rcpp::NumericVector& log_uniforms, int mcmc_steps,
      Rcpp::Nullable<Rcpp::NumericVector> current_values_input = R_NilValue,
      double proposal_df = R_PosInf,
      Rcpp::Nullable<Rcpp::NumericVector> proposal_scales_input = R_NilValue) {
    validate_sampler(
      theta, eta_input, sigma, omega, proposal_roots, normals,
      log_uniforms, mcmc_steps);
    if (proposal_modes.nrow() != static_cast<int>(tapes_.size()) ||
        proposal_modes.ncol() != n_eta_ ||
        proposal_precisions.size() != static_cast<int>(tapes_.size())) {
      throw std::invalid_argument(
        "Laplace independence proposals must match subjects and ETAs.");
    }
    const bool student_t = std::isfinite(proposal_df);
    if (student_t && proposal_df <= 2.0) {
      throw std::invalid_argument(
        "Student-t Laplace proposals require degrees of freedom above two.");
    }
    return metropolis(
      theta, eta_input, sigma, omega, proposal_roots, proposal_modes,
      proposal_precisions, normals, log_uniforms, mcmc_steps, 1.0,
      current_values_input, true, proposal_df, proposal_scales_input);
  }

  Rcpp::List bayes_sample(
      const Rcpp::List& map_config, int n_burn, int n_sample, int n_thin,
      double step_scale, double eta_step, bool adapt,
      const std::string& outer_kernel, int adaptive_start,
      int adaptive_interval, double target_acceptance,
      double delayed_rejection_scale,
      const std::string& eta_kernel, int eta_refresh, int eta_maxit,
      double eta_tolerance, double eta_df, double eta_rescue_probability,
      double eta_parameter_refresh, double eta_low_acceptance,
      bool gibbs_omega) {
    if (n_burn < 0 || n_sample < 1 || n_thin < 1 ||
        !std::isfinite(step_scale) || step_scale <= 0.0 ||
        !std::isfinite(eta_step) || eta_step <= 0.0 ||
        adaptive_start < 2 || adaptive_interval < 1 ||
        !std::isfinite(target_acceptance) || target_acceptance <= 0.0 ||
        target_acceptance >= 1.0 ||
        !std::isfinite(delayed_rejection_scale) ||
        delayed_rejection_scale < 0.0 || delayed_rejection_scale >= 1.0 ||
        (outer_kernel != "isotropic" &&
         outer_kernel != "adaptive_metropolis") ||
        (eta_kernel != "random_walk" && eta_kernel != "laplace" &&
         eta_kernel != "student_t") ||
        eta_refresh < 1 || eta_maxit < 1 || !std::isfinite(eta_tolerance) ||
        eta_tolerance <= 0.0 || !std::isfinite(eta_df) || eta_df <= 2.0 ||
        !std::isfinite(eta_rescue_probability) ||
        eta_rescue_probability < 0.0 || eta_rescue_probability >= 1.0 ||
        !std::isfinite(eta_parameter_refresh) || eta_parameter_refresh <= 0.0 ||
        !std::isfinite(eta_low_acceptance) || eta_low_acceptance < 0.0 ||
        eta_low_acceptance >= 1.0) {
      throw std::invalid_argument("Native BAYES controls are invalid.");
    }
    StochasticBayesMap map(map_config);
    const Rcpp::List mu_input = map_config.containsElementNamed("mu") ?
      Rcpp::List(map_config["mu"]) : Rcpp::List::create(
        Rcpp::Named("active") = false);
    const StochasticMuConfig mu(
      mu_input, static_cast<int>(tapes_.size()), n_eta_);
    Vector outer(static_cast<Eigen::Index>(map.start().size()));
    for (std::size_t index = 0; index < map.start().size(); ++index) {
      outer[static_cast<Eigen::Index>(index)] = map.start()[index];
    }
    StochasticBayesParameters parameters = map.decode(outer);
    Matrix eta = Matrix::Zero(
      static_cast<Eigen::Index>(tapes_.size()), n_eta_);
    std::vector<double> subject_values;
    double current = bayes_log_posterior(
      map, parameters, eta, subject_values);
    if (!std::isfinite(current)) {
      throw std::domain_error("Initial native BAYES posterior is not finite.");
    }

    const int total_iterations = n_burn + n_sample * n_thin;
    const int n_native = n_theta_ + n_sigma_ + n_omega_;
    const int output_columns = n_native +
      static_cast<int>(tapes_.size()) * n_eta_ + 1;
    Rcpp::NumericMatrix chain(n_sample, output_columns);
    std::vector<int> random_positions;
    random_positions.reserve(map.dimension());
    for (int position = 0; position < static_cast<int>(map.dimension()); ++position) {
      bool linked = false;
      if (mu.active) {
        for (int theta : mu.theta) {
          if (map.theta_outer_position(theta) == position) {
            linked = true;
            break;
          }
        }
      }
      if (!linked) random_positions.push_back(position);
    }
    struct GibbsOmega {
      int native_index;
      int outer_position;
      int effect;
      double shape;
      double rate;
    };
    std::vector<GibbsOmega> gibbs_omegas;
    if (gibbs_omega && n_eta_ == map.omega_covariance(parameters).rows()) {
      for (int native_index = 0;
           native_index < static_cast<int>(parameters.omega.size());
           ++native_index) {
        double shape = 0.0, rate = 0.0;
        const int outer_position = map.omega_outer_position(native_index);
        const int effect = map.omega_effect(native_index);
        if (outer_position >= 0 && effect >= 0 &&
            map.diagonal_omega_inverse_gamma(native_index, shape, rate)) {
          gibbs_omegas.push_back(GibbsOmega{
            native_index, outer_position, effect, shape, rate});
        }
      }
      for (const GibbsOmega& update : gibbs_omegas) {
        random_positions.erase(std::remove(
          random_positions.begin(), random_positions.end(),
          update.outer_position), random_positions.end());
      }
    }
    const int dimension = static_cast<int>(random_positions.size());
    Matrix proposal_root = Matrix::Identity(dimension, dimension) * step_scale;
    Matrix proposal_covariance = proposal_root * proposal_root.transpose();
    Vector adaptive_mean = Vector::Zero(dimension);
    Matrix adaptive_m2 = Matrix::Zero(dimension, dimension);
    double log_multiplier = 0.0;
    int adaptive_n = 0;
    int covariance_updates = 0;
    int covariance_regularizations = 0;
    int accepted_outer = 0, attempted_outer = 0;
    int accepted_delayed = 0, attempted_delayed = 0;
    int accepted_eta = 0, attempted_eta = 0;
    int accepted_mu = 0, attempted_mu = 0;
    int gibbs_omega_updates = 0, gibbs_omega_draws = 0;
    int keep = 0;
    int omega_factorizations = 0, omega_cache_hits = 0;
    std::vector<double> cached_omega;
    std::vector<Matrix> eta_roots(tapes_.size());
    Rcpp::NumericMatrix eta_proposal_modes;
    Rcpp::List eta_proposal_roots;
    Rcpp::List eta_proposal_precisions;
    std::vector<double> eta_proposal_anchor;
    bool eta_force_refresh = false;
    int eta_refreshes = 0, eta_refresh_failures = 0;
    int eta_rescue_iterations = 0, eta_parameter_refreshes = 0;
    int eta_acceptance_refreshes = 0, eta_fallback_iterations = 0;
    std::string eta_last_error;

    for (int iteration = 1; iteration <= total_iterations; ++iteration) {
      bool accepted_outer_iteration = false;
      if (dimension > 0) {
        Vector normals(dimension);
        for (int index = 0; index < dimension; ++index) {
          normals[index] = R::rnorm(0.0, 1.0);
        }
        Vector candidate_outer = outer;
        const Vector increment = proposal_root * normals;
        for (int index = 0; index < dimension; ++index) {
          candidate_outer[random_positions[static_cast<std::size_t>(index)]] +=
            increment[index];
        }
        double candidate_logp = -std::numeric_limits<double>::infinity();
        StochasticBayesParameters candidate_parameters;
        std::vector<double> candidate_values;
        if (map.in_bounds(candidate_outer)) {
          try {
            candidate_parameters = map.decode(candidate_outer);
            candidate_logp = bayes_log_posterior(
              map, candidate_parameters, eta, candidate_values);
          } catch (const std::exception&) {
            candidate_logp = -std::numeric_limits<double>::infinity();
          }
        }
        ++attempted_outer;
        const double first_log_alpha = std::min(0.0, candidate_logp - current);
        if (std::log(R::runif(0.0, 1.0)) < first_log_alpha) {
          outer.swap(candidate_outer);
          parameters = std::move(candidate_parameters);
          subject_values.swap(candidate_values);
          current = candidate_logp;
          ++accepted_outer;
          accepted_outer_iteration = true;
        } else if (delayed_rejection_scale > 0.0) {
          ++attempted_delayed;
          Vector second_normal(dimension);
          for (int index = 0; index < dimension; ++index) {
            second_normal[index] = R::rnorm(0.0, 1.0);
          }
          Vector second_outer = outer;
          const Vector second_increment = delayed_rejection_scale *
            proposal_root * second_normal;
          for (int index = 0; index < dimension; ++index) {
            second_outer[random_positions[static_cast<std::size_t>(index)]] +=
              second_increment[index];
          }
          double second_logp = -std::numeric_limits<double>::infinity();
          StochasticBayesParameters second_parameters;
          std::vector<double> second_values;
          if (map.in_bounds(second_outer)) {
            try {
              second_parameters = map.decode(second_outer);
              second_logp = bayes_log_posterior(
                map, second_parameters, eta, second_values);
            } catch (const std::exception&) {
              second_logp = -std::numeric_limits<double>::infinity();
            }
          }
          ++attempted_outer;
          if (std::isfinite(second_logp)) {
            Vector first_from_current(dimension);
            Vector first_from_second(dimension);
            for (int index = 0; index < dimension; ++index) {
              const int position = random_positions[
                static_cast<std::size_t>(index)];
              first_from_current[index] = candidate_outer[position] - outer[position];
              first_from_second[index] = candidate_outer[position] -
                second_outer[position];
            }
            const double reverse_first_alpha = std::min(
              0.0, candidate_logp - second_logp);
            const double correction = second_logp - current -
              0.5 * gaussian_quadratic(first_from_second, proposal_root) +
              0.5 * gaussian_quadratic(first_from_current, proposal_root) +
              log_one_minus_acceptance(reverse_first_alpha) -
              log_one_minus_acceptance(first_log_alpha);
            if (std::log(R::runif(0.0, 1.0)) < std::min(0.0, correction)) {
              outer.swap(second_outer);
              parameters = std::move(second_parameters);
              subject_values.swap(second_values);
              current = second_logp;
              ++accepted_outer;
              ++accepted_delayed;
              accepted_outer_iteration = true;
            }
          }
        }
      }

      if (mu.active) {
        ++attempted_mu;
        if (bayes_mu_step(
              mu, map, outer, parameters, eta, subject_values, current)) {
          ++accepted_mu;
        }
      }

      for (const GibbsOmega& update : gibbs_omegas) {
        double sum_squares = 0.0;
        for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
          const double value = eta(
            static_cast<Eigen::Index>(subject), update.effect);
          sum_squares += value * value;
        }
        const double posterior_shape = update.shape +
          0.5 * static_cast<double>(tapes_.size());
        const double posterior_rate = update.rate + 0.5 * sum_squares;
        ++gibbs_omega_updates;
        bool updated = false;
        for (int attempt = 0; attempt < 10000 && !updated; ++attempt) {
          ++gibbs_omega_draws;
          const double precision = R::rgamma(
            posterior_shape, 1.0 / posterior_rate);
          if (!(precision > 0.0) || !std::isfinite(precision)) continue;
          StochasticBayesParameters candidate_parameters = parameters;
          candidate_parameters.omega[static_cast<std::size_t>(
            update.native_index)] = 1.0 / precision;
          Vector candidate_outer;
          try {
            candidate_outer = map.encode(candidate_parameters);
          } catch (const std::exception&) {
            continue;
          }
          std::vector<double> candidate_values;
          const double candidate_logp = bayes_log_posterior(
            map, candidate_parameters, eta, candidate_values);
          if (!std::isfinite(candidate_logp)) continue;
          outer.swap(candidate_outer);
          parameters = std::move(candidate_parameters);
          subject_values.swap(candidate_values);
          current = candidate_logp;
          updated = true;
        }
        if (!updated) {
          throw std::runtime_error(
            "Conjugate OMEGA update could not draw inside parameter bounds.");
        }
      }

      bool eta_independence = false;
      if (n_eta_ > 0 && eta_kernel != "random_walk") {
        std::vector<double> anchor;
        anchor.reserve(parameters.theta.size() + parameters.sigma.size() +
                       parameters.omega.size());
        anchor.insert(anchor.end(), parameters.theta.begin(), parameters.theta.end());
        anchor.insert(anchor.end(), parameters.sigma.begin(), parameters.sigma.end());
        anchor.insert(anchor.end(), parameters.omega.begin(), parameters.omega.end());
        bool parameter_refresh = eta_proposal_anchor.size() != anchor.size();
        if (!parameter_refresh && !anchor.empty()) {
          double movement = 0.0;
          for (std::size_t index = 0; index < anchor.size(); ++index) {
            movement = std::max(
              movement, std::abs(anchor[index] - eta_proposal_anchor[index]) /
                (1.0 + std::abs(eta_proposal_anchor[index])));
          }
          parameter_refresh = movement > eta_parameter_refresh;
        }
        const bool scheduled = eta_proposal_anchor.empty() ||
          (iteration - 1) % eta_refresh == 0;
        if (scheduled || parameter_refresh || eta_force_refresh) {
          if (parameter_refresh && !eta_proposal_anchor.empty()) {
            ++eta_parameter_refreshes;
          }
          try {
            Rcpp::NumericMatrix starts = libertad::eigen_matrix_to_r(eta);
            const Rcpp::List proposal = laplace_proposal(
              Rcpp::wrap(parameters.theta), starts,
              Rcpp::wrap(parameters.sigma), Rcpp::wrap(parameters.omega),
              eta_maxit, eta_tolerance);
            eta_proposal_modes = Rcpp::as<Rcpp::NumericMatrix>(proposal["modes"]);
            eta_proposal_roots = Rcpp::as<Rcpp::List>(proposal["roots"]);
            eta_proposal_precisions = Rcpp::as<Rcpp::List>(proposal["precisions"]);
            eta_proposal_anchor = anchor;
            eta_force_refresh = false;
            ++eta_refreshes;
          } catch (const std::exception& error) {
            ++eta_refresh_failures;
            eta_last_error = error.what();
          }
        }
        eta_independence = eta_proposal_modes.nrow() ==
          static_cast<int>(tapes_.size());
        if (!eta_independence) {
          ++eta_fallback_iterations;
        } else if (eta_rescue_probability > 0.0 &&
                   R::runif(0.0, 1.0) < eta_rescue_probability) {
          eta_independence = false;
          ++eta_rescue_iterations;
        }
      }

      if (n_eta_ > 0) {
        const int accepted_eta_before = accepted_eta;
        const int attempted_eta_before = attempted_eta;
        if (!eta_independence && cached_omega != parameters.omega) {
          std::vector<Matrix> covariance_groups;
          std::vector<Matrix> root_groups;
          for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
            Matrix covariance = subject_eta_covariance(
              map, parameters, subject);
            regularize_curvature(covariance, "native BAYES ETA covariance");
            std::size_t group = covariance_groups.size();
            for (std::size_t candidate = 0;
                 candidate < covariance_groups.size(); ++candidate) {
              const Matrix& established = covariance_groups[candidate];
              if (established.rows() == covariance.rows() &&
                  established.cols() == covariance.cols() &&
                  (established.array() == covariance.array()).all()) {
                group = candidate;
                break;
              }
            }
            if (group < root_groups.size()) {
              eta_roots[subject] = root_groups[group];
              continue;
            }
            Eigen::LLT<Matrix> factor(covariance);
            if (factor.info() != Eigen::Success) {
              throw std::runtime_error(
                "Native BAYES ETA covariance factorization failed.");
            }
            Matrix root(factor.matrixL());
            covariance_groups.push_back(covariance);
            root_groups.push_back(root);
            eta_roots[subject] = std::move(root);
          }
          cached_omega = parameters.omega;
          omega_factorizations += static_cast<int>(root_groups.size());
        } else if (!eta_independence) {
          ++omega_cache_hits;
        }
        Matrix eta_normals(
          static_cast<Eigen::Index>(tapes_.size()), n_eta_);
        // Match R's matrix(rnorm(), n_subjects, n_eta) column ordering.
        for (int effect = 0; effect < n_eta_; ++effect) {
          for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
            eta_normals(static_cast<Eigen::Index>(subject), effect) =
              R::rnorm(0.0, 1.0);
          }
        }
        std::vector<double> eta_log_uniforms(tapes_.size());
        for (double& value : eta_log_uniforms) {
          value = std::log(R::runif(0.0, 1.0));
        }
        std::ostringstream messages;
        if (fused_enabled_ && native_threads_ > 1 && tapes_.size() > 1U) {
          struct EtaCandidate {
            Vector eta;
            std::vector<double> point;
            double current_quad = 0.0;
            double candidate_quad = 0.0;
            double value = std::numeric_limits<double>::infinity();
          };
          std::vector<EtaCandidate> candidates(tapes_.size());
          for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
            EtaCandidate& candidate = candidates[subject];
            candidate.eta.resize(n_eta_);
            if (eta_independence) {
              Rcpp::NumericVector proposal_mode(n_eta_);
              for (int effect = 0; effect < n_eta_; ++effect) {
                proposal_mode[effect] = eta_proposal_modes(
                  static_cast<int>(subject), effect);
              }
              const Rcpp::NumericMatrix root = eta_proposal_roots[
                static_cast<int>(subject)];
              const Rcpp::NumericMatrix proposal_precision =
                eta_proposal_precisions[static_cast<int>(subject)];
              Rcpp::NumericVector current_eta(n_eta_);
              const double proposal_scale = eta_kernel == "student_t" ?
                std::sqrt(eta_df / R::rchisq(eta_df)) : 1.0;
              for (int row = 0; row < n_eta_; ++row) {
                current_eta[row] = eta(
                  static_cast<Eigen::Index>(subject), row);
                double increment = 0.0;
                for (int column = 0; column < n_eta_; ++column) {
                  increment += root(row, column) * eta_normals(
                    static_cast<Eigen::Index>(subject), column);
                }
                candidate.eta[row] = proposal_mode[row] +
                  proposal_scale * increment;
              }
              Rcpp::NumericVector candidate_eta_r =
                libertad::eigen_vector_to_r(candidate.eta);
              candidate.current_quad = quadratic(
                current_eta, proposal_mode, proposal_precision);
              candidate.candidate_quad = quadratic(
                candidate_eta_r, proposal_mode, proposal_precision);
            } else {
              candidate.eta = eta.row(
                static_cast<Eigen::Index>(subject)).transpose() +
                eta_step * eta_roots[subject] * eta_normals.row(
                  static_cast<Eigen::Index>(subject)).transpose();
            }
            candidate.point = points_[subject];
            fill_point_native(candidate.point, parameters, candidate.eta);
          }
          parallel_subjects(tapes_.size(), [&](std::size_t subject) {
            candidates[subject].value = direct_value_raw(
              subject, parameters, candidates[subject].eta);
          });
          evaluations_ += static_cast<long long>(tapes_.size());
          fused_evaluations_ += static_cast<long long>(tapes_.size());
          for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
            EtaCandidate& candidate = candidates[subject];
            ++attempted_eta;
            double log_ratio = -0.5 *
              (candidate.value - subject_values[subject]);
            if (eta_independence && std::isfinite(candidate.value)) {
              if (eta_kernel == "student_t") {
                log_ratio += 0.5 * (eta_df + static_cast<double>(n_eta_)) *
                  (std::log1p(candidate.candidate_quad / eta_df) -
                   std::log1p(candidate.current_quad / eta_df));
              } else {
                log_ratio += 0.5 *
                  (candidate.candidate_quad - candidate.current_quad);
              }
            }
            if (std::isfinite(candidate.value) &&
                eta_log_uniforms[subject] < log_ratio) {
              eta.row(static_cast<Eigen::Index>(subject)) =
                candidate.eta.transpose();
              points_[subject].swap(candidate.point);
              current -= 0.5 *
                (candidate.value - subject_values[subject]);
              subject_values[subject] = candidate.value;
              ++accepted_eta;
            }
          }
        } else for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
          Vector candidate_eta(n_eta_);
          double current_quad = 0.0;
          Rcpp::NumericVector proposal_mode;
          Rcpp::NumericMatrix proposal_precision;
          if (eta_independence) {
            proposal_mode = Rcpp::NumericVector(n_eta_);
            for (int effect = 0; effect < n_eta_; ++effect) {
              proposal_mode[effect] = eta_proposal_modes(
                static_cast<int>(subject), effect);
            }
            const Rcpp::NumericMatrix root = eta_proposal_roots[
              static_cast<int>(subject)];
            proposal_precision = Rcpp::NumericMatrix(
              eta_proposal_precisions[static_cast<int>(subject)]);
            Rcpp::NumericVector current_eta(n_eta_);
            const double proposal_scale = eta_kernel == "student_t" ?
              std::sqrt(eta_df / R::rchisq(eta_df)) : 1.0;
            for (int row = 0; row < n_eta_; ++row) {
              current_eta[row] = eta(static_cast<Eigen::Index>(subject), row);
              double increment = 0.0;
              for (int column = 0; column < n_eta_; ++column) {
                increment += root(row, column) * eta_normals(
                  static_cast<Eigen::Index>(subject), column);
              }
              candidate_eta[row] = proposal_mode[row] +
                proposal_scale * increment;
            }
            current_quad = quadratic(
              current_eta, proposal_mode, proposal_precision);
          } else {
            candidate_eta = eta.row(
              static_cast<Eigen::Index>(subject)).transpose() +
              eta_step * eta_roots[subject] * eta_normals.row(
                static_cast<Eigen::Index>(subject)).transpose();
          }
          std::vector<double> candidate_point = points_[subject];
          fill_point_native(
            candidate_point, parameters, candidate_eta);
          const double candidate_value = guarded_value(
            subject, parameters, candidate_eta, candidate_point,
            "native BAYES stochastic objective");
          ++attempted_eta;
          double log_ratio = -0.5 *
            (candidate_value - subject_values[subject]);
          if (eta_independence && std::isfinite(candidate_value)) {
            Rcpp::NumericVector candidate_eta_r =
              libertad::eigen_vector_to_r(candidate_eta);
            const double candidate_quad = quadratic(
              candidate_eta_r, proposal_mode, proposal_precision);
            if (eta_kernel == "student_t") {
              log_ratio += 0.5 * (eta_df + static_cast<double>(n_eta_)) *
                (std::log1p(candidate_quad / eta_df) -
                 std::log1p(current_quad / eta_df));
            } else {
              log_ratio += 0.5 * (candidate_quad - current_quad);
            }
          }
          if (std::isfinite(candidate_value) &&
              eta_log_uniforms[subject] < log_ratio) {
            eta.row(static_cast<Eigen::Index>(subject)) =
              candidate_eta.transpose();
            points_[subject].swap(candidate_point);
            current -= 0.5 * (candidate_value - subject_values[subject]);
            subject_values[subject] = candidate_value;
            ++accepted_eta;
          }
        }
        if (eta_independence) {
          const double acceptance =
            static_cast<double>(accepted_eta - accepted_eta_before) /
            static_cast<double>(std::max(
              attempted_eta - attempted_eta_before, 1));
          if (acceptance < eta_low_acceptance) {
            eta_force_refresh = true;
            ++eta_acceptance_refreshes;
          }
        }
      }

      if (dimension > 0 && outer_kernel == "adaptive_metropolis") {
        ++adaptive_n;
        Vector adaptive_value(dimension);
        for (int index = 0; index < dimension; ++index) {
          adaptive_value[index] = outer[
            random_positions[static_cast<std::size_t>(index)]];
        }
        const Vector delta = adaptive_value - adaptive_mean;
        adaptive_mean += delta / static_cast<double>(adaptive_n);
        adaptive_m2 += delta * (adaptive_value - adaptive_mean).transpose();
        if (adapt && iteration <= n_burn) {
          const double gain = std::min(
            0.02, std::pow(static_cast<double>(adaptive_n + 10), -0.6));
          log_multiplier += gain *
            ((accepted_outer_iteration ? 1.0 : 0.0) - target_acceptance);
          if (adaptive_n >= adaptive_start &&
              adaptive_n % adaptive_interval == 0) {
            const Matrix empirical = adaptive_m2 /
              static_cast<double>(std::max(adaptive_n - 1, 1));
            const double optimal = 2.38 * 2.38 /
              static_cast<double>(std::max(dimension, 1));
            const double ridge = step_scale * step_scale * 1e-3;
            Matrix candidate = std::exp(2.0 * log_multiplier) *
              (optimal * empirical +
               Matrix::Identity(dimension, dimension) * ridge);
            const double jitter = regularize_curvature(
              candidate, "adaptive native BAYES population proposal");
            covariance_regularizations += jitter > 0.0 ? 1 : 0;
            Eigen::LLT<Matrix> factor(candidate);
            if (factor.info() != Eigen::Success) {
              throw std::runtime_error(
                "Adaptive native BAYES proposal factorization failed.");
            }
            proposal_covariance = candidate;
            proposal_root = Matrix(factor.matrixL());
            ++covariance_updates;
          }
        }
      } else if (dimension > 0 && adapt && iteration <= n_burn &&
                 iteration % 50 == 0) {
        const double rate = static_cast<double>(accepted_outer) /
          static_cast<double>(std::max(attempted_outer, 1));
        step_scale *= std::exp(rate > 0.3 ? 0.1 : -0.1);
        proposal_root = Matrix::Identity(dimension, dimension) * step_scale;
        proposal_covariance = proposal_root * proposal_root.transpose();
      }

      if (iteration > n_burn && (iteration - n_burn) % n_thin == 0) {
        int column = 0;
        for (double value : parameters.theta) chain(keep, column++) = value;
        for (double value : parameters.sigma) chain(keep, column++) = value;
        for (double value : parameters.omega) chain(keep, column++) = value;
        for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
          for (int effect = 0; effect < n_eta_; ++effect) {
            chain(keep, column++) = eta(
              static_cast<Eigen::Index>(subject), effect);
          }
        }
        chain(keep, column) = current;
        ++keep;
      }
      if (iteration % 50 == 0) Rcpp::checkUserInterrupt();
    }
    ++calls_;
    Rcpp::NumericVector final_theta = Rcpp::wrap(parameters.theta);
    Rcpp::NumericVector final_sigma = Rcpp::wrap(parameters.sigma);
    Rcpp::NumericVector final_omega = Rcpp::wrap(parameters.omega);
    return Rcpp::List::create(
      Rcpp::Named("chain") = chain,
      Rcpp::Named("final_theta") = final_theta,
      Rcpp::Named("final_sigma") = final_sigma,
      Rcpp::Named("final_omega") = final_omega,
      Rcpp::Named("final_eta") = libertad::eigen_matrix_to_r(eta),
      Rcpp::Named("final_log_posterior") = current,
      Rcpp::Named("outer_acceptance") = static_cast<double>(accepted_outer) /
        static_cast<double>(std::max(attempted_outer, 1)),
      Rcpp::Named("eta_acceptance") = static_cast<double>(accepted_eta) /
        static_cast<double>(std::max(attempted_eta, 1)),
      Rcpp::Named("accepted_outer") = accepted_outer,
      Rcpp::Named("attempted_outer") = attempted_outer,
      Rcpp::Named("accepted_delayed") = accepted_delayed,
      Rcpp::Named("attempted_delayed") = attempted_delayed,
      Rcpp::Named("delayed_rejection_scale") = delayed_rejection_scale,
      Rcpp::Named("accepted_eta") = accepted_eta,
      Rcpp::Named("attempted_eta") = attempted_eta,
      Rcpp::Named("accepted_mu") = accepted_mu,
      Rcpp::Named("attempted_mu") = attempted_mu,
      Rcpp::Named("mu_acceptance") = static_cast<double>(accepted_mu) /
        static_cast<double>(std::max(attempted_mu, 1)),
      Rcpp::Named("gibbs_omega_parameters") =
        static_cast<int>(gibbs_omegas.size()),
      Rcpp::Named("gibbs_omega_updates") = gibbs_omega_updates,
      Rcpp::Named("gibbs_omega_draws") = gibbs_omega_draws,
      Rcpp::Named("final_step_scale") = outer_kernel == "adaptive_metropolis" ?
        std::exp(log_multiplier) : step_scale,
      Rcpp::Named("covariance_updates") = covariance_updates,
      Rcpp::Named("covariance_regularizations") = covariance_regularizations,
      Rcpp::Named("covariance") =
        libertad::eigen_matrix_to_r(proposal_covariance),
      Rcpp::Named("multiplier") = std::exp(log_multiplier),
      Rcpp::Named("omega_factorizations") = omega_factorizations,
      Rcpp::Named("omega_cache_hits") = omega_cache_hits,
      Rcpp::Named("eta_kernel") = eta_kernel,
      Rcpp::Named("eta_refreshes") = eta_refreshes,
      Rcpp::Named("eta_refresh_failures") = eta_refresh_failures,
      Rcpp::Named("eta_rescue_iterations") = eta_rescue_iterations,
      Rcpp::Named("eta_parameter_refreshes") = eta_parameter_refreshes,
      Rcpp::Named("eta_acceptance_refreshes") = eta_acceptance_refreshes,
      Rcpp::Named("eta_fallback_iterations") = eta_fallback_iterations,
      Rcpp::Named("eta_last_error") = eta_last_error,
      Rcpp::Named("backend") = "persistent-native-cpp-bayes-coordinator");
  }

  Rcpp::List telemetry() const {
    return Rcpp::List::create(
      Rcpp::Named("subjects") = static_cast<int>(tapes_.size()),
      Rcpp::Named("evaluations") = static_cast<double>(evaluations_),
      Rcpp::Named("calls") = calls_,
      Rcpp::Named("laplace_refreshes") = laplace_refreshes_,
      Rcpp::Named("laplace_mode_evaluations") =
        static_cast<double>(laplace_mode_evaluations_),
      Rcpp::Named("tape_records") = tape_records_,
      Rcpp::Named("tape_retapes") = tape_retapes_,
      Rcpp::Named("ode_owned_tapes") = use_ode_,
      Rcpp::Named("fused_advan_requested") = fused_requested_,
      Rcpp::Named("fused_advan_enabled") = fused_enabled_,
      Rcpp::Named("fused_advan_evaluations") =
        static_cast<double>(fused_evaluations_),
      Rcpp::Named("recorded_tape_evaluations") =
        static_cast<double>(recorded_evaluations_),
      Rcpp::Named("fused_advan_guard_failures") = fused_guard_failures_,
      Rcpp::Named("fused_advan_fallback_reason") = fused_fallback_reason_,
      Rcpp::Named("native_subject_threads") = native_threads_,
      Rcpp::Named("persistent_worker_pool") = !is_null_pool(),
      Rcpp::Named("worker_pool_dispatches") = subject_pool_ ?
        static_cast<double>(subject_pool_->dispatches()) : 0.0);
  }

 private:
  ModelEngine* engine_ = nullptr;
  int n_theta_ = 0;
  int n_eta_ = 0;
  int n_sigma_ = 0;
  int n_omega_ = 0;
  std::size_t domain_ = 0U;
  bool use_ode_ = false;
  double guard_radius_ = 0.5;
  SEXP retained_engine_ = R_NilValue;
  SEXP retained_subject_data_ = R_NilValue;
  std::vector<SEXP> subject_data_;
  std::vector<std::shared_ptr<const OwnedEventTable>> owned_subject_data_;
  std::vector<ObjectiveTape*> tapes_;
  std::vector<std::unique_ptr<ObjectiveTape>> owned_tapes_;
  std::vector<std::vector<double>> points_;
  std::vector<std::vector<double>> anchors_;
  std::vector<std::size_t> eta_positions_;
  long long evaluations_ = 0;
  int calls_ = 0;
  int laplace_refreshes_ = 0;
  long long laplace_mode_evaluations_ = 0;
  int tape_records_ = 0;
  int tape_retapes_ = 0;
  bool fused_requested_ = false;
  bool fused_enabled_ = false;
  int native_threads_ = 1;
  std::unique_ptr<NativeSubjectPool> subject_pool_;
  std::vector<bool> fused_subject_;
  long long fused_evaluations_ = 0;
  long long recorded_evaluations_ = 0;
  int fused_guard_failures_ = 0;
  std::string fused_fallback_reason_;

  bool is_null_pool() const { return subject_pool_ == nullptr; }

  template <class Function>
  void parallel_subjects(std::size_t count, Function function) const {
    const int workers = std::min(
      native_threads_, static_cast<int>(std::max<std::size_t>(count, 1U)));
    if (workers <= 1 || count < 2U) {
      for (std::size_t subject = 0; subject < count; ++subject) {
        function(subject);
      }
      return;
    }
    if (!subject_pool_) {
      throw std::logic_error("The native subject worker pool is unavailable.");
    }
    subject_pool_->run(count, std::move(function));
  }

  double direct_value_raw(
      std::size_t subject, const StochasticBayesParameters& parameters,
      const Vector& eta) const {
    if (!fused_enabled_ &&
        (subject >= fused_subject_.size() || !fused_subject_[subject])) {
      // During construction fused_enabled_ is not final yet, so the per-subject
      // eligibility flag is intentionally not required until after the guard.
      if (subject >= owned_subject_data_.size() ||
          !owned_subject_data_[subject]) {
        throw std::logic_error("A fused ADVAN subject has no native event data.");
      }
    }
    if (!engine_ || subject >= owned_subject_data_.size() ||
        !owned_subject_data_[subject]) {
      throw std::logic_error("Fused ADVAN likelihood evaluation is unavailable.");
    }
    const EventDataView data(owned_subject_data_[subject]);
    std::vector<double> eta_values(static_cast<std::size_t>(eta.size()));
    for (Eigen::Index effect = 0; effect < eta.size(); ++effect) {
      eta_values[static_cast<std::size_t>(effect)] = eta[effect];
    }
    return population_joint_nll_t<double>(
      *engine_, data, parameters.theta, eta_values, parameters.sigma,
      parameters.omega, true);
  }

  std::vector<double> native_point(
      const StochasticBayesParameters& parameters, const Vector& eta) const {
    std::vector<double> point(domain_, 0.0);
    fill_point_native(point, parameters, eta);
    return point;
  }

  bool material_movement(
      std::size_t subject, const StochasticBayesParameters& parameters,
      const Vector& eta) const {
    if (!use_ode_) return false;
    const std::vector<double> point = native_point(parameters, eta);
    const std::vector<double>& anchor = anchors_[subject];
    if (point.size() != anchor.size()) return true;
    double distance = 0.0;
    for (std::size_t index = 0; index < point.size(); ++index) {
      distance = std::max(
        distance, std::abs(point[index] - anchor[index]) /
          std::max(std::abs(anchor[index]), 1.0));
    }
    return !std::isfinite(distance) || distance > guard_radius_;
  }

  void record_subject_tape(
      std::size_t subject, const StochasticBayesParameters& parameters,
      const Vector& eta, bool retape) {
    if (!engine_ || subject >= subject_data_.size()) {
      throw std::runtime_error(
        "A stochastic tape cannot be rebuilt without its model and data.");
    }
    const EventDataView data = event_data_view(subject_data_[subject]);
    Rcpp::NumericVector theta = Rcpp::wrap(parameters.theta);
    Rcpp::NumericVector sigma = Rcpp::wrap(parameters.sigma);
    Rcpp::NumericVector omega = Rcpp::wrap(parameters.omega);
    Rcpp::NumericMatrix eta_matrix(1, n_eta_);
    for (int effect = 0; effect < n_eta_; ++effect) {
      eta_matrix(0, effect) = eta[effect];
    }
    owned_tapes_[subject] = record_objective_tape(
      *engine_, data, theta, eta_matrix, sigma, omega, true);
    tapes_[subject] = owned_tapes_[subject].get();
    anchors_[subject] = native_point(parameters, eta);
    ++tape_records_;
    if (retape) ++tape_retapes_;
  }

  double guarded_value(
      std::size_t subject, const StochasticBayesParameters& parameters,
      const Vector& eta, std::vector<double>& point,
      const std::string& context) {
    if (fused_enabled_ && subject < fused_subject_.size() &&
        fused_subject_[subject]) {
      const double value = direct_value_raw(subject, parameters, eta);
      ++evaluations_;
      ++fused_evaluations_;
      return value;
    }
    if (material_movement(subject, parameters, eta)) {
      record_subject_tape(subject, parameters, eta, true);
    }
    for (int attempt = 0; attempt < 4; ++attempt) {
      try {
        std::ostringstream messages;
        const std::vector<double> value = tapes_[subject]->fun.Forward(
          0, point, messages);
        require_unchanged_path(tapes_[subject]->fun, context);
        ++evaluations_;
        ++recorded_evaluations_;
        return value.empty() ? std::numeric_limits<double>::infinity() :
          value[0];
      } catch (const TapePathChange&) {
        record_subject_tape(subject, parameters, eta, true);
      }
    }
    throw std::runtime_error(
      "A stochastic objective remained structurally unstable after retaping.");
  }

  Matrix subject_eta_covariance(
      const StochasticBayesMap& map,
      const StochasticBayesParameters& parameters,
      std::size_t subject) const {
    const Matrix base = map.omega_covariance(parameters);
    if (!engine_) {
      if (base.rows() != n_eta_) {
        throw std::invalid_argument(
          "Expanded random effects require a retained model engine.");
      }
      return base;
    }
    const EventDataView data = event_data_view(subject_data_[subject]);
    return expanded_omega_t<double>(*engine_, data, base, n_eta_);
  }

  static double gaussian_quadratic(
      const Vector& difference, const Matrix& root) {
    if (root.rows() != difference.size() || root.cols() != difference.size()) {
      return std::numeric_limits<double>::infinity();
    }
    const Vector standardized =
      root.triangularView<Eigen::Lower>().solve(difference);
    return standardized.allFinite() ? standardized.squaredNorm() :
      std::numeric_limits<double>::infinity();
  }

  static double log_one_minus_acceptance(double log_alpha) {
    if (log_alpha >= 0.0) {
      return -std::numeric_limits<double>::infinity();
    }
    if (!std::isfinite(log_alpha)) return 0.0;
    return log_alpha < -std::log(2.0) ?
      std::log1p(-std::exp(log_alpha)) : std::log(-std::expm1(log_alpha));
  }

  static double regularize_curvature(
      Matrix& matrix, const std::string& context) {
    matrix = 0.5 * (matrix + matrix.transpose()).eval();
    if (!matrix.allFinite()) {
      throw std::domain_error(context + " is not finite.");
    }
    const auto eigen = libertad::detail::self_adjoint_eigen(matrix, false);
    if (eigen.info != Eigen::Success || !eigen.values.allFinite()) {
      throw std::runtime_error(context + " decomposition failed.");
    }
    const double largest = std::max(
      eigen.values.cwiseAbs().maxCoeff(), 1.0);
    const double jitter = std::max(
      0.0, largest * 1e-9 - eigen.values.minCoeff());
    if (jitter > largest * 1e-2) {
      throw std::domain_error(
        context + " is not sufficiently positive definite.");
    }
    matrix.diagonal().array() += jitter;
    return jitter;
  }

  void validate_parameters(
      const Rcpp::NumericVector& theta, const Rcpp::NumericMatrix& eta,
      const Rcpp::NumericVector& sigma,
      const Rcpp::NumericVector& omega) const {
    if (theta.size() != n_theta_ || eta.nrow() != static_cast<int>(tapes_.size()) ||
        eta.ncol() != n_eta_ || sigma.size() != n_sigma_ ||
        omega.size() != n_omega_) {
      throw std::invalid_argument(
        "Persistent stochastic parameter dimensions changed.");
    }
    for (double value : theta) if (!std::isfinite(value)) {
      throw std::invalid_argument("Persistent stochastic THETAs must be finite.");
    }
    for (double value : eta) if (!std::isfinite(value)) {
      throw std::invalid_argument("Persistent stochastic ETAs must be finite.");
    }
    for (double value : sigma) if (!std::isfinite(value)) {
      throw std::invalid_argument("Persistent stochastic SIGMAs must be finite.");
    }
    for (double value : omega) if (!std::isfinite(value)) {
      throw std::invalid_argument("Persistent stochastic OMEGAs must be finite.");
    }
  }

  void validate_sampler(
      const Rcpp::NumericVector& theta, const Rcpp::NumericMatrix& eta,
      const Rcpp::NumericVector& sigma, const Rcpp::NumericVector& omega,
      const Rcpp::List& roots, const Rcpp::NumericMatrix& normals,
      const Rcpp::NumericVector& uniforms, int steps) const {
    validate_parameters(theta, eta, sigma, omega);
    if (roots.size() != static_cast<int>(tapes_.size()) || steps < 1 ||
        normals.nrow() != static_cast<int>(tapes_.size()) * steps ||
        normals.ncol() != n_eta_ || uniforms.size() != normals.nrow()) {
      throw std::invalid_argument(
        "Persistent stochastic Metropolis inputs are inconsistent.");
    }
  }

  void fill_point(
      std::size_t subject, const Rcpp::NumericVector& theta,
      const Rcpp::NumericMatrix& eta, const Rcpp::NumericVector& sigma,
      const Rcpp::NumericVector& omega) {
    std::vector<double>& point = points_[subject];
    std::copy(theta.begin(), theta.end(), point.begin());
    for (int effect = 0; effect < n_eta_; ++effect) {
      point[static_cast<std::size_t>(n_theta_ + effect)] =
        eta(static_cast<int>(subject), effect);
    }
    std::copy(
      sigma.begin(), sigma.end(), point.begin() + n_theta_ + n_eta_);
    std::copy(
      omega.begin(), omega.end(),
      point.begin() + n_theta_ + n_eta_ + n_sigma_);
  }

  void fill_point_native(
      std::vector<double>& point,
      const StochasticBayesParameters& parameters,
      const Vector& eta) const {
    std::copy(parameters.theta.begin(), parameters.theta.end(), point.begin());
    for (int effect = 0; effect < n_eta_; ++effect) {
      point[static_cast<std::size_t>(n_theta_ + effect)] = eta[effect];
    }
    std::copy(
      parameters.sigma.begin(), parameters.sigma.end(),
      point.begin() + n_theta_ + n_eta_);
    std::copy(
      parameters.omega.begin(), parameters.omega.end(),
      point.begin() + n_theta_ + n_eta_ + n_sigma_);
  }

  double bayes_log_posterior(
      const StochasticBayesMap& map,
      const StochasticBayesParameters& parameters,
      const Matrix& eta, std::vector<double>& subject_values) {
    const double prior = map.log_prior(parameters);
    if (!std::isfinite(prior)) {
      return -std::numeric_limits<double>::infinity();
    }
    subject_values.assign(tapes_.size(), 0.0);
    if (fused_enabled_) {
      parallel_subjects(tapes_.size(), [&](std::size_t subject) {
        subject_values[subject] = direct_value_raw(
          subject, parameters,
          eta.row(static_cast<Eigen::Index>(subject)).transpose());
      });
      evaluations_ += static_cast<long long>(tapes_.size());
      fused_evaluations_ += static_cast<long long>(tapes_.size());
      double total = 0.0;
      for (double value : subject_values) {
        if (!std::isfinite(value)) {
          return -std::numeric_limits<double>::infinity();
        }
        total += value;
      }
      return -0.5 * total + prior + parameters.log_jacobian;
    }
    double total = 0.0;
    for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
      Vector current_eta = eta.row(
        static_cast<Eigen::Index>(subject)).transpose();
      fill_point_native(points_[subject], parameters, current_eta);
      const double value = guarded_value(
        subject, parameters, current_eta, points_[subject],
        "native BAYES population objective");
      if (!std::isfinite(value)) {
        return -std::numeric_limits<double>::infinity();
      }
      subject_values[subject] = value;
      total += value;
    }
    return -0.5 * total + prior + parameters.log_jacobian;
  }

  struct MuSystem {
    Matrix hessian;
    Vector mean;
  };

  MuSystem bayes_mu_system(
      const StochasticMuConfig& mu, const StochasticBayesMap& map,
      const StochasticBayesParameters& parameters, const Matrix& eta) const {
    const int p = static_cast<int>(mu.theta.size());
    Matrix hessian = Matrix::Zero(p, p);
    Vector score = Vector::Zero(p);
    const Vector beta = mu.beta(parameters);
    for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
      Matrix covariance = subject_eta_covariance(
        map, parameters, subject);
      regularize_curvature(covariance, "native BAYES MU covariance");
      Eigen::LLT<Matrix> factor(covariance);
      if (factor.info() != Eigen::Success) {
        throw std::runtime_error("Native BAYES MU covariance factorization failed.");
      }
      const Matrix precision = factor.solve(
        Matrix::Identity(n_eta_, n_eta_));
      Matrix design(n_eta_, p);
      for (int column = 0; column < p; ++column) {
        design.col(column) = mu.design_columns[static_cast<std::size_t>(column)]
          .row(static_cast<Eigen::Index>(subject)).transpose();
      }
      const Vector centered = eta.row(
        static_cast<Eigen::Index>(subject)).transpose() + design * beta;
      hessian.noalias() += design.transpose() * precision * design;
      score.noalias() += design.transpose() * precision * centered;
    }
    regularize_curvature(hessian, "native BAYES MU information");
    Eigen::LLT<Matrix> factor(hessian);
    if (factor.info() != Eigen::Success) {
      throw std::runtime_error("Native BAYES MU information factorization failed.");
    }
    return MuSystem{hessian, factor.solve(score)};
  }

  static double bayes_mu_log_proposal(
      const StochasticMuConfig& mu,
      const StochasticBayesParameters& parameters, const Vector& beta,
      const MuSystem& system) {
    const Vector difference = beta - system.mean;
    Eigen::LLT<Matrix> factor(system.hessian);
    if (factor.info() != Eigen::Success) {
      return -std::numeric_limits<double>::infinity();
    }
    const Matrix lower = Matrix(factor.matrixL());
    const double logdet = 2.0 * lower.diagonal().array().log().sum();
    return 0.5 * logdet -
      0.5 * static_cast<double>(beta.size()) *
        std::log(2.0 * std::acos(-1.0)) -
      0.5 * difference.dot(system.hessian * difference) +
      mu.log_native_jacobian(parameters);
  }

  bool bayes_mu_step(
      const StochasticMuConfig& mu, const StochasticBayesMap& map,
      Vector& outer, StochasticBayesParameters& parameters, Matrix& eta,
      std::vector<double>& subject_values, double& current) {
    if (!mu.active || mu.theta.empty()) return false;
    const MuSystem system = bayes_mu_system(mu, map, parameters, eta);
    Eigen::LLT<Matrix> factor(system.hessian);
    if (factor.info() != Eigen::Success) return false;
    Vector normal(system.mean.size());
    for (Eigen::Index index = 0; index < normal.size(); ++index) {
      normal[index] = R::rnorm(0.0, 1.0);
    }
    const Matrix upper = Matrix(factor.matrixU());
    const Vector candidate_beta = system.mean +
      upper.triangularView<Eigen::Upper>().solve(normal);
    StochasticBayesParameters candidate_parameters = parameters;
    mu.set_beta(candidate_parameters, candidate_beta);
    Vector candidate_outer;
    try {
      candidate_outer = map.encode(candidate_parameters);
    } catch (const std::exception&) {
      return false;
    }
    const Vector current_beta = mu.beta(parameters);
    Matrix candidate_eta = mu.recenter(eta, current_beta, candidate_beta);
    std::vector<double> candidate_values;
    const double candidate_logp = bayes_log_posterior(
      map, candidate_parameters, candidate_eta, candidate_values);
    const double log_ratio = candidate_logp - current +
      bayes_mu_log_proposal(mu, parameters, current_beta, system) -
      bayes_mu_log_proposal(
        mu, candidate_parameters, candidate_beta, system);
    if (std::isfinite(candidate_logp) &&
        std::log(R::runif(0.0, 1.0)) < log_ratio) {
      outer.swap(candidate_outer);
      parameters = std::move(candidate_parameters);
      eta.swap(candidate_eta);
      subject_values.swap(candidate_values);
      current = candidate_logp;
      return true;
    }
    return false;
  }

  static double quadratic(
      const Rcpp::NumericVector& value, const Rcpp::NumericVector& center,
      const Rcpp::NumericMatrix& precision) {
    const int dimension = value.size();
    double result = 0.0;
    for (int row = 0; row < dimension; ++row) {
      const double left = value[row] - center[row];
      for (int column = 0; column < dimension; ++column) {
        result += left * precision(row, column) *
          (value[column] - center[column]);
      }
    }
    return result;
  }

  Rcpp::List metropolis(
      const Rcpp::NumericVector& theta, const Rcpp::NumericMatrix& eta_input,
      const Rcpp::NumericVector& sigma, const Rcpp::NumericVector& omega,
      const Rcpp::List& roots, SEXP modes_input, SEXP precisions_input,
      const Rcpp::NumericMatrix& normals,
      const Rcpp::NumericVector& log_uniforms, int steps, double scale,
      Rcpp::Nullable<Rcpp::NumericVector> current_values_input,
      bool independent, double proposal_df,
      Rcpp::Nullable<Rcpp::NumericVector> proposal_scales_input) {
    const bool use_current = current_values_input.isNotNull();
    Rcpp::NumericVector supplied;
    if (use_current) {
      supplied = Rcpp::NumericVector(current_values_input);
      if (supplied.size() != static_cast<int>(tapes_.size())) {
        throw std::invalid_argument(
          "Cached stochastic values must match the number of subjects.");
      }
    }
    Rcpp::NumericMatrix modes;
    Rcpp::List precisions;
    const bool student_t = independent && std::isfinite(proposal_df);
    Rcpp::NumericVector proposal_scales;
    if (independent) {
      modes = Rcpp::NumericMatrix(modes_input);
      precisions = Rcpp::List(precisions_input);
      if (student_t) {
        if (proposal_df <= 2.0 || proposal_scales_input.isNull()) {
          throw std::invalid_argument(
            "Student-t Laplace proposals require degrees of freedom and scales.");
        }
        proposal_scales = Rcpp::NumericVector(proposal_scales_input);
        if (proposal_scales.size() !=
            static_cast<int>(tapes_.size()) * steps) {
          throw std::invalid_argument(
            "Student-t Laplace proposal scales have invalid dimensions.");
        }
      }
    }
    StochasticBayesParameters native_parameters;
    native_parameters.theta = Rcpp::as<std::vector<double>>(theta);
    native_parameters.sigma = Rcpp::as<std::vector<double>>(sigma);
    native_parameters.omega = Rcpp::as<std::vector<double>>(omega);
    if (fused_enabled_) {
      const std::size_t subjects = tapes_.size();
      Matrix eta_native(eta_input.nrow(), eta_input.ncol());
      for (int subject = 0; subject < eta_input.nrow(); ++subject) {
        for (int effect = 0; effect < eta_input.ncol(); ++effect) {
          eta_native(subject, effect) = eta_input(subject, effect);
        }
      }
      std::vector<Matrix> roots_native(subjects);
      std::vector<Vector> modes_native(subjects);
      std::vector<Matrix> precisions_native(subjects);
      for (std::size_t subject = 0; subject < subjects; ++subject) {
        const Rcpp::NumericMatrix root = roots[static_cast<int>(subject)];
        if (root.nrow() != n_eta_ || root.ncol() != n_eta_) {
          throw std::invalid_argument(
            "A fused stochastic proposal root has invalid dimensions.");
        }
        roots_native[subject].resize(n_eta_, n_eta_);
        for (int row = 0; row < n_eta_; ++row) {
          for (int column = 0; column < n_eta_; ++column) {
            roots_native[subject](row, column) = root(row, column);
          }
        }
        if (independent) {
          modes_native[subject].resize(n_eta_);
          for (int effect = 0; effect < n_eta_; ++effect) {
            modes_native[subject][effect] = modes(
              static_cast<int>(subject), effect);
          }
          const Rcpp::NumericMatrix precision =
            precisions[static_cast<int>(subject)];
          if (precision.nrow() != n_eta_ || precision.ncol() != n_eta_) {
            throw std::invalid_argument(
              "A fused stochastic proposal precision has invalid dimensions.");
          }
          precisions_native[subject].resize(n_eta_, n_eta_);
          for (int row = 0; row < n_eta_; ++row) {
            for (int column = 0; column < n_eta_; ++column) {
              precisions_native[subject](row, column) = precision(row, column);
            }
          }
        }
      }
      Matrix normals_native(normals.nrow(), normals.ncol());
      for (int row = 0; row < normals.nrow(); ++row) {
        for (int column = 0; column < normals.ncol(); ++column) {
          normals_native(row, column) = normals(row, column);
        }
      }
      const std::vector<double> uniforms_native =
        Rcpp::as<std::vector<double>>(log_uniforms);
      const std::vector<double> supplied_native = use_current ?
        Rcpp::as<std::vector<double>>(supplied) : std::vector<double>();
      const std::vector<double> scales_native = student_t ?
        Rcpp::as<std::vector<double>>(proposal_scales) :
        std::vector<double>(subjects * static_cast<std::size_t>(steps), 1.0);
      struct FusedResult {
        Vector eta;
        std::vector<double> point;
        double value = std::numeric_limits<double>::infinity();
        int accepted = 0;
      };
      std::vector<FusedResult> results(subjects);
      parallel_subjects(subjects, [&](std::size_t subject) {
        FusedResult& result = results[subject];
        result.eta = eta_native.row(
          static_cast<Eigen::Index>(subject)).transpose();
        result.point = points_[subject];
        double current_value = use_current ? supplied_native[subject] :
          direct_value_raw(subject, native_parameters, result.eta);
        if (!std::isfinite(current_value)) {
          throw std::domain_error(
            "Current fused stochastic ETA objective is not finite.");
        }
        double current_quad = independent ?
          (result.eta - modes_native[subject]).dot(
            precisions_native[subject] *
            (result.eta - modes_native[subject])) : 0.0;
        for (int step = 0; step < steps; ++step) {
          const int draw = static_cast<int>(subject) * steps + step;
          Vector candidate_eta = independent ? modes_native[subject] : result.eta;
          const Vector increment = roots_native[subject] *
            normals_native.row(draw).transpose();
          if (independent) {
            candidate_eta +=
              scales_native[static_cast<std::size_t>(draw)] * increment;
          } else {
            candidate_eta += scale * increment;
          }
          const double candidate_value = direct_value_raw(
            subject, native_parameters, candidate_eta);
          double log_ratio = std::isfinite(candidate_value) ?
            -0.5 * (candidate_value - current_value) :
            -std::numeric_limits<double>::infinity();
          double candidate_quad = 0.0;
          if (independent && std::isfinite(candidate_value)) {
            const Vector difference = candidate_eta - modes_native[subject];
            candidate_quad = difference.dot(
              precisions_native[subject] * difference);
            if (student_t) {
              log_ratio += 0.5 *
                (proposal_df + static_cast<double>(n_eta_)) *
                (std::log1p(candidate_quad / proposal_df) -
                 std::log1p(current_quad / proposal_df));
            } else {
              log_ratio += 0.5 * (candidate_quad - current_quad);
            }
          }
          if (uniforms_native[static_cast<std::size_t>(draw)] < log_ratio) {
            result.eta = std::move(candidate_eta);
            current_value = candidate_value;
            current_quad = candidate_quad;
            ++result.accepted;
          }
        }
        fill_point_native(result.point, native_parameters, result.eta);
        result.value = current_value;
      });
      Rcpp::NumericMatrix eta = Rcpp::clone(eta_input);
      Rcpp::NumericVector values(static_cast<R_xlen_t>(subjects));
      int accepted = 0;
      for (std::size_t subject = 0; subject < subjects; ++subject) {
        for (int effect = 0; effect < n_eta_; ++effect) {
          eta(static_cast<int>(subject), effect) = results[subject].eta[effect];
        }
        points_[subject].swap(results[subject].point);
        values[static_cast<R_xlen_t>(subject)] = results[subject].value;
        accepted += results[subject].accepted;
      }
      const long long candidates = static_cast<long long>(subjects) * steps;
      const long long currents = use_current ? 0LL :
        static_cast<long long>(subjects);
      evaluations_ += candidates + currents;
      fused_evaluations_ += candidates + currents;
      ++calls_;
      return Rcpp::List::create(
        Rcpp::Named("eta") = eta,
        Rcpp::Named("value") = values,
        Rcpp::Named("accepted") = accepted,
        Rcpp::Named("attempted") = static_cast<int>(subjects) * steps,
        Rcpp::Named("current_evaluations") = static_cast<int>(currents),
        Rcpp::Named("current_cache_hits") = use_current ?
          static_cast<int>(subjects) : 0,
        Rcpp::Named("candidate_evaluations") =
          static_cast<int>(candidates),
        Rcpp::Named("kernel") = independent ?
          (student_t ? "fused-laplace-student-t-independence" :
           "fused-laplace-independence") : "fused-random-walk");
    }
    Rcpp::NumericMatrix eta = Rcpp::clone(eta_input);
    Rcpp::NumericVector values(static_cast<R_xlen_t>(tapes_.size()));
    int accepted = 0;
    int current_evaluations = 0;
    int current_cache_hits = 0;
    int candidate_evaluations = 0;
    std::ostringstream messages;
    for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
      fill_point(subject, theta, eta, sigma, omega);
      EtaEvaluation current;
      if (use_current) {
        current.value = supplied[static_cast<R_xlen_t>(subject)];
        current.finite = std::isfinite(current.value);
        ++current_cache_hits;
      } else {
        Vector subject_eta(n_eta_);
        for (int effect = 0; effect < n_eta_; ++effect) {
          subject_eta[effect] = eta(static_cast<int>(subject), effect);
        }
        current.value = guarded_value(
          subject, native_parameters, subject_eta, points_[subject],
          "persistent stochastic current ETA objective");
        current.finite = std::isfinite(current.value);
        ++current_evaluations;
      }
      if (!current.finite) {
        throw std::domain_error("Current stochastic ETA objective is not finite.");
      }
      Rcpp::NumericMatrix root = roots[static_cast<int>(subject)];
      if (root.nrow() != n_eta_ || root.ncol() != n_eta_) {
        throw std::invalid_argument(
          "A persistent stochastic proposal root has invalid dimensions.");
      }
      Rcpp::NumericVector mode;
      Rcpp::NumericMatrix precision;
      Rcpp::NumericVector current_eta(n_eta_);
      for (int effect = 0; effect < n_eta_; ++effect) {
        current_eta[effect] = eta(static_cast<int>(subject), effect);
      }
      double current_quad = 0.0;
      if (independent) {
        mode = Rcpp::NumericVector(n_eta_);
        for (int effect = 0; effect < n_eta_; ++effect) {
          mode[effect] = modes(static_cast<int>(subject), effect);
        }
        precision = Rcpp::NumericMatrix(precisions[static_cast<int>(subject)]);
        if (precision.nrow() != n_eta_ || precision.ncol() != n_eta_) {
          throw std::invalid_argument(
            "A Laplace proposal precision has invalid dimensions.");
        }
        current_quad = quadratic(current_eta, mode, precision);
      }
      for (int step = 0; step < steps; ++step) {
        const int draw = static_cast<int>(subject) * steps + step;
        Rcpp::NumericVector candidate_eta(n_eta_);
        std::vector<double> candidate_point = points_[subject];
        for (int row = 0; row < n_eta_; ++row) {
          double increment = 0.0;
          for (int column = 0; column < n_eta_; ++column) {
            increment += root(row, column) * normals(draw, column);
          }
          const double proposal_scale = student_t ? proposal_scales[draw] : 1.0;
          candidate_eta[row] = independent ?
            mode[row] + proposal_scale * increment :
            eta(static_cast<int>(subject), row) + scale * increment;
          candidate_point[eta_positions_[static_cast<std::size_t>(row)]] =
            candidate_eta[row];
        }
        Vector candidate_eta_eigen(n_eta_);
        for (int effect = 0; effect < n_eta_; ++effect) {
          candidate_eta_eigen[effect] = candidate_eta[effect];
        }
        EtaEvaluation candidate;
        candidate.value = guarded_value(
          subject, native_parameters, candidate_eta_eigen, candidate_point,
          "persistent stochastic candidate ETA objective");
        candidate.finite = std::isfinite(candidate.value);
        ++candidate_evaluations;
        double log_ratio = candidate.finite ?
          -0.5 * (candidate.value - current.value) :
          -std::numeric_limits<double>::infinity();
        double candidate_quad = 0.0;
        if (independent && candidate.finite) {
          candidate_quad = quadratic(candidate_eta, mode, precision);
          if (student_t) {
            log_ratio += 0.5 * (proposal_df + static_cast<double>(n_eta_)) *
              (std::log1p(candidate_quad / proposal_df) -
               std::log1p(current_quad / proposal_df));
          } else {
            log_ratio += 0.5 * (candidate_quad - current_quad);
          }
        }
        if (log_uniforms[draw] < log_ratio) {
          for (int effect = 0; effect < n_eta_; ++effect) {
            eta(static_cast<int>(subject), effect) = candidate_eta[effect];
            current_eta[effect] = candidate_eta[effect];
          }
          points_[subject].swap(candidate_point);
          current = std::move(candidate);
          current_quad = candidate_quad;
          ++accepted;
        }
      }
      values[static_cast<R_xlen_t>(subject)] = current.value;
      if ((subject + 1U) % 64U == 0U) Rcpp::checkUserInterrupt();
    }
    ++calls_;
    return Rcpp::List::create(
      Rcpp::Named("eta") = eta,
      Rcpp::Named("value") = values,
      Rcpp::Named("accepted") = accepted,
      Rcpp::Named("attempted") = static_cast<int>(tapes_.size()) * steps,
      Rcpp::Named("current_evaluations") = current_evaluations,
      Rcpp::Named("current_cache_hits") = current_cache_hits,
      Rcpp::Named("candidate_evaluations") = candidate_evaluations,
      Rcpp::Named("kernel") = independent ?
        (student_t ? "laplace-student-t-independence" :
         "laplace-independence") : "random-walk");
  }
};

// Persistent native coordinator for the LibeR-optimized Gaussian-quadrature
// estimator.  The NONMEM-compatible policy deliberately keeps its established
// R/L-BFGS-B coordinator; this class removes R callbacks from adaptive proposal
// construction, signed-grid evaluation, score search, and exact finite-grid
// refinement without changing the quadrature objective.
