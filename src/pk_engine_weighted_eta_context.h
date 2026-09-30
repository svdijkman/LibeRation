// Internal implementation fragment included inside namespace liberation by
// pk_engine.cpp. It intentionally has no include guard or namespace declaration.
// Contains the persistent weighted-ETA context shared by ITS, IMP, and SAEM.

class WeightedEtaCollection {
 public:
  WeightedEtaCollection(
      SEXP engine_pointer, const Rcpp::List& tape_pointers,
      const Rcpp::List& subject_data,
      int n_theta, int n_eta, int n_sigma, int n_omega,
      bool use_ode, bool reduced_population_tape, int native_threads,
      int ode_support_tape_limit)
      : n_theta_(n_theta), n_eta_(n_eta), n_sigma_(n_sigma),
        n_omega_(n_omega), retained_subject_data_(subject_data),
        use_ode_(use_ode),
        native_threads_(std::max(1, native_threads)),
        reduced_requested_(reduced_population_tape),
        ode_support_tape_limit_(ode_support_tape_limit) {
    if (tape_pointers.size() < 1 || subject_data.size() != tape_pointers.size() ||
        n_theta < 0 || n_eta < 0 || n_sigma < 0 || n_omega < 0 ||
        native_threads < 1 || ode_support_tape_limit < 1) {
      throw std::invalid_argument(
        "Persistent weighted-ETA inputs are inconsistent.");
    }
    domain_ = static_cast<std::size_t>(
      n_theta_ + n_eta_ + n_sigma_ + n_omega_);
    Rcpp::XPtr<ModelEngine> engine(engine_pointer);
    engine_ = engine.get();
    retained_engine_ = engine_pointer;
    R_PreserveObject(retained_engine_);
    tapes_.reserve(static_cast<std::size_t>(tape_pointers.size()));
    subject_data_.reserve(static_cast<std::size_t>(subject_data.size()));
    points_.resize(static_cast<std::size_t>(tape_pointers.size()));
    grids_.resize(static_cast<std::size_t>(tape_pointers.size()));
    weights_.resize(static_cast<std::size_t>(tape_pointers.size()));
    mean_eta_ = Matrix::Zero(tape_pointers.size(), n_eta_);
    second_moment_.assign(
      static_cast<std::size_t>(tape_pointers.size()),
      Matrix::Zero(n_eta_, n_eta_));
    ode_support_tapes_.resize(static_cast<std::size_t>(tape_pointers.size()));
    owned_tapes_.resize(static_cast<std::size_t>(tape_pointers.size()));
    std::unordered_set<ObjectiveTape*> unique_tapes;
    for (int subject = 0; subject < tape_pointers.size(); ++subject) {
      Rcpp::XPtr<ObjectiveTape> tape(tape_pointers[subject]);
      if (tape->domain_names.size() != domain_) {
        throw std::invalid_argument(
          "A persistent weighted-ETA tape has an inconsistent domain.");
      }
      tapes_.push_back(tape.get());
      unique_tapes.insert(tape.get());
      const SEXP subject_input = subject_data[subject];
      subject_data_.push_back(subject_input);
      dynamic_values_.push_back(objective_dynamic_values(
        *tape, subject_input));
      points_[static_cast<std::size_t>(subject)].assign(domain_, 0.0);
    }
    native_threads_ = std::min({
      native_threads_, static_cast<int>(tapes_.size()),
      static_cast<int>(CPPAD_MAX_NUM_THREADS) - 1});
    if (unique_tapes.size() != tapes_.size()) {
      native_threads_ = 1;
      native_parallel_fallback_reason_ =
        "one or more subjects share a mutable CppAD objective tape";
    }
    if (native_threads_ > 1) {
      subject_pool_ = std::make_unique<NativeSubjectPool>(native_threads_);
    }
    if (use_ode_) reduced_requested_ = false;
    R_PreserveObject(retained_subject_data_);
  }

  ~WeightedEtaCollection() {
    subject_pool_.reset();
    if (retained_engine_ != R_NilValue) R_ReleaseObject(retained_engine_);
    if (retained_subject_data_ != R_NilValue) {
      R_ReleaseObject(retained_subject_data_);
    }
  }

  WeightedEtaCollection(const WeightedEtaCollection&) = delete;
  WeightedEtaCollection& operator=(const WeightedEtaCollection&) = delete;

  void set_grids(const Rcpp::List& eta_input,
                 const Rcpp::List& weight_input) {
    if (eta_input.size() != static_cast<int>(tapes_.size()) ||
        weight_input.size() != static_cast<int>(tapes_.size())) {
      throw std::invalid_argument(
        "Weighted ETA grids require one matrix and weight vector per subject.");
    }
    common_support_ = true;
    int reference_support = -1;
    std::vector<double> reference_weights;
    for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
      const Rcpp::NumericMatrix eta(eta_input[static_cast<int>(subject)]);
      const Rcpp::NumericVector weights(weight_input[static_cast<int>(subject)]);
      if (eta.nrow() < 1 || eta.ncol() != n_eta_ ||
          weights.size() != eta.nrow()) {
        throw std::invalid_argument("A weighted ETA grid has invalid dimensions.");
      }
      Matrix grid(eta.nrow(), eta.ncol());
      std::vector<double> normalized(static_cast<std::size_t>(weights.size()));
      double total = 0.0;
      for (int row = 0; row < eta.nrow(); ++row) {
        const double weight = weights[row];
        if (!std::isfinite(weight) || weight < 0.0) {
          throw std::invalid_argument("Weighted ETA probabilities must be finite and non-negative.");
        }
        normalized[static_cast<std::size_t>(row)] = weight;
        total += weight;
        for (int effect = 0; effect < n_eta_; ++effect) {
          const double value = eta(row, effect);
          if (!std::isfinite(value)) {
            throw std::invalid_argument("Weighted ETA support points must be finite.");
          }
          grid(row, effect) = value;
        }
      }
      if (!(total > 0.0) || !std::isfinite(total)) {
        throw std::invalid_argument("Weighted ETA probabilities have zero mass.");
      }
      for (double& weight : normalized) weight /= total;
      grids_[subject] = std::move(grid);
      weights_[subject] = std::move(normalized);
      if (reference_support < 0) {
        reference_support = eta.nrow();
        reference_weights = weights_[subject];
      } else if (reference_support != eta.nrow() ||
                 reference_weights != weights_[subject]) {
        common_support_ = false;
      }
    }
    if (common_support_) {
      common_weights_ = std::move(reference_weights);
      for (std::vector<double>& subject_weights : weights_) {
        subject_weights.clear();
        subject_weights.shrink_to_fit();
      }
    } else {
      common_weights_.clear();
    }
    ++grid_updates_;
    recompute_moments();
    reduced_dynamic_dirty_ = true;
    observe_support_shape();
  }

  Rcpp::List set_importance(
      const Rcpp::NumericVector& theta,
      const Rcpp::NumericVector& sigma,
      const Rcpp::NumericVector& omega,
      const Rcpp::List& eta_input,
      const Rcpp::List& log_proposal_input) {
    validate_parameters(theta, sigma, omega);
    if (eta_input.size() != static_cast<int>(tapes_.size()) ||
        log_proposal_input.size() != static_cast<int>(tapes_.size())) {
      throw std::invalid_argument(
        "Native IMP installation requires one proposal per subject.");
    }
    Rcpp::List uniform(static_cast<int>(tapes_.size()));
    for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
      const Rcpp::NumericMatrix eta(eta_input[static_cast<int>(subject)]);
      const Rcpp::NumericVector log_proposal(
        log_proposal_input[static_cast<int>(subject)]);
      if (eta.nrow() < 1 || eta.ncol() != n_eta_ ||
          log_proposal.size() != eta.nrow()) {
        throw std::invalid_argument(
          "A native IMP proposal has invalid dimensions.");
      }
      uniform[static_cast<int>(subject)] = Rcpp::NumericVector(
        eta.nrow(), 1.0 / static_cast<double>(eta.nrow()));
    }
    set_grids(eta_input, uniform);
    common_support_ = false;
    common_weights_.clear();
    weights_.assign(tapes_.size(), std::vector<double>());
    const std::vector<double> theta_native =
      Rcpp::as<std::vector<double>>(theta);
    const std::vector<double> sigma_native =
      Rcpp::as<std::vector<double>>(sigma);
    const std::vector<double> omega_native =
      Rcpp::as<std::vector<double>>(omega);
    const int retape_limit = ode_retape_limit();
    std::vector<std::vector<double>> objective_values(tapes_.size());
    bool completed = false;
    for (int attempt = 0; attempt < retape_limit && !completed; ++attempt) {
      prepare_native_points(theta_native, sigma_native, omega_native);
      std::vector<int> retape_requested(tapes_.size(), 0);
      std::vector<Eigen::Index> retape_supports(tapes_.size(), -1);
      std::vector<std::vector<double>> retape_points(tapes_.size());
      const auto evaluate_subject = [&](std::size_t subject) {
        std::vector<double>& point = points_[subject];
        const Eigen::Index count = subject_support_count(subject);
        std::vector<double>& values = objective_values[subject];
        values.assign(static_cast<std::size_t>(count),
                      std::numeric_limits<double>::infinity());
        std::ostringstream messages;
        Eigen::Index active = -1;
        try {
          for (Eigen::Index support = 0; support < count; ++support) {
            active = support;
            for (int effect = 0; effect < n_eta_; ++effect) {
              point[static_cast<std::size_t>(n_theta_ + effect)] =
                grids_[subject](support, effect);
            }
            ObjectiveTape& tape = objective_tape(subject, support);
            const std::vector<double> evaluated =
              tape.fun.Forward(0, point, messages);
            require_unchanged_path(tape.fun, "native IMP expectation");
            if (evaluated.size() != 1U || !std::isfinite(evaluated[0])) {
              throw std::domain_error(
                "A native IMP support objective was non-finite.");
            }
            values[static_cast<std::size_t>(support)] = evaluated[0];
          }
        } catch (const TapePathChange& change) {
          retape_requested[subject] = 1;
          retape_supports[subject] = active;
          retape_points[subject] = change.point().empty() ? point : change.point();
        }
      };
      if (subject_pool_) {
        subject_pool_->run_cppad(tapes_.size(), evaluate_subject);
      } else {
        for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
          evaluate_subject(subject);
          if ((subject + 1U) % 64U == 0U) Rcpp::checkUserInterrupt();
        }
      }
      bool retape = false;
      for (int requested : retape_requested) retape = retape || requested != 0;
      if (!retape) {
        completed = true;
        break;
      }
      if (!use_ode_) {
        const auto found = std::find(
          retape_requested.begin(), retape_requested.end(), 1);
        const std::size_t subject = static_cast<std::size_t>(
          std::distance(retape_requested.begin(), found));
        throw TapePathChange(
          "native IMP expectation", retape_points[subject]);
      }
      for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
        if (retape_requested[subject]) {
          retape_subject(
            subject, retape_supports[subject], retape_points[subject], true);
        }
      }
    }
    if (!completed) {
      throw std::runtime_error(
        "Native IMP expectation exceeded the ODE retaping limit.");
    }
    Rcpp::NumericVector ess(static_cast<R_xlen_t>(tapes_.size()));
    for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
      const Rcpp::NumericVector log_proposal(
        log_proposal_input[static_cast<int>(subject)]);
      const std::vector<double>& values = objective_values[subject];
      double maximum = -std::numeric_limits<double>::infinity();
      std::vector<double> log_weight(values.size());
      for (std::size_t support = 0; support < values.size(); ++support) {
        log_weight[support] = -0.5 * values[support] -
          log_proposal[static_cast<R_xlen_t>(support)];
        maximum = std::max(maximum, log_weight[support]);
      }
      double total = 0.0;
      std::vector<double>& probability = weights_[subject];
      probability.resize(values.size());
      for (std::size_t support = 0; support < values.size(); ++support) {
        probability[support] = std::exp(log_weight[support] - maximum);
        total += probability[support];
      }
      if (!(total > 0.0) || !std::isfinite(total)) {
        throw std::domain_error("Native IMP importance weights have zero mass.");
      }
      double squares = 0.0;
      for (double& value : probability) {
        value /= total;
        squares += value * value;
      }
      ess[static_cast<R_xlen_t>(subject)] = 1.0 / squares;
    }
    recompute_moments();
    reduced_dynamic_dirty_ = true;
    ++importance_updates_;
    point_evaluations_ += support_point_count();
    last_point_evaluations_ = support_point_count();
    return Rcpp::List::create(
      Rcpp::Named("ess") = ess,
      Rcpp::Named("subjects") = static_cast<int>(tapes_.size()),
      Rcpp::Named("support_points") =
        static_cast<double>(support_point_count()),
      Rcpp::Named("backend") =
        "persistent-native-importance-installation");
  }

  void update_common(const Rcpp::NumericMatrix& eta, double gamma,
                     int max_support = 0, double prune_tolerance = 0.0) {
    if (eta.nrow() != static_cast<int>(tapes_.size()) ||
        eta.ncol() != n_eta_ || !std::isfinite(gamma) || gamma <= 0.0 ||
        gamma > 1.0 || max_support < 0 || !std::isfinite(prune_tolerance) ||
        prune_tolerance < 0.0 || prune_tolerance >= 1.0) {
      throw std::invalid_argument(
        "The stochastic weighted-ETA update is invalid.");
    }
    if (support_count() && !common_support_) {
      throw std::logic_error(
        "A common stochastic update cannot extend subject-specific weights.");
    }
    const bool reset = support_count() == 0 ||
      gamma >= 1.0 - std::numeric_limits<double>::epsilon();
    for (int subject = 0; subject < eta.nrow(); ++subject) {
      for (int effect = 0; effect < n_eta_; ++effect) {
        if (!std::isfinite(eta(subject, effect))) {
          throw std::invalid_argument("Stochastic ETA states must be finite.");
        }
      }
    }
    if (reset) {
      common_weights_.assign(1U, 1.0);
    } else {
      for (double& weight : common_weights_) weight *= 1.0 - gamma;
      common_weights_.push_back(gamma);
      normalize(common_weights_);
    }
    const Eigen::Index previous_support = reset ? 0 :
      static_cast<Eigen::Index>(common_weights_.size() - 1U);
    for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
      if (reset) {
        if (grids_[subject].rows() < 1) grids_[subject].resize(1, n_eta_);
        for (int effect = 0; effect < n_eta_; ++effect) {
          grids_[subject](0, effect) = eta(static_cast<int>(subject), effect);
        }
      } else {
        if (grids_[subject].rows() <= previous_support) {
          const Eigen::Index capacity = std::max<Eigen::Index>(
            previous_support + 1,
            std::max<Eigen::Index>(4, 2 * grids_[subject].rows()));
          Matrix next(capacity, n_eta_);
          if (previous_support) {
            next.topRows(previous_support) =
              grids_[subject].topRows(previous_support);
          }
          grids_[subject].swap(next);
          ++support_reallocations_;
        }
        for (int effect = 0; effect < n_eta_; ++effect) {
          grids_[subject](previous_support, effect) =
            eta(static_cast<int>(subject), effect);
        }
      }
      weights_[subject].clear();
    }
    common_support_ = true;
    if (!reset && (prune_tolerance > 0.0 || max_support > 0)) {
      compress_common(max_support, prune_tolerance);
    } else {
      if (reset) {
        for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
          const Vector current = grids_[subject].row(0).transpose();
          mean_eta_.row(static_cast<Eigen::Index>(subject)) = current.transpose();
          second_moment_[subject].noalias() = current * current.transpose();
        }
      } else {
        for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
          Vector current(n_eta_);
          for (int effect = 0; effect < n_eta_; ++effect) {
            current[effect] = eta(static_cast<int>(subject), effect);
          }
          mean_eta_.row(static_cast<Eigen::Index>(subject)) =
            (1.0 - gamma) * mean_eta_.row(static_cast<Eigen::Index>(subject)) +
            gamma * current.transpose();
          second_moment_[subject] =
            (1.0 - gamma) * second_moment_[subject] +
            gamma * current * current.transpose();
        }
      }
    }
    reduced_dynamic_dirty_ = true;
    observe_support_shape();
    ++stochastic_updates_;
  }

  Rcpp::List evaluate_aggregate(
      const Rcpp::NumericVector& theta, const Rcpp::NumericVector& sigma,
      const Rcpp::NumericVector& omega) {
    const SaemEvaluation native = evaluate_native(
      Rcpp::as<std::vector<double>>(theta),
      Rcpp::as<std::vector<double>>(sigma),
      Rcpp::as<std::vector<double>>(omega));
    Rcpp::NumericVector gradient = libertad::eigen_vector_to_r(native.gradient);
    if (!tapes_.empty()) gradient.attr("names") =
      Rcpp::wrap(tapes_.front()->domain_names);
    return Rcpp::List::create(
      Rcpp::Named("value") = native.value,
      Rcpp::Named("gradient") = gradient,
      Rcpp::Named("evaluations") = static_cast<double>(evaluations_),
      Rcpp::Named("point_evaluations") =
        static_cast<double>(last_point_evaluations_));
  }

  SaemEvaluation evaluate_native(
      const std::vector<double>& theta, const std::vector<double>& sigma,
      const std::vector<double>& omega) {
    validate_parameters_native(theta, sigma, omega);
    if (!support_count()) {
      throw std::logic_error("The weighted-ETA context has no support points.");
    }
    if (ensure_reduced_tape(theta, sigma, omega)) {
      return evaluate_reduced_native(theta, sigma, omega, true);
    }
    const std::vector<double> reverse_weight(1U, 1.0);
    const int retape_limit = ode_retape_limit();
    for (int attempt = 0; attempt < retape_limit; ++attempt) {
      prepare_native_points(theta, sigma, omega);
      std::vector<double> subject_values(tapes_.size(), 0.0);
      std::vector<std::vector<double>> subject_gradients(
        tapes_.size(), std::vector<double>(domain_, 0.0));
      std::vector<long long> subject_evaluations(tapes_.size(), 0);
      std::vector<int> retape_requested(tapes_.size(), 0);
      std::vector<Eigen::Index> active_support(tapes_.size(), -1);
      std::vector<Eigen::Index> retape_supports(tapes_.size(), -1);
      std::vector<std::vector<double>> retape_points(tapes_.size());
      const auto evaluate_subject = [&](std::size_t subject) {
        std::vector<double>& point = points_[subject];
        double subject_value = 0.0;
        std::vector<double>& subject_gradient = subject_gradients[subject];
        std::ostringstream messages;
        try {
          for (Eigen::Index support = 0;
               support < subject_support_count(subject); ++support) {
            active_support[subject] = support;
            ObjectiveTape& tape = objective_tape(subject, support);
            for (int effect = 0; effect < n_eta_; ++effect) {
              point[static_cast<std::size_t>(n_theta_ + effect)] =
                grids_[subject](support, effect);
            }
            const std::vector<double> value = tape.fun.Forward(0, point, messages);
            require_unchanged_path(tape.fun, "persistent weighted-ETA objective");
            if (value.size() != 1U || !std::isfinite(value[0])) {
              throw std::domain_error("A weighted-ETA objective was non-finite.");
            }
            const double probability = probability_at(subject, support);
            subject_value += probability * value[0];
            const std::vector<double> derivative =
              tape.fun.Reverse(1, reverse_weight);
            require_unchanged_path(tape.fun, "persistent weighted-ETA gradient");
            if (derivative.size() != domain_) {
              throw std::runtime_error(
                "A weighted-ETA tape returned an invalid gradient length.");
            }
            for (std::size_t column = 0; column < domain_; ++column) {
              subject_gradient[column] += probability * derivative[column];
            }
            ++subject_evaluations[subject];
          }
          subject_values[subject] = subject_value;
        } catch (const TapePathChange& change) {
          retape_requested[subject] = 1;
          retape_points[subject] = change.point().empty() ? point : change.point();
          retape_supports[subject] = active_support[subject];
        }
      };
      if (subject_pool_) {
        subject_pool_->run_cppad(tapes_.size(), evaluate_subject);
      } else {
        for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
          evaluate_subject(subject);
          if ((subject + 1U) % 64U == 0U) Rcpp::checkUserInterrupt();
        }
      }
      bool retape = false;
      for (int requested : retape_requested) retape = retape || requested != 0;
      if (retape) {
        if (!use_ode_) {
          const auto found = std::find(retape_requested.begin(),
                                       retape_requested.end(), 1);
          const std::size_t subject = static_cast<std::size_t>(
            std::distance(retape_requested.begin(), found));
          throw TapePathChange(
            "persistent weighted-ETA objective", retape_points[subject]);
        }
        for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
          if (retape_requested[subject]) {
            retape_subject(
              subject, retape_supports[subject], retape_points[subject], true);
          }
        }
        continue;
      }
      double value_total = 0.0;
      std::vector<double> gradient_total(domain_, 0.0);
      long long point_evaluations = 0;
      // Reduction order is deliberately independent of worker scheduling.
      for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
        value_total += subject_values[subject];
        for (std::size_t column = 0; column < domain_; ++column) {
          gradient_total[column] += subject_gradients[subject][column];
        }
        point_evaluations += subject_evaluations[subject];
      }
      ++evaluations_;
      point_evaluations_ += point_evaluations;
      last_point_evaluations_ = point_evaluations;
      SaemEvaluation result;
      result.value = value_total;
      result.gradient = Eigen::Map<Vector>(
        gradient_total.data(), static_cast<Eigen::Index>(gradient_total.size()));
      return result;
    }
    throw std::runtime_error(
      "Weighted ETA evaluation exceeded the ODE retaping limit.");
  }

  double value_native(
      const std::vector<double>& theta, const std::vector<double>& sigma,
      const std::vector<double>& omega) {
    validate_parameters_native(theta, sigma, omega);
    if (!support_count()) {
      throw std::logic_error("The weighted-ETA context has no support points.");
    }
    if (ensure_reduced_tape(theta, sigma, omega)) {
      return evaluate_reduced_native(theta, sigma, omega, false).value;
    }
    const int retape_limit = ode_retape_limit();
    for (int attempt = 0; attempt < retape_limit; ++attempt) {
      prepare_native_points(theta, sigma, omega);
      std::vector<double> subject_values(tapes_.size(), 0.0);
      std::vector<long long> subject_evaluations(tapes_.size(), 0);
      std::vector<int> retape_requested(tapes_.size(), 0);
      std::vector<Eigen::Index> active_support(tapes_.size(), -1);
      std::vector<Eigen::Index> retape_supports(tapes_.size(), -1);
      std::vector<std::vector<double>> retape_points(tapes_.size());
      const auto evaluate_subject = [&](std::size_t subject) {
        std::vector<double>& point = points_[subject];
        std::ostringstream messages;
        double subject_value = 0.0;
        try {
          for (Eigen::Index support = 0;
               support < subject_support_count(subject); ++support) {
            active_support[subject] = support;
            ObjectiveTape& tape = objective_tape(subject, support);
            for (int effect = 0; effect < n_eta_; ++effect) {
              point[static_cast<std::size_t>(n_theta_ + effect)] =
                grids_[subject](support, effect);
            }
            const std::vector<double> value = tape.fun.Forward(0, point, messages);
            require_unchanged_path(tape.fun, "persistent weighted-ETA value");
            if (value.size() != 1U || !std::isfinite(value[0])) {
              subject_value = 1e100;
              break;
            }
            subject_value += probability_at(subject, support) * value[0];
            ++subject_evaluations[subject];
          }
          subject_values[subject] = subject_value;
        } catch (const TapePathChange& change) {
          retape_requested[subject] = 1;
          retape_points[subject] = change.point().empty() ? point : change.point();
          retape_supports[subject] = active_support[subject];
        }
      };
      if (subject_pool_) {
        subject_pool_->run_cppad(tapes_.size(), evaluate_subject);
      } else {
        for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
          evaluate_subject(subject);
          if ((subject + 1U) % 64U == 0U) Rcpp::checkUserInterrupt();
        }
      }
      bool retape = false;
      for (int requested : retape_requested) retape = retape || requested != 0;
      if (retape) {
        if (!use_ode_) {
          const auto found = std::find(retape_requested.begin(),
                                       retape_requested.end(), 1);
          const std::size_t subject = static_cast<std::size_t>(
            std::distance(retape_requested.begin(), found));
          throw TapePathChange(
            "persistent weighted-ETA value", retape_points[subject]);
        }
        for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
          if (retape_requested[subject]) {
            retape_subject(
              subject, retape_supports[subject], retape_points[subject], true);
          }
        }
        continue;
      }
      double value_total = 0.0;
      long long point_evaluations = 0;
      for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
        value_total += subject_values[subject];
        point_evaluations += subject_evaluations[subject];
      }
      ++evaluations_;
      point_evaluations_ += point_evaluations;
      last_point_evaluations_ = point_evaluations;
      return std::isfinite(value_total) ? value_total : 1e100;
    }
    throw std::runtime_error(
      "Weighted ETA value evaluation exceeded the ODE retaping limit.");
  }

  Rcpp::NumericMatrix mean_eta() const {
    if (!support_count()) {
      throw std::logic_error("The weighted-ETA context has no support points.");
    }
    return libertad::eigen_matrix_to_r(mean_eta_);
  }

  int eta_dimension() const { return n_eta_; }

  Rcpp::NumericVector common_weights() const {
    if (!common_support_ || !support_count()) return Rcpp::NumericVector();
    return Rcpp::wrap(common_weights_);
  }

  void recenter(const Rcpp::NumericMatrix& adjustment) {
    if (adjustment.nrow() != static_cast<int>(tapes_.size()) ||
        adjustment.ncol() != n_eta_) {
      throw std::invalid_argument("ETA recentering adjustment has invalid dimensions.");
    }
    for (int subject = 0; subject < adjustment.nrow(); ++subject) {
      for (int effect = 0; effect < adjustment.ncol(); ++effect) {
        if (!std::isfinite(adjustment(subject, effect))) {
          throw std::invalid_argument("ETA recentering adjustment must be finite.");
        }
      }
    }
    for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
      Vector shift(n_eta_);
      for (int effect = 0; effect < n_eta_; ++effect) {
        shift[effect] = adjustment(static_cast<int>(subject), effect);
      }
      const Vector previous_mean =
        mean_eta_.row(static_cast<Eigen::Index>(subject)).transpose();
      second_moment_[subject].noalias() +=
        previous_mean * shift.transpose() +
        shift * previous_mean.transpose() + shift * shift.transpose();
      mean_eta_.row(static_cast<Eigen::Index>(subject)) += shift.transpose();
      for (int effect = 0; effect < n_eta_; ++effect) {
        const double value = adjustment(static_cast<int>(subject), effect);
        grids_[subject].topRows(subject_support_count(subject)).col(effect).array() +=
          value;
      }
    }
    reduced_dynamic_dirty_ = true;
    ++recenters_;
  }

  Rcpp::NumericVector omega_sufficient(
      int n_eta_base, int iov, const Rcpp::IntegerVector& omega_rows,
      const Rcpp::IntegerVector& omega_cols) const {
    if (!support_count() || n_eta_base < 1 || iov < 0 || iov > n_eta_base ||
        omega_rows.size() != omega_cols.size()) {
      throw std::invalid_argument("Weighted OMEGA sufficient-statistic inputs are invalid.");
    }
    const int between = n_eta_base - iov;
    const int occasions = iov ? (n_eta_ - between) / iov : 0;
    if ((!iov && n_eta_ != n_eta_base) ||
        (iov && (n_eta_ < between || (n_eta_ - between) % iov != 0 ||
                 occasions < 1))) {
      throw std::invalid_argument("Weighted ETA columns do not match the IOV layout.");
    }
    Matrix covariance = Matrix::Zero(n_eta_base, n_eta_base);
    if (!iov) {
      for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
        covariance.noalias() += second_moment_[subject];
      }
      covariance /= static_cast<double>(tapes_.size());
    } else {
      if (between > 0) {
        for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
          covariance.topLeftCorner(between, between).noalias() +=
            second_moment_[subject].topLeftCorner(between, between);
        }
        covariance.topLeftCorner(between, between) /=
          static_cast<double>(tapes_.size());
      }
      for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
        for (int occasion = 0; occasion < occasions; ++occasion) {
          const int offset = between + occasion * iov;
          covariance.bottomRightCorner(iov, iov).noalias() +=
            second_moment_[subject].block(offset, offset, iov, iov);
        }
      }
      covariance.bottomRightCorner(iov, iov) /=
        static_cast<double>(tapes_.size() * static_cast<std::size_t>(occasions));
    }
    covariance.diagonal().array() += 1e-8;
    Rcpp::NumericVector result(omega_rows.size());
    for (R_xlen_t entry = 0; entry < omega_rows.size(); ++entry) {
      const int row = omega_rows[entry] - 1;
      const int column = omega_cols[entry] - 1;
      if (row < 0 || column < 0 || row >= n_eta_base || column >= n_eta_base) {
        throw std::invalid_argument("An OMEGA sufficient-statistic coordinate is invalid.");
      }
      result[entry] = covariance(row, column);
    }
    return result;
  }

  Rcpp::NumericVector sigma_expectation(
      SEXP engine_pointer, const Rcpp::DataFrame& data,
      const Rcpp::NumericVector& theta,
      const Rcpp::NumericVector& sigma) const {
    if (!common_support_ || !support_count()) {
      throw std::logic_error(
        "Native SIGMA expectation requires a common SAEM support.");
    }
    Rcpp::XPtr<ModelEngine> engine(engine_pointer);
    require_materialized_addl(data);
    if (engine->error_type != "additive" &&
        engine->error_type != "proportional" &&
        engine->error_type != "exponential") {
      throw std::invalid_argument(
        "Native weighted SIGMA updates require a simple residual model.");
    }
    const Rcpp::IntegerVector evid = data["EVID"];
    const Rcpp::IntegerVector mdv = data["MDV"];
    const Rcpp::NumericVector dv = data["DV"];
    Rcpp::IntegerVector dvid(data.nrows(), 1);
    if (data.containsElementNamed("DVID")) dvid = data["DVID"];
    std::vector<double> expected_variance(
      static_cast<std::size_t>(sigma.size()), 0.0);
    const std::vector<double>& support_weights = common_weights_;
    for (std::size_t support = 0; support < support_weights.size(); ++support) {
      Rcpp::NumericMatrix eta(static_cast<int>(tapes_.size()), n_eta_);
      for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
        for (int effect = 0; effect < n_eta_; ++effect) {
          eta(static_cast<int>(subject), effect) =
            grids_[subject](static_cast<Eigen::Index>(support), effect);
        }
      }
      const Rcpp::List simulation = simulate(*engine, data, theta, eta, sigma);
      const Rcpp::NumericVector prediction = simulation["ipred"];
      std::vector<double> sum_square(static_cast<std::size_t>(sigma.size()), 0.0);
      std::vector<int> count(static_cast<std::size_t>(sigma.size()), 0);
      for (int row = 0; row < data.nrows(); ++row) {
        if (evid[row] != 0 || mdv[row] != 0 || !std::isfinite(dv[row]) ||
            !std::isfinite(prediction[row])) continue;
        const int response = std::max(dvid[row], 1) - 1;
        if (response < 0 || response >= sigma.size()) continue;
        double residual = 0.0;
        if (engine->error_type == "additive") {
          residual = dv[row] - prediction[row];
        } else if (engine->error_type == "proportional") {
          residual = (dv[row] - prediction[row]) /
            std::max(std::abs(prediction[row]), 1e-12);
        } else {
          if (!(dv[row] > 0.0) || !(prediction[row] > 0.0)) continue;
          residual = std::log(dv[row]) - std::log(prediction[row]);
        }
        if (!std::isfinite(residual)) continue;
        sum_square[static_cast<std::size_t>(response)] += residual * residual;
        ++count[static_cast<std::size_t>(response)];
      }
      for (R_xlen_t response = 0; response < sigma.size(); ++response) {
        double variance = engine->sigma_parameterization == "variance" ?
          sigma[response] : sigma[response] * sigma[response];
        if (count[static_cast<std::size_t>(response)] > 0) {
          variance = sum_square[static_cast<std::size_t>(response)] /
            count[static_cast<std::size_t>(response)];
        }
        expected_variance[static_cast<std::size_t>(response)] +=
          support_weights[support] * variance;
      }
      if ((support + 1U) % 16U == 0U) Rcpp::checkUserInterrupt();
    }
    Rcpp::NumericVector result(sigma.size());
    for (R_xlen_t response = 0; response < sigma.size(); ++response) {
      result[response] = engine->sigma_parameterization == "variance" ?
        expected_variance[static_cast<std::size_t>(response)] :
        std::sqrt(std::max(0.0, expected_variance[static_cast<std::size_t>(response)]));
    }
    return result;
  }

  Rcpp::List telemetry() const {
    return Rcpp::List::create(
      Rcpp::Named("backend") = "persistent-cpp-weighted-eta",
      Rcpp::Named("support") = support_count(),
      Rcpp::Named("common_support") = common_support_,
      Rcpp::Named("weight_vectors") = common_support_ ? 1 :
        static_cast<int>(weights_.size()),
      Rcpp::Named("grid_updates") = static_cast<double>(grid_updates_),
      Rcpp::Named("stochastic_updates") = static_cast<double>(stochastic_updates_),
      Rcpp::Named("importance_updates") = static_cast<double>(importance_updates_),
      Rcpp::Named("evaluations") = static_cast<double>(evaluations_),
      Rcpp::Named("point_evaluations") = static_cast<double>(point_evaluations_),
      Rcpp::Named("dynamic_updates") = static_cast<double>(dynamic_updates_),
      Rcpp::Named("dynamic_cache_hits") =
        static_cast<double>(dynamic_cache_hits_),
      Rcpp::Named("native_subject_threads") = native_threads_,
      Rcpp::Named("native_subject_parallel") =
        static_cast<bool>(subject_pool_),
      Rcpp::Named("native_subject_dispatches") = subject_pool_ ?
        static_cast<double>(subject_pool_->dispatches()) : 0.0,
      Rcpp::Named("cppad_subject_dispatches") = subject_pool_ ?
        static_cast<double>(subject_pool_->cppad_dispatches()) : 0.0,
      Rcpp::Named("tape_records") = static_cast<double>(tape_records_),
      Rcpp::Named("tape_retapes") = static_cast<double>(tape_retapes_),
      Rcpp::Named("ode_support_tape_limit") = ode_support_tape_limit_,
      Rcpp::Named("native_parallel_fallback_reason") =
        native_parallel_fallback_reason_,
      Rcpp::Named("reduced_population_requested") = reduced_requested_,
      Rcpp::Named("reduced_population_available") =
        static_cast<bool>(reduced_tape_),
      Rcpp::Named("reduced_population_records") =
        static_cast<double>(reduced_records_),
      Rcpp::Named("reduced_population_evaluations") =
        static_cast<double>(reduced_evaluations_),
      Rcpp::Named("reduced_population_dynamic_updates") =
        static_cast<double>(reduced_dynamic_updates_),
      Rcpp::Named("reduced_population_shape_stability") =
        reduced_shape_stability_,
      Rcpp::Named("reduced_population_fallback_reason") =
        reduced_fallback_reason_,
      Rcpp::Named("recenters") = static_cast<double>(recenters_),
      Rcpp::Named("compressed_points") = static_cast<double>(compressed_points_),
      Rcpp::Named("support_reallocations") =
        static_cast<double>(support_reallocations_),
      Rcpp::Named("pruned_mass") = pruned_mass_);
  }

 private:
  int n_theta_ = 0;
  int n_eta_ = 0;
  int n_sigma_ = 0;
  int n_omega_ = 0;
  std::size_t domain_ = 0U;
  ModelEngine* engine_ = nullptr;
  SEXP retained_engine_ = R_NilValue;
  SEXP retained_subject_data_ = R_NilValue;
  bool use_ode_ = false;
  std::vector<ObjectiveTape*> tapes_;
  std::vector<std::unique_ptr<ObjectiveTape>> owned_tapes_;
  std::vector<std::vector<std::unique_ptr<ObjectiveTape>>> ode_support_tapes_;
  std::vector<SEXP> subject_data_;
  std::vector<std::vector<double>> dynamic_values_;
  std::vector<std::vector<double>> points_;
  std::vector<Matrix> grids_;
  std::vector<std::vector<double>> weights_;
  std::vector<double> common_weights_;
  Matrix mean_eta_;
  std::vector<Matrix> second_moment_;
  bool common_support_ = false;
  long long grid_updates_ = 0;
  long long stochastic_updates_ = 0;
  long long importance_updates_ = 0;
  long long evaluations_ = 0;
  long long point_evaluations_ = 0;
  long long last_point_evaluations_ = 0;
  long long dynamic_updates_ = 0;
  long long dynamic_cache_hits_ = 0;
  long long tape_records_ = 0;
  long long tape_retapes_ = 0;
  int native_threads_ = 1;
  std::unique_ptr<NativeSubjectPool> subject_pool_;
  std::string native_parallel_fallback_reason_;
  long long recenters_ = 0;
  long long compressed_points_ = 0;
  long long support_reallocations_ = 0;
  double pruned_mass_ = 0.0;
  bool reduced_requested_ = false;
  bool reduced_dynamic_dirty_ = true;
  std::unique_ptr<ObjectiveTape> reduced_tape_;
  std::vector<int> reduced_shape_;
  std::vector<int> observed_shape_;
  int reduced_shape_stability_ = 0;
  long long reduced_records_ = 0;
  long long reduced_evaluations_ = 0;
  long long reduced_dynamic_updates_ = 0;
  std::string reduced_fallback_reason_;
  double reduced_max_operations_ = 2e6;
  int ode_support_tape_limit_ = 4096;

  ObjectiveTape& objective_tape(std::size_t subject,
                                Eigen::Index support) {
    if (!use_ode_) return *tapes_[subject];
    if (support < 0 ||
        support >= static_cast<Eigen::Index>(ode_support_tapes_[subject].size()) ||
        !ode_support_tapes_[subject][static_cast<std::size_t>(support)]) {
      throw std::logic_error("An ODE weighted-support tape is unavailable.");
    }
    return *ode_support_tapes_[subject][static_cast<std::size_t>(support)];
  }

  int ode_retape_limit() const {
    if (!use_ode_) return 1;
    Eigen::Index maximum = 0;
    for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
      maximum = std::max(maximum, subject_support_count(subject));
    }
    return static_cast<int>(std::max<Eigen::Index>(8, maximum + 2));
  }

  void retape_subject(std::size_t subject, Eigen::Index support,
                      const std::vector<double>& point, bool retape) {
    if (!use_ode_ || !engine_ || subject >= subject_data_.size() ||
        support < 0 || point.size() != domain_) {
      throw std::runtime_error(
        "A weighted ETA objective cannot be retaped for this model.");
    }
    Rcpp::NumericVector theta(n_theta_), sigma(n_sigma_), omega(n_omega_);
    Rcpp::NumericMatrix eta(1, n_eta_);
    for (int index = 0; index < n_theta_; ++index) theta[index] = point[index];
    for (int index = 0; index < n_eta_; ++index) {
      eta(0, index) = point[static_cast<std::size_t>(n_theta_ + index)];
    }
    const int sigma_offset = n_theta_ + n_eta_;
    for (int index = 0; index < n_sigma_; ++index) {
      sigma[index] = point[static_cast<std::size_t>(sigma_offset + index)];
    }
    const int omega_offset = sigma_offset + n_sigma_;
    for (int index = 0; index < n_omega_; ++index) {
      omega[index] = point[static_cast<std::size_t>(omega_offset + index)];
    }
    const EventDataView data = event_data_view(subject_data_[subject]);
    std::unique_ptr<ObjectiveTape> recorded = record_objective_tape(
      *engine_, data, theta, eta, sigma, omega, true);
    auto& support_tapes = ode_support_tapes_[subject];
    const std::size_t index = static_cast<std::size_t>(support);
    if (support_tapes.size() <= index) support_tapes.resize(index + 1U);
    support_tapes[index] = std::move(recorded);
    points_[subject] = point;
    ++tape_records_;
    if (retape) ++tape_retapes_;
  }

  void prepare_native_points(
      const std::vector<double>& theta, const std::vector<double>& sigma,
      const std::vector<double>& omega) {
    // Dynamic input mutation and all R-facing checks stay on the main thread.
    // Workers subsequently own one independent tape and point buffer each.
    if (use_ode_ && support_point_count() > ode_support_tape_limit_) {
      throw std::runtime_error(
        "The ODE weighted ETA support exceeds the native retained-tape limit; "
        "increase LibeRation.ode_weighted_tape_limit or use the fallback path.");
    }
    for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
      std::vector<double>& point = points_[subject];
      std::copy(theta.begin(), theta.end(), point.begin());
      std::copy(sigma.begin(), sigma.end(), point.begin() + n_theta_ + n_eta_);
      std::copy(omega.begin(), omega.end(),
                point.begin() + n_theta_ + n_eta_ + n_sigma_);
      if (!use_ode_) {
        ObjectiveTape& tape = *tapes_[subject];
        if (tape.dynamic_values != dynamic_values_[subject]) {
          set_tape_dynamic_values(
            tape, dynamic_values_[subject], "Weighted ETA objective tape");
          ++dynamic_updates_;
        } else {
          ++dynamic_cache_hits_;
        }
        continue;
      }
      auto& support_tapes = ode_support_tapes_[subject];
      const Eigen::Index count = subject_support_count(subject);
      if (support_tapes.size() < static_cast<std::size_t>(count)) {
        support_tapes.resize(static_cast<std::size_t>(count));
      }
      for (Eigen::Index support = 0; support < count; ++support) {
        for (int effect = 0; effect < n_eta_; ++effect) {
          point[static_cast<std::size_t>(n_theta_ + effect)] =
            grids_[subject](support, effect);
        }
        if (!support_tapes[static_cast<std::size_t>(support)]) {
          retape_subject(subject, support, point, false);
        }
        ObjectiveTape& tape =
          *support_tapes[static_cast<std::size_t>(support)];
        if (tape.dynamic_values != dynamic_values_[subject]) {
          set_tape_dynamic_values(
            tape, dynamic_values_[subject], "Weighted ODE support tape");
          ++dynamic_updates_;
        } else {
          ++dynamic_cache_hits_;
        }
      }
    }
  }

  int support_count() const {
    if (grids_.empty()) return 0;
    return common_support_ ? static_cast<int>(common_weights_.size()) :
      static_cast<int>(grids_.front().rows());
  }

  Eigen::Index subject_support_count(std::size_t subject) const {
    return common_support_ ? static_cast<Eigen::Index>(common_weights_.size()) :
      grids_[subject].rows();
  }

  double probability_at(std::size_t subject, Eigen::Index support) const {
    return common_support_ ?
      common_weights_[static_cast<std::size_t>(support)] :
      weights_[subject][static_cast<std::size_t>(support)];
  }

  std::vector<int> support_shape() const {
    std::vector<int> result;
    result.reserve(grids_.size());
    for (std::size_t subject = 0; subject < grids_.size(); ++subject) {
      result.push_back(static_cast<int>(subject_support_count(subject)));
    }
    return result;
  }

  long long support_point_count() const {
    long long result = 0;
    for (std::size_t subject = 0; subject < grids_.size(); ++subject) {
      result += subject_support_count(subject);
    }
    return result;
  }

  void observe_support_shape() {
    const std::vector<int> shape = support_shape();
    if (shape == observed_shape_) {
      if (reduced_shape_stability_ < std::numeric_limits<int>::max()) {
        ++reduced_shape_stability_;
      }
    } else {
      observed_shape_ = shape;
      reduced_shape_stability_ = 1;
      if (shape != reduced_shape_) {
        reduced_tape_.reset();
        reduced_shape_.clear();
      }
    }
  }

  std::vector<double> reduced_dynamic_values() const {
    std::vector<double> result;
    std::size_t reserve = 0U;
    for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
      reserve += dynamic_values_[subject].size() +
        static_cast<std::size_t>(subject_support_count(subject)) *
          static_cast<std::size_t>(n_eta_ + 1);
    }
    result.reserve(reserve);
    for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
      result.insert(result.end(), dynamic_values_[subject].begin(),
                    dynamic_values_[subject].end());
      for (Eigen::Index support = 0;
           support < subject_support_count(subject); ++support) {
        for (int effect = 0; effect < n_eta_; ++effect) {
          result.push_back(grids_[subject](support, effect));
        }
        result.push_back(probability_at(subject, support));
      }
    }
    return result;
  }

  std::vector<double> reduced_point(
      const std::vector<double>& theta, const std::vector<double>& sigma,
      const std::vector<double>& omega) const {
    std::vector<double> result;
    result.reserve(theta.size() + sigma.size() + omega.size());
    result.insert(result.end(), theta.begin(), theta.end());
    result.insert(result.end(), sigma.begin(), sigma.end());
    result.insert(result.end(), omega.begin(), omega.end());
    return result;
  }

  bool ensure_reduced_tape(
      const std::vector<double>& theta, const std::vector<double>& sigma,
      const std::vector<double>& omega) {
    if (!reduced_requested_) return false;
    const std::vector<int> shape = support_shape();
    // A population aggregate pays off only when its domain survives across
    // expectation updates. Progressive IMP and uncompressed SAEM change the
    // support dimensions almost every iteration; retaping those large graphs
    // is slower than evaluating the persistent subject tapes directly.
    if (shape != observed_shape_ || reduced_shape_stability_ < 2) {
      reduced_fallback_reason_ =
        "support shape has not remained stable across expectation updates";
      return false;
    }
    double estimated_operations = 0.0;
    for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
      estimated_operations += static_cast<double>(tapes_[subject]->fun.size_op()) *
        static_cast<double>(subject_support_count(subject));
    }
    if (!std::isfinite(estimated_operations) ||
        estimated_operations > reduced_max_operations_) {
      reduced_fallback_reason_ =
        "estimated reduced population tape exceeds the operation limit";
      reduced_tape_.reset();
      reduced_shape_.clear();
      return false;
    }
    if (!reduced_tape_ || shape != reduced_shape_) {
      record_reduced_tape(theta, sigma, omega, shape);
    } else if (reduced_dynamic_dirty_) {
      const std::vector<double> dynamic = reduced_dynamic_values();
      set_tape_dynamic_values(
        *reduced_tape_, dynamic, "Reduced weighted population tape");
      ++reduced_dynamic_updates_;
      reduced_dynamic_dirty_ = false;
    }
    return static_cast<bool>(reduced_tape_);
  }

  void record_reduced_tape(
      const std::vector<double>& theta, const std::vector<double>& sigma,
      const std::vector<double>& omega, const std::vector<int>& shape) {
    using AD = CppAD::AD<double>;
    const std::vector<double> point = reduced_point(theta, sigma, omega);
    const std::vector<double> dynamic_values = reduced_dynamic_values();
    std::vector<AD> independent(point.begin(), point.end());
    std::vector<AD> dynamic(dynamic_values.begin(), dynamic_values.end());
    if (dynamic.empty()) CppAD::Independent(independent);
    else CppAD::Independent(independent, dynamic);
    CppADRecordingGuard<double> recording;
    std::size_t cursor = 0U;
    AD total = AD(0.0);
    std::ostringstream messages;
    for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
      ObjectiveTape& source = *tapes_[subject];
      auto nested = source.fun.base2ad();
      const std::size_t dynamic_count = source.fun.size_dyn_ind();
      if (cursor + dynamic_count > dynamic.size()) {
        throw std::logic_error(
          "Reduced weighted population dynamic offsets are inconsistent.");
      }
      if (dynamic_count) {
        std::vector<AD> subject_dynamic(
          dynamic.begin() + static_cast<std::ptrdiff_t>(cursor),
          dynamic.begin() + static_cast<std::ptrdiff_t>(cursor + dynamic_count));
        nested.new_dynamic(subject_dynamic);
      }
      cursor += dynamic_count;
      for (Eigen::Index support = 0;
           support < subject_support_count(subject); ++support) {
        if (cursor + static_cast<std::size_t>(n_eta_ + 1) > dynamic.size()) {
          throw std::logic_error(
            "Reduced weighted ETA offsets are inconsistent.");
        }
        std::vector<AD> source_point;
        source_point.reserve(domain_);
        source_point.insert(source_point.end(), independent.begin(),
                            independent.begin() + n_theta_);
        for (int effect = 0; effect < n_eta_; ++effect) {
          source_point.push_back(dynamic[cursor++]);
        }
        source_point.insert(
          source_point.end(), independent.begin() + n_theta_,
          independent.begin() + n_theta_ + n_sigma_);
        source_point.insert(
          source_point.end(), independent.begin() + n_theta_ + n_sigma_,
          independent.end());
        const AD probability = dynamic[cursor++];
        const std::vector<AD> value = nested.Forward(0, source_point, messages);
        if (value.size() != 1U) {
          throw std::logic_error(
            "A reduced weighted subject tape returned an invalid range.");
        }
        total += probability * value[0];
      }
    }
    if (cursor != dynamic.size()) {
      throw std::logic_error(
        "Reduced weighted population dynamic data were not fully consumed.");
    }
    std::vector<AD> dependent(1U, total);
    auto reduced = std::make_unique<ObjectiveTape>();
    reduced->fun.Dependent(independent, dependent);
    recording.release();
    reduced->fun.optimize();
    reduced->dynamic_values = dynamic_values;
    reduced->domain_names.reserve(point.size());
    for (int index = 0; index < n_theta_; ++index) {
      reduced->domain_names.push_back("THETA_" + std::to_string(index + 1));
    }
    for (int index = 0; index < n_sigma_; ++index) {
      reduced->domain_names.push_back("SIGMA_" + std::to_string(index + 1));
    }
    for (int index = 0; index < n_omega_; ++index) {
      reduced->domain_names.push_back("OMEGA_" + std::to_string(index + 1));
    }
    reduced_tape_ = std::move(reduced);
    reduced_shape_ = shape;
    reduced_dynamic_dirty_ = false;
    reduced_fallback_reason_.clear();
    ++reduced_records_;
  }

  SaemEvaluation evaluate_reduced_native(
      const std::vector<double>& theta, const std::vector<double>& sigma,
      const std::vector<double>& omega, bool gradient) {
    if (!reduced_tape_) {
      throw std::logic_error("The reduced weighted population tape is unavailable.");
    }
    const std::vector<double> point = reduced_point(theta, sigma, omega);
    std::ostringstream messages;
    const std::vector<double> value = reduced_tape_->fun.Forward(
      0, point, messages);
    require_unchanged_path(
      reduced_tape_->fun, "reduced weighted population objective");
    if (value.size() != 1U || !std::isfinite(value[0])) {
      throw std::domain_error(
        "The reduced weighted population objective is non-finite.");
    }
    SaemEvaluation result;
    result.value = value[0];
    result.gradient = Vector::Zero(static_cast<Eigen::Index>(domain_));
    if (gradient) {
      const std::vector<double> derivative = reduced_tape_->fun.Reverse(
        1, std::vector<double>(1U, 1.0));
      require_unchanged_path(
        reduced_tape_->fun, "reduced weighted population gradient");
      if (derivative.size() != point.size()) {
        throw std::logic_error(
          "The reduced weighted population gradient has the wrong length.");
      }
      for (int index = 0; index < n_theta_; ++index) {
        result.gradient[index] = derivative[static_cast<std::size_t>(index)];
      }
      const int source_sigma = n_theta_;
      const int target_sigma = n_theta_ + n_eta_;
      for (int index = 0; index < n_sigma_; ++index) {
        result.gradient[target_sigma + index] =
          derivative[static_cast<std::size_t>(source_sigma + index)];
      }
      const int source_omega = n_theta_ + n_sigma_;
      const int target_omega = n_theta_ + n_eta_ + n_sigma_;
      for (int index = 0; index < n_omega_; ++index) {
        result.gradient[target_omega + index] =
          derivative[static_cast<std::size_t>(source_omega + index)];
      }
    }
    ++evaluations_;
    ++reduced_evaluations_;
    last_point_evaluations_ = support_point_count();
    point_evaluations_ += last_point_evaluations_;
    return result;
  }

  static void normalize(std::vector<double>& weights) {
    const double total = std::accumulate(weights.begin(), weights.end(), 0.0);
    if (!(total > 0.0) || !std::isfinite(total)) {
      throw std::domain_error("Weighted ETA probabilities have zero mass.");
    }
    for (double& weight : weights) weight /= total;
  }

  void compress_common(int max_support, double tolerance) {
    const std::vector<double>& weights = common_weights_;
    if (weights.empty()) return;
    std::vector<int> keep;
    keep.reserve(weights.size());
    for (int index = 0; index < static_cast<int>(weights.size()); ++index) {
      if (weights[static_cast<std::size_t>(index)] >= tolerance) {
        keep.push_back(index);
      }
    }
    if (keep.empty()) {
      keep.push_back(static_cast<int>(std::distance(
        weights.begin(), std::max_element(weights.begin(), weights.end()))));
    }
    if (max_support > 0 && static_cast<int>(keep.size()) > max_support) {
      std::vector<int> ranked = keep;
      std::stable_sort(ranked.begin(), ranked.end(), [&](int left, int right) {
        return weights[static_cast<std::size_t>(left)] >
          weights[static_cast<std::size_t>(right)];
      });
      ranked.resize(static_cast<std::size_t>(max_support));
      std::sort(ranked.begin(), ranked.end());
      keep.swap(ranked);
    }
    if (keep.size() == weights.size()) {
      recompute_moments();
      return;
    }
    double retained = 0.0;
    for (int index : keep) retained += weights[static_cast<std::size_t>(index)];
    pruned_mass_ += std::max(0.0, 1.0 - retained);
    compressed_points_ += static_cast<long long>(weights.size() - keep.size());
    std::vector<double> compact_weights;
    compact_weights.reserve(keep.size());
    for (int index : keep) {
      compact_weights.push_back(weights[static_cast<std::size_t>(index)]);
    }
    normalize(compact_weights);
    for (std::size_t subject = 0; subject < grids_.size(); ++subject) {
      Matrix compact(static_cast<Eigen::Index>(keep.size()), n_eta_);
      for (std::size_t row = 0; row < keep.size(); ++row) {
        compact.row(static_cast<Eigen::Index>(row)) =
          grids_[subject].row(keep[row]);
      }
      grids_[subject].swap(compact);
      weights_[subject].clear();
    }
    common_weights_.swap(compact_weights);
    recompute_moments();
  }

  void recompute_moments() {
    mean_eta_.setZero();
    for (Matrix& moment : second_moment_) moment.setZero();
    for (std::size_t subject = 0; subject < grids_.size(); ++subject) {
      for (Eigen::Index support = 0;
           support < subject_support_count(subject); ++support) {
        const double probability = probability_at(subject, support);
        const Vector eta = grids_[subject].row(support).transpose();
        mean_eta_.row(static_cast<Eigen::Index>(subject)) +=
          probability * eta.transpose();
        second_moment_[subject].noalias() +=
          probability * eta * eta.transpose();
      }
    }
  }

  void validate_parameters(
      const Rcpp::NumericVector& theta, const Rcpp::NumericVector& sigma,
      const Rcpp::NumericVector& omega) const {
    validate_parameters_native(
      Rcpp::as<std::vector<double>>(theta),
      Rcpp::as<std::vector<double>>(sigma),
      Rcpp::as<std::vector<double>>(omega));
  }

  void validate_parameters_native(
      const std::vector<double>& theta, const std::vector<double>& sigma,
      const std::vector<double>& omega) const {
    if (static_cast<int>(theta.size()) != n_theta_ ||
        static_cast<int>(sigma.size()) != n_sigma_ ||
        static_cast<int>(omega.size()) != n_omega_) {
      throw std::invalid_argument(
        "Persistent weighted-ETA population parameter dimensions changed.");
    }
    for (double value : theta) if (!std::isfinite(value)) {
      throw std::invalid_argument("Weighted-ETA THETAs must be finite.");
    }
    for (double value : sigma) if (!std::isfinite(value)) {
      throw std::invalid_argument("Weighted-ETA SIGMAs must be finite.");
    }
    for (double value : omega) if (!std::isfinite(value)) {
      throw std::invalid_argument("Weighted-ETA OMEGAs must be finite.");
    }
  }
};

