struct EtaEvaluation {
  double value = std::numeric_limits<double>::infinity();
  Vector gradient;
  bool finite = false;
};

EtaEvaluation objective_eta_evaluate(
    ObjectiveTape& tape, const std::vector<double>& point,
    const std::vector<std::size_t>& positions, bool gradient = true) {
  std::ostringstream messages;
  const std::vector<double> value = tape.fun.Forward(0, point, messages);
  EtaEvaluation result;
  // An invalid line-search trial must be rejected by the optimizer, not used
  // as a retaping anchor. At extreme ETAs an exponentiated rate can underflow
  // to zero; replay then becomes non-finite and comparison changes at that
  // same point are immaterial. Only finite trials are eligible for retaping.
  if (value.empty() || !std::isfinite(value[0])) return result;
  if (tape.fun.compare_change_number() != 0U) {
    throw TapePathChange("conditional objective", point);
  }
  result.value = value[0];
  result.finite = true;
  result.gradient = Vector::Zero(static_cast<Eigen::Index>(positions.size()));
  if (!gradient || positions.empty()) return result;
  const std::vector<double> weight(1, 1.0);
  const std::vector<double> full = tape.fun.Reverse(1, weight);
  for (std::size_t i = 0; i < positions.size(); ++i) {
    result.gradient[static_cast<Eigen::Index>(i)] = full[positions[i]];
    if (!std::isfinite(result.gradient[static_cast<Eigen::Index>(i)])) {
      result.finite = false;
    }
  }
  return result;
}

Matrix objective_eta_hessian(
    ObjectiveTape& tape, const std::vector<double>& point,
    const std::vector<std::size_t>& positions) {
  const std::size_t domain = tape.domain_names.size();
  const std::size_t dimension = positions.size();
  Matrix hessian = Matrix::Zero(
    static_cast<Eigen::Index>(dimension), static_cast<Eigen::Index>(dimension));
  std::ostringstream messages;
  tape.fun.Forward(0, point, messages);
  require_unchanged_path(tape.fun, "conditional objective Hessian");
  const std::vector<double> weight(1, 1.0);
  std::vector<double> direction(domain, 0.0);
  for (std::size_t column = 0; column < dimension; ++column) {
    direction[positions[column]] = 1.0;
    tape.fun.Forward(1, direction, messages);
    direction[positions[column]] = 0.0;
    const std::vector<double> reverse = tape.fun.Reverse(2, weight);
    for (std::size_t row = 0; row < dimension; ++row) {
      hessian(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) =
        reverse[positions[row] * 2U + 1U];
    }
  }
  return 0.5 * (hessian + hessian.transpose()).eval();
}

Rcpp::List objective_eta_mode(
    ObjectiveTape& tape, std::vector<double> point,
    const std::vector<std::size_t>& positions,
    const Rcpp::NumericVector& start, int maxit, double tolerance,
    bool exact_hessian) {
  const Eigen::Index dimension = static_cast<Eigen::Index>(positions.size());
  if (start.size() != dimension) {
    throw std::invalid_argument("ETA starting point has the wrong length.");
  }
  if (maxit < 1 || !std::isfinite(tolerance) || tolerance <= 0.0) {
    throw std::invalid_argument("ETA optimizer controls are invalid.");
  }
  Vector eta(dimension);
  for (Eigen::Index i = 0; i < dimension; ++i) {
    eta[i] = start[i];
    point[positions[static_cast<std::size_t>(i)]] = eta[i];
  }
  EtaEvaluation current = objective_eta_evaluate(tape, point, positions, true);
  int evaluations = 1;
  int gradient_evaluations = 1;
  int convergence = current.finite ? 1 : 52;
  int iterations = 0;
  Matrix inverse = Matrix::Identity(dimension, dimension);
  const Matrix identity = Matrix::Identity(dimension, dimension);
  // Require a gradient-based stop so that warm starts do not change the
  // conditional objective through an early relative-function-value stop.
  const double gradient_tolerance = std::max(1e-8, tolerance);

  for (int iteration = 0; current.finite && iteration < maxit; ++iteration) {
    iterations = iteration;
    if (current.gradient.lpNorm<Eigen::Infinity>() <= gradient_tolerance) {
      convergence = 0;
      break;
    }
    Vector direction = -inverse * current.gradient;
    double directional = current.gradient.dot(direction);
    if (!std::isfinite(directional) || directional >= -1e-14) {
      inverse.setIdentity();
      direction = -current.gradient;
      directional = -current.gradient.squaredNorm();
    }
    double step_scale = 1.0;
    EtaEvaluation candidate;
    Vector candidate_eta(dimension);
    std::vector<double> candidate_point;
    bool accepted = false;
    for (int line_search = 0; line_search < 32; ++line_search) {
      candidate_eta = eta + step_scale * direction;
      candidate_point = point;
      for (Eigen::Index i = 0; i < dimension; ++i) {
        candidate_point[positions[static_cast<std::size_t>(i)]] = candidate_eta[i];
      }
      candidate = objective_eta_evaluate(tape, candidate_point, positions, false);
      ++evaluations;
      if (candidate.finite &&
          candidate.value <= current.value + 1e-4 * step_scale * directional) {
        accepted = true;
        break;
      }
      step_scale *= 0.5;
    }
    if (!accepted) {
      convergence = 52;
      break;
    }
    candidate = objective_eta_evaluate(tape, candidate_point, positions, true);
    ++evaluations;
    ++gradient_evaluations;
    if (!candidate.finite) {
      convergence = 52;
      break;
    }
    const Vector displacement = candidate_eta - eta;
    const Vector gradient_change = candidate.gradient - current.gradient;
    const double curvature = gradient_change.dot(displacement);
    if (std::isfinite(curvature) &&
        curvature > 1e-12 * displacement.norm() * gradient_change.norm()) {
      const double rho = 1.0 / curvature;
      const Matrix left = identity - rho * displacement * gradient_change.transpose();
      inverse = left * inverse * left.transpose() +
        rho * displacement * displacement.transpose();
    } else {
      inverse.setIdentity();
    }
    eta = candidate_eta;
    point.swap(candidate_point);
    current = std::move(candidate);
    iterations = iteration + 1;
    if (current.gradient.lpNorm<Eigen::Infinity>() <= gradient_tolerance) {
      convergence = 0;
      break;
    }
    if ((iteration + 1) % 20 == 0) Rcpp::checkUserInterrupt();
  }

  if (current.finite && convergence != 0 &&
      current.gradient.lpNorm<Eigen::Infinity>() <= 10.0 * gradient_tolerance) {
    convergence = 0;
  }
  Matrix hessian = exact_hessian ? objective_eta_hessian(tape, point, positions) :
    Matrix::Zero(0, 0);
  Rcpp::NumericVector par(dimension);
  Rcpp::NumericVector gradient(dimension);
  for (Eigen::Index i = 0; i < dimension; ++i) {
    par[i] = eta[i];
    gradient[i] = current.gradient[i];
  }
  return Rcpp::List::create(
    Rcpp::Named("par") = par,
    Rcpp::Named("value") = current.value,
    Rcpp::Named("convergence") = convergence,
    Rcpp::Named("hessian") = libertad::eigen_matrix_to_r(hessian),
    Rcpp::Named("gradient") = gradient,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("evaluations") = evaluations,
    Rcpp::Named("gradient_evaluations") = gradient_evaluations
  );
}

struct PopulationParameters {
  std::vector<double> theta;
  std::vector<double> sigma;
  std::vector<double> omega;
  Matrix transform;
};

struct PopulationPrior {
  int native_index = -1;
  std::string family;
  double mean = 0.0;
  double sd = 1.0;
  double shape = std::numeric_limits<double>::quiet_NaN();
  double rate = std::numeric_limits<double>::quiet_NaN();
};

// Persistent population objective used by R's mature L-BFGS-B/BFGS driver.
// The R callbacks around this object only transfer one encoded parameter
// vector. Parameter decoding, conditional modes, curvature, priors, AD
// derivatives, and same-point caches all remain in this compiled object.
class PopulationObjective {
 public:
  PopulationObjective(
      SEXP engine_pointer, const Rcpp::List& subject_data,
      const Rcpp::List& primary_tape_pointers,
      const Rcpp::List& curvature_tape_pointers,
      const Rcpp::List& config) {
    Rcpp::XPtr<ModelEngine> engine(engine_pointer);
    engine_ = engine.get();
    approximation_ = Rcpp::as<std::string>(config["approximation"]);
    if (approximation_ != "fo" && approximation_ != "its" &&
        approximation_ != "foce" && approximation_ != "focei" &&
        approximation_ != "laplace") {
      throw std::invalid_argument("Unknown compiled population approximation.");
    }
    theta_base_ = Rcpp::as<std::vector<double>>(config["theta"]);
    sigma_base_ = Rcpp::as<std::vector<double>>(config["sigma"]);
    omega_base_ = Rcpp::as<std::vector<double>>(config["omega"]);
    theta_free_ = zero_based(Rcpp::as<std::vector<int>>(config["theta_free"]));
    sigma_free_ = zero_based(Rcpp::as<std::vector<int>>(config["sigma_free"]));
    omega_free_ = zero_based(Rcpp::as<std::vector<int>>(config["omega_free"]));
    omega_full_ = Rcpp::as<bool>(config["omega_full"]);
    omega_rows_ = zero_based(Rcpp::as<std::vector<int>>(config["omega_rows"]));
    omega_cols_ = zero_based(Rcpp::as<std::vector<int>>(config["omega_cols"]));
    n_eta_ = Rcpp::as<int>(config["n_eta"]);
    n_eta_base_ = Rcpp::as<int>(config["n_eta_base"]);
    eta_maxit_ = Rcpp::as<int>(config["eta_maxit"]);
    tolerance_ = Rcpp::as<double>(config["tolerance"]);
    use_ode_ = Rcpp::as<bool>(config["use_ode"]);
    fo_population_batch_requested_ = Rcpp::as<bool>(config["fo_population_batch"]);
    fo_population_max_operations_ =
      Rcpp::as<double>(config["fo_population_max_operations"]);
    guard_radius_ = Rcpp::as<double>(config["guard_radius"]);
    start_ = Rcpp::as<std::vector<double>>(config["start"]);
    if (eta_maxit_ < 1 || tolerance_ <= 0.0 || !std::isfinite(tolerance_) ||
        guard_radius_ <= 0.0 || !std::isfinite(guard_radius_) ||
        fo_population_max_operations_ <= 0.0 ||
        !std::isfinite(fo_population_max_operations_)) {
      throw std::invalid_argument("Compiled population controls are invalid.");
    }
    if (omega_rows_.size() != omega_base_.size() ||
        omega_cols_.size() != omega_base_.size()) {
      throw std::invalid_argument("Compiled OMEGA mapping is inconsistent.");
    }
    const std::vector<int> prior_index = zero_based(
      Rcpp::as<std::vector<int>>(config["prior_index"]));
    const std::vector<std::string> prior_family =
      Rcpp::as<std::vector<std::string>>(config["prior_family"]);
    const std::vector<double> prior_mean =
      Rcpp::as<std::vector<double>>(config["prior_mean"]);
    const std::vector<double> prior_sd =
      Rcpp::as<std::vector<double>>(config["prior_sd"]);
    const std::vector<double> prior_shape =
      Rcpp::as<std::vector<double>>(config["prior_shape"]);
    const std::vector<double> prior_rate =
      Rcpp::as<std::vector<double>>(config["prior_rate"]);
    const std::size_t prior_count = prior_index.size();
    if (prior_family.size() != prior_count || prior_mean.size() != prior_count ||
        prior_sd.size() != prior_count || prior_shape.size() != prior_count ||
        prior_rate.size() != prior_count) {
      throw std::invalid_argument("Compiled prior mapping is inconsistent.");
    }
    priors_.reserve(prior_count);
    for (std::size_t index = 0; index < prior_count; ++index) {
      priors_.push_back(PopulationPrior{
        prior_index[index], prior_family[index], prior_mean[index],
        prior_sd[index], prior_shape[index], prior_rate[index]
      });
    }

    const int subjects = subject_data.size();
    if (subjects < 1) throw std::invalid_argument("Population data have no subjects.");
    subject_data_.reserve(static_cast<std::size_t>(subjects));
    for (int subject = 0; subject < subjects; ++subject) {
      subject_data_.push_back(subject_data[subject]);
    }
    starts_ = Matrix::Zero(subjects, n_eta_);
    if (config.containsElementNamed("eta_start")) {
      Rcpp::NumericMatrix eta_start = config["eta_start"];
      if (eta_start.nrow() != subjects || eta_start.ncol() != n_eta_) {
        throw std::invalid_argument("Compiled ETA start matrix has the wrong dimensions.");
      }
      for (int subject = 0; subject < subjects; ++subject) {
        for (int effect = 0; effect < n_eta_; ++effect) {
          const double value = eta_start(subject, effect);
          if (!std::isfinite(value)) {
            throw std::invalid_argument("Compiled ETA starts must be finite.");
          }
          starts_(subject, effect) = value;
        }
      }
    }
    primary_.resize(static_cast<std::size_t>(subjects), nullptr);
    curvature_.resize(static_cast<std::size_t>(subjects), nullptr);
    owned_prediction_.resize(static_cast<std::size_t>(subjects));
    owned_primary_.resize(static_cast<std::size_t>(subjects));
    owned_curvature_.resize(static_cast<std::size_t>(subjects));
    anchors_.resize(static_cast<std::size_t>(subjects));

    if (use_ode_) {
      const PopulationParameters initial = decode(start_);
      for (int subject = 0; subject < subjects; ++subject) {
        record_subject(subject, initial, starts_.row(subject).transpose(), false);
      }
    } else {
      if (primary_tape_pointers.size() != subjects) {
        throw std::invalid_argument("Population tapes do not match subject data.");
      }
      for (int subject = 0; subject < subjects; ++subject) {
        SEXP source = primary_tape_pointers[subject];
        Rcpp::XPtr<ObjectiveTape> tape(source);
        primary_[static_cast<std::size_t>(subject)] = tape.get();
      }
      if (has_curvature()) {
        if (curvature_tape_pointers.size() != subjects) {
          throw std::invalid_argument("Curvature tapes do not match subject data.");
        }
        for (int subject = 0; subject < subjects; ++subject) {
          SEXP source = curvature_tape_pointers[subject];
          Rcpp::XPtr<ObjectiveTape> tape(source);
          curvature_[static_cast<std::size_t>(subject)] = tape.get();
        }
      }
    }
    if (is_fo()) {
      std::unordered_map<ObjectiveTape*, bool> unique;
      for (ObjectiveTape* tape : primary_) unique[tape] = true;
      fo_unique_subject_tapes_ = static_cast<int>(unique.size());
      if (!use_ode_ && fo_population_batch_requested_) {
        try {
          record_fo_population(decode(start_));
        } catch (const std::exception& error) {
          fo_population_error_ = error.what();
          fo_population_.reset();
        }
      }
    }
  }

  double value(const Rcpp::NumericVector& encoded) {
    ++value_requests_;
    const std::vector<double> point = Rcpp::as<std::vector<double>>(encoded);
    if (same_key(point)) {
      ++value_cache_hits_;
      return cache_value_;
    }
    evaluate_value(point);
    return cache_value_;
  }

  Rcpp::NumericVector gradient(const Rcpp::NumericVector& encoded) {
    ++gradient_requests_;
    const std::vector<double> point = Rcpp::as<std::vector<double>>(encoded);
    const bool reused_value = same_key(point);
    if (!reused_value) evaluate_value(point);
    else ++shared_state_hits_;
    if (cache_gradient_valid_) {
      ++gradient_cache_hits_;
      return Rcpp::wrap(cache_gradient_);
    }
    if (!std::isfinite(cache_value_) || cache_value_ >= penalty()) {
      throw std::runtime_error("Cannot differentiate a failed population objective.");
    }
    evaluate_gradient();
    return Rcpp::wrap(cache_gradient_);
  }

  Rcpp::NumericMatrix hessian(const Rcpp::NumericVector& encoded) {
    ++hessian_requests_;
    const std::vector<double> point = Rcpp::as<std::vector<double>>(encoded);
    if (!same_key(point)) evaluate_value(point);
    if (!std::isfinite(cache_value_) || cache_value_ >= penalty()) {
      throw std::runtime_error("Cannot differentiate a failed population objective.");
    }
    if (std::any_of(
          cache_mode_convergence_.begin(), cache_mode_convergence_.end(),
          [](int value) { return value != 0; })) {
      throw std::runtime_error(
        "An exact population Hessian requires converged conditional modes.");
    }

    const int n_theta = static_cast<int>(cache_parameters_.theta.size());
    const int n_sigma = static_cast<int>(cache_parameters_.sigma.size());
    const int n_omega = static_cast<int>(cache_parameters_.omega.size());
    const int n_native = n_theta + n_sigma + n_omega;
    const int n_outer = static_cast<int>(point.size());
    Vector native_gradient = Vector::Zero(n_native);
    Matrix native_hessian = Matrix::Zero(n_native, n_native);
    Vector prior_gradient;
    prior_nll(cache_parameters_, &prior_gradient);
    native_gradient += prior_gradient;
    native_hessian += prior_hessian(cache_parameters_);

    std::vector<int> population;
    population.reserve(static_cast<std::size_t>(n_native));
    for (int index = 0; index < n_theta; ++index) population.push_back(index);
    for (int index = 0; index < n_sigma; ++index) {
      population.push_back(n_theta + n_eta_ + index);
    }
    for (int index = 0; index < n_omega; ++index) {
      population.push_back(n_theta + n_eta_ + n_sigma + index);
    }

    if (is_fo()) {
      const std::vector<double> native_point = fo_point(cache_parameters_);
      for (std::size_t subject = 0; subject < primary_.size(); ++subject) {
        ObjectiveTape& objective = *primary_[subject];
        Rcpp::DataFrame data(subject_data_[subject]);
        set_fo_dynamic(objective, data);
        ++fo_dynamic_updates_;
        const Vector gradient = tape_gradient(
          objective, native_point, "FO exact population Hessian");
        const Matrix hessian = tape_hessian(
          objective, native_point, "FO exact population Hessian");
        native_gradient += gradient;
        native_hessian += hessian;
      }
    } else {
      std::vector<int> eta_positions(static_cast<std::size_t>(n_eta_));
      for (int effect = 0; effect < n_eta_; ++effect) {
        eta_positions[static_cast<std::size_t>(effect)] = n_theta + effect;
      }
      for (std::size_t subject = 0; subject < primary_.size(); ++subject) {
        exact_subject_marginal_curvature(
          *primary_[subject],
          has_curvature() ? curvature_[subject] : nullptr,
          cache_points_[subject], population, eta_positions,
          native_gradient, native_hessian);
      }
    }

    CppAD::ADFun<double> transform_fun = native_transform_tape(point);
    const std::vector<double> jacobian_values = transform_fun.Jacobian(point);
    if (jacobian_values.size() !=
        static_cast<std::size_t>(n_native * n_outer)) {
      throw std::logic_error("The population transform Jacobian has the wrong size.");
    }
    Matrix jacobian(n_native, n_outer);
    for (int row = 0; row < n_native; ++row) {
      for (int column = 0; column < n_outer; ++column) {
        jacobian(row, column) = jacobian_values[
          static_cast<std::size_t>(row * n_outer + column)];
      }
    }
    Matrix outer_hessian = jacobian.transpose() * native_hessian * jacobian;
    std::vector<double> weight(static_cast<std::size_t>(n_native), 0.0);
    for (int native = 0; native < n_native; ++native) {
      if (native_gradient[native] == 0.0) continue;
      weight[static_cast<std::size_t>(native)] = 1.0;
      const Matrix second = adfun_hessian(transform_fun, point, weight);
      outer_hessian += native_gradient[native] * second;
      weight[static_cast<std::size_t>(native)] = 0.0;
    }
    outer_hessian = 0.5 * (outer_hessian + outer_hessian.transpose()).eval();
    if (!outer_hessian.allFinite()) {
      throw std::runtime_error("The exact population Hessian is not finite.");
    }
    return libertad::eigen_matrix_to_r(outer_hessian);
  }

  Rcpp::List state(const Rcpp::NumericVector& encoded) {
    const std::vector<double> point = Rcpp::as<std::vector<double>>(encoded);
    if (!same_key(point)) evaluate_value(point);
    Rcpp::NumericMatrix eta(starts_.rows(), starts_.cols());
    for (Eigen::Index row = 0; row < starts_.rows(); ++row) {
      for (Eigen::Index column = 0; column < starts_.cols(); ++column) {
        eta(row, column) = starts_(row, column);
      }
    }
    Rcpp::List modes(starts_.rows());
    for (Eigen::Index subject = 0; subject < starts_.rows(); ++subject) {
      Rcpp::NumericVector par(starts_.cols());
      for (Eigen::Index effect = 0; effect < starts_.cols(); ++effect) {
        par[effect] = starts_(subject, effect);
      }
      modes[subject] = Rcpp::List::create(
        Rcpp::Named("par") = par,
        Rcpp::Named("value") = cache_subject_values_.empty() ? NA_REAL :
          cache_subject_values_[static_cast<std::size_t>(subject)],
        Rcpp::Named("convergence") = cache_mode_convergence_.empty() ? 0 :
          cache_mode_convergence_[static_cast<std::size_t>(subject)],
        Rcpp::Named("hessian") = Rcpp::NumericMatrix(0, 0),
        Rcpp::Named("logdet") = cache_curvature_values_.empty() ? 0.0 :
          cache_curvature_values_[static_cast<std::size_t>(subject)],
        Rcpp::Named("jitter") = 0.0,
        Rcpp::Named("gradient") = Rcpp::NumericVector(starts_.cols()),
        Rcpp::Named("iterations") = 0,
        Rcpp::Named("evaluations") = 0,
        Rcpp::Named("backend") = "cpp-population-cache"
      );
    }
    return Rcpp::List::create(
      Rcpp::Named("eta") = eta,
      Rcpp::Named("modes") = modes,
      Rcpp::Named("value") = cache_value_
    );
  }

  Rcpp::List telemetry() const {
    return Rcpp::List::create(
      Rcpp::Named("backend") = fo_population_ ?
        "persistent-cpp-batched-fo-population-objective" :
        "persistent-cpp-population-objective",
      Rcpp::Named("approximation") = approximation_,
      Rcpp::Named("value_requests") = value_requests_,
      Rcpp::Named("gradient_requests") = gradient_requests_,
      Rcpp::Named("hessian_requests") = hessian_requests_,
      Rcpp::Named("parameter_evaluations") = parameter_evaluations_,
      Rcpp::Named("value_cache_hits") = value_cache_hits_,
      Rcpp::Named("gradient_cache_hits") = gradient_cache_hits_,
      Rcpp::Named("shared_state_hits") = shared_state_hits_,
      Rcpp::Named("mode_iterations") = mode_iterations_,
      Rcpp::Named("mode_evaluations") = mode_evaluations_,
      Rcpp::Named("mode_recoveries") = mode_recoveries_,
      Rcpp::Named("tape_records") = tape_records_,
      Rcpp::Named("tape_retapes") = tape_retapes_,
      Rcpp::Named("ode_owned_tapes") = use_ode_,
      Rcpp::Named("fo_unique_subject_tapes") = fo_unique_subject_tapes_,
      Rcpp::Named("fo_shared_subject_tapes") = is_fo() ?
        static_cast<int>(primary_.size()) - fo_unique_subject_tapes_ : 0,
      Rcpp::Named("fo_population_batched") = static_cast<bool>(fo_population_),
      Rcpp::Named("fo_population_operations") = fo_population_ ?
        static_cast<double>(fo_population_->fun.size_op()) : 0.0,
      Rcpp::Named("fo_population_fallbacks") = fo_population_fallbacks_,
      Rcpp::Named("fo_dynamic_updates") = static_cast<double>(fo_dynamic_updates_),
      Rcpp::Named("fo_population_error") = fo_population_error_,
      Rcpp::Named("propagation_kernel") = propagation_kernel_name(*engine_)
    );
  }

 private:
  ModelEngine* engine_ = nullptr;
  std::string approximation_;
  std::vector<double> theta_base_, sigma_base_, omega_base_, start_;
  std::vector<int> theta_free_, sigma_free_, omega_free_;
  std::vector<int> omega_rows_, omega_cols_;
  bool omega_full_ = false;
  bool use_ode_ = false;
  bool fo_population_batch_requested_ = true;
  double fo_population_max_operations_ = 2e6;
  int n_eta_ = 0;
  int n_eta_base_ = 0;
  int eta_maxit_ = 100;
  double tolerance_ = 1e-7;
  double guard_radius_ = 0.5;
  std::vector<PopulationPrior> priors_;
  std::vector<SEXP> subject_data_;
  std::vector<ObjectiveTape*> primary_, curvature_;
  std::vector<std::unique_ptr<PredictionTape>> owned_prediction_;
  std::vector<std::unique_ptr<ObjectiveTape>> owned_primary_, owned_curvature_;
  std::unique_ptr<ObjectiveTape> fo_population_;
  std::vector<std::vector<double>> anchors_;
  Matrix starts_;

  bool cache_valid_ = false;
  bool cache_gradient_valid_ = false;
  std::vector<double> cache_key_, cache_gradient_, cache_subject_values_;
  std::vector<double> cache_curvature_values_;
  std::vector<int> cache_mode_convergence_;
  std::vector<std::vector<double>> cache_points_;
  PopulationParameters cache_parameters_;
  double cache_value_ = std::numeric_limits<double>::infinity();

  int value_requests_ = 0;
  int gradient_requests_ = 0;
  int hessian_requests_ = 0;
  int parameter_evaluations_ = 0;
  int value_cache_hits_ = 0;
  int gradient_cache_hits_ = 0;
  int shared_state_hits_ = 0;
  long long mode_iterations_ = 0;
  long long mode_evaluations_ = 0;
  int mode_recoveries_ = 0;
  int tape_records_ = 0;
  int tape_retapes_ = 0;
  int fo_unique_subject_tapes_ = 0;
  int fo_population_fallbacks_ = 0;
  long long fo_dynamic_updates_ = 0;
  std::string fo_population_error_;

  static double penalty() { return 1e100; }

  static std::vector<int> zero_based(std::vector<int> source) {
    for (int& value : source) {
      if (value < 1) throw std::invalid_argument("A parameter index is invalid.");
      --value;
    }
    return source;
  }

  bool is_fo() const { return approximation_ == "fo"; }
  bool has_curvature() const {
    return approximation_ == "foce" || approximation_ == "focei" ||
      approximation_ == "laplace";
  }
  bool interaction() const { return approximation_ != "foce"; }

  bool same_key(const std::vector<double>& point) const {
    return cache_valid_ && point.size() == cache_key_.size() &&
      std::equal(point.begin(), point.end(), cache_key_.begin());
  }

  PopulationParameters decode(const std::vector<double>& encoded) const {
    PopulationParameters result;
    result.theta = theta_base_;
    result.sigma = sigma_base_;
    result.omega = omega_base_;
    const int n_native = static_cast<int>(
      result.theta.size() + result.sigma.size() + result.omega.size());
    result.transform = Matrix::Zero(n_native, encoded.size());
    std::size_t cursor = 0;
    for (int index : theta_free_) {
      if (cursor >= encoded.size() || index < 0 ||
          index >= static_cast<int>(result.theta.size())) {
        throw std::invalid_argument("Encoded THETA mapping is invalid.");
      }
      result.theta[static_cast<std::size_t>(index)] = encoded[cursor];
      result.transform(index, static_cast<Eigen::Index>(cursor)) = 1.0;
      ++cursor;
    }
    const int sigma_offset = static_cast<int>(result.theta.size());
    for (int index : sigma_free_) {
      if (cursor >= encoded.size() || index < 0 ||
          index >= static_cast<int>(result.sigma.size())) {
        throw std::invalid_argument("Encoded SIGMA mapping is invalid.");
      }
      const double value = std::exp(encoded[cursor]);
      result.sigma[static_cast<std::size_t>(index)] = value;
      result.transform(sigma_offset + index, static_cast<Eigen::Index>(cursor)) = value;
      ++cursor;
    }
    const int omega_offset = sigma_offset + static_cast<int>(result.sigma.size());
    if (omega_full_ && !omega_free_.empty()) {
      if (n_eta_base_ < 1 || result.omega.size() != omega_rows_.size() ||
          cursor + result.omega.size() > encoded.size()) {
        throw std::invalid_argument("Encoded full OMEGA mapping is invalid.");
      }
      Matrix lower = Matrix::Zero(n_eta_base_, n_eta_base_);
      for (std::size_t entry = 0; entry < result.omega.size(); ++entry) {
        const int row = omega_rows_[entry];
        const int column = omega_cols_[entry];
        if (row < 0 || column < 0 || row >= n_eta_base_ || column > row) {
          throw std::invalid_argument("OMEGA Cholesky coordinates are invalid.");
        }
        lower(row, column) = row == column ? std::exp(encoded[cursor + entry]) :
          encoded[cursor + entry];
      }
      const Matrix covariance = lower * lower.transpose();
      for (std::size_t entry = 0; entry < result.omega.size(); ++entry) {
        result.omega[entry] = covariance(omega_rows_[entry], omega_cols_[entry]);
      }
      for (std::size_t encoded_entry = 0; encoded_entry < result.omega.size();
           ++encoded_entry) {
        Matrix derivative_lower = Matrix::Zero(n_eta_base_, n_eta_base_);
        const int row = omega_rows_[encoded_entry];
        const int column = omega_cols_[encoded_entry];
        derivative_lower(row, column) = row == column ? lower(row, column) : 1.0;
        const Matrix derivative = derivative_lower * lower.transpose() +
          lower * derivative_lower.transpose();
        for (std::size_t native = 0; native < result.omega.size(); ++native) {
          result.transform(
            omega_offset + static_cast<int>(native),
            static_cast<Eigen::Index>(cursor + encoded_entry)) = derivative(
              omega_rows_[native], omega_cols_[native]);
        }
      }
      cursor += result.omega.size();
    } else {
      for (int index : omega_free_) {
        if (cursor >= encoded.size() || index < 0 ||
            index >= static_cast<int>(result.omega.size())) {
          throw std::invalid_argument("Encoded OMEGA mapping is invalid.");
        }
        const double value = std::exp(encoded[cursor]);
        result.omega[static_cast<std::size_t>(index)] = value;
        result.transform(omega_offset + index, static_cast<Eigen::Index>(cursor)) = value;
        ++cursor;
      }
    }
    if (cursor != encoded.size()) {
      throw std::invalid_argument("Encoded population parameter length is invalid.");
    }
    return result;
  }

  std::vector<double> native_values(const PopulationParameters& parameters) const {
    std::vector<double> result;
    result.reserve(parameters.theta.size() + parameters.sigma.size() +
                   parameters.omega.size());
    result.insert(result.end(), parameters.theta.begin(), parameters.theta.end());
    result.insert(result.end(), parameters.sigma.begin(), parameters.sigma.end());
    result.insert(result.end(), parameters.omega.begin(), parameters.omega.end());
    return result;
  }

  double prior_nll(const PopulationParameters& parameters,
                   Vector* derivative = nullptr) const {
    const std::vector<double> native = native_values(parameters);
    if (derivative) derivative->setZero(static_cast<Eigen::Index>(native.size()));
    const double log_two_pi = std::log(2.0 * std::acos(-1.0));
    double log_density = 0.0;
    for (const PopulationPrior& prior : priors_) {
      if (prior.native_index < 0 ||
          prior.native_index >= static_cast<int>(native.size())) {
        throw std::invalid_argument("A prior refers to an invalid parameter.");
      }
      const double value = native[static_cast<std::size_t>(prior.native_index)];
      double density = -std::numeric_limits<double>::infinity();
      double gradient = std::numeric_limits<double>::quiet_NaN();
      if (prior.family == "normal" || prior.family == "half_normal") {
        if (prior.sd > 0.0 && std::isfinite(value) &&
            (prior.family != "half_normal" || value >= 0.0)) {
          const double z = (value - prior.mean) / prior.sd;
          density = -0.5 * log_two_pi - std::log(prior.sd) - 0.5 * z * z;
          if (prior.family == "half_normal") density += std::log(2.0);
          gradient = 2.0 * (value - prior.mean) / (prior.sd * prior.sd);
        }
      } else if (prior.family == "lognormal") {
        if (value > 0.0 && prior.sd > 0.0) {
          const double z = (std::log(value) - prior.mean) / prior.sd;
          density = -std::log(value) - 0.5 * log_two_pi -
            std::log(prior.sd) - 0.5 * z * z;
          gradient = 2.0 / value + 2.0 * (std::log(value) - prior.mean) /
            (prior.sd * prior.sd * value);
        }
      } else if (prior.family == "inverse_gamma") {
        if (value > 0.0 && prior.shape > 0.0 && prior.rate > 0.0) {
          density = prior.shape * std::log(prior.rate) - std::lgamma(prior.shape) -
            (prior.shape + 1.0) * std::log(value) - prior.rate / value;
          gradient = 2.0 * (prior.shape + 1.0) / value -
            2.0 * prior.rate / (value * value);
        }
      } else {
        throw std::invalid_argument("Unknown compiled prior family.");
      }
      if (!std::isfinite(density)) return penalty();
      log_density += density;
      if (derivative) (*derivative)[prior.native_index] += gradient;
    }
    return -2.0 * log_density;
  }

  Matrix prior_hessian(const PopulationParameters& parameters) const {
    const std::vector<double> native = native_values(parameters);
    Matrix result = Matrix::Zero(
      static_cast<Eigen::Index>(native.size()),
      static_cast<Eigen::Index>(native.size()));
    for (const PopulationPrior& prior : priors_) {
      if (prior.native_index < 0 ||
          prior.native_index >= static_cast<int>(native.size())) {
        throw std::invalid_argument("A prior refers to an invalid parameter.");
      }
      const double value = native[static_cast<std::size_t>(prior.native_index)];
      double second = std::numeric_limits<double>::quiet_NaN();
      if (prior.family == "normal" || prior.family == "half_normal") {
        second = 2.0 / (prior.sd * prior.sd);
      } else if (prior.family == "lognormal") {
        const double centered = std::log(value) - prior.mean;
        second = -2.0 / (value * value) +
          2.0 * (1.0 - centered) /
            (prior.sd * prior.sd * value * value);
      } else if (prior.family == "inverse_gamma") {
        second = -2.0 * (prior.shape + 1.0) / (value * value) +
          4.0 * prior.rate / (value * value * value);
      } else {
        throw std::invalid_argument("Unknown compiled prior family.");
      }
      if (!std::isfinite(second)) {
        throw std::runtime_error("A prior Hessian contribution is not finite.");
      }
      result(prior.native_index, prior.native_index) += second;
    }
    return result;
  }

  static Matrix adfun_hessian(
      CppAD::ADFun<double>& fun, const std::vector<double>& point,
      const std::vector<double>& weight) {
    const std::vector<double> values = fun.Hessian(point, weight);
    const int dimension = static_cast<int>(point.size());
    if (values.size() != static_cast<std::size_t>(dimension * dimension)) {
      throw std::logic_error("A CppAD Hessian has the wrong size.");
    }
    Matrix result(dimension, dimension);
    for (int row = 0; row < dimension; ++row) {
      for (int column = 0; column < dimension; ++column) {
        result(row, column) = values[
          static_cast<std::size_t>(row * dimension + column)];
      }
    }
    return 0.5 * (result + result.transpose()).eval();
  }

  static Vector tape_gradient(
      ObjectiveTape& tape, const std::vector<double>& point,
      const std::string& context) {
    if (tape.fun.Range() != 1U || point.size() != tape.fun.Domain()) {
      throw std::invalid_argument("A scalar objective tape has incompatible dimensions.");
    }
    std::ostringstream messages;
    tape.fun.Forward(0, point, messages);
    require_unchanged_path(tape.fun, context);
    const std::vector<double> weight(1, 1.0);
    const std::vector<double> derivative = tape.fun.Reverse(1, weight);
    require_unchanged_path(tape.fun, context);
    Vector result(static_cast<Eigen::Index>(derivative.size()));
    for (std::size_t index = 0; index < derivative.size(); ++index) {
      result[static_cast<Eigen::Index>(index)] = derivative[index];
    }
    return result;
  }

  static Matrix tape_hessian(
      ObjectiveTape& tape, const std::vector<double>& point,
      const std::string& context) {
    if (tape.fun.Range() != 1U || point.size() != tape.fun.Domain()) {
      throw std::invalid_argument("A scalar objective tape has incompatible dimensions.");
    }
    const std::vector<double> weight(1, 1.0);
    const Matrix result = adfun_hessian(tape.fun, point, weight);
    require_unchanged_path(tape.fun, context);
    return result;
  }

  static CppAD::ADFun<double> mode_gradient_tape(
      ObjectiveTape& objective, const std::vector<double>& point,
      const std::vector<int>& eta_positions) {
    using AD = CppAD::AD<double>;
    std::vector<AD> independent(point.begin(), point.end());
    CppAD::Independent(independent);
    CppADRecordingGuard<double> recording;
    auto nested = objective.fun.base2ad();
    if (objective.fun.size_dyn_ind()) {
      if (objective.dynamic_values.size() != objective.fun.size_dyn_ind()) {
        throw std::logic_error(
          "An objective tape's dynamic values are unavailable for nested AD.");
      }
      std::vector<AD> dynamic(
        objective.dynamic_values.begin(), objective.dynamic_values.end());
      nested.new_dynamic(dynamic);
    }
    std::ostringstream messages;
    nested.Forward(0, independent, messages);
    const std::vector<AD> weight(1, AD(1.0));
    const std::vector<AD> gradient = nested.Reverse(1, weight);
    std::vector<AD> dependent(eta_positions.size());
    for (std::size_t effect = 0; effect < eta_positions.size(); ++effect) {
      const int position = eta_positions[effect];
      if (position < 0 || position >= static_cast<int>(gradient.size())) {
        throw std::logic_error("A conditional-mode derivative position is invalid.");
      }
      dependent[effect] = gradient[static_cast<std::size_t>(position)];
    }
    CppAD::ADFun<double> result;
    result.Dependent(independent, dependent);
    recording.release();
    result.optimize();
    return result;
  }

  void exact_subject_marginal_curvature(
      ObjectiveTape& objective, ObjectiveTape* curvature,
      const std::vector<double>& point, const std::vector<int>& population,
      const std::vector<int>& eta_positions, Vector& native_gradient,
      Matrix& native_hessian) {
    const int n_native = static_cast<int>(population.size());
    const int n_eta = static_cast<int>(eta_positions.size());
    Vector joint_gradient = tape_gradient(
      objective, point, "exact marginal objective Hessian");
    Matrix joint_hessian = tape_hessian(
      objective, point, "exact marginal objective Hessian");
    Vector total_gradient = joint_gradient;
    Matrix total_hessian = joint_hessian;
    if (curvature) {
      total_gradient += tape_gradient(
        *curvature, point, "exact marginal curvature Hessian");
      total_hessian += tape_hessian(
        *curvature, point, "exact marginal curvature Hessian");
    }

    Matrix mode_hessian(n_eta, n_eta);
    Matrix mixed(n_eta, n_native);
    for (int row = 0; row < n_eta; ++row) {
      for (int column = 0; column < n_eta; ++column) {
        mode_hessian(row, column) = joint_hessian(
          eta_positions[static_cast<std::size_t>(row)],
          eta_positions[static_cast<std::size_t>(column)]);
      }
      for (int column = 0; column < n_native; ++column) {
        mixed(row, column) = joint_hessian(
          eta_positions[static_cast<std::size_t>(row)],
          population[static_cast<std::size_t>(column)]);
      }
    }
    if (n_eta) {
      auto eigen = libertad::detail::self_adjoint_eigen(mode_hessian, false);
      if (eigen.info != Eigen::Success) {
        throw std::runtime_error(
          "Conditional curvature decomposition failed for the exact Hessian.");
      }
      const double largest = std::max(1.0, eigen.values.cwiseAbs().maxCoeff());
      if (eigen.values.minCoeff() <= largest * 1e-10) {
        throw std::runtime_error(
          "The exact population Hessian requires positive-definite, unregularized "
          "conditional curvature.");
      }
    }
    Matrix sensitivity;
    if (n_eta) {
      sensitivity = -mode_hessian.ldlt().solve(mixed);
    } else {
      sensitivity = Matrix::Zero(0, n_native);
    }

    Vector population_gradient(n_native);
    Vector mode_gradient(n_eta);
    for (int index = 0; index < n_native; ++index) {
      population_gradient[index] = total_gradient[
        population[static_cast<std::size_t>(index)]];
    }
    for (int effect = 0; effect < n_eta; ++effect) {
      mode_gradient[effect] = total_gradient[
        eta_positions[static_cast<std::size_t>(effect)]];
    }
    native_gradient += population_gradient +
      sensitivity.transpose() * mode_gradient;

    Matrix pp(n_native, n_native);
    Matrix pe(n_native, n_eta);
    Matrix ep(n_eta, n_native);
    Matrix ee(n_eta, n_eta);
    for (int row = 0; row < n_native; ++row) {
      for (int column = 0; column < n_native; ++column) {
        pp(row, column) = total_hessian(
          population[static_cast<std::size_t>(row)],
          population[static_cast<std::size_t>(column)]);
      }
      for (int effect = 0; effect < n_eta; ++effect) {
        pe(row, effect) = total_hessian(
          population[static_cast<std::size_t>(row)],
          eta_positions[static_cast<std::size_t>(effect)]);
      }
    }
    for (int effect = 0; effect < n_eta; ++effect) {
      for (int column = 0; column < n_native; ++column) {
        ep(effect, column) = total_hessian(
          eta_positions[static_cast<std::size_t>(effect)],
          population[static_cast<std::size_t>(column)]);
      }
      for (int other = 0; other < n_eta; ++other) {
        ee(effect, other) = total_hessian(
          eta_positions[static_cast<std::size_t>(effect)],
          eta_positions[static_cast<std::size_t>(other)]);
      }
    }
    Matrix result = pp + pe * sensitivity + sensitivity.transpose() * ep +
      sensitivity.transpose() * ee * sensitivity;

    if (n_eta && mode_gradient.lpNorm<Eigen::Infinity>() > 1e-14) {
      CppAD::ADFun<double> gradient_tape = mode_gradient_tape(
        objective, point, eta_positions);
      const Vector multiplier = mode_hessian.ldlt().solve(mode_gradient);
      std::vector<double> weight(static_cast<std::size_t>(n_eta), 0.0);
      for (int equation = 0; equation < n_eta; ++equation) {
        if (multiplier[equation] == 0.0) continue;
        weight[static_cast<std::size_t>(equation)] = 1.0;
        const Matrix third = adfun_hessian(gradient_tape, point, weight);
        weight[static_cast<std::size_t>(equation)] = 0.0;
        Matrix contraction(n_native, n_native);
        for (int first = 0; first < n_native; ++first) {
          for (int second = 0; second < n_native; ++second) {
            double value = third(
              population[static_cast<std::size_t>(first)],
              population[static_cast<std::size_t>(second)]);
            for (int left = 0; left < n_eta; ++left) {
              value += third(
                population[static_cast<std::size_t>(first)],
                eta_positions[static_cast<std::size_t>(left)]) *
                  sensitivity(left, second);
              value += third(
                eta_positions[static_cast<std::size_t>(left)],
                population[static_cast<std::size_t>(second)]) *
                  sensitivity(left, first);
              for (int right = 0; right < n_eta; ++right) {
                value += sensitivity(left, first) * third(
                  eta_positions[static_cast<std::size_t>(left)],
                  eta_positions[static_cast<std::size_t>(right)]) *
                  sensitivity(right, second);
              }
            }
            contraction(first, second) = value;
          }
        }
        result -= multiplier[equation] * contraction;
      }
    }
    native_hessian += 0.5 * (result + result.transpose()).eval();
  }

  CppAD::ADFun<double> native_transform_tape(
      const std::vector<double>& encoded) const {
    using AD = CppAD::AD<double>;
    std::vector<AD> independent(encoded.begin(), encoded.end());
    CppAD::Independent(independent);
    CppADRecordingGuard<double> recording;
    std::vector<AD> theta(theta_base_.begin(), theta_base_.end());
    std::vector<AD> sigma(sigma_base_.begin(), sigma_base_.end());
    std::vector<AD> omega(omega_base_.begin(), omega_base_.end());
    std::size_t cursor = 0;
    for (int index : theta_free_) theta[static_cast<std::size_t>(index)] = independent[cursor++];
    for (int index : sigma_free_) sigma[static_cast<std::size_t>(index)] = CppAD::exp(independent[cursor++]);
    if (omega_full_ && !omega_free_.empty()) {
      MatrixT<AD> lower = MatrixT<AD>::Zero(n_eta_base_, n_eta_base_);
      for (std::size_t entry = 0; entry < omega.size(); ++entry) {
        const int row = omega_rows_[entry];
        const int column = omega_cols_[entry];
        lower(row, column) = row == column ?
          CppAD::exp(independent[cursor + entry]) : independent[cursor + entry];
      }
      const MatrixT<AD> covariance = lower * lower.transpose();
      for (std::size_t entry = 0; entry < omega.size(); ++entry) {
        omega[entry] = covariance(omega_rows_[entry], omega_cols_[entry]);
      }
      cursor += omega.size();
    } else {
      for (int index : omega_free_) {
        omega[static_cast<std::size_t>(index)] = CppAD::exp(independent[cursor++]);
      }
    }
    if (cursor != independent.size()) {
      throw std::logic_error("The population transform tape consumed the wrong input count.");
    }
    std::vector<AD> dependent;
    dependent.reserve(theta.size() + sigma.size() + omega.size());
    dependent.insert(dependent.end(), theta.begin(), theta.end());
    dependent.insert(dependent.end(), sigma.begin(), sigma.end());
    dependent.insert(dependent.end(), omega.begin(), omega.end());
    CppAD::ADFun<double> result;
    result.Dependent(independent, dependent);
    recording.release();
    result.optimize();
    return result;
  }

  std::vector<double> objective_point(
      const PopulationParameters& parameters, const Vector& eta) const {
    std::vector<double> result;
    result.reserve(parameters.theta.size() + static_cast<std::size_t>(eta.size()) +
                   parameters.sigma.size() + parameters.omega.size());
    result.insert(result.end(), parameters.theta.begin(), parameters.theta.end());
    for (Eigen::Index index = 0; index < eta.size(); ++index) result.push_back(eta[index]);
    result.insert(result.end(), parameters.sigma.begin(), parameters.sigma.end());
    result.insert(result.end(), parameters.omega.begin(), parameters.omega.end());
    return result;
  }

  std::vector<double> fo_point(const PopulationParameters& parameters) const {
    std::vector<double> result;
    result.reserve(parameters.theta.size() + parameters.sigma.size() +
                   parameters.omega.size());
    result.insert(result.end(), parameters.theta.begin(), parameters.theta.end());
    result.insert(result.end(), parameters.sigma.begin(), parameters.sigma.end());
    result.insert(result.end(), parameters.omega.begin(), parameters.omega.end());
    return result;
  }

  std::vector<double> fo_population_dynamic_values() {
    std::vector<double> result;
    for (std::size_t subject = 0; subject < primary_.size(); ++subject) {
      Rcpp::DataFrame data(subject_data_[subject]);
      const std::vector<double> values = fo_dynamic_values(*primary_[subject], data);
      result.insert(result.end(), values.begin(), values.end());
    }
    fo_dynamic_updates_ += static_cast<long long>(primary_.size());
    return result;
  }

  void record_fo_population(const PopulationParameters& parameters) {
    if (!is_fo() || primary_.empty()) return;
    double estimated_operations = 0.0;
    for (ObjectiveTape* tape : primary_) {
      estimated_operations += static_cast<double>(tape->fun.size_op());
    }
    if (estimated_operations > fo_population_max_operations_) {
      throw std::runtime_error(
        "The estimated fused FO tape exceeds LibeRation.fo_population_max_operations; "
        "using subject tapes.");
    }
    const std::vector<double> point = fo_point(parameters);
    std::vector<double> dynamic_values = fo_population_dynamic_values();
    std::vector<CppAD::AD<double>> independent(point.begin(), point.end());
    std::vector<CppAD::AD<double>> dynamic(
      dynamic_values.begin(), dynamic_values.end());
    if (dynamic.empty()) CppAD::Independent(independent);
    else CppAD::Independent(independent, dynamic);

    std::vector<CppAD::AD<double>> dependent(primary_.size());
    std::size_t cursor = 0;
    std::ostringstream messages;
    for (std::size_t subject = 0; subject < primary_.size(); ++subject) {
      ObjectiveTape& source = *primary_[subject];
      if (source.domain_names.size() != point.size()) {
        throw std::invalid_argument("FO subject tapes have inconsistent parameter domains.");
      }
      auto nested = source.fun.base2ad();
      const std::size_t count = source.fun.size_dyn_ind();
      if (cursor + count > dynamic.size()) {
        throw std::logic_error("FO population dynamic offsets are inconsistent.");
      }
      if (count) {
        std::vector<CppAD::AD<double>> current(
          dynamic.begin() + static_cast<std::ptrdiff_t>(cursor),
          dynamic.begin() + static_cast<std::ptrdiff_t>(cursor + count));
        nested.new_dynamic(current);
      }
      const std::vector<CppAD::AD<double>> value =
        nested.Forward(0, independent, messages);
      if (value.size() != 1U) {
        throw std::logic_error("An FO subject tape did not return one objective value.");
      }
      dependent[subject] = value[0];
      cursor += count;
    }
    if (cursor != dynamic.size()) {
      throw std::logic_error("FO population dynamic data were not fully consumed.");
    }
    auto population = std::make_unique<ObjectiveTape>();
    population->fun.Dependent(independent, dependent);
    population->fun.optimize();
    population->domain_names = primary_[0]->domain_names;
    population->dynamic_values = std::move(dynamic_values);
    fo_population_ = std::move(population);
    ++tape_records_;
  }

  void evaluate_fo_population(const PopulationParameters& parameters,
                              double prior) {
    if (!fo_population_) throw std::logic_error("FO population tape is unavailable.");
    const std::vector<double> dynamic_values = fo_population_dynamic_values();
    if (dynamic_values.size() != fo_population_->fun.size_dyn_ind()) {
      throw std::logic_error("FO population dynamic data have the wrong length.");
    }
    if (!dynamic_values.empty()) fo_population_->fun.new_dynamic(dynamic_values);
    fo_population_->dynamic_values = dynamic_values;
    const std::vector<double> point = fo_point(parameters);
    std::ostringstream messages;
    const std::vector<double> values = fo_population_->fun.Forward(0, point, messages);
    require_unchanged_path(fo_population_->fun, "batched FO population objective");
    if (values.size() != primary_.size()) {
      throw std::logic_error("Batched FO population output has the wrong length.");
    }
    cache_points_.assign(primary_.size(), point);
    cache_subject_values_ = values;
    double total = prior;
    for (double value : values) {
      if (!std::isfinite(value)) {
        total = penalty();
        break;
      }
      total += value;
    }
    cache_value_ = std::isfinite(total) ? total : penalty();
    cache_valid_ = true;
  }

  std::vector<double> anchor_point(
      const PopulationParameters& parameters, const Vector& eta) const {
    return objective_point(parameters, eta);
  }

  bool material_movement(int subject, const PopulationParameters& parameters,
                         const Vector& eta) const {
    if (!use_ode_) return false;
    const std::vector<double> point = anchor_point(parameters, eta);
    const std::vector<double>& anchor = anchors_[static_cast<std::size_t>(subject)];
    if (point.size() != anchor.size()) return true;
    double distance = 0.0;
    for (std::size_t index = 0; index < point.size(); ++index) {
      distance = std::max(distance, std::abs(point[index] - anchor[index]) /
        std::max(std::abs(anchor[index]), 1.0));
    }
    return std::isfinite(distance) && distance > guard_radius_;
  }

  void record_subject(int subject, const PopulationParameters& parameters,
                      const Vector& eta, bool retape) {
    Rcpp::DataFrame data(subject_data_[static_cast<std::size_t>(subject)]);
    Rcpp::NumericVector theta = Rcpp::wrap(parameters.theta);
    Rcpp::NumericVector sigma = Rcpp::wrap(parameters.sigma);
    Rcpp::NumericVector omega = Rcpp::wrap(parameters.omega);
    Rcpp::NumericMatrix eta_matrix(1, n_eta_);
    Rcpp::NumericVector eta_vector(n_eta_);
    for (int effect = 0; effect < n_eta_; ++effect) {
      eta_matrix(0, effect) = eta[effect];
      eta_vector[effect] = eta[effect];
    }
    const std::size_t index = static_cast<std::size_t>(subject);
    if (is_fo()) {
      owned_prediction_[index] = record_prediction_tape(
        *engine_, data, theta, eta_matrix, sigma);
      owned_primary_[index] = record_fo_tape(
        *engine_, *owned_prediction_[index], data, theta, sigma, omega);
    } else {
      owned_primary_[index] = record_objective_tape(
        *engine_, data, theta, eta_matrix, sigma, omega, interaction());
      if (has_curvature()) {
        owned_prediction_[index] = record_prediction_tape(
          *engine_, data, theta, eta_matrix, sigma);
        owned_curvature_[index] = record_curvature_tape(
          *engine_, *owned_prediction_[index], *owned_primary_[index], data,
          theta, eta_vector, sigma, omega, approximation_);
      }
    }
    primary_[index] = owned_primary_[index].get();
    curvature_[index] = has_curvature() ? owned_curvature_[index].get() : nullptr;
    anchors_[index] = anchor_point(parameters, eta);
    ++tape_records_;
    if (retape) ++tape_retapes_;
  }

  void record_primary_subject(int subject,
                              const PopulationParameters& parameters,
                              const Vector& eta, bool retape) {
    if (is_fo()) {
      record_subject(subject, parameters, eta, retape);
      return;
    }
    Rcpp::DataFrame data(subject_data_[static_cast<std::size_t>(subject)]);
    Rcpp::NumericVector theta = Rcpp::wrap(parameters.theta);
    Rcpp::NumericVector sigma = Rcpp::wrap(parameters.sigma);
    Rcpp::NumericVector omega = Rcpp::wrap(parameters.omega);
    Rcpp::NumericMatrix eta_matrix(1, n_eta_);
    for (int effect = 0; effect < n_eta_; ++effect) {
      eta_matrix(0, effect) = eta[effect];
    }
    const std::size_t index = static_cast<std::size_t>(subject);
    owned_primary_[index] = record_objective_tape(
      *engine_, data, theta, eta_matrix, sigma, omega, interaction());
    primary_[index] = owned_primary_[index].get();
    anchors_[index] = anchor_point(parameters, eta);
    ++tape_records_;
    if (retape) ++tape_retapes_;
  }

  bool ensure_tape(int subject, const PopulationParameters& parameters,
                   const Vector& eta) {
    if (!material_movement(subject, parameters, eta)) return false;
    record_subject(subject, parameters, eta, true);
    return true;
  }

  static double tape_value(ObjectiveTape& tape, const std::vector<double>& point) {
    if (point.size() != tape.domain_names.size()) {
      throw std::invalid_argument("A compiled population point has the wrong length.");
    }
    std::ostringstream messages;
    const std::vector<double> value = tape.fun.Forward(0, point, messages);
    if (value.empty() || !std::isfinite(value[0])) {
      return std::numeric_limits<double>::infinity();
    }
    require_unchanged_path(tape.fun, "compiled population objective");
    return value[0];
  }

  Rcpp::List mode_at_once(int subject, const PopulationParameters& parameters,
                          const Vector& start) {
    const std::size_t index = static_cast<std::size_t>(subject);
    std::vector<double> point = objective_point(parameters, start);
    std::vector<std::size_t> positions(static_cast<std::size_t>(n_eta_));
    for (int effect = 0; effect < n_eta_; ++effect) {
      positions[static_cast<std::size_t>(effect)] = parameters.theta.size() +
        static_cast<std::size_t>(effect);
    }
    Rcpp::NumericVector start_vector(n_eta_);
    for (int effect = 0; effect < n_eta_; ++effect) start_vector[effect] = start[effect];
    Rcpp::List result = objective_eta_mode(
      *primary_[index], std::move(point), positions, start_vector,
      eta_maxit_, tolerance_, false);
    int total_iterations = Rcpp::as<int>(result["iterations"]);
    int total_evaluations = Rcpp::as<int>(result["evaluations"]);
    auto accept_relative_mode = [this, &parameters, &positions, index](
        Rcpp::List& candidate) {
      if (Rcpp::as<int>(candidate["convergence"]) == 0) return;
      const Rcpp::NumericVector gradient = candidate["gradient"];
      const Rcpp::NumericVector eta = candidate["par"];
      const double value = Rcpp::as<double>(candidate["value"]);
      double norm = 0.0;
      for (double current : gradient) norm = std::max(norm, std::abs(current));
      const double relative_tolerance = tolerance_ * (1.0 + std::abs(value));
      if (std::isfinite(value) && std::isfinite(norm) &&
          norm <= std::max(10.0 * tolerance_, relative_tolerance)) {
        candidate["convergence"] = 0;
        return;
      }
      // A very small Newton displacement is a scale-aware convergence test
      // when residual variance makes the conditional gradient and Hessian
      // simultaneously large.
      Vector eta_eigen(n_eta_), gradient_eigen(n_eta_);
      for (int effect = 0; effect < n_eta_; ++effect) {
        eta_eigen[effect] = eta[effect];
        gradient_eigen[effect] = gradient[effect];
      }
      std::vector<double> current_point = objective_point(parameters, eta_eigen);
      Matrix hessian = objective_eta_hessian(
        *primary_[index], current_point, positions);
      auto eigen = libertad::detail::self_adjoint_eigen(hessian, false);
      if (eigen.info != Eigen::Success) return;
      const double largest = std::max(eigen.values.cwiseAbs().maxCoeff(), 1.0);
      const double jitter = std::max(0.0, largest * 1e-12 -
        eigen.values.minCoeff());
      hessian.diagonal().array() += jitter;
      const Vector displacement = hessian.ldlt().solve(gradient_eigen);
      if (displacement.allFinite() &&
          displacement.lpNorm<Eigen::Infinity>() <= std::sqrt(tolerance_) *
            (1.0 + eta_eigen.lpNorm<Eigen::Infinity>())) {
        candidate["convergence"] = 0;
      }
    };
    accept_relative_mode(result);
    // A fresh inverse-Hessian approximation is often enough to recover from
    // a conditional line-search failure at an extreme outer trial point. The
    // R implementation historically obtained the same robustness by falling
    // back to a second BFGS invocation; keep that recovery inside C++.
    for (int restart = 0;
         restart < 2 && Rcpp::as<int>(result["convergence"]) != 0;
         ++restart) {
      Rcpp::NumericVector restart_at = result["par"];
      std::vector<double> restart_point = objective_point(parameters, Vector::Zero(n_eta_));
      for (int effect = 0; effect < n_eta_; ++effect) {
        restart_point[parameters.theta.size() + static_cast<std::size_t>(effect)] =
          restart_at[effect];
      }
      result = objective_eta_mode(
        *primary_[index], std::move(restart_point), positions, restart_at,
        eta_maxit_, tolerance_, false);
      total_iterations += Rcpp::as<int>(result["iterations"]);
      total_evaluations += Rcpp::as<int>(result["evaluations"]);
      accept_relative_mode(result);
    }
    result["iterations"] = total_iterations;
    result["evaluations"] = total_evaluations;
    return result;
  }

  Rcpp::List mode_at(int subject, const PopulationParameters& parameters,
                     const Vector& start) {
    Vector anchor = start;
    for (int attempt = 0; attempt < 12; ++attempt) {
      try {
        return mode_at_once(subject, parameters, anchor);
      } catch (const TapePathChange& change) {
        Vector candidate = anchor;
        if (change.point().size() >= parameters.theta.size() +
                                   static_cast<std::size_t>(n_eta_)) {
          for (int effect = 0; effect < n_eta_; ++effect) {
            candidate[effect] = change.point()[parameters.theta.size() +
                                               static_cast<std::size_t>(effect)];
          }
        }
        // A Laplace curvature tape may not be positive definite away from the
        // conditional mode. Retape only the primary objective during mode
        // search; the curvature tape is refreshed once the mode is known.
        // A path-changing line-search trial can also sit outside the valid
        // pharmacological domain (for example exp(ETA) underflowing a rate to
        // zero). Backtrack toward the last valid anchor until direct recording
        // succeeds rather than promoting that invalid trial to a tape anchor.
        // Bound a single retape displacement as well: a BFGS line search can
        // briefly propose ETAs hundreds of units away, which may still pass a
        // simple positivity test but is not a numerically meaningful anchor.
        Vector displacement = candidate - anchor;
        const double maximum = displacement.allFinite() ?
          displacement.lpNorm<Eigen::Infinity>() :
          std::numeric_limits<double>::infinity();
        constexpr double maximum_retape_eta_step = 2.0;
        if (maximum > maximum_retape_eta_step) {
          candidate = anchor + displacement * (maximum_retape_eta_step / maximum);
        }
        bool recorded = false;
        for (int backtrack = 0; backtrack < 24 && !recorded; ++backtrack) {
          if (!candidate.allFinite()) candidate = anchor;
          try {
            record_primary_subject(subject, parameters, candidate, true);
            anchor = candidate;
            recorded = true;
          } catch (const std::domain_error&) {
            candidate = 0.5 * (anchor + candidate);
          }
        }
        if (!recorded) throw;
      }
    }
    return mode_at_once(subject, parameters, anchor);
  }

  double guarded_tape_value(int subject, const PopulationParameters& parameters,
                            const Vector& eta, bool curvature,
                            const std::vector<double>& point) {
    for (int attempt = 0; attempt < 3; ++attempt) {
      try {
        ObjectiveTape* tape = curvature ?
          curvature_[static_cast<std::size_t>(subject)] :
          primary_[static_cast<std::size_t>(subject)];
        return tape_value(*tape, point);
      } catch (const TapePathChange&) {
        if (curvature) record_subject(subject, parameters, eta, true);
        else record_primary_subject(subject, parameters, eta, true);
      }
    }
    ObjectiveTape* tape = curvature ?
      curvature_[static_cast<std::size_t>(subject)] :
      primary_[static_cast<std::size_t>(subject)];
    return tape_value(*tape, point);
  }

  void evaluate_value(const std::vector<double>& encoded) {
    cache_valid_ = false;
    cache_gradient_valid_ = false;
    cache_key_ = encoded;
    cache_parameters_ = decode(encoded);
    ++parameter_evaluations_;
    const int subjects = static_cast<int>(primary_.size());
    cache_points_.assign(static_cast<std::size_t>(subjects), std::vector<double>());
    cache_subject_values_.assign(static_cast<std::size_t>(subjects),
                                 std::numeric_limits<double>::infinity());
    cache_curvature_values_.assign(static_cast<std::size_t>(subjects), 0.0);
    cache_mode_convergence_.assign(static_cast<std::size_t>(subjects), 0);
    double total = prior_nll(cache_parameters_);
    if (!std::isfinite(total) || total >= penalty()) {
      cache_value_ = penalty();
      cache_valid_ = true;
      return;
    }
    if (is_fo() && fo_population_) {
      try {
        evaluate_fo_population(cache_parameters_, total);
        return;
      } catch (const TapePathChange&) {
        fo_population_.reset();
        ++fo_population_fallbacks_;
        fo_population_error_ = "A batched tape path changed; using subject tapes.";
      }
    }
    if (is_fo()) {
      const Vector zero_eta = Vector::Zero(n_eta_);
      const std::vector<double> point = fo_point(cache_parameters_);
      for (int subject = 0; subject < subjects; ++subject) {
        if (use_ode_) ensure_tape(subject, cache_parameters_, zero_eta);
        Rcpp::DataFrame data(subject_data_[static_cast<std::size_t>(subject)]);
        set_fo_dynamic(*primary_[static_cast<std::size_t>(subject)], data);
        ++fo_dynamic_updates_;
        const double current = guarded_tape_value(
          subject, cache_parameters_, zero_eta, false, point);
        if (!std::isfinite(current)) {
          total = penalty();
          break;
        }
        cache_points_[static_cast<std::size_t>(subject)] = point;
        cache_subject_values_[static_cast<std::size_t>(subject)] = current;
        total += current;
      }
    } else {
      for (int subject = 0; subject < subjects; ++subject) {
        Vector start = starts_.row(subject).transpose();
        if (use_ode_) ensure_tape(subject, cache_parameters_, start);
        Vector eta = start;
        double current = std::numeric_limits<double>::infinity();
        int convergence = 0;
        if (n_eta_) {
          Rcpp::List mode = mode_at(subject, cache_parameters_, start);
          convergence = Rcpp::as<int>(mode["convergence"]);
          mode_iterations_ += Rcpp::as<int>(mode["iterations"]);
          mode_evaluations_ += Rcpp::as<int>(mode["evaluations"]);
          Rcpp::NumericVector par = mode["par"];
          for (int effect = 0; effect < n_eta_; ++effect) eta[effect] = par[effect];
          current = Rcpp::as<double>(mode["value"]);
          if (convergence == 0 && use_ode_ &&
              ensure_tape(subject, cache_parameters_, eta)) {
            mode = mode_at(subject, cache_parameters_, eta);
            convergence = Rcpp::as<int>(mode["convergence"]);
            mode_iterations_ += Rcpp::as<int>(mode["iterations"]);
            mode_evaluations_ += Rcpp::as<int>(mode["evaluations"]);
            par = Rcpp::as<Rcpp::NumericVector>(mode["par"]);
            for (int effect = 0; effect < n_eta_; ++effect) eta[effect] = par[effect];
            current = Rcpp::as<double>(mode["value"]);
          }
        } else {
          const std::vector<double> point = objective_point(cache_parameters_, eta);
          current = guarded_tape_value(
            subject, cache_parameters_, eta, false, point);
          ++mode_evaluations_;
        }
        if (!std::isfinite(current)) {
          total = penalty();
          break;
        }
        if (convergence != 0) {
          // Keep a finite approximate mode at an extreme outer line-search
          // point so L-BFGS-B can obtain a gradient and reject that point.
          // Converged points remain exact; the recovery is reported and the
          // final cached state retains its convergence code.
          ++mode_recoveries_;
          cache_mode_convergence_[static_cast<std::size_t>(subject)] = convergence;
        }
        starts_.row(subject) = eta.transpose();
        const std::vector<double> point = objective_point(cache_parameters_, eta);
        cache_points_[static_cast<std::size_t>(subject)] = point;
        cache_subject_values_[static_cast<std::size_t>(subject)] = current;
        total += current;
        if (has_curvature()) {
          const double determinant = guarded_tape_value(
            subject, cache_parameters_, eta, true, point);
          if (!std::isfinite(determinant)) {
            total = penalty();
            break;
          }
          cache_curvature_values_[static_cast<std::size_t>(subject)] = determinant;
          total += determinant;
        }
      }
    }
    cache_value_ = std::isfinite(total) ? total : penalty();
    cache_valid_ = true;
  }

  void evaluate_gradient() {
    const int n_theta = static_cast<int>(cache_parameters_.theta.size());
    const int n_sigma = static_cast<int>(cache_parameters_.sigma.size());
    const int n_omega = static_cast<int>(cache_parameters_.omega.size());
    const int n_native = n_theta + n_sigma + n_omega;
    const int n_outer = static_cast<int>(cache_key_.size());
    Vector outer = Vector::Zero(n_outer);
    Vector native_prior;
    prior_nll(cache_parameters_, &native_prior);
    const std::vector<double> weight(1, 1.0);
    std::vector<int> population;
    population.reserve(static_cast<std::size_t>(n_native));
    for (int index = 0; index < n_theta; ++index) population.push_back(index);
    for (int index = 0; index < n_sigma; ++index) {
      population.push_back(n_theta + n_eta_ + index);
    }
    for (int index = 0; index < n_omega; ++index) {
      population.push_back(n_theta + n_eta_ + n_sigma + index);
    }

    if (is_fo() && fo_population_) {
      // The population tape has one output per subject.  Form its Jacobian and
      // add rows in subject order so the fused route retains the deterministic
      // summation order of the established subject-tape implementation.  With
      // far fewer population parameters than subjects CppAD uses a small
      // number of forward sweeps here, rather than one reverse sweep per
      // subject.
      const std::vector<double> point = fo_point(cache_parameters_);
      const std::vector<double> derivative = fo_population_->fun.Jacobian(point);
      require_unchanged_path(fo_population_->fun, "batched FO population gradient");
      Vector native = Vector::Zero(n_native);
      if (derivative.size() != primary_.size() * static_cast<std::size_t>(n_native)) {
        throw std::logic_error("Batched FO population gradient has the wrong length.");
      }
      for (std::size_t subject = 0; subject < primary_.size(); ++subject) {
        for (int index = 0; index < n_native; ++index) {
          native[index] += derivative[
            subject * static_cast<std::size_t>(n_native) +
            static_cast<std::size_t>(index)];
        }
      }
      native += native_prior;
      outer = cache_parameters_.transform.transpose() * native;
    } else if (is_fo() || approximation_ == "its") {
      Vector native = Vector::Zero(n_native);
      for (std::size_t subject = 0; subject < primary_.size(); ++subject) {
        ObjectiveTape& objective = *primary_[subject];
        if (is_fo()) {
          Rcpp::DataFrame data(subject_data_[subject]);
          set_fo_dynamic(objective, data);
          ++fo_dynamic_updates_;
        }
        std::ostringstream messages;
        objective.fun.Forward(0, cache_points_[subject], messages);
        if (is_fo()) require_unchanged_path(objective.fun, "shared FO objective gradient");
        const std::vector<double> derivative = objective.fun.Reverse(1, weight);
        if (is_fo()) {
          for (int index = 0; index < n_native; ++index) native[index] += derivative[index];
        } else {
          for (int index = 0; index < n_native; ++index) {
            native[index] += derivative[static_cast<std::size_t>(population[index])];
          }
        }
      }
      native += native_prior;
      outer = cache_parameters_.transform.transpose() * native;
    } else {
      std::vector<std::size_t> eta_positions(static_cast<std::size_t>(n_eta_));
      for (int effect = 0; effect < n_eta_; ++effect) {
        eta_positions[static_cast<std::size_t>(effect)] =
          static_cast<std::size_t>(n_theta + effect);
      }
      for (std::size_t subject = 0; subject < primary_.size(); ++subject) {
        ObjectiveTape& objective = *primary_[subject];
        ObjectiveTape& curvature = *curvature_[subject];
        const std::vector<double>& point = cache_points_[subject];
        std::ostringstream messages;
        objective.fun.Forward(0, point, messages);
        const std::vector<double> objective_derivative = objective.fun.Reverse(1, weight);
        Matrix mixed(n_eta_, n_eta_ + n_native);
        std::vector<double> direction(objective.domain_names.size(), 0.0);
        for (int column = 0; column < n_eta_ + n_native; ++column) {
          const std::size_t position = column < n_eta_ ?
            eta_positions[static_cast<std::size_t>(column)] :
            static_cast<std::size_t>(population[static_cast<std::size_t>(column - n_eta_)]);
          direction[position] = 1.0;
          objective.fun.Forward(1, direction, messages);
          direction[position] = 0.0;
          const std::vector<double> reverse = objective.fun.Reverse(2, weight);
          for (int row = 0; row < n_eta_; ++row) {
            mixed(row, column) = reverse[
              eta_positions[static_cast<std::size_t>(row)] * 2U + 1U];
          }
        }
        Matrix eta_hessian;
        if (n_eta_) {
          eta_hessian = 0.5 *
            (mixed.leftCols(n_eta_) + mixed.leftCols(n_eta_).transpose()).eval();
        } else {
          eta_hessian = Matrix::Zero(0, 0);
        }
        if (n_eta_) {
          auto eigen = libertad::detail::self_adjoint_eigen(eta_hessian, false);
          if (eigen.info != Eigen::Success) {
            throw std::runtime_error("Conditional ETA curvature decomposition failed.");
          }
          const double largest = std::max(eigen.values.cwiseAbs().maxCoeff(), 1.0);
          const double jitter = std::max(0.0, largest * 1e-9 -
            eigen.values.minCoeff());
          if (jitter > largest * 1e-2) {
            throw std::runtime_error("Conditional ETA curvature is not positive definite.");
          }
          eta_hessian.diagonal().array() += jitter;
        }
        Matrix sensitivity;
        if (n_eta_) {
          sensitivity = -eta_hessian.ldlt().solve(
            mixed.rightCols(n_native) * cache_parameters_.transform);
        } else {
          sensitivity = Matrix::Zero(0, n_outer);
        }
        curvature.fun.Forward(0, point, messages);
        const std::vector<double> curvature_derivative = curvature.fun.Reverse(1, weight);
        for (int encoded = 0; encoded < n_outer; ++encoded) {
          double current = 0.0;
          for (int native = 0; native < n_native; ++native) {
            current += (objective_derivative[static_cast<std::size_t>(population[native])] +
              curvature_derivative[static_cast<std::size_t>(population[native])]) *
              cache_parameters_.transform(native, encoded);
          }
          for (int effect = 0; effect < n_eta_; ++effect) {
            current += curvature_derivative[eta_positions[static_cast<std::size_t>(effect)]] *
              sensitivity(effect, encoded);
          }
          outer[encoded] += current;
        }
      }
      outer += cache_parameters_.transform.transpose() * native_prior;
    }
    cache_gradient_.resize(static_cast<std::size_t>(n_outer));
    for (int index = 0; index < n_outer; ++index) {
      cache_gradient_[static_cast<std::size_t>(index)] = outer[index];
      if (!std::isfinite(outer[index])) {
        throw std::runtime_error("Compiled population gradient is not finite.");
      }
    }
    cache_gradient_valid_ = true;
  }
};
