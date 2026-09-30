// Internal implementation fragment included inside namespace liberation by
// pk_engine.cpp. It intentionally has no include guard or namespace declaration.
// Contains the fixed-ETA SAEM objective, projected gradient, limited-memory optimiser, and line search.

class SaemFixedEtaObjective {
 public:
  SaemFixedEtaObjective(
      const Rcpp::List& tape_pointers, const Rcpp::NumericMatrix& eta,
      const Rcpp::NumericVector& theta, const Rcpp::NumericVector& sigma,
      const Rcpp::NumericVector& omega,
      const Rcpp::IntegerVector& theta_free,
      const Rcpp::IntegerVector& sigma_free,
      const Rcpp::List& prior_config)
      : theta_base_(Rcpp::as<std::vector<double>>(theta)),
        sigma_base_(Rcpp::as<std::vector<double>>(sigma)),
        omega_(Rcpp::as<std::vector<double>>(omega)),
        eta_(eta.nrow(), eta.ncol()), n_eta_(eta.ncol()) {
    if (tape_pointers.size() != eta.nrow() || tape_pointers.size() < 1) {
      throw std::invalid_argument(
        "SAEM fixed-ETA tapes and ETA rows must have equal non-zero length.");
    }
    for (int row = 0; row < eta.nrow(); ++row) {
      for (int column = 0; column < eta.ncol(); ++column) {
        if (!std::isfinite(eta(row, column))) {
          throw std::invalid_argument("SAEM fixed ETAs must be finite.");
        }
        eta_(row, column) = eta(row, column);
      }
    }
    tapes_.reserve(static_cast<std::size_t>(tape_pointers.size()));
    const std::size_t expected_domain = theta_base_.size() +
      static_cast<std::size_t>(eta.ncol()) + sigma_base_.size() + omega_.size();
    for (int subject = 0; subject < tape_pointers.size(); ++subject) {
      Rcpp::XPtr<ObjectiveTape> tape(tape_pointers[subject]);
      if (tape->domain_names.size() != expected_domain) {
        throw std::invalid_argument(
          "An SAEM objective tape has an inconsistent domain length.");
      }
      tapes_.push_back(tape.get());
    }
    theta_free_ = zero_based(Rcpp::as<std::vector<int>>(theta_free));
    sigma_free_ = zero_based(Rcpp::as<std::vector<int>>(sigma_free));
    for (int index : theta_free_) {
      if (index < 0 || index >= static_cast<int>(theta_base_.size())) {
        throw std::invalid_argument("An SAEM free THETA index is invalid.");
      }
    }
    for (int index : sigma_free_) {
      if (index < 0 || index >= static_cast<int>(sigma_base_.size())) {
        throw std::invalid_argument("An SAEM free SIGMA index is invalid.");
      }
    }
    parse_priors(prior_config);
  }

  SaemFixedEtaObjective(
      WeightedEtaCollection& weighted,
      const Rcpp::NumericVector& theta, const Rcpp::NumericVector& sigma,
      const Rcpp::NumericVector& omega,
      const Rcpp::IntegerVector& theta_free,
      const Rcpp::IntegerVector& sigma_free,
      const Rcpp::List& prior_config)
      : theta_base_(Rcpp::as<std::vector<double>>(theta)),
        sigma_base_(Rcpp::as<std::vector<double>>(sigma)),
        omega_(Rcpp::as<std::vector<double>>(omega)), eta_(0, 0),
        n_eta_(weighted.eta_dimension()), weighted_(&weighted) {
    theta_free_ = zero_based(Rcpp::as<std::vector<int>>(theta_free));
    sigma_free_ = zero_based(Rcpp::as<std::vector<int>>(sigma_free));
    for (int index : theta_free_) {
      if (index < 0 || index >= static_cast<int>(theta_base_.size())) {
        throw std::invalid_argument("An SAEM free THETA index is invalid.");
      }
    }
    for (int index : sigma_free_) {
      if (index < 0 || index >= static_cast<int>(sigma_base_.size())) {
        throw std::invalid_argument("An SAEM free SIGMA index is invalid.");
      }
    }
    parse_priors(prior_config);
  }

  std::size_t dimension() const {
    return theta_free_.size() + sigma_free_.size();
  }

  std::vector<double> start() const {
    std::vector<double> result;
    result.reserve(dimension());
    for (int index : theta_free_) {
      result.push_back(theta_base_[static_cast<std::size_t>(index)]);
    }
    for (int index : sigma_free_) {
      const double value = sigma_base_[static_cast<std::size_t>(index)];
      if (!(value > 0.0) || !std::isfinite(value)) {
        throw std::invalid_argument("A free SAEM SIGMA start is not positive.");
      }
      result.push_back(std::log(value));
    }
    return result;
  }

  void decode(const Vector& encoded, std::vector<double>& theta,
              std::vector<double>& sigma) const {
    if (encoded.size() != static_cast<Eigen::Index>(dimension())) {
      throw std::invalid_argument("The SAEM M-step parameter length is invalid.");
    }
    theta = theta_base_;
    sigma = sigma_base_;
    Eigen::Index cursor = 0;
    for (int index : theta_free_) {
      theta[static_cast<std::size_t>(index)] = encoded[cursor++];
    }
    for (int index : sigma_free_) {
      const double value = std::exp(encoded[cursor++]);
      if (!(value > 0.0) || !std::isfinite(value)) {
        throw std::domain_error("An encoded SAEM SIGMA is not finite and positive.");
      }
      sigma[static_cast<std::size_t>(index)] = value;
    }
  }

  SaemEvaluation evaluate(const Vector& encoded) {
    if (cache_valid_ && encoded.size() == cache_point_.size() &&
        encoded.isApprox(cache_point_, 0.0)) {
      return cache_;
    }
    std::vector<double> theta, sigma;
    decode(encoded, theta, sigma);
    const Eigen::Index native_size = static_cast<Eigen::Index>(
      theta.size() + sigma.size() + omega_.size());
    Vector native_gradient = Vector::Zero(native_size);
    const Eigen::Index full_native_size = static_cast<Eigen::Index>(
      theta.size() + static_cast<std::size_t>(n_eta_) + sigma.size() +
      omega_.size());
    Vector full_native_gradient = Vector::Zero(full_native_size);
    double value = 0.0;
    const std::vector<double> reverse_weight(1U, 1.0);
    std::ostringstream messages;
    if (weighted_) {
      const SaemEvaluation weighted = weighted_->evaluate_native(
        theta, sigma, omega_);
      value = weighted.value;
      if (weighted.gradient.size() != static_cast<Eigen::Index>(
          theta.size() + static_cast<std::size_t>(n_eta_) +
          sigma.size() + omega_.size())) {
        throw std::runtime_error(
          "A weighted SAEM objective returned an invalid gradient length.");
      }
      full_native_gradient = weighted.gradient;
      for (std::size_t index = 0; index < theta.size(); ++index) {
        native_gradient[static_cast<Eigen::Index>(index)] =
          weighted.gradient[static_cast<Eigen::Index>(index)];
      }
      const Eigen::Index source_sigma = static_cast<Eigen::Index>(
        theta.size() + static_cast<std::size_t>(n_eta_));
      const Eigen::Index target_sigma = static_cast<Eigen::Index>(theta.size());
      for (std::size_t index = 0; index < sigma.size(); ++index) {
        native_gradient[target_sigma + static_cast<Eigen::Index>(index)] =
          weighted.gradient[source_sigma + static_cast<Eigen::Index>(index)];
      }
    } else for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
      ObjectiveTape& tape = *tapes_[subject];
      std::vector<double> point;
      point.reserve(tape.domain_names.size());
      point.insert(point.end(), theta.begin(), theta.end());
      for (Eigen::Index effect = 0; effect < eta_.cols(); ++effect) {
        point.push_back(eta_(static_cast<Eigen::Index>(subject), effect));
      }
      point.insert(point.end(), sigma.begin(), sigma.end());
      point.insert(point.end(), omega_.begin(), omega_.end());
      const std::vector<double> forward = tape.fun.Forward(0, point, messages);
      require_unchanged_path(tape.fun, "native SAEM M-step objective");
      if (forward.size() != 1U || !std::isfinite(forward[0])) {
        return penalty_evaluation(encoded.size());
      }
      value += forward[0];
      const std::vector<double> derivative = tape.fun.Reverse(1, reverse_weight);
      require_unchanged_path(tape.fun, "native SAEM M-step gradient");
      if (derivative.size() != point.size()) {
        throw std::runtime_error("An SAEM tape returned an invalid gradient length.");
      }
      for (std::size_t index = 0; index < derivative.size(); ++index) {
        full_native_gradient[static_cast<Eigen::Index>(index)] +=
          derivative[index];
      }
      for (std::size_t index = 0; index < theta.size(); ++index) {
        native_gradient[static_cast<Eigen::Index>(index)] += derivative[index];
      }
      const std::size_t sigma_domain_offset = theta.size() +
        static_cast<std::size_t>(n_eta_);
      const Eigen::Index sigma_native_offset =
        static_cast<Eigen::Index>(theta.size());
      for (std::size_t index = 0; index < sigma.size(); ++index) {
        native_gradient[sigma_native_offset + static_cast<Eigen::Index>(index)] +=
          derivative[sigma_domain_offset + index];
      }
      if ((subject + 1U) % 64U == 0U) Rcpp::checkUserInterrupt();
    }
    Vector prior_gradient = Vector::Zero(native_size);
    const double prior = prior_nll(theta, sigma, prior_gradient);
    if (!std::isfinite(prior) || prior >= penalty()) {
      return penalty_evaluation(encoded.size());
    }
    native_gradient += prior_gradient;
    for (std::size_t index = 0; index < theta.size(); ++index) {
      full_native_gradient[static_cast<Eigen::Index>(index)] +=
        prior_gradient[static_cast<Eigen::Index>(index)];
    }
    const Eigen::Index full_sigma_offset = static_cast<Eigen::Index>(
      theta.size() + static_cast<std::size_t>(n_eta_));
    const Eigen::Index compact_sigma_offset = static_cast<Eigen::Index>(
      theta.size());
    for (std::size_t index = 0; index < sigma.size(); ++index) {
      full_native_gradient[full_sigma_offset + static_cast<Eigen::Index>(index)] +=
        prior_gradient[compact_sigma_offset + static_cast<Eigen::Index>(index)];
    }
    const Eigen::Index full_omega_offset = full_sigma_offset +
      static_cast<Eigen::Index>(sigma.size());
    const Eigen::Index compact_omega_offset = compact_sigma_offset +
      static_cast<Eigen::Index>(sigma.size());
    for (std::size_t index = 0; index < omega_.size(); ++index) {
      full_native_gradient[full_omega_offset + static_cast<Eigen::Index>(index)] +=
        prior_gradient[compact_omega_offset + static_cast<Eigen::Index>(index)];
    }
    value += prior;
    Vector gradient(static_cast<Eigen::Index>(dimension()));
    Eigen::Index cursor = 0;
    for (int index : theta_free_) gradient[cursor++] = native_gradient[index];
    const Eigen::Index sigma_native_offset =
      static_cast<Eigen::Index>(theta.size());
    for (int index : sigma_free_) {
      gradient[cursor++] = native_gradient[sigma_native_offset + index] *
        sigma[static_cast<std::size_t>(index)];
    }
    if (!std::isfinite(value) || !gradient.allFinite()) {
      return penalty_evaluation(encoded.size());
    }
    cache_point_ = encoded;
    cache_.value = value;
    cache_.gradient = gradient;
    cache_.native_gradient = full_native_gradient;
    cache_valid_ = true;
    return cache_;
  }

  double value(const Vector& encoded) {
    if (cache_valid_ && encoded.size() == cache_point_.size() &&
        encoded.isApprox(cache_point_, 0.0)) {
      return cache_.value;
    }
    std::vector<double> theta, sigma;
    decode(encoded, theta, sigma);
    double result = 0.0;
    if (weighted_) {
      result = weighted_->value_native(theta, sigma, omega_);
      Vector derivative = Vector::Zero(static_cast<Eigen::Index>(
        theta.size() + sigma.size() + omega_.size()));
      const double prior = prior_nll(theta, sigma, derivative);
      if (!std::isfinite(prior) || prior >= penalty()) return penalty();
      result += prior;
      return std::isfinite(result) ? result : penalty();
    }
    std::ostringstream messages;
    for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
      ObjectiveTape& tape = *tapes_[subject];
      std::vector<double> point;
      point.reserve(tape.domain_names.size());
      point.insert(point.end(), theta.begin(), theta.end());
      for (Eigen::Index effect = 0; effect < eta_.cols(); ++effect) {
        point.push_back(eta_(static_cast<Eigen::Index>(subject), effect));
      }
      point.insert(point.end(), sigma.begin(), sigma.end());
      point.insert(point.end(), omega_.begin(), omega_.end());
      const std::vector<double> forward = tape.fun.Forward(0, point, messages);
      require_unchanged_path(tape.fun, "native SAEM M-step value");
      if (forward.size() != 1U || !std::isfinite(forward[0])) return penalty();
      result += forward[0];
      if ((subject + 1U) % 256U == 0U) Rcpp::checkUserInterrupt();
    }
    Vector derivative = Vector::Zero(static_cast<Eigen::Index>(
      theta.size() + sigma.size() + omega_.size()));
    const double prior = prior_nll(theta, sigma, derivative);
    if (!std::isfinite(prior) || prior >= penalty()) return penalty();
    result += prior;
    return std::isfinite(result) ? result : penalty();
  }

  const std::vector<double>& omega() const { return omega_; }

 private:
  std::vector<ObjectiveTape*> tapes_;
  std::vector<double> theta_base_, sigma_base_, omega_;
  Matrix eta_;
  int n_eta_ = 0;
  WeightedEtaCollection* weighted_ = nullptr;
  std::vector<int> theta_free_, sigma_free_;
  std::vector<SaemPrior> priors_;
  bool cache_valid_ = false;
  Vector cache_point_;
  SaemEvaluation cache_;

  static double penalty() { return 1e100; }

  static std::vector<int> zero_based(std::vector<int> source) {
    for (int& value : source) {
      if (value < 1) {
        throw std::invalid_argument("An SAEM parameter index is invalid.");
      }
      --value;
    }
    return source;
  }

  static SaemEvaluation penalty_evaluation(Eigen::Index dimension) {
    SaemEvaluation result;
    result.value = penalty();
    result.gradient = Vector::Zero(dimension);
    result.native_gradient = Vector();
    return result;
  }

  void parse_priors(const Rcpp::List& config) {
    if (!config.containsElementNamed("index")) return;
    const std::vector<int> index = zero_based(
      Rcpp::as<std::vector<int>>(config["index"]));
    const std::vector<std::string> family =
      Rcpp::as<std::vector<std::string>>(config["family"]);
    const std::vector<double> mean =
      Rcpp::as<std::vector<double>>(config["mean"]);
    const std::vector<double> sd =
      Rcpp::as<std::vector<double>>(config["sd"]);
    const std::vector<double> shape =
      Rcpp::as<std::vector<double>>(config["shape"]);
    const std::vector<double> rate =
      Rcpp::as<std::vector<double>>(config["rate"]);
    if (family.size() != index.size() || mean.size() != index.size() ||
        sd.size() != index.size() || shape.size() != index.size() ||
        rate.size() != index.size()) {
      throw std::invalid_argument("The SAEM prior configuration is inconsistent.");
    }
    priors_.reserve(index.size());
    for (std::size_t prior = 0; prior < index.size(); ++prior) {
      priors_.push_back(SaemPrior{
        index[prior], family[prior], mean[prior], sd[prior],
        shape[prior], rate[prior]
      });
    }
  }

  double prior_nll(const std::vector<double>& theta,
                   const std::vector<double>& sigma,
                   Vector& derivative) const {
    std::vector<double> native;
    native.reserve(theta.size() + sigma.size() + omega_.size());
    native.insert(native.end(), theta.begin(), theta.end());
    native.insert(native.end(), sigma.begin(), sigma.end());
    native.insert(native.end(), omega_.begin(), omega_.end());
    const double log_two_pi = std::log(2.0 * std::acos(-1.0));
    double log_density = 0.0;
    for (const SaemPrior& prior : priors_) {
      if (prior.native_index < 0 ||
          prior.native_index >= static_cast<int>(native.size())) {
        throw std::invalid_argument("An SAEM prior refers to an invalid parameter.");
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
        throw std::invalid_argument("Unknown SAEM prior family.");
      }
      if (!std::isfinite(density) || !std::isfinite(gradient)) return penalty();
      log_density += density;
      derivative[prior.native_index] += gradient;
    }
    return -2.0 * log_density;
  }
};

inline Vector saem_projected_gradient(
    const Vector& point, const Vector& gradient,
    const Vector& lower, const Vector& upper) {
  Vector result = gradient;
  for (Eigen::Index index = 0; index < point.size(); ++index) {
    const double margin = 1e-12 * std::max(1.0, std::abs(point[index]));
    if ((point[index] <= lower[index] + margin && gradient[index] > 0.0) ||
        (point[index] >= upper[index] - margin && gradient[index] < 0.0)) {
      result[index] = 0.0;
    }
  }
  return result;
}

class SaemLbfgsState {
 public:
  explicit SaemLbfgsState(std::size_t memory = 5U) : memory_(memory) {}

  void prepare(Eigen::Index dimension, const Vector& proposed_scale) {
    if (dimension_ != dimension || scale_.size() != dimension) {
      clear();
      dimension_ = dimension;
      scale_ = proposed_scale;
      return;
    }
    // A persistent L-BFGS history is useful across adjacent stochastic
    // M-steps, but its coordinates cease to be meaningful after a large
    // parameter-scale change.  Refresh only in that exceptional case so
    // ordinary iterations retain their curvature information.
    bool scale_drifted = false;
    for (Eigen::Index index = 0; index < dimension; ++index) {
      const double ratio = proposed_scale[index] / scale_[index];
      if (!std::isfinite(ratio) || ratio < 0.25 || ratio > 4.0) {
        scale_drifted = true;
        break;
      }
    }
    if (scale_drifted) {
      clear();
      scale_ = proposed_scale;
    }
  }

  const Vector& scale() const { return scale_; }

  void clear() {
    s_.clear();
    y_.clear();
    rho_.clear();
  }

  Vector direction(const Vector& gradient) const {
    if (s_.empty()) return -gradient;
    Vector q = gradient;
    std::vector<double> alpha(s_.size(), 0.0);
    for (std::size_t offset = 0; offset < s_.size(); ++offset) {
      const std::size_t index = s_.size() - offset - 1U;
      alpha[index] = rho_[index] * s_[index].dot(q);
      q.noalias() -= alpha[index] * y_[index];
    }
    const double yy = y_.back().squaredNorm();
    const double gamma = yy > 0.0 ?
      std::max(1e-8, s_.back().dot(y_.back()) / yy) : 1.0;
    Vector result = gamma * q;
    for (std::size_t index = 0; index < s_.size(); ++index) {
      const double beta = rho_[index] * y_[index].dot(result);
      result.noalias() += (alpha[index] - beta) * s_[index];
    }
    return -result;
  }

  bool update(const Vector& displacement, Vector gradient_change) {
    const double ss = displacement.squaredNorm();
    if (!(ss > 0.0) || !std::isfinite(ss) || !gradient_change.allFinite()) {
      return false;
    }
    double curvature = displacement.dot(gradient_change);
    // Dampen weak or slightly negative curvature instead of either accepting
    // an unstable pair or discarding the entire memory.
    const double target = 1e-6 * ss;
    if (!std::isfinite(curvature)) return false;
    if (curvature < target) {
      gradient_change.noalias() +=
        ((target - curvature) / ss) * displacement;
      curvature = displacement.dot(gradient_change);
    }
    if (!(curvature > 0.0) || !std::isfinite(curvature)) return false;
    if (s_.size() == memory_) {
      s_.erase(s_.begin());
      y_.erase(y_.begin());
      rho_.erase(rho_.begin());
    }
    s_.push_back(displacement);
    y_.push_back(std::move(gradient_change));
    rho_.push_back(1.0 / curvature);
    return true;
  }

  std::size_t size() const { return s_.size(); }

 private:
  std::size_t memory_ = 5U;
  Eigen::Index dimension_ = 0;
  Vector scale_;
  std::vector<Vector> s_, y_;
  std::vector<double> rho_;
};

inline Rcpp::List optimize_saem_fixed_eta(
    SaemFixedEtaObjective& objective, const Rcpp::NumericVector& lower_source,
    const Rcpp::NumericVector& upper_source, int maxit,
    double tolerance, int trace, SaemLbfgsState& state) {
  const std::vector<double> start_source = objective.start();
  const Eigen::Index dimension = static_cast<Eigen::Index>(start_source.size());
  if (dimension < 1 || lower_source.size() != dimension ||
      upper_source.size() != dimension || maxit < 1 || tolerance <= 0.0 ||
      !std::isfinite(tolerance)) {
    throw std::invalid_argument("Native SAEM M-step controls are invalid.");
  }
  Vector proposed_scale(dimension);
  for (Eigen::Index index = 0; index < dimension; ++index) {
    proposed_scale[index] = std::max(
      std::abs(start_source[static_cast<std::size_t>(index)]), 1.0);
  }
  state.prepare(dimension, proposed_scale);
  const Vector scale = state.scale();
  Vector point(dimension), lower(dimension), upper(dimension);
  for (Eigen::Index index = 0; index < dimension; ++index) {
    point[index] = start_source[static_cast<std::size_t>(index)] / scale[index];
    lower[index] = lower_source[index] / scale[index];
    upper[index] = upper_source[index] / scale[index];
    if (lower[index] > upper[index] || point[index] < lower[index] ||
        point[index] > upper[index]) {
      throw std::invalid_argument("Native SAEM M-step start is outside its bounds.");
    }
  }
  auto evaluate_scaled = [&](const Vector& scaled) {
    const Vector encoded = scaled.cwiseProduct(scale);
    SaemEvaluation result = objective.evaluate(encoded);
    result.gradient = result.gradient.cwiseProduct(scale);
    return result;
  };
  SaemEvaluation current = evaluate_scaled(point);
  const double objective_scale = std::max(std::abs(current.value), 1.0);
  current.value /= objective_scale;
  current.gradient /= objective_scale;
  auto value_scaled = [&](const Vector& scaled) {
    return objective.value(scaled.cwiseProduct(scale)) / objective_scale;
  };
  int evaluations = 1;
  int gradient_evaluations = 1;
  int convergence = 1;
  int iterations = 0;
  std::string message = "iteration limit reached";
  std::vector<int> trace_iteration;
  std::vector<double> trace_value, trace_gradient, trace_step;
  for (int iteration = 0; iteration < maxit; ++iteration) {
    const Vector projected = saem_projected_gradient(
      point, current.gradient, lower, upper);
    const double norm = projected.lpNorm<Eigen::Infinity>();
    trace_iteration.push_back(iteration);
    trace_value.push_back(current.value);
    trace_gradient.push_back(norm);
    trace_step.push_back(0.0);
    if (trace > 0) {
      Rcpp::Rcout << "[LibeRation/SAEM-native] ITERATION " << iteration
                  << " OFV " << current.value
                  << " PROJECTED_GRADIENT " << norm << "\n";
    }
    if (norm <= std::max(tolerance, 1e-8) * (1.0 + std::abs(current.value))) {
      convergence = 0;
      message = "projected gradient tolerance reached";
      iterations = iteration;
      break;
    }
    Vector direction = state.direction(projected);
    for (Eigen::Index index = 0; index < dimension; ++index) {
      if (projected[index] == 0.0) direction[index] = 0.0;
    }
    double directional = current.gradient.dot(direction);
    if (!std::isfinite(directional) || directional >= -1e-14) {
      state.clear();
      direction = -projected;
      directional = current.gradient.dot(direction);
    }
    auto feasible_step = [&](const Vector& candidate_direction) {
      double maximum = 1.0;
      for (Eigen::Index index = 0; index < dimension; ++index) {
        if (candidate_direction[index] > 0.0 && std::isfinite(upper[index])) {
          maximum = std::min(
            maximum, (upper[index] - point[index]) / candidate_direction[index]);
        } else if (candidate_direction[index] < 0.0 &&
                   std::isfinite(lower[index])) {
          maximum = std::min(
            maximum, (lower[index] - point[index]) / candidate_direction[index]);
        }
      }
      return maximum;
    };
    double maximum_step = feasible_step(direction);
    double step = std::max(0.0, maximum_step);
    Vector candidate = point;
    SaemEvaluation next;
    bool accepted = false;
    int backtracks = 0;
    bool history_direction = state.size() > 0U;
    for (int line_search = 0; line_search < 20 && step > 1e-16; ++line_search) {
      if (line_search == 2 && history_direction) {
        // Adjacent SAEM M-steps change the fixed-ETA objective.  Retained
        // curvature is useful only while it predicts a readily acceptable
        // step; otherwise restart promptly from projected steepest descent.
        state.clear();
        direction = -projected;
        directional = current.gradient.dot(direction);
        maximum_step = feasible_step(direction);
        step = std::max(0.0, maximum_step);
        history_direction = false;
      }
      candidate = (point + step * direction).cwiseMax(lower).cwiseMin(upper);
      const Vector actual = candidate - point;
      const double armijo_directional = current.gradient.dot(actual);
      const double trial = value_scaled(candidate);
      ++evaluations;
      if (std::isfinite(trial) && trial < 1e100 &&
          trial <= current.value + 1e-4 * armijo_directional) {
        next = evaluate_scaled(candidate);
        next.value /= objective_scale;
        next.gradient /= objective_scale;
        ++gradient_evaluations;
        accepted = true;
        break;
      }
      ++backtracks;
      double interpolated = step * 0.5;
      const double denominator = 2.0 *
        (trial - current.value - step * directional);
      if (std::isfinite(denominator) && denominator > 0.0) {
        const double quadratic = -directional * step * step / denominator;
        if (std::isfinite(quadratic)) {
          interpolated = std::max(0.1 * step, std::min(0.5 * step, quadratic));
        }
      }
      step = interpolated;
    }
    if (!accepted) {
      convergence = 52;
      message = "line search failed";
      iterations = iteration;
      break;
    }
    const Vector displacement = candidate - point;
    const Vector change = next.gradient - current.gradient;
    const double accepted_directional = next.gradient.dot(direction);
    const bool strong_wolfe = std::isfinite(accepted_directional) &&
      std::abs(accepted_directional) <= 0.9 * std::abs(directional);
    if (strong_wolfe) state.update(displacement, change);
    else state.clear();
    if (backtracks > 12) state.clear();
    const double previous = current.value;
    point = candidate;
    current = next;
    iterations = iteration + 1;
    trace_step.back() = step;
    if (std::abs(previous - current.value) <=
        tolerance * (1.0 + std::abs(current.value))) {
      const Vector projected_next = saem_projected_gradient(
        point, current.gradient, lower, upper);
      if (projected_next.lpNorm<Eigen::Infinity>() <=
          std::sqrt(tolerance) * (1.0 + std::abs(current.value))) {
        convergence = 0;
        message = "relative objective and gradient tolerance reached";
        break;
      }
    }
    if ((iteration + 1) % 10 == 0) Rcpp::checkUserInterrupt();
  }
  const Vector encoded = point.cwiseProduct(scale);
  std::vector<double> theta, sigma;
  objective.decode(encoded, theta, sigma);
  for (double& value : trace_value) value *= objective_scale;
  Rcpp::NumericVector par(dimension), gradient(dimension);
  for (Eigen::Index index = 0; index < dimension; ++index) {
    par[index] = encoded[index];
    gradient[index] = current.gradient[index] * objective_scale / scale[index];
  }
  const Rcpp::IntegerVector counts = Rcpp::IntegerVector::create(
    Rcpp::Named("function") = evaluations,
    Rcpp::Named("gradient") = gradient_evaluations);
  Rcpp::NumericVector native_gradient =
    current.native_gradient.size() ?
      libertad::eigen_vector_to_r(current.native_gradient) :
      Rcpp::NumericVector();
  return Rcpp::List::create(
    Rcpp::Named("par") = par,
    Rcpp::Named("theta") = Rcpp::wrap(theta),
    Rcpp::Named("sigma") = Rcpp::wrap(sigma),
    Rcpp::Named("omega") = Rcpp::wrap(objective.omega()),
    Rcpp::Named("value") = current.value * objective_scale,
    Rcpp::Named("convergence") = convergence,
    Rcpp::Named("message") = message,
    Rcpp::Named("counts") = counts,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("objective_evaluations") = evaluations,
    Rcpp::Named("gradient_evaluations") = gradient_evaluations,
    Rcpp::Named("objective_scale") = objective_scale,
    Rcpp::Named("lbfgs_memory") = static_cast<int>(state.size()),
    Rcpp::Named("gradient") = gradient,
    Rcpp::Named("native_gradient") = native_gradient,
    Rcpp::Named("backend") = "native-cpp-fixed-eta-lbfgs",
    Rcpp::Named("telemetry") = Rcpp::DataFrame::create(
      Rcpp::Named("iteration") = trace_iteration,
      Rcpp::Named("objective") = trace_value,
      Rcpp::Named("projected_gradient") = trace_gradient,
      Rcpp::Named("step") = trace_step));
}
