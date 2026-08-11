#pragma once

// This header is included near the end of pk_engine.cpp, inside namespace
// liberation, after ObjectiveTape and the shared parameter helpers are defined.
// It deliberately contains no R callbacks: a complete HMC/NUTS trajectory
// stays in C++ and only checks for cancellation between sampler iterations.

struct HmcParameters {
  std::vector<double> theta;
  std::vector<double> sigma;
  std::vector<double> omega;
  Matrix transform;
  Matrix eta_factor;
  std::vector<Matrix> eta_factor_gradient;
  double eta_log_jacobian = 0.0;
  Vector eta_log_jacobian_gradient;
  double log_jacobian = 0.0;
  Vector log_jacobian_gradient;
};

struct HmcEvaluation {
  double logp = -std::numeric_limits<double>::infinity();
  Vector gradient;
  HmcParameters parameters;
  Vector eta;
  std::vector<double> outer;
  bool finite = false;
};

class HmcTarget {
 public:
  HmcTarget(ObjectiveTape& tape, const Rcpp::List& config)
      : tape_(tape) {
    theta_base_ = Rcpp::as<std::vector<double>>(config["theta"]);
    sigma_base_ = Rcpp::as<std::vector<double>>(config["sigma"]);
    omega_base_ = Rcpp::as<std::vector<double>>(config["omega"]);
    theta_free_ = zero_based(Rcpp::as<std::vector<int>>(config["theta_free"]));
    sigma_free_ = zero_based(Rcpp::as<std::vector<int>>(config["sigma_free"]));
    omega_free_ = zero_based(Rcpp::as<std::vector<int>>(config["omega_free"]));
    omega_full_ = Rcpp::as<bool>(config["omega_full"]);
    omega_rows_ = zero_based(Rcpp::as<std::vector<int>>(config["omega_rows"]));
    omega_cols_ = zero_based(Rcpp::as<std::vector<int>>(config["omega_cols"]));
    n_eta_base_ = Rcpp::as<int>(config["n_eta_base"]);
    n_subjects_ = Rcpp::as<int>(config["n_subjects"]);
    n_eta_ = Rcpp::as<int>(config["n_eta"]);
    noncentered_ = Rcpp::as<bool>(config["noncentered"]);
    lower_ = Rcpp::as<std::vector<double>>(config["lower"]);
    upper_ = Rcpp::as<std::vector<double>>(config["upper"]);
    initial_ = Rcpp::as<std::vector<double>>(config["initial"]);

    const std::vector<int> prior_index = zero_based_allow_empty(
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
      throw std::invalid_argument("Native HMC prior mapping is inconsistent.");
    }
    priors_.reserve(prior_count);
    for (std::size_t index = 0; index < prior_count; ++index) {
      priors_.push_back(PopulationPrior{
        prior_index[index], prior_family[index], prior_mean[index],
        prior_sd[index], prior_shape[index], prior_rate[index]
      });
    }

    n_outer_ = theta_free_.size() + sigma_free_.size() +
      (omega_full_ && !omega_free_.empty() ? omega_base_.size() :
       omega_free_.size());
    n_eta_total_ = static_cast<std::size_t>(n_subjects_) *
      static_cast<std::size_t>(n_eta_);
    if (lower_.size() != n_outer_ || upper_.size() != n_outer_ ||
        initial_.size() != n_outer_ + n_eta_total_) {
      throw std::invalid_argument("Native HMC parameter dimensions are inconsistent.");
    }
    if (omega_rows_.size() != omega_base_.size() ||
        omega_cols_.size() != omega_base_.size()) {
      throw std::invalid_argument("Native HMC OMEGA mapping is inconsistent.");
    }
    const std::size_t expected_domain = theta_base_.size() + n_eta_total_ +
      sigma_base_.size() + omega_base_.size();
    if (tape_.domain_names.size() != expected_domain) {
      throw std::invalid_argument("Native HMC objective tape has the wrong domain.");
    }
    if (noncentered_ && n_eta_ != n_eta_base_) {
      throw std::invalid_argument(
        "Native non-centred HMC requires the base ETA layout."
      );
    }
  }

  std::size_t dimension() const { return n_outer_ + n_eta_total_; }
  std::size_t n_outer() const { return n_outer_; }
  const std::vector<double>& initial() const { return initial_; }
  long long evaluations() const { return evaluations_; }

  HmcEvaluation evaluate(const Vector& q) {
    ++evaluations_;
    HmcEvaluation result;
    result.gradient = Vector::Constant(
      static_cast<Eigen::Index>(dimension()),
      std::numeric_limits<double>::quiet_NaN());
    if (static_cast<std::size_t>(q.size()) != dimension() ||
        !q.allFinite()) {
      return result;
    }

    BoundResult bounded = bound_forward(q.head(
      static_cast<Eigen::Index>(n_outer_)));
    if (!bounded.finite) return result;
    HmcParameters parameters;
    try {
      parameters = decode(bounded.value);
    } catch (...) {
      return result;
    }

    const Vector latent_eta = q.tail(static_cast<Eigen::Index>(n_eta_total_));
    result.eta = latent_eta;
    if (noncentered_) {
      for (int subject = 0; subject < n_subjects_; ++subject) {
        const Eigen::Index offset = static_cast<Eigen::Index>(subject * n_eta_);
        result.eta.segment(offset, n_eta_) = parameters.eta_factor *
          latent_eta.segment(offset, n_eta_);
      }
    }
    std::vector<double> point;
    point.reserve(tape_.domain_names.size());
    point.insert(point.end(), parameters.theta.begin(), parameters.theta.end());
    for (Eigen::Index index = 0; index < result.eta.size(); ++index) {
      point.push_back(result.eta[index]);
    }
    point.insert(point.end(), parameters.sigma.begin(), parameters.sigma.end());
    point.insert(point.end(), parameters.omega.begin(), parameters.omega.end());

    std::ostringstream messages;
    const std::vector<double> objective = tape_.fun.Forward(0, point, messages);
    if (objective.size() != 1U || !std::isfinite(objective[0]) ||
        tape_.fun.compare_change_number() != 0U) {
      return result;
    }
    const std::vector<double> weight(1, 1.0);
    const std::vector<double> objective_gradient = tape_.fun.Reverse(1, weight);
    if (objective_gradient.size() != point.size()) return result;

    Vector prior_gradient;
    const double prior = prior_nll(parameters, &prior_gradient);
    if (!std::isfinite(prior) || prior >= penalty()) return result;

    const std::size_t n_native = theta_base_.size() + sigma_base_.size() +
      omega_base_.size();
    Vector native_gradient(static_cast<Eigen::Index>(n_native));
    std::size_t native = 0U;
    for (std::size_t index = 0; index < theta_base_.size(); ++index) {
      native_gradient[static_cast<Eigen::Index>(native++)] =
        -0.5 * objective_gradient[index];
    }
    const std::size_t sigma_objective_offset = theta_base_.size() + n_eta_total_;
    for (std::size_t index = 0; index < sigma_base_.size(); ++index) {
      native_gradient[static_cast<Eigen::Index>(native++)] =
        -0.5 * objective_gradient[sigma_objective_offset + index];
    }
    const std::size_t omega_objective_offset =
      sigma_objective_offset + sigma_base_.size();
    for (std::size_t index = 0; index < omega_base_.size(); ++index) {
      native_gradient[static_cast<Eigen::Index>(native++)] =
        -0.5 * objective_gradient[omega_objective_offset + index];
    }
    native_gradient -= 0.5 * prior_gradient;

    Vector outer_gradient = parameters.transform.transpose() * native_gradient +
      parameters.log_jacobian_gradient;
    const std::size_t eta_offset = theta_base_.size();
    if (noncentered_) {
      for (std::size_t outer = 0; outer < n_outer_; ++outer) {
        const Matrix& derivative = parameters.eta_factor_gradient[outer];
        if (derivative.squaredNorm() == 0.0) continue;
        for (int subject = 0; subject < n_subjects_; ++subject) {
          const Eigen::Index offset = static_cast<Eigen::Index>(subject * n_eta_);
          const Vector eta_gradient = -0.5 * Eigen::Map<const Vector>(
            objective_gradient.data() + eta_offset +
              static_cast<std::size_t>(offset), n_eta_
          );
          outer_gradient[static_cast<Eigen::Index>(outer)] +=
            eta_gradient.dot(
              derivative * latent_eta.segment(offset, n_eta_)
            );
        }
      }
      outer_gradient += parameters.eta_log_jacobian_gradient;
    }
    outer_gradient = outer_gradient.cwiseProduct(bounded.derivative) +
      bounded.log_jacobian_gradient;
    result.gradient.head(static_cast<Eigen::Index>(n_outer_)) = outer_gradient;
    if (noncentered_) {
      for (int subject = 0; subject < n_subjects_; ++subject) {
        const Eigen::Index offset = static_cast<Eigen::Index>(subject * n_eta_);
        const Vector eta_gradient = -0.5 * Eigen::Map<const Vector>(
          objective_gradient.data() + eta_offset +
            static_cast<std::size_t>(offset), n_eta_
        );
        result.gradient.segment(
          static_cast<Eigen::Index>(n_outer_) + offset, n_eta_
        ) = parameters.eta_factor.transpose() * eta_gradient;
      }
    } else {
      for (std::size_t index = 0; index < n_eta_total_; ++index) {
        result.gradient[static_cast<Eigen::Index>(n_outer_ + index)] =
          -0.5 * objective_gradient[eta_offset + index];
      }
    }

    result.logp = -0.5 * objective[0] - 0.5 * prior +
      parameters.log_jacobian + bounded.log_jacobian +
      parameters.eta_log_jacobian;
    result.parameters = std::move(parameters);
    result.outer.assign(
      bounded.value.data(),
      bounded.value.data() + bounded.value.size());
    result.finite = std::isfinite(result.logp) && result.gradient.allFinite();
    if (!result.finite) {
      result.logp = -std::numeric_limits<double>::infinity();
    }
    return result;
  }

  Rcpp::NumericVector native_row(const HmcEvaluation& evaluated) const {
    const std::size_t size = theta_base_.size() + sigma_base_.size() +
      omega_base_.size() + n_eta_total_ + 1U;
    Rcpp::NumericVector result(static_cast<R_xlen_t>(size));
    std::size_t cursor = 0U;
    for (double value : evaluated.parameters.theta) result[cursor++] = value;
    for (double value : evaluated.parameters.sigma) result[cursor++] = value;
    for (double value : evaluated.parameters.omega) result[cursor++] = value;
    for (Eigen::Index index = 0; index < evaluated.eta.size(); ++index) {
      result[cursor++] = evaluated.eta[index];
    }
    result[cursor] = evaluated.logp;
    return result;
  }

 private:
  struct BoundResult {
    Vector value;
    Vector derivative;
    Vector log_jacobian_gradient;
    double log_jacobian = 0.0;
    bool finite = true;
  };

  ObjectiveTape& tape_;
  std::vector<double> theta_base_, sigma_base_, omega_base_;
  std::vector<int> theta_free_, sigma_free_, omega_free_;
  std::vector<int> omega_rows_, omega_cols_;
  std::vector<double> lower_, upper_, initial_;
  std::vector<PopulationPrior> priors_;
  bool omega_full_ = false;
  bool noncentered_ = false;
  int n_eta_base_ = 0;
  int n_subjects_ = 0;
  int n_eta_ = 0;
  std::size_t n_outer_ = 0U;
  std::size_t n_eta_total_ = 0U;
  long long evaluations_ = 0;

  static double penalty() { return 1e100; }

  static std::vector<int> zero_based(std::vector<int> source) {
    for (int& value : source) {
      if (value < 1) {
        throw std::invalid_argument("A native HMC parameter index is invalid.");
      }
      --value;
    }
    return source;
  }

  static std::vector<int> zero_based_allow_empty(std::vector<int> source) {
    if (source.empty()) return source;
    return zero_based(std::move(source));
  }

  static double softplus(double value) {
    if (value > 0.0) return value + std::log1p(std::exp(-value));
    return std::log1p(std::exp(value));
  }

  BoundResult bound_forward(const Vector& q) const {
    BoundResult result;
    result.value.resize(q.size());
    result.derivative.resize(q.size());
    result.log_jacobian_gradient.resize(q.size());
    for (Eigen::Index index = 0; index < q.size(); ++index) {
      const double current = q[index];
      const double lo = lower_[static_cast<std::size_t>(index)];
      const double hi = upper_[static_cast<std::size_t>(index)];
      if (std::isfinite(lo) && std::isfinite(hi)) {
        const double width = hi - lo;
        const double p = current >= 0.0 ?
          1.0 / (1.0 + std::exp(-current)) :
          std::exp(current) / (1.0 + std::exp(current));
        result.value[index] = lo + width * p;
        result.derivative[index] = width * p * (1.0 - p);
        result.log_jacobian += std::log(width) -
          softplus(-current) - softplus(current);
        result.log_jacobian_gradient[index] = 1.0 - 2.0 * p;
      } else if (std::isfinite(lo)) {
        const double scale = std::exp(current);
        result.value[index] = lo + scale;
        result.derivative[index] = scale;
        result.log_jacobian += current;
        result.log_jacobian_gradient[index] = 1.0;
      } else if (std::isfinite(hi)) {
        const double scale = std::exp(current);
        result.value[index] = hi - scale;
        result.derivative[index] = -scale;
        result.log_jacobian += current;
        result.log_jacobian_gradient[index] = 1.0;
      } else {
        result.value[index] = current;
        result.derivative[index] = 1.0;
        result.log_jacobian_gradient[index] = 0.0;
      }
    }
    result.finite = result.value.allFinite() && result.derivative.allFinite() &&
      result.log_jacobian_gradient.allFinite() &&
      std::isfinite(result.log_jacobian);
    return result;
  }

  HmcParameters decode(const Vector& encoded) const {
    HmcParameters result;
    result.theta = theta_base_;
    result.sigma = sigma_base_;
    result.omega = omega_base_;
    const int n_native = static_cast<int>(
      theta_base_.size() + sigma_base_.size() + omega_base_.size());
    result.transform = Matrix::Zero(n_native, encoded.size());
    result.log_jacobian_gradient = Vector::Zero(encoded.size());
    result.eta_factor_gradient.assign(
      static_cast<std::size_t>(encoded.size()),
      Matrix::Zero(n_eta_base_, n_eta_base_)
    );
    result.eta_log_jacobian_gradient = Vector::Zero(encoded.size());
    std::size_t cursor = 0U;
    for (int index : theta_free_) {
      result.theta[static_cast<std::size_t>(index)] =
        encoded[static_cast<Eigen::Index>(cursor)];
      result.transform(index, static_cast<Eigen::Index>(cursor)) = 1.0;
      ++cursor;
    }
    const int sigma_offset = static_cast<int>(theta_base_.size());
    for (int index : sigma_free_) {
      const double value = std::exp(encoded[static_cast<Eigen::Index>(cursor)]);
      result.sigma[static_cast<std::size_t>(index)] = value;
      result.transform(sigma_offset + index, static_cast<Eigen::Index>(cursor)) =
        value;
      result.log_jacobian += std::log(value);
      result.log_jacobian_gradient[static_cast<Eigen::Index>(cursor)] = 1.0;
      ++cursor;
    }
    const int omega_offset = sigma_offset + static_cast<int>(sigma_base_.size());
    const std::size_t omega_cursor = cursor;
    if (omega_full_ && !omega_free_.empty()) {
      Matrix lower = Matrix::Zero(n_eta_base_, n_eta_base_);
      for (std::size_t entry = 0; entry < omega_base_.size(); ++entry) {
        const int row = omega_rows_[entry];
        const int column = omega_cols_[entry];
        const double encoded_value =
          encoded[static_cast<Eigen::Index>(cursor + entry)];
        lower(row, column) = row == column ? std::exp(encoded_value) :
          encoded_value;
      }
      const Matrix covariance = lower * lower.transpose();
      for (std::size_t entry = 0; entry < omega_base_.size(); ++entry) {
        result.omega[entry] =
          covariance(omega_rows_[entry], omega_cols_[entry]);
      }
      result.log_jacobian += static_cast<double>(n_eta_base_) * std::log(2.0);
      for (int row = 0; row < n_eta_base_; ++row) {
        result.log_jacobian += static_cast<double>(n_eta_base_ + 1 - row) *
          std::log(lower(row, row));
      }
      for (std::size_t encoded_entry = 0;
           encoded_entry < omega_base_.size(); ++encoded_entry) {
        Matrix derivative_lower = Matrix::Zero(n_eta_base_, n_eta_base_);
        const int row = omega_rows_[encoded_entry];
        const int column = omega_cols_[encoded_entry];
        derivative_lower(row, column) =
          row == column ? lower(row, column) : 1.0;
        const Matrix derivative = derivative_lower * lower.transpose() +
          lower * derivative_lower.transpose();
        for (std::size_t native = 0; native < omega_base_.size(); ++native) {
          result.transform(
            omega_offset + static_cast<int>(native),
            static_cast<Eigen::Index>(cursor + encoded_entry)) =
              derivative(omega_rows_[native], omega_cols_[native]);
        }
        if (row == column) {
          result.log_jacobian_gradient[
            static_cast<Eigen::Index>(cursor + encoded_entry)] =
              static_cast<double>(n_eta_base_ + 1 - row);
        }
      }
      cursor += omega_base_.size();
    } else {
      for (int index : omega_free_) {
        const double value =
          std::exp(encoded[static_cast<Eigen::Index>(cursor)]);
        result.omega[static_cast<std::size_t>(index)] = value;
        result.transform(
          omega_offset + index, static_cast<Eigen::Index>(cursor)) = value;
        result.log_jacobian += std::log(value);
        result.log_jacobian_gradient[static_cast<Eigen::Index>(cursor)] = 1.0;
        ++cursor;
      }
    }
    if (cursor != static_cast<std::size_t>(encoded.size())) {
      throw std::invalid_argument("Native HMC encoded parameter length is invalid.");
    }
    if (noncentered_) {
      Matrix covariance = Matrix::Zero(n_eta_base_, n_eta_base_);
      for (std::size_t entry = 0; entry < result.omega.size(); ++entry) {
        const int row = omega_rows_[entry];
        const int column = omega_cols_[entry];
        covariance(row, column) = result.omega[entry];
        covariance(column, row) = result.omega[entry];
      }
      Eigen::LLT<Matrix> decomposition(covariance);
      if (decomposition.info() != Eigen::Success) {
        throw std::invalid_argument(
          "Native non-centred HMC requires positive-definite OMEGA."
        );
      }
      result.eta_factor = decomposition.matrixL();
      if (omega_full_ && !omega_free_.empty()) {
        for (std::size_t entry = 0; entry < omega_base_.size(); ++entry) {
          const int row = omega_rows_[entry];
          const int column = omega_cols_[entry];
          result.eta_factor_gradient[omega_cursor + entry](row, column) =
            row == column ? result.eta_factor(row, column) : 1.0;
        }
      } else if (!omega_full_ && !omega_free_.empty()) {
        for (std::size_t entry = 0; entry < omega_free_.size(); ++entry) {
          const int native = omega_free_[entry];
          const int row = omega_rows_[static_cast<std::size_t>(native)];
          result.eta_factor_gradient[omega_cursor + entry](row, row) =
            0.5 * result.eta_factor(row, row);
        }
      }
      result.eta_log_jacobian = static_cast<double>(n_subjects_) *
        result.eta_factor.diagonal().array().log().sum();
      for (std::size_t outer = 0; outer < n_outer_; ++outer) {
        result.eta_log_jacobian_gradient[
          static_cast<Eigen::Index>(outer)] =
            static_cast<double>(n_subjects_) *
            result.eta_factor.triangularView<Eigen::Lower>().solve(
              result.eta_factor_gradient[outer]
            ).trace();
      }
    }
    return result;
  }

  double prior_nll(const HmcParameters& parameters,
                   Vector* derivative = nullptr) const {
    std::vector<double> native;
    native.reserve(theta_base_.size() + sigma_base_.size() + omega_base_.size());
    native.insert(native.end(), parameters.theta.begin(), parameters.theta.end());
    native.insert(native.end(), parameters.sigma.begin(), parameters.sigma.end());
    native.insert(native.end(), parameters.omega.begin(), parameters.omega.end());
    if (derivative) {
      derivative->setZero(static_cast<Eigen::Index>(native.size()));
    }
    const double log_two_pi = std::log(2.0 * std::acos(-1.0));
    double log_density = 0.0;
    for (const PopulationPrior& prior : priors_) {
      if (prior.native_index < 0 ||
          prior.native_index >= static_cast<int>(native.size())) {
        return penalty();
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
          gradient = 2.0 / value +
            2.0 * (std::log(value) - prior.mean) /
            (prior.sd * prior.sd * value);
        }
      } else if (prior.family == "inverse_gamma") {
        if (value > 0.0 && prior.shape > 0.0 && prior.rate > 0.0) {
          density = prior.shape * std::log(prior.rate) -
            std::lgamma(prior.shape) -
            (prior.shape + 1.0) * std::log(value) - prior.rate / value;
          gradient = 2.0 * (prior.shape + 1.0) / value -
            2.0 * prior.rate / (value * value);
        }
      } else {
        throw std::invalid_argument("Unknown native HMC prior family.");
      }
      if (!std::isfinite(density) || !std::isfinite(gradient)) return penalty();
      log_density += density;
      if (derivative) (*derivative)[prior.native_index] += gradient;
    }
    return -2.0 * log_density;
  }
};

class HmcRandom {
 public:
  explicit HmcRandom(std::uint64_t seed) : engine_(seed) {}

  double uniform() {
    const long double numerator =
      static_cast<long double>(engine_()) + 0.5L;
    const long double denominator =
      static_cast<long double>(std::numeric_limits<std::uint64_t>::max()) + 1.0L;
    return static_cast<double>(numerator / denominator);
  }

  double normal() {
    if (has_spare_) {
      has_spare_ = false;
      return spare_;
    }
    const double radius = std::sqrt(-2.0 * std::log(uniform()));
    const double angle = 2.0 * std::acos(-1.0) * uniform();
    spare_ = radius * std::sin(angle);
    has_spare_ = true;
    return radius * std::cos(angle);
  }

  double exponential() { return -std::log(uniform()); }

  bool bernoulli(double probability) {
    return uniform() < std::min(std::max(probability, 0.0), 1.0);
  }

  int direction() { return uniform() < 0.5 ? -1 : 1; }

 private:
  std::mt19937_64 engine_;
  bool has_spare_ = false;
  double spare_ = 0.0;
};

struct HmcMetric {
  std::string kind = "diagonal";
  Matrix mass;
  Matrix inverse_mass;
  Matrix root;
  std::vector<std::vector<int>> blocks;

  HmcMetric() = default;
  HmcMetric(int dimension, const std::string& metric_kind,
            const std::vector<std::vector<int>>& metric_blocks)
      : kind(metric_kind), mass(Matrix::Identity(dimension, dimension)),
        inverse_mass(Matrix::Identity(dimension, dimension)),
        root(Matrix::Identity(dimension, dimension)), blocks(metric_blocks) {}
};

inline HmcMetric hmc_metric_from_mass(
    const Matrix& proposed_mass, const std::string& kind,
    const std::vector<std::vector<int>>& blocks) {
  const Matrix symmetric = 0.5 * (proposed_mass + proposed_mass.transpose());
  Eigen::SelfAdjointEigenSolver<Matrix> decomposition(symmetric);
  if (decomposition.info() != Eigen::Success) {
    throw std::runtime_error("HMC metric eigendecomposition failed.");
  }
  Vector values = decomposition.eigenvalues().array().max(1e-3).min(1e3);
  HmcMetric result;
  result.kind = kind;
  result.blocks = blocks;
  result.mass = decomposition.eigenvectors() * values.asDiagonal() *
    decomposition.eigenvectors().transpose();
  result.inverse_mass = decomposition.eigenvectors() *
    values.cwiseInverse().asDiagonal() * decomposition.eigenvectors().transpose();
  Eigen::LLT<Matrix> cholesky(result.mass);
  if (cholesky.info() != Eigen::Success) {
    throw std::runtime_error("Regularized HMC metric is not positive definite.");
  }
  result.root = cholesky.matrixL();
  return result;
}

inline Vector hmc_momentum(const HmcMetric& metric, HmcRandom& random) {
  Vector standard(metric.mass.rows());
  for (Eigen::Index index = 0; index < standard.size(); ++index) {
    standard[index] = random.normal();
  }
  return metric.root * standard;
}

inline Vector hmc_velocity(const Vector& momentum, const HmcMetric& metric) {
  return metric.inverse_mass * momentum;
}

struct HmcMove {
  Vector q;
  Vector momentum;
  HmcEvaluation evaluated;
  bool valid = false;
};

inline double hmc_kinetic(const Vector& momentum, const HmcMetric& metric) {
  return 0.5 * momentum.dot(hmc_velocity(momentum, metric));
}

HmcMove hmc_leapfrog(const Vector& q, const Vector& momentum,
                     const Vector& gradient, double epsilon,
                     const HmcMetric& metric, HmcTarget& target) {
  HmcMove result;
  result.momentum = momentum + 0.5 * epsilon * gradient;
  result.q = q + epsilon * hmc_velocity(result.momentum, metric);
  result.evaluated = target.evaluate(result.q);
  if (!result.evaluated.finite) return result;
  result.momentum += 0.5 * epsilon * result.evaluated.gradient;
  result.valid = result.momentum.allFinite();
  return result;
}

double hmc_find_step(const Vector& q, const HmcEvaluation& evaluated,
                     const HmcMetric& metric, HmcTarget& target,
                     HmcRandom& random) {
  double epsilon = 1.0;
  Vector momentum = hmc_momentum(metric, random);
  HmcMove proposal = hmc_leapfrog(
    q, momentum, evaluated.gradient, epsilon, metric, target);
  double log_accept = proposal.valid ?
    proposal.evaluated.logp - hmc_kinetic(proposal.momentum, metric) -
      evaluated.logp + hmc_kinetic(momentum, metric) :
    -std::numeric_limits<double>::infinity();
  const int direction =
    std::isfinite(log_accept) && log_accept > std::log(0.5) ? 1 : -1;
  for (int iteration = 0; iteration < 20; ++iteration) {
    const double candidate = epsilon * (direction > 0 ? 2.0 : 0.5);
    proposal = hmc_leapfrog(
      q, momentum, evaluated.gradient, candidate, metric, target);
    const double candidate_accept = proposal.valid ?
      proposal.evaluated.logp - hmc_kinetic(proposal.momentum, metric) -
        evaluated.logp + hmc_kinetic(momentum, metric) :
      -std::numeric_limits<double>::infinity();
    const bool keep_going = direction > 0 ?
      candidate_accept > std::log(0.5) :
      candidate_accept < std::log(0.5);
    epsilon = candidate;
    if (!keep_going || epsilon < 1e-8 || epsilon > 1e2) break;
  }
  return std::min(std::max(epsilon, 1e-8), 1e2);
}

class DualAverage {
 public:
  DualAverage(double initial, double target)
      : mu_(std::log(10.0 * initial)),
        log_step_(std::log(initial)),
        log_step_bar_(std::log(initial)),
        target_(target) {}

  double update(double acceptance) {
    ++iteration_;
    const double current = static_cast<double>(iteration_);
    hbar_ = (1.0 - 1.0 / (current + 10.0)) * hbar_ +
      (target_ - acceptance) / (current + 10.0);
    log_step_ = mu_ - std::sqrt(current) / 0.05 * hbar_;
    const double weight = std::pow(current, -0.75);
    log_step_bar_ = weight * log_step_ +
      (1.0 - weight) * log_step_bar_;
    return std::exp(log_step_);
  }

  double final() const { return std::exp(log_step_bar_); }

 private:
  double mu_;
  double log_step_;
  double log_step_bar_;
  double target_;
  double hbar_ = 0.0;
  int iteration_ = 0;
};

struct HmcTransition {
  Vector q;
  HmcEvaluation evaluated;
  double acceptance = 0.0;
  bool accepted = false;
  bool divergence = false;
  double energy_error = std::numeric_limits<double>::quiet_NaN();
  double energy = std::numeric_limits<double>::quiet_NaN();
  int tree_depth = -1;
  int leapfrog = 0;
  bool max_depth_reached = false;
};

HmcTransition hmc_transition(const Vector& q,
                             const HmcEvaluation& evaluated,
                             double step_size, const HmcMetric& metric,
                             int n_leapfrog, HmcTarget& target,
                             double divergence_threshold,
                             HmcRandom& random) {
  Vector momentum = hmc_momentum(metric, random);
  const Vector initial_momentum = momentum;
  const double energy = -evaluated.logp + hmc_kinetic(initial_momentum, metric);
  Vector proposal_q = q;
  HmcEvaluation proposal = evaluated;
  bool valid = true;
  int completed = 0;
  for (int step = 0; step < n_leapfrog; ++step) {
    HmcMove moved = hmc_leapfrog(
      proposal_q, momentum, proposal.gradient, step_size, metric, target);
    ++completed;
    if (!moved.valid) {
      valid = false;
      break;
    }
    proposal_q = std::move(moved.q);
    momentum = std::move(moved.momentum);
    proposal = std::move(moved.evaluated);
  }
  const double energy_error = valid ?
    -(proposal.logp - hmc_kinetic(momentum, metric)) +
      (evaluated.logp - hmc_kinetic(initial_momentum, metric)) :
    std::numeric_limits<double>::infinity();
  const double acceptance = std::isfinite(energy_error) ?
    std::min(1.0, std::exp(-energy_error)) : 0.0;
  const bool accepted = valid && random.bernoulli(acceptance);
  HmcTransition result;
  result.q = accepted ? proposal_q : q;
  result.evaluated = accepted ? proposal : evaluated;
  result.acceptance = acceptance;
  result.accepted = accepted;
  result.divergence = !std::isfinite(energy_error) ||
    std::abs(energy_error) > divergence_threshold;
  result.energy_error = energy_error;
  result.energy = energy;
  result.leapfrog = completed;
  return result;
}

inline bool nuts_stop(const Vector& q_minus, const Vector& q_plus,
                      const Vector& r_minus, const Vector& r_plus,
                      const HmcMetric& metric) {
  const Vector delta = q_plus - q_minus;
  return delta.dot(hmc_velocity(r_minus, metric)) >= 0.0 &&
    delta.dot(hmc_velocity(r_plus, metric)) >= 0.0;
}

struct NutsTree {
  Vector q_minus, r_minus, q_plus, r_plus, q_proposal;
  HmcEvaluation e_minus, e_plus, e_proposal;
  int n = 0;
  bool active = false;
  double alpha = 0.0;
  int n_alpha = 0;
  bool divergent = false;
  int leapfrog = 0;
};

NutsTree nuts_tree(const Vector& q, const Vector& r,
                   const HmcEvaluation& evaluated, double log_slice,
                   int direction, int depth, double step_size,
                   const HmcMetric& metric, HmcTarget& target, double joint0,
                   double divergence_threshold, HmcRandom& random) {
  if (depth == 0) {
    HmcMove moved = hmc_leapfrog(
      q, r, evaluated.gradient, direction * step_size, metric, target);
    NutsTree result;
    result.q_minus = result.q_plus = moved.q;
    result.r_minus = result.r_plus = moved.momentum;
    result.e_minus = result.e_plus = moved.evaluated;
    result.q_proposal = moved.valid ? moved.q : q;
    result.e_proposal = moved.valid ? moved.evaluated : evaluated;
    result.n_alpha = 1;
    result.leapfrog = 1;
    if (!moved.valid) {
      result.divergent = true;
      return result;
    }
    const double joint =
      moved.evaluated.logp - hmc_kinetic(moved.momentum, metric);
    const double error = joint0 - joint;
    result.n = log_slice <= joint ? 1 : 0;
    result.active = std::isfinite(joint) &&
      log_slice - divergence_threshold < joint;
    result.alpha = std::min(1.0, std::exp(std::min(0.0, joint - joint0)));
    result.divergent = !std::isfinite(error) ||
      std::abs(error) > divergence_threshold;
    return result;
  }

  NutsTree left = nuts_tree(
    q, r, evaluated, log_slice, direction, depth - 1, step_size, metric,
    target, joint0, divergence_threshold, random);
  if (!left.active) return left;
  NutsTree right = direction < 0 ?
    nuts_tree(
      left.q_minus, left.r_minus, left.e_minus, log_slice, direction,
      depth - 1, step_size, metric, target, joint0, divergence_threshold,
      random) :
    nuts_tree(
      left.q_plus, left.r_plus, left.e_plus, log_slice, direction,
      depth - 1, step_size, metric, target, joint0, divergence_threshold,
      random);

  NutsTree result;
  result.q_proposal = left.q_proposal;
  result.e_proposal = left.e_proposal;
  if (right.n > 0 &&
      random.bernoulli(
        static_cast<double>(right.n) /
        static_cast<double>(std::max(left.n + right.n, 1)))) {
    result.q_proposal = right.q_proposal;
    result.e_proposal = right.e_proposal;
  }
  result.q_minus = direction < 0 ? right.q_minus : left.q_minus;
  result.r_minus = direction < 0 ? right.r_minus : left.r_minus;
  result.e_minus = direction < 0 ? right.e_minus : left.e_minus;
  result.q_plus = direction < 0 ? left.q_plus : right.q_plus;
  result.r_plus = direction < 0 ? left.r_plus : right.r_plus;
  result.e_plus = direction < 0 ? left.e_plus : right.e_plus;
  result.n = left.n + right.n;
  result.active = right.active && nuts_stop(
    result.q_minus, result.q_plus, result.r_minus, result.r_plus, metric);
  result.alpha = left.alpha + right.alpha;
  result.n_alpha = left.n_alpha + right.n_alpha;
  result.divergent = left.divergent || right.divergent;
  result.leapfrog = left.leapfrog + right.leapfrog;
  return result;
}

HmcTransition nuts_transition(const Vector& q,
                              const HmcEvaluation& evaluated,
                              double step_size, const HmcMetric& metric,
                              int max_depth, HmcTarget& target,
                              double divergence_threshold,
                              HmcRandom& random) {
  Vector momentum = hmc_momentum(metric, random);
  const double joint0 = evaluated.logp - hmc_kinetic(momentum, metric);
  const double log_slice = joint0 - random.exponential();
  Vector q_minus = q, q_plus = q, proposal_q = q;
  Vector r_minus = momentum, r_plus = momentum;
  HmcEvaluation e_minus = evaluated, e_plus = evaluated;
  HmcEvaluation proposal = evaluated;
  int n = 1;
  bool active = true;
  double alpha = 0.0;
  int n_alpha = 0;
  bool divergent = false;
  int leapfrog = 0;
  int depth_reached = 0;
  for (int depth = 0; depth < max_depth && active; ++depth) {
    const int direction = random.direction();
    NutsTree tree = direction < 0 ?
      nuts_tree(
        q_minus, r_minus, e_minus, log_slice, direction, depth,
        step_size, metric, target, joint0, divergence_threshold, random) :
      nuts_tree(
        q_plus, r_plus, e_plus, log_slice, direction, depth,
        step_size, metric, target, joint0, divergence_threshold, random);
    if (tree.active && tree.n > 0 &&
        random.bernoulli(
          static_cast<double>(tree.n) /
          static_cast<double>(std::max(n + tree.n, 1)))) {
      proposal_q = tree.q_proposal;
      proposal = tree.e_proposal;
    }
    if (direction < 0) {
      q_minus = tree.q_minus;
      r_minus = tree.r_minus;
      e_minus = tree.e_minus;
    } else {
      q_plus = tree.q_plus;
      r_plus = tree.r_plus;
      e_plus = tree.e_plus;
    }
    n += tree.n;
    active = tree.active &&
      nuts_stop(q_minus, q_plus, r_minus, r_plus, metric);
    alpha += tree.alpha;
    n_alpha += tree.n_alpha;
    divergent = divergent || tree.divergent;
    leapfrog += tree.leapfrog;
    depth_reached = depth + 1;
  }
  HmcTransition result;
  result.q = proposal_q;
  result.evaluated = proposal;
  result.acceptance = alpha / static_cast<double>(std::max(n_alpha, 1));
  result.accepted = (proposal_q - q).squaredNorm() > 0.0;
  result.divergence = divergent;
  result.energy = -joint0;
  result.tree_depth = depth_reached;
  result.leapfrog = leapfrog;
  result.max_depth_reached = depth_reached >= max_depth && active;
  return result;
}

inline double hmc_log_add_exp(double left, double right) {
  if (!std::isfinite(left)) return right;
  if (!std::isfinite(right)) return left;
  const double maximum = std::max(left, right);
  return maximum + std::log(
    std::exp(left - maximum) + std::exp(right - maximum));
}

inline bool nuts_generalized_stop(const Vector& r_minus,
                                  const Vector& r_plus,
                                  const Vector& rho,
                                  const HmcMetric& metric) {
  return hmc_velocity(r_minus, metric).dot(rho) >= 0.0 &&
    hmc_velocity(r_plus, metric).dot(rho) >= 0.0;
}

struct MultinomialNutsTree {
  Vector q_minus, r_minus, q_plus, r_plus, q_proposal, rho;
  HmcEvaluation e_minus, e_plus, e_proposal;
  double log_weight = -std::numeric_limits<double>::infinity();
  bool active = false;
  double alpha = 0.0;
  int n_alpha = 0;
  bool divergent = false;
  int leapfrog = 0;
};

MultinomialNutsTree nuts_multinomial_tree(
    const Vector& q, const Vector& r, const HmcEvaluation& evaluated,
    int direction, int depth, double step_size, const HmcMetric& metric,
    HmcTarget& target, double joint0, double divergence_threshold,
    HmcRandom& random) {
  if (depth == 0) {
    HmcMove moved = hmc_leapfrog(
      q, r, evaluated.gradient, direction * step_size, metric, target);
    MultinomialNutsTree result;
    result.q_minus = result.q_plus = moved.q;
    result.r_minus = result.r_plus = moved.momentum;
    result.e_minus = result.e_plus = moved.evaluated;
    result.rho = moved.valid ? moved.momentum : Vector::Zero(r.size());
    result.n_alpha = 1;
    result.leapfrog = 1;
    const double joint = moved.valid ?
      moved.evaluated.logp - hmc_kinetic(moved.momentum, metric) :
      -std::numeric_limits<double>::infinity();
    const double error = joint0 - joint;
    const bool valid = moved.valid && std::isfinite(error) &&
      std::abs(error) <= divergence_threshold;
    result.q_proposal = valid ? moved.q : q;
    result.e_proposal = valid ? moved.evaluated : evaluated;
    result.log_weight = valid ? joint - joint0 :
      -std::numeric_limits<double>::infinity();
    result.active = valid;
    result.alpha = moved.valid ?
      std::min(1.0, std::exp(std::min(0.0, joint - joint0))) : 0.0;
    result.divergent = !valid;
    return result;
  }
  MultinomialNutsTree left = nuts_multinomial_tree(
    q, r, evaluated, direction, depth - 1, step_size, metric, target,
    joint0, divergence_threshold, random);
  if (!left.active) return left;
  MultinomialNutsTree right = direction < 0 ?
    nuts_multinomial_tree(
      left.q_minus, left.r_minus, left.e_minus, direction, depth - 1,
      step_size, metric, target, joint0, divergence_threshold, random) :
    nuts_multinomial_tree(
      left.q_plus, left.r_plus, left.e_plus, direction, depth - 1,
      step_size, metric, target, joint0, divergence_threshold, random);
  // A subtree that already terminated is not a reversible extension of its
  // parent trajectory. Its leapfrog/acceptance diagnostics still inform
  // adaptation, but its states must not enter the multinomial candidate set.
  if (!right.active) {
    left.active = false;
    left.alpha += right.alpha;
    left.n_alpha += right.n_alpha;
    left.divergent = left.divergent || right.divergent;
    left.leapfrog += right.leapfrog;
    return left;
  }
  MultinomialNutsTree result;
  result.log_weight = hmc_log_add_exp(left.log_weight, right.log_weight);
  const double right_probability = std::isfinite(right.log_weight) ?
    std::exp(right.log_weight - result.log_weight) : 0.0;
  const bool choose_right = random.bernoulli(right_probability);
  result.q_proposal = choose_right ? right.q_proposal : left.q_proposal;
  result.e_proposal = choose_right ? right.e_proposal : left.e_proposal;
  result.q_minus = direction < 0 ? right.q_minus : left.q_minus;
  result.r_minus = direction < 0 ? right.r_minus : left.r_minus;
  result.e_minus = direction < 0 ? right.e_minus : left.e_minus;
  result.q_plus = direction < 0 ? left.q_plus : right.q_plus;
  result.r_plus = direction < 0 ? left.r_plus : right.r_plus;
  result.e_plus = direction < 0 ? left.e_plus : right.e_plus;
  result.rho = left.rho + right.rho;
  const bool generalized_turn = nuts_generalized_stop(
    result.r_minus, result.r_plus, result.rho, metric);
  result.active = right.active && generalized_turn;
  result.alpha = left.alpha + right.alpha;
  result.n_alpha = left.n_alpha + right.n_alpha;
  result.divergent = left.divergent || right.divergent;
  result.leapfrog = left.leapfrog + right.leapfrog;
  return result;
}

HmcTransition nuts_multinomial_transition(
    const Vector& q, const HmcEvaluation& evaluated, double step_size,
    const HmcMetric& metric, int max_depth, HmcTarget& target,
    double divergence_threshold, HmcRandom& random) {
  Vector momentum = hmc_momentum(metric, random);
  const double joint0 = evaluated.logp - hmc_kinetic(momentum, metric);
  Vector q_minus = q, q_plus = q, proposal_q = q;
  Vector r_minus = momentum, r_plus = momentum, rho = momentum;
  HmcEvaluation e_minus = evaluated, e_plus = evaluated, proposal = evaluated;
  double log_weight = 0.0;
  bool active = true;
  double alpha = 0.0;
  int n_alpha = 0, leapfrog = 0, depth_reached = 0;
  bool divergent = false;
  for (int depth = 0; depth < max_depth && active; ++depth) {
    const int direction = random.direction();
    MultinomialNutsTree tree = direction < 0 ?
      nuts_multinomial_tree(
        q_minus, r_minus, e_minus, direction, depth, step_size, metric,
        target, joint0, divergence_threshold, random) :
      nuts_multinomial_tree(
        q_plus, r_plus, e_plus, direction, depth, step_size, metric,
        target, joint0, divergence_threshold, random);
    alpha += tree.alpha;
    n_alpha += tree.n_alpha;
    divergent = divergent || tree.divergent;
    leapfrog += tree.leapfrog;
    depth_reached = depth + 1;
    if (!tree.active) {
      active = false;
      break;
    }
    const double total_weight = hmc_log_add_exp(log_weight, tree.log_weight);
    const double tree_probability = std::isfinite(tree.log_weight) ?
      std::exp(tree.log_weight - total_weight) : 0.0;
    if (random.bernoulli(tree_probability)) {
      proposal_q = tree.q_proposal;
      proposal = tree.e_proposal;
    }
    if (direction < 0) {
      q_minus = tree.q_minus; r_minus = tree.r_minus; e_minus = tree.e_minus;
    } else {
      q_plus = tree.q_plus; r_plus = tree.r_plus; e_plus = tree.e_plus;
    }
    rho += tree.rho;
    log_weight = total_weight;
    active = tree.active && nuts_generalized_stop(
      r_minus, r_plus, rho, metric);
  }
  HmcTransition result;
  result.q = proposal_q;
  result.evaluated = proposal;
  result.acceptance = alpha / static_cast<double>(std::max(n_alpha, 1));
  result.accepted = (proposal_q - q).squaredNorm() > 0.0;
  result.divergence = divergent;
  result.energy = -joint0;
  result.tree_depth = depth_reached;
  result.leapfrog = leapfrog;
  result.max_depth_reached = depth_reached >= max_depth && active;
  return result;
}

HmcMetric adapted_metric(const std::vector<Vector>& warmup, int start, int count,
                         const std::string& kind,
                         const std::vector<std::vector<int>>& configured_blocks) {
  const int dimension = static_cast<int>(warmup.front().size());
  if (count - start < 2 || kind == "unit") {
    return HmcMetric(dimension, kind, configured_blocks);
  }
  Vector mean = Vector::Zero(warmup.front().size());
  for (int index = start; index < count; ++index) mean += warmup[index];
  mean /= static_cast<double>(count - start);
  std::vector<std::vector<int>> blocks;
  if (kind == "diagonal") {
    for (int index = 0; index < dimension; ++index) blocks.push_back({index});
  } else if (kind == "dense") {
    blocks.push_back(std::vector<int>(static_cast<std::size_t>(dimension)));
    std::iota(blocks.front().begin(), blocks.front().end(), 0);
  } else {
    blocks = configured_blocks;
  }
  const double shrink = static_cast<double>(count - start) /
    static_cast<double>(count - start + 5);
  Matrix mass = Matrix::Zero(dimension, dimension);
  for (const std::vector<int>& block : blocks) {
    const int width = static_cast<int>(block.size());
    Matrix covariance = Matrix::Zero(width, width);
    for (int sample = start; sample < count; ++sample) {
      for (int row = 0; row < width; ++row) {
        const double centered_row = warmup[sample][block[row]] - mean[block[row]];
        for (int column = 0; column < width; ++column) {
          covariance(row, column) += centered_row *
            (warmup[sample][block[column]] - mean[block[column]]);
        }
      }
    }
    covariance /= static_cast<double>(count - start - 1);
    covariance = shrink * covariance +
      (1.0 - shrink) * 1e-3 * Matrix::Identity(width, width);
    Eigen::SelfAdjointEigenSolver<Matrix> decomposition(
      0.5 * (covariance + covariance.transpose()));
    if (decomposition.info() != Eigen::Success) {
      throw std::runtime_error("HMC covariance adaptation failed.");
    }
    const Vector variance =
      decomposition.eigenvalues().array().max(1e-3).min(1e3);
    const Matrix precision = decomposition.eigenvectors() *
      variance.cwiseInverse().asDiagonal() *
      decomposition.eigenvectors().transpose();
    for (int row = 0; row < width; ++row) {
      for (int column = 0; column < width; ++column) {
        mass(block[row], block[column]) = precision(row, column);
      }
    }
  }
  return hmc_metric_from_mass(mass, kind, configured_blocks);
}

struct HmcWarmupSchedule {
  int initial_buffer = 0;
  int terminal_buffer = 0;
  int first_window = 25;
  std::vector<int> ends;
};

HmcWarmupSchedule hmc_warmup_schedule(int n_warmup, int initial_buffer,
                                      int terminal_buffer, int first_window) {
  HmcWarmupSchedule result;
  result.initial_buffer = initial_buffer;
  result.terminal_buffer = terminal_buffer;
  result.first_window = first_window;
  if (n_warmup < 20) {
    result.initial_buffer = n_warmup;
    result.terminal_buffer = 0;
    return result;
  }
  if (initial_buffer + terminal_buffer + first_window > n_warmup) {
    result.initial_buffer = std::max(1, static_cast<int>(std::floor(0.15 * n_warmup)));
    result.terminal_buffer = std::max(1, static_cast<int>(std::floor(0.10 * n_warmup)));
    result.first_window = std::max(
      1, n_warmup - result.initial_buffer - result.terminal_buffer);
  }
  const int slow = n_warmup - result.initial_buffer - result.terminal_buffer;
  if (slow <= 0) {
    result.initial_buffer = n_warmup;
    result.terminal_buffer = 0;
    return result;
  }
  int consumed = 0;
  int width = std::min(result.first_window, slow);
  while (consumed < slow) {
    const int remaining = slow - consumed;
    if (remaining < 2 * width) width = remaining;
    consumed += width;
    result.ends.push_back(result.initial_buffer + consumed);
    width = std::min(2 * width, slow - consumed);
    if (width <= 0) break;
  }
  return result;
}

Rcpp::List hmc_evaluation_to_r(const HmcEvaluation& evaluated) {
  return Rcpp::List::create(
    Rcpp::Named("logp") = evaluated.logp,
    Rcpp::Named("gradient") = libertad::eigen_vector_to_r(evaluated.gradient),
    Rcpp::Named("theta") = evaluated.parameters.theta,
    Rcpp::Named("sigma") = evaluated.parameters.sigma,
    Rcpp::Named("omega") = evaluated.parameters.omega,
    Rcpp::Named("eta") = libertad::eigen_vector_to_r(evaluated.eta),
    Rcpp::Named("outer") = evaluated.outer,
    Rcpp::Named("finite") = evaluated.finite
  );
}

Rcpp::List native_hmc_target_eval(ObjectiveTape& tape,
                                  const Rcpp::NumericVector& q,
                                  const Rcpp::List& config) {
  HmcTarget target(tape, config);
  const auto mapped = libertad::r_vector_map(q);
  HmcEvaluation evaluated = target.evaluate(mapped);
  return hmc_evaluation_to_r(evaluated);
}

Rcpp::List native_hmc_sample(ObjectiveTape& tape,
                             const Rcpp::List& config,
                             const std::string& method,
                             int n_warmup, int n_sample, int n_thin,
                             int n_chains, std::uint64_t seed,
                             double step_size, double target_acceptance,
                             bool adapt_mass, int n_leapfrog, int max_depth,
                             double divergence_threshold, int print_every) {
  if (method != "HMC" && method != "NUTS") {
    throw std::invalid_argument("Native sampler method must be HMC or NUTS.");
  }
  if (n_warmup < 0 || n_sample < 1 || n_thin < 1 || n_chains < 1 ||
      n_leapfrog < 1 || max_depth < 1 ||
      !(target_acceptance > 0.0 && target_acceptance < 1.0) ||
      !(divergence_threshold > 0.0)) {
    throw std::invalid_argument("Native HMC controls are invalid.");
  }
  HmcTarget target(tape, config);
  const std::string nuts_variant = config.containsElementNamed("nuts_variant") ?
    Rcpp::as<std::string>(config["nuts_variant"]) : "classic_slice";
  if (method == "NUTS" && nuts_variant != "classic_slice" &&
      nuts_variant != "multinomial") {
    throw std::invalid_argument("Native NUTS variant is invalid.");
  }
  const int dimension = static_cast<int>(target.dimension());
  const std::string metric_kind = config.containsElementNamed("hmc_metric") ?
    Rcpp::as<std::string>(config["hmc_metric"]) : "diagonal";
  if (metric_kind != "unit" && metric_kind != "diagonal" &&
      metric_kind != "dense" && metric_kind != "block") {
    throw std::invalid_argument("Native HMC metric is invalid.");
  }
  std::vector<std::vector<int>> metric_blocks;
  if (config.containsElementNamed("hmc_metric_blocks")) {
    Rcpp::List supplied = config["hmc_metric_blocks"];
    std::vector<bool> used(static_cast<std::size_t>(dimension), false);
    for (R_xlen_t block_index = 0; block_index < supplied.size(); ++block_index) {
      Rcpp::IntegerVector source = supplied[block_index];
      std::vector<int> block;
      block.reserve(static_cast<std::size_t>(source.size()));
      for (int index : source) {
        if (index < 0 || index >= dimension || used[static_cast<std::size_t>(index)]) {
          throw std::invalid_argument("Native HMC metric blocks are invalid or overlap.");
        }
        used[static_cast<std::size_t>(index)] = true;
        block.push_back(index);
      }
      if (block.empty()) throw std::invalid_argument("Native HMC metric blocks cannot be empty.");
      metric_blocks.push_back(std::move(block));
    }
    for (int index = 0; index < dimension; ++index) {
      if (!used[static_cast<std::size_t>(index)]) metric_blocks.push_back({index});
    }
  }
  const int initial_buffer = config.containsElementNamed("adapt_initial_buffer") ?
    Rcpp::as<int>(config["adapt_initial_buffer"]) : 75;
  const int terminal_buffer = config.containsElementNamed("adapt_terminal_buffer") ?
    Rcpp::as<int>(config["adapt_terminal_buffer"]) : 50;
  const int first_window = config.containsElementNamed("adapt_first_window") ?
    Rcpp::as<int>(config["adapt_first_window"]) : 25;
  const HmcWarmupSchedule warmup_schedule = hmc_warmup_schedule(
    n_warmup, initial_buffer, terminal_buffer, first_window);
  const std::string initialization = config.containsElementNamed("initialization") ?
    Rcpp::as<std::string>(config["initialization"]) : "auto";
  const double initialization_radius =
    config.containsElementNamed("initialization_radius") ?
      Rcpp::as<double>(config["initialization_radius"]) : 2.0;
  const int initialization_attempts =
    config.containsElementNamed("initialization_attempts") ?
      Rcpp::as<int>(config["initialization_attempts"]) : 100;
  if ((initialization != "auto" && initialization != "model" &&
       initialization != "jitter") || initialization_radius < 0.0 ||
      initialization_attempts < 1) {
    throw std::invalid_argument("Native HMC initialization controls are invalid.");
  }
  const int total = n_warmup + n_sample * n_thin;
  const int output_columns = static_cast<int>(
    config.containsElementNamed("output_columns") ?
      Rcpp::as<int>(config["output_columns"]) : dimension + 1);
  Rcpp::List chains(n_chains);
  Rcpp::List diagnostics(n_chains);
  const long long evaluations_before = target.evaluations();

  for (int chain_id = 0; chain_id < n_chains; ++chain_id) {
    HmcRandom random(seed + static_cast<std::uint64_t>(chain_id));
    Vector q(dimension);
    for (int index = 0; index < dimension; ++index) {
      q[index] = target.initial()[static_cast<std::size_t>(index)];
    }
    HmcEvaluation evaluated;
    bool initialized = false;
    int initialization_count = 0;
    std::string initialization_source;
    if (initialization == "auto" || initialization == "model") {
      evaluated = target.evaluate(q);
      initialization_count = 1;
      initialized = evaluated.finite;
      if (initialized) initialization_source = "model initial";
      if (!initialized && initialization == "model") {
        throw std::runtime_error(
          "The native HMC model initial point has non-finite posterior density or gradient.");
      }
    }
    for (int attempt = 0; !initialized && attempt < initialization_attempts; ++attempt) {
      const double radius = initialization_radius *
        std::pow(0.5, static_cast<double>(attempt / 25));
      for (int index = 0; index < dimension; ++index) {
        q[index] = target.initial()[static_cast<std::size_t>(index)] +
          radius * (2.0 * random.uniform() - 1.0);
      }
      evaluated = target.evaluate(q);
      ++initialization_count;
      initialized = evaluated.finite;
      if (initialized) initialization_source = "bounded random jitter";
    }
    if (!initialized) {
      throw std::runtime_error(
        "Unable to initialize native HMC at a finite posterior density after bounded attempts.");
    }

    HmcMetric metric(dimension, metric_kind, metric_blocks);
    double epsilon = std::isfinite(step_size) && step_size > 0.0 ?
      step_size : hmc_find_step(q, evaluated, metric, target, random);
    DualAverage dual(epsilon, target_acceptance);
    std::vector<Vector> warmup(
      static_cast<std::size_t>(std::max(n_warmup, 1)),
      Vector::Zero(dimension));
    std::size_t metric_window = 0;
    int metric_window_start = warmup_schedule.ends.empty() ? 0 :
      warmup_schedule.initial_buffer;
    Rcpp::NumericMatrix draws(n_sample, output_columns);
    Rcpp::NumericMatrix trace(total, 6);
    Rcpp::CharacterVector trace_names = Rcpp::CharacterVector::create(
      "acceptance", "divergence", "tree_depth", "leapfrog", "step_size",
      "energy");
    trace.attr("dimnames") =
      Rcpp::List::create(R_NilValue, trace_names);
    int keep = 0;
    int post_divergences = 0;
    int max_depth_hits = 0;
    double post_acceptance = 0.0;
    int post_iterations = 0;

    for (int iteration = 0; iteration < total; ++iteration) {
      HmcTransition transition = method == "NUTS" ?
        (nuts_variant == "multinomial" ?
          nuts_multinomial_transition(
            q, evaluated, epsilon, metric, max_depth, target,
            divergence_threshold, random) :
          nuts_transition(
            q, evaluated, epsilon, metric, max_depth, target,
            divergence_threshold, random)) :
        hmc_transition(
          q, evaluated, epsilon, metric, n_leapfrog, target,
          divergence_threshold, random);
      q = std::move(transition.q);
      evaluated = std::move(transition.evaluated);
      trace(iteration, 0) = transition.acceptance;
      trace(iteration, 1) = transition.divergence ? 1.0 : 0.0;
      trace(iteration, 2) =
        transition.tree_depth < 0 ? NA_REAL : transition.tree_depth;
      trace(iteration, 3) = transition.leapfrog;
      trace(iteration, 4) = epsilon;
      trace(iteration, 5) = transition.energy;

      if (iteration < n_warmup) {
        warmup[static_cast<std::size_t>(iteration)] = q;
        epsilon = dual.update(transition.acceptance);
        if (adapt_mass && metric_kind != "unit" &&
            metric_window < warmup_schedule.ends.size() &&
            iteration + 1 == warmup_schedule.ends[metric_window]) {
          metric = adapted_metric(
            warmup, metric_window_start, iteration + 1,
            metric_kind, metric_blocks);
          epsilon = hmc_find_step(q, evaluated, metric, target, random);
          dual = DualAverage(epsilon, target_acceptance);
          metric_window_start = iteration + 1;
          ++metric_window;
        }
        if (iteration + 1 == n_warmup) epsilon = dual.final();
      } else {
        ++post_iterations;
        post_acceptance += transition.acceptance;
        if (transition.divergence) ++post_divergences;
        if (transition.max_depth_reached) ++max_depth_hits;
        if ((iteration + 1 - n_warmup) % n_thin == 0) {
          Rcpp::NumericVector row = target.native_row(evaluated);
          if (row.size() != output_columns) {
            throw std::logic_error("Native HMC output layout is inconsistent.");
          }
          for (int column = 0; column < output_columns; ++column) {
            draws(keep, column) = row[column];
          }
          ++keep;
        }
      }

      if (print_every > 0 && (iteration + 1) % print_every == 0) {
        Rcpp::Rcout << "[LibeRation] " << method << " CHAIN "
                    << chain_id + 1 << " ITERATION " << iteration + 1
                    << " LOGPOST " << evaluated.logp
                    << " ACCEPT " << transition.acceptance
                    << " STEP " << epsilon
                    << " DIVERGENT "
                    << (transition.divergence ? "TRUE" : "FALSE") << "\n";
      }
      Rcpp::checkUserInterrupt();
    }
    chains[chain_id] = draws;
    double ebfmi = NA_REAL;
    if (post_iterations >= 3) {
      double mean_energy = 0.0;
      for (int index = n_warmup; index < total; ++index) {
        mean_energy += trace(index, 5);
      }
      mean_energy /= static_cast<double>(post_iterations);
      double numerator = 0.0;
      double denominator = 0.0;
      for (int index = n_warmup; index < total; ++index) {
        const double centered = trace(index, 5) - mean_energy;
        denominator += centered * centered;
        if (index > n_warmup) {
          const double delta = trace(index, 5) - trace(index - 1, 5);
          numerator += delta * delta;
        }
      }
      if (std::isfinite(numerator) && denominator > 0.0) {
        ebfmi = numerator / denominator;
      }
    }
    Rcpp::RObject mass_output = metric.kind == "diagonal" ?
      Rcpp::RObject(libertad::eigen_vector_to_r(metric.mass.diagonal())) :
      Rcpp::RObject(libertad::eigen_matrix_to_r(metric.mass));
    Rcpp::List schedule_output = Rcpp::List::create(
      Rcpp::Named("n_warmup") = n_warmup,
      Rcpp::Named("initial_buffer") = warmup_schedule.initial_buffer,
      Rcpp::Named("terminal_buffer") = warmup_schedule.terminal_buffer,
      Rcpp::Named("first_window") = warmup_schedule.first_window,
      Rcpp::Named("windows") = Rcpp::wrap(warmup_schedule.ends));
    diagnostics[chain_id] = Rcpp::List::create(
      Rcpp::Named("trace") = trace,
      Rcpp::Named("step_size") = epsilon,
      Rcpp::Named("mass") = mass_output,
      Rcpp::Named("metric") = metric.kind,
      Rcpp::Named("warmup_schedule") = schedule_output,
      Rcpp::Named("initialization") = Rcpp::List::create(
        Rcpp::Named("source") = initialization_source,
        Rcpp::Named("attempts") = initialization_count),
      Rcpp::Named("ebfmi") = ebfmi,
      Rcpp::Named("divergences") = post_divergences,
      Rcpp::Named("mean_acceptance") =
        post_iterations ? post_acceptance / post_iterations : NA_REAL,
      Rcpp::Named("max_depth_hits") = max_depth_hits
    );
  }
  return Rcpp::List::create(
    Rcpp::Named("chains") = chains,
    Rcpp::Named("diagnostics") = diagnostics,
    Rcpp::Named("objective_evaluations") =
      static_cast<double>(target.evaluations() - evaluations_before),
    Rcpp::Named("backend") = "native-cpp-cppad"
  );
}
