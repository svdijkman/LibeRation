// Internal implementation fragment included inside namespace liberation by
// pk_engine.cpp. It intentionally has no include guard or namespace declaration.
// Contains shared stochastic parameter maps, MU configuration, subject-parallel support, and the native subject pool.

// Native fixed-ETA SAEM M-step.  This file is included inside namespace
// liberation by pk_engine.cpp so that it can reuse ObjectiveTape without
// exposing implementation-only types in installed headers.

struct SaemPrior {
  int native_index = -1;
  std::string family;
  double mean = 0.0;
  double sd = 1.0;
  double shape = std::numeric_limits<double>::quiet_NaN();
  double rate = std::numeric_limits<double>::quiet_NaN();
};

struct SaemEvaluation {
  double value = 1e100;
  Vector gradient;
  // Gradient in the complete objective-tape domain
  // THETA, ETA, SIGMA, OMEGA.  The optimizer uses the encoded subset above;
  // retaining this vector lets exact sufficient-statistic updates reuse the
  // final Reverse(1) sweep instead of evaluating the full weighted Q again.
  Vector native_gradient;
};

struct StochasticBayesParameters {
  std::vector<double> theta;
  std::vector<double> sigma;
  std::vector<double> omega;
  double log_jacobian = 0.0;
};

// Parameter-map subset needed by the optimized random-walk BAYES sampler.
// This mirrors .nm_outer_map() exactly but intentionally omits derivatives:
// the Metropolis coordinator requires only decoded native values, bounds,
// priors, and the transformation Jacobian.
class StochasticBayesMap {
 public:
  explicit StochasticBayesMap(const Rcpp::List& config) {
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
    start_ = Rcpp::as<std::vector<double>>(config["start"]);
    lower_ = Rcpp::as<std::vector<double>>(config["lower"]);
    upper_ = Rcpp::as<std::vector<double>>(config["upper"]);
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
    if (prior_family.size() != prior_index.size() ||
        prior_mean.size() != prior_index.size() ||
        prior_sd.size() != prior_index.size() ||
        prior_shape.size() != prior_index.size() ||
        prior_rate.size() != prior_index.size()) {
      throw std::invalid_argument("Native BAYES prior mapping is inconsistent.");
    }
    for (std::size_t index = 0; index < prior_index.size(); ++index) {
      priors_.push_back(PopulationPrior{
        prior_index[index], prior_family[index], prior_mean[index],
        prior_sd[index], prior_shape[index], prior_rate[index]
      });
    }
    const std::size_t expected = theta_free_.size() + sigma_free_.size() +
      (omega_full_ && !omega_free_.empty() ? omega_base_.size() :
       omega_free_.size());
    if (start_.size() != expected || lower_.size() != expected ||
        upper_.size() != expected || n_eta_base_ < 0 ||
        omega_rows_.size() != omega_base_.size() ||
        omega_cols_.size() != omega_base_.size()) {
      throw std::invalid_argument("Native BAYES parameter map is inconsistent.");
    }
  }

  std::size_t dimension() const { return start_.size(); }
  const std::vector<double>& start() const { return start_; }
  const std::vector<double>& lower() const { return lower_; }
  const std::vector<double>& upper() const { return upper_; }

  bool in_bounds(const Vector& encoded) const {
    if (static_cast<std::size_t>(encoded.size()) != dimension() ||
        !encoded.allFinite()) return false;
    for (Eigen::Index index = 0; index < encoded.size(); ++index) {
      if (encoded[index] < lower_[static_cast<std::size_t>(index)] ||
          encoded[index] > upper_[static_cast<std::size_t>(index)]) return false;
    }
    return true;
  }

  StochasticBayesParameters decode(const Vector& encoded) const {
    if (!in_bounds(encoded)) {
      throw std::domain_error("Native BAYES proposal is outside its bounds.");
    }
    StochasticBayesParameters result;
    result.theta = theta_base_;
    result.sigma = sigma_base_;
    result.omega = omega_base_;
    std::size_t cursor = 0U;
    for (int index : theta_free_) {
      result.theta[static_cast<std::size_t>(index)] =
        encoded[static_cast<Eigen::Index>(cursor++)];
    }
    for (int index : sigma_free_) {
      const double value = std::exp(
        encoded[static_cast<Eigen::Index>(cursor++)]);
      result.sigma[static_cast<std::size_t>(index)] = value;
      result.log_jacobian += std::log(value);
    }
    if (omega_full_ && !omega_free_.empty()) {
      Matrix lower = Matrix::Zero(n_eta_base_, n_eta_base_);
      for (std::size_t entry = 0; entry < omega_base_.size(); ++entry) {
        const int row = omega_rows_[entry];
        const int column = omega_cols_[entry];
        if (row < 0 || column < 0 || row >= n_eta_base_ || column > row) {
          throw std::invalid_argument("Native BAYES OMEGA coordinates are invalid.");
        }
        const double value = encoded[
          static_cast<Eigen::Index>(cursor + entry)];
        lower(row, column) = row == column ? std::exp(value) : value;
      }
      const Matrix covariance = lower * lower.transpose();
      for (std::size_t entry = 0; entry < omega_base_.size(); ++entry) {
        result.omega[entry] = covariance(
          omega_rows_[entry], omega_cols_[entry]);
      }
      result.log_jacobian += static_cast<double>(n_eta_base_) * std::log(2.0);
      for (int row = 0; row < n_eta_base_; ++row) {
        result.log_jacobian += static_cast<double>(n_eta_base_ + 1 - row) *
          std::log(lower(row, row));
      }
      cursor += omega_base_.size();
    } else {
      for (int index : omega_free_) {
        const double value = std::exp(
          encoded[static_cast<Eigen::Index>(cursor++)]);
        result.omega[static_cast<std::size_t>(index)] = value;
        result.log_jacobian += std::log(value);
      }
    }
    if (cursor != dimension() || !std::isfinite(result.log_jacobian)) {
      throw std::invalid_argument("Native BAYES decoding failed.");
    }
    return result;
  }

  Vector encode(const StochasticBayesParameters& parameters) const {
    Vector encoded(static_cast<Eigen::Index>(dimension()));
    std::size_t cursor = 0U;
    for (int index : theta_free_) {
      encoded[static_cast<Eigen::Index>(cursor++)] =
        parameters.theta[static_cast<std::size_t>(index)];
    }
    for (int index : sigma_free_) {
      const double value = parameters.sigma[static_cast<std::size_t>(index)];
      if (!(value > 0.0)) {
        throw std::domain_error("Native BAYES SIGMA encoding is invalid.");
      }
      encoded[static_cast<Eigen::Index>(cursor++)] = std::log(value);
    }
    if (omega_full_ && !omega_free_.empty()) {
      Matrix covariance = omega_covariance(parameters);
      Eigen::LLT<Matrix> factor(covariance);
      if (factor.info() != Eigen::Success) {
        throw std::domain_error("Native BAYES OMEGA encoding is invalid.");
      }
      const Matrix lower = Matrix(factor.matrixL());
      for (std::size_t entry = 0; entry < omega_base_.size(); ++entry) {
        const int row = omega_rows_[entry];
        const int column = omega_cols_[entry];
        encoded[static_cast<Eigen::Index>(cursor++)] = row == column ?
          std::log(lower(row, column)) : lower(row, column);
      }
    } else {
      for (int index : omega_free_) {
        const double value = parameters.omega[static_cast<std::size_t>(index)];
        if (!(value > 0.0)) {
          throw std::domain_error("Native BAYES OMEGA encoding is invalid.");
        }
        encoded[static_cast<Eigen::Index>(cursor++)] = std::log(value);
      }
    }
    if (cursor != dimension() || !in_bounds(encoded)) {
      throw std::domain_error("Native BAYES encoded state is outside its bounds.");
    }
    return encoded;
  }

  int theta_outer_position(int native_index) const {
    const auto found = std::find(
      theta_free_.begin(), theta_free_.end(), native_index);
    return found == theta_free_.end() ? -1 :
      static_cast<int>(std::distance(theta_free_.begin(), found));
  }

  int omega_outer_position(int native_index) const {
    if (omega_full_) return -1;
    const auto found = std::find(
      omega_free_.begin(), omega_free_.end(), native_index);
    return found == omega_free_.end() ? -1 :
      static_cast<int>(theta_free_.size() + sigma_free_.size() +
        std::distance(omega_free_.begin(), found));
  }

  bool diagonal_omega_inverse_gamma(
      int native_index, double& shape, double& rate) const {
    if (omega_full_ || native_index < 0 ||
        native_index >= static_cast<int>(omega_rows_.size()) ||
        omega_rows_[static_cast<std::size_t>(native_index)] !=
          omega_cols_[static_cast<std::size_t>(native_index)] ||
        omega_outer_position(native_index) < 0) return false;
    const int prior_native = static_cast<int>(theta_base_.size() +
      sigma_base_.size()) + native_index;
    int matches = 0;
    for (const PopulationPrior& prior : priors_) {
      if (prior.native_index == prior_native) {
        if (prior.family != "inverse_gamma" || !(prior.shape > 0.0) ||
            !(prior.rate > 0.0)) return false;
        shape = prior.shape;
        rate = prior.rate;
        ++matches;
      }
    }
    return matches == 1;
  }

  int omega_effect(int native_index) const {
    return native_index >= 0 &&
      native_index < static_cast<int>(omega_rows_.size()) ?
      omega_rows_[static_cast<std::size_t>(native_index)] : -1;
  }

  double log_prior(const StochasticBayesParameters& parameters) const {
    std::vector<double> native;
    native.reserve(parameters.theta.size() + parameters.sigma.size() +
                   parameters.omega.size());
    native.insert(native.end(), parameters.theta.begin(), parameters.theta.end());
    native.insert(native.end(), parameters.sigma.begin(), parameters.sigma.end());
    native.insert(native.end(), parameters.omega.begin(), parameters.omega.end());
    const double log_two_pi = std::log(2.0 * std::acos(-1.0));
    double total = 0.0;
    for (const PopulationPrior& prior : priors_) {
      if (prior.native_index < 0 ||
          prior.native_index >= static_cast<int>(native.size())) {
        return -std::numeric_limits<double>::infinity();
      }
      const double value = native[static_cast<std::size_t>(prior.native_index)];
      double density = -std::numeric_limits<double>::infinity();
      if (prior.family == "normal" || prior.family == "half_normal") {
        if (prior.sd > 0.0 && std::isfinite(value) &&
            (prior.family != "half_normal" || value >= 0.0)) {
          const double z = (value - prior.mean) / prior.sd;
          density = -0.5 * log_two_pi - std::log(prior.sd) - 0.5 * z * z;
          if (prior.family == "half_normal") density += std::log(2.0);
        }
      } else if (prior.family == "lognormal") {
        if (value > 0.0 && prior.sd > 0.0) {
          const double z = (std::log(value) - prior.mean) / prior.sd;
          density = -std::log(value) - 0.5 * log_two_pi -
            std::log(prior.sd) - 0.5 * z * z;
        }
      } else if (prior.family == "inverse_gamma") {
        if (value > 0.0 && prior.shape > 0.0 && prior.rate > 0.0) {
          density = prior.shape * std::log(prior.rate) -
            std::lgamma(prior.shape) -
            (prior.shape + 1.0) * std::log(value) - prior.rate / value;
        }
      } else {
        throw std::invalid_argument("Unknown native BAYES prior family.");
      }
      if (!std::isfinite(density)) {
        return -std::numeric_limits<double>::infinity();
      }
      total += density;
    }
    return total;
  }

  // Return the same -2 log-prior contribution and native-parameter
  // derivative used by PopulationObjective.  Keeping this next to the map
  // makes stochastic marginal objectives independent of R callbacks without
  // changing their prior convention.
  double prior_nll(
      const StochasticBayesParameters& parameters,
      Vector* derivative = nullptr) const {
    std::vector<double> native;
    native.reserve(parameters.theta.size() + parameters.sigma.size() +
                   parameters.omega.size());
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
        throw std::invalid_argument("A native prior refers to an invalid parameter.");
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
          density = prior.shape * std::log(prior.rate) -
            std::lgamma(prior.shape) -
            (prior.shape + 1.0) * std::log(value) - prior.rate / value;
          gradient = 2.0 * (prior.shape + 1.0) / value -
            2.0 * prior.rate / (value * value);
        }
      } else {
        throw std::invalid_argument("Unknown native prior family.");
      }
      if (!std::isfinite(density) ||
          (derivative && !std::isfinite(gradient))) {
        return 1e100;
      }
      log_density += density;
      if (derivative) (*derivative)[prior.native_index] += gradient;
    }
    return -2.0 * log_density;
  }

  // Map a THETA/SIGMA/OMEGA derivative into the optimizer coordinates.  The
  // full-OMEGA branch differentiates Omega = L L' explicitly, matching
  // .nm_outer_map() rather than approximating the Cholesky chain rule.
  Vector outer_gradient(
      const Vector& encoded, const StochasticBayesParameters& parameters,
      const Vector& native_gradient) const {
    const std::size_t native_size = theta_base_.size() + sigma_base_.size() +
      omega_base_.size();
    if (static_cast<std::size_t>(native_gradient.size()) != native_size ||
        static_cast<std::size_t>(encoded.size()) != dimension()) {
      throw std::invalid_argument("Native stochastic gradient dimensions are inconsistent.");
    }
    Vector result = Vector::Zero(static_cast<Eigen::Index>(dimension()));
    std::size_t cursor = 0U;
    for (int index : theta_free_) {
      result[static_cast<Eigen::Index>(cursor++)] = native_gradient[index];
    }
    const int sigma_offset = static_cast<int>(theta_base_.size());
    for (int index : sigma_free_) {
      result[static_cast<Eigen::Index>(cursor++)] =
        native_gradient[sigma_offset + index] *
        parameters.sigma[static_cast<std::size_t>(index)];
    }
    const int omega_offset = sigma_offset + static_cast<int>(sigma_base_.size());
    if (omega_full_ && !omega_free_.empty()) {
      Matrix lower = Matrix::Zero(n_eta_base_, n_eta_base_);
      for (std::size_t entry = 0; entry < omega_base_.size(); ++entry) {
        const int row = omega_rows_[entry];
        const int column = omega_cols_[entry];
        const double value = encoded[static_cast<Eigen::Index>(cursor + entry)];
        lower(row, column) = row == column ? std::exp(value) : value;
      }
      for (std::size_t encoded_entry = 0; encoded_entry < omega_base_.size();
           ++encoded_entry) {
        Matrix derivative_lower = Matrix::Zero(n_eta_base_, n_eta_base_);
        const int row = omega_rows_[encoded_entry];
        const int column = omega_cols_[encoded_entry];
        derivative_lower(row, column) = row == column ?
          lower(row, column) : 1.0;
        const Matrix derivative = derivative_lower * lower.transpose() +
          lower * derivative_lower.transpose();
        double value = 0.0;
        for (std::size_t native = 0; native < omega_base_.size(); ++native) {
          value += native_gradient[omega_offset + static_cast<int>(native)] *
            derivative(omega_rows_[native], omega_cols_[native]);
        }
        result[static_cast<Eigen::Index>(cursor + encoded_entry)] = value;
      }
      cursor += omega_base_.size();
    } else {
      for (int index : omega_free_) {
        result[static_cast<Eigen::Index>(cursor++)] =
          native_gradient[omega_offset + index] *
          parameters.omega[static_cast<std::size_t>(index)];
      }
    }
    if (cursor != dimension() || !result.allFinite()) {
      throw std::runtime_error("Native stochastic gradient mapping failed.");
    }
    return result;
  }

  Matrix omega_covariance(
      const StochasticBayesParameters& parameters) const {
    Matrix covariance = Matrix::Zero(n_eta_base_, n_eta_base_);
    if (omega_rows_.size() != parameters.omega.size()) {
      throw std::invalid_argument("Native BAYES OMEGA dimensions changed.");
    }
    for (std::size_t entry = 0; entry < parameters.omega.size(); ++entry) {
      covariance(omega_rows_[entry], omega_cols_[entry]) =
        parameters.omega[entry];
      covariance(omega_cols_[entry], omega_rows_[entry]) =
        parameters.omega[entry];
    }
    return covariance;
  }

 private:
  std::vector<double> theta_base_, sigma_base_, omega_base_;
  std::vector<int> theta_free_, sigma_free_, omega_free_;
  std::vector<int> omega_rows_, omega_cols_;
  std::vector<double> start_, lower_, upper_;
  std::vector<PopulationPrior> priors_;
  bool omega_full_ = false;
  int n_eta_base_ = 0;

  static std::vector<int> zero_based(std::vector<int> source) {
    for (int& value : source) {
      if (value < 1) {
        throw std::invalid_argument("A native BAYES parameter index is invalid.");
      }
      --value;
    }
    return source;
  }
};

struct StochasticMuConfig {
  bool active = false;
  std::vector<int> theta;
  std::vector<std::string> links;
  std::vector<Matrix> design_columns;

  StochasticMuConfig() = default;

  StochasticMuConfig(
      const Rcpp::List& config, int subjects, int n_eta) {
    active = config.containsElementNamed("active") &&
      Rcpp::as<bool>(config["active"]);
    if (!active) return;
    theta = Rcpp::as<std::vector<int>>(config["theta"]);
    links = Rcpp::as<std::vector<std::string>>(config["links"]);
    const Rcpp::List source = config["design_columns"];
    if (theta.empty() || theta.size() != links.size() ||
        source.size() != static_cast<int>(theta.size())) {
      throw std::invalid_argument("Native BAYES MU configuration is inconsistent.");
    }
    design_columns.reserve(theta.size());
    for (std::size_t column = 0; column < theta.size(); ++column) {
      if (theta[column] < 1) {
        throw std::invalid_argument("A native BAYES MU THETA index is invalid.");
      }
      --theta[column];
      if (links[column] != "identity" && links[column] != "log") {
        throw std::invalid_argument("A native BAYES MU link is invalid.");
      }
      Rcpp::NumericMatrix input = source[static_cast<int>(column)];
      if (input.nrow() != subjects || input.ncol() != n_eta) {
        throw std::invalid_argument("A native BAYES MU design has invalid dimensions.");
      }
      Matrix design(subjects, n_eta);
      for (int row = 0; row < subjects; ++row) {
        for (int effect = 0; effect < n_eta; ++effect) {
          design(row, effect) = input(row, effect);
        }
      }
      design_columns.push_back(std::move(design));
    }
  }

  Vector beta(const StochasticBayesParameters& parameters) const {
    Vector result(static_cast<Eigen::Index>(theta.size()));
    for (std::size_t column = 0; column < theta.size(); ++column) {
      const double value = parameters.theta[static_cast<std::size_t>(theta[column])];
      if (links[column] == "log" && !(value > 0.0)) {
        throw std::domain_error("A log-linked MU THETA is not positive.");
      }
      result[static_cast<Eigen::Index>(column)] =
        links[column] == "log" ? std::log(value) : value;
    }
    return result;
  }

  void set_beta(StochasticBayesParameters& parameters, const Vector& value) const {
    for (std::size_t column = 0; column < theta.size(); ++column) {
      parameters.theta[static_cast<std::size_t>(theta[column])] =
        links[column] == "log" ?
          std::exp(value[static_cast<Eigen::Index>(column)]) :
          value[static_cast<Eigen::Index>(column)];
    }
  }

  Matrix recenter(
      const Matrix& eta, const Vector& old_beta,
      const Vector& new_beta) const {
    Matrix result = eta;
    for (std::size_t column = 0; column < design_columns.size(); ++column) {
      result += design_columns[column] *
        (old_beta[static_cast<Eigen::Index>(column)] -
         new_beta[static_cast<Eigen::Index>(column)]);
    }
    return result;
  }

  double log_native_jacobian(
      const StochasticBayesParameters& parameters) const {
    double result = 0.0;
    for (std::size_t column = 0; column < theta.size(); ++column) {
      if (links[column] == "log") {
        result -= std::log(
          parameters.theta[static_cast<std::size_t>(theta[column])]);
      }
    }
    return result;
  }
};

struct NativeGqEvaluation {
  double value = 1e100;
  Vector native_gradient;
  Matrix modes;
  std::vector<double> effective_points;
  std::vector<double> cancellation_ratio;
  bool valid = false;
  long long node_evaluations = 0;
};

// CppAD's allocator needs a stable thread identity whenever independent ADFun
// instances are evaluated concurrently.  Subject workers occupy slots 1..47;
// the R/main thread remains slot zero.  The callbacks are process-global, but
// parallel mode is enabled only for the bounded lifetime of a pool dispatch.
inline thread_local std::size_t cppad_subject_thread_number = 0U;
inline std::atomic<int> cppad_subject_parallel_dispatches{0};

inline bool cppad_subject_in_parallel() {
  return cppad_subject_parallel_dispatches.load(std::memory_order_acquire) > 0;
}

inline std::size_t cppad_subject_thread_num() {
  return cppad_subject_thread_number;
}

inline void configure_cppad_subject_parallelism() {
  static std::once_flag configured;
  std::call_once(configured, []() {
    CppAD::thread_alloc::parallel_setup(
      CPPAD_MAX_NUM_THREADS, cppad_subject_in_parallel,
      cppad_subject_thread_num);
    CppAD::thread_alloc::hold_memory(true);
  });
}

// A small persistent worker team for both direct-double stochastic kernels and
// independent subject CppAD tapes.  Work is statically partitioned and every
// result is reduced later on the main thread in subject order.  CppAD mode must
// only be used when no ADFun pointer is shared between subjects because Forward
// and Reverse mutate per-tape work state.
class NativeSubjectPool {
 public:
  explicit NativeSubjectPool(int workers)
      : workers_(std::max(1, workers)), errors_(static_cast<std::size_t>(workers_)) {
    threads_.reserve(static_cast<std::size_t>(workers_));
    for (int worker = 0; worker < workers_; ++worker) {
      threads_.emplace_back([this, worker]() { worker_loop(worker); });
    }
  }

  ~NativeSubjectPool() {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      stopping_ = true;
      ++generation_;
    }
    start_.notify_all();
    for (std::thread& worker : threads_) {
      if (worker.joinable()) worker.join();
    }
  }

  NativeSubjectPool(const NativeSubjectPool&) = delete;
  NativeSubjectPool& operator=(const NativeSubjectPool&) = delete;

  template <class Function>
  void run(std::size_t count, Function&& function) {
    run_impl(count, std::forward<Function>(function), false);
  }

  template <class Function>
  void run_cppad(std::size_t count, Function&& function) {
    configure_cppad_subject_parallelism();
    run_impl(count, std::forward<Function>(function), true);
  }

  long long dispatches() const { return dispatches_; }
  long long cppad_dispatches() const { return cppad_dispatches_; }
  int workers() const { return workers_; }

 private:
  template <class Function>
  void run_impl(std::size_t count, Function&& function, bool cppad) {
    if (!count) return;
    {
      std::lock_guard<std::mutex> lock(mutex_);
      if (active_ != 0) {
        throw std::logic_error("The native subject pool cannot be nested.");
      }
      count_ = count;
      task_ = std::forward<Function>(function);
      cppad_task_ = cppad;
      std::fill(errors_.begin(), errors_.end(), std::exception_ptr());
      active_ = workers_;
      ++generation_;
      ++dispatches_;
      if (cppad) {
        ++cppad_dispatches_;
        cppad_subject_parallel_dispatches.fetch_add(1, std::memory_order_release);
      }
    }
    start_.notify_all();
    {
      std::unique_lock<std::mutex> lock(mutex_);
      finished_.wait(lock, [this]() { return active_ == 0; });
    }
    if (cppad) {
      cppad_subject_parallel_dispatches.fetch_sub(1, std::memory_order_release);
    }
    for (const std::exception_ptr& error : errors_) {
      if (error) std::rethrow_exception(error);
    }
  }
  int workers_ = 1;
  std::vector<std::thread> threads_;
  mutable std::mutex mutex_;
  std::condition_variable start_;
  std::condition_variable finished_;
  bool stopping_ = false;
  std::size_t generation_ = 0U;
  std::size_t count_ = 0U;
  int active_ = 0;
  std::function<void(std::size_t)> task_;
  bool cppad_task_ = false;
  std::vector<std::exception_ptr> errors_;
  long long dispatches_ = 0;
  long long cppad_dispatches_ = 0;

  void worker_loop(int worker) {
    std::size_t observed_generation = 0U;
    for (;;) {
      std::function<void(std::size_t)> task;
      std::size_t count = 0U;
      bool cppad = false;
      {
        std::unique_lock<std::mutex> lock(mutex_);
        start_.wait(lock, [&]() {
          return stopping_ || generation_ != observed_generation;
        });
        if (stopping_) return;
        observed_generation = generation_;
        task = task_;
        count = count_;
        cppad = cppad_task_;
      }
      try {
        if (cppad) cppad_subject_thread_number =
          static_cast<std::size_t>(worker + 1);
        const std::size_t begin = count * static_cast<std::size_t>(worker) /
          static_cast<std::size_t>(workers_);
        const std::size_t end = count * static_cast<std::size_t>(worker + 1) /
          static_cast<std::size_t>(workers_);
        for (std::size_t subject = begin; subject < end; ++subject) {
          task(subject);
        }
      } catch (...) {
        errors_[static_cast<std::size_t>(worker)] = std::current_exception();
      }
      cppad_subject_thread_number = 0U;
      {
        std::lock_guard<std::mutex> lock(mutex_);
        --active_;
        if (active_ == 0) finished_.notify_one();
      }
    }
  }
};

// Persistent subject collection for stochastic estimators. Dynamic
// observations/covariates are installed once and parameter/ETA point buffers
// are reused across iterations. This is restricted to stable, non-retaping
// optimized contexts by the R coordinator.
