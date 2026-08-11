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
class SaemFixedEtaCollection {
 public:
  SaemFixedEtaCollection(
      const Rcpp::List& tape_pointers, const Rcpp::NumericMatrix& eta,
      int n_theta, int n_sigma, int n_omega)
      : n_theta_(n_theta), n_eta_(eta.ncol()), n_sigma_(n_sigma),
        n_omega_(n_omega) {
    if (tape_pointers.size() != eta.nrow() || tape_pointers.size() < 1 ||
        n_theta < 0 || n_sigma < 0 || n_omega < 0 || eta.ncol() < 0) {
      throw std::invalid_argument(
        "Persistent SAEM fixed-ETA inputs are inconsistent.");
    }
    domain_ = static_cast<std::size_t>(
      n_theta_ + n_eta_ + n_sigma_ + n_omega_);
    tapes_.reserve(static_cast<std::size_t>(tape_pointers.size()));
    points_.resize(static_cast<std::size_t>(tape_pointers.size()));
    for (int subject = 0; subject < tape_pointers.size(); ++subject) {
      Rcpp::XPtr<ObjectiveTape> tape(tape_pointers[subject]);
      if (tape->domain_names.size() != domain_) {
        throw std::invalid_argument(
          "A persistent SAEM objective tape has an inconsistent domain.");
      }
      tapes_.push_back(tape.get());
      std::vector<double>& point = points_[static_cast<std::size_t>(subject)];
      point.assign(domain_, 0.0);
      for (int effect = 0; effect < n_eta_; ++effect) {
        const double value = eta(subject, effect);
        if (!std::isfinite(value)) {
          throw std::invalid_argument("Persistent SAEM ETAs must be finite.");
        }
        point[static_cast<std::size_t>(n_theta_ + effect)] = value;
      }
    }
  }

  Rcpp::List evaluate(
      const Rcpp::NumericVector& theta, const Rcpp::NumericVector& sigma,
      const Rcpp::NumericVector& omega) {
    validate_parameters(theta, sigma, omega);
    Rcpp::NumericVector values(static_cast<R_xlen_t>(tapes_.size()));
    Rcpp::NumericMatrix gradients(
      static_cast<int>(tapes_.size()), static_cast<int>(domain_));
    const std::vector<double> weight(1U, 1.0);
    std::ostringstream messages;
    for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
      std::vector<double>& point = points_[subject];
      std::copy(theta.begin(), theta.end(), point.begin());
      std::copy(
        sigma.begin(), sigma.end(), point.begin() + n_theta_ + n_eta_);
      std::copy(
        omega.begin(), omega.end(),
        point.begin() + n_theta_ + n_eta_ + n_sigma_);
      ObjectiveTape& tape = *tapes_[subject];
      const std::vector<double> value = tape.fun.Forward(0, point, messages);
      require_unchanged_path(tape.fun, "persistent SAEM fixed-ETA objective");
      if (value.size() != 1U || !std::isfinite(value[0])) {
        throw std::domain_error(
          "A persistent SAEM fixed-ETA objective was non-finite.");
      }
      values[static_cast<R_xlen_t>(subject)] = value[0];
      const std::vector<double> derivative = tape.fun.Reverse(1, weight);
      require_unchanged_path(tape.fun, "persistent SAEM fixed-ETA gradient");
      if (derivative.size() != domain_) {
        throw std::runtime_error(
          "A persistent SAEM tape returned an invalid gradient length.");
      }
      for (std::size_t column = 0; column < domain_; ++column) {
        gradients(static_cast<int>(subject), static_cast<int>(column)) =
          derivative[column];
      }
      if ((subject + 1U) % 256U == 0U) Rcpp::checkUserInterrupt();
    }
    ++evaluations_;
    return Rcpp::List::create(
      Rcpp::Named("value") = values,
      Rcpp::Named("gradient") = gradients,
      Rcpp::Named("evaluations") = evaluations_);
  }

  Rcpp::List evaluate_aggregate(
      const Rcpp::NumericVector& theta, const Rcpp::NumericVector& sigma,
      const Rcpp::NumericVector& omega) {
    validate_parameters(theta, sigma, omega);
    double value_total = 0.0;
    Rcpp::NumericVector gradient_total(static_cast<R_xlen_t>(domain_));
    const std::vector<double> weight(1U, 1.0);
    std::ostringstream messages;
    for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
      std::vector<double>& point = points_[subject];
      std::copy(theta.begin(), theta.end(), point.begin());
      std::copy(
        sigma.begin(), sigma.end(), point.begin() + n_theta_ + n_eta_);
      std::copy(
        omega.begin(), omega.end(),
        point.begin() + n_theta_ + n_eta_ + n_sigma_);
      ObjectiveTape& tape = *tapes_[subject];
      const std::vector<double> value = tape.fun.Forward(0, point, messages);
      require_unchanged_path(
        tape.fun, "persistent aggregate SAEM fixed-ETA objective");
      if (value.size() != 1U || !std::isfinite(value[0])) {
        throw std::domain_error(
          "A persistent aggregate SAEM fixed-ETA objective was non-finite.");
      }
      value_total += value[0];
      const std::vector<double> derivative = tape.fun.Reverse(1, weight);
      require_unchanged_path(
        tape.fun, "persistent aggregate SAEM fixed-ETA gradient");
      if (derivative.size() != domain_) {
        throw std::runtime_error(
          "A persistent aggregate SAEM tape returned an invalid gradient length.");
      }
      for (std::size_t column = 0; column < domain_; ++column) {
        gradient_total[static_cast<R_xlen_t>(column)] += derivative[column];
      }
      if ((subject + 1U) % 256U == 0U) Rcpp::checkUserInterrupt();
    }
    ++evaluations_;
    if (!tapes_.empty()) {
      gradient_total.attr("names") = Rcpp::wrap(tapes_.front()->domain_names);
    }
    return Rcpp::List::create(
      Rcpp::Named("value") = value_total,
      Rcpp::Named("gradient") = gradient_total,
      Rcpp::Named("evaluations") = evaluations_);
  }

 private:
  int n_theta_ = 0;
  int n_eta_ = 0;
  int n_sigma_ = 0;
  int n_omega_ = 0;
  std::size_t domain_ = 0U;
  int evaluations_ = 0;
  std::vector<ObjectiveTape*> tapes_;
  std::vector<std::vector<double>> points_;

  void validate_parameters(
      const Rcpp::NumericVector& theta, const Rcpp::NumericVector& sigma,
      const Rcpp::NumericVector& omega) const {
    if (theta.size() != n_theta_ || sigma.size() != n_sigma_ ||
        omega.size() != n_omega_) {
      throw std::invalid_argument(
        "Persistent SAEM population parameter dimensions changed.");
    }
    for (double value : theta) if (!std::isfinite(value)) {
      throw std::invalid_argument("Persistent SAEM THETAs must be finite.");
    }
    for (double value : sigma) if (!std::isfinite(value)) {
      throw std::invalid_argument("Persistent SAEM SIGMAs must be finite.");
    }
    for (double value : omega) if (!std::isfinite(value)) {
      throw std::invalid_argument("Persistent SAEM OMEGAs must be finite.");
    }
  }
};

// Persistent complete-data expectation shared by ITS, IMP and SAEM.  ETA
// support points and their normalized weights may be replaced (ITS/IMP) or
// advanced through the exact Robbins--Monro recurrence (SAEM), while the
// subject tapes, dynamic inputs and point buffers remain resident.  All
// reductions deliberately retain subject-then-support order so the
// compatibility path does not acquire scheduler-dependent floating-point
// results.
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
