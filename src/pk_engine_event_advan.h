
template <class Base>
class CppADRecordingGuard {
 public:
  CppADRecordingGuard() = default;
  CppADRecordingGuard(const CppADRecordingGuard&) = delete;
  CppADRecordingGuard& operator=(const CppADRecordingGuard&) = delete;
  ~CppADRecordingGuard() {
    if (active_) CppAD::AD<Base>::abort_recording();
  }
  void release() noexcept { active_ = false; }

 private:
  bool active_ = true;
};

struct AffineMap {
  Matrix transition;
  Vector offset;
};

struct Topology {
  Matrix k;
  std::vector<std::string> state_names;
  int default_dose = 0;
  int default_observation = 0;
  std::vector<double> default_scales;
};

struct MatrixFlow {
  int from = 0;
  int to = -1;
  std::string type;
  std::string parameter;
  std::string volume_parameter;
};

struct MatrixGraph {
  std::vector<std::string> names;
  std::vector<std::string> scale_parameters;
  std::vector<MatrixFlow> flows;
  bool enabled = false;
};

struct ActiveInfusion {
  double end = 0.0;
  int compartment = 0;
  double rate = 0.0;
};

struct OdeControl {
  double rtol = 1e-8;
  double atol = 1e-10;
  int max_steps = 100000;
  double initial_step = 0.0;
};

using Parameters = std::unordered_map<std::string, double>;

inline bool finite_positive(double value) {
  return std::isfinite(value) && value > 0.0;
}

double get_parameter(const Parameters& parameters,
                     std::initializer_list<const char*> names,
                     double fallback = std::numeric_limits<double>::quiet_NaN()) {
  for (const char* name : names) {
    auto it = parameters.find(name);
    if (it != parameters.end() && std::isfinite(it->second)) return it->second;
  }
  return fallback;
}

double get_positive(const Parameters& parameters,
                    std::initializer_list<const char*> names,
                    double fallback = std::numeric_limits<double>::quiet_NaN()) {
  for (const char* name : names) {
    auto it = parameters.find(name);
    if (it != parameters.end() && finite_positive(it->second)) return it->second;
  }
  return fallback;
}

void require_positive(double value, const std::string& name, int advan) {
  if (!finite_positive(value)) {
    throw std::domain_error("ADVAN" + std::to_string(advan) +
                            " requires a positive " + name + ".");
  }
}

Matrix matrix_exp(const Matrix& matrix) {
  if (matrix.rows() != matrix.cols()) {
    throw std::invalid_argument("Matrix exponential requires a square matrix.");
  }
  if (!matrix.allFinite()) {
    throw std::domain_error("Matrix exponential input contains non-finite values.");
  }
  Matrix result = matrix.exp();
  if (!result.allFinite()) {
    throw std::domain_error("Matrix exponential produced non-finite values.");
  }
  return result;
}

AffineMap affine_map(const Matrix& k, const Vector& input, double dt) {
  const Eigen::Index n = k.rows();
  if (k.cols() != n || input.size() != n) {
    throw std::invalid_argument("Affine propagation dimensions are inconsistent.");
  }
  if (dt < -1e-12 || !std::isfinite(dt)) {
    throw std::domain_error("Propagation interval must be finite and non-negative.");
  }
  if (dt <= 0.0) return {Matrix::Identity(n, n), Vector::Zero(n)};
  Matrix augmented = Matrix::Zero(n + 1, n + 1);
  augmented.topLeftCorner(n, n) = k;
  augmented.topRightCorner(n, 1) = input;
  Matrix exponential = matrix_exp(augmented * dt);
  return {exponential.topLeftCorner(n, n), exponential.topRightCorner(n, 1)};
}

Vector propagate(const Matrix& k, const Vector& input, double dt, const Vector& state) {
  AffineMap map = affine_map(k, input, dt);
  return map.transition * state + map.offset;
}

Vector solve_periodic(const Matrix& transition, const Vector& offset,
                      const std::string& context) {
  Matrix system = Matrix::Identity(transition.rows(), transition.cols()) - transition;
  Eigen::FullPivLU<Matrix> lu(system);
  lu.setThreshold(1e-12);
  if (!lu.isInvertible()) {
    throw std::domain_error(context +
      " steady state does not exist or is numerically singular (I - Phi is not invertible).");
  }
  Vector solution = lu.solve(offset);
  const double scale = std::max(1.0, offset.norm());
  const double residual = (system * solution - offset).norm() / scale;
  if (!solution.allFinite() || !std::isfinite(residual) || residual > 1e-8) {
    throw std::domain_error(context + " steady-state solve is ill-conditioned.");
  }
  return solution;
}

Topology build_topology(int advan, const Parameters& p) {
  const bool oral = advan == 2 || advan == 4 || advan == 12;
  const int n = (advan == 1 ? 1 :
                 advan == 2 || advan == 3 ? 2 :
                 advan == 4 || advan == 11 ? 3 :
                 advan == 12 ? 4 : 0);
  if (n == 0) {
    throw std::invalid_argument("The current analytical engine supports ADVAN1-4/11/12.");
  }
  Topology topology;
  topology.k = Matrix::Zero(n, n);
  topology.default_scales.assign(static_cast<std::size_t>(n), 1.0);

  const double vc = get_positive(p, {"VC", "V1", "V"});
  const double vp1 = get_positive(p, {"VP", "VP1", "V2"});
  const double vp2 = get_positive(p, {"VP2", "V3"});
  double cl = get_positive(p, {"CL"});
  const double q1 = get_positive(p, {"Q2", "Q", "Q1"});
  const double q2 = get_positive(p, {"Q3", "Q4"});

  if (advan == 1) {
    double k10 = get_positive(p, {"K10", "K"});
    if (!finite_positive(k10) && finite_positive(cl) && finite_positive(vc)) k10 = cl / vc;
    require_positive(k10, "K10 or CL/V", advan);
    topology.k(0, 0) = -k10;
    topology.state_names = {"CENTRAL"};
    topology.default_scales[0] = finite_positive(vc) ? vc : 1.0;
  } else if (advan == 2) {
    const double ka = get_positive(p, {"KA"});
    double k20 = get_positive(p, {"K20", "K10", "K"});
    if (!finite_positive(k20) && finite_positive(cl) && finite_positive(vc)) k20 = cl / vc;
    require_positive(ka, "KA", advan);
    require_positive(k20, "K20 or CL/V", advan);
    topology.k(0, 0) = -ka;
    topology.k(1, 0) = ka;
    topology.k(1, 1) = -k20;
    topology.state_names = {"DEPOT", "CENTRAL"};
    topology.default_observation = 1;
    topology.default_scales[1] = finite_positive(vc) ? vc : 1.0;
  } else if (advan == 3) {
    double k10 = get_positive(p, {"K10"});
    double k12 = get_positive(p, {"K12"});
    double k21 = get_positive(p, {"K21"});
    if (!finite_positive(k10) && finite_positive(cl) && finite_positive(vc)) k10 = cl / vc;
    if (!finite_positive(k12) && finite_positive(q1) && finite_positive(vc)) k12 = q1 / vc;
    if (!finite_positive(k21) && finite_positive(q1) && finite_positive(vp1)) k21 = q1 / vp1;
    require_positive(k10, "K10 or CL/V1", advan);
    require_positive(k12, "K12 or Q/V1", advan);
    require_positive(k21, "K21 or Q/V2", advan);
    topology.k << -(k10 + k12), k21,
                   k12, -k21;
    topology.state_names = {"CENTRAL", "PERIPHERAL1"};
    topology.default_scales[0] = finite_positive(vc) ? vc : 1.0;
    topology.default_scales[1] = finite_positive(vp1) ? vp1 : 1.0;
  } else if (advan == 4) {
    const double ka = get_positive(p, {"KA"});
    double k20 = get_positive(p, {"K20", "K10"});
    double k23 = get_positive(p, {"K23", "K12"});
    double k32 = get_positive(p, {"K32", "K21"});
    if (!finite_positive(k20) && finite_positive(cl) && finite_positive(vc)) k20 = cl / vc;
    if (!finite_positive(k23) && finite_positive(q1) && finite_positive(vc)) k23 = q1 / vc;
    if (!finite_positive(k32) && finite_positive(q1) && finite_positive(vp1)) k32 = q1 / vp1;
    require_positive(ka, "KA", advan);
    require_positive(k20, "K20 or CL/VC", advan);
    require_positive(k23, "K23 or Q/VC", advan);
    require_positive(k32, "K32 or Q/VP", advan);
    topology.k(0, 0) = -ka;
    topology.k(1, 0) = ka;
    topology.k(1, 1) = -(k20 + k23);
    topology.k(1, 2) = k32;
    topology.k(2, 1) = k23;
    topology.k(2, 2) = -k32;
    topology.state_names = {"DEPOT", "CENTRAL", "PERIPHERAL1"};
    topology.default_observation = 1;
    topology.default_scales[1] = finite_positive(vc) ? vc : 1.0;
    topology.default_scales[2] = finite_positive(vp1) ? vp1 : 1.0;
  } else if (advan == 11) {
    double k10 = get_positive(p, {"K10"});
    double k12 = get_positive(p, {"K12"});
    double k21 = get_positive(p, {"K21"});
    double k13 = get_positive(p, {"K13"});
    double k31 = get_positive(p, {"K31"});
    if (!finite_positive(k10) && finite_positive(cl) && finite_positive(vc)) k10 = cl / vc;
    if (!finite_positive(k12) && finite_positive(q1) && finite_positive(vc)) k12 = q1 / vc;
    if (!finite_positive(k21) && finite_positive(q1) && finite_positive(vp1)) k21 = q1 / vp1;
    if (!finite_positive(k13) && finite_positive(q2) && finite_positive(vc)) k13 = q2 / vc;
    if (!finite_positive(k31) && finite_positive(q2) && finite_positive(vp2)) k31 = q2 / vp2;
    require_positive(k10, "K10 or CL/V1", advan);
    require_positive(k12, "K12 or Q2/V1", advan);
    require_positive(k21, "K21 or Q2/V2", advan);
    require_positive(k13, "K13 or Q3/V1", advan);
    require_positive(k31, "K31 or Q3/V3", advan);
    topology.k(0, 0) = -(k10 + k12 + k13);
    topology.k(0, 1) = k21;
    topology.k(0, 2) = k31;
    topology.k(1, 0) = k12;
    topology.k(1, 1) = -k21;
    topology.k(2, 0) = k13;
    topology.k(2, 2) = -k31;
    topology.state_names = {"CENTRAL", "PERIPHERAL1", "PERIPHERAL2"};
    topology.default_scales[0] = finite_positive(vc) ? vc : 1.0;
    topology.default_scales[1] = finite_positive(vp1) ? vp1 : 1.0;
    topology.default_scales[2] = finite_positive(vp2) ? vp2 : 1.0;
  } else if (advan == 12) {
    const double ka = get_positive(p, {"KA"});
    double k20 = get_positive(p, {"K20", "K10"});
    double k23 = get_positive(p, {"K23", "K12"});
    double k32 = get_positive(p, {"K32", "K21"});
    double k24 = get_positive(p, {"K24", "K13"});
    double k42 = get_positive(p, {"K42", "K31"});
    if (!finite_positive(k20) && finite_positive(cl) && finite_positive(vc)) k20 = cl / vc;
    if (!finite_positive(k23) && finite_positive(q1) && finite_positive(vc)) k23 = q1 / vc;
    if (!finite_positive(k32) && finite_positive(q1) && finite_positive(vp1)) k32 = q1 / vp1;
    if (!finite_positive(k24) && finite_positive(q2) && finite_positive(vc)) k24 = q2 / vc;
    if (!finite_positive(k42) && finite_positive(q2) && finite_positive(vp2)) k42 = q2 / vp2;
    require_positive(ka, "KA", advan);
    require_positive(k20, "K20 or CL/VC", advan);
    require_positive(k23, "K23 or Q2/VC", advan);
    require_positive(k32, "K32 or Q2/VP1", advan);
    require_positive(k24, "K24 or Q3/VC", advan);
    require_positive(k42, "K42 or Q3/VP2", advan);
    topology.k(0, 0) = -ka;
    topology.k(1, 0) = ka;
    topology.k(1, 1) = -(k20 + k23 + k24);
    topology.k(1, 2) = k32;
    topology.k(1, 3) = k42;
    topology.k(2, 1) = k23;
    topology.k(2, 2) = -k32;
    topology.k(3, 1) = k24;
    topology.k(3, 3) = -k42;
    topology.state_names = {"DEPOT", "CENTRAL", "PERIPHERAL1", "PERIPHERAL2"};
    topology.default_observation = 1;
    topology.default_scales[1] = finite_positive(vc) ? vc : 1.0;
    topology.default_scales[2] = finite_positive(vp1) ? vp1 : 1.0;
    topology.default_scales[3] = finite_positive(vp2) ? vp2 : 1.0;
  }
  (void)oral;
  return topology;
}

Topology build_graph_topology(const MatrixGraph& graph, const Parameters& p) {
  if (!graph.enabled || graph.names.empty()) {
    throw std::invalid_argument("Matrix graph is empty.");
  }
  const int n = static_cast<int>(graph.names.size());
  Topology topology;
  topology.k = Matrix::Zero(n, n);
  topology.state_names = graph.names;
  topology.default_scales.assign(static_cast<std::size_t>(n), 1.0);
  for (int i = 0; i < n; ++i) {
    if (i < static_cast<int>(graph.scale_parameters.size()) &&
        !graph.scale_parameters[static_cast<std::size_t>(i)].empty()) {
      const std::string& name = graph.scale_parameters[static_cast<std::size_t>(i)];
      auto it = p.find(name);
      if (it == p.end() || !finite_positive(it->second)) {
        throw std::domain_error("Matrix graph scale parameter '" + name + "' must be positive.");
      }
      topology.default_scales[static_cast<std::size_t>(i)] = it->second;
    }
  }
  for (const MatrixFlow& flow : graph.flows) {
    auto parameter = p.find(flow.parameter);
    if (parameter == p.end() || !finite_positive(parameter->second)) {
      throw std::domain_error("Matrix graph flow parameter '" + flow.parameter + "' must be positive.");
    }
    double rate = parameter->second;
    if (flow.type == "clearance") {
      auto volume = p.find(flow.volume_parameter);
      if (volume == p.end() || !finite_positive(volume->second)) {
        throw std::domain_error("Matrix graph volume parameter '" +
                                flow.volume_parameter + "' must be positive.");
      }
      rate /= volume->second;
    }
    if (flow.from < 0 || flow.from >= n || flow.to >= n) {
      throw std::logic_error("Matrix graph contains an invalid compiled compartment index.");
    }
    topology.k(flow.from, flow.from) -= rate;
    if (flow.to >= 0) topology.k(flow.to, flow.from) += rate;
  }
  return topology;
}

struct ResidualGroupSpec {
  std::string label;
  std::vector<double> dvid;
  std::vector<std::string> source;
  std::vector<int> index;
  std::vector<double> value;
  std::string transform = "tanh";
};

class ModelEngine {
 public:
  int advan;
  int trans;
  int model_ss;
  int dose_cmp;
  int obs_cmp;
  int n_theta;
  int n_eta;
  int n_state;
  bool direct_prediction = false;
  std::string solver;
  std::string error_type;
  std::string omega_type;
  std::string sigma_correlation;
  std::string sigma_parameterization = "sd";
  std::string blq_method;
  double ar1_rho = 0.0;
  std::string ar1_parameter_source = "fixed";
  int ar1_parameter_index = -1;
  std::string ar1_transform = "tanh";
  double lloq = std::numeric_limits<double>::quiet_NaN();
  int iov = 0;
  bool specialized_advan = true;
  std::vector<int> omega_rows;
  std::vector<int> omega_cols;
  bool re_enabled = false;
  std::vector<std::vector<int>> re_blocks;
  std::vector<double> mixture_probabilities;
  std::vector<ResidualGroupSpec> residual_groups;
  std::shared_ptr<const libertad::Program> pred;
  std::shared_ptr<const libertad::Program> post_pred;
  std::shared_ptr<const libertad::Program> des;
  std::shared_ptr<const libertad::Program> alg;
  std::shared_ptr<const libertad::Program> error;
  std::vector<std::size_t> all_outputs;
  std::vector<std::size_t> post_all_outputs;
  std::vector<std::size_t> likelihood_output;
  std::string likelihood_scale;
  bool hmm_enabled = false;
  bool hmm_by_dvid = true;
  bool hmm_continuous = false;
  int hmm_states = 0;
  std::vector<std::string> hmm_state_names;
  std::vector<std::size_t> hmm_outputs;
  std::string hmm_initial_scale;
  std::string hmm_transition_scale;
  std::string hmm_rate_scale;
  std::string hmm_emission_scale;
  bool kalman_enabled = false;
  bool kalman_by_dvid = true;
  bool kalman_prediction_baseline = true;
  std::string kalman_filter_type = "linear";
  int kalman_states = 0;
  std::vector<std::string> kalman_state_names;
  std::vector<std::string> kalman_state_inputs;
  double kalman_jacobian_step = 1e-5;
  double kalman_ukf_alpha = 0.5;
  double kalman_ukf_beta = 2.0;
  double kalman_ukf_kappa = 0.0;
  int kalman_particles = 256;
  double kalman_ess_threshold = 0.5;
  int kalman_seed = 20260721;
  std::string kalman_dynamics = "discrete";
  std::string kalman_sde_method = "euler";
  int kalman_sde_substeps = 8;
  std::vector<std::size_t> kalman_outputs;
  bool switching_enabled = false;
  int switching_regimes = 0;
  std::vector<std::string> switching_regime_names;
  std::string switching_initial_scale = "probability";
  std::string switching_transition_scale = "probability";
  std::vector<std::size_t> switching_outputs;
  std::vector<std::string> selected_output_names;
  std::vector<std::size_t> derivative_outputs;
  std::vector<std::size_t> algebraic_outputs;
  bool dde_enabled = false;
  double dde_step = 0.05;
  int dde_max_steps = 100000;
  double dde_minimum_delay = 0.0;
  std::vector<double> dde_history;
  std::vector<std::string> dde_lag_inputs;
  std::vector<int> dde_lag_states;
  std::vector<std::string> dde_lag_delays;
  bool dae_enabled = false;
  std::vector<std::string> dae_variables;
  std::vector<double> dae_initial;
  double dae_tolerance = 1e-9;
  int dae_maxit = 12;
  double dae_jacobian_step = 1e-6;
  std::vector<bool> dae_sparsity;
  std::vector<std::vector<int>> dae_block_rows;
  std::vector<std::vector<int>> dae_block_columns;
  std::vector<std::string> state_names;
  OdeControl ode_control;
  MatrixGraph matrix_graph;

  explicit ModelEngine(const Rcpp::List& spec)
      : advan(Rcpp::as<int>(spec["advan"])),
        trans(Rcpp::as<int>(spec["trans"])),
        model_ss(Rcpp::as<int>(spec["model_ss"])),
        dose_cmp(Rcpp::as<int>(spec["dose_cmp"])),
        obs_cmp(Rcpp::as<int>(spec["obs_cmp"])),
        n_theta(Rcpp::as<int>(spec["n_theta"])),
        n_eta(Rcpp::as<int>(spec["n_eta"])),
        n_state(Rcpp::as<int>(spec["n_state"])),
        direct_prediction(
          spec.containsElementNamed("pred_mode") &&
          Rcpp::as<std::string>(spec["pred_mode"]) == "pred"),
        solver(Rcpp::as<std::string>(spec["solver"])),
        error_type(Rcpp::as<std::string>(spec["error_type"])),
        pred(std::make_shared<const libertad::Program>(Rcpp::as<Rcpp::List>(spec["pred_ir"]))) {
    all_outputs.resize(pred->output_names.size());
    std::iota(all_outputs.begin(), all_outputs.end(), 0U);
    Rcpp::RObject post_pred_ir = spec.containsElementNamed("post_pred_ir") ?
      Rcpp::RObject(spec["post_pred_ir"]) : Rcpp::RObject(R_NilValue);
    if (!Rf_isNull(post_pred_ir)) {
      post_pred = std::make_shared<const libertad::Program>(
        Rcpp::as<Rcpp::List>(post_pred_ir));
      post_all_outputs.resize(post_pred->output_names.size());
      std::iota(post_all_outputs.begin(), post_all_outputs.end(), 0U);
    }
    if (spec.containsElementNamed("output_names")) {
      selected_output_names = Rcpp::as<std::vector<std::string>>(spec["output_names"]);
      // Validate the serialized selection at engine construction rather than
      // failing part-way through a long estimation or simulation.
      for (const std::string& name : selected_output_names) {
        const bool in_pk = std::find(
          pred->output_names.begin(), pred->output_names.end(), name
        ) != pred->output_names.end();
        const bool in_post = post_pred && std::find(
          post_pred->output_names.begin(), post_pred->output_names.end(), name
        ) != post_pred->output_names.end();
        if (!in_pk && !in_post) {
          throw std::invalid_argument(
            "Selected model output '" + name + "' is not compiled.");
        }
      }
    }
    Rcpp::RObject des_ir = spec["des_ir"];
    if (!Rf_isNull(des_ir)) {
      des = std::make_shared<const libertad::Program>(Rcpp::as<Rcpp::List>(des_ir));
      derivative_outputs.reserve(static_cast<std::size_t>(n_state));
      for (int i = 1; i <= n_state; ++i) {
        derivative_outputs.push_back(des->select_outputs({"DADT_" + std::to_string(i)}).front());
      }
    }
    Rcpp::RObject dde_object = spec.containsElementNamed("dde_config") ?
      Rcpp::RObject(spec["dde_config"]) : Rcpp::RObject(R_NilValue);
    if (!Rf_isNull(dde_object)) {
      Rcpp::List dde(dde_object);
      dde_enabled = true;
      dde_step = Rcpp::as<double>(dde["step"]);
      dde_max_steps = Rcpp::as<int>(dde["max_steps"]);
      dde_history = Rcpp::as<std::vector<double>>(dde["history"]);
      if (dde.containsElementNamed("minimum_delay") &&
          !Rf_isNull(dde["minimum_delay"])) {
        dde_minimum_delay = Rcpp::as<double>(dde["minimum_delay"]);
      } else {
        dde_minimum_delay = dde_step;
      }
      Rcpp::List lags(dde["lags"]);
      for (R_xlen_t lag = 0; lag < lags.size(); ++lag) {
        Rcpp::List value(lags[lag]);
        dde_lag_inputs.push_back(Rcpp::as<std::string>(value["input"]));
        dde_lag_states.push_back(Rcpp::as<int>(value["state"]) - 1);
        dde_lag_delays.push_back(Rcpp::as<std::string>(value["delay"]));
      }
      if (!des || dde_step <= 0.0 || dde_max_steps < 1 ||
          dde_history.size() != static_cast<std::size_t>(n_state) ||
          dde_lag_inputs.empty() || dde_lag_inputs.size() != dde_lag_states.size() ||
          dde_lag_inputs.size() != dde_lag_delays.size()) {
        throw std::invalid_argument("DDE configuration is inconsistent with the compiled model.");
      }
    }
    Rcpp::RObject alg_ir = spec.containsElementNamed("alg_ir") ?
      Rcpp::RObject(spec["alg_ir"]) : Rcpp::RObject(R_NilValue);
    Rcpp::RObject dae_object = spec.containsElementNamed("dae_config") ?
      Rcpp::RObject(spec["dae_config"]) : Rcpp::RObject(R_NilValue);
    if (!Rf_isNull(dae_object)) {
      if (Rf_isNull(alg_ir)) {
        throw std::invalid_argument("DAE configuration requires a compiled ALG residual program.");
      }
      alg = std::make_shared<const libertad::Program>(Rcpp::as<Rcpp::List>(alg_ir));
      Rcpp::List dae(dae_object);
      dae_enabled = true;
      dae_variables = Rcpp::as<std::vector<std::string>>(dae["variables"]);
      dae_initial = Rcpp::as<std::vector<double>>(dae["initial"]);
      dae_tolerance = Rcpp::as<double>(dae["tolerance"]);
      dae_maxit = Rcpp::as<int>(dae["maxit"]);
      dae_jacobian_step = Rcpp::as<double>(dae["jacobian_step"]);
      std::vector<std::string> names;
      names.reserve(dae_variables.size());
      for (std::size_t index = 0; index < dae_variables.size(); ++index) {
        names.push_back("DAE_RES_" + std::to_string(index + 1U));
      }
      algebraic_outputs = alg->select_outputs(names);
      if (dae.containsElementNamed("sparsity") && !Rf_isNull(dae["sparsity"])) {
        Rcpp::LogicalMatrix sparsity(dae["sparsity"]);
        if (sparsity.nrow() != static_cast<int>(dae_variables.size()) ||
            sparsity.ncol() != static_cast<int>(dae_variables.size())) {
          throw std::invalid_argument("DAE sparsity dimensions are inconsistent.");
        }
        for (int row = 0; row < sparsity.nrow(); ++row) {
          for (int column = 0; column < sparsity.ncol(); ++column) {
            dae_sparsity.push_back(sparsity(row, column) == TRUE);
          }
        }
      }
      if (!des || dae_variables.empty() || dae_initial.size() != dae_variables.size() ||
          dae_tolerance <= 0.0 || dae_maxit < 1 || dae_jacobian_step <= 0.0) {
        throw std::invalid_argument("DAE configuration is inconsistent with the compiled model.");
      }
      const int dimension = static_cast<int>(dae_variables.size());
      if (dae_sparsity.empty()) {
        std::vector<int> all(static_cast<std::size_t>(dimension));
        std::iota(all.begin(), all.end(), 0);
        dae_block_rows.push_back(all);
        dae_block_columns.push_back(all);
      } else {
        std::vector<bool> visited(static_cast<std::size_t>(2 * dimension), false);
        for (int start = 0; start < 2 * dimension; ++start) {
          if (visited[static_cast<std::size_t>(start)]) continue;
          std::queue<int> pending; pending.push(start);
          std::vector<int> rows, columns;
          visited[static_cast<std::size_t>(start)] = true;
          while (!pending.empty()) {
            const int node = pending.front(); pending.pop();
            if (node < dimension) {
              rows.push_back(node);
              for (int column = 0; column < dimension; ++column) {
                if (dae_sparsity[static_cast<std::size_t>(node * dimension + column)] &&
                    !visited[static_cast<std::size_t>(dimension + column)]) {
                  visited[static_cast<std::size_t>(dimension + column)] = true;
                  pending.push(dimension + column);
                }
              }
            } else {
              const int column = node - dimension;
              columns.push_back(column);
              for (int row = 0; row < dimension; ++row) {
                if (dae_sparsity[static_cast<std::size_t>(row * dimension + column)] &&
                    !visited[static_cast<std::size_t>(row)]) {
                  visited[static_cast<std::size_t>(row)] = true;
                  pending.push(row);
                }
              }
            }
          }
          if (rows.size() != columns.size() || rows.empty()) {
            throw std::invalid_argument(
              "Every DAE sparsity block must contain the same non-zero number of residuals and variables.");
          }
          dae_block_rows.push_back(std::move(rows));
          dae_block_columns.push_back(std::move(columns));
        }
      }
    }
    Rcpp::RObject error_ir = spec.containsElementNamed("error_ir") ?
      Rcpp::RObject(spec["error_ir"]) : Rcpp::RObject(R_NilValue);
    if (!Rf_isNull(error_ir)) {
      error = std::make_shared<const libertad::Program>(Rcpp::as<Rcpp::List>(error_ir));
      Rcpp::RObject hmm_object = spec.containsElementNamed("hmm_config") ?
        Rcpp::RObject(spec["hmm_config"]) : Rcpp::RObject(R_NilValue);
      if (!Rf_isNull(hmm_object)) {
        Rcpp::List hmm(hmm_object);
        hmm_enabled = true;
        hmm_state_names = Rcpp::as<std::vector<std::string>>(hmm["states"]);
        hmm_states = static_cast<int>(hmm_state_names.size());
        hmm_by_dvid = Rcpp::as<bool>(hmm["by_dvid"]);
        const std::string transition_type = hmm.containsElementNamed("transition_type") ?
          Rcpp::as<std::string>(hmm["transition_type"]) : "discrete";
        hmm_continuous = transition_type == "continuous";
        hmm_initial_scale = Rcpp::as<std::string>(hmm["initial_scale"]);
        if (hmm_continuous) {
          hmm_rate_scale = Rcpp::as<std::string>(hmm["rate_scale"]);
        } else {
          hmm_transition_scale = Rcpp::as<std::string>(hmm["transition_scale"]);
        }
        hmm_emission_scale = Rcpp::as<std::string>(hmm["emission_scale"]);
        if (hmm_states < 2 ||
            (hmm_initial_scale != "probability" && hmm_initial_scale != "log") ||
            (!hmm_continuous && hmm_transition_scale != "probability" &&
             hmm_transition_scale != "log") ||
            (hmm_continuous && hmm_rate_scale != "rate" && hmm_rate_scale != "log") ||
            (hmm_emission_scale != "likelihood" && hmm_emission_scale != "log")) {
          throw std::invalid_argument("Hidden Markov likelihood scales or state count are inconsistent.");
        }
        std::vector<std::string> names =
          Rcpp::as<std::vector<std::string>>(hmm["initial"]);
        const std::vector<std::string> emission =
          Rcpp::as<std::vector<std::string>>(hmm["emission"]);
        if (names.size() != static_cast<std::size_t>(hmm_states) ||
            emission.size() != static_cast<std::size_t>(hmm_states)) {
          throw std::invalid_argument("Hidden Markov output dimensions are inconsistent.");
        }
        if (hmm_continuous) {
          Rcpp::CharacterMatrix generator(hmm["generator"]);
          if (generator.nrow() != hmm_states || generator.ncol() != hmm_states) {
            throw std::invalid_argument("Continuous-time HMM generator dimensions are inconsistent.");
          }
          names.reserve(static_cast<std::size_t>(hmm_states * (hmm_states + 1)));
          for (int from = 0; from < hmm_states; ++from) {
            for (int to = 0; to < hmm_states; ++to) {
              if (from != to) {
                names.push_back(Rcpp::as<std::string>(generator(from, to)));
              }
            }
          }
        } else {
          Rcpp::CharacterMatrix transition(hmm["transition"]);
          if (transition.nrow() != hmm_states || transition.ncol() != hmm_states) {
            throw std::invalid_argument("Discrete-time HMM transition dimensions are inconsistent.");
          }
          names.reserve(static_cast<std::size_t>(hmm_states * (hmm_states + 2)));
          for (int from = 0; from < hmm_states; ++from) {
            for (int to = 0; to < hmm_states; ++to) {
              names.push_back(Rcpp::as<std::string>(transition(from, to)));
            }
          }
        }
        names.insert(names.end(), emission.begin(), emission.end());
        hmm_outputs = error->select_outputs(names);
      } else {
        Rcpp::RObject kalman_object = spec.containsElementNamed("kalman_config") ?
          Rcpp::RObject(spec["kalman_config"]) : Rcpp::RObject(R_NilValue);
        if (!Rf_isNull(kalman_object)) {
          Rcpp::List kalman(kalman_object);
          kalman_enabled = true;
          kalman_state_names = Rcpp::as<std::vector<std::string>>(kalman["states"]);
          kalman_states = static_cast<int>(kalman_state_names.size());
          kalman_by_dvid = Rcpp::as<bool>(kalman["by_dvid"]);
          kalman_prediction_baseline =
            Rcpp::as<std::string>(kalman["baseline"]) == "prediction";
          kalman_filter_type = kalman.containsElementNamed("filter") ?
            Rcpp::as<std::string>(kalman["filter"]) : "linear";
          kalman_state_inputs = kalman.containsElementNamed("state_inputs") ?
            Rcpp::as<std::vector<std::string>>(kalman["state_inputs"]) :
            std::vector<std::string>();
          kalman_jacobian_step = kalman.containsElementNamed("jacobian_step") ?
            Rcpp::as<double>(kalman["jacobian_step"]) : 1e-5;
          kalman_ukf_alpha = kalman.containsElementNamed("ukf_alpha") ?
            Rcpp::as<double>(kalman["ukf_alpha"]) : 0.5;
          kalman_ukf_beta = kalman.containsElementNamed("ukf_beta") ?
            Rcpp::as<double>(kalman["ukf_beta"]) : 2.0;
          kalman_ukf_kappa = kalman.containsElementNamed("ukf_kappa") ?
            Rcpp::as<double>(kalman["ukf_kappa"]) : 0.0;
          kalman_particles = kalman.containsElementNamed("particles") ?
            Rcpp::as<int>(kalman["particles"]) : 256;
          kalman_ess_threshold = kalman.containsElementNamed("ess_threshold") ?
            Rcpp::as<double>(kalman["ess_threshold"]) : 0.5;
          kalman_seed = kalman.containsElementNamed("seed") ?
            Rcpp::as<int>(kalman["seed"]) : 20260721;
          kalman_dynamics = kalman.containsElementNamed("dynamics") ?
            Rcpp::as<std::string>(kalman["dynamics"]) : "discrete";
          kalman_sde_method = kalman.containsElementNamed("sde_method") ?
            Rcpp::as<std::string>(kalman["sde_method"]) : "euler";
          kalman_sde_substeps = kalman.containsElementNamed("sde_substeps") ?
            Rcpp::as<int>(kalman["sde_substeps"]) : 8;
          if (kalman_filter_type != "linear" && kalman_filter_type != "ekf" &&
              kalman_filter_type != "ukf" && kalman_filter_type != "particle") {
            throw std::invalid_argument("Unknown state-space filter type.");
          }
          if ((kalman_dynamics != "discrete" && kalman_dynamics != "sde") ||
              (kalman_sde_method != "euler" && kalman_sde_method != "milstein") ||
              kalman_sde_substeps < 1) {
            throw std::invalid_argument("Invalid state-space dynamics controls.");
          }
          std::vector<std::string> names =
            Rcpp::as<std::vector<std::string>>(kalman["initial_mean"]);
          const std::vector<std::string> observation =
            Rcpp::as<std::vector<std::string>>(kalman["observation"]);
          const std::vector<std::string> observation_variance =
            Rcpp::as<std::vector<std::string>>(kalman["observation_variance"]);
          Rcpp::CharacterMatrix initial_covariance(kalman["initial_covariance"]);
          Rcpp::CharacterMatrix process_covariance(kalman["process_covariance"]);
          if (kalman_states < 1 || names.size() != static_cast<std::size_t>(kalman_states) ||
              observation_variance.size() != 1U ||
              initial_covariance.nrow() != kalman_states ||
              initial_covariance.ncol() != kalman_states ||
              process_covariance.nrow() != kalman_states ||
              process_covariance.ncol() != kalman_states) {
            throw std::invalid_argument("Kalman state-space output dimensions are inconsistent.");
          }
          names.reserve(static_cast<std::size_t>(
            3 * kalman_states * kalman_states + 2 * kalman_states + 1));
          auto append_matrix = [&](const Rcpp::CharacterMatrix& matrix) {
            for (int row = 0; row < kalman_states; ++row) {
              for (int column = 0; column < kalman_states; ++column) {
                names.push_back(Rcpp::as<std::string>(matrix(row, column)));
              }
            }
          };
          append_matrix(initial_covariance);
          if (kalman_filter_type == "linear") {
            Rcpp::CharacterMatrix transition(kalman["transition"]);
            if (transition.nrow() != kalman_states || transition.ncol() != kalman_states ||
                observation.size() != static_cast<std::size_t>(kalman_states)) {
              throw std::invalid_argument("Linear Kalman transition/observation dimensions are inconsistent.");
            }
            append_matrix(transition);
          } else {
            const std::vector<std::string> transition =
              Rcpp::as<std::vector<std::string>>(kalman["transition"]);
            if (transition.size() != static_cast<std::size_t>(kalman_states) ||
                observation.size() != 1U ||
                kalman_state_inputs.size() != static_cast<std::size_t>(kalman_states)) {
              throw std::invalid_argument("Nonlinear state-space transition/observation dimensions are inconsistent.");
            }
            names.insert(names.end(), transition.begin(), transition.end());
          }
          append_matrix(process_covariance);
          names.insert(names.end(), observation.begin(), observation.end());
          names.push_back(observation_variance.front());
          kalman_outputs = error->select_outputs(names);
          if (kalman.containsElementNamed("switching") &&
              !Rf_isNull(kalman["switching"])) {
            if (kalman_filter_type != "particle") {
              throw std::invalid_argument("Switching state-space inference requires the particle filter.");
            }
            Rcpp::List switching(kalman["switching"]);
            switching_regime_names =
              Rcpp::as<std::vector<std::string>>(switching["regimes"]);
            switching_regimes = static_cast<int>(switching_regime_names.size());
            switching_initial_scale =
              Rcpp::as<std::string>(switching["initial_scale"]);
            switching_transition_scale =
              Rcpp::as<std::string>(switching["transition_scale"]);
            std::vector<std::string> switching_names =
              Rcpp::as<std::vector<std::string>>(switching["initial"]);
            Rcpp::CharacterMatrix regime_transition(switching["transition"]);
            Rcpp::List state_transition(switching["state_transition"]);
            Rcpp::List process(switching["process"]);
            const std::vector<std::string> observation =
              Rcpp::as<std::vector<std::string>>(switching["observation"]);
            const std::vector<std::string> observation_variance =
              Rcpp::as<std::vector<std::string>>(switching["observation_variance"]);
            if (switching_regimes < 2 ||
                switching_names.size() != static_cast<std::size_t>(switching_regimes) ||
                regime_transition.nrow() != switching_regimes ||
                regime_transition.ncol() != switching_regimes ||
                state_transition.size() != switching_regimes ||
                process.size() != switching_regimes ||
                observation.size() != static_cast<std::size_t>(switching_regimes) ||
                observation_variance.size() != static_cast<std::size_t>(switching_regimes)) {
              throw std::invalid_argument("Switching state-space dimensions are inconsistent.");
            }
            for (int from = 0; from < switching_regimes; ++from) {
              for (int to = 0; to < switching_regimes; ++to) {
                switching_names.push_back(
                  Rcpp::as<std::string>(regime_transition(from, to)));
              }
            }
            for (int regime = 0; regime < switching_regimes; ++regime) {
              const std::vector<std::string> local_transition =
                Rcpp::as<std::vector<std::string>>(state_transition[regime]);
              Rcpp::CharacterMatrix local_process(process[regime]);
              if (local_transition.size() != static_cast<std::size_t>(kalman_states) ||
                  local_process.nrow() != kalman_states ||
                  local_process.ncol() != kalman_states) {
                throw std::invalid_argument("Switching regime state dimensions are inconsistent.");
              }
              switching_names.insert(
                switching_names.end(), local_transition.begin(), local_transition.end());
              for (int matrix_row = 0; matrix_row < kalman_states; ++matrix_row) {
                for (int matrix_column = 0; matrix_column < kalman_states; ++matrix_column) {
                  switching_names.push_back(
                    Rcpp::as<std::string>(local_process(matrix_row, matrix_column)));
                }
              }
              switching_names.push_back(observation[static_cast<std::size_t>(regime)]);
              switching_names.push_back(
                observation_variance[static_cast<std::size_t>(regime)]);
            }
            switching_outputs = error->select_outputs(switching_names);
            switching_enabled = true;
          }
        } else {
        const std::string output = Rcpp::as<std::string>(spec["likelihood_output"]);
        likelihood_scale = Rcpp::as<std::string>(spec["likelihood_scale"]);
        likelihood_output = error->select_outputs({output});
        if (likelihood_output.size() != 1U ||
            (likelihood_scale != "log" && likelihood_scale != "likelihood")) {
          throw std::invalid_argument("User likelihood specification is inconsistent.");
        }
        }
      }
    }
    state_names = Rcpp::as<std::vector<std::string>>(spec["state_names"]);
    if (state_names.size() != static_cast<std::size_t>(n_state)) {
      state_names.clear();
      for (int i = 1; i <= n_state; ++i) state_names.push_back("COMPARTMENT" + std::to_string(i));
    }
    Rcpp::List control = spec["ode_control"];
    ode_control.rtol = Rcpp::as<double>(control["rtol"]);
    ode_control.atol = Rcpp::as<double>(control["atol"]);
    ode_control.max_steps = Rcpp::as<int>(control["max_steps"]);
    ode_control.initial_step = Rcpp::as<double>(control["initial_step"]);
    Rcpp::List likelihood = spec["lik_config"];
    error_type = Rcpp::as<std::string>(likelihood["error"]);
    omega_type = Rcpp::as<std::string>(likelihood["omega"]);
    sigma_correlation = Rcpp::as<std::string>(likelihood["sigma_corr"]);
    if (likelihood.containsElementNamed("sigma_parameterization")) {
      sigma_parameterization = Rcpp::as<std::string>(likelihood["sigma_parameterization"]);
    }
    blq_method = Rcpp::as<std::string>(likelihood["blq_method"]);
    ar1_rho = Rcpp::as<double>(likelihood["ar1_rho"]);
    if (likelihood.containsElementNamed("ar1_source")) {
      ar1_parameter_source = Rcpp::as<std::string>(likelihood["ar1_source"]);
    }
    if (likelihood.containsElementNamed("ar1_index")) {
      ar1_parameter_index = Rcpp::as<int>(likelihood["ar1_index"]) - 1;
    }
    if (likelihood.containsElementNamed("ar1_transform")) {
      ar1_transform = Rcpp::as<std::string>(likelihood["ar1_transform"]);
    }
    lloq = Rcpp::as<double>(likelihood["lloq"]);
    iov = Rcpp::as<int>(likelihood["iov"]);
    if (likelihood.containsElementNamed("residual_groups") &&
        !Rf_isNull(likelihood["residual_groups"])) {
      Rcpp::List groups(likelihood["residual_groups"]);
      residual_groups.reserve(groups.size());
      for (R_xlen_t group_index = 0; group_index < groups.size(); ++group_index) {
        Rcpp::List input(groups[group_index]);
        ResidualGroupSpec group;
        group.label = Rcpp::as<std::string>(input["label"]);
        group.dvid = Rcpp::as<std::vector<double>>(input["dvid"]);
        group.transform = Rcpp::as<std::string>(input["parameter_transform"]);
        Rcpp::CharacterMatrix source(input["source"]);
        Rcpp::IntegerMatrix index(input["index"]);
        Rcpp::NumericMatrix value(input["value"]);
        const int dimension = static_cast<int>(group.dvid.size());
        if (dimension < 2 || source.nrow() != dimension || source.ncol() != dimension ||
            index.nrow() != dimension || index.ncol() != dimension ||
            value.nrow() != dimension || value.ncol() != dimension) {
          throw std::invalid_argument("Correlated residual group dimensions are inconsistent.");
        }
        group.source.reserve(static_cast<std::size_t>(dimension * dimension));
        group.index.reserve(static_cast<std::size_t>(dimension * dimension));
        group.value.reserve(static_cast<std::size_t>(dimension * dimension));
        for (int row = 0; row < dimension; ++row) {
          for (int column = 0; column < dimension; ++column) {
            group.source.push_back(Rcpp::as<std::string>(source(row, column)));
            group.index.push_back(index(row, column) - 1);
            group.value.push_back(value(row, column));
          }
        }
        residual_groups.push_back(std::move(group));
      }
    }
    if (spec.containsElementNamed("specialized_advan")) {
      specialized_advan = Rcpp::as<bool>(spec["specialized_advan"]);
    }
    Rcpp::IntegerVector omega_row = spec["omega_row"];
    Rcpp::IntegerVector omega_col = spec["omega_col"];
    if (omega_row.size() != omega_col.size()) {
      throw std::invalid_argument("OMEGA row/column vectors have different lengths.");
    }
    omega_rows.reserve(static_cast<std::size_t>(omega_row.size()));
    omega_cols.reserve(static_cast<std::size_t>(omega_col.size()));
    for (R_xlen_t i = 0; i < omega_row.size(); ++i) {
      omega_rows.push_back(omega_row[i] - 1);
      omega_cols.push_back(omega_col[i] - 1);
    }
    Rcpp::RObject re_object = spec.containsElementNamed("re_config") ?
      Rcpp::RObject(spec["re_config"]) : Rcpp::RObject(R_NilValue);
    if (!Rf_isNull(re_object)) {
      Rcpp::List re_config(re_object);
      Rcpp::List blocks(re_config["blocks"]);
      re_blocks.reserve(blocks.size());
      std::vector<bool> assigned(static_cast<std::size_t>(n_eta), false);
      for (R_xlen_t block_index = 0; block_index < blocks.size(); ++block_index) {
        Rcpp::List block(blocks[block_index]);
        std::vector<int> indices = Rcpp::as<std::vector<int>>(block["etas"]);
        for (int& index : indices) {
          --index;
          if (index < 0 || index >= n_eta || assigned[static_cast<std::size_t>(index)]) {
            throw std::invalid_argument("Random-effect block ETA indices are inconsistent.");
          }
          assigned[static_cast<std::size_t>(index)] = true;
        }
        re_blocks.push_back(std::move(indices));
      }
      if (std::find(assigned.begin(), assigned.end(), false) != assigned.end()) {
        throw std::invalid_argument("Random-effect design does not assign every ETA.");
      }
      re_enabled = true;
    }
    Rcpp::RObject mixture_object = likelihood.containsElementNamed("mixtures") ?
      Rcpp::RObject(likelihood["mixtures"]) : Rcpp::RObject(R_NilValue);
    if (!Rf_isNull(mixture_object)) {
      Rcpp::List mixture(mixture_object);
      mixture_probabilities = Rcpp::as<std::vector<double>>(mixture["probability"]);
      if (mixture_probabilities.size() < 2U) {
        throw std::invalid_argument("A finite mixture requires at least two components.");
      }
    }

    Rcpp::RObject graph_object = spec["matrix_graph"];
    if (!Rf_isNull(graph_object)) {
      Rcpp::List graph(graph_object);
      matrix_graph.names = Rcpp::as<std::vector<std::string>>(graph["names"]);
      matrix_graph.scale_parameters =
        Rcpp::as<std::vector<std::string>>(graph["scale_parameter"]);
      Rcpp::IntegerVector from = graph["from"];
      Rcpp::IntegerVector to = graph["to"];
      std::vector<std::string> type = Rcpp::as<std::vector<std::string>>(graph["type"]);
      std::vector<std::string> parameter =
        Rcpp::as<std::vector<std::string>>(graph["parameter"]);
      std::vector<std::string> volume =
        Rcpp::as<std::vector<std::string>>(graph["volume_parameter"]);
      const std::size_t count = static_cast<std::size_t>(from.size());
      if (to.size() != from.size() || type.size() != count ||
          parameter.size() != count || volume.size() != count) {
        throw std::invalid_argument("Matrix graph flow vectors have different lengths.");
      }
      matrix_graph.flows.reserve(count);
      for (std::size_t i = 0; i < count; ++i) {
        matrix_graph.flows.push_back({
          from[static_cast<R_xlen_t>(i)] - 1,
          to[static_cast<R_xlen_t>(i)] - 1,
          type[i], parameter[i], volume[i]
        });
      }
      matrix_graph.enabled = true;
      n_state = static_cast<int>(matrix_graph.names.size());
      state_names = matrix_graph.names;
    }
  }

  bool is_ode() const { return static_cast<bool>(des); }
};

inline bool implicit_ode_advan(int advan) {
  return advan == 8 || advan == 9 || advan == 13 || advan == 14;
}

inline std::string ode_kernel_name(int advan) {
  switch (advan) {
    case 8: return "advan8-stiff-implicit";
    case 9: return "advan9-dae-implicit";
    case 10: return "advan10-michaelis-menten-rk45";
    case 13: return "advan13-implicit";
    case 14: return "advan14-stiff-nonstiff-implicit";
    default: return "advan6-rk45";
  }
}

bool starts_with(const std::string& value, const char* prefix) {
  return value.rfind(prefix, 0) == 0;
}

int indexed_name(const std::string& name, const char* prefix) {
  if (!starts_with(name, prefix)) return -1;
  try {
    int index = std::stoi(name.substr(std::string(prefix).size()));
    return index - 1;
  } catch (...) {
    return -1;
  }
}

double data_value(const Rcpp::DataFrame& data, const std::string& name, int row) {
  if (!data.containsElementNamed(name.c_str())) {
    throw std::invalid_argument("PRED input '" + name + "' is not present in the dataset.");
  }
  Rcpp::RObject object = data[name];
  if (TYPEOF(object) == REALSXP) return Rcpp::NumericVector(object)[row];
  if (TYPEOF(object) == INTSXP || TYPEOF(object) == LGLSXP) {
    int value = Rcpp::IntegerVector(object)[row];
    return value == NA_INTEGER ? NA_REAL : static_cast<double>(value);
  }
  throw std::invalid_argument("PRED input '" + name + "' must be numeric.");
}

int eta_column(const ModelEngine& engine, const Rcpp::DataFrame& data,
               int row, int eta_index, int eta_columns) {
  if (eta_index < 0 || eta_index >= engine.n_eta) {
    throw std::out_of_range("ETA index exceeds the model ETA definitions.");
  }
  if (engine.re_enabled) {
    const std::string column_name = ".ETA_COLUMN_" + std::to_string(eta_index + 1);
    if (!data.containsElementNamed(column_name.c_str())) {
      throw std::invalid_argument("General random-effect execution requires compiled ETA mapping columns.");
    }
    const int column = static_cast<int>(data_value(data, column_name, row)) - 1;
    if (column < 0 || column >= eta_columns) {
      throw std::out_of_range("Mapped random-effect ETA index exceeds the supplied ETA matrix.");
    }
    return column;
  }
  if (engine.iov <= 0 || eta_index < engine.n_eta - engine.iov) return eta_index;
  if (!data.containsElementNamed(".OCC_INDEX")) {
    throw std::invalid_argument("IOV execution requires compiled .OCC_INDEX data.");
  }
  const int between = engine.n_eta - engine.iov;
  const int occasion = static_cast<int>(data_value(data, ".OCC_INDEX", row)) - 1;
  const int column = between + occasion * engine.iov + (eta_index - between);
  if (occasion < 0 || column < 0 || column >= eta_columns) {
    throw std::out_of_range("Occasion-specific ETA index exceeds the supplied ETA matrix.");
  }
  return column;
}

int required_eta_columns(const ModelEngine& engine,
                         const Rcpp::DataFrame& data) {
  if (engine.re_enabled) {
    int maximum = 0;
    for (int eta = 1; eta <= engine.n_eta; ++eta) {
      const std::string name = ".ETA_COLUMN_" + std::to_string(eta);
      if (!data.containsElementNamed(name.c_str())) {
        throw std::invalid_argument("Compiled random-effect ETA mapping is missing.");
      }
      Rcpp::IntegerVector values(data[name]);
      for (int value : values) maximum = std::max(maximum, value);
    }
    return maximum;
  }
  if (engine.iov <= 0) return engine.n_eta;
  if (!data.containsElementNamed(".OCC_INDEX")) {
    throw std::invalid_argument("IOV execution requires compiled .OCC_INDEX data.");
  }
  Rcpp::IntegerVector occasion(data[".OCC_INDEX"]);
  int count = 0;
  for (int value : occasion) count = std::max(count, value);
  return engine.n_eta - engine.iov + count * engine.iov;
}

Parameters evaluate_parameters(const ModelEngine& engine,
                               const Rcpp::DataFrame& data,
                               int row, int subject,
                               const Rcpp::NumericVector& theta,
                               const Rcpp::NumericMatrix& eta,
                               const Rcpp::NumericVector& sigma) {
  std::vector<double> inputs(engine.pred->input_names.size(), 0.0);
  for (std::size_t i = 0; i < engine.pred->input_names.size(); ++i) {
    const std::string& name = engine.pred->input_names[i];
    int index = indexed_name(name, "THETA_");
    if (index >= 0) {
      if (index >= theta.size()) throw std::out_of_range("THETA index exceeds supplied values.");
      inputs[i] = theta[index];
      continue;
    }
    index = indexed_name(name, "ETA_");
    if (index >= 0) {
      inputs[i] = eta(subject, eta_column(engine, data, row, index, eta.ncol()));
      continue;
    }
    index = indexed_name(name, "SIGMA_");
    if (index >= 0) {
      if (index >= sigma.size()) throw std::out_of_range("SIGMA index exceeds supplied values.");
      inputs[i] = sigma[index];
      continue;
    }
    if (starts_with(name, "ERR_")) {
      inputs[i] = 0.0;
      continue;
    }
    if (name == "F") {
      inputs[i] = 0.0;
      continue;
    }
    if (name == "MIXNUM") {
      inputs[i] = data.containsElementNamed("MIXNUM") ?
        data_value(data, "MIXNUM", row) : 1.0;
      continue;
    }
    inputs[i] = data_value(data, name, row);
    if (!std::isfinite(inputs[i])) {
      throw std::domain_error("PRED input '" + name + "' is non-finite at row " +
                              std::to_string(row + 1) + ".");
    }
  }
  std::vector<double> output = engine.pred->eval_outputs(inputs, engine.all_outputs);
  Parameters parameters;
  for (std::size_t i = 0; i < output.size(); ++i) {
    parameters[engine.pred->output_names[i]] = output[i];
  }
  return parameters;
}

double evaluate_post_prediction(
    const ModelEngine& engine, const Rcpp::DataFrame& data,
    int row, int subject, double time, const Vector& state,
    const Rcpp::NumericVector& theta, const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma, double advan_prediction,
    Parameters& parameters) {
  if (!engine.post_pred) return advan_prediction;
  std::vector<double> inputs(engine.post_pred->input_names.size(), 0.0);
  for (std::size_t i = 0; i < engine.post_pred->input_names.size(); ++i) {
    const std::string& name = engine.post_pred->input_names[i];
    int index = indexed_name(name, "THETA_");
    if (index >= 0) {
      if (index >= theta.size()) throw std::out_of_range("THETA index exceeds values.");
      inputs[i] = theta[index];
      continue;
    }
    index = indexed_name(name, "ETA_");
    if (index >= 0) {
      inputs[i] = eta(subject, eta_column(engine, data, row, index, eta.ncol()));
      continue;
    }
    index = indexed_name(name, "SIGMA_");
    if (index >= 0) {
      if (index >= sigma.size()) throw std::out_of_range("SIGMA index exceeds values.");
      inputs[i] = sigma[index];
      continue;
    }
    index = indexed_name(name, "A_");
    if (index >= 0) {
      if (index >= state.size()) throw std::out_of_range("$PRED A() index exceeds state dimension.");
      inputs[i] = state[index];
      continue;
    }
    if (name == "F_ADVAN") { inputs[i] = advan_prediction; continue; }
    if (name == "T" || name == "TIME") { inputs[i] = time; continue; }
    if (name == "MIXNUM") {
      inputs[i] = data.containsElementNamed("MIXNUM") ?
        data_value(data, "MIXNUM", row) : 1.0;
      continue;
    }
    const auto assigned = parameters.find(name);
    if (assigned != parameters.end()) {
      inputs[i] = assigned->second;
      continue;
    }
    inputs[i] = data_value(data, name, row);
    if (!std::isfinite(inputs[i])) {
      throw std::domain_error(
        "Post-ADVAN $PRED input '" + name + "' is non-finite at row " +
        std::to_string(row + 1) + ".");
    }
  }
  const std::vector<double> output =
    engine.post_pred->eval_outputs(inputs, engine.post_all_outputs);
  for (std::size_t i = 0; i < output.size(); ++i) {
    parameters[engine.post_pred->output_names[i]] = output[i];
  }
  const auto prediction = parameters.find("F");
  if (prediction == parameters.end() || !std::isfinite(prediction->second)) {
    throw std::domain_error("Post-ADVAN $PRED did not produce a finite F.");
  }
  return prediction->second;
}

Vector evaluate_algebraic_residuals(
    const ModelEngine& engine, const Rcpp::DataFrame& data,
    int row, int subject, double time, const Vector& state,
    const Parameters& parameters, const Rcpp::NumericVector& theta,
    const Rcpp::NumericMatrix& eta, const Rcpp::NumericVector& sigma,
    const Vector& algebraic) {
  if (!engine.alg) throw std::logic_error("DAE algebraic residual program is missing.");
  std::vector<double> inputs(engine.alg->input_names.size(), 0.0);
  for (std::size_t i = 0; i < engine.alg->input_names.size(); ++i) {
    const std::string& name = engine.alg->input_names[i];
    int index = indexed_name(name, "A_");
    if (index >= 0) {
      if (index >= state.size()) throw std::out_of_range("A() index exceeds the DAE state dimension.");
      inputs[i] = state[index];
      continue;
    }
    if (name == "T") { inputs[i] = time; continue; }
    auto variable = std::find(engine.dae_variables.begin(), engine.dae_variables.end(), name);
    if (variable != engine.dae_variables.end()) {
      inputs[i] = algebraic[std::distance(engine.dae_variables.begin(), variable)];
      continue;
    }
    auto parameter = parameters.find(name);
    if (parameter != parameters.end()) { inputs[i] = parameter->second; continue; }
    index = indexed_name(name, "THETA_");
    if (index >= 0) { inputs[i] = theta[index]; continue; }
    index = indexed_name(name, "ETA_");
    if (index >= 0) {
      inputs[i] = eta(subject, eta_column(engine, data, row, index, eta.ncol()));
      continue;
    }
    index = indexed_name(name, "SIGMA_");
    if (index >= 0) { inputs[i] = sigma[index]; continue; }
    if (starts_with(name, "ERR_") || name == "F") continue;
    inputs[i] = data_value(data, name, row);
  }
  const std::vector<double> output = engine.alg->eval_outputs(inputs, engine.algebraic_outputs);
  Vector residual(static_cast<Eigen::Index>(output.size()));
  for (std::size_t i = 0; i < output.size(); ++i) residual[static_cast<Eigen::Index>(i)] = output[i];
  return residual;
}

Vector solve_algebraic(
    const ModelEngine& engine, const Rcpp::DataFrame& data,
    int row, int subject, double time, const Vector& state,
    const Parameters& parameters, const Rcpp::NumericVector& theta,
    const Rcpp::NumericMatrix& eta, const Rcpp::NumericVector& sigma) {
  Vector value(static_cast<Eigen::Index>(engine.dae_initial.size()));
  for (std::size_t i = 0; i < engine.dae_initial.size(); ++i) value[static_cast<Eigen::Index>(i)] = engine.dae_initial[i];
  Vector residual;
  for (int iteration = 0; iteration < engine.dae_maxit; ++iteration) {
    residual = evaluate_algebraic_residuals(
      engine, data, row, subject, time, state, parameters, theta, eta, sigma, value);
    if (!residual.allFinite()) throw std::domain_error("DAE residual is non-finite.");
    if (residual.cwiseAbs().maxCoeff() <= engine.dae_tolerance) return value;
    Matrix jacobian = Matrix::Zero(value.size(), value.size());
    for (Eigen::Index column = 0; column < value.size(); ++column) {
      const double delta = engine.dae_jacobian_step * std::max(1.0, std::abs(value[column]));
      Vector plus = value; plus[column] += delta;
      Vector minus = value; minus[column] -= delta;
      const Vector upper = evaluate_algebraic_residuals(
        engine, data, row, subject, time, state, parameters, theta, eta, sigma, plus);
      const Vector lower = evaluate_algebraic_residuals(
        engine, data, row, subject, time, state, parameters, theta, eta, sigma, minus);
      jacobian.col(column) = (upper - lower) / (2.0 * delta);
      if (!engine.dae_sparsity.empty()) {
        for (Eigen::Index row_index = 0; row_index < value.size(); ++row_index) {
          if (!engine.dae_sparsity[static_cast<std::size_t>(row_index * value.size() + column)]) {
            jacobian(row_index, column) = 0.0;
          }
        }
      }
    }
    Vector update = Vector::Zero(value.size());
    for (std::size_t block = 0; block < engine.dae_block_rows.size(); ++block) {
      const auto& rows = engine.dae_block_rows[block];
      const auto& columns = engine.dae_block_columns[block];
      Matrix local(rows.size(), columns.size());
      Vector rhs(rows.size());
      for (std::size_t local_row = 0; local_row < rows.size(); ++local_row) {
        rhs[static_cast<Eigen::Index>(local_row)] = -residual[rows[local_row]];
        for (std::size_t local_column = 0; local_column < columns.size(); ++local_column) {
          local(static_cast<Eigen::Index>(local_row),
                static_cast<Eigen::Index>(local_column)) =
            jacobian(rows[local_row], columns[local_column]);
        }
      }
      Eigen::FullPivLU<Matrix> lu(local);
      if (!lu.isInvertible()) throw std::runtime_error("DAE Newton Jacobian block is singular.");
      const Vector local_update = lu.solve(rhs);
      for (std::size_t local_column = 0; local_column < columns.size(); ++local_column) {
        update[columns[local_column]] = local_update[static_cast<Eigen::Index>(local_column)];
      }
    }
    if (!update.allFinite()) throw std::runtime_error("DAE Newton update is non-finite.");
    value += update;
    if (update.cwiseAbs().maxCoeff() <= engine.dae_tolerance) {
      residual = evaluate_algebraic_residuals(
        engine, data, row, subject, time, state, parameters, theta, eta, sigma, value);
      if (residual.cwiseAbs().maxCoeff() <= 10.0 * engine.dae_tolerance) return value;
    }
  }
  throw std::runtime_error("DAE algebraic Newton solve did not converge.");
}
