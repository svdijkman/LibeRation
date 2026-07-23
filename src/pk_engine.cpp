// [[Rcpp::depends(LibeRtAD)]]
// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>
#include <LibeRtAD/eigen_r.hpp>
#include <LibeRtAD/sparse_hessian.hpp>
#include <unsupported/Eigen/MatrixFunctions>
#include <LibeRtAD/program.hpp>
#include "eigen_solver.h"

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <numeric>
#include <queue>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace liberation {

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

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
  std::shared_ptr<const libertad::Program> des;
  std::shared_ptr<const libertad::Program> alg;
  std::shared_ptr<const libertad::Program> error;
  std::vector<std::size_t> all_outputs;
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
        solver(Rcpp::as<std::string>(spec["solver"])),
        error_type(Rcpp::as<std::string>(spec["error_type"])),
        pred(std::make_shared<const libertad::Program>(Rcpp::as<Rcpp::List>(spec["pred_ir"]))) {
    all_outputs.resize(pred->output_names.size());
    std::iota(all_outputs.begin(), all_outputs.end(), 0U);
    if (spec.containsElementNamed("output_names")) {
      selected_output_names = Rcpp::as<std::vector<std::string>>(spec["output_names"]);
      // Validate the serialized selection at engine construction rather than
      // failing part-way through a long estimation or simulation.
      pred->select_outputs(selected_output_names);
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

struct DdeHistory {
  std::vector<double> time;
  std::vector<Vector> state;
  std::vector<double> baseline;

  void reset(double at, const Vector& value, const std::vector<double>& history) {
    time.assign(1U, at); state.assign(1U, value); baseline = history;
  }
  void append(double at, const Vector& value) {
    if (!time.empty() && std::abs(time.back() - at) <= 1e-12) {
      state.back() = value; return;
    }
    time.push_back(at); state.push_back(value);
  }
  double at(double target, int component) const {
    if (component < 0 || component >= static_cast<int>(baseline.size())) {
      throw std::out_of_range("DDE lag state is outside the state vector.");
    }
    if (time.empty() || target < time.front() - 1e-12) return baseline[static_cast<std::size_t>(component)];
    if (target >= time.back() - 1e-12) return state.back()[component];
    auto upper = std::upper_bound(time.begin(), time.end(), target);
    const std::size_t right = static_cast<std::size_t>(std::distance(time.begin(), upper));
    const std::size_t left = right - 1U;
    const double fraction = (target - time[left]) / (time[right] - time[left]);
    return state[left][component] + fraction * (state[right][component] - state[left][component]);
  }
};

Vector evaluate_derivatives(const ModelEngine& engine,
                            const Rcpp::DataFrame& data,
                            int row, int subject, double t,
                            const Vector& state,
                            const Parameters& parameters,
                            const Rcpp::NumericVector& theta,
                            const Rcpp::NumericMatrix& eta,
                            const Rcpp::NumericVector& sigma,
                            const Vector* lag_values = nullptr) {
  if (!engine.des) throw std::logic_error("ODE derivative program is missing.");
  const Vector algebraic = engine.dae_enabled ? solve_algebraic(
    engine, data, row, subject, t, state, parameters, theta, eta, sigma) : Vector();
  std::vector<double> inputs(engine.des->input_names.size(), 0.0);
  for (std::size_t i = 0; i < engine.des->input_names.size(); ++i) {
    const std::string& name = engine.des->input_names[i];
    int index = indexed_name(name, "A_");
    if (index >= 0) {
      if (index >= state.size()) throw std::out_of_range("A() index exceeds the ODE state dimension.");
      inputs[i] = state[index];
      continue;
    }
    if (name == "T") {
      inputs[i] = t;
      continue;
    }
    auto lag = std::find(engine.dde_lag_inputs.begin(), engine.dde_lag_inputs.end(), name);
    if (lag != engine.dde_lag_inputs.end()) {
      if (lag_values == nullptr) throw std::logic_error("DDE lag history is unavailable.");
      inputs[i] = (*lag_values)[std::distance(engine.dde_lag_inputs.begin(), lag)];
      continue;
    }
    auto variable = std::find(engine.dae_variables.begin(), engine.dae_variables.end(), name);
    if (variable != engine.dae_variables.end()) {
      inputs[i] = algebraic[std::distance(engine.dae_variables.begin(), variable)];
      continue;
    }
    auto parameter = parameters.find(name);
    if (parameter != parameters.end()) {
      inputs[i] = parameter->second;
      continue;
    }
    index = indexed_name(name, "THETA_");
    if (index >= 0) {
      if (index >= theta.size()) throw std::out_of_range("THETA index exceeds supplied values in DES.");
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
      if (index >= sigma.size()) throw std::out_of_range("SIGMA index exceeds supplied values in DES.");
      inputs[i] = sigma[index];
      continue;
    }
    if (starts_with(name, "ERR_") || name == "F") {
      inputs[i] = 0.0;
      continue;
    }
    inputs[i] = data_value(data, name, row);
    if (!std::isfinite(inputs[i])) {
      throw std::domain_error("DES input '" + name + "' is non-finite at row " +
                              std::to_string(row + 1) + ".");
    }
  }
  std::vector<double> values = engine.des->eval_outputs(inputs, engine.derivative_outputs);
  Vector derivative(static_cast<Eigen::Index>(values.size()));
  for (std::size_t i = 0; i < values.size(); ++i) derivative[static_cast<Eigen::Index>(i)] = values[i];
  if (!derivative.allFinite()) {
    throw std::domain_error("DES produced a non-finite derivative at time " + std::to_string(t) + ".");
  }
  return derivative;
}

using OdeRhs = std::function<Vector(double, const Vector&)>;

double scaled_error(const Vector& error, const Vector& before,
                    const Vector& after, const OdeControl& control) {
  double maximum = 0.0;
  for (Eigen::Index i = 0; i < error.size(); ++i) {
    const double scale = control.atol + control.rtol *
      std::max(std::abs(before[i]), std::abs(after[i]));
    maximum = std::max(maximum, std::abs(error[i]) / scale);
  }
  return maximum;
}

Vector integrate_dopri54(const OdeRhs& rhs, Vector state, double from, double to,
                         const OdeControl& control) {
  if (to <= from) return state;
  double t = from;
  const double span = to - from;
  double h = control.initial_step > 0.0 ? std::min(control.initial_step, span) : span / 10.0;
  h = std::max(h, std::min(span, 1e-8));
  int attempts = 0;
  while (t < to) {
    if (++attempts > control.max_steps) {
      throw std::runtime_error("ADVAN6 exceeded ODE_CONTROL$max_steps.");
    }
    h = std::min(h, to - t);
    const double minimum = 32.0 * std::numeric_limits<double>::epsilon() *
      std::max({1.0, std::abs(t), std::abs(to)});
    if (h < minimum) throw std::runtime_error("ADVAN6 ODE step size underflow.");

    const Vector k1 = rhs(t, state);
    const Vector k2 = rhs(t + h * (1.0 / 5.0),
      state + h * ((1.0 / 5.0) * k1));
    const Vector k3 = rhs(t + h * (3.0 / 10.0),
      state + h * ((3.0 / 40.0) * k1 + (9.0 / 40.0) * k2));
    const Vector k4 = rhs(t + h * (4.0 / 5.0),
      state + h * ((44.0 / 45.0) * k1 - (56.0 / 15.0) * k2 + (32.0 / 9.0) * k3));
    const Vector k5 = rhs(t + h * (8.0 / 9.0),
      state + h * ((19372.0 / 6561.0) * k1 - (25360.0 / 2187.0) * k2 +
                   (64448.0 / 6561.0) * k3 - (212.0 / 729.0) * k4));
    const Vector k6 = rhs(t + h,
      state + h * ((9017.0 / 3168.0) * k1 - (355.0 / 33.0) * k2 +
                   (46732.0 / 5247.0) * k3 + (49.0 / 176.0) * k4 -
                   (5103.0 / 18656.0) * k5));
    const Vector fifth = state + h * ((35.0 / 384.0) * k1 +
      (500.0 / 1113.0) * k3 + (125.0 / 192.0) * k4 -
      (2187.0 / 6784.0) * k5 + (11.0 / 84.0) * k6);
    const Vector k7 = rhs(t + h, fifth);
    const Vector fourth = state + h * ((5179.0 / 57600.0) * k1 +
      (7571.0 / 16695.0) * k3 + (393.0 / 640.0) * k4 -
      (92097.0 / 339200.0) * k5 + (187.0 / 2100.0) * k6 + (1.0 / 40.0) * k7);
    const double error = scaled_error(fifth - fourth, state, fifth, control);
    if (!std::isfinite(error)) throw std::domain_error("ADVAN6 ODE error estimate is non-finite.");
    if (error <= 1.0) {
      state = fifth;
      t += h;
    }
    const double factor = error == 0.0 ? 5.0 :
      std::clamp(0.9 * std::pow(error, -0.2), 0.1, 5.0);
    h *= factor;
  }
  return state;
}

bool implicit_trapezoid_step(const OdeRhs& rhs, const Vector& before,
                             double t, double h, const OdeControl& control,
                             Vector& after) {
  const Vector f0 = rhs(t, before);
  after = before + h * f0;
  const Eigen::Index n = before.size();
  for (int iteration = 0; iteration < 12; ++iteration) {
    const Vector f1 = rhs(t + h, after);
    const Vector residual = after - before - 0.5 * h * (f0 + f1);
    if (scaled_error(residual, before, after, control) < 0.03) return after.allFinite();
    Matrix jacobian(n, n);
    for (Eigen::Index j = 0; j < n; ++j) {
      Vector perturbed = after;
      const double delta = std::sqrt(std::numeric_limits<double>::epsilon()) *
        std::max(1.0, std::abs(after[j]));
      perturbed[j] += delta;
      jacobian.col(j) = (rhs(t + h, perturbed) - f1) / delta;
    }
    Matrix system = Matrix::Identity(n, n) - 0.5 * h * jacobian;
    Eigen::FullPivLU<Matrix> lu(system);
    if (!lu.isInvertible()) return false;
    Vector update = lu.solve(-residual);
    if (!update.allFinite()) return false;
    after += update;
    if (scaled_error(update, before, after, control) < 0.03) return after.allFinite();
  }
  return false;
}

Vector integrate_implicit_trapezoid(const OdeRhs& rhs, Vector state,
                                    double from, double to,
                                    const OdeControl& control) {
  if (to <= from) return state;
  double t = from;
  const double span = to - from;
  double h = control.initial_step > 0.0 ? std::min(control.initial_step, span) : span / 10.0;
  h = std::max(h, std::min(span, 1e-8));
  int attempts = 0;
  while (t < to) {
    if (++attempts > control.max_steps) {
      throw std::runtime_error("ADVAN13 exceeded ODE_CONTROL$max_steps.");
    }
    h = std::min(h, to - t);
    const double minimum = 32.0 * std::numeric_limits<double>::epsilon() *
      std::max({1.0, std::abs(t), std::abs(to)});
    if (h < minimum) throw std::runtime_error("ADVAN13 ODE step size underflow.");

    Vector full, half, two_half;
    const bool converged = implicit_trapezoid_step(rhs, state, t, h, control, full) &&
      implicit_trapezoid_step(rhs, state, t, h * 0.5, control, half) &&
      implicit_trapezoid_step(rhs, half, t + h * 0.5, h * 0.5, control, two_half);
    double error = std::numeric_limits<double>::infinity();
    if (converged) error = scaled_error((two_half - full) / 3.0, state, two_half, control);
    if (converged && std::isfinite(error) && error <= 1.0) {
      state = two_half + (two_half - full) / 3.0;
      t += h;
    }
    const double factor = converged && error == 0.0 ? 4.0 :
      (converged && std::isfinite(error) ?
        std::clamp(0.9 * std::pow(error, -1.0 / 3.0), 0.1, 4.0) : 0.25);
    h *= factor;
  }
  return state;
}

int compartment_index(int cmt, int fallback, int n) {
  int index = cmt > 0 ? cmt - 1 : fallback;
  if (index < 0 || index >= n) {
    throw std::out_of_range("Event compartment is outside the model state vector.");
  }
  return index;
}

double row_optional(const Rcpp::DataFrame& data, const std::string& name,
                    int row, double fallback) {
  if (!data.containsElementNamed(name.c_str())) return fallback;
  double value = data_value(data, name, row);
  return std::isfinite(value) ? value : fallback;
}

double bioavailability(const Parameters& p, const Rcpp::DataFrame& data,
                       int row, int cmt) {
  const std::string name = "F" + std::to_string(cmt);
  double value = get_positive(p, {name.c_str()});
  if (!finite_positive(value)) value = row_optional(data, name, row, 1.0);
  return finite_positive(value) ? value : 1.0;
}

double event_infusion_rate(const Parameters& p, const Rcpp::DataFrame& data,
                           int row, int cmt, double amount, double rate_code) {
  if (rate_code >= 0.0) return rate_code;
  if (rate_code != -1.0 && rate_code != -2.0) {
    throw std::domain_error("Negative RATE must be -1 (modelled Rn) or -2 (modelled Dn).");
  }
  const std::string name = std::string(rate_code == -1.0 ? "R" : "D") +
    std::to_string(cmt);
  double value = get_positive(p, {name.c_str()});
  if (!finite_positive(value)) value = row_optional(data, name, row, NA_REAL);
  if (!finite_positive(value)) {
    throw std::domain_error("RATE=" + std::to_string(static_cast<int>(rate_code)) +
                            " requires a positive " + name + " value.");
  }
  return rate_code == -1.0 ? value : amount / value;
}

double observation_scale(const Parameters& p, const Rcpp::DataFrame& data,
                         int row, int cmt, const Topology& topology) {
  const std::string name = "S" + std::to_string(cmt);
  double value = get_positive(p, {name.c_str()});
  if (!finite_positive(value)) value = row_optional(data, name, row, NA_REAL);
  const int index = cmt - 1;
  if (!finite_positive(value) && index >= 0 && index < static_cast<int>(topology.default_scales.size())) {
    value = topology.default_scales[static_cast<std::size_t>(index)];
  }
  return finite_positive(value) ? value : 1.0;
}

Vector infusion_input(int n, const std::vector<ActiveInfusion>& active) {
  Vector input = Vector::Zero(n);
  for (const ActiveInfusion& infusion : active) {
    input[infusion.compartment] += infusion.rate;
  }
  return input;
}

void remove_finished(std::vector<ActiveInfusion>& active, double time) {
  active.erase(
    std::remove_if(active.begin(), active.end(),
      [time](const ActiveInfusion& infusion) { return infusion.end <= time + 1e-12; }),
    active.end()
  );
}

Vector propagate_to(const Matrix& k, Vector state, double from, double to,
                    std::vector<ActiveInfusion>& active) {
  double cursor = from;
  remove_finished(active, cursor);
  while (cursor < to - 1e-12) {
    double segment_end = to;
    for (const ActiveInfusion& infusion : active) {
      if (infusion.end > cursor + 1e-12) segment_end = std::min(segment_end, infusion.end);
    }
    state = propagate(k, infusion_input(k.rows(), active), segment_end - cursor, state);
    cursor = segment_end;
    remove_finished(active, cursor);
  }
  return state;
}

Vector propagate_ode_to(const ModelEngine& engine,
                        const Rcpp::DataFrame& data,
                        int row, int subject,
                        const Rcpp::NumericVector& theta,
                        const Rcpp::NumericMatrix& eta,
                        const Rcpp::NumericVector& sigma,
                        const Parameters& parameters,
                        Vector state, double from, double to,
                        std::vector<ActiveInfusion>& active,
                        DdeHistory* dde_history = nullptr) {
  double cursor = from;
  remove_finished(active, cursor);
  while (cursor < to - 1e-12) {
    double segment_end = to;
    for (const ActiveInfusion& infusion : active) {
      if (infusion.end > cursor + 1e-12) segment_end = std::min(segment_end, infusion.end);
    }
    const Vector input = infusion_input(engine.n_state, active);
    if (engine.dde_enabled) {
      if (dde_history == nullptr) {
        throw std::logic_error("DDE propagation requires an initialized history.");
      }
      double time = cursor;
      int steps = 0;
      while (time < segment_end - 1e-12) {
        if (++steps > engine.dde_max_steps) {
          throw std::runtime_error("DDE method-of-steps exceeded max_steps.");
        }
        const double h = std::min(engine.dde_step, segment_end - time);
        auto rhs = [&](double stage_time, const Vector& value) {
          Vector lag_values(static_cast<Eigen::Index>(engine.dde_lag_inputs.size()));
          for (std::size_t lag = 0; lag < engine.dde_lag_inputs.size(); ++lag) {
            auto delay = parameters.find(engine.dde_lag_delays[lag]);
            if (delay == parameters.end() || !std::isfinite(delay->second) ||
                delay->second < engine.dde_minimum_delay || delay->second < h) {
              throw std::domain_error("DDE delay '" + engine.dde_lag_delays[lag] +
                                      "' must be finite and at least the integration step.");
            }
            lag_values[static_cast<Eigen::Index>(lag)] = dde_history->at(
              stage_time - delay->second, engine.dde_lag_states[lag]);
          }
          Vector derivative = evaluate_derivatives(
            engine, data, row, subject, stage_time, value, parameters,
            theta, eta, sigma, &lag_values);
          derivative += input;
          return derivative;
        };
        const Vector k1 = rhs(time, state);
        const Vector k2 = rhs(time + 0.5 * h, state + 0.5 * h * k1);
        const Vector k3 = rhs(time + 0.5 * h, state + 0.5 * h * k2);
        const Vector k4 = rhs(time + h, state + h * k3);
        state += (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        time += h;
        dde_history->append(time, state);
      }
      cursor = segment_end;
      remove_finished(active, cursor);
      continue;
    }
    OdeRhs rhs = [&](double t, const Vector& y) {
      Vector derivative = evaluate_derivatives(
        engine, data, row, subject, t, y, parameters, theta, eta, sigma
      );
      derivative += input;
      return derivative;
    };
    state = engine.advan == 13 ?
      integrate_implicit_trapezoid(rhs, state, cursor, segment_end, engine.ode_control) :
      integrate_dopri54(rhs, state, cursor, segment_end, engine.ode_control);
    cursor = segment_end;
    remove_finished(active, cursor);
  }
  return state;
}

Vector steady_bolus_post(const Matrix& k, const Vector& dose, double interval) {
  if (!(interval > 0.0) || !std::isfinite(interval)) {
    throw std::domain_error("Steady-state bolus requires II > 0.");
  }
  Matrix transition = matrix_exp(k * interval);
  return solve_periodic(transition, dose, "Bolus");
}

Vector steady_infusion_pre(const Matrix& k, const Vector& rate,
                           double duration, double interval) {
  if (!(duration > 0.0) || !(interval > 0.0) || !std::isfinite(duration) ||
      !std::isfinite(interval)) {
    throw std::domain_error("Steady-state infusion requires finite duration and II > 0.");
  }
  const int complete = static_cast<int>(std::floor(duration / interval + 1e-12));
  double remainder = duration - complete * interval;
  if (remainder < 1e-12) remainder = 0.0;
  const Vector baseline = static_cast<double>(complete) * rate;
  AffineMap first = affine_map(k, baseline + (remainder > 0.0 ? rate : Vector::Zero(k.rows())),
                               remainder);
  AffineMap second = affine_map(k, baseline, interval - remainder);
  Matrix period_transition = second.transition * first.transition;
  Vector period_offset = second.transition * first.offset + second.offset;
  return solve_periodic(period_transition, period_offset, "Infusion");
}

double relative_state_change(const Vector& before, const Vector& after) {
  return (after - before).norm() / std::max(1.0, after.norm());
}

std::vector<ActiveInfusion> periodic_infusions(double time, double duration,
                                               double interval, int compartment,
                                               double rate) {
  const int previous = std::max(0, static_cast<int>(std::ceil(duration / interval - 1e-12)) - 1);
  std::vector<ActiveInfusion> active;
  active.reserve(static_cast<std::size_t>(previous + 1));
  for (int dose = 0; dose <= previous; ++dose) {
    const double end = time + duration - dose * interval;
    if (end > time + 1e-12) active.push_back({end, compartment, rate});
  }
  return active;
}

Vector steady_ode_bolus_post(
    const ModelEngine& engine, const Rcpp::DataFrame& data,
    int row, int subject, const Rcpp::NumericVector& theta,
    const Rcpp::NumericMatrix& eta, const Rcpp::NumericVector& sigma,
    const Parameters& parameters, const Vector& dose,
    double time, double interval) {
  if (!(interval > 0.0)) throw std::domain_error("ODE steady-state bolus requires II > 0.");
  Vector current = dose;
  const double tolerance = std::max(1e-10, engine.ode_control.rtol * 5.0);
  for (int iteration = 0; iteration < 10000; ++iteration) {
    std::vector<ActiveInfusion> active;
    Vector next = propagate_ode_to(
      engine, data, row, subject, theta, eta, sigma, parameters,
      current, time, time + interval, active) + dose;
    if (relative_state_change(current, next) <= tolerance) return next;
    current = next;
  }
  throw std::runtime_error("ODE bolus periodic shooting did not converge.");
}

Vector steady_ode_infusion_pre(
    const ModelEngine& engine, const Rcpp::DataFrame& data,
    int row, int subject, const Rcpp::NumericVector& theta,
    const Rcpp::NumericMatrix& eta, const Rcpp::NumericVector& sigma,
    const Parameters& parameters, int compartment, double administered_rate,
    double duration, double time, double interval) {
  if (!(duration > 0.0) || !(interval > 0.0)) {
    throw std::domain_error("ODE steady-state infusion requires duration and II > 0.");
  }
  Vector current = Vector::Zero(engine.n_state);
  const double tolerance = std::max(1e-10, engine.ode_control.rtol * 5.0);
  for (int iteration = 0; iteration < 10000; ++iteration) {
    std::vector<ActiveInfusion> active = periodic_infusions(
      time, duration, interval, compartment, administered_rate);
    Vector next = propagate_ode_to(
      engine, data, row, subject, theta, eta, sigma, parameters,
      current, time, time + interval, active);
    if (relative_state_change(current, next) <= tolerance) return next;
    current = next;
  }
  throw std::runtime_error("ODE infusion periodic shooting did not converge.");
}

Rcpp::List simulate(ModelEngine& engine,
                    const Rcpp::DataFrame& data,
                    const Rcpp::NumericVector& theta,
                    const Rcpp::NumericMatrix& eta,
                    const Rcpp::NumericVector& sigma) {
  const int n_rows = data.nrows();
  if (theta.size() != engine.n_theta) Rcpp::stop("Theta vector has the wrong length.");
  const int minimum_eta_columns = required_eta_columns(engine, data);
  const int between_eta = engine.n_eta - engine.iov;
  if (eta.ncol() < minimum_eta_columns ||
      (!engine.re_enabled && engine.iov > 0 &&
       (eta.ncol() - between_eta) % engine.iov != 0)) {
    Rcpp::stop("ETA matrix has the wrong number of between-subject/occasion columns.");
  }

  Rcpp::NumericVector time = data["TIME"];
  Rcpp::NumericVector amount = data["AMT"];
  Rcpp::NumericVector rate = data["RATE"];
  Rcpp::NumericVector interval = data["II"];
  Rcpp::IntegerVector evid = data["EVID"];
  Rcpp::IntegerVector cmt = data["CMT"];
  Rcpp::IntegerVector ss = data["SS"];
  Rcpp::IntegerVector eta_subject_index = data[".ID_INDEX"];
  Rcpp::IntegerVector subject_index = data.containsElementNamed(".STRUCT_ID_INDEX") ?
    Rcpp::IntegerVector(data[".STRUCT_ID_INDEX"]) : eta_subject_index;
  int n_subjects = 0;
  for (int value : eta_subject_index) n_subjects = std::max(n_subjects, value);
  if (eta.nrow() != n_subjects) Rcpp::stop("ETA matrix has the wrong number of subject rows.");

  const int n_state = engine.n_state;
  if (n_state < 1) Rcpp::stop("Model state dimension must be positive.");
  Rcpp::NumericVector prediction(n_rows, NA_REAL);
  Rcpp::NumericMatrix amounts(n_rows, n_state);
  Rcpp::NumericMatrix generated(n_rows, engine.selected_output_names.size());
  Vector state = Vector::Zero(n_state);
  std::vector<ActiveInfusion> active;
  DdeHistory dde_history;
  Matrix previous_k = Matrix::Zero(n_state, n_state);
  Parameters previous_parameters;
  int previous_row = -1;
  bool have_previous = false;
  int previous_subject = -1;
  double previous_time = 0.0;
  std::vector<std::string> state_names;

  for (int row = 0; row < n_rows; ++row) {
    const int structural_subject = subject_index[row] - 1;
    const int subject = eta_subject_index[row] - 1;
    if (structural_subject != previous_subject) {
      state.setZero();
      active.clear();
      have_previous = false;
      previous_time = time[row];
      previous_subject = structural_subject;
      if (engine.dde_enabled) dde_history.reset(time[row], state, engine.dde_history);
    }
    if (have_previous) {
      if (time[row] < previous_time - 1e-12) Rcpp::stop("Subject event times are decreasing.");
      if (engine.is_ode()) {
        state = propagate_ode_to(
          engine, data, previous_row, subject, theta, eta, sigma,
          previous_parameters, state, previous_time, time[row], active
          , engine.dde_enabled ? &dde_history : nullptr
        );
      } else {
        state = propagate_to(previous_k, state, previous_time, time[row], active);
      }
    }

    Parameters parameters = evaluate_parameters(engine, data, row, subject, theta, eta, sigma);
    for (std::size_t output = 0; output < engine.selected_output_names.size(); ++output) {
      auto found = parameters.find(engine.selected_output_names[output]);
      generated(row, static_cast<int>(output)) = found == parameters.end() ? NA_REAL : found->second;
    }
    Topology topology;
    if (engine.is_ode()) {
      topology.k = Matrix::Zero(n_state, n_state);
      topology.state_names = engine.state_names;
      topology.default_scales.assign(static_cast<std::size_t>(n_state), 1.0);
    } else {
      topology = engine.matrix_graph.enabled ?
        build_graph_topology(engine.matrix_graph, parameters) :
        build_topology(engine.advan, parameters);
    }
    if (topology.k.rows() != n_state) Rcpp::stop("ADVAN state dimension changed unexpectedly.");
    state_names = topology.state_names;

    const bool reset = evid[row] == 3 || evid[row] == 4;
    if (reset) {
      state.setZero();
      active.clear();
      if (engine.dde_enabled) dde_history.reset(time[row], state, engine.dde_history);
    }
    const bool dosing = amount[row] > 0.0 && (evid[row] == 1 || evid[row] == 4 || evid[row] == 0);
    if (dosing) {
      const int dose_cmt = cmt[row] > 0 ? cmt[row] : engine.dose_cmp;
      const int dose_index = compartment_index(dose_cmt, topology.default_dose, n_state);
      const double f = bioavailability(parameters, data, row, dose_cmt);
      const double event_rate = event_infusion_rate(
        parameters, data, row, dose_cmt, amount[row], rate[row]);
      const int ss_flag = ss[row] != 0 ? ss[row] : engine.model_ss;
      if (ss_flag == 1) {
        state.setZero();
        active.clear();
      } else if (ss_flag != 0 && ss_flag != 2) {
        Rcpp::stop("Only SS=0, SS=1, and SS=2 are supported.");
      }
      if (engine.dde_enabled && ss_flag != 0) {
        Rcpp::stop("Experimental DDE models currently require SS=0; provide an explicit warm-up regimen.");
      }

      if (event_rate > 0.0) {
        const double duration = amount[row] / event_rate;
        Vector input = Vector::Zero(n_state);
        input[dose_index] = event_rate * f;
        if (ss_flag != 0) {
          Vector periodic = engine.is_ode() ? steady_ode_infusion_pre(
            engine, data, row, subject, theta, eta, sigma, parameters,
            dose_index, event_rate * f, duration, time[row], interval[row]) :
            steady_infusion_pre(topology.k, input, duration, interval[row]);
          if (ss_flag == 1) state = periodic;
          else state += periodic;
          std::vector<ActiveInfusion> periodic_active = periodic_infusions(
            time[row], duration, interval[row], dose_index, event_rate * f);
          active.insert(active.end(), periodic_active.begin(), periodic_active.end());
        } else {
          active.push_back({time[row] + duration, dose_index, event_rate * f});
        }
      } else {
        Vector dose = Vector::Zero(n_state);
        dose[dose_index] = amount[row] * f;
        if (ss_flag != 0) {
          Vector periodic = engine.is_ode() ? steady_ode_bolus_post(
            engine, data, row, subject, theta, eta, sigma, parameters,
            dose, time[row], interval[row]) :
            steady_bolus_post(topology.k, dose, interval[row]);
          if (ss_flag == 1) state = periodic;
          else state += periodic;
        } else {
          state += dose;
        }
      }
      if (engine.dde_enabled) dde_history.append(time[row], state);
    }

    const int observation_cmt = cmt[row] > 0 && evid[row] == 0 ? cmt[row] : engine.obs_cmp;
    const int observation_index = compartment_index(
      observation_cmt, topology.default_observation, n_state
    );
    const double scale = observation_scale(parameters, data, row, observation_cmt, topology);
    const auto direct_prediction = parameters.find("F");
    prediction[row] = direct_prediction == parameters.end() ?
      state[observation_index] / scale : direct_prediction->second;
    for (int j = 0; j < n_state; ++j) amounts(row, j) = state[j];

    previous_k = topology.k;
    previous_parameters = std::move(parameters);
    previous_row = row;
    previous_time = time[row];
    have_previous = true;
  }

  return Rcpp::List::create(
    Rcpp::Named("ipred") = prediction,
    Rcpp::Named("amounts") = amounts,
    Rcpp::Named("generated") = generated,
    Rcpp::Named("output_names") = Rcpp::wrap(engine.selected_output_names),
    Rcpp::Named("state_names") = state_names,
    Rcpp::Named("solver") = engine.solver == "auto" ?
      (engine.is_ode() ? (engine.advan == 13 ? "advan13-implicit" : "advan6-rk45") : "advan") :
      engine.solver
  );
}

// Scalar-generic analytical path used to record the *complete* event and
// compartment calculation with CppAD.  The ordinary simulation path above is
// deliberately retained as an independently testable double-precision
// implementation.  Keeping the two entry points separate also prevents an R
// callback or an adaptive solver decision from being hidden inside an AD
// recording.
template <class Scalar>
using MatrixT = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

template <class Scalar>
using VectorT = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

template <class Scalar>
using ParametersT = std::unordered_map<std::string, Scalar>;

template <class Scalar>
struct TopologyT {
  MatrixT<Scalar> k;
  int default_dose = 0;
  int default_observation = 0;
  std::vector<Scalar> default_scales;
};

template <class Scalar>
struct AffineMapT {
  MatrixT<Scalar> transition;
  VectorT<Scalar> offset;
};

template <class Scalar>
struct ActiveInfusionT {
  Scalar end = Scalar(0.0);
  int compartment = 0;
  Scalar rate = Scalar(0.0);
};

inline double scalar_value(double value) { return value; }

inline double scalar_value(const CppAD::AD<double>& value) {
  return CppAD::Value(CppAD::Var2Par(value));
}

template <class Scalar>
bool path_lt(const Scalar& left, const Scalar& right) { return left < right; }
template <class Scalar>
bool path_le(const Scalar& left, const Scalar& right) { return left <= right; }
template <class Scalar>
bool path_gt(const Scalar& left, const Scalar& right) { return left > right; }
template <class Scalar>
bool path_ne(const Scalar& left, const Scalar& right) { return left != right; }

template <class Scalar>
bool scalar_finite(const Scalar& value) {
  return std::isfinite(scalar_value(value));
}

template <class Scalar>
bool scalar_positive(const Scalar& value) {
  // Positivity is a parameter-domain check, not a valid alternate execution
  // path. Recording it as a CompareOp would attempt to retape at invalid
  // underflowed ETA trials instead of letting the optimizer reject them.
  return scalar_finite(value) && scalar_value(value) > 0.0;
}

template <class Scalar>
Scalar parameter_value(const ParametersT<Scalar>& parameters,
                       std::initializer_list<const char*> names,
                       double fallback = std::numeric_limits<double>::quiet_NaN()) {
  for (const char* name : names) {
    auto it = parameters.find(name);
    if (it != parameters.end() && scalar_finite(it->second)) return it->second;
  }
  return Scalar(fallback);
}

template <class Scalar>
Scalar positive_parameter(const ParametersT<Scalar>& parameters,
                          std::initializer_list<const char*> names,
                          double fallback = std::numeric_limits<double>::quiet_NaN()) {
  for (const char* name : names) {
    auto it = parameters.find(name);
    if (it != parameters.end() && scalar_positive(it->second)) return it->second;
  }
  return Scalar(fallback);
}

template <class Scalar>
void require_scalar_positive(const Scalar& value, const std::string& name, int advan) {
  if (!scalar_positive(value)) {
    throw std::domain_error("ADVAN" + std::to_string(advan) +
                            " requires a positive " + name + ".");
  }
}

template <class Scalar>
MatrixT<Scalar> solve_linear(MatrixT<Scalar> matrix, MatrixT<Scalar> rhs,
                             const std::string& context) {
  const Eigen::Index n = matrix.rows();
  if (matrix.cols() != n || rhs.rows() != n) {
    throw std::invalid_argument(context + " linear solve has inconsistent dimensions.");
  }
  // Pivot selection is made from the recording point and then becomes a fixed
  // tape structure. Arithmetic after the row choice remains fully
  // differentiable. Retaping is required when a materially different point
  // changes numerical pivoting, which the tape wrapper detects at construction.
  for (Eigen::Index column = 0; column < n; ++column) {
    Eigen::Index pivot = column;
    Scalar largest_value = libertad::scalar_abs(matrix(column, column));
    double largest = scalar_value(largest_value);
    for (Eigen::Index row = column + 1; row < n; ++row) {
      const Scalar candidate_value = libertad::scalar_abs(matrix(row, column));
      if (path_gt(candidate_value, largest_value)) {
        largest_value = candidate_value;
        largest = scalar_value(candidate_value);
        pivot = row;
      }
    }
    if (!std::isfinite(largest) || largest <= 1e-14) {
      throw std::domain_error(context + " linear system is singular at the recording point.");
    }
    if (pivot != column) {
      matrix.row(column).swap(matrix.row(pivot));
      rhs.row(column).swap(rhs.row(pivot));
    }
    const Scalar diagonal = matrix(column, column);
    matrix.row(column) /= diagonal;
    rhs.row(column) /= diagonal;
    for (Eigen::Index row = 0; row < n; ++row) {
      if (row == column) continue;
      const Scalar factor = matrix(row, column);
      matrix.row(row) -= factor * matrix.row(column);
      rhs.row(row) -= factor * rhs.row(column);
    }
  }
  return rhs;
}

template <class Scalar>
MatrixT<Scalar> matrix_exp_pade(const MatrixT<Scalar>& input) {
  if (input.rows() != input.cols()) {
    throw std::invalid_argument("Matrix exponential requires a square matrix.");
  }
  const Eigen::Index n = input.rows();
  Scalar norm_one_value = Scalar(0.0);
  for (Eigen::Index column = 0; column < n; ++column) {
    Scalar sum = Scalar(0.0);
    for (Eigen::Index row = 0; row < n; ++row) {
      const double value = scalar_value(input(row, column));
      if (!std::isfinite(value)) {
        throw std::domain_error("Matrix exponential input contains non-finite values.");
      }
      sum += libertad::scalar_abs(input(row, column));
    }
    // Which column attains the norm is not a structural execution choice.
    // Keep it as a conditional expression and guard only the later Pade
    // scaling boundary; otherwise harmless changes in the largest column
    // force needless retaping.
    norm_one_value = libertad::choose_gt(
      sum, norm_one_value, sum, norm_one_value);
  }
  const double norm_one = scalar_value(norm_one_value);
  constexpr double theta13 = 5.371920351148152;
  int scaling = 0;
  if (path_gt(norm_one_value, Scalar(theta13))) {
    scaling = std::max(0, static_cast<int>(std::ceil(std::log2(norm_one / theta13))));
  }
  const Scalar divisor = Scalar(std::ldexp(1.0, scaling));
  const MatrixT<Scalar> a = input / divisor;
  const MatrixT<Scalar> identity = MatrixT<Scalar>::Identity(n, n);
  const MatrixT<Scalar> a2 = a * a;
  const MatrixT<Scalar> a4 = a2 * a2;
  const MatrixT<Scalar> a6 = a4 * a2;
  const double b[] = {
    64764752532480000.0, 32382376266240000.0, 7771770303897600.0,
    1187353796428800.0, 129060195264000.0, 10559470521600.0,
    670442572800.0, 33522128640.0, 1323241920.0, 40840800.0,
    960960.0, 16380.0, 182.0, 1.0
  };
  const MatrixT<Scalar> u = a * (
    a6 * (Scalar(b[13]) * a6 + Scalar(b[11]) * a4 + Scalar(b[9]) * a2) +
    Scalar(b[7]) * a6 + Scalar(b[5]) * a4 + Scalar(b[3]) * a2 + Scalar(b[1]) * identity
  );
  const MatrixT<Scalar> v =
    a6 * (Scalar(b[12]) * a6 + Scalar(b[10]) * a4 + Scalar(b[8]) * a2) +
    Scalar(b[6]) * a6 + Scalar(b[4]) * a4 + Scalar(b[2]) * a2 + Scalar(b[0]) * identity;
  MatrixT<Scalar> result = solve_linear(
    MatrixT<Scalar>(v - u), MatrixT<Scalar>(v + u), "Matrix exponential");
  for (int i = 0; i < scaling; ++i) result = result * result;
  return result;
}

template <class Scalar>
Scalar scalar_sinh_t(const Scalar& value) {
  using std::sinh;
  return sinh(value);
}

template <class Scalar>
Scalar scalar_cosh_t(const Scalar& value) {
  using std::cosh;
  return cosh(value);
}

// sinh(z) / z with a smooth series around zero. The series branch is selected
// at tape-recording time, but remains accurate after a 1000-fold movement from
// its threshold, which covers the default NONMEM-style parameter bounds.
template <class Scalar>
Scalar scalar_sinhc_t(const Scalar& value) {
  if (path_gt(libertad::scalar_abs(value), Scalar(1e-4))) {
    return scalar_sinh_t(value) / value;
  }
  const Scalar square = value * value;
  return Scalar(1.0) + square * (
    Scalar(1.0 / 6.0) + square * (
      Scalar(1.0 / 120.0) + square * (
        Scalar(1.0 / 5040.0) + square * Scalar(1.0 / 362880.0)
      )
    )
  );
}

// Stable first divided difference of exp(lambda * time). This is the core of
// the absorption kernels and remains finite when two disposition rates meet.
template <class Scalar>
Scalar divided_exp_t(const Scalar& left, const Scalar& right,
                     const Scalar& time) {
  const Scalar mean = Scalar(0.5) * (left + right);
  const Scalar delta = Scalar(0.5) * (left - right) * time;
  return time * libertad::scalar_exp(mean * time) * scalar_sinhc_t(delta);
}

template <class Scalar>
struct TwoByTwoExponentialT {
  MatrixT<Scalar> transition;
  MatrixT<Scalar> centered;
  Scalar mean = Scalar(0.0);
  Scalar delta = Scalar(0.0);
};

template <class Scalar>
TwoByTwoExponentialT<Scalar> two_by_two_exp_t(
    const MatrixT<Scalar>& matrix, const Scalar& time) {
  if (matrix.rows() != 2 || matrix.cols() != 2) {
    throw std::invalid_argument("Two-by-two propagation requires a 2 x 2 matrix.");
  }
  TwoByTwoExponentialT<Scalar> result;
  result.mean = Scalar(0.5) * (matrix(0, 0) + matrix(1, 1));
  result.centered = matrix - result.mean * MatrixT<Scalar>::Identity(2, 2);
  const Scalar half_difference = Scalar(0.5) * (matrix(0, 0) - matrix(1, 1));
  const Scalar discriminant =
    half_difference * half_difference + matrix(0, 1) * matrix(1, 0);
  if (path_lt(discriminant, Scalar(-1e-12))) {
    throw std::domain_error("The ADVAN two-compartment transition has complex rates.");
  }
  result.delta = libertad::scalar_sqrt(
    path_lt(discriminant, Scalar(0.0)) ? Scalar(0.0) : discriminant);
  const Scalar scaled_delta = result.delta * time;
  const Scalar multiplier = libertad::scalar_exp(result.mean * time);
  result.transition = multiplier * (
    scalar_cosh_t(scaled_delta) * MatrixT<Scalar>::Identity(2, 2) +
    time * scalar_sinhc_t(scaled_delta) * result.centered
  );
  return result;
}

template <class Scalar>
MatrixT<Scalar> advan2_transition_t(const MatrixT<Scalar>& k,
                                    const Scalar& time) {
  MatrixT<Scalar> result = MatrixT<Scalar>::Zero(2, 2);
  result(0, 0) = libertad::scalar_exp(k(0, 0) * time);
  result(1, 1) = libertad::scalar_exp(k(1, 1) * time);
  result(1, 0) = k(1, 0) * divided_exp_t(k(0, 0), k(1, 1), time);
  return result;
}

template <class Scalar>
MatrixT<Scalar> advan4_transition_t(const MatrixT<Scalar>& k,
                                    const Scalar& time) {
  MatrixT<Scalar> result = MatrixT<Scalar>::Zero(3, 3);
  const Scalar depot_rate = k(0, 0);
  result(0, 0) = libertad::scalar_exp(depot_rate * time);
  const MatrixT<Scalar> central = k.block(1, 1, 2, 2);
  const TwoByTwoExponentialT<Scalar> disposition =
    two_by_two_exp_t(central, time);
  result.block(1, 1, 2, 2) = disposition.transition;

  // f(Kc)b, where f(lambda) is the divided exponential between a
  // disposition eigenvalue and the depot eigenvalue. This avoids the
  // singular (Kc + KA I)^-1 formula when KA equals a hybrid rate.
  const Scalar delta_value = libertad::scalar_abs(disposition.delta);
  const Scalar scale = libertad::choose_gt(
    libertad::scalar_abs(disposition.mean), Scalar(1.0),
    libertad::scalar_abs(disposition.mean), Scalar(1.0));
  if (path_le(delta_value, Scalar(1e-8) * scale)) {
    return matrix_exp_pade(MatrixT<Scalar>(k * time));
  }
  const Scalar lambda_plus = disposition.mean + disposition.delta;
  const Scalar lambda_minus = disposition.mean - disposition.delta;
  const Scalar f_plus = divided_exp_t(lambda_plus, depot_rate, time);
  const Scalar f_minus = divided_exp_t(lambda_minus, depot_rate, time);
  const VectorT<Scalar> coupling = k.block(1, 0, 2, 1);
  const VectorT<Scalar> cross =
    Scalar(0.5) * (f_plus + f_minus) * coupling +
    Scalar(0.5) * (f_plus - f_minus) / disposition.delta *
      (disposition.centered * coupling);
  result.block(1, 0, 2, 1) = cross;
  return result;
}

inline bool specialized_advan_number(int advan) {
  return advan == 1 || advan == 2 || advan == 3 || advan == 4 ||
    advan == 11 || advan == 12;
}

inline bool use_specialized_advan(const ModelEngine& engine) {
  return engine.specialized_advan && !engine.is_ode() &&
    !engine.matrix_graph.enabled && specialized_advan_number(engine.advan);
}

inline std::string propagation_kernel_name(const ModelEngine& engine) {
  if (engine.dde_enabled) return "dde-rk4-method-of-steps";
  if (engine.dae_enabled) {
    return engine.advan == 13 ? "dae-advan13-implicit-newton" :
      "dae-advan6-rk45-newton";
  }
  if (engine.is_ode()) {
    return engine.advan == 13 ? "advan13-implicit" : "advan6-rk45";
  }
  if (use_specialized_advan(engine)) {
    return "specialized-advan" + std::to_string(engine.advan);
  }
  return "general-matrix-exponential";
}

template <class Scalar>
MatrixT<Scalar> specialized_advan_transition_t(
    int advan, const MatrixT<Scalar>& k, const Scalar& time) {
  if (path_le(time, Scalar(0.0))) {
    return MatrixT<Scalar>::Identity(k.rows(), k.cols());
  }
  if (advan == 1) {
    MatrixT<Scalar> result(1, 1);
    result(0, 0) = libertad::scalar_exp(k(0, 0) * time);
    return result;
  }
  if (advan == 2) return advan2_transition_t(k, time);
  if (advan == 3) return two_by_two_exp_t(k, time).transition;
  if (advan == 4) return advan4_transition_t(k, time);
  // The three-compartment kernels retain the robust Padé transition on the
  // native 3 x 3 or 4 x 4 system. Their specialization comes from avoiding
  // the larger affine augmentation and using the exact phi-one offset below.
  if (advan == 11 || advan == 12) {
    return matrix_exp_pade(MatrixT<Scalar>(k * time));
  }
  throw std::invalid_argument("No specialized propagation kernel exists for this ADVAN.");
}

template <class Scalar>
Scalar matrix_one_norm_value(const MatrixT<Scalar>& matrix) {
  Scalar norm_value = Scalar(0.0);
  for (Eigen::Index column = 0; column < matrix.cols(); ++column) {
    Scalar sum = Scalar(0.0);
    for (Eigen::Index row = 0; row < matrix.rows(); ++row) {
      sum += libertad::scalar_abs(matrix(row, column));
    }
    norm_value = libertad::choose_gt(sum, norm_value, sum, norm_value);
  }
  return norm_value;
}

template <class Scalar>
bool zero_input_t(const VectorT<Scalar>& input) {
  for (Eigen::Index row = 0; row < input.rows(); ++row) {
    if (path_ne(input[row], Scalar(0.0))) return false;
  }
  return true;
}

template <class Scalar>
AffineMapT<Scalar> specialized_advan_affine_map_t(
    int advan, const MatrixT<Scalar>& k, const VectorT<Scalar>& input,
    const Scalar& time) {
  const Eigen::Index n = k.rows();
  if (k.cols() != n || input.size() != n) {
    throw std::invalid_argument("Specialized ADVAN propagation dimensions are inconsistent.");
  }
  if (path_lt(time, Scalar(-1e-12)) || !scalar_finite(time)) {
    throw std::domain_error("Propagation interval must be finite and non-negative.");
  }
  if (path_le(time, Scalar(0.0))) {
    return {MatrixT<Scalar>::Identity(n, n), VectorT<Scalar>::Zero(n)};
  }
  const MatrixT<Scalar> transition = specialized_advan_transition_t(advan, k, time);
  if (zero_input_t(input)) {
    return {transition, VectorT<Scalar>::Zero(n)};
  }

  VectorT<Scalar> offset(n);
  const MatrixT<Scalar> scaled = k * time;
  if (path_le(matrix_one_norm_value(scaled), Scalar(1e-4))) {
    // phi_1(A) = I + A/2! + A^2/3! + ... avoids cancellation in
    // K^-1(exp(Kt)-I) for very short intervals or very slow rates.
    const MatrixT<Scalar> identity = MatrixT<Scalar>::Identity(n, n);
    MatrixT<Scalar> phi = identity;
    MatrixT<Scalar> term = identity;
    for (int order = 1; order <= 18; ++order) {
      term = term * scaled / Scalar(static_cast<double>(order));
      phi += term / Scalar(static_cast<double>(order + 1));
    }
    offset = time * phi * input;
  } else {
    MatrixT<Scalar> rhs(n, 1);
    rhs.col(0) = (transition - MatrixT<Scalar>::Identity(n, n)) * input;
    offset = solve_linear(k, rhs, "Specialized ADVAN affine offset").col(0);
  }
  return {transition, offset};
}

template <class Scalar>
AffineMapT<Scalar> affine_map_t(const MatrixT<Scalar>& k,
                                const VectorT<Scalar>& input, const Scalar& dt) {
  const Eigen::Index n = k.rows();
  if (k.cols() != n || input.size() != n) {
    throw std::invalid_argument("Affine propagation dimensions are inconsistent.");
  }
  const double dt_value = scalar_value(dt);
  if (path_lt(dt, Scalar(-1e-12)) || !std::isfinite(dt_value)) {
    throw std::domain_error("Propagation interval must be finite and non-negative.");
  }
  if (path_le(dt, Scalar(0.0))) {
    return {MatrixT<Scalar>::Identity(n, n), VectorT<Scalar>::Zero(n)};
  }
  MatrixT<Scalar> augmented = MatrixT<Scalar>::Zero(n + 1, n + 1);
  augmented.topLeftCorner(n, n) = k;
  augmented.topRightCorner(n, 1) = input;
  MatrixT<Scalar> exponential = matrix_exp_pade(MatrixT<Scalar>(augmented * dt));
  return {exponential.topLeftCorner(n, n), exponential.topRightCorner(n, 1)};
}

template <class Scalar>
VectorT<Scalar> propagate_t(const MatrixT<Scalar>& k,
                            const VectorT<Scalar>& input, const Scalar& dt,
                            const VectorT<Scalar>& state) {
  AffineMapT<Scalar> map = affine_map_t(k, input, dt);
  return map.transition * state + map.offset;
}

template <class Scalar>
AffineMapT<Scalar> engine_affine_map_t(
    const ModelEngine& engine, const MatrixT<Scalar>& k,
    const VectorT<Scalar>& input, const Scalar& dt) {
  if (use_specialized_advan(engine)) {
    return specialized_advan_affine_map_t(engine.advan, k, input, dt);
  }
  return affine_map_t(k, input, dt);
}

template <class Scalar>
MatrixT<Scalar> engine_transition_t(
    const ModelEngine& engine, const MatrixT<Scalar>& k, const Scalar& dt) {
  if (use_specialized_advan(engine)) {
    return specialized_advan_transition_t(engine.advan, k, dt);
  }
  return matrix_exp_pade(MatrixT<Scalar>(k * dt));
}

template <class Scalar>
VectorT<Scalar> solve_periodic_t(const MatrixT<Scalar>& transition,
                                 const VectorT<Scalar>& offset,
                                 const std::string& context) {
  MatrixT<Scalar> system = MatrixT<Scalar>::Identity(
    transition.rows(), transition.cols()) - transition;
  MatrixT<Scalar> rhs(offset.rows(), 1);
  rhs.col(0) = offset;
  return solve_linear(system, rhs, context + " steady-state").col(0);
}

template <class Scalar>
TopologyT<Scalar> build_topology_t(int advan, const ParametersT<Scalar>& p) {
  const int n = (advan == 1 ? 1 :
                 advan == 2 || advan == 3 ? 2 :
                 advan == 4 || advan == 11 ? 3 :
                 advan == 12 ? 4 : 0);
  if (n == 0) {
    throw std::invalid_argument("Differentiable analytical engine supports ADVAN1-4/11/12.");
  }
  TopologyT<Scalar> topology;
  topology.k = MatrixT<Scalar>::Zero(n, n);
  topology.default_scales.assign(static_cast<std::size_t>(n), Scalar(1.0));
  const Scalar vc = positive_parameter(p, {"VC", "V1", "V"});
  const Scalar vp1 = positive_parameter(p, {"VP", "VP1", "V2"});
  const Scalar vp2 = positive_parameter(p, {"VP2", "V3"});
  const Scalar cl = positive_parameter(p, {"CL"});
  const Scalar q1 = positive_parameter(p, {"Q2", "Q", "Q1"});
  const Scalar q2 = positive_parameter(p, {"Q3", "Q4"});

  if (advan == 1) {
    Scalar k10 = positive_parameter(p, {"K10", "K"});
    if (!scalar_positive(k10) && scalar_positive(cl) && scalar_positive(vc)) k10 = cl / vc;
    require_scalar_positive(k10, "K10 or CL/V", advan);
    topology.k(0, 0) = -k10;
    topology.default_scales[0] = scalar_positive(vc) ? vc : Scalar(1.0);
  } else if (advan == 2) {
    const Scalar ka = positive_parameter(p, {"KA"});
    Scalar k20 = positive_parameter(p, {"K20", "K10", "K"});
    if (!scalar_positive(k20) && scalar_positive(cl) && scalar_positive(vc)) k20 = cl / vc;
    require_scalar_positive(ka, "KA", advan);
    require_scalar_positive(k20, "K20 or CL/V", advan);
    topology.k(0, 0) = -ka;
    topology.k(1, 0) = ka;
    topology.k(1, 1) = -k20;
    topology.default_observation = 1;
    topology.default_scales[1] = scalar_positive(vc) ? vc : Scalar(1.0);
  } else if (advan == 3) {
    Scalar k10 = positive_parameter(p, {"K10"});
    Scalar k12 = positive_parameter(p, {"K12"});
    Scalar k21 = positive_parameter(p, {"K21"});
    if (!scalar_positive(k10) && scalar_positive(cl) && scalar_positive(vc)) k10 = cl / vc;
    if (!scalar_positive(k12) && scalar_positive(q1) && scalar_positive(vc)) k12 = q1 / vc;
    if (!scalar_positive(k21) && scalar_positive(q1) && scalar_positive(vp1)) k21 = q1 / vp1;
    require_scalar_positive(k10, "K10 or CL/V1", advan);
    require_scalar_positive(k12, "K12 or Q/V1", advan);
    require_scalar_positive(k21, "K21 or Q/V2", advan);
    topology.k(0, 0) = -(k10 + k12);
    topology.k(0, 1) = k21;
    topology.k(1, 0) = k12;
    topology.k(1, 1) = -k21;
    topology.default_scales[0] = scalar_positive(vc) ? vc : Scalar(1.0);
    topology.default_scales[1] = scalar_positive(vp1) ? vp1 : Scalar(1.0);
  } else if (advan == 4) {
    const Scalar ka = positive_parameter(p, {"KA"});
    Scalar k20 = positive_parameter(p, {"K20", "K10"});
    Scalar k23 = positive_parameter(p, {"K23", "K12"});
    Scalar k32 = positive_parameter(p, {"K32", "K21"});
    if (!scalar_positive(k20) && scalar_positive(cl) && scalar_positive(vc)) k20 = cl / vc;
    if (!scalar_positive(k23) && scalar_positive(q1) && scalar_positive(vc)) k23 = q1 / vc;
    if (!scalar_positive(k32) && scalar_positive(q1) && scalar_positive(vp1)) k32 = q1 / vp1;
    require_scalar_positive(ka, "KA", advan);
    require_scalar_positive(k20, "K20 or CL/VC", advan);
    require_scalar_positive(k23, "K23 or Q/VC", advan);
    require_scalar_positive(k32, "K32 or Q/VP", advan);
    topology.k(0, 0) = -ka;
    topology.k(1, 0) = ka;
    topology.k(1, 1) = -(k20 + k23);
    topology.k(1, 2) = k32;
    topology.k(2, 1) = k23;
    topology.k(2, 2) = -k32;
    topology.default_observation = 1;
    topology.default_scales[1] = scalar_positive(vc) ? vc : Scalar(1.0);
    topology.default_scales[2] = scalar_positive(vp1) ? vp1 : Scalar(1.0);
  } else if (advan == 11) {
    Scalar k10 = positive_parameter(p, {"K10"});
    Scalar k12 = positive_parameter(p, {"K12"});
    Scalar k21 = positive_parameter(p, {"K21"});
    Scalar k13 = positive_parameter(p, {"K13"});
    Scalar k31 = positive_parameter(p, {"K31"});
    if (!scalar_positive(k10) && scalar_positive(cl) && scalar_positive(vc)) k10 = cl / vc;
    if (!scalar_positive(k12) && scalar_positive(q1) && scalar_positive(vc)) k12 = q1 / vc;
    if (!scalar_positive(k21) && scalar_positive(q1) && scalar_positive(vp1)) k21 = q1 / vp1;
    if (!scalar_positive(k13) && scalar_positive(q2) && scalar_positive(vc)) k13 = q2 / vc;
    if (!scalar_positive(k31) && scalar_positive(q2) && scalar_positive(vp2)) k31 = q2 / vp2;
    require_scalar_positive(k10, "K10 or CL/V1", advan);
    require_scalar_positive(k12, "K12 or Q2/V1", advan);
    require_scalar_positive(k21, "K21 or Q2/V2", advan);
    require_scalar_positive(k13, "K13 or Q3/V1", advan);
    require_scalar_positive(k31, "K31 or Q3/V3", advan);
    topology.k(0, 0) = -(k10 + k12 + k13);
    topology.k(0, 1) = k21;
    topology.k(0, 2) = k31;
    topology.k(1, 0) = k12;
    topology.k(1, 1) = -k21;
    topology.k(2, 0) = k13;
    topology.k(2, 2) = -k31;
    topology.default_scales[0] = scalar_positive(vc) ? vc : Scalar(1.0);
    topology.default_scales[1] = scalar_positive(vp1) ? vp1 : Scalar(1.0);
    topology.default_scales[2] = scalar_positive(vp2) ? vp2 : Scalar(1.0);
  } else {
    const Scalar ka = positive_parameter(p, {"KA"});
    Scalar k20 = positive_parameter(p, {"K20", "K10"});
    Scalar k23 = positive_parameter(p, {"K23", "K12"});
    Scalar k32 = positive_parameter(p, {"K32", "K21"});
    Scalar k24 = positive_parameter(p, {"K24", "K13"});
    Scalar k42 = positive_parameter(p, {"K42", "K31"});
    if (!scalar_positive(k20) && scalar_positive(cl) && scalar_positive(vc)) k20 = cl / vc;
    if (!scalar_positive(k23) && scalar_positive(q1) && scalar_positive(vc)) k23 = q1 / vc;
    if (!scalar_positive(k32) && scalar_positive(q1) && scalar_positive(vp1)) k32 = q1 / vp1;
    if (!scalar_positive(k24) && scalar_positive(q2) && scalar_positive(vc)) k24 = q2 / vc;
    if (!scalar_positive(k42) && scalar_positive(q2) && scalar_positive(vp2)) k42 = q2 / vp2;
    require_scalar_positive(ka, "KA", advan);
    require_scalar_positive(k20, "K20 or CL/VC", advan);
    require_scalar_positive(k23, "K23 or Q2/VC", advan);
    require_scalar_positive(k32, "K32 or Q2/VP1", advan);
    require_scalar_positive(k24, "K24 or Q3/VC", advan);
    require_scalar_positive(k42, "K42 or Q3/VP2", advan);
    topology.k(0, 0) = -ka;
    topology.k(1, 0) = ka;
    topology.k(1, 1) = -(k20 + k23 + k24);
    topology.k(1, 2) = k32;
    topology.k(1, 3) = k42;
    topology.k(2, 1) = k23;
    topology.k(2, 2) = -k32;
    topology.k(3, 1) = k24;
    topology.k(3, 3) = -k42;
    topology.default_observation = 1;
    topology.default_scales[1] = scalar_positive(vc) ? vc : Scalar(1.0);
    topology.default_scales[2] = scalar_positive(vp1) ? vp1 : Scalar(1.0);
    topology.default_scales[3] = scalar_positive(vp2) ? vp2 : Scalar(1.0);
  }
  return topology;
}

template <class Scalar>
TopologyT<Scalar> build_graph_topology_t(const MatrixGraph& graph,
                                         const ParametersT<Scalar>& p) {
  if (!graph.enabled || graph.names.empty()) {
    throw std::invalid_argument("Matrix graph is empty.");
  }
  const int n = static_cast<int>(graph.names.size());
  TopologyT<Scalar> topology;
  topology.k = MatrixT<Scalar>::Zero(n, n);
  topology.default_scales.assign(static_cast<std::size_t>(n), Scalar(1.0));
  for (int i = 0; i < n; ++i) {
    if (i < static_cast<int>(graph.scale_parameters.size()) &&
        !graph.scale_parameters[static_cast<std::size_t>(i)].empty()) {
      const std::string& name = graph.scale_parameters[static_cast<std::size_t>(i)];
      auto it = p.find(name);
      if (it == p.end() || !scalar_positive(it->second)) {
        throw std::domain_error("Matrix graph scale parameter '" + name + "' must be positive.");
      }
      topology.default_scales[static_cast<std::size_t>(i)] = it->second;
    }
  }
  for (const MatrixFlow& flow : graph.flows) {
    auto parameter = p.find(flow.parameter);
    if (parameter == p.end() || !scalar_positive(parameter->second)) {
      throw std::domain_error("Matrix graph flow parameter '" + flow.parameter + "' must be positive.");
    }
    Scalar rate = parameter->second;
    if (flow.type == "clearance") {
      auto volume = p.find(flow.volume_parameter);
      if (volume == p.end() || !scalar_positive(volume->second)) {
        throw std::domain_error("Matrix graph volume parameter '" +
                                flow.volume_parameter + "' must be positive.");
      }
      rate /= volume->second;
    }
    topology.k(flow.from, flow.from) -= rate;
    if (flow.to >= 0) topology.k(flow.to, flow.from) += rate;
  }
  return topology;
}

// Numeric covariates used by PRED/DES can change between subjects without
// changing the event topology.  During prediction-tape recording these values
// are CppAD dynamic parameters, allowing one structural tape to be reused.
// Event ordering, dosing-mode, compartment, mixture, and IOV fields remain
// ordinary constants because changes to them require a different tape.
template <class Scalar>
struct DynamicDataT {
  int n_rows = 0;
  std::unordered_map<std::string, std::size_t> column_positions;
  std::vector<Scalar> values;

  const Scalar* find(const std::string& name, int row) const {
    auto it = column_positions.find(name);
    if (it == column_positions.end()) return nullptr;
    const std::size_t position = it->second * static_cast<std::size_t>(n_rows) +
      static_cast<std::size_t>(row);
    return &values.at(position);
  }
};

template <class Scalar>
Scalar dynamic_row_value(const Rcpp::DataFrame& data, const std::string& name,
                         int row, const DynamicDataT<Scalar>* dynamic_data) {
  if (dynamic_data != nullptr) {
    const Scalar* value = dynamic_data->find(name, row);
    if (value != nullptr) return *value;
  }
  return Scalar(data_value(data, name, row));
}

template <class Scalar>
Scalar dynamic_row_optional(const Rcpp::DataFrame& data, const std::string& name,
                            int row, double fallback,
                            const DynamicDataT<Scalar>* dynamic_data) {
  if (dynamic_data != nullptr) {
    const Scalar* value = dynamic_data->find(name, row);
    if (value != nullptr) return *value;
  }
  return Scalar(row_optional(data, name, row, fallback));
}

template <class Scalar>
ParametersT<Scalar> evaluate_parameters_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data, int row, int subject,
    const std::vector<Scalar>& theta, const std::vector<Scalar>& eta,
    int eta_columns, const std::vector<Scalar>& sigma, int mixture_number,
    const DynamicDataT<Scalar>* dynamic_data = nullptr) {
  std::vector<Scalar> inputs(engine.pred->input_names.size(), Scalar(0.0));
  for (std::size_t i = 0; i < engine.pred->input_names.size(); ++i) {
    const std::string& name = engine.pred->input_names[i];
    int index = indexed_name(name, "THETA_");
    if (index >= 0) {
      if (index >= static_cast<int>(theta.size())) throw std::out_of_range("THETA index exceeds values.");
      inputs[i] = theta[static_cast<std::size_t>(index)];
      continue;
    }
    index = indexed_name(name, "ETA_");
    if (index >= 0) {
      const int column = eta_column(engine, data, row, index, eta_columns);
      const std::size_t position = static_cast<std::size_t>(subject * eta_columns + column);
      if (index >= engine.n_eta || position >= eta.size()) throw std::out_of_range("ETA index exceeds values.");
      inputs[i] = eta[position];
      continue;
    }
    index = indexed_name(name, "SIGMA_");
    if (index >= 0) {
      if (index >= static_cast<int>(sigma.size())) throw std::out_of_range("SIGMA index exceeds values.");
      inputs[i] = sigma[static_cast<std::size_t>(index)];
      continue;
    }
    if (starts_with(name, "ERR_") || name == "F") continue;
    if (name == "MIXNUM") {
      inputs[i] = Scalar(mixture_number);
      continue;
    }
    const Scalar value = dynamic_row_value(data, name, row, dynamic_data);
    if (!std::isfinite(scalar_value(value))) {
      throw std::domain_error("PRED input '" + name + "' is non-finite at row " +
                              std::to_string(row + 1) + ".");
    }
    inputs[i] = value;
  }
  std::vector<Scalar> output = engine.pred->eval_outputs(inputs, engine.all_outputs);
  ParametersT<Scalar> parameters;
  for (std::size_t i = 0; i < output.size(); ++i) {
    parameters[engine.pred->output_names[i]] = output[i];
  }
  return parameters;
}

template <class Scalar>
VectorT<Scalar> evaluate_algebraic_residuals_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data,
    int row, int subject, const Scalar& time, const VectorT<Scalar>& state,
    const ParametersT<Scalar>& parameters,
    const std::vector<Scalar>& theta, const std::vector<Scalar>& eta,
    int eta_columns, const std::vector<Scalar>& sigma, int mixture_number,
    const VectorT<Scalar>& algebraic,
    const DynamicDataT<Scalar>* dynamic_data = nullptr) {
  if (!engine.alg) throw std::logic_error("DAE algebraic residual program is missing.");
  std::vector<Scalar> inputs(engine.alg->input_names.size(), Scalar(0.0));
  for (std::size_t i = 0; i < engine.alg->input_names.size(); ++i) {
    const std::string& name = engine.alg->input_names[i];
    int index = indexed_name(name, "A_");
    if (index >= 0) { inputs[i] = state[index]; continue; }
    if (name == "T") { inputs[i] = time; continue; }
    auto variable = std::find(engine.dae_variables.begin(), engine.dae_variables.end(), name);
    if (variable != engine.dae_variables.end()) {
      inputs[i] = algebraic[std::distance(engine.dae_variables.begin(), variable)];
      continue;
    }
    auto parameter = parameters.find(name);
    if (parameter != parameters.end()) { inputs[i] = parameter->second; continue; }
    index = indexed_name(name, "THETA_");
    if (index >= 0) { inputs[i] = theta.at(static_cast<std::size_t>(index)); continue; }
    index = indexed_name(name, "ETA_");
    if (index >= 0) {
      const int column = eta_column(engine, data, row, index, eta_columns);
      inputs[i] = eta.at(static_cast<std::size_t>(subject * eta_columns + column));
      continue;
    }
    index = indexed_name(name, "SIGMA_");
    if (index >= 0) { inputs[i] = sigma.at(static_cast<std::size_t>(index)); continue; }
    if (starts_with(name, "ERR_") || name == "F") continue;
    if (name == "MIXNUM") { inputs[i] = Scalar(mixture_number); continue; }
    inputs[i] = dynamic_row_value(data, name, row, dynamic_data);
  }
  const std::vector<Scalar> output = engine.alg->eval_outputs(inputs, engine.algebraic_outputs);
  VectorT<Scalar> residual(static_cast<Eigen::Index>(output.size()));
  for (std::size_t i = 0; i < output.size(); ++i) residual[static_cast<Eigen::Index>(i)] = output[i];
  return residual;
}

template <class Scalar>
VectorT<Scalar> solve_algebraic_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data,
    int row, int subject, const Scalar& time, const VectorT<Scalar>& state,
    const ParametersT<Scalar>& parameters,
    const std::vector<Scalar>& theta, const std::vector<Scalar>& eta,
    int eta_columns, const std::vector<Scalar>& sigma, int mixture_number,
    const DynamicDataT<Scalar>* dynamic_data = nullptr) {
  VectorT<Scalar> value(static_cast<Eigen::Index>(engine.dae_initial.size()));
  for (std::size_t i = 0; i < engine.dae_initial.size(); ++i) value[static_cast<Eigen::Index>(i)] = Scalar(engine.dae_initial[i]);
  for (int iteration = 0; iteration < engine.dae_maxit; ++iteration) {
    const VectorT<Scalar> residual = evaluate_algebraic_residuals_t(
      engine, data, row, subject, time, state, parameters, theta, eta,
      eta_columns, sigma, mixture_number, value, dynamic_data);
    Scalar maximum = Scalar(0.0);
    for (Eigen::Index i = 0; i < residual.size(); ++i) {
      maximum = libertad::choose_gt(
        libertad::scalar_abs(residual[i]), maximum,
        libertad::scalar_abs(residual[i]), maximum);
    }
    if (path_le(maximum, Scalar(engine.dae_tolerance))) return value;
    MatrixT<Scalar> jacobian = MatrixT<Scalar>::Zero(value.size(), value.size());
    for (Eigen::Index column = 0; column < value.size(); ++column) {
      const double delta = engine.dae_jacobian_step *
        std::max(1.0, std::abs(scalar_value(value[column])));
      VectorT<Scalar> plus = value; plus[column] += Scalar(delta);
      VectorT<Scalar> minus = value; minus[column] -= Scalar(delta);
      const VectorT<Scalar> upper = evaluate_algebraic_residuals_t(
        engine, data, row, subject, time, state, parameters, theta, eta,
        eta_columns, sigma, mixture_number, plus, dynamic_data);
      const VectorT<Scalar> lower = evaluate_algebraic_residuals_t(
        engine, data, row, subject, time, state, parameters, theta, eta,
        eta_columns, sigma, mixture_number, minus, dynamic_data);
      jacobian.col(column) = (upper - lower) / Scalar(2.0 * delta);
      if (!engine.dae_sparsity.empty()) {
        for (Eigen::Index row_index = 0; row_index < value.size(); ++row_index) {
          if (!engine.dae_sparsity[static_cast<std::size_t>(row_index * value.size() + column)]) {
            jacobian(row_index, column) = Scalar(0.0);
          }
        }
      }
    }
    VectorT<Scalar> update = VectorT<Scalar>::Zero(value.size());
    for (std::size_t block = 0; block < engine.dae_block_rows.size(); ++block) {
      const auto& rows = engine.dae_block_rows[block];
      const auto& columns = engine.dae_block_columns[block];
      MatrixT<Scalar> local(rows.size(), columns.size());
      MatrixT<Scalar> rhs(rows.size(), 1);
      for (std::size_t local_row = 0; local_row < rows.size(); ++local_row) {
        rhs(static_cast<Eigen::Index>(local_row), 0) = -residual[rows[local_row]];
        for (std::size_t local_column = 0; local_column < columns.size(); ++local_column) {
          local(static_cast<Eigen::Index>(local_row),
                static_cast<Eigen::Index>(local_column)) =
            jacobian(rows[local_row], columns[local_column]);
        }
      }
      const MatrixT<Scalar> local_update = solve_linear(local, rhs, "DAE Newton block");
      for (std::size_t local_column = 0; local_column < columns.size(); ++local_column) {
        update[columns[local_column]] = local_update(static_cast<Eigen::Index>(local_column), 0);
      }
    }
    value += update;
  }
  const VectorT<Scalar> residual = evaluate_algebraic_residuals_t(
    engine, data, row, subject, time, state, parameters, theta, eta,
    eta_columns, sigma, mixture_number, value, dynamic_data);
  double maximum = 0.0;
  for (Eigen::Index i = 0; i < residual.size(); ++i) {
    maximum = std::max(maximum, std::abs(scalar_value(residual[i])));
  }
  if (maximum > 10.0 * engine.dae_tolerance) {
    throw std::runtime_error("DAE algebraic Newton solve did not converge.");
  }
  return value;
}

template <class Scalar>
struct DdeHistoryT {
  std::vector<Scalar> time;
  std::vector<VectorT<Scalar>> state;
  std::vector<double> baseline;

  void reset(const Scalar& at, const VectorT<Scalar>& value,
             const std::vector<double>& history) {
    time.assign(1U, at); state.assign(1U, value); baseline = history;
  }
  void append(const Scalar& at, const VectorT<Scalar>& value) {
    if (!time.empty() && std::abs(scalar_value(time.back() - at)) <= 1e-12) {
      state.back() = value; return;
    }
    time.push_back(at); state.push_back(value);
  }
  Scalar at(const Scalar& target, int component) const {
    if (component < 0 || component >= static_cast<int>(baseline.size())) {
      throw std::out_of_range("DDE lag state is outside the state vector.");
    }
    const double target_value = scalar_value(target);
    if (time.empty() || target_value < scalar_value(time.front()) - 1e-12) {
      return Scalar(baseline[static_cast<std::size_t>(component)]);
    }
    if (target_value >= scalar_value(time.back()) - 1e-12) return state.back()[component];
    std::size_t right = 1U;
    while (right < time.size() && scalar_value(time[right]) <= target_value) ++right;
    const std::size_t left = right - 1U;
    const Scalar fraction = (target - time[left]) / (time[right] - time[left]);
    return state[left][component] + fraction * (state[right][component] - state[left][component]);
  }
};

template <class Scalar>
VectorT<Scalar> evaluate_derivatives_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data,
    int row, int subject, const Scalar& time, const VectorT<Scalar>& state,
    const ParametersT<Scalar>& parameters,
    const std::vector<Scalar>& theta, const std::vector<Scalar>& eta,
    int eta_columns, const std::vector<Scalar>& sigma, int mixture_number,
    const DynamicDataT<Scalar>* dynamic_data = nullptr,
    const VectorT<Scalar>* lag_values = nullptr) {
  if (!engine.des) throw std::logic_error("ODE derivative program is missing.");
  const VectorT<Scalar> algebraic = engine.dae_enabled ? solve_algebraic_t(
    engine, data, row, subject, time, state, parameters, theta, eta,
    eta_columns, sigma, mixture_number, dynamic_data) : VectorT<Scalar>();
  std::vector<Scalar> inputs(engine.des->input_names.size(), Scalar(0.0));
  for (std::size_t i = 0; i < engine.des->input_names.size(); ++i) {
    const std::string& name = engine.des->input_names[i];
    int index = indexed_name(name, "A_");
    if (index >= 0) {
      if (index >= state.size()) throw std::out_of_range("A() index exceeds ODE state dimension.");
      inputs[i] = state[index];
      continue;
    }
    if (name == "T") {
      inputs[i] = time;
      continue;
    }
    auto lag = std::find(engine.dde_lag_inputs.begin(), engine.dde_lag_inputs.end(), name);
    if (lag != engine.dde_lag_inputs.end()) {
      if (lag_values == nullptr) throw std::logic_error("DDE lag history is unavailable.");
      inputs[i] = (*lag_values)[std::distance(engine.dde_lag_inputs.begin(), lag)];
      continue;
    }
    auto variable = std::find(engine.dae_variables.begin(), engine.dae_variables.end(), name);
    if (variable != engine.dae_variables.end()) {
      inputs[i] = algebraic[std::distance(engine.dae_variables.begin(), variable)];
      continue;
    }
    auto parameter = parameters.find(name);
    if (parameter != parameters.end()) {
      inputs[i] = parameter->second;
      continue;
    }
    index = indexed_name(name, "THETA_");
    if (index >= 0) {
      inputs[i] = theta.at(static_cast<std::size_t>(index));
      continue;
    }
    index = indexed_name(name, "ETA_");
    if (index >= 0) {
      const int column = eta_column(engine, data, row, index, eta_columns);
      inputs[i] = eta.at(static_cast<std::size_t>(subject * eta_columns + column));
      continue;
    }
    index = indexed_name(name, "SIGMA_");
    if (index >= 0) {
      inputs[i] = sigma.at(static_cast<std::size_t>(index));
      continue;
    }
    if (starts_with(name, "ERR_") || name == "F") continue;
    if (name == "MIXNUM") {
      inputs[i] = Scalar(mixture_number);
      continue;
    }
    const Scalar value = dynamic_row_value(data, name, row, dynamic_data);
    if (!std::isfinite(scalar_value(value))) {
      throw std::domain_error("DES input '" + name + "' is non-finite at row " +
                              std::to_string(row + 1) + ".");
    }
    inputs[i] = value;
  }
  std::vector<Scalar> values = engine.des->eval_outputs(inputs, engine.derivative_outputs);
  VectorT<Scalar> derivative(static_cast<Eigen::Index>(values.size()));
  for (std::size_t i = 0; i < values.size(); ++i) derivative[static_cast<Eigen::Index>(i)] = values[i];
  return derivative;
}

template <class Scalar>
Scalar scaled_error_t(const VectorT<Scalar>& error,
                      const VectorT<Scalar>& before,
                      const VectorT<Scalar>& after,
                      const OdeControl& control) {
  Scalar maximum = Scalar(0.0);
  for (Eigen::Index i = 0; i < error.size(); ++i) {
    const double scale = control.atol + control.rtol * std::max(
      std::abs(scalar_value(before[i])), std::abs(scalar_value(after[i])));
    const Scalar scaled = libertad::scalar_abs(error[i]) / Scalar(scale);
    maximum = libertad::choose_gt(scaled, maximum, scaled, maximum);
  }
  return maximum;
}

template <class Scalar, class Rhs>
VectorT<Scalar> integrate_dopri54_t(const Rhs& rhs, VectorT<Scalar> state,
                                    double from, double to,
                                    const OdeControl& control) {
  if (to <= from) return state;
  double time = from;
  const double span = to - from;
  double h = control.initial_step > 0.0 ? std::min(control.initial_step, span) : span / 10.0;
  h = std::max(h, std::min(span, 1e-8));
  int attempts = 0;
  while (time < to) {
    if (++attempts > control.max_steps) throw std::runtime_error("ADVAN6 exceeded ODE_CONTROL$max_steps.");
    h = std::min(h, to - time);
    const double minimum = 32.0 * std::numeric_limits<double>::epsilon() *
      std::max({1.0, std::abs(time), std::abs(to)});
    if (h < minimum) throw std::runtime_error("ADVAN6 ODE step size underflow.");
    const VectorT<Scalar> k1 = rhs(time, state);
    const VectorT<Scalar> k2 = rhs(time + h * (1.0 / 5.0),
      state + Scalar(h) * (Scalar(1.0 / 5.0) * k1));
    const VectorT<Scalar> k3 = rhs(time + h * (3.0 / 10.0),
      state + Scalar(h) * (Scalar(3.0 / 40.0) * k1 + Scalar(9.0 / 40.0) * k2));
    const VectorT<Scalar> k4 = rhs(time + h * (4.0 / 5.0),
      state + Scalar(h) * (Scalar(44.0 / 45.0) * k1 - Scalar(56.0 / 15.0) * k2 +
                           Scalar(32.0 / 9.0) * k3));
    const VectorT<Scalar> k5 = rhs(time + h * (8.0 / 9.0),
      state + Scalar(h) * (Scalar(19372.0 / 6561.0) * k1 -
        Scalar(25360.0 / 2187.0) * k2 + Scalar(64448.0 / 6561.0) * k3 -
        Scalar(212.0 / 729.0) * k4));
    const VectorT<Scalar> k6 = rhs(time + h,
      state + Scalar(h) * (Scalar(9017.0 / 3168.0) * k1 - Scalar(355.0 / 33.0) * k2 +
        Scalar(46732.0 / 5247.0) * k3 + Scalar(49.0 / 176.0) * k4 -
        Scalar(5103.0 / 18656.0) * k5));
    const VectorT<Scalar> fifth = state + Scalar(h) * (
      Scalar(35.0 / 384.0) * k1 + Scalar(500.0 / 1113.0) * k3 +
      Scalar(125.0 / 192.0) * k4 - Scalar(2187.0 / 6784.0) * k5 +
      Scalar(11.0 / 84.0) * k6);
    const VectorT<Scalar> k7 = rhs(time + h, fifth);
    const VectorT<Scalar> fourth = state + Scalar(h) * (
      Scalar(5179.0 / 57600.0) * k1 + Scalar(7571.0 / 16695.0) * k3 +
      Scalar(393.0 / 640.0) * k4 - Scalar(92097.0 / 339200.0) * k5 +
      Scalar(187.0 / 2100.0) * k6 + Scalar(1.0 / 40.0) * k7);
    const Scalar error_expression = scaled_error_t(
      VectorT<Scalar>(fifth - fourth), state, fifth, control);
    const double error = scalar_value(error_expression);
    if (!std::isfinite(error)) throw std::domain_error("ADVAN6 ODE error estimate is non-finite.");
    if (path_le(error_expression, Scalar(1.0))) {
      state = fifth;
      time += h;
    }
    const double factor = error == 0.0 ? 5.0 :
      std::clamp(0.9 * std::pow(error, -0.2), 0.1, 5.0);
    h *= factor;
  }
  return state;
}

template <class Scalar, class Rhs>
bool implicit_trapezoid_step_t(const Rhs& rhs, const VectorT<Scalar>& before,
                               double time, double h, const OdeControl& control,
                               VectorT<Scalar>& after) {
  const VectorT<Scalar> f0 = rhs(time, before);
  after = before + Scalar(h) * f0;
  const Eigen::Index n = before.size();
  for (int iteration = 0; iteration < 12; ++iteration) {
    const VectorT<Scalar> f1 = rhs(time + h, after);
    const VectorT<Scalar> residual = after - before - Scalar(0.5 * h) * (f0 + f1);
    if (path_lt(scaled_error_t(residual, before, after, control),
                Scalar(0.03))) return true;
    MatrixT<Scalar> jacobian(n, n);
    for (Eigen::Index column = 0; column < n; ++column) {
      VectorT<Scalar> perturbed = after;
      const double delta = std::sqrt(std::numeric_limits<double>::epsilon()) *
        std::max(1.0, std::abs(scalar_value(after[column])));
      perturbed[column] += Scalar(delta);
      jacobian.col(column) = (rhs(time + h, perturbed) - f1) / Scalar(delta);
    }
    MatrixT<Scalar> system = MatrixT<Scalar>::Identity(n, n) - Scalar(0.5 * h) * jacobian;
    MatrixT<Scalar> rhs_matrix(n, 1);
    rhs_matrix.col(0) = -residual;
    const VectorT<Scalar> update = solve_linear(system, rhs_matrix, "ADVAN13 Newton").col(0);
    after += update;
    if (path_lt(scaled_error_t(update, before, after, control),
                Scalar(0.03))) return true;
  }
  return false;
}

template <class Scalar, class Rhs>
VectorT<Scalar> integrate_implicit_trapezoid_t(
    const Rhs& rhs, VectorT<Scalar> state, double from, double to,
    const OdeControl& control) {
  if (to <= from) return state;
  double time = from;
  const double span = to - from;
  double h = control.initial_step > 0.0 ? std::min(control.initial_step, span) : span / 10.0;
  h = std::max(h, std::min(span, 1e-8));
  int attempts = 0;
  while (time < to) {
    if (++attempts > control.max_steps) throw std::runtime_error("ADVAN13 exceeded ODE_CONTROL$max_steps.");
    h = std::min(h, to - time);
    const double minimum = 32.0 * std::numeric_limits<double>::epsilon() *
      std::max({1.0, std::abs(time), std::abs(to)});
    if (h < minimum) throw std::runtime_error("ADVAN13 ODE step size underflow.");
    VectorT<Scalar> full, half, two_half;
    const bool converged = implicit_trapezoid_step_t(rhs, state, time, h, control, full) &&
      implicit_trapezoid_step_t(rhs, state, time, h * 0.5, control, half) &&
      implicit_trapezoid_step_t(rhs, half, time + h * 0.5, h * 0.5, control, two_half);
    Scalar error_expression = Scalar(std::numeric_limits<double>::infinity());
    if (converged) error_expression = scaled_error_t(
      VectorT<Scalar>((two_half - full) / Scalar(3.0)), state, two_half, control);
    const double error = scalar_value(error_expression);
    if (converged && std::isfinite(error) &&
        path_le(error_expression, Scalar(1.0))) {
      state = two_half + (two_half - full) / Scalar(3.0);
      time += h;
    }
    const double factor = converged && error == 0.0 ? 4.0 :
      (converged && std::isfinite(error) ?
        std::clamp(0.9 * std::pow(error, -1.0 / 3.0), 0.1, 4.0) : 0.25);
    h *= factor;
  }
  return state;
}

template <class Scalar, class Rhs>
VectorT<Scalar> integrate_dopri54_interval_t(
    const Rhs& rhs, VectorT<Scalar> state, const Scalar& from,
    const Scalar& to, const OdeControl& control) {
  const Scalar span = to - from;
  const double span_value = scalar_value(span);
  if (path_le(span, Scalar(0.0))) return state;
  OdeControl normalized = control;
  if (normalized.initial_step > 0.0) {
    normalized.initial_step = std::min(1.0, normalized.initial_step / span_value);
  }
  auto scaled_rhs = [&](double fraction, const VectorT<Scalar>& value) {
    const Scalar time = from + span * Scalar(fraction);
    return VectorT<Scalar>(span * rhs(time, value));
  };
  return integrate_dopri54_t<Scalar>(scaled_rhs, state, 0.0, 1.0, normalized);
}

template <class Scalar, class Rhs>
VectorT<Scalar> integrate_implicit_interval_t(
    const Rhs& rhs, VectorT<Scalar> state, const Scalar& from,
    const Scalar& to, const OdeControl& control) {
  const Scalar span = to - from;
  const double span_value = scalar_value(span);
  if (path_le(span, Scalar(0.0))) return state;
  OdeControl normalized = control;
  if (normalized.initial_step > 0.0) {
    normalized.initial_step = std::min(1.0, normalized.initial_step / span_value);
  }
  auto scaled_rhs = [&](double fraction, const VectorT<Scalar>& value) {
    const Scalar time = from + span * Scalar(fraction);
    return VectorT<Scalar>(span * rhs(time, value));
  };
  return integrate_implicit_trapezoid_t<Scalar>(
    scaled_rhs, state, 0.0, 1.0, normalized);
}

template <class Scalar>
Scalar bioavailability_t(const ParametersT<Scalar>& p, const Rcpp::DataFrame& data,
                         int row, int cmt,
                         const DynamicDataT<Scalar>* dynamic_data = nullptr) {
  const std::string name = "F" + std::to_string(cmt);
  Scalar value = positive_parameter(p, {name.c_str()});
  if (!scalar_positive(value)) {
    value = dynamic_row_optional(data, name, row, 1.0, dynamic_data);
  }
  return scalar_positive(value) ? value : Scalar(1.0);
}

template <class Scalar>
Scalar event_infusion_rate_t(const ParametersT<Scalar>& p,
                             const Rcpp::DataFrame& data, int row, int cmt,
                             double amount, double rate_code,
                             const DynamicDataT<Scalar>* dynamic_data = nullptr) {
  if (rate_code >= 0.0) return Scalar(rate_code);
  if (rate_code != -1.0 && rate_code != -2.0) {
    throw std::domain_error("Negative RATE must be -1 (modelled Rn) or -2 (modelled Dn).");
  }
  const std::string name = std::string(rate_code == -1.0 ? "R" : "D") +
    std::to_string(cmt);
  Scalar value = positive_parameter(p, {name.c_str()});
  if (!scalar_positive(value)) {
    value = dynamic_row_optional(data, name, row, NA_REAL, dynamic_data);
  }
  if (!scalar_positive(value)) {
    throw std::domain_error("RATE=" + std::to_string(static_cast<int>(rate_code)) +
                            " requires a positive " + name + " value.");
  }
  return rate_code == -1.0 ? value : Scalar(amount) / value;
}

template <class Scalar>
Scalar observation_scale_t(const ParametersT<Scalar>& p,
                           const Rcpp::DataFrame& data, int row, int cmt,
                           const TopologyT<Scalar>& topology,
                           const DynamicDataT<Scalar>* dynamic_data = nullptr) {
  const std::string name = "S" + std::to_string(cmt);
  Scalar value = positive_parameter(p, {name.c_str()});
  if (!scalar_positive(value)) {
    value = dynamic_row_optional(data, name, row, NA_REAL, dynamic_data);
  }
  const int index = cmt - 1;
  if (!scalar_positive(value) && index >= 0 &&
      index < static_cast<int>(topology.default_scales.size())) {
    value = topology.default_scales[static_cast<std::size_t>(index)];
  }
  return scalar_positive(value) ? value : Scalar(1.0);
}

template <class Scalar>
VectorT<Scalar> infusion_input_t(int n, const std::vector<ActiveInfusionT<Scalar>>& active) {
  VectorT<Scalar> input = VectorT<Scalar>::Zero(n);
  for (const auto& infusion : active) input[infusion.compartment] += infusion.rate;
  return input;
}

template <class Scalar>
void remove_finished_t(std::vector<ActiveInfusionT<Scalar>>& active,
                       const Scalar& time) {
  active.erase(std::remove_if(active.begin(), active.end(), [&time](const auto& infusion) {
    return path_le(infusion.end, time + Scalar(1e-12));
  }), active.end());
}

template <class Scalar>
VectorT<Scalar> propagate_to_t(const ModelEngine& engine,
                              const MatrixT<Scalar>& k, VectorT<Scalar> state,
                              const Scalar& from, const Scalar& to,
                              std::vector<ActiveInfusionT<Scalar>>& active) {
  Scalar cursor = from;
  remove_finished_t(active, cursor);
  while (path_lt(cursor, to - Scalar(1e-12))) {
    Scalar segment_end = to;
    for (const auto& infusion : active) {
      if (path_gt(infusion.end, cursor + Scalar(1e-12)) &&
          path_lt(infusion.end, segment_end)) {
        segment_end = infusion.end;
      }
    }
    const AffineMapT<Scalar> map = engine_affine_map_t(
      engine, k, infusion_input_t(k.rows(), active), segment_end - cursor);
    state = map.transition * state + map.offset;
    cursor = segment_end;
    remove_finished_t(active, cursor);
  }
  return state;
}

template <class Scalar>
VectorT<Scalar> propagate_ode_to_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data,
    int row, int subject, const std::vector<Scalar>& theta,
    const std::vector<Scalar>& eta, int eta_columns,
    const std::vector<Scalar>& sigma, const ParametersT<Scalar>& parameters,
    int mixture_number, VectorT<Scalar> state, double from, double to,
    std::vector<ActiveInfusionT<Scalar>>& active,
    const DynamicDataT<Scalar>* dynamic_data = nullptr,
    DdeHistoryT<Scalar>* dde_history = nullptr) {
  Scalar cursor = Scalar(from);
  const Scalar endpoint = Scalar(to);
  remove_finished_t(active, cursor);
  while (path_lt(cursor, endpoint - Scalar(1e-12))) {
    Scalar segment_end = endpoint;
    for (const auto& infusion : active) {
      if (path_gt(infusion.end, cursor + Scalar(1e-12)) &&
          path_lt(infusion.end, segment_end)) {
        segment_end = infusion.end;
      }
    }
    const VectorT<Scalar> input = infusion_input_t(engine.n_state, active);
    if (engine.dde_enabled) {
      if (dde_history == nullptr) {
        throw std::logic_error("DDE propagation requires an initialized history.");
      }
      Scalar time = cursor;
      int steps = 0;
      while (path_lt(time, segment_end - Scalar(1e-12))) {
        if (++steps > engine.dde_max_steps) {
          throw std::runtime_error("DDE method-of-steps exceeded max_steps.");
        }
        const double remaining = scalar_value(segment_end - time);
        const double h = std::min(engine.dde_step, remaining);
        auto rhs = [&](const Scalar& stage_time, const VectorT<Scalar>& value) {
          VectorT<Scalar> lag_values(static_cast<Eigen::Index>(engine.dde_lag_inputs.size()));
          for (std::size_t lag = 0; lag < engine.dde_lag_inputs.size(); ++lag) {
            auto delay = parameters.find(engine.dde_lag_delays[lag]);
            if (delay == parameters.end() ||
                !std::isfinite(scalar_value(delay->second)) ||
                path_lt(delay->second, Scalar(engine.dde_minimum_delay)) ||
                path_lt(delay->second, Scalar(h))) {
              throw std::domain_error("DDE delay '" + engine.dde_lag_delays[lag] +
                                      "' must be finite and at least the integration step.");
            }
            lag_values[static_cast<Eigen::Index>(lag)] = dde_history->at(
              stage_time - delay->second, engine.dde_lag_states[lag]);
          }
          VectorT<Scalar> derivative = evaluate_derivatives_t(
            engine, data, row, subject, stage_time, value, parameters,
            theta, eta, eta_columns, sigma, mixture_number, dynamic_data, &lag_values);
          derivative += input;
          return derivative;
        };
        const VectorT<Scalar> k1 = rhs(time, state);
        const VectorT<Scalar> k2 = rhs(time + Scalar(0.5 * h), state + Scalar(0.5 * h) * k1);
        const VectorT<Scalar> k3 = rhs(time + Scalar(0.5 * h), state + Scalar(0.5 * h) * k2);
        const VectorT<Scalar> k4 = rhs(time + Scalar(h), state + Scalar(h) * k3);
        state += Scalar(h / 6.0) * (k1 + Scalar(2.0) * k2 + Scalar(2.0) * k3 + k4);
        time += Scalar(h);
        dde_history->append(time, state);
      }
      cursor = segment_end;
      remove_finished_t(active, cursor);
      continue;
    }
    auto rhs = [&](const Scalar& time, const VectorT<Scalar>& value) {
      VectorT<Scalar> derivative = evaluate_derivatives_t(
        engine, data, row, subject, time, value, parameters,
        theta, eta, eta_columns, sigma, mixture_number, dynamic_data);
      derivative += input;
      return derivative;
    };
    state = engine.advan == 13 ?
      integrate_implicit_interval_t(rhs, state, cursor, segment_end, engine.ode_control) :
      integrate_dopri54_interval_t(rhs, state, cursor, segment_end, engine.ode_control);
    cursor = segment_end;
    remove_finished_t(active, cursor);
  }
  return state;
}

template <class Scalar>
Scalar relative_state_change_t(const VectorT<Scalar>& before,
                               const VectorT<Scalar>& after) {
  Scalar numerator = Scalar(0.0);
  Scalar denominator = Scalar(0.0);
  for (Eigen::Index i = 0; i < before.size(); ++i) {
    const Scalar difference = after[i] - before[i];
    numerator += difference * difference;
    denominator += after[i] * after[i];
  }
  const Scalar norm = libertad::scalar_sqrt(denominator);
  const Scalar scale = libertad::choose_gt(
    norm, Scalar(1.0), norm, Scalar(1.0));
  return libertad::scalar_sqrt(numerator) / scale;
}

template <class Scalar>
std::vector<ActiveInfusionT<Scalar>> periodic_infusions_t(
    double time, const Scalar& duration, double interval, int compartment,
    const Scalar& rate) {
  const double duration_value = scalar_value(duration);
  const int previous = std::max(
    0, static_cast<int>(std::ceil(duration_value / interval - 1e-12)) - 1);
  std::vector<ActiveInfusionT<Scalar>> active;
  active.reserve(static_cast<std::size_t>(previous + 1));
  for (int dose = 0; dose <= previous; ++dose) {
    const Scalar end = Scalar(time) + duration - Scalar(dose * interval);
    if (path_gt(end, Scalar(time + 1e-12))) {
      active.push_back({end, compartment, rate});
    }
  }
  return active;
}

template <class Scalar>
VectorT<Scalar> steady_ode_bolus_post_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data,
    int row, int subject, const std::vector<Scalar>& theta,
    const std::vector<Scalar>& eta, int eta_columns,
    const std::vector<Scalar>& sigma, const ParametersT<Scalar>& parameters,
    int mixture_number, const VectorT<Scalar>& dose,
    double time, double interval,
    const DynamicDataT<Scalar>* dynamic_data = nullptr) {
  if (!(interval > 0.0)) throw std::domain_error("ODE steady-state bolus requires II > 0.");
  VectorT<Scalar> current = dose;
  const double tolerance = std::max(1e-10, engine.ode_control.rtol * 5.0);
  for (int iteration = 0; iteration < 10000; ++iteration) {
    std::vector<ActiveInfusionT<Scalar>> active;
    VectorT<Scalar> next = propagate_ode_to_t(
      engine, data, row, subject, theta, eta, eta_columns, sigma,
      parameters, mixture_number, current, time, time + interval, active,
      dynamic_data) + dose;
    if (path_le(relative_state_change_t(current, next), Scalar(tolerance))) {
      return next;
    }
    current = next;
  }
  throw std::runtime_error("ODE bolus periodic shooting did not converge.");
}

template <class Scalar>
VectorT<Scalar> steady_ode_infusion_pre_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data,
    int row, int subject, const std::vector<Scalar>& theta,
    const std::vector<Scalar>& eta, int eta_columns,
    const std::vector<Scalar>& sigma, const ParametersT<Scalar>& parameters,
    int mixture_number, int compartment, const Scalar& administered_rate,
    const Scalar& duration, double time, double interval,
    const DynamicDataT<Scalar>* dynamic_data = nullptr) {
  if (!path_gt(duration, Scalar(0.0)) || !(interval > 0.0)) {
    throw std::domain_error("ODE steady-state infusion requires duration and II > 0.");
  }
  VectorT<Scalar> current = VectorT<Scalar>::Zero(engine.n_state);
  const double tolerance = std::max(1e-10, engine.ode_control.rtol * 5.0);
  for (int iteration = 0; iteration < 10000; ++iteration) {
    std::vector<ActiveInfusionT<Scalar>> active = periodic_infusions_t(
      time, duration, interval, compartment, administered_rate);
    VectorT<Scalar> next = propagate_ode_to_t(
      engine, data, row, subject, theta, eta, eta_columns, sigma,
      parameters, mixture_number, current, time, time + interval, active,
      dynamic_data);
    if (path_le(relative_state_change_t(current, next), Scalar(tolerance))) {
      return next;
    }
    current = next;
  }
  throw std::runtime_error("ODE infusion periodic shooting did not converge.");
}

template <class Scalar>
VectorT<Scalar> steady_bolus_post_t(const ModelEngine& engine,
                                    const MatrixT<Scalar>& k,
                                    const VectorT<Scalar>& dose, double interval) {
  if (!(interval > 0.0) || !std::isfinite(interval)) {
    throw std::domain_error("Steady-state bolus requires II > 0.");
  }
  return solve_periodic_t(
    engine_transition_t(engine, k, Scalar(interval)), dose, "Bolus");
}

template <class Scalar>
VectorT<Scalar> steady_infusion_pre_t(const ModelEngine& engine,
                                      const MatrixT<Scalar>& k,
                                      const VectorT<Scalar>& rate,
                                      const Scalar& duration, double interval) {
  const double duration_value = scalar_value(duration);
  if (!(duration_value > 0.0) || !(interval > 0.0)) {
    throw std::domain_error("Steady-state infusion requires duration and II > 0.");
  }
  const int complete = static_cast<int>(std::floor(duration_value / interval + 1e-12));
  Scalar remainder = duration - Scalar(complete * interval);
  if (path_lt(remainder, Scalar(1e-12))) remainder = Scalar(0.0);
  const VectorT<Scalar> baseline = Scalar(static_cast<double>(complete)) * rate;
  VectorT<Scalar> first_input = baseline;
  if (path_gt(remainder, Scalar(0.0))) first_input += rate;
  AffineMapT<Scalar> on = engine_affine_map_t(engine, k, first_input, remainder);
  AffineMapT<Scalar> off = engine_affine_map_t(
    engine, k, baseline, Scalar(interval) - remainder);
  return solve_periodic_t(
    MatrixT<Scalar>(off.transition * on.transition),
    VectorT<Scalar>(off.transition * on.offset + off.offset), "Infusion");
}

template <class Scalar>
std::vector<Scalar> simulate_analytical_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data,
    const std::vector<Scalar>& theta, const std::vector<Scalar>& eta,
    const std::vector<Scalar>& sigma,
    const std::vector<int>& mixture_assignment = std::vector<int>(),
    const DynamicDataT<Scalar>* dynamic_data = nullptr) {
  const int n_rows = data.nrows();
  Rcpp::NumericVector time = data["TIME"];
  Rcpp::NumericVector amount = data["AMT"];
  Rcpp::NumericVector rate = data["RATE"];
  Rcpp::NumericVector interval = data["II"];
  Rcpp::IntegerVector evid = data["EVID"];
  Rcpp::IntegerVector cmt = data["CMT"];
  Rcpp::IntegerVector ss = data["SS"];
  Rcpp::IntegerVector eta_subject_index = data[".ID_INDEX"];
  Rcpp::IntegerVector subject_index = data.containsElementNamed(".STRUCT_ID_INDEX") ?
    Rcpp::IntegerVector(data[".STRUCT_ID_INDEX"]) : eta_subject_index;
  int n_subjects = 0;
  for (int value : eta_subject_index) n_subjects = std::max(n_subjects, value);
  if (n_subjects < 1 || eta.size() % static_cast<std::size_t>(n_subjects) != 0U) {
    throw std::invalid_argument("ETA vector cannot be divided into subject rows.");
  }
  const int eta_columns = static_cast<int>(eta.size() / static_cast<std::size_t>(n_subjects));
  const int n_state = engine.n_state;
  std::vector<Scalar> prediction(static_cast<std::size_t>(n_rows), Scalar(0.0));
  VectorT<Scalar> state = VectorT<Scalar>::Zero(n_state);
  std::vector<ActiveInfusionT<Scalar>> active;
  DdeHistoryT<Scalar> dde_history;
  MatrixT<Scalar> previous_k = MatrixT<Scalar>::Zero(n_state, n_state);
  ParametersT<Scalar> previous_parameters;
  int previous_row = -1;
  int previous_mixture = 1;
  int previous_subject = -1;
  double previous_time = 0.0;
  bool have_previous = false;

  for (int row = 0; row < n_rows; ++row) {
    const int structural_subject = subject_index[row] - 1;
    const int subject = eta_subject_index[row] - 1;
    if (structural_subject != previous_subject) {
      state.setZero();
      active.clear();
      have_previous = false;
      previous_time = time[row];
      previous_subject = structural_subject;
      if (engine.dde_enabled) dde_history.reset(Scalar(time[row]), state, engine.dde_history);
    }
    if (have_previous) {
      if (time[row] < previous_time - 1e-12) throw std::domain_error("Subject event times decrease.");
      if (engine.is_ode()) {
        state = propagate_ode_to_t(
          engine, data, previous_row, subject, theta, eta, eta_columns, sigma,
          previous_parameters, previous_mixture, state,
          previous_time, time[row], active, dynamic_data,
          engine.dde_enabled ? &dde_history : nullptr);
      } else {
        state = propagate_to_t(
          engine, previous_k, state, Scalar(previous_time), Scalar(time[row]), active);
      }
    }
    const int mixture_number = mixture_assignment.empty() ?
      static_cast<int>(row_optional(data, "MIXNUM", row, 1.0)) :
      mixture_assignment.at(static_cast<std::size_t>(subject));
    ParametersT<Scalar> parameters = evaluate_parameters_t(
      engine, data, row, subject, theta, eta, eta_columns, sigma,
      mixture_number, dynamic_data);
    TopologyT<Scalar> topology;
    if (engine.is_ode()) {
      topology.k = MatrixT<Scalar>::Zero(n_state, n_state);
      topology.default_scales.assign(static_cast<std::size_t>(n_state), Scalar(1.0));
    } else {
      topology = engine.matrix_graph.enabled ?
        build_graph_topology_t(engine.matrix_graph, parameters) :
        build_topology_t(engine.advan, parameters);
    }
    if (topology.k.rows() != n_state) throw std::logic_error("ADVAN state dimension changed.");
    if (evid[row] == 3 || evid[row] == 4) {
      state.setZero();
      active.clear();
      if (engine.dde_enabled) dde_history.reset(Scalar(time[row]), state, engine.dde_history);
    }
    const bool dosing = amount[row] > 0.0 &&
      (evid[row] == 1 || evid[row] == 4 || evid[row] == 0);
    if (dosing) {
      const int dose_cmt = cmt[row] > 0 ? cmt[row] : engine.dose_cmp;
      const int dose_index = compartment_index(dose_cmt, topology.default_dose, n_state);
      const Scalar f = bioavailability_t(
        parameters, data, row, dose_cmt, dynamic_data);
      const Scalar event_rate = event_infusion_rate_t(
        parameters, data, row, dose_cmt, amount[row], rate[row], dynamic_data);
      const int ss_flag = ss[row] != 0 ? ss[row] : engine.model_ss;
      if (ss_flag == 1) {
        state.setZero();
        active.clear();
      } else if (ss_flag != 0 && ss_flag != 2) {
        throw std::domain_error("Only SS=0, SS=1, and SS=2 are supported.");
      }
      if (engine.dde_enabled && ss_flag != 0) {
        throw std::domain_error("Experimental DDE models currently require SS=0; provide an explicit warm-up regimen.");
      }
      if (scalar_positive(event_rate)) {
        const Scalar duration = Scalar(amount[row]) / event_rate;
        VectorT<Scalar> input = VectorT<Scalar>::Zero(n_state);
        input[dose_index] = event_rate * f;
        if (ss_flag != 0) {
          VectorT<Scalar> periodic = engine.is_ode() ? steady_ode_infusion_pre_t(
            engine, data, row, subject, theta, eta, eta_columns, sigma,
            parameters, mixture_number, dose_index, event_rate * f,
            duration, time[row], interval[row], dynamic_data) :
            steady_infusion_pre_t(engine, topology.k, input, duration, interval[row]);
          if (ss_flag == 1) state = periodic; else state += periodic;
          std::vector<ActiveInfusionT<Scalar>> periodic_active = periodic_infusions_t(
            time[row], duration, interval[row], dose_index, event_rate * f);
          active.insert(active.end(), periodic_active.begin(), periodic_active.end());
        } else {
          active.push_back({Scalar(time[row]) + duration, dose_index, event_rate * f});
        }
      } else {
        VectorT<Scalar> dose = VectorT<Scalar>::Zero(n_state);
        dose[dose_index] = Scalar(amount[row]) * f;
        if (ss_flag != 0) {
          VectorT<Scalar> periodic = engine.is_ode() ? steady_ode_bolus_post_t(
            engine, data, row, subject, theta, eta, eta_columns, sigma,
            parameters, mixture_number, dose, time[row], interval[row],
            dynamic_data) :
            steady_bolus_post_t(engine, topology.k, dose, interval[row]);
          if (ss_flag == 1) state = periodic; else state += periodic;
        } else {
          state += dose;
        }
      }
      if (engine.dde_enabled) dde_history.append(Scalar(time[row]), state);
    }
    const int observation_cmt = cmt[row] > 0 && evid[row] == 0 ? cmt[row] : engine.obs_cmp;
    const int observation_index = compartment_index(
      observation_cmt, topology.default_observation, n_state);
    const Scalar scale = observation_scale_t(
      parameters, data, row, observation_cmt, topology, dynamic_data);
    const auto direct_prediction = parameters.find("F");
    prediction[static_cast<std::size_t>(row)] = direct_prediction == parameters.end() ?
      state[observation_index] / scale : direct_prediction->second;
    previous_k = topology.k;
    previous_parameters = std::move(parameters);
    previous_row = row;
    previous_mixture = mixture_number;
    previous_time = time[row];
    have_previous = true;
  }
  return prediction;
}

struct ObjectiveTape {
  CppAD::ADFun<double> fun;
  std::vector<std::string> domain_names;
  std::vector<std::string> dynamic_columns;
  std::vector<int> dynamic_observed_rows;
  std::vector<int> structural_dvid;
  std::vector<double> dynamic_values;
  int n_rows = 0;
  libertad::SparseHessianCache hessian_cache;
};

class TapePathChange : public std::runtime_error {
 public:
  explicit TapePathChange(const std::string& context,
                          std::vector<double> point = std::vector<double>())
      : std::runtime_error("CppAD tape path changed in " + context +
                           "; automatic retaping is required."),
        point_(std::move(point)) {}

  const std::vector<double>& point() const { return point_; }

 private:
  std::vector<double> point_;
};

void require_unchanged_path(CppAD::ADFun<double>& fun,
                            const std::string& context) {
  if (fun.compare_change_number() != 0U) throw TapePathChange(context);
}

struct PredictionTape {
  CppAD::ADFun<double> fun;
  std::vector<std::string> domain_names;
  std::vector<std::string> dynamic_columns;
  std::vector<double> dynamic_values;
  int n_rows = 0;
  std::string propagation_kernel;
  std::size_t operation_count = 0;
  std::size_t variable_count = 0;
  std::string derivative_strategy = "not-evaluated";
  std::size_t jacobian_nonzeros = 0;
};

bool structural_data_input(const std::string& name) {
  static const std::unordered_map<std::string, bool> structural = {
    {"ID", true}, {"TIME", true}, {"AMT", true}, {"RATE", true},
    {"II", true}, {"ADDL", true}, {"EVID", true}, {"CMT", true},
    {"SS", true}, {"MIXNUM", true}, {"DVID", true}, {"DV", true},
    {"MDV", true}, {"LLOQ", true}, {"BLQ", true}, {"CENS", true},
    {".ID_INDEX", true}, {".OCC_INDEX", true}
  };
  return structural.find(name) != structural.end();
}

bool data_backed_model_input(const std::string& name) {
  if (structural_data_input(name) || name == "F" || name == "T" ||
      name == "MIXNUM" || starts_with(name, "THETA_") ||
      starts_with(name, "ETA_") || starts_with(name, "SIGMA_") ||
      starts_with(name, "ERR_") || starts_with(name, "A_")) {
    return false;
  }
  return true;
}

std::vector<std::string> prediction_dynamic_columns(
    const ModelEngine& engine, const Rcpp::DataFrame& data) {
  std::vector<std::string> columns;
  std::unordered_map<std::string, bool> seen;
  auto append = [&](const std::vector<std::string>& inputs) {
    for (const std::string& name : inputs) {
      if (!data_backed_model_input(name) ||
          !data.containsElementNamed(name.c_str()) || seen[name]) continue;
      for (int row = 0; row < data.nrows(); ++row) {
        if (!std::isfinite(data_value(data, name, row))) {
          throw std::domain_error("Dynamic model input '" + name +
                                  "' contains a non-finite value.");
        }
      }
      seen[name] = true;
      columns.push_back(name);
    }
  };
  append(engine.pred->input_names);
  if (engine.des) append(engine.des->input_names);
  return columns;
}

std::vector<double> prediction_dynamic_values(
    const std::vector<std::string>& columns, const Rcpp::DataFrame& data,
    int expected_rows = -1) {
  if (expected_rows >= 0 && data.nrows() != expected_rows) {
    throw std::invalid_argument("Dynamic prediction data has a different row count.");
  }
  std::vector<double> values;
  values.reserve(columns.size() * static_cast<std::size_t>(data.nrows()));
  for (const std::string& name : columns) {
    if (!data.containsElementNamed(name.c_str())) {
      throw std::invalid_argument("Dynamic prediction data is missing column '" + name + "'.");
    }
    for (int row = 0; row < data.nrows(); ++row) {
      const double value = data_value(data, name, row);
      if (!std::isfinite(value)) {
        throw std::domain_error("Dynamic prediction input '" + name +
                                "' contains a non-finite value.");
      }
      values.push_back(value);
    }
  }
  return values;
}

std::vector<int> fo_observed_rows(const Rcpp::DataFrame& data) {
  Rcpp::NumericVector dv = data["DV"];
  Rcpp::NumericVector evid = data["EVID"];
  Rcpp::NumericVector mdv = data["MDV"];
  std::vector<int> observed;
  for (int row = 0; row < data.nrows(); ++row) {
    if (evid[row] == 0.0 && mdv[row] == 0.0 && std::isfinite(dv[row])) {
      observed.push_back(row);
    }
  }
  return observed;
}

std::vector<int> fo_dvid_values(const Rcpp::DataFrame& data) {
  std::vector<int> result(static_cast<std::size_t>(data.nrows()), 1);
  if (!data.containsElementNamed("DVID")) return result;
  Rcpp::NumericVector dvid = data["DVID"];
  for (int row = 0; row < data.nrows(); ++row) {
    result[static_cast<std::size_t>(row)] =
      std::max(1, static_cast<int>(dvid[row]));
  }
  return result;
}

std::vector<double> fo_dynamic_values(const ObjectiveTape& tape,
                                      const Rcpp::DataFrame& data) {
  if (tape.n_rows != data.nrows()) {
    throw std::invalid_argument("A shared FO tape received a different number of rows.");
  }
  if (fo_observed_rows(data) != tape.dynamic_observed_rows) {
    throw std::invalid_argument("A shared FO tape received a different observation pattern.");
  }
  if (fo_dvid_values(data) != tape.structural_dvid) {
    throw std::invalid_argument("A shared FO tape received a different DVID pattern.");
  }
  std::vector<double> values = prediction_dynamic_values(
    tape.dynamic_columns, data, tape.n_rows);
  Rcpp::NumericVector dv = data["DV"];
  values.reserve(values.size() + tape.dynamic_observed_rows.size());
  for (int row : tape.dynamic_observed_rows) {
    const double value = dv[row];
    if (!std::isfinite(value)) {
      throw std::domain_error("A shared FO tape received a non-finite observation.");
    }
    values.push_back(value);
  }
  return values;
}

void set_fo_dynamic(ObjectiveTape& tape, const Rcpp::DataFrame& data) {
  const std::vector<double> values = fo_dynamic_values(tape, data);
  if (values.size() != tape.fun.size_dyn_ind()) {
    throw std::logic_error("Shared FO dynamic data do not match the recorded tape.");
  }
  if (!values.empty()) tape.fun.new_dynamic(values);
  tape.dynamic_values = values;
}

std::vector<double> flatten_parameters(const Rcpp::NumericVector& theta,
                                       const Rcpp::NumericMatrix& eta,
                                       const Rcpp::NumericVector& sigma) {
  std::vector<double> result;
  result.reserve(theta.size() + eta.size() + sigma.size());
  for (double value : theta) result.push_back(value);
  for (int row = 0; row < eta.nrow(); ++row) {
    for (int column = 0; column < eta.ncol(); ++column) result.push_back(eta(row, column));
  }
  for (double value : sigma) result.push_back(value);
  return result;
}

std::vector<std::string> parameter_names(int n_theta, int n_subjects,
                                         int n_eta, int n_sigma) {
  std::vector<std::string> names;
  for (int i = 0; i < n_theta; ++i) names.push_back("THETA_" + std::to_string(i + 1));
  for (int subject = 0; subject < n_subjects; ++subject) {
    for (int i = 0; i < n_eta; ++i) {
      names.push_back("ETA_" + std::to_string(subject + 1) + "_" + std::to_string(i + 1));
    }
  }
  for (int i = 0; i < n_sigma; ++i) names.push_back("SIGMA_" + std::to_string(i + 1));
  return names;
}

std::unique_ptr<PredictionTape> record_prediction_tape(
    const ModelEngine& engine, const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta, const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma) {
  const int minimum_eta_columns = required_eta_columns(engine, data);
  const int between_eta = engine.n_eta - engine.iov;
  if (theta.size() != engine.n_theta || eta.ncol() < minimum_eta_columns ||
      (!engine.re_enabled && engine.iov > 0 &&
       (eta.ncol() - between_eta) % engine.iov != 0)) {
    throw std::invalid_argument("Prediction tape parameter dimensions are inconsistent with the model.");
  }
  std::vector<double> point = flatten_parameters(theta, eta, sigma);
  std::vector<CppAD::AD<double>> independent(point.begin(), point.end());
  const std::vector<std::string> dynamic_columns =
    prediction_dynamic_columns(engine, data);
  const std::vector<double> dynamic_values =
    prediction_dynamic_values(dynamic_columns, data);
  std::vector<CppAD::AD<double>> dynamic(dynamic_values.begin(), dynamic_values.end());
  if (dynamic.empty()) CppAD::Independent(independent);
  else CppAD::Independent(independent, dynamic);
  DynamicDataT<CppAD::AD<double>> dynamic_data;
  dynamic_data.n_rows = data.nrows();
  dynamic_data.values = dynamic;
  for (std::size_t column = 0; column < dynamic_columns.size(); ++column) {
    dynamic_data.column_positions[dynamic_columns[column]] = column;
  }
  std::size_t cursor = 0;
  std::vector<CppAD::AD<double>> theta_ad(static_cast<std::size_t>(theta.size()));
  for (auto& value : theta_ad) value = independent[cursor++];
  std::vector<CppAD::AD<double>> eta_ad(static_cast<std::size_t>(eta.size()));
  for (auto& value : eta_ad) value = independent[cursor++];
  std::vector<CppAD::AD<double>> sigma_ad(static_cast<std::size_t>(sigma.size()));
  for (auto& value : sigma_ad) value = independent[cursor++];
  std::vector<CppAD::AD<double>> predictions = simulate_analytical_t(
    engine, data, theta_ad, eta_ad, sigma_ad, std::vector<int>(),
    dynamic.empty() ? nullptr : &dynamic_data);
  auto tape = std::make_unique<PredictionTape>();
  tape->fun.Dependent(independent, predictions);
  tape->fun.optimize();
  tape->operation_count = tape->fun.size_op();
  tape->variable_count = tape->fun.size_var();
  tape->domain_names = parameter_names(theta.size(), eta.nrow(), eta.ncol(), sigma.size());
  tape->dynamic_columns = dynamic_columns;
  tape->dynamic_values = dynamic_values;
  tape->n_rows = data.nrows();
  tape->propagation_kernel = propagation_kernel_name(engine);
  return tape;
}

std::vector<double> prediction_point(PredictionTape& tape,
                                     const Rcpp::NumericVector& point) {
  if (point.size() != static_cast<R_xlen_t>(tape.domain_names.size())) {
    throw std::invalid_argument("Prediction tape point has the wrong length.");
  }
  return Rcpp::as<std::vector<double>>(point);
}

template <class Scalar>
Scalar scalar_floor_t(const Scalar& value, double floor) {
  return libertad::choose_gt(value, Scalar(floor), value, Scalar(floor));
}

inline double scalar_erf_t(double value) { return std::erf(value); }
inline CppAD::AD<double> scalar_erf_t(const CppAD::AD<double>& value) {
  return CppAD::erf(value);
}

template <class Scalar>
Scalar normal_cdf_t(const Scalar& value) {
  const Scalar probability = Scalar(0.5) *
    (Scalar(1.0) + scalar_erf_t(value * Scalar(0.7071067811865475244)));
  return scalar_floor_t(probability, 1e-300);
}

template <class Scalar>
Scalar residual_variance_t(const ModelEngine& engine,
                           const Scalar& prediction,
                           const std::vector<Scalar>& sigma,
                           int dvid) {
  int per_response = engine.error_type == "combined" ? 2 :
    (engine.error_type == "power" ? 2 : 1);
  int offset = std::max(0, dvid - 1) * per_response;
  if (offset + per_response > static_cast<int>(sigma.size())) offset = 0;
  if (sigma.empty()) {
    throw std::domain_error("Likelihood evaluation requires residual SIGMA parameters.");
  }
  const Scalar s1 = sigma[static_cast<std::size_t>(offset)];
  const auto sigma_variance = [&](const Scalar& value) {
    return engine.sigma_parameterization == "variance" ? value : value * value;
  };
  Scalar variance;
  if (engine.error_type == "additive" || engine.error_type == "exponential") {
    variance = sigma_variance(s1);
  } else if (engine.error_type == "proportional") {
    variance = sigma_variance(s1) * prediction * prediction;
  } else if (engine.error_type == "power") {
    if (sigma.size() <= static_cast<std::size_t>(offset + 1)) {
      throw std::domain_error("Power residual error requires two SIGMA parameters.");
    }
    const Scalar magnitude = scalar_floor_t(libertad::scalar_abs(prediction), 1e-12);
    variance = sigma_variance(s1) * libertad::scalar_pow(
      magnitude, Scalar(2.0) * sigma[static_cast<std::size_t>(offset + 1)]);
  } else {
    if (sigma.size() <= static_cast<std::size_t>(offset + 1)) {
      throw std::domain_error("Combined residual error requires two SIGMA parameters.");
    }
    const Scalar s2 = sigma[static_cast<std::size_t>(offset + 1)];
    variance = sigma_variance(s1) * prediction * prediction + sigma_variance(s2);
  }
  return scalar_floor_t(variance, 1e-16);
}

template <class Scalar>
MatrixT<Scalar> omega_matrix_t(const ModelEngine& engine,
                               const std::vector<Scalar>& omega) {
  if (omega.size() != engine.omega_rows.size()) {
    throw std::invalid_argument("OMEGA parameter vector has the wrong length.");
  }
  MatrixT<Scalar> covariance = MatrixT<Scalar>::Zero(engine.n_eta, engine.n_eta);
  for (std::size_t i = 0; i < omega.size(); ++i) {
    const int row = engine.omega_rows[i];
    const int column = engine.omega_cols[i];
    if (row < 0 || column < 0 || row >= engine.n_eta || column >= engine.n_eta) {
      throw std::logic_error("OMEGA covariance index is outside the ETA dimension.");
    }
    covariance(row, column) = omega[i];
    covariance(column, row) = omega[i];
  }
  return covariance;
}

template <class Scalar>
MatrixT<Scalar> expanded_omega_t(const ModelEngine& engine,
                                 const Rcpp::DataFrame& data,
                                 const MatrixT<Scalar>& base,
                                 int expanded_dimension) {
  if (engine.re_enabled) {
    MatrixT<Scalar> output = MatrixT<Scalar>::Zero(
      expanded_dimension, expanded_dimension);
    int offset = 0;
    for (std::size_t block_index = 0; block_index < engine.re_blocks.size(); ++block_index) {
      const std::vector<int>& indices = engine.re_blocks[block_index];
      const std::string total_name = ".RE_TOTAL_" + std::to_string(block_index + 1U);
      const int units = static_cast<int>(data_value(data, total_name, 0));
      for (int unit = 0; unit < units; ++unit) {
        for (std::size_t row = 0; row < indices.size(); ++row) {
          for (std::size_t column = 0; column < indices.size(); ++column) {
            output(offset + unit * static_cast<int>(indices.size()) + static_cast<int>(row),
                   offset + unit * static_cast<int>(indices.size()) + static_cast<int>(column)) =
              base(indices[row], indices[column]);
          }
        }
      }
      offset += units * static_cast<int>(indices.size());
    }
    if (offset != expanded_dimension) {
      throw std::invalid_argument("Expanded random-effect covariance has the wrong dimension.");
    }
    return output;
  }
  if (engine.iov == 0) return base;
  const int between = engine.n_eta - engine.iov;
  if (expanded_dimension < between ||
      (expanded_dimension - between) % engine.iov != 0) {
    throw std::invalid_argument("Expanded IOV covariance has an invalid layout.");
  }
  const int occasions = (expanded_dimension - between) / engine.iov;
  MatrixT<Scalar> output = MatrixT<Scalar>::Zero(
    expanded_dimension, expanded_dimension);
  if (between) output.topLeftCorner(between, between) =
    base.topLeftCorner(between, between);
  for (int occasion = 0; occasion < occasions; ++occasion) {
    const int target = between + occasion * engine.iov;
    output.block(target, target, engine.iov, engine.iov) =
      base.bottomRightCorner(engine.iov, engine.iov);
  }
  return output;
}

template <class Scalar>
Scalar omega_subject_prior_t(const MatrixT<Scalar>& covariance,
                             const VectorT<Scalar>& eta) {
  const Eigen::Index n = covariance.rows();
  MatrixT<Scalar> lower = MatrixT<Scalar>::Zero(n, n);
  Scalar logdet = Scalar(0.0);
  for (Eigen::Index row = 0; row < n; ++row) {
    for (Eigen::Index column = 0; column <= row; ++column) {
      Scalar value = covariance(row, column);
      for (Eigen::Index k = 0; k < column; ++k) {
        value -= lower(row, k) * lower(column, k);
      }
      if (row == column) {
        if (!(scalar_value(value) > 1e-14)) {
          throw std::domain_error("OMEGA matrix is not positive definite at the recording point.");
        }
        lower(row, column) = libertad::scalar_sqrt(value);
        logdet += Scalar(2.0) * libertad::scalar_log(lower(row, column));
      } else {
        lower(row, column) = value / lower(column, column);
      }
    }
  }
  VectorT<Scalar> standardized(n);
  for (Eigen::Index row = 0; row < n; ++row) {
    Scalar value = eta[row];
    for (Eigen::Index column = 0; column < row; ++column) {
      value -= lower(row, column) * standardized[column];
    }
    standardized[row] = value / lower(row, row);
  }
  Scalar quadratic = Scalar(0.0);
  for (Eigen::Index i = 0; i < n; ++i) quadratic += standardized[i] * standardized[i];
  return logdet + quadratic;
}

template <class Scalar>
std::vector<Scalar> evaluate_error_outputs_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data, int row, int subject,
    const std::vector<Scalar>& theta, const std::vector<Scalar>& eta,
    int eta_columns, const std::vector<Scalar>& sigma, int mixture_number,
    const Scalar& prediction, bool first, double previous_dv,
    double previous_time, const std::vector<std::size_t>& selected,
    const char* context,
    const ParametersT<Scalar>* input_overrides = nullptr) {
  if (!engine.error || selected.empty()) {
    throw std::logic_error(std::string("Compiled ") + context + " outputs are missing.");
  }
  const double current_time = row_optional(data, "TIME", row, 0.0);
  const double current_dv = row_optional(data, "DV", row, NA_REAL);
  ParametersT<Scalar> parameters = evaluate_parameters_t(
    engine, data, row, subject, theta, eta, eta_columns, sigma,
    mixture_number);
  std::vector<Scalar> inputs(engine.error->input_names.size(), Scalar(0.0));
  for (std::size_t index = 0; index < engine.error->input_names.size(); ++index) {
    const std::string& name = engine.error->input_names[index];
    int parameter_index = indexed_name(name, "THETA_");
    if (parameter_index >= 0) {
      if (parameter_index >= static_cast<int>(theta.size())) {
        throw std::out_of_range("THETA index exceeds values in user likelihood.");
      }
      inputs[index] = theta[static_cast<std::size_t>(parameter_index)];
      continue;
    }
    parameter_index = indexed_name(name, "ETA_");
    if (parameter_index >= 0) {
      const int column = eta_column(
        engine, data, row, parameter_index, eta_columns);
      const std::size_t position =
        static_cast<std::size_t>(subject * eta_columns + column);
      if (position >= eta.size()) {
        throw std::out_of_range("ETA index exceeds values in user likelihood.");
      }
      inputs[index] = eta[position];
      continue;
    }
    parameter_index = indexed_name(name, "SIGMA_");
    if (parameter_index >= 0) {
      if (parameter_index >= static_cast<int>(sigma.size())) {
        throw std::out_of_range("SIGMA index exceeds values in user likelihood.");
      }
      inputs[index] = sigma[static_cast<std::size_t>(parameter_index)];
      continue;
    }
    if (name == "F" || name == "PRED" || name == "IPRED") {
      inputs[index] = prediction;
      continue;
    }
    if (name == "DV") {
      inputs[index] = Scalar(current_dv);
      continue;
    }
    if (name == "PREV_DV") {
      inputs[index] = Scalar(previous_dv);
      continue;
    }
    if (name == "PREV_TIME") {
      inputs[index] = Scalar(previous_time);
      continue;
    }
    if (name == "DT") {
      inputs[index] = Scalar(first ? 0.0 : current_time - previous_time);
      continue;
    }
    if (name == "FIRST") {
      inputs[index] = Scalar(first ? 1.0 : 0.0);
      continue;
    }
    if (name == "MIXNUM") {
      inputs[index] = Scalar(mixture_number);
      continue;
    }
    if (input_overrides != nullptr) {
      const auto overridden = input_overrides->find(name);
      if (overridden != input_overrides->end()) {
        inputs[index] = overridden->second;
        continue;
      }
    }
    const auto assigned = parameters.find(name);
    if (assigned != parameters.end()) {
      inputs[index] = assigned->second;
      continue;
    }
    inputs[index] = Scalar(data_value(data, name, row));
    if (!std::isfinite(scalar_value(inputs[index]))) {
      throw std::domain_error("User likelihood input '" + name +
                              "' is non-finite at row " +
                              std::to_string(row + 1) + ".");
    }
  }
  const std::vector<Scalar> output = engine.error->eval_outputs(inputs, selected);
  for (const Scalar& value : output) {
    if (!std::isfinite(scalar_value(value))) {
      throw std::domain_error(std::string(context) +
                              " returned a non-finite value at row " +
                              std::to_string(row + 1) + ".");
    }
  }
  return output;
}

template <class Scalar>
Scalar user_likelihood_nll_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data, int row, int subject,
    const std::vector<Scalar>& theta, const std::vector<Scalar>& eta,
    int eta_columns, const std::vector<Scalar>& sigma, int mixture_number,
    const Scalar& prediction, bool first, double previous_dv,
    double previous_time) {
  if (engine.likelihood_output.size() != 1U) {
    throw std::logic_error("Compiled user likelihood is missing.");
  }
  const std::vector<Scalar> output = evaluate_error_outputs_t(
    engine, data, row, subject, theta, eta, eta_columns, sigma,
    mixture_number, prediction, first, previous_dv, previous_time,
    engine.likelihood_output, "User likelihood");
  if (engine.likelihood_scale == "log") return Scalar(-2.0) * output[0];
  if (!(scalar_value(output[0]) > 0.0)) {
    throw std::domain_error("LIK must be positive at row " +
                            std::to_string(row + 1) + ".");
  }
  return Scalar(-2.0) * libertad::scalar_log(
    scalar_floor_t(output[0], 1e-300));
}

template <class Scalar>
Scalar log_sum_exp_t(const std::vector<Scalar>& values) {
  if (values.empty()) throw std::invalid_argument("log-sum-exp requires values.");
  Scalar maximum = values.front();
  for (std::size_t index = 1; index < values.size(); ++index) {
    maximum = libertad::choose_gt(values[index], maximum, values[index], maximum);
  }
  Scalar total = Scalar(0.0);
  for (const Scalar& value : values) total += libertad::scalar_exp(value - maximum);
  return maximum + libertad::scalar_log(total);
}

template <class Scalar>
std::vector<Scalar> hmm_normalize_weights_t(
    const std::vector<Scalar>& values, const std::string& scale,
    const char* label, int row) {
  if (values.empty()) throw std::invalid_argument("Hidden Markov weights are empty.");
  std::vector<Scalar> result(values.size());
  if (scale == "log") {
    const Scalar normalizer = log_sum_exp_t(values);
    for (std::size_t index = 0; index < values.size(); ++index) {
      result[index] = libertad::scalar_exp(values[index] - normalizer);
    }
    return result;
  }
  Scalar total = Scalar(0.0);
  for (const Scalar& value : values) {
    if (scalar_value(value) < 0.0) {
      throw std::domain_error(std::string("HMM ") + label +
                              " weights must be non-negative at row " +
                              std::to_string(row + 1) + ".");
    }
    total += value;
  }
  if (!(scalar_value(total) > 0.0)) {
    throw std::domain_error(std::string("HMM ") + label +
                            " weights sum to zero at row " +
                            std::to_string(row + 1) + ".");
  }
  for (std::size_t index = 0; index < values.size(); ++index) {
    result[index] = values[index] / total;
  }
  return result;
}

template <class Scalar>
struct HmmRowComponents {
  std::vector<Scalar> initial;
  std::vector<std::vector<Scalar>> transition;
  std::vector<Scalar> emission;
  std::vector<Scalar> log_emission;
};

template <class Scalar>
HmmRowComponents<Scalar> hmm_row_components_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data, int row, int subject,
    const std::vector<Scalar>& theta, const std::vector<Scalar>& eta,
    int eta_columns, const std::vector<Scalar>& sigma, int mixture_number,
    const Scalar& prediction, bool first, double previous_dv,
    double previous_time) {
  const int states = engine.hmm_states;
  const std::size_t state_count = static_cast<std::size_t>(states);
  const std::size_t transition_count = engine.hmm_continuous ?
    state_count * (state_count - 1U) : state_count * state_count;
  const std::size_t expected = state_count + transition_count + state_count;
  if (!engine.hmm_enabled || states < 2 || engine.hmm_outputs.size() != expected) {
    throw std::logic_error("Compiled hidden Markov likelihood is inconsistent.");
  }
  const std::vector<Scalar> output = evaluate_error_outputs_t(
    engine, data, row, subject, theta, eta, eta_columns, sigma,
    mixture_number, prediction, first, previous_dv, previous_time,
    engine.hmm_outputs, "Hidden Markov likelihood");
  HmmRowComponents<Scalar> result;
  if (first) {
    result.initial = hmm_normalize_weights_t(
      std::vector<Scalar>(output.begin(), output.begin() + states),
      engine.hmm_initial_scale, "initial", row);
  } else {
    const std::size_t transition_offset = state_count;
    result.transition.resize(state_count);
    if (engine.hmm_continuous) {
      MatrixT<Scalar> generator = MatrixT<Scalar>::Zero(states, states);
      std::size_t rate_index = transition_offset;
      for (int from = 0; from < states; ++from) {
        for (int to = 0; to < states; ++to) {
          if (from == to) continue;
          Scalar rate = output[rate_index++];
          if (engine.hmm_rate_scale == "log") rate = libertad::scalar_exp(rate);
          if (scalar_value(rate) < 0.0) {
            throw std::domain_error(
              "Continuous-time HMM rates must be non-negative at row " +
              std::to_string(row + 1) + ".");
          }
          generator(from, to) = rate;
          generator(from, from) -= rate;
        }
      }
      const double current_time = row_optional(data, "TIME", row, 0.0);
      const double dt = current_time - previous_time;
      if (!std::isfinite(dt) || dt < 0.0) {
        throw std::domain_error(
          "Continuous-time HMM observations must be ordered by non-decreasing TIME within sequence.");
      }
      const MatrixT<Scalar> transition = matrix_exp_pade(
        MatrixT<Scalar>(generator * Scalar(dt)));
      for (int from = 0; from < states; ++from) {
        result.transition[static_cast<std::size_t>(from)].resize(state_count);
        for (int to = 0; to < states; ++to) {
          result.transition[static_cast<std::size_t>(from)]
                           [static_cast<std::size_t>(to)] = transition(from, to);
        }
      }
    } else {
      for (int from = 0; from < states; ++from) {
        const auto begin = output.begin() + static_cast<std::ptrdiff_t>(
          transition_offset + static_cast<std::size_t>(from * states));
        result.transition[static_cast<std::size_t>(from)] =
          hmm_normalize_weights_t(
            std::vector<Scalar>(begin, begin + states),
            engine.hmm_transition_scale, "transition", row);
      }
    }
  }
  const std::size_t emission_offset = state_count + transition_count;
  result.emission.resize(state_count);
  result.log_emission.resize(state_count);
  for (std::size_t state = 0; state < state_count; ++state) {
    const Scalar value = output[emission_offset + state];
    if (engine.hmm_emission_scale == "log") {
      result.log_emission[state] = value;
      result.emission[state] = Scalar(0.0);
    } else {
      if (scalar_value(value) < 0.0) {
        throw std::domain_error(
          "HMM emission likelihoods must be non-negative at row " +
          std::to_string(row + 1) + ".");
      }
      result.emission[state] = value;
      result.log_emission[state] = libertad::scalar_log(
        scalar_floor_t(value, 1e-300));
    }
  }
  return result;
}

template <class Scalar>
Scalar hmm_row_nll_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data, int row, int subject,
    const std::vector<Scalar>& theta, const std::vector<Scalar>& eta,
    int eta_columns, const std::vector<Scalar>& sigma, int mixture_number,
    const Scalar& prediction, bool first, double previous_dv,
    double previous_time, std::vector<Scalar>& filtered,
    HmmRowComponents<Scalar>* components_output = nullptr) {
  const int states = engine.hmm_states;
  const std::size_t state_count = static_cast<std::size_t>(states);
  HmmRowComponents<Scalar> components = hmm_row_components_t(
    engine, data, row, subject, theta, eta, eta_columns, sigma,
    mixture_number, prediction, first, previous_dv, previous_time);
  std::vector<Scalar> prior(state_count);
  if (first) {
    prior = components.initial;
  } else {
    if (filtered.size() != state_count) {
      throw std::logic_error("Hidden Markov filter state has the wrong dimension.");
    }
    std::fill(prior.begin(), prior.end(), Scalar(0.0));
    for (int from = 0; from < states; ++from) {
      for (int to = 0; to < states; ++to) {
        prior[static_cast<std::size_t>(to)] +=
          filtered[static_cast<std::size_t>(from)] *
          components.transition[static_cast<std::size_t>(from)]
                               [static_cast<std::size_t>(to)];
      }
    }
  }
  if (components_output != nullptr) *components_output = components;
  filtered.assign(state_count, Scalar(0.0));
  if (engine.hmm_emission_scale == "log") {
    std::vector<Scalar> log_weight(state_count);
    for (std::size_t state = 0; state < state_count; ++state) {
      log_weight[state] = libertad::scalar_log(
        scalar_floor_t(prior[state], 1e-300)) + components.log_emission[state];
    }
    const Scalar log_likelihood = log_sum_exp_t(log_weight);
    for (std::size_t state = 0; state < state_count; ++state) {
      filtered[state] = libertad::scalar_exp(log_weight[state] - log_likelihood);
    }
    return Scalar(-2.0) * log_likelihood;
  }
  Scalar likelihood = Scalar(0.0);
  for (std::size_t state = 0; state < state_count; ++state) {
    filtered[state] = prior[state] * components.emission[state];
    likelihood += filtered[state];
  }
  if (!(scalar_value(likelihood) > 0.0)) {
    throw std::domain_error("HMM observation likelihood is zero at row " +
                            std::to_string(row + 1) + ".");
  }
  for (Scalar& value : filtered) value /= likelihood;
  return Scalar(-2.0) * libertad::scalar_log(
    scalar_floor_t(likelihood, 1e-300));
}

template <class Scalar>
Scalar positive_definite_gaussian_nll_t(
    const MatrixT<Scalar>& covariance, const VectorT<Scalar>& residual,
    const std::string& context);

inline int residual_group_endpoint(const ResidualGroupSpec& group, double dvid) {
  for (std::size_t index = 0; index < group.dvid.size(); ++index) {
    if (group.dvid[index] == dvid) return static_cast<int>(index);
  }
  return -1;
}

inline int residual_group_for_dvid(const ModelEngine& engine, double dvid) {
  for (std::size_t index = 0; index < engine.residual_groups.size(); ++index) {
    if (residual_group_endpoint(engine.residual_groups[index], dvid) >= 0) {
      return static_cast<int>(index);
    }
  }
  return -1;
}

template <class Scalar>
Scalar residual_group_correlation_t(
    const ResidualGroupSpec& group, int row, int column,
    const std::vector<Scalar>& theta, const std::vector<Scalar>& sigma) {
  const int dimension = static_cast<int>(group.dvid.size());
  if (row < 0 || column < 0 || row >= dimension || column >= dimension) {
    throw std::out_of_range("Residual-group endpoint index is outside its correlation matrix.");
  }
  const std::size_t position = static_cast<std::size_t>(row * dimension + column);
  const std::string& source = group.source[position];
  if (source == "fixed") return Scalar(group.value[position]);
  const std::vector<Scalar>& parameters = source == "theta" ? theta : sigma;
  const int index = group.index[position];
  if (index < 0 || index >= static_cast<int>(parameters.size())) {
    throw std::out_of_range("Residual correlation parameter index is outside its parameter vector.");
  }
  Scalar correlation = parameters[static_cast<std::size_t>(index)];
  if (group.transform == "tanh") correlation = libertad::scalar_tanh(correlation);
  if (!(std::abs(scalar_value(correlation)) < 1.0)) {
    throw std::domain_error("Cross-endpoint residual correlations must be strictly between -1 and 1.");
  }
  return correlation;
}

template <class Scalar>
std::vector<Scalar> residual_grouped_subject_nll_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data,
    const std::vector<Scalar>& prediction, const std::vector<Scalar>& theta,
    const std::vector<Scalar>& sigma,
    const std::vector<Scalar>* variance_prediction = nullptr) {
  Rcpp::NumericVector dv = data["DV"];
  Rcpp::IntegerVector evid = data["EVID"];
  Rcpp::IntegerVector mdv = data["MDV"];
  Rcpp::IntegerVector subjects = data[".ID_INDEX"];
  const bool has_dvid = data.containsElementNamed("DVID");
  int n_subjects = 0;
  for (int value : subjects) n_subjects = std::max(n_subjects, value);
  std::vector<Scalar> result(static_cast<std::size_t>(n_subjects), Scalar(0.0));
  std::vector<bool> processed(static_cast<std::size_t>(data.nrows()), false);
  auto observed = [&](int row) {
    return evid[row] == 0 && mdv[row] == 0 && std::isfinite(dv[row]);
  };
  for (int row = 0; row < data.nrows(); ++row) {
    if (!observed(row) || processed[static_cast<std::size_t>(row)]) continue;
    const int subject = subjects[row] - 1;
    const double dvid = has_dvid ? row_optional(data, "DVID", row, 1.0) : 1.0;
    const int group_index = residual_group_for_dvid(engine, dvid);
    std::vector<int> rows{row};
    if (group_index >= 0) {
      rows.clear();
      const double time = row_optional(data, "TIME", row, 0.0);
      std::vector<bool> endpoint_seen(
        engine.residual_groups[static_cast<std::size_t>(group_index)].dvid.size(), false);
      for (int candidate = 0; candidate < data.nrows(); ++candidate) {
        if (!observed(candidate) || subjects[candidate] - 1 != subject ||
            row_optional(data, "TIME", candidate, 0.0) != time) continue;
        const double candidate_dvid = has_dvid ?
          row_optional(data, "DVID", candidate, 1.0) : 1.0;
        const int endpoint = residual_group_endpoint(
          engine.residual_groups[static_cast<std::size_t>(group_index)], candidate_dvid);
        if (endpoint < 0) continue;
        if (endpoint_seen[static_cast<std::size_t>(endpoint)]) {
          throw std::domain_error(
            "A correlated residual group has duplicate DVID observations at one subject/time.");
        }
        endpoint_seen[static_cast<std::size_t>(endpoint)] = true;
        rows.push_back(candidate);
      }
    }
    const Eigen::Index dimension = static_cast<Eigen::Index>(rows.size());
    VectorT<Scalar> residual(dimension);
    VectorT<Scalar> variance(dimension);
    std::vector<int> endpoint(rows.size(), -1);
    for (Eigen::Index position = 0; position < dimension; ++position) {
      const int selected_row = rows[static_cast<std::size_t>(position)];
      const Scalar f = prediction[static_cast<std::size_t>(selected_row)];
      const Scalar scale_prediction = variance_prediction == nullptr ? f :
        variance_prediction->at(static_cast<std::size_t>(selected_row));
      const double selected_dvid = has_dvid ?
        row_optional(data, "DVID", selected_row, 1.0) : 1.0;
      variance[position] = residual_variance_t(
        engine, scale_prediction, sigma, std::max(1, static_cast<int>(selected_dvid)));
      if (engine.error_type == "exponential") {
        if (!(dv[selected_row] > 0.0)) {
          throw std::domain_error("Exponential residual likelihood requires positive DV.");
        }
        residual[position] = Scalar(std::log(dv[selected_row])) -
          libertad::scalar_log(scalar_floor_t(f, 1e-300));
      } else {
        residual[position] = Scalar(dv[selected_row]) - f;
      }
      if (group_index >= 0) {
        endpoint[static_cast<std::size_t>(position)] = residual_group_endpoint(
          engine.residual_groups[static_cast<std::size_t>(group_index)], selected_dvid);
      }
      processed[static_cast<std::size_t>(selected_row)] = true;
    }
    MatrixT<Scalar> covariance(dimension, dimension);
    for (Eigen::Index first = 0; first < dimension; ++first) {
      for (Eigen::Index second = 0; second < dimension; ++second) {
        Scalar correlation = first == second ? Scalar(1.0) : Scalar(0.0);
        if (group_index >= 0 && first != second) {
          correlation = residual_group_correlation_t(
            engine.residual_groups[static_cast<std::size_t>(group_index)],
            endpoint[static_cast<std::size_t>(first)],
            endpoint[static_cast<std::size_t>(second)], theta, sigma);
        }
        covariance(first, second) = correlation *
          libertad::scalar_sqrt(variance[first] * variance[second]);
      }
    }
    result[static_cast<std::size_t>(subject)] += positive_definite_gaussian_nll_t(
      covariance, residual, "Cross-endpoint residual covariance");
  }
  return result;
}

template <class Scalar>
struct KalmanFilterState {
  VectorT<Scalar> mean;
  MatrixT<Scalar> covariance;
};

template <class Scalar>
struct KalmanRowComponents {
  VectorT<Scalar> predicted_mean;
  MatrixT<Scalar> predicted_covariance;
  VectorT<Scalar> filtered_mean;
  MatrixT<Scalar> filtered_covariance;
  MatrixT<Scalar> transition;
  MatrixT<Scalar> smoother_cross_covariance;
  VectorT<Scalar> observation;
  Scalar observation_variance = Scalar(0.0);
  Scalar innovation = Scalar(0.0);
  Scalar innovation_variance = Scalar(0.0);
  std::vector<VectorT<Scalar>> particle_states;
  std::vector<Scalar> particle_weights;
  std::vector<int> particle_ancestors;
  std::vector<int> particle_regimes;
};

template <class Scalar>
MatrixT<Scalar> state_space_cholesky_t(const MatrixT<Scalar>& covariance,
                                       const char* context) {
  const Eigen::Index dimension = covariance.rows();
  MatrixT<Scalar> lower = MatrixT<Scalar>::Zero(dimension, dimension);
  for (Eigen::Index row = 0; row < dimension; ++row) {
    for (Eigen::Index column = 0; column <= row; ++column) {
      Scalar value = Scalar(0.5) *
        (covariance(row, column) + covariance(column, row));
      for (Eigen::Index k = 0; k < column; ++k) {
        value -= lower(row, k) * lower(column, k);
      }
      if (row == column) {
        if (!(scalar_value(value) > -1e-10)) {
          throw std::domain_error(std::string(context) + " is not positive semidefinite.");
        }
        lower(row, column) = libertad::scalar_sqrt(scalar_floor_t(value, 1e-12));
      } else {
        lower(row, column) = value / lower(column, column);
      }
    }
  }
  return lower;
}

inline std::uint64_t state_space_hash(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31U);
}

inline double state_space_uniform(int seed, int row, int particle,
                                  int dimension, int stream = 0) {
  std::uint64_t value = static_cast<std::uint64_t>(static_cast<std::uint32_t>(seed));
  value ^= static_cast<std::uint64_t>(row + 1) * 0x9e3779b97f4a7c15ULL;
  value ^= static_cast<std::uint64_t>(particle + 1) * 0xbf58476d1ce4e5b9ULL;
  value ^= static_cast<std::uint64_t>(dimension + 1) * 0x94d049bb133111ebULL;
  value ^= static_cast<std::uint64_t>(stream + 1) * 0xd6e8feb86659fd93ULL;
  const std::uint64_t hashed = state_space_hash(value);
  return (static_cast<double>((hashed >> 11U) + 1ULL)) /
    9007199254740993.0;
}

inline double state_space_normal(int seed, int row, int particle,
                                 int dimension, int stream = 0) {
  const double first = state_space_uniform(seed, row, particle, 2 * dimension, stream);
  const double second = state_space_uniform(seed, row, particle, 2 * dimension + 1, stream);
  return std::sqrt(-2.0 * std::log(std::max(first, 1e-15))) *
    std::cos(2.0 * 3.14159265358979323846 * second);
}

template <class Scalar>
struct NonlinearStateOutputs {
  VectorT<Scalar> initial_mean;
  MatrixT<Scalar> initial_covariance;
  VectorT<Scalar> transition;
  MatrixT<Scalar> process_covariance;
  Scalar observation = Scalar(0.0);
  Scalar observation_variance = Scalar(0.0);
};

template <class Scalar>
struct ParticleFilterState {
  std::vector<VectorT<Scalar>> particles;
  std::vector<Scalar> weights;
  std::vector<int> regimes;
};

template <class Scalar>
NonlinearStateOutputs<Scalar> nonlinear_state_raw_outputs_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data, int row, int subject,
    const std::vector<Scalar>& theta, const std::vector<Scalar>& eta,
    int eta_columns, const std::vector<Scalar>& sigma, int mixture_number,
    const Scalar& prediction, bool first, double previous_dv,
    double previous_time, const VectorT<Scalar>& state) {
  ParametersT<Scalar> overrides;
  for (int index = 0; index < engine.kalman_states; ++index) {
    overrides[engine.kalman_state_inputs[static_cast<std::size_t>(index)]] = state[index];
  }
  const std::vector<Scalar> output = evaluate_error_outputs_t(
    engine, data, row, subject, theta, eta, eta_columns, sigma,
    mixture_number, prediction, first, previous_dv, previous_time,
    engine.kalman_outputs, "Nonlinear state-space likelihood", &overrides);
  const int states = engine.kalman_states;
  const std::size_t state_count = static_cast<std::size_t>(states);
  const std::size_t expected = 2U * state_count * state_count +
    2U * state_count + 2U;
  if (output.size() != expected) {
    throw std::logic_error("Compiled nonlinear state-space outputs are inconsistent.");
  }
  std::size_t cursor = 0U;
  NonlinearStateOutputs<Scalar> result;
  result.initial_mean.resize(states);
  result.transition.resize(states);
  for (int index = 0; index < states; ++index) result.initial_mean[index] = output[cursor++];
  auto read_matrix = [&]() {
    MatrixT<Scalar> matrix(states, states);
    for (int matrix_row = 0; matrix_row < states; ++matrix_row) {
      for (int matrix_column = 0; matrix_column < states; ++matrix_column) {
        matrix(matrix_row, matrix_column) = output[cursor++];
      }
    }
    return matrix;
  };
  result.initial_covariance = read_matrix();
  for (int index = 0; index < states; ++index) result.transition[index] = output[cursor++];
  result.process_covariance = read_matrix();
  result.observation = output[cursor++];
  result.observation_variance = output[cursor++];
  return result;
}

template <class Scalar>
NonlinearStateOutputs<Scalar> nonlinear_state_outputs_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data, int row, int subject,
    const std::vector<Scalar>& theta, const std::vector<Scalar>& eta,
    int eta_columns, const std::vector<Scalar>& sigma, int mixture_number,
    const Scalar& prediction, bool first, double previous_dv,
    double previous_time, const VectorT<Scalar>& state) {
  NonlinearStateOutputs<Scalar> result = nonlinear_state_raw_outputs_t(
    engine, data, row, subject, theta, eta, eta_columns, sigma,
    mixture_number, prediction, first, previous_dv, previous_time, state);
  if (engine.kalman_dynamics != "sde" || first) return result;
  const double current_time = row_optional(data, "TIME", row, previous_time);
  const double interval = std::max(0.0, current_time - previous_time);
  const double step = interval / static_cast<double>(engine.kalman_sde_substeps);
  VectorT<Scalar> current = state;
  MatrixT<Scalar> covariance = MatrixT<Scalar>::Zero(
    engine.kalman_states, engine.kalman_states);
  if (step > 0.0) {
    for (int substep = 0; substep < engine.kalman_sde_substeps; ++substep) {
      const NonlinearStateOutputs<Scalar> local = nonlinear_state_raw_outputs_t(
        engine, data, row, subject, theta, eta, eta_columns, sigma,
        mixture_number, prediction, false, previous_dv,
        previous_time + substep * step, current);
      current += Scalar(step) * local.transition;
      covariance += Scalar(step) * local.process_covariance *
        local.process_covariance.transpose();
    }
  }
  result.transition = current;
  result.process_covariance = Scalar(0.5) *
    (covariance + covariance.transpose());
  return result;
}

template <class Scalar>
struct SwitchingStateOutputs {
  std::vector<Scalar> initial;
  MatrixT<Scalar> regime_transition;
  std::vector<NonlinearStateOutputs<Scalar>> regime;
};

template <class Scalar>
void normalize_switching_weights(std::vector<Scalar>& value,
                                 const std::string& scale,
                                 const char* context) {
  Scalar total = Scalar(0.0);
  if (scale == "log") {
    double maximum = -std::numeric_limits<double>::infinity();
    for (const Scalar& item : value) maximum = std::max(maximum, scalar_value(item));
    for (Scalar& item : value) {
      item = libertad::scalar_exp(item - Scalar(maximum));
      total += item;
    }
  } else {
    for (const Scalar& item : value) {
      if (scalar_value(item) < 0.0) {
        throw std::domain_error(std::string(context) + " contains a negative weight.");
      }
      total += item;
    }
  }
  if (!(scalar_value(total) > 0.0)) {
    throw std::domain_error(std::string(context) + " has zero total weight.");
  }
  for (Scalar& item : value) item /= total;
}

template <class Scalar>
SwitchingStateOutputs<Scalar> switching_state_raw_outputs_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data, int row, int subject,
    const std::vector<Scalar>& theta, const std::vector<Scalar>& eta,
    int eta_columns, const std::vector<Scalar>& sigma, int mixture_number,
    const Scalar& prediction, bool first, double previous_dv,
    double previous_time, const VectorT<Scalar>& state) {
  ParametersT<Scalar> overrides;
  for (int index = 0; index < engine.kalman_states; ++index) {
    overrides[engine.kalman_state_inputs[static_cast<std::size_t>(index)]] = state[index];
  }
  const std::vector<Scalar> output = evaluate_error_outputs_t(
    engine, data, row, subject, theta, eta, eta_columns, sigma,
    mixture_number, prediction, first, previous_dv, previous_time,
    engine.switching_outputs, "Switching state-space likelihood", &overrides);
  const int regimes = engine.switching_regimes;
  const int states = engine.kalman_states;
  const std::size_t expected = static_cast<std::size_t>(regimes + regimes * regimes +
    regimes * (states + states * states + 2));
  if (output.size() != expected) {
    throw std::logic_error("Compiled switching state-space outputs are inconsistent.");
  }
  std::size_t cursor = 0U;
  SwitchingStateOutputs<Scalar> result;
  result.initial.resize(static_cast<std::size_t>(regimes));
  for (int regime = 0; regime < regimes; ++regime) result.initial[regime] = output[cursor++];
  normalize_switching_weights(
    result.initial, engine.switching_initial_scale, "Initial regime distribution");
  result.regime_transition.resize(regimes, regimes);
  for (int from = 0; from < regimes; ++from) {
    std::vector<Scalar> local(static_cast<std::size_t>(regimes));
    for (int to = 0; to < regimes; ++to) local[to] = output[cursor++];
    normalize_switching_weights(
      local, engine.switching_transition_scale, "Regime transition row");
    for (int to = 0; to < regimes; ++to) result.regime_transition(from, to) = local[to];
  }
  result.regime.resize(static_cast<std::size_t>(regimes));
  for (int regime = 0; regime < regimes; ++regime) {
    NonlinearStateOutputs<Scalar>& local = result.regime[static_cast<std::size_t>(regime)];
    local.transition.resize(states);
    for (int index = 0; index < states; ++index) local.transition[index] = output[cursor++];
    local.process_covariance.resize(states, states);
    for (int matrix_row = 0; matrix_row < states; ++matrix_row) {
      for (int matrix_column = 0; matrix_column < states; ++matrix_column) {
        local.process_covariance(matrix_row, matrix_column) = output[cursor++];
      }
    }
    local.observation = output[cursor++];
    local.observation_variance = output[cursor++];
  }
  return result;
}

template <class Scalar>
int sample_switching_regime(const std::vector<Scalar>& probability,
                            double uniform) {
  Scalar cumulative = probability.front();
  int selected = 0;
  while (selected + 1 < static_cast<int>(probability.size()) &&
         path_lt(cumulative, Scalar(uniform))) {
    ++selected;
    cumulative += probability[static_cast<std::size_t>(selected)];
  }
  return selected;
}

template <class Scalar>
Scalar particle_row_nll_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data, int row, int subject,
    const std::vector<Scalar>& theta, const std::vector<Scalar>& eta,
    int eta_columns, const std::vector<Scalar>& sigma, int mixture_number,
    const Scalar& prediction, bool first, double previous_dv,
    double previous_time, ParticleFilterState<Scalar>& filter,
    KalmanRowComponents<Scalar>* components_output = nullptr) {
  const int states = engine.kalman_states;
  const int particles = engine.kalman_particles;
  const VectorT<Scalar> zero = VectorT<Scalar>::Zero(states);
  const NonlinearStateOutputs<Scalar> base = nonlinear_state_outputs_t(
    engine, data, row, subject, theta, eta, eta_columns, sigma,
    mixture_number, prediction, first, previous_dv, previous_time,
    first ? zero : filter.particles.front());
  MatrixT<Scalar> root = MatrixT<Scalar>::Zero(states, states);
  if (first || !engine.switching_enabled) {
    const MatrixT<Scalar> covariance = first ?
      base.initial_covariance : base.process_covariance;
    root = state_space_cholesky_t(
      covariance, first ? "Particle initial covariance" : "Particle process covariance");
  }
  std::vector<VectorT<Scalar>> propagated(static_cast<std::size_t>(particles));
  std::vector<int> propagated_regime(static_cast<std::size_t>(particles), 0);
  std::vector<Scalar> regime_importance(
    static_cast<std::size_t>(particles), Scalar(1.0));
  std::vector<Scalar> prior_weight(static_cast<std::size_t>(particles),
                                   Scalar(1.0 / static_cast<double>(particles)));
  if (!first) {
    if (filter.particles.size() != static_cast<std::size_t>(particles) ||
        filter.weights.size() != static_cast<std::size_t>(particles) ||
        (engine.switching_enabled &&
         filter.regimes.size() != static_cast<std::size_t>(particles))) {
      throw std::logic_error("Particle filter state dimension changed within a sequence.");
    }
    prior_weight = filter.weights;
  }
  for (int particle = 0; particle < particles; ++particle) {
    VectorT<Scalar> normal(states);
    for (int state = 0; state < states; ++state) {
      normal[state] = Scalar(state_space_normal(
        engine.kalman_seed + 7919 * subject, row, particle, state, first ? 1 : 2));
    }
    if (engine.switching_enabled) {
      const VectorT<Scalar>& current = first ? zero :
        filter.particles[static_cast<std::size_t>(particle)];
      const SwitchingStateOutputs<Scalar> switching = switching_state_raw_outputs_t(
        engine, data, row, subject, theta, eta, eta_columns, sigma,
        mixture_number, prediction, first, previous_dv, previous_time, current);
      std::vector<Scalar> regime_probability;
      if (first) {
        regime_probability = switching.initial;
      } else {
        regime_probability.resize(static_cast<std::size_t>(engine.switching_regimes));
        const int previous_regime = filter.regimes[static_cast<std::size_t>(particle)];
        for (int regime = 0; regime < engine.switching_regimes; ++regime) {
          regime_probability[static_cast<std::size_t>(regime)] =
            switching.regime_transition(previous_regime, regime);
        }
      }
      const int target_regime = particle % engine.switching_regimes;
      propagated_regime[static_cast<std::size_t>(particle)] = target_regime;
      const int target_count = particles / engine.switching_regimes +
        (target_regime < particles % engine.switching_regimes ? 1 : 0);
      const double proposal_probability =
        static_cast<double>(target_count) / static_cast<double>(particles);
      regime_importance[static_cast<std::size_t>(particle)] =
        regime_probability[static_cast<std::size_t>(target_regime)] /
        Scalar(proposal_probability);
    }
    if (first) {
      propagated[static_cast<std::size_t>(particle)] = base.initial_mean + root * normal;
    } else if (engine.kalman_dynamics == "sde") {
      VectorT<Scalar> current = filter.particles[static_cast<std::size_t>(particle)];
      const double current_time = row_optional(data, "TIME", row, previous_time);
      const double interval = std::max(0.0, current_time - previous_time);
      const double step = interval / static_cast<double>(engine.kalman_sde_substeps);
      for (int substep = 0; substep < engine.kalman_sde_substeps; ++substep) {
        const NonlinearStateOutputs<Scalar> local = engine.switching_enabled ?
          switching_state_raw_outputs_t(
            engine, data, row, subject, theta, eta, eta_columns, sigma,
            mixture_number, prediction, false, previous_dv,
            previous_time + substep * step, current).regime[
              static_cast<std::size_t>(propagated_regime[static_cast<std::size_t>(particle)])] :
          nonlinear_state_raw_outputs_t(
            engine, data, row, subject, theta, eta, eta_columns, sigma,
            mixture_number, prediction, false, previous_dv,
            previous_time + substep * step, current);
        VectorT<Scalar> increment(states);
        for (int state = 0; state < states; ++state) {
          increment[state] = Scalar(std::sqrt(step) * state_space_normal(
            engine.kalman_seed + 7919 * subject, row, particle,
            substep * states + state, 6));
        }
        VectorT<Scalar> next = current + Scalar(step) * local.transition +
          local.process_covariance * increment;
        if (engine.kalman_sde_method == "milstein") {
          for (int first_state = 0; first_state < states; ++first_state) {
            for (int second_state = 0; second_state < states; ++second_state) {
              if (first_state != second_state &&
                  std::abs(scalar_value(local.process_covariance(first_state, second_state))) > 1e-12) {
                throw std::domain_error(
                  "Milstein propagation currently requires diagonal diffusion.");
              }
            }
            const double derivative_step = engine.kalman_jacobian_step *
              std::max(1.0, std::abs(scalar_value(current[first_state])));
            VectorT<Scalar> plus = current;
            VectorT<Scalar> minus = current;
            plus[first_state] += Scalar(derivative_step);
            minus[first_state] -= Scalar(derivative_step);
            const int regime = propagated_regime[static_cast<std::size_t>(particle)];
            const Scalar upper = engine.switching_enabled ?
              switching_state_raw_outputs_t(
                engine, data, row, subject, theta, eta, eta_columns, sigma,
                mixture_number, prediction, false, previous_dv,
                previous_time + substep * step, plus).regime[
                  static_cast<std::size_t>(regime)].process_covariance(first_state, first_state) :
              nonlinear_state_raw_outputs_t(
                engine, data, row, subject, theta, eta, eta_columns, sigma,
                mixture_number, prediction, false, previous_dv,
                previous_time + substep * step, plus).process_covariance(first_state, first_state);
            const Scalar lower = engine.switching_enabled ?
              switching_state_raw_outputs_t(
                engine, data, row, subject, theta, eta, eta_columns, sigma,
                mixture_number, prediction, false, previous_dv,
                previous_time + substep * step, minus).regime[
                  static_cast<std::size_t>(regime)].process_covariance(first_state, first_state) :
              nonlinear_state_raw_outputs_t(
                engine, data, row, subject, theta, eta, eta_columns, sigma,
                mixture_number, prediction, false, previous_dv,
                previous_time + substep * step, minus).process_covariance(first_state, first_state);
            const Scalar derivative = (upper - lower) / Scalar(2.0 * derivative_step);
            const Scalar diffusion = local.process_covariance(first_state, first_state);
            next[first_state] += Scalar(0.5) * diffusion * derivative *
              (increment[first_state] * increment[first_state] - Scalar(step));
          }
        }
        current = next;
      }
      propagated[static_cast<std::size_t>(particle)] = current;
    } else {
      const VectorT<Scalar>& current = filter.particles[static_cast<std::size_t>(particle)];
      const NonlinearStateOutputs<Scalar> point = engine.switching_enabled ?
        switching_state_raw_outputs_t(
          engine, data, row, subject, theta, eta, eta_columns, sigma,
          mixture_number, prediction, first, previous_dv, previous_time, current).regime[
            static_cast<std::size_t>(propagated_regime[static_cast<std::size_t>(particle)])] :
        nonlinear_state_outputs_t(
          engine, data, row, subject, theta, eta, eta_columns, sigma,
          mixture_number, prediction, first, previous_dv, previous_time, current);
      const MatrixT<Scalar> local_root = engine.switching_enabled ?
        state_space_cholesky_t(point.process_covariance,
                               "Switching particle process covariance") : root;
      propagated[static_cast<std::size_t>(particle)] =
        point.transition + local_root * normal;
    }
  }
  KalmanRowComponents<Scalar> components;
  std::vector<Scalar> predictive_weight = prior_weight;
  if (engine.switching_enabled) {
    Scalar predictive_total = Scalar(0.0);
    for (int particle = 0; particle < particles; ++particle) {
      predictive_weight[static_cast<std::size_t>(particle)] *=
        regime_importance[static_cast<std::size_t>(particle)];
      predictive_total += predictive_weight[static_cast<std::size_t>(particle)];
    }
    if (!(scalar_value(predictive_total) > 0.0)) {
      throw std::domain_error("Switching-state predictive regime mass is zero.");
    }
    for (Scalar& weight : predictive_weight) weight /= predictive_total;
  }
  components.predicted_mean = VectorT<Scalar>::Zero(states);
  for (int particle = 0; particle < particles; ++particle) {
    components.predicted_mean += predictive_weight[static_cast<std::size_t>(particle)] *
      propagated[static_cast<std::size_t>(particle)];
  }
  components.predicted_covariance = MatrixT<Scalar>::Zero(states, states);
  for (int particle = 0; particle < particles; ++particle) {
    const VectorT<Scalar> centered = propagated[static_cast<std::size_t>(particle)] -
      components.predicted_mean;
    components.predicted_covariance += predictive_weight[static_cast<std::size_t>(particle)] *
      centered * centered.transpose();
  }
  std::vector<Scalar> log_weight(static_cast<std::size_t>(particles));
  std::vector<Scalar> observation_value(static_cast<std::size_t>(particles));
  double maximum = -std::numeric_limits<double>::infinity();
  const Scalar baseline = engine.kalman_prediction_baseline ? prediction : Scalar(0.0);
  const Scalar observed = Scalar(row_optional(data, "DV", row, NA_REAL));
  Scalar predicted_observation = Scalar(0.0);
  for (int particle = 0; particle < particles; ++particle) {
    const NonlinearStateOutputs<Scalar> point = engine.switching_enabled ?
      switching_state_raw_outputs_t(
        engine, data, row, subject, theta, eta, eta_columns, sigma,
        mixture_number, prediction, first, previous_dv, previous_time,
        propagated[static_cast<std::size_t>(particle)]).regime[
          static_cast<std::size_t>(propagated_regime[static_cast<std::size_t>(particle)])] :
      nonlinear_state_outputs_t(
        engine, data, row, subject, theta, eta, eta_columns, sigma,
        mixture_number, prediction, first, previous_dv, previous_time,
        propagated[static_cast<std::size_t>(particle)]);
    if (!(scalar_value(point.observation_variance) > 1e-14)) {
      throw std::domain_error("Particle observation variance must be positive.");
    }
    observation_value[static_cast<std::size_t>(particle)] = baseline + point.observation;
    predicted_observation += predictive_weight[static_cast<std::size_t>(particle)] *
      observation_value[static_cast<std::size_t>(particle)];
    const Scalar residual = observed - observation_value[static_cast<std::size_t>(particle)];
    log_weight[static_cast<std::size_t>(particle)] =
      libertad::scalar_log(scalar_floor_t(
        prior_weight[static_cast<std::size_t>(particle)] *
          regime_importance[static_cast<std::size_t>(particle)], 1e-300)) -
      Scalar(0.5) * (libertad::scalar_log(point.observation_variance) +
                     residual * residual / point.observation_variance);
    maximum = std::max(maximum, scalar_value(log_weight[static_cast<std::size_t>(particle)]));
  }
  Scalar normalizer = Scalar(0.0);
  for (const Scalar& value : log_weight) normalizer += libertad::scalar_exp(value - Scalar(maximum));
  const Scalar log_likelihood = Scalar(maximum) + libertad::scalar_log(normalizer);
  filter.weights.resize(static_cast<std::size_t>(particles));
  for (int particle = 0; particle < particles; ++particle) {
    filter.weights[static_cast<std::size_t>(particle)] =
      libertad::scalar_exp(log_weight[static_cast<std::size_t>(particle)] - log_likelihood);
  }
  components.innovation = observed - predicted_observation;
  components.innovation_variance = Scalar(0.0);
  for (int particle = 0; particle < particles; ++particle) {
    const Scalar centered = observation_value[static_cast<std::size_t>(particle)] -
      predicted_observation;
    components.innovation_variance += predictive_weight[static_cast<std::size_t>(particle)] *
      centered * centered;
  }
  components.transition = MatrixT<Scalar>::Identity(states, states);
  components.smoother_cross_covariance = MatrixT<Scalar>::Zero(states, states);
  VectorT<Scalar> filtered_mean = VectorT<Scalar>::Zero(states);
  for (int particle = 0; particle < particles; ++particle) {
    filtered_mean += filter.weights[static_cast<std::size_t>(particle)] *
      propagated[static_cast<std::size_t>(particle)];
  }
  MatrixT<Scalar> filtered_covariance = MatrixT<Scalar>::Zero(states, states);
  for (int particle = 0; particle < particles; ++particle) {
    const VectorT<Scalar> centered = propagated[static_cast<std::size_t>(particle)] - filtered_mean;
    filtered_covariance += filter.weights[static_cast<std::size_t>(particle)] *
      centered * centered.transpose();
  }
  Scalar squared_weight = Scalar(0.0);
  for (const Scalar& weight : filter.weights) squared_weight += weight * weight;
  std::vector<int> ancestors(static_cast<std::size_t>(particles));
  std::iota(ancestors.begin(), ancestors.end(), 0);
  // Keep both the resampling decision and systematic ancestry comparisons on
  // the CppAD comparison tape. A changed ancestry then raises TapePathChange
  // on replay and the population engine records the correct particle graph at
  // the new parameter point instead of silently reusing stale indices.
  const Scalar resampling_boundary = Scalar(
    1.0 / (engine.kalman_ess_threshold * static_cast<double>(particles)));
  if (squared_weight > resampling_boundary) {
    const double offset = state_space_uniform(
      engine.kalman_seed + 104729 * subject, row, 0, 0, 9) /
      static_cast<double>(particles);
    std::vector<VectorT<Scalar>> resampled(static_cast<std::size_t>(particles));
    std::vector<int> resampled_regime(static_cast<std::size_t>(particles));
    int selected = 0;
    Scalar cumulative = filter.weights[0];
    for (int particle = 0; particle < particles; ++particle) {
      const double target = offset + static_cast<double>(particle) /
        static_cast<double>(particles);
      while (selected + 1 < particles && cumulative < Scalar(target)) {
        ++selected;
        cumulative += filter.weights[static_cast<std::size_t>(selected)];
      }
      resampled[static_cast<std::size_t>(particle)] =
        propagated[static_cast<std::size_t>(selected)];
      ancestors[static_cast<std::size_t>(particle)] = selected;
      resampled_regime[static_cast<std::size_t>(particle)] =
        propagated_regime[static_cast<std::size_t>(selected)];
    }
    filter.particles = std::move(resampled);
    if (engine.switching_enabled) filter.regimes = std::move(resampled_regime);
    std::fill(filter.weights.begin(), filter.weights.end(),
              Scalar(1.0 / static_cast<double>(particles)));
  } else {
    filter.particles = std::move(propagated);
    if (engine.switching_enabled) filter.regimes = std::move(propagated_regime);
  }
  if (components_output != nullptr) {
    components.observation = VectorT<Scalar>::Zero(states);
    components.filtered_mean = filtered_mean;
    components.filtered_covariance = filtered_covariance;
    components.particle_states = filter.particles;
    components.particle_weights = filter.weights;
    components.particle_ancestors = ancestors;
    components.particle_regimes = filter.regimes;
    *components_output = components;
  }
  return Scalar(-2.0) * log_likelihood;
}

template <class Scalar>
Scalar nonlinear_kalman_row_nll_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data, int row, int subject,
    const std::vector<Scalar>& theta, const std::vector<Scalar>& eta,
    int eta_columns, const std::vector<Scalar>& sigma, int mixture_number,
    const Scalar& prediction, bool first, double previous_dv,
    double previous_time, KalmanFilterState<Scalar>& filtered,
    KalmanRowComponents<Scalar>* components_output) {
  const int states = engine.kalman_states;
  const VectorT<Scalar> zero = VectorT<Scalar>::Zero(states);
  NonlinearStateOutputs<Scalar> base = nonlinear_state_outputs_t(
    engine, data, row, subject, theta, eta, eta_columns, sigma,
    mixture_number, prediction, first, previous_dv, previous_time,
    first ? zero : filtered.mean);
  KalmanRowComponents<Scalar> components;
  components.transition = MatrixT<Scalar>::Identity(states, states);
  components.smoother_cross_covariance = MatrixT<Scalar>::Zero(states, states);
  if (first) {
    components.predicted_mean = base.initial_mean;
    components.predicted_covariance = Scalar(0.5) *
      (base.initial_covariance + base.initial_covariance.transpose());
  } else if (engine.kalman_filter_type == "ekf") {
    const VectorT<Scalar> previous_mean = filtered.mean;
    const MatrixT<Scalar> previous_covariance = filtered.covariance;
    components.predicted_mean = base.transition;
    MatrixT<Scalar> jacobian(states, states);
    for (int column = 0; column < states; ++column) {
      const double step = engine.kalman_jacobian_step *
        std::max(1.0, std::abs(scalar_value(previous_mean[column])));
      VectorT<Scalar> plus = previous_mean;
      VectorT<Scalar> minus = previous_mean;
      plus[column] += Scalar(step);
      minus[column] -= Scalar(step);
      const VectorT<Scalar> upper = nonlinear_state_outputs_t(
        engine, data, row, subject, theta, eta, eta_columns, sigma,
        mixture_number, prediction, first, previous_dv, previous_time, plus).transition;
      const VectorT<Scalar> lower = nonlinear_state_outputs_t(
        engine, data, row, subject, theta, eta, eta_columns, sigma,
        mixture_number, prediction, first, previous_dv, previous_time, minus).transition;
      jacobian.col(column) = (upper - lower) / Scalar(2.0 * step);
    }
    components.transition = jacobian;
    components.smoother_cross_covariance = previous_covariance * jacobian.transpose();
    components.predicted_covariance = jacobian * previous_covariance *
      jacobian.transpose() + base.process_covariance;
  } else if (engine.kalman_filter_type == "ukf") {
    const VectorT<Scalar> previous_mean = filtered.mean;
    const MatrixT<Scalar> previous_covariance = filtered.covariance;
    const double lambda = engine.kalman_ukf_alpha * engine.kalman_ukf_alpha *
      (static_cast<double>(states) + engine.kalman_ukf_kappa) -
      static_cast<double>(states);
    const double scale = static_cast<double>(states) + lambda;
    const int points = 2 * states + 1;
    std::vector<double> mean_weight(static_cast<std::size_t>(points),
                                    1.0 / (2.0 * scale));
    std::vector<double> covariance_weight = mean_weight;
    mean_weight[0] = lambda / scale;
    covariance_weight[0] = mean_weight[0] +
      (1.0 - engine.kalman_ukf_alpha * engine.kalman_ukf_alpha +
       engine.kalman_ukf_beta);
    const MatrixT<Scalar> root = state_space_cholesky_t(
      previous_covariance, "UKF filtered covariance") * Scalar(std::sqrt(scale));
    std::vector<VectorT<Scalar>> sigma_points(static_cast<std::size_t>(points));
    std::vector<VectorT<Scalar>> propagated(static_cast<std::size_t>(points));
    sigma_points[0] = previous_mean;
    for (int column = 0; column < states; ++column) {
      sigma_points[static_cast<std::size_t>(1 + column)] = previous_mean + root.col(column);
      sigma_points[static_cast<std::size_t>(1 + states + column)] = previous_mean - root.col(column);
    }
    components.predicted_mean = VectorT<Scalar>::Zero(states);
    for (int point = 0; point < points; ++point) {
      propagated[static_cast<std::size_t>(point)] = nonlinear_state_outputs_t(
        engine, data, row, subject, theta, eta, eta_columns, sigma,
        mixture_number, prediction, first, previous_dv, previous_time,
        sigma_points[static_cast<std::size_t>(point)]).transition;
      components.predicted_mean += Scalar(mean_weight[static_cast<std::size_t>(point)]) *
        propagated[static_cast<std::size_t>(point)];
    }
    components.predicted_covariance = base.process_covariance;
    for (int point = 0; point < points; ++point) {
      const VectorT<Scalar> before = sigma_points[static_cast<std::size_t>(point)] - previous_mean;
      const VectorT<Scalar> after = propagated[static_cast<std::size_t>(point)] -
        components.predicted_mean;
      components.predicted_covariance +=
        Scalar(covariance_weight[static_cast<std::size_t>(point)]) *
        after * after.transpose();
      components.smoother_cross_covariance +=
        Scalar(covariance_weight[static_cast<std::size_t>(point)]) *
        before * after.transpose();
    }
  } else {
    throw std::logic_error("Particle filtering uses the dedicated particle state path.");
  }
  components.predicted_covariance = Scalar(0.5) *
    (components.predicted_covariance + components.predicted_covariance.transpose());

  VectorT<Scalar> observation_gradient(states);
  Scalar observation_mean = Scalar(0.0);
  if (engine.kalman_filter_type == "ukf") {
    const double lambda = engine.kalman_ukf_alpha * engine.kalman_ukf_alpha *
      (static_cast<double>(states) + engine.kalman_ukf_kappa) -
      static_cast<double>(states);
    const double scale = static_cast<double>(states) + lambda;
    const int points = 2 * states + 1;
    std::vector<double> mean_weight(static_cast<std::size_t>(points),
                                    1.0 / (2.0 * scale));
    std::vector<double> covariance_weight = mean_weight;
    mean_weight[0] = lambda / scale;
    covariance_weight[0] = mean_weight[0] +
      (1.0 - engine.kalman_ukf_alpha * engine.kalman_ukf_alpha +
       engine.kalman_ukf_beta);
    const MatrixT<Scalar> root = state_space_cholesky_t(
      components.predicted_covariance, "UKF predicted covariance") *
      Scalar(std::sqrt(scale));
    std::vector<VectorT<Scalar>> points_state(static_cast<std::size_t>(points));
    std::vector<Scalar> values(static_cast<std::size_t>(points));
    points_state[0] = components.predicted_mean;
    for (int column = 0; column < states; ++column) {
      points_state[static_cast<std::size_t>(1 + column)] =
        components.predicted_mean + root.col(column);
      points_state[static_cast<std::size_t>(1 + states + column)] =
        components.predicted_mean - root.col(column);
    }
    for (int point = 0; point < points; ++point) {
      values[static_cast<std::size_t>(point)] = nonlinear_state_outputs_t(
        engine, data, row, subject, theta, eta, eta_columns, sigma,
        mixture_number, prediction, first, previous_dv, previous_time,
        points_state[static_cast<std::size_t>(point)]).observation;
      observation_mean += Scalar(mean_weight[static_cast<std::size_t>(point)]) *
        values[static_cast<std::size_t>(point)];
    }
    components.innovation_variance = base.observation_variance;
    VectorT<Scalar> cross = VectorT<Scalar>::Zero(states);
    for (int point = 0; point < points; ++point) {
      const Scalar centered = values[static_cast<std::size_t>(point)] - observation_mean;
      components.innovation_variance +=
        Scalar(covariance_weight[static_cast<std::size_t>(point)]) * centered * centered;
      cross += Scalar(covariance_weight[static_cast<std::size_t>(point)]) *
        (points_state[static_cast<std::size_t>(point)] - components.predicted_mean) * centered;
    }
    observation_gradient = cross / components.innovation_variance;
  } else {
    const NonlinearStateOutputs<Scalar> at_mean = nonlinear_state_outputs_t(
      engine, data, row, subject, theta, eta, eta_columns, sigma,
      mixture_number, prediction, first, previous_dv, previous_time,
      components.predicted_mean);
    observation_mean = at_mean.observation;
    for (int column = 0; column < states; ++column) {
      const double step = engine.kalman_jacobian_step *
        std::max(1.0, std::abs(scalar_value(components.predicted_mean[column])));
      VectorT<Scalar> plus = components.predicted_mean;
      VectorT<Scalar> minus = components.predicted_mean;
      plus[column] += Scalar(step);
      minus[column] -= Scalar(step);
      const Scalar upper = nonlinear_state_outputs_t(
        engine, data, row, subject, theta, eta, eta_columns, sigma,
        mixture_number, prediction, first, previous_dv, previous_time, plus).observation;
      const Scalar lower = nonlinear_state_outputs_t(
        engine, data, row, subject, theta, eta, eta_columns, sigma,
        mixture_number, prediction, first, previous_dv, previous_time, minus).observation;
      observation_gradient[column] = (upper - lower) / Scalar(2.0 * step);
    }
    components.innovation_variance =
      (observation_gradient.transpose() * components.predicted_covariance *
       observation_gradient)[0] + at_mean.observation_variance;
    observation_gradient = components.predicted_covariance *
      observation_gradient / components.innovation_variance;
  }
  const Scalar baseline = engine.kalman_prediction_baseline ? prediction : Scalar(0.0);
  components.innovation = Scalar(row_optional(data, "DV", row, NA_REAL)) -
    baseline - observation_mean;
  if (!(scalar_value(components.innovation_variance) > 1e-14)) {
    throw std::domain_error("Nonlinear state-space innovation variance must be positive.");
  }
  const VectorT<Scalar> gain = observation_gradient;
  filtered.mean = components.predicted_mean + gain * components.innovation;
  filtered.covariance = components.predicted_covariance -
    gain * components.innovation_variance * gain.transpose();
  filtered.covariance = Scalar(0.5) *
    (filtered.covariance + filtered.covariance.transpose());
  components.filtered_mean = filtered.mean;
  components.filtered_covariance = filtered.covariance;
  if (components_output != nullptr) *components_output = components;
  return libertad::scalar_log(components.innovation_variance) +
    components.innovation * components.innovation /
      components.innovation_variance;
}

template <class Scalar>
Scalar kalman_row_nll_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data, int row, int subject,
    const std::vector<Scalar>& theta, const std::vector<Scalar>& eta,
    int eta_columns, const std::vector<Scalar>& sigma, int mixture_number,
    const Scalar& prediction, bool first, double previous_dv,
    double previous_time, KalmanFilterState<Scalar>& filtered,
    KalmanRowComponents<Scalar>* components_output = nullptr) {
  if (engine.kalman_filter_type != "linear") {
    return nonlinear_kalman_row_nll_t(
      engine, data, row, subject, theta, eta, eta_columns, sigma,
      mixture_number, prediction, first, previous_dv, previous_time,
      filtered, components_output);
  }
  const int states = engine.kalman_states;
  const std::size_t state_count = static_cast<std::size_t>(states);
  const std::size_t matrix_count = state_count * state_count;
  const std::size_t expected = 2U * state_count + 3U * matrix_count + 1U;
  if (!engine.kalman_enabled || states < 1 || engine.kalman_outputs.size() != expected) {
    throw std::logic_error("Compiled Kalman state-space likelihood is inconsistent.");
  }
  const std::vector<Scalar> output = evaluate_error_outputs_t(
    engine, data, row, subject, theta, eta, eta_columns, sigma,
    mixture_number, prediction, first, previous_dv, previous_time,
    engine.kalman_outputs, "Kalman state-space likelihood");
  std::size_t cursor = 0U;
  VectorT<Scalar> initial_mean(states);
  for (int state = 0; state < states; ++state) initial_mean[state] = output[cursor++];
  auto read_matrix = [&]() {
    MatrixT<Scalar> matrix(states, states);
    for (int matrix_row = 0; matrix_row < states; ++matrix_row) {
      for (int matrix_column = 0; matrix_column < states; ++matrix_column) {
        matrix(matrix_row, matrix_column) = output[cursor++];
      }
    }
    return matrix;
  };
  const MatrixT<Scalar> initial_covariance = read_matrix();
  const MatrixT<Scalar> transition = read_matrix();
  const MatrixT<Scalar> process_covariance = read_matrix();
  VectorT<Scalar> observation(states);
  for (int state = 0; state < states; ++state) observation[state] = output[cursor++];
  const Scalar observation_variance = output[cursor++];
  if (!(scalar_value(observation_variance) >= 0.0)) {
    throw std::domain_error(
      "Kalman observation variance must be non-negative at row " +
      std::to_string(row + 1) + ".");
  }
  KalmanRowComponents<Scalar> components;
  components.transition = transition;
  components.smoother_cross_covariance = MatrixT<Scalar>::Zero(states, states);
  components.observation = observation;
  components.observation_variance = observation_variance;
  if (first) {
    components.predicted_mean = initial_mean;
    components.predicted_covariance = Scalar(0.5) *
      (initial_covariance + initial_covariance.transpose());
  } else {
    if (filtered.mean.size() != states || filtered.covariance.rows() != states ||
        filtered.covariance.cols() != states) {
      throw std::logic_error("Kalman filter state has the wrong dimension.");
    }
    components.smoother_cross_covariance = filtered.covariance * transition.transpose();
    components.predicted_mean = transition * filtered.mean;
    components.predicted_covariance = transition * filtered.covariance *
      transition.transpose() + process_covariance;
    components.predicted_covariance = Scalar(0.5) *
      (components.predicted_covariance + components.predicted_covariance.transpose());
  }
  const Scalar baseline = engine.kalman_prediction_baseline ? prediction : Scalar(0.0);
  const Scalar observation_mean = baseline +
    (observation.transpose() * components.predicted_mean)[0];
  components.innovation = Scalar(row_optional(data, "DV", row, NA_REAL)) -
    observation_mean;
  components.innovation_variance =
    (observation.transpose() * components.predicted_covariance * observation)[0] +
    observation_variance;
  if (!(scalar_value(components.innovation_variance) > 1e-14)) {
    throw std::domain_error(
      "Kalman innovation variance is not positive at row " +
      std::to_string(row + 1) + ".");
  }
  const VectorT<Scalar> gain =
    components.predicted_covariance * observation / components.innovation_variance;
  filtered.mean = components.predicted_mean + gain * components.innovation;
  const MatrixT<Scalar> identity = MatrixT<Scalar>::Identity(states, states);
  const MatrixT<Scalar> update = identity - gain * observation.transpose();
  filtered.covariance = update * components.predicted_covariance * update.transpose() +
    gain * observation_variance * gain.transpose();
  filtered.covariance = Scalar(0.5) *
    (filtered.covariance + filtered.covariance.transpose());
  components.filtered_mean = filtered.mean;
  components.filtered_covariance = filtered.covariance;
  if (components_output != nullptr) *components_output = components;
  return libertad::scalar_log(components.innovation_variance) +
    components.innovation * components.innovation /
      components.innovation_variance;
}

template <class Scalar>
Scalar ar1_rho_t(const ModelEngine& engine,
                 const std::vector<Scalar>& theta,
                 const std::vector<Scalar>& sigma) {
  if (engine.ar1_parameter_source == "fixed") return Scalar(engine.ar1_rho);
  const std::vector<Scalar>& parameters =
    engine.ar1_parameter_source == "theta" ? theta : sigma;
  if (engine.ar1_parameter_index < 0 ||
      engine.ar1_parameter_index >= static_cast<int>(parameters.size())) {
    throw std::out_of_range("Estimated AR(1) parameter index is outside the parameter vector.");
  }
  Scalar rho = parameters[static_cast<std::size_t>(engine.ar1_parameter_index)];
  if (engine.ar1_transform == "tanh") rho = libertad::scalar_tanh(rho);
  if (!(std::abs(scalar_value(rho)) < 1.0)) {
    throw std::domain_error("Estimated AR(1) correlation must be strictly between -1 and 1.");
  }
  return rho;
}

template <class Scalar>
std::vector<Scalar> residual_subject_nll_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data,
    const std::vector<Scalar>& prediction, const std::vector<Scalar>& theta,
    const std::vector<Scalar>& eta, const std::vector<Scalar>& sigma,
    const std::vector<int>& mixture_assignment = std::vector<int>(),
    const std::vector<Scalar>* variance_prediction = nullptr) {
  if (!engine.residual_groups.empty()) {
    if (engine.error_type == "likelihood" || engine.sigma_correlation != "independent" ||
        engine.blq_method != "none") {
      throw std::logic_error("Correlated residual groups were combined with an incompatible likelihood option.");
    }
    return residual_grouped_subject_nll_t(
      engine, data, prediction, theta, sigma, variance_prediction);
  }
  Rcpp::NumericVector dv = data["DV"];
  Rcpp::IntegerVector evid = data["EVID"];
  Rcpp::IntegerVector mdv = data["MDV"];
  Rcpp::IntegerVector subjects = data[".ID_INDEX"];
  int n_subjects = 0;
  for (int value : subjects) n_subjects = std::max(n_subjects, value);
  std::vector<Scalar> result(static_cast<std::size_t>(n_subjects), Scalar(0.0));
  const bool has_dvid = data.containsElementNamed("DVID");
  const bool has_lloq = data.containsElementNamed("LLOQ");
  const bool has_blq = data.containsElementNamed("BLQ");
  const bool has_cens = data.containsElementNamed("CENS");
  int previous_subject = -1;
  std::unordered_map<int, Scalar> previous_standardized_residual;
  std::unordered_map<int, double> previous_outcome;
  std::unordered_map<int, double> previous_outcome_time;
  std::unordered_map<int, std::vector<Scalar>> hmm_filtered;
  std::unordered_map<int, KalmanFilterState<Scalar>> kalman_filtered;
  std::unordered_map<int, ParticleFilterState<Scalar>> particle_filtered;
  const int eta_columns = n_subjects > 0 ?
    static_cast<int>(eta.size() / static_cast<std::size_t>(n_subjects)) : 0;
  const Scalar ar1_rho = ar1_rho_t(engine, theta, sigma);

  for (int row = 0; row < data.nrows(); ++row) {
    const int subject = subjects[row] - 1;
    const int dvid = has_dvid ?
      static_cast<int>(row_optional(data, "DVID", row, 1.0)) : 1;
    if (subject != previous_subject) {
      previous_outcome.clear();
      previous_outcome_time.clear();
      hmm_filtered.clear();
      kalman_filtered.clear();
      particle_filtered.clear();
      previous_standardized_residual.clear();
    }
    if (evid[row] == 0 && mdv[row] == 0 && std::isfinite(dv[row])) {
      const Scalar f = prediction[static_cast<std::size_t>(row)];
      if (engine.error_type == "likelihood") {
        const int sequence =
          (engine.hmm_enabled && !engine.hmm_by_dvid) ||
          (engine.kalman_enabled && !engine.kalman_by_dvid) ? 1 : dvid;
        const auto previous = previous_outcome.find(sequence);
        const bool first = previous == previous_outcome.end();
        const double previous_value = first ? dv[row] : previous->second;
        const double previous_time = first ? row_optional(data, "TIME", row, 0.0) :
          previous_outcome_time.at(sequence);
        const int mixture_number = mixture_assignment.empty() ?
          static_cast<int>(row_optional(data, "MIXNUM", row, 1.0)) :
          mixture_assignment.at(static_cast<std::size_t>(subject));
        if (engine.hmm_enabled) {
          result[static_cast<std::size_t>(subject)] += hmm_row_nll_t(
            engine, data, row, subject, theta, eta, eta_columns, sigma,
            mixture_number, f, first, previous_value, previous_time,
            hmm_filtered[sequence]);
        } else if (engine.kalman_enabled) {
          if (engine.kalman_filter_type == "particle") {
            result[static_cast<std::size_t>(subject)] += particle_row_nll_t(
              engine, data, row, subject, theta, eta, eta_columns, sigma,
              mixture_number, f, first, previous_value, previous_time,
              particle_filtered[sequence]);
          } else {
            result[static_cast<std::size_t>(subject)] += kalman_row_nll_t(
              engine, data, row, subject, theta, eta, eta_columns, sigma,
              mixture_number, f, first, previous_value, previous_time,
              kalman_filtered[sequence]);
          }
        } else {
          result[static_cast<std::size_t>(subject)] += user_likelihood_nll_t(
            engine, data, row, subject, theta, eta, eta_columns, sigma,
            mixture_number, f, first, previous_value, previous_time);
        }
        previous_standardized_residual.erase(dvid);
      } else {
      const Scalar scale_prediction = variance_prediction == nullptr ? f :
        variance_prediction->at(static_cast<std::size_t>(row));
      const Scalar variance = residual_variance_t(
        engine, scale_prediction, sigma, dvid);
      const Scalar sd = libertad::scalar_sqrt(variance);
      double limit = has_lloq ? row_optional(data, "LLOQ", row, engine.lloq) : engine.lloq;
      bool censored = false;
      if (engine.blq_method != "none" && std::isfinite(limit)) {
        if (has_blq) censored = row_optional(data, "BLQ", row, 0.0) == 1.0;
        if (has_cens) censored = censored || row_optional(data, "CENS", row, 0.0) == 1.0;
        if (!has_blq && !has_cens) censored = dv[row] < limit;
      }
      if (censored) {
        Scalar z;
        if (engine.error_type == "exponential") {
          z = (Scalar(std::log(std::max(limit, 1e-300))) -
            libertad::scalar_log(scalar_floor_t(f, 1e-300))) / sd;
        } else {
          z = (Scalar(limit) - f) / sd;
        }
        Scalar probability = normal_cdf_t(z);
        if (engine.blq_method == "m4" && engine.error_type != "exponential") {
          const Scalar below_zero = normal_cdf_t((-f) / sd);
          probability = scalar_floor_t(
            (probability - below_zero) /
              scalar_floor_t(Scalar(1.0) - below_zero, 1e-300),
            1e-300);
        }
        result[static_cast<std::size_t>(subject)] -=
          Scalar(2.0) * libertad::scalar_log(probability);
        previous_standardized_residual.erase(dvid);
      } else {
        Scalar residual = Scalar(dv[row]) - f;
        if (engine.error_type == "exponential") {
          if (!(dv[row] > 0.0)) {
            throw std::domain_error("Exponential residual likelihood requires positive DV.");
          }
          residual = Scalar(std::log(dv[row])) -
            libertad::scalar_log(scalar_floor_t(f, 1e-300));
        }
        const auto previous_residual = previous_standardized_residual.find(dvid);
        if (engine.sigma_correlation == "ar1" &&
            previous_residual != previous_standardized_residual.end()) {
          const Scalar standardized = residual / sd;
          const Scalar innovation = standardized -
            ar1_rho * previous_residual->second;
          const Scalar innovation_variance =
            Scalar(1.0) - ar1_rho * ar1_rho;
          result[static_cast<std::size_t>(subject)] +=
            libertad::scalar_log(variance) +
            libertad::scalar_log(innovation_variance) +
            innovation * innovation / innovation_variance;
        } else {
          result[static_cast<std::size_t>(subject)] +=
            libertad::scalar_log(variance) + residual * residual / variance;
        }
        previous_standardized_residual[dvid] = residual / sd;
      }
      }
    }
    if (evid[row] == 0 && std::isfinite(dv[row]) &&
        ((!engine.hmm_enabled && !engine.kalman_enabled) || mdv[row] == 0)) {
      const int sequence =
        (engine.hmm_enabled && !engine.hmm_by_dvid) ||
        (engine.kalman_enabled && !engine.kalman_by_dvid) ? 1 : dvid;
      previous_outcome[sequence] = dv[row];
      previous_outcome_time[sequence] = row_optional(data, "TIME", row, 0.0);
    }
    previous_subject = subject;
  }
  return result;
}

template <class Scalar>
Scalar population_joint_nll_t(const ModelEngine& engine,
                              const Rcpp::DataFrame& data,
                              const std::vector<Scalar>& theta,
                              const std::vector<Scalar>& eta,
                              const std::vector<Scalar>& sigma,
                              const std::vector<Scalar>& omega,
                              bool interaction = true) {
  Rcpp::IntegerVector subjects = data[".ID_INDEX"];
  int n_subjects = 0;
  for (int value : subjects) n_subjects = std::max(n_subjects, value);
  Scalar total = Scalar(0.0);
  if (engine.mixture_probabilities.empty()) {
    const std::vector<Scalar> prediction = simulate_analytical_t(
      engine, data, theta, eta, sigma);
    std::vector<Scalar> variance_prediction;
    if (!interaction && engine.error_type != "likelihood") {
      const std::vector<Scalar> zero_eta(eta.size(), Scalar(0.0));
      variance_prediction = simulate_analytical_t(
        engine, data, theta, zero_eta, sigma);
    }
    const std::vector<Scalar> residual = residual_subject_nll_t(
      engine, data, prediction, theta, eta, sigma, std::vector<int>(),
      interaction || engine.error_type == "likelihood" ? nullptr : &variance_prediction);
    for (const Scalar& value : residual) total += value;
  } else {
    std::vector<std::vector<Scalar>> component_nll;
    component_nll.reserve(engine.mixture_probabilities.size());
    for (std::size_t component = 0; component < engine.mixture_probabilities.size(); ++component) {
      std::vector<int> assignment(static_cast<std::size_t>(n_subjects),
                                  static_cast<int>(component + 1));
      const std::vector<Scalar> prediction = simulate_analytical_t(
        engine, data, theta, eta, sigma, assignment);
      std::vector<Scalar> variance_prediction;
      if (!interaction && engine.error_type != "likelihood") {
        const std::vector<Scalar> zero_eta(eta.size(), Scalar(0.0));
        variance_prediction = simulate_analytical_t(
          engine, data, theta, zero_eta, sigma, assignment);
      }
      component_nll.push_back(residual_subject_nll_t(
        engine, data, prediction, theta, eta, sigma, assignment,
        interaction || engine.error_type == "likelihood" ? nullptr : &variance_prediction));
    }
    for (int subject = 0; subject < n_subjects; ++subject) {
      std::vector<Scalar> log_component(engine.mixture_probabilities.size());
      std::size_t maximum_index = 0;
      for (std::size_t component = 0; component < log_component.size(); ++component) {
        log_component[component] = Scalar(std::log(engine.mixture_probabilities[component])) -
          Scalar(0.5) * component_nll[component][static_cast<std::size_t>(subject)];
        if (scalar_value(log_component[component]) > scalar_value(log_component[maximum_index])) {
          maximum_index = component;
        }
      }
      const Scalar maximum = log_component[maximum_index];
      Scalar sum = Scalar(0.0);
      for (const Scalar& value : log_component) {
        sum += libertad::scalar_exp(value - maximum);
      }
      total -= Scalar(2.0) * (maximum + libertad::scalar_log(sum));
    }
  }

  if (engine.n_eta > 0) {
    const MatrixT<Scalar> covariance = omega_matrix_t(engine, omega);
    int n_subjects = 0;
    for (int value : subjects) n_subjects = std::max(n_subjects, value);
    if (n_subjects < 1 || eta.size() % static_cast<std::size_t>(n_subjects) != 0U) {
      throw std::invalid_argument("ETA vector cannot be divided into subject rows for its prior.");
    }
    const int eta_columns = static_cast<int>(eta.size() / static_cast<std::size_t>(n_subjects));
    for (int subject = 0; subject < n_subjects; ++subject) {
      if (engine.re_enabled) {
        int offset = 0;
        for (std::size_t block_index = 0; block_index < engine.re_blocks.size(); ++block_index) {
          const std::vector<int>& indices = engine.re_blocks[block_index];
          const std::string total_name = ".RE_TOTAL_" + std::to_string(block_index + 1U);
          if (!data.containsElementNamed(total_name.c_str())) {
            throw std::invalid_argument("Random-effect unit totals are missing from compiled data.");
          }
          const int units = static_cast<int>(data_value(data, total_name, 0));
          MatrixT<Scalar> block_covariance(indices.size(), indices.size());
          for (std::size_t row = 0; row < indices.size(); ++row) {
            for (std::size_t column = 0; column < indices.size(); ++column) {
              block_covariance(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) =
                covariance(indices[row], indices[column]);
            }
          }
          for (int unit = 0; unit < units; ++unit) {
            VectorT<Scalar> effect(static_cast<Eigen::Index>(indices.size()));
            for (std::size_t index = 0; index < indices.size(); ++index) {
              const int column = offset + unit * static_cast<int>(indices.size()) +
                static_cast<int>(index);
              if (column >= eta_columns) {
                throw std::invalid_argument("Expanded random-effect columns exceed the ETA matrix.");
              }
              effect[static_cast<Eigen::Index>(index)] =
                eta[static_cast<std::size_t>(subject * eta_columns + column)];
            }
            total += omega_subject_prior_t(block_covariance, effect);
          }
          offset += units * static_cast<int>(indices.size());
        }
        if (offset != eta_columns) {
          throw std::invalid_argument("Expanded random-effect columns do not match the design.");
        }
        continue;
      }
      if (engine.iov <= 0) {
        if (eta_columns != engine.n_eta) {
          throw std::invalid_argument("ETA matrix columns do not match OMEGA dimension.");
        }
        VectorT<Scalar> effect(engine.n_eta);
        for (int index = 0; index < engine.n_eta; ++index) {
          effect[index] = eta[static_cast<std::size_t>(subject * eta_columns + index)];
        }
        total += omega_subject_prior_t(covariance, effect);
        continue;
      }
      const int between = engine.n_eta - engine.iov;
      if (eta_columns < between || (eta_columns - between) % engine.iov != 0) {
        throw std::invalid_argument("IOV ETA columns do not match the occasion layout.");
      }
      if (between > 0) {
        MatrixT<Scalar> between_covariance = covariance.topLeftCorner(between, between);
        VectorT<Scalar> between_effect(between);
        for (int index = 0; index < between; ++index) {
          between_effect[index] = eta[static_cast<std::size_t>(subject * eta_columns + index)];
        }
        total += omega_subject_prior_t(between_covariance, between_effect);
      }
      const MatrixT<Scalar> iov_covariance = covariance.bottomRightCorner(engine.iov, engine.iov);
      const int occasions = (eta_columns - between) / engine.iov;
      for (int occasion = 0; occasion < occasions; ++occasion) {
        VectorT<Scalar> occasion_effect(engine.iov);
        for (int index = 0; index < engine.iov; ++index) {
          const int column = between + occasion * engine.iov + index;
          occasion_effect[index] = eta[static_cast<std::size_t>(subject * eta_columns + column)];
        }
        total += omega_subject_prior_t(iov_covariance, occasion_effect);
      }
    }
  }
  return total;
}

template <class Scalar>
Scalar positive_definite_gaussian_nll_t(
    const MatrixT<Scalar>& covariance, const VectorT<Scalar>& residual,
    const std::string& context) {
  const Eigen::Index dimension = covariance.rows();
  if (covariance.cols() != dimension || residual.size() != dimension) {
    throw std::invalid_argument(context + " dimensions are inconsistent.");
  }
  if (!dimension) return Scalar(0.0);
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
  VectorT<Scalar> forward(dimension);
  for (Eigen::Index row = 0; row < dimension; ++row) {
    Scalar value = residual[row];
    for (Eigen::Index column = 0; column < row; ++column) {
      value -= lower(row, column) * forward[column];
    }
    forward[row] = value / lower(row, row);
  }
  VectorT<Scalar> solution(dimension);
  for (Eigen::Index offset = 0; offset < dimension; ++offset) {
    const Eigen::Index row = dimension - offset - 1;
    Scalar value = forward[row];
    for (Eigen::Index column = row + 1; column < dimension; ++column) {
      value -= lower(column, row) * solution[column];
    }
    solution[row] = value / lower(row, row);
  }
  Scalar quadratic = Scalar(0.0);
  for (Eigen::Index row = 0; row < dimension; ++row) {
    quadratic += residual[row] * solution[row];
  }
  return logdet + quadratic;
}

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

std::unique_ptr<ObjectiveTape> record_fo_tape(
    const ModelEngine& engine, PredictionTape& prediction_tape,
    const Rcpp::DataFrame& data, const Rcpp::NumericVector& theta,
    const Rcpp::NumericVector& sigma, const Rcpp::NumericVector& omega) {
  const int n_theta = theta.size();
  const int n_sigma = sigma.size();
  const int n_omega = omega.size();
  const int n_eta = static_cast<int>(prediction_tape.domain_names.size()) -
    n_theta - n_sigma;
  if (n_eta < 0 || n_omega != static_cast<int>(engine.omega_rows.size())) {
    throw std::invalid_argument("FO tape parameter dimensions are inconsistent with the model.");
  }
  Rcpp::NumericVector dv = data["DV"];
  Rcpp::NumericVector dvid = data.containsElementNamed("DVID") ?
    Rcpp::NumericVector(data["DVID"]) : Rcpp::NumericVector(data.nrows(), 1.0);
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

  MatrixT<CppAD::AD<double>> residual_covariance(n_observed, n_observed);
  for (Eigen::Index row = 0; row < n_observed; ++row) {
    for (Eigen::Index column = 0; column < n_observed; ++column) {
      CppAD::AD<double> correlation = row == column ?
        CppAD::AD<double>(1.0) : CppAD::AD<double>(0.0);
      if (engine.sigma_correlation == "ar1" && dvid[observed[row]] == dvid[observed[column]]) {
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
        const int group_index = residual_group_for_dvid(engine, dvid[observed[row]]);
        if (group_index >= 0 && residual_group_for_dvid(engine, dvid[observed[column]]) == group_index) {
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
  std::vector<CppAD::AD<double>> dependent(1);
  dependent[0] = positive_definite_gaussian_nll_t(
    marginal, residual, "FO marginal covariance");
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
  return tape;
}

std::unique_ptr<ObjectiveTape> record_curvature_tape(
    const ModelEngine& engine, PredictionTape& prediction_tape,
    ObjectiveTape& objective_tape, const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta, const Rcpp::NumericVector& eta,
    const Rcpp::NumericVector& sigma, const Rcpp::NumericVector& omega,
    const std::string& approximation) {
  if (approximation != "foce" && approximation != "focei" &&
      approximation != "laplace") {
    throw std::invalid_argument("Unknown conditional-curvature approximation.");
  }
  const int n_theta = theta.size();
  const int n_eta = eta.size();
  const int n_sigma = sigma.size();
  const int n_omega = omega.size();
  std::vector<double> point;
  point.reserve(static_cast<std::size_t>(n_theta + n_eta + n_sigma + n_omega));
  for (double value : theta) point.push_back(value);
  for (double value : eta) point.push_back(value);
  for (double value : sigma) point.push_back(value);
  for (double value : omega) point.push_back(value);
  if (point.size() != objective_tape.domain_names.size() ||
      prediction_tape.domain_names.size() !=
        static_cast<std::size_t>(n_theta + n_eta + n_sigma)) {
    throw std::invalid_argument("Curvature tape parameter dimensions are inconsistent.");
  }
  std::vector<CppAD::AD<double>> independent(point.begin(), point.end());
  CppAD::Independent(independent);
  std::ostringstream messages;
  MatrixT<CppAD::AD<double>> curvature(n_eta, n_eta);

  if (approximation == "laplace") {
    auto objective_ad = objective_tape.fun.base2ad();
    objective_ad.Forward(0, independent, messages);
    std::vector<CppAD::AD<double>> direction(
      independent.size(), CppAD::AD<double>(0.0));
    const std::vector<CppAD::AD<double>> weight(1, CppAD::AD<double>(1.0));
    for (int column = 0; column < n_eta; ++column) {
      const std::size_t position = static_cast<std::size_t>(n_theta + column);
      direction[position] = CppAD::AD<double>(1.0);
      objective_ad.Forward(1, direction, messages);
      direction[position] = CppAD::AD<double>(0.0);
      const std::vector<CppAD::AD<double>> reverse =
        objective_ad.Reverse(2, weight);
      for (int row = 0; row < n_eta; ++row) {
        const std::size_t row_position = static_cast<std::size_t>(n_theta + row);
        curvature(row, column) = reverse[row_position * 2U + 1U];
      }
    }
    curvature = CppAD::AD<double>(0.5) *
      MatrixT<CppAD::AD<double>>(curvature + curvature.transpose());
  } else {
    std::vector<CppAD::AD<double>> prediction_point(
      independent.begin(), independent.begin() + n_theta + n_eta + n_sigma);
    auto prediction_ad = prediction_tape.fun.base2ad();
    const std::vector<CppAD::AD<double>> prediction =
      prediction_ad.Forward(0, prediction_point, messages);
    MatrixT<CppAD::AD<double>> eta_jacobian(
      static_cast<Eigen::Index>(prediction.size()), n_eta);
    std::vector<CppAD::AD<double>> direction(
      prediction_point.size(), CppAD::AD<double>(0.0));
    for (int column = 0; column < n_eta; ++column) {
      direction[static_cast<std::size_t>(n_theta + column)] = CppAD::AD<double>(1.0);
      const std::vector<CppAD::AD<double>> derivative =
        prediction_ad.Forward(1, direction, messages);
      direction[static_cast<std::size_t>(n_theta + column)] = CppAD::AD<double>(0.0);
      for (std::size_t row = 0; row < derivative.size(); ++row) {
        eta_jacobian(static_cast<Eigen::Index>(row), column) = derivative[row];
      }
    }
    std::vector<CppAD::AD<double>> scale_prediction = prediction;
    if (approximation == "foce") {
      std::vector<CppAD::AD<double>> zero_eta_point = prediction_point;
      for (int column = 0; column < n_eta; ++column) {
        zero_eta_point[static_cast<std::size_t>(n_theta + column)] =
          CppAD::AD<double>(0.0);
      }
      scale_prediction = prediction_ad.Forward(0, zero_eta_point, messages);
    }
    std::vector<CppAD::AD<double>> sigma_ad(
      independent.begin() + n_theta + n_eta,
      independent.begin() + n_theta + n_eta + n_sigma);
    Rcpp::NumericVector dv = data["DV"];
    Rcpp::NumericVector evid = data["EVID"];
    Rcpp::NumericVector mdv = data["MDV"];
    Rcpp::NumericVector dvid = data.containsElementNamed("DVID") ?
      Rcpp::NumericVector(data["DVID"]) : Rcpp::NumericVector(data.nrows(), 1.0);
    curvature.setZero();
    for (int row = 0; row < data.nrows(); ++row) {
      if (evid[row] != 0.0 || mdv[row] != 0.0 || !std::isfinite(dv[row])) continue;
      const CppAD::AD<double> variance = residual_variance_t(
        engine, scale_prediction[static_cast<std::size_t>(row)], sigma_ad,
        std::max(1, static_cast<int>(dvid[row])));
      for (int first = 0; first < n_eta; ++first) {
        for (int second = 0; second < n_eta; ++second) {
          curvature(first, second) += CppAD::AD<double>(2.0) *
            eta_jacobian(row, first) * eta_jacobian(row, second) / variance;
        }
      }
    }

    MatrixT<CppAD::AD<double>> base_omega =
      MatrixT<CppAD::AD<double>>::Zero(engine.n_eta, engine.n_eta);
    const std::size_t omega_offset = static_cast<std::size_t>(n_theta + n_eta + n_sigma);
    for (int index = 0; index < n_omega; ++index) {
      const int row = engine.omega_rows[static_cast<std::size_t>(index)];
      const int column = engine.omega_cols[static_cast<std::size_t>(index)];
      const CppAD::AD<double> value = independent[omega_offset + index];
      base_omega(row, column) = value;
      base_omega(column, row) = value;
    }
    MatrixT<CppAD::AD<double>> effect_omega = expanded_omega_t(
      engine, data, base_omega, n_eta);
    MatrixT<CppAD::AD<double>> identity =
      MatrixT<CppAD::AD<double>>::Identity(n_eta, n_eta);
    const MatrixT<CppAD::AD<double>> omega_inverse = solve_linear(
      effect_omega, identity, "Conditional OMEGA curvature");
    curvature += CppAD::AD<double>(2.0) * omega_inverse;
  }

  std::vector<CppAD::AD<double>> dependent(1);
  dependent[0] = positive_definite_logdet_t(
    curvature, "Conditional curvature determinant");
  auto tape = std::make_unique<ObjectiveTape>();
  tape->fun.Dependent(independent, dependent);
  tape->fun.optimize();
  tape->domain_names = objective_tape.domain_names;
  return tape;
}

std::unique_ptr<ObjectiveTape> record_objective_tape(
    const ModelEngine& engine, const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta, const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma, const Rcpp::NumericVector& omega,
    bool interaction) {
  std::vector<double> point = flatten_parameters(theta, eta, sigma);
  for (double value : omega) point.push_back(value);
  std::vector<CppAD::AD<double>> independent(point.begin(), point.end());
  CppAD::Independent(independent);
  std::size_t cursor = 0;
  std::vector<CppAD::AD<double>> theta_ad(static_cast<std::size_t>(theta.size()));
  for (auto& value : theta_ad) value = independent[cursor++];
  std::vector<CppAD::AD<double>> eta_ad(static_cast<std::size_t>(eta.size()));
  for (auto& value : eta_ad) value = independent[cursor++];
  std::vector<CppAD::AD<double>> sigma_ad(static_cast<std::size_t>(sigma.size()));
  for (auto& value : sigma_ad) value = independent[cursor++];
  std::vector<CppAD::AD<double>> omega_ad(static_cast<std::size_t>(omega.size()));
  for (auto& value : omega_ad) value = independent[cursor++];
  std::vector<CppAD::AD<double>> dependent(1);
  dependent[0] = population_joint_nll_t(
    engine, data, theta_ad, eta_ad, sigma_ad, omega_ad, interaction);
  auto tape = std::make_unique<ObjectiveTape>();
  tape->fun.Dependent(independent, dependent);
  tape->fun.optimize();
  tape->domain_names = parameter_names(theta.size(), eta.nrow(), eta.ncol(), sigma.size());
  for (int i = 0; i < omega.size(); ++i) {
    tape->domain_names.push_back("OMEGA_" + std::to_string(i + 1));
  }
  return tape;
}

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
      auto eigen = liberation::detail::self_adjoint_eigen(hessian, false);
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
          auto eigen = liberation::detail::self_adjoint_eigen(eta_hessian, false);
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

struct HmmDecodeObservation {
  int row = -1;
  std::vector<double> initial;
  std::vector<std::vector<double>> transition;
  std::vector<double> log_emission;
  std::vector<double> filtered;
  double log_scale = 0.0;
};

struct HmmDecodeSequence {
  int subject = 0;
  int sequence = 1;
  std::vector<HmmDecodeObservation> observations;
};

double hmm_log_probability(double value) {
  return value > 0.0 ? std::log(value) :
    -std::numeric_limits<double>::infinity();
}

double hmm_log_sum_exp(const std::vector<double>& values) {
  if (values.empty()) {
    throw std::invalid_argument("HMM log-sum-exp requires values.");
  }
  const double maximum = *std::max_element(values.begin(), values.end());
  if (!std::isfinite(maximum)) return maximum;
  double total = 0.0;
  for (double value : values) total += std::exp(value - maximum);
  return maximum + std::log(total);
}

Rcpp::List hmm_filter(
    const ModelEngine& engine, const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta, const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma) {
  if (!engine.hmm_enabled) {
    throw std::invalid_argument("The model does not define a hidden Markov likelihood.");
  }
  if (!engine.mixture_probabilities.empty()) {
    throw std::invalid_argument(
      "Hidden-state filtering for an additional finite-mixture layer is not yet available.");
  }
  Rcpp::IntegerVector subjects = data[".ID_INDEX"];
  Rcpp::IntegerVector evid = data["EVID"];
  Rcpp::IntegerVector mdv = data["MDV"];
  Rcpp::NumericVector dv = data["DV"];
  int n_subjects = 0;
  for (int value : subjects) n_subjects = std::max(n_subjects, value);
  if (eta.nrow() != n_subjects) {
    throw std::invalid_argument("The ETA matrix does not match the HMM dataset subjects.");
  }
  const std::vector<double> theta_values = Rcpp::as<std::vector<double>>(theta);
  std::vector<double> eta_values;
  eta_values.reserve(static_cast<std::size_t>(eta.size()));
  for (int row = 0; row < eta.nrow(); ++row) {
    for (int column = 0; column < eta.ncol(); ++column) {
      eta_values.push_back(eta(row, column));
    }
  }
  const std::vector<double> sigma_values = Rcpp::as<std::vector<double>>(sigma);
  const std::vector<double> prediction = simulate_analytical_t(
    engine, data, theta_values, eta_values, sigma_values);
  const int rows = data.nrows();
  const int states = engine.hmm_states;
  Rcpp::NumericMatrix probability(rows, states);
  Rcpp::NumericMatrix smoothed(rows, states);
  std::fill(probability.begin(), probability.end(), NA_REAL);
  std::fill(smoothed.begin(), smoothed.end(), NA_REAL);
  Rcpp::IntegerVector filtered_state(rows, NA_INTEGER);
  Rcpp::IntegerVector smoothed_state(rows, NA_INTEGER);
  Rcpp::IntegerVector viterbi_state(rows, NA_INTEGER);
  Rcpp::NumericVector row_nll(rows, NA_REAL);
  const bool has_dvid = data.containsElementNamed("DVID");
  int previous_subject = -1;
  std::unordered_map<int, double> previous_outcome;
  std::unordered_map<int, double> previous_outcome_time;
  std::unordered_map<int, std::vector<double>> filtered;
  std::unordered_map<int, std::size_t> sequence_lookup;
  std::vector<HmmDecodeSequence> sequences;
  double total_nll = 0.0;
  for (int row = 0; row < rows; ++row) {
    const int subject = subjects[row] - 1;
    if (subject != previous_subject) {
      previous_outcome.clear();
      previous_outcome_time.clear();
      filtered.clear();
      sequence_lookup.clear();
    }
    const int dvid = has_dvid ?
      static_cast<int>(row_optional(data, "DVID", row, 1.0)) : 1;
    const int sequence = engine.hmm_by_dvid ? dvid : 1;
    if (evid[row] == 0 && mdv[row] == 0 && std::isfinite(dv[row])) {
      const auto previous = previous_outcome.find(sequence);
      const bool first = previous == previous_outcome.end();
      const double previous_value = first ? dv[row] : previous->second;
      const double previous_time = first ? row_optional(data, "TIME", row, 0.0) :
        previous_outcome_time.at(sequence);
      const int mixture_number = static_cast<int>(
        row_optional(data, "MIXNUM", row, 1.0));
      HmmRowComponents<double> components;
      const double contribution = hmm_row_nll_t(
        engine, data, row, subject, theta_values, eta_values, eta.ncol(),
        sigma_values, mixture_number, prediction[static_cast<std::size_t>(row)],
        first, previous_value, previous_time, filtered[sequence], &components);
      row_nll[row] = contribution;
      total_nll += contribution;
      int most_likely = 0;
      for (int hidden = 0; hidden < states; ++hidden) {
        const double value = filtered[sequence][static_cast<std::size_t>(hidden)];
        probability(row, hidden) = value;
        if (value > filtered[sequence][static_cast<std::size_t>(most_likely)]) {
          most_likely = hidden;
        }
      }
      filtered_state[row] = most_likely + 1;
      if (first) {
        sequence_lookup[sequence] = sequences.size();
        HmmDecodeSequence track;
        track.subject = subject + 1;
        track.sequence = sequence;
        sequences.push_back(std::move(track));
      }
      const auto track = sequence_lookup.find(sequence);
      if (track == sequence_lookup.end()) {
        throw std::logic_error("Hidden Markov decoding sequence was not initialized.");
      }
      HmmDecodeObservation observation;
      observation.row = row;
      observation.initial = std::move(components.initial);
      observation.transition = std::move(components.transition);
      observation.log_emission = std::move(components.log_emission);
      if (engine.hmm_emission_scale != "log") {
        for (int hidden = 0; hidden < states; ++hidden) {
          if (!(components.emission[static_cast<std::size_t>(hidden)] > 0.0)) {
            observation.log_emission[static_cast<std::size_t>(hidden)] =
              -std::numeric_limits<double>::infinity();
          }
        }
      }
      observation.filtered = filtered[sequence];
      observation.log_scale = -0.5 * contribution;
      sequences[track->second].observations.push_back(std::move(observation));
      previous_outcome[sequence] = dv[row];
      previous_outcome_time[sequence] = row_optional(data, "TIME", row, 0.0);
    }
    previous_subject = subject;
  }

  Rcpp::IntegerVector summary_subject(static_cast<R_xlen_t>(sequences.size()));
  Rcpp::IntegerVector summary_sequence(static_cast<R_xlen_t>(sequences.size()));
  Rcpp::IntegerVector summary_observations(static_cast<R_xlen_t>(sequences.size()));
  Rcpp::NumericVector summary_log_likelihood(static_cast<R_xlen_t>(sequences.size()));
  Rcpp::NumericVector summary_viterbi_joint(static_cast<R_xlen_t>(sequences.size()));
  Rcpp::NumericVector summary_viterbi_posterior(static_cast<R_xlen_t>(sequences.size()));

  for (std::size_t sequence_index = 0; sequence_index < sequences.size();
       ++sequence_index) {
    const HmmDecodeSequence& sequence = sequences[sequence_index];
    const std::size_t length = sequence.observations.size();
    if (length == 0U) continue;
    const std::size_t state_count = static_cast<std::size_t>(states);

    std::vector<std::vector<double>> log_beta(
      length, std::vector<double>(state_count, 0.0));
    for (std::size_t current = length - 1U; current > 0U; --current) {
      const HmmDecodeObservation& next = sequence.observations[current];
      for (std::size_t from = 0; from < state_count; ++from) {
        std::vector<double> terms(state_count);
        for (std::size_t to = 0; to < state_count; ++to) {
          terms[to] = hmm_log_probability(next.transition[from][to]) +
            next.log_emission[to] + log_beta[current][to];
        }
        log_beta[current - 1U][from] =
          hmm_log_sum_exp(terms) - next.log_scale;
      }
    }
    for (std::size_t current = 0; current < length; ++current) {
      const HmmDecodeObservation& observation = sequence.observations[current];
      std::vector<double> log_probability(state_count);
      for (std::size_t hidden = 0; hidden < state_count; ++hidden) {
        log_probability[hidden] = hmm_log_probability(observation.filtered[hidden]) +
          log_beta[current][hidden];
      }
      const double normalizer = hmm_log_sum_exp(log_probability);
      int most_likely = 0;
      for (std::size_t hidden = 0; hidden < state_count; ++hidden) {
        const double value = std::exp(log_probability[hidden] - normalizer);
        smoothed(observation.row, static_cast<int>(hidden)) = value;
        if (value > smoothed(observation.row, most_likely)) {
          most_likely = static_cast<int>(hidden);
        }
      }
      smoothed_state[observation.row] = most_likely + 1;
    }

    std::vector<std::vector<double>> delta(
      length, std::vector<double>(state_count,
        -std::numeric_limits<double>::infinity()));
    std::vector<std::vector<int>> back_pointer(
      length, std::vector<int>(state_count, -1));
    const HmmDecodeObservation& first = sequence.observations.front();
    for (std::size_t hidden = 0; hidden < state_count; ++hidden) {
      delta[0][hidden] = hmm_log_probability(first.initial[hidden]) +
        first.log_emission[hidden];
    }
    for (std::size_t current = 1; current < length; ++current) {
      const HmmDecodeObservation& observation = sequence.observations[current];
      for (std::size_t to = 0; to < state_count; ++to) {
        int best_from = 0;
        double best = delta[current - 1U][0] +
          hmm_log_probability(observation.transition[0][to]);
        for (std::size_t from = 1; from < state_count; ++from) {
          const double candidate = delta[current - 1U][from] +
            hmm_log_probability(observation.transition[from][to]);
          if (candidate > best) {
            best = candidate;
            best_from = static_cast<int>(from);
          }
        }
        delta[current][to] = best + observation.log_emission[to];
        back_pointer[current][to] = best_from;
      }
    }
    int final_state = 0;
    for (int hidden = 1; hidden < states; ++hidden) {
      if (delta[length - 1U][static_cast<std::size_t>(hidden)] >
          delta[length - 1U][static_cast<std::size_t>(final_state)]) {
        final_state = hidden;
      }
    }
    std::vector<int> path(length, final_state);
    for (std::size_t current = length - 1U; current > 0U; --current) {
      path[current - 1U] = back_pointer[current]
        [static_cast<std::size_t>(path[current])];
    }
    for (std::size_t current = 0; current < length; ++current) {
      viterbi_state[sequence.observations[current].row] = path[current] + 1;
    }

    double sequence_log_likelihood = 0.0;
    for (const HmmDecodeObservation& observation : sequence.observations) {
      sequence_log_likelihood += observation.log_scale;
    }
    const double viterbi_log_joint =
      delta[length - 1U][static_cast<std::size_t>(final_state)];
    summary_subject[static_cast<R_xlen_t>(sequence_index)] = sequence.subject;
    summary_sequence[static_cast<R_xlen_t>(sequence_index)] = sequence.sequence;
    summary_observations[static_cast<R_xlen_t>(sequence_index)] =
      static_cast<int>(length);
    summary_log_likelihood[static_cast<R_xlen_t>(sequence_index)] =
      sequence_log_likelihood;
    summary_viterbi_joint[static_cast<R_xlen_t>(sequence_index)] =
      viterbi_log_joint;
    summary_viterbi_posterior[static_cast<R_xlen_t>(sequence_index)] =
      viterbi_log_joint - sequence_log_likelihood;
  }

  probability.attr("dimnames") = Rcpp::List::create(
    R_NilValue, Rcpp::wrap(engine.hmm_state_names));
  smoothed.attr("dimnames") = Rcpp::List::create(
    R_NilValue, Rcpp::wrap(engine.hmm_state_names));
  auto state_labels = [&](const Rcpp::IntegerVector& state) {
    Rcpp::CharacterVector labels(rows, NA_STRING);
    for (int row = 0; row < rows; ++row) {
      if (state[row] != NA_INTEGER) {
        labels[row] = engine.hmm_state_names[
          static_cast<std::size_t>(state[row] - 1)];
      }
    }
    return labels;
  };
  const Rcpp::CharacterVector filtered_label = state_labels(filtered_state);
  const Rcpp::CharacterVector smoothed_label = state_labels(smoothed_state);
  const Rcpp::CharacterVector viterbi_label = state_labels(viterbi_state);
  const Rcpp::DataFrame sequence_summary = Rcpp::DataFrame::create(
    Rcpp::Named(".ID_INDEX") = summary_subject,
    Rcpp::Named("HMM_SEQUENCE") = summary_sequence,
    Rcpp::Named("OBSERVATIONS") = summary_observations,
    Rcpp::Named("LOG_LIKELIHOOD") = summary_log_likelihood,
    Rcpp::Named("VITERBI_LOG_JOINT") = summary_viterbi_joint,
    Rcpp::Named("VITERBI_LOG_POSTERIOR") = summary_viterbi_posterior);
  return Rcpp::List::create(
    Rcpp::Named("filtered") = probability,
    Rcpp::Named("smoothed") = smoothed,
    Rcpp::Named("state") = filtered_state,
    Rcpp::Named("state_label") = filtered_label,
    Rcpp::Named("filtered_state") = filtered_state,
    Rcpp::Named("filtered_state_label") = filtered_label,
    Rcpp::Named("smoothed_state") = smoothed_state,
    Rcpp::Named("smoothed_state_label") = smoothed_label,
    Rcpp::Named("viterbi_state") = viterbi_state,
    Rcpp::Named("viterbi_state_label") = viterbi_label,
    Rcpp::Named("row_nll") = row_nll,
    Rcpp::Named("nll") = total_nll,
    Rcpp::Named("log_likelihood") = -0.5 * total_nll,
    Rcpp::Named("sequence_summary") = sequence_summary,
    Rcpp::Named("states") = engine.hmm_state_names);
}

struct KalmanDecodeObservation {
  int row = -1;
  Vector predicted_mean;
  Matrix predicted_covariance;
  Vector filtered_mean;
  Matrix filtered_covariance;
  Matrix transition;
  Matrix smoother_cross_covariance;
  std::vector<Vector> particle_states;
  std::vector<double> particle_weights;
  std::vector<int> particle_ancestors;
  std::vector<int> particle_regimes;
  double innovation = NA_REAL;
  double innovation_variance = NA_REAL;
  double row_nll = NA_REAL;
};

struct KalmanDecodeSequence {
  std::vector<KalmanDecodeObservation> observations;
};

Rcpp::List kalman_filter(
    const ModelEngine& engine, const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta, const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma) {
  if (!engine.kalman_enabled) {
    throw std::invalid_argument("The model does not define a Kalman state-space likelihood.");
  }
  if (!engine.mixture_probabilities.empty()) {
    throw std::invalid_argument(
      "Kalman decoding for an additional finite-mixture layer is not yet available.");
  }
  Rcpp::IntegerVector subjects = data[".ID_INDEX"];
  Rcpp::IntegerVector evid = data["EVID"];
  Rcpp::IntegerVector mdv = data["MDV"];
  Rcpp::NumericVector dv = data["DV"];
  int n_subjects = 0;
  for (int value : subjects) n_subjects = std::max(n_subjects, value);
  if (eta.nrow() != n_subjects) {
    throw std::invalid_argument("The ETA matrix does not match the Kalman dataset subjects.");
  }
  const std::vector<double> theta_values = Rcpp::as<std::vector<double>>(theta);
  std::vector<double> eta_values;
  eta_values.reserve(static_cast<std::size_t>(eta.size()));
  for (int row = 0; row < eta.nrow(); ++row) {
    for (int column = 0; column < eta.ncol(); ++column) {
      eta_values.push_back(eta(row, column));
    }
  }
  const std::vector<double> sigma_values = Rcpp::as<std::vector<double>>(sigma);
  const std::vector<double> prediction = simulate_analytical_t(
    engine, data, theta_values, eta_values, sigma_values);
  const int rows = data.nrows();
  const int states = engine.kalman_states;
  Rcpp::NumericMatrix predicted_mean(rows, states);
  Rcpp::NumericMatrix filtered_mean(rows, states);
  Rcpp::NumericMatrix smoothed_mean(rows, states);
  Rcpp::NumericMatrix filtered_variance(rows, states);
  Rcpp::NumericMatrix smoothed_variance(rows, states);
  std::fill(predicted_mean.begin(), predicted_mean.end(), NA_REAL);
  std::fill(filtered_mean.begin(), filtered_mean.end(), NA_REAL);
  std::fill(smoothed_mean.begin(), smoothed_mean.end(), NA_REAL);
  std::fill(filtered_variance.begin(), filtered_variance.end(), NA_REAL);
  std::fill(smoothed_variance.begin(), smoothed_variance.end(), NA_REAL);
  Rcpp::NumericVector innovation(rows, NA_REAL);
  Rcpp::NumericVector innovation_variance(rows, NA_REAL);
  Rcpp::NumericVector row_nll(rows, NA_REAL);
  Rcpp::NumericMatrix filtered_regime(
    rows, engine.switching_enabled ? engine.switching_regimes : 0);
  Rcpp::NumericMatrix smoothed_regime(
    rows, engine.switching_enabled ? engine.switching_regimes : 0);
  std::fill(filtered_regime.begin(), filtered_regime.end(), NA_REAL);
  std::fill(smoothed_regime.begin(), smoothed_regime.end(), NA_REAL);
  const bool has_dvid = data.containsElementNamed("DVID");
  int previous_subject = -1;
  std::unordered_map<int, double> previous_outcome;
  std::unordered_map<int, double> previous_outcome_time;
  std::unordered_map<int, KalmanFilterState<double>> filtered;
  std::unordered_map<int, ParticleFilterState<double>> particle_filter;
  std::unordered_map<int, std::size_t> sequence_lookup;
  std::vector<KalmanDecodeSequence> sequences;
  double total_nll = 0.0;
  for (int row = 0; row < rows; ++row) {
    const int subject = subjects[row] - 1;
    if (subject != previous_subject) {
      previous_outcome.clear();
      previous_outcome_time.clear();
      filtered.clear();
      particle_filter.clear();
      sequence_lookup.clear();
    }
    const int dvid = has_dvid ?
      static_cast<int>(row_optional(data, "DVID", row, 1.0)) : 1;
    const int sequence = engine.kalman_by_dvid ? dvid : 1;
    if (evid[row] == 0 && mdv[row] == 0 && std::isfinite(dv[row])) {
      const auto previous = previous_outcome.find(sequence);
      const bool first = previous == previous_outcome.end();
      const double previous_value = first ? dv[row] : previous->second;
      const double previous_time = first ? row_optional(data, "TIME", row, 0.0) :
        previous_outcome_time.at(sequence);
      const int mixture_number = static_cast<int>(
        row_optional(data, "MIXNUM", row, 1.0));
      KalmanRowComponents<double> components;
      const double contribution = engine.kalman_filter_type == "particle" ?
        particle_row_nll_t(
          engine, data, row, subject, theta_values, eta_values, eta.ncol(),
          sigma_values, mixture_number, prediction[static_cast<std::size_t>(row)],
          first, previous_value, previous_time, particle_filter[sequence], &components) :
        kalman_row_nll_t(
          engine, data, row, subject, theta_values, eta_values, eta.ncol(),
          sigma_values, mixture_number, prediction[static_cast<std::size_t>(row)],
          first, previous_value, previous_time, filtered[sequence], &components);
      if (first) {
        sequence_lookup[sequence] = sequences.size();
        sequences.emplace_back();
      }
      const auto track = sequence_lookup.find(sequence);
      if (track == sequence_lookup.end()) {
        throw std::logic_error("Kalman decoding sequence was not initialized.");
      }
      KalmanDecodeObservation observation;
      observation.row = row;
      observation.predicted_mean = components.predicted_mean;
      observation.predicted_covariance = components.predicted_covariance;
      observation.filtered_mean = components.filtered_mean;
      observation.filtered_covariance = components.filtered_covariance;
      observation.transition = components.transition;
      observation.smoother_cross_covariance = components.smoother_cross_covariance;
      observation.particle_states = components.particle_states;
      observation.particle_weights = components.particle_weights;
      observation.particle_ancestors = components.particle_ancestors;
      observation.particle_regimes = components.particle_regimes;
      observation.innovation = components.innovation;
      observation.innovation_variance = components.innovation_variance;
      observation.row_nll = contribution;
      sequences[track->second].observations.push_back(std::move(observation));
      for (int state = 0; state < states; ++state) {
        predicted_mean(row, state) = components.predicted_mean[state];
        filtered_mean(row, state) = components.filtered_mean[state];
        filtered_variance(row, state) = components.filtered_covariance(state, state);
      }
      if (engine.switching_enabled) {
        for (int regime = 0; regime < engine.switching_regimes; ++regime) {
          double probability = 0.0;
          for (std::size_t particle = 0;
               particle < components.particle_regimes.size(); ++particle) {
            if (components.particle_regimes[particle] == regime &&
                particle < components.particle_weights.size()) {
              probability += components.particle_weights[particle];
            }
          }
          filtered_regime(row, regime) = probability;
        }
      }
      innovation[row] = components.innovation;
      innovation_variance[row] = components.innovation_variance;
      row_nll[row] = contribution;
      total_nll += contribution;
      previous_outcome[sequence] = dv[row];
      previous_outcome_time[sequence] = row_optional(data, "TIME", row, 0.0);
    }
    previous_subject = subject;
  }
  for (const KalmanDecodeSequence& sequence : sequences) {
    const std::size_t length = sequence.observations.size();
    if (!length) continue;
    std::vector<Vector> smooth_mean(length);
    std::vector<Matrix> smooth_covariance(length);
    if (engine.kalman_filter_type == "particle") {
      std::vector<std::vector<double>> smooth_weight(length);
      smooth_weight[length - 1U] = sequence.observations[length - 1U].particle_weights;
      for (std::size_t offset = 1U; offset < length; ++offset) {
        const std::size_t current = length - offset - 1U;
        const KalmanDecodeObservation& next = sequence.observations[current + 1U];
        const std::size_t particles = sequence.observations[current].particle_states.size();
        smooth_weight[current].assign(particles, 0.0);
        for (std::size_t particle = 0; particle < next.particle_ancestors.size(); ++particle) {
          const int ancestor = next.particle_ancestors[particle];
          if (ancestor >= 0 && static_cast<std::size_t>(ancestor) < particles &&
              particle < smooth_weight[current + 1U].size()) {
            smooth_weight[current][static_cast<std::size_t>(ancestor)] +=
              smooth_weight[current + 1U][particle];
          }
        }
        const double total = std::accumulate(
          smooth_weight[current].begin(), smooth_weight[current].end(), 0.0);
        if (total > 0.0) for (double& weight : smooth_weight[current]) weight /= total;
      }
      for (std::size_t current = 0; current < length; ++current) {
        const auto& cloud = sequence.observations[current].particle_states;
        const auto& weights = smooth_weight[current];
        smooth_mean[current] = Vector::Zero(states);
        for (std::size_t particle = 0; particle < cloud.size(); ++particle) {
          smooth_mean[current] += weights[particle] * cloud[particle];
        }
        smooth_covariance[current] = Matrix::Zero(states, states);
        for (std::size_t particle = 0; particle < cloud.size(); ++particle) {
          const Vector centered = cloud[particle] - smooth_mean[current];
          smooth_covariance[current] += weights[particle] * centered * centered.transpose();
        }
        if (engine.switching_enabled) {
          const int row = sequence.observations[current].row;
          const auto& regimes = sequence.observations[current].particle_regimes;
          for (int regime = 0; regime < engine.switching_regimes; ++regime) {
            double probability = 0.0;
            for (std::size_t particle = 0; particle < regimes.size() &&
                 particle < weights.size(); ++particle) {
              if (regimes[particle] == regime) probability += weights[particle];
            }
            smoothed_regime(row, regime) = probability;
          }
        }
      }
    } else {
      smooth_mean[length - 1U] = sequence.observations[length - 1U].filtered_mean;
      smooth_covariance[length - 1U] =
        sequence.observations[length - 1U].filtered_covariance;
      for (std::size_t offset = 1U; offset < length; ++offset) {
        const std::size_t current = length - offset - 1U;
        const KalmanDecodeObservation& observation = sequence.observations[current];
        const KalmanDecodeObservation& next = sequence.observations[current + 1U];
        Eigen::LDLT<Matrix> decomposition(next.predicted_covariance);
        if (decomposition.info() != Eigen::Success) {
          throw std::runtime_error("RTS smoother predicted covariance factorization failed.");
        }
        const Matrix smoother_gain = next.smoother_cross_covariance *
          decomposition.solve(Matrix::Identity(states, states));
        smooth_mean[current] = observation.filtered_mean + smoother_gain *
          (smooth_mean[current + 1U] - next.predicted_mean);
        smooth_covariance[current] = observation.filtered_covariance + smoother_gain *
          (smooth_covariance[current + 1U] - next.predicted_covariance) *
            smoother_gain.transpose();
        smooth_covariance[current] = 0.5 *
          (smooth_covariance[current] + smooth_covariance[current].transpose()).eval();
      }
    }
    for (std::size_t current = 0; current < length; ++current) {
      const int row = sequence.observations[current].row;
      for (int state = 0; state < states; ++state) {
        smoothed_mean(row, state) = smooth_mean[current][state];
        smoothed_variance(row, state) = smooth_covariance[current](state, state);
      }
    }
  }
  const Rcpp::List dimnames = Rcpp::List::create(
    R_NilValue, Rcpp::wrap(engine.kalman_state_names));
  predicted_mean.attr("dimnames") = dimnames;
  filtered_mean.attr("dimnames") = dimnames;
  smoothed_mean.attr("dimnames") = dimnames;
  filtered_variance.attr("dimnames") = dimnames;
  smoothed_variance.attr("dimnames") = dimnames;
  if (engine.switching_enabled) {
    const Rcpp::List regime_dimnames = Rcpp::List::create(
      R_NilValue, Rcpp::wrap(engine.switching_regime_names));
    filtered_regime.attr("dimnames") = regime_dimnames;
    smoothed_regime.attr("dimnames") = regime_dimnames;
  }
  return Rcpp::List::create(
    Rcpp::Named("predicted_mean") = predicted_mean,
    Rcpp::Named("filtered_mean") = filtered_mean,
    Rcpp::Named("smoothed_mean") = smoothed_mean,
    Rcpp::Named("filtered_variance") = filtered_variance,
    Rcpp::Named("smoothed_variance") = smoothed_variance,
    Rcpp::Named("innovation") = innovation,
    Rcpp::Named("innovation_variance") = innovation_variance,
    Rcpp::Named("row_nll") = row_nll,
    Rcpp::Named("filtered_regime") = filtered_regime,
    Rcpp::Named("smoothed_regime") = smoothed_regime,
    Rcpp::Named("regimes") = engine.switching_regime_names,
    Rcpp::Named("nll") = total_nll,
    Rcpp::Named("log_likelihood") = -0.5 * total_nll,
    Rcpp::Named("states") = engine.kalman_state_names,
    Rcpp::Named("filter") = engine.kalman_filter_type,
    Rcpp::Named("smoother") = engine.kalman_filter_type == "particle" ?
      "genealogical" : "RTS");
}

Matrix covariance_sampling_root(const Matrix& covariance,
                                const std::string& context) {
  const Matrix symmetric = 0.5 * (covariance + covariance.transpose()).eval();
  auto eigen = liberation::detail::self_adjoint_eigen(symmetric, true);
  if (eigen.info != Eigen::Success) {
    throw std::runtime_error(context + " eigen decomposition failed.");
  }
  const double scale = std::max(eigen.values.cwiseAbs().maxCoeff(), 1.0);
  if (eigen.values.minCoeff() < -1e-10 * scale) {
    throw std::domain_error(context + " is not positive semidefinite.");
  }
  return eigen.vectors * eigen.values.cwiseMax(0.0).cwiseSqrt().asDiagonal();
}

Rcpp::NumericVector kalman_simulate(
    const ModelEngine& engine, const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta, const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma,
    const Rcpp::NumericMatrix& process_normals,
    const Rcpp::NumericVector& observation_normals) {
  if (!engine.kalman_enabled) {
    throw std::invalid_argument("The model does not define a Kalman state-space likelihood.");
  }
  const int rows = data.nrows();
  const int states = engine.kalman_states;
  const int process_columns = states *
    (engine.kalman_dynamics == "sde" ? engine.kalman_sde_substeps : 1);
  if (process_normals.nrow() != rows || process_normals.ncol() != process_columns ||
      observation_normals.size() != rows) {
    throw std::invalid_argument("Kalman simulation normal draws have inconsistent dimensions.");
  }
  Rcpp::IntegerVector subjects = data[".ID_INDEX"];
  Rcpp::IntegerVector evid = data["EVID"];
  Rcpp::IntegerVector mdv = data["MDV"];
  Rcpp::NumericVector input_dv = data["DV"];
  Rcpp::NumericVector simulated = Rcpp::clone(input_dv);
  int n_subjects = 0;
  for (int value : subjects) n_subjects = std::max(n_subjects, value);
  if (eta.nrow() != n_subjects) {
    throw std::invalid_argument("The ETA matrix does not match the Kalman simulation subjects.");
  }
  const std::vector<double> theta_values = Rcpp::as<std::vector<double>>(theta);
  std::vector<double> eta_values;
  eta_values.reserve(static_cast<std::size_t>(eta.size()));
  for (int row = 0; row < eta.nrow(); ++row) {
    for (int column = 0; column < eta.ncol(); ++column) eta_values.push_back(eta(row, column));
  }
  const std::vector<double> sigma_values = Rcpp::as<std::vector<double>>(sigma);
  const std::vector<double> prediction = simulate_analytical_t(
    engine, data, theta_values, eta_values, sigma_values);
  const bool has_dvid = data.containsElementNamed("DVID");
  int previous_subject = -1;
  std::unordered_map<int, Vector> latent;
  std::unordered_map<int, int> latent_regime;
  std::unordered_map<int, double> previous_time;
  std::unordered_map<int, double> previous_dv;
  for (int row = 0; row < rows; ++row) {
    const int subject = subjects[row] - 1;
    if (subject != previous_subject) {
      latent.clear();
      latent_regime.clear();
      previous_time.clear();
      previous_dv.clear();
    }
    if (evid[row] != 0 || mdv[row] != 0) {
      previous_subject = subject;
      continue;
    }
    const int dvid = has_dvid ?
      static_cast<int>(row_optional(data, "DVID", row, 1.0)) : 1;
    const int sequence = engine.kalman_by_dvid ? dvid : 1;
    const bool first = latent.find(sequence) == latent.end();
    const double previous_value = first ? 0.0 : previous_dv.at(sequence);
    const double previous_observation_time = first ?
      row_optional(data, "TIME", row, 0.0) : previous_time.at(sequence);
    const int mixture_number = static_cast<int>(row_optional(data, "MIXNUM", row, 1.0));
    Vector normal(states);
    for (int state = 0; state < states; ++state) normal[state] = process_normals(row, state);
    if (engine.kalman_filter_type != "linear") {
      const Vector evaluation_state = first ? Vector::Zero(states) : latent[sequence];
      NonlinearStateOutputs<double> output = nonlinear_state_outputs_t(
        engine, data, row, subject, theta_values, eta_values, eta.ncol(),
        sigma_values, mixture_number, prediction[static_cast<std::size_t>(row)],
        first, previous_value, previous_observation_time, evaluation_state);
      if (engine.switching_enabled) {
        const SwitchingStateOutputs<double> switching = switching_state_raw_outputs_t(
          engine, data, row, subject, theta_values, eta_values, eta.ncol(),
          sigma_values, mixture_number, prediction[static_cast<std::size_t>(row)],
          first, previous_value, previous_observation_time, evaluation_state);
        std::vector<double> probability;
        if (first) {
          probability = switching.initial;
        } else {
          probability.resize(static_cast<std::size_t>(engine.switching_regimes));
          const int from = latent_regime.at(sequence);
          for (int to = 0; to < engine.switching_regimes; ++to) {
            probability[static_cast<std::size_t>(to)] = switching.regime_transition(from, to);
          }
        }
        latent_regime[sequence] = sample_switching_regime(
          probability, state_space_uniform(
            engine.kalman_seed + 65537 * subject, row, sequence, 0, first ? 31 : 32));
        const NonlinearStateOutputs<double>& local =
          switching.regime[static_cast<std::size_t>(latent_regime[sequence])];
        output.transition = local.transition;
        output.process_covariance = local.process_covariance;
        output.observation = local.observation;
        output.observation_variance = local.observation_variance;
      }
      if (!(output.observation_variance >= 0.0)) {
        throw std::domain_error("Nonlinear state-space simulation observation variance must be non-negative.");
      }
      if (first) {
        latent[sequence] = output.initial_mean + covariance_sampling_root(
          output.initial_covariance, "State-space initial covariance") * normal;
      } else if (engine.kalman_dynamics == "sde") {
        Vector current = latent[sequence];
        const double current_time = row_optional(data, "TIME", row, previous_observation_time);
        const double interval = std::max(0.0, current_time - previous_observation_time);
        const double step = interval / static_cast<double>(engine.kalman_sde_substeps);
        for (int substep = 0; substep < engine.kalman_sde_substeps; ++substep) {
          const NonlinearStateOutputs<double> local = engine.switching_enabled ?
            switching_state_raw_outputs_t(
              engine, data, row, subject, theta_values, eta_values, eta.ncol(),
              sigma_values, mixture_number, prediction[static_cast<std::size_t>(row)],
              false, previous_value, previous_observation_time + substep * step,
              current).regime[static_cast<std::size_t>(latent_regime[sequence])] :
            nonlinear_state_raw_outputs_t(
              engine, data, row, subject, theta_values, eta_values, eta.ncol(),
              sigma_values, mixture_number, prediction[static_cast<std::size_t>(row)],
              false, previous_value, previous_observation_time + substep * step, current);
          Vector increment(states);
          for (int state = 0; state < states; ++state) {
            increment[state] = std::sqrt(step) *
              process_normals(row, substep * states + state);
          }
          Vector next = current + step * local.transition +
            local.process_covariance * increment;
          if (engine.kalman_sde_method == "milstein") {
            for (int first_state = 0; first_state < states; ++first_state) {
              for (int second_state = 0; second_state < states; ++second_state) {
                if (first_state != second_state &&
                    std::abs(local.process_covariance(first_state, second_state)) > 1e-12) {
                  throw std::domain_error(
                    "Milstein simulation currently requires diagonal diffusion.");
                }
              }
              const double derivative_step = engine.kalman_jacobian_step *
                std::max(1.0, std::abs(current[first_state]));
              Vector plus = current;
              Vector minus = current;
              plus[first_state] += derivative_step;
              minus[first_state] -= derivative_step;
              const double upper = engine.switching_enabled ?
                switching_state_raw_outputs_t(
                  engine, data, row, subject, theta_values, eta_values, eta.ncol(),
                  sigma_values, mixture_number, prediction[static_cast<std::size_t>(row)],
                  false, previous_value, previous_observation_time + substep * step,
                  plus).regime[static_cast<std::size_t>(latent_regime[sequence])].
                    process_covariance(first_state, first_state) :
                nonlinear_state_raw_outputs_t(
                  engine, data, row, subject, theta_values, eta_values, eta.ncol(),
                  sigma_values, mixture_number, prediction[static_cast<std::size_t>(row)],
                  false, previous_value, previous_observation_time + substep * step,
                  plus).process_covariance(first_state, first_state);
              const double lower = engine.switching_enabled ?
                switching_state_raw_outputs_t(
                  engine, data, row, subject, theta_values, eta_values, eta.ncol(),
                  sigma_values, mixture_number, prediction[static_cast<std::size_t>(row)],
                  false, previous_value, previous_observation_time + substep * step,
                  minus).regime[static_cast<std::size_t>(latent_regime[sequence])].
                    process_covariance(first_state, first_state) :
                nonlinear_state_raw_outputs_t(
                  engine, data, row, subject, theta_values, eta_values, eta.ncol(),
                  sigma_values, mixture_number, prediction[static_cast<std::size_t>(row)],
                  false, previous_value, previous_observation_time + substep * step,
                  minus).process_covariance(first_state, first_state);
              const double derivative = (upper - lower) / (2.0 * derivative_step);
              const double diffusion = local.process_covariance(first_state, first_state);
              next[first_state] += 0.5 * diffusion * derivative *
                (increment[first_state] * increment[first_state] - step);
            }
          }
          current = next;
        }
        latent[sequence] = current;
      } else {
        latent[sequence] = output.transition + covariance_sampling_root(
          output.process_covariance, "State-space process covariance") * normal;
      }
      const NonlinearStateOutputs<double> observed = engine.switching_enabled ?
        switching_state_raw_outputs_t(
          engine, data, row, subject, theta_values, eta_values, eta.ncol(),
          sigma_values, mixture_number, prediction[static_cast<std::size_t>(row)],
          first, previous_value, previous_observation_time,
          latent[sequence]).regime[static_cast<std::size_t>(latent_regime[sequence])] :
        nonlinear_state_outputs_t(
          engine, data, row, subject, theta_values, eta_values, eta.ncol(),
          sigma_values, mixture_number, prediction[static_cast<std::size_t>(row)],
          first, previous_value, previous_observation_time, latent[sequence]);
      const double baseline = engine.kalman_prediction_baseline ?
        prediction[static_cast<std::size_t>(row)] : 0.0;
      simulated[row] = baseline + observed.observation +
        std::sqrt(observed.observation_variance) * observation_normals[row];
      previous_time[sequence] = row_optional(data, "TIME", row, 0.0);
      previous_dv[sequence] = simulated[row];
      previous_subject = subject;
      continue;
    }
    const std::vector<double> output = evaluate_error_outputs_t(
      engine, data, row, subject, theta_values, eta_values, eta.ncol(),
      sigma_values, mixture_number, prediction[static_cast<std::size_t>(row)],
      first, previous_value, previous_observation_time, engine.kalman_outputs,
      "Kalman state-space simulation");
    std::size_t cursor = 0U;
    Vector initial_mean(states);
    for (int state = 0; state < states; ++state) initial_mean[state] = output[cursor++];
    auto read_matrix = [&]() {
      Matrix matrix(states, states);
      for (int matrix_row = 0; matrix_row < states; ++matrix_row) {
        for (int matrix_column = 0; matrix_column < states; ++matrix_column) {
          matrix(matrix_row, matrix_column) = output[cursor++];
        }
      }
      return matrix;
    };
    const Matrix initial_covariance = read_matrix();
    const Matrix transition = read_matrix();
    const Matrix process_covariance = read_matrix();
    Vector observation(states);
    for (int state = 0; state < states; ++state) observation[state] = output[cursor++];
    const double observation_variance = output[cursor++];
    if (!(observation_variance >= 0.0)) {
      throw std::domain_error("Kalman simulation observation variance must be non-negative.");
    }
    if (first) {
      latent[sequence] = initial_mean + covariance_sampling_root(
        initial_covariance, "Kalman initial covariance") * normal;
    } else {
      latent[sequence] = transition * latent[sequence] + covariance_sampling_root(
        process_covariance, "Kalman process covariance") * normal;
    }
    const double baseline = engine.kalman_prediction_baseline ?
      prediction[static_cast<std::size_t>(row)] : 0.0;
    simulated[row] = baseline + observation.dot(latent[sequence]) +
      std::sqrt(observation_variance) * observation_normals[row];
    previous_time[sequence] = row_optional(data, "TIME", row, 0.0);
    previous_dv[sequence] = simulated[row];
    previous_subject = subject;
  }
  return simulated;
}

#include "hmc_sampler.h"

}  // namespace liberation

// [[Rcpp::export(name = ".liberation_population_objective_create")]]
SEXP liberation_population_objective_create(
    SEXP engine_pointer, const Rcpp::List& subject_data,
    const Rcpp::List& primary_tape_pointers,
    const Rcpp::List& curvature_tape_pointers,
    const Rcpp::List& config) {
  Rcpp::XPtr<liberation::PopulationObjective> pointer(
    new liberation::PopulationObjective(
      engine_pointer, subject_data, primary_tape_pointers,
      curvature_tape_pointers, config),
    true
  );
  pointer.attr("class") = Rcpp::CharacterVector::create(
    "liberation_population_objective_ptr", "externalptr");
  // Keep every non-owning analytical tape, the model engine, and the subject
  // data alive for exactly as long as the persistent evaluator is reachable.
  pointer.attr("keepers") = Rcpp::List::create(
    engine_pointer, subject_data, primary_tape_pointers,
    curvature_tape_pointers);
  return pointer;
}

// [[Rcpp::export(name = ".liberation_population_objective_value")]]
double liberation_population_objective_value(
    SEXP pointer, const Rcpp::NumericVector& encoded) {
  Rcpp::XPtr<liberation::PopulationObjective> objective(pointer);
  return objective->value(encoded);
}

// [[Rcpp::export(name = ".liberation_population_objective_gradient")]]
Rcpp::NumericVector liberation_population_objective_gradient(
    SEXP pointer, const Rcpp::NumericVector& encoded) {
  Rcpp::XPtr<liberation::PopulationObjective> objective(pointer);
  return objective->gradient(encoded);
}

// [[Rcpp::export(name = ".liberation_population_objective_state")]]
Rcpp::List liberation_population_objective_state(
    SEXP pointer, const Rcpp::NumericVector& encoded) {
  Rcpp::XPtr<liberation::PopulationObjective> objective(pointer);
  return objective->state(encoded);
}

// [[Rcpp::export(name = ".liberation_population_objective_telemetry")]]
Rcpp::List liberation_population_objective_telemetry(SEXP pointer) {
  Rcpp::XPtr<liberation::PopulationObjective> objective(pointer);
  return objective->telemetry();
}

// [[Rcpp::export(name = ".liberation_engine_create")]]
SEXP liberation_engine_create(const Rcpp::List& specification) {
  Rcpp::XPtr<liberation::ModelEngine> pointer(
    new liberation::ModelEngine(specification), true
  );
  pointer.attr("class") = Rcpp::CharacterVector::create("liberation_engine_ptr", "externalptr");
  return pointer;
}

// [[Rcpp::export(name = ".liberation_engine_simulate")]]
Rcpp::List liberation_engine_simulate(
    SEXP engine_pointer,
    const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta,
    const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma) {
  Rcpp::XPtr<liberation::ModelEngine> engine(engine_pointer);
  return liberation::simulate(*engine, data, theta, eta, sigma);
}

// [[Rcpp::export(name = ".liberation_engine_hmm_filter")]]
Rcpp::List liberation_engine_hmm_filter(
    SEXP engine_pointer,
    const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta,
    const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma) {
  Rcpp::XPtr<liberation::ModelEngine> engine(engine_pointer);
  return liberation::hmm_filter(*engine, data, theta, eta, sigma);
}

// [[Rcpp::export(name = ".liberation_engine_kalman_filter")]]
Rcpp::List liberation_engine_kalman_filter(
    SEXP engine_pointer,
    const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta,
    const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma) {
  Rcpp::XPtr<liberation::ModelEngine> engine(engine_pointer);
  return liberation::kalman_filter(*engine, data, theta, eta, sigma);
}

// [[Rcpp::export(name = ".liberation_engine_kalman_simulate")]]
Rcpp::NumericVector liberation_engine_kalman_simulate(
    SEXP engine_pointer, const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta, const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma,
    const Rcpp::NumericMatrix& process_normals,
    const Rcpp::NumericVector& observation_normals) {
  Rcpp::XPtr<liberation::ModelEngine> engine(engine_pointer);
  return liberation::kalman_simulate(
    *engine, data, theta, eta, sigma, process_normals, observation_normals);
}

// [[Rcpp::export(name = ".liberation_engine_derivative")]]
Rcpp::NumericVector liberation_engine_derivative(
    SEXP engine_pointer,
    const Rcpp::DataFrame& data,
    int row,
    int subject,
    double time,
    const Rcpp::NumericVector& state,
    const Rcpp::NumericVector& theta,
    const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma) {
  Rcpp::XPtr<liberation::ModelEngine> engine(engine_pointer);
  if (row < 1 || row > data.nrows()) Rcpp::stop("Derivative row is outside the dataset.");
  if (subject < 1 || subject > eta.nrow()) Rcpp::stop("Derivative subject is outside the ETA matrix.");
  if (state.size() != engine->n_state) Rcpp::stop("Derivative state has the wrong length.");
  liberation::Parameters parameters = liberation::evaluate_parameters(
    *engine, data, row - 1, subject - 1, theta, eta, sigma
  );
  const auto mapped = libertad::r_vector_map(state);
  return libertad::eigen_vector_to_r(liberation::evaluate_derivatives(
    *engine, data, row - 1, subject - 1, time, mapped, parameters, theta, eta, sigma
  ));
}

// [[Rcpp::export(name = ".liberation_matrix_exp")]]
Rcpp::NumericMatrix liberation_matrix_exp(const Rcpp::NumericMatrix& matrix,
                                           double dt = 1.0) {
  const auto mapped = libertad::r_matrix_map(matrix);
  return libertad::eigen_matrix_to_r(liberation::matrix_exp(mapped * dt));
}

// [[Rcpp::export(name = ".liberation_advan_matrix")]]
Rcpp::List liberation_advan_matrix(int advan, const Rcpp::List& parameters) {
  liberation::Parameters p;
  Rcpp::CharacterVector names = parameters.names();
  for (R_xlen_t i = 0; i < parameters.size(); ++i) {
    p[Rcpp::as<std::string>(names[i])] = Rcpp::as<double>(parameters[i]);
  }
  liberation::Topology topology = liberation::build_topology(advan, p);
  return Rcpp::List::create(
    Rcpp::Named("K") = libertad::eigen_matrix_to_r(topology.k),
    Rcpp::Named("states") = topology.state_names
  );
}

// [[Rcpp::export(name = ".liberation_prediction_tape_create")]]
SEXP liberation_prediction_tape_create(
    SEXP engine_pointer, const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta, const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma) {
  Rcpp::XPtr<liberation::ModelEngine> engine(engine_pointer);
  std::unique_ptr<liberation::PredictionTape> tape = liberation::record_prediction_tape(
    *engine, data, theta, eta, sigma);
  Rcpp::XPtr<liberation::PredictionTape> pointer(tape.release(), true);
  pointer.attr("class") = Rcpp::CharacterVector::create(
    "liberation_prediction_tape_ptr", "externalptr");
  pointer.attr("domain") = Rcpp::wrap(pointer->domain_names);
  pointer.attr("dynamic_columns") = Rcpp::wrap(pointer->dynamic_columns);
  pointer.attr("dynamic_parameters") =
    static_cast<double>(pointer->fun.size_dyn_ind());
  pointer.attr("propagation_kernel") = pointer->propagation_kernel;
  pointer.attr("operation_count") = static_cast<double>(pointer->operation_count);
  pointer.attr("variable_count") = static_cast<double>(pointer->variable_count);
  return pointer;
}

// [[Rcpp::export(name = ".liberation_prediction_tape_info")]]
Rcpp::List liberation_prediction_tape_info(SEXP tape_pointer) {
  Rcpp::XPtr<liberation::PredictionTape> tape(tape_pointer);
  const std::size_t taylor_bytes =
    tape->fun.size_var() * tape->fun.size_order() *
    std::max<std::size_t>(tape->fun.size_direction(), 1U) * sizeof(double);
  const std::size_t resident_proxy = tape->fun.size_op_seq() +
    tape->fun.size_random() + tape->fun.size_forward_bool() +
    tape->fun.size_forward_set() + taylor_bytes;
  return Rcpp::List::create(
    Rcpp::Named("operations") = static_cast<double>(tape->fun.size_op()),
    Rcpp::Named("operator_arguments") =
      static_cast<double>(tape->fun.size_op_arg()),
    Rcpp::Named("variables") = static_cast<double>(tape->fun.size_var()),
    Rcpp::Named("parameters") = static_cast<double>(tape->fun.size_par()),
    Rcpp::Named("dynamic_independent") =
      static_cast<double>(tape->fun.size_dyn_ind()),
    Rcpp::Named("dynamic_parameters") =
      static_cast<double>(tape->fun.size_dyn_par()),
    Rcpp::Named("dynamic_arguments") =
      static_cast<double>(tape->fun.size_dyn_arg()),
    Rcpp::Named("taylor_orders") =
      static_cast<double>(tape->fun.size_order()),
    Rcpp::Named("taylor_directions") =
      static_cast<double>(tape->fun.size_direction()),
    Rcpp::Named("operation_sequence_bytes") =
      static_cast<double>(tape->fun.size_op_seq()),
    Rcpp::Named("random_access_bytes") =
      static_cast<double>(tape->fun.size_random()),
    Rcpp::Named("forward_sparsity_bytes") = static_cast<double>(
      tape->fun.size_forward_bool() + tape->fun.size_forward_set()),
    Rcpp::Named("taylor_bytes_proxy") = static_cast<double>(taylor_bytes),
    Rcpp::Named("resident_bytes_proxy") =
      static_cast<double>(resident_proxy),
    Rcpp::Named("propagation_kernel") = tape->propagation_kernel,
    Rcpp::Named("derivative_strategy") = tape->derivative_strategy,
    Rcpp::Named("jacobian_nonzeros") =
      static_cast<double>(tape->jacobian_nonzeros)
  );
}

// [[Rcpp::export(name = ".liberation_prediction_tape_new_dynamic")]]
Rcpp::NumericVector liberation_prediction_tape_new_dynamic(
    SEXP tape_pointer, const Rcpp::DataFrame& data) {
  Rcpp::XPtr<liberation::PredictionTape> tape(tape_pointer);
  std::vector<double> values = liberation::prediction_dynamic_values(
    tape->dynamic_columns, data, tape->n_rows);
  tape->fun.new_dynamic(values);
  tape->dynamic_values = values;
  Rcpp::NumericVector result(values.begin(), values.end());
  result.attr("columns") = Rcpp::wrap(tape->dynamic_columns);
  return result;
}

// [[Rcpp::export(name = ".liberation_fo_tape_new_dynamic")]]
Rcpp::NumericVector liberation_fo_tape_new_dynamic(
    SEXP tape_pointer, const Rcpp::DataFrame& data) {
  Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointer);
  liberation::set_fo_dynamic(*tape, data);
  return Rcpp::wrap(tape->dynamic_values);
}

// [[Rcpp::export(name = ".liberation_prediction_tape_eval")]]
Rcpp::List liberation_prediction_tape_eval(
    SEXP tape_pointer, const Rcpp::NumericVector& point, bool jacobian = true) {
  Rcpp::XPtr<liberation::PredictionTape> tape(tape_pointer);
  std::vector<double> x = liberation::prediction_point(*tape, point);
  std::ostringstream messages;
  std::vector<double> value = tape->fun.Forward(0, x, messages);
  liberation::require_unchanged_path(tape->fun, "prediction evaluation");
  Rcpp::List result = Rcpp::List::create(Rcpp::Named("value") = Rcpp::wrap(value));
  result.attr("domain") = Rcpp::wrap(tape->domain_names);
  if (jacobian) {
    const std::size_t n = tape->domain_names.size();
    const std::size_t m = static_cast<std::size_t>(tape->n_rows);
    Rcpp::NumericMatrix derivative(m, n);
    std::size_t nonzeros = 0U;
    if (m * n >= 4096U && m >= 32U) {
      CppAD::vectorBool select_domain(n), select_range(m);
      for (std::size_t column = 0; column < n; ++column) select_domain[column] = true;
      for (std::size_t row = 0; row < m; ++row) select_range[row] = true;
      using SizeVector = CppAD::vector<std::size_t>;
      using BaseVector = CppAD::vector<double>;
      CppAD::sparse_rcv<SizeVector, BaseVector> sparse;
      BaseVector sparse_point(x.size());
      for (std::size_t index = 0; index < x.size(); ++index) sparse_point[index] = x[index];
      tape->fun.subgraph_jac_rev(
        select_domain, select_range, sparse_point, sparse);
      liberation::require_unchanged_path(
        tape->fun, "sparse prediction evaluation");
      for (std::size_t index = 0; index < sparse.nnz(); ++index) {
        derivative(sparse.row()[index], sparse.col()[index]) = sparse.val()[index];
      }
      nonzeros = sparse.nnz();
      tape->derivative_strategy = "subgraph-reverse";
    } else {
      constexpr std::size_t block_max = 16U;
      for (std::size_t first = 0; first < n; first += block_max) {
        const std::size_t directions = std::min(block_max, n - first);
        std::vector<double> seed(n * directions, 0.0);
        for (std::size_t direction = 0; direction < directions; ++direction) {
          seed[(first + direction) * directions + direction] = 1.0;
        }
        const std::vector<double> forward = directions == 1U ?
          tape->fun.Forward(1, seed) :
          tape->fun.Forward(1, directions, seed);
        for (std::size_t row = 0; row < m; ++row) {
          for (std::size_t direction = 0; direction < directions; ++direction) {
            const double current = forward[row * directions + direction];
            derivative(row, first + direction) = current;
            if (current != 0.0) ++nonzeros;
          }
        }
      }
      tape->derivative_strategy = n == 1U ? "forward" : "multi-forward";
    }
    tape->jacobian_nonzeros = nonzeros;
    derivative.attr("dimnames") = Rcpp::List::create(R_NilValue, Rcpp::wrap(tape->domain_names));
    result["jacobian"] = derivative;
  }
  result.attr("derivative_strategy") = tape->derivative_strategy;
  result.attr("jacobian_nonzeros") = static_cast<double>(tape->jacobian_nonzeros);
  return result;
}

// [[Rcpp::export(name = ".liberation_prediction_tape_eval_subset")]]
Rcpp::List liberation_prediction_tape_eval_subset(
    SEXP tape_pointer, const Rcpp::NumericVector& point,
    const Rcpp::IntegerVector& columns) {
  Rcpp::XPtr<liberation::PredictionTape> tape(tape_pointer);
  std::vector<double> x = liberation::prediction_point(*tape, point);
  std::ostringstream messages;
  const std::vector<double> value = tape->fun.Forward(0, x, messages);
  liberation::require_unchanged_path(tape->fun, "prediction subset evaluation");
  const std::size_t domain = tape->domain_names.size();
  const std::size_t range = static_cast<std::size_t>(tape->n_rows);
  Rcpp::NumericMatrix derivative(range, columns.size());
  Rcpp::CharacterVector names(columns.size());
  std::vector<std::size_t> selected_columns(static_cast<std::size_t>(columns.size()));
  for (R_xlen_t selected = 0; selected < columns.size(); ++selected) {
    const int column = columns[selected] - 1;
    if (column < 0 || static_cast<std::size_t>(column) >= domain) {
      Rcpp::stop("Prediction derivative column is outside the tape domain.");
    }
    selected_columns[static_cast<std::size_t>(selected)] =
      static_cast<std::size_t>(column);
    names[selected] = tape->domain_names[static_cast<std::size_t>(column)];
  }
  constexpr std::size_t block_max = 16U;
  for (std::size_t first = 0; first < selected_columns.size(); first += block_max) {
    const std::size_t directions = std::min(block_max, selected_columns.size() - first);
    std::vector<double> seed(domain * directions, 0.0);
    for (std::size_t direction = 0; direction < directions; ++direction) {
      seed[selected_columns[first + direction] * directions + direction] = 1.0;
    }
    const std::vector<double> forward = directions == 1U ?
      tape->fun.Forward(1, seed) :
      tape->fun.Forward(1, directions, seed);
    for (std::size_t row = 0; row < range; ++row) {
      for (std::size_t direction = 0; direction < directions; ++direction) {
        derivative(static_cast<int>(row), static_cast<int>(first + direction)) =
          forward[row * directions + direction];
      }
    }
  }
  tape->derivative_strategy = selected_columns.size() <= 1U ?
    "forward-subset" : "multi-forward-subset";
  derivative.attr("dimnames") = Rcpp::List::create(R_NilValue, names);
  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("value") = Rcpp::wrap(value),
    Rcpp::Named("jacobian") = derivative
  );
  result.attr("domain") = names;
  return result;
}

// [[Rcpp::export(name = ".liberation_matrix_exp_pade")]]
Rcpp::NumericMatrix liberation_matrix_exp_pade(const Rcpp::NumericMatrix& matrix,
                                                double dt = 1.0) {
  const auto mapped = libertad::r_matrix_map(matrix);
  return libertad::eigen_matrix_to_r(
    liberation::matrix_exp_pade(Eigen::MatrixXd(mapped * dt)));
}

// [[Rcpp::export(name = ".liberation_fo_tape_create")]]
SEXP liberation_fo_tape_create(
    SEXP engine_pointer, SEXP prediction_tape_pointer,
    const Rcpp::DataFrame& data, const Rcpp::NumericVector& theta,
    const Rcpp::NumericVector& sigma, const Rcpp::NumericVector& omega) {
  Rcpp::XPtr<liberation::ModelEngine> engine(engine_pointer);
  Rcpp::XPtr<liberation::PredictionTape> prediction_tape(prediction_tape_pointer);
  std::unique_ptr<liberation::ObjectiveTape> tape = liberation::record_fo_tape(
    *engine, *prediction_tape, data, theta, sigma, omega);
  Rcpp::XPtr<liberation::ObjectiveTape> pointer(tape.release(), true);
  pointer.attr("class") = Rcpp::CharacterVector::create(
    "liberation_fo_tape_ptr", "liberation_objective_tape_ptr", "externalptr");
  pointer.attr("domain") = Rcpp::wrap(pointer->domain_names);
  return pointer;
}

// [[Rcpp::export(name = ".liberation_curvature_tape_create")]]
SEXP liberation_curvature_tape_create(
    SEXP engine_pointer, SEXP prediction_tape_pointer,
    SEXP objective_tape_pointer, const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta, const Rcpp::NumericVector& eta,
    const Rcpp::NumericVector& sigma, const Rcpp::NumericVector& omega,
    const std::string& approximation) {
  Rcpp::XPtr<liberation::ModelEngine> engine(engine_pointer);
  Rcpp::XPtr<liberation::PredictionTape> prediction_tape(prediction_tape_pointer);
  Rcpp::XPtr<liberation::ObjectiveTape> objective_tape(objective_tape_pointer);
  std::unique_ptr<liberation::ObjectiveTape> tape = liberation::record_curvature_tape(
    *engine, *prediction_tape, *objective_tape, data,
    theta, eta, sigma, omega, approximation);
  Rcpp::XPtr<liberation::ObjectiveTape> pointer(tape.release(), true);
  pointer.attr("class") = Rcpp::CharacterVector::create(
    "liberation_curvature_tape_ptr", "liberation_objective_tape_ptr", "externalptr");
  pointer.attr("domain") = Rcpp::wrap(pointer->domain_names);
  return pointer;
}

// [[Rcpp::export(name = ".liberation_objective_tape_create")]]
SEXP liberation_objective_tape_create(
    SEXP engine_pointer, const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta, const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma, const Rcpp::NumericVector& omega,
    bool interaction = true) {
  Rcpp::XPtr<liberation::ModelEngine> engine(engine_pointer);
  std::unique_ptr<liberation::ObjectiveTape> tape = liberation::record_objective_tape(
    *engine, data, theta, eta, sigma, omega, interaction);
  Rcpp::XPtr<liberation::ObjectiveTape> pointer(tape.release(), true);
  pointer.attr("class") = Rcpp::CharacterVector::create(
    "liberation_objective_tape_ptr", "externalptr");
  pointer.attr("domain") = Rcpp::wrap(pointer->domain_names);
  return pointer;
}

// [[Rcpp::export(name = ".liberation_objective_tape_eval")]]
Rcpp::List liberation_objective_tape_eval(
    SEXP tape_pointer, const Rcpp::NumericVector& point,
    bool gradient = true, bool hessian = false) {
  Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointer);
  if (point.size() != static_cast<R_xlen_t>(tape->domain_names.size())) {
    Rcpp::stop("Objective tape point has the wrong length.");
  }
  std::vector<double> x = Rcpp::as<std::vector<double>>(point);
  std::ostringstream messages;
  std::vector<double> value = tape->fun.Forward(0, x, messages);
  liberation::require_unchanged_path(tape->fun, "objective evaluation");
  Rcpp::List result = Rcpp::List::create(Rcpp::Named("value") = value[0]);
  if (gradient || hessian) {
    std::vector<double> weight(1, 1.0);
    std::vector<double> derivative = tape->fun.Reverse(1, weight);
    Rcpp::NumericVector output(derivative.begin(), derivative.end());
    output.attr("names") = Rcpp::wrap(tape->domain_names);
    result["gradient"] = output;
  }
  if (hessian) {
    const std::size_t n = tape->domain_names.size();
    Rcpp::NumericMatrix output(n, n);
    libertad::analyse_hessian_sparsity(
      tape->fun, tape->hessian_cache);
    if (tape->hessian_cache.use_sparse) {
      const std::vector<double> values = libertad::sparse_hessian(
        tape->fun, x, tape->hessian_cache);
      liberation::require_unchanged_path(
        tape->fun, "sparse objective Hessian evaluation");
      for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t column = 0; column < n; ++column) {
          output(row, column) = values[row * n + column];
        }
      }
    } else {
      std::vector<double> direction(n, 0.0);
      std::vector<double> weight(1, 1.0);
      for (std::size_t column = 0; column < n; ++column) {
        direction[column] = 1.0;
        tape->fun.Forward(1, direction, messages);
        direction[column] = 0.0;
        std::vector<double> reverse = tape->fun.Reverse(2, weight);
        for (std::size_t row = 0; row < n; ++row) {
          output(row, column) = reverse[row * 2 + 1];
        }
      }
    }
    output.attr("dimnames") = Rcpp::List::create(
      Rcpp::wrap(tape->domain_names), Rcpp::wrap(tape->domain_names));
    result["hessian"] = output;
  }
  result.attr("domain") = Rcpp::wrap(tape->domain_names);
  return result;
}

// [[Rcpp::export(name = ".liberation_objective_tape_info")]]
Rcpp::List liberation_objective_tape_info(SEXP tape_pointer) {
  Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointer);
  const std::size_t taylor_bytes =
    tape->fun.size_var() * tape->fun.size_order() *
    std::max<std::size_t>(tape->fun.size_direction(), 1U) * sizeof(double);
  return Rcpp::List::create(
    Rcpp::Named("operations") = static_cast<double>(tape->fun.size_op()),
    Rcpp::Named("operator_arguments") =
      static_cast<double>(tape->fun.size_op_arg()),
    Rcpp::Named("variables") = static_cast<double>(tape->fun.size_var()),
    Rcpp::Named("parameters") = static_cast<double>(tape->fun.size_par()),
    Rcpp::Named("dynamic_independent") =
      static_cast<double>(tape->fun.size_dyn_ind()),
    Rcpp::Named("dynamic_parameters") =
      static_cast<double>(tape->fun.size_dyn_par()),
    Rcpp::Named("operation_sequence_bytes") =
      static_cast<double>(tape->fun.size_op_seq()),
    Rcpp::Named("random_access_bytes") =
      static_cast<double>(tape->fun.size_random()),
    Rcpp::Named("forward_sparsity_bytes") = static_cast<double>(
      tape->fun.size_forward_bool() + tape->fun.size_forward_set()),
    Rcpp::Named("taylor_bytes_proxy") = static_cast<double>(taylor_bytes),
    Rcpp::Named("hessian_strategy") = tape->hessian_cache.strategy,
    Rcpp::Named("hessian_nonzeros") =
      static_cast<double>(tape->hessian_cache.nonzeros),
    Rcpp::Named("hessian_density") = tape->hessian_cache.density,
    Rcpp::Named("hessian_sweeps") =
      static_cast<double>(tape->hessian_cache.sweeps),
    Rcpp::Named("resident_bytes_proxy") = static_cast<double>(
      tape->fun.size_op_seq() + tape->fun.size_random() +
      tape->fun.size_forward_bool() + tape->fun.size_forward_set() +
      taylor_bytes)
  );
}

// [[Rcpp::export(name = ".liberation_hmc_target_eval")]]
Rcpp::List liberation_hmc_target_eval(
    SEXP tape_pointer, const Rcpp::NumericVector& q,
    const Rcpp::List& config) {
  Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointer);
  return liberation::native_hmc_target_eval(*tape, q, config);
}

// [[Rcpp::export(name = ".liberation_hmc_sample")]]
Rcpp::List liberation_hmc_sample(
    SEXP tape_pointer, const Rcpp::List& config, const std::string& method,
    int n_warmup, int n_sample, int n_thin, int n_chains, double seed,
    double step_size, double target_acceptance, bool adapt_mass,
    int n_leapfrog, int max_depth, double divergence_threshold,
    int print_every) {
  if (!std::isfinite(seed) || seed < 0.0) {
    Rcpp::stop("Native HMC seed must be a non-negative finite number.");
  }
  Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointer);
  return liberation::native_hmc_sample(
    *tape, config, method, n_warmup, n_sample, n_thin, n_chains,
    static_cast<std::uint64_t>(seed), step_size, target_acceptance,
    adapt_mass, n_leapfrog, max_depth, divergence_threshold, print_every
  );
}

// [[Rcpp::export(name = ".liberation_objective_tape_eta_values")]]
Rcpp::NumericVector liberation_objective_tape_eta_values(
    SEXP tape_pointer, const Rcpp::NumericVector& point,
    const Rcpp::IntegerVector& eta_positions,
    const Rcpp::NumericMatrix& eta) {
  Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointer);
  const std::size_t domain = tape->domain_names.size();
  if (point.size() != static_cast<R_xlen_t>(domain)) {
    Rcpp::stop("Objective tape point has the wrong length.");
  }
  if (eta.ncol() != eta_positions.size()) {
    Rcpp::stop("ETA samples have the wrong number of columns.");
  }
  std::vector<std::size_t> positions;
  positions.reserve(static_cast<std::size_t>(eta_positions.size()));
  for (int value : eta_positions) {
    if (value < 1 || static_cast<std::size_t>(value) > domain) {
      Rcpp::stop("ETA position is outside the objective tape domain.");
    }
    positions.push_back(static_cast<std::size_t>(value - 1));
  }
  std::vector<double> x = Rcpp::as<std::vector<double>>(point);
  Rcpp::NumericVector values(eta.nrow());
  std::ostringstream messages;
  for (int sample = 0; sample < eta.nrow(); ++sample) {
    for (int column = 0; column < eta.ncol(); ++column) {
      x[positions[static_cast<std::size_t>(column)]] = eta(sample, column);
    }
    const std::vector<double> value = tape->fun.Forward(0, x, messages);
    liberation::require_unchanged_path(tape->fun, "objective ETA batch");
    values[sample] = value.empty() ? NA_REAL : value[0];
    if ((sample + 1) % 256 == 0) Rcpp::checkUserInterrupt();
  }
  return values;
}

// [[Rcpp::export(name = ".liberation_objective_tape_collection_values")]]
Rcpp::NumericVector liberation_objective_tape_collection_values(
    const Rcpp::List& tape_pointers, const Rcpp::NumericMatrix& points) {
  if (points.nrow() != tape_pointers.size()) {
    Rcpp::stop("Objective point rows must match the number of tapes.");
  }
  Rcpp::NumericVector values(points.nrow());
  for (int row = 0; row < points.nrow(); ++row) {
    Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointers[row]);
    if (points.ncol() != static_cast<int>(tape->domain_names.size())) {
      Rcpp::stop("An objective point has the wrong length.");
    }
    std::vector<double> point(static_cast<std::size_t>(points.ncol()));
    for (int column = 0; column < points.ncol(); ++column) {
      point[static_cast<std::size_t>(column)] = points(row, column);
    }
    std::ostringstream messages;
    const std::vector<double> value = tape->fun.Forward(0, point, messages);
    liberation::require_unchanged_path(tape->fun, "objective collection");
    values[row] = value.empty() ? NA_REAL : value[0];
    if ((row + 1) % 256 == 0) Rcpp::checkUserInterrupt();
  }
  return values;
}

// [[Rcpp::export(name = ".liberation_objective_tape_collection_gradients")]]
Rcpp::NumericMatrix liberation_objective_tape_collection_gradients(
    const Rcpp::List& tape_pointers, const Rcpp::NumericMatrix& points) {
  if (points.nrow() != tape_pointers.size()) {
    Rcpp::stop("Objective point rows must match the number of tapes.");
  }
  if (!points.nrow()) return Rcpp::NumericMatrix(0, points.ncol());
  Rcpp::NumericMatrix gradients(points.nrow(), points.ncol());
  const std::vector<double> weight(1, 1.0);
  for (int row = 0; row < points.nrow(); ++row) {
    Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointers[row]);
    if (points.ncol() != static_cast<int>(tape->domain_names.size())) {
      Rcpp::stop("An objective point has the wrong length.");
    }
    std::vector<double> point(static_cast<std::size_t>(points.ncol()));
    for (int column = 0; column < points.ncol(); ++column) {
      point[static_cast<std::size_t>(column)] = points(row, column);
    }
    std::ostringstream messages;
    tape->fun.Forward(0, point, messages);
    liberation::require_unchanged_path(tape->fun, "objective gradient collection");
    const std::vector<double> derivative = tape->fun.Reverse(1, weight);
    for (int column = 0; column < points.ncol(); ++column) {
      gradients(row, column) = derivative[static_cast<std::size_t>(column)];
    }
    if ((row + 1) % 256 == 0) Rcpp::checkUserInterrupt();
  }
  return gradients;
}

// [[Rcpp::export(name = ".liberation_objective_tape_hessian_subset")]]
Rcpp::NumericMatrix liberation_objective_tape_hessian_subset(
    SEXP tape_pointer, const Rcpp::NumericVector& point,
    const Rcpp::IntegerVector& row_positions,
    const Rcpp::IntegerVector& column_positions) {
  Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointer);
  const std::size_t domain = tape->domain_names.size();
  if (point.size() != static_cast<R_xlen_t>(domain)) {
    Rcpp::stop("Objective tape point has the wrong length.");
  }
  auto positions = [domain](const Rcpp::IntegerVector& source) {
    std::vector<std::size_t> result;
    result.reserve(static_cast<std::size_t>(source.size()));
    for (int value : source) {
      if (value < 1 || static_cast<std::size_t>(value) > domain) {
        Rcpp::stop("Hessian position is outside the objective tape domain.");
      }
      result.push_back(static_cast<std::size_t>(value - 1));
    }
    return result;
  };
  const std::vector<std::size_t> rows = positions(row_positions);
  const std::vector<std::size_t> columns = positions(column_positions);
  Rcpp::NumericMatrix result(rows.size(), columns.size());
  std::vector<double> x = Rcpp::as<std::vector<double>>(point);
  std::ostringstream messages;
  tape->fun.Forward(0, x, messages);
  liberation::require_unchanged_path(tape->fun, "objective Hessian subset");
  const std::vector<double> weight(1, 1.0);
  std::vector<double> direction(domain, 0.0);
  for (std::size_t column = 0; column < columns.size(); ++column) {
    direction[columns[column]] = 1.0;
    tape->fun.Forward(1, direction, messages);
    direction[columns[column]] = 0.0;
    const std::vector<double> reverse = tape->fun.Reverse(2, weight);
    for (std::size_t row = 0; row < rows.size(); ++row) {
      result(static_cast<int>(row), static_cast<int>(column)) =
        reverse[rows[row] * 2U + 1U];
    }
  }
  return result;
}

// [[Rcpp::export(name = ".liberation_nested_population_gradient")]]
Rcpp::List liberation_nested_population_gradient(
    const Rcpp::List& objective_tapes, const Rcpp::List& curvature_tapes,
    const Rcpp::NumericMatrix& points,
    const Rcpp::IntegerVector& eta_positions,
    const Rcpp::IntegerVector& population_positions,
    const Rcpp::NumericMatrix& transform) {
  const int subjects = points.nrow();
  const int n_eta = eta_positions.size();
  const int n_population = population_positions.size();
  const int n_outer = transform.ncol();
  if (objective_tapes.size() != subjects || curvature_tapes.size() != subjects ||
      transform.nrow() != n_population) {
    Rcpp::stop("Nested-gradient batch dimensions are inconsistent.");
  }
  Rcpp::NumericMatrix subject_gradients(subjects, n_outer);
  Rcpp::NumericVector jitters(subjects);
  const std::vector<double> weight(1, 1.0);
  for (int subject = 0; subject < subjects; ++subject) {
    Rcpp::XPtr<liberation::ObjectiveTape> objective(objective_tapes[subject]);
    Rcpp::XPtr<liberation::ObjectiveTape> curvature(curvature_tapes[subject]);
    const std::size_t domain = objective->domain_names.size();
    if (points.ncol() != static_cast<int>(domain) ||
        curvature->domain_names.size() != domain) {
      Rcpp::stop("A nested-gradient objective point has the wrong length.");
    }
    std::vector<std::size_t> eta, population;
    for (int value : eta_positions) {
      if (value < 1 || static_cast<std::size_t>(value) > domain) {
        Rcpp::stop("ETA position is outside a nested-gradient tape domain.");
      }
      eta.push_back(static_cast<std::size_t>(value - 1));
    }
    for (int value : population_positions) {
      if (value < 1 || static_cast<std::size_t>(value) > domain) {
        Rcpp::stop("Population position is outside a nested-gradient tape domain.");
      }
      population.push_back(static_cast<std::size_t>(value - 1));
    }
    std::vector<double> point(domain);
    for (int column = 0; column < points.ncol(); ++column) {
      point[static_cast<std::size_t>(column)] = points(subject, column);
    }
    std::ostringstream messages;
    objective->fun.Forward(0, point, messages);
    const std::vector<double> objective_derivative = objective->fun.Reverse(1, weight);
    Eigen::MatrixXd mixed(n_eta, n_eta + n_population);
    std::vector<double> direction(domain, 0.0);
    for (int column = 0; column < n_eta + n_population; ++column) {
      const std::size_t position = column < n_eta ?
        eta[static_cast<std::size_t>(column)] :
        population[static_cast<std::size_t>(column - n_eta)];
      direction[position] = 1.0;
      objective->fun.Forward(1, direction, messages);
      direction[position] = 0.0;
      const std::vector<double> reverse = objective->fun.Reverse(2, weight);
      for (int row = 0; row < n_eta; ++row) {
        mixed(row, column) = reverse[eta[static_cast<std::size_t>(row)] * 2U + 1U];
      }
    }
    Eigen::MatrixXd eta_hessian;
    if (n_eta) {
      eta_hessian = 0.5 *
        (mixed.leftCols(n_eta) + mixed.leftCols(n_eta).transpose()).eval();
    } else {
      eta_hessian = Eigen::MatrixXd::Zero(0, 0);
    }
    double jitter = 0.0;
    if (n_eta) {
      auto eigen = liberation::detail::self_adjoint_eigen(eta_hessian, false);
      if (eigen.info != Eigen::Success) {
        Rcpp::stop("Conditional ETA curvature eigen decomposition failed.");
      }
      const double largest = std::max(eigen.values.cwiseAbs().maxCoeff(), 1.0);
      jitter = std::max(0.0, largest * 1e-9 - eigen.values.minCoeff());
      if (jitter > largest * 1e-2) {
        Rcpp::stop("Conditional ETA curvature is not sufficiently positive definite.");
      }
      eta_hessian.diagonal().array() += jitter;
    }
    jitters[subject] = jitter;
    Eigen::MatrixXd mapped_transform(n_population, n_outer);
    for (int row = 0; row < n_population; ++row) {
      for (int column = 0; column < n_outer; ++column) {
        mapped_transform(row, column) = transform(row, column);
      }
    }
    Eigen::MatrixXd sensitivity;
    if (n_eta) {
      sensitivity = -eta_hessian.ldlt().solve(
        mixed.rightCols(n_population) * mapped_transform);
    } else {
      sensitivity = Eigen::MatrixXd::Zero(0, n_outer);
    }
    curvature->fun.Forward(0, point, messages);
    const std::vector<double> curvature_derivative = curvature->fun.Reverse(1, weight);
    for (int outer = 0; outer < n_outer; ++outer) {
      double derivative = 0.0;
      for (int native = 0; native < n_population; ++native) {
        const double chain = transform(native, outer);
        derivative += (objective_derivative[population[static_cast<std::size_t>(native)]] +
          curvature_derivative[population[static_cast<std::size_t>(native)]]) * chain;
      }
      for (int effect = 0; effect < n_eta; ++effect) {
        derivative += curvature_derivative[eta[static_cast<std::size_t>(effect)]] *
          sensitivity(effect, outer);
      }
      subject_gradients(subject, outer) = derivative;
    }
    if ((subject + 1) % 64 == 0) Rcpp::checkUserInterrupt();
  }
  Rcpp::NumericVector gradient(n_outer);
  for (int outer = 0; outer < n_outer; ++outer) {
    for (int subject = 0; subject < subjects; ++subject) {
      gradient[outer] += subject_gradients(subject, outer);
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("gradient") = gradient,
    Rcpp::Named("subject_gradients") = subject_gradients,
    Rcpp::Named("eta_jitter") = jitters);
}

// [[Rcpp::export(name = ".liberation_objective_tape_eta_mode")]]
Rcpp::List liberation_objective_tape_eta_mode(
    SEXP tape_pointer, const Rcpp::NumericVector& point,
    const Rcpp::IntegerVector& eta_positions,
    const Rcpp::NumericVector& start, int maxit = 100,
    double tolerance = 1e-7, bool exact_hessian = true) {
  Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointer);
  if (point.size() != static_cast<R_xlen_t>(tape->domain_names.size())) {
    Rcpp::stop("Objective tape point has the wrong length.");
  }
  std::vector<std::size_t> positions;
  positions.reserve(static_cast<std::size_t>(eta_positions.size()));
  for (int value : eta_positions) {
    if (value < 1 || static_cast<std::size_t>(value) > tape->domain_names.size()) {
      Rcpp::stop("ETA position is outside the objective tape domain.");
    }
    positions.push_back(static_cast<std::size_t>(value - 1));
  }
  return liberation::objective_eta_mode(
    *tape, Rcpp::as<std::vector<double>>(point), positions, start,
    maxit, tolerance, exact_hessian);
}

// [[Rcpp::export(name = ".liberation_objective_tape_eta_modes")]]
Rcpp::List liberation_objective_tape_eta_modes(
    const Rcpp::List& tape_pointers, const Rcpp::NumericMatrix& points,
    const Rcpp::IntegerVector& eta_positions,
    const Rcpp::NumericMatrix& starts, int maxit = 100,
    double tolerance = 1e-7, bool exact_hessian = true) {
  if (points.nrow() != tape_pointers.size() || starts.nrow() != points.nrow()) {
    Rcpp::stop("ETA-mode rows must match the number of objective tapes.");
  }
  if (starts.ncol() != eta_positions.size()) {
    Rcpp::stop("ETA starting values have the wrong number of columns.");
  }
  Rcpp::List result(points.nrow());
  for (int row = 0; row < points.nrow(); ++row) {
    Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointers[row]);
    const std::size_t domain = tape->domain_names.size();
    if (points.ncol() != static_cast<int>(domain)) {
      Rcpp::stop("An ETA-mode objective point has the wrong length.");
    }
    std::vector<std::size_t> positions;
    positions.reserve(static_cast<std::size_t>(eta_positions.size()));
    for (int value : eta_positions) {
      if (value < 1 || static_cast<std::size_t>(value) > domain) {
        Rcpp::stop("ETA position is outside an objective tape domain.");
      }
      positions.push_back(static_cast<std::size_t>(value - 1));
    }
    Rcpp::NumericVector start(starts.ncol());
    std::vector<double> point(static_cast<std::size_t>(points.ncol()));
    for (int column = 0; column < points.ncol(); ++column) {
      point[static_cast<std::size_t>(column)] = points(row, column);
    }
    for (int column = 0; column < starts.ncol(); ++column) {
      start[column] = starts(row, column);
    }
    result[row] = liberation::objective_eta_mode(
      *tape, std::move(point), positions, start, maxit, tolerance,
      exact_hessian);
    if ((row + 1) % 64 == 0) Rcpp::checkUserInterrupt();
  }
  return result;
}

// [[Rcpp::export(name = ".liberation_objective_tape_point_gradients")]]
Rcpp::List liberation_objective_tape_point_gradients(
    SEXP tape_pointer, const Rcpp::NumericMatrix& points) {
  Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointer);
  const std::size_t domain = tape->domain_names.size();
  if (points.ncol() != static_cast<int>(domain)) {
    Rcpp::stop("Objective sample points have the wrong number of columns.");
  }
  Rcpp::NumericVector values(points.nrow());
  Rcpp::NumericMatrix gradients(points.nrow(), points.ncol());
  const std::vector<double> weight(1, 1.0);
  std::ostringstream messages;
  for (int row = 0; row < points.nrow(); ++row) {
    std::vector<double> point(domain);
    for (int column = 0; column < points.ncol(); ++column) {
      point[static_cast<std::size_t>(column)] = points(row, column);
    }
    const std::vector<double> value = tape->fun.Forward(0, point, messages);
    values[row] = value.empty() ? NA_REAL : value[0];
    const std::vector<double> derivative = tape->fun.Reverse(1, weight);
    for (int column = 0; column < points.ncol(); ++column) {
      gradients(row, column) = derivative[static_cast<std::size_t>(column)];
    }
    if ((row + 1) % 256 == 0) Rcpp::checkUserInterrupt();
  }
  gradients.attr("dimnames") = Rcpp::List::create(
    R_NilValue, Rcpp::wrap(tape->domain_names));
  return Rcpp::List::create(
    Rcpp::Named("value") = values,
    Rcpp::Named("gradient") = gradients);
}

// [[Rcpp::export(name = ".liberation_objective_tape_eta_metropolis")]]
Rcpp::List liberation_objective_tape_eta_metropolis(
    const Rcpp::List& tape_pointers, const Rcpp::NumericMatrix& points,
    const Rcpp::IntegerVector& eta_positions,
    const Rcpp::NumericMatrix& current_eta,
    const Rcpp::List& proposal_roots, const Rcpp::NumericMatrix& normals,
    const Rcpp::NumericVector& log_uniforms, int mcmc_steps,
    double step_scale = 0.5) {
  const int subjects = points.nrow();
  const int dimension = eta_positions.size();
  if (subjects != tape_pointers.size() || current_eta.nrow() != subjects ||
      current_eta.ncol() != dimension || proposal_roots.size() != subjects ||
      mcmc_steps < 1 || normals.nrow() != subjects * mcmc_steps ||
      normals.ncol() != dimension || log_uniforms.size() != normals.nrow() ||
      !std::isfinite(step_scale) || step_scale <= 0.0) {
    Rcpp::stop("Batched ETA Metropolis inputs are inconsistent.");
  }
  Rcpp::NumericMatrix eta = Rcpp::clone(current_eta);
  Rcpp::NumericVector values(subjects);
  int accepted = 0;
  for (int subject = 0; subject < subjects; ++subject) {
    Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointers[subject]);
    const std::size_t domain = tape->domain_names.size();
    if (points.ncol() != static_cast<int>(domain)) {
      Rcpp::stop("A batched Metropolis objective point has the wrong length.");
    }
    std::vector<std::size_t> positions;
    positions.reserve(static_cast<std::size_t>(dimension));
    for (int value : eta_positions) {
      if (value < 1 || static_cast<std::size_t>(value) > domain) {
        Rcpp::stop("ETA position is outside a Metropolis objective tape domain.");
      }
      positions.push_back(static_cast<std::size_t>(value - 1));
    }
    Rcpp::NumericMatrix root = proposal_roots[subject];
    if (root.nrow() != dimension || root.ncol() != dimension) {
      Rcpp::stop("A Metropolis proposal root has the wrong dimensions.");
    }
    std::vector<double> point(domain);
    for (int column = 0; column < points.ncol(); ++column) {
      point[static_cast<std::size_t>(column)] = points(subject, column);
    }
    for (int column = 0; column < dimension; ++column) {
      point[positions[static_cast<std::size_t>(column)]] = eta(subject, column);
    }
    liberation::EtaEvaluation current = liberation::objective_eta_evaluate(
      *tape, point, positions, false);
    if (!current.finite) Rcpp::stop("Current ETA objective is not finite.");
    for (int step = 0; step < mcmc_steps; ++step) {
      const int draw = subject * mcmc_steps + step;
      std::vector<double> candidate_point = point;
      Rcpp::NumericVector candidate_eta(dimension);
      for (int row = 0; row < dimension; ++row) {
        double increment = 0.0;
        for (int column = 0; column < dimension; ++column) {
          increment += root(row, column) * normals(draw, column);
        }
        candidate_eta[row] = eta(subject, row) + step_scale * increment;
        candidate_point[positions[static_cast<std::size_t>(row)]] = candidate_eta[row];
      }
      liberation::EtaEvaluation candidate = liberation::objective_eta_evaluate(
        *tape, candidate_point, positions, false);
      if (candidate.finite && log_uniforms[draw] <
          -0.5 * (candidate.value - current.value)) {
        for (int row = 0; row < dimension; ++row) {
          eta(subject, row) = candidate_eta[row];
        }
        point.swap(candidate_point);
        current = std::move(candidate);
        ++accepted;
      }
    }
    values[subject] = current.value;
    if ((subject + 1) % 64 == 0) Rcpp::checkUserInterrupt();
  }
  return Rcpp::List::create(
    Rcpp::Named("eta") = eta,
    Rcpp::Named("value") = values,
    Rcpp::Named("accepted") = accepted,
    Rcpp::Named("attempted") = subjects * mcmc_steps);
}

// [[Rcpp::export(name = ".liberation_mixture_component_nll")]]
Rcpp::NumericMatrix liberation_mixture_component_nll(
    SEXP engine_pointer, const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta, const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma) {
  Rcpp::XPtr<liberation::ModelEngine> engine(engine_pointer);
  if (engine->mixture_probabilities.empty()) {
    Rcpp::stop("The model does not define a finite mixture.");
  }
  Rcpp::IntegerVector subject_index = data[".ID_INDEX"];
  int n_subjects = 0;
  for (int value : subject_index) n_subjects = std::max(n_subjects, value);
  std::vector<double> theta_values = Rcpp::as<std::vector<double>>(theta);
  std::vector<double> eta_values;
  eta_values.reserve(static_cast<std::size_t>(eta.size()));
  for (int row = 0; row < eta.nrow(); ++row) {
    for (int column = 0; column < eta.ncol(); ++column) eta_values.push_back(eta(row, column));
  }
  std::vector<double> sigma_values = Rcpp::as<std::vector<double>>(sigma);
  Rcpp::NumericMatrix result(n_subjects, engine->mixture_probabilities.size());
  for (std::size_t component = 0; component < engine->mixture_probabilities.size(); ++component) {
    std::vector<int> assignment(static_cast<std::size_t>(n_subjects),
                                static_cast<int>(component + 1));
    std::vector<double> prediction = liberation::simulate_analytical_t(
      *engine, data, theta_values, eta_values, sigma_values, assignment);
    std::vector<double> nll = liberation::residual_subject_nll_t(
      *engine, data, prediction, theta_values, eta_values, sigma_values,
      assignment);
    for (int subject = 0; subject < n_subjects; ++subject) {
      result(subject, static_cast<int>(component)) = nll[static_cast<std::size_t>(subject)];
    }
  }
  return result;
}
