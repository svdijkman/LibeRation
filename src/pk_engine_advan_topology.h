// Internal implementation fragment included inside namespace liberation by
// pk_engine.cpp. It intentionally has no include guard or namespace declaration.
// Contains ADVAN topology, matrix-exponential, and steady-state construction.

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

