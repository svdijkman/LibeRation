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
    return "dae-advan" + std::to_string(engine.advan) +
      (implicit_ode_advan(engine.advan) ? "-implicit-newton" : "-rk45-newton");
  }
  if (engine.is_ode()) return ode_kernel_name(engine.advan);
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
Scalar evaluate_post_prediction_t(
    const ModelEngine& engine, const Rcpp::DataFrame& data,
    int row, int subject, const Scalar& time, const VectorT<Scalar>& state,
    const std::vector<Scalar>& theta, const std::vector<Scalar>& eta,
    int eta_columns, const std::vector<Scalar>& sigma, int mixture_number,
    const Scalar& advan_prediction, ParametersT<Scalar>& parameters,
    const DynamicDataT<Scalar>* dynamic_data = nullptr) {
  if (!engine.post_pred) return advan_prediction;
  std::vector<Scalar> inputs(
    engine.post_pred->input_names.size(), Scalar(0.0));
  for (std::size_t i = 0; i < engine.post_pred->input_names.size(); ++i) {
    const std::string& name = engine.post_pred->input_names[i];
    int index = indexed_name(name, "THETA_");
    if (index >= 0) {
      inputs[i] = theta.at(static_cast<std::size_t>(index));
      continue;
    }
    index = indexed_name(name, "ETA_");
    if (index >= 0) {
      const int column = eta_column(engine, data, row, index, eta_columns);
      inputs[i] = eta.at(
        static_cast<std::size_t>(subject * eta_columns + column));
      continue;
    }
    index = indexed_name(name, "SIGMA_");
    if (index >= 0) {
      inputs[i] = sigma.at(static_cast<std::size_t>(index));
      continue;
    }
    index = indexed_name(name, "A_");
    if (index >= 0) {
      if (index >= state.size()) {
        throw std::out_of_range("$PRED A() index exceeds state dimension.");
      }
      inputs[i] = state[index];
      continue;
    }
    if (name == "F_ADVAN") { inputs[i] = advan_prediction; continue; }
    if (name == "T" || name == "TIME") { inputs[i] = time; continue; }
    if (name == "MIXNUM") { inputs[i] = Scalar(mixture_number); continue; }
    const auto assigned = parameters.find(name);
    if (assigned != parameters.end()) {
      inputs[i] = assigned->second;
      continue;
    }
    inputs[i] = dynamic_row_value(data, name, row, dynamic_data);
    if (!std::isfinite(scalar_value(inputs[i]))) {
      throw std::domain_error(
        "Post-ADVAN $PRED input '" + name + "' is non-finite at row " +
        std::to_string(row + 1) + ".");
    }
  }
  const std::vector<Scalar> output =
    engine.post_pred->eval_outputs(inputs, engine.post_all_outputs);
  for (std::size_t i = 0; i < output.size(); ++i) {
    parameters[engine.post_pred->output_names[i]] = output[i];
  }
  const auto prediction = parameters.find("F");
  if (prediction == parameters.end() ||
      !std::isfinite(scalar_value(prediction->second))) {
    throw std::domain_error("Post-ADVAN $PRED did not produce a finite F.");
  }
  return prediction->second;
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
  struct Jump {
    Scalar time;
    VectorT<Scalar> before;
    VectorT<Scalar> after;
  };

  std::vector<Scalar> time;
  std::vector<VectorT<Scalar>> state;
  std::vector<double> baseline;
  std::vector<Jump> jumps;

  void reset(const Scalar& at, const VectorT<Scalar>& value,
             const std::vector<double>& history) {
    time.assign(1U, at); state.assign(1U, value); baseline = history;
    jumps.clear();
  }
  void append(const Scalar& at, const VectorT<Scalar>& value) {
    if (!time.empty() && std::abs(scalar_value(time.back() - at)) <= 1e-12) {
      double maximum = 0.0;
      for (Eigen::Index index = 0; index < value.size(); ++index) {
        maximum = std::max(
          maximum, std::abs(scalar_value(state.back()[index] - value[index])));
      }
      if (maximum > 1e-14) {
        if (!jumps.empty() &&
            std::abs(scalar_value(jumps.back().time - at)) <= 1e-12) {
          jumps.back().after = value;
        } else {
          VectorT<Scalar> before = state.back();
          if (time.size() == 1U && jumps.empty() &&
              baseline.size() == static_cast<std::size_t>(before.size())) {
            for (Eigen::Index index = 0; index < before.size(); ++index) {
              before[index] = Scalar(
                baseline[static_cast<std::size_t>(index)]);
            }
          }
          jumps.push_back({at, before, value});
        }
      }
      state.back() = value; return;
    }
    time.push_back(at); state.push_back(value);
  }
  const Jump* jump_at(const Scalar& target) const {
    for (auto jump = jumps.rbegin(); jump != jumps.rend(); ++jump) {
      if (std::abs(scalar_value(jump->time - target)) <= 1e-12) return &*jump;
    }
    return nullptr;
  }
  const Jump* jump_at_index(std::size_t index) const {
    return jump_at(time[index]);
  }
  Scalar at(const Scalar& target, int component, bool left_limit = false) const {
    if (component < 0 || component >= static_cast<int>(baseline.size())) {
      throw std::out_of_range("DDE lag state is outside the state vector.");
    }
    const double target_value = scalar_value(target);
    if (time.empty() || target_value < scalar_value(time.front()) - 1e-12) {
      return Scalar(baseline[static_cast<std::size_t>(component)]);
    }
    if (const Jump* jump = jump_at(target)) {
      return (left_limit ? jump->before : jump->after)[component];
    }
    if (time.size() == 1U) return state.back()[component];
    if (target_value >= scalar_value(time.back()) - 1e-12) {
      const std::size_t right = time.size() - 1U;
      const std::size_t left = right - 1U;
      const Scalar fraction = (target - time[left]) / (time[right] - time[left]);
      return state[left][component] +
        fraction * (state[right][component] - state[left][component]);
    }
    std::size_t right = 1U;
    while (right < time.size() && scalar_value(time[right]) < target_value) ++right;
    const std::size_t left = right - 1U;
    const Scalar fraction = (target - time[left]) / (time[right] - time[left]);
    const Scalar left_state = state[left][component];
    const Jump* right_jump = jump_at_index(right);
    const Scalar right_state = right_jump == nullptr ?
      state[right][component] : right_jump->before[component];
    return left_state + fraction * (right_state - left_state);
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
    const VectorT<Scalar> update =
      solve_linear(system, rhs_matrix, "Implicit ADVAN Newton").col(0);
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
    if (++attempts > control.max_steps) throw std::runtime_error("Implicit ADVAN solver exceeded ODE_CONTROL$max_steps.");
    h = std::min(h, to - time);
    const double minimum = 32.0 * std::numeric_limits<double>::epsilon() *
      std::max({1.0, std::abs(time), std::abs(to)});
    if (h < minimum) throw std::runtime_error("Implicit ADVAN ODE step size underflow.");
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
        Scalar step_end = remaining <= engine.dde_step + 1e-12 ?
          segment_end : time + Scalar(engine.dde_step);
        bool ends_at_delayed_jump = false;
        for (const auto& jump : dde_history->jumps) {
          for (const std::string& delay_name : engine.dde_lag_delays) {
            auto delay = parameters.find(delay_name);
            if (delay == parameters.end()) continue;
            const Scalar boundary = jump.time + delay->second;
            if (path_gt(boundary, time + Scalar(1e-12)) &&
                path_lt(boundary, step_end - Scalar(1e-12))) {
              step_end = boundary;
              ends_at_delayed_jump = true;
            } else if (
              path_gt(boundary, time + Scalar(1e-12)) &&
              std::abs(scalar_value(boundary - step_end)) <= 1e-12
            ) {
              step_end = boundary;
              ends_at_delayed_jump = true;
            }
          }
        }
        const Scalar h = step_end - time;
        auto rhs = [&](const Scalar& stage_time,
                       const VectorT<Scalar>& value,
                       bool left_limit = false) {
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
              stage_time - delay->second, engine.dde_lag_states[lag],
              left_limit);
          }
          VectorT<Scalar> derivative = evaluate_derivatives_t(
            engine, data, row, subject, stage_time, value, parameters,
            theta, eta, eta_columns, sigma, mixture_number, dynamic_data, &lag_values);
          derivative += input;
          return derivative;
        };
        const VectorT<Scalar> k1 = rhs(time, state);
        const VectorT<Scalar> k2 = rhs(
          time + Scalar(0.5) * h, state + Scalar(0.5) * h * k1);
        const VectorT<Scalar> k3 = rhs(
          time + Scalar(0.5) * h, state + Scalar(0.5) * h * k2);
        const VectorT<Scalar> k4 = rhs(
          time + h, state + h * k3, ends_at_delayed_jump);
        state += (h / Scalar(6.0)) *
          (k1 + Scalar(2.0) * k2 + Scalar(2.0) * k3 + k4);
        if (ends_at_delayed_jump &&
            std::abs(scalar_value(step_end - segment_end)) <= 1e-12) {
          const VectorT<Scalar> right_derivative = rhs(
            step_end, state, false);
          state -= (
            step_end - Scalar(scalar_value(step_end))
          ) * right_derivative;
        }
        time = step_end;
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
    state = implicit_ode_advan(engine.advan) ?
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
    if (have_previous && !engine.direct_prediction) {
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
    if (engine.direct_prediction) {
      topology.k = MatrixT<Scalar>::Zero(n_state, n_state);
      topology.default_scales.assign(
        static_cast<std::size_t>(n_state), Scalar(1.0));
    } else if (engine.is_ode()) {
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
    const bool dosing = !engine.direct_prediction && amount[row] > 0.0 &&
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
    const auto direct_prediction = parameters.find("F");
    Scalar raw_prediction = Scalar(0.0);
    if (engine.direct_prediction) {
      if (direct_prediction == parameters.end()) {
        throw std::domain_error("Direct $PRED evaluation did not produce F.");
      }
      raw_prediction = direct_prediction->second;
    } else {
      const int observation_cmt =
        cmt[row] > 0 && evid[row] == 0 ? cmt[row] : engine.obs_cmp;
      const int observation_index = compartment_index(
        observation_cmt, topology.default_observation, n_state);
      const Scalar scale = observation_scale_t(
        parameters, data, row, observation_cmt, topology, dynamic_data);
      raw_prediction =
        direct_prediction == parameters.end() ?
        state[observation_index] / scale : direct_prediction->second;
    }
    prediction[static_cast<std::size_t>(row)] = evaluate_post_prediction_t(
      engine, data, row, subject, Scalar(time[row]), state, theta, eta,
      eta_columns, sigma, mixture_number, raw_prediction, parameters,
      dynamic_data);
    previous_k = topology.k;
    previous_parameters = std::move(parameters);
    previous_row = row;
    previous_mixture = mixture_number;
    previous_time = time[row];
    have_previous = true;
  }
  return prediction;
}
