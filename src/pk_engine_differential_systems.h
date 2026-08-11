struct DdeHistory {
  struct DenseSegment {
    double from;
    double to;
    Vector before;
    std::array<Vector, 3> derivative;
  };

  struct Jump {
    double time;
    Vector before;
    Vector after;
  };

  std::vector<double> time;
  std::vector<Vector> state;
  std::vector<double> baseline;
  std::vector<Jump> jumps;
  std::vector<DenseSegment> dense;

  void reset(double at, const Vector& value, const std::vector<double>& history) {
    time.assign(1U, at); state.assign(1U, value); baseline = history;
    jumps.clear(); dense.clear();
  }
  void append(double at, const Vector& value) {
    if (!time.empty() && std::abs(time.back() - at) <= 1e-12) {
      if ((state.back() - value).cwiseAbs().maxCoeff() > 1e-14) {
        if (!jumps.empty() && std::abs(jumps.back().time - at) <= 1e-12) {
          jumps.back().after = value;
        } else {
          Vector before = state.back();
          if (time.size() == 1U && jumps.empty() &&
              baseline.size() >= static_cast<std::size_t>(before.size())) {
            for (Eigen::Index index = 0; index < before.size(); ++index) {
              before[index] = baseline[static_cast<std::size_t>(index)];
            }
          }
          jumps.push_back({at, before, value});
        }
      }
      state.back() = value; return;
    }
    time.push_back(at); state.push_back(value);
  }
  void append_radau(double from, double to, const Vector& before,
                    const std::array<Vector, 3>& derivative,
                    const Vector& after) {
    dense.push_back({from, to, before, derivative});
    append(to, after);
  }
  const Jump* jump_at(double target) const {
    for (auto jump = jumps.rbegin(); jump != jumps.rend(); ++jump) {
      if (std::abs(jump->time - target) <= 1e-12) return &*jump;
    }
    return nullptr;
  }
  const Jump* jump_at_index(std::size_t index) const {
    return jump_at(time[index]);
  }
  double at(double target, int component, bool left_limit = false) const {
    if (component < 0 || component >= static_cast<int>(baseline.size())) {
      throw std::out_of_range("DDE lag state is outside the state vector.");
    }
    if (time.empty() || target < time.front() - 1e-12) return baseline[static_cast<std::size_t>(component)];
    if (const Jump* jump = jump_at(target)) {
      return (left_limit ? jump->before : jump->after)[component];
    }
    // At a shared continuous mesh point prefer the segment ending there.  Its
    // stiffly accurate endpoint derivative is the DDE left limit; explicit
    // jumps have already been handled above.
    for (auto segment = dense.begin(); segment != dense.end(); ++segment) {
      if (target < segment->from - 1e-12 || target > segment->to + 1e-12) continue;
      const double width = segment->to - segment->from;
      if (!(width > 0.0)) continue;
      const double theta = std::clamp((target - segment->from) / width, 0.0, 1.0);
      const double root = std::sqrt(6.0);
      const std::array<double, 3> nodes = {
        (4.0 - root) / 10.0, (4.0 + root) / 10.0, 1.0
      };
      double value = segment->before[component];
      for (std::size_t stage = 0; stage < 3U; ++stage) {
        const std::size_t first = (stage + 1U) % 3U;
        const std::size_t second = (stage + 2U) % 3U;
        const double denominator =
          (nodes[stage] - nodes[first]) * (nodes[stage] - nodes[second]);
        const double q2 = 1.0 / denominator;
        const double q1 = -(nodes[first] + nodes[second]) / denominator;
        const double q0 = nodes[first] * nodes[second] / denominator;
        const double integral = q2 * theta * theta * theta / 3.0 +
          q1 * theta * theta / 2.0 + q0 * theta;
        value += width * integral * segment->derivative[stage][component];
      }
      return value;
    }
    if (target >= time.back() - 1e-12) return state.back()[component];
    auto upper = std::upper_bound(time.begin(), time.end(), target);
    const std::size_t right = static_cast<std::size_t>(std::distance(time.begin(), upper));
    const std::size_t left = right - 1U;
    const double fraction = (target - time[left]) / (time[right] - time[left]);
    const double left_state = state[left][component];
    const Jump* right_jump = jump_at_index(right);
    const double right_state = right_jump == nullptr ?
      state[right][component] : right_jump->before[component];
    return left_state + fraction * (right_state - left_state);
  }
};

Vector evaluate_derivatives(const ModelEngine& engine,
                            const EventDataView& data,
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
    const CompiledInputBinding& binding = engine.des_inputs[i];
    const std::string& name = binding.name;
    int index = binding.state;
    if (index >= 0) {
      if (index >= state.size()) throw std::out_of_range("A() index exceeds the ODE state dimension.");
      inputs[i] = state[index];
      continue;
    }
    if (name == "T") {
      inputs[i] = t;
      continue;
    }
    if (binding.lag >= 0) {
      if (lag_values == nullptr) throw std::logic_error("DDE lag history is unavailable.");
      inputs[i] = (*lag_values)[binding.lag];
      continue;
    }
    if (binding.algebraic >= 0) {
      inputs[i] = algebraic[binding.algebraic];
      continue;
    }
    auto parameter = parameters.find(name);
    if (parameter != parameters.end()) {
      inputs[i] = parameter->second;
      continue;
    }
    index = binding.theta;
    if (index >= 0) {
      if (index >= theta.size()) throw std::out_of_range("THETA index exceeds supplied values in DES.");
      inputs[i] = theta[index];
      continue;
    }
    index = binding.eta;
    if (index >= 0) {
      inputs[i] = eta(subject, eta_column(engine, data, row, index, eta.ncol()));
      continue;
    }
    index = binding.sigma;
    if (index >= 0) {
      if (index >= sigma.size()) throw std::out_of_range("SIGMA index exceeds supplied values in DES.");
      inputs[i] = sigma[index];
      continue;
    }
    if (binding.zero_error) {
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
      throw std::runtime_error("Implicit ADVAN solver exceeded ODE_CONTROL$max_steps.");
    }
    h = std::min(h, to - t);
    const double minimum = 32.0 * std::numeric_limits<double>::epsilon() *
      std::max({1.0, std::abs(t), std::abs(to)});
    if (h < minimum) throw std::runtime_error("Implicit ADVAN ODE step size underflow.");

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

double row_optional(const EventDataView& data, const std::string& name,
                    int row, double fallback) {
  if (!data.containsElementNamed(name.c_str())) return fallback;
  double value = data_value(data, name, row);
  return std::isfinite(value) ? value : fallback;
}

double bioavailability(const Parameters& p, const EventDataView& data,
                       int row, int cmt) {
  const std::string name = "F" + std::to_string(cmt);
  double value = get_positive(p, {name.c_str()});
  if (!finite_positive(value)) value = row_optional(data, name, row, 1.0);
  return finite_positive(value) ? value : 1.0;
}

double event_infusion_rate(const Parameters& p, const EventDataView& data,
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

double observation_scale(const Parameters& p, const EventDataView& data,
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

struct RadauDdeStep {
  Vector after;
  std::array<Vector, 3> derivative;
  double error = std::numeric_limits<double>::infinity();
  bool converged = false;
};

template <class Rhs>
RadauDdeStep radau_iia5_dde_step(const Rhs& rhs, const Vector& before,
                                 double time, double h,
                                 const OdeControl& control,
                                 bool ends_at_delayed_jump,
                                 bool optimized = false) {
  // Clean-room implementation of the published three-stage, fifth-order
  // Radau IIA collocation formula.  It deliberately uses the full coupled
  // Newton matrix rather than any implementation-specific RADAR5 transform.
  const double root = std::sqrt(6.0);
  const std::array<double, 3> c = {
    (4.0 - root) / 10.0, (4.0 + root) / 10.0, 1.0
  };
  const double a[3][3] = {
    {(88.0 - 7.0 * root) / 360.0,
     (296.0 - 169.0 * root) / 1800.0,
     (-2.0 + 3.0 * root) / 225.0},
    {(296.0 + 169.0 * root) / 1800.0,
     (88.0 + 7.0 * root) / 360.0,
     (-2.0 - 3.0 * root) / 225.0},
    {(16.0 - root) / 36.0,
     (16.0 + root) / 36.0,
     1.0 / 9.0}
  };
  const Eigen::Index n = before.size();
  const Vector f0 = rhs(time, before, false);
  std::array<Vector, 3> stage = {
    before + c[0] * h * f0,
    before + c[1] * h * f0,
    before + h * f0
  };
  std::array<Vector, 3> derivative = {f0, f0, f0};
  std::array<Matrix, 3> jacobian;
  bool jacobian_ready = false;
  Eigen::FullPivLU<Matrix> optimized_lu;
  if (optimized) {
    Matrix step_jacobian(n, n);
    for (Eigen::Index column = 0; column < n; ++column) {
      Vector perturbed = before;
      const double delta = std::sqrt(std::numeric_limits<double>::epsilon()) *
        std::max(1.0, std::abs(before[column]));
      perturbed[column] += delta;
      step_jacobian.col(column) =
        (rhs(time, perturbed, false) - f0) / delta;
    }
    Matrix system = Matrix::Zero(3 * n, 3 * n);
    for (std::size_t i = 0; i < 3U; ++i) {
      for (std::size_t j = 0; j < 3U; ++j) {
        system.block(static_cast<Eigen::Index>(i) * n,
                     static_cast<Eigen::Index>(j) * n, n, n) =
          -h * a[i][j] * step_jacobian;
        if (i == j) {
          system.block(static_cast<Eigen::Index>(i) * n,
                       static_cast<Eigen::Index>(j) * n, n, n).diagonal().array() += 1.0;
        }
      }
    }
    optimized_lu.compute(system);
    if (!optimized_lu.isInvertible()) return RadauDdeStep();
  }
  RadauDdeStep result;
  for (int iteration = 0; iteration < 12; ++iteration) {
    Vector residual(3 * n);
    for (std::size_t j = 0; j < 3U; ++j) {
      const bool left = ends_at_delayed_jump && j == 2U;
      derivative[j] = rhs(time + c[j] * h, stage[j], left);
      if (!optimized && !jacobian_ready) {
        jacobian[j].resize(n, n);
        for (Eigen::Index column = 0; column < n; ++column) {
          Vector perturbed = stage[j];
          const double delta = std::sqrt(std::numeric_limits<double>::epsilon()) *
            std::max(1.0, std::abs(stage[j][column]));
          perturbed[column] += delta;
          jacobian[j].col(column) =
            (rhs(time + c[j] * h, perturbed, left) - derivative[j]) / delta;
        }
      }
    }
    jacobian_ready = true;
    for (std::size_t i = 0; i < 3U; ++i) {
      Vector local = stage[i] - before;
      for (std::size_t j = 0; j < 3U; ++j) local -= h * a[i][j] * derivative[j];
      residual.segment(static_cast<Eigen::Index>(i) * n, n) = local;
    }
    double residual_error = 0.0;
    for (std::size_t i = 0; i < 3U; ++i) {
      residual_error = std::max(
        residual_error,
        scaled_error(
          residual.segment(static_cast<Eigen::Index>(i) * n, n),
          before, stage[i], control));
    }
    if (residual_error < 0.03) {
      result.converged = true;
      break;
    }
    Vector update;
    if (optimized) {
      update = optimized_lu.solve(-residual);
    } else {
      Matrix system = Matrix::Zero(3 * n, 3 * n);
      for (std::size_t i = 0; i < 3U; ++i) {
        for (std::size_t j = 0; j < 3U; ++j) {
          system.block(static_cast<Eigen::Index>(i) * n,
                       static_cast<Eigen::Index>(j) * n, n, n) =
            -h * a[i][j] * jacobian[j];
          if (i == j) {
            system.block(static_cast<Eigen::Index>(i) * n,
                         static_cast<Eigen::Index>(j) * n, n, n).diagonal().array() += 1.0;
          }
        }
      }
      Eigen::FullPivLU<Matrix> lu(system);
      if (!lu.isInvertible()) return result;
      update = lu.solve(-residual);
    }
    if (!update.allFinite()) return result;
    double update_error = 0.0;
    for (std::size_t i = 0; i < 3U; ++i) {
      const Vector local = update.segment(static_cast<Eigen::Index>(i) * n, n);
      stage[i] += local;
      update_error = std::max(update_error,
        scaled_error(local, before, stage[i], control));
    }
    if (update_error < 0.03) {
      result.converged = true;
      break;
    }
  }
  if (!result.converged) return result;
  for (std::size_t j = 0; j < 3U; ++j) {
    derivative[j] = rhs(
      time + c[j] * h, stage[j], ends_at_delayed_jump && j == 2U);
  }
  result.after = stage[2];  // Radau IIA is stiffly accurate.
  result.derivative = derivative;

  // An independently constructed order-three embedded formula supplies a
  // local error estimate.  gamma is deliberately explicit so the estimator
  // follows from the B(3) moment equations, not from vendor implementation.
  constexpr double gamma = 0.1;
  const std::array<double, 3> embedded = {
    4.0 / 9.0 - root / 36.0 - gamma / 3.0 - root * gamma / 2.0,
    4.0 / 9.0 + root / 36.0 - gamma / 3.0 + root * gamma / 2.0,
    1.0 / 9.0 - gamma / 3.0
  };
  Vector lower = before + h * gamma * f0;
  for (std::size_t stage_index = 0; stage_index < 3U; ++stage_index) {
    lower += h * embedded[stage_index] * derivative[stage_index];
  }
  result.error = scaled_error(result.after - lower, before, result.after, control);
  if (!std::isfinite(result.error)) result.converged = false;
  return result;
}

struct PropagationCacheEntry {
  Matrix k;
  Vector input;
  double interval = 0.0;
  AffineMap map;
};

class PropagationCache {
 public:
  const AffineMap& get(const Matrix& k, const Vector& input, double interval) {
    for (const PropagationCacheEntry& entry : entries_) {
      if (entry.interval == interval && entry.k.rows() == k.rows() &&
          entry.k.cols() == k.cols() && entry.input.size() == input.size() &&
          (entry.k.array() == k.array()).all() &&
          (entry.input.array() == input.array()).all()) {
        ++hits_;
        return entry.map;
      }
    }
    if (entries_.size() >= 64U) entries_.clear();
    entries_.push_back({k, input, interval, affine_map(k, input, interval)});
    ++misses_;
    return entries_.back().map;
  }

  void clear() { entries_.clear(); }
  std::size_t hits() const { return hits_; }
  std::size_t misses() const { return misses_; }

 private:
  std::vector<PropagationCacheEntry> entries_;
  std::size_t hits_ = 0U;
  std::size_t misses_ = 0U;
};

Vector propagate_to(const Matrix& k, Vector state, double from, double to,
                    std::vector<ActiveInfusion>& active,
                    PropagationCache* cache = nullptr) {
  double cursor = from;
  remove_finished(active, cursor);
  while (cursor < to - 1e-12) {
    double segment_end = to;
    for (const ActiveInfusion& infusion : active) {
      if (infusion.end > cursor + 1e-12) segment_end = std::min(segment_end, infusion.end);
    }
    const Vector input = infusion_input(k.rows(), active);
    if (cache == nullptr) {
      state = propagate(k, input, segment_end - cursor, state);
    } else {
      const AffineMap& map = cache->get(k, input, segment_end - cursor);
      state = map.transition * state + map.offset;
    }
    cursor = segment_end;
    remove_finished(active, cursor);
  }
  return state;
}

Vector propagate_ode_to(const ModelEngine& engine,
                        const EventDataView& data,
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
      double suggested_step = engine.dde_step;
      int steps = 0;
      while (time < segment_end - 1e-12) {
        if (++steps > engine.dde_max_steps) {
          throw std::runtime_error("DDE method-of-steps exceeded max_steps.");
        }
        double step_end = time + std::min(suggested_step, segment_end - time);
        bool ends_at_delayed_jump = false;
        for (const auto& jump : dde_history->jumps) {
          for (const std::string& delay_name : engine.dde_lag_delays) {
            auto delay = parameters.find(delay_name);
            if (delay == parameters.end()) continue;
            const double boundary = jump.time + delay->second;
            if (boundary > time + 1e-12 &&
                boundary < step_end - 1e-12) {
              step_end = boundary;
              ends_at_delayed_jump = true;
            } else if (boundary > time + 1e-12 &&
                       std::abs(boundary - step_end) <= 1e-12) {
              step_end = boundary;
              ends_at_delayed_jump = true;
            }
          }
        }
        const double h = step_end - time;
        std::vector<double> lag_cache_time;
        std::vector<bool> lag_cache_left;
        std::vector<Vector> lag_cache_value;
        auto rhs = [&](double stage_time, const Vector& value,
                       bool left_limit = false) {
          // The minimum-delay contract keeps lagged values entirely in
          // accepted history.  Reusing the value at an identical collocation
          // time therefore cannot depend on the current Newton iterate and is
          // safe for both numerical policies.
          for (std::size_t cached = 0; cached < lag_cache_time.size(); ++cached) {
            if (lag_cache_left[cached] == left_limit &&
                std::abs(lag_cache_time[cached] - stage_time) <= 1e-14) {
              Vector derivative = evaluate_derivatives(
                engine, data, row, subject, stage_time, value, parameters,
                theta, eta, sigma, &lag_cache_value[cached]);
              derivative += input;
              return derivative;
            }
          }
          Vector lag_values(static_cast<Eigen::Index>(engine.dde_lag_inputs.size()));
          for (std::size_t lag = 0; lag < engine.dde_lag_inputs.size(); ++lag) {
            auto delay = parameters.find(engine.dde_lag_delays[lag]);
            if (delay == parameters.end() || !std::isfinite(delay->second) ||
                delay->second < engine.dde_minimum_delay || delay->second < h) {
              throw std::domain_error("DDE delay '" + engine.dde_lag_delays[lag] +
                                      "' must be finite and at least the integration step.");
            }
            const double target = stage_time - delay->second;
            const int lag_state = engine.dde_lag_states[lag];
            if (lag_state < engine.n_state ||
                target < dde_history->time.front() - 1e-12) {
              lag_values[static_cast<Eigen::Index>(lag)] = dde_history->at(
                target, lag_state, left_limit);
            } else {
              Vector delayed_state(engine.n_state);
              for (int component = 0; component < engine.n_state; ++component) {
                delayed_state[component] = dde_history->at(
                  target, component, left_limit);
              }
              const Vector delayed_algebraic = solve_algebraic(
                engine, data, row, subject, target, delayed_state, parameters,
                theta, eta, sigma);
              const int algebraic = lag_state - engine.n_state;
              if (algebraic < 0 || algebraic >= delayed_algebraic.size()) {
                throw std::out_of_range(
                  "DDE lag algebraic state is outside the DAE variable vector.");
              }
              lag_values[static_cast<Eigen::Index>(lag)] =
                delayed_algebraic[algebraic];
            }
          }
          lag_cache_time.push_back(stage_time);
          lag_cache_left.push_back(left_limit);
          lag_cache_value.push_back(lag_values);
          Vector derivative = evaluate_derivatives(
            engine, data, row, subject, stage_time, value, parameters,
            theta, eta, sigma, &lag_values);
          derivative += input;
          return derivative;
        };
        if (engine.advan == 16 || engine.advan == 17) {
          const RadauDdeStep trial = radau_iia5_dde_step(
            rhs, state, time, h, engine.ode_control, ends_at_delayed_jump,
            engine.liber_optimized);
          const double minimum = 32.0 * std::numeric_limits<double>::epsilon() *
            std::max({1.0, std::abs(time), std::abs(segment_end)});
          if (trial.converged && trial.error <= 1.0) {
            const Vector before = state;
            state = trial.after;
            dde_history->append_radau(time, step_end, before,
                                      trial.derivative, state);
            time = step_end;
          }
          const double factor = trial.converged && trial.error == 0.0 ? 3.0 :
            (trial.converged && std::isfinite(trial.error) ?
              std::clamp(0.9 * std::pow(trial.error, -0.25), 0.2, 3.0) : 0.25);
          suggested_step = std::min(engine.dde_step, h * factor);
          if (suggested_step < minimum) {
            throw std::runtime_error(
              "ADVAN" + std::to_string(engine.advan) +
              " Radau IIA DDE step size underflow.");
          }
          continue;
        }
        const Vector k1 = rhs(time, state);
        const Vector k2 = rhs(time + 0.5 * h, state + 0.5 * h * k1);
        const Vector k3 = rhs(time + 0.5 * h, state + 0.5 * h * k2);
        const Vector k4 = rhs(
          time + h, state + h * k3, ends_at_delayed_jump);
        state += (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        time = step_end;
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
    state = implicit_ode_advan(engine.advan) ?
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
    const ModelEngine& engine, const EventDataView& data,
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
    const ModelEngine& engine, const EventDataView& data,
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
                    const EventDataView& data,
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

  auto time = data.values("TIME");
  auto amount = data.values("AMT");
  auto rate = data.values("RATE");
  auto interval = data.values("II");
  auto evid = data.values("EVID");
  auto cmt = data.values("CMT");
  auto ss = data.values("SS");
  auto eta_subject_index = data.values(".ID_INDEX");
  auto subject_index = data.containsElementNamed(".STRUCT_ID_INDEX") ?
    data.values(".STRUCT_ID_INDEX") : eta_subject_index;
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
  PropagationCache propagation_cache;

  for (int row = 0; row < n_rows; ++row) {
    const int structural_subject = subject_index[row] - 1;
    const int subject = eta_subject_index[row] - 1;
    if (structural_subject != previous_subject) {
      state.setZero();
      active.clear();
      have_previous = false;
      previous_time = time[row];
      previous_subject = structural_subject;
      propagation_cache.clear();
      if (engine.dde_enabled) dde_history.reset(time[row], state, engine.dde_history);
    }
    if (have_previous && !engine.direct_prediction) {
      if (time[row] < previous_time - 1e-12) Rcpp::stop("Subject event times are decreasing.");
      if (engine.is_ode()) {
        state = propagate_ode_to(
          engine, data, previous_row, subject, theta, eta, sigma,
          previous_parameters, state, previous_time, time[row], active
          , engine.dde_enabled ? &dde_history : nullptr
        );
      } else {
        // Exact (K, input, interval) matches reuse the same affine map.  The
        // state multiplication and event order remain unchanged.
        state = propagate_to(
          previous_k, state, previous_time, time[row], active,
          &propagation_cache);
      }
    }

    Parameters parameters = evaluate_parameters(engine, data, row, subject, theta, eta, sigma);
    Topology topology;
    if (engine.direct_prediction) {
      topology.k = Matrix::Zero(n_state, n_state);
      topology.state_names = {"DIRECT_PRED"};
      topology.default_scales.assign(static_cast<std::size_t>(n_state), 1.0);
    } else if (engine.is_ode()) {
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
    const bool dosing = !engine.direct_prediction && amount[row] > 0.0 &&
      (evid[row] == 1 || evid[row] == 4 || evid[row] == 0);
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

    const auto direct_prediction = parameters.find("F");
    double raw_prediction = NA_REAL;
    if (engine.direct_prediction) {
      if (direct_prediction == parameters.end()) {
        Rcpp::stop("Direct $PRED evaluation did not produce F.");
      }
      raw_prediction = direct_prediction->second;
    } else {
      const int observation_cmt =
        cmt[row] > 0 && evid[row] == 0 ? cmt[row] : engine.obs_cmp;
      const int observation_index = compartment_index(
        observation_cmt, topology.default_observation, n_state
      );
      const double scale = observation_scale(
        parameters, data, row, observation_cmt, topology);
      raw_prediction = direct_prediction == parameters.end() ?
        state[observation_index] / scale : direct_prediction->second;
    }
    prediction[row] = evaluate_post_prediction(
      engine, data, row, subject, time[row], state, theta, eta, sigma,
      raw_prediction, parameters);
    for (std::size_t output = 0;
         output < engine.selected_output_names.size(); ++output) {
      const auto found = parameters.find(engine.selected_output_names[output]);
      generated(row, static_cast<int>(output)) =
        found == parameters.end() ? NA_REAL : found->second;
    }
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
      (engine.is_ode() ? ode_kernel_name(engine.advan) :
       (engine.matrix_graph.enabled ? "general-matrix-exponential" : "advan")) :
      engine.solver
  );
}

// Scalar-generic analytical path used to record the *complete* event and
// compartment calculation with CppAD.  The ordinary simulation path above is
// deliberately retained as an independently testable double-precision
// implementation.  Keeping the two entry points separate also prevents an R
// callback or an adaptive solver decision from being hidden inside an AD
// recording.
