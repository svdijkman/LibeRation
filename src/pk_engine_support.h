// Internal implementation fragment included inside namespace liberation by
// pk_engine.cpp. It intentionally has no include guard or namespace declaration.
// Contains CppAD recording guards and small shared engine structures.


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

