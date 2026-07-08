#pragma once

namespace nm_lik {

enum ErrorType : int {
  ERR_PROPADD = 0,
  ERR_ADD = 1,
  ERR_PROP = 2,
  ERR_LOG = 3,
  ERR_POWER = 4
};

enum OmegaType : int {
  OMEGA_DIAG = 0,
  OMEGA_BLOCK2 = 1
};

enum SigmaCorr : int {
  SIGMA_INDEP = 0,
  SIGMA_AR1 = 1
};

// Below-limit-of-quantification (censoring) likelihood method.
enum BlqMethod : int {
  BLQ_NONE = 0,
  BLQ_M3 = 1,
  BLQ_M4 = 2
};

struct LikConfig {
  int error_type = ERR_PROPADD;
  int omega_type = OMEGA_DIAG;
  int sigma_corr = SIGMA_INDEP;
  int iov = 0;
  double ar1_rho = 0.0;
  int blq_method = BLQ_NONE;
};

inline LikConfig& config() {
  static LikConfig cfg;
  return cfg;
}

inline void set_config(int error_type, int omega_type,
                       int sigma_corr = 0, int iov = 0,
                       double ar1_rho = 0.0, int blq_method = 0) {
  config().error_type = error_type;
  config().omega_type = omega_type;
  config().sigma_corr = sigma_corr;
  config().iov = iov;
  config().ar1_rho = ar1_rho;
  config().blq_method = blq_method;
}

}  // namespace nm_lik
