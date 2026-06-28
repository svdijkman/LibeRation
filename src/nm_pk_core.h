#pragma once
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <string>

namespace nm_pk {

constexpr int kMaxPkState = 7;

inline thread_local Rcpp::NumericMatrix* pk_route_amounts_ptr = nullptr;

inline void pk_write_amounts(int row, const std::vector<double>& st) {
  if (pk_route_amounts_ptr == nullptr || row < 0 ||
      row >= pk_route_amounts_ptr->nrow()) {
    return;
  }
  const int nc = std::min(static_cast<int>(st.size()), pk_route_amounts_ptr->ncol());
  for (int j = 0; j < nc; ++j) {
    (*pk_route_amounts_ptr)(row, j) = st[j];
  }
}

inline double safe_div(double a, double b, double eps = 1e-12) {
  return (std::fabs(b) < eps) ? 0.0 : a / b;
}

inline double lhopital_exp_diff(double ka, double ke, double dt) {
  if (std::fabs(ka - ke) < 1e-12) {
    return ka * dt * std::exp(-ka * dt);
  }
  return ka / (ka - ke) * (std::exp(-ke * dt) - std::exp(-ka * dt));
}

struct Hybrid2 {
  double lambda1, lambda2, E1, E2;
};

inline Hybrid2 hybrid_lambdas_2(double k10, double k12, double k21) {
  Hybrid2 h;
  h.E1 = k10 + k12;
  h.E2 = k21;
  const double disc = (h.E1 + h.E2) * (h.E1 + h.E2) - 4.0 * (h.E1 * h.E2 - k12 * k21);
  const double root = std::sqrt(std::max(disc, 0.0));
  h.lambda1 = 0.5 * ((h.E1 + h.E2) + root);
  h.lambda2 = 0.5 * ((h.E1 + h.E2) - root);
  return h;
}

struct Hybrid3 {
  double lambda1, lambda2, lambda3;
  double E1, E2, E3;
};

inline Hybrid3 hybrid_lambdas_3_iv(double k10, double k12, double k21, double k13, double k31) {
  Hybrid3 h;
  h.E1 = k10 + k12 + k13;
  h.E2 = k21;
  h.E3 = k31;
  const double a = h.E1 + h.E2 + h.E3;
  const double b = h.E1 * h.E2 + h.E3 * (h.E1 + h.E2) - k12 * k21 - k13 * k31;
  const double c = h.E1 * h.E2 * h.E3 - h.E3 * k12 * k21 - h.E2 * k13 * k31;
  const double m = (3.0 * b - a * a) / 3.0;
  const double n = (2.0 * a * a * a - 9.0 * a * b + 27.0 * c) / 27.0;
  const double Q = n * n / 4.0 + m * m * m / 27.0;
  const double alpha = std::sqrt(std::max(-Q, 0.0));
  const double beta = -n / 2.0;
  const double gamma = std::sqrt(beta * beta + alpha * alpha);
  const double theta = std::atan2(alpha, beta);
  const double g13 = std::cbrt(gamma);
  h.lambda1 = a / 3.0 + g13 * (std::cos(theta / 3.0) + std::sqrt(3.0) * std::sin(theta / 3.0));
  h.lambda2 = a / 3.0 + g13 * (std::cos(theta / 3.0) - std::sqrt(3.0) * std::sin(theta / 3.0));
  h.lambda3 = a / 3.0 - 2.0 * g13 * std::cos(theta / 3.0);
  return h;
}

inline Hybrid3 hybrid_lambdas_3_oral(double k20, double k23, double k32, double k24, double k42) {
  Hybrid3 h;
  h.E1 = k20 + k23 + k24;
  h.E2 = k32;
  h.E3 = k42;
  const double a = h.E1 + h.E2 + h.E3;
  const double b = h.E1 * h.E2 + h.E3 * (h.E1 + h.E2) - k23 * k32 - k24 * k42;
  const double c = h.E1 * h.E2 * h.E3 - h.E3 * k23 * k32 - h.E2 * k24 * k42;
  const double m = (3.0 * b - a * a) / 3.0;
  const double n = (2.0 * a * a * a - 9.0 * a * b + 27.0 * c) / 27.0;
  const double Q = n * n / 4.0 + m * m * m / 27.0;
  const double alpha = std::sqrt(std::max(-Q, 0.0));
  const double beta = -n / 2.0;
  const double gamma = std::sqrt(beta * beta + alpha * alpha);
  const double theta = std::atan2(alpha, beta);
  const double g13 = std::cbrt(gamma);
  h.lambda1 = a / 3.0 + g13 * (std::cos(theta / 3.0) + std::sqrt(3.0) * std::sin(theta / 3.0));
  h.lambda2 = a / 3.0 + g13 * (std::cos(theta / 3.0) - std::sqrt(3.0) * std::sin(theta / 3.0));
  h.lambda3 = a / 3.0 - 2.0 * g13 * std::cos(theta / 3.0);
  return h;
}

enum RouteType {
  ROUTE_IV_BOLUS = 1,
  ROUTE_IV_INFUSION = 2,
  ROUTE_ORAL = 3,
  ROUTE_TRANSIT = 4,
  ROUTE_ODE = 5
};

inline RouteType route_type(int advan, int trans, int cmt, double rate, int n_transit = 0) {
  if (n_transit > 0) return ROUTE_TRANSIT;
  if (rate > 0.0) return ROUTE_IV_INFUSION;
  if (rate <= -2.0) return ROUTE_IV_INFUSION;
  if (rate < 0.0) {
    if (cmt == 1 && (advan == 2 || advan == 4 || advan == 12 || advan == 6 || advan == 13)) {
      return ROUTE_ORAL;
    }
    return ROUTE_IV_BOLUS;
  }
  if (cmt == 1 && rate <= 0.0) {
    if (advan == 2 || advan == 4 || advan == 12) return ROUTE_ORAL;
    if (advan == 6 || advan == 13) return ROUTE_ORAL;
  }
  return ROUTE_IV_BOLUS;
}

static constexpr int kMaxScale = 10;

struct PkParams {
  double ka = 0.0;
  double k10 = 0.0, k12 = 0.0, k21 = 0.0, k13 = 0.0, k31 = 0.0;
  double k20 = 0.0, k23 = 0.0, k32 = 0.0, k24 = 0.0, k42 = 0.0, k30 = 0.0;
  double ktr = 0.0, kmf = 0.0, kme = 0.0;
  double cl = 0.0, vc = 0.0, vp = 0.0, vp2 = 0.0, q2 = 0.0, q3 = 0.0, q4 = 0.0;
  double v1 = 0.0, v2 = 0.0, v3 = 0.0, v4 = 0.0;
  double vss = 0.0;
  double aob = 0.0;
  double pk_alpha = 0.0;
  double pk_beta = 0.0;
  double pk_gamma = 0.0;
  double vmax = 0.0, km = 0.0, vm = 0.0;
  double scale[kMaxScale] = {0.0};
  double f[kMaxScale] = {0.0};
  double alag1 = 0.0;
  int n_transit = 0;
  bool has_metabolite = false;
};

inline double effective_ke(const PkParams& p) {
  return p.k20 > 0.0 ? p.k20 : p.k10;
}

inline int pk_ncomp(int advan) {
  if (advan == 10) return 1;
  if (advan == 1 || advan == 2) return 1;
  if (advan == 11 || advan == 12) return 3;
  if (advan == 3 || advan == 4) return 2;
  return 2;
}

inline void normalize_volume_aliases(PkParams& p) {
  if (p.vc <= 0.0 && p.v2 > 0.0) p.vc = p.v2;
  if (p.vp <= 0.0 && p.v3 > 0.0) p.vp = p.v3;
  if (p.v1 <= 0.0 && p.vc > 0.0) p.v1 = p.vc;
}

inline double pk_central_volume(const PkParams& p) {
  if (p.vc > 0.0) return p.vc;
  if (p.v1 > 0.0) return p.v1;
  if (p.v2 > 0.0) return p.v2;
  return 0.0;
}

inline double pk_oral_central_volume(const PkParams& p, int advan, int trans) {
  if (advan == 4 || advan == 12) {
    if (trans == 4 || trans == 3 || trans == 5 || trans == 6) {
      if (p.v2 > 0.0) return p.v2;
    }
  }
  return pk_central_volume(p);
}

inline double pk_peripheral_volume(const PkParams& p, int index = 1) {
  if (index == 2) {
    if (p.v4 > 0.0) return p.v4;
    if (p.vp2 > 0.0) return p.vp2;
    return 0.0;
  }
  if (p.v3 > 0.0) return p.v3;
  if (p.vp > 0.0) return p.vp;
  if (p.v2 > 0.0) return p.v2;
  return 0.0;
}

inline void ensure_vss(PkParams& p) {
  if (p.vss > 0.0) return;
  const double v = pk_central_volume(p);
  if (v > 0.0 && p.v2 > 0.0) p.vss = v + p.v2;
}

inline void set_elimination_rate(PkParams& p, double k, bool oral) {
  if (!(k > 0.0)) return;
  if (oral) {
    p.k20 = k;
    p.k10 = k;
  } else {
    p.k10 = k;
  }
}

inline void apply_trans1_oral_aliases(int advan, PkParams& p) {
  if (advan == 12) {
    if (p.k20 <= 0.0 && p.k10 > 0.0) p.k20 = p.k10;
    if (p.k12 > 0.0 && p.k21 > 0.0) {
      const double pk24 = p.k23 > 0.0 ? p.k23 : p.k24;
      const double pk42 = p.k32 > 0.0 ? p.k32 : p.k42;
      p.k23 = p.k12;
      p.k32 = p.k21;
      if (pk24 > 0.0) p.k24 = pk24;
      if (pk42 > 0.0) p.k42 = pk42;
    }
    return;
  }
  if (p.k20 <= 0.0 && p.k10 > 0.0) p.k20 = p.k10;
  if (p.k23 <= 0.0 && p.k12 > 0.0) p.k23 = p.k12;
  if (p.k32 <= 0.0 && p.k21 > 0.0) p.k32 = p.k21;
  if (p.k24 <= 0.0 && p.k13 > 0.0) p.k24 = p.k13;
  if (p.k42 <= 0.0 && p.k31 > 0.0) p.k42 = p.k31;
}

inline void apply_trans2(int advan, bool oral, PkParams& p) {
  const double vol = pk_central_volume(p);
  if (p.cl > 0.0 && vol > 0.0) {
    set_elimination_rate(p, p.cl / vol, oral);
  }
}

inline void apply_trans3(int advan, bool oral, PkParams& p) {
  ensure_vss(p);
  const double v = pk_central_volume(p);
  const double q = p.q2;
  const double v_per = p.vss > v ? (p.vss - v) : pk_peripheral_volume(p, 1);
  if (!(p.cl > 0.0 && v > 0.0)) return;
  const double k = p.cl / v;
  if (advan == 4 && oral) {
    set_elimination_rate(p, k, true);
    if (q > 0.0) {
      p.k23 = q / v;
      if (v_per > 0.0) p.k32 = q / v_per;
    }
    return;
  }
  if (advan == 3 && !oral) {
    set_elimination_rate(p, k, false);
    if (q > 0.0) {
      p.k12 = q / v;
      if (v_per > 0.0) p.k21 = q / v_per;
    }
  }
}

inline void apply_trans4(int advan, bool oral, PkParams& p) {
  const double q1 = p.q2;
  const double q2nd = p.q3 > 0.0 ? p.q3 : p.q4;
  if (advan == 3 && !oral) {
    const double v1 = p.v1 > 0.0 ? p.v1 : pk_central_volume(p);
    const double v2 = p.v2 > 0.0 ? p.v2 : pk_peripheral_volume(p, 1);
    if (p.cl > 0.0 && v1 > 0.0) {
      p.k10 = p.cl / v1;
      if (q1 > 0.0) {
        p.k12 = q1 / v1;
        if (v2 > 0.0) p.k21 = q1 / v2;
      }
    }
    return;
  }
  if (advan == 4 && oral) {
    const double vc = pk_oral_central_volume(p, advan, 4);
    const double vp = pk_peripheral_volume(p, 1);
    if (p.cl > 0.0 && vc > 0.0) {
      set_elimination_rate(p, p.cl / vc, true);
      if (q1 > 0.0) {
        p.k23 = q1 / vc;
        if (vp > 0.0) p.k32 = q1 / vp;
      }
    }
    return;
  }
  if (advan == 11 && !oral) {
    const double v1 = p.v1 > 0.0 ? p.v1 : pk_central_volume(p);
    const double v2 = p.v2 > 0.0 ? p.v2 : pk_peripheral_volume(p, 1);
    const double v3 = p.v3 > 0.0 ? p.v3 : pk_peripheral_volume(p, 2);
    if (p.cl > 0.0 && v1 > 0.0) {
      p.k10 = p.cl / v1;
      if (q1 > 0.0) {
        p.k12 = q1 / v1;
        if (v2 > 0.0) p.k21 = q1 / v2;
      }
      if (q2nd > 0.0) {
        p.k13 = q2nd / v1;
        if (v3 > 0.0) p.k31 = q2nd / v3;
      }
    }
    return;
  }
  if (advan == 12 && oral) {
    const double vc = pk_oral_central_volume(p, advan, 4);
    const double vp1 = pk_peripheral_volume(p, 1);
    const double vp2 = pk_peripheral_volume(p, 2);
    const double q_a = p.q3 > 0.0 ? p.q3 : q1;
    const double q_b = p.q4 > 0.0 ? p.q4 : q2nd;
    if (p.cl > 0.0 && vc > 0.0) {
      set_elimination_rate(p, p.cl / vc, true);
      if (q_a > 0.0) {
        p.k23 = q_a / vc;
        if (vp1 > 0.0) p.k32 = q_a / vp1;
      }
      if (q_b > 0.0) {
        p.k24 = q_b / vc;
        if (vp2 > 0.0) p.k42 = q_b / vp2;
      }
    }
  }
}

inline double nm_cardano_cbrt(double x) {
  if (x >= 0.0) return std::pow(x, 1.0 / 3.0);
  return -std::pow(-x, 1.0 / 3.0);
}

inline bool nm_cardano_3_roots(double tsum, double tb, double tc,
                               double& alpha, double& beta, double& gamma) {
  const double tq = tb / 3.0 - tsum * tsum / 9.0;
  const double tr = (tsum * tsum * tsum - 9.0 * tsum * tb + 27.0 * tc) / 108.0;
  const double td = tq * tq * tq + tr * tr;
  const double c1tr = nm_cardano_cbrt(tr + td);
  const double c2tr = nm_cardano_cbrt(tr - td);
  constexpr double rt3_2 = 0.866025403784;
  alpha = -tsum / 3.0 + (c1tr + c2tr) / 2.0;
  beta = -tsum / 3.0 + (c1tr * (-0.5 + rt3_2) + c2tr * (-0.5 - rt3_2)) / 2.0;
  gamma = -tsum / 3.0 + (c1tr * (-0.5 - rt3_2) + c2tr * (-0.5 + rt3_2)) / 2.0;
  return alpha > 0.0 && beta > 0.0 && gamma > 0.0;
}

inline void synthesize_hybrid_2_from_macro(double k10m, double k12m, double k21m,
                                            PkParams& p) {
  const double e1 = k10m + k12m;
  const double e2 = k21m;
  const double disc = (e1 - e2) * (e1 - e2) + 4.0 * k12m * k21m;
  if (disc < 0.0) return;
  p.pk_alpha = 0.5 * (e1 + e2 + std::sqrt(disc));
  p.pk_beta = 0.5 * (e1 + e2 - std::sqrt(disc));
}

inline void synthesize_trans6_hybrid(int advan, bool oral, PkParams& p) {
  if (advan == 11 && !oral &&
      p.k10 > 0.0 && p.k12 > 0.0 && p.k13 > 0.0 && p.k21 > 0.0 && p.k31 > 0.0) {
    const double e1 = p.k10 + p.k12 + p.k13;
    const double e2 = p.k21 + p.k31;
    const double e3 = p.k12 + p.k13;
    const double tsum = e1 + e2 + e3;
    const double tb = e1 * e2 + e1 * e3 + e2 * e3 - p.k12 * p.k21 - p.k13 * p.k31;
    const double tc = e1 * e2 * e3 - e3 * p.k12 * p.k21 - e2 * p.k13 * p.k31;
    nm_cardano_3_roots(tsum, tb, tc, p.pk_alpha, p.pk_beta, p.pk_gamma);
    return;
  }
  if (advan == 12 && oral &&
      p.k23 > 0.0 && p.k24 > 0.0 && p.k32 > 0.0 && p.k42 > 0.0) {
    const double k10m = p.k20 > 0.0 ? p.k20 : p.k10;
    if (!(k10m > 0.0)) return;
    const double e1 = k10m + p.k23 + p.k24;
    const double e2 = p.k32 + p.k42;
    const double e3 = p.k23 + p.k24;
    const double tsum = e1 + e2 + e3;
    const double tb = e1 * e2 + e1 * e3 + e2 * e3 - p.k23 * p.k32 - p.k24 * p.k42;
    const double tc = e1 * e2 * e3 - e3 * p.k23 * p.k32 - e2 * p.k24 * p.k42;
    nm_cardano_3_roots(tsum, tb, tc, p.pk_alpha, p.pk_beta, p.pk_gamma);
  }
}

inline void clear_trans5_oral_rates(PkParams& p) {
  p.k32 = 0.0;
  p.k23 = 0.0;
  p.k20 = 0.0;
  p.k10 = 0.0;
}

inline void clear_trans5_iv_rates(PkParams& p) {
  p.k21 = 0.0;
  p.k12 = 0.0;
  p.k10 = 0.0;
}

inline void apply_trans5(int advan, bool oral, PkParams& p) {
  if (!(p.pk_alpha > 0.0 && p.pk_beta > 0.0)) return;
  const double denom = p.aob + 1.0;
  if (advan == 3 && !oral) {
    if (p.k21 <= 0.0 && denom > 0.0) {
      p.k21 = (p.aob * p.pk_beta + p.pk_alpha) / denom;
    }
    if (p.k21 > 0.0) {
      p.k10 = (p.pk_alpha * p.pk_beta) / p.k21;
      p.k12 = p.pk_alpha + p.pk_beta - p.k21 - p.k10;
    }
    return;
  }
  if (advan == 4 && oral) {
    if (p.k32 <= 0.0 && denom > 0.0) {
      p.k32 = (p.aob * p.pk_beta + p.pk_alpha) / denom;
    }
    if (p.k32 > 0.0) {
      const double k = (p.pk_alpha * p.pk_beta) / p.k32;
      set_elimination_rate(p, k, true);
      p.k23 = p.pk_alpha + p.pk_beta - p.k32 - k;
    }
  }
}

inline void apply_trans6(int advan, bool oral, PkParams& p) {
  if (advan == 3 && !oral) {
    if (p.pk_alpha > 0.0 && p.pk_beta > 0.0 && p.k21 > 0.0) {
      p.k10 = (p.pk_alpha * p.pk_beta) / p.k21;
      p.k12 = p.pk_alpha + p.pk_beta - p.k21 - p.k10;
    }
    return;
  }
  if (advan == 4 && oral) {
    if (p.pk_alpha > 0.0 && p.pk_beta > 0.0 && p.k32 > 0.0) {
      const double k = (p.pk_alpha * p.pk_beta) / p.k32;
      set_elimination_rate(p, k, true);
      p.k23 = p.pk_alpha + p.pk_beta - p.k32 - k;
    }
    return;
  }
  if (advan == 11 && !oral) {
    if (p.pk_alpha > 0.0 && p.pk_beta > 0.0 && p.pk_gamma > 0.0 &&
        p.k21 > 0.0 && p.k31 > 0.0) {
      p.k10 = (p.pk_alpha * p.pk_beta * p.pk_gamma) / (p.k21 * p.k31);
      const double h_v1 = p.pk_alpha + p.pk_beta + p.pk_gamma;
      const double h_v2 = p.pk_alpha * p.pk_beta + p.pk_alpha * p.pk_gamma +
                          p.pk_beta * p.pk_gamma;
      const double den = p.k21 - p.k31;
      if (std::fabs(den) > 1e-12) {
        p.k13 = (h_v2 + p.k31 * p.k31 - p.k31 * h_v1 - p.k10 * p.k21) / den;
      }
      p.k12 = h_v1 - p.k10 - p.k13 - p.k21 - p.k31;
    }
    return;
  }
  if (advan == 12 && oral) {
    if (p.pk_alpha > 0.0 && p.pk_beta > 0.0 && p.pk_gamma > 0.0 &&
        p.k32 > 0.0 && p.k42 > 0.0) {
      const double k = (p.pk_alpha * p.pk_beta * p.pk_gamma) / (p.k32 * p.k42);
      set_elimination_rate(p, k, true);
      const double h_v2 = p.pk_alpha + p.pk_beta + p.pk_gamma;
      const double h_v3 = p.pk_alpha * p.pk_beta + p.pk_alpha * p.pk_gamma +
                          p.pk_beta * p.pk_gamma;
      const double den = p.k32 - p.k42;
      if (std::fabs(den) > 1e-12) {
        p.k24 = (h_v3 + p.k42 * p.k42 - p.k42 * h_v2 - k * p.k32) / den;
      }
      p.k23 = h_v2 - k - p.k24 - p.k32 - p.k42;
    }
  }
}

inline void finalize_volume_aliases(PkParams& p) {
  if (p.v1 <= 0.0 && p.vc > 0.0) p.v1 = p.vc;
  if (p.v2 <= 0.0 && p.vp > 0.0) p.v2 = p.vp;
  if (p.v3 <= 0.0 && p.vp2 > 0.0) p.v3 = p.vp2;
  if (p.vc <= 0.0 && p.v1 > 0.0) p.vc = p.v1;
}

inline void apply_trans(int trans, int advan, RouteType route, PkParams& p) {
  normalize_volume_aliases(p);
  ensure_vss(p);
  const bool oral = (route == ROUTE_ORAL || route == ROUTE_TRANSIT);
  const bool micro_trans1 = (trans == 1 && (advan == 3 || advan == 4 || advan == 11 || advan == 12));

  if (trans == 2) {
    apply_trans2(advan, oral, p);
  } else if (trans == 3) {
    apply_trans3(advan, oral, p);
  } else if (trans == 4) {
    apply_trans4(advan, oral, p);
  } else if (trans == 5) {
    if (!(p.pk_alpha > 0.0 && p.pk_beta > 0.0)) {
      apply_trans4(advan, oral, p);
      if (advan == 4 && oral) {
        const double k10m = p.k20 > 0.0 ? p.k20 : p.k10;
        if (k10m > 0.0 && p.k23 > 0.0 && p.k32 > 0.0) {
          synthesize_hybrid_2_from_macro(k10m, p.k23, p.k32, p);
        }
      } else if (advan == 3 && !oral && p.k10 > 0.0 && p.k12 > 0.0 && p.k21 > 0.0) {
        synthesize_hybrid_2_from_macro(p.k10, p.k12, p.k21, p);
      } else if (p.k10 > 0.0 && p.k12 > 0.0 && p.k21 > 0.0) {
        synthesize_hybrid_2_from_macro(p.k10, p.k12, p.k21, p);
      }
    }
    if (p.aob <= 0.0) p.aob = 1.0;
    if (advan == 4 && oral) clear_trans5_oral_rates(p);
    if (advan == 3 && !oral) clear_trans5_iv_rates(p);
    apply_trans5(advan, oral, p);
  } else if (trans == 6) {
    if (!(p.pk_alpha > 0.0 && p.pk_beta > 0.0 && p.pk_gamma > 0.0)) {
      synthesize_trans6_hybrid(advan, oral, p);
    }
    apply_trans6(advan, oral, p);
  } else if (!micro_trans1 && p.cl > 0.0) {
    apply_trans2(advan, oral, p);
  }

  if (oral && advan == 12 && (trans == 1 || (trans == 6 && p.pk_alpha <= 0.0))) {
    apply_trans1_oral_aliases(advan, p);
  } else if (oral && advan == 4 && trans == 6 && p.pk_alpha <= 0.0 &&
             p.k20 <= 0.0 && p.k10 > 0.0) {
    p.k20 = p.k10;
  } else if (oral && trans == 1 && !micro_trans1 && advan != 12) {
    apply_trans1_oral_aliases(advan, p);
  }

  // ADVAN4 TRANS1 micro: user $PK uses K10; oral central elimination is k20.
  if (oral && advan == 4 && trans == 1 && p.k20 <= 0.0 && p.k10 > 0.0) {
    p.k20 = p.k10;
  }

  // ADVAN6/13 ODE models use CL/V in $PK; SS init needs k10/k20 even when TRANS>1.
  if ((advan == 6 || advan == 13) && p.cl > 0.0 &&
      p.k10 <= 0.0 && p.k20 <= 0.0) {
    apply_trans2(advan, oral, p);
  }

  if (route == ROUTE_TRANSIT) {
    if (p.k20 <= 0.0 && p.cl > 0.0 && p.vc > 0.0) p.k20 = p.cl / p.vc;
    if (p.k20 <= 0.0 && p.k10 > 0.0) p.k20 = p.k10;
  }
  finalize_volume_aliases(p);
}

struct SubjectEvents {
  Rcpp::NumericVector time;
  Rcpp::NumericVector amt;
  Rcpp::NumericVector rate;
  Rcpp::NumericVector f1;
  Rcpp::IntegerVector cmt;
  Rcpp::IntegerVector evid;
  Rcpp::IntegerVector ss;
  Rcpp::NumericVector ii;
  Rcpp::NumericVector s1;
  Rcpp::NumericVector s2;
  Rcpp::NumericVector s3;
  Rcpp::NumericVector s4;
  Rcpp::NumericMatrix scale_mat;
  bool use_data_scale = false;
  Rcpp::NumericMatrix f_mat;
  bool use_data_f = false;
};

inline double pk_f_param(int cmt, const PkParams& p) {
  if (cmt < 1 || cmt > kMaxScale) return 0.0;
  return p.f[static_cast<size_t>(cmt - 1)];
}

inline double row_f_at(const SubjectEvents& ev, int i, int cmt) {
  if (i < 0 || cmt < 1 || cmt > kMaxScale) return 0.0;
  if (ev.f_mat.nrow() > 0 &&
      ev.f_mat.ncol() >= cmt &&
      i < ev.f_mat.nrow()) {
    const double v = ev.f_mat(i, cmt - 1);
    if (R_finite(v) && v > 0.0) return v;
  }
  if (cmt == 1 && ev.f1.size() > static_cast<R_xlen_t>(i) && ev.f1[i] > 0.0) {
    return ev.f1[i];
  }
  return 0.0;
}

inline double effective_f(int cmt, const PkParams& p, const SubjectEvents& ev, int i) {
  const double row_f = row_f_at(ev, i, cmt);
  const double pk_f = pk_f_param(cmt, p);
  if (ev.use_data_f && row_f > 0.0) return row_f;
  if (pk_f > 0.0) return pk_f;
  return 1.0;
}

inline double effective_rate(const SubjectEvents& ev, int i) {
  if (ev.rate.size() == 0 || i < 0 || i >= ev.rate.size()) return 0.0;
  const double r = ev.rate[i];
  if (r > 0.0) return r;
  return 0.0;
}

Rcpp::NumericVector solve_1comp_transit(
    const SubjectEvents& ev, const PkParams& p,
    int n_transit, int obs_cmp, int model_ss);

Rcpp::NumericVector solve_2comp_transit(
    const SubjectEvents& ev, const PkParams& p,
    int n_transit, int obs_cmp, int model_ss);

Rcpp::NumericVector solve_metab_route(
    const SubjectEvents& ev, const PkParams& p,
    int advan, RouteType route, int obs_cmp, int model_ss);

}  // namespace nm_pk

Rcpp::NumericVector nm_pk_route_cpp(
    int advan, int trans, int obs_cmp, int dose_cmp,
    int n_transit, bool use_ode, int model_ss,
    nm_pk::SubjectEvents ev, nm_pk::PkParams params);

Rcpp::List nm_pk_route_detail_cpp(
    int advan, int trans, int obs_cmp, int dose_cmp,
    int n_transit, bool use_ode, int model_ss,
    nm_pk::SubjectEvents ev, nm_pk::PkParams params, int n_state);

Rcpp::NumericVector nm_pk_route_r(
    int advan, int trans, int obs_cmp, int dose_cmp,
    int n_transit, bool use_ode, int model_ss,
    Rcpp::NumericVector time, Rcpp::NumericVector amt, Rcpp::NumericVector rate,
    Rcpp::NumericVector f1, Rcpp::IntegerVector cmt, Rcpp::IntegerVector evid,
    Rcpp::IntegerVector ss, Rcpp::NumericVector ii,
    Rcpp::List pk_params,
    Rcpp::NumericVector s1,
    Rcpp::NumericVector s2,
    Rcpp::NumericVector s3,
    Rcpp::NumericVector s4,
    Rcpp::NumericMatrix scale_mat,
    bool use_data_scale,
    Rcpp::NumericMatrix f_mat,
    bool use_data_f);

Rcpp::List nm_eval_pred_cpp(
    Rcpp::CharacterVector pred_lines,
    Rcpp::NumericVector theta,
    Rcpp::NumericVector eta,
    Rcpp::List covariates,
    Rcpp::CharacterVector des_lines);

bool nm_cpp_advan_supported(int advan, int trans);

bool nm_pred_expr_check_cpp(Rcpp::CharacterVector pred_lines);

