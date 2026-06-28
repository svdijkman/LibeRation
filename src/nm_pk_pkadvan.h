#pragma once
// Analytical 1–3 compartment PK step functions ported from PKADVAN
// (https://github.com/abuhelwa/PKADVAN_Rpackage, ADVAN-style non-superpositioning).
#include "nm_pk_core.h"
#include "nm_pk_ss.h"
#include <array>
#include <cmath>

namespace nm_pk {
namespace pkadvan {

inline void step_1_iv_bolus(double dt, double k10, double prev_a1, double& a1) {
  a1 = prev_a1 * std::exp(-dt * k10);
}

inline void step_1_iv_infusion(double dt, double k10, double rate, double prev_a1, double& a1) {
  a1 = rate / k10 * (1.0 - std::exp(-dt * k10)) + prev_a1 * std::exp(-dt * k10);
}

inline void step_1_oral(double dt, double ka, double k10,
                        double prev_gut, double prev_cen,
                        double& gut, double& cen) {
  cen = prev_gut * lhopital_exp_diff(ka, k10, dt) + prev_cen * std::exp(-dt * k10);
  gut = prev_gut * std::exp(-dt * ka);
}

inline void step_2_iv_bolus(double dt, double k10, double k12, double k21,
                            double prev_a1, double prev_a2,
                            double& a1, double& a2) {
  Hybrid2 h = hybrid_lambdas_2(k10, k12, k21);
  const double d = h.lambda2 - h.lambda1;
  if (std::fabs(d) < 1e-12) {
    const double e = std::exp(-h.lambda1 * dt);
    a1 = (prev_a1 + prev_a2 * k21 / h.lambda1) * e;
    a2 = prev_a2 * e;
    return;
  }
  a1 = (((prev_a1 * h.E2 + prev_a2 * k21) - prev_a1 * h.lambda1) * std::exp(-h.lambda1 * dt) -
        ((prev_a1 * h.E2 + prev_a2 * k21) - prev_a1 * h.lambda2) * std::exp(-h.lambda2 * dt)) / d;
  a2 = (((prev_a2 * h.E1 + prev_a1 * k12) - prev_a2 * h.lambda1) * std::exp(-h.lambda1 * dt) -
        ((prev_a2 * h.E1 + prev_a1 * k12) - prev_a2 * h.lambda2) * std::exp(-h.lambda2 * dt)) / d;
}

inline void step_2_iv_infusion(double dt, double k10, double k12, double k21, double rate,
                               double prev_a1, double prev_a2,
                               double& a1, double& a2) {
  Hybrid2 h = hybrid_lambdas_2(k10, k12, k21);
  const double d = h.lambda2 - h.lambda1;
  if (std::fabs(d) < 1e-12) {
    a1 = (prev_a1 + prev_a2 * k21 / h.lambda1 + rate) * std::exp(-h.lambda1 * dt);
    a2 = prev_a2 * std::exp(-h.lambda1 * dt);
    return;
  }
  a1 = (((prev_a1 * h.E2 + rate + prev_a2 * k21) - prev_a1 * h.lambda1) * std::exp(-h.lambda1 * dt) -
        ((prev_a1 * h.E2 + rate + prev_a2 * k21) - prev_a1 * h.lambda2) * std::exp(-h.lambda2 * dt)) / d;
  a1 += rate * h.E2 * (1.0 / (h.lambda1 * h.lambda2) +
         std::exp(-h.lambda1 * dt) / (h.lambda1 * (h.lambda1 - h.lambda2)) -
         std::exp(-h.lambda2 * dt) / (h.lambda2 * (h.lambda1 - h.lambda2)));
  a2 = (((prev_a2 * h.E1 + prev_a1 * k12) - prev_a2 * h.lambda1) * std::exp(-h.lambda1 * dt) -
        ((prev_a2 * h.E1 + prev_a1 * k12) - prev_a2 * h.lambda2) * std::exp(-h.lambda2 * dt)) / d;
  a2 += rate * k12 * (1.0 / (h.lambda1 * h.lambda2) +
         std::exp(-h.lambda1 * dt) / (h.lambda1 * (h.lambda1 - h.lambda2)) -
         std::exp(-h.lambda2 * dt) / (h.lambda2 * (h.lambda1 - h.lambda2)));
}

inline void step_2_oral(double dt, double ka, double k20, double k23, double k32,
                        double prev_gut, double prev_a2, double prev_a3,
                        double& gut, double& a2, double& a3) {
  Hybrid2 h = hybrid_lambdas_2(k20, k23, k32);
  const double E2 = h.E1;
  const double E3 = h.E2;
  const double d = h.lambda2 - h.lambda1;
  if (std::fabs(d) < 1e-12) {
    a2 = prev_a2 * std::exp(-h.lambda1 * dt);
    a3 = prev_a3 * std::exp(-h.lambda1 * dt);
  } else {
    a2 = (((prev_a2 * E3 + prev_a3 * k32) - prev_a2 * h.lambda1) * std::exp(-h.lambda1 * dt) -
          ((prev_a2 * E3 + prev_a3 * k32) - prev_a2 * h.lambda2) * std::exp(-h.lambda2 * dt)) / d;
    a3 = (((prev_a3 * E2 + prev_a2 * k23) - prev_a3 * h.lambda1) * std::exp(-h.lambda1 * dt) -
          ((prev_a3 * E2 + prev_a2 * k23) - prev_a3 * h.lambda2) * std::exp(-h.lambda2 * dt)) / d;
  }
  a2 += prev_gut * ka * (
      std::exp(-dt * ka) * (E3 - ka) / ((h.lambda1 - ka) * (h.lambda2 - ka)) +
      std::exp(-dt * h.lambda1) * (E3 - h.lambda1) / ((h.lambda2 - h.lambda1) * (ka - h.lambda1)) +
      std::exp(-dt * h.lambda2) * (E3 - h.lambda2) / ((h.lambda1 - h.lambda2) * (ka - h.lambda2)));
  a3 += prev_gut * ka * k23 * (
      std::exp(-dt * ka) / ((h.lambda1 - ka) * (h.lambda2 - ka)) +
      std::exp(-dt * h.lambda1) / ((h.lambda2 - h.lambda1) * (ka - h.lambda1)) +
      std::exp(-dt * h.lambda2) / ((h.lambda1 - h.lambda2) * (ka - h.lambda2)));
  gut = prev_gut * std::exp(-dt * ka);
}

inline void step_oral4_chain_rk4(double dt, double ka, double k10,
                                 double k12, double k21, double k23, double k32,
                                 double& gut, double& cen, double& p1, double& p2) {
  const int steps = std::max(1, static_cast<int>(std::ceil(dt * 20.0)));
  const double h = dt / steps;
  auto rhs = [](const std::array<double, 4>& y, std::array<double, 4>& dy,
                double ka0, double k100, double k120, double k210,
                double k230, double k320) {
    dy[0] = -ka0 * y[0];
    dy[1] = ka0 * y[0] - (k100 + k120) * y[1] + k210 * y[2];
    dy[2] = k120 * y[1] - (k210 + k230) * y[2] + k320 * y[3];
    dy[3] = k230 * y[2] - k320 * y[3];
  };
  for (int s = 0; s < steps; ++s) {
    std::array<double, 4> y{gut, cen, p1, p2};
    std::array<double, 4> k1{}, k2{}, k3{}, k4{}, tmp{};
    rhs(y, k1, ka, k10, k12, k21, k23, k32);
    tmp = {y[0] + 0.5 * h * k1[0], y[1] + 0.5 * h * k1[1],
           y[2] + 0.5 * h * k1[2], y[3] + 0.5 * h * k1[3]};
    rhs(tmp, k2, ka, k10, k12, k21, k23, k32);
    tmp = {y[0] + 0.5 * h * k2[0], y[1] + 0.5 * h * k2[1],
           y[2] + 0.5 * h * k2[2], y[3] + 0.5 * h * k2[3]};
    rhs(tmp, k3, ka, k10, k12, k21, k23, k32);
    tmp = {y[0] + h * k3[0], y[1] + h * k3[1], y[2] + h * k3[2], y[3] + h * k3[3]};
    rhs(tmp, k4, ka, k10, k12, k21, k23, k32);
    gut = y[0] + h / 6.0 * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]);
    cen = y[1] + h / 6.0 * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]);
    p1 = y[2] + h / 6.0 * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]);
    p2 = y[3] + h / 6.0 * (k1[3] + 2.0 * k2[3] + 2.0 * k3[3] + k4[3]);
  }
}

inline void step_3_iv_bolus(double dt, double k10, double k12, double k21, double k13, double k31,
                            double prev_a1, double prev_a2, double prev_a3,
                            double& a1, double& a2, double& a3) {
  Hybrid3 hy = hybrid_lambdas_3_iv(k10, k12, k21, k13, k31);
  const double E1 = hy.E1, E2 = hy.E2, E3 = hy.E3;
  const double l1 = hy.lambda1, l2 = hy.lambda2, l3 = hy.lambda3;
  const double B = prev_a2 * k21 + prev_a3 * k31;
  const double C = E3 * prev_a2 * k21 + E2 * prev_a3 * k31;
  const double I = prev_a1 * k12 * E3 - prev_a2 * k13 * k31 + prev_a3 * k12 * k31;
  const double J = prev_a1 * k13 * E2 + prev_a2 * k13 * k21 - prev_a3 * k12 * k21;
  a1 = prev_a1 * (std::exp(-dt * l1) * (E2 - l1) * (E3 - l1) / ((l2 - l1) * (l3 - l1)) +
                   std::exp(-dt * l2) * (E2 - l2) * (E3 - l2) / ((l1 - l2) * (l3 - l2)) +
                   std::exp(-dt * l3) * (E2 - l3) * (E3 - l3) / ((l1 - l3) * (l2 - l3)));
  a1 += std::exp(-dt * l1) * (C - B * l1) / ((l1 - l2) * (l1 - l3));
  a1 += std::exp(-dt * l2) * (B * l2 - C) / ((l1 - l2) * (l2 - l3));
  a1 += std::exp(-dt * l3) * (B * l3 - C) / ((l1 - l3) * (l3 - l2));
  a2 = prev_a2 * (std::exp(-dt * l1) * (E1 - l1) * (E3 - l1) / ((l2 - l1) * (l3 - l1)) +
                  std::exp(-dt * l2) * (E1 - l2) * (E3 - l2) / ((l1 - l2) * (l3 - l2)) +
                  std::exp(-dt * l3) * (E1 - l3) * (E3 - l3) / ((l1 - l3) * (l2 - l3)));
  a2 += std::exp(-dt * l1) * (I - prev_a1 * k12 * l1) / ((l1 - l2) * (l1 - l3));
  a2 += std::exp(-dt * l2) * (prev_a1 * k12 * l2 - I) / ((l1 - l2) * (l2 - l3));
  a2 += std::exp(-dt * l3) * (prev_a1 * k12 * l3 - I) / ((l1 - l3) * (l3 - l2));
  a3 = prev_a3 * (std::exp(-dt * l1) * (E1 - l1) * (E2 - l1) / ((l2 - l1) * (l3 - l1)) +
                  std::exp(-dt * l2) * (E1 - l2) * (E2 - l2) / ((l1 - l2) * (l3 - l2)) +
                  std::exp(-dt * l3) * (E1 - l3) * (E2 - l3) / ((l1 - l3) * (l2 - l3)));
  a3 += std::exp(-dt * l1) * (J - prev_a1 * k13 * l1) / ((l1 - l2) * (l1 - l3));
  a3 += std::exp(-dt * l2) * (prev_a1 * k13 * l2 - J) / ((l1 - l2) * (l2 - l3));
  a3 += std::exp(-dt * l3) * (prev_a1 * k13 * l3 - J) / ((l1 - l3) * (l3 - l2));
}

inline void step_3_iv_infusion(double dt, double k10, double k12, double k21, double k13, double k31,
                               double rate,
                               double prev_a1, double prev_a2, double prev_a3,
                               double& a1, double& a2, double& a3) {
  step_3_iv_bolus(dt, k10, k12, k21, k13, k31, prev_a1, prev_a2, prev_a3, a1, a2, a3);
  Hybrid3 hy = hybrid_lambdas_3_iv(k10, k12, k21, k13, k31);
  const double E2 = hy.E2, E3 = hy.E3;
  const double l1 = hy.lambda1, l2 = hy.lambda2, l3 = hy.lambda3;
  a1 += rate * ((E2 * E3) / (l1 * l2 * l3) -
                std::exp(-dt * l1) * (E2 - l1) * (E3 - l1) / (l1 * (l2 - l1) * (l3 - l1)) -
                std::exp(-dt * l2) * (E2 - l2) * (E3 - l2) / (l2 * (l1 - l2) * (l3 - l2)) -
                std::exp(-dt * l3) * (E2 - l3) * (E3 - l3) / (l3 * (l1 - l3) * (l2 - l3)));
  a2 += rate * k12 * (E3 / (l1 * l2 * l3) -
                      std::exp(-dt * l1) * (E3 - l1) / (l1 * (l2 - l1) * (l3 - l1)) -
                      std::exp(-dt * l2) * (E3 - l2) / (l2 * (l1 - l2) * (l3 - l2)) -
                      std::exp(-dt * l3) * (E3 - l3) / (l3 * (l1 - l3) * (l2 - l3)));
  a3 += rate * k13 * (E2 / (l1 * l2 * l3) -
                      std::exp(-dt * l1) * (E2 - l1) / (l1 * (l2 - l1) * (l3 - l1)) -
                      std::exp(-dt * l2) * (E2 - l2) / (l2 * (l1 - l2) * (l3 - l2)) -
                      std::exp(-dt * l3) * (E2 - l3) / (l3 * (l1 - l3) * (l2 - l3)));
}

inline void step_3_oral(double dt, double ka, double k20, double k23, double k32, double k24, double k42,
                        double prev_gut, double prev_a2, double prev_a3, double prev_a4,
                        double& gut, double& a2, double& a3, double& a4) {
  Hybrid3 hy = hybrid_lambdas_3_oral(k20, k23, k32, k24, k42);
  const double E2 = hy.E1, E3 = hy.E2, E4 = hy.E3;
  const double l1 = hy.lambda1, l2 = hy.lambda2, l3 = hy.lambda3;
  const double B = prev_a3 * k32 + prev_a4 * k42;
  const double C = E4 * prev_a3 * k32 + E3 * prev_a4 * k42;
  const double I = prev_a2 * k23 * E4 - prev_a3 * k24 * k42 + prev_a4 * k23 * k42;
  const double J = prev_a2 * k24 * E3 + prev_a3 * k24 * k32 - prev_a4 * k23 * k32;

  a2 = prev_a2 * (std::exp(-dt * l1) * (E3 - l1) * (E4 - l1) / ((l2 - l1) * (l3 - l1)) +
                  std::exp(-dt * l2) * (E3 - l2) * (E4 - l2) / ((l1 - l2) * (l3 - l2)) +
                  std::exp(-dt * l3) * (E3 - l3) * (E4 - l3) / ((l1 - l3) * (l2 - l3)));
  a2 += std::exp(-dt * l1) * (C - B * l1) / ((l1 - l2) * (l1 - l3));
  a2 += std::exp(-dt * l2) * (B * l2 - C) / ((l1 - l2) * (l2 - l3));
  a2 += std::exp(-dt * l3) * (B * l3 - C) / ((l1 - l3) * (l3 - l2));
  a2 += prev_gut * ka * (
      std::exp(-dt * l1) * (E3 - l1) * (E4 - l1) / ((l2 - l1) * (l3 - l1) * (ka - l1)) +
      std::exp(-dt * l2) * (E3 - l2) * (E4 - l2) / ((l1 - l2) * (l3 - l2) * (ka - l2)) +
      std::exp(-dt * l3) * (E3 - l3) * (E4 - l3) / ((l1 - l3) * (l2 - l3) * (ka - l3)) +
      std::exp(-dt * ka) * (E3 - ka) * (E4 - ka) / ((l1 - ka) * (l2 - ka) * (l3 - ka)));

  a3 = prev_a3 * (std::exp(-dt * l1) * (E2 - l1) * (E4 - l1) / ((l2 - l1) * (l3 - l1)) +
                  std::exp(-dt * l2) * (E2 - l2) * (E4 - l2) / ((l1 - l2) * (l3 - l2)) +
                  std::exp(-dt * l3) * (E2 - l3) * (E4 - l3) / ((l1 - l3) * (l2 - l3)));
  a3 += std::exp(-dt * l1) * (I - prev_a2 * k23 * l1) / ((l1 - l2) * (l1 - l3));
  a3 += std::exp(-dt * l2) * (prev_a2 * k23 * l2 - I) / ((l1 - l2) * (l2 - l3));
  a3 += std::exp(-dt * l3) * (prev_a2 * k23 * l3 - I) / ((l1 - l3) * (l3 - l2));
  a3 += prev_gut * ka * k23 * (
      std::exp(-dt * l1) * (E4 - l1) / ((l2 - l1) * (l3 - l1) * (ka - l1)) +
      std::exp(-dt * l2) * (E4 - l2) / ((l1 - l2) * (l3 - l2) * (ka - l2)) +
      std::exp(-dt * l3) * (E4 - l3) / ((l1 - l3) * (l2 - l3) * (ka - l3)) +
      std::exp(-dt * ka) * (E4 - ka) / ((l1 - ka) * (l2 - ka) * (l3 - ka)));

  a4 = prev_a4 * (std::exp(-dt * l1) * (E2 - l1) * (E3 - l1) / ((l2 - l1) * (l3 - l1)) +
                  std::exp(-dt * l2) * (E2 - l2) * (E3 - l2) / ((l1 - l2) * (l3 - l2)) +
                  std::exp(-dt * l3) * (E2 - l3) * (E3 - l3) / ((l1 - l3) * (l2 - l3)));
  a4 += std::exp(-dt * l1) * (J - prev_a2 * k24 * l1) / ((l1 - l2) * (l1 - l3));
  a4 += std::exp(-dt * l2) * (prev_a2 * k24 * l2 - J) / ((l1 - l2) * (l2 - l3));
  a4 += std::exp(-dt * l3) * (prev_a2 * k24 * l3 - J) / ((l1 - l3) * (l3 - l2));
  a4 += prev_gut * ka * k24 * (
      std::exp(-dt * l1) * (E3 - l1) / ((l2 - l1) * (l3 - l1) * (ka - l1)) +
      std::exp(-dt * l2) * (E3 - l2) / ((l1 - l2) * (l3 - l2) * (ka - l2)) +
      std::exp(-dt * l3) * (E3 - l3) / ((l1 - l3) * (l2 - l3) * (ka - l3)) +
      std::exp(-dt * ka) * (E3 - ka) / ((l1 - ka) * (l2 - ka) * (l3 - ka)));

  gut = prev_gut * std::exp(-dt * ka);
}

inline double row_f1(const SubjectEvents& ev, int i, const PkParams& p, int default_cmp = 1) {
  return effective_f(row_cmt(ev, i, default_cmp), p, ev, i);
}

inline void apply_iv_bolus_dose(const SubjectEvents& ev, int i, const PkParams& p,
                                double& a1, double& a2, double& a3) {
  double dose = ev.amt[i];
  if (dose <= 0.0) return;
  const int cmt = row_cmt(ev, i, pk_route_dose_cmp);
  dose *= effective_f(cmt, p, ev, i);
  if (cmt == 2) a2 += dose;
  else if (cmt == 3) a3 += dose;
  else a1 += dose;
}

inline void apply_oral_dose(const SubjectEvents& ev, int i, const PkParams& p,
                            std::vector<double>& amounts) {
  double dose = ev.amt[i];
  if (dose <= 0.0) return;
  const int cmt = row_cmt(ev, i, pk_route_dose_cmp);
  dose *= effective_f(cmt, p, ev, i);
  const int idx = cmt - 1;
  if (idx >= 0 && idx < static_cast<int>(amounts.size())) {
    amounts[static_cast<size_t>(idx)] += dose;
  } else if (!amounts.empty()) {
    amounts[0] += dose;
  }
}

inline double pk_scale_param(int cmt, const PkParams& p) {
  if (cmt < 1 || cmt > kMaxScale) return 0.0;
  return p.scale[static_cast<size_t>(cmt - 1)];
}

inline double row_scale_at(const SubjectEvents& ev, int i, int cmt) {
  if (i < 0 || cmt < 1 || cmt > kMaxScale) return 0.0;
  if (ev.scale_mat.nrow() > 0 &&
      ev.scale_mat.ncol() >= cmt &&
      i < ev.scale_mat.nrow()) {
    const double v = ev.scale_mat(i, cmt - 1);
    if (R_finite(v) && v > 0.0) return v;
  }
  switch (cmt) {
    case 1:
      return (ev.s1.size() > static_cast<R_xlen_t>(i)) ? ev.s1[i] : 0.0;
    case 2:
      return (ev.s2.size() > static_cast<R_xlen_t>(i)) ? ev.s2[i] : 0.0;
    case 3:
      return (ev.s3.size() > static_cast<R_xlen_t>(i)) ? ev.s3[i] : 0.0;
    case 4:
      return (ev.s4.size() > static_cast<R_xlen_t>(i)) ? ev.s4[i] : 0.0;
    default: return 0.0;
  }
}

inline double legacy_conc_iv(int obs_cmp, double a1, double a2, double a3, const PkParams& p) {
  if (obs_cmp == 2) return a2 / std::max(p.v2, 1e-8);
  if (obs_cmp == 3) return a3 / std::max(p.v3, 1e-8);
  return a1 / std::max(p.v1, 1e-8);
}

inline double legacy_conc_oral(int obs_cmp, double gut, double a2, double a3, double a4, const PkParams& p) {
  const double vc = p.vc > 0.0 ? p.vc : (p.v2 > 0.0 ? p.v2 : p.v1);
  if (obs_cmp == 1) return gut / std::max(vc, 1e-8);
  if (obs_cmp == 3) return a3 / std::max(p.v2 > 0.0 ? p.v2 : p.vp, 1e-8);
  if (obs_cmp == 4) return a4 / std::max(p.v3 > 0.0 ? p.v3 : p.vp2, 1e-8);
  return a2 / std::max(vc, 1e-8);
}

inline double effective_scale(int cmt, const PkParams& p, const SubjectEvents& ev, int i) {
  const double row_s = row_scale_at(ev, i, cmt);
  const double pk_s = pk_scale_param(cmt, p);
  if (ev.use_data_scale && row_s > 0.0) return row_s;
  if (pk_s > 0.0) return pk_s;
  return 0.0;
}

inline double amount_at(int obs_cmp, double a1, double a2, double a3, double a4 = 0.0) {
  if (obs_cmp == 2) return a2;
  if (obs_cmp == 3) return a3;
  if (obs_cmp == 4) return a4;
  return a1;
}

inline double conc_scaled(int obs_cmp, double amt, const PkParams& p,
                          const SubjectEvents& ev, int i) {
  const double s = effective_scale(obs_cmp, p, ev, i);
  if (s <= 0.0) return 0.0;
  return amt / std::max(s, 1e-8);
}

inline double conc_iv(int obs_cmp, double a1, double a2, double a3, const PkParams& p,
                      const SubjectEvents& ev = SubjectEvents(), int i = -1) {
  const double s = effective_scale(obs_cmp, p, ev, i);
  if (s <= 0.0) return legacy_conc_iv(obs_cmp, a1, a2, a3, p);
  return conc_scaled(obs_cmp, amount_at(obs_cmp, a1, a2, a3), p, ev, i);
}

inline double conc_oral(int obs_cmp, double gut, double a2, double a3, double a4, const PkParams& p,
                        const SubjectEvents& ev = SubjectEvents(), int i = -1) {
  const double s = effective_scale(obs_cmp, p, ev, i);
  if (s <= 0.0) return legacy_conc_oral(obs_cmp, gut, a2, a3, a4, p);
  const double amt = obs_cmp == 1 ? gut : amount_at(obs_cmp, gut, a2, a3, a4);
  return conc_scaled(obs_cmp, amt, p, ev, i);
}

inline void apply_dose_with_ss(int model_ss, int advan, RouteType route,
                               const SubjectEvents& ev, int i, const PkParams& p,
                               std::vector<double>& amounts, bool oral) {
  if (is_reset_row(ev, i)) {
    reset_amounts(amounts);
    return;
  }
  if (is_other_event_row(ev, i)) return;
  if (!is_dosing_row(ev, i)) return;
  if (apply_ss_at_dose(model_ss, row_ss(ev, i), advan, route, p, ev, i, amounts)) {
    return;
  }
  if (oral) {
    if (amounts.empty()) amounts.resize(1, 0.0);
    apply_oral_dose(ev, i, p, amounts);
  } else if (route != ROUTE_IV_INFUSION) {
    double a1 = amounts.size() > 0 ? amounts[0] : 0.0;
    double a2 = amounts.size() > 1 ? amounts[1] : 0.0;
    double a3 = amounts.size() > 2 ? amounts[2] : 0.0;
    apply_iv_bolus_dose(ev, i, p, a1, a2, a3);
    if (amounts.size() > 0) amounts[0] = a1;
    if (amounts.size() > 1) amounts[1] = a2;
    if (amounts.size() > 2) amounts[2] = a3;
  }
}

}  // namespace pkadvan
}  // namespace nm_pk
