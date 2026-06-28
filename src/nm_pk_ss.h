#pragma once
// Steady-state initial amounts via matrix exponential (SS_sol.R).
#include "nm_pk_core.h"
#include <vector>
#include <cmath>
#include <algorithm>

namespace nm_pk {

// Set in nm_pk_route_cpp before solving; used when row CMT is missing or zero.
extern thread_local int pk_route_dose_cmp;

namespace ss_mat {

inline const int MAXN = 4;

inline void mat_identity(int n, std::vector<double>& I) {
  I.assign(n * n, 0.0);
  for (int i = 0; i < n; ++i) I[i * n + i] = 1.0;
}

inline void mat_mul(const std::vector<double>& A, const std::vector<double>& B,
                    int n, std::vector<double>& C) {
  C.assign(n * n, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double s = 0.0;
      for (int k = 0; k < n; ++k) s += A[i * n + k] * B[k * n + j];
      C[i * n + j] = s;
    }
  }
}

inline void mat_vec(const std::vector<double>& A, const std::vector<double>& x,
                    int n, std::vector<double>& y) {
  y.assign(n, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) y[i] += A[i * n + j] * x[j];
  }
}

inline bool gauss_solve(std::vector<double> A, std::vector<double> b, int n,
                        std::vector<double>& x) {
  x = b;
  for (int col = 0; col < n; ++col) {
    int piv = col;
    double best = std::fabs(A[col * n + col]);
    for (int r = col + 1; r < n; ++r) {
      double v = std::fabs(A[r * n + col]);
      if (v > best) { best = v; piv = r; }
    }
    if (best < 1e-14) return false;
    if (piv != col) {
      for (int c = 0; c < n; ++c) std::swap(A[col * n + c], A[piv * n + c]);
      std::swap(x[col], x[piv]);
    }
    double d = A[col * n + col];
    for (int c = col; c < n; ++c) A[col * n + c] /= d;
    x[col] /= d;
    for (int r = 0; r < n; ++r) {
      if (r == col) continue;
      double f = A[r * n + col];
      for (int c = col; c < n; ++c) A[r * n + c] -= f * A[col * n + c];
      x[r] -= f * x[col];
    }
  }
  return true;
}

inline void mat_exp(const std::vector<double>& K, int n, double t,
                    std::vector<double>& Phi) {
  std::vector<double> A(n * n), B(n * n), C(n * n), T(n * n);
  for (int i = 0; i < n * n; ++i) A[i] = K[i] * t;
  double norm = 0.0;
  for (double v : A) norm = std::max(norm, std::fabs(v));
  int s = 0;
  while (norm > 0.5) {
    for (double& v : A) v *= 0.5;
    norm *= 0.5;
    ++s;
  }
  mat_identity(n, Phi);
  mat_identity(n, T);
  for (int k = 1; k <= 18; ++k) {
    mat_mul(T, A, n, C);
    for (double& v : C) v /= static_cast<double>(k);
    for (int i = 0; i < n * n; ++i) Phi[i] += C[i];
    T.swap(C);
  }
  for (int i = 0; i < s; ++i) {
    mat_mul(Phi, Phi, n, C);
    Phi.swap(C);
  }
}

inline void propagate(const std::vector<double>& K, const std::vector<double>& A0,
                      int n, double dt, std::vector<double>& out) {
  std::vector<double> Phi;
  mat_exp(K, n, dt, Phi);
  mat_vec(Phi, A0, n, out);
}

inline void propagate_infusion(const std::vector<double>& K,
                               const std::vector<double>& A0,
                               const std::vector<double>& rate,
                               int n, double dt, std::vector<double>& out) {
  const int m = n + 1;
  std::vector<double> M(m * m, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) M[i * m + j] = K[i * n + j];
    M[i * m + n] = rate[i];
  }
  std::vector<double> aug(m, 0.0);
  for (int i = 0; i < n; ++i) aug[i] = A0[i];
  aug[n] = 1.0;
  std::vector<double> Phi;
  mat_exp(M, m, dt, Phi);
  std::vector<double> res(m, 0.0);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < m; ++j) res[i] += Phi[i * m + j] * aug[j];
  }
  out.assign(n, 0.0);
  for (int i = 0; i < n; ++i) out[i] = res[i];
}

inline bool ss_bolus_or_oral_pre(const std::vector<double>& K,
                                 const std::vector<double>& D,
                                 int n, double tau, std::vector<double>& A_pre) {
  std::vector<double> Phi, I, M, rhs;
  mat_exp(K, n, tau, Phi);
  mat_identity(n, I);
  M.resize(n * n);
  for (int i = 0; i < n * n; ++i) M[i] = I[i] - Phi[i];
  mat_vec(Phi, D, n, rhs);
  return gauss_solve(M, rhs, n, A_pre);
}

inline bool ss_bolus_or_oral(const std::vector<double>& K,
                             const std::vector<double>& D,
                             int n, double tau, std::vector<double>& A_post) {
  std::vector<double> A_pre;
  if (!ss_bolus_or_oral_pre(K, D, n, tau, A_pre)) return false;
  A_post.resize(n);
  for (int i = 0; i < n; ++i) A_post[i] = A_pre[i] + D[i];
  return true;
}

inline bool ss_infusion(const std::vector<double>& K,
                        const std::vector<double>& rate,
                        int n, double dur, double tau,
                        std::vector<double>& A_pre) {
  if (dur <= 0.0 || tau <= dur) return false;
  std::vector<double> z(n, 0.0), after, end_cycle, Phi, I, M;
  propagate_infusion(K, z, rate, n, dur, after);
  propagate(K, after, n, tau - dur, end_cycle);
  mat_exp(K, n, tau, Phi);
  mat_identity(n, I);
  M.resize(n * n);
  for (int i = 0; i < n * n; ++i) M[i] = I[i] - Phi[i];
  return gauss_solve(M, end_cycle, n, A_pre);
}

}  // namespace ss_mat

inline bool ss_row_active(int model_ss, int row_ss) {
  return model_ss == 1 || row_ss == 1 || row_ss == 3;
}

inline bool ss_row_reset(int row_ss) {
  return row_ss == 2 || row_ss == 3;
}

inline double row_ii(const SubjectEvents& ev, int i) {
  if (ev.ii.size() == 0 || i < 0 || i >= ev.ii.size()) return 0.0;
  return ev.ii[i];
}

inline int row_ss(const SubjectEvents& ev, int i) {
  if (ev.ss.size() == 0 || i < 0 || i >= ev.ss.size()) return 0;
  return ev.ss[i];
}

inline bool is_dosing_row(const SubjectEvents& ev, int i) {
  int evid = ev.evid.size() > 0 ? ev.evid[i] : (ev.amt[i] > 0 ? 1 : 0);
  return evid == 1 || evid == 4;
}

inline int first_dosing_cmt(const SubjectEvents& ev, int dose_cmp) {
  for (int i = 0; i < ev.time.size(); ++i) {
    if (is_dosing_row(ev, i) && ev.cmt.size() > 0 && ev.cmt[i] > 0) {
      return ev.cmt[i];
    }
  }
  return dose_cmp > 0 ? dose_cmp : 1;
}

inline int row_cmt(const SubjectEvents& ev, int i, int default_cmp) {
  if (ev.cmt.size() > 0 && ev.cmt[i] > 0) return ev.cmt[i];
  return default_cmp > 0 ? default_cmp : 1;
}

inline int event_obs_cmp(int default_cmp, const SubjectEvents& ev, int i) {
  if (ev.evid.size() > 0 && ev.evid[i] != 0) return default_cmp;
  return row_cmt(ev, i, default_cmp);
}

inline bool is_reset_row(const SubjectEvents& ev, int i) {
  return ev.evid.size() > 0 && ev.evid[i] == 3;
}

inline bool is_other_event_row(const SubjectEvents& ev, int i) {
  return ev.evid.size() > 0 && ev.evid[i] == 2;
}

inline void reset_amounts(std::vector<double>& amounts) {
  std::fill(amounts.begin(), amounts.end(), 0.0);
}

inline void build_K_1_iv(const PkParams& p, std::vector<double>& K) {
  const double k10 = p.k10 > 0.0 ? p.k10 : effective_ke(p);
  K = { -k10 };
}

inline void build_K_1_oral(const PkParams& p, std::vector<double>& K) {
  const double k10 = p.k10 > 0.0 ? p.k10 : effective_ke(p);
  K = { -p.ka, 0.0, p.ka, -k10 };
}

inline void build_K_2_iv(const PkParams& p, std::vector<double>& K) {
  K = {
    -(p.k10 + p.k12), p.k21,
    p.k12, -p.k21
  };
}

inline void build_K_2_oral(const PkParams& p, std::vector<double>& K) {
  const double k10 = p.k20 > 0.0 ? p.k20 : p.k10;
  K = {
    -p.ka, 0.0, 0.0,
    p.ka, -(k10 + p.k23), p.k32,
    0.0, p.k23, -p.k32
  };
}

inline void build_K_3_iv(const PkParams& p, std::vector<double>& K) {
  K = {
    -(p.k10 + p.k12 + p.k13), p.k21, p.k31,
    p.k12, -p.k21, 0.0,
    p.k13, 0.0, -p.k31
  };
}

inline void build_K_4_oral_chain(const PkParams& p, std::vector<double>& K) {
  const double k10 = p.k20 > 0.0 ? p.k20 : p.k10;
  K = {
    -p.ka, 0.0, 0.0, 0.0,
    p.ka, -(k10 + p.k12), p.k21, 0.0,
    0.0, p.k12, -(p.k21 + p.k23), p.k32,
    0.0, 0.0, p.k23, -p.k32
  };
}

inline void build_K_3_oral(const PkParams& p, std::vector<double>& K) {
  const double k10 = p.k20 > 0.0 ? p.k20 : p.k10;
  K = {
    -p.ka, 0.0, 0.0, 0.0,
    p.ka, -(k10 + p.k23 + p.k24), p.k32, p.k42,
    0.0, p.k23, -p.k32, 0.0,
    0.0, p.k24, 0.0, -p.k42
  };
}

// Returns true if SS initialised amounts (dose already included for bolus/oral).
inline bool apply_ss_at_dose(int model_ss, int row_ss, int advan, RouteType route,
                             const PkParams& p, const SubjectEvents& ev, int i,
                             std::vector<double>& amounts,
                             bool lagged_oral_depot = false) {
  if (ss_row_reset(row_ss)) {
    for (double& a : amounts) a = 0.0;
  }
  const double tau = row_ii(ev, i);
  if (!ss_row_active(model_ss, row_ss) || tau <= 0.0) return false;

  const double dose = ev.amt[i];
  const int dose_cmt = row_cmt(ev, i, 1);
  const double ff = effective_f(dose_cmt, p, ev, i);
  const double rate0 = effective_rate(ev, i);
  const bool infusion = route == ROUTE_IV_INFUSION && rate0 > 0.0;
  const double dur = infusion ? dose / rate0 : 0.0;

  const int ncomp = pk_ncomp(advan);
  std::vector<double> K, D, rate_vec, A;

  if (ncomp == 1 && route == ROUTE_IV_BOLUS) {
    build_K_1_iv(p, K);
    D = { dose * ff };
    if (!ss_mat::ss_bolus_or_oral(K, D, 1, tau, A)) return false;
    amounts[0] = A[0];
    return true;
  }
  if (ncomp == 1 && route == ROUTE_ORAL) {
    build_K_1_oral(p, K);
    D = { lagged_oral_depot ? 0.0 : ff * dose, 0.0 };
    if (!ss_mat::ss_bolus_or_oral(K, D, 2, tau, A)) return false;
    if (amounts.size() < 2) amounts.resize(2, 0.0);
    amounts[0] = A[0];
    amounts[1] = A[1];
    return true;
  }
  if (ncomp == 1 && infusion) {
    build_K_1_iv(p, K);
    rate_vec = { rate0 };
    if (!ss_mat::ss_infusion(K, rate_vec, 1, dur, tau, A)) return false;
    amounts[0] = A[0];
    return true;
  }
  if (ncomp == 2 && route == ROUTE_IV_BOLUS) {
    build_K_2_iv(p, K);
    D = { 0.0, 0.0 };
    if (dose_cmt == 2) D[1] = dose * ff;
    else D[0] = dose * ff;
    if (!ss_mat::ss_bolus_or_oral(K, D, 2, tau, A)) return false;
    if (amounts.size() < 2) amounts.resize(2, 0.0);
    amounts[0] = A[0];
    amounts[1] = A[1];
    return true;
  }
  if (ncomp == 2 && route == ROUTE_ORAL) {
    if (advan == 4 || (p.k23 > 0.0 && p.k32 > 0.0)) {
      build_K_2_oral(p, K);
      D = { ff * dose, 0.0, 0.0 };
      if (!ss_mat::ss_bolus_or_oral(K, D, 3, tau, A)) return false;
      if (amounts.size() < 3) amounts.resize(3, 0.0);
      amounts[0] = A[0];
      amounts[1] = A[1];
      amounts[2] = A[2];
      return true;
    }
    build_K_1_oral(p, K);
    D = { ff * dose, 0.0 };
    if (!ss_mat::ss_bolus_or_oral(K, D, 2, tau, A)) return false;
    if (amounts.size() < 2) amounts.resize(2, 0.0);
    amounts[0] = A[0];
    amounts[1] = A[1];
    return true;
  }
  if (ncomp == 2 && infusion) {
    build_K_2_iv(p, K);
    rate_vec = { rate0, 0.0 };
    if (!ss_mat::ss_infusion(K, rate_vec, 2, dur, tau, A)) return false;
    if (amounts.size() < 2) amounts.resize(2, 0.0);
    amounts[0] = A[0];
    amounts[1] = A[1];
    return true;
  }
  if (ncomp == 3 && route == ROUTE_IV_BOLUS) {
    build_K_3_iv(p, K);
    D = { 0.0, 0.0, 0.0 };
    if (dose_cmt == 2) D[1] = dose * ff;
    else if (dose_cmt == 3) D[2] = dose * ff;
    else D[0] = dose * ff;
    if (!ss_mat::ss_bolus_or_oral(K, D, 3, tau, A)) return false;
    if (amounts.size() < 3) amounts.resize(3, 0.0);
    amounts[0] = A[0]; amounts[1] = A[1]; amounts[2] = A[2];
    return true;
  }
  if (ncomp == 3 && (route == ROUTE_ORAL || advan == 12)) {
    build_K_3_oral(p, K);
    D = { ff * dose, 0.0, 0.0, 0.0 };
    if (!ss_mat::ss_bolus_or_oral(K, D, 4, tau, A)) return false;
    if (amounts.size() < 4) amounts.resize(4, 0.0);
    amounts[0] = A[0]; amounts[1] = A[1]; amounts[2] = A[2]; amounts[3] = A[3];
    return true;
  }
  if (ncomp == 3 && infusion) {
    build_K_3_iv(p, K);
    rate_vec = { rate0, 0.0, 0.0 };
    if (!ss_mat::ss_infusion(K, rate_vec, 3, dur, tau, A)) return false;
    if (amounts.size() < 3) amounts.resize(3, 0.0);
    amounts[0] = A[0]; amounts[1] = A[1]; amounts[2] = A[2];
    return true;
  }
  return false;
}

}  // namespace nm_pk
