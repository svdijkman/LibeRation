#include <Rcpp.h>
#include <cmath>
#include <vector>
#include "nm_pk_core.h"
#include "nm_pk_ss.h"
#include "nm_pk_pkadvan.h"
using namespace Rcpp;
using namespace nm_pk;
using namespace pkadvan;

namespace {

inline double effective_ktr(const PkParams& p) {
  return p.ktr > 0.0 ? p.ktr : p.ka;
}

inline double conc_1comp(int obs_cmp, double a1, double a2, const PkParams& p,
                         const SubjectEvents& ev, int i) {
  return conc_oral(obs_cmp, a1, a2, 0.0, 0.0, p, ev, i);
}

inline double conc_2comp(int obs_cmp, double a1, double a2, double a3, const PkParams& p,
                         const SubjectEvents& ev, int i) {
  return conc_oral(obs_cmp, a1, a2, a3, 0.0, p, ev, i);
}

inline void hybrid_2_oral(double k20, double k23, double k32, double k30,
                        double& E2, double& E3, double& lambda1, double& lambda2) {
  E2 = k20 + k23;
  E3 = k32 + k30;
  const double sum = E2 + E3;
  const double root = std::sqrt(std::max(sum * sum - 4.0 * (E2 * E3 - k23 * k32), 0.0));
  lambda1 = 0.5 * (sum + root);
  lambda2 = 0.5 * (sum - root);
}

void apply_dose_1comp(int i, const SubjectEvents& ev, const PkParams& p, int model_ss,
                      int n_state, std::vector<double>& amounts) {
  if (!is_dosing_row(ev, i)) return;
  if (ss_row_reset(row_ss(ev, i))) {
    for (int k = 0; k < n_state; ++k) amounts[k] = 0.0;
  }
  const int cmt = row_cmt(ev, i, 1);
  const double ff = effective_f(cmt, p, ev, i);
  amounts[0] += ev.amt[i] * ff;
}

void apply_dose_2comp(int i, const SubjectEvents& ev, const PkParams& p, int model_ss,
                      int n_state, std::vector<double>& amounts) {
  if (!is_dosing_row(ev, i)) return;
  if (ss_row_reset(row_ss(ev, i))) {
    for (int k = 0; k < n_state; ++k) amounts[k] = 0.0;
  }
  const int cmt = row_cmt(ev, i, 1);
  const double ff = effective_f(cmt, p, ev, i);
  amounts[0] += ev.amt[i] * ff;
}

// amounts[0]=A1, amounts[1]=A2, amounts[2]=A3 (1st transit)
NumericVector solve_1comp_1transit(const SubjectEvents& ev, const PkParams& p,
                                   int obs_cmp, int model_ss) {
  const int n = ev.time.size();
  NumericVector ipred(n);
  const double ktr = effective_ktr(p);
  const double k20 = effective_ke(p);
  double a1 = 0.0, a2 = 0.0, a3 = 0.0;
  for (int i = 0; i < n; ++i) {
    if (i > 0) {
      const double dt = ev.time[i] - ev.time[i - 1];
      const double previousA1 = a1;
      const double previousA2 = a2;
      const double previousA3 = a3;
      a3 = previousA3 * std::exp(-dt * ktr) + ktr * previousA1 * dt * std::exp(-dt * ktr);
      a2 = previousA2 * std::exp(-dt * k20);
      a2 = a2 + ktr * (previousA3 * (std::exp(-dt * k20) / (ktr - k20) + std::exp(-dt * ktr) / (k20 - ktr)) +
                       previousA1 * ktr * (dt * std::exp(-dt * ktr) / (k20 - ktr) +
                                           std::exp(-dt * k20) / std::pow(k20 - ktr, 2.0) -
                                           std::exp(-dt * ktr) / std::pow(k20 - ktr, 2.0)));
      a1 = previousA1 * std::exp(-dt * ktr);
    }
    std::vector<double> amounts = {a1, a2, a3};
    apply_dose_1comp(i, ev, p, model_ss, 3, amounts);
    a1 = amounts[0]; a2 = amounts[1]; a3 = amounts[2];
    ipred[i] = conc_1comp(event_obs_cmp(obs_cmp, ev, i), a1, a2, p, ev, i);
    pk_write_amounts(i, amounts);
  }
  return ipred;
}

// amounts[0]=A1, [1]=A2, [2]=A3, [3]=A4
NumericVector solve_1comp_2transit(const SubjectEvents& ev, const PkParams& p,
                                   int obs_cmp, int model_ss) {
  const int n = ev.time.size();
  NumericVector ipred(n);
  const double ktr = effective_ktr(p);
  const double k20 = effective_ke(p);
  double a1 = 0.0, a2 = 0.0, a3 = 0.0, a4 = 0.0;
  for (int i = 0; i < n; ++i) {
    if (i > 0) {
      const double dt = ev.time[i] - ev.time[i - 1];
      const double previousA1 = a1;
      const double previousA2 = a2;
      const double previousA3 = a3;
      const double previousA4 = a4;
      a3 = previousA3 * std::exp(-dt * ktr) + ktr * previousA1 * dt * std::exp(-dt * ktr);
      a4 = previousA4 * std::exp(-dt * ktr) + ktr * previousA3 * dt * std::exp(-dt * ktr) +
           0.5 * std::pow(ktr, 2.0) * previousA1 * std::pow(dt, 2.0) * std::exp(-dt * ktr);
      a2 = previousA2 * std::exp(-dt * k20);
      a2 = a2 + ktr * (previousA4 * (std::exp(-dt * k20) / (ktr - k20) + std::exp(-dt * ktr) / (k20 - ktr)) +
                       previousA3 * ktr * (dt * std::exp(-dt * ktr) / (k20 - ktr) +
                                           std::exp(-dt * k20) / std::pow(k20 - ktr, 2.0) -
                                           std::exp(-dt * ktr) / std::pow(k20 - ktr, 2.0)) +
                       previousA1 * std::pow(ktr, 2.0) * (std::pow(dt, 2.0) * std::exp(-dt * ktr) / (2.0 * (k20 - ktr)) -
                                                          dt * std::exp(-dt * ktr) / std::pow(k20 - ktr, 2.0) -
                                                          std::exp(-dt * k20) / std::pow(k20 - ktr, 3.0) +
                                                          std::exp(-dt * ktr) / std::pow(k20 - ktr, 3.0)));
      a1 = previousA1 * std::exp(-dt * ktr);
    }
    std::vector<double> amounts = {a1, a2, a3, a4};
    apply_dose_1comp(i, ev, p, model_ss, 4, amounts);
    a1 = amounts[0]; a2 = amounts[1]; a3 = amounts[2]; a4 = amounts[3];
    ipred[i] = conc_1comp(event_obs_cmp(obs_cmp, ev, i), a1, a2, p, ev, i);
    pk_write_amounts(i, amounts);
  }
  return ipred;
}

// amounts[0]=A1, [1]=A2, [2]=A3, [3]=A4, [4]=A5
NumericVector solve_1comp_3transit(const SubjectEvents& ev, const PkParams& p,
                                   int obs_cmp, int model_ss) {
  const int n = ev.time.size();
  NumericVector ipred(n);
  const double ktr = effective_ktr(p);
  const double k20 = effective_ke(p);
  double a1 = 0.0, a2 = 0.0, a3 = 0.0, a4 = 0.0, a5 = 0.0;
  for (int i = 0; i < n; ++i) {
    if (i > 0) {
      const double dt = ev.time[i] - ev.time[i - 1];
      const double previousA1 = a1;
      const double previousA2 = a2;
      const double previousA3 = a3;
      const double previousA4 = a4;
      const double previousA5 = a5;
      a3 = previousA3 * std::exp(-dt * ktr) + ktr * previousA1 * dt * std::exp(-dt * ktr);
      a4 = previousA4 * std::exp(-dt * ktr) + ktr * previousA3 * dt * std::exp(-dt * ktr) +
           0.5 * std::pow(ktr, 2.0) * previousA1 * std::pow(dt, 2.0) * std::exp(-dt * ktr);
      a5 = previousA5 * std::exp(-dt * ktr) + ktr * previousA4 * dt * std::exp(-dt * ktr) +
           0.5 * std::pow(ktr, 2.0) * previousA3 * std::pow(dt, 2.0) * std::exp(-dt * ktr) +
           (1.0 / 6.0) * std::pow(ktr, 3.0) * previousA1 * std::pow(dt, 3.0) * std::exp(-dt * ktr);
      a2 = previousA2 * std::exp(-dt * k20);
      a2 = a2 + ktr * (previousA5 * (std::exp(-dt * k20) / (ktr - k20) + std::exp(-dt * ktr) / (k20 - ktr)) +
                       previousA4 * ktr * (dt * std::exp(-dt * ktr) / (k20 - ktr) +
                                           std::exp(-dt * k20) / std::pow(k20 - ktr, 2.0) -
                                           std::exp(-dt * ktr) / std::pow(k20 - ktr, 2.0)) +
                       previousA3 * std::pow(ktr, 2.0) * (std::pow(dt, 2.0) * std::exp(-dt * ktr) / (2.0 * (k20 - ktr)) -
                                                          dt * std::exp(-dt * ktr) / std::pow(k20 - ktr, 2.0) -
                                                          std::exp(-dt * k20) / std::pow(k20 - ktr, 3.0) +
                                                          std::exp(-dt * ktr) / std::pow(k20 - ktr, 3.0)) +
                       previousA1 * std::pow(ktr, 3.0) * (std::pow(dt, 3.0) * std::exp(-dt * ktr) / (6.0 * (k20 - ktr)) -
                                                          std::pow(dt, 2.0) * std::exp(-dt * ktr) / (2.0 * std::pow(k20 - ktr, 2.0)) +
                                                          dt * std::exp(-dt * ktr) / std::pow(k20 - ktr, 3.0) +
                                                          std::exp(-dt * k20) / std::pow(k20 - ktr, 4.0) -
                                                          std::exp(-dt * ktr) / std::pow(k20 - ktr, 4.0)));
      a1 = previousA1 * std::exp(-dt * ktr);
    }
    std::vector<double> amounts = {a1, a2, a3, a4, a5};
    apply_dose_1comp(i, ev, p, model_ss, 5, amounts);
    a1 = amounts[0]; a2 = amounts[1]; a3 = amounts[2]; a4 = amounts[3]; a5 = amounts[4];
    ipred[i] = conc_1comp(event_obs_cmp(obs_cmp, ev, i), a1, a2, p, ev, i);
    pk_write_amounts(i, amounts);
  }
  return ipred;
}

// amounts[0]=A1, [1]=A2, [2]=A3, [3]=A4, [4]=A5, [5]=A6
NumericVector solve_1comp_4transit(const SubjectEvents& ev, const PkParams& p,
                                   int obs_cmp, int model_ss) {
  const int n = ev.time.size();
  NumericVector ipred(n);
  const double ktr = effective_ktr(p);
  const double k20 = effective_ke(p);
  double a1 = 0.0, a2 = 0.0, a3 = 0.0, a4 = 0.0, a5 = 0.0, a6 = 0.0;
  for (int i = 0; i < n; ++i) {
    if (i > 0) {
      const double dt = ev.time[i] - ev.time[i - 1];
      const double previousA1 = a1;
      const double previousA2 = a2;
      const double previousA3 = a3;
      const double previousA4 = a4;
      const double previousA5 = a5;
      const double previousA6 = a6;
      a3 = previousA3 * std::exp(-dt * ktr) + ktr * previousA1 * dt * std::exp(-dt * ktr);
      a4 = previousA4 * std::exp(-dt * ktr) + ktr * previousA3 * dt * std::exp(-dt * ktr) +
           0.5 * std::pow(ktr, 2.0) * previousA1 * std::pow(dt, 2.0) * std::exp(-dt * ktr);
      a5 = previousA5 * std::exp(-dt * ktr) + ktr * previousA4 * dt * std::exp(-dt * ktr) +
           0.5 * std::pow(ktr, 2.0) * previousA3 * std::pow(dt, 2.0) * std::exp(-dt * ktr) +
           (1.0 / 6.0) * std::pow(ktr, 3.0) * previousA1 * std::pow(dt, 3.0) * std::exp(-dt * ktr);
      a6 = previousA6 * std::exp(-dt * ktr) + ktr * previousA5 * dt * std::exp(-dt * ktr) +
           0.5 * std::pow(ktr, 2.0) * previousA4 * std::pow(dt, 2.0) * std::exp(-dt * ktr) +
           (1.0 / 6.0) * std::pow(ktr, 3.0) * previousA3 * std::pow(dt, 3.0) * std::exp(-dt * ktr) +
           (1.0 / 24.0) * std::pow(ktr, 4.0) * previousA1 * std::pow(dt, 4.0) * std::exp(-dt * ktr);
      a2 = previousA2 * std::exp(-dt * k20);
      a2 = a2 + ktr * (previousA6 * (std::exp(-dt * k20) / (ktr - k20) + std::exp(-dt * ktr) / (k20 - ktr)) +
                       previousA5 * ktr * (dt * std::exp(-dt * ktr) / (k20 - ktr) +
                                           std::exp(-dt * k20) / std::pow(k20 - ktr, 2.0) -
                                           std::exp(-dt * ktr) / std::pow(k20 - ktr, 2.0)) +
                       previousA4 * std::pow(ktr, 2.0) * (std::pow(dt, 2.0) * std::exp(-dt * ktr) / (2.0 * (k20 - ktr)) -
                                                          dt * std::exp(-dt * ktr) / std::pow(k20 - ktr, 2.0) -
                                                          std::exp(-dt * k20) / std::pow(k20 - ktr, 3.0) +
                                                          std::exp(-dt * ktr) / std::pow(k20 - ktr, 3.0)) +
                       previousA3 * std::pow(ktr, 3.0) * (std::pow(dt, 3.0) * std::exp(-dt * ktr) / (6.0 * (k20 - ktr)) -
                                                          std::pow(dt, 2.0) * std::exp(-dt * ktr) / (2.0 * std::pow(k20 - ktr, 2.0)) +
                                                          dt * std::exp(-dt * ktr) / std::pow(k20 - ktr, 3.0) +
                                                          std::exp(-dt * k20) / std::pow(k20 - ktr, 4.0) -
                                                          std::exp(-dt * ktr) / std::pow(k20 - ktr, 4.0)) +
                       previousA1 * std::pow(ktr, 4.0) * (std::pow(dt, 4.0) * std::exp(-dt * ktr) / (24.0 * (k20 - ktr)) -
                                                          std::pow(dt, 3.0) * std::exp(-dt * ktr) / (6.0 * std::pow(k20 - ktr, 2.0)) +
                                                          std::pow(dt, 2.0) * std::exp(-dt * ktr) / (2.0 * std::pow(k20 - ktr, 3.0)) -
                                                          dt * std::exp(-dt * ktr) / std::pow(k20 - ktr, 4.0) -
                                                          std::exp(-dt * k20) / std::pow(k20 - ktr, 5.0) +
                                                          std::exp(-dt * ktr) / std::pow(k20 - ktr, 5.0)));
      a1 = previousA1 * std::exp(-dt * ktr);
    }
    std::vector<double> amounts = {a1, a2, a3, a4, a5, a6};
    apply_dose_1comp(i, ev, p, model_ss, 6, amounts);
    a1 = amounts[0]; a2 = amounts[1]; a3 = amounts[2]; a4 = amounts[3]; a5 = amounts[4]; a6 = amounts[5];
    ipred[i] = conc_1comp(event_obs_cmp(obs_cmp, ev, i), a1, a2, p, ev, i);
    pk_write_amounts(i, amounts);
  }
  return ipred;
}

// amounts[0]=A1, [1]=A2, [2]=A3, [3]=A4 (1st transit)
NumericVector solve_2comp_1transit(const SubjectEvents& ev, const PkParams& p,
                                   int obs_cmp, int model_ss) {
  const int n = ev.time.size();
  NumericVector ipred(n);
  const double ktr = effective_ktr(p);
  const double k20 = effective_ke(p);
  const double k23 = p.k23;
  const double k32 = p.k32;
  const double k30 = p.k30;
  double a1 = 0.0, a2 = 0.0, a3 = 0.0, a4 = 0.0;
  for (int i = 0; i < n; ++i) {
    if (i > 0) {
      const double dt = ev.time[i] - ev.time[i - 1];
      double E2, E3, lambda1, lambda2;
      hybrid_2_oral(k20, k23, k32, k30, E2, E3, lambda1, lambda2);
      const double previousA1 = a1;
      const double previousA2 = a2;
      const double previousA3 = a3;
      const double previousA4 = a4;
      a4 = previousA4 * std::exp(-dt * ktr) + ktr * previousA1 * dt * std::exp(-dt * ktr);
      a2 = (std::exp(-dt * lambda1) * ((previousA2 * E3 + previousA3 * k32) - previousA2 * lambda1) -
            std::exp(-dt * lambda2) * ((previousA2 * E3 + previousA3 * k32) - previousA2 * lambda2)) / (lambda2 - lambda1);
      a2 = a2 + ktr * E3 * (previousA4 * (std::exp(-dt * ktr) / ((lambda1 - ktr) * (lambda2 - ktr)) +
                                           std::exp(-dt * lambda1) / ((ktr - lambda1) * (lambda2 - lambda1)) +
                                           std::exp(-dt * lambda2) / ((ktr - lambda2) * (lambda1 - lambda2))) +
                              previousA1 * ktr * (std::exp(-dt * ktr) * (-lambda1 - lambda2 + 2.0 * ktr) / (std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) -
                                                  std::exp(-dt * lambda1) / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 2.0)) +
                                                  std::exp(-dt * lambda2) / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 2.0)) -
                                                  std::exp(-dt * ktr) * dt / ((lambda1 - ktr) * (ktr - lambda2))));
      a2 = a2 + ktr * (previousA4 * (std::exp(-dt * ktr) * ktr / ((lambda1 - ktr) * (ktr - lambda2)) +
                                     std::exp(-dt * lambda2) * lambda2 / ((lambda1 - lambda2) * (lambda2 - ktr)) -
                                     std::exp(-dt * lambda1) * lambda1 / ((lambda1 - lambda2) * (lambda1 - ktr))) +
                       previousA1 * ktr * (std::exp(-dt * ktr) * (lambda1 * lambda2 - std::pow(ktr, 2.0)) / (std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) +
                                           std::exp(-dt * lambda1) * lambda1 / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 2.0)) -
                                           std::exp(-dt * lambda2) * lambda2 / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 2.0)) +
                                           std::exp(-dt * ktr) * dt * ktr / ((lambda1 - ktr) * (ktr - lambda2))));
      a3 = (std::exp(-dt * lambda1) * ((previousA3 * E2 + k23 * previousA2) - previousA3 * lambda1) -
            std::exp(-dt * lambda2) * ((previousA3 * E2 + k23 * previousA2) - previousA3 * lambda2)) / (lambda2 - lambda1);
      a3 = a3 + ktr * k23 * (previousA4 * (std::exp(-dt * ktr) / ((lambda1 - ktr) * (lambda2 - ktr)) +
                                            std::exp(-dt * lambda1) / ((ktr - lambda1) * (lambda2 - lambda1)) +
                                            std::exp(-dt * lambda2) / ((ktr - lambda2) * (lambda1 - lambda2))) +
                             previousA1 * ktr * (std::exp(-dt * ktr) * (-lambda1 - lambda2 + 2.0 * ktr) / (std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) -
                                                 std::exp(-dt * lambda1) / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 2.0)) +
                                                 std::exp(-dt * lambda2) / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 2.0)) -
                                                 std::exp(-dt * ktr) * dt / ((lambda1 - ktr) * (ktr - lambda2))));
      a1 = previousA1 * std::exp(-dt * ktr);
    }
    std::vector<double> amounts = {a1, a2, a3, a4};
    apply_dose_2comp(i, ev, p, model_ss, 4, amounts);
    a1 = amounts[0]; a2 = amounts[1]; a3 = amounts[2]; a4 = amounts[3];
    ipred[i] = conc_2comp(event_obs_cmp(obs_cmp, ev, i), a1, a2, a3, p, ev, i);
    pk_write_amounts(i, amounts);
  }
  return ipred;
}

// amounts[0]=A1, [1]=A2, [2]=A3, [3]=A4, [4]=A5
NumericVector solve_2comp_2transit(const SubjectEvents& ev, const PkParams& p,
                                   int obs_cmp, int model_ss) {
  const int n = ev.time.size();
  NumericVector ipred(n);
  const double ktr = effective_ktr(p);
  const double k20 = effective_ke(p);
  const double k23 = p.k23;
  const double k32 = p.k32;
  const double k30 = p.k30;
  double a1 = 0.0, a2 = 0.0, a3 = 0.0, a4 = 0.0, a5 = 0.0;
  for (int i = 0; i < n; ++i) {
    if (i > 0) {
      const double dt = ev.time[i] - ev.time[i - 1];
      double E2, E3, lambda1, lambda2;
      hybrid_2_oral(k20, k23, k32, k30, E2, E3, lambda1, lambda2);
      const double previousA1 = a1;
      const double previousA2 = a2;
      const double previousA3 = a3;
      const double previousA4 = a4;
      const double previousA5 = a5;
      a4 = previousA4 * std::exp(-dt * ktr) + ktr * previousA1 * dt * std::exp(-dt * ktr);
      a5 = previousA5 * std::exp(-dt * ktr) + ktr * previousA4 * dt * std::exp(-dt * ktr) +
           0.5 * std::pow(ktr, 2.0) * previousA1 * std::pow(dt, 2.0) * std::exp(-dt * ktr);
      a2 = (std::exp(-dt * lambda1) * ((previousA2 * E3 + previousA3 * k32) - previousA2 * lambda1) -
            std::exp(-dt * lambda2) * ((previousA2 * E3 + previousA3 * k32) - previousA2 * lambda2)) / (lambda2 - lambda1);
      a2 = a2 + ktr * E3 * (previousA5 * (std::exp(-dt * ktr) / ((lambda1 - ktr) * (lambda2 - ktr)) +
                                           std::exp(-dt * lambda1) / ((ktr - lambda1) * (lambda2 - lambda1)) +
                                           std::exp(-dt * lambda2) / ((ktr - lambda2) * (lambda1 - lambda2))) +
                              previousA4 * ktr * (std::exp(-dt * ktr) * (-lambda1 - lambda2 + 2.0 * ktr) / (std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) -
                                                  std::exp(-dt * lambda1) / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 2.0)) +
                                                  std::exp(-dt * lambda2) / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 2.0)) -
                                                  std::exp(-dt * ktr) * dt / ((lambda1 - ktr) * (ktr - lambda2))) +
                              previousA1 * std::pow(ktr, 2.0) * ((std::exp(-dt * ktr) * (-std::pow(lambda1, 2.0) - lambda1 * lambda2 + 3.0 * lambda1 * ktr - std::pow(lambda2, 2.0) + 3.0 * lambda2 * ktr - 3.0 * std::pow(ktr, 2.0))) / (std::pow(lambda1 - ktr, 3.0) * std::pow(ktr - lambda2, 3.0)) -
                                                                 std::exp(-dt * ktr) * std::pow(dt, 2.0) / (2.0 * (lambda1 - ktr) * (ktr - lambda2)) +
                                                                 std::exp(-dt * ktr) * dt * (-lambda1 - lambda2 + 2.0 * ktr) / (std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) +
                                                                 std::exp(-dt * lambda1) / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 3.0)) -
                                                                 std::exp(-dt * lambda2) / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 3.0))));
      a2 = a2 + ktr * (previousA5 * (std::exp(-dt * ktr) * ktr / ((lambda1 - ktr) * (ktr - lambda2)) +
                                     std::exp(-dt * lambda2) * lambda2 / ((lambda1 - lambda2) * (lambda2 - ktr)) -
                                     std::exp(-dt * lambda1) * lambda1 / ((lambda1 - lambda2) * (lambda1 - ktr))) +
                       previousA4 * ktr * (std::exp(-dt * ktr) * (lambda1 * lambda2 - std::pow(ktr, 2.0)) / (std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) +
                                           std::exp(-dt * lambda1) * lambda1 / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 2.0)) -
                                           std::exp(-dt * lambda2) * lambda2 / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 2.0)) +
                                           std::exp(-dt * ktr) * dt * ktr / ((lambda1 - ktr) * (ktr - lambda2))) +
                       previousA1 * std::pow(ktr, 2.0) * (std::exp(-dt * ktr) * (std::pow(lambda1, 2.0) * lambda2 + lambda1 * std::pow(lambda2, 2.0) - 3.0 * lambda1 * lambda2 * ktr + std::pow(ktr, 3.0)) / (std::pow(lambda1 - ktr, 3.0) * std::pow(ktr - lambda2, 3.0)) +
                                                           std::exp(-dt * ktr) * dt * (lambda1 * lambda2 - std::pow(ktr, 2.0)) / (std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) +
                                                           std::exp(-dt * ktr) * ktr * std::pow(dt, 2.0) / (2.0 * (lambda1 - ktr) * (ktr - lambda2)) -
                                                           std::exp(-dt * lambda1) * lambda1 / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 3.0)) +
                                                           std::exp(-dt * lambda2) * lambda2 / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 3.0))));
      a3 = (std::exp(-dt * lambda1) * ((previousA3 * E2 + k23 * previousA2) - previousA3 * lambda1) -
            std::exp(-dt * lambda2) * ((previousA3 * E2 + k23 * previousA2) - previousA3 * lambda2)) / (lambda2 - lambda1);
      a3 = a3 + ktr * k23 * (previousA5 * (std::exp(-dt * ktr) / ((lambda1 - ktr) * (lambda2 - ktr)) +
                                            std::exp(-dt * lambda1) / ((ktr - lambda1) * (lambda2 - lambda1)) +
                                            std::exp(-dt * lambda2) / ((ktr - lambda2) * (lambda1 - lambda2))) +
                             previousA4 * ktr * (std::exp(-dt * ktr) * (-lambda1 - lambda2 + 2.0 * ktr) / (std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) -
                                                 std::exp(-dt * lambda1) / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 2.0)) +
                                                 std::exp(-dt * lambda2) / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 2.0)) -
                                                 std::exp(-dt * ktr) * dt / ((lambda1 - ktr) * (ktr - lambda2))) +
                             previousA1 * std::pow(ktr, 2.0) * ((std::exp(-dt * ktr) * (-std::pow(lambda1, 2.0) - lambda1 * lambda2 + 3.0 * lambda1 * ktr - std::pow(lambda2, 2.0) + 3.0 * lambda2 * ktr - 3.0 * std::pow(ktr, 2.0))) / (std::pow(lambda1 - ktr, 3.0) * std::pow(ktr - lambda2, 3.0)) -
                                                                 std::exp(-dt * ktr) * std::pow(dt, 2.0) / (2.0 * (lambda1 - ktr) * (ktr - lambda2)) +
                                                                 std::exp(-dt * ktr) * dt * (-lambda1 - lambda2 + 2.0 * ktr) / (std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) +
                                                                 std::exp(-dt * lambda1) / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 3.0)) -
                                                                 std::exp(-dt * lambda2) / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 3.0))));
      a1 = previousA1 * std::exp(-dt * ktr);
    }
    std::vector<double> amounts = {a1, a2, a3, a4, a5};
    apply_dose_2comp(i, ev, p, model_ss, 5, amounts);
    a1 = amounts[0]; a2 = amounts[1]; a3 = amounts[2]; a4 = amounts[3]; a5 = amounts[4];
    ipred[i] = conc_2comp(event_obs_cmp(obs_cmp, ev, i), a1, a2, a3, p, ev, i);
    pk_write_amounts(i, amounts);
  }
  return ipred;
}

// amounts[0]=A1, [1]=A2, [2]=A3, [3]=A4, [4]=A5, [5]=A6
NumericVector solve_2comp_3transit(const SubjectEvents& ev, const PkParams& p,
                                   int obs_cmp, int model_ss) {
  const int n = ev.time.size();
  NumericVector ipred(n);
  const double ktr = effective_ktr(p);
  const double k20 = effective_ke(p);
  const double k23 = p.k23;
  const double k32 = p.k32;
  const double k30 = p.k30;
  double a1 = 0.0, a2 = 0.0, a3 = 0.0, a4 = 0.0, a5 = 0.0, a6 = 0.0;
  for (int i = 0; i < n; ++i) {
    if (i > 0) {
      const double dt = ev.time[i] - ev.time[i - 1];
      double E2, E3, lambda1, lambda2;
      hybrid_2_oral(k20, k23, k32, k30, E2, E3, lambda1, lambda2);
      const double previousA1 = a1;
      const double previousA2 = a2;
      const double previousA3 = a3;
      const double previousA4 = a4;
      const double previousA5 = a5;
      const double previousA6 = a6;
      a4 = previousA4 * std::exp(-dt * ktr) + ktr * previousA1 * dt * std::exp(-dt * ktr);
      a5 = previousA5 * std::exp(-dt * ktr) + ktr * previousA4 * dt * std::exp(-dt * ktr) +
           0.5 * std::pow(ktr, 2.0) * previousA1 * std::pow(dt, 2.0) * std::exp(-dt * ktr);
      a6 = previousA6 * std::exp(-dt * ktr) + ktr * previousA5 * dt * std::exp(-dt * ktr) +
           0.5 * std::pow(ktr, 2.0) * previousA4 * std::pow(dt, 2.0) * std::exp(-dt * ktr) +
           (1.0 / 6.0) * std::pow(ktr, 3.0) * previousA1 * std::pow(dt, 3.0) * std::exp(-dt * ktr);
      a2 = (std::exp(-dt * lambda1) * ((previousA2 * E3 + previousA3 * k32) - previousA2 * lambda1) -
            std::exp(-dt * lambda2) * ((previousA2 * E3 + previousA3 * k32) - previousA2 * lambda2)) / (lambda2 - lambda1);
      a2 = a2 + ktr * k32 * (previousA6 * (std::exp(-dt * ktr) / ((lambda1 - ktr) * (lambda2 - ktr)) +
                                            std::exp(-dt * lambda1) / ((ktr - lambda1) * (lambda2 - lambda1)) +
                                            std::exp(-dt * lambda2) / ((ktr - lambda2) * (lambda1 - lambda2))) +
                              previousA5 * ktr * (std::exp(-dt * ktr) * (-lambda1 - lambda2 + 2.0 * ktr) / (std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) -
                                                  std::exp(-dt * lambda1) / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 2.0)) +
                                                  std::exp(-dt * lambda2) / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 2.0)) -
                                                  std::exp(-dt * ktr) * dt / ((lambda1 - ktr) * (ktr - lambda2))) +
                              previousA4 * std::pow(ktr, 2.0) * ((std::exp(-dt * ktr) * (-std::pow(lambda1, 2.0) - lambda1 * lambda2 + 3.0 * lambda1 * ktr - std::pow(lambda2, 2.0) + 3.0 * lambda2 * ktr - 3.0 * std::pow(ktr, 2.0))) / (std::pow(lambda1 - ktr, 3.0) * std::pow(ktr - lambda2, 3.0)) -
                                                                 std::exp(-dt * ktr) * std::pow(dt, 2.0) / (2.0 * (lambda1 - ktr) * (ktr - lambda2)) +
                                                                 std::exp(-dt * ktr) * dt * (-lambda1 - lambda2 + 2.0 * ktr) / (std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) +
                                                                 std::exp(-dt * lambda1) / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 3.0)) -
                                                                 std::exp(-dt * lambda2) / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 3.0))) +
                              previousA1 * std::pow(ktr, 3.0) * ((std::exp(-dt * ktr) * dt * (-std::pow(lambda1, 2.0) - lambda1 * lambda2 + 3.0 * lambda1 * ktr - std::pow(lambda2, 2.0) + 3.0 * lambda2 * ktr - 3.0 * std::pow(ktr, 2.0))) / (std::pow(lambda1 - ktr, 3.0) * std::pow(ktr - lambda2, 3.0)) +
                                                                 std::exp(-dt * ktr) * (-std::pow(lambda1, 3.0) - std::pow(lambda1, 2.0) * lambda2 + 4.0 * std::pow(lambda1, 2.0) * ktr - lambda1 * std::pow(lambda2, 2.0) + 4.0 * lambda1 * lambda2 * ktr - 6.0 * lambda1 * std::pow(ktr, 2.0) - std::pow(lambda2, 3.0) + 4.0 * std::pow(lambda2, 2.0) * ktr - 6.0 * lambda2 * std::pow(ktr, 2.0) + 4.0 * std::pow(ktr, 3.0)) / (std::pow(lambda1 - ktr, 4.0) * std::pow(ktr - lambda2, 4.0)) -
                                                                 std::exp(-dt * ktr) * std::pow(dt, 3.0) / (6.0 * (lambda1 - ktr) * (ktr - lambda2)) +
                                                                 std::exp(-dt * ktr) * std::pow(dt, 2.0) * (-lambda1 - lambda2 + 2.0 * ktr) / (2.0 * std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) -
                                                                 std::exp(-dt * lambda1) / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 4.0)) +
                                                                 std::exp(-dt * lambda2) / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 4.0))));
      a2 = a2 + ktr * (previousA6 * (std::exp(-dt * ktr) * ktr / ((lambda1 - ktr) * (ktr - lambda2)) +
                                     std::exp(-dt * lambda2) * lambda2 / ((lambda1 - lambda2) * (lambda2 - ktr)) -
                                     std::exp(-dt * lambda1) * lambda1 / ((lambda1 - lambda2) * (lambda1 - ktr))) +
                       previousA5 * ktr * (std::exp(-dt * ktr) * (lambda1 * lambda2 - std::pow(ktr, 2.0)) / (std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) +
                                           std::exp(-dt * lambda1) * lambda1 / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 2.0)) -
                                           std::exp(-dt * lambda2) * lambda2 / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 2.0)) +
                                           std::exp(-dt * ktr) * dt * ktr / ((lambda1 - ktr) * (ktr - lambda2))) +
                       previousA4 * std::pow(ktr, 2.0) * (std::exp(-dt * ktr) * (std::pow(lambda1, 2.0) * lambda2 + lambda1 * std::pow(lambda2, 2.0) - 3.0 * lambda1 * lambda2 * ktr + std::pow(ktr, 3.0)) / (std::pow(lambda1 - ktr, 3.0) * std::pow(ktr - lambda2, 3.0)) +
                                                           std::exp(-dt * ktr) * dt * (lambda1 * lambda2 - std::pow(ktr, 2.0)) / (std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) +
                                                           std::exp(-dt * ktr) * ktr * std::pow(dt, 2.0) / (2.0 * (lambda1 - ktr) * (ktr - lambda2)) -
                                                           std::exp(-dt * lambda1) * lambda1 / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 3.0)) +
                                                           std::exp(-dt * lambda2) * lambda2 / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 3.0))) +
                       previousA1 * std::pow(ktr, 3.0) * (std::exp(-dt * ktr) * dt * (std::pow(lambda1, 2.0) * lambda2 + lambda1 * std::pow(lambda2, 2.0) - 3.0 * lambda1 * lambda2 * ktr + std::pow(ktr, 3.0)) / (std::pow(lambda1 - ktr, 3.0) * std::pow(ktr - lambda2, 3.0)) +
                                                           std::exp(-dt * ktr) * (std::pow(lambda1, 3.0) * lambda2 + std::pow(lambda1, 2.0) * std::pow(lambda2, 2.0) - 4.0 * std::pow(lambda1, 2.0) * lambda2 * ktr + lambda1 * std::pow(lambda2, 3.0) - 4.0 * lambda1 * std::pow(lambda2, 2.0) * ktr + 6.0 * lambda1 * lambda2 * std::pow(ktr, 2.0) - std::pow(ktr, 4.0)) / (std::pow(lambda1 - ktr, 4.0) * std::pow(ktr - lambda2, 4.0)) +
                                                           std::exp(-dt * ktr) * std::pow(dt, 2.0) * (lambda1 * lambda2 - std::pow(ktr, 2.0)) / (2.0 * std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) +
                                                           std::exp(-dt * ktr) * ktr * std::pow(dt, 3.0) / (6.0 * (lambda1 - ktr) * (ktr - lambda2)) +
                                                           std::exp(-dt * lambda1) * lambda1 / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 4.0)) -
                                                           std::exp(-dt * lambda2) * lambda2 / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 4.0))));
      a3 = (std::exp(-dt * lambda1) * ((previousA3 * E2 + k23 * previousA2) - previousA3 * lambda1) -
            std::exp(-dt * lambda2) * ((previousA3 * E2 + k23 * previousA2) - previousA3 * lambda2)) / (lambda2 - lambda1);
      a3 = a3 + ktr * k23 * (previousA6 * (std::exp(-dt * ktr) / ((lambda1 - ktr) * (lambda2 - ktr)) +
                                            std::exp(-dt * lambda1) / ((ktr - lambda1) * (lambda2 - lambda1)) +
                                            std::exp(-dt * lambda2) / ((ktr - lambda2) * (lambda1 - lambda2))) +
                             previousA5 * ktr * (std::exp(-dt * ktr) * (-lambda1 - lambda2 + 2.0 * ktr) / (std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) -
                                                 std::exp(-dt * lambda1) / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 2.0)) +
                                                 std::exp(-dt * lambda2) / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 2.0)) -
                                                 std::exp(-dt * ktr) * dt / ((lambda1 - ktr) * (ktr - lambda2))) +
                             previousA4 * std::pow(ktr, 2.0) * ((std::exp(-dt * ktr) * (-std::pow(lambda1, 2.0) - lambda1 * lambda2 + 3.0 * lambda1 * ktr - std::pow(lambda2, 2.0) + 3.0 * lambda2 * ktr - 3.0 * std::pow(ktr, 2.0))) / (std::pow(lambda1 - ktr, 3.0) * std::pow(ktr - lambda2, 3.0)) -
                                                                 std::exp(-dt * ktr) * std::pow(dt, 2.0) / (2.0 * (lambda1 - ktr) * (ktr - lambda2)) +
                                                                 std::exp(-dt * ktr) * dt * (-lambda1 - lambda2 + 2.0 * ktr) / (std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) +
                                                                 std::exp(-dt * lambda1) / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 3.0)) -
                                                                 std::exp(-dt * lambda2) / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 3.0))) +
                             previousA1 * std::pow(ktr, 3.0) * ((std::exp(-dt * ktr) * dt * (-std::pow(lambda1, 2.0) - lambda1 * lambda2 + 3.0 * lambda1 * ktr - std::pow(lambda2, 2.0) + 3.0 * lambda2 * ktr - 3.0 * std::pow(ktr, 2.0))) / (std::pow(lambda1 - ktr, 3.0) * std::pow(ktr - lambda2, 3.0)) +
                                                                 std::exp(-dt * ktr) * (-std::pow(lambda1, 3.0) - std::pow(lambda1, 2.0) * lambda2 + 4.0 * std::pow(lambda1, 2.0) * ktr - lambda1 * std::pow(lambda2, 2.0) + 4.0 * lambda1 * lambda2 * ktr - 6.0 * lambda1 * std::pow(ktr, 2.0) - std::pow(lambda2, 3.0) + 4.0 * std::pow(lambda2, 2.0) * ktr - 6.0 * lambda2 * std::pow(ktr, 2.0) + 4.0 * std::pow(ktr, 3.0)) / (std::pow(lambda1 - ktr, 4.0) * std::pow(ktr - lambda2, 4.0)) -
                                                                 std::exp(-dt * ktr) * std::pow(dt, 3.0) / (6.0 * (lambda1 - ktr) * (ktr - lambda2)) +
                                                                 std::exp(-dt * ktr) * std::pow(dt, 2.0) * (-lambda1 - lambda2 + 2.0 * ktr) / (2.0 * std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) -
                                                                 std::exp(-dt * lambda1) / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 4.0)) +
                                                                 std::exp(-dt * lambda2) / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 4.0))));
      a1 = previousA1 * std::exp(-dt * ktr);
    }
    std::vector<double> amounts = {a1, a2, a3, a4, a5, a6};
    apply_dose_2comp(i, ev, p, model_ss, 6, amounts);
    a1 = amounts[0]; a2 = amounts[1]; a3 = amounts[2]; a4 = amounts[3]; a5 = amounts[4]; a6 = amounts[5];
    ipred[i] = conc_2comp(event_obs_cmp(obs_cmp, ev, i), a1, a2, a3, p, ev, i);
    pk_write_amounts(i, amounts);
  }
  return ipred;
}

// amounts[0]=A1, [1]=A2, [2]=A3, [3]=A4, [4]=A5, [5]=A6, [6]=A7
NumericVector solve_2comp_4transit(const SubjectEvents& ev, const PkParams& p,
                                   int obs_cmp, int model_ss) {
  const int n = ev.time.size();
  NumericVector ipred(n);
  const double ktr = effective_ktr(p);
  const double k20 = effective_ke(p);
  const double k23 = p.k23;
  const double k32 = p.k32;
  const double k30 = p.k30;
  double a1 = 0.0, a2 = 0.0, a3 = 0.0, a4 = 0.0, a5 = 0.0, a6 = 0.0, a7 = 0.0;
  for (int i = 0; i < n; ++i) {
    if (i > 0) {
      const double dt = ev.time[i] - ev.time[i - 1];
      double E2, E3, lambda1, lambda2;
      hybrid_2_oral(k20, k23, k32, k30, E2, E3, lambda1, lambda2);
      const double previousA1 = a1;
      const double previousA2 = a2;
      const double previousA3 = a3;
      const double previousA4 = a4;
      const double previousA5 = a5;
      const double previousA6 = a6;
      const double previousA7 = a7;
      a4 = previousA4 * std::exp(-dt * ktr) + ktr * previousA1 * dt * std::exp(-dt * ktr);
      a5 = previousA5 * std::exp(-dt * ktr) + ktr * previousA4 * dt * std::exp(-dt * ktr) +
           0.5 * std::pow(ktr, 2.0) * previousA1 * std::pow(dt, 2.0) * std::exp(-dt * ktr);
      a6 = previousA6 * std::exp(-dt * ktr) + ktr * previousA5 * dt * std::exp(-dt * ktr) +
           0.5 * std::pow(ktr, 2.0) * previousA4 * std::pow(dt, 2.0) * std::exp(-dt * ktr) +
           (1.0 / 6.0) * std::pow(ktr, 3.0) * previousA1 * std::pow(dt, 3.0) * std::exp(-dt * ktr);
      a7 = previousA7 * std::exp(-dt * ktr) + ktr * previousA6 * dt * std::exp(-dt * ktr) +
           0.5 * std::pow(ktr, 2.0) * previousA5 * std::pow(dt, 2.0) * std::exp(-dt * ktr) +
           (1.0 / 6.0) * std::pow(ktr, 3.0) * previousA4 * std::pow(dt, 3.0) * std::exp(-dt * ktr) +
           (1.0 / 24.0) * std::pow(ktr, 4.0) * previousA1 * std::pow(dt, 4.0) * std::exp(-dt * ktr);
      a2 = (std::exp(-dt * lambda1) * ((previousA2 * E3 + previousA3 * k32) - previousA2 * lambda1) -
            std::exp(-dt * lambda2) * ((previousA2 * E3 + previousA3 * k32) - previousA2 * lambda2)) / (lambda2 - lambda1);
      a2 = a2 + ktr * k32 * (previousA7 * (std::exp(-dt * ktr) / ((lambda1 - ktr) * (lambda2 - ktr)) +
                                            std::exp(-dt * lambda1) / ((ktr - lambda1) * (lambda2 - lambda1)) +
                                            std::exp(-dt * lambda2) / ((ktr - lambda2) * (lambda1 - lambda2))) +
                              previousA6 * ktr * (std::exp(-dt * ktr) * (-lambda1 - lambda2 + 2.0 * ktr) / (std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) -
                                                  std::exp(-dt * lambda1) / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 2.0)) +
                                                  std::exp(-dt * lambda2) / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 2.0)) -
                                                  std::exp(-dt * ktr) * dt / ((lambda1 - ktr) * (ktr - lambda2))) +
                              previousA5 * std::pow(ktr, 2.0) * ((std::exp(-dt * ktr) * (-std::pow(lambda1, 2.0) - lambda1 * lambda2 + 3.0 * lambda1 * ktr - std::pow(lambda2, 2.0) + 3.0 * lambda2 * ktr - 3.0 * std::pow(ktr, 2.0))) / (std::pow(lambda1 - ktr, 3.0) * std::pow(ktr - lambda2, 3.0)) -
                                                                 std::exp(-dt * ktr) * std::pow(dt, 2.0) / (2.0 * (lambda1 - ktr) * (ktr - lambda2)) +
                                                                 std::exp(-dt * ktr) * dt * (-lambda1 - lambda2 + 2.0 * ktr) / (std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) +
                                                                 std::exp(-dt * lambda1) / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 3.0)) -
                                                                 std::exp(-dt * lambda2) / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 3.0))) +
                              previousA4 * std::pow(ktr, 3.0) * ((std::exp(-dt * ktr) * dt * (-std::pow(lambda1, 2.0) - lambda1 * lambda2 + 3.0 * lambda1 * ktr - std::pow(lambda2, 2.0) + 3.0 * lambda2 * ktr - 3.0 * std::pow(ktr, 2.0))) / (std::pow(lambda1 - ktr, 3.0) * std::pow(ktr - lambda2, 3.0)) +
                                                                 std::exp(-dt * ktr) * (-std::pow(lambda1, 3.0) - std::pow(lambda1, 2.0) * lambda2 + 4.0 * std::pow(lambda1, 2.0) * ktr - lambda1 * std::pow(lambda2, 2.0) + 4.0 * lambda1 * lambda2 * ktr - 6.0 * lambda1 * std::pow(ktr, 2.0) - std::pow(lambda2, 3.0) + 4.0 * std::pow(lambda2, 2.0) * ktr - 6.0 * lambda2 * std::pow(ktr, 2.0) + 4.0 * std::pow(ktr, 3.0)) / (std::pow(lambda1 - ktr, 4.0) * std::pow(ktr - lambda2, 4.0)) -
                                                                 std::exp(-dt * ktr) * std::pow(dt, 3.0) / (6.0 * (lambda1 - ktr) * (ktr - lambda2)) +
                                                                 std::exp(-dt * ktr) * std::pow(dt, 2.0) * (-lambda1 - lambda2 + 2.0 * ktr) / (2.0 * std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) -
                                                                 std::exp(-dt * lambda1) / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 4.0)) +
                                                                 std::exp(-dt * lambda2) / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 4.0))) +
                              previousA1 * std::pow(ktr, 4.0) * ((std::exp(-dt * ktr) * std::pow(dt, 2.0) * (-std::pow(lambda1, 2.0) - lambda1 * lambda2 + 3.0 * lambda1 * ktr - std::pow(lambda2, 2.0) + 3.0 * lambda2 * ktr - 3.0 * std::pow(ktr, 2.0))) / (2.0 * std::pow(lambda1 - ktr, 3.0) * std::pow(ktr - lambda2, 3.0)) +
                                                                 std::exp(-dt * ktr) * dt * (-std::pow(lambda1, 3.0) - std::pow(lambda1, 2.0) * lambda2 + 4.0 * std::pow(lambda1, 2.0) * ktr - lambda1 * std::pow(lambda2, 2.0) + 4.0 * lambda1 * lambda2 * ktr - 6.0 * lambda1 * std::pow(ktr, 2.0) - std::pow(lambda2, 3.0) + 4.0 * std::pow(lambda2, 2.0) * ktr - 6.0 * lambda2 * std::pow(ktr, 2.0) + 4.0 * std::pow(ktr, 3.0)) / (std::pow(lambda1 - ktr, 4.0) * std::pow(ktr - lambda2, 4.0)) +
                                                                 (std::exp(-dt * ktr) / (std::pow(lambda1 - ktr, 5.0) * std::pow(ktr - lambda2, 5.0))) * (-std::pow(lambda1, 4.0) - std::pow(lambda1, 3.0) * lambda2 + 5.0 * std::pow(lambda1, 3.0) * ktr - std::pow(lambda1, 2.0) * std::pow(lambda2, 2.0) + 5.0 * std::pow(lambda1, 2.0) * lambda2 * ktr - 10.0 * std::pow(lambda1, 2.0) * std::pow(ktr, 2.0) - lambda1 * std::pow(lambda2, 3.0) + 5.0 * lambda1 * std::pow(lambda2, 2.0) * ktr - 10.0 * lambda1 * lambda2 * std::pow(ktr, 2.0) + 10.0 * lambda1 * std::pow(ktr, 3.0) - std::pow(lambda2, 4.0) + 5.0 * std::pow(lambda2, 3.0) * ktr - 10.0 * std::pow(lambda2, 2.0) * std::pow(ktr, 2.0) + 10.0 * lambda2 * std::pow(ktr, 3.0) - 5.0 * std::pow(ktr, 4.0)) -
                                                                 std::exp(-dt * ktr) * std::pow(dt, 4.0) / (24.0 * (lambda1 - ktr) * (ktr - lambda2)) +
                                                                 std::exp(-dt * ktr) * std::pow(dt, 3.0) * (-lambda1 - lambda2 + 2.0 * ktr) / (6.0 * std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) +
                                                                 std::exp(-dt * lambda1) / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 5.0)) -
                                                                 std::exp(-dt * lambda2) / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 5.0))));
      a2 = a2 + ktr * (previousA7 * (std::exp(-dt * ktr) * ktr / ((lambda1 - ktr) * (ktr - lambda2)) +
                                     std::exp(-dt * lambda2) * lambda2 / ((lambda1 - lambda2) * (lambda2 - ktr)) -
                                     std::exp(-dt * lambda1) * lambda1 / ((lambda1 - lambda2) * (lambda1 - ktr))) +
                       previousA6 * ktr * (std::exp(-dt * ktr) * (lambda1 * lambda2 - std::pow(ktr, 2.0)) / (std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) +
                                           std::exp(-dt * lambda1) * lambda1 / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 2.0)) -
                                           std::exp(-dt * lambda2) * lambda2 / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 2.0)) +
                                           std::exp(-dt * ktr) * dt * ktr / ((lambda1 - ktr) * (ktr - lambda2))) +
                       previousA5 * std::pow(ktr, 2.0) * (std::exp(-dt * ktr) * (std::pow(lambda1, 2.0) * lambda2 + lambda1 * std::pow(lambda2, 2.0) - 3.0 * lambda1 * lambda2 * ktr + std::pow(ktr, 3.0)) / (std::pow(lambda1 - ktr, 3.0) * std::pow(ktr - lambda2, 3.0)) +
                                                           std::exp(-dt * ktr) * dt * (lambda1 * lambda2 - std::pow(ktr, 2.0)) / (std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) +
                                                           std::exp(-dt * ktr) * ktr * std::pow(dt, 2.0) / (2.0 * (lambda1 - ktr) * (ktr - lambda2)) -
                                                           std::exp(-dt * lambda1) * lambda1 / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 3.0)) +
                                                           std::exp(-dt * lambda2) * lambda2 / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 3.0))) +
                       previousA4 * std::pow(ktr, 3.0) * (std::exp(-dt * ktr) * dt * (std::pow(lambda1, 2.0) * lambda2 + lambda1 * std::pow(lambda2, 2.0) - 3.0 * lambda1 * lambda2 * ktr + std::pow(ktr, 3.0)) / (std::pow(lambda1 - ktr, 3.0) * std::pow(ktr - lambda2, 3.0)) +
                                                           std::exp(-dt * ktr) * (std::pow(lambda1, 3.0) * lambda2 + std::pow(lambda1, 2.0) * std::pow(lambda2, 2.0) - 4.0 * std::pow(lambda1, 2.0) * lambda2 * ktr + lambda1 * std::pow(lambda2, 3.0) - 4.0 * lambda1 * std::pow(lambda2, 2.0) * ktr + 6.0 * lambda1 * lambda2 * std::pow(ktr, 2.0) - std::pow(ktr, 4.0)) / (std::pow(lambda1 - ktr, 4.0) * std::pow(ktr - lambda2, 4.0)) +
                                                           std::exp(-dt * ktr) * std::pow(dt, 2.0) * (lambda1 * lambda2 - std::pow(ktr, 2.0)) / (2.0 * std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) +
                                                           std::exp(-dt * ktr) * ktr * std::pow(dt, 3.0) / (6.0 * (lambda1 - ktr) * (ktr - lambda2)) +
                                                           std::exp(-dt * lambda1) * lambda1 / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 4.0)) -
                                                           std::exp(-dt * lambda2) * lambda2 / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 4.0))) +
                       previousA1 * std::pow(ktr, 4.0) * ((std::exp(-dt * ktr) * std::pow(dt, 2.0) * (std::pow(lambda1, 2.0) * lambda2 + lambda1 * std::pow(lambda2, 2.0) - 3.0 * lambda1 * lambda2 * ktr + std::pow(ktr, 3.0))) / (2.0 * std::pow(lambda1 - ktr, 3.0) * std::pow(ktr - lambda2, 3.0)) +
                                                           std::exp(-dt * ktr) * dt * (std::pow(lambda1, 3.0) * lambda2 + std::pow(lambda1, 2.0) * std::pow(lambda2, 2.0) - 4.0 * std::pow(lambda1, 2.0) * lambda2 * ktr + lambda1 * std::pow(lambda2, 3.0) - 4.0 * lambda1 * std::pow(lambda2, 2.0) * ktr + 6.0 * lambda1 * lambda2 * std::pow(ktr, 2.0) - std::pow(ktr, 4.0)) / (std::pow(lambda1 - ktr, 4.0) * std::pow(ktr - lambda2, 4.0)) +
                                                           (std::exp(-dt * ktr) / (std::pow(lambda1 - ktr, 5.0) * std::pow(ktr - lambda2, 5.0))) * (std::pow(lambda1, 4.0) * lambda2 + std::pow(lambda1, 3.0) * std::pow(lambda2, 2.0) - 5.0 * std::pow(lambda1, 3.0) * lambda2 * ktr + std::pow(lambda1, 2.0) * std::pow(lambda2, 3.0) - 5.0 * std::pow(lambda1, 2.0) * std::pow(lambda2, 2.0) * ktr + 10.0 * std::pow(lambda1, 2.0) * lambda2 * std::pow(ktr, 2.0) + lambda1 * std::pow(lambda2, 4.0) - 5.0 * lambda1 * std::pow(lambda2, 3.0) * ktr + 10.0 * lambda1 * std::pow(lambda2, 2.0) * std::pow(ktr, 2.0) - 10.0 * lambda1 * lambda2 * std::pow(ktr, 3.0) + std::pow(ktr, 5.0)) +
                                                           std::exp(-dt * ktr) * std::pow(dt, 3.0) * (lambda1 * lambda2 - std::pow(ktr, 2.0)) / (6.0 * std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) +
                                                           std::exp(-dt * ktr) * ktr * std::pow(dt, 4.0) / (24.0 * (lambda1 - ktr) * (ktr - lambda2)) -
                                                           std::exp(-dt * lambda1) * lambda1 / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 5.0)) +
                                                           std::exp(-dt * lambda2) * lambda2 / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 5.0))));
      a3 = (std::exp(-dt * lambda1) * ((previousA3 * E2 + k23 * previousA2) - previousA3 * lambda1) -
            std::exp(-dt * lambda2) * ((previousA3 * E2 + k23 * previousA2) - previousA3 * lambda2)) / (lambda2 - lambda1);
      a3 = a3 + ktr * k23 * (previousA7 * (std::exp(-dt * ktr) / ((lambda1 - ktr) * (lambda2 - ktr)) +
                                            std::exp(-dt * lambda1) / ((ktr - lambda1) * (lambda2 - lambda1)) +
                                            std::exp(-dt * lambda2) / ((ktr - lambda2) * (lambda1 - lambda2))) +
                             previousA6 * ktr * (std::exp(-dt * ktr) * (-lambda1 - lambda2 + 2.0 * ktr) / (std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) -
                                                 std::exp(-dt * lambda1) / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 2.0)) +
                                                 std::exp(-dt * lambda2) / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 2.0)) -
                                                 std::exp(-dt * ktr) * dt / ((lambda1 - ktr) * (ktr - lambda2))) +
                             previousA5 * std::pow(ktr, 2.0) * ((std::exp(-dt * ktr) * (-std::pow(lambda1, 2.0) - lambda1 * lambda2 + 3.0 * lambda1 * ktr - std::pow(lambda2, 2.0) + 3.0 * lambda2 * ktr - 3.0 * std::pow(ktr, 2.0))) / (std::pow(lambda1 - ktr, 3.0) * std::pow(ktr - lambda2, 3.0)) -
                                                                 std::exp(-dt * ktr) * std::pow(dt, 2.0) / (2.0 * (lambda1 - ktr) * (ktr - lambda2)) +
                                                                 std::exp(-dt * ktr) * dt * (-lambda1 - lambda2 + 2.0 * ktr) / (std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) +
                                                                 std::exp(-dt * lambda1) / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 3.0)) -
                                                                 std::exp(-dt * lambda2) / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 3.0))) +
                             previousA4 * std::pow(ktr, 3.0) * ((std::exp(-dt * ktr) * dt * (-std::pow(lambda1, 2.0) - lambda1 * lambda2 + 3.0 * lambda1 * ktr - std::pow(lambda2, 2.0) + 3.0 * lambda2 * ktr - 3.0 * std::pow(ktr, 2.0))) / (std::pow(lambda1 - ktr, 3.0) * std::pow(ktr - lambda2, 3.0)) +
                                                                 std::exp(-dt * ktr) * (-std::pow(lambda1, 3.0) - std::pow(lambda1, 2.0) * lambda2 + 4.0 * std::pow(lambda1, 2.0) * ktr - lambda1 * std::pow(lambda2, 2.0) + 4.0 * lambda1 * lambda2 * ktr - 6.0 * lambda1 * std::pow(ktr, 2.0) - std::pow(lambda2, 3.0) + 4.0 * std::pow(lambda2, 2.0) * ktr - 6.0 * lambda2 * std::pow(ktr, 2.0) + 4.0 * std::pow(ktr, 3.0)) / (std::pow(lambda1 - ktr, 4.0) * std::pow(ktr - lambda2, 4.0)) -
                                                                 std::exp(-dt * ktr) * std::pow(dt, 3.0) / (6.0 * (lambda1 - ktr) * (ktr - lambda2)) +
                                                                 std::exp(-dt * ktr) * std::pow(dt, 2.0) * (-lambda1 - lambda2 + 2.0 * ktr) / (2.0 * std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) -
                                                                 std::exp(-dt * lambda1) / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 4.0)) +
                                                                 std::exp(-dt * lambda2) / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 4.0))) +
                             previousA1 * std::pow(ktr, 4.0) * ((std::exp(-dt * ktr) * std::pow(dt, 2.0) * (-std::pow(lambda1, 2.0) - lambda1 * lambda2 + 3.0 * lambda1 * ktr - std::pow(lambda2, 2.0) + 3.0 * lambda2 * ktr - 3.0 * std::pow(ktr, 2.0))) / (2.0 * std::pow(lambda1 - ktr, 3.0) * std::pow(ktr - lambda2, 3.0)) +
                                                                 std::exp(-dt * ktr) * dt * (-std::pow(lambda1, 3.0) - std::pow(lambda1, 2.0) * lambda2 + 4.0 * std::pow(lambda1, 2.0) * ktr - lambda1 * std::pow(lambda2, 2.0) + 4.0 * lambda1 * lambda2 * ktr - 6.0 * lambda1 * std::pow(ktr, 2.0) - std::pow(lambda2, 3.0) + 4.0 * std::pow(lambda2, 2.0) * ktr - 6.0 * lambda2 * std::pow(ktr, 2.0) + 4.0 * std::pow(ktr, 3.0)) / (std::pow(lambda1 - ktr, 4.0) * std::pow(ktr - lambda2, 4.0)) +
                                                                 (std::exp(-dt * ktr) / (std::pow(lambda1 - ktr, 5.0) * std::pow(ktr - lambda2, 5.0))) * (-std::pow(lambda1, 4.0) - std::pow(lambda1, 3.0) * lambda2 + 5.0 * std::pow(lambda1, 3.0) * ktr - std::pow(lambda1, 2.0) * std::pow(lambda2, 2.0) + 5.0 * std::pow(lambda1, 2.0) * lambda2 * ktr - 10.0 * std::pow(lambda1, 2.0) * std::pow(ktr, 2.0) - lambda1 * std::pow(lambda2, 3.0) + 5.0 * lambda1 * std::pow(lambda2, 2.0) * ktr - 10.0 * lambda1 * lambda2 * std::pow(ktr, 2.0) + 10.0 * lambda1 * std::pow(ktr, 3.0) - std::pow(lambda2, 4.0) + 5.0 * std::pow(lambda2, 3.0) * ktr - 10.0 * std::pow(lambda2, 2.0) * std::pow(ktr, 2.0) + 10.0 * lambda2 * std::pow(ktr, 3.0) - 5.0 * std::pow(ktr, 4.0)) -
                                                                 std::exp(-dt * ktr) * std::pow(dt, 4.0) / (24.0 * (lambda1 - ktr) * (ktr - lambda2)) +
                                                                 std::exp(-dt * ktr) * std::pow(dt, 3.0) * (-lambda1 - lambda2 + 2.0 * ktr) / (6.0 * std::pow(lambda1 - ktr, 2.0) * std::pow(ktr - lambda2, 2.0)) +
                                                                 std::exp(-dt * lambda1) / ((lambda1 - lambda2) * std::pow(lambda1 - ktr, 5.0)) -
                                                                 std::exp(-dt * lambda2) / ((lambda1 - lambda2) * std::pow(lambda2 - ktr, 5.0))));
      a1 = previousA1 * std::exp(-dt * ktr);
    }
    std::vector<double> amounts = {a1, a2, a3, a4, a5, a6, a7};
    apply_dose_2comp(i, ev, p, model_ss, 7, amounts);
    a1 = amounts[0]; a2 = amounts[1]; a3 = amounts[2]; a4 = amounts[3]; a5 = amounts[4]; a6 = amounts[5]; a7 = amounts[6];
    ipred[i] = conc_2comp(event_obs_cmp(obs_cmp, ev, i), a1, a2, a3, p, ev, i);
    pk_write_amounts(i, amounts);
  }
  return ipred;
}

}  // namespace

namespace nm_pk {

NumericVector solve_1comp_transit(const SubjectEvents& ev, const PkParams& p,
                                  int n_transit, int obs_cmp, int model_ss) {
  switch (n_transit) {
    case 1: return solve_1comp_1transit(ev, p, obs_cmp, model_ss);
    case 2: return solve_1comp_2transit(ev, p, obs_cmp, model_ss);
    case 3: return solve_1comp_3transit(ev, p, obs_cmp, model_ss);
    case 4: return solve_1comp_4transit(ev, p, obs_cmp, model_ss);
    default: return solve_1comp_1transit(ev, p, obs_cmp, model_ss);
  }
}

NumericVector solve_2comp_transit(const SubjectEvents& ev, const PkParams& p,
                                  int n_transit, int obs_cmp, int model_ss) {
  switch (n_transit) {
    case 1: return solve_2comp_1transit(ev, p, obs_cmp, model_ss);
    case 2: return solve_2comp_2transit(ev, p, obs_cmp, model_ss);
    case 3: return solve_2comp_3transit(ev, p, obs_cmp, model_ss);
    case 4: return solve_2comp_4transit(ev, p, obs_cmp, model_ss);
    default: return solve_2comp_1transit(ev, p, obs_cmp, model_ss);
  }
}

}  // namespace nm_pk
