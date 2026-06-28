#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include "nm_pk_core.h"
#include "nm_pk_pkadvan.h"
using namespace Rcpp;
using namespace nm_pk;
using namespace pkadvan;

namespace nm_pk {
thread_local int pk_route_dose_cmp = 1;
}

namespace {

NumericVector solve_1_iv_bolus(const SubjectEvents& ev, const PkParams& p,
                               int obs_cmp, int model_ss, int advan) {
  const int n = ev.time.size();
  NumericVector ipred(n);
  double a1 = 0.0;
  const double k10 = p.k10 > 0.0 ? p.k10 : effective_ke(p);
  std::vector<double> amounts(1, 0.0);
  for (int i = 0; i < n; ++i) {
    if (i > 0) {
      double dt = ev.time[i] - ev.time[i - 1];
      step_1_iv_bolus(dt, k10, a1, a1);
    }
    amounts[0] = a1;
    apply_dose_with_ss(model_ss, advan, ROUTE_IV_BOLUS, ev, i, p, amounts, false);
    a1 = amounts[0];
    ipred[i] = conc_iv(event_obs_cmp(obs_cmp, ev, i), a1, 0.0, 0.0, p, ev, i);
    pk_write_amounts(i, amounts);
  }
  return ipred;
}

NumericVector solve_1_iv_inf(const SubjectEvents& ev, const PkParams& p,
                             int obs_cmp, int model_ss, int advan) {
  const int n = ev.time.size();
  NumericVector ipred(n);
  double a1 = 0.0;
  const double k10 = p.k10 > 0.0 ? p.k10 : effective_ke(p);
  std::vector<double> amounts(1, 0.0);
  for (int i = 0; i < n; ++i) {
    if (i > 0) {
      double dt = ev.time[i] - ev.time[i - 1];
      step_1_iv_infusion(dt, k10, effective_rate(ev, i), a1, a1);
    }
    amounts[0] = a1;
    apply_dose_with_ss(model_ss, advan, ROUTE_IV_INFUSION, ev, i, p, amounts, false);
    a1 = amounts[0];
    ipred[i] = conc_iv(event_obs_cmp(obs_cmp, ev, i), a1, 0.0, 0.0, p, ev, i);
    pk_write_amounts(i, amounts);
  }
  return ipred;
}

bool ss_oral_1_alag_cycle(const PkParams& p, double dose_ff, double tau, double lag,
                          double& gut, double& cen) {
  if (lag <= 0.0 || tau <= lag + 1e-12) return false;
  const double k10 = p.k10 > 0.0 ? p.k10 : effective_ke(p);
  gut = 0.0;
  cen = 0.0;
  for (int iter = 0; iter < 100; ++iter) {
    const double g0 = gut, c0 = cen;
    step_1_oral(lag, p.ka, k10, gut, cen, gut, cen);
    gut += dose_ff;
    step_1_oral(tau - lag, p.ka, k10, gut, cen, gut, cen);
    if (std::fabs(gut - g0) < 1e-8 && std::fabs(cen - c0) < 1e-8) return true;
  }
  return true;
}

NumericVector solve_1_oral_alag(const SubjectEvents& ev, const PkParams& p,
                              int obs_cmp, int model_ss, int advan) {
  const int n = ev.time.size();
  NumericVector ipred(n);
  double gut = 0.0, cen = 0.0;
  const double k10 = p.k10 > 0.0 ? p.k10 : effective_ke(p);
  struct PendingDose { double release_time; double amount; };
  std::vector<PendingDose> pending;
  std::vector<double> amounts(2, 0.0);
  int n_doses = 0;

  auto release_pending = [&](double t) {
    for (size_t j = 0; j < pending.size(); ) {
      if (pending[j].release_time <= t + 1e-9) {
        gut += pending[j].amount;
        pending.erase(pending.begin() + static_cast<long>(j));
      } else {
        ++j;
      }
    }
  };

  for (int i = 0; i < n; ++i) {
    const double t = ev.time[i];
    if (i > 0) {
      double t0 = ev.time[i - 1];
      double t1 = t;
      std::vector<double> breaks;
      breaks.push_back(t1);
      for (const auto& pd : pending) {
        if (pd.release_time > t0 + 1e-9 && pd.release_time < t1 - 1e-9) {
          breaks.push_back(pd.release_time);
        }
      }
      std::sort(breaks.begin(), breaks.end());
      double cur = t0;
      for (double tb : breaks) {
        if (tb > cur + 1e-9) {
          step_1_oral(tb - cur, p.ka, k10, gut, cen, gut, cen);
          cur = tb;
        }
        if (tb < t1 - 1e-9) {
          release_pending(cur);
        }
      }
    }
    amounts[0] = gut;
    amounts[1] = cen;
    if (is_reset_row(ev, i)) {
      gut = 0.0;
      cen = 0.0;
      pending.clear();
      n_doses = 0;
    } else if (ss_row_active(model_ss, row_ss(ev, i)) && row_ii(ev, i) > 0.0 &&
               p.alag1 > 0.0 && is_dosing_row(ev, i) && ev.amt[i] > 0.0) {
      if (ss_row_reset(row_ss(ev, i))) {
        gut = 0.0;
        cen = 0.0;
      }
      const double dose = ev.amt[i] * row_f1(ev, i, p);
      const double tau = row_ii(ev, i);
      ss_oral_1_alag_cycle(p, dose, tau, p.alag1, gut, cen);
      n_doses++;
      pending.push_back({t + p.alag1, dose});
    } else if (apply_ss_at_dose(model_ss, row_ss(ev, i), advan, ROUTE_ORAL, p, ev, i, amounts,
                               false)) {
      gut = amounts[0];
      cen = amounts[1];
      if (p.alag1 > 0.0 && is_dosing_row(ev, i) && ev.amt[i] > 0.0) {
        const double dose = ev.amt[i] * row_f1(ev, i, p);
        n_doses++;
        pending.push_back({t + p.alag1, dose});
      }
    } else if (is_dosing_row(ev, i) && ev.amt[i] > 0.0) {
      const double dose = ev.amt[i] * row_f1(ev, i, p);
      n_doses++;
      pending.push_back({t + p.alag1, dose});
    }
    ipred[i] = conc_oral(event_obs_cmp(obs_cmp, ev, i), gut, cen, 0.0, 0.0, p, ev, i);
    bool lag_release_row = false;
    if (p.alag1 > 0.0) {
      for (const auto& pd : pending) {
        if (std::fabs(pd.release_time - t) <= 1e-9) {
          lag_release_row = true;
          break;
        }
      }
    }
    const bool post_release_amt = lag_release_row && n_doses >= 2;
    if (post_release_amt) {
      release_pending(t);
      pk_write_amounts(i, std::vector<double>{gut, cen});
    } else {
      pk_write_amounts(i, std::vector<double>{gut, cen});
      release_pending(t);
    }
  }
  return ipred;
}

NumericVector solve_1_oral(const SubjectEvents& ev, const PkParams& p,
                           int obs_cmp, int model_ss, int advan) {
  if (p.alag1 > 0.0) {
    return solve_1_oral_alag(ev, p, obs_cmp, model_ss, advan);
  }
  const int n = ev.time.size();
  NumericVector ipred(n);
  double gut = 0.0, cen = 0.0;
  const double k10 = p.k10 > 0.0 ? p.k10 : effective_ke(p);
  std::vector<double> amounts(2, 0.0);
  for (int i = 0; i < n; ++i) {
    if (i > 0) {
      double dt = ev.time[i] - ev.time[i - 1];
      step_1_oral(dt, p.ka, k10, gut, cen, gut, cen);
    }
    amounts[0] = gut;
    amounts[1] = cen;
    apply_dose_with_ss(model_ss, advan, ROUTE_ORAL, ev, i, p, amounts, true);
    gut = amounts[0];
    cen = amounts[1];
    ipred[i] = conc_oral(event_obs_cmp(obs_cmp, ev, i), gut, cen, 0.0, 0.0, p, ev, i);
    pk_write_amounts(i, amounts);
  }
  return ipred;
}

NumericVector solve_2_iv_bolus(const SubjectEvents& ev, const PkParams& p,
                               int obs_cmp, int model_ss, int advan) {
  const int n = ev.time.size();
  NumericVector ipred(n);
  double a1 = 0.0, a2 = 0.0;
  std::vector<double> amounts(2, 0.0);
  for (int i = 0; i < n; ++i) {
    if (i > 0) {
      double dt = ev.time[i] - ev.time[i - 1];
      step_2_iv_bolus(dt, p.k10, p.k12, p.k21, a1, a2, a1, a2);
    }
    amounts[0] = a1;
    amounts[1] = a2;
    apply_dose_with_ss(model_ss, advan, ROUTE_IV_BOLUS, ev, i, p, amounts, false);
    a1 = amounts[0];
    a2 = amounts[1];
    ipred[i] = conc_iv(event_obs_cmp(obs_cmp, ev, i), a1, a2, 0.0, p, ev, i);
    pk_write_amounts(i, amounts);
  }
  return ipred;
}

NumericVector solve_2_iv_inf(const SubjectEvents& ev, const PkParams& p,
                             int obs_cmp, int model_ss, int advan) {
  const int n = ev.time.size();
  NumericVector ipred(n);
  double a1 = 0.0, a2 = 0.0;
  std::vector<double> amounts(2, 0.0);
  for (int i = 0; i < n; ++i) {
    if (i > 0) {
      double dt = ev.time[i] - ev.time[i - 1];
      step_2_iv_infusion(dt, p.k10, p.k12, p.k21, effective_rate(ev, i), a1, a2, a1, a2);
    }
    amounts[0] = a1;
    amounts[1] = a2;
    apply_dose_with_ss(model_ss, advan, ROUTE_IV_INFUSION, ev, i, p, amounts, false);
    a1 = amounts[0];
    a2 = amounts[1];
    ipred[i] = conc_iv(event_obs_cmp(obs_cmp, ev, i), a1, a2, 0.0, p, ev, i);
    pk_write_amounts(i, amounts);
  }
  return ipred;
}

NumericVector solve_2_oral(const SubjectEvents& ev, const PkParams& p,
                           int obs_cmp, int model_ss, int advan) {
  const int n = ev.time.size();
  NumericVector ipred(n);
  double gut = 0.0, a2 = 0.0, a3 = 0.0;
  std::vector<double> amounts(3, 0.0);
  for (int i = 0; i < n; ++i) {
    if (i > 0) {
      double dt = ev.time[i] - ev.time[i - 1];
      step_2_oral(dt, p.ka, p.k20, p.k23, p.k32, gut, a2, a3, gut, a2, a3);
    }
    amounts[0] = gut;
    amounts[1] = a2;
    amounts[2] = a3;
    apply_dose_with_ss(model_ss, advan, ROUTE_ORAL, ev, i, p, amounts, true);
    gut = amounts[0];
    a2 = amounts[1];
    a3 = amounts[2];
    ipred[i] = conc_oral(event_obs_cmp(obs_cmp, ev, i), gut, a2, a3, 0.0, p, ev, i);
    pk_write_amounts(i, amounts);
  }
  return ipred;
}

NumericVector solve_3_iv_bolus(const SubjectEvents& ev, const PkParams& p,
                               int obs_cmp, int model_ss, int advan) {
  const int n = ev.time.size();
  NumericVector ipred(n);
  double a1 = 0.0, a2 = 0.0, a3 = 0.0;
  std::vector<double> amounts(3, 0.0);
  for (int i = 0; i < n; ++i) {
    if (i > 0) {
      double dt = ev.time[i] - ev.time[i - 1];
      step_3_iv_bolus(dt, p.k10, p.k12, p.k21, p.k13, p.k31, a1, a2, a3, a1, a2, a3);
    }
    amounts[0] = a1;
    amounts[1] = a2;
    amounts[2] = a3;
    apply_dose_with_ss(model_ss, advan, ROUTE_IV_BOLUS, ev, i, p, amounts, false);
    a1 = amounts[0];
    a2 = amounts[1];
    a3 = amounts[2];
    ipred[i] = conc_iv(event_obs_cmp(obs_cmp, ev, i), a1, a2, a3, p, ev, i);
    pk_write_amounts(i, amounts);
  }
  return ipred;
}

NumericVector solve_3_iv_inf(const SubjectEvents& ev, const PkParams& p,
                             int obs_cmp, int model_ss, int advan) {
  const int n = ev.time.size();
  NumericVector ipred(n);
  double a1 = 0.0, a2 = 0.0, a3 = 0.0;
  std::vector<double> amounts(3, 0.0);
  for (int i = 0; i < n; ++i) {
    if (i > 0) {
      double dt = ev.time[i] - ev.time[i - 1];
      step_3_iv_infusion(dt, p.k10, p.k12, p.k21, p.k13, p.k31,
                         effective_rate(ev, i), a1, a2, a3, a1, a2, a3);
    }
    amounts[0] = a1;
    amounts[1] = a2;
    amounts[2] = a3;
    apply_dose_with_ss(model_ss, advan, ROUTE_IV_INFUSION, ev, i, p, amounts, false);
    a1 = amounts[0];
    a2 = amounts[1];
    a3 = amounts[2];
    ipred[i] = conc_iv(event_obs_cmp(obs_cmp, ev, i), a1, a2, a3, p, ev, i);
    pk_write_amounts(i, amounts);
  }
  return ipred;
}

NumericVector solve_3_oral(const SubjectEvents& ev, const PkParams& p,
                           int obs_cmp, int model_ss, int advan, int trans) {
  const int n = ev.time.size();
  NumericVector ipred(n);
  double gut = 0.0, a2 = 0.0, a3 = 0.0, a4 = 0.0;
  std::vector<double> amounts(4, 0.0);
  for (int i = 0; i < n; ++i) {
    if (i > 0) {
      double dt = ev.time[i] - ev.time[i - 1];
      step_3_oral(dt, p.ka, p.k20, p.k23, p.k32, p.k24, p.k42,
                  gut, a2, a3, a4, gut, a2, a3, a4);
    }
    amounts[0] = gut;
    amounts[1] = a2;
    amounts[2] = a3;
    amounts[3] = a4;
    apply_dose_with_ss(model_ss, advan, ROUTE_ORAL, ev, i, p, amounts, true);
    gut = amounts[0];
    a2 = amounts[1];
    a3 = amounts[2];
    a4 = amounts[3];
    ipred[i] = conc_oral(event_obs_cmp(obs_cmp, ev, i), gut, a2, a3, a4, p, ev, i);
    pk_write_amounts(i, amounts);
  }
  return ipred;
}

// Beal (1983): t = (A0-A)/Vm + (Km/Vm)*ln(A0/A) for dA/dt = -Vm*A/(Km+A).
double step_mm_iv_elim(double a0, double dt, double vm, double km) {
  if (dt <= 0.0 || a0 <= 0.0) return 0.0;
  if (vm <= 0.0) return a0;

  const double inv_vm = 1.0 / vm;
  const double km_inv_vm = km * inv_vm;
  auto residual = [&](double a) {
    if (a <= 0.0) return 1e100;
    return (a0 - a) * inv_vm + km_inv_vm * std::log(a0 / a) - dt;
  };

  double hi = a0;
  if (residual(hi) >= 0.0) return hi;

  double lo = std::min(a0 * 1e-12, km * 1e-12);
  lo = std::max(lo, 1e-300);
  while (lo > 0.0 && residual(lo) <= 0.0) {
    lo *= 0.1;
    if (lo < 1e-300) return 0.0;
  }

  for (int iter = 0; iter < 80; ++iter) {
    const double mid = 0.5 * (lo + hi);
    if (residual(mid) > 0.0) {
      lo = mid;
    } else {
      hi = mid;
    }
    if (hi - lo <= std::max(1e-15 * a0, 1e-300)) break;
  }
  return hi;
}

NumericVector solve_mm_iv(const SubjectEvents& ev, const PkParams& p,
                          int obs_cmp, int model_ss) {
  const int nt = ev.time.size();
  NumericVector ipred(nt);
  // ADVAN10: VM is maximum elimination rate; KM is in amount units (not concentration).
  const double mm_vm = p.vm > 0.0 ? p.vm : (p.vmax > 0.0 ? p.vmax : 1.0);
  const double mm_km = p.km > 0.0 ? p.km : 1.0;
  std::vector<double> a(1, 0.0);
  for (int i = 0; i < nt; ++i) {
    if (i > 0) {
      const double dt = ev.time[i] - ev.time[i - 1];
      a[0] = step_mm_iv_elim(a[0], dt, mm_vm, mm_km);
    }
    if (apply_ss_at_dose(model_ss, row_ss(ev, i), 1, ROUTE_IV_BOLUS, p, ev, i, a)) {
      // SS state initialized
    } else if (is_reset_row(ev, i)) {
      a[0] = 0.0;
    } else if (is_dosing_row(ev, i)) {
      const int cmt = row_cmt(ev, i, pk_route_dose_cmp);
      const double dose = ev.amt[i] * row_f1(ev, i, p);
      if (dose > 0.0 && cmt == 1) a[0] += dose;
    }
    ipred[i] = a[0];
    pk_write_amounts(i, a);
  }
  return ipred;
}

NumericVector solve_ode(int advan, RouteType route, int obs_cmp,
                        const SubjectEvents& ev, const PkParams& p,
                        int model_ss) {
  int n = 1;
  std::vector<double> kmat;
  if (advan == 1) {
    n = 1;
    kmat = { -p.k10 };
    if (route == ROUTE_ORAL) {
      n = 2;
      kmat = { -p.ka, 0.0, p.ka, -(p.k10 > 0 ? p.k10 : p.k20) };
    }
  } else if (advan == 2 || advan == 4) {
    if (route == ROUTE_ORAL) {
      n = 2;
      kmat = { -p.ka, 0.0, p.ka, -p.k20, 0.0, 0.0 };
    } else {
      n = 1;
      kmat = { -(p.k10 > 0 ? p.k10 : p.k20) };
    }
  } else if (advan == 3 || advan == 11 || advan == 12) {
    if (route == ROUTE_ORAL) {
      n = 3;
      kmat = {
        -p.ka, 0.0, 0.0,
        p.ka, -(p.k20 + p.k23), p.k32,
        0.0, p.k23, -p.k32
      };
      if (advan == 11 || advan == 12) {
        n = 4;
        kmat = {
          -p.ka, 0.0, 0.0, 0.0,
          p.ka, -(p.k20 + p.k23 + p.k24), p.k32, p.k42,
          0.0, p.k23, -p.k32, 0.0,
          0.0, p.k24, 0.0, -p.k42
        };
      }
    } else {
      n = (advan == 11) ? 3 : 2;
      if (n == 2) {
        kmat = { -(p.k10 + p.k12), p.k21, p.k12, -p.k21 };
      } else {
        kmat = {
          -(p.k10 + p.k12 + p.k13), p.k21, p.k31,
          p.k12, -p.k21, 0.0,
          p.k13, 0.0, -p.k31
        };
      }
    }
  } else if (advan >= 6) {
    double k20 = p.k20 > 0.0 ? p.k20 :
      (p.cl > 0.0 && (p.vc > 0.0 || p.v2 > 0.0) ?
         p.cl / (p.vc > 0.0 ? p.vc : p.v2) :
         (p.k10 > 0.0 ? p.k10 : 0.0));
    double k23 = p.k23 > 0.0 ? p.k23 :
      (p.q2 > 0.0 && (p.vc > 0.0 || p.v2 > 0.0) ?
         p.q2 / (p.vc > 0.0 ? p.vc : p.v2) : 0.0);
    double k32 = p.k32 > 0.0 ? p.k32 :
      (p.q2 > 0.0 && (p.vp > 0.0 || p.v3 > 0.0) ?
         p.q2 / (p.vp > 0.0 ? p.vp : p.v3) : 0.0);
    if (route == ROUTE_ORAL || p.ka > 0.0) {
      if (k23 > 0.0 && k32 > 0.0) {
        n = 3;
        kmat = {
          -p.ka, 0.0, 0.0,
          p.ka, -(k20 + k23), k32,
          0.0, k23, -k32
        };
      } else if (p.ka > 0.0) {
        n = 2;
        kmat = { -p.ka, 0.0, p.ka, -k20 };
      } else if (k23 > 0.0 && k32 > 0.0) {
        n = 2;
        kmat = { -(k20 + k23), k32, k23, -k32 };
      } else {
        n = 1;
        kmat = { -k20 };
      }
    } else if (k23 > 0.0 && k32 > 0.0) {
      n = 2;
      kmat = { -(k20 + k23), k32, k23, -k32 };
    } else {
      n = 1;
      kmat = { -k20 };
    }
  } else {
    n = 2;
    kmat = { -(p.k10 + p.k12), p.k21, p.k12, -p.k21 };
  }

  auto rk4_step = [](const std::vector<double>& km, int nc, double dt,
                     std::vector<double>& a, std::vector<double>& out) {
    auto f = [&](const std::vector<double>& y, std::vector<double>& dy) {
      for (int i = 0; i < nc; ++i) dy[i] = 0.0;
      for (int i = 0; i < nc; ++i) {
        for (int j = 0; j < nc; ++j) dy[i] += km[i * nc + j] * y[j];
      }
    };
    std::vector<double> k1(nc), k2(nc), k3(nc), k4(nc), tmp(nc);
    f(a, k1);
    for (int i = 0; i < nc; ++i) tmp[i] = a[i] + 0.5 * dt * k1[i];
    f(tmp, k2);
    for (int i = 0; i < nc; ++i) tmp[i] = a[i] + 0.5 * dt * k2[i];
    f(tmp, k3);
    for (int i = 0; i < nc; ++i) tmp[i] = a[i] + dt * k3[i];
    f(tmp, k4);
    for (int i = 0; i < nc; ++i) {
      out[i] = a[i] + dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
  };

  const int nt = ev.time.size();
  NumericVector ipred(nt);
  std::vector<double> a(n, 0.0), a2(n, 0.0);
  for (int i = 0; i < nt; ++i) {
    if (i > 0) {
      double dt = ev.time[i] - ev.time[i - 1];
      int steps = std::max(1, (int)std::ceil(dt * 20.0));
      double h = dt / steps;
      for (int s = 0; s < steps; ++s) {
        if (route == ROUTE_IV_INFUSION && effective_rate(ev, i) > 0.0) {
          a[0] += h * effective_rate(ev, i);
        }
        rk4_step(kmat, n, h, a, a2);
        a.swap(a2);
      }
    }
    if (apply_ss_at_dose(model_ss, row_ss(ev, i), advan, route, p, ev, i, a)) {
      // SS state set; bolus/oral dose already included where applicable
    } else if (is_reset_row(ev, i)) {
      std::fill(a.begin(), a.end(), 0.0);
    } else if (is_dosing_row(ev, i)) {
      const int cmt = row_cmt(ev, i, pk_route_dose_cmp);
      const double dose = ev.amt[i] * row_f1(ev, i, p);
      if (dose > 0.0) {
        const int idx = cmt - 1;
        if (idx >= 0 && idx < static_cast<int>(a.size())) {
          a[static_cast<size_t>(idx)] += dose;
        } else if (!a.empty()) {
          a[0] += dose;
        }
      }
    }
    if (route == ROUTE_ORAL) {
      ipred[i] = conc_oral(event_obs_cmp(obs_cmp, ev, i), a[0], n > 1 ? a[1] : 0.0, n > 2 ? a[2] : 0.0, n > 3 ? a[3] : 0.0, p, ev, i);
    } else {
      ipred[i] = conc_iv(event_obs_cmp(obs_cmp, ev, i), a[0], n > 1 ? a[1] : 0.0, n > 2 ? a[2] : 0.0, p, ev, i);
    }
    pk_write_amounts(i, a);
  }
  return ipred;
}

}  // namespace

NumericVector nm_pk_route_cpp(
    int advan, int trans, int obs_cmp, int dose_cmp,
    int n_transit, bool use_ode, int model_ss,
    SubjectEvents ev, PkParams params) {
  pk_route_dose_cmp = dose_cmp > 0 ? dose_cmp : 1;
  RouteType route = ROUTE_IV_BOLUS;
  int cmt0 = first_dosing_cmt(ev, dose_cmp);
  double rate0 = 0.0;
  for (int i = 0; i < ev.time.size(); ++i) {
    if (is_dosing_row(ev, i)) {
      if (ev.rate.size() > 0) rate0 = ev.rate[i];
      break;
    }
  }
  route = route_type(advan, trans, cmt0, rate0, n_transit);
  apply_trans(trans, advan, route, params);
  if (params.kmf > 0.0 || params.kme > 0.0) params.has_metabolite = true;

  if (advan == 10) {
    return solve_mm_iv(ev, params, obs_cmp, model_ss);
  }

  const int ncomp = pk_ncomp(advan);

  if (params.has_metabolite && !use_ode && advan < 6) {
    return solve_metab_route(ev, params, advan, route, obs_cmp, model_ss);
  }

  if (n_transit > 0 && !use_ode) {
    if (ncomp == 1) return solve_1comp_transit(ev, params, n_transit, obs_cmp, model_ss);
    if (ncomp == 2) return solve_2comp_transit(ev, params, n_transit, obs_cmp, model_ss);
    return solve_ode(advan, ROUTE_TRANSIT, obs_cmp, ev, params, model_ss);
  }

  if (use_ode || advan >= 6) {
    return solve_ode(advan, route, obs_cmp, ev, params, model_ss);
  }

  if (ncomp == 1) {
    if (route == ROUTE_IV_INFUSION) return solve_1_iv_inf(ev, params, obs_cmp, model_ss, advan);
    if (route == ROUTE_ORAL) return solve_1_oral(ev, params, obs_cmp, model_ss, advan);
    return solve_1_iv_bolus(ev, params, obs_cmp, model_ss, advan);
  }

  if (ncomp == 2) {
    if (route == ROUTE_ORAL) return solve_2_oral(ev, params, obs_cmp, model_ss, advan);
    if (route == ROUTE_IV_INFUSION) return solve_2_iv_inf(ev, params, obs_cmp, model_ss, advan);
    return solve_2_iv_bolus(ev, params, obs_cmp, model_ss, advan);
  }

  if (ncomp == 3) {
    if (route == ROUTE_ORAL || advan == 12) return solve_3_oral(ev, params, obs_cmp, model_ss, advan, trans);
    if (route == ROUTE_IV_INFUSION) return solve_3_iv_inf(ev, params, obs_cmp, model_ss, advan);
    return solve_3_iv_bolus(ev, params, obs_cmp, model_ss, advan);
  }

  return solve_ode(advan, route, obs_cmp, ev, params, model_ss);
}

List nm_pk_route_detail_cpp(
    int advan, int trans, int obs_cmp, int dose_cmp,
    int n_transit, bool use_ode, int model_ss,
    SubjectEvents ev, PkParams params, int n_state) {
  const int nt = ev.time.size();
  const int ncol = std::max(1, std::min(n_state, kMaxPkState));
  NumericMatrix amounts(nt, ncol);
  std::fill(amounts.begin(), amounts.end(), 0.0);
  CharacterVector cnames(ncol);
  for (int j = 0; j < ncol; ++j) {
    cnames[j] = "A" + std::to_string(j + 1);
  }
  colnames(amounts) = cnames;
  pk_route_amounts_ptr = &amounts;
  NumericVector ipred = nm_pk_route_cpp(
      advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss, ev, params);
  pk_route_amounts_ptr = nullptr;
  return List::create(
      Named("ipred") = ipred,
      Named("amounts") = amounts);
}
