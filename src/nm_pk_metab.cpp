#include <Rcpp.h>
#include <cmath>
#include "nm_pk_core.h"
#include "nm_pk_ss.h"
#include "nm_pk_pkadvan.h"
using namespace Rcpp;
using namespace nm_pk;
using namespace pkadvan;

namespace {

inline double ke_oral(const PkParams& p) {
  return p.k20 > 0 ? p.k20 : p.k10;
}

NumericVector solve_1_iv_metab(const SubjectEvents& ev, const PkParams& p, int obs_cmp, int model_ss) {
  const int n = ev.time.size();
  NumericVector ipred(n);
  double a1 = 0.0, am = 0.0;
  const double E1 = p.k10 + p.kmf;
  for (int i = 0; i < n; ++i) {
    if (i > 0) {
      double dt = ev.time[i] - ev.time[i - 1];
      double prev1 = a1, prevM = am;
      a1 = prev1 * std::exp(-dt * E1);
      am = prevM * std::exp(-dt * p.kme) +
          p.kmf * prev1 * (std::exp(-dt * p.kme) / (E1 - p.kme) +
                           std::exp(-dt * E1) / (p.kme - E1));
    }
    int evid = ev.evid.size() > 0 ? ev.evid[i] : (ev.amt[i] > 0 ? 1 : 0);
    if (evid == 1) {
      std::vector<double> st(1, a1);
      if (apply_ss_at_dose(model_ss, row_ss(ev, i), 1, ROUTE_IV_BOLUS, p, ev, i, st)) {
        a1 = st[0];
      } else {
        a1 += ev.amt[i];
      }
    }
    if (obs_cmp == 99 || obs_cmp == 5) ipred[i] = am / std::max(p.v1, 1e-8);
    else ipred[i] = conc_iv(event_obs_cmp(obs_cmp, ev, i), a1, 0.0, 0.0, p, ev, i);
    pk_write_amounts(i, std::vector<double>{a1, am});
  }
  return ipred;
}

NumericVector solve_2_iv_metab(const SubjectEvents& ev, const PkParams& p, int obs_cmp, int model_ss) {
  const int n = ev.time.size();
  NumericVector ipred(n);
  double a1 = 0.0, a2 = 0.0, am = 0.0;
  for (int i = 0; i < n; ++i) {
    if (i > 0) {
      double dt = ev.time[i] - ev.time[i - 1];
      Hybrid2 h = hybrid_lambdas_2(p.k10 + p.kmf, p.k12, p.k21);
      double d = h.lambda2 - h.lambda1;
      double prev1 = a1, prev2 = a2, prevM = am;
      if (std::fabs(d) < 1e-12) {
        a1 = (prev1 + prev2 * p.k21 / h.lambda1) * std::exp(-h.lambda1 * dt);
        a2 = prev2 * std::exp(-h.lambda1 * dt);
      } else {
        a1 = (((prev1 * h.E2 + prev2 * p.k21) - prev1 * h.lambda1) * std::exp(-h.lambda1 * dt) -
              ((prev1 * h.E2 + prev2 * p.k21) - prev1 * h.lambda2) * std::exp(-h.lambda2 * dt)) / d;
        a2 = (((prev2 * h.E1 + prev1 * p.k12) - prev2 * h.lambda1) * std::exp(-h.lambda1 * dt) -
              ((prev2 * h.E1 + prev1 * p.k12) - prev2 * h.lambda2) * std::exp(-h.lambda2 * dt)) / d;
      }
      am = prevM * std::exp(-dt * p.kme);
      am += p.kmf * (prev2 * p.k21 * (-std::exp(-dt * h.lambda2) / ((h.lambda1 - h.lambda2) * (h.lambda2 - p.kme)) +
                                      std::exp(-dt * h.lambda1) / ((h.lambda1 - h.lambda2) * (h.lambda1 - p.kme)) -
                                      std::exp(-dt * p.kme) / ((h.lambda1 - p.kme) * (p.kme - h.lambda2))) +
                       prev1 * (((h.lambda2 - h.E2) * std::exp(-dt * h.lambda2) / ((h.lambda1 - h.lambda2) * (h.lambda2 - p.kme))) +
                                ((p.kme - h.E2) * std::exp(-dt * p.kme) / ((h.lambda1 - p.kme) * (p.kme - h.lambda2))) +
                                ((h.E2 - h.lambda1) * std::exp(-dt * h.lambda1) / ((h.lambda1 - h.lambda2) * (h.lambda1 - p.kme)))));
    }
    int evid = ev.evid.size() > 0 ? ev.evid[i] : (ev.amt[i] > 0 ? 1 : 0);
    if (evid == 1) {
      int cmt = ev.cmt.size() > 0 ? ev.cmt[i] : 1;
      std::vector<double> st = {a1, a2};
      if (apply_ss_at_dose(model_ss, row_ss(ev, i), 3, ROUTE_IV_BOLUS, p, ev, i, st)) {
        a1 = st[0];
        a2 = st[1];
      } else {
        if (cmt == 1) a1 += ev.amt[i]; else if (cmt == 2) a2 += ev.amt[i]; else a1 += ev.amt[i];
      }
    }
    if (obs_cmp == 99 || obs_cmp == 5) ipred[i] = am / std::max(p.v1, 1e-8);
    else ipred[i] = conc_iv(event_obs_cmp(obs_cmp, ev, i), a1, a2, 0.0, p, ev, i);
    pk_write_amounts(i, std::vector<double>{a1, a2, am});
  }
  return ipred;
}

NumericVector solve_1_oral_metab(const SubjectEvents& ev, const PkParams& p, int obs_cmp, int model_ss) {
  const int n = ev.time.size();
  NumericVector ipred(n);
  double a_gut = 0.0, a_cen = 0.0, am = 0.0;
  const double E2 = ke_oral(p) + p.kmf;
  for (int i = 0; i < n; ++i) {
    if (i > 0) {
      double dt = ev.time[i] - ev.time[i - 1];
      double prevG = a_gut, prevC = a_cen, prevM = am;
      a_cen = (prevG * p.ka) / (p.ka - E2) * (std::exp(-E2 * dt) - std::exp(-p.ka * dt)) +
              prevC * std::exp(-E2 * dt);
      am = prevM * std::exp(-dt * p.kme);
      am += p.kmf * (prevG * p.ka * (-std::exp(-E2 * dt) / ((p.ka - E2) * (E2 - p.kme)) +
                                      std::exp(-p.ka * dt) / ((p.ka - E2) * (p.ka - p.kme)) -
                                      std::exp(-p.kme * dt) / ((p.ka - p.kme) * (p.kme - E2))) +
                       prevC * (std::exp(-E2 * dt) / (p.kme - E2) + std::exp(-p.kme * dt) / (E2 - p.kme)));
      a_gut = prevG * std::exp(-p.ka * dt);
    }
    int evid = ev.evid.size() > 0 ? ev.evid[i] : (ev.amt[i] > 0 ? 1 : 0);
    if (evid == 1) {
      const int cmt = row_cmt(ev, i, 1);
      a_gut += ev.amt[i] * effective_f(cmt, p, ev, i);
    }
    if (obs_cmp == 99 || obs_cmp == 5) ipred[i] = am / std::max(p.v1, 1e-8);
    else ipred[i] = conc_oral(event_obs_cmp(obs_cmp, ev, i), a_gut, a_cen, 0.0, 0.0, p, ev, i);
    pk_write_amounts(i, std::vector<double>{a_gut, a_cen, am});
  }
  return ipred;
}

}  // namespace

namespace nm_pk {

NumericVector solve_metab_route(const SubjectEvents& ev, const PkParams& p,
                                int advan, RouteType route, int obs_cmp, int model_ss) {
  const int ncomp = pk_ncomp(advan);
  if (ncomp == 1 && route == ROUTE_IV_BOLUS) return solve_1_iv_metab(ev, p, obs_cmp, model_ss);
  if (ncomp == 2 && route == ROUTE_IV_BOLUS) return solve_2_iv_metab(ev, p, obs_cmp, model_ss);
  if (ncomp == 1 && route == ROUTE_ORAL) return solve_1_oral_metab(ev, p, obs_cmp, model_ss);
  return solve_1_iv_metab(ev, p, obs_cmp, model_ss);
}

}  // namespace nm_pk
