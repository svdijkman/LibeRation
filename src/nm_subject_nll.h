#pragma once
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include "nm_pk_core.h"
#include "nm_lik_config.h"

namespace nm_subject {

struct SubjectEventExtras {
  Rcpp::NumericVector s1;
  Rcpp::NumericVector s2;
  Rcpp::NumericVector s3;
  Rcpp::NumericVector s4;
  Rcpp::NumericMatrix scale_mat;
  bool use_data_scale = false;
  Rcpp::NumericMatrix f_mat;
  bool use_data_f = false;
};

inline SubjectEventExtras extras_from_list(const Rcpp::List& subj) {
  SubjectEventExtras ex;
  if (subj.containsElementNamed("s1")) ex.s1 = subj["s1"];
  if (subj.containsElementNamed("s2")) ex.s2 = subj["s2"];
  if (subj.containsElementNamed("s3")) ex.s3 = subj["s3"];
  if (subj.containsElementNamed("s4")) ex.s4 = subj["s4"];
  if (subj.containsElementNamed("scale_mat")) {
    ex.scale_mat = Rcpp::as<Rcpp::NumericMatrix>(subj["scale_mat"]);
  }
  if (subj.containsElementNamed("use_data_scale")) {
    ex.use_data_scale = Rcpp::as<bool>(subj["use_data_scale"]);
  }
  if (subj.containsElementNamed("f_mat")) {
    ex.f_mat = Rcpp::as<Rcpp::NumericMatrix>(subj["f_mat"]);
  }
  if (subj.containsElementNamed("use_data_f")) {
    ex.use_data_f = Rcpp::as<bool>(subj["use_data_f"]);
  }
  return ex;
}

inline double residual_var(
    double yhat,
    double s1,
    double s2,
    int error_type) {
  switch (error_type) {
    case nm_lik::ERR_ADD:
      return std::max(s1 * s1, 1e-15);
    case nm_lik::ERR_PROP:
      return std::max((s1 * yhat) * (s1 * yhat), 1e-15);
    case nm_lik::ERR_LOG:
      return std::max(s1 * s1, 1e-15);
    case nm_lik::ERR_POWER: {
      const double fpos = std::max(std::abs(yhat), 1e-8);
      return std::max(s1 * s1 * std::pow(fpos, 2.0 * s2), 1e-15);
    }
    case nm_lik::ERR_PROPADD:
    default:
      return std::max((s1 * yhat) * (s1 * yhat) + s2 * s2, 1e-15);
  }
}

inline double residual_nll(
    const Rcpp::NumericVector& dv,
    const Rcpp::NumericVector& f,
    double s1,
    double s2,
    int error_type,
    const Rcpp::IntegerVector& dvid = Rcpp::IntegerVector(),
    const Rcpp::NumericVector& sigma = Rcpp::NumericVector()) {
  const int n = dv.size();
  if (n == 0) return 0.0;
  auto sigma_at = [&](int i) -> std::pair<double, double> {
    if (dvid.size() == n && sigma.size() >= 2) {
      const int k = dvid[i];
      if (k >= 1 && (k * 2) <= sigma.size()) {
        return {sigma[k * 2 - 2], sigma[k * 2 - 1]};
      }
    }
    return {s1, s2};
  };
  if (nm_lik::config().sigma_corr == nm_lik::SIGMA_AR1 && n > 1) {
    const double rho = nm_lik::config().ar1_rho;
    Rcpp::NumericVector r(n);
    Rcpp::NumericVector v(n);
    for (int i = 0; i < n; ++i) {
      const double yhat = f[i];
      if (error_type == nm_lik::ERR_LOG) {
        r[i] = std::log(std::max(dv[i], 1e-8)) - std::log(std::max(yhat, 1e-8));
      } else {
        r[i] = dv[i] - yhat;
      }
      const auto sg = sigma_at(i);
      v[i] = residual_var(yhat, sg.first, sg.second, error_type);
    }
    double nll = std::log(v[0]) + r[0] * r[0] / v[0];
    for (int i = 1; i < n; ++i) {
      const double var_eps = v[i] * (1.0 - rho * rho);
      const double e = r[i] - rho * r[i - 1];
      nll += std::log(std::max(var_eps, 1e-15)) + e * e / std::max(var_eps, 1e-15);
    }
    return nll;
  }
  double nll = 0.0;
  for (int i = 0; i < n; ++i) {
    const double yhat = f[i];
    double resid = dv[i] - yhat;
    if (error_type == nm_lik::ERR_LOG) {
      const double lf = std::log(std::max(yhat, 1e-8));
      const double ly = std::log(std::max(dv[i], 1e-8));
      resid = ly - lf;
    }
    const auto sg = sigma_at(i);
    const double var = residual_var(yhat, sg.first, sg.second, error_type);
    nll += std::log(var) + resid * resid / var;
  }
  return nll;
}

inline double prop_add_nll(
    const Rcpp::NumericVector& dv,
    const Rcpp::NumericVector& f,
    double s1,
    double s2) {
  return residual_nll(dv, f, s1, s2, nm_lik::ERR_PROPADD);
}

inline void prop_add_nll_grad_sigma(
    const Rcpp::NumericVector& dv,
    const Rcpp::NumericVector& f,
    double s1,
    double s2,
    double& g_s1,
    double& g_s2) {
  g_s1 = 0.0;
  g_s2 = 0.0;
  const int n = dv.size();
  for (int i = 0; i < n; ++i) {
    const double yhat = f[i];
    const double resid = dv[i] - yhat;
    const double var = std::max((s1 * yhat) * (s1 * yhat) + s2 * s2, 1e-15);
    const double inv_var = 1.0 / var;
    const double inv_var2 = inv_var * inv_var;
    const double dvar_ds1 = 2.0 * s1 * yhat * yhat;
    const double dvar_ds2 = 2.0 * s2;
    g_s1 += dvar_ds1 * inv_var - resid * resid * dvar_ds1 * inv_var2;
    g_s2 += dvar_ds2 * inv_var - resid * resid * dvar_ds2 * inv_var2;
  }
}

inline double logsumexp(const Rcpp::NumericVector& x) {
  const int n = x.size();
  if (n == 0) return R_NegInf;
  double m = x[0];
  for (int i = 1; i < n; ++i) {
    if (x[i] > m) m = x[i];
  }
  double s = 0.0;
  for (int i = 0; i < n; ++i) {
    s += std::exp(x[i] - m);
  }
  return m + std::log(s);
}

inline double omega_prior_nll_diag(
    const Rcpp::NumericVector& eta,
    const Rcpp::NumericVector& omega) {
  double nll = 0.0;
  for (int i = 0; i < eta.size(); ++i) {
    const double om = std::max(omega[i], 1e-15);
    nll += eta[i] * eta[i] / om + std::log(om);
  }
  return nll;
}

// Block 2x2 for first two etas (omega[0]=var1, omega[1]=var2, omega[2]=cov12)
inline double omega_prior_nll_block2(
    const Rcpp::NumericVector& eta,
    const Rcpp::NumericVector& omega) {
  if (eta.size() < 2 || omega.size() < 3) {
    return omega_prior_nll_diag(eta, omega);
  }
  const double v1 = std::max(omega[0], 1e-15);
  const double v2 = std::max(omega[1], 1e-15);
  const double c12 = omega[2];
  const double det = std::max(v1 * v2 - c12 * c12, 1e-15);
  const double e1 = eta[0];
  const double e2 = eta[1];
  const double quad = (v2 * e1 * e1 - 2.0 * c12 * e1 * e2 + v1 * e2 * e2) / det;
  double nll = std::log(det) + quad;
  if (eta.size() > 2) {
    Rcpp::NumericVector eta_rest(eta.size() - 2);
    Rcpp::NumericVector om_rest(eta.size() - 2);
    for (int i = 2; i < eta.size(); ++i) {
      eta_rest[i - 2] = eta[i];
      om_rest[i - 2] = (omega.size() > i) ? omega[i] : v1;
    }
    nll += omega_prior_nll_diag(eta_rest, om_rest);
  }
  return nll;
}

inline double omega_prior_nll(
    const Rcpp::NumericVector& eta,
    const Rcpp::NumericVector& omega) {
  const int iov = nm_lik::config().iov;
  if (iov > 0 && eta.size() > iov) {
    const int n_between = eta.size() - iov;
    Rcpp::NumericVector eta_b(n_between);
    Rcpp::NumericVector om_b(n_between);
    Rcpp::NumericVector eta_i(iov);
    Rcpp::NumericVector om_i(iov);
    for (int i = 0; i < n_between; ++i) {
      eta_b[i] = eta[i];
      om_b[i] = (i < omega.size()) ? omega[i] : omega[0];
    }
    for (int i = 0; i < iov; ++i) {
      eta_i[i] = eta[n_between + i];
      const int oi = n_between + i;
      om_i[i] = (oi < omega.size()) ? omega[oi] : omega[0];
    }
    double nll = 0.0;
    if (nm_lik::config().omega_type == nm_lik::OMEGA_BLOCK2 && n_between >= 2) {
      nll += omega_prior_nll_block2(eta_b, om_b);
    } else {
      nll += omega_prior_nll_diag(eta_b, om_b);
    }
    nll += omega_prior_nll_diag(eta_i, om_i);
    return nll;
  }
  if (nm_lik::config().omega_type == nm_lik::OMEGA_BLOCK2) {
    return omega_prior_nll_block2(eta, omega);
  }
  return omega_prior_nll_diag(eta, omega);
}

inline double omega_prior_grad_j(double eta_j, double omega_j) {
  const double om = std::max(omega_j, 1e-15);
  return -eta_j * eta_j / (om * om) + 1.0 / om;
}

inline Rcpp::NumericVector subject_ipred(
    const Rcpp::NumericVector& time,
    const Rcpp::NumericVector& amt,
    const Rcpp::NumericVector& rate,
    const Rcpp::NumericVector& f1,
    const Rcpp::IntegerVector& cmt,
    const Rcpp::IntegerVector& evid,
    const Rcpp::IntegerVector& ss,
    const Rcpp::NumericVector& ii,
    const Rcpp::CharacterVector& pred_lines,
    const Rcpp::NumericVector& eta,
    const Rcpp::NumericVector& theta,
    int advan,
    int trans,
    int obs_cmp,
    int dose_cmp,
    int n_transit,
    bool use_ode,
    int model_ss,
    const Rcpp::List& cov = Rcpp::List(),
    const Rcpp::CharacterVector& des_lines = Rcpp::CharacterVector(),
    const SubjectEventExtras& ex = SubjectEventExtras()) {
  Rcpp::List pk = nm_eval_pred_cpp(pred_lines, theta, eta, cov, des_lines);
  return nm_pk_route_r(
      advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
      time, amt, rate, f1, cmt, evid, ss, ii, pk,
      ex.s1, ex.s2, ex.s3, ex.s4, ex.scale_mat, ex.use_data_scale, ex.f_mat, ex.use_data_f);
}

inline void subject_f_at_obs(
    const Rcpp::NumericVector& time,
    const Rcpp::NumericVector& amt,
    const Rcpp::NumericVector& rate,
    const Rcpp::NumericVector& f1,
    const Rcpp::IntegerVector& cmt,
    const Rcpp::IntegerVector& evid,
    const Rcpp::IntegerVector& ss,
    const Rcpp::NumericVector& ii,
    const Rcpp::NumericVector& dv,
    const Rcpp::IntegerVector& obs_idx,
    const Rcpp::NumericVector& eta,
    const Rcpp::NumericVector& theta,
    const Rcpp::CharacterVector& pred_lines,
    int advan,
    int trans,
    int obs_cmp,
    int dose_cmp,
    int n_transit,
    bool use_ode,
    int model_ss,
    Rcpp::NumericVector& f_out,
    Rcpp::NumericVector& dv_out,
    const Rcpp::List& cov = Rcpp::List(),
    const Rcpp::CharacterVector& des_lines = Rcpp::CharacterVector(),
    const SubjectEventExtras& ex = SubjectEventExtras()) {
  Rcpp::NumericVector ipred = subject_ipred(
      time, amt, rate, f1, cmt, evid, ss, ii, pred_lines, eta, theta,
      advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss, cov, des_lines, ex);
  const int n = obs_idx.size();
  f_out = Rcpp::NumericVector(n);
  dv_out = Rcpp::NumericVector(n);
  for (int i = 0; i < n; ++i) {
    const int idx = obs_idx[i] - 1;
    f_out[i] = ipred[idx];
    dv_out[i] = dv[i];
  }
}

inline double subject_nll_from_f(
    const Rcpp::NumericVector& dv_obs,
    const Rcpp::NumericVector& f,
    const Rcpp::NumericVector& eta,
    const Rcpp::NumericVector& omega,
    const Rcpp::NumericVector& sigma,
    bool include_omega_prior,
    const Rcpp::IntegerVector& dvid = Rcpp::IntegerVector()) {
  const int err = nm_lik::config().error_type;
  double nll = residual_nll(
      dv_obs, f,
      sigma.size() > 0 ? sigma[0] : 0.0,
      sigma.size() > 1 ? sigma[1] : 0.0,
      err, dvid, sigma);
  if (include_omega_prior && eta.size() > 0) {
    nll += omega_prior_nll(eta, omega);
  }
  return nll;
}

inline double subject_nll(
    const Rcpp::NumericVector& time,
    const Rcpp::NumericVector& amt,
    const Rcpp::NumericVector& rate,
    const Rcpp::NumericVector& f1,
    const Rcpp::IntegerVector& cmt,
    const Rcpp::IntegerVector& evid,
    const Rcpp::IntegerVector& ss,
    const Rcpp::NumericVector& ii,
    const Rcpp::NumericVector& dv,
    const Rcpp::IntegerVector& obs_idx,
    const Rcpp::NumericVector& eta,
    const Rcpp::NumericVector& theta,
    const Rcpp::NumericVector& omega,
    const Rcpp::NumericVector& sigma,
    const Rcpp::CharacterVector& pred_lines,
    int advan,
    int trans,
    int obs_cmp,
    int dose_cmp,
    int n_transit,
    bool use_ode,
    int model_ss,
    bool include_omega_prior = true,
    const Rcpp::List& cov = Rcpp::List(),
    const Rcpp::IntegerVector& dvid = Rcpp::IntegerVector(),
    const Rcpp::CharacterVector& des_lines = Rcpp::CharacterVector(),
    const SubjectEventExtras& ex = SubjectEventExtras()) {
  Rcpp::NumericVector f;
  Rcpp::NumericVector dv_obs;
  subject_f_at_obs(
      time, amt, rate, f1, cmt, evid, ss, ii, dv, obs_idx, eta, theta,
      pred_lines, advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
      f, dv_obs, cov, des_lines, ex);
  return subject_nll_from_f(dv_obs, f, eta, omega, sigma, include_omega_prior, dvid);
}

inline void subject_nll_grad(
    const Rcpp::NumericVector& time,
    const Rcpp::NumericVector& amt,
    const Rcpp::NumericVector& rate,
    const Rcpp::NumericVector& f1,
    const Rcpp::IntegerVector& cmt,
    const Rcpp::IntegerVector& evid,
    const Rcpp::IntegerVector& ss,
    const Rcpp::NumericVector& ii,
    const Rcpp::NumericVector& dv,
    const Rcpp::IntegerVector& obs_idx,
    const Rcpp::NumericVector& eta,
    const Rcpp::NumericVector& theta,
    const Rcpp::NumericVector& omega,
    const Rcpp::NumericVector& sigma,
    const Rcpp::CharacterVector& pred_lines,
    int advan,
    int trans,
    int obs_cmp,
    int dose_cmp,
    int n_transit,
    bool use_ode,
    int model_ss,
    bool include_omega_prior,
    Rcpp::NumericVector& g_theta,
    Rcpp::NumericVector& g_omega,
    Rcpp::NumericVector& g_sigma,
    Rcpp::NumericVector& f_out,
    Rcpp::NumericVector& dv_out,
    const SubjectEventExtras& ex = SubjectEventExtras()) {
  subject_f_at_obs(
      time, amt, rate, f1, cmt, evid, ss, ii, dv, obs_idx, eta, theta,
      pred_lines, advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
      f_out, dv_out, Rcpp::List(), Rcpp::CharacterVector(), ex);
  const double s1 = sigma.size() > 0 ? sigma[0] : 0.0;
  const double s2 = sigma.size() > 1 ? sigma[1] : 0.0;
  const double eps = 1e-5;

  if (sigma.size() >= 1) {
    double gs1 = 0.0;
    double gs2 = 0.0;
    prop_add_nll_grad_sigma(dv_out, f_out, s1, s2, gs1, gs2);
    g_sigma[0] += gs1;
    if (sigma.size() >= 2) {
      g_sigma[1] += gs2;
    }
  }

  if (include_omega_prior && eta.size() > 0) {
    for (int j = 0; j < omega.size(); ++j) {
      g_omega[j] += omega_prior_grad_j(eta[j], omega[j]);
    }
  }

  const double prior_q = (include_omega_prior && eta.size() > 0)
      ? omega_prior_nll(eta, omega)
      : 0.0;
  for (int k = 0; k < theta.size(); ++k) {
    Rcpp::NumericVector thp = Rcpp::clone(theta);
    Rcpp::NumericVector thm = Rcpp::clone(theta);
    thp[k] += eps;
    thm[k] -= eps;
    Rcpp::NumericVector fp;
    Rcpp::NumericVector fm;
    Rcpp::NumericVector dv_obs = dv_out;
    subject_f_at_obs(
        time, amt, rate, f1, cmt, evid, ss, ii, dv, obs_idx, eta, thp,
        pred_lines, advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
        fp, dv_obs, Rcpp::List(), Rcpp::CharacterVector(), ex);
    subject_f_at_obs(
        time, amt, rate, f1, cmt, evid, ss, ii, dv, obs_idx, eta, thm,
        pred_lines, advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
        fm, dv_obs, Rcpp::List(), Rcpp::CharacterVector(), ex);
    const double nll_p = subject_nll_from_f(dv_obs, fp, eta, omega, sigma, false) + prior_q;
    const double nll_m = subject_nll_from_f(dv_obs, fm, eta, omega, sigma, false) + prior_q;
    g_theta[k] += (nll_p - nll_m) / (2.0 * eps);
  }
}

}  // namespace nm_subject
