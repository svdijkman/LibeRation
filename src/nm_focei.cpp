#include <Rcpp.h>
#include <cmath>
#include <limits>
#include <vector>
#include "nm_subject_nll.h"
#include "nm_lik_config.h"
using namespace Rcpp;
using nm_subject::subject_f_at_obs;
using nm_subject::residual_var;
using nm_subject::extras_from_list;

namespace {

double shi_fd_eps(double x, double f_ref = 0.0) {
  const double tau = 1e-3;
  const double base = std::sqrt(std::numeric_limits<double>::epsilon());
  double h = base * std::max(std::abs(x), tau);
  const double fs = std::abs(f_ref);
  if (fs > 0.0 && std::isfinite(fs)) {
    h = std::max(h, base * std::max(fs, tau));
  }
  if (h < 1e-7) h = 1e-7;
  if (h > 0.05) h = 0.05;
  return h;
}

double log_det_spd(NumericMatrix M) {
  const int p = M.nrow();
  if (p == 0) return 0.0;
  NumericMatrix L(p, p);
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j <= i; ++j) {
      double s = M(i, j);
      for (int k = 0; k < j; ++k) {
        s -= L(i, k) * L(j, k);
      }
      if (i == j) {
        if (s <= 1e-15) return R_NegInf;
        L(i, j) = std::sqrt(s);
      } else {
        L(i, j) = s / L(j, j);
      }
    }
  }
  double logdet = 0.0;
  for (int i = 0; i < p; ++i) {
    logdet += std::log(L(i, i));
  }
  return 2.0 * logdet;
}

NumericMatrix omega_matrix(const NumericVector& omega, int n_eta) {
  NumericMatrix OM(n_eta, n_eta);
  for (int i = 0; i < n_eta; ++i) {
    for (int j = 0; j < n_eta; ++j) {
      OM(i, j) = 0.0;
    }
  }
  if (n_eta == 0) return OM;
  if (nm_lik::config().omega_type == nm_lik::OMEGA_BLOCK2 && n_eta >= 2 &&
      omega.size() >= 3) {
    OM(0, 0) = std::max(omega[0], 1e-15);
    OM(1, 1) = std::max(omega[1], 1e-15);
    OM(0, 1) = OM(1, 0) = omega[2];
    for (int i = 2; i < n_eta; ++i) {
      const double v = (i < omega.size()) ? omega[i] : omega[0];
      OM(i, i) = std::max(v, 1e-15);
    }
    return OM;
  }
  for (int i = 0; i < n_eta; ++i) {
    const double v = (i < omega.size()) ? omega[i] : omega[0];
    OM(i, i) = std::max(v, 1e-15);
  }
  return OM;
}

NumericMatrix invert_spd(NumericMatrix OM) {
  const int p = OM.nrow();
  NumericMatrix inv(p, p);
  if (p == 1) {
    inv(0, 0) = 1.0 / OM(0, 0);
    return inv;
  }
  if (p == 2) {
    const double a = OM(0, 0);
    const double b = OM(0, 1);
    const double c = OM(1, 1);
    const double det = std::max(a * c - b * b, 1e-15);
    inv(0, 0) = c / det;
    inv(0, 1) = inv(1, 0) = -b / det;
    inv(1, 1) = a / det;
    return inv;
  }
  Environment base("package:base");
  Function solve = base["solve"];
  return NumericMatrix(solve(OM));
}

void subject_f_and_G(
    const List& subj,
    const NumericVector& eta,
    const NumericVector& theta,
    const NumericVector& omega,
    const NumericVector& sigma,
    const CharacterVector& pred_lines,
    int advan, int trans, int obs_cmp, int dose_cmp,
    int n_transit, bool use_ode, int model_ss,
    NumericVector& f_out,
    NumericVector& dv_out,
    NumericMatrix& G_out) {
  const int p = eta.size();
  List cov = subj.containsElementNamed("cov") ? subj["cov"] : List();
  CharacterVector des = CharacterVector();
  auto ex = nm_subject::extras_from_list(subj);
  subject_f_at_obs(
      subj["time"], subj["amt"], subj["rate"], subj["f1"], subj["cmt"],
      subj["evid"], subj["ss"], subj["ii"], subj["dv"], subj["obs_idx"],
      eta, theta, pred_lines, advan, trans, obs_cmp, dose_cmp,
      n_transit, use_ode, model_ss, f_out, dv_out, cov, des, ex);
  const int n = f_out.size();
  G_out = NumericMatrix(n, p);
  if (p == 0 || n == 0) return;
  double f_scale = 0.0;
  for (int j = 0; j < n; ++j) {
    f_scale = std::max(f_scale, std::abs(f_out[j]));
  }
  for (int k = 0; k < p; ++k) {
    const double eps = shi_fd_eps(eta[k], f_scale);
    NumericVector etap = clone(eta);
    NumericVector etam = clone(eta);
    etap[k] += eps;
    etam[k] -= eps;
    NumericVector fp, fm, dvp, dvm;
    subject_f_at_obs(
        subj["time"], subj["amt"], subj["rate"], subj["f1"], subj["cmt"],
        subj["evid"], subj["ss"], subj["ii"], subj["dv"], subj["obs_idx"],
        etap, theta, pred_lines, advan, trans, obs_cmp, dose_cmp,
        n_transit, use_ode, model_ss, fp, dvp, cov, des, ex);
    subject_f_at_obs(
        subj["time"], subj["amt"], subj["rate"], subj["f1"], subj["cmt"],
        subj["evid"], subj["ss"], subj["ii"], subj["dv"], subj["obs_idx"],
        etam, theta, pred_lines, advan, trans, obs_cmp, dose_cmp,
        n_transit, use_ode, model_ss, fm, dvm, cov, des, ex);
    for (int j = 0; j < n; ++j) {
      G_out(j, k) = (fp[j] - fm[j]) / (2.0 * eps);
    }
  }
}

double omega_prior_obj(const NumericVector& eta, const NumericMatrix& OM) {
  const int p = eta.size();
  if (p == 0) return 0.0;
  NumericMatrix invOM = invert_spd(OM);
  double quad = 0.0;
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < p; ++j) {
      quad += eta[i] * invOM(i, j) * eta[j];
    }
  }
  return quad + log_det_spd(OM);
}

double subject_focei_objective(
    const List& subj,
    const NumericVector& eta,
    const NumericVector& theta,
    const NumericVector& omega,
    const NumericVector& sigma,
    const CharacterVector& pred_lines,
    int advan, int trans, int obs_cmp, int dose_cmp,
    int n_transit, bool use_ode, int model_ss,
    double* interaction_out = nullptr) {
  const int p = eta.size();
  if (p == 0) {
    if (interaction_out) *interaction_out = 0.0;
    return 0.0;
  }
  NumericVector f, dv;
  NumericMatrix G;
  subject_f_and_G(
      subj, eta, theta, omega, sigma, pred_lines,
      advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
      f, dv, G);
  const int n = f.size();
  if (n == 0) {
    if (interaction_out) *interaction_out = 0.0;
    return 0.0;
  }
  const int err = nm_lik::config().error_type;
  const double s1 = sigma.size() > 0 ? sigma[0] : 0.0;
  const double s2 = sigma.size() > 1 ? sigma[1] : 0.0;
  NumericMatrix OM = omega_matrix(omega, p);
  NumericMatrix invOM = invert_spd(OM);
  double resid_obj = 0.0;
  NumericMatrix gvg(p, p);
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < p; ++j) gvg(i, j) = 0.0;
  }
  for (int j = 0; j < n; ++j) {
    const double yhat = f[j];
    const double rj = residual_var(yhat, s1, s2, err);
    const double inv_r = 1.0 / rj;
    const double resid = dv[j] - yhat;
    resid_obj += std::log(rj) + resid * resid * inv_r;
    for (int i = 0; i < p; ++i) {
      for (int k = 0; k < p; ++k) {
        gvg(i, k) += inv_r * G(j, i) * G(j, k);
      }
    }
  }
  NumericMatrix M(p, p);
  for (int i = 0; i < p; ++i) {
    for (int k = 0; k < p; ++k) {
      M(i, k) = invOM(i, k) + gvg(i, k);
    }
  }
  const double interaction = log_det_spd(M);
  if (interaction_out) *interaction_out = interaction;
  return resid_obj + omega_prior_obj(eta, OM) + interaction;
}

}  // namespace

// [[Rcpp::export]]
double nm_focei_objective_cpp(
    List subjects,
    NumericMatrix eta_modes,
    NumericVector theta,
    NumericVector omega,
    NumericVector sigma,
    CharacterVector pred_lines,
    int advan,
    int trans,
    int obs_cmp,
    int dose_cmp,
    int n_transit,
    bool use_ode,
    int model_ss) {
  const int n_sub = subjects.size();
  double total = 0.0;
  for (int s = 0; s < n_sub; ++s) {
    NumericVector eta_s(eta_modes.ncol());
    for (int j = 0; j < eta_modes.ncol(); ++j) {
      eta_s[j] = eta_modes(s, j);
    }
    total += subject_focei_objective(
        subjects[s], eta_s, theta, omega, sigma, pred_lines,
        advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss);
  }
  return total;
}

// [[Rcpp::export]]
double nm_focei_interaction_cpp(
    List subjects,
    NumericMatrix eta_modes,
    NumericVector theta,
    NumericVector omega,
    NumericVector sigma,
    CharacterVector pred_lines,
    int advan,
    int trans,
    int obs_cmp,
    int dose_cmp,
    int n_transit,
    bool use_ode,
    int model_ss) {
  const int n_sub = subjects.size();
  double penalty = 0.0;
  for (int s = 0; s < n_sub; ++s) {
    NumericVector eta_s(eta_modes.ncol());
    for (int j = 0; j < eta_modes.ncol(); ++j) {
      eta_s[j] = eta_modes(s, j);
    }
    if (eta_s.size() == 0) continue;
    double inter = 0.0;
    subject_focei_objective(
        subjects[s], eta_s, theta, omega, sigma, pred_lines,
        advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss, &inter);
    penalty += inter;
  }
  return penalty;
}
