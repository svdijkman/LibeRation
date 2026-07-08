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
    // Packing (matches .nm_focei_omega_matrix in R): omega[0..2] are the 2x2
    // block (var1, var2, cov12); the diagonal remainder for eta3..etaN starts
    // at omega[3], i.e. OM(i,i) = omega[i + 1] for i >= 2. The previous code
    // read omega[i] here, which for i == 2 is the covariance element, giving a
    // wrong Omega whenever n_eta > 2.
    for (int i = 2; i < n_eta; ++i) {
      const double v = (i + 1 < omega.size()) ? omega[i + 1] : omega[0];
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

// d R_j / d f_j divided by R_j, for the FOCE-INTER curvature term.
// R depends on eta only through f_j, so dR/deta = (dR/df) * G. add/log have
// no f-dependence (returns 0), leaving FOCEI identical to FOCE for those.
double residual_var_dlog_df(double yhat, double s1, double s2, double rj, int err) {
  switch (err) {
    case nm_lik::ERR_PROP:     // R = (s1 f)^2  -> R'/R = 2/f
      return (std::abs(yhat) > 1e-12) ? (2.0 / yhat) : 0.0;
    case nm_lik::ERR_PROPADD:  // R = (s1 f)^2 + s2^2 -> R' = 2 s1^2 f
      return (rj > 1e-300) ? (2.0 * s1 * s1 * yhat / rj) : 0.0;
    case nm_lik::ERR_POWER:    // R = s1^2 |f|^(2 s2) -> R'/R = 2 s2 / f
      return (std::abs(yhat) > 1e-12) ? (2.0 * s2 / yhat) : 0.0;
    case nm_lik::ERR_ADD:
    case nm_lik::ERR_LOG:
    default:
      return 0.0;
  }
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
    bool interaction = true,
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
    // FOCE curvature coefficient 1/R_j; FOCE-INTER adds 0.5 (R'_j/R_j)^2 to
    // account for the dependence of the residual variance on eta (the standard
    // NONMEM interaction term, cf. Wang 2007).
    double coef = inv_r;
    if (interaction) {
      const double rp_over_r = residual_var_dlog_df(yhat, s1, s2, rj, err);
      coef += 0.5 * rp_over_r * rp_over_r;
    }
    for (int i = 0; i < p; ++i) {
      for (int k = 0; k < p; ++k) {
        gvg(i, k) += coef * G(j, i) * G(j, k);
      }
    }
  }
  NumericMatrix M(p, p);
  for (int i = 0; i < p; ++i) {
    for (int k = 0; k < p; ++k) {
      M(i, k) = invOM(i, k) + gvg(i, k);
    }
  }
  const double curvature = log_det_spd(M);
  if (interaction_out) *interaction_out = curvature;
  return resid_obj + omega_prior_obj(eta, OM) + curvature;
}

// True NONMEM FO marginal -2LL for one subject.
//
// With f = f(eta = 0), G = d f / d eta at eta = 0, R = diag residual variance
// evaluated at eta = 0 predictions, and V = R + G Omega G', the per-subject
// contribution is  log|V| + r' V^{-1} r  with r = y - f. No etas are fitted.
double subject_fo_marginal(
    const List& subj,
    const NumericVector& theta,
    const NumericVector& omega,
    const NumericVector& sigma,
    int n_eta,
    const CharacterVector& pred_lines,
    int advan, int trans, int obs_cmp, int dose_cmp,
    int n_transit, bool use_ode, int model_ss) {
  NumericVector eta0(n_eta > 0 ? n_eta : 0);
  for (int k = 0; k < n_eta; ++k) eta0[k] = 0.0;
  NumericVector f, dv;
  NumericMatrix G;
  subject_f_and_G(
      subj, eta0, theta, omega, sigma, pred_lines,
      advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
      f, dv, G);
  const int n = f.size();
  if (n == 0) return 0.0;
  const int err = nm_lik::config().error_type;
  const double s1 = sigma.size() > 0 ? sigma[0] : 0.0;
  const double s2 = sigma.size() > 1 ? sigma[1] : 0.0;
  NumericMatrix OM = (n_eta > 0) ? omega_matrix(omega, n_eta) : NumericMatrix(0, 0);
  // V = diag(R) + G OM G'
  NumericMatrix V(n, n);
  NumericMatrix GO(n, n_eta > 0 ? n_eta : 0);
  for (int i = 0; i < n; ++i) {
    for (int a = 0; a < n_eta; ++a) {
      double acc = 0.0;
      for (int b = 0; b < n_eta; ++b) acc += G(i, b) * OM(b, a);
      GO(i, a) = acc;
    }
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double acc = 0.0;
      for (int a = 0; a < n_eta; ++a) acc += GO(i, a) * G(j, a);
      V(i, j) = acc;
    }
    V(i, i) += residual_var(f[i], s1, s2, err);
  }
  // Cholesky with a light jitter retry for numerical safety.
  NumericMatrix L(n, n);
  double jitter = 0.0;
  for (int attempt = 0; attempt < 3; ++attempt) {
    bool ok = true;
    for (int i = 0; i < n && ok; ++i) {
      for (int j = 0; j <= i; ++j) {
        double s = V(i, j) + (i == j ? jitter : 0.0);
        for (int k = 0; k < j; ++k) s -= L(i, k) * L(j, k);
        if (i == j) {
          if (s <= 0.0) { ok = false; break; }
          L(i, j) = std::sqrt(s);
        } else {
          L(i, j) = s / L(j, j);
        }
      }
    }
    if (ok) break;
    double scale = 0.0;
    for (int i = 0; i < n; ++i) scale = std::max(scale, std::abs(V(i, i)));
    jitter = (jitter > 0.0 ? jitter * 10.0 : 1e-10 * std::max(scale, 1.0));
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) L(i, j) = 0.0;
    if (attempt == 2) return R_PosInf;
  }
  double logdet = 0.0;
  for (int i = 0; i < n; ++i) logdet += std::log(L(i, i));
  logdet *= 2.0;
  // Solve V x = r via L L' x = r (r = dv - f).
  std::vector<double> w(n, 0.0), x(n, 0.0);
  for (int i = 0; i < n; ++i) {
    double s = dv[i] - f[i];
    for (int k = 0; k < i; ++k) s -= L(i, k) * w[k];
    w[i] = s / L(i, i);
  }
  for (int i = n - 1; i >= 0; --i) {
    double s = w[i];
    for (int k = i + 1; k < n; ++k) s -= L(k, i) * x[k];
    x[i] = s / L(i, i);
  }
  double quad = 0.0;
  for (int i = 0; i < n; ++i) quad += (dv[i] - f[i]) * x[i];
  return logdet + quad;
}

}  // namespace

// [[Rcpp::export]]
double nm_fo_marginal_cpp(
    List subjects,
    NumericVector theta,
    NumericVector omega,
    NumericVector sigma,
    int n_eta,
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
    total += subject_fo_marginal(
        subjects[s], theta, omega, sigma, n_eta, pred_lines,
        advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss);
  }
  return total;
}

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
    int model_ss,
    bool interaction = true) {
  const int n_sub = subjects.size();
  double total = 0.0;
  for (int s = 0; s < n_sub; ++s) {
    NumericVector eta_s(eta_modes.ncol());
    for (int j = 0; j < eta_modes.ncol(); ++j) {
      eta_s[j] = eta_modes(s, j);
    }
    total += subject_focei_objective(
        subjects[s], eta_s, theta, omega, sigma, pred_lines,
        advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
        interaction);
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
        advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
        true, &inter);
    penalty += inter;
  }
  return penalty;
}
