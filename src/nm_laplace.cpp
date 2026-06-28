#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include "nm_subject_nll.h"
#include "nm_eta_mode.h"
using namespace Rcpp;
using nm_subject::logsumexp;
using nm_subject::subject_nll;
using nm_subject::subject_nll_from_f;
using nm_subject::subject_f_at_obs;
using nm_subject::prop_add_nll_grad_sigma;
using nm_subject::omega_prior_nll;
using nm_subject::omega_prior_grad_j;
using nm_eta::SubjectData;
using nm_eta::list_to_subject;
using nm_eta::find_eta_mode;

namespace {

struct GhGrid {
  NumericMatrix z;
  NumericVector w;
};

struct ObsCache {
  NumericVector f;
  NumericVector dv;
};

struct LaplaceSubjectResult {
  double val = R_PosInf;
  std::vector<double> log_terms;
  std::vector<double> nll_q;
  std::vector<std::vector<double> > eta_q;
  std::vector<double> alpha;
  std::vector<ObsCache> obs_cache;
  NumericVector eta_mode;
};

GhGrid build_gh_product_grid(const NumericVector& nodes, const NumericVector& weights, int n_eta) {
  GhGrid grid;
  if (n_eta <= 0) {
    grid.z = NumericMatrix(1, 0);
    grid.w = NumericVector::create(1.0);
    return grid;
  }
  const int n = nodes.size();
  int n_q = 1;
  for (int i = 0; i < n_eta; ++i) n_q *= n;
  grid.z = NumericMatrix(n_q, n_eta);
  grid.w = NumericVector(n_q);
  std::vector<int> idx(n_eta, 0);
  for (int q = 0; q < n_q; ++q) {
    double w = 1.0;
    for (int j = 0; j < n_eta; ++j) {
      grid.z(q, j) = nodes[idx[j]];
      w *= weights[idx[j]];
    }
    grid.w[q] = w;
    int carry = n_eta - 1;
    while (carry >= 0) {
      idx[carry]++;
      if (idx[carry] < n) break;
      idx[carry] = 0;
      carry--;
    }
  }
  return grid;
}

void eta_at_quadrature(
    NumericVector& eta,
    const NumericVector& mode,
    const NumericVector& omega,
    const NumericVector& z_q,
    bool mode_centered) {
  for (int j = 0; j < eta.size(); ++j) {
    if (mode_centered) {
      eta[j] = mode[j] + std::sqrt(std::max(omega[j], 1e-15)) * z_q[j];
    } else {
      eta[j] = z_q[j];
    }
  }
}

List fwd_to_list(const LaplaceSubjectResult& fwd) {
  const int n_q = static_cast<int>(fwd.nll_q.size());
  List obs_list(n_q);
  for (int q = 0; q < n_q; ++q) {
    obs_list[q] = List::create(
      Named("f") = fwd.obs_cache[q].f,
      Named("dv") = fwd.obs_cache[q].dv
    );
  }
  return List::create(
    Named("val") = fwd.val,
    Named("log_terms") = wrap(fwd.log_terms),
    Named("nll_q") = wrap(fwd.nll_q),
    Named("eta_q") = wrap(fwd.eta_q),
    Named("alpha") = wrap(fwd.alpha),
    Named("eta_mode") = fwd.eta_mode,
    Named("obs_cache") = obs_list
  );
}

LaplaceSubjectResult fwd_from_list(const List& x) {
  LaplaceSubjectResult fwd;
  fwd.val = as<double>(x["val"]);
  fwd.log_terms = as<std::vector<double> >(x["log_terms"]);
  fwd.nll_q = as<std::vector<double> >(x["nll_q"]);
  fwd.eta_q = as<std::vector<std::vector<double> > >(x["eta_q"]);
  fwd.alpha = as<std::vector<double> >(x["alpha"]);
  fwd.eta_mode = x["eta_mode"];
  const List obs_list = x["obs_cache"];
  const int n_q = obs_list.size();
  fwd.obs_cache.resize(n_q);
  for (int q = 0; q < n_q; ++q) {
    List oc = obs_list[q];
    fwd.obs_cache[q].f = oc["f"];
    fwd.obs_cache[q].dv = oc["dv"];
  }
  return fwd;
}

LaplaceSubjectResult laplace_subject_forward(
    const SubjectData& subj,
    const NumericVector& theta,
    const NumericVector& omega,
    const NumericVector& sigma,
    const GhGrid& grid,
    const CharacterVector& pred_lines,
    int advan,
    int trans,
    int obs_cmp,
    int dose_cmp,
    int n_transit,
    bool use_ode,
    int model_ss,
    bool mode_centered,
    SEXP eta_init = R_NilValue) {
  LaplaceSubjectResult out;
  const int n_eta = omega.size();
  const int n_q = grid.w.size();
  out.log_terms.resize(n_q);
  out.nll_q.resize(n_q);
  out.eta_q.resize(n_q, std::vector<double>(n_eta));
  out.obs_cache.resize(n_q);
  if (n_eta <= 0) {
    ObsCache oc;
    subject_f_at_obs(
        subj.time, subj.amt, subj.rate, subj.f1, subj.cmt, subj.evid, subj.ss, subj.ii,
        subj.dv, subj.obs_idx, NumericVector(0), theta, pred_lines,
        advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss, oc.f, oc.dv,
        List(), CharacterVector(), subj.extras);
    const double nll = subject_nll_from_f(oc.dv, oc.f, NumericVector(0), omega, sigma, false);
    out.val = nll;
    out.log_terms = std::vector<double>(1, -0.5 * nll);
    out.nll_q = std::vector<double>(1, nll);
    out.alpha = std::vector<double>(1, 1.0);
    out.obs_cache[0] = oc;
    return out;
  }

  out.eta_mode = find_eta_mode(
      subj, theta, omega, sigma, pred_lines,
      advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss, eta_init);

  NumericVector eta(n_eta);
  NumericVector z_q(n_eta);

  for (int q = 0; q < n_q; ++q) {
    for (int j = 0; j < n_eta; ++j) {
      z_q[j] = grid.z(q, j);
    }
    eta_at_quadrature(eta, out.eta_mode, omega, z_q, mode_centered);
    for (int j = 0; j < n_eta; ++j) out.eta_q[q][j] = eta[j];
    subject_f_at_obs(
        subj.time, subj.amt, subj.rate, subj.f1, subj.cmt, subj.evid, subj.ss, subj.ii,
        subj.dv, subj.obs_idx, eta, theta, pred_lines,
        advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
        out.obs_cache[q].f, out.obs_cache[q].dv,
        List(), CharacterVector(), subj.extras);
    const double nll_q = subject_nll_from_f(
        out.obs_cache[q].dv, out.obs_cache[q].f, eta, omega, sigma, true);
    out.nll_q[q] = nll_q;
    out.log_terms[q] = std::log(std::max(grid.w[q], 1e-300)) - 0.5 * nll_q;
  }

  const NumericVector log_terms_r = wrap(out.log_terms);
  const double log_sum = logsumexp(log_terms_r);
  out.val = std::isfinite(log_sum) ? -2.0 * log_sum : R_PosInf;
  out.alpha.resize(n_q);
  for (int q = 0; q < n_q; ++q) {
    out.alpha[q] = std::exp(out.log_terms[q] - log_sum);
  }
  return out;
}

void laplace_subject_grad_full(
    const SubjectData& subj,
    const GhGrid& grid,
    const NumericVector& theta,
    const NumericVector& omega,
    const NumericVector& sigma,
    const CharacterVector& pred_lines,
    int advan,
    int trans,
    int obs_cmp,
    int dose_cmp,
    int n_transit,
    bool use_ode,
    int model_ss,
    bool mode_centered,
    SEXP eta_init,
    NumericVector& g_theta,
    NumericVector& g_omega,
    NumericVector& g_sigma) {
  const double eps = 1e-5;
  const LaplaceSubjectResult base = laplace_subject_forward(
      subj, theta, omega, sigma, grid, pred_lines,
      advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
      mode_centered, eta_init);

  auto subject_val = [&](const NumericVector& th,
                         const NumericVector& om,
                         const NumericVector& sg,
                         SEXP warm) {
    const LaplaceSubjectResult fwd = laplace_subject_forward(
        subj, th, om, sg, grid, pred_lines,
        advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
        mode_centered, warm);
    return fwd.val;
  };

  std::fill(g_theta.begin(), g_theta.end(), 0.0);
  for (int k = 0; k < theta.size(); ++k) {
    NumericVector thp = clone(theta);
    NumericVector thm = clone(theta);
    thp[k] += eps;
    thm[k] -= eps;
    const double vp = subject_val(thp, omega, sigma, base.eta_mode);
    const double vm = subject_val(thm, omega, sigma, base.eta_mode);
    if (std::isfinite(vp) && std::isfinite(vm)) {
      g_theta[k] = (vp - vm) / (2.0 * eps);
    }
  }

  std::fill(g_omega.begin(), g_omega.end(), 0.0);
  for (int j = 0; j < omega.size(); ++j) {
    NumericVector omp = clone(omega);
    NumericVector omm = clone(omega);
    omp[j] += eps;
    omm[j] -= eps;
    const double vp = subject_val(theta, omp, sigma, base.eta_mode);
    const double vm = subject_val(theta, omm, sigma, base.eta_mode);
    if (std::isfinite(vp) && std::isfinite(vm)) {
      g_omega[j] = (vp - vm) / (2.0 * eps);
    }
  }

  std::fill(g_sigma.begin(), g_sigma.end(), 0.0);
  for (int j = 0; j < sigma.size(); ++j) {
    NumericVector sgp = clone(sigma);
    NumericVector sgm = clone(sigma);
    sgp[j] += eps;
    sgm[j] -= eps;
    const double vp = subject_val(theta, omega, sgp, base.eta_mode);
    const double vm = subject_val(theta, omega, sgm, base.eta_mode);
    if (std::isfinite(vp) && std::isfinite(vm)) {
      g_sigma[j] = (vp - vm) / (2.0 * eps);
    }
  }
}

void laplace_subject_grad(
    const SubjectData& subj,
    const LaplaceSubjectResult& fwd,
    const NumericVector& theta,
    const NumericVector& omega,
    const NumericVector& sigma,
    const CharacterVector& pred_lines,
    int advan,
    int trans,
    int obs_cmp,
    int dose_cmp,
    int n_transit,
    bool use_ode,
    int model_ss,
    NumericVector& g_theta,
    NumericVector& g_omega,
    NumericVector& g_sigma) {
  const int n_q = static_cast<int>(fwd.nll_q.size());
  const int n_eta = omega.size();
  const double eps = 1e-5;
  const double s1 = sigma.size() > 0 ? sigma[0] : 0.0;
  const double s2 = sigma.size() > 1 ? sigma[1] : 0.0;

  std::fill(g_omega.begin(), g_omega.end(), 0.0);
  for (int j = 0; j < omega.size(); ++j) {
    double acc = 0.0;
    for (int q = 0; q < n_q; ++q) {
      const double dnll = omega_prior_grad_j(fwd.eta_q[q][j], omega[j]);
      acc += fwd.alpha[q] * dnll;
    }
    g_omega[j] = acc;
  }

  std::fill(g_theta.begin(), g_theta.end(), 0.0);
  for (int k = 0; k < theta.size(); ++k) {
    double acc = 0.0;
    for (int q = 0; q < n_q; ++q) {
      NumericVector eta(n_eta);
      for (int j = 0; j < n_eta; ++j) eta[j] = fwd.eta_q[q][j];
      const double prior_q = (n_eta > 0)
          ? omega_prior_nll(eta, omega)
          : 0.0;
      NumericVector thp = clone(theta);
      NumericVector thm = clone(theta);
      thp[k] += eps;
      thm[k] -= eps;
      NumericVector fp;
      NumericVector fm;
      NumericVector dv_obs = fwd.obs_cache[q].dv;
      subject_f_at_obs(
          subj.time, subj.amt, subj.rate, subj.f1, subj.cmt, subj.evid, subj.ss, subj.ii,
          subj.dv, subj.obs_idx, eta, thp, pred_lines,
          advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss, fp, dv_obs,
          List(), CharacterVector(), subj.extras);
      subject_f_at_obs(
          subj.time, subj.amt, subj.rate, subj.f1, subj.cmt, subj.evid, subj.ss, subj.ii,
          subj.dv, subj.obs_idx, eta, thm, pred_lines,
          advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss, fm, dv_obs,
          List(), CharacterVector(), subj.extras);
      const double nll_p = subject_nll_from_f(dv_obs, fp, eta, omega, sigma, false) + prior_q;
      const double nll_m = subject_nll_from_f(dv_obs, fm, eta, omega, sigma, false) + prior_q;
      acc += fwd.alpha[q] * (nll_p - nll_m) / (2.0 * eps);
    }
    g_theta[k] = acc;
  }

  std::fill(g_sigma.begin(), g_sigma.end(), 0.0);
  if (sigma.size() >= 1) {
    double g_s1 = 0.0;
    double g_s2 = 0.0;
    for (int q = 0; q < n_q; ++q) {
      double gs1_q = 0.0;
      double gs2_q = 0.0;
      prop_add_nll_grad_sigma(
          fwd.obs_cache[q].dv, fwd.obs_cache[q].f, s1, s2, gs1_q, gs2_q);
      g_s1 += fwd.alpha[q] * gs1_q;
      g_s2 += fwd.alpha[q] * gs2_q;
    }
    g_sigma[0] = g_s1;
    if (sigma.size() >= 2) {
      g_sigma[1] = g_s2;
    }
  }
}

List laplace_engine(
    List subjects,
    NumericVector theta,
    NumericVector omega,
    NumericVector sigma,
    NumericVector gh_nodes,
    NumericVector gh_weights,
    CharacterVector pred_lines,
    int advan,
    int trans,
    int obs_cmp,
    int dose_cmp,
    int n_transit,
    bool use_ode,
    int model_ss,
    bool mode_centered,
    bool need_grad,
    SEXP eta_modes = R_NilValue,
    SEXP fwd_cache = R_NilValue,
    bool grad_from_fwd = false) {
  const int n_eta = omega.size();
  const GhGrid grid = build_gh_product_grid(gh_nodes, gh_weights, n_eta);
  const int n_sub = subjects.size();
  const bool use_fwd_cache =
      grad_from_fwd && fwd_cache != R_NilValue && !mode_centered;

  double total = 0.0;
  NumericVector g_theta = clone(theta);
  NumericVector g_omega = clone(omega);
  NumericVector g_sigma = clone(sigma);
  std::fill(g_theta.begin(), g_theta.end(), 0.0);
  std::fill(g_omega.begin(), g_omega.end(), 0.0);
  std::fill(g_sigma.begin(), g_sigma.end(), 0.0);

  List fwd_out(n_sub);
  NumericMatrix eta_modes_out;
  if (n_sub > 0 && n_eta > 0) {
    eta_modes_out = NumericMatrix(n_sub, n_eta);
  }

  List fwd_list;
  if (use_fwd_cache) {
    fwd_list = as<List>(fwd_cache);
    if (fwd_list.size() != n_sub) {
      stop("fwd_cache length mismatch");
    }
  }

  NumericMatrix eta_warm;
  if (eta_modes != R_NilValue) {
    eta_warm = as<NumericMatrix>(eta_modes);
    if (eta_warm.nrow() != n_sub || eta_warm.ncol() != n_eta) {
      stop("eta_modes dimension mismatch");
    }
  }

  for (int s = 0; s < n_sub; ++s) {
    SubjectData subj = list_to_subject(subjects[s]);
    LaplaceSubjectResult fwd;
    if (use_fwd_cache) {
      fwd = fwd_from_list(fwd_list[s]);
      total += fwd.val;
    } else {
      SEXP eta_init = R_NilValue;
      if (eta_modes != R_NilValue) {
        NumericVector row(n_eta);
        for (int j = 0; j < n_eta; ++j) row[j] = eta_warm(s, j);
        eta_init = row;
      }
      fwd = laplace_subject_forward(
          subj, theta, omega, sigma, grid, pred_lines,
          advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
          mode_centered, eta_init);
      total += fwd.val;
      fwd_out[s] = fwd_to_list(fwd);
      if (n_eta > 0) {
        for (int j = 0; j < n_eta; ++j) eta_modes_out(s, j) = fwd.eta_mode[j];
      }
    }
    if (need_grad) {
      SEXP eta_init = R_NilValue;
      if (eta_modes != R_NilValue) {
        NumericVector row(n_eta);
        for (int j = 0; j < n_eta; ++j) row[j] = eta_warm(s, j);
        eta_init = row;
      } else if (!use_fwd_cache && n_eta > 0) {
        eta_init = fwd.eta_mode;
      }
      NumericVector gt = clone(theta);
      NumericVector go = clone(omega);
      NumericVector gs = clone(sigma);
      if (mode_centered) {
        laplace_subject_grad_full(
            subj, grid, theta, omega, sigma, pred_lines,
            advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
            mode_centered, eta_init, gt, go, gs);
      } else {
        if (use_fwd_cache) {
          fwd = fwd_from_list(fwd_list[s]);
        }
        laplace_subject_grad(
            subj, fwd, theta, omega, sigma, pred_lines,
            advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
            gt, go, gs);
      }
      for (int i = 0; i < gt.size(); ++i) g_theta[i] += gt[i];
      for (int i = 0; i < go.size(); ++i) g_omega[i] += go[i];
      for (int i = 0; i < gs.size(); ++i) g_sigma[i] += gs[i];
    }
  }

  return List::create(
    Named("objective") = total,
    Named("grad_theta") = g_theta,
    Named("grad_omega") = g_omega,
    Named("grad_sigma") = g_sigma,
    Named("fwd_subjects") = use_fwd_cache ? fwd_list : fwd_out,
    Named("eta_modes") = eta_modes_out
  );
}

}  // namespace

// [[Rcpp::export]]
double nm_laplace_nll_cpp(
    List subjects,
    NumericVector theta,
    NumericVector omega,
    NumericVector sigma,
    NumericVector gh_nodes,
    NumericVector gh_weights,
    CharacterVector pred_lines,
    int advan,
    int trans,
    int obs_cmp,
    int dose_cmp,
    int n_transit,
    bool use_ode,
    int model_ss,
    bool mode_centered = true,
    int n_threads = 0,
    SEXP eta_modes = R_NilValue) {
  (void)n_threads;
  List res = laplace_engine(
      subjects, theta, omega, sigma, gh_nodes, gh_weights, pred_lines,
      advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
      mode_centered, false, eta_modes, R_NilValue, false);
  return as<double>(res["objective"]);
}

// [[Rcpp::export]]
List nm_laplace_nll_grad_cpp(
    List subjects,
    NumericVector theta,
    NumericVector omega,
    NumericVector sigma,
    NumericVector gh_nodes,
    NumericVector gh_weights,
    CharacterVector pred_lines,
    int advan,
    int trans,
    int obs_cmp,
    int dose_cmp,
    int n_transit,
    bool use_ode,
    int model_ss,
    bool mode_centered = true,
    int n_threads = 0,
    SEXP eta_modes = R_NilValue,
    SEXP fwd_cache = R_NilValue,
    bool grad_from_fwd = false) {
  (void)n_threads;
  return laplace_engine(
      subjects, theta, omega, sigma, gh_nodes, gh_weights, pred_lines,
      advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
      mode_centered, true, eta_modes, fwd_cache, grad_from_fwd);
}

// [[Rcpp::export]]
List nm_laplace_nll_detailed_cpp(
    List subjects,
    NumericVector theta,
    NumericVector omega,
    NumericVector sigma,
    NumericVector gh_nodes,
    NumericVector gh_weights,
    CharacterVector pred_lines,
    int advan,
    int trans,
    int obs_cmp,
    int dose_cmp,
    int n_transit,
    bool use_ode,
    int model_ss,
    bool mode_centered = true,
    SEXP eta_modes = R_NilValue) {
  return laplace_engine(
      subjects, theta, omega, sigma, gh_nodes, gh_weights, pred_lines,
      advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
      mode_centered, false, eta_modes, R_NilValue, false);
}
