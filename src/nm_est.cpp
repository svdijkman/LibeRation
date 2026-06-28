#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include "nm_pk_core.h"
#include "nm_subject_nll.h"
using namespace Rcpp;
using nm_subject::subject_nll;
using nm_subject::subject_nll_grad;
using nm_subject::subject_ipred;
using nm_subject::subject_f_at_obs;
using nm_subject::subject_nll_from_f;
using nm_subject::SubjectEventExtras;
using nm_subject::extras_from_list;
using nm_subject::logsumexp;

namespace {

Rcpp::List fwd_to_list(const Rcpp::NumericVector& f, const Rcpp::NumericVector& dv) {
  return Rcpp::List::create(
    Rcpp::Named("f") = f,
    Rcpp::Named("dv") = dv
  );
}

void subject_nll_grad_from_fwd(
    const Rcpp::NumericVector& dv_obs,
    const Rcpp::NumericVector& f,
    const Rcpp::NumericVector& eta,
    const Rcpp::NumericVector& theta,
    const Rcpp::NumericVector& omega,
    const Rcpp::NumericVector& sigma,
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
    const SubjectEventExtras& ex = SubjectEventExtras()) {
  const double s1 = sigma.size() > 0 ? sigma[0] : 0.0;
  const double s2 = sigma.size() > 1 ? sigma[1] : 0.0;
  const double eps = 1e-5;
  if (sigma.size() >= 1) {
    double gs1 = 0.0;
    double gs2 = 0.0;
    nm_subject::prop_add_nll_grad_sigma(dv_obs, f, s1, s2, gs1, gs2);
    g_sigma[0] += gs1;
    if (sigma.size() >= 2) {
      g_sigma[1] += gs2;
    }
  }
  if (include_omega_prior && eta.size() > 0) {
    for (int j = 0; j < omega.size(); ++j) {
      g_omega[j] += nm_subject::omega_prior_grad_j(eta[j], omega[j]);
    }
  }
  const double prior_q = (include_omega_prior && eta.size() > 0)
      ? nm_subject::omega_prior_nll(eta, omega)
      : 0.0;
  for (int k = 0; k < theta.size(); ++k) {
    Rcpp::NumericVector thp = Rcpp::clone(theta);
    Rcpp::NumericVector thm = Rcpp::clone(theta);
    thp[k] += eps;
    thm[k] -= eps;
    Rcpp::NumericVector fp;
    Rcpp::NumericVector fm;
    Rcpp::NumericVector dv_copy = dv_obs;
    nm_subject::subject_f_at_obs(
        time, amt, rate, f1, cmt, evid, ss, ii, dv, obs_idx, eta, thp,
        pred_lines, advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
        fp, dv_copy, Rcpp::List(), Rcpp::CharacterVector(), ex);
    nm_subject::subject_f_at_obs(
        time, amt, rate, f1, cmt, evid, ss, ii, dv, obs_idx, eta, thm,
        pred_lines, advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
        fm, dv_copy, Rcpp::List(), Rcpp::CharacterVector(), ex);
    const double nll_p = nm_subject::subject_nll_from_f(
        dv_obs, fp, eta, omega, sigma, false) + prior_q;
    const double nll_m = nm_subject::subject_nll_from_f(
        dv_obs, fm, eta, omega, sigma, false) + prior_q;
    g_theta[k] += (nll_p - nll_m) / (2.0 * eps);
  }
}

Rcpp::List pop_nll_engine(
    Rcpp::List subjects,
    Rcpp::NumericVector theta,
    Rcpp::NumericVector omega,
    Rcpp::NumericVector sigma,
    Rcpp::NumericMatrix eta,
    Rcpp::CharacterVector pred_lines,
    int advan,
    int trans,
    int obs_cmp,
    int dose_cmp,
    int n_transit,
    bool use_ode,
    int model_ss,
    bool include_omega_prior,
    bool need_grad,
    SEXP fwd_cache = R_NilValue,
    bool grad_from_fwd = false) {
  const int n_sub = subjects.size();
  double total = 0.0;
  Rcpp::NumericVector g_theta = Rcpp::clone(theta);
  Rcpp::NumericVector g_omega = Rcpp::clone(omega);
  Rcpp::NumericVector g_sigma = Rcpp::clone(sigma);
  std::fill(g_theta.begin(), g_theta.end(), 0.0);
  std::fill(g_omega.begin(), g_omega.end(), 0.0);
  std::fill(g_sigma.begin(), g_sigma.end(), 0.0);
  Rcpp::List fwd_out(n_sub);
  Rcpp::List fwd_list;
  const bool use_fwd_cache = grad_from_fwd && fwd_cache != R_NilValue;
  if (use_fwd_cache) {
    fwd_list = Rcpp::as<Rcpp::List>(fwd_cache);
    if (fwd_list.size() != n_sub) {
      Rcpp::stop("fwd_cache length mismatch");
    }
  }

  for (int s = 0; s < n_sub; ++s) {
    Rcpp::List subj = subjects[s];
    SubjectEventExtras ex = extras_from_list(subj);
    Rcpp::NumericVector eta_s(eta.ncol());
    for (int j = 0; j < eta.ncol(); ++j) {
      eta_s[j] = eta(s, j);
    }
    Rcpp::NumericVector f;
    Rcpp::NumericVector dv;
    if (use_fwd_cache) {
      Rcpp::List fc = fwd_list[s];
      f = fc["f"];
      dv = fc["dv"];
    } else {
      nm_subject::subject_f_at_obs(
          subj["time"], subj["amt"], subj["rate"], subj["f1"], subj["cmt"],
          subj["evid"], subj["ss"], subj["ii"], subj["dv"], subj["obs_idx"],
          eta_s, theta, pred_lines, advan, trans, obs_cmp, dose_cmp, n_transit,
          use_ode, model_ss, f, dv, Rcpp::List(), Rcpp::CharacterVector(), ex);
      fwd_out[s] = fwd_to_list(f, dv);
    }
    total += nm_subject::subject_nll_from_f(
        dv, f, eta_s, omega, sigma, include_omega_prior);
    if (need_grad) {
      subject_nll_grad_from_fwd(
          dv, f, eta_s, theta, omega, sigma,
          subj["time"], subj["amt"], subj["rate"], subj["f1"], subj["cmt"],
          subj["evid"], subj["ss"], subj["ii"], subj["dv"], subj["obs_idx"],
          pred_lines, advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
          include_omega_prior, g_theta, g_omega, g_sigma, ex);
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("objective") = total,
    Rcpp::Named("grad_theta") = g_theta,
    Rcpp::Named("grad_omega") = g_omega,
    Rcpp::Named("grad_sigma") = g_sigma,
    Rcpp::Named("fwd_subjects") = use_fwd_cache ? fwd_list : fwd_out
  );
}

}  // namespace

namespace {

NumericVector saem_mh_one_subject(
    const List& subj,
    NumericVector cur,
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
    int n_mcmc,
    double step_scale) {
  const NumericVector time = subj["time"];
  const NumericVector amt = subj["amt"];
  const NumericVector rate = subj["rate"];
  const NumericVector f1 = subj["f1"];
  const IntegerVector cmt = subj["cmt"];
  const IntegerVector evid = subj["evid"];
  const IntegerVector ss = subj["ss"];
  const NumericVector ii = subj["ii"];
  const NumericVector dv = subj["dv"];
  const IntegerVector obs_idx = subj["obs_idx"];
  const SubjectEventExtras ex = extras_from_list(subj);
  const int n_eta = cur.size();

  NumericVector f;
  NumericVector dv_obs;
  subject_f_at_obs(
      time, amt, rate, f1, cmt, evid, ss, ii, dv, obs_idx, cur, theta,
      pred_lines, advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
      f, dv_obs, Rcpp::List(), Rcpp::CharacterVector(), ex);
  double n_old = subject_nll_from_f(dv_obs, f, cur, omega, sigma, true);

  for (int m = 0; m < n_mcmc; ++m) {
    NumericVector prop = clone(cur);
    for (int j = 0; j < n_eta; ++j) {
      const double sd = std::sqrt(std::max(omega[j], 1e-8)) * step_scale;
      prop[j] += R::rnorm(0.0, sd);
    }
    NumericVector f_prop;
    subject_f_at_obs(
        time, amt, rate, f1, cmt, evid, ss, ii, dv, obs_idx, prop, theta,
        pred_lines, advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
        f_prop, dv_obs, Rcpp::List(), Rcpp::CharacterVector(), ex);
    const double n_new = subject_nll_from_f(dv_obs, f_prop, prop, omega, sigma, true);
    if (std::log(std::max(R::runif(0.0, 1.0), 1e-300)) < n_old - n_new) {
      cur = prop;
      f = f_prop;
      n_old = n_new;
    }
  }
  return cur;
}

}  // namespace

// [[Rcpp::export]]
double nm_subject_nll_cpp(
    NumericVector time,
    NumericVector amt,
    NumericVector rate,
    NumericVector f1,
    IntegerVector cmt,
    IntegerVector evid,
    IntegerVector ss,
    NumericVector ii,
    NumericVector dv,
    IntegerVector obs_idx,
    NumericVector eta,
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
    bool include_omega_prior = true) {
  return subject_nll(
      time, amt, rate, f1, cmt, evid, ss, ii, dv, obs_idx, eta, theta, omega, sigma,
      pred_lines, advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
      include_omega_prior);
}

// [[Rcpp::export]]
NumericMatrix nm_saem_mh_cpp(
    NumericMatrix eta,
    List subjects,
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
    int n_mcmc = 1,
    double step_scale = 1.0) {
  const int n_sub = eta.nrow();
  const int n_eta = eta.ncol();
  NumericMatrix out = clone(eta);
  GetRNGstate();
  for (int s = 0; s < n_sub; ++s) {
    NumericVector cur(n_eta);
    for (int j = 0; j < n_eta; ++j) {
      cur[j] = out(s, j);
    }
    cur = saem_mh_one_subject(
        subjects[s], cur, theta, omega, sigma, pred_lines,
        advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
        n_mcmc, step_scale);
    for (int j = 0; j < n_eta; ++j) {
      out(s, j) = cur[j];
    }
  }
  PutRNGstate();
  return out;
}

// [[Rcpp::export]]
double nm_nll_cpp(
    List subjects,
    NumericVector theta,
    NumericVector omega,
    NumericVector sigma,
    NumericMatrix eta,
    CharacterVector pred_lines,
    int advan,
    int trans,
    int obs_cmp,
    int dose_cmp,
    int n_transit,
    bool use_ode,
    int model_ss,
    bool include_omega_prior = true) {
  double total = 0.0;
  for (int s = 0; s < subjects.size(); ++s) {
    List subj = subjects[s];
    NumericVector eta_s(eta.ncol());
    for (int j = 0; j < eta.ncol(); ++j) eta_s[j] = eta(s, j);
    total += subject_nll(
        subj["time"], subj["amt"], subj["rate"], subj["f1"], subj["cmt"], subj["evid"],
        subj["ss"], subj["ii"],
        subj["dv"], subj["obs_idx"],
        eta_s, theta, omega, sigma, pred_lines,
        advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
        include_omega_prior, Rcpp::List(), Rcpp::IntegerVector(),
        Rcpp::CharacterVector(), extras_from_list(subj));
  }
  return total;
}

// [[Rcpp::export]]
List nm_nll_detailed_cpp(
    List subjects,
    NumericVector theta,
    NumericVector omega,
    NumericVector sigma,
    NumericMatrix eta,
    CharacterVector pred_lines,
    int advan,
    int trans,
    int obs_cmp,
    int dose_cmp,
    int n_transit,
    bool use_ode,
    int model_ss,
    bool include_omega_prior = true) {
  return pop_nll_engine(
      subjects, theta, omega, sigma, eta, pred_lines,
      advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
      include_omega_prior, false, R_NilValue, false);
}

// [[Rcpp::export]]
List nm_nll_grad_cpp(
    List subjects,
    NumericVector theta,
    NumericVector omega,
    NumericVector sigma,
    NumericMatrix eta,
    CharacterVector pred_lines,
    int advan,
    int trans,
    int obs_cmp,
    int dose_cmp,
    int n_transit,
    bool use_ode,
    int model_ss,
    bool include_omega_prior = true,
    SEXP fwd_cache = R_NilValue,
    bool grad_from_fwd = false) {
  return pop_nll_engine(
      subjects, theta, omega, sigma, eta, pred_lines,
      advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
      include_omega_prior, true, fwd_cache, grad_from_fwd);
}
