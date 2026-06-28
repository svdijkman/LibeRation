#pragma once
#include <Rcpp.h>
#include "nm_subject_nll.h"

namespace nm_bayes {

struct SubjectCache {
  Rcpp::NumericVector f;
  Rcpp::NumericVector dv;
  double resid_nll = 0.0;
  double prior_nll = 0.0;
};

inline double log_prior_flat() {
  return 0.0;
}

// Log-normal prior on log(x); x > 0
inline double log_prior_lognormal(double x, double mu_log, double sd_log) {
  if (x <= 0.0 || sd_log <= 0.0) {
    return R_NegInf;
  }
  const double lx = std::log(x);
  const double z = (lx - mu_log) / sd_log;
  return -lx - std::log(sd_log) - 0.5 * z * z - 0.5 * std::log(2.0 * 3.14159265358979323846);
}

inline double pop_log_prior(
    const Rcpp::NumericVector& theta,
    const Rcpp::NumericVector& omega,
    const Rcpp::NumericVector& sigma,
    const Rcpp::IntegerVector& theta_prior_type,
    const Rcpp::NumericVector& theta_prior_mu,
    const Rcpp::NumericVector& theta_prior_sd,
    const Rcpp::IntegerVector& omega_prior_type,
    const Rcpp::NumericVector& omega_prior_mu,
    const Rcpp::NumericVector& omega_prior_sd,
    const Rcpp::IntegerVector& sigma_prior_type,
    const Rcpp::NumericVector& sigma_prior_mu,
    const Rcpp::NumericVector& sigma_prior_sd,
    const Rcpp::LogicalVector& theta_fix) {
  double lp = 0.0;
  for (int i = 0; i < theta.size(); ++i) {
    if (theta_fix[i]) continue;
    if (theta_prior_type[i] == 1) {
      lp += log_prior_lognormal(theta[i], theta_prior_mu[i], theta_prior_sd[i]);
    }
  }
  for (int j = 0; j < omega.size(); ++j) {
    if (omega_prior_type[j] == 1) {
      lp += log_prior_lognormal(omega[j], omega_prior_mu[j], omega_prior_sd[j]);
    }
  }
  for (int k = 0; k < sigma.size(); ++k) {
    if (sigma_prior_type[k] == 1) {
      lp += log_prior_lognormal(sigma[k], sigma_prior_mu[k], sigma_prior_sd[k]);
    }
  }
  return lp;
}

inline void refresh_subject_cache(
    SubjectCache& cache,
    const Rcpp::List& subj,
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
    const Rcpp::CharacterVector& des_lines = Rcpp::CharacterVector()) {
  Rcpp::List cov = subj.containsElementNamed("cov")
      ? Rcpp::as<Rcpp::List>(subj["cov"]) : Rcpp::List();
  nm_subject::subject_f_at_obs(
      subj["time"], subj["amt"], subj["rate"], subj["f1"], subj["cmt"],
      subj["evid"], subj["ss"], subj["ii"], subj["dv"], subj["obs_idx"],
      eta, theta, pred_lines, advan, trans, obs_cmp, dose_cmp, n_transit,
      use_ode, model_ss, cache.f, cache.dv, cov, des_lines,
      nm_subject::extras_from_list(subj));
  const double s1 = sigma.size() > 0 ? sigma[0] : 0.0;
  const double s2 = sigma.size() > 1 ? sigma[1] : 0.0;
  cache.resid_nll = nm_subject::residual_nll(
      cache.dv, cache.f, s1, s2, nm_lik::config().error_type,
      subj.containsElementNamed("dvid")
          ? Rcpp::as<Rcpp::IntegerVector>(subj["dvid"]) : Rcpp::IntegerVector(),
      sigma);
  cache.prior_nll = (eta.size() > 0)
      ? nm_subject::omega_prior_nll(eta, omega)
      : 0.0;
}

inline double subject_total_nll(const SubjectCache& cache) {
  return cache.resid_nll + cache.prior_nll;
}

inline Rcpp::NumericVector eta_mh_one_step(
    const Rcpp::List& subj,
    Rcpp::NumericVector cur,
    SubjectCache& cache,
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
    double step_scale,
    int& n_pk_eval,
    const Rcpp::CharacterVector& des_lines = Rcpp::CharacterVector()) {
  const int n_eta = cur.size();
  const double n_old = subject_total_nll(cache);
  Rcpp::NumericVector prop = Rcpp::clone(cur);
  for (int j = 0; j < n_eta; ++j) {
    const double sd = std::sqrt(std::max(omega[j], 1e-8)) * step_scale;
    prop[j] += R::rnorm(0.0, sd);
  }
  SubjectCache prop_cache;
  refresh_subject_cache(
      prop_cache, subj, prop, theta, omega, sigma, pred_lines,
      advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss, des_lines);
  n_pk_eval++;
  const double n_new = subject_total_nll(prop_cache);
  if (std::log(std::max(R::runif(0.0, 1.0), 1e-300)) < n_old - n_new) {
    cache = prop_cache;
    return prop;
  }
  return cur;
}

}  // namespace nm_bayes
