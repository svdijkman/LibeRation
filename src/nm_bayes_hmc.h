#pragma once
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include "nm_bayes.h"

namespace nm_bayes_hmc {

inline double subject_log_post_eta(
    const Rcpp::List& subj,
    const Rcpp::NumericVector& eta,
    const Rcpp::NumericVector& theta,
    const Rcpp::NumericVector& omega,
    const Rcpp::NumericVector& sigma,
    const Rcpp::CharacterVector& pred_lines,
    int advan, int trans, int obs_cmp, int dose_cmp,
    int n_transit, bool use_ode, int model_ss) {
  nm_bayes::SubjectCache cache;
  nm_bayes::refresh_subject_cache(
      cache, subj, eta, theta, omega, sigma, pred_lines,
      advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss);
  return -nm_bayes::subject_total_nll(cache);
}

inline void subject_grad_eta(
    const Rcpp::List& subj,
    const Rcpp::NumericVector& eta,
    const Rcpp::NumericVector& theta,
    const Rcpp::NumericVector& omega,
    const Rcpp::NumericVector& sigma,
    const Rcpp::CharacterVector& pred_lines,
    int advan, int trans, int obs_cmp, int dose_cmp,
    int n_transit, bool use_ode, int model_ss,
    Rcpp::NumericVector& grad) {
  const int p = eta.size();
  grad = Rcpp::NumericVector(p);
  const double eps = 1e-4;
  for (int j = 0; j < p; ++j) {
    Rcpp::NumericVector etap = Rcpp::clone(eta);
    Rcpp::NumericVector etam = Rcpp::clone(eta);
    etap[j] += eps;
    etam[j] -= eps;
    const double fp = subject_log_post_eta(
        subj, etap, theta, omega, sigma, pred_lines,
        advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss);
    const double fm = subject_log_post_eta(
        subj, etam, theta, omega, sigma, pred_lines,
        advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss);
    grad[j] = (fp - fm) / (2.0 * eps);
  }
}

inline Rcpp::NumericVector hmc_eta_step(
    const Rcpp::List& subj,
    const Rcpp::NumericVector& eta0,
    const Rcpp::NumericVector& theta,
    const Rcpp::NumericVector& omega,
    const Rcpp::NumericVector& sigma,
    const Rcpp::CharacterVector& pred_lines,
    int advan, int trans, int obs_cmp, int dose_cmp,
    int n_transit, bool use_ode, int model_ss,
    double epsilon,
    int n_leap,
    int& n_pk_eval) {
  const int p = eta0.size();
  Rcpp::NumericVector q = Rcpp::clone(eta0);
  Rcpp::NumericVector p_mom(p);
  for (int j = 0; j < p; ++j) {
    p_mom[j] = R::rnorm(0.0, 1.0);
  }
  const double H0 = subject_log_post_eta(
      subj, q, theta, omega, sigma, pred_lines,
      advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss)
      - 0.5 * Rcpp::sum(p_mom * p_mom);

  Rcpp::NumericVector grad(p);
  for (int step = 0; step < n_leap; ++step) {
    subject_grad_eta(
        subj, q, theta, omega, sigma, pred_lines,
        advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss, grad);
    n_pk_eval++;
    for (int j = 0; j < p; ++j) {
      p_mom[j] += 0.5 * epsilon * grad[j];
    }
    for (int j = 0; j < p; ++j) {
      q[j] += epsilon * p_mom[j];
    }
    subject_grad_eta(
        subj, q, theta, omega, sigma, pred_lines,
        advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss, grad);
    n_pk_eval++;
    for (int j = 0; j < p; ++j) {
      p_mom[j] += 0.5 * epsilon * grad[j];
    }
  }

  const double H1 = subject_log_post_eta(
      subj, q, theta, omega, sigma, pred_lines,
      advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss)
      - 0.5 * Rcpp::sum(p_mom * p_mom);
  if (std::log(std::max(R::runif(0.0, 1.0), 1e-300)) < H1 - H0) {
    return q;
  }
  return eta0;
}

inline Rcpp::NumericVector nuts_eta_step(
    const Rcpp::List& subj,
    const Rcpp::NumericVector& eta0,
    const Rcpp::NumericVector& theta,
    const Rcpp::NumericVector& omega,
    const Rcpp::NumericVector& sigma,
    const Rcpp::CharacterVector& pred_lines,
    int advan, int trans, int obs_cmp, int dose_cmp,
    int n_transit, bool use_ode, int model_ss,
    double epsilon,
    int max_depth,
    int& n_pk_eval) {
  const int p = eta0.size();
  Rcpp::NumericVector q0 = Rcpp::clone(eta0);
  Rcpp::NumericVector q = Rcpp::clone(eta0);
  Rcpp::NumericVector p_mom(p);
  for (int j = 0; j < p; ++j) {
    p_mom[j] = R::rnorm(0.0, 1.0);
  }
  Rcpp::NumericVector p_dir = Rcpp::clone(p_mom);
  const double H0 = subject_log_post_eta(
      subj, q, theta, omega, sigma, pred_lines,
      advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss)
      - 0.5 * Rcpp::sum(p_mom * p_mom);

  Rcpp::NumericVector grad(p);
  subject_grad_eta(
      subj, q, theta, omega, sigma, pred_lines,
      advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss, grad);
  n_pk_eval++;

  for (int j = 0; j < p; ++j) {
    p_mom[j] += 0.5 * epsilon * grad[j];
  }

  int n_steps = 1;
  for (int depth = 0; depth < max_depth; ++depth) {
    const int dir = (R::runif(0.0, 1.0) < 0.5) ? -1 : 1;
    for (int s = 0; s < n_steps; ++s) {
      for (int j = 0; j < p; ++j) {
        q[j] += dir * epsilon * p_dir[j];
      }
      subject_grad_eta(
          subj, q, theta, omega, sigma, pred_lines,
          advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss, grad);
      n_pk_eval++;
      for (int j = 0; j < p; ++j) {
        p_mom[j] += dir * epsilon * grad[j];
      }
    }
    n_steps *= 2;
    const double H1 = subject_log_post_eta(
        subj, q, theta, omega, sigma, pred_lines,
        advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss)
        - 0.5 * Rcpp::sum(p_mom * p_mom);
    if (H1 - H0 < -1e-10) {
      return q0;
    }
  }

  const double Hf = subject_log_post_eta(
      subj, q, theta, omega, sigma, pred_lines,
      advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss)
      - 0.5 * Rcpp::sum(p_mom * p_mom);
  if (std::log(std::max(R::runif(0.0, 1.0), 1e-300)) < Hf - H0) {
    return q;
  }
  return q0;
}

}  // namespace nm_bayes_hmc
