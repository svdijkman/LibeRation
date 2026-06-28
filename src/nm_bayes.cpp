#include <Rcpp.h>
#include <cmath>
#include <vector>
#include "nm_bayes.h"
#include "nm_bayes_hmc.h"

using namespace Rcpp;
using nm_bayes::SubjectCache;
using nm_bayes::pop_log_prior;
using nm_bayes::refresh_subject_cache;
using nm_bayes::subject_total_nll;
using nm_bayes::eta_mh_one_step;

namespace {

double pop_nll_from_cache(const std::vector<SubjectCache>& caches) {
  double total = 0.0;
  for (const auto& c : caches) {
    total += subject_total_nll(c);
  }
  return total;
}

void refresh_all_caches(
    std::vector<SubjectCache>& caches,
    const List& subjects,
    const NumericMatrix& eta,
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
    const CharacterVector& des_lines,
    int& n_pk_eval) {
  const int n_sub = subjects.size();
  for (int s = 0; s < n_sub; ++s) {
    NumericVector eta_s(eta.ncol());
    for (int j = 0; j < eta.ncol(); ++j) {
      eta_s[j] = eta(s, j);
    }
    refresh_subject_cache(
        caches[s], subjects[s], eta_s, theta, omega, sigma, pred_lines,
        advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss, des_lines);
    n_pk_eval++;
  }
}

void update_prior_nll_only(
    std::vector<SubjectCache>& caches,
    const NumericMatrix& eta,
    const NumericVector& omega) {
  const int n_sub = caches.size();
  const int n_eta = omega.size();
  for (int s = 0; s < n_sub; ++s) {
    NumericVector eta_s(n_eta);
    for (int j = 0; j < n_eta; ++j) {
      eta_s[j] = eta(s, j);
    }
    caches[s].prior_nll = (n_eta > 0)
        ? nm_subject::omega_prior_nll(eta_s, omega)
        : 0.0;
  }
}

void update_resid_nll_only(
    std::vector<SubjectCache>& caches,
    const List& subjects,
    const NumericVector& sigma) {
  const double s1 = sigma.size() > 0 ? sigma[0] : 0.0;
  const double s2 = sigma.size() > 1 ? sigma[1] : 0.0;
  for (size_t s = 0; s < caches.size(); ++s) {
    List subj = subjects[s];
    IntegerVector dvid;
    if (subj.containsElementNamed("dvid")) {
      dvid = subj["dvid"];
    }
    caches[s].resid_nll = nm_subject::residual_nll(
        caches[s].dv, caches[s].f, s1, s2, nm_lik::config().error_type,
        dvid, sigma);
  }
}

double log_posterior(
    const std::vector<SubjectCache>& caches,
    const NumericVector& theta,
    const NumericVector& omega,
    const NumericVector& sigma,
    const IntegerVector& theta_prior_type,
    const NumericVector& theta_prior_mu,
    const NumericVector& theta_prior_sd,
    const IntegerVector& omega_prior_type,
    const NumericVector& omega_prior_mu,
    const NumericVector& omega_prior_sd,
    const IntegerVector& sigma_prior_type,
    const NumericVector& sigma_prior_mu,
    const NumericVector& sigma_prior_sd,
    const LogicalVector& theta_fix) {
  return -0.5 * pop_nll_from_cache(caches) +
      pop_log_prior(
          theta, omega, sigma,
          theta_prior_type, theta_prior_mu, theta_prior_sd,
          omega_prior_type, omega_prior_mu, omega_prior_sd,
          sigma_prior_type, sigma_prior_mu, sigma_prior_sd,
          theta_fix);
}

}  // namespace

// [[Rcpp::export]]
List nm_bayes_mcmc_cpp(
    NumericMatrix eta,
    NumericVector theta,
    NumericVector omega,
    NumericVector sigma,
    List subjects,
    CharacterVector pred_lines,
    int advan,
    int trans,
    int obs_cmp,
    int dose_cmp,
    int n_transit,
    bool use_ode,
    int model_ss,
    int n_burn = 100,
    int n_sample = 500,
    int n_thin = 1,
    NumericVector step_theta = NumericVector(),
    NumericVector step_omega = NumericVector(),
    NumericVector step_sigma = NumericVector(),
    double step_eta = 1.0,
    IntegerVector theta_prior_type = IntegerVector(),
    NumericVector theta_prior_mu = NumericVector(),
    NumericVector theta_prior_sd = NumericVector(),
    IntegerVector omega_prior_type = IntegerVector(),
    NumericVector omega_prior_mu = NumericVector(),
    NumericVector omega_prior_sd = NumericVector(),
    IntegerVector sigma_prior_type = IntegerVector(),
    NumericVector sigma_prior_mu = NumericVector(),
    NumericVector sigma_prior_sd = NumericVector(),
    LogicalVector theta_fix = LogicalVector(),
    CharacterVector des_lines = CharacterVector(),
    std::string sampler = "mh",
    double hmc_epsilon = 0.05,
    int hmc_leap = 10,
    int nuts_depth = 5) {
  const int n_sub = eta.nrow();
  const int n_eta = eta.ncol();
  const int n_th = theta.size();
  const int n_om = omega.size();
  const int n_sg = sigma.size();

  if (step_theta.size() == 0) {
    step_theta = NumericVector(n_th, 0.05);
  }
  if (step_omega.size() == 0) {
    step_omega = NumericVector(n_om, 0.05);
  }
  if (step_sigma.size() == 0) {
    step_sigma = NumericVector(n_sg, 0.05);
  }
  if (theta_fix.size() == 0) {
    theta_fix = LogicalVector(n_th, false);
  }

  const int n_keep = n_sample / std::max(n_thin, 1);
  NumericMatrix theta_chain(n_keep, n_th);
  NumericMatrix omega_chain(n_keep, n_om);
  NumericMatrix sigma_chain(n_keep, n_sg);
  NumericMatrix eta_chain;
  if (n_eta > 0) {
    eta_chain = NumericMatrix(n_keep, n_sub * n_eta);
  }

  IntegerVector acc_theta(n_th);
  IntegerVector acc_omega(n_om);
  IntegerVector acc_sigma(n_sg);
  int acc_eta = 0;
  int eta_props = 0;
  int n_pk_eval = 0;

  const bool use_hmc = (sampler == "hmc");
  const bool use_nuts = (sampler == "nuts");

  std::vector<SubjectCache> caches(n_sub);
  refresh_all_caches(
      caches, subjects, eta, theta, omega, sigma, pred_lines,
      advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
      des_lines, n_pk_eval);

  GetRNGstate();
  int keep_idx = 0;
  const int total_iter = n_burn + n_sample;

  for (int iter = 1; iter <= total_iter; ++iter) {
    // eta block
    for (int s = 0; s < n_sub; ++s) {
      NumericVector cur(n_eta);
      for (int j = 0; j < n_eta; ++j) {
        cur[j] = eta(s, j);
      }
      eta_props++;
      NumericVector upd;
      if (use_hmc) {
        upd = nm_bayes_hmc::hmc_eta_step(
            subjects[s], cur, theta, omega, sigma, pred_lines,
            advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
            hmc_epsilon, hmc_leap, n_pk_eval);
      } else if (use_nuts) {
        upd = nm_bayes_hmc::nuts_eta_step(
            subjects[s], cur, theta, omega, sigma, pred_lines,
            advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
            hmc_epsilon, nuts_depth, n_pk_eval);
      } else {
        upd = eta_mh_one_step(
            subjects[s], cur, caches[s], theta, omega, sigma, pred_lines,
            advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
            step_eta, n_pk_eval, des_lines);
      }
      if (!use_hmc && !use_nuts) {
        // caches updated inside eta_mh_one_step
      } else {
        refresh_subject_cache(
            caches[s], subjects[s], upd, theta, omega, sigma, pred_lines,
            advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
            des_lines);
        n_pk_eval++;
      }
      bool accepted = false;
      for (int j = 0; j < n_eta; ++j) {
        if (upd[j] != cur[j]) {
          accepted = true;
        }
        eta(s, j) = upd[j];
      }
      if (accepted) {
        acc_eta++;
      }
    }

    double log_post = log_posterior(
        caches, theta, omega, sigma,
        theta_prior_type, theta_prior_mu, theta_prior_sd,
        omega_prior_type, omega_prior_mu, omega_prior_sd,
        sigma_prior_type, sigma_prior_mu, sigma_prior_sd,
        theta_fix);

    // theta block (log-scale RW; restore subject caches on reject)
    for (int k = 0; k < n_th; ++k) {
      if (theta_fix[k]) continue;
      const double log_post_old = log_post;
      const std::vector<SubjectCache> caches_backup = caches;
      NumericVector th_try = clone(theta);
      const double log_th = std::log(std::max(theta[k], 1e-12));
      const double prop_log = log_th + step_theta[k] * R::rnorm(0.0, 1.0);
      th_try[k] = std::exp(prop_log);
      if (!std::isfinite(th_try[k]) || th_try[k] <= 0.0) {
        continue;
      }
      refresh_all_caches(
          caches, subjects, eta, th_try, omega, sigma, pred_lines,
          advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
          des_lines, n_pk_eval);
      const double log_post_new = log_posterior(
          caches, th_try, omega, sigma,
          theta_prior_type, theta_prior_mu, theta_prior_sd,
          omega_prior_type, omega_prior_mu, omega_prior_sd,
          sigma_prior_type, sigma_prior_mu, sigma_prior_sd,
          theta_fix);
      const double log_acc = log_post_new - log_post_old + prop_log - log_th;
      if (std::log(std::max(R::runif(0.0, 1.0), 1e-300)) < log_acc) {
        theta = th_try;
        log_post = log_post_new;
        acc_theta[k]++;
      } else {
        caches = caches_backup;
      }
    }

    // omega block (log-scale RW; update prior only)
    for (int j = 0; j < n_om; ++j) {
      const double log_post_old = log_post;
      NumericVector om_try = clone(omega);
      const double log_om = std::log(std::max(omega[j], 1e-12));
      const double prop_log = log_om + step_omega[j] * R::rnorm(0.0, 1.0);
      om_try[j] = std::exp(prop_log);
      if (!std::isfinite(om_try[j]) || om_try[j] <= 0.0) {
        continue;
      }
      update_prior_nll_only(caches, eta, om_try);
      const double log_post_new = log_posterior(
          caches, theta, om_try, sigma,
          theta_prior_type, theta_prior_mu, theta_prior_sd,
          omega_prior_type, omega_prior_mu, omega_prior_sd,
          sigma_prior_type, sigma_prior_mu, sigma_prior_sd,
          theta_fix);
      const double log_acc = log_post_new - log_post_old + prop_log - log_om;
      if (std::log(std::max(R::runif(0.0, 1.0), 1e-300)) < log_acc) {
        omega = om_try;
        log_post = log_post_new;
        acc_omega[j]++;
      } else {
        update_prior_nll_only(caches, eta, omega);
      }
    }

    // sigma block (log-scale RW; update residual only)
    for (int k = 0; k < n_sg; ++k) {
      const double log_post_old = log_post;
      NumericVector sg_try = clone(sigma);
      const double log_sg = std::log(std::max(sigma[k], 1e-12));
      const double prop_log = log_sg + step_sigma[k] * R::rnorm(0.0, 1.0);
      sg_try[k] = std::exp(prop_log);
      if (!std::isfinite(sg_try[k]) || sg_try[k] <= 0.0) {
        continue;
      }
      update_resid_nll_only(caches, subjects, sg_try);
      const double log_post_new = log_posterior(
          caches, theta, omega, sg_try,
          theta_prior_type, theta_prior_mu, theta_prior_sd,
          omega_prior_type, omega_prior_mu, omega_prior_sd,
          sigma_prior_type, sigma_prior_mu, sigma_prior_sd,
          theta_fix);
      const double log_acc = log_post_new - log_post_old + prop_log - log_sg;
      if (std::log(std::max(R::runif(0.0, 1.0), 1e-300)) < log_acc) {
        sigma = sg_try;
        log_post = log_post_new;
        acc_sigma[k]++;
      } else {
        update_resid_nll_only(caches, subjects, sigma);
      }
    }

    if (iter > n_burn && ((iter - n_burn - 1) % std::max(n_thin, 1) == 0) &&
        keep_idx < n_keep) {
      for (int k = 0; k < n_th; ++k) {
        theta_chain(keep_idx, k) = theta[k];
      }
      for (int j = 0; j < n_om; ++j) {
        omega_chain(keep_idx, j) = omega[j];
      }
      for (int k = 0; k < n_sg; ++k) {
        sigma_chain(keep_idx, k) = sigma[k];
      }
      if (n_eta > 0) {
        int col = 0;
        for (int s = 0; s < n_sub; ++s) {
          for (int j = 0; j < n_eta; ++j) {
            eta_chain(keep_idx, col++) = eta(s, j);
          }
        }
      }
      keep_idx++;
    }
  }
  PutRNGstate();

  return List::create(
      Named("theta") = theta,
      Named("omega") = omega,
      Named("sigma") = sigma,
      Named("eta") = eta,
      Named("theta_chain") = theta_chain,
      Named("omega_chain") = omega_chain,
      Named("sigma_chain") = sigma_chain,
      Named("eta_chain") = eta_chain,
      Named("log_posterior") = log_posterior(
          caches, theta, omega, sigma,
          theta_prior_type, theta_prior_mu, theta_prior_sd,
          omega_prior_type, omega_prior_mu, omega_prior_sd,
          sigma_prior_type, sigma_prior_mu, sigma_prior_sd,
          theta_fix),
      Named("acceptance") = List::create(
          Named("theta") = acc_theta,
          Named("omega") = acc_omega,
          Named("sigma") = acc_sigma,
          Named("eta") = NumericVector::create(
              Rcpp::_["accepted"] = acc_eta,
              Rcpp::_["proposed"] = eta_props
          )
      ),
      Named("n_burn") = n_burn,
      Named("n_sample") = n_sample,
      Named("n_thin") = n_thin,
      Named("n_pk_eval") = n_pk_eval,
      Named("n_mcmc_iter") = total_iter,
      Named("n_keep") = keep_idx
  );
}
