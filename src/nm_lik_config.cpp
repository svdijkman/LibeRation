#include <Rcpp.h>
#include "nm_lik_config.h"

// [[Rcpp::export]]
void nm_lik_config_set(int error_type = 0, int omega_type = 0,
                       int sigma_corr = 0, int iov = 0,
                       double ar1_rho = 0.0, int blq_method = 0) {
  nm_lik::set_config(error_type, omega_type, sigma_corr, iov, ar1_rho,
                     blq_method);
}

// [[Rcpp::export]]
Rcpp::List nm_lik_config_get() {
  const nm_lik::LikConfig& c = nm_lik::config();
  return Rcpp::List::create(
      Rcpp::_["error_type"] = c.error_type,
      Rcpp::_["omega_type"] = c.omega_type,
      Rcpp::_["sigma_corr"] = c.sigma_corr,
      Rcpp::_["iov"] = c.iov,
      Rcpp::_["ar1_rho"] = c.ar1_rho,
      Rcpp::_["blq_method"] = c.blq_method
  );
}
