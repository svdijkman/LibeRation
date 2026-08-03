#ifndef LIBERATION_POPULATION_OBJECTIVE_API_HPP
#define LIBERATION_POPULATION_OBJECTIVE_API_HPP

#include <Rcpp.h>

namespace liberation {

SEXP population_objective_create_api(
    SEXP engine_pointer, const Rcpp::List& subject_data,
    const Rcpp::List& primary_tape_pointers,
    const Rcpp::List& curvature_tape_pointers,
    const Rcpp::List& config);
double population_objective_value_api(
    SEXP pointer, const Rcpp::NumericVector& encoded);
Rcpp::NumericVector population_objective_gradient_api(
    SEXP pointer, const Rcpp::NumericVector& encoded);
Rcpp::NumericMatrix population_objective_hessian_api(
    SEXP pointer, const Rcpp::NumericVector& encoded);
Rcpp::List population_objective_state_api(
    SEXP pointer, const Rcpp::NumericVector& encoded);
Rcpp::List population_objective_telemetry_api(SEXP pointer);

}  // namespace liberation

#endif
