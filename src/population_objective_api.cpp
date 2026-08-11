#include "population_objective_api.h"
#include "execution_contract.h"

// This translation unit is the stable R-facing seam for the persistent
// population objective. The CppAD-heavy implementation remains private to the
// engine and can be decomposed without changing generated R entry points.

// [[Rcpp::export(name = ".liberation_population_objective_create")]]
SEXP liberation_population_objective_create(
    SEXP engine_pointer, const Rcpp::List& subject_data,
    const Rcpp::List& primary_tape_pointers,
    const Rcpp::List& curvature_tape_pointers,
    const Rcpp::List& config) {
  liberation::require_materialized_addl(subject_data);
  return liberation::population_objective_create_api(
    engine_pointer, subject_data, primary_tape_pointers,
    curvature_tape_pointers, config);
}

// [[Rcpp::export(name = ".liberation_population_objective_value")]]
double liberation_population_objective_value(
    SEXP pointer, const Rcpp::NumericVector& encoded) {
  return liberation::population_objective_value_api(pointer, encoded);
}

// [[Rcpp::export(name = ".liberation_population_objective_gradient")]]
Rcpp::NumericVector liberation_population_objective_gradient(
    SEXP pointer, const Rcpp::NumericVector& encoded) {
  return liberation::population_objective_gradient_api(pointer, encoded);
}

// [[Rcpp::export(name = ".liberation_population_objective_hessian")]]
Rcpp::NumericMatrix liberation_population_objective_hessian(
    SEXP pointer, const Rcpp::NumericVector& encoded) {
  return liberation::population_objective_hessian_api(pointer, encoded);
}

// [[Rcpp::export(name = ".liberation_population_objective_state")]]
Rcpp::List liberation_population_objective_state(
    SEXP pointer, const Rcpp::NumericVector& encoded) {
  return liberation::population_objective_state_api(pointer, encoded);
}

// [[Rcpp::export(name = ".liberation_population_objective_telemetry")]]
Rcpp::List liberation_population_objective_telemetry(SEXP pointer) {
  return liberation::population_objective_telemetry_api(pointer);
}

// Optimise a persistent compiled population objective without returning to R
// for each value/gradient evaluation.
// [[Rcpp::export(name = ".liberation_population_objective_native_optimizer")]]
Rcpp::List liberation_population_objective_native_optimizer(
    SEXP pointer, const Rcpp::NumericVector& start,
    const Rcpp::NumericVector& lower, const Rcpp::NumericVector& upper,
    int maxit = 200, double tolerance = 1e-6, int trace = 0) {
  return liberation::population_objective_native_optimizer_api(
    pointer, start, lower, upper, maxit, tolerance, trace);
}
