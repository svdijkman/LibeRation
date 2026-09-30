// [[Rcpp::depends(LibeRtAD)]]
// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>
#include <LibeRtAD/eigen_r.hpp>
#include <Eigen/Eigenvalues>
#include "execution_contract.h"
#include "population_objective_api.h"
#include "native_optimizer_api.h"
#include <LibeRtAD/sparse_hessian.hpp>
#include <unsupported/Eigen/MatrixFunctions>
#include <LibeRtAD/program.hpp>
#include <LibeRtAD/eigen_solver.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <cmath>
#include <condition_variable>
#include <functional>
#include <iomanip>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <queue>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace liberation {

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

// Internal implementation units are included into one coordinator TU so
// CppAD/Eigen template definitions remain visible without duplication.
// The population-objective R boundary is separately compiled in
// population_objective_api.cpp.
#include "pk_engine_support.h"
#include "pk_engine_advan_topology.h"
#include "pk_engine_model_engine.h"
#include "pk_engine_event_data.h"
#include "pk_engine_model_evaluation.h"
#include "pk_engine_differential_systems.h"
#include "pk_engine_ad_propagation.h"
#include "pk_engine_objective_tape.h"
#include "pk_engine_residual_models.h"
#include "pk_engine_fo_marginal.h"
#include "pk_engine_conditional.h"
#include "pk_engine_population.h"
#include "pk_engine_stochastic_support.h"
#include "pk_engine_stochastic_eta_context.h"
#include "pk_engine_quadrature.h"
#include "pk_engine_saem_context.h"
#include "pk_engine_weighted_eta_context.h"
#include "pk_engine_saem.h"
#include "pk_engine_state_space.h"

}  // namespace liberation

// Retain stochastic subject tapes, dynamic inputs, and point buffers across
// optimized SAEM/BAYES iterations.
// [[Rcpp::export(name = ".liberation_stochastic_eta_context_create")]]
SEXP liberation_stochastic_eta_context_create(
    SEXP engine_pointer, const Rcpp::List& tape_pointers,
    const Rcpp::List& subject_data, int n_theta, int n_eta, int n_sigma,
    int n_omega, bool use_ode, const Rcpp::NumericVector& initial_theta,
    const Rcpp::NumericVector& initial_sigma,
    const Rcpp::NumericVector& initial_omega, double guard_radius,
    bool fused_values, int native_threads) {
  auto context = std::make_unique<liberation::StochasticEtaCollection>(
    engine_pointer, tape_pointers, subject_data, n_theta, n_eta, n_sigma,
    n_omega, use_ode, initial_theta, initial_sigma, initial_omega,
    guard_radius, fused_values, native_threads);
  Rcpp::XPtr<liberation::StochasticEtaCollection> pointer(
    context.release(), true);
  pointer.attr("class") = Rcpp::CharacterVector::create(
    "liberation_stochastic_eta_context_ptr", "externalptr");
  return pointer;
}

// [[Rcpp::export(name = ".liberation_stochastic_eta_context_eval")]]
Rcpp::NumericVector liberation_stochastic_eta_context_eval(
    SEXP context_pointer, const Rcpp::NumericVector& theta,
    const Rcpp::NumericMatrix& eta, const Rcpp::NumericVector& sigma,
    const Rcpp::NumericVector& omega) {
  Rcpp::XPtr<liberation::StochasticEtaCollection> context(context_pointer);
  return context->evaluate(theta, eta, sigma, omega);
}

// Build all subject-specific Gaussian Laplace proposals in the persistent
// stochastic context.  Conditional modes use the same scale-aware convergence,
// restart, and bounded-curvature repair policy as the optimized Laplace
// estimator, while curvature is evaluated only at the accepted mode.
// [[Rcpp::export(name = ".liberation_stochastic_eta_context_laplace_proposal")]]
Rcpp::List liberation_stochastic_eta_context_laplace_proposal(
    SEXP context_pointer, const Rcpp::NumericVector& theta,
    const Rcpp::NumericMatrix& starts, const Rcpp::NumericVector& sigma,
    const Rcpp::NumericVector& omega, int maxit, double tolerance) {
  Rcpp::XPtr<liberation::StochasticEtaCollection> context(context_pointer);
  return context->laplace_proposal(
    theta, starts, sigma, omega, maxit, tolerance);
}

// [[Rcpp::export(name = ".liberation_stochastic_eta_context_random_walk")]]
Rcpp::List liberation_stochastic_eta_context_random_walk(
    SEXP context_pointer, const Rcpp::NumericVector& theta,
    const Rcpp::NumericMatrix& eta, const Rcpp::NumericVector& sigma,
    const Rcpp::NumericVector& omega, const Rcpp::List& proposal_roots,
    const Rcpp::NumericMatrix& normals,
    const Rcpp::NumericVector& log_uniforms, int mcmc_steps,
    double step_scale,
    Rcpp::Nullable<Rcpp::NumericVector> current_values = R_NilValue) {
  Rcpp::XPtr<liberation::StochasticEtaCollection> context(context_pointer);
  return context->random_walk(
    theta, eta, sigma, omega, proposal_roots, normals, log_uniforms,
    mcmc_steps, step_scale, current_values);
}

// [[Rcpp::export(name = ".liberation_stochastic_eta_context_independence")]]
Rcpp::List liberation_stochastic_eta_context_independence(
    SEXP context_pointer, const Rcpp::NumericVector& theta,
    const Rcpp::NumericMatrix& eta, const Rcpp::NumericVector& sigma,
    const Rcpp::NumericVector& omega,
    const Rcpp::NumericMatrix& proposal_modes,
    const Rcpp::List& proposal_roots,
    const Rcpp::List& proposal_precisions,
    const Rcpp::NumericMatrix& normals,
    const Rcpp::NumericVector& log_uniforms, int mcmc_steps,
    Rcpp::Nullable<Rcpp::NumericVector> current_values = R_NilValue,
    double proposal_df = 1e300,
    Rcpp::Nullable<Rcpp::NumericVector> proposal_scales = R_NilValue) {
  Rcpp::XPtr<liberation::StochasticEtaCollection> context(context_pointer);
  return context->laplace_independence(
    theta, eta, sigma, omega, proposal_modes, proposal_roots,
    proposal_precisions, normals, log_uniforms, mcmc_steps, current_values,
    proposal_df >= 1e250 ? R_PosInf : proposal_df, proposal_scales);
}

// [[Rcpp::export(name = ".liberation_stochastic_eta_context_bayes")]]
Rcpp::List liberation_stochastic_eta_context_bayes(
    SEXP context_pointer, const Rcpp::List& map_config,
    int n_burn, int n_sample, int n_thin, double step_scale,
    double eta_step, bool adapt, const std::string& outer_kernel,
    int adaptive_start, int adaptive_interval, double target_acceptance,
    double delayed_rejection_scale,
    const std::string& eta_kernel, int eta_refresh, int eta_maxit,
    double eta_tolerance, double eta_df, double eta_rescue_probability,
    double eta_parameter_refresh, double eta_low_acceptance,
    bool gibbs_omega) {
  Rcpp::RNGScope scope;
  Rcpp::XPtr<liberation::StochasticEtaCollection> context(context_pointer);
  return context->bayes_sample(
    map_config, n_burn, n_sample, n_thin, step_scale, eta_step, adapt,
    outer_kernel, adaptive_start, adaptive_interval, target_acceptance,
    delayed_rejection_scale,
    eta_kernel, eta_refresh, eta_maxit, eta_tolerance,
    eta_df, eta_rescue_probability, eta_parameter_refresh,
    eta_low_acceptance, gibbs_omega);
}

// Retain the complete adaptive/fixed Gaussian-quadrature coordinator across
// native optimizer callbacks.  The referenced stochastic subject context is
// preserved by the coordinator for its full lifetime, including refinement
// and final diagnostic evaluation.
// [[Rcpp::export(name = ".liberation_gq_context_create")]]
SEXP liberation_gq_context_create(
    SEXP stochastic_context_pointer, const Rcpp::List& map_config,
    const Rcpp::NumericMatrix& nodes,
    const Rcpp::NumericVector& log_measure,
    const Rcpp::NumericVector& measure_sign, bool adaptive,
    int eta_maxit, double tolerance) {
  Rcpp::XPtr<liberation::StochasticEtaCollection> stochastic(
    stochastic_context_pointer);
  auto context = std::make_unique<liberation::NativeGqCoordinator>(
    stochastic.get(), stochastic_context_pointer, map_config, nodes,
    log_measure, measure_sign, adaptive, eta_maxit, tolerance);
  Rcpp::XPtr<liberation::NativeGqCoordinator> pointer(context.release(), true);
  pointer.attr("class") = Rcpp::CharacterVector::create(
    "liberation_gq_context_ptr", "externalptr");
  return pointer;
}

// [[Rcpp::export(name = ".liberation_gq_context_eval")]]
Rcpp::List liberation_gq_context_eval(
    SEXP context_pointer, const Rcpp::NumericVector& encoded,
    bool gradient = true) {
  Rcpp::XPtr<liberation::NativeGqCoordinator> context(context_pointer);
  return context->evaluate(encoded, gradient);
}

// [[Rcpp::export(name = ".liberation_gq_context_optimize")]]
Rcpp::List liberation_gq_context_optimize(
    SEXP context_pointer, int maxit, int trace = 0,
    bool exact_refinement = true) {
  Rcpp::XPtr<liberation::NativeGqCoordinator> context(context_pointer);
  return context->optimize(maxit, trace, exact_refinement);
}

// [[Rcpp::export(name = ".liberation_stochastic_eta_context_telemetry")]]
Rcpp::List liberation_stochastic_eta_context_telemetry(SEXP context_pointer) {
  Rcpp::XPtr<liberation::StochasticEtaCollection> context(context_pointer);
  return context->telemetry();
}

// Retain fixed ETAs and objective-tape references across R optim() callbacks.
// [[Rcpp::export(name = ".liberation_saem_fixed_eta_context_create")]]
SEXP liberation_saem_fixed_eta_context_create(
    const Rcpp::List& tape_pointers, const Rcpp::NumericMatrix& eta,
    int n_theta, int n_sigma, int n_omega) {
  auto context = std::make_unique<liberation::SaemFixedEtaCollection>(
    tape_pointers, eta, n_theta, n_sigma, n_omega);
  Rcpp::XPtr<liberation::SaemFixedEtaCollection> pointer(
    context.release(), true);
  pointer.attr("class") = Rcpp::CharacterVector::create(
    "liberation_saem_fixed_eta_context_ptr", "externalptr");
  return pointer;
}

// [[Rcpp::export(name = ".liberation_saem_fixed_eta_context_eval")]]
Rcpp::List liberation_saem_fixed_eta_context_eval(
    SEXP context_pointer, const Rcpp::NumericVector& theta,
    const Rcpp::NumericVector& sigma, const Rcpp::NumericVector& omega) {
  Rcpp::XPtr<liberation::SaemFixedEtaCollection> context(context_pointer);
  return context->evaluate(theta, sigma, omega);
}

// Aggregate the same subject-ordered objective and gradient on the native
// side.  This avoids returning an N-subject gradient matrix on every R
// L-BFGS-B callback while retaining the established optimizer coordinator.
// [[Rcpp::export(name = ".liberation_saem_fixed_eta_context_eval_aggregate")]]
Rcpp::List liberation_saem_fixed_eta_context_eval_aggregate(
    SEXP context_pointer, const Rcpp::NumericVector& theta,
    const Rcpp::NumericVector& sigma, const Rcpp::NumericVector& omega) {
  Rcpp::XPtr<liberation::SaemFixedEtaCollection> context(context_pointer);
  return context->evaluate_aggregate(theta, sigma, omega);
}

// Persistent weighted complete-data expectation for ITS, IMP and SAEM.
// [[Rcpp::export(name = ".liberation_weighted_eta_context_create")]]
SEXP liberation_weighted_eta_context_create(
    SEXP engine_pointer, const Rcpp::List& tape_pointers,
    const Rcpp::List& subject_data,
    int n_theta, int n_eta, int n_sigma, int n_omega,
    bool use_ode = false, bool reduced_population_tape = false,
    int native_threads = 1, int ode_support_tape_limit = 4096) {
  auto context = std::make_unique<liberation::WeightedEtaCollection>(
    engine_pointer, tape_pointers, subject_data, n_theta, n_eta, n_sigma,
    n_omega, use_ode, reduced_population_tape, native_threads,
    ode_support_tape_limit);
  Rcpp::XPtr<liberation::WeightedEtaCollection> pointer(
    context.release(), true);
  pointer.attr("class") = Rcpp::CharacterVector::create(
    "liberation_weighted_eta_context_ptr", "externalptr");
  return pointer;
}

// [[Rcpp::export(name = ".liberation_weighted_eta_context_set")]]
void liberation_weighted_eta_context_set(
    SEXP context_pointer, const Rcpp::List& eta,
    const Rcpp::List& weights) {
  Rcpp::XPtr<liberation::WeightedEtaCollection> context(context_pointer);
  context->set_grids(eta, weights);
}

// Evaluate and normalize all IMP proposal points and install the resulting
// subject-specific probabilities without returning ETA/weight lists to R.
// [[Rcpp::export(name = ".liberation_weighted_eta_context_set_importance")]]
Rcpp::List liberation_weighted_eta_context_set_importance(
    SEXP context_pointer,
    const Rcpp::NumericVector& theta,
    const Rcpp::NumericVector& sigma,
    const Rcpp::NumericVector& omega,
    const Rcpp::List& eta,
    const Rcpp::List& log_proposal) {
  Rcpp::XPtr<liberation::WeightedEtaCollection> context(context_pointer);
  return context->set_importance(
    theta, sigma, omega, eta, log_proposal);
}

// [[Rcpp::export(name = ".liberation_weighted_eta_context_update")]]
void liberation_weighted_eta_context_update(
    SEXP context_pointer, const Rcpp::NumericMatrix& eta, double gamma,
    int max_support = 0, double prune_tolerance = 0.0) {
  Rcpp::XPtr<liberation::WeightedEtaCollection> context(context_pointer);
  context->update_common(eta, gamma, max_support, prune_tolerance);
}

// [[Rcpp::export(name = ".liberation_weighted_eta_context_eval")]]
Rcpp::List liberation_weighted_eta_context_eval(
    SEXP context_pointer, const Rcpp::NumericVector& theta,
    const Rcpp::NumericVector& sigma, const Rcpp::NumericVector& omega) {
  Rcpp::XPtr<liberation::WeightedEtaCollection> context(context_pointer);
  return context->evaluate_aggregate(theta, sigma, omega);
}

// [[Rcpp::export(name = ".liberation_weighted_eta_context_mean")]]
Rcpp::NumericMatrix liberation_weighted_eta_context_mean(
    SEXP context_pointer) {
  Rcpp::XPtr<liberation::WeightedEtaCollection> context(context_pointer);
  return context->mean_eta();
}

// [[Rcpp::export(name = ".liberation_weighted_eta_context_weights")]]
Rcpp::NumericVector liberation_weighted_eta_context_weights(
    SEXP context_pointer) {
  Rcpp::XPtr<liberation::WeightedEtaCollection> context(context_pointer);
  return context->common_weights();
}

// [[Rcpp::export(name = ".liberation_weighted_eta_context_recenter")]]
void liberation_weighted_eta_context_recenter(
    SEXP context_pointer, const Rcpp::NumericMatrix& adjustment) {
  Rcpp::XPtr<liberation::WeightedEtaCollection> context(context_pointer);
  context->recenter(adjustment);
}

// [[Rcpp::export(name = ".liberation_weighted_eta_context_omega")]]
Rcpp::NumericVector liberation_weighted_eta_context_omega(
    SEXP context_pointer, int n_eta_base, int iov,
    const Rcpp::IntegerVector& omega_rows,
    const Rcpp::IntegerVector& omega_cols) {
  Rcpp::XPtr<liberation::WeightedEtaCollection> context(context_pointer);
  return context->omega_sufficient(
    n_eta_base, iov, omega_rows, omega_cols);
}

// [[Rcpp::export(name = ".liberation_weighted_eta_context_sigma")]]
Rcpp::NumericVector liberation_weighted_eta_context_sigma(
    SEXP context_pointer, SEXP engine_pointer, const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta, const Rcpp::NumericVector& sigma) {
  Rcpp::XPtr<liberation::WeightedEtaCollection> context(context_pointer);
  return context->sigma_expectation(engine_pointer, data, theta, sigma);
}

// [[Rcpp::export(name = ".liberation_weighted_eta_context_telemetry")]]
Rcpp::List liberation_weighted_eta_context_telemetry(SEXP context_pointer) {
  Rcpp::XPtr<liberation::WeightedEtaCollection> context(context_pointer);
  return context->telemetry();
}

// Run one fixed-ETA SAEM maximisation entirely in C++.  Objective tapes,
// exact gradients, line search, and BFGS state remain on the native side.
// [[Rcpp::export(name = ".liberation_saem_mstep")]]
Rcpp::List liberation_saem_mstep(
    const Rcpp::List& tape_pointers, const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& theta, const Rcpp::NumericVector& sigma,
    const Rcpp::NumericVector& omega,
    const Rcpp::IntegerVector& theta_free,
    const Rcpp::IntegerVector& sigma_free,
    const Rcpp::NumericVector& lower,
    const Rcpp::NumericVector& upper,
    const Rcpp::List& prior_config,
    int maxit = 20, double tolerance = 1e-6, int trace = 0,
    SEXP optimizer_state = R_NilValue) {
  liberation::SaemFixedEtaObjective objective(
    tape_pointers, eta, theta, sigma, omega, theta_free, sigma_free,
    prior_config);
  SEXP state_result = optimizer_state;
  liberation::SaemLbfgsState* state = nullptr;
  if (Rf_isNull(optimizer_state)) {
    auto created = std::make_unique<liberation::SaemLbfgsState>();
    Rcpp::XPtr<liberation::SaemLbfgsState> pointer(created.release(), true);
    pointer.attr("class") = Rcpp::CharacterVector::create(
      "liberation_saem_lbfgs_state_ptr", "externalptr");
    state = pointer.get();
    state_result = pointer;
  } else {
    Rcpp::XPtr<liberation::SaemLbfgsState> pointer(optimizer_state);
    state = pointer.get();
  }
  Rcpp::List result = liberation::optimize_saem_fixed_eta(
    objective, lower, upper, maxit, tolerance, trace, *state);
  result["optimizer_state"] = state_result;
  return result;
}

// Run the native SAEM maximizer against the complete retained weighted-Q
// support.  The weighted context and L-BFGS history both remain resident, so
// post-burn iterations no longer fall back to R optimizer callbacks.
// [[Rcpp::export(name = ".liberation_saem_weighted_mstep")]]
Rcpp::List liberation_saem_weighted_mstep(
    SEXP weighted_context,
    const Rcpp::NumericVector& theta, const Rcpp::NumericVector& sigma,
    const Rcpp::NumericVector& omega,
    const Rcpp::IntegerVector& theta_free,
    const Rcpp::IntegerVector& sigma_free,
    const Rcpp::NumericVector& lower,
    const Rcpp::NumericVector& upper,
    const Rcpp::List& prior_config,
    int maxit = 20, double tolerance = 1e-6, int trace = 0,
    SEXP optimizer_state = R_NilValue) {
  Rcpp::XPtr<liberation::WeightedEtaCollection> weighted(weighted_context);
  liberation::SaemFixedEtaObjective objective(
    *weighted, theta, sigma, omega, theta_free, sigma_free, prior_config);
  SEXP state_result = optimizer_state;
  liberation::SaemLbfgsState* state = nullptr;
  if (Rf_isNull(optimizer_state)) {
    auto created = std::make_unique<liberation::SaemLbfgsState>();
    Rcpp::XPtr<liberation::SaemLbfgsState> pointer(created.release(), true);
    pointer.attr("class") = Rcpp::CharacterVector::create(
      "liberation_saem_lbfgs_state_ptr", "externalptr");
    state = pointer.get();
    state_result = pointer;
  } else {
    Rcpp::XPtr<liberation::SaemLbfgsState> pointer(optimizer_state);
    state = pointer.get();
  }
  Rcpp::List result = liberation::optimize_saem_fixed_eta(
    objective, lower, upper, maxit, tolerance, trace, *state);
  result["optimizer_state"] = state_result;
  result["backend"] = "native-cpp-weighted-eta-lbfgs";
  return result;
}

// [[Rcpp::export(name = ".liberation_saem_omega_sufficient")]]
Rcpp::NumericVector liberation_saem_omega_sufficient(
    const Rcpp::NumericMatrix& eta, int n_eta_base, int iov,
    const Rcpp::IntegerVector& omega_rows,
    const Rcpp::IntegerVector& omega_cols) {
  if (eta.nrow() < 1 || n_eta_base < 1 || iov < 0 || iov > n_eta_base ||
      omega_rows.size() != omega_cols.size()) {
    Rcpp::stop("Native SAEM OMEGA sufficient-statistic inputs are invalid.");
  }
  const int between = n_eta_base - iov;
  const int occasions = iov ? (eta.ncol() - between) / iov : 0;
  if ((!iov && eta.ncol() != n_eta_base) ||
      (iov && (eta.ncol() < between ||
               (eta.ncol() - between) % iov != 0 || occasions < 1))) {
    Rcpp::stop("Native SAEM ETA columns do not match the IOV layout.");
  }
  liberation::Matrix covariance = liberation::Matrix::Zero(
    n_eta_base, n_eta_base);
  if (!iov) {
    for (int subject = 0; subject < eta.nrow(); ++subject) {
      for (int row = 0; row < n_eta_base; ++row) {
        for (int column = 0; column < n_eta_base; ++column) {
          covariance(row, column) += eta(subject, row) * eta(subject, column);
        }
      }
    }
    covariance /= static_cast<double>(eta.nrow());
  } else {
    if (between) {
      for (int subject = 0; subject < eta.nrow(); ++subject) {
        for (int row = 0; row < between; ++row) {
          for (int column = 0; column < between; ++column) {
            covariance(row, column) += eta(subject, row) * eta(subject, column);
          }
        }
      }
      covariance.topLeftCorner(between, between) /=
        static_cast<double>(eta.nrow());
    }
    const double denominator = static_cast<double>(eta.nrow() * occasions);
    for (int subject = 0; subject < eta.nrow(); ++subject) {
      for (int occasion = 0; occasion < occasions; ++occasion) {
        const int offset = between + occasion * iov;
        for (int row = 0; row < iov; ++row) {
          for (int column = 0; column < iov; ++column) {
            covariance(between + row, between + column) +=
              eta(subject, offset + row) * eta(subject, offset + column);
          }
        }
      }
    }
    covariance.bottomRightCorner(iov, iov) /= denominator;
  }
  covariance.diagonal().array() += 1e-8;
  Rcpp::NumericVector result(omega_rows.size());
  for (R_xlen_t entry = 0; entry < omega_rows.size(); ++entry) {
    const int row = omega_rows[entry] - 1;
    const int column = omega_cols[entry] - 1;
    if (row < 0 || column < 0 || row >= n_eta_base || column >= n_eta_base) {
      Rcpp::stop("An OMEGA sufficient-statistic coordinate is invalid.");
    }
    result[entry] = covariance(row, column);
  }
  return result;
}

// [[Rcpp::export(name = ".liberation_engine_create")]]
SEXP liberation_engine_create(const Rcpp::List& specification) {
  Rcpp::XPtr<liberation::ModelEngine> pointer(
    new liberation::ModelEngine(specification), true
  );
  pointer.attr("class") = Rcpp::CharacterVector::create("liberation_engine_ptr", "externalptr");
  return pointer;
}

// [[Rcpp::export(name = ".liberation_engine_simulate")]]
Rcpp::List liberation_engine_simulate(
    SEXP engine_pointer,
    const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta,
    const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma) {
  Rcpp::XPtr<liberation::ModelEngine> engine(engine_pointer);
  liberation::require_materialized_addl(data);
  return liberation::simulate(*engine, data, theta, eta, sigma);
}

// Retain one canonical event table and its contiguous subject ranges.  The
// returned external pointer owns only R column references and integer ranges;
// no subject columns are duplicated.
// [[Rcpp::export(name = ".liberation_subject_store_create")]]
SEXP liberation_subject_store_create(
    const Rcpp::DataFrame& data, const Rcpp::IntegerVector& starts,
    const Rcpp::IntegerVector& lengths) {
  liberation::require_materialized_addl(data);
  Rcpp::XPtr<liberation::NativeSubjectStore> pointer(
    new liberation::NativeSubjectStore(data, starts, lengths), true);
  pointer.attr("class") = Rcpp::CharacterVector::create(
    "liberation_subject_store_ptr", "externalptr");
  return pointer;
}

// [[Rcpp::export(name = ".liberation_subject_view_signature")]]
std::string liberation_subject_view_signature(
    SEXP data_input, const Rcpp::CharacterVector& ignored_columns,
    bool include_fo_layout = false) {
  const liberation::EventDataView data = liberation::event_data_view(data_input);
  std::unordered_map<std::string, bool> ignored;
  for (R_xlen_t index = 0; index < ignored_columns.size(); ++index) {
    ignored[Rcpp::as<std::string>(ignored_columns[index])] = true;
  }
  return liberation::event_data_signature(data, ignored, include_fo_layout);
}

// Return only the columns/rows still needed by an R-side calculation.  This is
// deliberately a named list rather than a data frame, avoiding row.names,
// class reconstruction, and unrelated column copies.
// [[Rcpp::export(name = ".liberation_subject_view_project")]]
Rcpp::List liberation_subject_view_project(
    SEXP data_input, const Rcpp::CharacterVector& columns,
    bool observed_only = false, bool first_only = false) {
  const liberation::EventDataView data = liberation::event_data_view(data_input);
  std::vector<int> rows;
  rows.reserve(static_cast<std::size_t>(data.nrows()));
  for (int row = 0; row < data.nrows(); ++row) {
    if (observed_only && !(
        liberation::data_value(data, "EVID", row) == 0.0 &&
        liberation::data_value(data, "MDV", row) == 0.0 &&
        std::isfinite(liberation::data_value(data, "DV", row)))) continue;
    rows.push_back(row);
    if (first_only) break;
  }
  Rcpp::List result(columns.size());
  Rcpp::CharacterVector names(columns.size());
  for (R_xlen_t column = 0; column < columns.size(); ++column) {
    const std::string name = Rcpp::as<std::string>(columns[column]);
    names[column] = name;
    if (!data.containsElementNamed(name.c_str())) {
      result[column] = R_NilValue;
      continue;
    }
    Rcpp::NumericVector values(rows.size());
    for (std::size_t row = 0; row < rows.size(); ++row) {
      values[static_cast<R_xlen_t>(row)] =
        liberation::data_value(data, name, rows[row]);
    }
    result[column] = values;
  }
  result.attr("names") = names;
  result.attr("rows") = Rcpp::wrap(rows);
  result.attr("source") = "native-subject-view-projection";
  return result;
}

// Build the subject-specific random-effect covariance directly from a native
// row view.  This replaces R-side projection of .OCC_INDEX/.RE_TOTAL_* and
// avoids constructing an intermediate subject data frame.
// [[Rcpp::export(name = ".liberation_subject_effect_covariance")]]
Rcpp::NumericMatrix liberation_subject_effect_covariance(
    SEXP engine_pointer, SEXP data_input,
    const Rcpp::NumericVector& omega, int expanded_dimension) {
  Rcpp::XPtr<liberation::ModelEngine> engine(engine_pointer);
  const liberation::EventDataView data = liberation::event_data_view(data_input);
  std::vector<double> values(omega.begin(), omega.end());
  const liberation::MatrixT<double> base =
    liberation::omega_matrix_t<double>(*engine, values);
  const liberation::MatrixT<double> expanded =
    liberation::expanded_omega_t<double>(
      *engine, data, base, expanded_dimension);
  return libertad::eigen_matrix_to_r(expanded);
}

// Evaluate the established first-order Gaussian conditional covariance used
// by ITS in one native subject sweep. This changes execution only; the
// conditional-mode approximation and curvature formula remain unchanged.
// [[Rcpp::export(name = ".liberation_its_gaussian_covariance")]]
Rcpp::List liberation_its_gaussian_covariance(
    SEXP engine_pointer, const Rcpp::List& prediction_tapes,
    const Rcpp::List& subject_data, const Rcpp::NumericVector& theta,
    const Rcpp::NumericMatrix& modes, const Rcpp::NumericVector& sigma,
    const Rcpp::NumericVector& omega) {
  Rcpp::XPtr<liberation::ModelEngine> engine(engine_pointer);
  const int subjects = prediction_tapes.size();
  const int n_eta = modes.ncol();
  if (subjects < 1 || subject_data.size() != subjects ||
      modes.nrow() != subjects || n_eta < 1 || sigma.size() < 1 ||
      omega.size() != static_cast<R_xlen_t>(engine->omega_rows.size())) {
    Rcpp::stop("Native ITS Gaussian covariance inputs are inconsistent.");
  }
  if (engine->error_type != "additive" &&
      engine->error_type != "proportional" &&
      engine->error_type != "combined" &&
      engine->error_type != "exponential" &&
      engine->error_type != "power") {
    Rcpp::stop("Native ITS covariance requires a Gaussian residual model.");
  }
  const std::vector<double> sigma_native =
    Rcpp::as<std::vector<double>>(sigma);
  const liberation::Matrix base_omega =
    liberation::omega_matrix_t<double>(
      *engine, Rcpp::as<std::vector<double>>(omega));
  Rcpp::List covariance(subjects);
  Rcpp::NumericVector jitters(subjects);
  std::ostringstream messages;
  for (int subject = 0; subject < subjects; ++subject) {
    const liberation::EventDataView data =
      liberation::event_data_view(subject_data[subject]);
    liberation::require_materialized_addl(data);
    Rcpp::XPtr<liberation::PredictionTape> tape(prediction_tapes[subject]);
    if (tape->n_rows != data.nrows() ||
        static_cast<int>(tape->domain_names.size()) !=
          theta.size() + n_eta + sigma.size()) {
      Rcpp::stop("An ITS prediction tape has an inconsistent domain.");
    }
    std::vector<double> dynamic = liberation::prediction_dynamic_values(
      tape->dynamic_columns, data, tape->n_rows);
    if (!dynamic.empty()) tape->fun.new_dynamic(dynamic);
    tape->dynamic_values = dynamic;
    std::vector<double> point(tape->domain_names.size(), 0.0);
    std::copy(theta.begin(), theta.end(), point.begin());
    for (int effect = 0; effect < n_eta; ++effect) {
      const double value = modes(subject, effect);
      if (!std::isfinite(value)) {
        Rcpp::stop("ITS conditional modes must be finite.");
      }
      point[static_cast<std::size_t>(theta.size() + effect)] = value;
    }
    std::copy(sigma.begin(), sigma.end(),
              point.begin() + theta.size() + n_eta);
    const std::vector<double> prediction = tape->fun.Forward(0, point, messages);
    liberation::require_unchanged_path(
      tape->fun, "native ITS prediction evaluation");
    const std::size_t domain = tape->domain_names.size();
    std::vector<double> seed(domain * static_cast<std::size_t>(n_eta), 0.0);
    for (int effect = 0; effect < n_eta; ++effect) {
      seed[(static_cast<std::size_t>(theta.size()) +
            static_cast<std::size_t>(effect)) *
             static_cast<std::size_t>(n_eta) +
           static_cast<std::size_t>(effect)] = 1.0;
    }
    const std::vector<double> forward = n_eta == 1 ?
      tape->fun.Forward(1, seed) : tape->fun.Forward(1, n_eta, seed);
    liberation::require_unchanged_path(
      tape->fun, "native ITS prediction derivatives");
    const std::vector<int> observed = liberation::fo_observed_rows(data);
    const std::vector<int> dvid = liberation::fo_dvid_values(data);
    liberation::Matrix curvature = liberation::Matrix::Zero(n_eta, n_eta);
    for (int row : observed) {
      const double current = prediction[static_cast<std::size_t>(row)];
      const double variance = liberation::residual_variance_t<double>(
        *engine, current, sigma_native, dvid[static_cast<std::size_t>(row)]);
      liberation::Vector derivative(n_eta);
      for (int effect = 0; effect < n_eta; ++effect) {
        derivative[effect] = forward[
          static_cast<std::size_t>(row) * static_cast<std::size_t>(n_eta) +
          static_cast<std::size_t>(effect)];
      }
      curvature.noalias() += (2.0 / variance) *
        derivative * derivative.transpose();
    }
    const liberation::Matrix expanded = liberation::expanded_omega_t<double>(
      *engine, data, base_omega, n_eta);
    Eigen::LDLT<liberation::Matrix> omega_factor(expanded);
    if (omega_factor.info() != Eigen::Success) {
      Rcpp::stop("ITS OMEGA factorization failed.");
    }
    curvature.noalias() += 2.0 * omega_factor.solve(
      liberation::Matrix::Identity(n_eta, n_eta));
    curvature = 0.5 * (curvature + curvature.transpose()).eval();
    const auto eigen = libertad::detail::self_adjoint_eigen(curvature, false);
    if (eigen.info != Eigen::Success || !eigen.values.allFinite()) {
      Rcpp::stop("ITS first-order curvature decomposition failed.");
    }
    const double largest = std::max(eigen.values.cwiseAbs().maxCoeff(), 1.0);
    const double jitter = std::max(0.0, largest * 1e-9 - eigen.values.minCoeff());
    if (jitter > largest * 1e-2) {
      Rcpp::stop("ITS first-order conditional curvature is not sufficiently positive definite.");
    }
    curvature.diagonal().array() += jitter;
    Eigen::LDLT<liberation::Matrix> factor(curvature);
    if (factor.info() != Eigen::Success) {
      Rcpp::stop("ITS first-order curvature factorization failed.");
    }
    liberation::Matrix current_covariance = 2.0 * factor.solve(
      liberation::Matrix::Identity(n_eta, n_eta));
    current_covariance = 0.5 *
      (current_covariance + current_covariance.transpose()).eval();
    covariance[subject] = libertad::eigen_matrix_to_r(current_covariance);
    jitters[subject] = jitter;
    if ((subject + 1) % 32 == 0) Rcpp::checkUserInterrupt();
  }
  return Rcpp::List::create(
    Rcpp::Named("covariance") = covariance,
    Rcpp::Named("jitter") = jitters,
    Rcpp::Named("backend") = "cpp-batched-first-order-gaussian");
}

// [[Rcpp::export(name = ".liberation_subject_observation_count")]]
int liberation_subject_observation_count(SEXP data_input) {
  const liberation::EventDataView data = liberation::event_data_view(data_input);
  int observations = 0;
  for (int row = 0; row < data.nrows(); ++row) {
    if (liberation::data_value(data, "EVID", row) == 0.0 &&
        liberation::data_value(data, "MDV", row) == 0.0 &&
        std::isfinite(liberation::data_value(data, "DV", row))) {
      ++observations;
    }
  }
  return observations;
}

// [[Rcpp::export(name = ".liberation_subject_observation_counts")]]
Rcpp::IntegerVector liberation_subject_observation_counts(SEXP subject_source) {
  auto count_view = [](const liberation::EventDataView& data) {
    int observations = 0;
    for (int row = 0; row < data.nrows(); ++row) {
      if (liberation::data_value(data, "EVID", row) == 0.0 &&
          liberation::data_value(data, "MDV", row) == 0.0 &&
          std::isfinite(liberation::data_value(data, "DV", row))) {
        ++observations;
      }
    }
    return observations;
  };
  if (Rf_inherits(subject_source, "liberation_subject_store_ptr")) {
    Rcpp::XPtr<liberation::NativeSubjectStore> store(subject_source);
    Rcpp::IntegerVector result(store->starts.size());
    for (R_xlen_t subject = 0; subject < result.size(); ++subject) {
      result[subject] = count_view(store->view(static_cast<int>(subject)));
    }
    return result;
  }
  if (TYPEOF(subject_source) != VECSXP) {
    Rcpp::stop(
      "Observation counting requires a native subject store or subject-data list.");
  }
  Rcpp::List subjects(subject_source);
  Rcpp::IntegerVector result(subjects.size());
  for (R_xlen_t subject = 0; subject < subjects.size(); ++subject) {
    result[subject] = count_view(
      liberation::event_data_view(subjects[subject]));
  }
  return result;
}

// Draw a batch from a positive-semidefinite multivariate Gaussian in one
// native allocation. Small round-off eigenvalues are clipped; materially
// indefinite covariance remains an error.
// [[Rcpp::export(name = ".liberation_mvn_draws")]]
Rcpp::NumericMatrix liberation_mvn_draws(
    const Rcpp::NumericVector& mean_input,
    const Rcpp::NumericMatrix& covariance_input, int draws) {
  if (draws < 1 || covariance_input.nrow() != mean_input.size() ||
      covariance_input.ncol() != mean_input.size()) {
    Rcpp::stop("Multivariate-normal draw dimensions are invalid.");
  }
  const Eigen::VectorXd mean = libertad::r_vector_map(mean_input);
  const Eigen::MatrixXd covariance_source =
    libertad::r_matrix_map(covariance_input);
  const Eigen::MatrixXd covariance = 0.5 * (
    covariance_source + covariance_source.transpose()).eval();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(covariance);
  if (solver.info() != Eigen::Success || !mean.allFinite() ||
      !covariance.allFinite()) {
    Rcpp::stop("ETA covariance eigendecomposition failed.");
  }
  Eigen::VectorXd values = solver.eigenvalues();
  const double scale = std::max(1.0, values.cwiseAbs().maxCoeff());
  if (values.minCoeff() < -std::sqrt(std::numeric_limits<double>::epsilon()) * scale) {
    Rcpp::stop("ETA covariance is not positive semidefinite.");
  }
  values = values.cwiseMax(0.0).cwiseSqrt();
  const Eigen::MatrixXd root = solver.eigenvectors() * values.asDiagonal();
  Rcpp::NumericMatrix output(draws, mean.size());
  Eigen::Map<Eigen::MatrixXd> result(output.begin(), draws, mean.size());
  for (int row = 0; row < draws; ++row) {
    Eigen::VectorXd normal(mean.size());
    for (Eigen::Index column = 0; column < mean.size(); ++column) {
      normal[column] = R::rnorm(0.0, 1.0);
    }
    result.row(row) = (mean + root * normal).transpose();
  }
  return output;
}

// Preserve R's replicate-wise, column-major normal stream while avoiding one
// R matrix allocation and BLAS dispatch per simulation replicate.
// [[Rcpp::export(name = ".liberation_eta_draw_pool")]]
Rcpp::NumericMatrix liberation_eta_draw_pool(
    const Rcpp::NumericMatrix& root_input, int subjects, int replicates) {
  if (subjects < 1 || replicates < 1 ||
      root_input.nrow() != root_input.ncol()) {
    Rcpp::stop("ETA draw-pool dimensions are invalid.");
  }
  const Eigen::MatrixXd root = libertad::r_matrix_map(root_input);
  const int effects = root.rows();
  Rcpp::NumericMatrix output(subjects * replicates, effects);
  Eigen::Map<Eigen::MatrixXd> result(
    output.begin(), subjects * replicates, effects);
  Eigen::MatrixXd normal(subjects, effects);
  for (int replicate = 0; replicate < replicates; ++replicate) {
    for (int effect = 0; effect < effects; ++effect) {
      for (int subject = 0; subject < subjects; ++subject) {
        normal(subject, effect) = R::rnorm(0.0, 1.0);
      }
    }
    result.middleRows(replicate * subjects, subjects).noalias() =
      normal * root.transpose();
  }
  return output;
}

// Build ADDL expansion indices once. R subsets every user column with this
// layout in one operation instead of constructing and rbinding one data frame
// per generated dose.
// [[Rcpp::export(name = ".liberation_expand_addl_layout")]]
Rcpp::List liberation_expand_addl_layout(
    const Rcpp::NumericVector& time, const Rcpp::IntegerVector& evid,
    const Rcpp::NumericVector& amount, const Rcpp::NumericVector& interval,
    const Rcpp::IntegerVector& addl,
    const Rcpp::IntegerVector& source_row) {
  const R_xlen_t rows = time.size();
  if (evid.size() != rows || amount.size() != rows || interval.size() != rows ||
      addl.size() != rows || source_row.size() != rows) {
    Rcpp::stop("ADDL expansion columns have inconsistent lengths.");
  }
  std::size_t output_rows = static_cast<std::size_t>(rows);
  for (R_xlen_t row = 0; row < rows; ++row) {
    const int count = addl[row];
    if (count == NA_INTEGER || count < 0) {
      Rcpp::stop("ADDL must be a non-negative integer at row " +
        std::to_string(row + 1) + ".");
    }
    if (count && (!std::isfinite(interval[row]) || interval[row] <= 0.0)) {
      Rcpp::stop("ADDL > 0 requires II > 0 at row " +
        std::to_string(row + 1) + ".");
    }
    if (count && (!((evid[row] == 1) || (evid[row] == 4)) || amount[row] <= 0.0)) {
      Rcpp::stop("ADDL is only valid on positive dosing records at row " +
        std::to_string(row + 1) + ".");
    }
    output_rows += static_cast<std::size_t>(count);
  }
  Rcpp::IntegerVector source(output_rows), generated_source(output_rows),
    priority(output_rows);
  Rcpp::LogicalVector generated(output_rows);
  Rcpp::NumericVector expanded_time(output_rows);
  std::size_t position = 0;
  for (R_xlen_t row = 0; row < rows; ++row) {
    source[position] = row + 1;
    expanded_time[position] = time[row];
    generated[position] = false;
    generated_source[position] = source_row[row];
    priority[position] = 0;
    ++position;
    for (int repeat = 1; repeat <= addl[row]; ++repeat) {
      source[position] = row + 1;
      expanded_time[position] = time[row] + repeat * interval[row];
      generated[position] = true;
      generated_source[position] = source_row[row];
      priority[position] = -1;
      ++position;
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("source") = source,
    Rcpp::Named("time") = expanded_time,
    Rcpp::Named("generated") = generated,
    Rcpp::Named("source_row") = generated_source,
    Rcpp::Named("sort_priority") = priority);
}

// [[Rcpp::export(name = ".liberation_event_order")]]
Rcpp::IntegerVector liberation_event_order(
    const Rcpp::IntegerVector& id, const Rcpp::NumericVector& time,
    const Rcpp::IntegerVector& priority,
    const Rcpp::IntegerVector& source_row) {
  const R_xlen_t rows = id.size();
  if (time.size() != rows || priority.size() != rows ||
      source_row.size() != rows) {
    Rcpp::stop("Event-order columns have inconsistent lengths.");
  }
  std::vector<int> order(static_cast<std::size_t>(rows));
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](int left, int right) {
    if (id[left] != id[right]) return id[left] < id[right];
    if (time[left] != time[right]) return time[left] < time[right];
    if (priority[left] != priority[right]) return priority[left] < priority[right];
    return source_row[left] < source_row[right];
  });
  Rcpp::IntegerVector result(rows);
  for (R_xlen_t row = 0; row < rows; ++row) result[row] = order[row] + 1;
  return result;
}

// [[Rcpp::export(name = ".liberation_mu_program_create")]]
SEXP liberation_mu_program_create(
    const Rcpp::List& ir, const Rcpp::IntegerVector& eta) {
  Rcpp::XPtr<liberation::NativeMuProgram> pointer(
    new liberation::NativeMuProgram(ir, eta), true);
  pointer.attr("class") = Rcpp::CharacterVector::create(
    "liberation_mu_program_ptr", "externalptr");
  return pointer;
}

namespace {

void liberation_mu_fill_row(
    const liberation::NativeMuProgram& program,
    const liberation::EventDataView& data,
    const Rcpp::NumericVector& theta, Rcpp::NumericMatrix& result,
    int subject) {
  const std::vector<double> values =
    liberation::evaluate_mu_program(program, data, theta);
  for (std::size_t output = 0; output < values.size(); ++output) {
    const int eta = program.eta[output];
    if (eta < 0 || eta >= result.ncol()) {
      Rcpp::stop("A compiled MU ETA destination exceeds the expanded ETA dimension.");
    }
    result(subject, eta) = values[output];
  }
}

Rcpp::NumericMatrix liberation_mu_evaluate_source(
    const liberation::NativeMuProgram& program, SEXP subject_source,
    const Rcpp::NumericVector& theta, int n_eta) {
  if (n_eta < 0) Rcpp::stop("The MU ETA dimension cannot be negative.");
  if (Rf_inherits(subject_source, "liberation_subject_store_ptr")) {
    Rcpp::XPtr<liberation::NativeSubjectStore> store(subject_source);
    Rcpp::NumericMatrix result(
      static_cast<int>(store->starts.size()), n_eta);
    for (int subject = 0; subject < result.nrow(); ++subject) {
      liberation_mu_fill_row(
        program, store->view(subject), theta, result, subject);
    }
    return result;
  }
  if (TYPEOF(subject_source) != VECSXP) {
    Rcpp::stop("MU evaluation requires a native subject store or subject-data list.");
  }
  Rcpp::List subjects(subject_source);
  Rcpp::NumericMatrix result(subjects.size(), n_eta);
  for (R_xlen_t subject = 0; subject < subjects.size(); ++subject) {
    liberation_mu_fill_row(
      program, liberation::event_data_view(subjects[subject]), theta,
      result, static_cast<int>(subject));
  }
  return result;
}

double liberation_mu_inverse_link(double beta, const std::string& link) {
  if (link == "identity") return beta;
  if (link == "log") return std::exp(beta);
  Rcpp::stop("Unsupported native MU link: " + link + ".");
  return beta;
}

}  // namespace

// Evaluate all subject MU values while reading covariates by row reference.
// [[Rcpp::export(name = ".liberation_mu_program_eval")]]
Rcpp::NumericMatrix liberation_mu_program_eval(
    SEXP program_pointer, SEXP subject_source,
    const Rcpp::NumericVector& theta, int n_eta) {
  Rcpp::XPtr<liberation::NativeMuProgram> program(program_pointer);
  return liberation_mu_evaluate_source(*program, subject_source, theta, n_eta);
}

// Derive the exact affine MU offset/design through basis evaluations of the
// compiled equations.  R establishes that the equations are affine; C++ then
// evaluates every subject without projecting covariate columns into R.
// [[Rcpp::export(name = ".liberation_mu_affine_design")]]
Rcpp::List liberation_mu_affine_design(
    SEXP program_pointer, SEXP subject_source,
    const Rcpp::NumericVector& theta,
    const Rcpp::IntegerVector& theta_indices,
    const Rcpp::CharacterVector& links, int n_eta) {
  if (theta_indices.size() != links.size()) {
    Rcpp::stop("MU THETA indices and links must have equal length.");
  }
  Rcpp::XPtr<liberation::NativeMuProgram> program(program_pointer);
  Rcpp::NumericVector baseline_theta = Rcpp::clone(theta);
  for (R_xlen_t column = 0; column < theta_indices.size(); ++column) {
    const int index = theta_indices[column] - 1;
    if (index < 0 || index >= baseline_theta.size()) {
      Rcpp::stop("A MU design THETA index exceeds supplied values.");
    }
    baseline_theta[index] = liberation_mu_inverse_link(
      0.0, Rcpp::as<std::string>(links[column]));
  }
  Rcpp::NumericMatrix offset = liberation_mu_evaluate_source(
    *program, subject_source, baseline_theta, n_eta);
  Rcpp::List design(theta_indices.size());
  for (R_xlen_t column = 0; column < theta_indices.size(); ++column) {
    Rcpp::NumericVector basis_theta = Rcpp::clone(baseline_theta);
    const int index = theta_indices[column] - 1;
    basis_theta[index] = liberation_mu_inverse_link(
      1.0, Rcpp::as<std::string>(links[column]));
    Rcpp::NumericMatrix basis = liberation_mu_evaluate_source(
      *program, subject_source, basis_theta, n_eta);
    Rcpp::NumericMatrix current(basis.nrow(), basis.ncol());
    for (R_xlen_t value = 0; value < basis.size(); ++value) {
      current[value] = basis[value] - offset[value];
    }
    design[column] = current;
  }
  return Rcpp::List::create(
    Rcpp::Named("offset") = offset,
    Rcpp::Named("design_columns") = design,
    Rcpp::Named("backend") = "cpp-native-row-view");
}

// [[Rcpp::export(name = ".liberation_saem_sigma_sufficient")]]
Rcpp::NumericVector liberation_saem_sigma_sufficient(
    SEXP engine_pointer, const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta, const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma) {
  Rcpp::XPtr<liberation::ModelEngine> engine(engine_pointer);
  liberation::require_materialized_addl(data);
  if (engine->error_type != "additive" &&
      engine->error_type != "proportional" &&
      engine->error_type != "exponential") {
    Rcpp::stop("Native SAEM SIGMA updates require a simple residual model.");
  }
  const Rcpp::List simulation = liberation::simulate(
    *engine, data, theta, eta, sigma);
  const Rcpp::NumericVector prediction = simulation["ipred"];
  const Rcpp::IntegerVector evid = data["EVID"];
  const Rcpp::IntegerVector mdv = data["MDV"];
  const Rcpp::NumericVector dv = data["DV"];
  Rcpp::IntegerVector dvid(data.nrows(), 1);
  if (data.containsElementNamed("DVID")) dvid = data["DVID"];
  Rcpp::NumericVector result = Rcpp::clone(sigma);
  std::vector<double> sum_square(static_cast<std::size_t>(sigma.size()), 0.0);
  std::vector<int> count(static_cast<std::size_t>(sigma.size()), 0);
  for (int row = 0; row < data.nrows(); ++row) {
    if (evid[row] != 0 || mdv[row] != 0 ||
        !std::isfinite(dv[row]) || !std::isfinite(prediction[row])) continue;
    const int response = std::max(dvid[row], 1) - 1;
    if (response < 0 || response >= sigma.size()) continue;
    double residual = 0.0;
    if (engine->error_type == "additive") {
      residual = dv[row] - prediction[row];
    } else if (engine->error_type == "proportional") {
      residual = (dv[row] - prediction[row]) /
        std::max(std::abs(prediction[row]), 1e-12);
    } else {
      if (!(dv[row] > 0.0) || !(prediction[row] > 0.0)) continue;
      residual = std::log(dv[row]) - std::log(prediction[row]);
    }
    if (!std::isfinite(residual)) continue;
    sum_square[static_cast<std::size_t>(response)] += residual * residual;
    ++count[static_cast<std::size_t>(response)];
  }
  for (R_xlen_t response = 0; response < sigma.size(); ++response) {
    const int observations = count[static_cast<std::size_t>(response)];
    if (!observations) continue;
    const double variance =
      sum_square[static_cast<std::size_t>(response)] / observations;
    if (std::isfinite(variance) && variance > 0.0) {
      result[response] = engine->sigma_parameterization == "variance" ?
        variance : std::sqrt(variance);
    }
  }
  return result;
}

// [[Rcpp::export(name = ".liberation_engine_simulate_batch")]]
Rcpp::List liberation_engine_simulate_batch(
    SEXP engine_pointer,
    const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta,
    const Rcpp::List& eta_values,
    const Rcpp::NumericVector& sigma) {
  Rcpp::XPtr<liberation::ModelEngine> engine(engine_pointer);
  liberation::require_materialized_addl(data);
  Rcpp::List result(eta_values.size());
  for (R_xlen_t replicate = 0; replicate < eta_values.size(); ++replicate) {
    Rcpp::NumericMatrix eta(eta_values[replicate]);
    result[replicate] = liberation::simulate(*engine, data, theta, eta, sigma);
    if ((replicate + 1) % 32 == 0) Rcpp::checkUserInterrupt();
  }
  return result;
}

// [[Rcpp::export(name = ".liberation_engine_hmm_filter")]]
Rcpp::List liberation_engine_hmm_filter(
    SEXP engine_pointer,
    const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta,
    const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma) {
  Rcpp::XPtr<liberation::ModelEngine> engine(engine_pointer);
  liberation::require_materialized_addl(data);
  return liberation::hmm_filter(*engine, data, theta, eta, sigma);
}

// [[Rcpp::export(name = ".liberation_engine_kalman_filter")]]
Rcpp::List liberation_engine_kalman_filter(
    SEXP engine_pointer,
    const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta,
    const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma) {
  Rcpp::XPtr<liberation::ModelEngine> engine(engine_pointer);
  liberation::require_materialized_addl(data);
  return liberation::kalman_filter(*engine, data, theta, eta, sigma);
}

// [[Rcpp::export(name = ".liberation_engine_kalman_simulate")]]
Rcpp::NumericVector liberation_engine_kalman_simulate(
    SEXP engine_pointer, const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta, const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma,
    const Rcpp::NumericMatrix& process_normals,
    const Rcpp::NumericVector& observation_normals) {
  Rcpp::XPtr<liberation::ModelEngine> engine(engine_pointer);
  liberation::require_materialized_addl(data);
  return liberation::kalman_simulate(
    *engine, data, theta, eta, sigma, process_normals, observation_normals);
}

// [[Rcpp::export(name = ".liberation_engine_derivative")]]
Rcpp::NumericVector liberation_engine_derivative(
    SEXP engine_pointer,
    const Rcpp::DataFrame& data,
    int row,
    int subject,
    double time,
    const Rcpp::NumericVector& state,
    const Rcpp::NumericVector& theta,
    const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma) {
  Rcpp::XPtr<liberation::ModelEngine> engine(engine_pointer);
  liberation::require_materialized_addl(data);
  if (row < 1 || row > data.nrows()) Rcpp::stop("Derivative row is outside the dataset.");
  if (subject < 1 || subject > eta.nrow()) Rcpp::stop("Derivative subject is outside the ETA matrix.");
  if (state.size() != engine->n_state) Rcpp::stop("Derivative state has the wrong length.");
  liberation::Parameters parameters = liberation::evaluate_parameters(
    *engine, data, row - 1, subject - 1, theta, eta, sigma
  );
  const auto mapped = libertad::r_vector_map(state);
  return libertad::eigen_vector_to_r(liberation::evaluate_derivatives(
    *engine, data, row - 1, subject - 1, time, mapped, parameters, theta, eta, sigma
  ));
}

// [[Rcpp::export(name = ".liberation_matrix_exp")]]
Rcpp::NumericMatrix liberation_matrix_exp(const Rcpp::NumericMatrix& matrix,
                                           double dt = 1.0) {
  const auto mapped = libertad::r_matrix_map(matrix);
  return libertad::eigen_matrix_to_r(liberation::matrix_exp(mapped * dt));
}

// [[Rcpp::export(name = ".liberation_advan_matrix")]]
Rcpp::List liberation_advan_matrix(int advan, const Rcpp::List& parameters) {
  liberation::Parameters p;
  Rcpp::CharacterVector names = parameters.names();
  for (R_xlen_t i = 0; i < parameters.size(); ++i) {
    p[Rcpp::as<std::string>(names[i])] = Rcpp::as<double>(parameters[i]);
  }
  liberation::Topology topology = liberation::build_topology(advan, p);
  return Rcpp::List::create(
    Rcpp::Named("K") = libertad::eigen_matrix_to_r(topology.k),
    Rcpp::Named("states") = topology.state_names
  );
}

// [[Rcpp::export(name = ".liberation_prediction_tape_create")]]
SEXP liberation_prediction_tape_create(
    SEXP engine_pointer, SEXP data_input,
    const Rcpp::NumericVector& theta, const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma) {
  Rcpp::XPtr<liberation::ModelEngine> engine(engine_pointer);
  const liberation::EventDataView data = liberation::event_data_view(data_input);
  liberation::require_materialized_addl(data);
  std::unique_ptr<liberation::PredictionTape> tape = liberation::record_prediction_tape(
    *engine, data, theta, eta, sigma);
  Rcpp::XPtr<liberation::PredictionTape> pointer(tape.release(), true);
  pointer.attr("class") = Rcpp::CharacterVector::create(
    "liberation_prediction_tape_ptr", "externalptr");
  pointer.attr("domain") = Rcpp::wrap(pointer->domain_names);
  pointer.attr("dynamic_columns") = Rcpp::wrap(pointer->dynamic_columns);
  pointer.attr("dynamic_parameters") =
    static_cast<double>(pointer->fun.size_dyn_ind());
  pointer.attr("propagation_kernel") = pointer->propagation_kernel;
  pointer.attr("operation_count") = static_cast<double>(pointer->operation_count);
  pointer.attr("variable_count") = static_cast<double>(pointer->variable_count);
  return pointer;
}

// [[Rcpp::export(name = ".liberation_prediction_tape_info")]]
Rcpp::List liberation_prediction_tape_info(SEXP tape_pointer) {
  Rcpp::XPtr<liberation::PredictionTape> tape(tape_pointer);
  const std::size_t taylor_bytes =
    tape->fun.size_var() * tape->fun.size_order() *
    std::max<std::size_t>(tape->fun.size_direction(), 1U) * sizeof(double);
  const std::size_t resident_proxy = tape->fun.size_op_seq() +
    tape->fun.size_random() + tape->fun.size_forward_bool() +
    tape->fun.size_forward_set() + taylor_bytes;
  return Rcpp::List::create(
    Rcpp::Named("operations") = static_cast<double>(tape->fun.size_op()),
    Rcpp::Named("operator_arguments") =
      static_cast<double>(tape->fun.size_op_arg()),
    Rcpp::Named("variables") = static_cast<double>(tape->fun.size_var()),
    Rcpp::Named("parameters") = static_cast<double>(tape->fun.size_par()),
    Rcpp::Named("dynamic_independent") =
      static_cast<double>(tape->fun.size_dyn_ind()),
    Rcpp::Named("dynamic_parameters") =
      static_cast<double>(tape->fun.size_dyn_par()),
    Rcpp::Named("dynamic_arguments") =
      static_cast<double>(tape->fun.size_dyn_arg()),
    Rcpp::Named("taylor_orders") =
      static_cast<double>(tape->fun.size_order()),
    Rcpp::Named("taylor_directions") =
      static_cast<double>(tape->fun.size_direction()),
    Rcpp::Named("operation_sequence_bytes") =
      static_cast<double>(tape->fun.size_op_seq()),
    Rcpp::Named("random_access_bytes") =
      static_cast<double>(tape->fun.size_random()),
    Rcpp::Named("forward_sparsity_bytes") = static_cast<double>(
      tape->fun.size_forward_bool() + tape->fun.size_forward_set()),
    Rcpp::Named("taylor_bytes_proxy") = static_cast<double>(taylor_bytes),
    Rcpp::Named("resident_bytes_proxy") =
      static_cast<double>(resident_proxy),
    Rcpp::Named("propagation_kernel") = tape->propagation_kernel,
    Rcpp::Named("derivative_strategy") = tape->derivative_strategy,
    Rcpp::Named("jacobian_nonzeros") =
      static_cast<double>(tape->jacobian_nonzeros)
  );
}

// [[Rcpp::export(name = ".liberation_prediction_tape_new_dynamic")]]
Rcpp::NumericVector liberation_prediction_tape_new_dynamic(
    SEXP tape_pointer, SEXP data_input) {
  const liberation::EventDataView data = liberation::event_data_view(data_input);
  liberation::require_materialized_addl(data);
  Rcpp::XPtr<liberation::PredictionTape> tape(tape_pointer);
  std::vector<double> values = liberation::prediction_dynamic_values(
    tape->dynamic_columns, data, tape->n_rows);
  tape->fun.new_dynamic(values);
  tape->dynamic_values = values;
  Rcpp::NumericVector result(values.begin(), values.end());
  result.attr("columns") = Rcpp::wrap(tape->dynamic_columns);
  return result;
}

// [[Rcpp::export(name = ".liberation_prediction_tape_set_dynamic")]]
void liberation_prediction_tape_set_dynamic(
    SEXP tape_pointer, const Rcpp::NumericVector& values) {
  Rcpp::XPtr<liberation::PredictionTape> tape(tape_pointer);
  liberation::set_tape_dynamic_values(
    *tape, Rcpp::as<std::vector<double>>(values), "Prediction tape");
}

// [[Rcpp::export(name = ".liberation_fo_tape_new_dynamic")]]
Rcpp::NumericVector liberation_fo_tape_new_dynamic(
    SEXP tape_pointer, SEXP data_input) {
  const liberation::EventDataView data = liberation::event_data_view(data_input);
  liberation::require_materialized_addl(data);
  Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointer);
  liberation::set_fo_dynamic(*tape, data);
  return Rcpp::wrap(tape->dynamic_values);
}

// [[Rcpp::export(name = ".liberation_prediction_tape_eval")]]
Rcpp::List liberation_prediction_tape_eval(
    SEXP tape_pointer, const Rcpp::NumericVector& point, bool jacobian = true) {
  Rcpp::XPtr<liberation::PredictionTape> tape(tape_pointer);
  std::vector<double> x = liberation::prediction_point(*tape, point);
  std::ostringstream messages;
  std::vector<double> value = tape->fun.Forward(0, x, messages);
  liberation::require_unchanged_path(tape->fun, "prediction evaluation");
  Rcpp::List result = Rcpp::List::create(Rcpp::Named("value") = Rcpp::wrap(value));
  result.attr("domain") = Rcpp::wrap(tape->domain_names);
  if (jacobian) {
    const std::size_t n = tape->domain_names.size();
    const std::size_t m = static_cast<std::size_t>(tape->n_rows);
    Rcpp::NumericMatrix derivative(m, n);
    std::size_t nonzeros = 0U;
    if (m * n >= 4096U && m >= 32U) {
      CppAD::vectorBool select_domain(n), select_range(m);
      for (std::size_t column = 0; column < n; ++column) select_domain[column] = true;
      for (std::size_t row = 0; row < m; ++row) select_range[row] = true;
      using SizeVector = CppAD::vector<std::size_t>;
      using BaseVector = CppAD::vector<double>;
      CppAD::sparse_rcv<SizeVector, BaseVector> sparse;
      BaseVector sparse_point(x.size());
      for (std::size_t index = 0; index < x.size(); ++index) sparse_point[index] = x[index];
      tape->fun.subgraph_jac_rev(
        select_domain, select_range, sparse_point, sparse);
      liberation::require_unchanged_path(
        tape->fun, "sparse prediction evaluation");
      for (std::size_t index = 0; index < sparse.nnz(); ++index) {
        derivative(sparse.row()[index], sparse.col()[index]) = sparse.val()[index];
      }
      nonzeros = sparse.nnz();
      tape->derivative_strategy = "subgraph-reverse";
    } else {
      constexpr std::size_t block_max = 16U;
      for (std::size_t first = 0; first < n; first += block_max) {
        const std::size_t directions = std::min(block_max, n - first);
        std::vector<double> seed(n * directions, 0.0);
        for (std::size_t direction = 0; direction < directions; ++direction) {
          seed[(first + direction) * directions + direction] = 1.0;
        }
        const std::vector<double> forward = directions == 1U ?
          tape->fun.Forward(1, seed) :
          tape->fun.Forward(1, directions, seed);
        for (std::size_t row = 0; row < m; ++row) {
          for (std::size_t direction = 0; direction < directions; ++direction) {
            const double current = forward[row * directions + direction];
            derivative(row, first + direction) = current;
            if (current != 0.0) ++nonzeros;
          }
        }
      }
      tape->derivative_strategy = n == 1U ? "forward" : "multi-forward";
    }
    tape->jacobian_nonzeros = nonzeros;
    derivative.attr("dimnames") = Rcpp::List::create(R_NilValue, Rcpp::wrap(tape->domain_names));
    result["jacobian"] = derivative;
  }
  result.attr("derivative_strategy") = tape->derivative_strategy;
  result.attr("jacobian_nonzeros") = static_cast<double>(tape->jacobian_nonzeros);
  return result;
}

// [[Rcpp::export(name = ".liberation_prediction_tape_eval_subset")]]
Rcpp::List liberation_prediction_tape_eval_subset(
    SEXP tape_pointer, const Rcpp::NumericVector& point,
    const Rcpp::IntegerVector& columns) {
  Rcpp::XPtr<liberation::PredictionTape> tape(tape_pointer);
  std::vector<double> x = liberation::prediction_point(*tape, point);
  std::ostringstream messages;
  const std::vector<double> value = tape->fun.Forward(0, x, messages);
  liberation::require_unchanged_path(tape->fun, "prediction subset evaluation");
  const std::size_t domain = tape->domain_names.size();
  const std::size_t range = static_cast<std::size_t>(tape->n_rows);
  Rcpp::NumericMatrix derivative(range, columns.size());
  Rcpp::CharacterVector names(columns.size());
  std::vector<std::size_t> selected_columns(static_cast<std::size_t>(columns.size()));
  for (R_xlen_t selected = 0; selected < columns.size(); ++selected) {
    const int column = columns[selected] - 1;
    if (column < 0 || static_cast<std::size_t>(column) >= domain) {
      Rcpp::stop("Prediction derivative column is outside the tape domain.");
    }
    selected_columns[static_cast<std::size_t>(selected)] =
      static_cast<std::size_t>(column);
    names[selected] = tape->domain_names[static_cast<std::size_t>(column)];
  }
  constexpr std::size_t block_max = 16U;
  for (std::size_t first = 0; first < selected_columns.size(); first += block_max) {
    const std::size_t directions = std::min(block_max, selected_columns.size() - first);
    std::vector<double> seed(domain * directions, 0.0);
    for (std::size_t direction = 0; direction < directions; ++direction) {
      seed[selected_columns[first + direction] * directions + direction] = 1.0;
    }
    const std::vector<double> forward = directions == 1U ?
      tape->fun.Forward(1, seed) :
      tape->fun.Forward(1, directions, seed);
    for (std::size_t row = 0; row < range; ++row) {
      for (std::size_t direction = 0; direction < directions; ++direction) {
        derivative(static_cast<int>(row), static_cast<int>(first + direction)) =
          forward[row * directions + direction];
      }
    }
  }
  tape->derivative_strategy = selected_columns.size() <= 1U ?
    "forward-subset" : "multi-forward-subset";
  derivative.attr("dimnames") = Rcpp::List::create(R_NilValue, names);
  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("value") = Rcpp::wrap(value),
    Rcpp::Named("jacobian") = derivative
  );
  result.attr("domain") = names;
  return result;
}

// [[Rcpp::export(name = ".liberation_matrix_exp_pade")]]
Rcpp::NumericMatrix liberation_matrix_exp_pade(const Rcpp::NumericMatrix& matrix,
                                                double dt = 1.0) {
  const auto mapped = libertad::r_matrix_map(matrix);
  return libertad::eigen_matrix_to_r(
    liberation::matrix_exp_pade(Eigen::MatrixXd(mapped * dt)));
}

// [[Rcpp::export(name = ".liberation_fo_tape_create")]]
SEXP liberation_fo_tape_create(
    SEXP engine_pointer, SEXP prediction_tape_pointer,
    SEXP data_input, const Rcpp::NumericVector& theta,
    const Rcpp::NumericVector& sigma, const Rcpp::NumericVector& omega,
    bool low_rank = false, double low_rank_tolerance = 1e-9,
    double low_rank_condition_tolerance = 1e-12) {
  const liberation::EventDataView data = liberation::event_data_view(data_input);
  liberation::require_materialized_addl(data);
  Rcpp::XPtr<liberation::ModelEngine> engine(engine_pointer);
  Rcpp::XPtr<liberation::PredictionTape> prediction_tape(prediction_tape_pointer);
  std::unique_ptr<liberation::ObjectiveTape> tape = liberation::record_fo_tape(
    *engine, *prediction_tape, data, theta, sigma, omega, low_rank,
    low_rank_tolerance, low_rank_condition_tolerance);
  Rcpp::XPtr<liberation::ObjectiveTape> pointer(tape.release(), true);
  pointer.attr("class") = Rcpp::CharacterVector::create(
    "liberation_fo_tape_ptr", "liberation_objective_tape_ptr", "externalptr");
  pointer.attr("domain") = Rcpp::wrap(pointer->domain_names);
  return pointer;
}

// [[Rcpp::export(name = ".liberation_curvature_tape_create")]]
SEXP liberation_curvature_tape_create(
    SEXP engine_pointer, SEXP prediction_tape_pointer,
    SEXP objective_tape_pointer, SEXP data_input,
    const Rcpp::NumericVector& theta, const Rcpp::NumericVector& eta,
    const Rcpp::NumericVector& sigma, const Rcpp::NumericVector& omega,
    const std::string& approximation) {
  const liberation::EventDataView data = liberation::event_data_view(data_input);
  liberation::require_materialized_addl(data);
  Rcpp::XPtr<liberation::ModelEngine> engine(engine_pointer);
  Rcpp::XPtr<liberation::PredictionTape> prediction_tape(prediction_tape_pointer);
  Rcpp::XPtr<liberation::ObjectiveTape> objective_tape(objective_tape_pointer);
  std::unique_ptr<liberation::ObjectiveTape> tape = liberation::record_curvature_tape(
    *engine, *prediction_tape, *objective_tape, data,
    theta, eta, sigma, omega, approximation);
  Rcpp::XPtr<liberation::ObjectiveTape> pointer(tape.release(), true);
  pointer.attr("class") = Rcpp::CharacterVector::create(
    "liberation_curvature_tape_ptr", "liberation_objective_tape_ptr", "externalptr");
  pointer.attr("domain") = Rcpp::wrap(pointer->domain_names);
  return pointer;
}

// [[Rcpp::export(name = ".liberation_objective_tape_create")]]
SEXP liberation_objective_tape_create(
    SEXP engine_pointer, SEXP data_input,
    const Rcpp::NumericVector& theta, const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma, const Rcpp::NumericVector& omega,
    bool interaction = true) {
  const liberation::EventDataView data = liberation::event_data_view(data_input);
  liberation::require_materialized_addl(data);
  Rcpp::XPtr<liberation::ModelEngine> engine(engine_pointer);
  std::unique_ptr<liberation::ObjectiveTape> tape = liberation::record_objective_tape(
    *engine, data, theta, eta, sigma, omega, interaction);
  Rcpp::XPtr<liberation::ObjectiveTape> pointer(tape.release(), true);
  pointer.attr("class") = Rcpp::CharacterVector::create(
    "liberation_objective_tape_ptr", "externalptr");
  pointer.attr("domain") = Rcpp::wrap(pointer->domain_names);
  return pointer;
}

// [[Rcpp::export(name = ".liberation_shared_fo_objective_tape_create")]]
SEXP liberation_shared_fo_objective_tape_create(
    SEXP engine_pointer, SEXP prediction_tape_pointer,
    SEXP data_input, const Rcpp::NumericVector& theta,
    const Rcpp::NumericMatrix& eta, const Rcpp::NumericVector& sigma,
    const Rcpp::NumericVector& omega) {
  const liberation::EventDataView data = liberation::event_data_view(data_input);
  liberation::require_materialized_addl(data);
  Rcpp::XPtr<liberation::ModelEngine> engine(engine_pointer);
  Rcpp::XPtr<liberation::PredictionTape> prediction(prediction_tape_pointer);
  std::unique_ptr<liberation::ObjectiveTape> tape =
    liberation::record_shared_fo_objective_tape(
      *engine, *prediction, data, theta, eta, sigma, omega);
  Rcpp::XPtr<liberation::ObjectiveTape> pointer(tape.release(), true);
  pointer.attr("class") = Rcpp::CharacterVector::create(
    "liberation_objective_tape_ptr", "externalptr");
  pointer.attr("domain") = Rcpp::wrap(pointer->domain_names);
  return pointer;
}

// [[Rcpp::export(name = ".liberation_objective_tape_new_dynamic")]]
Rcpp::NumericVector liberation_objective_tape_new_dynamic(
    SEXP tape_pointer, SEXP data_input) {
  const liberation::EventDataView data = liberation::event_data_view(data_input);
  Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointer);
  liberation::set_shared_objective_dynamic(*tape, data);
  Rcpp::NumericVector result = Rcpp::wrap(tape->dynamic_values);
  result.attr("columns") = Rcpp::wrap(tape->dynamic_columns);
  return result;
}

// [[Rcpp::export(name = ".liberation_objective_tape_set_dynamic")]]
void liberation_objective_tape_set_dynamic(
    SEXP tape_pointer, const Rcpp::NumericVector& values) {
  Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointer);
  liberation::set_tape_dynamic_values(
    *tape, Rcpp::as<std::vector<double>>(values), "Objective tape");
}

// [[Rcpp::export(name = ".liberation_objective_tape_eval")]]
Rcpp::List liberation_objective_tape_eval(
    SEXP tape_pointer, const Rcpp::NumericVector& point,
    bool gradient = true, bool hessian = false) {
  Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointer);
  if (point.size() != static_cast<R_xlen_t>(tape->domain_names.size())) {
    Rcpp::stop("Objective tape point has the wrong length.");
  }
  std::vector<double> x = Rcpp::as<std::vector<double>>(point);
  std::ostringstream messages;
  std::vector<double> value = tape->fun.Forward(0, x, messages);
  liberation::require_unchanged_path(tape->fun, "objective evaluation");
  Rcpp::List result = Rcpp::List::create(Rcpp::Named("value") = value[0]);
  if (gradient || hessian) {
    std::vector<double> weight(1, 1.0);
    std::vector<double> derivative = tape->fun.Reverse(1, weight);
    Rcpp::NumericVector output(derivative.begin(), derivative.end());
    output.attr("names") = Rcpp::wrap(tape->domain_names);
    result["gradient"] = output;
  }
  if (hessian) {
    const std::size_t n = tape->domain_names.size();
    Rcpp::NumericMatrix output(n, n);
    libertad::analyse_hessian_sparsity(
      tape->fun, tape->hessian_cache);
    if (tape->hessian_cache.use_sparse) {
      const std::vector<double> values = libertad::sparse_hessian(
        tape->fun, x, tape->hessian_cache);
      liberation::require_unchanged_path(
        tape->fun, "sparse objective Hessian evaluation");
      for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t column = 0; column < n; ++column) {
          output(row, column) = values[row * n + column];
        }
      }
    } else {
      std::vector<double> direction(n, 0.0);
      std::vector<double> weight(1, 1.0);
      for (std::size_t column = 0; column < n; ++column) {
        direction[column] = 1.0;
        tape->fun.Forward(1, direction, messages);
        direction[column] = 0.0;
        std::vector<double> reverse = tape->fun.Reverse(2, weight);
        for (std::size_t row = 0; row < n; ++row) {
          output(row, column) = reverse[row * 2 + 1];
        }
      }
    }
    output.attr("dimnames") = Rcpp::List::create(
      Rcpp::wrap(tape->domain_names), Rcpp::wrap(tape->domain_names));
    result["hessian"] = output;
  }
  result.attr("domain") = Rcpp::wrap(tape->domain_names);
  return result;
}

// [[Rcpp::export(name = ".liberation_objective_tape_info")]]
Rcpp::List liberation_objective_tape_info(SEXP tape_pointer) {
  Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointer);
  const std::size_t taylor_bytes =
    tape->fun.size_var() * tape->fun.size_order() *
    std::max<std::size_t>(tape->fun.size_direction(), 1U) * sizeof(double);
  return Rcpp::List::create(
    Rcpp::Named("operations") = static_cast<double>(tape->fun.size_op()),
    Rcpp::Named("operator_arguments") =
      static_cast<double>(tape->fun.size_op_arg()),
    Rcpp::Named("variables") = static_cast<double>(tape->fun.size_var()),
    Rcpp::Named("parameters") = static_cast<double>(tape->fun.size_par()),
    Rcpp::Named("dynamic_independent") =
      static_cast<double>(tape->fun.size_dyn_ind()),
    Rcpp::Named("dynamic_parameters") =
      static_cast<double>(tape->fun.size_dyn_par()),
    Rcpp::Named("operation_sequence_bytes") =
      static_cast<double>(tape->fun.size_op_seq()),
    Rcpp::Named("random_access_bytes") =
      static_cast<double>(tape->fun.size_random()),
    Rcpp::Named("forward_sparsity_bytes") = static_cast<double>(
      tape->fun.size_forward_bool() + tape->fun.size_forward_set()),
    Rcpp::Named("taylor_bytes_proxy") = static_cast<double>(taylor_bytes),
    Rcpp::Named("hessian_strategy") = tape->hessian_cache.strategy,
    Rcpp::Named("hessian_nonzeros") =
      static_cast<double>(tape->hessian_cache.nonzeros),
    Rcpp::Named("hessian_density") = tape->hessian_cache.density,
    Rcpp::Named("hessian_sweeps") =
      static_cast<double>(tape->hessian_cache.sweeps),
    Rcpp::Named("resident_bytes_proxy") = static_cast<double>(
      tape->fun.size_op_seq() + tape->fun.size_random() +
      tape->fun.size_forward_bool() + tape->fun.size_forward_set() +
      taylor_bytes)
  );
}

// [[Rcpp::export(name = ".liberation_hmc_target_eval")]]
Rcpp::List liberation_hmc_target_eval(
    SEXP tape_pointer, const Rcpp::NumericVector& q,
    const Rcpp::List& config) {
  Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointer);
  return liberation::native_hmc_target_eval(*tape, q, config);
}

// [[Rcpp::export(name = ".liberation_hmc_sample")]]
Rcpp::List liberation_hmc_sample(
    SEXP tape_pointer, const Rcpp::List& config, const std::string& method,
    int n_warmup, int n_sample, int n_thin, int n_chains, double seed,
    double step_size, double target_acceptance, bool adapt_mass,
    int n_leapfrog, int max_depth, double divergence_threshold,
    int print_every) {
  if (!std::isfinite(seed) || seed < 0.0) {
    Rcpp::stop("Native HMC seed must be a non-negative finite number.");
  }
  Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointer);
  return liberation::native_hmc_sample(
    *tape, config, method, n_warmup, n_sample, n_thin, n_chains,
    static_cast<std::uint64_t>(seed), step_size, target_acceptance,
    adapt_mass, n_leapfrog, max_depth, divergence_threshold, print_every
  );
}

// [[Rcpp::export(name = ".liberation_objective_tape_eta_values")]]
Rcpp::NumericVector liberation_objective_tape_eta_values(
    SEXP tape_pointer, const Rcpp::NumericVector& point,
    const Rcpp::IntegerVector& eta_positions,
    const Rcpp::NumericMatrix& eta) {
  Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointer);
  const std::size_t domain = tape->domain_names.size();
  if (point.size() != static_cast<R_xlen_t>(domain)) {
    Rcpp::stop("Objective tape point has the wrong length.");
  }
  if (eta.ncol() != eta_positions.size()) {
    Rcpp::stop("ETA samples have the wrong number of columns.");
  }
  std::vector<std::size_t> positions;
  positions.reserve(static_cast<std::size_t>(eta_positions.size()));
  for (int value : eta_positions) {
    if (value < 1 || static_cast<std::size_t>(value) > domain) {
      Rcpp::stop("ETA position is outside the objective tape domain.");
    }
    positions.push_back(static_cast<std::size_t>(value - 1));
  }
  std::vector<double> x = Rcpp::as<std::vector<double>>(point);
  Rcpp::NumericVector values(eta.nrow());
  std::ostringstream messages;
  for (int sample = 0; sample < eta.nrow(); ++sample) {
    for (int column = 0; column < eta.ncol(); ++column) {
      x[positions[static_cast<std::size_t>(column)]] = eta(sample, column);
    }
    const std::vector<double> value = tape->fun.Forward(0, x, messages);
    liberation::require_unchanged_path(tape->fun, "objective ETA batch");
    values[sample] = value.empty() ? NA_REAL : value[0];
    if ((sample + 1) % 256 == 0) Rcpp::checkUserInterrupt();
  }
  return values;
}

// [[Rcpp::export(name = ".liberation_objective_tape_collection_values")]]
Rcpp::NumericVector liberation_objective_tape_collection_values(
    const Rcpp::List& tape_pointers, const Rcpp::NumericMatrix& points,
    const Rcpp::List& subject_data) {
  if (points.nrow() != tape_pointers.size() ||
      subject_data.size() != tape_pointers.size()) {
    Rcpp::stop("Objective point rows must match the number of tapes.");
  }
  Rcpp::NumericVector values(points.nrow());
  for (int row = 0; row < points.nrow(); ++row) {
    Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointers[row]);
    liberation::set_objective_dynamic_input(*tape, subject_data[row]);
    if (points.ncol() != static_cast<int>(tape->domain_names.size())) {
      Rcpp::stop("An objective point has the wrong length.");
    }
    std::vector<double> point(static_cast<std::size_t>(points.ncol()));
    for (int column = 0; column < points.ncol(); ++column) {
      point[static_cast<std::size_t>(column)] = points(row, column);
    }
    std::ostringstream messages;
    const std::vector<double> value = tape->fun.Forward(0, point, messages);
    liberation::require_unchanged_path(tape->fun, "objective collection");
    values[row] = value.empty() ? NA_REAL : value[0];
    if ((row + 1) % 256 == 0) Rcpp::checkUserInterrupt();
  }
  return values;
}

// [[Rcpp::export(name = ".liberation_objective_tape_collection_gradients")]]
Rcpp::NumericMatrix liberation_objective_tape_collection_gradients(
    const Rcpp::List& tape_pointers, const Rcpp::NumericMatrix& points,
    const Rcpp::List& subject_data) {
  if (points.nrow() != tape_pointers.size() ||
      subject_data.size() != tape_pointers.size()) {
    Rcpp::stop("Objective point rows must match the number of tapes.");
  }
  if (!points.nrow()) return Rcpp::NumericMatrix(0, points.ncol());
  Rcpp::NumericMatrix gradients(points.nrow(), points.ncol());
  const std::vector<double> weight(1, 1.0);
  for (int row = 0; row < points.nrow(); ++row) {
    Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointers[row]);
    liberation::set_objective_dynamic_input(*tape, subject_data[row]);
    if (points.ncol() != static_cast<int>(tape->domain_names.size())) {
      Rcpp::stop("An objective point has the wrong length.");
    }
    std::vector<double> point(static_cast<std::size_t>(points.ncol()));
    for (int column = 0; column < points.ncol(); ++column) {
      point[static_cast<std::size_t>(column)] = points(row, column);
    }
    std::ostringstream messages;
    tape->fun.Forward(0, point, messages);
    liberation::require_unchanged_path(tape->fun, "objective gradient collection");
    const std::vector<double> derivative = tape->fun.Reverse(1, weight);
    for (int column = 0; column < points.ncol(); ++column) {
      gradients(row, column) = derivative[static_cast<std::size_t>(column)];
    }
    if ((row + 1) % 256 == 0) Rcpp::checkUserInterrupt();
  }
  return gradients;
}

// Evaluate fixed-ETA subject objectives and exact gradients in one pass.
// R's L-BFGS-B driver normally requests fn and gr consecutively at the same
// point; returning both avoids a duplicate Forward(0) replay while retaining
// R-side subject-order aggregation in compatibility-sensitive estimators.
// [[Rcpp::export(name = ".liberation_objective_tape_collection_value_gradients")]]
Rcpp::List liberation_objective_tape_collection_value_gradients(
    const Rcpp::List& tape_pointers, const Rcpp::NumericMatrix& points,
    const Rcpp::List& subject_data) {
  if (points.nrow() != tape_pointers.size() ||
      subject_data.size() != tape_pointers.size()) {
    Rcpp::stop("Objective point rows must match the number of tapes.");
  }
  Rcpp::NumericVector values(points.nrow());
  Rcpp::NumericMatrix gradients(points.nrow(), points.ncol());
  const std::vector<double> weight(1, 1.0);
  for (int row = 0; row < points.nrow(); ++row) {
    Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointers[row]);
    liberation::set_objective_dynamic_input(*tape, subject_data[row]);
    if (points.ncol() != static_cast<int>(tape->domain_names.size())) {
      Rcpp::stop("An objective point has the wrong length.");
    }
    std::vector<double> point(static_cast<std::size_t>(points.ncol()));
    for (int column = 0; column < points.ncol(); ++column) {
      point[static_cast<std::size_t>(column)] = points(row, column);
    }
    std::ostringstream messages;
    const std::vector<double> value = tape->fun.Forward(0, point, messages);
    liberation::require_unchanged_path(
      tape->fun, "objective value-gradient collection");
    values[row] = value.empty() ? NA_REAL : value[0];
    const std::vector<double> derivative = tape->fun.Reverse(1, weight);
    liberation::require_unchanged_path(
      tape->fun, "objective value-gradient collection");
    for (int column = 0; column < points.ncol(); ++column) {
      gradients(row, column) = derivative[static_cast<std::size_t>(column)];
    }
    if ((row + 1) % 256 == 0) Rcpp::checkUserInterrupt();
  }
  return Rcpp::List::create(
    Rcpp::Named("value") = values,
    Rcpp::Named("gradient") = gradients
  );
}

// [[Rcpp::export(name = ".liberation_objective_tape_hessian_subset")]]
Rcpp::NumericMatrix liberation_objective_tape_hessian_subset(
    SEXP tape_pointer, const Rcpp::NumericVector& point,
    const Rcpp::IntegerVector& row_positions,
    const Rcpp::IntegerVector& column_positions) {
  Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointer);
  const std::size_t domain = tape->domain_names.size();
  if (point.size() != static_cast<R_xlen_t>(domain)) {
    Rcpp::stop("Objective tape point has the wrong length.");
  }
  auto positions = [domain](const Rcpp::IntegerVector& source) {
    std::vector<std::size_t> result;
    result.reserve(static_cast<std::size_t>(source.size()));
    for (int value : source) {
      if (value < 1 || static_cast<std::size_t>(value) > domain) {
        Rcpp::stop("Hessian position is outside the objective tape domain.");
      }
      result.push_back(static_cast<std::size_t>(value - 1));
    }
    return result;
  };
  const std::vector<std::size_t> rows = positions(row_positions);
  const std::vector<std::size_t> columns = positions(column_positions);
  Rcpp::NumericMatrix result(rows.size(), columns.size());
  std::vector<double> x = Rcpp::as<std::vector<double>>(point);
  std::ostringstream messages;
  tape->fun.Forward(0, x, messages);
  liberation::require_unchanged_path(tape->fun, "objective Hessian subset");
  const std::vector<double> weight(1, 1.0);
  std::vector<double> direction(domain, 0.0);
  for (std::size_t column = 0; column < columns.size(); ++column) {
    direction[columns[column]] = 1.0;
    tape->fun.Forward(1, direction, messages);
    direction[columns[column]] = 0.0;
    const std::vector<double> reverse = tape->fun.Reverse(2, weight);
    for (std::size_t row = 0; row < rows.size(); ++row) {
      result(static_cast<int>(row), static_cast<int>(column)) =
        reverse[rows[row] * 2U + 1U];
    }
  }
  return result;
}

// [[Rcpp::export(name = ".liberation_nested_population_gradient")]]
Rcpp::List liberation_nested_population_gradient(
    const Rcpp::List& objective_tapes, const Rcpp::List& curvature_tapes,
    const Rcpp::NumericMatrix& points,
    const Rcpp::IntegerVector& eta_positions,
    const Rcpp::IntegerVector& population_positions,
    const Rcpp::NumericMatrix& transform) {
  const int subjects = points.nrow();
  const int n_eta = eta_positions.size();
  const int n_population = population_positions.size();
  const int n_outer = transform.ncol();
  if (objective_tapes.size() != subjects || curvature_tapes.size() != subjects ||
      transform.nrow() != n_population) {
    Rcpp::stop("Nested-gradient batch dimensions are inconsistent.");
  }
  Rcpp::NumericMatrix subject_gradients(subjects, n_outer);
  Rcpp::NumericVector jitters(subjects);
  const std::vector<double> weight(1, 1.0);
  for (int subject = 0; subject < subjects; ++subject) {
    Rcpp::XPtr<liberation::ObjectiveTape> objective(objective_tapes[subject]);
    Rcpp::XPtr<liberation::ObjectiveTape> curvature(curvature_tapes[subject]);
    const std::size_t domain = objective->domain_names.size();
    if (points.ncol() != static_cast<int>(domain) ||
        curvature->domain_names.size() != domain) {
      Rcpp::stop("A nested-gradient objective point has the wrong length.");
    }
    std::vector<std::size_t> eta, population;
    for (int value : eta_positions) {
      if (value < 1 || static_cast<std::size_t>(value) > domain) {
        Rcpp::stop("ETA position is outside a nested-gradient tape domain.");
      }
      eta.push_back(static_cast<std::size_t>(value - 1));
    }
    for (int value : population_positions) {
      if (value < 1 || static_cast<std::size_t>(value) > domain) {
        Rcpp::stop("Population position is outside a nested-gradient tape domain.");
      }
      population.push_back(static_cast<std::size_t>(value - 1));
    }
    std::vector<double> point(domain);
    for (int column = 0; column < points.ncol(); ++column) {
      point[static_cast<std::size_t>(column)] = points(subject, column);
    }
    std::ostringstream messages;
    objective->fun.Forward(0, point, messages);
    const std::vector<double> objective_derivative = objective->fun.Reverse(1, weight);
    Eigen::MatrixXd mixed(n_eta, n_eta + n_population);
    std::vector<double> direction(domain, 0.0);
    for (int column = 0; column < n_eta + n_population; ++column) {
      const std::size_t position = column < n_eta ?
        eta[static_cast<std::size_t>(column)] :
        population[static_cast<std::size_t>(column - n_eta)];
      direction[position] = 1.0;
      objective->fun.Forward(1, direction, messages);
      direction[position] = 0.0;
      const std::vector<double> reverse = objective->fun.Reverse(2, weight);
      for (int row = 0; row < n_eta; ++row) {
        mixed(row, column) = reverse[eta[static_cast<std::size_t>(row)] * 2U + 1U];
      }
    }
    Eigen::MatrixXd eta_hessian;
    if (n_eta) {
      eta_hessian = 0.5 *
        (mixed.leftCols(n_eta) + mixed.leftCols(n_eta).transpose()).eval();
    } else {
      eta_hessian = Eigen::MatrixXd::Zero(0, 0);
    }
    double jitter = 0.0;
    if (n_eta) {
      auto eigen = libertad::detail::self_adjoint_eigen(eta_hessian, false);
      if (eigen.info != Eigen::Success) {
        Rcpp::stop("Conditional ETA curvature eigen decomposition failed.");
      }
      const double largest = std::max(eigen.values.cwiseAbs().maxCoeff(), 1.0);
      jitter = std::max(0.0, largest * 1e-9 - eigen.values.minCoeff());
      if (jitter > largest * 1e-2) {
        Rcpp::stop("Conditional ETA curvature is not sufficiently positive definite.");
      }
      eta_hessian.diagonal().array() += jitter;
    }
    jitters[subject] = jitter;
    Eigen::MatrixXd mapped_transform(n_population, n_outer);
    for (int row = 0; row < n_population; ++row) {
      for (int column = 0; column < n_outer; ++column) {
        mapped_transform(row, column) = transform(row, column);
      }
    }
    Eigen::MatrixXd sensitivity;
    if (n_eta) {
      sensitivity = -eta_hessian.ldlt().solve(
        mixed.rightCols(n_population) * mapped_transform);
    } else {
      sensitivity = Eigen::MatrixXd::Zero(0, n_outer);
    }
    curvature->fun.Forward(0, point, messages);
    const std::vector<double> curvature_derivative = curvature->fun.Reverse(1, weight);
    for (int outer = 0; outer < n_outer; ++outer) {
      double derivative = 0.0;
      for (int native = 0; native < n_population; ++native) {
        const double chain = transform(native, outer);
        derivative += (objective_derivative[population[static_cast<std::size_t>(native)]] +
          curvature_derivative[population[static_cast<std::size_t>(native)]]) * chain;
      }
      for (int effect = 0; effect < n_eta; ++effect) {
        derivative += curvature_derivative[eta[static_cast<std::size_t>(effect)]] *
          sensitivity(effect, outer);
      }
      subject_gradients(subject, outer) = derivative;
    }
    if ((subject + 1) % 64 == 0) Rcpp::checkUserInterrupt();
  }
  Rcpp::NumericVector gradient(n_outer);
  for (int outer = 0; outer < n_outer; ++outer) {
    for (int subject = 0; subject < subjects; ++subject) {
      gradient[outer] += subject_gradients(subject, outer);
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("gradient") = gradient,
    Rcpp::Named("subject_gradients") = subject_gradients,
    Rcpp::Named("eta_jitter") = jitters);
}

// [[Rcpp::export(name = ".liberation_objective_tape_eta_mode")]]
Rcpp::List liberation_objective_tape_eta_mode(
    SEXP tape_pointer, const Rcpp::NumericVector& point,
    const Rcpp::IntegerVector& eta_positions,
    const Rcpp::NumericVector& start, int maxit = 100,
    double tolerance = 1e-7, bool exact_hessian = true) {
  Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointer);
  if (point.size() != static_cast<R_xlen_t>(tape->domain_names.size())) {
    Rcpp::stop("Objective tape point has the wrong length.");
  }
  std::vector<std::size_t> positions;
  positions.reserve(static_cast<std::size_t>(eta_positions.size()));
  for (int value : eta_positions) {
    if (value < 1 || static_cast<std::size_t>(value) > tape->domain_names.size()) {
      Rcpp::stop("ETA position is outside the objective tape domain.");
    }
    positions.push_back(static_cast<std::size_t>(value - 1));
  }
  return liberation::objective_eta_mode(
    *tape, Rcpp::as<std::vector<double>>(point), positions, start,
    maxit, tolerance, exact_hessian);
}

// Native conditional-mode fitting with a non-zero/custom Gaussian prior. The
// recorded objective already contains the population ETA prior; this applies
// the exact quadratic replacement without any R optimizer callbacks.
// [[Rcpp::export(name = ".liberation_objective_tape_eta_mode_prior")]]
Rcpp::List liberation_objective_tape_eta_mode_prior(
    SEXP tape_pointer, const Rcpp::NumericVector& point,
    const Rcpp::IntegerVector& eta_positions,
    const Rcpp::NumericVector& start,
    const Rcpp::NumericVector& prior_mean,
    const Rcpp::NumericMatrix& base_precision_input,
    const Rcpp::NumericMatrix& prior_precision_input,
    int maxit = 100, double tolerance = 1e-7,
    bool exact_hessian = true) {
  Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointer);
  std::vector<std::size_t> positions;
  positions.reserve(static_cast<std::size_t>(eta_positions.size()));
  for (int value : eta_positions) {
    if (value < 1 || value > static_cast<int>(tape->domain_names.size())) {
      Rcpp::stop("ETA position is outside the objective tape domain.");
    }
    positions.push_back(static_cast<std::size_t>(value - 1));
  }
  const int dimension = eta_positions.size();
  if (prior_mean.size() != dimension ||
      base_precision_input.nrow() != dimension ||
      base_precision_input.ncol() != dimension ||
      prior_precision_input.nrow() != dimension ||
      prior_precision_input.ncol() != dimension) {
    Rcpp::stop("Custom ETA prior dimensions do not match the ETA vector.");
  }
  liberation::EtaPriorAdjustment adjustment;
  adjustment.mean = libertad::r_vector_map(prior_mean);
  adjustment.base_precision = libertad::r_matrix_map(base_precision_input);
  adjustment.prior_precision = libertad::r_matrix_map(prior_precision_input);
  if (!adjustment.mean.allFinite() || !adjustment.base_precision.allFinite() ||
      !adjustment.prior_precision.allFinite()) {
    Rcpp::stop("Custom ETA prior inputs must be finite.");
  }
  return liberation::objective_eta_mode(
    *tape, Rcpp::as<std::vector<double>>(point), positions, start,
    maxit, tolerance, exact_hessian, &adjustment);
}

// [[Rcpp::export(name = ".liberation_objective_tape_eta_modes")]]
Rcpp::List liberation_objective_tape_eta_modes(
    const Rcpp::List& tape_pointers, const Rcpp::NumericMatrix& points,
    const Rcpp::IntegerVector& eta_positions,
    const Rcpp::NumericMatrix& starts, int maxit = 100,
    double tolerance = 1e-7, bool exact_hessian = true,
    SEXP subject_data_input = R_NilValue,
    bool reuse_optimizer_state = false) {
  const Rcpp::List subject_data = Rf_isNull(subject_data_input) ?
    Rcpp::List() : Rcpp::as<Rcpp::List>(subject_data_input);
  if (points.nrow() != tape_pointers.size() || starts.nrow() != points.nrow() ||
      (subject_data.size() && subject_data.size() != tape_pointers.size())) {
    Rcpp::stop("ETA-mode rows must match the number of objective tapes.");
  }
  if (starts.ncol() != eta_positions.size()) {
    Rcpp::stop("ETA starting values have the wrong number of columns.");
  }
  Rcpp::List result(points.nrow());
  for (int row = 0; row < points.nrow(); ++row) {
    Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointers[row]);
    if (subject_data.size()) {
      liberation::set_objective_dynamic_input(*tape, subject_data[row]);
    }
    const std::size_t domain = tape->domain_names.size();
    if (points.ncol() != static_cast<int>(domain)) {
      Rcpp::stop("An ETA-mode objective point has the wrong length.");
    }
    std::vector<std::size_t> positions;
    positions.reserve(static_cast<std::size_t>(eta_positions.size()));
    for (int value : eta_positions) {
      if (value < 1 || static_cast<std::size_t>(value) > domain) {
        Rcpp::stop("ETA position is outside an objective tape domain.");
      }
      positions.push_back(static_cast<std::size_t>(value - 1));
    }
    Rcpp::NumericVector start(starts.ncol());
    std::vector<double> point(static_cast<std::size_t>(points.ncol()));
    for (int column = 0; column < points.ncol(); ++column) {
      point[static_cast<std::size_t>(column)] = points(row, column);
    }
    for (int column = 0; column < starts.ncol(); ++column) {
      start[column] = starts(row, column);
    }
    result[row] = liberation::objective_eta_mode(
      *tape, std::move(point), positions, start, maxit, tolerance,
      exact_hessian, nullptr, reuse_optimizer_state);
    if ((row + 1) % 64 == 0) Rcpp::checkUserInterrupt();
  }
  return result;
}

// [[Rcpp::export(name = ".liberation_objective_tape_point_gradients")]]
Rcpp::List liberation_objective_tape_point_gradients(
    SEXP tape_pointer, const Rcpp::NumericMatrix& points) {
  Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointer);
  const std::size_t domain = tape->domain_names.size();
  if (points.ncol() != static_cast<int>(domain)) {
    Rcpp::stop("Objective sample points have the wrong number of columns.");
  }
  Rcpp::NumericVector values(points.nrow());
  Rcpp::NumericMatrix gradients(points.nrow(), points.ncol());
  const std::vector<double> weight(1, 1.0);
  std::ostringstream messages;
  for (int row = 0; row < points.nrow(); ++row) {
    std::vector<double> point(domain);
    for (int column = 0; column < points.ncol(); ++column) {
      point[static_cast<std::size_t>(column)] = points(row, column);
    }
    const std::vector<double> value = tape->fun.Forward(0, point, messages);
    values[row] = value.empty() ? NA_REAL : value[0];
    const std::vector<double> derivative = tape->fun.Reverse(1, weight);
    for (int column = 0; column < points.ncol(); ++column) {
      gradients(row, column) = derivative[static_cast<std::size_t>(column)];
    }
    if ((row + 1) % 256 == 0) Rcpp::checkUserInterrupt();
  }
  gradients.attr("dimnames") = Rcpp::List::create(
    R_NilValue, Rcpp::wrap(tape->domain_names));
  return Rcpp::List::create(
    Rcpp::Named("value") = values,
    Rcpp::Named("gradient") = gradients);
}

// [[Rcpp::export(name = ".liberation_objective_tape_eta_metropolis")]]
Rcpp::List liberation_objective_tape_eta_metropolis(
    const Rcpp::List& tape_pointers, const Rcpp::NumericMatrix& points,
    const Rcpp::IntegerVector& eta_positions,
    const Rcpp::NumericMatrix& current_eta,
    const Rcpp::List& proposal_roots, const Rcpp::NumericMatrix& normals,
    const Rcpp::NumericVector& log_uniforms, int mcmc_steps,
    double step_scale = 0.5,
    Rcpp::Nullable<Rcpp::NumericVector> current_values_input = R_NilValue) {
  const int subjects = points.nrow();
  const int dimension = eta_positions.size();
  if (subjects != tape_pointers.size() || current_eta.nrow() != subjects ||
      current_eta.ncol() != dimension || proposal_roots.size() != subjects ||
      mcmc_steps < 1 || normals.nrow() != subjects * mcmc_steps ||
      normals.ncol() != dimension || log_uniforms.size() != normals.nrow() ||
      !std::isfinite(step_scale) || step_scale <= 0.0) {
    Rcpp::stop("Batched ETA Metropolis inputs are inconsistent.");
  }
  const bool use_current_values = current_values_input.isNotNull();
  Rcpp::NumericVector supplied_values;
  if (use_current_values) {
    supplied_values = Rcpp::NumericVector(current_values_input);
    if (supplied_values.size() != subjects) {
      Rcpp::stop("Cached ETA objective values must contain one value per subject.");
    }
    for (double value : supplied_values) {
      if (!std::isfinite(value)) {
        Rcpp::stop("Cached ETA objective values must be finite.");
      }
    }
  }
  std::vector<std::size_t> positions;
  positions.reserve(static_cast<std::size_t>(dimension));
  for (int value : eta_positions) {
    if (value < 1 || value > points.ncol()) {
      Rcpp::stop("ETA position is outside a Metropolis objective tape domain.");
    }
    positions.push_back(static_cast<std::size_t>(value - 1));
  }
  Rcpp::NumericMatrix eta = Rcpp::clone(current_eta);
  Rcpp::NumericVector values(subjects);
  int accepted = 0;
  int current_evaluations = 0;
  int current_cache_hits = 0;
  int candidate_evaluations = 0;
  std::ostringstream messages;
  for (int subject = 0; subject < subjects; ++subject) {
    Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointers[subject]);
    const std::size_t domain = tape->domain_names.size();
    if (points.ncol() != static_cast<int>(domain)) {
      Rcpp::stop("A batched Metropolis objective point has the wrong length.");
    }
    Rcpp::NumericMatrix root = proposal_roots[subject];
    if (root.nrow() != dimension || root.ncol() != dimension) {
      Rcpp::stop("A Metropolis proposal root has the wrong dimensions.");
    }
    std::vector<double> point(domain);
    for (int column = 0; column < points.ncol(); ++column) {
      point[static_cast<std::size_t>(column)] = points(subject, column);
    }
    for (int column = 0; column < dimension; ++column) {
      point[positions[static_cast<std::size_t>(column)]] = eta(subject, column);
    }
    liberation::EtaEvaluation current;
    if (use_current_values) {
      current.value = supplied_values[subject];
      current.finite = true;
      ++current_cache_hits;
    } else {
      current = liberation::objective_eta_evaluate(
        *tape, point, positions, false, &messages);
      ++current_evaluations;
    }
    if (!current.finite) Rcpp::stop("Current ETA objective is not finite.");
    std::vector<double> candidate_eta(static_cast<std::size_t>(dimension));
    for (int step = 0; step < mcmc_steps; ++step) {
      const int draw = subject * mcmc_steps + step;
      std::vector<double> candidate_point = point;
      for (int row = 0; row < dimension; ++row) {
        double increment = 0.0;
        for (int column = 0; column < dimension; ++column) {
          increment += root(row, column) * normals(draw, column);
        }
        candidate_eta[static_cast<std::size_t>(row)] =
          eta(subject, row) + step_scale * increment;
        candidate_point[positions[static_cast<std::size_t>(row)]] = candidate_eta[row];
      }
      liberation::EtaEvaluation candidate = liberation::objective_eta_evaluate(
        *tape, candidate_point, positions, false, &messages);
      ++candidate_evaluations;
      if (candidate.finite && log_uniforms[draw] <
          -0.5 * (candidate.value - current.value)) {
        for (int row = 0; row < dimension; ++row) {
          eta(subject, row) = candidate_eta[row];
        }
        point.swap(candidate_point);
        current = std::move(candidate);
        ++accepted;
      }
    }
    values[subject] = current.value;
    if ((subject + 1) % 64 == 0) Rcpp::checkUserInterrupt();
  }
  return Rcpp::List::create(
    Rcpp::Named("eta") = eta,
    Rcpp::Named("value") = values,
    Rcpp::Named("accepted") = accepted,
    Rcpp::Named("attempted") = subjects * mcmc_steps,
    Rcpp::Named("current_evaluations") = current_evaluations,
    Rcpp::Named("current_cache_hits") = current_cache_hits,
    Rcpp::Named("candidate_evaluations") = candidate_evaluations);
}

// [[Rcpp::export(name = ".liberation_objective_tape_importance_collection")]]
Rcpp::List liberation_objective_tape_importance_collection(
    const Rcpp::List& tape_pointers, const Rcpp::NumericMatrix& base_points,
    const Rcpp::IntegerVector& eta_positions,
    const Rcpp::List& eta_samples, const Rcpp::List& log_proposals,
    const Rcpp::List& log_measures, const Rcpp::List& measure_signs,
    bool gradient = true, SEXP subject_data_input = R_NilValue) {
  const int subjects = tape_pointers.size();
  const Rcpp::List subject_data = Rf_isNull(subject_data_input) ?
    Rcpp::List() : Rcpp::as<Rcpp::List>(subject_data_input);
  if (base_points.nrow() != subjects || eta_samples.size() != subjects ||
      log_proposals.size() != subjects || log_measures.size() != subjects ||
      measure_signs.size() != subjects ||
      (subject_data.size() && subject_data.size() != subjects)) {
    Rcpp::stop("Importance-collection inputs must have one entry per subject.");
  }
  std::vector<std::size_t> positions;
  positions.reserve(static_cast<std::size_t>(eta_positions.size()));
  for (int value : eta_positions) {
    if (value < 1 || value > base_points.ncol()) {
      Rcpp::stop("An ETA position is outside the objective tape domain.");
    }
    positions.push_back(static_cast<std::size_t>(value - 1));
  }

  Rcpp::List states(subjects);
  Rcpp::NumericVector total_gradient(base_points.ncol());
  if (gradient && subjects > 0) {
    Rcpp::XPtr<liberation::ObjectiveTape> first_tape(tape_pointers[0]);
    total_gradient.attr("names") = Rcpp::wrap(first_tape->domain_names);
  }
  double total_value = 0.0;
  bool total_gradient_valid = gradient;
  const std::vector<double> reverse_weight(1U, 1.0);
  std::ostringstream messages;

  for (int subject = 0; subject < subjects; ++subject) {
    Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointers[subject]);
    if (subject_data.size()) {
      liberation::set_objective_dynamic_input(*tape, subject_data[subject]);
    }
    const std::size_t domain = tape->domain_names.size();
    if (base_points.ncol() != static_cast<int>(domain)) {
      Rcpp::stop("An importance base point has the wrong domain length.");
    }
    Rcpp::NumericMatrix samples(eta_samples[subject]);
    Rcpp::NumericVector log_proposal(log_proposals[subject]);
    Rcpp::NumericVector log_measure(log_measures[subject]);
    Rcpp::NumericVector signs(measure_signs[subject]);
    const int draws = samples.nrow();
    if (samples.ncol() != eta_positions.size() ||
        log_proposal.size() != draws || log_measure.size() != draws ||
        signs.size() != draws || draws < 1) {
      Rcpp::stop("Importance sample dimensions are inconsistent.");
    }
    std::vector<double> point(domain, 0.0);
    for (std::size_t column = 0; column < domain; ++column) {
      point[column] = base_points(subject, static_cast<int>(column));
    }
    std::vector<double> log_integrand(static_cast<std::size_t>(draws));
    std::vector<double> sample_gradient;
    if (gradient) sample_gradient.resize(
      static_cast<std::size_t>(draws) * domain, 0.0);
    double maximum = -std::numeric_limits<double>::infinity();
    for (int draw = 0; draw < draws; ++draw) {
      for (int effect = 0; effect < samples.ncol(); ++effect) {
        point[positions[static_cast<std::size_t>(effect)]] = samples(draw, effect);
      }
      const std::vector<double> value = tape->fun.Forward(0, point, messages);
      liberation::require_unchanged_path(
        tape->fun, "native importance objective");
      const double current = value.empty() ?
        std::numeric_limits<double>::infinity() : value[0];
      const double integrand = -0.5 * current - log_proposal[draw] +
        log_measure[draw];
      log_integrand[static_cast<std::size_t>(draw)] = integrand;
      if (std::isfinite(integrand) && std::isfinite(signs[draw]) &&
          signs[draw] != 0.0) {
        maximum = std::max(maximum, integrand);
      }
      if (gradient) {
        const std::vector<double> derivative =
          tape->fun.Reverse(1, reverse_weight);
        liberation::require_unchanged_path(
          tape->fun, "native importance gradient");
        for (std::size_t column = 0; column < domain; ++column) {
          sample_gradient[static_cast<std::size_t>(draw) * domain + column] =
            derivative[column];
        }
      }
    }

    double signed_total = 0.0;
    double absolute_total = 0.0;
    double squared_absolute = 0.0;
    std::vector<double> scaled(static_cast<std::size_t>(draws), 0.0);
    if (std::isfinite(maximum)) {
      for (int draw = 0; draw < draws; ++draw) {
        const double integrand = log_integrand[static_cast<std::size_t>(draw)];
        if (!std::isfinite(integrand) || !std::isfinite(signs[draw]) ||
            signs[draw] == 0.0) continue;
        const double current = signs[draw] * std::exp(integrand - maximum);
        scaled[static_cast<std::size_t>(draw)] = current;
        signed_total += current;
        absolute_total += std::abs(current);
      }
    }
    const bool valid = std::isfinite(signed_total) &&
      std::isfinite(absolute_total) && signed_total >
        std::numeric_limits<double>::epsilon() * std::max(1.0, absolute_total);
    if (!valid) {
      total_value = std::numeric_limits<double>::infinity();
      total_gradient_valid = false;
      states[subject] = Rcpp::List::create(
        Rcpp::Named("value") = R_PosInf,
        Rcpp::Named("native_gradient") = R_NilValue,
        Rcpp::Named("effective_sample_size") = 0.0,
        Rcpp::Named("cancellation_ratio") = 0.0,
        Rcpp::Named("quadrature_valid") = false);
      continue;
    }
    const double value = -2.0 * (maximum + std::log(signed_total));
    total_value += value;
    Rcpp::NumericVector native_gradient(domain);
    Rcpp::NumericVector posterior_weights(draws);
    for (int draw = 0; draw < draws; ++draw) {
      const double absolute_weight =
        std::abs(scaled[static_cast<std::size_t>(draw)]) / absolute_total;
      squared_absolute += absolute_weight * absolute_weight;
      if (!gradient) continue;
      const double normalized =
        scaled[static_cast<std::size_t>(draw)] / signed_total;
      posterior_weights[draw] = normalized;
      for (std::size_t column = 0; column < domain; ++column) {
        native_gradient[static_cast<R_xlen_t>(column)] += normalized *
          sample_gradient[static_cast<std::size_t>(draw) * domain + column];
      }
    }
    if (!gradient) {
      for (int draw = 0; draw < draws; ++draw) {
        posterior_weights[draw] =
          scaled[static_cast<std::size_t>(draw)] / signed_total;
      }
    }
    if (gradient) {
      for (std::size_t column = 0; column < domain; ++column) {
        total_gradient[static_cast<R_xlen_t>(column)] +=
          native_gradient[static_cast<R_xlen_t>(column)];
      }
      native_gradient.attr("names") = Rcpp::wrap(tape->domain_names);
    }
    states[subject] = Rcpp::List::create(
      Rcpp::Named("value") = value,
      Rcpp::Named("native_gradient") = gradient ?
        Rcpp::RObject(native_gradient) : Rcpp::RObject(R_NilValue),
      Rcpp::Named("weights") = posterior_weights,
      Rcpp::Named("effective_sample_size") =
        squared_absolute > 0.0 ? 1.0 / squared_absolute : 0.0,
      Rcpp::Named("cancellation_ratio") = signed_total / absolute_total,
      Rcpp::Named("quadrature_valid") = true);
    if ((subject + 1) % 32 == 0) Rcpp::checkUserInterrupt();
  }
  return Rcpp::List::create(
    Rcpp::Named("value") = total_value,
    Rcpp::Named("native_gradient") = total_gradient_valid ?
      Rcpp::RObject(total_gradient) : Rcpp::RObject(R_NilValue),
    Rcpp::Named("states") = states);
}

// Evaluate a common support grid against every subject tape without crossing
// the R/C++ boundary once per subject/support combination.  This is used by
// the optimized NPML/NPAG path; the compatibility path keeps the original R
// orchestration for paired reference work.
// [[Rcpp::export(name = ".liberation_objective_tape_eta_grid")]]
Rcpp::List liberation_objective_tape_eta_grid(
    const Rcpp::List& tape_pointers, const Rcpp::NumericMatrix& base_points,
    const Rcpp::IntegerVector& eta_positions,
    const Rcpp::NumericMatrix& eta_grid, bool gradient = false) {
  const int subjects = tape_pointers.size();
  const int supports = eta_grid.nrow();
  if (base_points.nrow() != subjects || supports < 1 ||
      eta_grid.ncol() != eta_positions.size()) {
    Rcpp::stop("ETA-grid dimensions are inconsistent.");
  }
  std::vector<std::size_t> positions;
  positions.reserve(static_cast<std::size_t>(eta_positions.size()));
  for (int value : eta_positions) {
    if (value < 1 || value > base_points.ncol()) {
      Rcpp::stop("An ETA position is outside the objective tape domain.");
    }
    positions.push_back(static_cast<std::size_t>(value - 1));
  }
  Rcpp::NumericMatrix values(subjects, supports);
  Rcpp::List gradients(gradient ? subjects : 0);
  const std::vector<double> reverse_weight(1U, 1.0);
  std::ostringstream messages;
  for (int subject = 0; subject < subjects; ++subject) {
    Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointers[subject]);
    const std::size_t domain = tape->domain_names.size();
    if (base_points.ncol() != static_cast<int>(domain)) {
      Rcpp::stop("An ETA-grid base point has the wrong domain length.");
    }
    std::vector<double> point(domain, 0.0);
    for (std::size_t column = 0; column < domain; ++column) {
      point[column] = base_points(subject, static_cast<int>(column));
    }
    Rcpp::NumericMatrix subject_gradient;
    if (gradient) subject_gradient = Rcpp::NumericMatrix(supports, domain);
    for (int support = 0; support < supports; ++support) {
      for (int effect = 0; effect < eta_grid.ncol(); ++effect) {
        point[positions[static_cast<std::size_t>(effect)]] =
          eta_grid(support, effect);
      }
      const std::vector<double> value = tape->fun.Forward(0, point, messages);
      values(subject, support) = value.empty() ? R_PosInf : value[0];
      if (gradient) {
        const std::vector<double> derivative =
          tape->fun.Reverse(1, reverse_weight);
        for (std::size_t column = 0; column < domain; ++column) {
          subject_gradient(support, static_cast<int>(column)) =
            derivative[column];
        }
      }
    }
    if (gradient) {
      Rcpp::colnames(subject_gradient) = Rcpp::wrap(tape->domain_names);
      gradients[subject] = subject_gradient;
    }
    if ((subject + 1) % 32 == 0) Rcpp::checkUserInterrupt();
  }
  return Rcpp::List::create(
    Rcpp::Named("value") = values,
    Rcpp::Named("gradient") = gradient ?
      Rcpp::RObject(gradients) : Rcpp::RObject(R_NilValue));
}

// [[Rcpp::export(name = ".liberation_np_weights")]]
Rcpp::List liberation_np_weights(
    const Rcpp::NumericMatrix& loglik,
    const Rcpp::NumericVector& initial, int maxit = 1000,
    double tolerance = 1e-8) {
  const int subjects = loglik.nrow();
  const int supports = loglik.ncol();
  if (subjects < 1 || supports < 1 || maxit < 1 ||
      !std::isfinite(tolerance) || tolerance < 0.0) {
    Rcpp::stop("Nonparametric weight inputs are invalid.");
  }
  Rcpp::NumericVector weights(supports);
  if (initial.size() == supports) {
    double total = 0.0;
    for (int support = 0; support < supports; ++support) {
      weights[support] = std::max(initial[support], 0.0);
      total += weights[support];
    }
    if (!(total > 0.0)) total = 1.0;
    for (double& weight : weights) weight /= total;
  } else {
    std::fill(weights.begin(), weights.end(), 1.0 / supports);
  }
  Rcpp::NumericMatrix responsibilities(subjects, supports);
  std::vector<double> history;
  history.reserve(static_cast<std::size_t>(maxit));
  double log_likelihood = -std::numeric_limits<double>::infinity();
  auto expectation = [&]() {
    log_likelihood = 0.0;
    for (int subject = 0; subject < subjects; ++subject) {
      double maximum = -std::numeric_limits<double>::infinity();
      for (int support = 0; support < supports; ++support) {
        const double score = std::log(std::max(
          weights[support], std::numeric_limits<double>::min())) +
          loglik(subject, support);
        responsibilities(subject, support) = score;
        maximum = std::max(maximum, score);
      }
      double total = 0.0;
      for (int support = 0; support < supports; ++support) {
        const double current = std::exp(
          responsibilities(subject, support) - maximum);
        responsibilities(subject, support) = current;
        total += current;
      }
      const double normalizer = maximum + std::log(total);
      log_likelihood += normalizer;
      for (int support = 0; support < supports; ++support) {
        responsibilities(subject, support) /= total;
      }
    }
  };
  int iterations = 0;
  for (; iterations < maxit; ++iterations) {
    expectation();
    Rcpp::NumericVector next(supports);
    for (int support = 0; support < supports; ++support) {
      double total = 0.0;
      for (int subject = 0; subject < subjects; ++subject) {
        total += responsibilities(subject, support);
      }
      next[support] = std::max(
        total / subjects, std::numeric_limits<double>::epsilon());
    }
    double total = std::accumulate(next.begin(), next.end(), 0.0);
    double maximum_change = 0.0;
    for (int support = 0; support < supports; ++support) {
      next[support] /= total;
      maximum_change = std::max(
        maximum_change, std::abs(next[support] - weights[support]));
    }
    history.push_back(log_likelihood);
    weights = next;
    if (maximum_change <= tolerance) {
      ++iterations;
      break;
    }
  }
  expectation();
  return Rcpp::List::create(
    Rcpp::Named("weights") = weights,
    Rcpp::Named("responsibilities") = responsibilities,
    Rcpp::Named("log_likelihood") = log_likelihood,
      Rcpp::Named("iterations") = iterations,
      Rcpp::Named("history") = Rcpp::wrap(history));
}

Rcpp::List liberation_np_responsibilities(
    const Rcpp::NumericMatrix& loglik,
    const Rcpp::NumericVector& weights);

// [[Rcpp::export(name = ".liberation_np_weights_interior")]]
Rcpp::List liberation_np_weights_interior(
    const Rcpp::NumericMatrix& loglik,
    const Rcpp::NumericVector& initial, int maxit = 1000,
    double tolerance = 1e-8) {
  const int subjects = loglik.nrow();
  const int supports = loglik.ncol();
  if (subjects < 1 || supports < 1 || maxit < 1 ||
      !std::isfinite(tolerance) || tolerance <= 0.0) {
    Rcpp::stop("Interior nonparametric weight inputs are invalid.");
  }
  Eigen::MatrixXd likelihood(subjects, supports);
  Eigen::VectorXd row_shift(subjects);
  for (int subject = 0; subject < subjects; ++subject) {
    double maximum = -std::numeric_limits<double>::infinity();
    for (int support = 0; support < supports; ++support) {
      maximum = std::max(maximum, loglik(subject, support));
    }
    if (!std::isfinite(maximum)) {
      Rcpp::stop("Nonparametric log likelihood contains a non-finite row.");
    }
    row_shift[subject] = maximum;
    for (int support = 0; support < supports; ++support) {
      likelihood(subject, support) = std::exp(
        loglik(subject, support) - maximum);
    }
  }
  Eigen::VectorXd weights(supports);
  if (initial.size() == supports) {
    for (int support = 0; support < supports; ++support) {
      weights[support] = std::max(initial[support], 1e-12);
    }
    weights /= weights.sum();
  } else {
    weights.setConstant(1.0 / static_cast<double>(supports));
  }
  double lambda = 0.0;
  double barrier = 0.1;
  int iterations = 0;
  std::vector<double> history;
  history.reserve(static_cast<std::size_t>(maxit));
  const double minimum = std::numeric_limits<double>::min();
  auto objective = [&](const Eigen::VectorXd& value, double mu) {
    if ((value.array() <= 0.0).any()) {
      return std::numeric_limits<double>::infinity();
    }
    const Eigen::VectorXd marginal = likelihood * value;
    if ((marginal.array() <= 0.0).any() || !marginal.allFinite()) {
      return std::numeric_limits<double>::infinity();
    }
    return -marginal.array().log().sum() - mu * value.array().log().sum();
  };
  while (barrier > std::max(tolerance * 0.1, 1e-12) &&
         iterations < maxit) {
    const int inner_limit = std::min(50, maxit - iterations);
    for (int inner = 0; inner < inner_limit; ++inner) {
      const Eigen::VectorXd marginal =
        (likelihood * weights).array().max(minimum);
      Eigen::MatrixXd scaled = likelihood;
      for (int subject = 0; subject < subjects; ++subject) {
        scaled.row(subject) /= marginal[subject];
      }
      const Eigen::VectorXd gradient =
        -scaled.colwise().sum().transpose() -
        barrier * weights.cwiseInverse();
      Eigen::MatrixXd hessian = scaled.transpose() * scaled;
      hessian.diagonal().array() +=
        barrier * weights.array().square().inverse();
      Eigen::MatrixXd kkt(supports + 1, supports + 1);
      kkt.setZero();
      kkt.topLeftCorner(supports, supports) = hessian;
      kkt.topRightCorner(supports, 1).setOnes();
      kkt.bottomLeftCorner(1, supports).setOnes();
      Eigen::VectorXd residual(supports + 1);
      residual.head(supports) = gradient.array() + lambda;
      residual[supports] = weights.sum() - 1.0;
      const Eigen::FullPivLU<Eigen::MatrixXd> factor(kkt);
      if (!factor.isInvertible()) break;
      const Eigen::VectorXd direction = factor.solve(-residual);
      if (!direction.allFinite()) break;
      const Eigen::VectorXd delta = direction.head(supports);
      const double decrement = residual.cwiseAbs().maxCoeff();
      if (decrement <= std::max(tolerance, barrier * 0.1)) break;
      double step = 1.0;
      for (int support = 0; support < supports; ++support) {
        if (delta[support] < 0.0) {
          step = std::min(step, 0.99 * -weights[support] / delta[support]);
        }
      }
      const double current = objective(weights, barrier);
      const double directional = gradient.dot(delta);
      bool accepted = false;
      for (int line = 0; line < 40; ++line) {
        Eigen::VectorXd candidate = weights + step * delta;
        candidate /= candidate.sum();
        if ((candidate.array() > 0.0).all() &&
            objective(candidate, barrier) <=
              current + 1e-4 * step * directional) {
          weights.swap(candidate);
          lambda += step * direction[supports];
          accepted = true;
          break;
        }
        step *= 0.5;
      }
      ++iterations;
      const Eigen::VectorXd fitted =
        (likelihood * weights).array().max(minimum);
      history.push_back((row_shift.array() + fitted.array().log()).sum());
      if (!accepted) break;
    }
    barrier *= 0.1;
    Rcpp::checkUserInterrupt();
  }
  Rcpp::NumericVector weights_r(supports);
  for (int support = 0; support < supports; ++support) {
    weights_r[support] = weights[support];
  }
  Rcpp::List state = liberation_np_responsibilities(loglik, weights_r);
  return Rcpp::List::create(
    Rcpp::Named("weights") = weights_r,
    Rcpp::Named("responsibilities") = state["responsibilities"],
    Rcpp::Named("log_likelihood") = state["log_likelihood"],
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("history") = Rcpp::wrap(history),
    Rcpp::Named("solver") = "primal-dual-barrier-newton",
    Rcpp::Named("backend") = "native-cpp");
}

// [[Rcpp::export(name = ".liberation_np_responsibilities")]]
Rcpp::List liberation_np_responsibilities(
    const Rcpp::NumericMatrix& loglik,
    const Rcpp::NumericVector& weights) {
  const int subjects = loglik.nrow();
  const int supports = loglik.ncol();
  if (subjects < 1 || supports < 1 || weights.size() != supports) {
    Rcpp::stop("Nonparametric responsibility dimensions are invalid.");
  }
  Rcpp::NumericMatrix responsibilities(subjects, supports);
  double log_likelihood = 0.0;
  for (int subject = 0; subject < subjects; ++subject) {
    double maximum = -std::numeric_limits<double>::infinity();
    for (int support = 0; support < supports; ++support) {
      const double score = std::log(std::max(
        weights[support], std::numeric_limits<double>::min())) +
        loglik(subject, support);
      responsibilities(subject, support) = score;
      maximum = std::max(maximum, score);
    }
    double total = 0.0;
    for (int support = 0; support < supports; ++support) {
      const double value = std::exp(responsibilities(subject, support) - maximum);
      responsibilities(subject, support) = value;
      total += value;
    }
    if (!(total > 0.0) || !std::isfinite(total)) {
      Rcpp::stop("Nonparametric responsibilities are not finite.");
    }
    log_likelihood += maximum + std::log(total);
    for (int support = 0; support < supports; ++support) {
      responsibilities(subject, support) /= total;
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("responsibilities") = responsibilities,
    Rcpp::Named("log_likelihood") = log_likelihood);
}

// Sum subject/support population gradients under posterior support weights.
// The gradient array is stored in R column-major order with dimensions
// subject x support x parameter.
// [[Rcpp::export(name = ".liberation_np_gradient_reduce")]]
Rcpp::NumericVector liberation_np_gradient_reduce(
    const Rcpp::NumericVector& gradient,
    const Rcpp::NumericMatrix& responsibilities) {
  Rcpp::IntegerVector dimensions = gradient.attr("dim");
  if (dimensions.size() != 3 ||
      dimensions[0] != responsibilities.nrow() ||
      dimensions[1] != responsibilities.ncol()) {
    Rcpp::stop("Nonparametric gradient array dimensions are invalid.");
  }
  const int subjects = dimensions[0];
  const int supports = dimensions[1];
  const int parameters = dimensions[2];
  Rcpp::NumericVector result(parameters);
  for (int parameter = 0; parameter < parameters; ++parameter) {
    double total = 0.0;
    const R_xlen_t parameter_offset =
      static_cast<R_xlen_t>(parameter) * subjects * supports;
    for (int support = 0; support < supports; ++support) {
      const R_xlen_t support_offset = parameter_offset +
        static_cast<R_xlen_t>(support) * subjects;
      for (int subject = 0; subject < subjects; ++subject) {
        total += gradient[support_offset + subject] *
          responsibilities(subject, support);
      }
    }
    result[parameter] = total;
  }
  return result;
}

// [[Rcpp::export(name = ".liberation_waic_components")]]
Rcpp::List liberation_waic_components(const Rcpp::NumericMatrix& loglik) {
  if (loglik.nrow() < 2 || loglik.ncol() < 1) {
    Rcpp::stop("WAIC requires at least two draws and one pointwise unit.");
  }
  Rcpp::NumericVector lppd(loglik.ncol()), variance(loglik.ncol()),
    elpd(loglik.ncol()), waic(loglik.ncol());
  for (int unit = 0; unit < loglik.ncol(); ++unit) {
    double maximum = -std::numeric_limits<double>::infinity();
    double mean = 0.0;
    for (int draw = 0; draw < loglik.nrow(); ++draw) {
      const double value = loglik(draw, unit);
      if (!std::isfinite(value)) Rcpp::stop("WAIC log likelihood must be finite.");
      maximum = std::max(maximum, value);
      mean += value;
    }
    mean /= loglik.nrow();
    double exponential = 0.0;
    double square = 0.0;
    for (int draw = 0; draw < loglik.nrow(); ++draw) {
      exponential += std::exp(loglik(draw, unit) - maximum);
      const double centered = loglik(draw, unit) - mean;
      square += centered * centered;
    }
    lppd[unit] = maximum + std::log(exponential / loglik.nrow());
    variance[unit] = square / (loglik.nrow() - 1.0);
    elpd[unit] = lppd[unit] - variance[unit];
    waic[unit] = -2.0 * elpd[unit];
  }
  return Rcpp::List::create(
    Rcpp::Named("lppd") = lppd, Rcpp::Named("p_waic") = variance,
    Rcpp::Named("elpd") = elpd, Rcpp::Named("waic") = waic);
}

// [[Rcpp::export(name = ".liberation_ar1_standardize")]]
Rcpp::NumericVector liberation_ar1_standardize(
    const Rcpp::LogicalVector& observed,
    const Rcpp::IntegerVector& groups, double rho, bool use_ar1) {
  if (groups.size() != observed.size()) {
    Rcpp::stop("Residual group and observation vectors must have equal length.");
  }
  int maximum_group = 0;
  for (int group : groups) {
    if (group != NA_INTEGER) maximum_group = std::max(maximum_group, group);
  }
  Rcpp::NumericVector standardized(observed.size());
  const double innovation_scale = std::sqrt(std::max(0.0, 1.0 - rho * rho));
  for (int group = 1; group <= maximum_group; ++group) {
    bool first = true;
    double previous = 0.0;
    for (R_xlen_t row = 0; row < observed.size(); ++row) {
      if (observed[row] != TRUE || groups[row] != group) continue;
      const double innovation = R::rnorm(0.0, 1.0);
      const double current = use_ar1 && !first ?
        rho * previous + innovation_scale * innovation : innovation;
      standardized[row] = current;
      previous = current;
      first = false;
    }
  }
  return standardized;
}

// [[Rcpp::export(name = ".liberation_mixture_component_nll")]]
Rcpp::NumericMatrix liberation_mixture_component_nll(
    SEXP engine_pointer, const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta, const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma) {
  liberation::require_materialized_addl(data);
  Rcpp::XPtr<liberation::ModelEngine> engine(engine_pointer);
  if (engine->mixture_probabilities.empty()) {
    Rcpp::stop("The model does not define a finite mixture.");
  }
  Rcpp::IntegerVector subject_index = data[".ID_INDEX"];
  int n_subjects = 0;
  for (int value : subject_index) n_subjects = std::max(n_subjects, value);
  std::vector<double> theta_values = Rcpp::as<std::vector<double>>(theta);
  std::vector<double> eta_values;
  eta_values.reserve(static_cast<std::size_t>(eta.size()));
  for (int row = 0; row < eta.nrow(); ++row) {
    for (int column = 0; column < eta.ncol(); ++column) eta_values.push_back(eta(row, column));
  }
  std::vector<double> sigma_values = Rcpp::as<std::vector<double>>(sigma);
  Rcpp::NumericMatrix result(n_subjects, engine->mixture_probabilities.size());
  for (std::size_t component = 0; component < engine->mixture_probabilities.size(); ++component) {
    std::vector<int> assignment(static_cast<std::size_t>(n_subjects),
                                static_cast<int>(component + 1));
    std::vector<double> prediction = liberation::simulate_analytical_t(
      *engine, data, theta_values, eta_values, sigma_values, assignment);
    std::vector<double> nll = liberation::residual_subject_nll_t(
      *engine, data, prediction, theta_values, eta_values, sigma_values,
      assignment);
    for (int subject = 0; subject < n_subjects; ++subject) {
      result(subject, static_cast<int>(component)) = nll[static_cast<std::size_t>(subject)];
    }
  }
  return result;
}
