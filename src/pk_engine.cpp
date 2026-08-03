// [[Rcpp::depends(LibeRtAD)]]
// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>
#include <LibeRtAD/eigen_r.hpp>
#include "execution_contract.h"
#include "population_objective_api.h"
#include <LibeRtAD/sparse_hessian.hpp>
#include <unsupported/Eigen/MatrixFunctions>
#include <LibeRtAD/program.hpp>
#include <LibeRtAD/eigen_solver.hpp>

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <numeric>
#include <queue>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace liberation {

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

// Internal implementation units are included into one coordinator TU so
// CppAD/Eigen template definitions remain visible without duplication.
// The population-objective R boundary is separately compiled in
// population_objective_api.cpp.
#include "pk_engine_event_advan.h"
#include "pk_engine_differential_systems.h"
#include "pk_engine_ad_propagation.h"
#include "pk_engine_likelihood.h"
#include "pk_engine_population.h"
#include "pk_engine_state_space.h"

}  // namespace liberation

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
    SEXP engine_pointer, const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta, const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma) {
  Rcpp::XPtr<liberation::ModelEngine> engine(engine_pointer);
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
    SEXP tape_pointer, const Rcpp::DataFrame& data) {
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

// [[Rcpp::export(name = ".liberation_fo_tape_new_dynamic")]]
Rcpp::NumericVector liberation_fo_tape_new_dynamic(
    SEXP tape_pointer, const Rcpp::DataFrame& data) {
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
    const Rcpp::DataFrame& data, const Rcpp::NumericVector& theta,
    const Rcpp::NumericVector& sigma, const Rcpp::NumericVector& omega) {
  liberation::require_materialized_addl(data);
  Rcpp::XPtr<liberation::ModelEngine> engine(engine_pointer);
  Rcpp::XPtr<liberation::PredictionTape> prediction_tape(prediction_tape_pointer);
  std::unique_ptr<liberation::ObjectiveTape> tape = liberation::record_fo_tape(
    *engine, *prediction_tape, data, theta, sigma, omega);
  Rcpp::XPtr<liberation::ObjectiveTape> pointer(tape.release(), true);
  pointer.attr("class") = Rcpp::CharacterVector::create(
    "liberation_fo_tape_ptr", "liberation_objective_tape_ptr", "externalptr");
  pointer.attr("domain") = Rcpp::wrap(pointer->domain_names);
  return pointer;
}

// [[Rcpp::export(name = ".liberation_curvature_tape_create")]]
SEXP liberation_curvature_tape_create(
    SEXP engine_pointer, SEXP prediction_tape_pointer,
    SEXP objective_tape_pointer, const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta, const Rcpp::NumericVector& eta,
    const Rcpp::NumericVector& sigma, const Rcpp::NumericVector& omega,
    const std::string& approximation) {
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
    SEXP engine_pointer, const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta, const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma, const Rcpp::NumericVector& omega,
    bool interaction = true) {
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
    const Rcpp::List& tape_pointers, const Rcpp::NumericMatrix& points) {
  if (points.nrow() != tape_pointers.size()) {
    Rcpp::stop("Objective point rows must match the number of tapes.");
  }
  Rcpp::NumericVector values(points.nrow());
  for (int row = 0; row < points.nrow(); ++row) {
    Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointers[row]);
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
    const Rcpp::List& tape_pointers, const Rcpp::NumericMatrix& points) {
  if (points.nrow() != tape_pointers.size()) {
    Rcpp::stop("Objective point rows must match the number of tapes.");
  }
  if (!points.nrow()) return Rcpp::NumericMatrix(0, points.ncol());
  Rcpp::NumericMatrix gradients(points.nrow(), points.ncol());
  const std::vector<double> weight(1, 1.0);
  for (int row = 0; row < points.nrow(); ++row) {
    Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointers[row]);
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

// [[Rcpp::export(name = ".liberation_objective_tape_eta_modes")]]
Rcpp::List liberation_objective_tape_eta_modes(
    const Rcpp::List& tape_pointers, const Rcpp::NumericMatrix& points,
    const Rcpp::IntegerVector& eta_positions,
    const Rcpp::NumericMatrix& starts, int maxit = 100,
    double tolerance = 1e-7, bool exact_hessian = true) {
  if (points.nrow() != tape_pointers.size() || starts.nrow() != points.nrow()) {
    Rcpp::stop("ETA-mode rows must match the number of objective tapes.");
  }
  if (starts.ncol() != eta_positions.size()) {
    Rcpp::stop("ETA starting values have the wrong number of columns.");
  }
  Rcpp::List result(points.nrow());
  for (int row = 0; row < points.nrow(); ++row) {
    Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointers[row]);
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
      exact_hessian);
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
    double step_scale = 0.5) {
  const int subjects = points.nrow();
  const int dimension = eta_positions.size();
  if (subjects != tape_pointers.size() || current_eta.nrow() != subjects ||
      current_eta.ncol() != dimension || proposal_roots.size() != subjects ||
      mcmc_steps < 1 || normals.nrow() != subjects * mcmc_steps ||
      normals.ncol() != dimension || log_uniforms.size() != normals.nrow() ||
      !std::isfinite(step_scale) || step_scale <= 0.0) {
    Rcpp::stop("Batched ETA Metropolis inputs are inconsistent.");
  }
  Rcpp::NumericMatrix eta = Rcpp::clone(current_eta);
  Rcpp::NumericVector values(subjects);
  int accepted = 0;
  for (int subject = 0; subject < subjects; ++subject) {
    Rcpp::XPtr<liberation::ObjectiveTape> tape(tape_pointers[subject]);
    const std::size_t domain = tape->domain_names.size();
    if (points.ncol() != static_cast<int>(domain)) {
      Rcpp::stop("A batched Metropolis objective point has the wrong length.");
    }
    std::vector<std::size_t> positions;
    positions.reserve(static_cast<std::size_t>(dimension));
    for (int value : eta_positions) {
      if (value < 1 || static_cast<std::size_t>(value) > domain) {
        Rcpp::stop("ETA position is outside a Metropolis objective tape domain.");
      }
      positions.push_back(static_cast<std::size_t>(value - 1));
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
    liberation::EtaEvaluation current = liberation::objective_eta_evaluate(
      *tape, point, positions, false);
    if (!current.finite) Rcpp::stop("Current ETA objective is not finite.");
    for (int step = 0; step < mcmc_steps; ++step) {
      const int draw = subject * mcmc_steps + step;
      std::vector<double> candidate_point = point;
      Rcpp::NumericVector candidate_eta(dimension);
      for (int row = 0; row < dimension; ++row) {
        double increment = 0.0;
        for (int column = 0; column < dimension; ++column) {
          increment += root(row, column) * normals(draw, column);
        }
        candidate_eta[row] = eta(subject, row) + step_scale * increment;
        candidate_point[positions[static_cast<std::size_t>(row)]] = candidate_eta[row];
      }
      liberation::EtaEvaluation candidate = liberation::objective_eta_evaluate(
        *tape, candidate_point, positions, false);
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
    Rcpp::Named("attempted") = subjects * mcmc_steps);
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
