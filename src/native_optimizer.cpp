// [[Rcpp::depends(LibeRtAD)]]
// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>
#include <LibeRtAD/eigen_r.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace {

double optimizer_value(const Rcpp::Function& function,
                       const Eigen::VectorXd& point,
                       const Eigen::VectorXd& scale) {
  Rcpp::NumericVector native(point.size());
  for (Eigen::Index i = 0; i < point.size(); ++i) native[i] = point[i] * scale[i];
  Rcpp::NumericVector value = function(native);
  if (value.size() != 1 || !std::isfinite(value[0])) return 1e100;
  return value[0];
}

Eigen::VectorXd optimizer_gradient(const Rcpp::Function& function,
                                   const Eigen::VectorXd& point,
                                   const Eigen::VectorXd& scale) {
  Rcpp::NumericVector native(point.size());
  for (Eigen::Index i = 0; i < point.size(); ++i) native[i] = point[i] * scale[i];
  Rcpp::NumericVector source = function(native);
  if (source.size() != point.size()) {
    Rcpp::stop("The native optimizer gradient has the wrong length.");
  }
  Eigen::VectorXd result(point.size());
  for (Eigen::Index i = 0; i < point.size(); ++i) {
    result[i] = source[i] * scale[i];
    if (!std::isfinite(result[i])) {
      Rcpp::stop("The native optimizer gradient is not finite.");
    }
  }
  return result;
}

Eigen::VectorXd projected_gradient(const Eigen::VectorXd& point,
                                   const Eigen::VectorXd& gradient,
                                   const Eigen::VectorXd& lower,
                                   const Eigen::VectorXd& upper) {
  Eigen::VectorXd result = gradient;
  for (Eigen::Index i = 0; i < point.size(); ++i) {
    const double margin = 1e-12 * std::max(1.0, std::abs(point[i]));
    if ((point[i] <= lower[i] + margin && gradient[i] > 0.0) ||
        (point[i] >= upper[i] - margin && gradient[i] < 0.0)) {
      result[i] = 0.0;
    }
  }
  return result;
}

}  // namespace

// A scaled, box-constrained BFGS implementation keeps all line-search and
// convergence bookkeeping in compiled code while accepting the population
// objective/gradient closures used by the R estimation layer.
// [[Rcpp::export(name = ".liberation_native_optimizer")]]
Rcpp::List liberation_native_optimizer(
    const Rcpp::Function& objective, const Rcpp::Function& gradient,
    const Rcpp::NumericVector& start, const Rcpp::NumericVector& lower,
    const Rcpp::NumericVector& upper, int maxit = 200,
    double tolerance = 1e-6, int trace = 0) {
  const Eigen::Index dimension = start.size();
  if (lower.size() != dimension || upper.size() != dimension || maxit < 1 ||
      !std::isfinite(tolerance) || tolerance <= 0.0) {
    Rcpp::stop("Native optimizer controls or bounds are invalid.");
  }
  Eigen::VectorXd scale(dimension), point(dimension), low(dimension), high(dimension);
  for (Eigen::Index i = 0; i < dimension; ++i) {
    scale[i] = std::max(std::abs(start[i]), 1.0);
    point[i] = start[i] / scale[i];
    low[i] = lower[i] / scale[i];
    high[i] = upper[i] / scale[i];
    if (low[i] > high[i] || point[i] < low[i] || point[i] > high[i]) {
      Rcpp::stop("Native optimizer start is outside its bounds.");
    }
  }
  double value = optimizer_value(objective, point, scale);
  Eigen::VectorXd derivative = optimizer_gradient(gradient, point, scale);
  int function_evaluations = 1;
  int gradient_evaluations = 1;
  int convergence = 1;
  int iterations = 0;
  std::string message = "iteration limit reached";
  const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(dimension, dimension);
  Eigen::MatrixXd initial_inverse = identity;
  for (Eigen::Index i = 0; i < dimension; ++i) {
    initial_inverse(i, i) = 1.0 / (scale[i] * scale[i]);
  }
  Eigen::MatrixXd inverse = initial_inverse;
  std::vector<int> trace_iteration;
  std::vector<double> trace_value, trace_gradient, trace_step;

  for (int iteration = 0; iteration < maxit; ++iteration) {
    Eigen::VectorXd projected = projected_gradient(point, derivative, low, high);
    const double norm = dimension ? projected.lpNorm<Eigen::Infinity>() : 0.0;
    trace_iteration.push_back(iteration);
    trace_value.push_back(value);
    trace_gradient.push_back(norm);
    trace_step.push_back(iteration ? trace_step.back() : 0.0);
    if (trace > 0) {
      Rcpp::Rcout << "[LibeRation/native] ITERATION " << iteration
                  << " OFV " << value << " PROJECTED_GRADIENT " << norm << "\n";
    }
    if (norm <= std::max(tolerance, 1e-8) * (1.0 + std::abs(value))) {
      convergence = 0;
      message = "projected gradient tolerance reached";
      iterations = iteration;
      break;
    }
    Eigen::VectorXd direction = -inverse * projected;
    for (Eigen::Index i = 0; i < dimension; ++i) {
      if (projected[i] == 0.0) direction[i] = 0.0;
    }
    double directional = derivative.dot(direction);
    if (!std::isfinite(directional) || directional >= -1e-14) {
      inverse = initial_inverse;
      direction = -inverse * projected;
      directional = derivative.dot(direction);
    }
    double maximum_step = 1.0;
    for (Eigen::Index i = 0; i < dimension; ++i) {
      if (direction[i] > 0.0 && std::isfinite(high[i])) {
        maximum_step = std::min(maximum_step, (high[i] - point[i]) / direction[i]);
      } else if (direction[i] < 0.0 && std::isfinite(low[i])) {
        maximum_step = std::min(maximum_step, (low[i] - point[i]) / direction[i]);
      }
    }
    double step = std::max(0.0, maximum_step);
    Eigen::VectorXd candidate = point;
    double candidate_value = 1e100;
    bool accepted = false;
    for (int line_search = 0; line_search < 40 && step > 1e-16; ++line_search) {
      candidate = point + step * direction;
      candidate = candidate.cwiseMax(low).cwiseMin(high);
      candidate_value = optimizer_value(objective, candidate, scale);
      ++function_evaluations;
      if (std::isfinite(candidate_value) &&
          candidate_value <= value + 1e-4 * step * directional) {
        accepted = true;
        break;
      }
      step *= 0.5;
    }
    if (!accepted) {
      convergence = 52;
      message = "line search failed";
      iterations = iteration;
      break;
    }
    Eigen::VectorXd candidate_derivative = optimizer_gradient(gradient, candidate, scale);
    ++gradient_evaluations;
    const Eigen::VectorXd displacement = candidate - point;
    const Eigen::VectorXd change = candidate_derivative - derivative;
    const double curvature = displacement.dot(change);
    if (std::isfinite(curvature) &&
        curvature > 1e-12 * displacement.norm() * change.norm()) {
      const double rho = 1.0 / curvature;
      const Eigen::MatrixXd left = identity - rho * displacement * change.transpose();
      inverse = left * inverse * left.transpose() +
        rho * displacement * displacement.transpose();
    } else {
      inverse = initial_inverse;
    }
    const double previous = value;
    point = candidate;
    value = candidate_value;
    derivative = candidate_derivative;
    iterations = iteration + 1;
    trace_step.back() = step;
    if (std::abs(previous - value) <= tolerance * (1.0 + std::abs(value))) {
      Eigen::VectorXd next_projected = projected_gradient(point, derivative, low, high);
      if (!dimension || next_projected.lpNorm<Eigen::Infinity>() <=
          std::sqrt(tolerance) * (1.0 + std::abs(value))) {
        convergence = 0;
        message = "relative objective and gradient tolerance reached";
        break;
      }
    }
    if ((iteration + 1) % 10 == 0) Rcpp::checkUserInterrupt();
  }
  Rcpp::NumericVector par(dimension), final_gradient(dimension);
  for (Eigen::Index i = 0; i < dimension; ++i) {
    par[i] = point[i] * scale[i];
    final_gradient[i] = derivative[i] / scale[i];
  }
  Rcpp::IntegerVector counts = Rcpp::IntegerVector::create(
    Rcpp::Named("function") = function_evaluations,
    Rcpp::Named("gradient") = gradient_evaluations);
  return Rcpp::List::create(
    Rcpp::Named("par") = par,
    Rcpp::Named("value") = value,
    Rcpp::Named("convergence") = convergence,
    Rcpp::Named("message") = message,
    Rcpp::Named("counts") = counts,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("objective_evaluations") = function_evaluations,
    Rcpp::Named("gradient_evaluations") = gradient_evaluations,
    Rcpp::Named("gradient") = final_gradient,
    Rcpp::Named("telemetry") = Rcpp::DataFrame::create(
      Rcpp::Named("iteration") = trace_iteration,
      Rcpp::Named("objective") = trace_value,
      Rcpp::Named("projected_gradient") = trace_gradient,
      Rcpp::Named("step") = trace_step));
}

