#ifndef LIBERATION_NATIVE_OPTIMIZER_API_HPP
#define LIBERATION_NATIVE_OPTIMIZER_API_HPP

#include <Rcpp.h>
#include <LibeRtAD/eigen_r.hpp>

#include <functional>

namespace liberation {

using NativeValueFunction = std::function<double(const Eigen::VectorXd&)>;
using NativeGradientFunction =
  std::function<Eigen::VectorXd(const Eigen::VectorXd&)>;

Rcpp::List native_optimizer_core(
    const NativeValueFunction& objective,
    const NativeGradientFunction& gradient,
    const Rcpp::NumericVector& start,
    const Rcpp::NumericVector& lower,
    const Rcpp::NumericVector& upper,
    int maxit, double tolerance, int trace);

}  // namespace liberation

#endif
