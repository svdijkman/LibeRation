#ifndef LIBERATION_NATIVE_EIGEN_SOLVER_HPP
#define LIBERATION_NATIVE_EIGEN_SOLVER_HPP

#include <LibeRtAD/eigen.hpp>

namespace liberation {
namespace detail {

struct SelfAdjointEigenResult {
  Eigen::ComputationInfo info = Eigen::InvalidInput;
  Eigen::VectorXd values;
  Eigen::MatrixXd vectors;
};

SelfAdjointEigenResult self_adjoint_eigen(
  const Eigen::MatrixXd& matrix, bool compute_vectors = true);

}  // namespace detail
}  // namespace liberation

#endif
