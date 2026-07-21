// GCC 14 reports false-positive uninitialized-pointer diagnostics inside
// Eigen's optimized self-adjoint matrix/vector kernels. Keep the
// suppression confined to this third-party adapter translation unit.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif

#include "eigen_solver.h"

namespace liberation {
namespace detail {

SelfAdjointEigenResult self_adjoint_eigen(
    const Eigen::MatrixXd& matrix, bool compute_vectors) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
  solver.compute(matrix, compute_vectors ? Eigen::ComputeEigenvectors :
                                           Eigen::EigenvaluesOnly);
  SelfAdjointEigenResult result;
  result.info = solver.info();
  if (result.info == Eigen::Success) {
    result.values = solver.eigenvalues();
    if (compute_vectors) result.vectors = solver.eigenvectors();
  }
  return result;
}

}  // namespace detail
}  // namespace liberation
