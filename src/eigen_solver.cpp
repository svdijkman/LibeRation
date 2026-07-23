#if defined(__GNUC__) && !defined(__clang__)
// Eigen tracks this optimizer false positive as upstream issue #2304 and
// suppresses it for GCC in its own test builds.
#define LIBERATION_GCC_DIAGNOSTIC(x) _Pragma(#x)
LIBERATION_GCC_DIAGNOSTIC(GCC diagnostic ignored "-Wmaybe-uninitialized")
#endif

#include "eigen_solver.h"

namespace liberation {
namespace detail {

SelfAdjointEigenResult self_adjoint_eigen(
    const Eigen::MatrixXd& matrix, bool compute_vectors) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
  solver.compute(matrix, compute_vectors ? Eigen::ComputeEigenvectors :
                                           Eigen::EigenvaluesOnly);
  SelfAdjointEigenResult output{};
  output.info = solver.info();
  if (output.info == Eigen::Success) {
    output.values = solver.eigenvalues();
    if (compute_vectors) output.vectors = solver.eigenvectors();
  }
  return output;
}

}  // namespace detail
}  // namespace liberation

#if defined(__GNUC__) && !defined(__clang__)
#undef LIBERATION_GCC_DIAGNOSTIC
#endif
