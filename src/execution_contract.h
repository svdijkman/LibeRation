#ifndef LIBERATION_EXECUTION_CONTRACT_H
#define LIBERATION_EXECUTION_CONTRACT_H

#include <Rcpp.h>
#include <cmath>

namespace liberation {

inline void require_materialized_addl(const Rcpp::DataFrame& data) {
  if (!data.containsElementNamed("ADDL")) return;
  const Rcpp::NumericVector addl = Rcpp::as<Rcpp::NumericVector>(data["ADDL"]);
  for (R_xlen_t row = 0; row < addl.size(); ++row) {
    if (!std::isfinite(addl[row]) || addl[row] != 0.0) {
      Rcpp::stop(
        "Native execution requires ADDL/II doses to be materialized by "
        "nm_dataset(); a non-zero or invalid ADDL value reached C++."
      );
    }
  }
}

inline void require_materialized_addl(const Rcpp::List& subject_data) {
  for (R_xlen_t index = 0; index < subject_data.size(); ++index) {
    SEXP input = subject_data[index];
    // Dynamic vectors and native subject-view descriptors do not carry an
    // event table. Their canonical parent table was validated when its native
    // store was created.
    if (!Rf_inherits(input, "data.frame")) continue;
    require_materialized_addl(Rcpp::as<Rcpp::DataFrame>(input));
  }
}

}  // namespace liberation

#endif
