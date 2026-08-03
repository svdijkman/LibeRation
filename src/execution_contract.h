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
    require_materialized_addl(Rcpp::as<Rcpp::DataFrame>(subject_data[index]));
  }
}

}  // namespace liberation

#endif
