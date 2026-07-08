#pragma once
#include <Rcpp.h>

namespace nm_pk_ad {

void reverse_pk_oral2_trans4(Rcpp::Environment myvar, SEXP mygrad);
void reverse_pk_oral1(Rcpp::Environment myvar, SEXP mygrad);
void reverse_pk_bolus1(Rcpp::Environment myvar, SEXP mygrad);
void replay_pk_oral2_trans4(Rcpp::Environment node);
void replay_pk_oral1(Rcpp::Environment node);
void replay_pk_bolus1(Rcpp::Environment node);

}  // namespace nm_pk_ad
