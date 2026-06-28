#pragma once
#include <Rcpp.h>
#include "nm_subject_nll.h"

namespace nm_eta {

struct SubjectData {
  Rcpp::NumericVector time, amt, rate, f1, ii, dv;
  Rcpp::IntegerVector cmt, evid, ss, obs_idx;
  nm_subject::SubjectEventExtras extras;
};

SubjectData list_to_subject(const Rcpp::List& subj);

Rcpp::NumericVector find_eta_mode(
    const SubjectData& subj,
    const Rcpp::NumericVector& theta,
    const Rcpp::NumericVector& omega,
    const Rcpp::NumericVector& sigma,
    const Rcpp::CharacterVector& pred_lines,
    int advan,
    int trans,
    int obs_cmp,
    int dose_cmp,
    int n_transit,
    bool use_ode,
    int model_ss,
    SEXP eta_init = R_NilValue,
    int max_iter = 40);

}  // namespace nm_eta
