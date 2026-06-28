#include <R.h>
#include <Rcpp.h>
#include "nm_pk_ad.h"

namespace {

typedef void (*register_custom_reverse_fn)(const char*, void (*)(SEXP, SEXP));

void wrap_reverse_pk_bolus1(SEXP node, SEXP grad) {
  nm_pk_ad::reverse_pk_bolus1(Rcpp::Environment(node), grad);
}

void wrap_reverse_pk_oral1(SEXP node, SEXP grad) {
  nm_pk_ad::reverse_pk_oral1(Rcpp::Environment(node), grad);
}

void wrap_reverse_pk_oral2_trans4(SEXP node, SEXP grad) {
  nm_pk_ad::reverse_pk_oral2_trans4(Rcpp::Environment(node), grad);
}

void register_pk_reverse_handlers() {
  register_custom_reverse_fn reg = (register_custom_reverse_fn)
      R_GetCCallable("LibeRtAD", "register_custom_reverse");
  if (reg == NULL) {
    Rf_warning("LibeRtAD custom reverse registry not found; PK AD gradients may fail");
    return;
  }
  reg("pk_bolus1", wrap_reverse_pk_bolus1);
  reg("pk_oral1", wrap_reverse_pk_oral1);
  reg("pk_oral2_trans4", wrap_reverse_pk_oral2_trans4);
}

}  // namespace

// [[Rcpp::init]]
void nm_register_pk_ad_reverse(DllInfo* dll) {
  register_pk_reverse_handlers();
}
