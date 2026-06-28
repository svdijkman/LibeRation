#' C++ PK route with optional row-level S/F scaling
#'
#' @param s1,s2,s3,s4 Optional row vectors overriding $PK scaling factors (legacy).
#' @param scale_mat Optional \code{n_events x n_comp} matrix of row-level S values.
#' @param use_data_scale When \code{TRUE}, use data \code{Sx} columns when present.
#' @param f_mat Optional \code{n_events x n_comp} matrix of row-level F values.
#' @param use_data_f When \code{TRUE}, use data \code{Fx} columns when present.
#' @keywords internal
nm_pk_route_r <- function(advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
                          time, amt, rate, f1, cmt, evid, ss, ii, pk_params,
                          s1 = numeric(0), s2 = numeric(0), s3 = numeric(0), s4 = numeric(0),
                          scale_mat = matrix(numeric(0), 0L, 0L),
                          use_data_scale = FALSE,
                          f_mat = matrix(numeric(0), 0L, 0L),
                          use_data_f = FALSE) {
  .Call(
    `_LibeRation_nm_pk_route_r`,
    advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
    time, amt, rate, f1, cmt, evid, ss, ii, pk_params,
    s1, s2, s3, s4, scale_mat, use_data_scale, f_mat, use_data_f
  )
}

#' @keywords internal
nm_pk_route_detail_r <- function(advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
                                 n_state, time, amt, rate, f1, cmt, evid, ss, ii, pk_params,
                                 s1 = numeric(0), s2 = numeric(0), s3 = numeric(0), s4 = numeric(0),
                                 scale_mat = matrix(numeric(0), 0L, 0L),
                                 use_data_scale = FALSE,
                                 f_mat = matrix(numeric(0), 0L, 0L),
                                 use_data_f = FALSE) {
  .Call(
    `_LibeRation_nm_pk_route_detail_r`,
    advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss, n_state,
    time, amt, rate, f1, cmt, evid, ss, ii, pk_params,
    s1, s2, s3, s4, scale_mat, use_data_scale, f_mat, use_data_f
  )
}
