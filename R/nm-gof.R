#' @keywords internal
.nm_eta_matrix_usable <- function(eta_mat, n_sub, n_eta) {
  if (n_eta == 0L) {
    return(is.null(eta_mat))
  }
  is.matrix(eta_mat) &&
    nrow(eta_mat) >= n_sub &&
    ncol(eta_mat) >= n_eta &&
    any(is.finite(eta_mat)) &&
    max(abs(eta_mat), na.rm = TRUE) > 1e-12
}

#' @keywords internal
.nm_fit_eta_matrix <- function(fit, data = fit$data, control = list()) {
  model <- fit$model
  n_eta <- .nm_n_eta(model)
  if (n_eta == 0L) {
    return(NULL)
  }
  dat <- .nm_prepare_data(data, model$INPUT, model)
  n_sub <- length(.nm_subject_ids(dat))
  eta_mat <- fit$eta
  if (.nm_eta_matrix_usable(eta_mat, n_sub, n_eta)) {
    return(eta_mat)
  }
  p <- .nm_unpack(model, fit$par)
  pk <- fit$pk_engine %||% "auto"
  .nm_fit_all_eta_modes(
    model, dat, p$theta, p$omega, p$sigma,
    if (is.matrix(eta_mat)) eta_mat else NULL,
    pk, control
  )
}

#' Predicted values from a fit
#'
#' @param object An \code{nm_fit} object.
#' @param type \code{"ipred"} (individual) or \code{"ppred"} (population, eta = 0).
#' @param ... Unused.
#' @return A \code{data.table} with predictions and residuals columns.
#' @examples
#' \dontrun{
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "FO",
#'               control = list(maxit = 5L, compute_inference = FALSE))
#' head(predict(fit))
#' }
#' @export
predict.nm_fit <- function(object, type = c("ipred", "ppred"), ...) {
  type <- match.arg(type)
  if (!requireNamespace("data.table", quietly = TRUE)) {
    .nm_stop("Package 'data.table' is required.")
  }
  model <- object$model
  dat <- .nm_prepare_data(object$data, model$INPUT, model)
  .nm_sync_lik_config(model)
  ids <- .nm_subject_ids(dat)
  n_eta <- .nm_n_eta(model)
  err <- .nm_lik_config(model)$error
  out <- data.table::copy(dat)
  out$IPRED <- 0
  out$PRED <- 0
  out$RES <- NA_real_
  out$WRES <- NA_real_
  out$IWRES <- NA_real_
  p <- .nm_unpack(model, object$par)
  pk_eng <- object$pk_engine %||% "auto"
  eta_mat <- if (type == "ipred" && n_eta > 0L) {
    .nm_fit_eta_matrix(object, data = object$data)
  } else {
    object$eta
  }
  for (j in seq_along(ids)) {
    id <- ids[j]
    sub <- .nm_subject_slice(dat, id)
    eta_j <- if (type == "ppred" || n_eta == 0L) {
      numeric(n_eta)
    } else if (is.matrix(eta_mat)) {
      eta_mat[j, ]
    } else {
      numeric(n_eta)
    }
    pred <- .nm_subject_ipred(model, sub, p$theta, p$omega, eta_j, pk_engine = pk_eng)
    pred_pop <- if (n_eta > 0L) {
      .nm_subject_ipred(
        model, sub, p$theta, p$omega, numeric(n_eta), pk_engine = pk_eng
      )
    } else {
      pred
    }
    idx_id <- out$ID == id
    out$IPRED[idx_id] <- pred$ipred
    if (length(pred$obs_idx) > 0L) {
      obs_rows <- which(idx_id & out$MDV == 0L & out$EVID == 0L)
      f <- pred$F
      f_pop <- pred_pop$F
      dv <- out$DV[obs_rows]
      out$PRED[obs_rows] <- f_pop
      out$RES[obs_rows] <- dv - f
      var <- .nm_residual_var(f, p$sigma, err)
      out$WRES[obs_rows] <- out$RES[obs_rows] / sqrt(var)
      out$IWRES[obs_rows] <- out$WRES[obs_rows]
    }
  }
  out
}

#' Extract residuals from a fit
#' @param object An \code{nm_fit} object.
#' @examples
#' \dontrun{
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "FO",
#'               control = list(maxit = 5L, compute_inference = FALSE))
#' head(residuals(fit))
#' }
#' @export
residuals.nm_fit <- function(object, ...) {
  pred <- predict(object, ...)
  pred[pred$MDV == 0L & pred$EVID == 0L,
       list(ID, TIME, DV, PRED, RES, WRES, IWRES)]
}

#' NONMEM-style ETAB (subject random effects and post-hoc summaries)
#'
#' @param fit An \code{nm_fit} object.
#' @return A \code{data.frame} with one row per subject.
#' @examples
#' \dontrun{
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "FO",
#'               control = list(maxit = 5L, compute_inference = FALSE))
#' nm_etab(fit)
#' }
#' @export
nm_etab <- function(fit) {
  model <- fit$model
  dat <- .nm_prepare_data(fit$data, model$INPUT, model)
  ids <- .nm_subject_ids(dat)
  n_eta <- .nm_n_eta(model)
  eta <- .nm_fit_eta_matrix(fit, data = fit$data)
  if (is.null(eta) || n_eta == 0L) {
    return(data.frame(ID = ids, stringsAsFactors = FALSE))
  }
  if (!is.matrix(eta)) {
    eta <- matrix(eta, nrow = length(ids), ncol = n_eta, byrow = TRUE)
  }
  cn <- paste0("ETA", seq_len(n_eta))
  tab <- data.frame(ID = ids, eta, stringsAsFactors = FALSE)
  names(tab) <- c("ID", cn)
  tab$OBJ <- NA_real_
  p <- .nm_unpack(model, fit$par)
  for (j in seq_along(ids)) {
    sub <- .nm_subject_slice(dat, ids[j])
    tab$OBJ[j] <- .nm_subject_nll_internal(
      model, sub, p$theta, p$omega, p$sigma, eta[j, ],
      include_omega_prior = TRUE, pk_engine = fit$pk_engine %||% "auto"
    )
  }
  tab
}

#' Shrinkage diagnostics for random effects
#' @param fit An \code{nm_fit} object.
#' @examples
#' \dontrun{
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "FO",
#'               control = list(maxit = 5L, compute_inference = FALSE))
#' nm_shrinkage(fit)
#' }
#' @export
nm_shrinkage <- function(fit) {
  model <- fit$model
  n_eta <- .nm_n_eta(model)
  if (n_eta == 0L) {
    return(NULL)
  }
  eta <- .nm_fit_eta_matrix(fit, data = fit$data)
  if (is.null(eta)) {
    return(NULL)
  }
  if (!is.matrix(eta)) {
    eta <- matrix(eta, ncol = n_eta)
  }
  omega <- fit$omega
  sd_eta <- apply(eta, 2, sd)
  sd_omega <- sqrt(pmax(omega[seq_len(n_eta)], 1e-8))
  shrink <- 1 - sd_eta / sd_omega
  data.frame(
    ETA = paste0("ETA", seq_len(n_eta)),
    sd_eta = sd_eta,
    sd_omega = sd_omega,
    shrinkage = shrink,
    row.names = NULL
  )
}
