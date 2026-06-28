#' FOCEI conditional residuals and objective on any fit (nlmixr2-style addCwres)
#'
#' Adds population predictions, conditional predictions (CPRED), conditional
#' residuals (CRES, CWRES), and a FOCEI objective evaluated at the fit's
#' post-hoc ETAs. Works for FO, FOCE, SAEM, etc., not only FOCEI fits.
#'
#' @param fit An \code{nm_fit} object.
#' @param data Optional dataset; defaults to \code{fit$data}.
#' @param refit_eta If \code{TRUE}, re-fit subject ETAs at \code{fit$par} before GOF.
#' @param pk_engine PK engine for predictions and optional eta refit.
#' @param control Passed to eta refit when \code{refit_eta = TRUE}.
#' @return The fit object with \code{gof} (row-level table), \code{cwres_obj},
#'   and \code{cwres_eta} attached.
#' @examples
#' \dontrun{
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "FO",
#'               control = list(maxit = 5L, compute_inference = FALSE))
#' fit <- nm_add_cwres(fit)
#' head(fit$gof[, c("ID", "TIME", "CWRES")])
#' }
#' @export
nm_add_cwres <- function(fit,
                         data = fit$data,
                         refit_eta = FALSE,
                         pk_engine = NULL,
                         control = list()) {
  if (is.null(fit) || is.null(fit$model)) {
    .nm_stop("fit must be an nm_fit object with a model.")
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    .nm_stop("Package 'data.table' is required.")
  }
  model <- fit$model
  pk <- pk_engine %||% fit$pk_engine %||% "auto"
  .nm_sync_lik_config(model)
  dat <- .nm_prepare_data(data, model$INPUT, model)
  ids <- .nm_subject_ids(dat)
  n_eta <- .nm_n_eta(model)
  p <- .nm_unpack(model, fit$par)
  eta_mat <- .nm_fit_eta_matrix(fit, data = data, control = control)
  if (isTRUE(refit_eta) && n_eta > 0L) {
    eta_init <- .nm_mceta_init(model, dat, p$theta, p$omega, p$sigma, eta_mat, control)
    eta_mat <- .nm_fit_all_eta_modes(
      model, dat, p$theta, p$omega, p$sigma, eta_init, pk, control
    )
    eta_mat <- .nm_sanitize_eta_mat(eta_mat, p$omega)
  }
  cfg <- .nm_lik_config(model)
  err <- cfg$error
  out <- data.table::copy(dat)
  out$IPRED <- 0
  out$PRED <- 0
  out$CPRED <- NA_real_
  out$CRES <- NA_real_
  out$CWRES <- NA_real_
  out$RES <- NA_real_
  out$WRES <- NA_real_
  out$IWRES <- NA_real_
  for (j in seq_along(ids)) {
    id <- ids[j]
    subj <- .nm_subject_slice(dat, id)
    eta_j <- if (n_eta == 0L) {
      numeric(0)
    } else if (is.matrix(eta_mat) && nrow(eta_mat) >= j) {
      eta_mat[j, ]
    } else {
      numeric(n_eta)
    }
    pred <- .nm_subject_ipred(model, subj, p$theta, p$omega, eta_j, pk_engine = pk)
    pred_pop <- if (n_eta > 0L) {
      .nm_subject_ipred(model, subj, p$theta, p$omega, numeric(n_eta), pk_engine = pk)
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
      var_pop <- .nm_residual_var(f, p$sigma, err)
      out$WRES[obs_rows] <- out$RES[obs_rows] / sqrt(var_pop)
      out$IWRES[obs_rows] <- out$WRES[obs_rows]
      gh <- .nm_focei_subject_G_sens(
        model, subj, p$theta, p$omega, eta_j, p$sigma, pk
      )
      G <- gh$G
      OM <- .nm_focei_omega_matrix(p$omega, n_eta, cfg$omega)
      for (k in seq_along(obs_rows)) {
        fj <- as.numeric(.nm_focei_pick_f(f, k))
        rj <- .nm_focei_residual_var(fj, p$sigma, err)
        g_row <- if (n_eta > 0L && ncol(G) > 0L) {
          vapply(seq_len(n_eta), function(ii) {
            as.numeric(.nm_focei_G_ij(G, k, ii, ad = FALSE))
          }, numeric(1))
        } else {
          numeric(0)
        }
        gvg <- if (length(g_row) > 0L) {
          sum((OM %*% g_row) * g_row)
        } else {
          0
        }
        v_cond <- max(as.numeric(rj) + gvg, .Machine$double.eps)
        out$CPRED[obs_rows[k]] <- fj
        out$CRES[obs_rows[k]] <- dv[k] - fj
        out$CWRES[obs_rows[k]] <- out$CRES[obs_rows[k]] / sqrt(v_cond)
      }
    }
  }
  cwres_obj <- if (n_eta > 0L && is.matrix(eta_mat)) {
    .nm_focei_objective(model, data, p$theta, p$omega, p$sigma, eta_mat, pk)
  } else {
    .nm_nll_internal(
      model, data, p$theta, p$omega, p$sigma,
      eta = eta_mat, include_omega_prior = TRUE, pk_engine = pk
    )
  }
  fit$gof <- out
  fit$cwres_obj <- cwres_obj
  fit$cwres_eta <- eta_mat
  fit
}
