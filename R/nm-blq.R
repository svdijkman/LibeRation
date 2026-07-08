#' Below-limit-of-quantification (BLQ / censoring) likelihood support
#'
#' Implements the Beal M3 and M4 methods for handling observations below the
#' lower limit of quantification (LLOQ). Censored observations contribute the
#' log-probability of being below the limit rather than a point density.
#'
#' The convention accepted for flagging a censored observation is any of:
#' \itemize{
#'   \item a \code{BLQ} column equal to 1,
#'   \item a \code{CENS} column equal to 1 (NONMEM convention), or
#'   \item \code{DV < LLOQ} when an \code{LLOQ} value is available.
#' }
#' The LLOQ threshold is taken from a per-row \code{LLOQ} column when present,
#' otherwise from \code{nm_lik_config(lloq = ...)}.
#' @name nm-blq
#' @keywords internal
NULL

#' Is the BLQ (M3/M4) likelihood active for a config?
#' @keywords internal
.nm_blq_active <- function(cfg) {
  if (is.null(cfg)) {
    return(FALSE)
  }
  m <- cfg$blq_method %||% "none"
  isTRUE(m %in% c("m3", "m4"))
}

#' Does a model use the BLQ likelihood?
#' @keywords internal
.nm_model_blq_active <- function(model) {
  .nm_blq_active(.nm_lik_config(model))
}

#' Does a model declare covariates?
#' @keywords internal
.nm_model_has_covariates <- function(model) {
  covs <- model$COVARIATES
  !is.null(covs) && length(as.character(covs)) > 0L
}

#' Does a model require the R likelihood path (rather than the monolithic C++
#' objective / eta fitter)?
#'
#' The C++ objective (nm_est.cpp) and C++ eta fitter (nm_eta_mode.cpp) do not
#' bind covariates or implement BLQ; the per-subject R likelihood loop does
#' (while still using the C++ PK solver for predictions). Models needing either
#' feature are routed through the R path.
#' @keywords internal
.nm_model_needs_r_lik <- function(model) {
  .nm_model_blq_active(model) || .nm_model_has_covariates(model)
}

#' Resolve LLOQ threshold and censoring flag for a set of observation rows.
#'
#' @param subj_ev Subject event table (post ADDL expansion).
#' @param obs_idx Integer indices of observation rows in \code{subj_ev}.
#' @param cfg Likelihood config (for the scalar \code{lloq} fallback).
#' @return List with numeric \code{lloq} and logical \code{cens} (per obs).
#' @keywords internal
.nm_blq_obs_info <- function(subj_ev, obs_idx, cfg) {
  n <- length(obs_idx)
  lloq <- rep(NA_real_, n)
  if (!is.null(cfg$lloq) && length(cfg$lloq) == 1L && is.finite(cfg$lloq)) {
    lloq[] <- as.numeric(cfg$lloq)
  }
  if ("LLOQ" %in% names(subj_ev)) {
    col <- suppressWarnings(as.numeric(subj_ev$LLOQ[obs_idx]))
    lloq[is.finite(col)] <- col[is.finite(col)]
  }
  dv <- suppressWarnings(as.numeric(subj_ev$DV[obs_idx]))
  cens <- rep(FALSE, n)
  if ("BLQ" %in% names(subj_ev)) {
    cens <- cens | (suppressWarnings(as.numeric(subj_ev$BLQ[obs_idx])) == 1)
  }
  if ("CENS" %in% names(subj_ev)) {
    cens <- cens | (suppressWarnings(as.numeric(subj_ev$CENS[obs_idx])) == 1)
  }
  if (!("BLQ" %in% names(subj_ev)) && !("CENS" %in% names(subj_ev))) {
    cens <- cens | (is.finite(lloq) & is.finite(dv) & dv < lloq)
  }
  cens[is.na(cens)] <- FALSE
  # A row can only be treated as censored when a finite LLOQ is available.
  cens <- cens & is.finite(lloq)
  list(lloq = lloq, cens = cens)
}

#' Standard-normal CDF clamped away from 0/1 for stable logs.
#' @keywords internal
.nm_blq_phi <- function(z) {
  p <- stats::pnorm(z)
  pmin(pmax(p, 1e-300), 1 - 1e-16)
}

#' Residual -2 log-likelihood with M3/M4 censoring (R path, numeric only).
#'
#' Uncensored observations use the same residual contribution as the standard
#' likelihood (\code{log(var) + resid^2 / var}); censored observations use the
#' Beal M3 or M4 log-probability of being below \code{lloq}. Consistent with the
#' package convention the additive \eqn{\log(2\pi)} constant is dropped for
#' uncensored observations; this constant is independent of the parameters so it
#' does not affect parameter estimates or likelihood-ratio comparisons.
#'
#' @keywords internal
.nm_residual_nll_blq <- function(dv, f, sigma, error = "propadd",
                                 lloq = NULL, cens = NULL, blq_method = "m3",
                                 dvid = NULL, ar1_rho = 0.0,
                                 sigma_corr = "indep") {
  n <- length(dv)
  if (n == 0L) {
    return(0)
  }
  if (is.null(cens) || !any(cens)) {
    return(.nm_residual_nll_scalar(
      dv, f, sigma, error, dvid, ar1_rho, sigma_corr
    ))
  }
  is_log <- identical(error, "log")
  nll <- 0
  for (i in seq_len(n)) {
    var <- .nm_residual_var(f[i], sigma, error, dvid, i)
    sdv <- sqrt(max(var, 1e-15))
    if (isTRUE(cens[i]) && is.finite(lloq[i])) {
      if (is_log) {
        lf <- log(max(f[i], 1e-8))
        llq <- log(max(lloq[i], 1e-8))
        z_lloq <- (llq - lf) / sdv
        p_lloq <- .nm_blq_phi(z_lloq)
        prob <- p_lloq
      } else {
        z_lloq <- (lloq[i] - f[i]) / sdv
        p_lloq <- .nm_blq_phi(z_lloq)
        if (identical(blq_method, "m4")) {
          z_zero <- (0 - f[i]) / sdv
          p_zero <- .nm_blq_phi(z_zero)
          prob <- (p_lloq - p_zero) / max(1 - p_zero, 1e-300)
        } else {
          prob <- p_lloq
        }
      }
      prob <- min(max(prob, 1e-300), 1)
      nll <- nll - 2 * log(prob)
    } else {
      resid <- if (is_log) {
        log(max(dv[i], 1e-8)) - log(max(f[i], 1e-8))
      } else {
        dv[i] - f[i]
      }
      nll <- nll + log(var) + resid^2 / var
    }
  }
  nll
}
