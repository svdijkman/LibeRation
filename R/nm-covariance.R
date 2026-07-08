#' Covariance estimators for \code{nm_fit} objects
#'
#' @param fit An \code{nm_fit} object.
#' @param data Optional dataset; defaults to \code{fit$data}.
#' @param type \code{"hessian"} (inverse Hessian), \code{"linfim"} (OPG / linearized FIM),
#'   or \code{"sandwich"} (\eqn{H^{-1} B H^{-1}} with OPG as \eqn{B}).
#' @param hessian Hessian source for \code{type = "hessian"} and sandwich.
#' @return Named covariance matrix aligned with \code{fit$par}, or \code{NULL}.
#' @examples
#' \dontrun{
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "FO",
#'               control = list(maxit = 5L, compute_inference = FALSE))
#' nm_fit_covariance(fit, type = "linfim")
#' }
#' @export
nm_fit_covariance <- function(fit,
                                data = fit$data,
                                type = c("hessian", "linfim", "sandwich"),
                                hessian = c("auto", "ad", "numeric")) {
  type <- match.arg(type)
  if (is.null(fit) || is.null(fit$par)) {
    return(NULL)
  }
  if (!is.null(fit$covariance) && !is.null(fit$covariance[[type]])) {
    return(fit$covariance[[type]])
  }
  labels <- .nm_par_labels(fit$model)
  free <- which(!.nm_fix_mask(fit$model))
  if (length(free) == 0L) {
    return(matrix(0, length(labels), length(labels), dimnames = list(labels, labels)))
  }
  if (identical(type, "linfim")) {
    opg <- .nm_fit_opg_matrix(fit, data = data)
    if (is.null(opg)) {
      return(NULL)
    }
    vcov <- .nm_stable_vcov(opg)
    if (is.null(vcov)) {
      return(NULL)
    }
    out <- matrix(NA_real_, length(labels), length(labels), dimnames = list(labels, labels))
    out[free, free] <- vcov
    return(out)
  }
  hres <- .nm_fit_hessian_at_fit(fit, data = data, hessian = hessian)
  if (is.null(hres) || is.null(hres$hessian)) {
    return(NULL)
  }
  H <- (hres$hessian + t(hres$hessian)) / 2
  free_h <- hres$free
  vcov_h <- .nm_stable_vcov(H)
  out <- matrix(NA_real_, length(labels), length(labels), dimnames = list(labels, labels))
  if (!is.null(vcov_h)) {
    out[free_h, free_h] <- vcov_h
  }
  if (identical(type, "sandwich")) {
    opg <- .nm_fit_opg_matrix(fit, data = data)
    if (is.null(opg) || is.null(vcov_h)) {
      return(out)
    }
    B <- opg[free_h, free_h, drop = FALSE]
    sw <- tryCatch(
      vcov_h %*% B %*% vcov_h,
      error = function(e) NULL
    )
    if (!is.null(sw)) {
      out[free_h, free_h] <- sw
    }
  }
  out
}

#' @keywords internal
.nm_fit_opg_matrix <- function(fit, data = fit$data) {
  model <- fit$model
  labels <- .nm_par_labels(model)
  free <- which(!.nm_fix_mask(model))
  if (length(free) == 0L) {
    return(NULL)
  }
  par <- fit$par
  pk <- fit$pk_engine %||% "auto"
  eta_mat <- fit$eta
  dat <- .nm_prepare_data(data, model$INPUT, model)
  ids <- .nm_subject_ids(dat)
  n_sub <- length(ids)
  G <- matrix(0, n_sub, length(free))
  p <- .nm_unpack(model, par)
  for (j in seq_along(ids)) {
    if (j == 1L || j == n_sub || j %% max(1L, min(5L, n_sub)) == 0L) {
      .nm_job_progress_event(
        "cov_opg",
        list(subject = j, n_sub = n_sub),
        log_msg = paste0("Covariance OPG: subject ", j, "/", n_sub)
      )
    }
    subj <- .nm_subject_slice(dat, ids[j])
    eta_j <- if (is.matrix(eta_mat) && nrow(eta_mat) >= j) {
      eta_mat[j, ]
    } else {
      numeric(.nm_n_eta(model))
    }
    g_j <- .nm_subject_score_numeric(
      fit, subj, par, eta_j, pk_engine = pk
    )
    if (!is.null(g_j)) {
      G[j, ] <- g_j[free]
    }
  }
  crossprod(G)
}

#' @keywords internal
.nm_subject_score_numeric <- function(fit, subj, par, eta, pk_engine = "auto") {
  model <- fit$model
  labels <- .nm_par_labels(model)
  fn <- function(p) {
    pp <- .nm_unpack(model, .nm_apply_fix(model, p))
    if (identical(fit$method, "FOCEI")) {
      return(.nm_focei_subject_nll_internal(
        model, subj, pp$theta, pp$omega, pp$sigma, eta, pk_engine = pk_engine
      ))
    }
    .nm_subject_nll_internal(
      model, subj, pp$theta, pp$omega, pp$sigma, eta,
      include_omega_prior = TRUE, pk_engine = pk_engine
    )
  }
  stats::setNames(.nm_num_grad(fn, par), labels)
}

#' @keywords internal
.nm_matrix_cond_number <- function(M) {
  if (is.null(M) || length(M) == 0L) {
    return(NA_real_)
  }
  ev <- tryCatch(eigen((M + t(M)) / 2, symmetric = TRUE, only.values = TRUE)$values,
                 error = function(e) NULL)
  if (is.null(ev) || length(ev) == 0L) {
    return(NA_real_)
  }
  ev <- ev[is.finite(ev) & ev > 0]
  if (length(ev) == 0L) {
    return(Inf)
  }
  max(ev) / min(ev)
}

#' @keywords internal
.nm_cov_default_method <- function(method) {
  if (method %in% c("FO", "FOCE", "FOCEI")) {
    "linfim"
  } else {
    "hessian"
  }
}

#' @keywords internal
.nm_cov_resolve_method <- function(method, fit_method) {
  method <- match.arg(method, c("auto", "hessian", "linfim", "sandwich"))
  if (identical(method, "auto")) {
    .nm_cov_default_method(fit_method)
  } else {
    method
  }
}

#' @keywords internal
.nm_cov_refit_eta <- function(fit, data = fit$data, control = list()) {
  if (is.null(fit) || identical(fit$method, "BAYES")) {
    return(fit)
  }
  model <- fit$model
  if (is.null(model)) {
    return(fit)
  }
  pk <- fit$pk_engine %||% "auto"
  dat <- .nm_prepare_data(data, model$INPUT, model)
  p <- .nm_unpack(model, fit$par)
  eta_hint <- fit$eta
  if (identical(fit$method, "FOCEI")) {
    setup <- fit$focei_setup
    if (is.null(setup)) {
      setup <- nm_focei_setup(model, data, pk_engine = pk)
      fit$focei_setup <- setup
    }
    nested <- .nm_focei_nested_objective(
      model, setup$dat, data, fit$par, eta_hint, pk, control, use_cache = FALSE
    )
    fit$eta <- nested$eta
    fit$objective <- nested$value
  } else if (fit$method %in% c("FO", "FOCE", "SAEM", "LAPLACE", "IMP")) {
    fit$eta <- .nm_fit_all_eta_modes(
      model, dat, p$theta, p$omega, p$sigma, eta_hint, pk, control
    )
    if (fit$method %in% c("FOCE", "FO")) {
      fit$objective <- .nm_fit_inference_objective(fit, fit$par, data = data)
    }
  }
  fit
}

#' @keywords internal
.nm_cov_se_from_vcov <- function(vcov, labels) {
  se <- rep(NA_real_, length(labels))
  if (is.null(vcov)) {
    return(stats::setNames(se, labels))
  }
  d <- diag(vcov)
  if (length(d) == length(labels)) {
    se <- sqrt(pmax(d, 0))
  } else {
    free <- intersect(labels, rownames(vcov))
    se[match(free, labels)] <- sqrt(pmax(diag(vcov)[match(free, rownames(vcov))], 0))
  }
  stats::setNames(se, labels)
}

#' Covariance step at final estimates (NONMEM-style)
#'
#' Re-fits subject ETAs at the final population parameters (when applicable),
#' then computes Hessian, linearized FIM (OPG), and sandwich covariance matrices.
#' Updates \code{fit$par_se} from the selected method.
#'
#' @param fit An \code{nm_fit} object.
#' @param data Optional dataset; defaults to \code{fit$data}.
#' @param method \code{"auto"} uses linearized FIM for FO/FOCE/FOCEI and inverse
#'   Hessian otherwise; or \code{"hessian"}, \code{"linfim"}, \code{"sandwich"}.
#' @param hessian Hessian source for \code{"hessian"} and \code{"sandwich"}.
#' @param refit_eta If \code{TRUE}, re-estimate ETAs at fixed THETA/OMEGA/SIGMA
#'   before computing covariance (recommended).
#' @param update_se If \code{TRUE}, set \code{fit$par_se} from \code{method}.
#' @param control Passed to ETA refit (e.g. \code{maxit} for inner eta fits).
#' @return Updated \code{nm_fit} object with \code{covariance}, \code{covariance_method},
#'   and optionally refreshed \code{par_se}.
#' @examples
#' \dontrun{
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "FO",
#'               control = list(maxit = 5L, compute_inference = FALSE))
#' nm_cov_step(fit, refit_eta = FALSE)
#' }
#' @export
nm_cov_step <- function(fit,
                        data = fit$data,
                        method = c("auto", "hessian", "linfim", "sandwich"),
                        hessian = c("auto", "ad", "numeric"),
                        refit_eta = TRUE,
                        update_se = TRUE,
                        control = list()) {
  if (is.null(fit)) {
    return(fit)
  }
  if (identical(fit$method, "BAYES")) {
    return(.nm_fit_attach_bayes_inference(fit))
  }
  hessian <- match.arg(hessian)
  cov_method <- .nm_cov_resolve_method(match.arg(method), fit$method)
  if (isTRUE(refit_eta)) {
    .nm_job_progress_event(
      "phase",
      list(step = "cov_refit_eta", method = fit$method),
      log_msg = paste0("Covariance: refitting individual ETAs (", fit$method, ")")
    )
    fit <- .nm_cov_refit_eta(fit, data = data, control = control)
  }
  .nm_seed_laplace_eta_modes(fit, data = data)
  fit <- .nm_fit_attach_covariance(fit, data = data, hessian = hessian)
  fit$covariance_method <- cov_method
  fit$cov_step <- list(
    method = cov_method,
    refit_eta = isTRUE(refit_eta),
    hessian = hessian
  )
  labels <- .nm_par_labels(fit$model)
  vc <- fit$covariance[[cov_method]]
  if (isTRUE(update_se)) {
    fit$par_se <- .nm_cov_se_from_vcov(vc, labels)
    if (!any(is.finite(unname(fit$par_se)))) {
      fit$inference_note <- paste0(
        "Covariance step (", cov_method, ") failed: matrix singular or non-finite."
      )
    } else {
      fit$inference_method <- cov_method
    }
  }
  fit
}

#' Correlation matrix from a fit covariance step
#'
#' @param fit An \code{nm_fit} object.
#' @param type Covariance type; defaults to \code{fit$covariance_method}.
#' @return Named correlation matrix, or \code{NULL}.
#' @examples
#' \dontrun{
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "FO",
#'               control = list(maxit = 5L))
#' nm_fit_correlation(fit)
#' }
#' @export
nm_fit_correlation <- function(fit, type = NULL) {
  if (is.null(fit) || is.null(fit$covariance)) {
    return(NULL)
  }
  type <- type %||% fit$covariance_method %||% "hessian"
  vc <- fit$covariance[[type]]
  if (is.null(vc)) {
    return(NULL)
  }
  sd <- sqrt(pmax(diag(vc), 0))
  if (any(!is.finite(sd)) || any(sd <= 0)) {
    return(NULL)
  }
  cor_mat <- vc / outer(sd, sd)
  diag(cor_mat) <- 1
  cor_mat
}

#' @keywords internal
.nm_fit_attach_covariance <- function(fit, data = fit$data, hessian = c("auto", "ad", "numeric")) {
  if (is.null(fit) || identical(fit$method, "BAYES")) {
    return(fit)
  }
  hessian <- match.arg(hessian)
  free <- which(!.nm_fix_mask(fit$model))
  .nm_job_progress_event(
    "phase",
    list(step = "cov_hessian", method = fit$method),
    log_msg = "Covariance: Hessian-based variance matrix"
  )
  vc_h <- nm_fit_covariance(fit, data = data, type = "hessian", hessian = hessian)
  .nm_job_progress_event(
    "phase",
    list(step = "cov_linfim", method = fit$method),
    log_msg = "Covariance: linearization (FIM) matrix"
  )
  vc_lf <- nm_fit_covariance(fit, data = data, type = "linfim")
  .nm_job_progress_event(
    "phase",
    list(step = "cov_sandwich", method = fit$method),
    log_msg = "Covariance: sandwich matrix"
  )
  vc_sw <- nm_fit_covariance(fit, data = data, type = "sandwich", hessian = hessian)
  fit$covariance <- list(
    hessian = vc_h,
    linfim = vc_lf,
    sandwich = vc_sw,
    cond_number_hessian = if (!is.null(vc_h) && length(free) > 0L) {
      .nm_matrix_cond_number(vc_h[free, free, drop = FALSE])
    } else {
      NA_real_
    },
    cond_number_linfim = if (!is.null(vc_lf) && length(free) > 0L) {
      .nm_matrix_cond_number(vc_lf[free, free, drop = FALSE])
    } else {
      NA_real_
    }
  )
  fit
}
