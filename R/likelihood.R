#' Configure population likelihood features
#'
#' @param error Residual model. `auto` uses the model's `ERROR` block.
#' @param omega Random-effect covariance structure.
#' @param sigma_corr Independent or within-subject AR(1) residual errors.
#' @param sigma_parameterization Interpret SIGMA values as residual standard
#'   deviations (the historical LibeRation convention) or variances (the
#'   native NONMEM `$SIGMA` convention).
#' @param ar1_rho Fixed AR(1) correlation in (-1, 1).
#' @param blq_method No censoring, Beal M3, or Beal M4.
#' @param lloq Optional scalar LLOQ used when a row-specific `LLOQ` column is
#'   absent.
#' @param iov Number of trailing ETA definitions that are occasion-specific.
#' @param occasion_col Dataset column identifying occasions.
#' @param priors Optional fixed-effect prior definitions.
#' @param mixtures Optional mixture-model definitions.
#' @return A validated serializable likelihood configuration.
#' @export
nm_lik_config <- function(
    error = c("auto", "none", "combined", "additive", "proportional", "exponential", "power"),
    omega = c("diagonal", "full"),
    sigma_corr = c("independent", "ar1"),
    sigma_parameterization = c("sd", "variance"),
    ar1_rho = 0,
    blq_method = c("none", "m3", "m4"),
    lloq = NA_real_,
    iov = 0L,
    occasion_col = "OCC",
    priors = NULL,
    mixtures = NULL) {
  error <- match.arg(error)
  omega <- match.arg(omega)
  sigma_corr <- match.arg(sigma_corr)
  sigma_parameterization <- match.arg(sigma_parameterization)
  blq_method <- match.arg(blq_method)
  ar1_rho <- as.numeric(ar1_rho)
  lloq <- as.numeric(lloq)
  iov <- as.integer(iov)
  if (length(ar1_rho) != 1L || !is.finite(ar1_rho) || abs(ar1_rho) >= 1) {
    .nm_stop("`ar1_rho` must be a finite scalar strictly between -1 and 1.")
  }
  if (length(lloq) != 1L || (!is.na(lloq) && (!is.finite(lloq) || lloq <= 0))) {
    .nm_stop("`lloq` must be NA or a positive finite scalar.")
  }
  if (length(iov) != 1L || is.na(iov) || iov < 0L) {
    .nm_stop("`iov` must be a non-negative integer.")
  }
  if (length(occasion_col) != 1L || !nzchar(occasion_col)) {
    .nm_stop("`occasion_col` must be one non-empty column name.")
  }
  if (!is.null(priors)) {
    priors <- as.data.frame(priors, stringsAsFactors = FALSE)
    required <- c("parameter", "distribution", "mean", "sd", "shape", "rate")
    if (!all(required %in% names(priors))) {
      .nm_stop("`priors` must be built by row-binding `nm_prior()` results.")
    }
    priors$parameter <- toupper(gsub("_", "", priors$parameter))
  }
  if (!is.null(mixtures)) {
    if (is.numeric(mixtures)) mixtures <- nm_mixture(mixtures)
    if (!inherits(mixtures, "nm_mixture")) {
      .nm_stop("`mixtures` must be created by `nm_mixture()` or be a probability vector.")
    }
  }
  structure(
    list(
      version = 1L, error = error, omega = omega,
      sigma_corr = sigma_corr, sigma_parameterization = sigma_parameterization,
      ar1_rho = ar1_rho,
      blq_method = blq_method, lloq = lloq, iov = iov,
      occasion_col = as.character(occasion_col), priors = priors,
      mixtures = mixtures
    ),
    class = "nm_lik_config"
  )
}

.nm_lik_config <- function(config, error_type, iov) {
  if (is.null(config)) config <- nm_lik_config(error = error_type, iov = iov)
  if (!inherits(config, "nm_lik_config")) {
    if (!is.list(config)) .nm_stop("LIK_CONFIG must be created by `nm_lik_config()` or be a named list.")
    config$version <- NULL
    if (is.list(config$mixtures) && !inherits(config$mixtures, "nm_mixture")) {
      probability <- config$mixtures$probability
      label <- config$mixtures$label %||% paste0("MIX", seq_along(probability))
      config$mixtures <- nm_mixture(probability, label)
    }
    config <- do.call(nm_lik_config, config)
  }
  if (config$error == "auto") config$error <- error_type
  if (iov > 0L && config$iov == 0L) config$iov <- as.integer(iov)
  config
}

#' @export
print.nm_lik_config <- function(x, ...) {
  cat("LibeRation likelihood\n")
  cat("  residual:", x$error, " omega:", x$omega,
      " correlation:", x$sigma_corr, " SIGMA:", x$sigma_parameterization,
      " BLQ:", x$blq_method,
      " IOV:", x$iov, "\n")
  invisible(x)
}
