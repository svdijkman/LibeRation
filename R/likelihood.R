#' Configure population likelihood features
#'
#' @param error Residual model. `auto` uses the model's `ERROR` block.
#'   `likelihood` compiles a user-defined row contribution from `LIK` (a
#'   probability or density) or `LOGLIK` (a log likelihood) in `$ERROR`.
#' @param omega Random-effect covariance structure.
#' @param sigma_corr Independent or within-subject AR(1) residual errors.
#' @param sigma_parameterization Interpret SIGMA values as residual standard
#'   deviations (the historical LibeRation convention) or variances (the
#'   native NONMEM `$SIGMA` convention).
#' @param ar1_rho Fixed AR(1) correlation in (-1, 1). Ignored when
#'   `ar1_parameter` is supplied.
#' @param ar1_parameter Optional estimated AR(1) parameter written as
#'   `"THETA(i)"` or `"SIGMA(i)"`.
#' @param ar1_transform Transform applied to `ar1_parameter`. The default
#'   hyperbolic tangent maps an unconstrained parameter exactly into `(-1, 1)`;
#'   `"identity"` expects the parameter itself to remain strictly inside that
#'   interval.
#' @param blq_method No censoring, Beal M3, or Beal M4.
#' @param lloq Optional scalar LLOQ used when a row-specific `LLOQ` column is
#'   absent.
#' @param iov Number of trailing ETA definitions that are occasion-specific.
#' @param occasion_col Dataset column identifying occasions.
#' @param residual_groups Optional correlated endpoint declarations created by
#'   [nm_residual_group()].
#' @param priors Optional fixed-effect prior definitions.
#' @param mixtures Optional mixture-model definitions.
#' @return A validated serializable likelihood configuration.
#' @export
nm_lik_config <- function(
    error = c("auto", "none", "combined", "additive", "proportional", "exponential", "power", "likelihood"),
    omega = c("diagonal", "full"),
    sigma_corr = c("independent", "ar1"),
    sigma_parameterization = c("sd", "variance"),
    ar1_rho = 0,
    ar1_parameter = NULL,
    ar1_transform = c("tanh", "identity"),
    blq_method = c("none", "m3", "m4"),
    lloq = NA_real_,
    iov = 0L,
    occasion_col = "OCC",
    residual_groups = NULL,
    priors = NULL,
    mixtures = NULL) {
  error <- match.arg(error)
  omega <- match.arg(omega)
  sigma_corr <- match.arg(sigma_corr)
  sigma_parameterization <- match.arg(sigma_parameterization)
  ar1_transform <- match.arg(ar1_transform)
  blq_method <- match.arg(blq_method)
  ar1_rho <- as.numeric(ar1_rho)
  lloq <- as.numeric(lloq)
  iov <- as.integer(iov)
  if (length(ar1_rho) != 1L || !is.finite(ar1_rho) || abs(ar1_rho) >= 1) {
    .nm_stop("`ar1_rho` must be a finite scalar strictly between -1 and 1.")
  }
  ar1_source <- "fixed"
  ar1_index <- 0L
  if (!is.null(ar1_parameter)) {
    ar1_parameter <- toupper(gsub("[[:space:]_]", "", as.character(ar1_parameter)))
    if (length(ar1_parameter) != 1L) {
      .nm_stop("`ar1_parameter` must be one THETA(i) or SIGMA(i) reference.")
    }
    matched <- regmatches(
      ar1_parameter, regexec("^(THETA|SIGMA)\\(([1-9][0-9]*)\\)$", ar1_parameter)
    )[[1L]]
    if (!length(matched)) {
      .nm_stop("`ar1_parameter` must be written as THETA(i) or SIGMA(i).")
    }
    ar1_source <- tolower(matched[[2L]])
    ar1_index <- as.integer(matched[[3L]])
  } else {
    ar1_parameter <- NULL
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
  residual_groups <- .nm_residual_groups(residual_groups)
  if (!is.null(residual_groups) && sigma_corr != "independent") {
    .nm_stop("Cross-endpoint residual groups currently require `sigma_corr = 'independent'`.")
  }
  if (!is.null(residual_groups) && blq_method != "none") {
    .nm_stop("Censored/BLQ residual likelihoods cannot currently be combined with correlated endpoint groups.")
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
      version = 2L, error = error, omega = omega,
      sigma_corr = sigma_corr, sigma_parameterization = sigma_parameterization,
      ar1_rho = ar1_rho, ar1_parameter = ar1_parameter,
      ar1_source = ar1_source, ar1_index = ar1_index,
      ar1_transform = ar1_transform,
      blq_method = blq_method, lloq = lloq, iov = iov,
      occasion_col = as.character(occasion_col), residual_groups = residual_groups,
      priors = priors,
      mixtures = mixtures
    ),
    class = "nm_lik_config"
  )
}

.nm_ar1_rho <- function(model, theta = model$THETAS$Value,
                        sigma = model$SIGMAS$Value) {
  config <- model$LIK_CONFIG
  source <- config$ar1_source %||% "fixed"
  if (identical(source, "fixed")) return(as.numeric(config$ar1_rho))
  values <- if (identical(source, "theta")) theta else sigma
  index <- as.integer(config$ar1_index)
  if (index < 1L || index > length(values)) {
    .nm_stop("Estimated AR(1) parameter index is outside the ", toupper(source), " table.")
  }
  value <- as.numeric(values[[index]])
  rho <- if (identical(config$ar1_transform %||% "tanh", "tanh")) tanh(value) else value
  if (!is.finite(rho) || abs(rho) >= 1) {
    .nm_stop("The estimated AR(1) correlation is not strictly between -1 and 1.")
  }
  rho
}

.nm_lik_config <- function(config, error_type, iov) {
  if (is.null(config)) config <- nm_lik_config(error = error_type, iov = iov)
  if (!inherits(config, "nm_lik_config")) {
    if (!is.list(config)) .nm_stop("LIK_CONFIG must be created by `nm_lik_config()` or be a named list.")
    config$version <- NULL
    config$ar1_source <- NULL
    config$ar1_index <- NULL
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
  if (!identical(x$ar1_source %||% "fixed", "fixed")) {
    cat("  AR(1) parameter:", x$ar1_parameter,
        " transform:", x$ar1_transform, "\n")
  }
  invisible(x)
}
