#' Resolve the baseline value of a covariate for one subject.
#'
#' The PK/PRED code is evaluated once per subject, so only a single (baseline)
#' covariate value can be used. This helper makes the two previously-silent
#' failure modes explicit:
#'   * a missing (\code{NA}) baseline value is no longer silently replaced with
#'     0 (which quietly corrupts estimates); we carry forward the first
#'     non-missing value with a warning, and error only if every value is
#'     missing;
#'   * a covariate that changes within a subject (time-varying) triggers a
#'     one-time warning that only the baseline value is used.
#' @keywords internal
.nm_cov_baseline_value <- function(values, cv) {
  values <- as.numeric(values)
  if (length(values) == 0L) {
    .nm_stop("Covariate '", cv, "' has no values for a subject.")
  }
  non_na <- values[!is.na(values)]
  if (length(non_na) == 0L) {
    .nm_stop(
      "Covariate '", cv, "' is missing (all NA) for a subject; ",
      "impute or remove this subject before fitting."
    )
  }
  if (length(non_na) > 1L && any(non_na[-1L] != non_na[1L])) {
    .nm_cov_warn_once(
      paste0("cov_tv_", cv),
      sprintf(
        paste0(
          "Covariate '%s' varies within subject(s); the current engine uses ",
          "the baseline value only. Time-varying covariate effects are not ",
          "yet applied per time segment."
        ),
        cv
      )
    )
  }
  val <- values[1L]
  if (is.na(val)) {
    .nm_cov_warn_once(
      paste0("cov_na_", cv),
      sprintf(
        "Covariate '%s' baseline value is NA; carrying forward first non-missing value.",
        cv
      )
    )
    val <- non_na[1L]
  }
  val
}

#' @keywords internal
.nm_cov_warn_once <- function(key, msg) {
  store <- getOption("LibeRation.cov_warned", NULL)
  if (is.null(store)) {
    store <- new.env(parent = emptyenv())
    options(LibeRation.cov_warned = store)
  }
  if (isTRUE(store[[key]])) {
    return(invisible(FALSE))
  }
  store[[key]] <- TRUE
  warning(msg, call. = FALSE)
  invisible(TRUE)
}

#' Bind covariate columns into a PRED evaluation environment
#'
#' @param model An \code{nm_model} object.
#' @param subj Subject event data (one ID).
#' @param env Environment to augment (typically from \code{.nm_make_env}).
#' @keywords internal
.nm_bind_covariates <- function(model, subj, env) {
  covs <- model$COVARIATES
  if (is.null(covs) || length(covs) == 0L) {
    return(invisible(env))
  }
  covs <- as.character(covs)
  missing <- setdiff(covs, names(subj))
  if (length(missing) > 0L) {
    .nm_stop("Covariate columns missing from data: ", paste(missing, collapse = ", "))
  }
  for (cv in covs) {
    env[[cv]] <- .nm_cov_baseline_value(subj[[cv]], cv)
  }
  invisible(env)
}

#' @keywords internal
.nm_covariate_values <- function(model, dat, id) {
  subj <- .nm_subject_slice(dat, id)
  covs <- model$COVARIATES
  if (is.null(covs) || length(covs) == 0L) {
    return(list())
  }
  stats::setNames(
    lapply(as.character(covs), function(cv) .nm_cov_baseline_value(subj[[cv]], cv)),
    covs
  )
}
