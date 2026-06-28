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
    val <- subj[[cv]][1L]
    if (is.na(val)) {
      val <- 0
    }
    env[[cv]] <- val
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
    lapply(as.character(covs), function(cv) subj[[cv]][1L]),
    covs
  )
}
