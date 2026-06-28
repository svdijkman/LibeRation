#' Profiling counters for estimation runs
#' @keywords internal
.nm_profile_env <- function() {
  if (is.null(.nm_state$profile)) {
    .nm_state$profile <- list()
  }
  .nm_state$profile
}

#' @keywords internal
.nm_profile_reset <- function() {
  .nm_state$profile <- list()
  invisible(NULL)
}

#' @keywords internal
.nm_profile_add <- function(label, time_sec, count = 1L) {
  pe <- .nm_profile_env()
  if (is.null(pe[[label]])) {
    pe[[label]] <- list(time = 0, count = 0L)
  }
  pe[[label]]$time <- pe[[label]]$time + time_sec
  pe[[label]]$count <- pe[[label]]$count + as.integer(count)
  .nm_state$profile <- pe
  invisible(NULL)
}

#' @keywords internal
.nm_profile_snapshot <- function() {
  pe <- .nm_profile_env()
  if (length(pe) == 0L) {
    return(NULL)
  }
  data.frame(
    label = names(pe),
    time_sec = vapply(pe, function(x) x$time, numeric(1)),
    count = vapply(pe, function(x) x$count, integer(1)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' @rdname profile.nm_fit
#' @method print nm_profile
#' @param x A profiling summary data frame.
#' @param ... Unused.
#' @examples
#' \dontrun{
#' options(LibeRation.profile = TRUE)
#' sim <- nm_synthetic_theo(n_sub = 2L)
#' fit <- nm_est(sim$model, sim$data, method = "FO",
#'               control = list(maxit = 3L, compute_inference = FALSE))
#' options(LibeRation.profile = FALSE)
#' if (!is.null(fit$profile)) print(profile(fit))
#' }
#' @export
print.nm_profile <- function(x, ...) {
  cat("LibeRation profiling summary\n")
  print(x[order(-x$time_sec), , drop = FALSE])
  invisible(x)
}

#' Summarise profiling data attached to a fit
#' @param object An \code{nm_fit} object.
#' @examples
#' \dontrun{
#' options(LibeRation.profile = TRUE)
#' sim <- nm_synthetic_theo(n_sub = 2L)
#' fit <- nm_est(sim$model, sim$data, method = "FO",
#'               control = list(maxit = 3L, compute_inference = FALSE))
#' options(LibeRation.profile = FALSE)
#' profile(fit)
#' }
#' @export
profile.nm_fit <- function(object, ...) {
  if (is.null(object$profile)) {
    return(invisible(NULL))
  }
  structure(object$profile, class = "nm_profile")
}
