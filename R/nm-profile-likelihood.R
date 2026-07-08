#' Likelihood profiling of a THETA parameter
#'
#' Profiles the objective function value (OFV, -2 log-likelihood) over a grid of
#' values for a chosen \code{THETA}, re-optimising all the other parameters at
#' each grid point (the profiled parameter is held fixed). Returns the profile
#' and a likelihood-based confidence interval defined by the change in OFV
#' (\eqn{\Delta OFV = 3.84} for a 95\% CI, the \eqn{\chi^2_1} quantile).
#'
#' This is likelihood-based \emph{parameter} profiling. For the runtime timing
#' profiler see \code{\link{nm_time_profile}}.
#'
#' @param fit A converged \code{nm_fit} object.
#' @param parameter THETA to profile: an integer THETA index, or a label such as
#'   \code{"THETA2"}.
#' @param grid Optional numeric vector of values to profile over. When
#'   \code{NULL} a symmetric grid of \code{n} points spanning
#'   \code{estimate * (1 +/- span)} is used.
#' @param n Number of grid points when \code{grid} is \code{NULL}.
#' @param span Fractional half-width of the default grid around the estimate.
#' @param delta_ofv Change in OFV defining the CI (default 3.84 = 95\%).
#' @param method Estimation method for the re-optimisations (default: the
#'   method of \code{fit}).
#' @param control Control list passed to \code{\link{nm_est}}
#'   (\code{compute_inference = FALSE} is forced).
#' @param ... Passed to \code{\link{nm_est}}.
#' @return A list with class \code{nm_profile_likelihood}: \code{profile}
#'   (data.frame of \code{value}, \code{ofv}, \code{dofv}), \code{estimate},
#'   \code{ci} (named \code{lower}/\code{upper}), \code{delta_ofv}, and
#'   \code{parameter}.
#' @examples
#' \dontrun{
#' sim <- nm_synthetic_theo(n_sub = 20L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "FO")
#' pl <- nm_profile_likelihood(fit, "THETA1")
#' pl$ci
#' }
#' @export
nm_profile_likelihood <- function(fit,
                                  parameter,
                                  grid = NULL,
                                  n = 7L,
                                  span = 0.3,
                                  delta_ofv = 3.84,
                                  method = NULL,
                                  control = list(),
                                  ...) {
  if (is.null(fit) || !inherits(fit, "nm_fit") || is.null(fit$model)) {
    .nm_stop("fit must be a fitted nm_fit object with a model.")
  }
  model <- fit$model
  data <- fit$data
  if (is.null(data)) {
    .nm_stop("fit must carry its data (fit$data) for profiling.")
  }
  idx <- .nm_profile_ll_theta_index(model, parameter)
  label <- paste0("THETA", model$THETAS$THETA[idx])
  est_val <- as.numeric(fit$theta[idx])
  method <- method %||% fit$method
  if (identical(method, "POSTHOC")) {
    .nm_stop("Cannot profile a POSTHOC fit; re-fit with an estimation method.")
  }
  base_ofv <- as.numeric(fit$objective)
  if (is.null(grid)) {
    half <- max(abs(est_val) * span, .Machine$double.eps)
    grid <- seq(est_val - half, est_val + half, length.out = as.integer(n))
  }
  grid <- sort(unique(as.numeric(grid)))

  p_fit <- .nm_pack(model, fit$theta, fit$omega, fit$sigma)
  ofv <- rep(NA_real_, length(grid))
  for (k in seq_along(grid)) {
    m2 <- model
    if (!"FIX" %in% names(m2$THETAS)) {
      m2$THETAS$FIX <- FALSE
    }
    m2$THETAS$FIX[idx] <- TRUE
    m2$THETAS$Value[idx] <- grid[k]
    start <- p_fit
    start[idx] <- grid[k]
    fk <- tryCatch(
      nm_est(
        m2, data, method = method, start = start, grad = "numeric",
        control = c(list(compute_inference = FALSE), control), ...
      ),
      error = function(e) NULL
    )
    if (!is.null(fk) && is.finite(fk$objective)) {
      ofv[k] <- as.numeric(fk$objective)
    }
  }
  ref_ofv <- min(c(base_ofv, ofv), na.rm = TRUE)
  profile <- data.frame(
    parameter = label,
    value = grid,
    ofv = ofv,
    dofv = ofv - ref_ofv,
    stringsAsFactors = FALSE
  )
  ci <- .nm_profile_ll_ci(grid, ofv, est_val, ref_ofv + delta_ofv)
  structure(
    list(
      parameter = label,
      estimate = est_val,
      profile = profile,
      ci = ci,
      delta_ofv = delta_ofv,
      ref_ofv = ref_ofv,
      method = method
    ),
    class = "nm_profile_likelihood"
  )
}

#' @keywords internal
.nm_profile_ll_theta_index <- function(model, parameter) {
  n_th <- nrow(model$THETAS)
  if (is.numeric(parameter)) {
    i <- as.integer(parameter[1L])
    if (i >= 1L && i <= n_th) {
      return(i)
    }
    .nm_stop("THETA index out of range: ", parameter)
  }
  p <- toupper(trimws(as.character(parameter[1L])))
  labs <- toupper(paste0("THETA", model$THETAS$THETA))
  hit <- match(p, labs)
  if (is.na(hit)) {
    hit <- match(sub("^THETA", "", p), as.character(model$THETAS$THETA))
  }
  if (is.na(hit)) {
    .nm_stop(
      "Could not resolve parameter '", parameter,
      "'. Use a THETA index or a label like 'THETA2'."
    )
  }
  hit
}

#' Linear-interpolation likelihood CI from a profile.
#' @keywords internal
.nm_profile_ll_ci <- function(grid, ofv, est, target) {
  cross <- function(xs, ys) {
    # find x where y crosses `target`, interpolating between bracketing points
    out <- NA_real_
    for (i in seq_len(length(xs) - 1L)) {
      y0 <- ys[i]
      y1 <- ys[i + 1L]
      if (is.na(y0) || is.na(y1)) {
        next
      }
      if ((y0 - target) * (y1 - target) <= 0 && y0 != y1) {
        frac <- (target - y0) / (y1 - y0)
        out <- xs[i] + frac * (xs[i + 1L] - xs[i])
      }
    }
    out
  }
  lower_mask <- grid <= est
  upper_mask <- grid >= est
  lower <- if (sum(lower_mask) >= 2L) {
    cross(rev(grid[lower_mask]), rev(ofv[lower_mask]))
  } else {
    NA_real_
  }
  upper <- if (sum(upper_mask) >= 2L) {
    cross(grid[upper_mask], ofv[upper_mask])
  } else {
    NA_real_
  }
  c(lower = lower, upper = upper)
}

#' @rdname nm_profile_likelihood
#' @method print nm_profile_likelihood
#' @param x An \code{nm_profile_likelihood} object.
#' @param ... Unused.
#' @export
print.nm_profile_likelihood <- function(x, ...) {
  cat("Likelihood profile:", x$parameter, "\n")
  cat("  estimate:", format(x$estimate, digits = 6), "\n")
  cat("  95% CI (dOFV =", x$delta_ofv, "): [",
      format(x$ci[["lower"]], digits = 6), ",",
      format(x$ci[["upper"]], digits = 6), "]\n")
  print(x$profile)
  invisible(x)
}
