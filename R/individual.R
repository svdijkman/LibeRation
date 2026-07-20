#' Fit one patient's individual random effects
#'
#' Estimates a patient's ETA vector while holding the population model and its
#' THETA, OMEGA, and SIGMA values fixed. The conditional objective and its
#' derivatives are evaluated by the same persistent C++/CppAD engine used by
#' population estimation. A non-zero Gaussian prior permits sequential and
#' dynamic-state updating without treating an earlier point estimate as exact.
#'
#' @param model An [nm_model()] or compiled [NMEngine].
#' @param data A one-subject NONMEM-style event dataset.
#' @param theta,sigma,omega Fixed population parameter values.
#' @param start Initial ETA vector.
#' @param prior_mean Optional Gaussian prior mean for the expanded ETA vector.
#' @param prior_covariance Optional Gaussian prior covariance. Supplying either
#'   custom-prior argument requires both.
#' @param maxit Maximum conditional-mode iterations.
#' @param tolerance Relative optimizer tolerance.
#' @param interaction Include ETA-residual interaction.
#' @param exact_hessian Use the exact conditional Hessian.
#' @return An `nm_individual_fit` containing the MAP ETA state, Laplace
#'   covariance, predictions, diagnostics, and the fixed population values.
#' @export
nm_individual_fit <- function(model, data, theta = NULL, sigma = NULL,
                              omega = NULL, start = NULL,
                              prior_mean = NULL, prior_covariance = NULL,
                              maxit = 100L, tolerance = 1e-7,
                              interaction = TRUE, exact_hessian = TRUE) {
  started <- proc.time()[["elapsed"]]
  engine <- if (inherits(model, "NMEngine")) model else nm_compile(model)
  specification <- engine$model
  theta <- as.numeric(theta %||% specification$THETAS$Value)
  sigma <- as.numeric(sigma %||% specification$SIGMAS$Value)
  omega <- as.numeric(omega %||% specification$OMEGAS$Value)
  normalized <- .nm_engine_data(specification, data)
  if (length(unique(normalized$.ID_INDEX)) != 1L) {
    .nm_stop("`nm_individual_fit()` requires data for exactly one patient.")
  }
  if (!any(normalized$EVID == 0L & normalized$MDV == 0L & is.finite(normalized$DV))) {
    .nm_stop("Individual fitting requires at least one observed DV record.")
  }
  n_eta <- .nm_eta_columns(specification, normalized)
  if (!is.null(prior_mean) != !is.null(prior_covariance)) {
    .nm_stop("Supply both `prior_mean` and `prior_covariance`, or neither.")
  }
  if (is.null(prior_mean)) prior_mean <- rep(0, n_eta)
  prior_mean <- as.numeric(prior_mean)
  if (length(prior_mean) != n_eta || anyNA(prior_mean) || any(!is.finite(prior_mean))) {
    .nm_stop("`prior_mean` must contain one finite value per expanded ETA.")
  }
  base_covariance <- .nm_effect_covariance(specification, normalized, omega)
  custom_prior <- !is.null(prior_covariance)
  if (is.null(prior_covariance)) prior_covariance <- base_covariance
  prior_covariance <- as.matrix(prior_covariance)
  if (!identical(dim(prior_covariance), c(n_eta, n_eta)) ||
      anyNA(prior_covariance) || any(!is.finite(prior_covariance))) {
    .nm_stop("`prior_covariance` must be a finite square matrix matching the expanded ETA vector.")
  }
  prior_covariance <- .nm_positive_definite(
    prior_covariance, "Individual ETA prior covariance"
  )$matrix
  start <- as.numeric(start %||% prior_mean)
  if (length(start) != n_eta || anyNA(start) || any(!is.finite(start))) {
    .nm_stop("`start` must contain one finite value per expanded ETA.")
  }
  evaluator <- .NMSubjectEvaluator$new(
    engine, normalized, theta, sigma, omega, n_eta = n_eta
  )

  if (!custom_prior && all(prior_mean == 0)) {
    mode <- evaluator$eta_mode(
      theta, sigma, omega, start = start, maxit = maxit,
      tolerance = tolerance, interaction = interaction,
      exact_hessian = exact_hessian
    )
  } else {
    base_precision <- solve(.nm_positive_definite(
      base_covariance, "Population ETA covariance"
    )$matrix)
    prior_precision <- solve(prior_covariance)
    eta_positions <- length(theta) + seq_len(n_eta)
    adjusted <- function(eta, gradient = FALSE, hessian = FALSE) {
      raw <- evaluator$objective(
        theta, eta, sigma, omega, gradient = gradient,
        hessian = hessian, interaction = interaction
      )
      centered <- eta - prior_mean
      raw$value <- as.numeric(
        raw$value - crossprod(eta, base_precision %*% eta) +
          crossprod(centered, prior_precision %*% centered)
      )
      if (gradient) {
        raw$gradient[eta_positions] <- raw$gradient[eta_positions] -
          2 * as.numeric(base_precision %*% eta) +
          2 * as.numeric(prior_precision %*% centered)
      }
      if (hessian) {
        raw$hessian[eta_positions, eta_positions] <-
          raw$hessian[eta_positions, eta_positions, drop = FALSE] -
          2 * base_precision + 2 * prior_precision
      }
      raw
    }
    objective <- function(eta) adjusted(eta, gradient = FALSE)$value
    gradient <- function(eta) {
      adjusted(eta, gradient = TRUE)$gradient[eta_positions]
    }
    optimized <- stats::optim(
      start, objective, gradient, method = "BFGS",
      control = list(maxit = as.integer(maxit), reltol = as.numeric(tolerance))
    )
    at_mode <- adjusted(
      optimized$par, gradient = TRUE, hessian = isTRUE(exact_hessian)
    )
    curvature <- if (isTRUE(exact_hessian)) {
      .nm_positive_definite(
        at_mode$hessian[eta_positions, eta_positions, drop = FALSE],
        "Individual ETA posterior curvature"
      )
    } else list(matrix = matrix(numeric(), 0L, 0L), jitter = 0)
    mode <- list(
      par = optimized$par, value = at_mode$value,
      convergence = optimized$convergence,
      hessian = curvature$matrix, jitter = curvature$jitter,
      gradient = at_mode$gradient[eta_positions],
      iterations = as.integer(optimized$counts[["gradient"]]),
      evaluations = as.integer(optimized$counts[["function"]]),
      backend = "r-bfgs-exact-cppad"
    )
  }

  posterior_covariance <- if (n_eta && length(mode$hessian)) {
    # The engine objective is -2 log posterior, so the Laplace covariance is
    # the inverse Hessian of objective/2.
    2 * solve(mode$hessian)
  } else matrix(numeric(), n_eta, n_eta)
  iov <- specification$LIK_CONFIG$iov
  eta_names <- if (iov > 0L) {
    between <- specification$n_eta - iov
    occasions <- (n_eta - between) / iov
    c(
      if (between) paste0("ETA", seq_len(between)) else character(),
      unlist(lapply(seq_len(occasions), function(occasion) {
        paste0("ETA", between + seq_len(iov), "_OCC", occasion)
      }), use.names = FALSE)
    )
  } else paste0("ETA", seq_len(n_eta))
  names(mode$par) <- eta_names
  dimnames(posterior_covariance) <- list(eta_names, eta_names)
  predictions <- engine$simulate(
    normalized, theta = theta, eta = matrix(mode$par, 1L, n_eta), sigma = sigma
  )
  structure(list(
    eta = mode$par,
    eta_covariance = posterior_covariance,
    eta_sd = stats::setNames(sqrt(pmax(diag(posterior_covariance), 0)), eta_names),
    objective = as.numeric(mode$value), convergence = as.integer(mode$convergence),
    gradient = as.numeric(mode$gradient %||% numeric()),
    hessian = mode$hessian, predictions = predictions,
    theta = theta, sigma = sigma, omega = omega,
    prior_mean = stats::setNames(prior_mean, eta_names),
    prior_covariance = prior_covariance,
    model = specification, data = normalized,
    diagnostics = list(
      backend = mode$backend %||% "cpp", iterations = mode$iterations %||% NA_integer_,
      evaluations = mode$evaluations %||% NA_integer_, jitter = mode$jitter %||% 0,
      interaction = isTRUE(interaction), custom_prior = custom_prior,
      elapsed_seconds = unname(proc.time()[["elapsed"]] - started)
    )
  ), class = "nm_individual_fit")
}

#' @export
print.nm_individual_fit <- function(x, ...) {
  cat("LibeRation individual fit\n")
  cat("  ETA:", paste(names(x$eta), format(x$eta, digits = 5), collapse = ", "), "\n")
  cat("  objective:", format(x$objective), " convergence:", x$convergence, "\n")
  invisible(x)
}
