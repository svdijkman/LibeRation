#' Configure a Gaussian state-space model
#'
#' Defines a continuous-observation linear or nonlinear Gaussian state-space
#' likelihood. Linear models use an exact Kalman filter. Nonlinear transition
#' and observation assignments can use an extended Kalman filter (EKF),
#' unscented Kalman filter (UKF), or a seeded bootstrap particle filter. All
#' deterministic filter operations remain in the compiled C++ objective.
#' Every name refers to an assignment in `$ERROR`, so transition, process,
#' observation, and covariance terms may depend on `THETA`, `ETA`, covariates,
#' `DT`, and ordinary `$PK/$PRED` outputs. The compiled C++ engine runs the
#' Kalman filter as part of the exact CppAD population objective and exposes a
#' Rauch--Tung--Striebel retrospective smoother.
#'
#' The measurement equation is `DV = F + H x + epsilon` when `baseline` is
#' `"prediction"`, and `DV = H x + epsilon` when it is `"zero"`.
#'
#' @param states Unique labels for latent continuous states.
#' @param initial_mean Character vector naming the initial-state means.
#' @param initial_covariance Square symmetric character matrix naming the
#'   initial covariance entries.
#' @param transition For `filter = "linear"`, a square character matrix naming
#'   the discrete transition matrix entries. For a nonlinear filter, a
#'   character vector naming the next-state function for every state.
#' @param process_covariance Square symmetric character matrix naming process
#'   covariance entries for the current interval.
#' @param observation For `filter = "linear"`, a character vector naming the
#'   observation-loading row `H`. For a nonlinear filter, one assignment naming
#'   the scalar observation function.
#' @param observation_variance One assignment naming the positive measurement
#'   variance.
#' @param baseline Add the structural prediction `F` to the latent-state
#'   observation mean, or use a zero baseline.
#' @param by_dvid Maintain a separate latent trajectory for each DVID within
#'   subject.
#' @param filter State inference algorithm: exact linear Kalman, EKF, UKF, or
#'   seeded bootstrap particle filtering.
#' @param state_inputs Input symbols used for latent states by nonlinear
#'   transition and observation assignments. The default is `STATE_<name>`.
#' @param jacobian_step Relative central-difference step used for EKF state
#'   Jacobians. Parameter derivatives are still propagated by CppAD.
#' @param ukf_alpha,ukf_beta,ukf_kappa Scaled unscented-transform controls.
#' @param particles Number of particles for the bootstrap filter.
#' @param ess_threshold Resample when effective sample size falls below this
#'   proportion of `particles`.
#' @param seed Reproducible common-random-number seed for particle inference.
#' @param dynamics Discrete transition assignments or continuous-discrete SDE
#'   drift/diffusion assignments. Prefer [nm_sde_config()] for SDE models.
#' @param sde_method Fixed-step Euler--Maruyama or diagonal Milstein
#'   propagation for SDE dynamics.
#' @param sde_substeps Number of fixed SDE steps per observation interval.
#' @return A validated serializable linear Gaussian state-space declaration.
#' @export
nm_kalman_config <- function(
    states, initial_mean, initial_covariance, transition,
    process_covariance, observation, observation_variance,
    baseline = c("prediction", "zero"), by_dvid = TRUE,
    filter = c("linear", "ekf", "ukf", "particle"), state_inputs = NULL,
    jacobian_step = 1e-5, ukf_alpha = 0.5, ukf_beta = 2,
    ukf_kappa = 0, particles = 256L, ess_threshold = 0.5,
    seed = 20260721L, dynamics = c("discrete", "sde"),
    sde_method = c("euler", "milstein"), sde_substeps = 8L) {
  states <- trimws(as.character(states))
  if (!length(states) || anyNA(states) || any(!nzchar(states)) || anyDuplicated(states)) {
    .nm_stop("`states` must contain unique, non-empty latent-state labels.")
  }
  names_vector <- function(value, expected, label) {
    value <- trimws(as.character(value))
    if (length(value) != expected || anyNA(value) || any(!nzchar(value))) {
      .nm_stop("`", label, "` must name exactly ", expected, " non-empty $ERROR assignment(s).")
    }
    unname(value)
  }
  names_matrix <- function(value, label, symmetric = FALSE) {
    value <- as.matrix(value)
    if (!identical(dim(value), rep(length(states), 2L))) {
      .nm_stop("`", label, "` must be square with one row and column per state.")
    }
    value <- matrix(
      names_vector(as.vector(value), length(states)^2L, label),
      length(states), length(states), dimnames = list(states, states)
    )
    if (symmetric && !identical(unname(value), unname(t(value)))) {
      .nm_stop("`", label, "` must be symmetric and reference matching assignments.")
    }
    value
  }
  filter <- match.arg(filter)
  dynamics <- match.arg(dynamics)
  sde_method <- match.arg(sde_method)
  if (dynamics == "sde" && filter == "linear") {
    .nm_stop("SDE dynamics require `filter = 'ekf'`, `'ukf'`, or `'particle'`.")
  }
  initial_mean <- names_vector(initial_mean, length(states), "initial_mean")
  if (identical(filter, "linear")) {
    observation <- names_vector(observation, length(states), "observation")
  } else {
    transition <- names_vector(transition, length(states), "transition")
    observation <- names_vector(observation, 1L, "observation")
  }
  observation_variance <- names_vector(
    observation_variance, 1L, "observation_variance"
  )
  initial_covariance <- names_matrix(
    initial_covariance, "initial_covariance", symmetric = TRUE
  )
  if (identical(filter, "linear")) {
    transition <- names_matrix(transition, "transition")
  }
  process_covariance <- names_matrix(
    process_covariance, "process_covariance", symmetric = dynamics != "sde"
  )
  baseline <- match.arg(baseline)
  if (length(by_dvid) != 1L || is.na(by_dvid)) {
    .nm_stop("`by_dvid` must be TRUE or FALSE.")
  }
  if (is.null(state_inputs)) state_inputs <- paste0("STATE_", make.names(states, unique = TRUE))
  state_inputs <- trimws(as.character(state_inputs))
  if (length(state_inputs) != length(states) || anyNA(state_inputs) ||
      any(!nzchar(state_inputs)) || anyDuplicated(state_inputs)) {
    .nm_stop("`state_inputs` must contain one unique, non-empty symbol per state.")
  }
  if (identical(filter, "linear")) state_inputs <- character()
  jacobian_step <- as.numeric(jacobian_step)
  ukf_alpha <- as.numeric(ukf_alpha)
  ukf_beta <- as.numeric(ukf_beta)
  ukf_kappa <- as.numeric(ukf_kappa)
  particles <- as.integer(particles)
  ess_threshold <- as.numeric(ess_threshold)
  seed <- as.integer(seed)
  sde_substeps <- as.integer(sde_substeps)
  if (length(jacobian_step) != 1L || !is.finite(jacobian_step) || jacobian_step <= 0 ||
      length(ukf_alpha) != 1L || !is.finite(ukf_alpha) || ukf_alpha <= 0 ||
      length(ukf_beta) != 1L || !is.finite(ukf_beta) ||
      length(ukf_kappa) != 1L || !is.finite(ukf_kappa) ||
      length(particles) != 1L || is.na(particles) || particles < 16L ||
      length(ess_threshold) != 1L || !is.finite(ess_threshold) ||
      ess_threshold <= 0 || ess_threshold > 1 ||
      length(seed) != 1L || is.na(seed) || length(sde_substeps) != 1L ||
      is.na(sde_substeps) || sde_substeps < 1L || sde_substeps > 10000L) {
    .nm_stop("Invalid EKF/UKF/particle filter control value.")
  }
  if (filter == "ukf" && length(states) + ukf_alpha^2 *
      (length(states) + ukf_kappa) - length(states) <= 0) {
    .nm_stop("UKF controls require `length(states) + lambda > 0`.")
  }
  structure(
    list(
      version = 2L, states = states, initial_mean = initial_mean,
      initial_covariance = initial_covariance, transition = transition,
      process_covariance = process_covariance, observation = observation,
      observation_variance = observation_variance, baseline = baseline,
      by_dvid = isTRUE(by_dvid), filter = filter,
      state_inputs = state_inputs, jacobian_step = jacobian_step,
      ukf_alpha = ukf_alpha, ukf_beta = ukf_beta, ukf_kappa = ukf_kappa,
      particles = particles, ess_threshold = ess_threshold, seed = seed,
      dynamics = dynamics, sde_method = sde_method,
      sde_substeps = sde_substeps
    ),
    class = "nm_kalman_config"
  )
}

#' Configure a continuous-discrete stochastic differential equation
#'
#' Creates a compiled state-space declaration for an Ito SDE. `drift` names
#' one `$ERROR` assignment per state and `diffusion` names the square diffusion
#' loading matrix. The C++ engine discretizes each observation interval and
#' applies EKF, UKF, or particle inference. Milstein is available for
#' componentwise/diagonal diffusion; Euler--Maruyama is general.
#'
#' @inheritParams nm_kalman_config
#' @param drift Character vector naming the SDE drift assignments.
#' @param diffusion Square character matrix naming the diffusion loading.
#' @param method Euler--Maruyama or Milstein propagation.
#' @param substeps Fixed substeps per observation interval.
#' @export
nm_sde_config <- function(
    states, initial_mean, initial_covariance, drift, diffusion,
    observation, observation_variance,
    baseline = c("prediction", "zero"), by_dvid = TRUE,
    filter = c("ukf", "ekf", "particle"), state_inputs = NULL,
    method = c("euler", "milstein"), substeps = 8L,
    jacobian_step = 1e-5, ukf_alpha = 0.5, ukf_beta = 2,
    ukf_kappa = 0, particles = 256L, ess_threshold = 0.5,
    seed = 20260721L) {
  config <- nm_kalman_config(
    states = states, initial_mean = initial_mean,
    initial_covariance = initial_covariance, transition = drift,
    process_covariance = diffusion, observation = observation,
    observation_variance = observation_variance,
    baseline = match.arg(baseline), by_dvid = by_dvid,
    filter = match.arg(filter), state_inputs = state_inputs,
    jacobian_step = jacobian_step, ukf_alpha = ukf_alpha,
    ukf_beta = ukf_beta, ukf_kappa = ukf_kappa,
    particles = particles, ess_threshold = ess_threshold, seed = seed,
    dynamics = "sde", sde_method = match.arg(method),
    sde_substeps = substeps
  )
  class(config) <- c("nm_sde_config", class(config))
  config
}

#' Configure an ARMA residual process
#'
#' Generates an exact linear state-space representation for residual
#' autocorrelation around the structural prediction. Coefficients and
#' variances are expressions, so they may reference THETA/SIGMA assignments.
#' The generated `$ERROR` component is inserted automatically by [nm_model()].
#'
#' @param ar Character expressions for AR coefficients.
#' @param ma Character expressions for MA coefficients.
#' @param innovation_variance Expression for innovation variance.
#' @param observation_variance Additional independent measurement variance;
#'   zero is valid.
#' @param initial_variance Diffuse initial variance expression.
#' @param prefix Prefix for generated assignment names.
#' @param by_dvid Maintain independent ARMA histories by DVID.
#' @return An `nm_kalman_config` carrying generated `$ERROR` assignments.
#' @export
nm_arma_config <- function(ar = character(), ma = character(),
                           innovation_variance = "SIGMA(1)",
                           observation_variance = "0",
                           initial_variance = "1000",
                           prefix = "ARMA", by_dvid = TRUE) {
  ar <- trimws(as.character(ar)); ma <- trimws(as.character(ma))
  if (anyNA(c(ar, ma)) || any(!nzchar(c(ar, ma)))) {
    .nm_stop("AR and MA coefficients must be non-empty expressions.")
  }
  prefix <- gsub("[^A-Za-z0-9_]", "_", toupper(as.character(prefix)[[1L]]))
  if (!nzchar(prefix)) prefix <- "ARMA"
  p <- length(ar); q <- length(ma)
  y_count <- max(1L, p)
  dimension <- y_count + q
  states <- c(paste0("residual_lag", seq_len(y_count)),
              if (q) paste0("innovation_lag", seq_len(q)) else character())
  zero <- paste0(prefix, "_ZERO")
  initial_mean_names <- paste0(prefix, "_M", seq_len(dimension))
  initial_covariance_names <- matrix(zero, dimension, dimension)
  diag(initial_covariance_names) <- paste0(prefix, "_P", seq_len(dimension))
  transition_names <- matrix(zero, dimension, dimension)
  transition_names[1L, seq_len(p)] <- if (p) paste0(prefix, "_AR", seq_len(p)) else character()
  if (q) transition_names[1L, y_count + seq_len(q)] <- paste0(prefix, "_MA", seq_len(q))
  if (y_count > 1L) {
    for (index in 2:y_count) transition_names[index, index - 1L] <- paste0(prefix, "_ONE")
  }
  if (q > 1L) {
    for (index in 2:q) {
      transition_names[y_count + index, y_count + index - 1L] <- paste0(prefix, "_ONE")
    }
  }
  loading <- numeric(dimension); loading[[1L]] <- 1
  if (q) loading[[y_count + 1L]] <- 1
  process_names <- matrix(zero, dimension, dimension)
  for (row in seq_len(dimension)) for (column in seq_len(dimension)) {
    if (loading[[row]] && loading[[column]]) {
      process_names[row, column] <- paste0(prefix, "_Q")
    }
  }
  observation_names <- rep(zero, dimension)
  observation_names[[1L]] <- paste0(prefix, "_ONE")
  observation_variance_name <- paste0(prefix, "_R")
  code <- c(
    paste0(zero, " = 0"), paste0(prefix, "_ONE = 1"),
    paste0(initial_mean_names, " = 0"),
    paste0(diag(initial_covariance_names), " = ", initial_variance),
    if (p) paste0(prefix, "_AR", seq_len(p), " = ", ar) else character(),
    if (q) paste0(prefix, "_MA", seq_len(q), " = ", ma) else character(),
    paste0(prefix, "_Q = ", innovation_variance),
    paste0(observation_variance_name, " = ", observation_variance)
  )
  config <- nm_kalman_config(
    states = states, initial_mean = initial_mean_names,
    initial_covariance = initial_covariance_names,
    transition = transition_names, process_covariance = process_names,
    observation = observation_names,
    observation_variance = observation_variance_name,
    baseline = "prediction", by_dvid = by_dvid, filter = "linear"
  )
  attr(config, "generated_error") <- paste(code, collapse = "\n")
  attr(config, "arma_order") <- c(p = p, q = q)
  class(config) <- c("nm_arma_config", class(config))
  config
}

#' @export
print.nm_arma_config <- function(x, ...) {
  order <- attr(x, "arma_order", exact = TRUE) %||% c(p = NA, q = NA)
  cat("LibeRation ARMA(", order[[1L]], ",", order[[2L]],
      ") residual state-space model\n", sep = "")
  invisible(x)
}

#' @export
print.nm_sde_config <- function(x, ...) {
  cat("LibeRation continuous-discrete Ito SDE\n")
  cat("  states:", paste(x$states, collapse = ", "),
      " filter:", x$filter, " method:", x$sde_method,
      " substeps:", x$sde_substeps, "\n")
  invisible(x)
}

.nm_kalman_config <- function(config) {
  if (is.null(config)) return(NULL)
  if (inherits(config, "nm_kalman_config")) return(config)
  if (!is.list(config) || is.null(names(config))) {
    .nm_stop("KALMAN_CONFIG must be created by `nm_kalman_config()` or be a named list.")
  }
  config$version <- NULL
  do.call(nm_kalman_config, config)
}

#' Filter and smooth latent Gaussian states
#'
#' @param object An [nm_model()], compiled [NMEngine], or fitted `nm_fit` with
#'   `KALMAN_CONFIG`.
#' @param data Observation data; defaults to fitted data for an `nm_fit`.
#' @param type Use individual or population ETAs for a fitted model.
#' @param theta,eta,sigma Optional parameters for an unfitted model.
#' @return Observation data augmented with filtered, predicted, and smoothed
#'   latent-state means and standard deviations, innovations, innovation
#'   variances, and row likelihood contributions.
#' @export
nm_kalman_decode <- function(object, data = NULL,
                             type = c("individual", "population"),
                             theta = NULL, eta = NULL, sigma = NULL) {
  type <- match.arg(type)
  fit <- if (inherits(object, "nm_fit")) object else NULL
  if (!is.null(fit)) {
    model <- fit$model
    data <- .nm_engine_data(model, data %||% fit$data)
    theta <- theta %||% fit$theta
    sigma <- sigma %||% fit$sigma
    eta <- eta %||% .nm_fit_eta_for_data(fit, data, type)
    engine <- nm_compile(model)
  } else {
    engine <- if (inherits(object, "NMEngine")) object else nm_compile(object)
    if (is.null(data)) .nm_stop("`data` is required when `object` is not an nm_fit.")
    data <- .nm_engine_data(engine$model, data)
    theta <- theta %||% engine$model$THETAS$Value
    sigma <- sigma %||% engine$model$SIGMAS$Value
  }
  decoded <- engine$kalman_filter(data, theta = theta, eta = eta, sigma = sigma)
  output <- as.data.frame(data)
  state_names <- make.names(engine$model$KALMAN_CONFIG$states, unique = TRUE)
  for (index in seq_along(state_names)) {
    suffix <- state_names[[index]]
    output[[paste0("KF_PRED_", suffix)]] <- decoded$predicted_mean[, index]
    output[[paste0("KF_FILTER_", suffix)]] <- decoded$filtered_mean[, index]
    output[[paste0("KF_SMOOTH_", suffix)]] <- decoded$smoothed_mean[, index]
    output[[paste0("KF_FILTER_SD_", suffix)]] <- sqrt(decoded$filtered_variance[, index])
    output[[paste0("KF_SMOOTH_SD_", suffix)]] <- sqrt(decoded$smoothed_variance[, index])
  }
  output$KF_INNOVATION <- decoded$innovation
  output$KF_INNOVATION_VARIANCE <- decoded$innovation_variance
  output$KF_ROW_NLL <- decoded$row_nll
  if (length(decoded$regimes %||% character())) {
    regime_names <- make.names(decoded$regimes, unique = TRUE)
    for (index in seq_along(regime_names)) {
      output[[paste0("REGIME_FILTER_PROB_", regime_names[[index]])]] <-
        decoded$filtered_regime[, index]
      output[[paste0("REGIME_SMOOTH_PROB_", regime_names[[index]])]] <-
        decoded$smoothed_regime[, index]
    }
    output$REGIME_FILTER <- decoded$regimes[
      max.col(decoded$filtered_regime, ties.method = "first")]
    output$REGIME_SMOOTH <- decoded$regimes[
      max.col(decoded$smoothed_regime, ties.method = "first")]
    attr(output, "regimes") <- decoded$regimes
  }
  attr(output, "log_likelihood") <- decoded$log_likelihood
  attr(output, "states") <- engine$model$KALMAN_CONFIG$states
  attr(output, "eta_type") <- type
  attr(output, "filter") <- decoded$filter %||%
    engine$model$KALMAN_CONFIG$filter %||% "linear"
  attr(output, "smoother") <- decoded$smoother %||% "RTS"
  class(output) <- c("nm_kalman_decode", class(output))
  output
}

.nm_simulate_kalman <- function(engine, result, theta, eta, sigma) {
  states <- length(engine$model$KALMAN_CONFIG$states)
  draw_columns <- states * if (identical(
    engine$model$KALMAN_CONFIG$dynamics %||% "discrete", "sde"
  )) engine$model$KALMAN_CONFIG$sde_substeps %||% 8L else 1L
  process_normals <- matrix(
    stats::rnorm(nrow(result) * draw_columns), nrow(result), draw_columns
  )
  observation_normals <- stats::rnorm(nrow(result))
  result$DV <- .liberation_engine_kalman_simulate(
    engine$pointer, result, as.numeric(theta), as.matrix(eta),
    as.numeric(sigma), process_normals, observation_normals
  )
  result
}

#' @export
print.nm_kalman_config <- function(x, ...) {
  cat("LibeRation Gaussian state-space model\n")
  cat("  states:", paste(x$states, collapse = ", "),
      " filter:", x$filter %||% "linear", " baseline:", x$baseline,
      " by DVID:", x$by_dvid, "\n")
  invisible(x)
}
