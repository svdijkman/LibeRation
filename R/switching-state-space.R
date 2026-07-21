#' Configure an experimental switching nonlinear state-space model
#'
#' A discrete latent regime is propagated jointly with the continuous state by
#' the compiled particle filter. A deterministic stratified regime proposal
#' and differentiable importance weights preserve transition-parameter
#' gradients, while each regime has its own nonlinear
#' transition (or SDE drift/diffusion) and observation model. Regime filtering
#' and genealogical smoothing are returned by [nm_kalman_decode()].
#'
#' @param regimes Unique regime labels.
#' @param initial_regime Assignment names for initial regime weights.
#' @param regime_transition Square matrix of assignment names for regime
#'   transition weights.
#' @param transition List with one next-state/drift assignment vector per regime.
#' @param process_covariance List with one covariance/diffusion matrix per regime.
#' @param observation,observation_variance One assignment per regime.
#' @param initial_scale,transition_scale Probability or log weights.
#' @inheritParams nm_kalman_config
#' @export
nm_switching_state_space_config <- function(
    regimes, initial_regime, regime_transition,
    states, initial_mean, initial_covariance,
    transition, process_covariance, observation, observation_variance,
    initial_scale = c("probability", "log"),
    transition_scale = c("probability", "log"),
    baseline = c("prediction", "zero"), by_dvid = TRUE,
    state_inputs = NULL, particles = 512L, ess_threshold = 0.5,
    seed = 20260721L, dynamics = c("discrete", "sde"),
    sde_method = c("euler", "milstein"), sde_substeps = 8L,
    jacobian_step = 1e-5) {
  regimes <- trimws(as.character(regimes))
  states <- trimws(as.character(states))
  if (length(regimes) < 2L || anyNA(regimes) || any(!nzchar(regimes)) ||
      anyDuplicated(regimes)) {
    .nm_stop("A switching state-space model requires at least two unique regimes.")
  }
  if (!length(states) || anyNA(states) || any(!nzchar(states)) || anyDuplicated(states)) {
    .nm_stop("Continuous states must be unique, non-empty labels.")
  }
  assignments <- function(value, count, label) {
    value <- trimws(as.character(value))
    if (length(value) != count || anyNA(value) || any(!nzchar(value))) {
      .nm_stop("`", label, "` requires exactly ", count, " assignment names.")
    }
    unname(value)
  }
  initial_regime <- assignments(initial_regime, length(regimes), "initial_regime")
  regime_transition <- as.matrix(regime_transition)
  if (!identical(dim(regime_transition), rep(length(regimes), 2L))) {
    .nm_stop("`regime_transition` must be square with one row/column per regime.")
  }
  regime_transition <- matrix(
    assignments(as.vector(regime_transition), length(regimes)^2L, "regime_transition"),
    length(regimes), length(regimes), dimnames = list(regimes, regimes))
  if (!is.list(transition) || length(transition) != length(regimes) ||
      !is.list(process_covariance) || length(process_covariance) != length(regimes)) {
    .nm_stop("Transition and process specifications require one list element per regime.")
  }
  transition <- lapply(transition, assignments, count = length(states), label = "transition")
  process_covariance <- lapply(process_covariance, function(value) {
    value <- as.matrix(value)
    if (!identical(dim(value), rep(length(states), 2L))) {
      .nm_stop("Every process covariance/diffusion matrix must match the state dimension.")
    }
    matrix(assignments(as.vector(value), length(states)^2L, "process_covariance"),
           length(states), length(states), dimnames = list(states, states))
  })
  observation <- assignments(observation, length(regimes), "observation")
  observation_variance <- assignments(
    observation_variance, length(regimes), "observation_variance")
  dynamics <- match.arg(dynamics)
  config <- nm_kalman_config(
    states = states, initial_mean = initial_mean,
    initial_covariance = initial_covariance,
    transition = transition[[1L]], process_covariance = process_covariance[[1L]],
    observation = observation[[1L]], observation_variance = observation_variance[[1L]],
    baseline = match.arg(baseline), by_dvid = by_dvid, filter = "particle",
    state_inputs = state_inputs, jacobian_step = jacobian_step,
    particles = particles, ess_threshold = ess_threshold, seed = seed,
    dynamics = dynamics, sde_method = match.arg(sde_method),
    sde_substeps = sde_substeps
  )
  config$switching <- list(
    version = 1L, regimes = regimes, initial = initial_regime,
    transition = regime_transition, initial_scale = match.arg(initial_scale),
    transition_scale = match.arg(transition_scale),
    state_transition = transition, process = process_covariance,
    observation = observation, observation_variance = observation_variance
  )
  class(config) <- c("nm_switching_state_space_config", class(config))
  config
}

#' @export
print.nm_switching_state_space_config <- function(x, ...) {
  cat("LibeRation experimental switching state-space model\n")
  cat("  regimes:", paste(x$switching$regimes, collapse = ", "),
      " continuous states:", paste(x$states, collapse = ", "), "\n")
  cat("  inference: compiled bootstrap particle filter (", x$particles,
      " particles)\n", sep = "")
  invisible(x)
}
