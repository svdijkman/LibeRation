#' Define parameter priors
#'
#' @param parameter Parameter name such as `THETA1`, `SIGMA1`, or `OMEGA1`.
#' @param distribution Prior family.
#' @param mean,sd Normal/log-normal parameters.
#' @param shape,rate Inverse-gamma parameters.
#' @return A one-row prior table suitable for `nm_lik_config(priors=...)`.
#' @export
nm_prior <- function(parameter,
                     distribution = c("normal", "lognormal", "half_normal", "inverse_gamma"),
                     mean = 0, sd = 1, shape = NA_real_, rate = NA_real_) {
  distribution <- match.arg(distribution)
  parameter <- toupper(gsub("_", "", as.character(parameter)))
  if (length(parameter) != 1L || !grepl("^(THETA|SIGMA|OMEGA)[0-9]+$", parameter)) {
    .nm_stop("Prior parameter must look like THETA1, SIGMA1, or OMEGA1.")
  }
  out <- data.frame(
    parameter = parameter, distribution = distribution,
    mean = as.numeric(mean), sd = as.numeric(sd),
    shape = as.numeric(shape), rate = as.numeric(rate),
    stringsAsFactors = FALSE
  )
  if (distribution %in% c("normal", "lognormal", "half_normal") &&
      (!is.finite(out$sd) || out$sd <= 0 || !is.finite(out$mean))) {
    .nm_stop("Normal-family priors require finite `mean` and positive `sd`.")
  }
  if (distribution == "inverse_gamma" &&
      (!is.finite(out$shape) || out$shape <= 0 || !is.finite(out$rate) || out$rate <= 0)) {
    .nm_stop("Inverse-gamma priors require positive finite `shape` and `rate`.")
  }
  out
}

.nm_parameter_value <- function(parameters, name) {
  index <- as.integer(sub("^[A-Z]+", "", name))
  if (startsWith(name, "THETA")) return(parameters$theta[[index]])
  if (startsWith(name, "SIGMA")) return(parameters$sigma[[index]])
  if (startsWith(name, "OMEGA")) return(parameters$omega[[index]])
  NA_real_
}

.nm_log_prior <- function(model, parameters) {
  priors <- model$LIK_CONFIG$priors
  if (is.null(priors) || !nrow(priors)) return(0)
  total <- 0
  for (row in seq_len(nrow(priors))) {
    value <- .nm_parameter_value(parameters, priors$parameter[[row]])
    if (!is.finite(value)) return(-Inf)
    density <- switch(
      priors$distribution[[row]],
      normal = stats::dnorm(value, priors$mean[[row]], priors$sd[[row]], log = TRUE),
      lognormal = if (value > 0) stats::dlnorm(
        value, priors$mean[[row]], priors$sd[[row]], log = TRUE
      ) else -Inf,
      half_normal = if (value >= 0) {
        log(2) + stats::dnorm(value, priors$mean[[row]], priors$sd[[row]], log = TRUE)
      } else -Inf,
      inverse_gamma = if (value > 0) {
        priors$shape[[row]] * log(priors$rate[[row]]) -
          lgamma(priors$shape[[row]]) -
          (priors$shape[[row]] + 1) * log(value) - priors$rate[[row]] / value
      } else -Inf
    )
    total <- total + density
  }
  total
}

.nm_prior_nll <- function(model, parameters) -2 * .nm_log_prior(model, parameters)

# Resolve prior targets and normalization constants once for iterative
# stochastic estimators. This is algebraically identical to `.nm_log_prior()`
# but avoids reparsing parameter names and dispatching density functions on
# every MCMC proposal.
.nm_prior_evaluator <- function(model) {
  priors <- model$LIK_CONFIG$priors
  if (is.null(priors) || !nrow(priors)) {
    return(list(
      log_density = function(parameters) 0,
      nll = function(parameters) 0,
      count = 0L
    ))
  }
  names <- toupper(as.character(priors$parameter))
  family <- sub("[0-9]+$", "", names)
  index <- as.integer(sub("^[A-Z]+", "", names))
  distribution <- as.character(priors$distribution)
  mean <- as.numeric(priors$mean)
  sd <- as.numeric(priors$sd)
  shape <- as.numeric(priors$shape)
  rate <- as.numeric(priors$rate)
  normal_constant <- -log(sd) - 0.5 * log(2 * pi)
  inverse_constant <- shape * log(rate) - lgamma(shape)
  values <- function(parameters) {
    result <- numeric(length(index))
    theta <- family == "THETA"
    sigma <- family == "SIGMA"
    omega <- family == "OMEGA"
    if (any(theta)) result[theta] <- parameters$theta[index[theta]]
    if (any(sigma)) result[sigma] <- parameters$sigma[index[sigma]]
    if (any(omega)) result[omega] <- parameters$omega[index[omega]]
    result
  }
  log_density <- function(parameters) {
    value <- values(parameters)
    if (any(!is.finite(value))) return(-Inf)
    density <- numeric(length(value))
    normal <- distribution == "normal"
    if (any(normal)) {
      standardized <- (value[normal] - mean[normal]) / sd[normal]
      density[normal] <- normal_constant[normal] - 0.5 * standardized^2
    }
    lognormal <- distribution == "lognormal"
    if (any(lognormal)) {
      if (any(value[lognormal] <= 0)) return(-Inf)
      logged <- log(value[lognormal])
      standardized <- (logged - mean[lognormal]) / sd[lognormal]
      density[lognormal] <- normal_constant[lognormal] - logged -
        0.5 * standardized^2
    }
    half_normal <- distribution == "half_normal"
    if (any(half_normal)) {
      if (any(value[half_normal] < 0)) return(-Inf)
      standardized <- (value[half_normal] - mean[half_normal]) /
        sd[half_normal]
      density[half_normal] <- log(2) + normal_constant[half_normal] -
        0.5 * standardized^2
    }
    inverse_gamma <- distribution == "inverse_gamma"
    if (any(inverse_gamma)) {
      if (any(value[inverse_gamma] <= 0)) return(-Inf)
      selected <- value[inverse_gamma]
      density[inverse_gamma] <- inverse_constant[inverse_gamma] -
        (shape[inverse_gamma] + 1) * log(selected) -
        rate[inverse_gamma] / selected
    }
    sum(density)
  }
  list(
    log_density = log_density,
    nll = function(parameters) -2 * log_density(parameters),
    count = length(index)
  )
}

.nm_prior_nll_native_gradient <- function(model, parameters) {
  gradient <- numeric(
    nrow(model$THETAS) + nrow(model$SIGMAS) + nrow(model$OMEGAS)
  )
  priors <- model$LIK_CONFIG$priors
  if (is.null(priors) || !nrow(priors)) return(gradient)
  offsets <- c(THETA = 0L, SIGMA = nrow(model$THETAS),
               OMEGA = nrow(model$THETAS) + nrow(model$SIGMAS))
  for (row in seq_len(nrow(priors))) {
    name <- toupper(priors$parameter[[row]])
    family <- sub("[0-9]+$", "", name)
    index <- as.integer(sub("^[A-Z]+", "", name))
    value <- .nm_parameter_value(parameters, name)
    derivative <- switch(
      priors$distribution[[row]],
      normal = 2 * (value - priors$mean[[row]]) / priors$sd[[row]]^2,
      half_normal = 2 * (value - priors$mean[[row]]) / priors$sd[[row]]^2,
      lognormal = 2 / value +
        2 * (log(value) - priors$mean[[row]]) /
          (priors$sd[[row]]^2 * value),
      inverse_gamma = 2 * (priors$shape[[row]] + 1) / value -
        2 * priors$rate[[row]] / value^2
    )
    gradient[offsets[[family]] + index] <-
      gradient[offsets[[family]] + index] + derivative
  }
  gradient
}
