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
