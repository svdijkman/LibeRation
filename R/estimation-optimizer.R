# Outer parameter maps, objectives, optimizers, and conditional modes.
# Split from estimation.R as a behaviour-preserving source move.

.nm_outer_map <- function(model) {
  theta_fixed <- model$THETAS$FIX
  theta_free <- which(!theta_fixed)
  sigma_free <- which(!model$SIGMAS$FIX)
  omega_full <- any(model$OMEGAS$ROW != model$OMEGAS$COL)
  if (omega_full && any(model$OMEGAS$FIX) && !all(model$OMEGAS$FIX)) {
    .nm_stop("A correlated OMEGA must currently be either entirely fixed or entirely estimated.")
  }
  omega_free <- if (omega_full && all(model$OMEGAS$FIX)) integer() else {
    which(!model$OMEGAS$FIX)
  }
  omega_encode <- function(values) {
    if (!length(omega_free)) return(numeric())
    if (!omega_full) return(log(values[omega_free]))
    lower <- t(chol(.nm_omega_matrix(model, values)))
    vapply(seq_len(nrow(model$OMEGAS)), function(i) {
      row <- model$OMEGAS$ROW[[i]]
      column <- model$OMEGAS$COL[[i]]
      if (row == column) log(lower[row, column]) else lower[row, column]
    }, numeric(1))
  }
  omega_decode <- function(encoded) {
    if (!length(omega_free)) return(model$OMEGAS$Value)
    if (!omega_full) {
      values <- model$OMEGAS$Value
      values[omega_free] <- exp(encoded)
      return(values)
    }
    lower <- matrix(0, model$n_eta, model$n_eta)
    for (i in seq_len(nrow(model$OMEGAS))) {
      row <- model$OMEGAS$ROW[[i]]
      column <- model$OMEGAS$COL[[i]]
      lower[row, column] <- if (row == column) exp(encoded[[i]]) else encoded[[i]]
    }
    covariance <- lower %*% t(lower)
    vapply(seq_len(nrow(model$OMEGAS)), function(i) {
      covariance[model$OMEGAS$ROW[[i]], model$OMEGAS$COL[[i]]]
    }, numeric(1))
  }
  start <- c(model$THETAS$Value[theta_free],
             log(model$SIGMAS$Value[sigma_free]),
             omega_encode(model$OMEGAS$Value))
  theta_lower <- model$THETAS$LOWER %||% rep(-Inf, nrow(model$THETAS))
  theta_upper <- model$THETAS$UPPER %||% rep(Inf, nrow(model$THETAS))
  omega_parameter_count <- if (omega_full && length(omega_free)) {
    nrow(model$OMEGAS)
  } else length(omega_free)
  lower <- c(theta_lower[theta_free], rep(-Inf, length(sigma_free) + omega_parameter_count))
  upper <- c(theta_upper[theta_free], rep(Inf, length(sigma_free) + omega_parameter_count))
  parameter_names <- c(
    if (length(theta_free)) paste0("THETA", theta_free) else character(),
    if (length(sigma_free)) paste0("log_SIGMA", sigma_free) else character(),
    if (omega_parameter_count) {
      if (omega_full) paste0("OMEGA_CHOL", seq_len(omega_parameter_count))
      else paste0("log_OMEGA", omega_free)
    } else character()
  )
  decode <- function(parameters) {
    cursor <- 0L
    theta <- model$THETAS$Value
    if (length(theta_free)) {
      theta[theta_free] <- parameters[seq_len(length(theta_free))]
      cursor <- length(theta_free)
    }
    sigma <- model$SIGMAS$Value
    if (length(sigma_free)) {
      sigma_index <- cursor + seq_len(length(sigma_free))
      sigma[sigma_free] <- exp(parameters[sigma_index])
      cursor <- cursor + length(sigma_free)
    }
    omega <- model$OMEGAS$Value
    omega_parameter_count <- if (omega_full && length(omega_free)) {
      nrow(model$OMEGAS)
    } else length(omega_free)
    if (omega_parameter_count) {
      omega_index <- cursor + seq_len(omega_parameter_count)
      omega <- omega_decode(parameters[omega_index])
    }
    list(theta = theta, sigma = sigma, omega = omega)
  }
  encode <- function(parameters) c(
    parameters$theta[theta_free],
    log(parameters$sigma[sigma_free]),
    omega_encode(parameters$omega)
  )
  log_jacobian <- function(parameters) {
    value <- sum(log(parameters$sigma[sigma_free]))
    if (!length(omega_free)) return(value)
    if (!omega_full) return(value + sum(log(parameters$omega[omega_free])))
    lower <- t(chol(.nm_omega_matrix(model, parameters$omega)))
    value + model$n_eta * log(2) + sum(
      (model$n_eta + 2L - seq_len(model$n_eta)) * log(diag(lower))
    )
  }
  log_jacobian_gradient <- function(parameters) {
    result <- numeric(length(start))
    cursor <- length(theta_free)
    if (length(sigma_free)) {
      result[cursor + seq_len(length(sigma_free))] <- 1
      cursor <- cursor + length(sigma_free)
    }
    if (!length(omega_free)) return(result)
    if (!omega_full) {
      result[cursor + seq_len(length(omega_free))] <- 1
      return(result)
    }
    for (encoded in seq_len(nrow(model$OMEGAS))) {
      row <- model$OMEGAS$ROW[[encoded]]
      column <- model$OMEGAS$COL[[encoded]]
      if (row == column) {
        result[cursor + encoded] <- model$n_eta + 2L - row
      }
    }
    result
  }
  in_bounds <- function(parameters) {
    length(parameters) == length(lower) &&
      all(parameters >= lower) && all(parameters <= upper)
  }
  jacobian <- function(parameters) {
    n_native <- nrow(model$THETAS) + nrow(model$SIGMAS) + nrow(model$OMEGAS)
    result <- matrix(0, n_native, length(start))
    cursor <- 0L
    for (index in theta_free) {
      cursor <- cursor + 1L
      result[index, cursor] <- 1
    }
    sigma_offset <- nrow(model$THETAS)
    for (index in sigma_free) {
      cursor <- cursor + 1L
      result[sigma_offset + index, cursor] <- parameters$sigma[[index]]
    }
    omega_offset <- sigma_offset + nrow(model$SIGMAS)
    if (!length(omega_free)) return(result)
    if (!omega_full) {
      for (index in omega_free) {
        cursor <- cursor + 1L
        result[omega_offset + index, cursor] <- parameters$omega[[index]]
      }
      return(result)
    }
    covariance <- .nm_omega_matrix(model, parameters$omega)
    lower_cholesky <- t(chol(covariance))
    for (encoded in seq_len(nrow(model$OMEGAS))) {
      cursor <- cursor + 1L
      row <- model$OMEGAS$ROW[[encoded]]
      column <- model$OMEGAS$COL[[encoded]]
      derivative_lower <- matrix(0, model$n_eta, model$n_eta)
      derivative_lower[row, column] <- if (row == column) {
        lower_cholesky[row, column]
      } else 1
      derivative <- derivative_lower %*% t(lower_cholesky) +
        lower_cholesky %*% t(derivative_lower)
      for (native in seq_len(nrow(model$OMEGAS))) {
        result[omega_offset + native, cursor] <- derivative[
          model$OMEGAS$ROW[[native]], model$OMEGAS$COL[[native]]
        ]
      }
    }
    result
  }
  list(start = start, lower = lower, upper = upper, names = parameter_names,
       decode = decode, encode = encode, in_bounds = in_bounds, theta_free = theta_free,
       sigma_free = sigma_free, omega_free = omega_free,
       omega_full = omega_full, log_jacobian = log_jacobian,
       omega_parameterization = if (omega_full) "log_cholesky" else "log_variance",
       log_jacobian_gradient = log_jacobian_gradient, jacobian = jacobian)
}

.nm_cpp_prior_config <- function(model) {
  priors <- model$LIK_CONFIG$priors
  if (is.null(priors) || !nrow(priors)) {
    return(list(
      index = integer(), family = character(), mean = numeric(), sd = numeric(),
      shape = numeric(), rate = numeric()
    ))
  }
  offsets <- c(
    THETA = 0L, SIGMA = nrow(model$THETAS),
    OMEGA = nrow(model$THETAS) + nrow(model$SIGMAS)
  )
  family <- sub("[0-9]+$", "", toupper(priors$parameter))
  index <- as.integer(sub("^[A-Z]+", "", priors$parameter))
  list(
    index = unname(offsets[family]) + index,
    family = as.character(priors$distribution),
    mean = as.numeric(priors$mean), sd = as.numeric(priors$sd),
    shape = as.numeric(priors$shape), rate = as.numeric(priors$rate)
  )
}

.nm_cpp_population_objective <- function(context, map, approximation,
                                          eta_maxit, tolerance,
                                          initial_eta = NULL) {
  started <- proc.time()[["elapsed"]]
  finish <- function(result) {
    result$initialization_seconds <- unname(
      proc.time()[["elapsed"]] - started
    )
    result
  }
  if (!isTRUE(getOption("LibeRation.cpp_population_objective", TRUE))) {
    return(finish(list(pointer = NULL, reason = "disabled by option")))
  }
  if (!is.null(context$parallel)) {
    return(finish(list(
      pointer = NULL, reason = "PSOCK workers require R coordination"
    )))
  }
  approximation <- match.arg(
    tolower(approximation), c("fo", "its", "foce", "focei", "laplace")
  )
  parameters <- map$decode(map$start)
  mu <- if (.nm_liber_optimized(context) && approximation != "fo" &&
            context$n_eta > 0L) {
    .nm_mu_specialization(context, map, enabled = TRUE)
  } else NULL
  primary <- curvature <- list()
  # Adaptive ODE tapes are owned and retaped by the C++ population object.
  # Analytical models retain the already-recorded subject tapes, including
  # structurally shared prediction tapes used to construct curvature tapes.
  if (!isTRUE(context$model$USE_ODE)) {
    if (approximation == "fo") {
      fo_pool <- new.env(parent = emptyenv())
      invisible(lapply(context$subjects, function(evaluator) {
        .nm_fo_pool_tape(
          fo_pool, evaluator, parameters$theta, parameters$sigma, parameters$omega
        )
      }))
      primary <- lapply(context$subjects, function(evaluator) evaluator$fo_tape$pointer)
    } else {
      interaction <- approximation != "foce"
      primary <- lapply(context$subjects, function(evaluator) {
        if (interaction) evaluator$objective_tape$pointer else
          evaluator$noninteraction_tape$pointer
      })
      if (approximation %in% c("foce", "focei", "laplace")) {
        curvature_anchors <- if (approximation == "laplace" && context$n_eta) {
          initial_modes <- .nm_subject_modes(
            context, parameters, starts = initial_eta,
            maxit = eta_maxit, tolerance = tolerance,
            interaction = TRUE, exact_hessian = TRUE
          )
          if (any(vapply(initial_modes, `[[`, integer(1), "convergence") != 0L)) {
            .nm_stop("Initial conditional modes did not converge for the compiled Laplace objective.")
          }
          lapply(initial_modes, `[[`, "par")
        } else {
          if (is.null(initial_eta)) {
            rep(list(rep(0, context$n_eta)), context$n_subjects)
          } else {
            lapply(seq_len(context$n_subjects), function(subject) initial_eta[subject, ])
          }
        }
        invisible(Map(function(evaluator, eta) {
          evaluator$ensure_curvature_tape(
            parameters$theta, eta, parameters$sigma,
            parameters$omega, approximation
          )
        }, context$subjects, curvature_anchors))
        curvature <- lapply(context$subjects, function(evaluator) {
          evaluator$curvature_tapes[[approximation]]$pointer
        })
      }
    }
  }
  priors <- .nm_cpp_prior_config(context$model)
  config <- list(
    approximation = approximation,
    theta = parameters$theta, sigma = parameters$sigma, omega = parameters$omega,
    theta_free = map$theta_free, sigma_free = map$sigma_free,
    omega_free = map$omega_free, omega_full = map$omega_full,
    omega_rows = context$model$OMEGAS$ROW,
    omega_cols = context$model$OMEGAS$COL,
    n_eta = context$n_eta, n_eta_base = context$model$n_eta,
    eta_maxit = as.integer(eta_maxit), tolerance = as.numeric(tolerance),
    use_ode = isTRUE(context$model$USE_ODE),
    guard_radius = as.numeric(getOption("LibeRation.tape_guard_radius", 0.5)),
    start = map$start,
    eta_start = initial_eta %||% matrix(0, context$n_subjects, context$n_eta),
    prior_index = priors$index, prior_family = priors$family,
    prior_mean = priors$mean, prior_sd = priors$sd,
    prior_shape = priors$shape, prior_rate = priors$rate,
    mu = if (isTRUE(mu$active) && length(mu$theta)) list(
      active = TRUE, theta = as.integer(mu$theta),
      links = unname(mu$links[as.character(mu$theta)]),
      design_columns = unname(mu$design_columns)
    ) else list(active = FALSE),
    fo_population_batch = isTRUE(getOption("LibeRation.fo_population_batch", TRUE)),
    fo_population_scalar = isTRUE(getOption(
      "LibeRation.fo_population_scalar", TRUE
    )),
    fo_low_rank = isTRUE(getOption("LibeRation.fo_low_rank", TRUE)),
    fo_low_rank_tolerance = as.numeric(getOption(
      "LibeRation.fo_low_rank_tolerance", 1e-9
    )),
    fo_low_rank_condition_tolerance = as.numeric(getOption(
      "LibeRation.fo_low_rank_condition_tolerance", 1e-12
    )),
    fo_population_max_operations = as.numeric(getOption(
      "LibeRation.fo_population_max_operations", 2e6
    )),
    subject_materializer = NULL,
    subject_store = if (is.null(context$subject_store)) NULL else
      context$subject_store$pointer
  )
  subject_inputs <- lapply(context$subjects, function(evaluator) {
    if (is.null(context$subject_store)) return(evaluator$data_frame())
    if (approximation == "fo") return(evaluator$fo_dynamic %||% numeric())
    if (approximation == "foce") evaluator$noninteraction_dynamic %||% numeric() else
      evaluator$objective_dynamic %||% numeric()
  })
  finish(tryCatch(
    list(
      pointer = .liberation_population_objective_create(
        context$engine$pointer,
        subject_inputs,
        primary, curvature, config
      ),
      reason = NULL
    ),
    error = function(error) list(
      pointer = NULL,
      reason = paste("compiled population initialization failed:", conditionMessage(error))
    )
  ))
}

.nm_log_gradient <- function(iteration, objective, parameters, map, value = NULL,
                             gradient_function = NULL) {
  baseline <- value %||% objective(parameters)
  gradient <- if (is.function(gradient_function)) {
    as.numeric(gradient_function(parameters))
  } else {
    result <- numeric(length(parameters))
    for (index in seq_along(parameters)) {
      step <- 1e-5 * max(abs(parameters[[index]]), 1)
      low <- high <- parameters
      low[[index]] <- max(map$lower[[index]], parameters[[index]] - step)
      high[[index]] <- min(map$upper[[index]], parameters[[index]] + step)
      low_value <- if (low[[index]] < parameters[[index]]) objective(low) else baseline
      high_value <- if (high[[index]] > parameters[[index]]) objective(high) else baseline
      width <- high[[index]] - low[[index]]
      result[[index]] <- if (width > 0 && is.finite(low_value) && is.finite(high_value)) {
        (high_value - low_value) / width
      } else NA_real_
    }
    result
  }
  names(gradient) <- map$names
  cat(sprintf(
    "[LibeRation] OUTER EVALUATION %d OFV %.10g SCALED GRADIENT %s\n",
    as.integer(iteration), as.numeric(baseline),
    paste(sprintf("%s=%.6g", names(gradient), gradient), collapse = " ")
  ))
  try(flush(stdout()), silent = TRUE)
  invisible(gradient)
}

.nm_outer_optim <- function(map, objective, maxit, tolerance, trace = 0L,
                            print_every = 0L, gradient = NULL,
                            optimizer_backend = c("auto", "native", "r"),
                            compiled_objective = NULL,
                            strict_convergence = FALSE,
                            allow_fd_gradient = getOption(
                              "LibeRation.allow_fd_gradient", FALSE
                            )) {
  optimizer_backend <- match.arg(optimizer_backend)
  if (optimizer_backend == "auto") optimizer_backend <- "r"
  compiled_pointer <- compiled_objective$pointer %||% NULL
  compiled <- !is.null(compiled_pointer)
  objective_initialization_seconds <- as.numeric(
    compiled_objective$initialization_seconds %||% 0
  )
  started <- proc.time()[["elapsed"]]
  print_every <- as.integer(print_every)
  evaluations <- 0L
  gradient_evaluations <- 0L
  gradient_fallbacks <- 0L
  gradient_fallback_evaluations <- 0L
  fd_warning_emitted <- FALSE
  pending_log <- NULL
  objective_scale <- 1
  raw <- function(parameters) {
    if (!map$in_bounds(parameters)) return(1e100)
    value <- tryCatch(
      if (compiled) {
        .liberation_population_objective_value(compiled_pointer, parameters)
      } else objective(map$decode(parameters)),
      error = function(e) Inf
    )
    if (length(value) != 1L || !is.finite(value)) 1e100 else value
  }
  safe <- function(parameters) {
    evaluations <<- evaluations + 1L
    value <- raw(parameters)
    if (print_every > 0L &&
        (evaluations == 1L || evaluations %% print_every == 0L)) {
      if (is.function(gradient) || compiled) {
        pending_log <<- list(
          iteration = evaluations, parameters = parameters, value = value
        )
      } else {
        .nm_log_gradient(evaluations, raw, parameters, map, value)
      }
    }
    value
  }
  finite_difference_gradient <- function(parameters, baseline = NULL) {
    baseline <- baseline %||% raw(parameters)
    valid <- function(value) {
      length(value) == 1L && is.finite(value) && value < 1e99
    }
    result <- rep(NA_real_, length(parameters))
    for (index in seq_along(parameters)) {
      step <- 1e-5 * max(abs(parameters[[index]]), 1)
      low <- high <- parameters
      low[[index]] <- max(map$lower[[index]], parameters[[index]] - step)
      high[[index]] <- min(map$upper[[index]], parameters[[index]] + step)
      low_value <- if (low[[index]] < parameters[[index]]) {
        gradient_fallback_evaluations <<- gradient_fallback_evaluations + 1L
        raw(low)
      } else baseline
      high_value <- if (high[[index]] > parameters[[index]]) {
        gradient_fallback_evaluations <<- gradient_fallback_evaluations + 1L
        raw(high)
      } else baseline
      result[[index]] <- if (
        low[[index]] < parameters[[index]] &&
          high[[index]] > parameters[[index]] &&
          valid(low_value) && valid(high_value)
      ) {
        (high_value - low_value) / (high[[index]] - low[[index]])
      } else if (
        high[[index]] > parameters[[index]] &&
          valid(baseline) && valid(high_value)
      ) {
        (high_value - baseline) / (high[[index]] - parameters[[index]])
      } else if (
        low[[index]] < parameters[[index]] &&
          valid(baseline) && valid(low_value)
      ) {
        (baseline - low_value) / (parameters[[index]] - low[[index]])
      } else NA_real_
    }
    result
  }
  safe_gradient <- if (is.function(gradient) || compiled) function(parameters) {
    gradient_evaluations <<- gradient_evaluations + 1L
    value <- tryCatch(
      if (compiled) {
        as.numeric(.liberation_population_objective_gradient(
          compiled_pointer, parameters
        ))
      } else as.numeric(gradient(map$decode(parameters))),
      error = function(error) rep(NA_real_, length(parameters))
    )
    if (length(value) != length(parameters) || any(!is.finite(value))) {
      # L-BFGS-B evaluates the gradient at its generalized Cauchy point before
      # line-searching back. A large likelihood gradient can put that point on
      # a numerically invalid boundary even when the starting point is sound.
      # Return a finite inward barrier direction there; retain the hard error
      # for a non-finite derivative at a finite objective value.
      point_value <- raw(parameters)
      if (is.finite(point_value) && point_value >= 1e99) {
        parameter_scale <- pmax(abs(map$start), 1)
        inward <- (parameters - map$start) / parameter_scale
        largest <- max(abs(inward))
        if (is.finite(largest) && largest > 0) {
          return(objective_scale * inward / largest)
        }
      }
      if (!isTRUE(allow_fd_gradient)) {
        .nm_stop(
          "The population objective gradient is not finite. Automatic finite-",
          "difference substitution is disabled because it changes the estimator's ",
          "derivative contract. Resolve the tape/path failure or rerun with ",
          "`allow_fd_gradient = TRUE` and inspect `gradient_fallbacks`."
        )
      }
      fallback <- finite_difference_gradient(parameters, point_value)
      if (length(fallback) == length(parameters) &&
          all(is.finite(fallback))) {
        gradient_fallbacks <<- gradient_fallbacks + 1L
        if (!fd_warning_emitted) {
          warning(
            "Using an explicitly enabled finite-difference population-gradient fallback; ",
            "the completed fit must be reviewed via `diagnostics$gradient_fallbacks`.",
            call. = FALSE
          )
          fd_warning_emitted <<- TRUE
        }
        value <- fallback
      } else {
        .nm_stop("The population objective gradient is not finite.")
      }
    }
    if (!is.null(pending_log) &&
        identical(as.numeric(parameters), as.numeric(pending_log$parameters))) {
      .nm_log_gradient(
        pending_log$iteration, raw, parameters, map, pending_log$value,
        gradient_function = function(ignored) value
      )
      pending_log <<- NULL
    }
    value
  } else NULL
  if (!length(map$start)) {
    result <- list(
      par = numeric(), value = safe(numeric()), convergence = 0L,
      counts = c(`function` = 1L, gradient = NA_integer_),
      iterations = 0L, objective_evaluations = 1L,
      gradient_evaluations = 0L, backend = "fixed-parameters",
      gradient_fallbacks = 0L, gradient_fallback_evaluations = 0L,
      elapsed_seconds = unname(proc.time()[["elapsed"]] - started),
      objective_initialization_seconds = objective_initialization_seconds,
      message = NULL,
      objective_backend = if (compiled) "persistent-cpp-population-objective" else
        "r-orchestrated-population-objective",
      population_objective = if (compiled) {
        .liberation_population_objective_telemetry(compiled_pointer)
      } else NULL
    )
    if (compiled) {
      result$objective_backend <- result$population_objective$backend %||%
        result$objective_backend
    }
    return(result)
  }
  if (optimizer_backend == "native" && compiled) {
    result <- .liberation_population_objective_native_optimizer(
      compiled_pointer, map$start, map$lower, map$upper,
      as.integer(maxit), as.numeric(tolerance), as.integer(trace)
    )
    result$backend <- "native-bfgs"
    result$coordinator <- "direct-cpp-population-objective"
    result$objective_backend <- "persistent-cpp-population-objective"
    result$population_objective <-
      .liberation_population_objective_telemetry(compiled_pointer)
    result$objective_backend <- result$population_objective$backend %||%
      result$objective_backend
    result$elapsed_seconds <- unname(proc.time()[["elapsed"]] - started)
    result$objective_initialization_seconds <- objective_initialization_seconds
    result$gradient_fallbacks <- 0L
    result$gradient_fallback_evaluations <- 0L
    return(result)
  }
  if (optimizer_backend == "native" && is.function(safe_gradient)) {
    result <- .liberation_native_optimizer(
      safe, safe_gradient, map$start, map$lower, map$upper,
      as.integer(maxit), as.numeric(tolerance), as.integer(trace)
    )
    result$backend <- "native-bfgs-r-callback"
    result$objective_backend <- if (compiled) {
      "persistent-cpp-population-objective"
    } else "r-orchestrated-population-objective"
    result$population_objective <- if (compiled) {
      .liberation_population_objective_telemetry(compiled_pointer)
    } else NULL
    if (compiled) {
      result$objective_backend <- result$population_objective$backend %||%
        result$objective_backend
    }
    result$elapsed_seconds <- unname(proc.time()[["elapsed"]] - started)
    result$objective_initialization_seconds <- objective_initialization_seconds
    result$gradient_fallbacks <- as.integer(gradient_fallbacks)
    result$gradient_fallback_evaluations <- as.integer(
      gradient_fallback_evaluations
    )
    return(result)
  }
  bounded <- any(is.finite(map$lower)) || any(is.finite(map$upper))
  initial_value <- raw(map$start)
  if (is.finite(initial_value) && initial_value < 1e99) {
    objective_scale <- max(abs(initial_value), 1)
  }
  arguments <- list(
    par = map$start, fn = safe,
    method = if (bounded) "L-BFGS-B" else if (is.function(safe_gradient) ||
      length(map$start) == 1L) "BFGS" else "Nelder-Mead",
    control = list(
      maxit = as.integer(maxit), reltol = tolerance, trace = trace,
      fnscale = objective_scale
    )
  )
  if (is.function(safe_gradient)) arguments$gr <- safe_gradient
  if (bounded) {
    arguments$lower <- map$lower
    arguments$upper <- map$upper
    if (isTRUE(strict_convergence)) {
      # Exact FO gradients can be very large at the initial point. A loose
      # function-reduction test stopped before OMEGA reached the same solution
      # under mathematically equivalent summation orders. Keep the function
      # test at machine accuracy and use a squared scaled-gradient target.
      arguments$control$factr <- 1
      arguments$control$pgtol <- tolerance^2
    } else {
      arguments$control$factr <- max(tolerance / .Machine$double.eps, 1)
    }
    arguments$control$reltol <- NULL
  }
  result <- do.call(stats::optim, arguments)
  iterations <- suppressWarnings(as.integer(result$counts[["gradient"]]))
  if (!length(iterations) || is.na(iterations)) {
    iterations <- suppressWarnings(as.integer(result$counts[["function"]]))
  }
  result$iterations <- if (!length(iterations) || is.na(iterations)) {
    as.integer(evaluations)
  } else iterations
  result$objective_evaluations <- as.integer(evaluations)
  result$gradient_evaluations <- as.integer(gradient_evaluations)
  result$gradient_fallbacks <- as.integer(gradient_fallbacks)
  result$gradient_fallback_evaluations <- as.integer(
    gradient_fallback_evaluations
  )
  result$backend <- if (compiled) {
    paste0("r-", tolower(arguments$method), "-cpp-objective")
  } else if (is.function(safe_gradient)) "r-optim-gradient" else "r-optim"
  result$population_objective <- if (compiled) {
    .liberation_population_objective_telemetry(compiled_pointer)
  } else NULL
  result$objective_backend <- if (compiled) {
    result$population_objective$backend %||% "persistent-cpp-population-objective"
  } else "r-orchestrated-population-objective"
  result$objective_scale <- objective_scale
  result$elapsed_seconds <- unname(proc.time()[["elapsed"]] - started)
  result$objective_initialization_seconds <- objective_initialization_seconds
  result
}

.nm_subject_modes <- function(context, parameters, starts = NULL,
                              maxit = 100L, tolerance = 1e-7,
                              interaction = TRUE, exact_hessian = TRUE) {
  if (is.null(starts)) starts <- matrix(0, context$n_subjects, context$n_eta)
  batch_modes <- function(evaluators, starts) {
    if (!length(evaluators)) return(list())
    native_eligible <- all(vapply(evaluators, function(evaluator) {
      tape <- if (isTRUE(interaction)) evaluator$objective_tape else
        evaluator$noninteraction_tape
      !is.null(tape) && !is.null(tape$pointer)
    }, logical(1)))
    if (!native_eligible) {
      return(lapply(seq_along(evaluators), function(subject) {
        evaluators[[subject]]$eta_mode(
          parameters$theta, parameters$sigma, parameters$omega,
          start = starts[subject, ], maxit = maxit, tolerance = tolerance
        )
      }))
    }
    ode_guard <- isTRUE(evaluators[[1L]]$engine$model$USE_ODE)
    if (ode_guard) invisible(Map(function(evaluator, subject) {
      evaluator$ensure_valid_tapes(
        parameters$theta, parameters$sigma, parameters$omega, starts[subject, ]
      )
    }, evaluators, seq_along(evaluators)))
    points <- cbind(
      matrix(parameters$theta, nrow(starts), length(parameters$theta), byrow = TRUE),
      starts,
      matrix(parameters$sigma, nrow(starts), length(parameters$sigma), byrow = TRUE),
      matrix(parameters$omega, nrow(starts), length(parameters$omega), byrow = TRUE)
    )
    tapes <- lapply(evaluators, function(evaluator) {
      if (isTRUE(interaction)) evaluator$objective_tape$pointer else
        evaluator$noninteraction_tape$pointer
    })
    raw <- tryCatch(
      .liberation_objective_tape_eta_modes(
        tapes, points, length(parameters$theta) + seq_len(context$n_eta), starts,
        as.integer(maxit), as.numeric(tolerance), isTRUE(exact_hessian),
        lapply(
          evaluators, .nm_subject_dynamic_input, interaction = interaction
        ),
        isTRUE(.nm_liber_optimized(context))
      ), error = identity
    )
    if (inherits(raw, "error")) {
      if (!grepl("CppAD tape path changed", conditionMessage(raw), fixed = TRUE)) {
        stop(raw)
      }
      return(lapply(seq_along(evaluators), function(subject) {
        evaluators[[subject]]$eta_mode(
          parameters$theta, parameters$sigma, parameters$omega,
          start = starts[subject, ], maxit = maxit, tolerance = tolerance,
          interaction = interaction, exact_hessian = exact_hessian
        )
      }))
    }
    lapply(seq_along(raw), function(subject) {
      mode <- raw[[subject]]
      if (!identical(as.integer(mode$convergence), 0L)) {
        return(evaluators[[subject]]$eta_mode(
          parameters$theta, parameters$sigma, parameters$omega,
          start = starts[subject, ], maxit = maxit, tolerance = tolerance,
          interaction = interaction, exact_hessian = exact_hessian
        ))
      }
      if (ode_guard && evaluators[[subject]]$ensure_valid_tapes(
        parameters$theta, parameters$sigma, parameters$omega, mode$par
      )) {
        return(evaluators[[subject]]$eta_mode(
          parameters$theta, parameters$sigma, parameters$omega,
          start = mode$par, maxit = maxit, tolerance = tolerance,
          interaction = interaction, exact_hessian = exact_hessian
        ))
      }
      curvature <- if (isTRUE(exact_hessian)) {
        .nm_positive_definite(mode$hessian, "Conditional ETA curvature")
      } else list(matrix = matrix(numeric(), 0L, 0L), logdet = 0, jitter = 0)
      list(
        par = as.numeric(mode$par), value = as.numeric(mode$value),
        convergence = 0L, hessian = curvature$matrix,
        logdet = curvature$logdet, jitter = curvature$jitter,
        gradient = as.numeric(mode$gradient),
        iterations = as.integer(mode$iterations),
        evaluations = as.integer(mode$evaluations),
        optimizer_state_reused = isTRUE(mode$optimizer_state_reused),
        backend = "cpp-batch"
      )
    })
  }
  if (!is.null(context$parallel)) {
    chunks <- context$parallel$chunks
    pieces <- parallel::clusterApply(
      context$parallel$cluster, seq_along(chunks),
       function(index, start_chunks, theta, sigma, omega, maxit, tolerance,
                interaction, exact_hessian, optimized) {
        worker_state <- get(
          ".nm_parallel_worker_state", envir = asNamespace("LibeRation")
        )
        evaluators <- worker_state()$subjects
        worker_starts <- start_chunks[[index]]
        context <- list(n_eta = ncol(worker_starts), optimized = optimized)
        parameters <- list(theta = theta, sigma = sigma, omega = omega)
        batch <- get(".nm_subject_modes_batch", envir = asNamespace("LibeRation"))
        batch(evaluators, context, parameters, worker_starts, maxit, tolerance,
              interaction, exact_hessian)
      },
      start_chunks = lapply(chunks, function(rows) starts[rows, , drop = FALSE]),
      theta = parameters$theta, sigma = parameters$sigma,
       omega = parameters$omega, maxit = maxit, tolerance = tolerance,
       interaction = interaction, exact_hessian = exact_hessian,
       optimized = isTRUE(.nm_liber_optimized(context))
    )
    return(unlist(pieces, recursive = FALSE))
  }
  batch_modes(context$subjects, starts)
}

.nm_subject_modes_batch <- function(evaluators, context, parameters, starts,
                                    maxit, tolerance, interaction,
                                    exact_hessian) {
  if (!length(evaluators)) return(list())
  ode_guard <- isTRUE(evaluators[[1L]]$engine$model$USE_ODE)
  if (ode_guard) invisible(Map(function(evaluator, subject) {
    evaluator$ensure_valid_tapes(
      parameters$theta, parameters$sigma, parameters$omega, starts[subject, ]
    )
  }, evaluators, seq_along(evaluators)))
  points <- cbind(
    matrix(parameters$theta, nrow(starts), length(parameters$theta), byrow = TRUE),
    starts,
    matrix(parameters$sigma, nrow(starts), length(parameters$sigma), byrow = TRUE),
    matrix(parameters$omega, nrow(starts), length(parameters$omega), byrow = TRUE)
  )
  tapes <- lapply(evaluators, function(evaluator) {
    if (isTRUE(interaction)) evaluator$objective_tape$pointer else
      evaluator$noninteraction_tape$pointer
  })
  raw <- tryCatch(
    .liberation_objective_tape_eta_modes(
      tapes, points, length(parameters$theta) + seq_len(context$n_eta), starts,
      as.integer(maxit), as.numeric(tolerance), isTRUE(exact_hessian),
      lapply(
        evaluators, .nm_subject_dynamic_input, interaction = interaction
      ),
      isTRUE(context$optimized)
    ), error = identity
  )
  if (inherits(raw, "error")) {
    if (!grepl("CppAD tape path changed", conditionMessage(raw), fixed = TRUE)) {
      stop(raw)
    }
    return(lapply(seq_along(evaluators), function(subject) {
      evaluators[[subject]]$eta_mode(
        parameters$theta, parameters$sigma, parameters$omega,
        start = starts[subject, ], maxit = maxit, tolerance = tolerance,
        interaction = interaction, exact_hessian = exact_hessian
      )
    }))
  }
  lapply(seq_along(raw), function(subject) {
    mode <- raw[[subject]]
    if (!identical(as.integer(mode$convergence), 0L) ||
        (ode_guard && evaluators[[subject]]$ensure_valid_tapes(
          parameters$theta, parameters$sigma, parameters$omega, mode$par
        ))) {
      return(evaluators[[subject]]$eta_mode(
        parameters$theta, parameters$sigma, parameters$omega,
        start = mode$par, maxit = maxit, tolerance = tolerance,
        interaction = interaction, exact_hessian = exact_hessian
      ))
    }
    curvature <- if (isTRUE(exact_hessian)) {
      .nm_positive_definite(mode$hessian, "Conditional ETA curvature")
    } else list(matrix = matrix(numeric(), 0L, 0L), logdet = 0, jitter = 0)
    list(
      par = as.numeric(mode$par), value = as.numeric(mode$value), convergence = 0L,
      hessian = curvature$matrix, logdet = curvature$logdet,
      jitter = curvature$jitter, gradient = as.numeric(mode$gradient),
      iterations = as.integer(mode$iterations), evaluations = as.integer(mode$evaluations),
      optimizer_state_reused = isTRUE(mode$optimizer_state_reused),
      backend = "cpp-batch"
    )
  })
}

.nm_subject_curvature_logdet_reference <- function(
    context, evaluator, parameters, eta, approximation) {
  eta_columns <- length(parameters$theta) + seq_len(context$n_eta)
  if (approximation == "laplace") {
    curvature <- evaluator$objective_hessian_subset(
      parameters$theta, eta, parameters$sigma, parameters$omega,
      rows = eta_columns, columns = eta_columns, interaction = TRUE
    )
    return(.nm_positive_definite(
      curvature, "Laplace conditional curvature"
    )$logdet)
  }
  prediction <- evaluator$prediction(
    parameters$theta, eta, parameters$sigma, jacobian = TRUE,
    columns = eta_columns
  )
  observations <- evaluator$observation_data()
  observed_rows <- as.integer(attr(observations, "rows")) + 1L
  jacobian <- prediction$jacobian[observed_rows, , drop = FALSE]
  f <- prediction$value[observed_rows]
  scale_f <- if (approximation == "foce") {
    evaluator$prediction(
      parameters$theta, rep(0, context$n_eta), parameters$sigma,
      jacobian = FALSE
    )$value[observed_rows]
  } else f
  dvid <- observations$DVID %||% rep(1L, length(observed_rows))
  variance <- .nm_residual_variance(
    context$model, scale_f, parameters$sigma, dvid
  )
  omega_inverse <- solve(.nm_effect_covariance_evaluator_reference(
    context$model, evaluator, parameters$omega
  ))
  curvature <- 2 * crossprod(jacobian / sqrt(variance)) + 2 * omega_inverse
  .nm_positive_definite(
    curvature, paste0(toupper(approximation), " Gauss-Newton curvature")
  )$logdet
}

.nm_subject_curvature_logdet <- function(context, evaluator, parameters, eta,
                                         approximation) {
  result <- evaluator$curvature(
    parameters$theta, eta, parameters$sigma, parameters$omega,
    approximation, gradient = FALSE
  )
  as.numeric(result$value)
}

.nm_conditional_native_gradient <- function(context, parameters, eta,
                                            interaction = TRUE) {
  if (is.null(context$parallel)) {
    gradients <- .nm_objective_collection_gradient(
      context$subjects, parameters, eta, interaction = interaction
    )
    total <- colSums(gradients)
  } else {
    eta_chunks <- lapply(
      context$parallel$chunks, function(rows) eta[rows, , drop = FALSE]
    )
    pieces <- parallel::clusterApply(
      context$parallel$cluster, seq_along(context$parallel$chunks),
      function(index, eta_chunks, parameters, interaction) {
        worker_state <- get(
          ".nm_parallel_worker_state", envir = asNamespace("LibeRation")
        )
        evaluators <- worker_state()$subjects
        collection <- get(
          ".nm_objective_collection_gradient", envir = asNamespace("LibeRation")
        )
        colSums(collection(
          evaluators, parameters, eta_chunks[[index]], interaction = interaction
        ))
      }, eta_chunks = eta_chunks, parameters = parameters,
      interaction = interaction
    )
    total <- Reduce(`+`, pieces)
  }
  n_theta <- length(parameters$theta)
  n_sigma <- length(parameters$sigma)
  n_omega <- length(parameters$omega)
  population_positions <- c(
    seq_len(n_theta),
    n_theta + context$n_eta + seq_len(n_sigma),
    n_theta + context$n_eta + n_sigma + seq_len(n_omega)
  )
  as.numeric(total[population_positions]) +
    .nm_prior_nll_native_gradient(context$model, parameters)
}

.nm_nested_outer_gradient <- function(context, map, objective, parameters,
                                      approximation, relative_step = 1e-5) {
  value <- objective(parameters)
  if (!is.finite(value)) .nm_stop("Cannot differentiate a non-finite objective.")
  state <- attr(objective, "state")
  modes <- state$modes
  if (is.null(modes)) .nm_stop("Conditional modes are unavailable for differentiation.")
  eta <- if (context$n_eta) {
    do.call(rbind, lapply(modes, `[[`, "par"))
  } else matrix(numeric(), context$n_subjects, 0L)
  interaction <- approximation != "foce"
  transform <- map$jacobian(parameters)
  if (approximation == "its" || !context$n_eta || !ncol(transform)) {
    native <- .nm_conditional_native_gradient(
      context, parameters, eta, interaction = interaction
    )
    return(as.vector(native %*% transform))
  }
  if (is.null(context$parallel)) {
    result <- .nm_nested_gradient_batch(
      context$subjects, context$n_eta, parameters, eta, approximation, transform
    )
    gradient <- result$gradient
  } else {
    eta_chunks <- lapply(
      context$parallel$chunks, function(rows) eta[rows, , drop = FALSE]
    )
    pieces <- parallel::clusterApply(
      context$parallel$cluster, seq_along(context$parallel$chunks),
      function(index, eta_chunks, n_eta, parameters, approximation, transform) {
        worker_state <- get(
          ".nm_parallel_worker_state", envir = asNamespace("LibeRation")
        )
        evaluators <- worker_state()$subjects
        batch <- get(".nm_nested_gradient_batch", envir = asNamespace("LibeRation"))
        batch(
          evaluators, n_eta, parameters, eta_chunks[[index]],
          approximation, transform
        )$gradient
      }, eta_chunks = eta_chunks, n_eta = context$n_eta,
      parameters = parameters, approximation = approximation,
      transform = transform
    )
    gradient <- Reduce(`+`, pieces)
  }
  prior <- .nm_prior_nll_native_gradient(context$model, parameters)
  as.numeric(gradient) + as.vector(prior %*% transform)
}

.nm_nested_gradient_batch <- function(evaluators, n_eta, parameters, eta,
                                      approximation, transform) {
  interaction <- approximation != "foce"
  ode_guard <- isTRUE(evaluators[[1L]]$engine$model$USE_ODE)
  invisible(Map(function(evaluator, subject) {
    if (ode_guard) evaluator$ensure_valid_tapes(
      parameters$theta, parameters$sigma, parameters$omega, eta[subject, ]
    )
    evaluator$ensure_curvature_tape(
      parameters$theta, eta[subject, ], parameters$sigma,
      parameters$omega, approximation
    )
  }, evaluators, seq_along(evaluators)))
  points <- cbind(
    matrix(parameters$theta, nrow(eta), length(parameters$theta), byrow = TRUE),
    eta,
    matrix(parameters$sigma, nrow(eta), length(parameters$sigma), byrow = TRUE),
    matrix(parameters$omega, nrow(eta), length(parameters$omega), byrow = TRUE)
  )
  n_theta <- length(parameters$theta)
  n_sigma <- length(parameters$sigma)
  n_omega <- length(parameters$omega)
  population_positions <- c(
    seq_len(n_theta), n_theta + n_eta + seq_len(n_sigma),
    n_theta + n_eta + n_sigma + seq_len(n_omega)
  )
  .liberation_nested_population_gradient(
    lapply(evaluators, function(evaluator) {
      if (interaction) evaluator$objective_tape$pointer else
        evaluator$noninteraction_tape$pointer
    }),
    lapply(evaluators, function(evaluator) {
      evaluator$curvature_tapes[[approximation]]$pointer
    }),
    points, n_theta + seq_len(n_eta), population_positions, transform
  )
}

.nm_nested_objective <- function(context, approximation, eta_maxit, tolerance,
                                 initial_eta = NULL) {
  force(context); force(approximation)
  state <- new.env(parent = emptyenv())
  state$starts <- initial_eta %||% matrix(0, context$n_subjects, context$n_eta)
  state$key <- NULL
  state$value <- NULL
  state$modes <- NULL
  state$objective_calls <- 0L
  state$cache_hits <- 0L
  state$mode_iterations <- 0L
  state$mode_evaluations <- 0L
  objective <- function(parameters) {
    state$objective_calls <- state$objective_calls + 1L
    key <- c(parameters$theta, parameters$sigma, parameters$omega)
    if (!is.null(state$key) && identical(key, state$key)) {
      state$cache_hits <- state$cache_hits + 1L
      return(state$value)
    }
    modes <- .nm_subject_modes(
      context, parameters, starts = state$starts, maxit = eta_maxit,
      tolerance = tolerance, interaction = approximation != "foce",
      exact_hessian = approximation == "laplace"
    )
    if (any(vapply(modes, `[[`, integer(1), "convergence") != 0L)) return(Inf)
    state$mode_iterations <- state$mode_iterations + sum(vapply(
      modes, function(mode) as.integer(mode$iterations %||% 0L), integer(1)
    ))
    state$mode_evaluations <- state$mode_evaluations + sum(vapply(
      modes, function(mode) as.integer(mode$evaluations %||% 0L), integer(1)
    ))
    if (context$n_eta) {
      state$starts <- do.call(rbind, lapply(modes, `[[`, "par"))
    }
    value <- sum(vapply(modes, `[[`, numeric(1), "value"))
    if (approximation == "laplace") {
      value <- value + sum(vapply(modes, `[[`, numeric(1), "logdet"))
    } else if (approximation %in% c("foce", "focei")) {
      for (subject in seq_along(modes)) {
        value <- value + .nm_subject_curvature_logdet(
          context, context$subjects[[subject]], parameters,
          modes[[subject]]$par, approximation
        )
      }
    }
    value <- value + .nm_prior_nll(context$model, parameters)
    state$key <- key
    state$value <- value
    state$modes <- modes
    state$parameters <- parameters
    value
  }
  attr(objective, "state") <- state
  objective
}

.nm_fo_subject <- function(evaluator, model, theta, sigma, omega) {
  evaluator$fo_objective(theta, sigma, omega)$value
}

.nm_fo_subject_reference <- function(evaluator, model, theta, sigma, omega) {
  eta_columns <- length(theta) + seq_len(evaluator$n_eta)
  prediction <- evaluator$prediction(
    theta, rep(0, evaluator$n_eta), sigma, jacobian = TRUE,
    columns = eta_columns
  )
  observed_data <- evaluator$observation_data()
  observed_rows <- as.integer(attr(observed_data, "rows")) + 1L
  f <- prediction$value[observed_rows]
  dv <- observed_data$DV
  jacobian <- prediction$jacobian[observed_rows, , drop = FALSE]
  if (model$LIK_CONFIG$error == "exponential") {
    if (any(dv <= 0) || any(f <= 0)) return(Inf)
    residual <- log(dv) - log(f)
    jacobian <- jacobian / f
  } else residual <- dv - f
  dvid <- observed_data$DVID %||% rep(1L, length(dv))
  variance <- .nm_residual_variance(model, f, sigma, dvid)
  correlation <- diag(length(f))
  if (model$LIK_CONFIG$sigma_corr == "ar1" && length(f) > 1L) {
    rho <- .nm_ar1_rho(model, theta = theta, sigma = sigma)
    correlation <- outer(seq_along(f), seq_along(f), function(i, j) {
      rho^abs(i - j)
    })
  }
  if (length(model$LIK_CONFIG$residual_groups) && length(f) > 1L) {
    observed_dvid <- observed_data$DVID %||% rep(1L, length(dv))
    for (group in model$LIK_CONFIG$residual_groups) {
      group_correlation <- .nm_residual_group_value(group, theta, sigma)
      for (row in seq_along(dv)) {
        if (!observed_dvid[[row]] %in% group$dvid) next
        for (column in seq_along(dv)) {
          if (row == column || observed_data$.ID_INDEX[[row]] != observed_data$.ID_INDEX[[column]] ||
              observed_data$TIME[[row]] != observed_data$TIME[[column]] ||
              !observed_dvid[[column]] %in% group$dvid) next
          correlation[row, column] <- group_correlation[
            match(observed_dvid[[row]], group$dvid),
            match(observed_dvid[[column]], group$dvid)
          ]
        }
      }
    }
  }
  residual_covariance <- correlation * outer(sqrt(variance), sqrt(variance))
  marginal <- residual_covariance +
    jacobian %*% .nm_effect_covariance_evaluator(
      model, evaluator, omega
    ) %*% t(jacobian)
  pd <- .nm_positive_definite(marginal, "FO marginal covariance")
  as.numeric(pd$logdet + crossprod(residual, solve(pd$matrix, residual)))
}

.nm_fo_objective <- function(context, parameters) {
  if (is.null(context$parallel)) {
    values <- vapply(
      context$subjects, .nm_fo_subject, numeric(1),
      model = context$model, theta = parameters$theta,
      sigma = parameters$sigma, omega = parameters$omega
    )
  } else {
    pieces <- parallel::clusterApply(
      context$parallel$cluster, seq_along(context$parallel$chunks),
      function(index, parameters) {
        namespace <- asNamespace("LibeRation")
        subject_objective <- get(".nm_fo_subject", envir = namespace)
        state <- get(
          ".nm_parallel_worker_state", envir = namespace
        )()
        evaluators <- state$subjects
        model <- state$model
        vapply(
          evaluators, subject_objective, numeric(1), model = model,
          theta = parameters$theta, sigma = parameters$sigma,
          omega = parameters$omega
        )
      }, parameters = parameters
    )
    values <- unlist(pieces, use.names = FALSE)
  }
  sum(values) + .nm_prior_nll(context$model, parameters)
}

.nm_fo_collection_gradient <- function(evaluators, parameters) {
  if (!length(evaluators)) return(matrix(numeric(), 0L, 0L))
  # Shared FO tapes hold subject data as CppAD dynamic parameters. Select each
  # subject immediately before differentiating; collecting duplicate pointers
  # first would leave all rows on the final subject's dynamic values.
  do.call(rbind, lapply(evaluators, function(evaluator) {
    unname(evaluator$fo_objective(
      parameters$theta, parameters$sigma, parameters$omega,
      gradient = TRUE
    )$gradient)
  }))
}

.nm_fo_native_gradient <- function(context, parameters) {
  if (is.null(context$parallel)) {
    total <- colSums(.nm_fo_collection_gradient(context$subjects, parameters))
  } else {
    pieces <- parallel::clusterCall(
      context$parallel$cluster,
      function(parameters) {
        worker_state <- get(
          ".nm_parallel_worker_state", envir = asNamespace("LibeRation")
        )
        evaluators <- worker_state()$subjects
        collection <- get(".nm_fo_collection_gradient", envir = asNamespace("LibeRation"))
        colSums(collection(evaluators, parameters))
      }, parameters = parameters
    )
    total <- Reduce(`+`, pieces)
  }
  as.numeric(total) + .nm_prior_nll_native_gradient(context$model, parameters)
}

.nm_fo_outer_gradient <- function(context, map, parameters) {
  as.vector(.nm_fo_native_gradient(context, parameters) %*% map$jacobian(parameters))
}

.nm_assert_final_conditional_modes <- function(modes, method) {
  convergence <- if (length(modes)) {
    vapply(modes, function(mode) as.integer(mode$convergence %||% 0L), integer(1))
  } else integer()
  failed <- which(is.na(convergence) | convergence != 0L)
  if (length(failed)) {
    .nm_stop(
      method, " did not produce converged final conditional modes for subject",
      if (length(failed) == 1L) " " else "s ",
      paste(failed, collapse = ", "),
      ". The outer estimate has not been accepted because its reported ",
      "objective would otherwise depend on approximate ETA modes."
    )
  }
  invisible(convergence)
}

