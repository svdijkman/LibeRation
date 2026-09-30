# Iterative two-stage (ITS) estimation.
# Split from estimation-stochastic.R as a behaviour-preserving source move.

.nm_its_distribution <- function(context, parameters, starts, eta_maxit,
                                 tolerance) {
  gaussian_first_order <- context$model$LIK_CONFIG$error %in%
    c("additive", "proportional", "combined", "exponential", "power") &&
    identical(context$model$LIK_CONFIG$sigma_corr %||% "independent",
              "independent") &&
    !length(context$model$LIK_CONFIG$residual_groups)
  # Gaussian ITS covariance is built from the prediction Jacobian and OMEGA.
  # The exact conditional-objective Hessian is used only by the non-Gaussian
  # fallback, so avoid a redundant Reverse(2) sweep in the common path.
  modes <- .nm_subject_modes(
    context, parameters, starts = starts, maxit = eta_maxit,
    tolerance = tolerance, interaction = TRUE,
    exact_hessian = !gaussian_first_order
  )
  if (any(vapply(modes, `[[`, integer(1), "convergence") != 0L)) {
    .nm_stop("ITS conditional-mode calculation failed.")
  }
  native_covariance <- NULL
  native_error <- NULL
  if (gaussian_first_order && is.null(context$parallel) && context$n_eta > 0L &&
      isTRUE(getOption("LibeRation.its_native_covariance", TRUE))) {
    native_covariance <- tryCatch(
      .liberation_its_gaussian_covariance(
        context$engine$pointer,
        lapply(context$subjects, function(evaluator) {
          evaluator$prediction_tape$pointer
        }),
        lapply(context$subjects, function(evaluator) evaluator$data_input()),
        parameters$theta,
        do.call(rbind, lapply(modes, `[[`, "par")),
        parameters$sigma, parameters$omega
      ),
      error = function(error) {
        native_error <<- conditionMessage(error)
        NULL
      }
    )
  }
  conditional_covariance <- if (!is.null(native_covariance)) {
    Map(function(covariance, mode, evaluator) {
      eta_names <- names(mode$par)
      if (!length(eta_names)) {
        eta_names <- rownames(.nm_effect_covariance_evaluator(
          context$model, evaluator, parameters$omega
        ))
      }
      if (!length(eta_names)) {
        eta_names <- paste0("ETA_1_", seq_len(nrow(covariance)))
      }
      if (length(eta_names) == nrow(covariance)) {
        dimnames(covariance) <- list(eta_names, eta_names)
      }
      covariance
    }, native_covariance$covariance, modes, context$subjects)
  } else Map(function(mode, evaluator) {
    dimension <- length(mode$par)
    if (!dimension) return(matrix(numeric(), 0L, 0L))
    if (!gaussian_first_order) {
      return(.nm_positive_definite(
        2 * solve(mode$hessian), "ITS approximate conditional covariance"
      )$matrix)
    }
    eta_columns <- length(parameters$theta) + seq_len(context$n_eta)
    prediction <- evaluator$prediction(
      parameters$theta, mode$par, parameters$sigma,
      jacobian = TRUE, columns = eta_columns
    )
    observations <- evaluator$observation_data()
    rows <- as.integer(attr(observations, "rows")) + 1L
    f <- prediction$value[rows]
    jacobian <- prediction$jacobian[rows, , drop = FALSE]
    dvid <- observations$DVID %||% rep(1L, length(rows))
    variance <- .nm_residual_variance(
      context$model, f, parameters$sigma, dvid
    )
    omega_inverse <- solve(.nm_effect_covariance_evaluator(
      context$model, evaluator, parameters$omega
    ))
    curvature <- 2 * crossprod(jacobian / sqrt(variance)) + 2 * omega_inverse
    curvature <- .nm_positive_definite(
      curvature, "ITS first-order conditional curvature"
    )$matrix
    2 * solve(curvature)
  }, modes, context$subjects)
  grids <- Map(function(mode, covariance) {
    dimension <- length(mode$par)
    if (!dimension) return(matrix(numeric(), 1L, 0L))
    covariance <- .nm_positive_definite(
      covariance, "ITS approximate conditional covariance"
    )$matrix
    root <- t(chol(covariance))
    # The symmetric sigma grid reproduces the conditional mean and covariance
    # exactly. Consequently it evaluates the first-order (quadratic)
    # conditional expectation used by ITS without Monte-Carlo noise.
    offsets <- sqrt(dimension) * rbind(t(root), -t(root))
    result <- sweep(offsets, 2L, mode$par, `+`)
    eta_names <- rownames(covariance)
    if (length(eta_names) == ncol(result)) {
      dimnames(result) <- list(rep(eta_names, 2L), eta_names)
    }
    result
  }, modes, conditional_covariance)
  weights <- lapply(grids, function(grid) rep(1 / nrow(grid), nrow(grid)))
  list(
    modes = modes, eta = grids, weights = weights,
    covariance = conditional_covariance,
    covariance_backend = if (!is.null(native_covariance)) {
      native_covariance$backend
    } else "r-subject-first-order",
    covariance_native_error = native_error
  )
}

.nm_its_omega_sufficient <- function(context, modes, conditional_covariance) {
  dimension <- as.integer(context$model$n_eta)
  iov <- as.integer(context$model$LIK_CONFIG$iov %||% 0L)
  covariance <- matrix(0, dimension, dimension)
  if (!iov) {
    for (subject in seq_along(modes)) {
      eta <- as.numeric(modes[[subject]]$par)
      covariance <- covariance + tcrossprod(eta) +
        conditional_covariance[[subject]]
    }
    covariance <- covariance / max(length(modes), 1L)
  } else {
    between <- dimension - iov
    if (between > 0L) {
      index <- seq_len(between)
      for (subject in seq_along(modes)) {
        eta <- as.numeric(modes[[subject]]$par[index])
        covariance[index, index] <- covariance[index, index] +
          tcrossprod(eta) + conditional_covariance[[subject]][index, index,
                                                              drop = FALSE]
      }
      covariance[index, index] <- covariance[index, index] /
        max(length(modes), 1L)
    }
    target <- between + seq_len(iov)
    occasions <- 0L
    for (subject in seq_along(modes)) {
      eta <- as.numeric(modes[[subject]]$par)
      n_occasions <- as.integer((length(eta) - between) / iov)
      for (occasion in seq_len(n_occasions)) {
        source <- between + (occasion - 1L) * iov + seq_len(iov)
        covariance[target, target] <- covariance[target, target] +
          tcrossprod(eta[source]) +
          conditional_covariance[[subject]][source, source, drop = FALSE]
        occasions <- occasions + 1L
      }
    }
    covariance[target, target] <- covariance[target, target] /
      max(occasions, 1L)
  }
  result <- vapply(seq_len(nrow(context$model$OMEGAS)), function(index) {
    covariance[
      context$model$OMEGAS$ROW[[index]], context$model$OMEGAS$COL[[index]]
    ]
  }, numeric(1))
  diagonal <- context$model$OMEGAS$ROW == context$model$OMEGAS$COL
  result[diagonal] <- pmax(result[diagonal], 1e-12)
  result
}

.nm_est_its <- function(context, map, maxit, eta_maxit, tolerance, trace,
                        print_every = 0L, optimizer_backend = "auto",
                        its_mstep_maxit = NULL,
                        its_mstep_schedule = c("auto", "fixed", "progressive"),
                        its_acceleration = c("auto", "none", "aitken"),
                        its_eta_schedule = c("auto", "fixed", "progressive"),
                        its_eta_tolerance_multiplier = 100) {
  its_mstep_schedule <- match.arg(its_mstep_schedule)
  its_acceleration <- match.arg(its_acceleration)
  its_eta_schedule <- match.arg(its_eta_schedule)
  resolved_mstep_schedule <- if (its_mstep_schedule == "auto") {
    if (.nm_liber_optimized(context)) "progressive" else "fixed"
  } else its_mstep_schedule
  resolved_acceleration <- if (its_acceleration == "auto") {
    if (.nm_liber_optimized(context)) "aitken" else "none"
  } else its_acceleration
  resolved_eta_schedule <- if (its_eta_schedule == "auto") {
    if (.nm_liber_optimized(context)) "progressive" else "fixed"
  } else its_eta_schedule
  if (!.nm_liber_optimized(context)) resolved_acceleration <- "none"
  if (!.nm_liber_optimized(context)) resolved_eta_schedule <- "fixed"
  n_iter <- max(1L, as.integer(maxit))
  if (is.null(its_mstep_maxit)) {
    its_mstep_maxit <- if (.nm_liber_optimized(context)) 10L else 1L
  }
  its_mstep_maxit <- as.integer(its_mstep_maxit)
  if (length(its_mstep_maxit) != 1L || is.na(its_mstep_maxit) ||
      its_mstep_maxit < 1L) {
    .nm_stop("`its_mstep_maxit` must be one positive integer.")
  }
  if (length(its_eta_tolerance_multiplier) != 1L ||
      !is.finite(its_eta_tolerance_multiplier) ||
      its_eta_tolerance_multiplier < 1) {
    .nm_stop("`its_eta_tolerance_multiplier` must be finite and at least one.")
  }
  parameters <- map$decode(map$start)
  starts <- matrix(0, context$n_subjects, context$n_eta)
  objective_trace <- numeric(n_iter)
  parameter_trace <- if (length(map$start)) {
    matrix(NA_real_, n_iter, length(map$start))
  } else matrix(numeric(), n_iter, 0L)
  total_evaluations <- total_gradient_evaluations <- total_mstep_iterations <- 0L
  optimizer <- NULL
  modes <- vector("list", context$n_subjects)
  completed <- 0L
  phases <- .nm_stochastic_phase_timer()
  weighted_context <- .nm_weighted_eta_context(context)
  its_native_model <- context$model
  its_native_model$THETAS$Value <- parameters$theta
  its_native_model$SIGMAS$Value <- parameters$sigma
  its_native_model$OMEGAS$Value <- parameters$omega
  if (length(map$omega_free)) its_native_model$OMEGAS$FIX[] <- TRUE
  its_native_map <- .nm_outer_map(its_native_model)
  its_native_enabled <- isTRUE(getOption(
    "LibeRation.its_native_mstep", TRUE
  )) && .nm_liber_optimized(context) && !is.null(weighted_context) &&
    is.null(context$parallel) && optimizer_backend %in% c("auto", "native")
  its_native_state <- NULL
  its_native_attempts <- 0L
  its_native_successes <- 0L
  its_native_fallbacks <- 0L
  its_native_fallback_reason <- NULL
  mstep_history <- list()
  acceleration_attempts <- 0L
  acceleration_accepts <- 0L
  eta_tolerance_trace <- numeric(n_iter)
  for (iteration in seq_len(n_iter)) {
    current_eta_tolerance <- if (resolved_eta_schedule == "progressive") {
      min(1e-3, tolerance *
        its_eta_tolerance_multiplier^((n_iter - iteration) / max(n_iter - 1L, 1L)))
    } else tolerance
    eta_tolerance_trace[[iteration]] <- current_eta_tolerance
    expectation_state <- phases$time("conditional_distribution", {
      .nm_its_distribution(
        context, parameters, starts, eta_maxit, current_eta_tolerance
      )
    })
    modes <- expectation_state$modes
    if (context$n_eta) {
      starts <- do.call(rbind, lapply(modes, `[[`, "par"))
    }
    # NONMEM-style ITS is not Gaussian quadrature EM. THETA and SIGMA are
    # updated at the conditional modes; conditional variances enter the
    # OMEGA second moment. This is the documented approximation that makes ITS
    # approach, but not equal, the FOCE linearized objective.
    mstep_model <- context$model
    mstep_model$THETAS$Value <- parameters$theta
    mstep_model$SIGMAS$Value <- parameters$sigma
    mstep_model$OMEGAS$Value <- parameters$omega
    if (length(map$omega_free)) mstep_model$OMEGAS$FIX[] <- TRUE
    iteration_map <- .nm_outer_map(mstep_model)
    mode_grid <- lapply(modes, function(mode) matrix(mode$par, nrow = 1L))
    expectation <- phases$time("expectation_setup", {
      .nm_complete_data_expectation(
        context, iteration_map, mode_grid, rep(list(1), context$n_subjects),
        native_context = weighted_context
      )
    })
    current_mstep_maxit <- if (resolved_mstep_schedule == "progressive") {
      min(its_mstep_maxit, max(1L, as.integer(ceiling(
        its_mstep_maxit * iteration / n_iter
      ))))
    } else its_mstep_maxit
    native_result <- NULL
    if (its_native_enabled && length(its_native_map$start)) {
      its_native_attempts <- its_native_attempts + 1L
      native_result <- phases$time("mstep", tryCatch(
        .nm_native_weighted_mstep(
          context, parameters, weighted_context, its_native_map,
          current_mstep_maxit, tolerance,
          if (trace > 1L) trace else 0L, its_native_state
        ),
        error = identity
      ))
      if (inherits(native_result, "error") ||
          !is.finite(native_result$value %||% NA_real_) ||
          any(!is.finite(native_result$theta %||% NA_real_)) ||
          any(!is.finite(native_result$sigma %||% NA_real_))) {
        its_native_fallbacks <- its_native_fallbacks + 1L
        its_native_fallback_reason <- if (inherits(native_result, "error")) {
          conditionMessage(native_result)
        } else "native ITS M-step returned non-finite values"
        native_result <- NULL
        its_native_enabled <- FALSE
      } else {
        its_native_successes <- its_native_successes + 1L
        its_native_state <- native_result$optimizer_state
        native_result$optimizer_state <- NULL
      }
    }
    if (!is.null(native_result)) {
      optimizer <- native_result
      candidate <- list(
        theta = as.numeric(native_result$theta),
        sigma = as.numeric(native_result$sigma),
        omega = as.numeric(native_result$omega)
      )
      optimizer$par <- iteration_map$encode(candidate)
    } else if (length(iteration_map$start)) {
      optimizer <- phases$time("mstep", {
        .nm_outer_optim(
          iteration_map, expectation$objective, current_mstep_maxit, tolerance,
          if (trace > 1L) trace else 0L, 0L,
          gradient = expectation$gradient, optimizer_backend = optimizer_backend
        )
      })
      candidate <- iteration_map$decode(optimizer$par)
      candidate_point <- iteration_map$encode(candidate)
      mstep_history[[length(mstep_history) + 1L]] <- candidate_point
      if (resolved_acceleration == "aitken" &&
          length(mstep_history) >= 3L && length(candidate_point)) {
        acceleration_attempts <- acceleration_attempts + 1L
        recent <- tail(mstep_history, 3L)
        first_difference <- recent[[2L]] - recent[[1L]]
        second_difference <- recent[[3L]] - recent[[2L]]
        curvature_difference <- second_difference - first_difference
        denominator <- sum(curvature_difference^2)
        if (is.finite(denominator) && denominator > 1e-16) {
          factor <- -sum(second_difference * curvature_difference) /
            denominator
          factor <- min(2, max(0, factor))
          accelerated_point <- candidate_point + factor * second_difference
          if (iteration_map$in_bounds(accelerated_point)) {
            accelerated <- iteration_map$decode(accelerated_point)
            base_value <- expectation$objective(candidate)
            accelerated_value <- tryCatch(
              expectation$objective(accelerated), error = function(error) Inf
            )
            if (is.finite(accelerated_value) &&
                accelerated_value <= base_value) {
              candidate <- accelerated
              optimizer$par <- accelerated_point
              optimizer$value <- accelerated_value
              mstep_history[[length(mstep_history)]] <- accelerated_point
              acceleration_accepts <- acceleration_accepts + 1L
            }
          }
        }
      }
    } else {
      optimizer <- list(
        par = numeric(), value = expectation$objective(parameters),
        convergence = 0L, message = "ITS moment-only M-step",
        counts = c(`function` = 1L, gradient = 0L), iterations = 0L,
        objective_evaluations = 1L, gradient_evaluations = 0L,
        backend = "its-moment-update"
      )
      candidate <- parameters
    }
    previous <- map$encode(parameters)
    parameters$theta <- candidate$theta
    parameters$sigma <- candidate$sigma
    if (length(map$omega_free)) {
      omega_prior <- context$model$LIK_CONFIG$priors
      omega_prior <- !is.null(omega_prior) && nrow(omega_prior) &&
        any(startsWith(omega_prior$parameter, "OMEGA"))
      if (!omega_prior) {
        sufficient <- .nm_its_omega_sufficient(
          context, modes, expectation_state$covariance
        )
        parameters$omega[map$omega_free] <- sufficient[map$omega_free]
      } else {
        omega_model <- context$model
        omega_model$THETAS$Value <- parameters$theta
        omega_model$SIGMAS$Value <- parameters$sigma
        omega_model$OMEGAS$Value <- parameters$omega
        omega_model$THETAS$FIX[] <- TRUE
        omega_model$SIGMAS$FIX[] <- TRUE
        omega_map <- .nm_outer_map(omega_model)
        omega_expectation <- .nm_complete_data_expectation(
          context, omega_map, expectation_state$eta, expectation_state$weights,
          native_context = weighted_context
        )
        omega_optimizer <- .nm_outer_optim(
          omega_map, omega_expectation$objective, its_mstep_maxit, tolerance,
          if (trace > 1L) trace else 0L, 0L,
          gradient = omega_expectation$gradient,
          optimizer_backend = optimizer_backend
        )
        parameters$omega <- omega_map$decode(omega_optimizer$par)$omega
        optimizer$objective_evaluations <-
          as.integer(optimizer$objective_evaluations %||% 0L) +
          as.integer(omega_optimizer$objective_evaluations %||% 0L)
        optimizer$gradient_evaluations <-
          as.integer(optimizer$gradient_evaluations %||% 0L) +
          as.integer(omega_optimizer$gradient_evaluations %||% 0L)
        optimizer$iterations <- as.integer(optimizer$iterations %||% 0L) +
          as.integer(omega_optimizer$iterations %||% 0L)
      }
    }
    point <- map$encode(parameters)
    objective_trace[[iteration]] <- expectation$objective(parameters)
    if (length(point)) parameter_trace[iteration, ] <- point
    total_evaluations <- total_evaluations +
      as.integer(optimizer$objective_evaluations %||% 0L)
    total_gradient_evaluations <- total_gradient_evaluations +
      as.integer(optimizer$gradient_evaluations %||% 0L)
    total_mstep_iterations <- total_mstep_iterations +
      as.integer(optimizer$iterations %||% 0L)
    completed <- iteration
    if (print_every > 0L && iteration %% print_every == 0L) {
      cat(sprintf("[LibeRation] ITS ITERATION %d Q %.10g\n", iteration,
                  objective_trace[[iteration]]))
      try(flush(stdout()), silent = TRUE)
    }
    if (length(point) && iteration > 1L &&
        max(abs(point - previous) / (1 + abs(previous))) <= tolerance) break
  }
  expectation_state <- .nm_its_distribution(
    context, parameters, starts, eta_maxit, tolerance
  )
  modes <- expectation_state$modes
  final_mode_grid <- lapply(modes, function(mode) matrix(mode$par, nrow = 1L))
  final_expectation <- .nm_complete_data_expectation(
    context, map, final_mode_grid, rep(list(1), context$n_subjects),
    native_context = weighted_context
  )
  final_objective <- final_expectation$objective(parameters)
  parameter_converged <- completed < n_iter
  optimizer$convergence <- 0L
  optimizer$message <- if (completed < n_iter) "ITS parameter convergence reached" else
    "ITS iterations completed"
  optimizer$iterations <- completed
  optimizer$objective_evaluations <- total_evaluations
  optimizer$gradient_evaluations <- total_gradient_evaluations
  optimizer$mstep_iterations <- total_mstep_iterations
  .nm_fit_result(
    context, "ITS", parameters, final_objective, modes, optimizer,
    diagnostics = list(
      eta_convergence = vapply(modes, `[[`, integer(1), "convergence"),
      description = paste0(
        "iterative two-stage EM using conditional modes and first-order ",
        "approximate conditional variances"
      ),
      estimator_identity = paste0(
        "ITS conditional-mode/variance approximation with iterative ",
        "single-step population updates"
      ),
      theta_sigma_update = "complete-data gradient at conditional modes",
      omega_update = "conditional-mode second moment plus conditional variance",
      mstep_maxit = its_mstep_maxit,
      mstep_schedule = resolved_mstep_schedule,
      native_mstep = list(
        eligible = .nm_liber_optimized(context) && !is.null(weighted_context) &&
          is.null(context$parallel),
        enabled_at_completion = its_native_enabled,
        attempts = its_native_attempts,
        successes = its_native_successes,
        fallbacks = its_native_fallbacks,
        fallback_reason = its_native_fallback_reason,
        persistent_optimizer_state = !is.null(its_native_state),
        compatibility_policy = if (.nm_liber_optimized(context)) {
          "native weighted M-step enabled without changing ITS E-step semantics"
        } else {
          "R coordinator retained for NONMEM-compatible numerical policy"
        }
      ),
      acceleration = list(
        method = resolved_acceleration, attempts = acceleration_attempts,
        accepted = acceleration_accepts,
        monotonic_safeguard = TRUE
      ),
      eta_tolerance_schedule = list(
        method = resolved_eta_schedule,
        multiplier = its_eta_tolerance_multiplier,
        trace = eta_tolerance_trace[seq_len(completed)],
        final_tolerance = tolerance
      ),
      phase_timing = phases$snapshot(),
      expectation_backend = final_expectation$telemetry(),
      conditional_covariance_backend = expectation_state$covariance_backend,
      conditional_covariance_native_error =
        expectation_state$covariance_native_error,
      parameter_converged = parameter_converged,
      iteration_limit_reached = !parameter_converged,
      objective_trace = objective_trace[seq_len(completed)],
      parameter_trace = parameter_trace[seq_len(completed), , drop = FALSE],
      population_gradient = "exact CppAD gradient of the fixed ITS expectation"
    )
  )
}

