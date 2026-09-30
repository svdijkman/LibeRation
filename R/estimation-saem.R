# Stochastic approximation expectation-maximisation (SAEM).
# Split from estimation-stochastic.R as a behaviour-preserving source move.

.nm_stochastic_eta_context <- function(context,
                                       allow_compatibility = FALSE) {
  if ((!isTRUE(allow_compatibility) && !.nm_liber_optimized(context)) ||
      !is.null(context$parallel) ||
      !isTRUE(getOption(
        "LibeRation.stochastic_persistent_context", TRUE
      ))) {
    return(NULL)
  }
  native_threads <- as.integer(context$native_subject_threads %||% 1L)
  fused_values <- .nm_liber_optimized(context) && isTRUE(getOption(
    "LibeRation.fused_advan_stochastic", TRUE
  )) && (
    native_threads > 1L ||
      isTRUE(getOption("LibeRation.fused_advan_serial", TRUE))
  )
  tryCatch(
    .liberation_stochastic_eta_context_create(
      context$engine$pointer,
      lapply(context$subjects, function(evaluator) {
        evaluator$objective_tape$pointer
      }),
      lapply(context$subjects, function(evaluator) evaluator$data_input()),
      length(context$model$THETAS$Value), context$n_eta,
      length(context$model$SIGMAS$Value),
      length(context$model$OMEGAS$Value), isTRUE(context$model$USE_ODE),
      as.numeric(context$model$THETAS$Value),
      as.numeric(context$model$SIGMAS$Value),
      as.numeric(context$model$OMEGAS$Value),
      as.numeric(getOption("LibeRation.tape_guard_radius", 0.5)),
      fused_values, native_threads
    ),
    error = function(error) {
      warning(
        "Persistent stochastic context unavailable: ",
        conditionMessage(error), call. = FALSE
      )
      NULL
    }
  )
}

.nm_saem_conditional_components <- function(context, parameters, eta,
                                             persistent = NULL) {
  if (!is.null(persistent)) {
    subject <- .liberation_stochastic_eta_context_eval(
      persistent, parameters$theta, eta, parameters$sigma, parameters$omega
    )
  } else if (is.null(context$parallel)) {
    subject <- .nm_objective_collection(context$subjects, parameters, eta)
  } else {
    eta_chunks <- lapply(
      context$parallel$chunks, function(rows) eta[rows, , drop = FALSE]
    )
    pieces <- parallel::clusterApply(
      context$parallel$cluster, seq_along(context$parallel$chunks),
      function(index, eta_chunks, parameters) {
        namespace <- asNamespace("LibeRation")
        evaluators <- get(".nm_parallel_worker_state", envir = namespace)()$subjects
        worker_eta <- eta_chunks[[index]]
        collection <- get(".nm_objective_collection", envir = asNamespace("LibeRation"))
        collection(evaluators, parameters, worker_eta)
      }, eta_chunks = eta_chunks, parameters = parameters
    )
    subject <- unlist(pieces, use.names = FALSE)
  }
  prior <- .nm_prior_nll(context$model, parameters)
  list(
    value = sum(subject) + prior,
    subject = as.numeric(subject),
    prior = prior
  )
}

.nm_saem_conditional <- function(context, parameters, eta,
                                 persistent = NULL) {
  .nm_saem_conditional_components(
    context, parameters, eta, persistent = persistent
  )$value
}

.nm_fsaem_proposal <- function(context, parameters, starts, eta_maxit,
                               tolerance, persistent = NULL) {
  if (!is.null(persistent)) {
    proposal <- .liberation_stochastic_eta_context_laplace_proposal(
      persistent, parameters$theta, starts, parameters$sigma,
      parameters$omega, as.integer(eta_maxit), tolerance
    )
    proposal$anchor <- c(
      parameters$theta, parameters$sigma, parameters$omega
    )
    return(proposal)
  }
  modes <- .nm_subject_modes(
    context, parameters, starts = starts, maxit = eta_maxit,
    tolerance = tolerance, interaction = TRUE, exact_hessian = TRUE
  )
  convergence <- vapply(modes, `[[`, integer(1), "convergence")
  if (any(convergence != 0L)) {
    .nm_stop(
      "Laplace independence proposal failed for ",
      sum(convergence != 0L), " subject(s)."
    )
  }
  proposal <- lapply(modes, function(mode) {
    hessian_state <- .nm_positive_definite(
      mode$hessian, "f-SAEM conditional curvature"
    )
    hessian <- hessian_state$matrix
    covariance_state <- .nm_positive_definite(
      2 * chol2inv(chol(hessian)), "f-SAEM proposal covariance"
    )
    covariance <- covariance_state$matrix
    root <- chol(covariance)
    list(
      root = t(root),
      precision = if ((covariance_state$jitter %||% 0) == 0) {
        hessian / 2
      } else chol2inv(root)
    )
  })
  list(
    modes = do.call(rbind, lapply(modes, `[[`, "par")),
    roots = lapply(proposal, `[[`, "root"),
    precisions = lapply(proposal, `[[`, "precision"),
    mode_iterations = sum(vapply(
      modes, function(mode) as.integer(mode$iterations %||% 0L), integer(1)
    )),
    mode_evaluations = sum(vapply(
      modes, function(mode) as.integer(mode$evaluations %||% 0L), integer(1)
    )),
    anchor = c(parameters$theta, parameters$sigma, parameters$omega),
    backend = "r-coordinated-laplace-proposal"
  )
}

.nm_saem_conditional_gradient <- function(context, map, parameters, eta) {
  native <- .nm_conditional_native_gradient(
    context, parameters, eta, interaction = TRUE
  )
  as.vector(native %*% map$jacobian(parameters))
}

.nm_saem_paired_conditional <- function(context, map, eta) {
  cache <- new.env(parent = emptyenv())
  cache$key <- NULL
  cache$result <- NULL
  cache$evaluations <- 0L
  cache$hits <- 0L
  paired <- is.null(context$parallel) && isTRUE(getOption(
    "LibeRation.saem_paired_value_gradient", TRUE
  ))
  persistent_error <- NULL
  persistent <- NULL
  persistent_requested <- paired && !isTRUE(context$model$USE_ODE) &&
    isTRUE(getOption(
    "LibeRation.saem_persistent_fixed_eta", TRUE
  ))
  if (persistent_requested) {
    persistent <- tryCatch(
      .liberation_saem_fixed_eta_context_create(
        lapply(context$subjects, function(evaluator) {
          evaluator$objective_tape$pointer
        }),
        eta, length(context$model$THETAS$Value),
        length(context$model$SIGMAS$Value),
        length(context$model$OMEGAS$Value)
      ),
      error = function(error) {
        persistent_error <<- conditionMessage(error)
        NULL
      }
    )
  }
  aggregate <- !is.null(persistent) &&
    isTRUE(getOption("LibeRation.saem_native_aggregate", TRUE))
  evaluate <- function(parameters) {
    key <- map$encode(parameters)
    if (!is.null(cache$key) && identical(key, cache$key)) {
      cache$hits <- cache$hits + 1L
      return(cache$result)
    }
    cache$evaluations <- cache$evaluations + 1L
    if (paired) {
      collection <- if (aggregate) {
        .liberation_saem_fixed_eta_context_eval_aggregate(
          persistent, parameters$theta, parameters$sigma, parameters$omega
        )
      } else if (!is.null(persistent)) {
        .liberation_saem_fixed_eta_context_eval(
          persistent, parameters$theta, parameters$sigma, parameters$omega
        )
      } else {
        .nm_objective_collection_value_gradient(
          context$subjects, parameters, eta, interaction = TRUE
        )
      }
      total <- if (aggregate) collection$gradient else
        colSums(collection$gradient)
      n_theta <- length(parameters$theta)
      n_sigma <- length(parameters$sigma)
      n_omega <- length(parameters$omega)
      population_positions <- c(
        seq_len(n_theta),
        n_theta + context$n_eta + seq_len(n_sigma),
        n_theta + context$n_eta + n_sigma + seq_len(n_omega)
      )
      native <- as.numeric(total[population_positions]) +
        .nm_prior_nll_native_gradient(context$model, parameters)
      result <- list(
        value = if (aggregate) collection$value else sum(collection$value),
        gradient = as.vector(native %*% map$jacobian(parameters))
      )
      result$value <- result$value + .nm_prior_nll(context$model, parameters)
    } else {
      result <- list(
        value = .nm_saem_conditional(context, parameters, eta),
        gradient = .nm_saem_conditional_gradient(context, map, parameters, eta)
      )
    }
    cache$key <- key
    cache$result <- result
    result
  }
  list(
    objective = function(parameters) evaluate(parameters)$value,
    gradient = function(parameters) evaluate(parameters)$gradient,
    telemetry = function() list(
      paired = paired, evaluations = cache$evaluations,
      cache_hits = cache$hits, persistent = !is.null(persistent),
      native_aggregate = aggregate,
      persistent_requested = persistent_requested,
      persistent_error = persistent_error
    )
  )
}

.nm_saem_q_state <- function(context, max_support = 0L,
                             prune_tolerance = 0) {
  state <- new.env(parent = emptyenv())
  state$samples <- list()
  state$weights <- numeric()
  state$native <- .nm_weighted_eta_context(
    context, reduced_population_tape = max_support > 0L
  )
  state$native_error <- NULL
  update <- function(eta, gamma) {
    eta <- as.matrix(eta)
    if (!identical(dim(eta), c(context$n_subjects, context$n_eta)) ||
        any(!is.finite(eta)) || !is.finite(gamma) || gamma <= 0 || gamma > 1) {
      .nm_stop("SAEM stochastic-approximation state is invalid.")
    }
    if (!is.null(state$native)) {
      update_error <- tryCatch({
        .liberation_weighted_eta_context_update(
          state$native, eta, gamma, as.integer(max_support),
          as.numeric(prune_tolerance)
        )
        NULL
      }, error = identity)
      if (inherits(update_error, "error")) {
        .nm_stop(
          "Native SAEM stochastic-approximation update failed: ",
          conditionMessage(update_error)
        )
      }
      return(invisible(NULL))
    }
    if (gamma >= 1 - .Machine$double.eps || !length(state$samples)) {
      state$samples <- list(eta)
      state$weights <- 1
    } else {
      state$weights <- (1 - gamma) * state$weights
      state$samples[[length(state$samples) + 1L]] <- eta
      state$weights <- c(state$weights, gamma)
      state$weights <- state$weights / sum(state$weights)
    }
    invisible(NULL)
  }
  mean_eta <- function() {
    if (!is.null(state$native)) {
      return(.liberation_weighted_eta_context_mean(state$native))
    }
    if (!length(state$samples)) {
      return(matrix(0, context$n_subjects, context$n_eta))
    }
    Reduce(`+`, Map(`*`, state$samples, state$weights))
  }
  grids <- function() {
    if (!length(state$samples)) .nm_stop("SAEM Q state is empty.")
    eta <- lapply(seq_len(context$n_subjects), function(subject) {
      do.call(rbind, lapply(state$samples, function(sample) {
        matrix(sample[subject, ], nrow = 1L)
      }))
    })
    list(eta = eta, weights = rep(list(state$weights), context$n_subjects))
  }
  recenter <- function(mu, previous, current) {
    if (!isTRUE(mu$active)) return(invisible(NULL))
    if (!is.null(state$native)) {
      adjustment <- .nm_mu_recenter_eta(
        mu, previous, current,
        matrix(0, context$n_subjects, context$n_eta)
      )
      .liberation_weighted_eta_context_recenter(state$native, adjustment)
      return(invisible(NULL))
    }
    if (!length(state$samples)) return(invisible(NULL))
    state$samples <- lapply(state$samples, function(sample) {
      .nm_mu_recenter_eta(mu, previous, current, sample)
    })
    invisible(NULL)
  }
  list(
    update = update, mean_eta = mean_eta, grids = grids, recenter = recenter,
    count = function() if (!is.null(state$native)) {
      length(.liberation_weighted_eta_context_weights(state$native))
    } else length(state$samples),
    weights = function() if (!is.null(state$native)) {
      as.numeric(.liberation_weighted_eta_context_weights(state$native))
    } else state$weights,
    samples = function() state$samples,
    native = function() state$native,
    sigma = function(parameters) if (!is.null(state$native)) {
      .liberation_weighted_eta_context_sigma(
        state$native, context$engine$pointer, context$data,
        parameters$theta, parameters$sigma
      )
    } else NULL,
    omega = function() if (!is.null(state$native)) {
      .liberation_weighted_eta_context_omega(
        state$native, as.integer(context$model$n_eta),
        as.integer(context$model$LIK_CONFIG$iov),
        as.integer(context$model$OMEGAS$ROW),
        as.integer(context$model$OMEGAS$COL)
      )
    } else NULL,
    telemetry = function() if (!is.null(state$native)) {
      .liberation_weighted_eta_context_telemetry(state$native)
    } else list(backend = "r-retained-support", support = length(state$samples))
  )
}

.nm_saem_q_expectation <- function(context, map, q_state) {
  if (!is.null(q_state$native())) {
    return(.nm_native_weighted_expectation(context, map, q_state$native()))
  }
  grid <- q_state$grids()
  .nm_complete_data_expectation(context, map, grid$eta, grid$weights)
}

.nm_saem_sigma_expectation <- function(context, parameters, q_state) {
  native <- q_state$sigma(parameters)
  if (!is.null(native)) return(native)
  values <- lapply(q_state$samples(), function(eta) {
    .nm_saem_sigma_sufficient(context, parameters, eta)
  })
  if (!length(values) || any(vapply(values, is.null, logical(1)))) return(NULL)
  weights <- q_state$weights()
  matrix_values <- do.call(rbind, values)
  if (!identical(context$model$LIK_CONFIG$sigma_parameterization, "variance")) {
    matrix_values <- matrix_values^2
  }
  sufficient <- colSums(matrix_values * weights)
  if (!identical(context$model$LIK_CONFIG$sigma_parameterization, "variance")) {
    sufficient <- sqrt(pmax(sufficient, 0))
  }
  sufficient
}

.nm_saem_sigma_from_q_gradient <- function(context, parameters,
                                           native_gradient) {
  if (is.null(native_gradient) ||
      !context$model$LIK_CONFIG$error %in%
        c("additive", "proportional", "exponential") ||
      !identical(context$model$LIK_CONFIG$sigma_corr %||% "independent",
                 "independent") ||
      length(context$model$LIK_CONFIG$residual_groups)) {
    return(NULL)
  }
  n_theta <- length(parameters$theta)
  n_sigma <- length(parameters$sigma)
  positions <- n_theta + context$n_eta + seq_len(n_sigma)
  if (length(native_gradient) < max(positions) ||
      any(!is.finite(native_gradient[positions]))) {
    return(NULL)
  }
  observed <- context$data$EVID == 0L & context$data$MDV == 0L &
    is.finite(context$data$DV)
  if (identical(context$model$LIK_CONFIG$error, "exponential")) {
    observed <- observed & context$data$DV > 0
  }
  dvid <- if ("DVID" %in% names(context$data)) {
    pmax(as.integer(context$data$DVID), 1L)
  } else rep(1L, nrow(context$data))
  count <- tabulate(dvid[observed], nbins = n_sigma)
  value <- as.numeric(parameters$sigma)
  variance_parameterization <- identical(
    context$model$LIK_CONFIG$sigma_parameterization, "variance"
  )
  for (response in seq_len(n_sigma)) {
    observations <- count[[response]]
    current <- parameters$sigma[[response]]
    gradient <- native_gradient[[positions[[response]]]]
    if (!observations || !is.finite(current) || current <= 0) next
    # For a simple Gaussian response contribution, the complete-data
    # objective is N log(v) + SSE/v (or 2N log(s) + SSE/s^2).  Solving its
    # exact derivative for SSE/N recovers the same closed-form SAEM update
    # without simulating every retained stochastic-approximation support.
    variance <- if (variance_parameterization) {
      current - gradient * current^2 / observations
    } else {
      current^2 - gradient * current^3 / (2 * observations)
    }
    if (is.finite(variance) && variance > 0) {
      value[[response]] <- if (variance_parameterization) {
        variance
      } else sqrt(variance)
    }
  }
  value
}

.nm_saem_omega_expectation <- function(context, q_state) {
  native <- q_state$omega()
  if (!is.null(native)) return(native)
  values <- lapply(q_state$samples(), function(eta) {
    .nm_saem_omega_sufficient(context, eta)
  })
  if (!length(values)) return(NULL)
  colSums(do.call(rbind, values) * q_state$weights())
}

.nm_saem_native_mstep <- function(context, parameters, eta, map,
                                  maxit, tolerance, trace,
                                  optimizer_state = NULL) {
  started <- proc.time()[["elapsed"]]
  result <- .liberation_saem_mstep(
    lapply(context$subjects, function(evaluator) {
      evaluator$objective_tape$pointer
    }),
    eta, parameters$theta, parameters$sigma, parameters$omega,
    as.integer(map$theta_free), as.integer(map$sigma_free),
    as.numeric(map$lower), as.numeric(map$upper),
    .nm_cpp_prior_config(context$model), as.integer(maxit),
    as.numeric(tolerance), as.integer(trace), optimizer_state
  )
  result$elapsed_seconds <- unname(proc.time()[["elapsed"]] - started)
  result$objective_backend <- "native-cpp-fixed-eta-population-objective"
  result
}

.nm_native_weighted_mstep <- function(context, parameters, native_context, map,
                                      maxit, tolerance, trace,
                                      optimizer_state = NULL) {
  started <- proc.time()[["elapsed"]]
  result <- .liberation_saem_weighted_mstep(
    native_context, parameters$theta, parameters$sigma, parameters$omega,
    as.integer(map$theta_free), as.integer(map$sigma_free),
    as.numeric(map$lower), as.numeric(map$upper),
    .nm_cpp_prior_config(context$model), as.integer(maxit),
    as.numeric(tolerance), as.integer(trace), optimizer_state
  )
  result$elapsed_seconds <- unname(proc.time()[["elapsed"]] - started)
  result$objective_backend <- "native-cpp-weighted-eta-population-objective"
  result
}

.nm_saem_native_weighted_mstep <- function(context, parameters, q_state, map,
                                            maxit, tolerance, trace,
                                            optimizer_state = NULL) {
  .nm_native_weighted_mstep(
    context, parameters, q_state$native(), map, maxit, tolerance, trace,
    optimizer_state
  )
}

.nm_saem_omega_sufficient <- function(context, eta) {
  iov <- context$model$LIK_CONFIG$iov
  if (isTRUE(getOption("LibeRation.saem_native_sufficient_statistics", TRUE))) {
    return(.liberation_saem_omega_sufficient(
      eta, as.integer(context$model$n_eta), as.integer(iov),
      as.integer(context$model$OMEGAS$ROW),
      as.integer(context$model$OMEGAS$COL)
    ))
  }
  covariance <- matrix(0, context$model$n_eta, context$model$n_eta)
  if (iov == 0L) {
    covariance <- crossprod(eta) / max(nrow(eta), 1L)
  } else {
    between <- context$model$n_eta - iov
    if (between) {
      covariance[seq_len(between), seq_len(between)] <-
        crossprod(eta[, seq_len(between), drop = FALSE]) / max(nrow(eta), 1L)
    }
    occasions <- (ncol(eta) - between) / iov
    occasion_effects <- do.call(rbind, lapply(seq_len(occasions), function(occasion) {
      index <- between + (occasion - 1L) * iov + seq_len(iov)
      eta[, index, drop = FALSE]
    }))
    source <- between + seq_len(iov)
    covariance[source, source] <- crossprod(occasion_effects) / max(nrow(occasion_effects), 1L)
  }
  covariance <- covariance + diag(1e-8, nrow(covariance))
  vapply(seq_len(nrow(context$model$OMEGAS)), function(i) {
    covariance[context$model$OMEGAS$ROW[[i]], context$model$OMEGAS$COL[[i]]]
  }, numeric(1))
}

.nm_saem_sigma_sufficient <- function(context, parameters, eta) {
  error <- context$model$LIK_CONFIG$error
  if (!error %in% c("additive", "proportional", "exponential")) return(NULL)
  if (isTRUE(getOption("LibeRation.saem_native_sufficient_statistics", TRUE))) {
    return(.liberation_saem_sigma_sufficient(
      context$engine$pointer, context$data, parameters$theta, eta,
      parameters$sigma
    ))
  }
  prediction <- context$engine$simulate(
    context$data, theta = parameters$theta, eta = eta,
    sigma = parameters$sigma
  )$IPRED
  observed <- context$data$EVID == 0L & context$data$MDV == 0L &
    is.finite(context$data$DV) & is.finite(prediction)
  if (!any(observed)) return(NULL)
  dvid <- if ("DVID" %in% names(context$data)) {
    pmax(as.integer(context$data$DVID), 1L)
  } else rep(1L, nrow(context$data))
  values <- parameters$sigma
  for (response in unique(dvid[observed])) {
    rows <- observed & dvid == response
    residual <- switch(
      error,
      additive = context$data$DV[rows] - prediction[rows],
      proportional = (context$data$DV[rows] - prediction[rows]) /
        pmax(abs(prediction[rows]), 1e-12),
      exponential = {
        valid <- context$data$DV[rows] > 0 & prediction[rows] > 0
        log(context$data$DV[rows][valid]) - log(prediction[rows][valid])
      }
    )
    variance <- mean(residual^2, na.rm = TRUE)
    if (is.finite(variance) && variance > 0 && response <= length(values)) {
      values[[response]] <- if (
        identical(context$model$LIK_CONFIG$sigma_parameterization, "variance")
      ) variance else sqrt(variance)
    }
  }
  values
}

.nm_saem_metropolis_chunk <- function(evaluators, parameters, eta,
                                      proposal_roots, normals, log_uniforms,
                                      mcmc_steps, step_scale,
                                      current_values = NULL) {
  if (!length(evaluators) || !ncol(eta)) {
    return(list(
      eta = eta, value = numeric(nrow(eta)), accepted = 0L, attempted = 0L,
      current_evaluations = 0L, current_cache_hits = 0L,
      candidate_evaluations = 0L
    ))
  }
  initial_eta <- eta
  ode_guard <- isTRUE(evaluators[[1L]]$engine$model$USE_ODE)
  if (ode_guard) invisible(Map(function(evaluator, subject) {
    evaluator$ensure_valid_tapes(
      parameters$theta, parameters$sigma, parameters$omega, eta[subject, ]
    )
  }, evaluators, seq_along(evaluators)))
  make_points <- function(eta) cbind(
    matrix(parameters$theta, nrow(eta), length(parameters$theta), byrow = TRUE),
    eta,
    matrix(parameters$sigma, nrow(eta), length(parameters$sigma), byrow = TRUE),
    matrix(parameters$omega, nrow(eta), length(parameters$omega), byrow = TRUE)
  )
  run <- function(eta, values = current_values) .liberation_objective_tape_eta_metropolis(
    lapply(evaluators, function(evaluator) evaluator$objective_tape$pointer),
    make_points(eta), length(parameters$theta) + seq_len(ncol(eta)), eta,
    proposal_roots, normals, log_uniforms, as.integer(mcmc_steps),
    as.numeric(step_scale), values
  )
  result <- run(initial_eta)
  retaped <- ode_guard && any(vapply(seq_along(evaluators), function(subject) {
    evaluators[[subject]]$ensure_valid_tapes(
      parameters$theta, parameters$sigma, parameters$omega,
      result$eta[subject, ]
    )
  }, logical(1)))
  if (retaped) result <- run(initial_eta, NULL)
  result
}

.nm_saem_independence_chunk <- function(
    evaluators, parameters, eta, proposal_modes, proposal_roots,
    proposal_precisions, normals, log_uniforms, mcmc_steps,
    current_values = NULL, proposal_df = Inf, proposal_scales = NULL) {
  if (!length(evaluators) || !ncol(eta)) {
    return(list(
      eta = eta, value = numeric(nrow(eta)), accepted = 0L, attempted = 0L,
      current_evaluations = 0L, current_cache_hits = 0L,
      candidate_evaluations = 0L, kernel = "laplace-independence"
    ))
  }
  if (nrow(proposal_modes) != length(evaluators) ||
      length(proposal_roots) != length(evaluators) ||
      length(proposal_precisions) != length(evaluators)) {
    .nm_stop("Laplace proposal chunks must contain one entry per subject.")
  }
  student_t <- is.finite(proposal_df)
  if (student_t && (length(proposal_df) != 1L || proposal_df <= 2)) {
    .nm_stop("Student-t Laplace proposals require `proposal_df > 2`.")
  }
  n_draws <- length(evaluators) * mcmc_steps
  proposal_scales <- proposal_scales %||% rep.int(1, n_draws)
  if (length(proposal_scales) != n_draws ||
      any(!is.finite(proposal_scales)) || any(proposal_scales <= 0)) {
    .nm_stop("Laplace proposal scales must contain one positive value per draw.")
  }
  ode_guard <- isTRUE(evaluators[[1L]]$engine$model$USE_ODE)
  evaluate <- function(evaluator, subject_eta) {
    if (ode_guard) evaluator$ensure_valid_tapes(
      parameters$theta, parameters$sigma, parameters$omega, subject_eta
    )
    evaluator$objective(
      parameters$theta, subject_eta, parameters$sigma, parameters$omega,
      gradient = FALSE
    )$value
  }
  values <- if (is.null(current_values)) {
    vapply(seq_along(evaluators), function(subject) {
      evaluate(evaluators[[subject]], eta[subject, ])
    }, numeric(1))
  } else as.numeric(current_values)
  if (length(values) != length(evaluators) || any(!is.finite(values))) {
    .nm_stop("Current f-SAEM subject objectives must be finite.")
  }
  current_evaluations <- if (is.null(current_values)) length(evaluators) else 0L
  current_cache_hits <- if (is.null(current_values)) 0L else length(evaluators)
  accepted <- candidate_evaluations <- 0L
  quadratic <- function(value, center, precision) {
    centered <- value - center
    drop(crossprod(centered, precision %*% centered))
  }
  for (subject in seq_along(evaluators)) {
    mode <- proposal_modes[subject, ]
    root <- proposal_roots[[subject]]
    precision <- proposal_precisions[[subject]]
    current_quad <- quadratic(eta[subject, ], mode, precision)
    for (step in seq_len(mcmc_steps)) {
      draw <- (subject - 1L) * mcmc_steps + step
      candidate <- as.vector(
        mode + proposal_scales[[draw]] * root %*% normals[draw, ]
      )
      candidate_value <- tryCatch(
        evaluate(evaluators[[subject]], candidate), error = function(error) Inf
      )
      candidate_evaluations <- candidate_evaluations + 1L
      candidate_quad <- quadratic(candidate, mode, precision)
      log_proposal_ratio <- if (student_t) {
        0.5 * (proposal_df + ncol(eta)) * (
          log1p(candidate_quad / proposal_df) -
            log1p(current_quad / proposal_df)
        )
      } else 0.5 * (candidate_quad - current_quad)
      log_ratio <- -0.5 * (candidate_value - values[[subject]]) +
        log_proposal_ratio
      if (is.finite(candidate_value) && log_uniforms[[draw]] < log_ratio) {
        eta[subject, ] <- candidate
        values[[subject]] <- candidate_value
        current_quad <- candidate_quad
        accepted <- accepted + 1L
      }
    }
  }
  list(
    eta = eta, value = values, accepted = accepted,
    attempted = length(evaluators) * mcmc_steps,
    current_evaluations = current_evaluations,
    current_cache_hits = current_cache_hits,
    candidate_evaluations = candidate_evaluations,
    kernel = if (student_t) "laplace-student-t-independence" else
      "laplace-independence"
  )
}

.nm_proposal_root_groups <- function(context) {
  if (is.null(context$model$RE_CONFIG)) {
    if (context$model$LIK_CONFIG$iov == 0L) {
      return(rep.int(1L, context$n_subjects))
    }
    return(match(
      vapply(context$subjects, `[[`, integer(1), "n_eta"),
      unique(vapply(context$subjects, `[[`, integer(1), "n_eta"))
    ))
  }
  total_names <- paste0(
    ".RE_TOTAL_", seq_along(context$model$RE_CONFIG$blocks)
  )
  keys <- vapply(context$subjects, function(evaluator) {
    totals <- evaluator$project(total_names, first_only = TRUE)
    paste(
      vapply(total_names, function(name) {
        as.integer(totals[[name]][[1L]])
      }, integer(1)),
      collapse = ":"
    )
  }, character(1))
  match(keys, unique(keys))
}

.nm_proposal_root_cache <- function(context) {
  cache <- new.env(parent = emptyenv())
  cache$omega <- NULL
  cache$roots <- NULL
  cache$groups <- .nm_proposal_root_groups(context)
  cache$hits <- 0L
  cache$misses <- 0L
  cache$factorizations <- 0L
  cache
}

.nm_proposal_roots <- function(context, omega, cache = NULL) {
  # Reusing an exactly identical factorisation is arithmetic-neutral and does
  # not alter proposal ordering or random-number consumption, so it is safe in
  # both numerical policies.
  if (!is.null(cache) && !is.null(cache$omega) &&
      identical(cache$omega, omega)) {
    cache$hits <- cache$hits + 1L
    return(cache$roots)
  }
  groups <- cache$groups %||% .nm_proposal_root_groups(context)
  representatives <- match(unique(groups), groups)
  unique_roots <- lapply(representatives, function(subject) {
    covariance <- .nm_effect_covariance_evaluator(
      context$model, context$subjects[[subject]], omega
    )
    t(chol(covariance))
  })
  roots <- lapply(groups, function(group) unique_roots[[group]])
  if (!is.null(cache)) {
    cache$omega <- omega
    cache$roots <- roots
    cache$misses <- cache$misses + 1L
    cache$factorizations <- cache$factorizations + length(unique_roots)
  }
  roots
}

.nm_saem_metropolis <- function(context, parameters, eta, mcmc_steps,
                                 step_scale, proposal_roots = NULL,
                                 proposal_cache = NULL,
                                 current_values = NULL,
                                 persistent = NULL,
                                 independence = NULL) {
  if (!context$n_eta) {
    return(list(
      eta = eta, value = rep(0, context$n_subjects),
      accepted = 0L, attempted = 0L
    ))
  }
  roots <- if (!is.null(independence)) independence$roots else
    proposal_roots %||% .nm_proposal_roots(
      context, parameters$omega, cache = proposal_cache
    )
  if (length(roots) != context$n_subjects) {
    .nm_stop("`proposal_roots` must contain one covariance root per subject.")
  }
  normals <- matrix(
    stats::rnorm(context$n_subjects * mcmc_steps * context$n_eta),
    context$n_subjects * mcmc_steps, context$n_eta
  )
  log_uniforms <- log(stats::runif(context$n_subjects * mcmc_steps))
  proposal_df <- if (is.null(independence)) Inf else
    as.numeric(independence$df %||% Inf)
  student_t <- !is.null(independence) && is.finite(proposal_df)
  proposal_scales <- if (student_t) {
    sqrt(proposal_df / stats::rchisq(
      context$n_subjects * mcmc_steps, df = proposal_df
    ))
  } else rep.int(1, context$n_subjects * mcmc_steps)
  if (is.null(context$parallel)) {
    if (!is.null(persistent)) {
      if (!is.null(independence)) {
        return(.liberation_stochastic_eta_context_independence(
          persistent, parameters$theta, eta, parameters$sigma,
          parameters$omega, independence$modes, independence$roots,
          independence$precisions, normals, log_uniforms,
          as.integer(mcmc_steps), current_values, proposal_df,
          proposal_scales
        ))
      }
      return(.liberation_stochastic_eta_context_random_walk(
        persistent, parameters$theta, eta, parameters$sigma,
        parameters$omega, roots, normals, log_uniforms,
        as.integer(mcmc_steps), as.numeric(step_scale), current_values
      ))
    }
    if (!is.null(independence)) {
      return(.nm_saem_independence_chunk(
        context$subjects, parameters, eta, independence$modes,
        independence$roots, independence$precisions, normals, log_uniforms,
        mcmc_steps, current_values, proposal_df, proposal_scales
      ))
    }
    return(.nm_saem_metropolis_chunk(
      context$subjects, parameters, eta, roots, normals, log_uniforms,
      mcmc_steps, step_scale, current_values
    ))
  }
  chunks <- context$parallel$chunks
  eta_chunks <- lapply(chunks, function(rows) eta[rows, , drop = FALSE])
  root_chunks <- lapply(chunks, function(rows) roots[rows])
  normal_chunks <- lapply(chunks, function(rows) {
    draws <- unlist(lapply(rows, function(subject) {
      (subject - 1L) * mcmc_steps + seq_len(mcmc_steps)
    }))
    normals[draws, , drop = FALSE]
  })
  uniform_chunks <- lapply(chunks, function(rows) {
    draws <- unlist(lapply(rows, function(subject) {
      (subject - 1L) * mcmc_steps + seq_len(mcmc_steps)
    }))
    log_uniforms[draws]
  })
  scale_chunks <- lapply(chunks, function(rows) {
    draws <- unlist(lapply(rows, function(subject) {
      (subject - 1L) * mcmc_steps + seq_len(mcmc_steps)
    }))
    proposal_scales[draws]
  })
  value_chunks <- if (is.null(current_values)) {
    rep(list(NULL), length(chunks))
  } else lapply(chunks, function(rows) current_values[rows])
  mode_chunks <- if (is.null(independence)) {
    rep(list(NULL), length(chunks))
  } else lapply(chunks, function(rows) {
    independence$modes[rows, , drop = FALSE]
  })
  precision_chunks <- if (is.null(independence)) {
    rep(list(NULL), length(chunks))
  } else lapply(chunks, function(rows) independence$precisions[rows])
  pieces <- parallel::clusterApply(
    context$parallel$cluster, seq_along(chunks),
    function(index, parameters, eta_chunks, root_chunks, normal_chunks,
             uniform_chunks, scale_chunks, value_chunks, mode_chunks,
             precision_chunks, mcmc_steps, step_scale, independent,
             proposal_df) {
      namespace <- asNamespace("LibeRation")
      evaluators <- get(".nm_parallel_worker_state", envir = namespace)()$subjects
      if (isTRUE(independent)) {
        sampler <- get(".nm_saem_independence_chunk", envir = namespace)
        sampler(
          evaluators, parameters, eta_chunks[[index]], mode_chunks[[index]],
          root_chunks[[index]], precision_chunks[[index]],
          normal_chunks[[index]], uniform_chunks[[index]], mcmc_steps,
          value_chunks[[index]], proposal_df, scale_chunks[[index]]
        )
      } else {
        sampler <- get(".nm_saem_metropolis_chunk", envir = namespace)
        sampler(
          evaluators, parameters, eta_chunks[[index]], root_chunks[[index]],
          normal_chunks[[index]], uniform_chunks[[index]], mcmc_steps,
          step_scale, value_chunks[[index]]
        )
      }
    }, parameters = parameters, eta_chunks = eta_chunks,
    root_chunks = root_chunks, normal_chunks = normal_chunks,
    uniform_chunks = uniform_chunks, scale_chunks = scale_chunks,
    value_chunks = value_chunks,
    mode_chunks = mode_chunks, precision_chunks = precision_chunks,
    mcmc_steps = mcmc_steps, step_scale = step_scale,
    independent = !is.null(independence), proposal_df = proposal_df
  )
  list(
    eta = do.call(rbind, lapply(pieces, `[[`, "eta")),
    value = unlist(lapply(pieces, `[[`, "value"), use.names = FALSE),
    accepted = sum(vapply(pieces, `[[`, integer(1), "accepted")),
    attempted = sum(vapply(pieces, `[[`, integer(1), "attempted")),
    current_evaluations = sum(vapply(
      pieces, `[[`, integer(1), "current_evaluations"
    )),
    current_cache_hits = sum(vapply(
      pieces, `[[`, integer(1), "current_cache_hits"
    )),
    candidate_evaluations = sum(vapply(
      pieces, `[[`, integer(1), "candidate_evaluations"
    ))
  )
}

.nm_saem_stationarity <- function(objective, parameters, burn,
                                  window = 20L, tolerance = 1e-3) {
  completed <- length(objective)
  first_post_burn <- as.integer(burn) + 1L
  available <- completed - first_post_burn + 1L
  required <- max(4L, min(as.integer(window), 10L))
  if (available < required || !ncol(parameters)) {
    return(list(
      ready = FALSE, converged = FALSE, window = max(available, 0L),
      parameter_drift = NA_real_, objective_drift = NA_real_,
      tolerance = tolerance
    ))
  }
  width <- min(as.integer(window), available)
  rows <- seq.int(completed - width + 1L, completed)
  parameter_window <- parameters[rows, , drop = FALSE]
  parameter_scale <- pmax(abs(colMeans(parameter_window)), 1)
  parameter_drift <- max(abs(
    parameter_window[nrow(parameter_window), ] - parameter_window[1L, ]
  ) / parameter_scale)
  objective_window <- objective[rows]
  objective_scale <- max(abs(mean(objective_window)), 1)
  index <- seq_along(objective_window)
  centered_index <- index - mean(index)
  slope <- sum(centered_index * (objective_window - mean(objective_window))) /
    sum(centered_index^2)
  objective_drift <- abs(slope) * max(width - 1L, 1L) / objective_scale
  list(
    ready = TRUE,
    converged = is.finite(parameter_drift) && is.finite(objective_drift) &&
      parameter_drift <= tolerance && objective_drift <= tolerance,
    window = width, parameter_drift = parameter_drift,
    objective_drift = objective_drift, tolerance = tolerance
  )
}

.nm_est_saem_single <- function(context, map, maxit, tolerance, trace,
                         n_iter = 200L, burn = NULL, mcmc_steps = 2L,
                         step_scale = 0.5, sa_power = 0.7,
                         mstep_maxit = 20L, seed = 20260713L,
                         print_every = 0L, adapt_proposal = TRUE,
                         target_acceptance = 0.3, closed_form_sigma = TRUE,
                         optimizer_backend = "auto", initial_eta = NULL,
                         mu_specialization = TRUE,
                         saem_kernel = c("auto", "random_walk", "fsaem"),
                         fsaem_distribution = c("auto", "gaussian", "student_t"),
                         fsaem_df = 7,
                         fsaem_refresh = 25L,
                         fsaem_eta_maxit = 50L,
                         fsaem_rescue_probability = 0.1,
                         fsaem_parameter_refresh = 0.15,
                         fsaem_low_acceptance = 0.1,
                         stationarity_window = 20L,
                         stationarity_tolerance = 1e-3,
                         auto_stop = NULL,
                         auto_stop_consecutive = 3L,
                         auto_stop_min_iterations = NULL,
                         saem_support_max = 0L,
                         saem_support_prune = 0,
                         saem_mstep_interval_burn = NULL,
                         saem_mstep_interval = NULL,
                         saem_parameter_averaging = c("auto", "none", "polyak"),
                         saem_average_start = NULL) {
  n_iter <- as.integer(n_iter)
  burn <- as.integer(burn %||% floor(n_iter / 3))
  mcmc_steps <- as.integer(mcmc_steps)
  saem_kernel <- match.arg(saem_kernel)
  fsaem_distribution <- match.arg(fsaem_distribution)
  saem_parameter_averaging <- match.arg(saem_parameter_averaging)
  fsaem_df <- as.numeric(fsaem_df)
  fsaem_refresh <- as.integer(fsaem_refresh)
  fsaem_eta_maxit <- as.integer(fsaem_eta_maxit)
  stationarity_window <- as.integer(stationarity_window)
  auto_stop_consecutive <- as.integer(auto_stop_consecutive)
  auto_stop_min_iterations <- as.integer(
    auto_stop_min_iterations %||% max(burn + stationarity_window, burn + 10L)
  )
  auto_stop <- isTRUE(auto_stop %||% .nm_liber_optimized(context))
  if (!.nm_liber_optimized(context)) auto_stop <- FALSE
  saem_support_max <- as.integer(saem_support_max)
  saem_support_prune <- as.numeric(saem_support_prune)
  saem_mstep_interval_burn <- as.integer(
    saem_mstep_interval_burn %||% if (.nm_liber_optimized(context)) 4L else 1L
  )
  saem_mstep_interval <- as.integer(
    saem_mstep_interval %||% if (.nm_liber_optimized(context)) 2L else 1L
  )
  if (!.nm_liber_optimized(context)) {
    saem_mstep_interval_burn <- 1L
    saem_mstep_interval <- 1L
  }
  resolved_parameter_averaging <- if (saem_parameter_averaging == "auto") {
    if (.nm_liber_optimized(context)) "polyak" else "none"
  } else saem_parameter_averaging
  if (!.nm_liber_optimized(context)) resolved_parameter_averaging <- "none"
  saem_average_start <- as.integer(saem_average_start %||% (burn + 1L))
  if (n_iter < 2L || burn < 0L || burn >= n_iter || mcmc_steps < 1L) {
    .nm_stop("SAEM requires n_iter >= 2, 0 <= burn < n_iter, and mcmc_steps >= 1.")
  }
  if (!is.finite(step_scale) || step_scale <= 0 ||
      !is.finite(target_acceptance) || target_acceptance <= 0 ||
      target_acceptance >= 1) {
    .nm_stop("SAEM proposal scale must be positive and target acceptance must lie in (0, 1).")
  }
  if (is.na(fsaem_refresh) || fsaem_refresh < 1L ||
      is.na(fsaem_eta_maxit) || fsaem_eta_maxit < 1L) {
    .nm_stop(
      "f-SAEM refresh and ETA-mode iteration controls must be positive integers."
    )
  }
  if (length(fsaem_df) != 1L || !is.finite(fsaem_df) || fsaem_df <= 2 ||
      !is.finite(fsaem_rescue_probability) ||
      fsaem_rescue_probability < 0 || fsaem_rescue_probability >= 1 ||
      !is.finite(fsaem_parameter_refresh) || fsaem_parameter_refresh <= 0 ||
      !is.finite(fsaem_low_acceptance) || fsaem_low_acceptance < 0 ||
      fsaem_low_acceptance >= 1) {
    .nm_stop(
      "f-SAEM rescue probability and low-acceptance threshold must lie in ",
      "[0, 1), and the parameter-refresh threshold must be positive."
    )
  }
  if (is.na(stationarity_window) || stationarity_window < 4L ||
      !is.finite(stationarity_tolerance) || stationarity_tolerance <= 0 ||
      is.na(auto_stop_consecutive) || auto_stop_consecutive < 1L ||
      is.na(auto_stop_min_iterations) || auto_stop_min_iterations < 2L ||
      is.na(saem_support_max) || saem_support_max < 0L ||
      length(saem_support_prune) != 1L || !is.finite(saem_support_prune) ||
      saem_support_prune < 0 || saem_support_prune >= 1 ||
      is.na(saem_mstep_interval_burn) || saem_mstep_interval_burn < 1L ||
      is.na(saem_mstep_interval) || saem_mstep_interval < 1L ||
      is.na(saem_average_start) || saem_average_start < 1L ||
      saem_average_start > n_iter) {
    .nm_stop(
      "SAEM stationarity requires a window >= 4, a positive tolerance, ",
      "positive consecutive count, and minimum iterations >= 2."
    )
  }
  fsaem_eligible <- .nm_liber_optimized(context) && context$n_eta > 0L
  resolved_kernel <- if (saem_kernel == "auto") {
    if (fsaem_eligible) "fsaem" else "random_walk"
  } else saem_kernel
  resolved_fsaem_distribution <- if (fsaem_distribution == "auto") {
    if (resolved_kernel == "fsaem") "student_t" else "gaussian"
  } else fsaem_distribution
  if (resolved_kernel == "fsaem" && !fsaem_eligible) {
    .nm_stop(
      "The f-SAEM Laplace-independence kernel requires liber_optimized ",
      "and at least one random effect."
    )
  }
  set.seed(seed)
  parameters <- map$decode(map$start)
  eta <- initial_eta %||% matrix(0, context$n_subjects, context$n_eta)
  accepted <- attempted <- 0L
  objective_trace <- numeric(n_iter)
  acceptance_trace <- numeric(n_iter)
  step_scale_trace <- numeric(n_iter)
  parameter_trace <- matrix(NA_real_, n_iter, length(map$start))
  colnames(parameter_trace) <- map$names
  completed_iterations <- 0L
  stationary_iterations <- 0L
  stationarity <- .nm_saem_stationarity(
    numeric(), parameter_trace[FALSE, , drop = FALSE], burn,
    stationarity_window, stationarity_tolerance
  )
  mstep_objective_evaluations <- 0L
  mstep_gradient_evaluations <- 0L
  mstep_iterations <- 0L
  mstep_elapsed <- 0
  mstep_backend <- "unknown"
  mu <- .nm_mu_specialization(context, map, enabled = mu_specialization)
  mu_updates <- 0L
  mu_fallbacks <- 0L
  mu_closed_form_only_iterations <- 0L
  proposal_cache <- .nm_proposal_root_cache(context)
  stochastic_context <- .nm_stochastic_eta_context(
    context, allow_compatibility = TRUE
  )
  phases <- .nm_stochastic_phase_timer()
  fsaem_proposal <- NULL
  fsaem_refreshes <- 0L
  fsaem_refresh_failures <- 0L
  fsaem_fallback_iterations <- 0L
  fsaem_mode_iterations <- 0L
  fsaem_mode_evaluations <- 0L
  fsaem_last_error <- NULL
  fsaem_force_refresh <- FALSE
  fsaem_rescue_iterations <- 0L
  fsaem_acceptance_refreshes <- 0L
  fsaem_parameter_refreshes <- 0L
  current_subject_values <- NULL
  eta_current_evaluations <- 0L
  eta_current_cache_hits <- 0L
  eta_candidate_evaluations <- 0L
  saem_priors <- context$model$LIK_CONFIG$priors
  saem_sigma_prior <- !is.null(saem_priors) && nrow(saem_priors) &&
    any(startsWith(toupper(saem_priors$parameter), "SIGMA"))
  simple_sigma <- isTRUE(closed_form_sigma) && !saem_sigma_prior &&
    context$model$LIK_CONFIG$error %in%
      c("additive", "proportional", "exponential") &&
    identical(context$model$LIK_CONFIG$sigma_corr %||% "independent",
              "independent") &&
    !length(context$model$LIK_CONFIG$residual_groups)
  native_mstep_model <- context$model
  native_mstep_model$THETAS$Value <- parameters$theta
  native_mstep_model$SIGMAS$Value <- parameters$sigma
  native_mstep_model$OMEGAS$Value <- parameters$omega
  if (length(map$omega_free)) native_mstep_model$OMEGAS$FIX[] <- TRUE
  if (simple_sigma && length(map$sigma_free)) {
    native_mstep_model$SIGMAS$FIX[] <- TRUE
  }
  native_mstep_map <- .nm_outer_map(native_mstep_model)
  native_mstep_enabled <- isTRUE(getOption(
    "LibeRation.saem_native_mstep", TRUE
  )) && .nm_liber_optimized(context) && is.null(context$parallel) &&
    optimizer_backend %in% c("auto", "native") &&
    !isTRUE(mu$active)
  native_mstep_attempts <- 0L
  native_mstep_successes <- 0L
  native_mstep_fallbacks <- 0L
  native_mstep_fallback_reason <- NULL
  native_mstep_state <- NULL
  sigma_gradient_updates <- 0L
  sigma_expectation_fallbacks <- 0L
  paired_objective_evaluations <- 0L
  paired_objective_cache_hits <- 0L
  paired_objective_iterations <- 0L
  persistent_objective_iterations <- 0L
  persistent_objective_fallbacks <- 0L
  persistent_objective_fallback_reason <- NULL
  q_state <- .nm_saem_q_state(
    context,
    max_support = if (.nm_liber_optimized(context)) saem_support_max else 0L,
    prune_tolerance = if (.nm_liber_optimized(context)) {
      saem_support_prune
    } else 0
  )
  mstep_performed_trace <- logical(n_iter)
  parameter_average <- numeric(length(map$start))
  parameter_average_count <- 0L
  update_parameter_average <- function(iteration, parameters) {
    if (resolved_parameter_averaging != "polyak" ||
        iteration < saem_average_start || !length(parameter_average)) {
      return(invisible(NULL))
    }
    parameter_average_count <<- parameter_average_count + 1L
    point <- map$encode(parameters)
    parameter_average <<- parameter_average +
      (point - parameter_average) / parameter_average_count
    invisible(NULL)
  }
  for (iteration in seq_len(n_iter)) {
    previous_parameters <- parameters
    if (context$n_eta) {
      independence <- NULL
      parameter_refresh <- FALSE
      if (resolved_kernel == "fsaem" && !is.null(fsaem_proposal)) {
        current_anchor <- c(
          parameters$theta, parameters$sigma, parameters$omega
        )
        parameter_refresh <- max(
          abs(current_anchor - fsaem_proposal$anchor) /
            (1 + abs(fsaem_proposal$anchor))
        ) > fsaem_parameter_refresh
      }
      if (resolved_kernel == "fsaem" &&
          (is.null(fsaem_proposal) || fsaem_force_refresh ||
           parameter_refresh ||
           (iteration - 1L) %% fsaem_refresh == 0L)) {
        if (parameter_refresh) {
          fsaem_parameter_refreshes <- fsaem_parameter_refreshes + 1L
        }
        refreshed <- tryCatch(
          .nm_fsaem_proposal(
            context, parameters,
            if (is.null(fsaem_proposal)) eta else fsaem_proposal$modes,
            fsaem_eta_maxit, tolerance, persistent = stochastic_context
          ),
          error = identity
        )
        if (inherits(refreshed, "error")) {
          fsaem_refresh_failures <- fsaem_refresh_failures + 1L
          fsaem_last_error <- conditionMessage(refreshed)
        } else {
          fsaem_proposal <- refreshed
          fsaem_refreshes <- fsaem_refreshes + 1L
          fsaem_mode_iterations <- fsaem_mode_iterations +
            refreshed$mode_iterations
          fsaem_mode_evaluations <- fsaem_mode_evaluations +
            refreshed$mode_evaluations
          fsaem_force_refresh <- FALSE
        }
      }
      if (resolved_kernel == "fsaem") {
        independence <- fsaem_proposal
        if (!is.null(independence)) {
          independence$df <- if (resolved_fsaem_distribution == "student_t") {
            fsaem_df
          } else Inf
        }
        if (is.null(independence)) {
          fsaem_fallback_iterations <- fsaem_fallback_iterations + 1L
        } else if (fsaem_rescue_probability > 0 &&
                   stats::runif(1) < fsaem_rescue_probability) {
          # A random-walk rescue kernel is individually invariant for the
          # conditional target. Randomly mixing it with the Laplace
          # independence kernel therefore preserves the exact target while
          # improving tail and secondary-mode exploration.
          independence <- NULL
          fsaem_rescue_iterations <- fsaem_rescue_iterations + 1L
        }
      }
      sampled <- phases$time("eta_sampling", {
        .nm_saem_metropolis(
          context, parameters, eta, mcmc_steps, step_scale,
          proposal_cache = proposal_cache,
          current_values = current_subject_values,
          persistent = stochastic_context,
          independence = independence
        )
      })
      eta <- sampled$eta
      current_subject_values <- NULL
      accepted <- accepted + sampled$accepted
      attempted <- attempted + sampled$attempted
      eta_current_evaluations <- eta_current_evaluations +
        as.integer(sampled$current_evaluations %||% 0L)
      eta_current_cache_hits <- eta_current_cache_hits +
        as.integer(sampled$current_cache_hits %||% 0L)
      eta_candidate_evaluations <- eta_candidate_evaluations +
        as.integer(sampled$candidate_evaluations %||% 0L)
      acceptance_trace[[iteration]] <- sampled$accepted / max(sampled$attempted, 1L)
      if (resolved_kernel == "fsaem" && !is.null(independence) &&
          acceptance_trace[[iteration]] < fsaem_low_acceptance) {
        fsaem_force_refresh <- TRUE
        fsaem_acceptance_refreshes <- fsaem_acceptance_refreshes + 1L
      }
      if (isTRUE(adapt_proposal) && iteration <= burn &&
          is.null(independence)) {
        gain <- min(0.1, 1 / sqrt(iteration))
        step_scale <- step_scale * exp(
          gain * (acceptance_trace[[iteration]] - target_acceptance)
        )
      }
    }
    step_scale_trace[[iteration]] <- step_scale
    gamma <- if (iteration <= burn) 1 else (iteration - burn)^(-sa_power)
    phases$time("stochastic_approximation", q_state$update(eta, gamma))
    mstep_interval_current <- if (iteration <= burn) {
      saem_mstep_interval_burn
    } else saem_mstep_interval
    mstep_due <- !.nm_liber_optimized(context) || iteration == 1L ||
      iteration == n_iter || iteration %% mstep_interval_current == 0L
    mstep_performed_trace[[iteration]] <- mstep_due
    if (!mstep_due) {
      objective_state <- .nm_saem_conditional_components(
        context, parameters, eta, persistent = stochastic_context
      )
      objective_trace[[iteration]] <- objective_state$value
      if (length(map$start)) {
        parameter_trace[iteration, ] <- map$encode(parameters)
      }
      completed_iterations <- iteration
      current_subject_values <- objective_state$subject
      update_parameter_average(iteration, parameters)
      next
    }
    mstep_model <- context$model
    mstep_model$THETAS$Value <- parameters$theta
    mstep_model$SIGMAS$Value <- parameters$sigma
    mstep_model$OMEGAS$Value <- parameters$omega
    if (length(map$omega_free)) mstep_model$OMEGAS$FIX[] <- TRUE
    if (simple_sigma && length(map$sigma_free)) mstep_model$SIGMAS$FIX[] <- TRUE
    eta_mstep <- q_state$mean_eta()
    mu_iteration_active <- FALSE
    if (isTRUE(mu$saem_eligible) && isTRUE(mu$active)) {
      mu_update <- .nm_mu_gls_update(mu, context, parameters, eta_mstep)
      if (isTRUE(mu_update$valid)) {
        mstep_model$THETAS$Value <- mu_update$parameters$theta
        mstep_model$THETAS$FIX[mu$theta] <- TRUE
        eta_mstep <- mu_update$eta
        mu_updates <- mu_updates + 1L
        mu_iteration_active <- TRUE
      } else {
        mu_fallbacks <- mu_fallbacks + 1L
        mu$runtime_reason <- mu_update$reason %||% "MU GLS update unavailable"
      }
    }
    paired_conditional <- NULL
    use_native_mstep <- native_mstep_enabled && !mu_iteration_active &&
      !is.null(q_state$native())
    native_result <- NULL
    if (use_native_mstep && length(native_mstep_map$start)) {
      native_mstep_attempts <- native_mstep_attempts + 1L
      native_result <- phases$time("mstep", tryCatch(
          .nm_saem_native_weighted_mstep(
            context, parameters, q_state, native_mstep_map,
            min(as.integer(mstep_maxit), as.integer(maxit)), tolerance,
            if (trace > 1L) trace else 0L, native_mstep_state
          ),
          error = identity
        ))
      if (inherits(native_result, "error") ||
          !is.finite(native_result$value %||% NA_real_) ||
          any(!is.finite(native_result$theta %||% NA_real_)) ||
          any(!is.finite(native_result$sigma %||% NA_real_))) {
        native_mstep_fallbacks <- native_mstep_fallbacks + 1L
        native_mstep_fallback_reason <- if (inherits(native_result, "error")) {
          conditionMessage(native_result)
        } else {
          "native SAEM M-step returned non-finite values"
        }
        native_result <- NULL
        native_mstep_enabled <- FALSE
      } else {
        native_mstep_successes <- native_mstep_successes + 1L
        native_mstep_state <- native_result$optimizer_state
        native_result$optimizer_state <- NULL
      }
    }
    iteration_map <- NULL
    if (!is.null(native_result)) {
      maximized <- native_result
      candidate <- list(
        theta = as.numeric(native_result$theta),
        sigma = as.numeric(native_result$sigma),
        omega = as.numeric(native_result$omega)
      )
    } else {
      iteration_map <- .nm_outer_map(mstep_model)
    }
    if (is.null(native_result) && !length(iteration_map$start)) {
      maximized <- list(
        par = numeric(), value = NA_real_, convergence = 0L,
        message = "All SAEM M-step parameters used closed-form updates",
        counts = c(`function` = 0L, gradient = 0L),
        iterations = 0L, objective_evaluations = 0L,
        gradient_evaluations = 0L, elapsed_seconds = 0,
        backend = "saem-closed-form"
      )
      if (mu_iteration_active) {
        mu_closed_form_only_iterations <- mu_closed_form_only_iterations + 1L
      }
      candidate <- iteration_map$decode(maximized$par)
    } else if (is.null(native_result)) {
      paired_conditional <- if (!is.null(q_state$native())) {
        .nm_saem_q_expectation(context, iteration_map, q_state)
      } else if (q_state$count() == 1L) {
        .nm_saem_paired_conditional(context, iteration_map, eta_mstep)
      } else {
        .nm_saem_q_expectation(context, iteration_map, q_state)
      }
      maximized <- phases$time("mstep", {
        .nm_outer_optim(
          iteration_map, paired_conditional$objective,
          min(as.integer(mstep_maxit), as.integer(maxit)),
          tolerance, if (trace > 1L) trace else 0L,
          gradient = paired_conditional$gradient,
          optimizer_backend = optimizer_backend
        )
      })
      paired_telemetry <- paired_conditional$telemetry()
      paired_objective_evaluations <- paired_objective_evaluations +
        as.integer(paired_telemetry$evaluations)
      paired_objective_cache_hits <- paired_objective_cache_hits +
        as.integer(paired_telemetry$cache_hits)
      paired_objective_iterations <- paired_objective_iterations +
        as.integer(isTRUE(paired_telemetry$paired %||% FALSE))
      persistent_objective_iterations <- persistent_objective_iterations +
        as.integer(isTRUE(paired_telemetry$persistent %||% FALSE))
      if (isTRUE(paired_telemetry$persistent_requested %||% FALSE) &&
          !isTRUE(paired_telemetry$persistent)) {
        persistent_objective_fallbacks <- persistent_objective_fallbacks + 1L
        persistent_objective_fallback_reason <-
          paired_telemetry$persistent_error %||%
          persistent_objective_fallback_reason
      }
      candidate <- iteration_map$decode(maximized$par)
    }
    mstep_objective_evaluations <- mstep_objective_evaluations +
      as.integer(maximized$objective_evaluations %||% 0L)
    mstep_gradient_evaluations <- mstep_gradient_evaluations +
      as.integer(maximized$gradient_evaluations %||% 0L)
    mstep_iterations <- mstep_iterations + as.integer(maximized$iterations %||% 0L)
    mstep_elapsed <- mstep_elapsed + as.numeric(maximized$elapsed_seconds %||% 0)
    mstep_backend <- maximized$backend %||% mstep_backend
    # Q_k already contains the Robbins--Monro stochastic approximation. The
    # M-step must maximize that accumulated auxiliary function; applying gamma
    # a second time to the resulting parameter estimate is not canonical SAEM.
    parameters$theta[map$theta_free] <- candidate$theta[map$theta_free]
    if (mu_iteration_active) {
      eta <- .nm_mu_recenter_eta(
        mu, previous_parameters, parameters, eta
      )
      q_state$recenter(mu, previous_parameters, parameters)
    }
    sigma_parameters <- candidate
    sigma_parameters$theta <- parameters$theta
    sigma_sufficient <- if (simple_sigma) {
      phases$time("sigma_sufficient", {
        # The native weighted M-step has already evaluated the complete Q
        # gradient at its accepted point. Reuse that exact Reverse(1) result
        # rather than replaying every retained ETA support solely to recover
        # the closed-form residual-variance update.
        gradient_sufficient <- if (!is.null(native_result)) {
          .nm_saem_sigma_from_q_gradient(
            context, sigma_parameters, native_result$native_gradient
          )
        } else NULL
        if (!is.null(q_state$native())) {
          if (is.null(gradient_sufficient)) {
            q_expectation <- paired_conditional
            if (is.null(q_expectation)) {
              gradient_map <- iteration_map %||% native_mstep_map
              q_expectation <- .nm_saem_q_expectation(
                context, gradient_map, q_state
              )
            }
            q_gradient <- tryCatch(
              q_expectation$native_gradient(sigma_parameters),
              error = function(error) NULL
            )
            gradient_sufficient <- .nm_saem_sigma_from_q_gradient(
              context, sigma_parameters, q_gradient
            )
          }
        }
        if (!is.null(gradient_sufficient)) {
          sigma_gradient_updates <- sigma_gradient_updates + 1L
          gradient_sufficient
        } else {
          sigma_expectation_fallbacks <- sigma_expectation_fallbacks + 1L
          .nm_saem_sigma_expectation(context, sigma_parameters, q_state)
        }
      })
    } else NULL
    if (!is.null(sigma_sufficient) && length(map$sigma_free)) {
      candidate$sigma[map$sigma_free] <- sigma_sufficient[map$sigma_free]
    }
    parameters$sigma[map$sigma_free] <- candidate$sigma[map$sigma_free]
    if (length(map$omega_free) && context$n_eta) {
      omega_sufficient <- phases$time("omega_sufficient", {
        .nm_saem_omega_expectation(context, q_state)
      })
      parameters$omega[map$omega_free] <- omega_sufficient[map$omega_free]
    }
    objective_state <- .nm_saem_conditional_components(
      context, parameters, eta, persistent = stochastic_context
    )
    objective_trace[[iteration]] <- objective_state$value
    if (length(map$start)) parameter_trace[iteration, ] <- map$encode(parameters)
    completed_iterations <- iteration
    current_subject_values <- objective_state$subject
    update_parameter_average(iteration, parameters)
    if (print_every > 0L && iteration %% print_every == 0L && length(map$start)) {
      point <- map$encode(parameters)
      objective_at <- function(value) .nm_saem_conditional(
        context, map$decode(value), eta, persistent = stochastic_context
      )
      gradient_at <- function(value) {
        .nm_saem_conditional_gradient(context, map, map$decode(value), eta)
      }
      .nm_log_gradient(
        iteration, objective_at, point, map, objective_trace[[iteration]],
        gradient_function = gradient_at
      )
    }
    stationarity <- .nm_saem_stationarity(
      objective_trace[seq_len(iteration)],
      parameter_trace[seq_len(iteration), , drop = FALSE], burn,
      stationarity_window, stationarity_tolerance
    )
    stationary_iterations <- if (isTRUE(stationarity$converged)) {
      stationary_iterations + 1L
    } else 0L
    if (isTRUE(auto_stop) && iteration >= auto_stop_min_iterations &&
        stationary_iterations >= auto_stop_consecutive) break
  }
  used <- seq_len(completed_iterations)
  objective_trace <- objective_trace[used]
  acceptance_trace <- acceptance_trace[used]
  step_scale_trace <- step_scale_trace[used]
  parameter_trace <- parameter_trace[used, , drop = FALSE]
  parameter_average_applied <- resolved_parameter_averaging == "polyak" &&
    parameter_average_count > 0L && length(parameter_average)
  if (parameter_average_applied) {
    parameters <- map$decode(parameter_average)
  }
  final_objective_state <- .nm_saem_conditional_components(
    context, parameters, eta, persistent = stochastic_context
  )
  modes <- lapply(seq_len(context$n_subjects), function(subject) {
    list(par = eta[subject, ], convergence = 0L, jitter = 0)
  })
  optimizer <- list(
    convergence = 0L,
    message = if (completed_iterations < n_iter) {
      "SAEM stationarity criterion reached"
    } else "SAEM iterations completed",
    counts = c(`function` = mstep_objective_evaluations,
               gradient = mstep_gradient_evaluations),
    iterations = completed_iterations,
    objective_evaluations = mstep_objective_evaluations,
    gradient_evaluations = mstep_gradient_evaluations,
    mstep_iterations = mstep_iterations, elapsed_seconds = mstep_elapsed,
    backend = paste0("saem+", mstep_backend),
    objective_backend = if (native_mstep_successes > 0L) {
      if (native_mstep_fallbacks > 0L) {
        "native-cpp-fixed-eta-population-objective-with-r-fallback"
      } else {
        "native-cpp-fixed-eta-population-objective"
      }
    } else {
      if (paired_objective_iterations > 0L) {
        "r-optimizer+native-paired-value-gradient"
      } else {
        "r-orchestrated-population-objective"
      }
    }
  )
  support_approximation <- .nm_liber_optimized(context) &&
    (saem_support_max > 0L || saem_support_prune > 0)
  accelerated_variant <- .nm_liber_optimized(context) &&
    (resolved_kernel == "fsaem" || saem_mstep_interval_burn > 1L ||
       saem_mstep_interval > 1L || resolved_parameter_averaging == "polyak")
  estimator_variant <- if (support_approximation) {
    "support-pruned approximate accelerated SAEM"
  } else if (accelerated_variant && resolved_kernel == "fsaem") {
    "accelerated f-SAEM"
  } else if (accelerated_variant) {
    "accelerated SAEM"
  } else {
    "canonical SAEM"
  }
  fit <- .nm_fit_result(
    context, "SAEM", parameters, final_objective_state$value, modes, optimizer,
    diagnostics = list(
      objective_trace = objective_trace, acceptance = accepted / max(attempted, 1L),
      acceptance_trace = acceptance_trace, step_scale_trace = step_scale_trace,
      parameter_trace = parameter_trace,
      estimator_identity = paste0(
        estimator_variant,
        " auxiliary-function stochastic approximation with ",
        if (resolved_kernel == "fsaem") "f-SAEM independence MH" else
          "random-walk MH"
      ),
      estimator_variant = estimator_variant,
      finite_iteration_schedule = list(
        mstep_every_iteration = saem_mstep_interval_burn == 1L &&
          saem_mstep_interval == 1L,
        parameter_averaging = resolved_parameter_averaging,
        support_approximation = support_approximation
      ),
      stochastic_approximation = list(
        state = "complete-data auxiliary function Q_k",
        retained_latent_support = q_state$count(),
        retained_weights = q_state$weights(),
        gain_power = sa_power,
        burn = burn,
        support_max = if (.nm_liber_optimized(context)) {
          saem_support_max
        } else 0L,
        support_prune_tolerance = if (.nm_liber_optimized(context)) {
          saem_support_prune
        } else 0,
        backend = q_state$telemetry()
      ),
      phase_timing = phases$snapshot(),
      stationarity = c(stationarity, list(
        auto_stop = isTRUE(auto_stop),
        stopped_early = completed_iterations < n_iter,
        consecutive_confirmations = stationary_iterations,
        required_confirmations = auto_stop_consecutive,
        minimum_iterations = auto_stop_min_iterations
      )),
      final_step_scale = step_scale, target_acceptance = target_acceptance,
      adaptive_proposal = isTRUE(adapt_proposal),
      eta_sampler = list(
        requested_kernel = saem_kernel,
        resolved_kernel = resolved_kernel,
        backend = if (!is.null(stochastic_context)) {
          "persistent cached conditional-subject C++ Metropolis"
        } else "cached batched conditional-subject C++ Metropolis",
        current_evaluations = eta_current_evaluations,
        current_cache_hits = eta_current_cache_hits,
        candidate_evaluations = eta_candidate_evaluations,
        proposal_root_cache_hits = as.integer(proposal_cache$hits),
        proposal_root_cache_misses = as.integer(proposal_cache$misses),
        proposal_root_factorizations = as.integer(
          proposal_cache$factorizations
        ),
        proposal_root_groups = length(unique(proposal_cache$groups)),
        fsaem = list(
          proposal = paste(
            "Laplace", resolved_fsaem_distribution,
            "independent Metropolis-Hastings"
          ),
          distribution = resolved_fsaem_distribution,
          degrees_of_freedom = if (
            resolved_fsaem_distribution == "student_t"
          ) fsaem_df else Inf,
          refresh_every = fsaem_refresh,
          refreshes = fsaem_refreshes,
          refresh_failures = fsaem_refresh_failures,
          fallback_iterations = fsaem_fallback_iterations,
          mode_iterations = fsaem_mode_iterations,
          mode_evaluations = fsaem_mode_evaluations,
          rescue_probability = fsaem_rescue_probability,
          rescue_iterations = fsaem_rescue_iterations,
          parameter_refresh_threshold = fsaem_parameter_refresh,
          parameter_triggered_refreshes = fsaem_parameter_refreshes,
          low_acceptance_threshold = fsaem_low_acceptance,
          acceptance_triggered_refreshes = fsaem_acceptance_refreshes,
          last_error = fsaem_last_error
        ),
        persistent_context = if (!is.null(stochastic_context)) {
          .liberation_stochastic_eta_context_telemetry(stochastic_context)
        } else NULL
      ),
      closed_form_sigma = simple_sigma,
      sigma_sufficient_statistics = list(
        gradient_recovery_updates = sigma_gradient_updates,
        simulation_fallbacks = sigma_expectation_fallbacks,
        backend = if (sigma_gradient_updates > 0L) {
          "exact-complete-data-gradient"
        } else if (simple_sigma) {
          "weighted-simulation"
        } else "not-applicable"
      ),
      closed_form_omega = length(map$omega_free) > 0L,
      native_sufficient_statistics = isTRUE(getOption(
        "LibeRation.saem_native_sufficient_statistics", TRUE
      )),
      paired_value_gradient = list(
        iterations = paired_objective_iterations,
        evaluations = paired_objective_evaluations,
        cache_hits = paired_objective_cache_hits,
        persistent_iterations = persistent_objective_iterations,
        persistent_fallbacks = persistent_objective_fallbacks,
        persistent_fallback_reason = persistent_objective_fallback_reason
      ),
      native_mstep = list(
        eligible = isTRUE(getOption("LibeRation.saem_native_mstep", TRUE)) &&
          .nm_liber_optimized(context) && is.null(context$parallel) &&
          optimizer_backend %in% c("auto", "native") && !isTRUE(mu$active),
        attempts = native_mstep_attempts,
        successes = native_mstep_successes,
        fallbacks = native_mstep_fallbacks,
        fallback_reason = native_mstep_fallback_reason,
        objective_backend = if (native_mstep_successes) {
          if (isTRUE(context$model$USE_ODE)) {
            "native-cpp-weighted-eta-ode-population-objective"
          } else {
            "native-cpp-weighted-eta-population-objective"
          }
        } else {
          if (paired_objective_iterations > 0L) {
            "r-optimizer+native-paired-value-gradient"
          } else {
            "r-orchestrated-population-objective"
          }
        }
      ),
      mu_specialization = c(
        .nm_mu_diagnostic(mu),
        list(
          closed_form_updates = mu_updates,
          closed_form_only_iterations = mu_closed_form_only_iterations,
          runtime_fallbacks = mu_fallbacks,
          runtime_reason = mu$runtime_reason %||% NULL
        )
      ),
      n_iter = completed_iterations, requested_n_iter = n_iter,
      mstep_schedule = list(
        burn_interval = saem_mstep_interval_burn,
        post_burn_interval = saem_mstep_interval,
        performed = mstep_performed_trace[used],
        count = sum(mstep_performed_trace[used])
      ),
      parameter_averaging = list(
        requested = saem_parameter_averaging,
        resolved = resolved_parameter_averaging,
        start_iteration = saem_average_start,
        samples = parameter_average_count,
        applied = parameter_average_applied
      ),
      burn = burn, seed = seed,
      population_gradient = "exact conditional CppAD gradient"
    )
  )
  fit
}

.nm_est_saem <- function(context, map, maxit, tolerance, trace,
                         n_replicates = 1L,
                         replicate_seed_stride = 100003L,
                         parallel_replicates = FALSE,
                         replicate_score_samples = 200L,
                         replicate_score_seed = 20260811L,
                         replicate_score_eta_maxit = 100L, ...) {
  n_replicates <- as.integer(n_replicates)
  replicate_seed_stride <- as.integer(replicate_seed_stride)
  replicate_score_samples <- as.integer(replicate_score_samples)
  replicate_score_seed <- as.integer(replicate_score_seed)
  replicate_score_eta_maxit <- as.integer(replicate_score_eta_maxit)
  if (is.na(n_replicates) || n_replicates < 1L ||
      is.na(replicate_seed_stride) || replicate_seed_stride < 1L) {
    .nm_stop("SAEM replicate count and seed stride must be positive integers.")
  }
  if (is.na(replicate_score_samples) || replicate_score_samples < 5L ||
      is.na(replicate_score_seed) ||
      is.na(replicate_score_eta_maxit) || replicate_score_eta_maxit < 1L) {
    .nm_stop(
      "SAEM replicate scoring requires at least five importance samples, ",
      "one integer seed, and a positive ETA-mode iteration limit."
    )
  }
  if (isTRUE(parallel_replicates) && n_replicates > 1L) {
    warning(
      "SAEM replicates are sequenced in this R session because compiled tape ",
      "contexts are not shared across processes. Submit independently seeded ",
      "jobs to a LibeR queue for replicate-level parallelism.", call. = FALSE
    )
  }
  controls <- list(...)
  base_seed <- as.integer(controls$seed %||% 20260713L)
  controls$seed <- NULL
  fits <- lapply(seq_len(n_replicates), function(replicate) {
    do.call(.nm_est_saem_single, c(list(
      context = context, map = map, maxit = maxit,
      tolerance = tolerance, trace = trace,
      seed = base_seed + (replicate - 1L) * replicate_seed_stride
    ), controls))
  })
  conditional_objectives <- vapply(fits, `[[`, numeric(1), "objective")
  selection_scores <- rep(Inf, n_replicates)
  selection_metric <- paste0(
    "common-random-number antithetic importance marginal objective (",
    replicate_score_samples, " samples per subject)"
  )
  score_errors <- rep(NA_character_, n_replicates)
  # All replicates are scored with the same base-normal draws. Each fit still
  # constructs its own correctly normalized conditional proposal, making the
  # resulting values comparable finite-sample marginal likelihood estimates
  # rather than last-ETA complete-data objectives. The single-replicate case
  # uses the same calculation so the public SAEM objective has clear likelihood
  # semantics too.
  score_normals <- .nm_imp_normals(
    context, replicate_score_samples, replicate_score_seed,
    sampling = "antithetic", proposal = "gaussian"
  )
  selection_scores <- vapply(seq_along(fits), function(index) {
    fit <- fits[[index]]
    parameters <- list(
      theta = fit$theta, sigma = fit$sigma, omega = fit$omega
    )
    score <- tryCatch(
      .nm_imp_evaluate(
        context, parameters, score_normals,
        replicate_score_eta_maxit, tolerance, gradient = FALSE,
        initial_eta = fit$eta
      )$value,
      error = identity
    )
    if (inherits(score, "error")) {
      score_errors[[index]] <<- conditionMessage(score)
      return(Inf)
    }
    as.numeric(score)
  }, numeric(1))
  best <- which.min(ifelse(is.finite(selection_scores), selection_scores, Inf))
  if (!length(best) || !is.finite(selection_scores[[best]])) {
    warning(
      "All independent SAEM replicate marginal scores failed; selecting the ",
      "replicate with the smallest reported conditional objective. Inspect ",
      "`diagnostics$replicates$score_errors`.", call. = FALSE
    )
    best <- which.min(ifelse(
      is.finite(conditional_objectives), conditional_objectives, Inf
    ))
    if (!length(best) || !is.finite(conditional_objectives[[best]])) best <- 1L
    selection_metric <- "conditional objective fallback after failed marginal scoring"
  }
  result <- fits[[best]]
  marginal_score_available <- is.finite(selection_scores[[best]]) &&
    !grepl("fallback", selection_metric, fixed = TRUE)
  if (marginal_score_available) {
    result$objective <- selection_scores[[best]]
    result$objective_type <-
      "negative_twice_importance_sampled_marginal_log_likelihood"
    result$objective_comparable <- TRUE
    result$diagnostics$objective_semantics <- list(
      type = result$objective_type,
      likelihood_comparable = TRUE,
      reported_point = "selected SAEM population estimate",
      Monte_Carlo_samples_per_subject = replicate_score_samples,
      Monte_Carlo_seed = replicate_score_seed
    )
  } else {
    result$objective_type <-
      "negative_twice_complete_data_objective_at_final_latent_state"
    result$objective_comparable <- FALSE
    result$diagnostics$objective_semantics <- list(
      type = result$objective_type,
      likelihood_comparable = FALSE,
      recommendation = "Repeat marginal scoring before likelihood comparison."
    )
  }
  estimates <- do.call(rbind, lapply(fits, function(value) {
    c(value$theta, value$sigma, value$omega)
  }))
  colnames(estimates) <- .nm_parameter_names(
    result$theta, result$sigma, result$omega
  )
  result$diagnostics$replicates <- list(
    count = n_replicates,
    selected = best,
    seeds = base_seed + (seq_len(n_replicates) - 1L) * replicate_seed_stride,
    objectives = conditional_objectives,
    conditional_objectives = conditional_objectives,
    selection_metric = selection_metric,
    selection_scores = selection_scores,
    score_samples = replicate_score_samples,
    score_seed = replicate_score_seed,
    score_errors = score_errors,
    estimates = estimates,
    between_replicate_sd = if (n_replicates > 1L) {
      apply(estimates, 2L, stats::sd)
    } else stats::setNames(rep(NA_real_, ncol(estimates)), colnames(estimates)),
    stationarity = lapply(fits, function(value) {
      value$diagnostics$stationarity
    })
  )
  result
}

