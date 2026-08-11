.nm_log_mean_exp <- function(values) {
  maximum <- max(values)
  maximum + log(mean(exp(values - maximum)))
}

.nm_log_sum_exp <- function(values) {
  maximum <- max(values)
  maximum + log(sum(exp(values - maximum)))
}

.nm_stochastic_phase_timer <- function() {
  state <- new.env(parent = emptyenv())
  state$seconds <- numeric()
  state$calls <- integer()
  time <- function(phase, expression) {
    started <- proc.time()[["elapsed"]]
    on.exit({
      elapsed <- unname(proc.time()[["elapsed"]] - started)
      previous_seconds <- if (phase %in% names(state$seconds)) {
        state$seconds[[phase]]
      } else 0
      previous_calls <- if (phase %in% names(state$calls)) {
        state$calls[[phase]]
      } else 0L
      state$seconds[[phase]] <- previous_seconds + elapsed
      state$calls[[phase]] <- previous_calls + 1L
    }, add = TRUE)
    force(expression)
  }
  list(
    time = time,
    snapshot = function() list(
      seconds = as.list(state$seconds), calls = as.list(state$calls),
      total_seconds = sum(state$seconds)
    )
  )
}

.nm_weighted_eta_context <- function(context, allow_compatibility = TRUE,
                                     reduced_population_tape = NULL) {
  if ((!isTRUE(allow_compatibility) && !.nm_liber_optimized(context)) ||
      !is.null(context$parallel) ||
      !isTRUE(getOption("LibeRation.weighted_eta_context", TRUE))) {
    return(NULL)
  }
  if (is.null(reduced_population_tape)) {
    reduced_population_tape <- .nm_liber_optimized(context)
  }
  tryCatch(
    .liberation_weighted_eta_context_create(
      context$engine$pointer,
      lapply(context$subjects, function(evaluator) {
        evaluator$objective_tape$pointer
      }),
      lapply(context$subjects, function(evaluator) evaluator$data_input()),
      length(context$model$THETAS$Value), context$n_eta,
      length(context$model$SIGMAS$Value),
      length(context$model$OMEGAS$Value),
      isTRUE(context$model$USE_ODE),
      isTRUE(reduced_population_tape) && .nm_liber_optimized(context),
      as.integer(context$native_subject_threads %||% 1L),
      as.integer(getOption("LibeRation.ode_weighted_tape_limit", 4096L))
    ),
    error = function(error) NULL
  )
}

.nm_gq_tensor_fits <- function(order, dimension, max_points) {
  points <- 1
  for (axis in seq_len(dimension)) {
    if (points > max_points / order) return(FALSE)
    points <- points * order
  }
  TRUE
}

.nm_gq_design <- function(context, order = 5L, max_points = 100000L,
                          grid = c("auto", "tensor", "smolyak"), level = 3L) {
  order <- as.integer(order)
  level <- as.integer(level)
  max_points <- as.integer(max_points)
  if (!length(grid)) .nm_stop("`gq_grid` must contain one grid strategy.")
  grid <- tolower(as.character(grid[[1L]] %||% "auto"))
  if (length(grid) != 1L || is.na(grid)) {
    .nm_stop("`gq_grid` must contain one grid strategy.")
  }
  if (identical(grid, "sparse")) grid <- "smolyak"
  if (!grid %in% c("auto", "tensor", "smolyak")) {
    .nm_stop("`gq_grid` must be one of auto, tensor, or smolyak.")
  }
  requested_grid <- grid
  if (grid == "auto") {
    tensor_fits <- .nm_gq_tensor_fits(order, context$n_eta, max_points)
    grid <- if (tensor_fits && (context$n_eta <= 3L || order == 1L)) {
      "tensor"
    } else "smolyak"
  }
  rule <- if (grid == "tensor") {
    LibeRtAD::ad_gauss_hermite(
      order = order, dimension = context$n_eta, max_points = max_points
    )
  } else {
    LibeRtAD::ad_smolyak_gauss_hermite(
      level = level, dimension = context$n_eta, max_points = max_points
    )
  }
  nodes <- rule$nodes
  attr(nodes, "log_measure") <- as.numeric(
    rule$log_abs_weights %||% rule$log_weights
  )
  attr(nodes, "measure_sign") <- as.numeric(rule$signs %||% sign(rule$weights))
  attr(nodes, "quadrature_method") <- paste0(grid, "-gauss-hermite")
  list(
    normals = rep(list(nodes), context$n_subjects),
    method = paste0(grid, "-gauss-hermite"),
    actual_samples = as.integer(rule$points),
    candidate_points = as.integer(rule$candidate_points %||% rule$points),
    quadrature_order = if (grid == "tensor") as.integer(rule$order) else NA_integer_,
    quadrature_level = if (grid == "smolyak") as.integer(rule$level) else NA_integer_,
    requested_grid = requested_grid,
    resolved_grid = grid,
    negative_weights = as.integer(rule$negative_weights %||% 0L),
    max_points = max_points
  )
}

.nm_radical_inverse <- function(index, base) {
  index <- as.integer(index)
  result <- numeric(length(index))
  factor <- 1 / base
  while (any(index > 0L)) {
    result <- result + factor * (index %% base)
    index <- index %/% base
    factor <- factor / base
  }
  result
}

.nm_imp_rqmc_uniform <- function(draws, dimension, offset = 0L) {
  if (!dimension) return(matrix(numeric(), draws, 0L))
  primes <- c(
    2L, 3L, 5L, 7L, 11L, 13L, 17L, 19L, 23L, 29L, 31L, 37L,
    41L, 43L, 47L, 53L, 59L, 61L, 67L, 71L, 73L, 79L, 83L,
    89L, 97L, 101L, 103L, 107L, 109L, 113L, 127L, 131L
  )
  if (dimension > length(primes)) {
    .nm_stop("Randomized quasi-Monte Carlo IMP currently supports at most ",
             length(primes), " ETA dimensions.")
  }
  sequence <- offset + seq_len(draws)
  result <- vapply(seq_len(dimension), function(axis) {
    (.nm_radical_inverse(sequence, primes[[axis]]) + stats::runif(1L)) %% 1
  }, numeric(draws))
  matrix(pmin(1 - 1e-12, pmax(1e-12, result)), draws, dimension)
}

.nm_imp_normals <- function(context, n_imp, seed,
                            sampling = c("random", "antithetic", "rqmc"),
                            proposal = c("gaussian", "student_t", "defensive"),
                            proposal_df = 7) {
  sampling <- match.arg(sampling)
  proposal <- match.arg(proposal)
  n_imp <- as.integer(n_imp)
  seed <- as.integer(seed)
  if (length(n_imp) == 1L) n_imp <- rep.int(n_imp, context$n_subjects)
  if (length(n_imp) != context$n_subjects || anyNA(n_imp) || any(n_imp < 5L)) {
    .nm_stop("Importance-sampling information requires at least 5 samples.")
  }
  if (length(seed) != 1L || is.na(seed)) .nm_stop("`seed` must be one integer.")
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) previous_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  lapply(seq_len(context$n_subjects), function(subject) {
    draws <- n_imp[[subject]]
    dimension <- context$n_eta
    draw_gaussian <- function(draws) {
      if (sampling == "rqmc") {
        uniforms <- .nm_imp_rqmc_uniform(
          draws, dimension, offset = (subject - 1L) * max(n_imp)
        )
        return(matrix(stats::qnorm(uniforms), draws, dimension))
      }
      if (sampling == "antithetic" && draws > 1L) {
        half <- as.integer(ceiling(draws / 2))
        base <- matrix(stats::rnorm(half * dimension), half, dimension)
        rbind(base, -base)[seq_len(draws), , drop = FALSE]
      } else matrix(stats::rnorm(draws * dimension), draws, dimension)
    }
    draw_student_t <- function(draws) {
      scale <- sqrt((proposal_df - 2) / proposal_df)
      if (sampling == "rqmc") {
        uniforms <- .nm_imp_rqmc_uniform(
          draws, dimension, offset = (subject - 1L) * max(n_imp)
        )
        return(matrix(stats::qt(uniforms, df = proposal_df), draws, dimension) *
                 scale)
      }
      if (sampling == "antithetic" && draws > 1L) {
        half <- as.integer(ceiling(draws / 2))
        base <- matrix(
          stats::rt(half * dimension, df = proposal_df),
          half, dimension
        ) * scale
        rbind(base, -base)[seq_len(draws), , drop = FALSE]
      } else {
        matrix(
          stats::rt(draws * dimension, df = proposal_df),
          draws, dimension
        ) * scale
      }
    }
    if (proposal == "gaussian") {
      z <- draw_gaussian(draws)
    } else if (proposal == "student_t") {
      z <- draw_student_t(draws)
    } else {
      gaussian_draws <- as.integer(ceiling(draws / 2))
      t_draws <- draws - gaussian_draws
      z <- rbind(
        draw_gaussian(gaussian_draws),
        if (t_draws) draw_student_t(t_draws) else
          matrix(numeric(), 0L, dimension)
      )
      # Deterministic-mixture importance sampling must evaluate the same
      # component allocation that was drawn. Odd budgets are not exactly
      # 50:50, so retain the realised Gaussian fraction for the balance
      # heuristic used when constructing the proposal density.
      attr(z, "imp_defensive_gaussian_weight") <- gaussian_draws / draws
    }
    attr(z, "imp_proposal") <- proposal
    attr(z, "imp_proposal_df") <- proposal_df
    attr(z, "imp_sampling") <- sampling
    z
  })
}

.nm_imp_covariance_design <- function(context, samples, seed) {
  samples <- as.integer(samples)
  dimension <- as.integer(context$n_eta)
  if (!dimension) {
    normals <- rep(list(matrix(numeric(), 1L, 0L)), context$n_subjects)
    return(list(
      normals = normals, method = "none", actual_samples = 1L,
      quadrature_order = 0L
    ))
  }
  order <- min(15L, max(3L, as.integer(ceiling(samples^(1 / dimension)))))
  nodes_required <- order^dimension
  use_quadrature <- nodes_required <= max(4L * samples, 1024L)
  if (!use_quadrature) {
    return(list(
      normals = .nm_imp_normals(context, samples, seed),
      method = "random-normal", actual_samples = samples,
      quadrature_order = NA_integer_
    ))
  }
  .nm_gq_design(context, order = order, max_points = nodes_required)
}

.nm_complete_data_expectation <- function(context, map, eta, weights,
                                          native_context = NULL) {
  if (length(eta) != context$n_subjects || length(weights) != context$n_subjects) {
    .nm_stop("Complete-data expectation requires one ETA grid and weight vector per subject.")
  }
  for (subject in seq_len(context$n_subjects)) {
    eta[[subject]] <- as.matrix(eta[[subject]])
    weights[[subject]] <- as.numeric(weights[[subject]])
    if (ncol(eta[[subject]]) != context$n_eta ||
        nrow(eta[[subject]]) != length(weights[[subject]]) ||
        any(!is.finite(eta[[subject]])) || any(!is.finite(weights[[subject]])) ||
        any(weights[[subject]] < 0) || sum(weights[[subject]]) <= 0) {
      .nm_stop("Complete-data ETA grids and probability weights are invalid.")
    }
    weights[[subject]] <- weights[[subject]] / sum(weights[[subject]])
  }
  native_error <- NULL
  native_ready <- FALSE
  if (!is.null(native_context)) {
    native_ready <- tryCatch({
      .liberation_weighted_eta_context_set(native_context, eta, weights)
      TRUE
    }, error = function(error) {
      native_error <<- conditionMessage(error)
      FALSE
    })
  }
  cache <- new.env(parent = emptyenv())
  cache$key <- NULL
  cache$result <- NULL
  cache$evaluations <- 0L
  cache$hits <- 0L
  evaluate <- function(parameters) {
    key <- map$encode(parameters)
    if (!is.null(cache$key) && identical(key, cache$key)) {
      cache$hits <- cache$hits + 1L
      return(cache$result)
    }
    if (native_ready) {
      evaluated <- tryCatch(
        .liberation_weighted_eta_context_eval(
          native_context, parameters$theta, parameters$sigma,
          parameters$omega
        ),
        error = identity
      )
      if (inherits(evaluated, "error")) {
        native_error <<- conditionMessage(evaluated)
        native_ready <<- FALSE
      }
    }
    if (native_ready) {
      value <- as.numeric(evaluated$value)
      full_gradient <- as.numeric(evaluated$gradient)
    } else {
      value <- 0
      full_gradient <- numeric(
        length(parameters$theta) + context$n_eta +
          length(parameters$sigma) + length(parameters$omega)
      )
      for (subject in seq_len(context$n_subjects)) {
        evaluated <- context$subjects[[subject]]$objective_eta_batch(
          parameters$theta, eta[[subject]], parameters$sigma, parameters$omega
        )
        value <- value + sum(weights[[subject]] * evaluated$value)
        full_gradient <- full_gradient + colSums(
          evaluated$gradient * weights[[subject]]
        )
      }
    }
    n_theta <- length(parameters$theta)
    n_sigma <- length(parameters$sigma)
    population_positions <- c(
      seq_len(n_theta), n_theta + context$n_eta + seq_len(n_sigma),
      n_theta + context$n_eta + n_sigma + seq_along(parameters$omega)
    )
    native_gradient <- as.numeric(full_gradient[population_positions]) +
      .nm_prior_nll_native_gradient(context$model, parameters)
    cache$result <- list(
      value = value + .nm_prior_nll(context$model, parameters),
      gradient = as.vector(native_gradient %*% map$jacobian(parameters)),
      native_gradient = as.numeric(full_gradient)
    )
    cache$key <- key
    cache$evaluations <- cache$evaluations + 1L
    cache$result
  }
  list(
    objective = function(parameters) evaluate(parameters)$value,
    gradient = function(parameters) evaluate(parameters)$gradient,
    native_gradient = function(parameters) evaluate(parameters)$native_gradient,
    telemetry = function() list(
      evaluations = cache$evaluations, cache_hits = cache$hits,
      native = native_ready,
      native_error = native_error,
      native_context = if (!is.null(native_context)) {
        tryCatch(
          .liberation_weighted_eta_context_telemetry(native_context),
          error = function(error) NULL
        )
      } else NULL
    )
  )
}

.nm_native_weighted_expectation <- function(context, map, native_context) {
  cache <- new.env(parent = emptyenv())
  cache$key <- NULL
  cache$result <- NULL
  cache$evaluations <- 0L
  cache$hits <- 0L
  evaluate <- function(parameters) {
    key <- map$encode(parameters)
    if (!is.null(cache$key) && identical(key, cache$key)) {
      cache$hits <- cache$hits + 1L
      return(cache$result)
    }
    evaluated <- .liberation_weighted_eta_context_eval(
      native_context, parameters$theta, parameters$sigma, parameters$omega
    )
    n_theta <- length(parameters$theta)
    n_sigma <- length(parameters$sigma)
    population_positions <- c(
      seq_len(n_theta), n_theta + context$n_eta + seq_len(n_sigma),
      n_theta + context$n_eta + n_sigma + seq_along(parameters$omega)
    )
    native_gradient <- as.numeric(evaluated$gradient[population_positions]) +
      .nm_prior_nll_native_gradient(context$model, parameters)
    cache$result <- list(
      value = as.numeric(evaluated$value) +
        .nm_prior_nll(context$model, parameters),
      gradient = as.vector(native_gradient %*% map$jacobian(parameters)),
      native_gradient = as.numeric(evaluated$gradient)
    )
    cache$key <- key
    cache$evaluations <- cache$evaluations + 1L
    cache$result
  }
  list(
    objective = function(parameters) evaluate(parameters)$value,
    gradient = function(parameters) evaluate(parameters)$gradient,
    native_gradient = function(parameters) evaluate(parameters)$native_gradient,
    telemetry = function() list(
      paired = TRUE, evaluations = cache$evaluations,
      cache_hits = cache$hits, persistent = TRUE,
      native_aggregate = TRUE, persistent_requested = TRUE,
      persistent_error = NULL,
      weighted_context = .liberation_weighted_eta_context_telemetry(
        native_context
      )
    )
  )
}

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

.nm_imp_subject_objective <- function(evaluator, parameters, normals,
                                      eta_maxit, tolerance) {
  .nm_imp_subject_state(
    evaluator, parameters, normals, eta_maxit, tolerance, gradient = FALSE
  )$value
}

.nm_imp_proposal_from_mode <- function(evaluator, parameters, normals, mode) {
  if (mode$convergence != 0L) {
    return(list(valid = FALSE, mode = mode, eta = NULL, log_proposal = NULL))
  }
  dimension <- length(mode$par)
  if (!dimension) {
    return(list(
      valid = TRUE, mode = mode, eta = matrix(numeric(), 1L, 0L),
      log_proposal = 0, log_measure = 0, sampling = "none"
    ))
  }
  covariance <- 2 * solve(mode$hessian)
  covariance <- .nm_positive_definite(covariance, "IMP proposal covariance")$matrix
  root <- t(chol(covariance))
  logdet <- as.numeric(determinant(covariance, logarithm = TRUE)$modulus)
  z <- normals
  proposal_family <- attr(normals, "imp_proposal", exact = TRUE) %||%
    "gaussian"
  proposal_df <- as.numeric(
    attr(normals, "imp_proposal_df", exact = TRUE) %||% 7
  )
  defensive_gaussian_weight <- as.numeric(
    attr(normals, "imp_defensive_gaussian_weight", exact = TRUE) %||% 0.5
  )
  log_measure <- attr(normals, "log_measure", exact = TRUE)
  measure_sign <- attr(normals, "measure_sign", exact = TRUE)
  sampling <- attr(normals, "quadrature_method", exact = TRUE)
  if (is.null(log_measure)) sampling <- "random-normal"
  if (is.null(sampling)) sampling <- "tensor-gauss-hermite"
  if (is.null(log_measure)) log_measure <- rep(-log(nrow(z)), nrow(z))
  if (is.null(measure_sign)) measure_sign <- rep(1, nrow(z))
  gaussian_log <- -0.5 * (
    dimension * log(2 * pi) + logdet + rowSums(z^2)
  )
  t_scale <- (proposal_df - 2) / proposal_df
  t_quadratic <- rowSums(z^2) / t_scale
  t_log <- lgamma((proposal_df + dimension) / 2) -
    lgamma(proposal_df / 2) - dimension * log(proposal_df * pi) / 2 -
    (logdet + dimension * log(t_scale)) / 2 -
    (proposal_df + dimension) * log1p(t_quadratic / proposal_df) / 2
  log_proposal <- switch(
    proposal_family,
    gaussian = gaussian_log,
    student_t = t_log,
    defensive = {
      if (length(defensive_gaussian_weight) != 1L ||
          !is.finite(defensive_gaussian_weight) ||
          defensive_gaussian_weight <= 0 || defensive_gaussian_weight >= 1) {
        .nm_stop("The defensive IMP proposal allocation is invalid.")
      }
      maximum <- pmax(gaussian_log, t_log)
      maximum + log(
        defensive_gaussian_weight * exp(gaussian_log - maximum) +
          (1 - defensive_gaussian_weight) * exp(t_log - maximum)
      )
    },
    .nm_stop("Unknown IMP proposal family: ", proposal_family)
  )
  list(
    valid = TRUE, mode = mode,
    eta = sweep(z %*% t(root), 2L, mode$par, `+`),
    log_proposal = log_proposal,
    log_measure = as.numeric(log_measure),
    measure_sign = as.numeric(measure_sign), sampling = sampling,
    proposal_family = proposal_family, proposal_df = proposal_df
  )
}

.nm_imp_subject_proposal <- function(evaluator, parameters, normals,
                                     eta_maxit, tolerance, start = NULL) {
  if (!is.null(start)) {
    start <- as.numeric(start)
    if (length(start) != evaluator$n_eta || any(!is.finite(start))) {
      .nm_stop("Adaptive proposal ETA starts must match the subject ETA dimension and be finite.")
    }
  }
  mode <- evaluator$eta_mode(
    parameters$theta, parameters$sigma, parameters$omega,
    start = start %||% rep(0, evaluator$n_eta),
    maxit = eta_maxit, tolerance = tolerance
  )
  .nm_imp_proposal_from_mode(evaluator, parameters, normals, mode)
}

.nm_gq_fixed_subject_proposal <- function(evaluator, parameters, normals) {
  dimension <- evaluator$n_eta
  mode <- list(
    par = rep(0, dimension), convergence = 0L, iterations = 0L,
    evaluations = 0L, backend = "fixed-omega"
  )
  if (!dimension) {
    return(list(
      valid = TRUE, mode = mode, eta = matrix(numeric(), 1L, 0L),
      log_proposal = 0, log_measure = 0,
      sampling = "fixed-tensor-gauss-hermite", measure_sign = 1
    ))
  }
  covariance <- .nm_positive_definite(
    .nm_omega_matrix(evaluator$engine$model, parameters$omega),
    "Fixed GQ OMEGA covariance"
  )$matrix
  root <- t(chol(covariance))
  logdet <- as.numeric(determinant(covariance, logarithm = TRUE)$modulus)
  z <- normals
  log_measure <- attr(normals, "log_measure", exact = TRUE)
  measure_sign <- attr(normals, "measure_sign", exact = TRUE)
  sampling <- attr(normals, "quadrature_method", exact = TRUE) %||%
    "tensor-gauss-hermite"
  if (is.null(log_measure)) {
    .nm_stop("Fixed Gaussian quadrature requires deterministic node weights.")
  }
  if (is.null(measure_sign)) measure_sign <- rep(1, nrow(z))
  list(
    valid = TRUE, mode = mode,
    eta = z %*% t(root),
    log_proposal = -0.5 * (
      dimension * log(2 * pi) + logdet + rowSums(z^2)
    ),
    log_measure = as.numeric(log_measure),
    measure_sign = as.numeric(measure_sign),
    sampling = paste0("fixed-", sampling)
  )
}

.nm_imp_subject_from_proposal <- function(evaluator, parameters, proposal,
                                          gradient = TRUE) {
  if (!isTRUE(proposal$valid)) {
    return(list(value = Inf, native_gradient = NULL, mode = proposal$mode))
  }
  dimension <- ncol(proposal$eta)
  if (!dimension) {
    evaluated <- evaluator$objective(
      parameters$theta, numeric(), parameters$sigma, parameters$omega,
      gradient = gradient
    )
    return(list(
      value = evaluated$value,
      native_gradient = if (isTRUE(gradient)) as.numeric(evaluated$gradient) else NULL,
      mode = proposal$mode, effective_sample_size = 1
    ))
  }
  evaluated <- if (isTRUE(gradient)) {
    evaluator$objective_eta_batch(
      parameters$theta, proposal$eta, parameters$sigma, parameters$omega
    )
  } else list(value = evaluator$objective_eta_values(
    parameters$theta, proposal$eta, parameters$sigma, parameters$omega
  ))
  log_weight <- -0.5 * evaluated$value - proposal$log_proposal
  log_integrand <- log_weight + proposal$log_measure
  measure_sign <- proposal$measure_sign %||% rep(1, length(log_integrand))
  finite <- is.finite(log_integrand) & is.finite(measure_sign) & measure_sign != 0
  if (!any(finite)) {
    return(list(
      value = Inf, native_gradient = NULL, mode = proposal$mode,
      effective_sample_size = 0, cancellation_ratio = 0,
      quadrature_valid = FALSE
    ))
  }
  maximum <- max(log_integrand[finite])
  scaled <- numeric(length(log_integrand))
  scaled[finite] <- measure_sign[finite] * exp(log_integrand[finite] - maximum)
  signed_total <- sum(scaled)
  absolute_total <- sum(abs(scaled))
  valid <- is.finite(signed_total) && is.finite(absolute_total) &&
    signed_total > .Machine$double.eps * max(1, absolute_total)
  if (!valid) {
    return(list(
      value = Inf, native_gradient = NULL, mode = proposal$mode,
      effective_sample_size = 0, cancellation_ratio = 0,
      quadrature_valid = FALSE
    ))
  }
  value <- -2 * (maximum + log(signed_total))
  native_gradient <- NULL
  absolute_weights <- abs(scaled) / absolute_total
  effective_sample_size <- 1 / sum(absolute_weights^2)
  cancellation_ratio <- signed_total / absolute_total
  if (isTRUE(gradient)) {
    weights <- scaled / signed_total
    native_gradient <- colSums(evaluated$gradient * weights)
  }
  list(
    value = value, native_gradient = native_gradient, mode = proposal$mode,
    effective_sample_size = effective_sample_size,
    cancellation_ratio = cancellation_ratio, quadrature_valid = TRUE
  )
}

.nm_imp_subject_state <- function(evaluator, parameters, normals,
                                  eta_maxit, tolerance, gradient = TRUE) {
  proposal <- .nm_imp_subject_proposal(
    evaluator, parameters, normals, eta_maxit, tolerance
  )
  .nm_imp_subject_from_proposal(evaluator, parameters, proposal, gradient)
}

.nm_imp_prepare_proposals <- function(context, parameters, normals,
                                      eta_maxit, tolerance, adaptive = TRUE,
                                      initial_eta = NULL,
                                      cached_modes = NULL,
                                      proposal_curvature = c("exact", "fisher")) {
  proposal_curvature <- match.arg(proposal_curvature)
  if (!is.null(initial_eta)) {
    initial_eta <- as.matrix(initial_eta)
    expected <- c(context$n_subjects, context$n_eta)
    if (!identical(dim(initial_eta), expected) || any(!is.finite(initial_eta))) {
      .nm_stop(
        "`initial_eta` must be a finite ", expected[[1L]], " x ",
        expected[[2L]], " subject-by-ETA matrix."
      )
    }
  }
  if (!context$n_eta) {
    # With no random effects the marginal subject contribution is identical
    # to the conditional contribution.  In particular, do not ask the
    # Fisher-proposal path to invert a 0 x 0 covariance matrix.
    return(Map(function(evaluator, normal) {
      .nm_gq_fixed_subject_proposal(evaluator, parameters, normal)
    }, context$subjects, normals))
  }
  prepare_chunk <- function(evaluators, chunk_normals, chunk_starts = NULL) {
    lapply(seq_along(evaluators), function(subject) {
      if (isTRUE(adaptive)) {
        .nm_imp_subject_proposal(
          evaluators[[subject]], parameters, chunk_normals[[subject]],
          eta_maxit, tolerance,
          start = if (is.null(chunk_starts)) NULL else chunk_starts[subject, ]
        )
      } else {
        .nm_gq_fixed_subject_proposal(
          evaluators[[subject]], parameters, chunk_normals[[subject]]
        )
      }
    })
  }
  if (!is.null(cached_modes)) {
    if (length(cached_modes) != context$n_subjects) {
      .nm_stop("Cached IMP modes must contain one mode per subject.")
    }
    return(Map(function(evaluator, normal, mode) {
      .nm_imp_proposal_from_mode(evaluator, parameters, normal, mode)
    }, context$subjects, normals, cached_modes))
  }
  if (isTRUE(adaptive) &&
      is.null(context$parallel) && !isTRUE(context$model$USE_ODE)) {
    starts <- initial_eta %||% matrix(0, context$n_subjects, context$n_eta)
    modes <- if (proposal_curvature == "fisher") {
      distribution <- .nm_its_distribution(
        context, parameters, starts, eta_maxit, tolerance
      )
      Map(function(mode, covariance) {
        # Importance weights remain exact; this Fisher/Gauss-Newton curvature
        # changes proposal efficiency only, not the MCEM target.
        mode$hessian <- 2 * solve(.nm_positive_definite(
          covariance, "IMP Fisher proposal covariance"
        )$matrix)
        mode
      }, distribution$modes, distribution$covariance)
    } else {
      .nm_subject_modes(
        context, parameters, starts = starts, maxit = eta_maxit,
        tolerance = tolerance, interaction = TRUE, exact_hessian = TRUE
      )
    }
    return(Map(function(evaluator, normal, mode) {
      .nm_imp_proposal_from_mode(evaluator, parameters, normal, mode)
    }, context$subjects, normals, modes))
  }
  if (is.null(context$parallel)) {
    return(prepare_chunk(context$subjects, normals, initial_eta))
  }
  chunks <- context$parallel$chunks
  normal_chunks <- lapply(chunks, function(rows) normals[rows])
  start_chunks <- if (is.null(initial_eta)) {
    rep(list(NULL), length(chunks))
  } else {
    lapply(chunks, function(rows) initial_eta[rows, , drop = FALSE])
  }
  pieces <- parallel::clusterApply(
    context$parallel$cluster, seq_along(chunks),
      function(index, normal_chunks, start_chunks, parameters,
               eta_maxit, tolerance, adaptive) {
        namespace <- asNamespace("LibeRation")
        evaluators <- get(".nm_parallel_worker_state", envir = namespace)()$subjects
        prepare <- get(
          if (isTRUE(adaptive)) ".nm_imp_subject_proposal" else
            ".nm_gq_fixed_subject_proposal",
          envir = asNamespace("LibeRation")
        )
        lapply(seq_along(evaluators), function(subject) {
          if (isTRUE(adaptive)) {
            prepare(
              evaluators[[subject]], parameters,
              normal_chunks[[index]][[subject]], eta_maxit, tolerance,
              start = if (is.null(start_chunks[[index]])) {
                NULL
              } else start_chunks[[index]][subject, ]
            )
          } else {
            prepare(
              evaluators[[subject]], parameters,
              normal_chunks[[index]][[subject]]
            )
          }
        })
      }, normal_chunks = normal_chunks, start_chunks = start_chunks,
      parameters = parameters,
      eta_maxit = eta_maxit, tolerance = tolerance, adaptive = adaptive
  )
  unlist(pieces, recursive = FALSE)
}

.nm_imp_evaluate_fixed <- function(context, parameters, proposals,
                                   gradient = TRUE) {
  evaluate_chunk <- function(evaluators, chunk_proposals) {
    lapply(seq_along(evaluators), function(subject) {
      .nm_imp_subject_from_proposal(
        evaluators[[subject]], parameters, chunk_proposals[[subject]], gradient
      )
    })
  }
  if (is.null(context$parallel) &&
      !isTRUE(context$model$USE_ODE) && length(proposals) &&
      all(vapply(proposals, function(proposal) isTRUE(proposal$valid), logical(1)))) {
    points <- cbind(
      matrix(parameters$theta, context$n_subjects, length(parameters$theta),
             byrow = TRUE),
      matrix(0, context$n_subjects, context$n_eta),
      matrix(parameters$sigma, context$n_subjects, length(parameters$sigma),
             byrow = TRUE),
      matrix(parameters$omega, context$n_subjects, length(parameters$omega),
             byrow = TRUE)
    )
    native <- .liberation_objective_tape_importance_collection(
      lapply(context$subjects, function(evaluator) evaluator$objective_tape$pointer),
      points, length(parameters$theta) + seq_len(context$n_eta),
      lapply(proposals, `[[`, "eta"),
      lapply(proposals, `[[`, "log_proposal"),
      lapply(proposals, `[[`, "log_measure"),
      lapply(proposals, function(proposal) {
        proposal$measure_sign %||% rep(1, nrow(proposal$eta))
      }),
      isTRUE(gradient),
      lapply(context$subjects, function(evaluator) evaluator$data_input())
    )
    native$states <- Map(function(state, proposal) {
      state$mode <- proposal$mode
      state
    }, native$states, proposals)
    return(native)
  }
  if (is.null(context$parallel)) {
    states <- evaluate_chunk(context$subjects, proposals)
  } else {
    chunks <- context$parallel$chunks
    proposal_chunks <- lapply(chunks, function(rows) proposals[rows])
    pieces <- parallel::clusterApply(
      context$parallel$cluster, seq_along(chunks),
      function(index, chunks, parameters, gradient) {
        namespace <- asNamespace("LibeRation")
        evaluators <- get(".nm_parallel_worker_state", envir = namespace)()$subjects
        evaluate <- get(
          ".nm_imp_subject_from_proposal", envir = asNamespace("LibeRation")
        )
        lapply(seq_along(evaluators), function(subject) {
          evaluate(
            evaluators[[subject]], parameters, chunks[[index]][[subject]], gradient
          )
        })
      }, chunks = proposal_chunks, parameters = parameters, gradient = gradient
    )
    states <- unlist(pieces, recursive = FALSE)
  }
  value <- sum(vapply(states, `[[`, numeric(1), "value"))
  if (!isTRUE(gradient)) return(list(value = value, states = states))
  gradients <- lapply(states, `[[`, "native_gradient")
  if (any(vapply(gradients, is.null, logical(1)))) {
    return(list(value = value, native_gradient = NULL, states = states))
  }
  list(
    value = value, native_gradient = Reduce(`+`, gradients), states = states
  )
}

.nm_imp_objective <- function(context, parameters, normals,
                              eta_maxit, tolerance) {
  if (is.null(context$parallel)) {
    subject_values <- vapply(seq_len(context$n_subjects), function(subject) {
      .nm_imp_subject_objective(
        context$subjects[[subject]], parameters, normals[[subject]],
        eta_maxit, tolerance
      )
    }, numeric(1))
  } else {
    normal_chunks <- lapply(context$parallel$chunks, function(rows) normals[rows])
    pieces <- parallel::clusterApply(
      context$parallel$cluster, seq_along(context$parallel$chunks),
      function(index, chunks, parameters, eta_maxit, tolerance) {
        objective <- get(".nm_imp_subject_objective", envir = asNamespace("LibeRation"))
        namespace <- asNamespace("LibeRation")
        evaluators <- get(".nm_parallel_worker_state", envir = namespace)()$subjects
        worker_normals <- chunks[[index]]
        vapply(seq_along(evaluators), function(subject) {
          objective(
            evaluators[[subject]], parameters, worker_normals[[subject]],
            eta_maxit, tolerance
          )
        }, numeric(1))
      }, chunks = normal_chunks, parameters = parameters,
      eta_maxit = eta_maxit, tolerance = tolerance
    )
    subject_values <- unlist(pieces, use.names = FALSE)
  }
  sum(subject_values) + .nm_prior_nll(context$model, parameters)
}

.nm_imp_evaluate <- function(context, parameters, normals, eta_maxit, tolerance,
                             gradient = TRUE, initial_eta = NULL) {
  proposals <- .nm_imp_prepare_proposals(
    context, parameters, normals, eta_maxit, tolerance,
    adaptive = TRUE, initial_eta = initial_eta
  )
  evaluated <- .nm_imp_evaluate_fixed(
    context, parameters, proposals, gradient = gradient
  )
  states <- evaluated$states
  value <- evaluated$value +
    .nm_prior_nll(context$model, parameters)
  if (!isTRUE(gradient)) return(list(value = value, states = states))
  if (is.null(evaluated$native_gradient)) {
    return(list(value = value, native_gradient = NULL, states = states))
  }
  full <- evaluated$native_gradient
  n_theta <- length(parameters$theta)
  n_sigma <- length(parameters$sigma)
  n_omega <- length(parameters$omega)
  population_positions <- c(
    seq_len(n_theta), n_theta + context$n_eta + seq_len(n_sigma),
    n_theta + context$n_eta + n_sigma + seq_len(n_omega)
  )
  list(
    value = value,
    native_gradient = as.numeric(full[population_positions]) +
      .nm_prior_nll_native_gradient(context$model, parameters),
    states = states
  )
}

.nm_conditional_state_cache <- function(context, evaluator, mu = NULL) {
  cache <- new.env(parent = emptyenv())
  cache$key <- NULL
  cache$parameters <- NULL
  cache$modes <- NULL
  cache$hits <- 0L
  cache$misses <- 0L
  cache$mu_recentered_starts <- 0L

  mode_matrix <- function(result) {
    if (!context$n_eta) return(matrix(numeric(), context$n_subjects, 0L))
    states <- result$states %||% list()
    if (length(states) != context$n_subjects) return(NULL)
    modes <- do.call(rbind, lapply(states, function(state) {
      mode <- state$mode$par %||% state$mode %||% NULL
      if (is.null(mode)) rep(NA_real_, context$n_eta) else as.numeric(mode)
    }))
    if (!identical(dim(modes), c(context$n_subjects, context$n_eta))) NULL else modes
  }

  evaluate <- function(parameters) {
    key <- c(parameters$theta, parameters$sigma, parameters$omega)
    if (!is.null(cache$key) && identical(cache$key, key)) {
      cache$hits <- cache$hits + 1L
      return(cache$result)
    }
    starts <- NULL
    if (!is.null(cache$parameters) && !is.null(cache$modes) &&
        !is.null(mu) && isTRUE(mu$mapped) && isTRUE(mu$enabled) &&
        identical(dim(cache$modes), c(context$n_subjects, context$n_eta)) &&
        all(is.finite(cache$modes))) {
      starts <- .nm_mu_recenter_eta(mu, cache$parameters, parameters, cache$modes)
      cache$mu_recentered_starts <- cache$mu_recentered_starts + 1L
    }
    cache$result <- evaluator(parameters, starts)
    cache$key <- key
    cache$parameters <- parameters
    cache$modes <- mode_matrix(cache$result)
    cache$misses <- cache$misses + 1L
    cache$result
  }

  telemetry <- function() list(
    hits = cache$hits,
    misses = cache$misses,
    recentered_mode_starts = cache$mu_recentered_starts
  )
  list(evaluate = evaluate, telemetry = telemetry)
}

.nm_est_imp_marginal <- function(context, map, maxit, eta_maxit, tolerance, trace,
                                 n_imp = 200L, seed = 20260713L,
                                 print_every = 0L,
                                 imp_gradient = c("score", "finite_crn"),
                                 optimizer_backend = "auto",
                                 mu_specialization = TRUE) {
  n_imp <- as.integer(n_imp)
  if (n_imp < 5L) .nm_stop("IMP requires `n_imp >= 5`.")
  imp_gradient <- match.arg(imp_gradient)
  normals <- .nm_imp_normals(context, n_imp, seed)
  mu <- .nm_mu_specialization(context, map, enabled = mu_specialization)
  cache <- .nm_conditional_state_cache(
    context,
    function(parameters, starts) .nm_imp_evaluate(
      context, parameters, normals, eta_maxit, tolerance,
      gradient = imp_gradient == "score", initial_eta = starts
    ),
    mu = mu
  )
  evaluate <- cache$evaluate
  objective <- function(parameters) evaluate(parameters)$value
  gradient <- if (imp_gradient == "score") function(parameters) {
    result <- evaluate(parameters)
    as.vector(result$native_gradient %*% map$jacobian(parameters))
  } else NULL
  optimizer <- .nm_outer_optim(
    map, objective, maxit, tolerance, trace, print_every,
    gradient = gradient, optimizer_backend = optimizer_backend
  )
  fallback <- NULL
  refinement_reason <- NULL
  if (imp_gradient == "score") {
    refinement_reason <- if (
      !identical(as.integer(optimizer$convergence), 0L)
    ) {
      "score search did not converge"
    } else if (isTRUE(mu$active) && isTRUE(mu$covariate_design)) {
      paste0(
        "subject-varying MU design requires an exact finite-CRN ",
        "objective refinement"
      )
    } else NULL
  }
  if (!is.null(refinement_reason)) {
    # The normalized importance-score gradient deliberately omits proposal
    # derivatives. It is an efficient search direction but is not the exact
    # derivative of the finite common-random-number objective, so L-BFGS-B can
    # report an abnormal line-search termination after reaching its vicinity.
    # Finish such runs against the exact finite-CRN objective without an
    # analytic gradient. This retains the practical fast path while ensuring
    # that a completed fit has an optimizer convergence result it can defend.
    fallback_map <- map
    # An abnormal line-search endpoint can also be below finite-difference
    # resolution, so that case restarts from the declared model values.
    # A converged score search is a useful warm start for the exact refinement.
    fallback_map$start <- if (
      identical(refinement_reason, "score search did not converge")
    ) map$start else optimizer$par
    fallback <- .nm_outer_optim(
      fallback_map, objective, maxit, tolerance, trace, print_every,
      gradient = NULL, optimizer_backend = "r"
    )
    if (identical(as.integer(fallback$convergence), 0L) &&
        is.finite(fallback$value) &&
        fallback$value <= optimizer$value +
          tolerance * max(abs(optimizer$value), 1)) {
      fallback$score_search <- list(
        convergence = optimizer$convergence, value = optimizer$value,
        par = optimizer$par, backend = optimizer$backend,
        objective_evaluations = optimizer$objective_evaluations,
        gradient_evaluations = optimizer$gradient_evaluations
      )
      fallback$backend <- paste0(
        optimizer$backend, "+", fallback$backend, "-finite-crn-refinement"
      )
      optimizer <- fallback
    }
  }
  parameters <- map$decode(optimizer$par)
  modes <- .nm_subject_modes(
    context, parameters, maxit = eta_maxit, tolerance = tolerance,
    exact_hessian = FALSE
  )
  .nm_fit_result(
    context, "IMP", parameters, optimizer$value, modes, optimizer,
    diagnostics = list(
      n_imp = n_imp, seed = seed, eta_maxit = eta_maxit,
      common_random_numbers = TRUE, imp_gradient = imp_gradient,
      mu_specialization = c(
        .nm_mu_diagnostic(mu),
        list(recentered_mode_starts = cache$telemetry()$recentered_mode_starts)
      ),
      conditional_state_cache = cache$telemetry(),
      finite_crn_fallback = !is.null(optimizer$score_search),
      finite_crn_refinement_reason = if (!is.null(optimizer$score_search)) {
        refinement_reason
      } else NULL,
      estimator_identity = "direct finite-sample importance marginal maximum likelihood",
      population_gradient = if (imp_gradient == "score") {
        paste0(
          "normalized importance-score CppAD gradient (proposal derivative ",
          "omitted)",
          if (!is.null(optimizer$score_search)) {
            " with exact finite-CRN derivative-free convergence fallback"
          } else ""
        )
      } else "finite common-random-number objective"
    )
  )
}

.nm_imp_expectation_state <- function(context, parameters, normals, eta_maxit,
                                      tolerance, starts = NULL,
                                      cached_modes = NULL,
                                      proposal_curvature = c("exact", "fisher"),
                                      native_context = NULL) {
  proposal_curvature <- match.arg(proposal_curvature)
  proposals <- .nm_imp_prepare_proposals(
    context, parameters, normals, eta_maxit, tolerance,
    adaptive = TRUE, initial_eta = starts, cached_modes = cached_modes,
    proposal_curvature = proposal_curvature
  )
  direct_native_error <- NULL
  if (!is.null(native_context) &&
      all(vapply(proposals, function(proposal) {
        isTRUE(proposal$valid)
      }, logical(1))) &&
      isTRUE(getOption("LibeRation.imp_native_expectation_install", TRUE))) {
    installed <- tryCatch(
      .liberation_weighted_eta_context_set_importance(
        native_context, parameters$theta, parameters$sigma, parameters$omega,
        lapply(proposals, `[[`, "eta"),
        lapply(proposals, `[[`, "log_proposal")
      ),
      error = identity
    )
    if (!inherits(installed, "error")) {
      return(list(
        eta = NULL, weights = NULL, ess = as.numeric(installed$ess),
        modes = lapply(proposals, `[[`, "mode"),
        backend = installed$backend,
        native_error = NULL, reused_modes = !is.null(cached_modes),
        native_context_installed = TRUE,
        support_points = as.numeric(installed$support_points)
      ))
    }
    direct_native_error <- conditionMessage(installed)
  }
  eta <- vector("list", context$n_subjects)
  weights <- vector("list", context$n_subjects)
  ess <- numeric(context$n_subjects)
  modes <- vector("list", context$n_subjects)
  native_states <- NULL
  native_error <- NULL
  if (is.null(context$parallel) && !isTRUE(context$model$USE_ODE) &&
      all(vapply(proposals, function(proposal) isTRUE(proposal$valid), logical(1))) &&
      isTRUE(getOption("LibeRation.imp_native_expectation", TRUE))) {
    points <- cbind(
      matrix(parameters$theta, context$n_subjects, length(parameters$theta),
             byrow = TRUE),
      matrix(0, context$n_subjects, context$n_eta),
      matrix(parameters$sigma, context$n_subjects, length(parameters$sigma),
             byrow = TRUE),
      matrix(parameters$omega, context$n_subjects, length(parameters$omega),
             byrow = TRUE)
    )
    native <- tryCatch(
      .liberation_objective_tape_importance_collection(
        lapply(context$subjects, function(evaluator) {
          evaluator$objective_tape$pointer
        }),
        points, length(parameters$theta) + seq_len(context$n_eta),
        lapply(proposals, `[[`, "eta"),
        lapply(proposals, `[[`, "log_proposal"),
        lapply(proposals, function(proposal) rep(0, nrow(proposal$eta))),
        lapply(proposals, function(proposal) rep(1, nrow(proposal$eta))),
        FALSE,
        lapply(context$subjects, function(evaluator) evaluator$data_input())
      ),
      error = identity
    )
    if (inherits(native, "error")) {
      native_error <- conditionMessage(native)
    } else {
      candidate_states <- native$states
      valid_states <- length(candidate_states) == context$n_subjects &&
        all(vapply(seq_len(context$n_subjects), function(subject) {
          state <- candidate_states[[subject]]
          isTRUE(state$quadrature_valid) &&
            length(state$weights) == nrow(proposals[[subject]]$eta) &&
            all(is.finite(state$weights)) && all(state$weights >= 0)
        }, logical(1)))
      if (valid_states) native_states <- candidate_states else
        native_error <- "native importance weights were invalid"
    }
  }
  for (subject in seq_len(context$n_subjects)) {
    proposal <- proposals[[subject]]
    if (!isTRUE(proposal$valid)) {
      .nm_stop("IMP conditional proposal failed for subject ", subject, ".")
    }
    probability <- if (!is.null(native_states)) {
      as.numeric(native_states[[subject]]$weights)
    } else {
      joint <- context$subjects[[subject]]$objective_eta_values(
        parameters$theta, proposal$eta, parameters$sigma, parameters$omega
      )
      log_weight <- -0.5 * joint - proposal$log_proposal
      normalizer <- .nm_log_sum_exp(log_weight)
      value <- exp(log_weight - normalizer)
      value / sum(value)
    }
    eta[[subject]] <- proposal$eta
    weights[[subject]] <- probability
    ess[[subject]] <- 1 / sum(probability^2)
    modes[[subject]] <- proposal$mode
  }
  list(
    eta = eta, weights = weights, ess = ess, modes = modes,
    backend = if (!is.null(native_states)) {
      "cpp-batched-importance-weights"
    } else "r-subject-importance-weights",
    native_error = paste(
      Filter(nzchar, c(direct_native_error %||% "", native_error %||% "")),
      collapse = "; "
    ),
    reused_modes = !is.null(cached_modes),
    native_context_installed = FALSE
  )
}

.nm_est_imp <- function(context, map, maxit, eta_maxit, tolerance, trace,
                        n_imp = 200L, seed = 20260713L, print_every = 0L,
                        imp_gradient = c("score", "finite_crn"),
                        optimizer_backend = "auto", mu_specialization = TRUE,
                        imp_algorithm = c("auto", "mcem", "marginal_ml"),
                        imp_mstep_maxit = NULL,
                        imp_mstep_schedule = c("auto", "fixed", "progressive"),
                        imp_sample_schedule = c("auto", "fixed", "progressive"),
                        imp_min_samples = NULL,
                        imp_sampling = c("auto", "random", "antithetic", "rqmc"),
                        imp_proposal = c("auto", "gaussian", "student_t", "defensive"),
                        imp_proposal_curvature = c("auto", "exact", "fisher"),
                        imp_proposal_df = 7,
                        imp_reuse_modes = NULL,
                        imp_mode_refresh_threshold = 0.025,
                        imp_mode_reuse_ess = 0.25,
                        imp_subject_allocation = c("auto", "fixed", "ess"),
                        imp_auto_stop = NULL,
                        imp_stationarity_window = 12L,
                        imp_stationarity_tolerance = 2e-3,
                        imp_stationarity_consecutive = 3L) {
  imp_algorithm <- match.arg(imp_algorithm)
  imp_mstep_schedule <- match.arg(imp_mstep_schedule)
  imp_sample_schedule <- match.arg(imp_sample_schedule)
  imp_sampling <- match.arg(imp_sampling)
  imp_proposal <- match.arg(imp_proposal)
  imp_proposal_curvature <- match.arg(imp_proposal_curvature)
  imp_subject_allocation <- match.arg(imp_subject_allocation)
  resolved_algorithm <- if (imp_algorithm == "auto") "mcem" else imp_algorithm
  if (resolved_algorithm == "marginal_ml") {
    fit <- .nm_est_imp_marginal(
      context, map, maxit, eta_maxit, tolerance, trace, n_imp, seed,
      print_every, imp_gradient, optimizer_backend, mu_specialization
    )
    fit$diagnostics$algorithm_requested <- imp_algorithm
    fit$diagnostics$algorithm_resolved <- resolved_algorithm
    return(fit)
  }
  n_imp <- as.integer(n_imp)
  n_iter <- max(1L, as.integer(maxit))
  if (n_imp < 5L) .nm_stop("IMP requires `n_imp >= 5`.")
  optimized <- .nm_liber_optimized(context)
  resolved_mstep_schedule <- if (imp_mstep_schedule == "auto") {
    if (optimized) "progressive" else "fixed"
  } else imp_mstep_schedule
  resolved_sample_schedule <- if (imp_sample_schedule == "auto") {
    if (optimized) "progressive" else "fixed"
  } else imp_sample_schedule
  resolved_sampling <- if (imp_sampling == "auto") {
    # Antithetic draws remain the optimized automatic policy: standard-profile
    # measurement showed lower generation cost and a smoother MCEM Q surface
    # than independently shifted RQMC at the same nominal draw budget. RQMC is
    # retained as an explicit variance-reduction option for accuracy studies.
    if (optimized) "antithetic" else "random"
  } else imp_sampling
  resolved_proposal <- if (imp_proposal == "auto") {
    if (optimized) "defensive" else "gaussian"
  } else imp_proposal
  resolved_proposal_curvature <- if (imp_proposal_curvature == "auto") {
    if (optimized) "fisher" else "exact"
  } else imp_proposal_curvature
  resolved_subject_allocation <- if (imp_subject_allocation == "auto") {
    # Reallocating a fixed total draw budget by the previous iteration's ESS
    # is valid importance sampling, but changes the finite MCEM surface enough
    # to increase optimizer evaluations on the standard benchmark. Keep it as
    # an explicit difficult-subject option; fixed allocation is the faster and
    # more stable automatic policy in both numerical modes.
    "fixed"
  } else imp_subject_allocation
  imp_min_samples <- as.integer(imp_min_samples %||%
    if (optimized) max(10L, min(n_imp, as.integer(ceiling(n_imp / 4)))) else n_imp)
  imp_reuse_modes <- isTRUE(imp_reuse_modes %||% optimized)
  imp_auto_stop <- isTRUE(imp_auto_stop %||% optimized)
  imp_stationarity_window <- as.integer(imp_stationarity_window)
  imp_stationarity_consecutive <- as.integer(imp_stationarity_consecutive)
  if (!optimized) {
    resolved_mstep_schedule <- "fixed"
    resolved_sample_schedule <- "fixed"
    resolved_sampling <- "random"
    resolved_proposal <- "gaussian"
    resolved_proposal_curvature <- "exact"
    imp_min_samples <- n_imp
    imp_reuse_modes <- FALSE
    imp_auto_stop <- FALSE
    resolved_subject_allocation <- "fixed"
  }
  if (is.na(imp_min_samples) || imp_min_samples < 5L ||
      imp_min_samples > n_imp || !is.finite(imp_proposal_df) ||
      imp_proposal_df <= 2 || !is.finite(imp_mode_refresh_threshold) ||
      imp_mode_refresh_threshold <= 0 || !is.finite(imp_mode_reuse_ess) ||
      imp_mode_reuse_ess <= 0 || imp_mode_reuse_ess >= 1) {
    .nm_stop("Adaptive IMP controls are invalid.")
  }
  if (is.na(imp_stationarity_window) || imp_stationarity_window < 4L ||
      !is.finite(imp_stationarity_tolerance) ||
      imp_stationarity_tolerance <= 0 ||
      is.na(imp_stationarity_consecutive) ||
      imp_stationarity_consecutive < 1L) {
    .nm_stop("IMP stationarity controls are invalid.")
  }
  if (is.null(imp_mstep_maxit)) {
    imp_mstep_maxit <- if (.nm_liber_optimized(context)) 10L else 1L
  }
  imp_mstep_maxit <- as.integer(imp_mstep_maxit)
  if (length(imp_mstep_maxit) != 1L || is.na(imp_mstep_maxit) ||
      imp_mstep_maxit < 1L) {
    .nm_stop("`imp_mstep_maxit` must be one positive integer.")
  }
  parameters <- map$decode(map$start)
  mu <- .nm_mu_specialization(context, map, enabled = mu_specialization)
  mu_recentered_starts <- 0L
  starts <- matrix(0, context$n_subjects, context$n_eta)
  objective_trace <- numeric(n_iter)
  ess_trace <- matrix(NA_real_, n_iter, context$n_subjects)
  parameter_trace <- if (length(map$start)) {
    matrix(NA_real_, n_iter, length(map$start))
  } else matrix(numeric(), n_iter, 0L)
  total_evaluations <- total_gradient_evaluations <- total_mstep_iterations <- 0L
  optimizer <- NULL
  completed <- 0L
  expectation_state <- NULL
  phases <- .nm_stochastic_phase_timer()
  weighted_context <- .nm_weighted_eta_context(
    context,
    reduced_population_tape = resolved_sample_schedule == "fixed" &&
      resolved_subject_allocation == "fixed"
  )
  imp_priors <- context$model$LIK_CONFIG$priors
  has_sigma_prior <- !is.null(imp_priors) && nrow(imp_priors) &&
    any(startsWith(toupper(imp_priors$parameter), "SIGMA"))
  has_omega_prior <- !is.null(imp_priors) && nrow(imp_priors) &&
    any(startsWith(toupper(imp_priors$parameter), "OMEGA"))
  imp_simple_sigma <- !has_sigma_prior &&
    context$model$LIK_CONFIG$error %in%
      c("additive", "proportional", "exponential") &&
    identical(context$model$LIK_CONFIG$sigma_corr %||% "independent",
              "independent") &&
    !length(context$model$LIK_CONFIG$residual_groups)
  imp_native_model <- context$model
  imp_native_model$THETAS$Value <- parameters$theta
  imp_native_model$SIGMAS$Value <- parameters$sigma
  imp_native_model$OMEGAS$Value <- parameters$omega
  if (length(map$omega_free)) imp_native_model$OMEGAS$FIX[] <- TRUE
  if (imp_simple_sigma && length(map$sigma_free)) {
    imp_native_model$SIGMAS$FIX[] <- TRUE
  }
  imp_native_map <- .nm_outer_map(imp_native_model)
  imp_native_enabled <- isTRUE(getOption(
    "LibeRation.imp_native_mstep", TRUE
  )) && optimized && !is.null(weighted_context) &&
    is.null(context$parallel) &&
    optimizer_backend %in% c("auto", "native") &&
    !has_omega_prior
  imp_native_state <- NULL
  imp_native_attempts <- 0L
  imp_native_successes <- 0L
  imp_native_fallbacks <- 0L
  imp_native_fallback_reason <- NULL
  imp_sigma_gradient_updates <- 0L
  imp_omega_moment_updates <- 0L
  sample_trace <- integer(n_iter)
  subject_sample_trace <- matrix(NA_integer_, n_iter, context$n_subjects)
  mstep_effort_trace <- integer(n_iter)
  proposal_reuse_trace <- logical(n_iter)
  cached_modes <- NULL
  proposal_anchor <- NULL
  previous_relative_ess <- 0
  previous_subject_relative_ess <- rep(NA_real_, context$n_subjects)
  stationary_iterations <- 0L
  stationarity <- .nm_saem_stationarity(
    numeric(), parameter_trace[FALSE, , drop = FALSE], 0L,
    imp_stationarity_window, imp_stationarity_tolerance
  )
  for (iteration in seq_len(n_iter)) {
    # Advancing the deterministic seed by iteration gives an independent,
    # reproducible Monte-Carlo E-step rather than silently reusing one finite
    # common-random-number objective as a surrogate for MCEM.
    scheduled_samples <- if (resolved_sample_schedule == "progressive") {
      imp_min_samples + as.integer(ceiling(
        (n_imp - imp_min_samples) * (iteration / n_iter)^1.5
      ))
    } else n_imp
    if (previous_relative_ess < imp_mode_reuse_ess / 2 && iteration > 1L) {
      scheduled_samples <- max(scheduled_samples, min(n_imp,
        max(imp_min_samples, 2L * sample_trace[[iteration - 1L]])))
    }
    scheduled_samples <- min(n_imp, max(5L, scheduled_samples))
    subject_samples <- rep.int(scheduled_samples, context$n_subjects)
    if (resolved_subject_allocation == "ess" && iteration > 1L &&
        all(is.finite(previous_subject_relative_ess))) {
      difficulty <- 1 / sqrt(pmax(previous_subject_relative_ess, 0.05))
      proposed <- as.integer(round(
        scheduled_samples * context$n_subjects * difficulty / sum(difficulty)
      ))
      subject_samples <- pmin(n_imp, pmax(5L, proposed))
    }
    sample_trace[[iteration]] <- as.integer(round(mean(subject_samples)))
    subject_sample_trace[iteration, ] <- subject_samples
    normals <- phases$time("normal_generation", {
      .nm_imp_normals(
        context, subject_samples, seed + iteration - 1L,
        sampling = resolved_sampling, proposal = resolved_proposal,
        proposal_df = imp_proposal_df
      )
    })
    current_anchor <- c(parameters$theta, parameters$sigma, parameters$omega)
    anchor_drift <- if (is.null(proposal_anchor)) Inf else max(
      abs(current_anchor - proposal_anchor) / (1 + abs(proposal_anchor))
    )
    reusable_modes <- optimized && imp_reuse_modes &&
      !is.null(cached_modes) && is.finite(anchor_drift) &&
      anchor_drift <= imp_mode_refresh_threshold &&
      previous_relative_ess >= imp_mode_reuse_ess
    expectation_state <- phases$time("expectation", {
      .nm_imp_expectation_state(
        context, parameters, normals, eta_maxit, tolerance, starts,
        cached_modes = if (reusable_modes) cached_modes else NULL,
        proposal_curvature = resolved_proposal_curvature,
        native_context = if (optimized) weighted_context else NULL
      )
    })
    proposal_reuse_trace[[iteration]] <- reusable_modes
    if (!reusable_modes) {
      cached_modes <- expectation_state$modes
      proposal_anchor <- current_anchor
    }
    if (context$n_eta) {
      starts <- do.call(rbind, lapply(expectation_state$modes, `[[`, "par"))
    }
    expectation <- phases$time("expectation_setup", {
      if (isTRUE(expectation_state$native_context_installed)) {
        .nm_native_weighted_expectation(context, map, weighted_context)
      } else {
        .nm_complete_data_expectation(
          context, map, expectation_state$eta, expectation_state$weights,
          native_context = weighted_context
        )
      }
    })
    iteration_map <- map
    iteration_map$start <- map$encode(parameters)
    current_mstep_maxit <- if (resolved_mstep_schedule == "progressive") {
      min(imp_mstep_maxit, max(1L, as.integer(ceiling(
        imp_mstep_maxit * iteration / n_iter
      ))))
    } else imp_mstep_maxit
    mstep_effort_trace[[iteration]] <- current_mstep_maxit
    native_result <- NULL
    if (imp_native_enabled && length(imp_native_map$start)) {
      imp_native_attempts <- imp_native_attempts + 1L
      native_result <- phases$time("mstep", tryCatch(
        .nm_native_weighted_mstep(
          context, parameters, weighted_context, imp_native_map,
          current_mstep_maxit, tolerance,
          if (trace > 1L) trace else 0L, imp_native_state
        ),
        error = identity
      ))
      if (inherits(native_result, "error") ||
          !is.finite(native_result$value %||% NA_real_) ||
          any(!is.finite(native_result$theta %||% NA_real_)) ||
          any(!is.finite(native_result$sigma %||% NA_real_))) {
        imp_native_fallbacks <- imp_native_fallbacks + 1L
        imp_native_fallback_reason <- if (inherits(native_result, "error")) {
          conditionMessage(native_result)
        } else "native IMP M-step returned non-finite values"
        native_result <- NULL
        imp_native_enabled <- FALSE
      } else {
        imp_native_successes <- imp_native_successes + 1L
        imp_native_state <- native_result$optimizer_state
        native_result$optimizer_state <- NULL
      }
    }
    optimizer <- if (!is.null(native_result)) {
      native_result
    } else phases$time("mstep", {
      .nm_outer_optim(
        iteration_map, expectation$objective, current_mstep_maxit, tolerance,
        if (trace > 1L) trace else 0L, 0L,
        gradient = expectation$gradient, optimizer_backend = optimizer_backend
      )
    })
    previous_parameters <- parameters
    previous <- map$encode(previous_parameters)
    parameters <- if (!is.null(native_result)) {
      list(
        theta = as.numeric(native_result$theta),
        sigma = as.numeric(native_result$sigma),
        omega = as.numeric(native_result$omega)
      )
    } else map$decode(optimizer$par)
    if (!is.null(native_result) && imp_simple_sigma && length(map$sigma_free)) {
      sigma_update <- .nm_saem_sigma_from_q_gradient(
        context, parameters, native_result$native_gradient
      )
      if (!is.null(sigma_update)) {
        parameters$sigma[map$sigma_free] <- sigma_update[map$sigma_free]
        imp_sigma_gradient_updates <- imp_sigma_gradient_updates + 1L
      }
    }
    if (!is.null(native_result) && length(map$omega_free) && context$n_eta) {
      omega_update <- .liberation_weighted_eta_context_omega(
        weighted_context, as.integer(context$model$n_eta),
        as.integer(context$model$LIK_CONFIG$iov),
        as.integer(context$model$OMEGAS$ROW),
        as.integer(context$model$OMEGAS$COL)
      )
      parameters$omega[map$omega_free] <- omega_update[map$omega_free]
      imp_omega_moment_updates <- imp_omega_moment_updates + 1L
    }
    if (context$n_eta && isTRUE(mu$enabled) && isTRUE(mu$mapped)) {
      starts <- .nm_mu_recenter_eta(
        mu, previous_parameters, parameters, starts
      )
      mu_recentered_starts <- mu_recentered_starts + 1L
    }
    point <- map$encode(parameters)
    objective_trace[[iteration]] <- expectation$objective(parameters)
    optimizer$par <- point
    optimizer$value <- objective_trace[[iteration]]
    ess_trace[iteration, ] <- expectation_state$ess
    previous_subject_relative_ess <- expectation_state$ess / subject_samples
    previous_relative_ess <- mean(previous_subject_relative_ess)
    if (length(point)) parameter_trace[iteration, ] <- point
    total_evaluations <- total_evaluations +
      as.integer(optimizer$objective_evaluations %||% 0L)
    total_gradient_evaluations <- total_gradient_evaluations +
      as.integer(optimizer$gradient_evaluations %||% 0L)
    total_mstep_iterations <- total_mstep_iterations +
      as.integer(optimizer$iterations %||% 0L)
    completed <- iteration
    if (print_every > 0L && iteration %% print_every == 0L) {
      cat(sprintf(
        "[LibeRation] IMP MCEM ITERATION %d Q %.10g ESS %.1f\n",
        iteration, objective_trace[[iteration]], mean(expectation_state$ess)
      ))
      try(flush(stdout()), silent = TRUE)
    }
    stationarity <- .nm_saem_stationarity(
      objective_trace[seq_len(iteration)],
      parameter_trace[seq_len(iteration), , drop = FALSE], 0L,
      imp_stationarity_window, imp_stationarity_tolerance
    )
    stationary_iterations <- if (
      isTRUE(stationarity$converged) && previous_relative_ess >= 0.1
    ) stationary_iterations + 1L else 0L
    if (optimized && imp_auto_stop &&
        stationary_iterations >= imp_stationarity_consecutive) break
    if (!optimized && length(point) && iteration > 1L &&
        max(abs(point - previous) / (1 + abs(previous))) <= tolerance) break
  }
  final_normals <- .nm_imp_normals(
    context, n_imp, seed + completed, sampling = resolved_sampling,
    proposal = resolved_proposal, proposal_df = imp_proposal_df
  )
  final <- .nm_imp_evaluate(
    context, parameters, final_normals, eta_maxit, tolerance,
    gradient = FALSE, initial_eta = starts
  )
  modes <- .nm_subject_modes(
    context, parameters, starts = starts, maxit = eta_maxit,
    tolerance = tolerance, exact_hessian = FALSE
  )
  parameter_converged <- completed < n_iter
  optimizer$convergence <- 0L
  optimizer$message <- if (completed < n_iter) "IMP MCEM parameter convergence reached" else
    "IMP MCEM iterations completed"
  optimizer$iterations <- completed
  optimizer$objective_evaluations <- total_evaluations
  optimizer$gradient_evaluations <- total_gradient_evaluations
  optimizer$mstep_iterations <- total_mstep_iterations
  .nm_fit_result(
    context, "IMP", parameters, final$value, modes, optimizer,
    diagnostics = list(
      n_imp = n_imp, seed = seed, eta_maxit = eta_maxit,
      algorithm_requested = imp_algorithm,
      algorithm_resolved = resolved_algorithm,
      estimator_identity = "importance-sampling Monte-Carlo EM",
      mstep_maxit = imp_mstep_maxit,
      mstep_schedule = resolved_mstep_schedule,
      mstep_effort = mstep_effort_trace[seq_len(completed)],
      sample_schedule = resolved_sample_schedule,
      sample_count = sample_trace[seq_len(completed)],
      subject_sample_allocation = list(
        strategy = resolved_subject_allocation,
        count = subject_sample_trace[seq_len(completed), , drop = FALSE]
      ),
      sampling = resolved_sampling,
      proposal = resolved_proposal,
      proposal_curvature = resolved_proposal_curvature,
      proposal_df = imp_proposal_df,
      proposal_mode_reuse = list(
        enabled = imp_reuse_modes,
        refresh_threshold = imp_mode_refresh_threshold,
        minimum_relative_ess = imp_mode_reuse_ess,
        reused_iterations = sum(proposal_reuse_trace[seq_len(completed)]),
        trace = proposal_reuse_trace[seq_len(completed)]
      ),
      stationarity = c(stationarity, list(
        auto_stop = imp_auto_stop,
        stopped_early = completed < n_iter,
        consecutive_confirmations = stationary_iterations,
        required_confirmations = imp_stationarity_consecutive
      )),
      expectation_backend = expectation_state$backend,
      expectation_native_error = expectation_state$native_error,
      weighted_expectation = expectation$telemetry(),
      native_weighted_mstep = list(
        eligible = !has_omega_prior && optimized &&
          is.null(context$parallel),
        enabled = imp_native_enabled,
        attempts = imp_native_attempts,
        successes = imp_native_successes,
        fallbacks = imp_native_fallbacks,
        fallback_reason = imp_native_fallback_reason,
        persistent_optimizer_state = !is.null(imp_native_state),
        sigma_gradient_updates = imp_sigma_gradient_updates,
        omega_moment_updates = imp_omega_moment_updates
      ),
      phase_timing = phases$snapshot(),
      parameter_converged = parameter_converged,
      iteration_limit_reached = !parameter_converged,
      independent_e_steps = TRUE,
      effective_sample_size = ess_trace[seq_len(completed), , drop = FALSE],
      objective_trace = objective_trace[seq_len(completed)],
      parameter_trace = parameter_trace[seq_len(completed), , drop = FALSE],
      imp_gradient = "complete-data-expectation",
      common_random_numbers = FALSE,
      finite_crn_fallback = FALSE,
      mu_specialization = c(
        .nm_mu_diagnostic(mu),
        list(recentered_mode_starts = mu_recentered_starts)
      ),
      population_gradient = "exact CppAD gradient of the fixed Monte-Carlo E-step"
    )
  )
}

.nm_gq_evaluate <- function(context, parameters, normals, eta_maxit, tolerance,
                            adaptive = TRUE, gradient = TRUE,
                            initial_eta = NULL) {
  proposals <- .nm_imp_prepare_proposals(
    context, parameters, normals, eta_maxit, tolerance, adaptive = adaptive,
    initial_eta = initial_eta
  )
  evaluated <- .nm_imp_evaluate_fixed(
    context, parameters, proposals, gradient = gradient
  )
  value <- evaluated$value + .nm_prior_nll(context$model, parameters)
  if (!isTRUE(gradient) || is.null(evaluated$native_gradient)) {
    return(list(value = value, native_gradient = NULL, states = evaluated$states,
                proposals = proposals))
  }
  n_theta <- length(parameters$theta)
  n_sigma <- length(parameters$sigma)
  n_omega <- length(parameters$omega)
  population_positions <- c(
    seq_len(n_theta), n_theta + context$n_eta + seq_len(n_sigma),
    n_theta + context$n_eta + n_sigma + seq_len(n_omega)
  )
  list(
    value = value,
    native_gradient = as.numeric(evaluated$native_gradient[population_positions]) +
      .nm_prior_nll_native_gradient(context$model, parameters),
    states = evaluated$states, proposals = proposals
  )
}

.nm_est_gq <- function(context, map, maxit, eta_maxit, tolerance, trace,
                       gq_order = 5L, gq_adaptive = TRUE,
                       gq_max_points = 100000L,
                       gq_grid = c("auto", "tensor", "smolyak"),
                       gq_level = 3L,
                       gq_gradient = c("auto", "score", "finite_grid"),
                       print_every = 0L, optimizer_backend = "auto",
                       mu_specialization = TRUE) {
  gq_order <- as.integer(gq_order)
  gq_level <- as.integer(gq_level)
  gq_max_points <- as.integer(gq_max_points)
  if (length(gq_order) != 1L || is.na(gq_order) || gq_order < 1L) {
    .nm_stop("`gq_order` must be one positive integer.")
  }
  if (length(gq_max_points) != 1L || is.na(gq_max_points) ||
      gq_max_points < 1L) {
    .nm_stop("`gq_max_points` must be one positive integer.")
  }
  if (length(gq_level) != 1L || is.na(gq_level) || gq_level < 1L) {
    .nm_stop("`gq_level` must be one positive integer.")
  }
  if (length(gq_adaptive) != 1L || is.na(gq_adaptive)) {
    .nm_stop("`gq_adaptive` must be TRUE or FALSE.")
  }
  gq_adaptive <- isTRUE(gq_adaptive)
  if (length(gq_grid) > 1L) gq_grid <- gq_grid[[1L]]
  if (length(gq_grid) != 1L || is.na(gq_grid)) {
    .nm_stop("`gq_grid` must contain one grid strategy.")
  }
  gq_grid <- tolower(as.character(gq_grid))
  if (identical(gq_grid, "sparse")) gq_grid <- "smolyak"
  if (!gq_grid %in% c("auto", "tensor", "smolyak")) {
    .nm_stop("`gq_grid` must be one of auto, tensor, or smolyak.")
  }
  gq_gradient <- match.arg(gq_gradient)
  resolved_gradient <- if (gq_gradient == "auto") {
    if (.nm_liber_optimized(context)) "score" else "finite_grid"
  } else gq_gradient
  design <- .nm_gq_design(
    context, order = gq_order, max_points = gq_max_points,
    grid = gq_grid, level = gq_level
  )
  mu <- .nm_mu_specialization(context, map, enabled = mu_specialization)
  native_status <- list(
    eligible = .nm_liber_optimized(context) &&
      identical(resolved_gradient, "score") &&
      !identical(optimizer_backend, "r") &&
      is.null(context$parallel) && !isTRUE(mu$mapped) &&
      isTRUE(getOption("LibeRation.gq_native_coordinator", TRUE)),
    used = FALSE, fallback_reason = NULL
  )
  if (isTRUE(native_status$eligible)) {
    started <- proc.time()[["elapsed"]]
    native <- tryCatch({
      stochastic <- .nm_stochastic_eta_context(context)
      if (is.null(stochastic)) {
        stop("persistent stochastic subject context was unavailable",
             call. = FALSE)
      }
      nodes <- design$normals[[1L]]
      pointer <- .liberation_gq_context_create(
        stochastic, .nm_bayes_cpp_map_config(context, map, mu), nodes,
        as.numeric(attr(nodes, "log_measure", exact = TRUE)),
        as.numeric(attr(nodes, "measure_sign", exact = TRUE)),
        gq_adaptive, as.integer(eta_maxit), as.numeric(tolerance)
      )
      result <- .liberation_gq_context_optimize(
        pointer, as.integer(maxit), as.integer(trace), TRUE
      )
      result$pointer <- pointer
      result
    }, error = identity)
    if (!inherits(native, "error")) {
      native_status$used <- TRUE
      optimizer <- native$optimizer
      optimizer$elapsed_seconds <- unname(
        proc.time()[["elapsed"]] - started
      )
      optimizer$objective_initialization_seconds <- 0
      parameters <- map$decode(optimizer$par)
      modes <- .nm_subject_modes(
        context, parameters, maxit = eta_maxit, tolerance = tolerance,
        exact_hessian = FALSE
      )
      return(.nm_fit_result(
        context, "GQ", parameters, optimizer$value, modes, optimizer,
        diagnostics = list(
          quadrature_order = design$quadrature_order,
          quadrature_level = design$quadrature_level,
          quadrature_points = design$actual_samples,
          quadrature_candidate_points = design$candidate_points,
          quadrature_max_points = design$max_points,
          quadrature_grid_requested = design$requested_grid,
          quadrature_grid = design$resolved_grid,
          quadrature_negative_weights = design$negative_weights,
          adaptive = gq_adaptive, gq_gradient = gq_gradient,
          gq_gradient_resolved = resolved_gradient,
          estimator_identity = paste0(
            if (gq_adaptive) "adaptive " else "fixed ",
            design$resolved_grid, " Gauss-Hermite quadrature"
          ),
          exact_finite_grid_refinement = isTRUE(
            native$exact_finite_grid_refinement
          ),
          mu_specialization = c(
            .nm_mu_diagnostic(mu), list(recentered_mode_starts = 0L)
          ),
          conditional_state_cache = list(
            hits = native$telemetry$cache_hits %||% 0L,
            misses = native$telemetry$parameter_evaluations %||% 0L,
            recentered_mode_starts = 0L
          ),
          native_gq_coordinator = c(native_status, list(
            telemetry = native$telemetry
          )),
          effective_quadrature_points = as.numeric(
            native$effective_quadrature_points
          ),
          quadrature_cancellation_ratio = as.numeric(
            native$quadrature_cancellation_ratio
          ),
          population_gradient = paste0(
            "normalized quadrature-score CppAD search direction (node ",
            "derivative omitted); final convergence against the complete ",
            "finite quadrature-grid objective in the persistent C++ coordinator"
          )
        )
      ))
    }
    native_status$fallback_reason <- conditionMessage(native)
  } else {
    native_status$fallback_reason <- if (!.nm_liber_optimized(context)) {
      "NONMEM-compatible policy retains the established R/L-BFGS-B coordinator"
    } else if (!identical(resolved_gradient, "score")) {
      "native coordination currently requires the score-search policy"
    } else if (identical(optimizer_backend, "r")) {
      "the R optimizer backend was requested explicitly"
    } else if (!is.null(context$parallel)) {
      "a PSOCK subject cluster is active"
    } else if (isTRUE(mu$mapped)) {
      "MU-aware proposal recentering remains on the established coordinator"
    } else "native GQ coordination was disabled"
  }
  cache <- .nm_conditional_state_cache(
    context,
    function(parameters, starts) .nm_gq_evaluate(
      context, parameters, design$normals, eta_maxit, tolerance,
      adaptive = gq_adaptive, gradient = resolved_gradient == "score",
      initial_eta = starts
    ),
    mu = mu
  )
  evaluate <- cache$evaluate
  objective <- function(parameters) evaluate(parameters)$value
  gradient <- if (resolved_gradient == "score") function(parameters) {
    result <- evaluate(parameters)
    if (is.null(result$native_gradient)) return(NULL)
    as.vector(result$native_gradient %*% map$jacobian(parameters))
  } else NULL
  optimizer <- .nm_outer_optim(
    map, objective, maxit, tolerance, trace, print_every,
    gradient = gradient, optimizer_backend = optimizer_backend
  )
  score_search <- NULL
  if (resolved_gradient == "score" && length(map$start)) {
    # Adaptive-node score derivatives are an excellent search direction but
    # omit movement of the conditional mode and quadrature transform. Finish
    # against the complete finite-grid objective so the returned estimate is a
    # stationary point of the objective actually reported by the estimator.
    score_search <- optimizer
    exact_map <- map
    exact_map$start <- optimizer$par
    refined <- .nm_outer_optim(
      exact_map, objective, maxit, tolerance, trace, print_every,
      gradient = NULL, optimizer_backend = "r"
    )
    if (is.finite(refined$value) &&
        refined$value <= optimizer$value + tolerance * max(1, abs(optimizer$value))) {
      refined$score_search <- list(
        value = optimizer$value, convergence = optimizer$convergence,
        objective_evaluations = optimizer$objective_evaluations,
        gradient_evaluations = optimizer$gradient_evaluations
      )
      refined$backend <- paste0(
        optimizer$backend, "+", refined$backend, "-exact-finite-grid-refinement"
      )
      optimizer <- refined
    }
  }
  parameters <- map$decode(optimizer$par)
  final <- evaluate(parameters)
  modes <- .nm_subject_modes(
    context, parameters, maxit = eta_maxit, tolerance = tolerance,
    exact_hessian = FALSE
  )
  effective <- vapply(
    final$states, function(state) state$effective_sample_size %||% NA_real_,
    numeric(1)
  )
  cancellation <- vapply(
    final$states, function(state) state$cancellation_ratio %||% 1,
    numeric(1)
  )
  .nm_fit_result(
    context, "GQ", parameters, optimizer$value, modes, optimizer,
    diagnostics = list(
      quadrature_order = design$quadrature_order,
      quadrature_level = design$quadrature_level,
      quadrature_points = design$actual_samples,
      quadrature_candidate_points = design$candidate_points,
      quadrature_max_points = design$max_points,
      quadrature_grid_requested = design$requested_grid,
      quadrature_grid = design$resolved_grid,
      quadrature_negative_weights = design$negative_weights,
      adaptive = gq_adaptive, gq_gradient = gq_gradient,
      gq_gradient_resolved = resolved_gradient,
      estimator_identity = paste0(
        if (gq_adaptive) "adaptive " else "fixed ",
        design$resolved_grid, " Gauss-Hermite quadrature"
      ),
      exact_finite_grid_refinement = !is.null(optimizer$score_search),
      mu_specialization = c(
        .nm_mu_diagnostic(mu),
        list(recentered_mode_starts = cache$telemetry()$recentered_mode_starts)
      ),
      conditional_state_cache = cache$telemetry(),
      native_gq_coordinator = native_status,
      effective_quadrature_points = effective,
      quadrature_cancellation_ratio = cancellation,
      population_gradient = if (resolved_gradient == "score") {
        paste0(
          "normalized quadrature-score CppAD search direction (node derivative omitted)",
          if (!is.null(optimizer$score_search))
            "; final convergence against complete finite quadrature-grid objective" else ""
        )
      } else {
        "finite quadrature-grid objective"
      }
    )
  )
}

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

.nm_bayes_state <- function(map, outer, eta, subject_values = NULL) {
  list(
    outer = outer, parameters = map$decode(outer), eta = eta,
    subject_values = subject_values
  )
}

.nm_bayes_adaptive_state <- function(dimension, initial_scale,
                                     adapt_start, adapt_interval,
                                     target_acceptance) {
  state <- new.env(parent = emptyenv())
  state$dimension <- as.integer(dimension)
  state$n <- 0L
  state$mean <- numeric(dimension)
  state$m2 <- matrix(0, dimension, dimension)
  state$root <- diag(initial_scale, dimension)
  state$covariance <- diag(initial_scale^2, dimension)
  state$log_multiplier <- 0
  state$adapt_start <- as.integer(adapt_start)
  state$adapt_interval <- as.integer(adapt_interval)
  state$target_acceptance <- as.numeric(target_acceptance)
  state$root_updates <- 0L
  state$regularizations <- 0L
  state$update <- function(value, accepted, adapting = TRUE) {
    value <- as.numeric(value)
    state$n <- state$n + 1L
    delta <- value - state$mean
    state$mean <- state$mean + delta / state$n
    state$m2 <- state$m2 + tcrossprod(delta, value - state$mean)
    if (!isTRUE(adapting)) return(invisible(state$root))
    gain <- min(0.02, (state$n + 10)^(-0.6))
    state$log_multiplier <- state$log_multiplier +
      gain * (as.numeric(isTRUE(accepted)) - state$target_acceptance)
    if (state$n < state$adapt_start ||
        state$n %% state$adapt_interval != 0L) {
      return(invisible(state$root))
    }
    empirical <- state$m2 / max(state$n - 1L, 1L)
    optimal <- 2.38^2 / max(state$dimension, 1L)
    ridge <- initial_scale^2 * 1e-3
    candidate <- exp(2 * state$log_multiplier) *
      (optimal * empirical + diag(ridge, state$dimension))
    repaired <- .nm_positive_definite(
      candidate, "adaptive BAYES population proposal"
    )
    state$regularizations <- state$regularizations +
      as.integer((repaired$jitter %||% 0) > 0)
    state$covariance <- repaired$matrix
    state$root <- t(chol(repaired$matrix))
    state$root_updates <- state$root_updates + 1L
    invisible(state$root)
  }
  state
}

.nm_bayes_cpp_map_config <- function(context, map, mu = NULL) {
  priors <- .nm_cpp_prior_config(context$model)
  list(
    theta = as.numeric(context$model$THETAS$Value),
    sigma = as.numeric(context$model$SIGMAS$Value),
    omega = as.numeric(context$model$OMEGAS$Value),
    theta_free = as.integer(map$theta_free),
    sigma_free = as.integer(map$sigma_free),
    omega_free = as.integer(map$omega_free),
    omega_full = isTRUE(map$omega_full),
    omega_rows = as.integer(context$model$OMEGAS$ROW),
    omega_cols = as.integer(context$model$OMEGAS$COL),
    n_eta_base = as.integer(context$model$n_eta),
    start = as.numeric(map$start), lower = as.numeric(map$lower),
    upper = as.numeric(map$upper),
    prior_index = priors$index, prior_family = priors$family,
    prior_mean = priors$mean, prior_sd = priors$sd,
    prior_shape = priors$shape, prior_rate = priors$rate,
    mu = if (isTRUE(mu$active) && length(mu$theta)) list(
      active = TRUE, theta = as.integer(mu$theta),
      links = unname(mu$links[as.character(mu$theta)]),
      design_columns = unname(mu$design_columns)
    ) else list(active = FALSE)
  )
}

.nm_est_bayes_single <- function(context, map, tolerance,
                          n_burn = 500L, n_sample = 1000L, n_thin = 1L,
                          step_scale = 0.03, eta_step = 0.35,
                          seed = 20260713L, adapt = TRUE,
                          print_every = 0L,
                          mu_specialization = TRUE,
                          outer_kernel = c(
                            "auto", "isotropic", "adaptive_metropolis"
                          ),
                          adaptive_start = 50L,
                          adaptive_interval = 10L,
                          adaptive_target = NULL,
                          delayed_rejection_scale = NULL,
                          eta_kernel = c(
                            "auto", "random_walk", "laplace", "student_t"
                          ),
                          bayes_eta_refresh = 25L,
                          bayes_eta_maxit = 50L,
                          bayes_eta_df = 7,
                          bayes_eta_rescue_probability = 0.05,
                          bayes_eta_parameter_refresh = 0.15,
                          bayes_eta_low_acceptance = 0.1,
                          bayes_gibbs_omega = TRUE) {
  n_burn <- as.integer(n_burn)
  n_sample <- as.integer(n_sample)
  n_thin <- as.integer(n_thin)
  outer_kernel <- match.arg(outer_kernel)
  eta_kernel <- match.arg(eta_kernel)
  adaptive_start <- as.integer(adaptive_start)
  adaptive_interval <- as.integer(adaptive_interval)
  bayes_eta_refresh <- as.integer(bayes_eta_refresh)
  bayes_eta_maxit <- as.integer(bayes_eta_maxit)
  bayes_eta_df <- as.numeric(bayes_eta_df)
  delayed_rejection_scale <- as.numeric(
    delayed_rejection_scale %||% if (.nm_liber_optimized(context)) 0.25 else 0
  )
  if (n_burn < 0L || n_sample < 1L || n_thin < 1L) {
    .nm_stop("BAYES requires n_burn >= 0, n_sample >= 1, and n_thin >= 1.")
  }
  if (is.na(adaptive_start) || adaptive_start < 2L ||
      is.na(adaptive_interval) || adaptive_interval < 1L) {
    .nm_stop("Adaptive BAYES start and interval must be positive integers.")
  }
  if (is.na(bayes_eta_refresh) || bayes_eta_refresh < 1L ||
      is.na(bayes_eta_maxit) || bayes_eta_maxit < 1L ||
      !is.finite(bayes_eta_rescue_probability) ||
      bayes_eta_rescue_probability < 0 || bayes_eta_rescue_probability >= 1 ||
      !is.finite(bayes_eta_parameter_refresh) ||
      bayes_eta_parameter_refresh <= 0 ||
      length(bayes_eta_df) != 1L || !is.finite(bayes_eta_df) ||
      bayes_eta_df <= 2 ||
      !is.finite(bayes_eta_low_acceptance) || bayes_eta_low_acceptance < 0 ||
      bayes_eta_low_acceptance >= 1) {
    .nm_stop("Optimized BAYES ETA proposal controls are invalid.")
  }
  if (length(delayed_rejection_scale) != 1L ||
      !is.finite(delayed_rejection_scale) || delayed_rejection_scale < 0 ||
      delayed_rejection_scale >= 1) {
    .nm_stop("`delayed_rejection_scale` must lie in [0, 1).")
  }
  if (length(bayes_gibbs_omega) != 1L || is.na(bayes_gibbs_omega)) {
    .nm_stop("`bayes_gibbs_omega` must be TRUE or FALSE.")
  }
  if (delayed_rejection_scale > 0 && !.nm_liber_optimized(context)) {
    .nm_stop("Delayed-rejection BAYES is available only in liber_optimized mode.")
  }
  set.seed(seed)
  stochastic_context <- .nm_stochastic_eta_context(
    context, allow_compatibility = TRUE
  )
  prior <- .nm_prior_evaluator(context$model)
  log_posterior <- function(state) {
    if (!map$in_bounds(state$outer)) {
      return(list(value = -Inf, subject_values = NULL))
    }
    parameters <- state$parameters
    subject_values <- if (!is.null(stochastic_context)) {
      .liberation_stochastic_eta_context_eval(
        stochastic_context, parameters$theta, state$eta,
        parameters$sigma, parameters$omega
      )
    } else tryCatch(
      .nm_saem_conditional_components(context, parameters, state$eta)$subject,
      error = function(e) NULL
    )
    prior_log_density <- prior$log_density(parameters)
    if (is.null(subject_values) || any(!is.finite(subject_values)) ||
        !is.finite(prior_log_density)) {
      return(list(value = -Inf, subject_values = NULL))
    }
    jacobian <- map$log_jacobian(parameters)
    list(
      value = -0.5 * sum(subject_values) + prior_log_density + jacobian,
      subject_values = as.numeric(subject_values)
    )
  }
  state <- .nm_bayes_state(
    map, map$start, matrix(0, context$n_subjects, context$n_eta)
  )
  initial <- log_posterior(state)
  current <- initial$value
  state$subject_values <- initial$subject_values
  mu <- .nm_mu_specialization(context, map, enabled = mu_specialization)
  mu_outer <- if (isTRUE(mu$active) && length(mu$theta)) {
    match(mu$theta, map$theta_free)
  } else integer()
  mu_outer <- mu_outer[!is.na(mu_outer)]
  random_walk_outer <- setdiff(seq_along(map$start), mu_outer)
  resolved_outer_kernel <- if (outer_kernel == "auto") {
    if (.nm_liber_optimized(context) && isTRUE(adapt) &&
        length(random_walk_outer)) {
      "adaptive_metropolis"
    } else "isotropic"
  } else outer_kernel
  resolved_eta_kernel <- if (eta_kernel == "auto") {
    if (.nm_liber_optimized(context) && context$n_eta >= 2L) {
      "laplace"
    } else "random_walk"
  } else eta_kernel
  if (resolved_eta_kernel %in% c("laplace", "student_t") &&
      !.nm_liber_optimized(context)) {
    .nm_stop(
      "The Laplace-independence BAYES ETA kernel is available only in ",
      "liber_optimized mode."
    )
  }
  if (resolved_outer_kernel == "adaptive_metropolis" &&
      !.nm_liber_optimized(context)) {
    .nm_stop(
      "The adaptive BAYES population proposal is available only in ",
      "liber_optimized mode."
    )
  }
  if (resolved_outer_kernel == "adaptive_metropolis" && !isTRUE(adapt)) {
    .nm_stop("Adaptive BAYES requires `adapt = TRUE`.")
  }
  target <- adaptive_target %||% if (length(random_walk_outer) == 1L) {
    0.44
  } else 0.234
  target <- as.numeric(target)
  if (length(target) != 1L || !is.finite(target) ||
      target <= 0 || target >= 1) {
    .nm_stop("`adaptive_target` must lie strictly between zero and one.")
  }
  adaptive_state <- if (
    resolved_outer_kernel == "adaptive_metropolis" &&
      length(random_walk_outer)
  ) {
    .nm_bayes_adaptive_state(
      length(random_walk_outer), step_scale, adaptive_start,
      adaptive_interval, target
    )
  } else NULL
  total_iterations <- n_burn + n_sample * n_thin
  kept <- vector("list", n_sample)
  accepted_outer <- attempted_outer <- accepted_eta <- attempted_eta <- 0L
  accepted_delayed <- attempted_delayed <- 0L
  accepted_mu <- attempted_mu <- 0L
  proposal_cache <- .nm_proposal_root_cache(context)
  eta_current_evaluations <- 0L
  eta_current_cache_hits <- 0L
  eta_candidate_evaluations <- 0L
  eta_proposal <- NULL
  eta_force_refresh <- FALSE
  eta_refreshes <- eta_refresh_failures <- 0L
  eta_rescue_iterations <- eta_acceptance_refreshes <- 0L
  eta_parameter_refreshes <- eta_fallback_iterations <- 0L
  eta_last_error <- NULL
  keep <- 0L
  optimized_native_policy <- .nm_liber_optimized(context)
  compatibility_native_policy <- !optimized_native_policy &&
    identical(resolved_outer_kernel, "isotropic") &&
    identical(resolved_eta_kernel, "random_walk") &&
    identical(delayed_rejection_scale, 0) &&
    !isTRUE(mu$active)
  native_bayes_option <- isTRUE(getOption(
    "LibeRation.bayes_native_coordinator", TRUE
  ))
  native_bayes_eligible <-
    (optimized_native_policy || compatibility_native_policy) &&
    !is.null(stochastic_context) && print_every == 0L &&
    native_bayes_option
  native_bayes_ineligibility <- if (native_bayes_eligible) {
    NULL
  } else if (!native_bayes_option) {
    "native BAYES coordinator disabled by option"
  } else if (print_every > 0L) {
    "iteration printing requires the R BAYES coordinator"
  } else if (is.null(stochastic_context)) {
    "a persistent serial stochastic context is unavailable"
  } else if (!optimized_native_policy && isTRUE(mu$active)) {
    "compatibility MU interweaving remains R-coordinated"
  } else {
    "the selected compatibility controls are not arithmetic-neutral"
  }
  native_gibbs_omega <- optimized_native_policy &&
    isTRUE(bayes_gibbs_omega) &&
    is.null(context$model$RE_CONFIG) &&
    identical(as.integer(context$model$LIK_CONFIG$iov %||% 0L), 0L)
  native_bayes_error <- NULL
  native_bayes_seed <- if (exists(
    ".Random.seed", envir = .GlobalEnv, inherits = FALSE
  )) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  native_bayes <- if (native_bayes_eligible) tryCatch(
    .liberation_stochastic_eta_context_bayes(
      stochastic_context, .nm_bayes_cpp_map_config(context, map, mu),
      n_burn, n_sample, n_thin, step_scale, eta_step, isTRUE(adapt),
      resolved_outer_kernel, adaptive_start, adaptive_interval, target,
      delayed_rejection_scale,
      resolved_eta_kernel, bayes_eta_refresh, bayes_eta_maxit, tolerance,
      bayes_eta_df,
      bayes_eta_rescue_probability, bayes_eta_parameter_refresh,
      bayes_eta_low_acceptance,
      native_gibbs_omega
    ),
    error = function(error) {
      native_bayes_error <<- conditionMessage(error)
      if (!is.null(native_bayes_seed)) {
        assign(".Random.seed", native_bayes_seed, envir = .GlobalEnv)
      }
      NULL
    }
  ) else NULL
  if (is.null(native_bayes)) for (iteration in seq_len(total_iterations)) {
    accepted_outer_iteration <- FALSE
    if (length(random_walk_outer)) {
      proposed_outer <- state$outer
      proposal_root <- if (!is.null(adaptive_state)) adaptive_state$root else
        diag(step_scale, length(random_walk_outer))
      increment <- as.vector(
        proposal_root %*% stats::rnorm(length(random_walk_outer))
      )
      proposed_outer[random_walk_outer] <-
        proposed_outer[random_walk_outer] + increment
      proposal <- .nm_bayes_state(
        map, proposed_outer, state$eta
      )
      proposed_state <- log_posterior(proposal)
      proposed <- proposed_state$value
      attempted_outer <- attempted_outer + 1L
      first_log_alpha <- min(0, proposed - current)
      if (log(stats::runif(1)) < first_log_alpha) {
        proposal$subject_values <- proposed_state$subject_values
        state <- proposal
        current <- proposed
        accepted_outer <- accepted_outer + 1L
        accepted_outer_iteration <- TRUE
      } else if (delayed_rejection_scale > 0) {
        attempted_delayed <- attempted_delayed + 1L
        second_outer <- state$outer
        second_outer[random_walk_outer] <-
          second_outer[random_walk_outer] + as.vector(
            delayed_rejection_scale * proposal_root %*%
              stats::rnorm(length(random_walk_outer))
          )
        second <- .nm_bayes_state(map, second_outer, state$eta)
        second_state <- log_posterior(second)
        second_logp <- second_state$value
        attempted_outer <- attempted_outer + 1L
        if (is.finite(second_logp)) {
          gaussian_quad <- function(value) {
            standardized <- forwardsolve(proposal_root, value)
            drop(crossprod(standardized))
          }
          log_one_minus <- function(log_alpha) {
            if (log_alpha >= 0) return(-Inf)
            if (!is.finite(log_alpha)) return(0)
            if (log_alpha < -log(2)) {
              log1p(-exp(log_alpha))
            } else log(-expm1(log_alpha))
          }
          first_from_current <- proposed_outer[random_walk_outer] -
            state$outer[random_walk_outer]
          first_from_second <- proposed_outer[random_walk_outer] -
            second_outer[random_walk_outer]
          reverse_first <- min(0, proposed - second_logp)
          correction <- second_logp - current -
            0.5 * gaussian_quad(first_from_second) +
            0.5 * gaussian_quad(first_from_current) +
            log_one_minus(reverse_first) - log_one_minus(first_log_alpha)
          if (log(stats::runif(1)) < min(0, correction)) {
            second$subject_values <- second_state$subject_values
            state <- second
            current <- second_logp
            accepted_outer <- accepted_outer + 1L
            accepted_delayed <- accepted_delayed + 1L
            accepted_outer_iteration <- TRUE
          }
        }
      }
    }
    if (isTRUE(mu$active) && length(mu$theta)) {
      proposed_mu <- .nm_mu_bayes_proposal(
        mu, context, state, map, log_posterior, current
      )
      if (isTRUE(proposed_mu$attempted)) {
        attempted_mu <- attempted_mu + 1L
      }
      if (isTRUE(proposed_mu$accepted)) {
        accepted_mu <- accepted_mu + 1L
        state <- proposed_mu$state
        current <- proposed_mu$log_posterior
      }
    }
    if (context$n_eta) {
      independence <- NULL
      parameter_refresh <- FALSE
      if (resolved_eta_kernel %in% c("laplace", "student_t") &&
          !is.null(eta_proposal)) {
        anchor <- c(
          state$parameters$theta, state$parameters$sigma,
          state$parameters$omega
        )
        parameter_refresh <- max(
          abs(anchor - eta_proposal$anchor) / (1 + abs(eta_proposal$anchor))
        ) > bayes_eta_parameter_refresh
      }
      if (resolved_eta_kernel %in% c("laplace", "student_t") &&
          (is.null(eta_proposal) || eta_force_refresh || parameter_refresh ||
           (iteration - 1L) %% bayes_eta_refresh == 0L)) {
        if (parameter_refresh) {
          eta_parameter_refreshes <- eta_parameter_refreshes + 1L
        }
        refreshed <- tryCatch(
          .nm_fsaem_proposal(
            context, state$parameters,
            if (is.null(eta_proposal)) state$eta else eta_proposal$modes,
            bayes_eta_maxit, tolerance,
            persistent = if (isTRUE(context$model$USE_ODE)) {
              NULL
            } else stochastic_context
          ),
          error = identity
        )
        if (inherits(refreshed, "error")) {
          eta_refresh_failures <- eta_refresh_failures + 1L
          eta_last_error <- conditionMessage(refreshed)
        } else {
          eta_proposal <- refreshed
          eta_refreshes <- eta_refreshes + 1L
          eta_force_refresh <- FALSE
        }
      }
      if (resolved_eta_kernel %in% c("laplace", "student_t")) {
        independence <- eta_proposal
        if (!is.null(independence)) {
          independence$df <- if (resolved_eta_kernel == "student_t") {
            bayes_eta_df
          } else Inf
        }
        if (is.null(independence)) {
          eta_fallback_iterations <- eta_fallback_iterations + 1L
        } else if (bayes_eta_rescue_probability > 0 &&
                   stats::runif(1) < bayes_eta_rescue_probability) {
          independence <- NULL
          eta_rescue_iterations <- eta_rescue_iterations + 1L
        }
      }
      # Conditional independence makes all subject ETA updates separable at a
      # fixed population parameter point.  The batched C++ kernel evaluates
      # each subject tape exactly once per proposal instead of re-evaluating a
      # full population tape N times during every sweep.
      previous_subject_total <- sum(state$subject_values)
      sampled <- .nm_saem_metropolis(
        context, state$parameters, state$eta, mcmc_steps = 1L,
        step_scale = eta_step, proposal_cache = proposal_cache,
        current_values = state$subject_values,
        persistent = if (isTRUE(context$model$USE_ODE)) NULL else
          stochastic_context,
        independence = independence
      )
      state$eta <- sampled$eta
      state$subject_values <- sampled$value
      accepted_eta <- accepted_eta + sampled$accepted
      attempted_eta <- attempted_eta + sampled$attempted
      eta_current_evaluations <- eta_current_evaluations +
        as.integer(sampled$current_evaluations %||% 0L)
      eta_current_cache_hits <- eta_current_cache_hits +
        as.integer(sampled$current_cache_hits %||% 0L)
      eta_candidate_evaluations <- eta_candidate_evaluations +
        as.integer(sampled$candidate_evaluations %||% 0L)
      eta_acceptance <- sampled$accepted / max(sampled$attempted, 1L)
      if (resolved_eta_kernel %in% c("laplace", "student_t") &&
          !is.null(independence) &&
          eta_acceptance < bayes_eta_low_acceptance) {
        eta_force_refresh <- TRUE
        eta_acceptance_refreshes <- eta_acceptance_refreshes + 1L
      }
      current <- current - 0.5 * (
        sum(sampled$value) - previous_subject_total
      )
    }
    if (!is.null(adaptive_state) && length(random_walk_outer)) {
      adaptive_state$update(
        state$outer[random_walk_outer], accepted_outer_iteration,
        adapting = iteration <= n_burn
      )
    } else if (isTRUE(adapt) && length(random_walk_outer) &&
        iteration <= n_burn && iteration %% 50L == 0L) {
      rate <- accepted_outer / max(attempted_outer, 1L)
      step_scale <- step_scale * exp(if (rate > 0.3) 0.1 else -0.1)
    }
    if (iteration > n_burn && (iteration - n_burn) %% n_thin == 0L) {
      keep <- keep + 1L
      kept[[keep]] <- c(
        state$parameters$theta, state$parameters$sigma, state$parameters$omega,
        as.vector(t(state$eta)), LOG_POSTERIOR = current
      )
    }
    if (print_every > 0L && iteration %% print_every == 0L) {
      population <- tryCatch(
        .nm_conditional_native_gradient(
          context, state$parameters, state$eta, interaction = TRUE
        ),
        error = function(error) rep(
          NA_real_, length(state$parameters$theta) +
            length(state$parameters$sigma) + length(state$parameters$omega)
        )
      )
      names(population) <- .nm_parameter_names(
        state$parameters$theta, state$parameters$sigma, state$parameters$omega
      )
      cat(sprintf(
        "[LibeRation] MCMC ITERATION %d -2LOGPOST %.10g GRADIENT %s\n",
        iteration, -2 * current,
        paste(sprintf("%s=%.6g", names(population), population), collapse = " ")
      ))
      try(flush(stdout()), silent = TRUE)
    }
  }
  chain <- if (!is.null(native_bayes)) {
    accepted_outer <- as.integer(native_bayes$accepted_outer)
    attempted_outer <- as.integer(native_bayes$attempted_outer)
    accepted_eta <- as.integer(native_bayes$accepted_eta)
    attempted_eta <- as.integer(native_bayes$attempted_eta)
    accepted_mu <- as.integer(native_bayes$accepted_mu %||% 0L)
    attempted_mu <- as.integer(native_bayes$attempted_mu %||% 0L)
    accepted_delayed <- as.integer(native_bayes$accepted_delayed %||% 0L)
    attempted_delayed <- as.integer(native_bayes$attempted_delayed %||% 0L)
    eta_current_evaluations <- 0L
    eta_current_cache_hits <- attempted_eta
    eta_candidate_evaluations <- attempted_eta
    as.matrix(native_bayes$chain)
  } else do.call(rbind, kept)
  n_theta <- nrow(context$model$THETAS)
  n_sigma <- nrow(context$model$SIGMAS)
  n_omega <- nrow(context$model$OMEGAS)
  colnames(chain) <- c(
    .nm_numbered_names("THETA", n_theta), .nm_numbered_names("SIGMA", n_sigma),
    .nm_numbered_names("OMEGA", n_omega),
    if (context$n_eta) unlist(lapply(seq_len(context$n_subjects), function(subject) {
      paste0("ETA", subject, "_", seq_len(context$n_eta))
    })) else character(), "LOG_POSTERIOR"
  )
  parameters <- list(
    theta = colMeans(chain[, seq_len(n_theta), drop = FALSE]),
    sigma = colMeans(chain[, n_theta + seq_len(n_sigma), drop = FALSE]),
    omega = colMeans(chain[, n_theta + n_sigma + seq_len(n_omega), drop = FALSE])
  )
  eta_start <- n_theta + n_sigma + n_omega
  eta <- if (context$n_eta) matrix(
    colMeans(chain[, eta_start + seq_len(context$n_subjects * context$n_eta), drop = FALSE]),
    context$n_subjects, context$n_eta, byrow = TRUE
  ) else matrix(numeric(), context$n_subjects, 0L)
  final_state <- .nm_bayes_state(map, map$encode(parameters), eta)
  final_objective <- -2 * log_posterior(final_state)$value
  modes <- lapply(seq_len(context$n_subjects), function(subject) {
    list(par = eta[subject, ], convergence = 0L, jitter = 0)
  })
  optimizer <- list(
    convergence = 0L, message = "Bayesian sampling completed",
    counts = c(`function` = total_iterations, gradient = NA_integer_),
    iterations = total_iterations, objective_evaluations = total_iterations,
    backend = native_bayes$backend %||% "r-coordinated-bayes"
  )
  fit <- .nm_fit_result(
    context, "BAYES", parameters, final_objective, modes, optimizer,
    diagnostics = list(
      outer_acceptance = native_bayes$outer_acceptance %||%
        (accepted_outer / max(attempted_outer, 1L)),
      eta_acceptance = native_bayes$eta_acceptance %||%
        (accepted_eta / max(attempted_eta, 1L)),
      mu_acceptance = accepted_mu / max(attempted_mu, 1L),
      mu_specialization = c(
        .nm_mu_diagnostic(mu),
        list(
          attempted_blocks = attempted_mu,
          accepted_blocks = accepted_mu,
          random_walk_outer_parameters = length(random_walk_outer)
        )
      ),
      n_burn = n_burn, n_sample = n_sample, n_thin = n_thin,
      seed = seed,
      final_step_scale = native_bayes$final_step_scale %||% if (!is.null(adaptive_state)) {
        exp(adaptive_state$log_multiplier)
      } else step_scale,
      outer_sampler = list(
        requested_kernel = outer_kernel,
        resolved_kernel = resolved_outer_kernel,
        target_acceptance = target,
        adaptive_start = adaptive_start,
        adaptive_interval = adaptive_interval,
        covariance_updates = native_bayes$covariance_updates %||%
          (adaptive_state$root_updates %||% 0L),
        covariance_regularizations = native_bayes$covariance_regularizations %||%
          (adaptive_state$regularizations %||% 0L),
        covariance = native_bayes$covariance %||%
          (adaptive_state$covariance %||% NULL),
        multiplier = native_bayes$multiplier %||% if (!is.null(adaptive_state)) {
          exp(adaptive_state$log_multiplier)
        } else 1,
        delayed_rejection = list(
          scale = delayed_rejection_scale,
          attempted = attempted_delayed,
          accepted = accepted_delayed,
          acceptance = accepted_delayed / max(attempted_delayed, 1L)
        )
      ),
      # Retain the established top-level diagnostic label for downstream
      # consumers; the selected algorithm is exposed in eta_sampler below.
      eta_kernel = "batched conditional-subject C++ Metropolis",
      eta_sampler = list(
        backend = if (!is.null(native_bayes)) {
          native_bayes$backend
        } else if (!is.null(stochastic_context)) {
          "persistent cached conditional-subject C++ Metropolis"
        } else "cached batched conditional-subject C++ Metropolis",
        current_evaluations = eta_current_evaluations,
        current_cache_hits = eta_current_cache_hits,
        candidate_evaluations = eta_candidate_evaluations,
        proposal_root_cache_hits = if (!is.null(native_bayes)) {
          as.integer(native_bayes$omega_cache_hits %||% 0L)
        } else as.integer(proposal_cache$hits),
        proposal_root_cache_misses = as.integer(proposal_cache$misses),
        proposal_root_factorizations = if (!is.null(native_bayes)) {
          as.integer(native_bayes$omega_factorizations %||% 0L)
        } else as.integer(proposal_cache$factorizations),
        proposal_root_groups = length(unique(proposal_cache$groups)),
        requested_kernel = eta_kernel,
        resolved_kernel = resolved_eta_kernel,
        laplace = list(
          refresh_every = bayes_eta_refresh,
          refreshes = native_bayes$eta_refreshes %||% eta_refreshes,
          refresh_failures = native_bayes$eta_refresh_failures %||%
            eta_refresh_failures,
          fallback_iterations = native_bayes$eta_fallback_iterations %||%
            eta_fallback_iterations,
          rescue_probability = bayes_eta_rescue_probability,
          rescue_iterations = native_bayes$eta_rescue_iterations %||%
            eta_rescue_iterations,
          parameter_refresh_threshold = bayes_eta_parameter_refresh,
          parameter_triggered_refreshes =
            native_bayes$eta_parameter_refreshes %||%
              eta_parameter_refreshes,
          low_acceptance_threshold = bayes_eta_low_acceptance,
          acceptance_triggered_refreshes =
            native_bayes$eta_acceptance_refreshes %||%
              eta_acceptance_refreshes,
          last_error = native_bayes$eta_last_error %||% eta_last_error
        ),
        native_coordinator = list(
          eligible = native_bayes_eligible,
          used = !is.null(native_bayes),
          policy = if (!is.null(native_bayes)) {
            if (compatibility_native_policy) {
              "compatibility-preserving-native"
            } else "optimized-native"
          } else "R-coordinated",
          compatibility_preserving = compatibility_native_policy,
          fallback_reason = native_bayes_error %||%
            native_bayes_ineligibility,
          omega_factorizations = native_bayes$omega_factorizations %||% 0L,
          omega_cache_hits = native_bayes$omega_cache_hits %||% 0L,
          conjugate_omega = list(
            requested = isTRUE(bayes_gibbs_omega),
            enabled = native_gibbs_omega,
            parameters = native_bayes$gibbs_omega_parameters %||% 0L,
            updates = native_bayes$gibbs_omega_updates %||% 0L,
            draws = native_bayes$gibbs_omega_draws %||% 0L
          )
        ),
        persistent_context = if (!is.null(stochastic_context)) {
          .liberation_stochastic_eta_context_telemetry(stochastic_context)
        } else NULL
      ),
      objective_semantics = list(
        type = "negative_twice_log_posterior_at_posterior_mean_sampling_coordinates",
        likelihood_comparable = FALSE,
        reported_point = "posterior mean population parameters and ETAs",
        recommendation = "Use WAIC or PSIS-LOO for predictive model comparison."
      )
    )
  )
  fit$chain <- chain
  population_names <- .nm_parameter_names(
    parameters$theta, parameters$sigma, parameters$omega
  )
  population_chain <- chain[, population_names, drop = FALSE]
  population_covariance <- if (nrow(population_chain) > 1L) {
    stats::cov(population_chain)
  } else {
    matrix(NA_real_, ncol(population_chain), ncol(population_chain),
           dimnames = list(population_names, population_names))
  }
  population_sd <- apply(population_chain, 2, stats::sd)
  population_correlation <- population_covariance / outer(population_sd, population_sd)
  diag(population_correlation) <- 1
  sampling_diagnostics <- .nm_mcmc_diagnostics(list(population_chain))
  sampling_diagnostics <- lapply(sampling_diagnostics, function(value) {
    stats::setNames(value, population_names)
  })
  fit$posterior <- list(
    mean = colMeans(chain),
    sd = apply(chain, 2, stats::sd),
    quantile = apply(chain, 2, stats::quantile, probs = c(0.025, 0.5, 0.975)),
    population = list(
      mean = colMeans(population_chain), sd = population_sd,
      quantile = apply(
        population_chain, 2, stats::quantile, probs = c(0.025, 0.5, 0.975)
      ),
      covariance = population_covariance,
      correlation = population_correlation,
      rhat = sampling_diagnostics$rhat,
      ess = sampling_diagnostics$bulk_ess,
      bulk_ess = sampling_diagnostics$bulk_ess,
      tail_ess = sampling_diagnostics$tail_ess,
      mcse_mean = sampling_diagnostics$mcse_mean,
      diagnostics_method = "rank-normalized split R-hat and bulk/tail ESS"
    )
  )
  fit
}

.nm_est_bayes <- function(context, map, tolerance, n_chains = 1L,
                          parallel_chains = FALSE,
                          seed = 20260713L, ...) {
  n_chains <- as.integer(n_chains)
  seed <- as.integer(seed)
  if (length(n_chains) != 1L || is.na(n_chains) || n_chains < 1L) {
    .nm_stop("`n_chains` must be a positive integer for BAYES.")
  }
  if (length(seed) != 1L || is.na(seed) ||
      seed > .Machine$integer.max - n_chains + 1L) {
    .nm_stop("`seed` must leave room for one independent seed per BAYES chain.")
  }
  if (length(parallel_chains) != 1L || is.na(parallel_chains)) {
    .nm_stop("`parallel_chains` must be TRUE or FALSE.")
  }
  if (isTRUE(parallel_chains) && n_chains > 1L) {
    # Subject-parallel contexts own live worker connections and compiled tape
    # pointers, so they cannot safely be serialized into another process.
    # Independent chains are deliberately sequenced here; an outer job queue
    # can still run them concurrently without sharing mutable CppAD state.
    warning(
      "BAYES chains are being run sequentially in this R session. For true ",
      "chain-level parallelism, submit independent seeded jobs through the ",
      "local or remote queue.", call. = FALSE
    )
  }
  chain_fits <- lapply(seq_len(n_chains), function(chain_id) {
    .nm_est_bayes_single(
      context = context, map = map, tolerance = tolerance,
      seed = seed + chain_id - 1L, ...
    )
  })
  if (n_chains == 1L) {
    chain_fits[[1L]]$chains <- list(chain_fits[[1L]]$chain)
    chain_fits[[1L]]$diagnostics$n_chains <- 1L
    return(chain_fits[[1L]])
  }
  chains <- lapply(chain_fits, `[[`, "chain")
  chain <- do.call(rbind, chains)
  n_theta <- nrow(context$model$THETAS)
  n_sigma <- nrow(context$model$SIGMAS)
  n_omega <- nrow(context$model$OMEGAS)
  theta_names <- .nm_numbered_names("THETA", n_theta)
  sigma_names <- .nm_numbered_names("SIGMA", n_sigma)
  omega_names <- .nm_numbered_names("OMEGA", n_omega)
  parameters <- list(
    theta = colMeans(chain[, theta_names, drop = FALSE]),
    sigma = colMeans(chain[, sigma_names, drop = FALSE]),
    omega = colMeans(chain[, omega_names, drop = FALSE])
  )
  eta_names <- if (context$n_eta) unlist(lapply(
    seq_len(context$n_subjects), function(subject) {
      paste0("ETA", subject, "_", seq_len(context$n_eta))
    }
  )) else character()
  eta <- if (length(eta_names)) matrix(
    colMeans(chain[, eta_names, drop = FALSE]),
    context$n_subjects, context$n_eta, byrow = TRUE
  ) else matrix(numeric(), context$n_subjects, 0L)
  result <- chain_fits[[1L]]
  result$theta <- parameters$theta
  result$sigma <- parameters$sigma
  result$omega <- parameters$omega
  result$eta <- eta
  mean_subject <- tryCatch(
    .nm_saem_conditional_components(context, parameters, eta)$subject,
    error = function(error) NULL
  )
  mean_prior <- .nm_prior_evaluator(context$model)$log_density(parameters)
  if (!is.null(mean_subject) && all(is.finite(mean_subject)) &&
      is.finite(mean_prior)) {
    result$objective <- -2 * (
      -0.5 * sum(mean_subject) + mean_prior + map$log_jacobian(parameters)
    )
  }
  result$chain <- chain
  result$chains <- chains
  population_names <- c(theta_names, sigma_names, omega_names)
  population_chains <- lapply(chains, function(value) {
    value[, population_names, drop = FALSE]
  })
  population_chain <- chain[, population_names, drop = FALSE]
  population_covariance <- stats::cov(population_chain)
  population_sd <- apply(population_chain, 2, stats::sd)
  population_correlation <- population_covariance /
    outer(population_sd, population_sd)
  diag(population_correlation) <- 1
  sampling_diagnostics <- .nm_mcmc_diagnostics(population_chains)
  sampling_diagnostics <- lapply(sampling_diagnostics, function(value) {
    stats::setNames(value, population_names)
  })
  result$posterior <- list(
    mean = colMeans(chain), sd = apply(chain, 2, stats::sd),
    quantile = apply(
      chain, 2, stats::quantile, probs = c(0.025, 0.5, 0.975)
    ),
    population = list(
      mean = colMeans(population_chain), sd = population_sd,
      quantile = apply(
        population_chain, 2, stats::quantile,
        probs = c(0.025, 0.5, 0.975)
      ),
      covariance = population_covariance,
      correlation = population_correlation,
      rhat = sampling_diagnostics$rhat,
      ess = sampling_diagnostics$bulk_ess,
      bulk_ess = sampling_diagnostics$bulk_ess,
      tail_ess = sampling_diagnostics$tail_ess,
      mcse_mean = sampling_diagnostics$mcse_mean,
      diagnostics_method = "rank-normalized split R-hat and bulk/tail ESS"
    )
  )
  result$diagnostics$n_chains <- n_chains
  result$diagnostics$chain_seeds <- seed + seq_len(n_chains) - 1L
  result$diagnostics$chains <- lapply(chain_fits, function(value) list(
    outer_acceptance = value$diagnostics$outer_acceptance,
    eta_acceptance = value$diagnostics$eta_acceptance,
    mu_acceptance = value$diagnostics$mu_acceptance,
    backend = value$diagnostics$optimizer$backend
  ))
  result$diagnostics$outer_acceptance <- mean(vapply(
    chain_fits, function(value) value$diagnostics$outer_acceptance, numeric(1)
  ))
  result$diagnostics$eta_acceptance <- mean(vapply(
    chain_fits, function(value) value$diagnostics$eta_acceptance, numeric(1)
  ))
  result$diagnostics$mu_acceptance <- mean(vapply(
    chain_fits, function(value) value$diagnostics$mu_acceptance, numeric(1)
  ))
  result
}

.nm_est_stochastic <- function(context, map, method, maxit, eta_maxit,
                               tolerance, trace, print_every = 0L,
                               optimizer_backend = "auto", initial_eta = NULL, ...) {
  controls <- list(...)
  if (method == "ITS") {
    return(do.call(.nm_est_its, c(list(
      context = context, map = map, maxit = maxit, eta_maxit = eta_maxit,
      tolerance = tolerance, trace = trace, print_every = print_every,
      optimizer_backend = optimizer_backend
    ), controls)))
  }
  if (method == "IMP") {
    return(do.call(.nm_est_imp, c(list(
      context = context, map = map, maxit = maxit, eta_maxit = eta_maxit,
      tolerance = tolerance, trace = trace, print_every = print_every,
      optimizer_backend = optimizer_backend
    ), controls)))
  }
  if (method == "GQ") {
    return(do.call(.nm_est_gq, c(list(
      context = context, map = map, maxit = maxit, eta_maxit = eta_maxit,
      tolerance = tolerance, trace = trace, print_every = print_every,
      optimizer_backend = optimizer_backend
    ), controls)))
  }
  if (method == "SAEM") {
    return(do.call(.nm_est_saem, c(list(
      context = context, map = map, maxit = maxit,
      tolerance = tolerance, trace = trace, print_every = print_every,
      optimizer_backend = optimizer_backend, initial_eta = initial_eta
    ), controls)))
  }
  if (method == "BAYES") {
    return(do.call(.nm_est_bayes, c(list(
      context = context, map = map, tolerance = tolerance,
      print_every = print_every
    ), controls)))
  }
  if (method %in% c("HMC", "NUTS")) {
    return(do.call(.nm_est_hmc, c(list(
      context = context, map = map, method = method,
      print_every = print_every
    ), controls)))
  }
  if (method %in% c("NPML", "NPAG")) {
    return(do.call(.nm_est_nonparametric, c(list(
      context = context, method = method, maxit = maxit,
      tolerance = tolerance, trace = trace, print_every = print_every,
      optimizer_backend = optimizer_backend, eta_maxit = eta_maxit
    ), controls)))
  }
  .nm_stop("Unknown stochastic estimation method: ", method)
}
