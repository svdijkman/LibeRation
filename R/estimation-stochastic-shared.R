# Shared stochastic-estimation helpers and method dispatch.
# Split from estimation-stochastic.R as a behaviour-preserving source move.

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
