.nm_log_mean_exp <- function(values) {
  maximum <- max(values)
  maximum + log(mean(exp(values - maximum)))
}

.nm_log_sum_exp <- function(values) {
  maximum <- max(values)
  maximum + log(sum(exp(values - maximum)))
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

.nm_imp_normals <- function(context, n_imp, seed) {
  n_imp <- as.integer(n_imp)
  seed <- as.integer(seed)
  if (length(n_imp) != 1L || is.na(n_imp) || n_imp < 5L) {
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
    matrix(stats::rnorm(n_imp * context$n_eta), n_imp, context$n_eta)
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

.nm_est_its <- function(context, map, maxit, eta_maxit, tolerance, trace,
                        print_every = 0L, optimizer_backend = "auto") {
  objective <- .nm_nested_objective(context, "its", eta_maxit, tolerance)
  gradient <- function(parameters) .nm_nested_outer_gradient(
    context, map, objective, parameters, "its"
  )
  compiled <- .nm_cpp_population_objective(
    context, map, "its", eta_maxit, tolerance
  )
  optimizer <- .nm_outer_optim(
    map, objective, maxit, tolerance, trace, print_every,
    gradient = gradient, optimizer_backend = optimizer_backend,
    compiled_objective = compiled
  )
  parameters <- map$decode(optimizer$par)
  cached <- if (!is.null(compiled$pointer)) tryCatch(
    .liberation_population_objective_state(compiled$pointer, optimizer$par),
    error = function(error) NULL
  ) else NULL
  modes <- cached$modes %||% .nm_subject_modes(
    context, parameters, maxit = eta_maxit, tolerance = tolerance,
    exact_hessian = FALSE
  )
  work <- optimizer$population_objective
  if (is.null(work)) {
    work <- list(
      value_requests = attr(objective, "state")$objective_calls,
      value_cache_hits = attr(objective, "state")$cache_hits,
      mode_iterations = attr(objective, "state")$mode_iterations,
      mode_evaluations = attr(objective, "state")$mode_evaluations
    )
  }
  .nm_fit_result(
    context, "ITS", parameters, optimizer$value, modes, optimizer,
    diagnostics = list(
      eta_convergence = vapply(modes, `[[`, integer(1), "convergence"),
      description = "iterative conditional modes without a Laplace determinant",
      conditional_mode_work = list(
        objective_calls = work$value_requests %||% work$parameter_evaluations,
        cache_hits = sum(
          work$value_cache_hits %||% 0L,
          work$gradient_cache_hits %||% 0L,
          work$shared_state_hits %||% 0L
        ),
        iterations = work$mode_iterations,
        evaluations = work$mode_evaluations
      ),
      population_gradient = "exact envelope gradient"
    )
  )
}

.nm_imp_subject_objective <- function(evaluator, parameters, normals,
                                      eta_maxit, tolerance) {
  .nm_imp_subject_state(
    evaluator, parameters, normals, eta_maxit, tolerance, gradient = FALSE
  )$value
}

.nm_imp_subject_proposal <- function(evaluator, parameters, normals,
                                     eta_maxit, tolerance) {
  mode <- evaluator$eta_mode(
    parameters$theta, parameters$sigma, parameters$omega,
    maxit = eta_maxit, tolerance = tolerance
  )
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
  log_measure <- attr(normals, "log_measure", exact = TRUE)
  measure_sign <- attr(normals, "measure_sign", exact = TRUE)
  sampling <- attr(normals, "quadrature_method", exact = TRUE)
  if (is.null(log_measure)) sampling <- "random-normal"
  if (is.null(sampling)) sampling <- "tensor-gauss-hermite"
  if (is.null(log_measure)) log_measure <- rep(-log(nrow(z)), nrow(z))
  if (is.null(measure_sign)) measure_sign <- rep(1, nrow(z))
  list(
    valid = TRUE, mode = mode,
    eta = sweep(z %*% t(root), 2L, mode$par, `+`),
    log_proposal = -0.5 * (
      dimension * log(2 * pi) + logdet + rowSums(z^2)
    ),
    log_measure = as.numeric(log_measure),
    measure_sign = as.numeric(measure_sign), sampling = sampling
  )
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
                                      eta_maxit, tolerance, adaptive = TRUE) {
  prepare_chunk <- function(evaluators, chunk_normals) {
    lapply(seq_along(evaluators), function(subject) {
      if (isTRUE(adaptive)) {
        .nm_imp_subject_proposal(
          evaluators[[subject]], parameters, chunk_normals[[subject]],
          eta_maxit, tolerance
        )
      } else {
        .nm_gq_fixed_subject_proposal(
          evaluators[[subject]], parameters, chunk_normals[[subject]]
        )
      }
    })
  }
  if (is.null(context$parallel)) {
    return(prepare_chunk(context$subjects, normals))
  }
  chunks <- context$parallel$chunks
  normal_chunks <- lapply(chunks, function(rows) normals[rows])
  pieces <- parallel::clusterApply(
    context$parallel$cluster, seq_along(chunks),
      function(index, chunks, parameters, eta_maxit, tolerance, adaptive) {
        evaluators <- get(".liber_parallel_subjects", envir = .GlobalEnv)
        prepare <- get(
          if (isTRUE(adaptive)) ".nm_imp_subject_proposal" else
            ".nm_gq_fixed_subject_proposal",
          envir = asNamespace("LibeRation")
        )
        lapply(seq_along(evaluators), function(subject) {
          if (isTRUE(adaptive)) {
            prepare(
              evaluators[[subject]], parameters, chunks[[index]][[subject]],
              eta_maxit, tolerance
            )
          } else {
            prepare(evaluators[[subject]], parameters, chunks[[index]][[subject]])
          }
        })
      }, chunks = normal_chunks, parameters = parameters,
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
  if (is.null(context$parallel)) {
    states <- evaluate_chunk(context$subjects, proposals)
  } else {
    chunks <- context$parallel$chunks
    proposal_chunks <- lapply(chunks, function(rows) proposals[rows])
    pieces <- parallel::clusterApply(
      context$parallel$cluster, seq_along(chunks),
      function(index, chunks, parameters, gradient) {
        evaluators <- get(".liber_parallel_subjects", envir = .GlobalEnv)
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
        evaluators <- get(".liber_parallel_subjects", envir = .GlobalEnv)
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
                             gradient = TRUE) {
  evaluate_chunk <- function(evaluators, chunk_normals) {
    lapply(seq_along(evaluators), function(subject) {
      .nm_imp_subject_state(
        evaluators[[subject]], parameters, chunk_normals[[subject]],
        eta_maxit, tolerance, gradient = gradient
      )
    })
  }
  if (is.null(context$parallel)) {
    states <- evaluate_chunk(context$subjects, normals)
  } else {
    normal_chunks <- lapply(context$parallel$chunks, function(rows) normals[rows])
    pieces <- parallel::clusterApply(
      context$parallel$cluster, seq_along(context$parallel$chunks),
      function(index, chunks, parameters, eta_maxit, tolerance, gradient) {
        evaluators <- get(".liber_parallel_subjects", envir = .GlobalEnv)
        state <- get(".nm_imp_subject_state", envir = asNamespace("LibeRation"))
        lapply(seq_along(evaluators), function(subject) {
          state(
            evaluators[[subject]], parameters, chunks[[index]][[subject]],
            eta_maxit, tolerance, gradient = gradient
          )
        })
      }, chunks = normal_chunks, parameters = parameters,
      eta_maxit = eta_maxit, tolerance = tolerance, gradient = gradient
    )
    states <- unlist(pieces, recursive = FALSE)
  }
  value <- sum(vapply(states, `[[`, numeric(1), "value")) +
    .nm_prior_nll(context$model, parameters)
  if (!isTRUE(gradient)) return(list(value = value, states = states))
  gradients <- lapply(states, `[[`, "native_gradient")
  if (any(vapply(gradients, is.null, logical(1)))) {
    return(list(value = value, native_gradient = NULL, states = states))
  }
  full <- Reduce(`+`, gradients)
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

.nm_est_imp <- function(context, map, maxit, eta_maxit, tolerance, trace,
                        n_imp = 200L, seed = 20260713L, print_every = 0L,
                        imp_gradient = c("score", "finite_crn"),
                        optimizer_backend = "auto") {
  n_imp <- as.integer(n_imp)
  if (n_imp < 5L) .nm_stop("IMP requires `n_imp >= 5`.")
  imp_gradient <- match.arg(imp_gradient)
  normals <- .nm_imp_normals(context, n_imp, seed)
  cache <- new.env(parent = emptyenv())
  cache$key <- NULL
  evaluate <- function(parameters) {
    key <- c(parameters$theta, parameters$sigma, parameters$omega)
    if (is.null(cache$key) || !identical(cache$key, key)) {
      cache$result <- .nm_imp_evaluate(
        context, parameters, normals, eta_maxit, tolerance,
        gradient = imp_gradient == "score"
      )
      cache$key <- key
    }
    cache$result
  }
  objective <- function(parameters) evaluate(parameters)$value
  gradient <- if (imp_gradient == "score") function(parameters) {
    result <- evaluate(parameters)
    as.vector(result$native_gradient %*% map$jacobian(parameters))
  } else NULL
  optimizer <- .nm_outer_optim(
    map, objective, maxit, tolerance, trace, print_every,
    gradient = gradient, optimizer_backend = optimizer_backend
  )
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
      population_gradient = if (imp_gradient == "score") {
        "normalized importance-score CppAD gradient (proposal derivative omitted)"
      } else "finite common-random-number objective"
    )
  )
}

.nm_gq_evaluate <- function(context, parameters, normals, eta_maxit, tolerance,
                            adaptive = TRUE, gradient = TRUE) {
  proposals <- .nm_imp_prepare_proposals(
    context, parameters, normals, eta_maxit, tolerance, adaptive = adaptive
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
                       gq_gradient = c("score", "finite_grid"),
                       print_every = 0L, optimizer_backend = "auto") {
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
  design <- .nm_gq_design(
    context, order = gq_order, max_points = gq_max_points,
    grid = gq_grid, level = gq_level
  )
  cache <- new.env(parent = emptyenv())
  cache$key <- NULL
  evaluate <- function(parameters) {
    key <- c(parameters$theta, parameters$sigma, parameters$omega)
    if (is.null(cache$key) || !identical(cache$key, key)) {
      cache$result <- .nm_gq_evaluate(
        context, parameters, design$normals, eta_maxit, tolerance,
        adaptive = gq_adaptive, gradient = gq_gradient == "score"
      )
      cache$key <- key
    }
    cache$result
  }
  objective <- function(parameters) evaluate(parameters)$value
  gradient <- if (gq_gradient == "score") function(parameters) {
    result <- evaluate(parameters)
    if (is.null(result$native_gradient)) return(NULL)
    as.vector(result$native_gradient %*% map$jacobian(parameters))
  } else NULL
  optimizer <- .nm_outer_optim(
    map, objective, maxit, tolerance, trace, print_every,
    gradient = gradient, optimizer_backend = optimizer_backend
  )
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
      effective_quadrature_points = effective,
      quadrature_cancellation_ratio = cancellation,
      population_gradient = if (gq_gradient == "score") {
        "normalized quadrature-score CppAD gradient (node derivative omitted)"
      } else {
        "finite adaptive-grid objective"
      }
    )
  )
}

.nm_saem_conditional <- function(context, parameters, eta) {
  if (is.null(context$parallel)) {
    value <- sum(.nm_objective_collection(context$subjects, parameters, eta))
  } else {
    eta_chunks <- lapply(
      context$parallel$chunks, function(rows) eta[rows, , drop = FALSE]
    )
    pieces <- parallel::clusterApply(
      context$parallel$cluster, seq_along(context$parallel$chunks),
      function(index, eta_chunks, parameters) {
        evaluators <- get(".liber_parallel_subjects", envir = .GlobalEnv)
        worker_eta <- eta_chunks[[index]]
        collection <- get(".nm_objective_collection", envir = asNamespace("LibeRation"))
        sum(collection(evaluators, parameters, worker_eta))
      }, eta_chunks = eta_chunks, parameters = parameters
    )
    value <- sum(unlist(pieces, use.names = FALSE))
  }
  value + .nm_prior_nll(context$model, parameters)
}

.nm_saem_conditional_gradient <- function(context, map, parameters, eta) {
  native <- .nm_conditional_native_gradient(
    context, parameters, eta, interaction = TRUE
  )
  as.vector(native %*% map$jacobian(parameters))
}

.nm_saem_omega_sufficient <- function(context, eta) {
  iov <- context$model$LIK_CONFIG$iov
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
                                      mcmc_steps, step_scale) {
  if (!length(evaluators) || !ncol(eta)) {
    return(list(eta = eta, accepted = 0L, attempted = 0L))
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
  run <- function(eta) .liberation_objective_tape_eta_metropolis(
    lapply(evaluators, function(evaluator) evaluator$objective_tape$pointer),
    make_points(eta), length(parameters$theta) + seq_len(ncol(eta)), eta,
    proposal_roots, normals, log_uniforms, as.integer(mcmc_steps),
    as.numeric(step_scale)
  )
  result <- run(initial_eta)
  retaped <- ode_guard && any(vapply(seq_along(evaluators), function(subject) {
    evaluators[[subject]]$ensure_valid_tapes(
      parameters$theta, parameters$sigma, parameters$omega,
      result$eta[subject, ]
    )
  }, logical(1)))
  if (retaped) result <- run(initial_eta)
  result
}

.nm_saem_metropolis <- function(context, parameters, eta, mcmc_steps,
                                step_scale) {
  if (!context$n_eta) return(list(eta = eta, accepted = 0L, attempted = 0L))
  roots <- lapply(context$subjects, function(evaluator) {
    t(chol(.nm_effect_covariance(
      context$model, evaluator$data, parameters$omega
    )))
  })
  normals <- matrix(
    stats::rnorm(context$n_subjects * mcmc_steps * context$n_eta),
    context$n_subjects * mcmc_steps, context$n_eta
  )
  log_uniforms <- log(stats::runif(context$n_subjects * mcmc_steps))
  if (is.null(context$parallel)) {
    return(.nm_saem_metropolis_chunk(
      context$subjects, parameters, eta, roots, normals, log_uniforms,
      mcmc_steps, step_scale
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
  pieces <- parallel::clusterApply(
    context$parallel$cluster, seq_along(chunks),
    function(index, parameters, eta_chunks, root_chunks, normal_chunks,
             uniform_chunks, mcmc_steps, step_scale) {
      evaluators <- get(".liber_parallel_subjects", envir = .GlobalEnv)
      sampler <- get(".nm_saem_metropolis_chunk", envir = asNamespace("LibeRation"))
      sampler(
        evaluators, parameters, eta_chunks[[index]], root_chunks[[index]],
        normal_chunks[[index]], uniform_chunks[[index]], mcmc_steps, step_scale
      )
    }, parameters = parameters, eta_chunks = eta_chunks,
    root_chunks = root_chunks, normal_chunks = normal_chunks,
    uniform_chunks = uniform_chunks, mcmc_steps = mcmc_steps,
    step_scale = step_scale
  )
  list(
    eta = do.call(rbind, lapply(pieces, `[[`, "eta")),
    accepted = sum(vapply(pieces, `[[`, integer(1), "accepted")),
    attempted = sum(vapply(pieces, `[[`, integer(1), "attempted"))
  )
}

.nm_est_saem <- function(context, map, maxit, tolerance, trace,
                         n_iter = 200L, burn = NULL, mcmc_steps = 2L,
                         step_scale = 0.5, sa_power = 0.7,
                         mstep_maxit = 20L, seed = 20260713L,
                         print_every = 0L, adapt_proposal = TRUE,
                         target_acceptance = 0.3, closed_form_sigma = TRUE,
                         optimizer_backend = "auto", initial_eta = NULL) {
  n_iter <- as.integer(n_iter)
  burn <- as.integer(burn %||% floor(n_iter / 3))
  mcmc_steps <- as.integer(mcmc_steps)
  if (n_iter < 2L || burn < 0L || burn >= n_iter || mcmc_steps < 1L) {
    .nm_stop("SAEM requires n_iter >= 2, 0 <= burn < n_iter, and mcmc_steps >= 1.")
  }
  if (!is.finite(step_scale) || step_scale <= 0 ||
      !is.finite(target_acceptance) || target_acceptance <= 0 ||
      target_acceptance >= 1) {
    .nm_stop("SAEM proposal scale must be positive and target acceptance must lie in (0, 1).")
  }
  set.seed(seed)
  parameters <- map$decode(map$start)
  eta <- initial_eta %||% matrix(0, context$n_subjects, context$n_eta)
  accepted <- attempted <- 0L
  objective_trace <- numeric(n_iter)
  acceptance_trace <- numeric(n_iter)
  step_scale_trace <- numeric(n_iter)
  mstep_objective_evaluations <- 0L
  mstep_gradient_evaluations <- 0L
  mstep_iterations <- 0L
  mstep_elapsed <- 0
  mstep_backend <- "unknown"
  for (iteration in seq_len(n_iter)) {
    if (context$n_eta) {
      sampled <- .nm_saem_metropolis(
        context, parameters, eta, mcmc_steps, step_scale
      )
      eta <- sampled$eta
      accepted <- accepted + sampled$accepted
      attempted <- attempted + sampled$attempted
      acceptance_trace[[iteration]] <- sampled$accepted / max(sampled$attempted, 1L)
      if (isTRUE(adapt_proposal) && iteration <= burn) {
        gain <- min(0.1, 1 / sqrt(iteration))
        step_scale <- step_scale * exp(
          gain * (acceptance_trace[[iteration]] - target_acceptance)
        )
      }
    }
    step_scale_trace[[iteration]] <- step_scale
    gamma <- if (iteration <= burn) 1 else (iteration - burn)^(-sa_power)
    simple_sigma <- isTRUE(closed_form_sigma) &&
      context$model$LIK_CONFIG$error %in% c("additive", "proportional", "exponential")
    mstep_model <- context$model
    mstep_model$THETAS$Value <- parameters$theta
    mstep_model$SIGMAS$Value <- parameters$sigma
    mstep_model$OMEGAS$Value <- parameters$omega
    if (length(map$omega_free)) mstep_model$OMEGAS$FIX[] <- TRUE
    if (simple_sigma && length(map$sigma_free)) mstep_model$SIGMAS$FIX[] <- TRUE
    iteration_map <- .nm_outer_map(mstep_model)
    conditional <- function(candidate) .nm_saem_conditional(context, candidate, eta)
    conditional_gradient <- function(candidate) {
      .nm_saem_conditional_gradient(context, iteration_map, candidate, eta)
    }
    maximized <- .nm_outer_optim(
      iteration_map, conditional, min(as.integer(mstep_maxit), as.integer(maxit)),
      tolerance, if (trace > 1L) trace else 0L,
      gradient = conditional_gradient, optimizer_backend = optimizer_backend
    )
    mstep_objective_evaluations <- mstep_objective_evaluations +
      as.integer(maximized$objective_evaluations %||% 0L)
    mstep_gradient_evaluations <- mstep_gradient_evaluations +
      as.integer(maximized$gradient_evaluations %||% 0L)
    mstep_iterations <- mstep_iterations + as.integer(maximized$iterations %||% 0L)
    mstep_elapsed <- mstep_elapsed + as.numeric(maximized$elapsed_seconds %||% 0)
    mstep_backend <- maximized$backend %||% mstep_backend
    candidate <- iteration_map$decode(maximized$par)
    sigma_sufficient <- if (simple_sigma) {
      .nm_saem_sigma_sufficient(context, candidate, eta)
    } else NULL
    if (!is.null(sigma_sufficient) && length(map$sigma_free)) {
      candidate$sigma[map$sigma_free] <- sigma_sufficient[map$sigma_free]
    }
    parameters$theta[map$theta_free] <-
      (1 - gamma) * parameters$theta[map$theta_free] +
      gamma * candidate$theta[map$theta_free]
    parameters$sigma[map$sigma_free] <-
      (1 - gamma) * parameters$sigma[map$sigma_free] +
      gamma * candidate$sigma[map$sigma_free]
    if (length(map$omega_free) && context$n_eta) {
      sufficient <- .nm_saem_omega_sufficient(context, eta)
      parameters$omega[map$omega_free] <-
        (1 - gamma) * parameters$omega[map$omega_free] +
        gamma * sufficient[map$omega_free]
    }
    objective_trace[[iteration]] <- .nm_saem_conditional(context, parameters, eta)
    if (print_every > 0L && iteration %% print_every == 0L && length(map$start)) {
      point <- map$encode(parameters)
      objective_at <- function(value) .nm_saem_conditional(context, map$decode(value), eta)
      gradient_at <- function(value) {
        .nm_saem_conditional_gradient(context, map, map$decode(value), eta)
      }
      .nm_log_gradient(
        iteration, objective_at, point, map, objective_trace[[iteration]],
        gradient_function = gradient_at
      )
    }
  }
  modes <- lapply(seq_len(context$n_subjects), function(subject) {
    list(par = eta[subject, ], convergence = 0L, jitter = 0)
  })
  optimizer <- list(
    convergence = 0L, message = "SAEM iterations completed",
    counts = c(`function` = mstep_objective_evaluations,
               gradient = mstep_gradient_evaluations),
    iterations = n_iter, objective_evaluations = mstep_objective_evaluations,
    gradient_evaluations = mstep_gradient_evaluations,
    mstep_iterations = mstep_iterations, elapsed_seconds = mstep_elapsed,
    backend = paste0("saem+", mstep_backend)
  )
  fit <- .nm_fit_result(
    context, "SAEM", parameters, tail(objective_trace, 1), modes, optimizer,
    diagnostics = list(
      objective_trace = objective_trace, acceptance = accepted / max(attempted, 1L),
      acceptance_trace = acceptance_trace, step_scale_trace = step_scale_trace,
      final_step_scale = step_scale, target_acceptance = target_acceptance,
      adaptive_proposal = isTRUE(adapt_proposal),
      closed_form_sigma = simple_sigma,
      closed_form_omega = length(map$omega_free) > 0L,
      n_iter = n_iter, burn = burn, seed = seed,
      population_gradient = "exact conditional CppAD gradient"
    )
  )
  fit
}

.nm_bayes_state <- function(map, outer, eta) {
  list(outer = outer, parameters = map$decode(outer), eta = eta)
}

.nm_est_bayes <- function(context, map, tolerance,
                          n_burn = 500L, n_sample = 1000L, n_thin = 1L,
                          step_scale = 0.03, eta_step = 0.35,
                          seed = 20260713L, adapt = TRUE,
                          print_every = 0L) {
  n_burn <- as.integer(n_burn)
  n_sample <- as.integer(n_sample)
  n_thin <- as.integer(n_thin)
  if (n_burn < 0L || n_sample < 1L || n_thin < 1L) {
    .nm_stop("BAYES requires n_burn >= 0, n_sample >= 1, and n_thin >= 1.")
  }
  set.seed(seed)
  full_tape <- context$engine$objective_tape(
    context$data, theta = context$model$THETAS$Value,
    eta = matrix(0, context$n_subjects, context$n_eta),
    sigma = context$model$SIGMAS$Value, omega = context$model$OMEGAS$Value
  )
  log_posterior <- function(state) {
    if (!map$in_bounds(state$outer)) return(-Inf)
    parameters <- state$parameters
    point <- c(parameters$theta, as.vector(t(state$eta)),
               parameters$sigma, parameters$omega)
    nll <- tryCatch(
      .liberation_objective_tape_eval(full_tape$pointer, point, FALSE, FALSE)$value,
      error = function(e) Inf
    )
    if (!is.finite(nll)) return(-Inf)
    jacobian <- map$log_jacobian(parameters)
    -0.5 * nll + .nm_log_prior(context$model, parameters) + jacobian
  }
  state <- .nm_bayes_state(
    map, map$start, matrix(0, context$n_subjects, context$n_eta)
  )
  current <- log_posterior(state)
  total_iterations <- n_burn + n_sample * n_thin
  kept <- vector("list", n_sample)
  accepted_outer <- attempted_outer <- accepted_eta <- attempted_eta <- 0L
  keep <- 0L
  for (iteration in seq_len(total_iterations)) {
    if (length(state$outer)) {
      proposal <- .nm_bayes_state(
        map, state$outer + stats::rnorm(length(state$outer), sd = step_scale), state$eta
      )
      proposed <- log_posterior(proposal)
      attempted_outer <- attempted_outer + 1L
      if (log(stats::runif(1)) < proposed - current) {
        state <- proposal
        current <- proposed
        accepted_outer <- accepted_outer + 1L
      }
    }
    if (context$n_eta) {
      for (subject in seq_len(context$n_subjects)) {
        proposal <- state
        root <- t(chol(.nm_effect_covariance(
          context$model, context$subjects[[subject]]$data, state$parameters$omega
        )))
        proposal$eta[subject, ] <- proposal$eta[subject, ] +
          as.vector(root %*% stats::rnorm(context$n_eta, sd = eta_step))
        proposed <- log_posterior(proposal)
        attempted_eta <- attempted_eta + 1L
        if (log(stats::runif(1)) < proposed - current) {
          state <- proposal
          current <- proposed
          accepted_eta <- accepted_eta + 1L
        }
      }
    }
    if (isTRUE(adapt) && iteration <= n_burn && iteration %% 50L == 0L) {
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
      point <- c(state$parameters$theta, as.vector(t(state$eta)),
                 state$parameters$sigma, state$parameters$omega)
      evaluated <- .liberation_objective_tape_eval(full_tape$pointer, point, TRUE, FALSE)
      population <- c(
        evaluated$gradient[seq_along(state$parameters$theta)],
        evaluated$gradient[length(state$parameters$theta) + length(state$eta) +
                             seq_along(state$parameters$sigma)],
        evaluated$gradient[length(state$parameters$theta) + length(state$eta) +
                             length(state$parameters$sigma) + seq_along(state$parameters$omega)]
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
  chain <- do.call(rbind, kept)
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
  final_state <- list(parameters = parameters, eta = eta)
  final_objective <- -2 * log_posterior(final_state)
  modes <- lapply(seq_len(context$n_subjects), function(subject) {
    list(par = eta[subject, ], convergence = 0L, jitter = 0)
  })
  optimizer <- list(
    convergence = 0L, message = "Bayesian sampling completed",
    counts = c(`function` = total_iterations, gradient = NA_integer_),
    iterations = total_iterations, objective_evaluations = total_iterations
  )
  fit <- .nm_fit_result(
    context, "BAYES", parameters, final_objective, modes, optimizer,
    diagnostics = list(
      outer_acceptance = accepted_outer / max(attempted_outer, 1L),
      eta_acceptance = accepted_eta / max(attempted_eta, 1L),
      n_burn = n_burn, n_sample = n_sample, n_thin = n_thin,
      seed = seed, final_step_scale = step_scale
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
      correlation = population_correlation
    )
  )
  fit
}

.nm_est_stochastic <- function(context, map, method, maxit, eta_maxit,
                               tolerance, trace, print_every = 0L,
                               optimizer_backend = "auto", initial_eta = NULL, ...) {
  controls <- list(...)
  if (method == "ITS") {
    return(.nm_est_its(context, map, maxit, eta_maxit, tolerance, trace,
                       print_every, optimizer_backend))
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
