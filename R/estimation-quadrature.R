# Gaussian-quadrature estimation.
# Split from estimation-stochastic.R as a behaviour-preserving source move.

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

