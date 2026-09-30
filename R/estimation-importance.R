# Importance-sampling and IMP estimation.
# Split from estimation-stochastic.R as a behaviour-preserving source move.

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

