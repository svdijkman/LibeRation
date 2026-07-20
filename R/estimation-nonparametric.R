.nm_np_loglik <- function(context, parameters, supports, gradient = FALSE) {
  supports <- as.matrix(supports)
  n_support <- nrow(supports)
  n_theta <- length(parameters$theta)
  n_sigma <- length(parameters$sigma)
  n_omega <- length(parameters$omega)
  population_positions <- c(
    seq_len(n_theta), n_theta + context$n_eta + seq_len(n_sigma),
    n_theta + context$n_eta + n_sigma + seq_len(n_omega)
  )
  values <- matrix(NA_real_, context$n_subjects, n_support)
  gradients <- if (isTRUE(gradient)) {
    array(0, c(context$n_subjects, n_support,
               n_theta + n_sigma + n_omega))
  } else NULL
  for (subject in seq_len(context$n_subjects)) {
    evaluator <- context$subjects[[subject]]
    covariance <- .nm_effect_covariance(context$model, evaluator$data, parameters$omega)
    positive <- .nm_positive_definite(covariance, "nonparametric reference OMEGA")
    inverse <- solve(positive$matrix)
    eta_prior_q <- context$n_eta * log(2 * pi) + positive$logdet +
      rowSums((supports %*% inverse) * supports)
    evaluated <- if (isTRUE(gradient)) {
      evaluator$objective_eta_batch(
        parameters$theta, supports, parameters$sigma, parameters$omega
      )
    } else list(value = evaluator$objective_eta_values(
      parameters$theta, supports, parameters$sigma, parameters$omega
    ))
    values[subject, ] <- -0.5 * (as.numeric(evaluated$value) - eta_prior_q)
    if (isTRUE(gradient)) {
      # The subtracted Gaussian ETA density has no THETA or SIGMA derivative.
      # OMEGA is fixed while fitting a discrete nonparametric distribution.
      gradients[subject, , ] <- evaluated$gradient[, population_positions, drop = FALSE]
    }
  }
  list(loglik = values, gradient = gradients)
}

.nm_np_responsibilities <- function(loglik, weights) {
  log_weights <- log(pmax(weights, .Machine$double.xmin))
  responsibilities <- matrix(0, nrow(loglik), ncol(loglik))
  marginal <- numeric(nrow(loglik))
  for (subject in seq_len(nrow(loglik))) {
    score <- log_weights + loglik[subject, ]
    normalizer <- .nm_log_sum_exp(score)
    marginal[[subject]] <- normalizer
    responsibilities[subject, ] <- exp(score - normalizer)
  }
  list(responsibilities = responsibilities, log_likelihood = sum(marginal))
}

.nm_np_weights <- function(loglik, initial = NULL, maxit = 1000L,
                           tolerance = 1e-8) {
  n_support <- ncol(loglik)
  weights <- if (is.null(initial) || length(initial) != n_support) {
    rep(1 / n_support, n_support)
  } else {
    initial <- pmax(as.numeric(initial), 0)
    initial / sum(initial)
  }
  history <- numeric()
  for (iteration in seq_len(as.integer(maxit))) {
    state <- .nm_np_responsibilities(loglik, weights)
    next_weights <- colMeans(state$responsibilities)
    next_weights <- pmax(next_weights, .Machine$double.eps)
    next_weights <- next_weights / sum(next_weights)
    history[[iteration]] <- state$log_likelihood
    if (max(abs(next_weights - weights)) <= tolerance) {
      weights <- next_weights
      break
    }
    weights <- next_weights
  }
  state <- .nm_np_responsibilities(loglik, weights)
  list(weights = weights, responsibilities = state$responsibilities,
       log_likelihood = state$log_likelihood, iterations = length(history),
       history = history)
}

.nm_np_unique_supports <- function(supports, tolerance = 1e-7) {
  supports <- as.matrix(supports)
  if (!nrow(supports)) return(supports)
  key <- apply(round(supports / tolerance) * tolerance, 1, paste, collapse = "|")
  supports[!duplicated(key), , drop = FALSE]
}

.nm_np_initial_supports <- function(context, parameters, supports = NULL,
                                    points = 25L, seed = 20260719L,
                                    eta_maxit = 100L, tolerance = 1e-6) {
  if (!context$n_eta) .nm_stop("NPML and NPAG require at least one ETA.")
  if (!is.null(supports)) {
    supports <- as.matrix(supports)
    storage.mode(supports) <- "double"
    if (ncol(supports) != context$n_eta || !nrow(supports) || any(!is.finite(supports))) {
      .nm_stop("`np_supports` must be a finite matrix with one column per ETA.")
    }
    return(.nm_np_unique_supports(supports))
  }
  modes <- .nm_subject_modes(
    context, parameters, maxit = eta_maxit, tolerance = tolerance,
    interaction = TRUE, exact_hessian = FALSE
  )
  mode_matrix <- do.call(rbind, lapply(modes, `[[`, "par"))
  candidate <- rbind(rep(0, context$n_eta), mode_matrix)
  candidate <- .nm_np_unique_supports(candidate)
  points <- max(1L, as.integer(points))
  if (nrow(candidate) > points) {
    selected <- unique(round(seq(1, nrow(candidate), length.out = points)))
    candidate <- candidate[selected, , drop = FALSE]
  }
  if (nrow(candidate) < points) {
    set.seed(seed)
    covariance <- .nm_effect_covariance(
      context$model, context$subjects[[1L]]$data, parameters$omega
    )
    root <- t(chol(.nm_positive_definite(covariance, "initial NP grid")$matrix))
    random <- matrix(stats::rnorm((points - nrow(candidate)) * context$n_eta),
                     ncol = context$n_eta) %*% t(root)
    candidate <- .nm_np_unique_supports(rbind(candidate, random))
  }
  candidate
}

.nm_np_population_fit <- function(context, model, map, parameters,
                                  supports, weights, maxit, tolerance,
                                  trace, print_every, optimizer_backend) {
  if (!length(map$start)) {
    evaluated <- .nm_np_loglik(context, parameters, supports, gradient = FALSE)
    state <- .nm_np_responsibilities(evaluated$loglik, weights)
    return(list(parameters = parameters, optimizer = list(
      par = numeric(), value = -2 * state$log_likelihood + .nm_prior_nll(model, parameters),
      convergence = 0L, message = "All population parameters fixed",
      counts = c(`function` = 1L, gradient = 0L), iterations = 0L,
      objective_evaluations = 1L, backend = "fixed-parameters"
    )))
  }
  objective <- function(candidate) {
    evaluated <- .nm_np_loglik(context, candidate, supports, gradient = FALSE)
    -2 * .nm_np_responsibilities(evaluated$loglik, weights)$log_likelihood +
      .nm_prior_nll(model, candidate)
  }
  gradient <- function(candidate) {
    evaluated <- .nm_np_loglik(context, candidate, supports, gradient = TRUE)
    state <- .nm_np_responsibilities(evaluated$loglik, weights)
    native <- numeric(
      length(candidate$theta) + length(candidate$sigma) + length(candidate$omega)
    )
    for (subject in seq_len(context$n_subjects)) {
      subject_gradient <- matrix(
        evaluated$gradient[subject, , , drop = FALSE],
        nrow = nrow(supports), ncol = length(native)
      )
      native <- native + colSums(
        subject_gradient * state$responsibilities[subject, ]
      )
    }
    native <- native + .nm_prior_nll_native_gradient(model, candidate)
    as.vector(native %*% map$jacobian(candidate))
  }
  local_map <- map
  local_map$start <- map$encode(parameters)
  optimizer <- .nm_outer_optim(
    local_map, objective, maxit = maxit, tolerance = tolerance, trace = trace,
    print_every = print_every, gradient = gradient,
    optimizer_backend = optimizer_backend
  )
  list(parameters = map$decode(optimizer$par), optimizer = optimizer)
}

.nm_np_prune <- function(supports, weights, minimum_weight, max_support) {
  keep <- which(weights >= minimum_weight)
  if (!length(keep)) keep <- which.max(weights)
  if (length(keep) > max_support) keep <- keep[order(weights[keep], decreasing = TRUE)[seq_len(max_support)]]
  supports <- supports[keep, , drop = FALSE]
  weights <- weights[keep]
  weights <- weights / sum(weights)
  list(supports = supports, weights = weights)
}

.nm_np_expand <- function(supports, root, step, max_candidates) {
  candidates <- list(supports)
  for (axis in seq_len(ncol(supports))) {
    delta <- step * root[, axis]
    candidates[[length(candidates) + 1L]] <- sweep(supports, 2, delta, "+")
    candidates[[length(candidates) + 1L]] <- sweep(supports, 2, delta, "-")
  }
  result <- .nm_np_unique_supports(do.call(rbind, candidates))
  if (nrow(result) > max_candidates) {
    selected <- unique(round(seq(1, nrow(result), length.out = max_candidates)))
    result <- result[selected, , drop = FALSE]
  }
  result
}

.nm_np_omega <- function(context, supports, weights, fallback) {
  dimension <- context$model$n_eta
  selected <- supports[, seq_len(min(dimension, ncol(supports))), drop = FALSE]
  if (ncol(selected) < dimension) return(fallback)
  mean <- colSums(selected * weights)
  centered <- sweep(selected, 2, mean, "-")
  covariance <- crossprod(centered * sqrt(weights))
  covariance <- covariance + diag(1e-8, dimension)
  vapply(seq_len(nrow(context$model$OMEGAS)), function(index) {
    covariance[context$model$OMEGAS$ROW[[index]], context$model$OMEGAS$COL[[index]]]
  }, numeric(1))
}

.nm_est_nonparametric <- function(context, method = c("NPML", "NPAG"),
                                  maxit = 100L, tolerance = 1e-6, trace = 0L,
                                  print_every = 0L, optimizer_backend = "auto",
                                  np_supports = NULL, np_points = 25L,
                                  np_max_support = 100L,
                                  np_min_weight = 1e-5,
                                  np_weight_maxit = 1000L,
                                  np_cycles = 3L, np_grid_step = 1,
                                  np_grid_decay = 0.5,
                                  np_max_candidates = 500L,
                                  np_estimate_population = TRUE,
                                  eta_maxit = 100L,
                                  seed = 20260719L, ...) {
  method <- match.arg(method)
  if (!context$n_eta) .nm_stop(method, " requires at least one ETA.")
  priors <- context$model$LIK_CONFIG$priors
  if (!is.null(priors) && nrow(priors) && any(startsWith(priors$parameter, "OMEGA"))) {
    .nm_stop(
      method, " replaces OMEGA with a discrete support distribution; OMEGA priors are not applicable."
    )
  }
  numeric_controls <- c(
    np_points = np_points, np_max_support = np_max_support,
    np_weight_maxit = np_weight_maxit, np_cycles = np_cycles,
    np_grid_step = np_grid_step, np_grid_decay = np_grid_decay,
    np_max_candidates = np_max_candidates
  )
  if (any(!is.finite(numeric_controls)) || any(numeric_controls <= 0) ||
      !is.finite(np_min_weight) || np_min_weight < 0 || np_min_weight >= 1) {
    .nm_stop("Nonparametric grid, cycle, support, and weight controls must be positive and finite.")
  }
  np_model <- context$model
  np_model$OMEGAS$FIX[] <- TRUE
  map <- .nm_outer_map(np_model)
  parameters <- list(
    theta = context$model$THETAS$Value,
    sigma = context$model$SIGMAS$Value,
    omega = context$model$OMEGAS$Value
  )
  supports <- .nm_np_initial_supports(
    context, parameters, np_supports, np_points, seed, eta_maxit, tolerance
  )
  weights <- rep(1 / nrow(supports), nrow(supports))
  cycles <- if (method == "NPAG") max(1L, as.integer(np_cycles)) else
    max(1L, as.integer(np_cycles))
  history <- data.frame(
    cycle = integer(), supports = integer(), log_likelihood = numeric(),
    grid_step = numeric(), stringsAsFactors = FALSE
  )
  optimizer <- NULL
  step <- as.numeric(np_grid_step)
  for (cycle in seq_len(cycles)) {
    evaluated <- .nm_np_loglik(context, parameters, supports, gradient = FALSE)
    weight_fit <- .nm_np_weights(
      evaluated$loglik, weights, np_weight_maxit, tolerance = tolerance * 0.1
    )
    weights <- weight_fit$weights
    pruned <- .nm_np_prune(
      supports, weights, np_min_weight, as.integer(np_max_support)
    )
    supports <- pruned$supports; weights <- pruned$weights
    if (isTRUE(np_estimate_population)) {
      population <- .nm_np_population_fit(
        context, np_model, map, parameters, supports, weights,
        maxit = max(1L, min(as.integer(maxit), 25L)), tolerance = tolerance,
        trace = trace, print_every = print_every,
        optimizer_backend = optimizer_backend
      )
      parameters <- population$parameters
      optimizer <- population$optimizer
    }
    evaluated <- .nm_np_loglik(context, parameters, supports, gradient = FALSE)
    weight_fit <- .nm_np_weights(
      evaluated$loglik, weights, np_weight_maxit, tolerance = tolerance * 0.1
    )
    weights <- weight_fit$weights
    if (method == "NPAG" && cycle < cycles && step > tolerance) {
      covariance <- .nm_effect_covariance(
        context$model, context$subjects[[1L]]$data, parameters$omega
      )
      root <- t(chol(.nm_positive_definite(covariance, "NPAG grid scale")$matrix))
      candidate <- .nm_np_expand(
        supports, root, step, as.integer(np_max_candidates)
      )
      candidate_eval <- .nm_np_loglik(context, parameters, candidate, gradient = FALSE)
      candidate_fit <- .nm_np_weights(
        candidate_eval$loglik, NULL, np_weight_maxit, tolerance = tolerance * 0.1
      )
      if (candidate_fit$log_likelihood >= weight_fit$log_likelihood - tolerance) {
        supports <- candidate
        weights <- candidate_fit$weights
        pruned <- .nm_np_prune(
          supports, weights, np_min_weight, as.integer(np_max_support)
        )
        supports <- pruned$supports; weights <- pruned$weights
      }
      step <- step * np_grid_decay
    }
    final_eval <- .nm_np_loglik(context, parameters, supports, gradient = FALSE)
    final_weights <- .nm_np_weights(
      final_eval$loglik, weights, np_weight_maxit, tolerance = tolerance * 0.1
    )
    weights <- final_weights$weights
    history <- rbind(history, data.frame(
      cycle = cycle, supports = nrow(supports),
      log_likelihood = final_weights$log_likelihood,
      grid_step = if (method == "NPAG") step else NA_real_
    ))
    if (print_every > 0L && (cycle == 1L || cycle %% print_every == 0L)) {
      cat(sprintf(
        "[LibeRation] %s CYCLE %d SUPPORTS %d -2LOGLIK %.10g\n",
        method, cycle, nrow(supports), -2 * final_weights$log_likelihood
      ))
      try(flush(stdout()), silent = TRUE)
    }
    if (nrow(history) > 1L &&
        abs(diff(tail(history$log_likelihood, 2L))) <= tolerance &&
        (method == "NPML" || step <= tolerance)) break
  }
  final_eval <- .nm_np_loglik(context, parameters, supports, gradient = FALSE)
  final_weights <- .nm_np_weights(
    final_eval$loglik, weights, np_weight_maxit, tolerance = tolerance * 0.1
  )
  weights <- final_weights$weights
  responsibilities <- final_weights$responsibilities
  eta <- responsibilities %*% supports
  parameters$omega <- .nm_np_omega(context, supports, weights, parameters$omega)
  modes <- lapply(seq_len(context$n_subjects), function(subject) {
    list(par = eta[subject, ], convergence = 0L, jitter = 0)
  })
  if (is.null(optimizer)) optimizer <- list(
    convergence = 0L, message = paste(method, "weight optimization completed"),
    counts = c(`function` = nrow(history), gradient = 0L),
    iterations = nrow(history), objective_evaluations = nrow(history),
    backend = if (method == "NPAG") "adaptive-grid-em" else "fixed-support-em"
  )
  optimizer$convergence <- 0L
  optimizer$message <- paste(method, "nonparametric estimation completed")
  optimizer$iterations <- nrow(history)
  objective <- -2 * final_weights$log_likelihood + .nm_prior_nll(np_model, parameters)
  fit <- .nm_fit_result(
    context, method, parameters, objective, modes, optimizer,
    diagnostics = list(
      nonparametric = list(
        distribution = "discrete ETA support",
        adaptive_grid = method == "NPAG", cycles = nrow(history),
        support_count = nrow(supports), weight_iterations = final_weights$iterations,
        history = history
      )
    )
  )
  colnames(supports) <- paste0("ETA", seq_len(context$n_eta))
  fit$nonparametric <- list(
    supports = supports, weights = weights,
    posterior_probabilities = responsibilities,
    posterior_eta = eta, log_likelihood = final_weights$log_likelihood,
    history = history,
    interpretation = "OMEGA is the weighted support covariance; inference is based on the discrete distribution."
  )
  fit
}
