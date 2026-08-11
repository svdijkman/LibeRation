.nm_hmc_bound_map <- function(map) {
  lower <- map$lower
  upper <- map$upper
  inverse <- function(x) {
    if (!length(x)) return(numeric())
    result <- numeric(length(x))
    for (i in seq_along(x)) {
      lo <- lower[[i]]
      hi <- upper[[i]]
      if (is.finite(lo) && is.finite(hi)) {
        width <- hi - lo
        p <- min(max((x[[i]] - lo) / width, 1e-8), 1 - 1e-8)
        result[[i]] <- stats::qlogis(p)
      } else if (is.finite(lo)) {
        result[[i]] <- log(max(x[[i]] - lo, 1e-8))
      } else if (is.finite(hi)) {
        result[[i]] <- log(max(hi - x[[i]], 1e-8))
      } else result[[i]] <- x[[i]]
    }
    result
  }
  forward <- function(q) {
    if (!length(q)) return(list(value = numeric(), derivative = numeric(),
                                log_jacobian = 0, log_jacobian_gradient = numeric()))
    value <- derivative <- log_gradient <- numeric(length(q))
    log_jacobian <- 0
    for (i in seq_along(q)) {
      lo <- lower[[i]]
      hi <- upper[[i]]
      if (is.finite(lo) && is.finite(hi)) {
        p <- stats::plogis(q[[i]])
        value[[i]] <- lo + (hi - lo) * p
        derivative[[i]] <- (hi - lo) * p * (1 - p)
        log_jacobian <- log_jacobian + log(hi - lo) +
          stats::plogis(q[[i]], log.p = TRUE) +
          stats::plogis(-q[[i]], log.p = TRUE)
        log_gradient[[i]] <- 1 - 2 * p
      } else if (is.finite(lo)) {
        scale <- exp(q[[i]])
        value[[i]] <- lo + scale
        derivative[[i]] <- scale
        log_jacobian <- log_jacobian + q[[i]]
        log_gradient[[i]] <- 1
      } else if (is.finite(hi)) {
        scale <- exp(q[[i]])
        value[[i]] <- hi - scale
        derivative[[i]] <- -scale
        log_jacobian <- log_jacobian + q[[i]]
        log_gradient[[i]] <- 1
      } else {
        value[[i]] <- q[[i]]
        derivative[[i]] <- 1
      }
    }
    list(value = value, derivative = derivative,
         log_jacobian = log_jacobian,
         log_jacobian_gradient = log_gradient)
  }
  list(inverse = inverse, forward = forward)
}

.nm_hmc_geometry_plan <- function(context, map,
                                  geometry = c("auto", "centered",
                                               "mu_noncentered")) {
  geometry <- match.arg(geometry)
  model <- context$model
  mu_eta <- sort(unique(as.integer(model$MU$ETA %||% integer())))
  expected_eta <- seq_len(context$n_eta)
  covariance <- tryCatch(.nm_omega_matrix(model), error = function(error) NULL)
  positive_definite <- !is.null(covariance) && nrow(covariance) == context$n_eta &&
    !inherits(tryCatch(chol(covariance), error = identity), "error")
  reason <- if (!.nm_liber_optimized(context)) {
    "MU non-centred geometry is restricted to liber_optimized mode"
  } else if (!context$n_eta) {
    "the model has no random effects"
  } else if (!identical(context$n_eta, model$n_eta)) {
    "expanded IOV random effects require centred geometry"
  } else if (!is.null(model$RE_CONFIG)) {
    "general random-effect mappings require centred geometry"
  } else if (!identical(mu_eta, expected_eta)) {
    "MU references do not cover every ETA"
  } else if (!positive_definite) {
    "OMEGA is not positive definite"
  } else NULL
  eligible <- is.null(reason)
  if (identical(geometry, "mu_noncentered") && !eligible) {
    .nm_stop("`geometry = \"mu_noncentered\"` is unavailable because ",
             reason, ".")
  }
  enabled <- eligible && !identical(geometry, "centered")
  list(
    requested = geometry, enabled = enabled,
    resolved = if (enabled) "mu_noncentered" else "centered",
    eligible = eligible,
    reason = if (enabled) NULL else reason %||%
      if (identical(geometry, "centered")) "centred geometry was requested" else NULL
  )
}

.nm_hmc_eta_transform <- function(context, map, parameters, encoded, plan) {
  n_eta <- context$n_eta
  n_outer <- length(encoded)
  factor <- t(chol(.nm_omega_matrix(context$model, parameters$omega)))
  derivatives <- replicate(
    n_outer, matrix(0, n_eta, n_eta), simplify = FALSE
  )
  omega_offset <- length(map$theta_free) + length(map$sigma_free)
  if (isTRUE(map$omega_full) && length(map$omega_free)) {
    for (entry in seq_len(nrow(context$model$OMEGAS))) {
      row <- context$model$OMEGAS$ROW[[entry]]
      column <- context$model$OMEGAS$COL[[entry]]
      derivatives[[omega_offset + entry]][row, column] <-
        if (row == column) factor[row, column] else 1
    }
  } else if (length(map$omega_free)) {
    for (entry in seq_along(map$omega_free)) {
      native <- map$omega_free[[entry]]
      row <- context$model$OMEGAS$ROW[[native]]
      column <- context$model$OMEGAS$COL[[native]]
      if (row != column) {
        .nm_stop("Non-centred diagonal OMEGA mapping encountered an off-diagonal entry.")
      }
      derivatives[[omega_offset + entry]][row, column] <- factor[row, column] / 2
    }
  }
  log_gradient <- vapply(derivatives, function(derivative) {
    context$n_subjects * sum(diag(solve(factor, derivative)))
  }, numeric(1))
  list(
    factor = factor, derivatives = derivatives,
    log_jacobian = context$n_subjects * sum(log(diag(factor))),
    log_jacobian_gradient = log_gradient
  )
}

.nm_hmc_target <- function(context, map, geometry_plan = NULL) {
  geometry_plan <- geometry_plan %||% .nm_hmc_geometry_plan(
    context, map, "centered"
  )
  full_tape <- context$engine$objective_tape(
    context$data, theta = context$model$THETAS$Value,
    eta = matrix(0, context$n_subjects, context$n_eta),
    sigma = context$model$SIGMAS$Value, omega = context$model$OMEGAS$Value
  )
  bounds <- .nm_hmc_bound_map(map)
  n_outer <- length(map$start)
  n_eta_total <- context$n_subjects * context$n_eta
  n_theta <- nrow(context$model$THETAS)
  n_sigma <- nrow(context$model$SIGMAS)
  n_omega <- nrow(context$model$OMEGAS)
  eta_positions <- n_theta + seq_len(n_eta_total)
  population_positions <- c(
    seq_len(n_theta),
    n_theta + n_eta_total + seq_len(n_sigma),
    n_theta + n_eta_total + n_sigma + seq_len(n_omega)
  )
  initial <- c(bounds$inverse(map$start), rep(0, n_eta_total))
  evaluate <- function(q) {
    if (length(q) != n_outer + n_eta_total || any(!is.finite(q))) {
      return(list(logp = -Inf, gradient = rep(NA_real_, length(q))))
    }
    transformed <- bounds$forward(q[seq_len(n_outer)])
    parameters <- tryCatch(map$decode(transformed$value), error = function(e) NULL)
    if (is.null(parameters)) {
      return(list(logp = -Inf, gradient = rep(NA_real_, length(q))))
    }
    latent_eta <- if (n_eta_total) matrix(
      q[n_outer + seq_len(n_eta_total)], context$n_subjects, context$n_eta,
      byrow = TRUE
    ) else matrix(numeric(), context$n_subjects, 0L)
    eta_transform <- if (isTRUE(geometry_plan$enabled)) {
      .nm_hmc_eta_transform(
        context, map, parameters, transformed$value, geometry_plan
      )
    } else NULL
    eta <- if (!is.null(eta_transform)) {
      latent_eta %*% t(eta_transform$factor)
    } else latent_eta
    point <- c(parameters$theta, as.vector(t(eta)), parameters$sigma, parameters$omega)
    evaluated <- tryCatch(
      .liberation_objective_tape_eval(full_tape$pointer, point, TRUE, FALSE),
      error = function(e) NULL
    )
    if (is.null(evaluated) || !is.finite(evaluated$value) ||
        any(!is.finite(evaluated$gradient))) {
      return(list(logp = -Inf, gradient = rep(NA_real_, length(q))))
    }
    prior <- .nm_log_prior(context$model, parameters)
    if (!is.finite(prior)) {
      return(list(logp = -Inf, gradient = rep(NA_real_, length(q))))
    }
    native_gradient <- -0.5 * as.numeric(evaluated$gradient[population_positions]) -
      0.5 * .nm_prior_nll_native_gradient(context$model, parameters)
    outer_gradient <- as.vector(native_gradient %*% map$jacobian(parameters)) +
      map$log_jacobian_gradient(parameters)
    eta_gradient_matrix <- if (n_eta_total) matrix(
      -0.5 * as.numeric(evaluated$gradient[eta_positions]),
      context$n_subjects, context$n_eta, byrow = TRUE
    ) else matrix(numeric(), context$n_subjects, 0L)
    if (!is.null(eta_transform)) {
      for (index in seq_len(n_outer)) {
        derivative <- eta_transform$derivatives[[index]]
        if (any(derivative != 0)) {
          outer_gradient[[index]] <- outer_gradient[[index]] + sum(
            eta_gradient_matrix * (latent_eta %*% t(derivative))
          )
        }
      }
      outer_gradient <- outer_gradient + eta_transform$log_jacobian_gradient
    }
    outer_gradient <- outer_gradient * transformed$derivative +
      transformed$log_jacobian_gradient
    eta_gradient <- if (n_eta_total) {
      as.vector(t(if (is.null(eta_transform)) eta_gradient_matrix else
        eta_gradient_matrix %*% eta_transform$factor))
    } else numeric()
    list(
      logp = -0.5 * evaluated$value + prior + map$log_jacobian(parameters) +
        transformed$log_jacobian + (eta_transform$log_jacobian %||% 0),
      gradient = c(outer_gradient, eta_gradient), parameters = parameters, eta = eta,
      outer = transformed$value
    )
  }
  list(evaluate = evaluate, initial = initial, bounds = bounds, map = map,
       n_outer = n_outer, n_eta_total = n_eta_total, full_tape = full_tape,
       geometry = geometry_plan)
}

.nm_hmc_cpp_config <- function(context, map, target) {
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
    n_subjects = as.integer(context$n_subjects),
    n_eta = as.integer(context$n_eta),
    noncentered = isTRUE(target$geometry$enabled),
    lower = as.numeric(map$lower),
    upper = as.numeric(map$upper),
    initial = as.numeric(target$initial),
    prior_index = priors$index,
    prior_family = priors$family,
    prior_mean = priors$mean,
    prior_sd = priors$sd,
    prior_shape = priors$shape,
    prior_rate = priors$rate,
    output_columns = as.integer(
      nrow(context$model$THETAS) + nrow(context$model$SIGMAS) +
        nrow(context$model$OMEGAS) +
        context$n_subjects * context$n_eta + 1L
    )
  )
}

.nm_hmc_native_target_eval <- function(target, context, map, q) {
  .liberation_hmc_target_eval(
    target$full_tape$pointer, as.numeric(q),
    .nm_hmc_cpp_config(context, map, target)
  )
}

.nm_hmc_metric_blocks <- function(dimension, context, blocks = NULL) {
  dimension <- as.integer(dimension)
  if (!is.null(blocks)) {
    if (!is.list(blocks)) .nm_stop("`hmc_metric_blocks` must be a list of integer index vectors.")
    blocks <- lapply(blocks, function(index) {
      index <- unique(as.integer(index))
      if (!length(index) || anyNA(index) || any(index < 1L | index > dimension)) {
        .nm_stop("Every HMC metric block must contain valid one-based parameter indices.")
      }
      index
    })
    used <- unlist(blocks, use.names = FALSE)
    if (anyDuplicated(used)) .nm_stop("HMC metric blocks must not overlap.")
    missing <- setdiff(seq_len(dimension), used)
    return(c(blocks, lapply(missing, function(index) index)))
  }
  n_eta_total <- as.integer(context$n_subjects * context$n_eta)
  n_outer <- dimension - n_eta_total
  result <- if (n_outer > 0L) list(seq_len(n_outer)) else list()
  if (context$n_eta > 0L) {
    result <- c(result, lapply(seq_len(context$n_subjects), function(subject) {
      n_outer + (subject - 1L) * context$n_eta + seq_len(context$n_eta)
    }))
  }
  result
}

.nm_hmc_metric <- function(dimension, kind = "diagonal", blocks = list()) {
  dimension <- as.integer(dimension)
  mass <- diag(1, dimension)
  structure(list(
    kind = kind, mass = mass, inverse_mass = mass, root = mass,
    blocks = blocks
  ), class = "nm_hmc_metric")
}

.nm_hmc_metric_from_mass <- function(mass, kind, blocks = list()) {
  mass <- (as.matrix(mass) + t(as.matrix(mass))) / 2
  decomposition <- eigen(mass, symmetric = TRUE)
  values <- pmin(pmax(decomposition$values, 1e-3), 1e3)
  mass <- decomposition$vectors %*% (values * t(decomposition$vectors))
  inverse_mass <- decomposition$vectors %*% ((1 / values) * t(decomposition$vectors))
  root <- t(chol(mass))
  structure(list(
    kind = kind, mass = mass, inverse_mass = inverse_mass, root = root,
    blocks = blocks
  ), class = "nm_hmc_metric")
}

.nm_hmc_metric_sample <- function(metric) {
  as.numeric(metric$root %*% stats::rnorm(nrow(metric$mass)))
}

.nm_hmc_velocity <- function(momentum, metric) {
  if (is.numeric(metric) && is.null(dim(metric))) return(momentum / metric)
  as.numeric(metric$inverse_mass %*% momentum)
}

.nm_hmc_kinetic <- function(momentum, metric) {
  as.numeric(crossprod(momentum, .nm_hmc_velocity(momentum, metric)) / 2)
}

.nm_hmc_leapfrog <- function(q, momentum, gradient, epsilon, metric, target) {
  next_momentum <- momentum + 0.5 * epsilon * gradient
  next_q <- q + epsilon * .nm_hmc_velocity(next_momentum, metric)
  evaluated <- target$evaluate(next_q)
  if (!is.finite(evaluated$logp) || any(!is.finite(evaluated$gradient))) {
    return(list(q = next_q, momentum = next_momentum,
                evaluated = evaluated, valid = FALSE))
  }
  next_momentum <- next_momentum + 0.5 * epsilon * evaluated$gradient
  list(q = next_q, momentum = next_momentum, evaluated = evaluated, valid = TRUE)
}

.nm_hmc_find_step <- function(q, evaluated, metric, target) {
  epsilon <- 1
  momentum <- .nm_hmc_metric_sample(metric)
  proposal <- .nm_hmc_leapfrog(q, momentum, evaluated$gradient, epsilon, metric, target)
  log_accept <- if (proposal$valid) {
    proposal$evaluated$logp - .nm_hmc_kinetic(proposal$momentum, metric) -
      evaluated$logp + .nm_hmc_kinetic(momentum, metric)
  } else -Inf
  direction <- if (is.finite(log_accept) && log_accept > log(0.5)) 1 else -1
  for (iteration in seq_len(20L)) {
    candidate <- epsilon * if (direction > 0) 2 else 0.5
    proposal <- .nm_hmc_leapfrog(q, momentum, evaluated$gradient, candidate, metric, target)
    candidate_accept <- if (proposal$valid) {
      proposal$evaluated$logp - .nm_hmc_kinetic(proposal$momentum, metric) -
        evaluated$logp + .nm_hmc_kinetic(momentum, metric)
    } else -Inf
    continue <- if (direction > 0) candidate_accept > log(0.5) else candidate_accept < log(0.5)
    epsilon <- candidate
    if (!continue || epsilon < 1e-8 || epsilon > 1e2) break
  }
  min(max(epsilon, 1e-8), 1e2)
}

.nm_dual_average <- function(initial, target_acceptance) {
  environment <- new.env(parent = emptyenv())
  environment$mu <- log(10 * initial)
  environment$log_step <- log(initial)
  environment$log_step_bar <- log(initial)
  environment$hbar <- 0
  environment$iteration <- 0L
  environment$update <- function(acceptance) {
    environment$iteration <- environment$iteration + 1L
    t <- environment$iteration
    environment$hbar <- (1 - 1 / (t + 10)) * environment$hbar +
      (target_acceptance - acceptance) / (t + 10)
    environment$log_step <- environment$mu - sqrt(t) / 0.05 * environment$hbar
    weight <- t^-0.75
    environment$log_step_bar <- weight * environment$log_step +
      (1 - weight) * environment$log_step_bar
    exp(environment$log_step)
  }
  environment$final <- function() exp(environment$log_step_bar)
  environment
}

.nm_hmc_transition <- function(q, evaluated, step_size, metric, n_leapfrog,
                               target, divergence_threshold) {
  momentum <- .nm_hmc_metric_sample(metric)
  initial_momentum <- momentum
  energy <- -evaluated$logp + .nm_hmc_kinetic(initial_momentum, metric)
  proposal_q <- q
  proposal <- evaluated
  valid <- TRUE
  for (step in seq_len(n_leapfrog)) {
    moved <- .nm_hmc_leapfrog(
      proposal_q, momentum, proposal$gradient, step_size, metric, target
    )
    if (!moved$valid) {
      valid <- FALSE
      break
    }
    proposal_q <- moved$q
    momentum <- moved$momentum
    proposal <- moved$evaluated
  }
  energy_error <- if (valid) {
    -(proposal$logp - .nm_hmc_kinetic(momentum, metric)) +
      (evaluated$logp - .nm_hmc_kinetic(initial_momentum, metric))
  } else Inf
  acceptance <- if (is.finite(energy_error)) min(1, exp(-energy_error)) else 0
  accepted <- isTRUE(valid) && stats::runif(1) < acceptance
  list(
    q = if (accepted) proposal_q else q,
    evaluated = if (accepted) proposal else evaluated,
    acceptance = acceptance, accepted = accepted,
    divergence = !is.finite(energy_error) || abs(energy_error) > divergence_threshold,
    energy_error = energy_error, energy = energy,
    tree_depth = NA_integer_, leapfrog = n_leapfrog
  )
}

.nm_nuts_stop <- function(q_minus, q_plus, r_minus, r_plus, metric) {
  delta <- q_plus - q_minus
  sum(delta * .nm_hmc_velocity(r_minus, metric)) >= 0 &&
    sum(delta * .nm_hmc_velocity(r_plus, metric)) >= 0
}

.nm_nuts_tree <- function(q, r, evaluated, log_slice, direction, depth,
                          step_size, metric, target, joint0,
                          divergence_threshold) {
  # Recursive R calls otherwise accumulate self-referential lazy promises
  # (notably `q = q` and `depth = depth - 1`) before a deep tree is forced.
  force(q); force(r); force(evaluated); force(log_slice); force(direction)
  force(depth); force(step_size); force(metric); force(target); force(joint0)
  force(divergence_threshold)
  if (depth == 0L) {
    moved <- .nm_hmc_leapfrog(
      q, r, evaluated$gradient, direction * step_size, metric, target
    )
    if (!moved$valid) {
      return(list(q_minus = moved$q, r_minus = moved$momentum, e_minus = moved$evaluated,
                  q_plus = moved$q, r_plus = moved$momentum, e_plus = moved$evaluated,
                  q_proposal = q, e_proposal = evaluated, n = 0L, s = FALSE,
                  alpha = 0, n_alpha = 1L, divergent = TRUE, leapfrog = 1L))
    }
    joint <- moved$evaluated$logp - .nm_hmc_kinetic(moved$momentum, metric)
    error <- joint0 - joint
    return(list(
      q_minus = moved$q, r_minus = moved$momentum, e_minus = moved$evaluated,
      q_plus = moved$q, r_plus = moved$momentum, e_plus = moved$evaluated,
      q_proposal = moved$q, e_proposal = moved$evaluated,
      n = as.integer(log_slice <= joint),
      s = is.finite(joint) && log_slice - divergence_threshold < joint,
      alpha = min(1, exp(min(0, joint - joint0))), n_alpha = 1L,
      divergent = !is.finite(error) || abs(error) > divergence_threshold,
      leapfrog = 1L
    ))
  }
  left <- .nm_nuts_tree(
    q, r, evaluated, log_slice, direction, depth - 1L, step_size, metric,
    target, joint0, divergence_threshold
  )
  if (!left$s) return(left)
  right <- if (direction < 0) {
    .nm_nuts_tree(
      left$q_minus, left$r_minus, left$e_minus, log_slice, direction, depth - 1L,
      step_size, metric, target, joint0, divergence_threshold
    )
  } else {
    .nm_nuts_tree(
      left$q_plus, left$r_plus, left$e_plus, log_slice, direction, depth - 1L,
      step_size, metric, target, joint0, divergence_threshold
    )
  }
  proposal_q <- left$q_proposal
  proposal_e <- left$e_proposal
  if (right$n > 0L && stats::runif(1) < right$n / max(left$n + right$n, 1L)) {
    proposal_q <- right$q_proposal
    proposal_e <- right$e_proposal
  }
  q_minus <- if (direction < 0) right$q_minus else left$q_minus
  r_minus <- if (direction < 0) right$r_minus else left$r_minus
  e_minus <- if (direction < 0) right$e_minus else left$e_minus
  q_plus <- if (direction < 0) left$q_plus else right$q_plus
  r_plus <- if (direction < 0) left$r_plus else right$r_plus
  e_plus <- if (direction < 0) left$e_plus else right$e_plus
  list(
    q_minus = q_minus, r_minus = r_minus, e_minus = e_minus,
    q_plus = q_plus, r_plus = r_plus, e_plus = e_plus,
    q_proposal = proposal_q, e_proposal = proposal_e,
    n = left$n + right$n,
    s = right$s && .nm_nuts_stop(q_minus, q_plus, r_minus, r_plus, metric),
    alpha = left$alpha + right$alpha,
    n_alpha = left$n_alpha + right$n_alpha,
    divergent = left$divergent || right$divergent,
    leapfrog = left$leapfrog + right$leapfrog
  )
}

.nm_nuts_transition <- function(q, evaluated, step_size, metric, max_depth,
                                target, divergence_threshold) {
  momentum <- .nm_hmc_metric_sample(metric)
  joint0 <- evaluated$logp - .nm_hmc_kinetic(momentum, metric)
  log_slice <- joint0 - stats::rexp(1)
  q_minus <- q_plus <- proposal_q <- q
  r_minus <- r_plus <- momentum
  e_minus <- e_plus <- proposal_e <- evaluated
  n <- 1L
  active <- TRUE
  alpha <- 0
  n_alpha <- 0L
  divergent <- FALSE
  leapfrog <- 0L
  depth_reached <- 0L
  for (depth in 0:(max_depth - 1L)) {
    if (!active) break
    direction <- sample(c(-1L, 1L), 1L)
    tree <- if (direction < 0) {
      .nm_nuts_tree(q_minus, r_minus, e_minus, log_slice, direction, depth,
                    step_size, metric, target, joint0, divergence_threshold)
    } else {
      .nm_nuts_tree(q_plus, r_plus, e_plus, log_slice, direction, depth,
                    step_size, metric, target, joint0, divergence_threshold)
    }
    if (tree$s && tree$n > 0L && stats::runif(1) < tree$n / max(n + tree$n, 1L)) {
      proposal_q <- tree$q_proposal
      proposal_e <- tree$e_proposal
    }
    if (direction < 0) {
      q_minus <- tree$q_minus; r_minus <- tree$r_minus; e_minus <- tree$e_minus
    } else {
      q_plus <- tree$q_plus; r_plus <- tree$r_plus; e_plus <- tree$e_plus
    }
    n <- n + tree$n
    active <- tree$s && .nm_nuts_stop(q_minus, q_plus, r_minus, r_plus, metric)
    alpha <- alpha + tree$alpha
    n_alpha <- n_alpha + tree$n_alpha
    divergent <- divergent || tree$divergent
    leapfrog <- leapfrog + tree$leapfrog
    depth_reached <- depth + 1L
  }
  list(
    q = proposal_q, evaluated = proposal_e,
    acceptance = alpha / max(n_alpha, 1L), accepted = !identical(proposal_q, q),
    divergence = divergent, energy_error = NA_real_, energy = -joint0,
    tree_depth = depth_reached, leapfrog = leapfrog,
    max_depth_reached = depth_reached >= max_depth && active
  )
}

.nm_log_add_exp <- function(left, right) {
  if (!is.finite(left)) return(right)
  if (!is.finite(right)) return(left)
  maximum <- max(left, right)
  maximum + log(exp(left - maximum) + exp(right - maximum))
}

.nm_nuts_generalized_stop <- function(r_minus, r_plus, rho, metric) {
  velocity_minus <- .nm_hmc_velocity(r_minus, metric)
  velocity_plus <- .nm_hmc_velocity(r_plus, metric)
  sum(velocity_minus * rho) >= 0 && sum(velocity_plus * rho) >= 0
}

.nm_nuts_multinomial_tree <- function(q, r, evaluated, direction, depth,
                                      step_size, metric, target, joint0,
                                      divergence_threshold) {
  if (depth == 0L) {
    moved <- .nm_hmc_leapfrog(
      q, r, evaluated$gradient, direction * step_size, metric, target
    )
    joint <- if (moved$valid) {
      moved$evaluated$logp - .nm_hmc_kinetic(moved$momentum, metric)
    } else -Inf
    energy_error <- joint0 - joint
    valid <- moved$valid && is.finite(energy_error) &&
      abs(energy_error) <= divergence_threshold
    return(list(
      q_minus = moved$q, r_minus = moved$momentum, e_minus = moved$evaluated,
      q_plus = moved$q, r_plus = moved$momentum, e_plus = moved$evaluated,
      q_proposal = if (valid) moved$q else q,
      e_proposal = if (valid) moved$evaluated else evaluated,
      rho = if (moved$valid) moved$momentum else rep(0, length(r)),
      log_weight = if (valid) joint - joint0 else -Inf,
      active = valid,
      alpha = if (moved$valid) min(1, exp(min(0, joint - joint0))) else 0,
      n_alpha = 1L, divergent = !valid, leapfrog = 1L
    ))
  }
  left <- .nm_nuts_multinomial_tree(
    q, r, evaluated, direction, depth - 1L, step_size, metric, target,
    joint0, divergence_threshold
  )
  if (!left$active) return(left)
  right <- if (direction < 0) {
    .nm_nuts_multinomial_tree(
      left$q_minus, left$r_minus, left$e_minus, direction, depth - 1L,
      step_size, metric, target, joint0, divergence_threshold
    )
  } else {
    .nm_nuts_multinomial_tree(
      left$q_plus, left$r_plus, left$e_plus, direction, depth - 1L,
      step_size, metric, target, joint0, divergence_threshold
    )
  }
  # A subtree that already contains a divergence or U-turn is not a valid
  # reversible extension of its parent trajectory. Keep its integration
  # statistics for adaptation and diagnostics, but exclude its states from
  # the parent's multinomial candidate set.
  if (!right$active) {
    left$active <- FALSE
    left$alpha <- left$alpha + right$alpha
    left$n_alpha <- left$n_alpha + right$n_alpha
    left$divergent <- left$divergent || right$divergent
    left$leapfrog <- left$leapfrog + right$leapfrog
    return(left)
  }
  total_weight <- .nm_log_add_exp(left$log_weight, right$log_weight)
  choose_right <- is.finite(right$log_weight) &&
    stats::runif(1) < exp(right$log_weight - total_weight)
  q_minus <- if (direction < 0) right$q_minus else left$q_minus
  r_minus <- if (direction < 0) right$r_minus else left$r_minus
  e_minus <- if (direction < 0) right$e_minus else left$e_minus
  q_plus <- if (direction < 0) left$q_plus else right$q_plus
  r_plus <- if (direction < 0) left$r_plus else right$r_plus
  e_plus <- if (direction < 0) left$e_plus else right$e_plus
  rho <- left$rho + right$rho
  generalized_turn <- .nm_nuts_generalized_stop(
    r_minus, r_plus, rho, metric
  )
  list(
    q_minus = q_minus, r_minus = r_minus, e_minus = e_minus,
    q_plus = q_plus, r_plus = r_plus, e_plus = e_plus,
    q_proposal = if (choose_right) right$q_proposal else left$q_proposal,
    e_proposal = if (choose_right) right$e_proposal else left$e_proposal,
    rho = rho, log_weight = total_weight,
    active = right$active && generalized_turn,
    alpha = left$alpha + right$alpha,
    n_alpha = left$n_alpha + right$n_alpha,
    divergent = left$divergent || right$divergent,
    leapfrog = left$leapfrog + right$leapfrog
  )
}

.nm_nuts_multinomial_transition <- function(q, evaluated, step_size, metric,
                                            max_depth, target,
                                            divergence_threshold) {
  momentum <- .nm_hmc_metric_sample(metric)
  joint0 <- evaluated$logp - .nm_hmc_kinetic(momentum, metric)
  q_minus <- q_plus <- proposal_q <- q
  r_minus <- r_plus <- momentum
  e_minus <- e_plus <- proposal_e <- evaluated
  rho <- momentum
  log_weight <- 0
  active <- TRUE
  alpha <- 0
  n_alpha <- leapfrog <- depth_reached <- 0L
  divergent <- FALSE
  for (depth in 0:(max_depth - 1L)) {
    if (!active) break
    direction <- sample(c(-1L, 1L), 1L)
    tree <- if (direction < 0) {
      .nm_nuts_multinomial_tree(
        q_minus, r_minus, e_minus, direction, depth, step_size, metric,
        target, joint0, divergence_threshold
      )
    } else {
      .nm_nuts_multinomial_tree(
        q_plus, r_plus, e_plus, direction, depth, step_size, metric,
        target, joint0, divergence_threshold
      )
    }
    alpha <- alpha + tree$alpha
    n_alpha <- n_alpha + tree$n_alpha
    divergent <- divergent || tree$divergent
    leapfrog <- leapfrog + tree$leapfrog
    depth_reached <- depth + 1L
    if (!tree$active) {
      active <- FALSE
      break
    }
    total_weight <- .nm_log_add_exp(log_weight, tree$log_weight)
    if (is.finite(tree$log_weight) &&
        stats::runif(1) < exp(tree$log_weight - total_weight)) {
      proposal_q <- tree$q_proposal
      proposal_e <- tree$e_proposal
    }
    if (direction < 0) {
      q_minus <- tree$q_minus; r_minus <- tree$r_minus; e_minus <- tree$e_minus
    } else {
      q_plus <- tree$q_plus; r_plus <- tree$r_plus; e_plus <- tree$e_plus
    }
    rho <- rho + tree$rho
    log_weight <- total_weight
    active <- tree$active && .nm_nuts_generalized_stop(
      r_minus, r_plus, rho, metric
    )
  }
  list(
    q = proposal_q, evaluated = proposal_e,
    acceptance = alpha / max(n_alpha, 1L), accepted = !identical(proposal_q, q),
    divergence = divergent, energy_error = NA_real_, energy = -joint0,
    tree_depth = depth_reached, leapfrog = leapfrog,
    max_depth_reached = depth_reached >= max_depth && active
  )
}

.nm_hmc_warmup_schedule <- function(n_warmup, initial_buffer = 75L,
                                    terminal_buffer = 50L,
                                    first_window = 25L) {
  n_warmup <- as.integer(n_warmup)
  initial_buffer <- as.integer(initial_buffer)
  terminal_buffer <- as.integer(terminal_buffer)
  first_window <- as.integer(first_window)
  if (anyNA(c(n_warmup, initial_buffer, terminal_buffer, first_window)) ||
      n_warmup < 0L || initial_buffer < 0L || terminal_buffer < 0L ||
      first_window < 1L) {
    .nm_stop("HMC warmup buffers must be non-negative and the first window positive.")
  }
  if (n_warmup < 20L) {
    return(list(
      n_warmup = n_warmup, initial_buffer = n_warmup,
      terminal_buffer = 0L, windows = integer(), first_window = first_window
    ))
  }
  if (initial_buffer + terminal_buffer + first_window > n_warmup) {
    initial_buffer <- max(1L, floor(0.15 * n_warmup))
    terminal_buffer <- max(1L, floor(0.10 * n_warmup))
    first_window <- max(1L, n_warmup - initial_buffer - terminal_buffer)
  }
  slow <- n_warmup - initial_buffer - terminal_buffer
  if (slow <= 0L) {
    return(list(
      n_warmup = n_warmup, initial_buffer = n_warmup,
      terminal_buffer = 0L, windows = integer(), first_window = first_window
    ))
  }
  ends <- integer()
  consumed <- 0L
  width <- min(first_window, slow)
  while (consumed < slow) {
    remaining <- slow - consumed
    if (remaining < 2L * width) width <- remaining
    consumed <- consumed + width
    ends <- c(ends, initial_buffer + consumed)
    width <- min(2L * width, slow - consumed)
    if (width <= 0L) break
  }
  list(
    n_warmup = n_warmup, initial_buffer = initial_buffer,
    terminal_buffer = terminal_buffer, windows = ends,
    first_window = first_window
  )
}

.nm_hmc_mass_windows <- function(n_warmup) {
  schedule <- .nm_hmc_warmup_schedule(n_warmup)
  result <- schedule$windows
  attr(result, "initial_buffer") <- schedule$initial_buffer
  attr(result, "terminal_buffer") <- schedule$terminal_buffer
  result
}

.nm_hmc_adapted_metric <- function(samples, kind = "diagonal", blocks = list()) {
  samples <- as.matrix(samples)
  if (nrow(samples) < 2L || kind == "unit") {
    return(.nm_hmc_metric(ncol(samples), kind, blocks))
  }
  shrink <- nrow(samples) / (nrow(samples) + 5)
  covariance <- diag(0, ncol(samples))
  metric_blocks <- switch(kind,
    diagonal = lapply(seq_len(ncol(samples)), function(index) index),
    dense = list(seq_len(ncol(samples))),
    block = blocks,
    .nm_stop("Unknown HMC metric `", kind, "`.")
  )
  for (index in metric_blocks) {
    block <- if (length(index) == 1L) {
      matrix(stats::var(samples[, index]), 1L, 1L)
    } else stats::cov(samples[, index, drop = FALSE])
    block[!is.finite(block)] <- 0
    block <- shrink * block + (1 - shrink) * diag(1e-3, length(index))
    decomposition <- eigen((block + t(block)) / 2, symmetric = TRUE)
    values <- pmin(pmax(decomposition$values, 1e-3), 1e3)
    block <- decomposition$vectors %*% (values * t(decomposition$vectors))
    covariance[index, index] <- block
  }
  # `mass` is the momentum covariance M; the inverse posterior covariance is
  # therefore the Euclidean metric used by the Hamiltonian dynamics.
  .nm_hmc_metric_from_mass(solve(covariance), kind, blocks)
}

.nm_hmc_initialize <- function(target, strategy = c("auto", "model", "jitter"),
                               radius = 2, attempts = 100L) {
  strategy <- match.arg(strategy)
  radius <- as.numeric(radius)
  attempts <- as.integer(attempts)
  if (!is.finite(radius) || radius < 0 || is.na(attempts) || attempts < 1L) {
    .nm_stop("HMC initialization requires a non-negative radius and positive attempts.")
  }
  evaluate <- function(q, attempt, source) {
    value <- target$evaluate(q)
    if (is.finite(value$logp) && all(is.finite(value$gradient))) {
      list(q = q, evaluated = value, attempts = attempt, source = source)
    } else NULL
  }
  if (strategy %in% c("auto", "model")) {
    exact <- evaluate(target$initial, 1L, "model initial")
    if (!is.null(exact)) return(exact)
    if (strategy == "model") .nm_stop("The model initial point has non-finite posterior density or gradient.")
  }
  for (attempt in seq_len(attempts)) {
    scale <- radius * 0.5^floor((attempt - 1L) / 25L)
    candidate <- target$initial + stats::runif(length(target$initial), -scale, scale)
    value <- evaluate(candidate, attempt + as.integer(strategy == "auto"), "bounded random jitter")
    if (!is.null(value)) return(value)
  }
  .nm_stop("Unable to initialize HMC at a finite posterior density after ", attempts, " bounded attempts.")
}

.nm_hmc_ebfmi <- function(energy) {
  energy <- as.numeric(energy)
  energy <- energy[is.finite(energy)]
  if (length(energy) < 3L) return(NA_real_)
  denominator <- sum((energy - mean(energy))^2)
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  sum(diff(energy)^2) / denominator
}

.nm_mcmc_split_chains <- function(chains) {
  if (!length(chains)) return(list())
  n <- min(vapply(chains, nrow, integer(1)))
  if (n < 4L) return(list())
  n_half <- floor(n / 2)
  unlist(lapply(chains, function(x) list(
    x[seq_len(n_half), , drop = FALSE],
    x[n - n_half + seq_len(n_half), , drop = FALSE]
  )), recursive = FALSE)
}

.nm_mcmc_rank_normalize <- function(chains) {
  if (!length(chains)) return(chains)
  sizes <- vapply(chains, nrow, integer(1))
  combined <- do.call(rbind, chains)
  normalized <- apply(combined, 2L, function(value) {
    rank <- rank(value, ties.method = "average")
    stats::qnorm((rank - 3 / 8) / (length(rank) + 1 / 4))
  })
  if (is.null(dim(normalized))) normalized <- matrix(normalized, ncol = 1L)
  ends <- cumsum(sizes)
  starts <- c(1L, utils::head(ends, -1L) + 1L)
  Map(function(first, last) normalized[first:last, , drop = FALSE], starts, ends)
}

.nm_mcmc_basic_rhat <- function(chains) {
  if (length(chains) < 2L) return(rep(NA_real_, ncol(chains[[1L]])))
  n_half <- min(vapply(chains, nrow, integer(1)))
  means <- do.call(rbind, lapply(chains, colMeans))
  variances <- do.call(rbind, lapply(chains, function(x) apply(x, 2, stats::var)))
  between <- n_half * apply(means, 2, stats::var)
  within <- colMeans(variances)
  result <- sqrt(((n_half - 1) / n_half * within + between / n_half) / within)
  result[!is.finite(result)] <- NA_real_
  result
}

.nm_mcmc_rhat <- function(chains) {
  if (length(chains) < 2L) return(rep(NA_real_, ncol(chains[[1L]])))
  split <- .nm_mcmc_split_chains(chains)
  if (!length(split)) return(rep(NA_real_, ncol(chains[[1L]])))
  rank_split <- .nm_mcmc_rank_normalize(split)
  rank_rhat <- .nm_mcmc_basic_rhat(rank_split)
  combined <- do.call(rbind, split)
  medians <- apply(combined, 2L, stats::median)
  folded <- lapply(split, function(value) {
    abs(sweep(value, 2L, medians, FUN = "-"))
  })
  folded_rhat <- .nm_mcmc_basic_rhat(.nm_mcmc_rank_normalize(folded))
  result <- pmax(rank_rhat, folded_rhat, na.rm = TRUE)
  result[!is.finite(result)] <- NA_real_
  result
}

.nm_mcmc_ess_transformed <- function(chains) {
  if (!length(chains)) return(numeric())
  n <- min(vapply(chains, nrow, integer(1)))
  m <- length(chains)
  if (n < 4L) return(rep(NA_real_, ncol(chains[[1L]])))
  vapply(seq_len(ncol(chains[[1L]])), function(column) {
    values <- lapply(chains, function(value) value[seq_len(n), column])
    within <- mean(vapply(values, stats::var, numeric(1)))
    if (!is.finite(within) || within <= 0) return(NA_real_)
    means <- vapply(values, mean, numeric(1))
    between <- if (m > 1L) n * stats::var(means) else 0
    variance <- (n - 1) / n * within + between / n
    autocovariance <- lapply(values, function(value) {
      as.numeric(stats::acf(
        value, lag.max = n - 1L, type = "covariance",
        plot = FALSE, demean = TRUE
      )$acf)
    })
    rho <- vapply(seq_len(n), function(index) {
      1 - (within - mean(vapply(
        autocovariance, `[`, numeric(1), index
      ))) / variance
    }, numeric(1))
    rho[[1L]] <- 1
    last <- if (length(rho) %% 2L) length(rho) - 1L else length(rho)
    pairs <- rho[seq.int(1L, last, by = 2L)] +
      rho[seq.int(2L, last, by = 2L)]
    negative <- which(!is.finite(pairs) | pairs < 0)
    if (length(negative)) pairs <- utils::head(pairs, negative[[1L]] - 1L)
    if (length(pairs) > 1L) {
      for (index in 2:length(pairs)) {
        pairs[[index]] <- min(pairs[[index]], pairs[[index - 1L]])
      }
    }
    tau <- max(-1 + 2 * sum(pairs), 1 / (m * n))
    min(m * n, m * n / tau)
  }, numeric(1))
}

.nm_mcmc_ess <- function(chains) {
  split <- .nm_mcmc_split_chains(chains)
  if (!length(split)) return(rep(NA_real_, ncol(chains[[1L]])))
  .nm_mcmc_ess_transformed(.nm_mcmc_rank_normalize(split))
}

.nm_mcmc_diagnostics <- function(chains) {
  split <- .nm_mcmc_split_chains(chains)
  columns <- ncol(chains[[1L]])
  if (!length(split)) {
    empty <- rep(NA_real_, columns)
    return(list(rhat = empty, bulk_ess = empty, tail_ess = empty,
                mcse_mean = empty))
  }
  combined <- do.call(rbind, split)
  bulk <- .nm_mcmc_ess_transformed(.nm_mcmc_rank_normalize(split))
  tail <- vapply(seq_len(columns), function(column) {
    quantiles <- stats::quantile(
      combined[, column], c(0.05, 0.95), names = FALSE
    )
    lower <- lapply(split, function(value) {
      matrix(as.numeric(value[, column] <= quantiles[[1L]]), ncol = 1L)
    })
    upper <- lapply(split, function(value) {
      matrix(as.numeric(value[, column] >= quantiles[[2L]]), ncol = 1L)
    })
    values <- c(
      .nm_mcmc_ess_transformed(lower),
      .nm_mcmc_ess_transformed(upper)
    )
    if (all(!is.finite(values))) NA_real_ else min(values, na.rm = TRUE)
  }, numeric(1))
  pooled_sd <- apply(combined, 2L, stats::sd)
  list(
    rhat = .nm_mcmc_rhat(chains), bulk_ess = bulk, tail_ess = tail,
    mcse_mean = pooled_sd / sqrt(bulk)
  )
}

.nm_mcmc_native_row <- function(evaluated, context) {
  c(
    evaluated$parameters$theta, evaluated$parameters$sigma,
    evaluated$parameters$omega, as.vector(t(evaluated$eta)),
    LOG_POSTERIOR = evaluated$logp
  )
}

.nm_est_hmc <- function(context, map, method = c("HMC", "NUTS"),
                        n_warmup = 500L, n_sample = 1000L, n_thin = 1L,
                        n_chains = 4L, seed = 20260719L,
                        step_size = NULL, target_acceptance = 0.8,
                        adapt_mass = TRUE, n_leapfrog = 10L,
                        max_depth = 10L, divergence_threshold = 1000,
                        print_every = 0L,
                        sampler_backend = c("native", "r"),
                        geometry = c("auto", "centered", "mu_noncentered"),
                        nuts_variant = c("auto", "classic_slice", "multinomial"),
                        hmc_metric = c("auto", "unit", "diagonal", "dense", "block"),
                        hmc_metric_blocks = NULL,
                        adapt_initial_buffer = 75L,
                        adapt_terminal_buffer = 50L,
                        adapt_first_window = 25L,
                        initialization = c("auto", "model", "jitter"),
                        initialization_radius = 2,
                        initialization_attempts = 100L,
                        ...) {
  method <- match.arg(method)
  sampler_backend <- match.arg(sampler_backend)
  geometry <- match.arg(geometry)
  nuts_variant <- match.arg(nuts_variant)
  hmc_metric <- match.arg(hmc_metric)
  initialization <- match.arg(initialization)
  resolved_metric <- if (hmc_metric == "auto") "diagonal" else hmc_metric
  resolved_nuts_variant <- if (method != "NUTS") {
    "not_applicable"
  } else if (nuts_variant == "auto") {
    if (.nm_liber_optimized(context)) "multinomial" else "classic_slice"
  } else nuts_variant
  n_warmup <- as.integer(n_warmup); n_sample <- as.integer(n_sample)
  n_thin <- as.integer(n_thin); n_chains <- as.integer(n_chains)
  if (anyNA(c(n_warmup, n_sample, n_thin, n_chains)) ||
      n_warmup < 0L || n_sample < 1L || n_thin < 1L || n_chains < 1L) {
    .nm_stop(method, " requires non-negative warmup and positive samples, thinning, and chains.")
  }
  n_leapfrog <- as.integer(n_leapfrog); max_depth <- as.integer(max_depth)
  if (is.na(n_leapfrog) || n_leapfrog < 1L || is.na(max_depth) || max_depth < 1L) {
    .nm_stop("`n_leapfrog` and `max_depth` must be positive integers.")
  }
  if (!is.null(step_size) &&
      (length(step_size) != 1L || !is.finite(step_size) || step_size <= 0)) {
    .nm_stop("`step_size` must be NULL or a positive finite number.")
  }
  if (length(divergence_threshold) != 1L || !is.finite(divergence_threshold) ||
      divergence_threshold <= 0) {
    .nm_stop("`divergence_threshold` must be a positive finite number.")
  }
  if (!is.finite(target_acceptance) || target_acceptance <= 0 || target_acceptance >= 1) {
    .nm_stop("`target_acceptance` must lie strictly between zero and one.")
  }
  initialization_attempts <- as.integer(initialization_attempts)
  initialization_radius <- as.numeric(initialization_radius)
  if (is.na(initialization_attempts) || initialization_attempts < 1L ||
      length(initialization_radius) != 1L || !is.finite(initialization_radius) ||
      initialization_radius < 0) {
    .nm_stop("HMC initialization requires a non-negative radius and positive attempts.")
  }
  maximum_seed <- .Machine$integer.max - n_chains + 1
  if (length(seed) != 1L || !is.finite(seed) || seed < 0 ||
      seed != floor(seed) || seed > maximum_seed) {
    .nm_stop(
      "`seed` must be a non-negative integer no larger than ",
      maximum_seed, " for the requested number of chains."
    )
  }
  geometry_plan <- .nm_hmc_geometry_plan(context, map, geometry)
  target <- .nm_hmc_target(context, map, geometry_plan)
  if (!length(target$initial)) .nm_stop(method, " requires at least one unknown parameter or ETA.")
  metric_blocks <- .nm_hmc_metric_blocks(
    length(target$initial), context,
    if (resolved_metric == "block") hmc_metric_blocks else NULL
  )
  warmup_schedule <- .nm_hmc_warmup_schedule(
    n_warmup, adapt_initial_buffer, adapt_terminal_buffer,
    adapt_first_window
  )
  objective_evaluations <- NA_integer_
  if (sampler_backend == "native") {
    native_config <- .nm_hmc_cpp_config(context, map, target)
    native_config$nuts_variant <- resolved_nuts_variant
    native_config$hmc_metric <- resolved_metric
    native_config$hmc_metric_blocks <- lapply(metric_blocks, function(index) {
      as.integer(index - 1L)
    })
    native_config$adapt_initial_buffer <- warmup_schedule$initial_buffer
    native_config$adapt_terminal_buffer <- warmup_schedule$terminal_buffer
    native_config$adapt_first_window <- warmup_schedule$first_window
    native_config$initialization <- initialization
    native_config$initialization_radius <- as.numeric(initialization_radius)
    native_config$initialization_attempts <- as.integer(initialization_attempts)
    native <- .liberation_hmc_sample(
      target$full_tape$pointer, native_config,
      method, n_warmup, n_sample, n_thin, n_chains, as.numeric(seed),
      if (is.null(step_size)) NaN else as.numeric(step_size),
      target_acceptance, isTRUE(adapt_mass), n_leapfrog, max_depth,
      divergence_threshold, as.integer(print_every)
    )
    chain_results <- native$chains
    chain_diagnostics <- lapply(native$diagnostics, function(value) {
      value$trace <- as.data.frame(value$trace, stringsAsFactors = FALSE)
      value$trace <- data.frame(
        iteration = seq_len(nrow(value$trace)), value$trace,
        check.names = FALSE
      )
      value$trace$divergence <- as.logical(value$trace$divergence)
      value
    })
    objective_evaluations <- as.integer(native$objective_evaluations)
  } else {
    chain_results <- vector("list", n_chains)
    chain_diagnostics <- vector("list", n_chains)
    for (chain_id in seq_len(n_chains)) {
      set.seed(as.integer(seed) + chain_id - 1L)
      initialized <- .nm_hmc_initialize(
        target, initialization, initialization_radius,
        initialization_attempts
      )
      q <- initialized$q
      evaluated <- initialized$evaluated
      metric <- .nm_hmc_metric(length(q), resolved_metric, metric_blocks)
      epsilon <- if (is.null(step_size)) {
        .nm_hmc_find_step(q, evaluated, metric, target)
      } else as.numeric(step_size)
      dual <- .nm_dual_average(epsilon, target_acceptance)
      warmup_q <- matrix(NA_real_, max(n_warmup, 1L), length(q))
      metric_window_ends <- warmup_schedule$windows
      metric_window_start <- if (length(metric_window_ends)) {
        warmup_schedule$initial_buffer + 1L
      } else 1L
      total <- n_warmup + n_sample * n_thin
      draws <- matrix(NA_real_, n_sample,
                      nrow(context$model$THETAS) + nrow(context$model$SIGMAS) +
                        nrow(context$model$OMEGAS) + context$n_subjects * context$n_eta + 1L)
      trace <- data.frame(
        iteration = seq_len(total), acceptance = NA_real_, divergence = FALSE,
        tree_depth = NA_integer_, leapfrog = NA_integer_, step_size = NA_real_,
        energy = NA_real_,
        stringsAsFactors = FALSE
      )
      keep <- 0L
      for (iteration in seq_len(total)) {
        transition <- if (method == "NUTS") {
          if (resolved_nuts_variant == "multinomial") {
            .nm_nuts_multinomial_transition(
              q, evaluated, epsilon, metric, max_depth, target,
              divergence_threshold
            )
          } else {
            .nm_nuts_transition(q, evaluated, epsilon, metric, max_depth, target,
                                divergence_threshold)
          }
        } else {
          .nm_hmc_transition(q, evaluated, epsilon, metric, n_leapfrog, target,
                             divergence_threshold)
        }
        q <- transition$q
        evaluated <- transition$evaluated
        trace$acceptance[[iteration]] <- transition$acceptance
        trace$divergence[[iteration]] <- transition$divergence
        trace$tree_depth[[iteration]] <- transition$tree_depth
        trace$leapfrog[[iteration]] <- transition$leapfrog
        trace$step_size[[iteration]] <- epsilon
        trace$energy[[iteration]] <- transition$energy
        if (iteration <= n_warmup) {
          warmup_q[iteration, ] <- q
          epsilon <- dual$update(transition$acceptance)
          if (isTRUE(adapt_mass) && resolved_metric != "unit" &&
              iteration %in% metric_window_ends) {
            metric <- .nm_hmc_adapted_metric(
              warmup_q[metric_window_start:iteration, , drop = FALSE],
              resolved_metric, metric_blocks
            )
            epsilon <- .nm_hmc_find_step(q, evaluated, metric, target)
            dual <- .nm_dual_average(epsilon, target_acceptance)
            metric_window_start <- iteration + 1L
          }
          if (iteration == n_warmup) epsilon <- dual$final()
        } else if ((iteration - n_warmup) %% n_thin == 0L) {
          keep <- keep + 1L
          draws[keep, ] <- .nm_mcmc_native_row(evaluated, context)
        }
        if (print_every > 0L && iteration %% print_every == 0L) {
          cat(sprintf(
            "[LibeRation] %s CHAIN %d ITERATION %d LOGPOST %.10g ACCEPT %.3f STEP %.5g DIVERGENT %s\n",
            method, chain_id, iteration, evaluated$logp, transition$acceptance,
            epsilon, transition$divergence
          ))
          try(flush(stdout()), silent = TRUE)
        }
      }
      chain_results[[chain_id]] <- draws
      post <- seq.int(n_warmup + 1L, total)
      chain_diagnostics[[chain_id]] <- list(
        trace = trace, step_size = epsilon,
        mass = if (resolved_metric == "diagonal") diag(metric$mass) else metric$mass,
        metric = resolved_metric, metric_blocks = metric_blocks,
        warmup_schedule = warmup_schedule,
        initialization = initialized[c("source", "attempts")],
        ebfmi = .nm_hmc_ebfmi(trace$energy[post]),
        divergences = sum(trace$divergence[post]),
        mean_acceptance = mean(trace$acceptance[post]),
        max_depth_hits = if (method == "NUTS") {
          sum(trace$tree_depth[post] >= max_depth, na.rm = TRUE)
        } else 0L
      )
    }
  }
  n_theta <- nrow(context$model$THETAS); n_sigma <- nrow(context$model$SIGMAS)
  n_omega <- nrow(context$model$OMEGAS)
  column_names <- c(
    .nm_numbered_names("THETA", n_theta), .nm_numbered_names("SIGMA", n_sigma),
    .nm_numbered_names("OMEGA", n_omega),
    if (context$n_eta) unlist(lapply(seq_len(context$n_subjects), function(subject) {
      paste0("ETA", subject, "_", seq_len(context$n_eta))
    })) else character(), "LOG_POSTERIOR"
  )
  chain_results <- lapply(chain_results, function(x) { colnames(x) <- column_names; x })
  chain <- do.call(rbind, chain_results)
  population_names <- .nm_parameter_names(
    context$model$THETAS$Value, context$model$SIGMAS$Value, context$model$OMEGAS$Value
  )
  population_chain <- chain[, population_names, drop = FALSE]
  parameters <- list(
    theta = colMeans(chain[, .nm_numbered_names("THETA", n_theta), drop = FALSE]),
    sigma = colMeans(chain[, .nm_numbered_names("SIGMA", n_sigma), drop = FALSE]),
    omega = colMeans(chain[, .nm_numbered_names("OMEGA", n_omega), drop = FALSE])
  )
  eta_start <- n_theta + n_sigma + n_omega
  eta <- if (context$n_eta) matrix(
    colMeans(chain[, eta_start + seq_len(context$n_subjects * context$n_eta), drop = FALSE]),
    context$n_subjects, context$n_eta, byrow = TRUE
  ) else matrix(numeric(), context$n_subjects, 0L)
  modes <- lapply(seq_len(context$n_subjects), function(subject) {
    list(par = eta[subject, ], convergence = 0L, jitter = 0)
  })
  optimizer <- list(
    convergence = 0L, message = paste(method, "sampling completed"),
    counts = if (is.finite(objective_evaluations)) {
      c(`function` = objective_evaluations, gradient = objective_evaluations)
    } else {
      c(
        `function` = sum(vapply(
          chain_diagnostics, function(x) nrow(x$trace), integer(1)
        )),
        gradient = NA_integer_
      )
    },
    iterations = n_warmup + n_sample * n_thin,
    objective_evaluations = objective_evaluations,
    backend = paste0(sampler_backend, "-cppad-", tolower(method))
  )
  ebfmi <- vapply(chain_diagnostics, function(value) {
    as.numeric(value$ebfmi %||% NA_real_)[[1L]]
  }, numeric(1))
  fit <- .nm_fit_result(
    context, method, parameters, -2 * max(chain[, "LOG_POSTERIOR"]), modes, optimizer,
    diagnostics = list(
      sampler = method, n_warmup = n_warmup, n_sample = n_sample,
      n_thin = n_thin, n_chains = n_chains, seed = seed,
      target_acceptance = target_acceptance,
      divergences = sum(vapply(chain_diagnostics, `[[`, numeric(1), "divergences")),
      mean_acceptance = mean(vapply(chain_diagnostics, `[[`, numeric(1), "mean_acceptance")),
      max_depth_hits = sum(vapply(chain_diagnostics, `[[`, numeric(1), "max_depth_hits")),
      ebfmi = ebfmi,
      low_ebfmi_chains = which(is.finite(ebfmi) & ebfmi < 0.3),
      chain = chain_diagnostics,
      gradient = "exact joint CppAD gradient",
      geometry = target$geometry,
      nuts = if (method == "NUTS") list(
        variant_requested = nuts_variant,
        variant_resolved = resolved_nuts_variant,
        trajectory_sampling = if (resolved_nuts_variant == "multinomial") {
          "multinomial progressive sampling"
        } else "slice sampling",
        metric_requested = hmc_metric,
        metric_resolved = resolved_metric,
        mass_adaptation = paste("three-stage windowed regularized", resolved_metric, "precision"),
        warmup_schedule = warmup_schedule,
        initialization = list(
          strategy = initialization, radius = initialization_radius,
          maximum_attempts = initialization_attempts
        ),
        ebfmi_threshold = 0.3
      ) else NULL,
      sampler_backend = sampler_backend,
      objective_evaluations = objective_evaluations,
      objective_semantics = list(
        type = "negative_twice_maximum_sampled_log_posterior_sampling_coordinates",
        likelihood_comparable = FALSE,
        reported_point = "highest retained sampled log-posterior",
        recommendation = "Use WAIC or PSIS-LOO for predictive model comparison."
      )
    )
  )
  population_covariance <- if (nrow(population_chain) > 1L) stats::cov(population_chain) else
    matrix(NA_real_, ncol(population_chain), ncol(population_chain),
           dimnames = list(population_names, population_names))
  population_sd <- apply(population_chain, 2, stats::sd)
  correlation <- population_covariance / outer(population_sd, population_sd)
  diag(correlation) <- 1
  sampling_diagnostics <- .nm_mcmc_diagnostics(lapply(
    chain_results, function(x) x[, population_names, drop = FALSE]
  ))
  sampling_diagnostics <- lapply(sampling_diagnostics, function(value) {
    stats::setNames(value, population_names)
  })
  fit$chain <- chain
  fit$chains <- chain_results
  fit$posterior <- list(
    mean = colMeans(chain), sd = apply(chain, 2, stats::sd),
    quantile = apply(chain, 2, stats::quantile, probs = c(0.025, 0.5, 0.975)),
    population = list(
      mean = colMeans(population_chain), sd = population_sd,
      quantile = apply(population_chain, 2, stats::quantile,
                       probs = c(0.025, 0.5, 0.975)),
      covariance = population_covariance, correlation = correlation,
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
