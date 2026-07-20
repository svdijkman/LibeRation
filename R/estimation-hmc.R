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

.nm_hmc_target <- function(context, map) {
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
    eta <- if (n_eta_total) matrix(
      q[n_outer + seq_len(n_eta_total)], context$n_subjects, context$n_eta,
      byrow = TRUE
    ) else matrix(numeric(), context$n_subjects, 0L)
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
    outer_gradient <- outer_gradient * transformed$derivative +
      transformed$log_jacobian_gradient
    eta_gradient <- if (n_eta_total) {
      -0.5 * as.numeric(evaluated$gradient[eta_positions])
    } else numeric()
    list(
      logp = -0.5 * evaluated$value + prior + map$log_jacobian(parameters) +
        transformed$log_jacobian,
      gradient = c(outer_gradient, eta_gradient), parameters = parameters, eta = eta,
      outer = transformed$value
    )
  }
  list(evaluate = evaluate, initial = initial, bounds = bounds, map = map,
       n_outer = n_outer, n_eta_total = n_eta_total, full_tape = full_tape)
}

.nm_hmc_leapfrog <- function(q, momentum, gradient, epsilon, mass, target) {
  next_momentum <- momentum + 0.5 * epsilon * gradient
  next_q <- q + epsilon * next_momentum / mass
  evaluated <- target$evaluate(next_q)
  if (!is.finite(evaluated$logp) || any(!is.finite(evaluated$gradient))) {
    return(list(q = next_q, momentum = next_momentum,
                evaluated = evaluated, valid = FALSE))
  }
  next_momentum <- next_momentum + 0.5 * epsilon * evaluated$gradient
  list(q = next_q, momentum = next_momentum, evaluated = evaluated, valid = TRUE)
}

.nm_hmc_find_step <- function(q, evaluated, mass, target) {
  epsilon <- 1
  momentum <- stats::rnorm(length(q), sd = sqrt(mass))
  proposal <- .nm_hmc_leapfrog(q, momentum, evaluated$gradient, epsilon, mass, target)
  log_accept <- if (proposal$valid) {
    proposal$evaluated$logp - sum(proposal$momentum^2 / mass) / 2 -
      evaluated$logp + sum(momentum^2 / mass) / 2
  } else -Inf
  direction <- if (is.finite(log_accept) && log_accept > log(0.5)) 1 else -1
  for (iteration in seq_len(20L)) {
    candidate <- epsilon * if (direction > 0) 2 else 0.5
    proposal <- .nm_hmc_leapfrog(q, momentum, evaluated$gradient, candidate, mass, target)
    candidate_accept <- if (proposal$valid) {
      proposal$evaluated$logp - sum(proposal$momentum^2 / mass) / 2 -
        evaluated$logp + sum(momentum^2 / mass) / 2
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

.nm_hmc_transition <- function(q, evaluated, step_size, mass, n_leapfrog,
                               target, divergence_threshold) {
  momentum <- stats::rnorm(length(q), sd = sqrt(mass))
  initial_momentum <- momentum
  proposal_q <- q
  proposal <- evaluated
  valid <- TRUE
  for (step in seq_len(n_leapfrog)) {
    moved <- .nm_hmc_leapfrog(
      proposal_q, momentum, proposal$gradient, step_size, mass, target
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
    -(proposal$logp - sum(momentum^2 / mass) / 2) +
      (evaluated$logp - sum(initial_momentum^2 / mass) / 2)
  } else Inf
  acceptance <- if (is.finite(energy_error)) min(1, exp(-energy_error)) else 0
  accepted <- isTRUE(valid) && stats::runif(1) < acceptance
  list(
    q = if (accepted) proposal_q else q,
    evaluated = if (accepted) proposal else evaluated,
    acceptance = acceptance, accepted = accepted,
    divergence = !is.finite(energy_error) || abs(energy_error) > divergence_threshold,
    energy_error = energy_error, tree_depth = NA_integer_, leapfrog = n_leapfrog
  )
}

.nm_nuts_stop <- function(q_minus, q_plus, r_minus, r_plus, mass) {
  delta <- q_plus - q_minus
  sum(delta * r_minus / mass) >= 0 && sum(delta * r_plus / mass) >= 0
}

.nm_nuts_tree <- function(q, r, evaluated, log_slice, direction, depth,
                          step_size, mass, target, joint0,
                          divergence_threshold) {
  # Recursive R calls otherwise accumulate self-referential lazy promises
  # (notably `q = q` and `depth = depth - 1`) before a deep tree is forced.
  force(q); force(r); force(evaluated); force(log_slice); force(direction)
  force(depth); force(step_size); force(mass); force(target); force(joint0)
  force(divergence_threshold)
  if (depth == 0L) {
    moved <- .nm_hmc_leapfrog(
      q, r, evaluated$gradient, direction * step_size, mass, target
    )
    if (!moved$valid) {
      return(list(q_minus = moved$q, r_minus = moved$momentum, e_minus = moved$evaluated,
                  q_plus = moved$q, r_plus = moved$momentum, e_plus = moved$evaluated,
                  q_proposal = q, e_proposal = evaluated, n = 0L, s = FALSE,
                  alpha = 0, n_alpha = 1L, divergent = TRUE, leapfrog = 1L))
    }
    joint <- moved$evaluated$logp - sum(moved$momentum^2 / mass) / 2
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
    q, r, evaluated, log_slice, direction, depth - 1L, step_size, mass,
    target, joint0, divergence_threshold
  )
  if (!left$s) return(left)
  right <- if (direction < 0) {
    .nm_nuts_tree(
      left$q_minus, left$r_minus, left$e_minus, log_slice, direction, depth - 1L,
      step_size, mass, target, joint0, divergence_threshold
    )
  } else {
    .nm_nuts_tree(
      left$q_plus, left$r_plus, left$e_plus, log_slice, direction, depth - 1L,
      step_size, mass, target, joint0, divergence_threshold
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
    s = right$s && .nm_nuts_stop(q_minus, q_plus, r_minus, r_plus, mass),
    alpha = left$alpha + right$alpha,
    n_alpha = left$n_alpha + right$n_alpha,
    divergent = left$divergent || right$divergent,
    leapfrog = left$leapfrog + right$leapfrog
  )
}

.nm_nuts_transition <- function(q, evaluated, step_size, mass, max_depth,
                                target, divergence_threshold) {
  momentum <- stats::rnorm(length(q), sd = sqrt(mass))
  joint0 <- evaluated$logp - sum(momentum^2 / mass) / 2
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
                    step_size, mass, target, joint0, divergence_threshold)
    } else {
      .nm_nuts_tree(q_plus, r_plus, e_plus, log_slice, direction, depth,
                    step_size, mass, target, joint0, divergence_threshold)
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
    active <- tree$s && .nm_nuts_stop(q_minus, q_plus, r_minus, r_plus, mass)
    alpha <- alpha + tree$alpha
    n_alpha <- n_alpha + tree$n_alpha
    divergent <- divergent || tree$divergent
    leapfrog <- leapfrog + tree$leapfrog
    depth_reached <- depth + 1L
  }
  list(
    q = proposal_q, evaluated = proposal_e,
    acceptance = alpha / max(n_alpha, 1L), accepted = !identical(proposal_q, q),
    divergence = divergent, energy_error = NA_real_,
    tree_depth = depth_reached, leapfrog = leapfrog,
    max_depth_reached = depth_reached >= max_depth && active
  )
}

.nm_mcmc_rhat <- function(chains) {
  if (length(chains) < 2L) return(rep(NA_real_, ncol(chains[[1L]])))
  n <- min(vapply(chains, nrow, integer(1)))
  if (n < 4L) return(rep(NA_real_, ncol(chains[[1L]])))
  n_half <- floor(n / 2)
  split <- unlist(lapply(chains, function(x) list(
    x[seq_len(n_half), , drop = FALSE],
    x[n - n_half + seq_len(n_half), , drop = FALSE]
  )), recursive = FALSE)
  means <- do.call(rbind, lapply(split, colMeans))
  variances <- do.call(rbind, lapply(split, function(x) apply(x, 2, stats::var)))
  between <- n_half * apply(means, 2, stats::var)
  within <- colMeans(variances)
  sqrt(((n_half - 1) / n_half * within + between / n_half) / within)
}

.nm_mcmc_ess <- function(chains) {
  combined <- do.call(rbind, chains)
  vapply(seq_len(ncol(combined)), function(column) {
    x <- combined[, column]
    n <- length(x)
    if (n < 4L || stats::var(x) <= 0) return(NA_real_)
    lag_max <- min(n - 1L, max(10L, floor(sqrt(n))))
    rho <- stats::acf(x, lag.max = lag_max, plot = FALSE, demean = TRUE)$acf[-1L]
    pair_sum <- numeric()
    if (length(rho) >= 2L) {
      pair_sum <- rho[seq(1L, length(rho) - 1L, by = 2L)] +
        rho[seq(2L, length(rho), by = 2L)]
      pair_sum <- pair_sum[seq_len(match(TRUE, pair_sum < 0, nomatch = length(pair_sum) + 1L) - 1L)]
    }
    min(n, n / max(1 + 2 * sum(pair_sum), 1 / n))
  }, numeric(1))
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
                        print_every = 0L, ...) {
  method <- match.arg(method)
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
  target <- .nm_hmc_target(context, map)
  if (!length(target$initial)) .nm_stop(method, " requires at least one unknown parameter or ETA.")
  chain_results <- vector("list", n_chains)
  chain_diagnostics <- vector("list", n_chains)
  final_evaluated <- NULL
  for (chain_id in seq_len(n_chains)) {
    set.seed(as.integer(seed) + chain_id - 1L)
    q <- target$initial + stats::rnorm(length(target$initial), sd = 0.02)
    evaluated <- target$evaluate(q)
    if (!is.finite(evaluated$logp)) {
      q <- target$initial
      evaluated <- target$evaluate(q)
    }
    if (!is.finite(evaluated$logp)) .nm_stop("Unable to initialize ", method, " at a finite posterior density.")
    mass <- rep(1, length(q))
    epsilon <- if (is.null(step_size)) .nm_hmc_find_step(q, evaluated, mass, target) else as.numeric(step_size)
    dual <- .nm_dual_average(epsilon, target_acceptance)
    warmup_q <- matrix(NA_real_, max(n_warmup, 1L), length(q))
    total <- n_warmup + n_sample * n_thin
    draws <- matrix(NA_real_, n_sample,
                    nrow(context$model$THETAS) + nrow(context$model$SIGMAS) +
                      nrow(context$model$OMEGAS) + context$n_subjects * context$n_eta + 1L)
    trace <- data.frame(
      iteration = seq_len(total), acceptance = NA_real_, divergence = FALSE,
      tree_depth = NA_integer_, leapfrog = NA_integer_, step_size = NA_real_,
      stringsAsFactors = FALSE
    )
    keep <- 0L
    for (iteration in seq_len(total)) {
      transition <- if (method == "NUTS") {
        .nm_nuts_transition(q, evaluated, epsilon, mass, max_depth, target,
                            divergence_threshold)
      } else {
        .nm_hmc_transition(q, evaluated, epsilon, mass, n_leapfrog, target,
                           divergence_threshold)
      }
      q <- transition$q
      evaluated <- transition$evaluated
      trace$acceptance[[iteration]] <- transition$acceptance
      trace$divergence[[iteration]] <- transition$divergence
      trace$tree_depth[[iteration]] <- transition$tree_depth
      trace$leapfrog[[iteration]] <- transition$leapfrog
      trace$step_size[[iteration]] <- epsilon
      if (iteration <= n_warmup) {
        warmup_q[iteration, ] <- q
        epsilon <- dual$update(transition$acceptance)
        if (isTRUE(adapt_mass) && n_warmup >= 20L && iteration == floor(n_warmup / 2)) {
          variance <- apply(warmup_q[seq_len(iteration), , drop = FALSE], 2, stats::var)
          mass <- pmin(pmax(variance + 1e-3, 1e-3), 1e3)
          epsilon <- .nm_hmc_find_step(q, evaluated, mass, target)
          dual <- .nm_dual_average(epsilon, target_acceptance)
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
    chain_diagnostics[[chain_id]] <- list(
      trace = trace, step_size = epsilon, mass = mass,
      divergences = sum(trace$divergence[(n_warmup + 1L):total]),
      mean_acceptance = mean(trace$acceptance[(n_warmup + 1L):total]),
      max_depth_hits = if (method == "NUTS") {
        sum(trace$tree_depth[(n_warmup + 1L):total] >= max_depth, na.rm = TRUE)
      } else 0L
    )
    final_evaluated <- evaluated
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
    counts = c(`function` = sum(vapply(chain_diagnostics, function(x) nrow(x$trace), integer(1))),
               gradient = NA_integer_),
    iterations = n_warmup + n_sample * n_thin,
    objective_evaluations = NA_integer_, backend = paste0("cppad-", tolower(method))
  )
  fit <- .nm_fit_result(
    context, method, parameters, -2 * max(chain[, "LOG_POSTERIOR"]), modes, optimizer,
    diagnostics = list(
      sampler = method, n_warmup = n_warmup, n_sample = n_sample,
      n_thin = n_thin, n_chains = n_chains, seed = seed,
      target_acceptance = target_acceptance,
      divergences = sum(vapply(chain_diagnostics, `[[`, numeric(1), "divergences")),
      mean_acceptance = mean(vapply(chain_diagnostics, `[[`, numeric(1), "mean_acceptance")),
      max_depth_hits = sum(vapply(chain_diagnostics, `[[`, numeric(1), "max_depth_hits")),
      chain = chain_diagnostics,
      gradient = "exact joint CppAD gradient"
    )
  )
  population_covariance <- if (nrow(population_chain) > 1L) stats::cov(population_chain) else
    matrix(NA_real_, ncol(population_chain), ncol(population_chain),
           dimnames = list(population_names, population_names))
  population_sd <- apply(population_chain, 2, stats::sd)
  correlation <- population_covariance / outer(population_sd, population_sd)
  diag(correlation) <- 1
  rhat <- .nm_mcmc_rhat(lapply(chain_results, function(x) x[, population_names, drop = FALSE]))
  ess <- .nm_mcmc_ess(lapply(chain_results, function(x) x[, population_names, drop = FALSE]))
  names(rhat) <- names(ess) <- population_names
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
      rhat = rhat, ess = ess
    )
  )
  fit
}
