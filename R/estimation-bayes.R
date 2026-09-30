# BAYES estimation.
# Split from estimation-stochastic.R as a behaviour-preserving source move.

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

