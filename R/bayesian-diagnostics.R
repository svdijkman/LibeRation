.nm_bayes_chain_indices <- function(fit, draws, seed) {
  if (!inherits(fit, "nm_fit") || !fit$method %in% c("BAYES", "HMC", "NUTS") ||
      is.null(fit$chain)) {
    .nm_stop("Bayesian predictive diagnostics require a BAYES, HMC, or NUTS fit with saved draws.")
  }
  available <- nrow(fit$chain)
  if (is.null(draws)) return(seq_len(available))
  if (length(draws) == 1L) {
    draws <- as.integer(draws)
    if (is.na(draws) || draws < 1L) .nm_stop("`draws` must be positive.")
    if (draws >= available) return(seq_len(available))
    set.seed(as.integer(seed))
    return(sort(sample.int(available, draws)))
  }
  draws <- as.integer(draws)
  if (anyNA(draws) || any(draws < 1L | draws > available)) {
    .nm_stop("Explicit posterior draw indices are outside the saved chain.")
  }
  unique(draws)
}

.nm_bayes_draw <- function(fit, row, n_subjects, n_eta) {
  chain <- fit$chain
  theta_names <- .nm_numbered_names("THETA", length(fit$theta))
  sigma_names <- .nm_numbered_names("SIGMA", length(fit$sigma))
  omega_names <- .nm_numbered_names("OMEGA", length(fit$omega))
  eta_names <- if (n_eta) unlist(lapply(seq_len(n_subjects), function(subject) {
    paste0("ETA", subject, "_", seq_len(n_eta))
  }), use.names = FALSE) else character()
  list(
    theta = as.numeric(chain[row, theta_names, drop = TRUE]),
    sigma = as.numeric(chain[row, sigma_names, drop = TRUE]),
    omega = as.numeric(chain[row, omega_names, drop = TRUE]),
    eta = if (n_eta) matrix(
      as.numeric(chain[row, eta_names, drop = TRUE]),
      n_subjects, n_eta, byrow = TRUE
    ) else matrix(numeric(), n_subjects, 0L)
  )
}

#' Pointwise Bayesian log-likelihood
#'
#' The pointwise unit is one independent subject/likelihood cluster. By
#' default ETAs are integrated from their fitted population distribution for
#' each posterior population draw, which is appropriate for subject-level
#' predictive criteria and avoids conditioning on held-out observations.
#'
#' @param fit A BAYES, HMC, or NUTS `nm_fit`.
#' @param draws Number of posterior draws, explicit draw indices, or `NULL` for
#'   all saved draws.
#' @param marginal Integrate over subject random effects. Set `FALSE` only for
#'   within-subject conditional prediction.
#' @param eta_samples Monte Carlo samples per subject and posterior draw for
#'   marginalization.
#' @param seed RNG seed.
#' @return Draw-by-subject log-likelihood matrix.
#' @export
nm_log_lik <- function(fit, draws = NULL, marginal = TRUE,
                       eta_samples = 64L, seed = 20260713L) {
  context <- .nm_estimation_context(fit$model, fit$data, method = fit$method)
  indices <- .nm_bayes_chain_indices(fit, draws, seed)
  eta_samples <- as.integer(eta_samples)
  if (isTRUE(marginal) && (length(eta_samples) != 1L ||
      is.na(eta_samples) || eta_samples < 2L)) {
    .nm_stop("`eta_samples` must be at least 2.")
  }
  result <- matrix(
    NA_real_, length(indices), context$n_subjects,
    dimnames = list(paste0("draw_", indices),
                    paste0("subject_", unique(context$data$ID)))
  )
  set.seed(as.integer(seed))
  for (draw_index in seq_along(indices)) {
    draw <- .nm_bayes_draw(
      fit, indices[[draw_index]], context$n_subjects, context$n_eta
    )
    for (subject in seq_len(context$n_subjects)) {
      evaluator <- context$subjects[[subject]]
      covariance <- .nm_effect_covariance(
        context$model, evaluator$data, draw$omega
      )
      positive <- .nm_positive_definite(
        covariance, "Bayesian random-effect covariance"
      )
      if (isTRUE(marginal) && context$n_eta) {
        eta <- matrix(
          stats::rnorm(eta_samples * context$n_eta),
          eta_samples, context$n_eta
        ) %*% chol(positive$matrix)
      } else {
        eta <- matrix(draw$eta[subject, ], 1L, context$n_eta)
      }
      joint <- evaluator$objective_eta_values(
        draw$theta, eta, draw$sigma, draw$omega
      )
      prior <- if (context$n_eta) {
        positive$logdet + rowSums(
          (eta %*% solve(positive$matrix)) * eta
        )
      } else 0
      # LibeRation's NONMEM-compatible Gaussian objective omits the common
      # log(2*pi) term. Predictive criteria require an actual log density, so
      # restore that term after removing the (also constant-free) ETA prior.
      observation_constant <- if (!identical(
        context$model$LIK_CONFIG$error, "likelihood"
      )) {
        sum(evaluator$data$EVID == 0L & evaluator$data$MDV == 0L &
              is.finite(evaluator$data$DV)) * log(2 * pi)
      } else 0
      conditional <- -0.5 * (joint - prior + observation_constant)
      result[draw_index, subject] <- if (isTRUE(marginal) && context$n_eta) {
        .nm_log_mean_exp(conditional)
      } else conditional[[1L]]
    }
  }
  attr(result, "pointwise_unit") <- "subject"
  attr(result, "marginal_random_effects") <- isTRUE(marginal)
  attr(result, "eta_samples") <- if (isTRUE(marginal)) eta_samples else 0L
  attr(result, "draw_indices") <- indices
  result
}

#' Widely applicable information criterion
#'
#' @param fit A Bayesian `nm_fit`, or a draw-by-subject log-likelihood matrix.
#' @param log_lik Optional precomputed [nm_log_lik()] matrix.
#' @param ... Passed to [nm_log_lik()] when needed.
#' @return An `nm_waic` object with total and pointwise contributions.
#' @export
nm_waic <- function(fit, log_lik = NULL, ...) {
  if (is.matrix(fit) && is.null(log_lik)) {
    log_lik <- fit
  } else if (is.null(log_lik)) log_lik <- nm_log_lik(fit, ...)
  log_lik <- as.matrix(log_lik)
  if (nrow(log_lik) < 2L || ncol(log_lik) < 1L || any(!is.finite(log_lik))) {
    .nm_stop("WAIC requires a finite log-likelihood matrix with at least two draws.")
  }
  lppd <- apply(log_lik, 2L, .nm_log_mean_exp)
  p_waic <- apply(log_lik, 2L, stats::var)
  elpd <- lppd - p_waic
  pointwise <- data.frame(
    unit = colnames(log_lik) %||% paste0("unit_", seq_len(ncol(log_lik))),
    lppd = lppd, p_waic = p_waic, elpd_waic = elpd,
    waic = -2 * elpd, stringsAsFactors = FALSE
  )
  structure(list(
    waic = sum(pointwise$waic), elpd_waic = sum(elpd),
    p_waic = sum(p_waic),
    se = sqrt(ncol(log_lik) * stats::var(pointwise$waic)),
    pointwise = pointwise, log_lik = log_lik,
    marginal_random_effects = attr(log_lik, "marginal_random_effects")
  ), class = "nm_waic")
}

#' @export
print.nm_waic <- function(x, ...) {
  cat("LibeRation WAIC\n")
  cat("  WAIC:", format(x$waic, digits = 7),
      " SE:", format(x$se, digits = 5),
      " p_WAIC:", format(x$p_waic, digits = 5), "\n")
  invisible(x)
}

#' Pareto-smoothed importance-sampling leave-one-subject-out validation
#'
#' @param fit A Bayesian `nm_fit`, or a log-likelihood matrix.
#' @param log_lik Optional precomputed [nm_log_lik()] matrix.
#' @param ... Passed to [nm_log_lik()] when needed.
#' @return An `nm_psis_loo` wrapper around `loo::loo()`.
#' @export
nm_psis_loo <- function(fit, log_lik = NULL, ...) {
  if (!requireNamespace("loo", quietly = TRUE)) {
    .nm_stop("PSIS-LOO requires the optional `loo` package; install it with install.packages('loo').")
  }
  if (is.matrix(fit) && is.null(log_lik)) {
    log_lik <- fit
  } else if (is.null(log_lik)) log_lik <- nm_log_lik(fit, ...)
  result <- loo::loo(as.matrix(log_lik))
  structure(list(
    loo = result,
    estimates = result$estimates,
    pointwise = result$pointwise,
    pareto_k = loo::pareto_k_values(result),
    log_lik = log_lik,
    marginal_random_effects = attr(log_lik, "marginal_random_effects")
  ), class = "nm_psis_loo")
}

#' @export
print.nm_psis_loo <- function(x, ...) {
  cat("LibeRation PSIS-LOO (subject level)\n")
  print(x$loo)
  invisible(x)
}

.nm_ppc_statistics <- function(value, statistics) {
  vapply(statistics, function(fun) {
    result <- as.numeric(fun(value))
    if (length(result) != 1L) {
      .nm_stop("Each posterior-predictive statistic must return one number.")
    }
    result[[1L]]
  }, numeric(1))
}

#' Bayesian posterior predictive checks
#'
#' @param fit A BAYES, HMC, or NUTS `nm_fit`.
#' @param draws Number of posterior draws or explicit draw indices.
#' @param predictive `population` samples new random effects; `conditional`
#'   reuses each posterior draw's subject effects.
#' @param statistics Named list of scalar functions applied to observed and
#'   replicated outcomes.
#' @param stratify Optional dataset column for separate checks by stratum.
#' @param seed RNG seed.
#' @param keep_data Retain replicated row-level data.
#' @return An `nm_ppc` object.
#' @export
nm_ppc <- function(fit, draws = 200L,
                   predictive = c("population", "conditional"),
                   statistics = NULL, stratify = NULL,
                   seed = 20260713L, keep_data = FALSE) {
  predictive <- match.arg(predictive)
  indices <- .nm_bayes_chain_indices(fit, draws, seed)
  context <- .nm_estimation_context(fit$model, fit$data, method = fit$method)
  statistics <- statistics %||% list(
    mean = mean, sd = stats::sd,
    q05 = function(x) stats::quantile(x, 0.05, names = FALSE),
    median = stats::median,
    q95 = function(x) stats::quantile(x, 0.95, names = FALSE)
  )
  if (!is.list(statistics) || !length(statistics) ||
      is.null(names(statistics)) || any(!vapply(statistics, is.function, logical(1)))) {
    .nm_stop("`statistics` must be a named list of scalar functions.")
  }
  stratify <- as.character(stratify %||% "")
  if (length(stratify) != 1L || (nzchar(stratify) &&
      !stratify %in% names(fit$data))) {
    .nm_stop("`stratify` must name one fitted-data column or be NULL.")
  }
  observed_rows <- fit$data$EVID == 0L & fit$data$MDV == 0L &
    is.finite(fit$data$DV)
  groups <- if (nzchar(stratify)) {
    as.character(fit$data[[stratify]])
  } else rep("overall", nrow(fit$data))
  groups[is.na(groups)] <- "<missing>"
  group_levels <- unique(groups[observed_rows])
  observed <- do.call(rbind, lapply(group_levels, function(group) {
    values <- fit$data$DV[observed_rows & groups == group]
    data.frame(group = group, t(.nm_ppc_statistics(values, statistics)),
               check.names = FALSE, stringsAsFactors = FALSE)
  }))
  replicated <- array(
    NA_real_, c(length(indices), length(group_levels), length(statistics)),
    dimnames = list(paste0("draw_", indices), group_levels, names(statistics))
  )
  data_replicates <- if (isTRUE(keep_data)) vector("list", length(indices)) else NULL
  set.seed(as.integer(seed))
  simulation_seeds <- sample.int(.Machine$integer.max, length(indices))
  for (index in seq_along(indices)) {
    draw <- .nm_bayes_draw(
      fit, indices[[index]], context$n_subjects, context$n_eta
    )
    simulated <- nm_simulate(
      fit$model, fit$data, theta = draw$theta, sigma = draw$sigma,
      omega = draw$omega,
      eta = if (predictive == "conditional") draw$eta else NULL,
      random_effects = predictive == "population", residual = TRUE,
      sample_mixture = TRUE, seed = simulation_seeds[[index]]
    )
    simulated_rows <- simulated$EVID == 0L & simulated$MDV == 0L &
      is.finite(simulated$DV)
    simulated_groups <- if (nzchar(stratify)) {
      as.character(simulated[[stratify]])
    } else rep("overall", nrow(simulated))
    simulated_groups[is.na(simulated_groups)] <- "<missing>"
    for (group in group_levels) {
      values <- simulated$DV[simulated_rows & simulated_groups == group]
      replicated[index, group, ] <- .nm_ppc_statistics(values, statistics)
    }
    if (isTRUE(keep_data)) data_replicates[[index]] <- simulated
  }
  checks <- do.call(rbind, lapply(seq_along(group_levels), function(group_index) {
    data.frame(
      group = group_levels[[group_index]],
      statistic = names(statistics),
      observed = unname(as.numeric(
        observed[group_index, names(statistics), drop = TRUE]
      )),
      replicated_mean = unname(apply(
        replicated[, group_index, , drop = FALSE], 3L, mean, na.rm = TRUE
      )),
      lower = unname(apply(replicated[, group_index, , drop = FALSE], 3L,
                    stats::quantile, 0.025, na.rm = TRUE)),
      upper = unname(apply(replicated[, group_index, , drop = FALSE], 3L,
                    stats::quantile, 0.975, na.rm = TRUE)),
      bayesian_p = vapply(seq_along(statistics), function(statistic) {
        mean(replicated[, group_index, statistic] >=
               observed[group_index, names(statistics)[[statistic]]], na.rm = TRUE)
      }, numeric(1)),
      stringsAsFactors = FALSE
    )
  }))
  structure(list(
    checks = checks, observed = observed, replicated = replicated,
    predictive = predictive, stratify = if (nzchar(stratify)) stratify else NULL,
    draw_indices = indices, seed = seed, data = data_replicates,
    fit_fingerprint = .nm_fit_fingerprint(fit)
  ), class = "nm_ppc")
}

#' @export
print.nm_ppc <- function(x, ...) {
  cat("LibeRation posterior predictive check\n")
  cat("  predictive target:", x$predictive,
      " draws:", length(x$draw_indices), "\n")
  print(x$checks, row.names = FALSE)
  invisible(x)
}
