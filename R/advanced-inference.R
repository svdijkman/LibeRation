.nm_log_mean_exp <- function(x) {
  maximum <- max(x)
  if (!is.finite(maximum)) return(maximum)
  maximum + log(mean(exp(x - maximum)))
}

.nm_fit_observations <- function(fit) {
  data <- fit$data
  sum(data$EVID == 0L & data$MDV == 0L & is.finite(data$DV))
}

.nm_fit_free_parameters <- function(fit) {
  sum(!fit$model$THETAS$FIX) + sum(!fit$model$SIGMAS$FIX) +
    if (any(fit$model$OMEGAS$ROW != fit$model$OMEGAS$COL) &&
        all(fit$model$OMEGAS$FIX)) 0L else sum(!fit$model$OMEGAS$FIX)
}

.nm_fit_fingerprint <- function(fit) {
  digest::digest(list(
    contract = nm_model_to_contract(fit$model),
    data = as.data.frame(fit$data),
    method = fit$method,
    objective = fit$objective,
    parameters = .nm_fit_native_parameters(fit)
  ), algo = "sha256", serialize = TRUE)
}

#' Create a reproducible, first-class model comparison
#'
#' The returned object retains model, dataset, fit and diagnostic fingerprints
#' together with all comparison rules. It is therefore safe to save as a
#' durable analysis artefact, feed into automated covariate workflows, or
#' summarize for the local Help AI.
#'
#' @param ... Fitted `nm_fit` objects, or one list of fits.
#' @param labels Optional unique display labels.
#' @param reference Reference-model index.
#' @param nested Logical flags declaring whether each non-reference candidate
#'   is a nested extension of the reference model. No likelihood-ratio
#'   inference is made unless nesting is explicitly declared.
#' @param diagnostics Optional named list of diagnostic objects for each fit.
#' @param alpha Decision threshold for declared nested comparisons.
#' @param bic_n Sample-size convention used by BIC: observed records or
#'   independent subjects.
#' @param notes Optional analysis notes retained in provenance.
#' @param parametric_bootstrap Number of parametric likelihood-ratio
#'   replications. Zero disables calibration.
#' @param seed Parametric-bootstrap seed.
#' @param refit_control Named list passed to [nm_est()] during bootstrap
#'   calibration.
#' @return An `nm_model_comparison`.
#' @export
nm_compare <- function(..., labels = NULL, reference = 1L, nested = FALSE,
                       diagnostics = NULL, alpha = 0.05,
                       bic_n = c("observations", "subjects"), notes = "",
                       parametric_bootstrap = 0L, seed = 20260713L,
                       refit_control = list()) {
  fits <- list(...)
  if (length(fits) == 1L && is.list(fits[[1L]]) &&
      !inherits(fits[[1L]], "nm_fit")) fits <- fits[[1L]]
  if (length(fits) < 2L || any(!vapply(fits, inherits, logical(1), "nm_fit"))) {
    .nm_stop("`nm_compare()` requires at least two `nm_fit` objects.")
  }
  reference <- as.integer(reference)
  if (length(reference) != 1L || is.na(reference) ||
      reference < 1L || reference > length(fits)) {
    .nm_stop("`reference` must identify one supplied fit.")
  }
  alpha <- as.numeric(alpha)
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    .nm_stop("`alpha` must lie strictly between zero and one.")
  }
  bic_n <- match.arg(bic_n)
  labels <- labels %||% names(fits)
  if (is.null(labels) || length(labels) != length(fits) ||
      any(!nzchar(labels))) labels <- paste0("Model ", seq_along(fits))
  labels <- make.unique(as.character(labels))
  names(fits) <- labels
  same_data <- vapply(fits, function(fit) {
    identical(
      digest::digest(as.data.frame(fit$data), algo = "sha256"),
      digest::digest(as.data.frame(fits[[reference]]$data), algo = "sha256")
    )
  }, logical(1))
  if (!length(nested) || !length(nested) %in% c(1L, length(fits)) ||
      anyNA(as.logical(nested))) {
    .nm_stop("`nested` must be one logical value or one flag per fit.")
  }
  nested <- rep_len(as.logical(nested), length(fits))
  nested[[reference]] <- TRUE
  n_parameters <- vapply(fits, .nm_fit_free_parameters, integer(1))
  n_observations <- vapply(fits, .nm_fit_observations, integer(1))
  objective <- vapply(fits, function(fit) as.numeric(fit$objective), numeric(1))
  methods <- vapply(fits, `[[`, character(1), "method")
  subjects <- vapply(
    fits, function(fit) length(unique(fit$data$ID)), integer(1)
  )
  error_models <- vapply(
    fits, function(fit) fit$model$LIK_CONFIG$error %||% "unknown", character(1)
  )
  comparable <- same_data & methods == methods[[reference]] &
    error_models == error_models[[reference]]
  bic_sample_size <- if (bic_n == "subjects") subjects else n_observations
  metrics <- data.frame(
    model = labels,
    method = methods,
    objective = objective,
    parameters = n_parameters,
    observations = n_observations,
    subjects = subjects,
    comparable = comparable,
    AIC = objective + 2 * n_parameters,
    BIC = objective + log(pmax(bic_sample_size, 1L)) * n_parameters,
    BIC_sample_size = bic_sample_size,
    convergence = vapply(fits, function(fit) as.integer(fit$convergence), integer(1)),
    same_data = same_data,
    fingerprint = vapply(fits, .nm_fit_fingerprint, character(1)),
    stringsAsFactors = FALSE
  )
  metrics$delta_AIC <- metrics$delta_BIC <-
    metrics$AIC_weight <- metrics$BIC_weight <- NA_real_
  metrics$delta_AIC[comparable] <- metrics$AIC[comparable] -
    min(metrics$AIC[comparable])
  metrics$delta_BIC[comparable] <- metrics$BIC[comparable] -
    min(metrics$BIC[comparable])
  metrics$AIC_weight[comparable] <- exp(
    -0.5 * metrics$delta_AIC[comparable]
  )
  metrics$AIC_weight[comparable] <- metrics$AIC_weight[comparable] /
    sum(metrics$AIC_weight[comparable])
  metrics$BIC_weight[comparable] <- exp(
    -0.5 * metrics$delta_BIC[comparable]
  )
  metrics$BIC_weight[comparable] <- metrics$BIC_weight[comparable] /
    sum(metrics$BIC_weight[comparable])
  ref <- fits[[reference]]
  free_names <- function(fit) {
    names(.nm_fit_native_parameters(fit))[c(
      !fit$model$THETAS$FIX, !fit$model$SIGMAS$FIX, !fit$model$OMEGAS$FIX
    )]
  }
  pairwise <- do.call(rbind, lapply(setdiff(seq_along(fits), reference), function(index) {
    compatible <- comparable[[index]]
    df <- n_parameters[[index]] - n_parameters[[reference]]
    delta <- objective[[reference]] - objective[[index]]
    valid_lrt <- compatible && nested[[index]] && df > 0L
    p <- if (valid_lrt) stats::pchisq(max(delta, 0), df = df, lower.tail = FALSE) else NA_real_
    added <- setdiff(free_names(fits[[index]]), free_names(ref))
    boundary <- any(grepl("^(OMEGA|SIGMA)", added)) ||
      !is.null(ref$model$LIK_CONFIG$mixtures) ||
      !is.null(fits[[index]]$model$LIK_CONFIG$mixtures)
    data.frame(
      reference = labels[[reference]], candidate = labels[[index]],
      nested = nested[[index]], compatible_likelihood = compatible,
      delta_objective = delta, delta_parameters = df,
      likelihood_ratio_p = p,
      boundary_warning = boundary,
      warning = if (valid_lrt && boundary) {
        "Ordinary chi-square LRT asymptotics may be invalid at a variance or mixture boundary; use parametric-bootstrap calibration."
      } else "",
      decision = if (!valid_lrt) {
        "descriptive-only"
      } else if (p < alpha) "prefer-candidate" else "retain-reference",
      stringsAsFactors = FALSE
    )
  }))
  parametric_bootstrap <- as.integer(parametric_bootstrap)
  if (length(parametric_bootstrap) != 1L || is.na(parametric_bootstrap) ||
      parametric_bootstrap < 0L) {
    .nm_stop("`parametric_bootstrap` must be a non-negative integer.")
  }
  if (!is.list(refit_control)) .nm_stop("`refit_control` must be a named list.")
  bootstrap <- list()
  if (parametric_bootstrap > 0L) {
    set.seed(as.integer(seed))
    seeds <- sample.int(.Machine$integer.max, parametric_bootstrap)
    fitted_model <- function(fit) .nm_model_rebuild(fit$model, list(
      THETAS = transform(fit$model$THETAS, Value = fit$theta),
      OMEGAS = transform(fit$model$OMEGAS, Value = fit$omega),
      SIGMAS = transform(fit$model$SIGMAS, Value = fit$sigma)
    ))
    reference_model <- fitted_model(ref)
    candidates <- match(
      pairwise$candidate[
        pairwise$nested & pairwise$compatible_likelihood &
          pairwise$delta_parameters > 0L
      ],
      labels
    )
    for (index in candidates) {
      observed_delta <- objective[[reference]] - objective[[index]]
      simulated_delta <- rep(NA_real_, parametric_bootstrap)
      errors <- character(parametric_bootstrap)
      for (replicate in seq_len(parametric_bootstrap)) {
        simulated <- tryCatch(nm_simulate(
          reference_model, ref$data, theta = ref$theta, sigma = ref$sigma,
          omega = ref$omega, random_effects = TRUE, residual = TRUE,
          sample_mixture = TRUE, seed = seeds[[replicate]]
        ), error = identity)
        if (inherits(simulated, "error")) {
          errors[[replicate]] <- conditionMessage(simulated)
          next
        }
        internal <- grep("^\\.", names(simulated), value = TRUE)
        simulated[internal] <- NULL
        base_fit <- tryCatch(do.call(nm_est, c(list(
          model = reference_model, data = simulated, method = ref$method
        ), refit_control)), error = identity)
        candidate_fit <- tryCatch(do.call(nm_est, c(list(
          model = fitted_model(fits[[index]]), data = simulated,
          method = fits[[index]]$method
        ), refit_control)), error = identity)
        if (inherits(base_fit, "error") || inherits(candidate_fit, "error")) {
          errors[[replicate]] <- paste(
            if (inherits(base_fit, "error")) conditionMessage(base_fit),
            if (inherits(candidate_fit, "error")) conditionMessage(candidate_fit),
            collapse = "; "
          )
        } else {
          simulated_delta[[replicate]] <- base_fit$objective - candidate_fit$objective
        }
      }
      finite <- is.finite(simulated_delta)
      bootstrap[[labels[[index]]]] <- list(
        reference = labels[[reference]], candidate = labels[[index]],
        observed_delta = observed_delta, simulated_delta = simulated_delta,
        successful = sum(finite), requested = parametric_bootstrap,
        p_value = if (any(finite)) {
          (1 + sum(simulated_delta[finite] >= observed_delta)) / (1 + sum(finite))
        } else NA_real_,
        errors = errors[nzchar(errors)], seeds = seeds
      )
    }
  }
  parameter_names <- unique(unlist(lapply(fits, function(fit) {
    c(
      .nm_numbered_names("THETA", length(fit$theta)),
      .nm_numbered_names("OMEGA", length(fit$omega)),
      .nm_numbered_names("SIGMA", length(fit$sigma))
    )
  }), use.names = FALSE))
  parameters <- data.frame(parameter = parameter_names, stringsAsFactors = FALSE)
  for (index in seq_along(fits)) {
    values <- .nm_fit_native_parameters(fits[[index]])
    parameters[[labels[[index]]]] <- unname(values[parameter_names])
  }
  if (is.null(diagnostics)) diagnostics <- stats::setNames(vector(
    "list", length(fits)
  ), labels)
  if (!is.list(diagnostics)) .nm_stop("`diagnostics` must be a list.")
  if (is.null(names(diagnostics))) names(diagnostics) <- labels[seq_along(diagnostics)]
  diagnostic_fingerprints <- stats::setNames(lapply(labels, function(label) {
    items <- diagnostics[[label]]
    if (is.null(items)) return(character())
    if (!is.list(items) || is.null(names(items))) items <- list(diagnostic = items)
    vapply(items, digest::digest, character(1), algo = "sha256",
           serialize = TRUE)
  }), labels)
  provenance <- list(
    schema = "liberation.model-comparison", version = 1L,
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    package_version = as.character(utils::packageVersion("LibeRation")),
    reference = labels[[reference]], alpha = alpha, bic_n = bic_n,
    nested = stats::setNames(nested, labels),
    parametric_bootstrap = list(
      requested = parametric_bootstrap, seed = as.integer(seed)
    ),
    rules = list(
      information_criteria = c("AIC", "BIC"),
      likelihood_ratio = "Only for explicitly declared nested models fitted to identical data with compatible likelihoods."
    ),
    notes = as.character(notes)
  )
  evidence <- list(
    fit_fingerprints = stats::setNames(metrics$fingerprint, labels),
    diagnostic_fingerprints = diagnostic_fingerprints,
    metrics = digest::digest(metrics, algo = "sha256", serialize = TRUE),
    pairwise = digest::digest(pairwise, algo = "sha256", serialize = TRUE),
    parameters = digest::digest(parameters, algo = "sha256", serialize = TRUE),
    reference = reference, nested = nested, alpha = alpha, bic_n = bic_n,
    parametric_bootstrap = digest::digest(
      bootstrap, algo = "sha256", serialize = TRUE
    ),
    notes = as.character(notes)
  )
  structure(list(
    id = digest::digest(evidence, algo = "sha256", serialize = TRUE),
    fits = fits, labels = labels, metrics = metrics, pairwise = pairwise,
    parametric_bootstrap = bootstrap,
    parameters = parameters, diagnostics = diagnostics,
    diagnostic_fingerprints = diagnostic_fingerprints, evidence = evidence,
    reference = reference, provenance = provenance
  ), class = "nm_model_comparison")
}

#' @export
print.nm_model_comparison <- function(x, ...) {
  cat("LibeRation reproducible model comparison\n")
  cat("  id:", substr(x$id, 1L, 12L), " reference:",
      x$labels[[x$reference]], "\n")
  print(x$metrics[, setdiff(names(x$metrics), "fingerprint"), drop = FALSE],
        row.names = FALSE)
  invisible(x)
}

#' @export
summary.nm_model_comparison <- function(object, ...) {
  list(
    id = object$id, reference = object$labels[[object$reference]],
    metrics = object$metrics, pairwise = object$pairwise,
    parameters = object$parameters, provenance = object$provenance
  )
}

#' Persist a reproducible model comparison
#'
#' @param comparison An [nm_compare()] result.
#' @param file RDS destination or source.
#' @return The comparison invisibly when saving, and the validated comparison
#'   when reading.
#' @export
nm_compare_save <- function(comparison, file) {
  if (!inherits(comparison, "nm_model_comparison")) {
    .nm_stop("`comparison` must be an nm_model_comparison.")
  }
  saveRDS(comparison, file, version = 3L)
  invisible(comparison)
}

#' @rdname nm_compare_save
#' @export
nm_compare_read <- function(file) {
  value <- readRDS(file)
  if (!inherits(value, "nm_model_comparison") ||
      !identical(value$provenance$schema, "liberation.model-comparison")) {
    .nm_stop("The file does not contain a LibeRation model comparison.")
  }
  current_fits <- stats::setNames(vapply(
    value$fits, .nm_fit_fingerprint, character(1)
  ), value$labels)
  current_diagnostics <- stats::setNames(lapply(value$labels, function(label) {
    items <- value$diagnostics[[label]]
    if (is.null(items)) return(character())
    if (!is.list(items) || is.null(names(items))) items <- list(diagnostic = items)
    vapply(items, digest::digest, character(1), algo = "sha256",
           serialize = TRUE)
  }), value$labels)
  valid <- identical(current_fits, value$evidence$fit_fingerprints) &&
    identical(current_diagnostics, value$evidence$diagnostic_fingerprints) &&
    identical(
      digest::digest(value$metrics, algo = "sha256", serialize = TRUE),
      value$evidence$metrics
    ) &&
    identical(
      digest::digest(value$pairwise, algo = "sha256", serialize = TRUE),
      value$evidence$pairwise
    ) &&
    identical(
      digest::digest(value$parameters, algo = "sha256", serialize = TRUE),
      value$evidence$parameters
    ) &&
    identical(
      digest::digest(
        value$parametric_bootstrap, algo = "sha256", serialize = TRUE
      ),
      value$evidence$parametric_bootstrap
    ) &&
    identical(
      digest::digest(value$evidence, algo = "sha256", serialize = TRUE),
      value$id
    )
  if (!isTRUE(valid)) {
    .nm_stop("The saved model-comparison evidence failed its integrity check.")
  }
  value
}

.nm_sir_mvt_draw <- function(n, center, scale, df) {
  root <- chol(scale)
  z <- matrix(stats::rnorm(n * length(center)), n, length(center)) %*% root
  sweep(z / sqrt(stats::rchisq(n, df = df) / df), 2L, center, "+")
}

.nm_sir_mvt_log_density <- function(x, center, scale, df) {
  x <- as.matrix(x)
  dimension <- length(center)
  root <- chol(scale)
  centered <- sweep(x, 2L, center, "-")
  q <- rowSums((centered %*% solve(root))^2)
  lgamma((df + dimension) / 2) - lgamma(df / 2) -
    0.5 * (dimension * log(df * pi) + 2 * sum(log(diag(root)))) -
    0.5 * (df + dimension) * log1p(q / df)
}

.nm_systematic_resample <- function(weights, n) {
  cumulative <- cumsum(weights)
  positions <- (stats::runif(1L) + 0:(n - 1L)) / n
  findInterval(positions, cumulative) + 1L
}

#' Sampling importance resampling uncertainty
#'
#' A heavy-tailed proposal is centred on the final estimate and scaled from
#' the fitted covariance matrix. The exact LibeRation population objective is
#' evaluated for each proposal, after which normalized importance weights and
#' reproducible systematic resampling produce the retained parameter sample.
#'
#' @param fit A frequentist `nm_fit`.
#' @param n_proposal Number of importance proposals.
#' @param n_resample Number of retained resamples.
#' @param inflation Proposal covariance inflation.
#' @param df Student-t proposal degrees of freedom.
#' @param seed RNG seed.
#' @param covariance Optional completed `nm_covariance`; by default the fit's
#'   covariance is used or calculated.
#' @param eta_maxit Conditional-mode iterations used by objective evaluations.
#' @return An `nm_sir` object.
#' @export
nm_sir <- function(fit, n_proposal = 2000L, n_resample = 1000L,
                   inflation = 1.5, df = 5, seed = 20260713L,
                   covariance = NULL, eta_maxit = NULL) {
  if (!inherits(fit, "nm_fit")) .nm_stop("`fit` must be an nm_fit.")
  if (fit$method %in% c("BAYES", "HMC", "NUTS", "NPML", "NPAG")) {
    .nm_stop("SIR requires a regular frequentist fit with a population objective.")
  }
  n_proposal <- as.integer(n_proposal)
  n_resample <- as.integer(n_resample)
  if (anyNA(c(n_proposal, n_resample)) || n_proposal < 10L ||
      n_resample < 1L || !is.finite(inflation) || inflation <= 0 ||
      !is.finite(df) || df <= 2) {
    .nm_stop("SIR requires at least 10 proposals, positive resamples/inflation, and `df > 2`.")
  }
  covariance <- covariance %||% fit$covariance
  if (is.null(covariance) || !identical(covariance$status, "completed")) {
    covariance <- nm_cov_step(fit, type = "auto", eta_maxit = eta_maxit)
  }
  context <- .nm_estimation_context(fit$model, fit$data, method = fit$method)
  map <- .nm_outer_map(fit$model)
  at <- map$encode(.nm_fit_parameters(fit))
  objective <- .nm_cov_objective(
    fit, context, map, anchor = at,
    eta_maxit = as.integer(eta_maxit %||% fit$diagnostics$eta_maxit %||% 100L),
    tolerance = fit$diagnostics$tolerance %||% 1e-7
  )
  parameters <- map$decode(at)
  transform <- .nm_native_transform_jacobian(fit$model, map, parameters)
  active <- which(rowSums(abs(transform)) > 0)
  jacobian <- transform[active, , drop = FALSE]
  native_covariance <- covariance$covariance
  native_names <- rownames(native_covariance)
  target_names <- .nm_parameter_names(fit$theta, fit$sigma, fit$omega)[active]
  native_covariance <- native_covariance[target_names, target_names, drop = FALSE]
  inverse <- if (nrow(jacobian) == ncol(jacobian)) solve(jacobian) else
    qr.solve(jacobian, diag(nrow(jacobian)))
  outer_covariance <- inverse %*% native_covariance %*% t(inverse)
  outer_covariance <- .nm_positive_definite(
    inflation * (outer_covariance + t(outer_covariance)) / 2,
    "SIR proposal covariance"
  )$matrix
  set.seed(as.integer(seed))
  proposals <- .nm_sir_mvt_draw(n_proposal, at, outer_covariance, df)
  valid <- apply(proposals, 1L, map$in_bounds)
  log_target <- rep(-Inf, n_proposal)
  log_target[valid] <- vapply(which(valid), function(index) {
    value <- tryCatch(objective(proposals[index, ]), error = function(error) Inf)
    if (is.finite(value)) -0.5 * value else -Inf
  }, numeric(1))
  log_proposal <- .nm_sir_mvt_log_density(proposals, at, outer_covariance, df)
  log_weight <- log_target - log_proposal
  maximum <- max(log_weight)
  if (!is.finite(maximum)) .nm_stop("No SIR proposal had finite target density.")
  weights <- exp(log_weight - maximum)
  weights <- weights / sum(weights)
  selected <- .nm_systematic_resample(weights, n_resample)
  natural <- t(vapply(selected, function(index) {
    decoded <- map$decode(proposals[index, ])
    c(decoded$theta, decoded$sigma, decoded$omega)
  }, numeric(length(fit$theta) + length(fit$sigma) + length(fit$omega))))
  colnames(natural) <- .nm_parameter_names(fit$theta, fit$sigma, fit$omega)
  display_names <- c(
    .nm_numbered_names("THETA", length(fit$theta)),
    .nm_numbered_names("OMEGA", length(fit$omega)),
    .nm_numbered_names("SIGMA", length(fit$sigma))
  )
  natural <- natural[, display_names, drop = FALSE]
  summary <- t(vapply(seq_len(ncol(natural)), function(index) {
    c(
      mean = mean(natural[, index]), sd = stats::sd(natural[, index]),
      lower_95 = stats::quantile(
        natural[, index], 0.025, names = FALSE
      ),
      median = stats::median(natural[, index]),
      upper_95 = stats::quantile(
        natural[, index], 0.975, names = FALSE
      )
    )
  }, numeric(5L)))
  rownames(summary) <- colnames(natural)
  structure(list(
    draws = natural, proposal = proposals, weights = weights,
    selected = selected,
    summary = summary,
    diagnostics = list(
      effective_sample_size = 1 / sum(weights^2),
      relative_effective_sample_size = 1 / sum(weights^2) / n_proposal,
      maximum_weight = max(weights),
      finite_proposals = sum(is.finite(log_target)),
      proposal_count = n_proposal, resample_count = n_resample,
      inflation = inflation, df = df, seed = seed
    ),
    covariance = covariance, fit_fingerprint = .nm_fit_fingerprint(fit)
  ), class = "nm_sir")
}

#' @export
print.nm_sir <- function(x, ...) {
  cat("LibeRation sampling importance resampling\n")
  cat("  proposals:", x$diagnostics$proposal_count,
      " retained:", x$diagnostics$resample_count,
      " ESS:", format(x$diagnostics$effective_sample_size, digits = 5), "\n")
  invisible(x)
}
