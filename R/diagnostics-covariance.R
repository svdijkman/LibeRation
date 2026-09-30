# Covariance-step construction and uncertainty diagnostics.
# Split from diagnostics.R as a behaviour-preserving source move.

.nm_outer_names <- function(model, map) {
  c(
    if (length(map$theta_free)) paste0("THETA", map$theta_free) else character(),
    if (length(map$sigma_free)) paste0("log_SIGMA", map$sigma_free) else character(),
    if (length(map$omega_free)) {
      if (map$omega_full) paste0("L_", model$OMEGAS$ROW, "_", model$OMEGAS$COL)
      else paste0("log_OMEGA", map$omega_free)
    } else character()
  )
}

.nm_native_transform_jacobian <- function(model, map, parameters) {
  if (is.function(map$jacobian)) return(map$jacobian(parameters))
  n_native <- nrow(model$THETAS) + nrow(model$SIGMAS) + nrow(model$OMEGAS)
  n_outer <- length(map$start)
  jacobian <- matrix(0, n_native, n_outer)
  cursor <- 0L
  if (length(map$theta_free)) {
    for (index in map$theta_free) {
      cursor <- cursor + 1L
      jacobian[index, cursor] <- 1
    }
  }
  sigma_offset <- nrow(model$THETAS)
  if (length(map$sigma_free)) {
    for (index in map$sigma_free) {
      cursor <- cursor + 1L
      jacobian[sigma_offset + index, cursor] <- parameters$sigma[[index]]
    }
  }
  omega_offset <- sigma_offset + nrow(model$SIGMAS)
  if (!length(map$omega_free)) return(jacobian)
  if (!map$omega_full) {
    for (index in map$omega_free) {
      cursor <- cursor + 1L
      jacobian[omega_offset + index, cursor] <- parameters$omega[[index]]
    }
    return(jacobian)
  }
  covariance <- .nm_omega_matrix(model, parameters$omega)
  lower <- t(chol(covariance))
  for (encoded in seq_len(nrow(model$OMEGAS))) {
    cursor <- cursor + 1L
    row <- model$OMEGAS$ROW[[encoded]]
    column <- model$OMEGAS$COL[[encoded]]
    derivative_lower <- matrix(0, model$n_eta, model$n_eta)
    derivative_lower[row, column] <- if (row == column) lower[row, column] else 1
    derivative <- derivative_lower %*% t(lower) + lower %*% t(derivative_lower)
    for (native in seq_len(nrow(model$OMEGAS))) {
      jacobian[omega_offset + native, cursor] <- derivative[
        model$OMEGAS$ROW[[native]], model$OMEGAS$COL[[native]]
      ]
    }
  }
  jacobian
}

.nm_covariance_eta_start <- function(fit, context) {
  if (!context$n_eta || !is.matrix(fit$eta)) return(NULL)
  expected <- c(context$n_subjects, context$n_eta)
  if (!identical(dim(fit$eta), expected) || any(!is.finite(fit$eta))) return(NULL)
  fit$eta
}

.nm_imp_information_objective <- function(context, map, normals, anchor,
                                          eta_maxit, tolerance,
                                          adaptive = TRUE,
                                          initial_eta = NULL) {
  anchor <- as.numeric(anchor)
  parameters <- map$decode(anchor)
  proposal_started <- proc.time()[["elapsed"]]
  proposals <- .nm_imp_prepare_proposals(
    context, parameters, normals, eta_maxit, tolerance,
    adaptive = adaptive, initial_eta = initial_eta
  )
  if (any(!vapply(proposals, function(proposal) isTRUE(proposal$valid), logical(1)))) {
    .nm_stop("Unable to construct finite importance proposals for covariance.")
  }
  telemetry <- new.env(parent = emptyenv())
  telemetry$proposal_seconds <- unname(proc.time()[["elapsed"]] - proposal_started)
  telemetry$parameter_evaluations <- 0L
  telemetry$cache_hits <- 0L
  telemetry$sample_evaluations <- 0
  telemetry$eta_warm_start <- !is.null(initial_eta) && isTRUE(adaptive)
  telemetry$proposal_mode_iterations <- sum(vapply(
    proposals, function(proposal) as.integer(proposal$mode$iterations %||% 0L),
    integer(1)
  ))
  telemetry$proposal_mode_evaluations <- sum(vapply(
    proposals, function(proposal) as.integer(proposal$mode$evaluations %||% 0L),
    integer(1)
  ))
  cache <- new.env(parent = emptyenv())
  cache$key <- NULL
  evaluate <- function(outer) {
    outer <- as.numeric(outer)
    if (!is.null(cache$key) && identical(cache$key, outer)) {
      telemetry$cache_hits <- telemetry$cache_hits + 1L
      return(cache$result)
    }
    candidate <- map$decode(outer)
    evaluated <- .nm_imp_evaluate_fixed(
      context, candidate, proposals, gradient = TRUE
    )
    native <- lapply(evaluated$states, `[[`, "native_gradient")
    if (any(vapply(native, is.null, logical(1)))) {
      .nm_stop("The fixed-proposal importance gradient is unavailable.")
    }
    native <- do.call(rbind, native)
    n_theta <- length(candidate$theta)
    n_sigma <- length(candidate$sigma)
    n_omega <- length(candidate$omega)
    population_positions <- c(
      seq_len(n_theta), n_theta + context$n_eta + seq_len(n_sigma),
      n_theta + context$n_eta + n_sigma + seq_len(n_omega)
    )
    transform <- map$jacobian(candidate)
    subject_gradient <- native[, population_positions, drop = FALSE] %*% transform
    prior_gradient <- as.vector(
      .nm_prior_nll_native_gradient(context$model, candidate) %*% transform
    )
    result <- list(
      value = evaluated$value + .nm_prior_nll(context$model, candidate),
      gradient = colSums(subject_gradient) + prior_gradient,
      scores = -0.5 * subject_gradient
    )
    telemetry$parameter_evaluations <- telemetry$parameter_evaluations + 1L
    telemetry$sample_evaluations <- telemetry$sample_evaluations + sum(vapply(
      proposals, function(proposal) nrow(proposal$eta), integer(1)
    ))
    cache$key <- outer
    cache$result <- result
    result
  }
  objective <- function(outer) evaluate(outer)$value
  attr(objective, "gradient") <- function(outer) evaluate(outer)$gradient
  attr(objective, "subject_scores") <- function(outer) evaluate(outer)$scores
  attr(objective, "objective_backend") <- "fixed-proposal-importance-score"
  attr(objective, "eta_warm_start") <- telemetry$eta_warm_start
  attr(objective, "telemetry") <- function() list(
    proposal_seconds = telemetry$proposal_seconds,
    eta_warm_start = telemetry$eta_warm_start,
    proposal_mode_iterations = telemetry$proposal_mode_iterations,
    proposal_mode_evaluations = telemetry$proposal_mode_evaluations,
    parameter_evaluations = telemetry$parameter_evaluations,
    cache_hits = telemetry$cache_hits,
    sample_evaluations = telemetry$sample_evaluations,
    proposals = length(proposals), samples = nrow(proposals[[1L]]$eta),
    sampling = proposals[[1L]]$sampling %||% "unknown"
  )
  objective
}

.nm_cov_objective <- function(fit, context, map, normals = NULL,
                              anchor = NULL, eta_maxit = 100L,
                              tolerance = 1e-7, adaptive = TRUE) {
  method <- fit$method
  fitted_eta <- .nm_covariance_eta_start(fit, context)
  deterministic <- switch(
    method, FO = "fo", FOCE = "foce", FOCEI = "focei",
    LAPLACE = "laplace", ITS = "its", NULL
  )
  if (!is.null(deterministic)) {
    # Covariance probes are local to the final fit. Record the persistent
    # population objective at that point, including any fitted conditional
    # modes, rather than at the model's original initials. This
    # keeps optimHess perturbations inside the relevant tape domain and leaves
    # the population object's automatic retaping as a genuine path-change
    # fallback instead of making it recover the whole fitted displacement.
    compiled_map <- map
    compiled_map$start <- as.numeric(anchor %||% map$start)
    compiled <- .nm_cpp_population_objective(
      context, compiled_map, deterministic, eta_maxit, tolerance,
      initial_eta = if (identical(method, "FO")) NULL else fitted_eta
    )
    if (!is.null(compiled$pointer)) {
      result <- function(outer) {
        .liberation_population_objective_value(compiled$pointer, outer)
      }
      attr(result, "gradient") <- function(outer) {
        .liberation_population_objective_gradient(compiled$pointer, outer)
      }
      attr(result, "compiled_objective") <- compiled
      attr(result, "eta_warm_start") <- !identical(method, "FO") &&
        !is.null(fitted_eta)
      return(result)
    }
  }
  if (method == "FO") {
    result <- function(outer) .nm_fo_objective(context, map$decode(outer))
    attr(result, "gradient") <- function(outer) {
      .nm_fo_outer_gradient(context, map, map$decode(outer))
    }
    return(result)
  }
  if (method %in% c("GQ", "IMP", "SAEM")) {
    return(.nm_imp_information_objective(
      context, map, normals, anchor %||% map$start, eta_maxit, tolerance,
      adaptive = adaptive,
      initial_eta = if (isTRUE(adaptive)) fitted_eta else NULL
    ))
  }
  approximation <- switch(
    method, FOCE = "foce", FOCEI = "focei", LAPLACE = "laplace",
    ITS = "its", NULL
  )
  if (is.null(approximation)) {
    .nm_stop("Covariance is available for FO, FOCE, FOCEI, LAPLACE, ITS, GQ, IMP, and SAEM fits.")
  }
  objective <- .nm_nested_objective(
    context, approximation, eta_maxit = eta_maxit, tolerance = tolerance,
    initial_eta = fitted_eta
  )
  result <- function(outer) objective(map$decode(outer))
  attr(result, "gradient") <- function(outer) {
    parameters <- map$decode(outer)
    .nm_nested_outer_gradient(
      context, map, objective, parameters, approximation
    )
  }
  attr(result, "eta_warm_start") <- !is.null(fitted_eta)
  result
}

.nm_numeric_gradient <- function(fn, at, relative_step = 1e-4) {
  at <- as.numeric(at)
  gradient <- numeric(length(at))
  baseline <- NULL
  for (index in seq_along(at)) {
    step <- relative_step * max(abs(at[[index]]), 1)
    upper <- lower <- at
    upper[[index]] <- upper[[index]] + step
    lower[[index]] <- lower[[index]] - step
    high <- fn(upper)
    low <- fn(lower)
    if (is.finite(high) && is.finite(low)) {
      gradient[[index]] <- (high - low) / (2 * step)
    } else {
      if (is.null(baseline)) baseline <- fn(at)
      if (is.finite(high) && is.finite(baseline)) {
        gradient[[index]] <- (high - baseline) / step
      } else if (is.finite(low) && is.finite(baseline)) {
        gradient[[index]] <- (baseline - low) / step
      } else {
        .nm_stop("Unable to evaluate a finite-difference marginal subject score.")
      }
    }
  }
  gradient
}

.nm_marginal_score_information <- function(objective, at,
                                           relative_step = 1e-4) {
  gradient <- attr(objective, "gradient", exact = TRUE)
  if (!is.function(gradient)) {
    .nm_stop("The marginal objective does not expose an estimator score.")
  }
  at <- as.numeric(at)
  dimension <- length(at)
  jacobian <- matrix(NA_real_, dimension, dimension)
  baseline <- NULL
  evaluations <- 0L
  one_sided <- integer()
  for (index in seq_len(dimension)) {
    step <- relative_step * max(abs(at[[index]]), 1)
    upper <- lower <- at
    upper[[index]] <- upper[[index]] + step
    lower[[index]] <- lower[[index]] - step
    high <- tryCatch(as.numeric(gradient(upper)), error = function(error) NULL)
    low <- tryCatch(as.numeric(gradient(lower)), error = function(error) NULL)
    evaluations <- evaluations + 2L
    high_ok <- length(high) == dimension && all(is.finite(high))
    low_ok <- length(low) == dimension && all(is.finite(low))
    if (high_ok && low_ok) {
      jacobian[, index] <- (high - low) / (2 * step)
      next
    }
    if (is.null(baseline)) {
      baseline <- tryCatch(as.numeric(gradient(at)), error = function(error) NULL)
      evaluations <- evaluations + 1L
    }
    baseline_ok <- length(baseline) == dimension && all(is.finite(baseline))
    if (high_ok && baseline_ok) {
      jacobian[, index] <- (high - baseline) / step
      one_sided <- c(one_sided, index)
    } else if (low_ok && baseline_ok) {
      jacobian[, index] <- (baseline - low) / step
      one_sided <- c(one_sided, index)
    } else {
      .nm_stop(
        "Unable to evaluate the marginal score around outer parameter ",
        index, "."
      )
    }
  }
  # Monte-Carlo and fixed-grid score Jacobians need not be exactly symmetric
  # at finite sample size. The observed-information estimate is its symmetric
  # part; retain asymmetry as a diagnostic rather than silently discarding it.
  asymmetry <- max(abs(jacobian - t(jacobian)))
  list(
    matrix = (jacobian + t(jacobian)) / 2,
    raw_jacobian = jacobian,
    evaluations = evaluations,
    one_sided_parameters = one_sided,
    maximum_asymmetry = asymmetry,
    relative_step = relative_step
  )
}

.nm_deterministic_subject_scores <- function(fit, context, map, parameters) {
  transform <- .nm_native_transform_jacobian(fit$model, map, parameters)
  scores <- matrix(0, context$n_subjects, ncol(transform))
  if (fit$method == "FO") {
    native <- .nm_fo_collection_gradient(context$subjects, parameters)
    return(-0.5 * native %*% transform)
  }
  approximation <- switch(
    fit$method, FOCE = "foce", FOCEI = "focei", LAPLACE = "laplace",
    ITS = "its", .nm_stop("Subject scores are unavailable for this method.")
  )
  interaction <- approximation != "foce"
  n_theta <- length(parameters$theta)
  n_sigma <- length(parameters$sigma)
  n_omega <- length(parameters$omega)
  eta_positions <- n_theta + seq_len(context$n_eta)
  population_positions <- c(
    seq_len(n_theta), n_theta + context$n_eta + seq_len(n_sigma),
    n_theta + context$n_eta + n_sigma + seq_len(n_omega)
  )
  for (subject in seq_len(context$n_subjects)) {
    evaluator <- context$subjects[[subject]]
    eta <- fit$eta[subject, ]
    derivative <- evaluator$objective(
      parameters$theta, eta, parameters$sigma, parameters$omega,
      gradient = TRUE, interaction = interaction
    )$gradient
    outer <- as.vector(derivative[population_positions] %*% transform)
    if (approximation != "its" && context$n_eta) {
      mixed <- evaluator$objective_hessian_subset(
        parameters$theta, eta, parameters$sigma, parameters$omega,
        rows = eta_positions,
        columns = c(eta_positions, population_positions),
        interaction = interaction
      )
      eta_hessian <- .nm_positive_definite(
        mixed[, seq_len(context$n_eta), drop = FALSE],
        "Conditional ETA curvature for subject score"
      )$matrix
      cross_native <- mixed[, context$n_eta + seq_along(population_positions),
                            drop = FALSE]
      sensitivity <- -solve(eta_hessian, cross_native %*% transform)
      curvature <- evaluator$curvature(
        parameters$theta, eta, parameters$sigma, parameters$omega,
        approximation, gradient = TRUE
      )$gradient
      outer <- outer + as.vector(
        curvature[population_positions] %*% transform +
          curvature[eta_positions] %*% sensitivity
      )
    }
    scores[subject, ] <- -0.5 * outer
  }
  scores
}

.nm_regularized_information <- function(matrix, tolerance) {
  matrix <- (matrix + t(matrix)) / 2
  eigenvalues <- eigen(matrix, symmetric = TRUE, only.values = TRUE)$values
  floor <- max(max(abs(eigenvalues)), 1) * tolerance
  regularization <- max(0, floor - min(eigenvalues))
  adjusted <- matrix + diag(regularization, nrow(matrix))
  list(
    matrix = matrix, adjusted = adjusted, eigenvalues = eigenvalues,
    floor = floor, regularization = regularization,
    condition = max(abs(eigenvalues)) / max(min(eigenvalues), floor),
    stable = all(is.finite(eigenvalues)) && min(eigenvalues) > floor &&
      max(abs(eigenvalues)) / min(eigenvalues) < 1 / tolerance
  )
}

#' Covariance step for a fitted model
#'
#' @param fit An `nm_fit`.
#' @param type `hessian`/`r` uses the population-objective Hessian,
#'   `opg`/`s` uses the subject-score matrix, and `sandwich` (the
#'   default) uses R-inverse S R-inverse. `auto` prefers a well-conditioned R
#'   matrix and falls back to sandwich or S.
#' @param hessian_backend Hessian implementation. `"auto"` uses the native
#'   CppAD/implicit-mode Hessian for deterministic compiled objectives. GQ,
#'   IMP, and SAEM use an estimator-specific numerical Jacobian of their fixed
#'   quadrature/importance score; this is deliberately distinct from a generic
#'   objective `optimHess` calculation.
#'   Deterministic objectives fall back explicitly to a numerical Jacobian of
#'   the population gradient.
#'   `"cppad"` requires the native path; `"numerical"` is retained as a
#'   comparator and for stochastic marginal objectives.
#' @param convention Information-matrix scaling. `"nonmem"` uses derivatives
#'   of the reported objective function for NONMEM-compatible R and S matrices.
#'   `"likelihood"` uses one-half of those derivatives and retains the former
#'   log-likelihood convention.
#' @param tolerance Positive-definite regularization tolerance.
#' @param samples Target integration budget for IMP/SAEM marginal information.
#'   Low-dimensional ETA integrals use a tensor Gauss--Hermite design near this
#'   budget; higher-dimensional integrals use this many random-normal samples.
#'   The default reuses the IMP fit sample count or uses 200 for SAEM. GQ fits
#'   reuse their estimation grid, tensor order or Smolyak level, and point
#'   limit.
#' @param seed Common-random-number seed used by the random-normal fallback for
#'   IMP/SAEM information.
#' @param eta_maxit Maximum conditional ETA iterations for GQ/IMP/SAEM
#'   information.
#' @return Covariance, correlation, standard errors, relative standard errors,
#'   eigenvalues, conditioning diagnostics, and whether fitted conditional
#'   modes warm-started the covariance calculation.
#' @export
nm_cov_step <- function(fit,
                        type = c("sandwich", "auto", "hessian", "opg", "r", "s"),
                        tolerance = 1e-8,
                        samples = NULL, seed = NULL, eta_maxit = NULL,
                        convention = c("nonmem", "likelihood"),
                        hessian_backend = c("auto", "cppad", "numerical")) {
  if (!inherits(fit, "nm_fit")) .nm_stop("`fit` must be an nm_fit.")
  type <- match.arg(type)
  convention <- match.arg(convention)
  hessian_backend <- match.arg(hessian_backend)
  information_scale <- if (convention == "nonmem") 1 else 0.5
  requested_type <- type
  if (type == "r") type <- "hessian"
  if (type == "s") type <- "opg"
  if (fit$method %in% c("BAYES", "HMC", "NUTS")) {
    .nm_stop(fit$method, " reports posterior uncertainty; a frequentist covariance step is not applicable.")
  }
  if (fit$method %in% c("NPML", "NPAG")) {
    .nm_stop("Nonparametric support uncertainty is non-regular; use bootstrap uncertainty for ", fit$method, ".")
  }
  if (!fit$method %in% c("FO", "FOCE", "FOCEI", "LAPLACE", "ITS", "GQ", "IMP", "SAEM")) {
    .nm_stop("The fitted method does not support a covariance step.")
  }
  tolerance <- as.numeric(tolerance)
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0) {
    .nm_stop("`tolerance` must be one positive finite number.")
  }
  context <- attr(fit, ".estimation_context", exact = TRUE)
  if (is.null(context)) context <- .nm_estimation_context(fit$model, fit$data)
  map <- .nm_outer_map(fit$model)
  parameters <- .nm_fit_parameters(fit)
  at <- map$encode(parameters)
  marginal <- fit$method %in% c("GQ", "IMP", "SAEM")
  normals <- NULL
  marginal_design <- NULL
  if (marginal) {
    eta_maxit <- as.integer(eta_maxit %||% fit$diagnostics$eta_maxit %||% 100L)
    if (length(eta_maxit) != 1L || is.na(eta_maxit) || eta_maxit < 1L) {
      .nm_stop("`eta_maxit` must be one positive integer.")
    }
    if (fit$method == "GQ") {
      marginal_design <- .nm_gq_design(
        context,
        order = fit$diagnostics$quadrature_order %||% 5L,
        max_points = fit$diagnostics$quadrature_max_points %||% 100000L,
        grid = fit$diagnostics$quadrature_grid %||% "tensor",
        level = fit$diagnostics$quadrature_level %||% 3L
      )
      samples <- marginal_design$actual_samples
      seed <- NULL
    } else {
      samples <- as.integer(samples %||%
        if (fit$method == "IMP") fit$diagnostics$n_imp %||% 200L else 200L)
      seed <- as.integer(seed %||% fit$diagnostics$seed %||% 20260713L)
      if (length(samples) != 1L || is.na(samples) || samples < 5L) {
        .nm_stop("`samples` must be one integer greater than or equal to 5.")
      }
      if (length(seed) != 1L || is.na(seed)) .nm_stop("`seed` must be one integer.")
      marginal_design <- .nm_imp_covariance_design(context, samples, seed)
    }
    normals <- marginal_design$normals
  }
  if (!length(at)) {
    empty <- matrix(numeric(), 0L, 0L, dimnames = list(character(), character()))
    return(structure(list(
      status = "completed", type = type, covariance = empty, correlation = empty,
      se = numeric(), rse = numeric(), eigenvalues = numeric(),
      condition = NA_real_, regularization = 0, convention = convention
    ), class = "nm_covariance"))
  }
  need_bread <- type %in% c("auto", "hessian", "sandwich")
  need_meat <- type %in% c("opg", "sandwich")
  bread <- meat <- scores <- NULL
  bread_compiled <- NULL
  bread_error <- meat_error <- NULL
  bread_source <- NULL
  bread_exact <- FALSE
  bread_fallback <- NULL
  marginal_information <- NULL
  objective <- NULL
  if (marginal && (need_bread || need_meat)) {
    objective <- .nm_cov_objective(
      fit, context, map, normals = normals, anchor = at,
      eta_maxit = eta_maxit, tolerance = fit$diagnostics$tolerance %||% 1e-7,
      adaptive = if (fit$method == "GQ") {
        isTRUE(fit$diagnostics$adaptive)
      } else TRUE
    )
  }
  if (need_bread) {
    if (is.null(objective)) {
      objective <- .nm_cov_objective(
        fit, context, map, normals = normals, anchor = at,
        eta_maxit = eta_maxit %||% 100L,
        tolerance = fit$diagnostics$tolerance %||% 1e-7,
        adaptive = if (fit$method == "GQ") {
          isTRUE(fit$diagnostics$adaptive)
        } else TRUE
      )
    }
    bread_compiled <- attr(objective, "compiled_objective", exact = TRUE)
    native_available <- !is.null(bread_compiled$pointer) && !marginal
    if (hessian_backend != "numerical" && native_available) {
      bread <- tryCatch(
        information_scale * .liberation_population_objective_hessian(
          bread_compiled$pointer, at
        ),
        error = function(error) {
          bread_error <<- conditionMessage(error)
          NULL
        }
      )
      if (!is.null(bread)) {
        bread_source <- paste(
          "native CppAD Hessian with implicit conditional-mode,",
          "curvature-determinant, and parameter-transform derivatives"
        )
        bread_exact <- TRUE
      } else if (hessian_backend == "cppad") {
        .nm_stop("Native CppAD R-matrix covariance failed: ", bread_error, ".")
      } else {
        bread_fallback <- bread_error
      }
    } else if (hessian_backend == "cppad") {
      .nm_stop(
        "A native CppAD Hessian is unavailable for this objective backend; ",
        "use `hessian_backend = \"auto\"` or `\"numerical\"`."
      )
    }
    if (is.null(bread) && marginal) {
      marginal_error <- NULL
      marginal_information <- tryCatch(
        .nm_marginal_score_information(objective, at),
        error = function(error) {
          marginal_error <<- conditionMessage(error)
          NULL
        }
      )
      if (!is.null(marginal_information)) {
        bread <- information_scale * marginal_information$matrix
        bread_source <- switch(
          fit$method,
          GQ = paste(
            "fixed-proposal quadrature-score Jacobian",
            "(CppAD node scores; proposal/node derivatives held at the fitted point)"
          ),
          IMP = paste(
            "fixed-proposal importance-score Jacobian",
            "(common integration design at the fitted point)"
          ),
          SAEM = paste(
            "post-fit fixed-proposal importance-score Jacobian",
            "(observed marginal-information approximation)"
          )
        )
      } else {
        bread_error <- paste(
          Filter(nzchar, c(bread_error %||% "", marginal_error %||% "")),
          collapse = "; "
        )
      }
    }
    if (is.null(bread) && !marginal) {
      numerical_error <- NULL
      bread <- tryCatch(
        information_scale * stats::optimHess(
          at, objective, gr = attr(objective, "gradient", exact = TRUE)
        ),
        error = function(error) {
          numerical_error <<- conditionMessage(error)
          NULL
        }
      )
      if (!is.null(bread)) {
        bread_source <- "stats::optimHess numerical Jacobian of the population gradient"
      } else {
        bread_error <- paste(
          Filter(nzchar, c(bread_error %||% "", numerical_error %||% "")),
          collapse = "; "
        )
      }
    }
  }
  if (type == "auto" && (is.null(bread) ||
      !.nm_regularized_information(bread, tolerance)$stable)) {
    need_meat <- TRUE
  }
  if (need_meat) {
    scores <- tryCatch({
      result <- matrix(0, context$n_subjects, length(at))
      if (marginal) {
        subject_scores <- attr(objective, "subject_scores", exact = TRUE)
        if (!is.function(subject_scores)) {
          .nm_stop("Marginal subject scores are unavailable.")
        }
        result <- subject_scores(at)
      } else {
        result <- .nm_deterministic_subject_scores(
          fit, context, map, parameters
        )
      }
      if (convention == "nonmem") result <- 2 * result
      result
    }, error = function(error) {
      meat_error <<- conditionMessage(error)
      NULL
    })
    if (!is.null(scores)) meat <- crossprod(scores)
  }
  if (type == "auto") {
    if (!is.null(bread) && .nm_regularized_information(bread, tolerance)$stable) {
      type <- "hessian"
    } else if (!is.null(bread) && !is.null(meat)) {
      type <- "sandwich"
    } else if (!is.null(meat)) {
      type <- "opg"
    } else {
      .nm_stop(
        "Automatic covariance failed. R: ", bread_error %||% "unavailable",
        "; S: ", meat_error %||% "unavailable", "."
      )
    }
  }
  if (type %in% c("hessian", "sandwich") && is.null(bread)) {
    .nm_stop("R-matrix covariance failed: ", bread_error %||% "unavailable", ".")
  }
  if (type %in% c("opg", "sandwich") && is.null(meat)) {
    .nm_stop("S-matrix covariance failed: ", meat_error %||% "unavailable", ".")
  }
  bread_info <- if (!is.null(bread)) .nm_regularized_information(bread, tolerance) else NULL
  meat_info <- if (!is.null(meat)) .nm_regularized_information(meat, tolerance) else NULL
  if (type == "hessian") {
    outer_covariance <- solve(bread_info$adjusted)
    selected <- bread_info
  } else if (type == "opg") {
    outer_covariance <- solve(meat_info$adjusted)
    selected <- meat_info
  } else {
    bread_inverse <- solve(bread_info$adjusted)
    outer_covariance <- bread_inverse %*% meat %*% bread_inverse
    selected <- bread_info
  }
  transform <- .nm_native_transform_jacobian(fit$model, map, parameters)
  active <- which(rowSums(abs(transform)) > 0)
  covariance <- transform[active, , drop = FALSE] %*%
    outer_covariance %*% t(transform[active, , drop = FALSE])
  native_names <- .nm_parameter_names(
    parameters$theta, parameters$sigma, parameters$omega
  )[active]
  native_estimates <- c(
    parameters$theta, parameters$sigma, parameters$omega
  )[active]
  dimnames(covariance) <- list(native_names, native_names)
  se <- sqrt(diag(covariance))
  correlation <- covariance / outer(se, se)
  structure(list(
    status = "completed", type = type, requested_type = requested_type,
    convention = convention,
    covariance = covariance, correlation = correlation,
    se = stats::setNames(se, native_names),
    rse = stats::setNames(100 * se / pmax(abs(native_estimates), 1e-12), native_names),
    eigenvalues = selected$eigenvalues,
    condition = selected$condition,
    regularization = selected$regularization,
    bread = bread, meat = meat, scores = scores,
    bread_source = bread_source,
    bread_exact = bread_exact,
    bread_fallback = bread_fallback,
    hessian_backend = if (bread_exact) "cppad" else if (!is.null(bread)) "numerical" else NULL,
    bread_condition = bread_info$condition %||% NA_real_,
    meat_condition = meat_info$condition %||% NA_real_,
    bread_regularization = bread_info$regularization %||% NA_real_,
    meat_regularization = meat_info$regularization %||% NA_real_,
    fallback = if (requested_type == "auto") type else NULL,
    objective_backend = attr(objective, "objective_backend", exact = TRUE) %||%
      if (!is.null(bread_compiled$pointer)) {
        "persistent-cpp-population-objective"
      } else "r-orchestrated-population-objective",
    eta_warm_start = isTRUE(attr(objective, "eta_warm_start", exact = TRUE)),
    objective_telemetry = {
      importance_telemetry <- attr(objective, "telemetry", exact = TRUE)
      if (is.function(importance_telemetry)) importance_telemetry()
      else if (!is.null(bread_compiled$pointer)) {
        .liberation_population_objective_telemetry(bread_compiled$pointer)
      } else NULL
    },
    marginal_information = if (marginal) c(
      list(
        estimator = fit$method,
        finite_sample = TRUE,
        exact_cppad_hessian = FALSE
      ),
      marginal_information[c(
        "evaluations", "one_sided_parameters", "maximum_asymmetry",
        "relative_step"
      )]
    ) else NULL,
    samples = if (marginal) samples else NULL,
    actual_samples = if (marginal) marginal_design$actual_samples else NULL,
    sampling = if (marginal) marginal_design$method else NULL,
    quadrature_order = if (marginal) marginal_design$quadrature_order else NULL,
    quadrature_level = if (marginal) marginal_design$quadrature_level else NULL,
    quadrature_grid = if (marginal) marginal_design$resolved_grid else NULL,
    seed = if (marginal) seed else NULL
  ), class = "nm_covariance")
}

