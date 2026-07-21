.nm_fit_parameters <- function(object) {
  list(theta = object$theta, sigma = object$sigma, omega = object$omega)
}

.nm_fit_eta_for_data <- function(object, data, type) {
  n_subjects <- length(unique(data$.ID_INDEX))
  n_eta <- .nm_eta_columns(object$model, data)
  if (type == "population") return(matrix(0, n_subjects, n_eta))
  if (nrow(object$eta) != n_subjects || ncol(object$eta) != n_eta) {
    .nm_stop("Individual predictions require the estimation dataset's subject/occasion layout.")
  }
  object$eta
}

.nm_np_predict <- function(object, data, type) {
  distribution <- object$nonparametric
  supports <- as.matrix(distribution$supports)
  n_subjects <- length(unique(data$.ID_INDEX))
  n_eta <- .nm_eta_columns(object$model, data)
  if (ncol(supports) != n_eta) {
    .nm_stop("The nonparametric support dimension does not match the requested occasion layout.")
  }
  probabilities <- if (type == "population") {
    matrix(distribution$weights, n_subjects, nrow(supports), byrow = TRUE)
  } else {
    value <- as.matrix(distribution$posterior_probabilities)
    if (!identical(dim(value), c(n_subjects, nrow(supports)))) {
      .nm_stop("Individual nonparametric predictions require the estimation dataset's subject layout.")
    }
    value
  }
  predictions <- lapply(seq_len(nrow(supports)), function(index) {
    nm_simulate(
      object$model, data, theta = object$theta,
      eta = matrix(supports[index, ], n_subjects, n_eta, byrow = TRUE),
      sigma = object$sigma, omega = object$omega
    )
  })
  result <- predictions[[1L]]
  generated <- setdiff(names(result), names(data))
  prediction_columns <- generated[
    vapply(result[generated], is.numeric, logical(1)) & !grepl("^ETA[0-9]+$", generated)
  ]
  row_probabilities <- probabilities[data$.ID_INDEX, , drop = FALSE]
  for (column in prediction_columns) {
    values <- do.call(cbind, lapply(predictions, `[[`, column))
    result[[column]] <- rowSums(values * row_probabilities)
  }
  eta <- probabilities %*% supports
  if (n_eta) {
    for (column in seq_len(n_eta)) {
      result[[paste0("ETA", column)]] <- eta[data$.ID_INDEX, column]
    }
  }
  attr(result, "solver") <- attr(predictions[[1L]], "solver")
  attr(result, "state_names") <- attr(predictions[[1L]], "state_names")
  result
}

#' Predictions from a fitted LibeRation model
#'
#' @param object An `nm_fit`.
#' @param newdata Optional event dataset. Individual predictions currently
#'   require the original subject and occasion layout.
#' @param type Individual or population predictions.
#' @param ... Reserved.
#' @return Event data augmented by fitted predictions and state amounts.
#' @export
predict.nm_fit <- function(object, newdata = NULL,
                           type = c("individual", "population"), ...) {
  type <- match.arg(type)
  data <- .nm_engine_data(object$model, newdata %||% object$data)
  if (object$method %in% c("NPML", "NPAG") && !is.null(object$nonparametric)) {
    return(.nm_np_predict(object, data, type))
  }
  eta <- .nm_fit_eta_for_data(object, data, type)
  nm_simulate(
    object$model, data, theta = object$theta, eta = eta,
    sigma = object$sigma, omega = object$omega
  )
}

#' Residual diagnostics for a fitted model
#'
#' @param object An `nm_fit`.
#' @param type Residual column to return.
#' @param ... Reserved.
#' @return A numeric residual vector aligned to event records.
#' @export
residuals.nm_fit <- function(object,
                             type = c("IWRES", "CWRES", "WRES", "IRES", "RES"), ...) {
  type <- match.arg(type)
  nm_gof(object)[[type]]
}

#' Goodness-of-fit table
#'
#' Computes population/individual predictions and residuals on the estimation
#' records. `WRES` and `IWRES` use the configured residual variance. `CWRES`
#' are decorrelated within subject using the exact AD prediction Jacobian and a
#' first-order conditional ETA covariance. Censored observations remain marked
#' and receive `NA` residual diagnostics.
#'
#' @param fit An `nm_fit`.
#' @return A data frame aligned to the fitted event data.
#' @export
nm_gof <- function(fit) {
  if (!inherits(fit, "nm_fit")) .nm_stop("`fit` must be an nm_fit.")
  individual <- predict(fit, type = "individual")
  population <- predict(fit, type = "population")
  output <- as.data.frame(fit$data)
  output$PRED <- population$IPRED
  output$IPRED <- individual$IPRED
  output$RES <- output$DV - output$PRED
  output$IRES <- output$DV - output$IPRED
  if (identical(fit$model$LIK_CONFIG$error, "likelihood")) {
    output[c("RES", "IRES", "WRES", "IWRES", "CWRES")] <- NA_real_
    selected <- fit$model$OUTPUT %||% character()
    generated <- setdiff(selected, names(output))
    for (name in intersect(generated, names(individual))) {
      output[[name]] <- individual[[name]]
    }
    if (!is.null(fit$model$KALMAN_CONFIG)) {
      decoded <- nm_kalman_decode(fit, data = fit$data, type = "individual")
      kalman_columns <- grep("^KF_", names(decoded), value = TRUE)
      for (name in kalman_columns) output[[name]] <- decoded[[name]]
      output$KF_STANDARDIZED_INNOVATION <- output$KF_INNOVATION /
        sqrt(output$KF_INNOVATION_VARIANCE)
      attr(output, "kalman_log_likelihood") <- attr(decoded, "log_likelihood")
      attr(output, "residual_note") <- paste(
        "Ordinary Gaussian WRES/IWRES/CWRES are replaced by one-step-ahead",
        "Kalman innovations; filtered and smoothed latent states are supplied in KF_* columns."
      )
    } else if (!is.null(fit$model$HMM_CONFIG)) {
      decoded <- nm_hmm_decode(fit, data = fit$data, type = "individual")
      hmm_columns <- grep("^HMM_", names(decoded), value = TRUE)
      for (name in hmm_columns) output[[name]] <- decoded[[name]]
      attr(output, "hmm_log_likelihood") <- attr(decoded, "log_likelihood")
      state_model <- if (isTRUE(fit$model$HMM_CONFIG$observed_states)) {
        "an observed continuous-time Markov likelihood"
      } else "a hidden Markov likelihood"
      attr(output, "residual_note") <- paste(
        "Gaussian WRES/IWRES/CWRES are not defined for", state_model,
        "; filtered state probabilities and classifications are supplied in HMM_* columns."
      )
    } else if (!is.null(fit$model$OUTCOMES)) {
      family <- nm_outcome_diagnostics(fit, predictions = individual)
      diagnostic_columns <- setdiff(names(family), names(output))
      for (name in diagnostic_columns) output[[name]] <- family[[name]]
      attr(output, "outcome_summary") <- attr(family, "summary", exact = TRUE)
      attr(output, "residual_note") <- paste(
        "Gaussian WRES/IWRES/CWRES are not defined for the compiled outcome likelihood;",
        "family-specific expected values, Pearson/deviance residuals and scores are supplied."
      )
    } else {
      attr(output, "residual_note") <- paste(
        "Gaussian WRES/IWRES/CWRES are not defined for a user likelihood;",
        "use likelihood-appropriate diagnostics such as categorical or Markov VPCs."
      )
    }
    return(output)
  }
  dvid <- if ("DVID" %in% names(output)) output$DVID else rep(1L, nrow(output))
  pop_variance <- .nm_residual_variance(fit$model, output$PRED, fit$sigma, dvid)
  ind_variance <- .nm_residual_variance(fit$model, output$IPRED, fit$sigma, dvid)
  output$WRES <- output$RES / sqrt(pop_variance)
  output$IWRES <- output$IRES / sqrt(ind_variance)
  unavailable <- output$EVID != 0L | output$MDV != 0L | !is.finite(output$DV)
  if ("CENS" %in% names(output)) unavailable <- unavailable | output$CENS == 1L
  if ("BLQ" %in% names(output)) unavailable <- unavailable | output$BLQ == 1L
  output$CWRES <- NA_real_
  available <- which(!unavailable & is.finite(ind_variance) & ind_variance > 0)
  if (length(available)) {
    n_eta <- ncol(fit$eta)
    if (!n_eta) {
      output$CWRES[available] <- output$IWRES[available]
    } else {
      derivative <- nm_prediction_derivatives(
        fit$model, fit$data, theta = fit$theta, eta = fit$eta,
        sigma = fit$sigma, jacobian = TRUE
      )
      omega <- .nm_effect_covariance(fit$model, fit$data, fit$omega)
      omega_inverse <- tryCatch(solve(omega), error = function(error) {
        solve(omega + diag(1e-10 * max(mean(diag(omega)), 1), nrow(omega)))
      })
      groups <- split(available, output$.ID_INDEX[available])
      for (subject_name in names(groups)) {
        rows <- groups[[subject_name]]
        subject <- as.integer(subject_name)
        columns <- match(paste0("ETA_", subject, "_", seq_len(n_eta)), derivative$domain)
        if (anyNA(columns)) {
          output$CWRES[rows] <- output$IWRES[rows]
          next
        }
        h <- derivative$jacobian[rows, columns, drop = FALSE]
        r <- pmax(ind_variance[rows], .Machine$double.eps)
        posterior <- tryCatch(
          solve(omega_inverse + crossprod(h, h / r)),
          error = function(error) NULL
        )
        if (is.null(posterior)) {
          output$CWRES[rows] <- output$IWRES[rows]
          next
        }
        covariance <- diag(r, nrow = length(rows)) + h %*% posterior %*% t(h)
        scale <- max(mean(diag(covariance)), 1)
        root <- tryCatch(chol(covariance), error = function(error) {
          chol(covariance + diag(1e-10 * scale, nrow(covariance)))
        })
        output$CWRES[rows] <- forwardsolve(t(root), output$IRES[rows])
      }
    }
  }
  output[unavailable, c("RES", "IRES", "WRES", "IWRES", "CWRES")] <- NA_real_
  selected <- fit$model$OUTPUT %||% character()
  generated <- setdiff(selected, names(output))
  for (name in intersect(generated, names(individual))) {
    output[[name]] <- individual[[name]]
  }
  output
}

#' Outcome-appropriate diagnostics
#'
#' Computes endpoint predictions, conditional variance, observed-category
#' probability or hazard, Pearson/deviance residuals, Brier/log scores, and
#' event-model cumulative hazard from a model declared with [nm_outcome()].
#' Unlike Gaussian CWRES, these quantities retain their natural interpretation
#' for categorical, count, joint, Markov, and event-time outcomes.
#'
#' @param fit An `nm_fit` with first-class `OUTCOMES`.
#' @param predictions Optional individual prediction table, used internally to
#'   avoid repeating model propagation.
#' @return A data frame aligned to the estimation records.
#' @export
nm_outcome_diagnostics <- function(fit, predictions = NULL) {
  if (!inherits(fit, "nm_fit")) .nm_stop("`fit` must be an nm_fit.")
  outcomes <- fit$model$OUTCOMES
  if (is.null(outcomes)) .nm_stop("The fitted model has no first-class OUTCOMES declaration.")
  predictions <- predictions %||% predict(fit, type = "individual")
  output <- as.data.frame(fit$data)
  n <- nrow(output)
  output$OUTCOME <- output$FAMILY <- rep(NA_character_, n)
  numeric_columns <- c(
    "EXPECTED", "CONDITIONAL_VARIANCE", "OBSERVED_PROBABILITY", "PRED_CATEGORY",
    "PEARSON_RESIDUAL", "DEVIANCE_RESIDUAL", "BRIER_SCORE", "LOG_SCORE",
    "HAZARD", "CUMULATIVE_HAZARD", "MARTINGALE_RESIDUAL"
  )
  for (name in numeric_columns) output[[name]] <- NA_real_
  summaries <- vector("list", length(outcomes))
  safe_log <- function(value) log(pmax(value, 1e-300))
  for (endpoint in seq_along(outcomes)) {
    outcome <- outcomes[[endpoint]]
    rows <- .nm_outcome_rows(output, outcome, include_mdv = TRUE)
    if (!length(rows)) next
    observed <- as.numeric(output$DV[rows])
    mu <- .nm_outcome_resolve(predictions, outcome$prediction, fit$theta, fit$sigma, rows)
    family <- outcome$family
    expected <- variance <- observed_probability <- pred_category <-
      pearson <- deviance <- brier <- log_score <- rep(NA_real_, length(rows))
    if (family %in% c("normal", "lognormal", "student_t")) {
      scale <- pmax(.nm_outcome_resolve(
        predictions, outcome$scale, fit$theta, fit$sigma, rows
      ), 1e-12)
      if (family == "normal") {
        expected <- mu
        variance <- scale^2
      } else if (family == "lognormal") {
        expected <- pmax(mu, 1e-300) * exp(scale^2 / 2)
        variance <- (exp(scale^2) - 1) * exp(2 * log(pmax(mu, 1e-300)) + scale^2)
      } else {
        expected <- mu
        variance <- if (outcome$df > 2) scale^2 * outcome$df / (outcome$df - 2) else NA_real_
      }
      pearson <- (observed - expected) / sqrt(variance)
      deviance <- pearson
    } else if (family == "bernoulli") {
      expected <- pmin(pmax(mu, 0), 1)
      variance <- expected * (1 - expected)
      observed_probability <- ifelse(observed == 1, expected, 1 - expected)
      pred_category <- as.numeric(expected >= 0.5)
      pearson <- (observed - expected) / sqrt(pmax(variance, 1e-12))
      brier <- (observed - expected)^2
      log_score <- -safe_log(observed_probability)
    } else if (family %in% c("categorical", "ordinal")) {
      probability <- vapply(outcome$probabilities, function(symbol) {
        .nm_outcome_resolve(predictions, symbol, fit$theta, fit$sigma, rows)
      }, numeric(length(rows)))
      probability <- pmax(probability, 0)
      probability <- probability / pmax(rowSums(probability), 1e-300)
      selected <- match(observed, outcome$categories)
      valid <- !is.na(selected)
      observed_probability[valid] <- probability[cbind(which(valid), selected[valid])]
      pred_category <- outcome$categories[max.col(probability, ties.method = "first")]
      brier <- rowSums((probability - vapply(outcome$categories, function(category) {
        as.numeric(observed == category)
      }, numeric(length(rows))))^2)
      log_score <- -safe_log(observed_probability)
    } else if (family %in% c("poisson", "negative_binomial", "binomial",
                             "zero_inflated_poisson", "hurdle_poisson")) {
      expected <- pmax(mu, 0)
      variance <- expected
      if (family == "negative_binomial") {
        size <- pmax(.nm_outcome_resolve(
          predictions, outcome$dispersion, fit$theta, fit$sigma, rows
        ), 1e-12)
        variance <- expected + expected^2 / size
      } else if (family == "binomial") {
        trials <- pmax(.nm_outcome_resolve(
          predictions, outcome$trials, fit$theta, fit$sigma, rows
        ), 0)
        probability <- pmin(pmax(mu, 0), 1)
        expected <- trials * probability
        variance <- trials * probability * (1 - probability)
      } else if (family %in% c("zero_inflated_poisson", "hurdle_poisson")) {
        zero <- pmin(pmax(.nm_outcome_resolve(
          predictions, outcome$zero_probability, fit$theta, fit$sigma, rows
        ), 0), 1)
        if (family == "zero_inflated_poisson") {
          expected <- (1 - zero) * expected
          variance <- (1 - zero) * mu * (1 + zero * mu)
        } else {
          positive_mean <- mu / pmax(1 - exp(-mu), 1e-12)
          expected <- (1 - zero) * positive_mean
          second <- (mu + mu^2) / pmax(1 - exp(-mu), 1e-12)
          variance <- (1 - zero) * second - expected^2
        }
      }
      pearson <- (observed - expected) / sqrt(pmax(variance, 1e-12))
      deviance <- sign(observed - expected) * sqrt(pmax(
        2 * (ifelse(observed > 0, observed * log(observed / pmax(expected, 1e-12)), 0) -
               (observed - expected)), 0
      ))
    } else if (family %in% c("tte", "recurrent_event", "competing_risks")) {
      hazard <- if (family == "competing_risks") {
        rowSums(vapply(outcome$cause_hazards, function(symbol) {
          .nm_outcome_resolve(predictions, symbol, fit$theta, fit$sigma, rows)
        }, numeric(length(rows))))
      } else mu
      hazard <- pmax(hazard, 0)
      output$HAZARD[rows] <- hazard
      groups <- split(seq_along(rows), interaction(
        output$.ID_INDEX[rows], if ("DVID" %in% names(output)) output$DVID[rows] else 1,
        drop = TRUE, lex.order = TRUE
      ))
      cumulative <- rep(NA_real_, length(rows))
      for (group in groups) {
        order_index <- group[order(output$TIME[rows[group]])]
        time <- output$TIME[rows[order_index]]
        cumulative[order_index] <- cumsum(hazard[order_index] * c(0, pmax(diff(time), 0)))
      }
      event <- if (family == "competing_risks") as.numeric(observed != 0) else
        as.numeric(observed == outcome$event)
      output$CUMULATIVE_HAZARD[rows] <- cumulative
      output$MARTINGALE_RESIDUAL[rows] <- event - cumulative
      expected <- hazard
    } else if (family %in% c("markov", "continuous_time_markov")) {
      groups <- split(seq_along(rows), interaction(
        output$.ID_INDEX[rows], if ("DVID" %in% names(output)) output$DVID[rows] else 1,
        drop = TRUE, lex.order = TRUE
      ))
      probability <- matrix(NA_real_, nrow = length(rows), ncol = length(outcome$categories))
      for (group in groups) {
        ordered <- group[order(output$TIME[rows[group]])]
        for (position in seq_along(ordered)) {
          local <- ordered[[position]]
          if (position == 1L) {
            probability[local, ] <- vapply(outcome$initial, function(symbol) {
              .nm_outcome_resolve(predictions, symbol, fit$theta, fit$sigma, rows[[local]])
            }, numeric(1))
          } else {
            previous <- match(observed[ordered[[position - 1L]]], outcome$categories)
            if (family == "markov") {
              probability[local, ] <- vapply(outcome$transition[previous, ], function(symbol) {
                .nm_outcome_resolve(predictions, symbol, fit$theta, fit$sigma, rows[[local]])
              }, numeric(1))
            } else {
              rates <- vapply(outcome$rates, function(symbol) {
                .nm_outcome_resolve(predictions, symbol, fit$theta, fit$sigma, rows[[local]])
              }, numeric(1))
              total <- max(sum(rates), 1e-12)
              dt <- max(output$TIME[rows[[local]]] -
                          output$TIME[rows[[ordered[[position - 1L]]]]], 0)
              p01 <- rates[[1L]] / total * (1 - exp(-total * dt))
              p10 <- rates[[2L]] / total * (1 - exp(-total * dt))
              probability[local, ] <- if (previous == 1L) c(1 - p01, p01) else c(p10, 1 - p10)
            }
          }
        }
      }
      probability <- pmax(probability, 0)
      probability <- probability / pmax(rowSums(probability), 1e-300)
      selected <- match(observed, outcome$categories)
      valid <- !is.na(selected)
      observed_probability[valid] <- probability[cbind(which(valid), selected[valid])]
      pred_category <- outcome$categories[max.col(probability, ties.method = "first")]
      brier <- rowSums((probability - vapply(outcome$categories, function(category) {
        as.numeric(observed == category)
      }, numeric(length(rows))))^2)
      log_score <- -safe_log(observed_probability)
    }
    output$OUTCOME[rows] <- outcome$name
    output$FAMILY[rows] <- family
    output$EXPECTED[rows] <- expected
    output$CONDITIONAL_VARIANCE[rows] <- variance
    output$OBSERVED_PROBABILITY[rows] <- observed_probability
    output$PRED_CATEGORY[rows] <- pred_category
    output$PEARSON_RESIDUAL[rows] <- pearson
    output$DEVIANCE_RESIDUAL[rows] <- deviance
    output$BRIER_SCORE[rows] <- brier
    output$LOG_SCORE[rows] <- log_score
    summaries[[endpoint]] <- data.frame(
      outcome = outcome$name, family = family, records = length(rows),
      mean_log_score = if (any(is.finite(log_score))) mean(log_score, na.rm = TRUE) else NA_real_,
      mean_brier_score = if (any(is.finite(brier))) mean(brier, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  attr(output, "summary") <- do.call(rbind, summaries[lengths(summaries) > 0L])
  class(output) <- c("nm_outcome_diagnostics", class(output))
  output
}

.nm_fit_selected_outputs <- function(fit) {
  selected <- fit$model$OUTPUT %||% character()
  if (!length(selected)) return(NULL)
  table <- nm_gof(fit)
  available <- intersect(selected, names(table))
  result <- data.frame(.ROW = seq_len(nrow(table)), check.names = FALSE)
  for (name in available) result[[name]] <- table[[name]]
  missing <- setdiff(selected, available)
  if (length(missing)) {
    for (name in missing) result[[name]] <- NA_real_
    attr(result, "unavailable") <- missing
  }
  result
}

#' Empirical Bayes estimates and shrinkage
#'
#' @param fit An `nm_fit`.
#' @return A list with subject ETA table and component-wise shrinkage.
#' @export
nm_etab <- function(fit) {
  if (!inherits(fit, "nm_fit")) .nm_stop("`fit` must be an nm_fit.")
  eta <- as.data.frame(fit$eta)
  eta$ID <- attr(fit$data, "id_levels") %||% unique(fit$data$ID)
  eta <- eta[c("ID", setdiff(names(eta), "ID"))]
  covariance <- .nm_effect_covariance(fit$model, fit$data, fit$omega)
  shrinkage <- if (ncol(fit$eta)) {
    1 - apply(fit$eta, 2, stats::sd) / sqrt(diag(covariance))
  } else numeric()
  names(shrinkage) <- colnames(fit$eta)
  list(eta = eta, shrinkage = shrinkage)
}

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

.nm_imp_information_objective <- function(context, map, normals, anchor,
                                          eta_maxit, tolerance,
                                          adaptive = TRUE) {
  anchor <- as.numeric(anchor)
  parameters <- map$decode(anchor)
  proposal_started <- proc.time()[["elapsed"]]
  proposals <- .nm_imp_prepare_proposals(
    context, parameters, normals, eta_maxit, tolerance,
    adaptive = adaptive
  )
  if (any(!vapply(proposals, function(proposal) isTRUE(proposal$valid), logical(1)))) {
    .nm_stop("Unable to construct finite importance proposals for covariance.")
  }
  telemetry <- new.env(parent = emptyenv())
  telemetry$proposal_seconds <- unname(proc.time()[["elapsed"]] - proposal_started)
  telemetry$parameter_evaluations <- 0L
  telemetry$cache_hits <- 0L
  telemetry$sample_evaluations <- 0
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
  attr(objective, "telemetry") <- function() list(
    proposal_seconds = telemetry$proposal_seconds,
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
  deterministic <- switch(
    method, FO = "fo", FOCE = "foce", FOCEI = "focei",
    LAPLACE = "laplace", ITS = "its", NULL
  )
  if (!is.null(deterministic)) {
    compiled <- .nm_cpp_population_objective(
      context, map, deterministic, eta_maxit, tolerance
    )
    if (!is.null(compiled$pointer)) {
      result <- function(outer) {
        .liberation_population_objective_value(compiled$pointer, outer)
      }
      attr(result, "gradient") <- function(outer) {
        .liberation_population_objective_gradient(compiled$pointer, outer)
      }
      attr(result, "compiled_objective") <- compiled
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
      adaptive = adaptive
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
    context, approximation, eta_maxit = eta_maxit, tolerance = tolerance
  )
  result <- function(outer) objective(map$decode(outer))
  attr(result, "gradient") <- function(outer) {
    parameters <- map$decode(outer)
    .nm_nested_outer_gradient(
      context, map, objective, parameters, approximation
    )
  }
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
#' @param type `hessian`/`r` uses objective curvature, `opg`/`s` uses the
#'   subject-score matrix, and `sandwich` uses R-inverse S R-inverse. `auto`
#'   prefers a well-conditioned R matrix and falls back to sandwich or S.
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
#'   eigenvalues, and conditioning diagnostics.
#' @export
nm_cov_step <- function(fit,
                        type = c("auto", "hessian", "opg", "sandwich", "r", "s"),
                        tolerance = 1e-8,
                        samples = NULL, seed = NULL, eta_maxit = NULL) {
  if (!inherits(fit, "nm_fit")) .nm_stop("`fit` must be an nm_fit.")
  type <- match.arg(type)
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
      condition = NA_real_, regularization = 0
    ), class = "nm_covariance"))
  }
  need_bread <- type %in% c("auto", "hessian", "sandwich")
  need_meat <- type %in% c("opg", "sandwich")
  bread <- meat <- scores <- NULL
  bread_compiled <- NULL
  bread_error <- meat_error <- NULL
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
    bread <- tryCatch(
      0.5 * stats::optimHess(
        at, objective, gr = attr(objective, "gradient", exact = TRUE)
      ),
      error = function(error) {
        bread_error <<- conditionMessage(error)
        NULL
      }
    )
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
    covariance = covariance, correlation = correlation,
    se = stats::setNames(se, native_names),
    rse = stats::setNames(100 * se / pmax(abs(native_estimates), 1e-12), native_names),
    eigenvalues = selected$eigenvalues,
    condition = selected$condition,
    regularization = selected$regularization,
    bread = bread, meat = meat, scores = scores,
    bread_condition = bread_info$condition %||% NA_real_,
    meat_condition = meat_info$condition %||% NA_real_,
    bread_regularization = bread_info$regularization %||% NA_real_,
    meat_regularization = meat_info$regularization %||% NA_real_,
    fallback = if (requested_type == "auto") type else NULL,
    objective_backend = attr(objective, "objective_backend", exact = TRUE) %||%
      if (!is.null(bread_compiled$pointer)) {
        "persistent-cpp-population-objective"
      } else "r-orchestrated-population-objective",
    objective_telemetry = {
      importance_telemetry <- attr(objective, "telemetry", exact = TRUE)
      if (is.function(importance_telemetry)) importance_telemetry()
      else if (!is.null(bread_compiled$pointer)) {
        .liberation_population_objective_telemetry(bread_compiled$pointer)
      } else NULL
    },
    samples = if (marginal) samples else NULL,
    actual_samples = if (marginal) marginal_design$actual_samples else NULL,
    sampling = if (marginal) marginal_design$method else NULL,
    quadrature_order = if (marginal) marginal_design$quadrature_order else NULL,
    quadrature_level = if (marginal) marginal_design$quadrature_level else NULL,
    quadrature_grid = if (marginal) marginal_design$resolved_grid else NULL,
    seed = if (marginal) seed else NULL
  ), class = "nm_covariance")
}

#' Visual predictive check summaries
#'
#' @param fit An `nm_fit`.
#' @param nsim Number of stochastic replicates.
#' @param breaks Optional time-bin breaks. The default uses quantile bins.
#' @param probs Observation quantiles.
#' @param level Simulation interval for each quantile.
#' @param seed RNG seed.
#' @param pc_correct Apply the legacy prediction correction `DV * PRED / IPRED`
#'   to observed and simulated values before bin summaries.
#' @param stratify Optional dataset column used to create additional
#'   stratum-specific VPC summaries. The unstratified VPC is always retained.
#' @return An `nm_vpc` list containing observed and simulated summaries.
#' @export
nm_vpc <- function(fit, nsim = 200L, breaks = NULL,
                   probs = c(0.05, 0.5, 0.95), level = 0.9,
                   seed = 20260713L, pc_correct = FALSE, stratify = NULL) {
  if (!inherits(fit, "nm_fit")) .nm_stop("`fit` must be an nm_fit.")
  pc_correct <- isTRUE(pc_correct)
  stratify <- as.character(stratify %||% "")
  if (length(stratify) != 1L || is.na(stratify)) {
    .nm_stop("`stratify` must be one dataset column name or `NULL`.")
  }
  if (!nzchar(stratify)) stratify <- NULL
  if (!is.null(stratify) && !stratify %in% names(fit$data)) {
    .nm_stop("VPC stratification column `", stratify, "` is not present in the estimation data.")
  }
  source_data <- if (pc_correct) nm_gof(fit) else as.data.frame(fit$data)
  observed <- source_data[source_data$EVID == 0L & source_data$MDV == 0L &
                            is.finite(source_data$DV), , drop = FALSE]
  if (pc_correct) {
    valid <- is.finite(observed$PRED) & is.finite(observed$IPRED) &
      abs(observed$IPRED) > sqrt(.Machine$double.eps)
    observed$DV[valid] <- observed$DV[valid] * observed$PRED[valid] / observed$IPRED[valid]
    observed$DV[!valid] <- NA_real_
  }
  if (!nrow(observed)) .nm_stop("VPC requires observed DV records.")
  if (is.null(breaks)) {
    breaks <- unique(stats::quantile(observed$TIME, probs = seq(0, 1, length.out = 7), na.rm = TRUE))
    if (length(breaks) < 2L) breaks <- range(observed$TIME) + c(-0.5, 0.5)
  }
  simulated <- nm_simulate(
    fit$model, fit$data, theta = fit$theta, sigma = fit$sigma,
    omega = fit$omega, nsim = nsim, random_effects = TRUE,
    residual = TRUE, sample_mixture = TRUE, seed = seed
  )
  if (pc_correct) {
    population <- predict(fit, type = "population")$IPRED
    reference <- rep(population, times = as.integer(nsim))
    valid <- is.finite(reference) & is.finite(simulated$IPRED) &
      abs(simulated$IPRED) > sqrt(.Machine$double.eps)
    simulated$DV[valid] <- simulated$DV[valid] * reference[valid] / simulated$IPRED[valid]
    simulated$DV[!valid] <- NA_real_
  }
  simulated <- simulated[simulated$EVID == 0L & simulated$MDV == 0L, , drop = FALSE]
  qnames <- paste0("Q", formatC(100 * probs, format = "fg"))
  summarize_set <- function(observed_set, simulated_set) {
    observed_set$BIN <- cut(observed_set$TIME, breaks = breaks, include.lowest = TRUE)
    simulated_set$BIN <- cut(simulated_set$TIME, breaks = breaks, include.lowest = TRUE)
    summarize <- function(frame) {
      if (!nrow(frame)) return(NULL)
      values <- stats::quantile(frame$DV, probs = probs, na.rm = TRUE, names = FALSE)
      data.frame(
        BIN = as.character(frame$BIN[[1L]]),
        TIME = stats::median(frame$TIME, na.rm = TRUE), N = nrow(frame),
        stats::setNames(as.list(values), qnames), check.names = FALSE
      )
    }
    observed_summary <- do.call(
      rbind, lapply(split(observed_set, observed_set$BIN, drop = TRUE), summarize)
    )
    if (!"SIM" %in% names(simulated_set)) simulated_set$SIM <- 1L
    per_sim <- do.call(rbind, lapply(split(
      simulated_set, list(simulated_set$SIM, simulated_set$BIN), drop = TRUE
    ), function(frame) {
      if (!nrow(frame)) return(NULL)
      cbind(SIM = frame$SIM[[1L]], summarize(frame))
    }))
    alpha <- (1 - level) / 2
    intervals <- do.call(rbind, lapply(split(per_sim, per_sim$BIN), function(frame) {
      values <- lapply(qnames, function(name) {
        stats::quantile(
          frame[[name]], c(alpha, 0.5, 1 - alpha), na.rm = TRUE, names = FALSE
        )
      })
      row <- data.frame(
        BIN = frame$BIN[[1L]], TIME = stats::median(frame$TIME, na.rm = TRUE)
      )
      for (i in seq_along(qnames)) {
        row[[paste0(qnames[[i]], "_lo")]] <- values[[i]][[1L]]
        row[[paste0(qnames[[i]], "_median")]] <- values[[i]][[2L]]
        row[[paste0(qnames[[i]], "_hi")]] <- values[[i]][[3L]]
      }
      row
    }))
    list(
      observed = observed_summary,
      simulated = intervals,
      points = observed_set[
        is.finite(observed_set$TIME) & is.finite(observed_set$DV),
        c("TIME", "DV"), drop = FALSE
      ],
      per_simulation = per_sim
    )
  }
  total <- summarize_set(observed, simulated)
  stratified <- list()
  if (!is.null(stratify)) {
    strata <- unique(as.character(observed[[stratify]]))
    strata <- strata[!is.na(strata) & nzchar(strata)]
    simulated_strata <- as.character(simulated[[stratify]])
    observed_strata <- as.character(observed[[stratify]])
    stratified <- unname(lapply(strata, function(stratum) {
      observed_set <- observed[!is.na(observed_strata) & observed_strata == stratum, , drop = FALSE]
      simulated_set <- simulated[!is.na(simulated_strata) & simulated_strata == stratum, , drop = FALSE]
      if (!nrow(observed_set) || !nrow(simulated_set)) return(NULL)
      c(list(level = stratum), summarize_set(observed_set, simulated_set))
    }))
    stratified <- Filter(Negate(is.null), stratified)
  }
  structure(c(
    total,
    list(
      breaks = breaks, probs = probs, level = level, nsim = nsim, seed = seed,
      pc_correct = pc_correct, stratify = stratify, stratified = stratified
    )
  ), class = "nm_vpc")
}

.nm_predictive_simulation_matrix <- function(fit, nsim, seed) {
  if (!inherits(fit, "nm_fit")) .nm_stop("`fit` must be an nm_fit.")
  nsim <- as.integer(nsim)
  if (length(nsim) != 1L || is.na(nsim) || nsim < 20L) {
    .nm_stop("`nsim` must be an integer of at least 20.")
  }
  simulated <- nm_simulate(
    fit$model, fit$data, theta = fit$theta, sigma = fit$sigma,
    omega = fit$omega, nsim = nsim, random_effects = TRUE,
    residual = TRUE, sample_mixture = TRUE, seed = seed
  )
  records <- nrow(fit$data)
  if (nrow(simulated) != records * nsim) {
    .nm_stop("Predictive simulation rows do not align with the fitted dataset.")
  }
  matrix(as.numeric(simulated$DV), nrow = records, ncol = nsim)
}

#' Numerical predictive check
#'
#' Computes the empirical predictive percentile of every uncensored
#' observation under repeated simulations from the fitted population model.
#'
#' @param fit An `nm_fit`.
#' @param nsim Number of predictive simulations; at least 20.
#' @param seed Reproducible RNG seed.
#' @return An `nm_npc` object with record-level percentiles and tail flags.
#' @export
nm_npc <- function(fit, nsim = 200L, seed = 20260713L) {
  simulations <- .nm_predictive_simulation_matrix(fit, nsim, seed)
  data <- as.data.frame(fit$data)
  observed <- data$EVID == 0L & data$MDV == 0L & is.finite(data$DV)
  if ("CENS" %in% names(data)) observed <- observed & data$CENS != 1L
  if ("BLQ" %in% names(data)) observed <- observed & data$BLQ != 1L
  rows <- which(observed)
  if (!length(rows)) .nm_stop("NPC requires uncensored observed DV records.")
  count <- vapply(rows, function(row) sum(simulations[row, ] <= data$DV[[row]], na.rm = TRUE), numeric(1))
  percentile <- (count + 0.5) / (ncol(simulations) + 1)
  table <- data.frame(
    ROW = rows, ID = data$ID[rows], TIME = data$TIME[rows], DV = data$DV[rows],
    PERCENTILE = percentile,
    TAIL_PROBABILITY = 2 * pmin(percentile, 1 - percentile),
    OUTSIDE_90 = percentile < 0.05 | percentile > 0.95,
    stringsAsFactors = FALSE
  )
  structure(list(
    table = table, nsim = ncol(simulations), seed = as.integer(seed),
    outside_90 = mean(table$OUTSIDE_90),
    histogram = graphics::hist(
      table$PERCENTILE, breaks = seq(0, 1, length.out = 11), plot = FALSE
    )
  ), class = "nm_npc")
}

#' Normalized prediction distribution errors
#'
#' Within each subject, simulations are centered and decorrelated with their
#' empirical predictive covariance before record-wise predictive ranks are
#' transformed to standard-normal NPDE values.
#'
#' @param fit An `nm_fit`.
#' @param nsim Number of predictive simulations; at least 20.
#' @param seed Reproducible RNG seed.
#' @param ridge Relative covariance ridge used for numerically singular blocks.
#' @return An `nm_npde` object with record-level NPDE and summary moments.
#' @export
nm_npde <- function(fit, nsim = 200L, seed = 20260713L, ridge = 1e-8) {
  simulations <- .nm_predictive_simulation_matrix(fit, nsim, seed)
  data <- as.data.frame(fit$data)
  observed <- data$EVID == 0L & data$MDV == 0L & is.finite(data$DV)
  if ("CENS" %in% names(data)) observed <- observed & data$CENS != 1L
  if ("BLQ" %in% names(data)) observed <- observed & data$BLQ != 1L
  rows <- which(observed)
  if (!length(rows)) .nm_stop("NPDE requires uncensored observed DV records.")
  ridge <- as.numeric(ridge)
  if (length(ridge) != 1L || !is.finite(ridge) || ridge <= 0) {
    .nm_stop("`ridge` must be a positive finite scalar.")
  }
  percentile <- rep(NA_real_, length(rows))
  groups <- split(seq_along(rows), data$.ID_INDEX[rows])
  for (indices in groups) {
    record_rows <- rows[indices]
    block <- simulations[record_rows, , drop = FALSE]
    center <- rowMeans(block)
    covariance <- if (nrow(block) == 1L) {
      matrix(stats::var(drop(block)), 1L, 1L)
    } else stats::cov(t(block))
    scale <- max(mean(diag(covariance)), 1)
    covariance <- covariance + diag(ridge * scale, nrow(covariance))
    root <- tryCatch(chol(covariance), error = function(e) NULL)
    if (is.null(root)) root <- chol(covariance + diag(sqrt(ridge) * scale, nrow(covariance)))
    observed_white <- forwardsolve(t(root), data$DV[record_rows] - center)
    simulated_white <- forwardsolve(t(root), sweep(block, 1L, center, "-"))
    percentile[indices] <- vapply(seq_along(indices), function(position) {
      (sum(simulated_white[position, ] <= observed_white[[position]]) + 0.5) /
        (ncol(simulated_white) + 1)
    }, numeric(1))
  }
  percentile <- pmin(pmax(percentile, 1e-12), 1 - 1e-12)
  npde <- stats::qnorm(percentile)
  centered <- npde - mean(npde)
  spread <- stats::sd(npde)
  table <- data.frame(
    ROW = rows, ID = data$ID[rows], TIME = data$TIME[rows], DV = data$DV[rows],
    PERCENTILE = percentile, NPDE = npde, stringsAsFactors = FALSE
  )
  structure(list(
    table = table, nsim = ncol(simulations), seed = as.integer(seed), ridge = ridge,
    summary = c(
      mean = mean(npde), sd = spread,
      skewness = if (is.finite(spread) && spread > 0) mean(centered^3) / spread^3 else NA_real_,
      kurtosis = if (is.finite(spread) && spread > 0) mean(centered^4) / spread^4 - 3 else NA_real_
    )
  ), class = "nm_npde")
}

.nm_vpc_breaks <- function(time, breaks = NULL) {
  if (!is.null(breaks)) return(as.numeric(breaks))
  breaks <- unique(stats::quantile(time, probs = seq(0, 1, length.out = 7), na.rm = TRUE))
  if (length(breaks) < 2L) breaks <- range(time, na.rm = TRUE) + c(-0.5, 0.5)
  breaks
}

#' Categorical visual predictive check
#'
#' First-class Bernoulli, categorical, ordinal, and Markov outcomes use their
#' declared probability vectors and can contain any number of categories. A
#' legacy binary user likelihood remains supported when `F` is the probability
#' of the non-reference category.
#'
#' @param fit An `nm_fit`.
#' @param outcome Binary observed outcome column.
#' @param nsim Number of simulations.
#' @param breaks Optional time bins.
#' @param level Simulation interval.
#' @param seed RNG seed.
#' @return An `nm_vpc_categorical`.
#' @export
nm_vpc_categorical <- function(fit, outcome = "DV", nsim = 200L, breaks = NULL,
                               level = 0.9, seed = 20260713L) {
  if (!inherits(fit, "nm_fit")) .nm_stop("`fit` must be an nm_fit.")
  outcome <- as.character(outcome)
  if (length(outcome) != 1L || !outcome %in% names(fit$data)) .nm_stop("Unknown categorical outcome column.")
  nsim <- as.integer(nsim)
  if (is.na(nsim) || nsim < 20L) .nm_stop("`nsim` must be at least 20.")
  if (!is.finite(level) || level <= 0 || level >= 1) .nm_stop("`level` must lie between zero and one.")
  data <- as.data.frame(fit$data)
  observed_rows <- data$EVID == 0L & data$MDV == 0L & !is.na(data[[outcome]])
  observed <- data[observed_rows, , drop = FALSE]
  categories <- sort(unique(observed[[outcome]]))
  declared <- NULL
  if (!is.null(fit$model$OUTCOMES)) {
    candidates <- Filter(function(value) {
      value$family %in% c("bernoulli", "categorical", "ordinal", "markov",
                          "continuous_time_markov")
    }, fit$model$OUTCOMES)
    if (length(candidates) == 1L) declared <- candidates[[1L]]
    if (length(candidates) > 1L && "DVID" %in% names(observed)) {
      dvid <- unique(observed$DVID)
      if (length(dvid) == 1L) {
        matching <- Filter(function(value) identical(as.numeric(value$dvid), as.numeric(dvid)), candidates)
        if (length(matching) == 1L) declared <- matching[[1L]]
      }
    }
  }
  if (!is.null(declared)) categories <- declared$categories %||% c(0, 1)
  if (length(categories) < 2L) .nm_stop("Categorical VPC requires at least two categories.")
  if (is.null(declared) && length(categories) != 2L) {
    .nm_stop("A multicategory VPC requires a first-class categorical OUTCOMES declaration.")
  }
  breaks <- .nm_vpc_breaks(observed$TIME, breaks)
  observed$BIN <- cut(observed$TIME, breaks, include.lowest = TRUE)
  observed_summary <- do.call(rbind, lapply(
    split(observed, observed$BIN, drop = TRUE), function(frame) do.call(rbind, lapply(
      categories, function(category) data.frame(
        BIN = as.character(frame$BIN[[1L]]), TIME = stats::median(frame$TIME),
        N = nrow(frame), CATEGORY = as.character(category),
        PROPORTION = mean(frame[[outcome]] == category), stringsAsFactors = FALSE
      )
    ))
  ))
  simulated <- nm_simulate(
    fit$model, fit$data, theta = fit$theta, sigma = fit$sigma, omega = fit$omega,
    nsim = nsim, random_effects = TRUE, residual = !is.null(declared),
    sample_mixture = TRUE, seed = seed
  )
  simulated <- simulated[rep(observed_rows, times = nsim), , drop = FALSE]
  if (is.null(declared)) {
    probability <- as.numeric(simulated$IPRED)
    if (mean(!is.finite(probability) | probability < -1e-6 | probability > 1 + 1e-6) > 0.01) {
      .nm_stop("Categorical VPC requires F/IPRED to represent a probability in [0, 1].")
    }
    probability <- pmin(pmax(probability, 0), 1)
    set.seed(as.integer(seed) + 1L)
    simulated$DV <- ifelse(stats::runif(nrow(simulated)) <= probability,
                           categories[[2L]], categories[[1L]])
  }
  simulated$CATEGORY <- simulated[[outcome]]
  simulated$BIN <- cut(simulated$TIME, breaks, include.lowest = TRUE)
  per_simulation <- do.call(rbind, lapply(
    split(simulated, list(simulated$SIM, simulated$BIN), drop = TRUE),
    function(frame) do.call(rbind, lapply(categories, function(category) data.frame(
      SIM = frame$SIM[[1L]], BIN = as.character(frame$BIN[[1L]]),
      TIME = stats::median(frame$TIME), CATEGORY = as.character(category),
      PROPORTION = mean(frame$CATEGORY == category), stringsAsFactors = FALSE
    )))
  ))
  alpha <- (1 - level) / 2
  intervals <- do.call(rbind, lapply(split(
    per_simulation, list(per_simulation$BIN, per_simulation$CATEGORY), drop = TRUE
  ), function(frame) {
    interval <- stats::quantile(frame$PROPORTION, c(alpha, 0.5, 1 - alpha),
                                names = FALSE, na.rm = TRUE)
    data.frame(
      BIN = frame$BIN[[1L]], TIME = stats::median(frame$TIME),
      CATEGORY = frame$CATEGORY[[1L]], lower = interval[[1L]],
      median = interval[[2L]], upper = interval[[3L]], stringsAsFactors = FALSE
    )
  }))
  structure(list(
    observed = observed_summary, simulated = intervals,
    per_simulation = per_simulation, categories = categories, outcome = outcome,
    breaks = breaks, level = level, nsim = nsim, seed = seed
  ), class = "nm_vpc_categorical")
}

#' Count visual predictive check
#'
#' Summarizes the mean, variance, zero fraction, median, and upper count
#' quantile in time bins for first-class Poisson, negative-binomial, binomial,
#' ZIP, and hurdle models.
#'
#' @param fit An `nm_fit`.
#' @param outcome Count column, normally `DV`.
#' @param dvid Optional endpoint `DVID` in a joint model.
#' @param nsim Number of simulations.
#' @param breaks Optional time bins.
#' @param level Simulation interval.
#' @param seed RNG seed.
#' @return An `nm_vpc_count`.
#' @export
nm_vpc_count <- function(fit, outcome = "DV", dvid = NULL, nsim = 200L,
                         breaks = NULL, level = 0.9, seed = 20260713L) {
  if (!inherits(fit, "nm_fit")) .nm_stop("`fit` must be an nm_fit.")
  if (is.null(fit$model$OUTCOMES)) .nm_stop("Count VPC requires first-class OUTCOMES.")
  if (!is.null(dvid) && (length(dvid) != 1L || !is.finite(dvid))) dvid <- NULL
  count_families <- c("poisson", "negative_binomial", "binomial",
                      "zero_inflated_poisson", "hurdle_poisson")
  candidates <- Filter(function(value) value$family %in% count_families,
                       fit$model$OUTCOMES)
  if (!is.null(dvid)) candidates <- Filter(function(value) {
    identical(as.numeric(value$dvid), as.numeric(dvid))
  }, candidates)
  if (length(candidates) != 1L) {
    .nm_stop("Select a unique count endpoint with `dvid`.")
  }
  endpoint <- candidates[[1L]]
  nsim <- as.integer(nsim)
  if (is.na(nsim) || nsim < 20L) .nm_stop("`nsim` must be at least 20.")
  if (!is.finite(level) || level <= 0 || level >= 1) .nm_stop("`level` must lie between zero and one.")
  data <- as.data.frame(fit$data)
  rows <- data$EVID == 0L & data$MDV == 0L & is.finite(data[[outcome]])
  if (!is.null(endpoint$dvid)) rows <- rows & data$DVID == endpoint$dvid
  observed <- data[rows, , drop = FALSE]
  if (!nrow(observed)) .nm_stop("No observed count records were found.")
  breaks <- .nm_vpc_breaks(observed$TIME, breaks)
  summarize <- function(frame) data.frame(
    BIN = as.character(frame$BIN[[1L]]), TIME = stats::median(frame$TIME),
    N = nrow(frame), MEAN = mean(frame[[outcome]]),
    VARIANCE = if (nrow(frame) > 1L) stats::var(frame[[outcome]]) else 0,
    ZERO = mean(frame[[outcome]] == 0),
    Q50 = unname(stats::quantile(frame[[outcome]], 0.5, type = 1)),
    Q90 = unname(stats::quantile(frame[[outcome]], 0.9, type = 1)),
    stringsAsFactors = FALSE
  )
  observed$BIN <- cut(observed$TIME, breaks, include.lowest = TRUE)
  observed_summary <- do.call(rbind, lapply(split(observed, observed$BIN, drop = TRUE), summarize))
  simulated <- nm_simulate(
    fit$model, fit$data, theta = fit$theta, sigma = fit$sigma, omega = fit$omega,
    nsim = nsim, random_effects = TRUE, residual = TRUE,
    sample_mixture = TRUE, seed = seed
  )
  simulated <- simulated[rep(rows, times = nsim), , drop = FALSE]
  simulated$BIN <- cut(simulated$TIME, breaks, include.lowest = TRUE)
  per_simulation <- do.call(rbind, lapply(split(
    simulated, list(simulated$SIM, simulated$BIN), drop = TRUE
  ), function(frame) cbind(SIM = frame$SIM[[1L]], summarize(frame))))
  measures <- c("MEAN", "VARIANCE", "ZERO", "Q50", "Q90")
  alpha <- (1 - level) / 2
  intervals <- do.call(rbind, lapply(split(per_simulation, per_simulation$BIN), function(frame) {
    row <- data.frame(BIN = frame$BIN[[1L]], TIME = stats::median(frame$TIME))
    for (measure in measures) {
      interval <- stats::quantile(frame[[measure]], c(alpha, 0.5, 1 - alpha),
                                  names = FALSE, na.rm = TRUE)
      row[[paste0(measure, "_lower")]] <- interval[[1L]]
      row[[paste0(measure, "_median")]] <- interval[[2L]]
      row[[paste0(measure, "_upper")]] <- interval[[3L]]
    }
    row
  }))
  structure(list(
    observed = observed_summary, simulated = intervals,
    per_simulation = per_simulation, outcome = outcome, dvid = endpoint$dvid,
    family = endpoint$family, breaks = breaks, level = level,
    nsim = nsim, seed = seed
  ), class = "nm_vpc_count")
}

.nm_km_at <- function(time, event, grid) {
  order <- order(time)
  time <- time[order]
  event <- event[order]
  unique_event <- sort(unique(time[event == 1L]))
  survival <- 1
  event_curve <- numeric(length(unique_event))
  for (index in seq_along(unique_event)) {
    current <- unique_event[[index]]
    at_risk <- sum(time >= current)
    events <- sum(time == current & event == 1L)
    if (at_risk > 0L) survival <- survival * (1 - events / at_risk)
    event_curve[[index]] <- survival
  }
  if (!length(unique_event)) return(rep(1, length(grid)))
  position <- findInterval(grid, unique_event)
  ifelse(position == 0L, 1, event_curve[pmax(position, 1L)])
}

#' Time-to-event visual predictive check
#'
#' `F`/`IPRED` is interpreted as a non-negative instantaneous hazard on each
#' subject's observation grid. The first event is simulated from the integrated
#' piecewise-constant hazard and Kaplan-Meier curves are summarized across
#' simulations.
#'
#' @param fit An `nm_fit`.
#' @param event Binary event-indicator column.
#' @param nsim Number of simulations.
#' @param level Simulation interval.
#' @param seed RNG seed.
#' @return An `nm_vpc_tte`.
#' @export
nm_vpc_tte <- function(fit, event = "DV", nsim = 200L, level = 0.9,
                       seed = 20260713L) {
  if (!inherits(fit, "nm_fit")) .nm_stop("`fit` must be an nm_fit.")
  event <- as.character(event)
  if (length(event) != 1L || !event %in% names(fit$data)) .nm_stop("Unknown event indicator column.")
  nsim <- as.integer(nsim)
  if (is.na(nsim) || nsim < 20L) .nm_stop("`nsim` must be at least 20.")
  if (!is.finite(level) || level <= 0 || level >= 1) .nm_stop("`level` must lie between zero and one.")
  data <- as.data.frame(fit$data)
  rows <- data$EVID == 0L & data$MDV == 0L & !is.na(data[[event]])
  observed <- data[rows, , drop = FALSE]
  if (!all(observed[[event]] %in% c(0, 1))) .nm_stop("The TTE event column must contain zero/one indicators.")
  subjects <- split(observed, observed$.ID_INDEX)
  subject_records <- do.call(rbind, lapply(subjects, function(frame) {
    event_rows <- which(frame[[event]] == 1L)
    data.frame(
      ID = frame$ID[[1L]],
      TIME = if (length(event_rows)) frame$TIME[[event_rows[[1L]]]] else max(frame$TIME),
      EVENT = as.integer(length(event_rows) > 0L)
    )
  }))
  grid <- sort(unique(observed$TIME))
  observed_curve <- data.frame(
    TIME = grid, SURVIVAL = .nm_km_at(subject_records$TIME, subject_records$EVENT, grid)
  )
  declared <- NULL
  if (!is.null(fit$model$OUTCOMES)) {
    candidates <- Filter(function(value) {
      value$family %in% c("tte", "competing_risks")
    }, fit$model$OUTCOMES)
    if (length(candidates) == 1L) declared <- candidates[[1L]]
  }
  simulated <- nm_simulate(
    fit$model, fit$data, theta = fit$theta, sigma = fit$sigma, omega = fit$omega,
    nsim = nsim, random_effects = TRUE, residual = !is.null(declared),
    sample_mixture = TRUE, seed = seed
  )
  simulated <- simulated[rep(rows, times = nsim), , drop = FALSE]
  if (is.null(declared) && any(!is.finite(simulated$IPRED) | simulated$IPRED < 0)) {
    .nm_stop("TTE VPC requires F/IPRED to be a finite non-negative hazard.")
  }
  set.seed(as.integer(seed) + 1L)
  curves <- matrix(1, nrow = length(grid), ncol = nsim)
  for (simulation in seq_len(nsim)) {
    sample <- simulated[simulated$SIM == simulation, , drop = FALSE]
    sample_subjects <- split(sample, sample$.ID_INDEX)
    records <- do.call(rbind, lapply(sample_subjects, function(frame) {
      frame <- frame[order(frame$TIME), , drop = FALSE]
      event_rows <- if (!is.null(declared)) {
        if (declared$family == "competing_risks") which(frame[[event]] != 0) else
          which(frame[[event]] == declared$event)
      } else {
        delta <- c(0, pmax(diff(frame$TIME), 0))
        probability <- 1 - exp(-pmax(frame$IPRED, 0) * delta)
        which(stats::runif(nrow(frame)) <= probability)
      }
      data.frame(
        TIME = if (length(event_rows)) frame$TIME[[event_rows[[1L]]]] else max(frame$TIME),
        EVENT = as.integer(length(event_rows) > 0L)
      )
    }))
    curves[, simulation] <- .nm_km_at(records$TIME, records$EVENT, grid)
  }
  alpha <- (1 - level) / 2
  intervals <- t(apply(curves, 1L, stats::quantile,
                       probs = c(alpha, 0.5, 1 - alpha), names = FALSE))
  simulated_curve <- data.frame(
    TIME = grid, lower = intervals[, 1L], median = intervals[, 2L],
    upper = intervals[, 3L]
  )
  structure(list(
    observed = observed_curve, simulated = simulated_curve,
    per_simulation = curves, event = event, level = level,
    nsim = nsim, seed = seed
  ), class = "nm_vpc_tte")
}

.nm_competing_curve <- function(records, causes, grid) {
  event_times <- sort(unique(records$TIME[records$CAUSE != 0]))
  survival <- 1
  cumulative <- stats::setNames(numeric(length(causes)), as.character(causes))
  history <- matrix(0, nrow = length(event_times), ncol = length(causes),
                    dimnames = list(NULL, as.character(causes)))
  for (index in seq_along(event_times)) {
    time <- event_times[[index]]
    at_risk <- sum(records$TIME >= time)
    if (at_risk > 0L) {
      events <- vapply(causes, function(cause) {
        sum(records$TIME == time & records$CAUSE == cause)
      }, numeric(1))
      cumulative <- cumulative + survival * events / at_risk
      survival <- survival * (1 - sum(events) / at_risk)
    }
    history[index, ] <- cumulative
  }
  if (!length(event_times)) history <- matrix(0, nrow = 1L, ncol = length(causes),
                                             dimnames = list(NULL, as.character(causes)))
  do.call(rbind, lapply(seq_along(causes), function(column) {
    position <- findInterval(grid, event_times)
    value <- if (!length(event_times)) rep(0, length(grid)) else
      ifelse(position == 0L, 0, history[pmax(position, 1L), column])
    data.frame(TIME = grid, CAUSE = as.character(causes[[column]]), CIF = value,
               stringsAsFactors = FALSE)
  }))
}

#' Competing-risk visual predictive check
#'
#' Uses Aalen-Johansen cumulative incidence curves for every declared cause and
#' compares them with simulation intervals from a first-class
#' `competing_risks` outcome.
#'
#' @param fit An `nm_fit`.
#' @param event Cause-code column (`0` denotes no event/censoring).
#' @param dvid Optional joint-endpoint `DVID`.
#' @param nsim Number of simulations.
#' @param level Simulation interval.
#' @param seed RNG seed.
#' @return An `nm_vpc_competing`.
#' @export
nm_vpc_competing <- function(fit, event = "DV", dvid = NULL, nsim = 200L,
                             level = 0.9, seed = 20260713L) {
  if (!inherits(fit, "nm_fit")) .nm_stop("`fit` must be an nm_fit.")
  if (!is.null(dvid) && (length(dvid) != 1L || !is.finite(dvid))) dvid <- NULL
  candidates <- Filter(function(value) value$family == "competing_risks",
                       fit$model$OUTCOMES %||% list())
  if (!is.null(dvid)) candidates <- Filter(function(value) {
    identical(as.numeric(value$dvid), as.numeric(dvid))
  }, candidates)
  if (length(candidates) != 1L) .nm_stop("Select a unique competing-risk endpoint with `dvid`.")
  endpoint <- candidates[[1L]]
  nsim <- as.integer(nsim)
  if (is.na(nsim) || nsim < 20L) .nm_stop("`nsim` must be at least 20.")
  data <- as.data.frame(fit$data)
  rows <- data$EVID == 0L & data$MDV == 0L & is.finite(data[[event]])
  if (!is.null(endpoint$dvid)) rows <- rows & data$DVID == endpoint$dvid
  observed <- data[rows, , drop = FALSE]
  records <- function(frame) do.call(rbind, lapply(split(frame, frame$.ID_INDEX), function(subject) {
    subject <- subject[order(subject$TIME), , drop = FALSE]
    event_row <- which(subject[[event]] != 0)[1L]
    if (is.na(event_row)) event_row <- integer()
    data.frame(
      TIME = if (length(event_row)) subject$TIME[[event_row]] else max(subject$TIME),
      CAUSE = if (length(event_row)) subject[[event]][[event_row]] else 0
    )
  }))
  grid <- sort(unique(observed$TIME))
  observed_curve <- .nm_competing_curve(records(observed), endpoint$categories, grid)
  simulated <- nm_simulate(
    fit$model, fit$data, theta = fit$theta, sigma = fit$sigma, omega = fit$omega,
    nsim = nsim, random_effects = TRUE, residual = TRUE,
    sample_mixture = TRUE, seed = seed
  )
  simulated <- simulated[rep(rows, times = nsim), , drop = FALSE]
  curves <- do.call(rbind, lapply(seq_len(nsim), function(index) {
    curve <- .nm_competing_curve(
      records(simulated[simulated$SIM == index, , drop = FALSE]),
      endpoint$categories, grid
    )
    curve$SIM <- index
    curve
  }))
  alpha <- (1 - level) / 2
  intervals <- do.call(rbind, lapply(split(
    curves, list(curves$CAUSE, curves$TIME), drop = TRUE
  ), function(frame) {
    interval <- stats::quantile(frame$CIF, c(alpha, 0.5, 1 - alpha), names = FALSE)
    data.frame(
      TIME = frame$TIME[[1L]], CAUSE = frame$CAUSE[[1L]],
      lower = interval[[1L]], median = interval[[2L]], upper = interval[[3L]],
      stringsAsFactors = FALSE
    )
  }))
  intervals <- intervals[order(as.numeric(intervals$CAUSE), intervals$TIME), , drop = FALSE]
  structure(list(
    observed = observed_curve, simulated = intervals, per_simulation = curves,
    event = event, dvid = endpoint$dvid, causes = endpoint$categories,
    level = level, nsim = nsim, seed = seed
  ), class = "nm_vpc_competing")
}

#' Recurrent-event visual predictive check
#'
#' Compares the observed mean cumulative event function with predictive
#' intervals from a first-class `recurrent_event` outcome.
#'
#' @param fit An `nm_fit`.
#' @param event Event-indicator column.
#' @param dvid Optional joint-endpoint `DVID`.
#' @param nsim Number of simulations.
#' @param level Simulation interval.
#' @param seed RNG seed.
#' @return An `nm_vpc_recurrent`.
#' @export
nm_vpc_recurrent <- function(fit, event = "DV", dvid = NULL, nsim = 200L,
                             level = 0.9, seed = 20260713L) {
  if (!inherits(fit, "nm_fit")) .nm_stop("`fit` must be an nm_fit.")
  if (!is.null(dvid) && (length(dvid) != 1L || !is.finite(dvid))) dvid <- NULL
  candidates <- Filter(function(value) value$family == "recurrent_event",
                       fit$model$OUTCOMES %||% list())
  if (!is.null(dvid)) candidates <- Filter(function(value) {
    identical(as.numeric(value$dvid), as.numeric(dvid))
  }, candidates)
  if (length(candidates) != 1L) .nm_stop("Select a unique recurrent-event endpoint with `dvid`.")
  endpoint <- candidates[[1L]]
  nsim <- as.integer(nsim)
  if (is.na(nsim) || nsim < 20L) .nm_stop("`nsim` must be at least 20.")
  data <- as.data.frame(fit$data)
  rows <- data$EVID == 0L & data$MDV == 0L & is.finite(data[[event]])
  if (!is.null(endpoint$dvid)) rows <- rows & data$DVID == endpoint$dvid
  observed <- data[rows, , drop = FALSE]
  grid <- sort(unique(observed$TIME))
  mean_cumulative <- function(frame) {
    subjects <- split(frame, frame$.ID_INDEX)
    curves <- vapply(subjects, function(subject) vapply(grid, function(time) {
      sum(subject[[event]][subject$TIME <= time] == endpoint$event)
    }, numeric(1)), numeric(length(grid)))
    if (is.null(dim(curves))) curves <- matrix(curves, ncol = 1L)
    rowMeans(curves)
  }
  observed_curve <- data.frame(TIME = grid, MEAN_CUMULATIVE = mean_cumulative(observed))
  simulated <- nm_simulate(
    fit$model, fit$data, theta = fit$theta, sigma = fit$sigma, omega = fit$omega,
    nsim = nsim, random_effects = TRUE, residual = TRUE,
    sample_mixture = TRUE, seed = seed
  )
  simulated <- simulated[rep(rows, times = nsim), , drop = FALSE]
  curves <- vapply(seq_len(nsim), function(index) {
    mean_cumulative(simulated[simulated$SIM == index, , drop = FALSE])
  }, numeric(length(grid)))
  alpha <- (1 - level) / 2
  intervals <- t(apply(curves, 1L, stats::quantile,
                       probs = c(alpha, 0.5, 1 - alpha), names = FALSE))
  simulated_curve <- data.frame(
    TIME = grid, lower = intervals[, 1L], median = intervals[, 2L], upper = intervals[, 3L]
  )
  structure(list(
    observed = observed_curve, simulated = simulated_curve,
    per_simulation = curves, event = event, dvid = endpoint$dvid,
    level = level, nsim = nsim, seed = seed
  ), class = "nm_vpc_recurrent")
}

#' Subject-level nonparametric bootstrap
#'
#' @param fit An `nm_fit` used as the model and estimation template.
#' @param n Number of bootstrap fits.
#' @param seed RNG seed.
#' @param level Percentile confidence level.
#' @param ... Controls passed to [nm_est()].
#' @return Bootstrap parameter matrix, convergence flags, and failed-run errors.
#' @export
nm_bootstrap <- function(fit, n = 100L, seed = 20260713L, level = 0.95, ...) {
  if (!inherits(fit, "nm_fit")) .nm_stop("`fit` must be an nm_fit.")
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) .nm_stop("`n` must be a positive integer.")
  if (!is.finite(level) || level <= 0 || level >= 1) .nm_stop("`level` must lie between zero and one.")
  set.seed(seed)
  ids <- unique(fit$data$ID)
  runs <- vector("list", n)
  errors <- character(n)
  for (iteration in seq_len(n)) {
    selected <- sample(ids, length(ids), replace = TRUE)
    pieces <- lapply(seq_along(selected), function(index) {
      block <- as.data.frame(fit$data[fit$data$ID == selected[[index]], , drop = FALSE])
      block$ID <- index
      internal <- grep("^\\.", names(block), value = TRUE)
      block[setdiff(internal, character())] <- NULL
      block
    })
    refit <- tryCatch(
      nm_est(fit$model, do.call(rbind, pieces), method = fit$method, ...),
      error = identity
    )
    if (inherits(refit, "error")) {
      errors[[iteration]] <- conditionMessage(refit)
    } else {
      runs[[iteration]] <- c(refit$theta, refit$sigma, refit$omega)
    }
  }
  successful <- Filter(Negate(is.null), runs)
  estimates <- if (length(successful)) do.call(rbind, successful) else matrix(numeric(), 0L, 0L)
  if (ncol(estimates)) {
    colnames(estimates) <- .nm_parameter_names(fit$theta, fit$sigma, fit$omega)
  }
  native <- .nm_fit_native_parameters(fit)
  alpha <- (1 - level) / 2
  summary <- if (nrow(estimates)) do.call(rbind, lapply(seq_along(native), function(index) {
    values <- estimates[, index]
    interval <- stats::quantile(values, c(alpha, 1 - alpha), na.rm = TRUE, names = FALSE)
    data.frame(
      parameter = names(native)[[index]], estimate = native[[index]],
      bootstrap_mean = mean(values, na.rm = TRUE), se = stats::sd(values, na.rm = TRUE),
      bias = mean(values, na.rm = TRUE) - native[[index]],
      lower = interval[[1L]], upper = interval[[2L]], stringsAsFactors = FALSE
    )
  })) else data.frame()
  structure(list(estimates = estimates, errors = errors[nzchar(errors)],
                 summary = summary, n = n, successful = nrow(estimates),
                 seed = seed, level = level),
            class = "nm_bootstrap")
}

#' @export
summary.nm_fit <- function(object, covariance = NULL, ...) {
  if (is.null(covariance) && !is.null(object$covariance)) covariance <- object$covariance
  parameter <- c(object$theta, object$sigma, object$omega)
  names(parameter) <- .nm_parameter_names(object$theta, object$sigma, object$omega)
  structure(list(
    method = object$method, objective = object$objective,
    convergence = object$convergence, parameters = parameter,
    covariance = covariance, posterior = object$posterior$population %||% NULL,
    eta = nm_etab(object)
  ), class = "summary.nm_fit")
}

#' @export
print.summary.nm_fit <- function(x, ...) {
  cat("LibeRation fit summary\n")
  cat("  method:", x$method, " objective:", format(x$objective),
      " convergence:", x$convergence, "\n\n")
  print(x$parameters)
  if (!is.null(x$covariance$se)) {
    cat("\nNative-scale standard errors\n")
    print(x$covariance$se)
  }
  if (!is.null(x$posterior$sd)) {
    cat("\nPosterior standard deviations\n")
    print(x$posterior$sd)
    cat("\nPosterior 95% credible intervals\n")
    print(x$posterior$quantile[c(1L, 3L), , drop = FALSE])
  }
  if (length(x$eta$shrinkage)) {
    cat("\nETA shrinkage\n")
    print(x$eta$shrinkage)
  }
  invisible(x)
}
