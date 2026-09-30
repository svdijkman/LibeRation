# Fit prediction, GOF, selected-output, ETAB, and summary methods.
# Split from diagnostics.R as a behaviour-preserving source move.

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


#' @export
summary.nm_fit <- function(object, covariance = NULL, ...) {
  if (is.null(covariance) && !is.null(object$covariance)) covariance <- object$covariance
  parameter <- c(object$theta, object$sigma, object$omega)
  names(parameter) <- .nm_parameter_names(object$theta, object$sigma, object$omega)
  gradient_description <- object$diagnostics$population_gradient %||%
    object$diagnostics$optimizer$population_gradient %||% "not reported"
  gradient_fallbacks <- object$diagnostics$optimizer$gradient_fallbacks %||% 0L
  gradient_class <- if (gradient_fallbacks > 0L) {
    "finite-difference fallback used"
  } else if (grepl(
    "omitted|finite common-random|finite adaptive-grid|derivative-free",
    gradient_description, ignore.case = TRUE
  )) {
    "score-incomplete or finite-grid derivative"
  } else if (grepl("exact|CppAD", gradient_description, ignore.case = TRUE)) {
    "CppAD-derived on the recorded smooth path"
  } else "estimator-specific; inspect description"
  structure(list(
    method = object$method, objective = object$objective,
    convergence = object$convergence, parameters = parameter,
    covariance = covariance, posterior = object$posterior$population %||% NULL,
    eta = nm_etab(object),
    derivative_provenance = list(
      gradient_class = gradient_class,
      gradient_description = gradient_description,
      gradient_fallbacks = as.integer(gradient_fallbacks),
      objective_backend = object$diagnostics$optimizer$objective_backend %||%
        object$diagnostics$optimizer$backend %||% "not reported",
      covariance_bread_source = covariance$bread_source %||% "not calculated",
      covariance_bread_exact = isTRUE(covariance$bread_exact)
    )
  ), class = "summary.nm_fit")
}

#' @export
print.summary.nm_fit <- function(x, ...) {
  cat("LibeRation fit summary\n")
  cat("  method:", x$method, " objective:", format(x$objective),
      " convergence:", x$convergence, "\n\n")
  cat("Derivative provenance\n")
  cat("  gradient:", x$derivative_provenance$gradient_class, "\n")
  cat("  detail:", x$derivative_provenance$gradient_description, "\n")
  cat("  objective backend:", x$derivative_provenance$objective_backend, "\n")
  if (!identical(x$derivative_provenance$covariance_bread_source, "not calculated")) {
    cat("  covariance bread:", x$derivative_provenance$covariance_bread_source,
        "(exact:", x$derivative_provenance$covariance_bread_exact, ")\n")
  }
  cat("\n")
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
