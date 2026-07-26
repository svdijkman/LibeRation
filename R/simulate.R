#' Simulate or predict a LibeRation model
#'
#' @param model An `nm_model` or compiled `NMEngine`.
#' @param data NONMEM-style event data.
#' @param theta Population fixed effects; defaults to the model table.
#' @param eta Subject ETA matrix ordered by first subject appearance.
#' @param sigma Residual parameters.
#' @param omega Random-effect covariance parameters.
#' @param nsim Number of simulation replicates.
#' @param random_effects Sample subject/occasion ETAs from OMEGA when `eta` is
#'   not supplied.
#' @param residual Sample observations from the configured residual model.
#' @param sample_mixture Sample finite-mixture membership per subject.
#' @param censor Mark simulated values below LLOQ as BLQ/CENS and replace their
#'   reported DV by the applicable quantification limit.
#' @param seed Optional reproducible RNG seed.
#' @param n_cores Number of parallel simulation workers. Replicates are
#'   distributed across persistent PSOCK workers on Windows, Linux, and macOS.
#' @return Input records augmented with `IPRED`, compartment amounts, and—when
#'   requested—simulated `DV`, ETAs, mixture membership, and replicate number.
#' @examples
#' \donttest{
#' model <- nm_model(
#'   INPUT = c("ID", "TIME", "EVID", "AMT", "CMT", "DV", "MDV"),
#'   ADVAN = 1,
#'   PRED = "CL=THETA(1); V=THETA(2); S1=V",
#'   ERROR = "Y=F",
#'   THETAS = data.frame(THETA = 1:2, Value = c(2, 20))
#' )
#' events <- data.frame(
#'   ID = 1, TIME = c(0, 1, 2), EVID = c(1, 0, 0),
#'   AMT = c(100, 0, 0), CMT = 1, DV = NA_real_, MDV = c(1, 0, 0)
#' )
#' nm_simulate(model, events)
#' }
#' @export
nm_simulate <- function(model, data, theta = NULL, eta = NULL, sigma = NULL,
                        omega = NULL, nsim = 1L, random_effects = FALSE,
                        residual = FALSE, sample_mixture = FALSE,
                        censor = FALSE, seed = NULL, n_cores = 1L) {
  engine <- if (inherits(model, "NMEngine")) model else nm_compile(model)
  if (isTRUE(residual) && identical(engine$model$LIK_CONFIG$error, "likelihood")) {
    if (is.null(engine$model$OUTCOMES) && is.null(engine$model$KALMAN_CONFIG)) {
      .nm_stop(
        "A free-form user likelihood does not define a generic outcome generator. ",
        "Declare OUTCOMES with `nm_outcome()` or use `residual = FALSE`."
      )
    }
  }
  theta <- theta %||% engine$model$THETAS$Value
  sigma <- sigma %||% engine$model$SIGMAS$Value
  omega <- omega %||% engine$model$OMEGAS$Value
  nsim <- as.integer(nsim)
  if (length(nsim) != 1L || is.na(nsim) || nsim < 1L) {
    .nm_stop("`nsim` must be a positive integer.")
  }
  n_cores <- as.integer(n_cores)
  if (length(n_cores) != 1L || is.na(n_cores) || n_cores < 1L) {
    .nm_stop("`n_cores` must be a positive integer.")
  }
  n_cores <- min(n_cores, nsim)
  if (!is.null(seed)) {
    seed <- as.integer(seed)
    if (length(seed) != 1L || is.na(seed)) .nm_stop("`seed` must be one integer.")
    set.seed(seed)
  }
  normalized <- .nm_engine_data(engine$model, data)
  if (isTRUE(random_effects) && !is.null(eta)) {
    .nm_stop("Supply `eta` or set `random_effects = TRUE`, not both.")
  }
  if (n_cores > 1L) {
    replicate_seeds <- sample.int(.Machine$integer.max, nsim)
    cluster <- parallel::makePSOCKcluster(n_cores, outfile = "")
    on.exit(try(parallel::stopCluster(cluster), silent = TRUE), add = TRUE)
    parallel::clusterCall(cluster, function(paths) {
      .libPaths(unique(c(paths, .libPaths())))
      if (!requireNamespace("LibeRation", quietly = TRUE)) {
        stop("LibeRation is not installed in the parallel worker library paths.")
      }
      TRUE
    }, .libPaths())
    pieces <- parallel::parLapply(
      cluster, seq_len(nsim),
      function(index, specification, dataset, theta, eta, sigma, omega,
               random_effects, residual, sample_mixture, censor, seeds) {
        result <- LibeRation::nm_simulate(
          specification, dataset, theta = theta, eta = eta, sigma = sigma,
          omega = omega, nsim = 1L, random_effects = random_effects,
          residual = residual, sample_mixture = sample_mixture,
          censor = censor, seed = seeds[[index]], n_cores = 1L
        )
        result$SIM <- index
        result
      },
      specification = engine$model, dataset = normalized, theta = theta,
      eta = eta, sigma = sigma, omega = omega,
      random_effects = random_effects, residual = residual,
      sample_mixture = sample_mixture, censor = censor,
      seeds = replicate_seeds
    )
    solver <- attr(pieces[[1L]], "solver")
    state_names <- attr(pieces[[1L]], "state_names")
    output <- do.call(rbind, pieces)
    rownames(output) <- NULL
    attr(output, "solver") <- solver
    attr(output, "state_names") <- state_names
    attr(output, "parallel_cores") <- n_cores
    return(output)
  }
  simulate_one <- function(index) {
    replicate_data <- normalized
    n_subjects <- length(unique(replicate_data$.ID_INDEX))
    n_eta <- .nm_eta_columns(engine$model, replicate_data)
    eta_value <- eta
    if (isTRUE(random_effects)) {
      covariance <- .nm_effect_covariance(engine$model, replicate_data, omega)
      if (n_eta) {
        root <- t(chol(covariance))
        eta_value <- matrix(stats::rnorm(n_subjects * n_eta), n_subjects, n_eta) %*%
          t(root)
      } else eta_value <- matrix(numeric(), n_subjects, 0L)
    }
    if (is.null(eta_value)) eta_value <- matrix(0, n_subjects, n_eta)
    eta_value <- as.matrix(eta_value)
    mixture <- engine$model$LIK_CONFIG$mixtures
    if (!is.null(mixture) && isTRUE(sample_mixture)) {
      assignment <- sample.int(
        length(mixture$probability), n_subjects, replace = TRUE,
        prob = mixture$probability
      )
      replicate_data$MIXNUM <- assignment[replicate_data$.ID_INDEX]
    }
    result <- engine$simulate(
      replicate_data, theta = theta, eta = eta_value, sigma = sigma
    )
    if ("PRED" %in% (engine$model$OUTPUT %||% character())) {
      if (!length(eta_value) || all(eta_value == 0)) {
        result$PRED <- result$IPRED
      } else {
        population_eta <- matrix(0, n_subjects, n_eta)
        result$PRED <- engine$simulate(
          replicate_data, theta = theta, eta = population_eta, sigma = sigma
        )$IPRED
      }
    }
    if (n_eta) {
      if (!is.null(engine$model$RE_CONFIG)) {
        for (effect in seq_len(engine$model$n_eta)) {
          mapped <- as.integer(result[[paste0(".ETA_COLUMN_", effect)]])
          result[[paste0("ETA", effect)]] <- eta_value[
            cbind(result$.ID_INDEX, mapped)
          ]
        }
      } else {
        for (column in seq_len(n_eta)) {
          result[[paste0("ETA", column)]] <- eta_value[result$.ID_INDEX, column]
        }
      }
    }
    if (!is.null(mixture)) {
      result$MIXNUM <- result$MIXNUM %||% rep(1L, nrow(result))
    }
    if (isTRUE(residual)) {
      result <- if (!is.null(engine$model$KALMAN_CONFIG)) {
        .nm_simulate_kalman(
          engine, result, theta = theta, eta = eta_value, sigma = sigma
        )
      } else if (identical(engine$model$LIK_CONFIG$error, "likelihood")) {
        .nm_simulate_first_class_outcomes(engine$model, result, theta, sigma)
      } else {
        .nm_simulate_residual(engine$model, result, sigma, censor)
      }
    }
    if (nsim > 1L) result$SIM <- index
    result
  }
  pieces <- lapply(seq_len(nsim), simulate_one)
  solver <- attr(pieces[[1L]], "solver")
  state_names <- attr(pieces[[1L]], "state_names")
  output <- do.call(rbind, pieces)
  rownames(output) <- NULL
  attr(output, "solver") <- solver
  attr(output, "state_names") <- state_names
  output
}

.nm_simulate_residual <- function(model, result, sigma, censor) {
  observed <- result$EVID == 0L & result$MDV == 0L & is.finite(result$IPRED)
  result$DV <- as.numeric(result$DV)
  if (!any(observed)) return(result)
  dvid <- if ("DVID" %in% names(result)) result$DVID else rep(1L, nrow(result))
  variance <- .nm_residual_variance(model, result$IPRED, sigma, dvid)
  standardized <- numeric(nrow(result))
  groups <- interaction(result$.ID_INDEX, dvid, drop = TRUE, lex.order = TRUE)
  rho <- .nm_ar1_rho(model, sigma = sigma)
  for (group in levels(groups)) {
    rows <- which(observed & groups == group)
    if (!length(rows)) next
    innovation <- stats::rnorm(length(rows))
    standardized[rows[[1L]]] <- innovation[[1L]]
    if (model$LIK_CONFIG$sigma_corr == "ar1" && length(rows) > 1L) {
      for (position in seq.int(2L, length(rows))) {
        standardized[rows[[position]]] <- rho * standardized[rows[[position - 1L]]] +
          sqrt(1 - rho^2) * innovation[[position]]
      }
    } else if (length(rows) > 1L) {
      standardized[rows[-1L]] <- innovation[-1L]
    }
  }
  if (length(model$LIK_CONFIG$residual_groups)) {
    time_groups <- interaction(
      result$.ID_INDEX, result$TIME, drop = TRUE, lex.order = TRUE
    )
    for (definition in model$LIK_CONFIG$residual_groups) {
      correlation <- .nm_residual_group_value(definition, model$THETAS$Value, sigma)
      eligible <- observed & dvid %in% definition$dvid
      for (time_group in levels(time_groups[eligible])) {
        rows <- which(eligible & time_groups == time_group)
        if (length(rows) < 2L) next
        endpoint <- match(dvid[rows], definition$dvid)
        if (anyDuplicated(endpoint)) {
          .nm_stop("A correlated residual group has duplicate DVID observations at one subject/time.")
        }
        submatrix <- correlation[endpoint, endpoint, drop = FALSE]
        standardized[rows] <- drop(stats::rnorm(length(rows)) %*% chol(submatrix))
      }
    }
  }
  error <- sqrt(variance) * standardized
  if (model$LIK_CONFIG$error == "exponential") {
    result$DV[observed] <- result$IPRED[observed] * exp(error[observed])
  } else {
    result$DV[observed] <- result$IPRED[observed] + error[observed]
  }
  if (isTRUE(censor) && model$LIK_CONFIG$blq_method != "none") {
    limit <- if ("LLOQ" %in% names(result)) {
      as.numeric(result$LLOQ)
    } else rep(model$LIK_CONFIG$lloq, nrow(result))
    censored <- observed & is.finite(limit) & result$DV < limit
    result$BLQ <- as.integer(censored)
    result$CENS <- as.integer(censored)
    result$DV[censored] <- limit[censored]
  }
  result
}

#' Exact derivatives of model predictions
#'
#' Records the complete analytical event/ADVAN/matrix path with CppAD and
#' returns derivatives with respect to THETA, every subject ETA, and SIGMA.
#' Dataset values and event times are treated as fixed inputs. General ODE,
#' Michaelis--Menten, and equilibrium-DAE ADVANs record their accepted
#' adaptive trajectories--including periodic steady-state shooting--on the
#' same persistent CppAD tape.
#'
#' @param model An `nm_model` or compiled `NMEngine`.
#' @param data NONMEM-style event data.
#' @param theta Population fixed effects.
#' @param eta Subject ETA matrix.
#' @param sigma Residual parameters.
#' @param jacobian Whether to return the full prediction Jacobian.
#' @return A list with `value`, `jacobian`, and differentiation `domain`.
#' @export
nm_prediction_derivatives <- function(model, data, theta = NULL, eta = NULL,
                                      sigma = NULL, jacobian = TRUE) {
  engine <- if (inherits(model, "NMEngine")) model else nm_compile(model)
  theta <- theta %||% engine$model$THETAS$Value
  sigma <- sigma %||% engine$model$SIGMAS$Value
  engine$prediction_derivatives(
    data, theta = theta, eta = eta, sigma = sigma, jacobian = jacobian
  )
}
