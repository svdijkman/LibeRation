.nm_engine_data <- function(model, data) {
  data <- if (inherits(data, "nm_dataset")) data else nm_dataset(data)
  iov <- model$LIK_CONFIG$iov
  if (iov > 0L) {
    column <- model$LIK_CONFIG$occasion_col
    if (!column %in% names(data)) {
      .nm_stop("IOV models require occasion column `", column, "` in the dataset.")
    }
    value <- as.character(data[[column]])
    if (anyNA(value) || any(!nzchar(value))) .nm_stop("Occasion identifiers must be non-missing.")
    if (".OCC_INDEX" %in% names(data) &&
        all(is.finite(data$.OCC_INDEX)) && all(data$.OCC_INDEX >= 1)) {
      data$.OCC_INDEX <- as.integer(data$.OCC_INDEX)
    } else {
      levels <- unique(value)
      data$.OCC_INDEX <- match(value, levels)
      attr(data, "occasion_levels") <- levels
    }
    total <- if (".OCC_TOTAL" %in% names(data) &&
                 all(is.finite(data$.OCC_TOTAL)) && all(data$.OCC_TOTAL >= 1)) {
      max(as.integer(data$.OCC_TOTAL))
    } else max(data$.OCC_INDEX)
    if (total < max(data$.OCC_INDEX)) {
      .nm_stop("Stored occasion layout is smaller than an observed occasion index.")
    }
    data$.OCC_TOTAL <- total
  }
  data
}

.nm_eta_columns <- function(model, data) {
  iov <- model$LIK_CONFIG$iov
  if (iov == 0L) return(model$n_eta)
  occasions <- if (".OCC_TOTAL" %in% names(data)) {
    max(data$.OCC_TOTAL)
  } else max(data$.OCC_INDEX)
  model$n_eta - iov + occasions * iov
}

.nm_effect_covariance <- function(model, data, values = model$OMEGAS$Value) {
  base <- .nm_omega_matrix(model, values)
  iov <- model$LIK_CONFIG$iov
  if (iov == 0L) return(base)
  between <- model$n_eta - iov
  occasions <- max(data$.OCC_INDEX)
  output <- matrix(0, between + occasions * iov, between + occasions * iov)
  if (between) output[seq_len(between), seq_len(between)] <- base[seq_len(between), seq_len(between)]
  for (occasion in seq_len(occasions)) {
    index <- between + (occasion - 1L) * iov + seq_len(iov)
    source <- between + seq_len(iov)
    output[index, index] <- base[source, source]
  }
  output
}

.nm_model_spec <- function(model) {
  list(
    version = model$version,
    advan = model$ADVAN,
    trans = model$TRANS,
    model_ss = model$SS,
    dose_cmp = model$DOSECMP,
    obs_cmp = model$OBSCMP,
    solver = model$SOLVER,
    error_type = model$ERROR_TYPE,
    lik_config = model$LIK_CONFIG,
    n_theta = nrow(model$THETAS),
    n_eta = model$n_eta,
    omega_row = model$OMEGAS$ROW,
    omega_col = model$OMEGAS$COL,
    pred_ir = model$pred_ir,
    output_names = intersect(
      model$OUTPUT %||% character(), model$pred_ir$output_names %||% character()
    ),
    des_ir = model$des_ir,
    n_state = model$n_state,
    ode_control = model$ODE_CONTROL,
    specialized_advan = isTRUE(getOption("LibeRation.specialized_advan", TRUE)),
    state_names = as.character(model$GRAPH$compartments$name %||% character())
    ,matrix_graph = .nm_matrix_graph_spec(model$GRAPH)
  )
}

#' Compiled LibeRation numerical engine
#'
#' The R6 object owns a disposable C++ pointer. Its `model` field remains fully
#' serializable and is the object sent to workers.
#'
#' @export
NMEngine <- R6::R6Class(
  "NMEngine",
  public = list(
    #' @field model Serializable validated [nm_model()] definition.
    model = NULL,
    #' @field pointer External pointer to the compiled C++ numerical engine.
    pointer = NULL,

    #' @description
    #' Compile a serializable model into a C++ engine.
    #' @param model A validated [nm_model()].
    #' @returns A new `NMEngine` object.
    initialize = function(model) {
      if (!inherits(model, "nm_model")) .nm_stop("`model` must be an nm_model.")
      self$model <- model
      self$pointer <- .liberation_engine_create(.nm_model_spec(model))
    },

    #' @description
    #' Run deterministic event-table prediction in C++.
    #' @param data NONMEM-style event data.
    #' @param theta Population fixed effects.
    #' @param eta Subject-by-effect random-effect matrix.
    #' @param sigma Residual parameters.
    #' @returns Event data augmented with `IPRED` and compartment amounts.
    simulate = function(data, theta = self$model$THETAS$Value,
                        eta = NULL, sigma = self$model$SIGMAS$Value) {
      data <- .nm_engine_data(self$model, data)
      n_subjects <- length(unique(data$.ID_INDEX))
      n_eta <- .nm_eta_columns(self$model, data)
      if (is.null(eta)) eta <- matrix(0, n_subjects, n_eta)
      eta <- as.matrix(eta)
      if (!identical(dim(eta), c(n_subjects, n_eta))) {
        .nm_stop("`eta` must have dimensions ", n_subjects, " x ", n_eta, ".")
      }
      raw <- .liberation_engine_simulate(
        self$pointer, data, as.numeric(theta), eta, as.numeric(sigma)
      )
      result <- as.data.frame(data)
      result$IPRED <- raw$ipred
      for (j in seq_len(ncol(raw$amounts))) {
        result[[paste0("A", j)]] <- raw$amounts[, j]
      }
      if (!is.null(raw$generated) && ncol(raw$generated)) {
        generated_names <- as.character(raw$output_names %||% colnames(raw$generated))
        for (j in seq_len(ncol(raw$generated))) {
          result[[generated_names[[j]]]] <- raw$generated[, j]
        }
      }
      attr(result, "solver") <- raw$solver
      attr(result, "state_names") <- raw$state_names
      result
    },

    #' @description
    #' Record the complete prediction calculation on a persistent AD tape.
    #' @param data NONMEM-style event data.
    #' @param theta Population fixed effects.
    #' @param eta Subject-by-effect random-effect matrix.
    #' @param sigma Residual parameters.
    #' @returns An internal `nm_prediction_tape` object.
    prediction_tape = function(data, theta = self$model$THETAS$Value,
                               eta = NULL, sigma = self$model$SIGMAS$Value) {
      data <- .nm_engine_data(self$model, data)
      n_subjects <- length(unique(data$.ID_INDEX))
      n_eta <- .nm_eta_columns(self$model, data)
      if (is.null(eta)) eta <- matrix(0, n_subjects, n_eta)
      eta <- as.matrix(eta)
      if (!identical(dim(eta), c(n_subjects, n_eta))) {
        .nm_stop("`eta` must have dimensions ", n_subjects, " x ", n_eta, ".")
      }
      theta <- as.numeric(theta)
      sigma <- as.numeric(sigma)
      pointer <- .liberation_prediction_tape_create(
        self$pointer, data, theta, eta, sigma
      )
      point <- c(theta, as.vector(t(eta)), sigma)
      structure(
        list(pointer = pointer, point = point, domain = attr(pointer, "domain"),
             data = data, n_subjects = n_subjects, n_eta = n_eta,
             dynamic_columns = attr(pointer, "dynamic_columns"),
             dynamic_parameters = attr(pointer, "dynamic_parameters"),
             propagation_kernel = attr(pointer, "propagation_kernel"),
             operation_count = attr(pointer, "operation_count"),
             variable_count = attr(pointer, "variable_count")),
        class = "nm_prediction_tape"
      )
    },

    #' @description
    #' Evaluate predictions and their exact automatic derivatives.
    #' @param data NONMEM-style event data.
    #' @param theta Population fixed effects.
    #' @param eta Subject-by-effect random-effect matrix.
    #' @param sigma Residual parameters.
    #' @param jacobian Whether to return the full Jacobian.
    #' @returns Prediction values, optional Jacobian, and AD domain names.
    prediction_derivatives = function(data, theta = self$model$THETAS$Value,
                                      eta = NULL, sigma = self$model$SIGMAS$Value,
                                      jacobian = TRUE) {
      tape <- self$prediction_tape(data, theta = theta, eta = eta, sigma = sigma)
      result <- .liberation_prediction_tape_eval(
        tape$pointer, tape$point, isTRUE(jacobian)
      )
      names(result$value) <- paste0("ROW_", seq_along(result$value))
      result$domain <- tape$domain
      result$propagation_kernel <- tape$propagation_kernel
      result$operation_count <- tape$operation_count
      result$variable_count <- tape$variable_count
      result$derivative_strategy <- attr(result, "derivative_strategy")
      result$jacobian_nonzeros <- attr(result, "jacobian_nonzeros")
      result
    },

    #' @description
    #' Record the joint population objective on a persistent AD tape.
    #' @param data NONMEM-style event data with observations.
    #' @param theta Population fixed effects.
    #' @param eta Subject-by-effect random-effect matrix.
    #' @param sigma Residual parameters.
    #' @param omega Random-effect covariance parameters.
    #' @param interaction Include ETA-residual interaction.
    #' @returns An internal `nm_objective_tape` object.
    objective_tape = function(data, theta = self$model$THETAS$Value,
                               eta = NULL, sigma = self$model$SIGMAS$Value,
                               omega = self$model$OMEGAS$Value,
                               interaction = TRUE) {
      data <- .nm_engine_data(self$model, data)
      n_subjects <- length(unique(data$.ID_INDEX))
      n_eta <- .nm_eta_columns(self$model, data)
      if (is.null(eta)) eta <- matrix(0, n_subjects, n_eta)
      eta <- as.matrix(eta)
      if (!identical(dim(eta), c(n_subjects, n_eta))) {
        .nm_stop("`eta` must have dimensions ", n_subjects, " x ", n_eta, ".")
      }
      theta <- as.numeric(theta)
      sigma <- as.numeric(sigma)
      omega <- as.numeric(omega)
      pointer <- .liberation_objective_tape_create(
        self$pointer, data, theta, eta, sigma, omega, isTRUE(interaction)
      )
      point <- c(theta, as.vector(t(eta)), sigma, omega)
      structure(
        list(pointer = pointer, point = point, domain = attr(pointer, "domain"),
             data = data, n_subjects = n_subjects, n_eta = n_eta),
        class = "nm_objective_tape"
      )
    },

    #' @description
    #' Evaluate the exact joint objective and requested derivatives.
    #' @param data NONMEM-style event data with observations.
    #' @param theta Population fixed effects.
    #' @param eta Subject-by-effect random-effect matrix.
    #' @param sigma Residual parameters.
    #' @param omega Random-effect covariance parameters.
    #' @param gradient Whether to return the exact gradient.
    #' @param hessian Whether to return the exact Hessian.
    #' @param interaction Include ETA-residual interaction.
    #' @returns Objective value, requested derivatives, and AD domain names.
    objective = function(data, theta = self$model$THETAS$Value,
                          eta = NULL, sigma = self$model$SIGMAS$Value,
                          omega = self$model$OMEGAS$Value,
                          gradient = TRUE, hessian = FALSE,
                          interaction = TRUE) {
      tape <- self$objective_tape(
        data, theta = theta, eta = eta, sigma = sigma, omega = omega,
        interaction = interaction
      )
      result <- .liberation_objective_tape_eval(
        tape$pointer, tape$point, isTRUE(gradient), isTRUE(hessian)
      )
      result$domain <- tape$domain
      result
    },

    #' @description
    #' Print a concise engine summary.
    #' @param ... Unused.
    #' @returns The engine, invisibly.
    print = function(...) {
      cat("Compiled LibeRation C++ engine\n")
      cat("  ADVAN", self$model$ADVAN, " solver:", self$model$SOLVER, "\n")
      cat("  pointer:", if (is.null(self$pointer)) "missing" else "ready", "\n")
      invisible(self)
    }
  )
)

#' Compile a serializable model into a pointer-backed engine
#' @param model An [nm_model()].
#' @export
nm_compile <- function(model) NMEngine$new(model)
