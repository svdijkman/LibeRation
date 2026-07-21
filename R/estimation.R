.nm_subject_data <- function(data, subject) {
  out <- as.data.frame(data[data$.ID_INDEX == subject, , drop = FALSE])
  internal <- intersect(c(".ID_INDEX", ".source_row", ".generated", ".sort_priority"), names(out))
  out[internal] <- NULL
  nm_dataset(out)
}

.nm_residual_variance <- function(model, prediction, sigma, dvid = 1L) {
  type <- model$LIK_CONFIG$error
  if (length(dvid) == 1L && length(prediction) > 1L) {
    dvid <- rep.int(dvid, length(prediction))
  }
  per_response <- if (type %in% c("combined", "power")) 2L else 1L
  offset <- (pmax(as.integer(dvid), 1L) - 1L) * per_response
  offset[offset + per_response > length(sigma)] <- 0L
  s1 <- sigma[offset + 1L]
  square <- function(value) {
    if (identical(model$LIK_CONFIG$sigma_parameterization, "variance")) value else value^2
  }
  variance <- switch(
    type,
    additive = square(s1),
    proportional = square(s1) * prediction^2,
    exponential = square(s1),
    power = square(s1) * pmax(abs(prediction), 1e-12)^(2 * sigma[offset + 2L]),
    combined = square(s1) * prediction^2 + square(sigma[offset + 2L]),
    .nm_stop("A residual likelihood is required for estimation.")
  )
  pmax(variance, 1e-16)
}

.nm_positive_definite <- function(matrix, context = "curvature") {
  matrix <- (matrix + t(matrix)) / 2
  if (!length(matrix)) return(list(matrix = matrix, logdet = 0, jitter = 0))
  eigenvalues <- eigen(matrix, symmetric = TRUE, only.values = TRUE)$values
  largest <- max(abs(eigenvalues), 1)
  jitter <- max(0, largest * 1e-9 - min(eigenvalues))
  if (jitter > largest * 1e-2) {
    .nm_stop(context, " is not sufficiently positive definite.")
  }
  adjusted <- matrix + diag(jitter, nrow(matrix))
  determinant <- determinant(adjusted, logarithm = TRUE)
  if (determinant$sign <= 0 || !is.finite(determinant$modulus)) {
    .nm_stop(context, " determinant is not positive and finite.")
  }
  list(matrix = adjusted, logdet = as.numeric(determinant$modulus), jitter = jitter)
}

.nm_prediction_dynamic_columns <- function(model, data) {
  inputs <- unique(c(
    model$pred_ir$input_names %||% character(),
    model$des_ir$input_names %||% character()
  ))
  structural <- c(
    "ID", "TIME", "AMT", "RATE", "II", "ADDL", "EVID", "CMT",
    "SS", "MIXNUM", "DVID", "DV", "MDV", "LLOQ", "BLQ", "CENS",
    ".ID_INDEX", ".OCC_INDEX", "F", "T"
  )
  parameter <- grepl("^(THETA|ETA|SIGMA|ERR|A)_", inputs)
  intersect(inputs[!parameter & !inputs %in% structural], names(data))
}

.nm_structure_key <- function(value) {
  raw <- serialize(value, NULL, version = 2L)
  bytes <- as.integer(raw)
  weight <- (seq_along(bytes) %% 65521L) + 1L
  paste(
    length(bytes),
    format(sum(bytes) %% 2147483647, scientific = FALSE, trim = TRUE),
    format(sum((bytes * weight) %% 2147483629) %% 2147483647,
           scientific = FALSE, trim = TRUE),
    sep = "-"
  )
}

.nm_prediction_structure <- function(model, data) {
  dynamic <- .nm_prediction_dynamic_columns(model, data)
  ignored <- c(
    "DV", "MDV", "LLOQ", "BLQ", "CENS", ".source_row", ".generated",
    ".sort_priority", ".ID_INDEX", dynamic
  )
  source <- paste(model$PRED %||% "", model$DES %||% "")
  if (!grepl("\\bID\\b", source, perl = TRUE)) ignored <- c(ignored, "ID")
  columns <- setdiff(names(data), ignored)
  value <- list(
    dynamic_columns = dynamic,
    structural_data = as.data.frame(data[columns], stringsAsFactors = FALSE),
    rows = nrow(data)
  )
  list(key = .nm_structure_key(value), value = value)
}

.nm_prediction_pool_tape <- function(pool, engine, data, theta, sigma, n_eta) {
  structure <- .nm_prediction_structure(engine$model, data)
  bucket <- pool[[structure$key]] %||% list()
  if (length(bucket)) {
    for (entry in bucket) {
      if (identical(entry$structure, structure$value)) return(entry$tape)
    }
  }
  tape <- engine$prediction_tape(
    data, theta = theta, eta = matrix(0, 1L, n_eta), sigma = sigma
  )
  pool[[structure$key]] <- c(
    bucket, list(list(structure = structure$value, tape = tape))
  )
  tape
}

.nm_fo_structure <- function(model, data) {
  prediction <- .nm_prediction_structure(model, data)
  observed <- data$EVID == 0L & data$MDV == 0L & is.finite(data$DV)
  value <- c(prediction$value, list(
    observed = as.logical(observed),
    dvid = if ("DVID" %in% names(data)) as.integer(data$DVID) else rep.int(1L, nrow(data))
  ))
  list(key = .nm_structure_key(value), value = value)
}

.nm_fo_pool_tape <- function(pool, evaluator, theta, sigma, omega) {
  structure <- .nm_fo_structure(evaluator$engine$model, evaluator$data)
  bucket <- pool[[structure$key]] %||% list()
  if (length(bucket)) {
    for (entry in bucket) {
      if (identical(entry$structure, structure$value)) {
        evaluator$fo_tape <- entry$tape
        return(invisible(entry$tape))
      }
    }
  }
  evaluator$ensure_fo_tape(theta, sigma, omega)
  pool[[structure$key]] <- c(
    bucket, list(list(structure = structure$value, tape = evaluator$fo_tape))
  )
  invisible(evaluator$fo_tape)
}

.NMSubjectEvaluator <- R6::R6Class(
  ".NMSubjectEvaluator",
  public = list(
    engine = NULL,
    data = NULL,
    objective_tape = NULL,
    noninteraction_tape = NULL,
    prediction_tape = NULL,
    fo_tape = NULL,
    curvature_tapes = NULL,
    tape_anchor = NULL,
    tape_signature = NULL,
    tape_records = 0L,
    tape_retapes = 0L,
    tape_checks = 0L,
    tape_profile = "full",
    n_theta = 0L,
    n_eta = 0L,
    n_sigma = 0L,

    initialize = function(engine, data, theta, sigma, omega, n_eta = NULL,
                          prediction_tape = NULL,
                          tape_profile = c("full", "fo")) {
      self$tape_profile <- match.arg(tape_profile)
      self$engine <- engine
      self$data <- .nm_engine_data(engine$model, data)
      self$n_theta <- length(theta)
      self$n_eta <- as.integer(n_eta %||% .nm_eta_columns(engine$model, self$data))
      self$n_sigma <- length(sigma)
      self$tape_signature <- .nm_prediction_structure(engine$model, self$data)$key
      self$record_tapes(
        theta, sigma, omega, rep(0, self$n_eta),
        prediction_tape = prediction_tape
      )
    },

    record_tapes = function(theta, sigma, omega, eta = rep(0, self$n_eta),
                            retape = FALSE, prediction_tape = NULL) {
      self$curvature_tapes <- list()
      self$fo_tape <- NULL
      eta <- matrix(as.numeric(eta), 1L, self$n_eta)
      self$prediction_tape <- prediction_tape %||% self$engine$prediction_tape(
        self$data, theta = theta, eta = eta, sigma = sigma
      )
      .liberation_prediction_tape_new_dynamic(
        self$prediction_tape$pointer, self$data
      )
      self$objective_tape <- self$engine$objective_tape(
        self$data, theta = theta, eta = eta, sigma = sigma, omega = omega
      )
      self$noninteraction_tape <- if (identical(self$tape_profile, "fo")) NULL else
        self$engine$objective_tape(
          self$data, theta = theta, eta = eta, sigma = sigma, omega = omega,
          interaction = FALSE
        )
      self$tape_anchor <- c(theta, as.numeric(eta), sigma, omega)
      self$tape_records <- self$tape_records + 1L
      if (isTRUE(retape)) self$tape_retapes <- self$tape_retapes + 1L
      invisible(self)
    },

    ensure_valid_tapes = function(theta, sigma, omega,
                                  eta = rep(0, self$n_eta)) {
      self$tape_checks <- self$tape_checks + 1L
      # CppAD conditional expressions remain valid without retaping. Adaptive
      # ODE solvers, however, record one accepted-step trajectory, so a
      # materially different population/ETA point is deliberately retaped.
      if (!isTRUE(self$engine$model$USE_ODE)) return(FALSE)
      point <- c(theta, as.numeric(eta), sigma, omega)
      anchor <- self$tape_anchor
      radius <- getOption("LibeRation.tape_guard_radius", 0.5)
      distance <- max(abs(point - anchor) / pmax(abs(anchor), 1), na.rm = TRUE)
      if (is.finite(distance) && distance > radius) {
        self$record_tapes(theta, sigma, omega, eta, retape = TRUE)
        return(TRUE)
      }
      FALSE
    },

    tape_telemetry = function() list(
      signature = self$tape_signature, records = self$tape_records,
      retapes = self$tape_retapes, validity_checks = self$tape_checks
    ),

    objective_point = function(theta, eta, sigma, omega) {
      c(theta, eta, sigma, omega)
    },

    prediction_point = function(theta, eta, sigma) c(theta, eta, sigma),

    fo_point = function(theta, sigma, omega) c(theta, sigma, omega),

    ensure_fo_tape = function(theta, sigma, omega) {
      if (isTRUE(self$engine$model$USE_ODE)) {
        self$ensure_valid_tapes(theta, sigma, omega)
      }
      if (is.null(self$fo_tape)) {
        .liberation_prediction_tape_new_dynamic(
          self$prediction_tape$pointer, self$data
        )
        self$fo_tape <- list(
          pointer = .liberation_fo_tape_create(
            self$engine$pointer, self$prediction_tape$pointer, self$data,
            as.numeric(theta), as.numeric(sigma), as.numeric(omega)
          )
        )
      }
      .liberation_fo_tape_new_dynamic(self$fo_tape$pointer, self$data)
      invisible(self$fo_tape)
    },

    fo_objective = function(theta, sigma, omega,
                            gradient = FALSE, hessian = FALSE) {
      private$with_retape(
        function() {
          self$ensure_fo_tape(theta, sigma, omega)
          .liberation_objective_tape_eval(
            self$fo_tape$pointer, self$fo_point(theta, sigma, omega),
            gradient, hessian
          )
        }, theta, sigma, omega, rep(0, self$n_eta)
      )
    },

    ensure_curvature_tape = function(theta, eta, sigma, omega, approximation) {
      if (isTRUE(self$engine$model$USE_ODE)) {
        self$ensure_valid_tapes(theta, sigma, omega, eta)
      }
      approximation <- match.arg(approximation, c("foce", "focei", "laplace"))
      if (is.null(self$curvature_tapes[[approximation]])) {
        # Validate the primary tape at the curvature anchor first. A changed
        # pivot or adaptive trajectory is retaped by objective() before the
        # nested base2ad curvature tape is recorded.
        self$objective(
          theta, eta, sigma, omega, gradient = FALSE,
          interaction = TRUE
        )
        .liberation_prediction_tape_new_dynamic(
          self$prediction_tape$pointer, self$data
        )
        self$curvature_tapes[[approximation]] <- list(
          pointer = .liberation_curvature_tape_create(
            self$engine$pointer, self$prediction_tape$pointer,
            self$objective_tape$pointer, self$data, as.numeric(theta),
            as.numeric(eta), as.numeric(sigma), as.numeric(omega),
            approximation
          )
        )
      }
      invisible(self$curvature_tapes[[approximation]])
    },

    curvature = function(theta, eta, sigma, omega, approximation,
                         gradient = TRUE) {
      private$with_retape(
        function() {
          self$ensure_curvature_tape(theta, eta, sigma, omega, approximation)
          .liberation_objective_tape_eval(
            self$curvature_tapes[[approximation]]$pointer,
            self$objective_point(theta, eta, sigma, omega),
            isTRUE(gradient), FALSE
          )
        }, theta, sigma, omega, eta
      )
    },

    objective = function(theta, eta, sigma, omega,
                         gradient = FALSE, hessian = FALSE,
                         interaction = TRUE) {
      if (isTRUE(self$engine$model$USE_ODE)) {
        self$ensure_valid_tapes(theta, sigma, omega, eta)
      }
      private$with_retape(
        function() {
          .liberation_objective_tape_eval(
            if (isTRUE(interaction)) self$objective_tape$pointer else self$noninteraction_tape$pointer,
            self$objective_point(theta, eta, sigma, omega),
            gradient, hessian
          )
        }, theta, sigma, omega, eta
      )
    },

    objective_eta_values = function(theta, eta, sigma, omega,
                                    interaction = TRUE) {
      eta <- as.matrix(eta)
      if (ncol(eta) != self$n_eta) {
        .nm_stop("ETA samples have the wrong number of columns.")
      }
      if (isTRUE(self$engine$model$USE_ODE)) {
        self$ensure_valid_tapes(theta, sigma, omega, eta[1L, ])
      }
      private$with_retape(
        function() {
          tape <- if (isTRUE(interaction)) self$objective_tape else self$noninteraction_tape
          .liberation_objective_tape_eta_values(
            tape$pointer,
            self$objective_point(theta, rep(0, self$n_eta), sigma, omega),
            self$n_theta + seq_len(self$n_eta), eta
          )
        }, theta, sigma, omega, eta[1L, ]
      )
    },

    objective_eta_batch = function(theta, eta, sigma, omega,
                                   interaction = TRUE) {
      eta <- as.matrix(eta)
      if (ncol(eta) != self$n_eta) {
        .nm_stop("ETA samples have the wrong number of columns.")
      }
      if (isTRUE(self$engine$model$USE_ODE)) {
        self$ensure_valid_tapes(theta, sigma, omega, eta[1L, ])
      }
      points <- cbind(
        matrix(theta, nrow(eta), length(theta), byrow = TRUE), eta,
        matrix(sigma, nrow(eta), length(sigma), byrow = TRUE),
        matrix(omega, nrow(eta), length(omega), byrow = TRUE)
      )
      private$with_retape(
        function() {
          tape <- if (isTRUE(interaction)) self$objective_tape else self$noninteraction_tape
          .liberation_objective_tape_point_gradients(tape$pointer, points)
        }, theta, sigma, omega, eta[1L, ]
      )
    },

    objective_hessian_subset = function(theta, eta, sigma, omega,
                                        rows, columns,
                                        interaction = TRUE) {
      if (isTRUE(self$engine$model$USE_ODE)) {
        self$ensure_valid_tapes(theta, sigma, omega, eta)
      }
      private$with_retape(
        function() {
          tape <- if (isTRUE(interaction)) self$objective_tape else self$noninteraction_tape
          .liberation_objective_tape_hessian_subset(
            tape$pointer, self$objective_point(theta, eta, sigma, omega),
            as.integer(rows), as.integer(columns)
          )
        }, theta, sigma, omega, eta
      )
    },

    prediction = function(theta, eta, sigma, jacobian = FALSE, columns = NULL) {
      omega_offset <- length(theta) + self$n_eta + length(sigma)
      omega <- self$tape_anchor[omega_offset + seq_len(
        length(self$tape_anchor) - omega_offset
      )]
      if (isTRUE(self$engine$model$USE_ODE)) {
        self$ensure_valid_tapes(theta, sigma, omega, eta)
      }
      private$with_retape(
        function() {
          .liberation_prediction_tape_new_dynamic(
            self$prediction_tape$pointer, self$data
          )
          point <- self$prediction_point(theta, eta, sigma)
          if (isTRUE(jacobian) && !is.null(columns)) {
            return(.liberation_prediction_tape_eval_subset(
              self$prediction_tape$pointer, point, as.integer(columns)
            ))
          }
          .liberation_prediction_tape_eval(
            self$prediction_tape$pointer, point, jacobian
          )
        }, theta, sigma, omega, eta
      )
    },

    eta_mode = function(theta, sigma, omega, start = rep(0, self$n_eta),
                         maxit = 100L, tolerance = 1e-7,
                         interaction = TRUE, exact_hessian = TRUE) {
      if (isTRUE(self$engine$model$USE_ODE)) {
        self$ensure_valid_tapes(theta, sigma, omega, start)
      }
      if (self$n_eta == 0L) {
        value <- self$objective(theta, numeric(), sigma, omega)$value
        return(list(par = numeric(), value = value, convergence = 0L,
                    hessian = matrix(numeric(), 0, 0), jitter = 0))
      }
      eta_positions <- self$n_theta + seq_len(self$n_eta)
      tape <- if (isTRUE(interaction)) self$objective_tape else self$noninteraction_tape
      native <- tryCatch(
        .liberation_objective_tape_eta_mode(
          tape$pointer,
          self$objective_point(theta, start, sigma, omega),
          eta_positions, as.numeric(start), as.integer(maxit),
          as.numeric(tolerance), isTRUE(exact_hessian)
        ),
        error = identity
      )
      if (!inherits(native, "error") && identical(as.integer(native$convergence), 0L)) {
        curvature <- if (isTRUE(exact_hessian)) {
          .nm_positive_definite(native$hessian, "Conditional ETA curvature")
        } else list(
          matrix = matrix(numeric(), 0L, 0L), logdet = 0, jitter = 0
        )
        return(list(
          par = as.numeric(native$par), value = as.numeric(native$value),
          convergence = 0L, hessian = curvature$matrix,
          logdet = curvature$logdet, jitter = curvature$jitter,
          gradient = as.numeric(native$gradient),
          iterations = as.integer(native$iterations),
          evaluations = as.integer(native$evaluations), backend = "cpp"
        ))
      }
      fn <- function(eta) {
        value <- tryCatch(
          self$objective(theta, eta, sigma, omega, gradient = FALSE,
                         interaction = interaction)$value,
          error = function(e) Inf
        )
        if (is.finite(value)) value else .Machine$double.xmax / 1e100
      }
      gr <- function(eta) {
        result <- self$objective(theta, eta, sigma, omega, gradient = TRUE,
                                 interaction = interaction)
        unname(result$gradient[eta_positions])
      }
      fit <- stats::optim(
        as.numeric(start), fn, gr, method = "BFGS",
        control = list(maxit = as.integer(maxit), reltol = tolerance)
      )
      at_mode <- self$objective(
        theta, fit$par, sigma, omega, gradient = TRUE,
        hessian = isTRUE(exact_hessian), interaction = interaction
      )
      curvature <- if (isTRUE(exact_hessian)) {
        .nm_positive_definite(
          at_mode$hessian[eta_positions, eta_positions, drop = FALSE],
          "Conditional ETA curvature"
        )
      } else list(
        matrix = matrix(numeric(), 0L, 0L), logdet = 0, jitter = 0
      )
      list(
        par = fit$par, value = at_mode$value,
        convergence = fit$convergence, hessian = curvature$matrix,
        logdet = curvature$logdet, jitter = curvature$jitter,
        gradient = at_mode$gradient[eta_positions],
        iterations = as.integer(fit$counts[["gradient"]]),
        evaluations = as.integer(fit$counts[["function"]]), backend = "r-fallback"
      )
    }
  ),
  private = list(
    with_retape = function(fun, theta, sigma, omega,
                           eta = rep(0, self$n_eta)) {
      tryCatch(
        fun(),
        error = function(error) {
          if (!grepl("CppAD tape path changed", conditionMessage(error),
                     fixed = TRUE)) stop(error)
          self$record_tapes(
            theta, sigma, omega, eta = eta, retape = TRUE
          )
          fun()
        }
      )
    }
  )
)

.nm_objective_collection <- function(evaluators, parameters, eta,
                                     interaction = TRUE) {
  eta <- as.matrix(eta)
  if (!length(evaluators)) return(numeric())
  if (nrow(eta) != length(evaluators)) {
    .nm_stop("ETA rows must match the number of subject evaluators.")
  }
  if (isTRUE(evaluators[[1L]]$engine$model$USE_ODE)) {
    invisible(Map(function(evaluator, subject) {
      evaluator$ensure_valid_tapes(
        parameters$theta, parameters$sigma, parameters$omega, eta[subject, ]
      )
    }, evaluators, seq_along(evaluators)))
  }
  points <- cbind(
    matrix(parameters$theta, nrow(eta), length(parameters$theta), byrow = TRUE),
    eta,
    matrix(parameters$sigma, nrow(eta), length(parameters$sigma), byrow = TRUE),
    matrix(parameters$omega, nrow(eta), length(parameters$omega), byrow = TRUE)
  )
  tapes <- lapply(evaluators, function(evaluator) {
    if (isTRUE(interaction)) evaluator$objective_tape$pointer else
      evaluator$noninteraction_tape$pointer
  })
  .liberation_objective_tape_collection_values(tapes, points)
}

.nm_objective_collection_gradient <- function(evaluators, parameters, eta,
                                              interaction = TRUE) {
  eta <- as.matrix(eta)
  if (!length(evaluators)) return(matrix(numeric(), 0L, 0L))
  if (nrow(eta) != length(evaluators)) {
    .nm_stop("ETA rows must match the number of subject evaluators.")
  }
  if (isTRUE(evaluators[[1L]]$engine$model$USE_ODE)) {
    invisible(Map(function(evaluator, subject) {
      evaluator$ensure_valid_tapes(
        parameters$theta, parameters$sigma, parameters$omega, eta[subject, ]
      )
    }, evaluators, seq_along(evaluators)))
  }
  points <- cbind(
    matrix(parameters$theta, nrow(eta), length(parameters$theta), byrow = TRUE),
    eta,
    matrix(parameters$sigma, nrow(eta), length(parameters$sigma), byrow = TRUE),
    matrix(parameters$omega, nrow(eta), length(parameters$omega), byrow = TRUE)
  )
  tapes <- lapply(evaluators, function(evaluator) {
    if (isTRUE(interaction)) evaluator$objective_tape$pointer else
      evaluator$noninteraction_tape$pointer
  })
  .liberation_objective_tape_collection_gradients(tapes, points)
}

.nm_estimation_context <- function(model, data, n_cores = 1L, method = NULL) {
  engine <- if (inherits(model, "NMEngine")) model else nm_compile(model)
  model <- engine$model
  data <- .nm_engine_data(model, data)
  user_likelihood <- identical(model$LIK_CONFIG$error, "likelihood")
  if (model$LIK_CONFIG$error == "none" || (!user_likelihood && !nrow(model$SIGMAS))) {
    .nm_stop("Estimation requires a residual error model or a compiled user likelihood.")
  }
  omega_diagonal <- model$OMEGAS$ROW == model$OMEGAS$COL
  if ((!user_likelihood && any(model$SIGMAS$Value <= 0)) ||
      any(model$OMEGAS$Value[omega_diagonal] <= 0)) {
    .nm_stop("Initial SIGMA and diagonal OMEGA values must be positive.")
  }
  n_subjects <- length(unique(data$.ID_INDEX))
  n_cores <- as.integer(n_cores)
  if (length(n_cores) != 1L || is.na(n_cores) || n_cores < 1L) {
    .nm_stop("`n_cores` must be a positive integer.")
  }
  n_cores <- min(n_cores, max(n_subjects, 1L))
  expanded_n_eta <- .nm_eta_columns(model, data)
  tape_profile <- if (!is.null(method) && identical(toupper(method), "FO")) "fo" else "full"
  subject_data <- lapply(seq_len(n_subjects), function(index) .nm_subject_data(data, index))
  prediction_pool <- new.env(parent = emptyenv())
  subjects <- lapply(seq_len(n_subjects), function(index) {
    prediction_tape <- .nm_prediction_pool_tape(
      prediction_pool, engine, subject_data[[index]], model$THETAS$Value,
      model$SIGMAS$Value, expanded_n_eta
    )
    .NMSubjectEvaluator$new(
      engine, subject_data[[index]], model$THETAS$Value,
      model$SIGMAS$Value, model$OMEGAS$Value, n_eta = expanded_n_eta,
      prediction_tape = prediction_tape, tape_profile = tape_profile
    )
  })
  parallel_state <- NULL
  if (n_cores > 1L) {
    starts <- floor((seq_len(n_cores) - 1L) * n_subjects / n_cores) + 1L
    ends <- floor(seq_len(n_cores) * n_subjects / n_cores)
    chunks <- Map(seq.int, starts, ends)
    cluster <- parallel::makePSOCKcluster(n_cores, outfile = "")
    # Configure the worker library search path before a closure from this
    # namespace is unserialized. Otherwise an older installed LibeRation can be
    # loaded from the default user library before the initialization body gets
    # a chance to update .libPaths().
    configure_library_paths <- function(library_paths) {
      .libPaths(unique(c(library_paths, .libPaths())))
      invisible(.libPaths())
    }
    environment(configure_library_paths) <- baseenv()
    parallel::clusterCall(cluster, configure_library_paths, .libPaths())
    initialized <- tryCatch({
      parallel::clusterApply(
        cluster, seq_along(chunks),
        function(index, data_chunks, specification, theta, sigma, omega,
                 n_eta, tape_profile, library_paths) {
          .libPaths(unique(c(library_paths, .libPaths())))
          namespace <- asNamespace("LibeRation")
          compiler <- get("nm_compile", envir = namespace)
          evaluator_class <- get(".NMSubjectEvaluator", envir = namespace)
          prediction_pool_tape <- get(".nm_prediction_pool_tape", envir = namespace)
          compiled <- compiler(specification)
          prediction_pool <- new.env(parent = emptyenv())
          evaluators <- lapply(data_chunks[[index]], function(subject_data) {
            prediction_tape <- prediction_pool_tape(
              prediction_pool, compiled, subject_data, theta, sigma, n_eta
            )
            evaluator_class$new(
              compiled, subject_data, theta, sigma, omega, n_eta = n_eta,
              prediction_tape = prediction_tape, tape_profile = tape_profile
            )
          })
          assign(".liber_parallel_subjects", evaluators, envir = .GlobalEnv)
          assign(".liber_parallel_model", specification, envir = .GlobalEnv)
          TRUE
        },
        data_chunks = lapply(chunks, function(rows) subject_data[rows]),
        specification = model, theta = model$THETAS$Value,
        sigma = model$SIGMAS$Value, omega = model$OMEGAS$Value,
        n_eta = expanded_n_eta, tape_profile = tape_profile,
        library_paths = .libPaths()
      )
      TRUE
    }, error = identity)
    if (inherits(initialized, "error")) {
      try(parallel::stopCluster(cluster), silent = TRUE)
      .nm_stop("Unable to initialize parallel estimation workers: ",
               conditionMessage(initialized))
    }
    parallel_state <- list(cluster = cluster, chunks = chunks, n_cores = n_cores)
  }
  list(engine = engine, model = model, data = data, subjects = subjects,
       n_subjects = n_subjects, n_eta = expanded_n_eta,
       parallel = parallel_state)
}

.nm_outer_map <- function(model) {
  theta_fixed <- model$THETAS$FIX
  theta_free <- which(!theta_fixed)
  sigma_free <- which(!model$SIGMAS$FIX)
  omega_full <- any(model$OMEGAS$ROW != model$OMEGAS$COL)
  if (omega_full && any(model$OMEGAS$FIX) && !all(model$OMEGAS$FIX)) {
    .nm_stop("A correlated OMEGA must currently be either entirely fixed or entirely estimated.")
  }
  omega_free <- if (omega_full && all(model$OMEGAS$FIX)) integer() else {
    which(!model$OMEGAS$FIX)
  }
  omega_encode <- function(values) {
    if (!length(omega_free)) return(numeric())
    if (!omega_full) return(log(values[omega_free]))
    lower <- t(chol(.nm_omega_matrix(model, values)))
    vapply(seq_len(nrow(model$OMEGAS)), function(i) {
      row <- model$OMEGAS$ROW[[i]]
      column <- model$OMEGAS$COL[[i]]
      if (row == column) log(lower[row, column]) else lower[row, column]
    }, numeric(1))
  }
  omega_decode <- function(encoded) {
    if (!length(omega_free)) return(model$OMEGAS$Value)
    if (!omega_full) {
      values <- model$OMEGAS$Value
      values[omega_free] <- exp(encoded)
      return(values)
    }
    lower <- matrix(0, model$n_eta, model$n_eta)
    for (i in seq_len(nrow(model$OMEGAS))) {
      row <- model$OMEGAS$ROW[[i]]
      column <- model$OMEGAS$COL[[i]]
      lower[row, column] <- if (row == column) exp(encoded[[i]]) else encoded[[i]]
    }
    covariance <- lower %*% t(lower)
    vapply(seq_len(nrow(model$OMEGAS)), function(i) {
      covariance[model$OMEGAS$ROW[[i]], model$OMEGAS$COL[[i]]]
    }, numeric(1))
  }
  start <- c(model$THETAS$Value[theta_free],
             log(model$SIGMAS$Value[sigma_free]),
             omega_encode(model$OMEGAS$Value))
  theta_lower <- model$THETAS$LOWER %||% rep(-Inf, nrow(model$THETAS))
  theta_upper <- model$THETAS$UPPER %||% rep(Inf, nrow(model$THETAS))
  omega_parameter_count <- if (omega_full && length(omega_free)) {
    nrow(model$OMEGAS)
  } else length(omega_free)
  lower <- c(theta_lower[theta_free], rep(-Inf, length(sigma_free) + omega_parameter_count))
  upper <- c(theta_upper[theta_free], rep(Inf, length(sigma_free) + omega_parameter_count))
  parameter_names <- c(
    if (length(theta_free)) paste0("THETA", theta_free) else character(),
    if (length(sigma_free)) paste0("log_SIGMA", sigma_free) else character(),
    if (omega_parameter_count) {
      if (omega_full) paste0("OMEGA_CHOL", seq_len(omega_parameter_count))
      else paste0("log_OMEGA", omega_free)
    } else character()
  )
  decode <- function(parameters) {
    cursor <- 0L
    theta <- model$THETAS$Value
    if (length(theta_free)) {
      theta[theta_free] <- parameters[seq_len(length(theta_free))]
      cursor <- length(theta_free)
    }
    sigma <- model$SIGMAS$Value
    if (length(sigma_free)) {
      sigma_index <- cursor + seq_len(length(sigma_free))
      sigma[sigma_free] <- exp(parameters[sigma_index])
      cursor <- cursor + length(sigma_free)
    }
    omega <- model$OMEGAS$Value
    omega_parameter_count <- if (omega_full && length(omega_free)) {
      nrow(model$OMEGAS)
    } else length(omega_free)
    if (omega_parameter_count) {
      omega_index <- cursor + seq_len(omega_parameter_count)
      omega <- omega_decode(parameters[omega_index])
    }
    list(theta = theta, sigma = sigma, omega = omega)
  }
  encode <- function(parameters) c(
    parameters$theta[theta_free],
    log(parameters$sigma[sigma_free]),
    omega_encode(parameters$omega)
  )
  log_jacobian <- function(parameters) {
    value <- sum(log(parameters$sigma[sigma_free]))
    if (!length(omega_free)) return(value)
    if (!omega_full) return(value + sum(log(parameters$omega[omega_free])))
    lower <- t(chol(.nm_omega_matrix(model, parameters$omega)))
    value + model$n_eta * log(2) + sum(
      (model$n_eta + 2L - seq_len(model$n_eta)) * log(diag(lower))
    )
  }
  log_jacobian_gradient <- function(parameters) {
    result <- numeric(length(start))
    cursor <- length(theta_free)
    if (length(sigma_free)) {
      result[cursor + seq_len(length(sigma_free))] <- 1
      cursor <- cursor + length(sigma_free)
    }
    if (!length(omega_free)) return(result)
    if (!omega_full) {
      result[cursor + seq_len(length(omega_free))] <- 1
      return(result)
    }
    for (encoded in seq_len(nrow(model$OMEGAS))) {
      row <- model$OMEGAS$ROW[[encoded]]
      column <- model$OMEGAS$COL[[encoded]]
      if (row == column) {
        result[cursor + encoded] <- model$n_eta + 2L - row
      }
    }
    result
  }
  in_bounds <- function(parameters) {
    length(parameters) == length(lower) &&
      all(parameters >= lower) && all(parameters <= upper)
  }
  jacobian <- function(parameters) {
    n_native <- nrow(model$THETAS) + nrow(model$SIGMAS) + nrow(model$OMEGAS)
    result <- matrix(0, n_native, length(start))
    cursor <- 0L
    for (index in theta_free) {
      cursor <- cursor + 1L
      result[index, cursor] <- 1
    }
    sigma_offset <- nrow(model$THETAS)
    for (index in sigma_free) {
      cursor <- cursor + 1L
      result[sigma_offset + index, cursor] <- parameters$sigma[[index]]
    }
    omega_offset <- sigma_offset + nrow(model$SIGMAS)
    if (!length(omega_free)) return(result)
    if (!omega_full) {
      for (index in omega_free) {
        cursor <- cursor + 1L
        result[omega_offset + index, cursor] <- parameters$omega[[index]]
      }
      return(result)
    }
    covariance <- .nm_omega_matrix(model, parameters$omega)
    lower_cholesky <- t(chol(covariance))
    for (encoded in seq_len(nrow(model$OMEGAS))) {
      cursor <- cursor + 1L
      row <- model$OMEGAS$ROW[[encoded]]
      column <- model$OMEGAS$COL[[encoded]]
      derivative_lower <- matrix(0, model$n_eta, model$n_eta)
      derivative_lower[row, column] <- if (row == column) {
        lower_cholesky[row, column]
      } else 1
      derivative <- derivative_lower %*% t(lower_cholesky) +
        lower_cholesky %*% t(derivative_lower)
      for (native in seq_len(nrow(model$OMEGAS))) {
        result[omega_offset + native, cursor] <- derivative[
          model$OMEGAS$ROW[[native]], model$OMEGAS$COL[[native]]
        ]
      }
    }
    result
  }
  list(start = start, lower = lower, upper = upper, names = parameter_names,
       decode = decode, encode = encode, in_bounds = in_bounds, theta_free = theta_free,
       sigma_free = sigma_free, omega_free = omega_free,
       omega_full = omega_full, log_jacobian = log_jacobian,
       log_jacobian_gradient = log_jacobian_gradient, jacobian = jacobian)
}

.nm_cpp_prior_config <- function(model) {
  priors <- model$LIK_CONFIG$priors
  if (is.null(priors) || !nrow(priors)) {
    return(list(
      index = integer(), family = character(), mean = numeric(), sd = numeric(),
      shape = numeric(), rate = numeric()
    ))
  }
  offsets <- c(
    THETA = 0L, SIGMA = nrow(model$THETAS),
    OMEGA = nrow(model$THETAS) + nrow(model$SIGMAS)
  )
  family <- sub("[0-9]+$", "", toupper(priors$parameter))
  index <- as.integer(sub("^[A-Z]+", "", priors$parameter))
  list(
    index = unname(offsets[family]) + index,
    family = as.character(priors$distribution),
    mean = as.numeric(priors$mean), sd = as.numeric(priors$sd),
    shape = as.numeric(priors$shape), rate = as.numeric(priors$rate)
  )
}

.nm_cpp_population_objective <- function(context, map, approximation,
                                          eta_maxit, tolerance,
                                          initial_eta = NULL) {
  if (!isTRUE(getOption("LibeRation.cpp_population_objective", TRUE))) {
    return(list(pointer = NULL, reason = "disabled by option"))
  }
  if (!is.null(context$parallel)) {
    return(list(pointer = NULL, reason = "PSOCK workers require R coordination"))
  }
  approximation <- match.arg(
    tolower(approximation), c("fo", "its", "foce", "focei", "laplace")
  )
  parameters <- map$decode(map$start)
  primary <- curvature <- list()
  # Adaptive ODE tapes are owned and retaped by the C++ population object.
  # Analytical models retain the already-recorded subject tapes, including
  # structurally shared prediction tapes used to construct curvature tapes.
  if (!isTRUE(context$model$USE_ODE)) {
    if (approximation == "fo") {
      fo_pool <- new.env(parent = emptyenv())
      invisible(lapply(context$subjects, function(evaluator) {
        .nm_fo_pool_tape(
          fo_pool, evaluator, parameters$theta, parameters$sigma, parameters$omega
        )
      }))
      primary <- lapply(context$subjects, function(evaluator) evaluator$fo_tape$pointer)
    } else {
      interaction <- approximation != "foce"
      primary <- lapply(context$subjects, function(evaluator) {
        if (interaction) evaluator$objective_tape$pointer else
          evaluator$noninteraction_tape$pointer
      })
      if (approximation %in% c("foce", "focei", "laplace")) {
        curvature_anchors <- if (approximation == "laplace" && context$n_eta) {
          initial_modes <- .nm_subject_modes(
            context, parameters, starts = initial_eta,
            maxit = eta_maxit, tolerance = tolerance,
            interaction = TRUE, exact_hessian = TRUE
          )
          if (any(vapply(initial_modes, `[[`, integer(1), "convergence") != 0L)) {
            .nm_stop("Initial conditional modes did not converge for the compiled Laplace objective.")
          }
          lapply(initial_modes, `[[`, "par")
        } else {
          if (is.null(initial_eta)) {
            rep(list(rep(0, context$n_eta)), context$n_subjects)
          } else {
            lapply(seq_len(context$n_subjects), function(subject) initial_eta[subject, ])
          }
        }
        invisible(Map(function(evaluator, eta) {
          evaluator$ensure_curvature_tape(
            parameters$theta, eta, parameters$sigma,
            parameters$omega, approximation
          )
        }, context$subjects, curvature_anchors))
        curvature <- lapply(context$subjects, function(evaluator) {
          evaluator$curvature_tapes[[approximation]]$pointer
        })
      }
    }
  }
  priors <- .nm_cpp_prior_config(context$model)
  config <- list(
    approximation = approximation,
    theta = parameters$theta, sigma = parameters$sigma, omega = parameters$omega,
    theta_free = map$theta_free, sigma_free = map$sigma_free,
    omega_free = map$omega_free, omega_full = map$omega_full,
    omega_rows = context$model$OMEGAS$ROW,
    omega_cols = context$model$OMEGAS$COL,
    n_eta = context$n_eta, n_eta_base = context$model$n_eta,
    eta_maxit = as.integer(eta_maxit), tolerance = as.numeric(tolerance),
    use_ode = isTRUE(context$model$USE_ODE),
    guard_radius = as.numeric(getOption("LibeRation.tape_guard_radius", 0.5)),
    start = map$start,
    eta_start = initial_eta %||% matrix(0, context$n_subjects, context$n_eta),
    prior_index = priors$index, prior_family = priors$family,
    prior_mean = priors$mean, prior_sd = priors$sd,
    prior_shape = priors$shape, prior_rate = priors$rate,
    fo_population_batch = isTRUE(getOption("LibeRation.fo_population_batch", TRUE)),
    fo_population_max_operations = as.numeric(getOption(
      "LibeRation.fo_population_max_operations", 2e6
    ))
  )
  tryCatch(
    list(
      pointer = .liberation_population_objective_create(
        context$engine$pointer,
        lapply(context$subjects, function(evaluator) evaluator$data),
        primary, curvature, config
      ),
      reason = NULL
    ),
    error = function(error) list(
      pointer = NULL,
      reason = paste("compiled population initialization failed:", conditionMessage(error))
    )
  )
}

.nm_log_gradient <- function(iteration, objective, parameters, map, value = NULL,
                             gradient_function = NULL) {
  baseline <- value %||% objective(parameters)
  gradient <- if (is.function(gradient_function)) {
    as.numeric(gradient_function(parameters))
  } else {
    result <- numeric(length(parameters))
    for (index in seq_along(parameters)) {
      step <- 1e-5 * max(abs(parameters[[index]]), 1)
      low <- high <- parameters
      low[[index]] <- max(map$lower[[index]], parameters[[index]] - step)
      high[[index]] <- min(map$upper[[index]], parameters[[index]] + step)
      low_value <- if (low[[index]] < parameters[[index]]) objective(low) else baseline
      high_value <- if (high[[index]] > parameters[[index]]) objective(high) else baseline
      width <- high[[index]] - low[[index]]
      result[[index]] <- if (width > 0 && is.finite(low_value) && is.finite(high_value)) {
        (high_value - low_value) / width
      } else NA_real_
    }
    result
  }
  names(gradient) <- map$names
  cat(sprintf(
    "[LibeRation] OUTER EVALUATION %d OFV %.10g SCALED GRADIENT %s\n",
    as.integer(iteration), as.numeric(baseline),
    paste(sprintf("%s=%.6g", names(gradient), gradient), collapse = " ")
  ))
  try(flush(stdout()), silent = TRUE)
  invisible(gradient)
}

.nm_outer_optim <- function(map, objective, maxit, tolerance, trace = 0L,
                            print_every = 0L, gradient = NULL,
                            optimizer_backend = c("auto", "native", "r"),
                            compiled_objective = NULL,
                            strict_convergence = FALSE) {
  optimizer_backend <- match.arg(optimizer_backend)
  if (optimizer_backend == "auto") optimizer_backend <- "r"
  compiled_pointer <- compiled_objective$pointer %||% NULL
  compiled <- !is.null(compiled_pointer)
  started <- proc.time()[["elapsed"]]
  print_every <- as.integer(print_every)
  evaluations <- 0L
  gradient_evaluations <- 0L
  pending_log <- NULL
  objective_scale <- 1
  raw <- function(parameters) {
    if (!map$in_bounds(parameters)) return(1e100)
    value <- tryCatch(
      if (compiled) {
        .liberation_population_objective_value(compiled_pointer, parameters)
      } else objective(map$decode(parameters)),
      error = function(e) Inf
    )
    if (length(value) != 1L || !is.finite(value)) 1e100 else value
  }
  safe <- function(parameters) {
    evaluations <<- evaluations + 1L
    value <- raw(parameters)
    if (print_every > 0L &&
        (evaluations == 1L || evaluations %% print_every == 0L)) {
      if (is.function(gradient) || compiled) {
        pending_log <<- list(
          iteration = evaluations, parameters = parameters, value = value
        )
      } else {
        .nm_log_gradient(evaluations, raw, parameters, map, value)
      }
    }
    value
  }
  safe_gradient <- if (is.function(gradient) || compiled) function(parameters) {
    gradient_evaluations <<- gradient_evaluations + 1L
    value <- tryCatch(
      if (compiled) {
        as.numeric(.liberation_population_objective_gradient(
          compiled_pointer, parameters
        ))
      } else as.numeric(gradient(map$decode(parameters))),
      error = function(error) rep(NA_real_, length(parameters))
    )
    if (length(value) != length(parameters) || any(!is.finite(value))) {
      # L-BFGS-B evaluates the gradient at its generalized Cauchy point before
      # line-searching back. A large likelihood gradient can put that point on
      # a numerically invalid boundary even when the starting point is sound.
      # Return a finite inward barrier direction there; retain the hard error
      # for a non-finite derivative at a finite objective value.
      point_value <- raw(parameters)
      if (is.finite(point_value) && point_value >= 1e99) {
        parameter_scale <- pmax(abs(map$start), 1)
        inward <- (parameters - map$start) / parameter_scale
        largest <- max(abs(inward))
        if (is.finite(largest) && largest > 0) {
          return(objective_scale * inward / largest)
        }
      }
      .nm_stop("The population objective gradient is not finite.")
    }
    if (!is.null(pending_log) &&
        identical(as.numeric(parameters), as.numeric(pending_log$parameters))) {
      .nm_log_gradient(
        pending_log$iteration, raw, parameters, map, pending_log$value,
        gradient_function = function(ignored) value
      )
      pending_log <<- NULL
    }
    value
  } else NULL
  if (!length(map$start)) {
    result <- list(
      par = numeric(), value = safe(numeric()), convergence = 0L,
      counts = c(`function` = 1L, gradient = NA_integer_),
      iterations = 0L, objective_evaluations = 1L,
      gradient_evaluations = 0L, backend = "fixed-parameters",
      elapsed_seconds = unname(proc.time()[["elapsed"]] - started),
      message = NULL,
      objective_backend = if (compiled) "persistent-cpp-population-objective" else
        "r-orchestrated-population-objective",
      population_objective = if (compiled) {
        .liberation_population_objective_telemetry(compiled_pointer)
      } else NULL
    )
    if (compiled) {
      result$objective_backend <- result$population_objective$backend %||%
        result$objective_backend
    }
    return(result)
  }
  if (optimizer_backend == "native" && is.function(safe_gradient)) {
    result <- .liberation_native_optimizer(
      safe, safe_gradient, map$start, map$lower, map$upper,
      as.integer(maxit), as.numeric(tolerance), as.integer(trace)
    )
    result$backend <- "native-bfgs"
    result$objective_backend <- if (compiled) {
      "persistent-cpp-population-objective"
    } else "r-orchestrated-population-objective"
    result$population_objective <- if (compiled) {
      .liberation_population_objective_telemetry(compiled_pointer)
    } else NULL
    if (compiled) {
      result$objective_backend <- result$population_objective$backend %||%
        result$objective_backend
    }
    result$elapsed_seconds <- unname(proc.time()[["elapsed"]] - started)
    return(result)
  }
  bounded <- any(is.finite(map$lower)) || any(is.finite(map$upper))
  initial_value <- raw(map$start)
  if (is.finite(initial_value) && initial_value < 1e99) {
    objective_scale <- max(abs(initial_value), 1)
  }
  arguments <- list(
    par = map$start, fn = safe,
    method = if (bounded) "L-BFGS-B" else if (is.function(safe_gradient) ||
      length(map$start) == 1L) "BFGS" else "Nelder-Mead",
    control = list(
      maxit = as.integer(maxit), reltol = tolerance, trace = trace,
      fnscale = objective_scale
    )
  )
  if (is.function(safe_gradient)) arguments$gr <- safe_gradient
  if (bounded) {
    arguments$lower <- map$lower
    arguments$upper <- map$upper
    if (isTRUE(strict_convergence)) {
      # Exact FO gradients can be very large at the initial point. A loose
      # function-reduction test stopped before OMEGA reached the same solution
      # under mathematically equivalent summation orders. Keep the function
      # test at machine accuracy and use a squared scaled-gradient target.
      arguments$control$factr <- 1
      arguments$control$pgtol <- tolerance^2
    } else {
      arguments$control$factr <- max(tolerance / .Machine$double.eps, 1)
    }
    arguments$control$reltol <- NULL
  }
  result <- do.call(stats::optim, arguments)
  iterations <- suppressWarnings(as.integer(result$counts[["gradient"]]))
  if (!length(iterations) || is.na(iterations)) {
    iterations <- suppressWarnings(as.integer(result$counts[["function"]]))
  }
  result$iterations <- if (!length(iterations) || is.na(iterations)) {
    as.integer(evaluations)
  } else iterations
  result$objective_evaluations <- as.integer(evaluations)
  result$gradient_evaluations <- as.integer(gradient_evaluations)
  result$backend <- if (compiled) {
    paste0("r-", tolower(arguments$method), "-cpp-objective")
  } else if (is.function(safe_gradient)) "r-optim-gradient" else "r-optim"
  result$population_objective <- if (compiled) {
    .liberation_population_objective_telemetry(compiled_pointer)
  } else NULL
  result$objective_backend <- if (compiled) {
    result$population_objective$backend %||% "persistent-cpp-population-objective"
  } else "r-orchestrated-population-objective"
  result$objective_scale <- objective_scale
  result$elapsed_seconds <- unname(proc.time()[["elapsed"]] - started)
  result
}

.nm_subject_modes <- function(context, parameters, starts = NULL,
                              maxit = 100L, tolerance = 1e-7,
                              interaction = TRUE, exact_hessian = TRUE) {
  if (is.null(starts)) starts <- matrix(0, context$n_subjects, context$n_eta)
  batch_modes <- function(evaluators, starts) {
    if (!length(evaluators)) return(list())
    ode_guard <- isTRUE(evaluators[[1L]]$engine$model$USE_ODE)
    if (ode_guard) invisible(Map(function(evaluator, subject) {
      evaluator$ensure_valid_tapes(
        parameters$theta, parameters$sigma, parameters$omega, starts[subject, ]
      )
    }, evaluators, seq_along(evaluators)))
    points <- cbind(
      matrix(parameters$theta, nrow(starts), length(parameters$theta), byrow = TRUE),
      starts,
      matrix(parameters$sigma, nrow(starts), length(parameters$sigma), byrow = TRUE),
      matrix(parameters$omega, nrow(starts), length(parameters$omega), byrow = TRUE)
    )
    tapes <- lapply(evaluators, function(evaluator) {
      if (isTRUE(interaction)) evaluator$objective_tape$pointer else
        evaluator$noninteraction_tape$pointer
    })
    raw <- tryCatch(
      .liberation_objective_tape_eta_modes(
        tapes, points, length(parameters$theta) + seq_len(context$n_eta), starts,
        as.integer(maxit), as.numeric(tolerance), isTRUE(exact_hessian)
      ), error = identity
    )
    if (inherits(raw, "error")) {
      if (!grepl("CppAD tape path changed", conditionMessage(raw), fixed = TRUE)) {
        stop(raw)
      }
      return(lapply(seq_along(evaluators), function(subject) {
        evaluators[[subject]]$eta_mode(
          parameters$theta, parameters$sigma, parameters$omega,
          start = starts[subject, ], maxit = maxit, tolerance = tolerance,
          interaction = interaction, exact_hessian = exact_hessian
        )
      }))
    }
    lapply(seq_along(raw), function(subject) {
      mode <- raw[[subject]]
      if (!identical(as.integer(mode$convergence), 0L)) {
        return(evaluators[[subject]]$eta_mode(
          parameters$theta, parameters$sigma, parameters$omega,
          start = starts[subject, ], maxit = maxit, tolerance = tolerance,
          interaction = interaction, exact_hessian = exact_hessian
        ))
      }
      if (ode_guard && evaluators[[subject]]$ensure_valid_tapes(
        parameters$theta, parameters$sigma, parameters$omega, mode$par
      )) {
        return(evaluators[[subject]]$eta_mode(
          parameters$theta, parameters$sigma, parameters$omega,
          start = mode$par, maxit = maxit, tolerance = tolerance,
          interaction = interaction, exact_hessian = exact_hessian
        ))
      }
      curvature <- if (isTRUE(exact_hessian)) {
        .nm_positive_definite(mode$hessian, "Conditional ETA curvature")
      } else list(matrix = matrix(numeric(), 0L, 0L), logdet = 0, jitter = 0)
      list(
        par = as.numeric(mode$par), value = as.numeric(mode$value),
        convergence = 0L, hessian = curvature$matrix,
        logdet = curvature$logdet, jitter = curvature$jitter,
        gradient = as.numeric(mode$gradient),
        iterations = as.integer(mode$iterations),
        evaluations = as.integer(mode$evaluations), backend = "cpp-batch"
      )
    })
  }
  if (!is.null(context$parallel)) {
    chunks <- context$parallel$chunks
    pieces <- parallel::clusterApply(
      context$parallel$cluster, seq_along(chunks),
       function(index, start_chunks, theta, sigma, omega, maxit, tolerance,
                interaction, exact_hessian) {
        evaluators <- get(".liber_parallel_subjects", envir = .GlobalEnv)
        worker_starts <- start_chunks[[index]]
        context <- list(n_eta = ncol(worker_starts))
        parameters <- list(theta = theta, sigma = sigma, omega = omega)
        batch <- get(".nm_subject_modes_batch", envir = asNamespace("LibeRation"))
        batch(evaluators, context, parameters, worker_starts, maxit, tolerance,
              interaction, exact_hessian)
      },
      start_chunks = lapply(chunks, function(rows) starts[rows, , drop = FALSE]),
      theta = parameters$theta, sigma = parameters$sigma,
       omega = parameters$omega, maxit = maxit, tolerance = tolerance,
       interaction = interaction, exact_hessian = exact_hessian
    )
    return(unlist(pieces, recursive = FALSE))
  }
  batch_modes(context$subjects, starts)
}

.nm_subject_modes_batch <- function(evaluators, context, parameters, starts,
                                    maxit, tolerance, interaction,
                                    exact_hessian) {
  if (!length(evaluators)) return(list())
  ode_guard <- isTRUE(evaluators[[1L]]$engine$model$USE_ODE)
  if (ode_guard) invisible(Map(function(evaluator, subject) {
    evaluator$ensure_valid_tapes(
      parameters$theta, parameters$sigma, parameters$omega, starts[subject, ]
    )
  }, evaluators, seq_along(evaluators)))
  points <- cbind(
    matrix(parameters$theta, nrow(starts), length(parameters$theta), byrow = TRUE),
    starts,
    matrix(parameters$sigma, nrow(starts), length(parameters$sigma), byrow = TRUE),
    matrix(parameters$omega, nrow(starts), length(parameters$omega), byrow = TRUE)
  )
  tapes <- lapply(evaluators, function(evaluator) {
    if (isTRUE(interaction)) evaluator$objective_tape$pointer else
      evaluator$noninteraction_tape$pointer
  })
  raw <- tryCatch(
    .liberation_objective_tape_eta_modes(
      tapes, points, length(parameters$theta) + seq_len(context$n_eta), starts,
      as.integer(maxit), as.numeric(tolerance), isTRUE(exact_hessian)
    ), error = identity
  )
  if (inherits(raw, "error")) {
    if (!grepl("CppAD tape path changed", conditionMessage(raw), fixed = TRUE)) {
      stop(raw)
    }
    return(lapply(seq_along(evaluators), function(subject) {
      evaluators[[subject]]$eta_mode(
        parameters$theta, parameters$sigma, parameters$omega,
        start = starts[subject, ], maxit = maxit, tolerance = tolerance,
        interaction = interaction, exact_hessian = exact_hessian
      )
    }))
  }
  lapply(seq_along(raw), function(subject) {
    mode <- raw[[subject]]
    if (!identical(as.integer(mode$convergence), 0L) ||
        (ode_guard && evaluators[[subject]]$ensure_valid_tapes(
          parameters$theta, parameters$sigma, parameters$omega, mode$par
        ))) {
      return(evaluators[[subject]]$eta_mode(
        parameters$theta, parameters$sigma, parameters$omega,
        start = mode$par, maxit = maxit, tolerance = tolerance,
        interaction = interaction, exact_hessian = exact_hessian
      ))
    }
    curvature <- if (isTRUE(exact_hessian)) {
      .nm_positive_definite(mode$hessian, "Conditional ETA curvature")
    } else list(matrix = matrix(numeric(), 0L, 0L), logdet = 0, jitter = 0)
    list(
      par = as.numeric(mode$par), value = as.numeric(mode$value), convergence = 0L,
      hessian = curvature$matrix, logdet = curvature$logdet,
      jitter = curvature$jitter, gradient = as.numeric(mode$gradient),
      iterations = as.integer(mode$iterations), evaluations = as.integer(mode$evaluations),
      backend = "cpp-batch"
    )
  })
}

.nm_subject_curvature_logdet <- function(context, evaluator, parameters, eta,
                                         approximation) {
  eta_columns <- length(parameters$theta) + seq_len(context$n_eta)
  if (approximation == "laplace") {
    curvature <- evaluator$objective_hessian_subset(
      parameters$theta, eta, parameters$sigma, parameters$omega,
      rows = eta_columns, columns = eta_columns, interaction = TRUE
    )
    return(.nm_positive_definite(
      curvature, "Laplace conditional curvature"
    )$logdet)
  }
  prediction <- evaluator$prediction(
    parameters$theta, eta, parameters$sigma, jacobian = TRUE,
    columns = eta_columns
  )
  observed <- evaluator$data$EVID == 0L & evaluator$data$MDV == 0L &
    is.finite(evaluator$data$DV)
  jacobian <- prediction$jacobian[observed, , drop = FALSE]
  f <- prediction$value[observed]
  scale_f <- if (approximation == "foce") {
    evaluator$prediction(
      parameters$theta, rep(0, context$n_eta), parameters$sigma,
      jacobian = FALSE
    )$value[observed]
  } else f
  dvid <- if ("DVID" %in% names(evaluator$data)) {
    evaluator$data$DVID[observed]
  } else rep(1L, sum(observed))
  variance <- .nm_residual_variance(
    context$model, scale_f, parameters$sigma, dvid
  )
  omega_inverse <- solve(.nm_effect_covariance(
    context$model, evaluator$data, parameters$omega
  ))
  curvature <- 2 * crossprod(jacobian / sqrt(variance)) + 2 * omega_inverse
  .nm_positive_definite(
    curvature, paste0(toupper(approximation), " Gauss-Newton curvature")
  )$logdet
}

.nm_conditional_native_gradient <- function(context, parameters, eta,
                                            interaction = TRUE) {
  if (is.null(context$parallel)) {
    gradients <- .nm_objective_collection_gradient(
      context$subjects, parameters, eta, interaction = interaction
    )
    total <- colSums(gradients)
  } else {
    eta_chunks <- lapply(
      context$parallel$chunks, function(rows) eta[rows, , drop = FALSE]
    )
    pieces <- parallel::clusterApply(
      context$parallel$cluster, seq_along(context$parallel$chunks),
      function(index, eta_chunks, parameters, interaction) {
        evaluators <- get(".liber_parallel_subjects", envir = .GlobalEnv)
        collection <- get(
          ".nm_objective_collection_gradient", envir = asNamespace("LibeRation")
        )
        colSums(collection(
          evaluators, parameters, eta_chunks[[index]], interaction = interaction
        ))
      }, eta_chunks = eta_chunks, parameters = parameters,
      interaction = interaction
    )
    total <- Reduce(`+`, pieces)
  }
  n_theta <- length(parameters$theta)
  n_sigma <- length(parameters$sigma)
  n_omega <- length(parameters$omega)
  population_positions <- c(
    seq_len(n_theta),
    n_theta + context$n_eta + seq_len(n_sigma),
    n_theta + context$n_eta + n_sigma + seq_len(n_omega)
  )
  as.numeric(total[population_positions]) +
    .nm_prior_nll_native_gradient(context$model, parameters)
}

.nm_nested_outer_gradient <- function(context, map, objective, parameters,
                                      approximation, relative_step = 1e-5) {
  value <- objective(parameters)
  if (!is.finite(value)) .nm_stop("Cannot differentiate a non-finite objective.")
  state <- attr(objective, "state")
  modes <- state$modes
  if (is.null(modes)) .nm_stop("Conditional modes are unavailable for differentiation.")
  eta <- if (context$n_eta) {
    do.call(rbind, lapply(modes, `[[`, "par"))
  } else matrix(numeric(), context$n_subjects, 0L)
  interaction <- approximation != "foce"
  transform <- map$jacobian(parameters)
  if (approximation == "its" || !context$n_eta || !ncol(transform)) {
    native <- .nm_conditional_native_gradient(
      context, parameters, eta, interaction = interaction
    )
    return(as.vector(native %*% transform))
  }
  if (is.null(context$parallel)) {
    result <- .nm_nested_gradient_batch(
      context$subjects, context$n_eta, parameters, eta, approximation, transform
    )
    gradient <- result$gradient
  } else {
    eta_chunks <- lapply(
      context$parallel$chunks, function(rows) eta[rows, , drop = FALSE]
    )
    pieces <- parallel::clusterApply(
      context$parallel$cluster, seq_along(context$parallel$chunks),
      function(index, eta_chunks, n_eta, parameters, approximation, transform) {
        evaluators <- get(".liber_parallel_subjects", envir = .GlobalEnv)
        batch <- get(".nm_nested_gradient_batch", envir = asNamespace("LibeRation"))
        batch(
          evaluators, n_eta, parameters, eta_chunks[[index]],
          approximation, transform
        )$gradient
      }, eta_chunks = eta_chunks, n_eta = context$n_eta,
      parameters = parameters, approximation = approximation,
      transform = transform
    )
    gradient <- Reduce(`+`, pieces)
  }
  prior <- .nm_prior_nll_native_gradient(context$model, parameters)
  as.numeric(gradient) + as.vector(prior %*% transform)
}

.nm_nested_gradient_batch <- function(evaluators, n_eta, parameters, eta,
                                      approximation, transform) {
  interaction <- approximation != "foce"
  ode_guard <- isTRUE(evaluators[[1L]]$engine$model$USE_ODE)
  invisible(Map(function(evaluator, subject) {
    if (ode_guard) evaluator$ensure_valid_tapes(
      parameters$theta, parameters$sigma, parameters$omega, eta[subject, ]
    )
    evaluator$ensure_curvature_tape(
      parameters$theta, eta[subject, ], parameters$sigma,
      parameters$omega, approximation
    )
  }, evaluators, seq_along(evaluators)))
  points <- cbind(
    matrix(parameters$theta, nrow(eta), length(parameters$theta), byrow = TRUE),
    eta,
    matrix(parameters$sigma, nrow(eta), length(parameters$sigma), byrow = TRUE),
    matrix(parameters$omega, nrow(eta), length(parameters$omega), byrow = TRUE)
  )
  n_theta <- length(parameters$theta)
  n_sigma <- length(parameters$sigma)
  n_omega <- length(parameters$omega)
  population_positions <- c(
    seq_len(n_theta), n_theta + n_eta + seq_len(n_sigma),
    n_theta + n_eta + n_sigma + seq_len(n_omega)
  )
  .liberation_nested_population_gradient(
    lapply(evaluators, function(evaluator) {
      if (interaction) evaluator$objective_tape$pointer else
        evaluator$noninteraction_tape$pointer
    }),
    lapply(evaluators, function(evaluator) {
      evaluator$curvature_tapes[[approximation]]$pointer
    }),
    points, n_theta + seq_len(n_eta), population_positions, transform
  )
}

.nm_nested_objective <- function(context, approximation, eta_maxit, tolerance,
                                 initial_eta = NULL) {
  force(context); force(approximation)
  state <- new.env(parent = emptyenv())
  state$starts <- initial_eta %||% matrix(0, context$n_subjects, context$n_eta)
  state$key <- NULL
  state$value <- NULL
  state$modes <- NULL
  state$objective_calls <- 0L
  state$cache_hits <- 0L
  state$mode_iterations <- 0L
  state$mode_evaluations <- 0L
  objective <- function(parameters) {
    state$objective_calls <- state$objective_calls + 1L
    key <- c(parameters$theta, parameters$sigma, parameters$omega)
    if (!is.null(state$key) && identical(key, state$key)) {
      state$cache_hits <- state$cache_hits + 1L
      return(state$value)
    }
    modes <- .nm_subject_modes(
      context, parameters, starts = state$starts, maxit = eta_maxit,
      tolerance = tolerance, interaction = approximation != "foce",
      exact_hessian = approximation == "laplace"
    )
    if (any(vapply(modes, `[[`, integer(1), "convergence") != 0L)) return(Inf)
    state$mode_iterations <- state$mode_iterations + sum(vapply(
      modes, function(mode) as.integer(mode$iterations %||% 0L), integer(1)
    ))
    state$mode_evaluations <- state$mode_evaluations + sum(vapply(
      modes, function(mode) as.integer(mode$evaluations %||% 0L), integer(1)
    ))
    if (context$n_eta) {
      state$starts <- do.call(rbind, lapply(modes, `[[`, "par"))
    }
    value <- sum(vapply(modes, `[[`, numeric(1), "value"))
    if (approximation == "laplace") {
      value <- value + sum(vapply(modes, `[[`, numeric(1), "logdet"))
    } else if (approximation %in% c("foce", "focei")) {
      for (subject in seq_along(modes)) {
        value <- value + .nm_subject_curvature_logdet(
          context, context$subjects[[subject]], parameters,
          modes[[subject]]$par, approximation
        )
      }
    }
    value <- value + .nm_prior_nll(context$model, parameters)
    state$key <- key
    state$value <- value
    state$modes <- modes
    state$parameters <- parameters
    value
  }
  attr(objective, "state") <- state
  objective
}

.nm_fo_subject <- function(evaluator, model, theta, sigma, omega) {
  evaluator$fo_objective(theta, sigma, omega)$value
}

.nm_fo_subject_reference <- function(evaluator, model, theta, sigma, omega) {
  eta_columns <- length(theta) + seq_len(evaluator$n_eta)
  prediction <- evaluator$prediction(
    theta, rep(0, evaluator$n_eta), sigma, jacobian = TRUE,
    columns = eta_columns
  )
  observed <- evaluator$data$EVID == 0L & evaluator$data$MDV == 0L & is.finite(evaluator$data$DV)
  f <- prediction$value[observed]
  dv <- evaluator$data$DV[observed]
  jacobian <- prediction$jacobian[observed, , drop = FALSE]
  if (model$LIK_CONFIG$error == "exponential") {
    if (any(dv <= 0) || any(f <= 0)) return(Inf)
    residual <- log(dv) - log(f)
    jacobian <- jacobian / f
  } else residual <- dv - f
  dvid <- if ("DVID" %in% names(evaluator$data)) evaluator$data$DVID[observed] else 1L
  variance <- .nm_residual_variance(model, f, sigma, dvid)
  correlation <- diag(length(f))
  if (model$LIK_CONFIG$sigma_corr == "ar1" && length(f) > 1L) {
    rho <- .nm_ar1_rho(model, theta = theta, sigma = sigma)
    correlation <- outer(seq_along(f), seq_along(f), function(i, j) {
      rho^abs(i - j)
    })
  }
  if (length(model$LIK_CONFIG$residual_groups) && length(f) > 1L) {
    observed_data <- evaluator$data[observed, , drop = FALSE]
    observed_dvid <- if ("DVID" %in% names(observed_data)) observed_data$DVID else
      rep(1L, nrow(observed_data))
    for (group in model$LIK_CONFIG$residual_groups) {
      group_correlation <- .nm_residual_group_value(group, theta, sigma)
      for (row in seq_len(nrow(observed_data))) {
        if (!observed_dvid[[row]] %in% group$dvid) next
        for (column in seq_len(nrow(observed_data))) {
          if (row == column || observed_data$.ID_INDEX[[row]] != observed_data$.ID_INDEX[[column]] ||
              observed_data$TIME[[row]] != observed_data$TIME[[column]] ||
              !observed_dvid[[column]] %in% group$dvid) next
          correlation[row, column] <- group_correlation[
            match(observed_dvid[[row]], group$dvid),
            match(observed_dvid[[column]], group$dvid)
          ]
        }
      }
    }
  }
  residual_covariance <- correlation * outer(sqrt(variance), sqrt(variance))
  marginal <- residual_covariance +
    jacobian %*% .nm_effect_covariance(model, evaluator$data, omega) %*% t(jacobian)
  pd <- .nm_positive_definite(marginal, "FO marginal covariance")
  as.numeric(pd$logdet + crossprod(residual, solve(pd$matrix, residual)))
}

.nm_fo_objective <- function(context, parameters) {
  if (is.null(context$parallel)) {
    values <- vapply(
      context$subjects, .nm_fo_subject, numeric(1),
      model = context$model, theta = parameters$theta,
      sigma = parameters$sigma, omega = parameters$omega
    )
  } else {
    pieces <- parallel::clusterApply(
      context$parallel$cluster, seq_along(context$parallel$chunks),
      function(index, parameters) {
        namespace <- asNamespace("LibeRation")
        subject_objective <- get(".nm_fo_subject", envir = namespace)
        evaluators <- get(".liber_parallel_subjects", envir = .GlobalEnv)
        model <- get(".liber_parallel_model", envir = .GlobalEnv)
        vapply(
          evaluators, subject_objective, numeric(1), model = model,
          theta = parameters$theta, sigma = parameters$sigma,
          omega = parameters$omega
        )
      }, parameters = parameters
    )
    values <- unlist(pieces, use.names = FALSE)
  }
  sum(values) + .nm_prior_nll(context$model, parameters)
}

.nm_fo_collection_gradient <- function(evaluators, parameters) {
  if (!length(evaluators)) return(matrix(numeric(), 0L, 0L))
  # Shared FO tapes hold subject data as CppAD dynamic parameters. Select each
  # subject immediately before differentiating; collecting duplicate pointers
  # first would leave all rows on the final subject's dynamic values.
  do.call(rbind, lapply(evaluators, function(evaluator) {
    unname(evaluator$fo_objective(
      parameters$theta, parameters$sigma, parameters$omega,
      gradient = TRUE
    )$gradient)
  }))
}

.nm_fo_native_gradient <- function(context, parameters) {
  if (is.null(context$parallel)) {
    total <- colSums(.nm_fo_collection_gradient(context$subjects, parameters))
  } else {
    pieces <- parallel::clusterCall(
      context$parallel$cluster,
      function(parameters) {
        evaluators <- get(".liber_parallel_subjects", envir = .GlobalEnv)
        collection <- get(".nm_fo_collection_gradient", envir = asNamespace("LibeRation"))
        colSums(collection(evaluators, parameters))
      }, parameters = parameters
    )
    total <- Reduce(`+`, pieces)
  }
  as.numeric(total) + .nm_prior_nll_native_gradient(context$model, parameters)
}

.nm_fo_outer_gradient <- function(context, map, parameters) {
  as.vector(.nm_fo_native_gradient(context, parameters) %*% map$jacobian(parameters))
}

.nm_fit_result <- function(context, method, parameters, objective, modes,
                           optimizer, diagnostics = list()) {
  eta <- if (context$n_eta) {
    do.call(rbind, lapply(modes, `[[`, "par"))
  } else matrix(numeric(), context$n_subjects, 0L)
  colnames(eta) <- if (context$n_eta) paste0("ETA", seq_len(context$n_eta)) else NULL
  iterations <- suppressWarnings(as.integer(optimizer$iterations %||% NA_integer_))
  if (!length(iterations) || is.na(iterations)) {
    iterations <- suppressWarnings(as.integer(optimizer$counts[["gradient"]] %||% NA_integer_))
  }
  if (!length(iterations) || is.na(iterations)) {
    iterations <- suppressWarnings(as.integer(optimizer$counts[["function"]] %||% NA_integer_))
  }
  objective_evaluations <- suppressWarnings(as.integer(
    optimizer$objective_evaluations %||% optimizer$counts[["function"]] %||% NA_integer_
  ))
  eta_iterations <- if (length(modes)) {
    sum(vapply(modes, function(mode) as.integer(mode$iterations %||% 0L), integer(1)))
  } else 0L
  eta_evaluations <- if (length(modes)) {
    sum(vapply(modes, function(mode) as.integer(mode$evaluations %||% 0L), integer(1)))
  } else 0L
  population_work <- optimizer$population_objective %||% list()
  if (eta_iterations == 0L && !is.null(population_work$mode_iterations)) {
    eta_iterations <- as.integer(population_work$mode_iterations)
  }
  if (eta_evaluations == 0L && !is.null(population_work$mode_evaluations)) {
    eta_evaluations <- as.integer(population_work$mode_evaluations)
  }
  tape <- lapply(context$subjects, function(evaluator) evaluator$tape_telemetry())
  diagnostics$optimizer <- list(
    backend = optimizer$backend %||% "unknown",
    objective_backend = optimizer$objective_backend %||%
      "r-orchestrated-population-objective",
    elapsed_seconds = optimizer$elapsed_seconds %||% NA_real_,
    objective_evaluations = objective_evaluations,
    gradient_evaluations = optimizer$gradient_evaluations %||% NA_integer_,
    trace = optimizer$telemetry %||% NULL,
    population_objective = optimizer$population_objective %||% NULL
  )
  diagnostics$conditional_modes <- list(
    iterations = eta_iterations, evaluations = eta_evaluations,
    backends = table(vapply(modes, function(mode) mode$backend %||% "unknown", character(1)))
  )
  diagnostics$tapes <- list(
    unique_structures = length(unique(vapply(tape, `[[`, character(1), "signature"))),
    subjects = length(tape),
    shared_prediction_tapes = max(
      0L, length(tape) - length(unique(vapply(tape, `[[`, character(1), "signature")))
    ),
    records = sum(vapply(tape, `[[`, integer(1), "records")) +
      as.integer(population_work$tape_records %||% 0L),
    retapes = sum(vapply(tape, `[[`, integer(1), "retapes")) +
      as.integer(population_work$tape_retapes %||% 0L),
    validity_checks = sum(vapply(tape, `[[`, integer(1), "validity_checks"))
  )
  structure(
    list(
      version = 1L, method = method, objective = as.numeric(objective),
      theta = parameters$theta, omega = parameters$omega,
      sigma = parameters$sigma, eta = eta,
      convergence = optimizer$convergence, message = optimizer$message,
      iterations = iterations,
      objective_evaluations = objective_evaluations,
      evaluations = optimizer$counts, model = context$model,
      data = context$data, diagnostics = diagnostics
    ),
    class = "nm_fit"
  )
}

.nm_est_fo <- function(context, map, maxit, tolerance, trace, print_every = 0L,
                       optimizer_backend = "auto") {
  if (context$model$LIK_CONFIG$blq_method != "none") {
    .nm_stop("FO does not support a censored Gaussian linearization; use FOCEI or LAPLACE for BLQ data.")
  }
  objective <- function(parameters) {
    .nm_fo_objective(context, parameters)
  }
  gradient <- function(parameters) .nm_fo_outer_gradient(context, map, parameters)
  compiled <- .nm_cpp_population_objective(
    context, map, "fo", eta_maxit = 100L, tolerance = tolerance
  )
  optimizer <- .nm_outer_optim(
    map, objective, maxit, tolerance, trace, print_every, gradient = gradient,
    optimizer_backend = optimizer_backend, compiled_objective = compiled,
    strict_convergence = TRUE
  )
  parameters <- map$decode(optimizer$par)
  modes <- .nm_subject_modes(
    context, parameters, maxit = 100L, tolerance = tolerance,
    exact_hessian = FALSE
  )
  .nm_fit_result(
    context, "FO", parameters, optimizer$value, modes, optimizer,
    diagnostics = list(population_gradient = "exact CppAD FO marginal gradient")
  )
}

.nm_est_nested <- function(context, map, method, maxit, eta_maxit, tolerance, trace,
                           print_every = 0L, optimizer_backend = "auto",
                           initial_eta = NULL) {
  approximation <- switch(
    method, FOCE = "foce", FOCEI = "focei", LAPLACE = "laplace", "laplace"
  )
  objective <- .nm_nested_objective(
    context, approximation, eta_maxit, tolerance, initial_eta = initial_eta
  )
  gradient <- function(parameters) .nm_nested_outer_gradient(
    context, map, objective, parameters, approximation
  )
  compiled <- .nm_cpp_population_objective(
    context, map, approximation, eta_maxit, tolerance, initial_eta = initial_eta
  )
  optimizer <- .nm_outer_optim(
    map, objective, maxit, tolerance, trace, print_every,
    gradient = gradient, optimizer_backend = optimizer_backend,
    compiled_objective = compiled
  )
  parameters <- map$decode(optimizer$par)
  cached <- if (!is.null(compiled$pointer)) tryCatch(
    .liberation_population_objective_state(compiled$pointer, optimizer$par),
    error = function(error) NULL
  ) else NULL
  modes <- cached$modes %||% .nm_subject_modes(
    context, parameters, starts = initial_eta,
    maxit = eta_maxit, tolerance = tolerance,
    interaction = approximation != "foce",
    exact_hessian = approximation == "laplace"
  )
  work <- optimizer$population_objective
  if (is.null(work)) {
    work <- list(
      value_requests = attr(objective, "state")$objective_calls,
      value_cache_hits = attr(objective, "state")$cache_hits,
      mode_iterations = attr(objective, "state")$mode_iterations,
      mode_evaluations = attr(objective, "state")$mode_evaluations
    )
  }
  .nm_fit_result(
    context, method, parameters, optimizer$value, modes, optimizer,
    diagnostics = list(
      eta_convergence = vapply(modes, `[[`, integer(1), "convergence"),
      eta_jitter = vapply(modes, `[[`, numeric(1), "jitter"),
      approximation = approximation,
      conditional_mode_work = list(
        objective_calls = work$value_requests %||% work$parameter_evaluations,
        cache_hits = sum(
          work$value_cache_hits %||% 0L,
          work$gradient_cache_hits %||% 0L,
          work$shared_state_hits %||% 0L
        ),
        iterations = work$mode_iterations,
        evaluations = work$mode_evaluations
      ),
      population_gradient = "exact CppAD curvature with implicit conditional-mode derivative"
    )
  )
}

#' Estimate a population pharmacometric model
#'
#' Deterministic conditional methods use exact CppAD gradients for ETA modes.
#' LAPLACE uses the exact conditional Hessian; FOCE/FOCEI use the
#' interaction-aware Gauss-Newton curvature; FO integrates the first-order
#' Gaussian linearization analytically. GQ integrates the exact joint
#' objective over ETAs with adaptive tensor or Smolyak sparse
#' Gauss--Hermite quadrature.
#' Stochastic methods use the same C++ joint objective and are implemented in
#' the stochastic estimation module. HMC and NUTS use exact joint CppAD
#' gradients with unconstrained parameter transforms, dual-averaged step-size
#' adaptation, and a diagonal mass matrix. NPML and NPAG replace the Gaussian
#' random-effect distribution with a discrete ETA support distribution; NPAG
#' adapts, expands, and prunes that support.
#' Serial FO, FOCE, FOCEI, Laplace, and ITS fits expose a persistent C++
#' population objective through thin callbacks to R's L-BFGS-B/BFGS optimizer;
#' PSOCK fits retain R coordination across their persistent C++ workers.
#'
#' @param model An [nm_model()] or compiled [NMEngine].
#' @param data A NONMEM-style event dataset.
#' @param method Estimation method.
#' @param maxit Maximum outer evaluations.
#' @param eta_maxit Maximum conditional ETA iterations.
#' @param tolerance Relative optimization tolerance.
#' @param trace Optimizer tracing level.
#' @param print_every Write the objective and scaled population gradient to
#'   stdout every N outer evaluations. Zero disables iteration logging.
#' @param n_cores Number of parallel subject workers. Parallel workers persist
#'   compiled C++ engines for the duration of the estimation.
#' @param optimizer_backend Use the experimental compiled scaled, bounded BFGS
#'   optimizer (`"native"`), R's mature L-BFGS-B/BFGS driver (`"r"`), or the
#'   calibrated production policy (`"auto"`). Auto selects R's driver while
#'   its objective, gradient, conditional-mode state, parameter transforms, and
#'   covariance derivatives remain in the persistent C++ population evaluator.
#' @param covariance Run a covariance step after estimation.
#' @param covariance_type Covariance estimator: automatic R/S fallback,
#'   objective Hessian (`"hessian"`), subject-score OPG (`"opg"`), or robust
#'   sandwich (`"sandwich"`). `"r"` and `"s"` are accepted aliases.
#' @param covariance_tolerance Positive-definite regularization tolerance for
#'   the covariance step.
#' @param covariance_samples Target integration budget for an IMP or SAEM
#'   covariance step. Low-dimensional ETA integrals use tensor Gauss--Hermite
#'   quadrature near this budget; higher-dimensional integrals use random-normal
#'   importance samples. The default reuses the IMP count or uses 200 for SAEM.
#' @param covariance_seed Common-random-number seed for the random-normal
#'   covariance fallback. The default reuses the estimation seed.
#' @param initial_eta Optional finite subject-by-ETA matrix used to warm-start
#'   compatible conditional or SAEM estimation steps.
#' @param collect_output Whether selected generated OUTPUT columns should be
#'   evaluated and retained with the completed fit.
#' @param ... Method-specific controls. For `method = "GQ"`, use `gq_grid`
#'   (`"auto"`, `"tensor"`, or `"smolyak"`), `gq_order` (tensor nodes per ETA
#'   dimension, default 5), `gq_level` (Smolyak level, default 3),
#'   `gq_adaptive` (default `TRUE`), `gq_max_points` (retained-grid allocation
#'   guard, default 100000), and `gq_gradient` (`"score"` or the slower
#'   `"finite_grid"`). Automatic selection uses tensor quadrature through
#'   three ETAs and Smolyak quadrature for higher-dimensional models.
#'   For `method = "HMC"` or `"NUTS"`, controls include `n_warmup` (500),
#'   `n_sample` (1000 per chain), `n_thin` (1), `n_chains` (4), `seed`,
#'   optional `step_size`, `target_acceptance` (0.8), `adapt_mass` (`TRUE`),
#'   `n_leapfrog` (10; HMC), `max_depth` (10; NUTS), and
#'   `divergence_threshold` (1000). For `method = "NPML"` or `"NPAG"`, use
#'   `np_supports` for an optional fixed starting matrix, `np_points` (25),
#'   `np_max_support` (100), `np_min_weight` (1e-5), `np_weight_maxit` (1000),
#'   `np_cycles` (3), and, for NPAG, `np_grid_step` (1), `np_grid_decay`
#'   (0.5), and `np_max_candidates` (500). `np_estimate_population` controls
#'   alternating THETA/SIGMA updates. Ordinary covariance is not regular for a
#'   discrete support distribution; use bootstrap uncertainty for NPML/NPAG.
#' @export
nm_est <- function(model, data,
                   method = c("FOCEI", "FOCE", "FO", "LAPLACE", "ITS",
                              "GQ", "IMP", "SAEM", "BAYES", "HMC", "NUTS",
                              "NPML", "NPAG"),
                   maxit = 200L, eta_maxit = 100L, tolerance = 1e-6,
                   trace = 0L, print_every = 0L, n_cores = 1L,
                   optimizer_backend = c("auto", "native", "r"),
                   covariance = FALSE,
                   covariance_type = c("auto", "hessian", "opg", "sandwich", "r", "s"),
                   covariance_tolerance = 1e-8,
                   covariance_samples = NULL, covariance_seed = NULL,
                   initial_eta = NULL, collect_output = TRUE, ...) {
  estimation_started <- proc.time()[["elapsed"]]
  method <- match.arg(method)
  optimizer_backend <- match.arg(optimizer_backend)
  covariance_type <- match.arg(covariance_type)
  print_every <- as.integer(print_every)
  if (length(print_every) != 1L || is.na(print_every) || print_every < 0L) {
    .nm_stop("`print_every` must be a non-negative integer.")
  }
  if (length(covariance) != 1L || is.na(covariance)) {
    .nm_stop("`covariance` must be TRUE or FALSE.")
  }
  covariance <- isTRUE(covariance)
  if (covariance && method %in% c("BAYES", "HMC", "NUTS")) {
    .nm_stop(method, " reports posterior SDs and credible intervals automatically; a frequentist covariance step is not applicable.")
  }
  if (covariance && !method %in% c("FO", "FOCE", "FOCEI", "LAPLACE", "ITS", "GQ", "IMP", "SAEM")) {
    .nm_stop("Covariance is available for FO, FOCE, FOCEI, LAPLACE, ITS, GQ, IMP, and SAEM fits.")
  }
  if (!inherits(model, c("nm_model", "NMEngine"))) {
    .nm_stop("`model` must be an nm_model or NMEngine.")
  }
  if (missing(data)) .nm_stop("`data` is required.")
  model_definition <- if (inherits(model, "NMEngine")) model$model else model
  if (identical(model_definition$LIK_CONFIG$error, "likelihood") &&
      method %in% c("FO", "FOCE", "FOCEI")) {
    .nm_stop(
      method, " assumes a Gaussian residual linearization and cannot be used ",
      "with a user-defined likelihood. Use LAPLACE for NONMEM-like conditional ",
      "likelihood estimation, or ITS/GQ/IMP/SAEM/BAYES/HMC/NUTS/NPML/NPAG."
    )
  }
  context <- .nm_estimation_context(model, data, n_cores = n_cores, method = method)
  if (!is.null(initial_eta)) {
    initial_eta <- as.matrix(initial_eta)
    expected <- c(context$n_subjects, context$n_eta)
    if (!identical(dim(initial_eta), expected) || any(!is.finite(initial_eta))) {
      .nm_stop("`initial_eta` must be a finite ", expected[[1L]], " x ",
               expected[[2L]], " matrix for this dataset.")
    }
  }
  if (!is.null(context$parallel)) {
    on.exit(try(parallel::stopCluster(context$parallel$cluster), silent = TRUE),
            add = TRUE)
  }
  map <- .nm_outer_map(context$model)
  fit <- if (method == "FO") {
    .nm_est_fo(
      context, map, maxit, tolerance, trace, print_every, optimizer_backend
    )
  } else if (method %in% c("FOCE", "FOCEI", "LAPLACE")) {
    .nm_est_nested(
      context, map, method, maxit, eta_maxit, tolerance, trace, print_every,
      optimizer_backend, initial_eta = initial_eta
    )
  } else {
    .nm_est_stochastic(
      context, map, method, maxit = maxit, eta_maxit = eta_maxit,
      tolerance = tolerance, trace = trace, print_every = print_every,
      optimizer_backend = optimizer_backend, initial_eta = initial_eta, ...
    )
  }
  model_fit_seconds <- unname(proc.time()[["elapsed"]] - estimation_started)
  fit$diagnostics$eta_maxit <- as.integer(eta_maxit)
  fit$diagnostics$tolerance <- as.numeric(tolerance)
  covariance_seconds <- NA_real_
  if (covariance) {
    covariance_started <- proc.time()[["elapsed"]]
    attr(fit, ".estimation_context") <- context
    fit$covariance <- tryCatch(
      nm_cov_step(
        fit, type = covariance_type, tolerance = covariance_tolerance,
        samples = covariance_samples, seed = covariance_seed,
        eta_maxit = eta_maxit
      ),
      error = function(error) structure(list(
        status = "failed", type = covariance_type,
        error = conditionMessage(error)
      ), class = "nm_covariance_error")
    )
    attr(fit, ".estimation_context") <- NULL
    covariance_seconds <- unname(proc.time()[["elapsed"]] - covariance_started)
  }
  fit$timing <- list(
    model_fit_seconds = as.numeric(model_fit_seconds),
    covariance_seconds = as.numeric(covariance_seconds),
    total_seconds = as.numeric(proc.time()[["elapsed"]] - estimation_started)
  )
  if (isTRUE(collect_output) && length(fit$model$OUTPUT %||% character())) {
    fit$output <- .nm_fit_selected_outputs(fit)
  }
  fit
}

#' @export
print.nm_fit <- function(x, ...) {
  cat("LibeRation fit\n")
  cat("  method:", .nm_fit_method_label(x), " objective:", format(x$objective),
      " convergence:", x$convergence, "\n")
  invisible(x)
}
