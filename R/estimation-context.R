# Subject-data views, tape pools, and estimation-context construction.
# Split from estimation.R as a behaviour-preserving source move.

.nm_parallel_registry <- new.env(parent = emptyenv())

.nm_parallel_worker_state <- function() {
  state <- .nm_parallel_registry$state
  if (!is.environment(state)) {
    .nm_stop("Parallel estimation worker state is not initialized.")
  }
  state
}

.nm_subject_data <- function(data, subject) {
  out <- as.data.frame(data[data$.ID_INDEX == subject, , drop = FALSE])
  internal <- intersect(c(".ID_INDEX", ".source_row", ".generated", ".sort_priority"), names(out))
  out[internal] <- NULL
  nm_dataset(out)
}

.nm_subject_store <- function(data) {
  index <- as.integer(data$.ID_INDEX)
  subject_ids <- unique(index)
  if (!length(subject_ids) || anyNA(subject_ids)) {
    .nm_stop("A shared subject-data store requires valid subject indices.")
  }
  starts <- c(1L, which(diff(index) != 0L) + 1L)
  lengths <- diff(c(starts, length(index) + 1L))
  if (!identical(index[starts], subject_ids) ||
      anyDuplicated(index[starts]) || length(starts) != length(subject_ids)) {
    .nm_stop("Subject rows must be contiguous in the normalized dataset.")
  }
  store <- new.env(parent = emptyenv())
  store$data <- data
  store$subject_ids <- subject_ids
  store$starts <- unname(starts)
  store$lengths <- unname(lengths)
  store$pointer <- .liberation_subject_store_create(
    data, as.integer(store$starts), as.integer(store$lengths)
  )
  store$metadata <- new.env(parent = emptyenv())
  store$metadata$columns <- names(data)
  store$materializations <- 0L
  store$projections <- 0L
  class(store) <- c("nm_subject_store", "environment")
  store
}

.nm_subject_store_view <- function(store, subject) {
  if (!inherits(store, "nm_subject_store")) {
    .nm_stop("Invalid shared subject-data store.")
  }
  subject <- as.integer(subject)
  if (length(subject) != 1L || is.na(subject) || subject < 1L ||
      subject > length(store$subject_ids)) {
    .nm_stop("Subject view is outside the shared data store.")
  }
  structure(
    list(
      pointer = store$pointer, subject = subject,
      rows = as.integer(store$lengths[[subject]]), metadata = store$metadata
    ),
    class = "nm_subject_view"
  )
}

.nm_subject_store_data <- function(store, subject) {
  if (!inherits(store, "nm_subject_store")) {
    .nm_stop("Invalid shared subject-data store.")
  }
  subject <- as.integer(subject)
  if (length(subject) != 1L || is.na(subject) || subject < 1L ||
      subject > length(store$subject_ids)) {
    .nm_stop("Subject view is outside the shared data store.")
  }
  store$materializations <- store$materializations + 1L
  rows <- seq.int(store$starts[[subject]], length.out = store$lengths[[subject]])
  out <- as.data.frame(store$data[rows, , drop = FALSE])
  internal <- intersect(
    c(".ID_INDEX", ".source_row", ".generated", ".sort_priority"), names(out)
  )
  out[internal] <- NULL
  nm_dataset(out)
}

.nm_subject_view_project <- function(view, columns, observed_only = FALSE,
                                     first_only = FALSE) {
  store <- attr(view, "store", exact = TRUE)
  if (inherits(store, "nm_subject_store")) {
    store$projections <- store$projections + 1L
  }
  .liberation_subject_view_project(
    view, unique(as.character(columns)), isTRUE(observed_only),
    isTRUE(first_only)
  )
}

.nm_subject_store_chunk <- function(store, subjects) {
  subjects <- as.integer(subjects)
  if (!length(subjects)) return(NULL)
  first <- store$starts[[subjects[[1L]]]]
  last_subject <- tail(subjects, 1L)
  last <- store$starts[[last_subject]] + store$lengths[[last_subject]] - 1L
  raw <- as.data.frame(store$data[seq.int(first, last), , drop = FALSE])
  raw$.ID_INDEX <- match(raw$.ID_INDEX, unique(raw$.ID_INDEX))
  rownames(raw) <- NULL
  attr(raw, "addl_materialized") <- TRUE
  class(raw) <- unique(c("nm_dataset", class(raw)))
  raw
}

.nm_context_observation_counts <- function(context) {
  source <- if (!is.null(context$subject_store)) {
    context$subject_store$pointer
  } else {
    lapply(context$subjects, function(evaluator) evaluator$data_input())
  }
  as.integer(.liberation_subject_observation_counts(source))
}

.nm_subject_dynamic_input <- function(evaluator, interaction = TRUE) {
  value <- if (isTRUE(interaction)) evaluator$objective_dynamic else
    evaluator$noninteraction_dynamic
  if (!is.null(value)) return(value)
  value <- evaluator$data_input
  if (is.function(value)) value <- value()
  if (is.null(value)) {
    .nm_stop("Subject evaluator does not expose dynamic input data.")
  }
  value
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
  repair <- nm_covariance_repair(
    matrix, method = "jitter", tolerance = 1e-9,
    preserve_diagonal = FALSE
  )
  jitter <- repair$diagnostics$diagonal_shift
  if (jitter > largest * 1e-2) {
    .nm_stop(context, " is not sufficiently positive definite.")
  }
  adjusted <- repair$matrix
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
  data_names <- if (inherits(data, "nm_subject_view")) {
    data$metadata$columns
  } else names(data)
  intersect(inputs[!parameter & !inputs %in% structural], data_names)
}

.nm_structure_key <- function(value) {
  digest::digest(value, algo = "sha256", serialize = TRUE)
}

.nm_prediction_structure <- function(model, data) {
  dynamic <- .nm_prediction_dynamic_columns(model, data)
  ignored <- c(
    "DV", "MDV", "LLOQ", "BLQ", "CENS", ".source_row", ".generated",
    ".sort_priority", ".ID_INDEX", dynamic
  )
  source <- paste(model$PRED %||% "", model$DES %||% "")
  if (!grepl("\\bID\\b", source, perl = TRUE)) ignored <- c(ignored, "ID")
  if (inherits(data, "nm_subject_view")) {
    key <- .liberation_subject_view_signature(
      data, unique(ignored), include_fo_layout = FALSE
    )
    return(list(
      key = key,
      value = list(
        signature = key, dynamic_columns = dynamic,
        rows = as.integer(data$rows), storage = "native-row-view"
      )
    ))
  }
  columns <- setdiff(names(data), ignored)
  structural <- list(
    dynamic_columns = dynamic,
    structural_data = as.data.frame(data[columns], stringsAsFactors = FALSE),
    rows = nrow(data)
  )
  key <- .nm_structure_key(structural)
  list(
    key = key,
    value = list(signature = key, dynamic_columns = dynamic, rows = nrow(data))
  )
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
  if (inherits(data, "nm_subject_view")) {
    dynamic <- .nm_prediction_dynamic_columns(model, data)
    ignored <- c(
      "DV", "MDV", "LLOQ", "BLQ", "CENS", ".source_row", ".generated",
      ".sort_priority", ".ID_INDEX", dynamic
    )
    source <- paste(model$PRED %||% "", model$DES %||% "")
    if (!grepl("\\bID\\b", source, perl = TRUE)) ignored <- c(ignored, "ID")
    key <- .liberation_subject_view_signature(
      data, unique(ignored), include_fo_layout = TRUE
    )
    return(list(
      key = key,
      value = list(
        signature = key, prediction_signature = prediction$key,
        rows = as.integer(data$rows), storage = "native-row-view"
      )
    ))
  }
  observed <- data$EVID == 0L & data$MDV == 0L & is.finite(data$DV)
  value <- list(
    prediction_signature = prediction$key,
    observed = as.logical(observed),
    dvid = if ("DVID" %in% names(data)) as.integer(data$DVID) else rep.int(1L, nrow(data))
  )
  list(key = .nm_structure_key(value), value = value)
}

.nm_fo_pool_tape <- function(pool, evaluator, theta, sigma, omega) {
  data <- evaluator$data_input()
  structure <- .nm_fo_structure(evaluator$engine$model, data)
  bucket <- pool[[structure$key]] %||% list()
  if (length(bucket)) {
    for (entry in bucket) {
      if (identical(entry$structure, structure$value)) {
        evaluator$fo_tape <- entry$tape
        evaluator$fo_dynamic <- .liberation_fo_tape_new_dynamic(
          entry$tape$pointer, data
        )
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

.nm_shared_fo_objective_eligible <- function(model) {
  likelihood <- model$LIK_CONFIG %||% list()
  !isTRUE(model$USE_ODE) &&
    !identical(likelihood$error, "likelihood") &&
    identical(likelihood$sigma_corr %||% "independent", "independent") &&
    identical(likelihood$blq_method %||% "none", "none") &&
    identical(as.integer(likelihood$iov %||% 0L), 0L) &&
    !length(likelihood$residual_groups %||% list()) &&
    !length(model$MIXTURE %||% list()) &&
    is.null(model$RE_CONFIG)
}

.nm_objective_pool_tape <- function(pool, engine, data, theta, sigma, omega,
                                    n_eta, prediction_tape) {
  structure <- .nm_fo_structure(engine$model, data)
  bucket <- pool[[structure$key]] %||% list()
  if (length(bucket)) {
    for (entry in bucket) {
      if (identical(entry$structure, structure$value)) {
        return(list(tape = entry$tape, recorded = FALSE))
      }
    }
  }
  tape <- list(pointer = .liberation_shared_fo_objective_tape_create(
    engine$pointer, prediction_tape$pointer, data, as.numeric(theta),
    matrix(0, 1L, n_eta), as.numeric(sigma), as.numeric(omega)
  ))
  pool[[structure$key]] <- c(
    bucket, list(list(structure = structure$value, tape = tape))
  )
  list(tape = tape, recorded = TRUE)
}

.NMSubjectEvaluator <- R6::R6Class(
  ".NMSubjectEvaluator",
  public = list(
    engine = NULL,
    data = NULL,
    data_store = NULL,
    data_subject = NA_integer_,
    objective_tape = NULL,
    noninteraction_tape = NULL,
    prediction_tape = NULL,
    fo_tape = NULL,
    curvature_tapes = NULL,
    prediction_dynamic = NULL,
    objective_dynamic = NULL,
    noninteraction_dynamic = NULL,
    fo_dynamic = NULL,
    tape_anchor = NULL,
    tape_signature = NULL,
    tape_records = 0L,
    tape_retapes = 0L,
    tape_checks = 0L,
    data_materializations = 0L,
    data_projections = 0L,
    tape_profile = "full",
    n_theta = 0L,
    n_eta = 0L,
    n_sigma = 0L,

    initialize = function(engine, data, theta, sigma, omega, n_eta = NULL,
                          prediction_tape = NULL,
                          objective_tape = NULL,
                          objective_tape_recorded = TRUE,
                          data_store = NULL, data_subject = NA_integer_,
                          tape_profile = c("full", "fo")) {
      self$tape_profile <- match.arg(tape_profile)
      self$engine <- engine
      normalized_data <- if (inherits(data, "nm_subject_view")) data else
        .nm_engine_data(engine$model, data)
      self$data_store <- data_store
      self$data_subject <- as.integer(data_subject)
      self$data <- if (is.null(data_store)) normalized_data else NULL
      self$n_theta <- length(theta)
      self$n_eta <- as.integer(n_eta %||% .nm_eta_columns(engine$model, normalized_data))
      self$n_sigma <- length(sigma)
      self$tape_signature <- .nm_prediction_structure(engine$model, normalized_data)$key
      self$record_tapes(
        theta, sigma, omega, rep(0, self$n_eta),
        prediction_tape = prediction_tape, objective_tape = objective_tape,
        objective_tape_recorded = objective_tape_recorded,
        recording_data = normalized_data
      )
    },

    data_frame = function() {
      if (!is.null(self$data)) return(self$data)
      self$data_materializations <- self$data_materializations + 1L
      .nm_subject_store_data(self$data_store, self$data_subject)
    },

    data_input = function() {
      if (!is.null(self$data)) return(self$data)
      .nm_subject_store_view(self$data_store, self$data_subject)
    },

    project = function(columns, observed_only = FALSE, first_only = FALSE) {
      columns <- unique(as.character(columns))
      self$data_projections <- self$data_projections + 1L
      if (!is.null(self$data)) {
        rows <- seq_len(nrow(self$data))
        if (isTRUE(observed_only)) {
          rows <- rows[
            self$data$EVID[rows] == 0L & self$data$MDV[rows] == 0L &
              is.finite(self$data$DV[rows])
          ]
        }
        if (isTRUE(first_only) && length(rows)) rows <- rows[[1L]]
        result <- lapply(columns, function(column) {
          if (column %in% names(self$data)) self$data[[column]][rows] else NULL
        })
        names(result) <- columns
        attr(result, "rows") <- as.integer(rows - 1L)
        return(result)
      }
      .nm_subject_view_project(
        self$data_input(), columns, observed_only = observed_only,
        first_only = first_only
      )
    },

    observation_data = function(columns = character()) {
      self$project(
        unique(c("DV", "DVID", "TIME", ".ID_INDEX", columns)),
        observed_only = TRUE
      )
    },

    observation_count = function() {
      .liberation_subject_observation_count(self$data_input())
    },

    record_tapes = function(theta, sigma, omega, eta = rep(0, self$n_eta),
                            retape = FALSE, prediction_tape = NULL,
                            objective_tape = NULL,
                            objective_tape_recorded = TRUE,
                            recording_data = NULL) {
      data <- recording_data %||% self$data_input()
      self$curvature_tapes <- list()
      self$fo_tape <- NULL
      self$fo_dynamic <- NULL
      eta <- matrix(as.numeric(eta), 1L, self$n_eta)
      self$prediction_tape <- prediction_tape %||% self$engine$prediction_tape(
        data, theta = theta, eta = eta, sigma = sigma
      )
      self$prediction_dynamic <- .liberation_prediction_tape_new_dynamic(
        self$prediction_tape$pointer, data
      )
      self$objective_tape <- objective_tape %||% self$engine$objective_tape(
        data, theta = theta, eta = eta, sigma = sigma, omega = omega
      )
      self$objective_dynamic <- .liberation_objective_tape_new_dynamic(
        self$objective_tape$pointer, data
      )
      self$noninteraction_tape <- if (identical(self$tape_profile, "fo")) NULL else
        self$engine$objective_tape(
          data, theta = theta, eta = eta, sigma = sigma, omega = omega,
          interaction = FALSE
        )
      self$noninteraction_dynamic <- if (is.null(self$noninteraction_tape)) NULL else
        .liberation_objective_tape_new_dynamic(
          self$noninteraction_tape$pointer, data
        )
      self$tape_anchor <- c(theta, as.numeric(eta), sigma, omega)
      self$tape_records <- self$tape_records +
        as.integer(is.null(objective_tape) || isTRUE(objective_tape_recorded))
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

    activate_objective_tape = function(tape = self$objective_tape,
                                       dynamic = self$objective_dynamic) {
      if (!is.null(tape)) {
        if (is.null(dynamic)) {
          dynamic <- .liberation_objective_tape_new_dynamic(
            tape$pointer, self$data_input()
          )
        } else {
          .liberation_objective_tape_set_dynamic(tape$pointer, dynamic)
        }
      }
      invisible(tape)
    },

    prediction_point = function(theta, eta, sigma) c(theta, eta, sigma),

    fo_point = function(theta, sigma, omega) c(theta, sigma, omega),

    ensure_fo_tape = function(theta, sigma, omega) {
      if (isTRUE(self$engine$model$USE_ODE)) {
        self$ensure_valid_tapes(theta, sigma, omega)
      }
      if (is.null(self$fo_tape)) {
        data <- self$data_input()
        .liberation_prediction_tape_set_dynamic(
          self$prediction_tape$pointer, self$prediction_dynamic
        )
        self$fo_tape <- list(
          pointer = .liberation_fo_tape_create(
            self$engine$pointer, self$prediction_tape$pointer, data,
            as.numeric(theta), as.numeric(sigma), as.numeric(omega),
            isTRUE(getOption("LibeRation.fo_low_rank", TRUE)),
            as.numeric(getOption("LibeRation.fo_low_rank_tolerance", 1e-9)),
            as.numeric(getOption(
              "LibeRation.fo_low_rank_condition_tolerance", 1e-12
            ))
          )
        )
        self$fo_dynamic <- .liberation_fo_tape_new_dynamic(
          self$fo_tape$pointer, data
        )
      }
      .liberation_objective_tape_set_dynamic(
        self$fo_tape$pointer, self$fo_dynamic
      )
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
        data <- self$data_input()
        .liberation_prediction_tape_set_dynamic(
          self$prediction_tape$pointer, self$prediction_dynamic
        )
        self$curvature_tapes[[approximation]] <- list(
          pointer = .liberation_curvature_tape_create(
            self$engine$pointer, self$prediction_tape$pointer,
            self$objective_tape$pointer, data, as.numeric(theta),
            as.numeric(eta), as.numeric(sigma), as.numeric(omega),
            approximation
          ),
          anchor = as.numeric(eta)
        )
        self$curvature_tapes[[approximation]]$dynamic <-
          .liberation_objective_tape_new_dynamic(
            self$curvature_tapes[[approximation]]$pointer, data
          )
      }
      invisible(self$curvature_tapes[[approximation]])
    },

    curvature = function(theta, eta, sigma, omega, approximation,
                         gradient = TRUE) {
      private$with_retape(
        function() {
          self$ensure_curvature_tape(theta, eta, sigma, omega, approximation)
          self$activate_objective_tape(
            self$curvature_tapes[[approximation]],
            self$curvature_tapes[[approximation]]$dynamic
          )
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
          self$activate_objective_tape(
            if (isTRUE(interaction)) self$objective_tape else self$noninteraction_tape,
            if (isTRUE(interaction)) self$objective_dynamic else
              self$noninteraction_dynamic
          )
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
        # Adaptive ODE trajectories may differ between proposal points.  A
        # single batched tape cannot safely represent all such paths, so let
        # the scalar evaluator validate/retape at every ETA point.
        return(vapply(seq_len(nrow(eta)), function(row) {
          self$objective(
            theta, eta[row, ], sigma, omega, gradient = FALSE,
            interaction = interaction
          )$value
        }, numeric(1)))
      }
      private$with_retape(
        function() {
          tape <- if (isTRUE(interaction)) self$objective_tape else self$noninteraction_tape
          self$activate_objective_tape(
            tape, if (isTRUE(interaction)) self$objective_dynamic else
              self$noninteraction_dynamic
          )
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
        evaluated <- lapply(seq_len(nrow(eta)), function(row) {
          self$objective(
            theta, eta[row, ], sigma, omega, gradient = TRUE,
            interaction = interaction
          )
        })
        return(list(
          value = vapply(evaluated, `[[`, numeric(1), "value"),
          gradient = do.call(rbind, lapply(evaluated, `[[`, "gradient"))
        ))
      }
      points <- cbind(
        matrix(theta, nrow(eta), length(theta), byrow = TRUE), eta,
        matrix(sigma, nrow(eta), length(sigma), byrow = TRUE),
        matrix(omega, nrow(eta), length(omega), byrow = TRUE)
      )
      private$with_retape(
        function() {
          tape <- if (isTRUE(interaction)) self$objective_tape else self$noninteraction_tape
          self$activate_objective_tape(
            tape, if (isTRUE(interaction)) self$objective_dynamic else
              self$noninteraction_dynamic
          )
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
          self$activate_objective_tape(
            tape, if (isTRUE(interaction)) self$objective_dynamic else
              self$noninteraction_dynamic
          )
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
          .liberation_prediction_tape_set_dynamic(
            self$prediction_tape$pointer, self$prediction_dynamic
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
                         interaction = TRUE, exact_hessian = TRUE,
                         prior_mean = NULL, base_precision = NULL,
                         prior_precision = NULL) {
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
      self$activate_objective_tape(
        tape, if (isTRUE(interaction)) self$objective_dynamic else
          self$noninteraction_dynamic
      )
      custom_prior <- !is.null(prior_mean) || !is.null(base_precision) ||
        !is.null(prior_precision)
      if (custom_prior &&
          (is.null(prior_mean) || is.null(base_precision) ||
           is.null(prior_precision))) {
        .nm_stop(
          "Custom ETA mode fitting requires prior mean, population precision, ",
          "and replacement prior precision."
        )
      }
      native <- tryCatch(
        if (custom_prior) {
          .liberation_objective_tape_eta_mode_prior(
            tape$pointer,
            self$objective_point(theta, start, sigma, omega),
            eta_positions, as.numeric(start), as.numeric(prior_mean),
            as.matrix(base_precision), as.matrix(prior_precision),
            as.integer(maxit), as.numeric(tolerance), isTRUE(exact_hessian)
          )
        } else {
          .liberation_objective_tape_eta_mode(
            tape$pointer,
            self$objective_point(theta, start, sigma, omega),
            eta_positions, as.numeric(start), as.integer(maxit),
            as.numeric(tolerance), isTRUE(exact_hessian)
          )
        },
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
        if (custom_prior && is.finite(value)) {
          centered <- eta - prior_mean
          value <- value - drop(crossprod(eta, base_precision %*% eta)) +
            drop(crossprod(centered, prior_precision %*% centered))
        }
        if (is.finite(value)) value else .Machine$double.xmax / 1e100
      }
      gr <- function(eta) {
        result <- self$objective(theta, eta, sigma, omega, gradient = TRUE,
                                 interaction = interaction)
        value <- unname(result$gradient[eta_positions])
        if (custom_prior) {
          value <- value - 2 * as.numeric(base_precision %*% eta) +
            2 * as.numeric(prior_precision %*% (eta - prior_mean))
        }
        value
      }
      fit <- stats::optim(
        as.numeric(start), fn, gr, method = "BFGS",
        control = list(maxit = as.integer(maxit), reltol = tolerance)
      )
      at_mode <- self$objective(
        theta, fit$par, sigma, omega, gradient = TRUE,
        hessian = isTRUE(exact_hessian), interaction = interaction
      )
      if (custom_prior) {
        centered <- fit$par - prior_mean
        at_mode$value <- at_mode$value -
          drop(crossprod(fit$par, base_precision %*% fit$par)) +
          drop(crossprod(centered, prior_precision %*% centered))
        at_mode$gradient[eta_positions] <-
          at_mode$gradient[eta_positions] -
          2 * as.numeric(base_precision %*% fit$par) +
          2 * as.numeric(prior_precision %*% centered)
        if (isTRUE(exact_hessian)) {
          at_mode$hessian[eta_positions, eta_positions] <-
            at_mode$hessian[eta_positions, eta_positions, drop = FALSE] -
            2 * base_precision + 2 * prior_precision
        }
      }
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
  .liberation_objective_tape_collection_values(
    tapes, points, lapply(
      evaluators, .nm_subject_dynamic_input, interaction = interaction
    )
  )
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
  .liberation_objective_tape_collection_gradients(
    tapes, points, lapply(
      evaluators, .nm_subject_dynamic_input, interaction = interaction
    )
  )
}

.nm_objective_collection_value_gradient <- function(
    evaluators, parameters, eta, interaction = TRUE) {
  eta <- as.matrix(eta)
  if (!length(evaluators)) {
    return(list(value = numeric(), gradient = matrix(numeric(), 0L, 0L)))
  }
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
  .liberation_objective_tape_collection_value_gradients(
    tapes, points, lapply(
      evaluators, .nm_subject_dynamic_input, interaction = interaction
    )
  )
}

.nm_estimation_context_build <- function(model, data, n_cores = 1L, method = NULL) {
  engine <- if (inherits(model, "NMEngine")) model else nm_compile(model)
  model <- engine$model
  data <- .nm_engine_data(model, data)
  if (!is.null(model$MU) && nrow(model$MU)) {
    declared <- unique(unlist(strsplit(
      paste(model$MU$COVARIATES, collapse = ";"), "\\s*[;,]\\s*", perl = TRUE
    )))
    declared <- declared[nzchar(declared)]
    for (covariate in declared) {
      if (!covariate %in% names(data)) {
        .nm_stop("MU covariate `", covariate, "` is absent from the dataset.")
      }
      varying <- vapply(split(data[[covariate]], data$.ID_INDEX), function(value) {
        length(unique(value[!is.na(value)])) > 1L
      }, logical(1))
      if (any(varying)) {
        .nm_stop(
          "MU covariate `", covariate,
          "` varies within subject; MU references require subject-level covariates."
        )
      }
    }
  }
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
  use_subject_views <- isTRUE(getOption("LibeRation.subject_data_views", TRUE))
  subject_store <- if (use_subject_views) .nm_subject_store(data) else NULL
  subject_data <- if (use_subject_views) NULL else
    lapply(seq_len(n_subjects), function(index) .nm_subject_data(data, index))
  subject_input <- function(index) {
    if (use_subject_views) .nm_subject_store_view(subject_store, index) else
      subject_data[[index]]
  }
  prediction_pool <- new.env(parent = emptyenv())
  objective_pool <- new.env(parent = emptyenv())
  shared_fo_objective <- identical(tape_profile, "fo") &&
    isTRUE(getOption("LibeRation.fo_objective_tape_sharing", TRUE)) &&
    .nm_shared_fo_objective_eligible(model)
  subjects <- lapply(seq_len(n_subjects), function(index) {
    current_data <- subject_input(index)
    prediction_tape <- .nm_prediction_pool_tape(
      prediction_pool, engine, current_data, model$THETAS$Value,
      model$SIGMAS$Value, expanded_n_eta
    )
    objective <- if (shared_fo_objective) {
      .nm_objective_pool_tape(
        objective_pool, engine, current_data, model$THETAS$Value,
        model$SIGMAS$Value, model$OMEGAS$Value, expanded_n_eta,
        prediction_tape
      )
    } else list(tape = NULL, recorded = TRUE)
    .NMSubjectEvaluator$new(
      engine, current_data, model$THETAS$Value,
      model$SIGMAS$Value, model$OMEGAS$Value, n_eta = expanded_n_eta,
      prediction_tape = prediction_tape, objective_tape = objective$tape,
      objective_tape_recorded = objective$recorded,
      data_store = subject_store,
      data_subject = if (use_subject_views) index else NA_integer_,
      tape_profile = tape_profile
    )
  })
  native_subject_threads <- 1L
  stochastic_method <- toupper(method %||% "")
  fused_native_thread_eligible <- n_cores > 1L &&
    isTRUE(getOption("LibeRation.native_subject_threads", TRUE)) &&
    identical(model$NUMERICAL_MODE %||% "nonmem_compatibility", "liber_optimized") &&
    stochastic_method %in% c("SAEM", "BAYES") &&
    all(vapply(subjects, function(evaluator) {
      startsWith(
        evaluator$prediction_tape$propagation_kernel %||% "",
        "specialized-advan"
      )
    }, logical(1)))
  cppad_native_thread_eligible <- n_cores > 1L &&
    isTRUE(getOption("LibeRation.cppad_subject_threads", TRUE)) &&
    identical(model$NUMERICAL_MODE %||% "nonmem_compatibility", "liber_optimized") &&
    stochastic_method %in% c("ITS", "IMP")
  native_thread_eligible <- fused_native_thread_eligible ||
    cppad_native_thread_eligible
  if (native_thread_eligible) native_subject_threads <- n_cores
  parallel_state <- NULL
  if (n_cores > 1L && native_subject_threads == 1L) {
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
                 n_eta, tape_profile, library_paths, use_subject_views) {
          .libPaths(unique(c(library_paths, .libPaths())))
          namespace <- asNamespace("LibeRation")
          compiler <- get("nm_compile", envir = namespace)
          evaluator_class <- get(".NMSubjectEvaluator", envir = namespace)
          prediction_pool_tape <- get(".nm_prediction_pool_tape", envir = namespace)
          objective_pool_tape <- get(".nm_objective_pool_tape", envir = namespace)
          subject_store_create <- get(".nm_subject_store", envir = namespace)
          subject_store_view <- get(".nm_subject_store_view", envir = namespace)
          shared_fo_eligible <- get(
            ".nm_shared_fo_objective_eligible", envir = namespace
          )
          compiled <- compiler(specification)
          prediction_pool <- new.env(parent = emptyenv())
          objective_pool <- new.env(parent = emptyenv())
          share_objective <- identical(tape_profile, "fo") &&
            shared_fo_eligible(specification)
          worker_store <- if (use_subject_views) {
            subject_store_create(data_chunks[[index]])
          } else NULL
          subject_count <- if (use_subject_views) {
            length(worker_store$subject_ids)
          } else length(data_chunks[[index]])
          evaluators <- lapply(seq_len(subject_count), function(subject) {
            subject_data <- if (use_subject_views) {
              subject_store_view(worker_store, subject)
            } else data_chunks[[index]][[subject]]
            prediction_tape <- prediction_pool_tape(
              prediction_pool, compiled, subject_data, theta, sigma, n_eta
            )
            objective <- if (share_objective) {
              objective_pool_tape(
                objective_pool, compiled, subject_data, theta, sigma, omega,
                n_eta, prediction_tape
              )
            } else list(tape = NULL, recorded = TRUE)
            evaluator_class$new(
              compiled, subject_data, theta, sigma, omega, n_eta = n_eta,
              prediction_tape = prediction_tape,
              objective_tape = objective$tape,
              objective_tape_recorded = objective$recorded,
              data_store = worker_store,
              data_subject = if (use_subject_views) subject else NA_integer_,
              tape_profile = tape_profile
            )
          })
          state <- new.env(parent = emptyenv())
          state$subjects <- evaluators
          state$model <- specification
          state$subject_store <- worker_store
          registry <- get(".nm_parallel_registry", envir = namespace)
          registry$state <- state
          TRUE
        },
        data_chunks = if (use_subject_views) {
          lapply(chunks, function(rows) .nm_subject_store_chunk(subject_store, rows))
        } else {
          lapply(chunks, function(rows) subject_data[rows])
        },
        specification = model, theta = model$THETAS$Value,
        sigma = model$SIGMAS$Value, omega = model$OMEGAS$Value,
        n_eta = expanded_n_eta, tape_profile = tape_profile,
        use_subject_views = use_subject_views,
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
       subject_store = subject_store,
       subject_data_layout = if (use_subject_views) "native-row-view" else "copied",
       n_subjects = n_subjects, n_eta = expanded_n_eta,
       parallel = parallel_state,
       native_subject_threads = native_subject_threads,
       shared_fo_objective = shared_fo_objective)
}

.nm_estimation_context_cache <- new.env(parent = emptyenv())
.nm_estimation_context_cache$.order <- character()

.nm_estimation_context_cache_clear <- function() {
  keys <- setdiff(ls(.nm_estimation_context_cache, all.names = TRUE), ".order")
  if (length(keys)) rm(list = keys, envir = .nm_estimation_context_cache)
  .nm_estimation_context_cache$.order <- character()
  invisible(NULL)
}

.nm_estimation_context <- function(model, data, n_cores = 1L, method = NULL) {
  definition <- if (inherits(model, "NMEngine")) model$model else model
  cacheable <- identical(toupper(method %||% ""), "FO") &&
    identical(as.integer(n_cores), 1L) &&
    isTRUE(getOption("LibeRation.fo_context_cache", TRUE)) &&
    .nm_shared_fo_objective_eligible(definition)
  if (!cacheable) {
    context <- .nm_estimation_context_build(model, data, n_cores, method)
    context$cache_hit <- FALSE
    return(context)
  }
  signature <- list(
    version = 2L, method = "FO", model = definition,
    data = as.data.frame(data, stringsAsFactors = FALSE),
    tape_options = list(
      sharing = getOption("LibeRation.fo_objective_tape_sharing", TRUE),
      low_rank = getOption("LibeRation.fo_low_rank", TRUE),
      low_rank_tolerance = getOption(
        "LibeRation.fo_low_rank_tolerance", 1e-9
      ),
      low_rank_condition_tolerance = getOption(
        "LibeRation.fo_low_rank_condition_tolerance", 1e-12
      ),
      tape_guard_radius = getOption("LibeRation.tape_guard_radius", 0.5),
      subject_data_views = getOption("LibeRation.subject_data_views", TRUE)
    )
  )
  key <- .nm_structure_key(signature)
  if (exists(key, envir = .nm_estimation_context_cache, inherits = FALSE)) {
    entry <- get(key, envir = .nm_estimation_context_cache, inherits = FALSE)
    if (identical(entry$signature, signature)) {
      context <- entry$context
      .nm_estimation_context_cache$.order <- c(
        setdiff(.nm_estimation_context_cache$.order, key), key
      )
      context$cache_hit <- TRUE
      return(context)
    }
  }
  context <- .nm_estimation_context_build(model, data, n_cores, method)
  context$cache_hit <- FALSE
  assign(key, list(signature = signature, context = context),
         envir = .nm_estimation_context_cache)
  .nm_estimation_context_cache$.order <- c(
    .nm_estimation_context_cache$.order, key
  )
  limit <- max(1L, as.integer(getOption("LibeRation.fo_context_cache_size", 4L)))
  while (length(.nm_estimation_context_cache$.order) > limit) {
    drop <- .nm_estimation_context_cache$.order[[1L]]
    rm(list = drop, envir = .nm_estimation_context_cache)
    .nm_estimation_context_cache$.order <-
      .nm_estimation_context_cache$.order[-1L]
  }
  context
}

