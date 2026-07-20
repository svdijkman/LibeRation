.nm_supported_estimation_methods <- c(
  "FOCEI", "FOCE", "FO", "LAPLACE", "ITS", "GQ", "IMP", "SAEM",
  "BAYES", "HMC", "NUTS", "NPML", "NPAG"
)

.nm_fit_method_label <- function(fit) {
  as.character(fit$sequence_label %||% fit$method %||% "")
}

#' Define one step in a sequential estimation
#'
#' @param method LibeRation estimation method.
#' @param ... Controls passed to [nm_est()] for this step.
#' @param label Optional human-readable step label.
#' @return A serializable `nm_est_stage` specification.
#' @export
nm_est_stage <- function(method, ..., label = NULL) {
  method <- toupper(as.character(method))
  if (length(method) != 1L || is.na(method) ||
      !method %in% .nm_supported_estimation_methods) {
    .nm_stop("Unknown estimation-stage method: ", paste(method, collapse = ", "), ".")
  }
  arguments <- list(...)
  if (length(arguments) && is.null(names(arguments))) {
    .nm_stop("Estimation-stage controls must be named.")
  }
  forbidden <- intersect(names(arguments), c("model", "data", "method"))
  if (length(forbidden)) {
    .nm_stop("Do not specify ", paste(forbidden, collapse = ", "),
             " inside an estimation stage.")
  }
  structure(
    list(method = method, label = as.character(label %||% method),
         arguments = arguments),
    class = "nm_est_stage"
  )
}

.nm_normalize_estimation_stages <- function(stages) {
  if (is.character(stages)) stages <- lapply(stages, nm_est_stage)
  if (!is.list(stages) || !length(stages)) {
    .nm_stop("`stages` must contain at least one estimation stage.")
  }
  lapply(stages, function(stage) {
    if (inherits(stage, "nm_est_stage")) return(stage)
    if (is.character(stage) && length(stage) == 1L) return(nm_est_stage(stage))
    if (!is.list(stage) || is.null(stage$method)) {
      .nm_stop("Each estimation stage must be created by nm_est_stage() or contain `method`.")
    }
    arguments <- stage$arguments %||% stage[setdiff(names(stage), c("method", "label"))]
    do.call(nm_est_stage, c(list(method = stage$method), arguments,
                            list(label = stage$label %||% stage$method)))
  })
}

.nm_model_with_estimates <- function(model, fit) {
  arguments <- model[intersect(names(model), names(formals(nm_model)))]
  theta <- model$THETAS
  omega <- model$OMEGAS
  sigma <- model$SIGMAS
  if (length(fit$theta) != nrow(theta) || length(fit$omega) != nrow(omega) ||
      length(fit$sigma) != nrow(sigma)) {
    .nm_stop("An estimation stage changed the population parameter dimensions.")
  }
  theta$Value <- as.numeric(fit$theta)
  omega$Value <- as.numeric(fit$omega)
  sigma$Value <- as.numeric(fit$sigma)
  arguments$THETAS <- theta
  arguments$OMEGAS <- omega
  arguments$SIGMAS <- sigma
  updated <- do.call(nm_model, arguments)
  attr(updated, "name") <- attr(model, "name", exact = TRUE)
  control <- attr(model, "nonmem_control", exact = TRUE)
  if (!is.null(control)) attr(updated, "nonmem_control") <- control
  updated
}

.nm_est_stage_record <- function(fit, specification, index, initial) {
  record <- fit
  record$model <- NULL
  record$data <- NULL
  record$output <- NULL
  record$stages <- NULL
  record$stage <- list(
    index = as.integer(index), label = specification$label,
    method = specification$method, arguments = specification$arguments,
    initial = initial
  )
  class(record) <- c("nm_est_stage_result", setdiff(class(record), "nm_est_stage_result"))
  record
}

#' Run multiple estimation methods sequentially
#'
#' Population estimates from each completed step become the initial THETA,
#' OMEGA, and SIGMA values for the next. Compatible conditional ETA estimates
#' are used as a warm start for FOCE/FOCEI/LAPLACE and SAEM. The returned object
#' remains an `nm_fit` for compatibility and contains compact per-step results
#' in `stages` plus `method_sequence` and `sequence_label` metadata.
#'
#' @param model An [nm_model()] or compiled [NMEngine()].
#' @param data NONMEM-style event data.
#' @param stages Ordered character vector or list of [nm_est_stage()] objects.
#' @param on_error Stop with an `nm_est_sequence_error`, or return the last
#'   successful fit with failure metadata.
#' @return The final `nm_fit`, augmented with sequential-estimation metadata.
#' @export
nm_est_sequence <- function(model, data, stages,
                            on_error = c("stop", "return")) {
  on_error <- match.arg(on_error)
  specifications <- .nm_normalize_estimation_stages(stages)
  current_model <- if (inherits(model, "NMEngine")) model$model else model
  if (!inherits(current_model, "nm_model")) {
    .nm_stop("`model` must be an nm_model or NMEngine.")
  }
  started <- proc.time()[["elapsed"]]
  records <- list()
  previous <- NULL
  for (index in seq_along(specifications)) {
    specification <- specifications[[index]]
    arguments <- specification$arguments
    arguments$collect_output <- index == length(specifications)
    warm_methods <- c("FOCE", "FOCEI", "LAPLACE", "SAEM")
    if (!is.null(previous) && specification$method %in% warm_methods &&
        is.null(arguments$initial_eta) &&
        identical(dim(previous$eta), c(length(unique(.nm_engine_data(current_model, data)$.ID_INDEX)),
                                        .nm_eta_columns(current_model, .nm_engine_data(current_model, data))))) {
      arguments$initial_eta <- previous$eta
    }
    initial <- list(
      theta = current_model$THETAS$Value,
      omega = current_model$OMEGAS$Value,
      sigma = current_model$SIGMAS$Value,
      eta_warm_start = !is.null(arguments$initial_eta)
    )
    fit <- tryCatch(
      do.call(nm_est, c(list(
        model = current_model, data = data, method = specification$method
      ), arguments)),
      error = identity
    )
    if (inherits(fit, "error")) {
      failure <- list(
        index = as.integer(index), method = specification$method,
        label = specification$label, message = conditionMessage(fit)
      )
      if (on_error == "return" && !is.null(previous)) {
        if (length(previous$model$OUTPUT %||% character())) {
          previous$output <- .nm_fit_selected_outputs(previous)
        }
        previous$stages <- records
        previous$sequence_failure <- failure
        previous$method_sequence <- vapply(specifications, `[[`, character(1), "method")
        previous$sequence_label <- paste(previous$method_sequence, collapse = " -> ")
        previous$sequence_complete <- FALSE
        return(previous)
      }
      condition <- structure(
        list(message = paste0("Estimation stage ", index, " (", specification$method,
                              ") failed: ", conditionMessage(fit)),
             call = NULL, stage = failure, completed = records),
        class = c("nm_est_sequence_error", "error", "condition")
      )
      stop(condition)
    }
    records[[index]] <- .nm_est_stage_record(fit, specification, index, initial)
    previous <- fit
    if (index < length(specifications)) {
      current_model <- .nm_model_with_estimates(current_model, fit)
    }
  }
  methods <- vapply(specifications, `[[`, character(1), "method")
  previous$stages <- records
  previous$method_sequence <- methods
  previous$sequence_label <- paste(methods, collapse = " -> ")
  previous$sequence_complete <- TRUE
  previous$timing$stage_total_seconds <- vapply(
    records, function(stage) as.numeric(stage$timing$total_seconds %||% NA_real_), numeric(1)
  )
  previous$timing$sequence_total_seconds <- as.numeric(proc.time()[["elapsed"]] - started)
  previous
}
