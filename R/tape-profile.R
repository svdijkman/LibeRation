.nm_profile_elapsed <- function(expression) {
  started <- proc.time()[["elapsed"]]
  value <- force(expression)
  list(value = value, seconds = unname(proc.time()[["elapsed"]] - started))
}

.nm_profile_repeated <- function(fun, repetitions) {
  fun()
  started <- proc.time()[["elapsed"]]
  for (iteration in seq_len(repetitions)) fun()
  elapsed <- unname(proc.time()[["elapsed"]] - started)
  list(
    repetitions = repetitions,
    total_seconds = elapsed,
    microseconds_per_call = elapsed * 1e6 / repetitions
  )
}

#' Profile persistent prediction-tape recording and evaluation
#'
#' This diagnostic separates engine compilation, CppAD recording, repeated
#' value evaluation, and repeated Jacobian evaluation. CppAD's operation
#' sequence and work-vector allocations are reported as a resident-memory
#' proxy; they are not a process-level peak-memory measurement.
#'
#' @param model A validated [nm_model()].
#' @param data A NONMEM-style event dataset.
#' @param repetitions Positive number of repeated tape evaluations.
#' @param jacobian Also profile exact Jacobian evaluation.
#' @return A one-row data frame with timing, tape-size, solver, and derivative
#'   strategy telemetry.
#' @export
nm_tape_profile <- function(model, data, repetitions = 20L, jacobian = TRUE) {
  if (!inherits(model, "nm_model")) .nm_stop("`model` must be an nm_model.")
  repetitions <- as.integer(repetitions)
  if (length(repetitions) != 1L || is.na(repetitions) || repetitions < 1L) {
    .nm_stop("`repetitions` must be a positive integer.")
  }
  data <- if (inherits(data, "nm_dataset")) data else nm_dataset(data)
  compiled <- .nm_profile_elapsed(nm_compile(model))
  engine <- compiled$value
  recorded <- .nm_profile_elapsed(engine$prediction_tape(data))
  tape <- recorded$value
  value_timing <- .nm_profile_repeated(function() {
    .liberation_prediction_tape_eval(tape$pointer, tape$point, FALSE)
  }, repetitions)
  jacobian_timing <- if (isTRUE(jacobian)) {
    .nm_profile_repeated(function() {
      .liberation_prediction_tape_eval(tape$pointer, tape$point, TRUE)
    }, repetitions)
  } else {
    list(microseconds_per_call = NA_real_)
  }
  info <- .liberation_prediction_tape_info(tape$pointer)
  data.frame(
    advan = model$ADVAN,
    solver = model$SOLVER,
    states = model$n_state,
    rows = nrow(data),
    subjects = length(unique(data$ID)),
    rtol = as.numeric(model$ODE_CONTROL$rtol %||% NA_real_),
    atol = as.numeric(model$ODE_CONTROL$atol %||% NA_real_),
    compile_seconds = compiled$seconds,
    record_seconds = recorded$seconds,
    value_microseconds = value_timing$microseconds_per_call,
    jacobian_microseconds = jacobian_timing$microseconds_per_call,
    operations = info$operations,
    operator_arguments = info$operator_arguments,
    variables = info$variables,
    parameters = info$parameters,
    dynamic_parameters = info$dynamic_parameters,
    operation_sequence_bytes = info$operation_sequence_bytes,
    random_access_bytes = info$random_access_bytes,
    taylor_bytes_proxy = info$taylor_bytes_proxy,
    resident_bytes_proxy = info$resident_bytes_proxy,
    propagation_kernel = as.character(info$propagation_kernel),
    derivative_strategy = as.character(info$derivative_strategy),
    jacobian_nonzeros = info$jacobian_nonzeros,
    stringsAsFactors = FALSE
  )
}

.nm_checkpoint_assessment <- function(probe) {
  cases <- unname(probe[c("advan1", "matrix2")])
  do.call(rbind, lapply(cases, function(value) {
    operation_ratio <- value$checkpoint_operations / value$direct_operations
    runtime_ratio <- value$checkpoint_microseconds / value$direct_microseconds
    data.frame(
      kernel = value$name,
      operation_ratio = operation_ratio,
      runtime_ratio = runtime_ratio,
      exact = value$max_value_difference <= 1e-12 &&
        value$max_jacobian_difference <= 1e-10,
      nested_ad_safe = isTRUE(value$nested_ad_safe),
      production_candidate = isTRUE(value$nested_ad_safe) &&
        operation_ratio < 0.8 && runtime_ratio < 1.15,
      stringsAsFactors = FALSE
    )
  }))
}

#' Benchmark nonlinear ODE tape scalability
#'
#' Recompiles an ADVAN6 or ADVAN13 model across solver tolerances and profiles
#' persistent tape size and repeated derivative evaluation. A guarded CppAD
#' checkpoint probe is returned as decision support; this function does not
#' silently switch production propagation kernels.
#'
#' @param model An ADVAN6 or ADVAN13 [nm_model()].
#' @param data A NONMEM-style event dataset.
#' @param tolerances Positive relative tolerances.
#' @param repetitions Repeated evaluation count per tolerance.
#' @param jacobian Profile the exact prediction Jacobian.
#' @param checkpoint_repetitions Repeated inner-kernel count for the checkpoint
#'   tape-size probe.
#' @param checkpoint_evaluations Timed evaluations in the checkpoint probe.
#' @return An `nm_ode_tape_benchmark` list containing `profiles`,
#'   `checkpoint`, and `assessment`.
#' @export
nm_ode_tape_benchmark <- function(
    model, data, tolerances = c(1e-5, 1e-7, 1e-9),
    repetitions = 10L, jacobian = TRUE,
    checkpoint_repetitions = 64L, checkpoint_evaluations = 1000L) {
  if (!inherits(model, "nm_model") || !isTRUE(model$USE_ODE)) {
    .nm_stop("`model` must be a validated ADVAN6 or ADVAN13 ODE model.")
  }
  tolerances <- as.numeric(tolerances)
  if (!length(tolerances) || any(!is.finite(tolerances)) ||
      any(tolerances <= 0)) {
    .nm_stop("`tolerances` must contain positive finite values.")
  }
  checkpoint_repetitions <- as.integer(checkpoint_repetitions)
  checkpoint_evaluations <- as.integer(checkpoint_evaluations)
  if (length(checkpoint_repetitions) != 1L ||
      is.na(checkpoint_repetitions) || checkpoint_repetitions < 1L ||
      length(checkpoint_evaluations) != 1L ||
      is.na(checkpoint_evaluations) || checkpoint_evaluations < 1L) {
    .nm_stop(
      "`checkpoint_repetitions` and `checkpoint_evaluations` must be ",
      "positive integers."
    )
  }
  profiles <- lapply(tolerances, function(tolerance) {
    candidate <- model
    candidate$ODE_CONTROL$rtol <- tolerance
    candidate$ODE_CONTROL$atol <- min(tolerance * 1e-2, tolerance)
    nm_tape_profile(
      candidate, data, repetitions = repetitions, jacobian = jacobian
    )
  })
  checkpoint <- LibeRtAD::ad_checkpoint_probe(
    repetitions = checkpoint_repetitions,
    evaluations = checkpoint_evaluations
  )
  structure(
    list(
      profiles = do.call(rbind, profiles),
      checkpoint = checkpoint,
      assessment = .nm_checkpoint_assessment(checkpoint),
      memory_note = paste(
        "`resident_bytes_proxy` covers CppAD operation sequences and current",
        "work vectors; it is not process peak memory."
      )
    ),
    class = "nm_ode_tape_benchmark"
  )
}

#' @export
print.nm_ode_tape_benchmark <- function(x, ...) {
  cat("LibeRation ODE tape scalability benchmark\n")
  print(x$profiles, row.names = FALSE)
  cat("\nCheckpoint decision support\n")
  print(x$assessment, row.names = FALSE)
  invisible(x)
}
