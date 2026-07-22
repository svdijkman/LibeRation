.nm_validation_values <- function(value, predictions = NULL) {
  if (inherits(value, "nm_fit")) {
    if (is.null(predictions)) {
      predictions <- tryCatch(predict(value), error = function(error) NULL)
    }
    return(list(
      objective = value$objective, theta = value$theta, omega = value$omega,
      sigma = value$sigma, eta = value$eta,
      predictions = if (is.data.frame(predictions) && "IPRED" %in% names(predictions))
        predictions$IPRED else predictions
    ))
  }
  if (is.list(value)) return(value)
  .nm_stop("Validation input must be an nm_fit or named list.")
}

.nm_validation_compare <- function(reference, candidate, tolerance, component) {
  if (is.null(reference) && is.null(candidate)) {
    return(data.frame(component = component, available = FALSE, passed = NA,
                      n = 0L, max_absolute = NA_real_, max_relative = NA_real_,
                      normalized_error = NA_real_, stringsAsFactors = FALSE))
  }
  if (is.null(reference) || is.null(candidate)) {
    return(data.frame(component = component, available = FALSE, passed = FALSE,
                      n = 0L, max_absolute = Inf, max_relative = Inf,
                      normalized_error = Inf, stringsAsFactors = FALSE))
  }
  reference <- as.numeric(reference)
  candidate <- as.numeric(candidate)
  if (length(reference) != length(candidate) || !length(reference)) {
    return(data.frame(component = component, available = TRUE, passed = FALSE,
                      n = min(length(reference), length(candidate)), max_absolute = Inf,
                      max_relative = Inf, normalized_error = Inf, stringsAsFactors = FALSE))
  }
  valid <- is.finite(reference) & is.finite(candidate)
  if (!all(valid) || !any(valid)) {
    return(data.frame(component = component, available = TRUE, passed = FALSE,
                      n = length(reference), max_absolute = Inf, max_relative = Inf,
                      normalized_error = Inf, stringsAsFactors = FALSE))
  }
  difference <- abs(candidate - reference)
  absolute <- as.numeric(tolerance$absolute %||% 0)
  relative <- as.numeric(tolerance$relative %||% 0)
  floor <- as.numeric(tolerance$floor %||% max(absolute, .Machine$double.eps))
  permitted <- absolute + relative * pmax(abs(reference), floor)
  normalized <- difference / pmax(permitted, .Machine$double.eps)
  data.frame(
    component = component, available = TRUE, passed = all(normalized <= 1),
    n = length(reference), max_absolute = max(difference),
    max_relative = max(difference / pmax(abs(reference), floor)),
    normalized_error = max(normalized), stringsAsFactors = FALSE
  )
}

#' Require numerical agreement before publishing a benchmark
#'
#' @param reference,candidate Reference-engine and LibeRation fit summaries or
#'   `nm_fit` objects.
#' @param reference_predictions,candidate_predictions Optional prediction
#'   vectors, avoiding a fresh `predict()` call.
#' @param tolerances Named per-component absolute/relative tolerances.
#' @param required Components that must be present and pass.
#' @return An `nm_validation_gate` report.
#' @export
nm_validation_gate <- function(
    reference, candidate, reference_predictions = NULL, candidate_predictions = NULL,
    tolerances = list(
      objective = list(absolute = 1e-3, relative = 1e-4),
      theta = list(absolute = 1e-5, relative = 0.02),
      omega = list(absolute = 1e-5, relative = 0.10),
      sigma = list(absolute = 1e-5, relative = 0.10),
      eta = list(absolute = 1e-4, relative = 0.10),
      predictions = list(absolute = 1e-6, relative = 0.02)
    ),
    required = c("objective", "theta", "omega", "sigma", "eta", "predictions")) {
  reference <- .nm_validation_values(reference, reference_predictions)
  candidate <- .nm_validation_values(candidate, candidate_predictions)
  components <- unique(c(names(tolerances), required))
  rows <- lapply(components, function(component) .nm_validation_compare(
    reference[[component]], candidate[[component]],
    tolerances[[component]] %||% list(absolute = 0, relative = 0), component
  ))
  table <- do.call(rbind, rows)
  required_rows <- match(required, table$component)
  passed <- !anyNA(required_rows) && all(table$available[required_rows]) &&
    all(table$passed[required_rows])
  structure(list(
    passed = passed, comparisons = table, required = required,
    tolerances = tolerances, checked = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")
  ), class = "nm_validation_gate")
}

#' Record reproducible benchmark provenance
#'
#' @param label Benchmark label.
#' @param repetitions Number of measured repetitions after warm-up.
#' @param warmup Number of discarded warm-up runs.
#' @param extra Additional named metadata.
#' @export
nm_benchmark_provenance <- function(label, repetitions, warmup = 1L, extra = list()) {
  packages <- c("LibeRtAD", "LibeRation", "LibeRties", "LibeRality")
  versions <- stats::setNames(vapply(packages, function(package) {
    if (requireNamespace(package, quietly = TRUE)) as.character(utils::packageVersion(package)) else NA_character_
  }, character(1)), packages)
  sha <- Sys.getenv("GITHUB_SHA", unset = "")
  if (!nzchar(sha)) sha <- NA_character_
  utils::modifyList(list(
    schema = "liber.benchmark/1", label = as.character(label),
    repetitions = as.integer(repetitions), warmup = as.integer(warmup),
    recorded = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC"),
    git_commit = sha, packages = versions, engine = LibeRtAD::ad_engine_info(),
    runtime = list(R = R.version.string, platform = R.version$platform,
                   machine = Sys.info()[["machine"]] %||% "unknown")
  ), extra)
}

#' Combine timing evidence with a correctness gate
#'
#' @param validation Result from [nm_validation_gate()].
#' @param timings Named timing summary.
#' @param provenance Result from [nm_benchmark_provenance()].
#' @return A benchmark record whose `publishable` flag is true only after a
#'   passing numerical gate and at least three measured repetitions.
#' @export
nm_benchmark_gate <- function(validation, timings, provenance) {
  if (!inherits(validation, "nm_validation_gate")) .nm_stop("`validation` must be an nm_validation_gate.")
  if (!is.list(timings) || is.null(names(timings))) .nm_stop("`timings` must be a named list.")
  if (!is.list(provenance) || !identical(provenance$schema %||% "", "liber.benchmark/1")) {
    .nm_stop("`provenance` must be created by nm_benchmark_provenance().")
  }
  structure(list(
    schema = "liber.benchmark-result/1", version = 1L,
    publishable = isTRUE(validation$passed) && provenance$repetitions >= 3L,
    validation = validation, timings = timings, provenance = provenance
  ), class = "nm_benchmark_result")
}

#' @export
print.nm_validation_gate <- function(x, ...) {
  cat("LibeR numerical validation gate:", if (x$passed) "PASS" else "FAIL", "\n")
  print(x$comparisons, row.names = FALSE)
  invisible(x)
}
