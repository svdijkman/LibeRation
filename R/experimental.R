#' Configure experimental engine behaviour
#'
#' Experimental models retain their approximation, reproducibility, and
#' provenance metadata in the serialized model and every queued job.
#'
#' @param enabled Explicitly acknowledge use of an experimental engine.
#' @param label Optional user-facing purpose or validation label.
#' @param strict Fail rather than fall back to a less exact algorithm.
#' @return Serializable experimental-engine metadata.
#' @export
nm_experimental_config <- function(enabled = TRUE, label = NULL, strict = TRUE) {
  if (length(enabled) != 1L || is.na(enabled) ||
      length(strict) != 1L || is.na(strict)) {
    .nm_stop("`enabled` and `strict` must be TRUE or FALSE.")
  }
  if (!is.null(label)) {
    label <- trimws(as.character(label))
    if (length(label) != 1L || is.na(label) || !nzchar(label)) {
      .nm_stop("`label` must be NULL or one non-empty string.")
    }
  }
  structure(
    list(version = 1L, enabled = isTRUE(enabled), label = label,
         strict = isTRUE(strict), features = character()),
    class = "nm_experimental_config"
  )
}

.nm_experimental_config <- function(config, features = character()) {
  if (is.null(config)) {
    if (length(features)) {
      .nm_stop(
        "Experimental model features require explicit acknowledgement with ",
        "`EXPERIMENTAL = nm_experimental_config(enabled = TRUE)`."
      )
    }
    config <- nm_experimental_config(enabled = FALSE)
  }
  if (!inherits(config, "nm_experimental_config")) {
    if (!is.list(config)) .nm_stop("EXPERIMENTAL must be created by `nm_experimental_config()`.")
    config$version <- NULL
    config <- do.call(nm_experimental_config, config)
  }
  if (length(features) && !isTRUE(config$enabled)) {
    .nm_stop(
      "Experimental model features require `EXPERIMENTAL = ",
      "nm_experimental_config(enabled = TRUE)`."
    )
  }
  config$features <- sort(unique(as.character(features)))
  config
}

#' @export
print.nm_experimental_config <- function(x, ...) {
  cat("LibeRation experimental engine\n")
  cat("  enabled:", x$enabled, " strict:", x$strict, "\n")
  if (length(x$features)) cat("  features:", paste(x$features, collapse = ", "), "\n")
  invisible(x)
}
