.liber_compatibility_fallback <- function() {
  list(
    schema = "liber.ecosystem/1", release = "0.9.0-research-beta.3",
    packages = list(
      LibeRtAD = list(version = "0.7.9", required = TRUE),
      LibeRation = list(version = "0.9.2", required = TRUE),
      LibeRary = list(version = "0.7.5", required = FALSE),
      LibeRator = list(version = "0.2.6", required = FALSE),
      LibeRality = list(version = "0.2.4", required = FALSE),
      LibeRties = list(version = "0.7.4", required = FALSE)
    ),
    contracts = list(model = 2L, job = 2L, result = 2L,
                     liberation_workspace = 2L, liberties_queue = 2L)
  )
}

.liber_compatibility <- function() {
  path <- system.file("ecosystem", "compatibility.json", package = "LibeRation")
  if (!nzchar(path) || !file.exists(path)) return(.liber_compatibility_fallback())
  tryCatch(jsonlite::fromJSON(path, simplifyVector = FALSE),
           error = function(error) .liber_compatibility_fallback())
}

.liber_version_status <- function(package, specification) {
  installed <- requireNamespace(package, quietly = TRUE)
  version <- if (installed) as.character(utils::packageVersion(package)) else NA_character_
  expected <- as.character(specification$version)
  compatible <- installed && utils::compareVersion(version, expected) >= 0L
  data.frame(
    package = package, installed = installed, version = version,
    expected = expected, required = isTRUE(specification$required),
    compatible = compatible, stringsAsFactors = FALSE
  )
}

#' Diagnose a LibeR ecosystem installation
#'
#' Reports package compatibility, compiled dependency provenance, wire-contract
#' versions, workspace health, and common deployment mistakes. It does not
#' modify the installation or workspace.
#'
#' @param workspace Optional [nm_workspace()] or workspace directory to inspect.
#' @param strict Throw an error if a required component is missing or
#'   incompatible.
#' @param verbose Print the report.
#' @return A `liber_diagnostics` list.
#' @export
liber_doctor <- function(workspace = NULL, strict = FALSE, verbose = interactive()) {
  compatibility <- .liber_compatibility()
  packages <- do.call(rbind, Map(
    .liber_version_status, names(compatibility$packages), compatibility$packages
  ))
  rownames(packages) <- NULL
  engine <- tryCatch(LibeRtAD::ad_engine_info(), error = function(error) {
    list(available = FALSE, error = conditionMessage(error))
  })
  workspace_status <- list(checked = FALSE)
  if (!is.null(workspace)) {
    workspace_status <- tryCatch({
      value <- if (inherits(workspace, "nm_workspace")) workspace else nm_workspace(workspace)
      probe <- tempfile("doctor-", tmpdir = value$path)
      writable <- tryCatch({
        writeLines("ok", probe, useBytes = TRUE); unlink(probe); TRUE
      }, error = function(error) FALSE)
      list(checked = TRUE, path = value$path, writable = writable,
           version = value$version %||% NA_integer_)
    }, error = function(error) list(checked = TRUE, healthy = FALSE,
                                    error = conditionMessage(error)))
  }
  queue <- if (requireNamespace("LibeRties", quietly = TRUE)) {
    tryCatch(LibeRties::ls_queue_capabilities(), error = function(error) {
      list(error = conditionMessage(error))
    })
  } else list(available = FALSE)
  problems <- character()
  required_bad <- packages$required & !packages$compatible
  optional_bad <- packages$installed & !packages$compatible
  if (any(required_bad)) {
    problems <- c(problems, paste0("Required package mismatch: ",
                                  paste(packages$package[required_bad], collapse = ", "), "."))
  }
  if (any(optional_bad)) {
    problems <- c(problems, paste0("Installed optional package mismatch: ",
                                  paste(packages$package[optional_bad], collapse = ", "), "."))
  }
  if (isTRUE(workspace_status$checked) && isFALSE(workspace_status$writable %||% FALSE)) {
    problems <- c(problems, "Workspace is unavailable or not writable.")
  }
  result <- structure(list(
    healthy = !length(problems), release = compatibility$release,
    packages = packages, contracts = compatibility$contracts,
    engine = engine, workspace = workspace_status, queue = queue,
    runtime = list(R = R.version.string, platform = R.version$platform,
                   make = unname(Sys.which("make")), compiler = unname(Sys.which("g++"))),
    problems = problems
  ), class = "liber_diagnostics")
  if (isTRUE(verbose)) print(result)
  if (isTRUE(strict) && length(problems)) .nm_stop(paste(problems, collapse = " "))
  invisible(result)
}

#' @export
print.liber_diagnostics <- function(x, ...) {
  cat("LibeR ecosystem diagnostics - ", x$release, "\n", sep = "")
  cat("Status: ", if (isTRUE(x$healthy)) "healthy" else "attention required", "\n", sep = "")
  print(x$packages, row.names = FALSE)
  if (length(x$problems)) cat("\n", paste0("- ", x$problems, collapse = "\n"), "\n", sep = "")
  invisible(x)
}
