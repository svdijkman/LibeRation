# Generated from tools/shared/liber-async.R.
# Run `Rscript tools/sync-shared-runtime.R`; do not edit this copy directly.

.liber_shared_task_registry <- function(session = NULL) {
  registry <- new.env(parent = emptyenv())
  registry$jobs <- list()
  registry$completed <- list()
  registry$session <- session
  class(registry) <- "liber_task_registry"
  if (!is.null(session) && is.function(session$onSessionEnded)) {
    session$onSessionEnded(function() {
      .liber_shared_task_cancel_all(registry)
    })
  }
  registry
}

.liber_shared_task_require_registry <- function(registry) {
  if (!inherits(registry, "liber_task_registry")) {
    stop("Invalid LibeR background-task registry.", call. = FALSE)
  }
  registry
}

.liber_shared_task_id <- function(prefix = "task") {
  paste0(
    gsub("[^A-Za-z0-9_-]", "-", as.character(prefix)[[1L]]), "-",
    format(Sys.time(), "%Y%m%d%H%M%OS6"), "-",
    sprintf("%08x", sample.int(.Machine$integer.max, 1L))
  )
}

.liber_shared_task_start <- function(
    registry, package, fun, args = list(), label = fun, metadata = list(),
    replace = FALSE) {
  registry <- .liber_shared_task_require_registry(registry)
  if (!requireNamespace("callr", quietly = TRUE)) {
    stop(
      "Non-blocking GUI tasks require the `callr` package.",
      call. = FALSE
    )
  }
  active <- names(registry$jobs)
  if (length(active)) {
    if (!isTRUE(replace)) {
      stop("Another background GUI task is already running.", call. = FALSE)
    }
    .liber_shared_task_cancel_all(registry)
  }
  package <- as.character(package)[[1L]]
  fun <- as.character(fun)[[1L]]
  if (!nzchar(package) || !nzchar(fun) || !is.list(args)) {
    stop("Invalid background GUI task specification.", call. = FALSE)
  }
  namespace <- asNamespace(package)
  specification <- getNamespaceInfo(namespace, "spec")
  expected_version <- as.character(specification[["version"]])
  if (!length(expected_version) || is.na(expected_version)) {
    expected_version <- ""
  }
  package_path <- getNamespaceInfo(namespace, "path")
  source_files <- if (
    length(package_path) == 1L && dir.exists(file.path(package_path, "R"))
  ) {
    list.files(file.path(package_path, "R"), pattern = "\\.[Rr]$")
  } else character()
  source_path <- if (length(source_files)) {
    normalizePath(package_path, winslash = "/", mustWork = TRUE)
  } else ""
  id <- .liber_shared_task_id(fun)
  process <- callr::r_bg(
    function(package, fun, args, expected_version, source_path) {
      if (nzchar(source_path)) {
        if (!requireNamespace("pkgload", quietly = TRUE)) {
          stop(
            "A source-loaded LibeR GUI requires `pkgload` in background ",
            "workers. Install the package or install `pkgload`.",
            call. = FALSE
          )
        }
        pkgload::load_all(
          source_path, quiet = TRUE, export_all = FALSE, helpers = FALSE
        )
      }
      namespace <- asNamespace(package)
      actual <- as.character(
        getNamespaceInfo(namespace, "spec")[["version"]]
      )
      if (nzchar(expected_version) && !identical(actual, expected_version)) {
        stop(
          "The active ", package, " GUI is version ", expected_version,
          " but its background worker loaded version ", actual,
          ". Restart R after reinstalling the package.",
          call. = FALSE
        )
      }
      target <- get(fun, envir = namespace, inherits = FALSE)
      do.call(target, args)
    },
    args = list(
      package = package, fun = fun, args = args,
      expected_version = expected_version, source_path = source_path
    ),
    libpath = .libPaths(), stdout = "|", stderr = "|",
    supervise = TRUE
  )
  registry$jobs[[id]] <- list(
    id = id, package = package, fun = fun, label = as.character(label)[[1L]],
    metadata = metadata, process = process, stdout = character(),
    stderr = character(), started_at = Sys.time(), status = "running"
  )
  id
}

.liber_shared_task_active <- function(registry) {
  registry <- .liber_shared_task_require_registry(registry)
  length(registry$jobs) > 0L
}

.liber_shared_task_snapshot <- function(registry) {
  registry <- .liber_shared_task_require_registry(registry)
  if (!length(registry$jobs)) {
    return(list(running = FALSE, id = "", label = "", cancellable = FALSE))
  }
  job <- registry$jobs[[1L]]
  list(
    running = TRUE, id = job$id, label = job$label,
    cancellable = TRUE, started_at = format(
      job$started_at, "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC"
    )
  )
}

.liber_shared_task_poll <- function(registry) {
  registry <- .liber_shared_task_require_registry(registry)
  if (!length(registry$jobs)) return(invisible(FALSE))
  for (id in names(registry$jobs)) {
    job <- registry$jobs[[id]]
    process <- job$process
    output <- tryCatch(process$read_output_lines(), error = function(error) character())
    errors <- tryCatch(process$read_error_lines(), error = function(error) character())
    if (length(output)) job$stdout <- c(job$stdout, output)
    if (length(errors)) job$stderr <- c(job$stderr, errors)
    if (isTRUE(process$is_alive())) {
      registry$jobs[[id]] <- job
      next
    }
    output <- tryCatch(process$read_all_output_lines(), error = function(error) character())
    errors <- tryCatch(process$read_all_error_lines(), error = function(error) character())
    if (length(output)) job$stdout <- c(job$stdout, output)
    if (length(errors)) job$stderr <- c(job$stderr, errors)
    value <- tryCatch(process$get_result(), error = identity)
    job$finished_at <- Sys.time()
    job$process <- NULL
    if (inherits(value, "error")) {
      job$status <- "failed"
      job$error <- conditionMessage(value)
      job$result <- NULL
    } else {
      job$status <- "completed"
      job$error <- NULL
      job$result <- value
    }
    registry$completed[[id]] <- job
    registry$jobs[[id]] <- NULL
  }
  invisible(TRUE)
}

.liber_shared_task_take_completed <- function(registry) {
  registry <- .liber_shared_task_require_registry(registry)
  completed <- registry$completed
  registry$completed <- list()
  unname(completed)
}

.liber_shared_task_cancel <- function(registry, id = NULL) {
  registry <- .liber_shared_task_require_registry(registry)
  ids <- names(registry$jobs)
  if (!is.null(id)) ids <- intersect(ids, as.character(id))
  for (current in ids) {
    job <- registry$jobs[[current]]
    process <- job$process
    if (!is.null(process) && isTRUE(process$is_alive())) {
      try(process$kill(), silent = TRUE)
    }
    job$process <- NULL
    job$status <- "cancelled"
    job$error <- NULL
    job$result <- NULL
    job$finished_at <- Sys.time()
    registry$completed[[current]] <- job
    registry$jobs[[current]] <- NULL
  }
  invisible(length(ids) > 0L)
}

.liber_shared_task_cancel_all <- function(registry) {
  .liber_shared_task_cancel(registry)
}

.liber_shared_task_notify <- function(session, input_id, registry) {
  if (is.null(session) || !is.function(session$sendCustomMessage)) {
    return(invisible(NULL))
  }
  state <- .liber_shared_task_snapshot(registry)
  state$inputId <- as.character(input_id)[[1L]]
  session$sendCustomMessage("liber-task-state", state)
  invisible(state)
}
