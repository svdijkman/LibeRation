#' @keywords internal
.nm_pkg_version_label <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    return("not installed")
  }
  as.character(utils::packageVersion(pkg))
}

#' Collect run environment metadata for workspace persistence
#'
#' @param job_id Optional local job id (including \code{remote_*} stubs).
#' @param job_root Job directory root.
#' @param started Run start time (\code{POSIXct} or character).
#' @param finished Run finish time (\code{POSIXct} or character).
#' @return List with package versions, queue info, and timing.
#' @keywords internal
.nm_run_info_collect <- function(job_id = NULL,
                                 job_root = nm_job_root(),
                                 started = NULL,
                                 finished = NULL) {
  vers <- list(
    LibeRation = .nm_pkg_version_label("LibeRation"),
    LibeRtAD = .nm_pkg_version_label("LibeRtAD")
  )
  queue <- list(type = "local", label = "Local job queue", job_id = job_id %||% "")
  if (!is.null(job_id) && nzchar(job_id)) {
    meta <- .nm_job_read_meta(job_id, job_root)
    if (!is.null(meta) && isTRUE(meta$remote)) {
      entry <- tryCatch(
        .nm_job_remote_entry(job_id, job_root),
        error = function(e) NULL
      )
      queue <- list(
        type = "remote",
        label = paste0(
          entry$name %||% entry$id %||% "Remote",
          if (nzchar(entry$base_url %||% "")) {
            paste0(" (", entry$base_url, ")")
          } else {
            ""
          }
        ),
        server_id = entry$id %||% "",
        server_name = entry$name %||% "",
        base_url = entry$base_url %||% "",
        remote_job_id = meta$remote_job_id %||% "",
        job_id = job_id
      )
      if (requireNamespace("LibeRties", quietly = TRUE)) {
        vers$LibeRties <- .nm_pkg_version_label("LibeRties")
      }
    }
  }
  dev_env <- getOption("LibeRation.job_dev_env", NULL)
  started_chr <- if (!is.null(started)) as.character(started) else ""
  finished_chr <- if (!is.null(finished)) as.character(finished) else ""
  duration_sec <- NA_real_
  if (!is.null(started) && !is.null(finished)) {
    duration_sec <- as.numeric(difftime(finished, started, units = "secs"))
    if (!is.finite(duration_sec) || duration_sec < 0) {
      duration_sec <- NA_real_
    }
  }
  list(
    versions = vers,
    queue = queue,
    dev_mode = identical(dev_env$mode, "dev"),
    started = started_chr,
    finished = finished_chr,
    duration_sec = duration_sec
  )
}

#' Format run duration for display
#'
#' @param seconds Duration in seconds.
#' @return Character label.
#' @keywords internal
.nm_run_info_format_duration <- function(seconds) {
  if (is.null(seconds) || length(seconds) != 1L || !is.finite(seconds)) {
    return("\u2014")
  }
  if (seconds < 60) {
    return(paste0(round(seconds, 1), " s"))
  }
  if (seconds < 3600) {
    mins <- floor(seconds / 60)
    secs <- round(seconds - mins * 60)
    return(paste0(mins, " min ", secs, " s"))
  }
  hrs <- floor(seconds / 3600)
  mins <- floor((seconds - hrs * 3600) / 60)
  paste0(hrs, " h ", mins, " min")
}

#' @keywords internal
.nm_run_info_as_rows <- function(run_info) {
  if (is.null(run_info) || length(run_info) == 0L) {
    return(data.frame(field = character(), value = character(), stringsAsFactors = FALSE))
  }
  vers <- run_info$versions %||% list()
  queue <- run_info$queue %||% list()
  rows <- list(
    c("Run started", run_info$started %||% ""),
    c("Run finished", run_info$finished %||% ""),
    c("Duration", .nm_run_info_format_duration(run_info$duration_sec)),
    c("Job queue", queue$label %||% if (identical(queue$type, "remote")) "Remote" else "Local"),
    c("LibeRation", vers$LibeRation %||% ""),
    c("LibeRtAD", vers$LibeRtAD %||% ""),
    c("LibeRties", vers$LibeRties %||% "")
  )
  if (isTRUE(run_info$dev_mode)) {
    rows <- c(rows, list(c("Package mode", "Development (load_all)")))
  }
  if (identical(queue$type, "remote") && nzchar(queue$remote_job_id %||% "")) {
    rows <- c(rows, list(c("Remote job id", queue$remote_job_id)))
  }
  if (nzchar(queue$job_id %||% "")) {
    rows <- c(rows, list(c("Local job id", queue$job_id)))
  }
  data.frame(
    field = vapply(rows, `[[`, character(1L), 1L),
    value = vapply(rows, function(x) {
      val <- x[[2L]]
      if (is.null(val) || (length(val) == 1L && !nzchar(as.character(val)))) {
        "\u2014"
      } else {
        as.character(val)
      }
    }, character(1L)),
    stringsAsFactors = FALSE
  )
}
