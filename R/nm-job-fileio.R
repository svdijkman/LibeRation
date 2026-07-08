#' Collision-tolerant file I/O helpers for concurrent job workers/readers.
#'
#' On Windows, a process reading a file (readLines/readRDS) while another
#' process writes it (append/saveRDS) can hit a sharing violation surfaced as
#' "cannot open the connection" / "Permission denied". These helpers retry
#' briefly so a transient collision never crashes a worker or a status poll.
#'
#' @keywords internal
.nm_io_retry_max <- function() {
  as.integer(getOption("LibeRation.io_retry_max", 100L))
}

#' @keywords internal
.nm_io_retry_sleep <- function() {
  as.numeric(getOption("LibeRation.io_retry_sleep", 0.015))
}

#' Append text to a file, retrying on transient open failures. Never throws.
#' @keywords internal
.nm_file_append <- function(path, text) {
  if (is.null(path) || !nzchar(path) || is.null(text) || length(text) == 0L) {
    return(invisible(FALSE))
  }
  retries <- .nm_io_retry_max()
  sleep <- .nm_io_retry_sleep()
  for (i in seq_len(retries)) {
    con <- suppressWarnings(tryCatch(file(path, open = "a"), error = function(e) NULL))
    if (!is.null(con)) {
      ok <- tryCatch({
        writeLines(text, con, useBytes = TRUE)
        flush(con)
        TRUE
      }, error = function(e) FALSE)
      try(close(con), silent = TRUE)
      if (isTRUE(ok)) {
        return(invisible(TRUE))
      }
    }
    Sys.sleep(sleep)
  }
  invisible(FALSE)
}

#' Overwrite a file with text lines, retrying on transient failures. Never throws.
#' @keywords internal
.nm_file_write_lines <- function(path, text) {
  if (is.null(path) || !nzchar(path)) {
    return(invisible(FALSE))
  }
  retries <- .nm_io_retry_max()
  sleep <- .nm_io_retry_sleep()
  for (i in seq_len(retries)) {
    con <- suppressWarnings(tryCatch(file(path, open = "w"), error = function(e) NULL))
    if (!is.null(con)) {
      ok <- tryCatch({
        writeLines(as.character(text), con, useBytes = TRUE)
        flush(con)
        TRUE
      }, error = function(e) FALSE)
      try(close(con), silent = TRUE)
      if (isTRUE(ok)) {
        return(invisible(TRUE))
      }
    }
    Sys.sleep(sleep)
  }
  invisible(FALSE)
}

#' Read all lines from a file, retrying on transient failures.
#' Returns character(0) if the file is missing or unreadable after retries.
#' @keywords internal
.nm_file_read_lines <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(character(0))
  }
  retries <- .nm_io_retry_max()
  sleep <- .nm_io_retry_sleep()
  sentinel <- structure(list(), class = ".nm_io_fail")
  for (i in seq_len(retries)) {
    res <- suppressWarnings(tryCatch(
      readLines(path, warn = FALSE),
      error = function(e) sentinel
    ))
    if (!inherits(res, ".nm_io_fail")) {
      return(res)
    }
    Sys.sleep(sleep)
  }
  character(0)
}

#' Read an RDS file, retrying on transient failures. Returns `default` on failure.
#' @keywords internal
.nm_read_rds_safe <- function(path, default = NULL) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(default)
  }
  retries <- .nm_io_retry_max()
  sleep <- .nm_io_retry_sleep()
  sentinel <- structure(list(), class = ".nm_io_fail")
  for (i in seq_len(retries)) {
    res <- suppressWarnings(tryCatch(readRDS(path), error = function(e) sentinel))
    if (!inherits(res, ".nm_io_fail")) {
      return(res)
    }
    Sys.sleep(sleep)
  }
  default
}

#' Write an RDS file, retrying on transient open failures. Never throws.
#'
#' A concurrent reader hitting the file mid-write gets a sharing violation
#' (failed open) or a truncated-stream error, both of which surface as an
#' error the reader retries away; readers never observe a partial object
#' silently. So a retrying direct write plus retrying reads is sufficient and
#' avoids Windows `file.rename` (which cannot overwrite an existing target).
#' @keywords internal
.nm_save_rds_safe <- function(obj, path) {
  if (is.null(path) || !nzchar(path)) {
    return(invisible(FALSE))
  }
  retries <- .nm_io_retry_max()
  sleep <- .nm_io_retry_sleep()
  for (i in seq_len(retries)) {
    ok <- suppressWarnings(tryCatch({
      saveRDS(obj, path)
      TRUE
    }, error = function(e) FALSE))
    if (isTRUE(ok)) {
      return(invisible(TRUE))
    }
    Sys.sleep(sleep)
  }
  invisible(FALSE)
}
