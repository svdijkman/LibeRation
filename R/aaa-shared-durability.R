# Generated from tools/shared/liber-durability.R.
# Run `Rscript tools/sync-shared-runtime.R`; do not edit this copy directly.

.liber_shared_fail <- function(error, message) {
  if (is.function(error)) return(error(message))
  stop(message, call. = FALSE)
}

.liber_shared_acl_warning <- local({
  warned <- FALSE
  function(message) {
    if (!warned) {
      warning(message, call. = FALSE)
      warned <<- TRUE
    }
    invisible(FALSE)
  }
})

.liber_shared_windows_acl <- function(path, error = NULL) {
  icacls <- unname(Sys.which("icacls.exe"))
  whoami <- unname(Sys.which("whoami.exe"))
  failure <- function(detail) {
    message <- paste0(
      "Unable to apply the owner/SYSTEM/Administrators Windows ACL to ",
      path, ": ", detail,
      ". Encrypted LibeR workspaces remain encrypted, but filesystem ACL ",
      "defence in depth is unavailable."
    )
    if (isTRUE(getOption("LibeR.strict_windows_acl", FALSE))) {
      .liber_shared_fail(error, message)
    }
    .liber_shared_acl_warning(message)
  }
  if (!nzchar(icacls) || !nzchar(whoami)) {
    return(failure("icacls.exe or whoami.exe was not found"))
  }
  identity <- tryCatch(
    suppressWarnings(system2(
      whoami, c("/user", "/fo", "csv", "/nh"),
      stdout = TRUE, stderr = TRUE
    )), error = identity
  )
  if (inherits(identity, "error") || !length(identity)) {
    return(failure("the current Windows security identifier could not be read"))
  }
  row <- tryCatch(
    utils::read.csv(
      text = paste(identity, collapse = "\n"), header = FALSE,
      stringsAsFactors = FALSE
    ), error = identity
  )
  sid <- if (!inherits(row, "error") && ncol(row) >= 2L) {
    trimws(as.character(row[[2L]][[1L]]))
  } else ""
  if (!grepl("^S-[0-9-]+$", sid)) {
    return(failure("the current Windows security identifier was invalid"))
  }
  output <- tryCatch(
    suppressWarnings(system2(
      icacls,
      c(
        shQuote(normalizePath(path, winslash = "\\", mustWork = TRUE)),
        "/inheritance:r", "/grant:r",
        shQuote(paste0("*", sid, ":(F)")),
        shQuote("*S-1-5-18:(F)"),
        shQuote("*S-1-5-32-544:(F)"), "/q"
      ),
      stdout = TRUE, stderr = TRUE
    )), error = identity
  )
  status <- if (inherits(output, "error")) 1L else attr(output, "status")
  if (is.null(status)) status <- 0L
  if (!identical(as.integer(status), 0L)) {
    detail <- if (inherits(output, "error")) conditionMessage(output) else
      paste(output, collapse = " ")
    return(failure(detail))
  }
  invisible(TRUE)
}

.liber_shared_publish_file <- function(
    temporary, path, mode = "0600", error = NULL) {
  backup <- paste0(path, ".previous")
  had_previous <- file.exists(path)
  if (file.exists(backup)) unlink(backup, force = TRUE)
  if (had_previous && !file.rename(path, backup)) {
    .liber_shared_fail(error, paste0("Unable to rotate file: ", path))
  }
  if (!file.rename(temporary, path)) {
    if (had_previous && file.exists(backup)) file.rename(backup, path)
    .liber_shared_fail(error, paste0("Unable to publish file: ", path))
  }
  if (file.exists(backup)) unlink(backup, force = TRUE)
  if (!is.null(mode)) {
    if (.Platform$OS.type == "windows") {
      .liber_shared_windows_acl(path, error = error)
    } else {
      try(Sys.chmod(path, mode = mode, use_umask = FALSE), silent = TRUE)
    }
  }
  invisible(path)
}

.liber_shared_atomic_publish <- function(
    path, writer, prefix = "write-", fileext = ".tmp", mode = "0600",
    error = NULL) {
  directory <- dirname(path)
  if (!dir.exists(directory) &&
      !dir.create(directory, recursive = TRUE, showWarnings = FALSE)) {
    .liber_shared_fail(error, paste0("Unable to create directory: ", directory))
  }
  temporary <- tempfile(prefix, tmpdir = directory, fileext = fileext)
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  writer(temporary)
  .liber_shared_publish_file(
    temporary, path, mode = mode, error = error
  )
}

.liber_shared_durable_read <- function(
    path, reader, attempts = 1L, delay = 0.01, missing = NULL,
    warn_recovery = TRUE, recovery_label = "durable", error = NULL) {
  attempts <- max(1L, as.integer(attempts))
  last <- NULL
  for (attempt in seq_len(attempts)) {
    candidates <- unique(c(path, paste0(path, ".previous")))
    candidates <- candidates[file.exists(candidates)]
    for (candidate in candidates) {
      outcome <- tryCatch(
        list(ok = TRUE, value = reader(candidate)),
        error = function(condition) {
          last <<- condition
          list(ok = FALSE, value = NULL)
        }
      )
      if (isTRUE(outcome$ok)) {
        if (isTRUE(warn_recovery) && !identical(candidate, path)) {
          warning(
            "Recovered interrupted ", recovery_label, " write from ",
            basename(candidate), ".",
            call. = FALSE
          )
        }
        return(outcome$value)
      }
    }
    if (attempt < attempts) Sys.sleep(max(0, as.numeric(delay)))
  }
  if (!is.null(missing) && is.null(last)) return(missing)
  detail <- if (inherits(last, "condition")) {
    conditionMessage(last)
  } else {
    "file is absent"
  }
  .liber_shared_fail(error, paste0("Unable to read ", path, ": ", detail))
}

.liber_shared_lock_acquire <- function(
    path, timeout = 5, stale_after = Inf, poll = 0.01, owner = NULL,
    error = NULL) {
  directory <- dirname(path)
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  }
  started <- proc.time()[["elapsed"]]
  repeat {
    if (dir.create(path, showWarnings = FALSE)) break
    info <- suppressWarnings(file.info(path))
    age <- if (nrow(info) && !is.na(info$mtime)) {
      as.numeric(difftime(Sys.time(), info$mtime, units = "secs"))
    } else {
      0
    }
    if (is.finite(age) && is.finite(stale_after) && age > stale_after) {
      unlink(path, recursive = TRUE, force = TRUE)
      next
    }
    if (proc.time()[["elapsed"]] - started >= as.numeric(timeout)) {
      .liber_shared_fail(error, paste0("Timed out acquiring lock: ", basename(path)))
    }
    Sys.sleep(max(0.001, as.numeric(poll)))
  }
  if (!is.null(owner)) {
    try(writeLines(enc2utf8(as.character(owner)), file.path(path, "owner"),
                   useBytes = TRUE), silent = TRUE)
  }
  path
}

.liber_shared_lock_release <- function(path) {
  if (!is.null(path) && dir.exists(path)) {
    unlink(path, recursive = TRUE, force = TRUE)
  }
  invisible(TRUE)
}

.liber_shared_with_lock <- function(
    path, operation, timeout = 5, stale_after = Inf, poll = 0.01,
    owner = NULL, error = NULL) {
  .liber_shared_lock_acquire(
    path, timeout = timeout, stale_after = stale_after, poll = poll,
    owner = owner, error = error
  )
  on.exit(.liber_shared_lock_release(path), add = TRUE)
  operation()
}

.liber_shared_component <- function(value, what = "path component",
                                    max_length = 128L, error = NULL) {
  value <- as.character(value)
  valid <- length(value) == 1L && !is.na(value) && nzchar(value) &&
    nchar(value, type = "bytes") <= as.integer(max_length) &&
    grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", value) &&
    !value %in% c(".", "..")
  if (!valid) .liber_shared_fail(error, paste0("Invalid ", what, "."))
  value
}
