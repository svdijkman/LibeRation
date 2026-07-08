#' @keywords internal
.nm_job_progress_env <- function() {
  getOption("LibeRation.job_progress")
}

#' @keywords internal
.nm_job_progress_init <- function(job_path) {
  if (is.null(job_path) || !nzchar(job_path)) {
    return(invisible(NULL))
  }
  env <- new.env(parent = emptyenv())
  env$job_path <- job_path
  env$progress_path <- file.path(job_path, "worker.progress")
  env$log_path <- file.path(job_path, "worker.log")
  env$heartbeat_path <- file.path(job_path, "worker.heartbeat")
  env$grad_step <- 0L
  options(LibeRation.job_progress = env)
  invisible(env)
}

#' @keywords internal
.nm_job_progress_clear <- function() {
  options(LibeRation.job_progress = NULL)
  invisible(NULL)
}

#' @keywords internal
.nm_job_progress_touch <- function() {
  env <- .nm_job_progress_env()
  if (is.null(env)) {
    return(invisible(NULL))
  }
  .nm_file_write_lines(env$heartbeat_path, as.character(Sys.time()))
  invisible(NULL)
}

#' @keywords internal
.nm_format_grad_log <- function(step, grad, prefix = "") {
  header <- if (!nzchar(prefix)) {
    paste0("Gradient evaluation ", step, ":")
  } else {
    paste0(prefix, " gradient evaluation ", step, ":")
  }
  nm <- names(grad)
  parts <- vapply(seq_along(grad), function(i) {
    lab <- if (!is.null(nm) && length(nm) >= i) nm[[i]] else paste0("p", i)
    sprintf("%s=%.6g", lab, grad[[i]])
  }, character(1L))
  paste0(header, "  ", paste(parts, collapse = "  "))
}

#' @keywords internal
.nm_job_progress_log <- function(text) {
  env <- .nm_job_progress_env()
  if (is.null(env) || is.null(text) || !nzchar(text)) {
    return(invisible(NULL))
  }
  .nm_file_append(env$log_path, text)
  invisible(NULL)
}

#' @keywords internal
.nm_job_progress_event <- function(phase, detail = list(), log_msg = NULL) {
  env <- .nm_job_progress_env()
  if (is.null(env)) {
    return(invisible(NULL))
  }
  .nm_job_progress_touch()
  payload <- c(
    list(
      time = as.character(Sys.time()),
      phase = as.character(phase)[1L]
    ),
    detail
  )
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    line <- jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null")
  } else {
    line <- paste0(
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      "\t",
      phase,
      "\t",
      paste(names(detail), detail, sep = "=", collapse = " ")
    )
  }
  .nm_file_append(env$progress_path, line)
  if (!is.null(log_msg) && nzchar(log_msg)) {
    .nm_job_progress_log(log_msg)
  }
  invisible(NULL)
}

#' @keywords internal
.nm_job_progress_grad <- function(step, grad, prefix = "", print_every = 0L) {
  step <- as.integer(step[1L])
  print_every <- as.integer(print_every[1L])
  grad_num <- as.numeric(grad)
  nm <- names(grad)
  if (length(nm) == 0L || length(nm) != length(grad_num)) {
    nm <- paste0("p", seq_along(grad_num))
    names(grad_num) <- nm
  }
  env <- .nm_job_progress_env()
  should_log <- print_every > 0L && step %% print_every == 0L
  if (is.null(env)) {
    if (should_log) {
      .nm_print_grad(step, grad_num, prefix = prefix)
    }
    return(invisible(grad))
  }
  .nm_job_progress_touch()
  env$grad_step <- step
  if (should_log) {
    detail <- list(
      step = step,
      prefix = prefix,
      grad = as.list(grad_num)
    )
    log_msg <- .nm_format_grad_log(step, grad_num, prefix = prefix)
    .nm_job_progress_event("gradient", detail, log_msg = log_msg)
  } else if (step %% 10L == 0L) {
    .nm_job_progress_event(
      "gradient",
      list(step = step, prefix = prefix),
      log_msg = NULL
    )
  }
  invisible(grad)
}

#' @keywords internal
.nm_est_progress_phase <- function(method, label, detail = list(), log_msg = NULL) {
  detail$method <- method
  detail$label <- label
  if (is.null(log_msg)) {
    log_msg <- paste0("[", method, "] ", label)
  }
  .nm_job_progress_event("phase", detail, log_msg = log_msg)
}

#' @keywords internal
.nm_est_progress_outer <- function(method, iter, max_iter, objective = NA_real_) {
  obj_txt <- if (is.finite(objective)) {
    format(objective, digits = 8)
  } else {
    "NA"
  }
  .nm_est_progress_phase(
    method,
    paste0("outer ", iter, "/", max_iter),
    list(iter = iter, max_iter = max_iter, objective = objective),
    log_msg = paste0(
      method, " outer iteration ", iter, "/", max_iter,
      "  objective = ", obj_txt
    )
  )
}

#' @keywords internal
.nm_est_progress_outer_done <- function(method, iter, max_iter, objective) {
  obj_txt <- if (is.finite(objective)) {
    format(objective, digits = 8)
  } else {
    "NA"
  }
  .nm_est_progress_phase(
    method,
    paste0("outer ", iter, "/", max_iter, " complete"),
    list(iter = iter, max_iter = max_iter, objective = objective),
    log_msg = paste0(
      method, " outer iteration ", iter, "/", max_iter,
      " finished  objective = ", obj_txt
    )
  )
}

#' @keywords internal
.nm_est_progress_covariance <- function(method, cov_method = "auto") {
  .nm_est_progress_phase(
    method,
    "covariance step",
    list(cov_method = cov_method),
    log_msg = paste0(method, " covariance / standard errors (", cov_method, ")")
  )
}

#' @keywords internal
.nm_est_progress_iter <- function(method, iter, max_iter, detail = list()) {
  detail$iter <- iter
  detail$max_iter <- max_iter
  .nm_est_progress_phase(
    method,
    paste0("iteration ", iter, "/", max_iter),
    detail,
    log_msg = paste0(method, " iteration ", iter, "/", max_iter)
  )
}

#' @keywords internal
.nm_job_progress_active <- function(job_path, idle_sec = 900) {
  if (is.null(job_path) || !nzchar(job_path) || !dir.exists(job_path)) {
    return(FALSE)
  }
  now <- Sys.time()
  prog <- file.path(job_path, "worker.progress")
  hb <- file.path(job_path, "worker.heartbeat")
  log_path <- file.path(job_path, "worker.log")
  if (file.exists(prog)) {
    if (isTRUE(file.info(prog)$mtime > now - idle_sec)) {
      return(TRUE)
    }
  }
  if (file.exists(hb)) {
    if (isTRUE(file.info(hb)$mtime > now - idle_sec)) {
      return(TRUE)
    }
  }
  if (file.exists(log_path)) {
    if (isTRUE(file.info(log_path)$mtime > now - idle_sec)) {
      return(TRUE)
    }
  }
  FALSE
}

#' @keywords internal
.nm_job_progress_signature <- function(job_path) {
  prog <- file.path(job_path, "worker.progress")
  if (!file.exists(prog)) {
    return("0:0")
  }
  info <- file.info(prog)
  paste0(info$mtime, ":", info$size)
}
