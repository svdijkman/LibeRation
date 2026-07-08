#' @keywords internal
.nm_job_submit_remote_sim <- function(model,
                                      data,
                                      project,
                                      version_id,
                                      sim_id,
                                      n_sim = 1L,
                                      seed = 1L,
                                      n_cores = 1L,
                                      pk_engine = "cpp",
                                      theta = NULL,
                                      omega = NULL,
                                      sigma = NULL,
                                      label = NULL,
                                      use_fit = FALSE,
                                      est_run_id = NULL,
                                      design = NULL,
                                      vpc = FALSE,
                                      sim_compute_npc = FALSE,
                                      sim_compute_npde = FALSE,
                                      diag_n_sim = 50L,
                                      diag_refit_eta = TRUE,
                                      diag_only = FALSE,
                                      workspace_root = nm_workspace_root(),
                                      job_root = nm_job_root(),
                                      server = NULL,
                                      fit = NULL) {
  entry <- .nm_remote_get_server(server)
  sim_args <- list(
    project = project,
    version_id = version_id,
    sim_id = sim_id,
    n_sim = as.integer(n_sim),
    seed = as.integer(seed),
    n_cores = as.integer(n_cores),
    pk_engine = pk_engine,
    theta = theta,
    omega = omega,
    sigma = sigma,
    label = label,
    use_fit = isTRUE(use_fit),
    est_run_id = est_run_id,
    design = design,
    vpc = isTRUE(vpc),
    sim_compute_npc = isTRUE(sim_compute_npc),
    sim_compute_npde = isTRUE(sim_compute_npde),
    diag_n_sim = as.integer(diag_n_sim),
    diag_refit_eta = isTRUE(diag_refit_eta),
    diag_only = isTRUE(diag_only),
    workspace_root = workspace_root,
    remote = TRUE
  )
  sim_args <- sim_args[!vapply(sim_args, is.null, logical(1L))]
  if (!is.null(fit)) {
    sim_args$fit <- fit
  }
  payload <- list(
    job_type = "sim",
    label = label,
    model_b64 = .nm_rds_b64(model),
    data_b64 = .nm_rds_b64(data),
    est_args = sim_args
  )
  resp <- .nm_remote_http("POST", "/v1/jobs", server = entry$id, body = payload)
  remote <- resp$job
  job_id <- paste0("remote_", remote$id)
  job_path <- .nm_job_path(job_id, job_root)
  dir.create(job_path, recursive = TRUE, showWarnings = FALSE)

  meta <- list(
    id = job_id,
    remote = TRUE,
    remote_job_id = remote$id,
    server_id = entry$id,
    server_name = entry$name %||% entry$id,
    label = remote$label %||% label,
    job_type = "sim",
    status = remote$status %||% "queued",
    method = "SIM",
    sim_id = sim_id,
    version_id = version_id,
    created = remote$created %||% as.character(Sys.time()),
    started = remote$started %||% "",
    finished = remote$finished %||% "",
    pid = NA_integer_,
    objective = NA_real_,
    error = "",
    warnings = character()
  )
  saveRDS(meta, file.path(job_path, "meta.rds"))
  saveRDS(entry, file.path(job_path, "remote_server.rds"))

  structure(
    list(
      id = job_id,
      path = job_path,
      process = NULL,
      remote = TRUE,
      status = meta$status
    ),
    class = "nm_job_handle"
  )
}

#' Submit a job to a remote scheduler
#'
#' @param model,method,data,... Passed to estimation; see \code{\link{nm_job_submit}}.
#' @param server Remote server id (NULL = default).
#' @param data_ref Optional list with \code{dataset_id} and \code{md5} for cluster data.
#' @param label Job label.
#' @param job_root Local stub directory for tracking.
#' @return Job handle (class \code{nm_job_handle}).
#' @keywords internal
.nm_job_submit_remote <- function(model,
                                  data = NULL,
                                  method = "FO",
                                  label = NULL,
                                  server = NULL,
                                  data_ref = NULL,
                                  job_root = nm_job_root(),
                                  workspace_project = NULL,
                                  workspace_version_id = NULL,
                                  est_run_id = NULL,
                                  workspace_root = NULL,
                                  ...) {
  entry <- .nm_remote_get_server(server)
  est_args <- list(...)
  payload <- list(
    job_type = "est",
    method = method,
    label = label,
    model_b64 = .nm_rds_b64(model),
    est_args = est_args
  )
  if (!is.null(data_ref)) {
    payload$data_ref <- data_ref
  } else if (!is.null(data)) {
    payload$data_b64 <- .nm_rds_b64(data)
    payload$data_md5 <- .nm_data_md5(data)
  } else {
    stop("data or data_ref is required for remote estimation.", call. = FALSE)
  }

  resp <- .nm_remote_http("POST", "/v1/jobs", server = entry$id, body = payload)
  remote <- resp$job
  job_id <- paste0("remote_", remote$id)
  job_path <- .nm_job_path(job_id, job_root)
  dir.create(job_path, recursive = TRUE, showWarnings = FALSE)

  meta <- list(
    id = job_id,
    remote = TRUE,
    remote_job_id = remote$id,
    server_id = entry$id,
    server_name = entry$name %||% entry$id,
    label = remote$label %||% label,
    job_type = "est",
    status = remote$status %||% "queued",
    method = method,
    created = remote$created %||% as.character(Sys.time()),
    started = remote$started %||% "",
    finished = remote$finished %||% "",
    pid = NA_integer_,
    objective = NA_real_,
    error = "",
    warnings = character(),
    workspace_project = workspace_project %||% "",
    workspace_version_id = workspace_version_id %||% "",
    est_run_id = est_run_id %||% "",
    workspace_root = workspace_root %||% "",
    workspace_saved = FALSE
  )
  saveRDS(meta, file.path(job_path, "meta.rds"))
  saveRDS(entry, file.path(job_path, "remote_server.rds"))

  structure(
    list(
      id = job_id,
      path = job_path,
      process = NULL,
      remote = TRUE,
      status = meta$status
    ),
    class = "nm_job_handle"
  )
}

#' @keywords internal
.nm_job_is_remote <- function(job_id, job_root = nm_job_root()) {
  meta <- .nm_job_read_meta(job_id, job_root)
  !is.null(meta) && isTRUE(meta$remote)
}

#' @keywords internal
.nm_job_remote_entry <- function(job_id, job_root = nm_job_root()) {
  path <- file.path(.nm_job_path(job_id, job_root), "remote_server.rds")
  if (!file.exists(path)) {
    entry <- .nm_remote_get_server(NULL)
    return(entry)
  }
  readRDS(path)
}

#' @keywords internal
.nm_job_remote_sync <- function(job_id, job_root = nm_job_root()) {
  meta <- .nm_job_read_meta(job_id, job_root)
  if (is.null(meta) || !isTRUE(meta$remote)) {
    return(meta)
  }
  entry <- .nm_job_remote_entry(job_id, job_root)
  rid <- meta$remote_job_id
  resp <- tryCatch(
    .nm_remote_http("GET", paste0("/v1/jobs/", rid), server = entry$id, timeout = 30L),
    error = function(e) {
      meta$remote_sync_error <- conditionMessage(e)
      return(NULL)
    }
  )
  if (is.null(resp)) {
    return(meta)
  }
  remote <- resp$job
  if (!is.null(resp$warnings) && length(resp$warnings) > 0L) {
    meta$warnings <- unique(c(meta$warnings %||% character(), resp$warnings))
  }
  for (fld in c("status", "started", "finished", "objective", "error")) {
    if (!is.null(remote[[fld]])) {
      meta[[fld]] <- remote[[fld]]
    }
  }
  saveRDS(meta, .nm_job_meta_path(job_id, job_root))
  meta
}

#' @keywords internal
.nm_job_status_remote <- function(job_id, job_root = nm_job_root()) {
  .nm_job_remote_sync(job_id, job_root)
}

#' @keywords internal
.nm_job_result_remote <- function(job_id, job_root = nm_job_root()) {
  meta <- .nm_job_remote_sync(job_id, job_root)
  if (!identical(meta$status, "success")) {
    stop("Remote job not successful: ", meta$status, call. = FALSE)
  }
  entry <- .nm_job_remote_entry(job_id, job_root)
  resp <- .nm_remote_http(
    "GET",
    paste0("/v1/jobs/", meta$remote_job_id, "/result"),
    server = entry$id,
    timeout = 600L
  )
  .nm_b64_rds(resp$result_b64)
}

#' @keywords internal
.nm_job_cache_remote_log <- function(stub_id, server_id, remote_id,
                                     job_root = nm_job_root(), tail = 500L) {
  if (is.null(remote_id) || !nzchar(remote_id)) {
    return(invisible(FALSE))
  }
  stub_path <- .nm_job_path(stub_id, job_root)
  if (!dir.exists(stub_path)) {
    return(invisible(FALSE))
  }
  log_txt <- tryCatch({
    resp <- .nm_remote_http(
      "GET",
      paste0("/v1/jobs/", remote_id, "/log?tail=", as.integer(tail)),
      server = server_id
    )
    resp$log %||% ""
  }, error = function(e) "")
  if (!nzchar(log_txt)) {
    return(invisible(FALSE))
  }
  log_path <- file.path(stub_path, "worker.log")
  lines <- strsplit(log_txt, "\n", fixed = TRUE)[[1L]]
  writeLines(lines, log_path, useBytes = TRUE)
  invisible(TRUE)
}

#' @keywords internal
.nm_job_log_remote <- function(job_id, tail = 100L, job_root = nm_job_root()) {
  meta <- .nm_job_read_meta(job_id, job_root)
  if (is.null(meta) || !isTRUE(meta$remote)) {
    return("")
  }
  cached_path <- file.path(.nm_job_path(job_id, job_root), "worker.log")
  log_txt <- tryCatch({
    entry <- .nm_job_remote_entry(job_id, job_root)
    resp <- .nm_remote_http(
      "GET",
      paste0("/v1/jobs/", meta$remote_job_id, "/log?tail=", as.integer(tail)),
      server = entry$id
    )
    resp$log %||% ""
  }, error = function(e) "")
  if (nzchar(log_txt)) {
    if (meta$status %in% c("error", "success", "cancelled")) {
      tryCatch({
        lines <- strsplit(log_txt, "\n", fixed = TRUE)[[1L]]
        writeLines(lines, cached_path, useBytes = TRUE)
      }, error = function(e) NULL)
    }
    return(log_txt)
  }
  if (file.exists(cached_path)) {
    lines <- readLines(cached_path, warn = FALSE)
    if (length(lines) <= tail) {
      return(paste(lines, collapse = "\n"))
    }
    return(paste(tail(lines, tail), collapse = "\n"))
  }
  ""
}

#' @keywords internal
.nm_job_cancel_remote <- function(job_id, job_root = nm_job_root()) {
  meta <- .nm_job_read_meta(job_id, job_root)
  if (is.null(meta) || !isTRUE(meta$remote)) {
    return(NULL)
  }
  entry <- .nm_job_remote_entry(job_id, job_root)
  resp <- .nm_remote_http(
    "DELETE",
    paste0("/v1/jobs/", meta$remote_job_id),
    server = entry$id
  )
  meta <- resp$job
  meta$id <- job_id
  meta$remote <- TRUE
  meta$remote_job_id <- resp$job$id
  saveRDS(meta, .nm_job_meta_path(job_id, job_root))
  meta
}

#' @keywords internal
.nm_remote_job_empty_df <- function() {
  data.frame(
    id = character(),
    label = character(),
    job_type = character(),
    status = character(),
    method = character(),
    sim_id = character(),
    n_sim = integer(),
    created = character(),
    started = character(),
    finished = character(),
    objective = numeric(),
    error = character(),
    server = character(),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.nm_remote_jobs_to_df <- function(jobs, server_id = NULL) {
  if (is.null(jobs) || length(jobs) == 0L) {
    return(.nm_remote_job_empty_df())
  }
  if (is.data.frame(jobs)) {
    df <- jobs
  } else {
    df <- do.call(rbind, lapply(jobs, function(x) {
      as.data.frame(x, stringsAsFactors = FALSE)
    }))
  }
  if (nrow(df) == 0L) {
    return(.nm_remote_job_empty_df())
  }
  for (col in c("sim_id", "n_sim", "error")) {
    if (!col %in% names(df)) {
      df[[col]] <- if (col == "n_sim") NA_integer_ else ""
    }
  }
  srv_label <- server_id
  if (!is.null(server_id) && nzchar(server_id)) {
    srv_label <- tryCatch(
      {
        entry <- .nm_remote_get_server(server_id)
        as.character(entry$name %||% server_id)
      },
      error = function(e) server_id
    )
  }
  df$server <- srv_label
  df[order(df$created, decreasing = TRUE), , drop = FALSE]
}

#' @keywords internal
.nm_job_sync_remote_stubs <- function(server_id, jobs_df, job_root = nm_job_root()) {
  if (is.null(jobs_df) || nrow(jobs_df) == 0L) {
    return(invisible(0L))
  }
  n <- 0L
  for (i in seq_len(nrow(jobs_df))) {
    rid <- as.character(jobs_df$id[[i]])
    stub_id <- paste0("remote_", rid)
    meta_path <- .nm_job_meta_path(stub_id, job_root)
    if (!file.exists(meta_path)) {
      next
    }
    meta <- .nm_job_read_meta(stub_id, job_root)
    if (is.null(meta)) {
      next
    }
    changed <- FALSE
    for (fld in c("status", "started", "finished", "objective", "error", "label", "method")) {
      if (!fld %in% names(jobs_df)) {
        next
      }
      val <- jobs_df[[fld]][[i]]
      if (is.null(val) || (length(val) == 1L && is.na(val))) {
        next
      }
      if (!identical(meta[[fld]], val)) {
        meta[[fld]] <- val
        changed <- TRUE
      }
    }
    if (changed) {
      saveRDS(meta, meta_path)
      n <- n + 1L
    }
    if (meta$status %in% c("error", "success", "cancelled")) {
      .nm_job_cache_remote_log(stub_id, server_id, rid, job_root = job_root)
    }
  }
  invisible(n)
}

#' List jobs on a remote scheduler
#'
#' @param server Server id from \code{\link{nm_remote_server_list}}.
#' @param timeout HTTP timeout in seconds. Background pollers (the Shiny job hub)
#'   pass a short value so a slow/overloaded server cannot stall the poller.
#' @return Job data frame (same columns as \code{\link{nm_job_list}}).
#' @export
nm_remote_job_list <- function(server, timeout = 30L) {
  resp <- tryCatch(
    .nm_remote_http("GET", "/v1/jobs", server = server, timeout = timeout),
    error = function(e) {
      df <- .nm_remote_job_empty_df()
      attr(df, "error") <- conditionMessage(e)
      return(df)
    }
  )
  if (is.data.frame(resp) && !is.null(attr(resp, "error"))) {
    return(resp)
  }
  .nm_remote_jobs_to_df(resp$jobs, server_id = server)
}

#' Remote job status
#'
#' @param server Server id.
#' @param job_id Remote job id.
#' @return Job metadata list.
#' @keywords internal
.nm_remote_job_status <- function(server, job_id) {
  resp <- .nm_remote_http("GET", paste0("/v1/jobs/", job_id), server = server)
  resp$job
}

#' Remote worker log
#'
#' @param server Server id.
#' @param job_id Remote job id.
#' @param tail Number of log lines.
#' @keywords internal
.nm_remote_job_log <- function(server, job_id, tail = 100L) {
  resp <- .nm_remote_http(
    "GET",
    paste0("/v1/jobs/", job_id, "/log?tail=", as.integer(tail)),
    server = server
  )
  resp$log %||% ""
}

#' @keywords internal
#' @param timeout Per-request HTTP timeout in seconds. Kept short for the
#'   background poller so one slow log fetch cannot stall the whole hub.
.nm_job_remote_log_signature <- function(job_root, server_ids, timeout = 8L) {
  if (is.null(server_ids) || length(server_ids) == 0L || !dir.exists(job_root)) {
    return("")
  }
  stubs <- list.dirs(job_root, full.names = FALSE, recursive = FALSE)
  stubs <- stubs[grepl("^remote_", stubs)]
  if (length(stubs) == 0L) {
    return("")
  }
  parts <- character()
  for (sid in unique(as.character(server_ids))) {
    if (!nzchar(sid)) {
      next
    }
    for (stub in stubs) {
      meta <- .nm_job_read_meta(stub, job_root)
      if (is.null(meta) || !isTRUE(meta$remote)) {
        next
      }
      if (!identical(as.character(meta$server_id), sid)) {
        next
      }
      # Only running jobs have a growing log worth polling. Queued jobs have no
      # progress yet; their queued->running transition is detected via the job
      # list + stub sync, so polling their log every tick just adds blocking
      # round-trips (which is what stalled the hub during the burst stress test).
      if (!identical(meta$status, "running")) {
        next
      }
      rid <- meta$remote_job_id
      if (is.null(rid) || !nzchar(rid)) {
        next
      }
      resp <- tryCatch(
        .nm_remote_http(
          "GET",
          paste0("/v1/jobs/", rid, "/log?tail=12"),
          server = sid,
          timeout = timeout
        ),
        error = function(e) NULL
      )
      log_txt <- if (is.null(resp)) "" else as.character(resp$log %||% "")
      log_tail <- if (nchar(log_txt) > 120L) {
        substr(log_txt, nchar(log_txt) - 119L, nchar(log_txt))
      } else {
        log_txt
      }
      parts <- c(
        parts,
        paste0(stub, ":", nchar(log_txt), ":", log_tail)
      )
    }
  }
  paste(parts, collapse = "|")
}

#' @keywords internal
.nm_remote_job_cleanup <- function(server, job_root = nm_job_root()) {
  n <- 0L
  resp <- tryCatch(
    .nm_remote_http("POST", "/v1/jobs/cleanup", server = server, body = list()),
    error = function(e) NULL
  )
  if (!is.null(resp)) {
    n <- as.integer(resp$removed %||% 0L)
  }
  remaining <- tryCatch(nm_remote_job_list(server), error = function(e) NULL)
  remaining_ids <- if (is.null(remaining) || nrow(remaining) == 0L) {
    character()
  } else {
    as.character(remaining$id)
  }
  if (!dir.exists(job_root)) {
    return(n)
  }
  ids <- list.dirs(job_root, full.names = FALSE, recursive = FALSE)
  ids <- ids[grepl("^remote_", ids)]
  for (stub in ids) {
    meta_path <- file.path(job_root, stub, "meta.rds")
    if (!file.exists(meta_path)) {
      next
    }
    meta <- tryCatch(readRDS(meta_path), error = function(e) NULL)
    if (is.null(meta) || !identical(meta$server_id, server)) {
      next
    }
    rid <- sub("^remote_", "", stub, fixed = TRUE)
    if (!rid %in% remaining_ids || meta$status %in% c("success", "error", "cancelled")) {
      unlink(file.path(job_root, stub), recursive = TRUE)
    }
  }
  n
}

#' Cancel a remote job
#'
#' @param server Server id.
#' @param job_id Remote job id.
#' @keywords internal
.nm_remote_job_cancel <- function(server, job_id) {
  resp <- .nm_remote_http("DELETE", paste0("/v1/jobs/", job_id), server = server)
  resp$job
}

#' Fetch a remote job result
#'
#' @param server Server id.
#' @param job_id Remote job id.
#' @keywords internal
.nm_remote_job_result <- function(server, job_id) {
  resp <- .nm_remote_http(
    "GET",
    paste0("/v1/jobs/", job_id, "/result"),
    server = server,
    timeout = 600L
  )
  .nm_b64_rds(resp$result_b64)
}

#' List datasets on a remote cluster
#'
#' @param server Server id or NULL for default.
#' @return Data frame of datasets with id, label, md5.
#' @export
nm_remote_dataset_list <- function(server = NULL) {
  resp <- .nm_remote_http("GET", "/v1/datasets", server = server)
  resp$datasets
}
