#' Root directory for async estimation jobs
#'
#' Defaults to a per-user cache directory. Override with
#' \code{options(LibeRation.job_dir = "/path/to/jobs")}.
#'
#' @return Character path (created if missing).
#' @examples NULL
#' @export
nm_job_root <- function() {
  root <- getOption("LibeRation.job_dir")
  if (is.null(root)) {
    root <- if (getRversion() >= "4.0") {
      file.path(tools::R_user_dir("LibeRation", which = "cache"), "jobs")
    } else {
      file.path(path.expand("~"), ".LibeRation", "jobs")
    }
  }
  if (!dir.exists(root)) {
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(root, winslash = "/", mustWork = FALSE)
}

#' @keywords internal
.nm_job_path <- function(job_id, root = nm_job_root()) {
  file.path(root, job_id)
}

#' @keywords internal
.nm_job_meta_path <- function(job_id, root = nm_job_root()) {
  file.path(.nm_job_path(job_id, root), "meta.rds")
}

#' @keywords internal
.nm_job_read_meta <- function(job_id, root = nm_job_root()) {
  path <- .nm_job_meta_path(job_id, root)
  if (!file.exists(path)) {
    return(NULL)
  }
  .nm_read_rds_safe(path)
}

#' @keywords internal
.nm_job_write_meta <- function(meta, root = nm_job_root()) {
  path <- .nm_job_meta_path(meta$id, root)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(meta, path)
  invisible(meta)
}

#' @keywords internal
.nm_job_new_id <- function() {
  paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", substr(digest::digest(Sys.time()), 1L, 6L))
}

#' @keywords internal
.nm_job_new_id_fallback <- function() {
  paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", sample.int(999999L, 1L))
}

#' @keywords internal
.nm_job_make_id <- function() {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(.nm_job_new_id())
  }
  .nm_job_new_id_fallback()
}

#' @keywords internal
.nm_pkg_source_root <- function(pkg = "LibeRation") {
  if (!pkg %in% loadedNamespaces()) {
    return("")
  }
  ip <- system.file("", package = pkg)
  if (!nzchar(ip)) {
    return("")
  }
  ip <- normalizePath(ip, winslash = "/", mustWork = FALSE)
  if (file.exists(file.path(ip, "DESCRIPTION"))) {
    return(ip)
  }
  parent <- normalizePath(file.path(ip, ".."), winslash = "/", mustWork = FALSE)
  if (file.exists(file.path(parent, "DESCRIPTION"))) {
    return(parent)
  }
  ip
}

#' @keywords internal
.nm_job_dev_env <- function() {
  snap <- getOption("LibeRation.job_dev_env")
  if (is.list(snap) && !is.null(snap$mode)) {
    return(snap)
  }

  nm_root <- ""
  ad_root <- Sys.getenv("LibeRation_LIBERTAD_ROOT", "")

  if ("LibeRation" %in% loadedNamespaces()) {
    nm_root <- .nm_pkg_source_root("LibeRation")
  }
  if (!nzchar(nm_root)) {
    nm_root <- Sys.getenv("LibeRation_ROOT", "")
    if (nzchar(nm_root)) {
      nm_root <- normalizePath(nm_root, winslash = "/", mustWork = FALSE)
    }
  }
  if (!nzchar(ad_root) && "LibeRtAD" %in% loadedNamespaces()) {
    ad_root <- .nm_pkg_source_root("LibeRtAD")
  }
  if (!nzchar(ad_root) && nzchar(nm_root)) {
    for (cand in c(
      file.path(nm_root, "..", "LibeRtAD"),
      file.path(dirname(nm_root), "LibeRtAD")
    )) {
      if (dir.exists(cand)) {
        ad_root <- normalizePath(cand, winslash = "/", mustWork = FALSE)
        break
      }
    }
  }

  is_dev <- FALSE
  if (nzchar(nm_root)) {
    if (requireNamespace("pkgload", quietly = TRUE)) {
      is_dev <- isTRUE(pkgload::is_dev_package("LibeRation"))
    }
    if (!is_dev) {
      is_dev <- !file.exists(file.path(nm_root, "Meta", "package.rds"))
    }
  }

  if (is_dev && nzchar(nm_root)) {
    return(list(mode = "dev", nm_root = nm_root, ad_root = ad_root))
  }
  if (requireNamespace("LibeRation", quietly = TRUE)) {
    return(list(mode = "installed", nm_root = nm_root, ad_root = ad_root))
  }
  if (nzchar(nm_root)) {
    return(list(mode = "dev", nm_root = nm_root, ad_root = ad_root))
  }
  list(mode = "installed", nm_root = nm_root, ad_root = ad_root)
}

#' @keywords internal
.nm_job_worker_body <- function(job_path) {
  .nm_job_worker_impl(job_path)
}

#' @keywords internal
.nm_job_run_worker <- function(job_path) {
  .nm_job_worker_impl(job_path)
}

#' @keywords internal
.nm_job_worker_load <- function(dev_env, log_path = NULL) {
  log <- function(...) {
    if (!is.null(log_path)) {
      .nm_file_append(log_path, paste0(...))
    }
  }
  if (identical(dev_env$mode, "dev") && !is.null(dev_env$nm_root) && nzchar(dev_env$nm_root)) {
    log("Loading dev packages from: ", dev_env$nm_root, "\n")
    if (!requireNamespace("pkgload", quietly = TRUE)) {
      stop("Package 'pkgload' is required to run jobs from a dev load_all session.")
    }
    if (!is.null(dev_env$ad_root) && nzchar(dev_env$ad_root)) {
      pkgload::load_all(dev_env$ad_root, quiet = TRUE, compile = FALSE)
    }
    pkgload::load_all(dev_env$nm_root, quiet = TRUE, recompile = FALSE)
    return(invisible(NULL))
  }
  log("Loading installed LibeRation\n")
  pkg <- "LibeRation"
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package LibeRation is not installed.")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  invisible(NULL)
}

#' @keywords internal
.nm_job_clean_text <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return("")
  }
  txt <- as.character(x)[1L]
  if (is.na(txt) || !nzchar(txt) || identical(txt, "NA")) {
    ""
  } else {
    txt
  }
}

#' @keywords internal
.nm_job_worker_load_env <- function(dev_env, log_path = NULL) {
  log <- function(...) {
    if (!is.null(log_path)) {
      .nm_file_append(log_path, paste0(...))
    }
  }
  if ("LibeRation" %in% loadedNamespaces() &&
      exists("nm_est", envir = asNamespace("LibeRation"), inherits = FALSE)) {
    log("LibeRation already loaded (", system.file("", package = "LibeRation"), ")\n")
    return(invisible(TRUE))
  }
  nm_root <- if (is.null(dev_env$nm_root)) "" else as.character(dev_env$nm_root)
  ad_root <- if (is.null(dev_env$ad_root)) "" else as.character(dev_env$ad_root)
  is_dev_root <- function(root) {
    nzchar(root) && dir.exists(root) &&
      file.exists(file.path(root, "DESCRIPTION")) &&
      !file.exists(file.path(root, "Meta", "package.rds"))
  }
  load_root <- function(root, pkg) {
    if (!nzchar(root)) {
      return(invisible(FALSE))
    }
    if (is_dev_root(root)) {
      if (!requireNamespace("pkgload", quietly = TRUE)) {
        stop("Package 'pkgload' is required to run jobs from a dev load_all session.", call. = FALSE)
      }
      log("Loading dev ", pkg, " from: ", root, "\n")
      pkgload::load_all(root, quiet = TRUE, compile = FALSE, recompile = FALSE)
      return(invisible(TRUE))
    }
    if (requireNamespace(pkg, quietly = TRUE)) {
      log("Loading installed ", pkg, "\n")
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
      return(invisible(TRUE))
    }
    stop("Cannot load ", pkg, ".", call. = FALSE)
  }
  if (identical(dev_env$mode, "dev") && nzchar(nm_root)) {
    load_root(ad_root, "LibeRtAD")
    load_root(nm_root, "LibeRation")
    log("LibeRation loaded from dev: ", nm_root, "\n")
    return(invisible(TRUE))
  }
  if (!requireNamespace("LibeRation", quietly = TRUE)) {
    stop("Package LibeRation is not installed.", call. = FALSE)
  }
  log("Loading installed package LibeRation (", system.file("", package = "LibeRation"), ")\n")
  tryCatch(
    suppressPackageStartupMessages(library(LibeRation)),
    error = function(e) {
      stop("Failed to load LibeRation: ", conditionMessage(e), call. = FALSE)
    }
  )
  log("LibeRtAD path: ", system.file("", package = "LibeRtAD"), "\n")
  invisible(TRUE)
}

#' @keywords internal
.nm_job_worker_impl <- function(job_path) {
  meta_path <- file.path(job_path, "meta.rds")
  log_path <- file.path(job_path, "worker.log")
  meta <- .nm_read_rds_safe(meta_path)
  meta$status <- "running"
  meta$started <- as.character(Sys.time())
  .nm_save_rds_safe(meta, meta_path)
  .nm_file_write_lines(file.path(job_path, "worker.pid"), as.character(Sys.getpid()))
  .nm_file_write_lines(file.path(job_path, "worker.heartbeat"), as.character(Sys.time()))
  job_type <- if (is.null(meta$job_type)) "est" else meta$job_type
  .nm_file_append(
    log_path,
    paste0("LibeRation worker\nJob started: ", meta$id, " ( ", job_type, " )")
  )

  dev_env <- .nm_read_rds_safe(file.path(job_path, "env.rds"))
  options(LibeRation.job_dev_env = dev_env)
  .nm_job_worker_load_env(dev_env, log_path = log_path)

  args <- .nm_job_read_args(job_path)
  err <- NULL
  if (identical(job_type, "sim")) {
    nm_simulate <- getExportedValue("LibeRation", "nm_simulate")
    nm_workspace_save_sim <- getExportedValue("LibeRation", "nm_workspace_save_sim")
    compute_npc <- isTRUE(args$sim_compute_npc)
    compute_npde <- isTRUE(args$sim_compute_npde)
    diag_n_sim <- as.integer(args$diag_n_sim %||% 50L)
    diag_refit_eta <- isTRUE(args$diag_refit_eta)
    run_diag <- (compute_npc || compute_npde) && isTRUE(args$use_fit)
    diag_only <- run_diag && !isTRUE(args$vpc) && as.integer(args$n_sim) <= 1L
    run_replicate_sim <- !diag_only
    meta <- .nm_read_rds_safe(meta_path)
    sim_ok <- FALSE
    sim_data <- NULL
    fit_diag <- NULL
    err <- NULL

    if (run_replicate_sim) {
      sim_args <- args[c(
        "model", "data", "n_sim", "seed", "n_cores", "pk_engine",
        "theta", "omega", "sigma", "design", "vpc"
      )]
      sim_args <- sim_args[!vapply(sim_args, is.null, logical(1L))]
      if (!is.null(args$design)) {
        tmpl <- LibeRation::nm_sim_template_data(args$model, args$data, args$design)
        sim_args$data <- tmpl
        sim_args$model <- args$model
      }
      sim_args$design <- NULL
      n_cores <- as.integer(sim_args$n_cores %||% 1L)
      max_cpu <- meta$limits$max_cpu %||% NA_integer_
      if (!is.na(max_cpu) && max_cpu > 0L) {
        n_cores <- min(n_cores, as.integer(max_cpu))
      }
      sim_args$n_cores <- max(1L, n_cores)
      if (identical(dev_env$mode, "dev") && n_cores > 1L) {
        .nm_file_append(
          log_path,
          "Parallel simulation replicates run sequentially in development package sessions."
        )
      }
      sim_out <- tryCatch(
        do.call(nm_simulate, sim_args),
        error = function(e) {
          err <<- conditionMessage(e)
          NULL
        }
      )
      if (is.null(sim_out)) {
        meta$status <- "error"
        meta$error <- err
        .nm_file_write_lines(file.path(job_path, "error.txt"), err)
        .nm_file_append(log_path, paste0("Simulation failed: ", err))
      } else {
        sim_data <- LibeRation:::.nm_sim_pack_output(sim_out)
        if (!isTRUE(args$remote)) {
          do.call(
            nm_workspace_save_sim,
            c(
              list(
                project = args$project,
                version_id = args$version_id,
                sim_id = args$sim_id,
                sim_data = sim_data,
                root = args$workspace_root,
                label = args$label,
                seed = args$seed,
                n_sim = args$n_sim,
                use_fit = args$use_fit,
                est_run_id = args$est_run_id,
                vpc = isTRUE(args$vpc)
              )
            )
          )
        }
        sim_ok <- TRUE
        .nm_file_append(
          log_path,
          paste0("Simulation completed: ", args$sim_id, " replicates=", args$n_sim)
        )
      }
    } else if (run_diag) {
      sim_ok <- TRUE
      .nm_file_append(
        log_path,
        "Diagnostic-only job (NPC/NPDE); skipping workspace simulation save."
      )
    } else {
      meta$status <- "error"
      meta$error <- "Nothing to run: enable VPC, replications, NPC, or NPDE."
      .nm_file_write_lines(file.path(job_path, "error.txt"), meta$error)
      .nm_file_append(log_path, paste0("Job failed: ", meta$error))
    }

    if (sim_ok && run_diag) {
      nm_workspace_load_run_fit <- getExportedValue("LibeRation", "nm_workspace_load_run_fit")
      nm_workspace_save_run <- getExportedValue("LibeRation", "nm_workspace_save_run")
      nm_workspace_list_runs <- getExportedValue("LibeRation", "nm_workspace_list_runs")
      nm_add_npc_npde <- getExportedValue("LibeRation", "nm_add_npc_npde")
      est_run_id <- args$est_run_id
      if ((is.null(est_run_id) || !nzchar(est_run_id)) && !isTRUE(args$remote)) {
        runs <- nm_workspace_list_runs(
          args$project, args$version_id, root = args$workspace_root
        )
        if (nrow(runs) > 0L) {
          est_run_id <- runs$run_id[[1L]]
        }
      }
      if (!is.null(est_run_id) && nzchar(est_run_id)) {
        what <- c(
          if (compute_npc) "NPC",
          if (compute_npde) "NPDE"
        )
        .nm_file_append(
          log_path,
          paste0(paste(what, collapse = "/"), " on fit: ", est_run_id,
                 " - ", diag_n_sim, " simulations")
        )
        fit_diag <- tryCatch(
          if (isTRUE(args$remote) && !is.null(args$fit)) {
            args$fit
          } else {
            nm_workspace_load_run_fit(
              args$project, args$version_id, est_run_id, root = args$workspace_root
            )
          },
          error = function(e) NULL
        )
        if (!is.null(fit_diag) && !identical(fit_diag$method, "BAYES")) {
          pk_diag <- fit_diag$pk_engine %||% args$pk_engine %||% "cpp"
          fit_diag <- tryCatch(
            nm_add_npc_npde(
              fit_diag,
              n_sim = diag_n_sim,
              seed = as.integer(args$seed %||% 1L),
              refit_eta = diag_refit_eta,
              compute_npc = compute_npc,
              compute_npde = compute_npde,
              pk_engine = pk_diag,
              n_cores = 1L
            ),
            error = function(e) {
              .nm_file_append(log_path, paste0("NPC/NPDE failed: ", conditionMessage(e)))
              err <<- conditionMessage(e)
              NULL
            }
          )
          if (!is.null(fit_diag)) {
            if (!isTRUE(args$remote)) {
              nm_workspace_save_run(
                args$project,
                args$version_id,
                est_run_id,
                fit_diag,
                root = args$workspace_root,
                job_id = meta$id
              )
            }
            n_ok <- fit_diag$npc_npde$n_ok %||% NA_integer_
            .nm_file_append(
              log_path,
              paste0(paste(what, collapse = "/"), " completed: ", n_ok,
                     " obs rows / ", diag_n_sim, " simulations")
            )
            meta$diag_n_sim <- diag_n_sim
            meta$diag_ok <- n_ok
            meta$diag_est_run <- est_run_id
            meta$diag_npc <- compute_npc
            meta$diag_npde <- compute_npde
          } else {
            sim_ok <- FALSE
            meta$status <- "error"
            meta$error <- err %||% "NPC/NPDE failed."
            .nm_file_write_lines(file.path(job_path, "error.txt"), meta$error)
          }
        } else if (!is.null(fit_diag) && identical(fit_diag$method, "BAYES")) {
          .nm_file_append(log_path, "NPC/NPDE skipped: BAYES fit")
        }
      } else {
        .nm_file_append(log_path, "NPC/NPDE skipped: no estimation run for fit")
      }
    }

    if (identical(meta$status, "queued") || identical(meta$status, "running")) {
      if (sim_ok) {
        result <- list(
          sim_data = sim_data,
          sim_id = args$sim_id,
          project = args$project,
          version_id = args$version_id,
          label = args$label,
          seed = args$seed,
          n_sim = args$n_sim,
          use_fit = args$use_fit,
          est_run_id = args$est_run_id,
          vpc = isTRUE(args$vpc)
        )
        if (!is.null(meta$diag_est_run)) {
          result$diag_est_run <- meta$diag_est_run
        }
        if (!is.null(fit_diag)) {
          result$fit_diag <- fit_diag
        }
        .nm_job_write_result(result, job_path)
        meta$status <- "success"
        meta$sim_id <- args$sim_id
        meta$version_id <- args$version_id
        meta$error <- ""
      }
    }
  } else {
    nm_est_fun <- getExportedValue("LibeRation", "nm_est")
    est_args <- args
    est_args$task <- NULL
    bootstrap_n <- as.integer(est_args$bootstrap_n %||% 0L)
    bootstrap_seed <- as.integer(est_args$bootstrap_seed %||% 1L)
    est_args$bootstrap_n <- NULL
    est_args$bootstrap_seed <- NULL
    boot_dots <- list()
    if (!is.null(est_args$bootstrap_control)) {
      boot_dots$control <- est_args$bootstrap_control
      est_args$bootstrap_control <- NULL
    }
    for (nm in c("backend", "grad", "pk_engine", "engine", "max_outer", "n_iter", "n_burn",
                 "n_mcmc", "n_quad", "n_imp")) {
      if (!is.null(est_args[[nm]])) {
        boot_dots[[nm]] <- est_args[[nm]]
      }
    }
    cat(
      "Running estimation:", est_args$method %||% "FO", "\n",
      file = log_path, append = TRUE
    )
    .nm_job_progress_init(job_path)
    on.exit(.nm_job_progress_clear(), add = TRUE)
    if (is.null(est_args$control)) {
      est_args$control <- list()
    }
    if (is.null(est_args$control$print_grad_every)) {
      est_args$control$print_grad_every <- 1L
    }
    .nm_job_progress_event(
      "estimation_start",
      list(method = est_args$method %||% "FO"),
      log_msg = paste("Running estimation:", est_args$method %||% "FO")
    )
    fit <- tryCatch(
      do.call(nm_est_fun, est_args),
      error = function(e) {
        err <<- conditionMessage(e)
        NULL
      }
    )
    meta <- .nm_read_rds_safe(meta_path)
    if (is.null(fit)) {
      meta$status <- "error"
      meta$error <- err
      .nm_file_write_lines(file.path(job_path, "error.txt"), err)
      .nm_file_append(log_path, paste0("Job failed: ", err))
    } else {
      if (bootstrap_n > 0L && !identical(fit$method, "BAYES")) {
        .nm_file_append(
          log_path,
          paste0("Bootstrap SE: ", bootstrap_n, " replicates (seed=", bootstrap_seed, ")")
        )
        fit <- tryCatch(
          do.call(
            LibeRation:::.nm_bootstrap_attach,
            c(
              list(fit = fit, n_boot = bootstrap_n, seed = bootstrap_seed),
              boot_dots
            )
          ),
          error = function(e) {
            .nm_file_append(log_path, paste0("Bootstrap failed: ", conditionMessage(e)))
            fit
          }
        )
        if (!is.null(fit$bootstrap)) {
          .nm_file_append(
            log_path,
            paste0("Bootstrap completed: ", fit$bootstrap$n_ok, " / ", bootstrap_n)
          )
          meta$bootstrap_n <- bootstrap_n
          meta$bootstrap_ok <- fit$bootstrap$n_ok
        }
      }
      .nm_job_write_result(fit, job_path)
      meta$status <- "success"
      meta$objective <- fit$objective
      meta$method <- fit$method
      meta$error <- ""
      .nm_file_append(log_path, paste0("Job completed. objective = ", fit$objective))
      .nm_job_progress_event(
        "estimation_done",
        list(method = fit$method, objective = fit$objective),
        log_msg = paste("Job completed. objective =", fit$objective)
      )
    }
  }
  meta$finished <- as.character(Sys.time())
  .nm_save_rds_safe(meta, meta_path)
  invisible(NULL)
}

#' Submit an asynchronous \code{nm_est} job
#'
#' Runs estimation in a background R process via \pkg{callr}. Poll with
#' \code{\link{nm_job_status}} and retrieve the fit with \code{\link{nm_job_result}}.
#'
#' @inheritParams nm_est
#' @param label Optional short label for the job queue UI.
#' @param job_root Job storage directory; defaults to \code{\link{nm_job_root}}.
#' @param server Remote server id from \code{\link{nm_remote_server_list}}; \code{NULL} runs locally.
#' @param data_ref Optional list with \code{dataset_id} and \code{md5} for cluster-hosted data.
#' @return A list with \code{id}, \code{path}, and \code{process} (a \code{callr} background process).
#' @examples
#' \dontrun{
#' if (requireNamespace("callr", quietly = TRUE)) {
#'   sim <- nm_synthetic_theo(n_sub = 2L)
#'   job <- nm_job_submit(sim$model, sim$data, method = "FO",
#'                        control = list(maxit = 3L, compute_inference = FALSE))
#'   nm_job_status(job$id)
#' }
#' }
#' @export
nm_job_submit <- function(model,
                          data,
                          method = "FO",
                          label = NULL,
                          job_root = nm_job_root(),
                          server = NULL,
                          data_ref = NULL,
                          workspace_project = NULL,
                          workspace_version_id = NULL,
                          est_run_id = NULL,
                          workspace_root = NULL,
                          ...) {
  if (!is.null(server) && nzchar(as.character(server))) {
    return(.nm_job_submit_remote(
      model = model,
      data = data,
      method = method,
      label = label,
      server = server,
      data_ref = data_ref,
      job_root = job_root,
      workspace_project = workspace_project,
      workspace_version_id = workspace_version_id,
      est_run_id = est_run_id,
      workspace_root = workspace_root,
      ...
    ))
  }
  if (!requireNamespace("callr", quietly = TRUE)) {
    stop("Package 'callr' is required for nm_job_submit(). Install with install.packages('callr').")
  }
  nm_validate_model(model, data = data, stop_on_error = TRUE)
  job_id <- .nm_job_make_id()
  job_path <- .nm_job_path(job_id, job_root)
  dir.create(job_path, recursive = TRUE, showWarnings = FALSE)

  est_args <- c(
    list(model = model, data = data, method = method),
    list(...)
  )
  saveRDS(est_args, file.path(job_path, "args.rds"))
  saveRDS(.nm_job_dev_env(), file.path(job_path, "env.rds"))

  meta <- list(
    id = job_id,
    label = if (is.null(label)) paste(method, "fit") else label,
    job_type = "est",
    status = "queued",
    method = method,
    created = as.character(Sys.time()),
    started = "",
    finished = "",
    pid = NA_integer_,
    objective = NA_real_,
    error = ""
  )
  saveRDS(meta, file.path(job_path, "meta.rds"))

  message(
    "[LibeRation] Estimation job submitted: id=", job_id,
    "  method=", method,
    "  label=", meta$label
  )

  log_path <- file.path(job_path, "worker.log")
  stderr_path <- file.path(job_path, "worker.stderr")
  stdout_path <- file.path(job_path, "worker.stdout")
  proc <- callr::r_bg(
    func = function(job_path) {
      dev_env <- readRDS(file.path(job_path, "env.rds"))
      nm_root <- if (is.null(dev_env$nm_root)) "" else as.character(dev_env$nm_root)
      ad_root <- if (is.null(dev_env$ad_root)) "" else as.character(dev_env$ad_root)
      use_dev <- identical(dev_env$mode, "dev") && nzchar(nm_root)
      if (use_dev) {
        if (!requireNamespace("pkgload", quietly = TRUE)) {
          stop("Package 'pkgload' is required to run jobs from a dev load_all session.", call. = FALSE)
        }
        if (nzchar(ad_root)) {
          pkgload::load_all(ad_root, quiet = TRUE, compile = FALSE)
        }
        pkgload::load_all(nm_root, quiet = TRUE, recompile = FALSE)
      } else {
        if (!requireNamespace("LibeRation", quietly = TRUE)) {
          stop("Package LibeRation is not installed.", call. = FALSE)
        }
        tryCatch(
          suppressPackageStartupMessages(library(LibeRation)),
          error = function(e) {
            stop("Failed to load LibeRation: ", conditionMessage(e), call. = FALSE)
          }
        )
      }
      LibeRation:::.nm_job_worker_impl(job_path)
    },
    args = list(job_path = job_path),
    libpath = .libPaths(),
    repos = getOption("repos"),
    stdout = stdout_path,
    stderr = stderr_path,
    supervise = TRUE
  )

  Sys.sleep(0.2)
  pid <- tryCatch(proc$get_pid(), error = function(e) NA_integer_)
  if (is.na(pid)) {
    Sys.sleep(0.2)
    pid <- tryCatch(proc$get_pid(), error = function(e) NA_integer_)
  }
  meta$status <- "running"
  meta$pid <- pid
  meta$started <- as.character(Sys.time())
  saveRDS(meta, file.path(job_path, "meta.rds"))
  tryCatch(saveRDS(proc, file.path(job_path, ".process.rds")), error = function(e) NULL)

  structure(
    list(id = job_id, path = job_path, process = proc),
    class = "nm_job_handle"
  )
}

#' Submit an asynchronous simulation job
#'
#' @inheritParams nm_simulate
#' @param project Workspace project name (simulation saved under this project).
#' @param version_id Model version id.
#' @param sim_id Pre-allocated simulation id from \code{\link{nm_workspace_new_sim_id}}.
#' @param workspace_root Workspace root passed to the worker for saving results.
#' @examples
#' \dontrun{
#' if (requireNamespace("callr", quietly = TRUE)) {
#'   ws <- tempfile()
#'   nm_workspace_init(ws, create_demo_project = FALSE)
#'   nm_workspace_create_project("demo", path = ws, template = "theo")
#'   ver <- nm_workspace_list_versions("demo", root = ws)[[1L]]
#'   sim <- nm_synthetic_theo(n_sub = 2L)
#'   sim_id <- nm_workspace_new_sim_id("demo", ver, root = ws)
#'   nm_job_submit_sim(sim$model, sim$data, "demo", ver, sim_id,
#'                     n_sim = 1L, workspace_root = ws)
#' }
#' }
#' @export
nm_job_submit_sim <- function(model,
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
  if (!is.null(server) && nzchar(as.character(server))) {
    return(.nm_job_submit_remote_sim(
      model = model,
      data = data,
      project = project,
      version_id = version_id,
      sim_id = sim_id,
      n_sim = n_sim,
      seed = seed,
      n_cores = n_cores,
      pk_engine = pk_engine,
      theta = theta,
      omega = omega,
      sigma = sigma,
      label = label,
      use_fit = use_fit,
      est_run_id = est_run_id,
      design = design,
      vpc = vpc,
      sim_compute_npc = sim_compute_npc,
      sim_compute_npde = sim_compute_npde,
      diag_n_sim = diag_n_sim,
      diag_refit_eta = diag_refit_eta,
      diag_only = diag_only,
      workspace_root = workspace_root,
      job_root = job_root,
      server = server,
      fit = fit
    ))
  }
  if (!requireNamespace("callr", quietly = TRUE)) {
    stop("Package 'callr' is required for nm_job_submit_sim().")
  }
  job_id <- .nm_job_make_id()
  job_path <- .nm_job_path(job_id, job_root)
  dir.create(job_path, recursive = TRUE, showWarnings = FALSE)

  sim_args <- list(
    model = model,
    data = data,
    n_sim = as.integer(n_sim),
    seed = as.integer(seed),
    n_cores = as.integer(n_cores),
    pk_engine = pk_engine,
    theta = theta,
    omega = omega,
    sigma = sigma,
    project = project,
    version_id = version_id,
    sim_id = sim_id,
    workspace_root = workspace_root,
    label = label,
    use_fit = isTRUE(use_fit),
    est_run_id = est_run_id,
    design = design,
    vpc = isTRUE(vpc),
    sim_compute_npc = isTRUE(sim_compute_npc),
    sim_compute_npde = isTRUE(sim_compute_npde),
    diag_n_sim = as.integer(diag_n_sim %||% 50L),
    diag_refit_eta = isTRUE(diag_refit_eta),
    diag_only = isTRUE(diag_only)
  )
  saveRDS(sim_args, file.path(job_path, "args.rds"))
  saveRDS(.nm_job_dev_env(), file.path(job_path, "env.rds"))

  meta <- list(
    id = job_id,
    label = if (is.null(label) || !nzchar(label)) {
      if (isTRUE(diag_only)) {
        "NPC/NPDE diagnostic"
      } else {
        paste("sim", sim_id)
      }
    } else {
      label
    },
    job_type = "sim",
    status = "queued",
    method = "SIM",
    project = project,
    version_id = version_id,
    sim_id = sim_id,
    n_sim = as.integer(n_sim),
    created = as.character(Sys.time()),
    started = "",
    finished = "",
    pid = NA_integer_,
    objective = NA_real_,
    error = ""
  )
  saveRDS(meta, file.path(job_path, "meta.rds"))

  message(
    "[LibeRation] Simulation job submitted: id=", job_id,
    "  sim_id=", sim_id,
    "  replicates=", n_sim,
    "  label=", meta$label
  )

  log_path <- file.path(job_path, "worker.log")
  stderr_path <- file.path(job_path, "worker.stderr")
  stdout_path <- file.path(job_path, "worker.stdout")
  proc <- callr::r_bg(
    func = function(job_path) {
      dev_env <- readRDS(file.path(job_path, "env.rds"))
      nm_root <- if (is.null(dev_env$nm_root)) "" else as.character(dev_env$nm_root)
      ad_root <- if (is.null(dev_env$ad_root)) "" else as.character(dev_env$ad_root)
      use_dev <- identical(dev_env$mode, "dev") && nzchar(nm_root)
      if (use_dev) {
        if (!requireNamespace("pkgload", quietly = TRUE)) {
          stop("Package 'pkgload' is required to run jobs from a dev load_all session.", call. = FALSE)
        }
        if (nzchar(ad_root)) {
          pkgload::load_all(ad_root, quiet = TRUE, compile = FALSE)
        }
        pkgload::load_all(nm_root, quiet = TRUE, recompile = FALSE)
      } else {
        if (!requireNamespace("LibeRation", quietly = TRUE)) {
          stop("Package LibeRation is not installed.", call. = FALSE)
        }
        tryCatch(
          suppressPackageStartupMessages(library(LibeRation)),
          error = function(e) {
            stop("Failed to load LibeRation: ", conditionMessage(e), call. = FALSE)
          }
        )
      }
      LibeRation:::.nm_job_worker_impl(job_path)
    },
    args = list(job_path = job_path),
    libpath = .libPaths(),
    repos = getOption("repos"),
    stdout = stdout_path,
    stderr = stderr_path,
    supervise = TRUE
  )

  Sys.sleep(0.2)
  pid <- tryCatch(proc$get_pid(), error = function(e) NA_integer_)
  if (is.na(pid)) {
    Sys.sleep(0.2)
    pid <- tryCatch(proc$get_pid(), error = function(e) NA_integer_)
  }
  meta$status <- "running"
  meta$pid <- pid
  meta$started <- as.character(Sys.time())
  saveRDS(meta, file.path(job_path, "meta.rds"))
  tryCatch(saveRDS(proc, file.path(job_path, ".process.rds")), error = function(e) NULL)

  structure(
    list(id = job_id, path = job_path, process = proc),
    class = "nm_job_handle"
  )
}

#' @rdname nm_job_submit
#' @method print nm_job_handle
#' @param x A job handle object.
#' @param ... Unused.
#' @examples
#' \dontrun{
#' if (requireNamespace("callr", quietly = TRUE)) {
#'   sim <- nm_synthetic_theo(n_sub = 2L)
#'   job <- nm_job_submit(sim$model, sim$data, method = "FO",
#'                        control = list(maxit = 3L, compute_inference = FALSE))
#'   print(job)
#' }
#' }
#' @export
print.nm_job_handle <- function(x, ...) {
  cat("LibeRation job:", x$id, "\n")
  st <- nm_job_status(x$id)
  if (!is.null(st)) {
    cat("  status:", st$status, "\n")
  }
  invisible(x)
}

#' Read job status
#'
#' @param job_id Job identifier returned by \code{\link{nm_job_submit}}.
#' @param job_root Job storage directory.
#' @return A list with status fields, or \code{NULL} if unknown.
#' @keywords internal
.nm_job_start_grace_sec <- function() {
  5
}

#' @keywords internal
.nm_job_pid_alive <- function(pid) {
  if (is.null(pid) || length(pid) != 1L || is.na(pid) || pid <= 0L) {
    return(NA)
  }
  if (.Platform$OS.type == "windows") {
    out <- suppressWarnings(tryCatch(
      system2("tasklist", c("/FI", paste0("PID eq ", pid)), stdout = TRUE, stderr = FALSE),
      error = function(e) character()
    ))
    if (length(out) == 0L) {
      out <- tryCatch(
        shell(
          sprintf('tasklist /FI "PID eq %d"', as.integer(pid)),
          intern = TRUE,
          mustWork = FALSE
        ),
        error = function(e) character()
      )
    }
    if (length(out) == 0L) {
      return(NA)
    }
    if (any(grepl(as.character(pid), out, fixed = TRUE))) {
      return(TRUE)
    }
    return(FALSE)
  }
  ps_ok <- suppressWarnings(tryCatch({
    system2("ps", c("-p", as.character(pid)), stdout = FALSE, stderr = FALSE) == 0L
  }, error = function(e) NA))
  if (!is.na(ps_ok)) {
    return(ps_ok)
  }
  file.exists(file.path("/proc", as.character(pid)))
}

#' @keywords internal
.nm_job_in_start_grace <- function(meta) {
  started <- meta$started %||% meta$created %||% ""
  if (!nzchar(started)) {
    return(FALSE)
  }
  t0 <- suppressWarnings(as.POSIXct(started, tz = ""))
  if (is.na(t0)) {
    return(FALSE)
  }
  difftime(Sys.time(), t0, units = "secs") < .nm_job_start_grace_sec()
}

#' @keywords internal
.nm_job_process_alive <- function(job_path) {
  proc_path <- file.path(job_path, ".process.rds")
  if (!file.exists(proc_path)) {
    return(NA)
  }
  if (!requireNamespace("callr", quietly = TRUE)) {
    return(NA)
  }
  proc <- tryCatch(readRDS(proc_path), error = function(e) NULL)
  if (is.null(proc) || !is.function(proc$is_alive)) {
    return(NA)
  }
  alive <- suppressWarnings(tryCatch(proc$is_alive(), error = function(e) NA))
  if (identical(alive, TRUE)) {
    return(TRUE)
  }
  if (identical(alive, FALSE)) {
    return(FALSE)
  }
  NA
}

#' @keywords internal
.nm_job_estimation_phase <- function(job_path) {
  log_path <- file.path(job_path, "worker.log")
  if (!file.exists(log_path)) {
    return(FALSE)
  }
  lines <- .nm_file_read_lines(log_path)
  txt <- paste(lines, collapse = "\n")
  grepl("Running estimation:|Running simulation:", txt) &&
    !grepl("Job completed\\.|Simulation completed:|Job failed:|Worker error:", txt)
}

#' @keywords internal
.nm_job_log_loading_phase <- function(job_path) {
  log_path <- file.path(job_path, "worker.log")
  if (!file.exists(log_path)) {
    return(TRUE)
  }
  lines <- .nm_file_read_lines(log_path)
  if (length(lines) == 0L) {
    return(TRUE)
  }
  txt <- paste(lines, collapse = "\n")
  if (grepl("Job started:|Running estimation:|Running simulation:", txt)) {
    return(FALSE)
  }
  grepl("Loading dev |Loading installed ", txt)
}

#' @keywords internal
.nm_job_elapsed_sec <- function(meta) {
  started <- meta$started %||% meta$created %||% ""
  if (!nzchar(started)) {
    return(NA_real_)
  }
  t0 <- suppressWarnings(as.POSIXct(started, tz = ""))
  if (is.na(t0)) {
    return(NA_real_)
  }
  as.numeric(difftime(Sys.time(), t0, units = "secs"))
}

#' @keywords internal
.nm_job_active_pid <- function(meta, job_path) {
  pid_path <- file.path(job_path, "worker.pid")
  if (file.exists(pid_path)) {
    pid <- suppressWarnings(as.integer(.nm_file_read_lines(pid_path)[1L]))
    if (!is.na(pid) && pid > 0L) {
      return(pid)
    }
  }
  meta$pid
}

#' @keywords internal
.nm_job_worker_likely_running <- function(meta, job_path) {
  log_hint <- .nm_job_log_status_hint(job_path)
  if (identical(log_hint, "success") || identical(log_hint, "error")) {
    return(FALSE)
  }
  proc_alive <- .nm_job_process_alive(job_path)
  if (identical(proc_alive, TRUE)) {
    return(TRUE)
  }
  if (identical(proc_alive, FALSE)) {
    return(FALSE)
  }
  pid <- .nm_job_active_pid(meta, job_path)
  alive <- .nm_job_pid_alive(pid)
  if (identical(alive, TRUE)) {
    return(TRUE)
  }
  if (identical(alive, FALSE)) {
    return(FALSE)
  }
  if (.nm_job_log_loading_phase(job_path)) {
    elapsed <- .nm_job_elapsed_sec(meta)
    if (!is.finite(elapsed) || elapsed < 900) {
      return(TRUE)
    }
  }
  elapsed <- .nm_job_elapsed_sec(meta)
  if (is.finite(elapsed) && elapsed < .nm_job_start_grace_sec()) {
    return(TRUE)
  }
  hb_path <- file.path(job_path, "worker.heartbeat")
  if (file.exists(hb_path)) {
    hb <- suppressWarnings(as.POSIXct(.nm_file_read_lines(hb_path)[1L], tz = ""))
    if (!is.na(hb) && difftime(Sys.time(), hb, units = "secs") < 120) {
      return(TRUE)
    }
  }
  log_path <- file.path(job_path, "worker.log")
  if (file.exists(log_path)) {
    if (isTRUE(file.info(log_path)$mtime > Sys.time() - 120)) {
      return(TRUE)
    }
  }
  FALSE
}

#' @keywords internal
.nm_job_log_status_hint <- function(job_path) {
  log_path <- file.path(job_path, "worker.log")
  if (!file.exists(log_path)) {
    return("")
  }
  lines <- .nm_file_read_lines(log_path)
  if (length(lines) == 0L) {
    return("")
  }
  tail <- lines[max(1L, length(lines) - 29L):length(lines)]
  txt <- paste(tail, collapse = "\n")
  if (grepl("Job completed\\.|Simulation completed:", txt)) {
    return("success")
  }
  if (grepl("Job failed:|Worker error:|Error in |Execution halted", txt, ignore.case = TRUE)) {
    return("error")
  }
  ""
}

#' @keywords internal
.nm_job_worker_log_hint <- function(job_path) {
  log_path <- file.path(job_path, "worker.log")
  if (!file.exists(log_path)) {
    return("")
  }
  lines <- .nm_file_read_lines(log_path)
  if (length(lines) == 0L) {
    return("")
  }
  tail <- lines[max(1L, length(lines) - 4L):length(lines)]
  paste(tail, collapse = "\n")
}

#' @keywords internal
.nm_job_reconcile_status <- function(meta, job_path) {
  if (is.null(meta)) {
    return(NULL)
  }
  result_path <- file.path(job_path, "result.rds")
  error_path <- file.path(job_path, "error.txt")
  status <- meta$status %||% "queued"

  if (file.exists(result_path)) {
    meta$status <- "success"
    meta$error <- ""
    if (meta$status %in% c("success", "error", "cancelled") &&
        (is.null(meta$finished) || !nzchar(as.character(meta$finished)))) {
      log_path <- file.path(job_path, "worker.log")
      if (file.exists(log_path)) {
        meta$finished <- as.character(file.info(log_path)$mtime)
      } else {
        meta$finished <- as.character(Sys.time())
      }
    }
    return(meta)
  }

  running_like <- status %in% c("running", "queued") ||
    (identical(status, "error") && .nm_job_worker_likely_running(meta, job_path))

  if (running_like) {
    if (identical(.nm_job_process_alive(job_path), TRUE)) {
      meta$status <- "running"
      meta$error <- ""
      if (file.exists(error_path)) {
        unlink(error_path)
      }
      return(meta)
    }
    if (.nm_job_in_start_grace(meta)) {
      meta$status <- "running"
      meta$error <- ""
      if (file.exists(error_path)) {
        unlink(error_path)
      }
      return(meta)
    }
    log_path <- file.path(job_path, "worker.log")
    if (file.exists(log_path)) {
      log_tail <- tail(.nm_file_read_lines(log_path), 20L)
      log_txt <- paste(log_tail, collapse = "\n")
      if (grepl("Job failed:|Error in|Execution halted|Worker error:", log_txt, ignore.case = TRUE)) {
        meta$status <- "error"
        if (file.exists(error_path)) {
          meta$error <- paste(.nm_file_read_lines(error_path), collapse = "\n")
        } else if (grepl("Job failed:", log_txt, fixed = TRUE)) {
          meta$error <- sub(".*Job failed:\\s*", "", log_tail[length(log_tail)])
        } else {
          meta$error <- log_tail[length(log_tail)]
        }
        if (meta$status %in% c("success", "error", "cancelled") &&
            (is.null(meta$finished) || !nzchar(as.character(meta$finished)))) {
          meta$finished <- as.character(file.info(log_path)$mtime)
        }
        return(meta)
      }
      if (grepl("Job completed\\.|Simulation completed:", log_txt)) {
        if (file.exists(result_path)) {
          meta$status <- "success"
          meta$error <- ""
          fit <- tryCatch(readRDS(result_path), error = function(e) NULL)
          if (!is.null(fit) && !is.null(fit$objective)) {
            meta$objective <- fit$objective
          }
          meta$finished <- as.character(file.info(log_path)$mtime)
          return(meta)
        }
        meta$status <- "error"
        meta$error <- "Worker log reports completion but result.rds is missing."
        meta$finished <- as.character(file.info(log_path)$mtime)
        return(meta)
      }
    }
    if (.nm_job_worker_likely_running(meta, job_path)) {
      meta$status <- "running"
      meta$error <- ""
      if (file.exists(error_path)) {
        unlink(error_path)
      }
      return(meta)
    }
    meta$status <- "error"
    if (is.null(meta$error) || !nzchar(meta$error)) {
      hint <- .nm_job_worker_log_hint(job_path)
      meta$error <- if (nzchar(hint)) {
        paste0("Worker exited without result.\n", hint)
      } else {
        "Worker exited without result."
      }
    }
    if (!file.exists(error_path)) {
      writeLines(meta$error, error_path)
    }
    return(meta)
  }

  if (file.exists(error_path) && !file.exists(result_path)) {
    meta$status <- "error"
    err_txt <- paste(.nm_file_read_lines(error_path), collapse = "\n")
    if (nzchar(trimws(err_txt))) {
      meta$error <- err_txt
    }
  }
  if (identical(meta$status, "error") && file.exists(result_path)) {
    meta$status <- "success"
    meta$error <- ""
    fit <- tryCatch(readRDS(result_path), error = function(e) NULL)
    if (!is.null(fit) && !is.null(fit$objective)) {
      meta$objective <- fit$objective
    }
  }
  if (meta$status %in% c("success", "error", "cancelled") &&
      (is.null(meta$finished) || !nzchar(as.character(meta$finished)))) {
    log_path <- file.path(job_path, "worker.log")
    if (file.exists(log_path)) {
      meta$finished <- as.character(file.info(log_path)$mtime)
    } else {
      meta$finished <- as.character(Sys.time())
    }
  }
  meta
}

#' Read status of a background estimation or simulation job
#'
#' @param job_id Job identifier returned by \code{\link{nm_job_submit}}.
#' @param job_root Job storage directory; defaults to \code{\link{nm_job_root}}.
#' @return Job metadata list, or \code{NULL} if the job does not exist.
#' @examples
#' nm_job_list()
#' @export
nm_job_status <- function(job_id, job_root = nm_job_root()) {
  meta <- .nm_job_read_meta(job_id, job_root)
  if (is.null(meta)) {
    return(NULL)
  }
  if (isTRUE(meta$remote)) {
    return(.nm_job_status_remote(job_id, job_root))
  }
  job_path <- .nm_job_path(job_id, job_root)
  old_status <- meta$status
  meta <- .nm_job_reconcile_status(meta, job_path)
  if (!identical(old_status, meta$status)) {
    .nm_save_rds_safe(meta, .nm_job_meta_path(job_id, job_root))
  }
  meta
}

#' Watch signature for job directory changes (for Shiny auto-refresh)
#'
#' @param job_root Job storage directory.
#' @return Character signature string.
#' @examples
#' nm_job_watch_signature()
#' @export
nm_job_watch_signature <- function(job_root = nm_job_root(), light = FALSE,
                                   local_only = FALSE) {
  .nm_job_watch_signature(job_root, light = light, local_only = local_only)
}

#' @keywords internal
.nm_job_watch_signature <- function(job_root = nm_job_root(), light = FALSE,
                                   local_only = FALSE) {
  if (!dir.exists(job_root)) {
    return("")
  }
  ids <- list.dirs(job_root, full.names = FALSE, recursive = FALSE)
  ids <- ids[nzchar(ids)]
  if (length(ids) == 0L) {
    return("")
  }
  if (isTRUE(local_only)) {
    ids <- ids[vapply(ids, function(id) {
      meta <- .nm_job_read_meta(id, job_root)
      !is.null(meta) && !isTRUE(meta$remote)
    }, logical(1L))]
    if (length(ids) == 0L) {
      return("")
    }
  }
  parts <- vapply(ids, function(id) {
    job_path <- file.path(job_root, id)
    meta <- file.path(job_path, "meta.rds")
    log <- file.path(job_path, "worker.log")
    prog <- file.path(job_path, "worker.progress")
    result <- file.path(job_path, "result.rds")
    err <- file.path(job_path, "error.txt")
    if (isTRUE(light)) {
      return(paste(
        id,
        if (file.exists(meta)) file.info(meta)$mtime else 0,
        if (file.exists(prog)) file.info(prog)$mtime else 0,
        if (file.exists(prog)) file.info(prog)$size else 0,
        if (file.exists(result)) file.info(result)$mtime else 0,
        if (file.exists(err)) file.info(err)$mtime else 0
      ))
    }
    paste(
      id,
      if (file.exists(meta)) file.info(meta)$mtime else 0,
      if (file.exists(log)) file.info(log)$mtime else 0,
      if (file.exists(log)) file.info(log)$size else 0,
      if (file.exists(prog)) file.info(prog)$mtime else 0,
      if (file.exists(prog)) file.info(prog)$size else 0,
      if (file.exists(result)) file.info(result)$mtime else 0,
      if (file.exists(err)) file.info(err)$mtime else 0
    )
  }, character(1L))
  paste(parts, collapse = "|")
}

#' List jobs in the job root
#'
#' @param job_root Job storage directory.
#' @return A data.frame sorted by creation time (newest first).
#' @examples
#' nm_job_list()
#' @export
nm_job_list <- function(job_root = nm_job_root(), local_only = FALSE,
                        remote_only = FALSE, sync_active = TRUE) {
  if (isTRUE(local_only) && isTRUE(remote_only)) {
    stop("Specify at most one of local_only and remote_only.", call. = FALSE)
  }
  if (!dir.exists(job_root)) {
    return(data.frame(
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
    ))
  }
  ids <- list.dirs(job_root, full.names = FALSE, recursive = FALSE)
  ids <- ids[nzchar(ids)]
  ids <- ids[file.exists(file.path(job_root, ids, "meta.rds"))]
  if (length(ids) == 0L) {
    return(data.frame(
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
    ))
  }
  rows <- lapply(ids, function(id) {
    meta <- .nm_job_read_meta(id, job_root)
    if (is.null(meta)) {
      return(NULL)
    }
    if (isTRUE(local_only) && isTRUE(meta$remote)) {
      return(NULL)
    }
    if (isTRUE(remote_only) && !isTRUE(meta$remote)) {
      return(NULL)
    }
    terminal <- meta$status %in% c("success", "error", "cancelled")
    job_path <- .nm_job_path(id, job_root)
    result_path <- file.path(job_path, "result.rds")
    refresh_terminal <- function() {
      if (isTRUE(meta$remote)) {
        return(nm_job_status(id, job_root))
      }
      if (file.exists(result_path)) {
        return(nm_job_status(id, job_root))
      }
      if (identical(meta$status, "error") &&
          .nm_job_worker_likely_running(meta, job_path)) {
        return(nm_job_status(id, job_root))
      }
      meta
    }
    if (isTRUE(meta$remote)) {
      st <- if (terminal || !isTRUE(sync_active)) {
        if (terminal) refresh_terminal() else meta
      } else {
        nm_job_status(id, job_root)
      }
    } else if (terminal || !isTRUE(sync_active)) {
      st <- if (terminal) refresh_terminal() else meta
    } else {
      st <- nm_job_status(id, job_root)
    }
    if (is.null(st)) {
      return(NULL)
    }
    data.frame(
      id = st$id,
      label = if (is.null(st$label)) "" else st$label,
      job_type = if (is.null(st$job_type)) "est" else st$job_type,
      status = st$status,
      method = if (is.null(st$method)) "" else st$method,
      sim_id = .nm_job_clean_text(st$sim_id),
      n_sim = as.integer(if (is.null(st$n_sim)) NA_integer_ else st$n_sim),
      created = if (is.null(st$created)) "" else st$created,
      started = if (is.null(st$started)) "" else as.character(st$started),
      finished = if (is.null(st$finished)) "" else as.character(st$finished),
      objective = as.numeric(if (is.null(st$objective)) NA_real_ else st$objective),
      error = as.character(st$error %||% ""),
      server = if (isTRUE(st$remote)) {
        as.character(st$server_name %||% st$server_id %||% "remote")
      } else {
        "local"
      },
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1L))]
  if (length(rows) == 0L) {
    return(data.frame(
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
    ))
  }
  out <- do.call(rbind, rows)
  out <- out[order(out$created, decreasing = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Retrieve a completed job result
#'
#' @param job_id Job identifier.
#' @param job_root Job storage directory.
#' @return An \code{nm_fit} object.
#' @examples
#' # After a successful nm_job_submit(), use nm_job_result(job$id)
#' @export
nm_job_result <- function(job_id, job_root = nm_job_root()) {
  st <- nm_job_status(job_id, job_root)
  if (is.null(st)) {
    stop("Unknown job: ", job_id)
  }
  if (!identical(st$status, "success")) {
    stop("Job ", job_id, " is not successful (status = ", st$status, ").")
  }
  if (isTRUE(st$remote)) {
    return(.nm_job_result_remote(job_id, job_root))
  }
  path <- file.path(.nm_job_path(job_id, job_root), "result.rds")
  if (!file.exists(path)) {
    stop("Result file missing for job ", job_id)
  }
  readRDS(path)
}

#' Cancel a running job
#'
#' @param job_id Job identifier.
#' @param job_root Job storage directory.
#' @return Updated status list.
#' @examples
#' # After nm_job_submit(), use nm_job_cancel(job$id)
#' @export
nm_job_cancel <- function(job_id, job_root = nm_job_root()) {
  meta <- .nm_job_read_meta(job_id, job_root)
  if (is.null(meta)) {
    stop("Unknown job: ", job_id)
  }
  if (isTRUE(meta$remote)) {
    return(.nm_job_cancel_remote(job_id, job_root))
  }
  if (meta$status %in% c("success", "error", "cancelled")) {
    return(meta)
  }
  if (!is.na(meta$pid) && meta$pid > 0L) {
    tryCatch(
      tools::pskill(meta$pid),
      error = function(e) NULL
    )
  }
  meta$status <- "cancelled"
  meta$finished <- as.character(Sys.time())
  meta$error <- "Cancelled by user."
  .nm_save_rds_safe(meta, .nm_job_meta_path(job_id, job_root))
  meta
}

#' Read worker log for a job
#'
#' @param job_id Job identifier.
#' @param tail Number of lines to return from the end of the log.
#' @examples
#' # After nm_job_submit(), use nm_job_log(job$id)
#' @export
nm_job_log <- function(job_id, tail = 40L, job_root = nm_job_root()) {
  meta <- .nm_job_read_meta(job_id, job_root)
  if (!is.null(meta) && isTRUE(meta$remote)) {
    return(.nm_job_log_remote(job_id, tail = tail, job_root = job_root))
  }
  path <- file.path(.nm_job_path(job_id, job_root), "worker.log")
  if (!file.exists(path)) {
    return("")
  }
  lines <- .nm_file_read_lines(path)
  if (length(lines) <= tail) {
    return(paste(lines, collapse = "\n"))
  }
  paste(tail(lines, tail), collapse = "\n")
}

#' Remove finished jobs from disk
#'
#' @param job_root Job storage directory.
#' @param status Status values to remove (\code{"success"}, \code{"error"}, \code{"cancelled"}).
#' @return Number of jobs removed.
#' @examples
#' nm_job_cleanup()
#' @export
nm_job_cleanup <- function(job_root = nm_job_root(),
                           status = c("success", "error", "cancelled")) {
  jobs <- nm_job_list(job_root)
  if (nrow(jobs) == 0L) {
    return(0L)
  }
  rm_ids <- jobs$id[jobs$status %in% status]
  for (id in rm_ids) {
    unlink(.nm_job_path(id, job_root), recursive = TRUE)
  }
  length(rm_ids)
}
