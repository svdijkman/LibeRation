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
  readRDS(path)
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
      cat(..., file = log_path, append = TRUE)
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
.nm_job_worker_impl <- function(job_path) {
  meta_path <- file.path(job_path, "meta.rds")
  log_path <- file.path(job_path, "worker.log")
  meta <- readRDS(meta_path)
  meta$status <- "running"
  meta$started <- as.character(Sys.time())
  saveRDS(meta, meta_path)
  job_type <- if (is.null(meta$job_type)) "est" else meta$job_type
  cat("LibeRation worker\nJob started:", meta$id, "(", job_type, ")\n", file = log_path)

  dev_env <- readRDS(file.path(job_path, "env.rds"))
  nm_root <- if (is.null(dev_env$nm_root)) "" else as.character(dev_env$nm_root)
  ad_root <- if (is.null(dev_env$ad_root)) "" else as.character(dev_env$ad_root)
  use_dev <- identical(dev_env$mode, "dev") && nzchar(nm_root)

  if (use_dev) {
    cat("Loading dev sources:", nm_root, "\n", file = log_path, append = TRUE)
    if (!requireNamespace("pkgload", quietly = TRUE)) {
      stop("Package 'pkgload' is required to run jobs from a dev load_all session.", call. = FALSE)
    }
    if (nzchar(ad_root)) {
      pkgload::load_all(ad_root, quiet = TRUE, compile = FALSE)
    }
    pkgload::load_all(nm_root, quiet = TRUE, recompile = FALSE)
  } else {
    cat("Loading installed package LibeRation\n", file = log_path, append = TRUE)
    pkg <- "LibeRation"
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package LibeRation is not installed.", call. = FALSE)
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }

  args <- readRDS(file.path(job_path, "args.rds"))
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
    meta <- readRDS(meta_path)
    sim_ok <- FALSE
    sim_data <- NULL
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
        writeLines(err, file.path(job_path, "error.txt"))
        cat("Simulation failed:", err, "\n", file = log_path, append = TRUE)
      } else {
        sim_data <- LibeRation:::.nm_sim_pack_output(sim_out)
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
        sim_ok <- TRUE
        cat(
          "Simulation completed:", args$sim_id,
          " replicates=", args$n_sim, "\n",
          file = log_path, append = TRUE
        )
      }
    } else if (run_diag) {
      sim_ok <- TRUE
      cat("Diagnostic-only job (NPC/NPDE); skipping workspace simulation save.\n",
          file = log_path, append = TRUE)
    } else {
      meta$status <- "error"
      meta$error <- "Nothing to run: enable VPC, replications, NPC, or NPDE."
      writeLines(meta$error, file.path(job_path, "error.txt"))
      cat("Job failed:", meta$error, "\n", file = log_path, append = TRUE)
    }

    if (sim_ok && run_diag) {
      nm_workspace_load_run_fit <- getExportedValue("LibeRation", "nm_workspace_load_run_fit")
      nm_workspace_save_run <- getExportedValue("LibeRation", "nm_workspace_save_run")
      nm_workspace_list_runs <- getExportedValue("LibeRation", "nm_workspace_list_runs")
      nm_add_npc_npde <- getExportedValue("LibeRation", "nm_add_npc_npde")
      est_run_id <- args$est_run_id
      if (is.null(est_run_id) || !nzchar(est_run_id)) {
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
        cat(
          paste(what, collapse = "/"), "on fit:", est_run_id,
          " — ", diag_n_sim, "simulations\n",
          file = log_path, append = TRUE
        )
        fit_diag <- tryCatch(
          nm_workspace_load_run_fit(
            args$project, args$version_id, est_run_id, root = args$workspace_root
          ),
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
              cat("NPC/NPDE failed:", conditionMessage(e), "\n", file = log_path, append = TRUE)
              err <<- conditionMessage(e)
              NULL
            }
          )
          if (!is.null(fit_diag)) {
            nm_workspace_save_run(
              args$project,
              args$version_id,
              est_run_id,
              fit_diag,
              root = args$workspace_root,
              job_id = meta$id
            )
            n_ok <- fit_diag$npc_npde$n_ok %||% NA_integer_
            cat(
              paste(what, collapse = "/"), "completed:", n_ok, "obs rows /",
              diag_n_sim, "simulations\n",
              file = log_path, append = TRUE
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
            writeLines(meta$error, file.path(job_path, "error.txt"))
          }
        } else if (!is.null(fit_diag) && identical(fit_diag$method, "BAYES")) {
          cat("NPC/NPDE skipped: BAYES fit\n", file = log_path, append = TRUE)
        }
      } else {
        cat("NPC/NPDE skipped: no estimation run for fit\n", file = log_path, append = TRUE)
      }
    }

    if (identical(meta$status, "queued") || identical(meta$status, "running")) {
      if (sim_ok) {
        result <- list(sim_data = sim_data, sim_id = args$sim_id)
        if (!is.null(meta$diag_est_run)) {
          result$diag_est_run <- meta$diag_est_run
        }
        saveRDS(result, file.path(job_path, "result.rds"))
        meta$status <- "success"
        meta$sim_id <- args$sim_id
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
    fit <- tryCatch(
      do.call(nm_est_fun, est_args),
      error = function(e) {
        err <<- conditionMessage(e)
        NULL
      }
    )
    meta <- readRDS(meta_path)
    if (is.null(fit)) {
      meta$status <- "error"
      meta$error <- err
      writeLines(err, file.path(job_path, "error.txt"))
      cat("Job failed:", err, "\n", file = log_path, append = TRUE)
    } else {
      if (bootstrap_n > 0L && !identical(fit$method, "BAYES")) {
        cat(
          "Bootstrap SE:", bootstrap_n, "replicates (seed=", bootstrap_seed, ")\n",
          file = log_path, append = TRUE
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
            cat("Bootstrap failed:", conditionMessage(e), "\n", file = log_path, append = TRUE)
            fit
          }
        )
        if (!is.null(fit$bootstrap)) {
          cat(
            "Bootstrap completed:", fit$bootstrap$n_ok, "/", bootstrap_n, "\n",
            file = log_path, append = TRUE
          )
          meta$bootstrap_n <- bootstrap_n
          meta$bootstrap_ok <- fit$bootstrap$n_ok
        }
      }
      saveRDS(fit, file.path(job_path, "result.rds"))
      meta$status <- "success"
      meta$objective <- fit$objective
      meta$method <- fit$method
      meta$error <- ""
      cat("Job completed. objective =", fit$objective, "\n", file = log_path, append = TRUE)
    }
  }
  meta$finished <- as.character(Sys.time())
  saveRDS(meta, meta_path)
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
                          ...) {
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
        suppressPackageStartupMessages(library(LibeRation))
      }
      LibeRation:::.nm_job_worker_impl(job_path)
    },
    args = list(job_path = job_path),
    libpath = .libPaths(),
    repos = getOption("repos"),
    stdout = log_path,
    stderr = log_path,
    supervise = TRUE
  )

  Sys.sleep(0.2)
  pid <- tryCatch(proc$get_pid(), error = function(e) NA_integer_)
  meta$status <- "running"
  meta$pid <- pid
  meta$started <- as.character(Sys.time())
  saveRDS(meta, file.path(job_path, "meta.rds"))

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
                              job_root = nm_job_root()) {
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
        suppressPackageStartupMessages(library(LibeRation))
      }
      LibeRation:::.nm_job_worker_impl(job_path)
    },
    args = list(job_path = job_path),
    libpath = .libPaths(),
    repos = getOption("repos"),
    stdout = log_path,
    stderr = log_path,
    supervise = TRUE
  )

  Sys.sleep(0.2)
  pid <- tryCatch(proc$get_pid(), error = function(e) NA_integer_)
  meta$status <- "running"
  meta$pid <- pid
  meta$started <- as.character(Sys.time())
  saveRDS(meta, file.path(job_path, "meta.rds"))

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
.nm_job_reconcile_status <- function(meta, job_path) {
  if (is.null(meta)) {
    return(NULL)
  }
  result_path <- file.path(job_path, "result.rds")
  error_path <- file.path(job_path, "error.txt")
  status <- meta$status %||% "queued"
  if (file.exists(error_path) && !file.exists(result_path)) {
    meta$status <- "error"
    err_txt <- paste(readLines(error_path, warn = FALSE), collapse = "\n")
    if (nzchar(trimws(err_txt))) {
      meta$error <- err_txt
    }
  } else if (file.exists(result_path)) {
    meta$status <- "success"
  } else if (status %in% c("running", "queued")) {
    log_path <- file.path(job_path, "worker.log")
    if (file.exists(log_path)) {
      log_tail <- tail(readLines(log_path, warn = FALSE), 20L)
      log_txt <- paste(log_tail, collapse = "\n")
      if (grepl("Job failed:|Error in|Execution halted", log_txt, ignore.case = TRUE)) {
        meta$status <- "error"
        if (file.exists(error_path)) {
          meta$error <- paste(readLines(error_path, warn = FALSE), collapse = "\n")
        } else if (grepl("Job failed:", log_txt, fixed = TRUE)) {
          meta$error <- sub(".*Job failed:\\s*", "", log_tail[length(log_tail)])
        } else {
          meta$error <- log_tail[length(log_tail)]
        }
      }
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
  job_path <- .nm_job_path(job_id, job_root)
  old_status <- meta$status
  meta <- .nm_job_reconcile_status(meta, job_path)
  if (!identical(old_status, meta$status)) {
    saveRDS(meta, .nm_job_meta_path(job_id, job_root))
  }
  meta
}

#' @keywords internal
.nm_job_pid_alive <- function(pid) {
  if (is.na(pid) || pid <= 0L) {
    return(FALSE)
  }
  if (.Platform$OS.type == "windows") {
    out <- tryCatch(
      system2("tasklist", c("/FI", paste0("PID eq ", pid)), stdout = TRUE, stderr = FALSE),
      error = function(e) character()
    )
    return(any(grepl(as.character(pid), out, fixed = TRUE)))
  }
  ps_ok <- tryCatch({
    system2("ps", c("-p", as.character(pid)), stdout = FALSE, stderr = FALSE) == 0L
  }, error = function(e) NA)
  if (!is.na(ps_ok)) {
    return(ps_ok)
  }
  file.exists(file.path("/proc", as.character(pid)))
}

#' Watch signature for job directory changes (for Shiny auto-refresh)
#'
#' @param job_root Job storage directory.
#' @return Character signature string.
#' @examples
#' nm_job_watch_signature()
#' @export
nm_job_watch_signature <- function(job_root = nm_job_root()) {
  .nm_job_watch_signature(job_root)
}

#' @keywords internal
.nm_job_watch_signature <- function(job_root = nm_job_root()) {
  if (!dir.exists(job_root)) {
    return("")
  }
  ids <- list.dirs(job_root, full.names = FALSE, recursive = FALSE)
  ids <- ids[nzchar(ids)]
  if (length(ids) == 0L) {
    return("")
  }
  parts <- vapply(ids, function(id) {
    job_path <- file.path(job_root, id)
    meta <- file.path(job_path, "meta.rds")
    log <- file.path(job_path, "worker.log")
    result <- file.path(job_path, "result.rds")
    err <- file.path(job_path, "error.txt")
    paste(
      id,
      if (file.exists(meta)) file.info(meta)$mtime else 0,
      if (file.exists(log)) file.info(log)$mtime else 0,
      if (file.exists(log)) file.info(log)$size else 0,
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
nm_job_list <- function(job_root = nm_job_root()) {
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
      stringsAsFactors = FALSE
    ))
  }
  rows <- lapply(ids, function(id) {
    st <- nm_job_status(id, job_root)
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
  saveRDS(meta, .nm_job_meta_path(job_id, job_root))
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
  path <- file.path(.nm_job_path(job_id, job_root), "worker.log")
  if (!file.exists(path)) {
    return("")
  }
  lines <- readLines(path, warn = FALSE)
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
