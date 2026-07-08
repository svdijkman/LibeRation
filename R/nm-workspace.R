#' Workspace root for Pirana-style project management
#'
#' All file access in the GUI is confined under this directory.
#' Set with \code{options(LibeRation.workspace = "/path")} or \code{\link{nm_workspace_init}}.
#'
#' @return Normalized path.
#' @examples NULL
#' @export
nm_workspace_root <- function() {
  root <- getOption("LibeRation.workspace")
  if (is.null(root) || !nzchar(root)) {
    root <- file.path(path.expand("~"), "LibeRation_workspace")
  }
  normalizePath(root, winslash = "/", mustWork = FALSE)
}

#' @keywords internal
.nm_ws_marker <- function(root = nm_workspace_root()) {
  file.path(root, ".rcppnm-workspace")
}

#' @keywords internal
.nm_ws_valid_name <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
}

#' @keywords internal
.nm_ws_within <- function(path, root = nm_workspace_root()) {
  if (!.nm_ws_valid_name(path)) {
    return(FALSE)
  }
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  root <- sub("/$", "", root)
  identical(path, root) || startsWith(path, paste0(root, "/"))
}

#' @keywords internal
.nm_ws_resolve <- function(..., root = nm_workspace_root(), must_exist = FALSE) {
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  if (!dir.exists(root)) {
    .nm_stop("Workspace root does not exist: ", root)
  }
  parts <- list(...)
  for (p in parts) {
    if (!.nm_ws_valid_name(p)) {
      .nm_stop("Invalid workspace path segment.")
    }
  }
  path <- file.path(root, ...)
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (!.nm_ws_within(path, root)) {
    .nm_stop("Access outside workspace is not allowed: ", path)
  }
  if (isTRUE(must_exist) && !file.exists(path)) {
    .nm_stop("Path not found in workspace: ", path)
  }
  path
}

#' Initialize a workspace directory
#'
#' @param path Directory to use as workspace root.
#' @param create_demo_project Create a THEO demo project when \code{TRUE}.
#' @return Invisibly, the workspace root path.
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' }
#' @export
nm_workspace_init <- function(path = nm_workspace_root(), create_demo_project = TRUE) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(path, "projects"), recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c(
      "LibeRation workspace",
      paste("created:", Sys.time()),
      "version: 1"
    ),
    .nm_ws_marker(path)
  )
  options(LibeRation.workspace = path)
  if (isTRUE(create_demo_project) && length(nm_workspace_list_projects(path)) == 0L) {
    nm_workspace_create_project("theo_demo", path = path, template = "theo")
  }
  invisible(path)
}

#' List projects in the workspace
#'
#' @param root Workspace root.
#' @return Character vector of project names.
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_list_projects(ws)
#' }
#' @export
nm_workspace_list_projects <- function(root = nm_workspace_root()) {
  proj_dir <- .nm_ws_resolve("projects", root = root)
  if (!dir.exists(proj_dir)) {
    return(character())
  }
  ids <- list.dirs(proj_dir, full.names = FALSE, recursive = FALSE)
  ids[nzchar(ids)]
}

#' @keywords internal
.nm_ws_project_dir <- function(project, root = nm_workspace_root()) {
  .nm_ws_resolve("projects", project, root = root)
}

#' @keywords internal
.nm_ws_read_json <- function(path) {
  if (!file.exists(path)) {
    return(list())
  }
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    return(list(raw = txt))
  }
  jsonlite::fromJSON(txt, simplifyVector = TRUE)
}

#' @keywords internal
.nm_ws_meta_scalar <- function(x, default = "") {
  if (is.null(x) || length(x) == 0L) {
    return(default)
  }
  if (length(x) >= 1L && is.null(x[[1L]])) {
    return(default)
  }
  out <- as.character(x[[1L]])
  if (length(out) != 1L || is.na(out)) {
    default
  } else {
    out
  }
}

#' @keywords internal
.nm_ws_write_json <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE)
  } else {
    writeLines(as.character(unlist(x)), path)
  }
  invisible(path)
}

#' Create a new project directory
#'
#' @param name Project name (directory name).
#' @param path Workspace root.
#' @param template \code{"empty"} or \code{"theo"}.
#' @param description Optional free-text description stored in \code{project.json}.
#' @return Project metadata list.
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("demo", path = ws, template = "empty")
#' }
#' @export
nm_workspace_create_project <- function(name,
                                        path = nm_workspace_root(),
                                        template = c("empty", "theo"),
                                        description = NULL) {
  template <- match.arg(template)
  name <- .nm_ws_sanitize_name(name)
  proj_dir <- .nm_ws_resolve("projects", name, root = path)
  if (dir.exists(proj_dir)) {
    .nm_stop("Project already exists: ", name)
  }
  dir.create(proj_dir, recursive = TRUE)
  dir.create(file.path(proj_dir, "data"), recursive = TRUE)
  dir.create(file.path(proj_dir, "versions"), recursive = TRUE)
  dir.create(file.path(proj_dir, "runs"), recursive = TRUE)
  dir.create(file.path(proj_dir, "reports"), recursive = TRUE)
  meta <- list(
    name = name,
    created = as.character(Sys.time()),
    template = template
  )
  if (!is.null(description)) {
    desc <- trimws(as.character(description)[1L])
    if (!is.na(desc) && nzchar(desc)) {
      meta$description <- desc
    }
  }
  .nm_ws_write_json(meta, file.path(proj_dir, "project.json"))
  if (identical(template, "theo")) {
    nm_workspace_add_theo_demo(name, root = path)
  }
  meta
}

#' Copy a dataset file into a project's \code{data/} folder
#'
#' @param project Project name.
#' @param source_path Path to an existing file on disk.
#' @param dest_name Optional basename for the copied file (defaults to source basename).
#' @param root Workspace root.
#' @return Relative path suitable for \code{$DATA} (e.g. \code{"data/theo.csv"}).
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("demo", path = ws, template = "empty")
#' sim <- nm_synthetic_theo(n_sub = 2L)
#' csv <- tempfile(fileext = ".csv")
#' write.csv(sim$data$data, csv, row.names = FALSE)
#' nm_workspace_import_dataset("demo", csv, root = ws)
#' }
#' @export
nm_workspace_import_dataset <- function(project,
                                          source_path,
                                          dest_name = NULL,
                                          root = nm_workspace_root()) {
  if (!file.exists(source_path)) {
    .nm_stop("Dataset file not found: ", source_path)
  }
  if (is.null(dest_name) || !nzchar(dest_name)) {
    dest_name <- basename(source_path)
  }
  dest_name <- basename(dest_name)
  if (!grepl("\\.(csv|dat|txt|tsv)$", dest_name, ignore.case = TRUE)) {
    dest_name <- paste0(tools::file_path_sans_ext(dest_name), ".csv")
  }
  proj_dir <- .nm_ws_project_dir(project, root)
  dest_path <- file.path(proj_dir, "data", dest_name)
  if (!file.copy(source_path, dest_path, overwrite = TRUE)) {
    .nm_stop("Could not copy dataset to project.")
  }
  file.path("data", dest_name)
}

#' @keywords internal
.nm_ws_sanitize_name <- function(name) {
  if (!.nm_ws_valid_name(name)) {
    .nm_stop("Invalid project or model name.")
  }
  name <- gsub("[^A-Za-z0-9._-]", "_", name)
  name <- gsub("_+", "_", name)
  if (!nzchar(name)) {
    .nm_stop("Invalid project or model name.")
  }
  name
}

#' @keywords internal
.nm_model_to_ctl <- function(model,
                             data_file = "data.csv",
                             prob = "LibeRation model",
                             method = NULL) {
  parts <- list(
    problem = prob,
    advan = model$ADVAN,
    trans = model$TRANS,
    use_ode = isTRUE(model$USE_ODE),
    subroutine = "",
    data_file = data_file,
    input_cols = model$INPUT,
    thetas = model$THETAS,
    omegas = model$OMEGAS,
    sigmas = model$SIGMAS,
    pk = model$PRED,
    error = model$ERROR
  )
  if (isTRUE(model$USE_ODE)) {
    parts$subroutine <- if (model$ADVAN == 13L) "ADVAN13 TRANS1" else "ADVAN6 TRANS4"
  }
  ctl <- nm_ctl_compose(parts)
  if (!is.null(method) && nzchar(method)) {
    ctl <- paste0(ctl, "\n", sprintf("$EST METHOD=%s", method))
  }
  ctl
}

#' List dataset files in a project
#'
#' @param project Project name.
#' @param root Workspace root.
#' @return Data frame with \code{file} (basename) and \code{path} (relative to project).
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("theo_demo", path = ws, template = "theo")
#' nm_workspace_list_datasets("theo_demo", root = ws)
#' }
#' @export
nm_workspace_list_datasets <- function(project, root = nm_workspace_root()) {
  empty <- data.frame(file = character(), path = character(), stringsAsFactors = FALSE)
  if (!.nm_ws_valid_name(project)) {
    return(empty)
  }
  proj_dir <- .nm_ws_project_dir(project, root)
  rows <- list()
  data_dir <- file.path(proj_dir, "data")
  if (dir.exists(data_dir)) {
    files <- list.files(
      data_dir,
      pattern = "\\.(csv|dat|txt|tsv)$",
      ignore.case = TRUE,
      full.names = FALSE
    )
    for (f in files) {
      rows[[length(rows) + 1L]] <- data.frame(
        file = f,
        path = file.path("data", f),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0L) {
    return(empty)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$file), , drop = FALSE]
}

#' Delete a project from the workspace
#'
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("tmp", path = ws, template = "empty")
#' nm_workspace_delete_project("tmp", root = ws)
#' }
#' @export
nm_workspace_delete_project <- function(project, root = nm_workspace_root()) {
  if (!.nm_ws_valid_name(project)) {
    .nm_stop("Invalid project name.")
  }
  proj_dir <- .nm_ws_project_dir(project, root)
  if (!dir.exists(proj_dir)) {
    .nm_stop("Project not found: ", project)
  }
  unlink(proj_dir, recursive = TRUE)
  if (dir.exists(proj_dir)) {
    .nm_stop(
      "Could not delete project folder '", project,
      "'. Close any open files in that folder and try again."
    )
  }
  invisible(TRUE)
}

#' Delete a model version and its estimation runs
#'
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("demo", path = ws, template = "theo")
#' ver <- nm_workspace_new_version("demo", root = ws)
#' nm_workspace_delete_version("demo", ver, root = ws)
#' }
#' @export
nm_workspace_delete_version <- function(project, version_id, root = nm_workspace_root()) {
  version_id <- .nm_ws_sanitize_name(version_id)
  version_dir <- .nm_ws_version_dir(project, version_id, root = root)
  if (!dir.exists(version_dir)) {
    .nm_stop("Model version not found: ", version_id)
  }
  runs <- nm_workspace_list_runs(project, version_id, root = root)
  if (nrow(runs) > 0L) {
    for (i in seq_len(nrow(runs))) {
      nm_workspace_delete_run(project, version_id, runs$run_id[[i]], root = root)
    }
  }
  sims <- nm_workspace_list_sims(project, version_id, root = root)
  if (nrow(sims) > 0L) {
    for (i in seq_len(nrow(sims))) {
      nm_workspace_delete_sim(project, version_id, sims$sim_id[[i]], root = root)
    }
  }
  unlink(version_dir, recursive = TRUE)
  legacy <- file.path(.nm_ws_project_dir(project, root), "models", version_id)
  if (dir.exists(legacy)) {
    unlink(legacy, recursive = TRUE)
  }
  legacy_fit <- nm_workspace_fit_path(project, version_id, root = root)
  if (file.exists(legacy_fit)) {
    unlink(legacy_fit)
  }
  invisible(TRUE)
}

#' Delete an estimation run
#'
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("demo", path = ws, template = "theo")
#' ver <- nm_workspace_new_version("demo", root = ws)
#' run_id <- nm_workspace_new_run_id("demo", ver, root = ws)
#' nm_workspace_delete_run("demo", ver, run_id, root = ws)
#' }
#' @export
nm_workspace_delete_run <- function(project, version_id, run_id, root = nm_workspace_root()) {
  run_id <- .nm_ws_sanitize_name(run_id)
  run_dir <- .nm_ws_find_run_dir(project, version_id, run_id, root = root)
  if (!is.null(run_dir)) {
    unlink(run_dir, recursive = TRUE)
    return(invisible(TRUE))
  }
  legacy_id <- paste0(.nm_ws_sanitize_name(version_id), "_fit")
  if (identical(run_id, legacy_id)) {
    legacy <- nm_workspace_fit_path(project, version_id, root = root)
    if (file.exists(legacy)) {
      unlink(legacy)
      return(invisible(TRUE))
    }
  }
  .nm_stop("Estimation run not found: ", run_id)
}

#' @keywords internal
nm_workspace_add_theo_demo <- function(project, root = nm_workspace_root()) {
  sim <- nm_synthetic_theo(n_sub = 10L, seed = 1L)
  proj_dir <- .nm_ws_project_dir(project, root)
  data_path <- file.path(proj_dir, "data", "theo.csv")
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(sim$data$data, data_path)
  } else {
    write.csv(sim$data$data, data_path, row.names = FALSE)
  }
  ctl <- .nm_model_to_ctl(
    sim$model,
    data_file = "data/theo.csv",
    prob = "THEO demo",
    method = "FOCE"
  )
  nm_workspace_write_model(
    project, "theo_base", ctl, root = root, data_file = "data/theo.csv",
    label = "THEO base model"
  )
  invisible(TRUE)
}

#' @keywords internal
.nm_ws_version_ctl_path <- function(project, version_id, root = nm_workspace_root()) {
  version_id <- .nm_ws_sanitize_name(version_id)
  proj_dir <- .nm_ws_project_dir(project, root)
  cand <- file.path(proj_dir, "versions", version_id, "model.ctl")
  if (file.exists(cand)) {
    return(normalizePath(cand, winslash = "/", mustWork = FALSE))
  }
  legacy <- file.path(proj_dir, "models", version_id, "model.ctl")
  if (file.exists(legacy)) {
    return(normalizePath(legacy, winslash = "/", mustWork = FALSE))
  }
  normalizePath(cand, winslash = "/", mustWork = FALSE)
}

#' @keywords internal
.nm_ws_version_dir <- function(project, version_id, root = nm_workspace_root(), create = FALSE) {
  version_id <- .nm_ws_sanitize_name(version_id)
  proj_dir <- .nm_ws_project_dir(project, root)
  vdir <- file.path(proj_dir, "versions", version_id)
  if (create) {
    dir.create(vdir, recursive = TRUE, showWarnings = FALSE)
    return(normalizePath(vdir, winslash = "/", mustWork = FALSE))
  }
  if (dir.exists(vdir)) {
    return(normalizePath(vdir, winslash = "/", mustWork = FALSE))
  }
  legacy <- file.path(proj_dir, "models", version_id)
  if (dir.exists(legacy)) {
    return(normalizePath(legacy, winslash = "/", mustWork = FALSE))
  }
  normalizePath(vdir, winslash = "/", mustWork = FALSE)
}

#' @keywords internal
.nm_ws_version_runs_dir <- function(project,
                                    version_id,
                                    root = nm_workspace_root(),
                                    create = FALSE) {
  version_id <- .nm_ws_sanitize_name(version_id)
  vdir <- .nm_ws_version_dir(project, version_id, root = root, create = create)
  runs_dir <- file.path(vdir, "runs")
  if (create) {
    dir.create(runs_dir, recursive = TRUE, showWarnings = FALSE)
  }
  runs_dir
}

#' @keywords internal
.nm_ws_run_dir <- function(project,
                           version_id,
                           run_id,
                           root = nm_workspace_root(),
                           create = FALSE) {
  run_id <- .nm_ws_sanitize_name(run_id)
  runs_dir <- .nm_ws_version_runs_dir(project, version_id, root = root, create = create)
  run_dir <- file.path(runs_dir, run_id)
  if (create) {
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  }
  run_dir
}

#' @keywords internal
.nm_ws_legacy_runs_dir <- function(project, root = nm_workspace_root()) {
  file.path(.nm_ws_project_dir(project, root), "runs")
}

#' @keywords internal
.nm_ws_run_meta_version <- function(meta) {
  if (is.null(meta$version_id)) meta$model_version else meta$version_id
}

#' @keywords internal
.nm_ws_collect_runs_from_dir <- function(runs_dir, version_id) {
  empty <- data.frame(
    run_id = character(),
    method = character(),
    objective = numeric(),
    label = character(),
    created = character(),
    has_bootstrap = logical(),
    has_npc = logical(),
    has_npde = logical(),
    stringsAsFactors = FALSE
  )
  if (!dir.exists(runs_dir)) {
    return(empty)
  }
  ids <- list.dirs(runs_dir, full.names = FALSE, recursive = FALSE)
  ids <- ids[nzchar(ids)]
  rows <- list()
  for (rid in ids) {
    meta <- .nm_ws_read_json(file.path(runs_dir, rid, "meta.json"))
    vid <- .nm_ws_run_meta_version(meta)
    if (!identical(vid, version_id)) {
      next
    }
    fit_path <- file.path(runs_dir, rid, "fit.rds")
    obj <- meta$objective
    fit_obj <- NULL
    if (file.exists(fit_path)) {
      fit_obj <- tryCatch(readRDS(fit_path), error = function(e) NULL)
      if ((is.null(obj) || !is.finite(obj)) && !is.null(fit_obj) && is.finite(fit_obj$objective)) {
        obj <- fit_obj$objective
      }
    }
    has_bootstrap <- isTRUE(meta$has_bootstrap)
    has_npc <- isTRUE(meta$has_npc)
    has_npde <- isTRUE(meta$has_npde)
    if (!is.null(fit_obj) && (!has_bootstrap || !has_npc || !has_npde)) {
      if (!has_bootstrap) {
        has_bootstrap <- !is.null(fit_obj$bootstrap)
      }
      if (!has_npc) {
        has_npc <- .nm_fit_has_npc(fit_obj)
      }
      if (!has_npde) {
        has_npde <- .nm_fit_has_npde(fit_obj)
      }
    }
    rows[[length(rows) + 1L]] <- data.frame(
      run_id = rid,
      method = .nm_ws_meta_scalar(meta$method),
      objective = as.numeric(if (is.null(obj) || !is.finite(obj)) NA_real_ else obj),
      label = .nm_ws_meta_scalar(meta$label),
      created = .nm_ws_meta_scalar(meta$created),
      has_bootstrap = has_bootstrap,
      has_npc = has_npc,
      has_npde = has_npde,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0L) {
    return(empty)
  }
  do.call(rbind, rows)
}

#' @keywords internal
.nm_ws_find_run_dir <- function(project, version_id, run_id, root = nm_workspace_root()) {
  version_id <- .nm_ws_sanitize_name(version_id)
  run_id <- .nm_ws_sanitize_name(run_id)
  version_run <- file.path(
    .nm_ws_version_runs_dir(project, version_id, root = root),
    run_id
  )
  if (dir.exists(version_run)) {
    return(version_run)
  }
  legacy_run <- file.path(.nm_ws_legacy_runs_dir(project, root), run_id)
  if (dir.exists(legacy_run)) {
    meta <- .nm_ws_read_json(file.path(legacy_run, "meta.json"))
    vid <- .nm_ws_run_meta_version(meta)
    if (identical(vid, version_id)) {
      return(legacy_run)
    }
  }
  NULL
}

#' List model versions in a project
#'
#' Model versions are distinct control streams (different model code) within a project.
#'
#' @param project Project name.
#' @param root Workspace root.
#' @return Character vector of version ids.
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("theo_demo", path = ws, template = "theo")
#' nm_workspace_list_versions("theo_demo", root = ws)
#' }
#' @export
nm_workspace_list_versions <- function(project, root = nm_workspace_root()) {
  if (!.nm_ws_valid_name(project)) {
    return(character())
  }
  proj_dir <- .nm_ws_project_dir(project, root)
  ids <- character()
  for (sub in c("versions", "models")) {
    base <- file.path(proj_dir, sub)
    if (!dir.exists(base)) {
      next
    }
    found <- list.dirs(base, full.names = FALSE, recursive = FALSE)
    found <- found[nzchar(found)]
    found <- found[file.exists(file.path(base, found, "model.ctl"))]
    ids <- c(ids, found)
  }
  sort(unique(ids))
}

#' @rdname nm_workspace_list_versions
#' @export
nm_workspace_list_models <- function(project, root = nm_workspace_root()) {
  nm_workspace_list_versions(project, root = root)
}

#' Path to a model version control file
#'
#' @param project Project name.
#' @param version_id Model version id.
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("theo_demo", path = ws, template = "theo")
#' ver <- nm_workspace_list_versions("theo_demo", root = ws)[[1L]]
#' nm_workspace_model_path("theo_demo", ver, root = ws)
#' }
#' @export
nm_workspace_model_path <- function(project, version_id, root = nm_workspace_root()) {
  path <- .nm_ws_version_ctl_path(project, version_id, root = root)
  if (!file.exists(path)) {
    .nm_stop("Model version not found: ", version_id)
  }
  path
}

#' Read model control text for a model version
#'
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("theo_demo", path = ws, template = "theo")
#' ver <- nm_workspace_list_versions("theo_demo", root = ws)[[1L]]
#' head(strsplit(nm_workspace_read_model("theo_demo", ver, root = ws), "\n")[[1L]])
#' }
#' @export
nm_workspace_read_model <- function(project, version_id, root = nm_workspace_root()) {
  path <- nm_workspace_model_path(project, version_id, root = root)
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

#' Write model control text for a model version
#'
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("demo", path = ws, template = "empty")
#' parts <- nm_ctl_template(2L, 1L)
#' nm_workspace_write_model("demo", "mod001", nm_ctl_compose(parts), root = ws)
#' }
#' @export
nm_workspace_write_model <- function(project,
                                     version_id,
                                     ctl_text,
                                     root = nm_workspace_root(),
                                     data_file = NULL,
                                     label = NULL) {
  version_id <- .nm_ws_sanitize_name(version_id)
  version_dir <- .nm_ws_version_dir(project, version_id, root = root, create = TRUE)
  writeLines(unlist(strsplit(ctl_text, "\n", fixed = TRUE)), file.path(version_dir, "model.ctl"))
  meta <- list(
    version_id = version_id,
    updated = as.character(Sys.time()),
    data_file = data_file
  )
  if (!is.null(label)) {
    lbl <- trimws(as.character(label)[1L])
    if (!is.na(lbl) && nzchar(lbl)) {
      meta$label <- lbl
    }
  }
  meta_path <- file.path(version_dir, "meta.json")
  if (file.exists(meta_path)) {
    old <- .nm_ws_read_json(meta_path)
    if (!is.null(old$data_file) && is.null(data_file)) {
      meta$data_file <- old$data_file
    }
    if (!is.null(old$label) && is.null(label)) {
      meta$label <- old$label
    }
    if (!is.null(old$run_id) && is.null(meta$version_id)) {
      meta$version_id <- old$run_id
    }
  }
  .nm_ws_write_json(meta, meta_path)
  invisible(version_dir)
}

#' Create a new model version with optional template control stream
#'
#' @param copy_from Optional existing version id to duplicate (control stream only).
#' @param fit_inits Optional \code{nm_fit} object; when copying, replace
#'   \code{$THETA}/\code{$OMEGA}/\code{$SIGMA} initials with \code{fit_inits$par}.
#' @param data_file Optional \code{$DATA} path stored in version metadata.
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("demo", path = ws, template = "empty")
#' nm_workspace_new_version("demo", root = ws)
#' }
#' @export
nm_workspace_new_version <- function(project,
                                     version_id = NULL,
                                     root = nm_workspace_root(),
                                     template_ctl = NULL,
                                     label = NULL,
                                     data_file = NULL,
                                     copy_from = NULL,
                                     fit_inits = NULL) {
  existing <- nm_workspace_list_versions(project, root = root)
  if (is.null(version_id)) {
    n <- length(existing) + 1L
    version_id <- sprintf("mod%03d", n)
    while (version_id %in% existing) {
      n <- n + 1L
      version_id <- sprintf("mod%03d", n)
    }
  }
  from_meta <- list()
  if (!is.null(copy_from)) {
    copy_from <- .nm_ws_sanitize_name(copy_from)
    if (!copy_from %in% existing) {
      .nm_stop("Model version not found: ", copy_from)
    }
    template_ctl <- nm_ctl_compose(
      nm_ctl_parse(nm_workspace_read_model(project, copy_from, root = root))
    )
    if (!is.null(fit_inits)) {
      parts <- nm_ctl_parse(template_ctl)
      parts <- nm_ctl_apply_fit_inits(parts, fit_inits)
      template_ctl <- nm_ctl_compose(parts)
    }
    from_meta <- .nm_ws_read_json(
      file.path(.nm_ws_version_dir(project, copy_from, root = root), "meta.json")
    )
    if (is.null(label) || !nzchar(trimws(as.character(label)))) {
      label <- paste("Copy of", copy_from)
    }
    if (is.null(data_file) && !is.null(from_meta$data_file)) {
      data_file <- from_meta$data_file
    }
  }
  if (is.null(template_ctl)) {
    proj_dir <- .nm_ws_project_dir(project, root)
    meta <- .nm_ws_read_json(file.path(proj_dir, "project.json"))
    if (identical(meta$template, "theo") || file.exists(file.path(proj_dir, "data", "theo.csv"))) {
      sim <- nm_synthetic_theo(n_sub = 10L, seed = 1L)
      template_ctl <- .nm_model_to_ctl(
        sim$model,
        data_file = "data/theo.csv",
        prob = "THEO demo",
        method = "FOCE"
      )
      if (is.null(label)) {
        label <- "THEO base model"
      }
      if (is.null(data_file)) {
        data_file <- "data/theo.csv"
      }
    } else {
      parts <- nm_ctl_template(
        advan = 2L,
        trans = 2L,
        data_file = if (is.null(data_file)) "data.csv" else data_file,
        problem = "New model"
      )
      template_ctl <- nm_ctl_compose(parts)
    }
  }
  if (is.null(data_file)) {
    proj_dir <- .nm_ws_project_dir(project, root)
    if (file.exists(file.path(proj_dir, "data", "theo.csv"))) {
      data_file <- "data/theo.csv"
    }
  }
  nm_workspace_write_model(
    project, version_id, template_ctl, root = root,
    data_file = data_file,
    label = label
  )
  version_id
}

#' Duplicate an existing model version
#'
#' @param from_version_id Source model version id.
#' @param fit_inits Optional \code{nm_fit} object; when copying, replace
#'   \code{$THETA}/\code{$OMEGA}/\code{$SIGMA} initials with \code{fit_inits$par}.
#' @inheritParams nm_workspace_new_version
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("theo_demo", path = ws, template = "theo")
#' ver <- nm_workspace_list_versions("theo_demo", root = ws)[[1L]]
#' nm_workspace_copy_version("theo_demo", ver, root = ws)
#' }
#' @export
nm_workspace_copy_version <- function(project,
                                      from_version_id,
                                      version_id = NULL,
                                      root = nm_workspace_root(),
                                      label = NULL,
                                      fit_inits = NULL) {
  nm_workspace_new_version(
    project,
    version_id = version_id,
    root = root,
    copy_from = from_version_id,
    label = label,
    fit_inits = fit_inits
  )
}

#' @rdname nm_workspace_new_version
#' @export
nm_workspace_new_model <- function(project,
                                   run_id = NULL,
                                   root = nm_workspace_root(),
                                   template_ctl = NULL) {
  nm_workspace_new_version(
    project,
    version_id = run_id,
    root = root,
    template_ctl = template_ctl
  )
}

#' List estimation runs for a model version
#'
#' Estimation runs are fits of the same model version (e.g. FO vs FOCE attempts).
#'
#' @param project Project name.
#' @param version_id Model version id.
#' @param root Workspace root.
#' @return Data frame with \code{run_id}, \code{method}, \code{objective}, \code{label}, \code{created}.
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("theo_demo", path = ws, template = "theo")
#' ver <- nm_workspace_list_versions("theo_demo", root = ws)[[1L]]
#' nm_workspace_list_runs("theo_demo", ver, root = ws)
#' }
#' @export
nm_workspace_list_runs <- function(project, version_id, root = nm_workspace_root()) {
  empty <- data.frame(
    run_id = character(),
    method = character(),
    objective = numeric(),
    label = character(),
    created = character(),
    has_bootstrap = logical(),
    has_npc = logical(),
    has_npde = logical(),
    stringsAsFactors = FALSE
  )
  if (!.nm_ws_valid_name(project) || !.nm_ws_valid_name(version_id)) {
    return(empty)
  }
  version_id <- .nm_ws_sanitize_name(version_id)
  version_rows <- .nm_ws_collect_runs_from_dir(
    .nm_ws_version_runs_dir(project, version_id, root = root),
    version_id
  )
  legacy_rows <- .nm_ws_collect_runs_from_dir(
    .nm_ws_legacy_runs_dir(project, root),
    version_id
  )
  rows <- list()
  if (nrow(version_rows) > 0L) {
    rows[[length(rows) + 1L]] <- version_rows
  }
  if (nrow(legacy_rows) > 0L) {
    legacy_only <- legacy_rows[!legacy_rows$run_id %in% version_rows$run_id, , drop = FALSE]
    if (nrow(legacy_only) > 0L) {
      rows[[length(rows) + 1L]] <- legacy_only
    }
  }
  if (length(rows) == 0L) {
    legacy <- nm_workspace_fit_path(project, version_id, root = root)
    if (file.exists(legacy)) {
      fit <- tryCatch(readRDS(legacy), error = function(e) NULL)
      if (!is.null(fit)) {
        rows[[1L]] <- data.frame(
          run_id = paste0(version_id, "_fit"),
          method = as.character(if (is.null(fit$method)) "" else fit$method),
          objective = as.numeric(if (is.null(fit$objective)) NA_real_ else fit$objective),
          label = "legacy fit",
          created = "",
          has_bootstrap = !is.null(fit$bootstrap),
          has_npc = .nm_fit_has_npc(fit),
          has_npde = .nm_fit_has_npde(fit),
          stringsAsFactors = FALSE
        )
      }
    }
    if (length(rows) == 0L) {
      return(empty)
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$created, out$run_id, decreasing = TRUE), , drop = FALSE]
}

#' Allocate a new estimation run id for a model version
#'
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("theo_demo", path = ws, template = "theo")
#' ver <- nm_workspace_list_versions("theo_demo", root = ws)[[1L]]
#' nm_workspace_new_run_id("theo_demo", ver, root = ws)
#' }
#' @export
nm_workspace_new_run_id <- function(project, version_id, root = nm_workspace_root()) {
  runs_dir <- .nm_ws_version_runs_dir(project, version_id, root = root, create = TRUE)
  # A run id must be unique per *submission*, not per saved run. Remote/async runs
  # are only written to disk when they finish, so several jobs submitted before
  # any of them completes would otherwise all be handed the same "estNNN" id and
  # overwrite each other on save (leaving only a couple of runs visible). Reserve
  # the id atomically by creating its directory: dir.create() succeeds for the
  # winner and fails for anyone racing on the same name, so every caller walks
  # forward to a distinct id. Consider both saved runs and already-reserved dirs.
  taken <- unique(c(
    nm_workspace_list_runs(project, version_id, root = root)$run_id,
    {
      d <- list.dirs(runs_dir, full.names = FALSE, recursive = FALSE)
      d[nzchar(d)]
    }
  ))
  n <- length(taken) + 1L
  repeat {
    run_id <- sprintf("est%03d", n)
    if (!(run_id %in% taken)) {
      created <- suppressWarnings(
        dir.create(file.path(runs_dir, run_id), recursive = TRUE, showWarnings = FALSE)
      )
      if (isTRUE(created)) {
        return(run_id)
      }
      taken <- c(taken, run_id)
    }
    n <- n + 1L
  }
}

#' Release a run id reserved by nm_workspace_new_run_id() but never saved
#'
#' Only removes the directory if it is still empty (a pure reservation), so a
#' genuinely saved run is never deleted.
#' @keywords internal
.nm_ws_release_reserved_run <- function(project, version_id, run_id,
                                        root = nm_workspace_root()) {
  if (!.nm_ws_valid_name(project) || !.nm_ws_valid_name(version_id) ||
      is.null(run_id) || !nzchar(run_id)) {
    return(invisible(FALSE))
  }
  run_dir <- .nm_ws_run_dir(project, version_id, run_id, root = root, create = FALSE)
  if (!dir.exists(run_dir)) {
    return(invisible(FALSE))
  }
  if (length(list.files(run_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
    return(invisible(FALSE))
  }
  unlink(run_dir, recursive = TRUE)
  invisible(TRUE)
}

#' Save an estimation run (fit + metadata)
#'
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("demo", path = ws, template = "theo")
#' ver <- nm_workspace_list_versions("demo", root = ws)[[1L]]
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "FO",
#'               control = list(maxit = 5L, compute_inference = FALSE))
#' run_id <- nm_workspace_new_run_id("demo", ver, root = ws)
#' nm_workspace_save_run("demo", ver, run_id, fit, root = ws)
#' }
#' @export
nm_workspace_save_run <- function(project,
                                  version_id,
                                  run_id,
                                  fit,
                                  root = nm_workspace_root(),
                                  label = NULL,
                                  job_id = NULL,
                                  run_info = NULL) {
  version_id <- .nm_ws_sanitize_name(version_id)
  run_id <- .nm_ws_sanitize_name(run_id)
  run_dir <- .nm_ws_run_dir(project, version_id, run_id, root = root, create = TRUE)
  if (is.null(run_info) && !is.null(fit$run_info)) {
    run_info <- fit$run_info
  }
  if (!is.null(run_info)) {
    fit$run_info <- run_info
  }
  saveRDS(fit, file.path(run_dir, "fit.rds"))
  meta <- list(
    run_id = run_id,
    version_id = version_id,
    method = fit$method,
    objective = fit$objective,
    label = label,
    job_id = job_id,
    created = as.character(Sys.time()),
    has_bootstrap = !is.null(fit$bootstrap),
    bootstrap_n = if (!is.null(fit$bootstrap)) as.integer(fit$bootstrap$n_boot %||% NA_integer_) else NA_integer_,
    bootstrap_ok = if (!is.null(fit$bootstrap)) as.integer(fit$bootstrap$n_ok %||% NA_integer_) else NA_integer_,
    has_npc = .nm_fit_has_npc(fit),
    has_npde = .nm_fit_has_npde(fit),
    npc_n_sim = if (!is.null(fit$npc)) as.integer(fit$npc$n_sim %||% NA_integer_) else NA_integer_,
    npde_n_sim = if (!is.null(fit$npde)) as.integer(fit$npde$n_sim %||% NA_integer_) else NA_integer_,
    run_info = run_info
  )
  .nm_ws_write_json(meta, file.path(run_dir, "meta.json"))
  invisible(run_dir)
}

#' Load estimation run metadata
#'
#' @param project Project name.
#' @param version_id Model version id.
#' @param run_id Estimation run id.
#' @param root Workspace root.
#' @return Metadata list from \code{meta.json}, or \code{NULL}.
#' @export
nm_workspace_load_run_meta <- function(project,
                                       version_id,
                                       run_id,
                                       root = nm_workspace_root()) {
  version_id <- .nm_ws_sanitize_name(version_id)
  run_id <- .nm_ws_sanitize_name(run_id)
  run_dir <- .nm_ws_find_run_dir(project, version_id, run_id, root = root)
  if (is.null(run_dir)) {
    return(NULL)
  }
  meta_path <- file.path(run_dir, "meta.json")
  if (!file.exists(meta_path)) {
    return(NULL)
  }
  .nm_ws_read_json(meta_path)
}

#' Load fit for an estimation run
#'
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("demo", path = ws, template = "theo")
#' ver <- nm_workspace_list_versions("demo", root = ws)[[1L]]
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "FO",
#'               control = list(maxit = 5L, compute_inference = FALSE))
#' run_id <- nm_workspace_new_run_id("demo", ver, root = ws)
#' nm_workspace_save_run("demo", ver, run_id, fit, root = ws)
#' nm_workspace_load_run_fit("demo", ver, run_id, root = ws)
#' }
#' @export
nm_workspace_load_run_fit <- function(project,
                                      version_id,
                                      run_id,
                                      root = nm_workspace_root()) {
  version_id <- .nm_ws_sanitize_name(version_id)
  run_id <- .nm_ws_sanitize_name(run_id)
  run_dir <- .nm_ws_find_run_dir(project, version_id, run_id, root = root)
  if (!is.null(run_dir)) {
    run_path <- file.path(run_dir, "fit.rds")
    if (file.exists(run_path)) {
      return(readRDS(run_path))
    }
  }
  legacy_id <- paste0(version_id, "_fit")
  if (identical(run_id, legacy_id)) {
    legacy <- nm_workspace_fit_path(project, version_id, root = root)
    if (file.exists(legacy)) {
      return(readRDS(legacy))
    }
  }
  NULL
}

#' Resolve dataset path for a model version
#'
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("theo_demo", path = ws, template = "theo")
#' ver <- nm_workspace_list_versions("theo_demo", root = ws)[[1L]]
#' nm_workspace_model_data_path("theo_demo", ver, root = ws)
#' }
#' @export
nm_workspace_model_data_path <- function(project, version_id, root = nm_workspace_root()) {
  version_dir <- .nm_ws_version_dir(project, version_id, root = root)
  meta <- .nm_ws_read_json(file.path(version_dir, "meta.json"))
  data_rel <- meta$data_file
  if (is.null(data_rel) || !nzchar(data_rel)) {
    imported <- tryCatch(
      nm_workspace_parse_model(project, version_id, root = root),
      error = function(e) NULL
    )
    if (!is.null(imported) && !is.null(imported$data_path)) {
      data_rel <- imported$data_path
    }
  }
  if (is.null(data_rel) || !nzchar(data_rel)) {
    return(NULL)
  }
  proj_dir <- .nm_ws_project_dir(project, root)
  cand <- file.path(proj_dir, data_rel)
  if (file.exists(cand)) {
    return(normalizePath(cand, winslash = "/", mustWork = FALSE))
  }
  cand2 <- file.path(proj_dir, "data", basename(data_rel))
  if (file.exists(cand2)) {
    return(normalizePath(cand2, winslash = "/", mustWork = FALSE))
  }
  NULL
}

#' Parse a workspace model version into LibeRation objects
#'
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("theo_demo", path = ws, template = "theo")
#' ver <- nm_workspace_list_versions("theo_demo", root = ws)[[1L]]
#' nm_workspace_parse_model("theo_demo", ver, root = ws)$model$ADVAN
#' }
#' @export
nm_workspace_parse_model <- function(project, version_id, root = nm_workspace_root()) {
  ctl_path <- nm_workspace_model_path(project, version_id, root = root)
  proj_dir <- .nm_ws_project_dir(project, root)
  version_dir <- .nm_ws_version_dir(project, version_id, root = root)
  meta <- .nm_ws_read_json(file.path(version_dir, "meta.json"))
  data_override <- meta$data_file
  imported <- nm_read_nonmem(ctl_path, data_path = data_override)
  if (!is.null(imported$data_path)) {
    dp <- file.path(proj_dir, imported$data_path)
    if (!file.exists(dp)) {
      dp <- file.path(proj_dir, "data", basename(imported$data_path))
    }
    if (file.exists(dp)) {
      imported$data_path <- dp
      imported$data <- nm_dataset(dp)
    }
  }
  imported
}

#' Path for legacy per-version fit RDS (\code{fits/})
#'
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("theo_demo", path = ws, template = "theo")
#' ver <- nm_workspace_list_versions("theo_demo", root = ws)[[1L]]
#' nm_workspace_fit_path("theo_demo", ver, root = ws)
#' }
#' @export
nm_workspace_fit_path <- function(project, version_id, root = nm_workspace_root()) {
  version_id <- .nm_ws_sanitize_name(version_id)
  file.path(.nm_ws_project_dir(project, root), "fits", paste0(version_id, ".rds"))
}

#' Save fit to project (legacy: also writes under \code{fits/})
#'
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("demo", path = ws, template = "theo")
#' ver <- nm_workspace_list_versions("demo", root = ws)[[1L]]
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "FO",
#'               control = list(maxit = 5L, compute_inference = FALSE))
#' nm_workspace_save_fit("demo", ver, fit, root = ws)
#' }
#' @export
nm_workspace_save_fit <- function(project, version_id, fit, root = nm_workspace_root(), ...) {
  nm_workspace_save_run(
    project, version_id,
    run_id = nm_workspace_new_run_id(project, version_id, root = root),
    fit = fit,
    root = root,
    ...
  )
}

#' Load most recent fit for a model version
#'
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("demo", path = ws, template = "theo")
#' ver <- nm_workspace_list_versions("demo", root = ws)[[1L]]
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "FO",
#'               control = list(maxit = 5L, compute_inference = FALSE))
#' nm_workspace_save_fit("demo", ver, fit, root = ws)
#' nm_workspace_load_fit("demo", ver, root = ws)$method
#' }
#' @export
nm_workspace_load_fit <- function(project, version_id, root = nm_workspace_root()) {
  runs <- nm_workspace_list_runs(project, version_id, root = root)
  if (nrow(runs) == 0L) {
    return(NULL)
  }
  nm_workspace_load_run_fit(project, version_id, runs$run_id[[1L]], root = root)
}

#' Reports directory for a project (confined to workspace)
#'
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("demo", path = ws, template = "empty")
#' nm_workspace_reports_dir("demo", root = ws)
#' }
#' @export
nm_workspace_reports_dir <- function(project, root = nm_workspace_root()) {
  .nm_ws_resolve("projects", project, "reports", root = root)
}

#' List simulations for a model version
#'
#' @param project Project name.
#' @param version_id Model version id.
#' @param root Workspace root.
#' @return Data frame with \code{sim_id}, \code{label}, \code{seed}, \code{n_sim},
#'   \code{created}.
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("theo_demo", path = ws, template = "theo")
#' ver <- nm_workspace_list_versions("theo_demo", root = ws)[[1L]]
#' nm_workspace_list_sims("theo_demo", ver, root = ws)
#' }
#' @export
nm_workspace_list_sims <- function(project, version_id, root = nm_workspace_root()) {
  empty <- data.frame(
    sim_id = character(),
    label = character(),
    seed = integer(),
    n_sim = integer(),
    est_run_id = character(),
    vpc = logical(),
    created = character(),
    stringsAsFactors = FALSE
  )
  if (!.nm_ws_valid_name(project) || !.nm_ws_valid_name(version_id)) {
    return(empty)
  }
  version_id <- .nm_ws_sanitize_name(version_id)
  sims_dir <- file.path(.nm_ws_project_dir(project, root), "simulations")
  rows <- list()
  if (dir.exists(sims_dir)) {
    ids <- list.dirs(sims_dir, full.names = FALSE, recursive = FALSE)
    ids <- ids[nzchar(ids)]
    for (sid in ids) {
      meta <- .nm_ws_read_json(file.path(sims_dir, sid, "meta.json"))
      vid <- if (is.null(meta$version_id)) meta$model_version else meta$version_id
      if (!identical(vid, version_id)) {
        next
      }
      rows[[length(rows) + 1L]] <- data.frame(
        sim_id = sid,
        label = .nm_ws_meta_scalar(meta$label),
        seed = {
          s <- meta$seed
          if (is.null(s) || length(s) == 0L) {
            NA_integer_
          } else {
            as.integer(s[[1L]])
          }
        },
        n_sim = {
          n <- meta$n_sim
          if (is.null(n) || length(n) == 0L) {
            1L
          } else {
            as.integer(n[[1L]])
          }
        },
        est_run_id = .nm_ws_meta_scalar(meta$est_run_id),
        vpc = isTRUE(meta$vpc),
        created = .nm_ws_meta_scalar(meta$created),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0L) {
    return(empty)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  for (i in seq_len(nrow(out))) {
    if (!isTRUE(out$vpc[[i]]) && nzchar(out$est_run_id[[i]] %||% "")) {
      sim_path <- file.path(
        .nm_ws_project_dir(project, root), "simulations", out$sim_id[[i]], "sim_data.rds"
      )
      if (file.exists(sim_path)) {
        obj <- tryCatch(readRDS(sim_path), error = function(e) NULL)
        if (is.list(obj) && isTRUE(obj$vpc_mode)) {
          out$vpc[[i]] <- TRUE
        }
      }
    }
  }
  out[order(out$created, out$sim_id, decreasing = TRUE), , drop = FALSE]
}

#' List all simulations in a project (across model versions)
#'
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("theo_demo", path = ws, template = "theo")
#' nm_workspace_list_project_sims("theo_demo", root = ws)
#' }
#' @export
nm_workspace_list_project_sims <- function(project, root = nm_workspace_root()) {
  empty <- data.frame(
    sim_id = character(),
    version_id = character(),
    label = character(),
    seed = integer(),
    n_sim = integer(),
    created = character(),
    stringsAsFactors = FALSE
  )
  if (!.nm_ws_valid_name(project)) {
    return(empty)
  }
  sims_dir <- file.path(.nm_ws_project_dir(project, root), "simulations")
  if (!dir.exists(sims_dir)) {
    return(empty)
  }
  ids <- list.dirs(sims_dir, full.names = FALSE, recursive = FALSE)
  ids <- ids[nzchar(ids)]
  rows <- list()
  for (sid in ids) {
    meta <- .nm_ws_read_json(file.path(sims_dir, sid, "meta.json"))
    if (length(meta) == 0L) {
      next
    }
    vid <- if (is.null(meta$version_id)) meta$model_version else meta$version_id
    rows[[length(rows) + 1L]] <- data.frame(
      sim_id = sid,
      version_id = as.character(vid),
      label = .nm_ws_meta_scalar(meta$label),
      seed = as.integer(if (is.null(meta$seed)) NA_integer_ else meta$seed),
      n_sim = as.integer(if (is.null(meta$n_sim)) 1L else meta$n_sim),
      created = .nm_ws_meta_scalar(meta$created),
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0L) {
    return(empty)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$created, out$sim_id, decreasing = TRUE), , drop = FALSE]
}

#' List all estimation runs in a project (across model versions)
#'
#' @param project Project name.
#' @param root Workspace root.
#' @return Data frame with \code{run_id}, \code{version_id}, \code{method},
#'   \code{objective}, \code{label}, and \code{created}.
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("theo_demo", path = ws, template = "theo")
#' nm_workspace_list_project_runs("theo_demo", root = ws)
#' }
#' @export
nm_workspace_list_project_runs <- function(project, root = nm_workspace_root()) {
  empty <- data.frame(
    run_id = character(),
    version_id = character(),
    method = character(),
    objective = numeric(),
    label = character(),
    created = character(),
    has_bootstrap = logical(),
    has_npc = logical(),
    has_npde = logical(),
    stringsAsFactors = FALSE
  )
  if (!.nm_ws_valid_name(project)) {
    return(empty)
  }
  vers <- nm_workspace_list_versions(project, root = root)
  if (length(vers) == 0L) {
    return(empty)
  }
  rows <- list()
  for (ver in vers) {
    runs <- nm_workspace_list_runs(project, ver, root = root)
    if (nrow(runs) == 0L) {
      next
    }
    runs$version_id <- ver
    rows[[length(rows) + 1L]] <- runs
  }
  if (length(rows) == 0L) {
    return(empty)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$created, out$run_id, decreasing = TRUE), , drop = FALSE]
}

#' Allocate a new simulation id for a model version
#'
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("theo_demo", path = ws, template = "theo")
#' ver <- nm_workspace_list_versions("theo_demo", root = ws)[[1L]]
#' nm_workspace_new_sim_id("theo_demo", ver, root = ws)
#' }
#' @export
nm_workspace_new_sim_id <- function(project, version_id, root = nm_workspace_root()) {
  existing <- nm_workspace_list_sims(project, version_id, root = root)
  n <- if (nrow(existing) > 0L) nrow(existing) + 1L else 1L
  sim_id <- sprintf("sim%03d", n)
  while (sim_id %in% existing$sim_id) {
    n <- n + 1L
    sim_id <- sprintf("sim%03d", n)
  }
  sim_id
}

#' Save a simulation result
#'
#' @param sim_data An \code{nm_dataset}, or a list with \code{replicates} for
#'   multi-replicate simulations.
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("demo", path = ws, template = "theo")
#' ver <- nm_workspace_list_versions("demo", root = ws)[[1L]]
#' sim <- nm_synthetic_theo(n_sub = 2L)
#' sim_id <- nm_workspace_new_sim_id("demo", ver, root = ws)
#' nm_workspace_save_sim("demo", ver, sim_id, sim$data, root = ws)
#' }
#' @export
nm_workspace_save_sim <- function(project,
                                  version_id,
                                  sim_id,
                                  sim_data,
                                  root = nm_workspace_root(),
                                  label = NULL,
                                  seed = NULL,
                                  n_sim = 1L,
                                  use_fit = FALSE,
                                  est_run_id = NULL,
                                  vpc = FALSE) {
  version_id <- .nm_ws_sanitize_name(version_id)
  sim_id <- .nm_ws_sanitize_name(sim_id)
  sim_dir <- file.path(.nm_ws_project_dir(project, root), "simulations", sim_id)
  dir.create(sim_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(sim_data, file.path(sim_dir, "sim_data.rds"))
  vpc_flag <- isTRUE(vpc) ||
    (is.list(sim_data) && isTRUE(sim_data$vpc_mode))
  meta <- list(
    sim_id = sim_id,
    version_id = version_id,
    label = label,
    seed = seed,
    n_sim = as.integer(n_sim),
    use_fit = isTRUE(use_fit),
    est_run_id = est_run_id,
    vpc = vpc_flag,
    created = as.character(Sys.time())
  )
  .nm_ws_write_json(meta, file.path(sim_dir, "meta.json"))
  invisible(sim_dir)
}

#' Load a saved simulation
#'
#' Returns an \code{nm_dataset} for single simulations, or a list with
#' \code{replicates} and \code{primary} when multiple were saved.
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("demo", path = ws, template = "theo")
#' ver <- nm_workspace_list_versions("demo", root = ws)[[1L]]
#' sim <- nm_synthetic_theo(n_sub = 2L)
#' sim_id <- nm_workspace_new_sim_id("demo", ver, root = ws)
#' nm_workspace_save_sim("demo", ver, sim_id, sim$data, root = ws)
#' nm_workspace_load_sim("demo", ver, sim_id, root = ws)
#' }
#' @export
nm_workspace_load_sim <- function(project, version_id, sim_id, root = nm_workspace_root()) {
  version_id <- .nm_ws_sanitize_name(version_id)
  sim_id <- .nm_ws_sanitize_name(sim_id)
  sim_path <- file.path(
    .nm_ws_project_dir(project, root), "simulations", sim_id, "sim_data.rds"
  )
  if (!file.exists(sim_path)) {
    return(NULL)
  }
  obj <- readRDS(sim_path)
  if (is.list(obj) && !inherits(obj, "nm_dataset") && !is.null(obj$primary)) {
    return(obj)
  }
  obj
}

#' Load VPC summary stored for an estimation run
#'
#' Finds a VPC simulation linked to \code{est_run_id} and returns its summary.
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("demo", path = ws, template = "theo")
#' ver <- nm_workspace_list_versions("demo", root = ws)[[1L]]
#' run_id <- nm_workspace_new_run_id("demo", ver, root = ws)
#' nm_workspace_load_run_vpc("demo", ver, run_id, root = ws)
#' }
#' @export
nm_workspace_load_run_vpc <- function(project,
                                      version_id,
                                      run_id,
                                      root = nm_workspace_root()) {
  if (!.nm_ws_valid_name(project) || !.nm_ws_valid_name(version_id) ||
      !.nm_ws_valid_name(run_id)) {
    return(NULL)
  }
  sims <- nm_workspace_list_sims(project, version_id, root = root)
  if (nrow(sims) == 0L) {
    return(NULL)
  }
  hits <- sims[sims$vpc %in% TRUE & sims$est_run_id == run_id, , drop = FALSE]
  if (nrow(hits) == 0L) {
    return(NULL)
  }
  hits <- hits[order(hits$created, hits$sim_id, decreasing = TRUE), , drop = FALSE]
  sim_id <- hits$sim_id[[1L]]
  sim_obj <- nm_workspace_load_sim(project, version_id, sim_id, root = root)
  if (is.null(sim_obj) || !is.list(sim_obj) || !isTRUE(sim_obj$vpc_mode)) {
    return(NULL)
  }
  list(
    sim_id = sim_id,
    vpc = sim_obj$vpc,
    vpc_obs = sim_obj$vpc_obs
  )
}

#' Delete a simulation
#'
#' @examples
#' \dontrun{
#' ws <- tempfile()
#' nm_workspace_init(ws, create_demo_project = FALSE)
#' nm_workspace_create_project("demo", path = ws, template = "theo")
#' ver <- nm_workspace_list_versions("demo", root = ws)[[1L]]
#' sim <- nm_synthetic_theo(n_sub = 2L)
#' sim_id <- nm_workspace_new_sim_id("demo", ver, root = ws)
#' nm_workspace_save_sim("demo", ver, sim_id, sim$data, root = ws)
#' nm_workspace_delete_sim("demo", ver, sim_id, root = ws)
#' }
#' @export
nm_workspace_delete_sim <- function(project, version_id, sim_id, root = nm_workspace_root()) {
  sim_id <- .nm_ws_sanitize_name(sim_id)
  sim_dir <- file.path(.nm_ws_project_dir(project, root), "simulations", sim_id)
  if (dir.exists(sim_dir)) {
    unlink(sim_dir, recursive = TRUE)
    return(invisible(TRUE))
  }
  invisible(FALSE)
}

#' Extract primary dataset from a loaded simulation object
#'
#' @keywords internal
.nm_workspace_sim_dataset <- function(sim_obj) {
  if (is.null(sim_obj)) {
    return(NULL)
  }
  if (inherits(sim_obj, "nm_dataset")) {
    return(sim_obj)
  }
  if (is.list(sim_obj) && !is.null(sim_obj$combined)) {
    return(sim_obj$combined)
  }
  if (is.list(sim_obj) && !is.null(sim_obj$primary)) {
    return(sim_obj$primary)
  }
  if (is.list(sim_obj) && !is.null(sim_obj$data)) {
    return(structure(list(data = sim_obj$data), class = "nm_dataset"))
  }
  NULL
}
