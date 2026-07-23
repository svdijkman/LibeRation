.nm_workspace_path <- function(workspace) {
  path <- if (inherits(workspace, "nm_workspace")) workspace$path else workspace
  path <- as.character(path)
  if (length(path) != 1L || is.na(path) || !nzchar(path)) {
    .nm_stop("`workspace` must be one filesystem path or an nm_workspace.")
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.nm_workspace_component <- function(value, what = "project id") {
  value <- as.character(value)
  if (length(value) != 1L || is.na(value) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$", value) || value %in% c(".", "..")) {
    .nm_stop("Invalid ", what, ".")
  }
  value
}

.nm_workspace_atomic_save <- function(value, path) {
  directory <- dirname(path)
  if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE, showWarnings = FALSE)) {
    .nm_stop("Unable to create workspace directory: ", directory)
  }
  temporary <- tempfile("save-", tmpdir = directory, fileext = ".rds")
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  saveRDS(value, temporary, version = 3)
  if (file.exists(path)) {
    backup <- paste0(path, ".previous")
    unlink(backup, force = TRUE)
    if (!file.rename(path, backup)) .nm_stop("Unable to rotate workspace file: ", path)
    if (!file.rename(temporary, path)) {
      file.rename(backup, path)
      .nm_stop("Unable to publish workspace file: ", path)
    }
    unlink(backup, force = TRUE)
  } else if (!file.rename(temporary, path)) {
    .nm_stop("Unable to publish workspace file: ", path)
  }
  if (.Platform$OS.type != "windows") Sys.chmod(path, mode = "0600")
  invisible(path)
}

.nm_workspace_read <- function(path) {
  candidates <- c(path, paste0(path, ".previous"))
  candidates <- unique(candidates[file.exists(candidates)])
  if (!length(candidates)) .nm_stop("Workspace file does not exist: ", path)
  last <- NULL
  for (candidate in candidates) {
    value <- tryCatch(readRDS(candidate), error = function(error) {
      last <<- error
      NULL
    })
    if (!is.null(value)) {
      if (!identical(candidate, path)) {
        warning("Recovered interrupted workspace write from ",
                basename(candidate), ".", call. = FALSE)
      }
      return(value)
    }
  }
  .nm_stop(
    "Unable to read workspace file ", path, ": ",
    if (inherits(last, "condition")) conditionMessage(last) else "unknown read failure"
  )
}

.nm_workspace_now <- function() {
  format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")
}

.nm_workspace_lock <- function(workspace, name, timeout = 10, stale = 300) {
  root <- .nm_workspace_path(workspace)
  lock_root <- file.path(root, ".locks")
  if (!dir.exists(lock_root)) dir.create(lock_root, recursive = TRUE, showWarnings = FALSE)
  name <- .nm_workspace_component(name, "lock name")
  path <- file.path(lock_root, paste0(name, ".lock"))
  deadline <- Sys.time() + as.numeric(timeout)
  repeat {
    if (dir.create(path, showWarnings = FALSE)) {
      writeLines(c(paste0("pid=", Sys.getpid()), paste0("created=", .nm_workspace_now())),
                 file.path(path, "owner"), useBytes = TRUE)
      return(path)
    }
    info <- suppressWarnings(file.info(path))
    if (nrow(info) && !is.na(info$mtime) &&
        as.numeric(difftime(Sys.time(), info$mtime, units = "secs")) > stale) {
      unlink(path, recursive = TRUE, force = TRUE)
      next
    }
    if (Sys.time() >= deadline) .nm_stop("Timed out waiting for workspace lock: ", name, ".")
    Sys.sleep(0.025)
  }
}

.nm_workspace_unlock <- function(path) {
  if (!is.null(path) && dir.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
  invisible(TRUE)
}

.nm_object_reference <- function(hash, kind) {
  list(schema = "liberation.object-ref", version = 1L, hash = hash, kind = kind)
}

.nm_object_path <- function(workspace, hash) {
  hash <- as.character(hash)
  if (length(hash) != 1L || !grepl("^[a-f0-9]{64}$", hash)) {
    .nm_stop("Invalid workspace object hash.")
  }
  file.path(.nm_workspace_path(workspace), "objects", substr(hash, 1L, 2L),
            paste0(hash, ".rds"))
}

.nm_object_put <- function(workspace, value, kind) {
  if (is.null(value)) return(NULL)
  stored <- if (identical(kind, "model")) nm_model_to_contract(value) else value
  hash <- digest::digest(stored, algo = "sha256", serialize = TRUE)
  path <- .nm_object_path(workspace, hash)
  if (!file.exists(path)) {
    envelope <- list(
      schema = "liberation.object", version = 1L, hash = hash, kind = kind,
      created = .nm_workspace_now(), value = stored
    )
    .nm_workspace_atomic_save(envelope, path)
  }
  .nm_object_reference(hash, kind)
}

.nm_object_get <- function(workspace, reference) {
  if (is.null(reference)) return(NULL)
  if (!is.list(reference) ||
      !identical(as.character(reference$schema), "liberation.object-ref") ||
      !identical(as.integer(reference$version), 1L)) {
    .nm_stop("Invalid workspace object reference.")
  }
  path <- .nm_object_path(workspace, reference$hash)
  envelope <- .nm_workspace_read(path)
  if (!is.list(envelope) || !identical(envelope$schema, "liberation.object") ||
      !identical(envelope$hash, reference$hash) ||
      !identical(envelope$kind, reference$kind) ||
      !identical(digest::digest(envelope$value, algo = "sha256", serialize = TRUE),
                 reference$hash)) {
    .nm_stop("Workspace object failed integrity validation: ", reference$hash, ".")
  }
  if (identical(reference$kind, "model")) nm_model_from_contract(envelope$value) else envelope$value
}

.nm_snapshot_pack <- function(workspace, value) {
  model <- value$model %||% NULL
  data <- value$data %||% NULL
  result <- value$result %||% NULL
  hydrate_fit <- inherits(result, c("nm_fit", "nm_individual_fit"))
  if (hydrate_fit) {
    result <- result
    result$model <- NULL
    result$data <- NULL
  }
  value$model <- value$data <- value$result <- NULL
  value$version <- 2L
  value$schema <- "liberation.snapshot"
  value$storage <- list(
    schema = "liberation.snapshot-storage", version = 1L,
    model = .nm_object_put(workspace, model, "model"),
    data = .nm_object_put(workspace, data, "data"),
    result = .nm_object_put(workspace, result, "result"),
    hydrate_fit = hydrate_fit
  )
  value
}

.nm_snapshot_unpack <- function(workspace, value) {
  if (identical(as.integer(value$version %||% 1L), 1L)) return(value)
  if (!identical(value$schema %||% "", "liberation.snapshot") ||
      !is.list(value$storage) ||
      !identical(value$storage$schema %||% "", "liberation.snapshot-storage")) {
    .nm_stop("Project snapshot storage contract is invalid.")
  }
  model <- .nm_object_get(workspace, value$storage$model)
  data <- .nm_object_get(workspace, value$storage$data)
  result <- .nm_object_get(workspace, value$storage$result)
  if (isTRUE(value$storage$hydrate_fit) && !is.null(result)) {
    result$model <- model
    result$data <- data
  }
  value$model <- model
  value$data <- data
  value$result <- result
  value
}

.nm_manifest_upgrade_v1 <- function(manifest) {
  manifest$schema <- "liberation.project"
  manifest$version <- 2L
  manifest
}

.nm_manifest_migrate <- function(manifest) {
  version <- as.integer(manifest$version %||% 1L)
  if (length(version) != 1L || is.na(version) || version < 1L || version > 2L) {
    .nm_stop("Unsupported project-manifest schema version.")
  }
  while (version < 2L) {
    manifest <- switch(as.character(version), `1` = .nm_manifest_upgrade_v1(manifest))
    version <- as.integer(manifest$version)
  }
  if (!identical(manifest$schema %||% "", "liberation.project")) {
    .nm_stop("Project manifest schema is invalid.")
  }
  manifest
}

.nm_project_path <- function(workspace, project) {
  root <- .nm_workspace_path(workspace)
  project <- .nm_workspace_component(project)
  path <- file.path(root, "projects", project)
  normalized_parent <- normalizePath(dirname(path), winslash = "/", mustWork = TRUE)
  candidate <- file.path(normalized_parent, basename(path))
  if (!startsWith(tolower(candidate), paste0(tolower(root), "/"))) {
    .nm_stop("Resolved project path escaped the workspace.")
  }
  candidate
}

.nm_project_manifest <- function(workspace, project) {
  manifest <- .nm_manifest_migrate(.nm_workspace_read(file.path(
    .nm_project_path(workspace, project), "manifest.rds"
  )))
  snapshots <- manifest$snapshots
  had_covariance <- "has_covariance" %in% names(snapshots)
  required <- list(
    result_type = character(nrow(snapshots)),
    method = character(nrow(snapshots)),
    entry_type = rep("version", nrow(snapshots)),
    parent_id = character(nrow(snapshots)),
    run_number = rep(NA_integer_, nrow(snapshots)),
    has_vpc = rep(FALSE, nrow(snapshots)),
    has_gof = rep(FALSE, nrow(snapshots)),
    has_npc = rep(FALSE, nrow(snapshots)),
    has_npde = rep(FALSE, nrow(snapshots)),
    has_vpc_categorical = rep(FALSE, nrow(snapshots)),
    has_vpc_count = rep(FALSE, nrow(snapshots)),
    has_vpc_tte = rep(FALSE, nrow(snapshots)),
    has_vpc_competing = rep(FALSE, nrow(snapshots)),
    has_vpc_recurrent = rep(FALSE, nrow(snapshots)),
    has_bootstrap = rep(FALSE, nrow(snapshots)),
    has_profile = rep(FALSE, nrow(snapshots)),
    has_scm = rep(FALSE, nrow(snapshots)),
    has_covariance = rep(FALSE, nrow(snapshots)),
    queue_id = character(nrow(snapshots)),
    queue_job_id = character(nrow(snapshots))
  )
  for (name in names(required)) {
    if (!name %in% names(snapshots)) snapshots[[name]] <- required[[name]]
  }
  if (!had_covariance && nrow(snapshots)) {
    candidates <- which(
      snapshots$entry_type == "run" & snapshots$result_type == "estimation"
    )
    for (index in candidates) {
      stored <- tryCatch(
        .nm_snapshot_unpack(workspace, .nm_workspace_read(file.path(
          .nm_project_path(workspace, project), "snapshots",
          paste0(snapshots$id[[index]], ".rds")
        ))),
        error = function(error) NULL
      )
      covariance <- stored$result$covariance %||% NULL
      snapshots$has_covariance[[index]] <- !is.null(covariance) &&
        !identical(covariance$status %||% "completed", "failed")
    }
  }
  # Workspaces written by the prototype stored completed runs as peer snapshots.
  # Upgrade those records in memory by attaching each run to the preceding model
  # version. The upgraded manifest is persisted on the next mutation.
  if (!"entry_type" %in% names(manifest$snapshots) && nrow(snapshots)) {
    parent <- ""
    run_counts <- list()
    for (index in seq_len(nrow(snapshots))) {
      is_run <- isTRUE(snapshots$has_result[[index]]) &&
        snapshots$result_type[[index]] %in% c("estimation", "simulation") && nzchar(parent)
      if (is_run) {
        snapshots$entry_type[[index]] <- "run"
        snapshots$parent_id[[index]] <- parent
        count <- (run_counts[[parent]] %||% 0L) + 1L
        run_counts[[parent]] <- count
        snapshots$run_number[[index]] <- count
      } else {
        parent <- snapshots$id[[index]]
      }
    }
  }
  manifest$snapshots <- snapshots
  manifest
}

.nm_snapshot_result_type <- function(result) {
  if (is.null(result)) return("")
  if (inherits(result, "nm_fit")) return("estimation")
  if (inherits(result, "nm_vpc")) return("vpc")
  if (inherits(result, "nm_npc")) return("npc")
  if (inherits(result, "nm_npde")) return("npde")
  if (inherits(result, "nm_bootstrap")) return("bootstrap")
  if (is.data.frame(result) && "IPRED" %in% names(result)) return("simulation")
  "result"
}

#' Create or open a modelling workspace
#'
#' A workspace stores serializable model versions, their runs, and saved
#' diagnostics. Compiled external pointers are deliberately excluded and rebuilt
#' when a model version or run is opened.
#'
#' @param path Workspace directory.
#' @return An `nm_workspace` object.
#' @examples
#' workspace <- nm_workspace(tempfile("liber-workspace-"))
#' workspace
#' @export
nm_workspace <- function(path = file.path(tools::R_user_dir("LibeRation", "data"), "workspace")) {
  path <- path.expand(as.character(path))
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    .nm_stop("Unable to create workspace: ", path)
  }
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  directories <- file.path(path, c("projects", "reports", "objects", "backups", ".locks"))
  for (directory in directories) {
    if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE, showWarnings = FALSE)) {
      .nm_stop("Unable to create workspace directory: ", directory, ".")
    }
  }
  projects <- directories[[1L]]
  reports <- directories[[2L]]
  metadata_path <- file.path(path, "workspace.rds")
  if (file.exists(metadata_path)) {
    metadata <- .nm_workspace_read(metadata_path)
    if (!is.list(metadata) || !identical(metadata$schema %||% "", "liberation.workspace") ||
        !as.integer(metadata$version %||% 0L) %in% c(1L, 2L)) {
      .nm_stop("Workspace metadata is invalid or from an unsupported version.")
    }
    if (as.integer(metadata$version) < 2L) {
      metadata$version <- 2L
      metadata$migrated <- .nm_workspace_now()
      .nm_workspace_atomic_save(metadata, metadata_path)
    }
  } else {
    metadata <- list(schema = "liberation.workspace", version = 2L,
                     created = .nm_workspace_now(), updated = .nm_workspace_now())
    .nm_workspace_atomic_save(metadata, metadata_path)
  }
  if (.Platform$OS.type != "windows") {
    Sys.chmod(path, mode = "0700")
    for (directory in directories) Sys.chmod(directory, mode = "0700")
  }
  structure(list(version = 2L, path = path), class = "nm_workspace")
}

#' Create a project in a LibeRation workspace
#'
#' @param workspace Workspace object or path.
#' @param name Human-readable project name.
#' @param id Optional stable project identifier.
#' @param description Optional project description.
#' @return Project metadata.
#' @examples
#' workspace <- nm_workspace(tempfile("liber-workspace-"))
#' nm_project_create(workspace, "Example project")
#' @export
nm_project_create <- function(workspace, name, id = NULL, description = NULL) {
  workspace <- if (inherits(workspace, "nm_workspace")) workspace else nm_workspace(workspace)
  lock <- .nm_workspace_lock(workspace, "projects")
  on.exit(.nm_workspace_unlock(lock), add = TRUE)
  name <- trimws(as.character(name))
  if (length(name) != 1L || is.na(name) || !nzchar(name)) .nm_stop("Project name must not be empty.")
  if (is.null(id)) {
    stem <- tolower(gsub("(^-+|-+$)", "", gsub("[^A-Za-z0-9]+", "-", iconv(name, to = "ASCII//TRANSLIT"))))
    if (!nzchar(stem)) stem <- "project"
    id <- substr(stem, 1L, 96L)
    if (dir.exists(file.path(workspace$path, "projects", id))) {
      id <- paste0(id, "-", format(Sys.time(), "%Y%m%d%H%M%S"))
    }
  }
  id <- .nm_workspace_component(id)
  description <- trimws(as.character(description %||% ""))
  if (length(description) != 1L || is.na(description)) {
    .nm_stop("`description` must be one character value.")
  }
  directory <- .nm_project_path(workspace, id)
  if (dir.exists(directory)) .nm_stop("Project already exists: ", id, ".")
  if (!dir.create(file.path(directory, "snapshots"), recursive = TRUE, showWarnings = FALSE)) {
    .nm_stop("Unable to create project: ", id, ".")
  }
  if (!dir.create(file.path(directory, "diagnostics"), showWarnings = FALSE)) {
    .nm_stop("Unable to create project diagnostics directory: ", id, ".")
  }
  now <- format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")
  manifest <- list(
    schema = "liberation.project", version = 2L,
    id = id, name = name, description = description,
    created = now, updated = now,
    snapshots = data.frame(
      id = character(), label = character(), created = character(),
      has_model = logical(), has_data = logical(), has_result = logical(),
      result_type = character(), method = character(),
      entry_type = character(), parent_id = character(), run_number = integer(),
      has_vpc = logical(), has_npc = logical(), has_npde = logical(),
      has_gof = logical(),
      has_vpc_categorical = logical(), has_vpc_count = logical(),
      has_vpc_tte = logical(), has_vpc_competing = logical(),
      has_vpc_recurrent = logical(),
      has_bootstrap = logical(), has_profile = logical(), has_scm = logical(),
      has_covariance = logical(),
      stringsAsFactors = FALSE
    )
  )
  .nm_workspace_atomic_save(manifest, file.path(directory, "manifest.rds"))
  if (.Platform$OS.type != "windows") Sys.chmod(directory, mode = "0700")
  manifest
}

#' List projects or stored model versions and runs
#'
#' @param workspace Workspace object or path.
#' @param project Optional project id. When supplied, model versions and runs are
#'   returned. Use the `entry_type` and `parent_id` columns to reconstruct the
#'   hierarchy.
#' @return A data frame.
#' @export
nm_project_list <- function(workspace, project = NULL) {
  workspace <- if (inherits(workspace, "nm_workspace")) workspace else nm_workspace(workspace)
  if (!is.null(project)) {
    snapshots <- .nm_project_manifest(workspace, project)$snapshots
    if (!nrow(snapshots)) return(snapshots)
    snapshots$version <- seq_len(nrow(snapshots))
    return(snapshots[rev(seq_len(nrow(snapshots))), , drop = FALSE])
  }
  directories <- list.dirs(file.path(workspace$path, "projects"), recursive = FALSE, full.names = TRUE)
  records <- lapply(directories, function(directory) {
    tryCatch(
      .nm_project_manifest(workspace, basename(directory)),
      error = function(e) NULL
    )
  })
  records <- Filter(Negate(is.null), records)
  if (!length(records)) return(data.frame(
    id = character(), name = character(), description = character(),
    created = character(), updated = character(),
    snapshots = integer(), stringsAsFactors = FALSE
  ))
  output <- do.call(rbind, lapply(records, function(record) data.frame(
    id = record$id, name = record$name, created = record$created,
    description = record$description %||% "", updated = record$updated,
    snapshots = nrow(record$snapshots),
    versions = sum((record$snapshots$entry_type %||% rep("version", nrow(record$snapshots))) == "version"),
    stringsAsFactors = FALSE
  )))
  output[order(output$updated, decreasing = TRUE), , drop = FALSE]
}

#' Save a model version
#'
#' @param workspace Workspace object or path.
#' @param project Project id.
#' @param model Optional `nm_model` or `NMEngine`.
#' @param data Optional event dataset.
#' @param result Optional legacy result. New code should store completed work with
#'   [nm_project_save_run()] so it is nested beneath a model version.
#' @param label Model-version label. Empty labels are assigned `Mod001`,
#'   `Mod002`, and so on.
#' @param provenance Optional named provenance fields supplied by an integrated
#'   LibeR package. Core runtime provenance cannot be overwritten.
#' @return Model-version id, invisibly.
#' @export
nm_project_save <- function(workspace, project, model = NULL, data = NULL,
                            result = NULL, label = NULL, provenance = NULL) {
  workspace <- if (inherits(workspace, "nm_workspace")) workspace else nm_workspace(workspace)
  project <- .nm_workspace_component(project)
  lock <- .nm_workspace_lock(workspace, paste0("project-", project))
  on.exit(.nm_workspace_unlock(lock), add = TRUE)
  directory <- .nm_project_path(workspace, project)
  manifest <- .nm_project_manifest(workspace, project)
  if (inherits(model, "NMEngine")) model <- model$model
  if (!is.null(model) && !inherits(model, "nm_model")) .nm_stop("`model` must be an nm_model or NMEngine.")
  if (!is.null(data)) data <- if (inherits(data, "nm_dataset")) data else nm_dataset(data)
  if (!is.null(provenance) && (!is.list(provenance) ||
      (length(provenance) && is.null(names(provenance))))) {
    .nm_stop("`provenance` must be a named list.")
  }
  label <- trimws(as.character(label %||% ""))
  if (length(label) != 1L || is.na(label) || !nzchar(label)) {
    label <- sprintf("Mod%03d", sum(manifest$snapshots$entry_type == "version") + 1L)
  }
  id <- paste0(format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"), "-",
               sprintf("%08x", sample.int(.Machine$integer.max, 1L)))
  id <- .nm_workspace_component(id, "snapshot id")
  now <- format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")
  snapshot <- list(
    version = 2L, id = id, project = manifest$id, label = label, created = now,
    model = model, data = data, result = result,
    provenance = utils::modifyList(provenance %||% list(), list(
      R = R.version.string, platform = R.version$platform,
      LibeRation = as.character(utils::packageVersion("LibeRation"))
    ))
  )
  snapshot <- .nm_snapshot_pack(workspace, snapshot)
  .nm_workspace_atomic_save(snapshot, file.path(directory, "snapshots", paste0(id, ".rds")))
  manifest$snapshots <- rbind(manifest$snapshots, data.frame(
    id = id, label = label, created = now,
    has_model = !is.null(model), has_data = !is.null(data), has_result = !is.null(result),
    result_type = .nm_snapshot_result_type(result),
    method = if (inherits(result, "nm_fit")) .nm_fit_method_label(result) else "",
    entry_type = "version", parent_id = "", run_number = NA_integer_,
    has_vpc = FALSE, has_npc = FALSE, has_npde = FALSE, has_gof = FALSE,
    has_vpc_categorical = FALSE, has_vpc_count = FALSE,
    has_vpc_tte = FALSE, has_vpc_competing = FALSE,
    has_vpc_recurrent = FALSE,
    has_bootstrap = FALSE, has_profile = FALSE, has_scm = FALSE,
    has_covariance = FALSE, queue_id = "", queue_job_id = "",
    stringsAsFactors = FALSE
  ))
  manifest$updated <- now
  .nm_workspace_atomic_save(manifest, file.path(directory, "manifest.rds"))
  invisible(id)
}

#' Save a numbered model run
#'
#' A run is an estimation or simulation result owned by one immutable model
#' version. The model and data are copied into the run record so old work remains
#' reproducible even if the active editor changes later.
#'
#' @param workspace Workspace object or path.
#' @param project Project id.
#' @param version Parent model-version id.
#' @param result An `nm_fit` or simulation data frame.
#' @param label Optional human-readable run label.
#' @param model,data Optional run inputs; by default they are loaded from `version`.
#' @param queue_id,queue_job_id Optional execution provenance used to reconcile
#'   persistent queued jobs without creating duplicate model runs.
#' @return Run id, invisibly.
#' @export
nm_project_save_run <- function(workspace, project, version, result, label = NULL,
                                model = NULL, data = NULL,
                                queue_id = NULL, queue_job_id = NULL) {
  workspace <- if (inherits(workspace, "nm_workspace")) workspace else nm_workspace(workspace)
  project <- .nm_workspace_component(project)
  lock <- .nm_workspace_lock(workspace, paste0("project-", project))
  on.exit(.nm_workspace_unlock(lock), add = TRUE)
  manifest <- .nm_project_manifest(workspace, project)
  version <- .nm_workspace_component(version, "model version id")
  parent_index <- match(version, manifest$snapshots$id)
  if (is.na(parent_index) || manifest$snapshots$entry_type[[parent_index]] != "version") {
    .nm_stop("Unknown model version: ", version, ".")
  }
  kind <- .nm_snapshot_result_type(result)
  if (!kind %in% c("estimation", "simulation")) {
    .nm_stop("A model run must contain an estimation or simulation result.")
  }
  source <- nm_project_load(workspace, project, version)
  model <- model %||% source$model
  data <- data %||% source$data
  if (inherits(model, "NMEngine")) model <- model$model
  if (!is.null(model) && !inherits(model, "nm_model")) .nm_stop("`model` must be an nm_model or NMEngine.")
  if (!is.null(data)) data <- if (inherits(data, "nm_dataset")) data else nm_dataset(data)
  existing <- manifest$snapshots$entry_type == "run" & manifest$snapshots$parent_id == version
  number <- if (any(existing)) max(manifest$snapshots$run_number[existing], na.rm = TRUE) + 1L else 1L
  if (!is.finite(number)) number <- sum(existing) + 1L
  label <- trimws(as.character(label %||% ""))
  if (length(label) != 1L || is.na(label) || !nzchar(label)) {
    label <- if (kind == "estimation") sprintf("Estimation %03d", number) else sprintf("Simulation %03d", number)
  }
  id <- paste0(format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"), "-",
               sprintf("%08x", sample.int(.Machine$integer.max, 1L)))
  id <- .nm_workspace_component(id, "run id")
  now <- format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")
  run <- list(
    version = 2L, id = id, project = manifest$id, label = label, created = now,
    entry_type = "run", parent_id = version, run_number = number,
    model = model, data = data, result = result,
    provenance = list(
      R = R.version.string, platform = R.version$platform,
      LibeRation = as.character(utils::packageVersion("LibeRation")),
      queue_id = as.character(queue_id %||% ""),
      queue_job_id = as.character(queue_job_id %||% "")
    )
  )
  directory <- .nm_project_path(workspace, project)
  run <- .nm_snapshot_pack(workspace, run)
  .nm_workspace_atomic_save(run, file.path(directory, "snapshots", paste0(id, ".rds")))
  manifest$snapshots <- rbind(manifest$snapshots, data.frame(
    id = id, label = label, created = now,
    has_model = !is.null(model), has_data = !is.null(data), has_result = TRUE,
    result_type = kind,
    method = if (inherits(result, "nm_fit")) .nm_fit_method_label(result) else "",
    entry_type = "run", parent_id = version, run_number = as.integer(number),
    has_vpc = FALSE, has_npc = FALSE, has_npde = FALSE, has_gof = FALSE,
    has_vpc_categorical = FALSE, has_vpc_count = FALSE,
    has_vpc_tte = FALSE, has_vpc_competing = FALSE,
    has_vpc_recurrent = FALSE,
    has_bootstrap = FALSE, has_profile = FALSE, has_scm = FALSE,
    has_covariance = inherits(result, "nm_fit") &&
      !is.null(result$covariance) &&
      !identical(result$covariance$status %||% "completed", "failed"),
    queue_id = as.character(queue_id %||% ""),
    queue_job_id = as.character(queue_job_id %||% ""),
    stringsAsFactors = FALSE
  ))
  manifest$updated <- now
  .nm_workspace_atomic_save(manifest, file.path(directory, "manifest.rds"))
  invisible(id)
}

#' Save diagnostics for an estimation run
#'
#' @param workspace Workspace object or path.
#' @param project Project id.
#' @param run Estimation-run id.
#' @param diagnostics Named list containing saved plot data or secondary
#'   analyses for the run.
#' @return The merged diagnostic list, invisibly.
#' @export
nm_project_save_diagnostics <- function(workspace, project, run, diagnostics) {
  workspace <- if (inherits(workspace, "nm_workspace")) workspace else nm_workspace(workspace)
  project <- .nm_workspace_component(project)
  lock <- .nm_workspace_lock(workspace, paste0("project-", project))
  on.exit(.nm_workspace_unlock(lock), add = TRUE)
  run <- .nm_workspace_component(run, "run id")
  manifest <- .nm_project_manifest(workspace, project)
  index <- match(run, manifest$snapshots$id)
  if (is.na(index) || manifest$snapshots$entry_type[[index]] != "run" ||
      manifest$snapshots$result_type[[index]] != "estimation") {
    .nm_stop("Diagnostics can only be attached to an estimation run.")
  }
  if (!is.list(diagnostics) || is.null(names(diagnostics))) {
    .nm_stop("`diagnostics` must be a named list.")
  }
  allowed <- c(
    "gof", "vpc", "npc", "npde", "vpc_categorical", "vpc_count",
    "vpc_tte", "vpc_competing", "vpc_recurrent",
    "bootstrap", "profile", "scm"
  )
  if (length(setdiff(names(diagnostics), allowed))) .nm_stop("Unknown diagnostic type.")
  current <- nm_project_load_diagnostics(workspace, project, run)
  merged <- utils::modifyList(current, diagnostics)
  directory <- file.path(.nm_project_path(workspace, project), "diagnostics")
  if (!dir.exists(directory) && !dir.create(directory, showWarnings = FALSE)) {
    .nm_stop("Unable to create project diagnostics directory.")
  }
  .nm_workspace_atomic_save(merged, file.path(directory, paste0(run, ".rds")))
  for (name in allowed) {
    manifest$snapshots[[paste0("has_", name)]][[index]] <- !is.null(merged[[name]])
  }
  manifest$updated <- format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")
  .nm_workspace_atomic_save(manifest, file.path(.nm_project_path(workspace, project), "manifest.rds"))
  invisible(merged)
}

#' Load saved diagnostics for a model run
#' @param workspace Workspace object or path.
#' @param project Project id.
#' @param run Run id.
#' @return A named list of saved diagnostic objects.
#' @export
nm_project_load_diagnostics <- function(workspace, project, run) {
  workspace <- if (inherits(workspace, "nm_workspace")) workspace else nm_workspace(workspace)
  run <- .nm_workspace_component(run, "run id")
  path <- file.path(.nm_project_path(workspace, project), "diagnostics", paste0(run, ".rds"))
  if (!file.exists(path)) return(list())
  value <- .nm_workspace_read(path)
  if (!is.list(value)) .nm_stop("Saved run diagnostics are invalid.")
  value
}

#' Load a project snapshot
#'
#' @param workspace Workspace object or path.
#' @param project Project id.
#' @param snapshot Snapshot id or `"latest"`.
#' @return A serializable snapshot list.
#' @export
nm_project_load <- function(workspace, project, snapshot = "latest") {
  workspace <- if (inherits(workspace, "nm_workspace")) workspace else nm_workspace(workspace)
  manifest <- .nm_project_manifest(workspace, project)
  if (!nrow(manifest$snapshots)) .nm_stop("Project has no snapshots: ", project, ".")
  if (identical(snapshot, "latest")) snapshot <- tail(manifest$snapshots$id, 1L)
  snapshot <- .nm_workspace_component(snapshot, "snapshot id")
  if (!snapshot %in% manifest$snapshots$id) .nm_stop("Unknown snapshot: ", snapshot, ".")
  value <- .nm_snapshot_unpack(workspace, .nm_workspace_read(file.path(
    .nm_project_path(workspace, project), "snapshots", paste0(snapshot, ".rds")
  )))
  if (!is.list(value) || !as.integer(value$version %||% 0L) %in% c(1L, 2L) ||
      !identical(value$project, manifest$id)) {
    .nm_stop("Project snapshot is invalid.")
  }
  value
}

#' Copy a project snapshot into a new model version
#'
#' @param workspace Workspace object or path.
#' @param project Source project id.
#' @param snapshot Source snapshot id or `"latest"`.
#' @param target_project Destination project id; defaults to `project`.
#' @param label Optional label for the copied version.
#' @return The new snapshot id, invisibly.
#' @export
nm_project_copy <- function(workspace, project, snapshot = "latest",
                            target_project = project, label = NULL) {
  source <- nm_project_load(workspace, project, snapshot)
  if (is.null(label)) label <- paste(source$label, "copy")
  nm_project_save(
    workspace, target_project, source$model, source$data, NULL,
    label = label
  )
}

#' Delete a project snapshot
#'
#' @param workspace Workspace object or path.
#' @param project Project id.
#' @param snapshot Snapshot id.
#' @return `TRUE`, invisibly.
#' @export
nm_project_delete_snapshot <- function(workspace, project, snapshot) {
  workspace <- if (inherits(workspace, "nm_workspace")) workspace else nm_workspace(workspace)
  project <- .nm_workspace_component(project)
  lock <- .nm_workspace_lock(workspace, paste0("project-", project))
  on.exit(.nm_workspace_unlock(lock), add = TRUE)
  snapshot <- .nm_workspace_component(snapshot, "snapshot id")
  directory <- .nm_project_path(workspace, project)
  manifest <- .nm_project_manifest(workspace, project)
  index <- match(snapshot, manifest$snapshots$id)
  if (is.na(index)) .nm_stop("Unknown snapshot: ", snapshot, ".")
  targets <- snapshot
  if (manifest$snapshots$entry_type[[index]] == "version") {
    targets <- c(targets, manifest$snapshots$id[
      manifest$snapshots$entry_type == "run" & manifest$snapshots$parent_id == snapshot
    ])
  }
  for (target in targets) {
    path <- file.path(directory, "snapshots", paste0(target, ".rds"))
    status <- if (file.exists(path)) unlink(path, force = TRUE) else 0L
    if (!identical(status, 0L) || file.exists(path)) {
      .nm_stop("Unable to delete model version or run: ", target, ".")
    }
    diagnostic <- file.path(directory, "diagnostics", paste0(target, ".rds"))
    if (file.exists(diagnostic)) unlink(diagnostic, force = TRUE)
  }
  manifest$snapshots <- manifest$snapshots[!manifest$snapshots$id %in% targets, , drop = FALSE]
  manifest$updated <- format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")
  .nm_workspace_atomic_save(manifest, file.path(directory, "manifest.rds"))
  invisible(TRUE)
}

#' Delete a modelling project
#'
#' @param workspace Workspace object or path.
#' @param project Project id.
#' @return `TRUE`, invisibly.
#' @export
nm_project_delete <- function(workspace, project) {
  workspace <- if (inherits(workspace, "nm_workspace")) workspace else nm_workspace(workspace)
  project <- .nm_workspace_component(project)
  lock <- .nm_workspace_lock(workspace, paste0("project-", project))
  on.exit(.nm_workspace_unlock(lock), add = TRUE)
  directory <- .nm_project_path(workspace, project)
  if (!dir.exists(directory)) .nm_stop("Unknown project: ", project, ".")
  root <- normalizePath(file.path(workspace$path, "projects"), winslash = "/", mustWork = TRUE)
  target <- normalizePath(directory, winslash = "/", mustWork = TRUE)
  if (!startsWith(tolower(target), paste0(tolower(root), "/"))) {
    .nm_stop("Resolved project path escaped the workspace.")
  }
  status <- unlink(target, recursive = TRUE, force = TRUE)
  if (!identical(status, 0L) || dir.exists(target)) {
    .nm_stop("Unable to delete project: ", project, ".")
  }
  invisible(TRUE)
}

.nm_workspace_references <- function(workspace) {
  workspace <- if (inherits(workspace, "nm_workspace")) workspace else nm_workspace(workspace)
  projects <- list.dirs(file.path(workspace$path, "projects"), recursive = FALSE,
                        full.names = TRUE)
  references <- list()
  for (project in projects) {
    snapshots <- list.files(file.path(project, "snapshots"), pattern = "\\.rds$",
                            full.names = TRUE)
    for (path in snapshots) {
      value <- tryCatch(.nm_workspace_read(path), error = function(error) NULL)
      storage <- value$storage %||% NULL
      if (!is.list(storage) || !identical(storage$schema %||% "",
                                          "liberation.snapshot-storage")) next
      for (name in c("model", "data", "result")) {
        reference <- storage[[name]]
        if (!is.null(reference)) references[[length(references) + 1L]] <- reference
      }
    }
  }
  references
}

#' Verify a LibeRation workspace
#'
#' Checks manifests, snapshot contracts, content-addressed objects, and hashes
#' without changing the workspace.
#'
#' @param workspace Workspace object or path.
#' @return A `liber_workspace_verification` report.
#' @export
nm_workspace_verify <- function(workspace) {
  workspace <- if (inherits(workspace, "nm_workspace")) workspace else nm_workspace(workspace)
  issues <- list()
  add <- function(type, path, message) {
    issues[[length(issues) + 1L]] <<- data.frame(
      type = type, path = path, message = message, stringsAsFactors = FALSE
    )
  }
  projects <- list.dirs(file.path(workspace$path, "projects"), recursive = FALSE,
                        full.names = TRUE)
  checked_snapshots <- 0L
  for (directory in projects) {
    project <- basename(directory)
    manifest <- tryCatch(.nm_project_manifest(workspace, project), error = identity)
    if (inherits(manifest, "error")) {
      add("manifest", directory, conditionMessage(manifest)); next
    }
    for (id in manifest$snapshots$id) {
      path <- file.path(directory, "snapshots", paste0(id, ".rds"))
      value <- tryCatch(.nm_snapshot_unpack(workspace, .nm_workspace_read(path)), error = identity)
      checked_snapshots <- checked_snapshots + 1L
      if (inherits(value, "error")) add("snapshot", path, conditionMessage(value))
    }
  }
  table <- if (length(issues)) do.call(rbind, issues) else data.frame(
    type = character(), path = character(), message = character(), stringsAsFactors = FALSE
  )
  structure(list(
    valid = !nrow(table), workspace = workspace$path, projects = length(projects),
    snapshots = checked_snapshots, objects = length(.nm_workspace_references(workspace)),
    issues = table, checked = .nm_workspace_now()
  ), class = "liber_workspace_verification")
}

#' Migrate an existing workspace to the current schema
#'
#' Manifest upgrades are lossless. With `convert_snapshots = TRUE`, legacy
#' snapshot payloads are rewritten to content-addressed v2 storage after a
#' backup has been created.
#'
#' @param workspace Workspace object or path.
#' @param convert_snapshots Convert legacy embedded snapshots.
#' @param backup Create a complete backup before conversion.
#' @return Migration counts, invisibly.
#' @export
nm_workspace_migrate <- function(workspace, convert_snapshots = FALSE, backup = TRUE) {
  workspace <- if (inherits(workspace, "nm_workspace")) workspace else nm_workspace(workspace)
  if (isTRUE(backup)) nm_workspace_backup(workspace)
  directories <- list.dirs(file.path(workspace$path, "projects"), recursive = FALSE,
                           full.names = TRUE)
  converted <- 0L
  for (directory in directories) {
    project <- basename(directory)
    lock <- .nm_workspace_lock(workspace, paste0("project-", project))
    tryCatch({
      manifest_path <- file.path(directory, "manifest.rds")
      manifest <- .nm_project_manifest(workspace, project)
      .nm_workspace_atomic_save(manifest, manifest_path)
      if (isTRUE(convert_snapshots)) {
        for (id in manifest$snapshots$id) {
          path <- file.path(directory, "snapshots", paste0(id, ".rds"))
          value <- .nm_workspace_read(path)
          if (identical(as.integer(value$version %||% 1L), 1L)) {
            .nm_workspace_atomic_save(.nm_snapshot_pack(workspace, value), path)
            converted <- converted + 1L
          }
        }
      }
    }, finally = .nm_workspace_unlock(lock))
  }
  invisible(list(projects = length(directories), snapshots_converted = converted,
                 version = 2L))
}

#' Create a consistent workspace backup
#'
#' @param workspace Workspace object or path.
#' @param path Destination `.tar.gz`; defaults to the workspace backup folder.
#' @return Normalized backup path, invisibly.
#' @export
nm_workspace_backup <- function(workspace, path = NULL) {
  workspace <- if (inherits(workspace, "nm_workspace")) workspace else nm_workspace(workspace)
  if (is.null(path)) {
    path <- file.path(workspace$path, "backups", paste0(
      "liberation-workspace-", format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"), ".tar.gz"
    ))
  }
  path <- normalizePath(path.expand(path), winslash = "/", mustWork = FALSE)
  if (!dir.exists(dirname(path))) dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  lock <- .nm_workspace_lock(workspace, "workspace-backup", timeout = 30)
  on.exit(.nm_workspace_unlock(lock), add = TRUE)
  old <- setwd(workspace$path)
  on.exit(setwd(old), add = TRUE)
  include <- c("workspace.rds", "projects", "objects", "reports")
  include <- include[file.exists(include) | dir.exists(include)]
  utils::tar(path, files = include, compression = "gzip", tar = "internal")
  if (!file.exists(path)) .nm_stop("Workspace backup was not created.")
  invisible(path)
}

#' Remove unreferenced content-addressed workspace objects
#'
#' @param workspace Workspace object or path.
#' @param dry_run Report candidates without deleting them. This defaults to
#'   `TRUE` deliberately.
#' @return A data frame of unreferenced objects.
#' @export
nm_workspace_gc <- function(workspace, dry_run = TRUE) {
  workspace <- if (inherits(workspace, "nm_workspace")) workspace else nm_workspace(workspace)
  referenced <- unique(vapply(.nm_workspace_references(workspace), `[[`, character(1), "hash"))
  files <- list.files(file.path(workspace$path, "objects"), pattern = "^[a-f0-9]{64}\\.rds$",
                      recursive = TRUE, full.names = TRUE)
  hashes <- sub("\\.rds$", "", basename(files))
  candidates <- files[!hashes %in% referenced]
  report <- data.frame(hash = hashes[!hashes %in% referenced], path = candidates,
                       stringsAsFactors = FALSE)
  if (!isTRUE(dry_run) && nrow(report)) {
    lock <- .nm_workspace_lock(workspace, "workspace-gc", timeout = 30)
    on.exit(.nm_workspace_unlock(lock), add = TRUE)
    failed <- vapply(report$path, function(path) unlink(path, force = TRUE), integer(1)) != 0L
    if (any(failed)) .nm_stop("Unable to delete one or more unreferenced workspace objects.")
  }
  report
}

#' @export
print.nm_workspace <- function(x, ...) {
  projects <- nm_project_list(x)
  cat("LibeRation workspace\n  path:", x$path, "\n  projects:", nrow(projects), "\n")
  invisible(x)
}

#' @export
print.liber_workspace_verification <- function(x, ...) {
  cat("LibeRation workspace verification\n")
  cat("  valid:", x$valid, " projects:", x$projects, " snapshots:", x$snapshots, "\n")
  if (nrow(x$issues)) print(x$issues, row.names = FALSE)
  invisible(x)
}
