test_that("workspaces preserve immutable serializable modelling snapshots", {
  root <- tempfile("liber-workspace-")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  workspace <- nm_workspace(root)
  expect_identical(workspace$version, 2L)
  project <- nm_project_create(
    workspace, "First PK project", description = "Population PK analysis"
  )
  expect_equal(project$id, "first-pk-project")
  expect_equal(nm_project_list(workspace)$description, "Population PK analysis")

  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT"), ADVAN = 1,
    PRED = "CL=THETA(1);V=THETA(2);S1=V", ERROR = "Y=F",
    THETAS = data.frame(THETA = 1:2, Value = c(2, 20))
  )
  data <- data.frame(
    ID = 1, TIME = c(0, 1), EVID = c(1, 0), AMT = c(100, 0)
  )
  first <- nm_project_save(workspace, project$id, nm_compile(model), data,
                           nm_simulate(model, data), label = "Baseline")
  Sys.sleep(0.01)
  second <- nm_project_save(workspace, project$id, model, data, label = "Revision")
  expect_false(identical(first, second))
  expect_equal(nm_project_list(workspace)$snapshots, 2)
  expect_equal(nrow(nm_project_list(workspace, project$id)), 2)

  snapshot <- nm_project_load(workspace, project$id)
  expect_equal(snapshot$id, second)
  expect_s3_class(snapshot$model, "nm_model")
  expect_false(inherits(snapshot$model, "NMEngine"))
  expect_s3_class(snapshot$data, "nm_dataset")
  expect_null(snapshot$result)
  expect_true(nm_workspace_verify(workspace)$valid)
  object_files <- list.files(file.path(root, "objects"), pattern = "\\.rds$", recursive = TRUE)
  expect_true(length(object_files) >= 2L)
  expect_error(nm_project_load(workspace, "../escape"), "Invalid project id")

  copied <- nm_project_copy(workspace, project$id, first, label = "Baseline copy")
  expect_equal(nm_project_load(workspace, project$id, copied)$label, "Baseline copy")
  expect_true(nm_project_delete_snapshot(workspace, project$id, copied))
  expect_false(copied %in% nm_project_list(workspace, project$id)$id)
})

test_that("workspace v2 deduplicates model and data objects", {
  root <- tempfile("liber-deduplicate-")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  workspace <- nm_workspace(root)
  project <- nm_project_create(workspace, "Deduplication")
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT"), ADVAN = 1,
    PRED = "CL=THETA(1);V=THETA(2);S1=V", ERROR = "Y=F",
    THETAS = data.frame(THETA = 1:2, Value = c(2, 20))
  )
  data <- data.frame(ID = 1, TIME = c(0, 1), EVID = c(1, 0), AMT = c(100, 0))
  first <- nm_project_save(workspace, project$id, model, data)
  second <- nm_project_save(workspace, project$id, model, data)
  first_raw <- readRDS(file.path(root, "projects", project$id, "snapshots", paste0(first, ".rds")))
  second_raw <- readRDS(file.path(root, "projects", project$id, "snapshots", paste0(second, ".rds")))
  expect_identical(first_raw$storage$model$hash, second_raw$storage$model$hash)
  expect_identical(first_raw$storage$data$hash, second_raw$storage$data$hash)
  expect_true(nm_workspace_verify(workspace)$valid)
  expect_equal(nrow(nm_workspace_gc(workspace)), 0L)
})

test_that("legacy workspace manifests and snapshots migrate explicitly", {
  root <- tempfile("liber-migrate-")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  workspace <- nm_workspace(root)
  project <- nm_project_create(workspace, "Legacy")
  model <- nm_model(
    INPUT = c("ID", "TIME"), ADVAN = 1, PRED = "CL=THETA(1);V=10;S1=V",
    THETAS = data.frame(THETA = 1, Value = 2)
  )
  id <- nm_project_save(workspace, project$id, model)
  snapshot_path <- file.path(root, "projects", project$id, "snapshots", paste0(id, ".rds"))
  value <- nm_project_load(workspace, project$id, id)
  value$version <- 1L; value$schema <- value$storage <- NULL
  saveRDS(value, snapshot_path)
  manifest_path <- file.path(root, "projects", project$id, "manifest.rds")
  manifest <- readRDS(manifest_path); manifest$version <- 1L; manifest$schema <- NULL
  saveRDS(manifest, manifest_path)
  migrated <- nm_workspace_migrate(workspace, convert_snapshots = TRUE, backup = FALSE)
  expect_identical(migrated$snapshots_converted, 1L)
  expect_identical(readRDS(manifest_path)$version, 2L)
  expect_identical(readRDS(snapshot_path)$version, 2L)
  expect_s3_class(nm_project_load(workspace, project$id, id)$model, "nm_model")
})

test_that("integrated package provenance is retained without replacing runtime provenance", {
  root <- tempfile("liber-provenance-")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  workspace <- nm_workspace(root)
  project <- nm_project_create(workspace, "Provenance")
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT"), ADVAN = 1,
    PRED = "CL=THETA(1);V=THETA(2);S1=V", ERROR = "Y=F",
    THETAS = data.frame(THETA = 1:2, Value = c(2, 20))
  )
  id <- nm_project_save(workspace, project$id, model,
                        provenance = list(LibeRary = list(library_id = "test-model", version = "1.0.0")))
  saved <- nm_project_load(workspace, project$id, id)
  expect_equal(saved$provenance$LibeRary$library_id, "test-model")
  expect_equal(saved$provenance$LibeRation, as.character(utils::packageVersion("LibeRation")))
  expect_match(saved$provenance$R, "R version")
})

test_that("projects can be removed without escaping their workspace", {
  root <- tempfile("liber-workspace-")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  workspace <- nm_workspace(root)
  project <- nm_project_create(workspace, "Disposable")
  expect_true(nm_project_delete(workspace, project$id))
  expect_false(project$id %in% nm_project_list(workspace)$id)
})

test_that("workspace payload exposes projects to the React workbench", {
  root <- tempfile("liber-workspace-")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  workspace <- nm_workspace(root)
  project <- nm_project_create(workspace, "GUI project")
  payload <- LibeRation:::.liber_gui_workspace(workspace, project$id)
  expect_true(payload$enabled)
  expect_equal(payload$current, project$id)
  expect_equal(payload$projects[[1]]$name, "GUI project")
})

test_that("model versions own numbered runs and persistent diagnostics", {
  root <- tempfile("liber-hierarchy-")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  workspace <- nm_workspace(root)
  project <- nm_project_create(workspace, "Hierarchy")
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT"), ADVAN = 1,
    PRED = "CL=THETA(1);V=THETA(2);S1=V", ERROR = "Y=F",
    THETAS = data.frame(THETA = 1:2, Value = c(2, 20))
  )
  data <- data.frame(ID = 1, TIME = c(0, 1), EVID = c(1, 0), AMT = c(100, 0))
  version <- nm_project_save(workspace, project$id, model, data)
  expect_equal(nm_project_load(workspace, project$id, version)$label, "Mod001")
  simulation <- nm_simulate(model, data)
  first <- nm_project_save_run(workspace, project$id, version, simulation)
  second <- nm_project_save_run(workspace, project$id, version, simulation)
  records <- nm_project_list(workspace, project$id)
  runs <- records[records$entry_type == "run", , drop = FALSE]
  expect_equal(sort(runs$run_number), 1:2)
  expect_true(all(runs$parent_id == version))

  # Use a classed sentinel here: persistence must not execute or transform a
  # diagnostic object while saving or loading it.
  diagnostic <- structure(list(nsim = 25L), class = "nm_vpc")
  # Diagnostics are only valid on estimation runs, so turn the test record into
  # a lightweight estimation-shaped run at the manifest boundary.
  manifest_path <- file.path(root, "projects", project$id, "manifest.rds")
  manifest <- readRDS(manifest_path)
  index <- match(first, manifest$snapshots$id)
  manifest$snapshots$result_type[[index]] <- "estimation"
  saveRDS(manifest, manifest_path)
  count <- structure(list(nsim = 25L), class = "nm_vpc_count")
  competing <- structure(list(nsim = 25L), class = "nm_vpc_competing")
  recurrent <- structure(list(nsim = 25L), class = "nm_vpc_recurrent")
  saved <- nm_project_save_diagnostics(workspace, project$id, first, list(
    vpc = diagnostic, vpc_count = count,
    vpc_competing = competing, vpc_recurrent = recurrent
  ))
  expect_identical(saved$vpc, diagnostic)
  expect_identical(nm_project_load_diagnostics(workspace, project$id, first)$vpc, diagnostic)
  expect_true(nm_project_list(workspace, project$id)$has_vpc[
    nm_project_list(workspace, project$id)$id == first
  ])
  record <- nm_project_list(workspace, project$id)
  record <- record[record$id == first, , drop = FALSE]
  expect_true(record$has_vpc_count)
  expect_true(record$has_vpc_competing)
  expect_true(record$has_vpc_recurrent)

  expect_true(nm_project_delete_snapshot(workspace, project$id, version))
  expect_equal(nrow(nm_project_list(workspace, project$id)), 0L)
})

test_that("GUI default workspace and empty projects are usable", {
  expected <- if (.Platform$OS.type == "windows") {
    file.path(Sys.getenv("USERPROFILE", unset = path.expand("~")),
              "Documents", "LibeR", "workspace")
  } else {
    file.path(path.expand("~"), "LibeR", "workspace")
  }
  expect_identical(LibeRation:::.liber_default_workspace(), expected)

  root <- tempfile("liber-empty-")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  workspace <- nm_workspace(root)
  project <- nm_project_create(workspace, "Empty project")
  payload <- LibeRation:::.liber_gui_workspace(workspace, project$id)
  expect_equal(payload$current, project$id)
  expect_length(payload$versions, 0L)
  expect_true(nm_project_delete(workspace, project$id))
})
