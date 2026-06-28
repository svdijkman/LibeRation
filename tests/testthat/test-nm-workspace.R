test_that("workspace sandbox rejects paths outside root", {
  root <- file.path(tempdir(), "LibeRation_ws_test", format(Sys.time(), "%s"))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  nm_workspace_init(root, create_demo_project = FALSE)
  expect_error(
    .nm_ws_resolve("..", "etc", "passwd", root = root),
    "outside workspace"
  )
})

test_that("nm_ctl_compose handles nm_model THETAs without bounds columns", {
  sim <- nm_synthetic_theo(n_sub = 3L, seed = 1L)
  expect_false("Lower" %in% names(sim$model$THETAS))
  ctl <- .nm_model_to_ctl(sim$model, data_file = "data/theo.csv", prob = "THEO")
  expect_true(grepl("\\$THETA", ctl))
  parts <- nm_ctl_parse(ctl)
  expect_equal(nrow(parts$thetas), 5L)
})

test_that("workspace new version on theo project after empty versions", {
  root <- file.path(tempdir(), "LibeRation_ws_newver", format(Sys.time(), "%s"))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  nm_workspace_init(root, create_demo_project = FALSE)
  nm_workspace_create_project("theo", path = root, template = "theo")
  nm_workspace_delete_version("theo", "theo_base", root = root)
  expect_length(nm_workspace_list_versions("theo", root = root), 0L)
  vid <- nm_workspace_new_version("theo", root = root)
  expect_true(nzchar(vid))
  expect_length(nm_workspace_list_versions("theo", root = root), 1L)
})

test_that("nm_ctl essential input depends on ADVAN/TRANS", {
  oral <- nm_ctl_essential_input_cols(4L, 4L)
  expect_true(all(c("ID", "TIME", "KA", "F1") %in% oral))
  iv <- nm_ctl_essential_input_cols(3L, 4L)
  expect_false("KA" %in% iv)
  expect_true(all(c("ID", "TIME", "F1", "S2") %in% iv))
  expect_true(nm_ctl_is_valid_pair(4L, 4L))
  expect_false(nm_ctl_is_valid_pair(3L, 2L))
  expect_equal(nm_ctl_trans_choices(1L), c("1", "2"))
  expect_equal(nm_ctl_default_obscmp(4L, 4L), 2L)
  expect_equal(nm_ctl_default_obscmp(3L, 4L), 1L)
})

test_that("nm_ctl_parse and compose roundtrip", {
  sim <- nm_synthetic_theo(n_sub = 3L, seed = 1L)
  ctl <- .nm_model_to_ctl(sim$model, data_file = "data/theo.csv", prob = "THEO")
  parts <- nm_ctl_parse(ctl)
  expect_equal(parts$pk, trimws(sim$model$PRED))
  parts$output_cols <- c("ID", "TIME", "DV", "IPRED", "WRES")
  expect_equal(parts$advan, 4L)
  expect_equal(parts$trans, 4L)
  expect_equal(nrow(parts$thetas), 5L)
  rebuilt <- nm_ctl_compose(parts)
  parts2 <- nm_ctl_parse(rebuilt)
  expect_equal(parts2$thetas$Value, parts$thetas$Value)
  expect_equal(parts2$output_cols, parts$output_cols)
  expect_true(grepl("\\$PK", rebuilt))
  expect_true(grepl("\\$OUTPUT", rebuilt))
  expect_false(grepl("\\$EST", rebuilt))
})

test_that("estimation runs are stored per model version", {
  root <- file.path(tempdir(), "LibeRation_ws_runs", format(Sys.time(), "%s"))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  nm_workspace_init(root, create_demo_project = FALSE)
  nm_workspace_create_project("p", path = root, template = "empty")
  mod001 <- nm_workspace_new_version("p", root = root, label = "v1")
  mod002 <- nm_workspace_new_version("p", root = root, label = "v2")
  fit1 <- list(method = "FO", objective = 1, theta = 1)
  fit2 <- list(method = "FOCE", objective = 2, theta = 2)
  nm_workspace_save_run("p", mod001, "est001", fit1, root = root)
  nm_workspace_save_run("p", mod002, "est001", fit2, root = root)
  expect_equal(nrow(nm_workspace_list_runs("p", mod001, root = root)), 1L)
  expect_equal(nrow(nm_workspace_list_runs("p", mod002, root = root)), 1L)
  loaded1 <- nm_workspace_load_run_fit("p", mod001, "est001", root = root)
  loaded2 <- nm_workspace_load_run_fit("p", mod002, "est001", root = root)
  expect_equal(loaded1$objective, 1)
  expect_equal(loaded2$objective, 2)
})

test_that("workspace delete project version run", {
  root <- file.path(tempdir(), "LibeRation_ws_del", format(Sys.time(), "%s"))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  nm_workspace_init(root, create_demo_project = FALSE)
  nm_workspace_create_project("p", path = root, template = "theo")
  vers <- nm_workspace_list_versions("p", root = root)
  expect_length(vers, 1L)
  nm_workspace_save_run(
    "p", vers[[1L]], "est001",
    list(method = "FO", objective = 1, theta = 1),
    root = root
  )
  expect_equal(nrow(nm_workspace_list_runs("p", vers[[1L]], root = root)), 1L)
  nm_workspace_delete_run("p", vers[[1L]], "est001", root = root)
  expect_equal(nrow(nm_workspace_list_runs("p", vers[[1L]], root = root)), 0L)
  nm_workspace_delete_version("p", vers[[1L]], root = root)
  expect_length(nm_workspace_list_versions("p", root = root), 0L)
  nm_workspace_delete_project("p", root = root)
  expect_false("p" %in% nm_workspace_list_projects(root))
})

test_that("THEO demo control stream matches nm_synthetic_theo model", {
  root <- file.path(tempdir(), "LibeRation_ws_theo", format(Sys.time(), "%s"))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  nm_workspace_init(root, create_demo_project = FALSE)
  nm_workspace_create_project("theo", path = root, template = "theo")
  vers <- nm_workspace_list_versions("theo", root = root)
  expect_length(vers, 1L)
  sim <- nm_synthetic_theo(n_sub = 10L, seed = 1L)
  parsed <- nm_workspace_parse_model("theo", vers[[1L]], root = root)
  expect_equal(parsed$model$ADVAN, sim$model$ADVAN)
  expect_equal(parsed$model$TRANS, sim$model$TRANS)
  expect_equal(parsed$model$OBSCMP, sim$model$OBSCMP)
  expect_equal(parsed$model$INPUT, sim$model$INPUT)
  expect_equal(parsed$model$THETAS$Value, sim$model$THETAS$Value)
  expect_equal(parsed$model$OMEGAS$Value, sim$model$OMEGAS$Value)
  expect_equal(parsed$model$SIGMAS$Value, sim$model$SIGMAS$Value)
})

test_that("workspace create project and model I/O", {
  root <- file.path(tempdir(), "LibeRation_ws_io", format(Sys.time(), "%s"))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  nm_workspace_init(root, create_demo_project = FALSE)
  nm_workspace_create_project("testproj", path = root, template = "theo")
  projs <- nm_workspace_list_projects(root)
  expect_true("testproj" %in% projs)
  runs <- nm_workspace_list_versions("testproj", root = root)
  expect_length(runs, 1L)
  txt <- nm_workspace_read_model("testproj", runs[[1L]], root = root)
  expect_true(grepl("\\$PROBLEM", txt))
  parsed <- nm_workspace_parse_model("testproj", runs[[1L]], root = root)
  expect_false(is.null(parsed$model))
  expect_false(is.null(parsed$data))
  pt <- nm_workspace_param_table(parsed$model)
  expect_gt(nrow(pt), 0L)
})

test_that("nm_report_pdf writes pdf and manifest", {
  skip_if_not_installed("ggplot2")
  sim <- nm_synthetic_theo(n_sub = 5L, seed = 1L)
  fit <- nm_est(sim$model, sim$data, method = "FO", control = list(maxit = 5L))
  root <- file.path(tempdir(), "LibeRation_ws_report", format(Sys.time(), "%s"))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  nm_workspace_init(root, create_demo_project = FALSE)
  pdf_path <- file.path(root, "projects", "p", "reports", "r.pdf")
  dir.create(dirname(pdf_path), recursive = TRUE)
  res <- nm_report_pdf(
    fit,
    pdf_path,
    project_meta = list(project = "p", run_id = "run001")
  )
  expect_true(file.exists(res$pdf))
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    expect_true(file.exists(res$manifest))
  }
})

test_that("nm_ctl_template builds standard oral 2-compartment model", {
  parts <- nm_ctl_template(4L, 4L, data_file = "data/theo.csv", problem = "Test")
  expect_equal(parts$advan, 4L)
  expect_equal(parts$trans, 4L)
  expect_equal(nrow(parts$thetas), 5L)
  expect_equal(nrow(parts$omegas), 3L)
  expect_equal(nrow(parts$thetas), length(parts$thetas$Label))
  ctl <- nm_ctl_compose(parts)
  expect_true(grepl("\\$PK", ctl))
  expect_true(grepl("ADVAN=4", ctl))
  expect_true(grepl("KA", ctl))
})

test_that("nm_ctl_template builds IV and oral 3-compartment models", {
  iv <- nm_ctl_template(11L, 4L, data_file = "data.csv", problem = "IV3")
  expect_equal(iv$advan, 11L)
  expect_equal(nrow(iv$thetas), 6L)
  expect_true(grepl("Q3", iv$pk))
  expect_false(grepl("KA", iv$pk))
  oral <- nm_ctl_template(12L, 4L, data_file = "data.csv", problem = "Oral3")
  expect_equal(oral$advan, 12L)
  expect_equal(nrow(oral$thetas), 7L)
  expect_true(grepl("KA", oral$pk))
  expect_true(grepl("Q3", oral$pk))
})

test_that("ADVAN 6 template uses ode_ncomp for $DES", {
  parts <- nm_ctl_template(6L, data_file = "data.csv", problem = "ODE", ode_ncomp = 4L)
  expect_equal(parts$advan, 6L)
  expect_true(grepl("DADT\\(4\\)", parts$des, fixed = FALSE))
  expect_false(nm_ctl_show_trans(6L))
  expect_equal(nm_ctl_effective_trans(6L, NULL), 1L)
})

test_that("template thetas always align Label column length", {
  th <- LibeRation:::.nm_ctl_template_thetas(c(1, 2, 3))
  expect_equal(nrow(th), 3L)
  expect_equal(length(th$Label), 3L)
  parts <- nm_ctl_template(2L, 2L, data_file = "data.csv", problem = "Oral")
  expect_equal(nrow(parts$thetas), 4L)
  expect_equal(length(parts$thetas$Label), 4L)
  expect_silent(nm_ctl_compose(parts))
})

test_that("nm_workspace_copy_version duplicates model version", {
  root <- file.path(tempdir(), "LibeRation_ws_copy", format(Sys.time(), "%s"))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  nm_workspace_init(root, create_demo_project = FALSE)
  nm_workspace_create_project("p", path = root, template = "theo")
  v1 <- nm_workspace_list_versions("p", root = root)[[1L]]
  v2 <- nm_workspace_copy_version("p", v1, root = root)
  ctl1 <- nm_workspace_read_model("p", v1, root = root)
  ctl2 <- nm_workspace_read_model("p", v2, root = root)
  expect_equal(nm_ctl_parse(ctl1)$thetas$Value, nm_ctl_parse(ctl2)$thetas$Value)
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    meta <- jsonlite::fromJSON(
      file.path(root, "projects", "p", "versions", v2, "meta.json"),
      simplifyVector = TRUE
    )
    expect_equal(meta$label, paste("Copy of", v1))
  }
})

test_that("nm_workspace_import_dataset copies file into project", {
  root <- file.path(tempdir(), "LibeRation_ws_import", format(Sys.time(), "%s"))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  nm_workspace_init(root, create_demo_project = FALSE)
  nm_workspace_create_project("imp", path = root, template = "empty")
  src <- tempfile(fileext = ".csv")
  writeLines("ID,TIME,DV", src)
  rel <- nm_workspace_import_dataset("imp", src, "study.csv", root = root)
  expect_equal(rel, "data/study.csv")
  expect_true(file.exists(file.path(root, "projects", "imp", rel)))
})

test_that("nm_workspace_create_project stores optional description", {
  root <- file.path(tempdir(), "LibeRation_ws_desc", format(Sys.time(), "%s"))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  nm_workspace_init(root, create_demo_project = FALSE)
  nm_workspace_create_project("p1", path = root, template = "empty", description = "Test run")
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    meta <- jsonlite::fromJSON(
      file.path(root, "projects", "p1", "project.json"),
      simplifyVector = TRUE
    )
    expect_equal(meta$description, "Test run")
  }
})

test_that("nm_ctl compose and parse round-trip includes DES for ODE models", {
  parts <- nm_ctl_template(6L, 4L, data_file = "data.csv", problem = "ODE test")
  expect_true(nzchar(parts$des))
  ctl <- nm_ctl_compose(parts)
  expect_true(grepl("\\$DES", ctl))
  parsed <- nm_ctl_parse(ctl)
  expect_true(nzchar(parsed$des))
  expect_equal(parsed$advan, 6L)
})

test_that("project names with spaces are sanitized on create", {
  root <- file.path(tempdir(), "LibeRation_ws_space", format(Sys.time(), "%s"))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  nm_workspace_init(root, create_demo_project = FALSE)
  meta <- nm_workspace_create_project("Test large", path = root, template = "empty")
  expect_equal(meta$name, "Test_large")
  expect_true(dir.exists(file.path(root, "projects", "Test_large")))
  tmp <- tempfile(fileext = ".csv")
  writeLines("ID,TIME,DV", tmp)
  rel <- nm_workspace_import_dataset("Test_large", tmp, "d.csv", root = root)
  parts <- nm_ctl_template(2L, 2L, data_file = rel, problem = "Test")
  nm_workspace_new_version("Test_large", root = root, template_ctl = nm_ctl_compose(parts), data_file = rel)
  expect_length(nm_workspace_list_versions("Test_large", root = root), 1L)
})

test_that("empty project new version omits null label in meta", {
  root <- file.path(tempdir(), "LibeRation_ws_empty", format(Sys.time(), "%s"))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  nm_workspace_init(root, create_demo_project = FALSE)
  nm_workspace_create_project("empty1", path = root, template = "empty")
  ver <- nm_workspace_new_version("empty1", root = root)
  meta_path <- file.path(root, "projects", "empty1", "versions", ver, "meta.json")
  expect_true(file.exists(meta_path))
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    meta <- jsonlite::fromJSON(meta_path, simplifyVector = TRUE)
    expect_false("label" %in% names(meta))
  }
})

test_that("report base graphics helpers render without ggplot2", {
  obs <- data.frame(TIME = 1:3, DV = c(1, 2, 3), WRES = c(-0.1, 0, 0.1), IPRED = c(1.1, 2, 2.9))
  ind <- data.frame(ID = 1L, TIME = 1:3, IPRED = c(1.1, 2, 2.9))
  pdf <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf)
  expect_silent(LibeRation:::.nm_report_plot_time(obs, ind))
  expect_silent(LibeRation:::.nm_report_plot_scatter(obs, "IPRED", "DV vs IPRED"))
  expect_silent(LibeRation:::.nm_report_plot_wres(obs))
  grDevices::dev.off()
  expect_true(file.size(pdf) > 0L)
})

test_that("workspace simulation save and list work", {
  skip_if_not_installed("data.table")
  root <- file.path(tempdir(), "LibeRation_ws_sim", format(Sys.time(), "%s"))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  nm_workspace_init(root, create_demo_project = FALSE)
  nm_workspace_create_project("p", path = root, template = "theo")
  ver <- nm_workspace_list_versions("p", root = root)[[1L]]
  parsed <- nm_workspace_parse_model("p", ver, root = root)
  sim_dat <- nm_task("sim", parsed$model, parsed$data, seed = 42L)
  sim_out <- structure(list(data = sim_dat), class = "nm_dataset")
  sim_id <- nm_workspace_new_sim_id("p", ver, root = root)
  nm_workspace_save_sim("p", ver, sim_id, sim_out, root = root, label = "test", seed = 42L)
  sims <- nm_workspace_list_sims("p", ver, root = root)
  expect_equal(nrow(sims), 1L)
  expect_equal(sims$sim_id[[1L]], sim_id)
  loaded <- nm_workspace_load_sim("p", ver, sim_id, root = root)
  expect_true(inherits(loaded, "nm_dataset"))
  nm_workspace_delete_sim("p", ver, sim_id, root = root)
  expect_equal(nrow(nm_workspace_list_sims("p", ver, root = root)), 0L)
})
