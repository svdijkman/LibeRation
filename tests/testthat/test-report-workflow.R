test_that("report workflows persist and render PDF provenance bundles", {
  skip_if_not(exists(
    "_LibeRation_liberation_prediction_tape_new_dynamic",
    envir = asNamespace("LibeRation"), inherits = FALSE
  ), "The installed native library does not match the current R sources")
  fixture <- estimation_fixture()
  fit <- nm_est(fixture$model, fixture$data, method = "FO", maxit = 1)
  root <- tempfile("report-workflow-")
  workspace <- nm_workspace(root)
  project <- nm_project_create(workspace, "Report workflow")
  version <- nm_project_save(workspace, project$id, fixture$model, fixture$data,
                             label = "Base model")
  run <- nm_project_save_run(workspace, project$id, version, fit, "FO run")
  design <- nm_report_design(list(
    nm_report_block("introduction", text = "A reproducible report."),
    nm_report_block(
      "run", run_ids = run,
      elements = c("summary", "parameters", "gof", "run_info")
    ),
    nm_report_block("conclusion", source = "ai", text = "Reviewed draft.")
  ), formats = "pdf")

  expect_equal(nm_report_design_save(workspace, project$id, design), design$id)
  expect_s3_class(nm_report_design_load(workspace, project$id, design$id),
                  "nm_report_design")
  expect_equal(nrow(nm_report_design_load(workspace, project$id)), 1L)

  bundle <- nm_report_design_render(
    design, workspace, project$id, directory = root,
    name = "workflow-report", formats = "pdf"
  )
  expect_s3_class(bundle, "nm_report_bundle")
  expect_true(file.exists(bundle$pdf))
  expect_true(file.exists(bundle$json))
  expect_gt(file.info(bundle$pdf)$size, 1000)
})

test_that("report workflows generate DOCX when Pandoc is available", {
  skip_if_not(exists(
    "_LibeRation_liberation_prediction_tape_new_dynamic",
    envir = asNamespace("LibeRation"), inherits = FALSE
  ), "The installed native library does not match the current R sources")
  skip_if(!nzchar(LibeRation:::.nm_report_pandoc()), "Pandoc is unavailable")
  fixture <- estimation_fixture()
  fit <- nm_est(fixture$model, fixture$data, method = "FO", maxit = 1)
  root <- tempfile("report-docx-")
  workspace <- nm_workspace(root)
  project <- nm_project_create(workspace, "DOCX workflow")
  version <- nm_project_save(workspace, project$id, fixture$model, fixture$data)
  run <- nm_project_save_run(workspace, project$id, version, fit)
  design <- nm_report_design(list(
    nm_report_block("methods", text = "FO estimation."),
    nm_report_block("run", run_ids = run, elements = c("parameters", "gof"))
  ), formats = "docx")
  bundle <- nm_report_design_render(
    design, workspace, project$id, root, "workflow-docx", formats = "docx"
  )
  expect_true(file.exists(bundle$docx))
  expect_gt(file.info(bundle$docx)$size, 1000)
})

test_that("specialized outcome diagnostics are selectable report evidence", {
  elements <- c("vpc_count", "vpc_competing", "vpc_recurrent")
  block <- nm_report_block("run", run_ids = "run-1", elements = elements)
  expect_equal(block$elements, elements)

  fixtures <- list(
    vpc_count = structure(list(
      observed = data.frame(TIME = 0:1, MEAN = c(1, 2)),
      simulated = data.frame(
        TIME = 0:1, MEAN_lower = c(.5, 1), MEAN_median = c(1, 2),
        MEAN_upper = c(1.5, 3)
      )
    ), class = "nm_vpc_count"),
    vpc_competing = structure(list(
      observed = data.frame(
        TIME = rep(0:1, 2), CAUSE = rep(c("A", "B"), each = 2),
        CIF = c(0, .2, 0, .1)
      ),
      simulated = data.frame(
        TIME = rep(0:1, 2), CAUSE = rep(c("A", "B"), each = 2),
        lower = c(0, .1, 0, .05), median = c(0, .2, 0, .1),
        upper = c(0, .3, 0, .2)
      )
    ), class = "nm_vpc_competing"),
    vpc_recurrent = structure(list(
      observed = data.frame(TIME = 0:1, MEAN_CUMULATIVE = c(0, 1)),
      simulated = data.frame(
        TIME = 0:1, lower = c(0, .5), median = c(0, 1), upper = c(0, 1.5)
      )
    ), class = "nm_vpc_recurrent")
  )
  for (kind in names(fixtures)) {
    file <- tempfile(fileext = ".png")
    LibeRation:::.nm_report_plot_diagnostic(fixtures[[kind]], kind, file)
    expect_true(file.exists(file))
    expect_gt(file.info(file)$size, 1000)
  }
})
