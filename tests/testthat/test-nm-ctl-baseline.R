test_that("model file canonical baseline round-trips", {
  ws <- tempfile("nm_ws_baseline_")
  on.exit(unlink(ws, recursive = TRUE), add = TRUE)
  nm_workspace_init(ws, create_demo_project = FALSE)
  nm_workspace_create_project("demo", path = ws, template = "theo")
  vers <- nm_workspace_list_versions("demo", root = ws)
  expect_gt(length(vers), 0L)
  txt <- nm_workspace_read_model("demo", vers[[1L]], root = ws)
  baseline <- LibeRation:::.nm_ctl_baseline_text(txt)
  expect_true(nzchar(baseline))
  parts <- nm_ctl_parse(txt)
  recomposed <- nm_ctl_compose(parts)
  expect_equal(nm_ctl_canonical(recomposed), baseline)
  expect_equal(nm_ctl_canonical(txt), baseline)
})

test_that("reload baseline matches compose from parsed parts", {
  ws <- tempfile("nm_ws_baseline2_")
  on.exit(unlink(ws, recursive = TRUE), add = TRUE)
  nm_workspace_init(ws, create_demo_project = FALSE)
  nm_workspace_create_project("demo", path = ws, template = "empty")
  parts <- nm_ctl_template(2L, 1L, data_file = "data.csv", problem = "Baseline test")
  version_id <- nm_workspace_new_version(
    "demo",
    root = ws,
    template_ctl = nm_ctl_compose(parts),
    data_file = NULL,
    label = "v1"
  )
  txt <- nm_workspace_read_model("demo", version_id, root = ws)
  saved <- LibeRation:::.nm_ctl_baseline_text(txt)
  current <- nm_ctl_canonical(nm_ctl_compose(nm_ctl_parse(txt)))
  expect_equal(current, saved)
})
