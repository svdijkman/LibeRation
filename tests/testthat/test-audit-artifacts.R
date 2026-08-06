test_that("audit artifacts are opt-in for simulation and estimation", {
  fixture <- estimation_fixture()

  ordinary <- nm_simulate(fixture$model, fixture$data)
  expect_null(attr(ordinary, "audit_artifacts", exact = TRUE))

  simulated <- nm_simulate(
    fixture$model, fixture$data, seed = 42, audit_artifacts = TRUE
  )
  bundle <- attr(simulated, "audit_artifacts", exact = TRUE)
  expect_identical(bundle$schema, "liberation.audit-artifacts")
  expect_identical(bundle$operation, "simulate")
  expect_true(all(c("model.lst", "model.ctl", "model.tab") %in%
                    vapply(bundle$files, `[[`, character(1), "name")))

  fit <- nm_est(
    fixture$model, fixture$data, method = "FOCEI", maxit = 1,
    eta_maxit = 5, audit_artifacts = TRUE
  )
  fit_bundle <- attr(fit, "audit_artifacts", exact = TRUE)
  expect_identical(fit_bundle$operation, "estimate")
  expect_true(all(c("model.lst", "model.ctl", "model.ext", "model.phi") %in%
                    vapply(fit_bundle$files, `[[`, character(1), "name")))

  covariance_fit <- fit
  attr(covariance_fit, "audit_artifacts") <- NULL
  parameter_names <- c("THETA1", "THETA2", "OMEGA1", "SIGMA1")
  covariance <- diag(seq_along(parameter_names))
  dimnames(covariance) <- list(parameter_names, parameter_names)
  correlation <- diag(length(parameter_names))
  dimnames(correlation) <- list(parameter_names, parameter_names)
  covariance_fit$covariance <- list(
    status = "completed", covariance = covariance, correlation = correlation
  )
  covariance_fit <- LibeRation:::.nm_attach_audit_artifacts(
    covariance_fit, fixture$model, fixture$data, "estimate"
  )
  covariance_files <- vapply(
    attr(covariance_fit, "audit_artifacts")$files,
    `[[`, character(1), "name"
  )
  expect_true(all(c("model.cov", "model.cor") %in% covariance_files))
})

test_that("workspace runs materialize and verify audit artifacts without duplication", {
  root <- tempfile("liber-audit-workspace-")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  workspace <- nm_workspace(root)
  project <- nm_project_create(workspace, "Audit artifacts")
  fixture <- estimation_fixture()
  version <- nm_project_save(
    workspace, project$id, fixture$model, fixture$data
  )
  result <- nm_simulate(
    fixture$model, fixture$data, seed = 17, audit_artifacts = TRUE
  )
  run <- nm_project_save_run(
    workspace, project$id, version, result, label = "Audited simulation"
  )

  record <- nm_project_list(workspace, project$id)
  record <- record[record$id == run, , drop = FALSE]
  expect_true(record$has_audit_artifacts)
  saved <- nm_project_load(workspace, project$id, run)
  expect_null(attr(saved$result, "audit_artifacts", exact = TRUE))
  expect_identical(saved$provenance$audit_artifacts$engine, "liber")

  artifacts <- nm_project_audit_artifacts(
    workspace, project$id, run, verify = TRUE
  )
  expect_true(dir.exists(artifacts$path))
  expect_true(all(file.exists(artifacts$files)))
  expect_true(file.exists(file.path(artifacts$path, "audit-manifest.json")))

  expect_true(nm_project_delete_snapshot(workspace, project$id, run))
  expect_false(dir.exists(artifacts$path))
})

test_that("audit artifact bundles reject unsafe filenames and damaged content", {
  fixture <- estimation_fixture()
  result <- nm_simulate(
    fixture$model, fixture$data, audit_artifacts = TRUE
  )
  bundle <- attr(result, "audit_artifacts", exact = TRUE)
  unsafe <- bundle
  unsafe$files[[1L]]$name <- "../model.lst"
  expect_error(
    LibeRation:::.nm_audit_bundle_validate(unsafe), "unsafe or duplicate"
  )
  damaged <- bundle
  damaged$files[[1L]]$content <- paste0(damaged$files[[1L]]$content, "changed")
  expect_error(
    LibeRation:::.nm_audit_bundle_validate(damaged), "SHA-256"
  )
})

test_that("workbench dialogs expose disabled-by-default audit controls", {
  script <- paste(
    readLines(system.file(
      "htmlwidgets", "liberWorkbench.js", package = "LibeRation"
    ), warn = FALSE),
    collapse = "\n"
  )
  expect_match(script, "estimationAuditArtifacts = React.useState\\(false\\)")
  expect_match(script, "simulationAuditArtifacts = React.useState\\(false\\)")
  expect_match(script, "auditArtifacts:!!estimationAuditArtifacts\\[0\\]")
  expect_match(script, "auditArtifacts:!!simulationAuditArtifacts\\[0\\]")
})
