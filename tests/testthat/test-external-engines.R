test_that("execution-engine discovery is explicit and stable", {
  engines <- nm_execution_engines()
  expect_identical(engines$id, c("liber", "nonmem", "nlmixr2"))
  expect_identical(engines$label, c("LibeR", "NONMEM", "nlmixr2"))
  expect_true(engines$available[[1L]])
  expect_true(all(c("version", "location") %in% names(engines)))
})

test_that("NONMEM export can suppress estimation for simulation", {
  fixture <- estimation_fixture()
  control <- nm_control_write(
    fixture$model, data = "data.csv IGNORE=@",
    estimation = FALSE, covariance = FALSE
  )
  expect_false(grepl("$ESTIMATION", control, fixed = TRUE))
  expect_false(grepl("$COVARIANCE", control, fixed = TRUE))
})

test_that("nlmixr2 translation covers solved ADVAN and rejects unsupported errors", {
  fixture <- estimation_fixture()
  source <- LibeRation:::.nm_nlmixr_ui_text(fixture$model)
  expect_match(source, "ini\\(\\{")
  expect_match(source, "theta1 <-")
  expect_match(source, "eta1 ~")
  expect_match(source, "linCmt\\(\\) ~ add\\(sigma1\\)")
  expect_false(grepl("system\\(|library\\(|source\\(", source))

  model <- nm_model_update(fixture$model, ERROR = "Y=runif(1)")
  expect_error(
    LibeRation:::.nm_nlmixr_ui_text(model),
    "could not translate this \\$ERROR"
  )
})

test_that("nlmixr2 translation rejects model families it cannot preserve", {
  fixture <- estimation_fixture()
  model <- fixture$model
  model$HMM_CONFIG <- list(states = c("A", "B"))
  expect_error(
    LibeRation:::.nm_nlmixr_ui_text(model),
    "semantics-preserving translation.*hidden/semi-Markov"
  )
})

test_that("NONMEM output parser normalizes estimates and ETAs", {
  fixture <- estimation_fixture()
  directory <- tempfile("nonmem-parser-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE, force = TRUE), add = TRUE)
  writeLines(c(
    "ITERATION THETA1 THETA2 SIGMA(1,1) OMEGA(1,1) OBJ",
    "0 2.1 19.8 0.21 0.08 123.5",
    "-1000000000 2.2 19.7 0.22 0.07 122.5"
  ), file.path(directory, "model.ext"))
  writeLines(c(
    "SUBJECT_NO ID ETA(1)", "1 1 -0.1", "2 2 0.0", "3 3 0.2"
  ), file.path(directory, "model.phi"))
  writeLines(c(
    "ID TIME DV PRED IPRED CWRES", "1 1 4.5 4.4 4.5 0.1"
  ), file.path(directory, "model.tab"))
  writeLines(c(
    "MINIMIZATION SUCCESSFUL",
    "NO. OF FUNCTION EVALUATIONS USED: 12",
    "Elapsed estimation time in seconds: 0.45",
    "Elapsed covariance time in seconds: 0.05"
  ), file.path(directory, "model.lst"))

  fit <- LibeRation:::.nm_nonmem_fit(
    directory, fixture$model, nm_dataset(fixture$data), "FOCEI", 0.6
  )
  expect_s3_class(fit, "nm_fit")
  expect_identical(fit$execution_engine, "nonmem")
  expect_equal(fit$theta, c(2.2, 19.7))
  expect_equal(fit$omega, 0.07)
  expect_equal(fit$sigma, 0.22)
  expect_equal(unname(fit$eta[, 1L]), c(-0.1, 0, 0.2))
  expect_equal(fit$objective, 122.5)
  expect_identical(fit$convergence, 0L)
  expect_equal(fit$timing$total_seconds, 0.6)
})

test_that("external engine selection is present in both workbench dialogs", {
  script <- paste(
    readLines(system.file(
      "htmlwidgets", "liberWorkbench.js", package = "LibeRation"
    ), warn = FALSE), collapse = "\n"
  )
  expect_match(script, "estimationEngine = React.useState\\(\"liber\"\\)")
  expect_match(script, "simulationEngine = React.useState\\(\"liber\"\\)")
  expect_match(script, "engine:estimationEngine\\[0\\]")
  expect_match(script, "engine:simulationEngine\\[0\\]")
  expect_match(script, "value:\"nonmem\"")
  expect_match(script, "value:\"nlmixr2\"")
})

test_that("a live NONMEM smoke test is opt-in", {
  skip_if_not(identical(Sys.getenv("LIBERATION_TEST_NONMEM"), "true"))
  skip_if_not(nm_execution_engines()$available[nm_execution_engines()$id == "nonmem"])
  fixture <- estimation_fixture(fix = FALSE)
  fit <- nm_external_run(
    fixture$model, fixture$data, engine = "nonmem", operation = "estimate",
    method = "FO", maxit = 1L, covariance = FALSE,
    audit_artifacts = TRUE
  )
  expect_s3_class(fit, "nm_fit")
  expect_identical(fit$execution_engine, "nonmem")
  expect_true(any(vapply(
    attr(fit, "audit_artifacts")$files, `[[`, character(1), "name"
  ) == "model.lst"))
  simulated <- nm_external_run(
    fixture$model, fixture$data, engine = "nonmem", operation = "simulate",
    nsim = 2L, seed = 20260806L, audit_artifacts = TRUE
  )
  expect_s3_class(simulated, "nm_dataset")
  expect_equal(nrow(simulated), 2L * nrow(fixture$data))
  expect_true(all(c("PRED", "IPRED") %in% names(simulated)))
})

test_that("a live nlmixr2 simulation smoke test is opt-in", {
  skip_if_not(identical(Sys.getenv("LIBERATION_TEST_NLMIXR2"), "true"))
  skip_if_not(nm_execution_engines()$available[nm_execution_engines()$id == "nlmixr2"])
  fixture <- estimation_fixture()
  simulated <- nm_external_run(
    fixture$model, fixture$data, engine = "nlmixr2", operation = "simulate",
    nsim = 1L, seed = 20260806L, random_effects = TRUE,
    residual = FALSE, audit_artifacts = TRUE
  )
  expect_s3_class(simulated, "nm_dataset")
  expect_identical(attr(simulated, "solver"), "nlmixr2/rxode2")
  expect_true(nrow(simulated) > 0L)
  expect_true(any(vapply(
    attr(simulated, "audit_artifacts")$files, `[[`, character(1), "name"
  ) == "model.R"))

  fit_fixture <- estimation_fixture(fix = FALSE)
  fit <- nm_external_run(
    fit_fixture$model, fit_fixture$data,
    engine = "nlmixr2", operation = "estimate", method = "FOCEI",
    maxit = 2L, eta_maxit = 5L, covariance = FALSE, n_cores = 1L
  )
  expect_s3_class(fit, "nm_fit")
  expect_identical(fit$execution_engine, "nlmixr2")
  expect_length(fit$theta, nrow(fit_fixture$model$THETAS))
  expect_equal(nrow(fit$output), nrow(fit_fixture$data))
})
