test_that("fit diagnostics expose predictions, residuals, ETAs, and shrinkage", {
  fixture <- estimation_fixture()
  fit <- nm_est(fixture$model, fixture$data, method = "FOCEI", maxit = 3)
  individual <- predict(fit)
  population <- predict(fit, type = "population")
  gof <- nm_gof(fit)
  etab <- nm_etab(fit)
  expect_equal(nrow(individual), nrow(fixture$data))
  expect_equal(nrow(population), nrow(fixture$data))
  expect_true(all(c("PRED", "IPRED", "RES", "IRES", "WRES", "IWRES") %in% names(gof)))
  expect_equal(residuals(fit), gof$IWRES)
  expect_equal(nrow(etab$eta), 3)
  expect_named(etab$shrinkage, "ETA1")
  expect_s3_class(summary(fit), "summary.nm_fit")
})

test_that("GOF population and individual predictions use independent ETA values", {
  fixture <- estimation_fixture()
  data <- nm_dataset(fixture$data)
  eta <- matrix(c(-0.35, 0.05, 0.4), ncol = 1L)
  fit <- structure(list(
    model = fixture$model, data = data,
    theta = fixture$model$THETAS$Value,
    omega = fixture$model$OMEGAS$Value,
    sigma = fixture$model$SIGMAS$Value,
    eta = eta
  ), class = "nm_fit")
  gof <- nm_gof(fit)
  population <- nm_simulate(fixture$model, data, eta = matrix(0, 3L, 1L))
  individual <- nm_simulate(fixture$model, data, eta = eta)
  expect_equal(gof$PRED, population$IPRED)
  expect_equal(gof$IPRED, individual$IPRED)
  expect_gt(max(abs(gof$IPRED - gof$PRED)), 0.1)
})

test_that("OPG covariance uses exact subject score tapes", {
  fixture <- estimation_fixture()
  fixture$model$THETAS$FIX <- c(FALSE, TRUE)
  fit <- nm_est(
    fixture$model, fixture$data, method = "FOCEI",
    maxit = 35, eta_maxit = 60, tolerance = 1e-6,
    covariance = TRUE, covariance_type = "opg"
  )
  covariance <- fit$covariance
  expect_s3_class(covariance, "nm_covariance")
  expect_equal(covariance$status, "completed")
  expect_equal(dim(covariance$covariance), c(1, 1))
  expect_named(covariance$se, "THETA1")
  expect_true(is.finite(covariance$se[[1]]) && covariance$se[[1]] > 0)
  expect_true(is.finite(fit$timing$covariance_seconds))
  expect_gte(fit$timing$total_seconds, fit$timing$model_fit_seconds)
  expect_gte(fit$timing$total_seconds, fit$timing$covariance_seconds)
  payload <- LibeRation:::.liber_gui_fit(fit)
  expect_true(payload$covariance$requested)
  expect_equal(payload$covariance$status, "completed")
  expect_true(nzchar(payload$run_info[["Model fit time"]]))
  expect_true(nzchar(payload$run_info[["Covariance step time"]]))
  expect_true(nzchar(payload$run_info[["Total estimation time"]]))
  expect_true(payload$run_info[["Iterations to convergence"]] >= 0L)
  expect_true(is.finite(payload$parameters[[1]]$se))
})

test_that("covariance results use the native variance-parameter scale", {
  fixture <- estimation_fixture()
  fixture$model$THETAS$FIX <- TRUE
  fixture$model$OMEGAS$FIX <- TRUE
  fixture$model$SIGMAS$FIX <- FALSE
  fit <- nm_est(
    fixture$model, fixture$data, method = "FO", maxit = 15,
    covariance = TRUE, covariance_type = "hessian"
  )
  expect_equal(fit$covariance$status, "completed")
  expect_equal(
    fit$covariance$objective_backend,
    "persistent-cpp-population-objective"
  )
  expect_gte(fit$covariance$objective_telemetry$parameter_evaluations, 1L)
  expect_named(fit$covariance$se, "SIGMA1")
  expect_true(is.finite(fit$covariance$se[[1]]) && fit$covariance$se[[1]] > 0)
  fo_opg <- nm_est(
    fixture$model, fixture$data, method = "FO",
    covariance = TRUE, covariance_type = "opg"
  )
  expect_equal(fo_opg$covariance$status, "completed")
  expect_equal(fo_opg$covariance$type, "opg")
  expect_true(all(is.finite(fo_opg$covariance$scores)))
})

test_that("automatic and sandwich covariance expose R and S diagnostics", {
  fixture <- estimation_fixture()
  fixture$model$THETAS$FIX <- c(FALSE, TRUE)
  fit <- nm_est(fixture$model, fixture$data, method = "FOCEI", maxit = 8)
  automatic <- nm_cov_step(fit, type = "auto")
  sandwich <- nm_cov_step(fit, type = "sandwich")
  expect_equal(automatic$status, "completed")
  expect_true(automatic$type %in% c("hessian", "sandwich", "opg"))
  expect_equal(sandwich$type, "sandwich")
  expect_equal(dim(sandwich$bread), c(1L, 1L))
  expect_equal(dim(sandwich$meat), c(1L, 1L))
  expect_true(is.finite(sandwich$se[[1L]]))
})

test_that("IMP and SAEM covariance use marginal importance information", {
  fixture <- estimation_fixture()
  fixture$model$THETAS$FIX <- c(FALSE, TRUE)

  imp <- nm_est(
    fixture$model, fixture$data, method = "IMP",
    maxit = 8, eta_maxit = 50, n_imp = 20, seed = 41,
    covariance = TRUE, covariance_type = "hessian",
    covariance_samples = 20, covariance_seed = 41
  )
  expect_equal(imp$covariance$status, "completed")
  expect_equal(imp$covariance$samples, 20L)
  expect_equal(imp$covariance$seed, 41L)
  expect_named(imp$covariance$se, "THETA1")
  expect_true(is.finite(imp$covariance$se[[1]]) && imp$covariance$se[[1]] > 0)
  expect_equal(
    imp$covariance$objective_backend, "fixed-proposal-importance-score"
  )
  expect_equal(imp$covariance$objective_telemetry$proposals, 3L)
  expect_equal(imp$covariance$sampling, "tensor-gauss-hermite")
  expect_equal(
    imp$covariance$actual_samples,
    imp$covariance$quadrature_order^fixture$model$n_eta
  )
  expect_gte(imp$covariance$objective_telemetry$parameter_evaluations, 1L)

  saem <- nm_est(
    fixture$model, fixture$data, method = "SAEM",
    maxit = 3, eta_maxit = 50, n_iter = 8, burn = 2,
    mcmc_steps = 1, mstep_maxit = 2, seed = 43,
    covariance = TRUE, covariance_type = "opg",
    covariance_samples = 20, covariance_seed = 43
  )
  expect_equal(saem$covariance$status, "completed")
  expect_named(saem$covariance$rse, "THETA1")
  expect_true(is.finite(saem$covariance$rse[[1]]))
  expect_equal(
    saem$covariance$objective_backend, "fixed-proposal-importance-score"
  )
})

test_that("fixed-proposal importance gradients match their marginal objective", {
  fixture <- estimation_fixture()
  fixture$model$THETAS$FIX <- c(FALSE, TRUE)
  context <- LibeRation:::.nm_estimation_context(fixture$model, fixture$data)
  map <- LibeRation:::.nm_outer_map(context$model)
  normals <- LibeRation:::.nm_imp_normals(context, 30L, 91L)
  objective <- LibeRation:::.nm_imp_information_objective(
    context, map, normals, map$start, eta_maxit = 60L, tolerance = 1e-7
  )
  analytic <- attr(objective, "gradient")(map$start)
  numerical <- LibeRation:::.nm_numeric_gradient(
    objective, map$start, relative_step = 1e-5
  )
  expect_equal(analytic, numerical, tolerance = 2e-4)
  scores <- attr(objective, "subject_scores")(map$start)
  expect_equal(colSums(scores), -0.5 * analytic, tolerance = 1e-10)
})

test_that("marginal covariance designs use bounded quadrature and random fallback", {
  low <- LibeRation:::.nm_imp_covariance_design(
    list(n_eta = 3L, n_subjects = 2L), 200L, 17L
  )
  expect_equal(low$method, "tensor-gauss-hermite")
  expect_equal(low$actual_samples, 216L)
  expect_equal(dim(low$normals[[1L]]), c(216L, 3L))
  expect_equal(sum(exp(attr(low$normals[[1L]], "log_measure"))), 1,
               tolerance = 1e-12)

  high <- LibeRation:::.nm_imp_covariance_design(
    list(n_eta = 7L, n_subjects = 2L), 20L, 17L
  )
  expect_equal(high$method, "random-normal")
  expect_equal(high$actual_samples, 20L)
  expect_equal(dim(high$normals[[1L]]), c(20L, 7L))
})

test_that("Bayesian posterior uncertainty is distinct from covariance", {
  fixture <- estimation_fixture()
  fixture$model$THETAS$FIX <- c(FALSE, TRUE)
  fit <- nm_est(
    fixture$model, fixture$data, method = "BAYES",
    n_burn = 5, n_sample = 20, n_thin = 1, seed = 47
  )
  population <- fit$posterior$population
  expect_named(population$sd, c("THETA1", "THETA2", "SIGMA1", "OMEGA1"))
  expect_equal(dim(population$quantile), c(3L, 4L))
  expect_equal(dim(population$covariance), c(4L, 4L))
  expect_error(nm_cov_step(fit), "posterior uncertainty")

  payload <- LibeRation:::.liber_gui_fit(fit)
  expect_true(payload$posterior$available)
  expect_equal(payload$posterior$samples, 20L)
  expect_true(is.finite(payload$parameters[[1]]$posterior_sd))
  expect_true(is.finite(payload$parameters[[1]]$lower_95))
})

test_that("VPC and subject bootstrap produce versionable workflow results", {
  fixture <- estimation_fixture()
  fit <- nm_est(fixture$model, fixture$data, method = "FOCEI", maxit = 2)
  vpc <- nm_vpc(fit, nsim = 5, seed = 17)
  pc_vpc <- nm_vpc(fit, nsim = 5, seed = 17, pc_correct = TRUE)
  expect_s3_class(vpc, "nm_vpc")
  expect_s3_class(pc_vpc, "nm_vpc")
  expect_true(pc_vpc$pc_correct)
  expect_true(LibeRation:::.liber_gui_result(pc_vpc)$pc_correct)
  expect_true(nrow(vpc$observed) > 0)
  expect_true(nrow(vpc$simulated) > 0)
  fit$data$SEX <- rep(c("F", "M", "M"), each = 4L)
  stratified <- nm_vpc(fit, nsim = 5, seed = 18, stratify = "SEX")
  expect_equal(stratified$stratify, "SEX")
  expect_setequal(vapply(stratified$stratified, `[[`, character(1), "level"), c("F", "M"))
  payload <- LibeRation:::.liber_gui_result(stratified)
  expect_equal(length(payload$stratified), 2L)
  expect_true(all(vapply(payload$stratified, function(item) length(item$points) > 0L, logical(1))))
  bootstrap <- nm_bootstrap(
    fit, n = 2, seed = 17, maxit = 2, eta_maxit = 30
  )
  expect_s3_class(bootstrap, "nm_bootstrap")
  expect_equal(bootstrap$n, 2)
  expect_equal(bootstrap$successful + length(bootstrap$errors), 2)
})

test_that("NPC, NPDE, and reports support the legacy diagnostic workflow", {
  fixture <- estimation_fixture()
  fit <- nm_est(fixture$model, fixture$data, method = "FOCEI", maxit = 2)
  npc <- nm_npc(fit, nsim = 20, seed = 19)
  npde <- nm_npde(fit, nsim = 20, seed = 19)
  expect_s3_class(npc, "nm_npc")
  expect_s3_class(npde, "nm_npde")
  expect_true(all(npc$table$PERCENTILE > 0 & npc$table$PERCENTILE < 1))
  expect_true(all(is.finite(npde$table$NPDE)))

  path <- tempfile(fileext = ".pdf")
  on.exit(unlink(c(path, sub("[.]pdf$", ".json", path))), add = TRUE)
  report <- nm_report(fit, path, sections = c("summary", "parameters", "gof", "eta"))
  expect_s3_class(report, "nm_report")
  expect_true(file.exists(report$pdf))
  expect_true(file.exists(report$json))
})
