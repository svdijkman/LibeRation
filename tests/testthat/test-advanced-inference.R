test_that("MU references generate, validate, and round-trip semantic code", {
  fixture <- estimation_fixture()
  model <- nm_model(
    INPUT = fixture$model$INPUT, ADVAN = 1,
    PRED = "V=THETA(2); S1=V", ERROR = "Y=F+ERR(1)",
    THETAS = fixture$model$THETAS, OMEGAS = fixture$model$OMEGAS,
    SIGMAS = fixture$model$SIGMAS,
    MU = nm_mu(1, "log(THETA(1))", parameter = "CL")
  )
  expect_s3_class(model$MU, "nm_mu")
  expect_match(model$PRED, "MU_1=log\\(THETA\\(1\\)\\)")
  expect_match(model$PRED, "CL=exp\\(MU_1\\+ETA\\(1\\)\\)")

  contract <- nm_model_to_contract(model)
  expect_equal(contract$version, 4L)
  rebuilt <- nm_model_from_contract(contract)
  expect_equal(rebuilt$MU, model$MU)

  control <- nm_control_write(model)
  expect_match(control, "MU_1")
  imported <- nm_control_read(control)$model
  expect_equal(imported$MU$MU, 1L)
  expect_equal(imported$MU$ETA, 1L)

  eta <- matrix(c(-0.2, 0.1, 0.3), ncol = 1L)
  conventional <- nm_simulate(
    fixture$model, fixture$data, eta = eta,
    random_effects = FALSE, residual = FALSE
  )
  referenced <- nm_simulate(
    model, fixture$data, eta = eta,
    random_effects = FALSE, residual = FALSE
  )
  expect_equal(referenced$IPRED, conventional$IPRED, tolerance = 1e-12)
  expect_equal(nm_mu(2, "log(THETA(2))", "V")$MU, 2L)
})

test_that("MU covariates are required to be subject-level", {
  fixture <- estimation_fixture()
  model <- nm_model(
    INPUT = c(fixture$model$INPUT, "WT"), ADVAN = 1,
    PRED = "V=THETA(2); S1=V", ERROR = "Y=F+ERR(1)",
    THETAS = fixture$model$THETAS, OMEGAS = fixture$model$OMEGAS,
    SIGMAS = fixture$model$SIGMAS, COVARIATES = "WT",
    MU = nm_mu(1, "log(THETA(1))+0.75*log(WT/70)", "CL")
  )
  expect_equal(model$MU$COVARIATES, "WT")
  data <- fixture$data
  data$WT <- rep(c(60, 61, 60, 60), 3)
  expect_error(
    nm_est(model, data, method = "FOCEI", maxit = 1),
    "varies within subject"
  )
  subject_level <- data
  subject_level$WT <- rep(c(60, 70, 80), each = 4L)
  covariate_theta <- rbind(
    fixture$model$THETAS,
    data.frame(
      THETA = 3L, Value = 0.75, FIX = FALSE,
      LOWER = -1, UPPER = 2
    )
  )
  covariate_model <- nm_model(
    INPUT = c(fixture$model$INPUT, "WT"), ADVAN = 1,
    PRED = "V=THETA(2); S1=V", ERROR = "Y=F+ERR(1)",
    THETAS = covariate_theta, OMEGAS = fixture$model$OMEGAS,
    SIGMAS = fixture$model$SIGMAS, COVARIATES = "WT",
    MU = nm_mu(
      1, "log(THETA(1))+THETA(3)*log(WT/70)", "CL",
      covariates = "WT"
    )
  )
  specialization <- .nm_mu_specialization(
    .nm_estimation_context(covariate_model, subject_level, method = "IMP"),
    .nm_outer_map(covariate_model)
  )
  expect_true(specialization$covariate_design)
  expect_true(.nm_mu_diagnostic(specialization)$covariate_design)
})

test_that("MU specialization classifies, re-centres, and improves the Gaussian block", {
  fixture <- estimation_fixture(fix = FALSE)
  model <- nm_model(
    INPUT = fixture$model$INPUT, ADVAN = 1,
    PRED = "V=THETA(2); S1=V", ERROR = "Y=F+ERR(1)",
    THETAS = fixture$model$THETAS, OMEGAS = fixture$model$OMEGAS,
    SIGMAS = fixture$model$SIGMAS,
    MU = nm_mu(1, "log(THETA(1))", parameter = "CL")
  )
  context <- .nm_estimation_context(model, fixture$data, method = "SAEM")
  map <- .nm_outer_map(model)
  specialization <- .nm_mu_specialization(context, map)
  expect_true(specialization$active)
  expect_true(specialization$saem_eligible)
  expect_equal(specialization$theta, 1L)
  expect_equal(unname(specialization$links), "log")

  parameters <- map$decode(map$start)
  eta <- matrix(c(-0.4, 0.1, 0.35), ncol = 1L)
  updated <- .nm_mu_gls_update(
    specialization, context, parameters, eta
  )
  expect_true(updated$valid)
  old_phi <- eta + .nm_mu_values(specialization, parameters$theta)
  new_phi <- updated$eta +
    .nm_mu_values(specialization, updated$parameters$theta)
  expect_equal(new_phi, old_phi, tolerance = 1e-12)
  old_q <- sum(eta^2 / parameters$omega[[1L]])
  new_q <- sum(updated$eta^2 / parameters$omega[[1L]])
  expect_lte(new_q, old_q + 1e-12)
  repeated <- .nm_mu_gls_update(
    specialization, context, parameters, eta
  )
  expect_true(repeated$valid)
  cache_diagnostic <- .nm_mu_diagnostic(specialization)
  expect_true(cache_diagnostic$gls_vectorized)
  expect_equal(cache_diagnostic$gls_system_calls, 2L)
  expect_equal(cache_diagnostic$gls_cache_misses, 1L)
  expect_equal(cache_diagnostic$gls_cache_hits, 1L)

  nonlinear <- nm_model(
    INPUT = fixture$model$INPUT, ADVAN = 1,
    PRED = "V=THETA(2); S1=V", ERROR = "Y=F+ERR(1)",
    THETAS = fixture$model$THETAS, OMEGAS = fixture$model$OMEGAS,
    SIGMAS = fixture$model$SIGMAS,
    MU = nm_mu(1, "log(THETA(1))^2", parameter = "CL")
  )
  nonlinear_plan <- .nm_mu_plan(nonlinear, .nm_outer_map(nonlinear))
  expect_false(nonlinear_plan$affine)
  expect_false(nonlinear_plan$saem_eligible)
  expect_match(nonlinear_plan$reason, "not affine")

  aliased_theta <- fixture$model$THETAS
  aliased_theta <- rbind(
    aliased_theta,
    transform(aliased_theta[1L, , drop = FALSE],
              THETA = 3L, Value = 0.5, FIX = FALSE)
  )
  aliased <- nm_model(
    INPUT = fixture$model$INPUT, ADVAN = 1,
    PRED = "V=THETA(2); S1=V", ERROR = "Y=F+ERR(1)",
    THETAS = aliased_theta, OMEGAS = fixture$model$OMEGAS,
    SIGMAS = fixture$model$SIGMAS,
    MU = nm_mu(1, "THETA(1)+THETA(3)", parameter = "CL")
  )
  aliased_context <- .nm_estimation_context(
    aliased, fixture$data, method = "SAEM"
  )
  aliased_specialization <- .nm_mu_specialization(
    aliased_context, .nm_outer_map(aliased)
  )
  expect_false(aliased_specialization$active)
  expect_false(aliased_specialization$saem_eligible)
  expect_match(aliased_specialization$reason, "rank deficient")
})

test_that("SAEM, IMP, and BAYES report their MU-aware estimator paths", {
  fixture <- estimation_fixture(fix = FALSE)
  theta <- fixture$model$THETAS
  theta$FIX <- c(FALSE, TRUE)
  model <- nm_model(
    INPUT = fixture$model$INPUT, ADVAN = 1,
    PRED = "V=THETA(2); S1=V", ERROR = "Y=F+ERR(1)",
    THETAS = theta, OMEGAS = fixture$model$OMEGAS,
    SIGMAS = fixture$model$SIGMAS,
    MU = nm_mu(1, "log(THETA(1))", parameter = "CL")
  )
  saem <- nm_est(
    model, fixture$data, method = "SAEM", maxit = 2L,
    n_iter = 3L, burn = 1L, mcmc_steps = 1L, mstep_maxit = 1L,
    seed = 31L
  )
  expect_true(saem$diagnostics$mu_specialization$active)
  expect_equal(saem$diagnostics$mu_specialization$closed_form_updates, 3L)
  expect_equal(
    saem$diagnostics$mu_specialization$closed_form_only_iterations, 3L
  )
  expect_true(saem$diagnostics$mu_specialization$gls_vectorized)
  expect_equal(saem$diagnostics$mu_specialization$runtime_fallbacks, 0L)
  expect_match(
    .liber_gui_fit(saem, include_gof = FALSE)$run_info[[
      "MU estimator specialization"
    ]],
    "GLS M-step"
  )
  reference_saem <- nm_est(
    model, fixture$data, method = "SAEM", maxit = 2L,
    n_iter = 3L, burn = 1L, mcmc_steps = 1L, mstep_maxit = 1L,
    seed = 31L, mu_specialization = FALSE
  )
  expect_lt(saem$objective_evaluations, reference_saem$objective_evaluations)

  imp <- nm_est(
    model, fixture$data, method = "IMP", maxit = 2L,
    eta_maxit = 15L, n_imp = 5L, seed = 32L
  )
  expect_true(imp$diagnostics$mu_specialization$active)
  expect_gt(imp$diagnostics$mu_specialization$recentered_mode_starts, 0L)

  bayes <- nm_est(
    model, fixture$data, method = "BAYES",
    n_burn = 2L, n_sample = 3L, n_thin = 1L, seed = 33L
  )
  expect_true(bayes$diagnostics$mu_specialization$active)
  expect_equal(bayes$diagnostics$mu_specialization$attempted_blocks, 5L)
  expect_true(is.finite(bayes$diagnostics$mu_acceptance))
})

test_that("model comparisons retain evidence and explicit nested inference", {
  fixture <- estimation_fixture()
  fit <- nm_est(fixture$model, fixture$data, method = "FOCEI", maxit = 1)
  candidate <- fit
  candidate$objective <- fit$objective - 4
  candidate$model$THETAS$FIX[[1L]] <- FALSE
  comparison <- nm_compare(
    fit, candidate, labels = c("base", "candidate"),
    nested = c(TRUE, TRUE), notes = "unit test"
  )
  expect_s3_class(comparison, "nm_model_comparison")
  expect_equal(comparison$pairwise$delta_objective, 4)
  expect_equal(comparison$pairwise$delta_parameters, 1L)
  expect_equal(comparison$pairwise$decision, "prefer-candidate")
  expect_length(comparison$metrics$fingerprint, 2L)
  path <- tempfile(fileext = ".rds")
  nm_compare_save(comparison, path)
  expect_equal(nm_compare_read(path)$id, comparison$id)
  calibrated <- nm_compare(
    fit, candidate, labels = c("base", "candidate"),
    nested = c(TRUE, TRUE), parametric_bootstrap = 1L, seed = 7L,
    refit_control = list(maxit = 1L)
  )
  expect_equal(calibrated$parametric_bootstrap$candidate$requested, 1L)
  expect_equal(calibrated$parametric_bootstrap$candidate$successful, 1L)
  modified <- comparison
  modified$fits[[2L]]$objective <- modified$fits[[2L]]$objective + 1
  saveRDS(modified, path)
  expect_error(nm_compare_read(path), "integrity")
})

test_that("bootstrap supports strata, clusters, and parametric simulation", {
  fixture <- estimation_fixture()
  fixture$data$SITE <- rep(c("A", "B", "A"), each = 4)
  fixture$data$CLUSTER <- rep(c(1, 2, 3), each = 4)
  fit <- nm_est(fixture$model, fixture$data, method = "FOCEI", maxit = 1)
  stratified <- nm_bootstrap(
    fit, n = 2, strata = "SITE", maxit = 1, seed = 10
  )
  clustered <- nm_bootstrap(
    fit, n = 2, unit = "cluster", cluster = "CLUSTER",
    strata = "SITE", maxit = 1, seed = 11
  )
  parametric <- nm_bootstrap(
    fit, n = 2, type = "parametric", maxit = 1, seed = 12
  )
  expect_s3_class(stratified, "nm_bootstrap")
  expect_equal(clustered$unit, "cluster")
  expect_equal(clustered$strata, "SITE")
  expect_equal(parametric$type, "parametric")
})

test_that("SIR returns weighted and resampled population uncertainty", {
  fixture <- estimation_fixture(fix = FALSE)
  fit <- nm_est(
    fixture$model, fixture$data, method = "FOCEI", maxit = 8,
    covariance = TRUE, covariance_type = "hessian"
  )
  sir <- nm_sir(fit, n_proposal = 10, n_resample = 5, seed = 3)
  expect_s3_class(sir, "nm_sir")
  expect_equal(nrow(sir$draws), 5L)
  expect_equal(sum(sir$weights), 1, tolerance = 1e-12)
  expect_true(sir$diagnostics$effective_sample_size > 1)
  expect_equal(rownames(sir$summary), colnames(sir$draws))
  expect_equal(
    rownames(sir$summary),
    c("THETA1", "THETA2", "OMEGA1", "SIGMA1")
  )
})

test_that("Bayesian pointwise likelihood, WAIC, and PPC use posterior draws", {
  fixture <- estimation_fixture()
  fit <- nm_est(
    fixture$model, fixture$data, method = "BAYES",
    n_burn = 2, n_sample = 5, n_thin = 1, seed = 2
  )
  log_lik <- nm_log_lik(fit, draws = 3, eta_samples = 3, seed = 4)
  expect_equal(dim(log_lik), c(3L, 3L))
  expect_true(isTRUE(attr(log_lik, "marginal_random_effects")))
  waic <- nm_waic(log_lik)
  expect_s3_class(waic, "nm_waic")
  expect_true(is.finite(waic$waic))
  expect_no_warning(ppc <- nm_ppc(fit, draws = 3, seed = 4))
  expect_s3_class(ppc, "nm_ppc")
  expect_equal(nrow(ppc$checks), 5L)
  if (requireNamespace("loo", quietly = TRUE)) {
    expect_s3_class(nm_psis_loo(log_lik), "nm_psis_loo")
  }
})
