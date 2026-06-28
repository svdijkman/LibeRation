test_that("BAYES MCMC runs on synthetic THEO data", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 3L, seed = 45L)
  fit <- nm_est(
    sim$model, sim$data,
    method = "BAYES",
    engine = "cpp",
    n_burn = 10L,
    n_sample = 20L,
    n_thin = 2L,
    seed = 7L
  )
  expect_s3_class(fit, "nm_fit")
  expect_equal(fit$method, "BAYES")
  expect_true(is.finite(fit$log_posterior))
  expect_true(all(is.finite(fit$theta)))
  expect_true(is.matrix(fit$chains$theta))
  expect_equal(nrow(fit$chains$theta), fit$n_keep)
  expect_true(is.matrix(fit$eta))
  expect_true(!is.null(fit$acceptance$eta))
  expect_true(is.integer(fit$eval$pk) && fit$eval$pk > 0L)
  sm <- summary(fit)
  expect_s3_class(sm, "summary.nm_fit_bayes")
  expect_true(nrow(sm$theta) >= 1L)
})

test_that("BAYES attaches posterior SD and quantile intervals", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 3L, seed = 45L)
  fit <- nm_est(
    sim$model, sim$data,
    method = "BAYES",
    engine = "cpp",
    n_burn = 10L,
    n_sample = 30L,
    n_thin = 1L,
    seed = 7L
  )
  expect_false(is.null(fit$par_se))
  expect_true(any(is.finite(unname(fit$par_se))))
  expect_false(is.null(fit$par_ci_low))
  expect_false(is.null(fit$par_ci_high))
  expect_equal(fit$inference_method, "posterior")
  pt <- nm_fit_param_table(sim$model, fit)
  expect_true(any(is.finite(pt$se)))
})

test_that("min_retries and tweak_inits re-run on failed convergence", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 3L, seed = 1L)
  fit <- nm_est(
    sim$model, sim$data,
    method = "FO",
    control = list(maxit = 1L, min_retries = 1L, tweak_inits = TRUE, compute_inference = FALSE)
  )
  expect_true(!is.null(fit$est_retries) || identical(fit$convergence, 0L))
})

test_that("VPC simulation packs vpc summary for workspace", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 3L, seed = 2L)
  out <- nm_simulate(sim$model, sim$data, n_sim = 5L, seed = 3L, vpc = TRUE)
  expect_true(is.list(out))
  expect_false(is.null(out$vpc))
  expect_false(is.null(out$obs))
  expect_true("TIME" %in% names(out$vpc))
  expect_true(all(c("obs_med", "obs_lo", "obs_hi") %in% names(out$vpc)))
  expect_true(all(c("sim_med_lo", "sim_med_hi", "sim_lo_lo", "sim_hi_hi") %in% names(out$vpc)))
  expect_true(all(c("sim_med_lo_pc", "sim_med_hi_pc") %in% names(out$vpc)))
  packed <- LibeRation:::.nm_sim_pack_output(out)
  expect_true(isTRUE(packed$vpc_mode))
  expect_false(is.null(packed$vpc))
})

test_that("VPC bins handle duplicate observation times", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 3L, seed = 2L)
  d <- sim$data$data
  d$TIME <- 1
  sim2 <- structure(list(model = sim$model, data = structure(list(data = d), class = "nm_dataset")))
  expect_error(
    nm_simulate(sim2$model, sim2$data, n_sim = 3L, seed = 1L, vpc = TRUE),
    NA
  )
})

test_that("template PK matches ADVAN 3 TRANS 1 parameterization", {
  tpl <- nm_ctl_template(3L, 1L)
  expect_true(grepl("K10", tpl$pk, fixed = TRUE))
  expect_false(grepl("VC", tpl$pk, fixed = TRUE))
})

test_that("template PK matches ADVAN 3 TRANS 4 parameterization", {
  tpl <- nm_ctl_template(3L, 4L)
  expect_true(grepl("VC", tpl$pk, fixed = TRUE))
  expect_true(grepl("S2 = VP", tpl$pk, fixed = TRUE))
})

test_that("template PK matches ADVAN 12 TRANS 4 parameterization", {
  tpl <- nm_ctl_template(12L, 4L)
  expect_true(grepl("VP2", tpl$pk, fixed = TRUE))
  expect_true(grepl("Q3", tpl$pk, fixed = TRUE))
  expect_true(grepl("KA", tpl$pk, fixed = TRUE))
})

test_that("BAYES prior spec accepts lognormal", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 2L, seed = 46L)
  p0 <- .nm_init_par(sim$model)
  prior <- list(
    theta = list(type = "lognormal", meanlog = log(sim$model$THETAS$Value), sdlog = 0.5)
  )
  fit <- nm_est(
    sim$model, sim$data,
    method = "BAYES",
    engine = "cpp",
    n_burn = 5L,
    n_sample = 10L,
    seed = 3L,
    prior = prior
  )
  expect_s3_class(fit, "nm_fit")
  expect_equal(fit$method, "BAYES")
})
