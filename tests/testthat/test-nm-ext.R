test_that("nm_lik_config and engine detail work", {
  cfg <- nm_lik_config(error = "log", omega = "block2")
  expect_s3_class(cfg, "nm_lik_config")
  expect_equal(cfg$error_code, 3L)
  det <- LibeRation:::.nm_engine_detail("SAEM", "numeric", "numeric", "cpp", "cpp", "cpp", cfg)
  expect_s3_class(det, "nm_engine_detail")
})

test_that("predict and nm_etab run on FO fit", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 3L, seed = 50L)
  fit <- nm_est(sim$model, sim$data, method = "FO", grad = "numeric",
                pk_engine = "cpp", control = list(maxit = 15))
  pred <- predict(fit)
  expect_true("IWRES" %in% names(pred))
  etab <- nm_etab(fit)
  expect_equal(nrow(etab), 3L)
  sh <- nm_shrinkage(fit)
  expect_true(is.data.frame(sh))
})

test_that("FOCEI and IMP methods run", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 3L, seed = 51L)
  fit_i <- nm_est(sim$model, sim$data, method = "FOCEI", grad = "numeric",
                  pk_engine = "cpp", max_outer = 2L, control = list(maxit = 10))
  expect_equal(fit_i$method, "FOCEI")
  expect_equal(fit_i$focei_eta, "nested")
  fit_imp <- nm_est(sim$model, sim$data, method = "IMP", grad = "numeric",
                    pk_engine = "cpp", n_imp = 5L, control = list(maxit = 10))
  expect_equal(fit_imp$method, "IMP")
})

test_that("FOCEI focei_eta=outer runs with fixed eta during inner optim", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 3L, seed = 55L)
  fit <- nm_est(
    sim$model, sim$data, method = "FOCEI", grad = "auto", pk_engine = "cpp",
    max_outer = 2L,
    control = list(maxit = 12, focei_eta = "outer", compute_inference = FALSE)
  )
  expect_equal(fit$method, "FOCEI")
  expect_equal(fit$focei_eta, "outer")
  expect_true(is.finite(fit$objective))
  expect_true(is.matrix(fit$eta))
})

test_that("nm_ctl_apply_fit_inits updates control stream values", {
  sim <- nm_synthetic_theo(n_sub = 2L, seed = 56L)
  fit <- nm_est(
    sim$model, sim$data, method = "FO", grad = "numeric", pk_engine = "cpp",
    control = list(maxit = 8, compute_inference = FALSE)
  )
  ctl <- nm_ctl_compose(
    list(
      problem = "test", advan = 2L, trans = 2L, data_file = "d.csv",
      thetas = sim$model$THETAS,
      omegas = sim$model$OMEGAS,
      sigmas = sim$model$SIGMAS,
      pk = " CL = THETA(1)\n V = THETA(2)",
      error = "Y = F + D"
    )
  )
  parts <- nm_ctl_parse(ctl)
  parts2 <- nm_ctl_apply_fit_inits(parts, fit)
  pp <- LibeRation:::.nm_unpack(fit$model, fit$par)
  expect_equal(parts2$thetas$Value, pp$theta)
})

test_that("nm_proj save/load and support matrix", {
  p <- nm_proj("test")
  p <- nm_proj_set(p)
  expect_s3_class(p, "nm_proj")
  sm <- nm_support_matrix()
  expect_true("FOCEI" %in% sm$method)
})

test_that("nm_read_nonmem parses minimal ctl", {
  ctl <- tempfile(fileext = ".ctl")
  writeLines(c(
    "$PROB test ADVAN=2 TRANS=2",
    "$INPUT ID TIME DV AMT",
    "$DATA theo.csv",
    "$THETA (0,1)",
    "$OMEGA 0.1",
    "$SIGMA 0.1",
    "$PK CL = THETA(1)",
    "$ERROR Y = F"
  ), ctl)
  x <- nm_read_nonmem(ctl, data_path = "dummy.csv")
  expect_equal(x$method, "FO")
  expect_equal(trimws(x$model$PRED), "CL = THETA(1)")
  expect_true(nrow(x$model$THETAS) >= 1L)
  unlink(ctl)
})

test_that("nm_mcmc_diagnostics on BAYES fit", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 2L, seed = 52L)
  fit <- nm_est(sim$model, sim$data, method = "BAYES", engine = "cpp",
                n_burn = 5L, n_sample = 10L, seed = 1L)
  d <- nm_mcmc_diagnostics(fit)
  expect_s3_class(d, "nm_mcmc_diagnostics")
})

test_that("nm_lik_config AR1 and IOV sync to C++", {
  cfg <- nm_lik_config(error = "propadd", sigma_corr = "ar1", ar1_rho = 0.3, iov = 1L)
  expect_equal(cfg$sigma_corr_code, 1L)
  expect_equal(cfg$iov, 1L)
})

test_that("nm_to_nlmixr2 and nm_from_nlmixr2 roundtrip basics", {
  sim <- nm_synthetic_theo(n_sub = 2L, seed = 53L)
  txt <- nm_to_nlmixr2(sim$model, name = "m1")
  expect_true(any(grepl("function", txt)))
  imp <- nm_from_nlmixr2(list(theta = sim$model$THETAS$Value), model = sim$model)
  expect_true(!is.null(imp$par))
})

test_that("LibeRtAD tape cache and sparse jacobian", {
  skip_if_not_installed("LibeRtAD")
  f <- function(x) x^2
  LibeRtAD::autodiff(f, at = list(x = 2), mode = "reverse", record_tape = TRUE)
  saved <- LibeRtAD::ad_tape_save("t1")
  expect_true(saved$n_nodes > 0L)
  loaded <- LibeRtAD::ad_tape_load("t1")
  expect_equal(loaded$n_nodes, saved$n_nodes)
  sj <- LibeRtAD::ad_sparse_jacobian(f, x = 2)
  expect_true(length(sj$x) >= 1L)
})

test_that("BAYES HMC sampler runs", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 2L, seed = 54L)
  fit <- nm_est(sim$model, sim$data, method = "BAYES", engine = "cpp",
                sampler = "hmc", n_burn = 3L, n_sample = 6L, seed = 1L)
  expect_equal(fit$sampler, "hmc")
})

test_that("LibeRtAD hessian runs", {
  skip_if_not_installed("LibeRtAD")
  if (!exists("autodiff_hessian", where = asNamespace("LibeRtAD"), inherits = FALSE)) {
    skip("autodiff_hessian not exported")
  }
  f <- function(x, y) x^2 + x * y + y^2
  h <- LibeRtAD::autodiff_hessian(f, x = 1, y = 2)
  expect_equal(h$hessian[1, 1], 2, tolerance = 0.01)
})
