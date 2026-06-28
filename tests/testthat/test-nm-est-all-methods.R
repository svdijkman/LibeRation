#' Full algorithm smoke checks on synthetic THEO (perturbed initials, C++ path).
#'
#' Verifies each supported \code{nm_est} method runs without error, returns finite
#' objectives/parameters, and moves at least one free parameter from the start.

.nm_perturb_start <- function(model, factor = 1.25) {
  par <- .nm_init_par(model)
  free <- !.nm_fix_mask(model)
  par[free] <- par[free] * factor
  .nm_apply_fix(model, par)
}

.nm_par_movement <- function(par0, par1) {
  free <- seq_along(par0)
  max(abs(par1[free] - par0[free]))
}

test_that("all estimation methods run on THEO at N=10 with perturbed initials", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 10L, seed = 42L)
  start <- .nm_perturb_start(sim$model, 1.25)
  ctl <- list(maxit = 60L, compute_inference = FALSE, n_cores = 1L)

  fit_fo <- nm_est(
    sim$model, sim$data, start = start, method = "FO",
    grad = "cpp", pk_engine = "cpp", control = ctl
  )
  expect_s3_class(fit_fo, "nm_fit")
  expect_true(is.finite(fit_fo$objective))
  expect_true(.nm_par_movement(start, fit_fo$par) > 1e-4)

  fit_foce <- nm_est(
    sim$model, sim$data, start = start, method = "FOCE",
    grad = "cpp", pk_engine = "cpp", engine = "cpp",
    control = ctl, max_outer = 4L
  )
  expect_s3_class(fit_foce, "nm_fit")
  expect_true(is.finite(fit_foce$objective))
  expect_true(is.matrix(fit_foce$eta))
  expect_true(.nm_par_movement(start, fit_foce$par) > 1e-4)

  fit_focei <- nm_est(
    sim$model, sim$data, start = start, method = "FOCEI",
    grad = "cpp", pk_engine = "cpp", engine = "cpp",
    control = ctl, max_outer = 5L
  )
  expect_s3_class(fit_focei, "nm_fit")
  expect_true(is.finite(fit_focei$objective))
  expect_true(
    .nm_par_movement(start, fit_focei$par) > 1e-4 ||
      (is.matrix(fit_focei$eta) && max(abs(fit_focei$eta), na.rm = TRUE) > 1e-4)
  )

  fit_saem <- nm_est(
    sim$model, sim$data, start = start, method = "SAEM",
    grad = "cpp", pk_engine = "cpp", engine = "cpp",
    n_iter = 30L, n_burn = 8L, seed = 1L, control = ctl
  )
  expect_s3_class(fit_saem, "nm_fit")
  expect_true(is.finite(fit_saem$objective))
  expect_true(.nm_par_movement(start, fit_saem$par) > 1e-3)

  fit_lap <- nm_est(
    sim$model, sim$data, start = start, method = "LAPLACE",
    grad = "cpp", pk_engine = "cpp", engine = "cpp",
    n_quad = 3L, control = ctl
  )
  expect_s3_class(fit_lap, "nm_fit")
  expect_true(is.finite(fit_lap$objective))
  expect_true(is.matrix(fit_lap$eta))
  expect_true(.nm_par_movement(start, fit_lap$par) > 1e-4)

  fit_imp <- nm_est(
    sim$model, sim$data, start = start, method = "IMP",
    grad = "cpp", pk_engine = "cpp", engine = "cpp",
    n_imp = 15L, n_quad = 3L, control = ctl
  )
  expect_s3_class(fit_imp, "nm_fit")
  expect_true(is.finite(fit_imp$objective))
  expect_true(.nm_par_movement(start, fit_imp$par) > 1e-4)

  fit_bayes <- nm_est(
    sim$model, sim$data, start = start, method = "BAYES",
    engine = "cpp", n_burn = 10L, n_sample = 20L, seed = 1L
  )
  expect_s3_class(fit_bayes, "nm_fit")
  expect_true(all(is.finite(fit_bayes$theta)))
  expect_true(is.matrix(fit_bayes$eta))
})

test_that("FOCE, SAEM, and LAPLACE run on THEO at N=20 with perturbed initials", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 20L, seed = 43L)
  start <- .nm_perturb_start(sim$model, 1.2)
  ctl <- list(maxit = 80L, compute_inference = FALSE, n_cores = 1L)

  fit_foce <- nm_est(
    sim$model, sim$data, start = start, method = "FOCE",
    grad = "cpp", pk_engine = "cpp", engine = "cpp",
    control = ctl, max_outer = 4L
  )
  expect_s3_class(fit_foce, "nm_fit")
  expect_true(is.finite(fit_foce$objective))
  expect_true(.nm_par_movement(start, fit_foce$par) > 1e-4)

  fit_saem <- nm_est(
    sim$model, sim$data, start = start, method = "SAEM",
    grad = "cpp", pk_engine = "cpp", engine = "cpp",
    n_iter = 35L, n_burn = 10L, seed = 2L, control = ctl
  )
  expect_s3_class(fit_saem, "nm_fit")
  expect_true(is.finite(fit_saem$objective))
  expect_true(.nm_par_movement(start, fit_saem$par) > 1e-3)

  fit_lap <- nm_est(
    sim$model, sim$data, start = start, method = "LAPLACE",
    grad = "cpp", pk_engine = "cpp", engine = "cpp",
    n_quad = 3L, control = ctl
  )
  expect_s3_class(fit_lap, "nm_fit")
  expect_true(is.finite(fit_lap$objective))
  expect_true(.nm_par_movement(start, fit_lap$par) > 1e-4)
})

test_that("inference objective matches fit objective for FOCE and LAPLACE", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 5L, seed = 44L)
  ctl <- list(maxit = 40L, compute_inference = FALSE, n_cores = 1L)

  fit_foce <- nm_est(
    sim$model, sim$data, method = "FOCE",
    grad = "cpp", pk_engine = "cpp", engine = "cpp",
    control = ctl, max_outer = 3L
  )
  obj_foce <- .nm_fit_inference_objective(fit_foce, fit_foce$par, data = sim$data)
  expect_true(is.finite(obj_foce))
  expect_equal(fit_foce$objective, obj_foce, tolerance = 1e-2)

  fit_lap <- nm_est(
    sim$model, sim$data, method = "LAPLACE",
    grad = "cpp", pk_engine = "cpp", engine = "cpp",
    n_quad = 3L, control = ctl
  )
  obj_lap <- .nm_fit_inference_objective(fit_lap, fit_lap$par, data = sim$data)
  expect_true(is.finite(obj_lap))
  expect_equal(fit_lap$objective, obj_lap, tolerance = 1e-1)
})
