test_that("nm_cov_step attaches covariance matrices for FOCE", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_warf(n_sub = 8L, seed = 42L)
  fit <- nm_est(
    sim$model, sim$data,
    method = "FOCE",
    grad = "numeric",
    pk_engine = "cpp",
    control = list(
      maxit = 40L,
      compute_inference = TRUE,
      cov_method = "linfim",
      infer_hessian = "numeric"
    )
  )
  expect_false(is.null(fit$covariance))
  expect_false(is.null(fit$covariance$linfim))
  expect_false(is.null(fit$covariance$hessian))
  expect_false(is.null(fit$covariance$sandwich))
  expect_equal(fit$covariance_method, "linfim")
  expect_false(is.null(fit$par_se))
  expect_true(any(is.finite(unname(fit$par_se))))
  cor <- nm_fit_correlation(fit, type = "linfim")
  expect_true(is.matrix(cor))
  expect_equal(unname(diag(cor)), rep(1, nrow(cor)), tolerance = 1e-8)
})

test_that("nm_cov_step can be run standalone after estimation", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_pheno(n_sub = 6L, seed = 1L)
  fit <- nm_est(
    sim$model, sim$data,
    method = "FOCE",
    grad = "numeric",
    pk_engine = "cpp",
    control = list(maxit = 30L, compute_inference = FALSE)
  )
  expect_null(fit$covariance)
  fit2 <- nm_cov_step(fit, method = "sandwich", hessian = "numeric")
  expect_false(is.null(fit2$covariance$sandwich))
  expect_equal(fit2$covariance_method, "sandwich")
  expect_true(any(is.finite(unname(fit2$par_se))))
})
