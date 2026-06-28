test_that("bootstrap SE is SD of bootstrap estimates, not mean of replicate SEs", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_warf(n_sub = 6L, seed = 2L)
  fit <- nm_est(
    sim$model, sim$data,
    method = "FOCE",
    grad = "numeric",
    pk_engine = "cpp",
    control = list(
      maxit = 30L,
      compute_inference = TRUE,
      infer_hessian = "numeric"
    )
  )
  boot <- nm_bootstrap_se(
    fit,
    n_boot = 6L,
    seed = 3L,
    control = list(maxit = 15L, compute_inference = TRUE, infer_hessian = "numeric")
  )
  expect_equal(boot$se_method, "sd_of_bootstrap_estimates")
  labels <- .nm_par_labels(fit$model)
  for (lbl in labels) {
    x <- boot$boot_pars[, lbl]
    x <- x[is.finite(x)]
    if (length(x) >= 2L) {
      expect_equal(boot$se[[lbl]], stats::sd(x), tolerance = 1e-10)
    }
  }
})

test_that("predict ipred refits missing etas", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_digox(n_sub = 4L, seed = 1L)
  fit <- nm_est(
    sim$model, sim$data,
    method = "FOCE",
    grad = "numeric",
    pk_engine = "cpp",
    control = list(maxit = 25L, compute_inference = FALSE)
  )
  fit$eta <- NULL
  pred <- predict(fit, type = "ipred")
  pop <- predict(fit, type = "ppred")
  obs <- pred[pred$MDV == 0L & pred$EVID == 0L, ]
  obs2 <- pop[pop$MDV == 0L & pop$EVID == 0L, ]
  expect_lt(sum(abs(obs$IPRED - obs2$PRED) < 1e-10), nrow(obs))
})

test_that("nm_ctl_canonical normalizes equivalent control text", {
  parts <- nm_ctl_template(advan = 2L, trans = 1L, data_file = "data/x.csv", problem = "Test")
  a <- nm_ctl_compose(parts)
  b <- paste0(a, "\n")
  expect_identical(nm_ctl_canonical(a), nm_ctl_canonical(b))
})
