test_that("oral synthetic examples estimate with ADVAN 2 TRANS 1 K/V parameterization", {
  skip_if_not_installed("data.table")
  ids <- c("warf", "pheno", "digox")
  for (id in ids) {
    sim <- nm_synthetic_dataset(id = id, n_sub = 12L, seed = 42L)
    expect_equal(sim$model$ADVAN, 2L)
    expect_equal(sim$model$TRANS, 1L)
    expect_true(grepl("\\bK\\s*=", sim$model$PRED, ignore.case = FALSE))
    expect_true(grepl("S2\\s*=\\s*V", sim$model$PRED))
    fit <- nm_est(
      sim$model, sim$data,
      method = "FOCE", grad = "cpp", pk_engine = "cpp",
      control = list(maxit = 200)
    )
    th <- sim$model$THETAS$Value
    rel <- abs(fit$theta[1:3] - th[1:3]) / pmax(abs(th[1:3]), 1e-8)
    expect_true(max(rel) < 0.35, info = paste(id, "theta rel err", max(rel)))
  }
})

test_that("digoxin synthetic shows non-flat PK concentrations", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_digox(n_sub = 3L, seed = 1L)
  dv <- sim$data$data[sim$data$data$MDV == 0L & sim$data$data$EVID == 0L, ]
  expect_true(diff(range(dv$DV, na.rm = TRUE)) > 0.05)
})

test_that("predict refits missing or all-zero subject ETAs", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_digox(n_sub = 8L, seed = 2L)
  fit <- nm_est(
    sim$model, sim$data,
    method = "FOCE", grad = "numeric", pk_engine = "cpp",
    control = list(maxit = 80L, compute_inference = FALSE)
  )
  fit$eta <- matrix(0, nrow = nrow(fit$eta), ncol = ncol(fit$eta))
  eta_ref <- LibeRation:::.nm_fit_eta_matrix(fit, data = fit$data)
  expect_true(max(abs(eta_ref), na.rm = TRUE) > 1e-6)
  pred <- predict(fit, type = "ipred")
  pop <- predict(fit, type = "ppred")
  obs <- pred[pred$MDV == 0L & pred$EVID == 0L, ]
  obs2 <- pop[pop$MDV == 0L & pop$EVID == 0L, ]
  expect_true(max(abs(obs$IPRED - obs2$PRED), na.rm = TRUE) > 1e-6)
})

test_that("digoxin FOCE works with many subjects (cache key length)", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_digox(n_sub = 100L, seed = 3L)
  key <- LibeRation:::.nm_data_subject_cache_key(
    LibeRation:::.nm_prepare_data(sim$data, sim$model$INPUT, sim$model)
  )
  expect_true(nchar(key) < 1000L)
  fit <- nm_est(
    sim$model, sim$data,
    method = "FOCE", grad = "numeric", pk_engine = "cpp",
    control = list(maxit = 15L, compute_inference = FALSE, n_cores = 1L)
  )
  expect_true(max(abs(fit$eta), na.rm = TRUE) > 1e-6)
})

test_that("digoxin synthetic recovers visible random effects at n=20", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_digox(n_sub = 20L, seed = 4L)
  fit <- nm_est(
    sim$model, sim$data,
    method = "FOCE", grad = "numeric", pk_engine = "cpp",
    control = list(maxit = 100L, compute_inference = FALSE)
  )
  expect_true(max(apply(fit$eta, 2, sd), na.rm = TRUE) > 0.15)
  pred <- predict(fit, type = "ipred")
  pop <- predict(fit, type = "ppred")
  obs <- pred[pred$MDV == 0L & pred$EVID == 0L, ]
  obs2 <- pop[pop$MDV == 0L & pop$EVID == 0L, ]
  expect_true(max(abs(obs$IPRED - obs2$PRED), na.rm = TRUE) > 0.05)
})
