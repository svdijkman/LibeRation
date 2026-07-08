test_that("timing profiler renamed to nm_time_profile", {
  expect_true(is.function(nm_time_profile))
  # The old S3 profile method for nm_fit must be gone (name freed).
  expect_false("nm_fit" %in% attr(utils::methods("profile"), "info")$to)
  fit_stub <- structure(list(profile = NULL), class = "nm_fit")
  expect_null(nm_time_profile(fit_stub))
  fit_stub2 <- structure(
    list(profile = data.frame(label = "x", time_sec = 1, count = 1L)),
    class = "nm_fit"
  )
  tp <- nm_time_profile(fit_stub2)
  expect_s3_class(tp, "nm_profile")
})

test_that(".nm_profile_ll_theta_index resolves labels and indices", {
  sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
  m <- sim$model
  f <- LibeRation:::.nm_profile_ll_theta_index
  expect_equal(f(m, 2L), 2L)
  expect_equal(f(m, "THETA3"), 3L)
  expect_error(f(m, "NOPE"))
})

test_that("likelihood CI interpolation brackets the estimate", {
  ci <- LibeRation:::.nm_profile_ll_ci
  grid <- seq(1, 5, by = 0.5)
  # parabola minimum at 3, min OFV 100
  ofv <- 100 + 4 * (grid - 3)^2
  est <- 3
  out <- ci(grid, ofv, est, target = 100 + 3.84)
  expect_true(is.finite(out[["lower"]]))
  expect_true(is.finite(out[["upper"]]))
  expect_lt(out[["lower"]], est)
  expect_gt(out[["upper"]], est)
  # analytic crossing: 4 (x-3)^2 = 3.84 -> x = 3 +/- sqrt(0.96)
  expect_equal(out[["lower"]], 3 - sqrt(0.96), tolerance = 0.05)
  expect_equal(out[["upper"]], 3 + sqrt(0.96), tolerance = 0.05)
})

test_that("nm_profile_likelihood profiles a THETA and brackets the estimate", {
  skip_on_cran()
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 8L, seed = 4L)
  fit <- nm_est(sim$model, sim$data, method = "FO", grad = "numeric",
                control = list(maxit = 25L, compute_inference = FALSE))
  pl <- nm_profile_likelihood(
    fit, "THETA1", n = 5L, span = 0.3,
    control = list(maxit = 25L)
  )
  expect_s3_class(pl, "nm_profile_likelihood")
  expect_equal(nrow(pl$profile), 5L)
  expect_true(all(c("value", "ofv", "dofv") %in% names(pl$profile)))
  expect_true(is.finite(pl$ci[["lower"]]))
  expect_true(is.finite(pl$ci[["upper"]]))
  expect_lt(pl$ci[["lower"]], pl$estimate)
  expect_gt(pl$ci[["upper"]], pl$estimate)
})
