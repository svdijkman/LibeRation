test_that("POSTHOC fixes population parameters and computes ETAs", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 4L, seed = 21L)
  fit <- nm_est(sim$model, sim$data, method = "FOCE", grad = "numeric",
                control = list(maxit = 15L, compute_inference = FALSE))
  m2 <- sim$model
  m2$THETAS$Value <- fit$theta
  m2$OMEGAS$Value <- fit$omega
  m2$SIGMAS$Value <- fit$sigma
  ph <- nm_est(m2, sim$data, method = "POSTHOC")
  expect_equal(ph$method, "POSTHOC")
  expect_equal(ph$theta, fit$theta)          # population unchanged
  expect_true(is.finite(ph$objective))
  expect_equal(dim(ph$eta), c(4L, 3L))
  expect_true(all(is.finite(ph$eta)))
  # ETAs should be close to the FOCE post-hoc ETAs at the same parameters.
  expect_lt(max(abs(ph$eta - fit$eta)), 0.05)
})

test_that("control maxeval = 0 and posthoc = TRUE route to POSTHOC", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 3L, seed = 22L)
  ph1 <- nm_est(sim$model, sim$data, control = list(maxeval = 0L))
  expect_equal(ph1$method, "POSTHOC")
  ph2 <- nm_est(sim$model, sim$data, method = "FOCE",
                control = list(posthoc = TRUE))
  expect_equal(ph2$method, "POSTHOC")
})

test_that("nm_write_table writes and round-trips a NONMEM-style table", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 3L, seed = 23L)
  fit <- nm_est(sim$model, sim$data, method = "FOCE", grad = "numeric",
                control = list(maxit = 12L, compute_inference = FALSE))
  f <- tempfile(fileext = ".tab")
  tab <- nm_write_table(fit, f, nonmem_header = TRUE)
  expect_true(all(c("ID", "TIME", "DV", "PRED", "IPRED", "RES", "WRES",
                    "CWRES", "ETA1", "ETA2", "ETA3") %in% names(tab)))
  back <- nm_read_table(f)
  expect_equal(nrow(back), nrow(tab))
  expect_equal(names(back), names(tab))
  expect_equal(as.numeric(back$IPRED), as.numeric(tab$IPRED), tolerance = 1e-4)
})

test_that("non-finite objective guard maps to a large finite penalty", {
  pen <- LibeRation:::.nm_finite_penalty
  expect_equal(pen(NaN), .Machine$double.xmax)
  expect_equal(pen(Inf), .Machine$double.xmax)
  expect_equal(pen(-Inf), .Machine$double.xmax)
  expect_equal(pen(numeric(0)), .Machine$double.xmax)
  expect_equal(pen(3.5), 3.5)
})

test_that("nm_check_pk_engines agrees for a supported model", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 3L, seed = 25L)
  chk <- nm_check_pk_engines(sim$model, sim$data)
  expect_true(chk$comparable)
  expect_false(chk$diverged)
  expect_lt(chk$max_rel_diff, 1e-4)
  expect_equal(nrow(chk$per_subject), 3L)
})

test_that("nm_write_table supports firstonly and custom columns", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 4L, seed = 24L)
  fit <- nm_est(sim$model, sim$data, method = "FO", grad = "numeric",
                control = list(maxit = 10L, compute_inference = FALSE))
  first <- nm_write_table(fit, NULL, firstonly = TRUE)
  expect_equal(nrow(first), 4L)
  custom <- suppressWarnings(
    nm_write_table(fit, NULL, columns = c("ID", "TIME", "DV", "CWRES", "NOPE"))
  )
  expect_equal(names(custom), c("ID", "TIME", "DV", "CWRES"))
})
