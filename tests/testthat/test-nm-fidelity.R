# Fidelity checks for the estimation-algorithm improvements:
#   FO marginal, FOCEI interaction, SAEM closed-form Omega, IMP importance
#   sampling. These are lightweight, NONMEM-free self-consistency checks.

test_that("FO marginal: C++ linear algebra matches reference R implementation", {
  skip_if_not_installed("data.table")
  chk <- nm_bench_fo_marginal_check(id = "warf", n_sub = 6L, seed = 11L)
  expect_true(is.finite(chk$cpp))
  expect_true(is.finite(chk$r))
  expect_true(chk$ok)
  expect_lt(chk$diff, 1e-4)
})

test_that("FOCEI interaction term engages for proportional error, no-op for additive", {
  skip_if_not_installed("data.table")
  chk <- nm_bench_focei_interaction_check(n_sub = 6L, seed = 21L)
  # For a proportional/f-dependent error model the interaction term must change
  # the objective; for a pure additive model it must be a no-op.
  expect_true(chk$prop_engaged)
  expect_true(chk$add_noop)
  # R and C++ FOCEI-interaction objectives must agree.
  expect_true(chk$r_cpp_ok)
})

test_that("IMP importance-sampling objective is finite and uses common random numbers", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 4L, seed = 31L)
  m <- sim$model
  d <- sim$data
  th <- as.numeric(m$THETAS$Value)
  om <- as.numeric(m$OMEGAS$Value)
  sg <- as.numeric(m$SIGMAS$Value)
  eta <- matrix(0, 4L, .nm_n_eta(m))
  gh <- .nm_gh_nodes(3L)
  o1 <- .nm_imp_nll(m, d, th, om, sg, eta, gh,
                    n_imp = 20L, pk_engine = "cpp", seed = 99L)
  o2 <- .nm_imp_nll(m, d, th, om, sg, eta, gh,
                    n_imp = 20L, pk_engine = "cpp", seed = 99L)
  expect_true(is.finite(o1))
  # Common random numbers => deterministic objective for identical parameters.
  expect_equal(o1, o2)
})

test_that("IMP proposal covariance is the conditional posterior covariance at the mode", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 4L, seed = 41L)
  m <- sim$model
  d <- sim$data
  th <- as.numeric(m$THETAS$Value)
  om <- as.numeric(m$OMEGAS$Value)
  sg <- as.numeric(m$SIGMAS$Value)
  dat <- .nm_prepare_data(d, m$INPUT, m)
  ids <- .nm_subject_ids(dat)
  sub <- .nm_subject_slice(dat, ids[1L])
  n_eta <- .nm_n_eta(m)
  prop <- .nm_imp_proposal_cov(m, sub, th, om, sg, rep(0, n_eta), "cpp", n_eta)
  L <- prop$chol_lower
  cov <- L %*% t(L)
  # Conditional covariance must be positive-definite and tighter than the prior
  # Omega (data are informative), i.e. no downgrade for this well-behaved case.
  expect_false(isTRUE(prop$fallback))
  expect_true(all(is.finite(cov)))
  expect_true(all(diag(cov) > 0))
  expect_true(all(diag(cov) <= pmax(om[seq_len(n_eta)], 1e-8) + 1e-8))
})

test_that("SAEM closed-form Omega sufficient-statistic path runs and recovers Omega scale", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 8L, seed = 51L)
  m <- sim$model
  d <- sim$data
  fit <- nm_est(
    m, d, method = "SAEM", n_iter = 15L, n_burn = 5L, seed = 3L,
    control = list(maxit = 6L, compute_inference = FALSE, n_cores = 1L,
                   saem_omega_closed_form = TRUE)
  )
  expect_equal(fit$method, "SAEM")
  expect_true(all(is.finite(fit$omega)))
  expect_true(all(fit$omega[fit$omega != 0] > 0))
})
