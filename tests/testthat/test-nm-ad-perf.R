test_that("AD tape reuse matches full rebuild for population objective", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 3L, seed = 71L)
  m <- sim$model
  d <- sim$data
  par <- .nm_init_par(m)
  eta <- matrix(0, 3L, .nm_n_eta(m))
  pop <- .nm_build_pop_objective(m, d, eta, pk_engine = "cpp")
  labels <- .nm_par_labels(m)
  tape_key <- .nm_ad_tape_key(m, d, eta, "pop")
  at1 <- stats::setNames(as.list(par), labels)
  par2 <- par * 1.05
  at2 <- stats::setNames(as.list(par2), labels)
  withr::local_options(list(LibeRtAD.tape_reuse = TRUE))
  .nm_state$optim_cache <- NULL
  g1 <- .nm_ad_eval_cached(
    pop$fn, at1, labels, "cpp", need_grad = TRUE, tape_key = tape_key
  )
  g2 <- .nm_ad_eval_cached(
    pop$fn, at2, labels, "cpp", need_grad = TRUE, tape_key = tape_key
  )
  expect_true(all(is.finite(g1)))
  expect_true(all(is.finite(g2)))
  expect_true(any(abs(g2) > 0))
  .nm_state$optim_cache <- NULL
  g2_fresh <- .nm_do_call_autodiff(pop$fn, at2, "cpp", tape_key = NULL)
  expect_equal(unname(g2), unname(g2_fresh[labels]), tolerance = 1e-4)
})

test_that("AD tape key invalidates optim cache when ETAs change", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 3L, seed = 71L)
  m <- sim$model
  d <- sim$data
  par <- .nm_init_par(m)
  eta1 <- matrix(0, 3L, .nm_n_eta(m))
  eta2 <- matrix(rnorm(3L * .nm_n_eta(m), sd = 0.05),
                 nrow = 3L, ncol = .nm_n_eta(m))
  pop <- .nm_build_pop_objective(m, d, eta1, pk_engine = "cpp")
  labels <- .nm_par_labels(m)
  at <- stats::setNames(as.list(par), labels)
  key1 <- .nm_ad_tape_key(m, d, eta1, "pop")
  key2 <- .nm_ad_tape_key(m, d, eta2, "pop")
  expect_false(identical(key1, key2))
  .nm_state$optim_cache <- NULL
  pop$ctx$eta_mat <- eta1
  v1 <- .nm_ad_eval_cached(
    pop$fn, at, labels, "cpp", need_grad = FALSE, tape_key = key1
  )
  pop$ctx$eta_mat <- eta2
  v2 <- .nm_ad_eval_cached(
    pop$fn, at, labels, "cpp", need_grad = FALSE, tape_key = key2
  )
  expect_true(is.finite(v1))
  expect_true(is.finite(v2))
  expect_false(isTRUE(all.equal(v1, v2)))
})

test_that("FOCEI sensitivity AD grad agrees with numeric", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 3L, seed = 72L)
  m <- sim$model
  d <- sim$data
  setup <- nm_focei_setup(m, d, pk_engine = "cpp")
  par <- .nm_init_par(m)
  nested <- .nm_focei_nested_objective(
    m, setup$dat, d, par, NULL, "cpp", list(maxit_eta = 20L), use_cache = FALSE
  )
  g_ad <- .nm_focei_sensitivity_grad(
    m, d, par, nested$eta, "cpp", backend = "cpp", grad = "ad"
  )
  g_num <- .nm_focei_sensitivity_grad(
    m, d, par, nested$eta, "cpp", backend = "cpp", grad = "numeric"
  )
  expect_true(all(is.finite(g_ad)))
  expect_equal(g_ad, g_num, tolerance = 0.05)
})
