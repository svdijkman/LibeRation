test_that("nm_lik_config carries BLQ method and code", {
  cfg <- nm_lik_config(error = "prop", blq_method = "m3", lloq = 0.1)
  expect_s3_class(cfg, "nm_lik_config")
  expect_equal(cfg$blq_method, "m3")
  expect_equal(cfg$blq_code, 1L)
  expect_equal(cfg$lloq, 0.1)
  expect_equal(nm_lik_config(blq_method = "m4")$blq_code, 2L)
  expect_equal(nm_lik_config()$blq_code, 0L)
})

test_that("C++ lik config accepts and reports blq_method", {
  on.exit(LibeRation:::nm_lik_config_set(0L, 0L, 0L, 0L, 0.0, 0L), add = TRUE)
  LibeRation:::nm_lik_config_set(2L, 0L, 0L, 0L, 0.0, 2L)
  got <- LibeRation:::nm_lik_config_get()
  expect_equal(got$error_type, 2L)
  expect_equal(got$blq_method, 2L)
})

test_that("M3 residual matches -2 log Phi((LLOQ - f)/sqrt(R))", {
  blq <- LibeRation:::.nm_residual_nll_blq
  f <- c(5, 2)
  dv <- c(5, 0)
  sigma <- c(0.2, 0.3)  # propadd: var = (s1 f)^2 + s2^2
  lloq <- c(NA, 1.5)
  cens <- c(FALSE, TRUE)
  got <- blq(dv, f, sigma, error = "propadd", lloq = lloq, cens = cens,
             blq_method = "m3")
  var1 <- (0.2 * f[1])^2 + 0.3^2
  term_unc <- log(var1) + (dv[1] - f[1])^2 / var1
  var2 <- (0.2 * f[2])^2 + 0.3^2
  term_cens <- -2 * log(pnorm((lloq[2] - f[2]) / sqrt(var2)))
  expect_equal(got, term_unc + term_cens, tolerance = 1e-8)
})

test_that("M4 conditions on y > 0 and reduces to M3 when Phi(-f/sd) ~ 0", {
  blq <- LibeRation:::.nm_residual_nll_blq
  f <- 1.0
  sigma <- c(0.5, 0.2)
  lloq <- 1.2
  var <- (0.5 * f)^2 + 0.2^2
  sd <- sqrt(var)
  p_lloq <- pnorm((lloq - f) / sd)
  p_zero <- pnorm((0 - f) / sd)
  exp_m4 <- -2 * log((p_lloq - p_zero) / (1 - p_zero))
  got_m4 <- blq(0, f, sigma, error = "propadd", lloq = lloq, cens = TRUE,
                blq_method = "m4")
  expect_equal(got_m4, exp_m4, tolerance = 1e-8)
  # With f large relative to sd, p_zero ~ 0 so M4 ~ M3.
  got_m3 <- blq(0, 8, c(0.1, 0.1), error = "propadd", lloq = 1, cens = TRUE,
                blq_method = "m3")
  got_m4b <- blq(0, 8, c(0.1, 0.1), error = "propadd", lloq = 1, cens = TRUE,
                 blq_method = "m4")
  expect_equal(got_m3, got_m4b, tolerance = 1e-6)
})

test_that("no censored rows reproduces the standard residual", {
  blq <- LibeRation:::.nm_residual_nll_blq
  scal <- LibeRation:::.nm_residual_nll_scalar
  dv <- c(3, 4, 5)
  f <- c(3.1, 3.9, 5.2)
  sigma <- c(0.15, 0.25)
  a <- blq(dv, f, sigma, error = "propadd",
           lloq = rep(NA, 3), cens = rep(FALSE, 3), blq_method = "m3")
  b <- scal(dv, f, sigma, error = "propadd")
  expect_equal(a, b, tolerance = 1e-10)
})

test_that(".nm_blq_obs_info honours BLQ, CENS and DV<LLOQ conventions", {
  info <- LibeRation:::.nm_blq_obs_info
  cfg0 <- nm_lik_config()
  ev1 <- data.frame(DV = c(1, 0.05), LLOQ = c(0.1, 0.1), BLQ = c(0L, 1L))
  r1 <- info(ev1, seq_len(2), cfg0)
  expect_equal(r1$cens, c(FALSE, TRUE))
  ev2 <- data.frame(DV = c(1, 0.05), LLOQ = c(0.1, 0.1), CENS = c(0L, 1L))
  r2 <- info(ev2, seq_len(2), cfg0)
  expect_equal(r2$cens, c(FALSE, TRUE))
  # No flag column -> derive from DV < LLOQ
  ev3 <- data.frame(DV = c(1, 0.05), LLOQ = c(0.1, 0.1))
  r3 <- info(ev3, seq_len(2), cfg0)
  expect_equal(r3$cens, c(FALSE, TRUE))
  # Scalar LLOQ fallback from config when no column present
  ev4 <- data.frame(DV = c(1, 0.05))
  r4 <- info(ev4, seq_len(2), nm_lik_config(lloq = 0.1))
  expect_equal(r4$cens, c(FALSE, TRUE))
})

test_that("nm_nll wires BLQ end-to-end through data prep", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_iv1(n_sub = 3L, seed = 3L)
  model <- sim$model
  d <- data.table::as.data.table(sim$data$data)
  d <- d[d$ID == 1L]
  obs <- d$EVID == 0L & d$MDV == 0L
  lloq <- as.numeric(stats::quantile(d$DV[obs], 0.5, na.rm = TRUE))
  d$LLOQ <- lloq
  d$BLQ <- as.integer(obs & d$DV < lloq)
  th <- model$THETAS$Value
  om <- model$OMEGAS$Value
  sg <- model$SIGMAS$Value
  err <- LibeRation:::.nm_lik_config(model)$error

  m3 <- model
  m3$LIK_CONFIG <- nm_lik_config(error = err, blq_method = "m3")

  # Manual expected residual contribution at eta = 0.
  dat <- LibeRation:::.nm_prepare_data(d, model$INPUT, model)
  subj <- LibeRation:::.nm_subject_slice(dat, 1L)
  n_eta <- LibeRation:::.nm_n_eta(model)
  pred <- LibeRation:::.nm_subject_ipred(
    model, subj, th, om, rep(0, n_eta), pk_engine = "cpp"
  )
  dv <- pred$subj_ev$DV[pred$obs_idx]
  binfo <- LibeRation:::.nm_blq_obs_info(
    pred$subj_ev, pred$obs_idx, LibeRation:::.nm_lik_config(m3)
  )
  expect_true(any(binfo$cens))
  expected <- LibeRation:::.nm_residual_nll_blq(
    dv, pred$F, sg, error = err, lloq = binfo$lloq, cens = binfo$cens,
    blq_method = "m3"
  )
  got <- nm_nll(m3, d, th, om, sg, include_omega_prior = FALSE, pk_engine = "cpp")
  expect_equal(as.numeric(got), as.numeric(expected), tolerance = 1e-6)

  # BLQ objective differs from the naive (no-censoring) objective.
  got_none <- nm_nll(model, d, th, om, sg, include_omega_prior = FALSE,
                     pk_engine = "cpp")
  expect_false(isTRUE(all.equal(as.numeric(got), as.numeric(got_none))))
  expect_true(is.finite(as.numeric(got)))
})

test_that("M3 estimation runs end-to-end and returns a finite objective", {
  skip_on_cran()
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_iv1(n_sub = 5L, seed = 5L)
  model <- sim$model
  d <- data.table::as.data.table(sim$data$data)
  obs <- d$EVID == 0L & d$MDV == 0L
  lloq <- as.numeric(stats::quantile(d$DV[obs], 0.3, na.rm = TRUE))
  d$LLOQ <- lloq
  d$BLQ <- as.integer(obs & d$DV < lloq)
  model$LIK_CONFIG <- nm_lik_config(
    error = LibeRation:::.nm_lik_config(model)$error, blq_method = "m3"
  )
  fit <- nm_est(
    model, d, method = "FOCE",
    control = list(maxit = 15L, compute_inference = FALSE, n_cores = 1L),
    max_outer = 2L, tol = 1e-2
  )
  expect_s3_class(fit, "nm_fit")
  expect_true(is.finite(fit$objective))
})
