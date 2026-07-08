test_that("NONMEM ** power is normalized for R PRED evaluation", {
  sim <- nm_synthetic_theo(n_sub = 2L)
  model <- sim$model
  model$PRED <- paste(
    "CL = THETA(1) * exp(ETA(1))",
    "V = THETA(2) ** 2",
    sep = "\n"
  )
  out <- LibeRation:::.nm_eval_pred(
    model,
    model$THETAS$Value,
    model$OMEGAS$Value,
    rep(0, LibeRation:::.nm_n_eta(model))
  )
  expect_equal(as.numeric(out$V), model$THETAS$Value[2L]^2, tolerance = 1e-10)
})

test_that("nm_ctl_pk_highlight_symbols returns ADVAN-specific parameters", {
  hl <- nm_ctl_pk_highlight_symbols(4L, 4L)
  expect_true("CL" %in% hl$pk)
  expect_true("VC" %in% hl$pk)
  expect_true("Q2" %in% hl$flows || "Q2" %in% hl$pk)
})

test_that("nm_ctl_error_highlight_symbols returns error block tokens", {
  hl <- nm_ctl_error_highlight_symbols()
  expect_equal(hl$variant, "error")
  expect_true("Y" %in% hl$pk)
  expect_true("WRES" %in% hl$pk)
})
