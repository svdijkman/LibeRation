test_that("nm_validate_model accepts THEO-style ADVAN 4 TRANS 4 PRED", {
  sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
  out <- nm_validate_model(sim$model, sim$data, stop_on_error = FALSE)
  expect_true(out$ok)
  expect_true(all(c("CL", "VC", "VP", "Q2", "KA") %in% out$pred_symbols))
})

test_that("ADVAN 4 TRANS 4 NLL depends on VP and Q2", {
  sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
  model <- sim$model
  th <- model$THETAS$Value
  om <- model$OMEGAS$Value
  sg <- model$SIGMAS$Value
  nll <- function(x) {
    LibeRation:::.nm_nll_internal(
      model, sim$data, x, om, sg,
      eta = NULL, include_omega_prior = FALSE, pk_engine = "cpp"
    )
  }
  f0 <- nll(th)
  tp <- th
  tp[3] <- tp[3] * 1.05
  tp2 <- th
  tp2[4] <- tp2[4] * 1.05
  expect_true(abs(nll(tp) - f0) > 1e-6)
  expect_true(abs(nll(tp2) - f0) > 1e-6)
})
