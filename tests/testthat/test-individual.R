test_that("individual fitting estimates static and custom-prior ETA states", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "CMT", "DV", "MDV"),
    ADVAN = 2, TRANS = 2,
    PRED = "CL=THETA(1)*exp(ETA(1)); V=THETA(2); KA=THETA(3); S2=V",
    ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(THETA = 1:3, Value = c(4, 40, 1.2)),
    OMEGAS = data.frame(OMEGA = 1, Value = 0.25),
    SIGMAS = data.frame(SIGMA = 1, Value = 0.1)
  )
  events <- data.frame(
    ID = 1, TIME = c(0, 2, 6, 12), EVID = c(1, 0, 0, 0),
    AMT = c(200, 0, 0, 0), CMT = c(1, 2, 2, 2),
    DV = c(NA, 2.6, 1.7, 0.8), MDV = c(1, 0, 0, 0)
  )
  fit <- nm_individual_fit(model, events)
  expect_s3_class(fit, "nm_individual_fit")
  expect_length(fit$eta, 1L)
  expect_true(is.finite(fit$objective))
  expect_equal(dim(fit$eta_covariance), c(1L, 1L))
  expect_true(all(is.finite(fit$predictions$IPRED)))

  updated <- nm_individual_fit(
    model, events, prior_mean = fit$eta,
    prior_covariance = fit$eta_covariance + 0.05
  )
  expect_s3_class(updated, "nm_individual_fit")
  expect_true(isTRUE(updated$diagnostics$custom_prior))
})

test_that("individual fitting expands all-IOV states without resetting amounts", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "CMT", "DV", "MDV", "OCC"),
    ADVAN = 1,
    PRED = "CL=THETA(1)*exp(ETA(1)); V=THETA(2); S1=V",
    ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(THETA = 1:2, Value = c(3, 30)),
    OMEGAS = data.frame(OMEGA = 1, Value = 0.2),
    SIGMAS = data.frame(SIGMA = 1, Value = 0.1), IOV = 1L
  )
  events <- data.frame(
    ID = 1, TIME = c(0, 2, 8, 12, 14, 20),
    EVID = c(1, 0, 0, 1, 0, 0), AMT = c(100, 0, 0, 100, 0, 0),
    CMT = 1, DV = c(NA, 2.7, 1.4, NA, 2.2, 1.0),
    MDV = c(1, 0, 0, 1, 0, 0), OCC = c(1, 1, 1, 2, 2, 2)
  )
  prior <- matrix(c(0.2, 0.16, 0.16, 0.24), 2, 2)
  fit <- nm_individual_fit(
    model, events, prior_mean = c(0, 0), prior_covariance = prior
  )
  expect_named(fit$eta, c("ETA1_OCC1", "ETA1_OCC2"))
  expect_equal(dim(fit$eta_covariance), c(2L, 2L))
})
