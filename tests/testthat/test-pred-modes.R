prediction_mode_data <- function() {
  data.frame(
    ID = 1L, TIME = c(0, 1), EVID = c(1L, 0L), AMT = c(100, 0),
    RATE = 0, II = 0, SS = 0L, CMT = 1L, DV = c(NA, 0),
    MDV = c(1L, 0L), WT = 70
  )
}

test_that("direct PRED is row-wise and bypasses ADVAN dosing", {
  model <- nm_model(
    INPUT = names(prediction_mode_data()), PRED_MODE = "pred",
    PRED_SOURCE = "F=THETA(1)*TIME+WT/100",
    THETAS = data.frame(THETA = 1, Value = 2)
  )
  simulated <- nm_simulate(model, prediction_mode_data())
  expect_identical(model$SOLVER, "direct")
  expect_equal(simulated$IPRED, c(0.7, 2.7), tolerance = 1e-12)
  derivatives <- nm_prediction_derivatives(model, prediction_mode_data())
  expect_equal(
    drop(derivatives$jacobian[, match("THETA_1", derivatives$domain)]),
    c(0, 1), tolerance = 1e-12
  )
})

test_that("combined mode differentiably transforms ADVAN predictions and states", {
  model <- nm_model(
    INPUT = names(prediction_mode_data()), OUTPUT = "RAW", ADVAN = 1,
    PRED_MODE = "pk_pred",
    PK_SOURCE = "CL=THETA(1);V=THETA(2);S1=V",
    PRED_SOURCE = "RAW=F_ADVAN;F=RAW*(WT/70)+A(1)/100",
    THETAS = data.frame(THETA = 1:2, Value = c(2, 10))
  )
  simulated <- nm_simulate(model, prediction_mode_data())
  expected_state <- c(100, 100 * exp(-2 / 10))
  expect_equal(simulated$RAW, expected_state / 10, tolerance = 1e-10)
  expect_equal(simulated$IPRED, expected_state / 10 + expected_state / 100,
               tolerance = 1e-10)
  derivatives <- nm_prediction_derivatives(model, prediction_mode_data())
  expect_true(all(is.finite(derivatives$jacobian)))
})

test_that("combined mode validates its explicit final prediction contract", {
  arguments <- list(
    INPUT = c("ID", "TIME"), ADVAN = 1, PRED_MODE = "pk_pred",
    PK_SOURCE = "CL=THETA(1);V=10;S1=V",
    THETAS = data.frame(THETA = 1, Value = 2)
  )
  expect_error(
    do.call(nm_model, c(arguments, list(PRED_SOURCE = "RAW=F_ADVAN"))),
    "assign `F`"
  )
  expect_error(
    do.call(nm_model, c(arguments, list(PRED_SOURCE = "CL=1;F=F_ADVAN"))),
    "cannot overwrite"
  )
})
