test_that("prediction tape profiling separates recording and reuse telemetry", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT"),
    ADVAN = 6, DOSECMP = 1, OBSCMP = 1,
    PRED = "K=THETA(1); S1=THETA(2)",
    DES = "DADT(1)=-K*A(1)",
    ERROR = "Y=F",
    THETAS = data.frame(
      THETA = 1:2, Value = c(0.4, 20), FIX = TRUE
    )
  )
  data <- data.frame(
    ID = 1, TIME = c(0, 0.2, 1, 4),
    EVID = c(1, 0, 0, 0), AMT = c(100, 0, 0, 0)
  )

  profile <- nm_tape_profile(
    model, data, repetitions = 2L, jacobian = TRUE
  )
  expect_s3_class(profile, "data.frame")
  expect_equal(profile$advan, 6)
  expect_gt(profile$operations, 0)
  expect_gt(profile$resident_bytes_proxy, 0)
  expect_true(is.finite(profile$value_microseconds))
  expect_true(is.finite(profile$jacobian_microseconds))
  expect_match(profile$derivative_strategy, "forward|reverse|subgraph")
})

test_that("ODE checkpoint assessment remains measured and opt-in", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT"),
    ADVAN = 6, DOSECMP = 1, OBSCMP = 1,
    PRED = "K=THETA(1); S1=THETA(2)",
    DES = "DADT(1)=-K*A(1)",
    ERROR = "Y=F",
    THETAS = data.frame(
      THETA = 1:2, Value = c(0.4, 20), FIX = TRUE
    )
  )
  data <- data.frame(
    ID = 1, TIME = c(0, 1), EVID = c(1, 0), AMT = c(100, 0)
  )

  benchmark <- nm_ode_tape_benchmark(
    model, data, tolerances = 1e-5, repetitions = 1L,
    checkpoint_repetitions = 4L, checkpoint_evaluations = 3L
  )
  expect_s3_class(benchmark, "nm_ode_tape_benchmark")
  expect_equal(nrow(benchmark$profiles), 1)
  expect_equal(
    benchmark$assessment$kernel,
    c("ADVAN1 interval", "2x2 matrix state update")
  )
  expect_true(all(benchmark$assessment$exact))
  expect_true(all(benchmark$assessment$nested_ad_safe))
  expect_type(benchmark$assessment$production_candidate, "logical")
})
