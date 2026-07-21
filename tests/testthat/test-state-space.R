kalman_test_model <- function() {
  nm_model(
    INPUT = c("ID", "TIME", "DV", "MDV"), ADVAN = 1,
    PRED = "CL=1; V=1; S1=V; F=0",
    ERROR = paste(
      "M0 = 0", "P0 = THETA(2)",
      "A11 = exp(-THETA(1) * DT)",
      "Q11 = THETA(2) * (1 - exp(-2 * THETA(1) * DT))",
      "H1 = 1", "R1 = THETA(3) * THETA(3)", sep = "\n"
    ),
    THETAS = data.frame(
      THETA = 1:3, Value = c(0.35, 0.8, 0.25),
      LOWER = c(0.001, 0.001, 0.001), UPPER = c(5, 10, 5)
    ),
    KALMAN_CONFIG = nm_kalman_config(
      states = "deviation", initial_mean = "M0",
      initial_covariance = matrix("P0", 1, 1),
      transition = matrix("A11", 1, 1),
      process_covariance = matrix("Q11", 1, 1),
      observation = "H1", observation_variance = "R1",
      baseline = "prediction", by_dvid = FALSE
    )
  )
}

kalman_manual <- function(time, observation, theta) {
  mean <- 0
  covariance <- theta[[2]]
  nll <- 0
  filtered_mean <- filtered_variance <- predicted_mean <- numeric(length(time))
  for (row in seq_along(time)) {
    if (row > 1L) {
      dt <- time[[row]] - time[[row - 1L]]
      transition <- exp(-theta[[1]] * dt)
      process <- theta[[2]] * (1 - exp(-2 * theta[[1]] * dt))
      mean <- transition * mean
      covariance <- transition^2 * covariance + process
    }
    predicted_mean[[row]] <- mean
    innovation <- observation[[row]] - mean
    innovation_variance <- covariance + theta[[3]]^2
    nll <- nll + log(innovation_variance) + innovation^2 / innovation_variance
    gain <- covariance / innovation_variance
    mean <- mean + gain * innovation
    covariance <- (1 - gain)^2 * covariance + gain^2 * theta[[3]]^2
    filtered_mean[[row]] <- mean
    filtered_variance[[row]] <- covariance
  }
  list(
    nll = nll, predicted_mean = predicted_mean,
    filtered_mean = filtered_mean, filtered_variance = filtered_variance
  )
}

test_that("linear Gaussian state-space likelihood and gradients match reference", {
  model <- kalman_test_model()
  data <- data.frame(
    ID = "A", TIME = c(0, 0.4, 1.5, 3.1, 5),
    DV = c(0.3, 0.1, -0.4, 0.2, 0.5), MDV = 0
  )
  reference <- kalman_manual(data$TIME, data$DV, model$THETAS$Value)
  objective <- nm_objective(model, data, gradient = TRUE)
  expect_equal(objective$value, reference$nll, tolerance = 2e-10)
  expect_true(all(is.finite(objective$gradient)))
  step <- 1e-5
  numerical <- vapply(seq_len(3), function(index) {
    plus <- minus <- model$THETAS$Value
    plus[[index]] <- plus[[index]] + step
    minus[[index]] <- minus[[index]] - step
    (nm_objective(model, data, theta = plus, gradient = FALSE)$value -
       nm_objective(model, data, theta = minus, gradient = FALSE)$value) /
      (2 * step)
  }, numeric(1))
  expect_equal(unname(objective$gradient), numerical, tolerance = 2e-5)
})

test_that("Kalman decode exposes filtering and retrospective RTS smoothing", {
  model <- kalman_test_model()
  data <- data.frame(
    ID = "A", TIME = c(0, 0.4, 1.5, 3.1, 5),
    DV = c(0.3, 0.1, -0.4, 0.2, 0.5), MDV = 0
  )
  reference <- kalman_manual(data$TIME, data$DV, model$THETAS$Value)
  decoded <- nm_kalman_decode(model, data)
  expect_s3_class(decoded, "nm_kalman_decode")
  expect_equal(decoded$KF_PRED_deviation, reference$predicted_mean, tolerance = 1e-11)
  expect_equal(decoded$KF_FILTER_deviation, reference$filtered_mean, tolerance = 1e-11)
  expect_equal(decoded$KF_FILTER_SD_deviation^2, reference$filtered_variance,
               tolerance = 1e-11)
  expect_true(all(is.finite(decoded$KF_SMOOTH_deviation)))
  expect_true(all(decoded$KF_SMOOTH_SD_deviation >= 0))
  expect_equal(attr(decoded, "log_likelihood"), -0.5 * reference$nll,
               tolerance = 1e-11)
})

test_that("linear state-space trajectories simulate reproducibly", {
  model <- kalman_test_model()
  data <- data.frame(
    ID = rep(1:3, each = 5), TIME = rep(c(0, 0.4, 1.5, 3.1, 5), 3),
    DV = NA_real_, MDV = 0
  )
  first <- nm_simulate(model, data, residual = TRUE, nsim = 3, seed = 71)
  second <- nm_simulate(model, data, residual = TRUE, nsim = 3, seed = 71)
  expect_equal(first$DV, second$DV, tolerance = 0)
  expect_true(all(is.finite(first$DV)))
  expect_gt(stats::sd(first$DV), 0)
})

test_that("Kalman configuration validates dimensions and reserved likelihood outputs", {
  expect_error(
    nm_kalman_config(
      states = c("x", "y"), initial_mean = c("M1", "M2"),
      initial_covariance = matrix("P", 1, 1),
      transition = matrix("A", 2, 2), process_covariance = matrix("Q", 2, 2),
      observation = c("H1", "H2"), observation_variance = "R"
    ), "initial_covariance"
  )
  model <- kalman_test_model()
  expect_identical(model$likelihood_scale, "kalman")
  expect_error(nm_control_write(model), "KALMAN_CONFIG")
})
