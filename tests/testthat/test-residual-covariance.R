correlated_endpoint_fixture <- function(correlation = "THETA(3)") {
  raw_correlation <- stats::qlogis(0.7) / 2 # atanh(0.4)
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "DV", "MDV", "DVID"),
    ADVAN = 1, DOSECMP = 1, OBSCMP = 1,
    PRED = "CL=THETA(1); V=THETA(2); S1=V",
    ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(
      THETA = 1:3, Value = c(2, 20, raw_correlation),
      LOWER = c(0.01, 0.1, -5), UPPER = c(20, 200, 5)
    ),
    SIGMAS = data.frame(SIGMA = 1:2, Value = c(0.5, 0.8)),
    LIK_CONFIG = nm_lik_config(
      error = "additive",
      residual_groups = nm_residual_group(
        dvid = c(1, 2),
        correlation = matrix(c("1", correlation, correlation, "1"), 2, 2),
        parameter_transform = "tanh", label = "PK/PD endpoints"
      )
    )
  )
  data <- do.call(rbind, lapply(1:2, function(id) data.frame(
    ID = id, TIME = c(0, 1, 1, 4, 4),
    EVID = c(1, 0, 0, 0, 0), AMT = c(100, 0, 0, 0, 0),
    MDV = c(1, 0, 0, 0, 0), DVID = c(1, 1, 2, 1, 2),
    DV = NA_real_
  )))
  prediction <- nm_simulate(model, data, residual = FALSE)$IPRED
  residual <- rep(c(NA, 0.2, -0.3, -0.1, 0.4), 2)
  data$DV <- prediction + residual
  data$DV[data$EVID == 1] <- NA_real_
  list(model = model, data = data)
}

test_that("estimated AR1 correlation is transformed and differentiated", {
  rho <- 0.35
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "DV", "MDV"), ADVAN = 1,
    PRED = "CL=THETA(1); V=THETA(2); S1=V", ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(
      THETA = 1:3, Value = c(2, 20, atanh(rho)),
      LOWER = c(0.01, 0.1, -5), UPPER = c(20, 200, 5)
    ),
    SIGMAS = data.frame(SIGMA = 1, Value = 0.3),
    LIK_CONFIG = nm_lik_config(
      error = "additive", sigma_corr = "ar1",
      ar1_parameter = "THETA(3)", ar1_transform = "tanh"
    )
  )
  data <- data.frame(
    ID = 1, TIME = c(0, 1, 2, 4), EVID = c(1, 0, 0, 0),
    AMT = c(100, 0, 0, 0), MDV = c(1, 0, 0, 0), DV = NA_real_
  )
  prediction <- nm_simulate(model, data, residual = FALSE)$IPRED
  data$DV <- prediction + c(NA, 0.1, -0.2, 0.15)
  objective <- nm_objective(model, data, gradient = TRUE)
  step <- 1e-5
  plus <- minus <- model$THETAS$Value
  plus[[3]] <- plus[[3]] + step
  minus[[3]] <- minus[[3]] - step
  numerical <- (
    nm_objective(model, data, theta = plus, gradient = FALSE)$value -
      nm_objective(model, data, theta = minus, gradient = FALSE)$value
  ) / (2 * step)
  expect_equal(unname(objective$gradient[["THETA_3"]]), numerical, tolerance = 2e-6)
  expect_equal(LibeRation:::.nm_ar1_rho(model), rho, tolerance = 1e-12)
})

test_that("cross-endpoint SIGMA groups use one multivariate contribution", {
  fixture <- correlated_endpoint_fixture()
  objective <- nm_objective(fixture$model, fixture$data, gradient = TRUE)
  rho <- tanh(fixture$model$THETAS$Value[[3]])
  covariance <- matrix(c(
    0.5^2, rho * 0.5 * 0.8, rho * 0.5 * 0.8, 0.8^2
  ), 2, 2)
  residuals <- list(c(0.2, -0.3), c(-0.1, 0.4))
  expected_subject <- sum(vapply(residuals, function(residual) {
    as.numeric(determinant(covariance, logarithm = TRUE)$modulus +
                 crossprod(residual, solve(covariance, residual)))
  }, numeric(1)))
  expect_equal(objective$value, 2 * expected_subject, tolerance = 2e-10)

  step <- 1e-5
  plus <- minus <- fixture$model$THETAS$Value
  plus[[3]] <- plus[[3]] + step
  minus[[3]] <- minus[[3]] - step
  numerical <- (
    nm_objective(fixture$model, fixture$data, theta = plus, gradient = FALSE)$value -
      nm_objective(fixture$model, fixture$data, theta = minus, gradient = FALSE)$value
  ) / (2 * step)
  expect_equal(unname(objective$gradient[["THETA_3"]]), numerical, tolerance = 3e-6)
})

test_that("cross-endpoint correlation is used by FO and simulation", {
  fixture <- correlated_endpoint_fixture()
  context <- LibeRation:::.nm_estimation_context(fixture$model, fixture$data)
  parameters <- list(
    theta = fixture$model$THETAS$Value,
    sigma = fixture$model$SIGMAS$Value, omega = numeric()
  )
  fo <- sum(vapply(context$subjects, function(subject) {
    subject$fo_objective(parameters$theta, parameters$sigma, parameters$omega)$value
  }, numeric(1)))
  exact <- nm_objective(fixture$model, fixture$data, gradient = FALSE)$value
  expect_equal(fo, exact, tolerance = 2e-9)

  simulated <- nm_simulate(
    fixture$model, fixture$data, residual = TRUE, nsim = 400, seed = 902
  )
  rows <- simulated$EVID == 0 & simulated$TIME == 1
  wide <- reshape(
    simulated[rows, c("ID", "SIM", "DVID", "DV", "IPRED")],
    idvar = c("ID", "SIM"), timevar = "DVID", direction = "wide"
  )
  empirical <- stats::cor(wide$DV.1 - wide$IPRED.1, wide$DV.2 - wide$IPRED.2)
  expect_equal(empirical, 0.4, tolerance = 0.10)
})

test_that("residual-group declarations reject invalid and overlapping matrices", {
  expect_error(
    nm_residual_group(c(1, 2), matrix(c(1, 0.2, 0.4, 1), 2)),
    "symmetric"
  )
  group <- nm_residual_group(c(1, 2), matrix(c(1, 0.2, 0.2, 1), 2))
  expect_error(
    nm_lik_config(residual_groups = list(group, group)),
    "only one correlated residual group"
  )
})
