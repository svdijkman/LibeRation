objective_fixture <- function(config = NULL, error = "Y=F*(1+ERR(1))+ERR(2)") {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "DV", "MDV"),
    ADVAN = 1, DOSECMP = 1, OBSCMP = 1,
    PRED = "CL=THETA(1)*exp(ETA(1)); V=THETA(2); S1=V",
    ERROR = error,
    THETAS = data.frame(THETA = 1:2, Value = c(2, 20)),
    OMEGAS = data.frame(OMEGA = 1, Value = 0.1),
    SIGMAS = if (grepl("ERR\\(2\\)", error)) {
      data.frame(SIGMA = 1:2, Value = c(0.1, 0.2))
    } else data.frame(SIGMA = 1, Value = 0.2),
    LIK_CONFIG = config
  )
  data <- data.frame(
    ID = rep(1:2, each = 4), TIME = rep(c(0, 1, 5, 12), 2),
    EVID = rep(c(1, 0, 0, 0), 2), AMT = rep(c(100, 0, 0, 0), 2),
    MDV = rep(c(1, 0, 0, 0), 2),
    DV = c(NA, 4.6, 3.1, 1.5, NA, 4.2, 2.8, 1.2)
  )
  list(model = model, data = data, eta = matrix(c(0.1, -0.1), ncol = 1))
}

objective_central_difference <- function(fn, at, index, step = 1e-5) {
  plus <- minus <- at
  plus[[index]] <- plus[[index]] + step
  minus[[index]] <- minus[[index]] - step
  (fn(plus) - fn(minus)) / (2 * step)
}

test_that("joint likelihood gradient and Hessian come from the full AD tape", {
  fixture <- objective_fixture()
  engine <- nm_compile(fixture$model)
  tape <- engine$objective_tape(fixture$data, eta = fixture$eta)
  exact <- .liberation_objective_tape_eval(tape$pointer, tape$point, TRUE, TRUE)
  fn <- function(point) {
    .liberation_objective_tape_eval(tape$pointer, point, FALSE, FALSE)$value
  }
  gr <- function(point) {
    .liberation_objective_tape_eval(tape$pointer, point, TRUE, FALSE)$gradient
  }
  numerical_gradient <- vapply(seq_along(tape$point), function(index) {
    objective_central_difference(fn, tape$point, index, step = 2e-6)
  }, numeric(1))
  numerical_hessian <- vapply(seq_along(tape$point), function(index) {
    objective_central_difference(gr, tape$point, index, step = 3e-5)
  }, numeric(length(tape$point)))
  expect_equal(unname(exact$gradient), numerical_gradient, tolerance = 3e-6)
  expect_equal(unname(exact$hessian), unname(numerical_hessian), tolerance = 3e-5)
  expect_equal(exact$hessian, t(exact$hessian), tolerance = 2e-10)
})

test_that("native conditional modes and batched objective values match scalar evaluation", {
  fixture <- estimation_fixture()
  context <- LibeRation:::.nm_estimation_context(fixture$model, fixture$data)
  parameters <- list(
    theta = fixture$model$THETAS$Value,
    sigma = fixture$model$SIGMAS$Value,
    omega = fixture$model$OMEGAS$Value
  )
  evaluator <- context$subjects[[1L]]
  samples <- matrix(c(-0.2, 0, 0.15, 0.3), ncol = context$n_eta)
  batched <- evaluator$objective_eta_values(
    parameters$theta, samples, parameters$sigma, parameters$omega
  )
  scalar <- vapply(seq_len(nrow(samples)), function(index) {
    evaluator$objective(
      parameters$theta, samples[index, ], parameters$sigma, parameters$omega
    )$value
  }, numeric(1))
  expect_equal(batched, scalar, tolerance = 1e-12)

  eta <- matrix(c(-0.1, 0.05, 0.2), ncol = context$n_eta)
  collection <- LibeRation:::.nm_objective_collection(
    context$subjects, parameters, eta
  )
  individual <- vapply(seq_along(context$subjects), function(subject) {
    context$subjects[[subject]]$objective(
      parameters$theta, eta[subject, ], parameters$sigma, parameters$omega
    )$value
  }, numeric(1))
  expect_equal(collection, individual, tolerance = 1e-12)

  mode <- evaluator$eta_mode(
    parameters$theta, parameters$sigma, parameters$omega,
    exact_hessian = FALSE
  )
  expect_identical(mode$backend, "cpp")
  expect_lt(max(abs(mode$gradient)), 1e-5)
})

test_that("batched gradients and Hessian subsets match scalar tape derivatives", {
  fixture <- estimation_fixture()
  context <- LibeRation:::.nm_estimation_context(fixture$model, fixture$data)
  parameters <- list(
    theta = fixture$model$THETAS$Value,
    sigma = fixture$model$SIGMAS$Value,
    omega = fixture$model$OMEGAS$Value
  )
  eta <- matrix(c(-0.1, 0.05, 0.2), ncol = context$n_eta)
  batched <- LibeRation:::.nm_objective_collection_gradient(
    context$subjects, parameters, eta
  )
  scalar <- do.call(rbind, lapply(seq_along(context$subjects), function(subject) {
    context$subjects[[subject]]$objective(
      parameters$theta, eta[subject, ], parameters$sigma, parameters$omega,
      gradient = TRUE
    )$gradient
  }))
  expect_equal(unname(batched), unname(scalar), tolerance = 1e-12)

  evaluator <- context$subjects[[1L]]
  full <- evaluator$objective(
    parameters$theta, eta[1L, ], parameters$sigma, parameters$omega,
    hessian = TRUE
  )$hessian
  rows <- c(2L, 3L)
  columns <- c(1L, 3L, 5L)
  subset <- evaluator$objective_hessian_subset(
    parameters$theta, eta[1L, ], parameters$sigma, parameters$omega,
    rows = rows, columns = columns
  )
  expect_equal(unname(subset), unname(full[rows, columns]), tolerance = 1e-12)
})

test_that("M3 BLQ likelihood is differentiated in C++", {
  config <- nm_lik_config(error = "additive", blq_method = "m3", lloq = 2)
  fixture <- objective_fixture(config, error = "Y=F+ERR(1)")
  fixture$data$BLQ <- as.integer(!is.na(fixture$data$DV) & fixture$data$DV < 2)
  objective <- nm_objective(fixture$model, fixture$data, eta = fixture$eta)
  expect_true(is.finite(objective$value))
  expect_true(all(is.finite(objective$gradient)))
  # The censored contribution is a probability, so changing the fixed-effect
  # prediction changes the exact objective and its THETA derivative.
  expect_false(isTRUE(all.equal(unname(objective$gradient[["THETA_1"]]), 0)))
})

test_that("AR1 likelihood follows innovation form within each subject", {
  config <- nm_lik_config(error = "additive", sigma_corr = "ar1", ar1_rho = 0.35)
  fixture <- objective_fixture(config, error = "Y=F+ERR(1)")
  result <- nm_objective(fixture$model, fixture$data, eta = fixture$eta)
  prediction <- nm_simulate(fixture$model, fixture$data, eta = fixture$eta)$IPRED
  observed <- fixture$data$EVID == 0 & fixture$data$MDV == 0
  variance <- fixture$model$SIGMAS$Value[[1]]^2
  manual <- 0
  for (id in unique(fixture$data$ID)) {
    rows <- which(fixture$data$ID == id & observed)
    residual <- fixture$data$DV[rows] - prediction[rows]
    manual <- manual + log(variance) + residual[[1]]^2 / variance
    innovation_variance <- variance * (1 - 0.35^2)
    for (i in seq.int(2, length(residual))) {
      innovation <- residual[[i]] - 0.35 * residual[[i - 1L]]
      manual <- manual + log(innovation_variance) + innovation^2 / innovation_variance
    }
  }
  eta <- as.vector(fixture$eta)
  omega <- fixture$model$OMEGAS$Value[[1]]
  manual <- manual + sum(log(omega) + eta^2 / omega)
  expect_equal(result$value, manual, tolerance = 1e-10)
})

test_that("heteroscedastic AR1 correlates standardized residuals", {
  rho <- 0.45
  config <- nm_lik_config(
    error = "proportional", sigma_corr = "ar1", ar1_rho = rho
  )
  fixture <- objective_fixture(config, error = "Y=F*(1+ERR(1))")
  result <- nm_objective(fixture$model, fixture$data, eta = fixture$eta)
  prediction <- nm_simulate(fixture$model, fixture$data, eta = fixture$eta)$IPRED
  observed <- fixture$data$EVID == 0 & fixture$data$MDV == 0
  manual <- 0
  for (id in unique(fixture$data$ID)) {
    rows <- which(fixture$data$ID == id & observed)
    variance <- (fixture$model$SIGMAS$Value[[1]] * prediction[rows])^2
    z <- (fixture$data$DV[rows] - prediction[rows]) / sqrt(variance)
    manual <- manual + log(variance[[1]]) + z[[1]]^2
    for (i in seq.int(2, length(rows))) {
      innovation <- z[[i]] - rho * z[[i - 1L]]
      manual <- manual + log(variance[[i]]) + log(1 - rho^2) +
        innovation^2 / (1 - rho^2)
    }
  }
  eta <- as.vector(fixture$eta)
  omega <- fixture$model$OMEGAS$Value[[1]]
  manual <- manual + sum(log(omega) + eta^2 / omega)
  expect_equal(result$value, manual, tolerance = 1e-10)
})

test_that("full OMEGA covariance uses an exact Cholesky prior", {
  omega <- data.frame(
    OMEGA = 1:3, ROW = c(1, 2, 2), COL = c(1, 1, 2),
    Value = c(0.09, 0.02, 0.16)
  )
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "DV", "MDV"),
    ADVAN = 1, PRED = "CL=THETA(1)*exp(ETA(1)); V=THETA(2)*exp(ETA(2)); S1=V",
    ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(THETA = 1:2, Value = c(2, 20)),
    OMEGAS = omega, SIGMAS = data.frame(SIGMA = 1, Value = 0.2)
  )
  data <- data.frame(ID = 1, TIME = 0, EVID = 1, AMT = 100, DV = NA, MDV = 1)
  eta <- matrix(c(0.1, -0.2), 1)
  result <- nm_objective(model, data, eta = eta)
  covariance <- matrix(c(0.09, 0.02, 0.02, 0.16), 2)
  expected <- as.numeric(determinant(covariance, logarithm = TRUE)$modulus +
                           eta %*% solve(covariance, t(eta)))
  expect_equal(result$value, expected, tolerance = 1e-12)
  expect_equal(model$n_eta, 2)
})

test_that("IOV prior repeats the trailing OMEGA block by occasion", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "DV", "MDV", "OCC"),
    ADVAN = 1, PRED = "CL=THETA(1)*exp(ETA(1)+ETA(2)); V=THETA(2); S1=V",
    ERROR = "Y=F+ERR(1)", THETAS = data.frame(THETA = 1:2, Value = c(2, 20)),
    OMEGAS = data.frame(OMEGA = 1:2, Value = c(0.1, 0.2)),
    SIGMAS = data.frame(SIGMA = 1, Value = 0.1), IOV = 1
  )
  data <- data.frame(
    ID = 1, TIME = c(0, 24), EVID = c(1, 4), AMT = 100,
    DV = NA, MDV = 1, OCC = c(1, 2)
  )
  eta <- matrix(c(0.1, 0.2, -0.3), 1)
  value <- nm_objective(model, data, eta = eta, gradient = FALSE)$value
  expected <- log(0.1) + 0.1^2 / 0.1 +
    2 * log(0.2) + (0.2^2 + (-0.3)^2) / 0.2
  expect_equal(value, expected, tolerance = 1e-12)
})

test_that("finite mixtures combine subject likelihoods by log-sum-exp", {
  mixture <- nm_mixture(c(0.7, 0.3), c("slow", "fast"))
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "DV", "MDV"),
    ADVAN = 1,
    PRED = paste(
      "CL=ifelse(MIXNUM==1,THETA(1),THETA(2))",
      "V=THETA(3)", "S1=V", sep = ";"
    ),
    ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(THETA = 1:3, Value = c(1, 4, 20)),
    SIGMAS = data.frame(SIGMA = 1, Value = 0.2),
    LIK_CONFIG = nm_lik_config(error = "additive", mixtures = mixture)
  )
  data <- data.frame(
    ID = 1, TIME = c(0, 2, 8), EVID = c(1, 0, 0), AMT = c(100, 0, 0),
    MDV = c(1, 0, 0), DV = c(NA, 4.5, 3.4)
  )
  value <- nm_objective(model, data)$value
  component_nll <- vapply(1:2, function(component) {
    component_data <- transform(data, MIXNUM = component)
    plain <- model
    plain$LIK_CONFIG$mixtures <- NULL
    # Recompile so MIXNUM remains a fixed dataset input while mixture
    # integration itself is disabled.
    plain$pred_ir <- LibeRtAD::ad_ir(
      plain$PRED, inputs = c("THETA_1", "THETA_2", "THETA_3", "MIXNUM")
    )
    nm_objective(plain, component_data)$value
  }, numeric(1))
  expected <- -2 * log(sum(c(0.7, 0.3) * exp(-0.5 * component_nll)))
  expect_equal(value, expected, tolerance = 1e-11)
  expect_true(all(is.finite(nm_objective(model, data)$gradient)))
  posterior <- nm_mixture_posterior(model, data)
  expect_equal(rowSums(posterior[c("P_slow", "P_fast")]), 1, tolerance = 1e-14)
  expect_equal(posterior$MIXTURE, c("slow"))
})
