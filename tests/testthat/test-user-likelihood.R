markov_model <- function(explicit = FALSE) {
  probability_code <- paste(
    "PBASE = 1 / (1 + exp(-THETA(1)))",
    "P01 = 1 / (1 + exp(-THETA(2)))",
    "P11 = 1 / (1 + exp(-THETA(3)))",
    "CL = 1; V = 1; S1 = V; F = PBASE",
    sep = "\n"
  )
  likelihood_code <- paste(
    "PCURRENT = ifelse(",
    "  FIRST == 1,",
    "  ifelse(DV == 1, PBASE, 1 - PBASE),",
    "  ifelse(PREV_DV == 0,",
    "    ifelse(DV == 1, P01, 1 - P01),",
    "    ifelse(DV == 1, P11, 1 - P11)",
    "  )",
    ")",
    if (explicit) "Y = PCURRENT" else "LOGLIK = log(pmax(PCURRENT, 1e-12))",
    sep = "\n"
  )
  nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "CMT", "DV", "MDV"),
    ADVAN = 1,
    PRED = probability_code,
    ERROR = likelihood_code,
    THETAS = data.frame(
      THETA = 1:3,
      Value = stats::qlogis(c(0.4, 0.8, 0.7)),
      LOWER = -10, UPPER = 10
    ),
    LIK_CONFIG = if (explicit) nm_lik_config(error = "likelihood") else NULL,
    OUTPUT = c("PBASE", "P01", "P11")
  )
}

test_that("LIK and LOGLIK select the compiled user likelihood", {
  auto <- markov_model()
  expect_identical(auto$LIK_CONFIG$error, "likelihood")
  expect_identical(auto$likelihood_output, "LOGLIK")
  expect_identical(auto$likelihood_scale, "log")
  expect_s3_class(auto$error_ir, "libertad_ir")

  explicit <- markov_model(explicit = TRUE)
  expect_identical(explicit$likelihood_output, "Y")
  expect_identical(explicit$likelihood_scale, "likelihood")

  lik <- nm_model(
    INPUT = c("ID", "TIME", "DV", "MDV"), ADVAN = 1,
    PRED = "CL=1; V=1; S1=V; F=1/(1+exp(-THETA(1)))",
    ERROR = "LIK=ifelse(DV==1,F,1-F)",
    THETAS = data.frame(THETA = 1, Value = 0, LOWER = -10, UPPER = 10)
  )
  expect_identical(lik$LIK_CONFIG$error, "likelihood")
  expect_identical(lik$likelihood_output, "LIK")
  prediction <- nm_simulate(
    lik, data.frame(ID = "A", TIME = 0, DV = 1, MDV = 0), residual = FALSE
  )
  expect_equal(prediction$IPRED, 0.5, tolerance = 1e-12)
})

test_that("Markov helpers condition on MDV baseline states", {
  model <- markov_model()
  data <- data.frame(
    ID = "subject-A", TIME = c(0, 1, 3), DV = c(0, 1, 1),
    MDV = c(1, 0, 0), EVID = 0, AMT = 0, CMT = 1
  )
  value <- nm_objective(model, data, gradient = TRUE)
  expected <- -2 * (log(0.8) + log(0.7))
  expect_equal(value$value, expected, tolerance = 1e-10)
  expect_length(value$gradient, 3L)
  expect_true(all(is.finite(value$gradient)))
  # The baseline initial-state parameter is not used when the first record is
  # MDV=1, but that record is retained as PREV_DV for the first transition.
  expect_equal(unname(value$gradient[[1L]]), 0, tolerance = 1e-10)
})

test_that("Markov time helpers follow each response sequence", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "DV", "MDV", "DVID"), ADVAN = 1,
    PRED = "CL=1; V=1; S1=V; F=0.5",
    ERROR = paste(
      "GAP = (DT + TIME - PREV_TIME) / 2",
      "PSTAY = exp(-THETA(1) * GAP)",
      "LIK = ifelse(FIRST == 1, 1, ifelse(DV == PREV_DV, PSTAY, 1 - PSTAY))",
      sep = "\n"
    ),
    THETAS = data.frame(THETA = 1, Value = 0.2, LOWER = 0.001, UPPER = 2)
  )
  data <- data.frame(
    ID = "A", TIME = c(0, 0, 2, 3, 5), DVID = c(1, 2, 1, 2, 1),
    DV = c(0, 1, 0, 1, 1), MDV = c(1, 1, 0, 0, 0)
  )
  value <- nm_objective(model, data, gradient = FALSE)$value
  expected <- -2 * (log(exp(-0.2 * 2)) + log(exp(-0.2 * 3)) +
                    log(1 - exp(-0.2 * 3)))
  expect_equal(value, expected, tolerance = 1e-10)
})

test_that("user likelihoods reject Gaussian linearization methods", {
  model <- markov_model()
  data <- data.frame(
    ID = rep(1:2, each = 3), TIME = rep(0:2, 2),
    DV = c(0, 1, 1, 1, 0, 1), MDV = rep(c(1, 0, 0), 2),
    EVID = 0, AMT = 0, CMT = 1
  )
  expect_error(
    nm_est(model, data, method = "FOCEI", maxit = 1),
    "Gaussian residual linearization"
  )
})

test_that("Gaussian residual diagnostics are omitted for user likelihoods", {
  model <- markov_model()
  data <- nm_dataset(data.frame(
    ID = 1, TIME = c(0, 1, 2), DV = c(0, 1, 1),
    MDV = c(1, 0, 0), EVID = 0, AMT = 0, CMT = 1
  ))
  fit <- structure(list(
    model = model, data = data, method = "LAPLACE",
    theta = model$THETAS$Value, sigma = numeric(), omega = numeric(),
    eta = matrix(numeric(), 1, 0)
  ), class = "nm_fit")
  gof <- nm_gof(fit)
  expect_true(all(is.na(gof$CWRES)))
  expect_match(attr(gof, "residual_note"), "not defined")
  expect_true(all(c("PBASE", "P01", "P11") %in% names(gof)))
})
