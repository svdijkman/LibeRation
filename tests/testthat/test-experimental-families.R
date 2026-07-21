experimental_theta <- function(value) {
  data.frame(
    THETA = seq_along(value), Value = value,
    LOWER = rep(1e-6, length(value)), UPPER = rep(100, length(value))
  )
}

experimental_data <- function() {
  data.frame(
    ID = 1, TIME = c(0, 0.5, 1, 2), EVID = c(1, 0, 0, 0),
    AMT = c(10, 0, 0, 0)
  )
}

test_that("DDE method-of-steps simulates and differentiates", {
  acknowledgement <- nm_experimental_config(TRUE, label = "test")
  expect_error(nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT"), ADVAN = 6,
    PRED = "K=THETA(1); TAU=THETA(2); S1=1",
    DES = "DADT(1)=-K*A(1)+0.05*LAG(A(1),TAU)",
    THETAS = experimental_theta(c(0.4, 0.2)),
    DDE_CONFIG = nm_dde_config(step = 0.05, minimum_delay = 0.2)
  ), "explicit acknowledgement")
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT"), ADVAN = 6,
    PRED = "K=THETA(1); TAU=THETA(2); S1=1",
    DES = "DADT(1)=-K*A(1)+0.05*LAG(A(1),TAU)",
    THETAS = experimental_theta(c(0.4, 0.2)),
    DDE_CONFIG = nm_dde_config(step = 0.05, minimum_delay = 0.2),
    EXPERIMENTAL = acknowledgement
  )
  simulation <- nm_simulate(model, experimental_data())
  derivative <- nm_prediction_derivatives(model, experimental_data())
  expect_true(all(is.finite(simulation$IPRED)))
  expect_true(all(is.finite(derivative$jacobian)))
  expect_identical(attr(simulation, "solver"), "dde")
})

test_that("block-sparse index-1 DAE remains on the AD path", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT"), ADVAN = 6,
    PRED = "K=THETA(1); S1=1", DES = "DADT(1)=-Z",
    ALG = "RES(1)=Z-K*A(1)", THETAS = experimental_theta(0.4),
    DAE_CONFIG = nm_dae_config(
      "Z", initial = 1, sparsity = matrix(TRUE, 1, 1)),
    EXPERIMENTAL = nm_experimental_config(TRUE)
  )
  simulation <- nm_simulate(model, experimental_data())
  derivative <- nm_prediction_derivatives(model, experimental_data())
  expect_true(all(is.finite(simulation$IPRED)))
  expect_true(all(is.finite(derivative$jacobian)))
})

test_that("QSP networks and DES-scoped learned components compile", {
  acknowledgement <- nm_experimental_config(TRUE)
  system <- nm_qsp_system(
    c("Drug", "Metabolite"), matrix(c(-1, 1), 2, 1),
    rates = "K*Drug", dose_species = "Drug", observation_species = "Drug"
  )
  qsp <- nm_qsp_model(
    system, INPUT = c("ID", "TIME", "EVID", "AMT"),
    PRED = "K=THETA(1);S1=1", THETAS = experimental_theta(0.4),
    EXPERIMENTAL = acknowledgement
  )
  component <- nm_component(
    "learned_loss", "dense_nn", scope = "des", inputs = "A_1",
    outputs = "LOSS", weights = list(matrix(0.4, 1, 1)), biases = list(0)
  )
  hybrid <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT"), ADVAN = 6,
    PRED = "S1=1", DES = "DADT(1)=-LOSS", THETAS = experimental_theta(1),
    COMPONENTS = component, EXPERIMENTAL = acknowledgement
  )
  expect_true(all(is.finite(nm_simulate(qsp, experimental_data())$IPRED)))
  expect_true(all(is.finite(nm_simulate(hybrid, experimental_data())$IPRED)))
  expect_identical(hybrid$DES, "DADT(1)=-LOSS")
  expect_identical(hybrid$COMPONENTS[[1]]$scope, "des")
})

test_that("factorial HMM configuration retains chain structure", {
  first <- nm_factorial_chain(
    "disease", c("low", "high"), c("I1", "I2"),
    matrix(c("T11", "T21", "T12", "T22"), 2, 2)
  )
  second <- nm_factorial_chain(
    "response", c("off", "on"), c("J1", "J2"),
    matrix(c("U11", "U21", "U12", "U22"), 2, 2)
  )
  config <- nm_factorial_hmm_config(first, second, emission = paste0("E", 1:4))
  expect_s3_class(config, "nm_factorial_hmm_config")
  expect_equal(nrow(attr(config, "factorial")$grid), 4)
  expect_match(attr(config, "generated_error"), "FHMM_T4_4")
  error <- paste(
    "I1=.6", "I2=.4", "J1=.7", "J2=.3",
    "T11=.9", "T12=.1", "T21=.2", "T22=.8",
    "U11=.8", "U12=.2", "U21=.1", "U22=.9",
    "E1=exp(-.5*(DV-0)^2)", "E2=exp(-.5*(DV-1)^2)",
    "E3=exp(-.5*(DV-2)^2)", "E4=exp(-.5*(DV-3)^2)", sep = "\n"
  )
  model <- nm_model(
    INPUT = c("ID", "TIME", "DV", "MDV"), ADVAN = 1,
    PRED = "CL=1;V=1;S1=1;F=0", ERROR = error,
    THETAS = experimental_theta(1), HMM_CONFIG = config,
    EXPERIMENTAL = nm_experimental_config(TRUE)
  )
  decoded <- nm_hmm_decode(
    model, data.frame(ID = 1, TIME = 0:2, DV = c(0.1, 1.2, 2.4), MDV = 0),
    method = "all"
  )
  expect_true(any(grepl("^FHMM_disease_FILTER_PROB_", names(decoded))))
  expect_true(any(grepl("^FHMM_response_SMOOTH_PROB_", names(decoded))))
})

test_that("switching state-space filtering exposes regime probabilities", {
  config <- nm_switching_state_space_config(
    regimes = c("stable", "flare"), initial_regime = c("RI1", "RI2"),
    regime_transition = matrix(c("RT11", "RT21", "RT12", "RT22"), 2, 2),
    states = "x", initial_mean = "M0", initial_covariance = matrix("P0", 1, 1),
    transition = list("NX1", "NX2"),
    process_covariance = list(matrix("Q1", 1, 1), matrix("Q2", 1, 1)),
    observation = c("O1", "O2"), observation_variance = c("R1", "R2"),
    by_dvid = FALSE, particles = 32
  )
  error <- paste(
    "RP=1/(1+exp(-THETA(1)))",
    "M0=0", "P0=1", "NX1=.9*STATE_x", "NX2=-.7*STATE_x",
    "Q1=.1", "Q2=.2", "O1=STATE_x", "O2=STATE_x+1", "R1=.1", "R2=.1",
    "RI1=RP", "RI2=1-RP", "RT11=RP", "RT12=1-RP", "RT21=.2", "RT22=.8",
    sep = "\n"
  )
  model <- nm_model(
    INPUT = c("ID", "TIME", "DV", "MDV"), ADVAN = 1,
    PRED = "CL=1;V=1;S1=1;F=0", ERROR = error,
    THETAS = experimental_theta(1), KALMAN_CONFIG = config,
    EXPERIMENTAL = nm_experimental_config(TRUE)
  )
  data <- data.frame(ID = 1, TIME = 0:3, DV = c(0.1, 0.2, 0.8, 0.4), MDV = 0)
  objective <- nm_objective(model, data, gradient = TRUE)
  decoded <- nm_kalman_decode(model, data)
  expect_true(is.finite(objective$value))
  expect_true(all(is.finite(objective$gradient)))
  step <- 1e-5
  numerical <- (nm_objective(model, data, theta = 1 + step, gradient = FALSE)$value -
    nm_objective(model, data, theta = 1 - step, gradient = FALSE)$value) / (2 * step)
  expect_equal(unname(objective$gradient), numerical, tolerance = 2e-4)
  expect_true(all(c("REGIME_FILTER", "REGIME_SMOOTH") %in% names(decoded)))
  expect_equal(rowSums(decoded[grep("REGIME_FILTER_PROB_", names(decoded))]),
               rep(1, nrow(data)), tolerance = 1e-12)
})
