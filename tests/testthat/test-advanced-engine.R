nonlinear_filter_fixture <- function(filter = "ekf", particles = 512L) {
  data <- data.frame(
    ID = 1L, TIME = 0:5, EVID = 0L, AMT = 0, CMT = 1L,
    MDV = 0L, DV = c(0.2, 0.4, 0.1, 0.3, 0.2, 0.25)
  )
  config <- nm_kalman_config(
    states = "x", initial_mean = "M0",
    initial_covariance = matrix("P0", 1L), transition = "XN",
    process_covariance = matrix("Q0", 1L), observation = "HX",
    observation_variance = "R0", baseline = "zero", filter = filter,
    particles = particles, seed = 71L
  )
  model <- nm_model(
    INPUT = names(data), ADVAN = 1, PRED = "CL=1;V=1;S1=1",
    ERROR = paste(
      "M0=0;P0=1", "XN=STATE_x*exp(-THETA(1)*DT)",
      "Q0=THETA(2)*(1-exp(-2*THETA(1)*DT))", "HX=STATE_x",
      "R0=THETA(3)", sep = ";"
    ),
    THETAS = data.frame(THETA = 1:3, Value = c(0.3, 0.5, 0.2)),
    KALMAN_CONFIG = config
  )
  list(model = model, data = data)
}

test_that("EKF and UKF use differentiable nonlinear state functions", {
  ekf <- nonlinear_filter_fixture("ekf")
  ukf <- nonlinear_filter_fixture("ukf")
  ekf_value <- nm_compile(ekf$model)$objective(ekf$data)
  ukf_value <- nm_compile(ukf$model)$objective(ukf$data)
  expect_equal(ekf_value$value, ukf_value$value, tolerance = 1e-8)
  expect_equal(ekf_value$gradient, ukf_value$gradient, tolerance = 1e-6)
  decoded <- nm_kalman_decode(ukf$model, ukf$data)
  expect_identical(attr(decoded, "filter"), "ukf")
  expect_identical(attr(decoded, "smoother"), "RTS")
  expect_true(all(is.finite(decoded$KF_SMOOTH_x)))
})

test_that("particle likelihood and genealogy are seeded and reproducible", {
  fixture <- nonlinear_filter_fixture("particle", particles = 256L)
  first <- nm_compile(fixture$model)$objective(fixture$data)
  second <- nm_compile(fixture$model)$objective(fixture$data)
  expect_identical(first$value, second$value)
  expect_identical(first$gradient, second$gradient)
  expect_true(all(is.finite(first$gradient)))
  decoded <- nm_kalman_decode(fixture$model, fixture$data)
  expect_identical(attr(decoded, "filter"), "particle")
  expect_identical(attr(decoded, "smoother"), "genealogical")
  expect_true(all(is.finite(decoded$KF_SMOOTH_x)))
})

test_that("continuous-discrete SDEs support Euler and Milstein", {
  data <- data.frame(
    ID = 1L, TIME = 0:4, EVID = 0L, AMT = 0, CMT = 1L,
    MDV = 0L, DV = c(0.1, 0.3, 0.2, 0.4, 0.1)
  )
  make_model <- function(method, filter = "ukf") nm_model(
    INPUT = names(data), ADVAN = 1, PRED = "CL=1;V=1;S1=1",
    ERROR = paste(
      "M0=0;P0=1", "DRIFT=-THETA(1)*STATE_x", "G0=THETA(2)",
      "HX=STATE_x", "R0=THETA(3)", sep = ";"
    ),
    THETAS = data.frame(THETA = 1:3, Value = c(0.2, 0.3, 0.1)),
    KALMAN_CONFIG = nm_sde_config(
      states = "x", initial_mean = "M0",
      initial_covariance = matrix("P0", 1L), drift = "DRIFT",
      diffusion = matrix("G0", 1L), observation = "HX",
      observation_variance = "R0", baseline = "zero", filter = filter,
      method = method, substeps = 4L, particles = 64L, seed = 19L
    )
  )
  objective <- nm_compile(make_model("euler"))$objective(data)
  expect_true(is.finite(objective$value))
  expect_true(all(is.finite(objective$gradient)))
  first <- nm_simulate(make_model("milstein", "particle"), data,
                       residual = TRUE, seed = 90L)
  second <- nm_simulate(make_model("milstein", "particle"), data,
                        residual = TRUE, seed = 90L)
  expect_equal(first$DV, second$DV)
  expect_false(isTRUE(all.equal(first$DV, data$DV)))
})

test_that("nested and crossed random-effect blocks map shared ETAs exactly", {
  data <- expand.grid(ID = 1:4, TIME = 1:2)
  data <- data[order(data$ID, data$TIME), ]
  data$SITE <- rep(c("A", "B"), each = 4L)
  data$READER <- rep(c("R1", "R2", "R2", "R1"), each = 2L)
  data$EVID <- 0L; data$AMT <- 0; data$CMT <- 1L; data$MDV <- 0L
  data$DV <- 10 + data$ID / 10
  design <- nm_re_config(
    nm_re_block("site", "SITE", 1L),
    nm_re_block("patient", "ID", 2L),
    nm_re_block("reader", "READER", 3L)
  )
  model <- nm_model(
    INPUT = names(data), ADVAN = 1,
    PRED = "CL=THETA(1)*exp(ETA(1)+ETA(2)+ETA(3));V=1;S1=1;F=CL",
    ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(THETA = 1L, Value = 10),
    OMEGAS = data.frame(OMEGA = 1:3, Value = c(0.1, 0.2, 0.15)),
    SIGMAS = data.frame(SIGMA = 1L, Value = 0.1), RE_CONFIG = design
  )
  compiled_data <- LibeRation:::.nm_engine_data(model, data)
  expect_equal(length(unique(compiled_data$.ID_INDEX)), 1L)
  expect_equal(length(unique(compiled_data$.STRUCT_ID_INDEX)), 4L)
  expect_equal(LibeRation:::.nm_eta_columns(model, compiled_data), 8L)
  eta <- matrix(seq(0.01, 0.08, by = 0.01), 1L)
  objective <- nm_compile(model)$objective(compiled_data, eta = eta)
  expect_true(is.finite(objective$value))
  expect_true(all(is.finite(objective$gradient)))
})

test_that("ARMA declarations generate a compiled temporal likelihood", {
  data <- data.frame(
    ID = 1L, TIME = 0:5, EVID = 0L, AMT = 0, CMT = 1L,
    MDV = 0L, DV = c(1, 0.4, -0.1, 0.3, 0.2, -0.2)
  )
  model <- nm_model(
    INPUT = names(data), ADVAN = 1, PRED = "CL=1;V=1;S1=1;F=0",
    THETAS = data.frame(THETA = 1:2, Value = c(0.4, 0.2)),
    SIGMAS = data.frame(SIGMA = 1L, Value = 0.5),
    KALMAN_CONFIG = nm_arma_config(
      ar = "THETA(1)", ma = "THETA(2)",
      innovation_variance = "SIGMA(1)", initial_variance = "10"
    )
  )
  objective <- nm_compile(model)$objective(data)
  expect_true(is.finite(objective$value))
  expect_true(all(is.finite(objective$gradient)))
  expect_match(model$ERROR, "ARMA_AR1", fixed = TRUE)
})

test_that("HSMM declarations retain original states after exact expansion", {
  data <- data.frame(
    ID = 1L, TIME = 0:5, EVID = 0L, AMT = 0, CMT = 1L,
    MDV = 0L, DV = c(0, 0, 1, 1, 1, 0)
  )
  config <- nm_hsmm_config(
    states = c("controlled", "active"),
    initial = c("I1", "I2"),
    transition = matrix(c("T11", "T12", "T21", "T22"), 2, 2, byrow = TRUE),
    dwell = matrix(c("D11", "D12", "D13", "D21", "D22", "D23"),
                   2, 3, byrow = TRUE),
    emission = c("E1", "E2"), by_dvid = FALSE
  )
  model <- nm_model(
    INPUT = names(data), ADVAN = 1, PRED = "CL=1;V=1;S1=1;F=0",
    THETAS = data.frame(THETA = 1L, Value = 1),
    HMM_CONFIG = config,
    ERROR = paste(
      "I1=.8;I2=.2", "T11=.85;T12=.15;T21=.25;T22=.75",
      "D11=.6;D12=.3;D13=.1;D21=.2;D22=.5;D23=.3",
      "E1=ifelse(DV==0,.9,.1);E2=ifelse(DV==0,.15,.85)", sep = ";"
    )
  )
  objective <- nm_compile(model)$objective(data)
  expect_true(is.finite(objective$value))
  expect_true(all(is.finite(objective$gradient)))
  decoded <- nm_hmm_decode(model, data, method = "all")
  expect_identical(attr(decoded, "states"), c("controlled", "active"))
  expect_true(all(decoded$HMM_VITERBI_STATE %in% c("controlled", "active")))
  expect_false(any(grepl("@", decoded$HMM_VITERBI_STATE, fixed = TRUE)))
})
