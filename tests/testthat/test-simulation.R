theta_table <- function(values) {
  data.frame(THETA = seq_along(values), Value = values)
}

test_that("ADVAN1 agrees with the closed-form IV bolus solution", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT"),
    ADVAN = 1, DOSECMP = 1, OBSCMP = 1,
    PRED = "CL = THETA(1)\nV = THETA(2)\nS1 = V",
    ERROR = "Y = F",
    THETAS = theta_table(c(2, 20))
  )
  times <- c(0, 1, 5, 12)
  data <- data.frame(
    ID = 1, TIME = times, EVID = c(1, 0, 0, 0), AMT = c(100, 0, 0, 0)
  )
  result <- nm_simulate(model, data)
  expected <- 100 / 20 * exp(-(2 / 20) * times)
  expect_equal(result$IPRED, expected, tolerance = 1e-11)
  expect_equal(result$A1, expected * 20, tolerance = 1e-11)
})

test_that("selected model assignments are collected in the C++ prediction pass", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT"), OUTPUT = c("PRED", "CL", "K"),
    ADVAN = 1, PRED = "CL=THETA(1)*exp(ETA(1)); V=THETA(2); K=CL/V; S1=V",
    ERROR = "Y=F", THETAS = theta_table(c(2, 20)),
    OMEGAS = data.frame(OMEGA = 1, Value = 0.1)
  )
  data <- data.frame(
    ID = 1, TIME = c(0, 1), EVID = c(1, 0), AMT = c(100, 0)
  )
  result <- nm_simulate(model, data, eta = matrix(log(1.5), 1, 1))
  expect_equal(result$CL, rep(3, 2), tolerance = 1e-12)
  expect_equal(result$K, rep(0.15, 2), tolerance = 1e-12)
  expect_equal(result$PRED[[2L]], 5 * exp(-0.1), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(result$PRED[[2L]], result$IPRED[[2L]])))
})

test_that("ADVAN2 uses model dosing and observation compartments when CMT is absent", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT"),
    ADVAN = 2, DOSECMP = 1, OBSCMP = 2,
    PRED = "KA = THETA(1)\nCL = THETA(2)\nV = THETA(3)\nS2 = V",
    ERROR = "Y = F",
    THETAS = theta_table(c(1.5, 2, 20))
  )
  times <- c(0, 0.5, 1, 5)
  result <- nm_simulate(model, data.frame(
    ID = 1, TIME = times, EVID = c(1, 0, 0, 0), AMT = c(100, 0, 0, 0)
  ))
  ka <- 1.5
  k <- 2 / 20
  expected <- 100 * ka / (20 * (ka - k)) * (exp(-k * times) - exp(-ka * times))
  expect_equal(result$IPRED, expected, tolerance = 1e-10)
})

test_that("periodic bolus steady state is the post-dose fixed point", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "SS", "II"),
    ADVAN = 1, DOSECMP = 1, OBSCMP = 1,
    PRED = "CL = THETA(1)\nV = THETA(2)\nS1 = V",
    ERROR = "Y = F",
    THETAS = theta_table(c(2, 20))
  )
  result <- nm_simulate(model, data.frame(
    ID = 1, TIME = c(0, 6), EVID = c(1, 0), AMT = c(100, 0),
    SS = c(1, 0), II = c(12, 0)
  ))
  post <- 100 / (1 - exp(-(2 / 20) * 12))
  expect_equal(result$A1[[1]], post, tolerance = 1e-10)
  expect_equal(result$IPRED[[2]], post * exp(-(2 / 20) * 6) / 20,
               tolerance = 1e-10)
})

test_that("periodic infusion steady state composes on and off affine maps", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "RATE", "SS", "II"),
    ADVAN = 1, DOSECMP = 1, OBSCMP = 1,
    PRED = "CL = THETA(1)\nV = THETA(2)\nS1 = V",
    ERROR = "Y = F",
    THETAS = theta_table(c(2, 20))
  )
  result <- nm_simulate(model, data.frame(
    ID = 1, TIME = c(0, 5, 10, 24), EVID = c(1, 0, 0, 0),
    AMT = c(100, 0, 0, 0), RATE = c(10, 0, 0, 0),
    SS = c(1, 0, 0, 0), II = c(24, 0, 0, 0)
  ))
  k <- 2 / 20
  duration <- 10
  interval <- 24
  pre <- (10 / k) * (1 - exp(-k * duration)) *
    exp(-k * (interval - duration)) / (1 - exp(-k * interval))
  during <- pre * exp(-k * 5) + (10 / k) * (1 - exp(-k * 5))
  expect_equal(result$A1[[1]], pre, tolerance = 1e-9)
  expect_equal(result$A1[[2]], during, tolerance = 1e-9)
  expect_equal(result$A1[[4]], pre, tolerance = 1e-9)
})

test_that("ADVAN3, 4, 11, and 12 preserve non-negative linear states", {
  cases <- list(
    list(advan = 3, dose = 1, obs = 1,
         pred = "CL=THETA(1)\nVC=THETA(2)\nQ=THETA(3)\nVP=THETA(4)\nS1=VC",
         theta = c(2, 20, 1, 10)),
    list(advan = 4, dose = 1, obs = 2,
         pred = "KA=THETA(1)\nCL=THETA(2)\nVC=THETA(3)\nQ=THETA(4)\nVP=THETA(5)\nS2=VC",
         theta = c(1.5, 2, 20, 1, 10)),
    list(advan = 11, dose = 1, obs = 1,
         pred = paste("CL=THETA(1)", "VC=THETA(2)", "Q2=THETA(3)",
                      "VP1=THETA(4)", "Q3=THETA(5)", "VP2=THETA(6)",
                      "S1=VC", sep = "\n"),
         theta = c(2, 20, 1, 10, 0.5, 30)),
    list(advan = 12, dose = 1, obs = 2,
         pred = paste("KA=THETA(1)", "CL=THETA(2)", "VC=THETA(3)",
                      "Q2=THETA(4)", "VP1=THETA(5)", "Q3=THETA(6)",
                      "VP2=THETA(7)", "S2=VC", sep = "\n"),
         theta = c(1.5, 2, 20, 1, 10, 0.5, 30))
  )
  for (case in cases) {
    model <- nm_model(
      INPUT = c("ID", "TIME", "EVID", "AMT"), ADVAN = case$advan,
      DOSECMP = case$dose, OBSCMP = case$obs, PRED = case$pred,
      ERROR = "Y=F", THETAS = theta_table(case$theta)
    )
    result <- nm_simulate(model, data.frame(
      ID = 1, TIME = c(0, 1, 12, 48), EVID = c(1, 0, 0, 0),
      AMT = c(100, 0, 0, 0)
    ))
    amount_columns <- grep("^A[0-9]+$", names(result), value = TRUE)
    expect_true(all(as.matrix(result[amount_columns]) >= -1e-10),
                info = paste("ADVAN", case$advan))
    expect_true(all(is.finite(result$IPRED)), info = paste("ADVAN", case$advan))
    expect_lte(max(rowSums(result[amount_columns])), 100 + 1e-9)
  }
})

test_that("subjects have independent state and ETA rows", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT"), ADVAN = 1,
    DOSECMP = 1, OBSCMP = 1,
    PRED = "CL=THETA(1)*exp(ETA(1))\nV=THETA(2)\nS1=V",
    ERROR = "Y=F", THETAS = theta_table(c(2, 20)),
    OMEGAS = data.frame(OMEGA = 1, Value = 0.1)
  )
  result <- nm_simulate(
    model,
    data.frame(ID = c(1, 1, 2, 2), TIME = c(0, 1, 0, 1),
               EVID = c(1, 0, 1, 0), AMT = c(100, 0, 100, 0)),
    eta = matrix(c(0, log(2)), ncol = 1)
  )
  expect_equal(result$IPRED[result$ID == 1 & result$TIME == 1], 5 * exp(-0.1),
               tolerance = 1e-10)
  expect_equal(result$IPRED[result$ID == 2 & result$TIME == 1], 5 * exp(-0.2),
               tolerance = 1e-10)
})

test_that("ADVAN6 compiles DES indexing and agrees with a one-state solution", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "WT"), ADVAN = 6,
    DOSECMP = 1, OBSCMP = 1,
    PRED = "K=THETA(1)\nS1=THETA(2)",
    DES = "DADT(1) = -K * (1 + 0*T + 0*WT) * A(1)",
    COVARIATES = "WT", ERROR = "Y=F",
    THETAS = theta_table(c(0.4, 20))
  )
  times <- c(0, 0.2, 1, 4)
  result <- nm_simulate(model, data.frame(
    ID = 1, TIME = times, EVID = c(1, 0, 0, 0),
    AMT = c(100, 0, 0, 0), WT = 70
  ))
  expect_equal(result$A1, 100 * exp(-0.4 * times), tolerance = 2e-7)
  expect_equal(result$IPRED, result$A1 / 20, tolerance = 1e-12)
  expect_equal(attr(result, "solver"), "advan6-rk45")
  expect_equal(model$des_ir$output_names, "DADT_1")
})

test_that("ADVAN6 handles finite infusions inside ODE event intervals", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "RATE"), ADVAN = 6,
    DOSECMP = 1, OBSCMP = 1,
    PRED = "K=THETA(1)\nS1=1", DES = "DADT(1)=-K*A(1)",
    ERROR = "Y=F", THETAS = theta_table(0.2)
  )
  times <- c(0, 5, 10, 12)
  result <- nm_simulate(model, data.frame(
    ID = 1, TIME = times, EVID = c(1, 0, 0, 0), AMT = c(100, 0, 0, 0),
    RATE = c(10, 0, 0, 0)
  ))
  during <- (10 / 0.2) * (1 - exp(-0.2 * pmin(times, 10)))
  expected <- during * exp(-0.2 * pmax(times - 10, 0))
  expect_equal(result$A1, expected, tolerance = 2e-6)
})

test_that("ADVAN13 implicit integration resolves a stiff two-state model", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT"), ADVAN = 13,
    DOSECMP = 1, OBSCMP = 2,
    PRED = "KFAST=THETA(1)\nKSLOW=THETA(2)\nS2=1",
    DES = paste(
      "DADT(1)=-KFAST*A(1)",
      "DADT(2)=KFAST*A(1)-KSLOW*A(2)", sep = "\n"
    ),
    ERROR = "Y=F", THETAS = theta_table(c(1000, 1)),
    ODE_CONTROL = list(rtol = 2e-6, atol = 1e-9)
  )
  times <- c(0, 0.001, 0.01, 0.1)
  result <- nm_simulate(model, data.frame(
    ID = 1, TIME = times, EVID = c(1, 0, 0, 0), AMT = c(1, 0, 0, 0)
  ))
  expected_a1 <- exp(-1000 * times)
  expected_a2 <- 1000 / (1 - 1000) * (exp(-1000 * times) - exp(-times))
  expect_equal(result$A1, expected_a1, tolerance = 2e-5)
  expect_equal(result$A2, expected_a2, tolerance = 2e-5)
  expect_equal(attr(result, "solver"), "advan13-implicit")
})

test_that("ODE models reject malformed DES blocks and solve periodic steady state", {
  expect_error(
    nm_model(
      INPUT = c("ID", "TIME"), ADVAN = 6, PRED = "K=THETA(1)",
      DES = "DADT(2)=-K*A(2)", ERROR = "Y=F", THETAS = theta_table(1)
    ),
    "DADT\\(1\\) through"
  )
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "SS", "II"), ADVAN = 6,
    PRED = "K=THETA(1)", DES = "DADT(1)=-K*A(1)",
    ERROR = "Y=F", THETAS = theta_table(1)
  )
  result <- nm_simulate(model, data.frame(
    ID = 1, TIME = c(0, 1, 12), EVID = c(1, 0, 0), AMT = c(1, 0, 0),
    SS = c(1, 0, 0), II = c(12, 0, 0)
  ))
  expected_post <- 1 / (1 - exp(-12))
  expect_equal(
    result$A1,
    expected_post * exp(-c(0, 1, 12)),
    tolerance = 2e-8
  )
})

test_that("ODE periodic shooting handles overlapping steady-state infusions", {
  make_model <- function(advan, des = NULL) nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "RATE", "SS", "II"),
    ADVAN = advan, DOSECMP = 1, OBSCMP = 1,
    PRED = "K=THETA(1); V=THETA(2); S1=V", DES = des,
    ERROR = "Y=F", THETAS = theta_table(c(0.2, 20))
  )
  data <- data.frame(
    ID = 1, TIME = c(0, 2, 4, 6, 8, 12), EVID = c(1, 0, 0, 0, 0, 0),
    AMT = c(120, 0, 0, 0, 0, 0), RATE = c(10, 0, 0, 0, 0, 0),
    SS = c(1, 0, 0, 0, 0, 0), II = c(8, 0, 0, 0, 0, 0)
  )
  analytical <- nm_simulate(make_model(1), data)
  ode <- nm_simulate(make_model(6, "DADT(1)=-K*A(1)"), data)
  expect_equal(ode$A1, analytical$A1, tolerance = 3e-6)
  expect_equal(ode$IPRED, analytical$IPRED, tolerance = 2e-7)
})

test_that("IOV maps trailing ETAs to the active occasion", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "OCC"), ADVAN = 1,
    DOSECMP = 1, OBSCMP = 1,
    PRED = "CL=THETA(1)*exp(ETA(1)+ETA(2)); V=THETA(2); S1=V",
    ERROR = "Y=F+ERR(1)", THETAS = theta_table(c(2, 20)),
    OMEGAS = data.frame(OMEGA = 1:2, Value = c(0.1, 0.2)),
    SIGMAS = data.frame(SIGMA = 1, Value = 0.1), IOV = 1
  )
  data <- data.frame(
    ID = 1, TIME = c(0, 1, 24, 25), EVID = c(1, 0, 4, 0),
    AMT = c(100, 0, 100, 0), OCC = c(1, 1, 2, 2)
  )
  eta <- matrix(c(0, log(2), log(0.5)), 1)
  result <- nm_simulate(model, data, eta = eta)
  expect_equal(result$IPRED[c(2, 4)], 5 * exp(-c(0.2, 0.05)), tolerance = 1e-11)
  derivative <- nm_prediction_derivatives(model, data, eta = eta)
  expect_equal(ncol(derivative$jacobian), 2 + 3 + 1)
  expect_error(nm_simulate(model, transform(data, OCC = NULL)), "occasion column")
})

test_that("stochastic simulation draws ETAs, mixtures, residuals, and replicates reproducibly", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT"), ADVAN = 1,
    PRED = paste(
      "CL=ifelse(MIXNUM==1,THETA(1),THETA(2))*exp(ETA(1))",
      "V=THETA(3); S1=V", sep = ";"
    ), ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(THETA = 1:3, Value = c(1, 3, 20)),
    OMEGAS = data.frame(OMEGA = 1, Value = 0.1),
    SIGMAS = data.frame(SIGMA = 1, Value = 0.2),
    LIK_CONFIG = nm_lik_config(
      error = "additive", mixtures = nm_mixture(c(0.6, 0.4))
    )
  )
  data <- data.frame(
    ID = rep(1:4, each = 2), TIME = rep(c(0, 2), 4),
    EVID = rep(c(1, 0), 4), AMT = rep(c(100, 0), 4)
  )
  first <- nm_simulate(
    model, data, nsim = 3, random_effects = TRUE, residual = TRUE,
    sample_mixture = TRUE, seed = 412
  )
  second <- nm_simulate(
    model, data, nsim = 3, random_effects = TRUE, residual = TRUE,
    sample_mixture = TRUE, seed = 412
  )
  expect_equal(first, second)
  expect_equal(unique(first$SIM), 1:3)
  expect_true(all(first$MIXNUM %in% 1:2))
  expect_true(any(abs(first$ETA1) > 0))
  expect_true(all(is.finite(first$DV[first$EVID == 0])))
})

test_that("NONMEM variance and historical SD SIGMA parameterizations agree", {
  make_model <- function(value, parameterization) nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT"), ADVAN = 1,
    PRED = "CL=THETA(1); V=THETA(2); S1=V", ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(THETA = 1:2, Value = c(2, 20)),
    SIGMAS = data.frame(SIGMA = 1, Value = value),
    LIK_CONFIG = nm_lik_config(
      error = "additive", sigma_parameterization = parameterization
    )
  )
  data <- data.frame(
    ID = 1, TIME = c(0, 1, 2), EVID = c(1, 0, 0), AMT = c(100, 0, 0)
  )
  sd_result <- nm_simulate(make_model(0.2, "sd"), data, residual = TRUE, seed = 91)
  variance_result <- nm_simulate(
    make_model(0.04, "variance"), data, residual = TRUE, seed = 91
  )
  expect_equal(sd_result$DV, variance_result$DV, tolerance = 1e-14)
})

test_that("RATE=-1 and RATE=-2 use modelled infusion rate and duration", {
  data <- data.frame(
    ID = 1, TIME = c(0, 2, 6), EVID = c(1, 0, 0),
    AMT = c(100, 0, 0), RATE = c(-1, 0, 0)
  )
  rate_model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "RATE"), ADVAN = 1,
    PRED = "CL=THETA(1);V=THETA(2);R1=THETA(3);S1=V", ERROR = "Y=F",
    THETAS = theta_table(c(2, 20, 20))
  )
  explicit <- transform(data, RATE = c(20, 0, 0))
  expect_equal(nm_simulate(rate_model, data)$IPRED,
               nm_simulate(rate_model, explicit)$IPRED, tolerance = 1e-12)

  duration_model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "RATE"), ADVAN = 1,
    PRED = "CL=THETA(1);V=THETA(2);D1=THETA(3);S1=V", ERROR = "Y=F",
    THETAS = theta_table(c(2, 20, 5))
  )
  duration_data <- transform(data, RATE = c(-2, 0, 0))
  expect_equal(nm_simulate(duration_model, duration_data)$IPRED,
               nm_simulate(duration_model, explicit)$IPRED, tolerance = 1e-12)
  missing_rate_model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "RATE"), ADVAN = 1,
    PRED = "CL=THETA(1);V=THETA(2);S1=V", ERROR = "Y=F",
    THETAS = theta_table(c(2, 20))
  )
  expect_error(
    nm_simulate(missing_rate_model, data),
    "positive R1"
  )
})

test_that("simulation replicates can run in parallel workers", {
  skip_on_cran()
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT"), ADVAN = 1,
    PRED = "CL=THETA(1);V=THETA(2);S1=V", ERROR = "Y=F+ERR(1)",
    THETAS = theta_table(c(2, 20)),
    SIGMAS = data.frame(SIGMA = 1, Value = 0.1)
  )
  data <- data.frame(
    ID = 1, TIME = c(0, 1), EVID = c(1, 0), AMT = c(100, 0)
  )
  result <- nm_simulate(
    model, data, nsim = 2L, residual = TRUE, seed = 44,
    n_cores = 2L
  )
  expect_equal(sort(unique(result$SIM)), 1:2)
  expect_equal(attr(result, "parallel_cores"), 2L)
  expect_equal(nrow(result), 4L)
})
