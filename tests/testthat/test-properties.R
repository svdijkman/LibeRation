test_that("ADVAN1 obeys analytic superposition for random event schedules", {
  set.seed(20260723)
  for (case in seq_len(20L)) {
    clearance <- stats::runif(1L, 0.2, 8)
    volume <- stats::runif(1L, 5, 100)
    dose_times <- sort(stats::runif(sample(1:5, 1L), 0, 24))
    doses <- stats::runif(length(dose_times), 10, 500)
    observation_times <- sort(stats::runif(20L, 0, 36))
    rows <- rbind(
      data.frame(TIME = dose_times, EVID = 1L, AMT = doses, DV = NA, MDV = 1L),
      data.frame(TIME = observation_times, EVID = 0L, AMT = 0, DV = NA, MDV = 1L)
    )
    rows <- rows[order(rows$TIME, -rows$EVID), ]
    rows$ID <- 1L
    model <- nm_model(
      INPUT = c("ID", "TIME", "EVID", "AMT", "DV", "MDV"),
      ADVAN = 1, PRED = "CL=THETA(1);V=THETA(2);S1=V", ERROR = "Y=F",
      THETAS = data.frame(THETA = 1:2, Value = c(clearance, volume)),
      OMEGAS = NULL, SIGMAS = NULL
    )
    simulated <- nm_simulate(model, rows, residual = FALSE)
    observations <- simulated[simulated$EVID == 0L, ]
    expected <- vapply(observations$TIME, function(time) {
      active <- dose_times <= time
      sum(doses[active] / volume *
            exp(-clearance / volume * (time - dose_times[active])))
    }, numeric(1))
    expect_equal(observations$IPRED, expected, tolerance = 2e-10)
  }
})

test_that("parameter bounds round-trip under random finite initial values", {
  set.seed(20260724)
  for (initial in exp(stats::runif(50L, log(1e-6), log(1e6)))) {
    model <- nm_model(
      INPUT = c("ID", "TIME", "EVID", "AMT", "DV", "MDV"),
      ADVAN = 1, PRED = "CL=THETA(1);V=1;S1=V", ERROR = "Y=F",
      THETAS = data.frame(THETA = 1, Value = initial),
      OMEGAS = NULL, SIGMAS = NULL
    )
    expect_lte(model$THETAS$LOWER[[1L]], initial)
    expect_gte(model$THETAS$UPPER[[1L]], initial)
    expect_true(is.finite(model$THETAS$LOWER[[1L]]))
    expect_true(is.finite(model$THETAS$UPPER[[1L]]))
  }
})
