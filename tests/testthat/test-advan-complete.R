advan_event_data <- function(times = c(0, 1, 4)) {
  data.frame(
    ID = 1, TIME = times, EVID = c(1L, rep(0L, length(times) - 1L)),
    AMT = c(100, rep(0, length(times) - 1L)),
    RATE = 0, II = 0, SS = 0L, CMT = 1L,
    DV = NA_real_, MDV = c(1L, rep(0L, length(times) - 1L))
  )
}

test_that("ADVAN5 and ADVAN7 use the arbitrary linear matrix engine", {
  graph <- nm_matrix_model(
    data.frame(id = 1:2, name = c("CENTRAL", "PERIPHERAL")),
    data.frame(
      from = c(1L, 1L, 2L), to = c(0L, 2L, 1L),
      type = "rate", parameter = c("K10", "K12", "K21")
    )
  )
  make_model <- function(advan) nm_model(
    INPUT = names(advan_event_data()), ADVAN = advan, TRANS = 1,
    PRED = paste(
      "K10=THETA(1)", "K12=THETA(2)", "K21=THETA(3)",
      "V=THETA(4)", "S1=V", sep = "\n"
    ),
    ERROR = "Y=F", THETAS = data.frame(THETA = 1:4, Value = c(.2, .1, .05, 20)),
    GRAPH = graph
  )
  five <- nm_simulate(make_model(5), advan_event_data())
  seven <- nm_simulate(make_model(7), advan_event_data())
  expect_equal(five$IPRED, seven$IPRED, tolerance = 1e-12)
  expect_true(all(five$A1 >= 0) && all(five$A2 >= 0))
  expect_match(attr(five, "solver"), "matrix", ignore.case = TRUE)
  derivatives <- nm_compile(make_model(5))$prediction_derivatives(
    advan_event_data()
  )
  expect_true(all(is.finite(derivatives$jacobian)))
  expect_equal(derivatives$propagation_kernel, "general-matrix-exponential")
})

test_that("ADVAN8, ADVAN9, and ADVAN14 expose stiff general-model paths", {
  make_model <- function(advan) nm_model(
    INPUT = names(advan_event_data()), ADVAN = advan, TRANS = 1,
    DOSECMP = 1, OBSCMP = 2,
    PRED = "KFAST=THETA(1);KSLOW=THETA(2);S2=1",
    DES = "DADT(1)=-KFAST*A(1)\nDADT(2)=KFAST*A(1)-KSLOW*A(2)",
    ERROR = "Y=F", THETAS = data.frame(THETA = 1:2, Value = c(1000, 1)),
    ODE_CONTROL = list(rtol = 1e-7, atol = 1e-10)
  )
  reference <- nm_simulate(make_model(13), advan_event_data(c(0, .01, .1, 1)))
  for (advan in c(8L, 9L, 14L)) {
    model <- make_model(advan)
    data <- advan_event_data(c(0, .01, .1, 1))
    result <- nm_simulate(model, data)
    expect_equal(result$IPRED, reference$IPRED, tolerance = 1e-7)
    expect_match(attr(result, "solver"), paste0("advan", advan))
    derivatives <- nm_compile(model)$prediction_derivatives(data)
    expect_true(all(is.finite(derivatives$jacobian)))
    expect_match(derivatives$propagation_kernel, paste0("advan", advan))
  }
})

test_that("ADVAN9 maps equilibrium constraints onto the differentiable DAE path", {
  model <- nm_model(
    INPUT = names(advan_event_data()), ADVAN = 9, TRANS = 1,
    PRED = "K=THETA(1);BIND=THETA(2);V=THETA(3);S1=V",
    DES = "DADT(1)=-K*FREE",
    ALG = "RES(1)=FREE-A(1)/(1+BIND)",
    DAE_CONFIG = nm_dae_config("FREE", initial = 50),
    EXPERIMENTAL = nm_experimental_config(TRUE, "ADVAN9 equilibrium validation"),
    ERROR = "Y=F", THETAS = data.frame(THETA = 1:3, Value = c(.4, 1, 20))
  )
  result <- nm_simulate(model, advan_event_data())
  expected_amount <- 100 * exp(-.4 / 2 * advan_event_data()$TIME)
  expect_equal(result$A1, expected_amount, tolerance = 2e-5)
  derivatives <- nm_compile(model)$prediction_derivatives(advan_event_data())
  expect_true(all(is.finite(derivatives$jacobian)))
  expect_match(derivatives$propagation_kernel, "dae-advan9")
})

test_that("ADVAN15-18 expose differentiable DAE, DDE, and combined DDAE paths", {
  data <- advan_event_data(c(0, .5, 1, 2, 3, 4))
  fifteen <- nm_advan_template(15)
  result15 <- nm_simulate(fifteen, data)
  expect_true(all(is.finite(result15$IPRED)))
  expect_identical(attr(result15, "solver"), "dae")
  derivative15 <- nm_compile(fifteen)$prediction_derivatives(data)
  expect_true(all(is.finite(derivative15$jacobian)))
  expect_match(derivative15$propagation_kernel, "dae-advan15")

  for (advan in c(16L, 18L)) {
    model <- nm_advan_template(advan)
    result <- nm_simulate(model, data)
    early <- data$TIME <= 2
    expect_equal(
      result$A1[early], 100 * exp(-.4 * data$TIME[early]),
      tolerance = 2e-6
    )
    derivatives <- nm_compile(model)$prediction_derivatives(data)
    expect_true(all(is.finite(derivatives$jacobian)))
    expect_match(derivatives$propagation_kernel, paste0("dde-advan", advan))
    if (advan == 16L) {
      expect_match(derivatives$propagation_kernel, "radau-iia5")
    } else {
      expect_match(derivatives$propagation_kernel, "rk4-method-of-steps")
      expect_false(grepl("radau", derivatives$propagation_kernel, ignore.case = TRUE))
    }
  }

  seventeen <- nm_advan_template(17)
  expect_identical(seventeen$DDE_CONFIG$lags[[1L]]$state, 2L)
  expect_identical(seventeen$DDE_CONFIG$lags[[1L]]$algebraic, "FREE")
  result17 <- nm_simulate(seventeen, data)
  early <- data$TIME <= 2
  expect_equal(
    result17$A1[early], 100 * exp(-.2 * data$TIME[early]),
    tolerance = 2e-6
  )
  derivative17 <- nm_compile(seventeen)$prediction_derivatives(data)
  expect_true(all(is.finite(derivative17$jacobian)))
  expect_match(derivative17$propagation_kernel, "ddae-advan17")
  expect_match(derivative17$propagation_kernel, "radau-iia5")
})

test_that("ADVAN16/17 Radau paths retain parameter and delay sensitivities", {
  data <- advan_event_data(c(0, .5, 1, 2, 3, 4))
  for (advan in c(16L, 17L)) {
    model <- nm_advan_template(advan)
    compiled <- nm_compile(model)$prediction_derivatives(data)
    row <- match(3, data$TIME)
    columns <- seq_len(nrow(model$THETAS))
    automatic <- compiled$jacobian[row, columns]
    finite_difference <- vapply(columns, function(column) {
      delta <- 1e-5 * max(1, abs(model$THETAS$Value[[column]]))
      upper <- lower <- model
      upper$THETAS$Value[[column]] <- upper$THETAS$Value[[column]] + delta
      lower$THETAS$Value[[column]] <- lower$THETAS$Value[[column]] - delta
      (nm_simulate(upper, data)$IPRED[[row]] -
         nm_simulate(lower, data)$IPRED[[row]]) / (2 * delta)
    }, numeric(1))
    expect_equal(
      unname(automatic), finite_difference,
      tolerance = if (advan == 17L) 3e-4 else 3e-5,
      info = paste0("ADVAN", advan)
    )
  }
})

test_that("ADVAN15-18 enforce their structural contracts and export NONMEM blocks", {
  expect_error(nm_model(
    INPUT = names(advan_event_data()), ADVAN = 15,
    PRED = "K=THETA(1);S1=1", DES = "DADT(1)=-K*A(1)",
    THETAS = data.frame(THETA = 1, Value = .4)
  ), "requires DAE_CONFIG")
  expect_error(nm_model(
    INPUT = names(advan_event_data()), ADVAN = 16,
    PRED = "K=THETA(1);S1=1", DES = "DADT(1)=-K*A(1)",
    THETAS = data.frame(THETA = 1, Value = .4)
  ), "requires DDE_CONFIG")
  expect_error(nm_model(
    INPUT = names(advan_event_data()), ADVAN = 17,
    PRED = "K=THETA(1);TAU1=THETA(2);S1=1",
    DES = "DADT(1)=-K*A(1)+LAG(A(1),TAU1)",
    THETAS = data.frame(THETA = 1:2, Value = c(.4, 2)),
    DDE_CONFIG = nm_dde_config(step = .05),
    EXPERIMENTAL = nm_experimental_config(TRUE)
  ), "both DDE_CONFIG and DAE_CONFIG")
  expect_error(nm_model(
    INPUT = names(advan_event_data()), ADVAN = 15,
    PRED = "K=THETA(1);TAU1=THETA(2);S1=1",
    DES = "DADT(1)=-K*FREE+LAG(FREE,TAU1)",
    ALG = "RES(1)=FREE-A(1)",
    THETAS = data.frame(THETA = 1:2, Value = c(.4, 2)),
    DDE_CONFIG = nm_dde_config(step = .05),
    DAE_CONFIG = nm_dae_config("FREE"),
    EXPERIMENTAL = nm_experimental_config(TRUE)
  ), "Combined DDE_CONFIG and DAE_CONFIG requires ADVAN17")

  dde <- nm_control_write(nm_advan_template(16))
  expect_match(dde, "(?m)^;DDE$", perl = TRUE)
  expect_match(dde, "\\$SUBROUTINES ADVAN16 TRANS1 TOL=9")
  expect_match(dde, "AP_1_1")
  expect_match(dde, "AD_1_1")
  expect_false(grepl("LAG\\s*\\(", dde, perl = TRUE))

  ddae <- nm_control_write(nm_advan_template(17))
  expect_match(ddae, "EQUILIBRIUM")
  expect_match(ddae, "\\$AESINIT")
  expect_match(ddae, "\\$AES")
  expect_match(ddae, "E\\(2\\)")
  expect_match(ddae, "AD_2_1")
})

test_that("ADVAN10 implements NONMEM VM/KM Michaelis-Menten elimination", {
  model <- nm_model(
    INPUT = names(advan_event_data()), ADVAN = 10, TRANS = 1,
    PRED = "VM=THETA(1);KM=THETA(2);V=THETA(3);S1=V",
    ERROR = "Y=F", THETAS = data.frame(THETA = 1:3, Value = c(20, 50, 20))
  )
  result <- nm_simulate(model, advan_event_data())
  amount <- result$A1[-1L]
  time <- result$TIME[-1L]
  integrated <- amount + 50 * log(amount)
  initial <- 100 + 50 * log(100) - 20 * time
  expect_equal(integrated, initial, tolerance = 2e-5)
  expect_match(attr(result, "solver"), "advan10")
  derivatives <- nm_compile(model)$prediction_derivatives(advan_event_data())
  expect_true(all(is.finite(derivatives$jacobian)))
  expect_match(derivatives$propagation_kernel, "advan10")
})

test_that("general ADVAN families enter population estimation and recover a rate", {
  recover <- function(model, data, truth, tolerance, label) {
    generating <- model
    generating$THETAS$Value[[1L]] <- truth
    simulated <- nm_simulate(generating, data, residual = FALSE)
    observed <- data$EVID == 0L
    data$DV[observed] <- simulated$IPRED[observed]
    fit <- nm_est(model, data, method = "FO", maxit = 80L, tolerance = 1e-7)
    expect_equal(fit$convergence, 0L, info = label)
    expect_equal(fit$theta[[1L]], truth, tolerance = tolerance, info = label)
  }

  graph <- nm_matrix_model(
    data.frame(id = 1:2, name = c("CENTRAL", "PERIPHERAL")),
    data.frame(
      from = c(1L, 1L, 2L), to = c(0L, 2L, 1L),
      type = "rate", parameter = c("K10", "K12", "K21")
    )
  )
  linear_data <- advan_event_data(c(0, 0.5, 1, 2, 4))
  for (advan in c(5L, 7L)) {
    model <- nm_model(
      INPUT = names(linear_data), ADVAN = advan, TRANS = 1,
      PRED = "K10=THETA(1);K12=THETA(2);K21=THETA(3);S1=THETA(4)",
      ERROR = "Y=F+ERR(1)",
      THETAS = data.frame(
        THETA = 1:4, Value = c(.12, .1, .05, 20),
        FIX = c(FALSE, TRUE, TRUE, TRUE)
      ),
      SIGMAS = data.frame(SIGMA = 1, Value = .01, FIX = TRUE), GRAPH = graph
    )
    recover(model, linear_data, .2, .01, paste0("ADVAN", advan))
  }

  stiff_data <- advan_event_data(c(0, .01, .05, .2, 1))
  for (advan in c(8L, 9L, 14L)) {
    model <- nm_model(
      INPUT = names(stiff_data), ADVAN = advan, TRANS = 1,
      DOSECMP = 1, OBSCMP = 2,
      PRED = "KFAST=THETA(1);KSLOW=THETA(2);S2=1",
      DES = "DADT(1)=-KFAST*A(1)\nDADT(2)=KFAST*A(1)-KSLOW*A(2)",
      ERROR = "Y=F+ERR(1)",
      THETAS = data.frame(
        THETA = 1:2, Value = c(30, 1), FIX = c(FALSE, TRUE)
      ),
      SIGMAS = data.frame(SIGMA = 1, Value = .01, FIX = TRUE),
      ODE_CONTROL = list(rtol = 1e-8, atol = 1e-11)
    )
    recover(model, stiff_data, 50, .5, paste0("ADVAN", advan))
  }

  nonlinear_data <- advan_event_data(c(0, .5, 1, 2, 4))
  nonlinear <- nm_model(
    INPUT = names(nonlinear_data), ADVAN = 10, TRANS = 1,
    PRED = "VM=THETA(1);KM=THETA(2);V=THETA(3);S1=V",
    ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(
      THETA = 1:3, Value = c(12, 50, 20), FIX = c(FALSE, TRUE, TRUE)
    ),
    SIGMAS = data.frame(SIGMA = 1, Value = .01, FIX = TRUE),
    ODE_CONTROL = list(rtol = 1e-9, atol = 1e-11)
  )
  recover(nonlinear, nonlinear_data, 20, .2, "ADVAN10")
})

test_that("general linear NONMEM streams retain their model graph", {
  control <- c(
    "$PROBLEM ADVAN5 graph",
    "$INPUT ID TIME EVID AMT CMT DV MDV",
    "$DATA data.csv",
    "$SUBROUTINES ADVAN5 TRANS1",
    "$MODEL",
    "  COMP=(CENTRAL,DEFDOSE,DEFOBSERVATION)",
    "  COMP=(PERIPHERAL)",
    "$PK",
    "  K10=THETA(1)",
    "  K12=THETA(2)",
    "  K21=THETA(3)",
    "  S1=THETA(4)",
    "$ERROR Y=F",
    "$THETA 0.2 0.1 0.05 20"
  )
  imported <- nm_control_read(control)
  expect_s3_class(imported$model$GRAPH, "nm_matrix_model")
  expect_equal(imported$model$DOSECMP, 1L)
  expect_equal(imported$model$OBSCMP, 1L)
  written <- nm_control_write(imported$model)
  expect_match(written, "\\$MODEL")
  expect_match(written, "COMP=\\(CENTRAL,DEFDOSE,DEFOBSERVATION\\)")

  programmatic <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "CMT", "DV", "MDV"),
    ADVAN = 5, TRANS = 1, PRED = "K10=THETA(1);S1=1", ERROR = "Y=F",
    THETAS = data.frame(THETA = 1, Value = .2),
    GRAPH = nm_matrix_model(
      data.frame(id = 1:2, name = c("PERIPHERAL", "PERIPHERALLY")),
      data.frame(from = 1L, to = 0L, type = "rate", parameter = "K10")
    )
  )
  generated <- nm_control_write(programmatic)
  names <- regmatches(
    generated,
    gregexpr("(?<=COMP=\\()[A-Z0-9_]+", generated, perl = TRUE)
  )[[1L]]
  expect_true(all(nchar(names) <= 8L))
  expect_equal(anyDuplicated(names), 0L)
})

test_that("ODE imports honour $MODEL dose and observation flags", {
  control <- c(
    "$PROBLEM nonlinear oral model",
    "$SUBROUTINES ADVAN13",
    "$INPUT ID TIME EVID AMT CMT DV MDV WT",
    "$MODEL",
    "  COMP=(DEPOT,DEFDOSE)",
    "  COMP=(CENTRAL,DEFOBSERVATION)",
    "$THETA 1",
    "$OMEGA 0.1",
    "$SIGMA 0.1",
    "$PK KA=THETA(1); V=10",
    "$DES",
    "  DADT(1)=-KA*A(1)",
    "  DADT(2)=KA*A(1)-A(2)",
    "$ERROR Y=F*(1+ERR(1))"
  )
  imported <- nm_control_read(control, strict = FALSE)
  expect_equal(imported$model$DOSECMP, 1L)
  expect_equal(imported$model$OBSCMP, 2L)
})

test_that("ADVAN9 AES imports require an explicit safe DAE translation", {
  control <- c(
    "$PROBLEM ADVAN9 equilibrium",
    "$INPUT ID TIME EVID AMT CMT DV MDV",
    "$DATA data.csv",
    "$SUBROUTINES ADVAN9 TRANS1 TOL=9",
    "$MODEL COMP=(CENTRAL DEFDOSE DEFOBS) COMP=(BOUND EQUILIBRIUM)",
    "$PK K=THETA(1); BIND=THETA(2); S1=1",
    "$DES DADT(1)=-K*A(2)",
    "$AESINITIAL A(2)=A(1)/(1+BIND)",
    "$AES E(2)=A(2)-A(1)/(1+BIND)",
    "$ERROR Y=F",
    "$THETA 0.4 1"
  )
  expect_error(
    nm_control_read(control),
    "translation to LibeRation ALG and DAE_CONFIG"
  )
  preserved <- nm_control_read(control, strict = FALSE)
  expect_null(preserved$model)
  expect_true(preserved$compatibility$requires_manual_translation)
  expect_true(all(c("AES", "AESINITIAL") %in%
                    preserved$compatibility$preserved_records))
})
