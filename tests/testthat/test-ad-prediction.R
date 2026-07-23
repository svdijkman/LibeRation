theta_table_ad <- function(values) {
  data.frame(THETA = seq_along(values), Value = values)
}

central_difference <- function(fn, at, index, step = 1e-5) {
  plus <- minus <- at
  plus[[index]] <- plus[[index]] + step
  minus[[index]] <- minus[[index]] - step
  (fn(plus) - fn(minus)) / (2 * step)
}

test_that("the scalar-generic Pade exponential agrees with Eigen", {
  set.seed(2104)
  for (n in 1:5) {
    matrix <- matrix(rnorm(n * n, sd = 0.3), n, n)
    expect_equal(
      LibeRation:::.liberation_matrix_exp_pade(matrix),
      LibeRation:::.liberation_matrix_exp(matrix),
      tolerance = 2e-12
    )
  }
})

specialized_advan_cases <- function() {
  list(
    list(advan = 1L, dose = 1L, obs = 1L,
         pred = "CL=THETA(1);V=THETA(2);S1=V",
         theta = c(2, 20)),
    list(advan = 2L, dose = 1L, obs = 2L,
         pred = "KA=THETA(1);CL=THETA(2);V=THETA(3);S2=V",
         theta = c(1.1, 2, 20)),
    list(advan = 3L, dose = 1L, obs = 1L,
         pred = "CL=THETA(1);VC=THETA(2);Q=THETA(3);VP=THETA(4);S1=VC",
         theta = c(2, 20, 1.5, 30)),
    list(advan = 4L, dose = 1L, obs = 2L,
         pred = paste("KA=THETA(1)", "CL=THETA(2)", "VC=THETA(3)",
                      "Q=THETA(4)", "VP=THETA(5)", "S2=VC", sep = ";"),
         theta = c(1.1, 2, 20, 1.5, 30)),
    list(advan = 11L, dose = 1L, obs = 1L,
         pred = paste("CL=THETA(1)", "VC=THETA(2)", "Q2=THETA(3)",
                      "VP1=THETA(4)", "Q3=THETA(5)", "VP2=THETA(6)",
                      "S1=VC", sep = ";"),
         theta = c(2, 20, 1.5, 30, 0.8, 50)),
    list(advan = 12L, dose = 1L, obs = 2L,
         pred = paste("KA=THETA(1)", "CL=THETA(2)", "VC=THETA(3)",
                      "Q2=THETA(4)", "VP1=THETA(5)", "Q3=THETA(6)",
                      "VP2=THETA(7)", "S2=VC", sep = ";"),
         theta = c(1.1, 2, 20, 1.5, 30, 0.8, 50))
  )
}

test_that("specialized ADVAN tapes preserve values and exact derivatives", {
  previous <- getOption("LibeRation.specialized_advan")
  on.exit(options(LibeRation.specialized_advan = previous), add = TRUE)
  data <- data.frame(
    ID = 1, TIME = c(0, 0.25, 1, 4, 12, 24),
    EVID = c(1, 0, 0, 0, 0, 0), AMT = c(100, 0, 0, 0, 0, 0)
  )
  operation_reductions <- numeric()

  for (case in specialized_advan_cases()) {
    model <- nm_model(
      INPUT = c("ID", "TIME", "EVID", "AMT"), ADVAN = case$advan,
      DOSECMP = case$dose, OBSCMP = case$obs, PRED = case$pred,
      ERROR = "Y=F", THETAS = theta_table_ad(case$theta)
    )
    options(LibeRation.specialized_advan = TRUE)
    specialized <- nm_prediction_derivatives(model, data)
    options(LibeRation.specialized_advan = FALSE)
    general <- nm_prediction_derivatives(model, data)

    expect_identical(specialized$propagation_kernel,
                     paste0("specialized-advan", case$advan))
    expect_identical(general$propagation_kernel,
                     "general-matrix-exponential")
    expect_equal(specialized$value, general$value, tolerance = 2e-10,
                 info = paste("ADVAN", case$advan, "values"))
    expect_equal(specialized$jacobian, general$jacobian, tolerance = 2e-8,
                 info = paste("ADVAN", case$advan, "Jacobian"))
    expect_equal(specialized$value, nm_simulate(model, data)$IPRED,
                 tolerance = 2e-10, ignore_attr = TRUE,
                 info = paste("ADVAN", case$advan, "simulation"))
    operation_reductions <- c(
      operation_reductions,
      general$operation_count - specialized$operation_count
    )
  }
  expect_true(all(operation_reductions >= 0))
  expect_true(any(operation_reductions > 0))
})

test_that("subjects with the same topology share a dynamic-covariate tape", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "WT"), ADVAN = 1,
    DOSECMP = 1, OBSCMP = 1,
    PRED = "CL=THETA(1)*(WT/70)^0.75; V=THETA(2); S1=V",
    ERROR = "Y=F", THETAS = theta_table_ad(c(2, 20))
  )
  data <- data.frame(
    ID = rep(1:2, each = 4), TIME = rep(c(0, 1, 4, 12), 2),
    EVID = rep(c(1, 0, 0, 0), 2), AMT = rep(c(100, 0, 0, 0), 2),
    WT = rep(c(50, 100), each = 4)
  )
  normalized <- LibeRation:::.nm_engine_data(model, data)
  first <- LibeRation:::.nm_subject_data(normalized, 1L)
  second <- LibeRation:::.nm_subject_data(normalized, 2L)
  engine <- LibeRation:::NMEngine$new(model)
  pool <- new.env(parent = emptyenv())
  first_tape <- LibeRation:::.nm_prediction_pool_tape(
    pool, engine, first, model$THETAS$Value, numeric(), 0L
  )
  second_tape <- LibeRation:::.nm_prediction_pool_tape(
    pool, engine, second, model$THETAS$Value, numeric(), 0L
  )

  expect_identical(first_tape$pointer, second_tape$pointer)
  expect_identical(first_tape$dynamic_columns, "WT")
  expect_equal(first_tape$dynamic_parameters, nrow(first))

  evaluate <- function(subject) {
    LibeRation:::.liberation_prediction_tape_new_dynamic(
      first_tape$pointer, subject
    )
    LibeRation:::.liberation_prediction_tape_eval(
      first_tape$pointer, first_tape$point, FALSE
    )$value
  }
  expect_equal(evaluate(first), nm_simulate(model, first)$IPRED,
               tolerance = 2e-12, ignore_attr = TRUE)
  expect_equal(evaluate(second), nm_simulate(model, second)$IPRED,
               tolerance = 2e-12, ignore_attr = TRUE)
  expect_false(isTRUE(all.equal(evaluate(first), evaluate(second))))
})

test_that("specialized affine propagation covers steady-state infusions", {
  previous <- getOption("LibeRation.specialized_advan")
  on.exit(options(LibeRation.specialized_advan = previous), add = TRUE)
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "RATE", "SS", "II"),
    ADVAN = 3, DOSECMP = 1, OBSCMP = 1,
    PRED = "CL=THETA(1);VC=THETA(2);Q=THETA(3);VP=THETA(4);S1=VC",
    ERROR = "Y=F", THETAS = theta_table_ad(c(2, 20, 1.5, 30))
  )
  data <- data.frame(
    ID = 1, TIME = c(0, 2, 8, 24), EVID = c(1, 0, 0, 0),
    AMT = c(100, 0, 0, 0), RATE = c(10, 0, 0, 0),
    SS = c(1, 0, 0, 0), II = c(24, 0, 0, 0)
  )
  options(LibeRation.specialized_advan = TRUE)
  specialized <- nm_prediction_derivatives(model, data)
  options(LibeRation.specialized_advan = FALSE)
  general <- nm_prediction_derivatives(model, data)
  expect_equal(specialized$value, general$value, tolerance = 2e-9)
  expect_equal(specialized$jacobian, general$jacobian, tolerance = 2e-7)
})

test_that("CppAD differentiates the complete ADVAN1 event path", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT"), ADVAN = 1,
    DOSECMP = 1, OBSCMP = 1,
    PRED = "CL=THETA(1)*exp(ETA(1)); V=THETA(2); S1=V",
    ERROR = "Y=F", THETAS = theta_table_ad(c(2, 20)),
    OMEGAS = data.frame(OMEGA = 1, Value = 0.1)
  )
  data <- data.frame(
    ID = rep(1:2, each = 4), TIME = rep(c(0, 1, 5, 12), 2),
    EVID = rep(c(1, 0, 0, 0), 2), AMT = rep(c(100, 0, 0, 0), 2)
  )
  eta <- matrix(c(0.1, -0.2), ncol = 1)
  derivative <- nm_prediction_derivatives(model, data, eta = eta)
  expect_equal(derivative$value, nm_simulate(model, data, eta = eta)$IPRED,
               tolerance = 2e-12, ignore_attr = TRUE)

  point <- c(model$THETAS$Value, as.vector(t(eta)))
  value_at <- function(x) {
    nm_simulate(model, data, theta = x[1:2], eta = matrix(x[3:4], ncol = 1))$IPRED
  }
  numerical <- vapply(seq_along(point), function(index) {
    central_difference(value_at, point, index)
  }, numeric(nrow(data)))
  expect_equal(unname(derivative$jacobian), numerical, tolerance = 2e-7)
  expect_equal(colnames(derivative$jacobian),
               c("THETA_1", "THETA_2", "ETA_1_1", "ETA_2_1"))
})

test_that("AD remains exact through infusion and periodic steady state", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "RATE", "SS", "II"),
    ADVAN = 1, DOSECMP = 1, OBSCMP = 1,
    PRED = "CL=THETA(1); V=THETA(2); S1=V", ERROR = "Y=F",
    THETAS = theta_table_ad(c(2, 20))
  )
  data <- data.frame(
    ID = 1, TIME = c(0, 5, 10, 24), EVID = c(1, 0, 0, 0),
    AMT = c(100, 0, 0, 0), RATE = c(10, 0, 0, 0),
    SS = c(1, 0, 0, 0), II = c(24, 0, 0, 0)
  )
  derivative <- nm_prediction_derivatives(model, data)
  numerical <- vapply(1:2, function(index) {
    central_difference(
      function(theta) nm_simulate(model, data, theta = theta)$IPRED,
      c(2, 20), index
    )
  }, numeric(nrow(data)))
  expect_equal(unname(derivative$jacobian), numerical, tolerance = 1e-6)
})

test_that("arbitrary linear graphs compile to the same C++ and AD path", {
  graph <- nm_matrix_model(
    compartments = data.frame(
      id = c(10, 20), name = c("CENTRAL", "PERIPHERAL1"),
      volume_parameter = c("VC", "VP"),
      scale_parameter = c("VC", "VP")
    ),
    flows = data.frame(
      from = c(10, 10, 20), to = c(0, 20, 10),
      type = "clearance", parameter = c("CL", "Q", "Q")
    )
  )
  pred <- "CL=THETA(1); VC=THETA(2); Q=THETA(3); VP=THETA(4)"
  theta <- c(2, 20, 1, 10)
  matrix_model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT"), ADVAN = 1,
    SOLVER = "matrix", GRAPH = graph, DOSECMP = 1, OBSCMP = 1,
    PRED = pred, ERROR = "Y=F", THETAS = theta_table_ad(theta)
  )
  advan_model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT"), ADVAN = 3,
    DOSECMP = 1, OBSCMP = 1, PRED = paste0(pred, "; S1=VC"),
    ERROR = "Y=F", THETAS = theta_table_ad(theta)
  )
  data <- data.frame(
    ID = 1, TIME = c(0, 1, 4, 12, 36), EVID = c(1, 0, 0, 0, 0),
    AMT = c(100, 0, 0, 0, 0)
  )
  matrix_prediction <- nm_simulate(matrix_model, data)
  expect_equal(matrix_prediction$IPRED, nm_simulate(advan_model, data)$IPRED,
               tolerance = 2e-12)
  derivative <- nm_prediction_derivatives(matrix_model, data)
  expect_identical(derivative$propagation_kernel,
                   "general-matrix-exponential")
  numerical <- vapply(seq_along(theta), function(index) {
    central_difference(
      function(value) nm_simulate(matrix_model, data, theta = value)$IPRED,
      theta, index
    )
  }, numeric(nrow(data)))
  expect_equal(unname(derivative$jacobian), numerical, tolerance = 2e-7)
})

test_that("ADVAN6 records the accepted adaptive ODE trajectory with AD", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT"), ADVAN = 6,
    DOSECMP = 1, OBSCMP = 1,
    PRED = "K=THETA(1); S1=THETA(2)", DES = "DADT(1)=-K*A(1)",
    ERROR = "Y=F", THETAS = theta_table_ad(c(0.4, 20))
  )
  data <- data.frame(
    ID = 1, TIME = c(0, 0.2, 1, 4), EVID = c(1, 0, 0, 0),
    AMT = c(100, 0, 0, 0)
  )
  derivative <- nm_prediction_derivatives(model, data)
  expect_equal(derivative$value, nm_simulate(model, data)$IPRED,
               tolerance = 2e-10, ignore_attr = TRUE)
  numerical <- vapply(1:2, function(index) central_difference(
    function(theta) nm_simulate(model, data, theta = theta)$IPRED,
    c(0.4, 20), index, step = 2e-5
  ), numeric(nrow(data)))
  expect_equal(unname(derivative$jacobian), numerical, tolerance = 3e-6)
})

test_that("ADVAN13 stiff integration exposes exact parameter derivatives", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT"), ADVAN = 13,
    DOSECMP = 1, OBSCMP = 2,
    PRED = "KFAST=THETA(1); KSLOW=THETA(2); S2=1",
    DES = "DADT(1)=-KFAST*A(1); DADT(2)=KFAST*A(1)-KSLOW*A(2)",
    ERROR = "Y=F", THETAS = theta_table_ad(c(100, 1)),
    ODE_CONTROL = list(rtol = 2e-6, atol = 1e-9)
  )
  data <- data.frame(
    ID = 1, TIME = c(0, 0.01, 0.1), EVID = c(1, 0, 0), AMT = c(1, 0, 0)
  )
  derivative <- nm_prediction_derivatives(model, data)
  expect_equal(derivative$value, nm_simulate(model, data)$IPRED,
               tolerance = 2e-7, ignore_attr = TRUE)
  numerical <- vapply(1:2, function(index) central_difference(
    function(theta) nm_simulate(model, data, theta = theta)$IPRED,
    c(100, 1), index, step = if (index == 1) 1e-3 else 1e-5
  ), numeric(nrow(data)))
  expect_equal(unname(derivative$jacobian), numerical, tolerance = 2e-4)
})

test_that("periodic ODE steady state remains on the exact AD tape", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "RATE", "SS", "II"),
    ADVAN = 6, DOSECMP = 1, OBSCMP = 1,
    PRED = "K=THETA(1); V=THETA(2); S1=V",
    DES = "DADT(1)=-K*A(1)", ERROR = "Y=F",
    THETAS = theta_table_ad(c(0.2, 20))
  )
  data <- data.frame(
    ID = 1, TIME = c(0, 3, 8), EVID = c(1, 0, 0),
    AMT = c(120, 0, 0), RATE = c(10, 0, 0),
    SS = c(1, 0, 0), II = c(8, 0, 0)
  )
  derivative <- nm_prediction_derivatives(model, data)
  numerical <- vapply(1:2, function(index) central_difference(
    function(theta) nm_simulate(model, data, theta = theta)$IPRED,
    c(0.2, 20), index, step = 2e-5
  ), numeric(nrow(data)))
  expect_equal(derivative$value, nm_simulate(model, data)$IPRED,
               tolerance = 3e-8, ignore_attr = TRUE)
  expect_equal(unname(derivative$jacobian), numerical, tolerance = 2e-5)
})

test_that("AD differentiates through modelled infusion duration", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "RATE"), ADVAN = 6,
    DOSECMP = 1, OBSCMP = 1,
    PRED = "K=THETA(1);V=THETA(2);D1=THETA(3);S1=V",
    DES = "DADT(1)=-K*A(1)", ERROR = "Y=F",
    THETAS = theta_table_ad(c(0.1, 20, 5))
  )
  data <- data.frame(
    ID = 1, TIME = c(0, 2, 6, 10), EVID = c(1, 0, 0, 0),
    AMT = c(100, 0, 0, 0), RATE = c(-2, 0, 0, 0)
  )
  derivative <- nm_prediction_derivatives(model, data)
  numerical <- vapply(1:3, function(index) central_difference(
    function(theta) nm_simulate(model, data, theta = theta)$IPRED,
    c(0.1, 20, 5), index, step = 2e-5
  ), numeric(nrow(data)))
  expect_equal(derivative$value, nm_simulate(model, data)$IPRED,
               tolerance = 2e-8, ignore_attr = TRUE)
  expect_equal(unname(derivative$jacobian), numerical, tolerance = 5e-5)
})
