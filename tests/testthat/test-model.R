test_that("legacy model syntax compiles to serializable IR", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "CMT", "DV"),
    ADVAN = 1, TRANS = 2, DOSECMP = 1, OBSCMP = 1,
    PRED = "CL = THETA(1) * exp(ETA(1))\nV = THETA(2)\nS1 = V",
    ERROR = "Y = F * (1 + ERR(1)) + ERR(2)",
    THETAS = data.frame(THETA = 1:2, Value = c(2, 20)),
    OMEGAS = data.frame(OMEGA = 1, Value = 0.1),
    SIGMAS = data.frame(SIGMA = 1:2, Value = c(0.04, 0.01))
  )
  expect_s3_class(model, "nm_model")
  expect_s3_class(model$pred_ir, "libertad_ir")
  expect_equal(model$ERROR_TYPE, "combined")
  roundtrip <- unserialize(serialize(model, NULL))
  expect_equal(roundtrip$PRED, model$PRED)
})

test_that("layout is separate from semantic graph", {
  model <- nm_model(
    INPUT = c("ID", "TIME"), ADVAN = 1,
    PRED = "CL = THETA(1)\nV = THETA(2)",
    THETAS = data.frame(THETA = 1:2, Value = c(2, 20)),
    ERROR = "Y = F", LAYOUT = list(x = 10, y = 20)
  )
  expect_equal(model$GRAPH$compartments$name, "CENTRAL")
  expect_equal(model$LAYOUT$x, 10)
})

test_that("matrix graphs validate stable ids and flow semantics", {
  graph <- nm_matrix_model(
    compartments = data.frame(
      id = c(10, 20), name = c("CENTRAL", "PERIPHERAL"),
      volume_parameter = c("VC", "VP"),
      scale_parameter = c("VC", "VP")
    ),
    flows = data.frame(
      from = c(10, 10, 20), to = c(0, 20, 10),
      type = "clearance", parameter = c("CL", "Q", "Q")
    )
  )
  expect_s3_class(graph, "nm_matrix_model")
  expect_error(
    nm_matrix_model(
      data.frame(id = 1, name = "A"),
      data.frame(from = 1, to = 2, type = "rate", parameter = "K")
    ),
    "to.*compartment"
  )
})

test_that("correlated OMEGA requires a complete positive-definite triangle", {
  model <- nm_model(
    INPUT = c("ID", "TIME"), ADVAN = 1,
    PRED = "K=THETA(1)*exp(ETA(1)); V=THETA(2)*exp(ETA(2))",
    ERROR = "Y=F+ERR(1)", THETAS = data.frame(THETA = 1:2, Value = c(1, 10)),
    OMEGAS = data.frame(
      OMEGA = 1:3, ROW = c(1, 2, 2), COL = c(1, 1, 2),
      Value = c(0.1, -0.02, 0.2)
    ),
    SIGMAS = data.frame(SIGMA = 1, Value = 0.1)
  )
  expect_equal(model$n_eta, 2)
  expect_error(
    nm_model(
      INPUT = c("ID", "TIME"), ADVAN = 1,
      PRED = "K=THETA(1)*exp(ETA(1)); V=THETA(2)*exp(ETA(2))",
      ERROR = "Y=F+ERR(1)", THETAS = data.frame(THETA = 1:2, Value = c(1, 10)),
      OMEGAS = data.frame(
        OMEGA = 1:2, ROW = c(1, 2), COL = c(1, 2), Value = c(0.1, 0.2)
      ), SIGMAS = data.frame(SIGMA = 1, Value = 0.1),
      LIK_CONFIG = nm_lik_config(omega = "full")
    ),
    "complete lower triangle"
  )
})

test_that("fixed zero OMEGA is available for deterministic validation only", {
  model <- nm_model(
    INPUT = c("ID", "TIME"), ADVAN = 1,
    PRED = "K=THETA(1)*exp(ETA(1)); V=THETA(2)", ERROR = "Y=F",
    THETAS = data.frame(THETA = 1:2, Value = c(1, 10)),
    OMEGAS = data.frame(OMEGA = 1, Value = 0, FIX = TRUE)
  )
  expect_equal(model$OMEGAS$Value, 0)
  expect_error(
    nm_model(
      INPUT = c("ID", "TIME"), ADVAN = 1,
      PRED = "K=THETA(1)*exp(ETA(1)); V=THETA(2)", ERROR = "Y=F",
      THETAS = data.frame(THETA = 1:2, Value = c(1, 10)),
      OMEGAS = data.frame(OMEGA = 1, Value = 0, FIX = FALSE)
    ),
    "zero variances must be FIXed"
  )
})

test_that("mixture definitions normalize probabilities and expose MIXNUM", {
  mixture <- nm_mixture(c(7, 3), c("A", "B"))
  expect_equal(mixture$probability, c(0.7, 0.3))
  model <- nm_model(
    INPUT = c("ID", "TIME"), ADVAN = 1,
    PRED = "K=ifelse(MIXNUM==1,THETA(1),THETA(2)); V=THETA(3)",
    ERROR = "Y=F+ERR(1)", THETAS = data.frame(THETA = 1:3, Value = c(1, 2, 10)),
    SIGMAS = data.frame(SIGMA = 1, Value = 0.1),
    LIK_CONFIG = nm_lik_config(mixtures = mixture)
  )
  expect_true("MIXNUM" %in% model$pred_ir$input_names)
})

test_that("generated outputs are discovered statically and validated", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT"), OUTPUT = c("PRED", "CL", "V", "A1"),
    ADVAN = 1, PRED = "CL=THETA(1)*exp(ETA(1)); V=THETA(2); K=CL/V; S1=V",
    ERROR = "Y=F", THETAS = data.frame(THETA = 1:2, Value = c(2, 20)),
    OMEGAS = data.frame(OMEGA = 1, Value = 0.1)
  )
  catalog <- nm_model_outputs(model)
  expect_true(all(c("PRED", "IPRED", "CWRES", "ETA1", "A1", "CL", "V", "K", "S1") %in%
                    catalog$name))
  expect_equal(catalog$name[catalog$selected], c("PRED", "A1", "CL", "V"))
  expect_error(
    nm_model(
      INPUT = c("ID", "TIME"), OUTPUT = "NOT_CREATED", ADVAN = 1,
      PRED = "CL=THETA(1); V=THETA(2)", ERROR = "Y=F",
      THETAS = data.frame(THETA = 1:2, Value = c(2, 20))
    ),
    "Unknown OUTPUT"
  )
})
