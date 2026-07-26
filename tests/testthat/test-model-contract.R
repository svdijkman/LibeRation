test_that("versioned model contracts rebuild semantic models", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT"), ADVAN = 1,
    PRED = "CL=THETA(1);V=THETA(2);S1=V", ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(THETA = 1:2, Value = c(2, 20)),
    SIGMAS = data.frame(SIGMA = 1, Value = 0.1)
  )
  contract <- nm_model_to_contract(model)
  rebuilt <- nm_model_from_contract(contract)
  expect_identical(contract$schema, "liberation.model")
  expect_identical(contract$version, 3L)
  expect_s3_class(rebuilt, "nm_model")
  expect_equal(rebuilt$PRED, model$PRED)
  expect_equal(rebuilt$THETAS, model$THETAS)
  expect_false(any(c("pred_ir", "des_ir", "error_ir") %in% names(contract$fields)))
})

test_that("contract v3 preserves combined PK and post-ADVAN PRED sources", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT"), ADVAN = 1,
    PRED_MODE = "pk_pred",
    PK_SOURCE = "CL=THETA(1);V=THETA(2);S1=V",
    PRED_SOURCE = "F=F_ADVAN+A(1)/100",
    THETAS = data.frame(THETA = 1:2, Value = c(2, 20))
  )
  contract <- nm_model_to_contract(model)
  rebuilt <- nm_model_from_contract(contract)
  expect_identical(rebuilt$PRED_MODE, "pk_pred")
  expect_identical(rebuilt$PK_SOURCE, model$PK_SOURCE)
  expect_identical(rebuilt$PRED_SOURCE, model$PRED_SOURCE)
  expect_null(contract$fields$post_pred_ir)
})

test_that("contracts retain every advanced semantic configuration", {
  acknowledgement <- nm_experimental_config(enabled = TRUE)
  component <- nm_component(
    "turnover", "linear_spline", inputs = "A_1", outputs = "TURN",
    knots = c(0, 100), values = c(0, 1), scope = "des"
  )
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "DV", "MDV"),
    ADVAN = 6, USE_ODE = TRUE,
    PRED = "K=THETA(1)",
    DES = "DADT(1)=-K*A(1)",
    THETAS = data.frame(THETA = 1:2, Value = c(0.1, 10)),
    COMPONENTS = component, EXPERIMENTAL = acknowledgement
  )
  rebuilt <- nm_model_from_contract(nm_model_to_contract(model))
  expect_s3_class(rebuilt$COMPONENTS[[1]], "nm_component")
  expect_s3_class(rebuilt$EXPERIMENTAL, "nm_experimental_config")
  expect_equal(rebuilt$DES, model$DES)
})

test_that("contracts reject executable and unknown class content", {
  model <- nm_model(
    INPUT = c("ID", "TIME"), ADVAN = 1, PRED = "CL=THETA(1);V=10;S1=V",
    THETAS = data.frame(THETA = 1, Value = 2)
  )
  contract <- nm_model_to_contract(model)
  contract$fields$PRED <- identity
  expect_error(nm_model_from_contract(contract), "executable")
  contract <- nm_model_to_contract(model)
  contract$fields$INPUT <- structure(list("ID"), class = "unknown_contract_class")
  expect_error(nm_model_from_contract(contract), "unsupported class")
})
