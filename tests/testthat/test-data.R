test_that("ADDL doses are expanded with stable source identity", {
  data <- nm_dataset(data.frame(
    ID = 1, TIME = c(0, 12, 24), EVID = c(1, 0, 0),
    AMT = c(100, 0, 0), ADDL = c(2, 0, 0), II = c(12, 0, 0),
    CMT = 1, MDV = c(1, 0, 0)
  ))
  expect_equal(sum(data$EVID == 1), 3)
  expect_equal(sort(data$TIME[data$EVID == 1]), c(0, 12, 24))
  expect_equal(sum(data$.generated), 2)
})

test_that("every execution path materializes ADDL and native boundaries fail closed", {
  fixture <- estimation_fixture()
  raw <- data.frame(
    ID = 1, TIME = c(0, 12, 24), EVID = c(1, 0, 0),
    AMT = c(100, 0, 0), ADDL = c(2, 0, 0), II = c(12, 0, 0),
    CMT = 1, MDV = c(1, 1, 1)
  )
  deferred <- nm_dataset(raw, expand_addl = FALSE)
  expect_false(isTRUE(attr(deferred, "addl_materialized")))
  engine <- nm_compile(fixture$model)
  pointer <- engine$pointer
  expect_error(engine$pointer <- NULL, "read-only")
  expect_identical(engine$pointer, pointer)
  result <- engine$simulate(deferred)
  expect_true(isTRUE(attr(result, "addl_materialized")))
  expect_equal(sum(result$EVID == 1L), 3L)
  expect_true(all(result$ADDL == 0L))
  expect_error(
    LibeRation:::.liberation_engine_simulate(
      engine$pointer, deferred, fixture$model$THETAS$Value,
      matrix(0, 1, fixture$model$n_eta), fixture$model$SIGMAS$Value
    ),
    "requires ADDL/II doses to be materialized"
  )
})

test_that("NONMEM modelled infusion RATE codes are normalized safely", {
  data <- nm_dataset(data.frame(
    ID = 1, TIME = c(0, 1), EVID = c(1, 0), AMT = c(100, 0),
    RATE = c(-1, 0)
  ))
  expect_equal(data$RATE, c(-1, 0))
  expect_error(
    nm_dataset(data.frame(ID = 1, TIME = 0, EVID = 1, AMT = 100, RATE = -3)),
    "RATE must be"
  )
})

test_that("missing event cells use NONMEM zero defaults while DV remains missing", {
  data <- nm_dataset(data.frame(
    ID = 1, TIME = c(0, 1), EVID = c(1, NA), AMT = c(100, NA),
    RATE = c(NA, NA), DV = c(NA, 4.5)
  ))
  expect_equal(data$EVID, c(1L, 0L))
  expect_equal(data$AMT, c(100, 0))
  expect_equal(data$RATE, c(0, 0))
  expect_true(is.na(data$DV[[1L]]))
})
