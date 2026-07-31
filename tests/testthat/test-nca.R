test_that("native NCA recovers a mono-exponential profile", {
  time <- seq(0, 24, by = 2)
  data <- data.frame(ID = 1L, TIME = time, DV = 10 * exp(-0.1 * time))
  result <- nm_nca(data, dose = 100, route = "bolus")

  expect_s3_class(result, "nm_nca")
  expect_equal(result$backend, "native")
  expect_equal(result$results$CMAX, 10, tolerance = 1e-10)
  expect_equal(result$results$LAMBDA_Z, 0.1, tolerance = 1e-10)
  expect_equal(result$results$HALF_LIFE, log(2) / 0.1, tolerance = 1e-9)
  expect_equal(result$results$AUCINF_OBS, 100, tolerance = 1e-8)
  expect_equal(result$results$CL, 1, tolerance = 1e-8)
})

test_that("native NCA groups profiles and calculates partial exposure", {
  data <- rbind(
    data.frame(ID = "A", TIME = 0:4, DV = c(0, 5, 4, 3, 2)),
    data.frame(ID = "B", TIME = 0:4, DV = c(0, 4, 3, 2, 1))
  )
  result <- nm_nca(
    data, partial_auc = data.frame(name = "zero_to_two", start = 0, end = 2),
    method = "linear"
  )
  expect_equal(nrow(result$results), 2L)
  expect_equal(result$results$PAUC_zero_to_two, c(7, 5.5))
})

test_that("NCA preprocessing makes duplicate and BLQ policies explicit", {
  data <- data.frame(ID = 1, TIME = c(0, 1, 1, 2), DV = c(-1, 4, 6, 3))
  result <- nm_nca(data, duplicate = "mean", blq = "zero", method = "linear")
  expect_equal(result$results$N, 3L)
  expect_equal(result$results$CMAX, 5)
  expect_error(nm_nca(data, duplicate = "error"), "Duplicate")
})

test_that("ncar validation is auditable when the reference packages are available", {
  skip_if_not_installed("ncar")
  skip_if_not_installed("NonCompart")
  data <- data.frame(ID = 1, TIME = seq(0, 12, by = 1), DV = 8 * exp(-0.2 * seq(0, 12, by = 1)))
  result <- nm_nca(data, dose = 100, route = "extravascular", validate = TRUE)
  expect_true(is.data.frame(result$validation))
  expect_true(all(c("metric", "native", "ncar", "relative_difference") %in% names(result$validation)))
  compared <- result$validation[
    result$validation$metric %in% c("CMAX", "TMAX", "AUCLAST", "AUCINF_OBS", "LAMBDA_Z", "HALF_LIFE"),
  ]
  expect_true(all(is.finite(compared$relative_difference)))
  expect_lt(max(abs(compared$relative_difference)), 1e-8)

  reference <- nm_nca(data, dose = 100, route = "extravascular", engine = "ncar")
  expect_match(reference$backend, "ncar")
  expect_equal(reference$results$AUCLAST, result$results$AUCLAST, tolerance = 1e-8)
})

test_that("native NCA agrees with ncar on the multi-subject Theoph fixture", {
  skip_if_not_installed("ncar")
  skip_if_not_installed("NonCompart")
  data("Theoph", package = "datasets")
  result <- nm_nca(
    Theoph, time = "Time", concentration = "conc", id = "Subject",
    dose = 320, route = "oral", validate = TRUE,
    dose_unit = "mg", concentration_unit = "mg/L", time_unit = "h"
  )
  comparison <- result$validation
  invariant <- comparison[comparison$metric %in% c(
    "CMAX", "TMAX", "AUCLAST", "LAMBDA_Z", "AUCINF_OBS", "HALF_LIFE"
  ), , drop = FALSE]
  expect_true(all(is.finite(invariant$relative_difference)))
  expect_lt(max(abs(invariant$relative_difference)), 1e-8)
  moment <- comparison[comparison$metric %in% c("AUMCLAST", "AUMCINF_OBS", "MRT"), , drop = FALSE]
  expect_lt(max(abs(moment$relative_difference), na.rm = TRUE), 0.002)
})
