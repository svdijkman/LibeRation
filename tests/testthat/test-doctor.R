test_that("ecosystem diagnostics report versions without mutating state", {
  report <- liber_doctor(verbose = FALSE)
  expect_s3_class(report, "liber_diagnostics")
  expect_true(all(c("package", "installed", "version", "expected", "compatible") %in%
                  names(report$packages)))
  expect_identical(as.integer(report$contracts$model), 2L)
  expect_true(is.list(report$engine))
})

test_that("support declarations include an evidence tier", {
  support <- nm_support_matrix()
  expect_true(all(c("status", "validation", "recommended_use") %in% names(support)))
  expect_identical(support$validation[support$feature == "FO"], "reference-validated")
  expect_true(all(support$recommended_use[support$status == "experimental"] ==
                  "experimental research only"))
})
