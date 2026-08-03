test_that("covariance repairs are explicit and auditable", {
  roundoff <- matrix(c(1, 1 + 1e-12, 1 + 1e-12, 1), 2)
  none <- nm_covariance_repair(roundoff, "none", tolerance = 1e-10)
  expect_gte(min(eigen(none$matrix, symmetric = TRUE)$values), -2e-12)
  expect_equal(none$diagnostics$method, "none")

  indefinite <- matrix(c(1, 2, 2, 1), 2)
  expect_error(nm_covariance_repair(indefinite, "none"), "materially indefinite")
  expect_warning(clipped <- nm_covariance_repair(indefinite, "clip"), "material")
  expect_gte(min(clipped$diagnostics$repaired_eigenvalues), -1e-12)
  expect_true(clipped$diagnostics$material_indefiniteness)

  jittered <- nm_covariance_repair(indefinite, "jitter")
  expect_gt(jittered$diagnostics$diagonal_shift, 0)
  expect_gte(min(jittered$diagnostics$repaired_eigenvalues), 0)

  higham <- nm_covariance_repair(
    indefinite, "higham", preserve_diagonal = TRUE,
    tolerance = 1e-9, max_iterations = 200L
  )
  expect_equal(diag(higham$matrix), diag(indefinite), tolerance = 1e-10)
  expect_gte(min(higham$diagnostics$repaired_eigenvalues), -1e-10)
  expect_true(higham$diagnostics$converged)
  expect_gt(higham$diagnostics$adjustment_norm, 0)

  modified <- nm_covariance_repair(indefinite, "modified_cholesky")
  expect_gt(min(modified$diagnostics$repaired_eigenvalues), 0)
  expect_true(is.finite(modified$diagnostics$condition_number))
  expect_equal(sort(modified$diagnostics$permutation), 1:2)
  expect_error(
    nm_covariance_repair(indefinite, "modified_cholesky", preserve_diagonal = TRUE),
    "cannot preserve"
  )
})

test_that("repair sensitivity reports downstream changes and failures", {
  indefinite <- matrix(c(1, 2, 2, 1), 2)
  sensitivity <- nm_covariance_repair_sensitivity(
    indefinite,
    statistic = function(covariance) c(trace = sum(diag(covariance))),
    methods = c("none", "jitter", "higham", "modified_cholesky"),
    reference_method = "higham",
    preserve_diagonal = FALSE
  )
  expect_s3_class(sensitivity, "nm_covariance_sensitivity")
  expect_equal(sensitivity$comparison$status[[1L]], "failed")
  expect_true(all(c("estimate", "difference", "relative_matrix_adjustment") %in%
                    names(sensitivity$comparison)))
  expect_equal(
    sensitivity$comparison$difference[
      sensitivity$comparison$method == "higham" &
        sensitivity$comparison$statistic == "trace"
    ],
    0
  )
})
