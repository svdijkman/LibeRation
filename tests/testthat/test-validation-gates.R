test_that("benchmark publication is correctness-gated", {
  reference <- list(
    objective = 100, theta = c(2, 20), omega = 0.1, sigma = 0.2,
    eta = c(-0.1, 0.1), predictions = c(1, 2, 3)
  )
  close <- reference
  close$theta <- close$theta * c(1.001, 0.999)
  validation <- nm_validation_gate(reference, close)
  expect_true(validation$passed)
  provenance <- nm_benchmark_provenance("unit", repetitions = 3, warmup = 1)
  record <- nm_benchmark_gate(validation, list(core_seconds = c(1, 1.1, 0.9)), provenance)
  expect_true(record$publishable)

  far <- close; far$omega <- 1
  failed <- nm_validation_gate(reference, far)
  expect_false(failed$passed)
  expect_false(nm_benchmark_gate(failed, list(core_seconds = 1), provenance)$publishable)
})

test_that("benchmark gates reject missing required evidence and single runs", {
  partial <- list(objective = 1, theta = 2)
  validation <- nm_validation_gate(partial, partial)
  expect_false(validation$passed)
  provenance <- nm_benchmark_provenance("single", repetitions = 1, warmup = 0)
  complete <- list(objective = 1, theta = 2, omega = 0.1, sigma = 0.2,
                   eta = 0, predictions = 1)
  passed <- nm_validation_gate(complete, complete)
  expect_false(nm_benchmark_gate(passed, list(total_seconds = 1), provenance)$publishable)
})
