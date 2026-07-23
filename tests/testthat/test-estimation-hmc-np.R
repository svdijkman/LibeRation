test_that("HMC and NUTS use the exact joint target and retain posterior diagnostics", {
  fixture <- estimation_fixture(TRUE)
  for (method in c("HMC", "NUTS")) {
    fit <- nm_est(
      fixture$model, fixture$data, method = method,
      n_warmup = 5, n_sample = 8, n_chains = 1, n_thin = 1,
      n_leapfrog = 3, max_depth = 3, seed = 11
    )
    expect_s3_class(fit, "nm_fit")
    expect_identical(fit$method, method)
    expect_equal(nrow(fit$chain), 8)
    expect_true(is.list(fit$posterior$population))
    expect_true(is.finite(fit$objective))
    expect_true(is.numeric(fit$diagnostics$divergences))
  }
})

test_that("HMC bounds mapping is invertible and includes finite Jacobians", {
  map <- list(lower = c(-Inf, 0, -2), upper = c(Inf, Inf, 3))
  bounds <- LibeRation:::.nm_hmc_bound_map(map)
  q <- bounds$inverse(c(1, 2, 0))
  transformed <- bounds$forward(q)
  expect_equal(transformed$value, c(1, 2, 0), tolerance = 1e-10)
  expect_true(is.finite(transformed$log_jacobian))
  expect_true(all(is.finite(transformed$log_jacobian_gradient)))
})

test_that("HMC joint gradient includes population transforms and ETA derivatives", {
  fixture <- estimation_fixture(FALSE)
  context <- LibeRation:::.nm_estimation_context(
    fixture$model, fixture$data, n_cores = 1, method = "HMC"
  )
  map <- LibeRation:::.nm_outer_map(context$model)
  target <- LibeRation:::.nm_hmc_target(context, map)
  q <- target$initial + seq_along(target$initial) * 0.001
  evaluated <- target$evaluate(q)
  step <- 1e-6
  numerical <- vapply(seq_along(q), function(index) {
    plus <- minus <- q
    plus[[index]] <- plus[[index]] + step
    minus[[index]] <- minus[[index]] - step
    (target$evaluate(plus)$logp - target$evaluate(minus)$logp) / (2 * step)
  }, numeric(1))
  expect_true(is.finite(evaluated$logp))
  expect_equal(evaluated$gradient, numerical, tolerance = 2e-4)
})

test_that("native HMC target exactly matches the retained R reference target", {
  fixture <- estimation_fixture(FALSE)
  fixture$model$LIK_CONFIG$priors <- do.call(rbind, list(
    nm_prior("THETA1", "normal", mean = 2, sd = 0.5),
    nm_prior("OMEGA1", "lognormal", mean = log(0.09), sd = 0.3)
  ))
  context <- LibeRation:::.nm_estimation_context(
    fixture$model, fixture$data, n_cores = 1, method = "HMC"
  )
  map <- LibeRation:::.nm_outer_map(context$model)
  reference <- LibeRation:::.nm_hmc_target(context, map)
  q <- reference$initial + seq_along(reference$initial) * 0.0007

  expected <- reference$evaluate(q)
  observed <- LibeRation:::.nm_hmc_native_target_eval(
    reference, context, map, q
  )

  expect_equal(observed$logp, expected$logp, tolerance = 2e-11)
  expect_equal(observed$gradient, expected$gradient, tolerance = 2e-10)
  expect_equal(observed$outer, expected$outer, tolerance = 2e-12)
  expect_equal(
    as.numeric(observed$eta), as.numeric(expected$eta), tolerance = 2e-12
  )
})

test_that("HMC differentiates the full-OMEGA Cholesky transform", {
  fixture <- estimation_fixture(TRUE)
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "DV", "MDV"),
    ADVAN = 1, DOSECMP = 1, OBSCMP = 1,
    PRED = "CL=THETA(1)*exp(ETA(1)); V=THETA(2)*exp(ETA(2)); S1=V",
    ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(THETA = 1:2, Value = c(2, 20), FIX = TRUE),
    OMEGAS = data.frame(
      OMEGA = 1:3, ROW = c(1, 2, 2), COL = c(1, 1, 2),
      Value = c(0.09, 0.02, 0.16), FIX = FALSE
    ),
    SIGMAS = data.frame(SIGMA = 1, Value = 0.2, FIX = TRUE)
  )
  context <- LibeRation:::.nm_estimation_context(
    model, fixture$data, n_cores = 1, method = "HMC"
  )
  map <- LibeRation:::.nm_outer_map(context$model)
  target <- LibeRation:::.nm_hmc_target(context, map)
  q <- target$initial + seq_along(target$initial) * 0.0005
  evaluated <- target$evaluate(q)
  step <- 1e-6
  numerical <- vapply(seq_along(q), function(index) {
    plus <- minus <- q
    plus[[index]] <- plus[[index]] + step
    minus[[index]] <- minus[[index]] - step
    (target$evaluate(plus)$logp - target$evaluate(minus)$logp) / (2 * step)
  }, numeric(1))
  expect_equal(evaluated$gradient, numerical, tolerance = 3e-4)
  native <- LibeRation:::.nm_hmc_native_target_eval(
    target, context, map, q
  )
  expect_equal(native$logp, evaluated$logp, tolerance = 2e-11)
  expect_equal(native$gradient, evaluated$gradient, tolerance = 2e-10)
})

test_that("NPML estimates fixed-support weights and NPAG adapts its support", {
  fixture <- estimation_fixture(TRUE)
  npml <- nm_est(
    fixture$model, fixture$data, method = "NPML",
    np_supports = matrix(c(-0.4, 0.4), ncol = 1),
    np_cycles = 1, np_estimate_population = FALSE, np_weight_maxit = 50
  )
  expect_s3_class(npml, "nm_fit")
  expect_equal(sum(npml$nonparametric$weights), 1, tolerance = 1e-10)
  expect_equal(nrow(npml$nonparametric$supports), 2)
  expect_equal(dim(npml$nonparametric$posterior_probabilities), c(3, 2))
  individual <- predict(npml, type = "individual")
  components <- lapply(seq_len(2), function(index) nm_simulate(
    npml$model, npml$data, theta = npml$theta,
    eta = matrix(npml$nonparametric$supports[index, ], 3, 1),
    sigma = npml$sigma, omega = npml$omega
  )$IPRED)
  component_matrix <- do.call(cbind, components)
  row_weights <- npml$nonparametric$posterior_probabilities[npml$data$.ID_INDEX, ]
  expect_equal(individual$IPRED, rowSums(component_matrix * row_weights), tolerance = 1e-10)
  expect_s3_class(nm_gof(npml), "data.frame")

  npag <- nm_est(
    fixture$model, fixture$data, method = "NPAG",
    np_points = 2, np_cycles = 2, np_estimate_population = FALSE,
    np_weight_maxit = 50, np_max_support = 10, np_max_candidates = 20
  )
  expect_s3_class(npag, "nm_fit")
  expect_identical(npag$method, "NPAG")
  expect_true(nrow(npag$nonparametric$supports) >= 1)
  expect_equal(sum(npag$nonparametric$weights), 1, tolerance = 1e-10)

  free_fixture <- estimation_fixture(FALSE)
  population_fit <- nm_est(
    free_fixture$model, free_fixture$data, method = "NPML", maxit = 3,
    np_supports = matrix(c(-0.3, 0.3), ncol = 1), np_cycles = 1,
    np_weight_maxit = 20
  )
  expect_s3_class(population_fit, "nm_fit")
  expect_true(is.finite(population_fit$objective))
})

test_that("nonparametric likelihood removes the Gaussian OMEGA density", {
  fixture <- estimation_fixture(TRUE)
  context <- LibeRation:::.nm_estimation_context(
    fixture$model, fixture$data, n_cores = 1, method = "NPML"
  )
  supports <- matrix(c(-0.4, 0.4), ncol = 1)
  parameters <- list(theta = c(2, 20), sigma = 0.2, omega = 0.09)
  first <- LibeRation:::.nm_np_loglik(context, parameters, supports)$loglik
  parameters$omega <- 0.25
  second <- LibeRation:::.nm_np_loglik(context, parameters, supports)$loglik
  expect_equal(first, second, tolerance = 1e-8)
})
