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
    expect_false(fit$objective_comparable)
    expect_match(fit$objective_type, "log_posterior")
    expect_true(is.numeric(fit$diagnostics$divergences))
  }
})

test_that("optimized NUTS resolves to multinomial trajectory sampling", {
  fixture <- estimation_fixture(FALSE)
  fit <- nm_est(
    fixture$model, fixture$data, method = "NUTS",
    n_warmup = 6, n_sample = 8, n_chains = 1, n_thin = 1,
    max_depth = 3, seed = 193, numerical_mode = "liber_optimized"
  )
  expect_identical(fit$diagnostics$nuts$variant_resolved, "multinomial")
  expect_identical(
    fit$diagnostics$nuts$trajectory_sampling,
    "multinomial progressive sampling"
  )
  expect_match(fit$diagnostics$nuts$mass_adaptation, "windowed")
  expect_true(all(is.finite(fit$chain[, "LOG_POSTERIOR"])))
})

test_that("generalized NUTS termination uses endpoint velocities and momentum sum", {
  stop_rule <- LibeRation:::.nm_nuts_generalized_stop
  expect_true(stop_rule(c(1, 0), c(0.5, 0.5), c(2, 1), c(1, 2)))
  expect_false(stop_rule(c(-1, 0), c(0.5, 0.5), c(2, 1), c(1, 2)))
  windows <- LibeRation:::.nm_hmc_mass_windows(100L)
  expect_equal(attr(windows, "initial_buffer"), 15L)
  expect_true(all(diff(windows) > 0))
  expect_lte(tail(windows, 1L), 90L)
})

test_that("Stan-aligned warmup uses fast, expanding slow, and terminal phases", {
  schedule <- LibeRation:::.nm_hmc_warmup_schedule(500L)
  expect_equal(schedule$initial_buffer, 75L)
  expect_equal(schedule$terminal_buffer, 50L)
  expect_equal(schedule$windows, c(100L, 150L, 250L, 450L))

  short <- LibeRation:::.nm_hmc_warmup_schedule(100L)
  expect_equal(short$initial_buffer, 15L)
  expect_equal(short$terminal_buffer, 10L)
  expect_equal(tail(short$windows, 1L), 90L)
})

test_that("dense and block HMC metrics retain the requested covariance structure", {
  set.seed(120)
  first <- stats::rnorm(100)
  samples <- cbind(first, first + stats::rnorm(100, sd = 0.1), stats::rnorm(100))
  dense <- LibeRation:::.nm_hmc_adapted_metric(samples, "dense")
  blocked <- LibeRation:::.nm_hmc_adapted_metric(
    samples, "block", list(1:2, 3L)
  )
  expect_s3_class(dense, "nm_hmc_metric")
  expect_gt(abs(dense$mass[1, 2]), 0.01)
  expect_gt(abs(blocked$mass[1, 2]), 0.01)
  expect_equal(blocked$mass[1:2, 3], c(0, 0), tolerance = 1e-12)
  expect_true(all(eigen(dense$mass, symmetric = TRUE)$values > 0))
})

test_that("native optimized NUTS reports dense metric, energy, and E-BFMI", {
  fixture <- estimation_fixture(FALSE)
  fit <- nm_est(
    fixture$model, fixture$data, method = "NUTS",
    numerical_mode = "liber_optimized", n_warmup = 24, n_sample = 6,
    n_chains = 1, max_depth = 3, seed = 277, hmc_metric = "dense"
  )
  expect_identical(fit$diagnostics$nuts$metric_resolved, "dense")
  expect_true("energy" %in% names(fit$diagnostics$chain[[1L]]$trace))
  expect_true(is.finite(fit$diagnostics$ebfmi[[1L]]))
  expect_identical(
    fit$diagnostics$chain[[1L]]$initialization$source,
    "model initial"
  )
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

test_that("optimized HMC and NUTS use guarded MU non-centred geometry", {
  fixture <- estimation_fixture(FALSE)
  model <- nm_model(
    INPUT = fixture$model$INPUT, ADVAN = 1,
    PRED = "V=THETA(2); S1=V", ERROR = "Y=F+ERR(1)",
    THETAS = fixture$model$THETAS, OMEGAS = fixture$model$OMEGAS,
    SIGMAS = fixture$model$SIGMAS,
    MU = nm_mu(1, "log(THETA(1))", parameter = "CL"),
    NUMERICAL_MODE = "liber_optimized"
  )
  context <- LibeRation:::.nm_estimation_context(
    model, fixture$data, n_cores = 1, method = "HMC"
  )
  map <- LibeRation:::.nm_outer_map(context$model)
  plan <- LibeRation:::.nm_hmc_geometry_plan(context, map, "auto")
  expect_true(plan$enabled)
  expect_identical(plan$resolved, "mu_noncentered")

  target <- LibeRation:::.nm_hmc_target(context, map, plan)
  q <- target$initial + seq_along(target$initial) * 0.0007
  evaluated <- target$evaluate(q)
  step <- 1e-6
  numerical <- vapply(seq_along(q), function(index) {
    plus <- minus <- q
    plus[[index]] <- plus[[index]] + step
    minus[[index]] <- minus[[index]] - step
    (target$evaluate(plus)$logp - target$evaluate(minus)$logp) / (2 * step)
  }, numeric(1))
  expect_equal(evaluated$gradient, numerical, tolerance = 3e-4)
  native <- LibeRation:::.nm_hmc_native_target_eval(target, context, map, q)
  expect_equal(native$logp, evaluated$logp, tolerance = 2e-11)
  expect_equal(native$gradient, evaluated$gradient, tolerance = 2e-10)
  expect_equal(as.numeric(native$eta), as.numeric(evaluated$eta), tolerance = 2e-12)

  fit <- nm_est(
    model, fixture$data, method = "HMC", n_warmup = 2, n_sample = 3,
    n_chains = 1, n_leapfrog = 2, seed = 19
  )
  expect_identical(fit$diagnostics$geometry$resolved, "mu_noncentered")

  centred_model <- model
  centred_model$NUMERICAL_MODE <- "nonmem_compatibility"
  centred_context <- LibeRation:::.nm_estimation_context(
    centred_model, fixture$data, n_cores = 1, method = "HMC"
  )
  centred_plan <- LibeRation:::.nm_hmc_geometry_plan(
    centred_context, LibeRation:::.nm_outer_map(centred_model), "auto"
  )
  expect_false(centred_plan$enabled)
  expect_match(centred_plan$reason, "liber_optimized")
  expect_error(
    LibeRation:::.nm_hmc_geometry_plan(
      centred_context, LibeRation:::.nm_outer_map(centred_model),
      "mu_noncentered"
    ),
    "unavailable"
  )
})

test_that("MU non-centred geometry differentiates correlated OMEGA", {
  fixture <- estimation_fixture(TRUE)
  model <- nm_model(
    INPUT = fixture$model$INPUT, ADVAN = 1,
    PRED = "S1=V", ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(
      THETA = 1:2, Value = c(2, 20), FIX = FALSE,
      LOWER = c(0.01, 0.1), UPPER = c(20, 200)
    ),
    OMEGAS = data.frame(
      OMEGA = 1:3, ROW = c(1, 2, 2), COL = c(1, 1, 2),
      Value = c(0.09, 0.02, 0.16), FIX = FALSE
    ),
    SIGMAS = data.frame(SIGMA = 1, Value = 0.2, FIX = TRUE),
    MU = rbind(
      nm_mu(1, "log(THETA(1))", parameter = "CL"),
      nm_mu(2, "log(THETA(2))", parameter = "V")
    ),
    NUMERICAL_MODE = "liber_optimized"
  )
  context <- LibeRation:::.nm_estimation_context(
    model, fixture$data, n_cores = 1, method = "NUTS"
  )
  map <- LibeRation:::.nm_outer_map(context$model)
  plan <- LibeRation:::.nm_hmc_geometry_plan(context, map, "auto")
  target <- LibeRation:::.nm_hmc_target(context, map, plan)
  q <- target$initial + seq_along(target$initial) * 0.0003
  evaluated <- target$evaluate(q)
  step <- 1e-6
  numerical <- vapply(seq_along(q), function(index) {
    plus <- minus <- q
    plus[[index]] <- plus[[index]] + step
    minus[[index]] <- minus[[index]] - step
    (target$evaluate(plus)$logp - target$evaluate(minus)$logp) / (2 * step)
  }, numeric(1))
  native <- LibeRation:::.nm_hmc_native_target_eval(target, context, map, q)
  expect_equal(evaluated$gradient, numerical, tolerance = 4e-4)
  expect_equal(native$logp, evaluated$logp, tolerance = 3e-11)
  expect_equal(native$gradient, evaluated$gradient, tolerance = 3e-10)
})

test_that("optimized nonparametric weights use the native interior solver", {
  loglik <- matrix(c(
    -0.1, -1.2, -2.0,
    -1.1, -0.2, -1.5,
    -1.7, -0.8, -0.3
  ), nrow = 3L, byrow = TRUE)
  fit <- .nm_np_weights_interior(
    loglik, initial = c(0.3, 0.4, 0.3), maxit = 100L,
    tolerance = 1e-9
  )
  expect_identical(fit$solver, "primal-dual-barrier-newton")
  expect_identical(fit$backend, "native-cpp")
  expect_equal(sum(fit$weights), 1, tolerance = 2e-12)
  expect_true(all(fit$weights > 0))
  expect_equal(rowSums(fit$responsibilities), rep(1, 3), tolerance = 2e-12)
})

test_that("compatibility NP weight EM is coordinated in one native call", {
  loglik <- matrix(c(
    -1.2, -0.5, -1.8,
    -0.4, -1.1, -1.5,
    -1.0, -0.8, -0.6
  ), nrow = 3L, byrow = TRUE)
  native <- .nm_np_weights(
    loglik, c(0.2, 0.5, 0.3), maxit = 200L,
    tolerance = 1e-12, optimized = FALSE
  )
  expect_identical(native$solver, "expectation-maximization")
  expect_equal(sum(native$weights), 1, tolerance = 2e-15)
  expect_equal(
    rowSums(native$responsibilities), rep(1, 3), tolerance = 2e-15
  )
  expect_true(is.finite(native$log_likelihood))
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
  expect_true(
    npml$diagnostics$nonparametric$pruning$fixed_support_preserved
  )
  expect_gt(
    npml$diagnostics$nonparametric$likelihood_grid$exact_state_cache_hits,
    0L
  )
  expect_gt(
    npml$diagnostics$nonparametric$likelihood_grid$evaluations,
    0L
  )
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
  expect_match(
    npag$diagnostics$nonparametric$estimator_identity,
    "adaptive-grid maximum likelihood"
  )
  expect_true(npag$diagnostics$nonparametric$pruning$likelihood_guard)
  expect_gt(
    npag$diagnostics$nonparametric$likelihood_grid$exact_state_cache_hits,
    0L
  )
  expect_true(all(diff(npag$nonparametric$history$log_likelihood) >= -1e-6))

  free_fixture <- estimation_fixture(FALSE)
  population_fit <- nm_est(
    free_fixture$model, free_fixture$data, method = "NPML", maxit = 3,
    np_supports = matrix(c(-0.3, 0.3), ncol = 1), np_cycles = 1,
    np_weight_maxit = 20
  )
  expect_s3_class(population_fit, "nm_fit")
  expect_true(is.finite(population_fit$objective))
})

test_that("compatibility NPML never applies the NPAG support ceiling by default", {
  fixture <- estimation_fixture(TRUE)
  supports <- matrix(c(-1.5, -0.4, 0, 0.4, 1.5), ncol = 1L)
  fit <- nm_est(
    fixture$model, fixture$data, method = "NPML",
    np_supports = supports, np_cycles = 1L,
    np_max_support = 2L, np_estimate_population = FALSE,
    np_weight_maxit = 100L, numerical_mode = "nonmem_compatibility"
  )
  expect_equal(nrow(fit$nonparametric$supports), nrow(supports))
  expect_equal(
    fit$diagnostics$nonparametric$pruning$resolved_minimum_weight, 0
  )
  expect_true(fit$diagnostics$nonparametric$pruning$fixed_support_preserved)

  expect_warning(
    threshold_fit <- nm_est(
      fixture$model, fixture$data, method = "NPML",
      np_supports = supports, np_cycles = 1L,
      np_min_weight = 0.25, np_estimate_population = FALSE,
      np_weight_maxit = 100L
    ),
    "ignored for NPML"
  )
  expect_equal(nrow(threshold_fit$nonparametric$supports), nrow(supports))
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
