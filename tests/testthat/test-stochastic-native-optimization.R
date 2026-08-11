test_that("persistent weighted ETA expectation preserves value and gradient", {
  fixture <- estimation_fixture(FALSE)
  context <- .nm_estimation_context(
    fixture$model, fixture$data, method = "IMP"
  )
  map <- .nm_outer_map(context$model)
  parameters <- map$decode(map$start)
  eta <- lapply(seq_len(context$n_subjects), function(subject) {
    matrix(c(-0.08, 0.03, 0.11) + 0.01 * subject, ncol = context$n_eta)
  })
  weights <- rep(list(c(0.2, 0.5, 0.3)), context$n_subjects)
  reference <- .nm_complete_data_expectation(
    context, map, eta, weights, native_context = NULL
  )
  native_pointer <- .nm_weighted_eta_context(context)
  expect_s3_class(native_pointer, "liberation_weighted_eta_context_ptr")
  native <- .nm_complete_data_expectation(
    context, map, eta, weights, native_context = native_pointer
  )

  expect_equal(
    native$objective(parameters), reference$objective(parameters),
    tolerance = 2e-12
  )
  expect_equal(
    native$gradient(parameters), reference$gradient(parameters),
    tolerance = 2e-11
  )
  expect_equal(
    native$native_gradient(parameters), reference$native_gradient(parameters),
    tolerance = 2e-11
  )
  expect_true(native$telemetry()$native)
})

test_that("native SAEM Q recurrence retains exact moments and SIGMA update", {
  fixture <- estimation_fixture(FALSE)
  context <- .nm_estimation_context(
    fixture$model, fixture$data, method = "SAEM"
  )
  q_state <- .nm_saem_q_state(context)
  samples <- list(
    matrix(c(-0.12, 0.03, 0.19), context$n_subjects, context$n_eta),
    matrix(c(-0.04, 0.07, 0.13), context$n_subjects, context$n_eta),
    matrix(c(0.02, 0.05, 0.09), context$n_subjects, context$n_eta)
  )
  q_state$update(samples[[1L]], 1)
  q_state$update(samples[[2L]], 0.4)
  q_state$update(samples[[3L]], 0.25)
  expected_weights <- c(0.45, 0.30, 0.25)
  expected_mean <- Reduce(`+`, Map(`*`, samples, expected_weights))
  expect_equal(q_state$weights(), expected_weights, tolerance = 2e-15)
  expect_equal(q_state$mean_eta(), expected_mean, tolerance = 2e-15)
  growth_state <- .nm_saem_q_state(context)
  growth_state$update(samples[[1L]], 1)
  for (iteration in 2:17) {
    growth_state$update(samples[[3L]] + iteration / 1000, 1 / iteration)
  }
  expect_identical(growth_state$telemetry()$support, 17L)
  expect_lt(
    growth_state$telemetry()$support_reallocations,
    context$n_subjects * 16
  )

  capped_state <- .nm_saem_q_state(context, max_support = 10L)
  capped_state$update(samples[[1L]], 1)
  capped_state$update(samples[[2L]], 0.4)
  capped_state$update(samples[[3L]], 0.25)
  expect_equal(capped_state$mean_eta(), expected_mean, tolerance = 2e-15)

  omega_reference <- Reduce(`+`, Map(function(eta, weight) {
    weight * .nm_saem_omega_sufficient(context, eta)
  }, samples, expected_weights))
  expect_equal(q_state$omega(), omega_reference, tolerance = 2e-14)

  parameters <- list(
    theta = fixture$model$THETAS$Value,
    sigma = fixture$model$SIGMAS$Value,
    omega = fixture$model$OMEGAS$Value
  )
  expectation <- .nm_saem_q_expectation(
    context, .nm_outer_map(context$model), q_state
  )
  gradient_update <- .nm_saem_sigma_from_q_gradient(
    context, parameters, expectation$native_gradient(parameters)
  )
  simulation_update <- q_state$sigma(parameters)
  expect_equal(gradient_update, simulation_update, tolerance = 2e-10)
  expect_identical(q_state$telemetry()$backend, "persistent-cpp-weighted-eta")
})

test_that("native ITS covariance is equivalent to the R reference", {
  fixture <- estimation_fixture()
  context <- .nm_estimation_context(
    fixture$model, fixture$data, method = "ITS"
  )
  parameters <- list(
    theta = fixture$model$THETAS$Value,
    sigma = fixture$model$SIGMAS$Value,
    omega = fixture$model$OMEGAS$Value
  )
  starts <- matrix(0, context$n_subjects, context$n_eta)
  withr::local_options(LibeRation.its_native_covariance = FALSE)
  reference <- .nm_its_distribution(
    context, parameters, starts, eta_maxit = 50L, tolerance = 1e-8
  )
  options(LibeRation.its_native_covariance = TRUE)
  native <- .nm_its_distribution(
    context, parameters, starts, eta_maxit = 50L, tolerance = 1e-8
  )
  expect_identical(
    native$covariance_backend, "cpp-batched-first-order-gaussian"
  )
  expect_equal(native$covariance, reference$covariance, tolerance = 2e-10)
  expect_equal(native$eta, reference$eta, tolerance = 2e-10)
})

test_that("optimized ITS uses the persistent native weighted M-step only", {
  fixture <- estimation_fixture(FALSE)
  optimized <- nm_est(
    fixture$model, fixture$data, method = "ITS", maxit = 2L,
    its_mstep_maxit = 1L, eta_maxit = 40L, collect_output = FALSE,
    numerical_mode = "liber_optimized"
  )
  native <- optimized$diagnostics$native_mstep
  expect_true(native$eligible)
  expect_gt(native$attempts, 0L)
  expect_gt(native$successes, 0L)
  expect_identical(native$fallbacks, 0L)

  compatible <- nm_est(
    fixture$model, fixture$data, method = "ITS", maxit = 2L,
    its_mstep_maxit = 1L, eta_maxit = 40L, collect_output = FALSE,
    numerical_mode = "nonmem_compatibility"
  )
  expect_identical(compatible$diagnostics$native_mstep$attempts, 0L)
  expect_match(
    compatible$diagnostics$native_mstep$compatibility_policy,
    "NONMEM-compatible"
  )
})

test_that("optimized SAEM keeps ODE Q evaluation inside the native M-step", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "DV", "MDV"),
    ADVAN = 6, DOSECMP = 1, OBSCMP = 1,
    PRED = "K=THETA(1)*exp(ETA(1));V=THETA(2);S1=V",
    DES = "DADT(1)=-K*A(1)", ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(
      THETA = 1:2, Value = c(0.1, 20), FIX = c(FALSE, TRUE),
      LOWER = c(0.01, 0.02), UPPER = c(1, 20000)
    ),
    OMEGAS = data.frame(OMEGA = 1, Value = 0.1, FIX = TRUE),
    SIGMAS = data.frame(SIGMA = 1, Value = 0.2, FIX = TRUE)
  )
  subject <- data.frame(
    TIME = c(0, 1, 4), EVID = c(1, 0, 0), AMT = c(100, 0, 0),
    DV = c(NA, 4.52, 3.35), MDV = c(1, 0, 0)
  )
  data <- rbind(transform(subject, ID = 1L), transform(subject, ID = 2L))
  data <- data[, c("ID", setdiff(names(data), "ID"))]
  fit <- nm_est(
    model, data, method = "SAEM", maxit = 1L, n_iter = 2L, burn = 1L,
    mcmc_steps = 1L, mstep_maxit = 1L, saem_kernel = "random_walk",
    collect_output = FALSE, numerical_mode = "liber_optimized", seed = 1941L
  )
  expect_true(fit$diagnostics$native_mstep$eligible)
  expect_gt(fit$diagnostics$native_mstep$attempts, 0L)
  expect_gt(fit$diagnostics$native_mstep$successes, 0L)
  expect_identical(
    fit$diagnostics$native_mstep$objective_backend,
    "native-cpp-weighted-eta-ode-population-objective"
  )
  expect_true(fit$objective_comparable)
  expect_true(all(is.na(fit$diagnostics$replicates$score_errors)))
  expect_match(
    fit$objective_type,
    "importance_sampled_marginal_log_likelihood",
    fixed = TRUE
  )
})

test_that("native IMP weights preserve the R subject calculation", {
  fixture <- estimation_fixture()
  context <- .nm_estimation_context(
    fixture$model, fixture$data, method = "IMP"
  )
  parameters <- list(
    theta = fixture$model$THETAS$Value,
    sigma = fixture$model$SIGMAS$Value,
    omega = fixture$model$OMEGAS$Value
  )
  normals <- .nm_imp_normals(
    context, 12L, 1917L, sampling = "random", proposal = "gaussian"
  )
  withr::local_options(LibeRation.imp_native_expectation = FALSE)
  reference <- .nm_imp_expectation_state(
    context, parameters, normals, 50L, 1e-8
  )
  options(LibeRation.imp_native_expectation = TRUE)
  native <- .nm_imp_expectation_state(
    context, parameters, normals, 50L, 1e-8
  )
  expect_identical(native$backend, "cpp-batched-importance-weights")
  expect_equal(native$weights, reference$weights, tolerance = 2e-12)
  expect_equal(native$ess, reference$ess, tolerance = 2e-11)
})

test_that("native GQ preserves adaptive and fixed finite-grid objectives", {
  fixture <- estimation_fixture(FALSE)
  fixture$model$THETAS$FIX <- c(FALSE, TRUE)
  fixture$model$OMEGAS$FIX <- TRUE
  fixture$model$SIGMAS$FIX <- TRUE
  fixture$model$NUMERICAL_MODE <- "liber_optimized"
  context <- .nm_estimation_context(
    fixture$model, fixture$data, method = "GQ"
  )
  map <- .nm_outer_map(context$model)
  parameters <- map$decode(map$start)
  design <- .nm_gq_design(context, order = 3L)
  nodes <- design$normals[[1L]]
  stochastic <- .nm_stochastic_eta_context(context)
  mu <- .nm_mu_specialization(context, map, enabled = TRUE)
  config <- .nm_bayes_cpp_map_config(context, map, mu)

  for (adaptive in c(TRUE, FALSE)) {
    reference <- .nm_gq_evaluate(
      context, parameters, design$normals, 40L, 1e-7,
      adaptive = adaptive, gradient = TRUE
    )
    pointer <- .liberation_gq_context_create(
      stochastic, config, nodes,
      attr(nodes, "log_measure"), attr(nodes, "measure_sign"),
      adaptive, 40L, 1e-7
    )
    native <- .liberation_gq_context_eval(pointer, map$start, TRUE)
    expect_true(native$valid)
    expect_equal(native$value, reference$value, tolerance = 2e-11)
    expect_equal(
      native$gradient,
      as.vector(reference$native_gradient %*% map$jacobian(parameters)),
      tolerance = 2e-10
    )
    expect_equal(
      native$effective_quadrature_points,
      vapply(
        reference$states, `[[`, numeric(1), "effective_sample_size"
      ),
      tolerance = 2e-11
    )
    expect_equal(
      native$quadrature_cancellation_ratio,
      vapply(reference$states, `[[`, numeric(1), "cancellation_ratio"),
      tolerance = 2e-12
    )
  }

  correlated <- nm_model(
    INPUT = fixture$model$INPUT, ADVAN = 1,
    PRED = paste(
      "CL=THETA(1)*exp(ETA(1))",
      "V=THETA(2)*exp(ETA(2))", "S1=V", sep = ";"
    ),
    ERROR = "Y=F+ERR(1)", THETAS = fixture$model$THETAS,
    OMEGAS = data.frame(
      OMEGA = 1:3, ROW = c(1, 2, 2), COL = c(1, 1, 2),
      Value = c(0.09, 0.02, 0.16), FIX = FALSE
    ),
    SIGMAS = fixture$model$SIGMAS, NUMERICAL_MODE = "liber_optimized"
  )
  correlated_context <- .nm_estimation_context(
    correlated, fixture$data, method = "GQ"
  )
  correlated_map <- .nm_outer_map(correlated_context$model)
  correlated_parameters <- correlated_map$decode(correlated_map$start)
  correlated_design <- .nm_gq_design(correlated_context, order = 1L)
  correlated_nodes <- correlated_design$normals[[1L]]
  correlated_reference <- .nm_gq_evaluate(
    correlated_context, correlated_parameters, correlated_design$normals,
    40L, 1e-7, adaptive = TRUE, gradient = TRUE
  )
  correlated_stochastic <- .nm_stochastic_eta_context(correlated_context)
  correlated_pointer <- .liberation_gq_context_create(
    correlated_stochastic,
    .nm_bayes_cpp_map_config(
      correlated_context, correlated_map,
      .nm_mu_specialization(correlated_context, correlated_map, enabled = TRUE)
    ),
    correlated_nodes, attr(correlated_nodes, "log_measure"),
    attr(correlated_nodes, "measure_sign"), TRUE, 40L, 1e-7
  )
  correlated_native <- .liberation_gq_context_eval(
    correlated_pointer, correlated_map$start, TRUE
  )
  correlated_expected_gradient <- as.vector(
    correlated_reference$native_gradient %*%
      correlated_map$jacobian(correlated_parameters)
  )
  # Both paths solve adaptive conditional modes to tolerance before evaluating
  # the same finite quadrature grid.  Apple Accelerate and the reference BLAS
  # can legitimately finish on adjacent inner-solver iterates, so assert an
  # explicit absolute error bound rather than near-bitwise equality.
  expect_lt(
    max(abs(correlated_native$gradient - correlated_expected_gradient)),
    5e-5
  )

  optimized <- nm_est(
    fixture$model, fixture$data, method = "GQ", maxit = 2L,
    eta_maxit = 40L, gq_order = 3L, collect_output = FALSE,
    numerical_mode = "liber_optimized"
  )
  expect_true(optimized$diagnostics$native_gq_coordinator$used)
  expect_match(
    optimized$diagnostics$optimizer$coordinator,
    "persistent-native-cpp-gq", fixed = TRUE
  )
  expect_true(optimized$diagnostics$exact_finite_grid_refinement)

  compatible <- nm_est(
    fixture$model, fixture$data, method = "GQ", maxit = 1L,
    eta_maxit = 40L, gq_order = 3L, collect_output = FALSE,
    numerical_mode = "nonmem_compatibility"
  )
  expect_false(compatible$diagnostics$native_gq_coordinator$used)
  expect_match(
    compatible$diagnostics$native_gq_coordinator$fallback_reason,
    "NONMEM-compatible"
  )
})

test_that("native GQ preserves signed sparse grids and ODE retaping", {
  fixture <- estimation_fixture(FALSE)
  sparse_model <- nm_model(
    INPUT = fixture$model$INPUT, ADVAN = 1,
    PRED = paste(
      "CL=THETA(1)*exp(ETA(1)+0.2*ETA(2)+0.1*ETA(3)+0.05*ETA(4))",
      "V=THETA(2)", "S1=V", sep = ";"
    ),
    ERROR = "Y=F+ERR(1)",
    THETAS = transform(fixture$model$THETAS, FIX = TRUE),
    OMEGAS = data.frame(
      OMEGA = 1:4, Value = c(0.09, 0.07, 0.05, 0.03), FIX = TRUE
    ),
    SIGMAS = transform(fixture$model$SIGMAS, FIX = TRUE),
    NUMERICAL_MODE = "liber_optimized"
  )
  sparse_context <- .nm_estimation_context(
    sparse_model, fixture$data, method = "GQ"
  )
  sparse_map <- .nm_outer_map(sparse_context$model)
  sparse_parameters <- sparse_map$decode(sparse_map$start)
  sparse_design <- .nm_gq_design(sparse_context, order = 3L)
  expect_true(any(attr(sparse_design$normals[[1L]], "measure_sign") < 0))
  sparse_reference <- .nm_gq_evaluate(
    sparse_context, sparse_parameters, sparse_design$normals,
    30L, 1e-7, adaptive = TRUE, gradient = TRUE
  )
  sparse_nodes <- sparse_design$normals[[1L]]
  sparse_pointer <- .liberation_gq_context_create(
    .nm_stochastic_eta_context(sparse_context),
    .nm_bayes_cpp_map_config(
      sparse_context, sparse_map,
      .nm_mu_specialization(sparse_context, sparse_map, enabled = TRUE)
    ),
    sparse_nodes, attr(sparse_nodes, "log_measure"),
    attr(sparse_nodes, "measure_sign"), TRUE, 30L, 1e-7
  )
  sparse_native <- .liberation_gq_context_eval(
    sparse_pointer, sparse_map$start, TRUE
  )
  expect_true(sparse_native$valid)
  expect_equal(sparse_native$value, sparse_reference$value, tolerance = 2e-10)
  expect_equal(
    sparse_native$quadrature_cancellation_ratio,
    vapply(
      sparse_reference$states, `[[`, numeric(1), "cancellation_ratio"
    ),
    tolerance = 2e-11
  )

  ode_model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "DV", "MDV"),
    ADVAN = 6, DOSECMP = 1, OBSCMP = 1,
    PRED = "K=THETA(1)*exp(ETA(1));V=THETA(2);S1=V",
    DES = "DADT(1)=-K*A(1)", ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(
      THETA = 1:2, Value = c(0.1, 20), FIX = c(FALSE, TRUE),
      LOWER = c(0.01, 1), UPPER = c(1, 100)
    ),
    OMEGAS = data.frame(OMEGA = 1, Value = 0.1, FIX = TRUE),
    SIGMAS = data.frame(SIGMA = 1, Value = 0.2, FIX = TRUE),
    NUMERICAL_MODE = "liber_optimized"
  )
  subject <- data.frame(
    TIME = c(0, 1, 4), EVID = c(1, 0, 0), AMT = c(100, 0, 0),
    DV = c(NA, 4.52, 3.35), MDV = c(1, 0, 0)
  )
  ode_data <- rbind(
    transform(subject, ID = 1L), transform(subject, ID = 2L)
  )
  ode_data <- ode_data[, c("ID", setdiff(names(ode_data), "ID"))]
  ode_fit <- nm_est(
    ode_model, ode_data, method = "GQ", maxit = 1L,
    eta_maxit = 30L, gq_order = 1L, gq_adaptive = FALSE,
    collect_output = FALSE,
    numerical_mode = "liber_optimized"
  )
  expect_true(ode_fit$diagnostics$native_gq_coordinator$used)
  expect_true(is.finite(ode_fit$objective))
  expect_gt(
    ode_fit$diagnostics$native_gq_coordinator$telemetry$
      stochastic_context$tape_records,
    0
  )
})

test_that("defensive IMP uses its realised deterministic-mixture allocation", {
  fixture <- estimation_fixture()
  context <- .nm_estimation_context(
    fixture$model, fixture$data, method = "IMP"
  )
  parameters <- list(
    theta = fixture$model$THETAS$Value,
    sigma = fixture$model$SIGMAS$Value,
    omega = fixture$model$OMEGAS$Value
  )
  normals <- .nm_imp_normals(
    context, 9L, seed = 1920L,
    sampling = "random", proposal = "defensive"
  )
  expect_equal(
    attr(normals[[1L]], "imp_defensive_gaussian_weight"), 5 / 9,
    tolerance = 0
  )
  mode <- .nm_subject_modes(
    context, parameters, maxit = 50L, tolerance = 1e-8,
    exact_hessian = TRUE
  )[[1L]]
  proposal <- .nm_imp_proposal_from_mode(
    context$subjects[[1L]], parameters, normals[[1L]], mode
  )
  covariance <- .nm_positive_definite(
    2 * solve(mode$hessian), "test defensive covariance"
  )$matrix
  logdet <- as.numeric(determinant(covariance, logarithm = TRUE)$modulus)
  z <- normals[[1L]]
  dimension <- ncol(z)
  gaussian <- -0.5 * (
    dimension * log(2 * pi) + logdet + rowSums(z^2)
  )
  df <- attr(z, "imp_proposal_df")
  scale <- (df - 2) / df
  student <- lgamma((df + dimension) / 2) - lgamma(df / 2) -
    dimension * log(df * pi) / 2 -
    (logdet + dimension * log(scale)) / 2 -
    (df + dimension) * log1p(rowSums(z^2) / scale / df) / 2
  maximum <- pmax(gaussian, student)
  expected <- maximum + log(
    (5 / 9) * exp(gaussian - maximum) +
      (4 / 9) * exp(student - maximum)
  )
  expect_equal(proposal$log_proposal, expected, tolerance = 2e-14)
})

test_that("optimized stochastic policy extensions do not alter compatibility defaults", {
  fixture <- estimation_fixture()
  compatible_imp <- nm_est(
    fixture$model, fixture$data, method = "IMP", maxit = 2L,
    n_imp = 10L, imp_mstep_maxit = 1L, eta_maxit = 40L,
    collect_output = FALSE, numerical_mode = "nonmem_compatibility",
    seed = 1921L
  )
  expect_identical(compatible_imp$diagnostics$sampling, "random")
  expect_identical(compatible_imp$diagnostics$proposal, "gaussian")
  expect_identical(compatible_imp$diagnostics$proposal_curvature, "exact")
  expect_identical(compatible_imp$diagnostics$mstep_schedule, "fixed")
  expect_false(compatible_imp$diagnostics$proposal_mode_reuse$enabled)

  optimized_imp <- nm_est(
    fixture$model, fixture$data, method = "IMP", maxit = 2L,
    n_imp = 10L, imp_mstep_maxit = 1L, eta_maxit = 40L,
    collect_output = FALSE, numerical_mode = "liber_optimized",
    seed = 1921L
  )
  expect_identical(optimized_imp$diagnostics$sampling, "antithetic")
  expect_identical(optimized_imp$diagnostics$proposal, "defensive")
  expect_identical(optimized_imp$diagnostics$proposal_curvature, "fisher")
  expect_identical(optimized_imp$diagnostics$mstep_schedule, "progressive")
  expect_true(optimized_imp$diagnostics$proposal_mode_reuse$enabled)
  expect_true(is.list(
    optimized_imp$diagnostics$native_weighted_mstep
  ))
  expect_match(
    optimized_imp$diagnostics$expectation_backend,
    "persistent-native-importance"
  )
  expect_gt(
    optimized_imp$diagnostics$weighted_expectation$weighted_context$
      importance_updates,
    0L
  )
  expect_false(
    optimized_imp$diagnostics$weighted_expectation$weighted_context$
      reduced_population_requested
  )

  fixed_imp <- nm_est(
    fixture$model, fixture$data, method = "IMP", maxit = 2L,
    n_imp = 10L, imp_mstep_maxit = 1L, eta_maxit = 40L,
    imp_sample_schedule = "fixed", imp_subject_allocation = "fixed",
    collect_output = FALSE, numerical_mode = "liber_optimized",
    seed = 1921L
  )
  expect_true(
    fixed_imp$diagnostics$weighted_expectation$weighted_context$
      reduced_population_requested
  )
  expect_gte(
    fixed_imp$diagnostics$weighted_expectation$weighted_context$
      reduced_population_records,
    1
  )

  compatible_saem <- nm_est(
    fixture$model, fixture$data, method = "SAEM", maxit = 2L,
    n_iter = 3L, burn = 1L, mcmc_steps = 1L, mstep_maxit = 1L,
    collect_output = FALSE, numerical_mode = "nonmem_compatibility",
    seed = 1922L
  )
  expect_false(compatible_saem$diagnostics$stationarity$auto_stop)
  expect_identical(
    compatible_saem$diagnostics$parameter_averaging$resolved, "none"
  )
  expect_true(all(compatible_saem$diagnostics$mstep_schedule$performed))
  expect_gt(
    compatible_saem$diagnostics$sigma_sufficient_statistics$
      gradient_recovery_updates,
    0L
  )
})

test_that("optimized IMP retains its weighted M-step in the native context", {
  fixture <- estimation_fixture(FALSE)
  fit <- nm_est(
    fixture$model, fixture$data, method = "IMP", maxit = 2L,
    n_imp = 10L, imp_mstep_maxit = 1L, eta_maxit = 40L,
    collect_output = FALSE, numerical_mode = "liber_optimized",
    seed = 1926L
  )
  native <- fit$diagnostics$native_weighted_mstep
  expect_gt(native$attempts, 0L)
  expect_gt(native$successes, 0L)
  expect_true(native$persistent_optimizer_state)
  expect_gt(native$sigma_gradient_updates, 0L)
  expect_gt(native$omega_moment_updates, 0L)
})

test_that("weighted CppAD subject workers are deterministic", {
  fixture <- estimation_fixture(FALSE)
  fixture$model$NUMERICAL_MODE <- "liber_optimized"
  serial_context <- .nm_estimation_context(
    fixture$model, fixture$data, n_cores = 1L, method = "IMP"
  )
  threaded_context <- .nm_estimation_context(
    fixture$model, fixture$data, n_cores = 2L, method = "IMP"
  )
  eta <- lapply(seq_len(serial_context$n_subjects), function(subject) {
    matrix(c(-0.08, 0.03, 0.11) + subject / 100, ncol = 1L)
  })
  weights <- rep(list(c(0.2, 0.5, 0.3)), serial_context$n_subjects)
  serial <- .nm_weighted_eta_context(
    serial_context, reduced_population_tape = FALSE
  )
  threaded <- .nm_weighted_eta_context(
    threaded_context, reduced_population_tape = FALSE
  )
  .liberation_weighted_eta_context_set(serial, eta, weights)
  .liberation_weighted_eta_context_set(threaded, eta, weights)
  theta <- fixture$model$THETAS$Value
  sigma <- fixture$model$SIGMAS$Value
  omega <- fixture$model$OMEGAS$Value
  serial_value <- .liberation_weighted_eta_context_eval(
    serial, theta, sigma, omega
  )
  threaded_value <- .liberation_weighted_eta_context_eval(
    threaded, theta, sigma, omega
  )
  repeated_threaded_value <- .liberation_weighted_eta_context_eval(
    threaded, theta, sigma, omega
  )
  telemetry <- .liberation_weighted_eta_context_telemetry(threaded)

  expect_equal(threaded_value$value, serial_value$value, tolerance = 2e-12)
  expect_equal(
    threaded_value$gradient, serial_value$gradient, tolerance = 2e-11
  )
  expect_identical(
    unname(repeated_threaded_value$value), unname(threaded_value$value)
  )
  expect_identical(
    unname(repeated_threaded_value$gradient),
    unname(threaded_value$gradient)
  )
  expect_identical(telemetry$native_subject_threads, 2L)
  expect_true(telemetry$native_subject_parallel)
  expect_gt(telemetry$cppad_subject_dispatches, 0)
})

test_that("optimized IMP native M-step supports MU referencing", {
  fixture <- estimation_fixture(fix = FALSE)
  theta <- fixture$model$THETAS
  theta$FIX <- c(FALSE, TRUE)
  model <- nm_model(
    INPUT = fixture$model$INPUT, ADVAN = 1,
    PRED = "V=THETA(2); S1=V", ERROR = "Y=F+ERR(1)",
    THETAS = theta, OMEGAS = fixture$model$OMEGAS,
    SIGMAS = fixture$model$SIGMAS,
    MU = nm_mu(1, "log(THETA(1))", parameter = "CL")
  )
  fit <- nm_est(
    model, fixture$data, method = "IMP", maxit = 2L,
    n_imp = 10L, imp_mstep_maxit = 1L, eta_maxit = 30L,
    collect_output = FALSE, numerical_mode = "liber_optimized", seed = 1928L
  )
  native <- fit$diagnostics$native_weighted_mstep
  expect_true(fit$diagnostics$mu_specialization$active)
  expect_gt(native$attempts, 0L)
  expect_gt(native$successes, 0L)
  expect_null(native$fallback_reason)
})

test_that("weighted ODE objectives retape on the native IMP path", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "DV", "MDV"),
    ADVAN = 6, DOSECMP = 1, OBSCMP = 1,
    PRED = "K=THETA(1)*exp(ETA(1)); V=THETA(2); S1=V",
    DES = "DADT(1)=-K*A(1)", ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(
      THETA = 1:2, Value = c(0.1, 20), FIX = c(FALSE, TRUE),
      LOWER = c(0.01, 1), UPPER = c(1, 100)
    ),
    OMEGAS = data.frame(OMEGA = 1, Value = 0.1, FIX = TRUE),
    SIGMAS = data.frame(SIGMA = 1, Value = 0.2, FIX = TRUE),
    NUMERICAL_MODE = "liber_optimized"
  )
  subject <- data.frame(
    TIME = c(0, 1, 4), EVID = c(1, 0, 0), AMT = c(100, 0, 0),
    DV = c(NA, 4.52, 3.35), MDV = c(1, 0, 0)
  )
  data <- rbind(transform(subject, ID = 1L), transform(subject, ID = 2L))
  data <- data[, c("ID", setdiff(names(data), "ID"))]
  context <- .nm_estimation_context(model, data, n_cores = 2L, method = "IMP")
  native <- .nm_weighted_eta_context(
    context, reduced_population_tape = FALSE
  )
  eta <- rep(list(matrix(c(-0.25, 0.35), ncol = 1L)), 2L)
  weights <- rep(list(c(0.5, 0.5)), 2L)
  .liberation_weighted_eta_context_set(native, eta, weights)
  value <- .liberation_weighted_eta_context_eval(
    native, model$THETAS$Value, model$SIGMAS$Value, model$OMEGAS$Value
  )
  telemetry <- .liberation_weighted_eta_context_telemetry(native)
  expect_true(is.finite(value$value))
  expect_true(all(is.finite(value$gradient)))
  expect_identical(telemetry$native_subject_threads, 2L)
  expect_true(telemetry$native_subject_parallel)
  expect_gte(telemetry$tape_records, 4)
  fit <- nm_est(
    model, data, method = "IMP", maxit = 2L, n_imp = 10L,
    imp_mstep_maxit = 1L, eta_maxit = 20L, n_cores = 2L,
    collect_output = FALSE, numerical_mode = "liber_optimized", seed = 1929L
  )
  expect_gt(fit$diagnostics$native_weighted_mstep$attempts, 0L)
  expect_gt(fit$diagnostics$native_weighted_mstep$successes, 0L)
})

test_that("optimized randomized QMC importance draws are seeded and finite", {
  context <- list(n_subjects = 2L, n_eta = 2L)
  first <- .nm_imp_normals(
    context, c(11L, 17L), 913L, sampling = "rqmc", proposal = "gaussian"
  )
  second <- .nm_imp_normals(
    context, c(11L, 17L), 913L, sampling = "rqmc", proposal = "gaussian"
  )
  expect_equal(first, second, tolerance = 0)
  expect_true(all(vapply(first, function(draws) all(is.finite(draws)), logical(1))))
  expect_identical(attr(first[[1L]], "imp_sampling"), "rqmc")
})

test_that("native SAEM execution preserves the compatibility trajectory", {
  fixture <- estimation_fixture(FALSE)
  controls <- list(
    model = fixture$model, data = fixture$data, method = "SAEM",
    maxit = 3L, n_iter = 5L, burn = 2L, mcmc_steps = 1L,
    mstep_maxit = 2L, collect_output = FALSE,
    numerical_mode = "nonmem_compatibility", seed = 1927L,
    optimizer_backend = "r"
  )
  previous <- options(LibeRation.weighted_eta_context = FALSE)
  on.exit(options(previous), add = TRUE)
  reference <- do.call(nm_est, controls)
  options(LibeRation.weighted_eta_context = TRUE)
  native <- do.call(nm_est, controls)

  expect_equal(native$theta, reference$theta, tolerance = 2e-10)
  expect_equal(native$sigma, reference$sigma, tolerance = 2e-10)
  expect_equal(native$omega, reference$omega, tolerance = 2e-10)
  expect_equal(native$objective, reference$objective, tolerance = 2e-10)
  expect_equal(
    native$diagnostics$objective_trace,
    reference$diagnostics$objective_trace,
    tolerance = 2e-10
  )
  expect_identical(
    native$diagnostics$stochastic_approximation$backend$backend,
    "persistent-cpp-weighted-eta"
  )
})
