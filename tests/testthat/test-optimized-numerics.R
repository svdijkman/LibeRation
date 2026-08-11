test_that("optimized residual simulation preserves the compatibility stream", {
  fixture <- estimation_fixture(FALSE)
  compatible <- nm_simulate(
    fixture$model, fixture$data, nsim = 3L, random_effects = TRUE,
    residual = TRUE, seed = 741L,
    numerical_mode = "nonmem_compatibility"
  )
  optimized <- nm_simulate(
    fixture$model, fixture$data, nsim = 3L, random_effects = TRUE,
    residual = TRUE, seed = 741L, numerical_mode = "liber_optimized"
  )
  expect_identical(names(optimized), names(compatible))
  attr(compatible, "numerical_mode") <- NULL
  attr(optimized, "numerical_mode") <- NULL
  expect_equal(as.data.frame(optimized), as.data.frame(compatible), tolerance = 0)

  compatible_batch <- nm_simulate(
    fixture$model, fixture$data, nsim = 4L, random_effects = TRUE,
    residual = FALSE, seed = 742L,
    numerical_mode = "nonmem_compatibility"
  )
  optimized_batch <- nm_simulate(
    fixture$model, fixture$data, nsim = 4L, random_effects = TRUE,
    residual = FALSE, seed = 742L, numerical_mode = "liber_optimized"
  )
  attr(compatible_batch, "numerical_mode") <- NULL
  attr(optimized_batch, "numerical_mode") <- NULL
  expect_equal(
    as.data.frame(optimized_batch), as.data.frame(compatible_batch),
    tolerance = 0
  )
})

test_that("compatibility mode safely reuses deterministic simulation work", {
  fixture <- estimation_fixture()
  model <- nm_model_update(fixture$model, OUTPUT = c("PRED", "IPRED"))
  eta <- matrix(c(-0.15, 0.05, 0.2), ncol = 1L)
  single <- nm_simulate(
    model, fixture$data, eta = eta, nsim = 1L, residual = FALSE,
    numerical_mode = "nonmem_compatibility"
  )
  repeated <- nm_simulate(
    model, fixture$data, eta = eta, nsim = 4L, residual = FALSE,
    numerical_mode = "nonmem_compatibility"
  )
  reference <- single[, names(single), drop = FALSE]
  rownames(reference) <- NULL
  for (replicate in seq_len(4L)) {
    observed <- repeated[repeated$SIM == replicate,
                         setdiff(names(repeated), "SIM"), drop = FALSE]
    rownames(observed) <- NULL
    expect_equal(observed, reference, tolerance = 0)
  }
})

test_that("identical proposal covariance roots are policy neutral", {
  fixture <- estimation_fixture()
  context <- .nm_estimation_context(fixture$model, fixture$data)
  omega <- fixture$model$OMEGAS$Value
  cache <- .nm_proposal_root_cache(context)
  cached <- .nm_proposal_roots(context, omega, cache = cache)
  direct <- lapply(context$subjects, function(evaluator) {
    t(chol(.nm_effect_covariance(context$model, evaluator$data, omega)))
  })
  expect_identical(cached, direct)
  expect_identical(
    .nm_proposal_roots(context, omega, cache = cache), cached
  )
  expect_equal(cache$factorizations, length(unique(cache$groups)))
  expect_equal(cache$misses, 1L)
  expect_equal(cache$hits, 1L)
})

test_that("conditional ETA Metropolis reuses exact current subject values", {
  fixture <- estimation_fixture()
  context <- .nm_estimation_context(fixture$model, fixture$data)
  parameters <- list(
    theta = fixture$model$THETAS$Value,
    sigma = fixture$model$SIGMAS$Value,
    omega = fixture$model$OMEGAS$Value
  )
  eta <- matrix(0, context$n_subjects, context$n_eta)
  roots <- .nm_proposal_roots(context, parameters$omega)
  set.seed(8461L)
  mcmc_steps <- 2L
  normals <- matrix(
    stats::rnorm(context$n_subjects * mcmc_steps * context$n_eta),
    context$n_subjects * mcmc_steps, context$n_eta
  )
  log_uniforms <- log(stats::runif(context$n_subjects * mcmc_steps))
  values <- .nm_objective_collection(context$subjects, parameters, eta)
  ordinary <- .nm_saem_metropolis_chunk(
    context$subjects, parameters, eta, roots, normals, log_uniforms,
    mcmc_steps, 0.35
  )
  cached <- .nm_saem_metropolis_chunk(
    context$subjects, parameters, eta, roots, normals, log_uniforms,
    mcmc_steps, 0.35, values
  )
  expect_equal(cached$eta, ordinary$eta, tolerance = 0)
  expect_equal(cached$value, ordinary$value, tolerance = 0)
  expect_identical(cached$accepted, ordinary$accepted)
  expect_equal(cached$current_evaluations, 0L)
  expect_equal(cached$current_cache_hits, context$n_subjects)
  expect_equal(ordinary$current_evaluations, context$n_subjects)
})

test_that("persistent stochastic context preserves the random-walk kernel", {
  fixture <- estimation_fixture()
  fixture$model$NUMERICAL_MODE <- "liber_optimized"
  context <- .nm_estimation_context(fixture$model, fixture$data)
  parameters <- list(
    theta = fixture$model$THETAS$Value,
    sigma = fixture$model$SIGMAS$Value,
    omega = fixture$model$OMEGAS$Value
  )
  eta <- matrix(c(-0.08, 0.03, 0.11), context$n_subjects, context$n_eta)
  roots <- .nm_proposal_roots(context, parameters$omega)
  steps <- 3L
  set.seed(8462L)
  normals <- matrix(
    stats::rnorm(context$n_subjects * steps * context$n_eta),
    context$n_subjects * steps, context$n_eta
  )
  log_uniforms <- log(stats::runif(context$n_subjects * steps))
  current <- .nm_objective_collection(context$subjects, parameters, eta)
  reference <- .nm_saem_metropolis_chunk(
    context$subjects, parameters, eta, roots, normals, log_uniforms,
    steps, 0.35, current
  )
  persistent <- .nm_stochastic_eta_context(context)
  expect_s3_class(persistent, "liberation_stochastic_eta_context_ptr")
  expect_equal(
    .liberation_stochastic_eta_context_eval(
      persistent, parameters$theta, eta, parameters$sigma, parameters$omega
    ),
    current, tolerance = 2e-12
  )
  observed <- .liberation_stochastic_eta_context_random_walk(
    persistent, parameters$theta, eta, parameters$sigma, parameters$omega,
    roots, normals, log_uniforms, steps, 0.35, current
  )
  expect_equal(observed$eta, reference$eta, tolerance = 0)
  expect_equal(observed$value, reference$value, tolerance = 2e-12)
  expect_identical(observed$accepted, reference$accepted)
  expect_equal(observed$current_evaluations, 0L)
  expect_equal(observed$current_cache_hits, context$n_subjects)
})

test_that("Laplace independence kernel applies the proposal-density ratio", {
  fixture <- estimation_fixture()
  fixture$model <- nm_model_update(
    fixture$model,
    PRED = paste(
      "CL=THETA(1)*exp(ETA(1));",
      "V=THETA(2)*exp(ETA(2)); S1=V"
    ),
    OMEGAS = data.frame(
      OMEGA = 1:2, Value = c(0.09, 0.04), FIX = TRUE
    ),
    NUMERICAL_MODE = "liber_optimized"
  )
  context <- .nm_estimation_context(fixture$model, fixture$data)
  parameters <- list(
    theta = fixture$model$THETAS$Value,
    sigma = fixture$model$SIGMAS$Value,
    omega = fixture$model$OMEGAS$Value
  )
  eta <- matrix(
    c(-0.08, 0.03, 0.11, -0.04, 0.02, 0.07),
    context$n_subjects, context$n_eta, byrow = TRUE
  )
  modes <- matrix(0, context$n_subjects, context$n_eta)
  roots <- replicate(
    context$n_subjects, diag(c(0.22, 0.14)), simplify = FALSE
  )
  precisions <- lapply(roots, function(root) {
    solve(root %*% t(root))
  })
  set.seed(8463L)
  normals <- matrix(
    stats::rnorm(context$n_subjects * context$n_eta),
    context$n_subjects, context$n_eta
  )
  log_uniforms <- log(stats::runif(context$n_subjects))
  persistent <- .nm_stochastic_eta_context(context)
  current <- .liberation_stochastic_eta_context_eval(
    persistent, parameters$theta, eta, parameters$sigma, parameters$omega
  )
  expected_eta <- eta
  expected_value <- current
  accepted <- 0L
  for (subject in seq_len(context$n_subjects)) {
    candidate <- as.vector(modes[subject, ] +
      roots[[subject]] %*% normals[subject, ])
    candidate_eta <- expected_eta
    candidate_eta[subject, ] <- candidate
    candidate_value <- .liberation_stochastic_eta_context_eval(
      persistent, parameters$theta, candidate_eta,
      parameters$sigma, parameters$omega
    )[[subject]]
    centered_current <- expected_eta[subject, ] - modes[subject, ]
    centered_candidate <- candidate - modes[subject, ]
    log_ratio <- -0.5 * (candidate_value - expected_value[[subject]]) +
      0.5 * (
        drop(crossprod(centered_candidate,
                       precisions[[subject]] %*% centered_candidate)) -
        drop(crossprod(centered_current,
                       precisions[[subject]] %*% centered_current))
      )
    if (log_uniforms[[subject]] < log_ratio) {
      expected_eta[subject, ] <- candidate
      expected_value[[subject]] <- candidate_value
      accepted <- accepted + 1L
    }
  }
  observed <- .liberation_stochastic_eta_context_independence(
    persistent, parameters$theta, eta, parameters$sigma, parameters$omega,
    modes, roots, precisions, normals, log_uniforms, 1L, current
  )
  expect_equal(observed$eta, expected_eta, tolerance = 0)
  expect_equal(observed$value, expected_value, tolerance = 0)
  expect_equal(observed$accepted, accepted)
  expect_match(observed$kernel, "laplace-independence$")

  fallback <- .nm_saem_independence_chunk(
    context$subjects, parameters, eta, modes, roots, precisions,
    normals, log_uniforms, 1L, current
  )
  expect_equal(fallback$eta, observed$eta, tolerance = 1e-14)
  expect_equal(fallback$value, observed$value, tolerance = 1e-14)
  expect_identical(fallback$accepted, observed$accepted)
  expect_identical(fallback$kernel, "laplace-independence")
})

test_that("persistent f-SAEM proposal reuses optimized Laplace mode machinery", {
  fixture <- estimation_fixture()
  fixture$model <- nm_model_update(
    fixture$model,
    PRED = paste(
      "CL=THETA(1)*exp(ETA(1));",
      "V=THETA(2)*exp(ETA(2)); S1=V"
    ),
    OMEGAS = data.frame(
      OMEGA = 1:2, Value = c(0.09, 0.04), FIX = TRUE
    ),
    NUMERICAL_MODE = "liber_optimized"
  )
  context <- .nm_estimation_context(fixture$model, fixture$data)
  parameters <- list(
    theta = fixture$model$THETAS$Value,
    sigma = fixture$model$SIGMAS$Value,
    omega = fixture$model$OMEGAS$Value
  )
  persistent <- .nm_stochastic_eta_context(context)
  proposal <- .nm_fsaem_proposal(
    context, parameters,
    matrix(0, context$n_subjects, context$n_eta),
    eta_maxit = 30L, tolerance = 1e-6, persistent = persistent
  )
  expect_identical(proposal$backend, "persistent-cpp-laplace-proposal")
  expect_equal(dim(proposal$modes), c(context$n_subjects, context$n_eta))
  expect_true(all(is.finite(proposal$modes)))
  expect_true(all(vapply(proposal$roots, function(root) {
    identical(dim(root), c(context$n_eta, context$n_eta)) &&
      all(is.finite(root))
  }, logical(1))))
  expect_equal(
    proposal$anchor,
    c(parameters$theta, parameters$sigma, parameters$omega), tolerance = 0
  )
  telemetry <- .liberation_stochastic_eta_context_telemetry(persistent)
  expect_equal(telemetry$laplace_refreshes, 1L)
  expect_gte(telemetry$laplace_mode_evaluations, context$n_subjects)
})

test_that("native nonparametric EM agrees with the compatibility implementation", {
  loglik <- matrix(c(
    -0.2, -1.1, -2.0,
    -1.4, -0.1, -0.8,
    -2.1, -0.7, -0.2,
    -0.6, -0.3, -1.5
  ), nrow = 4L, byrow = TRUE)
  initial <- c(0.2, 0.5, 0.3)
  compatible <- .nm_np_weights(
    loglik, initial, maxit = 200L, tolerance = 1e-12,
    optimized = FALSE
  )
  optimized <- .nm_np_weights(
    loglik, initial, maxit = 200L, tolerance = 1e-12,
    optimized = TRUE
  )
  expect_equal(optimized$weights, compatible$weights, tolerance = 1e-6)
  expect_equal(
    optimized$responsibilities, compatible$responsibilities,
    tolerance = 1e-6
  )
  expect_equal(
    optimized$log_likelihood, compatible$log_likelihood,
    tolerance = 1e-6
  )
  expect_identical(optimized$solver, "primal-dual-barrier-newton")
  expect_gte(
    optimized$log_likelihood,
    compatible$log_likelihood - 1e-7
  )

  fixture <- estimation_fixture()
  parameters <- list(
    theta = fixture$model$THETAS$Value,
    sigma = fixture$model$SIGMAS$Value,
    omega = fixture$model$OMEGAS$Value
  )
  compatible_context <- .nm_estimation_context(fixture$model, fixture$data)
  optimized_context <- .nm_estimation_context(
    nm_model_update(fixture$model, NUMERICAL_MODE = "liber_optimized"),
    fixture$data
  )
  supports <- matrix(c(-0.4, -0.1, 0.2, 0.5), ncol = 1L)
  compatible_grid <- .nm_np_loglik(
    compatible_context, parameters, supports, gradient = TRUE
  )
  optimized_grid <- .nm_np_loglik(
    optimized_context, parameters, supports, gradient = TRUE
  )
  expect_equal(optimized_grid$loglik, compatible_grid$loglik, tolerance = 1e-13)
  expect_equal(optimized_grid$gradient, compatible_grid$gradient, tolerance = 1e-12)
})

test_that("batched importance collection retains subject values and gradients", {
  fixture <- estimation_fixture()
  optimized_model <- nm_model_update(
    fixture$model, NUMERICAL_MODE = "liber_optimized"
  )
  context <- .nm_estimation_context(optimized_model, fixture$data)
  parameters <- list(
    theta = fixture$model$THETAS$Value,
    sigma = fixture$model$SIGMAS$Value,
    omega = fixture$model$OMEGAS$Value
  )
  set.seed(194L)
  normals <- lapply(seq_len(context$n_subjects), function(index) {
    matrix(stats::rnorm(17L * context$n_eta), 17L, context$n_eta)
  })
  proposals <- .nm_imp_prepare_proposals(
    context, parameters, normals, eta_maxit = 80L,
    tolerance = 1e-9, adaptive = TRUE
  )
  optimized <- .nm_imp_evaluate_fixed(
    context, parameters, proposals, gradient = TRUE
  )
  compatible_states <- Map(function(evaluator, proposal) {
    .nm_imp_subject_from_proposal(
      evaluator, parameters, proposal, gradient = TRUE
    )
  }, context$subjects, proposals)
  expect_equal(
    optimized$value,
    sum(vapply(compatible_states, `[[`, numeric(1), "value")),
    tolerance = 1e-11
  )
  expect_equal(
    optimized$native_gradient,
    Reduce(`+`, lapply(compatible_states, `[[`, "native_gradient")),
    tolerance = 1e-10
  )
  expect_equal(
    vapply(optimized$states, `[[`, numeric(1), "effective_sample_size"),
    vapply(compatible_states, `[[`, numeric(1), "effective_sample_size"),
    tolerance = 1e-12
  )
})

test_that("native SAEM M-step matches the fixed-ETA objective and gradient", {
  fixture <- estimation_fixture(FALSE)
  fixture$model$NUMERICAL_MODE <- "liber_optimized"
  context <- .nm_estimation_context(
    fixture$model, fixture$data, method = "SAEM"
  )
  parameters <- list(
    theta = fixture$model$THETAS$Value,
    sigma = fixture$model$SIGMAS$Value,
    omega = fixture$model$OMEGAS$Value
  )
  mstep_model <- fixture$model
  mstep_model$OMEGAS$FIX[] <- TRUE
  map <- .nm_outer_map(mstep_model)
  eta <- matrix(
    rep(c(0.04, -0.03), length.out = context$n_subjects * context$n_eta),
    context$n_subjects, context$n_eta
  )
  native <- .nm_saem_native_mstep(
    context, parameters, eta, map, maxit = 2L,
    tolerance = 1e-8, trace = 0L
  )
  candidate <- list(
    theta = native$theta, sigma = native$sigma, omega = native$omega
  )
  expect_equal(
    native$value,
    .nm_saem_conditional(context, candidate, eta),
    tolerance = 1e-10
  )
  expect_equal(
    native$gradient,
    .nm_saem_conditional_gradient(context, map, candidate, eta),
    tolerance = 1e-9
  )
  expect_equal(
    native$native_gradient,
    colSums(.nm_objective_collection_gradient(
      context$subjects, candidate, eta, interaction = TRUE
    )),
    tolerance = 1e-9
  )
})

test_that("FO shares eligible conditional tapes without changing posthoc modes", {
  fixture <- estimation_fixture(FALSE)
  .nm_estimation_context_cache_clear()
  shared <- .nm_estimation_context(fixture$model, fixture$data, method = "FO")
  ordinary <- .nm_estimation_context_build(
    fixture$model, fixture$data, method = NULL
  )
  parameters <- list(
    theta = fixture$model$THETAS$Value,
    sigma = fixture$model$SIGMAS$Value,
    omega = fixture$model$OMEGAS$Value
  )
  shared_modes <- .nm_subject_modes(
    shared, parameters, maxit = 100L, tolerance = 1e-9
  )
  ordinary_modes <- .nm_subject_modes(
    ordinary, parameters, maxit = 100L, tolerance = 1e-9
  )
  expect_true(shared$shared_fo_objective)
  expect_lt(
    sum(vapply(shared$subjects, function(x) x$tape_records, integer(1))),
    shared$n_subjects
  )
  expect_equal(
    vapply(shared_modes, `[[`, numeric(1), "value"),
    vapply(ordinary_modes, `[[`, numeric(1), "value"), tolerance = 1e-10
  )
  expect_equal(
    do.call(rbind, lapply(shared_modes, `[[`, "par")),
    do.call(rbind, lapply(ordinary_modes, `[[`, "par")), tolerance = 1e-9
  )
})

test_that("FO context cache rejects changed data and core timing excludes setup", {
  fixture <- estimation_fixture(FALSE)
  .nm_estimation_context_cache_clear()
  first <- .nm_estimation_context(fixture$model, fixture$data, method = "FO")
  second <- .nm_estimation_context(fixture$model, fixture$data, method = "FO")
  changed <- fixture$data
  changed$DV[[2L]] <- changed$DV[[2L]] + 0.001
  third <- .nm_estimation_context(fixture$model, changed, method = "FO")
  expect_false(first$cache_hit)
  expect_true(second$cache_hit)
  expect_false(third$cache_hit)

  fit <- nm_est(
    fixture$model, fixture$data, method = "FO", maxit = 5L,
    covariance = TRUE, covariance_type = "hessian"
  )
  expect_true(is.finite(fit$timing$initialization_seconds))
  expect_equal(
    fit$timing$total_seconds,
    fit$timing$model_fit_seconds + fit$timing$covariance_seconds,
    tolerance = 1e-12
  )
  expect_lte(fit$timing$total_seconds, fit$timing$wall_total_seconds + 1e-8)
})

test_that("optimized SAEM exposes optimizer policy without changing canonical Q", {
  fixture <- estimation_fixture(FALSE)
  model <- nm_model_update(
    fixture$model, NUMERICAL_MODE = "liber_optimized"
  )
  controls <- list(
    model = model, data = fixture$data, method = "SAEM", maxit = 5L,
    n_iter = 8L, burn = 3L, mcmc_steps = 1L, mstep_maxit = 3L,
    seed = 81L
  )
  automatic <- do.call(nm_est, c(controls, list(optimizer_backend = "auto")))
  native <- do.call(nm_est, c(controls, list(optimizer_backend = "native")))
  expect_match(
    automatic$diagnostics$optimizer$backend,
    "native-cpp-weighted-eta-lbfgs"
  )
  expect_match(
    native$diagnostics$optimizer$backend,
    "native-cpp-weighted-eta-lbfgs"
  )
  expect_match(
    native$diagnostics$estimator_identity,
    "accelerated f-SAEM auxiliary-function"
  )
  expect_true(all(is.finite(c(native$theta, native$sigma, native$omega))))
})

test_that("optimized SAEM uses native burn-in steps and canonical Q thereafter", {
  fixture <- estimation_fixture(FALSE)
  fit <- nm_est(
    fixture$model, fixture$data, method = "SAEM", maxit = 3L,
    n_iter = 6L, burn = 2L, mcmc_steps = 1L, mstep_maxit = 2L,
    seed = 915L, collect_output = FALSE,
    numerical_mode = "liber_optimized", optimizer_backend = "native"
  )
  expect_true(fit$diagnostics$native_mstep$eligible)
  expect_gte(fit$diagnostics$native_mstep$successes, 1L)
  expect_lt(fit$diagnostics$native_mstep$successes, 6L)
  expect_equal(fit$diagnostics$native_mstep$fallbacks, 0L)
  expect_identical(
    fit$diagnostics$parameter_averaging$resolved, "polyak"
  )
  expect_true(fit$diagnostics$parameter_averaging$applied)
  expect_lt(fit$diagnostics$mstep_schedule$count, 6L)
  expect_gt(
    fit$diagnostics$stochastic_approximation$retained_latent_support, 1L
  )
  expect_equal(sum(
    fit$diagnostics$stochastic_approximation$retained_weights
  ), 1, tolerance = 1e-14)
})

test_that("SAEM Q state follows the Robbins-Monro recurrence", {
  state <- LibeRation:::.nm_saem_q_state(list(n_subjects = 1L, n_eta = 1L))
  state$update(matrix(1), 1)
  state$update(matrix(2), 0.5)
  state$update(matrix(4), 0.25)
  expect_equal(state$weights(), c(0.375, 0.375, 0.25), tolerance = 1e-14)
  expect_equal(state$mean_eta(), matrix(2.125), tolerance = 1e-14)
})

test_that("native SAEM sufficient statistics match the compatibility path", {
  fixture <- estimation_fixture(FALSE)
  context <- .nm_estimation_context(
    fixture$model, fixture$data, method = "SAEM"
  )
  parameters <- list(
    theta = fixture$model$THETAS$Value,
    sigma = fixture$model$SIGMAS$Value,
    omega = fixture$model$OMEGAS$Value
  )
  eta <- matrix(
    rep(c(0.04, -0.03), length.out = context$n_subjects * context$n_eta),
    context$n_subjects, context$n_eta
  )
  old <- options(LibeRation.saem_native_sufficient_statistics = FALSE)
  on.exit(options(old), add = TRUE)
  omega_reference <- .nm_saem_omega_sufficient(context, eta)
  sigma_reference <- .nm_saem_sigma_sufficient(context, parameters, eta)
  options(LibeRation.saem_native_sufficient_statistics = TRUE)
  expect_equal(
    .nm_saem_omega_sufficient(context, eta), omega_reference,
    tolerance = 1e-12
  )
  expect_equal(
    .nm_saem_sigma_sufficient(context, parameters, eta), sigma_reference,
    tolerance = 1e-12
  )
})

test_that("guarded FO low-rank covariance matches dense in compatibility mode", {
  fixture <- estimation_fixture(FALSE)
  fixture$model$NUMERICAL_MODE <- "nonmem_compatibility"
  parameters <- list(
    theta = fixture$model$THETAS$Value,
    sigma = fixture$model$SIGMAS$Value,
    omega = fixture$model$OMEGAS$Value
  )
  old <- options(LibeRation.fo_low_rank = FALSE)
  on.exit(options(old), add = TRUE)
  dense_context <- .nm_estimation_context(
    fixture$model, fixture$data, method = "FO"
  )
  dense <- lapply(dense_context$subjects, function(subject) {
    subject$fo_objective(
      parameters$theta, parameters$sigma, parameters$omega,
      gradient = TRUE
    )
  })
  options(LibeRation.fo_low_rank = TRUE)
  low_rank_context <- .nm_estimation_context(
    fixture$model, fixture$data, method = "FO"
  )
  low_rank <- lapply(low_rank_context$subjects, function(subject) {
    subject$fo_objective(
      parameters$theta, parameters$sigma, parameters$omega,
      gradient = TRUE
    )
  })
  expect_equal(
    vapply(low_rank, `[[`, numeric(1), "value"),
    vapply(dense, `[[`, numeric(1), "value"), tolerance = 2e-11
  )
  expect_equal(
    do.call(rbind, lapply(low_rank, `[[`, "gradient")),
    do.call(rbind, lapply(dense, `[[`, "gradient")), tolerance = 2e-10
  )
})

test_that("guarded FO low-rank covariance covers IOV and excludes AR1", {
  iov_model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "DV", "MDV", "OCC"),
    ADVAN = 1,
    PRED = "CL=THETA(1)*exp(ETA(1)+ETA(2)); V=THETA(2); S1=V",
    ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(THETA = 1:2, Value = c(2, 20), FIX = TRUE),
    OMEGAS = data.frame(OMEGA = 1:2, Value = c(0.1, 0.05), FIX = TRUE),
    SIGMAS = data.frame(SIGMA = 1, Value = 0.2, FIX = TRUE), IOV = 1,
    NUMERICAL_MODE = "nonmem_compatibility"
  )
  subject <- data.frame(
    TIME = c(0, 1, 2, 24, 25, 26), EVID = c(1, 0, 0, 1, 0, 0),
    AMT = c(100, 0, 0, 100, 0, 0),
    DV = c(NA, 4.5, 4.0, NA, 4.4, 3.9),
    MDV = c(1, 0, 0, 1, 0, 0), OCC = c(1, 1, 1, 2, 2, 2)
  )
  iov_data <- rbind(transform(subject, ID = 1L), transform(subject, ID = 2L))
  iov_data <- iov_data[, c("ID", setdiff(names(iov_data), "ID"))]
  parameters <- list(
    theta = iov_model$THETAS$Value, sigma = iov_model$SIGMAS$Value,
    omega = iov_model$OMEGAS$Value
  )
  old <- options(LibeRation.fo_low_rank = FALSE)
  on.exit(options(old), add = TRUE)
  dense_context <- .nm_estimation_context(iov_model, iov_data, method = "FO")
  dense <- dense_context$subjects[[1L]]$fo_objective(
    parameters$theta, parameters$sigma, parameters$omega, gradient = TRUE
  )
  options(LibeRation.fo_low_rank = TRUE)
  low_rank_context <- .nm_estimation_context(iov_model, iov_data, method = "FO")
  low_rank <- low_rank_context$subjects[[1L]]$fo_objective(
    parameters$theta, parameters$sigma, parameters$omega, gradient = TRUE
  )
  expect_equal(low_rank$value, dense$value, tolerance = 2e-11)
  expect_equal(low_rank$gradient, dense$gradient, tolerance = 2e-10)

  ar1 <- estimation_fixture(FALSE)
  ar1$model$NUMERICAL_MODE <- "nonmem_compatibility"
  ar1$model$LIK_CONFIG$sigma_corr <- "ar1"
  ar1$model$LIK_CONFIG$ar1_rho <- 0.3
  ar1_parameters <- list(
    theta = ar1$model$THETAS$Value, sigma = ar1$model$SIGMAS$Value,
    omega = ar1$model$OMEGAS$Value
  )
  options(LibeRation.fo_low_rank = FALSE)
  ar1_dense_context <- .nm_estimation_context(
    ar1$model, ar1$data, method = "FO"
  )
  ar1_dense <- ar1_dense_context$subjects[[1L]]$fo_objective(
    ar1_parameters$theta, ar1_parameters$sigma, ar1_parameters$omega,
    gradient = TRUE
  )
  options(LibeRation.fo_low_rank = TRUE)
  ar1_guarded_context <- .nm_estimation_context(
    ar1$model, ar1$data, method = "FO"
  )
  ar1_guarded <- ar1_guarded_context$subjects[[1L]]$fo_objective(
    ar1_parameters$theta, ar1_parameters$sigma, ar1_parameters$omega,
    gradient = TRUE
  )
  expect_equal(ar1_guarded$value, ar1_dense$value, tolerance = 0)
  expect_equal(ar1_guarded$gradient, ar1_dense$gradient, tolerance = 0)
})

test_that("scalar FO population tape matches vector in compatibility mode", {
  fixture <- estimation_fixture(FALSE)
  fixture$model$NUMERICAL_MODE <- "nonmem_compatibility"
  fixture$model$THETAS$FIX <- c(FALSE, TRUE)
  fixture$model$SIGMAS$FIX <- FALSE
  fixture$model$OMEGAS$FIX <- FALSE
  context <- .nm_estimation_context(fixture$model, fixture$data, method = "FO")
  map <- .nm_outer_map(context$model)
  old <- options(
    LibeRation.fo_population_batch = TRUE,
    LibeRation.fo_population_scalar = FALSE
  )
  on.exit(options(old), add = TRUE)
  vector_route <- .nm_cpp_population_objective(
    context, map, "fo", eta_maxit = 80L, tolerance = 1e-8
  )
  vector_value <- .liberation_population_objective_value(
    vector_route$pointer, map$start
  )
  vector_gradient <- .liberation_population_objective_gradient(
    vector_route$pointer, map$start
  )
  vector_telemetry <- .liberation_population_objective_telemetry(
    vector_route$pointer
  )

  options(LibeRation.fo_population_scalar = TRUE)
  scalar_route <- .nm_cpp_population_objective(
    context, map, "fo", eta_maxit = 80L, tolerance = 1e-8
  )
  scalar_value <- .liberation_population_objective_value(
    scalar_route$pointer, map$start
  )
  scalar_gradient <- .liberation_population_objective_gradient(
    scalar_route$pointer, map$start
  )
  scalar_state <- .liberation_population_objective_state(
    scalar_route$pointer, map$start
  )
  scalar_telemetry <- .liberation_population_objective_telemetry(
    scalar_route$pointer
  )
  expect_equal(scalar_value, vector_value, tolerance = 2e-11)
  expect_equal(scalar_gradient, vector_gradient, tolerance = 2e-10)
  expect_false(vector_telemetry$fo_population_scalar)
  expect_true(scalar_telemetry$fo_population_scalar)
  expect_gt(scalar_telemetry$fo_low_rank_subjects, 0L)
  expect_equal(scalar_telemetry$fo_low_rank_fallbacks, 0L)
  expect_true(all(is.finite(vapply(
    scalar_state$modes, `[[`, numeric(1), "value"
  ))))
})

test_that("ill-conditioned FO covariance automatically uses the dense route", {
  fixture <- estimation_fixture(FALSE)
  fixture$model$NUMERICAL_MODE <- "nonmem_compatibility"
  fixture$model <- nm_model_update(
    fixture$model,
    PRED = paste(
      "CL=THETA(1)*exp(ETA(1));",
      "V=THETA(2)*exp(ETA(2)); S1=V"
    ),
    OMEGAS = data.frame(
      OMEGA = 1:2, Value = c(0.09, 2e-14), FIX = FALSE
    )
  )
  context <- .nm_estimation_context(fixture$model, fixture$data, method = "FO")
  map <- .nm_outer_map(context$model)
  old <- options(
    LibeRation.fo_population_batch = TRUE,
    LibeRation.fo_population_scalar = TRUE,
    LibeRation.fo_low_rank = TRUE
  )
  on.exit(options(old), add = TRUE)
  compiled <- .nm_cpp_population_objective(
    context, map, "fo", eta_maxit = 80L, tolerance = 1e-8
  )
  telemetry <- .liberation_population_objective_telemetry(compiled$pointer)
  expect_equal(telemetry$fo_low_rank_subjects, 0L)
  expect_gt(telemetry$fo_low_rank_fallbacks, 0L)
  expect_match(
    telemetry$fo_low_rank_fallback_reason, "ill-conditioned|singular"
  )
  expect_true(is.finite(.liberation_population_objective_value(
    compiled$pointer, map$start
  )))
})

test_that("compatible fused FO tape reuses invariant dynamic data", {
  fixture <- estimation_fixture(FALSE)
  context <- .nm_estimation_context(fixture$model, fixture$data, method = "FO")
  map <- .nm_outer_map(context$model)
  old <- options(LibeRation.fo_population_batch = TRUE)
  on.exit(options(old), add = TRUE)
  compiled <- .nm_cpp_population_objective(
    context, map, "fo", eta_maxit = 80L, tolerance = 1e-8
  )
  initial <- .liberation_population_objective_telemetry(compiled$pointer)
  invisible(.liberation_population_objective_value(compiled$pointer, map$start))
  invisible(.liberation_population_objective_gradient(compiled$pointer, map$start))
  moved <- map$start
  if (length(moved)) moved[[1L]] <- moved[[1L]] + 1e-4
  invisible(.liberation_population_objective_value(compiled$pointer, moved))
  final <- .liberation_population_objective_telemetry(compiled$pointer)
  expect_identical(final$fo_dynamic_updates, initial$fo_dynamic_updates)
})

test_that("paired SAEM callbacks retain R objective and gradient", {
  fixture <- estimation_fixture(FALSE)
  context <- .nm_estimation_context(fixture$model, fixture$data, method = "SAEM")
  mstep_model <- context$model
  mstep_model$OMEGAS$FIX[] <- TRUE
  map <- .nm_outer_map(mstep_model)
  parameters <- map$decode(map$start)
  eta <- matrix(
    rep(c(0.04, -0.03), length.out = context$n_subjects * context$n_eta),
    context$n_subjects, context$n_eta
  )
  paired <- .nm_saem_paired_conditional(context, map, eta)
  expect_equal(
    paired$objective(parameters),
    .nm_saem_conditional(context, parameters, eta), tolerance = 1e-12
  )
  expect_equal(
    paired$gradient(parameters),
    .nm_saem_conditional_gradient(context, map, parameters, eta),
    tolerance = 1e-11
  )
  telemetry <- paired$telemetry()
  expect_true(telemetry$paired)
  expect_true(telemetry$persistent)
  expect_equal(telemetry$evaluations, 1L)
  expect_gte(telemetry$cache_hits, 1L)

  optimized_model <- nm_model_update(
    fixture$model, NUMERICAL_MODE = "liber_optimized"
  )
  optimized_context <- .nm_estimation_context(
    optimized_model, fixture$data, method = "SAEM"
  )
  optimized_map_model <- optimized_context$model
  optimized_map_model$OMEGAS$FIX[] <- TRUE
  optimized_map <- .nm_outer_map(optimized_map_model)
  optimized <- .nm_saem_paired_conditional(
    optimized_context, optimized_map, eta
  )
  optimized_parameters <- optimized_map$decode(optimized_map$start)
  expect_equal(
    optimized$objective(optimized_parameters), paired$objective(parameters),
    tolerance = 1e-12
  )
  expect_equal(
    optimized$gradient(optimized_parameters), paired$gradient(parameters),
    tolerance = 1e-11
  )
  expect_true(optimized$telemetry()$native_aggregate)
})

test_that("compatible SAEM retains R optimizer with paired native tapes", {
  fixture <- estimation_fixture(FALSE)
  fit <- nm_est(
    fixture$model, fixture$data, method = "SAEM", maxit = 3L,
    n_iter = 4L, burn = 1L, mcmc_steps = 1L, mstep_maxit = 2L,
    seed = 916L, collect_output = FALSE,
    optimizer_backend = "r", numerical_mode = "nonmem_compatibility"
  )
  expect_gt(fit$diagnostics$paired_value_gradient$iterations, 0L)
  expect_gt(fit$diagnostics$paired_value_gradient$cache_hits, 0L)
  expect_gt(
    fit$diagnostics$paired_value_gradient$persistent_iterations, 0L
  )
  expect_equal(
    fit$diagnostics$paired_value_gradient$persistent_fallbacks, 0L
  )
  expect_identical(
    fit$diagnostics$optimizer$objective_backend,
    "r-optimizer+native-paired-value-gradient"
  )
})

test_that("optimized estimator kernels retain fixed-model objectives", {
  fixture <- estimation_fixture()
  controls <- list(
    FO = list(), FOCE = list(), FOCEI = list(), LAPLACE = list(),
    ITS = list(
      its_mstep_schedule = "fixed", its_acceleration = "none",
      its_eta_schedule = "fixed"
    ),
    GQ = list(gq_order = 3L),
    IMP = list(
      n_imp = 16L, seed = 912L, imp_sampling = "random",
      imp_proposal = "gaussian", imp_sample_schedule = "fixed",
      imp_mstep_schedule = "fixed", imp_mstep_maxit = 1L,
      imp_reuse_modes = FALSE,
      imp_auto_stop = FALSE, imp_subject_allocation = "fixed"
    ),
    SAEM = list(
      n_iter = 6L, burn = 2L, mcmc_steps = 1L, seed = 913L,
      saem_kernel = "random_walk", auto_stop = FALSE,
      saem_mstep_interval_burn = 1L, saem_mstep_interval = 1L,
      saem_parameter_averaging = "none"
    ),
    BAYES = list(
      n_burn = 2L, n_sample = 5L, n_thin = 1L, seed = 914L,
      outer_kernel = "isotropic", eta_kernel = "random_walk",
      delayed_rejection_scale = 0
    ),
    NPML = list(np_points = 5L, np_cycles = 1L, np_weight_maxit = 50L),
    NPAG = list(
      np_points = 5L, np_cycles = 1L, np_weight_maxit = 50L,
      np_max_candidates = 10L
    )
  )
  for (method in names(controls)) {
    common <- c(list(
      model = fixture$model, data = fixture$data, method = method,
      maxit = 2L, eta_maxit = 50L, tolerance = 1e-7,
      collect_output = FALSE
    ), controls[[method]])
    compatible <- do.call(
      nm_est, c(common, list(numerical_mode = "nonmem_compatibility"))
    )
    optimized <- do.call(
      nm_est, c(common, list(numerical_mode = "liber_optimized"))
    )
    expect_true(is.finite(optimized$objective), info = method)
    expect_equal(
      optimized$objective, compatible$objective,
      tolerance = if (method %in% c("BAYES", "SAEM")) 1e-8 else if (
        method %in% c("NPML", "NPAG")
      ) 1e-4 else 1e-10,
      info = method
    )
  }
})

test_that("optimized SAEM selects f-SAEM for multivariate random effects", {
  fixture <- estimation_fixture()
  model <- nm_model_update(
    fixture$model,
    PRED = paste(
      "CL=THETA(1)*exp(ETA(1));",
      "V=THETA(2)*exp(ETA(2)); S1=V"
    ),
    OMEGAS = data.frame(
      OMEGA = 1:2, Value = c(0.09, 0.04), FIX = TRUE
    )
  )
  fit <- nm_est(
    model, fixture$data, method = "SAEM", maxit = 2L,
    n_iter = 5L, burn = 2L, mcmc_steps = 1L, mstep_maxit = 1L,
    seed = 919L, collect_output = FALSE,
    numerical_mode = "liber_optimized", saem_kernel = "auto",
    fsaem_refresh = 3L
  )
  sampler <- fit$diagnostics$eta_sampler
  expect_identical(sampler$resolved_kernel, "fsaem")
  expect_gt(sampler$fsaem$refreshes, 0L)
  expect_equal(sampler$fsaem$refresh_failures, 0L)
  expect_equal(sampler$fsaem$fallback_iterations, 0L)
  expect_gt(sampler$persistent_context$calls, 0L)
  expect_true(is.finite(fit$objective))

  expect_error(
    nm_est(
      model, fixture$data, method = "SAEM", n_iter = 2L, burn = 1L,
      mcmc_steps = 1L, collect_output = FALSE,
      numerical_mode = "nonmem_compatibility", saem_kernel = "fsaem"
    ),
    "requires liber_optimized"
  )
})

test_that("optimized BAYES adapts a full population proposal covariance", {
  fixture <- estimation_fixture(FALSE)
  fit <- nm_est(
    fixture$model, fixture$data, method = "BAYES",
    n_burn = 20L, n_sample = 10L, n_thin = 1L, seed = 920L,
    collect_output = FALSE, numerical_mode = "liber_optimized",
    outer_kernel = "auto", adaptive_start = 5L, adaptive_interval = 5L
  )
  sampler <- fit$diagnostics$outer_sampler
  expect_identical(sampler$resolved_kernel, "adaptive_metropolis")
  expect_gt(sampler$covariance_updates, 0L)
  expect_equal(dim(sampler$covariance), c(4L, 4L))
  expect_true(all(is.finite(sampler$covariance)))
  expect_true(all(is.finite(fit$posterior$population$mean)))
  expect_true(all(is.finite(fit$posterior$population$ess)))
  expect_true(fit$diagnostics$eta_sampler$native_coordinator$used)
  expect_identical(
    fit$diagnostics$optimizer$backend,
    "persistent-native-cpp-bayes-coordinator"
  )

  withr::local_options(LibeRation.bayes_native_coordinator = FALSE)
  fallback <- nm_est(
    fixture$model, fixture$data, method = "BAYES",
    n_burn = 3L, n_sample = 4L, n_thin = 1L, seed = 920L,
    collect_output = FALSE, numerical_mode = "liber_optimized",
    outer_kernel = "adaptive_metropolis", adaptive_start = 2L,
    adaptive_interval = 2L
  )
  expect_false(fallback$diagnostics$eta_sampler$native_coordinator$used)
  expect_identical(
    fallback$diagnostics$optimizer$backend, "r-coordinated-bayes"
  )

  compatible <- nm_est(
    fixture$model, fixture$data, method = "BAYES",
    n_burn = 2L, n_sample = 3L, n_thin = 1L, seed = 921L,
    collect_output = FALSE, numerical_mode = "nonmem_compatibility",
    outer_kernel = "auto"
  )
  expect_identical(
    compatible$diagnostics$outer_sampler$resolved_kernel, "isotropic"
  )
  expect_error(
    nm_est(
      fixture$model, fixture$data, method = "BAYES",
      n_burn = 2L, n_sample = 3L, collect_output = FALSE,
      numerical_mode = "nonmem_compatibility",
      outer_kernel = "adaptive_metropolis"
    ),
    "only in liber_optimized"
  )
})

test_that("native compatibility BAYES preserves the established chain", {
  fixture <- estimation_fixture(FALSE)
  controls <- list(
    model = fixture$model, data = fixture$data, method = "BAYES",
    n_burn = 8L, n_sample = 12L, n_thin = 1L, seed = 929L,
    collect_output = FALSE, numerical_mode = "nonmem_compatibility",
    outer_kernel = "isotropic", eta_kernel = "random_walk",
    delayed_rejection_scale = 0
  )
  previous <- options(LibeRation.bayes_native_coordinator = FALSE)
  on.exit(options(previous), add = TRUE)
  reference <- do.call(nm_est, controls)
  options(LibeRation.bayes_native_coordinator = TRUE)
  native <- do.call(nm_est, controls)

  expect_false(
    reference$diagnostics$eta_sampler$native_coordinator$used
  )
  expect_true(native$diagnostics$eta_sampler$native_coordinator$used)
  expect_identical(
    native$diagnostics$eta_sampler$native_coordinator$policy,
    "compatibility-preserving-native"
  )
  expect_true(
    native$diagnostics$eta_sampler$native_coordinator$
      compatibility_preserving
  )
  expect_false(
    native$diagnostics$eta_sampler$native_coordinator$
      conjugate_omega$enabled
  )
  expect_equal(native$chain, reference$chain, tolerance = 2e-11)
  expect_equal(native$objective, reference$objective, tolerance = 2e-11)
  expect_equal(
    native$diagnostics$outer_acceptance,
    reference$diagnostics$outer_acceptance,
    tolerance = 1e-14
  )
  expect_equal(
    native$diagnostics$eta_acceptance,
    reference$diagnostics$eta_acceptance,
    tolerance = 1e-14
  )
})

test_that("optimized stochastic estimators expose robust Laplace proposals", {
  fixture <- estimation_fixture()
  model <- nm_model_update(
    fixture$model,
    PRED = paste(
      "CL=THETA(1)*exp(ETA(1));",
      "V=THETA(2)*exp(ETA(2)); S1=V"
    ),
    OMEGAS = data.frame(
      OMEGA = 1:2, Value = c(0.09, 0.04), FIX = TRUE
    )
  )
  saem <- nm_est(
    model, fixture$data, method = "SAEM", maxit = 2L,
    n_iter = 4L, burn = 2L, mcmc_steps = 1L, mstep_maxit = 1L,
    seed = 930L, collect_output = FALSE,
    numerical_mode = "liber_optimized", saem_kernel = "fsaem",
    fsaem_distribution = "student_t", fsaem_df = 6
  )
  expect_identical(
    saem$diagnostics$eta_sampler$fsaem$distribution, "student_t"
  )
  expect_equal(saem$diagnostics$eta_sampler$fsaem$degrees_of_freedom, 6)
  expect_true(is.finite(saem$objective))

  bayes <- nm_est(
    model, fixture$data, method = "BAYES", n_burn = 2L,
    n_sample = 3L, seed = 931L, collect_output = FALSE,
    numerical_mode = "liber_optimized", eta_kernel = "student_t",
    bayes_eta_df = 6
  )
  expect_identical(
    bayes$diagnostics$eta_sampler$resolved_kernel, "student_t"
  )
  expect_true(bayes$diagnostics$eta_sampler$native_coordinator$used)
  expect_true(all(is.finite(bayes$posterior$population$mean)))
})

test_that("optimized BAYES combines independent chains with diagnostics", {
  fixture <- estimation_fixture()
  fit <- nm_est(
    fixture$model, fixture$data, method = "BAYES",
    n_burn = 2L, n_sample = 4L, n_chains = 2L, seed = 932L,
    collect_output = FALSE, numerical_mode = "liber_optimized"
  )
  expect_length(fit$chains, 2L)
  expect_equal(nrow(fit$chain), 8L)
  expect_equal(fit$diagnostics$n_chains, 2L)
  expect_equal(fit$diagnostics$chain_seeds, c(932L, 933L))
  expect_named(
    fit$posterior$population$rhat,
    .nm_parameter_names(fit$theta, fit$sigma, fit$omega)
  )
  expect_true(all(is.finite(c(fit$theta, fit$sigma, fit$omega))))
})

test_that("native optimized BAYES supports expanded IOV ETAs", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "DV", "OCC"), ADVAN = 1,
    PRED = "CL=THETA(1)*exp(ETA(1)+ETA(2)); V=THETA(2); S1=V",
    ERROR = "Y=F+ERR(1)", THETAS = data.frame(
      THETA = 1:2, Value = c(2, 20), FIX = TRUE
    ), OMEGAS = data.frame(
      OMEGA = 1:2, Value = c(0.1, 0.05), FIX = TRUE
    ), SIGMAS = data.frame(SIGMA = 1, Value = 0.2, FIX = TRUE), IOV = 1
  )
  data <- data.frame(
    ID = c(1, 1, 1, 1, 2, 2), TIME = c(0, 1, 24, 25, 0, 1),
    EVID = c(1, 0, 1, 0, 1, 0), AMT = c(100, 0, 100, 0, 100, 0),
    DV = c(NA, 4.5, NA, 4.4, NA, 4.6), OCC = c(1, 1, 2, 2, 1, 1)
  )
  fit <- nm_est(
    model, data, method = "BAYES", n_burn = 1L, n_sample = 2L,
    seed = 934L, collect_output = FALSE,
    numerical_mode = "liber_optimized", eta_kernel = "student_t"
  )
  expect_equal(dim(fit$eta), c(2L, 3L))
  expect_true(fit$diagnostics$eta_sampler$native_coordinator$used)
  expect_true(is.finite(fit$objective))
})

test_that("shared conditional-state cache reuses results and reports telemetry", {
  context <- list(n_eta = 1L, n_subjects = 2L)
  calls <- 0L
  cache <- .nm_conditional_state_cache(
    context,
    function(parameters, starts) {
      calls <<- calls + 1L
      list(
        value = sum(parameters$theta),
        states = lapply(c(0.1, -0.1), function(value) {
          list(mode = list(par = value))
        })
      )
    }
  )
  parameters <- list(theta = 1, sigma = 0.2, omega = 0.1)
  expect_equal(cache$evaluate(parameters)$value, 1)
  expect_equal(cache$evaluate(parameters)$value, 1)
  expect_equal(calls, 1L)
  expect_equal(cache$telemetry()$hits, 1L)
  expect_equal(cache$telemetry()$misses, 1L)
})

test_that("optimized analytical ADVAN stochastic execution is thread stable", {
  fixture <- estimation_fixture(FALSE)
  withr::local_options(list(
    LibeRation.fused_advan_stochastic = TRUE,
    LibeRation.native_subject_threads = TRUE
  ))
  controls <- list(
    model = fixture$model, data = fixture$data, method = "BAYES",
    n_burn = 2L, n_sample = 4L, n_thin = 1L, seed = 12041L,
    collect_output = FALSE, numerical_mode = "liber_optimized",
    outer_kernel = "isotropic", eta_kernel = "random_walk",
    delayed_rejection_scale = 0
  )
  serial <- do.call(nm_est, c(controls, list(n_cores = 1L)))
  threaded <- do.call(nm_est, c(controls, list(n_cores = 2L)))
  telemetry <- threaded$diagnostics$eta_sampler$persistent_context
  expect_true(telemetry$fused_advan_enabled)
  expect_equal(telemetry$native_subject_threads, 2L)
  expect_true(telemetry$persistent_worker_pool)
  expect_gt(telemetry$worker_pool_dispatches, 0L)
  expect_gt(telemetry$fused_advan_evaluations, 0L)
  expect_equal(threaded$chain, serial$chain, tolerance = 2e-11)
  expect_equal(threaded$objective, serial$objective, tolerance = 2e-11)
})

test_that("MCMC diagnostics expose rank-normalized bulk and tail measures", {
  set.seed(12042L)
  chains <- lapply(seq_len(4L), function(index) {
    matrix(stats::rnorm(800L), ncol = 2L)
  })
  diagnostics <- .nm_mcmc_diagnostics(chains)
  expect_named(
    diagnostics, c("rhat", "bulk_ess", "tail_ess", "mcse_mean")
  )
  expect_true(all(is.finite(diagnostics$rhat)))
  expect_true(all(diagnostics$rhat > 0.9 & diagnostics$rhat < 1.1))
  expect_true(all(is.finite(diagnostics$bulk_ess)))
  expect_true(all(is.finite(diagnostics$tail_ess)))
  expect_true(all(is.finite(diagnostics$mcse_mean)))
})

test_that("SAEM reports stationarity and independent replicate evidence", {
  fixture <- estimation_fixture(FALSE)
  fit <- nm_est(
    fixture$model, fixture$data, method = "SAEM", maxit = 2L,
    n_iter = 5L, burn = 2L, mcmc_steps = 1L, mstep_maxit = 1L,
    n_replicates = 2L, replicate_seed_stride = 19L, seed = 12043L,
    replicate_score_samples = 20L,
    stationarity_window = 4L, collect_output = FALSE,
    numerical_mode = "liber_optimized"
  )
  expect_equal(fit$diagnostics$replicates$count, 2L)
  expect_equal(fit$diagnostics$replicates$seeds, c(12043L, 12062L))
  expect_length(fit$diagnostics$replicates$objectives, 2L)
  expect_length(fit$diagnostics$replicates$selection_scores, 2L)
  expect_true(all(is.finite(fit$diagnostics$replicates$selection_scores)))
  expect_match(
    fit$diagnostics$replicates$selection_metric,
    "importance marginal objective"
  )
  expect_identical(
    fit$objective,
    fit$diagnostics$replicates$selection_scores[
      fit$diagnostics$replicates$selected
    ]
  )
  expect_true(fit$objective_comparable)
  expect_match(fit$objective_type, "marginal_log_likelihood")
  expect_named(
    fit$diagnostics$stationarity,
    c(
      "ready", "converged", "window", "parameter_drift",
      "objective_drift", "tolerance", "auto_stop", "stopped_early",
      "consecutive_confirmations", "required_confirmations",
      "minimum_iterations"
    )
  )
  expect_true(is.matrix(fit$diagnostics$parameter_trace))
  expect_true(is.finite(fit$objective))
})

test_that("native optimized BAYES owns and retapes ODE objectives", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "DV", "MDV"),
    ADVAN = 6, DOSECMP = 1, OBSCMP = 1,
    PRED = "K=THETA(1)*exp(ETA(1)); V=THETA(2); S1=V",
    DES = "DADT(1)=-K*A(1)", ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(THETA = 1:2, Value = c(0.1, 20), FIX = TRUE),
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
    model, data, method = "BAYES", n_burn = 1L, n_sample = 2L,
    seed = 935L, collect_output = FALSE,
    numerical_mode = "liber_optimized", eta_kernel = "student_t",
    bayes_eta_refresh = 1L
  )
  native <- fit$diagnostics$eta_sampler$native_coordinator
  persistent <- fit$diagnostics$eta_sampler$persistent_context
  expect_true(native$used)
  expect_true(persistent$ode_owned_tapes)
  expect_gte(persistent$tape_records, 0L)
  expect_true(is.finite(fit$objective))
})

test_that("optimized BAYES uses conjugate diagonal OMEGA updates when valid", {
  fixture <- estimation_fixture()
  fixture$model$OMEGAS$FIX[] <- FALSE
  fixture$model$LIK_CONFIG$priors <- nm_prior(
    "OMEGA1", "inverse_gamma", shape = 3, rate = 0.2
  )
  fit <- nm_est(
    fixture$model, fixture$data, method = "BAYES",
    n_burn = 2L, n_sample = 3L, seed = 936L,
    collect_output = FALSE, numerical_mode = "liber_optimized",
    bayes_gibbs_omega = TRUE
  )
  conjugate <- fit$diagnostics$eta_sampler$native_coordinator$conjugate_omega
  expect_equal(conjugate$parameters, 1L)
  expect_equal(conjugate$updates, 5L)
  expect_gte(conjugate$draws, conjugate$updates)
  expect_true(is.finite(fit$omega[[1L]]))
  expect_gt(fit$omega[[1L]], 0)
})

test_that("native optimized BAYES supports general random-effect designs", {
  data <- expand.grid(ID = 1:4, TIME = 1:2)
  data <- data[order(data$ID, data$TIME), ]
  data$SITE <- rep(c("A", "B"), each = 4L)
  data$READER <- rep(c("R1", "R2", "R2", "R1"), each = 2L)
  data$EVID <- 0L
  data$AMT <- 0
  data$CMT <- 1L
  data$MDV <- 0L
  data$DV <- 10 + data$ID / 10
  model <- nm_model(
    INPUT = names(data), ADVAN = 1,
    PRED = paste0(
      "CL=THETA(1)*exp(ETA(1)+ETA(2)+ETA(3));",
      "V=1;S1=1;F=CL"
    ),
    ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(THETA = 1L, Value = 10, FIX = TRUE),
    OMEGAS = data.frame(
      OMEGA = 1:3, Value = c(0.1, 0.2, 0.15), FIX = TRUE
    ),
    SIGMAS = data.frame(SIGMA = 1L, Value = 0.1, FIX = TRUE),
    RE_CONFIG = nm_re_config(
      nm_re_block("site", "SITE", 1L),
      nm_re_block("patient", "ID", 2L),
      nm_re_block("reader", "READER", 3L)
    )
  )
  fit <- nm_est(
    model, data, method = "BAYES", n_burn = 1L, n_sample = 2L,
    seed = 937L, collect_output = FALSE,
    numerical_mode = "liber_optimized", eta_kernel = "student_t"
  )
  expect_equal(dim(fit$eta), c(1L, 8L))
  expect_true(fit$diagnostics$eta_sampler$native_coordinator$used)
  expect_true(is.finite(fit$objective))
})

test_that("optimized BAYES supports native subject-parallel sampling", {
  skip_on_cran()
  fixture <- estimation_fixture()
  fit <- nm_est(
    fixture$model, fixture$data, method = "BAYES",
    n_burn = 1L, n_sample = 2L, n_cores = 2L, seed = 938L,
    collect_output = FALSE, numerical_mode = "liber_optimized",
    eta_kernel = "student_t", bayes_eta_refresh = 1L
  )
  expect_true(fit$diagnostics$eta_sampler$native_coordinator$used)
  expect_identical(
    fit$diagnostics$optimizer$backend,
    "persistent-native-cpp-bayes-coordinator"
  )
  expect_true(all(is.finite(c(fit$theta, fit$sigma, fit$omega))))
  expect_true(is.finite(fit$objective))
})
