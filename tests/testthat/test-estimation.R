test_that("FO uses the first-order marginal covariance", {
  fixture <- estimation_fixture()
  fit <- nm_est(fixture$model, fixture$data, method = "FO", maxit = 5)
  expect_s3_class(fit, "nm_fit")
  expect_true(is.finite(fit$objective))
  expect_equal(dim(fit$eta), c(3, 1))
  expect_true(is.numeric(fit$iterations) && fit$iterations >= 0L)
  expect_true(is.numeric(fit$objective_evaluations) && fit$objective_evaluations >= 1L)
  expect_true(is.finite(fit$timing$model_fit_seconds))
  expect_true(is.finite(fit$timing$total_seconds))
  expect_gte(fit$timing$total_seconds, fit$timing$model_fit_seconds)
})

test_that("FO marginal values and population gradients are taped exactly", {
  fixture <- estimation_fixture()
  fixture$model$THETAS$FIX <- c(FALSE, TRUE)
  fixture$model$SIGMAS$FIX <- FALSE
  fixture$model$OMEGAS$FIX <- FALSE
  context <- LibeRation:::.nm_estimation_context(fixture$model, fixture$data)
  map <- LibeRation:::.nm_outer_map(context$model)
  parameters <- map$decode(map$start)
  native <- vapply(
    context$subjects, LibeRation:::.nm_fo_subject, numeric(1),
    model = context$model, theta = parameters$theta,
    sigma = parameters$sigma, omega = parameters$omega
  )
  reference <- vapply(
    context$subjects, LibeRation:::.nm_fo_subject_reference, numeric(1),
    model = context$model, theta = parameters$theta,
    sigma = parameters$sigma, omega = parameters$omega
  )
  expect_equal(native, reference, tolerance = 2e-11)

  objective <- function(outer) {
    LibeRation:::.nm_fo_objective(context, map$decode(outer))
}

  exact <- LibeRation:::.nm_fo_outer_gradient(context, map, parameters)
  numerical <- vapply(seq_along(map$start), function(index) {
    step <- 1e-5 * max(1, abs(map$start[[index]]))
    high <- low <- map$start
    high[[index]] <- high[[index]] + step
    low[[index]] <- low[[index]] - step
    (objective(high) - objective(low)) / (2 * step)
  }, numeric(1))
  expect_equal(exact, numerical, tolerance = 3e-6)
})

test_that("FO shares structural likelihood tapes and fuses the population safely", {
  fixture <- estimation_fixture(FALSE)
  fixture$model$THETAS$FIX <- c(FALSE, TRUE)
  fixture$model$SIGMAS$FIX <- FALSE
  fixture$model$OMEGAS$FIX <- FALSE
  context <- LibeRation:::.nm_estimation_context(
    fixture$model, fixture$data, method = "FO"
  )
  expect_true(all(vapply(
    context$subjects, function(subject) is.null(subject$noninteraction_tape),
    logical(1)
  )))
  map <- LibeRation:::.nm_outer_map(context$model)

  old <- options(LibeRation.fo_population_batch = FALSE)
  on.exit(options(old), add = TRUE)
  subject_route <- LibeRation:::.nm_cpp_population_objective(
    context, map, "fo", eta_maxit = 80L, tolerance = 1e-8
  )
  subject_value <- LibeRation:::.liberation_population_objective_value(
    subject_route$pointer, map$start
  )
  subject_gradient <- LibeRation:::.liberation_population_objective_gradient(
    subject_route$pointer, map$start
  )
  subject_telemetry <- LibeRation:::.liberation_population_objective_telemetry(
    subject_route$pointer
  )

  options(LibeRation.fo_population_batch = TRUE)
  fused_route <- LibeRation:::.nm_cpp_population_objective(
    context, map, "fo", eta_maxit = 80L, tolerance = 1e-8
  )
  fused_value <- LibeRation:::.liberation_population_objective_value(
    fused_route$pointer, map$start
  )
  fused_gradient <- LibeRation:::.liberation_population_objective_gradient(
    fused_route$pointer, map$start
  )
  fused_telemetry <- LibeRation:::.liberation_population_objective_telemetry(
    fused_route$pointer
  )

  expect_equal(fused_value, subject_value, tolerance = 2e-12)
  expect_equal(fused_gradient, subject_gradient, tolerance = 2e-11)
  expect_equal(subject_telemetry$fo_unique_subject_tapes, 1L)
  expect_equal(subject_telemetry$fo_shared_subject_tapes, 2L)
  expect_true(fused_telemetry$fo_population_batched)
  expect_gt(fused_telemetry$fo_population_operations, 0)

  options(LibeRation.fo_population_max_operations = 1)
  limited_route <- LibeRation:::.nm_cpp_population_objective(
    context, map, "fo", eta_maxit = 80L, tolerance = 1e-8
  )
  limited_telemetry <- LibeRation:::.liberation_population_objective_telemetry(
    limited_route$pointer
  )
  expect_false(limited_telemetry$fo_population_batched)
  expect_match(limited_telemetry$fo_population_error, "exceeds")
  options(LibeRation.fo_population_max_operations = NULL)

  # The R-coordinated fallback must also select each subject's dynamic data
  # correctly after the tape pool has been established.
  parameters <- map$decode(map$start)
  fallback_gradient <- LibeRation:::.nm_fo_outer_gradient(
    context, map, parameters
  )
  expect_equal(fallback_gradient, subject_gradient, tolerance = 2e-11)
})

test_that("FO tapes cover residual-error families and AR1 covariance", {
  for (variant in c("proportional", "exponential", "combined", "power", "ar1")) {
    fixture <- estimation_fixture()
    if (variant %in% c("combined", "power")) {
      fixture$model$SIGMAS <- data.frame(
        SIGMA = 1:2,
        Value = if (variant == "power") c(0.2, 0.7) else c(0.1, 0.2),
        FIX = FALSE
      )
    }
    if (variant %in% c("proportional", "exponential", "combined", "power")) {
      fixture$model$LIK_CONFIG$error <- variant
    }
    if (variant == "ar1") {
      fixture$model$LIK_CONFIG$sigma_corr <- "ar1"
      fixture$model$LIK_CONFIG$ar1_rho <- 0.3
    }
    context <- LibeRation:::.nm_estimation_context(fixture$model, fixture$data)
    parameters <- list(
      theta = context$model$THETAS$Value,
      sigma = context$model$SIGMAS$Value,
      omega = context$model$OMEGAS$Value
    )
    native <- LibeRation:::.nm_fo_subject(
      context$subjects[[1L]], context$model, parameters$theta,
      parameters$sigma, parameters$omega
    )
    reference <- LibeRation:::.nm_fo_subject_reference(
      context$subjects[[1L]], context$model, parameters$theta,
      parameters$sigma, parameters$omega
    )
    expect_equal(native, reference, tolerance = 2e-11, info = variant)
  }
})

test_that("conditional methods use exact ETA modes and curvature", {
  fixture <- estimation_fixture()
  for (method in c("FOCE", "FOCEI", "LAPLACE")) {
    fit <- nm_est(
      fixture$model, fixture$data, method = method,
      maxit = 5, eta_maxit = 80, tolerance = 1e-7
    )
    expect_s3_class(fit, "nm_fit")
    expect_true(is.finite(fit$objective), info = method)
    expect_true(all(fit$diagnostics$eta_convergence == 0), info = method)
    expect_lt(max(abs(fit$eta)), 1)
  }
})

test_that("Laplace curvature tapes are anchored at conditional modes", {
  model <- LibeRation:::.liber_model_template(2L, trans = 2L)
  data <- LibeRation:::.liber_builtin_dataset(
    model, "theophylline", n_subjects = 3L, seed = 1L
  )
  context <- LibeRation:::.nm_estimation_context(model, data)
  map <- LibeRation:::.nm_outer_map(context$model)

  expect_error(
    context$subjects[[2L]]$ensure_curvature_tape(
      model$THETAS$Value, rep(0, context$n_eta), model$SIGMAS$Value,
      model$OMEGAS$Value, "laplace"
    ),
    "not positive definite"
  )
  compiled <- LibeRation:::.nm_cpp_population_objective(
    context, map, "laplace", eta_maxit = 100L, tolerance = 1e-7
  )
  expect_false(is.null(compiled$pointer), info = compiled$reason)
  fit <- nm_est(
    model, data, method = "LAPLACE", maxit = 2L,
    eta_maxit = 100L, tolerance = 1e-7
  )
  expect_true(is.finite(fit$objective))
  expect_true(all(fit$diagnostics$eta_convergence == 0L))
})

test_that("FOCE freezes residual variance at ETA zero while FOCEI interacts", {
  fixture <- estimation_fixture()
  fixture$model$ERROR <- "Y=F+F*ERR(1)"
  fixture$model$ERROR_TYPE <- "proportional"
  fixture$model$LIK_CONFIG$error <- "proportional"
  foce <- nm_est(
    fixture$model, fixture$data, method = "FOCE",
    maxit = 2, eta_maxit = 80, tolerance = 1e-8
  )
  focei <- nm_est(
    fixture$model, fixture$data, method = "FOCEI",
    maxit = 2, eta_maxit = 80, tolerance = 1e-8
  )
  expect_equal(foce$diagnostics$approximation, "foce")
  expect_equal(focei$diagnostics$approximation, "focei")
  expect_gt(abs(foce$objective - focei$objective), 1e-4)
  expect_gt(max(abs(foce$eta - focei$eta)), 1e-5)
})

test_that("IOV estimation uses one stable expanded ETA layout for every subject", {
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
  fit <- nm_est(model, data, method = "FOCEI", maxit = 2, eta_maxit = 40)
  expect_equal(dim(fit$eta), c(2, 3))
  expect_true(is.finite(fit$objective))
})

test_that("ITS and IMP use the exact conditional objective", {
  fixture <- estimation_fixture()
  its <- nm_est(
    fixture$model, fixture$data, method = "ITS",
    maxit = 5, eta_maxit = 60
  )
  imp <- nm_est(
    fixture$model, fixture$data, method = "IMP",
    maxit = 5, eta_maxit = 60, n_imp = 20, seed = 42
  )
  expect_true(is.finite(its$objective))
  expect_true(is.finite(imp$objective))
  expect_equal(imp$diagnostics$n_imp, 20)
  expect_equal(imp$diagnostics$imp_gradient, "score")
})

test_that("GQ integrates the exact joint objective deterministically", {
  fixture <- estimation_fixture()
  context <- LibeRation:::.nm_estimation_context(fixture$model, fixture$data)
  parameters <- list(
    theta = fixture$model$THETAS$Value,
    sigma = fixture$model$SIGMAS$Value,
    omega = fixture$model$OMEGAS$Value
  )
  evaluator <- context$subjects[[1L]]

  design <- LibeRation:::.nm_gq_design(context, order = 15L)
  proposal <- LibeRation:::.nm_imp_subject_proposal(
    evaluator, parameters, design$normals[[1L]], 100L, 1e-9
  )
  quadrature <- LibeRation:::.nm_imp_subject_from_proposal(
    evaluator, parameters, proposal, gradient = TRUE
  )
  centered_integral <- stats::integrate(
    function(eta) vapply(eta, function(value) exp(-0.5 * (
      evaluator$objective(
        parameters$theta, value, parameters$sigma, parameters$omega,
        gradient = FALSE
      )$value - proposal$mode$value
    )), numeric(1)),
    lower = -8, upper = 8, rel.tol = 1e-10
  )$value
  reference <- proposal$mode$value - 2 * log(centered_integral)
  expect_equal(quadrature$value, reference, tolerance = 2e-6)
  expect_true(all(is.finite(quadrature$native_gradient)))

  one_node <- LibeRation:::.nm_gq_design(context, order = 1L)
  laplace_proposal <- LibeRation:::.nm_imp_subject_proposal(
    evaluator, parameters, one_node$normals[[1L]], 100L, 1e-9
  )
  one_node_value <- LibeRation:::.nm_imp_subject_from_proposal(
    evaluator, parameters, laplace_proposal, gradient = FALSE
  )$value
  expect_equal(
    one_node_value,
    laplace_proposal$mode$value + laplace_proposal$mode$logdet - log(4 * pi),
    tolerance = 1e-10
  )
})

test_that("GQ exposes adaptive and fixed integration with covariance support", {
  fixture <- estimation_fixture()
  adaptive <- nm_est(
    fixture$model, fixture$data, method = "GQ", maxit = 2L,
    eta_maxit = 80L, gq_order = 3L
  )
  fixed <- nm_est(
    fixture$model, fixture$data, method = "GQ", maxit = 2L,
    eta_maxit = 80L, gq_order = 3L, gq_adaptive = FALSE
  )
  expect_true(is.finite(adaptive$objective))
  expect_true(is.finite(fixed$objective))
  expect_equal(adaptive$diagnostics$quadrature_points, 3L)
  expect_true(adaptive$diagnostics$adaptive)
  expect_false(fixed$diagnostics$adaptive)
  covariance <- nm_cov_step(adaptive)
  expect_identical(covariance$status, "completed")
  expect_equal(dim(covariance$covariance), c(0L, 0L))

  free <- estimation_fixture(FALSE)
  optimized <- nm_est(
    free$model, free$data, method = "GQ", maxit = 4L,
    eta_maxit = 80L, gq_order = 5L, covariance = TRUE,
    covariance_type = "hessian"
  )
  expect_true(is.finite(optimized$objective))
  expect_true(all(is.finite(c(
    optimized$theta, optimized$omega, optimized$sigma
  ))))
  expect_identical(optimized$covariance$status, "completed")
  expect_true(all(is.finite(optimized$covariance$se)))
})

test_that("GQ selects and evaluates Smolyak sparse grids", {
  low <- LibeRation:::.nm_gq_design(
    list(n_eta = 3L, n_subjects = 2L), order = 5L,
    grid = "auto", level = 3L
  )
  high <- LibeRation:::.nm_gq_design(
    list(n_eta = 4L, n_subjects = 2L), order = 5L,
    grid = "auto", level = 3L
  )
  expect_identical(low$resolved_grid, "tensor")
  expect_identical(high$resolved_grid, "smolyak")
  expect_equal(high$actual_samples, 49L)
  expect_lt(high$actual_samples, 5^4)
  expect_equal(
    sum(attr(high$normals[[1L]], "measure_sign") *
      exp(attr(high$normals[[1L]], "log_measure"))),
    1, tolerance = 1e-12
  )

  fixture <- estimation_fixture()
  tensor <- nm_est(
    fixture$model, fixture$data, method = "GQ", maxit = 2L,
    eta_maxit = 80L, gq_grid = "tensor", gq_order = 5L
  )
  sparse <- nm_est(
    fixture$model, fixture$data, method = "GQ", maxit = 2L,
    eta_maxit = 80L, gq_grid = "smolyak", gq_level = 3L
  )
  expect_equal(sparse$objective, tensor$objective, tolerance = 1e-10)
  expect_identical(sparse$diagnostics$quadrature_grid, "smolyak")
  expect_equal(sparse$diagnostics$quadrature_level, 3L)
  expect_true(all(sparse$diagnostics$quadrature_cancellation_ratio > 0))

  free <- estimation_fixture(FALSE)
  free$model$THETAS$FIX <- c(FALSE, TRUE)
  free$model$OMEGAS$FIX <- TRUE
  free$model$SIGMAS$FIX <- TRUE
  sparse_covariance <- nm_est(
    free$model, free$data, method = "GQ", maxit = 4L,
    eta_maxit = 80L, gq_grid = "smolyak", gq_level = 3L,
    covariance = TRUE, covariance_type = "hessian"
  )
  expect_identical(sparse_covariance$covariance$status, "completed")
  expect_identical(sparse_covariance$covariance$quadrature_grid, "smolyak")
  expect_equal(sparse_covariance$covariance$quadrature_level, 3L)
})

test_that("four-ETA adaptive Smolyak GQ has a finite signed-grid likelihood", {
  fixture <- estimation_fixture()
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "DV", "MDV"),
    ADVAN = 1, DOSECMP = 1, OBSCMP = 1,
    PRED = paste(
      "CL=THETA(1)*exp(ETA(1)+0.2*ETA(3))",
      "V=THETA(2)*exp(ETA(2)+0.2*ETA(4))", "S1=V", sep = ";"
    ),
    ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(THETA = 1:2, Value = c(2, 20), FIX = TRUE),
    OMEGAS = data.frame(OMEGA = 1:4, Value = rep(0.09, 4L), FIX = TRUE),
    SIGMAS = data.frame(SIGMA = 1, Value = 0.2, FIX = TRUE)
  )
  fit <- nm_est(
    model, fixture$data, method = "GQ", maxit = 2L,
    eta_maxit = 80L, gq_grid = "smolyak", gq_level = 3L
  )
  expect_true(is.finite(fit$objective))
  expect_equal(fit$diagnostics$quadrature_points, 49L)
  expect_gt(fit$diagnostics$quadrature_negative_weights, 0L)
  expect_true(all(fit$diagnostics$quadrature_cancellation_ratio > 0))
})

test_that("SAEM and BAYES return reproducible stochastic diagnostics", {
  fixture <- estimation_fixture()
  saem <- nm_est(
    fixture$model, fixture$data, method = "SAEM",
    n_iter = 8, burn = 2, mcmc_steps = 1, mstep_maxit = 2, seed = 99
  )
  bayes <- nm_est(
    fixture$model, fixture$data, method = "BAYES",
    n_burn = 5, n_sample = 10, n_thin = 1, seed = 99
  )
  expect_true(is.finite(saem$objective))
  expect_true(saem$diagnostics$acceptance >= 0 && saem$diagnostics$acceptance <= 1)
  expect_length(saem$diagnostics$acceptance_trace, 8L)
  expect_true(all(saem$diagnostics$step_scale_trace > 0))
  expect_equal(nrow(bayes$chain), 10)
  expect_true(all(is.finite(bayes$posterior$mean)))
  expect_true(bayes$diagnostics$eta_acceptance >= 0 &&
                bayes$diagnostics$eta_acceptance <= 1)
})

test_that("SAEM scales steep M-steps before L-BFGS-B boundary searches", {
  model <- LibeRation:::.liber_model_template(2L, trans = 2L)
  data <- LibeRation:::.liber_builtin_dataset(
    model, "theophylline", n_subjects = 3L, seed = 1L
  )
  fit <- nm_est(
    model, data, method = "SAEM", n_iter = 3L, burn = 1L,
    mcmc_steps = 1L, mstep_maxit = 3L, seed = 20260713L
  )
  expect_true(is.finite(fit$objective))
  expect_true(all(is.finite(fit$theta)))
  expect_true(all(is.finite(fit$omega)))
  expect_true(all(is.finite(fit$sigma)))
})

test_that("native population optimizer and structural tape pool report telemetry", {
  fixture <- estimation_fixture(FALSE)
  fit <- nm_est(
    fixture$model, fixture$data, method = "FOCEI", maxit = 5,
    eta_maxit = 50, optimizer_backend = "native"
  )
  expect_equal(fit$diagnostics$optimizer$backend, "native-bfgs")
  expect_true(is.data.frame(fit$diagnostics$optimizer$trace))
  expect_gte(fit$diagnostics$optimizer$objective_evaluations, 1L)
  expect_equal(fit$diagnostics$tapes$unique_structures, 1L)
  expect_equal(fit$diagnostics$tapes$shared_prediction_tapes, 2L)
  expect_true(fit$diagnostics$conditional_modes$evaluations > 0L)
})

test_that("L-BFGS-B uses the persistent C++ population objective and same-point cache", {
  fixture <- estimation_fixture(FALSE)
  context <- LibeRation:::.nm_estimation_context(fixture$model, fixture$data)
  map <- LibeRation:::.nm_outer_map(context$model)
  compiled <- LibeRation:::.nm_cpp_population_objective(
    context, map, "focei", eta_maxit = 80L, tolerance = 1e-8
  )
  expect_false(is.null(compiled$pointer), info = compiled$reason)

  value <- LibeRation:::.liberation_population_objective_value(
    compiled$pointer, map$start
  )
  gradient <- LibeRation:::.liberation_population_objective_gradient(
    compiled$pointer, map$start
  )
  repeated <- LibeRation:::.liberation_population_objective_gradient(
    compiled$pointer, map$start
  )
  reference <- LibeRation:::.nm_nested_objective(
    context, "focei", eta_maxit = 80L, tolerance = 1e-8
  )
  parameters <- map$decode(map$start)
  expect_equal(value, reference(parameters), tolerance = 2e-8)
  expect_equal(
    gradient,
    LibeRation:::.nm_nested_outer_gradient(
      context, map, reference, parameters, "focei"
    ),
    tolerance = 2e-7
  )
  expect_equal(repeated, gradient)
  telemetry <- LibeRation:::.liberation_population_objective_telemetry(
    compiled$pointer
  )
  expect_equal(telemetry$parameter_evaluations, 1L)
  expect_gte(telemetry$shared_state_hits, 1L)
  expect_gte(telemetry$gradient_cache_hits, 1L)
  expect_identical(telemetry$propagation_kernel, "specialized-advan1")

  fit <- nm_est(
    fixture$model, fixture$data, method = "FOCEI", maxit = 5,
    eta_maxit = 80L, tolerance = 1e-7
  )
  expect_equal(fit$diagnostics$optimizer$backend, "r-l-bfgs-b-cpp-objective")
  expect_equal(
    fit$diagnostics$optimizer$objective_backend,
    "persistent-cpp-population-objective"
  )
  expect_gte(
    fit$diagnostics$optimizer$population_objective$shared_state_hits, 1L
  )
})

test_that("Laplace retaping rejects extreme conditional-mode trial anchors", {
  model <- LibeRation:::.liber_model_template(2L, trans = 2L)
  data <- LibeRation:::.liber_builtin_dataset(
    model, "theophylline", n_subjects = 3L, seed = 1L
  )
  context <- LibeRation:::.nm_estimation_context(model, data)
  map <- LibeRation:::.nm_outer_map(context$model)
  compiled <- LibeRation:::.nm_cpp_population_objective(
    context, map, "laplace", eta_maxit = 100L, tolerance = 1e-7
  )

  value <- LibeRation:::.liberation_population_objective_value(
    compiled$pointer, map$start
  )
  gradient <- LibeRation:::.liberation_population_objective_gradient(
    compiled$pointer, map$start
  )
  telemetry <- LibeRation:::.liberation_population_objective_telemetry(
    compiled$pointer
  )

  expect_true(is.finite(value))
  expect_true(all(is.finite(gradient)))
  expect_gt(telemetry$tape_retapes, 0L)
  expect_gte(telemetry$tape_records, telemetry$tape_retapes)
  expect_lt(telemetry$tape_retapes, 12L * context$n_subjects)
})

test_that("adaptive ODE objective tapes retape after material domain movement", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "DV", "MDV"),
    ADVAN = 6, DOSECMP = 1, OBSCMP = 1,
    PRED = "K=THETA(1)*exp(ETA(1)); V=THETA(2); S1=V",
    DES = "DADT(1)=-K*A(1)", ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(THETA = 1:2, Value = c(0.1, 20), FIX = TRUE),
    OMEGAS = data.frame(OMEGA = 1, Value = 0.1, FIX = TRUE),
    SIGMAS = data.frame(SIGMA = 1, Value = 0.2, FIX = TRUE)
  )
  data <- data.frame(
    ID = 1, TIME = c(0, 1, 4), EVID = c(1, 0, 0),
    AMT = c(100, 0, 0), DV = c(NA, 4.5, 3.3), MDV = c(1, 0, 0)
  )
  context <- LibeRation:::.nm_estimation_context(model, data)
  evaluator <- context$subjects[[1L]]
  before <- evaluator$tape_telemetry()
  value <- evaluator$objective(c(1, 20), 0, 0.2, 0.1)$value
  after <- evaluator$tape_telemetry()
  expect_true(is.finite(value))
  expect_equal(after$retapes, before$retapes + 1L)
  expect_equal(after$records, before$records + 1L)
})

test_that("fixed-effect priors contribute to deterministic and Bayesian fits", {
  fixture <- estimation_fixture()
  fixture$model$LIK_CONFIG$priors <- nm_prior("THETA1", "normal", mean = 2, sd = 0.1)
  fit <- nm_est(fixture$model, fixture$data, method = "LAPLACE", maxit = 2)
  expect_true(is.finite(fit$objective))
  expect_equal(fit$theta[[1]], 2)
})

test_that("population gradients include conditional-mode and curvature derivatives", {
  fixture <- estimation_fixture()
  fixture$model$THETAS$FIX <- c(FALSE, TRUE)
  fixture$model$SIGMAS$FIX <- FALSE
  fixture$model$OMEGAS$FIX <- FALSE
  context <- LibeRation:::.nm_estimation_context(fixture$model, fixture$data)
  map <- LibeRation:::.nm_outer_map(context$model)
  central <- function(fn, at, relative_step = 1e-4) {
    vapply(seq_along(at), function(index) {
      step <- relative_step * max(1, abs(at[[index]]))
      high <- low <- at
      high[[index]] <- high[[index]] + step
      low[[index]] <- low[[index]] - step
      (fn(high) - fn(low)) / (2 * step)
    }, numeric(1))
  }
  for (approximation in c("its", "foce", "focei", "laplace")) {
    objective <- LibeRation:::.nm_nested_objective(
      context, approximation, eta_maxit = 100L, tolerance = 1e-8
    )
    at <- map$start
    encoded_objective <- function(outer) objective(map$decode(outer))
    encoded_objective(at)
    if (approximation != "its") {
      mode <- attr(objective, "state")$modes[[1L]]$par
      taped_curvature <- context$subjects[[1L]]$curvature(
        map$decode(at)$theta, mode, map$decode(at)$sigma,
        map$decode(at)$omega, approximation, gradient = TRUE
      )
      reference_curvature <- LibeRation:::.nm_subject_curvature_logdet(
        context, context$subjects[[1L]], map$decode(at), mode, approximation
      )
      expect_equal(taped_curvature$value, reference_curvature, tolerance = 2e-11)
      expect_true(all(is.finite(taped_curvature$gradient)))
    }
    exact <- LibeRation:::.nm_nested_outer_gradient(
      context, map, objective, map$decode(at), approximation
    )
    numerical <- central(encoded_objective, at)
    expect_equal(
      unname(exact), unname(numerical), tolerance = 8e-5,
      info = approximation
    )
  }
})

test_that("native prior gradients and full-OMEGA transforms are exact", {
  fixture <- estimation_fixture()
  fixture$model$LIK_CONFIG$priors <- do.call(rbind, list(
    nm_prior("THETA1", "normal", mean = 1.5, sd = 0.7),
    nm_prior("SIGMA1", "lognormal", mean = -1, sd = 0.4),
    nm_prior("OMEGA1", "inverse_gamma", shape = 3, rate = 0.2)
  ))
  parameters <- list(theta = c(2, 20), sigma = 0.25, omega = 0.12)
  analytic <- LibeRation:::.nm_prior_nll_native_gradient(
    fixture$model, parameters
  )
  point <- c(parameters$theta, parameters$sigma, parameters$omega)
  fn <- function(value) LibeRation:::.nm_prior_nll(
    fixture$model,
    list(theta = value[1:2], sigma = value[3], omega = value[4])
  )
  numerical <- vapply(seq_along(point), function(index) {
    high <- low <- point
    high[[index]] <- high[[index]] + 1e-6
    low[[index]] <- low[[index]] - 1e-6
    (fn(high) - fn(low)) / 2e-6
  }, numeric(1))
  expect_equal(analytic, numerical, tolerance = 2e-6)

  omega <- data.frame(
    OMEGA = 1:3, ROW = c(1, 2, 2), COL = c(1, 1, 2),
    Value = c(0.09, 0.02, 0.16), FIX = FALSE
  )
  correlated <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "DV", "MDV"), ADVAN = 1,
    PRED = "CL=THETA(1)*exp(ETA(1)); V=THETA(2)*exp(ETA(2)); S1=V",
    ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(THETA = 1:2, Value = c(2, 20), FIX = TRUE),
    OMEGAS = omega,
    SIGMAS = data.frame(SIGMA = 1, Value = 0.2, FIX = TRUE)
  )
  map <- LibeRation:::.nm_outer_map(correlated)
  at <- map$start
  analytic <- map$jacobian(map$decode(at))
  omega_rows <- nrow(correlated$THETAS) + nrow(correlated$SIGMAS) +
    seq_len(nrow(correlated$OMEGAS))
  numerical <- vapply(seq_along(at), function(index) {
    high <- low <- at
    high[[index]] <- high[[index]] + 1e-6
    low[[index]] <- low[[index]] - 1e-6
    (map$decode(high)$omega - map$decode(low)$omega) / 2e-6
  }, numeric(nrow(correlated$OMEGAS)))
  expect_equal(
    unname(analytic[omega_rows, , drop = FALSE]), unname(numerical),
    tolerance = 2e-7
  )
})

test_that("THETA bounds use NONMEM-style defaults and constrain estimation", {
  fixture <- estimation_fixture()
  expect_equal(fixture$model$THETAS$LOWER, c(0.002, 0.02))
  expect_equal(fixture$model$THETAS$UPPER, c(2000, 20000))

  fixture$model$THETAS$FIX <- c(FALSE, TRUE)
  fixture$model$THETAS$LOWER[[1L]] <- 1.95
  fixture$model$THETAS$UPPER[[1L]] <- 2.05
  fit <- nm_est(fixture$model, fixture$data, method = "FOCEI", maxit = 8)
  expect_gte(fit$theta[[1L]], 1.95)
  expect_lte(fit$theta[[1L]], 2.05)
})

test_that("estimation can retain compiled subject workers and print gradients", {
  skip_on_cran()
  fixture <- estimation_fixture()
  fixture$model$THETAS$FIX <- c(FALSE, TRUE)
  output <- capture.output(fit <- nm_est(
    fixture$model, fixture$data, method = "FOCEI", maxit = 2,
    n_cores = 2L, print_every = 1L
  ))
  expect_s3_class(fit, "nm_fit")
  expect_true(is.finite(fit$objective))
  expect_true(any(grepl("SCALED GRADIENT", output, fixed = TRUE)))
})

test_that("estimation methods can run as one ordered sequence", {
  fixture <- estimation_fixture(fix = TRUE)
  fit <- nm_est_sequence(
    fixture$model, fixture$data,
    stages = list(
      nm_est_stage("FO", maxit = 1),
      nm_est_stage("FOCE", maxit = 1, eta_maxit = 10)
    )
  )
  expect_s3_class(fit, "nm_fit")
  expect_equal(fit$method_sequence, c("FO", "FOCE"))
  expect_equal(fit$sequence_label, "FO -> FOCE")
  expect_true(fit$sequence_complete)
  expect_length(fit$stages, 2L)
  expect_equal(fit$stages[[2L]]$stage$initial$theta, fit$stages[[1L]]$theta)
  expect_true(fit$stages[[2L]]$stage$initial$eta_warm_start)
  expect_true(is.finite(fit$timing$sequence_total_seconds))
})
