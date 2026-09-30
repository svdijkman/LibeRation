# Fit construction, deterministic estimators, and the public estimation API.
# Split from estimation.R as a behaviour-preserving source move.

.nm_fit_result <- function(context, method, parameters, objective, modes,
                           optimizer, diagnostics = list()) {
  eta <- if (context$n_eta) {
    do.call(rbind, lapply(modes, `[[`, "par"))
  } else matrix(numeric(), context$n_subjects, 0L)
  colnames(eta) <- if (context$n_eta) paste0("ETA", seq_len(context$n_eta)) else NULL
  iterations <- suppressWarnings(as.integer(optimizer$iterations %||% NA_integer_))
  if (!length(iterations) || is.na(iterations)) {
    iterations <- suppressWarnings(as.integer(optimizer$counts[["gradient"]] %||% NA_integer_))
  }
  if (!length(iterations) || is.na(iterations)) {
    iterations <- suppressWarnings(as.integer(optimizer$counts[["function"]] %||% NA_integer_))
  }
  objective_evaluations <- suppressWarnings(as.integer(
    optimizer$objective_evaluations %||% optimizer$counts[["function"]] %||% NA_integer_
  ))
  eta_iterations <- if (length(modes)) {
    sum(vapply(modes, function(mode) as.integer(mode$iterations %||% 0L), integer(1)))
  } else 0L
  eta_evaluations <- if (length(modes)) {
    sum(vapply(modes, function(mode) as.integer(mode$evaluations %||% 0L), integer(1)))
  } else 0L
  population_work <- optimizer$population_objective %||% list()
  if (eta_iterations == 0L && !is.null(population_work$mode_iterations)) {
    eta_iterations <- as.integer(population_work$mode_iterations)
  }
  if (eta_evaluations == 0L && !is.null(population_work$mode_evaluations)) {
    eta_evaluations <- as.integer(population_work$mode_evaluations)
  }
  tape <- lapply(context$subjects, function(evaluator) evaluator$tape_telemetry())
  diagnostics$optimizer <- list(
    backend = optimizer$backend %||% "unknown",
    coordinator = optimizer$coordinator %||% "r-callback",
    objective_backend = optimizer$objective_backend %||%
      "r-orchestrated-population-objective",
    elapsed_seconds = optimizer$elapsed_seconds %||% NA_real_,
    objective_initialization_seconds =
      optimizer$objective_initialization_seconds %||% 0,
    objective_evaluations = objective_evaluations,
    gradient_evaluations = optimizer$gradient_evaluations %||% NA_integer_,
    gradient_fallbacks = optimizer$gradient_fallbacks %||% 0L,
    gradient_fallback_evaluations =
      optimizer$gradient_fallback_evaluations %||% 0L,
    trace = optimizer$telemetry %||% NULL,
    population_objective = optimizer$population_objective %||% NULL
  )
  diagnostics$conditional_modes <- list(
    iterations = eta_iterations, evaluations = eta_evaluations,
    backends = table(vapply(modes, function(mode) mode$backend %||% "unknown", character(1)))
  )
  diagnostics$data_layout <- list(
    storage = context$subject_data_layout %||% "copied",
    subjects = context$n_subjects,
    persistent_subject_frames = sum(vapply(
      context$subjects, function(evaluator) !is.null(evaluator$data), logical(1)
    )),
    subject_frame_materializations = sum(vapply(
      context$subjects, function(evaluator) evaluator$data_materializations,
      integer(1)
    )),
    minimal_projections = sum(vapply(
      context$subjects, function(evaluator) evaluator$data_projections,
      integer(1)
    )),
    shared_rows = nrow(context$data)
  )
  diagnostics$tapes <- list(
    unique_structures = length(unique(vapply(tape, `[[`, character(1), "signature"))),
    subjects = length(tape),
    shared_prediction_tapes = max(
      0L, length(tape) - length(unique(vapply(tape, `[[`, character(1), "signature")))
    ),
    records = sum(vapply(tape, `[[`, integer(1), "records")) +
      as.integer(population_work$tape_records %||% 0L),
    retapes = sum(vapply(tape, `[[`, integer(1), "retapes")) +
      as.integer(population_work$tape_retapes %||% 0L),
    validity_checks = sum(vapply(tape, `[[`, integer(1), "validity_checks"))
  )
  objective_semantics <- diagnostics$objective_semantics %||% list(
    type = "negative_twice_estimator_objective",
    likelihood_comparable = !method %in% c("BAYES", "HMC", "NUTS")
  )
  objective_type <- as.character(
    objective_semantics$type %||% "negative_twice_estimator_objective"
  )[[1L]]
  objective_comparable <- isTRUE(
    objective_semantics$likelihood_comparable %||% FALSE
  )
  diagnostics$objective_semantics <- utils::modifyList(
    objective_semantics,
    list(type = objective_type, likelihood_comparable = objective_comparable)
  )
  structure(
    list(
      version = 1L, method = method, objective = as.numeric(objective),
      objective_type = objective_type,
      objective_comparable = objective_comparable,
      theta = parameters$theta, omega = parameters$omega,
      sigma = parameters$sigma, eta = eta,
      convergence = optimizer$convergence, message = optimizer$message,
      iterations = iterations,
      objective_evaluations = objective_evaluations,
      evaluations = optimizer$counts, model = context$model,
      data = context$data, diagnostics = diagnostics
    ),
    class = "nm_fit"
  )
}

.nm_est_fo <- function(context, map, maxit, tolerance, trace, print_every = 0L,
                       optimizer_backend = "auto") {
  if (context$model$LIK_CONFIG$blq_method != "none") {
    .nm_stop("FO does not support a censored Gaussian linearization; use FOCEI or LAPLACE for BLQ data.")
  }
  objective <- function(parameters) {
    .nm_fo_objective(context, parameters)
  }
  gradient <- function(parameters) .nm_fo_outer_gradient(context, map, parameters)
  compiled <- .nm_cpp_population_objective(
    context, map, "fo", eta_maxit = 100L, tolerance = tolerance
  )
  optimizer <- .nm_outer_optim(
    map, objective, maxit, tolerance, trace, print_every, gradient = gradient,
    optimizer_backend = optimizer_backend, compiled_objective = compiled,
    strict_convergence = TRUE
  )
  parameters <- map$decode(optimizer$par)
  modes <- .nm_subject_modes(
    context, parameters, maxit = 100L, tolerance = tolerance,
    exact_hessian = FALSE
  )
  .nm_fit_result(
    context, "FO", parameters, optimizer$value, modes, optimizer,
    diagnostics = list(population_gradient = "exact CppAD FO marginal gradient")
  )
}

.nm_est_nested <- function(context, map, method, maxit, eta_maxit, tolerance, trace,
                           print_every = 0L, optimizer_backend = "auto",
                           initial_eta = NULL) {
  approximation <- switch(
    method, FOCE = "foce", FOCEI = "focei", LAPLACE = "laplace", "laplace"
  )
  objective <- .nm_nested_objective(
    context, approximation, eta_maxit, tolerance, initial_eta = initial_eta
  )
  gradient <- function(parameters) .nm_nested_outer_gradient(
    context, map, objective, parameters, approximation
  )
  compiled <- .nm_cpp_population_objective(
    context, map, approximation, eta_maxit, tolerance, initial_eta = initial_eta
  )
  optimizer <- .nm_outer_optim(
    map, objective, maxit, tolerance, trace, print_every,
    gradient = gradient, optimizer_backend = optimizer_backend,
    compiled_objective = compiled
  )
  parameters <- map$decode(optimizer$par)
  cached <- if (!is.null(compiled$pointer)) tryCatch(
    .liberation_population_objective_state(compiled$pointer, optimizer$par),
    error = function(error) NULL
  ) else NULL
  modes <- cached$modes %||% .nm_subject_modes(
    context, parameters, starts = initial_eta,
    maxit = eta_maxit, tolerance = tolerance,
    interaction = approximation != "foce",
    exact_hessian = approximation == "laplace"
  )
  .nm_assert_final_conditional_modes(modes, method)
  work <- optimizer$population_objective
  if (is.null(work)) {
    work <- list(
      value_requests = attr(objective, "state")$objective_calls,
      value_cache_hits = attr(objective, "state")$cache_hits,
      mode_iterations = attr(objective, "state")$mode_iterations,
      mode_evaluations = attr(objective, "state")$mode_evaluations
    )
  }
  .nm_fit_result(
    context, method, parameters, optimizer$value, modes, optimizer,
    diagnostics = list(
      eta_convergence = vapply(modes, `[[`, integer(1), "convergence"),
      eta_jitter = vapply(modes, `[[`, numeric(1), "jitter"),
      approximation = approximation,
      conditional_mode_work = list(
        objective_calls = work$value_requests %||% work$parameter_evaluations,
        cache_hits = sum(
          work$value_cache_hits %||% 0L,
          work$gradient_cache_hits %||% 0L,
          work$shared_state_hits %||% 0L
        ),
        iterations = work$mode_iterations,
        evaluations = work$mode_evaluations
      ),
      population_gradient = "exact CppAD curvature with implicit conditional-mode derivative"
    )
  )
}

#' Estimate a population pharmacometric model
#'
#' Deterministic conditional methods use exact CppAD gradients for ETA modes.
#' LAPLACE uses the exact conditional Hessian; FOCE/FOCEI use the
#' interaction-aware Gauss-Newton curvature; FO integrates the first-order
#' Gaussian linearization analytically. ITS alternates conditional-mode and
#' first-order conditional-variance calculations with approximate single-step
#' population updates. GQ approximates the ETA
#' marginal integral with finite adaptive tensor or Smolyak sparse
#' Gauss--Hermite quadrature. IMP uses importance-sampling Monte-Carlo EM;
#' each E-step calculates normalized importance weights and its M-step optimizes
#' the corresponding complete-data expectation. SAEM applies the
#' Robbins--Monro recurrence to that complete-data auxiliary function.
#' HMC and NUTS use exact joint CppAD gradients with unconstrained parameter
#' transforms and dual-averaged step-size adaptation. Optimized NUTS uses
#' multinomial progressive trajectory sampling, a generalized U-turn check,
#' and windowed regularized diagonal metric adaptation; the retained classic
#' slice variant is available for comparison. NPML estimates a discrete
#' fixed-support maximum-likelihood distribution. NPAG adds explicit support
#' expansion and pruning around the same constrained mixture likelihood.
#' Serial FO, FOCE, FOCEI, Laplace, and ITS fits expose a persistent C++
#' population objective through thin callbacks to R's L-BFGS-B/BFGS optimizer;
#' PSOCK fits retain R coordination across their persistent C++ workers.
#'
#' @param model An [nm_model()] or compiled [NMEngine].
#' @param data A NONMEM-style event dataset.
#' @param method Estimation method.
#' @param maxit Maximum outer evaluations.
#' @param eta_maxit Maximum conditional ETA iterations.
#' @param tolerance Relative optimization tolerance.
#' @param trace Optimizer tracing level.
#' @param print_every Write the objective and scaled population gradient to
#'   stdout every N outer evaluations. Zero disables iteration logging.
#' @param n_cores Number of parallel subject workers. Parallel workers persist
#'   compiled C++ engines for the duration of the estimation.
#' @param optimizer_backend Use the experimental compiled scaled, bounded BFGS
#'   optimizer (`"native"`), R's mature L-BFGS-B/BFGS driver (`"r"`), or the
#'   calibrated production policy (`"auto"`). Auto selects R's driver while
#'   its objective, gradient, conditional-mode state, parameter transforms, and
#'   covariance derivatives remain in the persistent C++ population evaluator.
#' @param covariance Run a covariance step after estimation.
#' @param covariance_type Covariance estimator: robust sandwich (default),
#'   automatic R/S fallback,
#'   objective Hessian (`"hessian"`), subject-score OPG (`"opg"`), or robust
#'   sandwich (`"sandwich"`). `"r"` and `"s"` are accepted aliases.
#' @param covariance_tolerance Positive-definite regularization tolerance for
#'   the covariance step.
#' @param covariance_samples Target integration budget for an IMP or SAEM
#'   covariance step. Low-dimensional ETA integrals use tensor Gauss--Hermite
#'   quadrature near this budget; higher-dimensional integrals use random-normal
#'   importance samples. The default reuses the IMP count or uses 200 for SAEM.
#' @param covariance_seed Common-random-number seed for the random-normal
#'   covariance fallback. The default reuses the estimation seed.
#' @param initial_eta Optional finite subject-by-ETA matrix used to warm-start
#'   compatible conditional or SAEM estimation steps.
#' @param collect_output Whether selected generated OUTPUT columns should be
#'   evaluated and retained with the completed fit.
#' @param allow_fd_gradient Permit an explicit finite-difference replacement
#'   if an advertised population gradient is non-finite. The default is
#'   `FALSE`; enabling it records and warns about every fallback.
#' @param audit_artifacts Generate an opt-in NONMEM-style audit bundle containing
#'   a listing, control-stream representation, estimates, covariance/correlation
#'   matrices when available, subject ETAs, and saved output tables.
#' @param ... Method-specific controls. ITS accepts `its_mstep_maxit`,
#'   `its_mstep_schedule` (`"fixed"` or `"progressive"`), and
#'   `its_acceleration` (`"none"` or safeguarded `"aitken"`).
#'   `its_eta_schedule` optionally relaxes early conditional-mode tolerances
#'   by `its_eta_tolerance_multiplier` before returning to the requested
#'   tolerance for the final distribution. The
#'   compatibility policy fixes the established single-step schedule without
#'   acceleration or tolerance relaxation; the optimized policy progressively
#'   increases M-step effort and accepts an accelerated point only when its
#'   auxiliary objective improves.
#'   IMP accepts `imp_algorithm` (`"mcem"`, the default, or the explicitly
#'   labelled `"marginal_ml"` finite-common-random-number alternative),
#'   `n_imp`, `seed`, and `imp_mstep_maxit`. Optimized MCEM additionally
#'   supports progressive sample/M-step schedules, antithetic or randomized
#'   quasi-Monte-Carlo (`"rqmc"`) draws, optional ESS-directed per-subject allocation,
#'   Gaussian, Student-t, or defensive-mixture proposals, guarded proposal-mode
#'   reuse, and stationarity stopping through the `imp_*` controls.
#'   `imp_proposal_curvature = "auto"` uses a Fisher/Gauss--Newton proposal
#'   curvature in optimized mode and the exact conditional Hessian in
#'   compatibility mode; exact importance weights preserve the MCEM target.
#'   The optimized path can also keep the weighted M-step and closed-form
#'   SIGMA/OMEGA updates in its persistent native context. The
#'   compatibility policy retains fixed sample/M-step counts, random Gaussian
#'   proposals, exact proposal curvature, fresh conditional modes, and its
#'   established optimizer trajectory.
#'   For `method = "GQ"`, use `gq_grid`
#'   (`"auto"`, `"tensor"`, or `"smolyak"`), `gq_order` (tensor nodes per ETA
#'   dimension, default 5), `gq_level` (Smolyak level, default 3),
#'   `gq_adaptive` (default `TRUE`), `gq_max_points` (retained-grid allocation
#'   guard, default 100000), and `gq_gradient` (`"score"`, `"finite_grid"`, or
#'   `"auto"`). A score-based search is always finished against the complete
#'   finite adaptive-grid objective before the estimate is returned. Automatic
#'   selection uses tensor quadrature through three ETAs and Smolyak quadrature
#'   for higher-dimensional models.
#'   For `method = "HMC"` or `"NUTS"`, controls include `n_warmup` (500),
#'   `n_sample` (1000 per chain), `n_thin` (1), `n_chains` (4), `seed`,
#'   optional `step_size`, `target_acceptance` (0.8), `adapt_mass` (`TRUE`),
#'   `n_leapfrog` (10; HMC), `max_depth` (10; NUTS), `nuts_variant`
#'   (`"auto"`, `"classic_slice"`, or `"multinomial"`), and
#'   `divergence_threshold` (1000). `hmc_metric` selects `"unit"`,
#'   `"diagonal"`, `"dense"`, or `"block"` Euclidean geometry; optional
#'   one-based `hmc_metric_blocks` override the default population/subject
#'   partition. Warmup uses initial fast, expanding slow metric windows, and
#'   terminal fast phases controlled by `adapt_initial_buffer` (75),
#'   `adapt_first_window` (25), and `adapt_terminal_buffer` (50).
#'   `initialization` selects the model point, bounded random jitter, or the
#'   robust automatic sequence, with `initialization_radius` and
#'   `initialization_attempts` controlling the latter. Per-chain energy and
#'   E-BFMI are retained alongside divergences and tree-depth diagnostics.
#'   `sampler_backend = "native"` keeps complete
#'   trajectories in C++; `"r"` retains the slower reference implementation
#'   for numerical comparison. `geometry = "auto"` uses a whitened,
#'   MU-aware non-centred ETA parameterization in `liber_optimized` mode when
#'   every ETA is MU-referenced and the base positive-definite OMEGA layout is
#'   available. `"centered"` forces the exact established target;
#'   `"mu_noncentered"` requires eligibility rather than silently falling back.
#'   The likelihood tape is unchanged and the Cholesky transform, Jacobian, and
#'   population/random-effect chain rules are differentiated exactly in both
#'   sampler backends. For `method = "NPML"` or `"NPAG"`, use
#'   `np_supports` for an optional fixed starting matrix, `np_points` (25),
#'   `np_max_support` (100), method-specific `np_min_weight` (0 for fixed-support
#'   NPML and 1e-5 for adaptive-grid NPAG), `np_weight_maxit` (1000),
#'   `np_cycles` (3), and, for NPAG, `np_grid_step` (1), `np_grid_decay`
#'   (0.5), and `np_max_candidates` (500). `np_estimate_population` controls
#'   alternating THETA/SIGMA updates. NPML ignores a nonzero pruning threshold
#'   with a warning; select NPAG when support adaptation is intended.
#'   Ordinary covariance is not regular for a
#'   discrete support distribution; use bootstrap uncertainty for NPML/NPAG.
#'   For `method = "SAEM"`, stochastic approximation is applied to a retained
#'   complete-data auxiliary-function state rather than to a single sampled
#'   conditional objective. `saem_kernel = "auto"` retains random-walk
#'   Metropolis in compatibility mode and selects a robust Laplace/Student-t
#'   independence Metropolis kernel (f-SAEM) for eligible ETA
#'   models in `liber_optimized` mode. `fsaem_refresh` (25) controls how often
#'   conditional modes and curvature are refreshed; `fsaem_eta_maxit` (50)
#'   bounds that calculation. `fsaem_distribution` selects `"student_t"`
#'   (the optimized automatic default) or `"gaussian"`; `fsaem_df` (7)
#'   controls Student-t tail weight. The optimized kernel also refreshes when the
#'   population point moves materially (`fsaem_parameter_refresh`, 0.15) or
#'   independence acceptance is poor (`fsaem_low_acceptance`, 0.1), and mixes
#'   in an exact OMEGA-scaled random-walk rescue kernel with probability
#'   `fsaem_rescue_probability` (0.1). Set that probability to zero and both
#'   adaptive thresholds conservatively to recover the fixed-refresh form.
#'   `saem_kernel = "random_walk"` is the explicit established-kernel comparator.
#'   SAEM always reports parameter/objective stationarity over
#'   `stationarity_window` (20) iterations using `stationarity_tolerance`
#'   (1e-3). `auto_stop` defaults to `FALSE` in compatibility mode and `TRUE`
#'   in optimized mode;
#'   `auto_stop_consecutive` and `auto_stop_min_iterations` make confirmation
#'   explicit. `saem_support_max` and `saem_support_prune` optionally bound the
#'   retained optimized Q support; both default to zero so the exact
#'   Robbins--Monro mixture is retained. Optimized execution can thin expensive
#'   numerical M-steps with `saem_mstep_interval_burn` and
#'   `saem_mstep_interval` while continuing every stochastic-approximation
#'   update, and uses post-burn Polyak parameter averaging by default;
#'   `saem_parameter_averaging = "none"` disables it and
#'   `saem_average_start` sets its first iteration. Compatibility mode forces
#'   both M-step intervals to one and disables averaging. `n_replicates` runs
#'   independently seeded sequential replicates, separated by
#'   `replicate_seed_stride`. Replicates are ranked with a common-seed marginal
#'   importance-sampling score rather than their incompatible stochastic
#'   auxiliary-function histories; `replicate_score_samples` (200),
#'   `replicate_score_seed`, and `replicate_score_eta_maxit` control that score.
#'   Use separate queued jobs for replicate-level parallelism.
#'   For `method = "BAYES"`, `outer_kernel = "auto"` retains the isotropic
#'   population random walk in compatibility mode and learns a full proposal
#'   covariance during burn-in in `liber_optimized` mode. The adaptive controls
#'   are `adaptive_start` (50), `adaptive_interval` (10), and optional
#'   `adaptive_target` (0.44 in one dimension, 0.234 otherwise). Explicit
#'   `outer_kernel = "isotropic"` preserves the former proposal. Eligible
#'   serial compatibility models run that established isotropic/random-walk
#'   algorithm in the persistent C++ coordinator while retaining R's RNG and
#'   proposal order; compatibility MU interweaving, subject parallelism, and
#'   iteration printing retain the R coordinator. Eligible optimized serial
#'   models additionally keep adaptation and Laplace/Student-t ETA sweeps in
#'   C++, including MU, IOV/general random effects, and guarded ODE retaping;
#'   `options(LibeRation.bayes_native_coordinator = FALSE)` retains the
#'   reference R coordinator for comparison. `eta_kernel = "auto"` uses the
#'   inexpensive OMEGA-scaled random walk for one ETA and a Laplace-Gaussian
#'   ETA proposal for multivariate ETAs; `"student_t"` plus `bayes_eta_df`
#'   provides a robust heavy-tailed alternative. Optimized BAYES also supports
#'   `delayed_rejection_scale` (0.25), independent `n_chains`, rank-normalized
#'   split R-hat, bulk/tail ESS, and mean Monte Carlo SE.
#'   `bayes_gibbs_omega = TRUE` uses exact inverse-gamma Gibbs updates when a
#'   free diagonal OMEGA has an explicit conjugate prior and the random-effect
#'   layout is eligible; all other OMEGA structures retain Metropolis updates.
#'   Subject-parallel models retain the exact R coordinator with PSOCK workers;
#'   independent chains can be distributed as separate queue jobs.
#' @param numerical_mode Numerical policy. `"nonmem_compatibility"` is the
#'   conservative default and uses the defining NONMEM method semantics and
#'   matched control interpretation for methods with a NONMEM counterpart;
#'   implementation-specific random streams and undocumented internal stopping
#'   details are not claimed to be bitwise identical. `"liber_optimized"`
#'   preserves each estimator's defining target and update equations while
#'   enabling validated LibeR-specific proposals, caching, batching, and solver
#'   accelerations. Explicit estimator controls take precedence in either mode.
#'   IMP, GQ, SAEM, and BAYES accept `mu_specialization = TRUE` (the default).
#'   IMP then re-centres cached conditional-mode starts as MU values change.
#'   When an affine MU design varies between subjects (for example an estimated
#'   covariate coefficient), IMP uses the score path as a warm start and
#'   automatically refines it against the exact finite common-random-number
#'   objective.
#'   GQ reuses the same execution-local conditional-state cache. In optimized
#'   mode, compiled FOCE/FOCEI/Laplace objectives also preserve `MU + ETA`
#'   while warm-starting modes as MU values move; compatibility mode is
#'   unchanged.
#'   Eligible affine MU models use a generalized least-squares fixed-effect
#'   M-step in SAEM, with vectorized OMEGA-keyed normal-equation caching and a
#'   closed-form-only fast path, and a Metropolis-corrected Gaussian MU block
#'   in BAYES.
#'   Non-affine, rank-deficient, or otherwise ineligible models automatically
#'   retain the ordinary estimator path and record the reason in
#'   `fit$diagnostics$mu_specialization`.
#' @export
nm_est <- function(model, data,
                   method = c("FOCEI", "FOCE", "FO", "LAPLACE", "ITS",
                              "GQ", "IMP", "SAEM", "BAYES", "HMC", "NUTS",
                              "NPML", "NPAG"),
                   maxit = 200L, eta_maxit = 100L, tolerance = 1e-6,
                   trace = 0L, print_every = 0L, n_cores = 1L,
                   optimizer_backend = c("auto", "native", "r"),
                   covariance = FALSE,
                   covariance_type = c("sandwich", "auto", "hessian", "opg", "r", "s"),
                   covariance_tolerance = 1e-8,
                   covariance_samples = NULL, covariance_seed = NULL,
                   initial_eta = NULL, collect_output = TRUE,
                   allow_fd_gradient = FALSE, audit_artifacts = FALSE,
                   numerical_mode = NULL, ...) {
  request_started <- proc.time()[["elapsed"]]
  method <- match.arg(method)
  model <- .nm_model_with_numerical_mode(model, numerical_mode)
  optimizer_backend <- match.arg(optimizer_backend)
  covariance_type <- match.arg(covariance_type)
  if (length(allow_fd_gradient) != 1L || is.na(allow_fd_gradient)) {
    .nm_stop("`allow_fd_gradient` must be TRUE or FALSE.")
  }
  if (length(audit_artifacts) != 1L || is.na(audit_artifacts)) {
    .nm_stop("`audit_artifacts` must be TRUE or FALSE.")
  }
  audit_artifacts <- isTRUE(audit_artifacts)
  previous_fd_option <- options(
    LibeRation.allow_fd_gradient = isTRUE(allow_fd_gradient)
  )
  on.exit(options(previous_fd_option), add = TRUE)
  print_every <- as.integer(print_every)
  if (length(print_every) != 1L || is.na(print_every) || print_every < 0L) {
    .nm_stop("`print_every` must be a non-negative integer.")
  }
  if (length(covariance) != 1L || is.na(covariance)) {
    .nm_stop("`covariance` must be TRUE or FALSE.")
  }
  covariance <- isTRUE(covariance)
  if (covariance && method %in% c("BAYES", "HMC", "NUTS")) {
    .nm_stop(method, " reports posterior SDs and credible intervals automatically; a frequentist covariance step is not applicable.")
  }
  if (covariance && !method %in% c("FO", "FOCE", "FOCEI", "LAPLACE", "ITS", "GQ", "IMP", "SAEM")) {
    .nm_stop("Covariance is available for FO, FOCE, FOCEI, LAPLACE, ITS, GQ, IMP, and SAEM fits.")
  }
  if (!inherits(model, c("nm_model", "NMEngine"))) {
    .nm_stop("`model` must be an nm_model or NMEngine.")
  }
  if (missing(data)) .nm_stop("`data` is required.")
  model_definition <- if (inherits(model, "NMEngine")) model$model else model
  if (identical(model_definition$LIK_CONFIG$error, "likelihood") &&
      method %in% c("FO", "FOCE", "FOCEI")) {
    .nm_stop(
      method, " assumes a Gaussian residual linearization and cannot be used ",
      "with a user-defined likelihood. Use LAPLACE for NONMEM-like conditional ",
      "likelihood estimation, or ITS/GQ/IMP/SAEM/BAYES/HMC/NUTS/NPML/NPAG."
    )
  }
  context_started <- proc.time()[["elapsed"]]
  context <- .nm_estimation_context(model, data, n_cores = n_cores, method = method)
  if (!is.null(initial_eta)) {
    initial_eta <- as.matrix(initial_eta)
    expected <- c(context$n_subjects, context$n_eta)
    if (!identical(dim(initial_eta), expected) || any(!is.finite(initial_eta))) {
      .nm_stop("`initial_eta` must be a finite ", expected[[1L]], " x ",
               expected[[2L]], " matrix for this dataset.")
    }
  }
  if (!is.null(context$parallel)) {
    on.exit(try(parallel::stopCluster(context$parallel$cluster), silent = TRUE),
            add = TRUE)
  }
  map <- .nm_outer_map(context$model)
  context_initialization_seconds <- unname(
    proc.time()[["elapsed"]] - context_started
  )
  estimation_started <- proc.time()[["elapsed"]]
  fit <- if (method == "FO") {
    .nm_est_fo(
      context, map, maxit, tolerance, trace, print_every, optimizer_backend
    )
  } else if (method %in% c("FOCE", "FOCEI", "LAPLACE")) {
    .nm_est_nested(
      context, map, method, maxit, eta_maxit, tolerance, trace, print_every,
      optimizer_backend, initial_eta = initial_eta
    )
  } else {
    .nm_est_stochastic(
      context, map, method, maxit = maxit, eta_maxit = eta_maxit,
      tolerance = tolerance, trace = trace, print_every = print_every,
      optimizer_backend = optimizer_backend, initial_eta = initial_eta, ...
    )
  }
  model_fit_wall_seconds <- unname(
    proc.time()[["elapsed"]] - estimation_started
  )
  objective_initialization_seconds <- as.numeric(
    fit$diagnostics$optimizer$objective_initialization_seconds %||% 0
  )
  model_fit_seconds <- max(
    0, model_fit_wall_seconds - objective_initialization_seconds
  )
  initialization_seconds <- context_initialization_seconds +
    objective_initialization_seconds
  fit$diagnostics$eta_maxit <- as.integer(eta_maxit)
  fit$diagnostics$tolerance <- as.numeric(tolerance)
  fit$diagnostics$numerical_mode <- model_definition$NUMERICAL_MODE %||%
    "nonmem_compatibility"
  fit$numerical_mode <- fit$diagnostics$numerical_mode
  covariance_seconds <- NA_real_
  if (covariance) {
    covariance_started <- proc.time()[["elapsed"]]
    attr(fit, ".estimation_context") <- context
    fit$covariance <- tryCatch(
      nm_cov_step(
        fit, type = covariance_type, tolerance = covariance_tolerance,
        samples = covariance_samples, seed = covariance_seed,
        eta_maxit = eta_maxit
      ),
      error = function(error) structure(list(
        status = "failed", type = covariance_type,
        error = conditionMessage(error)
      ), class = "nm_covariance_error")
    )
    attr(fit, ".estimation_context") <- NULL
    covariance_seconds <- unname(proc.time()[["elapsed"]] - covariance_started)
  }
  fit$timing <- list(
    initialization_seconds = as.numeric(initialization_seconds),
    model_fit_seconds = as.numeric(model_fit_seconds),
    covariance_seconds = as.numeric(covariance_seconds),
    total_seconds = as.numeric(model_fit_seconds +
      if (is.finite(covariance_seconds)) covariance_seconds else 0),
    wall_total_seconds = as.numeric(proc.time()[["elapsed"]] - request_started),
    context_cache_hit = isTRUE(context$cache_hit)
  )
  if (isTRUE(collect_output) && length(fit$model$OUTPUT %||% character())) {
    fit$output <- .nm_fit_selected_outputs(fit)
  }
  if (audit_artifacts) {
    fit <- .nm_attach_audit_artifacts(
      fit, fit$model, fit$data, "estimate",
      details = list(
        method = method, maxit = maxit, eta_maxit = eta_maxit,
        covariance = covariance, n_cores = n_cores,
        numerical_mode = fit$numerical_mode
      )
    )
  }
  fit
}

#' @export
print.nm_fit <- function(x, ...) {
  cat("LibeRation fit\n")
  cat(
    "  method:", .nm_fit_method_label(x),
    " reported objective:", format(x$objective),
    if (isFALSE(x$objective_comparable)) " (not likelihood-comparable)" else "",
    " convergence:", x$convergence, "\n"
  )
  invisible(x)
}
