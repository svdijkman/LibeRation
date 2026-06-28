test_that("NONMEM theta (lower, initial) parsing uses initial estimate", {
  lines <- c(" (0, 3)   ; CL", " (0, 20)  ; VC")
  th <- .nm_parse_theta_block(lines)
  expect_equal(th$Value, c(3, 20))
  expect_equal(th$Lower, c(0, 0))
})

test_that("fix mask length matches parameter vector", {
  sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
  mask <- .nm_fix_mask(sim$model)
  par <- .nm_init_par(sim$model)
  expect_length(mask, length(par))
  expect_no_warning({
    mask | grepl("^OMEGA", .nm_par_labels(sim$model))
  })
})

test_that("nm_nll and FO estimation run on synthetic THEO data", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 3L, seed = 42L)
  nll0 <- nm_nll(
    sim$model, sim$data,
    sim$model$THETAS$Value,
    sim$model$OMEGAS$Value,
    sim$model$SIGMAS$Value
  )
  expect_true(is.finite(nll0) && nll0 > 0)

  fit_fo <- nm_est(sim$model, sim$data, method = "FO", backend = "R",
                   grad = "numeric", control = list(maxit = 30))
  expect_s3_class(fit_fo, "nm_fit")
  expect_equal(fit_fo$method, "FO")
  expect_true(all(is.finite(fit_fo$theta)))
  expect_true(any(abs(fit_fo$omega - sim$model$OMEGAS$Value) > 1e-6))
  expect_true(!is.null(fit_fo$fo_estimation$step2))
})

test_that("FOCE, SAEM, and Laplace estimators run", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 3L, seed = 43L)

  fit_foce <- nm_est(sim$model, sim$data, method = "FOCE", backend = "cpp",
                     grad = "numeric", control = list(maxit = 20), max_outer = 3L)
  expect_s3_class(fit_foce, "nm_fit")
  expect_equal(fit_foce$method, "FOCE")
  expect_true(is.matrix(fit_foce$eta))
  expect_true(is.finite(fit_foce$objective))
  expect_equal(
    fit_foce$objective,
    .nm_fit_inference_objective(fit_foce, fit_foce$par, data = sim$data),
    tolerance = 1e-3
  )

  fit_saem <- nm_est(sim$model, sim$data, method = "SAEM", backend = "R",
                     grad = "numeric", n_iter = 5L, n_burn = 1L, seed = 1L)
  expect_s3_class(fit_saem, "nm_fit")
  expect_equal(fit_saem$method, "SAEM")

  fit_lap <- nm_est(sim$model, sim$data, method = "LAPLACE", backend = "R",
                    grad = "numeric", control = list(maxit = 20), n_quad = 3L)
  expect_s3_class(fit_lap, "nm_fit")
  expect_equal(fit_lap$method, "LAPLACE")
  expect_equal(fit_lap$n_quad, 3L)
  expect_true(is.matrix(fit_lap$eta))
  expect_equal(nrow(fit_lap$eta), length(unique(sim$data$data$ID)))
})

test_that("FOCE grad=cpp uses C++ population gradients (Shiny path)", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 100L, seed = 1L)
  expect_true(LibeRation:::.nm_use_cpp_pop_grad(sim$model, "cpp"))
  fit <- nm_est(
    sim$model, sim$data, method = "FOCE",
    grad = "cpp", pk_engine = "cpp",
    control = list(maxit = 25, n_cores = 1L, compute_inference = TRUE),
    max_outer = 3L
  )
  expect_s3_class(fit, "nm_fit")
  expect_equal(fit$grad, "cpp")
  expect_equal(fit$grad_backend, "cpp")
  expect_true(is.finite(fit$objective))
})

test_that("Laplace fit stores post-hoc eta for predict", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 10L, seed = 43L)
  fit_foce <- nm_est(
    sim$model, sim$data, method = "FOCE", engine = "cpp",
    grad = "numeric", control = list(maxit = 25), max_outer = 3L
  )
  fit_lap <- nm_est(
    sim$model, sim$data, method = "LAPLACE", engine = "cpp",
    grad = "numeric", control = list(maxit = 25), n_quad = 5L
  )
  expect_false(is.null(fit_lap$eta))
  expect_true(is.matrix(fit_lap$eta))
  pred <- predict(fit_foce, type = "ipred")
  obs <- pred[pred$MDV == 0L & pred$EVID == 0L, ]
  expect_true(any(abs(obs$IPRED - obs$PRED) > 1e-6))
})

test_that("nm_task simulation produces observations", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 2L, seed = 44L)
  out <- nm_task("sim", sim$model, sim$data, seed = 2L)
  expect_true("IPRED" %in% names(out))
  expect_true("REP" %in% names(out))
  expect_equal(out$REP[[1L]], 1L)
  expect_true(any(is.finite(out$DV)))
})

test_that("nm_simulate stacks replicates with REP column", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 2L, seed = 46L)
  sims <- nm_simulate(sim$model, sim$data, n_sim = 3L, seed = 1L)
  expect_length(sims, 3L)
  expect_equal(sims[[2L]]$data$REP[[1L]], 2L)
  packed <- LibeRation:::.nm_sim_pack_output(sims)
  expect_true(!is.null(packed$combined))
  comb <- packed$combined$data
  expect_true("REP" %in% names(comb))
  expect_setequal(unique(comb$REP), 1:3)
})

test_that("C++ PK and likelihood match R path", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 2L, seed = 45L)
  expect_equal(LibeRation:::.nm_n_transit(sim$model), 0L)
  th <- sim$model$THETAS$Value
  om <- sim$model$OMEGAS$Value
  sg <- sim$model$SIGMAS$Value
  nll_r <- nm_nll(sim$model, sim$data, th, om, sg, pk_engine = "R")
  nll_cpp <- nm_nll(sim$model, sim$data, th, om, sg, pk_engine = "cpp")
  expect_true(is.finite(nll_cpp))
  expect_true(is.finite(nll_r))
  expect_equal(nll_r, nll_cpp, tolerance = 1e-4)
  dat <- sim$data$data
  sub <- dat[dat$ID == dat$ID[1]]
  pred_r <- .nm_subject_ipred(
    sim$model, sub, th, om, c(0, 0, 0), sg, pk_engine = "R"
  )
  pred_cpp <- .nm_subject_ipred(
    sim$model, sub, th, om, c(0, 0, 0), sg, pk_engine = "cpp"
  )
  expect_equal(length(pred_r$F), length(pred_cpp$F))
})

test_that("grad and engine options work", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 2L, seed = 46L)

  fit_ad <- nm_est(sim$model, sim$data, method = "FO", grad = "auto",
                   pk_engine = "cpp", control = list(maxit = 8))
  expect_equal(fit_ad$grad, "ad")
  expect_equal(fit_ad$grad_requested, "auto")
  expect_equal(fit_ad$grad_backend, "cpp")

  fit_cpp <- nm_est(sim$model, sim$data, method = "FO", grad = "cpp",
                    pk_engine = "cpp", control = list(maxit = 8))
  expect_equal(fit_cpp$grad, "cpp")
  expect_equal(fit_cpp$grad_backend, "cpp")
  expect_true(any(abs(fit_cpp$theta - sim$model$THETAS$Value) > 1e-6))

  fit_cpp <- nm_est(sim$model, sim$data, method = "LAPLACE", engine = "cpp",
                    grad = "numeric", n_quad = 3L, control = list(maxit = 10))
  expect_equal(fit_cpp$engine, "cpp")

  fit_lap_cpp <- nm_est(sim$model, sim$data, method = "LAPLACE",
                        grad = "auto", n_quad = 3L, control = list(maxit = 12))
  expect_equal(fit_lap_cpp$grad, "cpp")
  expect_equal(fit_lap_cpp$grad_backend, "cpp")
  expect_false(fit_lap_cpp$laplace_mode_centered)

  gcpp <- .nm_laplace_nll_grad_cpp(
    sim$model, sim$data,
    sim$model$THETAS$Value, sim$model$OMEGAS$Value, sim$model$SIGMAS$Value,
    .nm_gh_nodes(3L)
  )
  expect_true(is.finite(gcpp$objective))
  expect_equal(length(gcpp$gradient), length(.nm_init_par(sim$model)))
  expect_true(any(abs(gcpp$gradient) > 0))

  fit_lap_ad_grad <- local({
    gh <- .nm_gh_nodes(3L)
    par0 <- .nm_init_par(sim$model)
    pn <- .nm_par_labels(sim$model)
    obj <- .nm_build_laplace_objective(sim$model, sim$data, gh, "cpp")
    g <- .nm_grad_population(obj, par0, pn, "ad", "cpp")
    g
  })
  expect_true(all(is.finite(fit_lap_ad_grad)))
  expect_true(any(abs(fit_lap_ad_grad) > 0))

  fit_lap_num <- nm_est(sim$model, sim$data, method = "LAPLACE", engine = "cpp",
                        grad = "numeric", pk_engine = "cpp", n_quad = 3L,
                        control = list(maxit = 40))
  fit_lap_cpp <- nm_est(sim$model, sim$data, method = "LAPLACE",
                        grad = "cpp", pk_engine = "cpp", n_quad = 3L,
                        control = list(maxit = 40))
  expect_equal(fit_lap_num$engine, "cpp")
  expect_equal(fit_lap_cpp$grad, "cpp")
  expect_equal(fit_lap_num$objective, fit_lap_cpp$objective, tolerance = 1.0)

  fit_saem_cpp <- nm_est(sim$model, sim$data, method = "SAEM", engine = "cpp",
                         grad = "numeric", n_iter = 3L, n_burn = 1L, seed = 2L)
  expect_equal(fit_saem_cpp$engine, "cpp")
  expect_true(is.list(fit_saem_cpp$optim) && length(fit_saem_cpp$optim) == 3L)
  counts <- .nm_sum_optim_counts(fit_saem_cpp$optim)
  expect_true(is.finite(counts[["function"]]) && counts[["function"]] > 0L)
})

test_that("nm_est attaches par_se and FOCE eta at final estimates", {
  skip_if_not_installed("data.table")
  sim <- nm_synthetic_theo(n_sub = 3L, seed = 42L)
  fit <- nm_est(
    sim$model, sim$data,
    method = "FOCE",
    grad = "numeric",
    pk_engine = "cpp",
    control = list(maxit = 25L, compute_inference = TRUE, infer_hessian = "numeric")
  )
  expect_true(is.matrix(fit$eta) && nrow(fit$eta) == 3L)
  expect_true(any(abs(fit$eta) > 1e-6))
  expect_false(is.null(fit$par_se))
  expect_true(any(is.finite(unname(fit$par_se))))
  expect_false(is.null(fit$par_grad))
})
