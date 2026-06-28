#' @keywords internal
.nm_tweak_inits <- function(model, par, tweak = TRUE) {
  if (isFALSE(tweak) || is.null(tweak)) {
    return(.nm_apply_fix(model, par))
  }
  factor <- if (isTRUE(tweak)) 1.05 else as.numeric(tweak[1L])
  if (!is.finite(factor) || factor <= 0) {
    return(.nm_apply_fix(model, par))
  }
  free <- which(!.nm_fix_mask(model))
  if (length(free) == 0L) {
    return(.nm_apply_fix(model, par))
  }
  par2 <- par
  delta <- abs(factor - 1)
  for (i in free) {
    dir <- sample(c(-1, 1), 1L)
    par2[i] <- par[i] * (1 + dir * delta)
  }
  .nm_apply_fix(model, par2)
}

#' @keywords internal
.nm_est_needs_retry <- function(fit) {
  if (is.null(fit) || identical(fit$method, "BAYES")) {
    return(FALSE)
  }
  conv <- fit$convergence
  is.null(conv) || !identical(as.integer(conv), 0L)
}

#' @keywords internal
.nm_finite_obj <- function(val) {
  if (length(val) != 1L || !is.numeric(val) || !is.finite(val)) {
    .Machine$double.xmax
  } else {
    as.numeric(val)
  }
}

#' @keywords internal
.nm_finite_grad <- function(g) {
  if (is.null(g) || !is.numeric(g)) {
    return(g)
  }
  bad <- !is.finite(g)
  if (any(bad)) {
    if (isTRUE(getOption("LibeRtAD.warn_nonfinite_grad", TRUE))) {
      warning(
        "Non-finite gradient component(s) replaced with 0 (indices: ",
        paste(which(bad), collapse = ", "),
        "). Set options(LibeRtAD.warn_nonfinite_grad = FALSE) to silence.",
        call. = FALSE
      )
    }
    g[bad] <- 0
  }
  g
}

#' @keywords internal
.nm_sanitize_eta_mat <- function(eta_mat, omega = NULL) {
  if (is.null(eta_mat) || !is.matrix(eta_mat)) {
    return(eta_mat)
  }
  if (any(!is.finite(eta_mat))) {
    eta_mat[!is.finite(eta_mat)] <- 0
  }
  n_eta <- ncol(eta_mat)
  if (n_eta > 0L) {
    bound <- rep(10, n_eta)
    if (!is.null(omega) && length(omega) >= n_eta) {
      bound <- pmax(5, 5 * sqrt(pmax(omega[seq_len(n_eta)], 1e-8)))
    }
    for (j in seq_len(n_eta)) {
      eta_mat[, j] <- pmax(pmin(eta_mat[, j], bound[j]), -bound[j])
    }
  }
  eta_mat
}

#' @keywords internal
.nm_run_optim <- function(par, fn, gr = NULL, method = "L-BFGS-B",
                          lower = -Inf, upper = Inf, control = list()) {
  fn_safe <- function(x) .nm_finite_obj(fn(x))
  gr_safe <- if (!is.null(gr)) {
    function(x) .nm_finite_grad(gr(x))
  } else {
    NULL
  }
  args <- list(
    par = par,
    fn = fn_safe,
    method = method,
    control = control
  )
  if (!is.null(gr_safe)) {
    args$gr <- gr_safe
  }
  if (identical(method, "L-BFGS-B")) {
    args$lower <- lower
    args$upper <- upper
  }
  tryCatch(
    do.call(optim, args),
    error = function(e) {
      msg <- conditionMessage(e)
      if (!grepl("non-finite", msg, fixed = TRUE)) {
        stop(e)
      }
      args2 <- args
      args2$gr <- NULL
      args2$method <- "Nelder-Mead"
      args2$lower <- NULL
      args2$upper <- NULL
      do.call(optim, args2)
    }
  )
}

#' @keywords internal
.nm_opt_grad <- function(f, par, par_names = NULL, grad = "auto", backend = "cpp") {
  grad <- match.arg(grad, c("auto", "ad", "numeric"))
  if (is.null(par_names)) {
    par_names <- names(par)
  }
  if (is.null(par_names)) {
    par_names <- paste0("p", seq_along(par))
  }
  if (grad == "numeric") {
    return(.nm_num_grad(f, par))
  }
  if (is.function(f) && !is.null(formals(f)) && length(formals(f)) > 0L) {
    return(.nm_grad_population(f, par, par_names, grad, backend))
  }
  named <- stats::setNames(as.list(par), par_names)
  wrapper <- function(...) {
    args <- list(...)
    x <- unlist(args[par_names])
    f(x)
  }
  formals(wrapper) <- stats::setNames(
    rep(list(quote(expr = )), length(par_names)),
    par_names
  )
  .nm_grad_population(wrapper, par, par_names, grad, backend)
}

#' @keywords internal
.nm_optimize_par <- function(model,
                             data,
                             par0,
                             objective,
                             backend = "cpp",
                             grad = "auto",
                             pk_engine = "auto",
                             control = list(),
                             fixed = NULL,
                             ad_objective = NULL,
                             ad_presolve = NULL,
                             eta_mat = NULL,
                             include_omega_prior = TRUE,
                             cpp_pop_grad = NULL,
                             gr_fn = NULL) {
  fix_mask <- .nm_fix_mask(model)
  if (!is.null(fixed)) {
    fix_mask <- fix_mask | fixed
  }
  free <- which(!fix_mask)
  par_names <- .nm_par_labels(model)
  if (length(free) == 0L) {
    return(list(par = par0, value = objective(par0), convergence = 0L))
  }
  par0 <- .nm_apply_fix(model, par0)
  est_ctl <- .nm_split_est_control(control)
  use_cpp_pop_grad <- if (is.null(cpp_pop_grad)) {
    .nm_use_cpp_pop_grad(model, grad)
  } else {
    isTRUE(cpp_pop_grad) && .nm_cpp_capable(model)
  }
  if (!use_cpp_pop_grad) {
    .nm_set_cpp_pk_ad(model, backend)
  }
  on.exit({
    .nm_state$use_cpp_pk <- FALSE
    .nm_state$cpp_pk_ad_mode <- NA_character_
    .nm_state$optim_cache <- NULL
  }, add = TRUE)
  ad_obj <- if (!use_cpp_pop_grad) .nm_ad_objective_fn(ad_objective) else NULL
  pinned_obj <- function(par) .nm_finite_obj(objective(par))
  obj_fn <- if (!is.null(ad_obj)) {
    ad_obj
  } else {
    objective
  }
  if (use_cpp_pop_grad) {
    .nm_clear_pop_optim_cache()
    f_free <- function(x) {
      par <- par0
      par[free] <- x
      par <- .nm_apply_fix(model, par)
      p <- .nm_unpack(model, par)
      val <- .nm_pop_nll_cpp(
        model, data, p$theta, p$omega, p$sigma, eta_mat,
        include_omega_prior = include_omega_prior, cache_fwd = TRUE
      )
      .nm_finite_obj(val)
    }
    g_free <- function(x) {
      par <- par0
      par[free] <- x
      par <- .nm_apply_fix(model, par)
      p <- .nm_unpack(model, par)
      g <- .nm_pop_nll_grad_cpp(
        model, data, p$theta, p$omega, p$sigma, eta_mat,
        include_omega_prior = include_omega_prior, use_fwd_cache = TRUE
      )$gradient[free]
      .nm_finite_grad(g)
    }
  } else {
    use_ad_value <- !is.null(ad_obj) && is.null(ad_presolve)
    f_free <- function(x) {
      par <- par0
      par[free] <- x
      par <- .nm_apply_fix(model, par)
      if (!is.null(ad_presolve)) {
        ad_presolve(par)
      }
      if (use_ad_value) {
        at <- stats::setNames(as.list(par), par_names)
        return(.nm_finite_obj(.nm_ad_eval_cached(
          ad_obj, at, par_names, backend, need_grad = FALSE
        )))
      }
      pinned_obj(par)
    }
    use_ad <- .nm_grad_uses_ad(grad) && !is.null(ad_obj) && is.null(gr_fn)
    g_free <- if (!is.null(gr_fn)) {
      function(x) {
        par <- par0
        par[free] <- x
        par <- .nm_apply_fix(model, par)
        .nm_finite_grad(gr_fn(par))
      }
    } else if (use_ad) {
      function(x) {
        par <- par0
        par[free] <- x
        par <- .nm_apply_fix(model, par)
        if (!is.null(ad_presolve)) {
          ad_presolve(par)
        }
        at <- stats::setNames(as.list(par), par_names)
        grad <- .nm_ad_eval_cached(
          ad_obj, at, par_names, backend, need_grad = TRUE
        )
        .nm_finite_grad(grad[free])
      }
    } else {
      function(x) {
        par <- par0
        par[free] <- x
        par <- .nm_apply_fix(model, par)
        g <- .nm_opt_grad(pinned_obj, par, par_names, "numeric", backend)[free]
        .nm_finite_grad(g)
      }
    }
  }
  grad_names <- par_names[free]
  g_free <- .nm_wrap_grad_trace(
    g_free,
    est_ctl$print_grad_every,
    grad_names,
    prefix = "Population"
  )
  use_scale <- .nm_par_scale_enabled(control)
  par_scale <- if (use_scale) {
    .nm_par_scale_vector(model, par0)
  } else {
    NULL
  }
  if (use_scale && !is.null(par_scale)) {
    wrapped <- .nm_wrap_scaled_optim(model, par0, free, f_free, g_free, par_scale)
    opt <- .nm_run_optim(
      par = wrapped$par,
      fn = wrapped$fn,
      gr = wrapped$gr,
      method = "L-BFGS-B",
      lower = .nm_par_lower(model)[free] / par_scale[free],
      upper = .nm_par_upper(model)[free] / par_scale[free],
      control = c(list(factr = 1e7), est_ctl$optim_control)
    )
    par <- par0
    par[free] <- .nm_par_from_scaled(opt$par, par_scale[free])
  } else {
    opt <- .nm_run_optim(
      par = par0[free],
      fn = f_free,
      gr = g_free,
      method = "L-BFGS-B",
      lower = .nm_par_lower(model)[free],
      upper = .nm_par_upper(model)[free],
      control = c(list(factr = 1e7), est_ctl$optim_control)
    )
    par <- par0
    par[free] <- opt$par
  }
  par <- .nm_apply_fix(model, par)
  list(par = par, value = opt$value, convergence = opt$convergence, optim = opt)
}

#' @keywords internal
.nm_fit_eta_subject <- function(model, subj, theta, omega, sigma, eta0 = NULL,
                                backend = "cpp", grad = "auto", pk_engine = "auto") {
  n_eta <- .nm_n_eta(model)
  if (n_eta == 0L) {
    return(numeric())
  }
  pk_eff <- .nm_resolve_pk_engine(pk_engine, model, theta, omega, sigma, eta0)
  if (is.null(eta0)) {
    eta0 <- rep(0, n_eta)
  }
  eta_names <- .nm_eta_labels(model)
  obj_fn <- .nm_build_eta_objective(model, subj, theta, omega, sigma, pk_eff)
  # Per-subject eta fits use numeric gradients even when outer estimation uses AD.
  eta_grad <- if (.nm_grad_uses_ad(grad)) "numeric" else grad
  use_ad <- .nm_grad_uses_ad(eta_grad)
  .nm_set_cpp_pk_ad(model, backend)
  on.exit({
    .nm_state$use_cpp_pk <- FALSE
    .nm_state$cpp_pk_ad_mode <- NA_character_
    .nm_state$optim_cache <- NULL
  }, add = TRUE)
  fn <- function(eta) {
    names(eta) <- eta_names
    if (use_ad) {
      at <- stats::setNames(as.list(eta), eta_names)
      val <- .nm_ad_eval_cached(
        obj_fn, at, eta_names, backend, need_grad = FALSE
      )
      .nm_finite_obj(val)
    } else {
      .nm_finite_obj(do.call(obj_fn, as.list(eta)))
    }
  }
  gr <- if (use_ad) {
    function(eta) {
      names(eta) <- eta_names
      at <- stats::setNames(as.list(eta), eta_names)
      .nm_finite_grad(.nm_ad_eval_cached(obj_fn, at, eta_names, backend, need_grad = TRUE))
    }
  } else {
    function(eta) {
      .nm_finite_grad(.nm_grad_eta(obj_fn, eta, eta_names, eta_grad, backend))
    }
  }
  opt <- .nm_run_optim(par = eta0, fn = fn, gr = gr, method = "BFGS")
  opt$par
}

#' @keywords internal
.nm_ad_objective_fn <- function(ad_objective) {
  if (is.null(ad_objective)) {
    return(NULL)
  }
  if (is.list(ad_objective) && !is.null(ad_objective$fn)) {
    return(ad_objective$fn)
  }
  ad_objective
}

#' Fit post-hoc ETAs for all subjects (FOCE / FOCEI helper).
#' @keywords internal
.nm_fit_all_eta_modes <- function(model, dat, theta, omega, sigma, eta_mat = NULL,
                                  pk_engine = "auto", control = list()) {
  eta_init <- .nm_mceta_init(model, dat, theta, omega, sigma, eta_mat, control)
  eta_out <- .nm_fit_all_eta(
    model, dat, theta, omega, sigma, eta_init,
    backend = "cpp", grad = "numeric", pk_engine = pk_engine, control = control
  )
  if (is.matrix(eta_out) && nrow(eta_out) > 0L && all(abs(eta_out) < 1e-12)) {
    mceta <- control$mceta %||% getOption("LibeRation.mceta", "last")
    if (!identical(mceta, "zero")) {
      ctl1 <- control
      ctl1$n_cores <- 1L
      eta_out <- .nm_fit_all_eta(
        model, dat, theta, omega, sigma, NULL,
        backend = "cpp", grad = "numeric", pk_engine = pk_engine, control = ctl1
      )
    }
  }
  eta_out
}

#' @keywords internal
.nm_fit_all_eta <- function(model, dat, theta, omega, sigma, eta_mat = NULL,
                            backend = "cpp", grad = "auto", pk_engine = "auto",
                            control = list()) {
  ids <- .nm_subject_ids(dat)
  n_eta <- .nm_n_eta(model)
  if (n_eta == 0L) {
    return(matrix(nrow = 0, ncol = 0))
  }
  pk_eff <- .nm_resolve_pk_engine(pk_engine, model, theta, omega, sigma, eta_mat)
  warm <- !is.null(eta_mat) && is.matrix(eta_mat) && nrow(eta_mat) > 0L
  eta_maxit <- .nm_focei_eta_max_iter(control, warm = warm)
  if (.nm_cpp_capable(model) && pk_eff == "cpp") {
    eta_out <- .nm_fit_all_eta_cpp(
      model, dat, theta, omega, sigma, eta_mat, max_iter = eta_maxit
    )
    return(.nm_sanitize_eta_mat(eta_out))
  }
  eta_out <- matrix(0, length(ids), n_eta)
  n_cores <- .nm_resolve_n_cores(control)
  if (n_cores > 1L && length(ids) > 1L) {
    jobs <- lapply(seq_along(ids), function(j) {
      list(
        model = model,
        subj = .nm_subject_slice(dat, ids[j]),
        theta = theta,
        omega = omega,
        sigma = sigma,
        eta0 = if (is.null(eta_mat)) NULL else eta_mat[j, ],
        backend = backend,
        grad = grad,
        pk_engine = pk_engine
      )
    })
    fit_one <- function(job) {
      .nm_fit_eta_subject(
        job$model, job$subj, job$theta, job$omega, job$sigma, job$eta0,
        job$backend, job$grad, job$pk_engine
      )
    }
    res <- .nm_parallel_lapply(jobs, fit_one, n_cores = n_cores)
    for (j in seq_along(ids)) {
      eta_out[j, ] <- res[[j]]
    }
    return(eta_out)
  }
  for (j in seq_along(ids)) {
    subj <- .nm_subject_slice(dat, ids[j])
    eta0 <- if (is.null(eta_mat)) NULL else eta_mat[j, ]
    eta_out[j, ] <- .nm_fit_eta_subject(
      model, subj, theta, omega, sigma, eta0, backend, grad, pk_engine
    )
  }
  eta_out
}

#' @keywords internal
.nm_est_opts <- function(model, data, par0, backend, grad, pk_engine, engine,
                         control, ...) {
  list(
    model = model,
    data = data,
    par0 = par0,
    backend = match.arg(backend, c("cpp", "R")),
    grad = match.arg(grad, c("auto", "ad", "numeric")),
    pk_engine = match.arg(pk_engine, c("auto", "cpp", "R")),
    engine = match.arg(engine, c("auto", "R", "cpp")),
    control = control,
    dots = list(...)
  )
}

#' @keywords internal
.nm_use_cpp_engine <- function(engine, model, method) {
  engine <- match.arg(engine, c("auto", "R", "cpp"))
  if (engine == "R") {
    return(FALSE)
  }
  if (!.nm_cpp_capable(model)) {
    return(FALSE)
  }
  if (engine == "cpp") {
    return(TRUE)
  }
  method %in% c("SAEM", "LAPLACE", "BAYES")
}

#' First-order (FO) estimation
#'
#' Two-step NONMEM-style FO: (1) THETA and SIGMA at eta = 0; (2) fit subject
#' etas then OMEGA from the random-effects prior \eqn{\sum_i \eta_i^{\top} \Omega^{-1} \eta_i}.
#'
#' @keywords internal
.nm_est_fo <- function(model, data, par0, backend = "cpp", grad = "auto",
                     pk_engine = "auto", engine = "auto", control = list()) {
  pk_eff <- .nm_effective_pk_engine(pk_engine, grad)
  ginfo <- .nm_resolve_estimation_grad(model, grad)
  pop_grad <- ginfo$grad
  par_names <- .nm_par_labels(model)
  omega_fixed <- grepl("^OMEGA", par_names)
  dat <- .nm_prepare_data(data, model$INPUT, model)

  ad_obj1 <- if (.nm_grad_uses_ad(pop_grad)) {
    .nm_build_pop_objective(
      model, data, eta_mat = NULL, include_omega_prior = FALSE, pk_engine = pk_eff
    )$fn
  } else {
    NULL
  }
  objective1 <- function(par) {
    p <- .nm_unpack(model, par)
    .nm_nll_internal(
      model, data, p$theta, p$omega, p$sigma,
      eta = NULL, include_omega_prior = FALSE, pk_engine = pk_eff
    )
  }
  fit1 <- .nm_optimize_par(
    model, data, par0, objective1, backend, pop_grad, pk_eff, control,
    fixed = omega_fixed, ad_objective = ad_obj1, include_omega_prior = FALSE
  )

  n_eta <- .nm_n_eta(model)
  fit2 <- NULL
  eta_mat <- NULL
  par <- fit1$par
  if (n_eta > 0L) {
    p1 <- .nm_unpack(model, par)
    eta_mat <- .nm_fit_all_eta(
      model, dat, p1$theta, p1$omega, p1$sigma,
      backend = backend, grad = grad, pk_engine = pk_eff, control = control
    )
    theta_sigma_fixed <- grepl("^THETA|^SIGMA", par_names)
    fit2 <- .nm_optimize_fo_omega(
      model, par, eta_mat, fixed = theta_sigma_fixed, control = control
    )
    par <- fit2$par
    p2 <- .nm_unpack(model, par)
    eta_mat <- .nm_fit_all_eta(
      model, dat, p2$theta, p2$omega, p2$sigma, eta_mat,
      backend = backend, grad = grad, pk_engine = pk_eff, control = control
    )
  }

  p <- .nm_unpack(model, par)
  total_obj <- .nm_fo_report_objective(
    model, data, par, eta_mat, pk_engine = pk_eff
  )
  conv <- fit1$convergence
  if (!is.null(fit2)) {
    conv <- max(conv, fit2$convergence)
  }
  structure(
    list(
      method = "FO",
      par = par,
      theta = p$theta,
      omega = p$omega,
      sigma = p$sigma,
      eta = eta_mat,
      objective = total_obj,
      convergence = conv,
      grad = pop_grad,
      grad_requested = grad,
      grad_backend = .nm_report_grad_backend(model, pop_grad, backend),
      pk_engine = pk_eff,
      fo_estimation = list(step1 = fit1, step2 = fit2),
      optim = list(step1 = fit1$optim, step2 = if (!is.null(fit2)) fit2$optim else NULL)
    ),
    class = "nm_fit"
  )
}

#' First-order conditional estimation (FOCE)
#'
#' @keywords internal
.nm_est_foce <- function(model, data, par0, backend = "cpp", grad = "auto",
                         pk_engine = "auto", engine = "auto",
                         control = list(), max_outer = 20L, tol = 1e-4) {
  pk_eff <- .nm_effective_pk_engine(pk_engine, grad)
  ginfo <- .nm_resolve_estimation_grad(model, grad)
  est_grad <- ginfo$grad
  .nm_clear_laplace_optim_cache()
  .nm_clear_laplace_eta_modes()
  dat <- .nm_prepare_data(data, model$INPUT, model)
  par <- .nm_apply_fix(model, par0)
  eta_mat <- NULL
  outer <- list()
  optim_runs <- list()
  fit <- NULL
  pop_ad <- if (.nm_grad_uses_ad(est_grad)) {
    .nm_build_pop_objective(
      model, data, eta_mat, include_omega_prior = TRUE, pk_engine = pk_eff
    )
  } else {
    NULL
  }
  prev_val <- NULL
  for (iter in seq_len(max_outer)) {
    .nm_clear_pop_optim_cache()
    p <- .nm_unpack(model, par)
    eta_mat <- .nm_fit_all_eta_modes(
      model, dat, p$theta, p$omega, p$sigma, eta_mat,
      pk_eff, control
    )
    eta_mat <- .nm_sanitize_eta_mat(eta_mat)
    if (!is.null(pop_ad)) {
      pop_ad$ctx$eta_mat <- eta_mat
    }
    ad_obj <- if (!is.null(pop_ad)) pop_ad$fn else NULL
    objective <- function(par) {
      pp <- .nm_unpack(model, par)
      .nm_finite_obj(.nm_nll_internal(
        model, data, pp$theta, pp$omega, pp$sigma,
        eta = eta_mat, include_omega_prior = TRUE, pk_engine = pk_eff
      ))
    }
    fit <- .nm_optimize_par(
      model, data, par, objective, backend, est_grad, pk_eff, control,
      ad_objective = ad_obj, eta_mat = eta_mat, include_omega_prior = TRUE
    )
    optim_runs[[iter]] <- fit$optim
    par <- fit$par
    outer[[iter]] <- fit$value
    if (!is.null(prev_val) &&
        abs(prev_val - fit$value) < tol * (1 + abs(prev_val))) {
      break
    }
    prev_val <- fit$value
  }
  p <- .nm_unpack(model, par)
  eta_mat <- .nm_fit_all_eta_modes(
    model, dat, p$theta, p$omega, p$sigma, eta_mat,
    pk_eff, control
  )
  final_obj <- .nm_nll_internal(
    model, data, p$theta, p$omega, p$sigma,
    eta = eta_mat, include_omega_prior = TRUE, pk_engine = pk_eff
  )
  structure(
    list(
      method = "FOCE",
      par = par,
      theta = p$theta,
      omega = p$omega,
      sigma = p$sigma,
      eta = eta_mat,
      objective = final_obj,
      outer = outer,
      convergence = fit$convergence,
      grad = est_grad,
      grad_requested = grad,
      grad_backend = .nm_report_grad_backend(model, est_grad, backend),
      pk_engine = pk_eff,
      optim = optim_runs
    ),
    class = "nm_fit"
  )
}

#' Stochastic approximation EM (SAEM)
#'
#' @keywords internal
.nm_est_saem <- function(model, data, par0, backend = "cpp", grad = "auto",
                         pk_engine = "auto", engine = "auto",
                         control = list(),
                         n_iter = 100L, n_burn = 20L, n_mcmc = 1L,
                         sa_rate = 0.6, seed = 1L) {
  pk_eff <- .nm_effective_pk_engine(pk_engine, grad)
  ginfo <- .nm_resolve_estimation_grad(model, grad)
  est_grad <- ginfo$grad
  set.seed(seed)
  dat <- .nm_prepare_data(data, model$INPUT, model)
  ids <- .nm_subject_ids(dat)
  n_sub <- length(ids)
  n_eta <- .nm_n_eta(model)
  par <- .nm_apply_fix(model, par0)
  p <- .nm_unpack(model, par)
  eta_mat <- matrix(0, n_sub, n_eta)
  if (n_eta > 0L) {
    eta_mat <- .nm_fit_all_eta(
      model, dat, p$theta, p$omega, p$sigma, NULL,
      backend, "numeric", pk_eff, control
    )
  }
  theta_sa <- p$theta
  omega_sa <- p$omega
  sigma_sa <- p$sigma
  trace <- matrix(NA_real_, n_iter, length(par))
  optim_runs <- list()
  par_mstep <- par
  use_cpp <- .nm_use_cpp_engine(engine, model, "SAEM")
  use_cpp_mstep <- .nm_use_cpp_pop_grad(model, est_grad)
  if (use_cpp_mstep) {
    .nm_cpp_subjects_cached(model, data)
  }
  pop_ad <- if (.nm_grad_uses_ad(est_grad) && !use_cpp_mstep) {
    .nm_build_pop_objective(
      model, data, eta_mat, include_omega_prior = TRUE, pk_engine = pk_eff
    )
  } else {
    NULL
  }
  for (k in seq_len(n_iter)) {
    p <- .nm_unpack(model, par)
    if (n_eta > 0L) {
      if (use_cpp) {
        eta_mat <- .nm_saem_mh_cpp(
          model, data, eta_mat, p$theta, p$omega, p$sigma,
          n_mcmc = n_mcmc
        )
      } else {
        eta_mat <- .nm_saem_mh_r(
          model, dat, ids, eta_mat, p$theta, p$omega, p$sigma,
          n_mcmc = n_mcmc, pk_engine = pk_eff, control = control
        )
      }
      eta_mat <- .nm_sanitize_eta_mat(eta_mat, p$omega)
    }
    if (!is.null(pop_ad)) {
      pop_ad$ctx$eta_mat <- eta_mat
    }
    ad_obj <- if (!is.null(pop_ad)) pop_ad$fn else NULL
    if (use_cpp_mstep) {
      objective <- function(par) {
        pp <- .nm_unpack(model, par)
        .nm_finite_obj(.nm_pop_nll_cpp(
          model, data, pp$theta, pp$omega, pp$sigma, eta_mat,
          include_omega_prior = TRUE
        ))
      }
    } else {
      objective <- function(par) {
        pp <- .nm_unpack(model, par)
        .nm_finite_obj(.nm_nll_internal(
          model, data, pp$theta, pp$omega, pp$sigma,
          eta = eta_mat, include_omega_prior = TRUE, pk_engine = pk_eff
        ))
      }
    }
    fit <- .nm_optimize_par(
      model, data, par_mstep, objective, backend, est_grad, pk_eff,
      control = .nm_saem_mstep_control(control, k, n_burn),
      ad_objective = ad_obj, eta_mat = eta_mat, include_omega_prior = TRUE
    )
    optim_runs[[k]] <- fit$optim
    par_mstep <- fit$par
    par_new <- fit$par
    pp_new <- .nm_unpack(model, par_new)
    gamma_k <- if (k <= n_burn) 1 else (k - n_burn)^(-sa_rate)
    gamma_k <- min(1, gamma_k)
    theta_sa <- (1 - gamma_k) * theta_sa + gamma_k * pp_new$theta
    sigma_sa <- (1 - gamma_k) * sigma_sa + gamma_k * pp_new$sigma
    if (n_eta > 0L) {
      omega_sa <- (1 - gamma_k) * omega_sa + gamma_k * pp_new$omega
    }
    par <- .nm_pack(model, theta_sa, omega_sa, sigma_sa)
    par <- .nm_apply_fix(model, par)
    trace[k, ] <- par
  }
  p <- .nm_unpack(model, par)
  structure(
    list(
      method = "SAEM",
      par = par,
      theta = p$theta,
      omega = p$omega,
      sigma = p$sigma,
      eta = eta_mat,
      objective = objective(par),
      trace = trace,
      convergence = 0L,
      grad = est_grad,
      grad_requested = grad,
      grad_backend = if (use_cpp_mstep) {
        "cpp"
      } else if (.nm_grad_uses_ad(est_grad)) {
        backend
      } else {
        NULL
      },
      pk_engine = pk_eff,
      engine = if (use_cpp) "cpp" else "R",
      optim = optim_runs
    ),
    class = "nm_fit"
  )
}

#' Laplace approximation with Gaussian quadrature (NONMEM-style)
#'
#' @keywords internal
.nm_subject_laplace_nll <- function(model, subj, theta, omega, sigma,
                                    gh, pk_engine = "auto") {
  val <- .nm_subject_laplace_nll_internal(
    model, subj, theta, omega, sigma, gh, pk_engine
  )
  if (.nm_any_ad(val)) {
    return(.ad_scalar_value(val))
  }
  if (!is.finite(val)) .Machine$double.xmax else val
}

#' @keywords internal
.nm_est_laplace <- function(model, data, par0, backend = "cpp", grad = "auto",
                            pk_engine = "auto", engine = "auto",
                            control = list(), n_quad = 5L) {
  pk_eff <- .nm_effective_pk_engine(pk_engine, grad)
  ginfo <- .nm_resolve_estimation_grad(model, grad)
  est_grad <- ginfo$grad
  .nm_clear_laplace_optim_cache()
  .nm_clear_laplace_eta_modes()
  dat <- .nm_prepare_data(data, model$INPUT, model)
  n_quad_eff <- .nm_effective_n_quad(model, n_quad)
  gh <- .nm_gh_nodes(n_quad_eff)
  ids <- .nm_subject_ids(dat)
  lap_grad <- .nm_resolve_laplace_grad(model, grad, gh, length(ids))
  if (identical(lap_grad, "ad") && !.nm_ad_pk_supported(model)) {
    .nm_stop(
      "Tape AD Laplace gradients are not supported for ADVAN ", model$ADVAN,
      " TRANS ", model$TRANS, ". Use grad = \"cpp\" or grad = \"auto\"."
    )
  }
  use_cpp_grad <- lap_grad == "cpp"
  use_ad_outer <- lap_grad == "ad"
  use_cpp <- .nm_use_cpp_engine(engine, model, "LAPLACE") || use_cpp_grad
  use_cpp <- use_cpp && .nm_cpp_capable(model)
  est_ctl <- .nm_split_est_control(control)
  ad_obj <- if (use_ad_outer) {
    .nm_build_laplace_objective(model, data, gh, pk_engine = pk_eff)
  } else {
    NULL
  }
  laplace_eval <- function(theta, omega, sigma) {
    if (use_cpp) {
      return(.nm_laplace_nll_cpp(model, data, theta, omega, sigma, gh))
    }
    val <- 0
    for (id in ids) {
      subj <- .nm_subject_slice(dat, id)
      val <- val + .nm_subject_laplace_nll(
        model, subj, theta, omega, sigma, gh, pk_eff
      )
    }
    val
  }
  objective <- function(par) {
    p <- .nm_unpack(model, par)
    laplace_eval(p$theta, p$omega, p$sigma)
  }
  if (use_cpp_grad) {
    .nm_clear_laplace_optim_cache()
    fix_mask <- .nm_fix_mask(model)
    free <- which(!fix_mask)
    par0 <- .nm_apply_fix(model, par0)
    par_names <- .nm_par_labels(model)
    f_free <- function(x) {
      par <- par0
      par[free] <- x
      par <- .nm_apply_fix(model, par)
      p <- .nm_unpack(model, par)
      val <- .nm_laplace_nll_cpp(
        model, data, p$theta, p$omega, p$sigma, gh,
        control = control, cache_fwd = TRUE
      )
      if (!is.finite(val)) .nm_finite_obj(val) else val
    }
    g_free <- function(x) {
      par <- par0
      par[free] <- x
      par <- .nm_apply_fix(model, par)
      p <- .nm_unpack(model, par)
      .nm_finite_grad(.nm_laplace_nll_grad_cpp(
        model, data, p$theta, p$omega, p$sigma, gh,
        control = control, use_fwd_cache = TRUE
      )$gradient[free])
    }
    g_free <- .nm_wrap_grad_trace(
      g_free, est_ctl$print_grad_every, par_names[free], prefix = "Population"
    )
    run_laplace_opt <- function() {
      .nm_run_optim(
        par = par0[free],
        fn = f_free,
        gr = g_free,
        method = "L-BFGS-B",
        lower = .nm_par_lower(model)[free],
        upper = .nm_par_upper(model)[free],
        control = c(list(factr = 1e7), est_ctl$optim_control)
      )
    }
    opt <- run_laplace_opt()
    par_start <- par0[free]
    if (isTRUE(getOption("LibeRtAD.laplace_opt_reset", TRUE)) &&
        max(abs(opt$par - par_start)) < 1e-6 * (1 + max(abs(par_start)))) {
      .nm_clear_laplace_eta_modes()
      .nm_clear_laplace_optim_cache()
      opt <- run_laplace_opt()
    }
    if (opt$convergence != 0L) {
      .nm_clear_laplace_eta_modes()
      .nm_clear_laplace_optim_cache()
      opt <- run_laplace_opt()
    }
    par <- par0
    par[free] <- opt$par
    par <- .nm_apply_fix(model, par)
    fit <- list(par = par, value = opt$value, convergence = opt$convergence, optim = opt)
  } else {
    fit <- .nm_optimize_par(
      model, data, par0,
      function(par) {
        val <- objective(par)
        if (!is.finite(val)) .nm_finite_obj(val) else val
      },
      backend, lap_grad, pk_eff, control,
      ad_objective = ad_obj,
      cpp_pop_grad = FALSE
    )
  }
  p <- .nm_unpack(model, fit$par)
  eta_mat <- .nm_laplace_final_eta_modes(
    model, data, p$theta, p$omega, p$sigma, gh, pk_engine = pk_eff
  )
  structure(
    list(
      method = "LAPLACE",
      par = fit$par,
      theta = p$theta,
      omega = p$omega,
      sigma = p$sigma,
      eta = eta_mat,
      objective = fit$value,
      n_quad = n_quad_eff,
      n_quad_requested = n_quad,
      convergence = fit$convergence,
      grad = lap_grad,
      grad_requested = grad,
      grad_effective = lap_grad,
      grad_backend = if (use_cpp_grad) "cpp" else if (use_ad_outer) backend else NULL,
      pk_engine = pk_eff,
      engine = if (use_cpp) "cpp" else "R",
      laplace_mode_centered = .nm_laplace_mode_centered(),
      optim = fit$optim
    ),
    class = "nm_fit"
  )
}
