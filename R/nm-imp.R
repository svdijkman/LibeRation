#' Importance sampling (IMP) estimation
#'
#' Uses Laplace eta-modes and Gaussian importance weights around each mode.
#'
#' @keywords internal
.nm_est_imp <- function(model, data, par0, backend = "cpp", grad = "auto",
                        pk_engine = "auto", engine = "auto", control = list(),
                        n_imp = 50L, n_quad = 5L, ...) {
  .nm_sync_lik_config(model)
  pk_eff <- .nm_effective_pk_engine(pk_engine, grad)
  ginfo <- .nm_resolve_estimation_grad(model, grad)
  est_grad <- ginfo$grad
  dat <- .nm_prepare_data(data, model$INPUT, model)
  ids <- .nm_subject_ids(dat)
  n_eta <- .nm_n_eta(model)
  par <- .nm_apply_fix(model, par0)
  p <- .nm_unpack(model, par)
  eta_mat <- matrix(0, length(ids), n_eta)
  if (n_eta > 0L) {
    eta_mat <- .nm_fit_all_eta(
      model, dat, p$theta, p$omega, p$sigma, NULL,
      backend, "numeric", pk_eff, control
    )
  }
  gh <- .nm_gh_nodes(n_quad)
  imp_diag <- new.env(parent = emptyenv())
  imp_diag$min_ess_frac <- Inf
  imp_diag$fallback <- FALSE
  imp_seed <- .nm_imp_seed(control)
  objective <- function(par) {
    pp <- .nm_unpack(model, par)
    .nm_imp_nll(
      model, data, pp$theta, pp$omega, pp$sigma, eta_mat,
      gh, n_imp = n_imp, pk_engine = pk_eff,
      diag_env = imp_diag, seed = imp_seed
    )
  }
  # IMP objective differs from standard population NLL; tape AD not implemented yet.
  ad_obj <- NULL
  if (.nm_grad_uses_ad(est_grad)) {
    est_grad <- "numeric"
  }
  gb <- if (.nm_use_cpp_pop_grad(model, est_grad)) "cpp" else backend
  .nm_est_progress_phase("IMP", "start", list(n_imp = n_imp, n_quad = n_quad))
  control_opt <- control
  control_opt$imp_seed <- NULL
  fit <- .nm_optimize_par(
    model, data, par, objective, backend, est_grad, pk_eff,
    control = control_opt, ad_objective = ad_obj, eta_mat = eta_mat,
    include_omega_prior = TRUE
  )
  .nm_imp_emit_diagnostics(imp_diag, n_imp)
  .nm_est_progress_phase("IMP", "optimization complete", list(objective = fit$value))
  p <- .nm_unpack(model, fit$par)
  structure(
    list(
      method = "IMP",
      par = fit$par,
      theta = p$theta,
      omega = p$omega,
      sigma = p$sigma,
      eta = eta_mat,
      objective = fit$value,
      convergence = fit$convergence,
      grad = est_grad,
      grad_requested = grad,
      grad_backend = gb,
      pk_engine = pk_eff,
      n_imp = n_imp,
      n_quad = n_quad,
      optim = fit$optim,
      engine_detail = .nm_engine_detail(
        "IMP", est_grad, grad, gb,
        pk_eff, "R", .nm_lik_config(model),
        mstep = "optim", eta_fit = "mode+is"
      )
    ),
    class = "nm_fit"
  )
}

#' Resolve the IMP common-random-number seed.
#' @keywords internal
.nm_imp_seed <- function(control) {
  s <- control$imp_seed
  if (is.null(s)) s <- control$seed
  if (is.null(s)) 20240624L else as.integer(s)[1L]
}

#' Emit one-time IMP diagnostics (proposal downgrade / low ESS).
#' @keywords internal
.nm_imp_emit_diagnostics <- function(diag_env, n_imp) {
  if (isTRUE(diag_env$fallback)) {
    warning(
      "IMP: the conditional (posterior) covariance was not positive-definite ",
      "for at least one subject; the importance proposal fell back to the ",
      "prior Omega covariance there. Estimates remain valid but sampling may ",
      "be less efficient.",
      call. = FALSE
    )
  }
  if (is.finite(diag_env$min_ess_frac) && diag_env$min_ess_frac < 0.1) {
    warning(
      sprintf(
        paste0(
          "IMP: low importance-sampling effective sample size (min ESS = ",
          "%.1f%% of %d draws). Consider increasing n_imp for a more reliable ",
          "objective."
        ),
        100 * diag_env$min_ess_frac, n_imp
      ),
      call. = FALSE
    )
  }
}

#' Conditional (posterior) covariance of eta at the mode for the IMP proposal.
#'
#' Approximates the Laplace posterior covariance as the inverse Hessian of the
#' joint negative log density g(eta) = 0.5 * subject_(-2LL) at the mode. Falls
#' back to the prior Omega covariance when the Hessian is not positive-definite.
#' @keywords internal
.nm_imp_proposal_cov <- function(model, sub, theta, omega, sigma, eta0,
                                 pk_engine, n_eta) {
  g <- function(eta) {
    0.5 * .nm_subject_nll_internal(
      model, sub, theta, omega, sigma, eta,
      include_omega_prior = TRUE, pk_engine = pk_engine
    )
  }
  cov <- NULL
  h <- tryCatch(.nm_imp_num_hessian(g, eta0), error = function(e) NULL)
  if (!is.null(h) && all(is.finite(h))) {
    h <- (h + t(h)) / 2
    cov <- tryCatch(solve(h), error = function(e) NULL)
  }
  ok <- !is.null(cov) && all(is.finite(cov))
  ch <- NULL
  if (ok) {
    cov <- (cov + t(cov)) / 2
    ch <- tryCatch(chol(cov), error = function(e) NULL)
    ok <- !is.null(ch)
  }
  if (!ok) {
    v <- pmax(omega[seq_len(n_eta)], 1e-8)
    cov <- diag(v, n_eta)
    ch <- chol(cov)
    return(list(chol_lower = t(ch), logdet = sum(log(v)), fallback = TRUE))
  }
  list(chol_lower = t(ch), logdet = 2 * sum(log(diag(ch))), fallback = FALSE)
}

#' Central-difference Hessian for small eta dimension.
#' @keywords internal
.nm_imp_num_hessian <- function(fn, x, h = 1e-4) {
  n <- length(x)
  hess <- matrix(0, n, n)
  f0 <- fn(x)
  if (!is.finite(f0)) stop("non-finite objective at mode")
  for (i in seq_len(n)) {
    for (k in i:n) {
      xi <- x
      xi[i] <- xi[i] + h
      xi[k] <- xi[k] + h
      fpp <- fn(xi)
      xi <- x
      xi[i] <- xi[i] + h
      xi[k] <- xi[k] - h
      fpm <- fn(xi)
      xi <- x
      xi[i] <- xi[i] - h
      xi[k] <- xi[k] + h
      fmp <- fn(xi)
      xi <- x
      xi[i] <- xi[i] - h
      xi[k] <- xi[k] - h
      fmm <- fn(xi)
      val <- (fpp - fpm - fmp + fmm) / (4 * h * h)
      hess[i, k] <- val
      hess[k, i] <- val
    }
  }
  hess
}

#' NONMEM-style importance-sampling (IMP) marginal objective.
#'
#' For each subject the marginal likelihood L_i = int p(y_i | eta) p(eta) d eta
#' is estimated by importance sampling with a Gaussian proposal centred at the
#' conditional mode with the conditional (posterior) covariance:
#'   eta_m ~ q = N(mode_i, C_i),  w_m = p(y_i, eta_m) / q(eta_m),
#'   -2 log L_i ~= -2 [ logsumexp_m log(w_m) - log(M) ].
#' Common random numbers (a fixed standardized draw set) are reused across
#' objective evaluations so the objective is a smooth function of the
#' parameters, which stabilises the finite-difference gradient.
#' @keywords internal
.nm_imp_nll <- function(model, data, theta, omega, sigma, eta_modes, gh,
                        n_imp = 50L, pk_engine = "auto",
                        diag_env = NULL, seed = 20240624L) {
  dat <- .nm_prepare_data(data, model$INPUT, model)
  ids <- .nm_subject_ids(dat)
  n_sub <- length(ids)
  n_eta <- .nm_n_eta(model)
  if (n_eta == 0L) {
    total <- 0
    for (j in seq_len(n_sub)) {
      sub <- .nm_subject_slice(dat, ids[j])
      total <- total + .nm_subject_nll_internal(
        model, sub, theta, omega, sigma, numeric(0),
        include_omega_prior = TRUE, pk_engine = pk_engine
      )
    }
    return(total)
  }
  # Common random numbers: draw a fixed standardized set, restore RNG after.
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  }
  set.seed(seed)
  z_all <- array(rnorm(n_sub * n_imp * n_eta), dim = c(n_sub, n_imp, n_eta))
  log2pi <- log(2 * pi)
  total <- 0
  ess_frac_min <- Inf
  for (j in seq_len(n_sub)) {
    sub <- .nm_subject_slice(dat, ids[j])
    eta0 <- if (is.matrix(eta_modes)) eta_modes[j, ] else rep(0, n_eta)
    prop <- .nm_imp_proposal_cov(
      model, sub, theta, omega, sigma, eta0, pk_engine, n_eta
    )
    if (isTRUE(prop$fallback) && !is.null(diag_env)) {
      diag_env$fallback <- TRUE
    }
    L <- prop$chol_lower
    logdet <- prop$logdet
    logw <- numeric(n_imp)
    for (m in seq_len(n_imp)) {
      z <- z_all[j, m, ]
      eta_m <- eta0 + as.numeric(L %*% z)
      nll_m <- .nm_subject_nll_internal(
        model, sub, theta, omega, sigma, eta_m,
        include_omega_prior = TRUE, pk_engine = pk_engine
      )
      log_joint <- -0.5 * nll_m
      log_q <- -0.5 * (n_eta * log2pi + logdet + sum(z * z))
      logw[m] <- log_joint - log_q
    }
    mx <- max(logw)
    if (!is.finite(mx)) {
      total <- total + .Machine$double.xmax / n_sub
      next
    }
    w <- exp(logw - mx)
    sw <- sum(w)
    lse <- mx + log(sw)
    total <- total - 2 * (lse - log(n_imp))
    ess <- (sw * sw) / sum(w * w)
    ess_frac_min <- min(ess_frac_min, ess / n_imp)
  }
  if (!is.null(diag_env) && is.finite(ess_frac_min)) {
    diag_env$min_ess_frac <- min(diag_env$min_ess_frac, ess_frac_min)
  }
  total
}
