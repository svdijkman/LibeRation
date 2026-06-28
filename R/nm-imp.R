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
  objective <- function(par) {
    pp <- .nm_unpack(model, par)
    .nm_imp_nll(
      model, data, pp$theta, pp$omega, pp$sigma, eta_mat,
      gh, n_imp = n_imp, pk_engine = pk_eff
    )
  }
  # IMP objective differs from standard population NLL; tape AD not implemented yet.
  ad_obj <- NULL
  if (.nm_grad_uses_ad(est_grad)) {
    est_grad <- "numeric"
  }
  gb <- if (.nm_use_cpp_pop_grad(model, est_grad)) "cpp" else backend
  fit <- .nm_optimize_par(
    model, data, par, objective, backend, est_grad, pk_eff,
    control = control, ad_objective = ad_obj, eta_mat = eta_mat,
    include_omega_prior = TRUE
  )
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

#' @keywords internal
.nm_imp_nll <- function(model, data, theta, omega, sigma, eta_modes, gh,
                        n_imp = 50L, pk_engine = "auto") {
  dat <- .nm_prepare_data(data, model$INPUT, model)
  ids <- .nm_subject_ids(dat)
  n_sub <- length(ids)
  total <- 0
  for (j in seq_len(n_sub)) {
    sub <- .nm_subject_slice(dat, ids[j])
    eta0 <- if (is.matrix(eta_modes)) eta_modes[j, ] else rep(0, .nm_n_eta(model))
    draws <- matrix(NA_real_, n_imp, length(eta0))
    for (m in seq_len(n_imp)) {
      draws[m, ] <- eta0 + rnorm(length(eta0), sd = sqrt(pmax(omega, 1e-8)) / 2)
    }
    ll <- vapply(seq_len(n_imp), function(m) {
      -0.5 * .nm_subject_nll_internal(
        model, sub, theta, omega, sigma, draws[m, ],
        include_omega_prior = TRUE, pk_engine = pk_engine
      )
    }, numeric(1))
    total <- total - 2 * .nm_logsumexp(ll - max(ll)) - 2 * max(ll) + 2 * log(n_imp)
  }
  total
}
