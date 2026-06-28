#' @keywords internal
.nm_as_ad_constants <- function(x, prefix = "c") {
  if (length(x) == 0L || .nm_any_ad(x)) {
    return(x)
  }
  if (is.list(x)) {
    return(x)
  }
  lapply(seq_along(x), function(i) {
    newConstant(name = paste0(prefix, i), value = x[i])
  })
}

#' @keywords internal
.nm_zero_like <- function(...) {
  if (.nm_any_ad(...)) {
    return(newConstant(name = "nm_zero", value = 0))
  }
  0
}

#' @keywords internal
.nm_subject_nll_internal <- function(model,
                                     subj,
                                     theta,
                                     omega,
                                     sigma,
                                     eta,
                                     include_omega_prior = TRUE,
                                     pk_engine = "auto") {
  if (.nm_any_ad(theta, omega, sigma) && length(eta) > 0L && !.nm_any_ad(eta)) {
    eta <- .nm_as_ad_constants(eta, "eta")
  }
  pk_engine <- .nm_resolve_pk_engine(pk_engine, model, theta, omega, sigma, eta)
  pred <- .nm_subject_ipred(model, subj, theta, omega, eta, pk_engine = pk_engine)
  dv <- pred$subj_ev$DV[pred$obs_idx]
  f <- pred$F
  nll <- .nm_zero_like(theta, omega, sigma, eta, f, dv)
  cfg <- .nm_lik_config(model)
  dvid <- if ("DVID" %in% names(pred$subj_ev)) pred$subj_ev$DVID[pred$obs_idx] else NULL
  nll <- .ad_add(nll, .nm_residual_nll(
    dv, f, sigma,
    error = cfg$error,
    dvid = dvid,
    ar1_rho = cfg$ar1_rho %||% 0,
    sigma_corr = cfg$sigma_corr %||% "indep"
  ))
  if (include_omega_prior) {
    nll <- .ad_add(nll, .nm_omega_prior(eta, omega))
  }
  nll
}

#' @keywords internal
.nm_nll_internal <- function(model,
                             data,
                             theta,
                             omega,
                             sigma,
                             eta = NULL,
                             include_omega_prior = TRUE,
                             pk_engine = "auto") {
  pk_engine <- .nm_resolve_pk_engine(pk_engine, model, theta, omega, sigma, eta)
  if (pk_engine == "cpp" && !.nm_any_ad(theta, omega, sigma, eta)) {
    return(.nm_nll_cpp(
      model, data, theta, omega, sigma, eta, include_omega_prior
    ))
  }
  dat <- .nm_prepare_data(data, model$INPUT, model)
  ids <- .nm_subject_ids(dat)
  n_eta <- .nm_n_eta(model)
  if (is.null(eta)) {
    eta_mat <- matrix(0, length(ids), n_eta)
  } else {
    if (is.vector(eta)) {
      eta <- matrix(eta, nrow = 1L)
    }
    eta_mat <- eta
    if (nrow(eta_mat) != length(ids)) {
      .nm_stop("eta must have one row per subject.")
    }
    if (ncol(eta_mat) != n_eta) {
      .nm_stop("eta must have ", n_eta, " columns.")
    }
  }
  nll <- .nm_zero_like(theta, omega, sigma, eta_mat)
  for (j in seq_along(ids)) {
    id <- ids[j]
    subj <- .nm_subject_slice(dat, id)
    nll <- .ad_add(nll, .nm_subject_nll_internal(
      model, subj, theta, omega, sigma, eta_mat[j, ],
      include_omega_prior = include_omega_prior,
      pk_engine = pk_engine
    ))
  }
  nll
}

#' @keywords internal
.nm_laplace_eta_at_q <- function(mode, omega, z_q, mode_centered = TRUE) {
  if (!isTRUE(mode_centered)) {
    return(as.numeric(z_q))
  }
  as.numeric(mode + sqrt(pmax(omega, 1e-8)) * z_q)
}

#' @keywords internal
.nm_subject_laplace_nll_internal <- function(model,
                                           subj,
                                           theta,
                                           omega,
                                           sigma,
                                           gh,
                                           pk_engine = "R",
                                           mode_centered = .nm_laplace_mode_centered()) {
  n_eta <- .nm_n_eta(model)
  if (n_eta == 0L) {
    return(.nm_subject_nll_internal(
      model, subj, theta, omega, sigma, numeric(),
      include_omega_prior = TRUE, pk_engine = pk_engine
    ))
  }
  grid <- .nm_product_grid(gh$nodes, gh$weights, n_eta)
  log_w <- log(grid$weights)
  use_mode <- isTRUE(mode_centered) && !.nm_any_ad(theta, omega, sigma)
  eta_mode <- if (use_mode) {
    .nm_fit_eta_subject(
      model, subj, theta, omega, sigma, eta0 = rep(0, n_eta),
      backend = "cpp", grad = "numeric", pk_engine = pk_engine
    )
  } else {
    rep(0, n_eta)
  }
  log_terms <- vector("list", nrow(grid$nodes))
  for (q in seq_len(nrow(grid$nodes))) {
    eta_q <- .nm_laplace_eta_at_q(
      eta_mode, omega, grid$nodes[q, ], use_mode
    )
    nll_q <- .nm_subject_nll_internal(
      model, subj, theta, omega, sigma, eta_q,
      include_omega_prior = TRUE, pk_engine = pk_engine
    )
    log_terms[[q]] <- .ad_add(
      newConstant(name = paste0("logw_", q), value = log_w[q]),
      .ad_mul(-0.5, nll_q)
    )
  }
  .ad_mul(-2, .ad_logsumexp_scalars(log_terms))
}

#' @keywords internal
.nm_laplace_nll_internal <- function(model,
                                      data,
                                      theta,
                                      omega,
                                      sigma,
                                      gh,
                                      pk_engine = "R") {
  pk_engine <- .nm_resolve_pk_engine(pk_engine, model, theta, omega, sigma, NULL)
  dat <- .nm_prepare_data(data, model$INPUT, model)
  ids <- .nm_subject_ids(dat)
  nll <- .nm_zero_like(theta, omega, sigma)
  for (id in ids) {
    subj <- .nm_subject_slice(dat, id)
    nll <- .ad_add(nll, .nm_subject_laplace_nll_internal(
      model, subj, theta, omega, sigma, gh, pk_engine
    ))
  }
  nll
}

#' Negative log-likelihood for a NONMEM-style model
#'
#' @inheritParams nm_est
#' @param pk_engine PK solver: \code{"auto"}, \code{"cpp"}, or \code{"R"}.
#' @examples
#' sim <- nm_synthetic_theo(n_sub = 2L)
#' nm_nll(sim$model, sim$data, sim$model$THETAS$Value,
#'        sim$model$OMEGAS$Value, sim$model$SIGMAS$Value)
#' @export
nm_nll <- function(model,
                   data,
                   theta,
                   omega,
                   sigma,
                   eta = NULL,
                   include_omega_prior = TRUE,
                   pk_engine = c("auto", "cpp", "R")) {
  pk_engine <- match.arg(pk_engine)
  .nm_nll_internal(
    model, data, theta, omega, sigma, eta,
    include_omega_prior = include_omega_prior,
    pk_engine = pk_engine
  )
}

#' @keywords internal
.nm_conditional_nll <- function(model, subj, theta, omega, sigma, eta,
                              pk_engine = "auto") {
  .nm_subject_nll_internal(
    model, subj, theta, omega, sigma, eta,
    include_omega_prior = TRUE,
    pk_engine = pk_engine
  )
}

#' @keywords internal
.nm_subject_nll <- function(model, subj, theta, omega, sigma, eta,
                            include_omega_prior = TRUE,
                            pk_engine = "auto") {
  .nm_subject_nll_internal(
    model, subj, theta, omega, sigma, eta,
    include_omega_prior = include_omega_prior,
    pk_engine = pk_engine
  )
}
