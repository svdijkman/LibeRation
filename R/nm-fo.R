#' First-order (FO) likelihood helpers (NONMEM-style)
#'
#' FO uses eta = 0 in the data term and estimates OMEGA from the first-order
#' expansion: \eqn{\log|\Omega + G_i^{\top} R_i^{-1} G_i|} per subject, where \eqn{G_i} is
#' d f / d eta evaluated at eta = 0.
#'
#' @keywords internal
.nm_fo_eta_eps <- function() {
  1e-4
}

#' FO omega step: \eqn{\sum_i \eta_i^{\top} \Omega^{-1} \eta_i + \log|\Omega|} with fixed etas
#'
#' @keywords internal
.nm_fo_report_objective <- function(model, data, par, eta_mat, pk_engine = "auto") {
  p <- .nm_unpack(model, par)
  obj1 <- .nm_nll_internal(
    model, data, p$theta, p$omega, p$sigma,
    eta = NULL, include_omega_prior = FALSE, pk_engine = pk_engine
  )
  obj2 <- if (!is.null(eta_mat) && is.matrix(eta_mat) && nrow(eta_mat) > 0L) {
    .nm_fo_omega_prior_nll_scalar(eta_mat, p$omega)
  } else {
    0
  }
  obj1 + obj2
}

#' @keywords internal
.nm_fo_omega_prior_nll_scalar <- function(eta_mat, omega) {
  if (length(omega) == 0L || nrow(eta_mat) == 0L) {
    return(0)
  }
  n_sub <- nrow(eta_mat)
  total <- 0
  for (j in seq_along(omega)) {
    om <- max(omega[j], 1e-15)
    total <- total + sum(eta_mat[, j]^2) / om + n_sub * log(om)
  }
  total
}

#' @keywords internal
.nm_fo_omega_prior_grad <- function(eta_mat, omega) {
  if (length(omega) == 0L || nrow(eta_mat) == 0L) {
    return(numeric(0))
  }
  n_sub <- nrow(eta_mat)
  vapply(seq_along(omega), function(j) {
    om <- max(omega[j], 1e-15)
    -sum(eta_mat[, j]^2) / om^2 + n_sub / om
  }, numeric(1))
}

#' @keywords internal
.nm_optimize_fo_omega <- function(model, par0, eta_mat, fixed = NULL, control = list()) {
  fix_mask <- .nm_fix_mask(model)
  if (!is.null(fixed)) {
    fix_mask <- fix_mask | fixed
  }
  free <- which(!fix_mask)
  par0 <- .nm_apply_fix(model, par0)
  if (length(free) == 0L) {
    return(list(
      par = par0,
      value = .nm_fo_omega_prior_nll_scalar(eta_mat, .nm_unpack(model, par0)$omega),
      convergence = 0L
    ))
  }
  est_ctl <- .nm_split_est_control(control)
  par_names <- .nm_par_labels(model)
  omega_idx <- grepl("^OMEGA", par_names)
  omega_free <- which(omega_idx & !fix_mask)
  if (length(omega_free) == 0L) {
    p <- .nm_unpack(model, par0)
    return(list(
      par = par0,
      value = .nm_fo_omega_prior_nll_scalar(eta_mat, p$omega),
      convergence = 0L
    ))
  }
  f_free <- function(x) {
    par <- par0
    par[omega_free] <- x
    par <- .nm_apply_fix(model, par)
    p <- .nm_unpack(model, par)
    val <- .nm_fo_omega_prior_nll_scalar(eta_mat, p$omega)
    if (!is.finite(val)) .Machine$double.xmax else val
  }
  g_free <- function(x) {
    par <- par0
    par[omega_free] <- x
    par <- .nm_apply_fix(model, par)
    p <- .nm_unpack(model, par)
    g_omega <- .nm_fo_omega_prior_grad(eta_mat, p$omega)
    names(g_omega) <- paste0("OMEGA", model$OMEGAS$OMEGA)
    unname(g_omega[par_names[omega_free]])
  }
  opt <- .nm_run_optim(
    par = par0[omega_free],
    fn = f_free,
    gr = g_free,
    method = "L-BFGS-B",
    lower = .nm_par_lower(model)[omega_free],
    upper = .nm_par_upper(model)[omega_free],
    control = c(list(factr = 1e7), est_ctl$optim_control)
  )
  par <- par0
  par[omega_free] <- opt$par
  par <- .nm_apply_fix(model, par)
  list(par = par, value = opt$value, convergence = opt$convergence, optim = opt)
}

#' @keywords internal
.nm_fo_omega_prior_nll <- function(model, eta_mat, omega) {
  if (length(omega) == 0L || nrow(eta_mat) == 0L) {
    return(.nm_zero_like(omega))
  }
  nll <- .nm_zero_like(omega)
  for (i in seq_len(nrow(eta_mat))) {
    nll <- .ad_add(nll, .nm_omega_prior(eta_mat[i, ], omega))
  }
  nll
}

#' @keywords internal
.nm_fo_subject_G <- function(model, subj, theta, sigma, pk_engine = "cpp") {
  n_eta <- .nm_n_eta(model)
  if (n_eta == 0L) {
    return(NULL)
  }
  omega0 <- rep(0.01, n_eta)
  eta0 <- rep(0, n_eta)
  pred0 <- .nm_subject_ipred(
    model, subj, theta, omega0, eta0, sigma, pk_engine = pk_engine
  )
  f0 <- as.numeric(pred0$F)
  n_obs <- length(f0)
  if (n_obs == 0L) {
    return(matrix(0, nrow = 0L, ncol = n_eta))
  }
  G <- matrix(0, nrow = n_obs, ncol = n_eta)
  eps <- .nm_fo_eta_eps()
  for (k in seq_len(n_eta)) {
    etap <- rep(0, n_eta)
    etam <- rep(0, n_eta)
    etap[k] <- eps
    etam[k] <- -eps
    fp <- as.numeric(.nm_subject_ipred(
      model, subj, theta, omega0, etap, sigma, pk_engine = pk_engine
    )$F)
    fm <- as.numeric(.nm_subject_ipred(
      model, subj, theta, omega0, etam, sigma, pk_engine = pk_engine
    )$F)
    G[, k] <- (fp - fm) / (2 * eps)
  }
  G
}

#' @keywords internal
.nm_fo_subject_S <- function(G, f0, sigma) {
  if (is.null(G) || nrow(G) == 0L) {
    return(matrix(0, nrow = ncol(G), ncol = ncol(G)))
  }
  s1 <- .nm_sigma_el(sigma, 1L)
  s2 <- .nm_sigma_el(sigma, 2L)
  var <- pmax((s1 * f0)^2 + s2^2, .Machine$double.eps)
  scale <- sqrt(1 / var)
  Gs <- G * scale
  crossprod(Gs)
}

#' @keywords internal
.nm_fo_prep_subjects <- function(model, data, theta, sigma, pk_engine = "cpp") {
  n_eta <- .nm_n_eta(model)
  if (n_eta == 0L) {
    return(list())
  }
  pk_engine <- .nm_resolve_pk_engine(
    pk_engine, model, theta, rep(0.01, n_eta), sigma, rep(0, n_eta)
  )
  dat <- .nm_prepare_data(data, model$INPUT, model)
  ids <- .nm_subject_ids(dat)
  lapply(ids, function(id) {
    subj <- .nm_subject_slice(dat, id)
    pred0 <- .nm_subject_ipred(
      model, subj, theta, rep(0.01, n_eta), rep(0, n_eta), sigma,
      pk_engine = pk_engine
    )
    f0 <- as.numeric(pred0$F)
    G <- .nm_fo_subject_G(model, subj, theta, sigma, pk_engine)
    list(S = .nm_fo_subject_S(G, f0, sigma))
  })
}

#' @keywords internal
.nm_fo_logdet_sym <- function(M) {
  if (is.list(M) && !is.null(M[[1L]]) && is.list(M[[1L]])) {
    return(.nm_fo_logdet_sym_ad(M))
  }
  if (.nm_any_ad(M)) {
    return(.nm_fo_logdet_sym_ad(M))
  }
  n <- nrow(M)
  if (n == 1L) {
    return(log(max(M[1L, 1L], .Machine$double.eps)))
  }
  ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  sum(log(pmax(Re(ev), .Machine$double.eps)))
}

#' @keywords internal
.nm_fo_logdet_sym_ad <- function(M) {
  n <- length(M)
  z <- M[[1L]][[1L]]
  if (n == 1L) {
    return(.ad_log(.ad_pmax(M[[1L]][[1L]], newConstant(name = "ld_eps", value = 1e-15))))
  }
  L <- vector("list", n)
  for (i in seq_len(n)) {
    L[[i]] <- vector("list", n)
    for (j in seq_len(n)) {
      L[[i]][[j]] <- M[[i]][[j]]
    }
  }
  for (i in seq_len(n)) {
    for (j in seq_len(i)) {
      s <- z
      if (j > 1L) {
        for (k in seq_len(j - 1L)) {
          s <- .ad_add(s, .ad_mul(L[[i]][[k]], L[[j]][[k]]))
        }
      }
      if (i == j) {
        L[[i]][[j]] <- .ad_sqrt(.ad_pmax(
          .ad_sub(M[[i]][[j]], s),
          newConstant(name = "chol_eps", value = 1e-15)
        ))
      } else {
        L[[i]][[j]] <- .ad_div(.ad_sub(M[[i]][[j]], s), L[[j]][[j]])
      }
    }
  }
  logdet <- .ad_log(L[[1L]][[1L]])
  if (n > 1L) {
    for (i in 2L:n) {
      logdet <- .ad_add(logdet, .ad_log(L[[i]][[i]]))
    }
  }
  .ad_mul(2, logdet)
}

#' FO omega term: \eqn{\sum_i \log|\Omega + G_i^{\top} R_i^{-1} G_i|}
#'
#' @keywords internal
.nm_fo_omega_nll <- function(model, fo_cache, omega) {
  if (length(fo_cache) == 0L) {
    return(.nm_zero_like(omega))
  }
  n_om <- length(omega)
  nll <- .nm_zero_like(omega)
  ad <- .nm_any_ad(omega)
  for (sub in fo_cache) {
    S <- sub$S
    if (ad) {
      M <- vector("list", n_om)
      for (i in seq_len(n_om)) {
        M[[i]] <- vector("list", n_om)
        for (j in seq_len(n_om)) {
          M[[i]][[j]] <- if (i == j) {
            .ad_add(
              newConstant(name = paste0("S_", i, "_", j), value = S[i, j]),
              omega[[i]]
            )
          } else {
            newConstant(name = paste0("S_", i, "_", j), value = S[i, j])
          }
        }
      }
      nll <- .ad_add(nll, .nm_fo_logdet_sym(M))
    } else {
      M <- S + diag(as.numeric(omega), n_om)
      nll <- nll + .nm_fo_logdet_sym(M)
    }
  }
  nll
}

#' Full FO -2LL: residual at eta = 0 plus omega prior at fitted etas
#'
#' @keywords internal
.nm_fo_nll_internal <- function(model,
                                  data,
                                  theta,
                                  omega,
                                  sigma,
                                  eta_mat = NULL,
                                  pk_engine = "auto") {
  res <- .nm_nll_internal(
    model, data, theta, omega, sigma,
    eta = NULL, include_omega_prior = FALSE, pk_engine = pk_engine
  )
  if (is.null(eta_mat)) {
    return(res)
  }
  om <- .nm_fo_omega_prior_nll(model, eta_mat, omega)
  if (.nm_any_ad(res, om)) {
    return(.ad_add(res, om))
  }
  res + om
}
