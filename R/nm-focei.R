#' @keywords internal
.nm_focei_is_ad_mat <- function(M) {
  if (!is.list(M) || length(M) == 0L) {
    return(.nm_any_ad(M))
  }
  if (is.list(M[[1L]])) {
    return(inherits(M[[1L]][[1L]], c("Variable", "Constant")))
  }
  .nm_any_ad(M)
}

#' @keywords internal
.nm_focei_is_diagonal <- function(M) {
  p <- if (is.list(M) && is.list(M[[1L]])) length(M) else nrow(M)
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      if (i == j) next
      off <- if (is.list(M) && is.list(M[[1L]])) {
        M[[i]][[j]]
      } else {
        M[i, j]
      }
      if (.nm_focei_is_ad_mat(M)) {
        if (!identical(as.numeric(.ad_scalar_value(off)), 0)) {
          return(FALSE)
        }
      } else if (abs(off) > 0) {
        return(FALSE)
      }
    }
  }
  TRUE
}

#' @keywords internal
.nm_focei_eta_eps <- function() {
  1e-4
}

#' @keywords internal
.nm_focei_pick_f <- function(f, j) {
  if (is.list(f) && !inherits(f, "Variable")) {
    return(f[[j]])
  }
  f[j]
}

#' @keywords internal
.nm_focei_n_obs <- function(f) {
  if (is.list(f) && !inherits(f, "Variable")) {
    return(length(f))
  }
  length(f)
}

#' Build omega matrix for FOCEI (numeric or AD list-matrix).
#' @keywords internal
.nm_focei_omega_matrix <- function(omega, n_eta, omega_struct = "diag") {
  eps <- newConstant(name = "om_eps", value = 1e-15)
  if (!.nm_any_ad(omega)) {
    if (identical(omega_struct, "block2") && n_eta >= 2L && length(omega) >= 3L) {
      OM <- matrix(0, n_eta, n_eta)
      OM[1L, 1L] <- max(omega[1L], 1e-15)
      OM[2L, 2L] <- max(omega[2L], 1e-15)
      OM[1L, 2L] <- OM[2L, 1L] <- omega[3L]
      if (n_eta > 2L) {
        diag(OM)[3L:n_eta] <- pmax(omega[4L:n_eta], 1e-15)
      }
      return(OM)
    }
    v <- pmax(omega[seq_len(n_eta)], 1e-15)
    if (n_eta == 1L) {
      return(matrix(v[1L], 1L, 1L))
    }
    return(diag(v))
  }
  om <- function(i, j, val) {
    if (i == j) .ad_pmax(val, eps) else val
  }
  M <- vector("list", n_eta)
  for (i in seq_len(n_eta)) {
    M[[i]] <- vector("list", n_eta)
    for (j in seq_len(n_eta)) {
      M[[i]][[j]] <- newConstant(name = paste0("om0_", i, "_", j), value = 0)
    }
  }
  if (identical(omega_struct, "block2") && n_eta >= 2L && length(omega) >= 3L) {
    M[[1L]][[1L]] <- om(1L, 1L, omega[[1L]])
    M[[2L]][[2L]] <- om(2L, 2L, omega[[2L]])
    M[[1L]][[2L]] <- M[[2L]][[1L]] <- omega[[3L]]
    if (n_eta > 2L) {
      for (k in 3L:n_eta) {
        idx <- if (k <= length(omega)) k else 1L
        M[[k]][[k]] <- om(k, k, omega[[idx]])
      }
    }
    return(M)
  }
  for (k in seq_len(n_eta)) {
    idx <- if (k <= length(omega)) k else 1L
    M[[k]][[k]] <- om(k, k, omega[[idx]])
  }
  M
}

#' @keywords internal
.nm_focei_invert_spd <- function(M) {
  ad <- .nm_focei_is_ad_mat(M)
  p <- if (ad) {
    length(M)
  } else if (is.matrix(M)) {
    nrow(M)
  } else {
    length(M)
  }
  if (p == 0L) {
    return(M)
  }
  if (!ad) {
    if (p == 1L) {
      return(matrix(1 / M[1L, 1L], 1L, 1L))
    }
    if (p == 2L) {
      a <- M[1L, 1L]
      b <- M[1L, 2L]
      c <- M[2L, 2L]
      det <- max(a * c - b * b, 1e-15)
      return(matrix(c(c / det, -b / det, -b / det, a / det), 2L, 2L, byrow = TRUE))
    }
    return(solve(M))
  }
  eps <- newConstant(name = "inv_eps", value = 1e-15)
  if (p == 1L) {
    return(list(list(.ad_div(1, M[[1L]][[1L]]))))
  }
  if (p == 2L) {
    a <- M[[1L]][[1L]]
    b <- M[[1L]][[2L]]
    c <- M[[2L]][[2L]]
    det <- .ad_pmax(.ad_sub(.ad_mul(a, c), .ad_mul(b, b)), eps)
    inv <- vector("list", 2L)
    inv[[1L]] <- list(.ad_div(c, det), .ad_div(.ad_mul(-1, b), det))
    inv[[2L]] <- list(.ad_div(.ad_mul(-1, b), det), .ad_div(a, det))
    return(inv)
  }
  if (.nm_focei_is_diagonal(M)) {
    inv <- vector("list", p)
    for (i in seq_len(p)) {
      inv[[i]] <- vector("list", p)
      for (j in seq_len(p)) {
        inv[[i]][[j]] <- if (i == j) {
          .ad_div(1, M[[i]][[j]])
        } else {
          newConstant(name = paste0("inv0_", i, "_", j), value = 0)
        }
      }
    }
    return(inv)
  }
  .nm_stop(
    "AD FOCEI matrix inverse for non-diagonal Omega with >2 ETAs is not supported; ",
    "use grad = \"numeric\"."
  )
}

#' @keywords internal
.nm_focei_mat_add <- function(A, B) {
  if (!.nm_focei_is_ad_mat(A) && !.nm_focei_is_ad_mat(B)) {
    return(A + B)
  }
  p <- length(A)
  M <- vector("list", p)
  for (i in seq_len(p)) {
    M[[i]] <- vector("list", p)
    for (j in seq_len(p)) {
      M[[i]][[j]] <- .ad_add(A[[i]][[j]], B[[i]][[j]])
    }
  }
  M
}

#' @keywords internal
.nm_focei_gvg_init <- function(p, seed) {
  z <- .nm_zero_like(seed)
  if (!.nm_any_ad(seed)) {
    return(matrix(0, p, p))
  }
  M <- vector("list", p)
  for (i in seq_len(p)) {
    M[[i]] <- vector("list", p)
    for (j in seq_len(p)) {
      M[[i]][[j]] <- z
    }
  }
  M
}

#' @keywords internal
.nm_focei_omega_prior <- function(eta, OM, invOM) {
  eta <- as.numeric(eta)
  p <- length(eta)
  quad <- .nm_zero_like(OM, invOM)
  if (.nm_focei_is_ad_mat(invOM)) {
    for (i in seq_len(p)) {
      for (j in seq_len(p)) {
        quad <- .ad_add(
          quad,
          .ad_mul(.ad_mul(eta[i], invOM[[i]][[j]]), eta[j])
        )
      }
    }
  } else {
    for (i in seq_len(p)) {
      for (j in seq_len(p)) {
        quad <- quad + eta[i] * invOM[i, j] * eta[j]
      }
    }
  }
  .ad_add(quad, .nm_fo_logdet_sym(OM))
}

#' FOCE residual variance R_j (no eta variance in V).
#' @keywords internal
.nm_focei_residual_var <- function(f, sigma, error = "propadd") {
  s1 <- .nm_sigma_el(sigma, 1L)
  s2 <- .nm_sigma_el(sigma, 2L)
  eps <- newConstant(name = "focei_rv_eps", value = 1e-15)
  if (.nm_any_ad(f, s1, s2)) {
    s1sq <- .ad_mul(s1, s1)
    s2sq <- .ad_mul(s2, s2)
    return(switch(
      error,
      add = .ad_pmax(s1sq, eps),
      prop = .ad_pmax(.ad_mul(.ad_mul(f, f), s1sq), eps),
      log = .ad_pmax(s1sq, eps),
      power = .ad_pmax(
        .ad_mul(
          s1sq,
          .ad_exp(.ad_mul(.ad_mul(2, s2), .ad_log(.ad_pmax(f, eps))))
        ),
        eps
      ),
      .ad_pmax(.ad_add(.ad_mul(.ad_mul(f, f), s1sq), s2sq), eps)
    ))
  }
  .nm_residual_var(f, sigma, error)
}

#' FOCE-INTER curvature addition: \eqn{0.5 (R'_j / R_j)^2}.
#'
#' R depends on eta only through the individual prediction \eqn{f_j}, so
#' \eqn{\partial R_j/\partial\eta = (dR_j/df_j)\, G_j}. This returns the extra
#' factor (beyond the FOCE \eqn{1/R_j}) multiplying \eqn{G_{ji}G_{jk}} in the
#' curvature/information term. Zero for additive/log error (no interaction).
#' @keywords internal
.nm_focei_interaction_coef <- function(fj, sigma, rj, error, ad) {
  s1 <- .nm_sigma_el(sigma, 1L)
  s2 <- .nm_sigma_el(sigma, 2L)
  if (ad) {
    rp <- switch(
      error,
      prop = .ad_div(newConstant(name = "two_pr", value = 2), fj),
      propadd = .ad_div(
        .ad_mul(.ad_mul(newConstant(name = "two_pa", value = 2), .ad_mul(s1, s1)), fj),
        rj
      ),
      power = .ad_div(.ad_mul(newConstant(name = "two_pw", value = 2), s2), fj),
      newConstant(name = "zero_int", value = 0)
    )
    return(.ad_mul(newConstant(name = "half_int", value = 0.5), .ad_mul(rp, rp)))
  }
  fj <- as.numeric(fj)
  rj <- as.numeric(rj)
  rp <- switch(
    error,
    prop = if (abs(fj) > 1e-12) 2 / fj else 0,
    propadd = if (rj > 0) 2 * as.numeric(s1)^2 * fj / rj else 0,
    power = if (abs(fj) > 1e-12) 2 * as.numeric(s2) / fj else 0,
    0
  )
  0.5 * rp^2
}

#' @keywords internal
.nm_focei_subject_G <- function(model, subj, theta, omega, eta, sigma, pk_engine) {
  .nm_focei_subject_G_sens(
    model, subj, theta, omega, eta, sigma, pk_engine
  )
}

#' @keywords internal
.nm_focei_G_ij <- function(G, j, i, ad = NULL) {
  if (is.null(ad)) {
    ad <- is.list(G) && length(G) > 0L && is.list(G[[1L, 1L]])
  }
  if (ad) {
    g <- G[j, i]
    if (is.list(g) && length(g) == 1L) {
      return(g[[1L]])
    }
    return(g)
  }
  G[j, i]
}

#' Per-subject FOCEI -2LL at a fixed eta mode (R path, AD-capable).
#' @keywords internal
.nm_focei_subject_nll_internal <- function(model,
                                           subj,
                                           theta,
                                           omega,
                                           sigma,
                                           eta,
                                           pk_engine = "auto",
                                           interaction = TRUE) {
  eta <- as.numeric(eta)
  n_eta <- length(eta)
  if (n_eta == 0L || length(omega) == 0L) {
    return(.nm_zero_like(theta, omega, sigma))
  }
  cfg <- .nm_lik_config(model)
  gh <- .nm_focei_subject_G(
    model, subj, theta, omega, eta, sigma, pk_engine = pk_engine
  )
  f <- gh$F
  dv <- gh$dv
  G <- gh$G
  n_obs <- .nm_focei_n_obs(f)
  if (n_obs == 0L) {
    return(.nm_zero_like(theta, omega, sigma, f))
  }
  ad <- .nm_any_ad(theta, omega, sigma, f)
  OM <- .nm_focei_omega_matrix(omega, n_eta, cfg$omega)
  invOM <- .nm_focei_invert_spd(OM)
  nll <- .nm_focei_omega_prior(eta, OM, invOM)
  gvg <- .nm_focei_gvg_init(n_eta, if (ad) invOM[[1L]][[1L]] else nll)
  for (j in seq_len(n_obs)) {
    fj <- .nm_focei_pick_f(f, j)
    rj <- .nm_focei_residual_var(fj, sigma, cfg$error)
    inv_r <- if (ad) .ad_div(1, rj) else 1 / rj
    resid <- if (ad) .ad_sub(dv[j], fj) else dv[j] - fj
    nll <- if (ad) {
      .ad_add(nll, .ad_add(.ad_log(rj), .ad_div(.ad_mul(resid, resid), rj)))
    } else {
      nll + log(rj) + resid^2 / rj
    }
    # FOCE-INTER: augment the curvature coefficient 1/R_j with 0.5 (R'_j/R_j)^2.
    coef <- if (isTRUE(interaction)) {
      extra <- .nm_focei_interaction_coef(fj, sigma, rj, cfg$error, ad)
      if (ad) .ad_add(inv_r, extra) else inv_r + extra
    } else {
      inv_r
    }
    for (i in seq_len(n_eta)) {
      for (k in seq_len(n_eta)) {
        gij <- .nm_focei_G_ij(G, j, i, ad)
        gkj <- .nm_focei_G_ij(G, j, k, ad)
        if (ad) {
          gvg[[i]][[k]] <- .ad_add(
            gvg[[i]][[k]],
            .ad_mul(.ad_mul(coef, gij), gkj)
          )
        } else {
          gvg[i, k] <- gvg[i, k] + coef * gij * gkj
        }
      }
    }
  }
  M <- .nm_focei_mat_add(invOM, gvg)
  .ad_add(nll, .nm_fo_logdet_sym(M))
}

#' Population FOCEI -2LL (C++ when numeric, R when AD active).
#' @keywords internal
.nm_focei_nll_internal <- function(model,
                                   data,
                                   theta,
                                   omega,
                                   sigma,
                                   eta_mat,
                                   pk_engine = "auto") {
  if (is.null(eta_mat) || ncol(eta_mat) == 0L) {
    return(.nm_nll_internal(
      model, data, theta, omega, sigma,
      eta = eta_mat, include_omega_prior = TRUE, pk_engine = pk_engine
    ))
  }
  pk_engine <- .nm_resolve_pk_engine(pk_engine, model, theta, omega, sigma, eta_mat)
  if (pk_engine == "cpp" && !.nm_any_ad(theta, omega, sigma)) {
    meta <- .nm_cpp_meta(model)
    subs <- .nm_cpp_subjects_cached(model, data)
    return(nm_focei_objective_cpp(
      subjects = subs,
      eta_modes = eta_mat,
      theta = as.numeric(theta),
      omega = as.numeric(omega),
      sigma = as.numeric(sigma),
      pred_lines = meta$pred_lines,
      advan = meta$advan,
      trans = meta$trans,
      obs_cmp = meta$obs_cmp,
      dose_cmp = meta$dose_cmp,
      n_transit = meta$n_transit,
      use_ode = meta$use_ode,
      model_ss = meta$model_ss
    ))
  }
  dat <- .nm_prepare_data(data, model$INPUT, model)
  ids <- .nm_subject_ids(dat)
  if (nrow(eta_mat) != length(ids)) {
    .nm_stop("eta_mat must have one row per subject.")
  }
  nll <- .nm_zero_like(theta, omega, sigma, eta_mat)
  for (j in seq_along(ids)) {
    subj <- .nm_subject_slice(dat, ids[j])
    nll <- .ad_add(nll, .nm_focei_subject_nll_internal(
      model, subj, theta, omega, sigma, eta_mat[j, ], pk_engine = pk_engine
    ))
  }
  nll
}

#' NONMEM-style FOCEI objective at conditional eta modes.
#'
#' Uses FOCE residual variance (no eta variance in V) plus the log
#' determinant interaction term \eqn{\log\det(\Omega^{-1} + G^{\top} R^{-1} G)}.
#'
#' @keywords internal
.nm_focei_objective <- function(model, data, theta, omega, sigma, eta_mat, pk_engine = "auto") {
  .nm_focei_nll_internal(
    model, data, theta, omega, sigma, eta_mat, pk_engine = pk_engine
  )
}

#' @keywords internal
.nm_focei_nested_cache <- function() {
  .nm_env_cache("focei_nested")
}

#' @keywords internal
.nm_clear_focei_nested_cache <- function() {
  cache <- .nm_focei_nested_cache()
  rm(list = ls(envir = cache, all.names = TRUE), envir = cache)
}

#' Nested FOCEI objective: re-fit conditional eta modes at each evaluation.
#' @keywords internal
.nm_focei_nested_objective <- function(model, dat, data, par, eta_hint, pk_engine, control,
                                      setup = NULL, use_cache = TRUE) {
  if (isTRUE(use_cache)) {
    key <- .ad_optim_cache_key(par)
    cache <- .nm_focei_nested_cache()
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache))
    }
  }
  pp <- .nm_unpack(model, par)
  eta_init <- .nm_mceta_init(
    model, dat, pp$theta, pp$omega, pp$sigma, eta_hint, control
  )
  eta_mat <- .nm_fit_all_eta_modes(
    model, dat, pp$theta, pp$omega, pp$sigma, eta_init,
    pk_engine, control
  )
  eta_mat <- .nm_sanitize_eta_mat(eta_mat, pp$omega)
  out <- list(
    value = .nm_focei_objective(
      model, data, pp$theta, pp$omega, pp$sigma, eta_mat, pk_engine
    ),
    eta = eta_mat
  )
  if (isTRUE(use_cache)) {
    assign(key, out, envir = cache)
  }
  out
}

#' @keywords internal
.nm_focei_interaction <- function(model, data, theta, omega, sigma, eta_mat) {
  if (is.null(eta_mat) || ncol(eta_mat) == 0L) {
    return(0)
  }
  meta <- .nm_cpp_meta(model)
  subs <- .nm_cpp_subjects_cached(model, data)
  nm_focei_interaction_cpp(
    subjects = subs,
    eta_modes = eta_mat,
    theta = as.numeric(theta),
    omega = as.numeric(omega),
    sigma = as.numeric(sigma),
    pred_lines = meta$pred_lines,
    advan = meta$advan,
    trans = meta$trans,
    obs_cmp = meta$obs_cmp,
    dose_cmp = meta$dose_cmp,
    n_transit = meta$n_transit,
    use_ode = meta$use_ode,
    model_ss = meta$model_ss
  )
}

#' FOCE with interaction (FOCEI)
#'
#' @keywords internal
.nm_est_focei <- function(model, data, par0, backend = "cpp", grad = "auto",
                          pk_engine = "auto", engine = "auto", control = list(),
                          max_outer = 10L, tol = 1e-5, ...) {
  pk_eff <- .nm_effective_pk_engine(pk_engine, grad)
  focei_grad <- control$focei_grad %||% grad
  est_grad <- .nm_resolve_focei_est_grad(focei_grad, model)
  eta_mode <- .nm_resolve_focei_eta_mode(control)
  setup <- nm_focei_setup(model, data, pk_engine = pk_eff)
  .nm_clear_laplace_optim_cache()
  .nm_clear_laplace_eta_modes()
  dat <- setup$dat
  par <- .nm_apply_fix(model, par0)
  eta_mat <- NULL
  optim_runs <- list()
  fit <- NULL
  use_sens <- identical(est_grad, "sensitivity")
  use_outer_eta <- identical(eta_mode, "outer")
  opt_grad <- "numeric"
  prev_val <- NULL
  for (iter in seq_len(max_outer)) {
    .nm_est_progress_outer("FOCEI", iter, max_outer, if (!is.null(prev_val)) prev_val else NA_real_)
    .nm_clear_pop_optim_cache()
    .nm_clear_focei_nested_cache()
    eta_env <- new.env(parent = emptyenv())
    eta_env$eta <- eta_mat
    eta_fixed <- NULL
    if (use_outer_eta) {
      nested_start <- .nm_focei_nested_objective(
        model, dat, data, par, eta_mat, pk_eff, control, use_cache = FALSE
      )
      eta_fixed <- nested_start$eta
      eta_env$eta <- eta_fixed
    }
    nest_state <- new.env(parent = emptyenv())
    nest_state$key <- NULL
    nest_state$nested <- NULL
    fetch <- if (use_outer_eta) {
      function(par) {
        .nm_focei_eval_at_eta(model, data, par, eta_fixed, pk_eff)
      }
    } else {
      function(par) {
        nested <- .nm_focei_nested_fetch(
          nest_state, model, dat, data, par, eta_env$eta, pk_eff, control
        )
        eta_env$eta <- nested$eta
        nested
      }
    }
    objective <- function(par) {
      .nm_finite_obj(fetch(par)$value)
    }
    gr_fn <- if (use_sens) {
      function(par) {
        nested <- fetch(par)
        .nm_focei_sensitivity_grad(
          model, data, par, nested$eta, pk_eff, backend, focei_grad
        )
      }
    } else {
      NULL
    }
    fit <- .nm_optimize_par(
      model, data, par, objective, backend, opt_grad, pk_eff, control,
      ad_objective = NULL, ad_presolve = NULL, eta_mat = NULL,
      include_omega_prior = TRUE, cpp_pop_grad = FALSE, gr_fn = gr_fn
    )
    .nm_est_progress_outer_done("FOCEI", iter, max_outer, fit$value)
    optim_runs[[iter]] <- fit$optim
    par <- fit$par
    eta_mat <- eta_env$eta
    if (!is.null(prev_val) &&
        abs(prev_val - fit$value) < tol * (1 + abs(prev_val))) {
      break
    }
    prev_val <- fit$value
  }
  p <- .nm_unpack(model, par)
  nested <- .nm_focei_nested_objective(
    model, dat, data, par, eta_mat, pk_eff, control
  )
  eta_mat <- nested$eta
  final_obj <- nested$value
  structure(
    list(
      method = "FOCEI",
      par = par,
      theta = p$theta,
      omega = p$omega,
      sigma = p$sigma,
      eta = eta_mat,
      objective = final_obj,
      outer = vapply(optim_runs, function(o) o$value, numeric(1)),
      convergence = fit$convergence,
      grad = if (use_sens) {
        if (.nm_grad_uses_ad(focei_grad)) "ad" else "sensitivity"
      } else {
        est_grad
      },
      grad_requested = focei_grad,
      grad_backend = if (use_sens && .nm_grad_uses_ad(focei_grad)) {
        backend
      } else if (use_sens) {
        "sensitivity"
      } else {
        .nm_report_grad_backend(model, est_grad, backend)
      },
      pk_engine = pk_eff,
      optim = optim_runs,
      focei_eta = eta_mode,
      interaction = .nm_focei_interaction(
        model, data, p$theta, p$omega, p$sigma, eta_mat
      ),
      focei_setup = setup
    ),
    class = "nm_fit"
  )
}
