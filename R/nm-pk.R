#' @keywords internal
.nm_ke <- function(CL, V) {
  CL / V
}

#' @keywords internal
.nm_cmp1bol_ss <- function(dat, ke, V, obs_cmp = 1L) {
  n <- nrow(dat)
  A <- rep(0, n)
  for (i in seq_len(n)) {
    if (dat$EVID[i] == 1L) {
      A[i] <- dat$AMT[i]
    } else if (i > 1L) {
      dt <- dat$TIME[i] - dat$TIME[i - 1L]
      A[i] <- A[i - 1L] * exp(-ke * dt)
    }
  }
  A / V
}

#' @keywords internal
.nm_cmp2oral_ss <- function(dat, ka, ke, V, obs_cmp = 1L) {
  n <- nrow(dat)
  A_dep <- rep(0, n)
  A_cen <- rep(0, n)
  for (i in seq_len(n)) {
    if (dat$EVID[i] == 1L) {
      A_dep[i] <- dat$AMT[i]
    } else if (i > 1L) {
      dt <- dat$TIME[i] - dat$TIME[i - 1L]
      d <- exp(-ke * dt) - exp(-ka * dt)
      if (abs(ka - ke) < 1e-12) {
        fac <- ka * dt * exp(-ka * dt)
      } else {
        fac <- ka / (ka - ke) * d
      }
      A_cen[i] <- A_cen[i - 1L] * exp(-ke * dt) + A_dep[i - 1L] * fac
      A_dep[i] <- A_dep[i - 1L] * exp(-ka * dt)
    }
  }
  if (obs_cmp == 1L) A_dep / V else A_cen / V
}

#' Two-compartment IV bolus (analytical hybrid-eigenvalue solution).
#'
#' Mirrors the validated C++ core (`step_2_iv_bolus` in `nm_pk_pkadvan.h`):
#' the amount vector is propagated over each interval with the exact
#' matrix-exponential of the linear system rather than the previous crude
#' first-order approximation. Doses are applied additively after propagation
#' so multi-dose designs and residual drug are handled correctly.
#' @keywords internal
.nm_cmp2bol_ss <- function(dat, k10, k12, k21, V1, V2, obs_cmp = 1L) {
  n <- nrow(dat)
  A1 <- rep(0, n)
  A2 <- rep(0, n)
  E1 <- k10 + k12
  E2 <- k21
  s <- E1 + E2
  disc <- s * s - 4 * (E1 * E2 - k12 * k21)
  root <- sqrt(max(disc, 0))
  l1 <- 0.5 * (s + root)
  l2 <- 0.5 * (s - root)
  d <- l2 - l1
  has_cmt <- "CMT" %in% names(dat)
  a1 <- 0
  a2 <- 0
  for (i in seq_len(n)) {
    if (i > 1L) {
      dt <- dat$TIME[i] - dat$TIME[i - 1L]
      p1 <- a1
      p2 <- a2
      if (abs(d) < 1e-12) {
        e <- exp(-l1 * dt)
        a1 <- (p1 + p2 * k21 / l1) * e
        a2 <- p2 * e
      } else {
        a1 <- (((p1 * E2 + p2 * k21) - p1 * l1) * exp(-l1 * dt) -
          ((p1 * E2 + p2 * k21) - p1 * l2) * exp(-l2 * dt)) / d
        a2 <- (((p2 * E1 + p1 * k12) - p2 * l1) * exp(-l1 * dt) -
          ((p2 * E1 + p1 * k12) - p2 * l2) * exp(-l2 * dt)) / d
      }
    }
    if (isTRUE(dat$EVID[i] == 3L)) {
      a1 <- 0
      a2 <- 0
    } else if (isTRUE(dat$EVID[i] == 1L)) {
      cmt <- if (has_cmt) as.integer(dat$CMT[i]) else 1L
      if (identical(cmt, 2L)) {
        a2 <- a2 + dat$AMT[i]
      } else {
        a1 <- a1 + dat$AMT[i]
      }
    }
    A1[i] <- a1
    A2[i] <- a2
  }
  if (obs_cmp == 2L) A2 / V2 else A1 / V1
}

#' @keywords internal
.nm_cmp2oral_trans4 <- function(dat, ka, CL, VC, VP, Q2, obs_cmp = 1L) {
  k20 <- CL / VC
  k23 <- Q2 / VC
  k32 <- Q2 / VP
  ad_mode <- .nm_any_ad(ka, CL, VC, VP, Q2)
  st <- .nm_cmp2oral_trans4_state(dat, ka, k20, k23, k32, ad_mode)
  if (!ad_mode) {
    return(if (obs_cmp == 1L) {
      st$gut / VC
    } else if (obs_cmp == 3L) {
      st$a3 / VP
    } else {
      st$a2 / VC
    })
  }
  if (obs_cmp == 1L) {
    vals <- lapply(st$gut, function(x) .ad_div(x, VC))
  } else if (obs_cmp == 3L) {
    vals <- lapply(st$a3, function(x) .ad_div(x, VP))
  } else {
    vals <- lapply(st$a2, function(x) .ad_div(x, VC))
  }
  structure(list(ad = TRUE, values = vals), class = "nm_ipred_ad")
}

#' @keywords internal
.nm_predict_subject_cpp <- function(model, subj, pred_vals) {
  ev <- .nm_input_event_vectors(subj)
  nm_pk_route_r(
    advan = as.integer(model$ADVAN),
    trans = as.integer(model$TRANS),
    obs_cmp = as.integer(model$OBSCMP),
    dose_cmp = as.integer(model$DOSECMP),
    n_transit = .nm_n_transit(model),
    use_ode = isTRUE(model$USE_ODE),
    model_ss = as.integer(model$SS),
    time = as.numeric(subj$TIME),
    amt = as.numeric(subj$AMT),
    rate = as.numeric(if ("RATE" %in% names(subj)) subj$RATE else 0),
    f1 = as.numeric(subj$F1),
    cmt = as.integer(if ("CMT" %in% names(subj)) subj$CMT else model$DOSECMP),
    evid = as.integer(subj$EVID),
    ss = as.integer(if ("SS" %in% names(subj)) subj$SS else 0L),
    ii = as.numeric(if ("II" %in% names(subj)) subj$II else subj$TAU),
    pk_params = pred_vals,
    s1 = ev$s1,
    s2 = ev$s2,
    s3 = ev$s3,
    s4 = ev$s4,
    scale_mat = ev$scale_mat,
    use_data_scale = ev$use_data_scale,
    f_mat = ev$f_mat,
    use_data_f = ev$use_data_f
  )
}

#' @keywords internal
.nm_predict_subject_cpp_detail <- function(model, subj, pred_vals, n_state = 7L) {
  ev <- .nm_input_event_vectors(subj)
  nm_pk_route_detail_r(
    advan = as.integer(model$ADVAN),
    trans = as.integer(model$TRANS),
    obs_cmp = as.integer(model$OBSCMP),
    dose_cmp = as.integer(model$DOSECMP),
    n_transit = .nm_n_transit(model),
    use_ode = isTRUE(model$USE_ODE),
    model_ss = as.integer(model$SS),
    n_state = as.integer(n_state),
    time = as.numeric(subj$TIME),
    amt = as.numeric(subj$AMT),
    rate = as.numeric(if ("RATE" %in% names(subj)) subj$RATE else 0),
    f1 = as.numeric(subj$F1),
    cmt = as.integer(if ("CMT" %in% names(subj)) subj$CMT else model$DOSECMP),
    evid = as.integer(subj$EVID),
    ss = as.integer(if ("SS" %in% names(subj)) subj$SS else 0L),
    ii = as.numeric(if ("II" %in% names(subj)) subj$II else subj$TAU),
    pk_params = pred_vals,
    s1 = ev$s1,
    s2 = ev$s2,
    s3 = ev$s3,
    s4 = ev$s4,
    scale_mat = ev$scale_mat,
    use_data_scale = ev$use_data_scale,
    f_mat = ev$f_mat,
    use_data_f = ev$use_data_f
  )
}

#' @keywords internal
.nm_predict_subject <- function(model, subj, pred_vals, pk_engine = "R") {
  if (pk_engine == "cpp") {
    return(.nm_predict_subject_cpp(model, subj, pred_vals))
  }
  advan <- model$ADVAN
  trans <- model$TRANS
  obs <- model$OBSCMP
  if (advan == 1L && trans == 2L) {
    ke <- pred_vals$KE %||% .nm_ke(pred_vals$CL, pred_vals$V)
    V <- pred_vals$V %||% pred_vals$VC
    ad_mode <- .nm_any_ad(ke, V)
    if (ad_mode && .nm_use_cpp_pk_ad()) {
      conc <- pk_bolus1_block_var(
        as.numeric(subj$TIME),
        as.numeric(subj$AMT),
        as.integer(subj$EVID),
        ke, V,
        as.integer(obs)
      )
      return(structure(
        list(ad = TRUE, values = conc, vector = TRUE),
        class = "nm_ipred_ad"
      ))
    }
    return(.nm_cmp1bol_ss(subj, ke, V, obs))
  }
  if (advan == 2L && trans %in% c(1L, 2L)) {
    ka <- pred_vals$KA
    ke <- pred_vals$KE %||% .nm_ke(pred_vals$CL, pred_vals$V)
    V <- pred_vals$V %||% pred_vals$VC
    ad_mode <- .nm_any_ad(ka, ke, V)
    if (ad_mode && .nm_use_cpp_pk_ad()) {
      conc <- pk_oral1_block_var(
        as.numeric(subj$TIME),
        as.numeric(subj$AMT),
        as.numeric(subj$F1),
        as.integer(subj$EVID),
        ka, ke, V,
        as.integer(obs)
      )
      return(structure(
        list(ad = TRUE, values = conc, vector = TRUE),
        class = "nm_ipred_ad"
      ))
    }
    return(.nm_cmp2oral_ss(subj, ka, ke, V, obs))
  }
  if ((advan == 3L || advan == 4L) && trans == 4L) {
    ka <- pred_vals$KA
    CL <- pred_vals$CL
    VC <- pred_vals$VC
    VP <- pred_vals$VP
    Q2 <- pred_vals$Q2
    ad_mode <- .nm_any_ad(ka, CL, VC, VP, Q2)
    if (ad_mode && .nm_use_cpp_pk_ad()) {
      k20 <- .ad_div(CL, VC)
      k23 <- .ad_div(Q2, VC)
      k32 <- .ad_div(Q2, VP)
      conc <- pk_oral2_trans4_block_var(
        as.numeric(subj$TIME),
        as.numeric(subj$AMT),
        as.numeric(subj$F1),
        as.integer(subj$EVID),
        ka, k20, k23, k32, VC, VP,
        as.integer(obs)
      )
      return(structure(
        list(ad = TRUE, values = conc, vector = TRUE),
        class = "nm_ipred_ad"
      ))
    }
    return(.nm_cmp2oral_trans4(
      subj, ka, CL, VC, VP, Q2, obs
    ))
  }
  if (advan == 3L && trans == 2L) {
    k10 <- pred_vals$CL / pred_vals$VC
    k12 <- pred_vals$Q2 / pred_vals$VC
    k21 <- pred_vals$Q2 / pred_vals$VP
    return(.nm_cmp2bol_ss(subj, k10, k12, k21, pred_vals$VC, pred_vals$VP, obs))
  }
  .nm_stop("Unsupported ADVAN/TRANS combination: ADVAN=", advan, " TRANS=", trans)
}

#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' @keywords internal
.nm_ipred_values <- function(ipred) {
  if (inherits(ipred, "nm_ipred_ad")) {
    if (isTRUE(ipred$vector)) {
      return(as.numeric(ipred$values$value))
    }
    return(vapply(ipred$values, function(x) x$value, numeric(1)))
  }
  ipred
}

#' @keywords internal
.nm_ipred_at <- function(ipred, idx) {
  if (inherits(ipred, "nm_ipred_ad")) {
    if (isTRUE(ipred$vector)) {
      if (length(idx) == 1L) {
        return(LibeRtAD:::subset_var(ipred$values, as.integer(idx)))
      }
      return(lapply(idx, function(i) LibeRtAD:::subset_var(ipred$values, as.integer(i))))
    }
    return(ipred$values[idx])
  }
  ipred[idx]
}

#' @keywords internal
.nm_subject_ipred <- function(model, subj, theta, omega, eta,
                              sigma = NULL, pk_engine = "auto",
                              with_amounts = FALSE, n_state = 7L) {
  if (is.null(sigma)) {
    sigma <- model$SIGMAS$Value
  }
  pk_engine <- .nm_resolve_pk_engine(pk_engine, model, theta, omega, sigma, eta)
  subj_ev <- .nm_subject_events(subj)
  pred_vals <- if (pk_engine == "cpp") {
    cov <- list()
    if (!is.null(model$COVARIATES) && length(model$COVARIATES) > 0L) {
      cv <- as.character(model$COVARIATES)
      missing_cv <- setdiff(cv, names(subj))
      if (length(missing_cv) > 0L) {
        .nm_stop("Covariate columns missing from data: ", paste(missing_cv, collapse = ", "))
      }
      cov <- stats::setNames(
        lapply(cv, function(cn) .nm_cov_baseline_value(subj[[cn]], cn)),
        cv
      )
    }
    nm_eval_pred_cpp(
      .nm_split_lines(model$PRED), as.numeric(theta), as.numeric(eta),
      covariates = cov,
      des_lines = if (isTRUE(model$USE_ODE)) {
        character()
      } else {
        .nm_split_lines(model$DES)
      }
    )
  } else {
    .nm_eval_pred(model, theta, omega, eta)
  }
  ipred_ev <- if (isTRUE(with_amounts) && pk_engine == "cpp") {
    detail <- .nm_predict_subject_cpp_detail(model, subj_ev, pred_vals, n_state = n_state)
    list(
      ipred = as.numeric(detail$ipred),
      amounts = detail$amounts,
      plain = TRUE
    )
  } else {
    .nm_predict_subject(model, subj_ev, pred_vals, pk_engine = pk_engine)
  }
  obs_idx <- which(subj_ev$MDV == 0L & subj_ev$EVID == 0L)
  ipred_vals <- if (is.list(ipred_ev) && isTRUE(ipred_ev$plain)) {
    ipred_ev$ipred
  } else {
    .nm_ipred_values(ipred_ev)
  }
  f_vals <- if (is.list(ipred_ev) && isTRUE(ipred_ev$plain)) {
    ipred_vals[obs_idx]
  } else {
    .nm_ipred_at(ipred_ev, obs_idx)
  }
  out <- list(
    ipred = .nm_ipred_align(subj, subj_ev, ipred_vals),
    F = f_vals,
    obs_idx = obs_idx,
    subj_ev = subj_ev,
    pred_vals = pred_vals
  )
  if (is.list(ipred_ev) && isTRUE(ipred_ev$plain) && !is.null(ipred_ev$amounts)) {
    out$amounts <- ipred_ev$amounts
  }
  out
}

#' @keywords internal
.nm_sigma_el <- function(sigma, i) {
  if (length(sigma) < i) {
    return(0)
  }
  if (is.list(sigma)) {
    return(sigma[[i]])
  }
  sigma[i]
}

#' @keywords internal
.nm_residual_nll <- function(dv, f, sigma, error = "propadd", dvid = NULL,
                             ar1_rho = 0.0, sigma_corr = "indep") {
  s1 <- .nm_sigma_el(sigma, 1L)
  s2 <- .nm_sigma_el(sigma, 2L)
  n_obs <- length(dv)
  if (is.list(f)) {
    nll <- .nm_zero_like(f[[1]], s1, s2)
    s1sq <- .ad_mul(s1, s1)
    s2sq <- .ad_mul(s2, s2)
    for (k in seq_len(n_obs)) {
      fk <- f[[k]]
      resid <- .ad_sub(dv[k], fk)
      var <- .ad_pmax(
        .ad_add(.ad_mul(.ad_mul(fk, fk), s1sq), s2sq),
        newConstant(name = "var_eps", value = .Machine$double.eps)
      )
      nll <- .ad_add(nll, .ad_add(.ad_log(var), .ad_div(.ad_mul(resid, resid), var)))
    }
    return(nll)
  }
  resid <- dv - f
  if (.nm_any_ad(f, sigma, s1, s2)) {
    var <- f * f * s1 * s1 + s2 * s2
    nll <- .nm_zero_like(f, s1, s2)
    for (k in seq_len(n_obs)) {
      nll <- nll + .ad_log(var[k]) + resid[k]^2 / var[k]
    }
    return(nll)
  }
  .nm_residual_nll_scalar(dv, f, sigma, error, dvid, ar1_rho, sigma_corr)
}

#' @keywords internal
.nm_omega_prior <- function(eta, omega) {
  if (length(eta) == 0L) {
    return(0)
  }
  if (is.list(omega) || is.list(eta)) {
    n <- max(length(omega), if (is.list(eta)) length(eta) else length(eta))
    nll <- .nm_zero_like(omega, eta)
    for (i in seq_len(n)) {
      ei <- if (is.list(eta)) eta[[i]] else eta[i]
      omi <- if (is.list(omega)) omega[[i]] else omega[i]
      if (!.nm_any_ad(ei, omi) && !is.list(omega)) {
        omi <- max(omi, .Machine$double.eps)
      }
      if (.nm_any_ad(ei, omi)) {
        nll <- .ad_add(nll, .ad_add(.ad_div(.ad_mul(ei, ei), omi), log(omi)))
      } else {
        nll <- nll + ei^2 / omi + log(omi)
      }
    }
    return(nll)
  }
  om <- if (.nm_any_ad(eta, omega)) omega else pmax(omega, .Machine$double.eps)
  sum(eta^2 / om + log(om))
}
