#' @keywords internal
.nm_focei_g_method <- function(method = NULL) {
  method <- method %||% getOption("LibeRation.focei_G", "auto")
  match.arg(method, c("auto", "ad", "fd", "shi"))
}

#' Shi-style adaptive FD step for FOCEI G (eta or F scale).
#' @keywords internal
.nm_focei_shi_eps <- function(x, f_ref = NULL) {
  tau <- 1e-3
  base <- sqrt(.Machine$double.eps)
  h <- base * pmax(abs(x), tau)
  if (!is.null(f_ref) && length(f_ref) > 0L) {
    fs <- max(abs(f_ref), na.rm = TRUE)
    if (is.finite(fs) && fs > 0) {
      h <- pmax(h, base * pmax(fs, tau))
    }
  }
  pmin(pmax(h, 1e-7), 0.05)
}

#' Models where eta sensitivities are discontinuous or AD-unreliable.
#' @keywords internal
.nm_focei_needs_fd_g <- function(model, subj = NULL) {
  pred <- paste(model$PRED, collapse = "\n")
  if (grepl("ALAG", pred, ignore.case = TRUE)) {
    return(TRUE)
  }
  if (grepl("DUR\\s*=", pred, ignore.case = TRUE)) {
    return(TRUE)
  }
  if (isTRUE(model$USE_ODE)) {
    return(TRUE)
  }
  if (!is.null(subj)) {
    user_f <- attr(subj, "user_f_cols")
    user_scale <- attr(subj, "user_scale_cols")
    if (length(user_f) > 0L || length(user_scale) > 0L) {
      return(TRUE)
    }
    if ("RATE" %in% names(subj) && any(as.numeric(subj$RATE) != 0, na.rm = TRUE)) {
      return(TRUE)
    }
  }
  FALSE
}

#' @keywords internal
.nm_focei_pick_g_method <- function(model, subj, pk_engine, method = NULL) {
  method <- .nm_focei_g_method(method)
  if (identical(method, "ad")) {
    return("ad")
  }
  if (identical(method, "fd")) {
    return("fd")
  }
  if (identical(method, "shi")) {
    return("shi")
  }
  if (.nm_focei_needs_fd_g(model, subj)) {
    return("shi")
  }
  # Default: Shi/central FD (C++ fast path). AD G is opt-in via focei_G = "ad".
  "shi"
}

#' Central FD with fixed or Shi step sizes.
#' @keywords internal
.nm_focei_subject_G_ad <- function(model, subj, theta, omega, eta, sigma,
                                   pk_engine, f, n_obs, n_eta) {
  pk_engine <- .nm_resolve_pk_engine(pk_engine, model, theta, omega, sigma, eta)
  fn <- function(e) {
    pred <- .nm_subject_ipred(
      model, subj, theta, omega, e, sigma, pk_engine = pk_engine
    )
    as.numeric(pred$F)
  }
  if (requireNamespace("numDeriv", quietly = TRUE)) {
    return(numDeriv::jacobian(fn, eta))
  }
  if (!requireNamespace("LibeRtAD", quietly = TRUE)) {
    return(NULL)
  }
  eta_names <- paste0("ETA", seq_len(n_eta))
  at <- stats::setNames(as.list(as.numeric(eta)), eta_names)
  backend <- if (.nm_cpp_capable(model)) "cpp" else "R"
  G <- matrix(0, n_obs, n_eta)
  for (j in seq_len(n_obs)) {
    fj_src <- paste0(
      "function(", paste(eta_names, collapse = ", "), ") ",
      "as.numeric(.nm_subject_ipred(model, subj, theta, omega, ",
      "list(", paste(eta_names, collapse = ", "), "), sigma, ",
      "pk_engine = pk_engine)$F)[", j, "L]"
    )
    fj <- eval(parse(text = fj_src))
    environment(fj) <- list2env(
      list(
        model = model, subj = subj, theta = theta, omega = omega,
        sigma = sigma, pk_engine = pk_engine,
        .nm_subject_ipred = .nm_subject_ipred
      ),
      parent = parent.env(environment())
    )
    gres <- tryCatch(
      LibeRtAD::backdiff(fj, at = at, backend = backend),
      error = function(e) NULL
    )
    if (is.null(gres) || is.null(gres$partials_flat)) {
      return(NULL)
    }
    G[j, ] <- unname(gres$partials_flat[eta_names])
  }
  G
}

#' Central FD with fixed or Shi step sizes.
#' @keywords internal
.nm_focei_subject_G_fd <- function(model, subj, theta, omega, eta, sigma,
                                   pk_engine, f, n_obs, n_eta, use_shi = TRUE) {
  ad <- .nm_any_ad(theta, omega, sigma, f)
  G <- if (ad) {
    matrix(list(), n_obs, n_eta)
  } else {
    matrix(0, n_obs, n_eta)
  }
  f_num <- if (ad) {
    vapply(seq_len(n_obs), function(j) {
      x <- .nm_focei_pick_f(f, j)
      if (.nm_any_ad(x)) {
        as.numeric(.ad_scalar_value(x))
      } else {
        as.numeric(x)
      }
    }, numeric(1))
  } else {
    as.numeric(f)
  }
  for (k in seq_len(n_eta)) {
    eps <- if (use_shi) {
      .nm_focei_shi_eps(eta[k], f_num)
    } else {
      .nm_focei_eta_eps()
    }
    etap <- eta
    etam <- eta
    etap[k] <- etap[k] + eps
    etam[k] <- etam[k] - eps
    fp <- .nm_subject_ipred(
      model, subj, theta, omega, etap, sigma, pk_engine = pk_engine
    )$F
    fm <- .nm_subject_ipred(
      model, subj, theta, omega, etam, sigma, pk_engine = pk_engine
    )$F
    half <- 2 * eps
    for (j in seq_len(n_obs)) {
      if (ad) {
        G[j, k] <- list(.ad_div(
          .ad_sub(.nm_focei_pick_f(fp, j), .nm_focei_pick_f(fm, j)),
          newConstant(name = "focei_g_half", value = half)
        ))
      } else {
        G[j, k] <- (as.numeric(.nm_focei_pick_f(fp, j)) -
          as.numeric(.nm_focei_pick_f(fm, j))) / half
      }
    }
  }
  G
}

#' @keywords internal
.nm_focei_subject_G_sens <- function(model, subj, theta, omega, eta, sigma,
                                     pk_engine, method = NULL) {
  eta <- as.numeric(eta)
  n_eta <- length(eta)
  pk_engine <- .nm_resolve_pk_engine(pk_engine, model, theta, omega, sigma, eta)
  pred0 <- .nm_subject_ipred(
    model, subj, theta, omega, eta, sigma, pk_engine = pk_engine
  )
  f <- pred0$F
  dv <- pred0$subj_ev$DV[pred0$obs_idx]
  n_obs <- .nm_focei_n_obs(f)
  if (n_obs == 0L || n_eta == 0L) {
    return(list(F = f, dv = dv, G = matrix(0, n_obs, n_eta)))
  }
  gmeth <- .nm_focei_pick_g_method(model, subj, pk_engine, method)
  G <- NULL
  if (identical(gmeth, "ad")) {
    G <- .nm_focei_subject_G_ad(
      model, subj, theta, omega, eta, sigma, pk_engine, f, n_obs, n_eta
    )
    if (is.null(G)) {
      gmeth <- "shi"
    }
  }
  if (is.null(G)) {
    use_shi <- identical(gmeth, "shi")
    G <- .nm_focei_subject_G_fd(
      model, subj, theta, omega, eta, sigma, pk_engine,
      f, n_obs, n_eta, use_shi = use_shi
    )
  }
  list(F = f, dv = dv, G = G, g_method = gmeth)
}
