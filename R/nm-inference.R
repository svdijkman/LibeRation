#' @keywords internal
.nm_fit_laplace_gh <- function(fit, model) {
  n_quad <- fit$n_quad %||% fit$n_quad_requested %||% 5L
  .nm_gh_nodes(.nm_effective_n_quad(model, as.integer(n_quad)))
}

#' @keywords internal
.nm_seed_laplace_eta_modes <- function(fit, data = fit$data) {
  if (identical(fit$method, "LAPLACE") &&
      !is.null(fit$eta) &&
      is.matrix(fit$eta) &&
      nrow(fit$eta) > 0L) {
    p <- .nm_unpack(fit$model, fit$par)
    .nm_laplace_eta_modes_set(
      fit$model, data, fit$eta,
      p$theta, p$omega, p$sigma
    )
  }
  invisible(NULL)
}

#' @keywords internal
.nm_fit_laplace_objective <- function(fit, par, data = fit$data) {
  model <- fit$model
  par <- .nm_apply_fix(model, par)
  p <- .nm_unpack(model, par)
  gh <- .nm_fit_laplace_gh(fit, model)
  pk <- fit$pk_engine %||% "auto"
  mode_centered <- fit$laplace_mode_centered %||% .nm_laplace_mode_centered()
  if (.nm_cpp_capable(model)) {
    pk_eff <- .nm_resolve_pk_engine(pk, model, p$theta, p$omega, p$sigma, NULL)
    if (pk_eff == "cpp") {
      return(.nm_laplace_nll_cpp(
        model, data, p$theta, p$omega, p$sigma, gh,
        mode_centered = mode_centered
      ))
    }
  }
  .nm_laplace_nll_internal(
    model, data, p$theta, p$omega, p$sigma, gh, pk_engine = pk
  )
}

#' Evaluate the reported objective at a parameter vector
#'
#' Matches the likelihood definition used for each estimation method at fixed
#' post-hoc ETAs (where applicable).
#'
#' @param fit An \code{nm_fit} object.
#' @param par Packed parameter vector.
#' @param data Optional dataset; defaults to \code{fit$data}.
#' @keywords internal
.nm_fit_inference_objective <- function(fit, par, data = fit$data) {
  model <- fit$model
  par <- .nm_apply_fix(model, par)
  pk <- fit$pk_engine %||% "auto"
  meth <- fit$method
  eta_mat <- fit$eta
  if (identical(meth, "FO")) {
    return(.nm_fo_report_objective(model, data, par, eta_mat, pk_engine = pk))
  }
  if (identical(meth, "LAPLACE")) {
    return(.nm_fit_laplace_objective(fit, par, data = data))
  }
  p <- .nm_unpack(model, par)
  if (identical(meth, "FOCEI")) {
    return(.nm_focei_objective(
      model, data, p$theta, p$omega, p$sigma, eta_mat, pk_engine = pk
    ))
  }
  .nm_nll_internal(
    model, data, p$theta, p$omega, p$sigma,
    eta = eta_mat, include_omega_prior = TRUE, pk_engine = pk
  )
}

#' Build a tape-AD objective matching \code{.nm_fit_inference_objective}.
#' @keywords internal
.nm_build_inference_objective <- function(fit, data = fit$data) {
  model <- fit$model
  eta_mat <- fit$eta
  pk <- fit$pk_engine %||% "auto"
  meth <- fit$method
  if (identical(meth, "FO")) {
    th_n <- paste0("THETA", model$THETAS$THETA)
    om_n <- paste0("OMEGA", model$OMEGAS$OMEGA)
    sg_n <- paste0("SIGMA", model$SIGMAS$SIGMA)
    all_n <- c(th_n, om_n, sg_n)
    src <- paste0(
      "function(", paste(all_n, collapse = ", "), ") ",
      ".nm_fo_report_objective(model, data, c(",
      paste(all_n, collapse = ", "), "), eta_mat, pk_engine = pk)"
    )
    fn <- eval(parse(text = src))
    environment(fn) <- list2env(
      list(
        model = model,
        data = data,
        eta_mat = eta_mat,
        pk = pk,
        .nm_fo_report_objective = .nm_fo_report_objective
      ),
      parent = parent.env(environment())
    )
    return(list(fn = fn))
  }
  if (identical(meth, "LAPLACE")) {
    fit_env <- new.env(parent = emptyenv())
    fit_env$fit <- fit
    fit_env$data <- data
    fit_env$.nm_fit_laplace_objective <- .nm_fit_laplace_objective
    th_n <- paste0("THETA", model$THETAS$THETA)
    om_n <- paste0("OMEGA", model$OMEGAS$OMEGA)
    sg_n <- paste0("SIGMA", model$SIGMAS$SIGMA)
    all_n <- c(th_n, om_n, sg_n)
    src <- paste0(
      "function(", paste(all_n, collapse = ", "), ") ",
      ".nm_fit_laplace_objective(fit, c(",
      paste(all_n, collapse = ", "), "), data = data)"
    )
    fn <- eval(parse(text = src))
    environment(fn) <- fit_env
    return(list(fn = fn))
  }
  if (identical(meth, "FOCEI")) {
    return(.nm_build_focei_objective(
      model, data, eta_mat, pk_engine = pk
    ))
  }
  .nm_build_pop_objective(
    model, data, eta_mat, include_omega_prior = TRUE, pk_engine = pk
  )
}

#' @keywords internal
.nm_fit_hessian_at_fit <- function(fit, data = fit$data, hessian = c("auto", "ad", "numeric")) {
  hessian <- match.arg(hessian)
  model <- fit$model
  if (is.null(model) || is.null(fit$par)) {
    return(NULL)
  }
  free <- which(!.nm_fix_mask(model))
  if (length(free) == 0L) {
    return(NULL)
  }
  fn <- function(x) {
    par <- fit$par
    par[free] <- x
    .nm_fit_inference_objective(fit, .nm_apply_fix(model, par), data = data)
  }
  x0 <- fit$par[free]
  H <- NULL
  .nm_job_progress_event(
    "cov_hessian_compute",
    list(method = fit$method, hessian = hessian),
    log_msg = paste0("Covariance: computing Hessian (", hessian, ")")
  )
  if (!identical(hessian, "numeric") &&
      !identical(fit$method, "BAYES")) {
    H <- .nm_fit_ad_hessian(fit, data = data, free = free)
  }
  if (is.null(H)) {
    H <- tryCatch(
      .nm_num_hessian(fn, x0),
      error = function(e) NULL
    )
  }
  if (is.null(H)) {
    return(NULL)
  }
  list(hessian = H, free = free)
}

#' @keywords internal
.nm_fit_ad_hessian <- function(fit, data = fit$data, free = NULL) {
  if (!requireNamespace("LibeRtAD", quietly = TRUE)) {
    return(NULL)
  }
  if (!exists("autodiff_hessian", where = asNamespace("LibeRtAD"), inherits = FALSE)) {
    return(NULL)
  }
  if (identical(fit$method, "BAYES")) {
    return(NULL)
  }
  model <- fit$model
  pk <- fit$pk_engine %||% "auto"
  backend <- fit$grad_backend %||% "cpp"
  if (identical(backend, "cpp") && !.nm_cpp_capable(model)) {
    backend <- "R"
  }
  infer_ad <- tryCatch(
    .nm_build_inference_objective(fit, data = data),
    error = function(e) NULL
  )
  if (is.null(infer_ad) || is.null(infer_ad$fn)) {
    return(NULL)
  }
  labels <- .nm_par_labels(model)
  if (is.null(free)) {
    free <- which(!.nm_fix_mask(model))
  }
  if (length(free) == 0L) {
    return(NULL)
  }
  at <- stats::setNames(as.list(fit$par), labels)
  Hres <- tryCatch(
    do.call(
      LibeRtAD::autodiff_hessian,
      c(list(f = infer_ad$fn), at, list(backend = backend))
    ),
    error = function(e) NULL
  )
  if (is.null(Hres) || is.null(Hres$hessian)) {
    return(NULL)
  }
  H <- Hres$hessian
  if (length(free) < length(labels)) {
    H <- H[free, free, drop = FALSE]
  }
  H
}

#' Posterior SD and quantile intervals from a BAYES fit MCMC chains.
#'
#' @param fit An \code{nm_fit} from \code{method = "BAYES"}.
#' @return List with \code{se}, \code{ci_low}, \code{ci_high} (named like \code{fit$par}).
#' @keywords internal
.nm_bayes_posterior_intervals <- function(fit) {
  ch <- fit$chains
  if (is.null(ch) || is.null(ch$theta) || nrow(ch$theta) < 2L) {
    return(NULL)
  }
  model <- fit$model
  labels <- .nm_par_labels(model)
  mat <- cbind(
    ch$theta,
    if (!is.null(ch$omega) && is.matrix(ch$omega) && ncol(ch$omega) > 0L) ch$omega,
    if (!is.null(ch$sigma) && is.matrix(ch$sigma) && ncol(ch$sigma) > 0L) ch$sigma
  )
  if (is.null(mat) || ncol(mat) != length(labels)) {
    return(NULL)
  }
  colnames(mat) <- labels
  se <- apply(mat, 2, stats::sd)
  ci_lo <- apply(mat, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
  ci_hi <- apply(mat, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
  list(
    se = stats::setNames(as.numeric(se), labels),
    ci_low = stats::setNames(as.numeric(ci_lo), labels),
    ci_high = stats::setNames(as.numeric(ci_hi), labels)
  )
}

#' @keywords internal
.nm_fit_attach_bayes_inference <- function(fit) {
  if (is.null(fit) || !identical(fit$method, "BAYES")) {
    return(fit)
  }
  iv <- .nm_bayes_posterior_intervals(fit)
  if (is.null(iv)) {
    fit$inference_note <- "Posterior intervals unavailable (no MCMC chain)."
    return(fit)
  }
  fit$par_se <- iv$se
  fit$par_ci_low <- iv$ci_low
  fit$par_ci_high <- iv$ci_high
  fit$inference_method <- "posterior"
  fit
}

#' Attach post-hoc gradients and standard errors to a fit (called once after estimation).
#' @keywords internal
.nm_fit_attach_inference <- function(fit,
                                     data = fit$data,
                                     hessian = c("auto", "ad", "numeric"),
                                     compute_grad = TRUE,
                                     cov_method = c("auto", "hessian", "linfim", "sandwich"),
                                     compute_covariance = TRUE,
                                     refit_eta = TRUE,
                                     control = list()) {
  if (is.null(fit)) {
    return(fit)
  }
  if (identical(fit$method, "BAYES")) {
    return(.nm_fit_attach_bayes_inference(fit))
  }
  hessian <- match.arg(hessian)
  cov_method <- match.arg(cov_method)
  .nm_seed_laplace_eta_modes(fit, data = data)
  if (isTRUE(compute_grad) && is.null(fit$par_grad)) {
    fit$par_grad <- tryCatch(
      nm_fit_par_gradients(fit, data = data),
      error = function(e) NULL
    )
  }
  if (isTRUE(compute_covariance)) {
    fit <- nm_cov_step(
      fit,
      data = data,
      method = cov_method,
      hessian = hessian,
      refit_eta = refit_eta,
      update_se = TRUE,
      control = control
    )
    if (identical(fit$method, "LAPLACE") && !is.null(fit$par_se)) {
      n_zero <- sum(is.finite(fit$par_se) & fit$par_se == 0, na.rm = TRUE)
      if (n_zero > 0L) {
        fit$inference_note <- paste0(
          "Some SEs are 0: the Laplace objective Hessian was singular or flat ",
          "for those parameters (common for highly correlated THETAs such as VC/VP). ",
          "Try more data, fewer ETAs, or options(LibeRtAD.laplace_opt_reset = TRUE) ",
          "to reset internal Laplace modes on a failed optimizer step."
        )
      }
    }
    return(fit)
  }
  if (is.null(fit$par_se)) {
    fit$par_se <- tryCatch(
      nm_fit_standard_errors(fit, data = data, hessian = hessian),
      error = function(e) NULL
    )
    if (!is.null(fit$par_se) && !any(is.finite(unname(fit$par_se)))) {
      fit$inference_note <- "Standard errors could not be computed (Hessian singular or non-finite)."
    }
  }
  fit
}

#' Population-parameter gradients at the fit optimum
#'
#' Uses C++ adjoint gradients when supported, otherwise LibeRtAD tape AD.
#' Numeric fallback is used only when AD is unavailable.
#'
#' @param fit An \code{nm_fit} object.
#' @param data Optional dataset; defaults to \code{fit$data}.
#' @return Named numeric vector aligned with \code{fit$par}, or \code{NULL}.
#' @examples
#' \dontrun{
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "FO",
#'               control = list(maxit = 5L, compute_inference = FALSE))
#' nm_fit_par_gradients(fit)
#' }
#' @export
nm_fit_par_gradients <- function(fit, data = fit$data) {
  if (is.null(fit) || is.null(fit$par) || is.null(fit$model)) {
    return(NULL)
  }
  if (!is.null(fit$par_grad) && is.numeric(fit$par_grad)) {
    return(fit$par_grad)
  }
  model <- fit$model
  labels <- .nm_par_labels(model)
  par <- fit$par
  eta_mat <- fit$eta
  pk <- fit$pk_engine %||% "auto"
  grad_req <- fit$grad_requested %||% fit$grad %||% "auto"
  p <- .nm_unpack(model, par)

  if (identical(fit$method, "FO")) {
    fn <- function(p) {
      .nm_fo_report_objective(model, data, .nm_apply_fix(model, p), eta_mat, pk_engine = pk)
    }
    g <- .nm_num_grad(fn, par)
    return(stats::setNames(g, labels))
  }

  if (identical(fit$method, "FOCEI")) {
    if (.nm_grad_uses_ad(grad_req)) {
      focei_ad <- .nm_build_focei_objective(
        model, data, eta_mat, pk_engine = pk
      )
      at <- stats::setNames(as.list(par), labels)
      tape_key <- .nm_ad_tape_key(model, data, eta_mat, "focei")
      g <- .nm_ad_eval_cached(
        focei_ad$fn, at, labels, fit$grad_backend %||% "cpp",
        need_grad = TRUE, tape_key = tape_key
      )
      return(stats::setNames(unname(g[labels]), labels))
    }
    fn <- function(p) {
      .nm_fit_inference_objective(fit, p, data = data)
    }
    return(stats::setNames(.nm_num_grad(fn, par), labels))
  }

  if (identical(fit$method, "LAPLACE")) {
    gh <- .nm_fit_laplace_gh(fit, model)
    lap_grad <- fit$grad_effective %||% fit$grad %||% "auto"
    if (identical(lap_grad, "cpp") || identical(lap_grad, "ad")) {
      gres <- .nm_laplace_nll_grad_cpp(
        model, data, p$theta, p$omega, p$sigma, gh,
        mode_centered = fit$laplace_mode_centered %||% .nm_laplace_mode_centered()
      )
      return(stats::setNames(gres$gradient, labels))
    }
    fn <- function(p) {
      .nm_fit_laplace_objective(fit, p, data = data)
    }
    return(stats::setNames(.nm_num_grad(fn, par), labels))
  }

  if (.nm_use_cpp_pop_grad(model, grad_req)) {
    gres <- .nm_pop_nll_grad_cpp(
      model, data, p$theta, p$omega, p$sigma, eta_mat,
      include_omega_prior = TRUE
    )
    return(stats::setNames(gres$gradient, labels))
  }

  if (.nm_grad_uses_ad(grad_req)) {
    pop_ad <- .nm_build_pop_objective(
      model, data, eta_mat, include_omega_prior = TRUE, pk_engine = pk
    )
    backend <- fit$grad_backend %||% "cpp"
    at <- stats::setNames(as.list(par), labels)
    tape_key <- .nm_ad_tape_key(model, data, eta_mat, "pop")
    g <- .nm_ad_eval_cached(
      pop_ad$fn, at, labels, backend, need_grad = TRUE, tape_key = tape_key
    )
    return(stats::setNames(unname(g[labels]), labels))
  }

  fn <- function(p) {
    .nm_fit_inference_objective(fit, p, data = data)
  }
  stats::setNames(.nm_num_grad(fn, par), labels)
}

#' Standard errors from the Hessian of the objective at the fit optimum
#'
#' Computed once at the final parameter estimates (not during optimization).
#' Prefers LibeRtAD \code{autodiff_hessian} on the method-appropriate inference
#' objective when available; otherwise uses an internal numeric Hessian for
#' all applicable methods (including FO and FOCEI). BAYES uses posterior SD from MCMC.
#'
#' @param fit An \code{nm_fit} object.
#' @param data Optional dataset; defaults to \code{fit$data}.
#' @param pk_engine Ignored (kept for compatibility); uses \code{fit$pk_engine}.
#' @param hessian \code{"auto"} tries AD then numeric; \code{"ad"} or \code{"numeric"}.
#' @param method Covariance method: \code{"hessian"}, \code{"linfim"}, or \code{"sandwich"}.
#' @return Named numeric vector aligned with \code{fit$par}.
#' @examples
#' \dontrun{
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "FO",
#'               control = list(maxit = 5L))
#' nm_fit_standard_errors(fit)
#' }
#' @export
nm_fit_standard_errors <- function(fit,
                                   data = fit$data,
                                   pk_engine = NULL,
                                   hessian = c("auto", "ad", "numeric"),
                                   method = c("hessian", "linfim", "sandwich")) {
  method <- match.arg(method)
  if (is.null(fit) || is.null(fit$par)) {
    return(NULL)
  }
  if (identical(fit$method, "BAYES")) {
    iv <- .nm_bayes_posterior_intervals(fit)
    if (is.null(iv)) {
      return(stats::setNames(rep(NA_real_, length(fit$par)), .nm_par_labels(fit$model)))
    }
    return(iv$se)
  }
  if (!identical(method, "hessian") && !is.null(fit$covariance)) {
    vc <- fit$covariance[[method]]
    if (!is.null(vc)) {
      labels <- .nm_par_labels(fit$model)
      se <- sqrt(pmax(diag(vc), 0))
      return(stats::setNames(se, labels))
    }
  }
  if (!is.null(fit$par_se) && identical(method, "hessian")) {
    return(fit$par_se)
  }
  model <- fit$model
  labels <- .nm_par_labels(model)
  if (!identical(method, "hessian")) {
    vc <- nm_fit_covariance(fit, data = data, type = method, hessian = hessian)
    if (!is.null(vc)) {
      se <- sqrt(pmax(diag(vc), 0))
      return(stats::setNames(se, labels))
    }
    return(stats::setNames(rep(NA_real_, length(labels)), labels))
  }
  hres <- .nm_fit_hessian_at_fit(fit, data = data, hessian = hessian)
  if (is.null(hres) || is.null(hres$hessian)) {
    return(stats::setNames(rep(NA_real_, length(fit$par)), labels))
  }
  H <- (hres$hessian + t(hres$hessian)) / 2
  free <- hres$free
  vcov <- .nm_stable_vcov(H)
  se <- rep(NA_real_, length(labels))
  if (!is.null(vcov)) {
    se[free] <- sqrt(pmax(diag(vcov), 0))
  }
  stats::setNames(se, labels)
}

#' Parameter table with estimates, SE, and confidence intervals
#'
#' @param model An \code{nm_model} object.
#' @param fit An \code{nm_fit} object.
#' @param level Confidence level (two-sided).
#' @param compute_se If \code{TRUE}, compute SE when \code{fit$par_se} is missing
#'   (normally SE is attached at the end of \code{\link{nm_est}}).
#' @examples
#' \dontrun{
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "FO",
#'               control = list(maxit = 5L))
#' nm_fit_param_table(sim$model, fit)
#' }
#' @export
nm_fit_param_table <- function(model, fit = NULL, level = 0.95, compute_se = FALSE) {
  if (is.null(fit)) {
    return(nm_workspace_param_table(model, fit = NULL))
  }
  labels <- .nm_par_labels(model)
  se <- if (isTRUE(compute_se) || !is.null(fit$par_se)) {
    if (!is.null(fit$par_se)) {
      fit$par_se
    } else {
      nm_fit_standard_errors(fit, data = fit$data)
    }
  } else {
    NULL
  }
  grad_vec <- if (!is.null(fit$par_grad)) {
    fit$par_grad
  } else if (isTRUE(compute_se)) {
    tryCatch(nm_fit_par_gradients(fit, data = fit$data), error = function(e) NULL)
  } else {
    NULL
  }
  z <- stats::qnorm(1 - (1 - level) / 2)
  rows <- list()
  for (i in seq_along(labels)) {
    lbl <- labels[[i]]
    type <- if (grepl("^THETA", lbl)) {
      "THETA"
    } else if (grepl("^OMEGA", lbl)) {
      "OMEGA"
    } else {
      "SIGMA"
    }
    est <- fit$par[[i]]
    se_i <- if (!is.null(se)) se[[lbl]] else NA_real_
    if (!is.null(fit$par_ci_low) && !is.null(fit$par_ci_high) &&
        lbl %in% names(fit$par_ci_low) && lbl %in% names(fit$par_ci_high)) {
      ci_lo <- fit$par_ci_low[[lbl]]
      ci_hi <- fit$par_ci_high[[lbl]]
    } else {
      ci_lo <- if (is.finite(se_i)) est - z * se_i else NA_real_
      ci_hi <- if (is.finite(se_i)) est + z * se_i else NA_real_
    }
    g_i <- if (!is.null(grad_vec) && lbl %in% names(grad_vec)) {
      grad_vec[[lbl]]
    } else {
      NA_real_
    }
    rows[[length(rows) + 1L]] <- data.frame(
      type = type,
      name = lbl,
      estimate = est,
      se = se_i,
      ci_low = ci_lo,
      ci_high = ci_hi,
      gradient = g_i,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0L) {
    return(data.frame(
      type = character(),
      name = character(),
      estimate = numeric(),
      se = numeric(),
      ci_low = numeric(),
      ci_high = numeric(),
      gradient = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' @rdname nm_fit_param_table
#' @export
nm_workspace_param_table <- function(model, fit = NULL, level = 0.95, compute_se = FALSE) {
  if (is.null(fit)) {
    rows <- list()
    for (i in seq_len(nrow(model$THETAS))) {
      rows[[length(rows) + 1L]] <- data.frame(
        type = "THETA",
        name = paste0("THETA", model$THETAS$THETA[i]),
        initial = as.numeric(model$THETAS$Value[i]),
        estimate = NA_real_,
        se = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        gradient = NA_real_,
        stringsAsFactors = FALSE
      )
    }
    for (i in seq_len(nrow(model$OMEGAS))) {
      rows[[length(rows) + 1L]] <- data.frame(
        type = "OMEGA",
        name = paste0("OMEGA", model$OMEGAS$OMEGA[i]),
        initial = as.numeric(model$OMEGAS$Value[i]),
        estimate = NA_real_,
        se = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        gradient = NA_real_,
        stringsAsFactors = FALSE
      )
    }
    for (i in seq_len(nrow(model$SIGMAS))) {
      rows[[length(rows) + 1L]] <- data.frame(
        type = "SIGMA",
        name = paste0("SIGMA", model$SIGMAS$SIGMA[i]),
        initial = as.numeric(model$SIGMAS$Value[i]),
        estimate = NA_real_,
        se = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        gradient = NA_real_,
        stringsAsFactors = FALSE
      )
    }
    if (length(rows) == 0L) {
      return(data.frame(
        type = character(), name = character(), initial = numeric(),
        estimate = numeric(), se = numeric(), ci_low = numeric(),
        ci_high = numeric(), gradient = numeric(),
        stringsAsFactors = FALSE
      ))
    }
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    return(out)
  }
  ext <- nm_fit_param_table(model, fit, level = level, compute_se = compute_se)
  init <- nm_workspace_param_table(model, fit = NULL)
  init_map <- stats::setNames(init$initial, init$name)
  ext$initial <- unname(init_map[ext$name])
  ext[, c("type", "name", "initial", "estimate", "se", "ci_low", "ci_high", "gradient")]
}
