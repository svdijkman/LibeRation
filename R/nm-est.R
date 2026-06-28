#' Estimate NONMEM-style model parameters
#'
#' Supports first-order (FO), first-order conditional estimation (FOCE),
#' stochastic approximation EM (SAEM), Laplace approximation with
#' Gaussian quadrature, and full Bayesian MCMC (METHOD = BAYES).
#'
#' @param model An \code{nm_model} object.
#' @param data An \code{nm_dataset} object.
#' @param method Estimation method: \code{"FO"} (two-step: THETA/SIGMA then OMEGA),
#'   \code{"FOCE"}, \code{"SAEM"}, \code{"LAPLACE"}, or \code{"BAYES"} (full Bayesian MCMC).
#' @param start Optional starting parameter vector (packed THETA, OMEGA, SIGMA).
#' @param backend Automatic differentiation backend for \code{LibeRtAD::autodiff}
#'   when \code{grad} uses AD (\code{"cpp"} or \code{"R"}).
#' @param grad Gradient method: \code{"auto"} (C++ adjoint when supported, else numeric
#'   or tape AD for small grids), \code{"ad"} (tape AD or C++ adjoint when supported),
#'   \code{"cpp"} (C++ quadrature + adjoint gradient), or \code{"numeric"}.
#'   Laplace with tape AD differentiates the quadrature objective in R; prefer
#'   \code{grad = "auto"} or \code{"cpp"} with \code{engine = "cpp"} for speed.
#' @param pk_engine PK solver: \code{"auto"} (C++ when supported and not using AD),
#'   \code{"cpp"}, or \code{"R"}.
#' @param engine Estimation engine for SAEM/Laplace inner loops:
#'   \code{"auto"} (C++ when the model is supported), \code{"cpp"}, or \code{"R"}.
#' @param control List passed to \code{\link[stats]{optim}} plus LibeRation options such as
#'   \code{compute_inference}, \code{cov_method}, \code{cov_refit_eta}, \code{n_cores},
#'   \code{maxit_burn}, \code{factr}, \code{min_retries}, and \code{tweak_inits}.
#' @param ... Method-specific arguments (\code{n_iter}, \code{n_burn}, \code{n_mcmc},
#'   \code{sa_rate}, \code{seed} for SAEM; \code{n_quad} for Laplace; \code{n_burn},
#'   \code{n_sample}, \code{n_thin}, \code{prior}, \code{step_scale} for BAYES, etc.).
#'   For \code{method = "BAYES"}, \code{grad} and \code{backend} are ignored.
#'   Laplace tape AD: \code{options(LibeRtAD.laplace_ad = "auto")} and
#'   \code{options(LibeRtAD.laplace_ad_max_terms = 500)}.
#' @return An \code{nm_fit} object.
#' @examples
#' \dontrun{
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "FO",
#'               control = list(maxit = 5L, compute_inference = FALSE))
#' print(fit)
#' }
#' @export
nm_est <- function(model,
                   data,
                   method = c("FO", "FOCE", "FOCEI", "SAEM", "LAPLACE", "IMP", "BAYES"),
                   start = NULL,
                   backend = c("cpp", "R"),
                   grad = c("auto", "ad", "numeric", "cpp"),
                   pk_engine = c("auto", "cpp", "R"),
                   engine = c("auto", "R", "cpp"),
                   control = list(),
                   ...) {
  method <- match.arg(method)
  grad <- match.arg(grad)
  backend <- match.arg(backend)
  gopts <- .nm_resolve_grad_options(grad, backend)
  grad <- gopts$grad
  backend <- gopts$backend
  pk_engine <- match.arg(pk_engine)
  engine <- match.arg(engine)
  if (!requireNamespace("data.table", quietly = TRUE)) {
    .nm_stop("Package 'data.table' is required for nm_est().")
  }
  par0 <- if (is.null(start)) .nm_init_par(model) else start
  par0 <- .nm_apply_fix(model, par0)
  .nm_sync_lik_config(model)
  nm_validate_model(model, data = data, stop_on_error = TRUE)
  if (isTRUE(getOption("LibeRation.profile", FALSE))) {
    .nm_profile_reset()
  }
  dots <- list(...)
  min_retries <- control$min_retries
  if (is.null(min_retries)) {
    min_retries <- 0L
  }
  min_retries <- as.integer(max(0L, min_retries[1L]))
  tweak_inits <- control$tweak_inits
  if (is.null(tweak_inits)) {
    tweak_inits <- FALSE
  }
  est_fun <- switch(
    method,
    FO = .nm_est_fo,
    FOCE = .nm_est_foce,
    FOCEI = .nm_est_focei,
    SAEM = .nm_est_saem,
    LAPLACE = .nm_est_laplace,
    IMP = .nm_est_imp,
    BAYES = .nm_est_bayes
  )
  attempt_par <- par0
  fit <- NULL
  for (attempt in seq_len(min_retries + 1L)) {
    if (attempt > 1L && !identical(method, "BAYES") &&
        (isTRUE(tweak_inits) || is.numeric(tweak_inits))) {
      attempt_par <- .nm_tweak_inits(model, par0, tweak_inits)
    }
    common <- list(
      model = model, data = data, par0 = attempt_par,
      backend = backend, grad = grad,
      pk_engine = pk_engine, engine = engine,
      control = control
    )
    fit <- do.call(est_fun, c(common, dots))
    if (!.nm_est_needs_retry(fit) || attempt >= min_retries + 1L) {
      break
    }
  }
  if (!is.null(fit) && min_retries > 0L) {
    fit$est_retries <- as.integer(attempt - 1L)
  }
  if (isTRUE(getOption("LibeRation.profile", FALSE)) &&
      is.null(fit$profile) && !identical(method, "BAYES")) {
    fit$profile <- .nm_profile_snapshot()
  }
  if (is.null(fit$engine_detail)) {
    fit$engine_detail <- .nm_engine_detail(
      fit$method, fit$grad, fit$grad_requested, fit$grad_backend,
      fit$pk_engine, fit[["engine"]], .nm_lik_config(model)
    )
  }
  fit$model <- model
  fit$data <- data
  infer <- control$compute_inference
  if (is.null(infer)) {
    infer <- control$compute_covariance
  }
  if (is.null(infer)) {
    infer <- TRUE
  }
  if (isTRUE(infer)) {
    hess <- control$infer_hessian
    if (is.null(hess)) {
      hess <- "numeric"
    }
    cov_method <- control$cov_method
    if (is.null(cov_method)) {
      cov_method <- "auto"
    }
    refit_eta <- control$cov_refit_eta
    if (is.null(refit_eta)) {
      refit_eta <- TRUE
    }
    fit <- .nm_fit_attach_inference(
      fit,
      data = data,
      hessian = hess,
      compute_grad = TRUE,
      cov_method = cov_method,
      compute_covariance = TRUE,
      refit_eta = refit_eta,
      control = control
    )
  }
  fit
}

#' Print method for \code{nm_fit} objects
#'
#' @param x An \code{nm_fit} object.
#' @param ... Unused.
#' @examples
#' \dontrun{
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "FO",
#'               control = list(maxit = 5L, compute_inference = FALSE))
#' print(fit)
#' }
#' @export
print.nm_fit <- function(x, ...) {
  cat("NONMEM-style fit (", x$method, ")\n", sep = "")
  if (identical(x$method, "BAYES")) {
    cat("  log posterior:", format(x$log_posterior, digits = 6), "\n")
    cat("  kept samples:", x$n_keep, " (burn =", x$n_burn, ", thin =", x$n_thin, ")\n")
  } else {
    cat("  objective (-2LL):", format(x$objective, digits = 6), "\n")
  }
  if (!is.null(x$grad) && (length(x$grad) == 0L || !is.na(x$grad)[1L])) {
    cat("  grad:", x$grad)
    if (!is.null(x$grad_requested) && !identical(x$grad_requested, x$grad)) {
      cat("  (requested:", x$grad_requested, ")")
    }
    if (!is.null(x$grad_effective) && !identical(x$grad_effective, x$grad)) {
      cat("  effective:", x$grad_effective)
    }
    if (!is.null(x$grad_backend)) cat("  backend:", x$grad_backend)
    if (!is.null(x$pk_engine)) cat("  pk_engine:", x$pk_engine)
    eng <- x[["engine"]]
    if (!is.null(eng) && (length(eng) == 0L || !is.na(eng)[1L])) {
      cat("  engine:", eng)
    }
    cat("\n")
  }
  cat("  THETA:\n")
  print(stats::setNames(x$theta, paste0("THETA", x$model$THETAS$THETA)))
  if (length(x$omega) > 0L) {
    cat("  OMEGA:\n")
    print(stats::setNames(x$omega, paste0("OMEGA", x$model$OMEGAS$OMEGA)))
  }
  if (length(x$sigma) > 0L) {
    cat("  SIGMA:\n")
    print(stats::setNames(x$sigma, paste0("SIGMA", x$model$SIGMAS$SIGMA)))
  }
  invisible(x)
}

#' Synthetic THEO-like dataset for examples and tests
#'
#' @param n_sub Number of subjects.
#' @param seed Random seed.
#' @return List with \code{model} and \code{data}.
#' @examples
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' names(sim)
#' head(sim$data$data)
#' @export
nm_synthetic_theo <- function(n_sub = 10L, seed = 1L) {
  set.seed(seed)
  thetas <- c(CL = 3.5, VC = 20, VP = 50, Q2 = 10, KA = 1.2)
  omega <- c(0.09, 0.09, 0.09)
  sigma <- c(0.1, 0.5)
  times <- c(0, 0.25, 0.5, 1, 2, 4, 6, 8, 12, 24)
  rows <- list()
  etas <- vector("list", n_sub)
  for (id in seq_len(n_sub)) {
    eta <- rnorm(3, sd = sqrt(omega))
    etas[[id]] <- eta
    ka <- thetas["KA"] * exp(eta[3])
    rows[[length(rows) + 1L]] <- data.frame(
      ID = id, TIME = 0, EVID = 1L, CMT = 1L, AMT = 320, RATE = 0,
      MDV = 1L, DV = 0, F1 = 1, S1 = 1, KA = ka
    )
    for (tm in times[-1]) {
      rows[[length(rows) + 1L]] <- data.frame(
        ID = id, TIME = tm, EVID = 0L, CMT = 2L, AMT = 0, RATE = 0,
        MDV = 0L, DV = NA_real_, F1 = 1, S1 = 1, KA = ka
      )
    }
  }
  dat <- do.call(rbind, rows)
  model <- nm_model(
    INPUT = c("AMT", "TIME", "EVID", "CMT", "MDV", "DV", "F1", "S1", "KA"),
    ADVAN = 4L, TRANS = 4L, DOSECMP = 1L, OBSCMP = 2L,
    PRED = paste(
      "CL = THETA(1) * exp(ETA(1))",
      "VC = THETA(2) * exp(ETA(2))",
      "VP = THETA(3)",
      "Q2 = THETA(4)",
      "KA = THETA(5) * exp(ETA(3))",
      "S1 = 1",
      "S2 = VC",
      "S3 = VP",
      sep = "\n"
    ),
    ERROR = "Y = F * (1 + ERR(1)) + ERR(2)",
    THETAS = data.frame(
      THETA = 1:5,
      Value = c(3, 20, 50, 10, 1),
      stringsAsFactors = FALSE
    ),
    OMEGAS = data.frame(OMEGA = 1:3, Value = rep(0.05, 3)),
    SIGMAS = data.frame(SIGMA = 1:2, Value = c(0.1, 0.5))
  )
  dat <- data.table::as.data.table(dat)
  list(
    model = model,
    data = .nm_synthetic_fill_dv(model, dat, etas, sigma)
  )
}
