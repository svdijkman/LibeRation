#' One-time FOCEI setup: prepared data, subject cache, likelihood config.
#'
#' Mirrors nlmixr2-style \code{foceiSetup} precomputation for nested FOCEI.
#'
#' @param model An \code{nm_model} object.
#' @param data Dataset used for estimation.
#' @param pk_engine PK engine passed to nested fits.
#' @return An environment with \code{model}, \code{data}, \code{dat}, \code{ids},
#'   optional \code{meta}/\code{subs}, and \code{pk_engine}.
#' @examples
#' sim <- nm_synthetic_theo(n_sub = 2L)
#' nm_focei_setup(sim$model, sim$data)
#' @export
nm_focei_setup <- function(model, data, pk_engine = "auto") {
  pk_engine <- match.arg(pk_engine, c("auto", "cpp", "R"))
  .nm_sync_lik_config(model)
  dat <- .nm_prepare_data(data, model$INPUT, model)
  ids <- .nm_subject_ids(dat)
  meta <- if (.nm_cpp_capable(model)) {
    .nm_cpp_meta(model)
  } else {
    NULL
  }
  subs <- if (!is.null(meta)) {
    .nm_cpp_subjects_cached(model, data)
  } else {
    NULL
  }
  env <- new.env(parent = emptyenv())
  env$model <- model
  env$data <- data
  env$dat <- dat
  env$ids <- ids
  env$meta <- meta
  env$subs <- subs
  env$pk_engine <- pk_engine
  env$eta_mat <- NULL
  env
}

#' @keywords internal
.nm_mceta_init <- function(model, dat, theta, omega, sigma, eta_mat, control = list()) {
  n_eta <- .nm_n_eta(model)
  if (n_eta == 0L) {
    return(NULL)
  }
  ids <- .nm_subject_ids(dat)
  mceta <- control$mceta %||% getOption("LibeRation.mceta", "last")
  mceta <- match.arg(mceta, c("last", "zero", "random"))
  if (identical(mceta, "zero")) {
    return(matrix(0, length(ids), n_eta))
  }
  if (identical(mceta, "random")) {
    return(matrix(rnorm(length(ids) * n_eta), nrow = length(ids), ncol = n_eta))
  }
  if (!is.null(eta_mat) && is.matrix(eta_mat) &&
      nrow(eta_mat) == length(ids) && ncol(eta_mat) == n_eta) {
    return(eta_mat)
  }
  NULL
}

#' Eta mode iteration limits (cold vs warm-started).
#' @keywords internal
.nm_focei_eta_max_iter <- function(control = list(), warm = FALSE) {
  if (isTRUE(warm)) {
    return(as.integer(
      control$maxit_eta_warm %||% getOption("LibeRation.maxit_eta_warm", 40L)
    ))
  }
  as.integer(
    control$maxit_eta %||% control$maxit %||% getOption("LibeRation.maxit_eta", 200L)
  )
}

#' Resolve FOCEI outer gradient method.
#'
#' \code{sensitivity}: one eta refit per \code{par}, then numeric grad of the
#' population objective at fixed eta modes (equals total dL/dtheta when eta is
#' at the conditional mode).
#' @keywords internal
.nm_resolve_focei_est_grad <- function(grad, model) {
  grad <- grad %||% "auto"
  if (identical(grad, "numeric")) {
    return("numeric")
  }
  if (identical(grad, "cpp")) {
    return("numeric")
  }
  if (identical(grad, "sensitivity") || identical(grad, "ad")) {
    if (!.nm_ad_pk_supported(model)) {
      return("numeric")
    }
    return("sensitivity")
  }
  if (identical(grad, "auto")) {
    if (.nm_ad_pk_supported(model) && requireNamespace("LibeRtAD", quietly = TRUE)) {
      return("sensitivity")
    }
    return("numeric")
  }
  "numeric"
}

#' Resolve FOCEI eta refit strategy.
#'
#' \code{nested} (default): re-fit conditional eta modes at each population
#' objective evaluation (NONMEM-style FOCEI).
#' \code{outer}: re-fit eta once per outer iteration and hold it fixed during
#' inner population optimization (faster approximate FOCEI).
#' @keywords internal
.nm_resolve_focei_eta_mode <- function(control = list()) {
  mode <- control$focei_eta %||% getOption("LibeRation.focei_eta", "nested")
  match.arg(mode, c("nested", "outer"))
}

#' FOCEI objective at fixed conditional eta modes (no eta re-fit).
#' @keywords internal
.nm_focei_eval_at_eta <- function(model, data, par, eta_mat, pk_engine) {
  pp <- .nm_unpack(model, .nm_apply_fix(model, par))
  list(
    value = .nm_focei_objective(
      model, data, pp$theta, pp$omega, pp$sigma, eta_mat, pk_engine
    ),
    eta = eta_mat
  )
}

#' Cached nested FOCEI objective + eta modes for fn/gr reuse at the same par.
#' @keywords internal
.nm_focei_nested_fetch <- function(state,
                                    model,
                                    dat,
                                    data,
                                    par,
                                    eta_hint,
                                    pk_engine,
                                    control) {
  key <- .ad_optim_cache_key(par)
  if (identical(state$key, key) && !is.null(state$nested)) {
    return(state$nested)
  }
  nested <- .nm_focei_nested_objective(
    model, dat, data, par, eta_hint, pk_engine, control, use_cache = TRUE
  )
  state$key <- key
  state$nested <- nested
  nested
}

#' Gradient of FOCEI at fixed conditional eta modes (total grad when eta is at mode).
#' @keywords internal
.nm_focei_sensitivity_grad <- function(model, data, par, eta_mat, pk_engine) {
  fn <- function(p) {
    pp <- .nm_unpack(model, .nm_apply_fix(model, p))
    .nm_focei_objective(
      model, data, pp$theta, pp$omega, pp$sigma, eta_mat, pk_engine
    )
  }
  stats::setNames(.nm_num_grad(fn, par), .nm_par_labels(model))
}
