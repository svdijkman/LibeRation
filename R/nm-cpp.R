#' @keywords internal
.nm_pred_rhs <- function(line) {
  line <- trimws(line)
  if (grepl("<-", line, fixed = TRUE)) {
    return(trimws(sub(".*<-", "", line)))
  }
  if (grepl("=", line, fixed = TRUE)) {
    return(trimws(sub(".*?=", "", line)))
  }
  line
}

#' @keywords internal
.nm_cpp_pred_supported <- function(model) {
  pred <- .nm_split_lines(model$PRED)
  if (length(pred) == 0L) {
    return(isTRUE(model$USE_ODE))
  }
  if (isTRUE(model$USE_ODE)) {
    pred <- pred[!grepl("A\\s*\\(", pred, ignore.case = TRUE)]
    if (length(pred) == 0L) {
      return(TRUE)
    }
  }
  nm_pred_expr_check_cpp(pred)
}

#' @keywords internal
.nm_n_transit <- function(model) {
  if (grepl("MDTR", model$PRED, fixed = TRUE) || grepl("NMDTR", model$PRED, fixed = TRUE)) {
    m <- regmatches(model$PRED, regexpr("NMDTR\\s*=\\s*[0-9]+", model$PRED))
    if (length(m)) {
      return(as.integer(sub(".*=\\s*", "", m[1])))
    }
  }
  if (grepl("KTR", model$PRED, fixed = TRUE)) {
    return(1L)
  }
  0L
}

#' @keywords internal
.nm_cpp_capable <- function(model) {
  if (!nm_cpp_advan_supported(as.integer(model$ADVAN), as.integer(model$TRANS))) {
    return(FALSE)
  }
  if (isTRUE(model$USE_ODE)) {
    return(TRUE)
  }
  .nm_cpp_pred_supported(model)
}

#' @keywords internal
.nm_resolve_pk_engine <- function(pk_engine, model, theta, omega, sigma, eta = NULL) {
  pk_engine <- match.arg(pk_engine, c("auto", "cpp", "R"))
  if (pk_engine == "R") {
    return("R")
  }
  if (isTRUE(model$USE_ODE) && .nm_cpp_capable(model)) {
    return("cpp")
  }
  if (.nm_any_ad(theta, omega, sigma, eta)) {
    return("R")
  }
  if (pk_engine == "cpp" || (pk_engine == "auto" && .nm_cpp_capable(model))) {
    if (.nm_cpp_capable(model)) {
      return("cpp")
    }
  }
  "R"
}

#' @keywords internal
.nm_data_fingerprint <- function(dat) {
  cols <- intersect(
    c("ID", "TIME", "EVID", "MDV", "DV", "AMT", "RATE", "CMT", "SS", "II", "F1", "S1", "KA"),
    names(dat)
  )
  if (length(cols) == 0L) {
    return("empty")
  }
  payload <- as.data.frame(dat)[, cols, drop = FALSE]
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(payload, algo = "xxhash64"))
  }
  # Short fallback key (must stay under R's 10000-byte symbol limit).
  paste0(
    "d", nrow(dat), "_",
    paste(vapply(cols, function(col) {
      x <- dat[[col]]
      paste0(sum(as.numeric(x), na.rm = TRUE), ":", length(x))
    }, character(1L)), collapse = "-")
  )
}

#' @keywords internal
.nm_data_subject_cache_key <- function(dat) {
  .nm_data_fingerprint(dat)
}

#' @keywords internal
.nm_clear_cpp_subjects_cache <- function() {
  cache <- .nm_env_cache("cpp_subjects")
  rm(list = ls(envir = cache, all.names = TRUE), envir = cache)
}

#' @keywords internal
.nm_cpp_subjects_cached <- function(model, data) {
  dat <- .nm_prepare_data(data, model$INPUT, model)
  ids <- .nm_subject_ids(dat)
  key <- paste(
    .nm_data_subject_cache_key(dat),
    model$ADVAN,
    model$TRANS,
    sep = "_"
  )
  cache <- .nm_env_cache("cpp_subjects")
  if (!is.null(cache[[key]])) {
    return(cache[[key]])
  }
  subs <- .nm_cpp_subjects(model, data)
  cache[[key]] <- subs
  subs
}

#' @keywords internal
.nm_env_cache <- function(name) {
  if (!exists(".nm_cache", envir = .GlobalEnv, inherits = FALSE)) {
    assign(".nm_cache", new.env(parent = emptyenv()), envir = .GlobalEnv)
  }
  cache_root <- get(".nm_cache", envir = .GlobalEnv)
  if (!exists(name, envir = cache_root, inherits = FALSE)) {
    assign(name, new.env(parent = emptyenv()), envir = cache_root)
  }
  get(name, envir = cache_root)
}

#' @keywords internal
.nm_cpp_subjects <- function(model, data) {
  dat <- .nm_prepare_data(data, model$INPUT, model)
  ids <- .nm_subject_ids(dat)
  lapply(ids, function(id) {
    subj <- .nm_subject_slice(dat, id)
    attr(subj, "user_scale_cols") <- attr(dat, "user_scale_cols")
    attr(subj, "user_f_cols") <- attr(dat, "user_f_cols")
    subj_ev <- .nm_subject_events(subj)
    ev <- .nm_input_event_vectors(subj_ev)
    obs_idx <- which(subj_ev$MDV == 0L & subj_ev$EVID == 0L)
    subj_list <- list(
      time = as.numeric(subj_ev$TIME),
      amt = as.numeric(subj_ev$AMT),
      rate = as.numeric(subj_ev$RATE),
      f1 = as.numeric(subj_ev$F1),
      cmt = as.integer(subj_ev$CMT),
      evid = as.integer(subj_ev$EVID),
      ss = as.integer(subj_ev$SS),
      ii = as.numeric(if ("II" %in% names(subj_ev)) subj_ev$II else subj_ev$TAU),
      dv = as.numeric(subj_ev$DV[obs_idx]),
      obs_idx = as.integer(obs_idx),
      s1 = ev$s1,
      s2 = ev$s2,
      s3 = ev$s3,
      s4 = ev$s4,
      scale_mat = ev$scale_mat,
      use_data_scale = ev$use_data_scale,
      f_mat = ev$f_mat,
      use_data_f = ev$use_data_f
    )
    if (!is.null(model$COVARIATES) && length(model$COVARIATES) > 0L) {
      cv <- as.character(model$COVARIATES)
      subj_row <- subj[1L, ]
      cov_vals <- vapply(cv, function(cn) {
        if (cn %in% names(subj_row)) as.numeric(subj_row[[cn]][1L]) else 0
      }, numeric(1))
      subj_list$cov <- stats::setNames(as.list(cov_vals), cv)
    }
    if (!is.null(model$DVID) && "DVID" %in% names(subj_ev)) {
      subj_list$dvid <- as.integer(subj_ev$DVID[obs_idx])
    }
    subj_list
  })
}

#' @keywords internal
.nm_cpp_meta <- function(model) {
  list(
    pred_lines = .nm_split_lines(model$PRED),
    des_lines = .nm_split_lines(model$DES),
    advan = as.integer(model$ADVAN),
    trans = as.integer(model$TRANS),
    obs_cmp = as.integer(model$OBSCMP),
    dose_cmp = as.integer(model$DOSECMP),
    n_transit = .nm_n_transit(model),
    use_ode = isTRUE(model$USE_ODE),
    model_ss = as.integer(model$SS)
  )
}

#' @keywords internal
.nm_pop_nll_cache_key <- function(theta, omega, sigma, eta = NULL) {
  base <- .nm_laplace_cache_key(theta, omega, sigma)
  if (is.null(eta)) {
    return(base)
  }
  em <- as.numeric(eta)
  if (length(em) == 0L) {
    return(base)
  }
  paste(base, .ad_optim_cache_key(em), sep = "\003")
}

#' @keywords internal
.nm_laplace_cache_key <- function(theta, omega, sigma) {
  .ad_optim_cache_key(c(as.numeric(theta), as.numeric(omega), as.numeric(sigma)))
}

#' @keywords internal
.nm_laplace_optim_cache <- function() {
  .nm_env_cache("laplace_optim")
}

#' @keywords internal
.nm_clear_laplace_optim_cache <- function() {
  cache <- .nm_laplace_optim_cache()
  rm(list = ls(envir = cache, all.names = TRUE), envir = cache)
}

#' @keywords internal
.nm_laplace_eta_modes_cache <- function() {
  .nm_env_cache("laplace_eta_modes")
}

#' @keywords internal
.nm_laplace_eta_modes_key <- function(model, data) {
  dat <- .nm_prepare_data(data, model$INPUT, model)
  paste(nrow(dat), length(.nm_subject_ids(dat)), model$ADVAN, model$TRANS, sep = "_")
}

#' @keywords internal
.nm_clear_laplace_eta_modes <- function() {
  cache <- .nm_laplace_eta_modes_cache()
  rm(list = ls(envir = cache, all.names = TRUE), envir = cache)
  invisible(NULL)
}

#' @keywords internal
.nm_laplace_eta_modes_cache_key <- function(model, data, theta, omega, sigma) {
  paste(
    .nm_laplace_eta_modes_key(model, data),
    .nm_laplace_cache_key(theta, omega, sigma),
    sep = "|"
  )
}

#' @keywords internal
.nm_laplace_eta_modes_get <- function(model, data, theta = NULL, omega = NULL, sigma = NULL) {
  cache <- .nm_laplace_eta_modes_cache()
  if (is.null(theta)) {
    return(NULL)
  }
  cache[[.nm_laplace_eta_modes_cache_key(model, data, theta, omega, sigma)]]
}

#' @keywords internal
.nm_laplace_eta_modes_set <- function(model, data, eta_modes, theta, omega, sigma) {
  if (is.null(eta_modes) || !is.matrix(eta_modes) || nrow(eta_modes) == 0L) {
    return(invisible(NULL))
  }
  cache <- .nm_laplace_eta_modes_cache()
  cache[[.nm_laplace_eta_modes_cache_key(model, data, theta, omega, sigma)]] <- eta_modes
  invisible(NULL)
}

#' @keywords internal
.nm_nll_cpp <- function(model, data, theta, omega, sigma, eta = NULL,
                        include_omega_prior = TRUE) {
  if (!.nm_cpp_capable(model)) {
    .nm_stop("Model is not supported by the C++ likelihood engine.")
  }
  meta <- .nm_cpp_meta(model)
  subs <- .nm_cpp_subjects_cached(model, data)
  ids <- .nm_subject_ids(.nm_prepare_data(data, model$INPUT, model))
  n_eta <- .nm_n_eta(model)
  if (is.null(eta)) {
    eta_mat <- matrix(0, length(ids), n_eta)
  } else {
    eta_mat <- eta
  }
  nm_nll_cpp(
    subjects = subs,
    theta = as.numeric(theta),
    omega = as.numeric(omega),
    sigma = as.numeric(sigma),
    eta = eta_mat,
    pred_lines = meta$pred_lines,
    advan = meta$advan,
    trans = meta$trans,
    obs_cmp = meta$obs_cmp,
    dose_cmp = meta$dose_cmp,
    n_transit = meta$n_transit,
    use_ode = meta$use_ode,
    model_ss = meta$model_ss,
    include_omega_prior = include_omega_prior
  )
}

#' @keywords internal
.nm_pop_optim_cache <- function() {
  .nm_env_cache("pop_optim")
}

#' @keywords internal
.nm_clear_pop_optim_cache <- function() {
  cache <- .nm_pop_optim_cache()
  rm(list = ls(envir = cache, all.names = TRUE), envir = cache)
  if (exists(".nm_clear_focei_nested_cache", mode = "function")) {
    .nm_clear_focei_nested_cache()
  }
}

#' @keywords internal
.nm_use_cpp_pop_grad <- function(model, grad) {
  if (!isTRUE(getOption("LibeRtAD.cpp_pop_grad", TRUE)) || !.nm_cpp_capable(model)) {
    return(FALSE)
  }
  if (identical(grad, "cpp") || .nm_grad_uses_ad(grad)) {
    return(TRUE)
  }
  isTRUE(getOption("LibeRation.cpp_pop_grad_numeric", TRUE)) && identical(grad, "numeric")
}

#' @keywords internal
.nm_fit_all_eta_cpp <- function(model, data, theta, omega, sigma, eta_mat = NULL,
                                max_iter = 200L) {
  meta <- .nm_cpp_meta(model)
  subs <- .nm_cpp_subjects_cached(model, data)
  eta_init <- if (is.null(eta_mat)) NULL else eta_mat
  nm_fit_all_eta_cpp(
    subjects = subs,
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
    model_ss = meta$model_ss,
    eta_init = eta_init,
    max_iter = as.integer(max_iter)
  )
}

#' @keywords internal
.nm_pop_nll_cpp <- function(model, data, theta, omega, sigma, eta = NULL,
                           include_omega_prior = TRUE, cache_fwd = FALSE) {
  meta <- .nm_cpp_meta(model)
  subs <- .nm_cpp_subjects_cached(model, data)
  ids <- .nm_subject_ids(.nm_prepare_data(data, model$INPUT, model))
  n_eta <- .nm_n_eta(model)
  if (is.null(eta)) {
    eta_mat <- matrix(0, length(ids), n_eta)
  } else {
    eta_mat <- eta
  }
  cache <- .nm_pop_optim_cache()
  cache_key <- .nm_pop_nll_cache_key(theta, omega, sigma, eta_mat)
  if (isTRUE(cache_fwd) && !is.null(cache$key) && identical(cache$key, cache_key)) {
    return(cache$objective)
  }
  res <- nm_nll_detailed_cpp(
    subjects = subs,
    theta = as.numeric(theta),
    omega = as.numeric(omega),
    sigma = as.numeric(sigma),
    eta = eta_mat,
    pred_lines = meta$pred_lines,
    advan = meta$advan,
    trans = meta$trans,
    obs_cmp = meta$obs_cmp,
    dose_cmp = meta$dose_cmp,
    n_transit = meta$n_transit,
    use_ode = meta$use_ode,
    model_ss = meta$model_ss,
    include_omega_prior = include_omega_prior
  )
  obj <- res$objective
  if (!is.finite(obj)) {
    obj <- .Machine$double.xmax
  }
  if (isTRUE(cache_fwd)) {
    cache$key <- cache_key
    cache$objective <- obj
    cache$fwd <- res$fwd_subjects
    return(cache$objective)
  }
  obj
}

#' @keywords internal
.nm_pop_nll_grad_cpp <- function(model, data, theta, omega, sigma, eta = NULL,
                                include_omega_prior = TRUE,
                                use_fwd_cache = FALSE) {
  meta <- .nm_cpp_meta(model)
  subs <- .nm_cpp_subjects_cached(model, data)
  ids <- .nm_subject_ids(.nm_prepare_data(data, model$INPUT, model))
  n_eta <- .nm_n_eta(model)
  if (is.null(eta)) {
    eta_mat <- matrix(0, length(ids), n_eta)
  } else {
    eta_mat <- eta
  }
  cache <- .nm_pop_optim_cache()
  cache_key <- .nm_pop_nll_cache_key(theta, omega, sigma, eta_mat)
  fwd_cache <- if (isTRUE(use_fwd_cache) &&
    !is.null(cache$key) &&
    identical(cache$key, cache_key)) {
    cache$fwd
  } else {
    NULL
  }
  names <- c(
    paste0("THETA", model$THETAS$THETA),
    paste0("OMEGA", model$OMEGAS$OMEGA),
    paste0("SIGMA", model$SIGMAS$SIGMA)
  )
  res <- nm_nll_grad_cpp(
    subjects = subs,
    theta = as.numeric(theta),
    omega = as.numeric(omega),
    sigma = as.numeric(sigma),
    eta = eta_mat,
    pred_lines = meta$pred_lines,
    advan = meta$advan,
    trans = meta$trans,
    obs_cmp = meta$obs_cmp,
    dose_cmp = meta$dose_cmp,
    n_transit = meta$n_transit,
    use_ode = meta$use_ode,
    model_ss = meta$model_ss,
    include_omega_prior = include_omega_prior,
    fwd_cache = fwd_cache,
    grad_from_fwd = !is.null(fwd_cache)
  )
  if (isTRUE(use_fwd_cache)) {
    .nm_clear_pop_optim_cache()
  }
  list(
    objective = res$objective,
    gradient = unname(c(res$grad_theta, res$grad_omega, res$grad_sigma)),
    names = names
  )
}

#' @keywords internal
.nm_laplace_mode_centered <- function() {
  isTRUE(getOption("LibeRtAD.laplace_mode_centered", FALSE))
}

#' @keywords internal
.nm_laplace_nll_cpp <- function(model, data, theta, omega, sigma, gh,
                                mode_centered = .nm_laplace_mode_centered(),
                                control = list(), cache_fwd = FALSE) {
  meta <- .nm_cpp_meta(model)
  subs <- .nm_cpp_subjects_cached(model, data)
  eta_modes <- .nm_laplace_eta_modes_get(model, data, theta, omega, sigma)
  cache <- .nm_laplace_optim_cache()
  cache_key <- .nm_laplace_cache_key(theta, omega, sigma)
  if (isTRUE(cache_fwd) && !is.null(cache$key) && identical(cache$key, cache_key)) {
    return(cache$objective)
  }
  res <- nm_laplace_nll_detailed_cpp(
    subjects = subs,
    theta = as.numeric(theta),
    omega = as.numeric(omega),
    sigma = as.numeric(sigma),
    gh_nodes = gh$nodes,
    gh_weights = gh$weights,
    pred_lines = meta$pred_lines,
    advan = meta$advan,
    trans = meta$trans,
    obs_cmp = meta$obs_cmp,
    dose_cmp = meta$dose_cmp,
    n_transit = meta$n_transit,
    use_ode = meta$use_ode,
    model_ss = meta$model_ss,
    mode_centered = isTRUE(mode_centered),
    eta_modes = eta_modes
  )
  .nm_laplace_eta_modes_set(model, data, res$eta_modes, theta, omega, sigma)
  if (isTRUE(cache_fwd)) {
    cache$key <- cache_key
    cache$objective <- res$objective
    cache$fwd <- res$fwd_subjects
  }
  res$objective
}

#' @keywords internal
.nm_laplace_final_eta_modes <- function(model,
                                         data,
                                         theta,
                                         omega,
                                         sigma,
                                         gh,
                                         pk_engine = "auto") {
  dat <- .nm_prepare_data(data, model$INPUT, model)
  ids <- .nm_subject_ids(dat)
  n_eta <- .nm_n_eta(model)
  if (n_eta == 0L) {
    return(matrix(numeric(0), nrow = length(ids), ncol = 0L))
  }
  cached <- .nm_laplace_eta_modes_get(model, data, theta, omega, sigma)
  if (!is.null(cached) && is.matrix(cached) &&
      nrow(cached) == length(ids) && ncol(cached) == n_eta) {
    return(cached)
  }
  if (.nm_cpp_capable(model)) {
    pk_eff <- .nm_resolve_pk_engine(pk_engine, model, theta, omega, sigma, NULL)
    if (pk_eff == "cpp") {
      res <- nm_laplace_nll_detailed_cpp(
        subjects = .nm_cpp_subjects_cached(model, data),
        theta = as.numeric(theta),
        omega = as.numeric(omega),
        sigma = as.numeric(sigma),
        gh_nodes = gh$nodes,
        gh_weights = gh$weights,
        pred_lines = .nm_cpp_meta(model)$pred_lines,
        advan = as.integer(model$ADVAN),
        trans = as.integer(model$TRANS),
        obs_cmp = as.integer(model$OBSCMP),
        dose_cmp = as.integer(model$DOSECMP),
        n_transit = .nm_n_transit(model),
        use_ode = isTRUE(model$USE_ODE),
        model_ss = as.integer(model$SS),
        mode_centered = isTRUE(.nm_laplace_mode_centered()),
        eta_modes = NULL
      )
      if (!is.null(res$eta_modes)) {
        .nm_laplace_eta_modes_set(model, data, res$eta_modes, theta, omega, sigma)
        return(res$eta_modes)
      }
    }
  }
  .nm_fit_all_eta(
    model, dat, theta, omega, sigma, NULL,
    backend = "cpp", grad = "numeric", pk_engine = pk_engine
  )
}

#' @keywords internal
.nm_laplace_nll_grad_cpp <- function(model, data, theta, omega, sigma, gh,
                                     mode_centered = .nm_laplace_mode_centered(),
                                     control = list(), use_fwd_cache = FALSE) {
  meta <- .nm_cpp_meta(model)
  subs <- .nm_cpp_subjects_cached(model, data)
  eta_modes <- .nm_laplace_eta_modes_get(model, data, theta, omega, sigma)
  cache <- .nm_laplace_optim_cache()
  cache_key <- .nm_laplace_cache_key(theta, omega, sigma)
  fwd_cache <- if (isTRUE(use_fwd_cache) &&
    !is.null(cache$key) &&
    identical(cache$key, cache_key)) {
    cache$fwd
  } else {
    NULL
  }
  names <- c(
    paste0("THETA", model$THETAS$THETA),
    paste0("OMEGA", model$OMEGAS$OMEGA),
    paste0("SIGMA", model$SIGMAS$SIGMA)
  )
  res <- nm_laplace_nll_grad_cpp(
    subjects = subs,
    theta = as.numeric(theta),
    omega = as.numeric(omega),
    sigma = as.numeric(sigma),
    gh_nodes = gh$nodes,
    gh_weights = gh$weights,
    pred_lines = meta$pred_lines,
    advan = meta$advan,
    trans = meta$trans,
    obs_cmp = meta$obs_cmp,
    dose_cmp = meta$dose_cmp,
    n_transit = meta$n_transit,
    use_ode = meta$use_ode,
    model_ss = meta$model_ss,
    mode_centered = isTRUE(mode_centered),
    n_threads = 1L,
    eta_modes = eta_modes,
    fwd_cache = fwd_cache,
    grad_from_fwd = !is.null(fwd_cache) && !isTRUE(mode_centered)
  )
  if (is.null(fwd_cache)) {
    .nm_laplace_eta_modes_set(model, data, res$eta_modes, theta, omega, sigma)
  }
  list(
    objective = res$objective,
    gradient = unname(c(res$grad_theta, res$grad_omega, res$grad_sigma)),
    names = names
  )
}

#' @keywords internal
.nm_saem_mh_cpp <- function(model, data, eta_mat, theta, omega, sigma,
                            n_mcmc = 1L, step_scale = 1.0) {
  meta <- .nm_cpp_meta(model)
  subs <- .nm_cpp_subjects_cached(model, data)
  nm_saem_mh_cpp(
    eta = eta_mat,
    subjects = subs,
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
    model_ss = meta$model_ss,
    n_mcmc = as.integer(n_mcmc),
    step_scale = step_scale
  )
}
