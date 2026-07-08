#' @keywords internal
.nm_is_ad <- function(x) {
  inherits(x, "Variable") || inherits(x, "Constant")
}

#' @keywords internal
.nm_any_ad <- function(...) {
  xs <- list(...)
  any(vapply(xs, function(x) {
    if (is.list(x) && !.nm_is_ad(x)) {
      any(vapply(x, .nm_is_ad, logical(1)))
    } else if (is.matrix(x)) {
      any(vapply(c(x), .nm_is_ad, logical(1)))
    } else if (is.vector(x) && length(x) > 1L && !.nm_is_ad(x)) {
      any(vapply(x, .nm_is_ad, logical(1)))
    } else {
      .nm_is_ad(x)
    }
  }, logical(1)))
}

#' @keywords internal
.nm_par_labels <- function(model) {
  c(
    paste0("THETA", model$THETAS$THETA),
    paste0("OMEGA", model$OMEGAS$OMEGA),
    paste0("SIGMA", model$SIGMAS$SIGMA)
  )
}

#' @keywords internal
.nm_eta_labels <- function(model) {
  if (.nm_n_eta(model) == 0L) {
    return(character())
  }
  paste0("ETA", seq_len(.nm_n_eta(model)))
}

#' @keywords internal
.nm_split_par_named <- function(model, args) {
  th_n <- paste0("THETA", model$THETAS$THETA)
  om_n <- paste0("OMEGA", model$OMEGAS$OMEGA)
  sg_n <- paste0("SIGMA", model$SIGMAS$SIGMA)
  list(
    theta = unlist(args[th_n], use.names = FALSE),
    omega = unlist(args[om_n], use.names = FALSE),
    sigma = unlist(args[sg_n], use.names = FALSE)
  )
}

#' @keywords internal
.nm_build_pop_objective <- function(model,
                                    data,
                                    eta_mat = NULL,
                                    include_omega_prior = TRUE,
                                    pk_engine = c("auto", "cpp", "R")) {
  th_n <- paste0("THETA", model$THETAS$THETA)
  om_n <- paste0("OMEGA", model$OMEGAS$OMEGA)
  sg_n <- paste0("SIGMA", model$SIGMAS$SIGMA)
  all_n <- c(th_n, om_n, sg_n)
  src <- paste0(
    "function(", paste(all_n, collapse = ", "), ") ",
    ".nm_nll_internal(model, data, ",
    "list(", paste(th_n, collapse = ", "), "), ",
    "list(", paste(om_n, collapse = ", "), "), ",
    "list(", paste(sg_n, collapse = ", "), "), ",
    "eta = eta_mat, include_omega_prior = ",
    include_omega_prior, ", pk_engine = \"", pk_engine, "\")"
  )
  fn <- eval(parse(text = src))
  ctx <- list2env(
    list(
      model = model,
      data = data,
      eta_mat = eta_mat,
      include_omega_prior = include_omega_prior,
      .nm_nll_internal = .nm_nll_internal
    ),
    parent = parent.env(environment())
  )
  environment(fn) <- ctx
  list(fn = fn, ctx = ctx)
}

#' @keywords internal
.nm_build_eta_objective <- function(model,
                                    subj,
                                    theta,
                                    omega,
                                    sigma,
                                    pk_engine = c("auto", "cpp", "R")) {
  n_eta <- .nm_n_eta(model)
  if (n_eta == 0L) {
    return(function() 0)
  }
  eta_n <- paste0("ETA", seq_len(n_eta))
  eta_expr <- paste0("c(", paste(eta_n, collapse = ", "), ")")
  src <- paste0(
    "function(", paste(eta_n, collapse = ", "), ") ",
    ".nm_subject_nll_internal(model, subj, theta, omega, sigma, ",
    "list(", paste(eta_n, collapse = ", "), "), ",
    "include_omega_prior = TRUE, pk_engine = \"", pk_engine, "\")"
  )
  fn <- eval(parse(text = src))
  environment(fn) <- list2env(
    list(
      model = model,
      subj = subj,
      theta = theta,
      omega = omega,
      sigma = sigma,
      .nm_subject_nll_internal = .nm_subject_nll_internal
    ),
    parent = parent.env(environment())
  )
  fn
}

#' @keywords internal
.nm_build_focei_objective <- function(model,
                                      data,
                                      eta_mat = NULL,
                                      pk_engine = c("auto", "cpp", "R")) {
  th_n <- paste0("THETA", model$THETAS$THETA)
  om_n <- paste0("OMEGA", model$OMEGAS$OMEGA)
  sg_n <- paste0("SIGMA", model$SIGMAS$SIGMA)
  all_n <- c(th_n, om_n, sg_n)
  src <- paste0(
    "function(", paste(all_n, collapse = ", "), ") ",
    ".nm_focei_nll_internal(model, data, ",
    "list(", paste(th_n, collapse = ", "), "), ",
    "list(", paste(om_n, collapse = ", "), "), ",
    "list(", paste(sg_n, collapse = ", "), "), ",
    "eta_mat, pk_engine = \"", pk_engine, "\")"
  )
  fn <- eval(parse(text = src))
  ctx <- list2env(
    list(
      model = model,
      data = data,
      eta_mat = eta_mat,
      pk_engine = pk_engine,
      .nm_focei_nll_internal = .nm_focei_nll_internal
    ),
    parent = parent.env(environment())
  )
  environment(fn) <- ctx
  list(fn = fn, ctx = ctx)
}

#' @keywords internal
.nm_build_fo_omega_objective <- function(model, eta_mat) {
  th_n <- paste0("THETA", model$THETAS$THETA)
  om_n <- paste0("OMEGA", model$OMEGAS$OMEGA)
  sg_n <- paste0("SIGMA", model$SIGMAS$SIGMA)
  all_n <- c(th_n, om_n, sg_n)
  src <- paste0(
    "function(", paste(all_n, collapse = ", "), ") ",
    ".nm_fo_omega_prior_nll(model, eta_mat, c(", paste(om_n, collapse = ", "), "))"
  )
  fn <- eval(parse(text = src))
  environment(fn) <- list2env(
    list(
      model = model,
      eta_mat = eta_mat,
      .nm_fo_omega_prior_nll = .nm_fo_omega_prior_nll
    ),
    parent = parent.env(environment())
  )
  fn
}

#' @keywords internal
.ad_max2 <- function(a, b) {
  if (.nm_any_ad(a, b)) .ad_pmax(a, b) else max(a, b)
}

#' @keywords internal
.ad_logsumexp_scalars <- function(terms) {
  if (length(terms) == 0L) {
    return(newConstant(name = "logsumexp_empty", value = -Inf))
  }
  if (!.nm_any_ad(terms)) {
    return(.nm_logsumexp(vapply(terms, as.numeric, numeric(1))))
  }
  if (.ad_use_cpp() && exists("logsumexp_scalars_var", mode = "function")) {
    return(logsumexp_scalars_var(terms))
  }
  m <- terms[[1L]]
  if (length(terms) > 1L) {
    for (i in seq_along(terms)[-1L]) {
      m <- .ad_max2(m, terms[[i]])
    }
  }
  s <- newConstant(name = "logsumexp_zero", value = 0)
  for (i in seq_along(terms)) {
    s <- .ad_add(s, .ad_exp(.ad_sub(terms[[i]], m)))
  }
  .ad_add(m, .ad_log(s))
}

#' @keywords internal
.nm_build_laplace_objective <- function(model,
                                       data,
                                       gh,
                                       pk_engine = c("auto", "cpp", "R")) {
  th_n <- paste0("THETA", model$THETAS$THETA)
  om_n <- paste0("OMEGA", model$OMEGAS$OMEGA)
  sg_n <- paste0("SIGMA", model$SIGMAS$SIGMA)
  all_n <- c(th_n, om_n, sg_n)
  src <- paste0(
    "function(", paste(all_n, collapse = ", "), ") ",
    ".nm_laplace_nll_internal(model, data, ",
    "list(", paste(th_n, collapse = ", "), "), ",
    "list(", paste(om_n, collapse = ", "), "), ",
    "list(", paste(sg_n, collapse = ", "), "), ",
    "gh, pk_engine = \"", pk_engine, "\")"
  )
  fn <- eval(parse(text = src))
  environment(fn) <- list2env(
    list(
      model = model,
      data = data,
      gh = gh,
      pk_engine = pk_engine,
      .nm_laplace_nll_internal = .nm_laplace_nll_internal
    ),
    parent = parent.env(environment())
  )
  fn
}

#' @keywords internal
.ad_exp <- function(x) {
  if (.nm_is_ad(x)) .ad_dispatch("exp", x) else exp(x)
}

#' @keywords internal
.ad_log <- function(x) {
  if (.nm_is_ad(x)) .ad_dispatch("log", x) else log(x)
}

#' @keywords internal
.ad_sub <- function(a, b) {
  if (.nm_any_ad(a, b)) .ad_dispatch("-", a, b) else a - b
}

#' @keywords internal
.ad_mul <- function(a, b) {
  if (.nm_any_ad(a, b)) .ad_dispatch("*", a, b) else a * b
}

#' @keywords internal
.ad_add <- function(a, b) {
  if (.nm_any_ad(a, b)) .ad_dispatch("+", a, b) else a + b
}

#' @keywords internal
.ad_div <- function(a, b) {
  if (.nm_any_ad(a, b)) .ad_dispatch("/", a, b) else a / b
}

#' @keywords internal
.ad_pmax <- function(a, b) {
  if (.nm_any_ad(a, b)) .ad_dispatch("pmax", a, b) else pmax(a, b)
}

#' @keywords internal
.nm_use_cpp_pk_ad <- function() {
  isTRUE(.nm_state$use_cpp_pk)
}

#' @keywords internal
.nm_set_cpp_pk_ad <- function(model, backend) {
  mode <- .nm_cpp_pk_ad_mode(model)
  .nm_state$use_cpp_pk <- isTRUE(
    backend == "cpp" &&
      !is.null(model) &&
      .nm_ad_pk_supported(model) &&
      !is.na(mode)
  )
  .nm_state$cpp_pk_ad_mode <- if (isTRUE(.nm_state$use_cpp_pk)) mode else NA_character_
  invisible(.nm_state$use_cpp_pk)
}

#' @keywords internal
.nm_resolve_laplace_grad <- function(model, grad, gh, n_sub) {
  grad <- match.arg(grad, c("auto", "ad", "numeric", "cpp"))
  if (grad == "cpp") {
    if (!.nm_cpp_capable(model)) {
      .nm_stop("Laplace grad = \"cpp\" requires a C++-capable model.")
    }
    return("cpp")
  }
  if (grad == "numeric") {
    return("numeric")
  }
  if (grad == "ad") {
    if (.nm_cpp_capable(model)) {
      return("cpp")
    }
    return("ad")
  }
  # auto
  if (.nm_cpp_capable(model)) {
    return("cpp")
  }
  mode <- getOption("LibeRtAD.laplace_ad", "auto")
  if (identical(mode, "tape")) {
    return("ad")
  }
  if (identical(mode, "numeric")) {
    return("numeric")
  }
  n_eta <- .nm_n_eta(model)
  n_quad <- length(gh$nodes)
  n_terms <- (n_quad^n_eta) * n_sub
  max_terms <- getOption("LibeRtAD.laplace_ad_max_terms", 500L)
  if (n_terms <= max_terms) "ad" else "numeric"
}

#' @keywords internal
.ad_optim_cache_key <- function(x) {
  paste(format(x, digits = 12, scientific = FALSE), collapse = "\001")
}

#' @keywords internal
.nm_ad_tape_reuse_enabled <- function() {
  isTRUE(getOption("LibeRtAD.tape_reuse", TRUE))
}

#' Stable cache key for AD tape structure (model + fixed ETAs, not par values).
#' @keywords internal
.nm_ad_tape_key <- function(model, data = NULL, eta_mat = NULL, kind = "pop") {
  cfg <- .nm_lik_config(model)
  paste(
    c(
      kind,
      model$ADVAN %||% 0L,
      model$TRANS %||% 0L,
      length(model$THETAS$THETA),
      length(model$OMEGAS$OMEGA),
      length(model$SIGMAS$SIGMA),
      cfg$error_code %||% 0L,
      cfg$omega_code %||% 0L,
      if (!is.null(data)) nrow(data) else 0L,
      if (!is.null(eta_mat) && length(eta_mat) > 0L) {
        .ad_optim_cache_key(eta_mat)
      } else {
        "noeta"
      }
    ),
    collapse = "\003"
  )
}

#' @keywords internal
.nm_ad_tape_cached <- function(tape_key) {
  !is.null(tape_key) &&
    .nm_ad_tape_reuse_enabled() &&
    !is.null(ad_tape_load(tape_key))
}

#' @keywords internal
.nm_ad_run_forward <- function(f, at, backend = "cpp", tape_key = NULL) {
  .ad_check_formals(f, at)
  ad_st <- LibeRtAD:::.ad_state
  reuse <- .nm_ad_tape_cached(tape_key)
  if (reuse) {
    cached <- ad_tape_load(tape_key)
    ad_st$tape <- cached$tape
    ad_st$graph_gen <- cached$graph_gen
    ad_st$active <- FALSE
    set_ops(backend)
    on.exit({
      ad_st$active <- FALSE
      set_ops("R")
    }, add = TRUE)
    parameters <- .ad_bind_tape_parameters(cached$tape, at)
    replay_tape_values_cpp(cached$tape, at)
    result <- cached$tape[[length(cached$tape)]]
    return(list(
      result = result,
      parameters = parameters,
      at = at,
      tape_key = tape_key,
      reused = TRUE
    ))
  }
  reset_tape()
  ad_st$active <- TRUE
  set_ops(backend)
  on.exit({
    ad_st$active <- FALSE
    set_ops("R")
  }, add = TRUE)
  parameters <- .ad_values_to_parameters(at)
  .ad_reset_gradients(parameters)
  env <- .ad_make_eval_env(parameters)
  parent.env(env) <- environment(f)
  args <- lapply(names(at), function(nm) env[[nm]])
  result <- do.call(f, args)
  result <- .ad_as_ad_node(result)
  if (.ad_node_len(result$value) != 1L) {
    .nm_stop("Differentiation requires a scalar NLL output.")
  }
  if (!is.null(tape_key) && .nm_ad_tape_reuse_enabled()) {
    tryCatch(ad_tape_save(tape_key), error = function(e) NULL)
  }
  list(
    result = result,
    parameters = parameters,
    at = at,
    tape_key = tape_key,
    reused = FALSE
  )
}

#' @keywords internal
.nm_ad_run_reverse_only <- function(fwd, backend = "cpp") {
  if (isTRUE(fwd$reused)) {
    reset_tape_grads_cpp(LibeRtAD:::.ad_state$tape)
  }
  .ad_run_reverse(fwd$result, backend)
  .ad_collect_reverse_partials(fwd$parameters)
}

#' @keywords internal
.nm_do_call_autodiff <- function(f, at, backend = "cpp", tape_key = NULL) {
  fwd <- .nm_ad_run_forward(f, at, backend, tape_key = tape_key)
  .nm_ad_run_reverse_only(fwd, backend)
}

#' @keywords internal
.nm_ad_eval_cached <- function(objective_fn, at, par_names, backend = "cpp",
                               need_grad = FALSE, tape_key = NULL) {
  if (is.null(.nm_state$optim_cache)) {
    .nm_state$optim_cache <- list()
  }
  key <- .ad_optim_cache_key(unlist(at[par_names], use.names = FALSE))
  cache <- .nm_state$optim_cache
  if (!is.null(cache) &&
      identical(cache$key, key) &&
      identical(cache$tape_key, tape_key)) {
    if (need_grad) {
      if (is.null(cache$fwd)) {
        cache$fwd <- .nm_ad_run_forward(
          objective_fn, at, backend, tape_key = tape_key
        )
        cache$key <- key
        cache$tape_key <- tape_key
      }
      grad <- .nm_ad_run_reverse_only(cache$fwd, backend)
      return(unname(grad[par_names]))
    }
    if (!is.null(cache$value)) {
      return(cache$value)
    }
  }
  fwd <- .nm_ad_run_forward(objective_fn, at, backend, tape_key = tape_key)
  val <- .ad_scalar_value(fwd$result)
  if (need_grad) {
    grad <- .nm_ad_run_reverse_only(fwd, backend)
    return(unname(grad[par_names]))
  }
  .nm_state$optim_cache <- list(
    key = key, tape_key = tape_key, fwd = fwd, value = val
  )
  val
}

#' @keywords internal
.nm_grad_population <- function(objective_fn, par, par_names, grad, backend) {
  named <- stats::setNames(as.list(par), par_names)
  if (grad == "numeric") {
    return(.nm_num_grad(function(x) {
      names(x) <- par_names
      do.call(objective_fn, as.list(x))
    }, par))
  }
  result <- tryCatch(
    .nm_do_call_autodiff(objective_fn, named, backend),
    error = function(e) {
      if (grad == "ad") {
        .nm_stop("AD gradient failed: ", conditionMessage(e))
      }
      if (requireNamespace("numDeriv", quietly = TRUE)) {
        numDeriv::grad(function(x) {
          names(x) <- par_names
          do.call(objective_fn, as.list(x))
        }, par)
      } else {
        .nm_num_grad(function(x) {
          names(x) <- par_names
          do.call(objective_fn, as.list(x))
        }, par)
      }
    }
  )
  unname(result[par_names])
}

#' @keywords internal
.nm_grad_eta <- function(objective_fn, eta, eta_names, grad, backend) {
  named <- stats::setNames(as.list(eta), eta_names)
  if (grad == "numeric") {
    return(.nm_num_grad(function(x) {
      names(x) <- eta_names
      do.call(objective_fn, as.list(x))
    }, eta))
  }
  result <- tryCatch(
    {
      flat <- .nm_do_call_autodiff(objective_fn, named, backend)
      unname(flat[eta_names])
    },
    error = function(e) {
      if (grad == "ad") {
        .nm_stop("AD gradient failed: ", conditionMessage(e))
      }
      if (requireNamespace("numDeriv", quietly = TRUE)) {
        numDeriv::grad(function(x) {
          names(x) <- eta_names
          do.call(objective_fn, as.list(x))
        }, eta)
      } else {
        .nm_num_grad(function(x) {
          names(x) <- eta_names
          do.call(objective_fn, as.list(x))
        }, eta)
      }
    }
  )
  result
}
