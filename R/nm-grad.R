#' @keywords internal
.nm_resolve_grad_options <- function(grad, backend) {
  if (identical(grad, "cpp")) {
    return(list(grad = "cpp", backend = "cpp"))
  }
  list(
    grad = match.arg(grad, c("auto", "ad", "numeric")),
    backend = match.arg(backend, c("cpp", "R"))
  )
}

#' @keywords internal
.nm_grad_uses_ad <- function(grad) {
  grad %in% c("auto", "ad")
}

#' @keywords internal
.nm_ad_pk_supported <- function(model) {
  advan <- as.integer(model$ADVAN)
  trans <- as.integer(model$TRANS)
  if (advan == 1L && trans == 2L) {
    return(TRUE)
  }
  if (advan == 2L && trans %in% c(1L, 2L)) {
    return(TRUE)
  }
  if (advan == 3L && trans == 4L) {
    return(TRUE)
  }
  if (advan == 4L && trans == 4L) {
    return(TRUE)
  }
  FALSE
}

#' @keywords internal
.nm_assert_ad_grad_supported <- function(model, grad, context = "population") {
  if (!identical(grad, "ad")) {
    return(invisible(TRUE))
  }
  if (.nm_ad_pk_supported(model)) {
    return(invisible(TRUE))
  }
  .nm_stop(
    "AD ", context, " gradients are not supported for ADVAN ", model$ADVAN,
    " TRANS ", model$TRANS, " yet. Use grad = \"numeric\" with pk_engine = \"cpp\" ",
    "or grad = \"auto\" (numeric fallback)."
  )
}

#' @keywords internal
.nm_resolve_population_grad <- function(model, grad, objective_fn, par, par_names) {
  if (grad == "numeric") {
    return(list(grad = "numeric", fn = NULL))
  }
  if (!.nm_ad_pk_supported(model)) {
    if (grad == "ad") {
      .nm_stop(
        "AD population gradients are not supported for ADVAN ", model$ADVAN,
        " TRANS ", model$TRANS, ". Use grad = \"numeric\" with pk_engine = \"cpp\"."
      )
    }
    return(list(grad = "numeric", fn = NULL))
  }
  list(grad = grad, fn = NULL)
}

#' @keywords internal
.nm_cpp_pk_ad_mode <- function(model) {
  if (!isTRUE(getOption("LibeRtAD.cpp_pk_ad", TRUE)) || is.null(model)) {
    return(NA_character_)
  }
  advan <- as.integer(model$ADVAN)
  trans <- as.integer(model$TRANS)
  if (advan == 1L && trans == 2L) {
    return("bolus1")
  }
  if (advan == 2L && trans %in% c(1L, 2L)) {
    return("oral1")
  }
  if (advan == 3L && trans == 4L) {
    return("oral2_trans4")
  }
  if (advan == 4L && trans == 4L) {
    return("oral2_trans4")
  }
  NA_character_
}

#' @keywords internal
.nm_effective_pk_engine <- function(pk_engine, grad) {
  match.arg(pk_engine, c("auto", "cpp", "R"))
}

.nm_resolve_estimation_grad <- function(model, grad) {
  if (identical(grad, "cpp")) {
    if (!.nm_cpp_capable(model)) {
      .nm_stop(
        "C++ population gradients require a C++-capable model (ADVAN/TRANS/PRED)."
      )
    }
    return(list(grad = "cpp", ad_ok = FALSE))
  }
  if (grad == "auto") {
    if (.nm_ad_pk_supported(model)) {
      return(list(grad = "ad", ad_ok = TRUE))
    }
    return(list(grad = "numeric", ad_ok = FALSE))
  }
  if (grad == "numeric") {
    return(list(grad = "numeric", ad_ok = .nm_ad_pk_supported(model)))
  }
  if (grad == "ad") {
    if (!.nm_ad_pk_supported(model)) {
      .nm_stop(
        "AD gradients are not supported for ADVAN ", model$ADVAN,
        " TRANS ", model$TRANS, ". Use grad = \"numeric\" with pk_engine = \"cpp\"."
      )
    }
    return(list(grad = "ad", ad_ok = TRUE))
  }
  list(grad = "numeric", ad_ok = FALSE)
}

#' @keywords internal
.nm_report_grad_backend <- function(model, est_grad, backend) {
  if (.nm_use_cpp_pop_grad(model, est_grad)) {
    return("cpp")
  }
  if (.nm_grad_uses_ad(est_grad)) {
    return(backend)
  }
  NULL
}

#' @keywords internal
.nm_bind_ad_eval_env <- function(env, ...) {
  if (.nm_any_ad(...)) {
    .ad_bind_math_ops(env)
    .ad_bind_control_ops(env)
  }
  invisible(env)
}
