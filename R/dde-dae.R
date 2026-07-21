#' Configure a delay differential-equation model
#'
#' Delayed states are referenced in `$DES` as `LAG(A(i), delay_name)`, where
#' `delay_name` is a positive assignment produced by `$PK/$PRED`. The compiled
#' method-of-steps solver uses fixed steps and differentiable linear history
#' interpolation.
#'
#' @param history Scalar or state-length numeric history before the first event.
#' @param step Fixed method-of-steps integration step.
#' @param interpolation Currently `"linear"`.
#' @param max_steps Maximum fixed steps per propagation call.
#' @param minimum_delay Optional declared lower bound for every delay. When
#'   supplied it must be at least `step`.
#' @export
nm_dde_config <- function(history = 0, step = 0.05,
                          interpolation = c("linear"), max_steps = 100000L,
                          minimum_delay = NULL) {
  interpolation <- match.arg(interpolation)
  history <- as.numeric(history)
  step <- as.numeric(step)
  max_steps <- as.integer(max_steps)
  if (!length(history) || any(!is.finite(history)) || length(step) != 1L ||
      !is.finite(step) || step <= 0 || length(max_steps) != 1L ||
      is.na(max_steps) || max_steps < 1L) {
    .nm_stop("DDE history must be finite and `step`/`max_steps` must be positive.")
  }
  if (!is.null(minimum_delay)) {
    minimum_delay <- as.numeric(minimum_delay)
    if (length(minimum_delay) != 1L || !is.finite(minimum_delay) ||
        minimum_delay < step) {
      .nm_stop("`minimum_delay` must be finite and at least as large as `step`.")
    }
  }
  structure(
    list(version = 1L, history = history, step = step,
         interpolation = interpolation, max_steps = max_steps,
         minimum_delay = minimum_delay, lags = list()),
    class = "nm_dde_config"
  )
}

.nm_dde_config <- function(config) {
  if (is.null(config)) return(NULL)
  if (inherits(config, "nm_dde_config")) return(config)
  if (!is.list(config)) .nm_stop("DDE_CONFIG must be created by `nm_dde_config()`.")
  config$version <- config$lags <- NULL
  do.call(nm_dde_config, config)
}

.nm_rewrite_dde_lags <- function(code, config) {
  if (is.null(config)) return(list(code = code, lags = list()))
  pattern <- "LAG\\s*\\(\\s*A\\s*\\(\\s*([0-9]+)\\s*\\)\\s*,\\s*([A-Za-z][A-Za-z0-9_]*)\\s*\\)"
  found <- regmatches(code, gregexpr(pattern, code, perl = TRUE))[[1L]]
  found <- unique(found[nzchar(found)])
  if (!length(found)) {
    .nm_stop("DDE_CONFIG requires at least one `LAG(A(i), delay_name)` expression in DES.")
  }
  lags <- vector("list", length(found))
  for (index in seq_along(found)) {
    state <- as.integer(sub(pattern, "\\1", found[[index]], perl = TRUE))
    delay <- sub(pattern, "\\2", found[[index]], perl = TRUE)
    input <- paste0("DDE_LAG_", index)
    code <- gsub(found[[index]], input, code, fixed = TRUE)
    lags[[index]] <- list(input = input, state = state, delay = delay)
  }
  list(code = code, lags = lags)
}

#' Configure a semi-explicit index-1 DAE
#'
#' Algebraic variables are ordinary symbols in `$DES`. The `ALG` block assigns
#' `RES(1)`, ..., one residual equation per variable. A fixed-iteration Newton
#' solve remains on the CppAD tape, including implicit parameter derivatives.
#'
#' @param variables Algebraic variable names.
#' @param initial Initial Newton values, one per variable.
#' @param tolerance Residual convergence tolerance.
#' @param maxit Maximum Newton iterations per derivative evaluation.
#' @param jacobian_step Relative central-difference step for the algebraic
#'   Jacobian.
#' @param sparsity Optional logical residual-by-variable sparsity pattern.
#' @export
nm_dae_config <- function(variables, initial = 0, tolerance = 1e-9,
                          maxit = 12L, jacobian_step = 1e-6,
                          sparsity = NULL) {
  variables <- trimws(as.character(variables))
  if (!length(variables) || anyNA(variables) || any(!nzchar(variables)) ||
      anyDuplicated(variables) || any(!grepl("^[A-Za-z][A-Za-z0-9_]*$", variables))) {
    .nm_stop("DAE algebraic variables must be unique identifier names.")
  }
  initial <- as.numeric(initial)
  if (length(initial) == 1L) initial <- rep(initial, length(variables))
  tolerance <- as.numeric(tolerance); maxit <- as.integer(maxit)
  jacobian_step <- as.numeric(jacobian_step)
  if (length(initial) != length(variables) || any(!is.finite(initial)) ||
      length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0 ||
      length(maxit) != 1L || is.na(maxit) || maxit < 1L || maxit > 100L ||
      length(jacobian_step) != 1L || !is.finite(jacobian_step) || jacobian_step <= 0) {
    .nm_stop("Invalid DAE initial values or Newton controls.")
  }
  if (!is.null(sparsity)) {
    sparsity <- as.matrix(sparsity)
    if (!identical(dim(sparsity), rep(length(variables), 2L))) {
      .nm_stop("DAE sparsity must be square with one row/column per algebraic variable.")
    }
    storage.mode(sparsity) <- "logical"
  }
  structure(
    list(version = 1L, variables = variables, initial = initial,
         tolerance = tolerance, maxit = maxit,
         jacobian_step = jacobian_step, sparsity = sparsity),
    class = "nm_dae_config"
  )
}

.nm_dae_config <- function(config) {
  if (is.null(config)) return(NULL)
  if (inherits(config, "nm_dae_config")) return(config)
  if (!is.list(config)) .nm_stop("DAE_CONFIG must be created by `nm_dae_config()`.")
  config$version <- NULL
  do.call(nm_dae_config, config)
}

.nm_compile_alg_ir <- function(alg, config, pred_ir, n_state,
                               n_theta, n_eta, covariates) {
  if (is.null(config)) return(NULL)
  code <- paste(alg, collapse = "\n")
  code <- .nm_rewrite_ode_indexing(code)
  code <- gsub("\\bRES\\s*\\(\\s*([0-9]+)\\s*\\)", "DAE_RES_\\1", code, perl = TRUE)
  outputs <- paste0("DAE_RES_", seq_along(config$variables))
  declared <- unique(c(
    pred_ir$output_names, paste0("A_", seq_len(n_state)), "T",
    config$variables,
    if (n_theta > 0L) paste0("THETA_", seq_len(n_theta)) else character(),
    if (n_eta > 0L) paste0("ETA_", seq_len(n_eta)) else character(),
    as.character(covariates %||% character())
  ))
  ir <- LibeRtAD::ad_ir(code, inputs = declared, outputs = outputs)
  missing <- setdiff(outputs, ir$output_names)
  if (length(missing)) .nm_stop("ALG must assign RES(1) through RES(", length(outputs), ").")
  ir
}

#' @export
print.nm_dde_config <- function(x, ...) {
  cat("LibeRation experimental DDE configuration\n")
  cat("  step:", x$step, " interpolation:", x$interpolation,
      " lags:", length(x$lags), "\n")
  invisible(x)
}

#' @export
print.nm_dae_config <- function(x, ...) {
  cat("LibeRation experimental index-1 DAE configuration\n")
  cat("  algebraic variables:", paste(x$variables, collapse = ", "),
      " Newton iterations:", x$maxit, "\n")
  invisible(x)
}

