#' @keywords internal
.nm_make_env <- function(theta, omega, sigma, eta = NULL, err = NULL) {
  THETA <- function(i) {
    i <- as.integer(i)
    if (i < 1L || i > length(theta)) {
      .nm_stop("THETA(", i, ") out of range.")
    }
    if (is.list(theta)) {
      return(theta[[i]])
    }
    theta[i]
  }
  OMEGA <- function(i) {
    i <- as.integer(i)
    if (i < 1L || i > length(omega)) {
      .nm_stop("OMEGA(", i, ") out of range.")
    }
    omega[i]
  }
  SIGMA <- function(i) {
    i <- as.integer(i)
    if (i < 1L || i > length(sigma)) {
      .nm_stop("SIGMA(", i, ") out of range.")
    }
    sigma[i]
  }
  ETA <- function(i) {
    i <- as.integer(i)
    if (is.null(eta)) {
      return(0)
    }
    if (i < 1L || i > length(eta)) {
      .nm_stop("ETA(", i, ") out of range.")
    }
    if (is.list(eta)) {
      return(eta[[i]])
    }
    eta[i]
  }
  ERR <- function(i) {
    i <- as.integer(i)
    if (is.null(err)) {
      return(0)
    }
    if (i < 1L || i > length(err)) {
      .nm_stop("ERR(", i, ") out of range.")
    }
    err[i]
  }
  list2env(
    list(THETA = THETA, OMEGA = OMEGA, SIGMA = SIGMA, ETA = ETA, ERR = ERR),
    parent = baseenv()
  )
}

#' @keywords internal
.nm_eval_block <- function(code, env, extra = list()) {
  lines <- .nm_split_lines(code)
  if (length(lines) == 0L) {
    return(list())
  }
  for (nm in names(extra)) {
    env[[nm]] <- extra[[nm]]
  }
  out <- list()
  for (line in lines) {
    expr <- parse(text = line)
    if (length(expr) == 0L) {
      next
    }
    e1 <- expr[[1]][[1]]
    if (is.call(expr[[1]]) && as.character(e1) %in% c("<-", "=")) {
      nm <- as.character(expr[[1]][[2]])
      res <- eval(expr, envir = env)
      out[[nm]] <- res
      env[[nm]] <- res
    } else {
      eval(expr, envir = env)
    }
  }
  out
}

#' @keywords internal
.nm_eval_pred <- function(model, theta, omega, eta, covariates = list()) {
  env <- .nm_make_env(theta, omega, rep(0, length(omega)), eta = eta)
  for (nm in names(covariates)) {
    env[[nm]] <- covariates[[nm]]
  }
  if (.nm_any_ad(theta, omega, eta)) {
    .ad_bind_math_ops(env)
    .ad_bind_control_ops(env)
  }
  .nm_eval_block(model$PRED, env)
}

#' @keywords internal
.nm_eval_error <- function(model, theta, omega, sigma, eta, err, F_val) {
  env <- .nm_make_env(theta, omega, sigma, eta = eta, err = err)
  env$F <- F_val
  if (.nm_any_ad(theta, omega, sigma, eta, err, F_val)) {
    .ad_bind_math_ops(env)
    .ad_bind_control_ops(env)
  }
  .nm_eval_block(model$ERROR, env)
}

#' @keywords internal
.nm_extract_y <- function(error_out) {
  if ("Y" %in% names(error_out)) {
    return(error_out$Y)
  }
  .nm_stop("ERROR block must assign Y = ...")
}
