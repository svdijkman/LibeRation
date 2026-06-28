#' Model likelihood / residual configuration
#'
#' @param error Residual error model: \code{"propadd"} (default), \code{"add"},
#'   \code{"prop"}, \code{"log"}, or \code{"power"}.
#' @param omega Random-effect prior structure: \code{"diag"} or \code{"block2"}
#'   (2x2 block on first two ETAs plus diagonal remainder).
#' @param sigma_corr Residual correlation: \code{"indep"} or \code{"ar1"}.
#' @param iov Number of trailing ETAs treated as inter-occasion variability
#'   (separate omega diagonal entries).
#' @param ar1_rho AR(1) correlation for residual errors when \code{sigma_corr = "ar1"}.
#' @return A list with integer codes for C++.
#' @examples
#' nm_lik_config(error = "propadd", omega = "diag")
#' @export
nm_lik_config <- function(error = c("propadd", "add", "prop", "log", "power"),
                          omega = c("diag", "block2"),
                          sigma_corr = c("indep", "ar1"),
                          iov = 0L,
                          ar1_rho = 0.0) {
  error <- match.arg(error)
  omega <- match.arg(omega)
  sigma_corr <- match.arg(sigma_corr)
  err_code <- switch(
    error,
    propadd = 0L, add = 1L, prop = 2L, log = 3L, power = 4L
  )
  om_code <- switch(omega, diag = 0L, block2 = 1L)
  sc_code <- switch(sigma_corr, indep = 0L, ar1 = 1L)
  structure(
    list(
      error = error,
      omega = omega,
      sigma_corr = sigma_corr,
      iov = as.integer(iov),
      ar1_rho = as.numeric(ar1_rho),
      error_code = err_code,
      omega_code = om_code,
      sigma_corr_code = sc_code
    ),
    class = "nm_lik_config"
  )
}

#' @rdname nm_lik_config
#' @method print nm_lik_config
#' @param x An \code{nm_lik_config} object.
#' @param ... Unused.
#' @examples
#' print(nm_lik_config())
#' @export
print.nm_lik_config <- function(x, ...) {
  cat(
    "Likelihood config: error =", x$error,
    ", omega =", x$omega,
    ", sigma_corr =", x$sigma_corr,
    ", iov =", x$iov, "\n"
  )
  invisible(x)
}

#' @keywords internal
.nm_lik_config <- function(model) {
  if (!is.null(model$LIK_CONFIG)) {
    return(model$LIK_CONFIG)
  }
  nm_lik_config()
}

#' @keywords internal
.nm_sync_lik_config <- function(model) {
  cfg <- .nm_lik_config(model)
  if (exists("nm_lik_config_set", mode = "function")) {
    nm_lik_config_set(
      cfg$error_code, cfg$omega_code,
      cfg$sigma_corr_code %||% 0L,
      cfg$iov %||% 0L,
      cfg$ar1_rho %||% 0.0
    )
  }
  invisible(cfg)
}

#' @keywords internal
.nm_error_codes <- function() {
  c(propadd = 0L, add = 1L, prop = 2L, log = 3L, power = 4L)
}

#' Residual variance at observations (R path)
#' @keywords internal
.nm_residual_var <- function(f, sigma, error = "propadd", dvid = NULL, i = NULL) {
  s1 <- if (length(sigma) >= 1L) sigma[1] else 0
  s2 <- if (length(sigma) >= 2L) sigma[2] else 0
  if (!is.null(dvid) && !is.null(i) && length(sigma) >= 2L) {
    k <- as.integer(dvid[i])
    if (k >= 1L && (k * 2L) <= length(sigma)) {
      s1 <- sigma[k * 2L - 1L]
      s2 <- sigma[k * 2L]
    }
  }
  switch(
    error,
    add = pmax(s1 * s1, 1e-15),
    prop = pmax((s1 * f)^2, 1e-15),
    log = pmax(s1 * s1, 1e-15),
    power = pmax(s1^2 * pmax(abs(f), 1e-8)^(2 * s2), 1e-15),
    propadd = pmax((s1 * f)^2 + s2^2, 1e-15)
  )
}

#' @keywords internal
.nm_residual_nll_scalar <- function(dv, f, sigma, error = "propadd", dvid = NULL,
                                    ar1_rho = 0.0, sigma_corr = "indep") {
  n <- length(dv)
  if (n == 0L) {
    return(0)
  }
  if (identical(sigma_corr, "ar1") && n > 1L) {
    r <- if (identical(error, "log")) {
      log(pmax(dv, 1e-8)) - log(pmax(f, 1e-8))
    } else {
      dv - f
    }
    v <- vapply(seq_len(n), function(i) {
      .nm_residual_var(f[i], sigma, error, dvid, i)
    }, numeric(1))
    nll <- log(v[1L]) + r[1L]^2 / v[1L]
    for (i in seq_len(n)[-1L]) {
      ve <- v[i] * (1 - ar1_rho^2)
      e <- r[i] - ar1_rho * r[i - 1L]
      nll <- nll + log(max(ve, 1e-15)) + e^2 / max(ve, 1e-15)
    }
    return(nll)
  }
  nll <- 0
  for (i in seq_len(n)) {
    var <- .nm_residual_var(f[i], sigma, error, dvid, i)
    resid <- dv[i] - f[i]
    if (identical(error, "log")) {
      resid <- log(pmax(dv[i], 1e-8)) - log(pmax(f[i], 1e-8))
    }
    nll <- nll + log(var) + resid^2 / var
  }
  nll
}

#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x
