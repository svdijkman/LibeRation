#' Build effective engine detail for a fit
#' @keywords internal
.nm_engine_detail <- function(method, grad, grad_requested, grad_backend,
                              pk_engine, engine, lik_config = NULL,
                              mstep = NULL, eta_fit = NULL) {
  structure(
    list(
      method = method,
      grad_requested = grad_requested %||% grad,
      grad_effective = grad,
      grad_backend = grad_backend,
      pk_engine = pk_engine,
      engine = engine,
      lik_error = lik_config$error %||% "propadd",
      lik_omega = lik_config$omega %||% "diag",
      mstep = mstep,
      eta_fit = eta_fit
    ),
    class = "nm_engine_detail"
  )
}

#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Print estimation engine diagnostics
#' @param x An \code{nm_engine_detail} object.
#' @param ... Unused.
#' @keywords internal
#' @examples
#' \dontrun{
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "FO",
#'               control = list(maxit = 5L, compute_inference = FALSE))
#' if (!is.null(fit$engine_detail)) print(fit$engine_detail)
#' }
#' @export
print.nm_engine_detail <- function(x, ...) {
  cat("Estimation engine detail\n")
  cat("  method:", x$method, "\n")
  cat("  grad (requested / effective / backend):",
      x$grad_requested, "/", x$grad_effective, "/", x$grad_backend %||% "-", "\n")
  cat("  pk_engine:", x$pk_engine, "  engine:", x$engine %||% "-", "\n")
  cat("  likelihood: error =", x$lik_error, ", omega =", x$lik_omega, "\n")
  if (!is.null(x$mstep)) cat("  M-step:", x$mstep, "\n")
  if (!is.null(x$eta_fit)) cat("  eta fit:", x$eta_fit, "\n")
  invisible(x)
}
