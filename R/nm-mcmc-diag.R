#' MCMC convergence diagnostics
#'
#' @param fit A BAYES \code{nm_fit} object.
#' @param probs Quantiles for posterior intervals.
#' @return List with \code{rhat}, \code{ess}, and \code{summary} tables.
#' @examples
#' \dontrun{
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "BAYES", n_burn = 5L, n_sample = 20L,
#'               control = list(compute_inference = FALSE))
#' nm_mcmc_diagnostics(fit)
#' }
#' @export
nm_mcmc_diagnostics <- function(fit, probs = c(0.025, 0.5, 0.975)) {
  if (!identical(fit$method, "BAYES")) {
    .nm_stop("MCMC diagnostics require method = BAYES.")
  }
  ch <- fit$chains
  diag_one <- function(mat, name_prefix) {
    if (is.null(mat) || nrow(mat) < 4L) {
      return(NULL)
    }
    n <- nrow(mat)
    half <- floor(n / 2)
    rhat <- vapply(seq_len(ncol(mat)), function(j) {
      x <- mat[, j]
      m <- mean(x)
      w <- var(x)
      x2 <- x[seq(half + 1L, n)]
      x1 <- x[seq_len(half)]
      b <- (half / 2) * (mean(x1) - mean(x2))^2
      var_plus <- ((half - 1) * w + b) / max(half - 1, 1)
      sqrt(var_plus / max(w, 1e-15))
    }, numeric(1))
    ess <- vapply(seq_len(ncol(mat)), function(j) {
      x <- mat[, j]
      acf1 <- stats::acf(x, plot = FALSE, lag.max = 1)$acf[2]
      rho <- max(0, acf1)
      n / (1 + 2 * rho)
    }, numeric(1))
    data.frame(
      param = paste0(name_prefix, seq_len(ncol(mat))),
      rhat = rhat,
      ess = ess,
      mean = colMeans(mat),
      t(vapply(seq_len(ncol(mat)), function(j) stats::quantile(mat[, j], probs), numeric(length(probs)))),
      row.names = NULL
    )
  }
  th <- diag_one(ch$theta, "THETA")
  om <- diag_one(ch$omega, "OMEGA")
  sg <- diag_one(ch$sigma, "SIGMA")
  structure(
    list(theta = th, omega = om, sigma = sg, n_keep = fit$n_keep),
    class = "nm_mcmc_diagnostics"
  )
}

#' @rdname nm_mcmc_diagnostics
#' @method print nm_mcmc_diagnostics
#' @param x An \code{nm_mcmc_diagnostics} object.
#' @param ... Unused.
#' @examples
#' \dontrun{
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "BAYES", n_burn = 5L, n_sample = 20L,
#'               control = list(compute_inference = FALSE))
#' print(nm_mcmc_diagnostics(fit))
#' }
#' @export
print.nm_mcmc_diagnostics <- function(x, ...) {
  cat("MCMC diagnostics (kept =", x$n_keep, ")\n")
  if (!is.null(x$theta)) {
    cat("THETA:\n")
    print(x$theta[, c("param", "rhat", "ess", "mean")])
  }
  invisible(x)
}
