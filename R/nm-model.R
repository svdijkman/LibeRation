#' Define a NONMEM-style pharmacometric model
#'
#' Models are specified using \code{PRED} and \code{ERROR} code blocks with
#' \code{THETA()}, \code{ETA()}, and \code{ERR()} helpers, similar to NONMEM.
#'
#' @param INPUT Character vector of required input column names.
#' @param ADVAN ADVAN compartment structure (1, 2, 3, 4, 11, 12 supported).
#' @param TRANS TRANS parameterization.
#' @param SS Steady-state flag (0/1).
#' @param DOSECMP Dosing compartment number.
#' @param OBSCMP Observation compartment number.
#' @param PRED Character scalar: PK / parameter code block.
#' @param ERROR Character scalar: residual error model (\code{Y = ...}).
#' @param DES Optional additional code block (reserved).
#' @param THETAS Data frame with columns \code{THETA}, \code{Value}, optional \code{FIX}.
#' @param OMEGAS Data frame with columns \code{OMEGA}, \code{Value}.
#' @param COVARIATES Optional character vector of covariate column names in the dataset.
#' @param USE_ODE Use ODE solver (ADVAN 6/13) when \code{TRUE}.
#' @param SIGMAS Data frame with columns \code{SIGMA}, \code{Value}.
#' @param IOV Number of trailing ETAs for inter-occasion variability.
#' @param LIK_CONFIG Likelihood configuration from \code{\link{nm_lik_config}}.
#' @return An \code{nm_model} object.
#' @examples
#' sim <- nm_synthetic_theo(n_sub = 2L)
#' sim$model
#' @export
nm_model <- function(INPUT,
                     ADVAN = 2L,
                     TRANS = 2L,
                     SS = 0L,
                     DOSECMP = 1L,
                     OBSCMP = 1L,
                     PRED = "",
                     ERROR = "Y = F",
                     DES = "",
                     THETAS,
                     OMEGAS = data.frame(OMEGA = integer(), Value = numeric()),
                     SIGMAS = data.frame(SIGMA = integer(), Value = numeric()),
                     COVARIATES = NULL,
                     USE_ODE = FALSE,
                     IOV = 0L,
                     LIK_CONFIG = NULL) {
  lik <- if (is.null(LIK_CONFIG)) nm_lik_config(iov = IOV) else LIK_CONFIG
  if (IOV > 0L && (is.null(lik$iov) || lik$iov == 0L)) {
    lik$iov <- as.integer(IOV)
  }
  structure(
    list(
      INPUT = INPUT,
      ADVAN = as.integer(ADVAN),
      TRANS = as.integer(TRANS),
      SS = as.integer(SS),
      DOSECMP = as.integer(DOSECMP),
      OBSCMP = as.integer(OBSCMP),
      PRED = PRED,
      ERROR = ERROR,
      DES = DES,
      THETAS = THETAS,
      OMEGAS = OMEGAS,
      SIGMAS = SIGMAS,
      COVARIATES = COVARIATES,
      USE_ODE = isTRUE(USE_ODE),
      IOV = as.integer(IOV),
      LIK_CONFIG = lik
    ),
    class = "nm_model"
  )
}

#' @rdname nm_model
#' @method print nm_model
#' @param x An \code{nm_model} object.
#' @param ... Unused.
#' @examples
#' sim <- nm_synthetic_theo(n_sub = 2L)
#' print(sim$model)
#' @export
print.nm_model <- function(x, ...) {
  cat("NONMEM-style model\n")
  cat("  ADVAN:", x$ADVAN, " TRANS:", x$TRANS, " SS:", x$SS, "\n")
  cat("  THETAS:", nrow(x$THETAS), " OMEGAS:", nrow(x$OMEGAS),
      " SIGMAS:", nrow(x$SIGMAS), "\n")
  invisible(x)
}
