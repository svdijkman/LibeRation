#' Evaluate a C++ population joint objective and its taped derivatives
#'
#' The objective is minus twice the observation log likelihood plus the ETA
#' Gaussian prior. The complete prediction and likelihood calculation is
#' recorded by CppAD, so gradients and Hessians of this joint objective do not
#' use finite differences on a valid smooth tape path. This statement does not
#' extend to outer marginal estimation, conditional-mode sensitivity, or the
#' covariance-step bread, whose provenance is reported separately.
#'
#' @param model An `nm_model` or compiled `NMEngine`.
#' @param data NONMEM-style event data containing `DV` and `MDV`.
#' @param theta,eta,sigma,omega Parameter values.
#' @param gradient Return the taped joint-objective gradient.
#' @param hessian Return the taped joint-objective Hessian.
#' @return Objective value and requested derivatives.
#' @export
nm_objective <- function(model, data, theta = NULL, eta = NULL,
                         sigma = NULL, omega = NULL,
                         gradient = TRUE, hessian = FALSE) {
  engine <- if (inherits(model, "NMEngine")) model else nm_compile(model)
  theta <- theta %||% engine$model$THETAS$Value
  sigma <- sigma %||% engine$model$SIGMAS$Value
  omega <- omega %||% engine$model$OMEGAS$Value
  engine$objective(
    data, theta = theta, eta = eta, sigma = sigma, omega = omega,
    gradient = gradient, hessian = hessian
  )
}
