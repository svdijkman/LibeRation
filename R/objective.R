#' Evaluate the exact C++ population joint objective
#'
#' The objective is minus twice the observation log likelihood plus the ETA
#' Gaussian prior. The complete prediction and likelihood calculation is
#' recorded by CppAD, so gradients and Hessians do not use finite differences.
#'
#' @param model An `nm_model` or compiled `NMEngine`.
#' @param data NONMEM-style event data containing `DV` and `MDV`.
#' @param theta,eta,sigma,omega Parameter values.
#' @param gradient Return the exact gradient.
#' @param hessian Return the exact Hessian.
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
