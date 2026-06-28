#' Export an \code{nm_model} to nlmixr2-style model text
#'
#' @param model An \code{nm_model} object.
#' @param name Model function name.
#' @return Character vector of R code suitable for \code{nlmixr2}.
#' @examples
#' sim <- nm_synthetic_theo(n_sub = 2L)
#' cat(nm_to_nlmixr2(sim$model), sep = "\n")
#' @export
nm_to_nlmixr2 <- function(model, name = "nm_export") {
  th <- model$THETAS
  th_lines <- paste0("  ", paste0("THETA", th$THETA, " <- fix(", th$Value, ")"), collapse = "\n")
  om <- model$OMEGAS
  om_lines <- if (nrow(om) > 0L) {
    paste0("  eta ~ ", paste0("c(", om$Value, ")", collapse = " + "))
  } else {
    "  eta ~ 0"
  }
  sg <- model$SIGMAS
  sg_lines <- if (nrow(sg) > 0L) {
    paste0("  ", paste0("SIGMA", sg$SIGMA, " <- ", sg$Value), collapse = "\n")
  } else {
    "  sigma <- 0.1"
  }
  pred <- gsub("THETA\\((\\d+)\\)", "THETA\\1", model$PRED)
  pred <- gsub("ETA\\((\\d+)\\)", "ETA.\\1", pred)
  c(
    paste0(name, " <- function() {"),
    "  ini({",
    th_lines,
    paste0("  ", om_lines),
    paste0("  ", sg_lines),
    "  })",
    "  model({",
    paste0("    ", strsplit(trimws(pred), "\n")[[1]]),
    paste0("    ", gsub("Y\\s*=\\s*", "DV ~ ", model$ERROR)),
    "  })",
    "}"
  )
}

#' Import a minimal nlmixr2 fit object into LibeRation format
#'
#' Requires \pkg{nlmixr2} for full parsing; falls back to theta/omega/sigma vectors.
#'
#' @param fit An \code{nlmixr2} fit object or list with \code{theta}, \code{omega}, \code{sigma}.
#' @param model Optional existing \code{nm_model} template.
#' @return List with \code{model} and \code{par}.
#' @examples
#' sim <- nm_synthetic_theo(n_sub = 2L)
#' nm_from_nlmixr2(list(theta = sim$model$THETAS$Value), model = sim$model)
#' @export
nm_from_nlmixr2 <- function(fit, model = NULL) {
  if (is.null(model)) {
    model <- nm_model(
      INPUT = c("ID", "TIME", "DV", "AMT", "MDV", "EVID", "CMT"),
      ADVAN = 2L,
      TRANS = 2L,
      PRED = "CL = THETA(1)\nVC = THETA(2)",
      THETAS = data.frame(THETA = 1:2, Value = c(1, 10), FIX = FALSE)
    )
  }
  th <- fit$theta
  if (is.null(th) && !is.null(fit$parFixed)) {
    th <- fit$parFixed$theta
  }
  om <- fit$omega
  if (is.null(om) && !is.null(fit$omega)) {
    om <- fit$omega
  }
  sg <- fit$sigma
  if (is.null(sg) && !is.null(fit$sigma)) {
    sg <- fit$sigma
  }
  if (is.matrix(om)) {
    om <- sqrt(diag(om))
  }
  if (length(th) == 0L) {
    th <- model$THETAS$Value
  }
  if (length(om) == 0L) {
    om <- if (nrow(model$OMEGAS) > 0L) model$OMEGAS$Value else numeric()
  }
  if (length(sg) == 0L) {
    sg <- if (nrow(model$SIGMAS) > 0L) model$SIGMAS$Value else 0.1
  }
  par <- .nm_pack(model, as.numeric(th), as.numeric(om), as.numeric(sg))
  list(model = model, par = par, source = "nlmixr2")
}
