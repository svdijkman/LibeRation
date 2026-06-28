#' @keywords internal
.nm_par_scale_enabled <- function(control) {
  if (!is.null(control$par_scale)) {
    return(isTRUE(control$par_scale))
  }
  isTRUE(getOption("LibeRation.par_scale", FALSE))
}

#' @keywords internal
.nm_par_scale_vector <- function(model, par) {
  labels <- .nm_par_labels(model)
  scale <- pmax(abs(par), ifelse(grepl("^OMEGA|^SIGMA", labels), 0.01, 0.1))
  fix <- .nm_fix_mask(model)
  scale[fix] <- 1
  scale
}

#' @keywords internal
.nm_par_to_scaled <- function(par, scale) {
  par / scale
}

#' @keywords internal
.nm_par_from_scaled <- function(x, scale) {
  x * scale
}

#' Wrap optim objective/gradient for internal parameter scaling.
#' @keywords internal
.nm_wrap_scaled_optim <- function(model, par0, free, f_free, g_free, scale) {
  x0 <- .nm_par_to_scaled(par0[free], scale[free])
  fn <- function(x) {
    f_free(.nm_par_from_scaled(x, scale[free]))
  }
  gr <- if (!is.null(g_free)) {
    function(x) {
      g_free(.nm_par_from_scaled(x, scale[free])) * scale[free]
    }
  } else {
    NULL
  }
  list(par = x0, fn = fn, gr = gr)
}
