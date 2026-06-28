#' @keywords internal
.nm_param_names <- function(model) {
  list(
    theta = paste0("THETA", model$THETAS$THETA),
    omega = paste0("OMEGA", model$OMEGAS$OMEGA),
    sigma = paste0("SIGMA", model$SIGMAS$SIGMA)
  )
}

#' @keywords internal
.nm_n_eta <- function(model) {
  nrow(model$OMEGAS)
}

#' @keywords internal
.nm_pack <- function(model, theta, omega, sigma) {
  c(theta, omega, sigma)
}

#' @keywords internal
.nm_unpack <- function(model, par) {
  n_th <- nrow(model$THETAS)
  n_om <- nrow(model$OMEGAS)
  n_sg <- nrow(model$SIGMAS)
  if (length(par) != n_th + n_om + n_sg) {
    .nm_stop("Parameter vector length mismatch.")
  }
  list(
    theta = par[seq_len(n_th)],
    omega = par[(n_th + 1L):(n_th + n_om)],
    sigma = par[(n_th + n_om + 1L):(n_th + n_om + n_sg)]
  )
}

#' @keywords internal
.nm_init_par <- function(model) {
  .nm_pack(
    model,
    model$THETAS$Value,
    model$OMEGAS$Value,
    model$SIGMAS$Value
  )
}

#' Default THETA lower/upper from an initial estimate (0.1x and 10x).
#' @keywords internal
.nm_theta_default_bounds <- function(value) {
  v <- as.numeric(value)
  if (!is.finite(v)) {
    v <- 1
  }
  av <- abs(v)
  if (av < 1e-15) {
    av <- 1
  }
  c(lower = max(1e-8, 0.1 * av), upper = 10 * av)
}

#' @keywords internal
.nm_omega_sigma_default_bounds <- function(value) {
  v <- as.numeric(value)
  if (!is.finite(v)) {
    v <- 0.1
  }
  av <- abs(v)
  if (av < 1e-15) {
    av <- 0.1
  }
  c(lower = max(1e-8, 0.1 * av), upper = max(av * 100, 1))
}

#' @keywords internal
.nm_par_lower <- function(model) {
  init <- .nm_init_par(model)
  lower <- rep(1e-8, length(init))
  n_th <- nrow(model$THETAS)
  n_om <- nrow(model$OMEGAS)
  n_sg <- nrow(model$SIGMAS)
  for (i in seq_len(n_th)) {
    val <- model$THETAS$Value[[i]]
    if ("Lower" %in% names(model$THETAS)) {
      lo <- model$THETAS$Lower[[i]]
      if (length(lo) == 1L && is.finite(lo) && !is.na(lo)) {
        lower[[i]] <- lo
        next
      }
    }
    lower[[i]] <- .nm_theta_default_bounds(val)[["lower"]]
  }
  for (i in seq_len(n_om)) {
    idx <- n_th + i
    bd <- .nm_omega_sigma_default_bounds(model$OMEGAS$Value[[i]])
    lower[[idx]] <- bd[["lower"]]
  }
  for (i in seq_len(n_sg)) {
    idx <- n_th + n_om + i
    bd <- .nm_omega_sigma_default_bounds(model$SIGMAS$Value[[i]])
    lower[[idx]] <- bd[["lower"]]
  }
  lower
}

#' @keywords internal
.nm_par_upper <- function(model) {
  init <- .nm_init_par(model)
  upper <- rep(Inf, length(init))
  n_th <- nrow(model$THETAS)
  n_om <- nrow(model$OMEGAS)
  n_sg <- nrow(model$SIGMAS)
  for (i in seq_len(n_th)) {
    val <- model$THETAS$Value[[i]]
    if ("Upper" %in% names(model$THETAS)) {
      up <- model$THETAS$Upper[[i]]
      if (length(up) == 1L && is.finite(up) && !is.na(up)) {
        upper[[i]] <- up
        next
      }
    }
    upper[[i]] <- .nm_theta_default_bounds(val)[["upper"]]
  }
  for (i in seq_len(n_om)) {
    idx <- n_th + i
    upper[[idx]] <- .nm_omega_sigma_default_bounds(model$OMEGAS$Value[[i]])[["upper"]]
  }
  for (i in seq_len(n_sg)) {
    idx <- n_th + n_om + i
    upper[[idx]] <- .nm_omega_sigma_default_bounds(model$SIGMAS$Value[[i]])[["upper"]]
  }
  upper
}

#' @keywords internal
.nm_fix_mask <- function(model) {
  n_th <- nrow(model$THETAS)
  n_om <- nrow(model$OMEGAS)
  n_sg <- nrow(model$SIGMAS)
  th_fix <- if ("FIX" %in% names(model$THETAS)) {
    as.logical(model$THETAS$FIX)
  } else {
    rep(FALSE, n_th)
  }
  if (length(th_fix) != n_th) {
    th_fix <- rep_len(th_fix, n_th)
  }
  c(th_fix, rep(FALSE, n_om + n_sg))
}

#' @keywords internal
.nm_apply_fix <- function(model, par) {
  init <- .nm_init_par(model)
  mask <- .nm_fix_mask(model)
  par[mask] <- init[mask]
  par
}

#' @keywords internal
.nm_pos_par <- function(x, eps = 1e-8) {
  pmax(x, eps)
}
