#' NPC and NPDE via simulation at fixed population parameters
#'
#' Monte Carlo predictive check (NPC) and normalized prediction distribution
#' errors (NPDE) using simulated datasets under the fitted model. Each replicate
#' uses the same design as the estimation data; conditional residuals (CWRES)
#' are compared to the observed data.
#'
#' @param fit An \code{nm_fit} object.
#' @param data Optional dataset; defaults to \code{fit$data}.
#' @param n_sim Number of simulated datasets.
#' @param seed Random seed for simulation.
#' @param refit_eta Re-fit subject ETAs on each replicate at fixed \code{fit$par}.
#' @param compute_npc If \code{TRUE}, compute NPC values.
#' @param compute_npde If \code{TRUE}, compute NPDE values.
#' @param pk_engine PK engine for simulation and residual calculation.
#' @param n_cores Parallel workers for simulation replicates.
#' @param control Passed to eta refit when \code{refit_eta = TRUE}.
#' @return The fit object with \code{gof} columns \code{NPC} and/or \code{NPDE},
#'   and metadata in \code{npc} / \code{npde}.
#' @examples
#' \dontrun{
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "FO",
#'               control = list(maxit = 5L, compute_inference = FALSE))
#' fit <- nm_add_npc_npde(fit, n_sim = 5L, refit_eta = FALSE)
#' }
#' @export
nm_add_npc_npde <- function(fit,
                            data = fit$data,
                            n_sim = 100L,
                            seed = 1L,
                            refit_eta = TRUE,
                            compute_npc = TRUE,
                            compute_npde = TRUE,
                            pk_engine = NULL,
                            n_cores = 1L,
                            control = list()) {
  if (is.null(fit) || is.null(fit$model)) {
    .nm_stop("fit must be an nm_fit object with a model.")
  }
  if (identical(fit$method, "BAYES")) {
    .nm_stop("NPC/NPDE is not supported for BAYES fits.")
  }
  if (!isTRUE(compute_npc) && !isTRUE(compute_npde)) {
    .nm_stop("At least one of compute_npc or compute_npde must be TRUE.")
  }
  n_sim <- max(1L, as.integer(n_sim))
  pk <- pk_engine %||% fit$pk_engine %||% "cpp"
  if (is.null(fit$gof) || !all(c("CWRES", "CRES") %in% names(fit$gof))) {
    fit <- nm_add_cwres(
      fit, data = data, refit_eta = isTRUE(refit_eta),
      pk_engine = pk, control = control
    )
  }
  gof <- fit$gof
  obs_idx <- which(gof$MDV == 0L & gof$EVID == 0L)
  if (length(obs_idx) == 0L) {
    .nm_stop("No observation rows (MDV=0, EVID=0) for NPC/NPDE.")
  }
  obs_cw <- gof$CWRES[obs_idx]
  na_cw <- !is.finite(obs_cw)
  if (any(na_cw)) {
    obs_cw[na_cw] <- gof$WRES[obs_idx][na_cw]
  }
  model <- fit$model
  sims <- nm_simulate(
    model, data, n_sim = n_sim, seed = seed,
    n_cores = n_cores, pk_engine = pk,
    theta = fit$theta, omega = fit$omega, sigma = fit$sigma
  )
  sim_mat <- matrix(NA_real_, length(obs_idx), n_sim)
  fit_stub <- fit
  fit_stub$data <- data
  for (k in seq_len(n_sim)) {
    fit_k <- fit_stub
    fit_k$data <- sims[[k]]
    fit_k <- tryCatch(
      nm_add_cwres(
        fit_k, data = sims[[k]], refit_eta = isTRUE(refit_eta),
        pk_engine = pk, control = control
      ),
      error = function(e) NULL
    )
    if (is.null(fit_k)) {
      next
    }
    gk <- fit_k$gof
    idx_k <- which(gk$MDV == 0L & gk$EVID == 0L)
    if (length(idx_k) != length(obs_idx)) {
      next
    }
    cw_k <- gk$CWRES[idx_k]
    bad <- !is.finite(cw_k)
    if (any(bad)) {
      cw_k[bad] <- gk$WRES[idx_k][bad]
    }
    sim_mat[, k] <- cw_k
  }
  npc_npde <- .nm_npc_npde_from_sim(obs_cw, sim_mat)
  if (isTRUE(compute_npc)) {
    if (!"NPC" %in% names(gof)) {
      gof$NPC <- NA_real_
    }
    gof$NPC[obs_idx] <- npc_npde$npc
    fit$npc <- list(
      n_sim = n_sim,
      n_ok = npc_npde$n_ok,
      refit_eta = isTRUE(refit_eta),
      seed = seed
    )
  }
  if (isTRUE(compute_npde)) {
    if (!"NPDE" %in% names(gof)) {
      gof$NPDE <- NA_real_
    }
    gof$NPDE[obs_idx] <- npc_npde$npde
    fit$npde <- list(
      n_sim = n_sim,
      n_ok = npc_npde$n_ok,
      refit_eta = isTRUE(refit_eta),
      seed = seed
    )
  }
  fit$gof <- gof
  fit$npc_npde <- list(
    n_sim = n_sim,
    n_ok = npc_npde$n_ok,
    refit_eta = isTRUE(refit_eta),
    seed = seed,
    compute_npc = isTRUE(compute_npc),
    compute_npde = isTRUE(compute_npde)
  )
  fit
}

#' @keywords internal
.nm_fit_has_npc <- function(fit) {
  if (is.null(fit)) {
    return(FALSE)
  }
  if (!is.null(fit$npc)) {
    return(TRUE)
  }
  if (!is.null(fit$gof) && "NPC" %in% names(fit$gof)) {
    any(is.finite(fit$gof$NPC))
  } else {
    FALSE
  }
}

#' @keywords internal
.nm_fit_has_npde <- function(fit) {
  if (is.null(fit)) {
    return(FALSE)
  }
  if (!is.null(fit$npde)) {
    return(TRUE)
  }
  if (!is.null(fit$gof) && "NPDE" %in% names(fit$gof)) {
    any(is.finite(fit$gof$NPDE))
  } else {
    FALSE
  }
}

#' @keywords internal
.nm_npc_npde_from_sim <- function(obs_cw, sim_mat) {
  n_obs <- length(obs_cw)
  npc <- rep(NA_real_, n_obs)
  npde <- rep(NA_real_, n_obs)
  n_ok <- 0L
  for (i in seq_len(n_obs)) {
    s <- sim_mat[i, ]
    s <- s[is.finite(s)]
    if (length(s) == 0L) {
      next
    }
    n_ok <- n_ok + 1L
    r <- rank(c(obs_cw[[i]], s), ties.method = "average")[[1L]]
    u <- (r - 0.5) / (length(s) + 1L)
    u <- max(min(u, 1 - 1e-12), 1e-12)
    npde[[i]] <- stats::qnorm(u)
    frac_hi <- mean(s >= obs_cw[[i]])
    frac_lo <- mean(s <= obs_cw[[i]])
    npc[[i]] <- 2 * min(frac_hi, frac_lo)
  }
  list(npc = npc, npde = npde, n_ok = n_ok)
}
