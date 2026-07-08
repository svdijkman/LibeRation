#' Build prior specification vectors for C++ Bayesian MCMC
#'
#' @keywords internal
.nm_bayes_prior_spec <- function(model, par0, prior = NULL) {
  p <- .nm_unpack(model, par0)
  n_th <- length(p$theta)
  n_om <- length(p$omega)
  n_sg <- length(p$sigma)
  if (is.null(prior)) {
    prior <- list()
  }
  parse_block <- function(block, n, default_meanlog = 0, default_sdlog = 1) {
    if (is.null(block)) {
      return(list(
        type = rep(0L, n),
        meanlog = rep(default_meanlog, n),
        sdlog = rep(default_sdlog, n)
      ))
    }
    tp <- block$type
    if (is.null(tp)) {
      tp <- rep("flat", n)
    }
    tp <- rep(as.character(tp), length.out = n)
    type <- ifelse(tp == "lognormal", 1L, 0L)
    mu <- block$meanlog
    if (is.null(mu)) {
      mu <- rep(default_meanlog, n)
    }
    sd <- block$sdlog
    if (is.null(sd)) {
      sd <- rep(default_sdlog, n)
    }
    list(
      type = as.integer(type),
      meanlog = as.numeric(rep(mu, length.out = n)),
      sdlog = as.numeric(rep(sd, length.out = n))
    )
  }
  th <- parse_block(
    prior$theta, n_th,
    default_meanlog = log(pmax(p$theta, 1e-8)),
    default_sdlog = 1
  )
  om <- parse_block(
    prior$omega, n_om,
    default_meanlog = log(pmax(sqrt(pmax(p$omega, 1e-8)), 1e-8)),
    default_sdlog = 1
  )
  sg <- parse_block(
    prior$sigma, n_sg,
    default_meanlog = log(pmax(p$sigma, 1e-8)),
    default_sdlog = 1
  )
  list(
    theta = th,
    omega = om,
    sigma = sg,
    theta_fix = as.logical(.nm_fix_mask(model))
  )
}

#' @keywords internal
.nm_bayes_step_scales <- function(model, par0, step_scale = NULL) {
  p <- .nm_unpack(model, par0)
  n_th <- length(p$theta)
  n_om <- length(p$omega)
  n_sg <- length(p$sigma)
  if (is.null(step_scale)) {
    return(list(
      theta = rep(0.05, n_th),
      omega = rep(0.05, n_om),
      sigma = rep(0.05, n_sg),
      eta = 1.0
    ))
  }
  if (is.numeric(step_scale) && length(step_scale) == 1L) {
    return(list(
      theta = rep(step_scale, n_th),
      omega = rep(step_scale, n_om),
      sigma = rep(step_scale, n_sg),
      eta = 1.0
    ))
  }
  list(
    theta = as.numeric(rep(step_scale$theta %||% 0.05, length.out = n_th)),
    omega = as.numeric(rep(step_scale$omega %||% 0.05, length.out = n_om)),
    sigma = as.numeric(rep(step_scale$sigma %||% 0.05, length.out = n_sg)),
    eta = as.numeric(step_scale$eta %||% 1)[1]
  )
}

#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x

#' @keywords internal
.nm_bayes_chain_mean <- function(chain) {
  if (is.null(chain) || !is.matrix(chain) || nrow(chain) == 0L) {
    return(NULL)
  }
  colMeans(chain)
}

#' @keywords internal
.nm_bayes_mcmc_cpp <- function(model, data, eta_mat, theta, omega, sigma,
                               n_burn, n_sample, n_thin, step_scale, prior,
                               par0, sampler = "mh", hmc_epsilon = 0.05,
                               hmc_leap = 10L, nuts_depth = 5L) {
  meta <- .nm_cpp_meta(model)
  subs <- .nm_cpp_subjects_cached(model, data)
  ps <- .nm_bayes_prior_spec(model, par0, prior)
  ss <- .nm_bayes_step_scales(model, par0, step_scale)
  nm_bayes_mcmc_cpp(
    eta = eta_mat,
    theta = as.numeric(theta),
    omega = as.numeric(omega),
    sigma = as.numeric(sigma),
    subjects = subs,
    pred_lines = meta$pred_lines,
    advan = meta$advan,
    trans = meta$trans,
    obs_cmp = meta$obs_cmp,
    dose_cmp = meta$dose_cmp,
    n_transit = meta$n_transit,
    use_ode = meta$use_ode,
    model_ss = meta$model_ss,
    n_burn = as.integer(n_burn),
    n_sample = as.integer(n_sample),
    n_thin = as.integer(n_thin),
    step_theta = ss$theta,
    step_omega = ss$omega,
    step_sigma = ss$sigma,
    step_eta = ss$eta,
    theta_prior_type = ps$theta$type,
    theta_prior_mu = ps$theta$meanlog,
    theta_prior_sd = ps$theta$sdlog,
    omega_prior_type = ps$omega$type,
    omega_prior_mu = ps$omega$meanlog,
    omega_prior_sd = ps$omega$sdlog,
    sigma_prior_type = ps$sigma$type,
    sigma_prior_mu = ps$sigma$meanlog,
    sigma_prior_sd = ps$sigma$sdlog,
    theta_fix = ps$theta_fix,
    des_lines = meta$des_lines,
    sampler = sampler,
    hmc_epsilon = hmc_epsilon,
    hmc_leap = as.integer(hmc_leap),
    nuts_depth = as.integer(nuts_depth)
  )
}

#' NONMEM-style full Bayesian MCMC (block Gibbs / Metropolis-Hastings)
#'
#' @keywords internal
.nm_est_bayes <- function(model, data, par0, backend = "cpp", grad = "auto",
                          pk_engine = "auto", engine = "auto",
                          control = list(),
                          n_burn = 100L, n_sample = 500L, n_thin = 1L,
                          seed = 1L, step_scale = NULL, prior = NULL,
                          store_eta_chain = TRUE,
                          sampler = c("mh", "hmc", "nuts"),
                          hmc_epsilon = 0.05, hmc_leap = 10L, nuts_depth = 5L) {
  sampler <- match.arg(sampler)
  pk_eff <- .nm_effective_pk_engine(pk_engine, grad)
  set.seed(seed)
  dat <- .nm_prepare_data(data, model$INPUT, model)
  ids <- .nm_subject_ids(dat)
  n_sub <- length(ids)
  n_eta <- .nm_n_eta(model)
  p0 <- .nm_unpack(model, par0)
  eta_mat <- matrix(0, n_sub, n_eta)
  if (n_eta > 0L) {
    eta_mat <- .nm_fit_all_eta(
      model, dat, p0$theta, p0$omega, p0$sigma, NULL,
      backend = "cpp", grad = "numeric", pk_engine = pk_eff, control = control
    )
  }
  use_cpp <- .nm_use_cpp_engine(engine, model, "BAYES") && .nm_cpp_capable(model)
  if (!use_cpp) {
    .nm_stop(
      "METHOD = BAYES requires a C++-capable model and engine = \"cpp\" or \"auto\"."
    )
  }
  .nm_est_progress_phase(
    "BAYES", "start",
    list(n_burn = n_burn, n_sample = n_sample, sampler = sampler, n_sub = n_sub)
  )
  .nm_est_progress_phase(
    "BAYES", "MCMC",
    list(n_burn = n_burn, n_sample = n_sample, n_thin = n_thin, sampler = sampler),
    log_msg = paste0(
      "BAYES MCMC: burn=", n_burn, " sample=", n_sample,
      " thin=", n_thin, " sampler=", sampler, " subjects=", n_sub
    )
  )
  res <- .nm_bayes_mcmc_cpp(
    model, data, eta_mat, p0$theta, p0$omega, p0$sigma,
    n_burn = n_burn, n_sample = n_sample, n_thin = n_thin,
    step_scale = step_scale, prior = prior, par0 = par0,
    sampler = sampler, hmc_epsilon = hmc_epsilon,
    hmc_leap = hmc_leap, nuts_depth = nuts_depth
  )
  acc_txt <- if (!is.null(res$acceptance)) {
    paste(format(unlist(res$acceptance), digits = 3), collapse = " ")
  } else {
    "NA"
  }
  .nm_est_progress_phase(
    "BAYES", "MCMC complete",
    list(
      log_posterior = res$log_posterior,
      n_burn = res$n_burn,
      n_sample = res$n_sample,
      acceptance = res$acceptance
    ),
    log_msg = paste0(
      "BAYES MCMC complete: log posterior = ",
      format(res$log_posterior, digits = 6),
      "  acceptance = ", acc_txt
    )
  )
  th_mean <- .nm_bayes_chain_mean(res$theta_chain)
  om_mean <- .nm_bayes_chain_mean(res$omega_chain)
  sg_mean <- .nm_bayes_chain_mean(res$sigma_chain)
  if (is.null(th_mean)) {
    th_mean <- res$theta
    om_mean <- res$omega
    sg_mean <- res$sigma
  }
  par_mean <- .nm_pack(model, th_mean, om_mean, sg_mean)
  par_mean <- .nm_apply_fix(model, par_mean)
  p_mean <- .nm_unpack(model, par_mean)
  eta_out <- res$eta
  eta_chain <- res$eta_chain
  if (!isTRUE(store_eta_chain)) {
    eta_chain <- NULL
  }
  structure(
    list(
      method = "BAYES",
      par = par_mean,
      theta = p_mean$theta,
      omega = p_mean$omega,
      sigma = p_mean$sigma,
      eta = eta_out,
      objective = -2 * res$log_posterior,
      log_posterior = res$log_posterior,
      convergence = 0L,
      grad = NA_character_,
      grad_requested = NULL,
      pk_engine = pk_eff,
      engine = "cpp",
      chains = list(
        theta = res$theta_chain,
        omega = res$omega_chain,
        sigma = res$sigma_chain,
        eta = eta_chain
      ),
      acceptance = res$acceptance,
      n_burn = res$n_burn,
      n_sample = res$n_sample,
      n_thin = res$n_thin,
      n_keep = res$n_keep,
      prior = prior,
      step_scale = step_scale,
      sampler = sampler,
      eval = list(
        pk = as.integer(res$n_pk_eval),
        mcmc_iter = as.integer(res$n_mcmc_iter),
        eta_proposed = as.integer(res$acceptance$eta[["proposed"]])
      )
    ),
    class = "nm_fit"
  )
}

#' Summary method for Bayesian \code{nm_fit} objects
#'
#' For non-Bayesian fits, returns the object invisibly. For \code{method = "BAYES"},
#' returns posterior mean and credible intervals for THETA, OMEGA, and SIGMA.
#'
#' @param object An \code{nm_fit} object.
#' @param ... Unused.
#' @examples
#' \dontrun{
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "BAYES", n_burn = 5L, n_sample = 10L,
#'               control = list(compute_inference = FALSE))
#' summary(fit)
#' }
#' @export
summary.nm_fit <- function(object, ...) {
  if (!identical(object$method, "BAYES")) {
    return(invisible(object))
  }
  ch <- object$chains
  if (is.null(ch$theta) || nrow(ch$theta) < 1L) {
    return(object)
  }
  pct <- function(x) stats::quantile(x, probs = c(0.025, 0.5, 0.975))
  thn <- paste0("THETA", object$model$THETAS$THETA)
  omn <- paste0("OMEGA", object$model$OMEGAS$OMEGA)
  sgn <- paste0("SIGMA", object$model$SIGMAS$SIGMA)
  theta_sum <- t(vapply(seq_along(thn), function(j) {
    c(mean = mean(ch$theta[, j]), pct(ch$theta[, j]))
  }, numeric(4)))
  rownames(theta_sum) <- thn
  colnames(theta_sum) <- c("mean", "2.5%", "50%", "97.5%")
  omega_sum <- if (length(omn) > 0L) {
    t(vapply(seq_along(omn), function(j) {
      c(mean = mean(ch$omega[, j]), pct(ch$omega[, j]))
    }, numeric(4)))
  } else {
    NULL
  }
  if (!is.null(omega_sum)) {
    rownames(omega_sum) <- omn
    colnames(omega_sum) <- colnames(theta_sum)
  }
  sigma_sum <- if (length(sgn) > 0L) {
    t(vapply(seq_along(sgn), function(j) {
      c(mean = mean(ch$sigma[, j]), pct(ch$sigma[, j]))
    }, numeric(4)))
  } else {
    NULL
  }
  if (!is.null(sigma_sum)) {
    rownames(sigma_sum) <- sgn
    colnames(sigma_sum) <- colnames(theta_sum)
  }
  structure(
    list(
      method = "BAYES",
      theta = theta_sum,
      omega = omega_sum,
      sigma = sigma_sum,
      acceptance = object$acceptance,
      n_keep = object$n_keep
    ),
    class = "summary.nm_fit_bayes"
  )
}

#' Print method for Bayesian MCMC summary objects
#'
#' @param x A \code{summary.nm_fit_bayes} object.
#' @param ... Unused.
#' @examples
#' \dontrun{
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "BAYES", n_burn = 5L, n_sample = 10L,
#'               control = list(compute_inference = FALSE))
#' print(summary(fit))
#' }
#' @export
print.summary.nm_fit_bayes <- function(x, ...) {
  cat("Bayesian MCMC summary (", x$method, ")\n", sep = "")
  cat("  kept samples:", x$n_keep, "\n")
  cat("  THETA:\n")
  print(x$theta)
  if (!is.null(x$omega) && nrow(x$omega) > 0L) {
    cat("  OMEGA:\n")
    print(x$omega)
  }
  if (!is.null(x$sigma) && nrow(x$sigma) > 0L) {
    cat("  SIGMA:\n")
    print(x$sigma)
  }
  if (!is.null(x$acceptance)) {
    cat("  acceptance (theta):", paste(round(x$acceptance$theta, 1), collapse = ", "), "\n")
    cat("  acceptance (omega):", paste(round(x$acceptance$omega, 1), collapse = ", "), "\n")
    cat("  acceptance (sigma):", paste(round(x$acceptance$sigma, 1), collapse = ", "), "\n")
    if (!is.null(x$acceptance$eta)) {
      cat("  acceptance (eta):",
          round(100 * x$acceptance$eta[["accepted"]] / max(1L, x$acceptance$eta[["proposed"]]), 1),
          "%\n")
    }
  }
  invisible(x)
}
