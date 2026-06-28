#' Bootstrap standard errors for population parameters
#'
#' Nonparametric bootstrap: resample subjects with replacement, refit, summarise.
#'
#' @param fit An \code{nm_fit} object (used for model, data, method, controls).
#' @param n_boot Number of bootstrap replicates.
#' @param seed Random seed.
#' @param warm_start If \code{TRUE} (default), start each replicate from the
#'   original fit parameters; otherwise use control-file initials.
#' @param ... Passed to \code{nm_est} (excluding \code{data}).
#' @return A list with \code{se}, \code{bias}, \code{par_mean}, and \code{boot_pars}.
#' @examples
#' \dontrun{
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "FO",
#'               control = list(maxit = 5L, compute_inference = FALSE))
#' nm_bootstrap_se(fit, n_boot = 3L, control = list(maxit = 3L))
#' }
#' @export
nm_bootstrap_se <- function(fit, n_boot = 50L, seed = 1L, warm_start = TRUE, ...) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    .nm_stop("Package 'data.table' is required.")
  }
  set.seed(seed)
  model <- fit$model
  dat <- fit$data
  ids <- .nm_subject_ids(.nm_prepare_data(dat, model$INPUT, model))
  n_sub <- length(ids)
  dots <- list(...)
  method <- fit$method
  if (method == "BAYES") {
    .nm_stop("Bootstrap SE is not defined for BAYES fits; use posterior intervals.")
  }
  labels <- .nm_par_labels(model)
  n_boot <- max(1L, as.integer(n_boot))
  boot_pars <- matrix(NA_real_, n_boot, length(labels))
  colnames(boot_pars) <- labels
  start_par <- if (isTRUE(warm_start)) fit$par else NULL
  for (b in seq_len(n_boot)) {
    samp <- sample(ids, n_sub, replace = TRUE)
    boot_dat <- .nm_bootstrap_dataset(dat, samp)
    .nm_clear_cpp_subjects_cache()
    .nm_clear_pop_optim_cache()
    boot_args <- c(
      list(
        model = model,
        data = boot_dat,
        method = method,
        start = start_par
      ),
      dots
    )
    boot_fit <- tryCatch(
      do.call(nm_est, boot_args),
      error = function(e) NULL
    )
    row <- .nm_bootstrap_par_row(boot_fit, labels)
    if (!is.null(row)) {
      boot_pars[b, ] <- row
    }
  }
  ok <- apply(boot_pars, 2, function(x) sum(is.finite(x)))
  se <- vapply(
    labels,
    function(lbl) {
      x <- boot_pars[, lbl]
      x <- x[is.finite(x)]
      if (length(x) < 2L) {
        NA_real_
      } else {
        stats::sd(x)
      }
    },
    numeric(1)
  )
  se <- stats::setNames(se, labels)
  n_ok <- if (length(ok)) as.integer(max(ok)) else 0L
  par_mean <- if (any(ok >= 1L)) {
    stats::setNames(
      vapply(labels, function(lbl) {
        x <- boot_pars[, lbl]
        x <- x[is.finite(x)]
        if (length(x) == 0L) NA_real_ else mean(x)
      }, numeric(1)),
      labels
    )
  } else {
    stats::setNames(rep(NA_real_, length(labels)), labels)
  }
  bias <- par_mean - stats::setNames(as.numeric(fit$par), labels)
  list(
    se = se,
    bias = bias,
    par_mean = par_mean,
    boot_pars = boot_pars,
    n_ok = n_ok,
    n_ok_col = stats::setNames(as.integer(ok), labels),
    se_method = "sd_of_bootstrap_estimates",
    n_boot = n_boot,
    seed = seed
  )
}

#' @keywords internal
.nm_bootstrap_par_row <- function(boot_fit, labels) {
  if (is.null(boot_fit) || is.null(boot_fit$par)) {
    return(NULL)
  }
  p <- as.numeric(boot_fit$par)
  if (length(p) != length(labels)) {
    return(NULL)
  }
  if (!is.null(names(boot_fit$par))) {
    out <- rep(NA_real_, length(labels))
    names(out) <- labels
    nm <- intersect(names(boot_fit$par), labels)
    out[nm] <- as.numeric(boot_fit$par[nm])
    p <- as.numeric(out)
  }
  stats::setNames(p, labels)
}

#' Run bootstrap SE and attach results to a fit object
#'
#' @param fit An \code{nm_fit} object.
#' @param n_boot Number of bootstrap replicates; \code{0} skips bootstrap.
#' @param seed Random seed.
#' @param ... Passed to \code{\link{nm_bootstrap_se}} / \code{\link{nm_est}}.
#' @return The fit with a \code{bootstrap} element when \code{n_boot >= 1}.
#' @keywords internal
.nm_bootstrap_attach <- function(fit, n_boot = 0L, seed = 1L, ...) {
  n_boot <- as.integer(n_boot)
  if (is.null(fit) || n_boot < 1L) {
    return(fit)
  }
  fit$bootstrap <- nm_bootstrap_se(fit, n_boot = n_boot, seed = seed, ...)
  fit
}

#' @keywords internal
.nm_bootstrap_dataset <- function(data, ids_sample) {
  dat <- data.table::copy(data$data)
  pieces <- lapply(seq_along(ids_sample), function(j) {
    sub <- dat[dat$ID == ids_sample[j]]
    sub$ID <- j
    sub
  })
  new_dat <- data.table::rbindlist(pieces)
  structure(list(data = new_dat), class = "nm_dataset")
}
