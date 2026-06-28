#' Run a NONMEM-style task (estimate or simulate)
#'
#' @param type Task type: \code{"est"} or \code{"sim"}.
#' @param model An \code{nm_model} object.
#' @param data An \code{nm_dataset} object.
#' @param method Estimation method when \code{type = "est"}.
#' @param ... Passed to \code{\link{nm_est}} or simulation helpers.
#' @return Task result (fit or simulated dataset).
#' @examples
#' \dontrun{
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_task("est", sim$model, sim$data, method = "FO",
#'                control = list(maxit = 5L, compute_inference = FALSE))
#' print(fit)
#' }
#' @export
nm_task <- function(type = c("est", "sim"),
                    model,
                    data,
                    method = "FO",
                    ...) {
  type <- match.arg(type)
  switch(
    type,
    est = nm_est(model, data, method = method, ...),
    sim = .nm_task_sim(model, data, ...)
  )
}

#' @keywords internal
.nm_task_sim <- function(model, data, theta = NULL, omega = NULL, sigma = NULL,
                         eta = NULL, seed = 1L, pk_engine = "cpp", rep = 1L,
                         with_amounts = FALSE, n_state = 7L) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    .nm_stop("Package 'data.table' is required for simulation.")
  }
  set.seed(seed)
  dat <- data.table::as.data.table(.nm_prepare_data(data, model$INPUT, model))
  if (is.null(theta)) {
    theta <- model$THETAS$Value
  }
  if (is.null(omega)) {
    omega <- model$OMEGAS$Value
  }
  if (is.null(sigma)) {
    sigma <- model$SIGMAS$Value
  }
  ids <- .nm_subject_ids(dat)
  n_eta <- .nm_n_eta(model)
  n_sub <- length(ids)
  eta_mat <- NULL
  if (is.null(eta) && n_eta > 0L && n_sub > 0L) {
    eta_mat <- matrix(
      stats::rnorm(n_sub * n_eta, sd = sqrt(pmax(omega, 0))),
      nrow = n_sub,
      ncol = n_eta,
      byrow = TRUE
    )
  }
  prop_sd <- if (length(sigma) >= 1L) sigma[[1L]] else 0
  add_sd <- if (length(sigma) >= 2L) sigma[[2L]] else 0
  ipred_vec <- numeric(nrow(dat))
  y_vec <- rep(NA_real_, nrow(dat))
  amt_cols <- character(0L)
  for (j in seq_len(n_sub)) {
    id <- ids[[j]]
    sub <- .nm_subject_slice(dat, id)
    eta_j <- if (is.null(eta)) {
      if (n_eta > 0L) eta_mat[j, , drop = TRUE] else numeric()
    } else if (is.matrix(eta)) {
      eta[j, , drop = TRUE]
    } else {
      eta
    }
    pred <- .nm_subject_ipred(
      model, sub, theta, omega, eta_j, sigma = sigma, pk_engine = pk_engine,
      with_amounts = with_amounts, n_state = n_state
    )
    idx_id <- which(dat$ID == id)
    ipred_vec[idx_id] <- pred$ipred
    if (isTRUE(with_amounts) && !is.null(pred$amounts)) {
      amt_mat <- pred$amounts
      if (length(amt_cols) == 0L) {
        amt_cols <- colnames(amt_mat)
        if (is.null(amt_cols)) {
          amt_cols <- sprintf("A%d", seq_len(ncol(amt_mat)))
        }
        for (cn in amt_cols) {
          data.table::set(dat, j = cn, value = NA_real_)
        }
      }
      for (k in seq_len(ncol(amt_mat))) {
        cn <- amt_cols[[k]]
        data.table::set(dat, i = idx_id, j = cn, value = amt_mat[, k])
      }
    }
    if (length(pred$obs_idx) > 0L) {
      f <- pred$F
      y_sim <- f * (1 + stats::rnorm(length(f), sd = prop_sd)) +
        stats::rnorm(length(f), sd = add_sd)
      idx_obs <- which(dat$ID == id & dat$MDV == 0L & dat$EVID == 0L)
      if (length(idx_obs) == length(y_sim)) {
        y_vec[idx_obs] <- y_sim
      }
    }
  }
  data.table::set(dat, j = "IPRED", value = ipred_vec)
  data.table::set(dat, j = "Y", value = y_vec)
  obs_ok <- !is.na(y_vec)
  if (any(obs_ok)) {
    data.table::set(dat, i = which(obs_ok), j = "DV", value = y_vec[obs_ok])
  }
  out <- as.data.frame(dat)
  out$REP <- as.integer(rep)
  out
}

#' Project container for NM workflows
#'
#' Supports save/load and fit chaining.
#'
#' @param name Project name.
#' @return An \code{nm_proj} object.
#' @examples
#' proj <- nm_proj("demo")
#' sim <- nm_synthetic_theo(n_sub = 2L)
#' nm_proj_set(proj, model = sim$model, data = sim$data)
#' @export
nm_proj <- function(name = "project") {
  structure(
    list(
      name = name,
      model = NULL,
      data = NULL,
      fit = NULL,
      history = list()
    ),
    class = "nm_proj"
  )
}

#' @rdname nm_proj
#' @param proj An \code{nm_proj} object.
#' @param model Optional \code{nm_model}.
#' @param data Optional \code{nm_dataset}.
#' @param fit Optional \code{nm_fit} (appended to history).
#' @examples
#' proj <- nm_proj("demo")
#' sim <- nm_synthetic_theo(n_sub = 2L)
#' nm_proj_set(proj, model = sim$model, data = sim$data)
#' @export
nm_proj_set <- function(proj, model = NULL, data = NULL, fit = NULL) {
  if (!is.null(model)) proj$model <- model
  if (!is.null(data)) proj$data <- data
  if (!is.null(fit)) {
    proj$fit <- fit
    proj$history[[length(proj$history) + 1L]] <- fit
  }
  proj
}

#' @rdname nm_proj
#' @param path File path for save/load.
#' @examples
#' proj <- nm_proj("demo")
#' path <- tempfile(fileext = ".rds")
#' nm_proj_save(proj, path)
#' @export
nm_proj_save <- function(proj, path) {
  saveRDS(proj, path)
  invisible(path)
}

#' @rdname nm_proj
#' @param path File path saved by \code{\link{nm_proj_save}}.
#' @examples
#' proj <- nm_proj("demo")
#' path <- tempfile(fileext = ".rds")
#' nm_proj_save(proj, path)
#' nm_proj_load(path)
#' @export
nm_proj_load <- function(path) {
  readRDS(path)
}

#' @rdname nm_proj
#' @method print nm_proj
#' @param x An \code{nm_proj} object.
#' @param ... Unused.
#' @examples
#' print(nm_proj("demo"))
#' @export
print.nm_proj <- function(x, ...) {
  cat("NM project:", x$name, "\n")
  if (!is.null(x$model)) cat("  model: defined\n")
  if (!is.null(x$data)) cat("  data:", nrow(x$data$data), "rows\n")
  if (!is.null(x$fit)) cat("  fit:", x$fit$method, "\n")
  if (length(x$history) > 0L) {
    cat("  history:", length(x$history), "fit(s)\n")
  }
  invisible(x)
}
