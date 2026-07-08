#' Stepwise covariate modelling (SCM)
#'
#' Automated forward selection followed by backward elimination of
#' covariate-parameter relationships, driven by the change in objective function
#' value (\eqn{\Delta OFV}) against \eqn{\chi^2} thresholds. Each candidate adds
#' one covariate effect (a new THETA) to a model parameter; the relationship can
#' be exponential, power, or linear.
#'
#' @param model A base \code{nm_model} (no candidate covariate effects).
#' @param data An \code{nm_dataset} or data.frame containing the covariate
#'   columns.
#' @param candidates A list of candidate relationships. Each element is a list
#'   with \code{parameter} (the model parameter to modify, e.g. \code{"CL"}),
#'   \code{covariate} (a data column name), and optional \code{relationship}
#'   (\code{"exponential"} default, \code{"power"}, or \code{"linear"}). You may
#'   also pass a data.frame with columns \code{parameter}, \code{covariate},
#'   \code{relationship}.
#' @param method Estimation method for each fit (default \code{"FOCE"}).
#' @param forward_alpha Significance level for forward inclusion (default 0.05).
#' @param backward_alpha Significance level for backward elimination
#'   (default 0.01, i.e. stricter than forward).
#' @param control Control list passed to \code{\link{nm_est}}.
#' @param verbose Logical; print per-step progress.
#' @param ... Passed to \code{\link{nm_est}}.
#' @return An object of class \code{nm_scm} with elements \code{steps}
#'   (per-step data.frame with dOFV, AIC, BIC), \code{retained} (list of the
#'   selected effects), \code{forest} (tidy data.frame of retained effect
#'   estimates + CI, suitable for \code{\link{nm_forest_plot}}), \code{base_ofv},
#'   \code{final_ofv}, and \code{final_fit}.
#' @examples
#' \dontrun{
#' scm <- nm_scm(model, data, candidates = list(
#'   list(parameter = "CL", covariate = "WT", relationship = "power"),
#'   list(parameter = "CL", covariate = "AGE")
#' ))
#' scm$steps
#' nm_forest_plot(scm)
#' }
#' @export
nm_scm <- function(model,
                   data,
                   candidates,
                   method = "FOCE",
                   forward_alpha = 0.05,
                   backward_alpha = 0.01,
                   control = list(),
                   verbose = TRUE,
                   ...) {
  candidates <- .nm_scm_normalize_candidates(candidates)
  if (length(candidates) == 0L) {
    .nm_stop("No candidates supplied.")
  }
  cand_covs <- unique(vapply(candidates, function(c) c$covariate, character(1L)))
  input_aug <- unique(c(as.character(model$INPUT), cand_covs))
  dat <- .nm_prepare_data(data, input_aug, model)
  refs <- .nm_scm_covariate_refs(dat, candidates)
  n_obs <- sum(dat$EVID == 0L & dat$MDV == 0L, na.rm = TRUE)
  ctl <- c(list(compute_inference = FALSE), control)

  n_base_th <- nrow(model$THETAS)

  # Build a warm-start (packed) vector for a trial model from a parent fit.
  # Covariate models are optimised on the slower numeric R likelihood path,
  # which can stall at a poor point when started from scratch; warm-starting
  # from the parent estimates (new/removed effects handled by index) begins the
  # optimiser at the parent optimum so a nested fit can only improve. This is
  # essential for valid (non-negative) dOFV comparisons.
  build_start <- function(effs, warm_fit, warm_effs) {
    if (is.null(warm_fit)) {
      return(NULL)
    }
    th <- as.numeric(warm_fit$theta)
    if (length(th) < n_base_th) {
      return(NULL)
    }
    th_base <- th[seq_len(n_base_th)]
    th_eff <- vapply(effs, function(e) {
      j <- 0L
      for (w in seq_along(warm_effs)) {
        if (.nm_scm_eff_eq(warm_effs[[w]], e)) {
          j <- w
          break
        }
      }
      if (j > 0L && length(th) >= n_base_th + j) th[n_base_th + j] else 0
    }, numeric(1L))
    m <- .nm_scm_build_model(model, effs, refs)
    tryCatch(
      .nm_pack(m, c(th_base, th_eff), warm_fit$omega, warm_fit$sigma),
      error = function(e) NULL
    )
  }

  fit_effects <- function(effs, warm_fit = NULL, warm_effs = list()) {
    m <- .nm_scm_build_model(model, effs, refs)
    # The C++ subject cache is keyed on core data columns only and does NOT
    # distinguish covariate sets, so a fit with a different covariate structure
    # on the same data could otherwise reuse stale (covariate-free) subjects.
    # Clear it between SCM fits so each model rebuilds its own bound subjects.
    if (exists(".nm_clear_cpp_subjects_cache", mode = "function")) {
      .nm_clear_cpp_subjects_cache()
    }
    start <- build_start(effs, warm_fit, warm_effs)
    f <- nm_est(m, data, method = method, start = start, control = ctl, ...)
    list(model = m, fit = f, ofv = as.numeric(f$objective),
         npar = length(f$theta) + length(f$omega) + length(f$sigma))
  }

  base <- fit_effects(list())
  base_ofv <- base$ofv
  steps <- list()
  add_step <- function(action, cand, dofv, ofv, npar, thr, sig) {
    steps[[length(steps) + 1L]] <<- data.frame(
      step = length(steps) + 1L,
      action = action,
      parameter = cand$parameter %||% NA_character_,
      covariate = cand$covariate %||% NA_character_,
      relationship = cand$relationship %||% NA_character_,
      dOFV = round(dofv, 4),
      ofv = round(ofv, 4),
      npar = npar,
      AIC = round(ofv + 2 * npar, 4),
      BIC = round(ofv + npar * log(max(n_obs, 1)), 4),
      threshold = round(thr, 4),
      significant = sig,
      stringsAsFactors = FALSE
    )
  }
  add_step("base", list(), 0, base_ofv, base$npar,
           0, TRUE)

  # ---- Forward selection ----
  thr_fwd <- stats::qchisq(1 - forward_alpha, df = 1L)
  included <- list()
  current_ofv <- base_ofv
  current_fit <- base$fit
  repeat {
    remaining <- candidates[!vapply(candidates, function(c) {
      .nm_scm_in(c, included)
    }, logical(1L))]
    if (length(remaining) == 0L) {
      break
    }
    best <- NULL
    for (cand in remaining) {
      trial <- fit_effects(c(included, list(cand)),
                           warm_fit = current_fit, warm_effs = included)
      dofv <- current_ofv - trial$ofv
      if (isTRUE(verbose)) {
        message(sprintf("  [forward] %s~%s (%s): dOFV = %.3f",
                        cand$parameter, cand$covariate,
                        cand$relationship, dofv))
      }
      if (is.null(best) || dofv > best$dofv) {
        best <- list(cand = cand, dofv = dofv, trial = trial)
      }
    }
    sig <- is.finite(best$dofv) && best$dofv > thr_fwd
    add_step("forward", best$cand, best$dofv, best$trial$ofv,
             best$trial$npar, thr_fwd, sig)
    if (!sig) {
      break
    }
    included <- c(included, list(best$cand))
    current_ofv <- best$trial$ofv
    current_fit <- best$trial$fit
  }

  # ---- Backward elimination ----
  thr_bwd <- stats::qchisq(1 - backward_alpha, df = 1L)
  repeat {
    if (length(included) <= 1L) {
      break
    }
    worst <- NULL
    for (i in seq_along(included)) {
      reduced <- included[-i]
      trial <- fit_effects(reduced,
                           warm_fit = current_fit, warm_effs = included)
      # dOFV of removing effect i: OFV increases when a real effect is dropped.
      dofv <- trial$ofv - current_ofv
      if (isTRUE(verbose)) {
        message(sprintf("  [backward] drop %s~%s: dOFV = %.3f",
                        included[[i]]$parameter, included[[i]]$covariate, dofv))
      }
      if (is.null(worst) || dofv < worst$dofv) {
        worst <- list(idx = i, dofv = dofv, trial = trial)
      }
    }
    # Retain the effect only if its removal is significant (dOFV >= threshold).
    retain <- is.finite(worst$dofv) && worst$dofv > thr_bwd
    cand_i <- included[[worst$idx]]
    add_step("backward", cand_i, worst$dofv, worst$trial$ofv,
             worst$trial$npar, thr_bwd, !retain)
    if (retain) {
      break
    }
    included <- included[-worst$idx]
    current_ofv <- worst$trial$ofv
    current_fit <- worst$trial$fit
  }

  # ---- Final model + inference for retained effects ----
  final_model <- .nm_scm_build_model(model, included, refs)
  if (exists(".nm_clear_cpp_subjects_cache", mode = "function")) {
    .nm_clear_cpp_subjects_cache()
  }
  final_start <- build_start(included, current_fit, included)
  final_fit <- nm_est(
    final_model, data, method = method, start = final_start,
    control = c(list(compute_inference = length(included) > 0L), control), ...
  )
  forest <- .nm_scm_forest_table(final_model, final_fit, model, included)

  structure(
    list(
      steps = do.call(rbind, steps),
      retained = included,
      forest = forest,
      base_ofv = base_ofv,
      final_ofv = as.numeric(final_fit$objective),
      final_model = final_model,
      final_fit = final_fit,
      forward_alpha = forward_alpha,
      backward_alpha = backward_alpha
    ),
    class = "nm_scm"
  )
}

#' @keywords internal
.nm_scm_normalize_candidates <- function(candidates) {
  if (is.data.frame(candidates)) {
    candidates <- lapply(seq_len(nrow(candidates)), function(i) {
      as.list(candidates[i, , drop = FALSE])
    })
  }
  lapply(candidates, function(c) {
    if (is.null(c$relationship) || is.na(c$relationship) ||
        !nzchar(c$relationship)) {
      c$relationship <- "exponential"
    }
    c$relationship <- match.arg(
      tolower(c$relationship), c("exponential", "power", "linear")
    )
    c$parameter <- as.character(c$parameter)
    c$covariate <- as.character(c$covariate)
    c
  })
}

#' @keywords internal
.nm_scm_eff_eq <- function(a, b) {
  identical(a$parameter, b$parameter) &&
    identical(a$covariate, b$covariate) &&
    identical(a$relationship, b$relationship)
}

#' @keywords internal
.nm_scm_in <- function(cand, lst) {
  any(vapply(lst, function(x) .nm_scm_eff_eq(x, cand), logical(1L)))
}

#' @keywords internal
.nm_scm_covariate_refs <- function(dat, candidates) {
  covs <- unique(vapply(candidates, function(c) c$covariate, character(1L)))
  refs <- list()
  for (cv in covs) {
    if (!cv %in% names(dat)) {
      .nm_stop("Covariate column not found in data: ", cv)
    }
    # Reference = median of the per-subject baseline covariate value.
    ids <- .nm_subject_ids(dat)
    vals <- vapply(ids, function(id) {
      .nm_cov_baseline_value(dat[dat$ID == id][[cv]], cv)
    }, numeric(1))
    refs[[cv]] <- stats::median(vals, na.rm = TRUE)
  }
  refs
}

#' @keywords internal
.nm_scm_build_model <- function(model, effects, refs) {
  m <- model
  for (eff in effects) {
    m <- .nm_scm_add_effect(m, eff$parameter, eff$covariate,
                            eff$relationship, refs[[eff$covariate]])
  }
  m
}

#' @keywords internal
.nm_scm_add_effect <- function(model, parameter, covariate, relationship, ref) {
  lines <- .nm_split_lines(model$PRED)
  lhs <- toupper(trimws(vapply(strsplit(lines, "=", fixed = TRUE),
                               function(x) x[[1L]], character(1L))))
  hit <- which(lhs == toupper(parameter))
  if (length(hit) == 0L) {
    .nm_stop("Parameter '", parameter, "' not found in $PRED for covariate effect.")
  }
  i <- hit[[1L]]
  k <- nrow(model$THETAS) + 1L
  refc <- formatC(ref, digits = 12L, format = "g")
  factor <- switch(
    relationship,
    exponential = sprintf("exp(THETA(%d) * (%s - %s))", k, covariate, refc),
    power = sprintf("exp(THETA(%d) * log(%s / %s))", k, covariate, refc),
    linear = sprintf("(1 + THETA(%d) * (%s - %s))", k, covariate, refc)
  )
  parts <- strsplit(lines[i], "=", fixed = TRUE)[[1L]]
  lhs_txt <- parts[[1L]]
  rhs_txt <- trimws(paste(parts[-1L], collapse = "="))
  lines[i] <- sprintf("%s= (%s) * %s", lhs_txt, rhs_txt, factor)
  model$PRED <- paste(lines, collapse = "\n")

  th <- model$THETAS
  if (!"Lower" %in% names(th)) th$Lower <- NA_real_
  if (!"Upper" %in% names(th)) th$Upper <- NA_real_
  new <- th[1L, , drop = FALSE]
  new[] <- lapply(new, function(x) if (is.character(x)) NA_character_ else NA)
  new$THETA <- k
  new$Value <- 0
  new$Lower <- -100
  new$Upper <- 100
  if ("FIX" %in% names(new)) new$FIX <- FALSE
  if ("Label" %in% names(new)) new$Label <- paste0(parameter, "_", covariate)
  model$THETAS <- rbind(th, new)
  model$COVARIATES <- unique(c(as.character(model$COVARIATES), covariate))
  if (!covariate %in% model$INPUT) {
    model$INPUT <- c(model$INPUT, covariate)
  }
  model
}

#' @keywords internal
.nm_scm_forest_table <- function(final_model, final_fit, base_model, effects) {
  if (length(effects) == 0L) {
    return(data.frame(
      parameter = character(), covariate = character(),
      relationship = character(),       estimate = numeric(),
      se = numeric(), lower = numeric(), upper = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  n_base_th <- nrow(base_model$THETAS)
  se_all <- tryCatch(.nm_scm_theta_se(final_fit), error = function(e) NULL)
  rows <- vector("list", length(effects))
  for (j in seq_along(effects)) {
    k <- n_base_th + j
    est <- as.numeric(final_fit$theta[k])
    se <- if (!is.null(se_all) && length(se_all) >= k) se_all[k] else NA_real_
    lower <- if (is.finite(se)) est - 1.96 * se else NA_real_
    upper <- if (is.finite(se)) est + 1.96 * se else NA_real_
    rows[[j]] <- data.frame(
      parameter = effects[[j]]$parameter,
      covariate = effects[[j]]$covariate,
      relationship = effects[[j]]$relationship,
      estimate = est,
      se = se,
      lower = lower,
      upper = upper,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

#' @keywords internal
.nm_scm_theta_se <- function(fit) {
  n_th <- length(fit$theta)
  # nm_est stores THETA/OMEGA/SIGMA standard errors (labelled) in fit$par_se;
  # the leading n_th entries are the THETA SEs. Fall back to other conventions.
  se <- fit$par_se %||% fit$standard_errors %||% fit$se
  if (!is.null(se) && length(se) >= n_th) {
    return(as.numeric(se)[seq_len(n_th)])
  }
  cov <- fit$covariance %||% fit$cov
  if (is.list(cov) && !is.matrix(cov)) {
    cov <- cov[[fit$covariance_method %||% names(cov)[1L]]]
  }
  if (!is.null(cov) && is.matrix(cov) && nrow(cov) >= n_th) {
    return(sqrt(pmax(diag(cov)[seq_len(n_th)], 0)))
  }
  rep(NA_real_, n_th)
}

#' @rdname nm_scm
#' @method print nm_scm
#' @param x An \code{nm_scm} object.
#' @param ... Unused.
#' @export
print.nm_scm <- function(x, ...) {
  cat("Stepwise covariate modelling (SCM)\n")
  cat("  base OFV:", format(x$base_ofv, digits = 6),
      " final OFV:", format(x$final_ofv, digits = 6), "\n")
  cat("  retained effects:", length(x$retained), "\n")
  print(x$steps)
  if (nrow(x$forest) > 0L) {
    cat("\nRetained effect estimates:\n")
    print(x$forest)
  }
  invisible(x)
}

#' Forest plot of retained SCM covariate effects
#'
#' @param scm An \code{nm_scm} object from \code{\link{nm_scm}}.
#' @param title Plot title.
#' @return A \code{ggplot} object (requires the suggested \pkg{ggplot2}).
#' @examples
#' \dontrun{
#' nm_forest_plot(nm_scm(model, data, candidates))
#' }
#' @export
nm_forest_plot <- function(scm, title = "SCM covariate effects") {
  if (!inherits(scm, "nm_scm")) {
    .nm_stop("scm must be an nm_scm object.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    .nm_stop("Package 'ggplot2' is required for nm_forest_plot() (Suggests).")
  }
  df <- scm$forest
  if (nrow(df) == 0L) {
    .nm_stop("No retained covariate effects to plot.")
  }
  df$label <- paste0(df$parameter, " ~ ", df$covariate)
  ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$estimate, y = .data$label)
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "grey50") +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data$lower, xmax = .data$upper), height = 0.2
    ) +
    ggplot2::labs(x = "Effect estimate (95% CI)", y = NULL, title = title) +
    ggplot2::theme_bw()
}
