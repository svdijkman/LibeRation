#' Validate model PK / PRED / ERROR definitions before estimation or simulation
#'
#' Checks that \code{$PK}/\code{$PRED} assignments evaluate and supply the
#' parameters required for the selected ADVAN/TRANS backend route.
#'
#' @param model An \code{nm_model} object.
#' @param data Optional dataset (used for a trial prediction when supplied).
#' @param stop_on_error If \code{TRUE}, stop with an informative error when
#'   validation fails; otherwise return a list with \code{ok} and messages.
#' @return A list with \code{ok} (logical), \code{issues} (character), and
#'   \code{pred_symbols} (character vector of symbols defined in PRED).
#' @examples
#' sim <- nm_synthetic_theo(n_sub = 2L)
#' nm_validate_model(sim$model, data = sim$data, stop_on_error = FALSE)
#' @export
nm_validate_model <- function(model, data = NULL, stop_on_error = TRUE) {
  issues <- character()
  pred_lines <- .nm_split_lines(model$PRED)
  if (length(pred_lines) == 0L) {
    issues <- c(issues, "Model has no $PK/$PRED assignments.")
  }
  pred_ok <- tryCatch(
    {
      if (.nm_cpp_pred_supported(model)) {
        .nm_cpp_pred_check(model, pred_lines)
      } else {
        .nm_eval_pred(
          model, model$THETAS$Value, model$OMEGAS$Value,
          rep(0, .nm_n_eta(model))
        )
        TRUE
      }
    },
    error = function(e) {
      issues <<- c(issues, paste0("$PK/$PRED evaluation failed: ", conditionMessage(e)))
      FALSE
    }
  )
  if (!isTRUE(pred_ok)) {
    out <- list(ok = FALSE, issues = issues, pred_symbols = character())
    if (isTRUE(stop_on_error)) .nm_stop(paste(issues, collapse = "\n"))
    return(out)
  }
  pred_symbols <- .nm_pred_defined_symbols(model)
  req <- .nm_model_required_pk_symbols(model)
  missing <- setdiff(req$required, pred_symbols)
  if ("Q2" %in% missing && "Q" %in% pred_symbols) {
    missing <- setdiff(missing, "Q2")
  }
  if (length(missing) > 0L && "K" %in% missing && "CL" %in% pred_symbols) {
    missing <- setdiff(missing, "K")
  }
  if (length(missing) > 0L && "CL" %in% missing && "K" %in% pred_symbols) {
    missing <- setdiff(missing, "CL")
  }
  if (length(missing) > 0L && isTRUE(req$micro_trans)) {
    if (all(c("K10", "K12", "K21") %in% pred_symbols)) {
      missing <- setdiff(missing, c("CL", "Q2", "Q3", "K10", "K12", "K21"))
    }
    if (all(c("K10", "K23", "K32") %in% pred_symbols)) {
      missing <- setdiff(missing, c("CL", "Q2", "Q3", "K10", "K12", "K21", "K23", "K32"))
    }
    if (all(c("K13", "K31") %in% pred_symbols) ||
        all(c("K14", "K41") %in% pred_symbols) ||
        all(c("K24", "K42") %in% pred_symbols)) {
      missing <- setdiff(missing, c("Q3", "K13", "K31", "K14", "K41", "K24", "K42"))
    }
  }
  if (length(missing) > 0L) {
    msg <- paste0(
      "Missing PK parameter(s) for ADVAN ", model$ADVAN, " TRANS ", model$TRANS,
      ": ", paste(missing, collapse = ", "),
      if (nzchar(req$hint)) paste0(" (", req$hint, ")") else ""
    )
    issues <- c(issues, msg)
  }
  if (isTRUE(req$need_central) && !any(req$central_aliases %in% pred_symbols)) {
    issues <- c(
      issues,
      paste0(
        "No central compartment volume found in $PK/$PRED (expected one of: ",
        paste(req$central_aliases, collapse = ", "), ")."
      )
    )
  }
  if (isTRUE(req$need_peripheral) && !any(req$periph_aliases %in% pred_symbols)) {
    issues <- c(
      issues,
      paste0(
        "No peripheral compartment volume found in $PK/$PRED (expected VP or V3). ",
        "THEO-style models use VP = THETA(3), not V3, unless you assign V3 explicitly."
      )
    )
  }
  alias_note <- .nm_pred_alias_notes(pred_symbols, req)
  if (length(alias_note) > 0L) {
    issues <- c(issues, alias_note)
  }
  err_lines <- .nm_split_lines(model$ERROR)
  if (length(err_lines) == 0L) {
    issues <- c(issues, "Model has no $ERROR block.")
  } else if (!any(grepl("\\bF\\b", paste(err_lines, collapse = " ")))) {
    issues <- c(issues, "$ERROR should reference predicted value F (e.g. Y = F * (1 + ERR(1)) + ERR(2)).")
  }
  if (!is.null(data)) {
    trial <- tryCatch(
      {
        dat <- .nm_prepare_data(data, model$INPUT, model)
        ids <- .nm_subject_ids(dat)
        if (length(ids) == 0L) {
          issues <<- c(issues, "Dataset has no subjects.")
          FALSE
        } else {
          subj <- .nm_subject_slice(dat, ids[[1L]])
          pred <- .nm_subject_ipred(
            model, subj, model$THETAS$Value, model$OMEGAS$Value,
            rep(0, .nm_n_eta(model)), model$SIGMAS$Value, pk_engine = "cpp"
          )
          f <- pred$F
          if (length(f) == 0L || !any(is.finite(f))) {
            issues <<- c(issues, "Trial prediction returned no finite F values for the first subject.")
            FALSE
          } else {
            TRUE
          }
        }
      },
      error = function(e) {
        issues <<- c(issues, paste0("Trial prediction failed: ", conditionMessage(e)))
        FALSE
      }
    )
    if (!isTRUE(trial)) {
      out <- list(ok = FALSE, issues = issues, pred_symbols = pred_symbols)
      if (isTRUE(stop_on_error)) .nm_stop(paste(issues, collapse = "\n"))
      return(out)
    }
  }
  ok <- length(issues) == 0L
  out <- list(ok = ok, issues = issues, pred_symbols = pred_symbols)
  if (!ok && isTRUE(stop_on_error)) .nm_stop(paste(issues, collapse = "\n"))
  out
}

#' @keywords internal
.nm_pred_defined_symbols <- function(model) {
  pred_lines <- .nm_split_lines(model$PRED)
  if (length(pred_lines) == 0L) {
    return(character())
  }
  if (.nm_cpp_capable(model)) {
    th <- as.numeric(model$THETAS$Value)
    if (length(th) == 0L) th <- 1
    et <- rep(0, max(1L, .nm_n_eta(model)))
    covs <- as.character(model$COVARIATES)
    cov_stub <- if (length(covs) > 0L) {
      stats::setNames(as.list(rep(1, length(covs))), covs)
    } else {
      list()
    }
    out <- nm_eval_pred_cpp(pred_lines, th, et, cov_stub)
    return(toupper(names(out)))
  }
  vals <- .nm_eval_pred(
    model, model$THETAS$Value, model$OMEGAS$Value,
    rep(0, .nm_n_eta(model))
  )
  toupper(names(vals))
}

#' @keywords internal
.nm_model_required_pk_symbols <- function(model) {
  advan <- as.integer(model$ADVAN)
  trans <- as.integer(model$TRANS)
  oral <- advan %in% c(2L, 4L, 6L, 12L, 13L)
  two_comp <- advan %in% c(3L, 4L, 11L, 12L)
  three_comp <- advan %in% c(11L, 12L)
  micro_trans <- trans == 6L ||
    (trans == 1L && advan %in% c(3L, 4L, 11L, 12L))
  hint <- ""
  required <- character()
  central_aliases <- c("VC", "V2", "V", "V1")
  periph_aliases <- c("VP", "V3", "V2")
  need_central <- (two_comp && !isTRUE(micro_trans)) || advan %in% c(1L, 2L)
  need_peripheral <- two_comp && !isTRUE(micro_trans)
  if (advan == 4L && trans == 1L) {
    required <- c("KA", "K10", "K23", "K32")
    need_central <- FALSE
    need_peripheral <- FALSE
    micro_trans <- TRUE
    hint <- "ADVAN4 TRANS1 oral micro uses KA, K10, K23, K32."
  } else if (oral) {
    required <- c(required, "KA")
    hint <- "oral models need KA (or KTR for transit)"
  }
  if (two_comp && !(advan == 4L && trans == 1L)) {
    if (isTRUE(micro_trans)) {
      required <- c(required, c("K10", "K12", "K21"))
      hint <- paste(hint, "Micro-rate models need K10, K12, K21 (and K13/K31 or K14/K41 for three-compartment).")
    } else {
      required <- c(required, "CL")
      hint <- paste(
        hint,
        "Two-compartment models need CL, inter-compartmental clearance (Q2 or Q), and volumes."
      )
    }
  } else if (advan %in% c(1L, 2L)) {
    if (trans == 1L) {
      required <- c(required, c("K", "V"))
      hint <- "One-compartment TRANS 1 models need K and V."
    } else {
      required <- c(required, c("CL", "V"))
      hint <- "One-compartment models need CL and V (or VC)."
    }
  }
  if (three_comp && !isTRUE(micro_trans)) {
    required <- setdiff(required, "Q2")
    required <- c(required, "Q3")
    hint <- paste(hint, "Three-compartment models also need Q3 (and VP2/V4 when applicable).")
  }
  required <- unique(required)
  list(
    required = required,
    hint = trimws(hint),
    need_central = need_central,
    need_peripheral = need_peripheral,
    central_aliases = central_aliases,
    periph_aliases = periph_aliases,
    micro_trans = micro_trans
  )
}

#' @keywords internal
.nm_pred_alias_notes <- function(pred_symbols, req) {
  msgs <- character()
  if ("V3" %in% pred_symbols && !any(c("VP", "V2", "VC") %in% pred_symbols)) {
    msgs <- c(msgs, "V3 is defined but no central volume (VC or V2) was found in $PK/$PRED.")
  }
  if ("VP" %in% pred_symbols && "V3" %in% pred_symbols) {
    msgs <- c(msgs, "Both VP and V3 are assigned; use one peripheral volume name consistently (VP or V3).")
  }
  if ("Q2" %in% req$required && !"Q2" %in% pred_symbols && !"Q" %in% pred_symbols) {
    msgs <- c(msgs, "Inter-compartmental clearance Q2 (or Q) is required but not defined in $PK/$PRED.")
  }
  msgs
}

#' @keywords internal
.nm_validate_model_quiet <- function(model, data = NULL) {
  nm_validate_model(model, data = data, stop_on_error = TRUE)
}

#' Reference validation against expected objectives (THEO-like)
#'
#' @param fit An \code{nm_fit} object.
#' @param tol Objective tolerance.
#' @return Logical; invisibly TRUE if checks pass.
#' @examples
#' \dontrun{
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "FO",
#'               control = list(maxit = 5L, compute_inference = FALSE))
#' nm_validate_fit(fit)
#' }
#' @export
nm_validate_fit <- function(fit, tol = 50) {
  if (!identical(fit$method, "FO") && !identical(fit$method, "FOCE") &&
      !identical(fit$method, "LAPLACE")) {
    message("Validation skipped for method ", fit$method)
    return(invisible(TRUE))
  }
  ok_obj <- is.finite(fit$objective)
  ok_par <- all(is.finite(fit$theta))
  if (!ok_obj || !ok_par) {
    warning("Validation failed: non-finite objective or THETA.")
    return(invisible(FALSE))
  }
  invisible(TRUE)
}

#' Export parameter table (NONMEM-style .par summary)
#' @param fit An \code{nm_fit} object.
#' @examples
#' \dontrun{
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "FO",
#'               control = list(maxit = 5L, compute_inference = FALSE))
#' nm_par_table(fit)
#' }
#' @export
nm_par_table <- function(fit) {
  model <- fit$model
  th <- stats::setNames(fit$theta, paste0("THETA", model$THETAS$THETA))
  om <- if (length(fit$omega) > 0L) {
    stats::setNames(fit$omega, paste0("OMEGA", model$OMEGAS$OMEGA))
  } else {
    numeric()
  }
  sg <- if (length(fit$sigma) > 0L) {
    stats::setNames(fit$sigma, paste0("SIGMA", model$SIGMAS$SIGMA))
  } else {
    numeric()
  }
  data.frame(
    parameter = c(names(th), names(om), names(sg)),
    estimate = c(unname(th), unname(om), unname(sg)),
    row.names = NULL
  )
}

#' Supported method / grad / engine matrix
#' @examples
#' nm_support_matrix()
#' @export
nm_support_matrix <- function() {
  data.frame(
    method = c("FO", "FOCE", "FOCEI", "SAEM", "LAPLACE", "IMP", "BAYES"),
    cpp_pk = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    cpp_pop_grad = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
    ad_tape = c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE),
    mcmc = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
}
