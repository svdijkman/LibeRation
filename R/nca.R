.nm_nca_stop <- function(...) stop(..., call. = FALSE)

.nm_nca_scalar <- function(value, group, default = NA_real_) {
  if (is.null(value)) return(default)
  if (is.character(value) && length(value) == 1L && value %in% names(group)) {
    available <- suppressWarnings(as.numeric(group[[value]]))
    available <- available[is.finite(available)]
    return(if (length(available)) available[[1L]] else default)
  }
  value <- suppressWarnings(as.numeric(value))
  if (!length(value) || !is.finite(value[[1L]])) default else value[[1L]]
}

.nm_nca_preprocess <- function(data, time, concentration,
                               duplicate = c("mean", "last", "error"),
                               blq = c("zero", "omit", "half_lloq"), lloq = NULL) {
  duplicate <- match.arg(duplicate); blq <- match.arg(blq)
  x <- suppressWarnings(as.numeric(data[[time]]))
  y <- suppressWarnings(as.numeric(data[[concentration]]))
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]; y <- y[keep]
  if (!length(x)) .nm_nca_stop("No finite time-concentration pairs remain after preprocessing.")
  if (blq == "omit") {
    keep <- y > 0; x <- x[keep]; y <- y[keep]
  } else if (blq == "zero") {
    y[y < 0] <- 0
  } else {
    limit <- suppressWarnings(as.numeric(lloq))
    if (length(limit) != 1L || !is.finite(limit) || limit <= 0) {
      .nm_nca_stop("`lloq` must be a positive scalar when `blq = \"half_lloq\"`.")
    }
    y[y <= 0] <- limit / 2
  }
  order <- order(x, seq_along(x)); x <- x[order]; y <- y[order]
  if (anyDuplicated(x)) {
    if (duplicate == "error") .nm_nca_stop("Duplicate sampling times require an explicit aggregation policy.")
    groups <- split(seq_along(x), x)
    x <- as.numeric(names(groups))
    y <- vapply(groups, function(index) {
      if (duplicate == "last") y[utils::tail(index, 1L)] else mean(y[index])
    }, numeric(1))
    order <- order(x); x <- x[order]; y <- y[order]
  }
  if (length(x) < 2L) .nm_nca_stop("At least two distinct sampling times are required for NCA.")
  data.frame(time = x, concentration = y)
}

.nm_nca_native_one <- function(profile, method, dose, tau, route,
                               terminal, partial_auc) {
  starts <- if (is.null(partial_auc)) numeric() else as.numeric(partial_auc$start)
  ends <- if (is.null(partial_auc)) numeric() else as.numeric(partial_auc$end)
  answer <- .liberation_nca_profile(
    profile$time, profile$concentration, method, dose, tau, route,
    as.integer(terminal %||% integer()), starts, ends
  )
  partial <- answer$partial_auc; answer$partial_auc <- NULL
  answer$terminal_indices <- I(list(answer$terminal_indices))
  if (length(partial)) {
    names(partial) <- as.character(partial_auc$name)
    for (name in names(partial)) answer[[paste0("PAUC_", make.names(name))]] <- partial[[name]]
  }
  as.data.frame(answer, check.names = FALSE, stringsAsFactors = FALSE)
}

.nm_nca_reference_one <- function(profile, method, dose, route, terminal,
                                  partial_auc, steady_state = FALSE,
                                  infusion_duration = 0, dose_unit = "mg",
                                  time_unit = "h", concentration_unit = "mg/L") {
  if (!requireNamespace("ncar", quietly = TRUE) ||
      !requireNamespace("NonCompart", quietly = TRUE)) {
    .nm_nca_stop("The optional `ncar` and `NonCompart` packages are required for the ncar reference backend.")
  }
  administration <- switch(route,
    oral =, extravascular =, ev = "Extravascular",
    bolus =, iv_bolus = "Bolus",
    infusion =, iv_infusion = "Infusion",
    "Extravascular"
  )
  result <- NonCompart::sNCA(
    profile$time, profile$concentration,
    dose = if (is.finite(dose)) dose else 0,
    adm = administration, dur = infusion_duration,
    doseUnit = dose_unit, timeUnit = time_unit, concUnit = concentration_unit,
    iAUC = if (is.null(partial_auc)) "" else data.frame(
      Name = partial_auc$name, Start = partial_auc$start, End = partial_auc$end
    ),
    down = if (method == "lin_up_log_down") "Log" else "Linear",
    SS = isTRUE(steady_state), UsePoints = terminal
  )
  values <- as.numeric(result)
  names(values) <- names(result)
  mapping <- c(
    CMAX = "CMAX", TMAX = "TMAX", CLST = "CLAST", CLSTP = "CLAST_PRED",
    TLST = "TLAST", AUCLST = "AUCLAST", AUMCLST = "AUMCLAST",
    LAMZ = "LAMBDA_Z", LAMZNPT = "LAMBDA_Z_N", LAMZLL = "LAMBDA_Z_LOWER",
    LAMZUL = "LAMBDA_Z_UPPER", R2 = "R2", R2ADJ = "R2_ADJ",
    LAMZHL = "HALF_LIFE", AUCIFO = "AUCINF_OBS", AUCPEO = "AUC_EXTRAP_PERCENT",
    AUMCIFO = "AUMCINF_OBS", CMAXD = "CMAX_DOSE_NORM", AUCIFOD = "AUCINF_DOSE_NORM",
    CLO = "CL", CLFO = "CL_F", VZO = "VZ", VZFO = "VZ_F",
    MRTIVIFO = "MRT", MRTEVIFO = "MRT"
  )
  output <- as.list(values)
  names(output) <- ifelse(names(values) %in% names(mapping), mapping[names(values)], names(values))
  output$N <- nrow(profile)
  output$terminal_indices <- I(list(attr(result, "UsedPoints") %||% terminal %||% integer()))
  as.data.frame(output, check.names = FALSE, stringsAsFactors = FALSE)
}

.nm_nca_compare <- function(native, reference) {
  metrics <- intersect(names(native), names(reference))
  metrics <- metrics[vapply(metrics, function(name) {
    is.numeric(native[[name]]) && is.numeric(reference[[name]])
  }, logical(1))]
  if (!length(metrics)) return(data.frame())
  do.call(rbind, lapply(metrics, function(metric) {
    observed <- as.numeric(native[[metric]])
    expected <- as.numeric(reference[[metric]])
    difference <- observed - expected
    data.frame(
      metric = metric, native = observed, ncar = expected,
      difference = difference, absolute_difference = abs(difference),
      relative_difference = difference / pmax(abs(expected), sqrt(.Machine$double.eps)),
      stringsAsFactors = FALSE
    )
  }))
}

#' Native noncompartmental pharmacokinetic analysis
#'
#' Calculates noncompartmental exposure and disposition summaries with a
#' native C++ implementation. The optional `ncar`/`NonCompart` route can be
#' selected as a reference backend, used automatically as a fallback, or run
#' alongside the native engine for validation.
#'
#' @param data A data frame containing one or more concentration-time profiles.
#' @param time,concentration Column names for time and concentration.
#' @param id Optional grouping columns. With `NULL`, the data form one profile.
#' @param dose Dose scalar or column name.
#' @param tau Dosing interval scalar or column name for steady-state summaries.
#' @param route Administration route: extravascular/oral, bolus, or infusion.
#' @param method Integration method: linear-up/log-down or linear trapezoidal.
#' @param terminal Optional terminal-phase row indices, or a named list by group.
#' @param partial_auc Optional data frame with `name`, `start`, and `end`.
#' @param engine Native engine, ncar reference backend, or automatic native with
#'   ncar fallback.
#' @param validate Compare native values with ncar when available.
#' @param duplicate Duplicate-time policy.
#' @param blq BLQ/non-positive concentration policy.
#' @param lloq LLOQ used by the half-LLOQ policy.
#' @param steady_state Passed to the ncar reference backend.
#' @param infusion_duration Infusion duration passed to the reference backend.
#' @param dose_unit,time_unit,concentration_unit Units passed to the ncar
#'   reference backend. Native results retain the input units.
#' @return An `nm_nca` object containing a wide results table, provenance, and
#'   optional validation differences.
#' @export
nm_nca <- function(data, time = "TIME", concentration = "DV", id = "ID",
                   dose = NULL, tau = NULL,
                   route = c("extravascular", "oral", "bolus", "infusion"),
                   method = c("lin_up_log_down", "linear"), terminal = NULL,
                   partial_auc = NULL, engine = c("native", "auto", "ncar"),
                   validate = FALSE, duplicate = c("mean", "last", "error"),
                   blq = c("zero", "omit", "half_lloq"), lloq = NULL,
                   steady_state = FALSE, infusion_duration = 0,
                   dose_unit = "mg", time_unit = "h",
                   concentration_unit = "mg/L") {
  started <- proc.time()[[3L]]
  data <- as.data.frame(data); route <- match.arg(route); method <- match.arg(method)
  engine <- match.arg(engine); duplicate <- match.arg(duplicate); blq <- match.arg(blq)
  required <- c(time, concentration, id %||% character())
  missing <- setdiff(required, names(data))
  if (length(missing)) .nm_nca_stop("NCA data are missing column(s): ", paste(missing, collapse = ", "), ".")
  if (!is.null(partial_auc)) {
    partial_auc <- as.data.frame(partial_auc)
    if (!all(c("name", "start", "end") %in% names(partial_auc))) {
      .nm_nca_stop("`partial_auc` must contain `name`, `start`, and `end`.")
    }
    if (any(!is.finite(partial_auc$start)) || any(!is.finite(partial_auc$end)) ||
        any(partial_auc$end <= partial_auc$start)) .nm_nca_stop("Partial-AUC intervals are invalid.")
  }
  key <- if (is.null(id) || !length(id)) rep("profile", nrow(data)) else {
    interaction(data[id], drop = TRUE, lex.order = TRUE, sep = " | ")
  }
  groups <- split(seq_len(nrow(data)), key, drop = TRUE)
  run_one <- function(index, backend, terminal_override = NULL) {
    group <- data[index, , drop = FALSE]
    profile <- .nm_nca_preprocess(group, time, concentration, duplicate, blq, lloq)
    group_name <- as.character(key[index[[1L]]])
    selected_terminal <- terminal_override %||%
      if (is.list(terminal)) terminal[[group_name]] %||% NULL else terminal
    group_dose <- .nm_nca_scalar(dose, group)
    group_tau <- .nm_nca_scalar(tau, group)
    if (backend == "native") {
      .nm_nca_native_one(profile, method, group_dose, group_tau, route,
                         selected_terminal, partial_auc)
    } else {
      .nm_nca_reference_one(profile, method, group_dose, route,
                            selected_terminal, partial_auc, steady_state,
                            infusion_duration, dose_unit, time_unit,
                            concentration_unit)
    }
  }
  backend_used <- engine
  results <- lapply(groups, function(index) {
    if (engine == "auto") {
      tryCatch(run_one(index, "native"), error = function(native_error) {
        backend_used <<- "ncar/NonCompart fallback"
        run_one(index, "ncar")
      })
    } else run_one(index, engine)
  })
  keys <- lapply(groups, function(index) {
    if (is.null(id) || !length(id)) data.frame(PROFILE = "profile")
    else data[index[[1L]], id, drop = FALSE]
  })
  table <- do.call(rbind, Map(function(key_row, result) cbind(key_row, result), keys, results))
  rownames(table) <- NULL
  validation <- NULL
  if (isTRUE(validate) && engine != "ncar") {
    if (requireNamespace("ncar", quietly = TRUE) && requireNamespace("NonCompart", quietly = TRUE)) {
      reference <- lapply(seq_along(groups), function(i) {
        selected <- results[[i]]$terminal_indices[[1L]] %||% NULL
        run_one(groups[[i]], "ncar", selected)
      })
      validation <- do.call(rbind, Map(function(key_row, native, expected) {
        comparison <- .nm_nca_compare(native, expected)
        if (!nrow(comparison)) return(comparison)
        cbind(key_row[rep(1L, nrow(comparison)), , drop = FALSE], comparison)
      }, keys, results, reference))
      rownames(validation) <- NULL
    } else {
      validation <- structure(data.frame(), unavailable = "Install ncar and NonCompart to run reference validation.")
    }
  }
  structure(list(
    schema = "liberation.nca", version = 1L, results = table,
    method = method, route = route, backend = backend_used,
    validation = validation, profiles = length(groups),
    elapsed_seconds = proc.time()[[3L]] - started,
    provenance = list(package = "LibeRation", package_version = as.character(utils::packageVersion("LibeRation")),
                      native = engine != "ncar", reference = "ncar/NonCompart"),
    call = match.call()
  ), class = "nm_nca")
}

#' Print an NCA result
#'
#' @param x An `nm_nca` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.nm_nca <- function(x, ...) {
  cat("LibeRation noncompartmental analysis\n")
  cat("  profiles:", x$profiles, " backend:", x$backend,
      " method:", x$method, "\n")
  print(x$results, row.names = FALSE)
  invisible(x)
}
