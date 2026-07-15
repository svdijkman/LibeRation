.nm_stop <- function(..., call. = FALSE) stop(..., call. = call.)

.nm_numbered_names <- function(prefix, count) {
  count <- as.integer(count)
  if (count <= 0L) character() else paste0(prefix, seq_len(count))
}

.nm_parameter_names <- function(theta, sigma, omega) {
  c(
    .nm_numbered_names("THETA", length(theta)),
    .nm_numbered_names("SIGMA", length(sigma)),
    .nm_numbered_names("OMEGA", length(omega))
  )
}
`%||%` <- function(x, y) if (is.null(x)) y else x

.nm_parameter_table <- function(x, kind, required = FALSE) {
  index <- toupper(kind)
  if (is.null(x)) {
    if (required) .nm_stop(kind, " table is required.")
    return(data.frame(stats::setNames(list(integer()), index), Value = numeric(),
                      FIX = logical(), stringsAsFactors = FALSE))
  }
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!index %in% names(x) || !"Value" %in% names(x)) {
    .nm_stop(kind, " table requires columns `", index, "` and `Value`.")
  }
  x[[index]] <- as.integer(x[[index]])
  x$Value <- as.numeric(x$Value)
  if (!"FIX" %in% names(x)) x$FIX <- FALSE
  x$FIX <- as.logical(x$FIX)
  if (anyNA(x[[index]]) || any(x[[index]] < 1L) || anyDuplicated(x[[index]])) {
    .nm_stop(kind, " indices must be unique positive integers.")
  }
  if (any(!is.finite(x$Value))) .nm_stop(kind, " values must be finite.")
  if (identical(index, "THETA")) {
    default_bounds <- function(initial) {
      if (initial > 0) return(c(initial / 1000, initial * 1000))
      if (initial < 0) return(c(initial * 1000, initial / 1000))
      c(-1000, 1000)
    }
    if (!"LOWER" %in% names(x)) x$LOWER <- NA_real_
    if (!"UPPER" %in% names(x)) x$UPPER <- NA_real_
    x$LOWER <- suppressWarnings(as.numeric(x$LOWER))
    x$UPPER <- suppressWarnings(as.numeric(x$UPPER))
    for (row in seq_len(nrow(x))) {
      fallback <- default_bounds(x$Value[[row]])
      lower <- x$LOWER[[row]]
      upper <- x$UPPER[[row]]
      if (!is.finite(lower)) lower <- fallback[[1L]]
      if (!is.finite(upper)) upper <- fallback[[2L]]
      if (!is.finite(lower) || !is.finite(upper) || lower >= upper ||
          x$Value[[row]] < lower || x$Value[[row]] > upper) {
        lower <- fallback[[1L]]
        upper <- fallback[[2L]]
      }
      x$LOWER[[row]] <- lower
      x$UPPER[[row]] <- upper
    }
  }
  x[order(x[[index]]), , drop = FALSE]
}

.nm_omega_table <- function(x) {
  table <- .nm_parameter_table(x, "OMEGA")
  if (!nrow(table)) {
    table$ROW <- table$COL <- integer()
    return(table)
  }
  has_row <- "ROW" %in% names(table)
  has_col <- "COL" %in% names(table)
  if (xor(has_row, has_col)) .nm_stop("OMEGAS must supply both `ROW` and `COL`.")
  if (!has_row) {
    table$ROW <- table$OMEGA
    table$COL <- table$OMEGA
  }
  table$ROW <- as.integer(table$ROW)
  table$COL <- as.integer(table$COL)
  if (anyNA(table$ROW) || anyNA(table$COL) || any(table$ROW < 1L) ||
      any(table$COL < 1L)) {
    .nm_stop("OMEGA ROW/COL indices must be positive integers.")
  }
  swap <- table$ROW < table$COL
  temporary <- table$ROW[swap]
  table$ROW[swap] <- table$COL[swap]
  table$COL[swap] <- temporary
  if (anyDuplicated(paste(table$ROW, table$COL, sep = ":"))) {
    .nm_stop("OMEGA covariance positions must be unique.")
  }
  n_eta <- max(table$ROW, table$COL)
  diagonal <- match(paste(seq_len(n_eta), seq_len(n_eta), sep = ":"),
                    paste(table$ROW, table$COL, sep = ":"))
  if (anyNA(diagonal) || any(table$Value[diagonal] < 0) ||
      any(table$Value[diagonal] == 0 & !table$FIX[diagonal])) {
    .nm_stop("OMEGAS requires a non-negative diagonal value for every ETA; zero variances must be FIXed.")
  }
  matrix <- matrix(0, n_eta, n_eta)
  for (i in seq_len(nrow(table))) {
    matrix[table$ROW[[i]], table$COL[[i]]] <- table$Value[[i]]
    matrix[table$COL[[i]], table$ROW[[i]]] <- table$Value[[i]]
  }
  eigenvalues <- eigen(matrix, symmetric = TRUE, only.values = TRUE)$values
  if (min(eigenvalues) < -1e-12 ||
      (min(eigenvalues) <= 0 && any(table$ROW != table$COL))) {
    .nm_stop("Initial correlated OMEGA covariance matrix must be positive definite.")
  }
  if (any(table$ROW != table$COL) && nrow(table) != n_eta * (n_eta + 1L) / 2L) {
    .nm_stop("A correlated OMEGA must supply the complete lower triangle.")
  }
  table
}

.nm_n_eta <- function(omega) {
  if (!nrow(omega)) 0L else max(as.integer(omega$ROW), as.integer(omega$COL))
}

.nm_omega_matrix <- function(model, values = model$OMEGAS$Value) {
  n_eta <- .nm_n_eta(model$OMEGAS)
  matrix <- matrix(0, n_eta, n_eta)
  for (i in seq_len(nrow(model$OMEGAS))) {
    row <- model$OMEGAS$ROW[[i]]
    column <- model$OMEGAS$COL[[i]]
    matrix[row, column] <- matrix[column, row] <- values[[i]]
  }
  matrix
}

.nm_error_type <- function(code, requested = "auto") {
  requested <- match.arg(requested, c("auto", "none", "additive", "proportional", "combined", "exponential"))
  if (requested != "auto") return(requested)
  compact <- toupper(gsub("[[:space:]]+", "", code %||% ""))
  if (!grepl("ERR\\(", compact)) return("none")
  has1 <- grepl("ERR\\(1\\)", compact)
  has2 <- grepl("ERR\\(2\\)", compact)
  if (has1 && has2) return("combined")
  if (grepl("EXP\\(ERR\\(1\\)\\)", compact)) return("exponential")
  if (grepl("F\\*.*ERR\\(1\\)|ERR\\(1\\).*\\*F", compact)) return("proportional")
  "additive"
}

.nm_known_graph <- function(advan) {
  advan <- as.integer(advan)
  names <- switch(
    as.character(advan),
    `1` = "CENTRAL",
    `2` = c("DEPOT", "CENTRAL"),
    `3` = c("CENTRAL", "PERIPHERAL1"),
    `4` = c("DEPOT", "CENTRAL", "PERIPHERAL1"),
    `11` = c("CENTRAL", "PERIPHERAL1", "PERIPHERAL2"),
    `12` = c("DEPOT", "CENTRAL", "PERIPHERAL1", "PERIPHERAL2"),
    character()
  )
  list(
    compartments = data.frame(
      id = seq_along(names), name = names, state = paste0("A", seq_along(names)),
      stringsAsFactors = FALSE
    ),
    source = if (length(names)) paste0("ADVAN", advan) else "user"
  )
}

.nm_compile_pred_ir <- function(pred, n_theta, n_eta, covariates) {
  declared <- c(
    if (n_theta > 0L) paste0("THETA_", seq_len(n_theta)) else character(),
    if (n_eta > 0L) paste0("ETA_", seq_len(n_eta)) else character(),
    as.character(covariates %||% character())
  )
  LibeRtAD::ad_ir(pred, inputs = declared)
}

.nm_rewrite_ode_indexing <- function(code) {
  code <- paste(code, collapse = "\n")
  for (name in c("A", "DADT")) {
    pattern <- paste0("\\b", name, "\\s*\\(\\s*([0-9]+)\\s*\\)")
    code <- gsub(pattern, paste0(name, "_\\1"), code, perl = TRUE)
  }
  code
}

.nm_compile_des_ir <- function(des, pred_ir, n_theta, n_eta, covariates) {
  rewritten <- .nm_rewrite_ode_indexing(des)
  parsed <- tryCatch(
    parse(text = rewritten),
    error = function(e) .nm_stop("Unable to parse DES block: ", conditionMessage(e))
  )
  lhs <- vapply(parsed, function(expr) {
    if (is.call(expr) && as.character(expr[[1L]]) %in% c("<-", "=") &&
        length(expr) == 3L && is.symbol(expr[[2L]])) {
      as.character(expr[[2L]])
    } else ""
  }, character(1))
  derivative_names <- grep("^DADT_[0-9]+$", lhs, value = TRUE)
  derivative_index <- suppressWarnings(as.integer(sub("^DADT_", "", derivative_names)))
  if (!length(derivative_index) || anyNA(derivative_index) ||
      !identical(sort(unique(derivative_index)), seq_len(max(derivative_index)))) {
    .nm_stop("DES must assign each derivative DADT(1) through DADT(n) exactly once.")
  }
  if (anyDuplicated(derivative_index)) {
    .nm_stop("Each DADT(i) derivative may only be assigned once in DES.")
  }
  n_state <- max(derivative_index)
  declared <- unique(c(
    pred_ir$output_names,
    paste0("A_", seq_len(n_state)), "T",
    if (n_theta > 0L) paste0("THETA_", seq_len(n_theta)) else character(),
    if (n_eta > 0L) paste0("ETA_", seq_len(n_eta)) else character(),
    as.character(covariates %||% character())
  ))
  list(
    ir = LibeRtAD::ad_ir(
      rewritten, inputs = declared, outputs = paste0("DADT_", seq_len(n_state))
    ),
    n_state = n_state
  )
}

.nm_ode_graph <- function(n_state) {
  list(
    compartments = data.frame(
      id = seq_len(n_state), name = paste0("COMPARTMENT", seq_len(n_state)),
      state = paste0("A", seq_len(n_state)), stringsAsFactors = FALSE
    ),
    source = "DES"
  )
}

.nm_ode_control <- function(control = NULL, advan) {
  defaults <- list(
    rtol = if (advan == 13L) 1e-7 else 1e-8,
    atol = if (advan == 13L) 1e-9 else 1e-10,
    max_steps = 100000L,
    initial_step = 0
  )
  if (!is.null(control)) {
    if (!is.list(control) || is.null(names(control))) {
      .nm_stop("ODE_CONTROL must be a named list.")
    }
    unknown <- setdiff(names(control), names(defaults))
    if (length(unknown)) .nm_stop("Unknown ODE_CONTROL setting(s): ", paste(unknown, collapse = ", "), ".")
    defaults[names(control)] <- control
  }
  defaults$rtol <- as.numeric(defaults$rtol)
  defaults$atol <- as.numeric(defaults$atol)
  defaults$max_steps <- as.integer(defaults$max_steps)
  defaults$initial_step <- as.numeric(defaults$initial_step)
  if (length(defaults$rtol) != 1L || !is.finite(defaults$rtol) || defaults$rtol <= 0 ||
      length(defaults$atol) != 1L || !is.finite(defaults$atol) || defaults$atol <= 0 ||
      length(defaults$max_steps) != 1L || is.na(defaults$max_steps) || defaults$max_steps < 1L ||
      length(defaults$initial_step) != 1L || !is.finite(defaults$initial_step) || defaults$initial_step < 0) {
    .nm_stop("ODE_CONTROL requires positive finite rtol/atol, max_steps >= 1, and initial_step >= 0.")
  }
  defaults
}

.nm_matrix_graph_spec <- function(graph) {
  if (!inherits(graph, "nm_matrix_model")) return(NULL)
  compartments <- graph$compartments
  flows <- graph$flows
  ids <- as.integer(compartments$id)
  position <- match(as.integer(flows$from), ids)
  target_id <- suppressWarnings(as.integer(flows$to))
  target_id[is.na(target_id)] <- 0L
  target <- ifelse(target_id == 0L, 0L, match(target_id, ids))
  volume_parameter <- if ("volume_parameter" %in% names(flows)) {
    as.character(flows$volume_parameter)
  } else {
    rep("", nrow(flows))
  }
  compartment_volume <- if ("volume_parameter" %in% names(compartments)) {
    as.character(compartments$volume_parameter)
  } else rep("", nrow(compartments))
  missing_volume <- !nzchar(volume_parameter)
  volume_parameter[missing_volume] <- compartment_volume[position[missing_volume]]
  scale_parameter <- if ("scale_parameter" %in% names(compartments)) {
    as.character(compartments$scale_parameter)
  } else rep("", nrow(compartments))
  list(
    names = as.character(compartments$name),
    scale_parameter = scale_parameter,
    from = as.integer(position),
    to = as.integer(target),
    type = tolower(as.character(flows$type)),
    parameter = as.character(flows$parameter),
    volume_parameter = volume_parameter
  )
}

#' Define a NONMEM-style pharmacometric model
#'
#' The established LibeRation `PRED` and `ERROR` strings are retained. They are
#' compiled to C++ expression IR; they are not evaluated inside the numerical
#' loop by R.
#'
#' @param INPUT Required dataset column names.
#' @param ADVAN ADVAN number.
#' @param TRANS TRANS parameterization number.
#' @param SS Model-level steady-state default.
#' @param DOSECMP Default dosing compartment.
#' @param OBSCMP Default observation compartment.
#' @param PRED Parameter/model assignment code.
#' @param ERROR Residual-error assignment code.
#' @param DES ODE derivative code for ADVAN6/13.
#' @param THETAS,OMEGAS,SIGMAS Parameter tables.
#' @param COVARIATES Dataset covariates exposed to `PRED`.
#' @param USE_ODE Whether an ODE solver is explicitly requested.
#' @param ODE_CONTROL Named list with `rtol`, `atol`, `max_steps`, and optional
#'   `initial_step`. ADVAN6 uses adaptive Dormand-Prince 5(4); ADVAN13 uses an
#'   A-stable adaptive implicit trapezoidal method.
#' @param IOV Number of trailing inter-occasion ETAs.
#' @param LIK_CONFIG Reserved likelihood configuration.
#' @param SOLVER `auto`, `advan`, `matrix`, or `ode`.
#' @param ERROR_TYPE Standard residual-error form, or `auto`.
#' @param GRAPH Optional semantic compartment graph.
#' @param LAYOUT Optional graphical layout, stored separately from `GRAPH`.
#' @param LANGUAGE Model source language (`R` or restricted `C++`).
#' @return A serializable `nm_model`.
#' @examples
#' model <- nm_model(
#'   INPUT = c("ID", "TIME", "EVID", "AMT", "CMT", "DV", "MDV"),
#'   ADVAN = 1,
#'   PRED = "CL=THETA(1)*exp(ETA(1)); V=THETA(2); S1=V",
#'   ERROR = "Y=F*(1+ERR(1))",
#'   THETAS = data.frame(THETA = 1:2, Value = c(2, 20)),
#'   OMEGAS = data.frame(OMEGA = 1, Value = 0.1),
#'   SIGMAS = data.frame(SIGMA = 1, Value = 0.1)
#' )
#' model
#' @export
nm_model <- function(INPUT,
                     ADVAN = 2L,
                     TRANS = 2L,
                     SS = 0L,
                     DOSECMP = 1L,
                     OBSCMP = 1L,
                     PRED = "",
                     ERROR = "Y = F",
                     DES = "",
                     THETAS,
                     OMEGAS = NULL,
                     SIGMAS = NULL,
                     COVARIATES = NULL,
                     USE_ODE = FALSE,
                     ODE_CONTROL = NULL,
                     IOV = 0L,
                     LIK_CONFIG = NULL,
                     SOLVER = c("auto", "advan", "matrix", "ode"),
                     ERROR_TYPE = c("auto", "none", "additive", "proportional", "combined", "exponential"),
                     GRAPH = NULL,
                     LAYOUT = NULL,
                     LANGUAGE = c("R", "C++")) {
  advan <- as.integer(ADVAN)
  supported <- c(1L, 2L, 3L, 4L, 6L, 11L, 12L, 13L)
  if (length(advan) != 1L || is.na(advan) || !advan %in% supported) {
    .nm_stop("ADVAN must be one of: ", paste(supported, collapse = ", "), ".")
  }
  solver <- match.arg(SOLVER)
  language <- match.arg(LANGUAGE)
  error_type <- .nm_error_type(ERROR, match.arg(ERROR_TYPE))
  lik_config <- .nm_lik_config(LIK_CONFIG, error_type, as.integer(IOV))
  theta <- .nm_parameter_table(THETAS, "THETA", required = TRUE)
  omega <- .nm_omega_table(OMEGAS)
  sigma <- .nm_parameter_table(SIGMAS, "SIGMA")
  if (lik_config$omega == "full" && .nm_n_eta(omega) > 1L &&
      nrow(omega) != .nm_n_eta(omega) * (.nm_n_eta(omega) + 1L) / 2L) {
    .nm_stop("A full OMEGA requires the complete lower triangle.")
  }
  if (error_type != "none" && nrow(sigma) == 0L) {
    .nm_stop("Residual error type '", error_type, "' requires SIGMAS.")
  }
  if (advan %in% c(6L, 13L) && !nzchar(trimws(DES))) {
    .nm_stop("ADVAN", advan, " requires a DES block.")
  }
  if (!nzchar(trimws(PRED))) .nm_stop("PRED code must not be empty.")
  n_eta <- .nm_n_eta(omega)
  if (lik_config$iov > n_eta) .nm_stop("IOV cannot exceed the number of ETA definitions.")
  if (lik_config$iov > 0L && n_eta > lik_config$iov) {
    between <- n_eta - lik_config$iov
    cross_level <- (omega$ROW <= between & omega$COL > between) |
      (omega$COL <= between & omega$ROW > between)
    if (any(cross_level & omega$Value != 0)) {
      .nm_stop("Between-subject and inter-occasion OMEGA blocks cannot be correlated.")
    }
  }
  compiler_covariates <- unique(c(
    COVARIATES, if (!is.null(lik_config$mixtures)) "MIXNUM" else character()
  ))
  pred_ir <- .nm_compile_pred_ir(PRED, nrow(theta), n_eta, compiler_covariates)
  des_info <- if (advan %in% c(6L, 13L)) {
    .nm_compile_des_ir(DES, pred_ir, nrow(theta), n_eta, compiler_covariates)
  } else NULL
  ode_control <- .nm_ode_control(ODE_CONTROL, advan)
  graph <- GRAPH %||% if (is.null(des_info)) .nm_known_graph(advan) else .nm_ode_graph(des_info$n_state)
  structure(
    list(
      version = 1L,
      INPUT = unique(as.character(INPUT)),
      ADVAN = advan,
      TRANS = as.integer(TRANS),
      SS = as.integer(SS),
      DOSECMP = as.integer(DOSECMP),
      OBSCMP = as.integer(OBSCMP),
      PRED = paste(PRED, collapse = "\n"),
      ERROR = paste(ERROR, collapse = "\n"),
      DES = paste(DES, collapse = "\n"),
      THETAS = theta,
      OMEGAS = omega,
      SIGMAS = sigma,
      COVARIATES = unique(as.character(COVARIATES %||% character())),
      USE_ODE = isTRUE(USE_ODE) || advan %in% c(6L, 13L),
      ODE_CONTROL = ode_control,
      IOV = as.integer(IOV),
      LIK_CONFIG = lik_config,
      SOLVER = solver,
      ERROR_TYPE = error_type,
      GRAPH = graph,
      LAYOUT = LAYOUT,
      LANGUAGE = language,
      pred_ir = pred_ir,
      des_ir = des_info$ir %||% NULL,
      n_eta = n_eta,
      n_state = des_info$n_state %||% nrow(graph$compartments %||% data.frame())
    ),
    class = "nm_model"
  )
}

#' Define a model using the restricted C++ expression form
#'
#' The first compiler version accepts scalar C++-style assignment expressions
#' using the same `THETA(i)`/`ETA(i)` calls as the R form. Arbitrary headers,
#' allocation, I/O, and side effects are intentionally not accepted.
#'
#' @param CPP Parameter/model assignment code.
#' @param ... Other arguments passed to [nm_model()].
#' @export
nm_model_cpp <- function(CPP, ...) {
  nm_model(PRED = CPP, LANGUAGE = "C++", ...)
}

#' Define an arbitrary linear compartment graph
#'
#' @param compartments Data frame with stable `id` and `name` columns.
#' @param flows Data frame describing directed flows and their parameterization.
#' @param observations Observation definitions.
#' @param inputs Dosing/input definitions.
#' @param layout Optional diagram layout kept separate from semantics.
#' @export
nm_matrix_model <- function(compartments, flows, observations = NULL,
                            inputs = NULL, layout = NULL) {
  compartments <- as.data.frame(compartments, stringsAsFactors = FALSE)
  flows <- as.data.frame(flows, stringsAsFactors = FALSE)
  if (!all(c("id", "name") %in% names(compartments))) {
    .nm_stop("`compartments` requires `id` and `name` columns.")
  }
  if (!all(c("from", "to", "type", "parameter") %in% names(flows))) {
    .nm_stop("`flows` requires `from`, `to`, `type`, and `parameter` columns.")
  }
  compartments$id <- as.integer(compartments$id)
  compartments$name <- as.character(compartments$name)
  if (anyNA(compartments$id) || any(compartments$id < 1L) ||
      anyDuplicated(compartments$id) || anyDuplicated(compartments$name) ||
      any(!nzchar(compartments$name))) {
    .nm_stop("Compartment `id` and `name` values must be unique, non-missing, and non-empty.")
  }
  flows$from <- as.integer(flows$from)
  flows$to <- suppressWarnings(as.integer(flows$to))
  flows$to[is.na(flows$to)] <- 0L
  flows$type <- tolower(as.character(flows$type))
  flows$parameter <- as.character(flows$parameter)
  allowed <- c("rate", "elimination", "clearance")
  if (anyNA(flows$from) || any(!flows$from %in% compartments$id)) {
    .nm_stop("Every flow `from` value must name a compartment id.")
  }
  if (any(!flows$to %in% c(0L, compartments$id))) {
    .nm_stop("Every flow `to` value must be 0 (elimination) or a compartment id.")
  }
  if (any(flows$from == flows$to)) .nm_stop("A flow cannot return to its source compartment.")
  if (any(!flows$type %in% allowed)) {
    .nm_stop("Flow `type` must be one of: ", paste(allowed, collapse = ", "), ".")
  }
  if (any(!nzchar(flows$parameter))) .nm_stop("Every flow requires a parameter name.")
  if (!"volume_parameter" %in% names(compartments)) {
    compartments$volume_parameter <- ""
  }
  if (!"scale_parameter" %in% names(compartments)) {
    compartments$scale_parameter <- compartments$volume_parameter
  }
  if (!"volume_parameter" %in% names(flows)) flows$volume_parameter <- ""
  clearance <- flows$type == "clearance"
  inherited_volume <- as.character(compartments$volume_parameter[match(flows$from, compartments$id)])
  resolved_volume <- as.character(flows$volume_parameter)
  resolved_volume[!nzchar(resolved_volume)] <- inherited_volume[!nzchar(resolved_volume)]
  if (any(clearance & !nzchar(resolved_volume))) {
    .nm_stop("Clearance flows require `volume_parameter` on the flow or source compartment.")
  }
  structure(
    list(
      version = 1L,
      compartments = compartments,
      flows = flows,
      observations = observations,
      inputs = inputs,
      layout = layout
    ),
    class = "nm_matrix_model"
  )
}

#' @export
print.nm_model <- function(x, ...) {
  cat("LibeRation model\n")
  cat("  ADVAN", x$ADVAN, " TRANS", x$TRANS, " solver:", x$SOLVER, "\n")
  cat("  source:", x$LANGUAGE, "  THETA:", nrow(x$THETAS),
      " ETA:", x$n_eta, " OMEGA parameters:", nrow(x$OMEGAS),
      " SIGMA:", nrow(x$SIGMAS), "\n")
  cat("  graph compartments:", nrow(x$GRAPH$compartments %||% data.frame()),
      " layout:", if (is.null(x$LAYOUT)) "none" else "present", "\n")
  invisible(x)
}

#' Report implemented feature support
#' @export
nm_support_matrix <- function() {
  feature <- c("ADVAN1", "ADVAN2", "ADVAN3", "ADVAN4", "ADVAN11", "ADVAN12",
               "matrix exponential", "steady-state bolus", "steady-state infusion",
               "nonlinear ODE steady state", "ADVAN6", "ADVAN13",
               "FO", "FOCE", "FOCEI", "LAPLACE", "ITS", "IMP", "SAEM", "BAYES",
               "full OMEGA", "IOV", "M3/M4 BLQ", "finite mixtures", "priors",
               "AR1 residuals", "stochastic simulation", "VPC", "categorical VPC",
               "time-to-event VPC", "bootstrap", "profile likelihood", "SCM",
               "NONMEM control-stream round-trip", "local queue", "remote queue",
               "React workbench", "model versions", "THETA bounds",
               "parallel estimation", "parallel simulation", "iteration gradients")
  data.frame(feature = feature, status = rep("implemented", length(feature)),
             stringsAsFactors = FALSE)
}
