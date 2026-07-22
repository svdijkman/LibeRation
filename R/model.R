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
  requested <- match.arg(requested, c(
    "auto", "none", "additive", "proportional", "combined",
    "exponential", "power", "likelihood"
  ))
  if (requested != "auto") return(requested)
  compact <- toupper(gsub("[[:space:]]+", "", code %||% ""))
  if (grepl("\\b(?:LOGLIK|LIK)(?:<-|=)", compact, perl = TRUE)) {
    return("likelihood")
  }
  if (!grepl("ERR\\(", compact)) return("none")
  has1 <- grepl("ERR\\(1\\)", compact)
  has2 <- grepl("ERR\\(2\\)", compact)
  if (has1 && has2) return("combined")
  if (grepl("EXP\\(ERR\\(1\\)\\)", compact)) return("exponential")
  if (grepl("F\\*.*ERR\\(1\\)|ERR\\(1\\).*\\*F", compact)) return("proportional")
  "additive"
}

.nm_compile_error_ir <- function(error, pred_ir, n_theta, n_eta, input,
                                 covariates, likelihood_type, hmm_config = NULL,
                                 kalman_config = NULL) {
  if (!identical(likelihood_type, "likelihood")) return(NULL)
  # Let the IR discover only symbols that the likelihood actually uses. This
  # avoids treating unused character identifiers or unrelated dataset columns
  # as numeric runtime inputs while still resolving referenced $PRED outputs
  # in the C++ likelihood evaluator.
  ir <- LibeRtAD::ad_ir(error)
  if (!is.null(kalman_config)) {
    required <- unique(c(
      kalman_config$initial_mean, as.vector(kalman_config$initial_covariance),
      as.vector(kalman_config$transition),
      as.vector(kalman_config$process_covariance), kalman_config$observation,
      kalman_config$observation_variance
    ))
    missing <- setdiff(required, ir$output_names)
    if (length(missing)) {
      .nm_stop(
        "KALMAN_CONFIG references $ERROR assignment(s) that were not found: ",
        paste(missing, collapse = ", "), "."
      )
    }
    if (any(c("LOGLIK", "LIK", "Y") %in% ir$output_names)) {
      .nm_stop(
        "A Kalman $ERROR block supplies state-space component outputs; ",
        "do not also assign LIK, LOGLIK, or Y."
      )
    }
    if (any(grepl("^ERR_", ir$input_names))) {
      .nm_stop("A state-space likelihood cannot use ERR(); define its observation variance explicitly.")
    }
    return(list(ir = ir, output = NULL, scale = "kalman"))
  }
  if (!is.null(hmm_config)) {
    transition_type <- hmm_config$transition_type %||% "discrete"
    transition_names <- if (identical(transition_type, "continuous")) {
      generator <- hmm_config$generator
      unname(generator[row(generator) != col(generator)])
    } else {
      as.vector(hmm_config$transition)
    }
    required <- unique(c(
      hmm_config$initial, transition_names, hmm_config$emission
    ))
    missing <- setdiff(required, ir$output_names)
    if (length(missing)) {
      .nm_stop(
        "HMM_CONFIG references $ERROR assignment(s) that were not found: ",
        paste(missing, collapse = ", "), "."
      )
    }
    if (any(c("LOGLIK", "LIK", "Y") %in% ir$output_names)) {
      .nm_stop(
        "An HMM $ERROR block supplies initial, transition/rate, and emission outputs; ",
        "do not also assign LIK, LOGLIK, or Y."
      )
    }
    if (any(grepl("^ERR_", ir$input_names))) {
      .nm_stop("A hidden Markov likelihood cannot use ERR(); define emission likelihoods explicitly.")
    }
    return(list(ir = ir, output = NULL, scale = "hmm"))
  }
  candidates <- intersect(c("LOGLIK", "LIK", "Y"), ir$output_names)
  if (!length(candidates)) {
    .nm_stop(
      "A user-defined likelihood `$ERROR` block must assign `LOGLIK` ",
      "(log likelihood) or `LIK` (probability/density). `Y` is also accepted ",
      "when LIK_CONFIG explicitly selects `error = \"likelihood\"`."
    )
  }
  preferred <- intersect(c("LOGLIK", "LIK"), candidates)
  if (length(preferred) > 1L) {
    .nm_stop("Assign only one of `LOGLIK` or `LIK` in a user-defined likelihood.")
  }
  output <- if (length(preferred)) preferred[[1L]] else "Y"
  if (any(grepl("^ERR_", ir$input_names))) {
    .nm_stop("A user-defined likelihood cannot use ERR(); put all randomness in LIK/LOGLIK.")
  }
  list(
    ir = ir,
    output = output,
    scale = if (identical(output, "LOGLIK")) "log" else "likelihood"
  )
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

.nm_compile_des_ir <- function(des, pred_ir, n_theta, n_eta, covariates,
                               dde_config = NULL,
                               algebraic_variables = character()) {
  lagged <- .nm_rewrite_dde_lags(paste(des, collapse = "\n"), dde_config)
  rewritten <- .nm_rewrite_ode_indexing(lagged$code)
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
  if (length(lagged$lags)) {
    lag_states <- vapply(lagged$lags, `[[`, integer(1), "state")
    lag_delays <- vapply(lagged$lags, `[[`, character(1), "delay")
    if (any(lag_states < 1L | lag_states > n_state)) {
      .nm_stop("A DDE LAG() state index exceeds the DES state dimension.")
    }
    missing_delay <- setdiff(lag_delays, pred_ir$output_names)
    if (length(missing_delay)) {
      .nm_stop("DDE delay assignment(s) missing from $PK/$PRED: ",
               paste(missing_delay, collapse = ", "), ".")
    }
  }
  declared <- unique(c(
    pred_ir$output_names,
    paste0("A_", seq_len(n_state)), "T",
    vapply(lagged$lags, `[[`, character(1), "input"),
    as.character(algebraic_variables),
    if (n_theta > 0L) paste0("THETA_", seq_len(n_theta)) else character(),
    if (n_eta > 0L) paste0("ETA_", seq_len(n_eta)) else character(),
    as.character(covariates %||% character())
  ))
  list(
    ir = LibeRtAD::ad_ir(
      rewritten, inputs = declared, outputs = paste0("DADT_", seq_len(n_state))
    ),
    n_state = n_state,
    lags = lagged$lags
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

.nm_standard_output_catalog <- function(n_state, n_eta) {
  rows <- list(
    data.frame(
      name = c("PRED", "IPRED", "RES", "IRES", "WRES", "IWRES", "CWRES"),
      source = "engine",
      availability = c("all", "all", rep("estimation", 5L)),
      description = c(
        "Population prediction", "Individual prediction",
        "Population residual", "Individual residual",
        "Population weighted residual", "Individual weighted residual",
        "Conditional weighted residual"
      ), stringsAsFactors = FALSE
    )
  )
  if (n_eta > 0L) {
    rows[[length(rows) + 1L]] <- data.frame(
      name = paste0("ETA", seq_len(n_eta)), source = "random effect",
      availability = "all", description = paste("Individual random effect", seq_len(n_eta)),
      stringsAsFactors = FALSE
    )
  }
  if (n_state > 0L) {
    rows[[length(rows) + 1L]] <- data.frame(
      name = paste0("A", seq_len(n_state)), source = "compartment state",
      availability = "all", description = paste("Amount in compartment", seq_len(n_state)),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

.nm_output_catalog <- function(pred_ir, n_state, n_eta, input = character()) {
  standard <- .nm_standard_output_catalog(n_state, n_eta)
  assigned <- setdiff(
    as.character(pred_ir$output_names %||% character()),
    c(standard$name, grep("^[.]value[0-9]+$", pred_ir$output_names %||% character(), value = TRUE))
  )
  generated <- data.frame(
    name = assigned, source = "model assignment", availability = "all",
    description = if (length(assigned)) paste("Assigned in $PK/$PRED:", assigned) else character(),
    stringsAsFactors = FALSE
  )
  output <- rbind(standard, generated)
  output$collision <- output$name %in% as.character(input)
  output$selectable <- !output$collision
  rownames(output) <- NULL
  output
}

.nm_validate_outputs <- function(output, catalog) {
  output <- unique(as.character(output %||% character()))
  if (anyNA(output) || any(!nzchar(trimws(output)))) {
    .nm_stop("OUTPUT names must be non-missing, non-empty strings.")
  }
  unknown <- setdiff(output, catalog$name)
  if (length(unknown)) {
    .nm_stop("Unknown OUTPUT column(s): ", paste(unknown, collapse = ", "), ".")
  }
  collisions <- output[match(output, catalog$name) %in% which(catalog$collision)]
  if (length(collisions)) {
    .nm_stop(
      "OUTPUT column(s) also occur in INPUT: ", paste(collisions, collapse = ", "),
      ". Rename the generated variable to keep the run table unambiguous."
    )
  }
  output
}

#' Define a NONMEM-style pharmacometric model
#'
#' The established LibeRation `PRED` and `ERROR` strings are retained. They are
#' compiled to C++ expression IR; they are not evaluated inside the numerical
#' loop by R.
#'
#' @param INPUT Required dataset column names.
#' @param OUTPUT Optional generated run columns. Available names are discovered
#'   from `$PK/$PRED` assignments and combined with standard engine outputs; see
#'   [nm_model_outputs()].
#' @param ADVAN ADVAN number.
#' @param TRANS TRANS parameterization number.
#' @param SS Model-level steady-state default.
#' @param DOSECMP Default dosing compartment.
#' @param OBSCMP Default observation compartment.
#' @param PRED Parameter/model assignment code.
#' @param ERROR Residual-error or user-likelihood assignment code. With
#'   `LIK_CONFIG = nm_lik_config(error = "likelihood")`, assign either `LIK`
#'   (a positive row probability/density) or `LOGLIK` (its logarithm). The
#'   compiled block may use `DV`, `F`/`PRED`/`IPRED`, model assignments,
#'   parameters, numeric input columns, and the Markov helpers `PREV_DV`,
#'   `PREV_TIME`, `DT`, and `FIRST`.
#' @param DES ODE derivative code for ADVAN6/13.
#' @param ALG Algebraic residual equations for an experimental index-1 DAE.
#' @param THETAS,OMEGAS,SIGMAS Parameter tables.
#' @param COVARIATES Dataset covariates exposed to `PRED`.
#' @param USE_ODE Whether an ODE solver is explicitly requested.
#' @param ODE_CONTROL Named list with `rtol`, `atol`, `max_steps`, and optional
#'   `initial_step`. ADVAN6 uses adaptive Dormand-Prince 5(4); ADVAN13 uses an
#'   A-stable adaptive implicit trapezoidal method.
#' @param IOV Number of trailing inter-occasion ETAs.
#' @param LIK_CONFIG Reserved likelihood configuration.
#' @param HMM_CONFIG Optional finite-state hidden Markov configuration created
#'   by [nm_hmm_config()] or [nm_cthmm_config()]. Its initial,
#'   transition/rate, and emission names refer to assignments in `ERROR`; the
#'   complete scaled forward likelihood runs in C++ and remains differentiable
#'   by CppAD.
#' @param KALMAN_CONFIG Optional linear Gaussian state-space configuration
#'   created by [nm_kalman_config()]. Filtering is part of the compiled exact
#'   likelihood; retrospective smoothing is available through
#'   [nm_kalman_decode()].
#' @param DDE_CONFIG Optional method-of-steps declaration created by
#'   [nm_dde_config()]. Delayed states use `LAG(A(i), delay_name)` in `DES`.
#' @param DAE_CONFIG Optional semi-explicit index-1 declaration created by
#'   [nm_dae_config()]. Algebraic residuals are supplied through `ALG`.
#' @param COMPONENTS Optional list of immutable offline [nm_component()]
#'   declarations expanded into compiled `$PK/$PRED` or `$DES` IR according
#'   to each component's scope.
#' @param EXPERIMENTAL Explicit experimental-engine acknowledgement and policy.
#' @param RE_CONFIG Optional nested/crossed random-effect design created by
#'   [nm_re_config()]. Every ETA is assigned to a grouping block; connected
#'   components are evaluated as independent conditional-likelihood units.
#' @param OUTCOMES Optional first-class observation declaration created by
#'   [nm_outcome()] or [nm_outcomes()]. When `ERROR` is omitted, LibeRation
#'   generates an editable normalized likelihood block. The declaration is
#'   also used for stochastic outcome simulation and family diagnostics.
#' @param SOLVER `auto`, `advan`, `matrix`, `ode`, `dde`, or `dae`.
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
                     OUTPUT = NULL,
                     ADVAN = 2L,
                     TRANS = 2L,
                     SS = 0L,
                     DOSECMP = 1L,
                     OBSCMP = 1L,
                     PRED = "",
                     ERROR = "Y = F",
                     DES = "",
                     ALG = "",
                     THETAS,
                     OMEGAS = NULL,
                     SIGMAS = NULL,
                     COVARIATES = NULL,
                     USE_ODE = FALSE,
                     ODE_CONTROL = NULL,
                     IOV = 0L,
                     LIK_CONFIG = NULL,
                     HMM_CONFIG = NULL,
                     KALMAN_CONFIG = NULL,
                     DDE_CONFIG = NULL,
                     DAE_CONFIG = NULL,
                     RE_CONFIG = NULL,
                     COMPONENTS = NULL,
                     EXPERIMENTAL = NULL,
                     OUTCOMES = NULL,
                     SOLVER = c("auto", "advan", "matrix", "ode", "dde", "dae"),
                     ERROR_TYPE = c("auto", "none", "additive", "proportional", "combined", "exponential", "power", "likelihood"),
                     GRAPH = NULL,
                     LAYOUT = NULL,
                     LANGUAGE = c("R", "C++")) {
  error_was_missing <- missing(ERROR)
  advan <- as.integer(ADVAN)
  supported <- c(1L, 2L, 3L, 4L, 6L, 11L, 12L, 13L)
  if (length(advan) != 1L || is.na(advan) || !advan %in% supported) {
    .nm_stop("ADVAN must be one of: ", paste(supported, collapse = ", "), ".")
  }
  solver <- match.arg(SOLVER)
  language <- match.arg(LANGUAGE)
  pred_source <- paste(PRED, collapse = "\n")
  des_source <- paste(DES, collapse = "\n")
  components <- .nm_components(COMPONENTS)
  pred_components <- Filter(function(value) identical(value$scope %||% "pred", "pred"), components)
  des_components <- Filter(function(value) identical(value$scope %||% "pred", "des"), components)
  if (length(pred_components)) {
    PRED <- paste(c(pred_source,
                    vapply(pred_components, nm_component_code, character(1))),
                  collapse = "\n")
  } else PRED <- pred_source
  if (length(des_components)) {
    DES <- paste(c(vapply(des_components, nm_component_code, character(1)), des_source),
                 collapse = "\n")
  } else DES <- des_source
  dde_config <- .nm_dde_config(DDE_CONFIG)
  dae_config <- .nm_dae_config(DAE_CONFIG)
  if (!is.null(dde_config) && !is.null(dae_config)) {
    .nm_stop("A model cannot currently combine DDE_CONFIG and DAE_CONFIG.")
  }
  if ((!is.null(dde_config) || !is.null(dae_config)) && !advan %in% c(6L, 13L)) {
    .nm_stop("DDE/DAE models require ADVAN6 or ADVAN13.")
  }
  hmm_config <- .nm_hmm_config(HMM_CONFIG)
  generated_hmm_error <- attr(hmm_config, "generated_error", exact = TRUE)
  if (!is.null(generated_hmm_error)) {
    ERROR <- if (isTRUE(error_was_missing)) generated_hmm_error else
      paste(ERROR, generated_hmm_error, sep = "\n")
  }
  kalman_config <- .nm_kalman_config(KALMAN_CONFIG)
  generated_kalman_error <- attr(kalman_config, "generated_error", exact = TRUE)
  if (!is.null(generated_kalman_error)) {
    ERROR <- if (isTRUE(error_was_missing)) generated_kalman_error else
      paste(ERROR, generated_kalman_error, sep = "\n")
  }
  outcomes <- .nm_outcomes(OUTCOMES)
  sequence_count <- sum(!vapply(
    list(outcomes, hmm_config, kalman_config), is.null, logical(1)
  ))
  if (sequence_count > 1L) {
    .nm_stop("Use only one of OUTCOMES, HMM_CONFIG, or KALMAN_CONFIG in a model.")
  }
  general_ctmc <- !is.null(outcomes) && any(vapply(outcomes, function(outcome) {
    identical(outcome$family, "continuous_time_markov") && length(outcome$initial) > 2L
  }, logical(1)))
  if (general_ctmc) {
    if (length(outcomes) != 1L) {
      .nm_stop("A general multi-state continuous-time Markov outcome must currently be the model's sole OUTCOME.")
    }
    if (!isTRUE(error_was_missing)) {
      .nm_stop("Omit ERROR for a general continuous-time Markov OUTCOME; LibeRation generates its matrix-exponential likelihood outputs.")
    }
    compiled_ctmc <- .nm_outcome_ctmc_hmm(outcomes[[1L]])
    hmm_config <- compiled_ctmc$config
    ERROR <- compiled_ctmc$error
  } else if (!is.null(outcomes) && isTRUE(error_was_missing)) {
    ERROR <- .nm_outcome_error(outcomes)
  }
  error_type <- .nm_error_type(ERROR, match.arg(ERROR_TYPE))
  if (!is.null(outcomes)) error_type <- "likelihood"
  if (!is.null(hmm_config)) error_type <- "likelihood"
  if (!is.null(kalman_config)) error_type <- "likelihood"
  lik_config <- .nm_lik_config(LIK_CONFIG, error_type, as.integer(IOV))
  if (!is.null(hmm_config)) {
    if (!lik_config$error %in% c("auto", "likelihood")) {
      .nm_stop("HMM_CONFIG requires LIK_CONFIG$error = 'likelihood'.")
    }
    lik_config$error <- "likelihood"
  }
  if (!is.null(kalman_config)) {
    if (!lik_config$error %in% c("auto", "likelihood")) {
      .nm_stop("KALMAN_CONFIG requires LIK_CONFIG$error = 'likelihood'.")
    }
    if (lik_config$blq_method != "none") {
      .nm_stop("KALMAN_CONFIG does not currently support censored/BLQ observations.")
    }
    lik_config$error <- "likelihood"
  }
  error_type <- lik_config$error
  theta <- .nm_parameter_table(THETAS, "THETA", required = TRUE)
  omega <- .nm_omega_table(OMEGAS)
  sigma <- .nm_parameter_table(SIGMAS, "SIGMA")
  if (!identical(lik_config$ar1_source %||% "fixed", "fixed")) {
    available <- if (identical(lik_config$ar1_source, "theta")) nrow(theta) else nrow(sigma)
    if (lik_config$ar1_index < 1L || lik_config$ar1_index > available) {
      .nm_stop(
        "Estimated AR(1) parameter ", lik_config$ar1_parameter,
        " is outside the model's ", toupper(lik_config$ar1_source), " table."
      )
    }
  }
  if (length(lik_config$residual_groups)) {
    for (group in lik_config$residual_groups) {
      for (source in c("theta", "sigma")) {
        indices <- group$index[group$source == source]
        available <- if (source == "theta") nrow(theta) else nrow(sigma)
        if (length(indices) && any(indices < 1L | indices > available)) {
          .nm_stop("Residual group '", group$label, "' references an unavailable ",
                   toupper(source), " parameter.")
        }
      }
      .nm_residual_group_value(group, theta$Value, sigma$Value)
    }
  }
  if (lik_config$omega == "full" && .nm_n_eta(omega) > 1L &&
      nrow(omega) != .nm_n_eta(omega) * (.nm_n_eta(omega) + 1L) / 2L) {
    .nm_stop("A full OMEGA requires the complete lower triangle.")
  }
  if (!error_type %in% c("none", "likelihood") && nrow(sigma) == 0L) {
    .nm_stop("Residual error type '", error_type, "' requires SIGMAS.")
  }
  if (advan %in% c(6L, 13L) && !nzchar(trimws(DES))) {
    .nm_stop("ADVAN", advan, " requires a DES block.")
  }
  if (!nzchar(trimws(PRED))) .nm_stop("PRED code must not be empty.")
  n_eta <- .nm_n_eta(omega)
  re_config <- .nm_re_config(RE_CONFIG, n_eta)
  if (!is.null(re_config) && lik_config$iov > 0L) {
    .nm_stop("Use either RE_CONFIG or legacy IOV expansion, not both.")
  }
  if (!is.null(re_config) && (!is.null(hmm_config) || !is.null(kalman_config) ||
      lik_config$sigma_corr != "independent" || length(lik_config$residual_groups))) {
    .nm_stop(
      "RE_CONFIG currently requires independent row residuals and cannot be combined ",
      "with HMM/KALMAN_CONFIG; sequence state belongs to structural subjects rather than ",
      "connected random-effect components."
    )
  }
  if (!is.null(re_config)) {
    eta_block <- integer(n_eta)
    for (block_index in seq_along(re_config$blocks)) {
      eta_block[re_config$blocks[[block_index]]$etas] <- block_index
    }
    cross_block <- omega$ROW != omega$COL &
      eta_block[omega$ROW] != eta_block[omega$COL] & omega$Value != 0
    if (any(cross_block)) {
      .nm_stop("OMEGA correlation is allowed within, but not between, random-effect blocks.")
    }
  }
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
  error_info <- .nm_compile_error_ir(
    ERROR, pred_ir, nrow(theta), n_eta, INPUT, compiler_covariates, error_type,
    hmm_config, kalman_config
  )
  des_info <- if (advan %in% c(6L, 13L)) {
    .nm_compile_des_ir(
      DES, pred_ir, nrow(theta), n_eta, compiler_covariates,
      dde_config = dde_config,
      algebraic_variables = dae_config$variables %||% character()
    )
  } else NULL
  if (!is.null(dde_config)) {
    history <- dde_config$history
    if (length(history) == 1L) history <- rep(history, des_info$n_state)
    if (length(history) != des_info$n_state) {
      .nm_stop("DDE history must be scalar or contain one value per DES state.")
    }
    dde_config$history <- history
    dde_config$lags <- des_info$lags
    if (identical(solver, "auto")) solver <- "dde"
  }
  alg_ir <- .nm_compile_alg_ir(
    ALG, dae_config, pred_ir, des_info$n_state %||% 0L,
    nrow(theta), n_eta, compiler_covariates
  )
  if (!is.null(dae_config) && identical(solver, "auto")) solver <- "dae"
  experimental_features <- c(
    if (!is.null(dde_config)) "delay differential equations",
    if (!is.null(dae_config)) "index-1 differential-algebraic equations",
    if (inherits(hmm_config, "nm_factorial_hmm_config")) "factorial hidden Markov models",
    if (inherits(kalman_config, "nm_switching_state_space_config")) "switching state-space models",
    if (length(components)) "offline hybrid components"
  )
  experimental <- .nm_experimental_config(EXPERIMENTAL, experimental_features)
  ode_control <- .nm_ode_control(ODE_CONTROL, advan)
  graph <- GRAPH %||% if (is.null(des_info)) .nm_known_graph(advan) else .nm_ode_graph(des_info$n_state)
  n_state <- des_info$n_state %||% nrow(graph$compartments %||% data.frame())
  output_catalog <- .nm_output_catalog(pred_ir, n_state, n_eta, INPUT)
  selected_output <- .nm_validate_outputs(OUTPUT, output_catalog)
  outcome_outputs <- .nm_outcome_symbols(outcomes)
  missing_outcome_outputs <- setdiff(outcome_outputs, output_catalog$name)
  if (length(missing_outcome_outputs)) {
    .nm_stop(
      "OUTCOMES references model assignment(s) not found in $PK/$PRED: ",
      paste(missing_outcome_outputs, collapse = ", "), "."
    )
  }
  selected_output <- unique(c(selected_output, outcome_outputs))
  structure(
    list(
      version = 2L,
      INPUT = unique(as.character(INPUT)),
      OUTPUT = selected_output,
      ADVAN = advan,
      TRANS = as.integer(TRANS),
      SS = as.integer(SS),
      DOSECMP = as.integer(DOSECMP),
      OBSCMP = as.integer(OBSCMP),
      PRED = pred_source,
      ERROR = paste(ERROR, collapse = "\n"),
      DES = des_source,
      ALG = paste(ALG, collapse = "\n"),
      THETAS = theta,
      OMEGAS = omega,
      SIGMAS = sigma,
      COVARIATES = unique(as.character(COVARIATES %||% character())),
      USE_ODE = isTRUE(USE_ODE) || advan %in% c(6L, 13L),
      ODE_CONTROL = ode_control,
      IOV = as.integer(IOV),
      LIK_CONFIG = lik_config,
      HMM_CONFIG = hmm_config,
      KALMAN_CONFIG = kalman_config,
      DDE_CONFIG = dde_config,
      DAE_CONFIG = dae_config,
      RE_CONFIG = re_config,
      COMPONENTS = components,
      EXPERIMENTAL = experimental,
      OUTCOMES = outcomes,
      outcome_error_generated = !is.null(outcomes) && isTRUE(error_was_missing),
      SOLVER = solver,
      ERROR_TYPE = error_type,
      GRAPH = graph,
      LAYOUT = LAYOUT,
      LANGUAGE = language,
      pred_ir = pred_ir,
      error_ir = error_info$ir %||% NULL,
      likelihood_output = error_info$output %||% NULL,
      likelihood_scale = error_info$scale %||% NULL,
      des_ir = des_info$ir %||% NULL,
      alg_ir = alg_ir,
      n_eta = n_eta,
      n_state = n_state,
      output_catalog = output_catalog
    ),
    class = "nm_model"
  )
}

#' List selectable generated model-run columns
#'
#' This is a static catalogue produced while `$PK/$PRED` is compiled. It does
#' not execute the model or require a dataset. Standard engine outputs are
#' included alongside named assignment results.
#'
#' @param model An [nm_model()] or compiled [NMEngine()].
#' @return A data frame with output name, source, run availability, description,
#'   and selection status.
#' @export
nm_model_outputs <- function(model) {
  if (inherits(model, "NMEngine")) model <- model$model
  if (!inherits(model, "nm_model")) .nm_stop("`model` must be an nm_model or NMEngine.")
  catalog <- model$output_catalog %||% .nm_output_catalog(
    model$pred_ir, model$n_state, model$n_eta, model$INPUT
  )
  catalog$selected <- catalog$name %in% (model$OUTPUT %||% character())
  catalog
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
  if (!is.null(x$HMM_CONFIG)) {
    cat("  hidden states:", paste(x$HMM_CONFIG$states, collapse = ", "),
        " sequence:", if (isTRUE(x$HMM_CONFIG$by_dvid)) "subject + DVID" else "subject",
        "\n")
  }
  if (!is.null(x$KALMAN_CONFIG)) {
    cat("  Gaussian states:", paste(x$KALMAN_CONFIG$states, collapse = ", "),
        " sequence:", if (isTRUE(x$KALMAN_CONFIG$by_dvid)) "subject + DVID" else "subject",
        "\n")
  }
  if (isTRUE(x$EXPERIMENTAL$enabled)) {
    cat("  experimental:", paste(x$EXPERIMENTAL$features, collapse = "; "),
        if (isTRUE(x$EXPERIMENTAL$strict)) "[strict]" else "[fallback allowed]", "\n")
  }
  if (!is.null(x$OUTCOMES)) {
    cat("  outcomes:", paste(vapply(x$OUTCOMES, function(value) {
      paste0(value$name, " [", value$family, "]")
    }, character(1)), collapse = "; "), "\n")
  }
  invisible(x)
}

#' Report implemented feature support
#' @export
nm_support_matrix <- function() {
  feature <- c("ADVAN1", "ADVAN2", "ADVAN3", "ADVAN4", "ADVAN11", "ADVAN12",
               "matrix exponential", "steady-state bolus", "steady-state infusion",
               "nonlinear ODE steady state", "ADVAN6", "ADVAN13",
               "FO", "FOCE", "FOCEI", "LAPLACE", "ITS", "GQ", "IMP", "SAEM", "BAYES",
               "HMC", "NUTS", "NPML", "NPAG",
               "full OMEGA", "IOV", "M3/M4 BLQ", "finite mixtures", "priors",
               "user-defined likelihood", "first-class outcome families",
               "joint DVID endpoints", "Markov likelihood helpers",
               "finite-state hidden Markov models", "continuous-time HMM",
               "hidden semi-Markov models",
               "arbitrary-state continuous-time Markov", "Kalman filter",
               "extended Kalman filter", "unscented Kalman filter",
               "particle filter", "genealogical particle smoothing",
               "RTS state smoothing", "state-space simulation",
               "continuous-discrete SDE", "Euler-Maruyama", "Milstein",
               "AR1 residuals", "estimated AR1 correlation",
               "ARMA residual processes", "nested random effects",
               "crossed random effects",
               "cross-endpoint residual covariance",
               "stochastic simulation", "VPC", "multicategory VPC",
               "count VPC", "time-to-event VPC", "recurrent-event VPC",
               "competing-risk VPC", "bootstrap", "profile likelihood", "SCM",
               "NONMEM control-stream round-trip", "local queue", "remote queue",
               "React workbench", "model versions", "THETA bounds",
               "parallel estimation", "parallel simulation", "iteration gradients")
  experimental <- c(
    "delay differential equations", "index-1 DAE",
    "block-sparse algebraic solves", "stoichiometric QSP networks",
    "exact factorial HMM", "switching nonlinear state-space models",
    "filtered and smoothed regime probabilities",
    "offline dense neural components", "offline spline components",
    "offline Gaussian-process components", "learned DES dynamics"
  )
  values <- c(feature, experimental)
  reference_validated <- c(
    "ADVAN1", "ADVAN2", "ADVAN3", "ADVAN4", "ADVAN11", "ADVAN12",
    "matrix exponential", "FO"
  )
  data.frame(
    feature = values,
    status = c(rep("implemented", length(feature)), rep("experimental", length(experimental))),
    validation = ifelse(
      values %in% experimental, "experimental-smoke-tested",
      ifelse(values %in% reference_validated, "reference-validated", "unit/integration-tested")
    ),
    recommended_use = ifelse(values %in% experimental,
                             "experimental research only", "research and teaching"),
    stringsAsFactors = FALSE
  )
}
