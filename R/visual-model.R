.nm_diagram_empty_compartments <- function() {
  data.frame(
    id = integer(), name = character(), kind = character(),
    volume_parameter = character(), scale_parameter = character(),
    dose = logical(), observe = logical(), x = numeric(), y = numeric(),
    stringsAsFactors = FALSE
  )
}

.nm_diagram_empty_flows <- function() {
  data.frame(
    id = character(), from = integer(), to = integer(), type = character(),
    parameter = character(), secondary_parameter = character(),
    expression = character(), label = character(),
    stringsAsFactors = FALSE
  )
}

.nm_diagram_empty_parameters <- function() {
  data.frame(
    name = character(), initial = numeric(), lower = numeric(), upper = numeric(),
    fixed = logical(), iiv = logical(), eta_variance = numeric(),
    stringsAsFactors = FALSE
  )
}

.nm_diagram_rows <- function(value, empty) {
  if (is.null(value) || !length(value)) return(empty)
  if (is.data.frame(value)) return(as.data.frame(value, stringsAsFactors = FALSE))
  if (!is.list(value)) .nm_stop("Diagram tables must be data frames or lists of rows.")
  rows <- lapply(value, function(row) {
    row <- as.list(row)
    for (name in names(empty)) {
      cell <- row[[name]]
      if (is.null(cell) || !length(cell)) cell <- NA
      if (length(cell) > 1L) cell <- cell[[1L]]
      row[[name]] <- cell
    }
    as.data.frame(row[names(empty)], stringsAsFactors = FALSE)
  })
  if (!length(rows)) return(empty)
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output
}

.nm_diagram_symbol <- function(value, what) {
  value <- trimws(as.character(value))
  if (length(value) != 1L || is.na(value) ||
      !grepl("^[A-Za-z][A-Za-z0-9_]*$", value)) {
    .nm_stop(what, " must be one valid model symbol.")
  }
  toupper(value)
}

.nm_diagram_compartments <- function(value) {
  output <- .nm_diagram_rows(value, .nm_diagram_empty_compartments())
  if (!nrow(output)) .nm_stop("A model diagram requires at least one compartment.")
  output$id <- as.integer(output$id)
  output$name <- trimws(as.character(output$name))
  output$kind <- tolower(trimws(as.character(output$kind)))
  output$kind[is.na(output$kind) | !nzchar(output$kind)] <- "amount"
  output$volume_parameter <- toupper(trimws(as.character(output$volume_parameter)))
  output$volume_parameter[is.na(output$volume_parameter)] <- ""
  output$scale_parameter <- toupper(trimws(as.character(output$scale_parameter)))
  output$scale_parameter[is.na(output$scale_parameter)] <- ""
  output$scale_parameter[!nzchar(output$scale_parameter)] <-
    output$volume_parameter[!nzchar(output$scale_parameter)]
  output$dose <- as.logical(output$dose)
  output$dose[is.na(output$dose)] <- FALSE
  output$observe <- as.logical(output$observe)
  output$observe[is.na(output$observe)] <- FALSE
  output$x <- suppressWarnings(as.numeric(output$x))
  output$y <- suppressWarnings(as.numeric(output$y))
  output$x[!is.finite(output$x)] <- seq(120, by = 180, length.out = nrow(output))[!is.finite(output$x)]
  output$y[!is.finite(output$y)] <- 180
  if (anyNA(output$id) || any(output$id < 1L) || anyDuplicated(output$id)) {
    .nm_stop("Diagram compartment ids must be unique positive integers.")
  }
  if (any(!nzchar(output$name)) || anyDuplicated(toupper(output$name))) {
    .nm_stop("Diagram compartment names must be unique and non-empty.")
  }
  if (any(!output$kind %in% c("amount", "response"))) {
    .nm_stop("Diagram compartment kind must be `amount` or `response`.")
  }
  symbols <- c(output$volume_parameter[nzchar(output$volume_parameter)],
               output$scale_parameter[nzchar(output$scale_parameter)])
  if (length(symbols)) invisible(lapply(symbols, .nm_diagram_symbol, "Compartment parameter"))
  if (!any(output$dose)) output$dose[[1L]] <- TRUE
  if (!any(output$observe)) output$observe[[min(2L, nrow(output))]] <- TRUE
  output
}

.nm_diagram_flows <- function(value, compartment_ids) {
  output <- .nm_diagram_rows(value, .nm_diagram_empty_flows())
  if (!nrow(output)) return(output)
  output$id <- as.character(output$id)
  missing_id <- is.na(output$id) | !nzchar(trimws(output$id))
  output$id[missing_id] <- paste0("flow-", which(missing_id))
  output$from <- suppressWarnings(as.integer(output$from))
  output$to <- suppressWarnings(as.integer(output$to))
  output$from[is.na(output$from)] <- 0L
  output$to[is.na(output$to)] <- 0L
  output$type <- tolower(trimws(as.character(output$type)))
  output$type[is.na(output$type) | !nzchar(output$type)] <- "rate"
  output$parameter <- toupper(trimws(as.character(output$parameter)))
  output$secondary_parameter <- toupper(trimws(as.character(output$secondary_parameter)))
  output$expression <- trimws(as.character(output$expression))
  output$label <- trimws(as.character(output$label))
  for (name in c("parameter", "secondary_parameter", "expression", "label")) {
    output[[name]][is.na(output[[name]])] <- ""
  }
  allowed <- c(
    "rate", "clearance", "bidirectional_clearance", "michaelis_menten",
    "zero_order", "custom"
  )
  if (any(!output$type %in% allowed)) {
    .nm_stop("Diagram flow type must be one of: ", paste(allowed, collapse = ", "), ".")
  }
  if (any(!output$from %in% c(0L, compartment_ids)) ||
      any(!output$to %in% c(0L, compartment_ids)) ||
      any(output$from == output$to)) {
    .nm_stop("Every diagram flow must connect two distinct compartments or one external source/sink.")
  }
  if (any(output$from == 0L & output$to == 0L)) {
    .nm_stop("A diagram flow cannot connect the external source to itself.")
  }
  needs_primary <- output$type != "custom"
  if (any(needs_primary & !nzchar(output$parameter))) {
    .nm_stop("Every non-custom diagram flow requires a primary parameter.")
  }
  needs_secondary <- output$type == "michaelis_menten"
  if (any(needs_secondary & !nzchar(output$secondary_parameter))) {
    .nm_stop("Michaelis-Menten flows require VMAX and KM parameter names.")
  }
  if (any(output$type == "custom" & !nzchar(output$expression))) {
    .nm_stop("Custom nonlinear flows require a flux expression.")
  }
  if (any(output$type == "bidirectional_clearance" &
          (output$from == 0L | output$to == 0L))) {
    .nm_stop("Bidirectional clearance must connect two compartments.")
  }
  primary <- output$parameter[nzchar(output$parameter)]
  secondary <- output$secondary_parameter[nzchar(output$secondary_parameter)]
  if (length(c(primary, secondary))) {
    invisible(lapply(c(primary, secondary), .nm_diagram_symbol, "Flow parameter"))
  }
  if (anyDuplicated(output$id)) .nm_stop("Diagram flow ids must be unique.")
  output
}

.nm_diagram_parameters <- function(value, required) {
  output <- .nm_diagram_rows(value, .nm_diagram_empty_parameters())
  if (nrow(output)) {
    output$name <- toupper(trimws(as.character(output$name)))
    invisible(lapply(output$name, .nm_diagram_symbol, "Parameter name"))
    if (anyDuplicated(output$name)) .nm_stop("Diagram parameter names must be unique.")
  }
  names_all <- unique(c(output$name, required))
  missing <- setdiff(names_all, output$name)
  if (length(missing)) {
    output <- rbind(output, data.frame(
      name = missing, initial = 1, lower = NA_real_, upper = NA_real_,
      fixed = FALSE, iiv = TRUE, eta_variance = 0.1,
      stringsAsFactors = FALSE
    ))
  }
  output <- output[match(names_all, output$name), , drop = FALSE]
  output$initial <- suppressWarnings(as.numeric(output$initial))
  output$lower <- suppressWarnings(as.numeric(output$lower))
  output$upper <- suppressWarnings(as.numeric(output$upper))
  output$fixed <- as.logical(output$fixed)
  output$iiv <- as.logical(output$iiv)
  output$eta_variance <- suppressWarnings(as.numeric(output$eta_variance))
  output$initial[!is.finite(output$initial) | output$initial <= 0] <- 1
  output$lower[!is.finite(output$lower)] <- NA_real_
  output$upper[!is.finite(output$upper)] <- NA_real_
  output$fixed[is.na(output$fixed)] <- FALSE
  output$iiv[is.na(output$iiv)] <- TRUE
  output$eta_variance[!is.finite(output$eta_variance) | output$eta_variance <= 0] <- 0.1
  if (any(is.finite(output$lower) & output$lower >= output$initial) ||
      any(is.finite(output$upper) & output$upper <= output$initial)) {
    .nm_stop("Diagram parameter bounds must contain their positive initial value.")
  }
  rownames(output) <- NULL
  output
}

.nm_diagram_required_parameters <- function(compartments, flows) {
  unique(c(
    compartments$volume_parameter[nzchar(compartments$volume_parameter)],
    compartments$scale_parameter[nzchar(compartments$scale_parameter)],
    flows$parameter[nzchar(flows$parameter)],
    flows$secondary_parameter[nzchar(flows$secondary_parameter)]
  ))
}

#' Define a visual PK/PD model diagram
#'
#' The diagram is a serializable semantic representation. It supports linear
#' rate and clearance flows, Michaelis-Menten elimination, zero-order source
#' or sink terms, and arbitrary nonlinear flux expressions. Layout is kept
#' alongside, but separate from, the generated model code.
#'
#' @param compartments Data frame or list of compartment rows.
#' @param flows Data frame or list of flow rows.
#' @param parameters Optional parameter settings. Missing referenced parameters
#'   are added with positive initials and log-normal IIV.
#' @param advan ODE solver, ADVAN6 or ADVAN13.
#' @param residual Residual-error scaffold: additive, proportional, or combined.
#' @param covariates Optional dataset covariates available to custom fluxes.
#' @param title Diagram/model title.
#' @return A serializable `nm_model_diagram`.
#' @export
nm_model_diagram <- function(compartments, flows = NULL, parameters = NULL,
                             advan = 6L,
                             residual = c("additive", "proportional", "combined"),
                             covariates = NULL,
                             title = "Diagram model") {
  compartments <- .nm_diagram_compartments(compartments)
  flows <- .nm_diagram_flows(flows, compartments$id)
  required <- .nm_diagram_required_parameters(compartments, flows)
  parameters <- .nm_diagram_parameters(parameters, required)
  if (!nrow(parameters)) {
    .nm_stop("A model diagram requires at least one structural parameter.")
  }
  advan <- as.integer(advan)
  if (length(advan) != 1L || is.na(advan) || !advan %in% c(6L, 13L)) {
    .nm_stop("A model diagram must use ADVAN6 or ADVAN13.")
  }
  covariates <- unique(toupper(trimws(as.character(covariates %||% character()))))
  covariates <- covariates[nzchar(covariates)]
  if (length(covariates)) invisible(lapply(covariates, .nm_diagram_symbol, "Covariate"))
  title <- trimws(as.character(title))
  if (length(title) != 1L || is.na(title) || !nzchar(title)) title <- "Diagram model"
  structure(list(
    schema = "liber.model-diagram/1", version = 1L,
    title = title, advan = advan, residual = match.arg(residual),
    covariates = covariates, compartments = compartments, flows = flows,
    parameters = parameters, generated = NULL
  ), class = "nm_model_diagram")
}

.nm_diagram_concentration <- function(compartment, compartments) {
  row <- match(compartment, compartments$id)
  volume <- compartments$volume_parameter[[row]]
  if (nzchar(volume)) paste0("(A(", compartment, ") / ", volume, ")") else paste0("A(", compartment, ")")
}

.nm_diagram_expand_concentrations <- function(expression, compartments) {
  for (id in compartments$id) {
    expression <- gsub(
      paste0("\\bC\\s*\\(\\s*", id, "\\s*\\)"),
      .nm_diagram_concentration(id, compartments), expression, perl = TRUE
    )
  }
  expression
}

.nm_diagram_code <- function(diagram) {
  compartments <- diagram$compartments
  parameters <- diagram$parameters
  eta_index <- cumsum(parameters$iiv)
  pred <- vapply(seq_len(nrow(parameters)), function(index) {
    eta <- if (parameters$iiv[[index]]) paste0(" * exp(ETA(", eta_index[[index]], "))") else ""
    paste0(parameters$name[[index]], " = THETA(", index, ")", eta)
  }, character(1))
  for (index in seq_len(nrow(compartments))) {
    scale <- compartments$scale_parameter[[index]]
    pred <- c(pred, paste0("S", compartments$id[[index]], " = ", if (nzchar(scale)) scale else "1"))
  }

  terms <- stats::setNames(vector("list", nrow(compartments)), as.character(compartments$id))
  add <- function(id, sign, expression) {
    if (!id) return(invisible(NULL))
    terms[[as.character(id)]] <<- c(
      terms[[as.character(id)]], paste0(if (sign > 0) "+ (" else "- (", expression, ")")
    )
    invisible(NULL)
  }
  flows <- diagram$flows
  if (nrow(flows)) for (index in seq_len(nrow(flows))) {
    flow <- flows[index, , drop = FALSE]
    from <- flow$from[[1L]]
    to <- flow$to[[1L]]
    type <- flow$type[[1L]]
    if (type == "rate") {
      flux <- if (from == 0L) flow$parameter[[1L]] else
        paste0(flow$parameter[[1L]], " * A(", from, ")")
      add(from, -1, flux); add(to, 1, flux)
    } else if (type == "clearance") {
      if (from == 0L) .nm_stop("Clearance flows require a source compartment.")
      flux <- paste0(flow$parameter[[1L]], " * ", .nm_diagram_concentration(from, compartments))
      add(from, -1, flux); add(to, 1, flux)
    } else if (type == "bidirectional_clearance") {
      forward <- paste0(flow$parameter[[1L]], " * ", .nm_diagram_concentration(from, compartments))
      reverse <- paste0(flow$parameter[[1L]], " * ", .nm_diagram_concentration(to, compartments))
      add(from, -1, forward); add(to, 1, forward)
      add(to, -1, reverse); add(from, 1, reverse)
    } else if (type == "michaelis_menten") {
      if (from == 0L) .nm_stop("Michaelis-Menten flows require a source compartment.")
      concentration <- .nm_diagram_concentration(from, compartments)
      flux <- paste0(
        flow$parameter[[1L]], " * ", concentration, " / (",
        flow$secondary_parameter[[1L]], " + ", concentration, ")"
      )
      add(from, -1, flux); add(to, 1, flux)
    } else if (type == "zero_order") {
      flux <- flow$parameter[[1L]]
      add(from, -1, flux); add(to, 1, flux)
    } else {
      flux <- .nm_diagram_expand_concentrations(flow$expression[[1L]], compartments)
      add(from, -1, flux); add(to, 1, flux)
    }
  }
  des <- vapply(seq_len(nrow(compartments)), function(index) {
    id <- compartments$id[[index]]
    rhs <- terms[[as.character(id)]]
    if (!length(rhs)) rhs <- "0"
    rhs <- sub("^[+] ", "", paste(rhs, collapse = " "))
    paste0("DADT(", id, ") = ", rhs)
  }, character(1))
  list(pred = paste(pred, collapse = "\n"), des = paste(des, collapse = "\n"))
}

#' Preview code generated by a model diagram
#'
#' @param diagram An [nm_model_diagram()] object.
#' @return Generated `$PK/$PRED`, `$DES`, THETA, OMEGA, and SIGMA scaffolds.
#' @export
nm_diagram_preview <- function(diagram) {
  if (!inherits(diagram, "nm_model_diagram")) .nm_stop("`diagram` must be an nm_model_diagram.")
  code <- .nm_diagram_code(diagram)
  parameters <- diagram$parameters
  theta <- data.frame(
    THETA = seq_len(nrow(parameters)), Value = parameters$initial,
    FIX = parameters$fixed, LOWER = parameters$lower, UPPER = parameters$upper,
    stringsAsFactors = FALSE
  )
  effects <- which(parameters$iiv)
  omega <- if (length(effects)) data.frame(
    OMEGA = seq_along(effects), Value = parameters$eta_variance[effects],
    FIX = FALSE, ROW = seq_along(effects), COL = seq_along(effects),
    stringsAsFactors = FALSE
  ) else NULL
  residual <- switch(
    diagram$residual,
    additive = list(error = "Y = F + ERR(1)", sigma = c(0.1)),
    proportional = list(error = "Y = F * (1 + ERR(1))", sigma = c(0.05)),
    combined = list(error = "Y = F * (1 + ERR(1)) + ERR(2)", sigma = c(0.05, 0.1))
  )
  sigma <- data.frame(
    SIGMA = seq_along(residual$sigma), Value = residual$sigma, FIX = FALSE,
    stringsAsFactors = FALSE
  )
  list(PRED = code$pred, DES = code$des, ERROR = residual$error,
       THETAS = theta, OMEGAS = omega, SIGMAS = sigma)
}

#' Generate a LibeRation model from a visual diagram
#'
#' Generated `$DES` remains ordinary editable model code. The semantic diagram
#' and last generated code are retained in `GRAPH`, so later visual regeneration
#' can be previewed explicitly without silently overwriting manual edits.
#'
#' @param diagram An [nm_model_diagram()] object.
#' @param input Optional NONMEM-style input-column names.
#' @param output Optional generated run columns.
#' @return A serializable [nm_model()].
#' @export
nm_diagram_generate <- function(diagram,
                                input = c("ID", "TIME", "EVID", "AMT", "RATE", "II",
                                          "SS", "CMT", "DV", "MDV"),
                                output = c("PRED", "IPRED", "CWRES")) {
  preview <- nm_diagram_preview(diagram)
  compartments <- diagram$compartments
  graph <- diagram
  graph$generated <- list(
    PRED = preview$PRED, DES = preview$DES,
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")
  )
  class(graph) <- c("nm_visual_model", class(diagram))
  model <- nm_model(
    INPUT = unique(c(input, diagram$covariates)), OUTPUT = output,
    ADVAN = diagram$advan, TRANS = 1L,
    DOSECMP = compartments$id[which(compartments$dose)[[1L]]],
    OBSCMP = compartments$id[which(compartments$observe)[[1L]]],
    PRED = preview$PRED, DES = preview$DES, ERROR = preview$ERROR,
    THETAS = preview$THETAS, OMEGAS = preview$OMEGAS, SIGMAS = preview$SIGMAS,
    COVARIATES = diagram$covariates, USE_ODE = TRUE,
    SOLVER = "ode", GRAPH = graph,
    LAYOUT = compartments[, c("id", "x", "y"), drop = FALSE]
  )
  attr(model, "name") <- diagram$title
  model
}

#' Recover a visual diagram from a generated model
#'
#' @param model An [nm_model()] or [NMEngine].
#' @return The stored `nm_model_diagram`, or `NULL` for a code-only model.
#' @export
nm_model_diagram_get <- function(model) {
  if (inherits(model, "NMEngine")) model <- model$model
  if (!inherits(model, "nm_model")) .nm_stop("`model` must be an nm_model or NMEngine.")
  if (!inherits(model$GRAPH, "nm_model_diagram")) return(NULL)
  diagram <- model$GRAPH
  class(diagram) <- "nm_model_diagram"
  diagram
}
