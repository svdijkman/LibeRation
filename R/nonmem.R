.nm_model_arguments <- function(model, overrides = list()) {
  if (inherits(model, "NMEngine")) model <- model$model
  if (!inherits(model, "nm_model")) .nm_stop("`model` must be an nm_model or NMEngine.")
  arguments <- model[intersect(names(model), names(formals(nm_model)))]
  arguments[names(overrides)] <- overrides
  arguments
}

.nm_model_rebuild <- function(model, overrides = list()) {
  rebuilt <- do.call(nm_model, .nm_model_arguments(model, overrides))
  attr(rebuilt, "name") <- attr(model, "name", exact = TRUE)
  control <- attr(model, "nonmem_control", exact = TRUE)
  if (!is.null(control)) attr(rebuilt, "nonmem_control") <- control
  rebuilt
}

.nm_control_sections <- function(text) {
  lines <- strsplit(paste(text, collapse = "\n"), "\r?\n", perl = TRUE)[[1L]]
  sections <- list()
  current <- NULL
  for (line in lines) {
    if (grepl("^\\s*\\$[A-Za-z]+", line)) {
      match <- regexec("^\\s*\\$([A-Za-z]+)\\s*(.*)$", line, perl = TRUE)
      parts <- regmatches(line, match)[[1L]]
      current <- toupper(parts[[2L]])
      sections[[length(sections) + 1L]] <- list(
        name = current, header = trimws(parts[[3L]]), lines = character()
      )
    } else if (!is.null(current)) {
      sections[[length(sections)]]$lines <- c(sections[[length(sections)]]$lines, line)
    }
  }
  sections
}

.nm_control_body <- function(section, comments = FALSE) {
  lines <- c(section$header, section$lines)
  if (!comments) lines <- sub(";.*$", "", lines)
  trimws(lines[nzchar(trimws(lines))])
}

.nm_control_parameter_records <- function(sections, kind) {
  selected <- Filter(function(section) identical(section$name, kind), sections)
  if (!length(selected)) return(list())
  records <- list()
  for (section in selected) {
    text <- paste(.nm_control_body(section), collapse = " ")
    if (kind %in% c("OMEGA", "SIGMA")) {
      text <- sub("^\\s*(?:BLOCK|DIAGONAL)\\s*\\(\\s*[0-9]+\\s*\\)\\s*", "", text,
                  ignore.case = TRUE, perl = TRUE)
    }
    matches <- gregexpr(
      "\\([^)]*\\)\\s*(?:FIXED|FIX)?|[-+]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[Ee][-+]?[0-9]+)?\\s*(?:FIXED|FIX)?",
      text, perl = TRUE, ignore.case = TRUE
    )
    tokens <- regmatches(text, matches)[[1L]]
    tokens <- tokens[nzchar(tokens)]
    for (token in tokens) {
      fixed <- grepl("\\bFIX(?:ED)?\\b", token, ignore.case = TRUE)
      clean <- gsub("\\bFIX(?:ED)?\\b", "", token, ignore.case = TRUE)
      clean <- gsub("[()]", "", clean)
      values <- suppressWarnings(as.numeric(strsplit(trimws(clean), "\\s*,\\s*|\\s+", perl = TRUE)[[1L]]))
      values <- values[is.finite(values)]
      if (!length(values)) next
      initial <- if (length(values) == 1L) values[[1L]] else values[[2L]]
      records[[length(records) + 1L]] <- list(
        value = initial, fixed = fixed,
        lower = if (length(values) >= 2L) values[[1L]] else -Inf,
        upper = if (length(values) >= 3L) values[[3L]] else Inf
      )
    }
  }
  records
}

.nm_control_parameter_table <- function(records, kind) {
  index <- toupper(kind)
  if (!length(records)) {
    return(data.frame(stats::setNames(list(integer()), index), Value = numeric(),
                      FIX = logical(), LOWER = numeric(), UPPER = numeric()))
  }
  data.frame(
    stats::setNames(list(seq_along(records)), index),
    Value = vapply(records, `[[`, numeric(1), "value"),
    FIX = vapply(records, `[[`, logical(1), "fixed"),
    LOWER = vapply(records, `[[`, numeric(1), "lower"),
    UPPER = vapply(records, `[[`, numeric(1), "upper"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

.nm_control_omega <- function(sections) {
  selected <- Filter(function(section) identical(section$name, "OMEGA"), sections)
  rows <- list()
  eta_offset <- 0L
  parameter_index <- 0L
  warnings <- character()
  for (section in selected) {
    header <- paste(section$header, collapse = " ")
    block_match <- regexec("BLOCK\\s*\\(\\s*([0-9]+)\\s*\\)", header,
                           ignore.case = TRUE, perl = TRUE)
    block_parts <- regmatches(header, block_match)[[1L]]
    records <- .nm_control_parameter_records(list(section), "OMEGA")
    if (length(block_parts)) {
      size <- as.integer(block_parts[[2L]])
      needed <- size * (size + 1L) / 2L
      if (length(records) < needed) .nm_stop("$OMEGA BLOCK(", size, ") is incomplete.")
      if (length(records) > needed) warnings <- c(warnings, "$OMEGA BLOCK contains trailing values.")
      cursor <- 0L
      for (row in seq_len(size)) for (column in seq_len(row)) {
        cursor <- cursor + 1L
        parameter_index <- parameter_index + 1L
        record <- records[[cursor]]
        rows[[length(rows) + 1L]] <- data.frame(
          OMEGA = parameter_index, Value = record$value, FIX = record$fixed,
          ROW = eta_offset + row, COL = eta_offset + column,
          LOWER = record$lower, UPPER = record$upper
        )
      }
      eta_offset <- eta_offset + size
    } else {
      for (record in records) {
        eta_offset <- eta_offset + 1L
        parameter_index <- parameter_index + 1L
        rows[[length(rows) + 1L]] <- data.frame(
          OMEGA = parameter_index, Value = record$value, FIX = record$fixed,
          ROW = eta_offset, COL = eta_offset,
          LOWER = record$lower, UPPER = record$upper
        )
      }
    }
  }
  list(
    table = if (length(rows)) do.call(rbind, rows) else NULL,
    warnings = warnings
  )
}

.nm_control_input <- function(sections) {
  section <- Filter(function(section) identical(section$name, "INPUT"), sections)
  if (!length(section)) return(character())
  tokens <- unlist(strsplit(paste(.nm_control_body(section[[1L]]), collapse = " "),
                            "[[:space:],]+", perl = TRUE), use.names = FALSE)
  tokens <- tokens[nzchar(tokens)]
  keep <- vapply(tokens, function(token) {
    parts <- strsplit(token, "=", fixed = TRUE)[[1L]]
    !toupper(tail(parts, 1L)) %in% c("DROP", "SKIP")
  }, logical(1))
  tokens <- tokens[keep]
  toupper(vapply(strsplit(tokens, "=", fixed = TRUE), `[[`, character(1), 1L))
}

.nm_control_first <- function(sections, name) {
  selected <- Filter(function(section) identical(section$name, name), sections)
  if (!length(selected)) return(character())
  .nm_control_body(selected[[1L]])
}

.nm_control_table_columns <- function(sections) {
  tables <- Filter(function(section) identical(section$name, "TABLE"), sections)
  if (!length(tables)) return(character())
  tokens <- unlist(lapply(tables, function(section) {
    strsplit(
      paste(.nm_control_body(section), collapse = " "),
      "[[:space:],]+", perl = TRUE
    )[[1L]]
  }), use.names = FALSE)
  tokens <- gsub("^[('\\\"]+|[)'\\\"]+$", "", tokens, perl = TRUE)
  options <- c(
    "APPEND", "NOAPPEND", "HEADER", "ONEHEADER", "NOHEADER", "PRINT",
    "NOPRINT", "FIRSTONLY", "NOFORWARD", "NOSUB", "PARAFILE"
  )
  keep <- nzchar(tokens) &
    grepl("^[A-Za-z][A-Za-z0-9_.]*$", tokens) &
    !toupper(tokens) %in% options &
    !grepl("=", tokens, fixed = TRUE)
  unique(tokens[keep])
}

.nm_control_math_functions <- function(code) {
  for (name in c("EXP", "LOG", "SQRT", "SIN", "COS", "TAN", "TANH",
                 "ABS", "EXPM1", "LOG1P", "IFELSE", "MIN", "MAX")) {
    code <- gsub(
      paste0("\\b", name, "(?=\\s*\\()"), tolower(name), code,
      ignore.case = TRUE, perl = TRUE
    )
  }
  code
}

#' Read and translate a NONMEM control stream
#'
#' The importer preserves unsupported records and reports them rather than
#' silently discarding them. The supported semantic subset includes INPUT,
#' DATA, SUBROUTINES, MODEL, PK/PRED, DES, ERROR, THETA, diagonal or BLOCK
#' OMEGA, SIGMA, ESTIMATION, COVARIANCE and TABLE records.
#'
#' @param x A control-stream path or character text.
#' @param strict If `TRUE`, fail when the translated model cannot be compiled.
#' @return An `nm_control_stream` with the translated `model`, execution
#'   metadata, preserved raw text, and a compatibility report.
#' @export
nm_control_read <- function(x, strict = TRUE) {
  if (length(x) == 1L && file.exists(x)) {
    source_path <- normalizePath(x, winslash = "/", mustWork = TRUE)
    text <- readLines(source_path, warn = FALSE)
  } else {
    source_path <- NULL
    text <- strsplit(paste(x, collapse = "\n"), "\r?\n", perl = TRUE)[[1L]]
  }
  sections <- .nm_control_sections(text)
  if (!length(sections)) .nm_stop("No NONMEM $RECORD sections were found.")
  names_present <- vapply(sections, `[[`, character(1), "name")
  supported <- c("PROBLEM", "INPUT", "DATA", "SUBROUTINES", "MODEL", "PK", "PRED",
                 "DES", "ERROR", "THETA", "OMEGA", "SIGMA", "ESTIMATION",
                 "COVARIANCE", "TABLE")
  unknown <- unique(setdiff(names_present, supported))
  subroutines <- paste(.nm_control_first(sections, "SUBROUTINES"), collapse = " ")
  advan_match <- regexec("ADVAN\\s*([0-9]+)", subroutines, ignore.case = TRUE, perl = TRUE)
  trans_match <- regexec("TRANS\\s*([0-9]+)", subroutines, ignore.case = TRUE, perl = TRUE)
  advan_parts <- regmatches(subroutines, advan_match)[[1L]]
  trans_parts <- regmatches(subroutines, trans_match)[[1L]]
  advan <- if (length(advan_parts)) as.integer(advan_parts[[2L]]) else 2L
  trans <- if (length(trans_parts)) as.integer(trans_parts[[2L]]) else 2L
  if (!advan %in% c(1L, 2L, 3L, 4L, 6L, 11L, 12L, 13L)) {
    .nm_stop("The control stream requests unsupported ADVAN", advan, ".")
  }
  input <- .nm_control_input(sections)
  tables <- lapply(Filter(function(section) identical(section$name, "TABLE"), sections),
                   function(section) paste(.nm_control_body(section), collapse = " "))
  table_columns <- .nm_control_table_columns(sections)
  theta <- .nm_control_parameter_table(.nm_control_parameter_records(sections, "THETA"), "THETA")
  if (!nrow(theta)) .nm_stop("The control stream does not define $THETA values.")
  omega_result <- .nm_control_omega(sections)
  sigma <- .nm_control_parameter_table(.nm_control_parameter_records(sections, "SIGMA"), "SIGMA")
  pred <- c(.nm_control_first(sections, "PK"), .nm_control_first(sections, "PRED"))
  pred <- .nm_control_math_functions(paste(pred, collapse = "\n"))
  des <- .nm_control_math_functions(paste(.nm_control_first(sections, "DES"), collapse = "\n"))
  error <- .nm_control_math_functions(paste(.nm_control_first(sections, "ERROR"), collapse = "\n"))
  error <- gsub("\\bEPS\\s*\\(", "ERR(", error, ignore.case = TRUE, perl = TRUE)
  error <- gsub("\\.EQ\\.", "==", error, ignore.case = TRUE, perl = TRUE)
  error <- gsub("\\.NE\\.", "!=", error, ignore.case = TRUE, perl = TRUE)
  error <- gsub("\\.GE\\.", ">=", error, ignore.case = TRUE, perl = TRUE)
  error <- gsub("\\.GT\\.", ">", error, ignore.case = TRUE, perl = TRUE)
  error <- gsub("\\.LE\\.", "<=", error, ignore.case = TRUE, perl = TRUE)
  error <- gsub("\\.LT\\.", "<", error, ignore.case = TRUE, perl = TRUE)
  if (!nzchar(trimws(error))) error <- "Y = F"
  if (!nzchar(trimws(pred))) .nm_stop("The control stream requires a $PK or $PRED block.")
  covariates <- setdiff(input, c("ID", "TIME", "EVID", "AMT", "RATE", "II", "SS",
                                  "CMT", "DV", "MDV", "CENS", "LLOQ", "DVID"))
  obs_cmp <- if (advan %in% c(2L, 4L, 12L)) 2L else 1L
  warnings <- omega_result$warnings
  estimation <- paste(.nm_control_first(sections, "ESTIMATION"), collapse = " ")
  user_likelihood <- grepl("\\b(?:LIKE|LIKELIHOOD)\\b", estimation,
                           ignore.case = TRUE, perl = TRUE)
  if (grepl("\\bIF\\s*\\(", paste(pred, des, error), ignore.case = TRUE)) {
    warnings <- c(warnings, "Runtime IF statements require manual conversion to tape-safe ifelse().")
  }
  arguments <- list(
    INPUT = unique(c(input, "ID", "TIME", "EVID", "AMT", "RATE", "II", "SS", "CMT", "DV", "MDV")),
    ADVAN = advan, TRANS = trans, OBSCMP = obs_cmp, PRED = pred, ERROR = error,
    DES = des, THETAS = theta, OMEGAS = omega_result$table, SIGMAS = sigma,
    COVARIATES = covariates,
    LIK_CONFIG = if (user_likelihood) nm_lik_config(error = "likelihood") else NULL
  )
  model_error <- NULL
  model <- tryCatch({
    translated <- do.call(nm_model, arguments)
    requested <- table_columns[!toupper(table_columns) %in% toupper(input)]
    if (length(requested)) {
      catalog <- nm_model_outputs(translated)
      matched <- match(toupper(requested), toupper(catalog$name))
      unsupported_table <- requested[is.na(matched)]
      if (length(unsupported_table)) {
        warnings <- c(
          warnings,
          paste0(
            "The following $TABLE columns are preserved but cannot be generated by this model: ",
            paste(unique(unsupported_table), collapse = ", "), "."
          )
        )
      }
      selected <- unique(catalog$name[stats::na.omit(matched)])
      if (length(selected)) {
        arguments$OUTPUT <- selected
        translated <- do.call(nm_model, arguments)
      }
    }
    translated
  }, error = function(error) {
    model_error <<- conditionMessage(error)
    NULL
  })
  if (is.null(model) && isTRUE(strict)) {
    .nm_stop("The NONMEM stream was parsed but its model needs manual translation: ", model_error)
  }
  data_record <- paste(.nm_control_first(sections, "DATA"), collapse = " ")
  data_path <- if (nzchar(data_record)) strsplit(data_record, "\\s+", perl = TRUE)[[1L]][[1L]] else NULL
  problem <- paste(.nm_control_first(sections, "PROBLEM"), collapse = " ")
  covariance <- paste(.nm_control_first(sections, "COVARIANCE"), collapse = " ")
  extra <- Filter(function(section) section$name %in% unknown, sections)
  compatibility <- list(
    translated = is.null(model_error), supported_records = intersect(unique(names_present), supported),
    preserved_records = unknown, warnings = unique(c(warnings, if (!is.null(model_error)) model_error)),
    requires_manual_translation = !is.null(model_error) || length(unknown) > 0L
  )
  result <- structure(list(
    model = model, problem = problem, data_path = data_path,
    data_record = data_record,
    input_record = paste(.nm_control_first(sections, "INPUT"), collapse = " "),
    model_record = paste(.nm_control_first(sections, "MODEL"), collapse = "\n"),
    estimation = estimation, estimation_present = "ESTIMATION" %in% names_present,
    covariance = covariance, covariance_present = "COVARIANCE" %in% names_present,
    tables = tables,
    extra_sections = extra, compatibility = compatibility,
    raw = paste(text, collapse = "\n"), source = source_path
  ), class = "nm_control_stream")
  if (!is.null(model)) {
    attr(model, "name") <- if (nzchar(problem)) problem else "Imported NONMEM model"
    attr(model, "nonmem_control") <- result[names(result) != "model"]
    result$model <- model
  }
  result
}

.nm_control_format_value <- function(value) {
  format(as.numeric(value), digits = 15L, scientific = FALSE, trim = TRUE)
}

.nm_control_write_parameter <- function(table, kind) {
  if (is.null(table) || !nrow(table)) return(character())
  vapply(seq_len(nrow(table)), function(index) {
    lower <- if ("LOWER" %in% names(table)) table$LOWER[[index]] else -Inf
    upper <- if ("UPPER" %in% names(table)) table$UPPER[[index]] else Inf
    value <- .nm_control_format_value(table$Value[[index]])
    record <- if (is.finite(lower) || is.finite(upper)) {
      paste0("(", if (is.finite(lower)) .nm_control_format_value(lower) else "-INF",
             ",", value, ",", if (is.finite(upper)) .nm_control_format_value(upper) else "INF", ")")
    } else value
    if (isTRUE(table$FIX[[index]])) paste(record, "FIX") else record
  }, character(1))
}

#' Write a NONMEM control stream
#'
#' @param x An `nm_model`, `NMEngine`, or object returned by
#'   [nm_control_read()].
#' @param file Optional destination `.ctl` path. If omitted, text is returned.
#' @param data Dataset path and options written to `$DATA`. The imported record
#'   is reused by default for an `nm_control_stream`.
#' @param estimation Optional `$ESTIMATION` options.
#' @param covariance Optional `$COVARIANCE` options; use `FALSE` to omit it.
#' @return Control-stream text, invisibly when written to `file`.
#' @export
nm_control_write <- function(x, file = NULL, data = NULL,
                             estimation = NULL, covariance = NULL) {
  control <- if (inherits(x, "nm_control_stream")) x else NULL
  model <- if (!is.null(control)) control$model else if (inherits(x, "NMEngine")) x$model else x
  if (!inherits(model, "nm_model")) .nm_stop("`x` must contain an nm_model.")
  if (!is.null(model$HMM_CONFIG)) {
    .nm_stop(
      "Automatic NONMEM control-stream export of HMM_CONFIG is not available; ",
      "the scaled forward recursion is a LibeRation engine feature. Export the ",
      "component expressions manually only if the target NONMEM implementation ",
      "provides an equivalent recursive likelihood."
    )
  }
  if (!is.null(model$KALMAN_CONFIG)) {
    .nm_stop(
      "Automatic NONMEM control-stream export of KALMAN_CONFIG is not available; ",
      "the compiled filter/smoother is a LibeRation engine feature."
    )
  }
  metadata <- if (!is.null(control)) control else attr(model, "nonmem_control", exact = TRUE) %||% list()
  problem <- metadata$problem %||% attr(model, "name", exact = TRUE) %||% "LibeRation model"
  input_record <- metadata$input_record %||% paste(model$INPUT, collapse = " ")
  if (is.null(data)) data <- metadata$data_record %||% metadata$data_path %||% "data.csv"
  lines <- c(
    paste("$PROBLEM", problem),
    paste("$INPUT", input_record),
    paste("$DATA", data),
    paste0("$SUBROUTINES ADVAN", model$ADVAN, " TRANS", model$TRANS)
  )
  if (nzchar(trimws(metadata$model_record %||% ""))) {
    lines <- c(lines, "$MODEL", paste0("  ", strsplit(metadata$model_record, "\n", fixed = TRUE)[[1L]]))
  }
  lines <- c(lines, "$PK", paste0("  ", strsplit(model$PRED, "\n", fixed = TRUE)[[1L]]))
  if (nzchar(trimws(model$DES))) lines <- c(lines, "$DES", paste0("  ", strsplit(model$DES, "\n", fixed = TRUE)[[1L]]))
  nonmem_error <- gsub("\\bERR\\s*\\(", "EPS(", model$ERROR,
                       ignore.case = TRUE, perl = TRUE)
  if (identical(model$LIK_CONFIG$error, "likelihood")) {
    if (identical(model$likelihood_output, "LIK")) {
      nonmem_error <- paste(nonmem_error, "Y = LIK", sep = "\n")
    } else if (identical(model$likelihood_output, "LOGLIK")) {
      nonmem_error <- paste(nonmem_error, "Y = exp(LOGLIK)", sep = "\n")
    }
  }
  lines <- c(lines, "$ERROR", paste0("  ", strsplit(nonmem_error, "\n", fixed = TRUE)[[1L]]),
             "$THETA", paste0("  ", .nm_control_write_parameter(model$THETAS, "THETA")))
  omega <- model$OMEGAS
  if (nrow(omega)) {
    correlated <- any(omega$ROW != omega$COL)
    if (correlated) {
      n_eta <- max(omega$ROW, omega$COL)
      lines <- c(lines, paste0("$OMEGA BLOCK(", n_eta, ")"),
                 paste0("  ", .nm_control_write_parameter(omega, "OMEGA")))
    } else {
      lines <- c(lines, "$OMEGA", paste0("  ", .nm_control_write_parameter(omega, "OMEGA")))
    }
  }
  if (nrow(model$SIGMAS)) lines <- c(lines, "$SIGMA", paste0("  ", .nm_control_write_parameter(model$SIGMAS, "SIGMA")))
  estimation <- estimation %||% metadata$estimation
  if (identical(model$LIK_CONFIG$error, "likelihood")) {
    if (is.null(estimation) || !nzchar(trimws(estimation))) {
      estimation <- "METHOD=COND LAPLACE LIKELIHOOD"
    } else if (!grepl("\\b(?:LIKE|LIKELIHOOD)\\b", estimation,
                      ignore.case = TRUE, perl = TRUE)) {
      estimation <- paste(estimation, "LIKELIHOOD")
    }
  }
  if ((!is.null(estimation) && nzchar(trimws(estimation))) || isTRUE(metadata$estimation_present)) {
    lines <- c(lines, trimws(paste("$ESTIMATION", estimation %||% "")))
  }
  covariance <- covariance %||% metadata$covariance
  if (!identical(covariance, FALSE) && ((!is.null(covariance) && nzchar(trimws(covariance))) ||
                                        isTRUE(metadata$covariance_present))) {
    lines <- c(lines, trimws(paste("$COVARIANCE", covariance %||% "")))
  }
  tables <- metadata$tables %||% list()
  if (length(tables)) {
    for (table in tables) lines <- c(lines, paste("$TABLE", table))
  } else if (length(model$OUTPUT)) {
    identifiers <- intersect(c("ID", "TIME", "DV"), model$INPUT)
    lines <- c(lines, paste("$TABLE", paste(unique(c(identifiers, model$OUTPUT)), collapse = " ")))
  }
  for (section in metadata$extra_sections %||% list()) {
    lines <- c(lines, paste0("$", section$name, if (nzchar(section$header)) paste0(" ", section$header) else ""), section$lines)
  }
  text <- paste(lines, collapse = "\n")
  if (!is.null(file)) {
    writeLines(text, file, useBytes = TRUE)
    return(invisible(text))
  }
  text
}
