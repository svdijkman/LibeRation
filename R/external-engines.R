.nm_execution_engine <- function(engine, include_liber = TRUE) {
  choices <- c(if (isTRUE(include_liber)) "liber", "nonmem", "nlmixr2")
  match.arg(tolower(as.character(engine)), choices)
}

.nm_nonmem_command <- function() {
  configured <- Sys.getenv("LIBERATION_NONMEM_COMMAND", unset = "")
  if (!nzchar(configured)) configured <- getOption("LibeRation.nonmem_command", "")
  if (length(configured) != 1L || is.na(configured)) configured <- ""
  candidate <- if (nzchar(configured)) configured else "execute"
  resolved <- unname(Sys.which(candidate))
  if (!nzchar(resolved) && file.exists(candidate)) {
    resolved <- normalizePath(candidate, winslash = "/", mustWork = TRUE)
  }
  resolved
}

.nm_nonmem_process_arguments <- function(command, arguments, directory) {
  result <- list(
    command = command, args = arguments, wd = directory, echo = FALSE,
    error_on_status = FALSE, stderr = "2>&1"
  )
  if (.Platform$OS.type == "windows") {
    # Rtools prepends its own Perl to PATH inside R sessions. PsN's Windows
    # launcher must instead find the Perl installation beside `execute.bat`.
    paths <- strsplit(Sys.getenv("PATH"), .Platform$path.sep, fixed = TRUE)[[1L]]
    installation_root <- normalizePath(
      dirname(dirname(dirname(command))), winslash = "/", mustWork = TRUE
    )
    canonical_paths <- vapply(paths, function(path) {
      if (!nzchar(path) || !file.exists(path)) return(path)
      normalizePath(path, winslash = "/", mustWork = TRUE)
    }, character(1))
    bundled <- startsWith(
      tolower(canonical_paths), paste0(tolower(installation_root), "/")
    )
    paths <- unique(c(dirname(command), paths[bundled], paths[!bundled]))
    quote <- function(value) {
      value <- as.character(value)
      if (grepl("\"", value, fixed = TRUE) ||
          grepl("\r", value, fixed = TRUE) ||
          grepl("\n", value, fixed = TRUE)) {
        .nm_stop("The configured NONMEM command contains unsafe command-shell text.")
      }
      if (grepl("[[:space:],]", value)) paste0("\"", value, "\"") else value
    }
    command_line <- paste(
      vapply(c(basename(command), arguments), quote, character(1)),
      collapse = " "
    )
    result$command <- Sys.getenv("COMSPEC", unset = "cmd.exe")
    result$args <- c("/d", "/s", "/c", command_line)
    result$env <- c(PATH = paste(paths, collapse = .Platform$path.sep))
  }
  result
}

.nm_nonmem_log <- function(process) {
  lines <- strsplit(
    paste(process$stdout %||% "", process$stderr %||% "", sep = "\n"),
    "\n", fixed = TRUE
  )[[1L]]
  private <- grepl(
    "License Registered to:|Expiration Date:|Days until program expires|^Current Date:",
    trimws(lines), ignore.case = TRUE
  )
  lines <- lines[!private & nzchar(trimws(lines))]
  if (length(lines)) cat(paste0(lines, collapse = "\n"), "\n")
  invisible(lines)
}

#' Report available model-execution engines
#'
#' Availability is evaluated in the current R process. Remote LibeRties
#' workers perform the same check inside their isolated execution environment.
#'
#' @return A data frame describing native LibeR, NONMEM/PsN, and nlmixr2.
#' @export
nm_execution_engines <- function() {
  nonmem <- .nm_nonmem_command()
  nlmixr <- requireNamespace("nlmixr2", quietly = TRUE) &&
    requireNamespace("nlmixr2est", quietly = TRUE) &&
    requireNamespace("rxode2", quietly = TRUE)
  data.frame(
    id = c("liber", "nonmem", "nlmixr2"),
    label = c("LibeR", "NONMEM", "nlmixr2"),
    available = c(TRUE, nzchar(nonmem), nlmixr),
    version = c(
      as.character(utils::packageVersion("LibeRation")),
      if (nzchar(nonmem)) "PsN execute" else "",
      if (nlmixr) as.character(utils::packageVersion("nlmixr2")) else ""
    ),
    location = c(
      getNamespaceInfo(asNamespace("LibeRation"), "path"), nonmem,
      if (nlmixr) getNamespaceInfo(asNamespace("nlmixr2"), "path") else ""
    ),
    stringsAsFactors = FALSE
  )
}

.nm_external_source_lines <- function(source) {
  source <- enc2utf8(as.character(source %||% ""))
  source <- gsub(";[[:space:]]*", "\n", source, perl = TRUE)
  lines <- trimws(strsplit(source, "\n", fixed = TRUE)[[1L]])
  lines[nzchar(lines)]
}

.nm_nonmem_prepared_model <- function(model, table_columns) {
  prepared <- model
  for (name in intersect(
    c("PRED", "ERROR", "DES", "PK_SOURCE", "PRED_SOURCE"), names(prepared)
  )) {
    prepared[[name]] <- paste(.nm_external_source_lines(prepared[[name]]), collapse = "\n")
  }
  if (!grepl("(^|[^A-Za-z0-9_])IPRED[[:space:]]*=", prepared$ERROR,
             ignore.case = TRUE, perl = TRUE)) {
    prepared$ERROR <- paste(prepared$ERROR, "IPRED = F", sep = "\n")
  }
  metadata <- attr(prepared, "nonmem_control", exact = TRUE) %||% list()
  metadata$estimation <- NULL
  metadata$estimation_present <- FALSE
  metadata$covariance <- NULL
  metadata$covariance_present <- FALSE
  metadata$tables <- list(paste(
    paste(unique(table_columns), collapse = " "),
    "FILE=model.tab ONEHEADER NOAPPEND NOPRINT"
  ))
  attr(prepared, "nonmem_control") <- metadata
  prepared
}

.nm_nonmem_data_write <- function(data, path) {
  data <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  connection <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)
  writeLines(paste0("@", paste(names(data), collapse = ",")), connection,
             useBytes = TRUE)
  utils::write.table(
    data, connection, sep = ",", row.names = FALSE, col.names = FALSE,
    quote = FALSE, na = "-99", qmethod = "double"
  )
  invisible(path)
}

.nm_nonmem_estimation_record <- function(method, arguments) {
  method <- toupper(as.character(method %||% "FOCEI"))
  maxit <- max(0L, as.integer(arguments$maxit %||% 200L))
  print <- max(0L, as.integer(arguments$print_every %||% 0L))
  common <- c(paste0("MAXEVAL=", maxit), paste0("PRINT=", print))
  record <- switch(
    method,
    FO = "METHOD=0",
    FOCE = "METHOD=COND NOINTERACTION",
    FOCEI = "METHOD=COND INTERACTION",
    LAPLACE = "METHOD=COND LAPLACE INTERACTION",
    ITS = paste("METHOD=ITS INTERACTION", paste0("NITER=", maxit)),
    IMP = paste(
      "METHOD=IMP INTERACTION",
      paste0("NITER=", maxit),
      paste0("ISAMPLE=", max(5L, as.integer(arguments$n_imp %||% 200L)))
    ),
    SAEM = paste(
      "METHOD=SAEM INTERACTION",
      paste0("NBURN=", max(0L, as.integer(arguments$burn %||% 60L))),
      paste0("NITER=", max(2L, as.integer(arguments$n_iter %||% maxit)))
    ),
    .nm_stop(
      "NONMEM execution currently supports FO, FOCE, FOCEI, LAPLACE, ITS, IMP, and SAEM; received ",
      method, "."
    )
  )
  paste(c(record, common), collapse = " ")
}

.nm_nonmem_covariance_record <- function(arguments) {
  if (!isTRUE(arguments$covariance)) return(FALSE)
  type <- tolower(as.character(arguments$covariance_type %||% "sandwich"))
  matrix <- switch(type, sandwich = "RS", hessian = "R", opg = "S", "RS")
  paste0("MATRIX=", matrix, " PRINT=E")
}

.nm_nonmem_read_table <- function(path) {
  if (!file.exists(path) || !file.info(path)$size) return(data.frame())
  lines <- readLines(path, warn = FALSE, encoding = "unknown")
  lines <- lines[nzchar(trimws(lines)) & !grepl("^[[:space:]]*TABLE NO[.]?", lines)]
  if (length(lines) < 2L) return(data.frame())
  header <- strsplit(trimws(lines[[1L]]), "[[:space:]]+")[[1L]]
  body_lines <- lines[-1L]
  body_lines <- body_lines[trimws(body_lines) != trimws(lines[[1L]])]
  if (!length(body_lines)) return(data.frame())
  body <- paste(body_lines, collapse = "\n")
  result <- tryCatch(
    utils::read.table(
      text = body, header = FALSE, col.names = make.unique(header),
      check.names = FALSE, stringsAsFactors = FALSE, fill = TRUE,
      na.strings = c("-99", "-9.9E+01", "***********")
    ),
    error = function(error) data.frame()
  )
  rownames(result) <- NULL
  result
}

.nm_nonmem_parameter <- function(final, pattern, fallback) {
  hit <- grep(pattern, names(final), value = TRUE, ignore.case = TRUE, perl = TRUE)
  if (!length(hit)) return(as.numeric(fallback))
  value <- suppressWarnings(as.numeric(final[[hit[[1L]]]]))
  if (!length(value) || !is.finite(value)) as.numeric(fallback) else value
}

.nm_nonmem_matrix_file <- function(path) {
  table <- .nm_nonmem_read_table(path)
  if (!nrow(table) || !ncol(table)) return(NULL)
  if (!all(vapply(table, is.numeric, logical(1)))) {
    labels <- as.character(table[[1L]])
    table <- table[-1L]
  } else labels <- names(table)
  value <- suppressWarnings(as.matrix(data.frame(lapply(table, as.numeric))))
  if (!nrow(value) || nrow(value) != ncol(value)) return(NULL)
  names <- labels[seq_len(nrow(value))]
  dimnames(value) <- list(names, names)
  value
}

.nm_nonmem_listing_info <- function(path) {
  lines <- if (file.exists(path)) readLines(path, warn = FALSE, encoding = "unknown") else character()
  successful <- any(grepl("MINIMIZATION SUCCESSFUL", lines, fixed = TRUE))
  evaluations <- grep("NO[.] OF FUNCTION EVALUATIONS USED:", lines, value = TRUE)
  evaluations <- if (length(evaluations)) {
    suppressWarnings(as.integer(sub(".*USED:[[:space:]]*([0-9]+).*", "\\1", tail(evaluations, 1L))))
  } else NA_integer_
  estimate_time <- grep("Elapsed estimation time in seconds:", lines, value = TRUE)
  covariance_time <- grep("Elapsed covariance time in seconds:", lines, value = TRUE)
  number <- function(value) {
    if (!length(value)) return(NA_real_)
    suppressWarnings(as.numeric(sub(".*seconds:[[:space:]]*", "", tail(value, 1L))))
  }
  list(
    successful = successful, evaluations = evaluations,
    estimate_seconds = number(estimate_time),
    covariance_seconds = number(covariance_time),
    message = if (successful) "NONMEM minimization successful" else
      "NONMEM returned estimates without a successful-minimization marker"
  )
}

.nm_nonmem_artifact_bundle <- function(directory, operation, model, data) {
  candidates <- c(
    "model.lst", "model.ctl", "model.ext", "model.cov", "model.cor",
    "model.phi", "model.tab"
  )
  files <- lapply(candidates[file.exists(file.path(directory, candidates))], function(name) {
    lines <- readLines(file.path(directory, name), warn = FALSE, encoding = "unknown")
    .nm_audit_file(
      name, c(
        if (name == "model.ctl") "; Executed NONMEM control stream generated by LibeRation." else character(),
        lines
      ),
      if (name == "model.tab") "text/tab-separated-values" else "text/plain"
    )
  })
  bundle <- list(
    schema = "liberation.audit-artifacts", version = 1L,
    engine = "nonmem", source = "original", operation = operation,
    created = .nm_workspace_now(),
    model_sha256 = digest::digest(nm_model_to_contract(model), algo = "sha256"),
    data_sha256 = digest::digest(as.data.frame(data), algo = "sha256"),
    files = files
  )
  .nm_audit_bundle_validate(bundle)
  bundle
}

.nm_nonmem_fit <- function(directory, model, data, method, elapsed) {
  ext <- .nm_nonmem_read_table(file.path(directory, "model.ext"))
  if (!nrow(ext)) .nm_stop("NONMEM did not produce a readable .ext parameter table.")
  final_rows <- which(suppressWarnings(as.numeric(ext$ITERATION)) == -1000000000)
  final <- ext[if (length(final_rows)) tail(final_rows, 1L) else nrow(ext), , drop = FALSE]
  theta <- vapply(seq_len(nrow(model$THETAS)), function(index) {
    .nm_nonmem_parameter(final, paste0("^THETA\\(?", index, "\\)?$"), model$THETAS$Value[[index]])
  }, numeric(1))
  sigma <- vapply(seq_len(nrow(model$SIGMAS)), function(index) {
    row <- if ("ROW" %in% names(model$SIGMAS)) model$SIGMAS$ROW[[index]] else index
    column <- if ("COL" %in% names(model$SIGMAS)) model$SIGMAS$COL[[index]] else index
    .nm_nonmem_parameter(
      final, paste0("^SIGMA\\(", row, "[,]", column, "\\)$"),
      model$SIGMAS$Value[[index]]
    )
  }, numeric(1))
  omega <- vapply(seq_len(nrow(model$OMEGAS)), function(index) {
    row <- model$OMEGAS$ROW[[index]]
    column <- model$OMEGAS$COL[[index]]
    .nm_nonmem_parameter(
      final, paste0("^OMEGA\\(", row, "[,]", column, "\\)$"),
      model$OMEGAS$Value[[index]]
    )
  }, numeric(1))
  phi <- .nm_nonmem_read_table(file.path(directory, "model.phi"))
  eta_columns <- grep("^ETA[.(]", names(phi), value = TRUE, ignore.case = TRUE)
  eta <- if (length(eta_columns)) as.matrix(phi[, eta_columns, drop = FALSE]) else
    matrix(numeric(), length(unique(data$ID)), 0L)
  storage.mode(eta) <- "double"
  if (ncol(eta)) colnames(eta) <- paste0("ETA", seq_len(ncol(eta)))
  if (nrow(eta) != length(unique(data$ID))) {
    .nm_stop("NONMEM .phi subject count does not match the submitted dataset.")
  }
  info <- .nm_nonmem_listing_info(file.path(directory, "model.lst"))
  objective <- suppressWarnings(as.numeric(final$OBJ[[1L]] %||% NA_real_))
  output <- .nm_nonmem_read_table(file.path(directory, "model.tab"))
  covariance <- .nm_nonmem_matrix_file(file.path(directory, "model.cov"))
  correlation <- .nm_nonmem_matrix_file(file.path(directory, "model.cor"))
  covariance_result <- if (is.null(covariance)) NULL else list(
    status = "completed", type = "NONMEM", covariance = covariance,
    correlation = correlation, bread_source = "NONMEM output",
    bread_exact = NA
  )
  iterations <- suppressWarnings(max(
    as.integer(ext$ITERATION[ext$ITERATION >= 0]), na.rm = TRUE
  ))
  if (!is.finite(iterations)) iterations <- NA_integer_
  structure(list(
    version = 1L, method = toupper(method), execution_engine = "nonmem",
    objective = objective, theta = theta, omega = omega, sigma = sigma,
    eta = eta, convergence = if (info$successful) 0L else 1L,
    message = info$message, iterations = as.integer(iterations),
    objective_evaluations = info$evaluations,
    evaluations = stats::setNames(
      c(info$evaluations, NA_integer_), c("function", "gradient")
    ),
    model = model, data = data,
    diagnostics = list(
      population_gradient = "NONMEM engine; inspect the original .grd/.lst output when retained",
      optimizer = list(backend = "NONMEM via PsN", objective_backend = "NONMEM")
    ),
    covariance = covariance_result,
    timing = list(
      model_fit_seconds = info$estimate_seconds,
      covariance_seconds = info$covariance_seconds,
      total_seconds = as.numeric(elapsed)
    ),
    output = if (nrow(output)) output else NULL
  ), class = "nm_fit")
}

.nm_nonmem_run <- function(model, data, operation, arguments) {
  command <- .nm_nonmem_command()
  if (!nzchar(command)) {
    .nm_stop(
      "NONMEM execution requires PsN's `execute` command in PATH or ",
      "LIBERATION_NONMEM_COMMAND configured by the machine administrator."
    )
  }
  data <- nm_dataset(data)
  audit <- isTRUE(arguments$audit_artifacts)
  seed <- as.integer(arguments$seed %||% 20260713L)
  nsim <- max(1L, as.integer(arguments$nsim %||% 1L))
  table_columns <- unique(c(
    model$INPUT, "PRED", "IPRED",
    if (identical(operation, "estimate")) "CWRES", model$OUTPUT
  ))
  prepared <- .nm_nonmem_prepared_model(model, table_columns)
  directory <- tempfile("liberation-nonmem-")
  if (!dir.create(directory, recursive = TRUE, showWarnings = FALSE)) {
    .nm_stop("Unable to create the isolated NONMEM working directory.")
  }
  on.exit(unlink(directory, recursive = TRUE, force = TRUE), add = TRUE)
  missing_input <- setdiff(model$INPUT, names(data))
  if (length(missing_input)) {
    .nm_stop(
      "NONMEM input data are missing model column(s): ",
      paste(missing_input, collapse = ", "), "."
    )
  }
  .nm_nonmem_data_write(
    data[, model$INPUT, drop = FALSE], file.path(directory, "data.csv")
  )
  if (identical(operation, "estimate")) {
    stages <- arguments$stages %||% list()
    final_arguments <- if (length(stages)) {
      utils::modifyList(
        tail(stages, 1L)[[1L]]$arguments %||% list(),
        arguments[setdiff(names(arguments), "stages")]
      )
    } else arguments
    records <- if (length(stages)) {
      vapply(stages, function(stage) {
        .nm_nonmem_estimation_record(stage$method, stage$arguments %||% list())
      }, character(1))
    } else .nm_nonmem_estimation_record(arguments$method, arguments)
    estimation <- paste(records, collapse = "\n$ESTIMATION ")
    covariance <- .nm_nonmem_covariance_record(final_arguments)
  } else {
    estimation <- FALSE
    covariance <- FALSE
  }
  control <- nm_control_write(
    prepared, data = "data.csv IGNORE=@", estimation = estimation,
    covariance = covariance
  )
  if (identical(operation, "simulate")) {
    lines <- strsplit(control, "\n", fixed = TRUE)[[1L]]
    table_at <- grep("^[[:space:]]*\\$TABLE\\b", lines, ignore.case = TRUE)
    simulation_record <- paste0(
      "$SIMULATION (", seed, ") ONLYSIM NSUBPROBLEMS=", nsim
    )
    if (length(table_at)) {
      position <- table_at[[1L]]
      lines <- append(lines, simulation_record, after = position - 1L)
    } else {
      lines <- c(lines, simulation_record)
    }
    control <- paste(lines, collapse = "\n")
  }
  writeLines(control, file.path(directory, "model.ctl"), useBytes = TRUE)
  file.copy(file.path(directory, "model.ctl"), file.path(directory, "model.mod"))
  started <- proc.time()[["elapsed"]]
  process_arguments <- .nm_nonmem_process_arguments(
    command,
    c(
      "model.mod", "-directory=psn-run", "-clean=0",
      "-nm_output=ext,cov,cor,phi", "-no-prepend_options_to_lst",
      "-display_iterations"
    ),
    directory
  )
  timeout <- as.numeric(arguments$timeout_seconds %||% Inf)
  if (length(timeout) != 1L || is.na(timeout) || timeout <= 0) {
    .nm_stop("`timeout_seconds` must be a positive number when supplied.")
  }
  if (is.finite(timeout)) process_arguments$timeout <- timeout * 1000
  process <- do.call(processx::run, process_arguments)
  .nm_nonmem_log(process)
  elapsed <- unname(proc.time()[["elapsed"]] - started)
  listing <- file.path(directory, "model.lst")
  if (!file.exists(listing)) {
    errors <- list.files(
      file.path(directory, "psn-run"), pattern = "(error|FREPORT)",
      recursive = TRUE, full.names = TRUE, ignore.case = TRUE
    )
    detail <- paste(unlist(lapply(errors, function(path) {
      paste(readLines(path, warn = FALSE), collapse = "\n")
    })), collapse = "\n")
    if (!nzchar(trimws(detail))) detail <- paste(process$stderr, process$stdout)
    .nm_stop("NONMEM/PsN did not produce a listing file. ", trimws(detail))
  }
  result <- if (identical(operation, "estimate")) {
    final_method <- final_arguments$method %||% if (length(arguments$stages %||% list())) {
      tail(arguments$stages, 1L)[[1L]]$method
    } else "FOCEI"
    .nm_nonmem_fit(
      directory, model, data, final_method, elapsed
    )
  } else {
    output <- .nm_nonmem_read_table(file.path(directory, "model.tab"))
    if (!nrow(output)) .nm_stop("NONMEM simulation did not produce model.tab.")
    class(output) <- unique(c("nm_dataset", class(output)))
    attr(output, "solver") <- "NONMEM via PsN"
    output
  }
  if (audit) {
    attr(result, "audit_artifacts") <- .nm_nonmem_artifact_bundle(
      directory, operation, model, data
    )
  }
  result
}

.nm_nlmixr_replace <- function(source) {
  lines <- .nm_external_source_lines(source)
  lines <- gsub("\\bTHETA[[:space:]]*\\(([0-9]+)\\)", "theta\\1", lines,
                ignore.case = TRUE, perl = TRUE)
  lines <- gsub("\\bETA[[:space:]]*\\(([0-9]+)\\)", "eta\\1", lines,
                ignore.case = TRUE, perl = TRUE)
  lines <- gsub("\\bA[[:space:]]*\\(([0-9]+)\\)", "A\\1", lines,
                ignore.case = TRUE, perl = TRUE)
  lines <- gsub("(?<![<>=])=(?!=)", "<-", lines, perl = TRUE)
  lines
}

.nm_nlmixr_preflight <- function(model) {
  unsupported <- c(
    HMM_CONFIG = "hidden/semi-Markov models",
    KALMAN_CONFIG = "state-space/Kalman models",
    DDE_CONFIG = "delay differential equations",
    DAE_CONFIG = "differential-algebraic equations",
    RE_CONFIG = "nested random-effect declarations",
    COMPONENTS = "composed component models",
    OUTCOMES = "multi-outcome likelihood declarations",
    EXPERIMENTAL = "experimental engine families"
  )
  present <- names(unsupported)[vapply(names(unsupported), function(name) {
    value <- model[[name]]
    if (identical(name, "EXPERIMENTAL")) return(isTRUE(value$enabled))
    !is.null(value) && length(value) > 0L
  }, logical(1))]
  if (length(present)) {
    labels <- unname(unsupported[present])
    .nm_stop(
      "The nlmixr2 adapter does not yet provide a semantics-preserving translation for ",
      paste(labels, collapse = ", "), ". Use the LibeR engine for this model."
    )
  }
  if (identical(model$LIK_CONFIG$error %||% "residual", "likelihood")) {
    .nm_stop(
      "The nlmixr2 adapter does not yet translate custom likelihood models. ",
      "Use the LibeR engine for this model."
    )
  }
  invisible(model)
}

.nm_nlmixr_toolchain_guard <- function() {
  variables <- c("PATH", "BINPREF", "rxBINPREF")
  previous <- Sys.getenv(variables, unset = NA_character_)
  restore <- function() {
    for (name in variables) {
      if (is.na(previous[[name]])) {
        Sys.unsetenv(name)
      } else {
        do.call(Sys.setenv, stats::setNames(list(previous[[name]]), name))
      }
    }
    invisible(NULL)
  }
  if (.Platform$OS.type != "windows") return(restore)
  compiler <- unname(Sys.which("gcc"))
  make <- unname(Sys.which("make"))
  if (!nzchar(compiler) || !nzchar(make)) return(restore)
  compiler <- normalizePath(compiler, winslash = "/", mustWork = TRUE)
  make <- normalizePath(make, winslash = "/", mustWork = TRUE)
  compiler_directory <- dirname(compiler)
  make_directory <- dirname(make)
  paths <- strsplit(Sys.getenv("PATH"), .Platform$path.sep, fixed = TRUE)[[1L]]
  paths <- paths[nzchar(paths)]
  normalized <- vapply(paths, function(path) {
    if (!dir.exists(path)) return(path)
    normalizePath(path, winslash = "/", mustWork = TRUE)
  }, character(1))
  # An older compiler's cc1/ld on PATH can be selected by rxode2's nested
  # R CMD SHLIB invocation even while the visible gcc executable is Rtools.
  # Remove only such conflicting helper directories for this execution.
  conflicts <- vapply(normalized, function(path) {
    !identical(tolower(path), tolower(compiler_directory)) &&
      any(file.exists(file.path(path, c("cc1.exe", "ld.exe"))))
  }, logical(1))
  paths <- unique(c(compiler_directory, make_directory, paths[!conflicts]))
  prefix <- paste0(compiler_directory, "/")
  Sys.setenv(
    PATH = paste(paths, collapse = .Platform$path.sep),
    BINPREF = prefix,
    rxBINPREF = prefix
  )
  restore
}

.nm_nlmixr_error <- function(model) {
  error <- toupper(gsub("[[:space:]]+", "", model$ERROR))
  rhs <- sub(".*Y=", "", error)
  terms <- character()
  for (index in seq_len(nrow(model$SIGMAS))) {
    token <- paste0("ERR\\(", index, "\\)")
    if (!grepl(token, rhs, perl = TRUE)) next
    name <- paste0("sigma", index)
    if (grepl(paste0("EXP\\(", token, "\\)"), rhs, perl = TRUE)) {
      terms <- c(terms, paste0("lnorm(", name, ")"))
    } else if (grepl(paste0("F\\*", token, "|", token, "\\*F|F\\*\\(1\\+", token),
                     rhs, perl = TRUE)) {
      terms <- c(terms, paste0("prop(", name, ")"))
    } else if (grepl(token, rhs, perl = TRUE)) {
      terms <- c(terms, paste0("add(", name, ")"))
    }
  }
  if (!length(terms)) {
    if (rhs %in% c("F", "IPRED")) return("add(sigma_fixed_zero)")
    .nm_stop(
      "The nlmixr2 adapter could not translate this $ERROR expression. ",
      "Use an additive, proportional, combined, or log-normal ERR model."
    )
  }
  paste(terms, collapse = " + ")
}

.nm_nlmixr_ini <- function(model) {
  theta <- vapply(seq_len(nrow(model$THETAS)), function(index) {
    row <- model$THETAS[index, , drop = FALSE]
    value <- format(row$Value, digits = 17, scientific = TRUE)
    if (isTRUE(row$FIX)) paste0("theta", index, " <- fixed(", value, ")") else
      paste0(
        "theta", index, " <- c(",
        format(row$LOWER, digits = 17, scientific = TRUE), ",", value, ",",
        format(row$UPPER, digits = 17, scientific = TRUE), ")"
      )
  }, character(1))
  sigma <- vapply(seq_len(nrow(model$SIGMAS)), function(index) {
    row <- model$SIGMAS[index, , drop = FALSE]
    value <- sqrt(max(0, row$Value))
    paste0(
      "sigma", index, " <- ",
      if (isTRUE(row$FIX)) paste0("fixed(", format(value, digits = 17), ")") else
        format(value, digits = 17)
    )
  }, character(1))
  if (!length(sigma) && toupper(gsub("[[:space:]]", "", model$ERROR)) %in%
      c("Y=F", "Y=IPRED")) sigma <- "sigma_fixed_zero <- fixed(0)"
  omega <- character()
  n_eta <- .nm_n_eta(model$OMEGAS)
  if (n_eta) {
    matrix <- .nm_omega_matrix(model)
    correlated <- any(model$OMEGAS$ROW != model$OMEGAS$COL)
    if (correlated) {
      if (any(model$OMEGAS$FIX)) {
        .nm_stop("The nlmixr2 adapter does not yet translate partially fixed correlated OMEGA blocks.")
      }
      values <- unlist(lapply(seq_len(n_eta), function(row) matrix[row, seq_len(row)]))
      omega <- paste0(
        paste(paste0("eta", seq_len(n_eta)), collapse = " + "), " ~ c(",
        paste(format(values, digits = 17, scientific = TRUE), collapse = ","), ")"
      )
    } else {
      omega <- vapply(seq_len(n_eta), function(index) {
        row <- model$OMEGAS[
          model$OMEGAS$ROW == index & model$OMEGAS$COL == index, , drop = FALSE
        ]
        value <- format(row$Value[[1L]], digits = 17, scientific = TRUE)
        paste0("eta", index, " ~ ", if (isTRUE(row$FIX[[1L]])) "fixed(" else "",
               value, if (isTRUE(row$FIX[[1L]])) ")" else "")
      }, character(1))
    }
  }
  c(theta, sigma, omega)
}

.nm_nlmixr_model_lines <- function(model) {
  pred <- .nm_nlmixr_replace(model$PRED)
  pred <- pred[!grepl("^[[:space:]]*F[[:space:]]*<-", pred, ignore.case = TRUE)]
  residual <- .nm_nlmixr_error(model)
  if (identical(model$PRED_MODE %||% "pk", "pred")) {
    direct <- .nm_nlmixr_replace(model$PRED_SOURCE %||% model$PRED)
    direct <- sub("^[[:space:]]*F[[:space:]]*<-", "ipred <-", direct,
                  ignore.case = TRUE)
    return(c(direct, paste0("ipred ~ ", residual)))
  }
  standard <- model$ADVAN %in% c(1L, 2L, 3L, 4L, 11L, 12L) &&
    !nzchar(trimws(model$DES %||% ""))
  if (standard) return(c(pred, paste0("linCmt() ~ ", residual)))
  if (!nzchar(trimws(model$DES %||% ""))) {
    .nm_stop("The nlmixr2 adapter requires a supported solved ADVAN or a $DES block.")
  }
  des <- .nm_nlmixr_replace(model$DES)
  des <- sub("^[[:space:]]*DADT\\(([0-9]+)\\)[[:space:]]*<-", "d/dt(A\\1) <-", des,
             ignore.case = TRUE, perl = TRUE)
  scale <- paste0("S", model$OBSCMP)
  scale_defined <- any(grepl(paste0("^[[:space:]]*", scale, "[[:space:]]*<-"), pred,
                             ignore.case = TRUE, perl = TRUE))
  endpoint <- paste0(
    "ipred <- A", model$OBSCMP, if (scale_defined) paste0(" / ", scale) else ""
  )
  c(pred, des, endpoint, paste0("ipred ~ ", residual))
}

.nm_nlmixr_ui_text <- function(model) {
  .nm_nlmixr_preflight(model)
  ini <- .nm_nlmixr_ini(model)
  body <- .nm_nlmixr_model_lines(model)
  paste(c(
    "function() {", "  ini({", paste0("    ", ini), "  })",
    "  model({", paste0("    ", body), "  })", "}"
  ), collapse = "\n")
}

.nm_nlmixr_ui <- function(model) {
  text <- .nm_nlmixr_ui_text(model)
  model_function <- tryCatch(
    eval(parse(text = text, keep.source = TRUE), envir = asNamespace("nlmixr2")),
    error = function(error) .nm_stop(
      "Unable to construct the translated nlmixr2 model: ", conditionMessage(error)
    )
  )
  attr(model_function, "liberation_source") <- text
  model_function
}

.nm_nlmixr_control <- function(method, arguments) {
  method <- toupper(method)
  covariance <- isTRUE(arguments$covariance)
  cov_method <- if (!covariance) "" else switch(
    tolower(arguments$covariance_type %||% "sandwich"),
    sandwich = "r,s", hessian = "r", opg = "s", "r,s"
  )
  if (method == "SAEM") {
    return(list(
      est = "saem",
        control = nlmixr2est::saemControl(
        seed = as.integer(arguments$method_seed %||% 20260713L),
        nBurn = max(0L, as.integer(arguments$burn %||% 60L)),
        nEm = max(1L, as.integer(arguments$n_iter %||% 200L) -
                    max(0L, as.integer(arguments$burn %||% 60L))),
        print = max(0L, as.integer(arguments$print_every %||% 0L)),
        covMethod = cov_method
      )
    ))
  }
  if (method %in% c("FOCE", "FOCEI")) {
    return(list(
      est = "focei",
      control = nlmixr2est::foceiControl(
        maxOuterIterations = max(1L, as.integer(arguments$maxit %||% 200L)),
        maxInnerIterations = max(1L, as.integer(arguments$eta_maxit %||% 100L)),
        print = max(0L, as.integer(arguments$print_every %||% 0L)),
        covMethod = cov_method
      )
    ))
  }
  if (method == "LAPLACE" && exists("laplaceControl", asNamespace("nlmixr2est"),
                                     inherits = FALSE)) {
    constructor <- get("laplaceControl", asNamespace("nlmixr2est"), inherits = FALSE)
    return(list(est = "laplace", control = constructor(
      maxOuterIterations = max(1L, as.integer(arguments$maxit %||% 200L)),
      print = max(0L, as.integer(arguments$print_every %||% 0L)),
      covMethod = cov_method
    )))
  }
  .nm_stop("nlmixr2 execution currently supports FOCE/FOCEI, SAEM, and available LAPLACE builds.")
}

.nm_nlmixr_fit <- function(fit, model, data, method, elapsed) {
  fixed <- tryCatch(as.numeric(nlmixr2est::fixef(fit)), error = function(error) numeric())
  fixed_names <- tryCatch(names(nlmixr2est::fixef(fit)), error = function(error) character())
  names(fixed) <- fixed_names
  theta <- vapply(seq_len(nrow(model$THETAS)), function(index) {
    fixed[[paste0("theta", index)]] %||% model$THETAS$Value[[index]]
  }, numeric(1))
  sigma <- vapply(seq_len(nrow(model$SIGMAS)), function(index) {
    value <- fixed[[paste0("sigma", index)]] %||% sqrt(model$SIGMAS$Value[[index]])
    as.numeric(value)^2
  }, numeric(1))
  omega_matrix <- tryCatch(as.matrix(fit$omega), error = function(error) NULL)
  omega <- vapply(seq_len(nrow(model$OMEGAS)), function(index) {
    row <- model$OMEGAS$ROW[[index]]
    column <- model$OMEGAS$COL[[index]]
    if (!is.null(omega_matrix) && nrow(omega_matrix) >= row && ncol(omega_matrix) >= column) {
      omega_matrix[row, column]
    } else model$OMEGAS$Value[[index]]
  }, numeric(1))
  eta <- tryCatch(as.matrix(fit$eta), error = function(error) NULL)
  if (is.null(eta)) eta <- matrix(numeric(), length(unique(data$ID)), 0L)
  if (ncol(eta)) {
    eta_columns <- grep("^eta", colnames(eta), ignore.case = TRUE)
    if (length(eta_columns)) eta <- eta[, eta_columns, drop = FALSE]
    colnames(eta) <- paste0("ETA", seq_len(ncol(eta)))
  }
  if (nrow(eta) != length(unique(data$ID))) {
    .nm_stop("nlmixr2 ETA subject count does not match the submitted dataset.")
  }
  output <- as.data.frame(fit, stringsAsFactors = FALSE, check.names = FALSE)
  names(output) <- toupper(names(output))
  objective <- tryCatch(as.numeric(fit$objf), error = function(error) NA_real_)
  if (!length(objective) || !is.finite(objective[[1L]])) objective <- NA_real_
  raw_covariance <- tryCatch(as.matrix(fit$cov), error = function(error) NULL)
  covariance <- if (is.null(raw_covariance) || !length(raw_covariance)) NULL else list(
    status = "completed", type = "nlmixr2", covariance = raw_covariance,
    correlation = tryCatch(stats::cov2cor(raw_covariance), error = function(error) NULL),
    bread_source = "nlmixr2 output", bread_exact = NA
  )
  iterations <- tryCatch(nrow(as.data.frame(fit$parHist)), error = function(error) NA_integer_)
  structure(list(
    version = 1L, method = toupper(method), execution_engine = "nlmixr2",
    objective = objective[[1L]], theta = theta, omega = omega, sigma = sigma,
    eta = eta, convergence = 0L, message = "nlmixr2 estimation completed",
    iterations = as.integer(iterations), objective_evaluations = NA_integer_,
    evaluations = stats::setNames(
      c(NA_integer_, NA_integer_), c("function", "gradient")
    ),
    model = model, data = data,
    diagnostics = list(
      population_gradient = "nlmixr2 estimator-specific derivative path",
      optimizer = list(backend = "nlmixr2", objective_backend = "nlmixr2")
    ), covariance = covariance,
    timing = list(
      model_fit_seconds = as.numeric(elapsed),
      covariance_seconds = NA_real_, total_seconds = as.numeric(elapsed)
    ), output = output
  ), class = "nm_fit")
}

.nm_nlmixr_audit <- function(result, model, data, operation, source) {
  details <- list(engine_label = "nlmixr2")
  bundle <- .nm_audit_bundle(
    result, model, data, operation, details = details,
    engine = "nlmixr2", source = "generated"
  )
  bundle$files[[length(bundle$files) + 1L]] <- .nm_audit_file("model.R", c(
    "# Generated nlmixr2 model executed by LibeRation.", source
  ), "text/x-r-source")
  .nm_audit_bundle_validate(bundle)
  bundle
}

.nm_nlmixr_run <- function(model, data, operation, arguments) {
  if (!requireNamespace("nlmixr2", quietly = TRUE) ||
      !requireNamespace("nlmixr2est", quietly = TRUE) ||
      !requireNamespace("rxode2", quietly = TRUE)) {
    .nm_stop(
      "nlmixr2 execution requires the optional `nlmixr2`, `nlmixr2est`, ",
      "and `rxode2` packages."
    )
  }
  restore_toolchain <- .nm_nlmixr_toolchain_guard()
  on.exit(restore_toolchain(), add = TRUE)
  data <- nm_dataset(data)
  ui <- .nm_nlmixr_ui(model)
  source <- attr(ui, "liberation_source", exact = TRUE)
  started <- proc.time()[["elapsed"]]
  if (identical(operation, "estimate")) {
    if (length(arguments$stages %||% list()) > 1L) {
      .nm_stop("Sequential nlmixr2 stages are not yet supported; submit each stage as a separate run.")
    }
    method <- toupper(arguments$method %||% "FOCEI")
    controls <- .nm_nlmixr_control(method, arguments)
    fit <- nlmixr2est::nlmixr2(
      ui, as.data.frame(data), est = controls$est, control = controls$control,
      table = nlmixr2est::tableControl(
        cwres = TRUE, addDosing = TRUE,
        cores = max(1L, as.integer(arguments$n_cores %||% 1L))
      )
    )
    elapsed <- unname(proc.time()[["elapsed"]] - started)
    result <- .nm_nlmixr_fit(fit, model, data, method, elapsed)
  } else {
    if (!is.null(arguments$seed)) set.seed(as.integer(arguments$seed))
    simulation_ui <- ui
    zero_re <- get0("zeroRe", envir = asNamespace("rxode2"), inherits = FALSE)
    random_effects <- isTRUE(arguments$random_effects %||% TRUE)
    residual <- isTRUE(arguments$residual %||% TRUE)
    zero <- c(if (!random_effects) "omega", if (!residual) "sigma")
    if (is.function(zero_re) && length(zero)) {
      simulation_ui <- zero_re(simulation_ui, which = zero, fix = TRUE)
    }
    solved <- rxode2::rxSolve(
      simulation_ui, as.data.frame(data),
      nStud = max(1L, as.integer(arguments$nsim %||% 1L)),
      cores = max(1L, as.integer(arguments$n_cores %||% 1L))
    )
    result <- as.data.frame(solved, stringsAsFactors = FALSE, check.names = FALSE)
    names(result) <- toupper(names(result))
    if ("SIM" %in% names(result) && !"DV" %in% names(result)) result$DV <- result$SIM
    if ("SIM.ID" %in% names(result)) result$SIM <- result$SIM.ID
    class(result) <- unique(c("nm_dataset", class(result)))
    attr(result, "solver") <- "nlmixr2/rxode2"
  }
  if (isTRUE(arguments$audit_artifacts)) {
    attr(result, "audit_artifacts") <- .nm_nlmixr_audit(
      result, model, data, operation, source
    )
  }
  result
}

.nm_engine_run <- function(engine, operation, model, data, arguments) {
  engine <- .nm_execution_engine(engine)
  operation <- match.arg(operation, c("estimate", "simulate"))
  if (engine == "liber") {
    return(if (operation == "estimate") {
      if (!is.null(arguments$stages)) {
        do.call(nm_est_sequence, c(list(model = model, data = data), arguments))
      } else do.call(nm_est, c(list(model = model, data = data), arguments))
    } else do.call(nm_simulate, c(list(model = model, data = data), arguments)))
  }
  if (!is.null(arguments$theta)) {
    if (length(arguments$theta) != nrow(model$THETAS)) {
      .nm_stop("External simulation THETA values have the wrong length.")
    }
    model$THETAS$Value <- as.numeric(arguments$theta)
    arguments$theta <- NULL
  }
  if (!is.null(arguments$omega)) {
    if (length(arguments$omega) != nrow(model$OMEGAS)) {
      .nm_stop("External simulation OMEGA values have the wrong length.")
    }
    model$OMEGAS$Value <- as.numeric(arguments$omega)
    arguments$omega <- NULL
  }
  if (!is.null(arguments$sigma)) {
    if (length(arguments$sigma) != nrow(model$SIGMAS)) {
      .nm_stop("External simulation SIGMA values have the wrong length.")
    }
    model$SIGMAS$Value <- as.numeric(arguments$sigma)
    arguments$sigma <- NULL
  }
  if (engine == "nonmem") return(.nm_nonmem_run(model, data, operation, arguments))
  .nm_nlmixr_run(model, data, operation, arguments)
}

#' Execute a model using NONMEM or nlmixr2
#'
#' External commands are selected from a fixed allow-list. NONMEM is invoked
#' through the administrator-configured PsN `execute` command; no executable
#' path or shell fragment is accepted from a remote job. nlmixr2 model code is
#' generated only from the already validated `nm_model` semantic contract.
#'
#' @param model An [nm_model()].
#' @param data NONMEM-style event data.
#' @param engine `"nonmem"` or `"nlmixr2"`.
#' @param operation `"estimate"` or `"simulate"`.
#' @param ... Typed estimation or simulation controls accepted by the adapter.
#' @return An `nm_fit` for estimation or an `nm_dataset` for simulation.
#' @export
nm_external_run <- function(model, data,
                            engine = c("nonmem", "nlmixr2"),
                            operation = c("estimate", "simulate"), ...) {
  engine <- match.arg(engine)
  operation <- match.arg(operation)
  if (!inherits(model, "nm_model")) .nm_stop("`model` must be an nm_model.")
  .nm_engine_run(engine, operation, model, data, list(...))
}
