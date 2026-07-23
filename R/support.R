.liber_redact_text <- function(value) {
  value <- enc2utf8(as.character(value))
  value <- gsub(
    "(?i)\\b[A-Z0-9._%+-]+@[A-Z0-9.-]+\\.[A-Z]{2,}\\b",
    "[REDACTED_EMAIL]", value, perl = TRUE
  )
  value <- gsub(
    "(?i)\\b(Bearer)\\s+[A-Z0-9._~+/-]+=*",
    "\\1 [REDACTED_TOKEN]", value, perl = TRUE
  )
  value <- gsub(
    "(?i)\\b(token|password|passwd|secret|api[_-]?key)\\b(\\s*[:=]\\s*)[^\\s,;]+",
    "\\1\\2[REDACTED_SECRET]", value, perl = TRUE
  )
  value <- gsub(
    "(?i)([A-Z]:[/\\\\]Users[/\\\\])[^/\\\\[:space:]]+",
    "\\1[REDACTED_USER]", value, perl = TRUE
  )
  value <- gsub(
    "(/(?:home|Users)/)[^/[:space:]]+",
    "\\1[REDACTED_USER]", value, perl = TRUE
  )
  value
}

.liber_redact_object <- function(value) {
  if (is.character(value)) return(.liber_redact_text(value))
  if (is.factor(value)) return(factor(.liber_redact_text(as.character(value))))
  if (is.data.frame(value)) {
    value[] <- lapply(value, .liber_redact_object)
    return(value)
  }
  if (is.list(value)) return(lapply(value, .liber_redact_object))
  value
}

.liber_dataset_summary <- function(data) {
  if (is.null(data)) return(NULL)
  if (!is.data.frame(data)) .nm_stop("`data` must be a data frame or `nm_dataset`.")
  classes <- vapply(data, function(value) paste(class(value), collapse = "/"), character(1))
  missing <- vapply(data, function(value) sum(is.na(value)), integer(1))
  list(
    rows = nrow(data),
    columns = ncol(data),
    column_names = names(data),
    column_classes = as.list(classes),
    missing_values = as.list(missing),
    subjects = if ("ID" %in% names(data)) length(unique(data$ID)) else NA_integer_,
    observations = if (all(c("EVID", "MDV") %in% names(data))) {
      sum(data$EVID == 0L & data$MDV == 0L, na.rm = TRUE)
    } else {
      NA_integer_
    },
    generated_events = if (".generated" %in% names(data)) {
      sum(data$.generated, na.rm = TRUE)
    } else {
      NA_integer_
    },
    content_sha256 = digest::digest(data, algo = "sha256", serialize = TRUE)
  )
}

.liber_model_summary <- function(model, include_code = FALSE) {
  if (is.null(model)) return(NULL)
  if (!inherits(model, "nm_model")) .nm_stop("`model` must be an `nm_model`.")
  result <- list(
    ADVAN = model$ADVAN,
    TRANS = model$TRANS,
    solver = model$SOLVER,
    language = model$LANGUAGE,
    theta_count = nrow(model$THETAS %||% data.frame()),
    eta_count = as.integer(model$n_eta %||% 0L),
    omega_parameter_count = nrow(model$OMEGAS %||% data.frame()),
    sigma_count = nrow(model$SIGMAS %||% data.frame()),
    outputs = nm_model_outputs(model),
    experimental = model$EXPERIMENTAL %||% list(enabled = FALSE),
    content_sha256 = digest::digest(
      nm_model_to_contract(model), algo = "sha256", serialize = TRUE
    )
  )
  if (isTRUE(include_code)) {
    result$code <- list(PRED = model$PRED %||% "", DES = model$DES %||% "",
                        ERROR = model$ERROR %||% "")
  }
  result
}

.liber_fit_summary <- function(fit) {
  if (is.null(fit)) return(NULL)
  if (!inherits(fit, "nm_fit")) .nm_stop("`fit` must be an `nm_fit`.")
  list(
    method = fit$method %||% "",
    convergence = fit$convergence %||% NA_integer_,
    message = fit$message %||% "",
    iterations = fit$iterations %||% NA_integer_,
    objective_evaluations = fit$objective_evaluations %||% NA_integer_,
    timing = fit$timing %||% list(),
    diagnostic_sections = names(fit$diagnostics %||% list()),
    model_sha256 = if (inherits(fit$model, "nm_model")) {
      digest::digest(nm_model_to_contract(fit$model), algo = "sha256", serialize = TRUE)
    } else {
      NA_character_
    }
  )
}

#' Report ecosystem capability and evidence tiers
#'
#' Returns the versioned support declaration for all LibeR packages. A
#' `validated` capability has an independent analytic or external-software
#' comparison, `verified` denotes deterministic internal and recovery testing,
#' and `experimental` denotes a research implementation without sufficient
#' qualification for routine use.
#'
#' @param package Optional package name or names to retain.
#' @param evidence_tier Optional tier or tiers to retain.
#' @return A data frame with capability status, evidence, gate, and recommended
#'   use.
#' @export
liber_support_matrix <- function(package = NULL, evidence_tier = NULL) {
  path <- system.file("ecosystem", "support-matrix.csv", package = "LibeRation")
  if (!nzchar(path) || !file.exists(path)) {
    path <- file.path("inst", "ecosystem", "support-matrix.csv")
  }
  if (!file.exists(path)) .nm_stop("The installed LibeR support matrix is unavailable.")
  result <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!is.null(package)) result <- result[result$package %in% package, , drop = FALSE]
  if (!is.null(evidence_tier)) {
    result <- result[result$evidence_tier %in% evidence_tier, , drop = FALSE]
  }
  rownames(result) <- NULL
  result
}

#' Create a redacted LibeR diagnostic support bundle
#'
#' The default archive contains compatibility, compiled-engine, runtime,
#' structural model/data, and fit metadata. It never contains data values,
#' parameter estimates, ETAs, environment variables, credentials, workspace
#' contents, or model code. Model code and redacted logs are opt-in and should
#' still be inspected before sharing.
#'
#' @param path Output `.zip` path.
#' @param workspace Optional [nm_workspace()] or workspace directory to verify.
#' @param model Optional `nm_model` to describe structurally.
#' @param data Optional data frame to describe without retaining values.
#' @param fit Optional `nm_fit` to describe without estimates or embedded data.
#' @param logs Optional character vector of log file paths or log text.
#' @param include_model_code Include `$PK/$PRED`, `$DES`, and `$ERROR` source.
#' @param include_logs Include supplied logs after secret/path redaction.
#' @param overwrite Replace an existing archive.
#' @return The normalized archive path, invisibly.
#' @export
liber_support_bundle <- function(
    path = file.path(getwd(), paste0("LibeR-support-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".zip")),
    workspace = NULL, model = NULL, data = NULL, fit = NULL, logs = NULL,
    include_model_code = FALSE, include_logs = FALSE, overwrite = FALSE) {
  path <- path.expand(path)
  if (!grepl("\\.zip$", path, ignore.case = TRUE)) path <- paste0(path, ".zip")
  parent <- dirname(path)
  if (!dir.exists(parent)) dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(parent)) .nm_stop("Unable to create support-bundle directory.")
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (file.exists(path) && !isTRUE(overwrite)) {
    .nm_stop("Support bundle already exists; set `overwrite = TRUE` to replace it.")
  }

  if (is.null(model) && inherits(fit, "nm_fit")) model <- fit$model
  if (is.null(data) && inherits(fit, "nm_fit")) data <- fit$data
  diagnostics <- liber_doctor(workspace = workspace, verbose = FALSE)
  diagnostics$workspace$path <- NULL
  diagnostics$runtime$make <- basename(diagnostics$runtime$make %||% "")
  diagnostics$runtime$compiler <- basename(diagnostics$runtime$compiler %||% "")
  workspace_verification <- if (!is.null(workspace)) {
    tryCatch({
      checked <- nm_workspace_verify(workspace)
      issues <- checked$issues %||% data.frame()
      if (nrow(issues)) issues$path <- NULL
      list(
        valid = checked$valid,
        projects = checked$projects,
        snapshots = checked$snapshots,
        objects = checked$objects,
        issues = issues
      )
    }, error = function(error) list(valid = FALSE, errors = conditionMessage(error)))
  } else {
    NULL
  }
  manifest <- list(
    schema = "liber.support-bundle/1",
    created_at_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC"),
    redaction = list(
      data_values_included = FALSE,
      parameter_estimates_included = FALSE,
      environment_variables_included = FALSE,
      model_code_included = isTRUE(include_model_code),
      logs_included = isTRUE(include_logs) && length(logs) > 0L
    ),
    diagnostics = unclass(diagnostics),
    workspace = workspace_verification,
    model = .liber_model_summary(model, include_code = include_model_code),
    data = .liber_dataset_summary(data),
    fit = .liber_fit_summary(fit),
    support_matrix_sha256 = digest::digest(
      liber_support_matrix(), algo = "sha256", serialize = TRUE
    )
  )
  manifest <- .liber_redact_object(manifest)

  stage <- tempfile("liber-support-")
  dir.create(stage, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(stage, recursive = TRUE, force = TRUE), add = TRUE)
  jsonlite::write_json(
    manifest, file.path(stage, "manifest.json"), auto_unbox = TRUE,
    pretty = TRUE, null = "null", na = "null", digits = 17
  )
  writeLines(
    .liber_redact_text(utils::capture.output(utils::sessionInfo())),
    file.path(stage, "session-info.txt"), useBytes = TRUE
  )
  utils::write.csv(
    liber_support_matrix(), file.path(stage, "support-matrix.csv"),
    row.names = FALSE, na = ""
  )
  writeLines(
    c(
      "LibeR redacted support bundle",
      "",
      "Inspect every file before sharing.",
      "No dataset values, parameter estimates, ETAs, environment variables,",
      "credentials, workspace contents, or model code are included by default."
    ),
    file.path(stage, "README.txt"), useBytes = TRUE
  )
  if (isTRUE(include_logs) && length(logs)) {
    log_lines <- unlist(lapply(logs, function(value) {
      if (length(value) == 1L && file.exists(value)) {
        readLines(value, warn = FALSE, encoding = "UTF-8")
      } else {
        as.character(value)
      }
    }), use.names = FALSE)
    writeLines(.liber_redact_text(log_lines), file.path(stage, "redacted-logs.txt"),
               useBytes = TRUE)
  }
  if (file.exists(path)) unlink(path)
  zip::zipr(path, list.files(stage, full.names = TRUE), root = stage)
  if (!file.exists(path)) .nm_stop("Support bundle could not be created.")
  invisible(path)
}
