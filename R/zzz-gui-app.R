.liber_report_default_directory <- function(workspace, project = NULL) {
  root <- file.path(.nm_workspace_path(workspace), "reports")
  if (!is.null(project) && length(project) && nzchar(as.character(project)[[1L]])) {
    root <- file.path(root, as.character(project)[[1L]])
  }
  normalizePath(root, winslash = "/", mustWork = FALSE)
}

.liber_report_choose_directory <- function(initial = getwd()) {
  initial <- path.expand(as.character(initial %||% getwd())[[1L]])
  candidate <- initial
  while (!dir.exists(candidate) && !identical(dirname(candidate), candidate)) {
    candidate <- dirname(candidate)
  }
  if (!dir.exists(candidate)) candidate <- getwd()
  selected <- if (.Platform$OS.type == "windows") {
    utils::choose.dir(default = candidate, caption = "Choose report output folder")
  } else if (capabilities("tcltk") && requireNamespace("tcltk", quietly = TRUE)) {
    as.character(tcltk::tk_choose.dir(default = candidate,
                                      caption = "Choose report output folder"))
  } else {
    .nm_stop("A native folder chooser is unavailable; enter the report directory manually.")
  }
  if (!length(selected) || is.na(selected) || !nzchar(selected)) return(NULL)
  normalizePath(selected, winslash = "/", mustWork = TRUE)
}

#' Launch the LibeR React modelling application
#'
#' This is the interactive workbench used to manage model versions, datasets,
#' simulations, estimations, diagnostics, reports, and local or remote jobs.
#'
#' @param model Optional model loaded into the workbench.
#' @param data Optional dataset loaded into the workbench.
#' @param queue Optional `LibeRQueue`; simulation and estimation jobs are
#'   submitted to it. `NULL` creates a persistent local queue under the
#'   workspace when LibeRties is installed. Use `FALSE` for in-process work.
#' @param workspace Workspace object or directory. By default the GUI uses
#'   `C:/Users/<username>/Documents/LibeR/workspace` on Windows and
#'   `~/LibeR/workspace` on Linux and macOS.
#' @param session_workspace Create a separate ephemeral workspace beneath
#'   `workspace` for every browser session. This is intended for hosted
#'   demonstrations where application users must not share project files.
#' @param project Optional project id to open when the application starts.
#' @param launch.browser Passed to [shiny::runApp()]. Use `NULL` to return the
#'   Shiny app object without launching it.
#' @param ... Additional arguments passed to [shiny::runApp()].
#' @examples
#' if (interactive()) liber_gui()
#' @export
liber_gui <- function(model = NULL, data = NULL, queue = NULL,
                      workspace = NULL, project = NULL, session_workspace = FALSE,
                      launch.browser = getOption("shiny.launch.browser", interactive()), ...) {
  dots <- list(...)
  model_input <- model
  data_input <- data
  queue_input <- queue
  workspace_input <- workspace
  project_input <- project
  favicon <- system.file("assets", "favicon.svg", package = "LibeRation")
  favicon_href <- if (nzchar(favicon) && file.exists(favicon)) {
    prefix <- paste0(
      "liberation-assets-",
      gsub("[^A-Za-z0-9_-]", "-", as.character(utils::packageVersion("LibeRation")))
    )
    if (!prefix %in% names(shiny::resourcePaths())) {
      shiny::addResourcePath(prefix, dirname(favicon))
    }
    paste0(prefix, "/favicon.svg")
  } else ""
  widget_directory <- system.file("htmlwidgets", package = "LibeRation")
  ai_worker_href <- ""
  if (nzchar(widget_directory) && dir.exists(widget_directory)) {
    widget_prefix <- paste0(
      "liberation-widget-",
      gsub("[^A-Za-z0-9_-]", "-", as.character(utils::packageVersion("LibeRation")))
    )
    if (!widget_prefix %in% names(shiny::resourcePaths())) {
      shiny::addResourcePath(widget_prefix, widget_directory)
    }
    ai_worker_href <- paste0(widget_prefix, "/liber-ai-worker.js")
  }
  ai_models <- .liber_ai_models()
  allowed_ai_models <- vapply(ai_models, `[[`, character(1), "id")
  ui <- .liber_full_page_ui(
    htmltools::tags$head(
      htmltools::tags$title("LibeRation"),
      htmltools::tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
      htmltools::tags$script(htmltools::HTML(
        "(function(){try{var t=localStorage.getItem('liber.theme');if(t!=='dark'&&t!=='light'){var l=localStorage.getItem('liberationDarkTheme');t=l==='1'?'dark':l==='0'?'light':(matchMedia('(prefers-color-scheme: dark)').matches?'dark':'light');}document.documentElement.setAttribute('data-liber-theme',t);}catch(e){}})();"
      )),
      if (nzchar(favicon_href)) {
        htmltools::tags$link(rel = "icon", type = "image/svg+xml", href = favicon_href)
      }
    ),
    liberWorkbenchOutput("workbench", height = "100vh")
  )

  server <- function(input, output, session) {
    workspace_path <- if (isTRUE(session_workspace)) {
      base <- if (inherits(workspace_input, "nm_workspace")) {
        workspace_input$path
      } else {
        workspace_input %||% file.path(tempdir(), "LibeRation-cloud")
      }
      file.path(base, "sessions", gsub("[^A-Za-z0-9_-]", "-", session$token))
    } else if (inherits(workspace_input, "nm_workspace")) {
      workspace_input$path
    } else {
      workspace_input %||% .liber_default_workspace()
    }
    workspace <- nm_workspace(workspace_path)
    queue <- queue_input
    if (identical(queue, FALSE)) {
      queue <- NULL
    } else if (is.null(queue) && requireNamespace("LibeRties", quietly = TRUE)) {
      queue <- LibeRties::ls_local_queue(
        root = file.path(workspace$path, ".jobs"), user = "local", max_workers = 1L
      )
    }
    client_settings <- .liber_client_settings_read(workspace)
    saved_remote_config <- client_settings$remotes
    saved_remote_queues <- list()
    saved_remote_meta <- list()
    if (length(saved_remote_config) && requireNamespace("LibeRties", quietly = TRUE)) {
      for (id in names(saved_remote_config)) {
        config <- saved_remote_config[[id]]
        remote <- tryCatch(
          LibeRties::ls_remote(config$url, config$token,
                               timeout = config$timeout %||% 30),
          error = function(error) NULL
        )
        if (is.null(remote)) next
        saved_remote_queues[[id]] <- remote
        saved_remote_meta[[id]] <- list(
          name = config$name %||% "Remote server", url = remote$url,
          user = config$user %||% ""
        )
      }
    }
    selected_queue <- as.character(client_settings$selected_queue %||% "local")[[1L]]
    if (!identical(selected_queue, "local") &&
        is.null(saved_remote_queues[[selected_queue]])) {
      selected_queue <- "local"
    }
    initial_jobs <- if (identical(selected_queue, "local") && !is.null(queue)) {
      tryCatch(
        as.data.frame(queue$list(), stringsAsFactors = FALSE),
        error = function(error) data.frame()
      )
    } else data.frame()
    initial_result <- NULL
    initial_snapshot <- NULL
    initial_run <- NULL
    initial_diagnostics <- list()
    session_model <- model_input
    session_data <- data_input
    if (!is.null(project_input)) {
      records <- nm_project_list(workspace, project_input)
      versions <- records[records$entry_type == "version", , drop = FALSE]
      if (nrow(versions)) {
        opened <- nm_project_load(workspace, project_input, versions$id[[1L]])
        session_model <- session_model %||% opened$model
        session_data <- session_data %||% opened$data
        initial_result <- opened$result
        initial_snapshot <- opened$id
      }
    }
    initial_fit <- if (inherits(initial_result, "nm_fit")) initial_result else NULL
    if (!client_settings$ai$help_model %in% allowed_ai_models) {
      client_settings$ai$help_model <- .liber_ai_default_help_model()
    }
    if (!client_settings$ai$report_model %in%
        c("same_as_help", allowed_ai_models)) {
      client_settings$ai$report_model <- .liber_ai_default_report_model()
    }
    client_settings$ai$model <- client_settings$ai$help_model
    saved_report_design <- function(project_id) {
      if (is.null(project_id) || !length(project_id) ||
          !nzchar(as.character(project_id)[[1L]])) return(NULL)
      metadata <- tryCatch(
        nm_report_design_load(workspace, project_id),
        error = function(error) NULL
      )
      if (is.null(metadata) || !nrow(metadata)) return(NULL)
      tryCatch(
        nm_report_design_load(workspace, project_id, metadata$id[[1L]]),
        error = function(error) NULL
      )
    }
    state <- shiny::reactiveValues(
      model = if (inherits(session_model, "NMEngine")) session_model$model else session_model,
      data = session_data, result = initial_result, fit = initial_fit,
      fit_payload = .liber_gui_fit(initial_fit, include_gof = FALSE),
      hmm_payload = NULL, kalman_payload = NULL,
      diagnostics = initial_diagnostics, report = NULL,
      report_design = saved_report_design(project_input),
      report_ai_context = list(
        available = FALSE, project = "", project_name = "", request_id = "",
        run_ids = character(), message = "Report evidence has not been requested.",
        runs = list()
      ),
      ai_config = client_settings$ai,
      ai_context = list(
        available = FALSE, project = "", project_name = "", request_id = "",
        scope = "index", message = "Project result summaries have not been requested.",
        runs = list()
      ),
      diagram_candidate = NULL,
      jobs = initial_jobs, job_log = character(), project = project_input,
      snapshot = initial_snapshot, run = initial_run,
      queue_id = selected_queue, remote_queues = saved_remote_queues,
      remote_meta = saved_remote_meta, remote_config = saved_remote_config,
      job_context = client_settings$pending_jobs %||% list(),
      hidden_jobs = character(), selected_job = NULL,
      active_page = "home", comparison_open = FALSE,
      draft_outputs = NULL, data_payload = FALSE, gof_payload = FALSE,
      diagnostic_payload = character(), refreshed = "Queue ready",
      log_level = "info", log_current = "Workbench ready",
      log_history = "Workbench ready"
    )
    poll_backoff <- new.env(parent = emptyenv())
    poll_backoff$until <- 0
    poll_backoff$ready <- FALSE
    active_queue <- function(queue_id = state$queue_id) {
      if (identical(queue_id, "local")) return(queue)
      state$remote_queues[[queue_id]]
    }
    with_ui_remote_timeout <- function(q, background = FALSE, operation) {
      if (!inherits(q, "LibeRRemote")) return(force(operation))
      original <- q$timeout
      on.exit(q$timeout <- original, add = TRUE)
      q$timeout <- min(original, if (isTRUE(background)) 2 else 5)
      force(operation)
    }
    job_context_key <- function(id, queue_id = state$queue_id) {
      paste(queue_id, as.character(id), sep = "::")
    }
    save_client_settings <- function() {
      persistent_jobs <- lapply(state$job_context, function(context) {
        context$model <- NULL
        context$data <- NULL
        context
      })
      .liber_client_settings_write(
        workspace, selected_queue = state$queue_id,
        remotes = state$remote_config, pending_jobs = persistent_jobs,
        ai = state$ai_config
      )
    }
    append_log <- function(message, level = "info") {
      message <- as.character(message)
      line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", message)
      state$log_level <- level
      state$log_current <- message
      state$log_history <- utils::tail(c(state$log_history, line), 500L)
      invisible(NULL)
    }
    update_fit <- function(fit) {
      state$fit <- if (inherits(fit, "nm_fit")) fit else NULL
      state$gof_payload <- FALSE
      state$hmm_payload <- NULL
      state$kalman_payload <- NULL
      state$fit_payload <- .liber_gui_fit(state$fit, include_gof = FALSE)
      invisible(state$fit)
    }
    diagram_from_event <- function(graph) {
      if (!is.list(graph)) .nm_stop("The visual model payload is invalid.")
      nm_model_diagram(
        compartments = graph$compartments,
        flows = graph$flows,
        parameters = graph$parameters,
        advan = as.integer(graph$advan %||% 6L),
        residual = as.character(graph$residual %||% "additive"),
        covariates = as.character(unlist(graph$covariates %||% character())),
        title = as.character(graph$title %||%
          attr(state$model, "name", exact = TRUE) %||% "Diagram model")
      )
    }
    report_design_from_event <- function(event) {
      blocks <- lapply(event$blocks %||% list(), function(block) {
        nm_report_block(
          type = as.character(block$type),
          title = as.character(block$title %||% NULL),
          source = as.character(block$source %||% "user"),
          text = as.character(block$text %||% ""),
          run_ids = as.character(unlist(block$run_ids %||% character())),
          elements = as.character(unlist(block$elements %||% character())),
          options = as.list(block$options %||% list()),
          id = as.character(block$id)
        )
      })
      nm_report_design(
        blocks = blocks,
        title = as.character(event$title %||% "LibeRation modelling report"),
        formats = as.character(unlist(event$formats %||% c("docx", "pdf"))),
        style = list(
          filename = as.character(event$name %||% "liberation-report"),
          output_directory = as.character(event$directory %||% "")
        ),
        id = if (nzchar(as.character(event$id %||% ""))) as.character(event$id) else NULL
      )
    }
    report_uploaded_text <- function(name, data) {
      name <- basename(as.character(name %||% "source.txt"))
      raw <- jsonlite::base64_dec(sub("^data:[^,]*,", "", as.character(data %||% "")))
      if (!length(raw)) .nm_stop("The uploaded report source is empty.")
      path <- tempfile(fileext = paste0(".", tolower(tools::file_ext(name))))
      on.exit(unlink(path, force = TRUE), add = TRUE)
      connection <- file(path, open = "wb")
      tryCatch(writeBin(raw, connection), finally = close(connection))
      extension <- tolower(tools::file_ext(path))
      if (extension == "pdf") {
        if (!requireNamespace("pdftools", quietly = TRUE)) {
          .nm_stop("Reading PDF source documents requires the optional `pdftools` package.")
        }
        return(paste(pdftools::pdf_text(path), collapse = "\n\n"))
      }
      .nm_report_document_text(path)
    }
    reset_lazy_payloads <- function(data = TRUE, diagnostics = TRUE) {
      if (isTRUE(data)) state$data_payload <- FALSE
      state$gof_payload <- FALSE
      state$hmm_payload <- NULL
      state$kalman_payload <- NULL
      if (isTRUE(diagnostics)) state$diagnostic_payload <- character()
      invisible(NULL)
    }
    invalidate_ai_context <- function() {
      state$ai_context <- list(
        available = FALSE, project = "", project_name = "", request_id = "",
        scope = "index", message = "Project result summaries have not been requested.",
        runs = list()
      )
      invisible(NULL)
    }
    record <- function(success, value, update_result = TRUE) {
      value <- tryCatch(force(value), error = identity)
      if (isTRUE(update_result) || inherits(value, "error")) state$result <- value
      if (inherits(value, "nm_fit")) update_fit(value)
      if (inherits(value, "nm_report") || inherits(value, "nm_report_bundle")) {
        state$report <- value
      }
      if (inherits(value, "error")) {
        append_log(conditionMessage(value), "error")
      } else {
        append_log(success, "info")
      }
      invisible(value)
    }
    ensure_parent_version <- function() {
      if (is.null(state$project)) return(NULL)
      if (is.null(state$snapshot)) {
        state$snapshot <- nm_project_save(
          workspace, state$project, state$model, state$data, NULL, label = NULL
        )
      }
      state$snapshot
    }
    materialize_job_result <- function(context, q, id) {
      records <- nm_project_list(workspace, context$project)
      existing <- if (nrow(records) && "queue_job_id" %in% names(records)) {
        which(records$entry_type == "run" &
                records$queue_id == context$queue_id &
                records$queue_job_id == context$job_id)
      } else integer()
      if (length(existing)) {
        context$run_id <- as.character(records$id[[existing[[1L]]]])
        context$materialize_error <- NULL
        context$materialize_attempted_at <- NULL
        return(context)
      }
      result <- q$result(id)
      source <- nm_project_load(
        workspace, context$project, snapshot = context$version
      )
      run_model <- context$model %||% if (inherits(result, "nm_fit")) {
        result$model
      } else source$model
      run_data <- context$data %||% if (inherits(result, "nm_fit")) {
        result$data
      } else source$data
      context$run_id <- nm_project_save_run(
        workspace, context$project, context$version, result,
        label = context$label, model = run_model, data = run_data,
        queue_id = context$queue_id, queue_job_id = context$job_id
      )
      if (identical(as.character(context$project), as.character(state$project))) {
        invalidate_ai_context()
      }
      context$materialize_error <- NULL
      context$materialize_attempted_at <- NULL
      context
    }
    materialization_due <- function(context, now = Sys.time()) {
      if (!identical(as.character(context$status %||% ""), "completed") ||
          nzchar(as.character(context$run_id %||% ""))) return(FALSE)
      attempted <- suppressWarnings(as.POSIXct(
        as.character(context$materialize_attempted_at %||% NA_character_),
        format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC"
      ))
      is.na(attempted) || as.numeric(difftime(now, attempted, units = "secs")) >= 30
    }
    reconcile_job_results <- function(jobs, q) {
      if (is.null(q) || !nrow(jobs) ||
          !all(c("id", "status") %in% names(jobs))) return(invisible(FALSE))
      contexts <- state$job_context
      changed <- FALSE
      for (row in seq_len(nrow(jobs))) {
        id <- as.character(jobs$id[[row]])
        key <- job_context_key(id)
        context <- contexts[[key]] %||% contexts[[id]]
        if (is.null(context)) next
        status <- as.character(jobs$status[[row]])
        if (!identical(as.character(context$status %||% ""), status)) {
          context$status <- status
          changed <- TRUE
        }
        if (materialization_due(context) &&
            nzchar(as.character(context$project %||% "")) &&
            nzchar(as.character(context$version %||% ""))) {
          context$materialize_attempted_at <- format(
            Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC"
          )
          materialized <- tryCatch(
            materialize_job_result(context, q, id), error = identity
          )
          if (inherits(materialized, "error")) {
            message <- conditionMessage(materialized)
            if (!identical(as.character(context$materialize_error %||% ""), message)) {
              append_log(
                paste("Completed job could not yet be saved as a model run:", message),
                "error"
              )
            }
            context$materialize_error <- message
            changed <- TRUE
          } else {
            context <- materialized
            changed <- TRUE
            append_log(paste("Saved completed queued run", context$label), "info")
          }
        }
        contexts[[key]] <- context
        if (!identical(key, id)) contexts[[id]] <- NULL
      }
      if (changed) {
        state$job_context <- contexts
        save_client_settings()
      }
      invisible(changed)
    }
    refresh_jobs <- function(start = FALSE, background = FALSE) {
      q <- active_queue()
      jobs <- if (is.null(q)) {
        data.frame()
      } else if (inherits(q, "LibeRQueue")) {
        q$poll(start = isTRUE(start))
        q$list()
      } else {
        with_ui_remote_timeout(q, background, q$list())
      }
      jobs <- as.data.frame(jobs %||% data.frame(), stringsAsFactors = FALSE)
      if (nrow(jobs) && length(state$hidden_jobs) && "id" %in% names(jobs)) {
        jobs <- jobs[!jobs$id %in% state$hidden_jobs, , drop = FALSE]
      }
      state$jobs <- jobs
      reconcile_job_results(jobs, q)
      if (!isTRUE(background) && !is.null(state$selected_job) &&
          nzchar(state$selected_job)) {
        selected <- as.character(state$selected_job)
        state$job_log <- tryCatch(
          with_ui_remote_timeout(
            q, FALSE,
            c("--- stdout ---", q$logs(selected, stream = "stdout"),
              "--- stderr ---", q$logs(selected, stream = "stderr"))
          ),
          error = function(error) c("Unable to read worker log:", conditionMessage(error))
        )
      }
      state$refreshed <- paste("Last refreshed", format(Sys.time(), "%H:%M:%S"))
      invisible(jobs)
    }
    load_snapshot <- function(project_id, snapshot_id = "latest") {
      opened <- nm_project_load(workspace, project_id, snapshot = snapshot_id)
      state$model <- opened$model
      state$draft_outputs <- NULL
      state$data <- opened$data
      state$project <- opened$project
      state$report_design <- saved_report_design(opened$project)
      state$snapshot <- opened$id
      state$run <- NULL
      state$diagnostics <- list()
      state$result <- structure(list(), class = "liber_gui_validation")
      reset_lazy_payloads()
      invalidate_ai_context()
      update_fit(NULL)
      opened
    }
    load_run <- function(project_id, run_id) {
      opened <- nm_project_load(workspace, project_id, snapshot = run_id)
      records <- nm_project_list(workspace, project_id)
      metadata <- records[records$id == run_id & records$entry_type == "run", , drop = FALSE]
      if (!nrow(metadata)) .nm_stop("The selected item is not a model run.")
      opened$entry_type <- "run"
      opened$parent_id <- opened$parent_id %||% metadata$parent_id[[1L]]
      opened$run_number <- opened$run_number %||% metadata$run_number[[1L]]
      state$model <- opened$model
      state$draft_outputs <- NULL
      state$data <- opened$data
      state$project <- opened$project
      state$report_design <- saved_report_design(opened$project)
      state$snapshot <- opened$parent_id
      state$run <- opened$id
      state$diagnostics <- nm_project_load_diagnostics(workspace, project_id, opened$id)
      state$result <- opened$result
      reset_lazy_payloads()
      update_fit(if (inherits(opened$result, "nm_fit")) opened$result else NULL)
      opened
    }
    persist_run <- function(result, label, project_id = state$project,
                            run_model = state$model, run_data = state$data,
                            parent_version = state$snapshot) {
      is_simulation <- is.data.frame(result) && "IPRED" %in% names(result)
      if (is.null(workspace) || is.null(project_id) ||
          (!inherits(result, "nm_fit") && !is_simulation)) return(NULL)
      if (is.null(parent_version)) {
        parent_version <- nm_project_save(
          workspace, project_id, run_model, run_data, NULL, label = NULL
        )
      }
      id <- nm_project_save_run(
        workspace, project_id, parent_version, result, label = label,
        model = run_model, data = run_data
      )
      invalidate_ai_context()
      if (identical(project_id, state$project)) {
        state$snapshot <- parent_version
        state$run <- id
        state$diagnostics <- list()
        reset_lazy_payloads(data = FALSE)
        if (is_simulation) {
          state$result <- result
          update_fit(NULL)
        }
      }
      id
    }

    shiny::observe({
      q <- active_queue()
      shiny::invalidateLater(if (inherits(q, "LibeRQueue")) 1000 else 3000, session)
      contexts <- state$job_context
      needs_reconciliation <- any(vapply(contexts, function(context) {
        status <- as.character(context$status %||% "queued")
        status %in% c("queued", "running")
      }, logical(1)))
      should_poll <- identical(state$active_page, "home") &&
        !isTRUE(state$comparison_open) && needs_reconciliation &&
        isTRUE(poll_backoff$ready) &&
        as.numeric(Sys.time()) >= poll_backoff$until
      if (should_poll && !is.null(q)) {
        tryCatch({
          refresh_jobs(
            start = inherits(q, "LibeRQueue"), background = !inherits(q, "LibeRQueue")
          )
          poll_backoff$until <- 0
        }, error = function(error) {
          poll_backoff$until <- as.numeric(Sys.time()) + 30
          append_log(paste("Queue polling paused for 30 seconds:", conditionMessage(error)),
                     "error")
        })
      }
    })

    session$onFlushed(function() {
      shiny::isolate({
        q <- active_queue()
        if (inherits(q, "LibeRQueue")) {
          tryCatch(refresh_jobs(start = TRUE), error = function(error) {
            append_log(paste("Initial queue refresh failed:", conditionMessage(error)), "error")
          })
        } else if (inherits(q, "LibeRRemote")) {
          state$refreshed <- "Remote queue restored; open Jobs to connect"
        }
        poll_backoff$ready <- TRUE
      })
    }, once = TRUE)

    output$workbench <- renderLiberWorkbench({
      remote_ids <- names(state$remote_meta)
      queues <- c(
        list(list(
          id = "local",
          name = if (is.null(queue)) "Local (in process)" else "Local queue",
          url = ""
        )),
        unname(lapply(remote_ids, function(id) c(list(id = id), state$remote_meta[[id]])))
      )
      q <- active_queue()
      selected_meta <- if (identical(state$queue_id, "local")) NULL else {
        state$remote_meta[[state$queue_id]]
      }
      server_info <- list(
        mode = if (identical(state$queue_id, "local")) {
          if (is.null(queue)) "local" else "local queue"
        } else "remote queue",
        connected = TRUE, platform = R.version$platform,
        worker = if (is.null(q)) "in-process" else if (inherits(q, "LibeRQueue")) {
          paste(q$max_workers, "callr worker(s)")
        } else paste0("remote user ", selected_meta$user %||% "authenticated"),
        isolation = if (is.null(q)) "current R session" else if (inherits(q, "LibeRQueue")) {
          paste("tenant", q$user)
        } else "server-managed tenant",
        queue_id = state$queue_id, queues = queues, refreshed = state$refreshed,
        queue_root = if (inherits(q, "LibeRQueue")) q$root else "",
        package_version = as.character(utils::packageVersion("LibeRation")),
        icon = favicon_href,
        job_count = nrow(as.data.frame(state$jobs %||% data.frame()))
      )
      run_output <- NULL
      if (isTRUE(state$data_payload)) {
        if (inherits(state$fit, "nm_fit")) {
          run_output <- state$fit$output
        } else if (is.data.frame(state$result) && length(state$model$OUTPUT %||% character())) {
          selected <- intersect(state$model$OUTPUT, names(state$result))
          if (length(selected)) {
            run_output <- state$result[, selected, drop = FALSE]
            run_output$.ROW <- seq_len(nrow(run_output))
          }
        }
      }
      liber_workbench(
        model = state$model, data = state$data, jobs = state$jobs,
        result = state$result, fit = state$fit_payload, report = state$report,
        report_design = state$report_design,
        diagnostics = state$diagnostics, hmm = state$hmm_payload,
        kalman = state$kalman_payload,
        log = list(level = state$log_level, current = state$log_current,
                   history = state$log_history),
        job_log = state$job_log, server = server_info,
        workspace = .liber_gui_workspace(
          workspace, state$project, state$snapshot, state$run,
          pending_jobs = state$job_context, jobs = state$jobs,
          queue_id = state$queue_id
        ),
        library = .liber_gui_library(),
        ai = c(state$ai_config, list(
          worker_url = ai_worker_href, models = ai_models,
          secure_context = TRUE,
          privacy = "Inference runs in a dedicated browser worker with no tools. Network APIs are disabled after the selected model has loaded."
        )),
        ai_context = state$ai_context,
        report_ai_context = state$report_ai_context,
        report_directory = state$report_design$style$output_directory %||%
          .liber_report_default_directory(workspace, state$project),
        output_catalog = state$draft_outputs,
        run_output = run_output,
        data_payload = state$data_payload,
        gof_payload = state$gof_payload,
        diagnostic_payload = state$diagnostic_payload,
        input_id = "liber_workbench", height = "100vh"
      )
    })

    shiny::observeEvent(input$liber_workbench_event, {
      event <- input$liber_workbench_event
      action <- as.character(event$action %||% "")

      if (identical(action, "ai_settings")) {
        allowed_models <- vapply(ai_models, `[[`, character(1), "id")
        requested_help <- as.character(event$help_model %||% event$model %||%
          state$ai_config$help_model %||% state$ai_config$model)[[1L]]
        requested_report <- as.character(event$report_model %||%
          state$ai_config$report_model)[[1L]]
        requested_help_context <- .liber_ai_context_setting(
          event$help_context %||% state$ai_config$help_context
        )
        requested_report_context <- .liber_ai_context_setting(
          event$report_context %||% state$ai_config$report_context
        )
        if (!requested_help %in% allowed_models) {
          requested_help <- .liber_ai_default_help_model()
        }
        if (!requested_report %in% c("same_as_help", allowed_models)) {
          requested_report <- .liber_ai_default_report_model()
        }
        state$ai_config <- list(
          activated = isTRUE(event$activated),
          consented = isTRUE(event$consented),
          help_model = requested_help,
          report_model = requested_report,
          help_context = requested_help_context,
          report_context = requested_report_context,
          model = requested_help
        )
        save_client_settings()
        append_log(if (state$ai_config$activated) {
          paste0(
            "Browser-local AI enabled; Help and Report models remain unloaded ",
            "until first use"
          )
        } else {
          "Browser-local AI settings saved; AI remains disabled"
        }, "info")
        return(invisible(NULL))
      }
      if (identical(action, "ai_context_request")) {
        requested_project <- as.character(event$project %||% "")[[1L]]
        request_id <- as.character(event$requestId %||% "")[[1L]]
        requested_scope <- as.character(event$scope %||% "index")[[1L]]
        if (!requested_scope %in% c("index", "results")) requested_scope <- "index"
        context <- tryCatch({
          if (!nzchar(requested_project) ||
              !identical(requested_project, as.character(state$project %||% ""))) {
            .nm_stop("The requested Help context is no longer the selected project.")
          }
          .liber_gui_ai_context(
            workspace, requested_project, selected_run = state$run,
            max_runs = if (identical(requested_scope, "results")) 12L else 30L,
            detail = requested_scope
          )
        }, error = function(error) list(
          available = FALSE, project = requested_project, project_name = "",
          scope = requested_scope, message = conditionMessage(error), run_count = 0L,
          included_runs = 0L, omitted_runs = 0L, runs = list()
        ))
        context$request_id <- request_id
        state$ai_context <- context
        append_log(if (isTRUE(context$available)) {
          paste("Loaded", context$included_runs, "saved run summaries for local Help AI")
        } else {
          paste("Help AI project summaries unavailable:", context$message)
        }, if (isTRUE(context$available)) "info" else "error")
        return(invisible(NULL))
      }
      if (identical(action, "diagram_preview")) {
        return(record("Diagram code preview generated", {
          diagram <- diagram_from_event(event$graph)
          preview <- nm_diagram_preview(diagram)
          state$diagram_candidate <- diagram
          baseline <- if (is.null(state$model)) "" else {
            as.character(state$model$GRAPH$generated$DES %||% "")
          }
          structure(list(
            nonce = as.character(event$nonce %||% format(Sys.time(), "%OS6")),
            preview = list(
              pred = preview$PRED, des = preview$DES, error = preview$ERROR,
              theta = .liber_gui_rows(preview$THETAS),
              omega = .liber_gui_rows(preview$OMEGAS),
              sigma = .liber_gui_rows(preview$SIGMAS)
            ),
            code_changed = nzchar(baseline) && !identical(trimws(baseline), trimws(state$model$DES))
          ), class = "liber_gui_diagram_preview")
        }))
      }
      if (identical(action, "diagram_apply")) {
        return(record("Visual model code generated and applied", {
          diagram <- if (is.list(event$graph)) diagram_from_event(event$graph) else state$diagram_candidate
          if (!inherits(diagram, "nm_model_diagram")) {
            .nm_stop("Preview the diagram before applying generated code.")
          }
          generated <- nm_diagram_generate(
            diagram,
            input = state$model$INPUT %||% c("ID", "TIME", "EVID", "AMT", "CMT", "DV", "MDV"),
            output = state$model$OUTPUT %||% c("PRED", "IPRED", "CWRES")
          )
          nm_compile(generated)
          state$model <- generated
          state$diagram_candidate <- NULL
          state$draft_outputs <- NULL
          state$snapshot <- NULL
          state$run <- NULL
          state$diagnostics <- list()
          reset_lazy_payloads(data = FALSE)
          update_fit(NULL)
          structure(list(), class = "liber_gui_validation")
        }))
      }
      if (identical(action, "report_document")) {
        value <- tryCatch(report_uploaded_text(event$name, event$data), error = identity)
        if (inherits(value, "error")) {
          append_log(conditionMessage(value), "error")
        } else {
          session$sendCustomMessage("liber-report-document", list(
            input_id = "liber_workbench", nonce = as.character(event$nonce %||% ""),
            block_id = as.character(event$blockId),
            name = basename(as.character(event$name)), text = value
          ))
          append_log(paste("Loaded report source", basename(as.character(event$name))), "info")
        }
        return(invisible(NULL))
      }
      if (identical(action, "report_directory_choose")) {
        selected <- tryCatch(
          .liber_report_choose_directory(
            as.character(event$directory %||%
              .liber_report_default_directory(workspace, state$project))
          ),
          error = identity
        )
        if (inherits(selected, "error")) {
          append_log(conditionMessage(selected), "error")
        } else if (!is.null(selected)) {
          session$sendCustomMessage("liber-report-directory", list(path = selected))
          append_log(paste("Report output folder selected:", selected), "info")
        }
        return(invisible(NULL))
      }
      if (identical(action, "report_ai_context_request")) {
        requested_project <- as.character(event$project %||% state$project %||% "")[[1L]]
        request_id <- as.character(event$requestId %||% "")[[1L]]
        run_ids <- unique(as.character(unlist(event$runs %||% character())))
        context <- tryCatch({
          if (!nzchar(requested_project) ||
              !identical(requested_project, as.character(state$project %||% ""))) {
            .nm_stop("The report evidence request is no longer for the selected project.")
          }
          .liber_gui_report_ai_context(workspace, requested_project, run_ids)
        }, error = function(error) list(
          available = FALSE, project = requested_project, project_name = "",
          run_ids = run_ids, message = conditionMessage(error), runs = list()
        ))
        context$request_id <- request_id
        state$report_ai_context <- context
        append_log(if (isTRUE(context$available)) {
          paste("Loaded", length(context$runs), "selected run(s) for local Report AI")
        } else {
          paste("Report AI evidence unavailable:", context$message)
        }, if (isTRUE(context$available)) "info" else "error")
        return(invisible(NULL))
      }
      if (identical(action, "report_design_save")) {
        return(record("Report workflow saved", {
          if (is.null(state$project)) .nm_stop("Open a project before saving a report workflow.")
          design <- report_design_from_event(event)
          nm_report_design_save(workspace, state$project, design)
          state$report_design <- design
          structure(list(), class = "liber_gui_validation")
        }))
      }
      if (identical(action, "report_design_render")) {
        return(record("DOCX/PDF report generated", {
          if (is.null(state$project)) .nm_stop("Open a project before generating a report.")
          design <- report_design_from_event(event)
          nm_report_design_save(workspace, state$project, design)
          state$report_design <- design
          name <- as.character(event$name %||% "liberation-report")
          directory <- trimws(as.character(event$directory %||% "")[[1L]])
          if (!nzchar(directory)) {
            directory <- .liber_report_default_directory(workspace, state$project)
          }
          directory <- path.expand(directory)
          nm_report_design_render(
            design, workspace, state$project, directory = directory,
            name = name, formats = design$formats
          )
        }))
      }

      if (identical(action, "load_payload")) {
        kind <- tolower(as.character(event$kind %||% ""))
        if (identical(kind, "data")) {
          state$data_payload <- !is.null(state$data)
        } else if (identical(kind, "gof")) {
          if (inherits(state$fit, "nm_fit")) {
            gof <- state$diagnostics$gof %||% nm_gof(state$fit)
            if (is.null(state$diagnostics$gof) && !is.null(state$run)) {
              state$diagnostics <- nm_project_save_diagnostics(
                workspace, state$project, state$run, list(gof = gof)
              )
              invalidate_ai_context()
            }
            state$gof_payload <- TRUE
            state$fit_payload <- .liber_gui_fit(state$fit, include_gof = TRUE, gof = gof)
          }
        } else if (identical(kind, "hmm")) {
          if (!inherits(state$fit, "nm_fit") || is.null(state$model$HMM_CONFIG)) {
            append_log("Open a fitted hidden Markov model before loading HMM results.", "error")
            return(invisible(NULL))
          }
          decoded <- tryCatch(
            nm_hmm_decode(state$fit, method = "all"),
            error = function(error) {
              append_log(paste("HMM decoding failed:", conditionMessage(error)), "error")
              NULL
            }
          )
          if (is.null(decoded)) return(invisible(NULL))
          state$hmm_payload <- .liber_gui_hmm(decoded, available = TRUE)
        } else if (identical(kind, "kalman")) {
          if (!inherits(state$fit, "nm_fit") || is.null(state$model$KALMAN_CONFIG)) {
            append_log("Open a fitted linear state-space model before loading state estimates.", "error")
            return(invisible(NULL))
          }
          decoded <- tryCatch(
            nm_kalman_decode(state$fit, type = "individual"),
            error = function(error) {
              append_log(paste("State-space decoding failed:", conditionMessage(error)), "error")
              NULL
            }
          )
          if (is.null(decoded)) return(invisible(NULL))
          state$kalman_payload <- .liber_gui_kalman(decoded, available = TRUE)
        } else if (kind %in% c("vpc", "npc", "npde", "vpc_categorical", "vpc_count",
                              "vpc_tte", "vpc_competing", "vpc_recurrent",
                              "bootstrap", "profile", "scm") &&
                   !is.null(state$diagnostics[[kind]])) {
          state$diagnostic_payload <- unique(c(state$diagnostic_payload, kind))
        }
        append_log(paste("Loaded", toupper(kind), "view data"), "info")
        return(invisible(NULL))
      }
      if (identical(action, "page_change")) {
        state$active_page <- as.character(event$page %||% "home")
        if (identical(state$active_page, "jobs")) {
          tryCatch(refresh_jobs(start = TRUE), error = function(error) {
            append_log(paste("Queue refresh failed:", conditionMessage(error)), "error")
          })
        }
        return(invisible(NULL))
      }
      if (identical(action, "comparison_close")) {
        state$comparison_open <- FALSE
        if (inherits(state$result, "liber_gui_comparison")) state$result <- NULL
        return(invisible(NULL))
      }

      if (identical(action, "clear_log")) {
        state$log_history <- character()
        state$log_current <- "Log cleared"
        state$log_level <- "info"
        return(invisible(NULL))
      }
      if (identical(action, "validate")) {
        return(record("Model validation completed", {
          if (is.null(state$model)) .nm_stop("Load a model before validation.")
          draft <- .liber_model_from_event(state$model, event)
          nm_compile(draft)
          state$draft_outputs <- nm_model_outputs(draft)
          structure(
            list(
              kind = "model_validation",
              nonce = format(Sys.time(), "%Y%m%d%H%M%OS6"),
              outputs = state$draft_outputs,
              parameters = list(
                theta = .liber_gui_rows(draft$THETAS),
                omega = .liber_gui_rows(draft$OMEGAS),
                sigma = .liber_gui_rows(draft$SIGMAS),
                priors = .liber_gui_rows(draft$LIK_CONFIG$priors %||% data.frame())
              )
            ),
            class = "liber_gui_validation"
          )
        }))
      }
      if (identical(action, "load_csv")) {
        return(record("Dataset imported", {
          connection <- textConnection(as.character(event$text %||% ""))
          on.exit(close(connection), add = TRUE)
          dataset <- nm_dataset(utils::read.csv(connection, check.names = FALSE))
          attr(dataset, "name") <- as.character(event$name %||% "Imported dataset")
          state$data <- dataset
          state$data_payload <- FALSE
          structure(list(), class = "liber_gui_validation")
        }))
      }
      if (identical(action, "control_import")) {
        return(record("NONMEM control stream imported", {
          text <- as.character(event$text %||% "")
          if (!nzchar(text)) .nm_stop("Select a NONMEM control stream.")
          control <- nm_control_read(text, strict = TRUE)
          project_id <- state$project
          if (is.null(project_id) || isTRUE(event$newProject)) {
            created <- nm_project_create(
              workspace,
              as.character(event$projectName %||% control$problem %||% "Imported NONMEM project")
            )
            project_id <- created$id
          }
          imported_data <- state$data
          data_text <- as.character(event$dataText %||% "")
          if (nzchar(data_text)) {
            connection <- textConnection(data_text)
            on.exit(close(connection), add = TRUE)
            filename <- tolower(as.character(event$dataName %||% "data.csv"))
            frame <- if (grepl("[.]csv$", filename)) {
              utils::read.csv(connection, check.names = FALSE)
            } else utils::read.table(connection, header = TRUE, check.names = FALSE)
            imported_data <- nm_dataset(frame)
            attr(imported_data, "name") <- basename(filename)
          }
          state$project <- project_id
          state$report_design <- saved_report_design(project_id)
          state$model <- control$model
          state$draft_outputs <- NULL
          state$data <- imported_data
          state$snapshot <- nm_project_save(
            workspace, project_id, control$model, imported_data, NULL,
            label = as.character(event$label %||% control$problem %||% "NONMEM import")
          )
          state$run <- NULL
          state$diagnostics <- list()
          reset_lazy_payloads()
          update_fit(NULL)
          if (length(control$compatibility$warnings)) {
            append_log(paste(control$compatibility$warnings, collapse = " | "), "warning")
          }
          structure(list(), class = "liber_gui_validation")
        }))
      }
      if (identical(action, "library_import")) {
        return(record("LibeRary model imported", {
          if (!requireNamespace("LibeRary", quietly = TRUE)) {
            .nm_stop("Install LibeRary to import catalogue models.")
          }
          library_id <- as.character(event$libraryId %||% "")
          if (!nzchar(library_id)) .nm_stop("Select a LibeRary model.")
          target <- state$project
          if (isTRUE(event$newProject) || is.null(target)) {
            requested_name <- trimws(as.character(event$projectName %||% ""))
            if (nzchar(requested_name)) {
              target <- nm_project_create(workspace, requested_name)$id
            } else target <- paste0("lib_", library_id)
          }
          imported <- LibeRary::library_use_in_workspace(
            library_id, project = target, workspace = workspace,
            version_label = as.character(event$label %||% "")
          )
          opened <- nm_project_load(workspace, imported$project, imported$version_id)
          state$project <- imported$project
          state$report_design <- saved_report_design(imported$project)
          state$snapshot <- imported$version_id
          state$run <- NULL
          state$model <- opened$model
          state$draft_outputs <- NULL
          state$data <- opened$data
          state$diagnostics <- list()
          reset_lazy_payloads()
          update_fit(NULL)
          structure(list(), class = "liber_gui_validation")
        }))
      }
      if (identical(action, "control_export")) {
        return(record("NONMEM control stream exported", {
          if (is.null(state$model)) .nm_stop("Load a model before exporting a control stream.")
          directory <- file.path(workspace$path, "exports")
          if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE, showWarnings = FALSE)) {
            .nm_stop("Unable to create the workspace exports directory.")
          }
          name <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(event$name %||% "model.ctl"))
          if (!grepl("[.]ctl$", name, ignore.case = TRUE)) name <- paste0(name, ".ctl")
          path <- file.path(directory, basename(name))
          nm_control_write(state$model, path, data = as.character(event$dataPath %||% "data.csv"))
          structure(list(path = normalizePath(path, winslash = "/", mustWork = TRUE)),
                    class = "liber_gui_control_export")
        }))
      }
      if (identical(action, "model_template")) {
        return(record("Model template created", {
          structural <- as.character(event$structuralTemplate %||% "standard")
          state$model <- if (!identical(structural, "standard")) {
            nm_model_template(structural)
          } else {
            .liber_model_template(
              as.integer(event$advan %||% 1L), trans = event$trans,
              n_state = event$nState,
              problem = as.character(event$problem %||% "Template model")
            )
          }
          if (!identical(structural, "standard")) {
            requested_name <- trimws(as.character(event$problem %||% ""))
            if (nzchar(requested_name)) attr(state$model, "name") <- requested_name
          }
          state$draft_outputs <- NULL
          state$snapshot <- if (!is.null(workspace) && !is.null(state$project)) {
            nm_project_save(
              workspace, state$project, state$model, state$data, NULL,
              label = as.character(event$label %||% attr(state$model, "name"))
            )
          } else NULL
          state$run <- NULL
          state$diagnostics <- list()
          reset_lazy_payloads(data = FALSE)
          update_fit(NULL)
          structure(list(), class = "liber_gui_validation")
        }))
      }
      if (identical(action, "update_model") || identical(action, "update_parameters")) {
        return(record("Model changes applied", {
          if (is.null(state$model)) .nm_stop("Load a model before editing it.")
          edited <- .liber_model_from_event(state$model, event)
          nm_compile(edited)
          state$model <- edited
          state$draft_outputs <- NULL
          state$snapshot <- NULL
          state$run <- NULL
          state$diagnostics <- list()
          reset_lazy_payloads(data = FALSE)
          structure(list(), class = "liber_gui_validation")
        }))
      }
      if (identical(action, "simulate")) {
        return(record("Simulation completed or submitted", {
          if (is.null(state$model) || is.null(state$data)) {
            .nm_stop("Load both a model and dataset before simulation.")
          }
          simulation_data <- if (is.null(event$nSubjects) && !isTRUE(event$useDesign)) {
            state$data
          } else {
            .liber_simulation_dataset(state$data, state$model, event)
          }
          arguments <- list(
            nsim = max(1L, as.integer(event$replicates %||% 1L)),
            random_effects = TRUE, residual = TRUE,
            seed = as.integer(event$seed %||% sample.int(.Machine$integer.max, 1L)),
            n_cores = max(1L, as.integer(event$nCores %||% 1L))
          )
          if (isTRUE(event$useFit) && inherits(state$fit, "nm_fit")) {
            arguments$theta <- state$fit$theta
            arguments$sigma <- state$fit$sigma
            arguments$omega <- state$fit$omega
          }
          q <- active_queue()
          if (is.null(q)) {
            result <- do.call(nm_simulate, c(list(model = state$model, data = simulation_data), arguments))
            persist_run(result, as.character(event$label %||% "Simulation"),
                        run_data = simulation_data)
            result
          } else {
            if (!requireNamespace("LibeRties", quietly = TRUE)) {
              .nm_stop("LibeRties is required for queued execution.")
            }
            parent_version <- ensure_parent_version()
            id <- q$submit(LibeRties::ls_job(
              "simulate", state$model, simulation_data, arguments = arguments,
              label = as.character(event$label %||% "Simulation")
            ))
            contexts <- state$job_context
            contexts[[job_context_key(id)]] <- list(
              project = state$project, model = state$model, data = simulation_data,
              label = as.character(event$label %||% "Simulation"),
              version = parent_version, type = "simulate", method = "simulation",
              queue_id = state$queue_id, job_id = as.character(id), status = "queued"
            )
            state$job_context <- contexts
            save_client_settings()
            refresh_jobs(start = TRUE)
            structure(list(id = id), class = "liber_gui_queued")
          }
        }))
      }
      if (identical(action, "estimate")) {
        return(record("Estimation completed or submitted", {
          if (is.null(state$model) || is.null(state$data)) {
            .nm_stop("Load both a model and dataset before estimation.")
          }
          stages <- .liber_estimation_stages(event)
          sequential <- length(stages) > 1L
          methods <- vapply(stages, `[[`, character(1), "method")
          arguments <- if (sequential) {
            list(stages = stages)
          } else {
            c(list(method = stages[[1L]]$method), stages[[1L]]$arguments)
          }
          q <- active_queue()
          run_label <- trimws(as.character(event$label %||% ""))
          method_label <- paste(methods, collapse = " -> ")
          if (!nzchar(run_label)) run_label <- paste(method_label, "estimation")
          if (is.null(q)) {
            result <- do.call(
              if (sequential) nm_est_sequence else nm_est,
              c(list(model = state$model, data = state$data), arguments)
            )
            persist_run(result, run_label)
            result
          } else {
            if (!requireNamespace("LibeRties", quietly = TRUE)) {
              .nm_stop("LibeRties is required for queued execution.")
            }
            parent_version <- ensure_parent_version()
            id <- q$submit(LibeRties::ls_job(
              if (sequential) "estimate_sequence" else "estimate",
              state$model, state$data, arguments = arguments,
              label = run_label
            ))
            contexts <- state$job_context
            contexts[[job_context_key(id)]] <- list(
              project = state$project, model = state$model, data = state$data,
              label = run_label, version = parent_version,
              type = "estimate", method = method_label,
              queue_id = state$queue_id, job_id = as.character(id), status = "queued"
            )
            state$job_context <- contexts
            save_client_settings()
            refresh_jobs(start = TRUE)
            structure(list(id = id), class = "liber_gui_queued")
          }
        }))
      }
      if (identical(action, "run_diagnostic")) {
        return(record("Selected diagnostics completed and saved", {
          if (!inherits(state$fit, "nm_fit") || is.null(state$run)) {
            .nm_stop("Open a saved estimation run before running diagnostics.")
          }
          selected <- unique(tolower(as.character(unlist(event$types %||% character()))))
          selected <- intersect(selected, c(
            "vpc", "npc", "npde", "vpc_categorical", "vpc_count",
            "vpc_tte", "vpc_competing", "vpc_recurrent"
          ))
          if (!length(selected)) .nm_stop("Select at least one diagnostic.")
          nsim <- max(20L, as.integer(event$nsim %||% 200L))
          seed <- as.integer(event$seed %||% 20260713L)
          created <- list()
          if ("vpc" %in% selected) {
            created$vpc <- nm_vpc(
              state$fit, nsim = nsim, seed = seed,
              pc_correct = isTRUE(event$pcCorrect),
              stratify = event$stratify %||% NULL
            )
          }
          if ("npc" %in% selected) created$npc <- nm_npc(state$fit, nsim = nsim, seed = seed)
          if ("npde" %in% selected) created$npde <- nm_npde(state$fit, nsim = nsim, seed = seed)
          if ("vpc_categorical" %in% selected) {
            created$vpc_categorical <- nm_vpc_categorical(
              state$fit, outcome = as.character(event$categoricalOutcome %||% "DV"),
              nsim = nsim, seed = seed
            )
          }
          if ("vpc_count" %in% selected) {
            created$vpc_count <- nm_vpc_count(
              state$fit, outcome = as.character(event$countOutcome %||% "DV"),
              dvid = suppressWarnings(as.numeric(event$countDvid %||% NA_real_)),
              nsim = nsim, seed = seed
            )
          }
          if ("vpc_tte" %in% selected) {
            created$vpc_tte <- nm_vpc_tte(
              state$fit, event = as.character(event$tteEvent %||% "DV"),
              nsim = nsim, seed = seed
            )
          }
          if ("vpc_competing" %in% selected) {
            created$vpc_competing <- nm_vpc_competing(
              state$fit, event = as.character(event$tteEvent %||% "DV"),
              dvid = suppressWarnings(as.numeric(event$competingDvid %||% NA_real_)),
              nsim = nsim, seed = seed
            )
          }
          if ("vpc_recurrent" %in% selected) {
            created$vpc_recurrent <- nm_vpc_recurrent(
              state$fit, event = as.character(event$tteEvent %||% "DV"),
              dvid = suppressWarnings(as.numeric(event$recurrentDvid %||% NA_real_)),
              nsim = nsim, seed = seed
            )
          }
          state$diagnostics <- nm_project_save_diagnostics(
            workspace, state$project, state$run, created
          )
          invalidate_ai_context()
          state$diagnostic_payload <- character()
          state$diagnostics
        }, update_result = FALSE))
      }
      if (identical(action, "run_uncertainty")) {
        return(record("Uncertainty analyses completed and saved", {
          if (!inherits(state$fit, "nm_fit") || is.null(state$run)) {
            .nm_stop("Open a saved estimation run before running uncertainty analyses.")
          }
          selected <- unique(tolower(as.character(unlist(event$types %||% character()))))
          selected <- intersect(selected, c("bootstrap", "profile"))
          if (!length(selected)) .nm_stop("Select bootstrap or profile likelihood.")
          created <- list()
          if ("bootstrap" %in% selected) {
            created$bootstrap <- nm_bootstrap(
              state$fit, n = as.integer(event$replicates %||% 100L),
              seed = as.integer(event$seed %||% 20260713L),
              level = as.numeric(event$level %||% 0.95),
              maxit = as.integer(event$maxit %||% 100L)
            )
          }
          if ("profile" %in% selected) {
            parameters <- trimws(unlist(strsplit(as.character(event$parameters %||% ""), "[,;[:space:]]+")))
            parameters <- parameters[nzchar(parameters)]
            created$profile <- nm_profile(
              state$fit, parameters = if (length(parameters)) parameters else NULL,
              points = as.integer(event$points %||% 9L),
              span = as.numeric(event$span %||% 3),
              level = as.numeric(event$level %||% 0.95),
              maxit = as.integer(event$maxit %||% 100L)
            )
          }
          state$diagnostics <- nm_project_save_diagnostics(
            workspace, state$project, state$run, created
          )
          invalidate_ai_context()
          state$diagnostic_payload <- character()
          state$diagnostics
        }, update_result = FALSE))
      }
      if (identical(action, "run_scm")) {
        return(record("Stepwise covariate modelling completed", {
          if (!inherits(state$fit, "nm_fit") || is.null(state$project)) {
            .nm_stop("Open an estimation run before starting SCM.")
          }
          lines <- strsplit(as.character(event$candidates %||% ""), "\r?\n", perl = TRUE)[[1L]]
          lines <- trimws(lines[nzchar(trimws(lines))])
          if (!length(lines)) .nm_stop("Enter at least one parameter,covariate candidate.")
          fields <- lapply(lines, function(line) trimws(strsplit(line, ",", fixed = TRUE)[[1L]]))
          item <- function(field, index, default = "") {
            if (length(field) >= index && nzchar(field[[index]])) field[[index]] else default
          }
          candidates <- do.call(rbind, lapply(fields, function(field) data.frame(
            parameter = item(field, 1L), covariate = item(field, 2L),
            form = item(field, 3L, "continuous"),
            reference = item(field, 4L, NA_character_),
            category = item(field, 5L, NA_character_), stringsAsFactors = FALSE
          )))
          scm <- nm_scm(
            state$fit, candidates,
            direction = as.character(event$direction %||% "both"),
            p_forward = as.numeric(event$pForward %||% 0.05),
            p_backward = as.numeric(event$pBackward %||% 0.01),
            max_steps = as.integer(event$maxSteps %||% 20L),
            maxit = as.integer(event$maxit %||% 100L)
          )
          version <- nm_project_save(
            workspace, state$project, scm$final_model, state$data, NULL,
            label = as.character(event$label %||% "SCM model")
          )
          run <- nm_project_save_run(
            workspace, state$project, version, scm$final_fit,
            label = paste("SCM", scm$final_fit$method), model = scm$final_model,
            data = state$data
          )
          state$diagnostics <- nm_project_save_diagnostics(
            workspace, state$project, run, list(scm = scm)
          )
          invalidate_ai_context()
          state$model <- scm$final_model
          state$draft_outputs <- NULL
          state$snapshot <- version
          state$run <- run
          reset_lazy_payloads(data = FALSE)
          update_fit(scm$final_fit)
          scm
        }, update_result = FALSE))
      }
      if (identical(action, "report")) {
        return(record("Report generated", {
          if (!inherits(state$fit, "nm_fit")) .nm_stop("Run or retrieve an estimation before reporting.")
          name <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(event$name %||% "report"))
          name <- sub("[.]pdf$", "", basename(name), ignore.case = TRUE)
          if (!nzchar(name)) name <- "report"
          directory <- if (is.null(workspace)) tempdir() else file.path(workspace$path, "reports")
          nm_report(
            state$fit, file.path(directory, paste0(name, ".pdf")),
            sections = as.character(unlist(event$sections %||% c(
              "summary", "parameters", "gof", "eta", "narrative_stub"
            ))), vpc = state$diagnostics$vpc
          )
        }))
      }
      if (identical(action, "project_create")) {
        return(record("Project created", {
          if (is.null(workspace)) .nm_stop("Launch with `workspace` to manage projects.")
          created <- nm_project_create(
            workspace, as.character(event$name %||% "New project"),
            description = as.character(event$description %||% "")
          )
          initialize <- identical(as.character(event$mode %||% "empty"), "template")
          initial <- tryCatch({
            if (!initialize) {
              NULL
            } else {
              model <- .liber_model_template(
                as.integer(event$advan %||% 4L), trans = event$trans,
                n_state = event$nState,
                problem = as.character(event$problem %||% "Template model")
              )
              source <- as.character(event$dataSource %||% "synthetic")
              data <- if (identical(source, "upload")) {
                text <- as.character(event$text %||% "")
                if (!nzchar(text)) .nm_stop("Select a dataset file for the initial model version.")
                connection <- textConnection(text)
                on.exit(close(connection), add = TRUE)
                filename <- tolower(as.character(event$fileName %||% "dataset.csv"))
                frame <- if (grepl("[.]csv$", filename)) {
                  utils::read.csv(connection, check.names = FALSE)
                } else {
                  utils::read.table(connection, header = TRUE, check.names = FALSE)
                }
                dataset <- nm_dataset(frame)
                attr(dataset, "name") <- basename(as.character(event$fileName %||% "Uploaded dataset"))
                dataset
              } else {
                .liber_builtin_dataset(
                  model, as.character(event$example %||% "theophylline"),
                  as.integer(event$nSubjects %||% 10L)
                )
              }
              snapshot <- nm_project_save(
                workspace, created$id, model, data, NULL,
                label = as.character(event$label %||% attr(model, "name"))
              )
              list(model = model, data = data, snapshot = snapshot)
            }
          }, error = function(error) {
            try(nm_project_delete(workspace, created$id), silent = TRUE)
            stop(error)
          })
          state$project <- created$id
          invalidate_ai_context()
          state$report_design <- NULL
          state$snapshot <- initial$snapshot %||% NULL
          state$run <- NULL
          state$diagnostics <- list()
          state$model <- initial$model %||% NULL
          state$draft_outputs <- NULL
          state$data <- initial$data %||% NULL
          reset_lazy_payloads()
          update_fit(NULL)
          structure(list(), class = "liber_gui_validation")
        }))
      }
      if (identical(action, "project_save")) {
        return(record("New model version saved", {
          if (is.null(workspace) || is.null(state$project)) {
            .nm_stop("Create or open a project before saving a version.")
          }
          if (is.null(state$model)) .nm_stop("Create or load a model before saving a version.")
          state$snapshot <- nm_project_save(
            workspace, state$project, state$model, state$data, NULL,
            label = as.character(event$label %||% "")
          )
          state$run <- NULL
          state$diagnostics <- list()
          reset_lazy_payloads(data = FALSE)
          structure(list(), class = "liber_gui_validation")
        }))
      }
      if (identical(action, "project_open")) {
        return(record("Project version opened", {
          if (is.null(workspace)) .nm_stop("Launch with `workspace` to open projects.")
          project_id <- as.character(event$id)
          records <- nm_project_list(workspace, project_id)
          versions <- records[records$entry_type == "version", , drop = FALSE]
          if (!nrow(versions)) {
            state$project <- project_id
            state$report_design <- saved_report_design(project_id)
            state$snapshot <- NULL
            state$run <- NULL
            state$model <- NULL
            state$draft_outputs <- NULL
            state$data <- NULL
            state$diagnostics <- list()
            reset_lazy_payloads()
            invalidate_ai_context()
            update_fit(NULL)
            structure(list(), class = "liber_gui_validation")
          } else {
            opened <- load_snapshot(
              project_id, as.character(event$snapshot %||% versions$id[[1L]])
            )
            opened$result %||% structure(list(), class = "liber_gui_validation")
          }
        }))
      }
      if (identical(action, "run_open")) {
        return(record("Model run opened", {
          if (is.null(workspace)) .nm_stop("Launch with `workspace` to open runs.")
          opened <- load_run(as.character(event$id %||% state$project), as.character(event$run))
          opened$result
        }))
      }
      if (identical(action, "project_copy")) {
        return(record("Model version copied", {
          if (is.null(workspace)) .nm_stop("Launch with `workspace` to copy versions.")
          project_id <- as.character(event$id %||% state$project)
          source <- nm_project_load(
            workspace, project_id,
            as.character(event$snapshot %||% state$snapshot %||% "latest")
          )
          copied_model <- source$model
          if (isTRUE(event$updateInits)) {
            fitted <- if (inherits(source$result, "nm_fit")) source$result else state$fit
            if (!inherits(fitted, "nm_fit")) {
              .nm_stop("No estimation result is available to update initial values.")
            }
            if (length(fitted$theta) != nrow(copied_model$THETAS) ||
                length(fitted$sigma) != nrow(copied_model$SIGMAS) ||
                length(fitted$omega) != nrow(copied_model$OMEGAS)) {
              .nm_stop("The fitted parameters do not match the source model version.")
            }
            copied_model$THETAS$Value <- fitted$theta
            copied_model$SIGMAS$Value <- fitted$sigma
            copied_model$OMEGAS$Value <- fitted$omega
          }
          state$snapshot <- nm_project_save(
            workspace, project_id, copied_model, source$data, NULL,
            label = paste(source$label, "copy")
          )
          load_snapshot(project_id, state$snapshot)
          structure(list(), class = "liber_gui_validation")
        }))
      }
      if (identical(action, "project_compare")) {
        return(record("Estimation runs compared", {
          if (is.null(workspace) || is.null(state$project)) {
            .nm_stop("Open a workspace project before comparing runs.")
          }
          ids <- unique(as.character(unlist(event$runs %||% event$snapshots %||% character())))
          if (length(ids) != 2L) .nm_stop("Select exactly two estimation runs to compare.")
          entries <- lapply(ids, function(id) nm_project_load(workspace, state$project, id))
          fits <- lapply(entries, `[[`, "result")
          if (any(!vapply(fits, inherits, logical(1), "nm_fit"))) {
            .nm_stop("Both selected runs must contain estimation results.")
          }
          labels <- make.unique(vapply(entries, function(entry) entry$label, character(1)))
                  gof_frames <- lapply(fits, nm_gof)
                  parameter_vector <- function(fit) {
                    .liber_gui_parameter_values(fit)
                  }
          vectors <- lapply(fits, parameter_vector)
          parameter_names <- unique(unlist(lapply(vectors, names), use.names = FALSE))
          parameters <- data.frame(Parameter = parameter_names, stringsAsFactors = FALSE)
          for (i in seq_along(vectors)) {
            parameters[[labels[[i]]]] <- unname(vectors[[i]][parameter_names])
            covariance <- fits[[i]]$covariance
            if (!is.null(covariance$se)) {
              parameters[[paste(labels[[i]], "SE")]] <- unname(covariance$se[parameter_names])
              parameters[[paste(labels[[i]], "RSE")]] <- unname(covariance$rse[parameter_names])
            }
          }
          gof_summary <- function(fit, gof) {
            observed <- gof$EVID == 0L & gof$MDV == 0L & is.finite(gof$DV)
            gof <- gof[observed, , drop = FALSE]
            n_parameters <- sum(!fit$model$THETAS$FIX) +
              sum(!fit$model$SIGMAS$FIX) + sum(!fit$model$OMEGAS$FIX)
            n_observations <- nrow(gof)
            c(
              OFV = fit$objective,
              AIC = fit$objective + 2 * n_parameters,
              BIC = fit$objective + log(max(1, n_observations)) * n_parameters,
              `Free parameters` = n_parameters,
              Observations = n_observations,
              `Population RMSE` = sqrt(mean((gof$DV - gof$PRED)^2, na.rm = TRUE)),
              `Individual RMSE` = sqrt(mean((gof$DV - gof$IPRED)^2, na.rm = TRUE)),
              `Mean WRES` = mean(gof$WRES, na.rm = TRUE),
              `SD WRES` = stats::sd(gof$WRES, na.rm = TRUE),
              `Mean IWRES` = mean(gof$IWRES, na.rm = TRUE),
              `SD IWRES` = stats::sd(gof$IWRES, na.rm = TRUE),
              `Mean CWRES` = mean(gof$CWRES, na.rm = TRUE),
              `SD CWRES` = stats::sd(gof$CWRES, na.rm = TRUE)
            )
          }
          summaries <- Map(gof_summary, fits, gof_frames)
          metric_names <- unique(unlist(lapply(summaries, names), use.names = FALSE))
          gof <- data.frame(Metric = metric_names, stringsAsFactors = FALSE)
          for (i in seq_along(summaries)) {
            gof[[labels[[i]]]] <- unname(summaries[[i]][metric_names])
          }
          runs <- data.frame(
            Run = labels,
            Method = vapply(fits, function(fit) fit$method, character(1)),
            Objective = vapply(fits, function(fit) fit$objective, numeric(1)),
            Convergence = vapply(fits, function(fit) fit$convergence, integer(1)),
            stringsAsFactors = FALSE
          )
          plots <- list(gof = unname(Map(function(label, frame) {
            list(
              label = label,
              fit = list(
                available = TRUE, gof_loaded = TRUE,
                gof = .liber_gui_rows(frame, 5000L)
              )
            )
          }, labels, gof_frames)))
          run_diagnostics <- lapply(ids, function(id) {
            nm_project_load_diagnostics(workspace, state$project, id)
          })
          for (kind in c("vpc", "vpc_categorical", "vpc_count", "vpc_tte",
                         "vpc_competing",
                         "vpc_recurrent", "npde", "npc")) {
            if (all(vapply(run_diagnostics, function(item) !is.null(item[[kind]]), logical(1)))) {
              plots[[kind]] <- unname(Map(function(label, item) {
                list(label = label, result = .liber_gui_result(item[[kind]]))
              }, labels, run_diagnostics))
            }
          }
          state$comparison_open <- TRUE
          structure(
            list(
              id = paste0(format(Sys.time(), "%Y%m%d%H%M%OS3"), "-", sample.int(999999L, 1L)),
              parameters = parameters, gof = gof, runs = runs, plots = plots
            ),
            class = "liber_gui_comparison"
          )
        }))
      }
      if (identical(action, "project_delete_snapshot")) {
        return(record("Model version deleted", {
          if (is.null(workspace)) .nm_stop("Launch with `workspace` to delete versions.")
          project_id <- as.character(event$id %||% state$project)
          nm_project_delete_snapshot(workspace, project_id, as.character(event$snapshot))
          remaining <- nm_project_list(workspace, project_id)
          remaining <- remaining[remaining$entry_type == "version", , drop = FALSE]
          if (nrow(remaining)) {
            load_snapshot(project_id, remaining$id[[1L]])
          } else {
            state$snapshot <- NULL
            state$run <- NULL
            state$model <- NULL
            state$draft_outputs <- NULL
            state$data <- NULL
            state$diagnostics <- list()
            reset_lazy_payloads()
            update_fit(NULL)
          }
          structure(list(), class = "liber_gui_validation")
        }))
      }
      if (identical(action, "project_delete_run")) {
        return(record("Model run deleted", {
          if (is.null(workspace)) .nm_stop("Launch with `workspace` to delete runs.")
          project_id <- as.character(event$id %||% state$project)
          nm_project_delete_snapshot(workspace, project_id, as.character(event$run))
          state$run <- NULL
          state$diagnostics <- list()
          reset_lazy_payloads(data = FALSE)
          if (!is.null(state$snapshot)) load_snapshot(project_id, state$snapshot)
          structure(list(), class = "liber_gui_validation")
        }))
      }
      if (identical(action, "project_delete")) {
        return(record("Project deleted", {
          if (is.null(workspace)) .nm_stop("Launch with `workspace` to delete projects.")
          if (!identical(as.character(event$confirmation %||% ""), "YES")) {
            .nm_stop('Type "YES" in the confirmation field before deleting a project.')
          }
          nm_project_delete(workspace, as.character(event$id %||% state$project))
          invalidate_ai_context()
          state$project <- NULL
          state$report_design <- NULL
          state$snapshot <- NULL
          state$run <- NULL
          state$model <- NULL
          state$draft_outputs <- NULL
          state$data <- NULL
          state$diagnostics <- list()
          reset_lazy_payloads()
          update_fit(NULL)
          structure(list(), class = "liber_gui_validation")
        }))
      }
      if (identical(action, "queue_select")) {
        return(record("Execution queue selected", {
          id <- as.character(event$id %||% "local")
          if (!identical(id, "local") && is.null(state$remote_queues[[id]])) {
            .nm_stop("Unknown execution queue.")
          }
          state$queue_id <- id
          state$hidden_jobs <- character()
          state$job_log <- character()
          state$selected_job <- NULL
          save_client_settings()
          refresh_jobs(start = TRUE)
          structure(list(), class = "liber_gui_validation")
        }, update_result = FALSE))
      }
      if (identical(action, "queue_add")) {
        return(record("Remote server connected", {
          if (!requireNamespace("LibeRties", quietly = TRUE)) {
            .nm_stop("LibeRties is required for remote execution.")
          }
          id <- as.character(event$id %||% "")
          if (!nzchar(id)) id <- paste0(
            "remote-", format(Sys.time(), "%Y%m%d%H%M%S"), "-",
            sprintf("%04d", sample.int(9999L, 1L))
          )
          if (identical(id, "local")) .nm_stop("Remote queue id cannot be `local`.")
          previous <- state$remote_config[[id]]
          token <- as.character(event$token %||% "")
          if (!nzchar(token) && !is.null(previous)) token <- as.character(previous$token %||% "")
          if (!nzchar(token)) .nm_stop("Enter the server bearer token.")
          remote <- LibeRties::ls_remote(
            as.character(event$url), token, timeout = 30
          )
          authentication <- remote$authenticate()
          server_name <- as.character(event$name %||% "Remote server")
          queues <- state$remote_queues
          queues[[id]] <- remote
          state$remote_queues <- queues
          metadata <- state$remote_meta
          metadata[[id]] <- list(
            name = server_name,
            url = remote$url, user = as.character(authentication$username %||% "")
          )
          state$remote_meta <- metadata
          config <- state$remote_config
          config[[id]] <- list(
            name = server_name, url = remote$url, token = token, timeout = 30,
            user = as.character(authentication$username %||% "")
          )
          state$remote_config <- config
          state$queue_id <- id
          state$hidden_jobs <- character()
          save_client_settings()
          refresh_jobs()
          structure(list(), class = "liber_gui_validation")
        }, update_result = FALSE))
      }
      if (identical(action, "queue_remove")) {
        return(record("Remote server removed", {
          id <- as.character(event$id)
          if (identical(id, "local")) .nm_stop("The local execution target cannot be removed.")
          queues <- state$remote_queues
          queues[[id]] <- NULL
          state$remote_queues <- queues
          metadata <- state$remote_meta
          metadata[[id]] <- NULL
          state$remote_meta <- metadata
          config <- state$remote_config
          config[[id]] <- NULL
          state$remote_config <- config
          state$queue_id <- "local"
          state$jobs <- data.frame()
          state$job_log <- character()
          state$selected_job <- NULL
          save_client_settings()
          refresh_jobs(start = TRUE)
          structure(list(), class = "liber_gui_validation")
        }, update_result = FALSE))
      }
      if (identical(action, "jobs_refresh")) {
        return(record("Job queue refreshed", {
          refresh_jobs(start = TRUE)
          structure(list(), class = "liber_gui_validation")
        }, update_result = FALSE))
      }
      if (identical(action, "jobs_clear")) {
        return(record("Finished jobs hidden from this view", {
          jobs <- as.data.frame(state$jobs, stringsAsFactors = FALSE)
          terminal <- if (nrow(jobs) && all(c("id", "status") %in% names(jobs))) {
            jobs$id[jobs$status %in% c("completed", "failed", "cancelled")]
          } else character()
          state$hidden_jobs <- unique(c(state$hidden_jobs, terminal))
          refresh_jobs()
          structure(list(), class = "liber_gui_validation")
        }, update_result = FALSE))
      }
      if (identical(action, "job_cancel")) {
        return(record("Job cancellation requested", {
          q <- active_queue()
          if (is.null(q)) .nm_stop("No queued execution target is selected.")
          q$cancel(as.character(event$id))
          refresh_jobs(start = TRUE)
          structure(list(), class = "liber_gui_validation")
        }, update_result = FALSE))
      }
      if (identical(action, "job_result")) {
        return(record("Saved model run opened", {
          requested_queue <- as.character(event$queueId %||% state$queue_id)
          if (!identical(requested_queue, state$queue_id)) {
            state$queue_id <- requested_queue
            state$hidden_jobs <- character()
          }
          q <- active_queue(requested_queue)
          if (is.null(q)) .nm_stop("No queued execution target is selected.")
          id <- as.character(event$id)
          key <- job_context_key(id, requested_queue)
          context <- state$job_context[[key]] %||% state$job_context[[id]]
          if (!is.null(context) &&
              !nzchar(as.character(context$run_id %||% ""))) {
            context <- materialize_job_result(context, q, id)
            contexts <- state$job_context
            contexts[[key]] <- context
            contexts[[id]] <- NULL
            state$job_context <- contexts
            save_client_settings()
          }
          if (!is.null(context) &&
              nzchar(as.character(context$run_id %||% ""))) {
            opened <- load_run(context$project, context$run_id)
            return(opened$result)
          }
          q$result(id)
        }))
      }
      if (identical(action, "job_select")) {
        value <- tryCatch({
          requested_queue <- as.character(event$queueId %||% state$queue_id)
          if (!identical(requested_queue, "local") &&
              is.null(state$remote_queues[[requested_queue]])) {
            .nm_stop("The execution queue for this model run is not configured.")
          }
          if (!identical(requested_queue, state$queue_id)) {
            state$queue_id <- requested_queue
            state$hidden_jobs <- character()
            save_client_settings()
          }
          q <- active_queue(requested_queue)
          if (is.null(q)) .nm_stop("No queued execution target is selected.")
          id <- as.character(event$id)
          state$selected_job <- id
          refresh_jobs(start = inherits(q, "LibeRQueue"))
          state$job_log <- with_ui_remote_timeout(
            q, FALSE,
            c("--- stdout ---", q$logs(id, stream = "stdout"),
              "--- stderr ---", q$logs(id, stream = "stderr"))
          )
          NULL
        }, error = identity)
        if (inherits(value, "error")) append_log(conditionMessage(value), "error")
        return(invisible(NULL))
      }
      invisible(NULL)
    }, ignoreInit = TRUE)
  }
  app <- shiny::shinyApp(ui, server)
  if (is.null(launch.browser)) return(app)
  if (is.null(dots$port) && is.null(getOption("shiny.port"))) {
    # Browser model storage is origin-scoped. A stable default port lets the
    # downloaded WebLLM weights persist across ordinary GUI sessions.
    dots$port <- 38764L
  }
  do.call(shiny::runApp, c(list(appDir = app, launch.browser = launch.browser), dots))
}
