test_that("react workbench serializes model and data without live pointers", {
  model <- nm_model(
    INPUT = c("ID", "TIME"), ADVAN = 1,
    PRED = "CL=THETA(1)\nV=THETA(2)\nS1=V", ERROR = "Y=F",
    THETAS = data.frame(THETA = 1:2, Value = c(2, 20))
  )
  widget <- liber_workbench(
    model, data.frame(ID = 1, TIME = 0), height = "600px"
  )
  expect_s3_class(widget, "htmlwidget")
  expect_equal(widget$x$tag$attribs$model$advan, 1)
  expect_equal(widget$x$tag$attribs$model$name, "ADVAN1 model")
  expect_equal(widget$x$tag$attribs$dataset$subjects, 1)
  expect_true("TIME" %in% widget$x$tag$attribs$dataset$numeric_columns)
  expect_equal(length(widget$x$tag$attribs$dataset$plot_rows), 1)
  # reactR omits zero-length R lists; the JavaScript boundary treats this as [].
  expect_null(widget$x$tag$attribs$jobs)
})

test_that("legacy workbench layout and workflow controls are present", {
  source <- paste(readLines(
    system.file("htmlwidgets", "liberWorkbench.js", package = "LibeRation"),
    warn = FALSE
  ), collapse = "\n")
  style <- paste(readLines(
    system.file("htmlwidgets", "liberWorkbench.css", package = "LibeRation"),
    warn = FALSE
  ), collapse = "\n")
  expect_match(source, 'label:"Home"', fixed = TRUE)
  expect_match(source, 'label:"Jobs"', fixed = TRUE)
  expect_match(source, 'label:"Data"', fixed = TRUE)
  expect_match(source, '"project_compare"', fixed = TRUE)
  expect_match(source, '"Simulate"', fixed = TRUE)
  expect_match(source, '"Estimate"', fixed = TRUE)
  expect_match(source, 'Run covariance step after estimation', fixed = TRUE)
  expect_match(source, 'label:"Covariance"', fixed = TRUE)
  expect_match(source, 'label:"Posterior"', fixed = TRUE)
  expect_match(source, 'Posterior SDs, posterior CVs and 95% credible intervals', fixed = TRUE)
  expect_match(source, 'Observed marginal information uses deterministic Gauss-Hermite', fixed = TRUE)
  expect_match(source, '"Run diagnostic"', fixed = TRUE)
  expect_match(source, '"New project"', fixed = TRUE)
  expect_match(source, '"Empty project"', fixed = TRUE)
  expect_match(source, '"Initial model version"', fixed = TRUE)
  expect_match(source, '"Dose amounts (TIME AMT per line, or AMT only)"', fixed = TRUE)
  expect_match(source, '"Edit server"', fixed = TRUE)
  expect_match(source, 'Manual X breaks', fixed = TRUE)
  expect_match(source, 'function CodeEditor', fixed = TRUE)
  expect_match(style, 'lw-syntax-parameter', fixed = TRUE)
  expect_match(style, 'lw-syntax-definition', fixed = TRUE)
  expect_match(style, 'lw-syntax-function', fixed = TRUE)
  expect_match(source, 'OMEGA matrix', fixed = TRUE)
  expect_match(source, '"OMEGA("', fixed = TRUE)
  expect_match(source, 'Sequential estimation', fixed = TRUE)
  expect_match(source, 'Generated run columns', fixed = TRUE)
  expect_match(source, 'Print gradients every N (0 = off)', fixed = TRUE)
  expect_match(source, 'GQ (adaptive Gauss-Hermite)', fixed = TRUE)
  expect_match(source, 'gqOrder', fixed = TRUE)
  expect_match(source, 'gqGrid', fixed = TRUE)
  expect_match(source, 'gqLevel', fixed = TRUE)
  expect_match(source, 'gqAdaptive', fixed = TRUE)
  expect_match(source, 'gqMaxPoints', fixed = TRUE)
  expect_match(source, 'Smolyak sparse grid', fixed = TRUE)
  expect_match(source, 'Open the saved model run and its results', fixed = TRUE)
  expect_match(source, 'comparison_close', fixed = TRUE)
  expect_match(source, 'page_change', fixed = TRUE)
  expect_match(source, 'lw-app-icon', fixed = TRUE)
  expect_match(source, 'lw-app-version', fixed = TRUE)
  expect_match(source, 'Estimation priors', fixed = TRUE)
  expect_match(source, 'lw-modal-comparison', fixed = TRUE)
  expect_match(source, 'lw-run-flag-cov', fixed = TRUE)
  expect_match(source, 'load_payload', fixed = TRUE)
  expect_match(source, 'function HmmPane', fixed = TRUE)
  expect_match(source, 'Retrospective smoothed', fixed = TRUE)
  expect_match(source, 'kind:"hmm"', fixed = TRUE)
  expect_match(source, 'function KalmanPane', fixed = TRUE)
  expect_match(source, 'kind:"kalman"', fixed = TRUE)
  expect_match(source, 'Stratify by', fixed = TRUE)
  expect_match(source, 'function ComparisonPlots', fixed = TRUE)
  expect_match(source, 'lw-button-danger-ghost', fixed = TRUE)
  expect_match(source, 'Type "YES" to confirm', fixed = TRUE)
  expect_match(source, 'confirmation:deleteConfirmation[0]', fixed = TRUE)
  expect_match(source, 'deleteConfirmation[0]!=="YES"', fixed = TRUE)
  expect_match(source, 'Selected project:', fixed = TRUE)
  expect_match(source, 'Selected model version:', fixed = TRUE)
  expect_match(source, 'Selected model run:', fixed = TRUE)
  expect_match(source, 'Dataset metadata:', fixed = TRUE)
  expect_match(source, 'ai_context_request', fixed = TRUE)
  expect_match(source, 'Saved project evidence', fixed = TRUE)
  expect_match(source, 'localAICancelPurpose', fixed = TRUE)
  expect_match(source, 'localAIRetryPending', fixed = TRUE)
  expect_match(source, 'DXGI_ERROR_DEVICE_', fixed = TRUE)
  expect_match(source, 'No usable WebGPU adapter is available', fixed = TRUE)
  expect_match(source, 'Reset local AI', fixed = TRUE)
  expect_match(source, 'localAIBudgetMessages', fixed = TRUE)
  expect_match(source, 'function AIContextSelect', fixed = TRUE)
  expect_match(source, 'help_context', fixed = TRUE)
  expect_match(source, 'report_context', fixed = TRUE)
  expect_match(source, 'context compacted', fixed = TRUE)
  expect_match(source, 'helpQuestionScope', fixed = TRUE)
  expect_match(source, 'scope:projectScope', fixed = TRUE)
  expect_match(source, 'Evidence scope selected for this question', fixed = TRUE)
  expect_match(source, 'report_ai_context_request', fixed = TRUE)
  expect_match(source, 'reportSelections', fixed = TRUE)
  expect_match(source, 'reportEvidenceText', fixed = TRUE)
  expect_match(source, 'Treat every listed saved run as available evidence', fixed = TRUE)
  expect_match(source, 'Save location', fixed = TRUE)
  expect_match(source, 'report_directory_choose', fixed = TRUE)
  expect_match(source, 'Generation stopped before a response was produced.', fixed = TRUE)
  worker_source <- paste(readLines(
    system.file("htmlwidgets", "liber-ai-worker.js", package = "LibeRation"),
    warn = FALSE
  ), collapse = "\n")
  expect_match(worker_source, 'isDisposedFailure', fixed = TRUE)
  expect_match(worker_source, 'recoverEngine', fixed = TRUE)
  expect_match(worker_source, 'Refreshing the local GPU session', fixed = TRUE)
  expect_match(worker_source, 'DXGI_ERROR_DEVICE_', fixed = TRUE)
  expect_match(worker_source, 'requestAdapter', fixed = TRUE)
  expect_match(worker_source, 'powerPreference: "high-performance"', fixed = TRUE)
  expect_match(worker_source, 'isContextFailure', fixed = TRUE)
  expect_match(worker_source, 'compactContextRequest', fixed = TRUE)
  expect_match(worker_source, 'context_compacted', fixed = TRUE)
  expect_match(worker_source, 'context_window_size: contextWindow', fixed = TRUE)
  expect_match(worker_source, 'contextFallbacks', fixed = TRUE)
  expect_match(worker_source, 'isMemoryAllocationFailure', fixed = TRUE)
  server_source <- paste(deparse(body(liber_gui)), collapse = "\n")
  expect_match(server_source, 'q\\$logs\\(selected,\\s+stream = "stdout"', perl = TRUE)
  expect_match(server_source, 'q\\$logs\\(id,\\s+stream = "stdout"', perl = TRUE)
  expect_match(server_source, 'event$confirmation', fixed = TRUE)
  estimator_source <- paste(deparse(body(LibeRation:::.liber_estimation_arguments)), collapse = "\n")
  expect_match(estimator_source, 'event$gqGrid', fixed = TRUE)
  expect_match(estimator_source, 'event$gqLevel', fixed = TRUE)
})

test_that("Help AI receives compact saved-run evidence on demand", {
  fixture <- estimation_fixture()
  fit <- nm_est(fixture$model, fixture$data, method = "FO", maxit = 1)
  root <- tempfile("gui-ai-context-")
  workspace <- nm_workspace(root)
  project <- nm_project_create(workspace, "AI evidence")
  version <- nm_project_save(
    workspace, project$id, fixture$model, fixture$data, label = "Base model"
  )
  run <- nm_project_save_run(
    workspace, project$id, version, fit, label = "FO estimation"
  )
  nm_project_save_diagnostics(
    workspace, project$id, run, list(gof = nm_gof(fit))
  )

  context <- LibeRation:::.liber_gui_ai_context(
    workspace, project$id, selected_run = run
  )
  expect_true(context$available)
  expect_equal(context$run_count, 1L)
  expect_equal(context$runs[[1L]]$label, "FO estimation")
  expect_equal(context$runs[[1L]]$objective, fit$objective)
  expect_true(context$runs[[1L]]$diagnostics$gof)
  expect_true(length(context$runs[[1L]]$parameters) > 0L)
  expect_false(any(c("data", "output", "gof") %in% names(context$runs[[1L]])))

  index_context <- LibeRation:::.liber_gui_ai_context(
    workspace, project$id, selected_run = run, detail = "index"
  )
  expect_equal(index_context$scope, "index")
  expect_false("parameters" %in% names(index_context$runs[[1L]]))
  expect_lt(as.numeric(object.size(index_context)), as.numeric(object.size(context)))

  app <- liber_gui(
    workspace = workspace, project = project$id, snapshot = version,
    launch.browser = NULL
  )
  server_function <- app[["serverFuncSource"]]()
  supplied <- NULL
  shiny::testServer(server_function, {
    session$setInputs(liber_workbench_event = list(action = "noop", nonce = 0))
    session$flushReact()
    session$setInputs(liber_workbench_event = list(
      action = "ai_context_request", project = project$id,
      requestId = "help-request-1", scope = "index", nonce = 1
    ))
    session$flushReact()
    supplied <<- shiny::isolate(state$ai_context)
  })
  expect_equal(supplied$request_id, "help-request-1")
  expect_true(supplied$available)
  expect_equal(supplied$scope, "index")
  expect_equal(supplied$runs[[1L]]$id, run)
  expect_false("parameters" %in% names(supplied$runs[[1L]]))

  second_run <- nm_project_save_run(
    workspace, project$id, version, fit, label = "Second FO estimation"
  )
  report_context <- LibeRation:::.liber_gui_report_ai_context(
    workspace, project$id, c(run, second_run)
  )
  expect_true(report_context$available)
  expect_equal(report_context$run_ids, c(run, second_run))
  expect_length(report_context$runs, 2L)
  expect_equal(report_context$runs[[1L]]$model$advan, fixture$model$ADVAN)
  expect_equal(report_context$runs[[1L]]$data$subjects,
               length(unique(fit$data$.ID_INDEX)))
  expect_true(length(report_context$runs[[1L]]$parameters) > 0L)
  expect_true(is.list(report_context$runs[[1L]]$gof_summary))
})

test_that("large data and diagnostic payloads are lazy", {
  data <- data.frame(
    ID = rep(1:100, each = 20), TIME = rep(seq(0, 19), 100),
    EVID = 0L, AMT = 0, DV = 1, MDV = 0L
  )
  metadata <- LibeRation:::.liber_gui_data(data, include_rows = FALSE)
  expect_true(metadata$loaded)
  expect_false(metadata$payload_loaded)
  expect_length(metadata$plot_rows, 0L)
  expect_equal(metadata$records, 2000L)
  expect_true(all(c("TIME", "DV") %in% metadata$numeric_columns))

  vpc <- structure(list(nsim = 20L), class = "nm_vpc")
  diagnostics <- LibeRation:::.liber_gui_diagnostics(list(vpc = vpc), payload = character())
  expect_true(diagnostics$available$vpc)
  expect_null(diagnostics$vpc)
})

test_that("HMM GUI payload exposes decoder trajectories and sequence evidence", {
  decoded <- data.frame(
    ID = c("A", "A"), TIME = 0:1, DVID = 1L,
    HMM_ROW_NLL = c(0.2, 0.3),
    HMM_FILTER_STATE_INDEX = c(1L, 2L),
    HMM_FILTER_STATE = c("low", "high"),
    HMM_SMOOTH_STATE_INDEX = c(1L, 1L),
    HMM_SMOOTH_STATE = c("low", "low"),
    HMM_VITERBI_STATE_INDEX = c(1L, 1L),
    HMM_VITERBI_STATE = c("low", "low"),
    HMM_FILTER_PROB_low = c(0.8, 0.4),
    HMM_FILTER_PROB_high = c(0.2, 0.6),
    HMM_SMOOTH_PROB_low = c(0.9, 0.7),
    HMM_SMOOTH_PROB_high = c(0.1, 0.3)
  )
  attr(decoded, "states") <- c("low", "high")
  attr(decoded, "method") <- "all"
  attr(decoded, "eta_type") <- "individual"
  attr(decoded, "log_likelihood") <- -0.5
  attr(decoded, "sequence_summary") <- data.frame(
    ID = "A", DVID = 1L, LOG_LIKELIHOOD = -0.5,
    VITERBI_LOG_JOINT = -0.7, VITERBI_LOG_POSTERIOR = -0.2
  )
  class(decoded) <- c("nm_hmm_decode", class(decoded))

  payload <- LibeRation:::.liber_gui_hmm(decoded)
  expect_true(payload$available)
  expect_true(payload$loaded)
  expect_equal(payload$observations, 2L)
  expect_equal(payload$sequences, 1L)
  expect_equal(vapply(payload$states, `[[`, character(1), "label"),
               c("low", "high"))
  expect_equal(payload$rows[[1L]]$SUBJECT, "A")
  expect_equal(payload$rows[[2L]]$HMM_SMOOTH_STATE, "low")
  expect_equal(payload$sequence_summary[[1L]]$VITERBI_LOG_POSTERIOR, -0.2)

  unloaded <- LibeRation:::.liber_gui_hmm(available = TRUE)
  expect_true(unloaded$available)
  expect_false(unloaded$loaded)
  expect_length(unloaded$rows, 0L)
})

test_that("linear state-space results have a lazy GUI payload", {
  decoded <- data.frame(
    ID = c("A", "A"), TIME = c(0, 1), DVID = 1L, DV = c(1.2, 1.4),
    KF_PRED_exposure = c(0, 0.8),
    KF_FILTER_exposure = c(0.7, 1.0),
    KF_SMOOTH_exposure = c(0.9, 1.0),
    KF_FILTER_SD_exposure = c(0.5, 0.4),
    KF_SMOOTH_SD_exposure = c(0.3, 0.4),
    KF_INNOVATION = c(1.2, 0.6),
    KF_INNOVATION_VARIANCE = c(1.44, 0.36),
    KF_ROW_NLL = c(1.1, 0.8)
  )
  attr(decoded, "states") <- "exposure"
  attr(decoded, "eta_type") <- "individual"
  attr(decoded, "log_likelihood") <- -1.9
  class(decoded) <- c("nm_kalman_decode", class(decoded))

  payload <- LibeRation:::.liber_gui_kalman(decoded)
  expect_true(payload$available)
  expect_true(payload$loaded)
  expect_equal(payload$observations, 2L)
  expect_equal(payload$sequences, 1L)
  expect_equal(payload$states[[1L]]$label, "exposure")
  expect_equal(payload$rows[[1L]]$KF_STANDARDIZED_INNOVATION, 1)

  unloaded <- LibeRation:::.liber_gui_kalman(available = TRUE)
  expect_true(unloaded$available)
  expect_false(unloaded$loaded)
  expect_length(unloaded$rows, 0L)
})

test_that("selected run outputs are loaded lazily into Data explorer", {
  data <- data.frame(
    ID = c(1, 1), TIME = c(0, 1), EVID = c(1, 0), AMT = c(100, 0),
    DV = c(NA, 4.5), MDV = c(1, 0)
  )
  output <- data.frame(.ROW = 1:2, CL = c(2, 2), CWRES = c(NA, 0.1))
  metadata <- LibeRation:::.liber_gui_data(
    data, include_rows = FALSE, run_output = output
  )
  expect_false("CL" %in% metadata$columns)

  payload <- LibeRation:::.liber_gui_data(
    data, include_rows = TRUE, run_output = output
  )
  expect_true(all(c("CL", "CWRES") %in% payload$columns))
  expect_equal(vapply(payload$plot_rows, `[[`, numeric(1), "CL"), c(2, 2))
})

test_that("draft model validation sees unsaved code and generated outputs", {
  model <- LibeRation:::.liber_model_template(1L)
  event <- list(
    pred = paste(model$PRED, "K = CL / V", sep = "\n"),
    output = c("PRED", "K")
  )
  draft <- LibeRation:::.liber_model_from_event(model, event)
  expect_true("K" %in% nm_model_outputs(draft)$name)
  expect_equal(draft$OUTPUT, c("PRED", "K"))
  expect_s3_class(nm_compile(draft), "NMEngine")
})

test_that("draft validation synchronizes THETA ETA and residual parameter tables", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "CMT", "DV", "MDV"),
    ADVAN = 1,
    PRED = "CL=THETA(1)*exp(ETA(1)); V=THETA(2); S1=V",
    ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(THETA = 1:2, Value = c(2, 20)),
    OMEGAS = data.frame(OMEGA = 1, Value = 0.1),
    SIGMAS = data.frame(SIGMA = 1, Value = 0.1)
  )
  draft <- LibeRation:::.liber_model_from_event(model, list(
    pred = paste(model$PRED, "BIO=THETA(3)*exp(ETA(2))", sep = "\n"),
    error = "Y=F*(1+ERR(1))+ERR(2)"
  ))
  expect_equal(nrow(draft$THETAS), 3L)
  expect_equal(draft$n_eta, 2L)
  expect_equal(nrow(draft$SIGMAS), 2L)
  expect_equal(draft$THETAS$Value[[3L]], 1)
  expect_equal(draft$OMEGAS$Value[[2L]], 0.1)
  expect_equal(draft$SIGMAS$Value[[2L]], 0.1)
  expect_s3_class(nm_compile(draft), "NMEngine")

  full <- LibeRation:::.liber_model_from_event(model, list(
    pred = paste(model$PRED, "BIO=THETA(3)*exp(ETA(2))", sep = "\n"),
    omega_structure = "full"
  ))
  expect_equal(nrow(full$OMEGAS), 3L)
  expect_equal(full$OMEGAS[, c("ROW", "COL")],
               data.frame(ROW = c(1L, 2L, 2L), COL = c(1L, 1L, 2L)))
})

test_that("GUI parameter estimates use THETA OMEGA SIGMA presentation order", {
  values <- LibeRation:::.liber_gui_parameter_values(list(
    theta = c(1, 2), omega = c(0.1, 0.2), sigma = 0.3
  ))
  expect_equal(
    names(values),
    c("THETA1", "THETA2", "OMEGA1", "OMEGA2", "SIGMA1")
  )
})

test_that("validation payload returns synchronized draft parameter tables", {
  validation <- structure(list(
    kind = "model_validation", nonce = "validation-1",
    parameters = list(theta = list(list(THETA = 1, Value = 1)),
                      omega = list(), sigma = list())
  ), class = "liber_gui_validation")
  payload <- LibeRation:::.liber_gui_result(validation)
  expect_equal(payload$status, "validated")
  expect_equal(payload$kind, "model_validation")
  expect_equal(payload$parameters$theta[[1L]]$THETA, 1)
})

test_that("GUI detects user likelihood edits and exposes compatible methods", {
  model <- LibeRation:::.liber_model_template(1L)
  draft <- LibeRation:::.liber_model_from_event(model, list(
    pred = paste(model$PRED, "P = 1 / (1 + exp(-THETA(1))); F = P", sep = "\n"),
    error = "LOGLIK = log(pmax(ifelse(DV == 1, F, 1 - F), 1e-12))"
  ))
  payload <- LibeRation:::.liber_gui_model(draft)
  expect_identical(payload$likelihood_type, "likelihood")
  expect_identical(draft$likelihood_output, "LOGLIK")

  source <- paste(readLines(
    system.file("htmlwidgets", "liberWorkbench.js", package = "LibeRation"),
    warn = FALSE
  ), collapse = "\n")
  expect_match(source, "User-defined likelihood detected", fixed = TRUE)
  expect_match(source, 'estimationMethod[1]("LAPLACE")', fixed = TRUE)
})

test_that("GUI exposes first-class outcomes, templates, and predictive checks", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "DV", "MDV"), ADVAN = 1,
    PRED = "MU=exp(THETA(1)); CL=1; V=1; S1=V; F=MU",
    THETAS = data.frame(THETA = 1, Value = 0, LOWER = -10, UPPER = 10),
    OUTCOMES = nm_outcome("poisson", name = "seizures", prediction = "MU")
  )
  payload <- LibeRation:::.liber_gui_model(model)
  expect_equal(payload$outcomes[[1L]]$family, "poisson")
  expect_equal(payload$outcomes[[1L]]$name, "seizures")
  expect_true(payload$outcomes[[1L]]$generated_error)

  count <- structure(list(
    nsim = 20L, outcome = "seizures", family = "poisson", dvid = NULL,
    observed = data.frame(TIME = 0, MEAN = 1),
    simulated = data.frame(TIME = 0, MEAN_median = 1,
                           MEAN_lower = 0.5, MEAN_upper = 1.5)
  ), class = "nm_vpc_count")
  diagnostics <- LibeRation:::.liber_gui_diagnostics(
    list(vpc_count = count), payload = "vpc_count"
  )
  expect_true(diagnostics$available$vpc_count)
  expect_equal(diagnostics$vpc_count$kind, "vpc_count")

  source <- paste(readLines(
    system.file("htmlwidgets", "liberWorkbench.js", package = "LibeRation"),
    warn = FALSE
  ), collapse = "\n")
  expect_match(source, "First-class outcomes", fixed = TRUE)
  expect_match(source, "Nonlinear elimination", fixed = TRUE)
  expect_match(source, "Count VPC", fixed = TRUE)
  expect_match(source, "Competing-risk VPC", fixed = TRUE)
  expect_match(source, "Recurrent-event VPC", fixed = TRUE)
})

test_that("the initial queue refresh runs inside a reactive isolate", {
  app <- liber_gui(workspace = tempfile("gui-flush-"), queue = FALSE,
                   launch.browser = NULL)
  server_function <- app[["serverFuncSource"]]()
  source <- paste(deparse(body(server_function)), collapse = "\n")
  expect_match(
    source,
    "session\\$onFlushed\\(function\\(\\) \\{\\s+shiny::isolate\\(\\{",
    perl = TRUE
  )
  expect_match(source, "needs_reconciliation", fixed = TRUE)
  expect_match(source, 'identical(state$active_page, "home")', fixed = TRUE)
  expect_match(source, "!isTRUE(state$comparison_open)", fixed = TRUE)
  expect_match(source, 'status %in% c("queued", "running")', fixed = TRUE)
  expect_match(source, "poll_backoff$until", fixed = TRUE)
  expect_match(source, "background = !inherits(q, \"LibeRQueue\")", fixed = TRUE)
})

test_that("client queue settings live in the workspace and survive GUI recreation", {
  skip_if_not_installed("LibeRties")
  root <- tempfile("gui-client-settings-")
  workspace <- nm_workspace(root)
  remotes <- list(team = list(
    name = "Team server", url = "https://server.example", token = "secret-token",
    timeout = 30, user = "scientist"
  ))
  pending <- list(`team::job-1` = list(
    queue_id = "team", job_id = "job-1", project = "demo",
    version = "version-1", type = "estimate", status = "running"
  ))
  path <- LibeRation:::.liber_client_settings_write(
    workspace, "team", remotes, pending_jobs = pending
  )
  expect_true(startsWith(normalizePath(path, winslash = "/"),
                         normalizePath(root, winslash = "/")))
  restored <- LibeRation:::.liber_client_settings_read(workspace)
  expect_equal(restored$selected_queue, "team")
  expect_equal(restored$remotes$team$token, "secret-token")
  expect_equal(restored$pending_jobs[[1L]]$status, "running")

  app <- liber_gui(workspace = workspace, launch.browser = NULL)
  server_function <- app[["serverFuncSource"]]()
  favicon_href <- get("favicon_href", envir = environment(server_function), inherits = TRUE)
  restored_state <- NULL
  shiny::testServer(server_function, {
    restored_state <<- shiny::isolate(list(
      queue_id = state$queue_id,
      remote_meta = state$remote_meta,
      remote_queues = state$remote_queues
    ))
  })
  expect_equal(restored_state$queue_id, "team")
  expect_equal(restored_state$remote_meta$team$name, "Team server")
  expect_s3_class(restored_state$remote_queues$team, "LibeRRemote")
  expect_match(favicon_href, "^liberation-assets-")
  expect_false(startsWith(favicon_href, "data:"))
  expect_lt(nchar(favicon_href), 100L)
})

test_that("GUI model editing supports full OMEGA matrices and priors", {
  model <- LibeRation:::.liber_model_template(2L)
  full <- list(
    list(OMEGA = 1, ROW = 1, COL = 1, Value = 0.1, FIX = FALSE),
    list(OMEGA = 2, ROW = 2, COL = 1, Value = 0.01, FIX = FALSE),
    list(OMEGA = 3, ROW = 2, COL = 2, Value = 0.1, FIX = FALSE),
    list(OMEGA = 4, ROW = 3, COL = 1, Value = 0.01, FIX = FALSE),
    list(OMEGA = 5, ROW = 3, COL = 2, Value = 0.01, FIX = FALSE),
    list(OMEGA = 6, ROW = 3, COL = 3, Value = 0.1, FIX = FALSE)
  )
  omega <- LibeRation:::.liber_parameter_table_update(model$OMEGAS, full, "OMEGA")
  priors <- LibeRation:::.liber_prior_table_update(list(
    list(parameter = "THETA1", distribution = "lognormal", mean = 0, sd = 0.5)
  ))
  arguments <- model[intersect(names(model), names(formals(nm_model)))]
  arguments$OMEGAS <- omega
  arguments$LIK_CONFIG <- nm_lik_config(omega = "full", priors = priors)
  edited <- do.call(nm_model, arguments)
  expect_equal(nrow(edited$OMEGAS), 6L)
  expect_equal(edited$LIK_CONFIG$omega, "full")
  expect_equal(edited$LIK_CONFIG$priors$parameter, "THETA1")
})

test_that("theophylline GUI example is oral ADVAN2 TRANS2 with a delayed peak", {
  model <- LibeRation:::.liber_model_template(2L, trans = 2L)
  expect_equal(model$ADVAN, 2L)
  expect_equal(model$TRANS, 2L)
  expect_equal(model$DOSECMP, 1L)
  expect_equal(model$OBSCMP, 2L)
  data <- LibeRation:::.liber_builtin_dataset(model, "theophylline", n_subjects = 1L)
  simulated <- nm_simulate(model, data, random_effects = FALSE, residual = FALSE)
  observations <- simulated[simulated$EVID == 0L, , drop = FALSE]
  expect_gt(observations$TIME[[which.max(observations$IPRED)]], min(observations$TIME))
})

test_that("GUI example datasets contain reproducible between-subject variability", {
  model <- LibeRation:::.liber_model_template(2L, trans = 2L)
  first <- LibeRation:::.liber_builtin_dataset(
    model, "theophylline", n_subjects = 12L, seed = 417L
  )
  second <- LibeRation:::.liber_builtin_dataset(
    model, "theophylline", n_subjects = 12L, seed = 417L
  )
  truth <- attr(first, "simulation_eta")
  expect_s3_class(truth, "data.frame")
  expect_equal(nrow(truth), 12L)
  expect_true(any(abs(as.matrix(truth[grep("^ETA", names(truth))])) > 0.05))
  expect_equal(first$DV, second$DV)
  expect_equal(attr(first, "simulation_eta"), attr(second, "simulation_eta"))
})

test_that("GUI uses the high-resolution blue LibeRation favicon", {
  favicon_path <- system.file("assets", "favicon.svg", package = "LibeRation")
  expect_true(file.exists(favicon_path))
  favicon <- paste(readLines(favicon_path, warn = FALSE), collapse = "\n")
  expect_match(favicon, 'width="1000"', fixed = TRUE)
  expect_match(favicon, 'height="1000"', fixed = TRUE)
  expect_match(favicon, 'href="data:image/png;base64,', fixed = TRUE)
  expect_match(favicon, 'id="liberation-blue"', fixed = TRUE)
  expect_false(grepl('id="liberties-red"', favicon, fixed = TRUE))
  expect_false(grepl('<circle cx="16"', favicon, fixed = TRUE))
})

test_that("GUI parameter names support models without random effects", {
  expect_equal(
    LibeRation:::.nm_parameter_names(c(2, 20), 0.05, numeric()),
    c("THETA1", "THETA2", "SIGMA1")
  )
  expect_equal(
    LibeRation:::.nm_parameter_names(numeric(), numeric(), numeric()),
    character()
  )
})

test_that("all supported ADVAN templates are valid models", {
  for (advan in c(1L, 2L, 3L, 4L, 6L, 11L, 12L, 13L)) {
    expect_s3_class(LibeRation:::.liber_model_template(advan), "nm_model")
  }
  ode <- LibeRation:::.liber_model_template(13L, n_state = 4L, problem = "Four-state ODE")
  expect_equal(ode$n_state, 4L)
  expect_equal(attr(ode, "name"), "Four-state ODE")
  expect_s3_class(liber_gui(workspace = tempfile("gui-workspace-"), launch.browser = NULL),
                  "shiny.appobj")
})

test_that("GUI page shell fills the viewport without Bootstrap gutters", {
  page <- LibeRation:::.liber_full_page_ui(
    htmltools::tags$head(htmltools::tags$title("LibeRation")),
    htmltools::tags$div(id = "workbench")
  )
  rendered <- htmltools::renderTags(page)
  markup <- paste(rendered$head, rendered$html)

  expect_match(markup, "body > .container-fluid", fixed = TRUE)
  expect_match(markup, "class=\"liberation-app-root\"", fixed = TRUE)
  expect_match(markup, "margin: 0; padding: 0; overflow: hidden;", fixed = TRUE)
})

test_that("workbench content remains scrollable in short viewports", {
  css_path <- system.file(
    "htmlwidgets", "liberWorkbench.css", package = "LibeRation"
  )
  css <- paste(readLines(css_path, warn = FALSE), collapse = "\n")

  expect_match(
    css, "width: 100%; height: 100%; min-height: 0;", fixed = TRUE
  )
  expect_match(
    css,
    ".lw-legacy-shell { min-height: 0; overflow: hidden; }.lw-page-host { overflow: auto; }",
    fixed = TRUE
  )
  expect_false(grepl("min-height: 720px", css, fixed = TRUE))
  expect_false(grepl("min-height: 760px", css, fixed = TRUE))

  dependency_path <- system.file(
    "htmlwidgets", "liberWorkbench.yaml", package = "LibeRation"
  )
  dependency <- paste(readLines(dependency_path, warn = FALSE), collapse = "\n")
  expect_match(dependency, "version: 0.8.0", fixed = TRUE)
})

test_that("hosted GUI sessions receive different ephemeral workspaces", {
  root <- tempfile("gui-hosted-")
  app <- liber_gui(
    workspace = root, queue = FALSE, session_workspace = TRUE,
    launch.browser = NULL
  )
  server <- app[["serverFuncSource"]]()
  paths <- character()
  shiny::testServer(server, {
    paths <<- c(paths, workspace$path)
  })
  shiny::testServer(server, {
    paths <<- c(paths, workspace$path)
  })
  expect_length(unique(paths), 2L)
  expect_true(all(startsWith(
    normalizePath(paths, winslash = "/"),
    paste0(normalizePath(root, winslash = "/", mustWork = FALSE), "/sessions/")
  )))
})

test_that("template project creation produces a linked synthetic dataset and version", {
  root <- tempfile("gui-template-")
  dir.create(root)
  app <- liber_gui(workspace = root, launch.browser = NULL)

  shiny::testServer(app[["serverFuncSource"]](), {
    session$setInputs(liber_workbench_event = list(action = "noop", nonce = 0))
    session$flushReact()
    session$setInputs(liber_workbench_event = list(
      action = "project_create", name = "Synthetic project",
      description = "Created from the legacy-style dialog", mode = "template",
      dataSource = "synthetic", example = "sparse", nSubjects = 3L,
      advan = 2L, trans = 2L, label = "Initial oral model",
      problem = "Oral population PK", nonce = 1
    ))
    session$flushReact()
  })

  workspace <- nm_workspace(root)
  expect_equal(nm_project_list(workspace)$description, "Created from the legacy-style dialog")
  versions <- nm_project_list(workspace, "synthetic-project")
  expect_equal(versions$label, "Initial oral model")
  opened <- nm_project_load(workspace, "synthetic-project")
  expect_equal(opened$model$ADVAN, 2L)
  expect_equal(length(unique(opened$data$ID)), 3L)
})

test_that("GUI result states distinguish validation, queueing, and failure", {
  validation <- structure(list(), class = "liber_gui_validation")
  queued <- structure(list(id = "job-1"), class = "liber_gui_queued")

  expect_equal(LibeRation:::.liber_gui_result(validation)$status, "validated")
  expect_equal(LibeRation:::.liber_gui_result(queued)$status, "queued")
  expect_equal(LibeRation:::.liber_gui_result(queued)$job_id, "job-1")
  expect_equal(LibeRation:::.liber_gui_result(simpleError("failed"))$status, "error")
})

test_that("GUI events create projects and immutable model versions", {
  root <- tempfile("gui-workflow-")
  dir.create(root)
  app <- liber_gui(workspace = root, launch.browser = NULL)

  shiny::testServer(app[["serverFuncSource"]](), {
    send <- function(event) {
      session$setInputs(liber_workbench_event = event)
      session$flushReact()
    }
    # The widget initializes this input in a browser. Prime it explicitly in
    # the server test because the event observer intentionally ignores init.
    send(list(action = "noop", nonce = 0))
    send(list(action = "project_create", name = "Demo", nonce = 1))
    send(list(action = "model_template", advan = 1L, nonce = 2))
    send(list(
      action = "load_csv", name = "demo.csv",
      text = paste(
        "ID,TIME,AMT,DV,EVID,CMT",
        "1,0,100,0,1,1",
        "1,1,0,5,0,1",
        sep = "\n"
      ),
      nonce = 3
    ))
    send(list(action = "project_save", label = "Initial model", nonce = 4))
  })

  workspace <- nm_workspace(root)
  versions <- nm_project_list(workspace, "demo")
  expect_equal(nrow(versions), 2L)
  expect_setequal(versions$label, c("Template model", "Initial model"))
  opened <- nm_project_load(workspace, "demo", versions$id[[1L]])
  expect_s3_class(opened$model, "nm_model")
  expect_s3_class(opened$data, "nm_dataset")
})

test_that("model comparison includes GOF and only diagnostics common to both runs", {
  fixture <- estimation_fixture()
  fit <- nm_est(fixture$model, fixture$data, method = "FOCEI", maxit = 1)
  root <- tempfile("gui-compare-plots-")
  workspace <- nm_workspace(root)
  project <- nm_project_create(workspace, "Comparison plots")
  version <- nm_project_save(workspace, project$id, fixture$model, fixture$data,
                             label = "Base model")
  first <- nm_project_save_run(workspace, project$id, version, fit, "Run A")
  second_fit <- fit
  second_fit$objective <- fit$objective + 1
  second <- nm_project_save_run(workspace, project$id, version, second_fit, "Run B")
  vpc <- nm_vpc(fit, nsim = 2, seed = 91)
  npc <- nm_npc(fit, nsim = 20, seed = 92)
  nm_project_save_diagnostics(workspace, project$id, first, list(vpc = vpc, npc = npc))
  nm_project_save_diagnostics(workspace, project$id, second, list(vpc = vpc))

  app <- liber_gui(workspace = workspace, project = project$id, launch.browser = NULL)
  server_function <- app[["serverFuncSource"]]()
  comparison <- NULL
  closed_result <- NULL
  comparison_open <- NULL
  shiny::testServer(server_function, {
    session$setInputs(liber_workbench_event = list(action = "noop", nonce = 0))
    session$flushReact()
    session$setInputs(liber_workbench_event = list(
      action = "project_compare", runs = c(first, second), nonce = 1
    ))
    session$flushReact()
    comparison <<- shiny::isolate(state$result)
    session$setInputs(liber_workbench_event = list(
      action = "comparison_close", nonce = 2
    ))
    session$flushReact()
    closed_result <<- shiny::isolate(state$result)
    comparison_open <<- shiny::isolate(state$comparison_open)
  })
  expect_s3_class(comparison, "liber_gui_comparison")
  expect_named(comparison$plots, c("gof", "vpc"))
  expect_length(comparison$plots$gof, 2L)
  expect_length(comparison$plots$vpc, 2L)
  expect_null(comparison$plots$npc)
  expect_null(closed_result)
  expect_false(comparison_open)
})

test_that("local queued jobs are durable and opening a result is idempotent", {
  skip_if_not_installed("LibeRties")
  root <- tempfile("gui-local-queue-")
  queue_root <- tempfile("gui-jobs-")
  dir.create(root)
  queue <- LibeRties::ls_local_queue(queue_root, "local", max_workers = 1L)
  app <- liber_gui(workspace = root, queue = queue, launch.browser = NULL)
  server_function <- app[["serverFuncSource"]]()

  shiny::testServer(server_function, {
    send <- function(event) {
      session$setInputs(liber_workbench_event = event)
      session$flushReact()
    }
    send(list(action = "noop", nonce = 0))
    send(list(
      action = "project_create", name = "Queued project", mode = "template",
      dataSource = "synthetic", example = "sparse", nSubjects = 1L,
      advan = 1L, trans = 2L, label = "Queued model", problem = "Queue test",
      nonce = 1
    ))
    send(list(action = "simulate", label = "Queued simulation", seed = 9L,
              replicates = 1L, useDesign = FALSE, nonce = 2))
    listed <- queue$list()
    expect_equal(nrow(listed), 1L)
    expect_true(listed$status[[1L]] %in% c("queued", "running", "completed"))
    expect_equal(nrow(state$jobs), 1L)
    expect_match(paste(output$workbench, collapse = ""), "Queued simulation", fixed = TRUE)
    id <- listed$id[[1L]]
    expect_equal(queue$wait(id, timeout = 30)$status, "completed")
    send(list(action = "jobs_refresh", nonce = 3))
    records <- nm_project_list(nm_workspace(root), "queued-project")
    expect_equal(sum(records$entry_type == "run"), 1L)
    expect_equal(records$queue_job_id[records$entry_type == "run"], id)
    send(list(action = "job_result", id = id, nonce = 4))
    expect_true(nzchar(shiny::isolate(state$run)))
    send(list(action = "job_result", id = id, nonce = 5))
  })

  records <- nm_project_list(nm_workspace(root), "queued-project")
  expect_equal(sum(records$entry_type == "run"), 1L)

  restarted <- liber_gui(workspace = root, queue = queue, launch.browser = NULL)
  restarted_server <- restarted[["serverFuncSource"]]()
  restarted_jobs <- NULL
  shiny::testServer(restarted_server, {
    restarted_jobs <<- shiny::isolate(state$jobs)
  })
  expect_equal(nrow(restarted_jobs), 1L)
  expect_equal(restarted_jobs$status[[1L]], "completed")
})

test_that("default liber_gui creates and exposes its persistent local queue", {
  skip_if_not_installed("LibeRties")
  root <- tempfile("gui-default-queue-")
  dir.create(root)
  app <- liber_gui(workspace = root, queue = NULL, launch.browser = NULL)
  server_function <- app[["serverFuncSource"]]()
  created_queue <- NULL

  shiny::testServer(server_function, {
    created_queue <<- queue
    send <- function(event) {
      session$setInputs(liber_workbench_event = event)
      session$flushReact()
    }
    send(list(action = "noop", nonce = 0))
    send(list(
      action = "project_create", name = "Default queue", mode = "template",
      dataSource = "synthetic", example = "sparse", nSubjects = 1L,
      advan = 1L, trans = 2L, label = "Model", problem = "Queue", nonce = 1
    ))
    send(list(action = "simulate", label = "Visible local job", seed = 5L,
              replicates = 1L, useDesign = FALSE, nonce = 2))
    expect_equal(nrow(shiny::isolate(state$jobs)), 1L)
    expect_equal(nrow(queue$list()), 1L)
    expect_match(paste(output$workbench, collapse = ""), "Visible local job", fixed = TRUE)
  })
  expect_s3_class(created_queue, "LibeRQueue")
  expect_equal(normalizePath(created_queue$root, winslash = "/"),
               normalizePath(file.path(root, ".jobs"), winslash = "/", mustWork = FALSE))
})
