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
  expect_match(source, 'Print gradients every N (0 = off)', fixed = TRUE)
  expect_match(source, 'Open the saved model run and its results', fixed = TRUE)
  expect_match(source, 'comparison_close', fixed = TRUE)
  expect_match(source, 'page_change', fixed = TRUE)
  expect_match(source, 'lw-app-icon', fixed = TRUE)
  expect_match(source, 'lw-app-version', fixed = TRUE)
  expect_match(source, 'Estimation priors', fixed = TRUE)
  expect_match(source, 'lw-modal-comparison', fixed = TRUE)
  expect_match(source, 'lw-run-flag-cov', fixed = TRUE)
  expect_match(source, 'load_payload', fixed = TRUE)
  expect_match(source, 'Stratify by', fixed = TRUE)
  expect_match(source, 'function ComparisonPlots', fixed = TRUE)
  expect_match(source, 'lw-button-danger-ghost', fixed = TRUE)
  server_source <- paste(deparse(body(liber_gui)), collapse = "\n")
  expect_match(server_source, 'q\\$logs\\(selected,\\s+stream = "stdout"', perl = TRUE)
  expect_match(server_source, 'q\\$logs\\(id,\\s+stream = "stdout"', perl = TRUE)
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
  state <- get("state", envir = environment(server_function), inherits = TRUE)
  favicon_href <- get("favicon_href", envir = environment(server_function), inherits = TRUE)
  expect_equal(shiny::isolate(state$queue_id), "team")
  expect_equal(shiny::isolate(state$remote_meta)$team$name, "Team server")
  expect_s3_class(shiny::isolate(state$remote_queues)$team, "LibeRRemote")
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
  state <- get("state", envir = environment(server_function), inherits = TRUE)
  comparison <- NULL
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
  })
  expect_s3_class(comparison, "liber_gui_comparison")
  expect_named(comparison$plots, c("gof", "vpc"))
  expect_length(comparison$plots$gof, 2L)
  expect_length(comparison$plots$vpc, 2L)
  expect_null(comparison$plots$npc)
  expect_null(shiny::isolate(state$result))
  expect_false(shiny::isolate(state$comparison_open))
})

test_that("local queued jobs are durable and opening a result is idempotent", {
  skip_if_not_installed("LibeRties")
  root <- tempfile("gui-local-queue-")
  queue_root <- tempfile("gui-jobs-")
  dir.create(root)
  queue <- LibeRties::ls_local_queue(queue_root, "local", max_workers = 1L)
  app <- liber_gui(workspace = root, queue = queue, launch.browser = NULL)
  server_function <- app[["serverFuncSource"]]()
  state <- get("state", envir = environment(server_function), inherits = TRUE)

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
  restarted_state <- get("state", envir = environment(restarted_server), inherits = TRUE)
  restarted_jobs <- shiny::isolate(restarted_state$jobs)
  expect_equal(nrow(restarted_jobs), 1L)
  expect_equal(restarted_jobs$status[[1L]], "completed")
})

test_that("default liber_gui creates and exposes its persistent local queue", {
  skip_if_not_installed("LibeRties")
  root <- tempfile("gui-default-queue-")
  dir.create(root)
  app <- liber_gui(workspace = root, queue = NULL, launch.browser = NULL)
  server_function <- app[["serverFuncSource"]]()
  state <- get("state", envir = environment(server_function), inherits = TRUE)
  queue <- get("queue", envir = environment(server_function), inherits = TRUE)
  expect_s3_class(queue, "LibeRQueue")
  expect_equal(normalizePath(queue$root, winslash = "/"),
               normalizePath(file.path(root, ".jobs"), winslash = "/", mustWork = FALSE))

  shiny::testServer(server_function, {
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
})
