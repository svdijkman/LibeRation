test_that("GUI persists lazy AI settings and applies nonlinear diagrams", {
  model <- LibeRation:::.liber_model_template(
    6L, n_state = 2L, problem = "GUI visual model"
  )
  data <- nm_dataset(data.frame(
    ID = c(1, 1), TIME = c(0, 1), EVID = c(1, 0), AMT = c(100, 0),
    RATE = 0, II = 0, SS = 0, CMT = 1, DV = c(NA, 4), MDV = c(1, 0)
  ))
  workspace <- nm_workspace(tempfile("visual-gui-"))
  project <- nm_project_create(workspace, "Visual GUI")
  nm_project_save(workspace, project$id, model, data, label = "Base")
  app <- liber_gui(workspace = workspace, project = project$id,
                   launch.browser = NULL)
  server <- app[["serverFuncSource"]]()
  graph <- list(
    title = "Nonlinear diagram", advan = 6L, residual = "additive",
    covariates = list(),
    compartments = list(list(
      id = 1L, name = "CENTRAL", kind = "amount",
      volume_parameter = "V", scale_parameter = "V",
      dose = TRUE, observe = TRUE, x = 250, y = 180
    )),
    flows = list(list(
      id = "elim", from = 1L, to = 0L, type = "michaelis_menten",
      parameter = "VMAX", secondary_parameter = "KM",
      expression = "", label = ""
    )), parameters = list()
  )
  graph$parameters <- Map(function(name, initial) list(
    name = name, initial = initial, lower = NULL, upper = NULL,
    fixed = FALSE, iiv = TRUE, eta_variance = 0.1
  ), c("V", "VMAX", "KM"), c(20, 4, 2))

  shiny::testServer(server, {
    session$setInputs(liber_workbench_event = list(action = "noop", nonce = 0))
    session$flushReact()
    session$setInputs(liber_workbench_event = list(
      action = "ai_settings", activated = TRUE, consented = TRUE,
      help_model = "SmolLM2-360M-Instruct-q4f16_1-MLC",
      report_model = "same_as_help", help_context = "8192",
      report_context = "6144", nonce = 1
    ))
    session$flushReact()
    expect_true(state$ai_config$activated)
    expect_true(state$ai_config$consented)

    session$setInputs(liber_workbench_event = list(
      action = "diagram_preview", graph = graph, nonce = 2
    ))
    session$flushReact()
    expect_s3_class(state$result, "liber_gui_diagram_preview")
    session$setInputs(liber_workbench_event = list(
      action = "diagram_apply", graph = graph, nonce = 3
    ))
    session$flushReact()
    expect_s3_class(state$model$GRAPH, "nm_model_diagram")
    expect_match(state$model$DES, "VMAX")
  })

  settings <- LibeRation:::.liber_client_settings_read(workspace)
  expect_true(settings$ai$activated)
  expect_equal(settings$ai$help_model, "SmolLM2-360M-Instruct-q4f16_1-MLC")
  expect_equal(settings$ai$report_model, "same_as_help")
  expect_equal(settings$ai$help_context, "8192")
  expect_equal(settings$ai$report_context, "6144")
  expect_equal(settings$ai$model, settings$ai$help_model)
})

test_that("browser-local AI exposes distinct quality and memory choices", {
  models <- LibeRation:::.liber_ai_models()
  ids <- vapply(models, `[[`, character(1), "id")
  expect_equal(
    LibeRation:::.liber_ai_default_help_model(),
    "Qwen2.5-Coder-3B-Instruct-q4f16_1-MLC"
  )
  expect_equal(
    LibeRation:::.liber_ai_default_report_model(),
    "Qwen2.5-7B-Instruct-q4f16_1-MLC"
  )
  expect_true(LibeRation:::.liber_ai_default_help_model() %in% ids)
  expect_true(LibeRation:::.liber_ai_default_report_model() %in% ids)
  expect_equal(length(ids), length(unique(ids)))
  expect_true("Qwen2.5-7B-Instruct-q4f16_1-MLC" %in% ids)
  expect_true("Qwen2.5-Coder-3B-Instruct-q4f16_1-MLC" %in% ids)
  expect_true(all(vapply(models, function(model) {
    nzchar(model$label) && nzchar(model$description) && model$vram_mb > 0
  }, logical(1))))
})

test_that("AI context settings are bounded and migration-safe", {
  expect_equal(LibeRation:::.liber_ai_context_setting("auto"), "auto")
  expect_equal(LibeRation:::.liber_ai_context_setting(8192), "8192")
  expect_equal(LibeRation:::.liber_ai_context_setting(10000), "10000")
  expect_equal(LibeRation:::.liber_ai_context_setting(512), "auto")
  expect_equal(LibeRation:::.liber_ai_context_setting(32768), "auto")
})

test_that("single-model AI settings migrate without losing the old choice", {
  workspace <- nm_workspace(tempfile("legacy-ai-settings-"))
  path <- LibeRation:::.liber_client_settings_path(workspace)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(
    version = 3L, selected_queue = "local", remotes = list(),
    pending_jobs = list(), ai = list(
      activated = TRUE, consented = TRUE,
      model = "Qwen2.5-1.5B-Instruct-q4f16_1-MLC"
    )
  ), path)
  settings <- LibeRation:::.liber_client_settings_read(workspace)
  expect_equal(settings$version, 5L)
  expect_equal(
    settings$ai$help_model,
    "Qwen2.5-1.5B-Instruct-q4f16_1-MLC"
  )
  expect_equal(
    settings$ai$report_model,
    "Qwen2.5-7B-Instruct-q4f16_1-MLC"
  )
  expect_equal(settings$ai$help_context, "auto")
  expect_equal(settings$ai$report_context, "auto")
})

test_that("Help and Report generation route through separate lazy models", {
  script <- readLines(
    system.file("htmlwidgets", "liberWorkbench.js", package = "LibeRation"),
    warn = FALSE
  )
  script <- paste(script, collapse = "\n")
  expect_match(script, 'purpose:"help"', fixed = TRUE)
  expect_match(script, 'purpose:"report"', fixed = TRUE)
  expect_match(script, "Only one model is held in GPU memory", fixed = TRUE)
  expect_match(script, "same_as_help", fixed = TRUE)
  expect_match(script, "Local AI settings", fixed = TRUE)
  expect_match(script, "lw-ai-settings-button", fixed = TRUE)
  expect_false(grepl('className:"lw-ai-model-header"', script, fixed = TRUE))
})

test_that("visual builder renders safe structural-parameter deletion", {
  script <- readLines(
    system.file("htmlwidgets", "liberWorkbench.js", package = "LibeRation"),
    warn = FALSE
  )
  script <- paste(script, collapse = "\n")
  expect_match(script, "function removeParameter\\(index\\)")
  expect_match(script, "This parameter is used by a compartment or flow", fixed = TRUE)
  expect_match(script, "lw-diagram-param-remove", fixed = TRUE)
  expect_match(script, "function diagramRenameParameter", fixed = TRUE)
  expect_match(script, "function renumberCompartment", fixed = TRUE)
  expect_match(script, 'label:"Compartment number"', fixed = TRUE)
})

test_that("code editor exposes resizable synchronized parameter tables", {
  script <- paste(readLines(
    system.file("htmlwidgets", "liberWorkbench.js", package = "LibeRation"),
    warn = FALSE
  ), collapse = "\n")
  expect_match(script, "function maximumCodeReference", fixed = TRUE)
  expect_match(script, "function synchronizedParameters", fixed = TRUE)
  expect_match(script, 'unitLabel:"THETA"', fixed = TRUE)
  expect_match(script, 'unitLabel:"ETA"', fixed = TRUE)
  expect_match(script, 'unitLabel:"SIGMA"', fixed = TRUE)
})

test_that("GUI report workflow renders without a selected estimation", {
  model <- LibeRation:::.liber_model_template(6L, n_state = 1L)
  data <- nm_dataset(data.frame(
    ID = c(1, 1), TIME = c(0, 1), EVID = c(1, 0), AMT = c(100, 0),
    RATE = 0, II = 0, SS = 0, CMT = 1, DV = c(NA, 4), MDV = c(1, 0)
  ))
  workspace <- nm_workspace(tempfile("report-gui-"))
  project <- nm_project_create(workspace, "Report GUI")
  nm_project_save(workspace, project$id, model, data)
  app <- liber_gui(workspace = workspace, project = project$id,
                   launch.browser = NULL)
  server <- app[["serverFuncSource"]]()
  blocks <- list(
    list(id = "intro", type = "introduction", title = "Introduction",
         source = "user", text = "A report.", run_ids = list(),
         elements = list(), options = list()),
    list(id = "conclusion", type = "conclusion", title = "Conclusion",
         source = "user", text = "Complete.", run_ids = list(),
         elements = list(), options = list())
  )
  shiny::testServer(server, {
    session$setInputs(liber_workbench_event = list(action = "noop", nonce = 0))
    session$flushReact()
    output_directory <- file.path(tempdir(), "liber-gui-selected-report-directory")
    session$setInputs(liber_workbench_event = list(
      action = "report_design_render", id = "gui-report", title = "GUI report",
      name = "gui-report", directory = output_directory,
      formats = list("pdf"), blocks = blocks, nonce = 1
    ))
    session$flushReact()
    expect_s3_class(state$report, "nm_report_bundle")
    expect_s3_class(state$report_design, "nm_report_design")
    expect_true(file.exists(state$report$pdf))
    expect_equal(normalizePath(dirname(state$report$pdf), winslash = "/"),
                 normalizePath(output_directory, winslash = "/"))
    expect_equal(state$report_design$style$output_directory, output_directory)
  })

  reopened <- liber_gui(workspace = workspace, project = project$id,
                        launch.browser = NULL)
  restored_design <- NULL
  shiny::testServer(reopened[["serverFuncSource"]](), {
    restored_design <<- shiny::isolate(state$report_design)
  })
  expect_s3_class(restored_design, "nm_report_design")
  expect_equal(restored_design$title, "GUI report")
})
