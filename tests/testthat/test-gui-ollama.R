test_that("Ollama URLs are restricted to loopback interfaces", {
  expect_equal(
    LibeRation:::.liber_ollama_url("http://localhost:11434/"),
    "http://localhost:11434"
  )
  expect_equal(
    LibeRation:::.liber_ollama_url("http://127.0.0.1:11434"),
    "http://127.0.0.1:11434"
  )
  expect_error(
    LibeRation:::.liber_ollama_url("http://192.168.1.10:11434"),
    "loopback|localhost"
  )
  expect_error(
    LibeRation:::.liber_ollama_url("https://example.com"),
    "loopback|localhost"
  )
})

test_that("Ollama is disabled for hosted and session-isolated applications", {
  expect_true(LibeRation:::.liber_ollama_runtime_allowed(
    host = "127.0.0.1", session_workspace = FALSE
  ))
  expect_false(LibeRation:::.liber_ollama_runtime_allowed(
    host = "0.0.0.0", session_workspace = FALSE
  ))
  expect_false(LibeRation:::.liber_ollama_runtime_allowed(
    allow_ollama = TRUE, host = "0.0.0.0", session_workspace = FALSE
  ))
  expect_false(LibeRation:::.liber_ollama_runtime_allowed(
    allow_ollama = TRUE, session_workspace = TRUE, host = "127.0.0.1"
  ))
  expect_false(LibeRation:::.liber_ai_backend_setting(
    "ollama", ollama_allowed = FALSE
  ) == "ollama")
})

test_that("Ollama discovery normalizes the installed model catalogue", {
  fake_models <- function(base_url) {
    expect_equal(base_url, "http://127.0.0.1:11434")
    data.frame(
      id = c("qwen3:8b", "qwen3:8b", "qwen2.5:3b"),
      size = c(5 * 1024^3, 5 * 1024^3, 2 * 1024^3)
    )
  }
  models <- LibeRation:::.liber_ollama_models(
    "http://127.0.0.1:11434", models = fake_models
  )
  expect_equal(vapply(models, `[[`, character(1), "id"),
               c("qwen3:8b", "qwen2.5:3b"))
  expect_match(models[[1L]]$label, "5.0 GB", fixed = TRUE)
  expect_equal(LibeRation:::.liber_ollama_default_model(models), "qwen3:8b")
  expect_equal(
    LibeRation:::.liber_ollama_default_model(list(
      list(id = "qwen3.6:27b"), list(id = "qwen3.5:9b"),
      list(id = "qwen2.5:7b-instruct")
    )),
    "qwen3.5:9b"
  )
})

test_that("Ollama messages preserve system prompt and conversation turns", {
  spec <- LibeRation:::.liber_ollama_message_spec(list(
    list(role = "system", content = "Use supplied evidence."),
    list(role = "user", content = "What is FOCEI?"),
    list(role = "assistant", content = "An estimation method."),
    list(role = "user", content = "Explain its fit.")
  ))
  expect_equal(spec$system_prompt, "Use supplied evidence.")
  expect_equal(spec$prompt, "Explain its fit.")
  expect_equal(vapply(spec$history, `[[`, character(1), "role"),
               c("user", "assistant"))
})

test_that("AI settings persist independent WebLLM and Ollama selections", {
  workspace <- nm_workspace(tempfile("ollama-settings-"))
  LibeRation:::.liber_client_settings_write(
    workspace,
    ai = list(
      activated = TRUE, consented = TRUE, backend = "ollama",
      help_model = LibeRation:::.liber_ai_default_help_model(),
      report_model = "same_as_help", help_context = "8192",
      report_context = "12288",
      ollama_url = "http://localhost:11434",
      ollama_help_model = "qwen3:8b",
      ollama_report_model = "qwen3:14b"
    )
  )
  restored <- LibeRation:::.liber_client_settings_read(workspace)
  expect_equal(restored$version, 6L)
  expect_equal(restored$ai$backend, "ollama")
  expect_equal(restored$ai$ollama_url, "http://localhost:11434")
  expect_equal(restored$ai$ollama_help_model, "qwen3:8b")
  expect_equal(restored$ai$ollama_report_model, "qwen3:14b")
})

test_that("hosted sessions reject a crafted Ollama settings event", {
  root <- tempfile("hosted-ollama-")
  app <- liber_gui(
    workspace = root, queue = FALSE, session_workspace = TRUE,
    allow_ollama = FALSE, launch.browser = NULL
  )
  server <- app[["serverFuncSource"]]()
  shiny::testServer(server, {
    session$setInputs(liber_workbench_event = list(
      action = "ai_settings", backend = "ollama",
      activated = TRUE, consented = TRUE,
      ollama_url = "http://127.0.0.1:11434",
      ollama_help_model = "qwen3:8b",
      ollama_report_model = "same_as_help",
      nonce = 1
    ))
    session$flushReact()
    expect_equal(state$ai_config$backend, "webllm")
  })
})

test_that("hosted deployment explicitly disables Ollama", {
  launcher <- file.path(
    testthat::test_path(), "..", "..", "..", "deploy", "shinyapps",
    "liberation", "app.R"
  )
  skip_if_not(
    file.exists(launcher),
    "The shinyapps launcher is validated by the monorepo CI preflight."
  )
  app <- readLines(launcher, warn = FALSE)
  expect_match(paste(app, collapse = "\n"), "allow_ollama = FALSE", fixed = TRUE)
})
