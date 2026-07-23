test_that("ecosystem support matrix uses explicit evidence tiers", {
  support <- liber_support_matrix()
  expect_true(all(c(
    "package", "capability", "status", "evidence_tier", "reference",
    "gate", "last_verified", "recommended_use"
  ) %in% names(support)))
  expect_setequal(
    unique(support$package),
    c("LibeRtAD", "LibeRation", "LibeRties", "LibeRary", "LibeRator", "LibeRality")
  )
  expect_true(all(support$evidence_tier %in% c("validated", "verified", "experimental")))
  expect_true(all(nzchar(support$gate)))
  expect_true(all(nm_support_matrix()$evidence_tier %in%
                    c("validated", "verified", "experimental")))
})

test_that("support bundles omit values and redact opt-in logs", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "DV", "MDV", "NOTE"),
    ADVAN = 1,
    PRED = "CL=THETA(1)*exp(ETA(1));V=THETA(2);S1=V",
    ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(THETA = 1:2, Value = c(2, 20), FIX = TRUE),
    OMEGAS = data.frame(OMEGA = 1, Value = 0.1, FIX = TRUE),
    SIGMAS = data.frame(SIGMA = 1, Value = 0.1, FIX = TRUE)
  )
  data <- data.frame(
    ID = 1, TIME = c(0, 1), EVID = c(1, 0), AMT = c(100, 0),
    DV = c(NA, 3.14159265), MDV = c(1, 0),
    NOTE = c("patient@example.org", "private-value")
  )
  output <- tempfile(fileext = ".zip")
  bundle <- liber_support_bundle(
    output, model = model, data = data,
    logs = c(
      "password=hunter2 token=secret-token patient@example.org",
      "C:\\Users\\Alice\\private\\worker.log"
    ),
    include_logs = TRUE
  )
  expect_true(file.exists(bundle))
  extracted <- tempfile("support-extract-")
  dir.create(extracted)
  utils::unzip(bundle, exdir = extracted)
  content <- paste(unlist(lapply(
    list.files(extracted, recursive = TRUE, full.names = TRUE),
    readLines, warn = FALSE
  )), collapse = "\n")
  expect_false(grepl("hunter2|secret-token|patient@example.org|private-value|3\\.14159265",
                     content))
  expect_false(grepl("Users\\\\Alice", content, fixed = TRUE))
  expect_match(content, "REDACTED_SECRET")
  expect_match(content, "REDACTED_EMAIL")

  manifest <- jsonlite::fromJSON(file.path(extracted, "manifest.json"))
  expect_false(manifest$redaction$data_values_included)
  expect_false(manifest$redaction$parameter_estimates_included)
  expect_false(manifest$redaction$model_code_included)
  expect_equal(manifest$data$rows, 2L)
  expect_equal(manifest$model$ADVAN, 1L)
})
