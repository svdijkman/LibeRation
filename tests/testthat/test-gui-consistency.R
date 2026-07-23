test_that("modelling GUI retains shared theme and accessible dialogs", {
  script <- paste(readLines(
    system.file("htmlwidgets", "liberWorkbench.js", package = "LibeRation"),
    warn = FALSE
  ), collapse = "\n")
  css <- paste(readLines(
    system.file("htmlwidgets", "liberWorkbench.css", package = "LibeRation"),
    warn = FALSE
  ), collapse = "\n")

  expect_match(script, 'localStorage\\.getItem\\("liber\\.theme"\\)')
  expect_match(script, "useDialogFocus", fixed = TRUE)
  expect_match(script, 'event\\.key === "Escape"')
  expect_match(script, '"aria-label": props.title', fixed = TRUE)
  expect_match(css, "focus-visible", fixed = TRUE)
  expect_match(css, "grid-template-rows: 58px 42px 32px", fixed = TRUE)
  expect_match(css, ".lw-app-icon { width: 42px; height: 42px", fixed = TRUE)
  expect_match(css, ".lw-button { min-height: 32px", fixed = TRUE)
  expect_match(css, "border-radius: 10px", fixed = TRUE)
})
