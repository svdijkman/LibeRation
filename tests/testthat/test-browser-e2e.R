test_that("real browser preserves full-page layout and modal reachability", {
  skip_if_not_installed("shinytest2")
  skip_if(Sys.getenv("LIBER_RUN_BROWSER_TESTS") != "true")
  root <- tempfile("liber-browser-")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  app <- LibeRation::liber_gui(workspace = root, queue = FALSE, launch.browser = NULL)
  driver <- shinytest2::AppDriver$new(
    app, name = "liberation-layout", width = 1366, height = 768,
    load_timeout = 120000, seed = 20260722
  )
  on.exit(driver$stop(), add = TRUE)
  driver$wait_for_idle()
  expect_identical(driver$get_js("document.title"), "LibeRation")
  expect_true(driver$get_js("!!document.querySelector('.liberation-app-root')"))
  expect_false(driver$get_js(
    "document.documentElement.scrollWidth > document.documentElement.clientWidth + 2"
  ))

  driver$run_js(paste(
    "Array.from(document.querySelectorAll('button'))",
    ".find(x => x.title === 'Create a new project').click()"
  ))
  Sys.sleep(0.25)
  expect_true(driver$get_js(
    "Array.from(document.querySelectorAll('.lw-modal-header strong')).some(x => x.textContent === 'New project')"
  ))
  expect_true(driver$get_js(paste(
    "(() => { const modal=document.querySelector('.lw-modal');",
    "if(!modal) return false; const r=modal.getBoundingClientRect();",
    "return r.top >= 0 && r.bottom <= window.innerHeight + 1; })()"
  )))
  driver$run_js("document.querySelector('.lw-modal-close').click()")

  expect_true(driver$get_js(paste(
    "(() => { const panel=document.querySelector('.lw-page-host');",
    "if(!panel) return false; panel.scrollTop=panel.scrollHeight;",
    "return panel.scrollTop + panel.clientHeight >= panel.scrollHeight - 2; })()"
  )))
})

test_that("large dataset remains metadata-only until explicitly requested", {
  model <- LibeRation::nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT"), ADVAN = 1,
    PRED = "CL=THETA(1);V=THETA(2);S1=V",
    THETAS = data.frame(THETA = 1:2, Value = c(2, 20))
  )
  data <- data.frame(
    ID = rep(seq_len(1000), each = 100), TIME = rep(seq(0, 24, length.out = 100), 1000),
    EVID = 0L, AMT = 0
  )
  payload <- LibeRation::liber_workbench(model, data, data_payload = FALSE)$x$tag$attribs
  expect_false(payload$dataset$payload_loaded)
  expect_length(payload$dataset$preview, 0L)
  expect_length(payload$dataset$plot_rows, 0L)
  expect_lt(as.numeric(object.size(payload)), 2e6)
})
