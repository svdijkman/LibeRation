test_that("nm_remote_job_list returns error attribute when server is down", {
  skip_if_not_installed("curl")
  skip_if_not_installed("jsonlite")
  srv_id <- paste0("test_down_", format(Sys.time(), "%s"))
  on.exit(tryCatch(nm_remote_server_remove(srv_id), error = function(e) NULL), add = TRUE)
  nm_remote_server_add(
    name = "Unreachable test",
    base_url = "http://127.0.0.1:59999",
    username = "tester",
    token = "test-token",
    id = srv_id,
    set_default = FALSE
  )
  df <- nm_remote_job_list(srv_id)
  expect_equal(nrow(df), 0L)
  expect_true(nzchar(attr(df, "error") %||% ""))
  expect_match(attr(df, "error"), "unreachable|connect|Could not", ignore.case = TRUE)
})
