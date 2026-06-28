test_that("nm_job_submit completes FO fit", {
  skip_if_not_installed("callr")
  skip_if_not(
    requireNamespace("LibeRation", quietly = TRUE),
    "LibeRation must be installed for background jobs"
  )
  sim <- nm_synthetic_theo(n_sub = 3L, seed = 99L)
  root <- file.path(tempdir(), "LibeRation_jobs_test", format(Sys.time(), "%s"))
  on.exit(unlink(dirname(root), recursive = TRUE), add = TRUE)
  old_env <- getOption("LibeRation.job_dev_env")
  on.exit(
    if (is.null(old_env)) {
      options(LibeRation.job_dev_env = NULL)
    } else {
      options(LibeRation.job_dev_env = old_env)
    },
    add = TRUE
  )
  options(LibeRation.job_dev_env = list(mode = "installed", nm_root = "", ad_root = ""))
  job <- nm_job_submit(
    sim$model, sim$data,
    method = "FO",
    grad = "numeric",
    pk_engine = "cpp",
    control = list(maxit = 8L),
    job_root = root,
    label = "test FO"
  )
  expect_s3_class(job, "nm_job_handle")
  job$process$wait()
  st <- nm_job_status(job$id, root)
  expect_equal(st$status, "success")
  fit <- nm_job_result(job$id, root)
  expect_s3_class(fit, "nm_fit")
  expect_equal(fit$method, "FO")
  expect_true(is.finite(fit$objective))
})

test_that("nm_job_list and cleanup work", {
  skip_if_not_installed("callr")
  root <- file.path(tempdir(), "LibeRation_jobs_list", format(Sys.time(), "%s"))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  df <- nm_job_list(root)
  expect_equal(nrow(df), 0L)
})

test_that("nm_job_cancel stops a running job", {
  skip_if_not_installed("callr")
  sim <- nm_synthetic_theo(n_sub = 2L, seed = 7L)
  root <- file.path(tempdir(), "LibeRation_jobs_cancel", format(Sys.time(), "%s"))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  job <- nm_job_submit(
    sim$model, sim$data,
    method = "FOCE",
    grad = "numeric",
    pk_engine = "cpp",
    control = list(maxit = 200L),
    max_outer = 50L,
    job_root = root,
    label = "cancel me"
  )
  st <- nm_job_cancel(job$id, root)
  expect_equal(st$status, "cancelled")
})
