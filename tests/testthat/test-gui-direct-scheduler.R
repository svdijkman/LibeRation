test_that("direct scheduler settings support a jump host without exposing shell input", {
  key <- paste(rep("ab", 32L), collapse = "")
  config <- .liber_direct_scheduler_normalize(list(
    backend = "slurm", queue_name = "teaching", storage_key = key,
    remote_rscript = "/opt/R/bin/Rscript", max_workers = 4L,
    max_cores_per_job = 16L, max_memory_mb = 8192L,
    partition = "compute", account = "group-a",
    ssh = list(
      host = "login.example.org", user = "student", port = 22L,
      proxy_host = "gateway.example.org", proxy_user = "student",
      proxy_port = 22L, identity_file = "", accept_new_host_key = TRUE
    )
  ), check_identity = FALSE)
  expect_identical(config$backend, "slurm")
  expect_identical(config$ssh$proxy_host, "gateway.example.org")
  expect_identical(config$limits$max_concurrent_jobs, 4L)
  expect_identical(.liber_direct_scheduler_normalize(
    config, check_identity = FALSE
  )$limits$max_memory_mb, 8192L)
  command <- .liber_ssh_scheduler_command(config, ssh = "ssh")
  expect_true("-J" %in% command$args)
  expect_true(any(command$args == "student@gateway.example.org"))
  expect_true(grepl("ls_direct_scheduler_cli", utils::tail(command$args, 1L), fixed = TRUE))
  expect_false(any(grepl(key, command$args, fixed = TRUE)))
})

test_that("direct scheduler response verification fails closed", {
  inner <- jsonlite::toJSON(
    list(ok = TRUE, payload = list(username = "local")),
    auto_unbox = TRUE, null = "null", force = TRUE
  )
  envelope <- jsonlite::toJSON(list(
    schema = "liberties.ssh.response", version = 1L,
    payload_json = as.character(inner), sha256 = paste(rep("0", 64L), collapse = "")
  ), auto_unbox = TRUE, null = "null", force = TRUE)
  expect_error(.liber_ssh_scheduler_response(paste(
    "LIBERTIES_SSH_RESPONSE_BEGIN", envelope, "LIBERTIES_SSH_RESPONSE_END", sep = "\n"
  )), "checksum mismatch")
})

test_that("saved direct scheduler secrets are not rendered into queue metadata", {
  public <- .liber_direct_scheduler_public_config(list(
    backend = "slurm", queue_name = "teaching",
    storage_key = paste(rep("ab", 32L), collapse = "")
  ))
  expect_identical(public$backend, "slurm")
  expect_identical(public$queue_name, "teaching")
  expect_false("storage_key" %in% names(public))
})

test_that("direct scheduler definitions survive client-settings upgrades", {
  workspace <- tempfile("direct-scheduler-settings-")
  dir.create(workspace, recursive = TRUE)
  key <- paste(rep("12", 32L), collapse = "")
  remotes <- list(cluster = list(
    name = "Teaching cluster", connection_mode = "ssh_scheduler", timeout = 30,
    scheduler = list(
      backend = "slurm", queue_name = "teaching", storage_key = key,
      remote_rscript = "Rscript", max_workers = 4L, max_cores_per_job = 16L,
      ssh = list(host = "login.example.org", user = "student", port = 22L)
    )
  ))
  .liber_client_settings_write(workspace, selected_queue = "cluster", remotes = remotes)
  restored <- .liber_client_settings_read(workspace)
  expect_identical(restored$version, 9L)
  expect_identical(restored$selected_queue, "cluster")
  expect_identical(restored$remotes$cluster$scheduler$storage_key, key)
})
