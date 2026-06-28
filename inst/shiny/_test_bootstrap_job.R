suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
if (!requireNamespace("callr", quietly = TRUE)) {
  stop("callr not installed")
}
sim <- nm_synthetic_theo(n_sub = 3L, seed = 99L)
root <- file.path(tempdir(), "LibeRation_boot_job", format(Sys.time(), "%s"))
job <- nm_job_submit(
  sim$model, sim$data,
  method = "FO",
  grad = "numeric",
  pk_engine = "cpp",
  control = list(maxit = 8L),
  bootstrap_n = 3L,
  bootstrap_seed = 1L,
  bootstrap_control = list(maxit = 5L, n_cores = 1L, compute_inference = FALSE),
  job_root = root,
  label = "boot test"
)
deadline <- Sys.time() + 180
repeat {
  st <- nm_job_status(job$id, root)
  if (st$status %in% c("success", "error", "cancelled")) break
  Sys.sleep(0.5)
}
cat("status:", st$status, "\n")
if (st$status == "success") {
  fit <- nm_job_result(job$id, root)
  cat("bootstrap n_ok:", fit$bootstrap$n_ok, "n_boot:", fit$bootstrap$n_boot, "\n")
  cat("se names:", paste(names(fit$bootstrap$se), collapse = ", "), "\n")
} else {
  cat("error:", st$error, "\n")
  log_path <- file.path(root, job$id, "worker.log")
  if (file.exists(log_path)) cat(readLines(log_path), sep = "\n")
}
