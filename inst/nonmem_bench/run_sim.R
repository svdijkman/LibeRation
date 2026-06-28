#!/usr/bin/env Rscript
# Run structural simulation benchmarks only (no estimation).
# Usage: Rscript inst/nonmem_bench/run_sim.R [output_dir]

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) {
  normalizePath(args[[1L]], winslash = "/", mustWork = FALSE)
} else {
  NA_character_
}

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg)) {
  pkg_root <- normalizePath(
    file.path(dirname(sub("^--file=", "", file_arg)), "..", ".."),
    winslash = "/",
    mustWork = FALSE
  )
} else {
  pkg_root <- getwd()
}

if (is.na(out_dir)) {
  out_dir <- file.path(pkg_root, "inst", "nonmem_bench", "sim_run")
}

if (requireNamespace("devtools", quietly = TRUE) &&
    file.exists(file.path(pkg_root, "DESCRIPTION"))) {
  devtools::load_all(pkg_root, quiet = TRUE, compile = TRUE)
}

if (!requireNamespace("LibeRation", quietly = TRUE)) {
  stop("Load LibeRation first, e.g. devtools::load_all('.')")
}

if (!LibeRation::nm_nonmem_available()) {
  stop("NONMEM not found on PATH (need nmfe73).")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("Simulation-only benchmark output: ", out_dir)

res <- LibeRation::nm_bench_run_all(
  work_dir = out_dir,
  n_per = 2L,
  seed = 2024L,
  run_nonmem = TRUE,
  run_est = FALSE,
  sim_ipred_rtol = 0.03
)

summary_path <- file.path(out_dir, "summary.csv")
sim_summary_path <- file.path(out_dir, "sim_summary.csv")
utils::write.csv(res$summary, summary_path, row.names = FALSE)
utils::write.csv(res$sim_summary, sim_summary_path, row.names = FALSE)

cat("\n=== Simulation summary ===\n")
print(res$sim_summary[, c("tag", "sim_ok", "max_ipred_rel", "max_amt_rel", "message")])

sim_reg_cols <- c(
  "tag",
  "sim_ok_single", "max_ipred_rel_single",
  "sim_ok_multiple", "max_ipred_rel_multiple",
  "sim_ok_steady_state", "max_ipred_rel_steady_state"
)
if (all(sim_reg_cols %in% names(res$sim_summary))) {
  cat("\n=== Simulation by regimen ===\n")
  print(res$sim_summary[, sim_reg_cols])
}

n_sim_pass <- sum(res$sim_summary$sim_ok, na.rm = TRUE)
n_sim_total <- sum(!is.na(res$sim_summary$sim_ok))

cat(sprintf("\nSim pass: %d / %d\n", n_sim_pass, n_sim_total))
cat("Summaries written to:\n")
cat(" ", summary_path, "\n")
cat(" ", sim_summary_path, "\n")

invisible(res)
