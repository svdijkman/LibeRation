#!/usr/bin/env Rscript
# Pilot NONMEM vs LibeRation benchmark (ADVAN 2 TRANS 2).
# Usage: Rscript bench/nonmem_bench/run_pilot.R

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg)) {
  pkg_root <- normalizePath(
    file.path(dirname(sub("^--file=", "", file_arg)), "..", ".."),
    winslash = "/",
    mustWork = FALSE
  )
} else {
  pkg_root <- getwd()
}
if (file.exists(file.path(pkg_root, "DESCRIPTION")) &&
    requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(pkg_root, quiet = TRUE)
}

if (!requireNamespace("LibeRation", quietly = TRUE)) {
  stop("Load LibeRation first, e.g. devtools::load_all('.')")
}

if (!LibeRation::nm_nonmem_available()) {
  stop("NONMEM not found on PATH (need nmfe73).")
}

res <- LibeRation::nm_bench_pilot(
  advan = 2L,
  trans = 2L,
  n_per = 2L,
  seed = 2024L,
  method = "FOCEI",
  run_nonmem = TRUE
)

cat("\nPilot benchmark complete.\n")
cat("Simulation OK:", res$sim_cmp$ok, "\n")
if (!is.null(res$est_cmp)) {
  cat("Estimation OK:", res$est_cmp$ok, "\n")
}

invisible(res)
