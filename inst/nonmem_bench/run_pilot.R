#!/usr/bin/env Rscript
# Pilot NONMEM vs LibeRation benchmark (ADVAN 2 TRANS 2).
# Usage: Rscript inst/nonmem_bench/run_pilot.R

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
if (dir.exists(file.path(pkg_root, "DESCRIPTION")) &&
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
  n_per = 3L,
  seed = 2024L,
  method = "FOCEI",
  run_nonmem = TRUE
)

cat("Work directory:", res$work_dir, "\n")
cat("LibeRation OBJ:", res$rcpp_fit$objective, "\n")
if (!is.null(res$nm_run) && isFALSE(res$nm_run$license_ok)) {
  cat("NONMEM: license expired or invalid (no .ext output)\n")
} else {
  cat("NONMEM OBJ:", res$nm_ext$obj, "\n")
  cat("Compare OK:", res$compare$ok, "\n")
  print(res$compare)
}
