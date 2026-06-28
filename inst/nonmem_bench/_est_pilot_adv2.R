pkgload::load_all("C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation", compile = FALSE, quiet = TRUE)
if (!nm_nonmem_available()) stop("no nm")
res <- nm_bench_pilot(advan = 2L, trans = 2L, n_per = 3L, seed = 2024L, method = "FOCEI", run_nonmem = TRUE)
cmp <- res$compare
cat("sim_ok:", res$status$sim_ok, " est_ok:", res$status$est_ok, "\n")
if (!is.null(cmp)) {
  cat("obj rel:", abs(cmp$obj_rcpp - cmp$obj_nm) / max(abs(cmp$obj_nm), 1e-8), "\n")
  cat("theta rel max:", max(cmp$d_theta), " omega max:", max(cmp$d_omega), " sigma max:", max(cmp$d_sigma), "\n")
  cat("ok:", cmp$ok, "\n")
}
