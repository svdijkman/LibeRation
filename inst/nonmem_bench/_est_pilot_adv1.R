pkgload::load_all("C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation", compile = FALSE, quiet = TRUE)
if (!nm_nonmem_available()) stop("NONMEM not available")
res <- nm_bench_pilot(advan = 1L, trans = 1L, n_per = 2L, seed = 2024L,
                      method = "FOCEI", run_nonmem = TRUE)
cmp <- res$compare
cat("est_ok:", res$status$est_ok, " outer:", length(res$rcpp_fit$outer), "\n")
cat("obj rel:", abs(cmp$obj_rcpp - cmp$obj_nm) / cmp$obj_nm, "\n")
cat("theta rel:", cmp$d_theta, " omega:", cmp$d_omega, " sigma:", cmp$d_sigma, "\n")
cat("ok:", cmp$ok, "\n")
