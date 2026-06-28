pkgload::load_all("C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation", compile = FALSE)
p <- nm_ctl_template(2L, 2L)
d <- nm_bench_mixed_design(2L, 2L, n_per = 1, seed = 2024L)
r <- nm_bench_simulate_rcppnm(p, d, seed = 2024L)
cat("adv2_trans2 rcpp max IPRED:", max(r$data$IPRED, na.rm = TRUE), "\n")

p12 <- nm_ctl_template(12L, 4L)
d12 <- nm_bench_mixed_design(12L, 4L, n_per = 1, seed = 2024L)
r12 <- nm_bench_simulate_rcppnm(p12, d12, seed = 2024L)
cat("adv12_trans4 rcpp max IPRED:", max(r12$data$IPRED, na.rm = TRUE), "\n")

pk <- .nm_bench_pk_nonmem_lines(p12)
cat("adv12 pk extras:\n", paste(pk, collapse = "\n"), "\n")
