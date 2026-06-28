pkgload::load_all("C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation", compile = FALSE, quiet = TRUE)
res <- nm_bench_pilot(advan = 1L, trans = 1L, n_per = 2L, seed = 2024L, method = "FOCEI", run_nonmem = TRUE)
fit <- res$rcpp_fit
ext <- res$nm_ext
cat("Rcpp theta:", fit$theta, "\n")
cat("Rcpp omega:", fit$omega, "\n")
cat("Rcpp sigma:", fit$sigma, "\n")
cat("Rcpp obj:", fit$objective, "\n")
cat("NM ext row:\n")
print(ext$table)
cat("NM parsed omega:", ext$omega, " sigma:", ext$sigma, "\n")
cmp <- res$compare
cat("sigma rel:", cmp$d_sigma, " omega rel:", cmp$d_omega, "\n")
