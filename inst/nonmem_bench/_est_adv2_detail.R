pkgload::load_all("C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation", compile = FALSE, quiet = TRUE)
res <- nm_bench_pilot(advan = 2L, trans = 2L, n_per = 3L, seed = 2024L, method = "FOCEI", run_nonmem = TRUE)
fit <- res$rcpp_fit; ext <- res$nm_ext
cat("rcpp sigma:", fit$sigma, "\n")
cat("nm sigma (SD):", ext$sigma, "\n")
cat("nm table sigma cols:\n")
print(ext$table[, grep("SIGMA", names(ext$table), value = TRUE)])
cat("theta rcpp:", fit$theta, "\n nm:", ext$theta, "\n")
cat("omega rcpp:", fit$omega, "\n nm:", ext$omega, "\n")
