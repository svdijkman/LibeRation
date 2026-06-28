pkgload::load_all("C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation", compile = FALSE, quiet = TRUE)
wd <- tempfile("adv1_est_"); dir.create(wd)
res <- nm_bench_pilot(advan = 1L, trans = 1L, n_per = 2L, seed = 2024L,
                      method = "FOCEI", run_nonmem = TRUE, work_dir = wd)
fit <- res$rcpp_fit
cat("convergence:", fit$convergence, " outer len:", length(fit$outer), "\n")
cat("outer vals:", unlist(fit$outer), "\n")
cmp <- res$compare
cat("theta rcpp:", fit$theta, " nm:", cmp$theta_nm, "\n")
cat("omega rcpp:", fit$omega, " nm:", res$nm_ext$omega, "\n")
cat("sigma rcpp:", fit$sigma, " nm sd:", res$nm_ext$sigma, "\n")
# obj at NM params
imp <- nm_read_nonmem(res$est_ctl_path, data_path = file.path(wd, "data.csv"))
imp$model$INPUT <- res$parts$input_cols
ds <- structure(list(data = res$data), class = "nm_dataset")
phi <- read.table(file.path(wd, "est.phi"), skip = 2, header = FALSE)
eta_nm <- matrix(as.numeric(phi[[3]]), ncol = 1)
obj_at_nm <- .nm_focei_objective(
  imp$model, ds, res$nm_ext$theta, res$nm_ext$omega, res$nm_ext$sigma,
  eta_nm, pk_engine = "cpp"
)
cat("obj at NM par:", obj_at_nm, " NM:", res$nm_ext$obj, "\n")
cat("obj at rcpp par:", fit$objective, "\n")
