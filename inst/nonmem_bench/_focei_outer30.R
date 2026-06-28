pkgload::load_all("C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation", compile = FALSE, quiet = TRUE)
wd <- list.files(tempdir(), pattern = "^adv1_est_", full.names = TRUE)
wd <- wd[length(wd)]
dat <- read.csv(file.path(wd, "data.csv"))
imp <- nm_read_nonmem(file.path(wd, "est.ctl"), data_path = file.path(wd, "data.csv"))
parts <- LibeRation:::.nm_bench_parts_for_mode(nm_ctl_template(1L, 1L), "est")
imp$model$INPUT <- parts$input_cols
fit30 <- nm_est(
  imp$model, structure(list(data = dat), class = "nm_dataset"),
  method = "FOCEI", grad = "auto", pk_engine = "cpp",
  max_outer = 30L, tol = 1e-5, control = list(maxit = 500, factr = 1e8)
)
ext <- nm_bench_read_ext(file.path(wd, "est.ext"))
cmp <- nm_bench_compare(fit30, ext, rtol = 0.15)
cat("outer n:", length(fit30$outer), " final obj:", fit30$objective, "\n")
cat("theta rel:", cmp$d_theta, " omega rel:", cmp$d_omega, " sigma rel:", cmp$d_sigma, "\n")
cat("obj rel:", abs(cmp$obj_rcpp - cmp$obj_nm) / cmp$obj_nm, " ok:", cmp$ok, "\n")
