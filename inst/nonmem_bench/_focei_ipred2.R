pkgload::load_all("C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation", compile = FALSE, quiet = TRUE)
wd <- "C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation/inst/nonmem_bench/_est_adv1_work"
ext <- nm_bench_read_ext(file.path(wd, "est.ext"))
dat <- read.csv(file.path(wd, "data.csv"))
parts <- LibeRation:::.nm_bench_parts_for_mode(nm_ctl_template(1L, 1L), "est")
tmp <- tempfile(fileext = ".ctl"); writeLines(nm_bench_ctl(parts, mode = "est"), tmp)
imp <- nm_read_nonmem(tmp, data_path = file.path(wd, "data.csv"))
model <- imp$model
phi <- read.table(file.path(wd, "est.phi"), skip = 2, header = FALSE)
for (i in 1:2) {
  sub <- dat[dat$ID == i, ]
  eta <- c(as.numeric(phi[[3]][i]))
  pr <- .nm_subject_ipred(model, sub, ext$theta, ext$omega, eta, ext$sigma, pk_engine = "cpp")
  cat("ID", i, "eta", eta, "F", pr$F, "DV", sub$DV[pr$obs_idx], "\n")
}
