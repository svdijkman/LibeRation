pkgload::load_all("C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation", compile = FALSE, quiet = TRUE)
wd <- "C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation/inst/nonmem_bench/_est_adv1_work"
ext <- nm_bench_read_ext(file.path(wd, "est.ext"))
th <- ext$theta; om <- ext$omega; sg <- ext$sigma
dat <- read.csv(file.path(wd, "data.csv"))
parts <- LibeRation:::.nm_bench_parts_for_mode(nm_ctl_template(1L, 1L), "est")
tmp <- tempfile(fileext = ".ctl"); writeLines(nm_bench_ctl(parts, mode = "est"), tmp)
imp <- nm_read_nonmem(tmp, data_path = file.path(wd, "data.csv"))
model <- imp$model
cat("error:", model$LIK_CONFIG$error, " sigma:", sg, " omega:", om, "\n")
phi <- utils::read.table(file.path(wd, "est.phi"), skip = 2, header = FALSE)
eta_nm <- as.numeric(phi[[3]])
for (i in 1:2) {
  sub <- dat[dat$ID == i, ]
  pr <- .nm_subject_ipred(model, sub, th, om, c(eta_nm[i]), sg, pk_engine = "cpp")
  obs <- sub$DV[sub$MDV == 0]
  cat("\nID", i, " eta=", eta_nm[i], "\n")
  print(data.frame(obs = obs, F = pr$F, rel = (obs - pr$F) / pmax(abs(pr$F), 1e-8)))
  # manual prop residual nll
  s1 <- sg[1]
  nll_m <- sum(log((s1 * pr$F)^2) + (obs - pr$F)^2 / (s1 * pr$F)^2)
  prior <- eta_nm[i]^2 / om[1] + log(om[1])
  cat("manual resid=", nll_m, " prior=", prior, " total=", nll_m + prior, "\n")
}
