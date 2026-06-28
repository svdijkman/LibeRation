pkgload::load_all("C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation", compile = FALSE, quiet = TRUE)
wd <- "C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation/inst/nonmem_bench/_est_adv1_work"
ext <- nm_bench_read_ext(file.path(wd, "est.ext"))
th <- ext$theta; om <- ext$omega; sg <- sqrt(ext$sigma)
dat <- read.csv(file.path(wd, "data.csv"))
parts <- LibeRation:::.nm_bench_parts_for_mode(nm_ctl_template(1L, 1L), "est")
tmp <- tempfile(fileext = ".ctl"); writeLines(nm_bench_ctl(parts, mode = "est"), tmp)
imp <- nm_read_nonmem(tmp, data_path = file.path(wd, "data.csv"))
model <- imp$model
phi <- utils::read.table(file.path(wd, "est.phi"), skip = 2, header = FALSE)
eta_nm <- as.numeric(phi[[3]])
obj_nm_sub <- as.numeric(phi[[5]])
ds <- structure(list(data = dat), class = "nm_dataset")
subs <- .nm_cpp_subjects_cached(model, ds)
meta <- .nm_cpp_meta(model)
for (i in 1:2) {
  sub <- dat[dat$ID == i, ]
  pr <- .nm_subject_ipred(model, sub, th, om, c(eta_nm[i]), sg, pk_engine = "cpp")
  obs <- sub$DV[sub$MDV == 0]
  cat("\nID", i, " n_obs=", length(obs), "\n")
  cat("F rcpp:", pr$F, "\n")
  cat("DV:", obs, "\n")
  eta_mat <- matrix(eta_nm[i], nrow = 1)
  inter_i <- nm_focei_interaction_cpp(
    subs[i], eta_mat, th, om, sg, meta$pred_lines, meta$advan, meta$trans,
    meta$obs_cmp, meta$dose_cmp, meta$n_transit, meta$use_ode, meta$model_ss
  )
  cat("interaction_i=", inter_i, " nm_obj=", obj_nm_sub[i], "\n")
}
