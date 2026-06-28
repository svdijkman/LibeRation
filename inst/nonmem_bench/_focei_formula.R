pkgload::load_all("C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation", compile = FALSE, quiet = TRUE)
wd <- "C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation/inst/nonmem_bench/_est_adv1_work"
ext <- nm_bench_read_ext(file.path(wd, "est.ext"))
th <- ext$theta; om <- ext$omega; sg <- ext$sigma
dat <- read.csv(file.path(wd, "data.csv"))
parts <- LibeRation:::.nm_bench_parts_for_mode(nm_ctl_template(1L, 1L), "est")
tmp <- tempfile(fileext = ".ctl"); writeLines(nm_bench_ctl(parts, mode = "est"), tmp)
imp <- nm_read_nonmem(tmp, data_path = file.path(wd, "data.csv"))
model <- imp$model
ds <- structure(list(data = dat), class = "nm_dataset")
subs <- .nm_cpp_subjects_cached(model, ds)
meta <- .nm_cpp_meta(model)
phi <- utils::read.table(file.path(wd, "est.phi"), skip = 2, header = FALSE)
eta_nm <- as.numeric(phi[[3]]); obj_nm <- as.numeric(phi[[5]])
log2pi <- log(2 * pi)

total_nm <- sum(obj_nm)
total_rcpp <- 0
total_f2 <- 0

for (i in seq_along(subs)) {
  sub <- subs[[i]]
  eta_i <- c(eta_nm[i])
  pr <- .nm_subject_ipred(model, dat[dat$ID == i, ], th, om, eta_i, sg, pk_engine = "cpp")
  obs <- dat$DV[dat$ID == i & dat$MDV == 0]
  n_obs <- length(obs)
  nll <- nm_subject_nll_cpp(
    sub$time, sub$amt, sub$rate, sub$f1, sub$cmt, sub$evid, sub$ss, sub$ii,
    sub$dv, sub$obs_idx, eta_i, th, om, sg,
    meta$pred_lines, meta$advan, meta$trans, meta$obs_cmp, meta$dose_cmp,
    meta$n_transit, meta$use_ode, meta$model_ss, TRUE
  )
  nll_np <- nm_subject_nll_cpp(
    sub$time, sub$amt, sub$rate, sub$f1, sub$cmt, sub$evid, sub$ss, sub$ii,
    sub$dv, sub$obs_idx, eta_i, th, om, sg,
    meta$pred_lines, meta$advan, meta$trans, meta$obs_cmp, meta$dose_cmp,
    meta$n_transit, meta$use_ode, meta$model_ss, FALSE
  )
  prior <- nll - nll_np
  eta_mat <- matrix(eta_i, nrow = 1)
  inter <- nm_focei_interaction_cpp(
    subs[i], eta_mat, th, om, sg, meta$pred_lines, meta$advan, meta$trans,
    meta$obs_cmp, meta$dose_cmp, meta$n_transit, meta$use_ode, meta$model_ss
  )
  n_eta <- length(eta_i)
  # candidate formulas
  f1 <- nll + inter
  f2 <- 2 * (nll + inter)
  f3 <- nll + inter + n_obs * log2pi + n_eta * log2pi
  f4 <- 2 * nll_np + prior + 2 * inter
  f5 <- nll + 2 * inter
  f6 <- 2 * nll_np + prior + inter + n_obs * log2pi + n_eta * log2pi
  total_rcpp <- total_rcpp + f1
  total_f2 <- total_f2 + f2
  cat(sprintf(
    "ID%d nm=%.4f f1=%.4f f2=%.4f f3=%.4f f4=%.4f f5=%.4f f6=%.4f resid=%.4f prior=%.4f inter=%.4f\n",
    i, obj_nm[i], f1, f2, f3, f4, f5, f6, nll_np, prior, inter
  ))
}
cat(sprintf("\nTOTAL nm=%.4f f1=%.4f f2=%.4f ext_obj=%.4f\n", total_nm, total_rcpp, total_f2, ext$obj))
