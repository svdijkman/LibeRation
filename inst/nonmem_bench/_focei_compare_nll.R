pkgload::load_all("C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation", compile = FALSE, quiet = TRUE)
source("C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation/inst/nonmem_bench/_focei_nm_formula.R", local = FALSE)
wd <- "C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation/inst/nonmem_bench/_est_adv1_work"
ext <- nm_bench_read_ext(file.path(wd, "est.ext"))
dat <- read.csv(file.path(wd, "data.csv"))
parts <- LibeRation:::.nm_bench_parts_for_mode(nm_ctl_template(1L, 1L), "est")
tmp <- tempfile(fileext = ".ctl"); writeLines(nm_bench_ctl(parts, mode = "est"), tmp)
imp <- nm_read_nonmem(tmp, data_path = file.path(wd, "data.csv"))
model <- imp$model
phi <- read.table(file.path(wd, "est.phi"), skip = 2, header = FALSE)
obj_nm <- as.numeric(phi[[5]])

foce_nll <- function(i) {
  sub <- dat[dat$ID == i, ]
  eta <- c(as.numeric(phi[[3]][i]))
  gh <- .nm_focei_subject_G(model, sub, ext$theta, ext$omega, eta, ext$sigma, pk_engine = "cpp")
  nll_np <- nm_subject_nll_cpp(
    subs[[i]]$time, subs[[i]]$amt, subs[[i]]$rate, subs[[i]]$f1, subs[[i]]$cmt,
    subs[[i]]$evid, subs[[i]]$ss, subs[[i]]$ii, subs[[i]]$dv, subs[[i]]$obs_idx,
    eta, ext$theta, ext$omega, ext$sigma, meta$pred_lines, meta$advan, meta$trans,
    meta$obs_cmp, meta$dose_cmp, meta$n_transit, meta$use_ode, meta$model_ss, TRUE
  )
  obj_i <- .nm_focei_subject_obj(model, sub, ext$theta, ext$omega, eta, ext$sigma, pk_engine = "cpp")
  cat(sprintf("ID%d focei=%.4f nll=%.4f nm=%.4f\n", i, obj_i, nll_np, obj_nm[i]))
}

ds <- structure(list(data = dat), class = "nm_dataset")
subs <- .nm_cpp_subjects_cached(model, ds)
meta <- .nm_cpp_meta(model)
for (i in 1:6) foce_nll(i)
