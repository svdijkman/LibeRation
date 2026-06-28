pkgload::load_all("C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation", compile = FALSE, quiet = TRUE)
wd <- "C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation/inst/nonmem_bench/_est_adv1_work"
ext <- nm_bench_read_ext(file.path(wd, "est.ext"))
th <- ext$theta; om <- ext$omega; sg <- sqrt(ext$sigma)
dat <- read.csv(file.path(wd, "data.csv"))
parts <- LibeRation:::.nm_bench_parts_for_mode(nm_ctl_template(1L, 1L), "est")
parts$data_file <- "data.csv"
tmp <- tempfile(fileext = ".ctl"); writeLines(nm_bench_ctl(parts, mode = "est"), tmp)
imp <- nm_read_nonmem(tmp, data_path = file.path(wd, "data.csv"))
model <- imp$model
phi <- utils::read.table(file.path(wd, "est.phi"), skip = 2, header = FALSE)
eta_nm <- as.numeric(phi[[5]])  # wait check column
# Actually OBJ is col 5, ETA col 3
eta_nm <- as.numeric(phi[[3]]); obj_nm_sub <- as.numeric(phi[[5]])
subs <- .nm_cpp_subjects_cached(model, structure(list(data = dat), class = "nm_dataset"))
for (i in seq_along(subs)) {
  eta_i <- c(eta_nm[i])
  nll_i <- nm_subject_nll_cpp(
    subs[[i]]$time, subs[[i]]$amt, subs[[i]]$rate, subs[[i]]$f1,
    subs[[i]]$cmt, subs[[i]]$evid, subs[[i]]$ss, subs[[i]]$ii,
    subs[[i]]$dv, subs[[i]]$obs_idx, eta_i, th, om, sg,
    .nm_cpp_meta(model)$pred_lines, .nm_cpp_meta(model)$advan,
    .nm_cpp_meta(model)$trans, .nm_cpp_meta(model)$obs_cmp,
    .nm_cpp_meta(model)$dose_cmp, .nm_cpp_meta(model)$n_transit,
    .nm_cpp_meta(model)$use_ode, .nm_cpp_meta(model)$model_ss,
    TRUE
  )
  cat("sub", i, " rcpp_nll=", nll_i, " nm_obj=", obj_nm_sub[i], " diff=", nll_i - obj_nm_sub[i], "\n")
}
