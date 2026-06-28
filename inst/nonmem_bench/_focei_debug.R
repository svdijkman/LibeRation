pkgload::load_all("C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation", compile = FALSE, quiet = TRUE)
wd <- "C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation/inst/nonmem_bench/_est_adv1_work"
parts <- LibeRation:::.nm_bench_parts_for_mode(nm_ctl_template(1L, 1L), "est")
tmp <- tempfile(fileext = ".ctl"); writeLines(nm_bench_ctl(parts, mode = "est"), tmp)
imp <- nm_read_nonmem(tmp, data_path = file.path(wd, "data.csv"))
model <- imp$model
cat("n_eta:", LibeRation:::.nm_n_eta(model), "\n")
cat("omegas:\n"); print(model$OMEGAS)
phi <- utils::read.table(file.path(wd, "est.phi"), skip = 2, header = FALSE)
cat("eta1:", phi[[3]][1], " len:", length(c(as.numeric(phi[[3]][1]))), "\n")
ext <- nm_bench_read_ext(file.path(wd, "est.ext"))
cat("ext omega:", ext$omega, "\n")
