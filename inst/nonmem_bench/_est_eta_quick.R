pkgload::load_all("C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation", compile = FALSE, quiet = TRUE)
wd <- "C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation/inst/nonmem_bench/_est_adv1_work"
ext <- nm_bench_read_ext(file.path(wd, "est.ext"))
th <- ext$theta; om <- ext$omega; sg <- sqrt(ext$sigma)
dat <- read.csv(file.path(wd, "data.csv"))
parts <- LibeRation:::.nm_bench_parts_for_mode(nm_ctl_template(1L, 1L), "est")
parts$data_file <- "data.csv"
ctl <- nm_bench_ctl(parts, mode = "est")
tmp <- tempfile(fileext = ".ctl"); writeLines(ctl, tmp)
imp <- nm_read_nonmem(tmp, data_path = file.path(wd, "data.csv"))
model <- imp$model
dat2 <- .nm_prepare_data(structure(list(data = dat), class = "nm_dataset"), model$INPUT, model)
eta_rcpp <- .nm_fit_all_eta_modes(model, dat2, th, om, sg, NULL, "cpp", list(maxit = 400))
phi <- utils::read.table(file.path(wd, "est.phi"), skip = 2, header = FALSE)
eta_nm <- as.numeric(phi[[3]])
cat("eta NM:", eta_nm, "\n")
cat("eta rcpp:\n"); print(eta_rcpp[, 1])
cat("max abs diff:", max(abs(eta_rcpp[, 1] - eta_nm)), "\n")

ds <- structure(list(data = dat), class = "nm_dataset")
eta_mat <- matrix(eta_nm, ncol = 1)
nll_nm_eta <- nm_nll(model, ds, th, om, sg, eta = eta_mat, pk_engine = "cpp")
inter <- .nm_focei_interaction(model, ds, th, om, sg, eta_mat)
cat("obj with NM eta:", nll_nm_eta + inter, " NM total:", ext$obj, "\n")
