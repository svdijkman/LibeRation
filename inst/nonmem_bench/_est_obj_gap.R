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
phi <- utils::read.table(file.path(wd, "est.phi"), skip = 2, header = FALSE)
eta_nm <- as.numeric(phi[[3]]); obj_nm <- as.numeric(phi[[5]])
eta_mat <- matrix(eta_nm, ncol = 1)
obj <- .nm_focei_objective(model, ds, th, om, sg, eta_mat, pk_engine = "cpp")
cat("FOCEI obj=", obj, " NM=", ext$obj, " gap=", ext$obj - obj, "\n")
cat("phi sum=", sum(obj_nm), " rel=", abs(ext$obj - obj) / ext$obj, "\n")
