pkgload::load_all("C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation", compile = FALSE, quiet = TRUE)
wd <- "C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation/inst/nonmem_bench/_est_adv1_work"
res <- nm_bench_pilot(advan = 1L, trans = 1L, n_per = 2L, seed = 2024L,
                      method = "FOCEI", run_nonmem = TRUE, work_dir = wd)
ext <- res$nm_ext
th <- ext$theta; om <- ext$omega; sg <- sqrt(ext$sigma)
imp <- nm_read_nonmem(res$est_ctl_path, data_path = res$data)
model <- imp$model
dat <- .nm_prepare_data(structure(list(data = res$data), class = "nm_dataset"), model$INPUT, model)
eta <- .nm_fit_all_eta_modes(model, dat, th, om, sg, NULL, "cpp", list(maxit = 400))
cat("Rcpp eta modes at NM pop params:\n")
print(eta)
# read NM posthoc eta from .phi if present
phi_path <- file.path(wd, "est.phi")
if (file.exists(phi_path)) {
  phi <- read.table(phi_path, skip = 1, header = FALSE)
  cat("NM phi file rows:", nrow(phi), " cols:", ncol(phi), "\n")
  print(head(phi))
}
