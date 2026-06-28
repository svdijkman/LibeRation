pkgload::load_all("C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation", compile = FALSE, quiet = TRUE)
res <- nm_bench_pilot(advan = 1L, trans = 1L, n_per = 2L, seed = 2024L,
                      method = "FOCEI", run_nonmem = TRUE,
                      work_dir = "C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation/inst/nonmem_bench/_est_adv1_work")
imp <- nm_read_nonmem(res$est_ctl_path, data_path = res$data)
model <- imp$model
dat <- .nm_prepare_data(structure(list(data = res$data), class = "nm_dataset"), model$INPUT, model)
ext <- res$nm_ext

th <- ext$theta
om <- ext$omega
sg <- sqrt(ext$sigma)  # variance -> SD for rcpp prop error

eta <- .nm_fit_all_eta_modes(model, dat, th, om, sg, NULL, "cpp", list(maxit = 200))
nll <- nm_nll(model, structure(list(data = res$data), class = "nm_dataset"), th, om, sg, eta = eta, pk_engine = "cpp")
inter <- .nm_focei_interaction(model, structure(list(data = res$data), class = "nm_dataset"), th, om, sg, eta)
cat("At NM theta/omega, sigma=sqrt(var): obj=", nll + inter, " nll=", nll, " inter=", inter, " NM=", ext$obj, "\n")

# grid search sigma SD around sqrt(var)
for (s in seq(0.25, 0.40, length.out = 8)) {
  eta2 <- .nm_fit_all_eta_modes(model, dat, th, om, s, NULL, "cpp", list(maxit = 200))
  nll2 <- nm_nll(model, structure(list(data = res$data), class = "nm_dataset"), th, om, s, eta = eta2, pk_engine = "cpp")
  inter2 <- .nm_focei_interaction(model, structure(list(data = res$data), class = "nm_dataset"), th, om, s, eta2)
  cat("sigma_sd=", round(s, 4), " obj=", nll2 + inter2, "\n")
}
