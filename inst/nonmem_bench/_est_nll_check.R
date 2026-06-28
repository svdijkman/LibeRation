pkgload::load_all("C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation", compile = FALSE, quiet = TRUE)
res <- nm_bench_pilot(advan = 1L, trans = 1L, n_per = 2L, seed = 2024L, method = "FOCEI", run_nonmem = TRUE)
fit <- res$rcpp_fit
ext <- res$nm_ext
imp <- nm_read_nonmem(res$est_ctl_path, data_path = res$data)
model <- imp$model
dat <- structure(list(data = res$data), class = "nm_dataset")

# NM params from .ext (sigma is variance on EPS for prop error)
th_nm <- ext$theta
om_nm <- ext$omega
sg_nm_var <- ext$sigma
sg_nm_sd <- sqrt(sg_nm_var)

cat("NM obj:", ext$obj, "\n")
for (label in c("rcpp_at_rcpp", "rcpp_at_nm_sd", "rcpp_at_nm_var")) {
  th <- th_nm; om <- om_nm
  sg <- switch(label,
    rcpp_at_rcpp = fit$sigma,
    rcpp_at_nm_sd = sg_nm_sd,
    rcpp_at_nm_var = sg_nm_var
  )
  th <- if (label == "rcpp_at_rcpp") fit$theta else th_nm
  eta <- res$rcpp_fit$eta
  nll <- nm_nll(model, dat, th, om, sg, eta = eta, pk_engine = "cpp")
  inter <- .nm_focei_interaction(model, dat, th, om, sg, eta)
  obj <- nll + inter
  cat(label, " sigma=", paste(signif(sg, 4), collapse = ","),
      " obj=", obj, " nll=", nll, " inter=", inter, "\n")
}

# error type
cat("lik config:", paste(model$LIK_CONFIG$error, model$LIK_CONFIG$omega), "\n")
