pkgload::load_all("C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation", compile = FALSE, quiet = TRUE)
wd <- "C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation/inst/nonmem_bench/_est_adv1_work"
ext <- nm_bench_read_ext(file.path(wd, "est.ext"))
om <- ext$omega
eta <- 0.577304
n_eta <- length(eta)
cat("n_eta", n_eta, "om", om, "class", class(om), "len", length(om), "\n")
OM <- diag(pmax(om[seq_len(n_eta)], 1e-15))
cat("dim OM", dim(OM), "\n")

.nm_omega_mat <- function(omega, n_eta) {
  if (length(omega) >= 3L && n_eta >= 2L) {
    om <- matrix(0, n_eta, n_eta)
    om[1, 1] <- omega[1]
    om[2, 2] <- omega[2]
    om[1, 2] <- om[2, 1] <- omega[3]
    if (n_eta > 2L) diag(om)[3:n_eta] <- omega[4:n_eta]
    return(om)
  }
  diag(pmax(omega[seq_len(n_eta)], 1e-15))
}
OM2 <- .nm_omega_mat(om, n_eta)
cat("dim OM2", dim(OM2), "\n")

parts <- LibeRation:::.nm_bench_parts_for_mode(nm_ctl_template(1L, 1L), "est")
tmp <- tempfile(fileext = ".ctl"); writeLines(nm_bench_ctl(parts, mode = "est"), tmp)
imp <- nm_read_nonmem(tmp, data_path = file.path(wd, "data.csv"))
model <- imp$model
dat <- read.csv(file.path(wd, "data.csv"))
sub <- dat[dat$ID == 1, ]
pred0 <- .nm_subject_ipred(model, sub, ext$theta, om, eta, ext$sigma, pk_engine = "cpp")
cat("F len", length(pred0$F), "\n")
