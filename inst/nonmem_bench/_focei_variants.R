pkgload::load_all("C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation", compile = FALSE, quiet = TRUE)
source("C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation/inst/nonmem_bench/_focei_nm_formula.R", local = FALSE)

.nm_focei_subject_obj_v <- function(model, subj, theta, omega, sigma, eta, include_gog = TRUE, include_t5 = TRUE, pk_engine = "cpp") {
  err <- model$LIK_CONFIG$error %||% "propadd"
  eta <- as.numeric(eta)
  n_eta <- length(eta)
  if (n_eta == 0L || length(omega) == 0L) return(0)
  gh <- .nm_focei_subject_G(model, subj, theta, omega, eta, sigma, pk_engine = pk_engine)
  F <- as.numeric(gh$F); G <- gh$G; obs <- as.numeric(gh$obs)
  n_obs <- length(F)
  if (n_obs == 0L) return(0)
  OM <- .nm_omega_mat(omega, n_eta)
  invOM <- solve(OM)
  logdetOM <- as.numeric(determinant(OM, logarithm = TRUE)$modulus)
  s1 <- sigma[1]; s2 <- if (length(sigma) >= 2) sigma[2] else 0
  term12 <- 0; gvg <- matrix(0, n_eta, n_eta)
  for (j in seq_len(n_obs)) {
    gj <- G[j, , drop = FALSE]
    gog <- if (include_gog) as.numeric(gj %*% OM %*% t(gj)) else 0
    rj <- if (err == "prop") (s1 * F[j])^2 else if (err == "add") s1^2 else (s1 * F[j])^2 + s2^2
    vj <- max(gog + rj, 1e-15)
    term12 <- term12 + log(vj) + (obs[j] - F[j])^2 / vj
    gvg <- gvg + (1 / vj) * (t(gj) %*% gj)
  }
  term3 <- as.numeric(t(eta) %*% invOM %*% eta)
  term5 <- if (include_t5) as.numeric(determinant(invOM + gvg, logarithm = TRUE)$modulus) else 0
  term12 + term3 + logdetOM + term5
}

wd <- "C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation/inst/nonmem_bench/_est_adv1_work"
ext <- nm_bench_read_ext(file.path(wd, "est.ext"))
dat <- read.csv(file.path(wd, "data.csv"))
parts <- LibeRation:::.nm_bench_parts_for_mode(nm_ctl_template(1L, 1L), "est")
tmp <- tempfile(fileext = ".ctl"); writeLines(nm_bench_ctl(parts, mode = "est"), tmp)
imp <- nm_read_nonmem(tmp, data_path = file.path(wd, "data.csv"))
model <- imp$model
phi <- read.table(file.path(wd, "est.phi"), skip = 2, header = FALSE)
obj_nm <- as.numeric(phi[[5]])

for (i in 1:6) {
  sub <- dat[dat$ID == i, ]
  eta <- c(as.numeric(phi[[3]][i]))
  full <- .nm_focei_subject_obj(model, sub, ext$theta, ext$omega, ext$sigma, eta, pk_engine = "cpp")
  no_gog <- .nm_focei_subject_obj_v(model, sub, ext$theta, ext$omega, ext$sigma, eta, FALSE, TRUE, "cpp")
  no_t5 <- .nm_focei_subject_obj_v(model, sub, ext$theta, ext$omega, ext$sigma, eta, TRUE, FALSE, "cpp")
  cat(sprintf("ID%d nm=%.3f full=%.3f no_gog=%.3f no_t5=%.3f\n", i, obj_nm[i], full, no_gog, no_t5))
}
