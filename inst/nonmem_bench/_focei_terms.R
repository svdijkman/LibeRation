source("C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation/inst/nonmem_bench/_focei_nm_formula.R", local = TRUE)
wd <- "C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation/inst/nonmem_bench/_est_adv1_work"
ext <- nm_bench_read_ext(file.path(wd, "est.ext"))
th <- ext$theta; om <- ext$omega; sg <- ext$sigma
dat <- read.csv(file.path(wd, "data.csv"))
parts <- LibeRation:::.nm_bench_parts_for_mode(nm_ctl_template(1L, 1L), "est")
tmp <- tempfile(fileext = ".ctl"); writeLines(nm_bench_ctl(parts, mode = "est"), tmp)
imp <- nm_read_nonmem(tmp, data_path = file.path(wd, "data.csv"))
model <- imp$model
phi <- utils::read.table(file.path(wd, "est.phi"), skip = 2, header = FALSE)
eta_nm <- as.numeric(phi[[3]]); obj_nm <- as.numeric(phi[[5]])

decomp <- function(i) {
  sub <- dat[dat$ID == i, ]
  eta <- c(eta_nm[i])
  n_eta <- length(eta)
  gh <- .nm_focei_subject_G(model, sub, th, om, sg, eta, pk_engine = "cpp")
  F <- gh$F; G <- gh$G; obs <- gh$obs
  OM <- .nm_omega_mat(om, n_eta)
  invOM <- solve(OM)
  logdetOM <- as.numeric(determinant(OM, logarithm = TRUE)$modulus)
  s1 <- sg[1]
  term12 <- 0; gvg <- matrix(0, n_eta, n_eta)
  for (j in seq_along(F)) {
    gj <- G[j, , drop = FALSE]
    gog <- as.numeric(gj %*% OM %*% t(gj))
    rj <- (s1 * F[j])^2
    vj <- max(gog + rj, 1e-15)
    term12 <- term12 + log(vj) + (obs[j] - F[j])^2 / vj
    gvg <- gvg + (1 / vj) * (t(gj) %*% gj)
  }
  term3 <- as.numeric(t(eta) %*% invOM %*% eta)
  term5 <- as.numeric(determinant(invOM + gvg, logarithm = TRUE)$modulus)
  cat(sprintf("ID%d eta=%.3f nm=%.4f\n  t12=%.4f t3=%.4f logdetOM=%.4f t5=%.4f sum=%.4f\n  G range: %.3e %.3e\n",
              i, eta, obj_nm[i], term12, term3, logdetOM, term5, term12+term3+logdetOM+term5,
              min(G), max(G)))
}

for (i in c(1, 2, 3, 6)) decomp(i)
