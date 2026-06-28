pkgload::load_all("C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation", compile = FALSE, quiet = TRUE)

.nm_focei_subject_G <- function(model, subj, theta, omega, eta, sigma, pk_engine = "cpp") {
  n_eta <- length(eta)
  pred0 <- .nm_subject_ipred(model, subj, theta, omega, eta, sigma, pk_engine = pk_engine)
  f0 <- as.numeric(pred0$F)
  obs <- as.numeric(subj$DV[pred0$obs_idx])
  if (length(f0) == 0L) {
    return(list(F = f0, obs = obs, G = matrix(0, 0, n_eta)))
  }
  n_obs <- length(f0)
  G <- matrix(0, nrow = n_obs, ncol = n_eta)
  eps <- 1e-4
  for (k in seq_len(n_eta)) {
    etap <- eta; etam <- eta
    etap[k] <- etap[k] + eps
    etam[k] <- etam[k] - eps
    fp <- as.numeric(.nm_subject_ipred(
      model, subj, theta, omega, etap, sigma, pk_engine = pk_engine
    )$F)
    fm <- as.numeric(.nm_subject_ipred(
      model, subj, theta, omega, etam, sigma, pk_engine = pk_engine
    )$F)
    G[, k] <- (fp - fm) / (2 * eps)
  }
  list(F = f0, obs = obs, G = G)
}

.nm_omega_mat <- function(omega, n_eta) {
  if (length(omega) >= 3L && n_eta >= 2L) {
    om <- matrix(0, n_eta, n_eta)
    om[1, 1] <- omega[1]
    om[2, 2] <- omega[2]
    om[1, 2] <- om[2, 1] <- omega[3]
    if (n_eta > 2L) diag(om)[3:n_eta] <- omega[4:n_eta]
    return(om)
  }
  v <- pmax(omega[seq_len(n_eta)], 1e-15)
  if (n_eta == 1L) {
    return(matrix(v[1L], 1L, 1L))
  }
  diag(v)
}

.nm_focei_subject_obj <- function(model, subj, theta, omega, sigma, eta, pk_engine = "cpp") {
  err <- model$LIK_CONFIG$error %||% "propadd"
  eta <- as.numeric(eta)
  n_eta <- length(eta)
  if (n_eta == 0L || length(omega) == 0L) {
    return(0)
  }
  gh <- .nm_focei_subject_G(model, subj, theta, omega, eta, sigma, pk_engine = pk_engine)
  F <- as.numeric(gh$F)
  G <- gh$G
  obs <- as.numeric(gh$obs)
  n_obs <- length(F)
  if (n_obs == 0L) return(0)
  OM <- .nm_omega_mat(omega, n_eta)
  invOM <- solve(OM)
  logdetOM <- as.numeric(determinant(OM, logarithm = TRUE)$modulus)
  s1 <- sigma[1]; s2 <- if (length(sigma) >= 2) sigma[2] else 0
  term12 <- 0
  gvg <- matrix(0, n_eta, n_eta)
  for (j in seq_len(n_obs)) {
    gj <- G[j, , drop = FALSE]
    gog <- gj %*% OM %*% t(gj)
    if (err == "prop") {
      rj <- (s1 * F[j])^2
    } else if (err == "add") {
      rj <- s1^2
    } else {
      rj <- (s1 * F[j])^2 + s2^2
    }
    vj <- max(as.numeric(gog) + rj, 1e-15)
    term12 <- term12 + log(vj) + (obs[j] - F[j])^2 / vj
    gvg <- gvg + (1 / vj) * (t(gj) %*% gj)
  }
  term3 <- as.numeric(t(eta) %*% invOM %*% eta)
  term5 <- as.numeric(determinant(invOM + gvg, logarithm = TRUE)$modulus)
  term12 + term3 + logdetOM + term5
}

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
total <- 0
for (i in seq_along(unique(dat$ID))) {
  sub <- dat[dat$ID == i, ]
  obj_i <- .nm_focei_subject_obj(model, sub, th, om, sg, c(eta_nm[i]), pk_engine = "cpp")
  total <- total + obj_i
  cat(sprintf("ID%d rcpp=%.4f nm=%.4f diff=%.4f\n", i, obj_i, obj_nm[i], obj_nm[i] - obj_i))
}
cat(sprintf("TOTAL rcpp=%.4f nm=%.4f ext=%.4f\n", total, sum(obj_nm), ext$obj))
