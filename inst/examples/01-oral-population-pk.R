# Oral population PK: simulation followed by a FOCEI fit.
library(LibeRation)

oral_model <- nm_model(
  INPUT = c("ID", "TIME", "EVID", "AMT", "CMT", "DV", "MDV"),
  ADVAN = 2, TRANS = 2, DOSECMP = 1, OBSCMP = 2,
  PRED = paste(
    "KA=THETA(1)*exp(ETA(1))",
    "CL=THETA(2)*exp(ETA(2))",
    "V=THETA(3)",
    "S2=V",
    sep = "\n"
  ),
  ERROR = "Y=F*(1+ERR(1))",
  THETAS = data.frame(
    THETA = 1:3, Value = c(1.1, 2, 20),
    Lower = c(0.05, 0.1, 5), Upper = c(5, 10, 80)
  ),
  OMEGAS = data.frame(OMEGA = 1:2, Value = c(0.12, 0.1)),
  SIGMAS = data.frame(SIGMA = 1, Value = 0.02),
  LIK_CONFIG = nm_lik_config(
    error = "proportional", sigma_parameterization = "variance"
  )
)

oral_design <- do.call(rbind, lapply(1:20, function(id) {
  data.frame(
    ID = id, TIME = c(0, 0.5, 1, 2, 4, 8, 12, 24),
    EVID = c(1L, rep(0L, 7)), AMT = c(100, rep(0, 7)),
    CMT = c(1L, rep(2L, 7)), DV = NA_real_,
    MDV = c(1L, rep(0L, 7))
  )
}))

oral_data <- nm_simulate(
  oral_model, oral_design, random_effects = TRUE, residual = TRUE,
  seed = 20260723
)

fit_oral_example <- function(covariance = TRUE) {
  nm_est(
    oral_model, oral_data, method = "FOCEI", covariance = covariance,
    maxit = 150, eta_maxit = 100
  )
}
