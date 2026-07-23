# Intermittent infusion at analytical periodic steady state.
library(LibeRation)

infusion_model <- nm_model(
  INPUT = c("ID", "TIME", "EVID", "AMT", "RATE", "II", "SS", "DV", "MDV"),
  ADVAN = 1, TRANS = 2, DOSECMP = 1, OBSCMP = 1,
  PRED = "CL=THETA(1)*exp(ETA(1));V=THETA(2);S1=V",
  ERROR = "Y=F*(1+ERR(1))",
  THETAS = data.frame(THETA = 1:2, Value = c(2, 20)),
  OMEGAS = data.frame(OMEGA = 1, Value = 0.1),
  SIGMAS = data.frame(SIGMA = 1, Value = 0.02),
  LIK_CONFIG = nm_lik_config(
    error = "proportional", sigma_parameterization = "variance"
  )
)

infusion_design <- data.frame(
  ID = 1, TIME = c(0, 1, 4, 8, 12),
  EVID = c(1L, rep(0L, 4)), AMT = c(100, rep(0, 4)),
  RATE = c(25, rep(0, 4)), II = c(12, rep(0, 4)),
  SS = c(1L, rep(0L, 4)), DV = NA_real_,
  MDV = c(1L, rep(0L, 4))
)

infusion_prediction <- nm_simulate(
  infusion_model, infusion_design, random_effects = FALSE, residual = FALSE
)
