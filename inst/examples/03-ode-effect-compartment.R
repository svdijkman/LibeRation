# Linked PK/effect-site ODE model using the stiff-capable ADVAN13 path.
library(LibeRation)

effect_model <- nm_model(
  INPUT = c("ID", "TIME", "EVID", "AMT", "CMT", "DV", "MDV"),
  ADVAN = 13, TRANS = 1, DOSECMP = 1, OBSCMP = 2,
  PRED = paste(
    "K=THETA(1)*exp(ETA(1))",
    "V=THETA(2)",
    "KE0=THETA(3)*exp(ETA(2))",
    "S1=V",
    "S2=1",
    sep = "\n"
  ),
  DES = paste(
    "DADT(1)=-K*A(1)",
    "DADT(2)=KE0*(A(1)/V-A(2))",
    sep = "\n"
  ),
  ERROR = "Y=F+ERR(1)",
  THETAS = data.frame(
    THETA = 1:3, Value = c(0.1, 20, 0.5),
    Lower = c(0.001, 1, 0.001), Upper = c(2, 100, 10)
  ),
  OMEGAS = data.frame(OMEGA = 1:2, Value = c(0.08, 0.08)),
  SIGMAS = data.frame(SIGMA = 1, Value = 0.02)
)

effect_design <- data.frame(
  ID = 1, TIME = c(0, 0.25, 0.5, 1, 2, 4, 8, 12),
  EVID = c(1L, rep(0L, 7)), AMT = c(100, rep(0, 7)),
  CMT = c(1L, rep(2L, 7)), DV = NA_real_,
  MDV = c(1L, rep(0L, 7))
)

effect_prediction <- nm_simulate(
  effect_model, effect_design, random_effects = FALSE, residual = FALSE
)
