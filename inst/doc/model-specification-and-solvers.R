## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(LibeRation)


## ----output-columns, eval=FALSE-----------------------------------------------
# model <- nm_model(
#   INPUT = c("ID", "TIME", "EVID", "AMT", "DV", "MDV"),
#   OUTPUT = c("PRED", "IPRED", "CWRES", "CL", "V", "A1"),
#   ADVAN = 1,
#   PRED = "CL=THETA(1)*exp(ETA(1)); V=THETA(2); K=CL/V; S1=V",
#   ERROR = "Y=F+ERR(1)",
#   THETAS = data.frame(THETA = 1:2, Value = c(2, 20)),
#   OMEGAS = data.frame(OMEGA = 1, Value = 0.1),
#   SIGMAS = data.frame(SIGMA = 1, Value = 0.1)
# )
# 
# nm_model_outputs(model)


## ----ode-example, eval=FALSE--------------------------------------------------
# ode_model <- nm_model(
#   INPUT = c("ID", "TIME", "EVID", "AMT", "CMT", "DV", "MDV"),
#   ADVAN = 13,
#   DOSECMP = 1,
#   OBSCMP = 2,
#   PRED = "KFAST=THETA(1); KSLOW=THETA(2); S2=THETA(3)",
#   DES = paste(
#     "DADT(1) = -KFAST * A(1)",
#     "DADT(2) = KFAST * A(1) - KSLOW * A(2)",
#     sep = "\n"
#   ),
#   ERROR = "Y = F + ERR(1)",
#   THETAS = data.frame(THETA = 1:3, Value = c(100, 1, 20)),
#   SIGMAS = data.frame(SIGMA = 1, Value = 0.1),
#   ODE_CONTROL = list(rtol = 1e-7, atol = 1e-9, max_steps = 100000)
# )


## ----parameter-tables---------------------------------------------------------
theta <- data.frame(
  THETA = 1:2,
  Lower = c(0.01, 1),
  Value = c(2, 20),
  Upper = c(20, 200),
  FIX = c(FALSE, FALSE)
)

omega_full <- data.frame(
  OMEGA = 1:3,
  ROW = c(1, 2, 2),
  COL = c(1, 1, 2),
  Value = c(0.10, -0.02, 0.20)
)


## ----likelihood-configuration-------------------------------------------------
priors <- rbind(
  nm_prior("THETA1", "normal", mean = 2, sd = 0.5),
  nm_prior("SIGMA1", "lognormal", mean = -2, sd = 0.5),
  nm_prior("OMEGA1", "inverse_gamma", shape = 3, rate = 0.2)
)

likelihood <- nm_lik_config(
  omega = "full",
  sigma_parameterization = "variance",
  blq_method = "m3",
  lloq = 0.05,
  iov = 1,
  occasion_col = "OCC",
  priors = priors,
  mixtures = nm_mixture(c(0.8, 0.2), c("typical", "slow"))
)


## ----cpp-model, eval=FALSE----------------------------------------------------
# cpp_model <- nm_model_cpp(
#   CPP = "CL = THETA(1) * exp(ETA(1)); V = THETA(2); S1 = V;",
#   INPUT = c("ID", "TIME", "EVID", "AMT", "CMT", "DV", "MDV"),
#   ADVAN = 1,
#   ERROR = "Y = F * (1 + ERR(1))",
#   THETAS = theta,
#   OMEGAS = data.frame(OMEGA = 1, Value = 0.1),
#   SIGMAS = data.frame(SIGMA = 1, Value = 0.1)
# )

