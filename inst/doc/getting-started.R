## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----model--------------------------------------------------------------------
library(LibeRation)

model <- nm_model(
  INPUT = c("ID", "TIME", "EVID", "AMT", "CMT", "DV", "MDV"),
  ADVAN = 2,
  TRANS = 2,
  PRED = paste(
    "CL = THETA(1) * exp(ETA(1))",
    "V = THETA(2)",
    "KA = THETA(3)",
    "S2 = V",
    sep = "\n"
  ),
  ERROR = "Y = F * (1 + ERR(1))",
  THETAS = data.frame(
    THETA = 1:3,
    Value = c(2.8, 32, 1.4),
    Lower = c(0.01, 0.1, 0.01),
    Upper = c(20, 200, 10)
  ),
  OMEGAS = data.frame(OMEGA = 1, Value = 0.09),
  SIGMAS = data.frame(SIGMA = 1, Value = 0.15)
)
model

## ----data---------------------------------------------------------------------
event_data <- data.frame(
  ID = 1L,
  TIME = c(0, 0.5, 1, 2, 4, 8, 12),
  EVID = c(1L, rep(0L, 6)),
  AMT = c(320, rep(0, 6)),
  CMT = c(1L, rep(2L, 6)),
  DV = NA_real_,
  MDV = c(1L, rep(0L, 6))
)

validated <- nm_dataset(event_data)
validated

## ----simulate-----------------------------------------------------------------
simulated <- nm_simulate(
  model,
  validated,
  nsim = 1,
  random_effects = TRUE,
  residual = TRUE,
  seed = 20260715
)
simulated[c("ID", "TIME", "EVID", "DV", "IPRED")]

## ----estimate, eval=FALSE-----------------------------------------------------
# fit <- nm_est(
#   model,
#   simulated,
#   method = "FOCEI",
#   maxit = 200,
#   eta_maxit = 100,
#   covariance = TRUE,
#   covariance_type = "auto",
#   print_every = 10,
#   n_cores = 1
# )
# 
# summary(fit)
# fit$covariance

## ----next, eval=FALSE---------------------------------------------------------
# gof <- nm_gof(fit)
# vpc <- nm_vpc(fit, nsim = 200, stratify = "SEX")
# npde <- nm_npde(fit, nsim = 200)
# npc <- nm_npc(fit, nsim = 200)
# 
# liber_gui(model = model, data = simulated)

