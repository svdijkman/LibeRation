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

## ----first-class-outcomes, eval=FALSE-----------------------------------------
# count_model <- nm_model(
#   INPUT = c("ID", "TIME", "DV", "MDV"),
#   ADVAN = 1,
#   PRED = "MU=exp(THETA(1)+ETA(1)); CL=1; V=1; S1=V; F=MU",
#   THETAS = data.frame(THETA=1:2, Value=log(c(2, 3)), LOWER=-10, UPPER=10),
#   OMEGAS = data.frame(OMEGA=1, Value=0.2),
#   OUTCOMES = nm_outcome(
#     "negative_binomial", prediction="MU", dispersion="THETA(2)",
#     max_count=200
#   )
# )
# 
# simulated <- nm_simulate(count_model, count_design, residual=TRUE, seed=1)
# fit <- nm_est(count_model, count_data, method="LAPLACE")
# nm_outcome_diagnostics(fit)
# nm_vpc_count(fit)

## ----joint-outcomes, eval=FALSE-----------------------------------------------
# joint <- nm_outcomes(
#   nm_outcome("normal", name="concentration", dvid=1,
#              prediction="CP", scale="SIGMA(1)"),
#   nm_outcome("bernoulli", name="response", dvid=2,
#              prediction="PRESPONSE")
# )

## ----advanced-templates, eval=FALSE-------------------------------------------
# tmdd <- nm_model_template("tmdd")
# transit <- nm_model_template("transit_absorption", n_transit=5)
# 
# hazard <- nm_piecewise("TIME", c(24, 72), c("H1", "H2", "H3"))
# trend <- nm_spline("TIME", c(0, 12, 48, 96),
#                    c("B1", "B2", "B3"), intercept="B0")

## ----markov-likelihood, eval=FALSE--------------------------------------------
# markov <- nm_model(
#   INPUT = c("ID", "TIME", "DV", "MDV"),
#   ADVAN = 1,
#   PRED = "
#     PBASE = 1 / (1 + exp(-THETA(1)))
#     P01 = 1 / (1 + exp(-(THETA(2) + ETA(1))))
#     P11 = 1 / (1 + exp(-(THETA(3) + ETA(1))))
#     CL = 1; V = 1; S1 = V; F = PBASE
#   ",
#   ERROR = "
#     PCURRENT = ifelse(
#       FIRST == 1,
#       ifelse(DV == 1, PBASE, 1 - PBASE),
#       ifelse(
#         PREV_DV == 0,
#         ifelse(DV == 1, P01, 1 - P01),
#         ifelse(DV == 1, P11, 1 - P11)
#       )
#     )
#     LOGLIK = log(pmax(PCURRENT, 1e-12))
#   ",
#   THETAS = data.frame(
#     THETA = 1:3, Value = qlogis(c(0.3, 0.2, 0.8)),
#     LOWER = -10, UPPER = 10
#   ),
#   OMEGAS = data.frame(OMEGA = 1, Value = 0.2),
#   OUTPUT = c("PBASE", "P01", "P11")
# )
# 
# fit <- nm_est(markov, transitions, method = "LAPLACE")

## ----hidden-markov-model, eval=FALSE------------------------------------------
# hmm <- nm_model(
#   INPUT = c("ID", "TIME", "DV", "MDV"),
#   ADVAN = 1,
#   PRED = "CL = 1; V = 1; S1 = V; F = 0",
#   ERROR = "
#     PI_C = 1 / (1 + exp(-THETA(1)))
#     PI_A = 1 - PI_C
#     P_CC = 1 / (1 + exp(-THETA(2)))
#     P_CA = 1 - P_CC
#     P_AC = 1 / (1 + exp(-THETA(3)))
#     P_AA = 1 - P_AC
#     EMIT_C = ifelse(DV == 0, 1 - THETA(4), THETA(4))
#     EMIT_A = ifelse(DV == 0, 1 - THETA(5), THETA(5))
#   ",
#   THETAS = data.frame(
#     THETA = 1:5,
#     Value = c(qlogis(0.8), qlogis(0.9), qlogis(0.3), 0.1, 0.8),
#     LOWER = c(rep(-10, 3), 1e-6, 1e-6),
#     UPPER = c(rep(10, 3), 1 - 1e-6, 1 - 1e-6)
#   ),
#   HMM_CONFIG = nm_hmm_config(
#     states = c("controlled", "active"),
#     initial = c("PI_C", "PI_A"),
#     transition = matrix(
#       c("P_CC", "P_CA", "P_AC", "P_AA"), 2, 2, byrow = TRUE
#     ),
#     emission = c("EMIT_C", "EMIT_A")
#   )
# )
# 
# fit <- nm_est(hmm, observations, method = "LAPLACE")
# states <- nm_hmm_decode(fit, method = "all")
# gof <- nm_gof(fit) # includes HMM_STATE and HMM_PROB_* columns

## ----hsmm-model, eval=FALSE---------------------------------------------------
# semi_markov <- nm_hsmm_config(
#   states = c("controlled", "active"),
#   initial = c("PI_C", "PI_A"),
#   transition = matrix(c("P_CC", "P_CA", "P_AC", "P_AA"), 2, 2,
#                       byrow = TRUE),
#   dwell = matrix(c("D_C1", "D_C2", "D_C3",
#                    "D_A1", "D_A2", "D_A3"), 2, 3, byrow = TRUE),
#   emission = c("EMIT_C", "EMIT_A")
# )

## ----random-effect-design, eval=FALSE-----------------------------------------
# design <- nm_re_config(
#   nm_re_block("site", "SITE", 1),
#   nm_re_block("patient", "ID", 2:3),
#   nm_re_block("reader", "READER", 4)
# )
# 
# model <- nm_model(
#   # ordinary INPUT/PRED/ERROR/parameter declarations ...
#   RE_CONFIG = design
# )

## ----kalman-model, eval=FALSE-------------------------------------------------
# ou <- nm_model(
#   INPUT = c("ID", "TIME", "DV", "MDV"), ADVAN = 1,
#   PRED = "CL=1; V=1; S1=V; F=0",
#   ERROR = "
#     M0=0; P0=THETA(2)
#     A11=exp(-THETA(1)*DT)
#     Q11=THETA(2)*(1-exp(-2*THETA(1)*DT))
#     H1=1; R1=THETA(3)*THETA(3)
#   ",
#   THETAS = data.frame(THETA=1:3, Value=c(0.3, 0.8, 0.2)),
#   KALMAN_CONFIG = nm_kalman_config(
#     states = "deviation", initial_mean = "M0",
#     initial_covariance = matrix("P0", 1),
#     transition = matrix("A11", 1),
#     process_covariance = matrix("Q11", 1),
#     observation = "H1", observation_variance = "R1",
#     by_dvid = FALSE
#   )
# )
# 
# fit <- nm_est(ou, observations, method = "LAPLACE")
# states <- nm_kalman_decode(fit)
# simulated <- nm_simulate(ou, schedule, residual = TRUE, nsim = 100)

## ----sde-model, eval=FALSE----------------------------------------------------
# sde <- nm_sde_config(
#   states = "deviation", initial_mean = "M0",
#   initial_covariance = matrix("P0", 1),
#   drift = "DRIFT", diffusion = matrix("G", 1),
#   observation = "HX", observation_variance = "R",
#   filter = "ukf", method = "euler", substeps = 8
# )

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

