## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")


## ----gq-template, eval=FALSE--------------------------------------------------
# gq_fit <- nm_est(
#   model, data, method = "GQ",
#   gq_grid = "auto", gq_order = 5, gq_level = 3,
#   covariance = TRUE
# )


## ----estimation-template, eval=FALSE------------------------------------------
# fit <- nm_est(
#   model,
#   data,
#   method = "FOCEI",
#   maxit = 200,
#   eta_maxit = 100,
#   tolerance = 1e-6,
#   optimizer_backend = "auto",
#   n_cores = 4,
#   print_every = 10,
#   covariance = TRUE,
#   covariance_type = "auto"
# )


## ----estimation-sequence, eval=FALSE------------------------------------------
# fit <- nm_est_sequence(
#   model, data,
#   stages = list(
#     nm_est_stage("FOCE", maxit = 100, eta_maxit = 50),
#     nm_est_stage("SAEM", n_iter = 500, seed = 20260719)
#   )
# )
# 
# fit$method_sequence
# fit$stages


## ----covariance-template, eval=FALSE------------------------------------------
# covariance <- nm_cov_step(
#   fit,
#   type = "auto",
#   tolerance = 1e-8,
#   samples = 200,
#   seed = 20260715
# )


## ----gof-template, eval=FALSE-------------------------------------------------
# gof <- nm_gof(fit)
# subset(gof, EVID == 0, c(ID, TIME, DV, PRED, IPRED, CWRES))


## ----predictive-template, eval=FALSE------------------------------------------
# vpc <- nm_vpc(
#   fit,
#   nsim = 500,
#   probs = c(0.05, 0.5, 0.95),
#   level = 0.90,
#   stratify = "SEX",
#   seed = 20260715
# )
# 
# npde <- nm_npde(fit, nsim = 500, seed = 20260715)
# npc <- nm_npc(fit, nsim = 500, seed = 20260715)


## ----special-vpc-template, eval=FALSE-----------------------------------------
# categorical_vpc <- nm_vpc_categorical(
#   categorical_fit, outcome = "DV", nsim = 500
# )
# 
# tte_vpc <- nm_vpc_tte(
#   tte_fit, event = "EVENT", nsim = 500
# )


## ----uncertainty-template, eval=FALSE-----------------------------------------
# bootstrap <- nm_bootstrap(
#   fit,
#   n = 1000,
#   seed = 20260715,
#   n_cores = 8
# )
# 
# profile <- nm_profile(
#   fit,
#   parameters = c("THETA1", "OMEGA1"),
#   points = 11,
#   span = 3,
#   level = 0.95
# )


## ----scm-template, eval=FALSE-------------------------------------------------
# candidates <- data.frame(
#   parameter = c("CL", "V"),
#   covariate = c("WT", "SEX"),
#   form = c("power", "categorical"),
#   reference = c(70, NA),
#   category = c(NA, 1),
#   initial = c(0.75, 0.1)
# )
# 
# scm <- nm_scm(
#   fit,
#   candidates,
#   direction = "both",
#   p_forward = 0.05,
#   p_backward = 0.01,
#   max_steps = 20
# )

