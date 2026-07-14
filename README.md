# LibeRation

LibeRation provides NONMEM-compatible population PK/PD model specification,
simulation, estimation, diagnostics, and a React-based graphical workflow.
The numerical model, event, ADVAN, ODE, likelihood, and automatic-
differentiation paths run in C++; models can be specified with established
LibeRation R syntax or the restricted C++ expression form.

Implemented model paths include ADVAN1-4/11/12, ADVAN6, ADVAN13, arbitrary
linear matrix propagation, infusions, analytical and nonlinear periodic
steady state, correlated OMEGA, IOV, priors, mixtures, BLQ likelihoods, and
time-varying covariates. Estimation methods include FO, FOCE, FOCEI, Laplace,
ITS, IMP, SAEM, and Bayesian estimation. Diagnostics and uncertainty workflows
include covariance, GOF/CWRES, VPC/NPDE/NPC, categorical and time-to-event
VPCs, bootstrap, profile likelihood, and SCM.

## Quick start

```r
library(LibeRation)

model <- nm_model(
  INPUT = c("ID", "TIME", "EVID", "AMT", "CMT", "DV", "MDV"),
  ADVAN = 1,
  PRED = "CL=THETA(1)*exp(ETA(1))\nV=THETA(2)\nS1=V",
  ERROR = "Y=F*(1+ERR(1))",
  THETAS = data.frame(THETA = 1:2, Value = c(2, 20)),
  OMEGAS = data.frame(OMEGA = 1, Value = 0.1),
  SIGMAS = data.frame(SIGMA = 1, Value = 0.1)
)

data <- data.frame(
  ID = 1, TIME = c(0, 1, 2, 4), EVID = c(1, 0, 0, 0),
  AMT = c(100, 0, 0, 0), CMT = 1,
  DV = c(NA, 4.5, 4.1, 3.4), MDV = c(1, 0, 0, 0)
)

fit <- nm_est(model, data, method = "FOCEI")
summary(fit)
liber_gui(model, data)
```

Install LibeRtAD first, then install LibeRation with R 4.1 or newer and a
C++17 toolchain. Install LibeRties as well to enable persistent local and
remote job queues.

LibeRation is MIT licensed. The remaining engineering work is tracked in
[TODO.md](TODO.md).
