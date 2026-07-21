# LibeRation

LibeRation provides NONMEM-compatible population PK/PD model specification,
simulation, estimation, diagnostics, and a React-based graphical workflow.
The numerical model, event, ADVAN, ODE, likelihood, and automatic-
differentiation paths run in C++; models can be specified with established
LibeRation R syntax or the restricted C++ expression form.

Implemented model paths include ADVAN1-4/11/12, ADVAN6, ADVAN13, arbitrary
linear matrix propagation, infusions, analytical and nonlinear periodic
steady state, correlated OMEGA, IOV, priors, mixtures, BLQ likelihoods,
compiled user-defined likelihoods plus declarative continuous, categorical,
ordinal/IRT, count, event-time, recurrent-event, competing-risk, observed
Markov, two-state continuous-time Markov, and joint DVID outcomes,
finite-state hidden Markov likelihoods with filtering, retrospective
smoothing, Viterbi decoding, hidden semi-Markov dwell distributions,
continuous-time HMMs, exact linear and nonlinear EKF/UKF/particle state-space
models, ARMA residual processes, continuous-discrete SDEs, generalized
nested/crossed random effects, and time-varying covariates. Estimation methods
include FO, FOCE, FOCEI, Laplace,
ITS, adaptive Gaussian quadrature (GQ) with automatic tensor/Smolyak sparse
grids, IMP, SAEM, and Bayesian estimation.
Bayesian workflows include random-walk BAYES, static HMC, and adaptive NUTS;
discrete nonparametric population distributions are available through
fixed-support NPML and adaptive-grid NPAG.
Diagnostics and uncertainty workflows
include covariance, GOF/CWRES, family-specific outcome residuals/scores,
VPC/NPDE/NPC, multicategory, count, time-to-event, recurrent-event, and
competing-risk VPCs, bootstrap, profile likelihood, and SCM.

The React workbench includes an explicit visual ADVAN6/13 model builder for
linear and nonlinear compartment systems. It generates previewable
`$PK/$PRED` and `$DES` code with log-normal ETA scaffolding while retaining the
normal editable code windows. A separate drag-and-drop report workflow renders
DOCX/PDF from user text, selected immutable model runs, comparisons, and saved
diagnostics. Fitted hidden Markov models gain a lazy HMM results tab with
filtered, retrospectively smoothed, Viterbi, and combined probability/path
views by subject, sequence, and hidden state. State-space models gain a lazy
States tab with filter/smoother trajectories, uncertainty bands, innovations,
and likelihood contributions. Optional modelling help and report drafting use a consented,
lazy-loaded WebGPU language model in a dedicated browser worker; inference is
local and the worker's network APIs are disabled after its weights load.

Advanced nonlinear model starts are available with `nm_model_template()` for
nonlinear elimination, transit/dual absorption, parent-metabolite,
effect-compartment, indirect-response, tumour-growth, and TMDD systems. They
are regular editable ADVAN13 models, not a separate execution path.

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

For deterministic marginal integration, `method = "GQ"` automatically uses
tensor Gauss--Hermite quadrature for up to three ETAs and a Smolyak sparse
grid above that. The grid can be selected explicitly with
`gq_grid = "tensor"` or `gq_grid = "smolyak"`; assess convergence by increasing
`gq_order` or `gq_level`, respectively.

`method = "HMC"` and `method = "NUTS"` sample the exact joint CppAD target
and retain divergence, acceptance, R-hat, and effective-sample-size
diagnostics. `method = "NPML"` estimates weights on a fixed ETA support;
`method = "NPAG"` expands and prunes the support grid. Bootstrap is the
recommended uncertainty procedure for the nonparametric methods.

Install LibeRtAD first, then install LibeRation with R 4.1 or newer and a
C++17 toolchain. Install LibeRties as well to enable persistent local and
remote job queues.

## AI-assisted development

GPT-5.6 was used as an AI engineering collaborator to help implement and review
the modelling engines, estimation and simulation workflows, GUI, regression tests, and documentation.
Scientific direction, architecture, validation criteria, and release decisions remain the responsibility of the project owner.

LibeRation is MIT licensed. The remaining engineering work is tracked in
[TODO.md](TODO.md).
Substantial new engine families are planned in
[ENGINE_MODEL_ROADMAP.md](ENGINE_MODEL_ROADMAP.md).
