# LibeRation

LibeRation provides NONMEM-compatible population PK/PD model specification,
simulation, estimation, diagnostics, and a React-based graphical workflow.
The numerical model, event, ADVAN, ODE, likelihood, and automatic-
differentiation paths run in C++; models can be specified with established
LibeRation R syntax or the restricted C++ expression form.

Noncompartmental analysis is available through `nm_nca()`. Its independent C++
engine supports linear or linear-up/log-down integration, terminal-slope and
partial-AUC analysis; `engine = "ncar"` and `validate = TRUE` provide an optional
reference path through `ncar`/`NonCompart`.

LibeRation 0.9 is a research beta. Install the exact
[ecosystem compatibility set](../docs/INSTALL.md), run `liber_doctor()`, and
use `liber_support_matrix("LibeRation")` to distinguish externally validated,
internally verified, and experimental workflows.

Implemented model paths include ADVAN1-14: analytical compartment kernels,
arbitrary linear matrix propagation (ADVAN5/7), explicit and stiff general
ODEs (ADVAN6/8/9/13/14), Michaelis--Menten elimination (ADVAN10), and
equilibrium DAE constraints (ADVAN9), plus infusions and analytical/nonlinear periodic
steady state, correlated OMEGA, IOV, priors, mixtures, BLQ likelihoods,
explicit NONMEM-compatible MU referencing,
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
Affine MU equations are exploited directly by SAEM and BAYES fixed-effect
blocks; IMP uses the same mapping to re-centre conditional-mode proposals.
The classifier falls back to the conventional estimator for nonlinear or
non-identifiable MU structures and records the reason with the fit.
Diagnostics and uncertainty workflows
include covariance, GOF/CWRES, family-specific outcome residuals/scores,
VPC/NPDE/NPC, multicategory, count, time-to-event, recurrent-event, and
competing-risk VPCs, bootstrap, profile likelihood, and SCM.
The uncertainty and comparison layer additionally provides sampling importance
resampling, stratified subject/cluster and parametric bootstrap, Bayesian
posterior predictive checks, WAIC, PSIS-LOO, and durable fingerprinted model
comparisons with explicit nested-model rules.

Paired NONMEM reference validation currently covers ADVAN1-13. ADVAN14 is
internally verified against its stiff-system reference but remains explicitly
pending external comparison with NONMEM 7.4 or newer.

The model editor exposes three definition routes:

- **ADVAN/PREDPP (`$PK`)** for the conventional compartment workflow.
- **Direct prediction (`$PRED`)** for row-wise models that assign `F` without
  dose-event propagation.
- **ADVAN + prediction layer (`$PK + $PRED`)**, a LibeRation extension in
  which `$PK` and ADVAN/`$DES` run first and a post-ADVAN `$PRED` transforms
  `F_ADVAN`, `A(i)`, model assignments, and row covariates into the final `F`.

All three routes remain inside the C++/CppAD objective. `$ERROR` continues to
define residual variability or a user likelihood. The combined route can be
exported to NONMEM by folding its marked prediction layer into `$ERROR`;
the marker lets LibeRation recover the two editable sources on re-import.

The React workbench includes an explicit visual general-ODE model builder for
linear and nonlinear compartment systems. It generates previewable
`$PK` and `$DES` code with log-normal ETA scaffolding while retaining the
normal editable code windows. A separate drag-and-drop report workflow renders
DOCX/PDF from user text, selected immutable model runs, comparisons, and saved
diagnostics. Fitted hidden Markov models gain a lazy HMM results tab with
filtered, retrospectively smoothed, Viterbi, and combined probability/path
views by subject, sequence, and hidden state. State-space models gain a lazy
States tab with filter/smoother trajectories, uncertainty bands, innovations,
and likelihood contributions. Optional modelling help and report drafting use
either a consented, lazy-loaded WebLLM model in a dedicated browser worker or
an installed Ollama model streamed asynchronously through `ellmer`. Both are
local-only: Ollama is loopback-restricted and unavailable to hosted or remote
browser sessions.

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

Deterministic FO, ITS, FOCE, FOCEI, and Laplace covariance steps use the native
CppAD marginal-objective Hessian, including implicit conditional-mode and
curvature derivatives. `nm_covariance_repair()` is a separate, explicit
recovery/diagnostic API for indefinite covariance matrices; it never changes a
fitted matrix merely because the matrix has a negative eigenvalue.

## Validation

Validation is not limited to continuous PK models. The standalone
[`validation/nonpk`](../validation/nonpk/README.md) campaign checks normalized
categorical/count and event-time likelihoods, observed Markov/CTMC models,
exact HMM/CT-HMM likelihood and decoding, and linear Gaussian state-space
inference against independent mathematical references. Bernoulli, interval-TTE,
and observed Markov likelihoods are additionally compared row-for-row and by
total objective with generated NONMEM 7.3 `LIKELIHOOD` fixtures.

The separate
[`validation/experimental-families`](../validation/experimental-families/README.md)
campaign checks canonical SDE filtering/simulation, smooth-history DDE
method-of-steps and delay sensitivities, nonlinear index-1 DAE constraints,
QSP reaction conservation, and immutable hybrid components against independent
analytic, convergence, finite-difference, metamorphic, and Monte Carlo
references. These remain experimental research interfaces; the validated
claims apply to the named canonical contracts, not every nonlinear system.
The companion
[`validation/edge-families`](../validation/edge-families/README.md) campaign
adds multiplicative/nonlinear SDE simulation, particle convergence,
delayed-dose discontinuities, stiff and larger DDE/DAE/QSP fixtures, compact
QSP recovery, and hybrid numerical/immutability edge cases.

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
