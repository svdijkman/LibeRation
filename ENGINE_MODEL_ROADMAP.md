# LibeRation model-engine roadmap

This roadmap covers model classes that cannot be made reliable merely by
generating a compiled row likelihood or an editable ODE template. Each phase
adds a reusable engine abstraction rather than a one-off model implementation.

## Design constraints

- Preserve `nm_model()` and editable `$PK/$PRED`, `$DES`, and `$ERROR` code.
- Keep propagation, likelihood, differentiation, conditional modes, and
  population objectives inside C++.
- Make new engines serializable for local/remote queues and deterministic under
  a saved seed.
- Expose exact/approximate method labels, numerical diagnostics, and validation
  provenance in fitted objects and the GUI.
- Retain current ADVAN and first-class outcome performance when new facilities
  are not used.

## Phase 1 — general sequence and covariance layer

**Implementation status (0.7.0)**

- [x] Arbitrary-state `Q(t) -> exp(Q * DT)` transitions for observed and hidden
  continuous-time Markov models, with filtering, smoothing, Viterbi, simulation,
  and exact CppAD gradients.
- [x] Full cross-endpoint residual correlation blocks at coincident DVID
  observation times, including FO, exact conditional likelihood, and simulation.
- [x] Estimated AR(1) parameters with bounded differentiable transforms and
  independent histories for interleaved endpoints.
- [x] General ARMA(p,q) residual processes through an exact generated
  state-space declaration. Custom temporal kernels remain expressible through
  `KALMAN_CONFIG`.
- [x] General nested/crossed random-effect designs with exact connected
  components and repeated block OMEGA structures. Sparse factorization for
  exceptionally large connected components remains a scalability refinement.

**Scope**

- General continuous-time multi-state models with an arbitrary generator
  matrix `Q(t)` and differentiable matrix-exponential transitions.
- Full cross-endpoint SIGMA blocks at coincident observation times.
- Estimated AR(1), ARMA(p,q), and user-defined residual correlation kernels.
- Study/site/subject/occasion and arbitrary nested or crossed random-effect
  indices with sparse block covariance.

**Engine work**

1. Replace the fixed subject/occasion ETA index calculation with a serializable
   random-effect design table (`level`, `unit`, `ETA block`, `row membership`).
2. Add grouped observation likelihood calls so several DVID records can share
   a multivariate residual contribution.
3. Add differentiable correlation parameter transforms and sparse Cholesky
   operations to the population objective.
4. Generalize the existing HMM transition evaluator to a `Q -> exp(Q*DT)`
   path using the matrix engine already used by linear PK propagation.

**Acceptance gates**

- Exact-gradient checks against finite differences for every covariance block.
- Agreement with analytic two-state CTMC results and independent multi-state
  reference calculations.
- Recovery and coverage simulations for nested/crossed effects and correlated
  endpoints.
- No measurable regression for ordinary diagonal OMEGA/SIGMA models.

## Phase 2 — state-space inference

**Implementation status (0.7.0)**

- [x] Linear Gaussian Kalman filtering in the compiled population objective.
- [x] Rauch--Tung--Striebel retrospective smoothing and state diagnostics.
- [x] Reproducible process/observation simulation with time-varying matrices.
- [x] Extended and unscented Kalman filters with CppAD parameter derivatives.
- [x] Seeded bootstrap particle filtering, genealogical smoothing, and a
  particle-marginal likelihood usable by the population estimators.

**Scope**

- Linear Gaussian Kalman filtering and Rauch–Tung–Striebel smoothing.
- Extended and unscented Kalman filters for nonlinear process/measurement
  models.
- Bootstrap and auxiliary particle filters, particle smoothing, and PMMH or
  particle-marginal likelihood estimation.
- Discrete-time process noise and continuous-discrete stochastic models.

**Engine work**

1. Define a process-model interface separate from the observation outcome:
   transition, process covariance, observation, and observation covariance.
2. Add per-subject latent-state caches and batched filter evaluation in C++.
3. Differentiate deterministic filters directly; use reparameterized or score
   gradients for particle methods with common random numbers.
4. Add filter degeneracy, effective sample size, smoothing, and prediction
   diagnostics to the result schema and GUI.

**Acceptance gates**

- Kalman likelihoods match closed-form/reference implementations.
- EKF/UKF recovery tested across nonlinear PK/PD fixtures.
- Particle estimates converge to the Kalman solution in linear Gaussian cases.
- Repeated seeded runs are bitwise reproducible on one platform and
  statistically reproducible across platforms.

## Phase 3 — stochastic differential equations

**Implementation status (0.7.0)**

- [x] Continuous-discrete Itô SDE declarations on the nonlinear state-space
  engine.
- [x] Fixed-step Euler--Maruyama and diagonal Milstein propagation.
- [x] EKF/UKF moment likelihoods, bootstrap particle likelihoods, and seeded
  stochastic simulation.
- [ ] Adaptive and higher-order SDE solvers remain future refinements.

**Scope**

- Itô SDEs with Euler–Maruyama and Milstein propagation.
- Continuous-discrete likelihood via EKF/UKF and particle filters.
- Optional higher-order or adaptive solvers where differentiability and
  reproducibility can be guaranteed.

**Engine work**

- Add `$DIFFUSION` alongside `$DES`, with a typed drift/diffusion IR.
- Record deterministic approximation paths and random-number streams as dynamic
  tape inputs; retape when adaptive paths change.
- Separate process variability from ETA/IOV and residual SIGMA in parameter
  reporting and covariance calculations.

**Acceptance gates**

- Ornstein–Uhlenbeck likelihood and moment agreement.
- Strong/weak convergence tests under time-step refinement.
- Parameter-recovery and uncertainty coverage simulations.

## Phase 4 — delayed, semi-Markov, and richer latent-state models

**Implementation status (0.7.0)**

- [x] Hidden semi-Markov models with explicit discrete dwell distributions,
  forward filtering, retrospective smoothing, and Viterbi decoding.
- [x] Fixed-step method-of-steps DDEs with parameterized delays,
  differentiable history interpolation, and automatic tape-path guards.
- [x] Exact factorial HMM enumeration with joint decoding and per-chain
  filtered/smoothed/Viterbi marginal output.
- [x] Joint discrete-regime/nonlinear-continuous switching state-space models
  with compiled particle likelihoods and genealogical regime smoothing.

**Scope**

- Delay differential equations with fixed or parameterized delays.
- Hidden semi-Markov models with explicit dwell-time distributions.
- Coupled/factorial HMMs and switching state-space models.
- General recurrent-event intensity models with history-dependent covariates.

**Engine work**

- A history-buffer abstraction shared by DDEs, semi-Markov durations, and
  recurrent-event history features.
- Forward algorithms over duration/state or factored state spaces, with sparse
  pruning and exact log-domain normalization.
- Retaping guards for history interpolation and state-space pruning decisions.

**Acceptance gates**

- DDE convergence against high-accuracy reference solvers.
- Exact enumeration agreement for small HSMM/factorial-HMM fixtures.
- Simulation-estimation recovery for duration and switching parameters.

## Phase 5 — large mechanistic and hybrid models

**Implementation status (0.7.0 experimental)**

- [x] Semi-explicit index-1 DAEs with Newton solves inside the CppAD objective.
- [x] User-declared residual/variable sparsity decomposes algebraic systems into
  independent solve blocks for both simulation and differentiation.
- [x] Stoichiometric QSP reaction-network builder with named species, dosing,
  observation, and optional algebraic constraints.
- [x] Immutable, hashed offline dense-network, linear-spline, and
  Gaussian-process components for `$PK/$PRED` or state-dependent `$DES` use.
- [ ] Large monolithic sparse factorizations, PDE/spatial models, and
  agent-based simulation remain separate future research interfaces.

**Scope**

- DAE support for conservation-constrained QSP models.
- Sparse PBPK/QSP Jacobians, scalable sensitivity methods, and modular organ or
  pathway components.
- Optional Gaussian-process, neural-ODE, or learned surrogate components.
- PDE/spatial or agent-based models only through a separate experimental
  interface; they should not complicate the core population PK/PD engine.

**Engine work**

- Sparse linear algebra and sparse automatic-differentiation Jacobians.
- Implicit DAE stepping with consistent initialization.
- A versioned external-component ABI with explicit offline/network and
  reproducibility boundaries for hybrid models.

**Acceptance gates**

- Mass-balance and conservation tests.
- Scaling benchmarks over state and parameter dimension.
- Gradient agreement and end-to-end parameter recovery.
- Hybrid components are opt-in and do not alter ordinary package dependencies.

## Recommended implementation order

1. Phase 1 general sequence/covariance layer.
2. Linear Kalman filter/smoother, then EKF/UKF.
3. Particle filter/smoother foundation.
4. SDE support on the state-space foundation.
5. DDE and HSMM history layer.
6. Sparse DAE/PBPK/QSP scaling.
7. Experimental hybrid-learning interface.

This order maximizes reuse: the grouped likelihood and random-effect design
from phase 1 support every later phase, while particle infrastructure supports
nonlinear state-space, SDE, switching, and history-dependent models.
