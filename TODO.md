# LibeR engineering TODO

Last reviewed: 2026-07-23

Legend: completed work is checked and struck through; partially completed work
retains an unchecked parent item with completed sub-items checked.

## 1. Structural tape sharing - partial

- [x] ~~Share prediction tapes when subject event, dose, time, and covariate
  values are identical.~~
- [x] ~~Exclude observation-specific DV/MDV/BLQ/LLOQ fields from the prediction
  structure key.~~
- [x] ~~Represent numeric PRED/DES covariates as CppAD dynamic parameters and
  share a prediction tape across subjects with the same event topology.~~
- [x] ~~Permit one prediction tape to serve heterogeneous covariate values
  without changing active derivatives or event ordering.~~
- [x] ~~Extend dynamic-parameter sharing to analytical Gaussian FO
  observation/objective tapes where it gives a measured benefit.~~
- [ ] Extend objective-tape sharing beyond the eligible FO surface only where
  complex likelihood, mixture, IOV, hierarchical random-effect, ODE and
  residual-correlation equivalence can be proved.
- [ ] Investigate dynamic dose magnitudes and event times only for layouts whose
  event order, dose compartment, ADDL expansion, and infusion topology are
  provably unchanged; retain separate tapes when topology differs.
- [x] ~~Exercise a 1,000-subject heterogeneous-WT pool and verify that all
  subjects use one dynamic prediction tape.~~
- [ ] Add peak-memory and worker-payload measurements for complex studies with
  thousands of subjects and heterogeneous covariates/dosing.

## 2. Batched population kernels - complete

- [x] ~~Evaluate conditional ETA modes in subject batches.~~
- [x] ~~Calculate FOCE/FOCEI/Laplace gradients, mode sensitivities, and
  curvature contributions in compiled population calls.~~
- [x] ~~Use a persistent C++ population objective for single-process outer
  optimization.~~
- [x] ~~Retain R only as the required cross-process coordinator for PSOCK
  execution.~~

## 3. Native outer optimizer - partial / experimental

- [x] ~~Provide a scaled, box-constrained native C++ BFGS optimizer.~~
- [x] ~~Implement bound projection, line search, projected-gradient convergence,
  iteration logging, and cancellation checks.~~
- [x] ~~Benchmark the native optimizer against R L-BFGS-B/BFGS.~~
- [x] ~~Use R's mature optimizer with thin persistent-C++ objective/gradient
  callbacks by default because it required fewer objective evaluations.~~
- [ ] Revisit a fully native L-BFGS-B or trust-region implementation only if it
  can beat the current callback path on the scenario matrix.

## 4. Tape validity and automatic retaping - complete

- [x] ~~Compile model `ifelse()` expressions as CppAD conditional expressions so
  branch outcomes can change without retaping.~~
- [x] ~~Guard ADVAN6/13 tapes against material parameter movement and retape the
  accepted adaptive ODE trajectory.~~
- [x] ~~Report tape records, validity checks, and retapes.~~
- [x] ~~Record and compare matrix-pivot and matrix-exponential scaling paths at
  later parameter points.~~
- [x] ~~Guard adaptive ODE accepted-step and implicit-Newton convergence paths,
  while retaining the parameter-distance guard as a proactive heuristic.~~
- [x] ~~Guard nonlinear steady-state convergence paths.~~
- [x] ~~Reject invalid or extreme conditional-mode trials as retape anchors and
  move incrementally to a finite, pharmacologically valid anchor.~~

## 5. IMP gradients - partial

- [x] ~~Provide the normalized importance-score gradient as the practical
  default.~~
- [x] ~~Provide the exact finite common-random-number objective through
  `imp_gradient = "finite_crn"`.~~
- [ ] Differentiate the finite-sample proposal completely, including mode,
  curvature, Cholesky, and reparameterized-sample derivatives.
- [ ] Benchmark the fully differentiated proposal against the score and
  numerically differentiated finite-CRN options.

## 6. SAEM engine - partial

- [x] ~~Batch subject Metropolis proposals and objective evaluations in C++.~~
- [x] ~~Adapt proposal scale towards a configurable acceptance target.~~
- [x] ~~Use closed-form OMEGA sufficient-statistic updates.~~
- [x] ~~Use closed-form SIGMA updates for additive, proportional, and
  exponential residual models.~~
- [x] ~~Evaluate fixed-ETA objective values and exact gradients in one native
  sweep for R-optimizer SAEM paths.~~
- [x] ~~Retain fixed ETAs, tape references, and reusable point buffers in a
  persistent native context across R-optimizer callbacks.~~
- [x] ~~Run eligible OMEGA and SIGMA sufficient-statistic kernels in C++ under
  both numerical policies.~~
- [x] ~~Provide an explicit native fixed-ETA M-step with value-only trials,
  objective scaling, safeguarded interpolation and damped L-BFGS memory for
  eligible LibeR-optimized analytical models, with fallback telemetry.~~
- [x] ~~Retain objective, acceptance, and proposal-scale traces.~~
- [ ] Move the remaining outer stochastic-approximation loop into C++ only if
  native profiling demonstrates a material benefit without changing the
  seeded stochastic trajectory.
- [ ] Cache/reuse OMEGA factorization within stable parameter states.
- [ ] Add formal stochastic-approximation stationarity and convergence
  diagnostics.

## 7. Covariance and uncertainty - mostly complete

- [x] ~~Implement NONMEM-style R/Hessian information.~~
- [x] ~~Implement S/OPG information.~~
- [x] ~~Implement R-inverse S R-inverse sandwich covariance and robust SEs.~~
- [x] ~~Report eigenvalues, condition numbers, and regularization.~~
- [x] ~~Automatically fall back from R to sandwich or S when appropriate.~~
- [x] ~~Provide profile-likelihood uncertainty with `nm_profile()`.~~
- [ ] Add specialized boundary-aware asymptotics/interval reporting for
  variance parameters.
- [ ] Support profile likelihood for individual elements of correlated OMEGA
  matrices.

## 8. Optimizer diagnostics - partial

- [x] ~~Retain objective and gradient evaluation counts.~~
- [x] ~~Retain outer and conditional-mode iterations/evaluations and convergence
  codes.~~
- [x] ~~Retain tape sharing, validity, record, and retape telemetry.~~
- [x] ~~Support periodic scaled-gradient output in worker logs.~~
- [x] ~~Report initialization, post-initialization model-fit, covariance,
  core fit-plus-covariance, and literal wall time.~~
- [x] ~~Retain objective, projected-gradient, and step-size traces for the native
  optimizer.~~
- [ ] Capture comparable per-iteration step, scaled/unscaled gradient, and
  parameter-scaling telemetry from the default R optimizer path.
- [ ] Retain deeper prediction/ETA/curvature/gradient phase timing as an
  opt-in developer profiler rather than expanding routine user reports.

## 9. Validation and benchmark matrix - partial

- [x] ~~Implement the complete NONMEM 7.6 ADVAN1-18 model surface and GUI
  templates.~~
- [x] ~~Compare predictions with NONMEM for ADVAN1-15 and ADVAN18, including
  ADVAN5/7 general linear, ADVAN8/9 stiff, ADVAN9 equilibrium DAE, and
  ADVAN10 Michaelis--Menten, ADVAN15 equilibrium DAE, and ADVAN18 DDE
  fixtures.~~
- [x] ~~Validate ADVAN16/17 delay dynamics against equivalent NONMEM ADVAN18
  controls and the ADVAN17 equilibrium component against ADVAN15.~~
- [ ] Repeat the direct ADVAN16/17 controls on a NONMEM installation licensed
  for the optional RADAR5NM extension. The current NONMEM 7.6 runtime rejects
  both before evaluation, and this remains `not-run` rather than a pass.
- [x] ~~Validate bolus, infusion, steady state, and modelled infusion
  rate/duration.~~
- [x] ~~Measure paired core and end-to-end NONMEM/LibeRation times.~~
- [x] ~~Provide oral, two-compartment, three-compartment, full-OMEGA, IOV, and
  ODE scenarios.~~
- [x] ~~Provide smoke, quick, and standard dataset profiles with configurable
  subject counts.~~
- [x] ~~Unit-test BLQ, mixtures, IOV, priors, correlated OMEGA, covariates, and
  residual-error variants.~~
- [x] ~~Add independent-reference validation campaigns for non-PK likelihood,
  Markov/HMM/state-space, SDE, DDE, nonlinear DAE, QSP, and offline hybrid
  model contracts.~~
- [x] ~~Extend the experimental-family campaign to named
  multiplicative/nonlinear SDEs, particle convergence, event-sensitive DDE
  delays, stiff and larger DDE/DAE/QSP fixtures, compact QSP recovery, and
  hybrid numerical/immutability edges.~~
- [ ] Add broad simulation-estimation coverage studies and scaling campaigns
  over substantially larger state, parameter, subject, and event dimensions.
- [ ] Extend paired NONMEM timing scenarios to ADVAN4/12, BLQ, mixtures,
  covariates, combined error, and time-varying covariates.
- [ ] Add a scenario-specific NONMEM control stream for paired IOV timing.
- [ ] Add a named very-large/thousands-of-subject benchmark profile with memory
  and startup measurements.

## 10. MU modelling / MU referencing - implemented / monitored

- [x] ~~Add an explicit serializable mapping from each eligible ETA to
  `MU_i = mu_i(theta, subject-level covariates)` and
  `phi_i = MU_i + ETA(i)`.~~
- [x] ~~Accept and round-trip NONMEM `MU_1`, `MU_2`, ... assignments without
  requiring users to rewrite imported control streams.~~
- [x] ~~Validate MU definitions: one ETA mapping, compatible scale, and no
  observation-varying covariates inside MU.~~
- [x] ~~Add symbolic affine/link classification and specific fallback
  diagnostics for nonlinear, aliased, rank-deficient, or otherwise
  ineligible `mu(theta)` relationships.~~
- [x] ~~Use eligible MU structure for bounded generalized least-squares SAEM
  M-steps and Metropolis-corrected Gaussian BAYES blocks; retain the existing
  optimizer or random-walk path otherwise.~~
- [x] ~~Reuse the mapping to preserve individual parameters when warm-starting
  IMP conditional-mode proposals.~~
- [ ] Reuse MU mappings when SCM generates common log-normal covariate models.
- [ ] Add GUI support that can generate MU definitions automatically from
  common log-normal parameter/covariate relationships while leaving expert code
  editable.
- [x] ~~Cross-check algebraic individual-parameter invariance, Gaussian-block
  improvement, safe nonlinear fallback, and estimator-path telemetry against
  equivalent LibeRation fixtures.~~
- [x] ~~Add paired runtime/estimate benchmarks against MU-referenced NONMEM
  FOCEI/IMP/SAEM models to the external validation matrix, retaining only
  scrubbed derived evidence.~~
- [x] ~~Vectorize and cache the MU GLS normal equations and bypass the generic
  SAEM M-step optimizer when all remaining parameters have closed-form
  updates.~~
- [x] ~~Extend the paired matrix to two ETAs, estimated diagonal/full OMEGA,
  estimated SIGMA, and a subject-level covariate coefficient.~~
- [x] ~~Refine IMP score solutions against the exact finite-CRN objective when
  a subject-varying MU design makes the practical score path insufficient.~~

MU modelling should remain optional. Existing expressions such as
`CL = THETA(1) * exp(ETA(1))` are already valid; explicit MU metadata is useful
only when the estimator exploits the additive individual-parameter structure or
when NONMEM control-stream compatibility requires it.

## 11. Advanced CppAD execution - implemented / monitored

- [x] ~~Use multi-direction Forward sweeps for suitable dense Jacobians.~~
- [x] ~~Use subgraph Reverse for sufficiently large sparse Jacobians and report
  the selected derivative strategy and nonzero count.~~
- [x] ~~Measure Hessian structure once, select sparse symmetric coloring only
  for sufficiently large sparse objectives, and cache its pattern/work object.~~
- [x] ~~Persist optimized CppAD graphs with exact version/commit provenance and
  reconstruct worker tapes without parsing or retaping.~~
- [x] ~~Prototype nested-AD-safe `chkpoint_two` ADVAN1 and 2x2 matrix kernels and
  verify values and Jacobians against their direct forms.~~
- [x] ~~Expose ODE recording/evaluation/tape-memory profiles and checkpoint
  decision support across solver tolerances. Current checkpoint prototypes are
  exact but materially slower and larger than direct small-kernel tapes, so
  they remain outside production.~~
- [ ] Promote checkpoint/atomic kernels to the production path only if a
  representative benchmark overcomes their current small-kernel overhead.
- [ ] Revisit sparsity thresholds and cache policy using large population
  models, ODE systems, and remote-worker startup benchmarks.

## 12. Native Bayesian trajectories - complete

- [x] ~~Move HMC leapfrog trajectories and NUTS tree construction into C++ with
  no R callback inside a trajectory.~~
- [x] ~~Keep bounded/native parameter transformations, exact Jacobians, priors,
  full OMEGA Cholesky derivatives, and all subject ETAs on the exact CppAD
  target.~~
- [x] ~~Implement dual-averaged step-size and diagonal mass adaptation with
  cancellation and iteration logging between sampler iterations.~~
- [x] ~~Retain the R sampler as an explicit reference backend and verify native
  target values/gradients against it.~~

## 13. Advanced uncertainty and reproducible comparison - implemented

- [x] ~~Add sampling importance resampling with a heavy-tailed proposal, exact
  population objective, normalized weights, ESS diagnostics, and systematic
  resampling.~~
- [x] ~~Add first-class, serializable model-comparison objects with immutable
  fit/data/model fingerprints, AIC/BIC weights, explicitly declared nested
  LRTs, boundary warnings, and optional parametric-bootstrap calibration.~~
- [x] ~~Use the model-comparison contract in SCM results, LibeRary dual-lane
  extraction reconciliation, the LibeRation comparison GUI, and Help-AI
  evidence.~~
- [x] ~~Add subject-marginal pointwise posterior log likelihood, posterior
  predictive checks, WAIC, and PSIS-LOO with Pareto-k diagnostics.~~
- [x] ~~Add subject- and cluster-level nonparametric bootstrap with optional
  strata, plus fitted-model parametric bootstrap.~~
- [ ] Add adaptive multi-round SIR proposals and proposal sources based on
  bootstrap or sandwich samples.
- [ ] Add BCa bootstrap intervals, explicit occasion-level resampling, and
  durable checkpoint/resume for large distributed uncertainty campaigns.

## 14. Cross-estimator execution architecture - implemented / qualifying

- [x] ~~Fuse eligible analytical ADVAN1--4/11/12 prediction and likelihood
  value sweeps inside the persistent stochastic context, with tape-equivalence
  guards and automatic generic fallback.~~
- [x] ~~Add a bounded native shared-memory subject executor for fused value-only
  work. Keep R-owned data, mutable CppAD tapes, RNG and retaping on the owning R
  thread and retain PSOCK for all ineligible cases.~~
- [x] ~~Use a common conditional-state cache contract for parameter anchors,
  modes and objective states. IMP and adaptive GQ share one implementation;
  deterministic conditional methods and SAEM/BAYES use their persistent C++
  objects because their curvature and seeded-transition semantics differ.~~
- [x] ~~Extend MU-preserving conditional-mode starts to adaptive GQ and the
  optimized compiled FOCE/FOCEI/Laplace objective, with recenter telemetry and
  no change to the NONMEM-compatible path.~~
- [x] ~~Add a unified stochastic diagnostic contract: SAEM
  phase/stationarity evidence and independent replicate summaries, plus
  rank-normalized split R-hat, bulk/tail ESS and MCSE for MCMC methods. Keep
  SAEM automatic stopping opt-in.~~
- [x] ~~Add guarded MU/non-centred geometry for joint HMC/NUTS. Eligible
  optimized models now whiten ETA through OMEGA's Cholesky factor with exact
  Jacobian and chain-rule gradients; incomplete MU, IOV and general-effect
  layouts retain explicit centred fallback.~~
- [x] ~~Run equivalence, fallback and complete standard matched-control
  benchmark gates before enabling any new fast path by default. The 10 August
  2026 results are retained under
  `validation/benchmark/results/20260810-full-cross-method-final/`.~~
