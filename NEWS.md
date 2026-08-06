# LibeRation 0.10.5

- Adds a direct SSH scheduler target beside HTTPS and SSH-tunnelled LibeRties.
  Local desktop sessions can submit typed jobs straight to Slurm or Grid
  Engine without an HTTP daemon, including through an optional SSH gateway.
- Persists scheduler identifiers and an encrypted remote queue so status,
  logs, cancellation and completed results can be reconciled after either the
  SSH connection or LibeRation session disappears.
- Extends the existing SSH readiness and key-management workflow to the direct
  scheduler route. Scheduler routing fields are validated and no client shell
  fragment, private key, or queue encryption key is exposed to browser state.

# LibeRation 0.10.4

- Removes institution-specific cluster and gateway host names from the
  SSH-tunnel wizard's initial values, edit fallbacks, placeholders, and server
  messages. New server definitions now start empty and use neutral
  `example.org` guidance.
- Replaces institution-specific SSH test fixtures with generic cluster hosts
  and adds regression assertions that prevent those addresses returning to
  the distributed GUI or server controller.

# LibeRation 0.10.3

- Adds allow-listed `LibeR`, `NONMEM`, and `nlmixr2` execution choices to the
  estimation and simulation dialogs. NONMEM jobs are translated to a control
  stream and run through administrator-configured PsN; nlmixr2 jobs use a
  restricted semantic translation of supported solved-ADVAN and `$DES` models.
- Normalizes supported external-engine estimates, ETAs, tables, timings, and
  covariance output into ordinary LibeRation runs. When audit artifacts are
  requested, original NONMEM outputs are preserved while nlmixr2 equivalents
  are generated with explicit provenance.
- Adds disabled-by-default audit-artifact controls to estimation and simulation
  dialogs. Selected runs generate checksum-verified NONMEM-style listings,
  control streams, estimate/ETA/covariance files where applicable, and output
  tables; unchecked runs retain the existing lightweight result path.
- Materializes optional audit bundles once under each saved run, records their
  provenance in the workspace manifest, verifies them on demand with
  `nm_project_audit_artifacts()`, and removes them with their owning run.
- Extends the remote-server dialog with an SSH readiness assistant. On Windows
  it can install the Microsoft OpenSSH Client through the optional-feature
  mechanism, enable `ssh-agent`, discover or generate protected Ed25519 keys,
  load a key into the agent, install its public half on either SSH hop, and
  test the gateway and destination independently. Passwords, key passphrases,
  MFA challenges, and host-policy decisions remain in visible OpenSSH
  terminals and are never passed through or retained by LibeRation.
- Adds a guided SSH-tunnel connection mode to the remote LibeRties dialog.
  Local desktop sessions can validate and manage a loopback-only OpenSSH port
  forward, including optional jump-host routing, ssh-agent or identity-file
  authentication, connection testing, safe host-key handling and automatic
  session cleanup. Hosted deployments continue to permit direct HTTPS only.
- Applies a selected identity file to both the destination and jump host using
  a restrictive ephemeral OpenSSH configuration, avoiding ProxyJump's failure
  to inherit `-i`. Automatic local-port selection now excludes ports that the
  operating system reserves even though they do not accept connections.
- Journals every queued submission before contacting LibeRties and replays an
  uncertain acknowledgement with a stable, payload-bound idempotency key.
  A connection loss after server acceptance can therefore recover the original
  job identifier without creating a duplicate run.
- Continues automatic reconciliation for completed jobs whose result has not
  yet been downloaded and saved under its model version. Transient result
  transfer failures are retried every 30 seconds and survive GUI/R restarts.

# LibeRation 0.10.1

- Requires LibeRtAD 0.8.1 so native builds inherit the complete CppAD
  R-console adapter and do not link direct `puts`/`putchar` calls.
- Refreshes the embedded beta.12 compatibility contract without rewriting the
  immutable LibeRation 0.10.0 release.

# LibeRation 0.10.0

- Makes the compiled engine pointer read-only and routes every engine input
  through one canonical ADDL materializer. Native entry points now fail closed
  if an unexpanded ADDL record bypasses that boundary.
- Disables silent finite-difference recovery for a non-finite outer gradient by
  default. `allow_fd_gradient = TRUE` is an explicit, warned, and telemetered
  compatibility escape hatch.
- Uses the complete CppAD population-objective Hessian for supported
  deterministic estimators. GQ, IMP, and SAEM covariance now use
  estimator-specific observed marginal score information at the final fit,
  with proposal adaptation held fixed and full provenance; they no longer use
  `optimHess` as their ordinary covariance bread.
- Keeps covariance repair opt-in, records spectral and adjustment diagnostics,
  and never silently changes a materially indefinite ETA covariance matrix.
- Completes the seam-aligned engine source split, adds regression coverage for
  parameter recovery and native-boundary invariants, and publishes the related
  NONMEM and covariance validation evidence.

# LibeRation 0.9.8

- Extends explicit covariance recovery with a symmetric pivoted
  modified-Cholesky option and adds downstream repair-sensitivity reporting.
  Estimated diagonal and correlated OMEGA structures continue to use intrinsic
  log-variance and log-Cholesky parameterizations rather than post-hoc repair.
- Adds `nm_nca()`, an independent native C++ noncompartmental-analysis engine
  with linear and linear-up/log-down integration, automatic or manual terminal
  slope selection, partial AUC, dose-normalised and dosing-interval summaries.
- Adds explicit BLQ and duplicate-time policies, typed provenance, optional
  `ncar`/`NonCompart` reference validation, and an automatic reference fallback.
- Adds an exact CppAD population-objective Hessian for deterministic FO, ITS,
  FOCE, FOCEI, and Laplace covariance calculations. Conditional-mode,
  curvature, log-determinant, transform, and prior derivatives are included;
  the numerical Hessian remains an explicitly labelled comparator/fallback.
- Adds `nm_covariance_repair()` with opt-in round-off clipping, diagonal jitter,
  and Higham nearest-PSD repair plus complete spectral and adjustment
  diagnostics. Material indefiniteness is never repaired silently.
- Adds local parameter-recovery estimation fixtures for FOCEI/Laplace and
  estimation coverage for ADVAN5/7/8/9/10/14. The former monolithic engine is
  partitioned into named event/ADVAN, differential-system, differentiable-
  propagation, likelihood, population, and state-space implementation units,
  with a separately compiled stable population-objective API seam.

# LibeRation 0.9.7

- Adds `nm_model_update()` as a shared, validated model-editing boundary for
  LibeRation and other ecosystem interfaces. Code, parameter tables, bounds,
  covariance structure, model route, and advanced outcome declarations are
  rebuilt through the canonical `nm_model()` validator.
- Sequential estimation now uses the same update boundary when propagating
  fitted THETA, OMEGA, and SIGMA values between stages.

# LibeRation 0.9.6

- Exposes `nm_advan_template()` as the supported public entry point for the
  same ADVAN 1--14 templates used by the New Model Version workflow, allowing
  other LibeR packages to reuse them without duplicating model definitions.
- NONMEM control-stream imports now honour `$MODEL` `DEFDOSE` and
  `DEFOBSERVATION` declarations for ODE and other non-general-linear ADVAN
  models. This is required for safe oral depot/central translations in the
  AEDapt catalogue migration.
- Opens the browser by default when `liber_gui()` is started by a
  non-interactive desktop launcher; hosted callers continue to request and
  receive the application object with `launch.browser = NULL`.

# LibeRation 0.9.5

- Adds explicit MU-reference metadata and editable `MU_n` code generation,
  subject-level covariate validation, NONMEM round-tripping, and model contract
  version 4.
- Specializes MU-aware estimation: symbolic affine/link classification drives
  bounded generalized least-squares SAEM M-steps, IMP preserves individual
  parameters when warm-starting conditional proposals, and BAYES uses a
  Metropolis-corrected Gaussian MU block. Ineligible models fall back safely
  and retain the reason in fit diagnostics.
- Vectorizes and caches the SAEM MU GLS system and skips the generic optimizer
  when MU, OMEGA, and SIGMA closed forms exhaust the M-step. On the paired
  100-subject reference model this reduced specialized SAEM core time by about
  half. IMP now detects subject-varying MU designs and refines its practical
  score solution against the exact finite common-random-number objective.
- Adds a scrubbed external comparison campaign for conventional and
  MU-referenced NONMEM 7.3 FOCEI, IMP, and SAEM controls, including matched
  estimates, individual ETAs, fresh-process/core runtime, and specialization
  telemetry. Scenarios cover one and two ETAs, fixed and estimated diagonal or
  correlated OMEGA, estimated SIGMA, and an estimated weight relationship.
- Adds sampling importance resampling, subject- or cluster-level stratified
  bootstrap, and fitted-model parametric bootstrap.
- Adds durable model-comparison objects with model/data/fit fingerprints,
  AIC/BIC weights, explicitly declared nested likelihood-ratio tests, boundary
  warnings, and optional parametric-bootstrap calibration. SCM, the comparison
  GUI, and Help AI use the same comparison evidence.
- Adds subject-marginal Bayesian pointwise likelihoods, population and
  conditional posterior predictive checks, WAIC, and optional PSIS-LOO with
  Pareto-k diagnostics.

# LibeRation 0.9.4

- Updates GOF, saved diagnostics, comparison, and queue views through targeted
  React messages while leaving the surrounding Shiny interface interactive.
- Automatically refreshes the active GOF view when switching completed runs
  and removes the global faded-busy overlay from routine output updates.
- Recovers an isolated non-finite AD gradient with a bounded, telemetered
  finite-difference step instead of aborting an otherwise valid LAPLACE fit.
- Increments the workbench asset version to invalidate browser caches from
  earlier research-beta installations.

# LibeRation 0.9.3

- Adds three explicit model-definition routes: conventional ADVAN/PREDPP
  `$PK`, row-wise direct `$PRED`, and a LibeRation extension that runs
  `$PK -> ADVAN/$DES -> $PRED -> $ERROR` in one differentiable C++ path.
  The combined layer can consume `F_ADVAN`, `A(i)`, `$PK` assignments, and
  row covariates before assigning the final `F`.
- Preserves separate `$PK` and `$PRED` editor drafts in model contract v3 and
  remote jobs, discovers generated columns from both programs, and folds the
  combined extension into a marked NONMEM `$ERROR` block for round-tripping.
- Completes the ADVAN1-14 model surface: arbitrary linear ADVAN5/7,
  stiff ADVAN8, equilibrium-capable ADVAN9, integrated Michaelis--Menten
  ADVAN10, and stiff/nonstiff ADVAN14 now have native model contracts,
  C++ dispatch, GUI templates, and exact CppAD prediction paths.
- Adds direct generated NONMEM fixtures for ADVAN5/7/8/9/10, including an
  ADVAN9 equilibrium compartment. ADVAN14 remains internally verified because
  the available NONMEM 7.3 installation predates that subroutine.
- Adds a standalone non-PK validation campaign for categorical/count,
  event-time, observed/hidden Markov, continuous-time Markov, and Gaussian
  state-space models. Exact/analytic references are supplemented by paired
  NONMEM 7.3 likelihood fixtures where a faithful counterpart exists.
- Adds a 19-check numerical validation campaign for canonical SDE, DDE,
  nonlinear index-1 DAE, QSP reaction-network, and offline hybrid-component
  contracts, with analytic, convergence, conservation, derivative, and seeded
  Monte Carlo evidence.
- Adds a complementary 27-check edge campaign covering multiplicative and
  nonlinear SDE simulation, particle convergence, delayed bolus boundaries,
  stiff DDEs, larger/coupled DAEs, larger/stiff QSP systems, compact parameter
  recovery, and hybrid numerical/immutability edges.
- Corrects continuous-discrete SDE EKF/UKF covariance propagation so process
  noise introduced during a substep is transported through subsequent drift,
  and retains differentiable DDE history interpolation at stored time nodes.
- Retains left/right DDE event history, splits integration at parameterized
  delayed discontinuities, and propagates the right-limit sensitivity across
  observation records.
- Compiles prediction-scoped hybrid components before user `$PK/$PRED`, uses a
  stable softplus expression for extreme inputs, and rejects modified component
  payloads whose immutable hash no longer matches.
- Imports and round-trips ADVAN5/7 `$MODEL` graphs and emits valid unique
  eight-character compartment identifiers for generated NONMEM streams.

# LibeRation 0.9.2

- Updates the exact ecosystem compatibility manifest for the LibeRties 0.7.4
  durable queue race correction in research beta 3.

# LibeRation 0.9.1

- Updates the ecosystem compatibility manifest for the cross-platform
  LibeRality source-publication correction in research beta 2.

# LibeRation 0.9.0

- Adds a machine-readable ecosystem support matrix with validated, verified,
  and experimental evidence tiers exposed through `liber_support_matrix()`.
- Adds `liber_support_bundle()` and a GUI Support action that create
  inspectable, privacy-preserving diagnostic archives without dataset values,
  estimates, ETAs, environment variables, workspace contents, or model code
  by default.
- Recovers interrupted workspace writes from the durable previous generation,
  with regression coverage for corrupt primary records.
- Adds three runnable onboarding workflows and scheduled 1,000-subject
  benchmarks that record cold/core time, payload sizes, and peak R heap.

# LibeRation 0.8.3

- Restores the established high-resolution LibeR dove and visibly aligns the
  workbench header, controls, panels, and dialog geometry with the ecosystem.
- Adds keyboard-safe dialogs with Escape handling, focus containment, and
  focus restoration, plus consistent focus indicators throughout the GUI.
- Uses the ecosystem-wide theme preference and a transparent package-coloured
  dove favicon.

# LibeRation 0.8.2

- Moves HMC and NUTS trajectory generation, dual averaging, diagonal mass
  adaptation, and posterior transformation into C++, with an explicit R
  reference fallback and target/gradient equivalence tests covering priors and
  full OMEGA Cholesky parameters.
- Adds persistent prediction/objective tape-size telemetry, nonlinear ODE
  tolerance profiling, and measured checkpoint decision support without
  silently changing the production solver.
- Reuses LibeRtAD's cached sparse-Hessian path for sufficiently large sparse
  full objectives while retaining the dense route elsewhere.
- Defines and tests process-isolated parallel execution; PSOCK workers now
  retain one private worker-state object instead of multiple global bindings.

# LibeRation 0.8.1

- Corrects covariance scaling for NONMEM and likelihood Hessian conventions and
  validates covariance standard errors against paired NONMEM runs.
- Extracts the native optimizer into a dedicated compilation unit, enables
  large Windows object files, and adds randomized ADVAN and parameter-bound
  property tests.
- Adds release-qualified capability metadata and reproducible validation
  provenance for analytical ADVAN, steady-state, ODE, FO, and FOCEI paths.

# LibeRation 0.8.0

- Adds a versioned semantic model contract used by LibeRties so advanced HMM,
  state-space, DDE/DAE, random-effect, component, experimental, and outcome
  configuration survives remote job round trips without executable R objects.
- Introduces workspace schema v2 with content-addressed model/data/result
  objects, atomic project locks, legacy migration, integrity verification,
  backup, and dry-run garbage collection. Repeated model versions no longer
  duplicate large datasets and fit payloads.
- Adds `liber_doctor()` and an explicit capability/validation tier to separate
  implemented, integration-tested, reference-validated, and experimental
  model families.
- Adds provenance-bearing numerical validation gates for objectives,
  parameters, variability, ETAs, and predictions. Benchmark summaries are not
  publishable unless correctness tolerances pass before timing is considered.
- Adds gated real-browser regression tests for responsive layout, modal
  reachability, scrolling, and deferred large-data transfer.

# LibeRation 0.7.3

- Keeps the workbench within the available viewport height and assigns vertical
  scrolling to its main content area, so lower sidebar actions remain reachable
  in short desktop windows and on mobile without restoring page gutters.

# LibeRation 0.7.2

- Removes Shiny/Bootstrap's default 15-pixel page gutters so the workbench
  reaches both viewport edges at desktop, laptop, tablet, and mobile widths,
  while preserving the workbench's own responsive scrolling.

# LibeRation 0.7.1

- Restores exact HMM emission-parameter gradients for `$ERROR` expressions
  that branch on fixed observations such as `ifelse(DV == 0, ...)`, via the
  corrected LibeRtAD 0.7.4 conditional-expression implementation.
- Moves the browser-local AI model and context controls into a compact settings
  dialog opened beside the `Activate AI` switch.
- Synchronizes THETA, OMEGA, and SIGMA tables with parameter references found
  in the editable model code during validation and application, while retaining
  explicit add/remove-row controls for manual model construction.
- Presents estimates consistently in pharmacometric order: THETA, OMEGA,
  SIGMA.
- Corrects visual-model parameter renaming so referenced parameter rows are
  renamed instead of duplicated, and permits compartment numbers to be swapped
  within the current diagram.

# LibeRation 0.7.0

- Adds opt-in per-browser-session workspaces for hosted demonstrations and
  moves all reactive workbench state inside the Shiny session boundary.
- Adds an explicitly gated experimental-engine layer with serialized feature,
  strictness, purpose, and provenance metadata retained by local/remote jobs.
- Adds delay differential equations through `nm_dde_config()`. Parameterized
  `LAG(A(i), delay)` expressions use a fixed-step method-of-steps RK4 kernel,
  differentiable linear history interpolation, and the same C++/CppAD path for
  simulation, prediction derivatives, and estimation.
- Adds semi-explicit index-1 DAEs through `$ALG` and `nm_dae_config()`.
  Algebraic Newton solves remain on the AD tape; declared sparsity decomposes
  independent residual/variable blocks.
- Adds `nm_qsp_system()` and `nm_qsp_model()` for named-species stoichiometric
  QSP reaction networks with optional algebraic constraints.
- Adds exact factorial HMMs with joint-state likelihoods and per-chain
  filtered, retrospectively smoothed, and Viterbi output.
- Adds switching nonlinear/SDE state-space models using a stratified,
  differentiably importance-weighted joint regime/continuous-state particle
  likelihood and genealogical regime smoothing.
- Adds immutable offline dense-network, spline, and Gaussian-process
  components. Components can augment `$PK/$PRED` or run inside `$DES` for
  learned state-dependent dynamics without network access at execution time.
- The React model editor shows experimental provenance and solver summaries,
  exposes `$ALG`, and identifies QSP, factorial, switching, and hybrid models.

- Adds compiled nonlinear state-space inference. `nm_kalman_config()` now
  selects an extended Kalman filter, unscented Kalman filter, or reproducibly
  seeded bootstrap particle filter in addition to the exact linear filter.
  Decoding supplies EKF/UKF RTS smoothing or genealogical particle smoothing.
- Adds continuous-discrete Itô stochastic differential equations through
  `nm_sde_config()`, with Euler--Maruyama and diagonal Milstein propagation,
  EKF/UKF moment likelihoods, particle likelihoods, and seeded simulation.
- Adds `nm_re_block()` and `nm_re_config()` for exact nested or crossed
  site/study/subject/reader-style random-effect designs. Independent connected
  components, repeated within-block OMEGA structures, structural subject
  propagation, and conditional-mode dimensions are generated automatically.
- Adds `nm_arma_config()` as an exact state-space declaration for general
  ARMA(p,q) residual processes around the structural prediction.
- Adds `nm_hsmm_config()` for explicit discrete dwell-time distributions.
  Sparse duration-state expansion reuses the compiled HMM forward,
  retrospective smoothing, Viterbi, and exact-gradient paths while decoded
  output is aggregated back to the original clinical states.
- Adds arbitrary-state continuous-time hidden Markov models through
  `nm_cthmm_config()`. The C++/CppAD engine constructs each generator diagonal
  from non-negative off-diagonal rates and evaluates `exp(Q * DT)` at irregular
  observation times. Filtering, retrospective smoothing, Viterbi decoding,
  and exact rate gradients share the same transition implementation.
- Generalizes first-class continuous-time Markov outcomes beyond two states.
  Observed states use deterministic emissions on the same matrix-exponential
  sequence engine, while simulation draws from the general transition matrix.
- Adds `nm_residual_group()` for full cross-endpoint residual correlation at
  coincident DVID observations. Fixed or transformed THETA/SIGMA correlations
  are supported in exact conditional objectives, FO marginal objectives, and
  stochastic simulation.
- AR(1) correlation may now be estimated from a THETA or SIGMA with a safe
  hyperbolic-tangent transform. Interleaved DVID histories are tracked
  independently in both the exact and FO paths.
- Adds a differentiable linear Gaussian state-space engine through
  `nm_kalman_config()`. It provides exact C++ Kalman likelihoods, process and
  observation simulation, `nm_kalman_decode()` filtering, and retrospective
  Rauch--Tung--Striebel smoothing with irregular/time-varying matrices defined
  in editable `$ERROR` code.
- Adds `nm_outcome()`/`nm_outcomes()` as a declarative outcome layer. Normal,
  log-normal, fixed-df Student t, Bernoulli, multinomial/ordinal, Poisson,
  negative-binomial, binomial, ZIP/hurdle, first/recurrent event,
  competing-risk, observed Markov, and exact two-state continuous-time Markov
  endpoints generate editable compiled `$ERROR` likelihoods and share one
  C++/CppAD population objective. Joint endpoints are selected by `DVID`.
- First-class outcomes now provide stochastic outcome generation and
  `nm_outcome_diagnostics()` with expected values, conditional variance,
  observed-category probabilities, Pearson/deviance residuals, Brier/log
  scores, hazard, cumulative hazard, and martingale residuals as applicable.
- Generalizes categorical VPCs beyond binary outcomes and adds count,
  recurrent-event, and competing-risk VPCs. Competing risks use
  Aalen--Johansen cumulative-incidence curves. These diagnostics persist as
  first-class project run metadata and can be selected in DOCX/PDF reports.
- Adds `nm_irt_outcomes()` for multi-item ordinal/IRT endpoint declarations.
- Adds eight editable advanced ADVAN13 templates through
  `nm_model_template()`: nonlinear elimination, transit/dual absorption,
  parent-metabolite, effect compartment, indirect response, tumour growth,
  and full TMDD. The GUI exposes these under New version from template.
- Adds tape-safe piecewise and restricted-cubic-spline expression generators.
- Documents the staged C++ engine plan for general multi-state/covariance,
  state-space, SDE, DDE/semi-Markov, sparse DAE/QSP, and hybrid models in
  `ENGINE_MODEL_ROADMAP.md`.

# LibeRation 0.6.9

- Extends `nm_hmm_decode()` with a scaled retrospective forward-backward
  smoother and log-domain Viterbi decoding. `method = "all"` returns explicit
  filtered, smoothed, and Viterbi state columns plus per-sequence likelihood
  and Viterbi log-posterior summaries.
- Adds a lazy HMM results tab to the React workbench. Filtered, retrospective
  smoothed, Viterbi, and combined views provide subject/sequence and state
  selectors, probability trajectories, classified paths, and per-sequence
  likelihood evidence without loading HMM rows until the tab is opened.
- Adds finite-state hidden Markov models through `nm_hmm_config()`. Initial,
  transition, and state-conditional emission expressions remain on the C++
  CppAD population-objective tape and are combined with a numerically scaled
  forward algorithm independently by subject and, optionally, `DVID`.
- Adds `nm_hmm_decode()` for record-level filtered state probabilities and
  classifications. HMM columns are also included in `nm_gof()` while
  Gaussian residual diagnostics are correctly reported as undefined.
- Added compiled user-defined observation likelihoods. `$ERROR` may assign a
  positive `LIK` probability/density or a `LOGLIK` contribution; the complete
  likelihood remains in C++ and on the CppAD population-objective tape.
- Added `PREV_DV`, `PREV_TIME`, `DT`, and `FIRST` helpers for first-order
  Markov models, tracked separately by subject and `DVID`. Finite MDV baseline
  outcomes seed the state without contributing to the objective.
- Direct `F` assignments in `$PK/$PRED` now override compartment-derived
  predictions, enabling general direct-prediction categorical/likelihood
  models. The GUI detects likelihood models, defaults them to LAPLACE, and
  hides incompatible Gaussian FO/FOCE/FOCEI choices.

- Added independent, persistent Help and Report context-window controls with
  model-aware Auto defaults, 1K--16K presets, and a guarded custom setting.
- Help now retains the complete conversation for dynamic budgeting, allowing
  an 8K Help context to preserve substantially more useful history than the
  previous fixed three-turn/4K limit.
- WebLLM now receives the selected context at model-load time. GPU-memory
  allocation and device-reset failures retry with smaller contexts, while the
  UI reports the approximate prompt budget, retained messages, output reserve,
  and any compaction.

# LibeRation 0.6.8

- AI-authored report sections now load every saved run selected anywhere in
  the visual workflow and receive the selected evidence types, including
  estimates and uncertainty, model/data metadata, timings, compact GOF
  statistics, covariance status, diagnostics, and model code.
- Report drafting now asks the local model to synthesize a connected account
  across selected runs and avoids the previous erroneous "no selected run
  evidence" fallback and generic missing-facts checklist.
- Added a persistent report save-location field beside the filename plus a
  native folder chooser, with manual path entry retained for platforms where a
  native chooser is unavailable.

# LibeRation 0.6.7

- Help AI now selects compact project-index, model-code, or result-detail
  evidence according to the question instead of attaching every available
  payload to every prompt.
- Added model-aware prompt budgeting and shorter conversation retention so
  browser-local requests stay inside the selected model's context window.
- Added a tokenizer-error fallback in the WebGPU worker that retries once with
  a smaller evidence payload and otherwise reports an actionable error.
- Project-index requests no longer deserialize complete saved fit objects,
  making simple project and run-availability questions faster.

# LibeRation 0.6.6

- Added an explicit high-performance WebGPU-adapter preflight before the local
  AI runtime or model artifacts are loaded.
- Distinguishes an unavailable browser GPU adapter from a missing physical GPU
  and explains browser/driver recovery while keeping the non-AI GUI usable.
- Added a `Reset local AI` control and clears stale AI errors when a new
  workbench page mounts, preventing an earlier browser GPU failure from looking
  like a LibeRation startup failure.

# LibeRation 0.6.5

- Help AI now loads compact summaries of the 20 most recent completed
  estimation and simulation runs in the selected project on demand.
  Objectives, parameter estimates and uncertainty, timings, convergence,
  output-column names, and
  saved diagnostic availability are supplied without loading row-level result
  data into the browser.
- Added a request/response context handshake so the first question can wait for
  saved project evidence and continue automatically; selecting a run is no
  longer required before asking about existing results.
- Recognises Windows D3D12 `DXGI_ERROR_DEVICE_REMOVED`, adapter, command-queue,
  and related WebGPU device-loss failures. A failed GPU worker is replaced once
  and the original request is retried from cached model artifacts, with a clear
  recovery message if the replacement also fails.
- Project/run selection changes now interrupt an in-flight answer without
  abruptly terminating the resident GPU worker, and model switches briefly
  wait for the previous GPU session to unload.

# LibeRation 0.6.4

- Added the selected project, model version, model run, and compact dataset
  metadata to the browser-local Help context. Help history is reset when that
  selection changes so answers from one project cannot bleed into another.
- Hardened the WebGPU worker lifecycle. Disposed/device-lost GPU sessions are
  rebuilt once from the browser cache before a request fails, worker crashes
  reject pending requests, and failed generations no longer remain visually
  stuck at `Generating...`.

# LibeRation 0.6.3

- Compiles against Eigen supplied directly by LibeRtAD and removes the
  RcppEigen build dependency without changing the numerical matrix backend.
- Uses LibeRtAD's controlled CppAD--Eigen compatibility header and explicit
  dense R/Eigen conversion helpers.

# LibeRation 0.6.2

- Added separate persistent Help and Report browser-local LLM selectors in the
  header and relevant panels. New workspaces use Qwen 2.5 Coder 3B for Help and
  Qwen 2.5 7B for reports; existing single-model settings migrate to Help. Only
  one lazy WebGPU worker remains resident, and cached models are switched on
  demand. Reports can alternatively use `Same as Help model`.
- Expanded the selectors from minimal models through Qwen 2.5 7B, Llama 3.1
  8B, and Gemma 2 9B choices with approximate memory and use-case labels.
- Tightened local-assistant and report-drafting evidence rules and sampling so
  missing project facts are reported as unavailable instead of being inferred.
- Added per-row deletion to the visual model builder's structural-parameter
  table. Parameters referenced by a compartment or flow are protected with an
  explanatory tooltip until that reference is removed or renamed.
- Requires LibeRtAD 0.7.2 so the complete bundled CppAD public-header tree is
  present when compiling LibeRation from source.
- Added a drag-and-drop visual structural-model builder for general nonlinear
  ADVAN6/13 systems. Compartment and flow semantics generate a previewable
  `$PK/$PRED`, `$DES`, error model, THETA, and log-normal ETA/OMEGA scaffold;
  applying is explicit and `$DES` always remains manually editable.
- Added opt-in browser-local WebGPU assistance in a dedicated lazy worker for
  modelling help and report drafting. Activation/consent/model choice persist,
  weights load only on first use, and worker network APIs are disabled before
  any model/run context is supplied.
- Added a linear drag-and-drop report workflow with user or local-AI narrative
  blocks, immutable model-run evidence, model comparisons and saved diagnostic
  plots. Workflows persist with projects and render DOCX/PDF plus a provenance
  manifest.
- Added ordered multi-stage estimation (`nm_est_sequence()` and
  `nm_est_stage()`), including population-parameter hand-off, compatible ETA
  warm starts, stage telemetry, local/remote queue execution, and GUI controls.
- Added static discovery and explicit selection of model-generated output
  columns. Selected PRED/IPRED/CWRES, ETA, compartment, and `$PK`/`$PRED`
  assignment columns are retained with fitted runs and loaded lazily into Data
  explorer.
- Validation now compiles the unsaved editor draft and refreshes its available
  output catalogue. NONMEM `$TABLE` columns round-trip through the same model
  output selection.
- Clarified full covariance labels as `OMEGA(row,col)` throughout the editor;
  ETA remains the random-effect vector governed by OMEGA.
- Added adaptive generalized Gaussian quadrature (`method = "GQ"`) with C++
  Gauss--Hermite grids from LibeRtAD, fixed-node reference integration,
  batched CppAD scores, covariance support, point-count safeguards, and GUI
  controls. GQ now also supports signed-weight Smolyak sparse grids and
  automatically selects them above three ETA dimensions.

- Compiles against the CppAD 20260000.0 headers supplied directly by LibeRtAD;
  RcppEigenAD and BH are no longer build dependencies.
- Shares prediction tapes across subjects with identical event topology while
  updating heterogeneous numeric PRED/DES covariates as CppAD dynamic
  parameters.
- Selects multi-direction Forward or sparse subgraph-Reverse prediction
  Jacobians according to graph dimensions and records strategy telemetry.
- Detects changed matrix pivots, matrix-exponential regimes, adaptive ODE
  trajectories, and steady-state convergence paths and retapes automatically.
- Bounds automatic conditional-mode retaping to finite, pharmacologically
  valid ETA anchors so extreme line-search trials cannot become tape anchors.
- Reuses FO likelihood tapes across structurally equivalent subjects by moving
  observations and covariates to CppAD dynamic parameters, and fuses eligible
  analytical FO objectives into one population tape with a safe subject-tape
  fallback.
- Uses a method-aware FO evaluator profile and a stricter exact-gradient
  convergence test, avoiding unused tape construction and premature OMEGA
  convergence under large rescaled objectives.

# LibeRation 0.6.0

- Rebuilt the numerical runtime around LibeRtAD and compiled C++ population
  objectives, gradients, event processing, and specialized ADVAN kernels.
- Added ADVAN1-4/11/12, arbitrary linear matrix propagation, ADVAN6 RK45,
  ADVAN13 implicit integration, infusions, and steady-state handling.
- Added FO, FOCE, FOCEI, Laplace, ITS, IMP, SAEM, and Bayesian estimation with
  bounds, priors, covariance diagnostics, parallel execution, and detailed
  run telemetry.
- Added bootstrap, profile likelihood, SCM, VPC/NPDE/NPC, categorical and
  time-to-event VPCs, CWRES GOF, and NONMEM control-stream round-tripping.
- Rebuilt the React workbench with named model versions, nested numbered runs,
  persistent queues/settings, lazy data and diagnostic loading, comparison
  views, syntax highlighting, and light/dark themes.
- Added persistent-C++ callbacks to R's mature optimizers and batched
  population kernels for substantially lower callback and compilation costs.

This release is an architectural and API break from the 0.4.x series.
