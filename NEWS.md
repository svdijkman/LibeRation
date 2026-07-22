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
