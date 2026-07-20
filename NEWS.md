# LibeRation 0.6.4

- Added the selected project, model version, model run, and compact dataset
  metadata to the browser-local Help context. Help history is reset when that
  selection changes so answers from one project cannot bleed into another.
- Hardened the WebGPU worker lifecycle. Disposed/device-lost GPU sessions are
  rebuilt once from the browser cache before a request fails, worker crashes
  reject pending requests, and failed generations no longer remain visually
  stuck at `Generating...`.

# LibeRation 0.6.3

- Compiles against Eigen 3.4.0 supplied directly by LibeRtAD and removes the
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
