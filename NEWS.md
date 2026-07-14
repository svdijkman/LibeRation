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
