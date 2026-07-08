## LibeRation 0.4.0

* Shiny client GUI (Pirana-style workspace, ribbon, dark mode), remote jobs via LibeRties,
  VPC/pcVPC, SCM, likelihood profiling, BLQ/censoring, POSTHOC, and expanded test coverage.
* Requires LibeRtAD >= 0.4.0.

## LibeRation 0.1.0

### Scientist-facing features

* BLQ / censoring (Beal M3 and M4): `nm_lik_config(blq_method = "m3"/"m4",
  lloq = )` handles below-limit-of-quantification observations. Censoring is
  taken from a `BLQ`/`CENS` flag or `DV < LLOQ`, with per-row `LLOQ` supported.
  Implemented in the R residual likelihood and in `src/nm_subject_nll.h`;
  estimation automatically uses the numeric R likelihood path when BLQ is on.
* POSTHOC-only estimation: `nm_est(method = "POSTHOC")` (or
  `control$maxeval == 0` / `control$posthoc = TRUE`) fixes the population
  parameters and computes empirical-Bayes ETAs and diagnostics.
* `$TABLE` export: `nm_write_table()` writes NONMEM-style output tables
  (ID TIME DV PRED IPRED RES WRES CWRES ETAs) and `nm_read_table()` reads them
  back.
* Stepwise covariate modelling: `nm_scm()` (forward selection + backward
  elimination on dOFV against chi-square thresholds) with per-step dOFV, AIC,
  BIC, a tidy forest data.frame, and `nm_forest_plot()` (needs ggplot2).
  Trials are warm-started from the parent fit for valid nested dOFV.
* Likelihood profiling: `nm_profile_likelihood()` profiles the OFV over a THETA
  grid (re-optimising the others) and returns a delta-OFV confidence interval.
* VPC / pcVPC: `nm_vpc()` and `nm_pcvpc()` with `strata` and the Bergstrand
  prediction correction; `nm_vpc_plot()` for a ggplot summary.

### Bug fixes

* The timing profiler `profile()` S3 method is renamed to `nm_time_profile()`
  (frees the name for likelihood profiling).
* Non-finite (`NaN`/`Inf`) objective values returned to the optimiser are now
  replaced with a large finite penalty so optimisation degrades gracefully.
* `nm_check_pk_engines()` compares the R and C++ PK engines and warns on
  divergence beyond tolerance.
* pcVPC (B12): when the prediction correction cannot be computed for a
  bin/stratum, the affected points are set to `NA` and a warning is issued
  instead of silently falling back to the uncorrected value.
* Covariate models are now recognised by the C++ prediction checker and routed
  through a covariate-aware likelihood path (model validation and prediction
  symbol detection declare covariates before the static C++ check).

### Infrastructure

* Remote cluster jobs: `nm_remote_server_*`, `nm_job_submit(..., server = )`,
  dataset MD5 verification via LibeRties API.
* Shiny Jobs tab: add/remove remote servers; estimation modal cluster selector.
