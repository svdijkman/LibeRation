# LibeRation

**NONMEM-style** nonlinear mixed-effects modelling for population PK/PD: FO, FOCE, FOCEI, SAEM, Laplace, IMP, and Bayesian MCMC.

Depends on [LibeRtAD](../LibeRtAD) for automatic differentiation.

## Installation

Install **LibeRtAD** first, then LibeRation:

```r
devtools::install("path/to/LibeRtAD")
devtools::install("path/to/LibeRation")
```

Suggested packages for the Shiny GUI and async jobs: `shiny`, `callr`, `ggplot2`, `DT`.

## Quick start

```r
library(LibeRation)

sim <- nm_synthetic_theo(n_sub = 12, seed = 1)
fit <- nm_est(
  sim$model, sim$data,
  method = "FOCE",
  grad = "cpp",
  pk_engine = "cpp",
  control = list(maxit = 200)
)

predict(fit)
nm_etab(fit)
summary(fit)
```

## Shiny GUI

```r
liberation_shiny()
# options(LibeRation.workspace = "path/to/projects")
```

## Project workspace

```r
nm_workspace_init("~/LibeRation_projects", create_demo_project = TRUE)
nm_workspace_create_project("my_study", template = "theo")
```

## Main function groups

| Area | Functions |
|------|-----------|
| Estimation | `nm_est()`, `nm_focei_setup()`, `nm_cov_step()` |
| Data / model | `nm_dataset()`, `nm_model()`, `nm_read_nonmem()` |
| Control stream | `nm_ctl_parse()`, `nm_ctl_compose()`, `nm_ctl_template()` |
| GOF | `predict()`, `nm_etab()`, `nm_add_cwres()`, `nm_add_npc_npde()` |
| Simulation | `nm_simulate()`, `nm_synthetic_*()` |
| Workspace | `nm_workspace_*()` |
| Jobs | `nm_job_submit()`, `nm_job_result()`, `nm_remote_server_*()` |
| Diagnostics | `nm_scm()`, `nm_vpc()`, `nm_pcvpc()`, `nm_profile_likelihood()` |
| Tables | `nm_write_table()`, `nm_read_table()` |

## Documentation

```r
?LibeRation
?nm_est
?liberation_shiny
```

Regenerate `.Rd` manuals:

```r
roxygen2::roxygenise("path/to/LibeRation")
```

## Vignette

```r
vignette("getting-started", package = "LibeRation")
```

Requires **knitr**, **rmarkdown**, and Pandoc.

## License

MIT — see `LICENSE`.
