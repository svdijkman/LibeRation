# NONMEM vs LibeRation benchmark library

Cross-validation workflow comparing LibeRation simulation and estimation to NONMEM on identical inputs.

## Workflow (per ADVAN/TRANS pair)

1. **Shared input** — `design.csv` from `nm_bench_mixed_design()` (single-dose, multi-dose, steady-state subjects).
2. **Structural simulation compare** — NONMEM and LibeRation both simulate with negligible IIV/residual (`sim_struct.ctl`); IPRED/PRED compared at observation times.
3. **Stochastic simulation** — NONMEM `$SIM ONLYSIM` fills `data.csv` with DV (full OMEGA/SIGMA).
4. **Re-estimation** — FOCEI in NONMEM and LibeRation on the same `data.csv`; THETA/OMEGA/SIGMA/OBJ compared via `.ext`.

## Run one pair

```r
devtools::load_all("LibeRation")
LibeRation::nm_bench_case(2L, 2L, n_per = 3L, seed = 2024L)
```

## Run all 23 ADVAN/TRANS pairs

```r
# From a dev checkout (scripts are excluded from the installed package):
source("inst/nonmem_bench/run_all.R")
# or call the exported API with a work directory outside the package:
LibeRation::nm_bench_run_all(work_dir = file.path(tempdir(), "nm_bench"), n_per = 2L)
```

Output: one folder per pair under the chosen `work_dir`, plus `summary.csv`.

> **Note:** `inst/nonmem_bench/` is for development only and is not included when
> LibeRation is installed (`R CMD INSTALL`). Do not store benchmark run outputs here;
> use a temp or project folder instead.

## Requirements

- `nmfe73` on `PATH`
- Valid NONMEM license in `nm_7.3.0_g/license/nonmem.lic` as a **single line** (if pasted one character per line, the runner auto-normalizes to a temp file)
- Set `NONMEM_LIC` or `NM_LIC` to override the license path
- Windows: `-prdefault` and portable GCC 4.6 prepended to `PATH` automatically

## API

| Function | Purpose |
|---|---|
| `nm_bench_pairs()` | List all valid ADVAN/TRANS combinations |
| `nm_bench_case()` | Full sim + est pipeline for one pair |
| `nm_bench_run_all()` | Loop all pairs |
| `nm_bench_summarize()` | Aggregate results data frame |
| `nm_bench_compare_sim()` | IPRED/PRED comparison |
| `nm_bench_compare()` | Parameter/OBJ comparison |

## Tests

`tests/testthat/test-nm-nonmem-bench.R`
