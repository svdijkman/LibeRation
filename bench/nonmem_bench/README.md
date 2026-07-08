# NONMEM vs LibeRation benchmark scripts (development only — not installed with the package)

Cross-validation workflow comparing LibeRation simulation and estimation to NONMEM on identical inputs.

## Workflow (per ADVAN/TRANS pair)

1. **Shared input** — `design.csv` from `nm_bench_mixed_design()` (single-dose, multi-dose, steady-state subjects).
2. **Structural simulation compare** — NONMEM and LibeRation both simulate with negligible IIV/residual; IPRED/PRED compared at observation times.
3. **Stochastic simulation** — NONMEM `$SIM ONLYSIM` fills `data.csv` with DV (full OMEGA/SIGMA).
4. **Re-estimation** — FOCEI in NONMEM and LibeRation on the same `data.csv`; THETA/OMEGA/SIGMA/OBJ compared via `.ext`.

Control streams and datasets are **generated at run time** by the exported API (`nm_bench_case()`, `nm_bench_ctl()`). This folder only contains helper scripts; it is not copied into the installed package.

## Run one pair

```r
devtools::load_all("LibeRation")
LibeRation::nm_bench_case(2L, 2L, n_per = 3L, seed = 2024L)
```

## Run all ADVAN/TRANS pairs

```r
# From a dev checkout:
source("bench/nonmem_bench/run_all.R")
# or call the exported API with a work directory outside the package:
LibeRation::nm_bench_run_all(work_dir = file.path(tempdir(), "nm_bench"), n_per = 2L)
```

Output: one folder per pair under the chosen `work_dir`, plus `summary.csv`.

## Requirements

- `nmfe73` on `PATH` (user-provided NONMEM installation; not bundled)
- Valid NONMEM license
- Set `NONMEM_LIC` or `NM_LIC` to override the license path
- Windows: `-prdefault` and portable GCC prepended to `PATH` automatically when detected

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
