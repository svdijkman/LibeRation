# nlmixr2 Ecosystem Analysis for LibeRtAD / LibeRation

**Status:** Analysis only — no implementation yet.  
**Date:** 2025-06-19  
**Constraint:** Learnings only; do not copy nlmixr2 source code. New implementations must preserve current LibeRation functionality (opt-in controls, existing defaults unchanged).

## Source material

Git clone was unavailable in the analysis environment. Review used a local clone at:

`C:\Users\svdijkman.DESKTOP-4OG10M4\Desktop\AD\_nlmixr2_review\`

Packages reviewed:

| Package | Role |
|---------|------|
| **nlmixr2est** | Core estimation (FO/FOCE/FOCEI/Laplace/AGQ/SAEM/nlm/posthoc) |
| **rxode2** | ODE solver, sensitivity ODE generation, event handling |
| **nlmixr2data** | Example datasets |
| **nlmixr2extra** | SCM, profile likelihood, covariate search |
| **nlmixr2plot** | Plotting helpers |
| **lotri** | Block-diagonal Ω DSL |

Dependencies used indirectly via nlmixr2est: **n1qn1**, **lbfgsb3c**, **PreciseSums**. Not cloned separately: **nlmixr2** (meta), **babelmixr2**, **rxode2parse**.

Key nlmixr2est files:

- `src/inner.cpp` — unified FOCE-family likelihood (~7.6k lines)
- `src/saem.cpp` — SAEM C++ engine
- `src/nlm.cpp` — alternate Laplace/NONMEM-like path via `stats::nlm`
- `src/shi21.cpp`, `R/nlmixrGrad.R` — Gill83 / Shi21 finite differences
- `src/cholse.cpp`, `src/nearPD.cpp`, `R/cov.R` — covariance machinery
- `src/cwres.cpp`, `src/npde.cpp`, `src/shrink.cpp` — post-hoc diagnostics
- `R/foceiControl.R`, `R/saemControl.R` — control objects

---

## 1. Architectural comparison

| Aspect | nlmixr2 | LibeRation (current) |
|--------|---------|------------------|
| Model language | rxode2 DSL (compiled ODE + sensitivity) | NONMEM-style `$PK` / `$ERROR` + ADVAN/TRANS |
| Estimation core | Single C++ engine (`inner.cpp`) with algorithm flags | Separate R/C++ paths per method |
| PK engine | rxode2 with auto-generated sensitivity ODEs | Custom C++ solvers + LibeRtAD tapes |
| Parallelism | OpenMP over subjects in C++; rxode2 thread pool | R `parallel::parLapply` for η refits; Laplace `n_threads` stub unused |
| Gradients | Analytical ∂F/∂η via sensitivity ODEs; Gill83/Shi21 FD for outer/cov | LibeRtAD population grad where supported; Shi FD for ∂F/∂η; FOCEI sensitivity grad at fixed η |

**Main nlmixr2 speed advantage:** compile-time sensitivity ODEs + unified C++ inner loop with OpenMP subject parallelism — not a single FOCEI trick.

---

## 2. Estimation methods

### 2.1 FO-family (FO, FOI, FOCE, FOCEI, posthoc)

nlmixr2 uses one engine toggled by `interaction`, `fo`, `nAGQ`, `maxOuterIterations=0` (posthoc).

LibeRation: FO, FOCE, FOCEI implemented; recent additions include `nm_focei_setup()`, warm η, fn/gr cache, sensitivity outer grad, `control$focei_eta = "outer"`.

| Feature | nlmixr2 | LibeRation gap | Priority |
|---------|---------|------------|----------|
| Inner η optimizer | **n1qn1** quasi-Newton with **Hessian warm-start** (curvature carry-over between outer steps) | BFGS/Newton without persistent per-subject curvature | **High** |
| η warm-start | `mceta`: last / zero / compare last+zero+N random Ω draws | `mceta = "last"` only | Medium |
| η reset policy | Z-test on standardized η (`resetEtaP`), nudge ladder on failure | Basic sanitize only | Medium |
| ODE solve skip | `recalc` flag — skip integration if η unchanged | No η-keyed solve memoization | **High** |
| Posthoc-only | `posthocControl()` = FOCEI with 0 outer iterations | No dedicated posthoc method | Low–Medium |
| Generalized LL | `needOptimHess` forces `interaction=0` | Not present | Low |

### 2.2 Laplace & quadrature

- nlmixr2: Laplace = `nAGQ=1`; AGQ = `nAGQ>1` with adaptive GH nodes; zero-point caching for odd quadrature order.
- LibeRation: `LAPLACE` with GH quadrature; `IMP` (importance sampling around Laplace modes).

Learnings:

- Unified Laplace/AGQ flag in one engine reduces code duplication.
- AGQ as first-class method (not only Laplace/IMP).
- Post-mode Laplace correction in inner likelihood (½ log|H| + log|Ω⁻¹|) — verify LibeRation matches NONMEM decomposition.
- nlmixr2 documents extra `+dnorm()` terms may be needed for closer NONMEM OFV on Laplace.

### 2.3 SAEM

nlmixr2 (`saem.cpp`): multi-kernel HMC (multivariate Gibbs, random walk, bootstrap), simulated annealing, Nelder–Mead/Newuoa for complex residual models, burn-in phases (`perNoCor`, `perFixOmega`, `perFixResid`), default **`linFim`** covariance.

LibeRation: SAEM with MH η updates + L-BFGS-B M-step — simpler.

| SAEM feature | Benefit for LibeRation |
|--------------|-------------------|
| Multi-kernel MCMC for η (`nu` vector) | Better mixing on correlated η |
| Simulated annealing during burn-in | Escapes poor basins |
| Correlation-free Ω phase (`perNoCor`) | Stabilizes early SAEM |
| C++ residual-error regression | Complex `$ERROR` without nested R optim |
| SAEM FIM (`fim`) vs `linFim` | Two SE options |
| `adjObf` OFV adjustment | NONMEM bench parity |

### 2.4 Methods comparison matrix

| Method | nlmixr2 | LibeRation |
|--------|---------|--------|
| FO / FOCE / FOCEI | ✓ unified engine | ✓ separate paths |
| Laplace | ✓ (`nAGQ=1`) | ✓ |
| AGQ | ✓ (`nAGQ>1`) | Partial (IMP only) |
| SAEM | ✓ rich C++ | ✓ basic |
| IMP / IMPREC | Not separate | ✓ IMP |
| BAYES / HMC | Limited | ✓ |
| nlm / nlme / nls | ✓ | — |
| Derivative-free outer (bobyqa, newuoa, uobyqa) | ✓ | L-BFGS-B only |
| Posthoc-only | ✓ | — |

---

## 3. Speed improvements (ranked)

### Tier 1 — largest wins

1. **Forward-sensitivity ODEs for ∂F/∂η** (rxode2 pattern)  
   Augment ODE system so sensitivities integrate alongside states. LibeRation currently uses Shi FD or LibeRtAD per observation for G — expensive for ODE models. LibeRtAD could generate sensitivity equations at compile time from `$PK` ODE text.

2. **OpenMP parallel η optimization across subjects**  
   nlmixr2: dynamic OpenMP in `innerOpt()`, per-subject buffers, no R API in workers, serial pre-draw of `mcetaSamples`. LibeRation uses R-level `parLapply` — process overhead. Move parallelism into `nm_fit_all_eta_cpp`.

3. **n1qn1-style inner optimizer with Hessian carry-over**  
   Avoid cold-starting η at every outer θ. Invalidate with sentinel between outer iterations but restore curvature when η is stable. Wrap **n1qn1** package or port minimal interface.

4. **ODE solve memoization (`recalc`)**  
   Skip PK integration when η (and relevant θ) unchanged. Valuable for fn/gr at same θ and line-search revisits.

### Tier 2 — solid wins, lower risk

5. **Gill83 / Shi21 adaptive FD for outer & covariance** — `derivMethod = "switch"` (forward → central near convergence); separate `covDerivMethod`.

6. **Unified `sigdig`-driven tolerance ladder** — one knob for inner/outer/ODE/boundary tolerances.

7. **Subject workload sorting before parallel loop** — balance OpenMP dynamic scheduling for stiff profiles.

8. **Per-subject sticky ODE tolerances** — `indTolRelax`, `maxOdeRecalc`, `odeRecalcFactor` for stiff subjects without global loosening.

### Tier 3 — situational / partially done in LibeRation

9. `focei_eta = "outer"` — LibeRation has this; nlmixr2 does not. Keep optional.

10. fn/gr cache at same θ — LibeRation: `.nm_focei_nested_fetch`.

11. Parameter scaling — LibeRation: optional `par_scale` (off by default).

---

## 4. Gradients & AD

### nlmixr2

- **Inner (η):** rxode2 sensitivity ODEs; FD fallback on event parameters and failed sensitivity (`fallbackFD`, dual `rxInner` / `rxPred` models).
- **Outer (θ, Ω):** Gill83 or Shi21; per-subject `thetaGrad` accumulated during inner eval for S-matrix.
- **Covariance:** separate eps (`hessEps`, `hessEpsLlik`); `covTryHarder` on unscaled space.

### LibeRation

- LibeRtAD tapes for population objective (FO/FOCE where supported).
- FOCEI outer: sensitivity grad at fixed η; full AD through interaction path still broken.
- ∂F/∂η: Shi FD in C++ or numDeriv/LibeRtAD fallback.

### Learnings for LibeRtAD/LibeRation

| Learning | Future action |
|----------|---------------|
| Dual-model FD fallback | Full model for sensitivity + lighter pred model for FD |
| Per-subject outer score during η fit | Avoid expensive post-hoc OPG pass |
| `derivMethod = "switch"` | Forward early, central near convergence |
| Event-parameter FD routing | Extend `.nm_focei_needs_fd_g()` to outer grad |
| Tape reuse across subjects | Cache LibeRtAD tapes per model signature |

---

## 5. Covariance, SE, post-hoc

### nlmixr2

- `covMethod`: `"r,s"` sandwich (default), `"r"`, `"s"`, or skip.
- Robustness: generalized Cholesky (`cholse.cpp`), `nearPD`, auto-fallback, `smatPer` threshold.
- `setCov()`: recalculate covariance without re-optimization.
- Modules: CWRES, IWRES, NPDE, shrinkage, `augPred`.

### LibeRation (current)

- `nm_fit_covariance(type = "hessian"|"linfim"|"sandwich")` — on-demand ✓
- `nm_add_cwres()` ✓
- `nm_shrinkage()` ✓

### Gaps

- NPDE pipeline
- `setCov()` in-place update
- cholSE / nearPD robustness chain
- Automatic R vs S vs sandwich selection
- SAEM-default `linFim` path
- Stage timing (`nlmixrWithTiming`, `benchmarking.R`)

---

## 6. rxode2 / PK engine

| Technique | Relevance |
|-----------|-----------|
| Sensitivity ODE generation (symbolic expansion) | **Critical** for ∂F/∂η speed |
| Compiled C callables (`R_GetCCallable`) | No R eval per ODE step |
| Parallel ODE solves (`setRxThreads`) | Subject/event batching |
| Memory pre-estimation | Avoid realloc in hot loop |
| Stiff solver retry per subject | `maxOdeRecalc`, tolerance relaxation |
| linCmt analytical shortcuts | LibeRation has ADVAN 1–4; extend detection for custom `$PK` |
| Event table (`etTran`, lag, bolus split) | Compare with LibeRation ALAG/infusion FD requirements |

---

## 7. Model & data features

| Feature | nlmixr2 location | LibeRation status |
|---------|------------------|---------------|
| IOV | `iov.R` | Partial (`nm_lik_config`) |
| Mu-referencing | `mu2.R`, rxode2 `mu.R` | Limited |
| Bounded param auto-transform | `preProcessBoundedTransform.R` | Manual bounds |
| lotri Ω DSL | lotri package | `$OMEGA` blocks |
| Censored data (M1/M3) | `censEst.cpp` | — |
| Complex residuals (Box-Cox, etc.) | SAEM C++ | Partial via `nm_lik_config` |
| SCM / covariate search | nlmixr2extra | — |
| Profile likelihood | nlmixr2extra | — |
| Simulation / VPC / NPDE | nlmixr2est | `nm_simulate` + basic VPC |
| NONMEM import | babelmixr2 | `nm_read_nonmem`, bench ✓ |

---

## 8. Control & workflow patterns

| Pattern | Benefit |
|---------|---------|
| Typed control objects (`foceiControl`, `saemControl`) | Prevents invalid combos |
| `getValidNlmixrControl()` S3 validation | Safer defaults |
| Persistent fit environment (`etaMat`, `phiC`, `phiH`) | Fast post-hoc — LibeRation started with `nm_focei_setup()` |
| `sigdig` master tolerance | Easier tuning |
| Benchmark harness + stage timings | Perf regression tests |
| Thread cap (avoid OpenMP nesting) | Prevent 2× oversubscription with R parallel |

---

## 9. Recommended implementation roadmap

All phases: **opt-in controls**, **defaults unchanged**, **existing tests must pass**.

### Phase A — Speed (no algorithm change)

1. OpenMP parallel `nm_fit_all_eta_cpp`
2. ODE / η solve memoization in C++ PK layer
3. Gill83 / Shi21 outer grad with forward/central switch
4. Stage timing + benchmark cases (mirror nlmixr2 `benchmarking.R`)

### Phase B — Algorithm quality (optional modes)

5. n1qn1 inner η optimizer with Hessian warm-start
6. Enhanced `mceta` (multi-candidate warm start)
7. Per-subject sticky ODE tolerance retry
8. Forward-sensitivity ODE generation via LibeRtAD for `$PK` ODE models

### Phase C — Functionality parity

9. NPDE + improved VPC on fit objects
10. `setCov()` / in-place covariance refresh
11. cholSE + nearPD + sandwich auto-fallback
12. AGQ as explicit method; unify with Laplace engine
13. SAEM: multi-kernel MCMC, residual regression, linFim default

### Phase D — Modeling & workflow

14. IOV in C++ likelihood
15. Bounded parameter auto-transform with cov Jacobian
16. Posthoc-only estimation method
17. SCM / profile likelihood (optional separate module)

---

## 10. Adopt carefully (regression risk)

| nlmixr2 behavior | Risk |
|------------------|------|
| Single outer L-BFGS-B pass | LibeRation `max_outer` loop is intentional |
| `adjObf` OFV adjustment | Changes objective constant — bench-only |
| `focei_eta = "outer"` | nlmixr2 doesn't do this; keep LibeRation as explicit approximate mode |
| rxode2-only syntax | NONMEM ctl compatibility must remain primary |
| Nested OpenMP + R parallel | Coordinate thread counts |
| Generalized log-likelihood | Disables interaction; niche |

---

## 11. Gaps vs standard NONMEM

| Area | Notes |
|------|-------|
| Inner η optimizer | nlmixr2 uses n1qn1, not NONMEM Newton |
| FOCEI interaction | Explicit a/B/c matrices; generalized LL disables interaction |
| Laplace setup | INT / YLO / LAPLACE controls differ |
| SAEM SE default | nlmixr2 `linFim` vs NONMEM reporting |
| OFV constant | `adjObf` / `adjLik` adjustments |
| Eta reset | `resetEtaP` vs ETABOUND semantics |
| Parallelism | OpenMP subjects vs NONMEM MPI |
| IMP / IMPREC | LibeRation has IMP; nlmixr2 uses AGQ instead |

---

## 12. Summary

**Highest leverage for LibeRtAD/LibeRation:**

1. Sensitivity ODE generation + C++ parallel η solves + n1qn1 inner optimizer  
2. Robust sandwich covariance (cholSE, nearPD, auto-fallback)  
3. SAEM / MCMC enrichment  

**LibeRation already strong in:** NONMEM-native interface, BAYES/HMC, IMP, NONMEM bench tooling, recent FOCEI speed options (warm η, nested cache, sensitivity grad, `focei_eta` modes).

**Suggested first implementation when ready:** Phase A — OpenMP η refits + ODE recalc cache (largest speed gain, least algorithmic change).

---

## 13. Useful nlmixr2 control reference (for parity discussions)

### `foceiControl()` highlights

- `sigdig`, `interaction`, `outerOpt`, `innerOpt` (n1qn1 vs L-BFGS-B)
- `maxInnerIterations`, `mceta`, `resetEtaP`, `resetHessianAndEta`
- `nAGQ` (0 / 1=Laplace / >1=AGQ)
- `derivMethod`, `covDerivMethod`, `covMethod`, `covTryHarder`
- `shi21maxOuter`, `shi21maxInner`, `fallbackFD`
- `indTolRelax`, `maxOdeRecalc`, `odeRecalcFactor`
- `scaleTo`, `scaleObjective`, `needOptimHess`

### `saemControl()` highlights

- `nBurn`, `nEm`, `nmc`, `nu`, `covMethod` (`linFim` default)
- `adjObf`, `nnodesGq`, `nsdGq`
- `perSa`, `perNoCor`, `perFixOmega`, `perFixResid`

### LibeRation existing FOCEI knobs (for cross-reference)

```r
control = list(
  mceta = "last",           # or "zero", "random"
  maxit_eta_warm = 40,
  maxit_eta = 200,
  focei_grad = "auto",      # or "numeric", "sensitivity"
  focei_eta = "nested",     # or "outer" (approximate fast mode)
  par_scale = FALSE,
  n_cores = 4
)
options(LibeRation.focei_G = "shi")
options(LibeRation.maxit_eta_warm = 40)
options(LibeRation.focei_eta = "nested")
```
