# Moved

NONMEM benchmark scripts and documentation are in `bench/nonmem_bench/` at the package root.

Benchmark datasets and control streams are generated at run time by `nm_bench_case()` and related functions. Run outputs should be written to a temp or project directory, not under `inst/`.

If this folder still contains old run artifacts from a previous checkout, you can delete `inst/nonmem_bench/` entirely — it is not used by the installed package.
