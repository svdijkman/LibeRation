setwd("C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation")
devtools::load_all(".", quiet = TRUE)

wd <- tempfile("nm_amt_")
dir.create(wd)
parts <- nm_ctl_template(6L, 1L, data_file = "design.csv")
parts$input_cols <- .nm_bench_nonmem_input_cols(6L, 1L)
parts_struct <- .nm_bench_parts_for_mode(parts, "sim_struct")
design <- nm_bench_mixed_design(6L, 1L, n_per = 1L, seed = 2024L)
.nm_bench_write_csv(design, file.path(wd, "design.csv"), input_cols = parts$input_cols)
ctl <- nm_bench_ctl(parts_struct, mode = "sim_struct", sim_seed = 2024L, sim_table = "simtab")
writeLines(ctl, file.path(wd, "run.ctl"))
cat("TABLE line:\n", tail(strsplit(ctl, "\n")[[1L]], 1L), "\n\n")
res <- nm_bench_run_nonmem(file.path(wd, "run.ctl"), mod_path = "run", work_dir = wd)
cat("NM status:", res$status, "\n")
if (file.exists(file.path(wd, "simtab"))) {
  tab <- nm_bench_read_simtab(file.path(wd, "simtab"))
  cat("Columns:", paste(names(tab), collapse = ", "), "\n")
} else {
  tab <- NULL
  cat("No simtab; log tail:\n")
  cat(paste(tail(res$log, 8L), collapse = "\n"), "\n")
}

rcpp <- nm_bench_simulate_rcppnm(parts_struct, design, seed = 2024L)
if (!is.null(tab)) {
  cmp <- nm_bench_compare_sim(rcpp$data, tab)
  print(cmp[c("ok", "ipred_ok", "amt_ok", "max_ipred_rel", "max_amt_rel")])
  print(cmp$amt_cols)
}
