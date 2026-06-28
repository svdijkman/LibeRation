setwd("C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation")
pkgload::load_all(compile = TRUE)

spot <- function(advan, trans, tag) {
  parts <- nm_ctl_template(advan, trans = trans)
  design <- nm_bench_mixed_design(advan, trans, n_per = 1, seed = 2024)
  rcpp <- nm_bench_simulate_rcppnm(parts, design, seed = 2024)$data
  wd <- file.path(tempdir(), paste0("spot_", tag))
  dir.create(wd, showWarnings = FALSE)
  ctl <- nm_bench_ctl(parts, mode = "sim_struct", sim_seed = 2024L, sim_table = "simtab_struct")
  writeLines(ctl, file.path(wd, "sim_struct.ctl"))
  write.csv(design, file.path(wd, "design.csv"), row.names = FALSE)
  nm_run <- nm_bench_run_nonmem(file.path(wd, "sim_struct.ctl"), mod_path = "sim_struct", work_dir = wd)
  nm_tab <- nm_bench_read_simtab(file.path(wd, "simtab_struct"))
  cmp <- if (!is.null(nm_tab)) nm_bench_compare_sim(rcpp, nm_tab) else list(ok = NA, reason = "no nm tab")
  cat(sprintf("%s sim_ok=%s max_ipred=%s max_amt=%s nm_status=%s\n",
              tag, cmp$ok %||% NA, cmp$max_ipred_rel %||% NA, cmp$max_amt_rel %||% NA,
              nm_run$status %||% "?"))
  invisible(cmp)
}

for (x in list(
  c(2, 2, "adv2_trans2"),
  c(3, 3, "adv3_trans3"),
  c(4, 1, "adv4_trans1"),
  c(4, 5, "adv4_trans5"),
  c(10, 1, "adv10_trans1"),
  c(12, 1, "adv12_trans1"),
  c(12, 4, "adv12_trans4"),
  c(12, 6, "adv12_trans6")
)) spot(as.integer(x[1]), as.integer(x[2]), x[3])
