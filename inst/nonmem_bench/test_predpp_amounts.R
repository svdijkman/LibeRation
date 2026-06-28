setwd("C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation")
devtools::load_all(".", quiet = TRUE)

test_case <- function(advan, trans) {
  wd <- tempfile(sprintf("adv%d_trans%d_", advan, trans))
  dir.create(wd)
  parts <- nm_ctl_template(advan, trans, data_file = "design.csv")
  parts$input_cols <- .nm_bench_nonmem_input_cols(advan, trans)
  parts_struct <- .nm_bench_parts_for_mode(parts, "sim_struct")
  design <- nm_bench_mixed_design(advan, trans, n_per = 1L, seed = 2024L)
  .nm_bench_write_csv(design, file.path(wd, "design.csv"), input_cols = parts$input_cols)
  ctl <- nm_bench_ctl(parts_struct, mode = "sim_struct", sim_seed = 2024L, sim_table = "simtab")
  writeLines(ctl, file.path(wd, "run.ctl"))
  res <- nm_bench_run_nonmem(file.path(wd, "run.ctl"), mod_path = "run", work_dir = wd)
  tab <- nm_bench_read_simtab(file.path(wd, "simtab"))
  cat(sprintf("adv%d_trans%d status=%s simtab=%s cols=%s\n",
              advan, trans, res$status, !is.null(tab),
              if (is.null(tab)) "" else paste(intersect(names(tab), c("A1","A2","A3","A4")), collapse=",")))
  invisible(tab)
}

test_case(1L, 1L)
test_case(2L, 1L)
test_case(2L, 2L)
