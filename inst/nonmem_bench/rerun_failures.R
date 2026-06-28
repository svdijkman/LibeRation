#!/usr/bin/env Rscript
setwd("C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation")
pkgload::load_all(compile = TRUE)

root <- "inst/nonmem_bench/all_run_v10"
tags <- c("adv2_trans2", "adv10_trans1", "adv12_trans4")
summary_path <- file.path(root, "summary.csv")
summary <- utils::read.csv(summary_path, stringsAsFactors = FALSE)

for (tag in tags) {
  m <- regexec("^adv(\\d+)_trans(\\d+)$", tag)
  g <- regmatches(tag, m)[[1]]
  adv <- as.integer(g[2])
  tr <- as.integer(g[3])
  message("Re-running ", tag)
  r <- nm_bench_case(
    advan = adv, trans = tr,
    work_dir = file.path(root, tag),
    n_per = 2L, seed = 2024L,
    run_nonmem = TRUE
  )
  row <- nm_bench_summarize(list(r))$combined
  summary[summary$tag == tag, ] <- row
}

utils::write.csv(summary, summary_path, row.names = FALSE)
cat("\nUpdated summary:\n")
print(summary[summary$tag %in% tags, c("tag", "sim_ok", "max_ipred_rel", "message")])
cat("\nSim pass:", sum(summary$sim_ok %in% TRUE, na.rm = TRUE), "/", nrow(summary), "\n")
cat("Overall pass:", sum(summary$ok, na.rm = TRUE), "/", nrow(summary), "\n")
