setwd("C:/Users/svdijkman.DESKTOP-4OG10M4/Desktop/AD/LibeRation")
pkgload::load_all()

diag_case <- function(advan, trans, tag) {
  parts <- nm_ctl_template(advan, trans = trans)
  design <- nm_bench_mixed_design(advan, trans, n_per = 1, seed = 12345)
  rcpp <- nm_bench_simulate_rcppnm(parts, design, seed = 12345)$data
  nm_path <- sprintf("inst/nonmem_bench/all_run_v9/%s/simtab_struct", tag)
  nm <- nm_bench_read_simtab(nm_path)
  cmp <- nm_bench_compare_sim(rcpp, nm)
  cat(sprintf("\n=== %s OBSCMP=%d design obs CMT=%d ===\n",
              tag, parts$advan, unique(design$CMT[design$EVID == 0])))
  imp <- nm_read_nonmem(
    sprintf("inst/nonmem_bench/all_run_v9/%s/sim_struct.ctl", tag),
    data_path = "design.csv"
  )
  cat("model OBSCMP:", imp$model$OBSCMP, "USE_ODE:", imp$model$USE_ODE, "\n")
  cat("max_ipred_rel:", cmp$max_ipred_rel, "ok:", cmp$ok, "n_obs:", cmp$n_obs, "\n")
  if (!is.null(cmp$merged) && nrow(cmp$merged) > 0) {
    m <- cmp$merged
    idx <- which.max(abs(m$IPRED_rcpp - m$IPRED_nm) / pmax(abs(m$IPRED_nm), 1e-8))
    cat("worst row TIME:", m$TIME[idx], "CMT:", m$CMT[idx],
        "rcpp:", m$IPRED_rcpp[idx], "nm:", m$IPRED_nm[idx], "\n")
  }
  invisible(cmp)
}

diag_case(6L, 1L, "adv6_trans1")
diag_case(13L, 1L, "adv13_trans1")
diag_case(12L, 4L, "adv12_trans4")
diag_case(2L, 2L, "adv2_trans2")
