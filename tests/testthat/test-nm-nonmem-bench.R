test_that("nm_ctl_compose includes PREDPP SUBROUTINE", {
  parts <- nm_ctl_template(2L, 2L, data_file = "data.csv")
  ctl <- nm_ctl_compose(parts)
  expect_true(grepl("\\$SUBROUTINE ADVAN2 TRANS2", ctl))
})

test_that("nm_bench_mixed_design covers single, multi, and SS regimens", {
  d <- nm_bench_mixed_design(2L, 2L, n_per = 2L)
  expect_true(all(c("ID", "TIME", "EVID", "SS", "II") %in% names(d)))
  ids <- unique(d$ID)
  expect_equal(length(ids), 6L)
  ss_ids <- unique(d$ID[d$SS == 1L & d$EVID == 1L])
  expect_equal(length(ss_ids), 2L)
  multi_ids <- ids[!ids %in% ss_ids][3:4]
  dose_counts <- table(d$EVID[d$ID %in% multi_ids & d$EVID == 1L])
  expect_true(all(dose_counts >= 2L))
  reg <- .nm_bench_regimen_by_id(d)
  expect_equal(unname(reg[as.character(ids[1:2])]), rep("single", 2L))
  expect_equal(unname(reg[as.character(ids[3:4])]), rep("multiple", 2L))
  expect_equal(unname(reg[as.character(ss_ids)]), rep("steady_state", 2L))
})

test_that("nm_bench_compare_sim reports regimen-specific metrics", {
  design <- nm_bench_mixed_design(2L, 2L, n_per = 1L)
  rcpp <- data.frame(
    ID = rep(1:3, each = 2L),
    TIME = rep(c(0.5, 1), 3L),
    CMT = 2L,
    MDV = 0L,
    EVID = 0L,
    IPRED = c(10, 9, 10, 9, 8, 7),
    A1 = c(100, 90, 100, 90, 80, 70),
    stringsAsFactors = FALSE
  )
  nm <- data.frame(
    ID = rep(1:3, each = 2L),
    TIME = rep(c(0.5, 1), 3L),
    CMT = 2L,
    MDV = 0L,
    EVID = 0L,
    PRED = c(10, 9, 10, 9, 8, 7),
    A1 = c(100, 90, 100, 90, 80, 70),
    stringsAsFactors = FALSE
  )
  cmp <- nm_bench_compare_sim(rcpp, nm, design = design)
  expect_true(cmp$ok)
  expect_equal(cmp$by_regimen$single$ok, TRUE)
  expect_equal(cmp$by_regimen$multiple$ok, TRUE)
  expect_equal(cmp$by_regimen$steady_state$ok, TRUE)
  expect_equal(cmp$by_regimen$single$n_obs, 2L)
  expect_equal(cmp$by_regimen$multiple$n_obs, 2L)
  expect_equal(cmp$by_regimen$steady_state$n_obs, 2L)

  rcpp_bad <- rcpp
  rcpp_bad$IPRED[rcpp_bad$ID == 3L] <- rcpp_bad$IPRED[rcpp_bad$ID == 3L] * 2
  cmp_bad <- nm_bench_compare_sim(rcpp_bad, nm, design = design)
  expect_false(cmp_bad$ok)
  expect_true(cmp_bad$by_regimen$single$ok)
  expect_true(cmp_bad$by_regimen$multiple$ok)
  expect_false(cmp_bad$by_regimen$steady_state$ok)
})

skip_nonmem <- function() {
  if (!nm_nonmem_available()) {
    skip("NONMEM (nmfe73) not on PATH")
  }
}

# --- NONMEM-free fidelity benchmarks (one per estimation method) ------------

test_that("FO marginal benchmark: C++ matches reference R implementation", {
  skip_if_not_installed("data.table")
  chk <- nm_bench_fo_marginal_check(id = "warf", n_sub = 8L, seed = 3L)
  expect_true(is.finite(chk$cpp))
  expect_true(is.finite(chk$r))
  expect_true(chk$ok, info = sprintf("cpp=%.6f r=%.6f diff=%.2e", chk$cpp, chk$r, chk$diff))
})

test_that("FOCEI interaction benchmark: engaged for prop, no-op for add", {
  skip_if_not_installed("data.table")
  chk <- nm_bench_focei_interaction_check(n_sub = 20L, seed = 11L)
  expect_true(chk$prop_engaged,
    info = sprintf("prop on=%.4f off=%.4f", chk$prop_on, chk$prop_off))
  expect_true(chk$add_noop,
    info = sprintf("add on=%.6f off=%.6f", chk$add_on, chk$add_off))
  expect_true(chk$r_cpp_ok,
    info = sprintf("R vs cpp FOCEI diff=%.2e", chk$r_cpp_diff))
})

test_that("FO benchmark: finite marginal OFV and structural theta recovery", {
  skip_if_not_installed("data.table")
  rec <- nm_bench_method_recovery(
    id = "warf", method = "FO", n_sub = 20L, seed = 5L, rtol = 0.6
  )
  expect_true(rec$obj_finite)
  expect_true(rec$theta_ok,
    info = sprintf("max theta rel err = %.3f", rec$max_theta_rel))
})

test_that("LAPLACE benchmark: theta recovery and FOCEI/LAPLACE agreement", {
  skip_if_not_installed("data.table")
  rec <- nm_bench_method_recovery(
    id = "iv1", method = "LAPLACE", n_sub = 20L, seed = 5L, rtol = 0.6,
    n_quad = 5L
  )
  expect_true(rec$obj_finite)
  expect_true(rec$theta_ok,
    info = sprintf("max theta rel err = %.3f", rec$max_theta_rel))

  cross <- nm_bench_cross_method(
    id = "iv1", methods = c("FOCEI", "LAPLACE"), n_sub = 20L, seed = 5L,
    theta_rtol = 0.15, max_outer = 4L, n_quad = 5L
  )
  expect_true(all(is.finite(cross$objective)))
  expect_true(cross$theta_agree,
    info = sprintf("FOCEI vs LAPLACE max theta rel = %.3f", cross$max_theta_rel))
})

test_that("SAEM benchmark: finite objective and structural theta recovery", {
  skip_if_not_installed("data.table")
  rec <- nm_bench_method_recovery(
    id = "warf", method = "SAEM", n_sub = 20L, seed = 5L, rtol = 0.6,
    n_iter = 60L, n_burn = 20L
  )
  expect_true(rec$obj_finite)
  expect_true(rec$theta_ok,
    info = sprintf("max theta rel err = %.3f", rec$max_theta_rel))
})

test_that("IMP benchmark: finite objective and structural theta recovery", {
  skip_if_not_installed("data.table")
  rec <- nm_bench_method_recovery(
    id = "warf", method = "IMP", n_sub = 20L, seed = 5L, rtol = 0.6,
    n_imp = 30L, n_quad = 5L
  )
  expect_true(rec$obj_finite)
  expect_true(rec$theta_ok,
    info = sprintf("max theta rel err = %.3f", rec$max_theta_rel))
})

test_that("BAYES benchmark: posterior-mean structural theta recovery", {
  skip_if_not_installed("data.table")
  rec <- nm_bench_method_recovery(
    id = "warf", method = "BAYES", n_sub = 20L, seed = 5L, rtol = 0.6,
    n_burn = 60L, n_sample = 120L
  )
  expect_true(rec$obj_finite)
  expect_true(rec$theta_ok,
    info = sprintf("max theta rel err = %.3f", rec$max_theta_rel))
})

test_that("nm_bench_pilot LibeRation arm runs without NONMEM", {
  skip_if_not_installed("data.table")
  res <- nm_bench_pilot(
    advan = 2L,
    trans = 2L,
    n_per = 2L,
    seed = 99L,
    method = "FOCE",
    run_nonmem = FALSE
  )
  expect_true(is.finite(res$rcpp_fit$objective))
  expect_true(nrow(res$data) > 20L)
  expect_true(file.exists(res$est_ctl_path))
})

test_that("pilot ADVAN2 TRANS2 benchmark vs NONMEM FOCEI", {
  skip_if_not_installed("data.table")
  skip_nonmem()
  res <- nm_bench_pilot(
    advan = 2L,
    trans = 2L,
    n_per = 3L,
    seed = 2024L,
    method = "FOCEI",
    run_nonmem = TRUE
  )
  expect_true(!is.null(res$rcpp_fit))
  expect_equal(res$rcpp_fit$method, "FOCEI")
  expect_true(is.finite(res$rcpp_fit$objective))
  if (is.null(res$nm_ext) || !is.finite(res$nm_ext$obj %||% NA_real_)) {
    if (!is.null(res$nm_run) && isFALSE(res$nm_run$license_ok)) {
      skip("NONMEM license expired or invalid")
    }
    skip("NONMEM estimation did not complete (check gfortran/NONMEM setup)")
  }
  cmp <- res$compare
  expect_false(is.null(cmp))
  expect_true(cmp$ok, info = paste(
    "theta rel diff:", paste(round(cmp$d_theta, 3), collapse = ", "),
    "| obj nm:", cmp$obj_nm, "rcpp:", cmp$obj_rcpp
  ))
})
