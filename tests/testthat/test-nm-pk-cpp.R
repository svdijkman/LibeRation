test_that("compartment scaling Sx divides amount by S(obs)", {
  time <- c(0, 1, 4, 8)
  amt <- c(100, 0, 0, 0)
  f1 <- rep(1, 4)
  rate <- rep(0, 4)
  cmt <- c(1L, 2L, 2L, 2L)
  evid <- c(1L, 0L, 0L, 0L)
  ss <- rep(0L, 4)
  ii <- rep(0, 4)
  pred_scaled <- list(KA = 1, CL = 3, V = 20, S1 = 1, S2 = 20)
  pred_legacy <- list(KA = 1, CL = 3, V = 20, VC = 20, V1 = 20)
  ip_s <- nm_pk_route_r(
    2L, 2L, 2L, 1L, 0L, FALSE, 0L,
    time, amt, rate, f1, cmt, evid, ss, ii, pred_scaled
  )
  ip_l <- nm_pk_route_r(
    2L, 2L, 2L, 1L, 0L, FALSE, 0L,
    time, amt, rate, f1, cmt, evid, ss, ii, pred_legacy
  )
  expect_true(all(is.finite(ip_s)))
  expect_equal(ip_s, ip_l, tolerance = 1e-8)
})

test_that("C++ PK router supports ADVAN/TRANS combinations", {
  skip_if_not_installed("data.table")
  mk_subj <- function(n = 5L) {
    data.frame(
      TIME = seq(0, 24, length.out = n),
      AMT = c(100, rep(0, n - 1L)),
      F1 = rep(1, n),
      RATE = rep(0, n),
      CMT = rep(1L, n),
      EVID = c(1L, rep(0L, n - 1L)),
      MDV = c(1L, rep(0L, n - 1L)),
      DV = c(NA, rep(1, n - 1L)),
      stringsAsFactors = FALSE
    )
  }
  time <- mk_subj()$TIME
  amt <- mk_subj()$AMT
  f1 <- mk_subj()$F1
  rate <- mk_subj()$RATE
  cmt <- mk_subj()$CMT
  evid <- mk_subj()$EVID
  ss <- rep(0L, length(time))
  ii <- rep(0, length(time))

  pred1 <- list(CL = 3, VC = 20, V = 20)
  ip1 <- nm_pk_route_r(1L, 2L, 1L, 1L, 0L, FALSE, 0L,
    time, amt, rate, f1, cmt, evid, ss, ii, pred1)
  expect_true(all(is.finite(ip1)))
  expect_true(length(ip1) == length(time))

  pred2 <- list(KA = 1, CL = 3, VC = 20)
  ip2 <- nm_pk_route_r(2L, 2L, 1L, 1L, 0L, FALSE, 0L,
    time, amt, rate, f1, cmt, evid, ss, ii, pred2)
  expect_true(all(is.finite(ip2)))

  pred3 <- list(CL = 3, VC = 20, VP = 50, Q2 = 10)
  ip3 <- nm_pk_route_r(3L, 2L, 1L, 1L, 0L, FALSE, 0L,
    time, amt, rate, f1, cmt, evid, ss, ii, pred3)
  expect_true(all(is.finite(ip3)))

  pred4 <- list(KA = 1, CL = 3, VC = 20, VP = 50, Q2 = 10)
  ip4 <- nm_pk_route_r(4L, 4L, 1L, 1L, 0L, FALSE, 0L,
    time, amt, rate, f1, cmt, evid, ss, ii, pred4)
  expect_true(all(is.finite(ip4)))

  pred11 <- list(CL = 3, VC = 20, VP = 50, Q2 = 10, VP2 = 80, Q3 = 5)
  ip11 <- nm_pk_route_r(11L, 4L, 1L, 1L, 0L, FALSE, 0L,
    time, amt, rate, f1, cmt, evid, ss, ii, pred11)
  expect_true(all(is.finite(ip11)))

  pred12 <- list(KA = 0.8, CL = 3, VC = 20, VP = 50, Q2 = 10, VP2 = 80, Q3 = 5)
  ip12 <- nm_pk_route_r(12L, 4L, 2L, 1L, 0L, FALSE, 0L,
    time, amt, rate, f1, cmt, evid, ss, ii, pred12)
  expect_true(all(is.finite(ip12)))

  inf_time <- c(0, 2, 4, 8, 24)
  inf_amt <- c(0, 0, 0, 0, 0)
  inf_rate <- c(50, 50, 0, 0, 0)
  inf_cmt <- rep(1L, length(inf_time))
  inf_evid <- c(1L, rep(0L, length(inf_time) - 1L))
  pred_inf <- list(CL = 3, VC = 20, VP = 50, Q2 = 10, S1 = 20, S2 = 50)
  ip_inf <- nm_pk_route_r(3L, 4L, 1L, 1L, 0L, FALSE, 0L,
    inf_time, inf_amt, inf_rate, f1[seq_along(inf_time)], inf_cmt, inf_evid,
    ss[seq_along(inf_time)], ii[seq_along(inf_time)], pred_inf)
  expect_true(all(is.finite(ip_inf)))
  expect_true(ip_inf[2] > ip_inf[1])

  pred_tr <- list(KTR = 1, CL = 3, VC = 20)
  ip_tr <- nm_pk_route_r(4L, 4L, 2L, 1L, 1L, FALSE, 0L,
    time, amt, rate, f1, cmt, evid, ss, ii, pred_tr)
  expect_true(all(is.finite(ip_tr)))

  expect_true(nm_cpp_advan_supported(3L, 4L))
  expect_false(nm_cpp_advan_supported(99L, 2L))
})

test_that("C++ PRED evaluator parses expressions with arithmetic", {
  lines <- c(
    "CL = THETA(1) * exp(ETA(1))",
    "VC = THETA(2) + 1",
    "Q2 = CL / VC"
  )
  out <- nm_eval_pred_cpp(lines, theta = c(3, 20), eta = 0.1, covariates = list())
  expect_equal(as.numeric(out$Q2), 3 * exp(0.1) / 21, tolerance = 1e-10)
  expect_true(nm_pred_expr_check_cpp(lines))
})

test_that("PRED evaluator skips ODE $DES lines with compartment amounts", {
  pk <- c(
    "CL = THETA(1) * exp(ETA(1))",
    "V2 = THETA(2) * exp(ETA(2))",
    "V3 = THETA(3)",
    "Q2 = THETA(4)",
    "K23 = Q2/V2",
    "K32 = Q2/V3",
    "KA = THETA(5) * exp(ETA(3))"
  )
  des <- c(
    "DADT(1) = -KA*A(1)",
    "DADT(2) = KA*A(1) - (CL/V2)*A(2) - K23*A(2) + K32*A(3)",
    "DADT(3) = K23*A(2) - K32*A(3)",
    "F = A(2)/V2"
  )
  out <- nm_eval_pred_cpp(pk, theta = 1:5, eta = c(0.1, -0.05, 0.02), covariates = list(), des_lines = des)
  expect_equal(as.numeric(out$K23), 4 / (2 * exp(-0.05)), tolerance = 1e-8)
  expect_equal(as.numeric(out$CL), exp(0.1), tolerance = 1e-8)
})

test_that("ADVAN 6 oral 2-compartment ODE route produces finite predictions", {
  skip_if_not_installed("data.table")
  mk_subj <- function(n = 8L) {
    data.frame(
      TIME = seq(0, 24, length.out = n),
      AMT = c(320, rep(0, n - 1L)),
      F1 = rep(1, n),
      RATE = rep(0, n),
      CMT = rep(1L, n),
      EVID = c(1L, rep(0L, n - 1L)),
      MDV = c(1L, rep(0L, n - 1L)),
      DV = c(NA, rep(5, n - 1L)),
      stringsAsFactors = FALSE
    )
  }
  subj <- mk_subj()
  pred <- list(
    KA = 1.2, CL = 3, V2 = 20, V3 = 50, Q2 = 10,
    K23 = 10 / 20, K32 = 10 / 50
  )
  ip <- nm_pk_route_r(
    6L, 1L, 2L, 1L, 0L, TRUE, 0L,
    subj$TIME, subj$AMT, subj$RATE, subj$F1, subj$CMT, subj$EVID,
    rep(0L, nrow(subj)), rep(0, nrow(subj)), pred
  )
  expect_true(all(is.finite(ip)))
  expect_true(any(ip[-1] > 0))
})

test_that("C++ PRED evaluator parses THEO-style assignments", {
  lines <- c(
    "CL = THETA(1) * exp(ETA(1))",
    "VC = THETA(2) * exp(ETA(2))",
    "VP = THETA(3)",
    "Q2 = THETA(4)",
    "KA = THETA(5) * exp(ETA(3))"
  )
  out <- nm_eval_pred_cpp(lines, theta = 1:5, eta = c(0.1, -0.05, 0.02), covariates = list())
  expect_equal(as.numeric(out$CL), 1 * exp(0.1), tolerance = 1e-10)
  expect_equal(as.numeric(out$VP), 3)
})

test_that("compartment scaling supports S1..S10 in C++ PK", {
  subj <- data.frame(
    TIME = c(0, 1, 2),
    AMT = c(100, 0, 0),
    F1 = 1, RATE = 0, CMT = 1L, EVID = c(1L, 0L, 0L),
    MDV = c(1L, 0L, 0L), stringsAsFactors = FALSE
  )
  pred <- list(CL = 2, V1 = 20, K10 = 0.1)
  pred$S5 <- 25
  ip <- nm_pk_route_r(
    1L, 1L, 1L, 1L, 0L, FALSE, 0L,
    subj$TIME, subj$AMT, subj$RATE, subj$F1, subj$CMT, subj$EVID,
    rep(0L, 3L), rep(0, 3L), pred,
    s1 = numeric(0), s2 = numeric(0), s3 = numeric(0), s4 = numeric(0),
    scale_mat = matrix(numeric(0), 0L, 0L), use_data_scale = FALSE,
    f_mat = matrix(numeric(0), 0L, 0L), use_data_f = FALSE
  )
  expect_true(all(is.finite(ip)))
  out <- nm_eval_pred_cpp(c("S5 = 25", "S1 = 10"), theta = numeric(0), eta = numeric(0), covariates = list())
  expect_equal(as.numeric(out$S5), 25)
})

test_that("compartment dose scaling uses Fx by CMT", {
  subj <- data.frame(
    TIME = c(0, 2),
    AMT = c(100, 0),
    F1 = 1, F2 = 1, RATE = 0,
    CMT = 2L,
    EVID = c(1L, 0L),
    MDV = c(1L, 0L),
    stringsAsFactors = FALSE
  )
  pred <- list(CL = 2, VC = 20, VP = 40, Q2 = 4, K10 = 0.1, K12 = 0.2, K21 = 0.1, F2 = 1)
  f_mat <- matrix(c(1, 1), ncol = 2)
  colnames(f_mat) <- c("F1", "F2")
  ip_full <- nm_pk_route_r(
    3L, 4L, 2L, 1L, 0L, FALSE, 0L,
    subj$TIME, subj$AMT, subj$RATE, subj$F1,
    subj$CMT, subj$EVID, rep(0L, 2L), rep(0, 2L), pred,
    s1 = numeric(0), s2 = numeric(0), s3 = numeric(0), s4 = numeric(0),
    scale_mat = matrix(numeric(0), 0L, 0L), use_data_scale = FALSE,
    f_mat = f_mat, use_data_f = TRUE
  )
  f_mat[, 2] <- 0.5
  pred$F2 <- 0.5
  ip_half <- nm_pk_route_r(
    3L, 4L, 2L, 1L, 0L, FALSE, 0L,
    subj$TIME, subj$AMT, subj$RATE, subj$F1,
    subj$CMT, subj$EVID, rep(0L, 2L), rep(0, 2L), pred,
    s1 = numeric(0), s2 = numeric(0), s3 = numeric(0), s4 = numeric(0),
    scale_mat = matrix(numeric(0), 0L, 0L), use_data_scale = FALSE,
    f_mat = f_mat, use_data_f = TRUE
  )
  expect_equal(ip_half[2] / ip_full[2], 0.5, tolerance = 1e-6)
  out <- nm_eval_pred_cpp(c("F2 = 0.25"), theta = numeric(0), eta = numeric(0), covariates = list())
  expect_equal(as.numeric(out$F2), 0.25)
})

test_that("ODE template emits S1..Sn scaling for ncomp > 2", {
  parts <- nm_ctl_template(6L, data_file = "data.csv", problem = "ODE", ode_ncomp = 4L)
  expect_true(grepl("S4 = 1", parts$pk, fixed = TRUE))
  expect_equal(nm_max_scale_n(), 10L)
})

test_that("nm_pk_route_detail_r returns compartment amounts", {
  pred <- list(KA = 1, CL = 1, V = 10, VC = 10, k20 = 0.1)
  detail <- nm_pk_route_detail_r(
    2L, 2L, 2L, 1L, 0L, FALSE, 0L, 2L,
    c(0, 1, 2), c(100, 0, 0), c(0, 0, 0), c(1, 1, 1),
    c(1, 2, 2), c(1, 0, 0), c(0, 0, 0), c(0, 0, 0), pred
  )
  expect_equal(length(detail$ipred), 3L)
  expect_equal(dim(detail$amounts), c(3L, 2L))
  expect_equal(colnames(detail$amounts), c("A1", "A2"))
  expect_gt(detail$amounts[2L, 2L], 0)
})

test_that("ADVAN4 TRANS3 macro maps to NONMEM microconstants (VSS peripheral)", {
  time <- c(0, 0.5, 1, 2, 4, 8)
  amt <- c(320, 0, 0, 0, 0, 0)
  rate <- rep(0, length(time))
  f1 <- rep(1, length(time))
  cmt <- rep(2L, length(time))
  cmt[1L] <- 1L
  evid <- c(1L, rep(0L, length(time) - 1L))
  ss <- rep(0L, length(time))
  ii <- rep(0, length(time))
  pred_macro <- list(CL = 1, V = 10, V2 = 20, Q = 2, VSS = 30, KA = 1.2, S2 = 10)
  pred_micro <- list(KA = 1.2, K10 = 0.1, K20 = 0.1, K23 = 0.2, K32 = 0.1, S2 = 10)
  ip_macro <- nm_pk_route_r(
    4L, 3L, 2L, 1L, 0L, FALSE, 0L,
    time, amt, rate, f1, cmt, evid, ss, ii, pred_macro
  )
  ip_micro <- nm_pk_route_r(
    4L, 1L, 2L, 1L, 0L, FALSE, 0L,
    time, amt, rate, f1, cmt, evid, ss, ii, pred_micro
  )
  expect_true(all(is.finite(c(ip_macro, ip_micro))))
  expect_equal(ip_macro, ip_micro, tolerance = 1e-6)
})

test_that("ADVAN4 TRANS1 oral micro aliases K10 to k20", {
  time <- c(0, 0.5, 1, 2)
  amt <- c(320, 0, 0, 0)
  rate <- rep(0, length(time))
  f1 <- rep(1, length(time))
  cmt <- rep(2L, length(time))
  cmt[1L] <- 1L
  evid <- c(1L, rep(0L, length(time) - 1L))
  ss <- rep(0L, length(time))
  ii <- rep(0, length(time))
  pred_k10 <- list(KA = 1.2, K10 = 0.1, K23 = 0.5, K32 = 0.3, S2 = 1)
  pred_both <- list(KA = 1.2, K10 = 0.1, K20 = 0.1, K23 = 0.5, K32 = 0.3, S2 = 1)
  ip_k10 <- nm_pk_route_r(
    4L, 1L, 2L, 1L, 0L, FALSE, 0L,
    time, amt, rate, f1, cmt, evid, ss, ii, pred_k10
  )
  ip_both <- nm_pk_route_r(
    4L, 1L, 2L, 1L, 0L, FALSE, 0L,
    time, amt, rate, f1, cmt, evid, ss, ii, pred_both
  )
  expect_equal(ip_k10, ip_both, tolerance = 1e-8)
})

test_that("ADVAN13 ODE steady-state uses CL/V when TRANS=4", {
  time <- c(0, 0.5, 1.5)
  amt <- c(320, 0, 0)
  rate <- rep(0, length(time))
  f1 <- rep(1, length(time))
  cmt <- c(1L, 2L, 2L)
  evid <- c(1L, 0L, 0L)
  ss <- c(1L, 0L, 0L)
  ii <- c(12, 0, 0)
  pred <- list(CL = 1, V = 10, KA = 1, S1 = 1, S2 = 10)
  ip_t4 <- nm_pk_route_r(
    13L, 4L, 2L, 1L, 0L, TRUE, 1L,
    time, amt, rate, f1, cmt, evid, ss, ii, pred
  )
  ip_t1 <- nm_pk_route_r(
    13L, 1L, 2L, 1L, 0L, TRUE, 1L,
    time, amt, rate, f1, cmt, evid, ss, ii, pred
  )
  expect_equal(ip_t4, ip_t1, tolerance = 1e-8)
})

test_that("ADVAN3 TRANS3 IV uses K21 = Q/(VSS - V)", {
  time <- c(0, 1, 2, 4, 8)
  amt <- c(100, 0, 0, 0, 0)
  rate <- rep(0, length(time))
  f1 <- rep(1, length(time))
  cmt <- rep(1L, length(time))
  evid <- c(1L, rep(0L, length(time) - 1L))
  ss <- rep(0L, length(time))
  ii <- rep(0, length(time))
  pred_macro <- list(CL = 1, V = 10, V2 = 20, Q = 2, VSS = 30, S1 = 10)
  pred_micro <- list(K10 = 0.1, K12 = 0.2, K21 = 0.1, S1 = 10)
  ip_macro <- nm_pk_route_r(
    3L, 3L, 1L, 1L, 0L, FALSE, 0L,
    time, amt, rate, f1, cmt, evid, ss, ii, pred_macro
  )
  ip_micro <- nm_pk_route_r(
    3L, 1L, 1L, 1L, 0L, FALSE, 0L,
    time, amt, rate, f1, cmt, evid, ss, ii, pred_micro
  )
  expect_equal(ip_macro, ip_micro, tolerance = 1e-6)
})
