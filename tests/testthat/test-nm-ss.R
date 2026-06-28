test_that("SS initial conditions match SS_sol.R for 1-comp IV bolus", {
  skip_if_not_installed("data.table")
  ss_path <- system.file("SS_initial_cond_sol", "SS_sol.R", package = "LibeRation")
  skip_if_not(nzchar(ss_path), "SS_sol.R not installed")
  source(ss_path, local = TRUE)

  tau <- 12
  cl <- 3
  vc <- 20
  dose <- 100
  times <- c(0, 3, 6, 9, 11.9)
  ref <- ss_1_bolus(dose, tau, times, CL = cl, Vc = vc)

  subj_ev <- data.frame(
    TIME = times,
    EVID = c(1L, rep(0L, length(times) - 1L)),
    AMT = c(dose, rep(0, length(times) - 1L)),
    RATE = 0,
    CMT = 1L,
    F1 = 1,
    SS = c(1L, rep(0L, length(times) - 1L)),
    II = c(tau, rep(0, length(times) - 1L)),
    stringsAsFactors = FALSE
  )
  cpp <- nm_pk_route_r(
    1L, 2L, 1L, 1L, 0L, FALSE, 0L,
    subj_ev$TIME, subj_ev$AMT, subj_ev$RATE, subj_ev$F1,
    subj_ev$CMT, subj_ev$EVID, subj_ev$SS, subj_ev$II,
    list(CL = cl, VC = vc, V = vc)
  )
  expect_equal(as.numeric(cpp), ref$conc, tolerance = 1e-6)
})

test_that("SS initial conditions match SS_sol.R for 2-comp oral", {
  skip_if_not_installed("data.table")
  ss_path <- system.file("SS_initial_cond_sol", "SS_sol.R", package = "LibeRation")
  skip_if_not(nzchar(ss_path), "SS_sol.R not installed")
  source(ss_path, local = TRUE)

  tau <- 24
  times <- seq(0, 24, by = 2)
  ref <- ss_2_oral(
    dose = 320, tau = tau, times = times,
    CL = 3, Vc = 20, Q = 10, Vp = 50, KA = 1, F = 1
  )

  subj_ev <- data.frame(
    TIME = times,
    EVID = c(1L, rep(0L, length(times) - 1L)),
    AMT = c(320, rep(0, length(times) - 1L)),
    RATE = 0,
    CMT = c(1L, rep(2L, length(times) - 1L)),
    F1 = 1,
    SS = c(1L, rep(0L, length(times) - 1L)),
    II = c(tau, rep(0, length(times) - 1L)),
    stringsAsFactors = FALSE
  )
  cpp <- nm_pk_route_r(
    4L, 4L, 2L, 1L, 0L, FALSE, 0L,
    subj_ev$TIME, subj_ev$AMT, subj_ev$RATE, subj_ev$F1,
    subj_ev$CMT, subj_ev$EVID, subj_ev$SS, subj_ev$II,
    list(KA = 1, CL = 3, VC = 20, VP = 50, Q2 = 10)
  )
  expect_equal(as.numeric(cpp), ref$conc, tolerance = 1e-5)
})

test_that("SS trough exceeds single-dose trough for repeated regimen", {
  skip_if_not_installed("data.table")
  times <- c(0, 11.9, 12, 23.9)
  mk <- function(ss_flag) {
    data.frame(
      TIME = times,
      EVID = c(1L, 0L, 1L, 0L),
      AMT = c(100, 0, 100, 0),
      RATE = 0, CMT = 1L, F1 = 1,
      SS = c(ss_flag, 0L, 0L, 0L),
      II = c(12, 0, 12, 0),
      stringsAsFactors = FALSE
    )
  }
  ss0 <- nm_pk_route_r(
    1L, 2L, 1L, 1L, 0L, FALSE, 0L,
    mk(0L)$TIME, mk(0L)$AMT, mk(0L)$RATE, mk(0L)$F1,
    mk(0L)$CMT, mk(0L)$EVID, mk(0L)$SS, mk(0L)$II,
    list(CL = 3, VC = 20)
  )
  ss1 <- nm_pk_route_r(
    1L, 2L, 1L, 1L, 0L, FALSE, 0L,
    mk(1L)$TIME, mk(1L)$AMT, mk(1L)$RATE, mk(1L)$F1,
    mk(1L)$CMT, mk(1L)$EVID, mk(1L)$SS, mk(1L)$II,
    list(CL = 3, VC = 20)
  )
  expect_true(ss1[2] > ss0[2])
  expect_true(ss1[4] > ss0[4])
})
