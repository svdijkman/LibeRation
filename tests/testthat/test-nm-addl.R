test_that("ADDL/II expands implied doses", {
  skip_if_not_installed("data.table")
  subj <- data.frame(
    TIME = c(0, 12, 24),
    EVID = c(1L, 0L, 0L),
    AMT = c(100, 0, 0),
    ADDL = c(2L, 0L, 0L),
    II = c(12, 0, 0),
    CMT = 1L,
    MDV = c(1L, 0L, 0L),
    stringsAsFactors = FALSE
  )
  ev <- .nm_expand_addl(subj)
  expect_equal(nrow(ev), 5L)
  expect_equal(sort(ev$TIME[ev$EVID == 1L]), c(0, 12, 24))
  expect_true(all(ev$ADDL == 0L))
})

test_that("ADDL/II superposition matches 1-comp analytical", {
  skip_if_not_installed("data.table")
  ke <- 0.1
  v <- 10
  amt <- 100
  ii <- 10
  times <- c(0, 5, 10, 15, 20, 25)
  subj <- data.frame(
    TIME = times,
    EVID = c(1L, rep(0L, length(times) - 1L)),
    AMT = c(amt, rep(0, length(times) - 1L)),
    ADDL = c(2L, rep(0L, length(times) - 1L)),
    II = c(ii, rep(0, length(times) - 1L)),
    RATE = 0,
    CMT = 1L,
    F1 = 1,
    MDV = c(1L, rep(0L, length(times) - 1L)),
    DV = 0,
    SS = 0L,
    stringsAsFactors = FALSE
  )
  subj_ev <- .nm_subject_events(subj)
  ss <- rep(0L, nrow(subj_ev))
  ii_vec <- rep(0, nrow(subj_ev))
  cpp <- nm_pk_route_r(
    1L, 2L, 1L, 1L, 0L, FALSE, 0L,
    subj_ev$TIME, subj_ev$AMT, subj_ev$RATE, subj_ev$F1,
    subj_ev$CMT, subj_ev$EVID, ss, ii_vec,
    list(CL = ke * v, VC = v, V = v)
  )
  dose_times <- c(0, ii, 2 * ii)
  anal <- vapply(subj_ev$TIME, function(t) {
    sum(amt * exp(-ke * (t - dose_times[dose_times <= t]))) / v
  }, numeric(1))
  expect_equal(cpp, anal, tolerance = 1e-10)
})

test_that("ADDL/II increases C++ exposure vs single dose", {
  skip_if_not_installed("data.table")
  times <- c(0, 6, 12, 18, 24, 30)
  pred_obs <- function(addl) {
    subj <- data.frame(
      TIME = times,
      EVID = c(1L, rep(0L, length(times) - 1L)),
      AMT = c(320, rep(0, length(times) - 1L)),
      ADDL = c(addl, rep(0L, length(times) - 1L)),
      II = c(12, rep(0, length(times) - 1L)),
      RATE = 0, CMT = 1L, F1 = 1,
      MDV = c(1L, rep(0L, length(times) - 1L)),
      stringsAsFactors = FALSE
    )
    ev <- .nm_subject_events(subj)
    ss <- rep(0L, nrow(ev))
    ii <- rep(0, nrow(ev))
    ipred <- nm_pk_route_r(
      2L, 2L, 2L, 1L, 0L, FALSE, 0L,
      ev$TIME, ev$AMT, ev$RATE, ev$F1, ev$CMT, ev$EVID, ss, ii,
      list(KA = 1, CL = 3, VC = 20)
    )
    obs <- ev$EVID == 0L
    stats::setNames(as.numeric(ipred[obs]), ev$TIME[obs])
  }
  single <- pred_obs(0L)
  repeat_dose <- pred_obs(1L)
  expect_true(all(repeat_dose >= single - 1e-12))
  expect_true(repeat_dose["30"] > single["30"])
})
