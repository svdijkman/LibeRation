test_that("candidate normalization fills defaults and validates relationship", {
  cands <- .nm_scm_normalize_candidates(list(
    list(parameter = "CL", covariate = "WT"),
    list(parameter = "V", covariate = "AGE", relationship = "power")
  ))
  expect_length(cands, 2L)
  expect_identical(cands[[1]]$relationship, "exponential")
  expect_identical(cands[[2]]$relationship, "power")
  expect_error(
    .nm_scm_normalize_candidates(list(list(
      parameter = "CL", covariate = "WT", relationship = "bogus"
    ))),
    "arg"
  )
})

test_that("candidate data.frame input is accepted", {
  df <- data.frame(
    parameter = c("CL", "V"),
    covariate = c("WT", "AGE"),
    relationship = c("power", "linear"),
    stringsAsFactors = FALSE
  )
  cands <- .nm_scm_normalize_candidates(df)
  expect_length(cands, 2L)
  expect_identical(cands[[1]]$covariate, "WT")
  expect_identical(cands[[2]]$relationship, "linear")
})

test_that("effect equality and membership helpers behave", {
  a <- list(parameter = "CL", covariate = "WT", relationship = "power")
  b <- list(parameter = "CL", covariate = "WT", relationship = "power")
  d <- list(parameter = "CL", covariate = "AGE", relationship = "power")
  expect_true(.nm_scm_eff_eq(a, b))
  expect_false(.nm_scm_eff_eq(a, d))
  expect_true(.nm_scm_in(a, list(d, b)))
  expect_false(.nm_scm_in(a, list(d)))
})

test_that("effect injection adds a THETA and covariate to the model", {
  sim <- nm_synthetic_iv1(n_sub = 2L, seed = 1L)
  m <- sim$model
  n_th0 <- nrow(m$THETAS)
  m2 <- .nm_scm_add_effect(m, "V", "WT", "power", 70)
  expect_equal(nrow(m2$THETAS), n_th0 + 1L)
  expect_true("WT" %in% as.character(m2$COVARIATES))
  expect_true("WT" %in% m2$INPUT)
  expect_true(grepl("WT", m2$PRED))
  # new THETA initialised at 0 (null effect)
  expect_equal(m2$THETAS$Value[nrow(m2$THETAS)], 0)
})

test_that("forest table is empty for no retained effects", {
  sim <- nm_synthetic_iv1(n_sub = 2L, seed = 1L)
  ft <- .nm_scm_forest_table(sim$model, NULL, sim$model, list())
  expect_s3_class(ft, "data.frame")
  expect_equal(nrow(ft), 0L)
  expect_true(all(c("parameter", "covariate", "estimate", "lower", "upper") %in%
                    names(ft)))
})

test_that("SCM selects a true covariate effect and rejects a null one", {
  skip_on_cran()
  set.seed(4321)
  sim <- nm_synthetic_iv1(n_sub = 16L, seed = 11L)
  m <- sim$model
  d <- data.table::as.data.table(sim$data$data)
  ids <- sort(unique(d$ID))
  wt  <- stats::setNames(round(stats::rnorm(length(ids), 70, 15), 1), ids)
  age <- stats::setNames(round(stats::rnorm(length(ids), 50, 12), 1), ids)
  d$WT  <- wt[as.character(d$ID)]
  d$AGE <- age[as.character(d$ID)]

  # simulate DV from a model where V depends on WT (power, exponent 0.9)
  true_m <- .nm_scm_add_effect(m, "V", "WT", "power", stats::median(wt))
  true_m$THETAS$Value[nrow(true_m$THETAS)] <- 0.9
  simres <- nm_simulate(true_m, nm_dataset_from_table(d), seed = 99L)
  simdt <- data.table::as.data.table(simres[[1L]]$data)
  d$DV <- simdt$DV[match(paste(d$ID, d$TIME), paste(simdt$ID, simdt$TIME))]

  scm <- nm_scm(
    m, nm_dataset_from_table(d),
    candidates = list(
      list(parameter = "V", covariate = "WT",  relationship = "power"),
      list(parameter = "V", covariate = "AGE", relationship = "power")
    ),
    method = "FO", control = list(n_cores = 1L), verbose = FALSE
  )

  expect_s3_class(scm, "nm_scm")
  # exactly the true effect retained
  expect_equal(length(scm$retained), 1L)
  expect_identical(scm$retained[[1]]$covariate, "WT")
  # steps report dOFV / AIC / BIC columns
  expect_true(all(c("dOFV", "AIC", "BIC", "significant") %in% names(scm$steps)))
  # the WT forward step is significant and large
  wt_step <- scm$steps[scm$steps$covariate == "WT" & scm$steps$action == "forward", ]
  expect_true(any(wt_step$significant))
  expect_gt(max(wt_step$dOFV, na.rm = TRUE), 10)
  # forest carries the estimate near the truth with a finite CI
  expect_equal(nrow(scm$forest), 1L)
  expect_true(is.finite(scm$forest$estimate[1]))
  expect_gt(scm$forest$estimate[1], 0.4)
})

test_that("nm_forest_plot errors without ggplot2 or with no effects", {
  sim <- nm_synthetic_iv1(n_sub = 2L, seed = 1L)
  empty <- structure(
    list(forest = .nm_scm_forest_table(sim$model, NULL, sim$model, list())),
    class = "nm_scm"
  )
  expect_error(nm_forest_plot(empty), "retained")
  expect_error(nm_forest_plot(list()), "nm_scm")
})
