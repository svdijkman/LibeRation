test_that("NONMEM control streams round-trip supported semantic records", {
  control <- paste(
    "$PROBLEM Round-trip model",
    "$INPUT ID TIME EVID AMT RATE II SS CMT DV MDV WT",
    "$DATA input.csv IGNORE=@",
    "$SUBROUTINES ADVAN2 TRANS2",
    "$PK",
    "KA=THETA(1)",
    "CL=THETA(2)*exp(ETA(1))",
    "V=THETA(3)",
    "S2=V",
    "$ERROR",
    "Y=F+F*EPS(1)",
    "$THETA",
    "(0,1,5)", "(0,5,20)", "50 FIX",
    "$OMEGA BLOCK(1)", "0.1",
    "$SIGMA", "0.05",
    "$ESTIMATION METHOD=1 INTERACTION",
    "$COVARIANCE MATRIX=S",
    "$TABLE ID TIME DV PRED",
    sep = "\n"
  )
  imported <- nm_control_read(control)
  expect_s3_class(imported, "nm_control_stream")
  expect_s3_class(imported$model, "nm_model")
  expect_equal(imported$model$ADVAN, 2L)
  expect_equal(imported$model$TRANS, 2L)
  expect_equal(imported$model$THETAS$Value, c(1, 5, 50))
  expect_equal(imported$model$THETAS$FIX, c(FALSE, FALSE, TRUE))
  expect_equal(imported$model$OMEGAS$Value, 0.1)
  expect_match(imported$model$ERROR, "ERR(1)", fixed = TRUE)
  expect_equal(imported$model$OUTPUT, "PRED")

  round_trip <- nm_control_read(nm_control_write(imported))
  expect_equal(round_trip$model$ADVAN, imported$model$ADVAN)
  expect_equal(round_trip$model$TRANS, imported$model$TRANS)
  expect_equal(round_trip$model$THETAS$Value, imported$model$THETAS$Value)
  expect_equal(round_trip$model$OMEGAS[, c("ROW", "COL", "Value")],
               imported$model$OMEGAS[, c("ROW", "COL", "Value")])
  expect_equal(round_trip$model$SIGMAS$Value, imported$model$SIGMAS$Value)
})

test_that("selected generated columns create a NONMEM TABLE record", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "DV"),
    OUTPUT = c("PRED", "CL", "CWRES"),
    ADVAN = 1,
    PRED = "CL=THETA(1); V=THETA(2); S1=V",
    ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(THETA = 1:2, Value = c(2, 20)),
    SIGMAS = data.frame(SIGMA = 1, Value = 0.1)
  )
  text <- nm_control_write(model)
  expect_match(text, "$TABLE ID TIME DV PRED CL CWRES", fixed = TRUE)
  expect_equal(nm_control_read(text)$model$OUTPUT, c("PRED", "CL", "CWRES"))
})

test_that("unknown NONMEM records are preserved and reported", {
  control <- paste(
    "$PROBLEM Preserved",
    "$INPUT ID TIME EVID AMT DV",
    "$SUBROUTINES ADVAN1 TRANS2",
    "$PK", "CL=THETA(1); V=THETA(2); S1=V",
    "$ERROR", "Y=F",
    "$THETA", "1", "10",
    "$MSFI old.msf",
    sep = "\n"
  )
  imported <- nm_control_read(control)
  expect_true(imported$compatibility$requires_manual_translation)
  expect_true("MSFI" %in% imported$compatibility$preserved_records)
  expect_match(nm_control_write(imported), "$MSFI", fixed = TRUE)
})

test_that("NONMEM likelihood records round-trip through the compiled path", {
  control <- paste(
    "$PROBLEM Binary likelihood",
    "$INPUT ID TIME DV MDV",
    "$DATA binary.csv",
    "$SUBROUTINES ADVAN1 TRANS2",
    "$PRED",
    "P=1/(1+EXP(-THETA(1)))",
    "CL=1; V=1; S1=V; F=P",
    "$ERROR",
    "Y=ifelse(DV.EQ.1,F,1-F)",
    "$THETA (-10,0,10)",
    "$ESTIMATION METHOD=COND LAPLACE LIKELIHOOD",
    sep = "\n"
  )
  imported <- nm_control_read(control)
  expect_identical(imported$model$LIK_CONFIG$error, "likelihood")
  expect_identical(imported$model$likelihood_output, "Y")
  written <- nm_control_write(imported$model)
  expect_match(written, "LIKELIHOOD")
  expect_identical(nm_control_read(written)$model$LIK_CONFIG$error, "likelihood")
})

test_that("direct PRED and combined PK-PRED modes round-trip explicitly", {
  direct <- nm_model(
    INPUT = c("ID", "TIME", "DV", "MDV"), PRED_MODE = "pred",
    PRED_SOURCE = "F=THETA(1)*TIME", ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(THETA = 1, Value = 2),
    SIGMAS = data.frame(SIGMA = 1, Value = 0.1)
  )
  direct_text <- nm_control_write(direct)
  expect_false(grepl("$SUBROUTINES", direct_text, fixed = TRUE))
  direct_rebuilt <- nm_control_read(direct_text)$model
  expect_identical(direct_rebuilt$PRED_MODE, "pred")
  expect_identical(direct_rebuilt$PRED_SOURCE, direct$PRED_SOURCE)

  combined <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "DV", "MDV"), ADVAN = 1,
    PRED_MODE = "pk_pred",
    PK_SOURCE = "CL=THETA(1);V=THETA(2);S1=V",
    PRED_SOURCE = "RAW=F_ADVAN;F=RAW+A(1)/100",
    ERROR = "Y=F",
    THETAS = data.frame(THETA = 1:2, Value = c(2, 20))
  )
  combined_text <- nm_control_write(combined)
  expect_match(combined_text, "LIBERATION_POST_PRED_BEGIN", fixed = TRUE)
  expect_false(grepl("\n$PRED", combined_text, fixed = TRUE))
  combined_rebuilt <- nm_control_read(combined_text)$model
  expect_identical(combined_rebuilt$PRED_MODE, "pk_pred")
  expect_identical(combined_rebuilt$PK_SOURCE, combined$PK_SOURCE)
  expect_identical(combined_rebuilt$PRED_SOURCE, combined$PRED_SOURCE)
})

test_that("a conventional monolithic NONMEM PRED is split at residual code", {
  control <- paste(
    "$PROBLEM Monolithic direct model",
    "$INPUT ID TIME DV MDV",
    "$PRED",
    "F=THETA(1)*TIME",
    "Y=F+EPS(1)",
    "$THETA 2",
    "$SIGMA 1",
    sep = "\n"
  )
  imported <- nm_control_read(control)$model
  expect_identical(imported$PRED_MODE, "pred")
  expect_match(imported$PRED_SOURCE, "F=THETA", fixed = TRUE)
  expect_false(grepl("Y=", imported$PRED_SOURCE, fixed = TRUE))
  expect_match(imported$ERROR, "Y=F+ERR(1)", fixed = TRUE)
  expect_s3_class(nm_compile(imported), "NMEngine")
})

test_that("CWRES are generated and remain aligned with fitted records", {
  fixture <- estimation_fixture()
  fit <- nm_est(fixture$model, fixture$data, method = "FOCEI", maxit = 2, eta_maxit = 5)
  gof <- nm_gof(fit)
  expect_true("CWRES" %in% names(gof))
  observed <- gof$EVID == 0L & gof$MDV == 0L
  expect_true(all(is.finite(gof$CWRES[observed])))
  expect_equal(residuals(fit), gof$IWRES)
  expect_equal(residuals(fit, "CWRES"), gof$CWRES)
})

test_that("binary categorical and hazard VPCs produce saved summaries", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "DV", "MDV"), ADVAN = 1,
    PRED = "CL=THETA(1); V=THETA(2); S1=V", ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(THETA = 1:2, Value = c(0.2, 2), FIX = TRUE),
    SIGMAS = data.frame(SIGMA = 1, Value = 0.01, FIX = TRUE)
  )
  data <- do.call(rbind, lapply(1:5, function(id) data.frame(
    ID = id, TIME = c(0, 1, 2, 4), EVID = c(1, 0, 0, 0),
    AMT = c(1, 0, 0, 0), MDV = c(1, 0, 0, 0),
    DV = c(NA, 0, as.integer(id %% 2 == 0), as.integer(id %% 3 == 0))
  )))
  fit <- nm_est(model, data, method = "FO", maxit = 1)
  categorical <- nm_vpc_categorical(fit, nsim = 20, seed = 12)
  tte <- nm_vpc_tte(fit, nsim = 20, seed = 12)
  expect_s3_class(categorical, "nm_vpc_categorical")
  expect_true(nrow(categorical$observed) > 0L)
  expect_true(all(c("lower", "median", "upper") %in% names(categorical$simulated)))
  expect_s3_class(tte, "nm_vpc_tte")
  expect_equal(tte$observed$SURVIVAL[[1L]], 1)
  expect_true(all(tte$simulated$lower <= tte$simulated$upper))
})

test_that("profile likelihood refits fixed grids and SCM evaluates candidates", {
  fixture <- estimation_fixture()
  fixture$data$WT <- rep(c(55, 70, 90), each = 4)
  model <- LibeRation:::.nm_model_rebuild(
    fixture$model,
    list(INPUT = c(fixture$model$INPUT, "WT"), COVARIATES = "WT")
  )
  fit <- nm_est(model, fixture$data, method = "FO", maxit = 1)
  profile <- nm_profile(fit, parameters = "THETA1", points = 3, span = 0.2, maxit = 1)
  expect_s3_class(profile, "nm_profile")
  expect_equal(nrow(profile$grid), 3L)
  expect_true(all(is.finite(profile$grid$objective)))

  scm <- nm_scm(
    fit,
    data.frame(parameter = "CL", covariate = "WT", form = "power"),
    direction = "forward", max_steps = 1, maxit = 1
  )
  expect_s3_class(scm, "nm_scm")
  expect_true(nrow(scm$steps) >= 1L)
  expect_s3_class(scm$final_fit, "nm_fit")
})
