test_that("advanced structural templates remain ordinary editable models", {
  catalogue <- nm_structural_templates()
  expect_true(all(c("template", "model", "initial_state", "notes") %in% names(catalogue)))
  expect_gte(nrow(catalogue), 8L)
  for (template in catalogue$template) {
    model <- nm_model_template(template, iiv = FALSE)
    expect_s3_class(model, "nm_model")
    expect_identical(attr(model, "template"), template)
    expect_true(nzchar(model$DES))
    expect_s3_class(nm_compile(model), "NMEngine")
  }
})

test_that("piecewise expressions are editable and tape-safe", {
  code <- nm_piecewise("TIME", c(2, 5), c("THETA(1)", "THETA(2)", "THETA(3)"))
  expect_match(code, "ifelse\\(TIME < 2")
  expect_match(code, "ifelse\\(TIME < 5")
  model <- nm_model(
    INPUT = c("ID", "TIME", "DV", "MDV"), ADVAN = 1,
    PRED = paste0("BASE=", code, "; CL=1; V=1; S1=V; F=BASE"),
    THETAS = data.frame(THETA = 1:3, Value = 1:3)
  )
  expect_s3_class(nm_compile(model), "NMEngine")

  spline <- nm_spline(
    "TIME", knots = c(0, 2, 5, 8),
    coefficients = c("THETA(1)", "THETA(2)", "THETA(3)"),
    intercept = "THETA(4)"
  )
  expect_match(spline, "pmax")
  spline_model <- nm_model(
    INPUT = c("ID", "TIME", "DV", "MDV"), ADVAN = 1,
    PRED = paste0("BASE=", spline, "; CL=1; V=1; S1=V; F=BASE"),
    THETAS = data.frame(THETA = 1:4, Value = c(0.1, 0.01, -0.01, 1))
  )
  expect_s3_class(nm_compile(spline_model), "NMEngine")
})
