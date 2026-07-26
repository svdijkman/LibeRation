test_that("visual diagrams generate editable PK and DES scaffolds", {
  diagram <- nm_model_diagram(
    compartments = data.frame(
      id = 1:2, name = c("CENTRAL", "PERIPHERAL"), kind = "amount",
      volume_parameter = c("VC", "VP"), scale_parameter = c("VC", "VP"),
      dose = c(TRUE, FALSE), observe = c(TRUE, FALSE),
      x = c(120, 360), y = c(180, 180)
    ),
    flows = data.frame(
      id = c("distribution", "elimination"), from = c(1, 1), to = c(2, 0),
      type = c("bidirectional_clearance", "clearance"),
      parameter = c("Q", "CL"), secondary_parameter = "",
      expression = "", label = c("Q", "CL")
    ),
    parameters = data.frame(
      name = c("VC", "VP", "Q", "CL"), initial = c(20, 30, 2, 1),
      lower = NA_real_, upper = NA_real_, fixed = FALSE,
      iiv = TRUE, eta_variance = 0.1
    )
  )
  preview <- nm_diagram_preview(diagram)
  expect_match(preview$PRED, "VC = THETA(1) * exp(ETA(1))", fixed = TRUE)
  expect_match(preview$DES, "Q * (A(1) / VC)", fixed = TRUE)
  expect_match(preview$DES, "Q * (A(2) / VP)", fixed = TRUE)

  model <- nm_diagram_generate(diagram)
  expect_s3_class(model, "nm_model")
  expect_equal(model$ADVAN, 6L)
  expect_equal(model$n_eta, 4L)
  expect_s3_class(nm_compile(model), "NMEngine")
  expect_s3_class(nm_model_diagram_get(model), "nm_model_diagram")

  manual <- model
  manual$DES <- paste(model$DES, "# user edit", sep = "\n")
  expect_equal(nm_model_diagram_get(manual)$generated$DES, preview$DES)
})

test_that("nonlinear and custom flows share the diagram code path", {
  diagram <- nm_model_diagram(
    compartments = data.frame(
      id = 1:2, name = c("CENTRAL", "RESPONSE"),
      kind = c("amount", "response"),
      volume_parameter = c("V", ""), scale_parameter = c("V", ""),
      dose = c(TRUE, FALSE), observe = c(TRUE, FALSE),
      x = c(100, 340), y = c(150, 220)
    ),
    flows = data.frame(
      id = c("mm", "response-in", "response-out"),
      from = c(1, 0, 2), to = c(0, 2, 0),
      type = c("michaelis_menten", "custom", "rate"),
      parameter = c("VMAX", "", "KOUT"),
      secondary_parameter = c("KM", "", ""),
      expression = c("", "KIN * (1 - IMAX * C(1) / (IC50 + C(1)))", ""),
      label = c("Saturable elimination", "Inhibition", "Loss")
    ),
    parameters = data.frame(
      name = c("V", "VMAX", "KM", "KIN", "IMAX", "IC50", "KOUT"),
      initial = c(20, 10, 2, 1, 0.8, 3, 0.1), lower = NA_real_, upper = NA_real_,
      fixed = FALSE, iiv = TRUE, eta_variance = 0.1
    )
  )
  preview <- nm_diagram_preview(diagram)
  expect_match(preview$DES, "VMAX * (A(1) / V) / (KM + (A(1) / V))", fixed = TRUE)
  expect_match(preview$DES, "KIN * (1 - IMAX * (A(1) / V)", fixed = TRUE)
  expect_s3_class(nm_compile(nm_diagram_generate(diagram)), "NMEngine")
})

test_that("visual diagrams can target every general ODE ADVAN", {
  for (advan in c(6L, 8L, 9L, 13L, 14L)) {
    diagram <- nm_model_diagram(
      compartments = data.frame(
        id = 1L, name = "CENTRAL", kind = "amount",
        volume_parameter = "V", scale_parameter = "V",
        dose = TRUE, observe = TRUE, x = 100, y = 100
      ),
      flows = data.frame(
        id = "elimination", from = 1L, to = 0L, type = "clearance",
        parameter = "CL", secondary_parameter = "", expression = "",
        label = "Elimination"
      ),
      parameters = data.frame(
        name = c("V", "CL"), initial = c(20, 2),
        lower = NA_real_, upper = NA_real_, fixed = FALSE,
        iiv = TRUE, eta_variance = .1
      ),
      advan = advan
    )
    model <- nm_diagram_generate(diagram)
    expect_equal(model$ADVAN, advan)
    expect_s3_class(model, "nm_model")
  }
})
