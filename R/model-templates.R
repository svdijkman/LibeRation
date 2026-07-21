.nm_template_residual <- function(type) {
  type <- match.arg(type, c("proportional", "additive", "combined", "lognormal", "none"))
  switch(type,
    proportional = list(code = "Y = F * (1 + ERR(1))",
                        sigma = data.frame(SIGMA = 1L, Value = 0.1)),
    additive = list(code = "Y = F + ERR(1)",
                    sigma = data.frame(SIGMA = 1L, Value = 0.1)),
    combined = list(code = "Y = F * (1 + ERR(1)) + ERR(2)",
                    sigma = data.frame(SIGMA = 1:2, Value = c(0.1, 0.1))),
    lognormal = list(code = "Y = F * exp(ERR(1))",
                     sigma = data.frame(SIGMA = 1L, Value = 0.1)),
    none = list(code = "Y = F", sigma = NULL)
  )
}

.nm_template_parameters <- function(names, values, iiv, bounded = character()) {
  theta <- data.frame(
    THETA = seq_along(names), Value = as.numeric(values),
    LOWER = ifelse(values > 0, values / 1000, -1000),
    UPPER = ifelse(values > 0, values * 1000, 1000), FIX = FALSE,
    stringsAsFactors = FALSE
  )
  code <- vapply(seq_along(names), function(index) {
    eta <- if (isTRUE(iiv)) paste0(" + ETA(", index, ")") else ""
    if (names[[index]] %in% bounded) {
      paste0(names[[index]], " = 1 / (1 + exp(-(THETA(", index, ")", eta, ")))")
    } else {
      paste0(names[[index]], " = THETA(", index, ") * exp(",
             if (isTRUE(iiv)) paste0("ETA(", index, ")") else "0", ")")
    }
  }, character(1))
  omega <- if (isTRUE(iiv)) data.frame(
    OMEGA = seq_along(names), Value = rep(0.1, length(names)), FIX = FALSE
  ) else NULL
  list(theta = theta, omega = omega, code = code)
}

#' Catalogue of advanced structural model templates
#'
#' @return A data frame describing templates accepted by
#'   [nm_model_template()].
#' @export
nm_structural_templates <- function() {
  data.frame(
    template = c(
      "nonlinear_elimination", "transit_absorption", "dual_absorption",
      "parent_metabolite", "effect_compartment", "indirect_response",
      "tumour_growth", "tmdd"
    ),
    model = c(
      "Michaelis-Menten elimination", "Transit-compartment absorption",
      "Parallel first-order absorption", "Parent-metabolite PK",
      "PK with effect compartment", "Indirect-response turnover",
      "PK-tumour growth/inhibition", "Full target-mediated disposition"
    ),
    initial_state = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE),
    notes = c(
      "One-compartment IV model", "Configurable number of transit compartments",
      "Dose fractions use F1/F2", "Observe parent/metabolite using CMT 1/2",
      "Observe plasma/effect-site using CMT 1/2",
      "Initialize response compartment to KIN/KOUT",
      "Initialize tumour compartment with baseline size",
      "Initialize free target compartment to KSYN/KDEG"
    ),
    stringsAsFactors = FALSE
  )
}

#' Create an editable advanced structural model
#'
#' These are ordinary [nm_model()] objects, not a second modelling language.
#' The generated `$PK/$PRED` and `$DES` blocks remain fully editable and run
#' through the same C++/CppAD ODE, likelihood, simulation, and estimation paths.
#'
#' @param template Template identifier from [nm_structural_templates()].
#' @param iiv Add log-normal ETA variability to generated parameters.
#' @param residual Residual model.
#' @param n_transit Number of transit compartments for `transit_absorption`.
#' @param ode_control Optional ADVAN13 solver controls.
#' @return An editable `nm_model` with template notes attached.
#' @export
nm_model_template <- function(
    template = c(
      "nonlinear_elimination", "transit_absorption", "dual_absorption",
      "parent_metabolite", "effect_compartment", "indirect_response",
      "tumour_growth", "tmdd"
    ),
    iiv = TRUE, residual = c("proportional", "additive", "combined", "lognormal", "none"),
    n_transit = 3L, ode_control = NULL) {
  template <- match.arg(template)
  residual <- .nm_template_residual(match.arg(residual))
  input <- c("ID", "TIME", "EVID", "AMT", "RATE", "CMT", "DV", "MDV", "DVID")
  build <- function(names, values, des, obs, dose = 1L, bounded = character(),
                    notes = character(), output = names, scale = NULL) {
    parameters <- .nm_template_parameters(names, values, iiv, bounded)
    scale <- scale %||% if (obs == 1L && "V" %in% names) "V" else "1"
    pred <- paste(c(parameters$code, paste0("S", obs, " = ", scale)), collapse = "\n")
    model <- nm_model(
      INPUT = input, OUTPUT = output, ADVAN = 13L, TRANS = 1L,
      DOSECMP = dose, OBSCMP = obs, PRED = pred, DES = des,
      ERROR = residual$code, THETAS = parameters$theta,
      OMEGAS = parameters$omega, SIGMAS = residual$sigma,
      ODE_CONTROL = ode_control, SOLVER = "ode"
    )
    attr(model, "name") <- nm_structural_templates()$model[
      match(template, nm_structural_templates()$template)
    ]
    attr(model, "template") <- template
    attr(model, "template_notes") <- notes
    model
  }
  if (template == "nonlinear_elimination") {
    return(build(
      c("VMAX", "KM", "V"), c(20, 2, 20),
      "DADT(1) = -VMAX * A(1) / (KM * V + A(1))", 1L,
      notes = "IV Michaelis-Menten elimination; concentration is A(1)/V."
    ))
  }
  if (template == "transit_absorption") {
    n_transit <- as.integer(n_transit)
    if (length(n_transit) != 1L || is.na(n_transit) || n_transit < 1L || n_transit > 20L) {
      .nm_stop("`n_transit` must be between 1 and 20.")
    }
    central <- n_transit + 1L
    derivatives <- c("DADT(1) = -KTR * A(1)")
    if (n_transit > 1L) derivatives <- c(derivatives, vapply(2:n_transit, function(index) {
      paste0("DADT(", index, ") = KTR * (A(", index - 1L, ") - A(", index, "))")
    }, character(1)))
    derivatives <- c(derivatives, paste0(
      "DADT(", central, ") = KTR * A(", n_transit,
      ") - CL / V * A(", central, ")"
    ))
    model <- build(
      c("KTR", "CL", "V"), c(1, 2, 20), paste(derivatives, collapse = "\n"),
      central, notes = paste(n_transit, "serial transit compartments; dose into CMT 1."),
      scale = "V"
    )
    return(model)
  }
  if (template == "dual_absorption") {
    parameters <- .nm_template_parameters(
      c("KA1", "KA2", "CL", "V", "FRAC"), c(1.5, 0.2, 2, 20, 0), iiv,
      bounded = "FRAC"
    )
    pred <- paste(c(parameters$code, "F1 = FRAC", "F2 = 1 - FRAC", "S3 = V"), collapse = "\n")
    model <- nm_model(
      INPUT = input, OUTPUT = c("KA1", "KA2", "CL", "V", "FRAC"),
      ADVAN = 13L, TRANS = 1L, DOSECMP = 1L, OBSCMP = 3L,
      PRED = pred,
      DES = paste("DADT(1) = -KA1 * A(1)", "DADT(2) = -KA2 * A(2)",
                  "DADT(3) = KA1 * A(1) + KA2 * A(2) - CL / V * A(3)", sep = "\n"),
      ERROR = residual$code, THETAS = parameters$theta, OMEGAS = parameters$omega,
      SIGMAS = residual$sigma, ODE_CONTROL = ode_control, SOLVER = "ode"
    )
    attr(model, "name") <- "Parallel first-order absorption"
    attr(model, "template") <- template
    attr(model, "template_notes") <- "Use matching dose records in CMT 1 and CMT 2; F1/F2 split the bioavailable amount."
    return(model)
  }
  if (template == "parent_metabolite") {
    parameters <- .nm_template_parameters(
      c("CLP", "VP", "FM", "CLM", "VM"), c(2, 20, 0, 1.5, 30), iiv,
      bounded = "FM"
    )
    pred <- paste(c(parameters$code, "S1 = VP", "S2 = VM"), collapse = "\n")
    model <- nm_model(
      INPUT = input, OUTPUT = c("CLP", "VP", "FM", "CLM", "VM"),
      ADVAN = 13L, TRANS = 1L, DOSECMP = 1L, OBSCMP = 1L, PRED = pred,
      DES = paste(
        "DADT(1) = -CLP / VP * A(1)",
        "DADT(2) = FM * CLP / VP * A(1) - CLM / VM * A(2)", sep = "\n"
      ), ERROR = residual$code, THETAS = parameters$theta, OMEGAS = parameters$omega,
      SIGMAS = residual$sigma, ODE_CONTROL = ode_control, SOLVER = "ode"
    )
    attr(model, "name") <- "Parent-metabolite PK"
    attr(model, "template") <- template
    attr(model, "template_notes") <- "Use CMT 1 for parent and CMT 2 for metabolite observations."
    return(model)
  }
  if (template == "effect_compartment") {
    return(build(
      c("CL", "V", "KE0"), c(2, 20, 0.5),
      paste("DADT(1) = -CL / V * A(1)",
            "DADT(2) = KE0 * (A(1) / V - A(2))", sep = "\n"),
      2L, notes = "A(2) is effect-site concentration; use CMT 1/2 for plasma/effect observations."
    ))
  }
  if (template == "indirect_response") {
    return(build(
      c("CL", "V", "KIN", "KOUT", "IC50"), c(2, 20, 10, 0.2, 2),
      paste(
        "DADT(1) = -CL / V * A(1)",
        "DADT(2) = KIN * (1 - (A(1) / V) / (IC50 + A(1) / V)) - KOUT * A(2)",
        sep = "\n"
      ), 2L, notes = "Before dosing, initialize CMT 2 with AMT = KIN/KOUT on an EVID=1 record."
    ))
  }
  if (template == "tumour_growth") {
    return(build(
      c("CL", "V", "KG", "KCAP", "KILL"), c(2, 20, 0.03, 100, 0.01),
      paste(
        "DADT(1) = -CL / V * A(1)",
        "DADT(2) = KG * A(2) * (1 - A(2) / KCAP) - KILL * A(1) / V * A(2)",
        sep = "\n"
      ), 2L, notes = "Initialize CMT 2 to baseline tumour size with an EVID=1 record."
    ))
  }
  build(
    c("CL", "V", "KON", "KOFF", "KINT", "KSYN", "KDEG"),
    c(2, 20, 0.1, 0.05, 0.02, 1, 0.1),
    paste(
      "DADT(1) = -CL / V * A(1) - KON / V * A(1) * A(2) + KOFF * A(3)",
      "DADT(2) = KSYN - KDEG * A(2) - KON / V * A(1) * A(2) + KOFF * A(3)",
      "DADT(3) = KON / V * A(1) * A(2) - (KOFF + KINT) * A(3)",
      sep = "\n"
    ), 1L, notes = paste(
      "Full TMDD. Initialize free target in CMT 2 to KSYN/KDEG before drug dosing;",
      "CMT 3 contains drug-target complex."
    )
  )
}

#' Generate a piecewise-constant model expression
#'
#' @param time Symbol containing time, normally `TIME`.
#' @param knots Increasing interval boundaries.
#' @param values One value/expression per interval (`length(knots) + 1`).
#' @return Editable nested `ifelse()` code for `$PK/$PRED` or `$ERROR`.
#' @export
nm_piecewise <- function(time = "TIME", knots, values) {
  time <- trimws(as.character(time))
  knots <- as.numeric(knots)
  values <- as.character(values)
  if (length(time) != 1L || !grepl("^[A-Za-z][A-Za-z0-9_.]*$", time) ||
      any(!is.finite(knots)) || is.unsorted(knots, strictly = TRUE) ||
      length(values) != length(knots) + 1L || any(!nzchar(trimws(values)))) {
    .nm_stop("Piecewise expressions require a time symbol, increasing knots, and one more value than knots.")
  }
  expression <- values[[length(values)]]
  for (index in rev(seq_along(knots))) {
    expression <- paste0("ifelse(", time, " < ", .nm_code_number(knots[[index]]),
                         ", ", values[[index]], ", ", expression, ")")
  }
  expression
}

#' Generate a restricted cubic spline expression
#'
#' @param x Predictor symbol, usually `TIME` or a covariate.
#' @param knots At least three strictly increasing knots.
#' @param coefficients Coefficients for the linear term followed by the
#'   `length(knots) - 2` nonlinear basis terms.
#' @param intercept Optional intercept expression.
#' @return Editable restricted-cubic-spline expression using only tape-safe
#'   arithmetic and `pmax()`.
#' @export
nm_spline <- function(x = "TIME", knots, coefficients, intercept = "0") {
  x <- trimws(as.character(x))
  knots <- as.numeric(knots)
  coefficients <- as.character(coefficients)
  intercept <- as.character(intercept)
  if (length(x) != 1L || !grepl("^[A-Za-z][A-Za-z0-9_.]*$", x) ||
      length(knots) < 3L || any(!is.finite(knots)) ||
      is.unsorted(knots, strictly = TRUE) ||
      length(coefficients) != length(knots) - 1L ||
      length(intercept) != 1L || !nzchar(trimws(intercept))) {
    .nm_stop(
      "Restricted cubic splines require a predictor, at least three increasing knots, ",
      "an intercept, and one linear plus K-2 nonlinear coefficients."
    )
  }
  last <- knots[[length(knots)]]
  penultimate <- knots[[length(knots) - 1L]]
  truncated <- function(knot) paste0("pmax(", x, " - ", .nm_code_number(knot), ", 0)^3")
  basis <- vapply(seq_len(length(knots) - 2L), function(index) {
    knot <- knots[[index]]
    paste0(
      "(", truncated(knot), " - ", truncated(penultimate), " * ",
      .nm_code_number((last - knot) / (last - penultimate)), " + ",
      truncated(last), " * ",
      .nm_code_number((penultimate - knot) / (last - penultimate)), ")"
    )
  }, character(1))
  terms <- c(paste0("(", coefficients[[1L]], ") * ", x),
             paste0("(", coefficients[-1L], ") * ", basis))
  paste(c(paste0("(", intercept, ")"), terms), collapse = " + ")
}
