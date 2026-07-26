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

#' Catalogue of editable model templates
#'
#' @return A data frame describing templates accepted by
#'   [nm_model_template()].
#' @export
nm_structural_templates <- function() {
  data.frame(
    template = c(
      "nonlinear_elimination", "transit_absorption", "dual_absorption",
      "parent_metabolite", "effect_compartment", "indirect_response",
      "tumour_growth", "tmdd", "bernoulli", "categorical", "ordinal",
      "poisson", "negative_binomial", "time_to_event", "recurrent_event",
      "competing_risks", "markov", "continuous_time_markov",
      "hidden_markov", "continuous_time_hidden_markov"
    ),
    model = c(
      "Michaelis-Menten elimination", "Transit-compartment absorption",
      "Parallel first-order absorption", "Parent-metabolite PK",
      "PK with effect compartment", "Indirect-response turnover",
      "PK-tumour growth/inhibition", "Full target-mediated disposition",
      "Binary outcome (Bernoulli)", "Multicategory outcome",
      "Ordered categorical outcome", "Count outcome (Poisson)",
      "Overdispersed count outcome", "Time to first event",
      "Recurrent-event intensity", "Competing-risks hazards",
      "Discrete-time Markov model", "Continuous-time Markov model",
      "Hidden Markov model", "Continuous-time hidden Markov model"
    ),
    initial_state = c(
      FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE,
      rep(FALSE, 12L)
    ),
    notes = c(
      "One-compartment IV model", "Configurable number of transit compartments",
      "Dose fractions use F1/F2", "Observe parent/metabolite using CMT 1/2",
      "Observe plasma/effect-site using CMT 1/2",
      "Initialize response compartment to KIN/KOUT",
      "Initialize tumour compartment with baseline size",
      "Initialize free target compartment to KSYN/KDEG",
      "DV must be coded 0/1", "DV defaults to categories 0, 1, and 2",
      "DV defaults to ordered categories 0, 1, and 2",
      "DV must be a non-negative integer count",
      "DV must be a non-negative integer count",
      "DV=1 denotes the event; rows define risk intervals",
      "DV=1 denotes each recurrent event",
      "DV=0 is no event; DV=1/2 are the competing causes",
      "DV is the observed state coded 0/1",
      "DV is the observed state coded 0/1; irregular times are supported",
      "A two-state binary-emission HMM with smoothing/Viterbi support",
      "A two-state binary-emission CT-HMM with irregular observation times"
    ),
    stringsAsFactors = FALSE
  )
}

.nm_likelihood_model_template <- function(template, iiv = TRUE) {
  input <- c("ID", "TIME", "DV", "MDV", "DVID", "EVID", "AMT", "CMT")
  eta <- if (isTRUE(iiv)) " + ETA(1)" else ""
  omega <- if (isTRUE(iiv)) {
    data.frame(OMEGA = 1L, Value = 0.1, FIX = FALSE)
  } else NULL
  theta <- function(values) data.frame(
    THETA = seq_along(values), Value = as.numeric(values),
    LOWER = rep(-10, length(values)), UPPER = rep(10, length(values)),
    FIX = FALSE
  )
  build <- function(pred, values, outcomes = NULL, error = NULL,
                    hmm = NULL, output = NULL) {
    arguments <- list(
      INPUT = input, OUTPUT = output, ADVAN = 1L, TRANS = 1L,
      DOSECMP = 1L, OBSCMP = 1L, PRED = pred,
      THETAS = theta(values), OMEGAS = omega
    )
    if (!is.null(outcomes)) arguments$OUTCOMES <- outcomes
    if (!is.null(error)) arguments$ERROR <- error
    if (!is.null(hmm)) arguments$HMM_CONFIG <- hmm
    model <- do.call(nm_model, arguments)
    catalogue <- nm_structural_templates()
    row <- match(template, catalogue$template)
    attr(model, "name") <- catalogue$model[[row]]
    attr(model, "template") <- template
    attr(model, "template_notes") <- catalogue$notes[[row]]
    model
  }
  base <- "CL=1\nV=1\nS1=V"
  if (template == "bernoulli") {
    return(build(
      paste0("LP=THETA(1)", eta, "\nP=1/(1+exp(-LP))\nF=P\n", base),
      0, nm_outcome("bernoulli", prediction = "P"), output = c("LP", "P")
    ))
  }
  if (template %in% c("categorical", "ordinal")) {
    pred <- paste0(
      "L1=THETA(1)", eta, "\nL2=THETA(2)\n",
      "E1=exp(L1)\nE2=exp(L2)\nDEN=1+E1+E2\n",
      "P0=1/DEN\nP1=E1/DEN\nP2=E2/DEN\nF=P1\n", base
    )
    return(build(
      pred, c(0, 0),
      nm_outcome(
        template, prediction = "P1",
        probabilities = c("P0", "P1", "P2"), categories = 0:2
      ),
      output = c("P0", "P1", "P2")
    ))
  }
  if (template == "poisson") {
    return(build(
      paste0("LOGMU=THETA(1)", eta, "\nMU=exp(LOGMU)\nF=MU\n", base),
      log(2), nm_outcome("poisson", prediction = "MU", max_count = 100L),
      output = c("LOGMU", "MU")
    ))
  }
  if (template == "negative_binomial") {
    return(build(
      paste0(
        "LOGMU=THETA(1)", eta,
        "\nMU=exp(LOGMU)\nSIZE=exp(THETA(2))\nF=MU\n", base
      ),
      log(c(2, 3)),
      nm_outcome(
        "negative_binomial", prediction = "MU",
        dispersion = "SIZE", max_count = 100L
      ),
      output = c("LOGMU", "MU", "SIZE")
    ))
  }
  if (template %in% c("time_to_event", "recurrent_event")) {
    family <- if (template == "time_to_event") "tte" else "recurrent_event"
    return(build(
      paste0("LOGHAZ=THETA(1)", eta, "\nHAZ=exp(LOGHAZ)\nF=HAZ\n", base),
      log(0.1), nm_outcome(family, prediction = "HAZ"),
      output = c("LOGHAZ", "HAZ")
    ))
  }
  if (template == "competing_risks") {
    return(build(
      paste0(
        "LOGH1=THETA(1)", eta,
        "\nH1=exp(LOGH1)\nH2=exp(THETA(2))\nF=H1+H2\n", base
      ),
      log(c(0.1, 0.05)),
      nm_outcome(
        "competing_risks", prediction = "H1",
        cause_hazards = c(`1` = "H1", `2` = "H2")
      ),
      output = c("H1", "H2")
    ))
  }
  if (template == "markov") {
    pred <- paste0(
      "P0=1/(1+exp(-THETA(1)))\nP1=1-P0\n",
      "T00=1/(1+exp(-(THETA(2)", eta, ")))\nT01=1-T00\n",
      "T10=1/(1+exp(-THETA(3)))\nT11=1-T10\nF=P1\n", base
    )
    return(build(
      pred, stats::qlogis(c(0.6, 0.8, 0.3)),
      nm_outcome(
        "markov", prediction = "P1", categories = c(0, 1),
        initial = c("P0", "P1"),
        transition = matrix(
          c("T00", "T01", "T10", "T11"), 2, 2, byrow = TRUE
        )
      ),
      output = c("P0", "P1", "T00", "T01", "T10", "T11")
    ))
  }
  if (template == "continuous_time_markov") {
    pred <- paste0(
      "PI0=1/(1+exp(-THETA(1)))\nPI1=1-PI0\n",
      "Q01=exp(THETA(2)", eta, ")\nQ10=exp(THETA(3))\nF=PI1\n", base
    )
    return(build(
      pred, c(stats::qlogis(0.7), log(0.2), log(0.4)),
      nm_outcome(
        "continuous_time_markov", prediction = "PI1",
        categories = c(0, 1), initial = c("PI0", "PI1"),
        rates = c("Q01", "Q10")
      ),
      output = c("PI0", "PI1", "Q01", "Q10")
    ))
  }
  if (template == "hidden_markov") {
    error <- paste0(
      "I1=1/(1+exp(-THETA(1)))\nI2=1-I1\n",
      "T11=1/(1+exp(-(THETA(2)", eta, ")))\nT12=1-T11\n",
      "T21=1/(1+exp(-THETA(3)))\nT22=1-T21\n",
      "E10=1/(1+exp(-THETA(4)))\nE20=1/(1+exp(-THETA(5)))\n",
      "E1=ifelse(DV==0,E10,1-E10)\nE2=ifelse(DV==0,E20,1-E20)"
    )
    return(build(
      paste("F=0", base, sep = "\n"),
      stats::qlogis(c(0.6, 0.8, 0.3, 0.9, 0.2)),
      error = error,
      hmm = nm_hmm_config(
        states = c("low", "high"), initial = c("I1", "I2"),
        transition = matrix(
          c("T11", "T12", "T21", "T22"), 2, 2, byrow = TRUE
        ),
        emission = c("E1", "E2"), by_dvid = TRUE
      )
    ))
  }
  if (template == "continuous_time_hidden_markov") {
    error <- paste0(
      "I1=1/(1+exp(-THETA(1)))\nI2=1-I1\n",
      "Q12=exp(THETA(2)", eta, ")\nQ21=exp(THETA(3))\n",
      "E10=1/(1+exp(-THETA(4)))\nE20=1/(1+exp(-THETA(5)))\n",
      "E1=ifelse(DV==0,E10,1-E10)\nE2=ifelse(DV==0,E20,1-E20)"
    )
    return(build(
      paste("F=0", base, sep = "\n"),
      c(stats::qlogis(0.6), log(0.2), log(0.4),
        stats::qlogis(0.9), stats::qlogis(0.2)),
      error = error,
      hmm = nm_cthmm_config(
        states = c("low", "high"), initial = c("I1", "I2"),
        generator = matrix(c("", "Q12", "Q21", ""), 2, 2, byrow = TRUE),
        emission = c("E1", "E2"), by_dvid = TRUE
      )
    ))
  }
  .nm_stop("Unknown likelihood model template: ", template, ".")
}

#' Create an editable model from a first-class template
#'
#' These are ordinary [nm_model()] objects, not a second modelling language.
#' Structural templates generate editable `$PK/$PRED` and `$DES` blocks.
#' Outcome, event, Markov, and hidden-Markov templates add their corresponding
#' likelihood configuration. Every template runs through the same C++/CppAD
#' compilation, likelihood, simulation, and estimation paths.
#'
#' @param template Template identifier from [nm_structural_templates()].
#' @param iiv Add ETA variability to generated parameters or linear predictors.
#' @param residual Residual model for structural PK/PD templates.
#' @param n_transit Number of transit compartments for `transit_absorption`.
#' @param ode_control Optional ADVAN13 solver controls.
#' @return An editable `nm_model` with template notes attached.
#' @export
nm_model_template <- function(
    template = "nonlinear_elimination",
    iiv = TRUE, residual = c("proportional", "additive", "combined", "lognormal", "none"),
    n_transit = 3L, ode_control = NULL) {
  template <- match.arg(template, nm_structural_templates()$template)
  likelihood_templates <- c(
    "bernoulli", "categorical", "ordinal", "poisson",
    "negative_binomial", "time_to_event", "recurrent_event",
    "competing_risks", "markov", "continuous_time_markov",
    "hidden_markov", "continuous_time_hidden_markov"
  )
  if (template %in% likelihood_templates) {
    return(.nm_likelihood_model_template(template, iiv = iiv))
  }
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
