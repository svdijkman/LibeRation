#' Declare a first-class pharmacometric outcome
#'
#' `nm_outcome()` describes an observation distribution independently of the
#' structural PK/PD model. LibeRation translates the declaration into editable
#' `$ERROR` code, compiles it with CppAD, and retains enough metadata for
#' stochastic simulation and outcome-appropriate diagnostics. Several
#' declarations can be combined with [nm_outcomes()] and selected by `DVID`.
#'
#' @param family Observation family. Supported families are Gaussian,
#'   log-normal, fixed-degrees-of-freedom Student t, Bernoulli, categorical,
#'   ordinal, Poisson, negative binomial, binomial, zero-inflated Poisson,
#'   hurdle Poisson, first/recurrent time-to-event, competing risks, and
#'   observed finite-state Markov models, and exact two-state continuous-time
#'   Markov models.
#' @param name Human-readable endpoint name.
#' @param dvid Optional numeric `DVID` selecting this endpoint. It is required
#'   when several outcomes are combined.
#' @param prediction Name of a `$PK/$PRED` assignment containing the mean,
#'   probability, intensity, or hazard. `F` uses the engine prediction.
#' @param scale Residual SD for continuous families, as a number, a named
#'   assignment, or `THETA(i)`/`SIGMA(i)`.
#' @param probabilities Named `$PK/$PRED` probability assignments for
#'   categorical or ordinal outcomes.
#' @param categories Numeric observed category codes corresponding to
#'   `probabilities`.
#' @param trials Number of binomial trials, either fixed or supplied by a
#'   numeric dataset/model symbol.
#' @param dispersion Negative-binomial size parameter.
#' @param zero_probability Structural-zero probability for ZIP/hurdle models.
#' @param event Numeric event code for TTE models.
#' @param cause_hazards Named hazard assignments for competing risks; names
#'   are the numeric cause codes.
#' @param initial Initial-state probability assignments for Markov models.
#' @param transition Square matrix of transition-probability assignments for
#'   Markov models; rows are previous states and columns are current states.
#' @param rates Two named off-diagonal intensity assignments (`q01`, `q10`)
#'   for a two-state continuous-time Markov model.
#' @param df Fixed Student-t degrees of freedom.
#' @param max_count Maximum supported count used to unroll exact factorial and
#'   negative-binomial normalizing terms on the AD tape.
#' @param risk Optional at-risk indicator/expression for event models.
#' @return An `nm_outcome` declaration.
#' @export
nm_outcome <- function(
    family = c(
      "normal", "lognormal", "student_t", "bernoulli", "categorical",
      "ordinal", "poisson", "negative_binomial", "binomial",
      "zero_inflated_poisson", "hurdle_poisson", "tte",
      "recurrent_event", "competing_risks", "markov", "continuous_time_markov"
    ),
    name = NULL, dvid = NULL, prediction = "F", scale = NULL,
    probabilities = NULL, categories = NULL, trials = NULL,
    dispersion = NULL, zero_probability = NULL, event = 1,
    cause_hazards = NULL, initial = NULL, transition = NULL, rates = NULL,
    df = 4, max_count = 100L, risk = NULL) {
  family <- match.arg(family)
  scalar_code <- function(value, label, required = TRUE) {
    if (is.null(value)) {
      if (required) .nm_stop("`", label, "` is required for ", family, " outcomes.")
      return(NULL)
    }
    if (is.numeric(value) && length(value) == 1L && is.finite(value)) {
      return(format(value, scientific = TRUE, digits = 17, trim = TRUE))
    }
    value <- trimws(as.character(value))
    if (length(value) != 1L || is.na(value) || !nzchar(value) ||
        !grepl("^(?:[A-Za-z][A-Za-z0-9_.]*|(?:THETA|SIGMA)\\([0-9]+\\))$", value)) {
      .nm_stop(
        "`", label, "` must be a finite number, a simple model/data symbol, ",
        "or THETA(i)/SIGMA(i)."
      )
    }
    value
  }
  prediction <- scalar_code(prediction, "prediction")
  if (is.null(name)) name <- gsub("_", " ", family, fixed = TRUE)
  name <- trimws(as.character(name))
  if (length(name) != 1L || is.na(name) || !nzchar(name)) {
    .nm_stop("`name` must be one non-empty string.")
  }
  if (!is.null(dvid)) {
    dvid <- as.numeric(dvid)
    if (length(dvid) != 1L || !is.finite(dvid)) .nm_stop("`dvid` must be one finite number.")
  }
  max_count <- as.integer(max_count)
  if (length(max_count) != 1L || is.na(max_count) || max_count < 1L || max_count > 10000L) {
    .nm_stop("`max_count` must be an integer between 1 and 10000.")
  }
  df <- as.numeric(df)
  if (family == "student_t" && (length(df) != 1L || !is.finite(df) || df <= 0)) {
    .nm_stop("Student-t `df` must be a positive fixed number.")
  }
  continuous <- family %in% c("normal", "lognormal", "student_t")
  if (continuous) scale <- scalar_code(scale %||% "SIGMA(1)", "scale")
  if (family %in% c("categorical", "ordinal")) {
    probabilities <- trimws(as.character(probabilities))
    if (length(probabilities) < 2L || anyNA(probabilities) || any(!nzchar(probabilities)) ||
        any(!grepl("^[A-Za-z][A-Za-z0-9_.]*$", probabilities))) {
      .nm_stop("Categorical/ordinal outcomes require at least two named probability assignments.")
    }
    categories <- categories %||% seq.int(0, length(probabilities) - 1L)
    categories <- as.numeric(categories)
    if (length(categories) != length(probabilities) || any(!is.finite(categories)) ||
        anyDuplicated(categories)) {
      .nm_stop("`categories` must contain one unique finite code per probability.")
    }
  }
  if (family == "negative_binomial") {
    dispersion <- scalar_code(dispersion, "dispersion")
  }
  if (family == "binomial") trials <- scalar_code(trials, "trials")
  if (family %in% c("zero_inflated_poisson", "hurdle_poisson")) {
    zero_probability <- scalar_code(zero_probability, "zero_probability")
  }
  event <- as.numeric(event)
  if (family %in% c("tte", "recurrent_event") &&
      (length(event) != 1L || !is.finite(event))) {
    .nm_stop("`event` must be one finite numeric event code.")
  }
  risk <- scalar_code(risk, "risk", required = FALSE)
  if (family == "competing_risks") {
    cause_names <- names(cause_hazards)
    cause_hazards <- trimws(as.character(cause_hazards))
    names(cause_hazards) <- cause_names
    causes <- suppressWarnings(as.numeric(cause_names))
    if (length(cause_hazards) < 2L || is.null(cause_names) ||
        anyNA(causes) || any(causes == 0) || anyDuplicated(causes) ||
        any(!grepl("^[A-Za-z][A-Za-z0-9_.]*$", cause_hazards))) {
      .nm_stop(
        "`cause_hazards` must be a named vector of at least two model assignments; ",
        "names are distinct non-zero numeric cause codes."
      )
    }
    categories <- causes
  }
  if (family == "markov") {
    initial <- trimws(as.character(initial))
    if (is.null(transition)) .nm_stop("Markov outcomes require a transition matrix.")
    transition <- as.matrix(transition)
    transition <- matrix(trimws(as.character(transition)), nrow = nrow(transition),
                         ncol = ncol(transition), dimnames = dimnames(transition))
    if (length(initial) < 2L || nrow(transition) != length(initial) ||
        ncol(transition) != length(initial) ||
        any(!grepl("^[A-Za-z][A-Za-z0-9_.]*$", c(initial, transition)))) {
      .nm_stop("Markov outcomes require named initial probabilities and a matching square transition matrix.")
    }
    categories <- categories %||% seq.int(0, length(initial) - 1L)
    categories <- as.numeric(categories)
    if (length(categories) != length(initial) || any(!is.finite(categories)) ||
        anyDuplicated(categories)) {
      .nm_stop("Markov `categories` must contain one unique finite state code per state.")
    }
  }
  if (family == "continuous_time_markov") {
    initial <- trimws(as.character(initial))
    if (length(initial) < 2L ||
        any(!grepl("^[A-Za-z][A-Za-z0-9_.]*$", initial))) {
      .nm_stop(
        "Continuous-time Markov outcomes require at least two named initial-probability assignments."
      )
    }
    if (is.null(dim(rates)) && length(rates) == 2L && length(initial) == 2L) {
      rates <- matrix(c("", rates[[1L]], rates[[2L]], ""), 2, 2, byrow = TRUE)
    } else {
      rates <- as.matrix(rates)
    }
    if (!identical(dim(rates), rep(length(initial), 2L))) {
      .nm_stop("Continuous-time Markov `rates` must be a square generator matrix matching `initial`.")
    }
    rates <- matrix(trimws(as.character(rates)), nrow = length(initial),
                    ncol = length(initial), dimnames = dimnames(rates))
    diagonal <- row(rates) == col(rates)
    rates[diagonal] <- ""
    if (anyNA(rates[!diagonal]) ||
        any(!grepl("^[A-Za-z][A-Za-z0-9_.]*$", rates[!diagonal]))) {
      .nm_stop("Every off-diagonal continuous-time Markov rate must name a model assignment.")
    }
    categories <- as.numeric(categories %||% seq.int(0, length(initial) - 1L))
    if (length(categories) != length(initial) || any(!is.finite(categories)) ||
        anyDuplicated(categories)) {
      .nm_stop("Continuous-time Markov `categories` must contain one unique finite code per state.")
    }
  }
  structure(list(
    version = 1L, name = name, family = family, dvid = dvid,
    prediction = prediction, scale = scale, probabilities = probabilities,
    categories = categories, trials = trials, dispersion = dispersion,
    zero_probability = zero_probability, event = event,
    cause_hazards = cause_hazards, initial = initial,
    transition = transition, rates = rates, df = df, max_count = max_count, risk = risk
  ), class = "nm_outcome")
}

#' Combine endpoint declarations
#'
#' @param ... [nm_outcome()] objects or lists containing them.
#' @return An `nm_outcomes` list.
#' @export
nm_outcomes <- function(...) {
  values <- list(...)
  if (length(values) == 1L && is.list(values[[1L]]) &&
      !inherits(values[[1L]], "nm_outcome")) values <- values[[1L]]
  if (!length(values) || any(!vapply(values, inherits, logical(1), "nm_outcome"))) {
    .nm_stop("`nm_outcomes()` requires one or more `nm_outcome()` declarations.")
  }
  if (length(values) > 1L) {
    dvid <- vapply(values, function(value) value$dvid %||% NA_real_, numeric(1))
    if (anyNA(dvid) || anyDuplicated(dvid)) {
      .nm_stop("Every endpoint in a joint model requires a unique `dvid`.")
    }
  }
  names(values) <- make.unique(vapply(values, `[[`, character(1), "name"))
  structure(values, class = c("nm_outcomes", "list"))
}

#' Declare an item-response/ordered-categorical endpoint set
#'
#' This convenience constructor creates one ordinal [nm_outcome()] per item and
#' combines them by `DVID`. Item characteristic curves remain ordinary named
#' `$PK/$PRED` probability assignments, so graded-response, partial-credit, or
#' custom IRT parameterizations can all use the same compiled engine.
#'
#' @param item_dvid Numeric DVID code for each item.
#' @param probabilities List containing one character probability vector per
#'   item.
#' @param categories Shared category codes or a list of codes per item.
#' @param names Optional item labels.
#' @return An `nm_outcomes` declaration.
#' @export
nm_irt_outcomes <- function(item_dvid, probabilities, categories = NULL,
                            names = NULL) {
  item_dvid <- as.numeric(item_dvid)
  probabilities <- as.list(probabilities)
  if (!length(item_dvid) || length(probabilities) != length(item_dvid) ||
      any(!is.finite(item_dvid)) || anyDuplicated(item_dvid)) {
    .nm_stop("IRT items require one unique finite DVID and one probability vector per item.")
  }
  names <- as.character(names %||% paste("Item", item_dvid))
  if (length(names) != length(item_dvid) || any(!nzchar(trimws(names)))) {
    .nm_stop("`names` must contain one non-empty label per IRT item.")
  }
  category_list <- if (is.list(categories)) categories else
    rep(list(categories), length(item_dvid))
  do.call(nm_outcomes, lapply(seq_along(item_dvid), function(index) {
    nm_outcome(
      "ordinal", name = names[[index]], dvid = item_dvid[[index]],
      prediction = probabilities[[index]][[1L]],
      probabilities = probabilities[[index]], categories = category_list[[index]]
    )
  }))
}

.nm_outcomes <- function(value) {
  if (is.null(value)) return(NULL)
  if (inherits(value, "nm_outcome")) value <- nm_outcomes(value)
  if (!inherits(value, "nm_outcomes")) value <- do.call(nm_outcomes, value)
  value
}

.nm_code_number <- function(value) {
  format(as.numeric(value), scientific = TRUE, digits = 17, trim = TRUE)
}

.nm_outcome_logfactorial <- function(symbol, maximum, prefix) {
  terms <- if (maximum >= 2L) vapply(seq.int(2L, maximum), function(index) {
    paste0("ifelse(", symbol, " >= ", index, ", ", .nm_code_number(log(index)), ", 0)")
  }, character(1)) else character()
  expression <- if (length(terms)) paste(terms, collapse = " + ") else "0"
  c(paste0(prefix, "LOGFACT = ", expression), paste0(prefix, "LOGFACT"))
}

.nm_outcome_probability_choice <- function(value, categories, probabilities) {
  expression <- "1e-300"
  for (index in rev(seq_along(categories))) {
    expression <- paste0(
      "ifelse(", value, " == ", .nm_code_number(categories[[index]]),
      ", ", probabilities[[index]], ", ", expression, ")"
    )
  }
  expression
}

.nm_outcome_component_code <- function(outcome, index) {
  prefix <- paste0("NM_F", index, "_")
  y <- "DV"
  mu <- outcome$prediction
  safe_mu <- paste0("pmax(", mu, ", 1e-12)")
  family <- outcome$family
  lines <- character()
  expression <- switch(
    family,
    normal = {
      variance <- paste0("pmax(", outcome$scale, " * ", outcome$scale, ", 1e-24)")
      paste0("-0.5 * (log(2 * pi * ", variance, ") + (", y, " - ", mu,
             ") * (", y, " - ", mu, ") / ", variance, ")")
    },
    lognormal = {
      variance <- paste0("pmax(", outcome$scale, " * ", outcome$scale, ", 1e-24)")
      paste0("-0.5 * (log(2 * pi * ", variance, ") + (log(pmax(", y,
             ", 1e-300)) - log(", safe_mu, ")) * (log(pmax(", y,
             ", 1e-300)) - log(", safe_mu, ")) / ", variance,
             ") - log(pmax(", y, ", 1e-300))")
    },
    student_t = {
      constant <- lgamma((outcome$df + 1) / 2) - lgamma(outcome$df / 2) -
        0.5 * log(outcome$df * pi)
      paste0(.nm_code_number(constant), " - log(pmax(", outcome$scale,
             ", 1e-12)) - ", .nm_code_number((outcome$df + 1) / 2),
             " * log1p((", y, " - ", mu, ") * (", y, " - ", mu,
             ") / (", .nm_code_number(outcome$df), " * pmax(", outcome$scale,
             " * ", outcome$scale, ", 1e-24)))")
    },
    bernoulli = paste0(
      "log(pmax(ifelse(", y, " == 1, pmin(pmax(", mu,
      ", 0), 1), 1 - pmin(pmax(", mu, ", 0), 1)), 1e-300))"
    ),
    categorical = paste0(
      "log(pmax(", .nm_outcome_probability_choice(y, outcome$categories,
                                                   outcome$probabilities), ", 1e-300))"
    ),
    ordinal = paste0(
      "log(pmax(", .nm_outcome_probability_choice(y, outcome$categories,
                                                   outcome$probabilities), ", 1e-300))"
    ),
    poisson = {
      factorial <- .nm_outcome_logfactorial(y, outcome$max_count, prefix)
      lines <- c(lines, factorial[[1L]])
      paste0(y, " * log(", safe_mu, ") - ", safe_mu, " - ", factorial[[2L]])
    },
    negative_binomial = {
      factorial <- .nm_outcome_logfactorial(y, outcome$max_count, prefix)
      rising_terms <- vapply(seq_len(outcome$max_count), function(k) {
        paste0("ifelse(", y, " >= ", k, ", log(pmax(", outcome$dispersion,
               " + ", k - 1L, ", 1e-12)), 0)")
      }, character(1))
      lines <- c(lines, factorial[[1L]],
                 paste0(prefix, "RISING = ", paste(rising_terms, collapse = " + ")))
      paste0(prefix, "RISING - ", factorial[[2L]], " + ", outcome$dispersion,
             " * log(pmax(", outcome$dispersion, " / (", outcome$dispersion,
             " + ", safe_mu, "), 1e-300)) + ", y,
             " * log(pmax(", safe_mu, " / (", outcome$dispersion, " + ",
             safe_mu, "), 1e-300))")
    },
    binomial = {
      terms <- vapply(seq_len(outcome$max_count), function(k) {
        paste0("ifelse(", y, " >= ", k, ", log(pmax((", outcome$trials,
               " - ", k - 1L, ") / ", k, ", 1e-300)), 0)")
      }, character(1))
      lines <- c(lines, paste0(prefix, "LOGCHOOSE = ", paste(terms, collapse = " + ")))
      probability <- paste0("pmin(pmax(", mu, ", 1e-12), 1 - 1e-12)")
      paste0(prefix, "LOGCHOOSE + ", y, " * log(", probability, ") + (",
             outcome$trials, " - ", y, ") * log(1 - ", probability, ")")
    },
    zero_inflated_poisson = {
      factorial <- .nm_outcome_logfactorial(y, outcome$max_count, prefix)
      lines <- c(lines, factorial[[1L]])
      zero <- paste0("pmin(pmax(", outcome$zero_probability, ", 0), 1)")
      paste0("ifelse(", y, " == 0, log(pmax(", zero, " + (1 - ", zero,
             ") * exp(-", safe_mu, "), 1e-300)), log(pmax(1 - ", zero,
             ", 1e-300)) + ", y, " * log(", safe_mu, ") - ", safe_mu,
             " - ", factorial[[2L]], ")")
    },
    hurdle_poisson = {
      factorial <- .nm_outcome_logfactorial(y, outcome$max_count, prefix)
      lines <- c(lines, factorial[[1L]])
      zero <- paste0("pmin(pmax(", outcome$zero_probability, ", 0), 1)")
      paste0("ifelse(", y, " == 0, log(pmax(", zero,
             ", 1e-300)), log(pmax(1 - ", zero, ", 1e-300)) + ", y,
             " * log(", safe_mu, ") - ", safe_mu, " - ", factorial[[2L]],
             " - log(pmax(1 - exp(-", safe_mu, "), 1e-300)))")
    },
    tte = paste0("ifelse(", y, " == ", .nm_code_number(outcome$event),
                 ", log(", safe_mu, "), 0) - ", safe_mu, " * DT"),
    recurrent_event = paste0("ifelse(", y, " == ", .nm_code_number(outcome$event),
                             ", log(", safe_mu, "), 0) - ", safe_mu, " * DT"),
    competing_risks = {
      total <- paste(outcome$cause_hazards, collapse = " + ")
      selected <- .nm_outcome_probability_choice(y, outcome$categories,
                                                  outcome$cause_hazards)
      paste0("ifelse(", y, " == 0, 0, log(pmax(", selected,
             ", 1e-300))) - (", total, ") * DT")
    },
    markov = {
      initial <- .nm_outcome_probability_choice(y, outcome$categories, outcome$initial)
      row_choice <- vapply(seq_along(outcome$categories), function(row) {
        .nm_outcome_probability_choice(y, outcome$categories,
                                       outcome$transition[row, ])
      }, character(1))
      transition <- .nm_outcome_probability_choice("PREV_DV", outcome$categories,
                                                    row_choice)
      paste0("log(pmax(ifelse(FIRST == 1, ", initial, ", ", transition,
             "), 1e-300))")
    },
    continuous_time_markov = {
      if (length(outcome$initial) != 2L) {
        .nm_stop("General continuous-time Markov outcomes use the compiled matrix-exponential path.")
      }
      initial <- .nm_outcome_probability_choice(y, outcome$categories, outcome$initial)
      q01 <- outcome$rates[1L, 2L]
      q10 <- outcome$rates[2L, 1L]
      total <- paste0("pmax(", q01, " + ", q10, ", 1e-12)")
      decay <- paste0("exp(-", total, " * DT)")
      p01 <- paste0(q01, " / ", total, " * (1 - ", decay, ")")
      p10 <- paste0(q10, " / ", total, " * (1 - ", decay, ")")
      p00 <- paste0("1 - (", p01, ")")
      p11 <- paste0("1 - (", p10, ")")
      from0 <- .nm_outcome_probability_choice(y, outcome$categories, c(p00, p01))
      from1 <- .nm_outcome_probability_choice(y, outcome$categories, c(p10, p11))
      transition <- .nm_outcome_probability_choice(
        "PREV_DV", outcome$categories, c(from0, from1)
      )
      paste0("log(pmax(ifelse(FIRST == 1, ", initial, ", ", transition,
             "), 1e-300))")
    }
  )
  if (!is.null(outcome$risk)) expression <- paste0("(", outcome$risk, ") * (", expression, ")")
  lines <- c(lines, paste0(prefix, "LOGLIK = ", expression))
  list(lines = lines, value = paste0(prefix, "LOGLIK"))
}

.nm_outcome_ctmc_hmm <- function(outcome) {
  states <- paste0("state_", format(outcome$categories, trim = TRUE, scientific = FALSE))
  states <- make.unique(states, sep = "_")
  n_state <- length(states)
  initial_names <- paste0("NM_CTM_I_", seq_len(n_state))
  generator_names <- matrix("", n_state, n_state)
  lines <- paste0(initial_names, " = ", outcome$initial)
  for (from in seq_len(n_state)) {
    for (to in seq_len(n_state)) {
      if (from == to) next
      name <- paste0("NM_CTM_Q_", from, "_", to)
      generator_names[from, to] <- name
      lines <- c(lines, paste0(name, " = ", outcome$rates[from, to]))
    }
  }
  emission_names <- paste0("NM_CTM_E_", seq_len(n_state))
  lines <- c(lines, vapply(seq_len(n_state), function(state) {
    paste0(emission_names[[state]], " = ifelse(DV == ",
           .nm_code_number(outcome$categories[[state]]), ", 1, 0)")
  }, character(1)))
  config <- nm_cthmm_config(
    states = states, initial = initial_names, generator = generator_names,
    emission = emission_names, initial_scale = "probability",
    rate_scale = "rate", emission_scale = "likelihood",
    by_dvid = !is.null(outcome$dvid)
  )
  config$observed_states <- TRUE
  config$categories <- outcome$categories
  list(config = config, error = paste(lines, collapse = "\n"))
}

.nm_outcome_error <- function(outcomes) {
  pieces <- lapply(seq_along(outcomes), function(index) {
    .nm_outcome_component_code(outcomes[[index]], index)
  })
  lines <- unlist(lapply(pieces, `[[`, "lines"), use.names = FALSE)
  if (length(outcomes) == 1L && is.null(outcomes[[1L]]$dvid)) {
    return(paste(c(lines, paste0("LOGLIK = ", pieces[[1L]]$value)), collapse = "\n"))
  }
  expression <- "-1e100"
  for (index in rev(seq_along(outcomes))) {
    expression <- paste0(
      "ifelse(DVID == ", .nm_code_number(outcomes[[index]]$dvid), ", ",
      pieces[[index]]$value, ", ", expression, ")"
    )
  }
  paste(c(lines, paste0("LOGLIK = ", expression)), collapse = "\n")
}

.nm_outcome_symbols <- function(outcomes) {
  if (is.null(outcomes)) return(character())
  candidates <- unlist(lapply(outcomes, function(value) c(
    value$prediction, value$scale, value$probabilities, value$dispersion,
    value$zero_probability, value$cause_hazards, value$initial,
    as.vector(value$transition), value$rates
  )), use.names = FALSE)
  candidates <- unique(candidates[grepl("^[A-Za-z][A-Za-z0-9_.]*$", candidates)])
  setdiff(candidates, c("F", "PRED", "IPRED", "DV", "DVID", "DT", "FIRST",
                        "PREV_DV", "PREV_TIME"))
}

.nm_outcome_resolve <- function(result, expression, theta, sigma, rows = seq_len(nrow(result))) {
  if (is.null(expression)) return(NULL)
  numeric_value <- suppressWarnings(as.numeric(expression))
  if (length(numeric_value) == 1L && is.finite(numeric_value)) {
    return(rep(numeric_value, length(rows)))
  }
  indexed <- regexec("^(THETA|SIGMA)\\(([0-9]+)\\)$", expression)
  matched <- regmatches(expression, indexed)[[1L]]
  if (length(matched)) {
    source <- if (matched[[2L]] == "THETA") theta else sigma
    index <- as.integer(matched[[3L]])
    if (index > length(source)) .nm_stop(expression, " is outside the fitted parameter vector.")
    return(rep(source[[index]], length(rows)))
  }
  column <- if (expression %in% c("F", "PRED", "IPRED")) "IPRED" else expression
  if (!column %in% names(result)) {
    .nm_stop("Outcome simulation requires generated model output `", expression, "`.")
  }
  as.numeric(result[[column]][rows])
}

.nm_outcome_rows <- function(result, outcome, include_mdv = TRUE) {
  rows <- result$EVID == 0L
  if (isTRUE(include_mdv)) rows <- rows & result$MDV == 0L
  if (!is.null(outcome$dvid)) {
    if (!"DVID" %in% names(result)) .nm_stop("Joint outcomes require a DVID column.")
    rows <- rows & as.numeric(result$DVID) == outcome$dvid
  }
  which(rows)
}

.nm_draw_categories <- function(probability) {
  probability <- as.matrix(probability)
  probability[!is.finite(probability) | probability < 0] <- 0
  totals <- rowSums(probability)
  if (any(totals <= 0)) .nm_stop("Outcome probabilities must have a positive row sum.")
  probability <- probability / totals
  uniform <- stats::runif(nrow(probability))
  cumulative <- t(apply(probability, 1L, cumsum))
  vapply(seq_len(nrow(probability)), function(row) {
    which(uniform[[row]] <= cumulative[row, ])[1L]
  }, integer(1))
}

.nm_simulate_first_class_outcomes <- function(model, result, theta, sigma) {
  outcomes <- model$OUTCOMES
  if (is.null(outcomes)) .nm_stop("The model has no first-class outcome generator.")
  for (outcome in outcomes) {
    rows <- .nm_outcome_rows(result, outcome, include_mdv = TRUE)
    if (!length(rows)) next
    family <- outcome$family
    prediction <- .nm_outcome_resolve(result, outcome$prediction, theta, sigma, rows)
    value <- switch(
      family,
      normal = stats::rnorm(length(rows), prediction,
                            pmax(.nm_outcome_resolve(result, outcome$scale, theta, sigma, rows), 0)),
      lognormal = stats::rlnorm(length(rows), log(pmax(prediction, 1e-300)),
                               pmax(.nm_outcome_resolve(result, outcome$scale, theta, sigma, rows), 0)),
      student_t = prediction +
        .nm_outcome_resolve(result, outcome$scale, theta, sigma, rows) *
        stats::rt(length(rows), df = outcome$df),
      bernoulli = stats::rbinom(length(rows), 1L, pmin(pmax(prediction, 0), 1)),
      categorical = {
        probability <- vapply(outcome$probabilities, function(symbol) {
          .nm_outcome_resolve(result, symbol, theta, sigma, rows)
        }, numeric(length(rows)))
        outcome$categories[.nm_draw_categories(probability)]
      },
      ordinal = {
        probability <- vapply(outcome$probabilities, function(symbol) {
          .nm_outcome_resolve(result, symbol, theta, sigma, rows)
        }, numeric(length(rows)))
        outcome$categories[.nm_draw_categories(probability)]
      },
      poisson = stats::rpois(length(rows), pmax(prediction, 0)),
      negative_binomial = stats::rnbinom(
        length(rows), mu = pmax(prediction, 0),
        size = pmax(.nm_outcome_resolve(result, outcome$dispersion, theta, sigma, rows), 1e-12)
      ),
      binomial = stats::rbinom(
        length(rows),
        size = pmax(round(.nm_outcome_resolve(result, outcome$trials, theta, sigma, rows)), 0),
        prob = pmin(pmax(prediction, 0), 1)
      ),
      zero_inflated_poisson = {
        structural <- stats::runif(length(rows)) <= pmin(pmax(
          .nm_outcome_resolve(result, outcome$zero_probability, theta, sigma, rows), 0), 1)
        count <- stats::rpois(length(rows), pmax(prediction, 0))
        count[structural] <- 0
        count
      },
      hurdle_poisson = {
        structural <- stats::runif(length(rows)) <= pmin(pmax(
          .nm_outcome_resolve(result, outcome$zero_probability, theta, sigma, rows), 0), 1)
        count <- stats::rpois(length(rows), pmax(prediction, 0))
        positive <- !structural
        while (any(positive & count == 0L)) {
          redraw <- positive & count == 0L
          count[redraw] <- stats::rpois(sum(redraw), pmax(prediction[redraw], 0))
          if (all(prediction[redraw] <= 0)) break
        }
        count[structural] <- 0
        count
      },
      tte = NULL,
      recurrent_event = NULL,
      competing_risks = NULL,
      markov = NULL,
      continuous_time_markov = NULL
    )
    if (!family %in% c("tte", "recurrent_event", "competing_risks", "markov",
                       "continuous_time_markov")) {
      result$DV[rows] <- value
      next
    }
    groups <- split(rows, interaction(result$.ID_INDEX[rows],
                                      if ("DVID" %in% names(result)) result$DVID[rows] else 1,
                                      drop = TRUE, lex.order = TRUE))
    if (family %in% c("tte", "recurrent_event", "competing_risks")) {
      result$DV[rows] <- 0
      for (group in groups) {
        group <- group[order(result$TIME[group])]
        active <- TRUE
        for (position in seq_along(group)) {
          row <- group[[position]]
          if (!active) next
          dt <- if (position == 1L) 0 else max(result$TIME[row] - result$TIME[group[[position - 1L]]], 0)
          if (family == "competing_risks") {
            hazards <- vapply(outcome$cause_hazards, function(symbol) {
              .nm_outcome_resolve(result, symbol, theta, sigma, row)
            }, numeric(1))
            hazards <- pmax(hazards, 0)
            total <- sum(hazards)
            if (total > 0 && stats::runif(1) <= 1 - exp(-total * dt)) {
              result$DV[row] <- outcome$categories[
                .nm_draw_categories(matrix(hazards, nrow = 1L))
              ]
              active <- FALSE
            }
          } else {
            hazard <- pmax(.nm_outcome_resolve(
              result, outcome$prediction, theta, sigma, row
            ), 0)
            if (stats::runif(1) <= 1 - exp(-hazard * dt)) {
              result$DV[row] <- outcome$event
              if (family == "tte") active <- FALSE
            }
          }
        }
      }
    } else {
      result$DV[rows] <- NA_real_
      for (group in groups) {
        group <- group[order(result$TIME[group])]
        previous <- NA_integer_
        for (position in seq_along(group)) {
          row <- group[[position]]
          if (position == 1L) {
            probability <- vapply(outcome$initial, function(symbol) {
              .nm_outcome_resolve(result, symbol, theta, sigma, row)
            }, numeric(1))
          } else if (family == "markov") {
            probability <- vapply(outcome$transition[previous, ], function(symbol) {
              .nm_outcome_resolve(result, symbol, theta, sigma, row)
            }, numeric(1))
          } else {
            n_state <- length(outcome$initial)
            generator <- matrix(0, n_state, n_state)
            for (from in seq_len(n_state)) {
              for (to in seq_len(n_state)) {
                if (from == to) next
                generator[from, to] <- pmax(.nm_outcome_resolve(
                  result, outcome$rates[from, to], theta, sigma, row
                ), 0)
              }
            }
            diag(generator) <- -rowSums(generator)
            dt <- max(result$TIME[row] - result$TIME[group[[position - 1L]]], 0)
            transition <- .liberation_matrix_exp_pade(generator, dt)
            probability <- pmax(transition[previous, ], 0)
          }
          previous <- .nm_draw_categories(matrix(probability, nrow = 1L))
          result$DV[row] <- outcome$categories[[previous]]
        }
      }
    }
  }
  result
}

.nm_validate_outcome_data <- function(model, data) {
  outcomes <- model$OUTCOMES
  if (is.null(outcomes)) return(invisible(data))
  observed <- data$EVID == 0L & data$MDV == 0L & is.finite(data$DV)
  for (outcome in outcomes) {
    rows <- observed
    if (!is.null(outcome$dvid)) {
      if (!"DVID" %in% names(data)) .nm_stop("Joint outcomes require a DVID column.")
      rows <- rows & as.numeric(data$DVID) == outcome$dvid
    }
    value <- data$DV[rows]
    if (!length(value)) next
    family <- outcome$family
    if (family == "bernoulli" && any(!value %in% c(0, 1))) {
      .nm_stop(outcome$name, " requires DV coded as zero/one.")
    }
    if (family %in% c("categorical", "ordinal", "markov", "continuous_time_markov") &&
        any(!value %in% outcome$categories)) {
      .nm_stop(outcome$name, " contains an undeclared category/state code.")
    }
    if (family %in% c("poisson", "negative_binomial", "binomial",
                      "zero_inflated_poisson", "hurdle_poisson") &&
        any(value < 0 | value != floor(value) | value > outcome$max_count)) {
      .nm_stop(outcome$name, " counts must be non-negative integers no larger than max_count.")
    }
    if (family %in% c("tte", "recurrent_event") && any(!value %in% c(0, outcome$event))) {
      .nm_stop(outcome$name, " must use zero for no event and the declared event code.")
    }
    if (family == "competing_risks" && any(!value %in% c(0, outcome$categories))) {
      .nm_stop(outcome$name, " contains an undeclared competing-risk cause code.")
    }
    if (family == "lognormal" && any(value <= 0)) {
      .nm_stop(outcome$name, " requires strictly positive observations.")
    }
  }
  invisible(data)
}

#' @export
print.nm_outcome <- function(x, ...) {
  cat("LibeRation outcome\n")
  cat("  ", x$name, ": ", x$family,
      if (!is.null(x$dvid)) paste0(" (DVID ", x$dvid, ")") else "", "\n", sep = "")
  invisible(x)
}
