#' Configure a finite-state hidden Markov model
#'
#' `nm_hmm_config()` maps assignments in a model's `$ERROR` block to the
#' initial-state weights, transition weights, and state-conditional emission
#' likelihoods used by LibeRation's compiled forward algorithm. Initial and
#' transition weights are normalized internally, so users may supply either
#' probabilities or positive unnormalised weights. Log weights are normalized
#' with a stable softmax. Emissions may be probabilities, densities, or log
#' likelihoods.
#'
#' @param states Unique labels for the hidden states.
#' @param initial Character vector naming one `$ERROR` assignment per state.
#' @param transition Character matrix naming `$ERROR` assignments for each
#'   from-state (row) to-state (column) transition.
#' @param emission Character vector naming one state-conditional `$ERROR`
#'   emission assignment per state.
#' @param initial_scale Whether `initial` assignments contain non-negative
#'   weights or log weights.
#' @param transition_scale Whether `transition` assignments contain
#'   non-negative weights or log weights.
#' @param emission_scale Whether `emission` assignments contain positive
#'   likelihoods/densities or log likelihoods.
#' @param by_dvid Maintain independent hidden-state sequences for each `DVID`
#'   within subject. If false, all observations within subject form one
#'   sequence.
#' @return A validated, serializable hidden Markov configuration.
#' @export
nm_hmm_config <- function(
    states,
    initial,
    transition,
    emission,
    initial_scale = c("probability", "log"),
    transition_scale = c("probability", "log"),
    emission_scale = c("likelihood", "log"),
    by_dvid = TRUE) {
  states <- trimws(as.character(states))
  if (length(states) < 2L || anyNA(states) || any(!nzchar(states)) ||
      anyDuplicated(states)) {
    .nm_stop("`states` must contain at least two unique, non-empty labels.")
  }
  names_vector <- function(value, expected, label) {
    value <- trimws(as.character(value))
    if (length(value) != expected || anyNA(value) || any(!nzchar(value))) {
      .nm_stop("`", label, "` must name exactly one non-empty $ERROR assignment per state.")
    }
    unname(value)
  }
  initial <- names_vector(initial, length(states), "initial")
  emission <- names_vector(emission, length(states), "emission")
  transition <- as.matrix(transition)
  if (!identical(dim(transition), rep(length(states), 2L))) {
    .nm_stop("`transition` must be a square character matrix with one row and column per state.")
  }
  transition <- matrix(
    names_vector(as.vector(transition), length(states)^2L, "transition"),
    nrow = length(states), ncol = length(states),
    dimnames = list(from = states, to = states)
  )
  initial_scale <- match.arg(initial_scale)
  transition_scale <- match.arg(transition_scale)
  emission_scale <- match.arg(emission_scale)
  if (length(by_dvid) != 1L || is.na(by_dvid)) {
    .nm_stop("`by_dvid` must be TRUE or FALSE.")
  }
  structure(
    list(
      version = 2L,
      transition_type = "discrete",
      states = states,
      initial = initial,
      transition = transition,
      generator = NULL,
      emission = emission,
      initial_scale = initial_scale,
      transition_scale = transition_scale,
      rate_scale = NULL,
      emission_scale = emission_scale,
      by_dvid = isTRUE(by_dvid)
    ),
    class = "nm_hmm_config"
  )
}

#' Configure a continuous-time hidden Markov model
#'
#' `nm_cthmm_config()` is the continuous-time counterpart of
#' [nm_hmm_config()]. Off-diagonal entries of `generator` name transition-rate
#' assignments in `$ERROR`; diagonal entries are generated as minus the sum of
#' the other entries in their row. At every observation the compiled engine
#' calculates `exp(Q * DT)` with a CppAD-compatible matrix exponential. This
#' supports irregular observation times, time-varying rates, exact gradients,
#' retrospective smoothing, and Viterbi decoding.
#'
#' @param states Unique labels for the hidden states.
#' @param initial Character vector naming one `$ERROR` initial-state assignment
#'   per state.
#' @param generator Square character matrix. Every off-diagonal element names
#'   the `$ERROR` assignment for the corresponding from-state to to-state rate.
#'   Diagonal elements are ignored and may be empty or `NA`.
#' @param emission Character vector naming one state-conditional `$ERROR`
#'   emission assignment per state.
#' @param initial_scale Whether `initial` contains probability weights or log
#'   weights.
#' @param rate_scale Whether off-diagonal generator assignments contain
#'   non-negative rates or unconstrained log rates.
#' @param emission_scale Whether `emission` contains likelihoods/densities or
#'   log likelihoods.
#' @param by_dvid Maintain independent state sequences for each `DVID` within
#'   subject.
#' @return A validated, serializable continuous-time hidden Markov
#'   configuration.
#' @export
nm_cthmm_config <- function(
    states,
    initial,
    generator,
    emission,
    initial_scale = c("probability", "log"),
    rate_scale = c("rate", "log"),
    emission_scale = c("likelihood", "log"),
    by_dvid = TRUE) {
  states <- trimws(as.character(states))
  if (length(states) < 2L || anyNA(states) || any(!nzchar(states)) ||
      anyDuplicated(states)) {
    .nm_stop("`states` must contain at least two unique, non-empty labels.")
  }
  names_vector <- function(value, expected, label) {
    value <- trimws(as.character(value))
    if (length(value) != expected || anyNA(value) || any(!nzchar(value))) {
      .nm_stop("`", label, "` must name exactly one non-empty $ERROR assignment per state.")
    }
    unname(value)
  }
  initial <- names_vector(initial, length(states), "initial")
  emission <- names_vector(emission, length(states), "emission")
  generator <- as.matrix(generator)
  if (!identical(dim(generator), rep(length(states), 2L))) {
    .nm_stop("`generator` must be a square character matrix with one row and column per state.")
  }
  generator <- matrix(
    trimws(as.character(generator)), nrow = length(states),
    ncol = length(states), dimnames = list(from = states, to = states)
  )
  diagonal <- row(generator) == col(generator)
  generator[diagonal] <- ""
  off_diagonal <- generator[!diagonal]
  if (anyNA(off_diagonal) || any(!nzchar(off_diagonal))) {
    .nm_stop("Every off-diagonal `generator` element must name a non-empty $ERROR rate assignment.")
  }
  initial_scale <- match.arg(initial_scale)
  rate_scale <- match.arg(rate_scale)
  emission_scale <- match.arg(emission_scale)
  if (length(by_dvid) != 1L || is.na(by_dvid)) {
    .nm_stop("`by_dvid` must be TRUE or FALSE.")
  }
  structure(
    list(
      version = 2L,
      transition_type = "continuous",
      states = states,
      initial = initial,
      transition = NULL,
      generator = generator,
      emission = emission,
      initial_scale = initial_scale,
      transition_scale = NULL,
      rate_scale = rate_scale,
      emission_scale = emission_scale,
      by_dvid = isTRUE(by_dvid)
    ),
    class = c("nm_cthmm_config", "nm_hmm_config")
  )
}

#' Configure a discrete hidden semi-Markov model
#'
#' Expands explicit state-duration distributions into a sparse finite-state
#' representation consumed by the compiled HMM forward, smoother, Viterbi, and
#' automatic-differentiation paths. Duration is measured in observation steps.
#'
#' @param states Original semi-Markov state labels.
#' @param initial,emission One `$ERROR` assignment per state.
#' @param transition Square matrix of between-state transition assignments.
#' @param dwell Matrix with one row per state and one column per possible dwell
#'   duration.
#' @param initial_scale,transition_scale,dwell_scale,emission_scale Probability
#'   or log scales.
#' @param by_dvid Maintain independent sequences by DVID.
#' @param prefix Prefix for generated expanded-state assignments.
#' @export
nm_hsmm_config <- function(
    states, initial, transition, dwell, emission,
    initial_scale = c("probability", "log"),
    transition_scale = c("probability", "log"),
    dwell_scale = c("probability", "log"),
    emission_scale = c("likelihood", "log"), by_dvid = TRUE,
    prefix = "HSMM") {
  states <- trimws(as.character(states))
  n_states <- length(states)
  if (n_states < 2L || anyNA(states) || any(!nzchar(states)) || anyDuplicated(states)) {
    .nm_stop("`states` must contain at least two unique labels.")
  }
  initial <- trimws(as.character(initial)); emission <- trimws(as.character(emission))
  transition <- as.matrix(transition); dwell <- as.matrix(dwell)
  if (length(initial) != n_states || length(emission) != n_states ||
      !identical(dim(transition), c(n_states, n_states)) ||
      nrow(dwell) != n_states || ncol(dwell) < 1L ||
      anyNA(c(initial, transition, dwell, emission)) ||
      any(!nzchar(trimws(c(initial, transition, dwell, emission))))) {
    .nm_stop("HSMM initial, transition, dwell, and emission assignments have inconsistent dimensions.")
  }
  initial_scale <- match.arg(initial_scale)
  transition_scale <- match.arg(transition_scale)
  dwell_scale <- match.arg(dwell_scale)
  emission_scale <- match.arg(emission_scale)
  prefix <- gsub("[^A-Za-z0-9_]", "_", toupper(as.character(prefix)[[1L]]))
  durations <- ncol(dwell)
  expanded <- as.vector(outer(states, seq_len(durations), paste, sep = "@"))
  mapping <- expand.grid(state = seq_len(n_states), duration = seq_len(durations))
  # `outer` varies the first argument fastest; align the explicit map.
  mapping <- mapping[order(mapping$duration, mapping$state), , drop = FALSE]
  count <- nrow(mapping)
  initial_names <- paste0(prefix, "_I", seq_len(count))
  emission_names <- paste0(prefix, "_E", seq_len(count))
  transition_names <- matrix(paste0(prefix, "_ZERO"), count, count)
  code <- paste0(prefix, "_ZERO = ",
                 if (transition_scale == "log") "-1e100" else "0")
  combine <- function(...) {
    values <- c(...)
    if (all(c(initial_scale, dwell_scale) == "log") ||
        all(c(transition_scale, dwell_scale) == "log")) {
      paste(values, collapse = " + ")
    } else paste(values, collapse = " * ")
  }
  for (index in seq_len(count)) {
    state <- mapping$state[[index]]; duration <- mapping$duration[[index]]
    initial_expression <- if (initial_scale == dwell_scale) {
      if (initial_scale == "log") paste(initial[[state]], dwell[state, duration], sep = " + ")
      else paste(initial[[state]], dwell[state, duration], sep = " * ")
    } else {
      left <- if (initial_scale == "log") paste0("exp(", initial[[state]], ")") else initial[[state]]
      right <- if (dwell_scale == "log") paste0("exp(", dwell[state, duration], ")") else dwell[state, duration]
      paste(left, right, sep = " * ")
    }
    if (initial_scale == "log" && dwell_scale != "log") initial_expression <- paste0("log(", initial_expression, ")")
    code <- c(code, paste0(initial_names[[index]], " = ", initial_expression),
              paste0(emission_names[[index]], " = ", emission[[state]]))
    if (duration > 1L) {
      target <- which(mapping$state == state & mapping$duration == duration - 1L)
      name <- paste0(prefix, "_T", index, "_", target)
      transition_names[index, target] <- name
      code <- c(code, paste0(name, " = ", if (transition_scale == "log") "0" else "1"))
    } else {
      for (next_state in seq_len(n_states)) for (next_duration in seq_len(durations)) {
        target <- which(mapping$state == next_state & mapping$duration == next_duration)
        transition_value <- transition[state, next_state]
        dwell_value <- dwell[next_state, next_duration]
        expression <- if (transition_scale == dwell_scale) {
          if (transition_scale == "log") paste(transition_value, dwell_value, sep = " + ")
          else paste(transition_value, dwell_value, sep = " * ")
        } else {
          left <- if (transition_scale == "log") paste0("exp(", transition_value, ")") else transition_value
          right <- if (dwell_scale == "log") paste0("exp(", dwell_value, ")") else dwell_value
          paste(left, right, sep = " * ")
        }
        if (transition_scale == "log" && dwell_scale != "log") expression <- paste0("log(", expression, ")")
        name <- paste0(prefix, "_T", index, "_", target)
        transition_names[index, target] <- name
        code <- c(code, paste0(name, " = ", expression))
      }
    }
  }
  config <- nm_hmm_config(
    states = expanded, initial = initial_names, transition = transition_names,
    emission = emission_names, initial_scale = initial_scale,
    transition_scale = transition_scale, emission_scale = emission_scale,
    by_dvid = by_dvid
  )
  attr(config, "generated_error") <- paste(code, collapse = "\n")
  attr(config, "semi_markov") <- list(
    states = states, mapping = mapping, max_duration = durations
  )
  class(config) <- c("nm_hsmm_config", class(config))
  config
}

#' @export
print.nm_hsmm_config <- function(x, ...) {
  metadata <- attr(x, "semi_markov", exact = TRUE)
  cat("LibeRation hidden semi-Markov model\n")
  cat("  states:", paste(metadata$states, collapse = ", "),
      " maximum duration:", metadata$max_duration, "observation steps\n")
  invisible(x)
}

.nm_hmm_config <- function(config) {
  if (is.null(config)) return(NULL)
  if (inherits(config, "nm_hmm_config")) return(config)
  if (!is.list(config) || is.null(names(config))) {
    .nm_stop("HMM_CONFIG must be created by `nm_hmm_config()` or be a named list.")
  }
  config$version <- NULL
  transition_type <- config$transition_type %||% if (!is.null(config$generator)) {
    "continuous"
  } else "discrete"
  config$transition_type <- NULL
  if (identical(transition_type, "continuous")) {
    config$transition <- NULL
    config$transition_scale <- NULL
    do.call(nm_cthmm_config, config)
  } else {
    config$generator <- NULL
    config$rate_scale <- NULL
    do.call(nm_hmm_config, config)
  }
}

#' @export
print.nm_hmm_config <- function(x, ...) {
  continuous <- identical(x$transition_type %||% "discrete", "continuous")
  observed <- isTRUE(x$observed_states)
  cat("LibeRation ", if (continuous) "continuous-time " else "",
      if (observed) "observed Markov model\n" else "hidden Markov model\n", sep = "")
  cat("  states:", paste(x$states, collapse = ", "), "\n")
  cat("  initial:", x$initial_scale)
  if (continuous) cat(" rates:", x$rate_scale) else
    cat(" transition:", x$transition_scale)
  cat(" emission:", x$emission_scale, " by DVID:", x$by_dvid, "\n")
  invisible(x)
}

#' Decode hidden states
#'
#' Runs the same scaled forward algorithm used by estimation and optionally a
#' retrospective scaled forward-backward smoother or log-domain Viterbi
#' decoder. Filtering conditions on observations up to the current record,
#' smoothing conditions on the complete sequence, and Viterbi returns the
#' single most probable joint state path. All results are conditional on the
#' supplied or fitted parameters and ETAs.
#'
#' @param object An [nm_model()], compiled [NMEngine], or fitted `nm_fit`.
#' @param data Observation data. May be omitted for an `nm_fit` to use its
#'   estimation dataset.
#' @param type For fitted models, use individual ETAs or population ETAs.
#' @param theta,eta,sigma Optional parameters for an unfitted model.
#' @param method Decoding result to return. `"filtered"`, `"smoothed"`, and
#'   `"viterbi"` expose the selected state through `HMM_STATE`; `"all"`
#'   returns explicitly prefixed columns for all three decoders.
#' @return The input data augmented with state, probability, and
#'   `HMM_ROW_NLL` columns. The total log likelihood and per-sequence Viterbi
#'   summaries are stored in the `log_likelihood` and `sequence_summary`
#'   attributes.
#' @export
nm_hmm_decode <- function(object, data = NULL,
                          type = c("individual", "population"),
                          theta = NULL, eta = NULL, sigma = NULL,
                          method = c("filtered", "smoothed", "viterbi", "all")) {
  type <- match.arg(type)
  method <- match.arg(method)
  fit <- if (inherits(object, "nm_fit")) object else NULL
  if (!is.null(fit)) {
    model <- fit$model
    data <- .nm_engine_data(model, data %||% fit$data)
    theta <- theta %||% fit$theta
    sigma <- sigma %||% fit$sigma
    eta <- eta %||% .nm_fit_eta_for_data(fit, data, type)
    engine <- nm_compile(model)
  } else {
    engine <- if (inherits(object, "NMEngine")) object else nm_compile(object)
    if (is.null(data)) .nm_stop("`data` is required when `object` is not an nm_fit.")
    data <- .nm_engine_data(engine$model, data)
    theta <- theta %||% engine$model$THETAS$Value
    sigma <- sigma %||% engine$model$SIGMAS$Value
  }
  decoded <- engine$hmm_filter(
    data, theta = theta, eta = eta, sigma = sigma
  )
  semi_markov <- attr(engine$model$HMM_CONFIG, "semi_markov", exact = TRUE)
  if (!is.null(semi_markov)) {
    mapping <- semi_markov$mapping
    aggregate_probabilities <- function(value) {
      output <- matrix(0, nrow(value), length(semi_markov$states))
      for (state in seq_along(semi_markov$states)) {
        output[, state] <- rowSums(value[, mapping$state == state, drop = FALSE])
      }
      colnames(output) <- semi_markov$states
      output
    }
    decoded$filtered <- aggregate_probabilities(decoded$filtered)
    decoded$smoothed <- aggregate_probabilities(decoded$smoothed)
    decoded$filtered_state <- max.col(decoded$filtered, ties.method = "first")
    decoded$smoothed_state <- max.col(decoded$smoothed, ties.method = "first")
    decoded$viterbi_state <- mapping$state[decoded$viterbi_state]
    decoded$filtered_state_label <- semi_markov$states[decoded$filtered_state]
    decoded$smoothed_state_label <- semi_markov$states[decoded$smoothed_state]
    decoded$viterbi_state_label <- semi_markov$states[decoded$viterbi_state]
    decoded$states <- semi_markov$states
  }
  factorial <- attr(engine$model$HMM_CONFIG, "factorial", exact = TRUE)
  factorial_result <- NULL
  if (!is.null(factorial)) {
    marginalize <- function(probabilities, chain_index, state_index) {
      rowSums(probabilities[, factorial$grid[[chain_index]] == state_index, drop = FALSE])
    }
    factorial_result <- lapply(seq_along(factorial$chains), function(chain_index) {
      chain <- factorial$chains[[chain_index]]
      filtered <- vapply(seq_along(chain$states), function(state_index) {
        marginalize(decoded$filtered, chain_index, state_index)
      }, numeric(nrow(decoded$filtered)))
      smoothed <- vapply(seq_along(chain$states), function(state_index) {
        marginalize(decoded$smoothed, chain_index, state_index)
      }, numeric(nrow(decoded$smoothed)))
      if (!is.matrix(filtered)) filtered <- matrix(filtered, ncol = length(chain$states))
      if (!is.matrix(smoothed)) smoothed <- matrix(smoothed, ncol = length(chain$states))
      colnames(filtered) <- colnames(smoothed) <- chain$states
      list(
        name = chain$name, states = chain$states,
        filtered = filtered, smoothed = smoothed,
        filtered_state = max.col(filtered, ties.method = "first"),
        smoothed_state = max.col(smoothed, ties.method = "first"),
        viterbi_state = factorial$grid[[chain_index]][decoded$viterbi_state]
      )
    })
    names(factorial_result) <- vapply(factorial_result, `[[`, character(1), "name")
  }
  output <- as.data.frame(data)
  output$HMM_ROW_NLL <- as.numeric(decoded$row_nll)
  state_names <- make.names(as.character(decoded$states), unique = TRUE)
  add_probabilities <- function(output, values, prefix) {
    for (index in seq_along(state_names)) {
      output[[paste0(prefix, state_names[[index]])]] <- values[, index]
    }
    output
  }
  if (method == "all") {
    output$HMM_FILTER_STATE_INDEX <- as.integer(decoded$filtered_state)
    output$HMM_FILTER_STATE <- as.character(decoded$filtered_state_label)
    output$HMM_SMOOTH_STATE_INDEX <- as.integer(decoded$smoothed_state)
    output$HMM_SMOOTH_STATE <- as.character(decoded$smoothed_state_label)
    output$HMM_VITERBI_STATE_INDEX <- as.integer(decoded$viterbi_state)
    output$HMM_VITERBI_STATE <- as.character(decoded$viterbi_state_label)
    output <- add_probabilities(output, decoded$filtered, "HMM_FILTER_PROB_")
    output <- add_probabilities(output, decoded$smoothed, "HMM_SMOOTH_PROB_")
    if (!is.null(factorial_result)) {
      for (chain in factorial_result) {
        prefix <- paste0("FHMM_", make.names(chain$name), "_")
        output[[paste0(prefix, "FILTER_STATE")]] <- chain$states[chain$filtered_state]
        output[[paste0(prefix, "SMOOTH_STATE")]] <- chain$states[chain$smoothed_state]
        output[[paste0(prefix, "VITERBI_STATE")]] <- chain$states[chain$viterbi_state]
        for (state_index in seq_along(chain$states)) {
          state <- make.names(chain$states[[state_index]])
          output[[paste0(prefix, "FILTER_PROB_", state)]] <- chain$filtered[, state_index]
          output[[paste0(prefix, "SMOOTH_PROB_", state)]] <- chain$smoothed[, state_index]
        }
      }
    }
  } else {
    state <- switch(
      method,
      filtered = list(index = decoded$filtered_state,
                      label = decoded$filtered_state_label,
                      probability = decoded$filtered),
      smoothed = list(index = decoded$smoothed_state,
                      label = decoded$smoothed_state_label,
                      probability = decoded$smoothed),
      viterbi = list(index = decoded$viterbi_state,
                     label = decoded$viterbi_state_label,
                     probability = NULL)
    )
    output$HMM_STATE_INDEX <- as.integer(state$index)
    output$HMM_STATE <- as.character(state$label)
    if (!is.null(state$probability)) {
      output <- add_probabilities(output, state$probability, "HMM_PROB_")
    }
    if (!is.null(factorial_result)) {
      for (chain in factorial_result) {
        prefix <- paste0("FHMM_", make.names(chain$name), "_")
        chain_index <- switch(method, filtered = chain$filtered_state,
                              smoothed = chain$smoothed_state,
                              viterbi = chain$viterbi_state)
        output[[paste0(prefix, "STATE")]] <- chain$states[chain_index]
        probabilities <- switch(method, filtered = chain$filtered,
                                smoothed = chain$smoothed, viterbi = NULL)
        if (!is.null(probabilities)) {
          for (state_index in seq_along(chain$states)) {
            output[[paste0(prefix, "PROB_", make.names(chain$states[[state_index]]))]] <-
              probabilities[, state_index]
          }
        }
      }
    }
  }
  attr(output, "states") <- as.character(decoded$states)
  attr(output, "method") <- method
  attr(output, "eta_type") <- type
  attr(output, "log_likelihood") <- as.numeric(decoded$log_likelihood)
  attr(output, "factorial_chains") <- factorial_result
  sequence_summary <- as.data.frame(decoded$sequence_summary)
  if (nrow(sequence_summary) && all(c(".ID_INDEX", "ID") %in% names(data))) {
    sequence_summary$ID <- data$ID[
      match(sequence_summary$.ID_INDEX, data$.ID_INDEX)
    ]
  }
  if (isTRUE(engine$model$HMM_CONFIG$by_dvid)) {
    sequence_summary$DVID <- sequence_summary$HMM_SEQUENCE
  }
  attr(output, "sequence_summary") <- sequence_summary
  class(output) <- c("nm_hmm_decode", class(output))
  output
}
