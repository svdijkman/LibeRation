#' Define one latent chain in a factorial HMM
#'
#' @param name Chain label.
#' @param states State labels.
#' @param initial Initial-weight assignment names.
#' @param transition Row-by-column transition assignment names.
#' @param initial_scale,transition_scale Probability or log scale.
#' @export
nm_factorial_chain <- function(name, states, initial, transition,
                               initial_scale = c("probability", "log"),
                               transition_scale = c("probability", "log")) {
  name <- trimws(as.character(name)); states <- trimws(as.character(states))
  transition <- as.matrix(transition); initial <- trimws(as.character(initial))
  if (length(name) != 1L || is.na(name) || !nzchar(name) || length(states) < 2L ||
      anyNA(states) || any(!nzchar(states)) || anyDuplicated(states) ||
      length(initial) != length(states) ||
      !identical(dim(transition), rep(length(states), 2L))) {
    .nm_stop("A factorial chain requires a name, at least two states, and matching initial/transition assignments.")
  }
  structure(list(
    name = name, states = states, initial = initial,
    transition = matrix(as.character(transition), nrow(transition), ncol(transition)),
    initial_scale = match.arg(initial_scale),
    transition_scale = match.arg(transition_scale)
  ), class = "nm_factorial_chain")
}

#' Configure an exact factorial hidden Markov model
#'
#' Independent chain transitions are combined lazily into one sparse-compatible
#' joint-state HMM. State-conditional emission assignments are supplied for the
#' Cartesian joint states. Decoding retains joint probabilities and adds chain
#' marginals.
#'
#' @param ... [nm_factorial_chain()] declarations.
#' @param emission One emission assignment per joint state.
#' @param emission_scale Likelihood or log scale.
#' @param by_dvid Maintain independent sequences by DVID.
#' @param max_joint_states Safety limit for exact enumeration.
#' @param prefix Generated assignment prefix.
#' @export
nm_factorial_hmm_config <- function(..., emission,
                                    emission_scale = c("likelihood", "log"),
                                    by_dvid = TRUE, max_joint_states = 256L,
                                    prefix = "FHMM") {
  chains <- list(...)
  if (length(chains) == 1L && is.list(chains[[1L]]) &&
      !inherits(chains[[1L]], "nm_factorial_chain")) chains <- chains[[1L]]
  if (length(chains) < 2L || any(!vapply(chains, inherits, logical(1), "nm_factorial_chain"))) {
    .nm_stop("A factorial HMM requires at least two `nm_factorial_chain()` declarations.")
  }
  chain_names <- vapply(chains, `[[`, character(1), "name")
  if (anyDuplicated(chain_names)) .nm_stop("Factorial chain names must be unique.")
  dimensions <- vapply(chains, function(chain) length(chain$states), integer(1))
  count <- prod(dimensions)
  max_joint_states <- as.integer(max_joint_states)
  if (length(max_joint_states) != 1L || is.na(max_joint_states) ||
      max_joint_states < 2L || count > max_joint_states) {
    .nm_stop("The factorial HMM has ", count, " joint states; increase `max_joint_states` only after checking memory/runtime.")
  }
  grid <- do.call(expand.grid, c(lapply(dimensions, seq_len),
                                 KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
  names(grid) <- chain_names
  labels <- vapply(seq_len(nrow(grid)), function(row) {
    paste(vapply(seq_along(chains), function(index) {
      paste0(chain_names[[index]], "=", chains[[index]]$states[[grid[row, index]]])
    }, character(1)), collapse = "|")
  }, character(1))
  emission <- trimws(as.character(emission))
  if (length(emission) != count || anyNA(emission) || any(!nzchar(emission))) {
    .nm_stop("`emission` must contain one assignment per factorial joint state.")
  }
  prefix <- gsub("[^A-Za-z0-9_]", "_", toupper(as.character(prefix)[[1L]]))
  prob <- function(symbol, scale) if (scale == "log") paste0("exp(", symbol, ")") else symbol
  initial_names <- paste0(prefix, "_I", seq_len(count))
  transition_names <- matrix("", count, count)
  code <- character()
  for (joint in seq_len(count)) {
    factors <- vapply(seq_along(chains), function(index) {
      prob(chains[[index]]$initial[[grid[joint, index]]], chains[[index]]$initial_scale)
    }, character(1))
    code <- c(code, paste0(initial_names[[joint]], " = ", paste(factors, collapse = " * ")))
    for (target in seq_len(count)) {
      name <- paste0(prefix, "_T", joint, "_", target)
      factors <- vapply(seq_along(chains), function(index) {
        prob(chains[[index]]$transition[grid[joint, index], grid[target, index]],
             chains[[index]]$transition_scale)
      }, character(1))
      transition_names[joint, target] <- name
      code <- c(code, paste0(name, " = ", paste(factors, collapse = " * ")))
    }
  }
  config <- nm_hmm_config(
    states = labels, initial = initial_names, transition = transition_names,
    emission = emission, initial_scale = "probability",
    transition_scale = "probability", emission_scale = match.arg(emission_scale),
    by_dvid = by_dvid
  )
  attr(config, "generated_error") <- paste(code, collapse = "\n")
  attr(config, "factorial") <- list(chains = chains, grid = grid, labels = labels)
  class(config) <- c("nm_factorial_hmm_config", class(config))
  config
}

#' @export
print.nm_factorial_hmm_config <- function(x, ...) {
  metadata <- attr(x, "factorial", exact = TRUE)
  cat("LibeRation exact factorial hidden Markov model\n")
  cat("  chains:", paste(vapply(metadata$chains, `[[`, character(1), "name"), collapse = ", "),
      " joint states:", nrow(metadata$grid), "\n")
  invisible(x)
}

