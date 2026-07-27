.nm_mu_assignments <- function(model) {
  source <- paste(c(
    model$PRED %||% "", model$PRED_SOURCE %||% ""
  ), collapse = "\n")
  lines <- trimws(unlist(strsplit(
    gsub(";", "\n", source, fixed = TRUE), "\n", fixed = TRUE
  ), use.names = FALSE))
  pattern <- "^([A-Za-z][A-Za-z0-9_]*)\\s*(?:<-|=)\\s*(.+)$"
  parts <- Filter(length, regmatches(lines, regexec(pattern, lines, perl = TRUE)))
  if (!length(parts)) return(list())
  result <- lapply(parts, function(part) {
    parsed <- tryCatch(parse(text = part[[3L]])[[1L]], error = function(error) NULL)
    parsed
  })
  names(result) <- toupper(vapply(parts, `[[`, character(1), 2L))
  result[!vapply(result, is.null, logical(1))]
}

.nm_mu_inline <- function(expression, assignments, seen = character()) {
  if (is.symbol(expression)) {
    key <- toupper(as.character(expression))
    if (key %in% names(assignments) && !grepl("^MU_?[0-9]+$", key)) {
      if (key %in% seen) .nm_stop("Circular assignment in MU expression: ", key, ".")
      return(.nm_mu_inline(assignments[[key]], assignments, c(seen, key)))
    }
    return(expression)
  }
  if (!is.call(expression)) return(expression)
  result <- as.list(expression)
  if (length(result) > 1L) {
    result[-1L] <- lapply(result[-1L], .nm_mu_inline,
                          assignments = assignments, seen = seen)
  }
  as.call(result)
}

.nm_mu_theta_leaf <- function(expression) {
  if (is.symbol(expression)) {
    match <- regexec("^THETA_?([0-9]+)$", toupper(as.character(expression)),
                     perl = TRUE)
    part <- regmatches(toupper(as.character(expression)), match)[[1L]]
    if (length(part)) {
      return(list(index = as.integer(part[[2L]]), link = "identity"))
    }
    return(NULL)
  }
  if (!is.call(expression)) return(NULL)
  head <- toupper(as.character(expression[[1L]]))
  if (head == "THETA" && length(expression) == 2L) {
    index <- suppressWarnings(tryCatch(
      as.integer(expression[[2L]]), error = function(error) NA_integer_
    ))
    if (length(index) == 1L && !is.na(index) && index > 0L) {
      return(list(index = index, link = "identity"))
    }
  }
  if (head == "LOG" && length(expression) == 2L) {
    inner <- .nm_mu_theta_leaf(expression[[2L]])
    if (!is.null(inner) && identical(inner$link, "identity")) {
      inner$link <- "log"
      return(inner)
    }
  }
  NULL
}

.nm_mu_theta_indices <- function(expression) {
  leaf <- .nm_mu_theta_leaf(expression)
  if (!is.null(leaf)) return(leaf$index)
  if (!is.call(expression) || length(expression) == 1L) return(integer())
  unique(unlist(lapply(as.list(expression)[-1L], .nm_mu_theta_indices),
                use.names = FALSE))
}

.nm_mu_language <- function(operator, left, right = NULL) {
  if (is.null(right)) call(operator, left) else call(operator, left, right)
}

.nm_mu_affine_expression <- function(expression) {
  leaf <- .nm_mu_theta_leaf(expression)
  if (!is.null(leaf)) {
    return(list(
      valid = TRUE, intercept = 0,
      terms = list(list(
        index = leaf$index, link = leaf$link, coefficient = 1
      ))
    ))
  }
  theta <- .nm_mu_theta_indices(expression)
  if (!length(theta)) {
    return(list(valid = TRUE, intercept = expression, terms = list()))
  }
  if (!is.call(expression)) {
    return(list(valid = FALSE, reason = "unsupported THETA expression"))
  }
  operator <- as.character(expression[[1L]])
  upper <- toupper(operator)
  if (operator == "(" && length(expression) == 2L) {
    return(.nm_mu_affine_expression(expression[[2L]]))
  }
  if (operator %in% c("+", "-")) {
    left <- .nm_mu_affine_expression(expression[[2L]])
    if (!isTRUE(left$valid)) return(left)
    if (length(expression) == 2L) {
      if (operator == "+") return(left)
      left$intercept <- .nm_mu_language("-", left$intercept)
      left$terms <- lapply(left$terms, function(term) {
        term$coefficient <- .nm_mu_language("-", term$coefficient)
        term
      })
      return(left)
    }
    right <- .nm_mu_affine_expression(expression[[3L]])
    if (!isTRUE(right$valid)) return(right)
    if (operator == "-") {
      right$intercept <- .nm_mu_language("-", right$intercept)
      right$terms <- lapply(right$terms, function(term) {
        term$coefficient <- .nm_mu_language("-", term$coefficient)
        term
      })
    }
    return(list(
      valid = TRUE,
      intercept = .nm_mu_language("+", left$intercept, right$intercept),
      terms = c(left$terms, right$terms)
    ))
  }
  if (operator %in% c("*", "/") && length(expression) == 3L) {
    left <- .nm_mu_affine_expression(expression[[2L]])
    right <- .nm_mu_affine_expression(expression[[3L]])
    if (!isTRUE(left$valid) || !isTRUE(right$valid)) {
      return(list(valid = FALSE, reason = "unsupported product or ratio"))
    }
    if (operator == "/" && length(right$terms)) {
      return(list(
        valid = FALSE,
        reason = "a denominator depending on THETA is nonlinear"
      ))
    }
    if (length(left$terms) && length(right$terms)) {
      return(list(
        valid = FALSE,
        reason = "a product of THETA-dependent terms is nonlinear"
      ))
    }
    if (operator == "/") {
      scale <- right$intercept
      left$intercept <- .nm_mu_language("/", left$intercept, scale)
      left$terms <- lapply(left$terms, function(term) {
        term$coefficient <- .nm_mu_language("/", term$coefficient, scale)
        term
      })
      return(left)
    }
    variable <- if (length(left$terms)) left else right
    constant <- if (length(left$terms)) right$intercept else left$intercept
    variable$intercept <- .nm_mu_language(
      "*", variable$intercept, constant
    )
    variable$terms <- lapply(variable$terms, function(term) {
      term$coefficient <- .nm_mu_language(
        "*", term$coefficient, constant
      )
      term
    })
    return(variable)
  }
  list(
    valid = FALSE,
    reason = paste0("`", upper, "` is nonlinear in a THETA-dependent term")
  )
}

.nm_mu_eval_environment <- function(data, theta = numeric()) {
  environment <- list2env(as.list(data[1L, , drop = FALSE]),
                          parent = baseenv())
  environment$THETA <- function(index) theta[[as.integer(index)]]
  if (length(theta)) {
    for (index in seq_along(theta)) {
      assign(paste0("THETA_", index), theta[[index]], envir = environment)
      assign(paste0("THETA", index), theta[[index]], envir = environment)
    }
  }
  aliases <- list(
    LOG = base::log, EXP = base::exp, SQRT = base::sqrt, ABS = base::abs,
    MIN = base::min, MAX = base::max, LOG10 = base::log10,
    LOGIT = stats::qlogis, INVLOGIT = stats::plogis
  )
  list2env(aliases, envir = environment)
  environment
}

.nm_mu_eval_scalar <- function(expression, data, theta = numeric(),
                               label = "MU expression") {
  value <- tryCatch(
    eval(expression, envir = .nm_mu_eval_environment(data, theta)),
    error = identity
  )
  if (inherits(value, "error") || length(value) != 1L ||
      !is.numeric(value) || !is.finite(value)) {
    detail <- if (inherits(value, "error")) conditionMessage(value) else
      "the result was not one finite numeric scalar"
    .nm_stop("Unable to evaluate ", label, ": ", detail, ".")
  }
  as.numeric(value)
}

.nm_mu_theta_outside_references <- function(model, indices) {
  if (!length(indices)) return(integer())
  source <- paste(c(
    model$PRED %||% "", model$PRED_SOURCE %||% "",
    model$DES %||% "", model$ERROR %||% ""
  ), collapse = "\n")
  lines <- unlist(strsplit(gsub(";", "\n", source, fixed = TRUE), "\n",
                           fixed = TRUE), use.names = FALSE)
  lines <- lines[!grepl(
    "^\\s*MU_?[0-9]+\\s*(?:<-|=)", lines, ignore.case = TRUE, perl = TRUE
  )]
  vapply(indices, function(index) {
    pattern <- paste0(
      "\\bTHETA\\s*\\(\\s*", index, "\\s*\\)|\\bTHETA_?", index, "\\b"
    )
    any(grepl(pattern, lines, ignore.case = TRUE, perl = TRUE))
  }, logical(1))
}

.nm_mu_plan <- function(model, map = .nm_outer_map(model)) {
  mu <- model$MU
  if (is.null(mu) || !nrow(mu)) {
    return(list(
      mapped = FALSE, affine = FALSE, saem_eligible = FALSE,
      reason = "the model has no MU references", expressions = list(),
      theta = integer(), links = character(), outside = integer()
    ))
  }
  assignments <- .nm_mu_assignments(model)
  parsed <- lapply(mu$EXPRESSION, function(value) {
    expression <- tryCatch(parse(text = value)[[1L]], error = identity)
    if (inherits(expression, "error")) {
      .nm_stop("Unable to parse MU expression `", value, "`: ",
               conditionMessage(expression), ".")
    }
    .nm_mu_inline(expression, assignments)
  })
  affine <- lapply(parsed, .nm_mu_affine_expression)
  invalid <- which(!vapply(affine, function(value) isTRUE(value$valid), logical(1)))
  if (length(invalid)) {
    return(list(
      mapped = TRUE, affine = FALSE, saem_eligible = FALSE,
      reason = paste0(
        "MU_", mu$MU[[invalid[[1L]]]], " is not affine: ",
        affine[[invalid[[1L]]]]$reason
      ),
      expressions = parsed, affine_expressions = affine,
      theta = integer(), links = character(), outside = integer()
    ))
  }
  terms <- unlist(lapply(affine, `[[`, "terms"), recursive = FALSE)
  theta <- if (length(terms)) {
    unique(vapply(terms, `[[`, integer(1), "index"))
  } else integer()
  free <- intersect(theta, map$theta_free)
  links <- if (length(free)) {
    stats::setNames(vapply(free, function(index) {
      values <- unique(vapply(
        Filter(function(term) term$index == index, terms),
        `[[`, character(1), "link"
      ))
      if (length(values) != 1L) NA_character_ else values
    }, character(1)), as.character(free))
  } else character()
  inconsistent <- anyNA(links)
  outside_flag <- .nm_mu_theta_outside_references(model, free)
  outside <- free[outside_flag]
  priors <- model$LIK_CONFIG$priors
  prior_theta <- if (is.null(priors) || !nrow(priors)) integer() else {
    parameters <- toupper(as.character(priors$parameter))
    as.integer(sub("^THETA", "", parameters[grepl("^THETA[0-9]+$", parameters)]))
  }
  prior_linked <- intersect(free, prior_theta)
  iov <- as.integer(model$LIK_CONFIG$iov %||% 0L)
  reason <- if (!length(free)) {
    "no free THETA is represented by the MU equations"
  } else if (inconsistent) {
    "one THETA uses inconsistent MU link functions"
  } else if (length(outside)) {
    paste0(
      "THETA", paste(outside, collapse = ", THETA"),
      " is also used outside its MU equation"
    )
  } else if (length(prior_linked)) {
    paste0(
      "THETA", paste(prior_linked, collapse = ", THETA"),
      " has an explicit prior"
    )
  } else if (iov > 0L) {
    "legacy IOV expansion is not yet eligible for the closed-form MU M-step"
  } else if (!is.null(model$RE_CONFIG)) {
    "general random-effect mappings are not eligible for the closed-form MU M-step"
  } else {
    NULL
  }
  list(
    mapped = TRUE, affine = !inconsistent,
    saem_eligible = is.null(reason), reason = reason,
    expressions = parsed, affine_expressions = affine,
    theta = free, links = links, outside = outside,
    prior_theta = prior_linked
  )
}

.nm_mu_link <- function(theta, link) {
  if (link == "identity") return(theta)
  if (link == "log") {
    if (theta <= 0) .nm_stop("A log-linked MU THETA must remain positive.")
    return(log(theta))
  }
  .nm_stop("Unsupported MU link: ", link, ".")
}

.nm_mu_inverse_link <- function(beta, link) {
  if (link == "identity") return(beta)
  if (link == "log") return(exp(beta))
  .nm_stop("Unsupported MU link: ", link, ".")
}

.nm_mu_beta <- function(specialization, theta) {
  vapply(seq_along(specialization$theta), function(column) {
    index <- specialization$theta[[column]]
    .nm_mu_link(theta[[index]], specialization$links[[as.character(index)]])
  }, numeric(1))
}

.nm_mu_specialization <- function(context, map, enabled = TRUE) {
  plan <- .nm_mu_plan(context$model, map)
  plan$enabled <- isTRUE(enabled)
  plan$active <- FALSE
  plan$design <- list()
  plan$design_columns <- list()
  plan$covariate_design <- FALSE
  plan$offset <- matrix(0, context$n_subjects, context$n_eta)
  plan$cache <- new.env(parent = emptyenv())
  plan$cache$omega <- NULL
  plan$cache$precision <- NULL
  plan$cache$weighted_design <- NULL
  plan$cache$hessian <- NULL
  plan$cache$chol_hessian <- NULL
  plan$cache$system_calls <- 0L
  plan$cache$hits <- 0L
  plan$cache$misses <- 0L
  plan$subject_data <- lapply(context$subjects, function(evaluator) {
    evaluator$data[1L, , drop = FALSE]
  })
  plan$eta <- if (isTRUE(plan$mapped)) as.integer(context$model$MU$ETA) else integer()
  if (!isTRUE(enabled)) {
    plan$reason <- "MU estimator specialization was disabled"
    return(plan)
  }
  if (!isTRUE(plan$mapped)) return(plan)
  # General MU re-centring only requires evaluable expressions, not affinity.
  invisible(.nm_mu_values(plan, context$model$THETAS$Value))
  if (!isTRUE(plan$affine) || !length(plan$theta)) return(plan)
  p <- length(plan$theta)
  plan$design <- vector("list", context$n_subjects)
  for (subject in seq_len(context$n_subjects)) {
    data <- plan$subject_data[[subject]]
    design <- matrix(0, context$n_eta, p)
    offset <- numeric(context$n_eta)
    for (row in seq_len(nrow(context$model$MU))) {
      eta <- context$model$MU$ETA[[row]]
      decomposition <- plan$affine_expressions[[row]]
      offset[[eta]] <- .nm_mu_eval_scalar(
        decomposition$intercept, data, label = "MU intercept"
      )
      for (term in decomposition$terms) {
        coefficient <- .nm_mu_eval_scalar(
          term$coefficient, data, label = "MU coefficient"
        )
        column <- match(term$index, plan$theta)
        if (is.na(column)) {
          offset[[eta]] <- offset[[eta]] + coefficient * .nm_mu_link(
            context$model$THETAS$Value[[term$index]], term$link
          )
        } else {
          design[eta, column] <- design[eta, column] + coefficient
        }
      }
    }
    plan$design[[subject]] <- design
    plan$offset[subject, ] <- offset
  }
  stacked <- do.call(rbind, plan$design)
  if (qr(stacked, tol = 1e-10)$rank < p) {
    plan$reason <- "the MU design matrix is rank deficient"
    plan$saem_eligible <- FALSE
    return(plan)
  }
  plan$design_columns <- lapply(seq_len(p), function(column) {
    output <- matrix(0, context$n_subjects, context$n_eta)
    for (subject in seq_len(context$n_subjects)) {
      output[subject, ] <- plan$design[[subject]][, column]
    }
    output
  })
  plan$covariate_design <- any(vapply(
    plan$design_columns,
    function(design) {
      if (nrow(design) < 2L) return(FALSE)
      any(abs(sweep(design, 2L, design[1L, ], FUN = "-")) > 1e-12)
    },
    logical(1)
  ))
  plan$beta_bounds <- .nm_mu_beta_bounds(plan, context$model)
  plan$active <- TRUE
  plan
}

.nm_mu_values <- function(specialization, theta) {
  result <- matrix(
    0, length(specialization$subject_data),
    ncol(specialization$offset)
  )
  if (!isTRUE(specialization$mapped)) return(result)
  if (isTRUE(specialization$active) && length(specialization$theta)) {
    beta <- .nm_mu_beta(specialization, theta)
    result <- specialization$offset
    for (column in seq_along(beta)) {
      result <- result +
        specialization$design_columns[[column]] * beta[[column]]
    }
    return(result)
  }
  for (subject in seq_len(nrow(result))) {
    for (row in seq_along(specialization$expressions)) {
      result[subject, specialization$eta[[row]]] <- .nm_mu_eval_scalar(
        specialization$expressions[[row]],
        specialization$subject_data[[subject]], theta,
        label = paste0("MU_", specialization$eta[[row]])
      )
    }
  }
  result
}

.nm_mu_recenter_eta <- function(specialization, old_parameters,
                                new_parameters, eta) {
  if (!isTRUE(specialization$enabled) || !isTRUE(specialization$mapped) ||
      !ncol(eta)) return(eta)
  if (isTRUE(specialization$active) && length(specialization$theta)) {
    difference <- .nm_mu_beta(specialization, old_parameters$theta) -
      .nm_mu_beta(specialization, new_parameters$theta)
    adjustment <- matrix(0, nrow(eta), ncol(eta))
    for (column in seq_along(difference)) {
      adjustment <- adjustment +
        specialization$design_columns[[column]] * difference[[column]]
    }
    return(eta + adjustment)
  }
  old_mu <- .nm_mu_values(specialization, old_parameters$theta)
  new_mu <- .nm_mu_values(specialization, new_parameters$theta)
  eta + old_mu - new_mu
}

.nm_mu_gls_system <- function(specialization, context, parameters, eta) {
  if (!isTRUE(specialization$active) || !length(specialization$theta)) {
    return(NULL)
  }
  p <- length(specialization$theta)
  beta <- .nm_mu_beta(specialization, parameters$theta)
  centered <- eta
  for (column in seq_along(beta)) {
    centered <- centered +
      specialization$design_columns[[column]] * beta[[column]]
  }
  cache <- specialization$cache
  cache$system_calls <- cache$system_calls + 1L
  if (!is.null(cache$omega) && identical(cache$omega, parameters$omega)) {
    cache$hits <- cache$hits + 1L
  } else {
    cache$misses <- cache$misses + 1L
    covariance <- .nm_effect_covariance(
      context$model, context$subjects[[1L]]$data, parameters$omega
    )
    covariance_pd <- tryCatch(
      .nm_positive_definite(covariance, "MU random-effect covariance"),
      error = identity
    )
    if (inherits(covariance_pd, "error")) {
      return(list(valid = FALSE, reason = conditionMessage(covariance_pd)))
    }
    covariance_chol <- chol(covariance_pd$matrix)
    precision <- chol2inv(covariance_chol)
    weighted_design <- lapply(
      specialization$design_columns,
      function(design) design %*% precision
    )
    hessian <- matrix(0, p, p)
    for (left in seq_len(p)) {
      for (right in seq_len(left)) {
        value <- sum(
          weighted_design[[left]] *
            specialization$design_columns[[right]]
        )
        hessian[left, right] <- hessian[right, left] <- value
      }
    }
    hessian <- (hessian + t(hessian)) / 2
    if (qr(hessian, tol = 1e-10)$rank < p) {
      return(list(
        valid = FALSE, reason = "the MU design matrix is rank deficient"
      ))
    }
    hessian_pd <- tryCatch(
      .nm_positive_definite(hessian, "MU fixed-effect information"),
      error = identity
    )
    if (inherits(hessian_pd, "error")) {
      return(list(valid = FALSE, reason = conditionMessage(hessian_pd)))
    }
    cache$omega <- parameters$omega
    cache$precision <- precision
    cache$weighted_design <- weighted_design
    cache$hessian <- hessian_pd$matrix
    cache$chol_hessian <- chol(hessian_pd$matrix)
  }
  score <- vapply(seq_len(p), function(column) {
    sum(cache$weighted_design[[column]] * centered)
  }, numeric(1))
  mean <- backsolve(
    cache$chol_hessian,
    forwardsolve(t(cache$chol_hessian), score)
  )
  list(
    valid = TRUE, hessian = cache$hessian, score = score,
    mean = as.vector(mean),
    phi = centered + specialization$offset
  )
}

.nm_mu_beta_bounds <- function(specialization, model) {
  lower <- model$THETAS$LOWER %||% rep(-Inf, nrow(model$THETAS))
  upper <- model$THETAS$UPPER %||% rep(Inf, nrow(model$THETAS))
  result_lower <- result_upper <- numeric(length(specialization$theta))
  for (column in seq_along(specialization$theta)) {
    index <- specialization$theta[[column]]
    link <- specialization$links[[as.character(index)]]
    if (link == "identity") {
      result_lower[[column]] <- lower[[index]]
      result_upper[[column]] <- upper[[index]]
    } else {
      result_lower[[column]] <- if (is.finite(lower[[index]]) &&
          lower[[index]] > 0) log(lower[[index]]) else -Inf
      result_upper[[column]] <- if (is.finite(upper[[index]]) &&
          upper[[index]] > 0) log(upper[[index]]) else Inf
    }
  }
  list(lower = result_lower, upper = result_upper)
}

.nm_mu_gls_update <- function(specialization, context, parameters, eta) {
  system <- .nm_mu_gls_system(specialization, context, parameters, eta)
  if (is.null(system) || !isTRUE(system$valid)) return(system)
  bounds <- specialization$beta_bounds %||%
    .nm_mu_beta_bounds(specialization, context$model)
  beta <- system$mean
  if (any(beta < bounds$lower) || any(beta > bounds$upper)) {
    start <- pmin(pmax(.nm_mu_beta(specialization, parameters$theta),
                       bounds$lower), bounds$upper)
    optimized <- stats::optim(
      start,
      function(value) {
        as.numeric(crossprod(value, system$hessian %*% value) / 2 -
                     crossprod(system$score, value))
      },
      function(value) as.vector(system$hessian %*% value - system$score),
      method = "L-BFGS-B", lower = bounds$lower, upper = bounds$upper
    )
    if (optimized$convergence != 0L) {
      return(list(valid = FALSE, reason = "the bounded MU GLS update failed"))
    }
    beta <- optimized$par
  }
  candidate <- parameters
  for (column in seq_along(specialization$theta)) {
    index <- specialization$theta[[column]]
    candidate$theta[[index]] <- .nm_mu_inverse_link(
      beta[[column]], specialization$links[[as.character(index)]]
    )
  }
  list(
    valid = TRUE, parameters = candidate,
    eta = .nm_mu_recenter_eta(specialization, parameters, candidate, eta),
    beta = beta, hessian = system$hessian
  )
}

.nm_mu_log_beta_proposal <- function(beta, mean, hessian,
                                      specialization, theta) {
  difference <- beta - mean
  logdet <- as.numeric(determinant(hessian, logarithm = TRUE)$modulus)
  value <- 0.5 * logdet - length(beta) * log(2 * pi) / 2 -
    as.numeric(crossprod(difference, hessian %*% difference)) / 2
  # Proposal density is constructed in beta coordinates but the BAYES state
  # stores native THETAs.
  for (index in specialization$theta) {
    if (specialization$links[[as.character(index)]] == "log") {
      value <- value - log(theta[[index]])
    }
  }
  value
}

.nm_mu_bayes_proposal <- function(specialization, context, state, map,
                                  log_posterior, current_logp) {
  system <- .nm_mu_gls_system(
    specialization, context, state$parameters, state$eta
  )
  if (is.null(system) || !isTRUE(system$valid)) {
    return(list(state = state, accepted = FALSE, attempted = FALSE,
                reason = system$reason %||% "MU GLS system unavailable"))
  }
  covariance <- solve(system$hessian)
  beta <- as.vector(
    system$mean + t(chol(covariance)) %*%
      stats::rnorm(length(system$mean))
  )
  candidate_parameters <- state$parameters
  for (column in seq_along(specialization$theta)) {
    index <- specialization$theta[[column]]
    candidate_parameters$theta[[index]] <- .nm_mu_inverse_link(
      beta[[column]], specialization$links[[as.character(index)]]
    )
  }
  candidate_outer <- tryCatch(
    map$encode(candidate_parameters), error = function(error) NULL
  )
  if (is.null(candidate_outer) || !map$in_bounds(candidate_outer)) {
    return(list(state = state, accepted = FALSE, attempted = TRUE))
  }
  candidate <- .nm_bayes_state(
    map, candidate_outer,
    .nm_mu_recenter_eta(
      specialization, state$parameters, candidate_parameters, state$eta
    )
  )
  proposed_logp <- log_posterior(candidate)
  current_beta <- .nm_mu_beta(specialization, state$parameters$theta)
  log_q_current <- .nm_mu_log_beta_proposal(
    current_beta, system$mean, system$hessian,
    specialization, state$parameters$theta
  )
  log_q_candidate <- .nm_mu_log_beta_proposal(
    beta, system$mean, system$hessian,
    specialization, candidate_parameters$theta
  )
  accepted <- is.finite(proposed_logp) &&
    log(stats::runif(1)) < proposed_logp - current_logp +
      log_q_current - log_q_candidate
  list(
    state = if (accepted) candidate else state,
    accepted = accepted, attempted = TRUE,
    log_posterior = if (accepted) proposed_logp else current_logp
  )
}

.nm_mu_diagnostic <- function(specialization) {
  list(
    enabled = isTRUE(specialization$enabled),
    mapped = isTRUE(specialization$mapped),
    affine = isTRUE(specialization$affine),
    active = isTRUE(specialization$active),
    recenter_available = isTRUE(specialization$enabled) &&
      isTRUE(specialization$mapped),
    saem_eligible = isTRUE(specialization$saem_eligible),
    theta = as.integer(specialization$theta %||% integer()),
    links = unname(specialization$links %||% character()),
    reason = specialization$reason %||% NULL,
    gls_system_calls = as.integer(specialization$cache$system_calls %||% 0L),
    gls_cache_hits = as.integer(specialization$cache$hits %||% 0L),
    gls_cache_misses = as.integer(specialization$cache$misses %||% 0L),
    gls_vectorized = isTRUE(specialization$active),
    covariate_design = isTRUE(specialization$covariate_design)
  )
}
