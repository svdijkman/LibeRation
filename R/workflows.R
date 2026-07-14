.nm_fit_native_parameters <- function(fit) {
  values <- c(fit$theta, fit$sigma, fit$omega)
  names(values) <- .nm_parameter_names(fit$theta, fit$sigma, fit$omega)
  values
}

.nm_profile_fixed_model <- function(fit, parameter, value) {
  model <- fit$model
  model$THETAS$Value <- fit$theta
  model$SIGMAS$Value <- fit$sigma
  model$OMEGAS$Value <- fit$omega
  match <- regexec("^(THETA|SIGMA|OMEGA)([0-9]+)$", toupper(parameter), perl = TRUE)
  parts <- regmatches(toupper(parameter), match)[[1L]]
  if (!length(parts)) .nm_stop("Unknown profile parameter: ", parameter, ".")
  table_name <- paste0(parts[[2L]], "S")
  index <- as.integer(parts[[3L]])
  table <- model[[table_name]]
  if (index < 1L || index > nrow(table)) .nm_stop("Unknown profile parameter: ", parameter, ".")
  if (parts[[2L]] == "OMEGA" && any(model$OMEGAS$ROW != model$OMEGAS$COL)) {
    .nm_stop("Profiling individual elements of a correlated OMEGA is not yet supported.")
  }
  table$Value[[index]] <- as.numeric(value)
  table$FIX[[index]] <- TRUE
  .nm_model_rebuild(model, stats::setNames(list(table), table_name))
}

#' Profile-likelihood parameter uncertainty
#'
#' Each grid point fixes one native-scale parameter and re-estimates every
#' remaining free parameter with the fit's estimation method. Confidence limits
#' use the likelihood-ratio cutoff on the `-2 log likelihood` objective scale.
#'
#' @param fit An `nm_fit`.
#' @param parameters Native parameter names such as `THETA1` or `SIGMA1`.
#' @param points Odd number of grid points per parameter.
#' @param span Half-width in standard errors (or fallback parameter scales).
#' @param level Confidence level.
#' @param ... Controls passed to [nm_est()].
#' @return An `nm_profile` containing grid fits and confidence intervals.
#' @export
nm_profile <- function(fit, parameters = NULL, points = 9L,
                       span = 3, level = 0.95, ...) {
  if (!inherits(fit, "nm_fit")) .nm_stop("`fit` must be an nm_fit.")
  points <- as.integer(points)
  if (length(points) != 1L || is.na(points) || points < 3L) .nm_stop("`points` must be at least 3.")
  if (points %% 2L == 0L) points <- points + 1L
  span <- as.numeric(span)
  level <- as.numeric(level)
  if (!is.finite(span) || span <= 0 || !is.finite(level) || level <= 0 || level >= 1) {
    .nm_stop("`span` must be positive and `level` must lie between zero and one.")
  }
  estimates <- .nm_fit_native_parameters(fit)
  free <- c(!fit$model$THETAS$FIX, !fit$model$SIGMAS$FIX, !fit$model$OMEGAS$FIX)
  names(free) <- names(estimates)
  if (any(fit$model$OMEGAS$ROW != fit$model$OMEGAS$COL)) {
    free[grepl("^OMEGA", names(free))] <- FALSE
  }
  if (is.null(parameters)) parameters <- names(estimates)[free]
  parameters <- unique(toupper(as.character(parameters)))
  if (!length(parameters) || length(setdiff(parameters, names(estimates)))) {
    .nm_stop("`parameters` contains unknown native parameter names.")
  }
  standard_errors <- fit$covariance$se %||% numeric()
  rows <- list()
  cutoff <- stats::qchisq(level, df = 1L)
  for (parameter in parameters) {
    estimate <- unname(estimates[[parameter]])
    se <- if (parameter %in% names(standard_errors)) unname(standard_errors[[parameter]]) else NA_real_
    if (!is.finite(se) || se <= 0) se <- max(abs(estimate) * 0.2, 0.05)
    grid <- seq(estimate - span * se, estimate + span * se, length.out = points)
    if (grepl("^(SIGMA|OMEGA)", parameter)) grid <- pmax(grid, .Machine$double.eps)
    grid <- sort(unique(c(grid, estimate)))
    for (value in grid) {
      if (isTRUE(all.equal(value, estimate, tolerance = 1e-12))) {
        profiled <- fit
      } else {
        fixed_model <- .nm_profile_fixed_model(fit, parameter, value)
        profiled <- tryCatch(
          nm_est(fixed_model, fit$data, method = fit$method, ...),
          error = identity
        )
      }
      rows[[length(rows) + 1L]] <- data.frame(
        parameter = parameter, value = value,
        objective = if (inherits(profiled, "nm_fit")) profiled$objective else NA_real_,
        convergence = if (inherits(profiled, "nm_fit")) profiled$convergence else NA_integer_,
        error = if (inherits(profiled, "error")) conditionMessage(profiled) else "",
        stringsAsFactors = FALSE
      )
    }
  }
  grid <- do.call(rbind, rows)
  intervals <- do.call(rbind, lapply(split(grid, grid$parameter), function(frame) {
    frame$delta <- frame$objective - min(frame$objective, na.rm = TRUE)
    accepted <- frame[is.finite(frame$delta) & frame$delta <= cutoff, , drop = FALSE]
    data.frame(
      parameter = frame$parameter[[1L]], estimate = estimates[[frame$parameter[[1L]]]],
      lower = if (nrow(accepted)) min(accepted$value) else NA_real_,
      upper = if (nrow(accepted)) max(accepted$value) else NA_real_,
      level = level, cutoff = cutoff,
      lower_truncated = nrow(accepted) && min(accepted$value) == min(frame$value),
      upper_truncated = nrow(accepted) && max(accepted$value) == max(frame$value),
      stringsAsFactors = FALSE
    )
  }))
  grid$delta <- stats::ave(grid$objective, grid$parameter,
                           FUN = function(value) value - min(value, na.rm = TRUE))
  structure(list(grid = grid, intervals = intervals, level = level, cutoff = cutoff,
                 method = fit$method), class = "nm_profile")
}

.nm_scm_candidates <- function(candidates, data) {
  candidates <- as.data.frame(candidates, stringsAsFactors = FALSE)
  names(candidates) <- tolower(names(candidates))
  required <- c("parameter", "covariate")
  if (length(setdiff(required, names(candidates))) || !nrow(candidates)) {
    .nm_stop("`candidates` requires non-empty `parameter` and `covariate` columns.")
  }
  if (!"form" %in% names(candidates)) candidates$form <- "continuous"
  if (!"reference" %in% names(candidates)) candidates$reference <- NA_character_
  if (!"category" %in% names(candidates)) candidates$category <- NA_character_
  if (!"scale" %in% names(candidates)) candidates$scale <- NA_real_
  if (!"initial" %in% names(candidates)) candidates$initial <- 0
  candidates$parameter <- toupper(trimws(candidates$parameter))
  candidates$covariate <- trimws(candidates$covariate)
  candidates$form <- tolower(trimws(candidates$form))
  if (any(!candidates$covariate %in% names(data))) .nm_stop("Every SCM covariate must be present in the fitted data.")
  if (any(!candidates$form %in% c("continuous", "power", "categorical"))) {
    .nm_stop("SCM forms are `continuous`, `power`, or `categorical`.")
  }
  candidates$id <- make.unique(paste(candidates$parameter, candidates$covariate,
                                     candidates$form, sep = ":"))
  candidates
}

.nm_scm_add_relationship <- function(model, candidate, data) {
  parameter <- candidate$parameter[[1L]]
  covariate <- candidate$covariate[[1L]]
  lines <- unlist(strsplit(gsub(";", "\n", model$PRED, fixed = TRUE), "\n", fixed = TRUE))
  pattern <- paste0("^\\s*", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", parameter),
                    "\\s*(?:<-|=)\\s*(.*)$")
  index <- grep(pattern, lines, ignore.case = TRUE, perl = TRUE)
  if (length(index) != 1L) {
    .nm_stop("SCM requires exactly one assignment to parameter `", parameter, "` in PRED.")
  }
  match <- regexec(pattern, lines[[index]], ignore.case = TRUE, perl = TRUE)
  rhs <- trimws(regmatches(lines[[index]], match)[[1L]][[2L]])
  theta_number <- nrow(model$THETAS) + 1L
  values <- data[[covariate]]
  form <- candidate$form[[1L]]
  reference <- suppressWarnings(as.numeric(candidate$reference[[1L]]))
  if (!is.finite(reference)) reference <- stats::median(as.numeric(values), na.rm = TRUE)
  if (form == "categorical") {
    levels <- unique(as.character(values[!is.na(values)]))
    category <- as.character(candidate$category[[1L]])
    if (is.na(category) || !nzchar(category)) {
      reference_level <- as.character(candidate$reference[[1L]])
      alternatives <- setdiff(levels, reference_level)
      category <- if (length(alternatives)) alternatives[[1L]] else levels[[1L]]
    }
    numeric_category <- suppressWarnings(as.numeric(category))
    if (!is.finite(numeric_category)) .nm_stop("Categorical SCM currently requires numeric category values.")
    effect <- paste0("exp(THETA(", theta_number, ")*ifelse(", covariate,
                     "==", format(numeric_category, digits = 15), ",1,0))")
  } else if (form == "power") {
    if (any(as.numeric(values) <= 0, na.rm = TRUE) || !is.finite(reference) || reference <= 0) {
      .nm_stop("Power covariate relationships require positive covariate and reference values.")
    }
    effect <- paste0("(", covariate, "/", format(reference, digits = 15), ")^THETA(", theta_number, ")")
  } else {
    scale <- suppressWarnings(as.numeric(candidate$scale[[1L]]))
    if (!is.finite(scale) || scale <= 0) scale <- stats::sd(as.numeric(values), na.rm = TRUE)
    if (!is.finite(scale) || scale <= 0) scale <- 1
    effect <- paste0("exp(THETA(", theta_number, ")*((", covariate, "-",
                     format(reference, digits = 15), ")/", format(scale, digits = 15), "))")
  }
  lines[[index]] <- paste0(parameter, "=(`", rhs, "`)*", effect)
  lines[[index]] <- gsub("`", "", lines[[index]], fixed = TRUE)
  theta <- model$THETAS
  new_row <- as.list(rep(NA, ncol(theta)))
  names(new_row) <- names(theta)
  new_row$THETA <- theta_number
  new_row$Value <- as.numeric(candidate$initial[[1L]] %||% 0)
  new_row$FIX <- FALSE
  if ("LOWER" %in% names(new_row)) new_row$LOWER <- -Inf
  if ("UPPER" %in% names(new_row)) new_row$UPPER <- Inf
  theta <- rbind(theta, as.data.frame(new_row, stringsAsFactors = FALSE))
  .nm_model_rebuild(model, list(
    PRED = paste(lines, collapse = "\n"), THETAS = theta,
    INPUT = unique(c(model$INPUT, covariate)),
    COVARIATES = unique(c(model$COVARIATES, covariate))
  ))
}

.nm_scm_model <- function(base_model, selected, candidates, data) {
  model <- base_model
  for (id in selected) {
    model <- .nm_scm_add_relationship(model, candidates[candidates$id == id, , drop = FALSE], data)
  }
  model
}

#' Stepwise covariate modelling
#'
#' Candidate relationships are tested with nested re-estimation and
#' likelihood-ratio thresholds. Continuous effects are exponential-linear,
#' `power` effects are normalized power functions, and categorical effects use
#' one indicator coefficient.
#'
#' @param fit Base `nm_fit`.
#' @param candidates Data frame with `parameter`, `covariate`, and optional
#'   `form`, `reference`, `category`, `scale`, and `initial` columns.
#' @param direction Forward selection, backward elimination, or both.
#' @param p_forward,p_backward Entry and retention significance levels.
#' @param max_steps Maximum accepted changes in each phase.
#' @param ... Controls passed to [nm_est()].
#' @return An `nm_scm` with the final model/fit and complete selection log.
#' @export
nm_scm <- function(fit, candidates, direction = c("both", "forward", "backward"),
                   p_forward = 0.05, p_backward = 0.01, max_steps = 20L, ...) {
  if (!inherits(fit, "nm_fit")) .nm_stop("`fit` must be an nm_fit.")
  direction <- match.arg(direction)
  candidates <- .nm_scm_candidates(candidates, fit$data)
  max_steps <- as.integer(max_steps)
  if (!is.finite(p_forward) || p_forward <= 0 || p_forward >= 1 ||
      !is.finite(p_backward) || p_backward <= 0 || p_backward >= 1 ||
      is.na(max_steps) || max_steps < 1L) {
    .nm_stop("SCM probabilities must lie between zero and one and `max_steps` must be positive.")
  }
  base <- fit$model
  base$THETAS$Value <- fit$theta
  base$SIGMAS$Value <- fit$sigma
  base$OMEGAS$Value <- fit$omega
  selected <- if (direction == "backward") candidates$id else character()
  current <- if (length(selected)) {
    nm_est(.nm_scm_model(base, selected, candidates, fit$data), fit$data,
           method = fit$method, ...)
  } else fit
  log <- list()
  step <- 0L
  if (direction %in% c("forward", "both")) {
    repeat {
      remaining <- setdiff(candidates$id, selected)
      if (!length(remaining) || step >= max_steps) break
      trials <- lapply(remaining, function(id) tryCatch(
        nm_est(.nm_scm_model(base, c(selected, id), candidates, fit$data), fit$data,
               method = fit$method, ...), error = identity
      ))
      drops <- vapply(trials, function(trial) if (inherits(trial, "nm_fit"))
        current$objective - trial$objective else -Inf, numeric(1))
      best <- which.max(drops)
      threshold <- stats::qchisq(1 - p_forward, df = 1L)
      accepted <- is.finite(drops[[best]]) && drops[[best]] >= threshold
      log[[length(log) + 1L]] <- data.frame(
        phase = "forward", step = step + 1L, candidate = remaining[[best]],
        objective = if (inherits(trials[[best]], "nm_fit")) trials[[best]]$objective else NA_real_,
        change = drops[[best]], threshold = threshold, accepted = accepted,
        error = if (inherits(trials[[best]], "error")) conditionMessage(trials[[best]]) else "",
        stringsAsFactors = FALSE
      )
      if (!accepted) break
      selected <- c(selected, remaining[[best]])
      current <- trials[[best]]
      step <- step + 1L
    }
  }
  if (direction %in% c("backward", "both") && length(selected)) {
    step <- 0L
    repeat {
      if (!length(selected) || step >= max_steps) break
      trials <- lapply(selected, function(id) tryCatch(
        nm_est(.nm_scm_model(base, setdiff(selected, id), candidates, fit$data), fit$data,
               method = fit$method, ...), error = identity
      ))
      increases <- vapply(trials, function(trial) if (inherits(trial, "nm_fit"))
        trial$objective - current$objective else Inf, numeric(1))
      best <- which.min(increases)
      threshold <- stats::qchisq(1 - p_backward, df = 1L)
      accepted <- is.finite(increases[[best]]) && increases[[best]] < threshold
      log[[length(log) + 1L]] <- data.frame(
        phase = "backward", step = step + 1L, candidate = selected[[best]],
        objective = if (inherits(trials[[best]], "nm_fit")) trials[[best]]$objective else NA_real_,
        change = increases[[best]], threshold = threshold, accepted = accepted,
        error = if (inherits(trials[[best]], "error")) conditionMessage(trials[[best]]) else "",
        stringsAsFactors = FALSE
      )
      if (!accepted) break
      selected <- setdiff(selected, selected[[best]])
      current <- trials[[best]]
      step <- step + 1L
    }
  }
  structure(list(
    selected = candidates[candidates$id %in% selected, , drop = FALSE],
    steps = if (length(log)) do.call(rbind, log) else data.frame(),
    final_model = current$model, final_fit = current, base_objective = fit$objective,
    final_objective = current$objective, direction = direction,
    p_forward = p_forward, p_backward = p_backward
  ), class = "nm_scm")
}
