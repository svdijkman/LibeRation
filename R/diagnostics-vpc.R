# VPC, NPC, NPDE, outcome-specific predictive checks, and bootstrap.
# Split from diagnostics.R as a behaviour-preserving source move.

#' Visual predictive check summaries
#'
#' @param fit An `nm_fit`.
#' @param nsim Number of stochastic replicates.
#' @param breaks Optional time-bin breaks. The default uses quantile bins.
#' @param probs Observation quantiles.
#' @param level Simulation interval for each quantile.
#' @param seed RNG seed.
#' @param pc_correct Apply the legacy prediction correction `DV * PRED / IPRED`
#'   to observed and simulated values before bin summaries.
#' @param stratify Optional dataset column used to create additional
#'   stratum-specific VPC summaries. The unstratified VPC is always retained.
#' @return An `nm_vpc` list containing observed and simulated summaries.
#' @export
nm_vpc <- function(fit, nsim = 200L, breaks = NULL,
                   probs = c(0.05, 0.5, 0.95), level = 0.9,
                   seed = 20260713L, pc_correct = FALSE, stratify = NULL) {
  if (!inherits(fit, "nm_fit")) .nm_stop("`fit` must be an nm_fit.")
  pc_correct <- isTRUE(pc_correct)
  stratify <- as.character(stratify %||% "")
  if (length(stratify) != 1L || is.na(stratify)) {
    .nm_stop("`stratify` must be one dataset column name or `NULL`.")
  }
  if (!nzchar(stratify)) stratify <- NULL
  if (!is.null(stratify) && !stratify %in% names(fit$data)) {
    .nm_stop("VPC stratification column `", stratify, "` is not present in the estimation data.")
  }
  source_data <- if (pc_correct) nm_gof(fit) else as.data.frame(fit$data)
  observed <- source_data[source_data$EVID == 0L & source_data$MDV == 0L &
                            is.finite(source_data$DV), , drop = FALSE]
  if (pc_correct) {
    valid <- is.finite(observed$PRED) & is.finite(observed$IPRED) &
      abs(observed$IPRED) > sqrt(.Machine$double.eps)
    observed$DV[valid] <- observed$DV[valid] * observed$PRED[valid] / observed$IPRED[valid]
    observed$DV[!valid] <- NA_real_
  }
  if (!nrow(observed)) .nm_stop("VPC requires observed DV records.")
  if (is.null(breaks)) {
    breaks <- unique(stats::quantile(observed$TIME, probs = seq(0, 1, length.out = 7), na.rm = TRUE))
    if (length(breaks) < 2L) breaks <- range(observed$TIME) + c(-0.5, 0.5)
  }
  simulated <- nm_simulate(
    fit$model, fit$data, theta = fit$theta, sigma = fit$sigma,
    omega = fit$omega, nsim = nsim, random_effects = TRUE,
    residual = TRUE, sample_mixture = TRUE, seed = seed
  )
  if (pc_correct) {
    population <- predict(fit, type = "population")$IPRED
    reference <- rep(population, times = as.integer(nsim))
    valid <- is.finite(reference) & is.finite(simulated$IPRED) &
      abs(simulated$IPRED) > sqrt(.Machine$double.eps)
    simulated$DV[valid] <- simulated$DV[valid] * reference[valid] / simulated$IPRED[valid]
    simulated$DV[!valid] <- NA_real_
  }
  simulated <- simulated[simulated$EVID == 0L & simulated$MDV == 0L, , drop = FALSE]
  qnames <- paste0("Q", formatC(100 * probs, format = "fg"))
  summarize_set <- function(observed_set, simulated_set) {
    observed_set$BIN <- cut(observed_set$TIME, breaks = breaks, include.lowest = TRUE)
    simulated_set$BIN <- cut(simulated_set$TIME, breaks = breaks, include.lowest = TRUE)
    summarize <- function(frame) {
      if (!nrow(frame)) return(NULL)
      values <- stats::quantile(frame$DV, probs = probs, na.rm = TRUE, names = FALSE)
      data.frame(
        BIN = as.character(frame$BIN[[1L]]),
        TIME = stats::median(frame$TIME, na.rm = TRUE), N = nrow(frame),
        stats::setNames(as.list(values), qnames), check.names = FALSE
      )
    }
    observed_summary <- do.call(
      rbind, lapply(split(observed_set, observed_set$BIN, drop = TRUE), summarize)
    )
    if (!"SIM" %in% names(simulated_set)) simulated_set$SIM <- 1L
    per_sim <- do.call(rbind, lapply(split(
      simulated_set, list(simulated_set$SIM, simulated_set$BIN), drop = TRUE
    ), function(frame) {
      if (!nrow(frame)) return(NULL)
      cbind(SIM = frame$SIM[[1L]], summarize(frame))
    }))
    alpha <- (1 - level) / 2
    intervals <- do.call(rbind, lapply(split(per_sim, per_sim$BIN), function(frame) {
      values <- lapply(qnames, function(name) {
        stats::quantile(
          frame[[name]], c(alpha, 0.5, 1 - alpha), na.rm = TRUE, names = FALSE
        )
      })
      row <- data.frame(
        BIN = frame$BIN[[1L]], TIME = stats::median(frame$TIME, na.rm = TRUE)
      )
      for (i in seq_along(qnames)) {
        row[[paste0(qnames[[i]], "_lo")]] <- values[[i]][[1L]]
        row[[paste0(qnames[[i]], "_median")]] <- values[[i]][[2L]]
        row[[paste0(qnames[[i]], "_hi")]] <- values[[i]][[3L]]
      }
      row
    }))
    list(
      observed = observed_summary,
      simulated = intervals,
      points = observed_set[
        is.finite(observed_set$TIME) & is.finite(observed_set$DV),
        c("TIME", "DV"), drop = FALSE
      ],
      per_simulation = per_sim
    )
  }
  total <- summarize_set(observed, simulated)
  stratified <- list()
  if (!is.null(stratify)) {
    strata <- unique(as.character(observed[[stratify]]))
    strata <- strata[!is.na(strata) & nzchar(strata)]
    simulated_strata <- as.character(simulated[[stratify]])
    observed_strata <- as.character(observed[[stratify]])
    stratified <- unname(lapply(strata, function(stratum) {
      observed_set <- observed[!is.na(observed_strata) & observed_strata == stratum, , drop = FALSE]
      simulated_set <- simulated[!is.na(simulated_strata) & simulated_strata == stratum, , drop = FALSE]
      if (!nrow(observed_set) || !nrow(simulated_set)) return(NULL)
      c(list(level = stratum), summarize_set(observed_set, simulated_set))
    }))
    stratified <- Filter(Negate(is.null), stratified)
  }
  structure(c(
    total,
    list(
      breaks = breaks, probs = probs, level = level, nsim = nsim, seed = seed,
      pc_correct = pc_correct, stratify = stratify, stratified = stratified
    )
  ), class = "nm_vpc")
}

.nm_predictive_simulation_matrix <- function(fit, nsim, seed) {
  if (!inherits(fit, "nm_fit")) .nm_stop("`fit` must be an nm_fit.")
  nsim <- as.integer(nsim)
  if (length(nsim) != 1L || is.na(nsim) || nsim < 20L) {
    .nm_stop("`nsim` must be an integer of at least 20.")
  }
  simulated <- nm_simulate(
    fit$model, fit$data, theta = fit$theta, sigma = fit$sigma,
    omega = fit$omega, nsim = nsim, random_effects = TRUE,
    residual = TRUE, sample_mixture = TRUE, seed = seed
  )
  records <- nrow(fit$data)
  if (nrow(simulated) != records * nsim) {
    .nm_stop("Predictive simulation rows do not align with the fitted dataset.")
  }
  matrix(as.numeric(simulated$DV), nrow = records, ncol = nsim)
}

#' Numerical predictive check
#'
#' Computes the empirical predictive percentile of every uncensored
#' observation under repeated simulations from the fitted population model.
#'
#' @param fit An `nm_fit`.
#' @param nsim Number of predictive simulations; at least 20.
#' @param seed Reproducible RNG seed.
#' @return An `nm_npc` object with record-level percentiles and tail flags.
#' @export
nm_npc <- function(fit, nsim = 200L, seed = 20260713L) {
  simulations <- .nm_predictive_simulation_matrix(fit, nsim, seed)
  data <- as.data.frame(fit$data)
  observed <- data$EVID == 0L & data$MDV == 0L & is.finite(data$DV)
  if ("CENS" %in% names(data)) observed <- observed & data$CENS != 1L
  if ("BLQ" %in% names(data)) observed <- observed & data$BLQ != 1L
  rows <- which(observed)
  if (!length(rows)) .nm_stop("NPC requires uncensored observed DV records.")
  count <- vapply(rows, function(row) sum(simulations[row, ] <= data$DV[[row]], na.rm = TRUE), numeric(1))
  percentile <- (count + 0.5) / (ncol(simulations) + 1)
  table <- data.frame(
    ROW = rows, ID = data$ID[rows], TIME = data$TIME[rows], DV = data$DV[rows],
    PERCENTILE = percentile,
    TAIL_PROBABILITY = 2 * pmin(percentile, 1 - percentile),
    OUTSIDE_90 = percentile < 0.05 | percentile > 0.95,
    stringsAsFactors = FALSE
  )
  structure(list(
    table = table, nsim = ncol(simulations), seed = as.integer(seed),
    outside_90 = mean(table$OUTSIDE_90),
    histogram = graphics::hist(
      table$PERCENTILE, breaks = seq(0, 1, length.out = 11), plot = FALSE
    )
  ), class = "nm_npc")
}

#' Normalized prediction distribution errors
#'
#' Within each subject, simulations are centered and decorrelated with their
#' empirical predictive covariance before record-wise predictive ranks are
#' transformed to standard-normal NPDE values.
#'
#' @param fit An `nm_fit`.
#' @param nsim Number of predictive simulations; at least 20.
#' @param seed Reproducible RNG seed.
#' @param ridge Relative covariance ridge used for numerically singular blocks.
#' @return An `nm_npde` object with record-level NPDE and summary moments.
#' @export
nm_npde <- function(fit, nsim = 200L, seed = 20260713L, ridge = 1e-8) {
  simulations <- .nm_predictive_simulation_matrix(fit, nsim, seed)
  data <- as.data.frame(fit$data)
  observed <- data$EVID == 0L & data$MDV == 0L & is.finite(data$DV)
  if ("CENS" %in% names(data)) observed <- observed & data$CENS != 1L
  if ("BLQ" %in% names(data)) observed <- observed & data$BLQ != 1L
  rows <- which(observed)
  if (!length(rows)) .nm_stop("NPDE requires uncensored observed DV records.")
  ridge <- as.numeric(ridge)
  if (length(ridge) != 1L || !is.finite(ridge) || ridge <= 0) {
    .nm_stop("`ridge` must be a positive finite scalar.")
  }
  percentile <- rep(NA_real_, length(rows))
  groups <- split(seq_along(rows), data$.ID_INDEX[rows])
  for (indices in groups) {
    record_rows <- rows[indices]
    block <- simulations[record_rows, , drop = FALSE]
    center <- rowMeans(block)
    covariance <- if (nrow(block) == 1L) {
      matrix(stats::var(drop(block)), 1L, 1L)
    } else stats::cov(t(block))
    scale <- max(mean(diag(covariance)), 1)
    covariance <- covariance + diag(ridge * scale, nrow(covariance))
    root <- tryCatch(chol(covariance), error = function(e) NULL)
    if (is.null(root)) root <- chol(covariance + diag(sqrt(ridge) * scale, nrow(covariance)))
    observed_white <- forwardsolve(t(root), data$DV[record_rows] - center)
    simulated_white <- forwardsolve(t(root), sweep(block, 1L, center, "-"))
    percentile[indices] <- vapply(seq_along(indices), function(position) {
      (sum(simulated_white[position, ] <= observed_white[[position]]) + 0.5) /
        (ncol(simulated_white) + 1)
    }, numeric(1))
  }
  percentile <- pmin(pmax(percentile, 1e-12), 1 - 1e-12)
  npde <- stats::qnorm(percentile)
  centered <- npde - mean(npde)
  spread <- stats::sd(npde)
  table <- data.frame(
    ROW = rows, ID = data$ID[rows], TIME = data$TIME[rows], DV = data$DV[rows],
    PERCENTILE = percentile, NPDE = npde, stringsAsFactors = FALSE
  )
  structure(list(
    table = table, nsim = ncol(simulations), seed = as.integer(seed), ridge = ridge,
    summary = c(
      mean = mean(npde), sd = spread,
      skewness = if (is.finite(spread) && spread > 0) mean(centered^3) / spread^3 else NA_real_,
      kurtosis = if (is.finite(spread) && spread > 0) mean(centered^4) / spread^4 - 3 else NA_real_
    )
  ), class = "nm_npde")
}

.nm_vpc_breaks <- function(time, breaks = NULL) {
  if (!is.null(breaks)) return(as.numeric(breaks))
  breaks <- unique(stats::quantile(time, probs = seq(0, 1, length.out = 7), na.rm = TRUE))
  if (length(breaks) < 2L) breaks <- range(time, na.rm = TRUE) + c(-0.5, 0.5)
  breaks
}

#' Categorical visual predictive check
#'
#' First-class Bernoulli, categorical, ordinal, and Markov outcomes use their
#' declared probability vectors and can contain any number of categories. A
#' legacy binary user likelihood remains supported when `F` is the probability
#' of the non-reference category.
#'
#' @param fit An `nm_fit`.
#' @param outcome Binary observed outcome column.
#' @param nsim Number of simulations.
#' @param breaks Optional time bins.
#' @param level Simulation interval.
#' @param seed RNG seed.
#' @return An `nm_vpc_categorical`.
#' @export
nm_vpc_categorical <- function(fit, outcome = "DV", nsim = 200L, breaks = NULL,
                               level = 0.9, seed = 20260713L) {
  if (!inherits(fit, "nm_fit")) .nm_stop("`fit` must be an nm_fit.")
  outcome <- as.character(outcome)
  if (length(outcome) != 1L || !outcome %in% names(fit$data)) .nm_stop("Unknown categorical outcome column.")
  nsim <- as.integer(nsim)
  if (is.na(nsim) || nsim < 20L) .nm_stop("`nsim` must be at least 20.")
  if (!is.finite(level) || level <= 0 || level >= 1) .nm_stop("`level` must lie between zero and one.")
  data <- as.data.frame(fit$data)
  observed_rows <- data$EVID == 0L & data$MDV == 0L & !is.na(data[[outcome]])
  observed <- data[observed_rows, , drop = FALSE]
  categories <- sort(unique(observed[[outcome]]))
  declared <- NULL
  if (!is.null(fit$model$OUTCOMES)) {
    candidates <- Filter(function(value) {
      value$family %in% c("bernoulli", "categorical", "ordinal", "markov",
                          "continuous_time_markov")
    }, fit$model$OUTCOMES)
    if (length(candidates) == 1L) declared <- candidates[[1L]]
    if (length(candidates) > 1L && "DVID" %in% names(observed)) {
      dvid <- unique(observed$DVID)
      if (length(dvid) == 1L) {
        matching <- Filter(function(value) identical(as.numeric(value$dvid), as.numeric(dvid)), candidates)
        if (length(matching) == 1L) declared <- matching[[1L]]
      }
    }
  }
  if (!is.null(declared)) categories <- declared$categories %||% c(0, 1)
  if (length(categories) < 2L) .nm_stop("Categorical VPC requires at least two categories.")
  if (is.null(declared) && length(categories) != 2L) {
    .nm_stop("A multicategory VPC requires a first-class categorical OUTCOMES declaration.")
  }
  breaks <- .nm_vpc_breaks(observed$TIME, breaks)
  observed$BIN <- cut(observed$TIME, breaks, include.lowest = TRUE)
  observed_summary <- do.call(rbind, lapply(
    split(observed, observed$BIN, drop = TRUE), function(frame) do.call(rbind, lapply(
      categories, function(category) data.frame(
        BIN = as.character(frame$BIN[[1L]]), TIME = stats::median(frame$TIME),
        N = nrow(frame), CATEGORY = as.character(category),
        PROPORTION = mean(frame[[outcome]] == category), stringsAsFactors = FALSE
      )
    ))
  ))
  simulated <- nm_simulate(
    fit$model, fit$data, theta = fit$theta, sigma = fit$sigma, omega = fit$omega,
    nsim = nsim, random_effects = TRUE, residual = !is.null(declared),
    sample_mixture = TRUE, seed = seed
  )
  simulated <- simulated[rep(observed_rows, times = nsim), , drop = FALSE]
  if (is.null(declared)) {
    probability <- as.numeric(simulated$IPRED)
    if (mean(!is.finite(probability) | probability < -1e-6 | probability > 1 + 1e-6) > 0.01) {
      .nm_stop("Categorical VPC requires F/IPRED to represent a probability in [0, 1].")
    }
    probability <- pmin(pmax(probability, 0), 1)
    set.seed(as.integer(seed) + 1L)
    simulated$DV <- ifelse(stats::runif(nrow(simulated)) <= probability,
                           categories[[2L]], categories[[1L]])
  }
  simulated$CATEGORY <- simulated[[outcome]]
  simulated$BIN <- cut(simulated$TIME, breaks, include.lowest = TRUE)
  per_simulation <- do.call(rbind, lapply(
    split(simulated, list(simulated$SIM, simulated$BIN), drop = TRUE),
    function(frame) do.call(rbind, lapply(categories, function(category) data.frame(
      SIM = frame$SIM[[1L]], BIN = as.character(frame$BIN[[1L]]),
      TIME = stats::median(frame$TIME), CATEGORY = as.character(category),
      PROPORTION = mean(frame$CATEGORY == category), stringsAsFactors = FALSE
    )))
  ))
  alpha <- (1 - level) / 2
  intervals <- do.call(rbind, lapply(split(
    per_simulation, list(per_simulation$BIN, per_simulation$CATEGORY), drop = TRUE
  ), function(frame) {
    interval <- stats::quantile(frame$PROPORTION, c(alpha, 0.5, 1 - alpha),
                                names = FALSE, na.rm = TRUE)
    data.frame(
      BIN = frame$BIN[[1L]], TIME = stats::median(frame$TIME),
      CATEGORY = frame$CATEGORY[[1L]], lower = interval[[1L]],
      median = interval[[2L]], upper = interval[[3L]], stringsAsFactors = FALSE
    )
  }))
  structure(list(
    observed = observed_summary, simulated = intervals,
    per_simulation = per_simulation, categories = categories, outcome = outcome,
    breaks = breaks, level = level, nsim = nsim, seed = seed
  ), class = "nm_vpc_categorical")
}

#' Count visual predictive check
#'
#' Summarizes the mean, variance, zero fraction, median, and upper count
#' quantile in time bins for first-class Poisson, negative-binomial, binomial,
#' ZIP, and hurdle models.
#'
#' @param fit An `nm_fit`.
#' @param outcome Count column, normally `DV`.
#' @param dvid Optional endpoint `DVID` in a joint model.
#' @param nsim Number of simulations.
#' @param breaks Optional time bins.
#' @param level Simulation interval.
#' @param seed RNG seed.
#' @return An `nm_vpc_count`.
#' @export
nm_vpc_count <- function(fit, outcome = "DV", dvid = NULL, nsim = 200L,
                         breaks = NULL, level = 0.9, seed = 20260713L) {
  if (!inherits(fit, "nm_fit")) .nm_stop("`fit` must be an nm_fit.")
  if (is.null(fit$model$OUTCOMES)) .nm_stop("Count VPC requires first-class OUTCOMES.")
  if (!is.null(dvid) && (length(dvid) != 1L || !is.finite(dvid))) dvid <- NULL
  count_families <- c("poisson", "negative_binomial", "binomial",
                      "zero_inflated_poisson", "hurdle_poisson")
  candidates <- Filter(function(value) value$family %in% count_families,
                       fit$model$OUTCOMES)
  if (!is.null(dvid)) candidates <- Filter(function(value) {
    identical(as.numeric(value$dvid), as.numeric(dvid))
  }, candidates)
  if (length(candidates) != 1L) {
    .nm_stop("Select a unique count endpoint with `dvid`.")
  }
  endpoint <- candidates[[1L]]
  nsim <- as.integer(nsim)
  if (is.na(nsim) || nsim < 20L) .nm_stop("`nsim` must be at least 20.")
  if (!is.finite(level) || level <= 0 || level >= 1) .nm_stop("`level` must lie between zero and one.")
  data <- as.data.frame(fit$data)
  rows <- data$EVID == 0L & data$MDV == 0L & is.finite(data[[outcome]])
  if (!is.null(endpoint$dvid)) rows <- rows & data$DVID == endpoint$dvid
  observed <- data[rows, , drop = FALSE]
  if (!nrow(observed)) .nm_stop("No observed count records were found.")
  breaks <- .nm_vpc_breaks(observed$TIME, breaks)
  summarize <- function(frame) data.frame(
    BIN = as.character(frame$BIN[[1L]]), TIME = stats::median(frame$TIME),
    N = nrow(frame), MEAN = mean(frame[[outcome]]),
    VARIANCE = if (nrow(frame) > 1L) stats::var(frame[[outcome]]) else 0,
    ZERO = mean(frame[[outcome]] == 0),
    Q50 = unname(stats::quantile(frame[[outcome]], 0.5, type = 1)),
    Q90 = unname(stats::quantile(frame[[outcome]], 0.9, type = 1)),
    stringsAsFactors = FALSE
  )
  observed$BIN <- cut(observed$TIME, breaks, include.lowest = TRUE)
  observed_summary <- do.call(rbind, lapply(split(observed, observed$BIN, drop = TRUE), summarize))
  simulated <- nm_simulate(
    fit$model, fit$data, theta = fit$theta, sigma = fit$sigma, omega = fit$omega,
    nsim = nsim, random_effects = TRUE, residual = TRUE,
    sample_mixture = TRUE, seed = seed
  )
  simulated <- simulated[rep(rows, times = nsim), , drop = FALSE]
  simulated$BIN <- cut(simulated$TIME, breaks, include.lowest = TRUE)
  per_simulation <- do.call(rbind, lapply(split(
    simulated, list(simulated$SIM, simulated$BIN), drop = TRUE
  ), function(frame) cbind(SIM = frame$SIM[[1L]], summarize(frame))))
  measures <- c("MEAN", "VARIANCE", "ZERO", "Q50", "Q90")
  alpha <- (1 - level) / 2
  intervals <- do.call(rbind, lapply(split(per_simulation, per_simulation$BIN), function(frame) {
    row <- data.frame(BIN = frame$BIN[[1L]], TIME = stats::median(frame$TIME))
    for (measure in measures) {
      interval <- stats::quantile(frame[[measure]], c(alpha, 0.5, 1 - alpha),
                                  names = FALSE, na.rm = TRUE)
      row[[paste0(measure, "_lower")]] <- interval[[1L]]
      row[[paste0(measure, "_median")]] <- interval[[2L]]
      row[[paste0(measure, "_upper")]] <- interval[[3L]]
    }
    row
  }))
  structure(list(
    observed = observed_summary, simulated = intervals,
    per_simulation = per_simulation, outcome = outcome, dvid = endpoint$dvid,
    family = endpoint$family, breaks = breaks, level = level,
    nsim = nsim, seed = seed
  ), class = "nm_vpc_count")
}

.nm_km_at <- function(time, event, grid) {
  order <- order(time)
  time <- time[order]
  event <- event[order]
  unique_event <- sort(unique(time[event == 1L]))
  survival <- 1
  event_curve <- numeric(length(unique_event))
  for (index in seq_along(unique_event)) {
    current <- unique_event[[index]]
    at_risk <- sum(time >= current)
    events <- sum(time == current & event == 1L)
    if (at_risk > 0L) survival <- survival * (1 - events / at_risk)
    event_curve[[index]] <- survival
  }
  if (!length(unique_event)) return(rep(1, length(grid)))
  position <- findInterval(grid, unique_event)
  ifelse(position == 0L, 1, event_curve[pmax(position, 1L)])
}

#' Time-to-event visual predictive check
#'
#' `F`/`IPRED` is interpreted as a non-negative instantaneous hazard on each
#' subject's observation grid. The first event is simulated from the integrated
#' piecewise-constant hazard and Kaplan-Meier curves are summarized across
#' simulations.
#'
#' @param fit An `nm_fit`.
#' @param event Binary event-indicator column.
#' @param nsim Number of simulations.
#' @param level Simulation interval.
#' @param seed RNG seed.
#' @return An `nm_vpc_tte`.
#' @export
nm_vpc_tte <- function(fit, event = "DV", nsim = 200L, level = 0.9,
                       seed = 20260713L) {
  if (!inherits(fit, "nm_fit")) .nm_stop("`fit` must be an nm_fit.")
  event <- as.character(event)
  if (length(event) != 1L || !event %in% names(fit$data)) .nm_stop("Unknown event indicator column.")
  nsim <- as.integer(nsim)
  if (is.na(nsim) || nsim < 20L) .nm_stop("`nsim` must be at least 20.")
  if (!is.finite(level) || level <= 0 || level >= 1) .nm_stop("`level` must lie between zero and one.")
  data <- as.data.frame(fit$data)
  rows <- data$EVID == 0L & data$MDV == 0L & !is.na(data[[event]])
  observed <- data[rows, , drop = FALSE]
  if (!all(observed[[event]] %in% c(0, 1))) .nm_stop("The TTE event column must contain zero/one indicators.")
  subjects <- split(observed, observed$.ID_INDEX)
  subject_records <- do.call(rbind, lapply(subjects, function(frame) {
    event_rows <- which(frame[[event]] == 1L)
    data.frame(
      ID = frame$ID[[1L]],
      TIME = if (length(event_rows)) frame$TIME[[event_rows[[1L]]]] else max(frame$TIME),
      EVENT = as.integer(length(event_rows) > 0L)
    )
  }))
  grid <- sort(unique(observed$TIME))
  observed_curve <- data.frame(
    TIME = grid, SURVIVAL = .nm_km_at(subject_records$TIME, subject_records$EVENT, grid)
  )
  declared <- NULL
  if (!is.null(fit$model$OUTCOMES)) {
    candidates <- Filter(function(value) {
      value$family %in% c("tte", "competing_risks")
    }, fit$model$OUTCOMES)
    if (length(candidates) == 1L) declared <- candidates[[1L]]
  }
  simulated <- nm_simulate(
    fit$model, fit$data, theta = fit$theta, sigma = fit$sigma, omega = fit$omega,
    nsim = nsim, random_effects = TRUE, residual = !is.null(declared),
    sample_mixture = TRUE, seed = seed
  )
  simulated <- simulated[rep(rows, times = nsim), , drop = FALSE]
  if (is.null(declared) && any(!is.finite(simulated$IPRED) | simulated$IPRED < 0)) {
    .nm_stop("TTE VPC requires F/IPRED to be a finite non-negative hazard.")
  }
  set.seed(as.integer(seed) + 1L)
  curves <- matrix(1, nrow = length(grid), ncol = nsim)
  for (simulation in seq_len(nsim)) {
    sample <- simulated[simulated$SIM == simulation, , drop = FALSE]
    sample_subjects <- split(sample, sample$.ID_INDEX)
    records <- do.call(rbind, lapply(sample_subjects, function(frame) {
      frame <- frame[order(frame$TIME), , drop = FALSE]
      event_rows <- if (!is.null(declared)) {
        if (declared$family == "competing_risks") which(frame[[event]] != 0) else
          which(frame[[event]] == declared$event)
      } else {
        delta <- c(0, pmax(diff(frame$TIME), 0))
        probability <- 1 - exp(-pmax(frame$IPRED, 0) * delta)
        which(stats::runif(nrow(frame)) <= probability)
      }
      data.frame(
        TIME = if (length(event_rows)) frame$TIME[[event_rows[[1L]]]] else max(frame$TIME),
        EVENT = as.integer(length(event_rows) > 0L)
      )
    }))
    curves[, simulation] <- .nm_km_at(records$TIME, records$EVENT, grid)
  }
  alpha <- (1 - level) / 2
  intervals <- t(apply(curves, 1L, stats::quantile,
                       probs = c(alpha, 0.5, 1 - alpha), names = FALSE))
  simulated_curve <- data.frame(
    TIME = grid, lower = intervals[, 1L], median = intervals[, 2L],
    upper = intervals[, 3L]
  )
  structure(list(
    observed = observed_curve, simulated = simulated_curve,
    per_simulation = curves, event = event, level = level,
    nsim = nsim, seed = seed
  ), class = "nm_vpc_tte")
}

.nm_competing_curve <- function(records, causes, grid) {
  event_times <- sort(unique(records$TIME[records$CAUSE != 0]))
  survival <- 1
  cumulative <- stats::setNames(numeric(length(causes)), as.character(causes))
  history <- matrix(0, nrow = length(event_times), ncol = length(causes),
                    dimnames = list(NULL, as.character(causes)))
  for (index in seq_along(event_times)) {
    time <- event_times[[index]]
    at_risk <- sum(records$TIME >= time)
    if (at_risk > 0L) {
      events <- vapply(causes, function(cause) {
        sum(records$TIME == time & records$CAUSE == cause)
      }, numeric(1))
      cumulative <- cumulative + survival * events / at_risk
      survival <- survival * (1 - sum(events) / at_risk)
    }
    history[index, ] <- cumulative
  }
  if (!length(event_times)) history <- matrix(0, nrow = 1L, ncol = length(causes),
                                             dimnames = list(NULL, as.character(causes)))
  do.call(rbind, lapply(seq_along(causes), function(column) {
    position <- findInterval(grid, event_times)
    value <- if (!length(event_times)) rep(0, length(grid)) else
      ifelse(position == 0L, 0, history[pmax(position, 1L), column])
    data.frame(TIME = grid, CAUSE = as.character(causes[[column]]), CIF = value,
               stringsAsFactors = FALSE)
  }))
}

#' Competing-risk visual predictive check
#'
#' Uses Aalen-Johansen cumulative incidence curves for every declared cause and
#' compares them with simulation intervals from a first-class
#' `competing_risks` outcome.
#'
#' @param fit An `nm_fit`.
#' @param event Cause-code column (`0` denotes no event/censoring).
#' @param dvid Optional joint-endpoint `DVID`.
#' @param nsim Number of simulations.
#' @param level Simulation interval.
#' @param seed RNG seed.
#' @return An `nm_vpc_competing`.
#' @export
nm_vpc_competing <- function(fit, event = "DV", dvid = NULL, nsim = 200L,
                             level = 0.9, seed = 20260713L) {
  if (!inherits(fit, "nm_fit")) .nm_stop("`fit` must be an nm_fit.")
  if (!is.null(dvid) && (length(dvid) != 1L || !is.finite(dvid))) dvid <- NULL
  candidates <- Filter(function(value) value$family == "competing_risks",
                       fit$model$OUTCOMES %||% list())
  if (!is.null(dvid)) candidates <- Filter(function(value) {
    identical(as.numeric(value$dvid), as.numeric(dvid))
  }, candidates)
  if (length(candidates) != 1L) .nm_stop("Select a unique competing-risk endpoint with `dvid`.")
  endpoint <- candidates[[1L]]
  nsim <- as.integer(nsim)
  if (is.na(nsim) || nsim < 20L) .nm_stop("`nsim` must be at least 20.")
  data <- as.data.frame(fit$data)
  rows <- data$EVID == 0L & data$MDV == 0L & is.finite(data[[event]])
  if (!is.null(endpoint$dvid)) rows <- rows & data$DVID == endpoint$dvid
  observed <- data[rows, , drop = FALSE]
  records <- function(frame) do.call(rbind, lapply(split(frame, frame$.ID_INDEX), function(subject) {
    subject <- subject[order(subject$TIME), , drop = FALSE]
    event_row <- which(subject[[event]] != 0)[1L]
    if (is.na(event_row)) event_row <- integer()
    data.frame(
      TIME = if (length(event_row)) subject$TIME[[event_row]] else max(subject$TIME),
      CAUSE = if (length(event_row)) subject[[event]][[event_row]] else 0
    )
  }))
  grid <- sort(unique(observed$TIME))
  observed_curve <- .nm_competing_curve(records(observed), endpoint$categories, grid)
  simulated <- nm_simulate(
    fit$model, fit$data, theta = fit$theta, sigma = fit$sigma, omega = fit$omega,
    nsim = nsim, random_effects = TRUE, residual = TRUE,
    sample_mixture = TRUE, seed = seed
  )
  simulated <- simulated[rep(rows, times = nsim), , drop = FALSE]
  curves <- do.call(rbind, lapply(seq_len(nsim), function(index) {
    curve <- .nm_competing_curve(
      records(simulated[simulated$SIM == index, , drop = FALSE]),
      endpoint$categories, grid
    )
    curve$SIM <- index
    curve
  }))
  alpha <- (1 - level) / 2
  intervals <- do.call(rbind, lapply(split(
    curves, list(curves$CAUSE, curves$TIME), drop = TRUE
  ), function(frame) {
    interval <- stats::quantile(frame$CIF, c(alpha, 0.5, 1 - alpha), names = FALSE)
    data.frame(
      TIME = frame$TIME[[1L]], CAUSE = frame$CAUSE[[1L]],
      lower = interval[[1L]], median = interval[[2L]], upper = interval[[3L]],
      stringsAsFactors = FALSE
    )
  }))
  intervals <- intervals[order(as.numeric(intervals$CAUSE), intervals$TIME), , drop = FALSE]
  structure(list(
    observed = observed_curve, simulated = intervals, per_simulation = curves,
    event = event, dvid = endpoint$dvid, causes = endpoint$categories,
    level = level, nsim = nsim, seed = seed
  ), class = "nm_vpc_competing")
}

#' Recurrent-event visual predictive check
#'
#' Compares the observed mean cumulative event function with predictive
#' intervals from a first-class `recurrent_event` outcome.
#'
#' @param fit An `nm_fit`.
#' @param event Event-indicator column.
#' @param dvid Optional joint-endpoint `DVID`.
#' @param nsim Number of simulations.
#' @param level Simulation interval.
#' @param seed RNG seed.
#' @return An `nm_vpc_recurrent`.
#' @export
nm_vpc_recurrent <- function(fit, event = "DV", dvid = NULL, nsim = 200L,
                             level = 0.9, seed = 20260713L) {
  if (!inherits(fit, "nm_fit")) .nm_stop("`fit` must be an nm_fit.")
  if (!is.null(dvid) && (length(dvid) != 1L || !is.finite(dvid))) dvid <- NULL
  candidates <- Filter(function(value) value$family == "recurrent_event",
                       fit$model$OUTCOMES %||% list())
  if (!is.null(dvid)) candidates <- Filter(function(value) {
    identical(as.numeric(value$dvid), as.numeric(dvid))
  }, candidates)
  if (length(candidates) != 1L) .nm_stop("Select a unique recurrent-event endpoint with `dvid`.")
  endpoint <- candidates[[1L]]
  nsim <- as.integer(nsim)
  if (is.na(nsim) || nsim < 20L) .nm_stop("`nsim` must be at least 20.")
  data <- as.data.frame(fit$data)
  rows <- data$EVID == 0L & data$MDV == 0L & is.finite(data[[event]])
  if (!is.null(endpoint$dvid)) rows <- rows & data$DVID == endpoint$dvid
  observed <- data[rows, , drop = FALSE]
  grid <- sort(unique(observed$TIME))
  mean_cumulative <- function(frame) {
    subjects <- split(frame, frame$.ID_INDEX)
    curves <- vapply(subjects, function(subject) vapply(grid, function(time) {
      sum(subject[[event]][subject$TIME <= time] == endpoint$event)
    }, numeric(1)), numeric(length(grid)))
    if (is.null(dim(curves))) curves <- matrix(curves, ncol = 1L)
    rowMeans(curves)
  }
  observed_curve <- data.frame(TIME = grid, MEAN_CUMULATIVE = mean_cumulative(observed))
  simulated <- nm_simulate(
    fit$model, fit$data, theta = fit$theta, sigma = fit$sigma, omega = fit$omega,
    nsim = nsim, random_effects = TRUE, residual = TRUE,
    sample_mixture = TRUE, seed = seed
  )
  simulated <- simulated[rep(rows, times = nsim), , drop = FALSE]
  curves <- vapply(seq_len(nsim), function(index) {
    mean_cumulative(simulated[simulated$SIM == index, , drop = FALSE])
  }, numeric(length(grid)))
  alpha <- (1 - level) / 2
  intervals <- t(apply(curves, 1L, stats::quantile,
                       probs = c(alpha, 0.5, 1 - alpha), names = FALSE))
  simulated_curve <- data.frame(
    TIME = grid, lower = intervals[, 1L], median = intervals[, 2L], upper = intervals[, 3L]
  )
  structure(list(
    observed = observed_curve, simulated = simulated_curve,
    per_simulation = curves, event = event, dvid = endpoint$dvid,
    level = level, nsim = nsim, seed = seed
  ), class = "nm_vpc_recurrent")
}

#' Nonparametric and parametric bootstrap uncertainty
#'
#' @param fit An `nm_fit` used as the model and estimation template.
#' @param n Number of bootstrap fits.
#' @param seed RNG seed.
#' @param level Percentile confidence level.
#' @param type `nonparametric` resamples observed units; `parametric`
#'   simulates new outcomes and random effects from the fitted model.
#' @param unit Resampling unit for a nonparametric bootstrap: subjects or a
#'   user-supplied cluster column.
#' @param strata Optional column whose subject/cluster-level values define
#'   independent resampling strata.
#' @param cluster Cluster column required for `unit = "cluster"`.
#' @param ... Controls passed to [nm_est()].
#' @return Bootstrap parameter matrix, convergence flags, and failed-run errors.
#' @export
nm_bootstrap <- function(fit, n = 100L, seed = 20260713L, level = 0.95,
                         type = c("nonparametric", "parametric"),
                         unit = c("subject", "cluster"),
                         strata = NULL, cluster = NULL, ...) {
  if (!inherits(fit, "nm_fit")) .nm_stop("`fit` must be an nm_fit.")
  type <- match.arg(type)
  unit <- match.arg(unit)
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) .nm_stop("`n` must be a positive integer.")
  if (!is.finite(level) || level <= 0 || level >= 1) .nm_stop("`level` must lie between zero and one.")
  strata <- as.character(strata %||% "")
  cluster <- as.character(cluster %||% "")
  if (length(strata) != 1L || length(cluster) != 1L) {
    .nm_stop("`strata` and `cluster` must each be one column name or NULL.")
  }
  if (type == "parametric") {
    unit <- "subject"
    strata <- cluster <- ""
  }
  if (nzchar(strata) && !strata %in% names(fit$data)) {
    .nm_stop("Bootstrap stratum column `", strata, "` is absent from the fitted data.")
  }
  if (unit == "cluster" && (!nzchar(cluster) || !cluster %in% names(fit$data))) {
    .nm_stop("Cluster bootstrap requires a valid `cluster` column.")
  }
  strip_internal <- function(data) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    data[grep("^\\.", names(data), value = TRUE)] <- NULL
    data
  }
  resample_nonparametric <- function() {
    data <- as.data.frame(fit$data, stringsAsFactors = FALSE)
    key <- if (unit == "subject") data$ID else data[[cluster]]
    units <- unique(key)
    unit_strata <- if (nzchar(strata)) vapply(units, function(value) {
      observed <- unique(data[[strata]][key == value])
      observed <- observed[!is.na(observed)]
      if (length(observed) > 1L) {
        .nm_stop("Bootstrap strata must be constant within each resampling unit.")
      }
      if (length(observed)) as.character(observed[[1L]]) else "<missing>"
    }, character(1)) else rep("all", length(units))
    selected <- unlist(lapply(split(units, unit_strata), function(values) {
      sample(values, length(values), replace = TRUE)
    }), use.names = FALSE)
    pieces <- lapply(seq_along(selected), function(index) {
      block <- data[key == selected[[index]], , drop = FALSE]
      if (unit == "subject") {
        block$ID <- index
      } else {
        original_ids <- unique(block$ID)
        block$ID <- match(block$ID, original_ids) +
          sum(vapply(seq_len(index - 1L), function(previous) {
            length(unique(data$ID[key == selected[[previous]]]))
          }, integer(1)))
        block[[cluster]] <- index
      }
      strip_internal(block)
    })
    do.call(rbind, pieces)
  }
  fitted_model <- .nm_model_rebuild(fit$model, list(
    THETAS = transform(fit$model$THETAS, Value = fit$theta),
    OMEGAS = transform(fit$model$OMEGAS, Value = fit$omega),
    SIGMAS = transform(fit$model$SIGMAS, Value = fit$sigma)
  ))
  set.seed(seed)
  runs <- vector("list", n)
  errors <- character(n)
  convergence <- rep(NA_integer_, n)
  replicate_seeds <- sample.int(.Machine$integer.max, n)
  for (iteration in seq_len(n)) {
    set.seed(replicate_seeds[[iteration]])
    dataset <- if (type == "parametric") {
      strip_internal(nm_simulate(
        fitted_model, fit$data, theta = fit$theta, sigma = fit$sigma,
        omega = fit$omega, random_effects = TRUE, residual = TRUE,
        sample_mixture = TRUE, seed = replicate_seeds[[iteration]]
      ))
    } else resample_nonparametric()
    refit <- tryCatch(
      nm_est(fitted_model, dataset, method = fit$method, ...),
      error = identity
    )
    if (inherits(refit, "error")) {
      errors[[iteration]] <- conditionMessage(refit)
    } else {
      convergence[[iteration]] <- refit$convergence
      runs[[iteration]] <- c(refit$theta, refit$sigma, refit$omega)
    }
  }
  successful <- Filter(Negate(is.null), runs)
  estimates <- if (length(successful)) do.call(rbind, successful) else matrix(numeric(), 0L, 0L)
  if (ncol(estimates)) {
    colnames(estimates) <- .nm_parameter_names(fit$theta, fit$sigma, fit$omega)
  }
  native <- .nm_fit_native_parameters(fit)
  alpha <- (1 - level) / 2
  summary <- if (nrow(estimates)) do.call(rbind, lapply(seq_along(native), function(index) {
    values <- estimates[, index]
    interval <- stats::quantile(values, c(alpha, 1 - alpha), na.rm = TRUE, names = FALSE)
    data.frame(
      parameter = names(native)[[index]], estimate = native[[index]],
      bootstrap_mean = mean(values, na.rm = TRUE), se = stats::sd(values, na.rm = TRUE),
      bias = mean(values, na.rm = TRUE) - native[[index]],
      lower = interval[[1L]], upper = interval[[2L]], stringsAsFactors = FALSE
    )
  })) else data.frame()
  if (nrow(summary)) {
    display_order <- c(
      .nm_numbered_names("THETA", length(fit$theta)),
      .nm_numbered_names("OMEGA", length(fit$omega)),
      .nm_numbered_names("SIGMA", length(fit$sigma))
    )
    summary <- summary[match(display_order, summary$parameter), , drop = FALSE]
    rownames(summary) <- NULL
  }
  structure(list(estimates = estimates, errors = errors[nzchar(errors)],
                 summary = summary, n = n, successful = nrow(estimates),
                 convergence = convergence, replicate_seeds = replicate_seeds,
                 seed = seed, level = level, type = type, unit = unit,
                 strata = if (nzchar(strata)) strata else NULL,
                 cluster = if (nzchar(cluster)) cluster else NULL,
                 fit_fingerprint = .nm_fit_fingerprint(fit)),
            class = "nm_bootstrap")
}
