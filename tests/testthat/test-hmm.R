hmm_test_model <- function(log_scale = FALSE, by_dvid = TRUE) {
  probability <- paste(
    "I1 = 1 / (1 + exp(-THETA(1)))",
    "I2 = 1 - I1",
    "T11 = 1 / (1 + exp(-THETA(2)))",
    "T12 = 1 - T11",
    "T21 = 1 / (1 + exp(-THETA(3)))",
    "T22 = 1 - T21",
    "E10 = 1 / (1 + exp(-THETA(4)))",
    "E20 = 1 / (1 + exp(-THETA(5)))",
    "E1 = ifelse(DV == 0, E10, 1 - E10)",
    "E2 = ifelse(DV == 0, E20, 1 - E20)",
    sep = "\n"
  )
  if (log_scale) {
    probability <- paste(
      probability,
      "LI1 = log(I1)", "LI2 = log(I2)",
      "LT11 = log(T11)", "LT12 = log(T12)",
      "LT21 = log(T21)", "LT22 = log(T22)",
      "LE1 = log(E1)", "LE2 = log(E2)", sep = "\n"
    )
  }
  nm_model(
    INPUT = c("ID", "TIME", "DV", "MDV", "DVID"),
    ADVAN = 1,
    PRED = "CL=1; V=1; S1=V; F=0",
    ERROR = probability,
    THETAS = data.frame(
      THETA = 1:5,
      Value = stats::qlogis(c(0.6, 0.8, 0.3, 0.9, 0.2)),
      LOWER = -10, UPPER = 10
    ),
    HMM_CONFIG = nm_hmm_config(
      states = c("low", "high"),
      initial = if (log_scale) c("LI1", "LI2") else c("I1", "I2"),
      transition = matrix(
        if (log_scale) c("LT11", "LT12", "LT21", "LT22") else
          c("T11", "T12", "T21", "T22"),
        2, 2, byrow = TRUE
      ),
      emission = if (log_scale) c("LE1", "LE2") else c("E1", "E2"),
      initial_scale = if (log_scale) "log" else "probability",
      transition_scale = if (log_scale) "log" else "probability",
      emission_scale = if (log_scale) "log" else "likelihood",
      by_dvid = by_dvid
    )
  )
}

hmm_manual_loglik <- function(outcome) {
  initial <- c(0.6, 0.4)
  transition <- matrix(c(0.8, 0.2, 0.3, 0.7), 2, byrow = TRUE)
  emission_zero <- c(0.9, 0.2)
  alpha <- NULL
  value <- 0
  for (index in seq_along(outcome)) {
    prior <- if (index == 1L) initial else drop(alpha %*% transition)
    emission <- if (outcome[[index]] == 0) emission_zero else 1 - emission_zero
    weight <- prior * emission
    value <- value + log(sum(weight))
    alpha <- weight / sum(weight)
  }
  value
}

hmm_exact_paths <- function(outcome) {
  initial <- c(0.6, 0.4)
  transition <- matrix(c(0.8, 0.2, 0.3, 0.7), 2, byrow = TRUE)
  emission_zero <- c(0.9, 0.2)
  paths <- as.matrix(expand.grid(rep(list(1:2), length(outcome))))
  weight <- apply(paths, 1L, function(path) {
    value <- initial[path[[1L]]]
    for (time in seq_along(outcome)) {
      emission <- if (outcome[[time]] == 0) emission_zero else 1 - emission_zero
      value <- value * emission[path[[time]]]
      if (time < length(outcome)) {
        value <- value * transition[path[[time]], path[[time + 1L]]]
      }
    }
    value
  })
  likelihood <- sum(weight)
  posterior <- weight / likelihood
  smoothed <- vapply(seq_along(outcome), function(time) {
    vapply(1:2, function(state) sum(posterior[paths[, time] == state]), numeric(1))
  }, numeric(2))
  list(
    log_likelihood = log(likelihood),
    smoothed = t(smoothed),
    viterbi = as.integer(paths[which.max(weight), ]),
    viterbi_log_joint = log(max(weight))
  )
}

cthmm_test_model <- function(log_rate = FALSE) {
  rates <- c(0.18, 0.05, 0.11, 0.09, 0.04, 0.16)
  rate_code <- if (log_rate) {
    paste0("Q", c("12", "13", "21", "23", "31", "32"),
           " = THETA(", seq_along(rates), ")")
  } else {
    paste0("Q", c("12", "13", "21", "23", "31", "32"),
           " = exp(THETA(", seq_along(rates), "))")
  }
  error <- paste(c(
    "I1 = 0.5", "I2 = 0.3", "I3 = 0.2", rate_code,
    "E1 = ifelse(DV == 0, 0.90, ifelse(DV == 1, 0.08, 0.02))",
    "E2 = ifelse(DV == 0, 0.10, ifelse(DV == 1, 0.80, 0.10))",
    "E3 = ifelse(DV == 0, 0.04, ifelse(DV == 1, 0.16, 0.80))"
  ), collapse = "\n")
  nm_model(
    INPUT = c("ID", "TIME", "DV", "MDV"), ADVAN = 1,
    PRED = "CL=1; V=1; S1=V; F=0", ERROR = error,
    THETAS = data.frame(
      THETA = seq_along(rates), Value = log(rates), LOWER = -10, UPPER = 2
    ),
    HMM_CONFIG = nm_cthmm_config(
      states = c("mild", "moderate", "severe"),
      initial = c("I1", "I2", "I3"),
      generator = matrix(c(
        "", "Q12", "Q13", "Q21", "", "Q23", "Q31", "Q32", ""
      ), 3, 3, byrow = TRUE),
      emission = c("E1", "E2", "E3"),
      rate_scale = if (log_rate) "log" else "rate",
      by_dvid = FALSE
    )
  )
}

cthmm_manual_loglik <- function(time, outcome, theta) {
  rates <- exp(theta)
  generator <- matrix(c(
    0, rates[[1]], rates[[2]], rates[[3]], 0, rates[[4]],
    rates[[5]], rates[[6]], 0
  ), 3, 3, byrow = TRUE)
  diag(generator) <- -rowSums(generator)
  emission <- matrix(c(
    0.90, 0.08, 0.02, 0.10, 0.80, 0.10, 0.04, 0.16, 0.80
  ), 3, 3, byrow = TRUE)
  alpha <- NULL
  value <- 0
  for (index in seq_along(outcome)) {
    prior <- if (index == 1L) c(0.5, 0.3, 0.2) else {
      transition <- LibeRation:::.liberation_matrix_exp_pade(
        generator, time[[index]] - time[[index - 1L]]
      )
      drop(alpha %*% transition)
    }
    weight <- prior * emission[, outcome[[index]] + 1L]
    value <- value + log(sum(weight))
    alpha <- weight / sum(weight)
  }
  value
}

test_that("HMM configuration validates named forward-algorithm outputs", {
  model <- hmm_test_model()
  expect_s3_class(model$HMM_CONFIG, "nm_hmm_config")
  expect_identical(model$ERROR_TYPE, "likelihood")
  expect_identical(model$likelihood_scale, "hmm")
  expect_null(model$likelihood_output)
  expect_error(
    nm_hmm_config(
      states = c("a", "b"), initial = c("I1", "I2"),
      transition = matrix("T", 3, 3), emission = c("E1", "E2")
    ),
    "square character matrix"
  )
})

test_that("continuous-time HMM uses an arbitrary-state differentiable generator", {
  model <- cthmm_test_model()
  expect_s3_class(model$HMM_CONFIG, "nm_cthmm_config")
  expect_identical(model$HMM_CONFIG$transition_type, "continuous")
  data <- data.frame(
    ID = "A", TIME = c(0, 0.4, 1.7, 3.2, 5.8),
    DV = c(0, 0, 1, 2, 1), MDV = 0
  )
  objective <- nm_objective(model, data, gradient = TRUE)
  expected <- -2 * cthmm_manual_loglik(
    data$TIME, data$DV, model$THETAS$Value
  )
  expect_equal(objective$value, expected, tolerance = 2e-10)
  expect_true(all(is.finite(objective$gradient)))

  step <- 1e-5
  numeric_gradient <- vapply(seq_len(6), function(index) {
    plus <- minus <- model$THETAS$Value
    plus[[index]] <- plus[[index]] + step
    minus[[index]] <- minus[[index]] - step
    (nm_objective(model, data, theta = plus, gradient = FALSE)$value -
       nm_objective(model, data, theta = minus, gradient = FALSE)$value) /
      (2 * step)
  }, numeric(1))
  expect_equal(unname(objective$gradient), numeric_gradient, tolerance = 3e-5)

  decoded <- nm_hmm_decode(model, data, method = "all")
  expect_equal(
    rowSums(decoded[c(
      "HMM_FILTER_PROB_mild", "HMM_FILTER_PROB_moderate",
      "HMM_FILTER_PROB_severe"
    )]), rep(1, nrow(data)), tolerance = 1e-11
  )
  expect_true(all(decoded$HMM_VITERBI_STATE %in% model$HMM_CONFIG$states))

  log_rate_model <- cthmm_test_model(log_rate = TRUE)
  expect_equal(
    nm_objective(log_rate_model, data, gradient = FALSE)$value,
    objective$value, tolerance = 2e-10
  )
})

test_that("continuous-time HMM validates generator dimensions", {
  expect_error(
    nm_cthmm_config(
      states = c("a", "b", "c"), initial = c("I1", "I2", "I3"),
      generator = matrix("Q", 2, 2), emission = c("E1", "E2", "E3")
    ), "square character matrix"
  )
})

test_that("scaled HMM objective and exact gradient match the forward likelihood", {
  data <- data.frame(
    ID = "A", TIME = 0:4, DV = c(0, 0, 1, 1, 0), MDV = 0, DVID = 1
  )
  for (log_scale in c(FALSE, TRUE)) {
    model <- hmm_test_model(log_scale = log_scale)
    objective <- nm_objective(model, data, gradient = TRUE)
    expected <- -2 * hmm_manual_loglik(data$DV)
    expect_equal(objective$value, expected, tolerance = 1e-10)
    expect_length(objective$gradient, 5L)
    expect_true(all(is.finite(objective$gradient)))

    step <- 1e-5
    numeric_gradient <- vapply(seq_len(5), function(index) {
      plus <- minus <- model$THETAS$Value
      plus[[index]] <- plus[[index]] + step
      minus[[index]] <- minus[[index]] - step
      (nm_objective(model, data, theta = plus, gradient = FALSE)$value -
         nm_objective(model, data, theta = minus, gradient = FALSE)$value) / (2 * step)
    }, numeric(1))
    expect_equal(unname(objective$gradient), numeric_gradient, tolerance = 2e-5)
  }
})

test_that("log-weight HMM is equivalent and exposes filtered state probabilities", {
  data <- data.frame(
    ID = "A", TIME = 0:4, DV = c(0, 0, 1, 1, 0), MDV = 0, DVID = 1
  )
  probability_model <- hmm_test_model()
  log_model <- hmm_test_model(log_scale = TRUE)
  probability_value <- nm_objective(probability_model, data, gradient = FALSE)$value
  log_value <- nm_objective(log_model, data, gradient = FALSE)$value
  expect_equal(log_value, probability_value, tolerance = 1e-10)

  decoded <- nm_hmm_decode(log_model, data)
  expect_s3_class(decoded, "nm_hmm_decode")
  expect_true(all(c("HMM_STATE", "HMM_PROB_low", "HMM_PROB_high") %in% names(decoded)))
  expect_equal(
    rowSums(decoded[c("HMM_PROB_low", "HMM_PROB_high")]),
    rep(1, nrow(decoded)), tolerance = 1e-12
  )
  expect_equal(attr(decoded, "log_likelihood"), -0.5 * log_value, tolerance = 1e-10)
})

test_that("forward-backward smoothing and Viterbi agree with exact path enumeration", {
  data <- data.frame(
    ID = "A", TIME = 0:4, DV = c(0, 0, 1, 1, 0), MDV = 0, DVID = 1
  )
  exact <- hmm_exact_paths(data$DV)
  decoded <- nm_hmm_decode(hmm_test_model(), data, method = "all")
  expect_true(all(c(
    "HMM_FILTER_STATE", "HMM_SMOOTH_STATE", "HMM_VITERBI_STATE",
    "HMM_FILTER_PROB_low", "HMM_FILTER_PROB_high",
    "HMM_SMOOTH_PROB_low", "HMM_SMOOTH_PROB_high"
  ) %in% names(decoded)))
  expect_equal(
    unname(as.matrix(decoded[c("HMM_SMOOTH_PROB_low", "HMM_SMOOTH_PROB_high")])),
    exact$smoothed, tolerance = 1e-12
  )
  expect_equal(decoded$HMM_VITERBI_STATE_INDEX, exact$viterbi)
  expect_equal(attr(decoded, "log_likelihood"), exact$log_likelihood,
               tolerance = 1e-12)
  summary <- attr(decoded, "sequence_summary")
  expect_equal(summary$VITERBI_LOG_JOINT, exact$viterbi_log_joint,
               tolerance = 1e-12)
  expect_equal(
    summary$VITERBI_LOG_POSTERIOR,
    exact$viterbi_log_joint - exact$log_likelihood,
    tolerance = 1e-12
  )
  expect_identical(summary$ID, "A")
  expect_identical(summary$DVID, 1L)

  log_decoded <- nm_hmm_decode(
    hmm_test_model(log_scale = TRUE), data, method = "all"
  )
  expect_equal(
    unname(as.matrix(log_decoded[c(
      "HMM_SMOOTH_PROB_low", "HMM_SMOOTH_PROB_high"
    )])),
    exact$smoothed, tolerance = 1e-12
  )
  expect_equal(log_decoded$HMM_VITERBI_STATE_INDEX, exact$viterbi)

  smooth_only <- nm_hmm_decode(hmm_test_model(), data, method = "smoothed")
  expect_equal(
    unname(as.matrix(smooth_only[c("HMM_PROB_low", "HMM_PROB_high")])),
    exact$smoothed, tolerance = 1e-12
  )
  expect_identical(attr(smooth_only, "method"), "smoothed")
  path_only <- nm_hmm_decode(hmm_test_model(), data, method = "viterbi")
  expect_equal(path_only$HMM_STATE_INDEX, exact$viterbi)
  expect_false(any(grepl("^HMM_PROB_", names(path_only))))
})

test_that("HMM filters keep DVID sequences independent", {
  model <- hmm_test_model(by_dvid = TRUE)
  data <- data.frame(
    ID = "A", TIME = c(0, 0, 1, 1), DVID = c(1, 2, 1, 2),
    DV = c(0, 1, 0, 1), MDV = 0
  )
  decoded <- nm_hmm_decode(model, data)
  initial_low <- 0.6 * 0.9 / (0.6 * 0.9 + 0.4 * 0.2)
  initial_high_outcome <- 0.6 * 0.1 / (0.6 * 0.1 + 0.4 * 0.8)
  expect_equal(decoded$HMM_PROB_low[1:2], c(initial_low, initial_high_outcome),
               tolerance = 1e-12)
  expect_false(isTRUE(all.equal(decoded$HMM_PROB_low[[3]], initial_low)))
  expect_false(isTRUE(all.equal(decoded$HMM_PROB_low[[4]], initial_high_outcome)))

  all_decoders <- nm_hmm_decode(model, data, method = "all")
  first <- nm_hmm_decode(model, data[data$DVID == 1, ], method = "all")
  second <- nm_hmm_decode(model, data[data$DVID == 2, ], method = "all")
  expect_equal(
    all_decoders$HMM_SMOOTH_PROB_low[data$DVID == 1],
    first$HMM_SMOOTH_PROB_low, tolerance = 1e-12
  )
  expect_equal(
    all_decoders$HMM_SMOOTH_PROB_low[data$DVID == 2],
    second$HMM_SMOOTH_PROB_low, tolerance = 1e-12
  )
  expect_equal(
    all_decoders$HMM_VITERBI_STATE_INDEX[data$DVID == 1],
    first$HMM_VITERBI_STATE_INDEX
  )
  expect_equal(nrow(attr(all_decoders, "sequence_summary")), 2L)
})

test_that("HMMs fit through the population engine and expose state diagnostics", {
  outcomes <- rbind(
    c(0, 0, 0, 1, 1, 1), c(0, 0, 1, 1, 1, 0),
    c(1, 1, 1, 0, 0, 0), c(1, 1, 0, 0, 0, 1)
  )
  data <- data.frame(
    ID = rep(seq_len(nrow(outcomes)), each = ncol(outcomes)),
    TIME = rep(seq_len(ncol(outcomes)) - 1L, nrow(outcomes)),
    DV = as.numeric(t(outcomes)), MDV = 0L, DVID = 1L
  )
  fit <- nm_est(hmm_test_model(), data, method = "LAPLACE", maxit = 12L)
  expect_s3_class(fit, "nm_fit")
  expect_true(is.finite(fit$objective))
  expect_true(all(is.finite(fit$theta)))

  diagnostics <- nm_gof(fit)
  expect_true(all(c("HMM_STATE", "HMM_PROB_low", "HMM_PROB_high") %in%
                    names(diagnostics)))
  expect_true(all(is.na(diagnostics$CWRES)))
  expect_match(attr(diagnostics, "residual_note"), "hidden Markov")
})

test_that("HMM export fails instead of writing an invalid NONMEM stream", {
  expect_error(nm_control_write(hmm_test_model()), "HMM_CONFIG")
})
