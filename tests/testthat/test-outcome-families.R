outcome_fit <- function(model, data) {
  data <- nm_dataset(data)
  structure(list(
    model = model, data = data, method = "LAPLACE",
    theta = model$THETAS$Value, sigma = model$SIGMAS$Value,
    omega = model$OMEGAS$Value,
    eta = matrix(numeric(), length(unique(data$.ID_INDEX)), model$n_eta),
    objective = NA_real_, convergence = 0L
  ), class = "nm_fit")
}

test_that("first-class Bernoulli and multicategory outcomes compile and simulate", {
  binary <- nm_model(
    INPUT = c("ID", "TIME", "DV", "MDV"), ADVAN = 1,
    PRED = "P=1/(1+exp(-THETA(1))); CL=1; V=1; S1=V; F=P",
    THETAS = data.frame(THETA = 1, Value = 0, LOWER = -10, UPPER = 10),
    OUTCOMES = nm_outcome("bernoulli", prediction = "P")
  )
  data <- data.frame(ID = rep(1:2, each = 3), TIME = rep(0:2, 2),
                     DV = c(0, 1, 1, 1, 0, 1), MDV = 0)
  objective <- nm_objective(binary, data, gradient = TRUE)
  expect_equal(objective$value, -2 * 6 * log(0.5), tolerance = 1e-10)
  expect_true(all(is.finite(objective$gradient)))
  simulated <- nm_simulate(binary, data, residual = TRUE, nsim = 2, seed = 11)
  expect_true(all(simulated$DV %in% 0:1))

  categorical <- nm_model(
    INPUT = c("ID", "TIME", "DV", "MDV"), ADVAN = 1,
    PRED = paste(
      "E1=exp(THETA(1)); E2=exp(THETA(2)); DEN=1+E1+E2",
      "P0=1/DEN; P1=E1/DEN; P2=E2/DEN",
      "CL=1; V=1; S1=V; F=P1", sep = "\n"
    ),
    THETAS = data.frame(THETA = 1:2, Value = c(0, 0), LOWER = -10, UPPER = 10),
    OUTCOMES = nm_outcome(
      "categorical", prediction = "P1",
      probabilities = c("P0", "P1", "P2"), categories = 0:2
    )
  )
  cat_data <- data.frame(ID = rep(1:3, each = 4), TIME = rep(0:3, 3),
                         DV = rep(0:2, length.out = 12), MDV = 0)
  expect_equal(nm_objective(categorical, cat_data, gradient = FALSE)$value,
               -2 * 12 * log(1 / 3), tolerance = 1e-10)
  cat_sim <- nm_simulate(categorical, cat_data, residual = TRUE, seed = 12)
  expect_true(all(cat_sim$DV %in% 0:2))
  fit <- outcome_fit(categorical, cat_data)
  diagnostic <- nm_outcome_diagnostics(fit)
  expect_true(all(diagnostic$PRED_CATEGORY %in% 0:2))
  expect_true(all(is.finite(diagnostic$OBSERVED_PROBABILITY)))
  vpc <- nm_vpc_categorical(fit, nsim = 20, seed = 13)
  expect_setequal(vpc$categories, 0:2)
  expect_setequal(unique(vpc$simulated$CATEGORY), as.character(0:2))

  irt <- nm_irt_outcomes(
    item_dvid = c(1, 2),
    probabilities = list(c("I1P0", "I1P1", "I1P2"),
                         c("I2P0", "I2P1", "I2P2")),
    categories = 0:2, names = c("Pain", "Function")
  )
  expect_s3_class(irt, "nm_outcomes")
  expect_equal(unname(vapply(irt, `[[`, character(1), "family")),
               c("ordinal", "ordinal"))
  expect_equal(unname(vapply(irt, `[[`, numeric(1), "dvid")), c(1, 2))
})

test_that("count families use normalized likelihoods and predictive checks", {
  poisson <- nm_model(
    INPUT = c("ID", "TIME", "DV", "MDV"), ADVAN = 1,
    PRED = "MU=exp(THETA(1)); CL=1; V=1; S1=V; F=MU",
    THETAS = data.frame(THETA = 1, Value = log(2), LOWER = -10, UPPER = 10),
    OUTCOMES = nm_outcome("poisson", prediction = "MU", max_count = 20)
  )
  data <- data.frame(ID = rep(1:4, each = 4), TIME = rep(0:3, 4),
                     DV = rep(c(0, 1, 2, 3), 4), MDV = 0)
  expected <- -2 * sum(data$DV * log(2) - 2 - lgamma(data$DV + 1))
  expect_equal(nm_objective(poisson, data, gradient = FALSE)$value,
               expected, tolerance = 1e-10)
  simulated <- nm_simulate(poisson, data, residual = TRUE, seed = 21)
  expect_true(all(simulated$DV >= 0 & simulated$DV == floor(simulated$DV)))
  fit <- outcome_fit(poisson, data)
  diagnostic <- nm_outcome_diagnostics(fit)
  expect_equal(diagnostic$EXPECTED, rep(2, nrow(data)))
  vpc <- nm_vpc_count(fit, nsim = 20, seed = 22)
  expect_s3_class(vpc, "nm_vpc_count")
  expect_true(all(c("MEAN_median", "ZERO_median") %in% names(vpc$simulated)))

  negative_binomial <- nm_model(
    INPUT = c("ID", "TIME", "DV", "MDV"), ADVAN = 1,
    PRED = "MU=exp(THETA(1)); SIZE=exp(THETA(2)); CL=1; V=1; S1=V; F=MU",
    THETAS = data.frame(THETA = 1:2, Value = log(c(2, 3)), LOWER = -10, UPPER = 10),
    OUTCOMES = nm_outcome(
      "negative_binomial", prediction = "MU", dispersion = "SIZE", max_count = 20
    )
  )
  score <- nm_objective(negative_binomial, data, gradient = TRUE)
  expect_true(is.finite(score$value))
  expect_true(all(is.finite(score$gradient)))
})

test_that("joint DVID outcomes share one compiled population likelihood", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "DV", "MDV", "DVID"), ADVAN = 1,
    PRED = paste(
      "MEAN=THETA(1); PBIN=1/(1+exp(-THETA(2)))",
      "CL=1; V=1; S1=V; F=MEAN", sep = "\n"
    ),
    THETAS = data.frame(THETA = 1:2, Value = c(2, 0), LOWER = -10, UPPER = 10),
    SIGMAS = data.frame(SIGMA = 1, Value = 0.5),
    OUTCOMES = nm_outcomes(
      nm_outcome("normal", name = "continuous", dvid = 1,
                 prediction = "MEAN", scale = "SIGMA(1)"),
      nm_outcome("bernoulli", name = "binary", dvid = 2, prediction = "PBIN")
    )
  )
  data <- data.frame(
    ID = rep(1:2, each = 4), TIME = rep(c(0, 0, 1, 1), 2),
    DVID = rep(c(1, 2), 4), DV = c(2.1, 1, 1.8, 0, 2.2, 1, 2, 1), MDV = 0
  )
  score <- nm_objective(model, data, gradient = TRUE)
  expect_true(is.finite(score$value))
  expect_true(all(is.finite(score$gradient)))
  simulated <- nm_simulate(model, data, residual = TRUE, seed = 31)
  expect_true(all(simulated$DV[simulated$DVID == 2] %in% 0:1))
})

test_that("event and continuous-time Markov outcomes fit and simulate", {
  tte <- nm_model(
    INPUT = c("ID", "TIME", "DV", "MDV"), ADVAN = 1,
    PRED = "HAZ=exp(THETA(1)); CL=1; V=1; S1=V; F=HAZ",
    THETAS = data.frame(THETA = 1, Value = log(0.2), LOWER = -10, UPPER = 10),
    OUTCOMES = nm_outcome("tte", prediction = "HAZ")
  )
  event_data <- data.frame(
    ID = rep(1:4, each = 5), TIME = rep(0:4, 4),
    DV = c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0),
    MDV = 0
  )
  score <- nm_objective(tte, event_data, gradient = TRUE)
  expect_true(is.finite(score$value))
  expect_true(all(is.finite(score$gradient)))
  tte_fit <- outcome_fit(tte, event_data)
  expect_s3_class(nm_vpc_tte(tte_fit, nsim = 20, seed = 41), "nm_vpc_tte")

  recurrent <- tte
  recurrent$OUTCOMES[[1]]$family <- "recurrent_event"
  recurrent_fit <- outcome_fit(recurrent, event_data)
  expect_s3_class(
    nm_vpc_recurrent(recurrent_fit, nsim = 20, seed = 42),
    "nm_vpc_recurrent"
  )

  ctmc <- nm_model(
    INPUT = c("ID", "TIME", "DV", "MDV"), ADVAN = 1,
    PRED = paste(
      "PI0=0.7; PI1=0.3; Q01=exp(THETA(1)); Q10=exp(THETA(2))",
      "CL=1; V=1; S1=V; F=PI1", sep = "\n"
    ),
    THETAS = data.frame(THETA = 1:2, Value = log(c(0.2, 0.4)), LOWER = -10, UPPER = 10),
    OUTCOMES = nm_outcome(
      "continuous_time_markov", prediction = "PI1", categories = c(0, 1),
      initial = c("PI0", "PI1"), rates = c("Q01", "Q10")
    )
  )
  states <- data.frame(ID = rep(1:3, each = 4), TIME = rep(0:3, 3),
                       DV = c(0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0), MDV = 0)
  ctmc_score <- nm_objective(ctmc, states, gradient = TRUE)
  expect_true(is.finite(ctmc_score$value))
  expect_true(all(is.finite(ctmc_score$gradient)))
  simulated <- nm_simulate(ctmc, states, residual = TRUE, seed = 43)
  expect_true(all(simulated$DV %in% 0:1))

  general_ctmc <- nm_model(
    INPUT = c("ID", "TIME", "DV", "MDV"), ADVAN = 1,
    PRED = paste(
      "PI0=0.5; PI1=0.3; PI2=0.2",
      "Q01=exp(THETA(1)); Q02=exp(THETA(2))",
      "Q10=exp(THETA(3)); Q12=exp(THETA(4))",
      "Q20=exp(THETA(5)); Q21=exp(THETA(6))",
      "CL=1; V=1; S1=V; F=0", sep = "\n"
    ),
    THETAS = data.frame(
      THETA = 1:6, Value = log(c(0.18, 0.05, 0.11, 0.09, 0.04, 0.16)),
      LOWER = -10, UPPER = 3
    ),
    OUTCOMES = nm_outcome(
      "continuous_time_markov", prediction = "PI0", categories = c(0, 1, 2),
      initial = c("PI0", "PI1", "PI2"),
      rates = matrix(c(
        "", "Q01", "Q02", "Q10", "", "Q12", "Q20", "Q21", ""
      ), 3, 3, byrow = TRUE)
    )
  )
  general_data <- data.frame(
    ID = "A", TIME = c(0, 0.4, 1.7, 3.2, 5.8),
    DV = c(0, 0, 1, 2, 1), MDV = 0
  )
  general_score <- nm_objective(general_ctmc, general_data, gradient = TRUE)
  rates <- exp(general_ctmc$THETAS$Value)
  generator <- matrix(c(
    0, rates[[1]], rates[[2]], rates[[3]], 0, rates[[4]],
    rates[[5]], rates[[6]], 0
  ), 3, 3, byrow = TRUE)
  diag(generator) <- -rowSums(generator)
  expected <- log(c(0.5, 0.3, 0.2)[general_data$DV[[1]] + 1L])
  for (row in 2:nrow(general_data)) {
    transition <- LibeRation:::.liberation_matrix_exp_pade(
      generator, general_data$TIME[[row]] - general_data$TIME[[row - 1L]]
    )
    expected <- expected + log(transition[
      general_data$DV[[row - 1L]] + 1L, general_data$DV[[row]] + 1L
    ])
  }
  expect_equal(general_score$value, -2 * expected, tolerance = 3e-10)
  expect_true(all(is.finite(general_score$gradient)))
  expect_true(isTRUE(general_ctmc$HMM_CONFIG$observed_states))
  decoded <- nm_hmm_decode(general_ctmc, general_data)
  expect_equal(decoded$HMM_STATE_INDEX, general_data$DV + 1L)
  simulated_general <- nm_simulate(
    general_ctmc, general_data, residual = TRUE, nsim = 20, seed = 44
  )
  expect_true(all(simulated_general$DV %in% 0:2))
})

test_that("competing-risk predictive checks use all declared causes", {
  model <- nm_model(
    INPUT = c("ID", "TIME", "DV", "MDV"), ADVAN = 1,
    PRED = "H1=exp(THETA(1)); H2=exp(THETA(2)); CL=1; V=1; S1=V; F=H1+H2",
    THETAS = data.frame(THETA = 1:2, Value = log(c(0.1, 0.05)), LOWER = -10, UPPER = 10),
    OUTCOMES = nm_outcome(
      "competing_risks", prediction = "H1", cause_hazards = c(`1` = "H1", `2` = "H2")
    )
  )
  data <- data.frame(
    ID = rep(1:6, each = 5), TIME = rep(0:4, 6), MDV = 0,
    DV = c(
      0, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0
    )
  )
  fit <- outcome_fit(model, data)
  result <- nm_vpc_competing(fit, nsim = 20, seed = 51)
  expect_s3_class(result, "nm_vpc_competing")
  expect_setequal(unique(result$observed$CAUSE), c("1", "2"))
  expect_setequal(unique(result$simulated$CAUSE), c("1", "2"))
})
