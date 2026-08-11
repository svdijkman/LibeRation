fidelity_correlated_fixture <- function() {
  fixture <- estimation_fixture(TRUE)
  fixture$data$WT <- rep(c(55, 70, 90), each = 4L)
  model <- nm_model(
    INPUT = c(names(fixture$data)), ADVAN = 1, DOSECMP = 1, OBSCMP = 1,
    PRED = paste(
      "CL=THETA(1)*(WT/70)^THETA(3)*exp(ETA(1))",
      "V=THETA(2)*exp(ETA(2))", "S1=V", sep = ";"
    ),
    ERROR = "Y=F+ERR(1)+F*ERR(2)",
    THETAS = data.frame(
      THETA = 1:3, Value = c(2, 20, 0.75), FIX = TRUE
    ),
    OMEGAS = data.frame(
      OMEGA = 1:3, ROW = c(1, 2, 2), COL = c(1, 1, 2),
      Value = c(0.09, 0.02, 0.16), FIX = TRUE
    ),
    SIGMAS = data.frame(
      SIGMA = 1:2, Value = c(0.04, 0.01), FIX = TRUE
    ),
    LIK_CONFIG = nm_lik_config(
      error = "combined", sigma_parameterization = "variance"
    )
  )
  list(model = model, data = fixture$data)
}

fidelity_scenario_gate <- function(scenario) {
  selected <- Sys.getenv("LIBER_FIDELITY_SCENARIO", "")
  if (nzchar(selected) && !identical(selected, scenario)) {
    skip(paste("fidelity runner selected", selected))
  }
}

test_that("dual numerical policies retain estimator targets with correlated ETAs", {
  fidelity_scenario_gate("correlated-two-eta")
  fixture <- fidelity_correlated_fixture()
  controls <- list(
    FO = list(), FOCE = list(), FOCEI = list(), LAPLACE = list(),
    ITS = list(
      its_mstep_schedule = "fixed", its_acceleration = "none",
      its_eta_schedule = "fixed"
    ),
    GQ = list(gq_order = 3L, gq_grid = "tensor"),
    IMP = list(
      n_imp = 12L, seed = 8101L, imp_sampling = "random",
      imp_proposal = "gaussian", imp_proposal_curvature = "exact",
      imp_sample_schedule = "fixed", imp_mstep_schedule = "fixed",
      imp_mstep_maxit = 1L, imp_reuse_modes = FALSE,
      imp_auto_stop = FALSE, imp_subject_allocation = "fixed"
    ),
    SAEM = list(
      n_iter = 4L, burn = 1L, mcmc_steps = 1L, seed = 8102L,
      saem_kernel = "random_walk", auto_stop = FALSE,
      saem_mstep_interval_burn = 1L, saem_mstep_interval = 1L,
      saem_parameter_averaging = "none", replicate_score_samples = 12L
    )
  )
  for (method in names(controls)) {
    common <- c(list(
      model = fixture$model, data = fixture$data, method = method,
      maxit = 2L, eta_maxit = 40L, tolerance = 1e-7,
      collect_output = FALSE
    ), controls[[method]])
    compatible <- do.call(
      nm_est, c(common, list(numerical_mode = "nonmem_compatibility"))
    )
    optimized <- do.call(
      nm_est, c(common, list(numerical_mode = "liber_optimized"))
    )
    expect_equal(
      optimized$objective, compatible$objective,
      tolerance = if (method %in% c("IMP", "SAEM")) 2e-5 else 2e-8,
      info = method
    )
    expect_equal(optimized$eta, compatible$eta, tolerance = 2e-5, info = method)
  }
})

test_that("posterior and nonparametric families cover correlated ETA structure", {
  selected <- Sys.getenv("LIBER_FIDELITY_SCENARIO", "")
  if (nzchar(selected) && !selected %in%
      c("correlated-posterior", "correlated-nonparametric")) {
    skip(paste("fidelity runner selected", selected))
  }
  fixture <- fidelity_correlated_fixture()
  if (!identical(selected, "correlated-nonparametric")) {
    fixture$model$THETAS$FIX[[1L]] <- FALSE
    for (method in c("HMC", "NUTS")) {
      for (mode in c("nonmem_compatibility", "liber_optimized")) {
        fit <- nm_est(
          fixture$model, fixture$data, method = method,
          n_warmup = 5L, n_sample = 6L, n_chains = 1L,
          n_leapfrog = 2L, max_depth = 2L, seed = 8103L,
          collect_output = FALSE, numerical_mode = mode
        )
        expect_true(all(is.finite(fit$chain[, "LOG_POSTERIOR"])),
                    info = paste(method, mode))
        expect_false(fit$objective_comparable)
      }
    }
  }

  if (!identical(selected, "correlated-posterior")) {
    supports <- rbind(c(-0.3, -0.3), c(-0.3, 0.3),
                      c(0.3, -0.3), c(0.3, 0.3))
    for (method in c("NPML", "NPAG")) {
      for (mode in c("nonmem_compatibility", "liber_optimized")) {
        fit <- nm_est(
          fixture$model, fixture$data, method = method,
          np_supports = supports, np_points = nrow(supports), np_cycles = 1L,
          np_max_support = 20L, np_weight_maxit = 100L,
          np_estimate_population = FALSE, seed = 8104L,
          collect_output = FALSE, numerical_mode = mode
        )
        expect_true(is.finite(fit$objective), info = paste(method, mode))
        expect_equal(sum(fit$nonparametric$weights), 1, tolerance = 1e-8)
        if (method == "NPML") {
          expect_equal(nrow(fit$nonparametric$supports), nrow(supports))
        }
      }
    }
  }
})

test_that("IOV and ODE estimator families retain finite accepted modes", {
  selected <- Sys.getenv("LIBER_FIDELITY_SCENARIO", "")
  if (nzchar(selected) && !selected %in% c("iov-two-occasion", "ode-advan6")) {
    skip(paste("fidelity runner selected", selected))
  }
  iov_model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "DV", "OCC"), ADVAN = 1,
    PRED = "CL=THETA(1)*exp(ETA(1)+ETA(2));V=THETA(2);S1=V",
    ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(THETA = 1:2, Value = c(2, 20), FIX = TRUE),
    OMEGAS = data.frame(OMEGA = 1:2, Value = c(0.1, 0.05), FIX = TRUE),
    SIGMAS = data.frame(SIGMA = 1, Value = 0.2, FIX = TRUE), IOV = 1
  )
  iov_data <- data.frame(
    ID = c(1, 1, 1, 1, 2, 2, 2, 2),
    TIME = rep(c(0, 1, 24, 25), 2), EVID = rep(c(1, 0, 1, 0), 2),
    AMT = rep(c(100, 0, 100, 0), 2),
    DV = c(NA, 4.5, NA, 4.4, NA, 4.7, NA, 4.6),
    OCC = rep(c(1, 1, 2, 2), 2)
  )
  if (!identical(selected, "ode-advan6")) for (method in c("FOCEI", "LAPLACE", "IMP", "SAEM")) {
    extra <- switch(
      method,
      IMP = list(n_imp = 10L, imp_mstep_maxit = 1L, seed = 8110L),
      SAEM = list(
        n_iter = 4L, burn = 1L, mcmc_steps = 1L, mstep_maxit = 1L,
        seed = 8111L, replicate_score_samples = 10L
      ),
      list()
    )
    fit <- do.call(nm_est, c(list(
      model = iov_model, data = iov_data, method = method,
      maxit = 2L, eta_maxit = 30L, collect_output = FALSE,
      numerical_mode = "liber_optimized"
    ), extra))
    expect_true(is.finite(fit$objective), info = paste("IOV", method))
    expect_true(all(is.finite(fit$eta)), info = paste("IOV", method))
  }

  ode_model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "DV", "MDV"),
    ADVAN = 6, DOSECMP = 1, OBSCMP = 1,
    PRED = "K=THETA(1)*exp(ETA(1));V=THETA(2);S1=V",
    DES = "DADT(1)=-K*A(1)", ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(THETA = 1:2, Value = c(0.1, 20), FIX = TRUE),
    OMEGAS = data.frame(OMEGA = 1, Value = 0.1, FIX = TRUE),
    SIGMAS = data.frame(SIGMA = 1, Value = 0.2, FIX = TRUE)
  )
  subject <- data.frame(
    TIME = c(0, 1, 4), EVID = c(1, 0, 0), AMT = c(100, 0, 0),
    DV = c(NA, 4.52, 3.35), MDV = c(1, 0, 0)
  )
  ode_data <- rbind(transform(subject, ID = 1L), transform(subject, ID = 2L))
  ode_data <- ode_data[, c("ID", setdiff(names(ode_data), "ID"))]
  if (!identical(selected, "iov-two-occasion")) for (method in c("FOCEI", "LAPLACE", "IMP")) {
    extra <- if (method == "IMP") {
      list(n_imp = 10L, imp_mstep_maxit = 1L, seed = 8120L)
    } else list()
    fit <- do.call(nm_est, c(list(
      model = ode_model, data = ode_data, method = method,
      maxit = 2L, eta_maxit = 30L, collect_output = FALSE,
      numerical_mode = "liber_optimized"
    ), extra))
    expect_true(is.finite(fit$objective), info = paste("ODE", method))
    expect_true(all(is.finite(fit$eta)), info = paste("ODE", method))
  }
})

test_that("BLQ and compiled user-likelihood estimators retain their targets", {
  selected <- Sys.getenv("LIBER_FIDELITY_SCENARIO", "")
  if (nzchar(selected) && !selected %in% c("blq-m3", "user-likelihood")) {
    skip(paste("fidelity runner selected", selected))
  }
  fixture <- estimation_fixture(TRUE)
  fixture$data$BLQ <- as.integer(
    !is.na(fixture$data$DV) & fixture$data$DV < 4.5
  )
  fixture$model$INPUT <- unique(c(fixture$model$INPUT, "BLQ"))
  fixture$model$LIK_CONFIG <- nm_lik_config(
    error = "additive", blq_method = "m3", lloq = 4.5,
    sigma_parameterization = "variance"
  )
  if (!identical(selected, "user-likelihood")) for (method in c("LAPLACE", "IMP", "SAEM")) {
    extra <- switch(
      method,
      IMP = list(n_imp = 10L, imp_mstep_maxit = 1L, seed = 8130L),
      SAEM = list(
        n_iter = 4L, burn = 1L, mcmc_steps = 1L, mstep_maxit = 1L,
        seed = 8131L, replicate_score_samples = 10L
      ),
      list()
    )
    fits <- lapply(c("nonmem_compatibility", "liber_optimized"), function(mode) {
      do.call(nm_est, c(list(
        model = fixture$model, data = fixture$data, method = method,
        maxit = 2L, eta_maxit = 30L, collect_output = FALSE,
        numerical_mode = mode
      ), extra))
    })
    expect_true(all(vapply(fits, function(x) is.finite(x$objective), logical(1))))
  }

  likelihood_model <- nm_model(
    INPUT = c("ID", "TIME", "DV", "MDV", "EVID", "AMT", "CMT"),
    ADVAN = 1,
    PRED = paste(
      "PBASE=1/(1+exp(-THETA(1)))",
      "P01=1/(1+exp(-THETA(2)))", "P11=1/(1+exp(-THETA(3)))",
      "CL=1", "V=1", "S1=V", "F=PBASE", sep = ";"
    ),
    ERROR = paste(
      "PCURRENT=ifelse(",
      "FIRST==1,ifelse(DV==1,PBASE,1-PBASE),",
      "ifelse(PREV_DV==0,ifelse(DV==1,P01,1-P01),",
      "ifelse(DV==1,P11,1-P11)))",
      "LOGLIK=log(pmax(PCURRENT,1e-12))", sep = "\n"
    ),
    THETAS = data.frame(
      THETA = 1:3, Value = qlogis(c(0.4, 0.8, 0.7)),
      LOWER = -10, UPPER = 10, FIX = FALSE
    ),
    LIK_CONFIG = nm_lik_config(error = "likelihood")
  )
  likelihood_data <- data.frame(
    ID = rep(1:2, each = 3), TIME = rep(0:2, 2),
    DV = c(0, 1, 1, 1, 0, 1), MDV = rep(c(1, 0, 0), 2),
    EVID = 0, AMT = 0, CMT = 1
  )
  if (!identical(selected, "blq-m3")) for (method in c("LAPLACE", "IMP", "SAEM", "BAYES")) {
    extra <- switch(
      method,
      IMP = list(n_imp = 10L, imp_mstep_maxit = 1L, seed = 8140L),
      SAEM = list(
        n_iter = 4L, burn = 1L, mcmc_steps = 1L, mstep_maxit = 1L,
        seed = 8141L, replicate_score_samples = 10L
      ),
      BAYES = list(n_burn = 2L, n_sample = 5L, seed = 8142L),
      list()
    )
    fit <- do.call(nm_est, c(list(
      model = likelihood_model, data = likelihood_data, method = method,
      maxit = 2L, eta_maxit = 20L, collect_output = FALSE,
      numerical_mode = "liber_optimized"
    ), extra))
    expect_true(is.finite(fit$objective), info = method)
  }
})

test_that("near-boundary variance models keep explicit estimator diagnostics", {
  fidelity_scenario_gate("boundary-variance")
  fixture <- estimation_fixture(FALSE)
  fixture$model$OMEGAS$Value[] <- 1e-4
  fixture$model$SIGMAS$Value[] <- 5e-3
  for (method in c("FO", "FOCEI")) {
    for (mode in c("nonmem_compatibility", "liber_optimized")) {
      fit <- nm_est(
        fixture$model, fixture$data, method = method, maxit = 2L,
        eta_maxit = 40L, collect_output = FALSE, numerical_mode = mode
      )
      expect_true(is.finite(fit$objective), info = paste(method, mode))
      expect_identical(fit$objective_comparable, TRUE)
      if (method == "FOCEI") {
        expect_true(all(fit$diagnostics$eta_convergence == 0L))
      }
    }
  }
})
