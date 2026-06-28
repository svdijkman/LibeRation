#' Catalog of built-in synthetic pharmacometric examples
#'
#' @return Named list with \code{id}, \code{label}, and \code{description} per entry.
#' @examples
#' nm_synthetic_catalog()
#' @export
nm_synthetic_catalog <- function() {
  list(
    theo = list(
      id = "theo",
      label = "THEO — oral 2-compartment",
      description = paste0(
        "Theophylline-style oral 2-compartment PK (ADVAN 4, TRANS 4). ",
        "Classic teaching dataset (Schumitzky & Jelliffe; NONMEM THEO demo)."
      ),
      csv = "theo.csv"
    ),
    warf = list(
      id = "warf",
      label = "Warfarin — oral 1-compartment",
      description = paste0(
        "Oral 1-compartment warfarin-like PK (ADVAN 2, TRANS 1). ",
        "Inspired by Sheiner et al. (1979) population analyses."
      ),
      csv = "warf.csv"
    ),
    pheno = list(
      id = "pheno",
      label = "Phenobarbital — oral 1-compartment",
      description = paste0(
        "Neonatal phenobarbital-style oral PK (ADVAN 2, TRANS 1). ",
        "Simplified from Boeckmann et al. (2005) without covariates."
      ),
      csv = "pheno.csv"
    ),
    digox = list(
      id = "digox",
      label = "Digoxin — oral 1-compartment",
      description = paste0(
        "Oral digoxin-like PK with slow absorption (ADVAN 2, TRANS 1). ",
        "Teaching-scale dose (5 mg) and volume (V=10 L) for visible concentrations over 72 h."
      ),
      csv = "digox.csv"
    ),
    iv1 = list(
      id = "iv1",
      label = "IV bolus — 1-compartment",
      description = "Intravenous bolus, 1-compartment (ADVAN 1, TRANS 1). Single-dose PK.",
      csv = "iv1.csv"
    ),
    iv2 = list(
      id = "iv2",
      label = "IV bolus — 2-compartment",
      description = "Intravenous bolus, 2-compartment (ADVAN 3, TRANS 4).",
      csv = "iv2.csv"
    ),
    iv3 = list(
      id = "iv3",
      label = "IV bolus — 3-compartment",
      description = "Intravenous bolus, 3-compartment (ADVAN 11, TRANS 4).",
      csv = "iv3.csv"
    ),
    oral3 = list(
      id = "oral3",
      label = "Oral — 3-compartment",
      description = "Oral absorption, 3-compartment disposition (ADVAN 12, TRANS 4).",
      csv = "oral3.csv"
    ),
    multidose = list(
      id = "multidose",
      label = "Oral — multiple doses",
      description = "Oral 1-compartment with two bolus doses per subject (ADVAN 2, TRANS 1).",
      csv = "multidose.csv"
    ),
    steady = list(
      id = "steady",
      label = "Oral — steady-state dosing",
      description = "Oral 1-compartment at steady state (SS=1, II=12 h; ADVAN 2, TRANS 1).",
      csv = "steady.csv"
    ),
    infusion = list(
      id = "infusion",
      label = "IV infusion — 1-compartment",
      description = "Intravenous infusion (RATE > 0) with 1-compartment PK (ADVAN 1, TRANS 1).",
      csv = "infusion.csv"
    ),
    pkpd_emax = list(
      id = "pkpd_emax",
      label = "PK/PD — Emax link",
      description = "Oral PK with an Emax pharmacodynamic response linked to concentration (ADVAN 2).",
      csv = "pkpd_emax.csv"
    )
  )
}

#' @keywords internal
.nm_synthetic_model_from_spec <- function(advan,
                                            trans,
                                            dose_cmp,
                                            obs_cmp,
                                            input = NULL,
                                            use_ode = FALSE,
                                            des = "",
                                            ss = 0L,
                                            error = NULL) {
  spec <- .nm_ctl_template_spec(advan, trans)
  spec$pk <- .nm_ctl_append_pk_scaling(spec$pk, advan, trans)
  if (is.null(input)) {
    input <- c("AMT", "TIME", "EVID", "CMT", "MDV", "DV")
    if (advan %in% c(2L, 4L, 12L)) {
      input <- c(input, "F1", "S1", "KA")
    }
    if (advan == 4L) {
      input <- c(input, "S2")
    }
  }
  nm_model(
    INPUT = input,
    ADVAN = advan,
    TRANS = trans,
    DOSECMP = dose_cmp,
    OBSCMP = obs_cmp,
    PRED = spec$pk,
    ERROR = error %||% spec$error,
    DES = des,
    USE_ODE = isTRUE(use_ode),
    SS = as.integer(ss),
    THETAS = spec$thetas,
    OMEGAS = spec$omegas,
    SIGMAS = spec$sigmas
  )
}

#' @keywords internal
.nm_synthetic_fill_dv <- function(model, dat, etas, sigma, pk_engine = "cpp",
                                    use_error = FALSE) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    .nm_stop("Package 'data.table' is required for synthetic datasets.")
  }
  dat <- data.table::as.data.table(dat)
  n_eta <- if (is.matrix(etas)) ncol(etas) else length(etas[[1L]])
  s1 <- sigma[1]
  s2 <- if (length(sigma) >= 2L) sigma[2] else 0
  if (isTRUE(model$USE_ODE)) {
    pk_engine <- "cpp"
  }
  for (id in unique(dat$ID)) {
    sub <- dat[ID == id]
    eta_j <- if (is.matrix(etas)) {
      etas[id, , drop = TRUE]
    } else {
      etas[[id]]
    }
    if (length(eta_j) != n_eta) {
      eta_j <- rep(0, n_eta)
    }
    pred <- .nm_subject_ipred(
      model, sub,
      model$THETAS$Value,
      model$OMEGAS$Value,
      eta_j,
      model$SIGMAS$Value,
      pk_engine = pk_engine
    )
    f <- pred$F
    if (isTRUE(use_error) && nzchar(model$ERROR %||% "")) {
      mu <- vapply(seq_along(f), function(k) {
        fk <- if (is.list(f)) f[[k]] else f[k]
        err <- rep(0, length(model$SIGMAS$Value))
        y_out <- .nm_eval_error(
          model,
          model$THETAS$Value,
          model$OMEGAS$Value,
          model$SIGMAS$Value,
          eta_j,
          err,
          fk
        )
        as.numeric(.nm_extract_y(y_out))
      }, numeric(1L))
      dv <- mu * (1 + rnorm(length(mu), sd = s1)) + rnorm(length(mu), sd = s2)
    } else {
      dv <- f * (1 + rnorm(length(f), sd = s1)) + rnorm(length(f), sd = s2)
    }
    idx <- which(dat$ID == id & dat$MDV == 0L)
    dat$DV[idx] <- dv
  }
  nm_dataset_from_table(dat)
}

#' Generate a synthetic dataset by catalog id
#'
#' @param id Catalog id (see \code{\link{nm_synthetic_catalog}}).
#' @param n_sub Number of subjects.
#' @param seed Random seed.
#' @return List with \code{model} and \code{data}, like \code{\link{nm_synthetic_theo}}.
#' @examples
#' nm_synthetic_dataset("theo", n_sub = 2L)
#' @export
nm_synthetic_dataset <- function(id = c(
                                   "theo", "warf", "pheno", "digox",
                                   "iv1", "iv2", "iv3", "oral3",
                                   "multidose", "steady", "infusion", "pkpd_emax"
                                 ),
                                 n_sub = 10L,
                                 seed = 1L) {
  id <- match.arg(id)
  switch(
    id,
    theo = nm_synthetic_theo(n_sub = n_sub, seed = seed),
    warf = nm_synthetic_warf(n_sub = n_sub, seed = seed),
    pheno = nm_synthetic_pheno(n_sub = n_sub, seed = seed),
    digox = nm_synthetic_digox(n_sub = n_sub, seed = seed),
    iv1 = nm_synthetic_iv1(n_sub = n_sub, seed = seed),
    iv2 = nm_synthetic_iv2(n_sub = n_sub, seed = seed),
    iv3 = nm_synthetic_iv3(n_sub = n_sub, seed = seed),
    oral3 = nm_synthetic_oral3(n_sub = n_sub, seed = seed),
    multidose = nm_synthetic_multidose(n_sub = n_sub, seed = seed),
    steady = nm_synthetic_steady(n_sub = n_sub, seed = seed),
    infusion = nm_synthetic_infusion(n_sub = n_sub, seed = seed),
    pkpd_emax = nm_synthetic_pkpd_emax(n_sub = n_sub, seed = seed)
  )
}

#' @keywords internal
.nm_synthetic_oral1_cl <- function(n_sub,
                                    seed,
                                    cl,
                                    v,
                                    ka,
                                    omega,
                                    sigma,
                                    amt,
                                    times) {
  set.seed(seed)
  rows <- list()
  etas <- vector("list", n_sub)
  for (id in seq_len(n_sub)) {
    eta <- rnorm(length(omega), sd = sqrt(omega))
    etas[[id]] <- eta
    rows[[length(rows) + 1L]] <- data.frame(
      ID = id, TIME = 0, EVID = 1L, CMT = 1L, AMT = amt, RATE = 0,
      MDV = 1L, DV = 0, F1 = 1, S1 = 1
    )
    for (tm in times[-1]) {
      rows[[length(rows) + 1L]] <- data.frame(
        ID = id, TIME = tm, EVID = 0L, CMT = 2L, AMT = 0, RATE = 0,
        MDV = 0L, DV = NA_real_, F1 = 1, S1 = 1
      )
    }
  }
  dat <- do.call(rbind, rows)
  model <- .nm_synthetic_model_from_spec(2L, 1L, dose_cmp = 1L, obs_cmp = 2L)
  model$INPUT <- setdiff(model$INPUT, "KA")
  model$THETAS <- data.frame(
    THETA = 1:3,
    Value = c(cl / v, v, ka),
    stringsAsFactors = FALSE
  )
  model$OMEGAS <- data.frame(
    OMEGA = seq_along(omega),
    Value = omega,
    stringsAsFactors = FALSE
  )
  model$SIGMAS <- data.frame(
    SIGMA = seq_along(sigma),
    Value = sigma,
    stringsAsFactors = FALSE
  )
  list(
    model = model,
    data = .nm_synthetic_fill_dv(model, dat, etas, sigma)
  )
}

#' Synthetic warfarin-like oral 1-compartment dataset
#'
#' @inheritParams nm_synthetic_theo
#' @examples
#' nm_synthetic_warf(n_sub = 2L)
#' @export
nm_synthetic_warf <- function(n_sub = 10L, seed = 1L) {
  .nm_synthetic_oral1_cl(
    n_sub = n_sub,
    seed = seed,
    cl = 0.18,
    v = 11,
    ka = 0.45,
    omega = c(0.04, 0.04),
    sigma = c(0.08, 0.15),
    amt = 5,
    times = c(0, 1, 2, 4, 8, 12, 24, 36, 48, 72)
  )
}

#' Synthetic phenobarbital-like oral 1-compartment dataset
#'
#' @inheritParams nm_synthetic_theo
#' @examples
#' nm_synthetic_pheno(n_sub = 2L)
#' @export
nm_synthetic_pheno <- function(n_sub = 10L, seed = 1L) {
  .nm_synthetic_oral1_cl(
    n_sub = n_sub,
    seed = seed,
    cl = 0.35,
    v = 2.5,
    ka = 1.4,
    omega = c(0.06, 0.09),
    sigma = c(0.1, 0.4),
    amt = 30,
    times = c(0, 0.5, 1, 2, 4, 8, 12, 24, 36)
  )
}

#' Synthetic digoxin-like oral 1-compartment dataset
#'
#' @inheritParams nm_synthetic_theo
#' @examples
#' nm_synthetic_digox(n_sub = 2L)
#' @export
nm_synthetic_digox <- function(n_sub = 10L, seed = 1L) {
  .nm_synthetic_oral1_cl(
    n_sub = n_sub,
    seed = seed,
    cl = 0.09,
    v = 10,
    ka = 0.65,
    omega = c(0.25, 0.35),
    sigma = c(0.04, 0.06),
    amt = 5,
    times = c(0, 0.25, 0.5, 1, 2, 4, 8, 12, 24, 48, 72)
  )
}

#' Synthetic IV bolus 1-compartment dataset
#'
#' @inheritParams nm_synthetic_theo
#' @examples
#' nm_synthetic_iv1(n_sub = 2L)
#' @export
nm_synthetic_iv1 <- function(n_sub = 10L, seed = 1L) {
  set.seed(seed)
  omega <- c(0.05)
  sigma <- c(0.08)
  times <- c(0, 0.5, 1, 2, 4, 8, 12, 24)
  rows <- list()
  etas <- vector("list", n_sub)
  for (id in seq_len(n_sub)) {
    eta <- rnorm(1, sd = sqrt(omega))
    etas[[id]] <- eta
    rows[[length(rows) + 1L]] <- data.frame(
      ID = id, TIME = 0, EVID = 1L, CMT = 1L, AMT = 100, RATE = 0,
      MDV = 1L, DV = 0, SS = 0L, II = 0
    )
    for (tm in times[-1]) {
      rows[[length(rows) + 1L]] <- data.frame(
        ID = id, TIME = tm, EVID = 0L, CMT = 1L, AMT = 0, RATE = 0,
        MDV = 0L, DV = NA_real_, SS = 0L, II = 0
      )
    }
  }
  dat <- do.call(rbind, rows)
  model <- .nm_synthetic_model_from_spec(1L, 1L, dose_cmp = 1L, obs_cmp = 1L)
  list(model = model, data = .nm_synthetic_fill_dv(model, dat, etas, sigma))
}

#' Synthetic IV bolus 2-compartment dataset
#'
#' @inheritParams nm_synthetic_theo
#' @examples
#' nm_synthetic_iv2(n_sub = 2L)
#' @export
nm_synthetic_iv2 <- function(n_sub = 10L, seed = 1L) {
  set.seed(seed)
  omega <- c(0.06, 0.06)
  sigma <- c(0.1)
  times <- c(0, 0.25, 0.5, 1, 2, 4, 8, 12, 24, 48)
  rows <- list()
  etas <- vector("list", n_sub)
  for (id in seq_len(n_sub)) {
    eta <- rnorm(2, sd = sqrt(omega))
    etas[[id]] <- eta
    rows[[length(rows) + 1L]] <- data.frame(
      ID = id, TIME = 0, EVID = 1L, CMT = 1L, AMT = 500, RATE = 0,
      MDV = 1L, DV = 0, SS = 0L, II = 0
    )
    for (tm in times[-1]) {
      rows[[length(rows) + 1L]] <- data.frame(
        ID = id, TIME = tm, EVID = 0L, CMT = 1L, AMT = 0, RATE = 0,
        MDV = 0L, DV = NA_real_, SS = 0L, II = 0
      )
    }
  }
  dat <- do.call(rbind, rows)
  model <- .nm_synthetic_model_from_spec(3L, 4L, dose_cmp = 1L, obs_cmp = 1L)
  list(model = model, data = .nm_synthetic_fill_dv(model, dat, etas, sigma))
}

#' Synthetic IV bolus 3-compartment dataset
#'
#' @inheritParams nm_synthetic_theo
#' @examples
#' nm_synthetic_iv3(n_sub = 2L)
#' @export
nm_synthetic_iv3 <- function(n_sub = 10L, seed = 1L) {
  set.seed(seed)
  omega <- c(0.05, 0.05)
  sigma <- c(0.12)
  times <- c(0, 0.5, 1, 2, 4, 8, 12, 24, 48, 72)
  rows <- list()
  etas <- vector("list", n_sub)
  for (id in seq_len(n_sub)) {
    eta <- rnorm(2, sd = sqrt(omega))
    etas[[id]] <- eta
    rows[[length(rows) + 1L]] <- data.frame(
      ID = id, TIME = 0, EVID = 1L, CMT = 1L, AMT = 250, RATE = 0,
      MDV = 1L, DV = 0, SS = 0L, II = 0
    )
    for (tm in times[-1]) {
      rows[[length(rows) + 1L]] <- data.frame(
        ID = id, TIME = tm, EVID = 0L, CMT = 1L, AMT = 0, RATE = 0,
        MDV = 0L, DV = NA_real_, SS = 0L, II = 0
      )
    }
  }
  dat <- do.call(rbind, rows)
  model <- .nm_synthetic_model_from_spec(11L, 4L, dose_cmp = 1L, obs_cmp = 1L)
  list(model = model, data = .nm_synthetic_fill_dv(model, dat, etas, sigma))
}

#' Synthetic oral 3-compartment dataset
#'
#' @inheritParams nm_synthetic_theo
#' @examples
#' nm_synthetic_oral3(n_sub = 2L)
#' @export
nm_synthetic_oral3 <- function(n_sub = 10L, seed = 1L) {
  set.seed(seed)
  omega <- c(0.07, 0.07, 0.07)
  sigma <- c(0.1, 0.4)
  times <- c(0, 0.5, 1, 2, 4, 8, 12, 24, 48)
  rows <- list()
  etas <- vector("list", n_sub)
  for (id in seq_len(n_sub)) {
    eta <- rnorm(3, sd = sqrt(omega))
    etas[[id]] <- eta
    ka <- 0.9 * exp(eta[3])
    rows[[length(rows) + 1L]] <- data.frame(
      ID = id, TIME = 0, EVID = 1L, CMT = 1L, AMT = 200, RATE = 0,
      MDV = 1L, DV = 0, F1 = 1, S1 = 1, KA = ka, SS = 0L, II = 0
    )
    for (tm in times[-1]) {
      rows[[length(rows) + 1L]] <- data.frame(
        ID = id, TIME = tm, EVID = 0L, CMT = 2L, AMT = 0, RATE = 0,
        MDV = 0L, DV = NA_real_, F1 = 1, S1 = 1, KA = ka, SS = 0L, II = 0
      )
    }
  }
  dat <- do.call(rbind, rows)
  model <- .nm_synthetic_model_from_spec(12L, 4L, dose_cmp = 1L, obs_cmp = 2L)
  list(model = model, data = .nm_synthetic_fill_dv(model, dat, etas, sigma))
}

#' Synthetic oral multiple-dose dataset
#'
#' @inheritParams nm_synthetic_theo
#' @examples
#' nm_synthetic_multidose(n_sub = 2L)
#' @export
nm_synthetic_multidose <- function(n_sub = 10L, seed = 1L) {
  set.seed(seed)
  omega <- c(0.05, 0.05)
  sigma <- c(0.1, 0.3)
  dose_times <- c(0, 24)
  obs_times <- c(0.5, 1, 2, 4, 8, 12, 18, 24, 25, 30, 36, 48)
  rows <- list()
  etas <- vector("list", n_sub)
  for (id in seq_len(n_sub)) {
    eta <- rnorm(2, sd = sqrt(omega))
    etas[[id]] <- eta
    ka <- 0.5 * exp(eta[2])
    for (dt in dose_times) {
      rows[[length(rows) + 1L]] <- data.frame(
        ID = id, TIME = dt, EVID = 1L, CMT = 1L, AMT = 10, RATE = 0,
        MDV = 1L, DV = 0, F1 = 1, S1 = 1, KA = ka, SS = 0L, II = 0
      )
    }
    for (tm in obs_times) {
      rows[[length(rows) + 1L]] <- data.frame(
        ID = id, TIME = tm, EVID = 0L, CMT = 2L, AMT = 0, RATE = 0,
        MDV = 0L, DV = NA_real_, F1 = 1, S1 = 1, KA = ka, SS = 0L, II = 0
      )
    }
  }
  dat <- do.call(rbind, rows)
  dat <- dat[order(dat$ID, dat$TIME, -dat$EVID), , drop = FALSE]
  model <- .nm_synthetic_model_from_spec(2L, 1L, dose_cmp = 1L, obs_cmp = 2L)
  list(model = model, data = .nm_synthetic_fill_dv(model, dat, etas, sigma))
}

#' Synthetic oral steady-state dataset
#'
#' @inheritParams nm_synthetic_theo
#' @examples
#' nm_synthetic_steady(n_sub = 2L)
#' @export
nm_synthetic_steady <- function(n_sub = 10L, seed = 1L) {
  set.seed(seed)
  omega <- c(0.05, 0.05)
  sigma <- c(0.1, 0.25)
  dose_ii <- 12
  obs_times <- seq(0.5, 11.5, by = 1)
  rows <- list()
  etas <- vector("list", n_sub)
  for (id in seq_len(n_sub)) {
    eta <- rnorm(2, sd = sqrt(omega))
    etas[[id]] <- eta
    ka <- 0.55 * exp(eta[2])
    rows[[length(rows) + 1L]] <- data.frame(
      ID = id, TIME = 0, EVID = 1L, CMT = 1L, AMT = 8, RATE = 0,
      MDV = 1L, DV = 0, F1 = 1, S1 = 1, KA = ka, SS = 1L, II = dose_ii
    )
    for (tm in obs_times) {
      rows[[length(rows) + 1L]] <- data.frame(
        ID = id, TIME = tm, EVID = 0L, CMT = 2L, AMT = 0, RATE = 0,
        MDV = 0L, DV = NA_real_, F1 = 1, S1 = 1, KA = ka, SS = 0L, II = 0
      )
    }
  }
  dat <- do.call(rbind, rows)
  model <- .nm_synthetic_model_from_spec(
    2L, 1L, dose_cmp = 1L, obs_cmp = 2L, ss = 1L
  )
  list(model = model, data = .nm_synthetic_fill_dv(model, dat, etas, sigma))
}

#' Synthetic IV infusion 1-compartment dataset
#'
#' @inheritParams nm_synthetic_theo
#' @examples
#' nm_synthetic_infusion(n_sub = 2L)
#' @export
nm_synthetic_infusion <- function(n_sub = 10L, seed = 1L) {
  set.seed(seed)
  omega <- c(0.05)
  sigma <- c(0.1)
  inf_rate <- 50
  inf_dur <- 2
  times <- c(0, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 24)
  rows <- list()
  etas <- vector("list", n_sub)
  for (id in seq_len(n_sub)) {
    eta <- rnorm(1, sd = sqrt(omega))
    etas[[id]] <- eta
    rows[[length(rows) + 1L]] <- data.frame(
      ID = id, TIME = 0, EVID = 1L, CMT = 1L,
      AMT = inf_rate * inf_dur, RATE = inf_rate,
      MDV = 1L, DV = 0, SS = 0L, II = 0
    )
    for (tm in times[-1]) {
      rows[[length(rows) + 1L]] <- data.frame(
        ID = id, TIME = tm, EVID = 0L, CMT = 1L, AMT = 0, RATE = 0,
        MDV = 0L, DV = NA_real_, SS = 0L, II = 0
      )
    }
  }
  dat <- do.call(rbind, rows)
  dat <- dat[order(dat$ID, dat$TIME, -dat$EVID), , drop = FALSE]
  model <- .nm_synthetic_model_from_spec(1L, 1L, dose_cmp = 1L, obs_cmp = 1L)
  list(model = model, data = .nm_synthetic_fill_dv(model, dat, etas, sigma))
}

#' Synthetic PK/PD Emax-linked dataset
#'
#' @inheritParams nm_synthetic_theo
#' @examples
#' nm_synthetic_pkpd_emax(n_sub = 2L)
#' @export
nm_synthetic_pkpd_emax <- function(n_sub = 10L, seed = 1L) {
  set.seed(seed)
  omega <- c(0.05, 0.05)
  sigma <- c(0.08, 0.5)
  times <- c(0, 0.5, 1, 2, 4, 8, 12, 24)
  rows <- list()
  etas <- vector("list", n_sub)
  for (id in seq_len(n_sub)) {
    eta <- rnorm(2, sd = sqrt(omega))
    etas[[id]] <- eta
    ka <- 0.7 * exp(eta[2])
    rows[[length(rows) + 1L]] <- data.frame(
      ID = id, TIME = 0, EVID = 1L, CMT = 1L, AMT = 100, RATE = 0,
      MDV = 1L, DV = 0, F1 = 1, S1 = 1, KA = ka, SS = 0L, II = 0
    )
    for (tm in times[-1]) {
      rows[[length(rows) + 1L]] <- data.frame(
        ID = id, TIME = tm, EVID = 0L, CMT = 2L, AMT = 0, RATE = 0,
        MDV = 0L, DV = NA_real_, F1 = 1, S1 = 1, KA = ka, SS = 0L, II = 0
      )
    }
  }
  dat <- do.call(rbind, rows)
  pk_base <- .nm_synthetic_model_from_spec(2L, 1L, dose_cmp = 1L, obs_cmp = 2L)
  cl <- 2.5
  v <- 20
  ka <- 0.7
  model <- nm_model(
    INPUT = pk_base$INPUT,
    ADVAN = pk_base$ADVAN,
    TRANS = pk_base$TRANS,
    DOSECMP = pk_base$DOSECMP,
    OBSCMP = pk_base$OBSCMP,
    PRED = pk_base$PRED,
    ERROR = paste(
      "E0 = THETA(4)",
      "EMAX = THETA(5)",
      "EC50 = THETA(6)",
      "CP = F",
      "Y = E0 + EMAX * CP / (EC50 + CP)",
      sep = "\n"
    ),
    THETAS = data.frame(
      THETA = 1:6,
      Value = c(cl / v, v, ka, 10, 40, 1.5),
      stringsAsFactors = FALSE
    ),
    OMEGAS = data.frame(OMEGA = 1:2, Value = c(0.05, 0.05)),
    SIGMAS = data.frame(SIGMA = 1:2, Value = sigma)
  )
  list(
    model = model,
    data = .nm_synthetic_fill_dv(model, dat, etas, sigma, use_error = TRUE)
  )
}
