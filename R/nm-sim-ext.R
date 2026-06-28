#' Extended simulation helpers
#'
#' @param model An \code{nm_model} object.
#' @param data Template \code{nm_dataset}.
#' @param n_sim Number of simulated datasets.
#' @param seed Random seed.
#' @param n_cores Number of parallel workers for replicates (\code{1L} = sequential).
#' @param pk_engine PK engine passed to \code{.nm_task_sim} (\code{"cpp"} recommended).
#' @param ... Passed to \code{.nm_task_sim}.
#' @return List of simulated \code{nm_dataset} objects, or VPC summary if
#'   \code{vpc = TRUE}.
#' @examples
#' \dontrun{
#' sim <- nm_synthetic_theo(n_sub = 2L)
#' nm_simulate(sim$model, sim$data, n_sim = 2L, seed = 1L)
#' }
#' @export
nm_simulate <- function(model, data, n_sim = 1L, seed = 1L, vpc = FALSE,
                        n_cores = 1L, pk_engine = "cpp", ...) {
  nm_validate_model(model, data = data, stop_on_error = TRUE)
  n_sim <- max(1L, as.integer(n_sim))
  n_cores <- max(1L, as.integer(n_cores))
  sim_args <- list(...)
  sim_one <- function(k) {
    structure(
      list(data = do.call(
        .nm_task_sim,
        c(
          list(
            model = model,
            data = data,
            seed = as.integer(seed) + as.integer(k) - 1L,
            pk_engine = pk_engine,
            rep = as.integer(k)
          ),
          sim_args
        )
      )),
      class = "nm_dataset"
    )
  }
  sims <- if (n_cores > 1L && n_sim > 1L) {
    .nm_parallel_lapply(seq_len(n_sim), sim_one, n_cores = min(n_cores, n_sim))
  } else {
    lapply(seq_len(n_sim), sim_one)
  }
  if (!isTRUE(vpc)) {
    return(sims)
  }
  .nm_vpc_summary(model, data, sims)
}

#' @keywords internal
.nm_parallel_lapply <- function(x, fun, n_cores = 2L) {
  n_cores <- max(1L, as.integer(n_cores))
  if (n_cores <= 1L || length(x) <= 1L) {
    return(lapply(x, fun))
  }
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makePSOCKcluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, c("fun"), envir = environment())
    parallel::clusterEvalQ(cl, {
      if (!requireNamespace("LibeRation", quietly = TRUE)) {
        stop("LibeRation is required for parallel simulation.", call. = FALSE)
      }
      suppressPackageStartupMessages(library(LibeRation))
    })
    parallel::parLapply(cl, x, fun)
  } else {
    parallel::mclapply(x, fun, mc.cores = n_cores)
  }
}

#' @keywords internal
.nm_vpc_compute_breaks <- function(times, n_bins = 10L) {
  times <- as.numeric(times)
  times <- times[is.finite(times)]
  if (length(times) == 0L) {
    return(c(0, 1))
  }
  uniq <- sort(unique(times))
  if (length(uniq) == 1L) {
    u <- uniq[[1L]]
    eps <- max(abs(u) * 1e-8, .Machine$double.eps)
    return(c(u - eps, u + eps))
  }
  n_bins <- max(2L, min(as.integer(n_bins), length(uniq)))
  probs <- seq(0, 1, length.out = n_bins + 1L)
  br <- stats::quantile(uniq, probs = probs, na.rm = TRUE, type = 7L)
  br <- sort(unique(as.numeric(br)))
  if (length(br) < 2L) {
    br <- range(uniq, na.rm = TRUE)
  }
  if (length(br) < 2L || br[[1L]] == br[[length(br)]]) {
    eps <- max(abs(br[[1L]]) * 1e-8, .Machine$double.eps)
    br <- c(br[[1L]] - eps, br[[1L]] + eps)
  }
  br
}

#' @keywords internal
.nm_vpc_assign_bins <- function(times, breaks) {
  times <- as.numeric(times)
  out <- rep(NA_character_, length(times))
  ok <- is.finite(times)
  if (any(ok)) {
    out[ok] <- as.character(base::cut(
      times[ok],
      breaks = breaks,
      include.lowest = TRUE,
      right = TRUE
    ))
  }
  lvls <- sort(unique(out[ok]))
  factor(out, levels = lvls)
}

#' @keywords internal
.nm_vpc_pc_dv <- function(obs, fit) {
  obs <- as.data.frame(obs)
  dv <- obs$DV
  if (is.null(fit)) {
    return(dv)
  }
  pr <- tryCatch(predict(fit, type = "ipred"), error = function(e) NULL)
  if (is.null(pr)) {
    return(dv)
  }
  pr <- as.data.frame(pr)
  pr_obs <- pr[pr$MDV == 0L & pr$EVID == 0L, c("ID", "TIME", "IPRED", "PRED"), drop = FALSE]
  if (nrow(pr_obs) == 0L) {
    return(dv)
  }
  key <- paste(obs$ID, obs$TIME)
  pr_key <- paste(pr_obs$ID, pr_obs$TIME)
  ip <- pr_obs$IPRED[match(key, pr_key)]
  pp <- pr_obs$PRED[match(key, pr_key)]
  dv_pc <- dv
  ok <- is.finite(dv) & is.finite(ip) & ip > 0 & is.finite(pp)
  dv_pc[ok] <- dv[ok] * pp[ok] / ip[ok]
  dv_pc
}

#' @keywords internal
.nm_vpc_aggregate_obs <- function(obs, fit = NULL, pc_correct = FALSE) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    .nm_stop("Package 'data.table' is required for VPC.")
  }
  obs <- as.data.frame(obs)
  obs <- obs[obs$MDV == 0L & obs$EVID == 0L & is.finite(obs$DV), , drop = FALSE]
  if (nrow(obs) == 0L || !"TIME_BIN" %in% names(obs)) {
    return(NULL)
  }
  dv_plot <- if (isTRUE(pc_correct)) {
    .nm_vpc_pc_dv(obs, fit)
  } else {
    obs$DV
  }
  dt <- data.table::as.data.table(obs)
  dt[, DV_plot := dv_plot]
  out <- dt[, .(
    TIME = stats::median(TIME, na.rm = TRUE),
    obs_med = stats::median(DV_plot, na.rm = TRUE),
    obs_lo = as.numeric(stats::quantile(DV_plot, 0.1, na.rm = TRUE)),
    obs_hi = as.numeric(stats::quantile(DV_plot, 0.9, na.rm = TRUE)),
    n_obs = .N
  ), by = TIME_BIN]
  as.data.frame(out)
}

#' @keywords internal
.nm_vpc_sim_pc_rows <- function(model, sim_data, pk_engine = "cpp") {
  full <- .nm_prepare_data(sim_data, model$INPUT, model)
  if (!"IPRED" %in% names(full)) {
    return(NULL)
  }
  theta <- model$THETAS$Value
  omega <- model$OMEGAS$Value
  n_eta <- .nm_n_eta(model)
  ids <- .nm_subject_ids(full)
  parts <- vector("list", length(ids))
  for (j in seq_along(ids)) {
    id <- ids[[j]]
    sub <- .nm_subject_slice(full, id)
    obs_mask <- sub$MDV == 0L & sub$EVID == 0L
    if (!any(obs_mask)) {
      next
    }
    ip <- as.numeric(sub$IPRED[obs_mask])
    pred_pop <- .nm_subject_ipred(
      model, sub, theta, omega, numeric(n_eta), pk_engine = pk_engine
    )
    pp <- as.numeric(pred_pop$F)
    if (length(pp) != sum(obs_mask)) {
      next
    }
    dv <- as.numeric(sub$DV[obs_mask])
    dv_pc <- dv
    ok <- is.finite(dv) & is.finite(ip) & ip > 0 & is.finite(pp)
    dv_pc[ok] <- dv[ok] * pp[ok] / ip[ok]
    parts[[j]] <- data.frame(
      ID = id,
      TIME = sub$TIME[obs_mask],
      DV_pc = dv_pc,
      stringsAsFactors = FALSE
    )
  }
  parts <- parts[!vapply(parts, is.null, logical(1L))]
  if (length(parts) == 0L) {
    return(NULL)
  }
  do.call(rbind, parts)
}

#' @keywords internal
.nm_vpc_sim_quantile_summary <- function(all_sim) {
  per_rep <- all_sim[, .(
    sim_med = stats::median(DV, na.rm = TRUE),
    sim_lo = as.numeric(stats::quantile(DV, 0.1, na.rm = TRUE)),
    sim_hi = as.numeric(stats::quantile(DV, 0.9, na.rm = TRUE))
  ), by = .(sim, TIME_BIN)]
  per_rep[, .(
    sim_med_lo = as.numeric(stats::quantile(sim_med, 0.05, na.rm = TRUE)),
    sim_med_hi = as.numeric(stats::quantile(sim_med, 0.95, na.rm = TRUE)),
    sim_lo_lo = as.numeric(stats::quantile(sim_lo, 0.05, na.rm = TRUE)),
    sim_lo_hi = as.numeric(stats::quantile(sim_lo, 0.95, na.rm = TRUE)),
    sim_hi_lo = as.numeric(stats::quantile(sim_hi, 0.05, na.rm = TRUE)),
    sim_hi_hi = as.numeric(stats::quantile(sim_hi, 0.95, na.rm = TRUE)),
    sim_med = stats::median(sim_med, na.rm = TRUE),
    sim_lo = stats::median(sim_lo, na.rm = TRUE),
    sim_hi = stats::median(sim_hi, na.rm = TRUE),
    n_sim = .N
  ), by = TIME_BIN]
}

#' @keywords internal
.nm_vpc_summary <- function(model, data, sims) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    .nm_stop("Package 'data.table' is required for VPC.")
  }
  obs <- .nm_prepare_data(data, model$INPUT, model)
  obs <- obs[obs$MDV == 0L & obs$EVID == 0L, , drop = FALSE]
  if (nrow(obs) == 0L) {
    .nm_stop("VPC requires at least one observation row (MDV=0, EVID=0).")
  }
  breaks <- .nm_vpc_compute_breaks(obs$TIME, n_bins = 10L)
  obs$TIME_BIN <- .nm_vpc_assign_bins(obs$TIME, breaks)
  sim_dv <- lapply(sims, function(s) {
    d <- s$data
    d <- d[d$MDV == 0L & d$EVID == 0L, c("ID", "TIME", "DV"), drop = FALSE]
    d$TIME_BIN <- .nm_vpc_assign_bins(d$TIME, breaks)
    d
  })
  all_sim <- data.table::rbindlist(sim_dv, idcol = "sim")
  vpc_sim <- .nm_vpc_sim_quantile_summary(all_sim)
  sim_pc <- lapply(sims, function(s) {
    pc <- .nm_vpc_sim_pc_rows(model, s$data)
    if (is.null(pc)) {
      return(NULL)
    }
    pc$TIME_BIN <- .nm_vpc_assign_bins(pc$TIME, breaks)
    pc[, c("ID", "TIME", "TIME_BIN", "DV_pc"), drop = FALSE]
  })
  sim_pc <- sim_pc[!vapply(sim_pc, is.null, logical(1L))]
  if (length(sim_pc) > 0L) {
    all_sim_pc <- data.table::rbindlist(
      lapply(seq_along(sim_pc), function(k) {
        d <- data.table::as.data.table(sim_pc[[k]])
        d[, sim := k]
        d[, DV := DV_pc]
        d[, DV_pc := NULL]
        d
      }),
      use.names = TRUE
    )
    vpc_sim_pc <- .nm_vpc_sim_quantile_summary(all_sim_pc)
    vpc_sim_pc <- as.data.frame(vpc_sim_pc)
    names(vpc_sim_pc) <- c(
      "TIME_BIN",
      "sim_med_lo_pc", "sim_med_hi_pc",
      "sim_lo_lo_pc", "sim_lo_hi_pc",
      "sim_hi_lo_pc", "sim_hi_hi_pc",
      "sim_med_pc", "sim_lo_pc", "sim_hi_pc",
      "n_sim_pc"
    )
    vpc_sim <- merge(vpc_sim, vpc_sim_pc, by = "TIME_BIN", all.x = TRUE)
  }
  obs_sum <- .nm_vpc_aggregate_obs(obs, fit = NULL, pc_correct = FALSE)
  vpc <- merge(as.data.frame(vpc_sim), obs_sum, by = "TIME_BIN", all.x = TRUE, sort = FALSE)
  vpc <- vpc[order(vpc$TIME, na.last = TRUE), , drop = FALSE]
  rownames(vpc) <- NULL
  # Legacy column names for downstream compatibility
  vpc$med <- vpc$sim_med
  vpc$lo <- vpc$sim_lo
  vpc$hi <- vpc$sim_hi
  list(obs = obs, vpc = vpc, sims = sims, breaks = breaks)
}

#' Estimate simulation workload (subjects x observation rows x replicates)
#'
#' @examples
#' sim <- nm_synthetic_theo(n_sub = 2L)
#' nm_sim_workload(sim$model, sim$data, n_sim = 2L)
#' @export
nm_sim_workload <- function(model, data, n_sim = 1L) {
  dat <- .nm_prepare_data(data, model$INPUT, model)
  n_sub <- length(.nm_subject_ids(dat))
  n_obs <- sum(dat$MDV == 0L & dat$EVID == 0L, na.rm = TRUE)
  n_sim <- max(1L, as.integer(n_sim))
  total <- as.numeric(n_sub) * as.numeric(n_obs) * as.numeric(n_sim)
  list(
    n_subjects = n_sub,
    n_obs_rows = as.integer(n_obs),
    n_replicates = n_sim,
    total_points = total
  )
}

#' Combine simulated replicates into one dataset with a REP column
#'
#' @param sims List of \code{nm_dataset} objects from \code{\link{nm_simulate}}.
#' @return A single \code{nm_dataset} with all rows stacked.
#' @keywords internal
.nm_sim_combine_replicates <- function(sims) {
  if (length(sims) == 0L) {
    return(NULL)
  }
  if (length(sims) == 1L) {
    return(sims[[1L]])
  }
  parts <- lapply(seq_along(sims), function(k) {
    d <- sims[[k]]$data
    if (is.null(d)) {
      return(NULL)
    }
    d <- as.data.frame(d)
    if (!"REP" %in% names(d)) {
      d$REP <- as.integer(k)
    }
    d
  })
  parts <- parts[!vapply(parts, is.null, logical(1L))]
  if (length(parts) == 0L) {
    return(NULL)
  }
  structure(
    list(data = do.call(rbind, parts)),
    class = "nm_dataset"
  )
}

#' Package simulation output for workspace persistence
#'
#' @keywords internal
.nm_sim_pack_output <- function(sims) {
  if (is.null(sims)) {
    return(NULL)
  }
  if (is.list(sims) && !is.null(sims$vpc)) {
    combined <- .nm_sim_combine_replicates(sims$sims)
    return(list(
      replicates = sims$sims,
      primary = combined,
      combined = combined,
      vpc = sims$vpc,
      vpc_obs = sims$obs,
      vpc_mode = TRUE
    ))
  }
  if (length(sims) == 0L) {
    return(NULL)
  }
  if (length(sims) == 1L) {
    return(sims[[1L]])
  }
  combined <- .nm_sim_combine_replicates(sims)
  list(replicates = sims, primary = combined, combined = combined)
}

#' Parse dose amount lines (\code{"TIME AMT"} or \code{"AMT"})
#'
#' @keywords internal
.nm_sim_parse_dose_table <- function(text, default_time = 0) {
  if (is.null(text) || !nzchar(trimws(as.character(text)[1L]))) {
    return(data.frame(TIME = default_time, AMT = 320, stringsAsFactors = FALSE))
  }
  lines <- strsplit(trimws(as.character(text)[1L]), "\n", fixed = FALSE)[[1L]]
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  rows <- lapply(lines, function(ln) {
    parts <- strsplit(ln, "[\\s,;]+", perl = TRUE)[[1L]]
    parts <- parts[nzchar(parts)]
    nums <- suppressWarnings(as.numeric(parts))
    nums <- nums[is.finite(nums)]
    if (length(nums) >= 2L) {
      data.frame(TIME = nums[[1L]], AMT = nums[[2L]], stringsAsFactors = FALSE)
    } else if (length(nums) == 1L) {
      data.frame(TIME = default_time, AMT = nums[[1L]], stringsAsFactors = FALSE)
    } else {
      NULL
    }
  })
  rows <- rows[!vapply(rows, is.null, logical(1L))]
  if (length(rows) == 0L) {
    return(data.frame(TIME = default_time, AMT = 320, stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$TIME), , drop = FALSE]
}

#' Resample or expand subjects in a template dataset
#'
#' @keywords internal
.nm_sim_resize_subjects <- function(data, model, n_sub, seed = NULL) {
  n_sub <- max(1L, as.integer(n_sub))
  dat <- .nm_prepare_data(data, model$INPUT, model)
  ids <- .nm_subject_ids(dat)
  if (length(ids) == n_sub) {
    return(nm_dataset_from_table(as.data.frame(dat)))
  }
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }
  if (length(ids) >= n_sub) {
    keep <- sample(ids, n_sub)
    sub <- dat[ID %in% keep]
    return(nm_dataset_from_table(as.data.frame(sub)))
  }
  parts <- lapply(seq_len(n_sub), function(j) {
    src <- ids[((j - 1L) %% length(ids)) + 1L]
    sub <- dat[ID == src]
    sub <- data.table::copy(sub)
    sub[, ID := j]
    sub
  })
  nm_dataset_from_table(as.data.frame(data.table::rbindlist(parts)))
}

#' Build a simulation template dataset from a dosing / sampling design
#'
#' @param model An \code{nm_model} object.
#' @param n_sub Number of subjects.
#' @param n_days Simulation horizon in days (TIME assumed in hours; 24 h per day).
#' @param obs_times Observation times (hours); default from \code{template_data} or a PK grid.
#' @param dose_mode \code{"single"}, \code{"repeat"}, or \code{"steady_state"}.
#' @param dose_amt Scalar dose amount when \code{dose_table} is empty.
#' @param dose_table Optional data frame or parsed text lines with \code{TIME} and \code{AMT}.
#' @param dose_n Number of doses for repeat dosing.
#' @param dose_ii Inter-dose interval (hours).
#' @param dose_cmt Dosing compartment.
#' @param template_data Optional \code{nm_dataset} to inherit observation times / amounts.
#' @examples
#' sim <- nm_synthetic_theo(n_sub = 2L)
#' nm_sim_design_dataset(sim$model, sim$data, n_sub = 2L, seed = 1L)
#' @export
nm_sim_design_dataset <- function(model,
                                  n_sub = 10L,
                                  n_days = 1L,
                                  obs_times = NULL,
                                  obs_per_day = 8L,
                                  dose_mode = c("single", "repeat", "steady_state"),
                                  dose_amt = 320,
                                  dose_table = NULL,
                                  dose_n = 3L,
                                  dose_ii = 12,
                                  dose_cmt = 1L,
                                  template_data = NULL,
                                  seed = NULL) {
  dose_mode <- match.arg(dose_mode)
  n_sub <- max(1L, as.integer(n_sub))
  n_days <- max(1L, as.integer(n_days))
  t_end <- n_days * 24
  dose_cmt <- as.integer(dose_cmt)
  obs_cmt <- as.integer(if (!is.null(model$OBSCMP)) model$OBSCMP else dose_cmt)
  dose_cmt <- as.integer(if (!is.null(model$DOSECMP)) model$DOSECMP else dose_cmt)
  if (is.character(dose_table)) {
    dose_table <- .nm_sim_parse_dose_table(dose_table)
  }
  if (is.null(dose_table) || !is.data.frame(dose_table) || nrow(dose_table) == 0L) {
    dose_table <- data.frame(TIME = 0, AMT = as.numeric(dose_amt), stringsAsFactors = FALSE)
  }
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }
  if (!is.null(template_data)) {
    tdat <- .nm_prepare_data(template_data, model$INPUT, model)
    obs_template <- tdat[tdat$MDV == 0L & tdat$EVID == 0L]
    if (is.null(obs_times) && nrow(obs_template) > 0L) {
      obs_times <- sort(unique(obs_template$TIME))
    }
    dose_template <- tdat[tdat$EVID %in% c(1L, 4L)]
    if (nrow(dose_template) > 0L && nrow(dose_table) == 1L && dose_table$TIME[[1L]] == 0) {
      dose_amt <- dose_template$AMT[[1L]]
      dose_table$AMT[[1L]] <- dose_amt
      if (dose_mode == "repeat" && "II" %in% names(dose_template)) {
        dose_ii <- dose_template$II[[1L]]
      }
    }
  }
  if (is.null(obs_times) || length(obs_times) == 0L) {
    obs_per_day <- max(3L, as.integer(obs_per_day))
    obs_times <- unique(c(
      seq(0.25, min(24, t_end), length.out = obs_per_day),
      if (t_end > 24) seq(24, t_end, length.out = min(obs_per_day, 8L)) else numeric()
    ))
  }
  obs_times <- sort(unique(as.numeric(obs_times)))
  obs_times <- obs_times[is.finite(obs_times) & obs_times > 0 & obs_times <= t_end]

  build_subject <- function(id) {
    rows <- list()
    mk_row <- function(time, evid, cmt, amt, mdv, dv = NA_real_, extra = list()) {
      row <- c(
        list(
          ID = id, TIME = time, EVID = evid, CMT = cmt,
          AMT = amt, MDV = mdv, DV = dv,
          RATE = 0, F1 = 1, SS = 0L, II = 0, ADDL = 0L
        ),
        extra
      )
      as.data.frame(row, stringsAsFactors = FALSE)
    }
    if (dose_mode == "single") {
      for (i in seq_len(nrow(dose_table))) {
        rows[[length(rows) + 1L]] <- mk_row(
          dose_table$TIME[[i]], 1L, dose_cmt, dose_table$AMT[[i]], 1L, 0
        )
      }
    } else if (dose_mode == "repeat") {
      dose_n <- max(1L, as.integer(dose_n))
      dose_ii <- as.numeric(dose_ii)
      base_amt <- dose_table$AMT[[1L]]
      for (k in seq_len(dose_n)) {
        t_dose <- (k - 1L) * dose_ii
        if (t_dose > t_end) {
          break
        }
        rows[[length(rows) + 1L]] <- mk_row(
          t_dose, 1L, dose_cmt, base_amt, 1L, 0,
          extra = list(II = dose_ii, ADDL = 0L, SS = 0L)
        )
      }
    } else {
      base_amt <- dose_table$AMT[[1L]]
      rows[[1L]] <- mk_row(
        0, 1L, dose_cmt, base_amt, 1L, 0,
        extra = list(II = as.numeric(dose_ii), ADDL = 0L, SS = 1L)
      )
    }
    for (tm in obs_times) {
      rows[[length(rows) + 1L]] <- mk_row(
        tm, 0L, obs_cmt, 0, 0L, NA_real_
      )
    }
    out <- if (requireNamespace("data.table", quietly = TRUE)) {
      as.data.frame(data.table::rbindlist(rows, fill = TRUE, use.names = TRUE))
    } else {
      do.call(rbind, rows)
    }
    out[order(out$TIME, -out$EVID), , drop = FALSE]
  }
  dat <- do.call(rbind, lapply(seq_len(n_sub), build_subject))
  rownames(dat) <- NULL
  nm_dataset_from_table(dat)
}

#' Resolve simulation template data from design options
#'
#' @param model Model object.
#' @param data Linked version dataset.
#' @param design List with design fields from the Shiny simulation dialog.
#' @examples
#' sim <- nm_synthetic_theo(n_sub = 2L)
#' nm_sim_template_data(sim$model, sim$data)
#' @export
nm_sim_template_data <- function(model, data, design = list()) {
  use_design <- isTRUE(design$use_design)
  n_sub <- max(1L, as.integer(if (is.null(design$n_sub)) 1L else design$n_sub))
  seed <- design$seed
  if (isTRUE(use_design)) {
    return(nm_sim_design_dataset(
      model = model,
      n_sub = n_sub,
      n_days = as.integer(if (is.null(design$n_days)) 1L else design$n_days),
      obs_times = design$obs_times,
      obs_per_day = as.integer(if (is.null(design$obs_per_day)) 8L else design$obs_per_day),
      dose_mode = if (is.null(design$dose_mode)) "single" else design$dose_mode,
      dose_amt = as.numeric(if (is.null(design$dose_amt)) 320 else design$dose_amt),
      dose_table = design$dose_table,
      dose_n = as.integer(if (is.null(design$dose_n)) 1L else design$dose_n),
      dose_ii = as.numeric(if (is.null(design$dose_ii)) 12 else design$dose_ii),
      dose_cmt = as.integer(if (is.null(design$dose_cmt)) 1L else design$dose_cmt),
      template_data = data,
      seed = seed
    ))
  }
  .nm_sim_resize_subjects(data, model, n_sub = n_sub, seed = seed)
}
