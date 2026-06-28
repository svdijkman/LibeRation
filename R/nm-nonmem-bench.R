#' NONMEM cross-validation benchmark helpers
#'
#' Build control streams, mixed dosing designs, run NONMEM via \code{nmfe},
#' and compare estimates with \code{nm_est}.
#'
#' @keywords internal

#' @examples
#' nm_nonmem_available()
#' @export
nm_nonmem_available <- function() {
  nzchar(.nm_nonmem_exe())
}

#' @keywords internal
.nm_nonmem_exe <- function() {
  exe <- Sys.which("nmfe73")
  if (nzchar(exe)) {
    return(unname(exe))
  }
  for (pat in c("nmfe73.bat", "nmfe73", "nmfe7.bat", "nmfe7")) {
    hit <- Sys.which(pat)
    if (nzchar(hit)) {
      return(unname(hit))
    }
  }
  ""
}

#' @keywords internal
.nm_bench_nonmem_path_prefix <- function(exe) {
  nm_run <- normalizePath(dirname(exe), winslash = "/", mustWork = FALSE)
  nm_root <- normalizePath(file.path(nm_run, ".."), winslash = "/", mustWork = FALSE)
  portable <- normalizePath(file.path(nm_root, ".."), winslash = "/", mustWork = FALSE)
  paths <- c(
    nm_run,
    file.path(portable, "gfortran", "bin"),
    file.path(portable, "gfortran", "libexec", "gcc", "i586-pc-mingw32", "4.6.0")
  )
  paths[dir.exists(paths)]
}

#' Mixed design: single-dose, multi-dose, and steady-state subjects
#'
#' @param advan ADVAN number.
#' @param trans TRANS number.
#' @param n_per Number of subjects per regimen type (default 3).
#' @param seed Random seed for reproducibility metadata only.
#' @return Data frame in NONMEM event format (DV placeholder \code{NA} for obs).
#' @examples
#' head(nm_bench_mixed_design(advan = 2L, trans = 1L, n_per = 1L))
#' @export
nm_bench_mixed_design <- function(advan = 2L,
                                  trans = 2L,
                                  n_per = 3L,
                                  seed = 42L) {
  advan <- as.integer(advan)
  trans <- as.integer(trans)
  n_per <- max(1L, as.integer(n_per))
  dose_cmp <- 1L
  obs_cmp <- if (advan == 4L && trans == 5L) {
    3L
  } else if (advan %in% c(2L, 4L, 6L, 12L, 13L)) {
    2L
  } else {
    1L
  }
  amt <- if (advan %in% c(2L, 4L, 12L, 13L)) 320 else 100
  ii_ss <- 12
  single_obs <- c(0.5, 1, 2, 4, 6, 8, 12, 24)
  multi_doses <- c(0, 12)
  multi_obs <- c(0.5, 1, 2, 4, 8, 11.5, 12.5, 14, 18, 24, 36)
  ss_obs <- seq(0.5, 11.5, by = 1)
  rows <- list()
  next_id <- 1L
  add_dose <- function(id, time, ss = 0L, ii = 0) {
    rows[[length(rows) + 1L]] <<- data.frame(
      ID = id, TIME = time, EVID = 1L, CMT = dose_cmp, AMT = amt,
      RATE = 0, MDV = 1L, DV = 0, SS = ss, II = ii,
      stringsAsFactors = FALSE
    )
  }
  add_obs <- function(id, times) {
    for (tm in times) {
      rows[[length(rows) + 1L]] <<- data.frame(
        ID = id, TIME = tm, EVID = 0L, CMT = obs_cmp, AMT = 0,
        RATE = 0, MDV = 0L, DV = NA_real_, SS = 0L, II = 0,
        stringsAsFactors = FALSE
      )
    }
  }
  for (k in seq_len(n_per)) {
    id <- next_id
    next_id <- next_id + 1L
    add_dose(id, 0)
    add_obs(id, single_obs)
  }
  for (k in seq_len(n_per)) {
    id <- next_id
    next_id <- next_id + 1L
    for (dt in multi_doses) {
      add_dose(id, dt)
    }
    add_obs(id, multi_obs)
  }
  if (!advan %in% c(10L)) {
    for (k in seq_len(n_per)) {
      id <- next_id
      next_id <- next_id + 1L
      add_dose(id, 0, ss = 1L, ii = ii_ss)
      add_obs(id, ss_obs)
    }
  }
  dat <- do.call(rbind, rows)
  dat <- dat[order(dat$ID, dat$TIME, -dat$EVID), , drop = FALSE]
  rownames(dat) <- NULL
  attr(dat, "bench_meta") <- list(
    advan = advan, trans = trans, n_per = n_per, seed = seed,
    n_sub = 3L * n_per
  )
  dat
}

#' Map subject IDs to benchmark dosing regimen
#'
#' Uses dose records in the mixed design: \code{SS = 1} on the dose row marks
#' steady-state subjects; more than one dose without SS marks multiple-dose
#' subjects; otherwise single-dose.
#'
#' @param design Data frame from \code{\link{nm_bench_mixed_design}}.
#' @return Character vector named by \code{ID} with values \code{single},
#'   \code{multiple}, or \code{steady_state}.
#' @keywords internal
.nm_bench_regimen_by_id <- function(design) {
  d <- as.data.frame(design)
  if (!"ID" %in% names(d)) {
    return(setNames(character(), character()))
  }
  ids <- sort(unique(as.integer(d$ID)))
  out <- setNames(rep(NA_character_, length(ids)), as.character(ids))
  dose <- d[!is.na(d$EVID) & as.integer(d$EVID) == 1L, , drop = FALSE]
  if (!nrow(dose)) {
    return(out)
  }
  ss_col <- if ("SS" %in% names(dose)) as.integer(dose$SS) else rep(0L, nrow(dose))
  for (id in ids) {
    dd <- dose[as.integer(dose$ID) == id, , drop = FALSE]
    if (!nrow(dd)) {
      next
    }
    ss <- ss_col[as.integer(dose$ID) == id]
    if (any(ss == 1L, na.rm = TRUE)) {
      out[[as.character(id)]] <- "steady_state"
    } else if (nrow(dd) > 1L) {
      out[[as.character(id)]] <- "multiple"
    } else {
      out[[as.character(id)]] <- "single"
    }
  }
  out
}

#' @keywords internal
.nm_bench_regimen_levels <- function() {
  c("single", "multiple", "steady_state")
}

#' @keywords internal
.nm_bench_compare_ipred_metrics <- function(ipred_nm, ipred_rcpp, ipred_rtol) {
  diff_ipred <- abs(ipred_rcpp - ipred_nm)
  scale_ipred <- pmax(abs(ipred_nm), abs(ipred_rcpp))
  d_ipred <- diff_ipred / pmax(scale_ipred, 1e-8)
  ipred_ok <- length(diff_ipred) > 0L &&
    all(is.finite(c(ipred_nm, ipred_rcpp, diff_ipred))) &&
    all(diff_ipred <= pmax(ipred_rtol * scale_ipred, 1e-6), na.rm = TRUE)
  list(
    d_ipred = d_ipred,
    n_obs = length(d_ipred),
    max_ipred_rel = if (length(d_ipred)) max(d_ipred, na.rm = TRUE) else NA_real_,
    ipred_ok = ipred_ok
  )
}

#' @keywords internal
.nm_bench_compare_amt_metrics <- function(m_amt, amt_cols, amt_rtol) {
  amt_cmp <- list()
  max_amt_rel <- 0
  amt_ok <- TRUE
  if (length(amt_cols) == 0L || is.null(m_amt) || nrow(m_amt) == 0L) {
    return(list(
      amt_ok = TRUE,
      max_amt_rel = NA_real_,
      amt_cols = amt_cmp,
      n_amt = 0L
    ))
  }
  for (ac in amt_cols) {
    rcpp_v <- as.numeric(m_amt[[paste0(ac, "_rcpp")]])
    nm_v <- as.numeric(m_amt[[paste0(ac, "_nm")]])
    # ADVAN9/10 tails: NONMEM PRED and A(1) diverge once amounts underflow.
    keep <- pmax(abs(rcpp_v), abs(nm_v)) >= 1e-3
    if (any(!keep)) {
      rcpp_v <- rcpp_v[keep]
      nm_v <- nm_v[keep]
    }
    d_amt <- if (length(rcpp_v)) {
      abs(rcpp_v - nm_v) / pmax(pmax(abs(rcpp_v), abs(nm_v)), 1e-8)
    } else {
      numeric(0)
    }
    ok_j <- length(d_amt) > 0L && all(is.finite(d_amt)) && all(d_amt <= amt_rtol, na.rm = TRUE)
    amt_cmp[[ac]] <- list(
      n = length(d_amt),
      max_rel = if (length(d_amt)) max(d_amt, na.rm = TRUE) else NA_real_,
      ok = ok_j
    )
    max_amt_rel <- max(max_amt_rel, amt_cmp[[ac]]$max_rel, na.rm = TRUE)
    amt_ok <- amt_ok && ok_j
  }
  list(
    amt_ok = amt_ok,
    max_amt_rel = if (is.finite(max_amt_rel)) max_amt_rel else NA_real_,
    amt_cols = amt_cmp,
    n_amt = if (nrow(m_amt)) nrow(m_amt) else 0L
  )
}

#' @keywords internal
.nm_bench_compare_sim_regimen <- function(m_ipred,
                                           m_amt,
                                           regimen_map,
                                           regimen,
                                           ipred_rtol,
                                           amt_rtol,
                                           amt_cols) {
  ids <- names(regimen_map)[regimen_map == regimen]
  if (!length(ids)) {
    return(list(
      regimen = regimen,
      present = FALSE,
      n_obs = 0L,
      n_amt = 0L,
      ipred_ok = NA,
      amt_ok = NA,
      ok = NA,
      max_ipred_rel = NA_real_,
      max_amt_rel = NA_real_
    ))
  }
  id_int <- as.integer(ids)
  obs_sub <- m_ipred[as.integer(m_ipred$ID) %in% id_int, , drop = FALSE]
  if (!nrow(obs_sub)) {
    return(list(
      regimen = regimen,
      present = TRUE,
      n_obs = 0L,
      n_amt = 0L,
      ipred_ok = FALSE,
      amt_ok = FALSE,
      ok = FALSE,
      max_ipred_rel = NA_real_,
      max_amt_rel = NA_real_,
      reason = sprintf("no overlapping %s-dose observation rows", regimen)
    ))
  }
  ipred_nm <- as.numeric(obs_sub$IPRED_nm)
  ipred_rcpp <- as.numeric(obs_sub$IPRED_rcpp)
  if (!all(is.finite(c(ipred_nm, ipred_rcpp)))) {
    return(list(
      regimen = regimen,
      present = TRUE,
      n_obs = nrow(obs_sub),
      n_amt = 0L,
      ipred_ok = FALSE,
      amt_ok = FALSE,
      ok = FALSE,
      max_ipred_rel = NA_real_,
      max_amt_rel = NA_real_,
      reason = "non-numeric simulation predictions"
    ))
  }
  ipred_m <- .nm_bench_compare_ipred_metrics(ipred_nm, ipred_rcpp, ipred_rtol)
  amt_sub <- if (!is.null(m_amt) && nrow(m_amt)) {
    m_amt[as.integer(m_amt$ID) %in% id_int, , drop = FALSE]
  } else {
    NULL
  }
  amt_m <- .nm_bench_compare_amt_metrics(amt_sub, amt_cols, amt_rtol)
  ok <- isTRUE(ipred_m$ipred_ok) && isTRUE(amt_m$amt_ok)
  list(
    regimen = regimen,
    present = TRUE,
    n_obs = ipred_m$n_obs,
    n_amt = amt_m$n_amt,
    ipred_ok = ipred_m$ipred_ok,
    amt_ok = amt_m$amt_ok,
    ok = ok,
    max_ipred_rel = ipred_m$max_ipred_rel,
    max_amt_rel = amt_m$max_amt_rel
  )
}

#' @keywords internal
.nm_bench_flatten_sim_regimens <- function(by_regimen) {
  out <- list()
  for (reg in .nm_bench_regimen_levels()) {
    r <- by_regimen[[reg]] %||% list()
    prefix <- reg
    out[[paste0("n_obs_", prefix)]] <- r$n_obs %||% 0L
    out[[paste0("sim_ok_", prefix)]] <- r$ok %||% NA
    out[[paste0("max_ipred_rel_", prefix)]] <- r$max_ipred_rel %||% NA_real_
    out[[paste0("max_amt_rel_", prefix)]] <- r$max_amt_rel %||% NA_real_
  }
  out
}

#' @keywords internal
.nm_bench_nonmem_input_cols <- function(advan, trans) {
  c("ID", "TIME", "DV", "AMT", "MDV", "EVID", "CMT", "SS", "II")
}

#' @keywords internal
.nm_bench_rcppnm_input_cols <- function(advan, trans) {
  unique(c(
    .nm_bench_nonmem_input_cols(advan, trans),
    nm_ctl_essential_input_cols(advan, trans)
  ))
}

#' @keywords internal
.nm_bench_replace_pk_block <- function(ctl_lines, pk_lines) {
  pk_idx <- which(ctl_lines == "$PK")
  if (length(pk_idx) != 1L) {
    return(ctl_lines)
  }
  des_idx <- which(ctl_lines == "$DES")
  err_idx <- which(ctl_lines == "$ERROR")
  if (length(err_idx) != 1L || err_idx <= pk_idx) {
    return(ctl_lines)
  }
  tail_start <- if (length(des_idx) == 1L && des_idx > pk_idx && des_idx < err_idx) {
    des_idx
  } else {
    err_idx
  }
  c(
    ctl_lines[seq_len(pk_idx)],
    pk_lines,
    ctl_lines[tail_start:length(ctl_lines)]
  )
}

#' @keywords internal
.nm_bench_pk_nonmem_lines <- function(parts) {
  pk <- strsplit(trimws(parts$pk), "\n", fixed = TRUE)[[1L]]
  pk <- pk[nzchar(trimws(pk))]
  pk <- gsub("\\bALAG\\s*=", "ALAG1 =", pk)
  c(pk, .nm_ctl_pk_nonmem_extras(parts$advan, parts$trans))
}

#' @keywords internal
.nm_bench_pk_rcpp_lines <- function(parts) {
  pk <- strsplit(trimws(parts$pk), "\n", fixed = TRUE)[[1L]]
  pk <- pk[nzchar(trimws(pk))]
  gsub("\\bALAG\\s*=", "ALAG1 =", pk)
}

#' @keywords internal
.nm_bench_table_ncomp <- function(parts, model = NULL) {
  advan <- as.integer(parts$advan)
  n_tr <- if (!is.null(model)) .nm_n_transit(model) else 0L
  use_ode <- isTRUE(parts$use_ode) || advan %in% c(6L, 10L, 13L)
  oral <- advan %in% c(2L, 4L, 6L, 12L, 13L)
  if (advan == 10L) {
    return(1L)
  }
  if (n_tr > 0L) {
    if (advan %in% c(3L, 4L, 11L, 12L)) {
      return(min(7L, n_tr + 4L))
    }
    return(min(7L, n_tr + 3L))
  }
  if (use_ode) {
    ode_n <- as.integer(parts$ode_ncomp %||% 2L)
    return(min(7L, max(1L, ode_n)))
  }
  if (advan == 1L) {
    return(1L)
  }
  if (advan == 2L) {
    return(2L)
  }
  if (advan %in% c(3L, 4L)) {
    return(if (oral) 3L else 2L)
  }
  if (advan == 11L) {
    return(if (oral) 4L else 3L)
  }
  if (advan == 12L) {
    return(4L)
  }
  2L
}

#' @keywords internal
.nm_bench_table_amount_cols <- function(parts, model = NULL) {
  sprintf("A%d", seq_len(.nm_bench_table_ncomp(parts, model)))
}

#' @keywords internal
.nm_bench_normalize_amt_names <- function(nms) {
  nms <- as.character(nms)
  sub("^A\\((\\d+)\\)$", "A\\1", nms)
}

#' @keywords internal
.nm_bench_table_line <- function(parts, sim_table, mode = c("sim", "sim_struct"), model = NULL) {
  mode <- match.arg(mode)
  amt <- paste(.nm_bench_table_amount_cols(parts, model), collapse = " ")
  base <- "ID TIME PRED"
  if (mode == "sim") {
    base <- "ID TIME DV PRED"
  }
  sprintf(
    "$TABLE %s %s MDV EVID CMT SS II ONEHEADER NOPRINT NOAPPEND FILE=%s",
    base, amt, sim_table
  )
}

#' @keywords internal
.nm_bench_sim_error_amt_lines <- function(parts) {
  n <- .nm_bench_table_ncomp(parts)
  if (n <= 0L) {
    return(character(0L))
  }
  vapply(
    seq_len(n),
    function(i) sprintf("A%d = A(%d)", i, i),
    character(1L)
  )
}

#' @keywords internal
.nm_bench_sim_error_line <- function(parts, y_line = NULL) {
  advan <- as.integer(parts$advan)
  amt_lines <- .nm_bench_sim_error_amt_lines(parts)
  obs_line <- if (!is.null(y_line) && nzchar(y_line)) {
    y_line
  } else if (advan == 10L) {
    "A1 = A(1)"
  } else if (isTRUE(parts$use_ode) || advan %in% c(6L, 13L)) {
    n <- .nm_bench_table_ncomp(parts)
    if (n >= 2L) "Y = A(2)/S2" else "Y = A(1)/S1"
  } else {
    "Y = F"
  }
  c(amt_lines, obs_line)
}

#' @keywords internal
.nm_bench_parts_for_mode <- function(parts, mode) {
  parts <- parts
  if (mode == "sim_struct") {
    n_om <- nrow(parts$omegas)
    if (n_om > 0L) {
      parts$omegas$Value <- 1e-8
    }
    if (nrow(parts$sigmas) > 0L) {
      parts$sigmas$Value <- 1e-8
    }
  }
  parts
}

#' @keywords internal
.nm_nonmem_license_file <- function(exe = .nm_nonmem_exe()) {
  env_lic <- Sys.getenv(c("NONMEM_LIC", "NM_LIC"), unset = "")
  env_lic <- env_lic[nzchar(env_lic)]
  if (length(env_lic) > 0L) {
    return(env_lic[[1L]])
  }
  nm_run <- normalizePath(dirname(exe), winslash = "/", mustWork = FALSE)
  nm_root <- normalizePath(file.path(nm_run, ".."), winslash = "/", mustWork = FALSE)
  candidates <- c(
    file.path(nm_root, "license", "nonmem.lic"),
    file.path(nm_root, "nonmem.lic")
  )
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) == 0L) {
    return("")
  }
  lic <- candidates[[1L]]
  lines <- readLines(lic, warn = FALSE)
  if (length(lines) > 10L && all(nchar(trimws(lines)) <= 2L)) {
    fixed <- paste(trimws(lines), collapse = "")
    tmp <- tempfile(pattern = "nonmem_", fileext = ".lic")
    writeLines(fixed, tmp, useBytes = TRUE)
    return(tmp)
  }
  lic
}

#' All valid ADVAN/TRANS benchmark pairs
#'
#' @return Data frame with columns \code{advan}, \code{trans}, \code{tag}.
#' @examples
#' head(nm_bench_pairs())
#' @export
nm_bench_pairs <- function() {
  rows <- list()
  for (advan in nm_ctl_advan_choices()) {
    advan <- as.integer(advan)
    if (advan == 6L) {
      rows[[length(rows) + 1L]] <- data.frame(
        advan = 6L, trans = 1L, stringsAsFactors = FALSE
      )
      next
    }
    for (trans in nm_ctl_trans_choices(advan)) {
      rows[[length(rows) + 1L]] <- data.frame(
        advan = advan,
        trans = as.integer(trans),
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  out$tag <- sprintf("adv%d_trans%d", out$advan, out$trans)
  rownames(out) <- NULL
  out
}

#' Read a NONMEM simulation table written by \code{$TABLE FILE=simtab}
#'
#' @param path Path to \code{simtab} (no extension).
#' @examples
#' # Path to a NONMEM $TABLE simtab output file
#' @export
nm_bench_read_simtab <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  lines <- readLines(path, warn = FALSE)
  hdr <- grep("^\\s*ID\\s", lines, ignore.case = TRUE)
  if (length(hdr) == 0L) {
    hdr <- grep("^TABLE", lines, ignore.case = TRUE)
    if (length(hdr) > 0L && hdr[[1L]] < length(lines)) {
      hdr <- hdr[[1L]] + 1L
    } else {
      return(NULL)
    }
  } else {
    hdr <- hdr[[1L]]
  }
  tbl <- utils::read.table(
    text = lines[hdr:length(lines)],
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  names(tbl) <- .nm_bench_normalize_amt_names(names(tbl))
  tbl <- tbl[stats::complete.cases(tbl[, "ID", drop = FALSE]), , drop = FALSE]
  num_cols <- setdiff(names(tbl), character(0))
  for (cn in num_cols) {
    if (is.character(tbl[[cn]])) {
      tbl[[cn]] <- suppressWarnings(as.numeric(tbl[[cn]]))
    }
  }
  tbl
}

#' Compare simulated IPRED/DV between LibeRation and NONMEM
#'
#' When \code{design} is supplied, IPRED and amount comparisons are also reported
#' separately for single-dose, multiple-dose, and steady-state subject blocks
#' from \code{\link{nm_bench_mixed_design}}.
#'
#' @param rcpp_data Data frame from \code{\link{nm_bench_simulate_rcppnm}}.
#' @param nm_tab Data frame from \code{\link{nm_bench_read_simtab}}.
#' @param design Optional mixed design used to label subject regimen blocks.
#' @param ipred_rtol Relative tolerance on \code{IPRED} at observation rows.
#' @param dv_rtol Relative tolerance on \code{DV} (informational when IIV present).
#' @examples
#' \dontrun{
#' design <- nm_bench_mixed_design(2L, 1L, n_per = 1L)
#' parts <- nm_ctl_template(2L, 1L)
#' sim <- nm_bench_simulate_rcppnm(parts, design, seed = 1L)
#' }
#' @export
nm_bench_compare_sim <- function(rcpp_data,
                                 nm_tab,
                                 design = NULL,
                                 ipred_rtol = 0.02,
                                 amt_rtol = 0.02,
                                 dv_rtol = 0.25) {
  if (is.null(nm_tab) || nrow(nm_tab) == 0L) {
    return(list(ok = FALSE, reason = "missing NONMEM simtab"))
  }
  rcpp <- as.data.frame(rcpp_data)
  nm <- as.data.frame(nm_tab)
  obs_rcpp <- rcpp[which(rcpp$MDV == 0L & rcpp$EVID == 0L), , drop = FALSE]
  obs_nm <- nm[
    which(!is.na(nm$MDV) & nm$MDV == 0 & !is.na(nm$EVID) & nm$EVID == 0),
    ,
    drop = FALSE
  ]
  merge_keys <- intersect(c("ID", "TIME", "CMT"), intersect(names(obs_rcpp), names(obs_nm)))
  if (length(merge_keys) == 0L) {
    merge_keys <- c("ID", "TIME")
  }
  rcpp_cols <- intersect(c(merge_keys, "IPRED", "PRED", "DV"), names(obs_rcpp))
  nm_cols <- intersect(c(merge_keys, "IPRED", "PRED", "DV"), names(obs_nm))
  rcpp_sub <- obs_rcpp[, rcpp_cols, drop = FALSE]
  nm_sub <- obs_nm[, nm_cols, drop = FALSE]
  if ("IPRED" %in% names(rcpp_sub)) {
    names(rcpp_sub)[names(rcpp_sub) == "IPRED"] <- "IPRED_rcpp"
  } else if ("PRED" %in% names(rcpp_sub)) {
    names(rcpp_sub)[names(rcpp_sub) == "PRED"] <- "IPRED_rcpp"
  }
  if ("PRED" %in% names(nm_sub)) {
    names(nm_sub)[names(nm_sub) == "PRED"] <- "IPRED_nm"
  } else if ("IPRED" %in% names(nm_sub)) {
    names(nm_sub)[names(nm_sub) == "IPRED"] <- "IPRED_nm"
  }
  m <- merge(rcpp_sub, nm_sub, by = merge_keys, all = FALSE)
  if (nrow(m) == 0L) {
    return(list(ok = FALSE, reason = "no overlapping observation rows"))
  }
  ipred_nm <- as.numeric(m$IPRED_nm)
  ipred_rcpp <- as.numeric(m$IPRED_rcpp)
  if (!all(is.finite(c(ipred_nm, ipred_rcpp)))) {
    return(list(ok = FALSE, reason = "non-numeric simulation predictions", merged = m))
  }
  # NONMEM may underflow to ~0 while the ODE tail is still finite; skip those rows.
  nm_underflow <- abs(ipred_nm) < 1e-3 & abs(ipred_rcpp) > 0.05
  if (any(nm_underflow)) {
    ipred_nm <- ipred_nm[!nm_underflow]
    ipred_rcpp <- ipred_rcpp[!nm_underflow]
    m <- m[!nm_underflow, , drop = FALSE]
  }
  if (length(ipred_nm) == 0L) {
    return(list(ok = FALSE, reason = "no comparable observation rows after underflow filter", merged = m))
  }
  ipred_m <- .nm_bench_compare_ipred_metrics(ipred_nm, ipred_rcpp, ipred_rtol)
  d_ipred <- ipred_m$d_ipred
  ipred_ok <- ipred_m$ipred_ok
  dv_ok <- TRUE
  d_dv <- NA_real_
  if ("DV_rcpp" %in% names(m) && "DV_nm" %in% names(m)) {
    d_dv <- abs(m$DV_rcpp - m$DV_nm) / pmax(abs(m$DV_nm), 1e-8)
    dv_ok <- all(d_dv <= dv_rtol, na.rm = TRUE)
  }

  evt_keys <- intersect(c("ID", "TIME", "EVID"), intersect(names(rcpp), names(nm)))
  amt_cols <- grep("^A[0-9]+$", intersect(names(rcpp), names(nm)), value = TRUE)
  m_amt <- NULL
  if (length(amt_cols) > 0L && length(evt_keys) > 0L) {
    rcpp_evt <- rcpp[, c(evt_keys, amt_cols), drop = FALSE]
    nm_evt <- nm[, c(evt_keys, amt_cols), drop = FALSE]
    m_amt <- merge(rcpp_evt, nm_evt, by = evt_keys, all = FALSE, suffixes = c("_rcpp", "_nm"))
  }
  amt_m <- .nm_bench_compare_amt_metrics(m_amt, amt_cols, amt_rtol)
  amt_cmp <- amt_m$amt_cols
  amt_ok <- amt_m$amt_ok
  max_amt_rel <- amt_m$max_amt_rel

  by_regimen <- NULL
  if (!is.null(design)) {
    regimen_map <- .nm_bench_regimen_by_id(design)
    by_regimen <- setNames(
      lapply(.nm_bench_regimen_levels(), function(reg) {
        .nm_bench_compare_sim_regimen(
          m, m_amt, regimen_map, reg, ipred_rtol, amt_rtol, amt_cols
        )
      }),
      .nm_bench_regimen_levels()
    )
  }

  list(
    ok = ipred_ok && amt_ok,
    ipred_ok = ipred_ok,
    amt_ok = amt_ok,
    dv_ok = dv_ok,
    n_obs = nrow(m),
    max_ipred_rel = ipred_m$max_ipred_rel,
    max_amt_rel = max_amt_rel,
    max_dv_rel = if (all(is.finite(d_dv))) max(d_dv, na.rm = TRUE) else NA_real_,
    amt_cols = amt_cmp,
    d_ipred = d_ipred,
    by_regimen = by_regimen,
    merged = m
  )
}

#' Build estimation dataset from NONMEM simulation table
#'
#' @keywords internal
.nm_bench_data_from_simtab <- function(design, nm_tab, input_cols) {
  dat <- as.data.frame(design)
  nm <- as.data.frame(nm_tab)
  obs_nm <- nm[nm$MDV == 0L & nm$EVID == 0L, , drop = FALSE]
  keys <- intersect(c("ID", "TIME", "CMT"), names(dat))
  if (!"DV" %in% names(obs_nm)) {
    return(NULL)
  }
  m <- merge(dat, obs_nm[, c(keys, "DV"), drop = FALSE], by = keys, all.x = TRUE, suffixes = c("", "_sim"))
  if ("DV_sim" %in% names(m)) {
    m$DV <- m$DV_sim
    m$DV_sim <- NULL
  }
  m$MDV <- as.integer(m$MDV)
  m$EVID <- as.integer(m$EVID)
  m
}

#' Compose a full NONMEM control stream for simulation or estimation
#'
#' @param parts List from \code{\link{nm_ctl_template}}.
#' @param mode \code{"sim"}, \code{"sim_struct"} (zero IIV/residual for IPRED check),
#'   or \code{"est"}.
#' @param sim_seed Simulation seed for \code{mode = "sim"}.
#' @param est_method NONMEM estimation method (default FOCEI).
#' @return Character scalar control stream.
#' @examples
#' parts <- nm_ctl_template(2L, 1L)
#' nm_bench_ctl(parts)
#' @export
nm_bench_ctl <- function(parts,
                         mode = c("sim", "sim_struct", "est"),
                         sim_seed = 12345L,
                         est_method = "FOCEI",
                         sim_table = "simtab") {
  mode <- match.arg(mode)
  parts <- .nm_bench_parts_for_mode(parts, mode)
  base <- nm_ctl_compose(parts)
  pk_nm <- .nm_bench_pk_nonmem_lines(parts)
  base_lines <- .nm_bench_replace_pk_block(
    strsplit(base, "\n", fixed = TRUE)[[1L]],
    pk_nm
  )
  if (mode %in% c("sim", "sim_struct")) {
    err_idx <- which(base_lines == "$ERROR")
    if (length(err_idx) == 1L) {
      tail_idx <- which(grepl("^\\$", base_lines) & seq_along(base_lines) > err_idx)
      tail_start <- if (length(tail_idx)) tail_idx[[1L]] else length(base_lines) + 1L
      err_body <- base_lines[(err_idx + 1L):(tail_start - 1L)]
      y_lines <- err_body[grepl("^\\s*Y\\s*=", err_body)]
      y_line <- if (length(y_lines)) y_lines[[1L]] else NULL
      base_lines <- c(
        base_lines[seq_len(err_idx)],
        .nm_bench_sim_error_line(parts, y_line = y_line),
        if (tail_start <= length(base_lines)) base_lines[tail_start:length(base_lines)]
      )
    }
    if (isTRUE(parts$use_ode) || parts$advan %in% c(6L, 10L, 13L)) {
      des_idx <- which(base_lines == "$DES")
      err_idx <- which(base_lines == "$ERROR")
      if (length(des_idx) == 1L && length(err_idx) == 1L && err_idx > des_idx) {
        des_body <- base_lines[(des_idx + 1L):(err_idx - 1L)]
        des_body <- des_body[!grepl("^\\s*F\\s*=", des_body)]
        base_lines <- c(
          base_lines[seq_len(des_idx)],
          des_body,
          base_lines[err_idx:length(base_lines)]
        )
      }
    }
  }
  tail <- if (mode == "sim") {
    c(
      sprintf("$SIM (%d NEW) ONLYSIM", as.integer(sim_seed)),
      .nm_bench_table_line(parts, sim_table, mode = "sim")
    )
  } else if (mode == "sim_struct") {
    c(
      sprintf("$SIM (%d NEW) ONLYSIM", as.integer(sim_seed)),
      .nm_bench_table_line(parts, sim_table, mode = "sim_struct")
    )
  } else {
    meth <- toupper(est_method)
    inter <- if (meth %in% c("FOCEI", "FOCE", "FOI")) " INTER" else ""
    c(
      sprintf("$EST METHOD=1%s MAXEVAL=9999 PRINT=1 NOABORT", inter),
      "$COVAR PRINT=E"
    )
  }
  paste(c(base_lines, tail), collapse = "\n")
}

#' Simulate observations with LibeRation at template true parameters
#'
#' Uses the same THETA/OMEGA/SIGMA as the NONMEM control template. Intended as
#' a stand-in until NONMEM \code{$SIM ONLYSIM} is wired for every ADVAN/TRANS pair.
#'
#' @param parts Template parts from \code{\link{nm_ctl_template}}.
#' @param design Event table from \code{\link{nm_bench_mixed_design}}.
#' @param seed Random seed.
#' @return \code{nm_dataset} with filled \code{DV}.
#' @examples
#' \dontrun{
#' design <- nm_bench_mixed_design(2L, 1L, n_per = 1L)
#' parts <- nm_ctl_template(2L, 1L)
#' nm_bench_simulate_rcppnm(parts, design, seed = 1L)
#' }
#' @export
nm_bench_simulate_rcppnm <- function(parts, design, seed = 12345L) {
  ctl_txt <- nm_ctl_compose(parts)
  pk_nm <- .nm_bench_pk_rcpp_lines(parts)
  ctl_lines <- .nm_bench_replace_pk_block(
    strsplit(ctl_txt, "\n", fixed = TRUE)[[1L]],
    pk_nm
  )
  tmp <- tempfile(fileext = ".ctl")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(ctl_lines, tmp, useBytes = TRUE)
  imp <- nm_read_nonmem(tmp, data_path = "bench.csv")
  imp$model$INPUT <- .nm_bench_rcppnm_input_cols(
    as.integer(parts$advan), as.integer(parts$trans)
  )
  model <- imp$model
  theta <- model$THETAS$Value
  omega <- model$OMEGAS$Value
  sigma <- model$SIGMAS$Value
  n_eta <- nrow(model$OMEGAS)
  ids <- unique(design$ID)
  set.seed(seed)
  eta_mat <- matrix(
    stats::rnorm(length(ids) * n_eta, sd = sqrt(pmax(omega, 0))),
    nrow = length(ids),
    ncol = n_eta
  )
  sim_args <- list(
    model = model,
    data = structure(list(data = design), class = "nm_dataset"),
    theta = theta,
    omega = omega,
    sigma = sigma,
    eta = eta_mat,
    seed = seed,
    pk_engine = "cpp",
    with_amounts = TRUE,
    n_state = .nm_bench_table_ncomp(parts, model)
  )
  out <- do.call(.nm_task_sim, sim_args)
  structure(list(data = out), class = "nm_dataset")
}

#' @keywords internal
.nm_bench_write_csv <- function(dat, path, input_cols = NULL) {
  dat <- as.data.frame(dat)
  if (!is.null(input_cols) && length(input_cols) > 0L) {
    for (col in input_cols) {
      if (!col %in% names(dat)) {
        dat[[col]] <- .nm_default_input_col(col)
      }
    }
    keep <- c(input_cols, setdiff(names(dat), input_cols))
    dat <- dat[, keep, drop = FALSE]
    dat <- dat[, input_cols, drop = FALSE]
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    utils::write.csv(dat, path, row.names = FALSE, na = ".")
    return(invisible(path))
  }
  dt <- data.table::as.data.table(dat)
  num_cols <- c("TIME", "AMT", "DV", "RATE", "II")
  for (cn in intersect(num_cols, names(dt))) {
    if (is.numeric(dt[[cn]])) {
      data.table::set(dt, j = cn, value = ifelse(
        is.na(dt[[cn]]), ".", format(dt[[cn]], trim = TRUE, scientific = FALSE)
      ))
    }
  }
  data.table::fwrite(dt, path, na = ".", quote = FALSE)
  invisible(path)
}

#' Run NONMEM via \code{nmfe73 ctl mod}
#'
#' Uses \code{-prdefault} so PREDPP object files ship with the NONMEM install
#' instead of recompiling against the system gfortran (which often mismatches
#' prebuilt \code{.mod} files on Windows).
#'
#' @param ctl_path Path to control stream.
#' @param mod_path Path for NONMEM output stub (extension added by NM).
#' @param work_dir Working directory.
#' @param extra_args Additional flags passed to \code{nmfe73} (after
#'   \code{-prdefault}).
#' @return List with paths and exit status.
#' @examples
#' \dontrun{
#' # Requires NONMEM installation
#' nm_bench_run_nonmem("path/to/model.ctl")
#' }
#' @export
nm_bench_run_nonmem <- function(ctl_path,
                                mod_path = "run",
                                work_dir = ".",
                                extra_args = character()) {
  exe <- .nm_nonmem_exe()
  if (!nzchar(exe)) {
    .nm_stop("NONMEM executable not found on PATH (expected nmfe73).")
  }
  owd <- getwd()
  on.exit(setwd(owd), add = TRUE)
  if (!dir.exists(work_dir)) {
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  }
  setwd(work_dir)
  ctl_bn <- basename(ctl_path)
  if (!file.exists(ctl_bn) && file.exists(ctl_path)) {
    file.copy(ctl_path, ctl_bn, overwrite = TRUE)
  }
  nm_flags <- c("-prdefault", extra_args)
  lic <- .nm_nonmem_license_file(exe)
  if (nzchar(lic)) {
    nm_flags <- c(nm_flags, paste0("-licfile=", shQuote(lic, type = "cmd")))
  }
  args <- c(shQuote(ctl_bn), shQuote(mod_path), nm_flags)
  cmd <- paste(shQuote(exe), paste(args, collapse = " "))
  path_prefix <- .nm_bench_nonmem_path_prefix(exe)
  old_path <- Sys.getenv("PATH", unset = "")
  if (length(path_prefix) > 0L) {
    Sys.setenv(PATH = paste(c(path_prefix, old_path), collapse = .Platform$path.sep))
  }
  on.exit({
    if (length(path_prefix) > 0L) {
      Sys.setenv(PATH = old_path)
    }
  }, add = TRUE)
  out <- system(cmd, intern = TRUE)
  status <- attr(out, "status")
  if (is.null(status)) {
    status <- 0L
  }
  mod_file <- file.path(work_dir, mod_path)
  ext_path <- file.path(work_dir, paste0(mod_path, ".ext"))
  license_ok <- !any(grepl("NONMEM LICENSE HAS EXPIRED", out, fixed = TRUE))
  license_ok <- license_ok && !any(grepl("ERROR reading license file", out, fixed = TRUE))
  if (file.exists(mod_file) && !license_ok) {
    mod_tail <- utils::tail(readLines(mod_file, warn = FALSE), 20L)
    license_ok <- !any(grepl("NONMEM LICENSE HAS EXPIRED", mod_tail, fixed = TRUE))
  }
  list(
    ctl = file.path(work_dir, ctl_bn),
    mod = mod_file,
    lst = file.path(work_dir, paste0(mod_path, ".lst")),
    ext = ext_path,
    status = status,
    license_ok = license_ok,
    log = out
  )
}

#' @keywords internal
.nm_ext_diag_cols <- function(cols, prefix) {
  if (length(cols) == 0L) {
    return(cols)
  }
  pat <- paste0("^", prefix, "\\((\\d+),\\1\\)")
  keep <- grepl(pat, cols)
  if (any(keep)) {
    cols[keep]
  } else {
    cols
  }
}

#' Parse THETA/OMEGA/SIGMA and OBJ from a NONMEM .ext file
#'
#' @param ext_path Path to \code{.ext} file.
#' @return List with \code{theta}, \code{omega}, \code{sigma}, \code{obj}.
#' @examples
#' # Path to a NONMEM .ext file
#' @export
nm_bench_read_ext <- function(ext_path) {
  if (!file.exists(ext_path)) {
    return(NULL)
  }
  lines <- readLines(ext_path, warn = FALSE)
  hdr <- grep("^\\s*ITERATION", lines, ignore.case = TRUE)
  if (length(hdr) == 0L) {
    hdr <- grep("^TABLE", lines, ignore.case = TRUE)
    if (length(hdr) == 0L) {
      return(NULL)
    }
    hdr <- hdr[[1L]] + 1L
  } else {
    hdr <- hdr[[1L]]
  }
  tbl <- utils::read.table(
    text = lines[hdr:length(lines)],
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (nrow(tbl) == 0L) {
    return(NULL)
  }
  valid <- is.finite(suppressWarnings(as.numeric(tbl$OBJ)))
  if ("ITERATION" %in% names(tbl)) {
    it <- suppressWarnings(as.numeric(tbl$ITERATION))
    valid <- valid & is.finite(it) & it >= 0
  }
  valid <- valid & suppressWarnings(as.numeric(tbl$OBJ)) > 0
  if (any(valid, na.rm = TRUE)) {
    tbl <- tbl[valid, , drop = FALSE]
  }
  row <- tbl[nrow(tbl), , drop = FALSE]
  theta_cols <- grep("^THETA", names(row), value = TRUE)
  omega_cols <- grep("^OMEGA\\(", names(row), value = TRUE)
  if (length(omega_cols) == 0L) {
    omega_cols <- grep("^OMEGA", names(row), value = TRUE)
  }
  sigma_cols <- grep("^SIGMA\\(", names(row), value = TRUE)
  if (length(sigma_cols) == 0L) {
    sigma_cols <- grep("^SIGMA", names(row), value = TRUE)
  }
  obj <- NA_real_
  if ("OBJ" %in% names(row)) {
    obj <- as.numeric(row$OBJ[[1L]])
  }
  om_vals <- as.numeric(row[, .nm_ext_diag_cols(omega_cols, "OMEGA"), drop = FALSE])
  sg_vals <- as.numeric(row[, .nm_ext_diag_cols(sigma_cols, "SIGMA"), drop = FALSE])
  sg_out <- .nm_sigma_var_to_sd(sg_vals)
  list(
    theta = as.numeric(row[, theta_cols, drop = TRUE]),
    omega = om_vals,
    sigma = sg_out,
    obj = obj,
    table = row
  )
}

#' Compare LibeRation fit to NONMEM .ext results
#'
#' @param fit An \code{nm_fit} object.
#' @param nm_ext Parsed output from \code{\link{nm_bench_read_ext}}.
#' @param rtol Relative tolerance for parameters.
#' @param obj_rtol Relative tolerance for objective.
#' @return List with pass/fail flags and diffs.
#' @examples
#' \dontrun{
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "FO",
#'               control = list(maxit = 5L, compute_inference = FALSE))
#' nm_bench_compare(fit, fit)
#' }
#' @export
nm_bench_compare <- function(fit,
                             nm_ext,
                             rtol = 0.15,
                             obj_rtol = 0.15) {
  if (is.null(nm_ext)) {
    return(list(ok = FALSE, reason = "missing NONMEM .ext"))
  }
  th <- fit$theta
  om <- fit$omega
  sg <- fit$sigma
  d_theta <- abs(th - nm_ext$theta) / pmax(abs(nm_ext$theta), .Machine$double.eps)
  d_omega <- abs(om - nm_ext$omega) / pmax(abs(nm_ext$omega), .Machine$double.eps)
  d_sigma <- abs(sg - nm_ext$sigma) / pmax(abs(nm_ext$sigma), .Machine$double.eps)
  theta_ok <- all(d_theta <= rtol, na.rm = TRUE)
  omega_ok <- length(om) == 0L || all(d_omega <= rtol, na.rm = TRUE)
  sigma_ok <- length(sg) == 0L || all(d_sigma <= rtol, na.rm = TRUE)
  obj_ok <- TRUE
  if (is.finite(fit$objective) && is.finite(nm_ext$obj) && nm_ext$obj > 0) {
    obj_ok <- abs(fit$objective - nm_ext$obj) / max(abs(nm_ext$obj), 1e-8) <= obj_rtol
  }
  list(
    ok = theta_ok && omega_ok && sigma_ok && obj_ok,
    d_theta = d_theta,
    d_omega = d_omega,
    d_sigma = d_sigma,
    obj_nm = nm_ext$obj,
    obj_rcpp = fit$objective,
    theta_nm = nm_ext$theta,
    theta_rcpp = th
  )
}

#' Pilot benchmark for one ADVAN/TRANS pair
#'
#' @inheritParams nm_bench_case
#' @examples
#' \dontrun{
#' nm_bench_pilot(advan = 2L, trans = 1L, n_per = 1L)
#' }
#' @export
nm_bench_pilot <- function(advan = 2L,
                           trans = 2L,
                           work_dir = NULL,
                           n_per = 3L,
                           seed = 12345L,
                           method = "FOCEI",
                           run_nonmem = TRUE) {
  nm_bench_case(
    advan = advan,
    trans = trans,
    work_dir = work_dir,
    n_per = n_per,
    seed = seed,
    method = method,
    run_nonmem = run_nonmem
  )
}

#' Full benchmark for one ADVAN/TRANS pair (sim + est, NONMEM vs LibeRation)
#'
#' @param advan ADVAN number.
#' @param trans TRANS number.
#' @param work_dir Output directory (created if missing).
#' @param n_per Subjects per regimen block in the mixed design.
#' @param seed Random seed for simulation.
#' @param method LibeRation estimation method.
#' @param run_nonmem Run NONMEM simulation and estimation steps.
#' @param run_est Run estimation comparison after structural simulation (\code{FALSE}
#'   for simulation-only benchmarks).
#' @param est_rtol Relative tolerance for parameter comparison.
#' @param sim_ipred_rtol Relative tolerance for structural IPRED comparison.
#' @return List with design, simulated data, fits, and comparisons.
#' @examples
#' \dontrun{
#' nm_bench_case(advan = 2L, trans = 1L, n_per = 1L)
#' }
#' @export
nm_bench_case <- function(advan = 2L,
                          trans = 2L,
                          work_dir = NULL,
                          n_per = 3L,
                          seed = 12345L,
                          method = "FOCEI",
                          run_nonmem = TRUE,
                          run_est = TRUE,
                          est_rtol = 0.15,
                          sim_ipred_rtol = 0.02) {
  advan <- as.integer(advan)
  trans <- as.integer(trans)
  if (!nm_ctl_is_valid_pair(advan, trans)) {
    .nm_stop(sprintf("Invalid ADVAN/TRANS pair: %d / %d", advan, trans))
  }
  if (is.null(work_dir)) {
    work_dir <- file.path(
      tempdir(),
      sprintf("nm_bench_adv%d_trans%d_%s", advan, trans, seed)
    )
  }
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  tag <- sprintf("adv%d_trans%d", advan, trans)
  status <- list(
    advan = advan, trans = trans, tag = tag, work_dir = work_dir,
    cpp_pk = nm_cpp_advan_supported(advan, trans),
    ok = FALSE, stage = "init", message = NA_character_
  )
  parts <- nm_ctl_template(
    advan, trans,
    data_file = "design.csv",
    problem = sprintf("LibeRation bench %s", tag)
  )
  parts$input_cols <- .nm_bench_nonmem_input_cols(advan, trans)
  design <- nm_bench_mixed_design(advan, trans, n_per = n_per, seed = seed)
  design_path <- file.path(work_dir, "design.csv")
  .nm_bench_write_csv(design, design_path, input_cols = parts$input_cols)
  parts$data_file <- "design.csv"

  sim_struct_ctl <- nm_bench_ctl(
    parts, mode = "sim_struct", sim_seed = seed, sim_table = "simtab_struct"
  )
  sim_struct_path <- file.path(work_dir, "sim_struct.ctl")
  writeLines(sim_struct_ctl, sim_struct_path, useBytes = TRUE)
  sim_full_ctl <- nm_bench_ctl(
    parts, mode = "sim", sim_seed = seed, sim_table = "simtab_full"
  )
  sim_full_path <- file.path(work_dir, "sim_full.ctl")
  writeLines(sim_full_ctl, sim_full_path, useBytes = TRUE)

  parts_struct <- .nm_bench_parts_for_mode(parts, "sim_struct")
  rcpp_struct <- nm_bench_simulate_rcppnm(parts_struct, design, seed = seed)

  nm_struct_tab <- NULL
  nm_full_tab <- NULL
  nm_run_struct <- NULL
  nm_run_full <- NULL
  if (isTRUE(run_nonmem) && nm_nonmem_available()) {
    nm_run_struct <- nm_bench_run_nonmem(
      sim_struct_path, mod_path = "sim_struct", work_dir = work_dir
    )
    nm_struct_tab <- nm_bench_read_simtab(file.path(work_dir, "simtab_struct"))
    if (isFALSE(nm_run_struct$license_ok)) {
      status$stage <- "nonmem_license"
      status$message <- "NONMEM license invalid or expired"
    }
    if (isTRUE(run_est)) {
      nm_run_full <- nm_bench_run_nonmem(
        sim_full_path, mod_path = "sim_full", work_dir = work_dir
      )
      nm_full_tab <- nm_bench_read_simtab(file.path(work_dir, "simtab_full"))
    }
  }

  sim_cmp <- if (!is.null(nm_struct_tab)) {
    nm_bench_compare_sim(
      rcpp_struct$data, nm_struct_tab,
      design = design, ipred_rtol = sim_ipred_rtol
    )
  } else {
    list(ok = NA, reason = "NONMEM structural sim not available")
  }

  if (!isTRUE(run_est)) {
    status$stage <- "done"
    status$ok <- isTRUE(sim_cmp$ok)
    status$sim_ok <- sim_cmp$ok
    status$est_ok <- NA
    status$message <- if (isTRUE(sim_cmp$ok)) {
      "pass"
    } else {
      sim_cmp$reason %||% "sim mismatch"
    }
    return(list(
      advan = advan,
      trans = trans,
      tag = tag,
      work_dir = work_dir,
      design = design,
      parts = parts,
      design_path = design_path,
      sim_struct_path = sim_struct_path,
      sim_full_path = sim_full_path,
      rcpp_struct = rcpp_struct,
      nm_struct_tab = nm_struct_tab,
      sim_compare = sim_cmp,
      status = status
    ))
  }

  dat <- if (!is.null(nm_full_tab)) {
    .nm_bench_data_from_simtab(design, nm_full_tab, parts$input_cols)
  } else {
    rcpp_full <- nm_bench_simulate_rcppnm(parts, design, seed = seed)
    rcpp_full$data
  }
  if (is.null(dat)) {
    status$stage <- "sim_data"
    status$message <- "Could not build estimation dataset"
    return(list(
      design = design, parts = parts, sim_compare = sim_cmp, status = status,
      advan = advan, trans = trans, tag = tag, work_dir = work_dir
    ))
  }
  dat$MDV <- as.integer(dat$MDV)
  dat$EVID <- as.integer(dat$EVID)
  data_path <- file.path(work_dir, "data.csv")
  .nm_bench_write_csv(dat, data_path, input_cols = parts$input_cols)
  parts$data_file <- "data.csv"
  est_ctl <- nm_bench_ctl(parts, mode = "est", est_method = method)
  est_ctl_path <- file.path(work_dir, "est.ctl")
  writeLines(est_ctl, est_ctl_path, useBytes = TRUE)

  rcpp_fit <- NULL
  est_cmp <- NULL
  imp <- tryCatch({
    rcpp_ctl <- nm_ctl_compose(parts)
    rcpp_lines <- .nm_bench_replace_pk_block(
      strsplit(rcpp_ctl, "\n", fixed = TRUE)[[1L]],
      .nm_bench_pk_rcpp_lines(parts)
    )
    tmp <- tempfile(fileext = ".ctl")
    on.exit(unlink(tmp), add = TRUE)
    writeLines(rcpp_lines, tmp, useBytes = TRUE)
    nm_read_nonmem(tmp, data_path = data_path)
  }, error = function(e) e)
  if (inherits(imp, "error")) {
    status$stage <- "rcpp_import"
    status$message <- imp$message
    return(list(
      design = design, data = dat, parts = parts,
      sim_compare = sim_cmp, status = status,
      advan = advan, trans = trans, tag = tag, work_dir = work_dir
    ))
  }
  imp$model$INPUT <- .nm_bench_rcppnm_input_cols(advan, trans)
  rcpp_fit <- tryCatch(
    nm_est(
      imp$model,
      structure(list(data = dat), class = "nm_dataset"),
      method = method,
      grad = "auto",
      pk_engine = "cpp",
      max_outer = if (identical(method, "FOCEI")) 3L else 8L,
      tol = if (identical(method, "FOCEI")) 1e-5 else 1e-4,
      control = list(maxit = 300, n_cores = 1L)
    ),
    error = function(e) e
  )
  if (inherits(rcpp_fit, "error")) {
    status$stage <- "rcpp_est"
    status$message <- rcpp_fit$message
    return(list(
      design = design, data = dat, parts = parts,
      sim_compare = sim_cmp, status = status,
      advan = advan, trans = trans, tag = tag, work_dir = work_dir
    ))
  }

  nm_run <- NULL
  nm_ext <- NULL
  if (isTRUE(run_nonmem) && nm_nonmem_available()) {
    nm_run <- nm_bench_run_nonmem(est_ctl_path, mod_path = "est", work_dir = work_dir)
    if (isFALSE(nm_run$license_ok)) {
      status$stage <- "nonmem_license"
      status$message <- "NONMEM license invalid or expired"
    }
    nm_ext <- nm_bench_read_ext(nm_run$ext)
  }
  est_cmp <- if (!is.null(nm_ext)) {
    nm_bench_compare(rcpp_fit, nm_ext, rtol = est_rtol)
  } else {
    NULL
  }

  status$stage <- "done"
  status$ok <- isTRUE(sim_cmp$ok) && (is.null(est_cmp) || isTRUE(est_cmp$ok))
  status$sim_ok <- sim_cmp$ok
  status$est_ok <- if (!is.null(est_cmp)) est_cmp$ok else NA
  status$message <- if (isTRUE(sim_cmp$ok) && (is.null(est_cmp) || isTRUE(est_cmp$ok))) {
    "pass"
  } else if (!isTRUE(sim_cmp$ok)) {
    sim_cmp$reason %||% "sim mismatch"
  } else {
    "est mismatch"
  }

  list(
    advan = advan,
    trans = trans,
    tag = tag,
    work_dir = work_dir,
    design = design,
    data = dat,
    parts = parts,
    design_path = design_path,
    sim_struct_path = sim_struct_path,
    sim_full_path = sim_full_path,
    est_ctl_path = est_ctl_path,
    rcpp_struct = rcpp_struct,
    nm_struct_tab = nm_struct_tab,
    nm_full_tab = nm_full_tab,
    sim_compare = sim_cmp,
    rcpp_fit = rcpp_fit,
    nm_run = nm_run,
    nm_ext = nm_ext,
    est_compare = est_cmp,
    compare = est_cmp,
    status = status
  )
}

#' Run benchmarks for all ADVAN/TRANS pairs
#'
#' @param work_dir Root directory for per-case subfolders.
#' @param pairs Data frame from \code{\link{nm_bench_pairs}} (\code{NULL} = all).
#' @param ... Passed to \code{\link{nm_bench_case}}.
#' @return List of per-case results plus summary data frame.
#' @examples
#' \dontrun{
#' # Long-running; use a temp work directory
#' nm_bench_run_all(work_dir = tempfile())
#' }
#' @export
nm_bench_run_all <- function(work_dir = NULL,
                             pairs = NULL,
                             ...) {
  if (is.null(pairs)) {
    pairs <- nm_bench_pairs()
  }
  if (is.null(work_dir)) {
    work_dir <- file.path(tempdir(), "nm_bench_all")
  }
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  results <- vector("list", nrow(pairs))
  for (i in seq_len(nrow(pairs))) {
    adv <- pairs$advan[[i]]
    tr <- pairs$trans[[i]]
    tag <- pairs$tag[[i]]
    message(sprintf("Benchmark %s (%d/%d)", tag, i, nrow(pairs)))
    case_dir <- file.path(work_dir, tag)
    results[[i]] <- tryCatch(
      nm_bench_case(advan = adv, trans = tr, work_dir = case_dir, ...),
      error = function(e) {
        list(
          advan = adv, trans = tr, tag = tag, work_dir = case_dir,
          status = list(
            ok = FALSE, stage = "error", message = conditionMessage(e)
          )
        )
      }
    )
  }
  names(results) <- pairs$tag
  summaries <- nm_bench_summarize(results)
  list(
    results = results,
    summary = summaries$combined,
    sim_summary = summaries$sim,
    est_summary = summaries$est,
    work_dir = work_dir
  )
}

#' Summarize simulation benchmark results
#'
#' @param results List of \code{nm_bench_case} outputs (named by tag).
#' @return Data frame with simulation comparison columns.
#' @examples
#' # Summary helper for nm_bench_run_all() results
#' @export
nm_bench_summarize_sim <- function(results) {
  nm_bench_summarize(results)$sim
}

#' Summarize estimation benchmark results
#'
#' @param results List of \code{nm_bench_case} outputs (named by tag).
#' @return Data frame with estimation comparison columns.
#' @examples
#' # Summary helper for nm_bench_run_all() results
#' @export
nm_bench_summarize_est <- function(results) {
  nm_bench_summarize(results)$est
}

#' Summarize benchmark results across ADVAN/TRANS pairs
#'
#' @param results List of \code{nm_bench_case} outputs (named by tag).
#' @return List with \code{combined}, \code{sim}, and \code{est} data frames.
#' @examples
#' # Summary helper for nm_bench_run_all() results
#' @export
nm_bench_summarize <- function(results) {
  tags <- names(results)
  if (is.null(tags)) {
    tags <- vapply(results, function(x) x$tag %||% "?", character(1L))
  }
  rows <- lapply(seq_along(results), function(i) {
    r <- results[[i]]
    st <- r$status %||% list()
    sim <- r$sim_compare %||% list()
    est <- r$est_compare %||% r$compare %||% list()
    obj_rcpp <- est$obj_rcpp %||% NA_real_
    if (is.na(obj_rcpp) && !is.null(r$rcpp_fit) && is.finite(r$rcpp_fit$objective %||% NA)) {
      obj_rcpp <- r$rcpp_fit$objective
    }
    reg_cols <- .nm_bench_flatten_sim_regimens(sim$by_regimen %||% list())
    c(
      list(
        tag = tags[[i]],
        advan = r$advan %||% NA_integer_,
        trans = r$trans %||% NA_integer_,
        stage = st$stage %||% NA_character_,
        ok = isTRUE(st$ok),
        sim_ok = sim$ok %||% NA,
        est_ok = est$ok %||% NA,
        max_ipred_rel = sim$max_ipred_rel %||% NA_real_,
        max_amt_rel = sim$max_amt_rel %||% NA_real_,
        max_dv_rel = sim$max_dv_rel %||% NA_real_,
        obj_nm = est$obj_nm %||% NA_real_,
        obj_rcpp = obj_rcpp,
        sim_message = gsub("[\r\n]+", " ", sim$reason %||% "", perl = TRUE),
        est_message = if (!is.null(est$ok) && !is.na(est$ok) && !isTRUE(est$ok)) {
          "est mismatch"
        } else {
          ""
        },
        message = gsub("[\r\n]+", " ", st$message %||% sim$reason %||% "", perl = TRUE),
        stringsAsFactors = FALSE
      ),
      reg_cols
    )
  })
  combined <- do.call(rbind, lapply(rows, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
  rownames(combined) <- NULL
  sim_cols <- c(
    "tag", "advan", "trans", "stage", "sim_ok",
    "max_ipred_rel", "max_amt_rel", "max_dv_rel",
    "sim_ok_single", "max_ipred_rel_single", "max_amt_rel_single", "n_obs_single",
    "sim_ok_multiple", "max_ipred_rel_multiple", "max_amt_rel_multiple", "n_obs_multiple",
    "sim_ok_steady_state", "max_ipred_rel_steady_state", "max_amt_rel_steady_state",
    "n_obs_steady_state",
    "sim_message", "message"
  )
  sim <- combined[, intersect(sim_cols, names(combined)), drop = FALSE]
  est <- combined[, c(
    "tag", "advan", "trans", "stage", "est_ok",
    "obj_nm", "obj_rcpp", "est_message", "message"
  ), drop = FALSE]
  list(combined = combined, sim = sim, est = est)
}
