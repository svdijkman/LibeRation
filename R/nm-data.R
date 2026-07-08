#' Load a NONMEM-style dataset
#'
#' @param path Path to a delimited file.
#' @param nm_table Logical; if \code{TRUE}, skip the first row (TABLE header).
#' @param ... Passed to \code{\link[data.table:fread]{data.table::fread}}.
#' @return An \code{nm_dataset} object.
#' @examples
#' sim <- nm_synthetic_theo(n_sub = 2L)
#' csv <- tempfile(fileext = ".csv")
#' write.csv(sim$data$data, csv, row.names = FALSE)
#' nm_dataset(csv)
#' @export
nm_dataset <- function(path, nm_table = FALSE, ...) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    .nm_stop("Package 'data.table' is required for nm_dataset().")
  }
  skip <- if (isTRUE(nm_table)) 1L else 0L
  dat <- data.table::fread(path, skip = skip, ...)
  nm_dataset_from_table(dat, path = path)
}

#' @rdname nm_dataset
#' @param dat A \code{data.table} or \code{data.frame}.
#' @param path Optional source path metadata.
#' @export
nm_dataset_from_table <- function(dat, path = NULL) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    .nm_stop("Package 'data.table' is required.")
  }
  dat <- data.table::as.data.table(dat)
  structure(
    list(data = dat, path = path),
    class = "nm_dataset"
  )
}

#' Default value when an $INPUT column is absent from the dataset
#' @keywords internal
.nm_default_input_col <- function(col) {
  col <- toupper(col)
  if (grepl("^F[0-9]+$", col) || grepl("^S[0-9]+$", col)) {
    return(1)
  }
  switch(
    col,
    F1 = 1,
    S1 = 1,
    S2 = 1,
    S3 = 1,
    S4 = 1,
    RATE = 0,
    AMT = 0,
    DV = 0,
    MDV = 0L,
    SS = 0L,
    II = 0,
    ADDL = 0L,
    TAU = 12,
    0
  )
}

#' Normalize NONMEM RATE=-1 (bolus) and RATE=-2 (duration in AMT) dose rows
#' @keywords internal
.nm_normalize_dose_records <- function(dat) {
  if (!"RATE" %in% names(dat)) {
    return(dat)
  }
  dose_idx <- which(dat$EVID %in% c(1L, 4L))
  if (length(dose_idx) == 0L) {
    return(dat)
  }
  for (i in dose_idx) {
    r <- as.numeric(dat$RATE[[i]])
    if (!is.finite(r)) {
      next
    }
    if (r == -1) {
      dat$RATE[[i]] <- 0
    } else if (r <= -2) {
      dur <- as.numeric(dat$AMT[[i]])
      dose_amt <- dur
      if ("DOSE" %in% names(dat)) {
        dose_amt <- as.numeric(dat$DOSE[[i]])
      }
      if (is.finite(dur) && dur > 0 && is.finite(dose_amt) && dose_amt > 0) {
        dat$RATE[[i]] <- dose_amt / dur
        dat$AMT[[i]] <- dose_amt
      }
    }
  }
  dat
}

#' Effective observation compartment for a row (CMT when EVID=0, else model default)
#' @keywords internal
.nm_row_obs_cmp <- function(model, subj, i = seq_len(nrow(subj))) {
  default <- as.integer(if (!is.null(model$OBSCMP)) model$OBSCMP else 1L)
  if (!"CMT" %in% names(subj) || !"EVID" %in% names(subj)) {
    return(default)
  }
  out <- rep(default, length(i))
  obs <- subj$EVID[i] == 0L & subj$MDV[i] == 0L
  cmt <- as.integer(subj$CMT[i])
  use <- obs & is.finite(cmt) & cmt > 0L
  out[use] <- cmt[use]
  out
}

#' Effective dose compartment for a row (CMT when dosing, else model default)
#' @keywords internal
.nm_row_dose_cmp <- function(model, subj, i = seq_len(nrow(subj))) {
  default <- as.integer(if (!is.null(model$DOSECMP)) model$DOSECMP else 1L)
  if (!"CMT" %in% names(subj) || !"EVID" %in% names(subj)) {
    return(default)
  }
  out <- rep(default, length(i))
  dose <- subj$EVID[i] %in% c(1L, 4L)
  cmt <- as.integer(subj$CMT[i])
  use <- dose & is.finite(cmt) & cmt > 0L
  out[use] <- cmt[use]
  out
}

#' @keywords internal
.nm_max_scale_n <- function() {
  10L
}

#' Maximum number of compartment scaling factors (S1..Sn) supported
#' @return Integer scalar (10, matching NONMEM $MODEL NCOMPS).
#' @examples
#' nm_max_scale_n()
#' @export
nm_max_scale_n <- function() {
  .nm_max_scale_n()
}

#' Sx column names S1..Sn
#' @keywords internal
.nm_scale_col_names <- function(n = .nm_max_scale_n()) {
  paste0("S", seq_len(as.integer(n)))
}

#' Sx symbols assigned in $PK (exclude from data defaults when set in model)
#' @keywords internal
.nm_pred_scale_cols <- function(model) {
  if (is.null(model) || is.null(model$PRED) || !nzchar(model$PRED)) {
    return(character())
  }
  lines <- .nm_split_lines(model$PRED)
  lhs <- toupper(trimws(vapply(
    strsplit(lines, "=", fixed = TRUE),
    function(x) x[[1L]],
    character(1L)
  )))
  intersect(.nm_scale_col_names(), lhs)
}

.nm_f_col_names <- function(n = .nm_max_scale_n()) {
  paste0("F", seq_len(as.integer(n)))
}

#' Fx symbols assigned in $PK (exclude from data defaults when set in model)
#' @keywords internal
.nm_pred_f_cols <- function(model) {
  if (is.null(model) || is.null(model$PRED) || !nzchar(model$PRED)) {
    return(character())
  }
  lines <- .nm_split_lines(model$PRED)
  lhs <- toupper(trimws(vapply(
    strsplit(lines, "=", fixed = TRUE),
    function(x) x[[1L]],
    character(1L)
  )))
  intersect(.nm_f_col_names(), lhs)
}

#' Row-level S/F vectors for C++ PK (data override when user supplied columns)
#' @keywords internal
.nm_input_event_vectors <- function(subj) {
  empty <- numeric(0)
  user_scale_cols <- attr(subj, "user_scale_cols")
  if (is.null(user_scale_cols)) {
    user_scale_cols <- character()
  }
  user_f_cols <- attr(subj, "user_f_cols")
  if (is.null(user_f_cols)) {
    user_f_cols <- character()
  }
  scale_names <- .nm_scale_col_names()
  f_names <- .nm_f_col_names()
  scale_mat <- matrix(NA_real_, nrow = nrow(subj), ncol = length(scale_names))
  colnames(scale_mat) <- scale_names
  f_mat <- matrix(NA_real_, nrow = nrow(subj), ncol = length(f_names))
  colnames(f_mat) <- f_names
  has_scale <- FALSE
  has_f <- FALSE
  for (j in seq_along(scale_names)) {
    sc <- scale_names[[j]]
    if (sc %in% names(subj)) {
      scale_mat[, j] <- as.numeric(subj[[sc]])
      has_scale <- TRUE
    }
  }
  for (j in seq_along(f_names)) {
    fc <- f_names[[j]]
    if (fc %in% names(subj)) {
      f_mat[, j] <- as.numeric(subj[[fc]])
      has_f <- TRUE
    }
  }
  list(
    s1 = if ("S1" %in% names(subj)) as.numeric(subj$S1) else empty,
    s2 = if ("S2" %in% names(subj)) as.numeric(subj$S2) else empty,
    s3 = if ("S3" %in% names(subj)) as.numeric(subj$S3) else empty,
    s4 = if ("S4" %in% names(subj)) as.numeric(subj$S4) else empty,
    scale_mat = if (has_scale) scale_mat else matrix(numeric(0), 0L, 0L),
    use_data_scale = length(user_scale_cols) > 0L,
    f_mat = if (has_f) f_mat else matrix(numeric(0), 0L, 0L),
    use_data_f = length(user_f_cols) > 0L
  )
}

#' Row-level S vectors for C++ PK (data override when user supplied S columns)
#' @keywords internal
.nm_scale_event_vectors <- function(subj) {
  .nm_input_event_vectors(subj)
}

#' @keywords internal
.nm_prepare_data <- function(data, input_cols, model = NULL) {
  if (inherits(data, "nm_dataset")) {
    dat <- data.table::as.data.table(data$data)
  } else if (is.data.frame(data)) {
    dat <- data.table::as.data.table(data)
  } else {
    .nm_stop("data must be an nm_dataset or data.frame/data.table.")
  }
  .nm_require_cols(dat, c("ID", "TIME", "EVID"), "dataset")
  pk_scale <- .nm_pred_scale_cols(model)
  pk_f <- .nm_pred_f_cols(model)
  user_scale_cols <- intersect(.nm_scale_col_names(), names(dat))
  user_f_cols <- intersect(.nm_f_col_names(), names(dat))
  for (col in input_cols) {
    col_u <- toupper(col)
    if (col_u %in% pk_scale || col_u %in% pk_f) {
      next
    }
    if (!col %in% names(dat)) {
      dat[[col]] <- .nm_default_input_col(col)
    }
  }
  if (!"MDV" %in% names(dat)) {
    dat$MDV <- as.integer(dat$EVID != 0L)
  }
  if (!"DV" %in% names(dat)) {
    dat$DV <- 0
  }
  if (!"AMT" %in% names(dat)) {
    dat$AMT <- 0
  }
  if (!"CMT" %in% names(dat)) {
    dat$CMT <- 1L
  }
  if (!"RATE" %in% names(dat)) {
    dat$RATE <- 0
  }
  if (!"SS" %in% names(dat)) {
    dat$SS <- 0L
  }
  if (!"II" %in% names(dat)) {
    dat$II <- 0
  }
  if (!"ADDL" %in% names(dat)) {
    dat$ADDL <- 0L
  }
  if (!"cEVID" %in% names(dat)) {
    dat$cEVID <- ave(seq_len(nrow(dat)), dat$ID, FUN = seq_along)
  }
  if (!"TAU" %in% names(dat)) {
    dat$TAU <- 12
  }
  if (!"F1" %in% names(dat)) {
    dat$F1 <- 1
  } else if (all(as.numeric(dat$F1) == 0, na.rm = TRUE)) {
    dat$F1 <- 1
  }
  dat <- .nm_normalize_dose_records(dat)
  for (sc in .nm_scale_col_names()) {
    if (sc %in% pk_scale) {
      next
    }
    if (!sc %in% names(dat)) {
      dat[[sc]] <- 1
    }
  }
  for (fc in .nm_f_col_names()) {
    if (fc %in% pk_f) {
      next
    }
    if (!fc %in% names(dat)) {
      dat[[fc]] <- 1
    }
  }
  if (!"KA" %in% names(dat)) {
    dat$KA <- 1
  }
  attr(dat, "user_scale_cols") <- user_scale_cols
  attr(dat, "user_f_cols") <- user_f_cols
  dat[]
}

#' @keywords internal
.nm_subject_ids <- function(dat) {
  sort(unique(dat$ID))
}

#' @keywords internal
.nm_subject_slice <- function(dat, id) {
  dat[ID == id][order(TIME, -EVID)]
}

#' Effective inter-dose interval for a dosing row (II, else TAU).
#' @keywords internal
.nm_effective_ii <- function(row) {
  ii <- if ("II" %in% names(row)) row$II else 0
  if (length(ii) == 0L || is.na(ii) || ii <= 0) {
    tau <- if ("TAU" %in% names(row)) row$TAU else 0
    if (length(tau) == 0L || is.na(tau) || tau <= 0) {
      return(0)
    }
    return(as.numeric(tau))
  }
  as.numeric(ii)
}

#' Expand NONMEM ADDL/II implied doses into explicit dosing events.
#' ADDL is the number of additional doses; II (or TAU) is the interval.
#' @keywords internal
.nm_expand_addl <- function(subj) {
  if (!"ADDL" %in% names(subj)) {
    return(subj)
  }
  subj <- data.table::as.data.table(subj)
  dose_idx <- which(subj$EVID %in% c(1L, 4L) & subj$ADDL > 0L)
  if (length(dose_idx) == 0L) {
    return(subj[])
  }
  extra_rows <- list()
  for (i in dose_idx) {
    addl <- as.integer(subj$ADDL[i])
    ii <- .nm_effective_ii(subj[i])
    if (addl <= 0L || ii <= 0) {
      next
    }
    for (k in seq_len(addl)) {
      new_row <- data.table::copy(subj[i])
      new_row[, TIME := TIME + k * ii]
      new_row[, ADDL := 0L]
      if ("SS" %in% names(new_row)) {
        new_row[, SS := 0L]
      }
      extra_rows[[length(extra_rows) + 1L]] <- new_row
    }
    subj[i, ADDL := 0L]
  }
  if (length(extra_rows) == 0L) {
    return(subj[])
  }
  out <- data.table::rbindlist(c(list(subj), extra_rows), use.names = TRUE, fill = TRUE)
  out[order(TIME, -EVID)][]
}

#' Map ipred from an expanded event table back to the original subject rows.
#'
#' When several rows share the same (TIME, EVID, CMT) key (duplicate/replicate
#' timepoints), a plain \code{match()} maps every duplicate original row to the
#' \emph{first} matching event row. That is fixed here by matching positionally:
#' the k-th original occurrence of a key maps to the k-th expanded occurrence
#' (falling back to the last available occurrence).
#' @keywords internal
.nm_ipred_align <- function(subj, subj_ev, ipred_ev) {
  ipred_ev <- as.numeric(ipred_ev)
  if (nrow(subj) == nrow(subj_ev)) {
    return(ipred_ev)
  }
  ev_key <- paste(subj_ev$TIME, subj_ev$EVID,
                  if ("CMT" %in% names(subj_ev)) subj_ev$CMT else 1L)
  orig_key <- paste(subj$TIME, subj$EVID,
                    if ("CMT" %in% names(subj)) subj$CMT else 1L)
  if (!anyDuplicated(orig_key)) {
    return(ipred_ev[match(orig_key, ev_key)])
  }
  ev_by_key <- split(seq_along(ev_key), ev_key)
  orig_occ <- stats::ave(seq_along(orig_key), orig_key, FUN = seq_along)
  out <- rep(NA_real_, length(orig_key))
  for (i in seq_along(orig_key)) {
    idxs <- ev_by_key[[orig_key[i]]]
    if (is.null(idxs)) {
      next
    }
    k <- min(orig_occ[i], length(idxs))
    out[i] <- ipred_ev[idxs[k]]
  }
  out
}

#' Subject event table with ADDL/II doses expanded for PK simulation.
#' @keywords internal
.nm_subject_events <- function(subj) {
  .nm_expand_addl(subj)
}
