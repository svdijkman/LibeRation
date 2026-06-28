#' Import a NONMEM control stream
#'
#' Parses \code{$PROB}, \code{$INPUT}, \code{$DATA}, \code{$THETA},
#' \code{$OMEGA}, \code{$SIGMA}, \code{$PK}/\code{$PRED}, \code{$ERROR},
#' \code{$DES}, \code{$COVAR}, \code{$EST}, \code{$BAYES}, and FIX/bounds.
#'
#' @param path Path to a \code{.ctl} file.
#' @param data_path Optional override for \code{$DATA} filename.
#' @return List with \code{model}, \code{data_path}, \code{method}, and parsed blocks.
#' @examples
#' ctl <- system.file("extdata", "theo.ctl", package = "LibeRation")
#' if (nzchar(ctl)) nm_read_nonmem(ctl)
#' @export
nm_read_nonmem <- function(path, data_path = NULL) {
  if (!file.exists(path)) {
    .nm_stop("Control file not found: ", path)
  }
  txt <- readLines(path, warn = FALSE)
  blocks <- .nm_parse_ctl_blocks(txt)
  prob <- blocks$PROB
  advan <- 2L
  trans <- 2L
  use_ode <- FALSE
  if (length(prob) > 0L) {
    prob_line <- paste(prob, collapse = " ")
    m <- regmatches(prob_line, gregexpr("ADVAN\\s*=\\s*[0-9]+", prob_line, ignore.case = TRUE))[[1]]
    if (length(m)) advan <- as.integer(sub(".*=\\s*", "", m))
    m2 <- regmatches(prob_line, gregexpr("TRANS\\s*=\\s*[0-9]+", prob_line, ignore.case = TRUE))[[1]]
    if (length(m2)) trans <- as.integer(sub(".*=\\s*", "", m2))
    if (grepl("ADVAN\\s*=\\s*(6|10|13)", prob_line, ignore.case = TRUE)) {
      use_ode <- TRUE
    }
  }
  input_cols <- .nm_parse_input_block(blocks$INPUT)
  if (is.null(data_path) && length(blocks$DATA) > 0L) {
    data_path <- .nm_parse_data_line(blocks$DATA[1])
  }
  thetas <- .nm_parse_theta_block(blocks$THETA)
  omegas <- .nm_parse_omega_block(blocks$OMEGA)
  sigmas <- .nm_parse_sigma_block(blocks$SIGMA)
  covariates <- .nm_parse_covar_block(blocks$COVAR)
  dosecmp <- 1L
  obscmp <- nm_ctl_default_obscmp(advan, trans)
  pred <- .nm_ctl_merge_pk_pred(blocks$PK, blocks$PRED)
  error <- paste(blocks$ERROR, collapse = "\n")
  des <- paste(blocks$DES, collapse = "\n")
  if (!nzchar(error)) {
    error <- "Y = F"
  }
  lik <- .nm_lik_from_error_block(error)
  sigmas <- .nm_sigmas_nm_to_rcpp(sigmas, lik)
  model <- nm_model(
    INPUT = input_cols,
    ADVAN = advan,
    TRANS = trans,
    DOSECMP = dosecmp,
    OBSCMP = obscmp,
    USE_ODE = use_ode,
    PRED = pred,
    ERROR = error,
    DES = des,
    THETAS = thetas,
    OMEGAS = omegas,
    SIGMAS = sigmas,
    COVARIATES = covariates,
    LIK_CONFIG = lik
  )
  method <- .nm_parse_est_method(blocks$EST)
  bayes_opts <- .nm_parse_bayes_block(blocks$BAYES)
  list(
    model = model,
    data_path = data_path,
    method = method,
    blocks = blocks,
    bayes = bayes_opts
  )
}

#' @keywords internal
.nm_parse_ctl_blocks <- function(lines) {
  blocks <- list()
  cur <- NULL
  for (ln in lines) {
    if (grepl("^\\s*;", ln)) next
    if (grepl("^\\s*\\$", ln)) {
      hdr <- sub("^\\s*\\$", "", ln)
      parts <- strsplit(trimws(hdr), "\\s+")[[1]]
      cur <- parts[[1L]]
      if (identical(cur, "PROBLEM")) {
        cur <- "PROB"
      }
      blocks[[cur]] <- character()
      if (length(parts) > 1L) {
        blocks[[cur]] <- paste(parts[-1L], collapse = " ")
      }
      next
    }
    if (!is.null(cur)) {
      blocks[[cur]] <- c(blocks[[cur]], ln)
    }
  }
  if (!is.null(blocks$PK) && is.null(blocks$PRED)) {
    blocks$PRED <- blocks$PK
  }
  blocks
}

#' @keywords internal
.nm_ctl_merge_pk_pred <- function(pk, pred) {
  pk <- pk %||% character()
  pred <- pred %||% character()
  pk <- pk[nzchar(trimws(pk))]
  pred <- pred[nzchar(trimws(pred))]
  pk_txt <- trimws(paste(pk, collapse = "\n"))
  pred_txt <- trimws(paste(pred, collapse = "\n"))
  if (!nzchar(pk_txt)) {
    return(pred_txt)
  }
  if (!nzchar(pred_txt)) {
    return(pk_txt)
  }
  if (identical(pk_txt, pred_txt)) {
    return(pk_txt)
  }
  trimws(paste(c(pk_txt, pred_txt), collapse = "\n"))
}

#' @keywords internal
.nm_parse_input_block <- function(lines) {
  if (length(lines) == 0L) {
    return(c("ID", "TIME", "DV", "AMT", "MDV", "EVID", "CMT", "RATE"))
  }
  cols <- character()
  for (ln in lines) {
    ln <- trimws(ln)
    if (!nzchar(ln)) next
    toks <- strsplit(ln, "\\s+")[[1]]
    toks <- toks[nzchar(toks)]
    cols <- c(cols, toupper(toks))
  }
  unique(cols)
}

#' @keywords internal
.nm_parse_theta_block <- function(lines) {
  if (length(lines) == 0L) {
    return(data.frame(THETA = 1L, Value = 1, FIX = FALSE, Lower = NA_real_, Upper = NA_real_))
  }
  rows <- list()
  k <- 0L
  for (ln in lines) {
    if (grepl("FIX", ln, ignore.case = TRUE)) next
    nums <- as.numeric(regmatches(ln, gregexpr("[0-9]+\\.?[0-9]*([eE][+-]?[0-9]+)?", ln))[[1]])
    nums <- nums[!is.na(nums)]
    if (length(nums) == 0L) next
    k <- k + 1L
    fix <- grepl("\\bFIX\\b", ln, ignore.case = TRUE)
    if (grepl("\\(", ln)) {
      lower <- if (length(nums) >= 1L) nums[1] else NA_real_
      val <- if (length(nums) >= 2L) nums[2] else nums[1]
      upper <- if (length(nums) >= 3L) nums[3] else NA_real_
    } else {
      val <- nums[1]
      lower <- if (length(nums) >= 2L) nums[2] else NA_real_
      upper <- if (length(nums) >= 3L) nums[3] else NA_real_
    }
    rows[[k]] <- data.frame(
      THETA = k, Value = val, FIX = fix,
      Lower = lower, Upper = upper,
      row.names = NULL
    )
  }
  if (length(rows) == 0L) {
    return(data.frame(THETA = 1L, Value = 1, FIX = FALSE, Lower = NA_real_, Upper = NA_real_))
  }
  do.call(rbind, rows)
}

#' @keywords internal
.nm_parse_omega_block <- function(lines) {
  .nm_parse_matrix_block(lines, "OMEGA")
}

#' @keywords internal
.nm_parse_sigma_block <- function(lines) {
  .nm_parse_matrix_block(lines, "SIGMA")
}

#' @keywords internal
.nm_parse_matrix_block <- function(lines, prefix) {
  if (length(lines) == 0L) {
    return(data.frame(
      stats::setNames(data.frame(integer(), numeric()), c(prefix, "Value"))
    ))
  }
  vals <- vapply(lines, function(ln) {
    m <- regmatches(ln, regexpr("[0-9]+\\.?[0-9]*([eE][+-]?[0-9]+)?", ln))
    if (length(m) == 0L || !nzchar(m)) {
      return(NA_real_)
    }
    as.numeric(m)
  }, numeric(1))
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0L) {
    return(data.frame(
      stats::setNames(data.frame(integer(), numeric()), c(prefix, "Value"))
    ))
  }
  data.frame(
    stats::setNames(
      data.frame(seq_along(vals), vals),
      c(prefix, "Value")
    ),
    row.names = NULL
  )
}

#' @keywords internal
.nm_parse_covar_block <- function(lines) {
  if (length(lines) == 0L) {
    return(NULL)
  }
  cols <- character()
  for (ln in lines) {
    m <- regmatches(ln, gregexpr("[A-Za-z][A-Za-z0-9_]*", ln))[[1]]
    cols <- c(cols, toupper(m))
  }
  unique(cols[nzchar(cols)])
}

#' @keywords internal
.nm_sigma_var_to_sd <- function(values) {
  sqrt(pmax(as.numeric(values), 0))
}

#' @keywords internal
.nm_sigma_sd_to_var <- function(values) {
  pmax(as.numeric(values), 0)^2
}

#' NONMEM stores diagonal SIGMA as variances; Rcpp residual formulas use SDs.
#' @keywords internal
.nm_sigmas_nm_to_rcpp <- function(sigmas, lik = NULL) {
  if (is.null(sigmas) || nrow(sigmas) == 0L) {
    return(sigmas)
  }
  sigmas$Value <- .nm_sigma_var_to_sd(sigmas$Value)
  sigmas
}

#' @keywords internal
.nm_lik_from_error_block <- function(error) {
  err <- "propadd"
  ar1 <- FALSE
  err_uc <- toupper(error)
  if (grepl("F\\s*\\*\\s*\\(\\s*1\\s*[+-]\\s*ERR", err_uc)) {
    err <- if (grepl("\\+\\s*ERR\\(\\s*2\\s*\\)", err_uc)) "propadd" else "prop"
  } else if (grepl("IPRED\\s*/\\s*DV", error, ignore.case = TRUE)) {
    err <- "prop"
  } else if (grepl("DV\\s*-\\s*IPRED", error, ignore.case = TRUE)) {
    err <- "add"
  } else if (grepl("LOG\\(", error, ignore.case = TRUE)) {
    err <- "log"
  }
  if (grepl("CORR|AR1|EPS\\(", error, ignore.case = TRUE)) {
    ar1 <- TRUE
  }
  nm_lik_config(
    error = err,
    sigma_corr = if (ar1) "ar1" else "indep",
    ar1_rho = if (ar1) 0.5 else 0.0
  )
}

#' @keywords internal
.nm_parse_est_method <- function(est_lines) {
  method <- "FO"
  if (length(est_lines) == 0L) {
    return(method)
  }
  est_line <- toupper(paste(est_lines, collapse = " "))
  if (grepl("FOCEI|FOCE\\s*INTER", est_line)) method <- "FOCEI"
  else if (grepl("FOCE", est_line)) method <- "FOCE"
  else if (grepl("SAEM", est_line)) method <- "SAEM"
  else if (grepl("LAPLACE", est_line)) method <- "LAPLACE"
  else if (grepl("IMP", est_line)) method <- "IMP"
  else if (grepl("BAYES", est_line)) method <- "BAYES"
  method
}

#' @keywords internal
.nm_parse_bayes_block <- function(lines) {
  if (length(lines) == 0L) {
    return(NULL)
  }
  txt <- toupper(paste(lines, collapse = " "))
  sampler <- if (grepl("NUTS", txt)) "nuts" else if (grepl("HMC", txt)) "hmc" else "mh"
  list(sampler = tolower(sampler))
}

#' @keywords internal
.nm_parse_data_line <- function(line) {
  line <- trimws(line)
  if (!nzchar(line)) {
    return(NULL)
  }
  parts <- strsplit(line, "\\s+")[[1]]
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0L) {
    return(NULL)
  }
  parts[[1L]]
}

`%||%` <- function(x, y) if (is.null(x)) y else x
