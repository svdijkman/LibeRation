#' Parse a NONMEM control stream into editable GUI parts
#'
#' @param ctl_text Full control stream text.
#' @return List of model components (no \code{$EST} block).
#' @examples
#' parts <- nm_ctl_template(2L, 1L)
#' ctl <- nm_ctl_compose(parts)
#' nm_ctl_parse(ctl)
#' @export
nm_ctl_parse <- function(ctl_text) {
  lines <- unlist(strsplit(ctl_text, "\n", fixed = TRUE))
  blocks <- .nm_parse_ctl_blocks(lines)
  prob_txt <- if (length(blocks$PROB) > 0L) trimws(paste(blocks$PROB, collapse = " ")) else ""
  advan <- 2L
  trans <- 2L
  use_ode <- FALSE
  if (nzchar(prob_txt)) {
    m <- regmatches(prob_txt, gregexpr("ADVAN\\s*=\\s*[0-9]+", prob_txt, ignore.case = TRUE))[[1]]
    if (length(m)) advan <- as.integer(sub(".*=\\s*", "", m))
    m2 <- regmatches(prob_txt, gregexpr("TRANS\\s*=\\s*[0-9]+", prob_txt, ignore.case = TRUE))[[1]]
    if (length(m2)) trans <- as.integer(sub(".*=\\s*", "", m2))
    prob_txt <- gsub("ADVAN\\s*=\\s*[0-9]+", "", prob_txt, ignore.case = TRUE)
    prob_txt <- gsub("TRANS\\s*=\\s*[0-9]+", "", prob_txt, ignore.case = TRUE)
    prob_txt <- trimws(gsub("\\s+", " ", prob_txt))
  }
  if (advan %in% c(6L, 13L)) {
    use_ode <- TRUE
  }
  subroutine <- ""
  if (length(blocks$SUBROUTINE) > 0L) {
    subroutine <- toupper(trimws(paste(blocks$SUBROUTINE, collapse = " ")))
  }
  data_file <- NULL
  if (length(blocks$DATA) > 0L) {
    data_file <- .nm_parse_data_line(blocks$DATA[1])
  }
  input_cols <- .nm_parse_input_block(blocks$INPUT)
  output_cols <- .nm_parse_input_block(blocks$OUTPUT)
  pk <- .nm_ctl_merge_pk_pred(blocks$PK, blocks$PRED)
  error <- trimws(paste(blocks$ERROR, collapse = "\n"))
  des <- trimws(paste(blocks$DES, collapse = "\n"))
  list(
    problem = prob_txt,
    advan = advan,
    trans = trans,
    use_ode = use_ode,
    subroutine = subroutine,
    data_file = data_file,
    input_cols = input_cols,
    output_cols = output_cols,
    thetas = .nm_parse_theta_block(blocks$THETA),
    omegas = .nm_parse_omega_block(blocks$OMEGA),
    sigmas = {
      sg <- .nm_parse_sigma_block(blocks$SIGMA)
      if (!is.null(sg) && nrow(sg) > 0L) {
        .nm_sigmas_nm_to_rcpp(sg)
      } else {
        sg
      }
    },
    pk = pk,
    des = des,
    error = if (nzchar(error)) error else "Y = F"
  )
}

#' Compose a NONMEM control stream from GUI parts (no \code{$EST})
#'
#' @param parts List as returned by \code{\link{nm_ctl_parse}}.
#' @return Character scalar control stream text.
#' @examples
#' parts <- nm_ctl_template(2L, 1L)
#' nm_ctl_compose(parts)
#' @export
nm_ctl_compose <- function(parts) {
  .nm_coalesce <- function(x, y) if (is.null(x)) y else x
  prob <- trimws(.nm_coalesce(parts$problem, ""))
  advan <- as.integer(.nm_coalesce(parts$advan, 2L))
  trans <- as.integer(.nm_coalesce(parts$trans, nm_ctl_default_trans(advan)))
  trans <- nm_ctl_effective_trans(advan, trans)
  prob_line <- if (nzchar(prob)) {
    sprintf("$PROBLEM %s ADVAN=%d TRANS=%d", prob, advan, trans)
  } else {
    sprintf("$PROBLEM ADVAN=%d TRANS=%d", advan, trans)
  }
  input_cols <- parts$input_cols
  if (is.null(input_cols) || length(input_cols) == 0L) {
    input_cols <- c("ID", "TIME", "DV", "AMT", "MDV", "EVID", "CMT")
  }
  data_rel <- .nm_coalesce(parts$data_file, "data.csv")
  data_line <- sprintf("$DATA %s IGNORE=@", basename(data_rel))
  thetas <- .nm_ctl_normalize_thetas(parts$thetas)
  theta_lines <- vapply(seq_len(nrow(thetas)), function(i) {
    lo <- thetas$Lower[[i]]
    val <- thetas$Value[[i]]
    up <- thetas$Upper[[i]]
    if (length(lo) != 1L || is.na(lo)) lo <- 0
    if (length(val) != 1L || is.na(val)) val <- 1
    line <- if (length(up) == 1L && !is.na(up) && is.finite(up)) {
      sprintf(" (%g, %g, %g)", lo, val, up)
    } else {
      sprintf(" (%g, %g)", lo, val)
    }
    if (isTRUE(thetas$FIX[[i]])) {
      line <- paste0(line, " FIX")
    }
    lbl <- if ("Label" %in% names(thetas)) thetas$Label[[i]] else ""
    if (length(lbl) == 1L && !is.na(lbl) && nzchar(lbl)) {
      line <- paste0(line, " ; ", lbl)
    } else {
      line <- paste0(line, sprintf(" ; THETA%d", i))
    }
    line
  }, character(1L))
  omegas <- parts$omegas
  if (is.null(omegas) || nrow(omegas) == 0L) {
    omegas <- data.frame(OMEGA = 1L, Value = 0.1)
  }
  omega_lines <- vapply(seq_len(nrow(omegas)), function(i) {
    sprintf(" %g", as.numeric(omegas$Value[i]))
  }, character(1L))
  sigmas <- parts$sigmas
  if (is.null(sigmas) || nrow(sigmas) == 0L) {
    sigmas <- data.frame(SIGMA = 1L, Value = 0.1)
  }
  sigma_lines <- vapply(seq_len(nrow(sigmas)), function(i) {
    sprintf(" %g", .nm_sigma_sd_to_var(as.numeric(sigmas$Value[i])))
  }, character(1L))
  pk_lines <- strsplit(trimws(.nm_coalesce(parts$pk, "")), "\n", fixed = TRUE)[[1L]]
  pk_lines <- pk_lines[nzchar(trimws(pk_lines))]
  if (length(pk_lines) == 0L) {
    pk_lines <- " CL = THETA(1)"
  }
  if (isTRUE(parts$nonmem_pk)) {
    pk_lines <- c(pk_lines, .nm_ctl_pk_nonmem_extras(advan, trans))
  }
  des_txt <- trimws(.nm_coalesce(parts$des, ""))
  des_lines <- character()
  if (nzchar(des_txt)) {
    des_lines <- strsplit(des_txt, "\n", fixed = TRUE)[[1L]]
    des_lines <- des_lines[nzchar(trimws(des_lines))]
  }
  err_lines <- strsplit(trimws(.nm_coalesce(parts$error, "Y = F")), "\n", fixed = TRUE)[[1L]]
  err_lines <- err_lines[nzchar(trimws(err_lines))]
  out <- c(
    prob_line,
    paste("$INPUT", paste(input_cols, collapse = " ")),
    data_line
  )
  output_cols <- parts$output_cols
  if (!is.null(output_cols) && length(output_cols) > 0L) {
    out <- c(out, paste("$OUTPUT", paste(output_cols, collapse = " ")))
  }
  out <- c(
    out,
    "$THETA",
    theta_lines,
    "$OMEGA",
    omega_lines,
    "$SIGMA",
    sigma_lines,
    "$PK",
    pk_lines
  )
  if (length(des_lines) > 0L) {
    out <- c(out, "$DES", des_lines)
  }
  out <- c(out, "$ERROR", err_lines)
  sub_line <- .nm_coalesce(parts$subroutine, nm_ctl_subroutine_text(advan, trans))
  if (nzchar(sub_line)) {
    if (!grepl("^\\$", sub_line)) {
      sub_line <- paste0("$SUBROUTINE ", sub_line)
    }
    out <- c(out[1L], sub_line, out[-1L])
  }
  model_txt <- trimws(.nm_coalesce(parts$model, ""))
  if (!nzchar(model_txt) && advan != 10L &&
      (isTRUE(parts$use_ode) || nm_ctl_use_ode(advan))) {
    ncomp <- as.integer(.nm_coalesce(parts$ode_ncomp, 2L))
    model_txt <- sprintf("NCOMPS=%d", max(1L, ncomp))
  }
  if (nzchar(model_txt)) {
    if (!grepl("^\\$", model_txt)) {
      model_txt <- paste0("$MODEL ", model_txt)
    }
    sub_idx <- which(grepl("^\\$SUBROUTINE", out))
    if (length(sub_idx) == 1L) {
      out <- c(out[seq_len(sub_idx)], model_txt, out[-seq_len(sub_idx)])
    } else {
      out <- c(out[1L], model_txt, out[-1L])
    }
  }
  paste(out, collapse = "\n")
}

#' Canonical control text for stable dirty-state comparison
#'
#' @param ctl_text Control stream text.
#' @return Normalized control text.
#' @examples
#' parts <- nm_ctl_template(2L, 1L)
#' ctl <- nm_ctl_compose(parts)
#' nm_ctl_canonical(ctl)
#' @export
nm_ctl_canonical <- function(ctl_text) {
  if (is.null(ctl_text) || !nzchar(trimws(ctl_text))) {
    return("")
  }
  nm_ctl_compose(nm_ctl_parse(ctl_text))
}

#' Read column names from a dataset file
#'
#' @param path Path to CSV or similar table file.
#' @return Character vector of column names (uppercase).
#' @examples
#' sim <- nm_synthetic_theo(n_sub = 2L)
#' csv <- tempfile(fileext = ".csv")
#' write.csv(sim$data$data, csv, row.names = FALSE)
#' nm_ctl_read_columns(csv)
#' @export
nm_ctl_read_columns <- function(path) {
  if (!file.exists(path)) {
    return(character())
  }
  if (requireNamespace("data.table", quietly = TRUE)) {
    hdr <- names(data.table::fread(path, nrows = 0L, showProgress = FALSE))
  } else {
    hdr <- names(utils::read.csv(path, nrows = 0L, check.names = FALSE))
  }
  toupper(trimws(hdr))
}

#' Valid ADVAN values for the GUI (standard PREDPP models)
#'
#' @return Character vector of ADVAN numbers.
#' @examples
#' nm_ctl_advan_choices()
#' @export
nm_ctl_advan_choices <- function() {
  c("1", "2", "3", "4", "6", "10", "11", "12", "13")
}

#' Valid TRANS values for a given ADVAN (NONMEM PREDPP pairings)
#'
#' @param advan ADVAN number.
#' @return Character vector of valid TRANS values.
#' @examples
#' nm_ctl_trans_choices(2L)
#' @export
nm_ctl_trans_choices <- function(advan) {
  advan <- as.character(as.integer(advan))
  if (as.integer(advan) == 6L) {
    return(character())
  }
  m <- nm_ctl_trans_map()
  if (advan %in% names(m)) {
    m[[advan]]
  } else {
    character()
  }
}

#' Default TRANS for an ADVAN (CL/V-style parameterization when available)
#'
#' @examples
#' nm_ctl_default_trans(4L)
#' @export
nm_ctl_default_trans <- function(advan) {
  advan <- as.character(as.integer(advan))
  if (as.integer(advan) == 6L) {
    return("1")
  }
  choices <- nm_ctl_trans_choices(advan)
  if (length(choices) == 0L) {
    return("2")
  }
  prefer <- switch(
    advan,
    "1" = "2", "2" = "2", "3" = "4", "4" = "4",
    "10" = "1", "11" = "4", "12" = "4", "13" = "1",
    choices[[1L]]
  )
  if (prefer %in% choices) prefer else choices[[1L]]
}

#' Check whether an ADVAN/TRANS pair is valid in NONMEM
#'
#' @examples
#' nm_ctl_is_valid_pair(4L, 4L)
#' @export
nm_ctl_is_valid_pair <- function(advan, trans) {
  if (as.integer(advan) == 6L) {
    return(TRUE)
  }
  as.character(as.integer(trans)) %in% nm_ctl_trans_choices(advan)
}

#' Default \code{$DES} text for ADVAN 6/13 templates
#'
#' @param advan ADVAN number.
#' @param trans TRANS number (ADVAN 6 uses TRANS 1 internally).
#' @param ncomp Number of compartments for ADVAN 6 (1–10).
#' @param oral When \code{TRUE}, first compartment is an absorption depot.
#' @return Character scalar (newline-separated ODE lines).
#' @examples
#' nm_ctl_default_des(6L, 1L, ncomp = 2L, oral = TRUE)
#' @export
nm_ctl_default_des <- function(advan = 6L, trans = 1L, ncomp = 2L, oral = TRUE) {
  advan <- as.integer(advan)
  trans <- as.integer(trans)
  ncomp <- max(1L, min(10L, as.integer(ncomp)))
  if (advan == 13L) {
    return(paste(
      "DADT(1) = -KA*A(1)",
      "DADT(2) = KA*A(1) - CL/V*A(2)",
      "F = A(2)/S2",
      sep = "\n"
    ))
  }
  lines <- character()
  if (isTRUE(oral) && ncomp >= 2L) {
    lines <- c(lines, "DADT(1) = -KA*A(1)")
    lines <- c(lines, "DADT(2) = KA*A(1) - (CL/V)*A(2)")
    if (ncomp >= 3L) {
      for (i in 3L:ncomp) {
        lines <- c(lines, sprintf("DADT(%d) = 0", i))
      }
    }
    lines <- c(lines, "F = A(2)/S2")
  } else if (ncomp == 1L) {
    lines <- c(lines, "DADT(1) = -(CL/V)*A(1)", "F = A(1)/S1")
  } else {
    lines <- c(lines, "DADT(1) = -(CL/V)*A(1)")
    if (ncomp >= 2L) {
      for (i in 2L:ncomp) {
        lines <- c(lines, sprintf("DADT(%d) = 0", i))
      }
    }
    lines <- c(lines, sprintf("F = A(%d)/S%d", ncomp, ncomp))
  }
  paste(lines, collapse = "\n")
}

#' Whether the GUI should show TRANS selection for an ADVAN
#'
#' ADVAN 6 uses a fixed TRANS in \code{$SUBROUTINE} and has no user-facing TRANS choice.
#'
#' @examples
#' nm_ctl_show_trans(6L)
#' @export
nm_ctl_show_trans <- function(advan) {
  as.integer(advan) != 6L
}

#' Default observation compartment for ADVAN/TRANS
#'
#' @examples
#' nm_ctl_default_obscmp(4L, 4L)
#' @export
nm_ctl_default_obscmp <- function(advan, trans = NULL) {
  advan <- as.integer(advan)
  trans <- as.integer(trans %||% nm_ctl_default_trans(advan))
  if (advan == 4L && trans == 5L) {
    return(3L)
  }
  if (advan %in% c(2L, 4L, 6L, 12L, 13L)) {
    return(2L)
  }
  if (advan == 4L && trans %in% c(5L, 6L)) {
    return(2L)
  }
  1L
}

#' Effective TRANS value passed to the backend
#'
#' @examples
#' nm_ctl_effective_trans(4L, 4L)
#' @export
nm_ctl_effective_trans <- function(advan, trans = NULL) {
  advan <- as.integer(advan)
  if (advan == 6L) {
    return(1L)
  }
  if (is.null(trans)) {
    return(as.integer(nm_ctl_default_trans(advan)))
  }
  as.integer(trans)
}

#' Whether ADVAN uses the general ODE (\code{$SUBROUTINE}) path
#'
#' @examples
#' nm_ctl_use_ode(6L)
#' @export
nm_ctl_use_ode <- function(advan) {
  as.integer(advan) %in% c(6L, 10L, 13L)
}

#' $SUBROUTINE text for ODE ADVAN models
#'
#' @examples
#' nm_ctl_subroutine_text(4L, 4L)
#' @export
nm_ctl_subroutine_text <- function(advan, trans) {
  advan <- as.integer(advan)
  trans <- as.integer(trans)
  if (advan == 13L) {
    return("ADVAN13 TRANS1 TOL=9")
  }
  if (advan == 6L) {
    return("ADVAN6 TRANS1 TOL=9")
  }
  if (advan == 10L) {
    return("ADVAN10 TRANS1 TOL=9")
  }
  if (advan %in% c(1L, 2L, 3L, 4L, 11L, 12L)) {
    if (is.na(trans)) {
      trans <- as.integer(nm_ctl_default_trans(advan))
    }
    trans <- nm_ctl_effective_trans(advan, trans)
    return(sprintf("ADVAN%d TRANS%d", advan, trans))
  }
  ""
}

#' Required dataset columns for \code{$INPUT} given ADVAN/TRANS
#'
#' @param advan ADVAN number; when \code{NULL} returns core identifiers only.
#' @param trans TRANS number.
#' @return Character vector of required column names.
#' @examples
#' nm_ctl_essential_input_cols(4L, 4L)
#' @export
nm_ctl_essential_input_cols <- function(advan = NULL, trans = NULL) {
  base <- c("ID", "TIME", "EVID", "CMT", "MDV", "DV", "AMT")
  if (is.null(advan)) {
    return(base)
  }
  advan <- as.integer(advan)
  if (is.null(trans)) {
    trans <- as.integer(nm_ctl_default_trans(advan))
  } else {
    trans <- as.integer(trans)
  }
  extra <- switch(
    as.character(advan),
    "1" = c("RATE", "F1", "S1"),
    "2" = c("F1", "S1", "S2", "RATE"),
    "3" = c("F1", "F2", "S1", "S2", "RATE"),
    "4" = c("F1", "F2", "S1", "S2", "KA", "RATE"),
    "6" = c("F1", "F2", "S1", "S2", "KA", "RATE"),
    "10" = c("F1", "S1"),
    "11" = c("F1", "F2", "F3", "S1", "S2", "S3", "RATE"),
    "12" = c("F1", "F2", "F3", "F4", "S1", "S2", "S3", "S4", "KA", "RATE"),
    "13" = c("F1", "F2", "S1", "S2", "KA", "RATE"),
    character()
  )
  unique(c(base, extra))
}

#' @keywords internal
nm_ctl_trans_map <- function() {
  list(
    "1" = c("1", "2"),
    "2" = c("1", "2"),
    "3" = c("1", "3", "4", "6"),
    "4" = c("1", "3", "4", "5", "6"),
    "6" = character(),
    "10" = c("1"),
    "11" = c("1", "4", "6"),
    "12" = c("1", "4", "6"),
    "13" = c("1", "4")
  )
}

#' Model structure and parameter help for ADVAN/TRANS pairs
#'
#' Used by the Shiny GUI help popup. Returns compartment layout and brief
#' parameter definitions for the selected PREDPP model.
#'
#' @param advan ADVAN number.
#' @param trans TRANS number (\code{NULL} uses the ADVAN default).
#' @return List with \code{title}, \code{summary}, \code{route},
#'   \code{compartments}, \code{flows}, and \code{parameters} (data frame).
#' @examples
#' nm_ctl_model_info(4L, 4L)
#' @export
nm_ctl_model_info <- function(advan, trans = NULL) {
  advan <- as.integer(advan)
  if (nm_ctl_use_ode(advan)) {
    trans <- 1L
  } else if (is.null(trans)) {
    trans <- as.integer(nm_ctl_default_trans(advan))
  } else {
    trans <- as.integer(trans)
  }
  key <- paste0("A", advan, "_T", trans)
  catalog <- list(
    A1_T1 = list(
      title = "ADVAN 1, TRANS 1 — IV bolus, 1-compartment",
      route = "Intravenous bolus",
      summary = "Single central compartment with first-order elimination.",
      compartments = c("Central (CMT 1, dose & obs)"),
      flows = list(),
      parameters = data.frame(
        symbol = c("CL", "V"),
        meaning = c("Clearance", "Central volume"),
        stringsAsFactors = FALSE
      )
    ),
    A1_T2 = list(
      title = "ADVAN 1, TRANS 2 — IV infusion, 1-compartment",
      route = "Intravenous infusion",
      summary = "Single central compartment; infusion rate handled via RATE on dosing records.",
      compartments = c("Central (CMT 1, dose & obs)"),
      flows = list(),
      parameters = data.frame(
        symbol = c("CL", "V"),
        meaning = c("Clearance", "Central volume"),
        stringsAsFactors = FALSE
      )
    ),
    A2_T1 = list(
      title = "ADVAN 2, TRANS 1 — Oral, 1-compartment",
      route = "Oral",
      summary = "First-order absorption into a single central compartment.",
      compartments = c("Depot / gut (CMT 1, dose)", "Central (CMT 2, obs)"),
      flows = list(list(from = 1L, to = 2L, label = "KA")),
      parameters = data.frame(
        symbol = c("CL", "V", "KA"),
        meaning = c("Clearance", "Central volume", "Absorption rate constant"),
        stringsAsFactors = FALSE
      )
    ),
    A2_T2 = list(
      title = "ADVAN 2, TRANS 2 — Oral with lag, 1-compartment",
      route = "Oral (lag time)",
      summary = "Absorption depot with optional lag (ALAG) before KA.",
      compartments = c("Depot / gut (CMT 1, dose)", "Central (CMT 2, obs)"),
      flows = list(list(from = 1L, to = 2L, label = "KA")),
      parameters = data.frame(
        symbol = c("CL", "V", "KA", "ALAG"),
        meaning = c("Clearance", "Central volume", "Absorption rate constant", "Absorption lag time"),
        stringsAsFactors = FALSE
      )
    ),
    A3_T1 = list(
      title = "ADVAN 3, TRANS 1 — IV bolus, 2-compartment (CL/V)",
      route = "Intravenous bolus",
      summary = "Central and peripheral compartments; macro CL, V1, V2, Q parameterization.",
      compartments = c("Central (CMT 1, dose & obs)", "Peripheral (CMT 2)"),
      flows = list(
        list(from = 1L, to = 2L, label = "Q"),
        list(from = 2L, to = 1L, label = "Q"),
        list(from = 1L, to = NULL, label = "CL")
      ),
      parameters = data.frame(
        symbol = c("CL", "V1", "V2", "Q"),
        meaning = c("Clearance", "Central volume", "Peripheral volume", "Inter-compartmental clearance"),
        stringsAsFactors = FALSE
      )
    ),
    A3_T3 = list(
      title = "ADVAN 3, TRANS 3 — IV bolus, 2-compartment (CL/V/V2/Q)",
      route = "Intravenous bolus",
      summary = "Micro-rate constants between central and peripheral compartments.",
      compartments = c("Central (CMT 1, dose & obs)", "Peripheral (CMT 2)"),
      flows = list(
        list(from = 1L, to = 2L, label = "K12"),
        list(from = 2L, to = 1L, label = "K21"),
        list(from = 1L, to = NULL, label = "K10")
      ),
      parameters = data.frame(
        symbol = c("K10", "K12", "K21"),
        meaning = c("Elimination from central", "Central to peripheral", "Peripheral to central"),
        stringsAsFactors = FALSE
      )
    ),
    A3_T4 = list(
      title = "ADVAN 3, TRANS 4 — IV bolus, 2-compartment (CL/VC/VP/Q)",
      route = "Intravenous bolus",
      summary = "Standard NONMEM macro parameterization: CL, VC, VP, Q2.",
      compartments = c("Central (CMT 1, dose & obs)", "Peripheral (CMT 2)"),
      flows = list(
        list(from = 1L, to = 2L, label = "Q2"),
        list(from = 2L, to = 1L, label = "Q2"),
        list(from = 1L, to = NULL, label = "CL")
      ),
      parameters = data.frame(
        symbol = c("CL", "VC", "VP", "Q2"),
        meaning = c("Clearance", "Central volume", "Peripheral volume", "Inter-compartmental clearance"),
        stringsAsFactors = FALSE
      )
    ),
    A3_T6 = list(
      title = "ADVAN 3, TRANS 6 — IV bolus, 2-compartment (K/V)",
      route = "Intravenous bolus",
      summary = "Elimination and distribution rate constants with volumes.",
      compartments = c("Central (CMT 1, dose & obs)", "Peripheral (CMT 2)"),
      flows = list(
        list(from = 1L, to = 2L, label = "K12"),
        list(from = 2L, to = 1L, label = "K21"),
        list(from = 1L, to = NULL, label = "K10")
      ),
      parameters = data.frame(
        symbol = c("K10", "K12", "K21", "V1", "V2"),
        meaning = c("Elimination rate", "Central to peripheral", "Peripheral to central", "Central volume", "Peripheral volume"),
        stringsAsFactors = FALSE
      )
    ),
    A4_T1 = list(
      title = "ADVAN 4, TRANS 1 — Oral, 2-compartment (CL/V)",
      route = "Oral",
      summary = "Oral absorption into a 2-compartment body model.",
      compartments = c("Depot (CMT 1, dose)", "Central (CMT 2, obs)", "Peripheral (CMT 3)"),
      flows = list(
        list(from = 1L, to = 2L, label = "KA"),
        list(from = 2L, to = 3L, label = "Q"),
        list(from = 3L, to = 2L, label = "Q"),
        list(from = 2L, to = NULL, label = "CL")
      ),
      parameters = data.frame(
        symbol = c("CL", "V1", "V2", "Q", "KA"),
        meaning = c("Clearance", "Central volume", "Peripheral volume", "Inter-compartmental clearance", "Absorption rate constant"),
        stringsAsFactors = FALSE
      )
    ),
    A4_T3 = list(
      title = "ADVAN 4, TRANS 3 — Oral, 2-compartment (K10/K12/K21)",
      route = "Oral",
      summary = "Oral route with micro-rate 2-compartment disposition.",
      compartments = c("Depot (CMT 1, dose)", "Central (CMT 2, obs)", "Peripheral (CMT 3)"),
      flows = list(
        list(from = 1L, to = 2L, label = "KA"),
        list(from = 2L, to = 3L, label = "K12"),
        list(from = 3L, to = 2L, label = "K21"),
        list(from = 2L, to = NULL, label = "K10")
      ),
      parameters = data.frame(
        symbol = c("KA", "K10", "K12", "K21"),
        meaning = c("Absorption rate constant", "Elimination from central", "Central to peripheral", "Peripheral to central"),
        stringsAsFactors = FALSE
      )
    ),
    A4_T4 = list(
      title = "ADVAN 4, TRANS 4 — Oral, 2-compartment (CL/VC/VP/Q2)",
      route = "Oral",
      summary = "Common oral 2-compartment model (e.g. THEO): depot, central, peripheral.",
      compartments = c("Depot (CMT 1, dose)", "Central (CMT 2, obs)", "Peripheral (CMT 3)"),
      flows = list(
        list(from = 1L, to = 2L, label = "KA"),
        list(from = 2L, to = 3L, label = "Q2"),
        list(from = 3L, to = 2L, label = "Q2"),
        list(from = 2L, to = NULL, label = "CL")
      ),
      parameters = data.frame(
        symbol = c("CL", "VC", "VP", "Q2", "KA"),
        meaning = c("Clearance", "Central volume", "Peripheral volume", "Inter-compartmental clearance", "Absorption rate constant"),
        stringsAsFactors = FALSE
      )
    ),
    A4_T5 = list(
      title = "ADVAN 4, TRANS 5 — Oral, 2-compartment (CL/VC/VP/Q2, obs CMT 3)",
      route = "Oral",
      summary = "Like TRANS 4 but observation compartment is peripheral (CMT 3).",
      compartments = c("Depot (CMT 1, dose)", "Central (CMT 2)", "Peripheral (CMT 3, obs)"),
      flows = list(
        list(from = 1L, to = 2L, label = "KA"),
        list(from = 2L, to = 3L, label = "Q2"),
        list(from = 3L, to = 2L, label = "Q2")
      ),
      parameters = data.frame(
        symbol = c("CL", "VC", "VP", "Q2", "KA"),
        meaning = c("Clearance", "Central volume", "Peripheral volume", "Inter-compartmental clearance", "Absorption rate constant"),
        stringsAsFactors = FALSE
      )
    ),
    A4_T6 = list(
      title = "ADVAN 4, TRANS 6 — Oral, 2-compartment (K/V)",
      route = "Oral",
      summary = "Oral absorption with micro-rate disposition and explicit volumes.",
      compartments = c("Depot (CMT 1, dose)", "Central (CMT 2, obs)", "Peripheral (CMT 3)"),
      flows = list(
        list(from = 1L, to = 2L, label = "KA"),
        list(from = 2L, to = 3L, label = "K12"),
        list(from = 3L, to = 2L, label = "K21"),
        list(from = 2L, to = NULL, label = "K10")
      ),
      parameters = data.frame(
        symbol = c("KA", "K10", "K12", "K21", "V1", "V2"),
        meaning = c("Absorption rate constant", "Elimination rate", "Central to peripheral", "Peripheral to central", "Central volume", "Peripheral volume"),
        stringsAsFactors = FALSE
      )
    ),
    A6_T1 = list(
      title = "ADVAN 6 — General ODE model",
      route = "User-defined ($DES)",
      summary = paste0(
        "Compartment structure and flows are defined in $DES. ",
        "Use CMT on dosing/observation records to target compartments."
      ),
      compartments = c("User-defined in $DES (e.g. depot, central, peripheral)"),
      flows = list(),
      parameters = data.frame(
        symbol = c("THETA(*)", "DADT(*)", "CMT"),
        meaning = c("Structural parameters in $PK", "ODE right-hand sides in $DES", "Dataset compartment index for dose/obs"),
        stringsAsFactors = FALSE
      )
    ),
    A10_T1 = list(
      title = "ADVAN 10 — Michaelis-Menten elimination",
      route = "Intravenous",
      summary = "Non-linear elimination with Vmax and Km; single central compartment.",
      compartments = c("Central (CMT 1, dose & obs)"),
      flows = list(list(from = 1L, to = NULL, label = "Vmax/Km")),
      parameters = data.frame(
        symbol = c("VMAX", "KM", "V"),
        meaning = c("Maximum elimination rate", "Michaelis constant", "Central volume"),
        stringsAsFactors = FALSE
      )
    ),
    A11_T1 = list(
      title = "ADVAN 11, TRANS 1 — IV bolus, 3-compartment (CL/V)",
      route = "Intravenous bolus",
      summary = "Three disposition compartments plus central elimination.",
      compartments = c("Central (CMT 1, dose & obs)", "Periph. 1 (CMT 2)", "Periph. 2 (CMT 3)"),
      flows = list(
        list(from = 1L, to = 2L, label = "Q3"),
        list(from = 2L, to = 1L, label = "Q3"),
        list(from = 1L, to = 3L, label = "Q4"),
        list(from = 3L, to = 1L, label = "Q4"),
        list(from = 1L, to = NULL, label = "CL")
      ),
      parameters = data.frame(
        symbol = c("CL", "V1", "V2", "V3", "Q3", "Q4"),
        meaning = c("Clearance", "Central volume", "Periph. 1 volume", "Periph. 2 volume", "Q to periph. 1", "Q to periph. 2"),
        stringsAsFactors = FALSE
      )
    ),
    A11_T4 = list(
      title = "ADVAN 11, TRANS 4 — IV bolus, 3-compartment (CL/VC/VP/VP2/Q2/Q3)",
      route = "Intravenous bolus",
      summary = "Macro CL/VC/VP parameterization with two peripheral compartments.",
      compartments = c("Central (CMT 1, dose & obs)", "Periph. 1 (CMT 2)", "Periph. 2 (CMT 3)"),
      flows = list(
        list(from = 1L, to = 2L, label = "Q2"),
        list(from = 2L, to = 1L, label = "Q2"),
        list(from = 1L, to = 3L, label = "Q3"),
        list(from = 3L, to = 1L, label = "Q3"),
        list(from = 1L, to = NULL, label = "CL")
      ),
      parameters = data.frame(
        symbol = c("CL", "VC", "VP", "Q2", "VP2", "Q3"),
        meaning = c("Clearance", "Central volume", "Periph. 1 volume", "Q to periph. 1", "Periph. 2 volume", "Q to periph. 2"),
        stringsAsFactors = FALSE
      )
    ),
    A11_T6 = list(
      title = "ADVAN 11, TRANS 6 — IV bolus, 3-compartment (K/V)",
      route = "Intravenous bolus",
      summary = "Micro-rate constants with explicit volumes for three compartments.",
      compartments = c("Central (CMT 1, dose & obs)", "Periph. 1 (CMT 2)", "Periph. 2 (CMT 3)"),
      flows = list(
        list(from = 1L, to = 2L, label = "K12"),
        list(from = 2L, to = 1L, label = "K21"),
        list(from = 1L, to = 3L, label = "K13"),
        list(from = 3L, to = 1L, label = "K31"),
        list(from = 1L, to = NULL, label = "K10")
      ),
      parameters = data.frame(
        symbol = c("K10", "K12", "K21", "K13", "K31", "V1", "V2", "V3"),
        meaning = c("Elimination from central", "Central to periph. 1", "Periph. 1 to central", "Central to periph. 2", "Periph. 2 to central", "Central volume", "Periph. 1 volume", "Periph. 2 volume"),
        stringsAsFactors = FALSE
      )
    ),
    A12_T1 = list(
      title = "ADVAN 12, TRANS 1 — Oral, 3-compartment (CL/V)",
      route = "Oral",
      summary = "Oral absorption into a three-compartment body model.",
      compartments = c("Depot (CMT 1, dose)", "Central (CMT 2, obs)", "Periph. 1 (CMT 3)", "Periph. 2 (CMT 4)"),
      flows = list(
        list(from = 1L, to = 2L, label = "KA"),
        list(from = 2L, to = 3L, label = "Q3"),
        list(from = 3L, to = 2L, label = "Q3"),
        list(from = 2L, to = 4L, label = "Q4"),
        list(from = 4L, to = 2L, label = "Q4"),
        list(from = 2L, to = NULL, label = "CL")
      ),
      parameters = data.frame(
        symbol = c("CL", "V1", "V2", "V3", "Q3", "Q4", "KA"),
        meaning = c("Clearance", "Central volume", "Periph. 1 volume", "Periph. 2 volume", "Q to periph. 1", "Q to periph. 2", "Absorption rate constant"),
        stringsAsFactors = FALSE
      )
    ),
    A12_T4 = list(
      title = "ADVAN 12, TRANS 4 — Oral, 3-compartment (CL/VC/VP/VP2/Q2/Q3)",
      route = "Oral",
      summary = "Oral route with central and two peripheral compartments.",
      compartments = c("Depot (CMT 1, dose)", "Central (CMT 2, obs)", "Periph. 1 (CMT 3)", "Periph. 2 (CMT 4)"),
      flows = list(
        list(from = 1L, to = 2L, label = "KA"),
        list(from = 2L, to = 3L, label = "Q2"),
        list(from = 3L, to = 2L, label = "Q2"),
        list(from = 2L, to = 4L, label = "Q3"),
        list(from = 4L, to = 2L, label = "Q3"),
        list(from = 2L, to = NULL, label = "CL")
      ),
      parameters = data.frame(
        symbol = c("CL", "VC", "VP", "Q2", "VP2", "Q3", "KA"),
        meaning = c("Clearance", "Central volume", "Periph. 1 volume", "Q to periph. 1", "Periph. 2 volume", "Q to periph. 2", "Absorption rate constant"),
        stringsAsFactors = FALSE
      )
    ),
    A12_T6 = list(
      title = "ADVAN 12, TRANS 6 — Oral, 3-compartment (K/V)",
      route = "Oral",
      summary = "Oral absorption with micro-rate three-compartment disposition.",
      compartments = c("Depot (CMT 1, dose)", "Central (CMT 2, obs)", "Periph. 1 (CMT 3)", "Periph. 2 (CMT 4)"),
      flows = list(
        list(from = 1L, to = 2L, label = "KA"),
        list(from = 2L, to = 3L, label = "K12"),
        list(from = 3L, to = 2L, label = "K21"),
        list(from = 2L, to = 4L, label = "K14"),
        list(from = 4L, to = 2L, label = "K41"),
        list(from = 2L, to = NULL, label = "K10")
      ),
      parameters = data.frame(
        symbol = c("KA", "K10", "K12", "K21", "K14", "K41", "V1", "V2", "V3"),
        meaning = c("Absorption rate constant", "Elimination from central", "Central to periph. 1", "Periph. 1 to central", "Central to periph. 2", "Periph. 2 to central", "Central volume", "Periph. 1 volume", "Periph. 2 volume"),
        stringsAsFactors = FALSE
      )
    ),
    A13_T1 = list(
      title = "ADVAN 13 — General ODE model (NM7+)",
      route = "User-defined ($DES)",
      summary = "Like ADVAN 6 but uses the ADVAN13 subroutine interface.",
      compartments = c("User-defined in $DES"),
      flows = list(),
      parameters = data.frame(
        symbol = c("THETA(*)", "DADT(*)", "CMT"),
        meaning = c("Structural parameters in $PK", "ODE right-hand sides in $DES", "Dataset compartment index for dose/obs"),
        stringsAsFactors = FALSE
      )
    )
  )
  info <- catalog[[key]]
  if (is.null(info)) {
    adv_lab <- paste0("ADVAN ", advan)
    tr_lab <- if (nm_ctl_show_trans(advan)) paste0(", TRANS ", trans) else ""
    info <- list(
      title = paste0(adv_lab, tr_lab),
      route = "See NONMEM PREDPP documentation",
      summary = "No detailed help entry for this pair yet. Check $PK parameter names in the template.",
      compartments = character(),
      flows = list(),
      parameters = data.frame(
        symbol = character(),
        meaning = character(),
        stringsAsFactors = FALSE
      )
    )
  }
  if (advan == 13L) {
    info <- catalog[["A13_T1"]]
  }
  if (advan == 6L) {
    info <- catalog[["A6_T1"]]
  }
  structure(
    c(info, list(advan = advan, trans = trans)),
    class = "nm_ctl_model_info"
  )
}

#' Common derived columns offered for \code{$OUTPUT} selection
#'
#' @return Character vector of column names.
#' @examples
#' nm_ctl_output_extra_cols()
#' @export
nm_ctl_output_extra_cols <- function() {
  c("IPRED", "PRED", "IWRES", "WRES", "RES", "CWRES", "CPRED", "CRES")
}

#' All selectable names for the column picker (dataset + output extras)
#'
#' @param dataset_cols Column names from the dataset file.
#' @return Character vector.
#' @examples
#' nm_ctl_picker_columns(c("ID", "TIME", "DV"))
#' @export
nm_ctl_picker_columns <- function(dataset_cols) {
  dataset_cols <- toupper(trimws(dataset_cols))
  dataset_cols <- dataset_cols[nzchar(dataset_cols)]
  extras <- nm_ctl_output_extra_cols()
  unique(c(dataset_cols, setdiff(extras, dataset_cols)))
}

#' @keywords internal
.nm_ctl_normalize_thetas <- function(thetas) {
  if (is.null(thetas) || nrow(thetas) == 0L) {
    return(data.frame(
      THETA = 1L, Value = 1, FIX = FALSE, Lower = 0, Upper = NA_real_,
      Label = "", stringsAsFactors = FALSE
    ))
  }
  if (!"FIX" %in% names(thetas)) {
    thetas$FIX <- FALSE
  }
  if (!"Lower" %in% names(thetas)) {
    thetas$Lower <- 0
  }
  if (!"Upper" %in% names(thetas)) {
    thetas$Upper <- NA_real_
  }
  if (!"Label" %in% names(thetas)) {
    thetas$Label <- ""
  }
  for (i in seq_len(nrow(thetas))) {
    val <- thetas$Value[[i]]
    bd <- .nm_theta_default_bounds(val)
    lo <- thetas$Lower[[i]]
    up <- thetas$Upper[[i]]
    if (length(lo) != 1L || is.na(lo)) {
      thetas$Lower[[i]] <- bd[["lower"]]
    }
    if (length(up) != 1L || is.na(up) || !is.finite(up)) {
      thetas$Upper[[i]] <- bd[["upper"]]
    }
  }
  thetas
}

#' @keywords internal
.nm_ctl_coalesce <- function(x, y) if (is.null(x)) y else x

#' @keywords internal
.nm_ctl_template_thetas <- function(values, labels = NULL) {
  nm <- names(values)
  values <- as.numeric(values)
  n <- length(values)
  if (is.null(labels) || length(labels) != n) {
    if (length(nm) == n) {
      labels <- as.character(nm)
    } else {
      labels <- paste0("THETA", seq_len(n))
    }
  } else {
    labels <- as.character(labels)
  }
  lowers <- vapply(values, function(v) .nm_theta_default_bounds(v)[["lower"]], numeric(1))
  uppers <- vapply(values, function(v) .nm_theta_default_bounds(v)[["upper"]], numeric(1))
  data.frame(
    THETA = seq_len(n),
    Value = values,
    FIX = rep(FALSE, n),
    Lower = lowers,
    Upper = uppers,
    Label = labels,
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.nm_ctl_template_omegas <- function(n, value = 0.09) {
  data.frame(
    OMEGA = seq_len(n),
    Value = rep(as.numeric(value), n),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.nm_ctl_template_sigmas <- function(values) {
  values <- as.numeric(values)
  data.frame(
    SIGMA = seq_along(values),
    Value = values,
    stringsAsFactors = FALSE
  )
}

#' Parameter name aliases for NM-TRAN (macro models use Q not Q2, etc.)
#' @keywords internal
.nm_ctl_pk_trans6_cardano_roots <- function() {
  c(
    "IF (TR+TD .GE. 0) THEN",
    " C1TR = (TR+TD)**0.333333333333",
    "ELSE",
    " C1TR = -((-TR-TD)**0.333333333333)",
    "ENDIF",
    "IF (TR-TD .GE. 0) THEN",
    " C2TR = (TR-TD)**0.333333333333",
    "ELSE",
    " C2TR = -((-TR+TD)**0.333333333333)",
    "ENDIF",
    "ALPHA = -TSUM/3 + (C1TR + C2TR)/2",
    "BETA = -TSUM/3 + (C1TR*(-0.5 + 0.866025403784) + C2TR*(-0.5 - 0.866025403784))/2",
    "GAMMA = -TSUM/3 + (C1TR*(-0.5 - 0.866025403784) + C2TR*(-0.5 + 0.866025403784))/2"
  )
}

.nm_ctl_pk_trans6_hybrid_lines <- function(advan, oral = FALSE) {
  advan <- as.integer(advan)
  oral <- isTRUE(oral)
  cardano <- .nm_ctl_pk_trans6_cardano_roots()
  if (advan %in% c(11L, 12L)) {
    if (advan == 12L && oral) {
      return(c(
        "E1 = K10 + K23 + K24",
        "E2 = K32 + K42",
        "E3 = K23 + K24",
        "TSUM = E1 + E2 + E3",
        "TB = E1*E2 + E1*E3 + E2*E3 - K23*K32 - K24*K42",
        "TC = E1*E2*E3 - E3*K23*K32 - E2*K24*K42",
        "TQ = TB/3 - TSUM*TSUM/9",
        "TR = (TSUM*TSUM*TSUM - 9*TSUM*TB + 27*TC)/108",
        "TD = (TQ*TQ*TQ + TR*TR)**0.5",
        cardano
      ))
    }
    return(c(
      "E1 = K10 + K12 + K13",
      "E2 = K21 + K31",
      "E3 = K12 + K13",
      "TSUM = E1 + E2 + E3",
      "TB = E1*E2 + E1*E3 + E2*E3 - K12*K21 - K13*K31",
      "TC = E1*E2*E3 - E3*K12*K21 - E2*K13*K31",
      "TQ = TB/3 - TSUM*TSUM/9",
      "TR = (TSUM*TSUM*TSUM - 9*TSUM*TB + 27*TC)/108",
      "TD = (TQ*TQ*TQ + TR*TR)**0.5",
      cardano
    ))
  }
  if (oral && advan == 4L) {
    return(c(
      "E1 = K10 + K23",
      "E2 = K32",
      "DISC = (E1-E2)*(E1-E2) + 4*K23*K32",
      "ALPHA = 0.5*(E1 + E2 + (DISC)**0.5)",
      "BETA = 0.5*(E1 + E2 - (DISC)**0.5)"
    ))
  }
  c(
    "E1 = K10 + K12",
    "E2 = K21",
    "DISC = (E1-E2)*(E1-E2) + 4*K12*K21",
    "ALPHA = 0.5*(E1 + E2 + (DISC)**0.5)",
    "BETA = 0.5*(E1 + E2 - (DISC)**0.5)"
  )
}

#' @keywords internal
.nm_ctl_pk_nm_param_aliases <- function(advan, trans) {
  advan <- as.integer(advan)
  trans <- as.integer(trans)
  lines <- character()
  if (advan == 3L && trans == 3L) {
    lines <- c(lines, "VSS = V + V2")
  }
  if (advan == 4L && trans == 3L) {
    lines <- c(lines, "VSS = V + V2")
  }
  if (advan %in% c(3L, 4L) && trans %in% c(4L, 2L, 5L)) {
    lines <- c(lines, "Q = Q2")
  }
  if (advan == 11L && trans == 4L) {
    lines <- c(lines, "Q = Q2")
  }
  if (advan == 3L && trans %in% c(1L, 6L)) {
    lines <- c(lines, "K = K10")
  }
  if (advan == 4L && trans %in% c(1L, 6L)) {
    lines <- c(lines, "K = K10")
  }
  if (advan == 11L && trans %in% c(1L, 6L)) {
    lines <- c(lines, "K = K10")
  }
  if (advan == 12L && trans %in% c(1L, 6L)) {
    lines <- c(lines, "K = K10")
  }
  if (advan == 4L && trans == 5L) {
    lines <- c(
      lines,
      "K10 = CL/VC",
      "K12 = Q2/VC",
      "K21 = Q2/VP",
      "K = CL/VC",
      "AOB = 1",
      .nm_ctl_pk_trans6_hybrid_lines(advan)
    )
  }
  if (advan == 12L && trans == 4L) {
    lines <- c(lines, "Q = Q2", "Q4 = Q3")
  }
  if (trans == 6L) {
    if (advan == 4L) {
      lines <- c(lines, .nm_ctl_pk_trans6_hybrid_lines(advan, oral = TRUE))
    } else if (advan == 3L) {
      lines <- c(lines, .nm_ctl_pk_trans6_hybrid_lines(advan, oral = FALSE))
    }
  }
  lines
}

#' Volume aliases required by NM-TRAN before Sx assignments (macro CL/VC/VP models)
#' @keywords internal
.nm_ctl_pk_volume_aliases <- function(advan, trans) {
  advan <- as.integer(advan)
  trans <- as.integer(trans)
  oral <- advan %in% c(2L, 4L, 12L, 13L) || nm_ctl_use_ode(advan)
  if (advan == 3L && trans %in% c(4L, 2L)) {
    return(c("V1 = VC", "V2 = VP"))
  }
  if (advan == 4L && trans %in% c(4L, 5L)) {
    return(c("V2 = VC", "V3 = VP"))
  }
  if (advan == 11L && trans == 4L) {
    return(c("V1 = VC", "V2 = VP", "V3 = VP2"))
  }
  if (advan == 12L && trans == 4L) {
    return(c("V2 = VC", "V3 = VP", "V4 = VP2"))
  }
  if (advan == 11L && trans == 2L) {
    return(character())
  }
  if (advan == 3L && trans == 2L) {
    return(character())
  }
  if (advan == 4L && !trans %in% c(4L, 5L, 6L, 1L, 3L)) {
    return(character())
  }
  character()
}

#' Compartment scaling lines for $PK (NONMEM Sx convention: F = A(obs)/S(obs))
#' @keywords internal
.nm_ctl_pk_scaling <- function(advan, trans, ode_ncomp = 2L) {
  advan <- as.integer(advan)
  trans <- as.integer(trans)
  if (advan == 10L) {
    return(character())
  }
  if (advan == 1L) {
    return("S1 = V")
  }
  if (advan == 2L) {
    return(c("S1 = 1", "S2 = V"))
  }
  if (advan == 3L) {
    if (trans %in% c(4L, 2L)) {
      if (trans == 4L) {
        return(c("S1 = VC", "S2 = VP"))
      }
      return(c("S1 = V1", "S2 = V2"))
    }
    if (trans == 6L) {
      return(c("S1 = V1", "S2 = V2"))
    }
    if (trans %in% c(1L, 3L)) {
      return(c("S1 = 1", "S2 = 1"))
    }
    return(c("S1 = V1", "S2 = V2"))
  }
  if (advan == 4L) {
    if (trans %in% c(4L, 5L)) {
      return(c("S1 = 1", "S2 = VC", "S3 = VP"))
    }
    if (trans == 6L) {
      return(c("S1 = 1", "S2 = V1", "S3 = V2"))
    }
    if (trans == 3L) {
      return(c("S1 = 1", "S2 = V", "S3 = V2"))
    }
    if (trans == 1L) {
      return(c("S1 = 1", "S2 = 1"))
    }
    return(c("S1 = 1", "S2 = V2"))
  }
  if (advan == 11L) {
    if (trans == 4L) {
      return(c("S1 = VC", "S2 = VP", "S3 = VP2"))
    }
    if (trans == 6L) {
      return(c("S1 = V1", "S2 = V2", "S3 = V3"))
    }
    if (trans %in% c(1L, 3L)) {
      return(c("S1 = 1", "S2 = 1", "S3 = 1"))
    }
    return(c("S1 = V1", "S2 = V2", "S3 = V3"))
  }
  if (advan == 12L) {
    if (trans == 4L) {
      return(c("S1 = 1", "S2 = VC", "S3 = VP", "S4 = VP2"))
    }
    if (trans == 6L) {
      return(c("S1 = 1", "S2 = V1", "S3 = V2", "S4 = V3"))
    }
    if (trans %in% c(1L, 3L)) {
      return(c("S1 = 1", "S2 = 1", "S3 = 1", "S4 = 1"))
    }
    return(c("S1 = 1", "S2 = V1", "S3 = V2", "S4 = V3"))
  }
  if (nm_ctl_use_ode(advan)) {
    return(.nm_ctl_pk_scaling_ode(ode_ncomp, oral = TRUE, advan = advan))
  }
  character()
}

#' ODE model S1..Sn scaling for $PK (supports NCOMPS up to 10)
#' @keywords internal
.nm_ctl_pk_scaling_ode <- function(ncomp = 2L, oral = TRUE, advan = 6L) {
  ncomp <- max(1L, min(.nm_max_scale_n(), as.integer(ncomp)))
  advan <- as.integer(advan)
  if (advan == 13L) {
    return(c("S1 = 1", "S2 = V"))
  }
  if (isTRUE(oral)) {
    lines <- c("S1 = 1", "S2 = V")
    if (ncomp > 2L) {
      for (i in 3:ncomp) {
        lines <- c(lines, sprintf("S%d = 1", i))
      }
    }
    return(lines)
  }
  if (ncomp == 1L) {
    return("S1 = V")
  }
  lines <- vapply(seq_len(ncomp - 1L), function(i) sprintf("S%d = 1", i), character(1L))
  c(lines, sprintf("S%d = V", ncomp))
}

#' @keywords internal
.nm_ctl_pk_nonmem_extras <- function(advan, trans) {
  c(
    .nm_ctl_pk_volume_aliases(advan, trans),
    .nm_ctl_pk_nm_param_aliases(advan, trans)
  )
}

#' @keywords internal
.nm_ctl_append_pk_scaling <- function(pk_body, advan, trans, ode_ncomp = 2L) {
  scale <- .nm_ctl_pk_scaling(advan, trans, ode_ncomp = ode_ncomp)
  if (length(scale) == 0L) {
    return(pk_body)
  }
  paste(c(pk_body, scale), collapse = "\n")
}

#' @keywords internal
.nm_ctl_template_spec <- function(advan, trans) {
  advan <- as.integer(advan)
  trans <- as.integer(trans)
  prop_add_err <- "Y = F * (1 + ERR(1)) + ERR(2)"
  prop_err <- "Y = F * (1 + ERR(1))"
  prop_add_err_ode <- "Y = F * (1 + ERR(1)) + ERR(2)"
  prop_err_ode <- "Y = F * (1 + ERR(1))"

  pk_k_v <- function(n_eta = 1L) {
    paste(
      "K = THETA(1) * exp(ETA(1))",
      "V = THETA(2)",
      sep = "\n"
    )
  }
  pk_k_v_ka <- function(n_eta = 2L) {
    paste(
      "K = THETA(1) * exp(ETA(1))",
      "V = THETA(2)",
      "KA = THETA(3) * exp(ETA(2))",
      sep = "\n"
    )
  }
  pk_cl_v <- function(n_eta = 1L) {
    paste(
      "CL = THETA(1) * exp(ETA(1))",
      "V = THETA(2)",
      sep = "\n"
    )
  }
  pk_cl_v_ka <- function(n_eta = 2L) {
    paste(
      "CL = THETA(1) * exp(ETA(1))",
      "V = THETA(2)",
      "KA = THETA(3) * exp(ETA(2))",
      sep = "\n"
    )
  }
  pk_cl_v1v2q <- function(n_eta = 2L) {
    paste(
      "CL = THETA(1) * exp(ETA(1))",
      "V1 = THETA(2) * exp(ETA(2))",
      "V2 = THETA(3)",
      "Q = THETA(4)",
      sep = "\n"
    )
  }
  pk_cl_vc_vp_q2 <- function(n_eta = 2L) {
    paste(
      "CL = THETA(1) * exp(ETA(1))",
      "VC = THETA(2) * exp(ETA(2))",
      "VP = THETA(3)",
      "Q2 = THETA(4)",
      sep = "\n"
    )
  }
  pk_micro_2c <- function() {
    paste(
      "K10 = THETA(1) * exp(ETA(1))",
      "K12 = THETA(2)",
      "K21 = THETA(3)",
      sep = "\n"
    )
  }
  pk_micro_2c_v <- function() {
    paste(
      "K10 = THETA(1) * exp(ETA(1))",
      "K12 = THETA(2)",
      "K21 = THETA(3)",
      "V1 = THETA(4) * exp(ETA(2))",
      "V2 = THETA(5)",
      sep = "\n"
    )
  }
  pk_micro_3c <- function() {
    paste(
      "K10 = THETA(1) * exp(ETA(1))",
      "K12 = THETA(2)",
      "K21 = THETA(3)",
      "K13 = THETA(4)",
      "K31 = THETA(5)",
      sep = "\n"
    )
  }
  pk_oral_micro_2c_adv4 <- function(n_eta = 3L) {
    paste(
      "KA = THETA(1) * exp(ETA(1))",
      "K10 = THETA(2) * exp(ETA(2))",
      "K23 = THETA(3)",
      "K32 = THETA(4)",
      sep = "\n"
    )
  }
  pk_oral_micro_4c <- function(n_eta = 2L) {
    paste(
      "KA = THETA(1) * exp(ETA(1))",
      "K10 = THETA(2) * exp(ETA(2))",
      "K12 = THETA(3)",
      "K21 = THETA(4)",
      "K23 = THETA(5)",
      "K32 = THETA(6)",
      sep = "\n"
    )
  }
  pk_oral_micro_4c_adv12 <- function(n_eta = 2L) {
    paste(
      "KA = THETA(1) * exp(ETA(1))",
      "K10 = THETA(2) * exp(ETA(2))",
      "K23 = THETA(3)",
      "K32 = THETA(4)",
      "K24 = THETA(5)",
      "K42 = THETA(6)",
      sep = "\n"
    )
  }
  pk_oral_micro_2c_hub_v <- function(n_eta = 3L) {
    paste(
      "KA = THETA(1) * exp(ETA(1))",
      "K10 = THETA(2) * exp(ETA(2))",
      "K23 = THETA(3)",
      "K32 = THETA(4)",
      "V1 = THETA(5) * exp(ETA(3))",
      "V2 = THETA(6)",
      sep = "\n"
    )
  }
  pk_oral_cl_v_v2q <- function(n_eta = 3L) {
    paste(
      "CL = THETA(1) * exp(ETA(1))",
      "V = THETA(2) * exp(ETA(2))",
      "V2 = THETA(3)",
      "Q = THETA(4)",
      "KA = THETA(5) * exp(ETA(3))",
      sep = "\n"
    )
  }
  pk_oral_micro_3c <- function(n_eta = 2L) {
    paste(
      "KA = THETA(1) * exp(ETA(1))",
      "K10 = THETA(2) * exp(ETA(2))",
      "K12 = THETA(3)",
      "K21 = THETA(4)",
      "K13 = THETA(5)",
      "K31 = THETA(6)",
      sep = "\n"
    )
  }
  pk_oral_cl_v1v2q <- function(n_eta = 3L) {
    paste(
      pk_cl_v1v2q(n_eta),
      "KA = THETA(5) * exp(ETA(3))",
      sep = "\n"
    )
  }
  pk_oral_cl_vc_vp_q2 <- function(n_eta = 3L) {
    paste(
      pk_cl_vc_vp_q2(n_eta),
      "KA = THETA(5) * exp(ETA(3))",
      sep = "\n"
    )
  }
  pk_oral_micro_2c <- function(n_eta = 3L) {
    paste(
      "KA = THETA(1) * exp(ETA(1))",
      "K10 = THETA(2) * exp(ETA(2))",
      "K12 = THETA(3)",
      "K21 = THETA(4)",
      sep = "\n"
    )
  }
  pk_oral_micro_2c_v <- function(n_eta = 3L) {
    paste(
      "KA = THETA(1) * exp(ETA(1))",
      "K10 = THETA(2) * exp(ETA(2))",
      "K12 = THETA(3)",
      "K21 = THETA(4)",
      "V1 = THETA(5) * exp(ETA(3))",
      "V2 = THETA(6)",
      sep = "\n"
    )
  }
  pk_3c_cl_v <- function(n_eta = 2L) {
    paste(
      "CL = THETA(1) * exp(ETA(1))",
      "V1 = THETA(2) * exp(ETA(2))",
      "V2 = THETA(3)",
      "V3 = THETA(4)",
      "Q3 = THETA(5)",
      "Q4 = THETA(6)",
      sep = "\n"
    )
  }
  pk_3c_cl_vc <- function(n_eta = 2L) {
    paste(
      "CL = THETA(1) * exp(ETA(1))",
      "VC = THETA(2) * exp(ETA(2))",
      "VP = THETA(3)",
      "Q2 = THETA(4)",
      "VP2 = THETA(5)",
      "Q3 = THETA(6)",
      sep = "\n"
    )
  }
  pk_oral_3c_cl_v <- function(n_eta = 3L) {
    paste(
      pk_3c_cl_v(n_eta),
      "KA = THETA(7) * exp(ETA(3))",
      sep = "\n"
    )
  }
  pk_oral_3c_cl_vc <- function(n_eta = 3L) {
    paste(
      pk_3c_cl_vc(n_eta),
      "KA = THETA(7) * exp(ETA(3))",
      sep = "\n"
    )
  }

  if (advan == 1L) {
    pk <- if (trans == 1L) pk_k_v() else pk_cl_v()
    return(list(
      thetas = .nm_ctl_template_thetas(c(if (trans == 1L) c(K = 0.1, V = 10) else c(CL = 1, V = 10))),
      omegas = .nm_ctl_template_omegas(1L),
      sigmas = .nm_ctl_template_sigmas(0.1),
      pk = pk,
      error = prop_err
    ))
  }
  if (advan == 2L && trans == 2L) {
    return(list(
      thetas = .nm_ctl_template_thetas(c(CL = 1, V = 10, KA = 1, ALAG = 0.5)),
      omegas = .nm_ctl_template_omegas(2L),
      sigmas = .nm_ctl_template_sigmas(c(0.1, 0.5)),
      pk = paste(
        pk_cl_v_ka(),
        "ALAG = THETA(4)",
        sep = "\n"
      ),
      error = prop_add_err
    ))
  }
  if (advan == 2L) {
    pk <- if (trans == 1L) pk_k_v_ka() else pk_cl_v_ka()
    return(list(
      thetas = .nm_ctl_template_thetas(
        if (trans == 1L) c(K = 0.1, V = 10, KA = 1) else c(CL = 1, V = 10, KA = 1)
      ),
      omegas = .nm_ctl_template_omegas(2L),
      sigmas = .nm_ctl_template_sigmas(c(0.1, 0.5)),
      pk = pk,
      error = prop_add_err
    ))
  }
  if (advan == 3L && trans == 1L) {
    return(list(
      thetas = .nm_ctl_template_thetas(c(K10 = 0.1, K12 = 0.5, K21 = 0.3)),
      omegas = .nm_ctl_template_omegas(1L),
      sigmas = .nm_ctl_template_sigmas(0.1),
      pk = pk_micro_2c(),
      error = prop_err
    ))
  }
  if (advan == 3L && trans == 3L) {
    return(list(
      thetas = .nm_ctl_template_thetas(c(CL = 1, V = 10, V2 = 20, Q = 2)),
      omegas = .nm_ctl_template_omegas(2L),
      sigmas = .nm_ctl_template_sigmas(0.1),
      pk = paste(
        "CL = THETA(1) * exp(ETA(1))",
        "V = THETA(2) * exp(ETA(2))",
        "V2 = THETA(3)",
        "Q = THETA(4)",
        sep = "\n"
      ),
      error = prop_err
    ))
  }
  if (advan == 3L && trans == 4L) {
    return(list(
      thetas = .nm_ctl_template_thetas(c(CL = 3, VC = 20, VP = 50, Q2 = 10)),
      omegas = .nm_ctl_template_omegas(2L),
      sigmas = .nm_ctl_template_sigmas(0.1),
      pk = pk_cl_vc_vp_q2(),
      error = prop_err
    ))
  }
  if (advan == 3L && trans == 6L) {
    return(list(
      thetas = .nm_ctl_template_thetas(c(K10 = 0.1, K12 = 0.5, K21 = 0.3, V1 = 10, V2 = 20)),
      omegas = .nm_ctl_template_omegas(2L),
      sigmas = .nm_ctl_template_sigmas(0.1),
      pk = pk_micro_2c_v(),
      error = prop_err
    ))
  }
  if (advan == 3L) {
    return(list(
      thetas = .nm_ctl_template_thetas(c(CL = 1, V1 = 10, V2 = 20, Q = 2)),
      omegas = .nm_ctl_template_omegas(2L),
      sigmas = .nm_ctl_template_sigmas(0.1),
      pk = pk_cl_v1v2q(),
      error = prop_err
    ))
  }
  if (advan == 4L && trans == 1L) {
    return(list(
      thetas = .nm_ctl_template_thetas(c(KA = 1.2, K10 = 0.1, K23 = 0.5, K32 = 0.3)),
      omegas = .nm_ctl_template_omegas(3L),
      sigmas = .nm_ctl_template_sigmas(c(0.1, 0.5)),
      pk = pk_oral_micro_2c_adv4(),
      error = prop_add_err
    ))
  }
  if (advan == 4L && trans == 3L) {
    return(list(
      thetas = .nm_ctl_template_thetas(c(CL = 1, V = 10, V2 = 20, Q = 2, KA = 1.2)),
      omegas = .nm_ctl_template_omegas(3L),
      sigmas = .nm_ctl_template_sigmas(c(0.1, 0.5)),
      pk = pk_oral_cl_v_v2q(),
      error = prop_add_err
    ))
  }
  if (advan == 4L && trans == 4L) {
    return(list(
      thetas = .nm_ctl_template_thetas(c(CL = 3, VC = 20, VP = 50, Q2 = 10, KA = 1.2)),
      omegas = .nm_ctl_template_omegas(3L),
      sigmas = .nm_ctl_template_sigmas(c(0.1, 0.5)),
      pk = pk_oral_cl_vc_vp_q2(),
      error = prop_add_err
    ))
  }
  if (advan == 4L && trans %in% c(5L, 6L)) {
    pk <- if (trans == 6L) pk_oral_micro_2c_hub_v() else pk_oral_cl_vc_vp_q2()
    th <- if (trans == 6L) {
      c(KA = 1, K10 = 0.1, K23 = 0.5, K32 = 0.3, V1 = 10, V2 = 20)
    } else {
      c(CL = 1, VC = 10, VP = 20, Q2 = 2, KA = 1)
    }
    return(list(
      thetas = .nm_ctl_template_thetas(th),
      omegas = .nm_ctl_template_omegas(3L),
      sigmas = .nm_ctl_template_sigmas(c(0.1, 0.5)),
      pk = pk,
      error = prop_add_err
    ))
  }
  if (advan == 4L) {
    return(list(
      thetas = .nm_ctl_template_thetas(c(CL = 1, V1 = 10, V2 = 20, Q = 2, KA = 1)),
      omegas = .nm_ctl_template_omegas(3L),
      sigmas = .nm_ctl_template_sigmas(c(0.1, 0.5)),
      pk = pk_oral_cl_v1v2q(),
      error = prop_add_err
    ))
  }
  if (advan == 11L && trans == 1L) {
    return(list(
      thetas = .nm_ctl_template_thetas(c(K10 = 0.1, K12 = 0.5, K21 = 0.3, K13 = 0.2, K31 = 0.15)),
      omegas = .nm_ctl_template_omegas(1L),
      sigmas = .nm_ctl_template_sigmas(0.1),
      pk = pk_micro_3c(),
      error = prop_err
    ))
  }
  if (advan == 11L && trans == 4L) {
    return(list(
      thetas = .nm_ctl_template_thetas(c(CL = 1, VC = 10, VP = 20, Q2 = 2, VP2 = 40, Q3 = 1)),
      omegas = .nm_ctl_template_omegas(2L),
      sigmas = .nm_ctl_template_sigmas(0.1),
      pk = pk_3c_cl_vc(),
      error = prop_err
    ))
  }
  if (advan == 11L && trans == 6L) {
    return(list(
      thetas = .nm_ctl_template_thetas(c(K10 = 0.1, K12 = 0.5, K21 = 0.3, K13 = 0.2, K31 = 0.15, V1 = 10, V2 = 20, V3 = 40)),
      omegas = .nm_ctl_template_omegas(2L),
      sigmas = .nm_ctl_template_sigmas(0.1),
      pk = paste(
        "K10 = THETA(1) * exp(ETA(1))",
        "K12 = THETA(2)",
        "K21 = THETA(3)",
        "K13 = THETA(4)",
        "K31 = THETA(5)",
        "V1 = THETA(6) * exp(ETA(2))",
        "V2 = THETA(7)",
        "V3 = THETA(8)",
        "ALPHA = 0.45",
        "BETA = 0.20",
        "GAMMA = 0.08",
        sep = "\n"
      ),
      error = prop_err
    ))
  }
  if (advan == 11L) {
    return(list(
      thetas = .nm_ctl_template_thetas(c(CL = 1, V1 = 10, V2 = 20, V3 = 40, Q3 = 1, Q4 = 0.5)),
      omegas = .nm_ctl_template_omegas(2L),
      sigmas = .nm_ctl_template_sigmas(0.1),
      pk = pk_3c_cl_v(),
      error = prop_err
    ))
  }
  if (advan == 12L && trans == 1L) {
    return(list(
      thetas = .nm_ctl_template_thetas(c(KA = 0.8, K10 = 0.1, K23 = 0.2, K32 = 0.15, K24 = 0.1, K42 = 0.08)),
      omegas = .nm_ctl_template_omegas(2L),
      sigmas = .nm_ctl_template_sigmas(c(0.1, 0.5)),
      pk = pk_oral_micro_4c_adv12(),
      error = prop_add_err
    ))
  }
  if (advan == 12L && trans == 4L) {
    return(list(
      thetas = .nm_ctl_template_thetas(c(CL = 1, VC = 10, VP = 20, Q2 = 2, VP2 = 40, Q3 = 1, KA = 0.8)),
      omegas = .nm_ctl_template_omegas(3L),
      sigmas = .nm_ctl_template_sigmas(c(0.1, 0.5)),
      pk = pk_oral_3c_cl_vc(),
      error = prop_add_err
    ))
  }
  if (advan == 12L && trans == 6L) {
    return(list(
      thetas = .nm_ctl_template_thetas(c(KA = 0.8, K10 = 0.1, V1 = 10, V2 = 20, V3 = 40)),
      omegas = .nm_ctl_template_omegas(3L),
      sigmas = .nm_ctl_template_sigmas(c(0.1, 0.5)),
      pk = paste(
        "KA = THETA(1) * exp(ETA(1))",
        "K10 = THETA(2) * exp(ETA(2))",
        "V1 = THETA(3) * exp(ETA(3))",
        "V2 = THETA(4)",
        "V3 = THETA(5)",
        "ALPHA = 0.45",
        "BETA = 0.20",
        "GAMMA = 0.08",
        "K32 = 0.30",
        "K42 = 0.12",
        sep = "\n"
      ),
      error = prop_add_err
    ))
  }
  if (advan == 12L) {
    return(list(
      thetas = .nm_ctl_template_thetas(c(CL = 1, V1 = 10, V2 = 20, V3 = 40, Q3 = 1, Q4 = 0.5, KA = 0.8)),
      omegas = .nm_ctl_template_omegas(3L),
      sigmas = .nm_ctl_template_sigmas(c(0.1, 0.5)),
      pk = pk_oral_3c_cl_v(),
      error = prop_add_err
    ))
  }
  if (advan == 10L) {
    return(list(
      thetas = .nm_ctl_template_thetas(c(VM = 10, KM = 1)),
      omegas = .nm_ctl_template_omegas(1L),
      sigmas = .nm_ctl_template_sigmas(0.1),
      pk = paste(
        "VM = THETA(1) * exp(ETA(1))",
        "KM = THETA(2)",
        sep = "\n"
      ),
      des = "DADT(1) = -(VM*A(1)/(KM + A(1)))",
      error = "A1 = A(1)\nY = F * (1 + ERR(1))"
    ))
  }
  if (nm_ctl_use_ode(advan)) {
    return(list(
      thetas = .nm_ctl_template_thetas(c(CL = 1, V = 10, KA = 1)),
      omegas = .nm_ctl_template_omegas(2L),
      sigmas = .nm_ctl_template_sigmas(c(0.1, 0.5)),
      pk = paste(
        "CL = THETA(1) * exp(ETA(1))",
        "V = THETA(2)",
        "KA = THETA(3) * exp(ETA(2))",
        sep = "\n"
      ),
      des = nm_ctl_default_des(advan, 1L, ncomp = 2L, oral = TRUE),
      error = prop_add_err_ode
    ))
  }
  list(
    thetas = .nm_ctl_template_thetas(c(CL = 1, V = 10)),
    omegas = .nm_ctl_template_omegas(1L),
    sigmas = .nm_ctl_template_sigmas(0.1),
    pk = pk_cl_v(),
    error = prop_err
  )
}

#' Update control-stream THETA / OMEGA / SIGMA initials from a fit
#'
#' @param parts List as returned by \code{\link{nm_ctl_parse}}.
#' @param fit An \code{nm_fit} object with \code{par} and \code{model}.
#' @return Updated \code{parts} list.
#' @examples
#' \dontrun{
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "FO",
#'               control = list(maxit = 5L, compute_inference = FALSE))
#' parts <- nm_ctl_template(4L, 4L)
#' nm_ctl_apply_fit_inits(parts, fit)
#' }
#' @export
nm_ctl_apply_fit_inits <- function(parts, fit) {
  if (is.null(parts) || is.null(fit) || is.null(fit$par) || is.null(fit$model)) {
    return(parts)
  }
  pp <- .nm_unpack(fit$model, fit$par)
  th <- parts$thetas
  if (!is.null(th) && nrow(th) > 0L) {
    n <- min(nrow(th), length(pp$theta))
    th$Value[seq_len(n)] <- as.numeric(pp$theta[seq_len(n)])
    parts$thetas <- th
  }
  om <- parts$omegas
  if (!is.null(om) && nrow(om) > 0L) {
    n <- min(nrow(om), length(pp$omega))
    om$Value[seq_len(n)] <- as.numeric(pp$omega[seq_len(n)])
    parts$omegas <- om
  }
  sg <- parts$sigmas
  if (!is.null(sg) && nrow(sg) > 0L) {
    n <- min(nrow(sg), length(pp$sigma))
    sg$Value[seq_len(n)] <- as.numeric(pp$sigma[seq_len(n)])
    parts$sigmas <- sg
  }
  parts
}

#' Build editable control-stream parts for a standard ADVAN/TRANS template
#'
#' Returns a list compatible with \code{\link{nm_ctl_compose}} with default
#' THETA/OMEGA/SIGMA, $PK, and $ERROR blocks for the requested model type.
#'
#' @param advan ADVAN number.
#' @param trans TRANS number.
#' @param data_file Relative path for \code{$DATA} (e.g. \code{"data/theo.csv"}).
#' @param problem Text for \code{$PROBLEM} (ADVAN/TRANS appended automatically).
#' @return List of control-stream parts (no \code{$EST}).
#' @examples
#' nm_ctl_template(4L, 4L, data_file = "data/theo.csv")
#' @export
nm_ctl_template <- function(advan,
                            trans = NULL,
                            data_file = "data.csv",
                            problem = "Template model",
                            ode_ncomp = 2L) {
  advan <- as.integer(advan)
  if (advan == 10L) {
    ode_ncomp <- 1L
  }
  if (is.null(trans)) {
    trans <- as.integer(nm_ctl_default_trans(advan))
  } else {
    trans <- as.integer(trans)
  }
  if (!nm_ctl_is_valid_pair(advan, trans)) {
    trans <- as.integer(nm_ctl_default_trans(advan))
  }
  trans <- nm_ctl_effective_trans(advan, trans)
  spec <- .nm_ctl_template_spec(advan, trans)
  spec$pk <- .nm_ctl_append_pk_scaling(spec$pk, advan, trans, ode_ncomp = ode_ncomp)
  if (nm_ctl_use_ode(advan)) {
    if (advan == 10L) {
      if (!nzchar(spec$des %||% "")) {
        spec$des <- "DADT(1) = -(VM*A(1)/(KM + A(1)))"
      }
    } else {
      spec$des <- nm_ctl_default_des(advan, trans, ncomp = ode_ncomp, oral = TRUE)
    }
  }
  list(
    problem = problem,
    advan = advan,
    trans = trans,
    use_ode = nm_ctl_use_ode(advan),
    ode_ncomp = ode_ncomp,
    subroutine = nm_ctl_subroutine_text(advan, trans),
    data_file = data_file,
    input_cols = nm_ctl_essential_input_cols(advan, trans),
    output_cols = character(),
    thetas = spec$thetas,
    omegas = spec$omegas,
    sigmas = spec$sigmas,
    pk = spec$pk,
    des = .nm_ctl_coalesce(spec$des, ""),
    error = spec$error
  )
}
