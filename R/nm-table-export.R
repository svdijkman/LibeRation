#' Write a NONMEM-style $TABLE output file for a fit
#'
#' Exports a whitespace-delimited table (with a header line) of the standard
#' NONMEM diagnostic quantities for an \code{nm_fit}: \code{ID TIME DV PRED
#' IPRED RES WRES CWRES} plus the empirical-Bayes \code{ETA1..ETAn}. Additional
#' quantities produced by \code{\link{nm_add_cwres}} (\code{CPRED}, \code{CRES},
#' \code{IWRES}) are also available.
#'
#' @param fit An \code{nm_fit} object.
#' @param file Output path. If \code{NULL} the table is returned but not written.
#' @param columns Character vector of columns to export. Defaults to the NONMEM
#'   standard set (\code{ID TIME DV PRED IPRED RES WRES CWRES}) followed by the
#'   ETA columns. Unknown columns are dropped with a warning.
#' @param data Optional dataset (defaults to \code{fit$data}).
#' @param firstonly Logical; if \code{TRUE} write one row per subject (the first
#'   record), mimicking NONMEM's \code{FIRSTONLY}.
#' @param nonmem_header Logical; if \code{TRUE} prepend a \code{TABLE NO. 1}
#'   banner line (NONMEM convention). Default \code{FALSE}.
#' @param digits Number of significant digits for numeric formatting.
#' @return Invisibly, the exported \code{data.frame}.
#' @examples
#' \dontrun{
#' sim <- nm_synthetic_theo(n_sub = 3L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "FOCE",
#'               control = list(maxit = 10L, compute_inference = FALSE))
#' tab <- nm_write_table(fit, tempfile(fileext = ".tab"))
#' head(tab)
#' }
#' @export
nm_write_table <- function(fit,
                           file = NULL,
                           columns = NULL,
                           data = fit$data,
                           firstonly = FALSE,
                           nonmem_header = FALSE,
                           digits = 6L) {
  if (is.null(fit) || !inherits(fit, "nm_fit")) {
    .nm_stop("fit must be an nm_fit object.")
  }
  tab <- .nm_table_build(fit, data = data, columns = columns)
  if (isTRUE(firstonly)) {
    keep <- !duplicated(tab$ID)
    tab <- tab[keep, , drop = FALSE]
  }
  rownames(tab) <- NULL
  if (!is.null(file)) {
    .nm_table_write_file(tab, file, nonmem_header = nonmem_header, digits = digits)
  }
  invisible(tab)
}

#' Read a table written by \code{\link{nm_write_table}} back into a data.frame
#'
#' @param file Path to a table file (whitespace-delimited, optional NONMEM
#'   banner line).
#' @return A \code{data.frame}.
#' @examples
#' \dontrun{
#' nm_read_table("run1.tab")
#' }
#' @export
nm_read_table <- function(file) {
  if (!file.exists(file)) {
    .nm_stop("File not found: ", file)
  }
  lines <- readLines(file, warn = FALSE)
  skip <- 0L
  if (length(lines) > 0L && grepl("^TABLE NO", lines[[1L]])) {
    skip <- 1L
  }
  utils::read.table(
    file, header = TRUE, skip = skip, stringsAsFactors = FALSE,
    check.names = TRUE
  )
}

#' @keywords internal
.nm_table_default_columns <- function() {
  c("ID", "TIME", "DV", "PRED", "IPRED", "RES", "WRES", "CWRES")
}

#' @keywords internal
.nm_table_build <- function(fit, data = fit$data, columns = NULL) {
  if (is.null(fit$gof)) {
    fit <- nm_add_cwres(fit, data = data)
  }
  gof <- as.data.frame(fit$gof)
  n_eta <- .nm_n_eta(fit$model)
  eta_names <- character()
  if (n_eta > 0L) {
    eta_mat <- fit$cwres_eta %||% fit$eta
    if (is.matrix(eta_mat) && nrow(eta_mat) > 0L) {
      ids <- .nm_subject_ids(data.table::as.data.table(gof))
      idx <- match(gof$ID, ids)
      eta_names <- paste0("ETA", seq_len(ncol(eta_mat)))
      for (k in seq_len(ncol(eta_mat))) {
        gof[[eta_names[k]]] <- eta_mat[idx, k]
      }
    }
  }
  if (is.null(columns)) {
    columns <- c(.nm_table_default_columns(), eta_names)
  }
  missing_cols <- setdiff(columns, names(gof))
  if (length(missing_cols) > 0L) {
    warning(
      "Dropping unavailable table column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
    columns <- intersect(columns, names(gof))
  }
  if (length(columns) == 0L) {
    .nm_stop("No valid columns available for the table export.")
  }
  gof[, columns, drop = FALSE]
}

#' @keywords internal
.nm_table_write_file <- function(tab, file, nonmem_header = FALSE, digits = 6L) {
  fmt <- as.data.frame(lapply(tab, function(col) {
    if (is.numeric(col)) {
      formatC(col, digits = digits, format = "g")
    } else {
      as.character(col)
    }
  }), stringsAsFactors = FALSE)
  names(fmt) <- names(tab)
  con <- file(file, open = "wt")
  on.exit(close(con), add = TRUE)
  if (isTRUE(nonmem_header)) {
    writeLines("TABLE NO.  1", con)
  }
  writeLines(paste(names(fmt), collapse = " "), con)
  utils::write.table(
    fmt, con, row.names = FALSE, col.names = FALSE,
    quote = FALSE, sep = " "
  )
  invisible(file)
}
