#' Define a correlated residual endpoint group
#'
#' Declares endpoint records observed at the same subject and time as one
#' multivariate Gaussian residual contribution. Marginal variances still come
#' from the model's residual-error definition and `SIGMAS`; this object defines
#' only the cross-endpoint correlation matrix.
#'
#' @param dvid Unique DVID values in matrix order.
#' @param correlation Square symmetric matrix. Diagonal values are set to one.
#'   Off-diagonal entries may be fixed numbers in `(-1, 1)` or parameter
#'   references written as `THETA(i)` or `SIGMA(i)`. One triangle may be left
#'   empty and is mirrored automatically.
#' @param parameter_transform Transform for parameter references. `"tanh"`
#'   maps unconstrained parameters into `(-1, 1)`; `"identity"` uses the
#'   parameter value directly.
#' @param label Optional human-readable group label.
#' @return A serializable residual-group declaration.
#' @export
nm_residual_group <- function(
    dvid, correlation,
    parameter_transform = c("tanh", "identity"),
    label = NULL) {
  dvid <- as.numeric(dvid)
  if (length(dvid) < 2L || any(!is.finite(dvid)) || anyDuplicated(dvid)) {
    .nm_stop("`dvid` must contain at least two unique finite endpoint codes.")
  }
  correlation <- as.matrix(correlation)
  if (!identical(dim(correlation), rep(length(dvid), 2L))) {
    .nm_stop("`correlation` must have one row and column per DVID.")
  }
  entries <- matrix(
    trimws(as.character(correlation)), nrow = length(dvid), ncol = length(dvid),
    dimnames = list(dvid, dvid)
  )
  entries[is.na(entries)] <- ""
  diag(entries) <- "1"
  for (row in seq_along(dvid)) {
    if (row == 1L) next
    for (column in seq_len(row - 1L)) {
      lower <- entries[row, column]
      upper <- entries[column, row]
      if (!nzchar(lower) && !nzchar(upper)) {
        .nm_stop("Every cross-endpoint correlation pair requires a value or parameter reference.")
      }
      if (nzchar(lower) && nzchar(upper) && !identical(lower, upper)) {
        .nm_stop("`correlation` must be symmetric; conflicting entries were supplied for DVID ",
                 dvid[[row]], " and ", dvid[[column]], ".")
      }
      value <- if (nzchar(lower)) lower else upper
      entries[row, column] <- entries[column, row] <- value
    }
  }
  source <- matrix("fixed", length(dvid), length(dvid))
  index <- matrix(0L, length(dvid), length(dvid))
  value <- matrix(0, length(dvid), length(dvid))
  diag(value) <- 1
  for (row in seq_along(dvid)) {
    for (column in seq_along(dvid)) {
      if (row == column) next
      entry <- toupper(gsub("[[:space:]_]", "", entries[row, column]))
      matched <- regmatches(
        entry, regexec("^(THETA|SIGMA)\\(([1-9][0-9]*)\\)$", entry)
      )[[1L]]
      if (length(matched)) {
        source[row, column] <- tolower(matched[[2L]])
        index[row, column] <- as.integer(matched[[3L]])
      } else {
        fixed <- suppressWarnings(as.numeric(entry))
        if (length(fixed) != 1L || !is.finite(fixed) || abs(fixed) >= 1) {
          .nm_stop(
            "Correlation entries must be fixed values in (-1, 1), THETA(i), or SIGMA(i)."
          )
        }
        value[row, column] <- fixed
      }
    }
  }
  parameter_transform <- match.arg(parameter_transform)
  if (is.null(label)) label <- paste0("DVID ", paste(dvid, collapse = "/"))
  if (length(label) != 1L || is.na(label) || !nzchar(trimws(label))) {
    .nm_stop("`label` must be one non-empty string.")
  }
  structure(
    list(
      version = 1L, label = as.character(label), dvid = dvid,
      correlation = entries, source = source, index = index, value = value,
      parameter_transform = parameter_transform
    ),
    class = "nm_residual_group"
  )
}

.nm_residual_groups <- function(value) {
  if (is.null(value)) return(NULL)
  if (inherits(value, "nm_residual_group")) value <- list(value)
  if (!is.list(value) || !length(value)) {
    .nm_stop("`residual_groups` must contain one or more `nm_residual_group()` objects.")
  }
  value <- lapply(value, function(group) {
    if (inherits(group, "nm_residual_group")) return(group)
    if (!is.list(group)) .nm_stop("Every residual group must be created by `nm_residual_group()`.")
    do.call(nm_residual_group, group[c(
      "dvid", "correlation", "parameter_transform", "label"
    )])
  })
  all_dvid <- unlist(lapply(value, `[[`, "dvid"), use.names = FALSE)
  if (anyDuplicated(all_dvid)) {
    .nm_stop("A DVID may occur in only one correlated residual group.")
  }
  value
}

.nm_residual_group_value <- function(group, theta, sigma) {
  result <- group$value
  for (row in seq_len(nrow(result))) {
    for (column in seq_len(ncol(result))) {
      source <- group$source[row, column]
      if (identical(source, "fixed")) next
      parameters <- if (identical(source, "theta")) theta else sigma
      index <- group$index[row, column]
      if (index < 1L || index > length(parameters)) {
        .nm_stop("Residual correlation parameter index is outside the ",
                 toupper(source), " table.")
      }
      parameter <- parameters[[index]]
      result[row, column] <- if (group$parameter_transform == "tanh") {
        tanh(parameter)
      } else parameter
    }
  }
  if (any(!is.finite(result)) || any(abs(result[row(result) != col(result)]) >= 1)) {
    .nm_stop("Resolved residual correlations must be finite and strictly between -1 and 1.")
  }
  if (isTRUE(tryCatch({ chol(result); TRUE }, error = function(e) FALSE))) return(result)
  .nm_stop("Resolved cross-endpoint residual correlation matrix is not positive definite.")
}

#' @export
print.nm_residual_group <- function(x, ...) {
  cat("LibeRation correlated residual group\n")
  cat("  label:", x$label, " DVID:", paste(x$dvid, collapse = ", "), "\n")
  cat("  parameter transform:", x$parameter_transform, "\n")
  invisible(x)
}
