.nm_event_defaults <- list(
  EVID = 0L, AMT = 0, RATE = 0, CMT = 0L, SS = 0L, II = 0,
  ADDL = 0L, MDV = 0L, DV = NA_real_
)

.nm_expand_addl <- function(data) {
  generated <- vector("list", 0L)
  for (i in seq_len(nrow(data))) {
    n <- as.integer(data$ADDL[[i]])
    if (is.na(n) || n < 0L) .nm_stop("ADDL must be a non-negative integer at row ", i, ".")
    if (n == 0L) next
    if (!is.finite(data$II[[i]]) || data$II[[i]] <= 0) {
      .nm_stop("ADDL > 0 requires II > 0 at row ", i, ".")
    }
    if (!(data$EVID[[i]] %in% c(1L, 4L)) || data$AMT[[i]] <= 0) {
      .nm_stop("ADDL is only valid on positive dosing records (row ", i, ").")
    }
    for (k in seq_len(n)) {
      row <- data[i, , drop = FALSE]
      row$TIME <- row$TIME + k * row$II
      row$ADDL <- 0L
      row$.generated <- TRUE
      row$.source_row <- data$.source_row[[i]]
      row$.sort_priority <- -1L
      generated[[length(generated) + 1L]] <- row
    }
  }
  data$ADDL <- 0L
  if (length(generated)) {
    data <- rbind(data, do.call(rbind, generated))
  }
  data
}

#' Normalize a NONMEM event dataset
#'
#' @param data Data frame containing at least `ID` and `TIME`.
#' @param expand_addl Materialize ADDL/II doses before C++ execution.
#' @return An `nm_dataset` data frame with stable event ordering.
#' @examples
#' events <- data.frame(
#'   ID = 1, TIME = c(0, 1, 2), EVID = c(1, 0, 0),
#'   AMT = c(100, 0, 0), DV = c(NA, 4.2, 3.1), MDV = c(1, 0, 0)
#' )
#' nm_dataset(events)
#' @export
nm_dataset <- function(data, expand_addl = TRUE) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  if (!all(c("ID", "TIME") %in% names(data))) {
    .nm_stop("Dataset requires `ID` and `TIME` columns.")
  }
  if (!nrow(data)) .nm_stop("Dataset must contain at least one row.")
  data$.source_row <- seq_len(nrow(data))
  data$.generated <- FALSE
  data$.sort_priority <- 0L
  for (name in names(.nm_event_defaults)) {
    if (!name %in% names(data)) data[[name]] <- .nm_event_defaults[[name]]
  }
  numeric_cols <- c("TIME", "AMT", "RATE", "II", "DV")
  integer_cols <- c("EVID", "CMT", "SS", "ADDL", "MDV")
  for (name in numeric_cols) data[[name]] <- as.numeric(data[[name]])
  for (name in integer_cols) data[[name]] <- as.integer(data[[name]])
  for (name in setdiff(names(.nm_event_defaults), c("TIME", "DV"))) {
    data[[name]][is.na(data[[name]])] <- .nm_event_defaults[[name]]
  }
  if (anyNA(data$ID) || any(!nzchar(as.character(data$ID)))) {
    .nm_stop("ID must be non-missing and non-empty.")
  }
  if (any(!is.finite(data$TIME))) .nm_stop("TIME must be finite.")
  if (any(!is.finite(data$AMT)) || any(data$AMT < 0)) {
    .nm_stop("AMT must be finite and non-negative.")
  }
  if (any(!is.finite(data$RATE)) || any(data$RATE < 0 & !data$RATE %in% c(-1, -2))) {
    .nm_stop("RATE must be non-negative, -1 (modelled Rn), or -2 (modelled Dn).")
  }
  if (any(!is.finite(data$II)) || any(data$II < 0)) .nm_stop("II must be finite and non-negative.")
  if (anyNA(data$EVID) || any(!data$EVID %in% 0:4)) .nm_stop("EVID must be one of 0, 1, 2, 3, or 4.")
  if (anyNA(data$CMT) || any(data$CMT < 0L)) .nm_stop("CMT must be a non-negative integer.")
  if (anyNA(data$SS) || any(!data$SS %in% 0:2)) .nm_stop("SS must be 0, 1, or 2.")
  if (anyNA(data$MDV) || any(!data$MDV %in% 0:1)) .nm_stop("MDV must be 0 or 1.")
  if (isTRUE(expand_addl)) data <- .nm_expand_addl(data)
  id_levels <- unique(as.character(data$ID))
  data$.ID_INDEX <- match(as.character(data$ID), id_levels)
  order_index <- order(data$.ID_INDEX, data$TIME, data$.sort_priority,
                       data$.source_row, method = "radix")
  data <- data[order_index, , drop = FALSE]
  rownames(data) <- NULL
  for (id in unique(data$.ID_INDEX)) {
    times <- data$TIME[data$.ID_INDEX == id]
    if (is.unsorted(times, strictly = FALSE)) .nm_stop("TIME ordering failed for subject index ", id, ".")
  }
  attr(data, "id_levels") <- id_levels
  class(data) <- c("nm_dataset", "data.frame")
  data
}

#' @export
print.nm_dataset <- function(x, ...) {
  cat("LibeRation event dataset\n")
  cat("  rows:", nrow(x), " subjects:", length(unique(x$.ID_INDEX)),
      " generated events:", sum(x$.generated), "\n")
  print(utils::head(as.data.frame(x), ...))
  invisible(x)
}
