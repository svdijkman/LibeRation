#' Define a random-effect grouping block
#'
#' @param name Stable block label.
#' @param column Dataset column identifying the unit receiving the effect.
#' @param etas ETA indices belonging to the block. Correlation is permitted
#'   within a block; different blocks are independent.
#' @return A serializable random-effect block declaration.
#' @export
nm_re_block <- function(name, column, etas) {
  name <- trimws(as.character(name))
  column <- trimws(as.character(column))
  etas <- sort(unique(as.integer(etas)))
  if (length(name) != 1L || is.na(name) || !nzchar(name) ||
      length(column) != 1L || is.na(column) || !nzchar(column) ||
      !length(etas) || anyNA(etas) || any(etas < 1L)) {
    .nm_stop("A random-effect block requires a name, grouping column, and positive ETA indices.")
  }
  structure(list(name = name, column = column, etas = etas), class = "nm_re_block")
}

#' Configure nested or crossed random effects
#'
#' Blocks may be nested or crossed. When `cluster` is omitted, LibeRation
#' constructs connected components from subject and block memberships; each
#' component becomes one conditional-likelihood unit. This is exact but a fully
#' crossed design may form one large conditional mode. Supplying `cluster`
#' avoids graph discovery, provided no grouping unit spans clusters.
#'
#' @param ... [nm_re_block()] declarations or one list of them.
#' @param cluster Optional dataset column defining independent top-level
#'   clusters.
#' @return A serializable random-effect design.
#' @export
nm_re_config <- function(..., cluster = NULL) {
  blocks <- list(...)
  if (length(blocks) == 1L && is.list(blocks[[1L]]) &&
      !inherits(blocks[[1L]], "nm_re_block")) blocks <- blocks[[1L]]
  if (!length(blocks) || any(!vapply(blocks, inherits, logical(1), "nm_re_block"))) {
    .nm_stop("`nm_re_config()` requires one or more `nm_re_block()` declarations.")
  }
  names <- vapply(blocks, `[[`, character(1), "name")
  etas <- unlist(lapply(blocks, `[[`, "etas"), use.names = FALSE)
  if (anyDuplicated(names)) .nm_stop("Random-effect block names must be unique.")
  if (anyDuplicated(etas)) .nm_stop("Each ETA may belong to only one random-effect block.")
  if (!is.null(cluster)) {
    cluster <- trimws(as.character(cluster))
    if (length(cluster) != 1L || is.na(cluster) || !nzchar(cluster)) {
      .nm_stop("`cluster` must be NULL or one dataset-column name.")
    }
  }
  structure(
    list(version = 1L, blocks = unname(blocks), cluster = cluster),
    class = "nm_re_config"
  )
}

.nm_re_config <- function(config, n_eta = NULL) {
  if (is.null(config)) return(NULL)
  if (!inherits(config, "nm_re_config")) {
    if (!is.list(config)) .nm_stop("RE_CONFIG must be created by `nm_re_config()`.")
    config$version <- NULL
    config <- do.call(nm_re_config, config)
  }
  if (!is.null(n_eta)) {
    assigned <- sort(unlist(lapply(config$blocks, `[[`, "etas"), use.names = FALSE))
    if (!identical(assigned, seq_len(as.integer(n_eta)))) {
      .nm_stop("RE_CONFIG must assign every ETA exactly once; expected ETA(1)-ETA(",
               as.integer(n_eta), ").")
    }
  }
  config
}

.nm_re_union_components <- function(data, columns) {
  rows <- nrow(data)
  parent <- seq_len(rows)
  find <- function(index) {
    while (parent[[index]] != index) {
      parent[[index]] <<- parent[[parent[[index]]]]
      index <- parent[[index]]
    }
    index
  }
  unite <- function(first, second) {
    first <- find(first); second <- find(second)
    if (first != second) parent[[second]] <<- first
  }
  for (column in columns) {
    values <- as.character(data[[column]])
    first <- new.env(parent = emptyenv(), hash = TRUE)
    for (row in seq_len(rows)) {
      key <- values[[row]]
      if (exists(key, first, inherits = FALSE)) unite(row, get(key, first))
      else assign(key, row, first)
    }
  }
  roots <- vapply(seq_len(rows), find, integer(1))
  match(roots, unique(roots))
}

.nm_re_engine_data <- function(model, data) {
  config <- model$RE_CONFIG
  if (is.null(config)) return(data)
  columns <- unique(c(
    vapply(config$blocks, `[[`, character(1), "column"), config$cluster
  ))
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    .nm_stop("Random-effect design column(s) missing from the dataset: ",
             paste(missing, collapse = ", "), ".")
  }
  for (column in columns) {
    value <- as.character(data[[column]])
    if (anyNA(value) || any(!nzchar(value))) {
      .nm_stop("Random-effect grouping column `", column, "` contains missing/empty identifiers.")
    }
  }
  data$.STRUCT_ID_INDEX <- data$.ID_INDEX
  if (is.null(config$cluster)) {
    data$.ID_INDEX <- .nm_re_union_components(
      data, unique(c(".STRUCT_ID_INDEX", vapply(config$blocks, `[[`, character(1), "column")))
    )
  } else {
    clusters <- as.character(data[[config$cluster]])
    data$.ID_INDEX <- match(clusters, unique(clusters))
    for (block in config$blocks) {
      membership <- split(data$.ID_INDEX, as.character(data[[block$column]]))
      if (any(vapply(membership, function(value) length(unique(value)) > 1L, logical(1)))) {
        .nm_stop("Random-effect units in `", block$column,
                 "` span `", config$cluster, "` clusters.")
      }
    }
  }
  offset <- 0L
  for (block_index in seq_along(config$blocks)) {
    block <- config$blocks[[block_index]]
    local <- integer(nrow(data))
    totals <- integer(max(data$.ID_INDEX))
    for (cluster in seq_len(max(data$.ID_INDEX))) {
      rows <- which(data$.ID_INDEX == cluster)
      values <- as.character(data[[block$column]][rows])
      local[rows] <- match(values, unique(values))
      totals[[cluster]] <- length(unique(values))
    }
    total_column <- paste0(".RE_TOTAL_", block_index)
    stored_total <- if (total_column %in% names(data)) {
      suppressWarnings(max(as.integer(data[[total_column]]), na.rm = TRUE))
    } else 0L
    if (!is.finite(stored_total)) stored_total <- 0L
    maximum <- max(max(totals), stored_total)
    data[[paste0(".RE_UNIT_", block_index)]] <- local
    data[[total_column]] <- maximum
    for (within in seq_along(block$etas)) {
      eta <- block$etas[[within]]
      data[[paste0(".ETA_COLUMN_", eta)]] <-
        offset + (local - 1L) * length(block$etas) + within
    }
    offset <- offset + maximum * length(block$etas)
  }
  order_index <- order(
    data$.ID_INDEX, data$.STRUCT_ID_INDEX, data$TIME,
    data$.sort_priority, data$.source_row, method = "radix"
  )
  data <- data[order_index, , drop = FALSE]
  rownames(data) <- NULL
  attr(data, "re_eta_columns") <- offset
  attr(data, "re_config") <- config
  class(data) <- unique(c("nm_dataset", class(data)))
  data
}

#' @export
print.nm_re_config <- function(x, ...) {
  cat("LibeRation random-effect design\n")
  cat("  independent cluster:", x$cluster %||% "auto-connected components", "\n")
  for (block in x$blocks) {
    cat("  ", block$name, ": ", block$column, " -> ETA(",
        paste(block$etas, collapse = ","), ")\n", sep = "")
  }
  invisible(x)
}
