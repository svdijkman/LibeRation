# LibeRation Shiny GUI — Pirana-style workspace + ribbon

if (file.exists("../../DESCRIPTION") &&
    requireNamespace("pkgload", quietly = TRUE)) {
  desc <- read.dcf("../../DESCRIPTION")
  if (desc[1L, "Package"] == "LibeRation") {
    if (dir.exists("../../R")) {
      ad <- Sys.getenv("LibeRation_LIBERTAD_ROOT", "")
      if (!nzchar(ad) && dir.exists("../../../LibeRtAD")) {
        ad <- normalizePath("../../../LibeRtAD")
      }
      if (nzchar(ad)) {
        pkgload::load_all(ad, quiet = TRUE, compile = FALSE)
      }
      pkgload::load_all("../..", quiet = TRUE, recompile = FALSE)
    }
  }
}

if (!requireNamespace("LibeRation", quietly = TRUE)) {
  stop("Package LibeRation must be loaded or installed.")
}
if (!requireNamespace("shiny", quietly = TRUE)) {
  stop("Package 'shiny' is required for the LibeRation GUI.")
}
if (!requireNamespace("DT", quietly = TRUE)) {
  stop("Package 'DT' is required for the LibeRation GUI.")
}

library(shiny, warn.conflicts = FALSE)

.valid_id <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
}

if (is.null(getOption("LibeRation.workspace"))) {
  nm_workspace_init(create_demo_project = TRUE)
}

.shiny_gof_obs <- function(pred) {
  pred[pred$MDV == 0L & pred$EVID == 0L, ]
}

.shiny_profile_lines <- function(pred, yvar = "IPRED") {
  lines <- pred[pred$MDV == 0L & pred$EVID == 0L, ]
  if (nrow(lines) == 0L) {
    return(lines)
  }
  lines <- lines[order(lines$ID, lines$TIME), ]
  lines[[yvar]] <- as.numeric(lines[[yvar]])
  lines
}

.shiny_col_dv <- "#222222"
.shiny_col_ipred <- "#D95F02"
.shiny_col_pred <- "#1B9E77"

.shiny_has_ggplot <- function() {
  requireNamespace("ggplot2", quietly = TRUE)
}

.shiny_named_trans_choices <- function(advan) {
  if (!nm_ctl_show_trans(as.integer(advan))) {
    return(NULL)
  }
  choices <- nm_ctl_trans_choices(advan)
  if (length(choices) == 0L) {
    return(NULL)
  }
  stats::setNames(choices, paste0("TRANS", choices))
}

.shiny_numeric_advan_choices <- function() {
  ch <- nm_ctl_advan_choices()
  stats::setNames(ch, ch)
}

.shiny_numeric_trans_choices <- function(advan) {
  if (!nm_ctl_show_trans(as.integer(advan))) {
    return(NULL)
  }
  ch <- nm_ctl_trans_choices(advan)
  if (length(ch) == 0L) {
    return(NULL)
  }
  stats::setNames(ch, ch)
}

.shiny_ctl_param_labels <- function(df, prefix) {
  if (is.null(df) || nrow(df) == 0L) {
    return(df)
  }
  idx_col <- prefix
  if (!idx_col %in% names(df)) {
    return(df)
  }
  if (!"Label" %in% names(df)) {
    df$Label <- paste0(prefix, df[[idx_col]])
  }
  for (i in seq_len(nrow(df))) {
    lbl <- trimws(as.character(df$Label[[i]]))
    if (!nzchar(lbl)) {
      df$Label[[i]] <- paste0(prefix, df[[idx_col]][[i]])
    }
  }
  df
}

.shiny_columns_picker_ui <- function(df) {
  if (nrow(df) == 0L) {
    return(tags$p(
      class = "text-muted",
      style = "font-size: 12px;",
      "Select a dataset to choose $INPUT / $OUTPUT columns."
    ))
  }
  rows <- lapply(seq_len(nrow(df)), function(i) {
    col <- df$Column[[i]]
    is_req <- identical(df$Required[[i]], "req")
    in_ds <- isTRUE(df$in_dataset[[i]])
    checked_in <- isTRUE(df$IN[[i]])
    checked_out <- isTRUE(df$OUT[[i]])
    tags$tr(
      tags$td(
        col,
        if (is_req) tags$span(class = "req-tag", "(required)")
      ),
      tags$td(
        if (in_ds) {
          tags$input(
            type = "checkbox",
            class = "ctl-col-pick",
            `data-col` = col,
            `data-kind` = "in",
            checked = if (checked_in) NA else NULL,
            disabled = if (is_req) NA else NULL
          )
        } else {
          tags$span("—", class = "text-muted")
        }
      ),
      tags$td(
        tags$input(
          type = "checkbox",
          class = "ctl-col-pick",
          `data-col` = col,
          `data-kind` = "out",
          checked = if (checked_out) NA else NULL
        )
      )
    )
  })
  tags$table(
    class = "ctl-col-picker-table",
    tags$thead(
      tags$tr(
        tags$th("Column"),
        tags$th("$INPUT"),
        tags$th("$OUTPUT")
      )
    ),
    tags$tbody(rows)
  )
}

.shiny_default_des <- function(advan = 6L, trans = 1L, ncomp = 2L, oral = TRUE) {
  fn <- get0("nm_ctl_default_des", envir = asNamespace("LibeRation"), inherits = FALSE)
  if (is.function(fn)) {
    return(fn(advan = advan, trans = trans, ncomp = ncomp, oral = oral))
  }
  advan <- as.integer(advan)
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
    lines <- c(lines, "F = A(2)/V")
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

.shiny_bin_midpoints <- function(levels_chr) {
  vapply(levels_chr, function(lbl) {
    nums <- suppressWarnings(as.numeric(unlist(
      regmatches(lbl, gregexpr("[+-]?[0-9]*\\.?[0-9]+([eE][+-]?[0-9]+)?", lbl, perl = TRUE))
    )))
    if (length(nums) >= 2L) {
      mean(nums[1:2])
    } else if (length(nums) == 1L) {
      nums
    } else {
      NA_real_
    }
  }, numeric(1))
}

.shiny_parse_breaks <- function(text) {
  if (is.null(text) || length(text) == 0L) {
    return(NULL)
  }
  txt <- trimws(as.character(text)[1L])
  if (!nzchar(txt)) {
    return(NULL)
  }
  txt <- gsub("[;\n\r]+", ",", txt)
  parts <- strsplit(txt, ",", fixed = TRUE)[[1L]]
  vals <- suppressWarnings(as.numeric(trimws(parts)))
  vals <- vals[is.finite(vals)]
  vals <- sort(unique(vals))
  if (length(vals) < 2L) {
    NULL
  } else {
    vals
  }
}

.shiny_bin_var <- function(x, nbins, breaks = NULL) {
  nbins <- max(3L, min(50L, as.integer(nbins)))
  if (!is.numeric(x)) {
    return(factor(as.character(x), levels = unique(as.character(x)), ordered = TRUE))
  }
  ok <- !is.na(x)
  if (sum(ok) < 2L) {
    return(factor(as.character(x), ordered = TRUE))
  }
  br <- if (!is.null(breaks) && length(breaks) >= 2L) {
    sort(unique(as.numeric(breaks)))
  } else {
    br <- unique(stats::quantile(x[ok], probs = seq(0, 1, length.out = nbins + 1L), na.rm = TRUE))
    if (length(br) < 3L) {
      pretty(range(x, na.rm = TRUE), n = nbins)
    } else {
      br
    }
  }
  b <- cut(x, breaks = br, include.lowest = TRUE, dig.lab = 4)
  lv <- levels(b)
  mids <- .shiny_bin_midpoints(lv)
  ord <- order(mids, lv)
  factor(b, levels = lv[ord], ordered = TRUE)
}

.shiny_explore_x_levels <- function(x) {
  if (is.factor(x)) {
    levels(x)
  } else {
    unique(as.character(x))
  }
}

.shiny_explore_x_numeric <- function(x) {
  if (is.factor(x)) {
    lv <- levels(x)
    mids <- setNames(.shiny_bin_midpoints(lv), lv)
    mids[as.character(x)]
  } else {
    as.numeric(x)
  }
}

.shiny_explore_bin_mean <- function(d) {
  lv <- .shiny_explore_x_levels(d$x)
  mids <- setNames(.shiny_bin_midpoints(lv), lv)
  has_panel <- "panel" %in% names(d) && length(unique(d$panel)) > 1L
  panels <- if (has_panel) levels(d$panel) else "All"
  groups <- unique(d$group)
  rows <- vector("list", length(lv) * length(groups) * length(panels))
  k <- 0L
  for (pn in panels) {
    dp <- if (has_panel) d[d$panel == pn, , drop = FALSE] else d
    for (g in groups) {
      dg <- dp[dp$group == g, , drop = FALSE]
      for (xl in lv) {
        z <- dg$y[as.character(dg$x) == xl]
        z <- z[!is.na(z)]
        k <- k + 1L
        rows[[k]] <- data.frame(
          x = xl,
          group = g,
          panel = pn,
          y = if (length(z)) mean(z) else NA_real_,
          n = length(z),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  out <- do.call(rbind, rows)
  out$x <- factor(out$x, levels = lv, ordered = TRUE)
  out$x_num <- mids[as.character(out$x)]
  if (has_panel) {
    out$panel <- factor(out$panel, levels = panels)
  }
  out[!is.na(out$y), , drop = FALSE]
}

.shiny_explore_bin_quantiles <- function(d, q_interval = 95, by_group = TRUE) {
  q_interval <- max(50, min(99.9, as.numeric(q_interval)))
  q_lo <- (100 - q_interval) / 200
  q_hi <- 1 - q_lo
  lv <- .shiny_explore_x_levels(d$x)
  mids <- setNames(.shiny_bin_midpoints(lv), lv)
  has_panel <- "panel" %in% names(d) && length(unique(d$panel)) > 1L
  panels <- if (has_panel) levels(d$panel) else "All"
  rows <- list()
  k <- 0L
  for (pn in panels) {
    dp <- if (has_panel) d[d$panel == pn, , drop = FALSE] else d
    groups <- if (isTRUE(by_group)) unique(dp$group) else "All"
    for (g in groups) {
      dg <- if (isTRUE(by_group)) dp[dp$group == g, , drop = FALSE] else dp
      for (xl in lv) {
        z <- dg$y[as.character(dg$x) == xl]
        z <- z[!is.na(z)]
        k <- k + 1L
        if (length(z) == 0L) {
          rows[[k]] <- data.frame(
            x = xl, group = g, panel = pn,
            y = NA_real_, y_lo = NA_real_, y_hi = NA_real_, n = 0L,
            stringsAsFactors = FALSE
          )
        } else {
          rows[[k]] <- data.frame(
            x = xl,
            group = g,
            panel = pn,
            y = stats::median(z),
            y_lo = unname(stats::quantile(z, q_lo)),
            y_hi = unname(stats::quantile(z, q_hi)),
            n = length(z),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  out <- do.call(rbind, rows)
  out$x <- factor(out$x, levels = lv, ordered = TRUE)
  out$x_num <- mids[as.character(out$x)]
  if (has_panel) {
    out$panel <- factor(out$panel, levels = panels)
  }
  out[!is.na(out$y), , drop = FALSE]
}

.shiny_explore_add_points <- function(p, d, use_strat, dodge,
                                       alpha_pt = 0.65, size = 0.85,
                                       scatter = 0.25, dodge_width = 0.75,
                                       dodge_points = use_strat, shape = 16L) {
  scatter <- max(0, min(1, as.numeric(scatter %||% 0)))
  jitter_width <- scatter * 0.45
  pt_size <- .shiny_explore_point_size(size)
  pt_shape <- as.integer(shape %||% 16L)
  if (scatter > 0) {
    if (use_strat && dodge_points) {
      p + ggplot2::geom_point(
        data = d,
        ggplot2::aes(x = .data[["x_plot"]], y = .data[["y_plot"]], colour = .data[["group"]]),
        alpha = alpha_pt, size = pt_size, shape = pt_shape,
        position = ggplot2::position_jitterdodge(
          jitter.width = jitter_width,
          jitter.height = 0,
          dodge.width = dodge_width
        ),
        inherit.aes = FALSE
      )
    } else if (use_strat) {
      p + ggplot2::geom_jitter(
        data = d,
        ggplot2::aes(x = .data[["x_plot"]], y = .data[["y_plot"]], colour = .data[["group"]]),
        alpha = alpha_pt, size = pt_size, shape = pt_shape, width = jitter_width, height = 0,
        inherit.aes = FALSE
      )
    } else {
      p + ggplot2::geom_jitter(
        data = d,
        ggplot2::aes(x = .data[["x_plot"]], y = .data[["y_plot"]]),
        alpha = alpha_pt, size = pt_size, shape = pt_shape, width = jitter_width, height = 0,
        colour = "black", inherit.aes = FALSE
      )
    }
  } else if (use_strat && dodge_points) {
    p + ggplot2::geom_point(
      data = d,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]], colour = .data[["group"]]),
      alpha = alpha_pt, size = pt_size, shape = pt_shape, position = dodge, inherit.aes = FALSE
    )
  } else if (use_strat) {
    p + ggplot2::geom_point(
      data = d,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]], colour = .data[["group"]]),
      alpha = alpha_pt, size = pt_size, shape = pt_shape, inherit.aes = FALSE
    )
  } else {
    p + ggplot2::geom_point(
      data = d,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
      alpha = alpha_pt, size = pt_size, shape = pt_shape, colour = "black", inherit.aes = FALSE
    )
  }
}

.shiny_explore_sort <- function(d, sort_cols = NULL) {
  if (is.factor(d$x)) {
    d$x <- factor(d$x, levels = levels(d$x), ordered = TRUE)
  }
  if (is.null(sort_cols)) {
    sort_cols <- if ("panel" %in% names(d)) {
      c("panel", "x", "group")
    } else {
      c("x", "group")
    }
  }
  sort_cols <- intersect(sort_cols, names(d))
  if (length(sort_cols) == 0L) {
    return(d)
  }
  ord <- do.call(order, d[sort_cols])
  d <- d[ord, , drop = FALSE]
  rownames(d) <- NULL
  d
}

.shiny_explore_line_id_col <- function(df, strat = NULL) {
  if (nzchar(strat) && strat %in% names(df)) {
    return(strat)
  }
  if ("ID" %in% names(df)) {
    return("ID")
  }
  NULL
}

.shiny_explore_agg_line <- function(d, how = c("mean", "median"), by_group = FALSE) {
  how <- match.arg(how)
  agg_fun <- if (how == "mean") {
    function(z) mean(z, na.rm = TRUE)
  } else {
    function(z) stats::median(z, na.rm = TRUE)
  }
  keys <- c("x", "x_plot")
  if (isTRUE(by_group) && "group" %in% names(d) && length(unique(d$group)) > 1L) {
    keys <- c(keys, "group")
  }
  if ("panel" %in% names(d) && nlevels(d$panel) > 1L) {
    keys <- c(keys, "panel")
  }
  keys <- unique(keys[keys %in% names(d)])
  rows <- lapply(split(d, interaction(d[keys], drop = TRUE)), function(chunk) {
    if (nrow(chunk) == 0L) {
      return(NULL)
    }
    y_val <- agg_fun(chunk$y_plot)
    out <- chunk[1L, keys, drop = FALSE]
    out$y_plot <- y_val
    out$y <- y_val
    out
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) {
    return(d[0L, , drop = FALSE])
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  sort_keys <- c("panel", "group", "x_plot")
  sort_keys <- intersect(sort_keys, names(out))
  out[do.call(order, out[sort_keys]), , drop = FALSE]
}

.shiny_explore_point_size <- function(point_size) {
  as.numeric(point_size %||% 0.85) * 2.2
}

.shiny_apply_facet <- function(p, use_split) {
  if (is.null(p) || !isTRUE(use_split) || !.shiny_has_ggplot()) {
    return(p)
  }
  p + ggplot2::facet_wrap(~ panel, scales = "free")
}

.shiny_data_explore_plot <- function(df, xvar, yvar, strat = NULL, split_by = NULL,
                                    bin_x = FALSE, bin_y = FALSE, nbins = 10L,
                                    bin_x_breaks = NULL, bin_y_breaks = NULL,
                                    bin_x_pos = "equal", bin_y_pos = "equal",
                                    plot_type = "points", show_points = FALSE,
                                    q_interval = 95, shade_alpha = 0.25,
                                    plot_title = NULL, xlab_custom = NULL, ylab_custom = NULL,
                                    point_scatter = 0.25, point_size = 0.85,
                                    point_shape = 16L,
                                    line_mode = "individual") {
  if (is.null(df) || nrow(df) == 0L || !nzchar(xvar) || !nzchar(yvar)) {
    return(NULL)
  }
  if (!(xvar %in% names(df)) || !(yvar %in% names(df))) {
    return(NULL)
  }
  plot_type <- match.arg(plot_type, c(
    "points", "lines", "both", "smooth", "regression",
    "boxplot", "violin", "mean_se", "median_q", "jitter"
  ))
  line_mode <- match.arg(line_mode, c("individual", "mean", "median", "none"))
  group_plots <- plot_type %in% c("boxplot", "violin", "mean_se", "median_q")
  line_id_col <- if (plot_type %in% c("lines", "both")) {
    .shiny_explore_line_id_col(df, strat)
  } else {
    NULL
  }
  cols <- unique(c(
    xvar, yvar,
    if (nzchar(strat)) strat else character(),
    if (nzchar(split_by)) split_by else character(),
    if (!is.null(line_id_col)) line_id_col else character()
  ))
  d <- as.data.frame(df[, cols, drop = FALSE])
  x_raw <- d[[xvar]]
  y_raw <- d[[yvar]]
  force_bin_x <- isTRUE(bin_x) || (group_plots && is.numeric(x_raw))
  x_binned <- FALSE
  y_binned <- FALSE
  bin_x_pos <- match.arg(bin_x_pos, c("equal", "midpoint"))
  bin_y_pos <- match.arg(bin_y_pos, c("equal", "midpoint"))
  if (force_bin_x && is.numeric(x_raw)) {
    d$x <- .shiny_bin_var(x_raw, nbins, breaks = bin_x_breaks)
    x_binned <- TRUE
  } else {
    d$x <- x_raw
  }
  if (isTRUE(bin_y) && is.numeric(y_raw)) {
    d$y <- .shiny_bin_var(y_raw, nbins, breaks = bin_y_breaks)
    y_binned <- TRUE
  } else {
    d$y <- y_raw
  }
  if (x_binned && bin_x_pos == "midpoint") {
    d$x_plot <- .shiny_explore_x_numeric(d$x)
  } else {
    d$x_plot <- d$x
  }
  if (y_binned && bin_y_pos == "midpoint") {
    d$y_plot <- .shiny_explore_x_numeric(d$y)
  } else {
    d$y_plot <- d$y
  }
  if (nzchar(strat) && strat %in% names(d)) {
    d$group <- as.character(d[[strat]])
  } else {
    d$group <- "All"
  }
  if (nzchar(split_by) && split_by %in% names(d)) {
    d$panel <- factor(d[[split_by]])
  } else {
    d$panel <- factor("All", levels = "All")
  }
  if (!is.null(line_id_col) && line_id_col %in% names(d)) {
    d$line_id <- as.character(d[[line_id_col]])
  } else {
    d$line_id <- "All"
  }
  d <- d[!is.na(d$x) & !is.na(d$y), , drop = FALSE]
  if (nrow(d) == 0L) {
    return(NULL)
  }
  title <- if (!is.null(plot_title) && nzchar(trimws(plot_title))) {
    trimws(plot_title)
  } else {
    paste(yvar, "vs", xvar)
  }
  xlab <- if (!is.null(xlab_custom) && nzchar(trimws(xlab_custom))) {
    trimws(xlab_custom)
  } else if (x_binned) {
    paste(xvar, "(binned)")
  } else {
    xvar
  }
  ylab <- if (!is.null(ylab_custom) && nzchar(trimws(ylab_custom))) {
    trimws(ylab_custom)
  } else if (isTRUE(bin_y) && is.numeric(y_raw)) {
    paste(yvar, "(binned)")
  } else {
    yvar
  }
  use_strat <- nzchar(strat) && length(unique(d$group)) > 1L
  use_split <- nzchar(split_by) && split_by %in% names(d) && nlevels(d$panel) > 1L
  x_discrete <- is.factor(d$x) || is.character(d$x)
  pt_size <- .shiny_explore_point_size(point_size)
  pt_shape <- as.integer(point_shape %||% 16L)
  xp <- if (x_binned && bin_x_pos == "midpoint") "x_plot" else "x"
  yp <- if (y_binned && bin_y_pos == "midpoint") "y_plot" else "y"
  use_dodge <- use_strat && (x_discrete || x_binned) &&
    !plot_type %in% c("lines", "both")
  colour_col <- if (use_strat) {
    "group"
  } else if (
    plot_type %in% c("lines", "both") &&
    line_mode == "individual" &&
    length(unique(d$line_id)) > 1L
  ) {
    "line_id"
  } else {
    NULL
  }
  use_colour <- !is.null(colour_col)
  colour_lbl <- if (use_strat) {
    strat
  } else if (identical(colour_col, "line_id") && !is.null(line_id_col)) {
    line_id_col
  } else {
    NULL
  }
  if (is.character(d$x)) {
    d$x <- factor(d$x, levels = unique(d$x), ordered = TRUE)
  }
  sum_d <- if (x_binned) .shiny_explore_bin_mean(d) else NULL
  q_d <- if (plot_type == "median_q") {
    .shiny_explore_bin_quantiles(d, q_interval = q_interval, by_group = FALSE)
  } else {
    NULL
  }

  if (.shiny_has_ggplot()) {
    p <- ggplot2::ggplot(d, ggplot2::aes(x = .data[[xp]], y = .data[[yp]])) +
      ggplot2::labs(title = title, x = xlab, y = ylab) +
      ggplot2::theme_bw(base_size = 11)
    if (use_strat && plot_type %in% c("boxplot", "violin")) {
      p <- p + ggplot2::aes(colour = .data[["group"]], fill = .data[["group"]])
    } else if (use_colour) {
      p <- p + ggplot2::aes(colour = .data[[colour_col]])
    }
    dodge <- ggplot2::position_dodge(width = 0.75)
    alpha_pt <- 0.65

    switch(
      plot_type,
      points = {
        if (use_colour && use_dodge) {
          p <- p + ggplot2::geom_point(
            alpha = alpha_pt, size = pt_size, shape = pt_shape, position = dodge
          )
        } else if (use_colour) {
          p <- p + ggplot2::geom_point(alpha = alpha_pt, size = pt_size, shape = pt_shape)
        } else {
          p <- p + ggplot2::geom_point(
            alpha = alpha_pt, size = pt_size, shape = pt_shape, colour = "steelblue"
          )
        }
      },
      jitter = {
        if (use_colour) {
          p <- p + ggplot2::geom_jitter(
            alpha = alpha_pt, size = pt_size, shape = pt_shape,
            width = 0.15, height = 0
          )
        } else {
          p <- p + ggplot2::geom_jitter(
            alpha = alpha_pt, size = pt_size, shape = pt_shape,
            width = 0.15, height = 0, colour = "steelblue"
          )
        }
      },
      lines = {
        if (x_binned && !is.null(sum_d) && nrow(sum_d) > 0L) {
          x_col <- if (bin_x_pos == "midpoint") "x_num" else "x"
          p <- ggplot2::ggplot(sum_d, ggplot2::aes(x = .data[[x_col]], y = .data[["y"]])) +
            ggplot2::labs(title = title, x = xlab, y = ylab) +
            ggplot2::theme_bw(base_size = 11)
          if (use_strat) {
            p <- p + ggplot2::aes(
              colour = .data[["group"]], group = .data[["group"]]
            ) +
              ggplot2::geom_line(linewidth = 0.85, alpha = 0.9) +
              ggplot2::geom_point(size = pt_size, shape = pt_shape, alpha = 0.9)
          } else {
            p <- p + ggplot2::geom_line(linewidth = 0.85, alpha = 0.9, colour = "steelblue") +
              ggplot2::geom_point(size = pt_size, shape = pt_shape, alpha = 0.9, colour = "steelblue")
          }
        } else if (line_mode == "none") {
          if (use_colour) {
            p <- p + ggplot2::geom_point(alpha = alpha_pt, size = pt_size, shape = pt_shape)
          } else {
            p <- p + ggplot2::geom_point(
              alpha = alpha_pt, size = pt_size, shape = pt_shape, colour = "steelblue"
            )
          }
        } else if (line_mode %in% c("mean", "median")) {
          line_df <- .shiny_explore_agg_line(d, how = line_mode, by_group = use_strat)
          p <- ggplot2::ggplot(line_df, ggplot2::aes(x = .data[[xp]], y = .data[[yp]])) +
            ggplot2::labs(title = title, x = xlab, y = ylab) +
            ggplot2::theme_bw(base_size = 11)
          if (use_strat) {
            p <- p + ggplot2::aes(
              colour = .data[["group"]], group = .data[["group"]]
            ) +
              ggplot2::geom_line(linewidth = 0.85, alpha = 0.9) +
              ggplot2::geom_point(size = pt_size, shape = pt_shape, alpha = 0.9)
          } else {
            p <- p + ggplot2::geom_line(linewidth = 0.85, alpha = 0.9, colour = "steelblue") +
              ggplot2::geom_point(size = pt_size, shape = pt_shape, alpha = 0.9, colour = "steelblue")
          }
        } else {
          line_df <- .shiny_explore_sort(
            d,
            sort_cols = c("panel", "line_id", xp)
          )
          p <- ggplot2::ggplot(line_df, ggplot2::aes(x = .data[[xp]], y = .data[[yp]])) +
            ggplot2::labs(title = title, x = xlab, y = ylab) +
            ggplot2::theme_bw(base_size = 11)
          if (use_colour) {
            p <- p + ggplot2::aes(
              colour = .data[[colour_col]], group = .data[["line_id"]]
            ) +
              ggplot2::geom_line(linewidth = 0.7, alpha = 0.85)
          } else {
            p <- p + ggplot2::geom_line(linewidth = 0.7, alpha = 0.85, colour = "steelblue")
          }
        }
      },
      both = {
        if (line_mode == "none") {
          line_df <- NULL
        } else if (line_mode %in% c("mean", "median")) {
          line_df <- .shiny_explore_agg_line(d, how = line_mode, by_group = use_strat)
        } else {
          line_df <- .shiny_explore_sort(
            d,
            sort_cols = c("panel", "line_id", xp)
          )
        }
        show_pts <- TRUE
        if (is.null(line_df) || nrow(line_df) == 0L) {
          line_df <- d
          show_pts <- line_mode == "none"
        }
        p <- ggplot2::ggplot(d, ggplot2::aes(x = .data[[xp]], y = .data[[yp]])) +
          ggplot2::labs(title = title, x = xlab, y = ylab) +
          ggplot2::theme_bw(base_size = 11)
        if (use_colour) {
          p <- p + ggplot2::aes(colour = .data[[colour_col]])
        }
        if (line_mode != "none" && !is.null(line_df) && nrow(line_df) > 0L) {
          if (line_mode == "individual" && use_colour) {
            p <- p + ggplot2::geom_line(
              data = line_df,
              ggplot2::aes(
                x = .data[[xp]], y = .data[[yp]],
                colour = .data[[colour_col]], group = .data[["line_id"]]
              ),
              linewidth = 0.65, alpha = 0.85, inherit.aes = FALSE
            )
          } else if (use_strat) {
            p <- p + ggplot2::geom_line(
              data = line_df,
              ggplot2::aes(
                x = .data[[xp]], y = .data[[yp]],
                colour = .data[["group"]], group = .data[["group"]]
              ),
              linewidth = 0.65, alpha = 0.85, inherit.aes = FALSE
            )
          } else {
            p <- p + ggplot2::geom_line(
              data = line_df,
              ggplot2::aes(x = .data[[xp]], y = .data[[yp]]),
              linewidth = 0.65, alpha = 0.85, colour = "steelblue",
              inherit.aes = FALSE
            )
          }
        }
        if (show_pts) {
          if (use_colour) {
            p <- p + ggplot2::geom_point(alpha = alpha_pt, size = pt_size, shape = pt_shape)
          } else {
            p <- p + ggplot2::geom_point(
              alpha = alpha_pt, size = pt_size, shape = pt_shape, colour = "steelblue"
            )
          }
        }
      },
      smooth = {
        if (x_binned && !is.null(sum_d) && nrow(sum_d) > 0L) {
          p <- ggplot2::ggplot(sum_d, ggplot2::aes(x = .data[["x_num"]], y = .data[["y"]])) +
            ggplot2::labs(title = paste(title, "(bin means)"), x = xlab, y = ylab) +
            ggplot2::theme_bw(base_size = 11)
          if (use_strat) {
            p <- p + ggplot2::aes(colour = .data[["group"]], group = .data[["group"]]) +
              ggplot2::geom_point(size = 2.2, alpha = 0.85) +
              ggplot2::geom_smooth(se = TRUE, method = "loess", formula = y ~ x,
                                   linewidth = 0.9, alpha = 0.15)
          } else {
            p <- p + ggplot2::geom_point(size = 2.2, alpha = 0.85, colour = "steelblue") +
              ggplot2::geom_smooth(se = TRUE, method = "loess", formula = y ~ x,
                                   colour = "firebrick", fill = "firebrick",
                                   linewidth = 0.9, alpha = 0.15)
          }
        } else {
          if (use_strat) {
            p <- p + ggplot2::geom_point(alpha = 0.35, size = 1.2) +
              ggplot2::geom_smooth(se = TRUE, method = "loess", formula = y ~ x,
                                   linewidth = 0.9, alpha = 0.15)
          } else {
            p <- p + ggplot2::geom_point(alpha = 0.35, size = 1.2, colour = "steelblue") +
              ggplot2::geom_smooth(se = TRUE, method = "loess", formula = y ~ x,
                                   colour = "firebrick", fill = "firebrick",
                                   linewidth = 0.9, alpha = 0.15)
          }
        }
      },
      regression = {
        if (x_binned && !is.null(sum_d) && nrow(sum_d) > 0L) {
          p <- ggplot2::ggplot(sum_d, ggplot2::aes(x = .data[["x_num"]], y = .data[["y"]])) +
            ggplot2::labs(title = paste(title, "(bin means)"), x = xlab, y = ylab) +
            ggplot2::theme_bw(base_size = 11)
          if (use_strat) {
            p <- p + ggplot2::aes(colour = .data[["group"]], group = .data[["group"]]) +
              ggplot2::geom_point(size = 2.2, alpha = 0.85) +
              ggplot2::geom_smooth(se = TRUE, method = "lm", formula = y ~ x,
                                   linewidth = 0.9, alpha = 0.15)
          } else {
            p <- p + ggplot2::geom_point(size = 2.2, alpha = 0.85, colour = "steelblue") +
              ggplot2::geom_smooth(se = TRUE, method = "lm", formula = y ~ x,
                                   colour = "firebrick", fill = "firebrick",
                                   linewidth = 0.9, alpha = 0.15)
          }
        } else {
          if (use_strat) {
            p <- p + ggplot2::geom_point(alpha = 0.35, size = 1.2) +
              ggplot2::geom_smooth(se = TRUE, method = "lm", formula = y ~ x,
                                   linewidth = 0.9, alpha = 0.15)
          } else {
            p <- p + ggplot2::geom_point(alpha = 0.35, size = 1.2, colour = "steelblue") +
              ggplot2::geom_smooth(se = TRUE, method = "lm", formula = y ~ x,
                                   colour = "firebrick", fill = "firebrick",
                                   linewidth = 0.9, alpha = 0.15)
          }
        }
      },
      boxplot = {
        if (use_strat) {
          p <- p + ggplot2::geom_boxplot(
            alpha = 0.55, outlier.size = 1.2, linewidth = 0.35,
            position = dodge, colour = "gray30"
          )
        } else {
          p <- p + ggplot2::geom_boxplot(
            alpha = 0.55, outlier.size = 1.2, linewidth = 0.35,
            fill = "steelblue", colour = "steelblue"
          )
        }
      },
      violin = {
        if (use_strat) {
          p <- p + ggplot2::geom_violin(
            alpha = 0.55, linewidth = 0.35, trim = FALSE,
            position = dodge, colour = "gray30"
          )
        } else {
          p <- p + ggplot2::geom_violin(
            alpha = 0.55, linewidth = 0.35, trim = FALSE,
            fill = "steelblue", colour = "steelblue"
          )
        }
      },
      mean_se = {
        se_fun <- function(z) {
          n <- sum(!is.na(z))
          m <- mean(z, na.rm = TRUE)
          se <- if (n > 1L) stats::sd(z, na.rm = TRUE) / sqrt(n) else 0
          data.frame(y = m, ymin = m - se, ymax = m + se)
        }
        if (use_strat) {
          p <- p + ggplot2::stat_summary(
            ggplot2::aes(group = .data[["group"]]),
            fun.data = se_fun, geom = "pointrange", size = 0.35,
            position = dodge, colour = "gray20"
          )
        } else {
          p <- p + ggplot2::stat_summary(
            fun.data = se_fun, geom = "pointrange", size = 0.35,
            colour = "steelblue"
          )
        }
        if (isTRUE(show_points)) {
          p <- .shiny_explore_add_points(
            p, d, use_strat, dodge, scatter = point_scatter, size = point_size,
            shape = point_shape
          )
        }
      },
      median_q = {
        if (is.null(q_d) || nrow(q_d) == 0L) {
          return(NULL)
        }
        q_x_col <- if (x_binned && bin_x_pos == "midpoint") "x_num" else "x"
        q_d <- q_d[order(q_d[[q_x_col]]), , drop = FALSE]
        q_linewidth <- 0.85
        med_linewidth <- 0.9
        ribbon_alpha <- max(0.05, min(0.75, as.numeric(shade_alpha %||% 0.25)))
        ribbon_fill <- "#4A90C8"
        p <- ggplot2::ggplot(q_d, ggplot2::aes(x = .data[[q_x_col]])) +
          ggplot2::labs(
            title = paste0(title, " (median + ", q_interval, "% quantiles)"),
            x = xlab, y = ylab
          ) +
          ggplot2::theme_bw(base_size = 11) +
          ggplot2::geom_ribbon(
            ggplot2::aes(ymin = .data[["y_lo"]], ymax = .data[["y_hi"]], group = 1),
            fill = ribbon_fill, alpha = ribbon_alpha, linewidth = 0, colour = NA
          ) +
          ggplot2::geom_line(
            ggplot2::aes(y = .data[["y_lo"]], group = 1),
            linewidth = q_linewidth, linetype = "dashed", colour = ribbon_fill, alpha = 0.95
          ) +
          ggplot2::geom_line(
            ggplot2::aes(y = .data[["y_hi"]], group = 1),
            linewidth = q_linewidth, linetype = "dashed", colour = ribbon_fill, alpha = 0.95
          ) +
          ggplot2::geom_line(
            ggplot2::aes(y = .data[["y"]], group = 1),
            linewidth = med_linewidth, linetype = "solid", colour = ribbon_fill
          ) +
          ggplot2::geom_point(
            ggplot2::aes(y = .data[["y"]]), size = 2, alpha = 0.9, colour = ribbon_fill
          )
        if (isTRUE(show_points)) {
          p <- .shiny_explore_add_points(
            p, d, use_strat, dodge, scatter = point_scatter,
            dodge_points = FALSE, size = point_size, shape = point_shape
          )
        }
      }
    )
    if (!is.null(colour_lbl) && plot_type %in% c(
      "points", "jitter", "smooth", "regression", "lines", "both", "mean_se"
    )) {
      p <- p + ggplot2::labs(colour = colour_lbl)
    } else if (use_strat && plot_type %in% c("boxplot", "violin")) {
      p <- p + ggplot2::labs(colour = strat, fill = strat)
    } else if (use_strat && plot_type == "median_q" && isTRUE(show_points)) {
      p <- p + ggplot2::labs(colour = strat)
    }
    return(.shiny_apply_facet(p, use_split))
  }

  # base R fallback
  if (plot_type == "boxplot" && (x_discrete || x_binned)) {
    if (use_strat) {
      graphics::boxplot(
        d$y ~ interaction(d$x, d$group, drop = TRUE),
        main = title, xlab = xlab, ylab = ylab, col = "steelblue"
      )
    } else {
      graphics::boxplot(d$y ~ d$x, main = title, xlab = xlab, ylab = ylab, col = "steelblue")
    }
    return(invisible(NULL))
  }
  xplot <- if (is.numeric(d$x)) d$x else as.numeric(d$x)
  if (length(unique(d$group)) > 1L) {
    groups <- unique(d$group)
    cols <- grDevices::rainbow(length(groups))
    graphics::plot(
      xplot, d$y, type = "n",
      main = title, xlab = xlab, ylab = ylab
    )
    for (i in seq_along(groups)) {
      sub <- d[d$group == groups[[i]], , drop = FALSE]
      xp <- if (is.numeric(sub$x)) sub$x else as.numeric(sub$x)
      if (plot_type == "lines") {
        ord <- order(xp)
        graphics::lines(xp[ord], sub$y[ord], col = cols[[i]])
      } else {
        graphics::points(xp, sub$y, col = cols[[i]], pch = 16, cex = 0.7)
      }
    }
    graphics::legend(
      "topright", legend = groups, col = cols, pch = 16, cex = 0.75,
      bty = "n", inset = c(0.02, 0.02)
    )
  } else {
    if (plot_type == "lines") {
      ord <- order(xplot)
      graphics::plot(
        xplot[ord], d$y[ord], type = "l", col = "steelblue",
        main = title, xlab = xlab, ylab = ylab
      )
    } else {
      graphics::plot(
        xplot, d$y, pch = 16, col = "steelblue", cex = 0.7,
        main = title, xlab = xlab, ylab = ylab
      )
    }
  }
  invisible(NULL)
}

.shiny_scatter_gof <- function(obs, xvar, title, xlab) {
  if (nrow(obs) == 0L) {
    return(NULL)
  }
  if (.shiny_has_ggplot()) {
    lims <- range(c(obs[[xvar]], obs$DV), finite = TRUE)
    ggplot2::ggplot(obs, ggplot2::aes(x = .data[[xvar]], y = .data[["DV"]])) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray50") +
      ggplot2::geom_point(alpha = 0.65, color = "steelblue", size = 1.5) +
      ggplot2::coord_equal(xlim = lims, ylim = lims) +
      ggplot2::labs(title = title, x = xlab, y = "DV") +
      ggplot2::theme_bw(base_size = 11)
  } else {
    lims <- range(c(obs[[xvar]], obs$DV), finite = TRUE)
    plot(obs[[xvar]], obs$DV, pch = 16, col = "steelblue", main = title,
         xlab = xlab, ylab = "DV", xlim = lims, ylim = lims)
    abline(0, 1, lty = 2, col = "gray50")
  }
}

.shiny_wres_vs <- function(obs, xvar, title, xlab) {
  ok <- !is.na(obs$WRES) & is.finite(obs$WRES)
  obs <- obs[ok, , drop = FALSE]
  if (nrow(obs) == 0L) {
    return(NULL)
  }
  if (.shiny_has_ggplot()) {
    ggplot2::ggplot(obs, ggplot2::aes(x = .data[[xvar]], y = WRES)) +
      ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "gray50") +
      ggplot2::geom_point(alpha = 0.65, color = "steelblue", size = 1.5) +
      ggplot2::labs(title = title, x = xlab, y = "WRES") +
      ggplot2::theme_bw(base_size = 11)
  } else {
    plot(obs[[xvar]], obs$WRES, pch = 16, col = "steelblue", main = title,
         xlab = xlab, ylab = "WRES")
    abline(h = 0, lty = 2, col = "gray50")
  }
}

.shiny_qq_plot <- function(x, title = "Normal Q-Q plot", xlab = "Theoretical quantiles") {
  x <- x[is.finite(x)]
  if (length(x) < 2L) {
    return(.shiny_empty_plot("Not enough finite values for Q-Q plot"))
  }
  if (.shiny_has_ggplot()) {
    d <- data.frame(sample = x)
    ggplot2::ggplot(d, ggplot2::aes(sample = sample)) +
      ggplot2::stat_qq() +
      ggplot2::stat_qq_line(linetype = 2, color = "gray40") +
      ggplot2::labs(title = title, x = xlab, y = "Sample quantiles") +
      ggplot2::theme_bw(base_size = 11)
  } else {
    stats::qqnorm(x, main = title, xlab = xlab)
    stats::qqline(x, lty = 2, col = "gray40")
    invisible(NULL)
  }
}

.shiny_gof_obs_table <- function(fit) {
  if (!is.null(fit$gof)) {
    return(as.data.frame(fit$gof))
  }
  pred <- predict(fit, type = "ipred")
  .shiny_gof_obs(pred)
}

.shiny_residual_vs <- function(obs, yvar, title, ylab) {
  ok <- !is.na(obs[[yvar]]) & is.finite(obs[[yvar]])
  obs <- obs[ok, , drop = FALSE]
  if (nrow(obs) == 0L) {
    return(NULL)
  }
  if (.shiny_has_ggplot()) {
    ggplot2::ggplot(obs, ggplot2::aes(x = TIME, y = .data[[yvar]])) +
      ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "gray50") +
      ggplot2::geom_point(alpha = 0.65, color = "steelblue", size = 1.5) +
      ggplot2::labs(title = title, x = "Time", y = ylab) +
      ggplot2::theme_bw(base_size = 11)
  } else {
    plot(obs$TIME, obs[[yvar]], pch = 16, col = "steelblue", main = title,
         xlab = "Time", ylab = ylab)
    abline(h = 0, lty = 2, col = "gray50")
  }
}

.shiny_copy_version_modal <- function(has_fit = FALSE) {
  fit_note <- if (isTRUE(has_fit)) {
    tags$p(
      class = "text-muted", style = "font-size: 11px;",
      "When a fit is loaded, you can copy its final estimates into the new version's initials."
    )
  } else {
    tags$p(
      class = "text-muted", style = "font-size: 11px;",
      "No estimation fit loaded — initials will match the source version."
    )
  }
  shiny::modalDialog(
    title = "Copy to new model version",
    fit_note,
    checkboxInput(
      "copy_update_inits",
      "Update THETA / OMEGA / SIGMA initials from current fit",
      value = isTRUE(has_fit)
    ),
    footer = tagList(
      modalButton("Cancel"),
      actionButton("confirm_copy_version", "Copy", class = "btn-primary")
    ),
    easyClose = TRUE
  )
}

.shiny_vpc_plot <- function(vpc,
                            obs = NULL,
                            fit = NULL,
                            pc_correct = FALSE,
                            title = "Visual predictive check") {
  if (is.null(vpc) || nrow(vpc) == 0L) {
    return(.shiny_empty_plot("VPC summary unavailable"))
  }
  vpc_agg_obs <- get0(
    ".nm_vpc_aggregate_obs",
    envir = asNamespace("LibeRation"),
    inherits = FALSE
  )
  vpc_df <- as.data.frame(vpc)
  if (!is.null(obs) && "TIME_BIN" %in% names(obs) && is.function(vpc_agg_obs)) {
    if (isTRUE(pc_correct) && !is.null(fit)) {
      obs_stats <- vpc_agg_obs(obs, fit = fit, pc_correct = TRUE)
      if (!is.null(obs_stats)) {
        keep <- setdiff(names(vpc_df), c("TIME", "obs_med", "obs_lo", "obs_hi", "n_obs"))
        vpc_df <- merge(vpc_df[, keep, drop = FALSE], obs_stats, by = "TIME_BIN", all.x = TRUE)
      }
    } else if (!all(c("obs_med", "obs_lo", "obs_hi") %in% names(vpc_df))) {
      obs_stats <- vpc_agg_obs(obs, fit = fit, pc_correct = FALSE)
      if (!is.null(obs_stats)) {
        vpc_df <- merge(vpc_df, obs_stats, by = "TIME_BIN", all.x = TRUE, sort = FALSE)
      }
    }
  }
  sim_suffix <- if (isTRUE(pc_correct) &&
      all(c("sim_med_lo_pc", "sim_med_hi_pc", "sim_lo_lo_pc", "sim_lo_hi_pc",
            "sim_hi_lo_pc", "sim_hi_hi_pc") %in% names(vpc_df))) {
    "_pc"
  } else {
    ""
  }
  req_cols <- c(
    "TIME", "obs_med", "obs_lo", "obs_hi",
    paste0("sim_med_lo", sim_suffix),
    paste0("sim_med_hi", sim_suffix),
    paste0("sim_lo_lo", sim_suffix),
    paste0("sim_lo_hi", sim_suffix),
    paste0("sim_hi_lo", sim_suffix),
    paste0("sim_hi_hi", sim_suffix)
  )
  if (!all(req_cols %in% names(vpc_df))) {
    msg <- if (isTRUE(pc_correct)) {
      "Prediction-corrected VPC unavailable (re-run VPC simulation or load linked fit)"
    } else {
      "VPC summary missing required columns (re-run VPC simulation)"
    }
    return(.shiny_empty_plot(msg))
  }
  vpc_df$sim_med_lo_plot <- vpc_df[[paste0("sim_med_lo", sim_suffix)]]
  vpc_df$sim_med_hi_plot <- vpc_df[[paste0("sim_med_hi", sim_suffix)]]
  vpc_df$sim_lo_lo_plot <- vpc_df[[paste0("sim_lo_lo", sim_suffix)]]
  vpc_df$sim_lo_hi_plot <- vpc_df[[paste0("sim_lo_hi", sim_suffix)]]
  vpc_df$sim_hi_lo_plot <- vpc_df[[paste0("sim_hi_lo", sim_suffix)]]
  vpc_df$sim_hi_hi_plot <- vpc_df[[paste0("sim_hi_hi", sim_suffix)]]
  vpc_df <- vpc_df[order(vpc_df$TIME, na.last = TRUE), , drop = FALSE]
  xvar <- "TIME"
  ylab <- if (isTRUE(pc_correct)) "Prediction-corrected DV" else "DV"
  if (!.shiny_has_ggplot()) {
    yrng <- range(
      c(vpc_df$obs_med, vpc_df$obs_lo, vpc_df$obs_hi,
        vpc_df$sim_med_lo_plot, vpc_df$sim_med_hi_plot,
        vpc_df$sim_lo_lo_plot, vpc_df$sim_hi_hi_plot),
      na.rm = TRUE
    )
    plot(vpc_df[[xvar]], vpc_df$obs_med, type = "n",
         ylim = yrng, main = title, xlab = "Time (binned)", ylab = ylab)
    polygon(
      c(vpc_df[[xvar]], rev(vpc_df[[xvar]])),
      c(vpc_df$sim_lo_lo_plot, rev(vpc_df$sim_lo_hi_plot)),
      col = adjustcolor("steelblue", 0.2), border = NA
    )
    polygon(
      c(vpc_df[[xvar]], rev(vpc_df[[xvar]])),
      c(vpc_df$sim_hi_lo_plot, rev(vpc_df$sim_hi_hi_plot)),
      col = adjustcolor("steelblue", 0.2), border = NA
    )
    polygon(
      c(vpc_df[[xvar]], rev(vpc_df[[xvar]])),
      c(vpc_df$sim_med_lo_plot, rev(vpc_df$sim_med_hi_plot)),
      col = adjustcolor("firebrick", 0.15), border = NA
    )
    lines(vpc_df[[xvar]], vpc_df$obs_med, col = "steelblue", lwd = 2)
    lines(vpc_df[[xvar]], vpc_df$obs_lo, lty = 2, col = "firebrick")
    lines(vpc_df[[xvar]], vpc_df$obs_hi, lty = 2, col = "firebrick")
    return(invisible(NULL))
  }
  col_blue <- "#3182BD"
  col_red <- "#DE2D26"
  ggplot2::ggplot(vpc_df, ggplot2::aes(x = .data[[xvar]])) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data[["sim_lo_lo_plot"]], ymax = .data[["sim_lo_hi_plot"]]),
      fill = col_blue, alpha = 0.22
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data[["sim_hi_lo_plot"]], ymax = .data[["sim_hi_hi_plot"]]),
      fill = col_blue, alpha = 0.22
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data[["sim_med_lo_plot"]], ymax = .data[["sim_med_hi_plot"]]),
      fill = col_red, alpha = 0.18
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data[["obs_med"]]),
      colour = col_blue, linewidth = 1.1
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data[["obs_lo"]]),
      colour = col_red, linewidth = 0.9, linetype = "dashed"
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data[["obs_hi"]]),
      colour = col_red, linewidth = 0.9, linetype = "dashed"
    ) +
    ggplot2::labs(
      title = title,
      x = "Time (binned)",
      y = ylab,
      caption = paste(
        "Blue solid: observed median; red dashed: observed 10th/90th;",
        "blue bands: 90% CI of simulated 10th/90th; red band: 90% CI of simulated median"
      )
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(plot.caption = ggplot2::element_text(size = 9, colour = "gray40"))
}

.shiny_log_max_entries <- 200L

.shiny_event_label <- function(event_expr) {
  lbl <- paste(deparse(event_expr, width.cutoff = 500L), collapse = " ")
  lbl <- sub("^input\\$", "", lbl)
  if (grepl("_rows_selected", lbl)) {
    lbl <- sub("_rows_selected.*", "", lbl)
    lbl <- paste0(lbl, " selection")
  }
  lbl
}

.shiny_empty_plot <- function(msg = "Unable to display plot") {
  graphics::plot.new()
  graphics::title(main = msg, cex.main = 0.9)
  invisible(NULL)
}

.shiny_sim_dataset <- function(sim_obj) {
  if (is.null(sim_obj)) {
    return(NULL)
  }
  if (inherits(sim_obj, "nm_dataset")) {
    return(sim_obj)
  }
  if (is.list(sim_obj) && !is.null(sim_obj$combined)) {
    return(sim_obj$combined)
  }
  if (is.list(sim_obj) && !is.null(sim_obj$primary)) {
    return(sim_obj$primary)
  }
  if (is.list(sim_obj) && !is.null(sim_obj$data)) {
    return(structure(list(data = sim_obj$data), class = "nm_dataset"))
  }
  NULL
}

.shiny_n_subjects <- function(data_obj) {
  tbl <- .shiny_dataset_table(data_obj)
  if (is.null(tbl) || !"ID" %in% names(tbl)) {
    return(10L)
  }
  as.integer(length(unique(tbl[["ID"]])))
}

.shiny_projects_tree_ui <- function(df, selected_project) {
  header <- tags$div(
    class = "panel-box-header",
    style = "padding: 4px 0 6px; border: none; background: transparent;",
    "Projects"
  )
  if (nrow(df) == 0L) {
    return(tagList(
      header,
      tags$p(
        class = "text-muted",
        style = "font-size: 11px; padding: 6px;",
        "No projects yet."
      )
    ))
  }
  rows <- lapply(seq_len(nrow(df)), function(i) {
    proj <- df$project[[i]]
    tags$div(
      class = paste("project-row", if (identical(proj, selected_project)) "selected" else NULL),
      `data-project` = proj,
      tags$span(class = "project-id", proj)
    )
  })
  tagList(header, tags$div(class = "version-tree project-tree", rows))
}

.shiny_run_estimation_modal <- function() {
  modalDialog(
    title = "Run estimation",
    size = "l",
    fluidRow(
      column(
        3L,
        selectInput(
          "method", "Method",
          choices = c("FO", "FOCE", "FOCEI", "SAEM", "LAPLACE", "IMP", "BAYES"),
          width = "100%"
        )
      ),
      column(3L, selectInput("grad", "Gradient", c("auto", "numeric", "ad", "cpp"), width = "100%")),
      column(3L, selectInput("pk_engine", "PK engine", c("auto", "cpp", "R"), width = "100%")),
      column(3L, numericInput("maxit", "maxit", value = 30L, min = 5L, step = 5L, width = "100%"))
    ),
    fluidRow(
      column(
        4L,
        conditionalPanel(
          "input.method == 'FOCE' || input.method == 'FOCEI'",
          numericInput("max_outer", "FOCE outer iterations", value = 5L, min = 1L, step = 1L, width = "100%")
        ),
        conditionalPanel(
          "input.method == 'SAEM'",
          numericInput("n_iter", "SAEM iterations", value = 30L, min = 5L, step = 5L, width = "100%"),
          numericInput("n_burn_saem", "SAEM burn-in", value = 10L, min = 0L, step = 5L, width = "100%"),
          numericInput("n_mcmc", "SAEM MCMC / subject", value = 1L, min = 1L, step = 1L, width = "100%")
        ),
        conditionalPanel(
          "input.method == 'LAPLACE'",
          numericInput("n_quad", "Gaussian quadrature points", value = 5L, min = 1L, max = 11L, step = 1L, width = "100%")
        ),
        conditionalPanel(
          "input.method == 'IMP'",
          numericInput("n_imp", "Importance samples", value = 50L, min = 5L, step = 5L, width = "100%"),
          numericInput("imp_n_quad", "IMP quadrature points", value = 5L, min = 1L, max = 11L, step = 1L, width = "100%")
        )
      ),
      column(4L, textInput("job_label", "Job label", value = "", width = "100%"))
    ),
    fluidRow(
      column(3L, numericInput("min_retries", "min_retries", value = 0L, min = 0L, max = 10L, step = 1L, width = "100%")),
      column(3L, checkboxInput("tweak_inits", "tweak_inits on retry", value = FALSE)),
      column(
        3L,
        conditionalPanel(
          "input.method != 'BAYES'",
          checkboxInput("compute_inference", "Compute standard errors", value = FALSE)
        )
      ),
      column(
        3L,
        conditionalPanel(
          "input.compute_inference && input.method != 'BAYES'",
          selectInput(
            "inference_method",
            "SE method",
            choices = c(
              "Covariance step" = "cov_step",
              "Hessian (numeric)" = "hessian_numeric",
              "Hessian (auto)" = "hessian_auto",
              "Hessian (ad)" = "hessian_ad"
            ),
            selected = "cov_step",
            width = "100%"
          )
        )
      )
    ),
    conditionalPanel(
      "input.method == 'BAYES'",
      fluidRow(
        column(4L, selectInput("sampler", "MCMC sampler", c("mh", "hmc", "nuts"), width = "100%")),
        column(4L, numericInput("n_burn", "Burn-in", value = 50L, min = 0L, step = 10L, width = "100%")),
        column(4L, numericInput("n_sample", "Samples", value = 100L, min = 10L, step = 10L, width = "100%"))
      )
    ),
    conditionalPanel(
      "input.method != 'BAYES'",
      fluidRow(
        column(
          4L,
          checkboxInput(
            "est_bootstrap", "Run bootstrap SE after estimation",
            value = FALSE
          )
        ),
        column(
          4L,
          conditionalPanel(
            "input.est_bootstrap",
            numericInput(
              "bootstrap_n", "Bootstrap replicates",
              value = 30L, min = 10L, max = 500L, step = 10L, width = "100%"
            )
          )
        ),
        column(
          4L,
          conditionalPanel(
            "input.est_bootstrap",
            numericInput(
              "bootstrap_seed", "Bootstrap seed",
              value = 1L, min = 1L, step = 1L, width = "100%"
            )
          )
        )
      )
    ),
    tags$p(
      class = "text-muted", style = "font-size: 11px;",
      "Estimation runs as a background job. ",
      "min_retries / tweak_inits follow PsN: re-start with perturbed initials when optim fails to converge. ",
      "Optional standard errors: covariance step (refit ETAs; linearized FIM or Hessian by method) or Hessian-only SE. ",
      "Bootstrap refits resampled subjects in the same worker (can take much longer than estimation alone). ",
      "BAYES SE uses posterior SD and 2.5\u201397.5% quantiles from MCMC."
    ),
    footer = tagList(
      modalButton("Cancel"),
      actionButton("submit_job", "Submit estimation", class = "btn-primary")
    ),
    easyClose = TRUE
  )
}

.shiny_sim_build_design <- function(input, seed) {
  dose_mode <- input$sim_dose_mode %||% "single"
  list(
    use_design = isTRUE(input$sim_use_design),
    n_sub = max(1L, as.integer(input$sim_n_subjects)),
    n_days = max(1L, as.integer(input$sim_n_days)),
    obs_per_day = max(3L, as.integer(input$sim_obs_per_day)),
    dose_mode = dose_mode,
    dose_amt = as.numeric(input$sim_dose_amt),
    dose_table = input$sim_dose_table,
    dose_n = max(1L, as.integer(input$sim_dose_n)),
    dose_ii = as.numeric(input$sim_dose_ii),
    dose_cmt = max(1L, as.integer(input$sim_dose_cmt)),
    seed = seed
  )
}

.shiny_versions_tree_ui <- function(df, proj, ws_root, selected_version, selected_sim,
                                     selected_est_run, expanded_versions) {
  if (nrow(df) == 0L) {
    return(tags$p(
      class = "text-muted",
      style = "font-size: 11px; padding: 6px;",
      "No model versions for this project."
    ))
  }
  rows <- lapply(seq_len(nrow(df)), function(i) {
    ver <- df$version[[i]]
    label <- df$label[[i]]
    runs_df <- nm_workspace_list_runs(proj, ver, root = ws_root)
    n_runs <- nrow(runs_df)
    expanded <- ver %in% expanded_versions
    sims <- nm_workspace_list_sims(proj, ver, root = ws_root)
    n_sims <- nrow(sims)
    standalone_sims <- if (n_sims > 0L) {
      sims[!(sims$vpc %in% TRUE & nzchar(sims$est_run_id %||% "")), , drop = FALSE]
    } else {
      sims
    }
    n_standalone <- nrow(standalone_sims)
    vpc_run_ids <- if (n_sims > 0L) {
      unique(sims$est_run_id[sims$vpc %in% TRUE & nzchar(sims$est_run_id)])
    } else {
      character()
    }
    arrow <- if (expanded) "\u25BC" else "\u25B6"
    ver_sel <- identical(ver, selected_version) &&
      !nzchar(selected_sim %||% "") &&
      !nzchar(selected_est_run %||% "")
    run_children <- if (expanded && n_runs > 0L) {
      lapply(seq_len(nrow(runs_df)), function(j) {
        rid <- runs_df$run_id[[j]]
        method <- runs_df$method[[j]]
        obj <- runs_df$objective[[j]]
        run_sel <- identical(rid, selected_est_run) && identical(ver, selected_version)
        has_vpc <- rid %in% vpc_run_ids
        has_bootstrap <- isTRUE(runs_df$has_bootstrap[[j]])
        has_npc <- isTRUE(runs_df$has_npc[[j]])
        has_npde <- isTRUE(runs_df$has_npde[[j]])
        obj_txt <- if (!is.na(obj) && is.finite(obj)) {
          paste0(" OFV=", format(round(obj, 2), nsmall = 2, trim = TRUE))
        } else {
          ""
        }
        tags$div(
          class = paste("run-row", if (run_sel) "selected" else NULL),
          `data-version` = ver,
          `data-run` = rid,
          tags$span(class = "run-id", rid),
          tags$span(class = "run-label", paste0(method, obj_txt)),
          if (has_vpc) {
            tags$span(class = "run-vpc-badge run-flag-badge", "VPC")
          },
          if (has_bootstrap) {
            tags$span(class = "run-flag-badge run-flag-bootstrap", "Bootstrap")
          },
          if (has_npc) {
            tags$span(class = "run-flag-badge run-flag-npc", "NPC")
          },
          if (has_npde) {
            tags$span(class = "run-flag-badge run-flag-npde", "NPDE")
          }
        )
      })
    } else if (expanded) {
      list(tags$div(class = "run-row run-empty", "No estimation runs yet"))
    } else {
      NULL
    }
    sim_children <- if (expanded && n_standalone > 0L) {
      lapply(seq_len(nrow(standalone_sims)), function(j) {
        sid <- standalone_sims$sim_id[[j]]
        slabel <- .shiny_meta_label(standalone_sims$label[[j]], fallback = sid)
        sim_sel <- identical(sid, selected_sim) && identical(ver, selected_version)
        tags$div(
          class = paste("sim-row", if (sim_sel) "selected" else NULL),
          `data-version` = ver,
          `data-sim` = sid,
          tags$span(class = "sim-id", sid),
          tags$span(class = "sim-label", slabel)
        )
      })
    } else if (expanded && n_runs == 0L && n_standalone == 0L) {
      list(tags$div(class = "sim-row sim-empty", "No simulations yet"))
    } else {
      NULL
    }
    child_rows <- c(
      if (!is.null(run_children)) run_children else list(),
      if (!is.null(sim_children)) sim_children else list()
    )
    tagList(
      tags$div(
        class = paste("version-row", if (ver_sel) "selected" else NULL),
        `data-version` = ver,
        tags$span(class = "version-toggle", `data-version` = ver, arrow),
        tags$span(class = "version-id", ver),
        tags$span(class = "version-label", label),
        tags$span(
          class = "version-meta",
          paste0(
            n_runs, " run", if (n_runs != 1L) "s" else "",
            if (length(vpc_run_ids) > 0L) {
              paste0(" (", length(vpc_run_ids), " VPC)")
            } else {
              ""
            },
            if (n_standalone > 0L) {
              paste0(", ", n_standalone, " sim", if (n_standalone != 1L) "s" else "")
            } else {
              ""
            }
          )
        )
      ),
      if (length(child_rows) > 0L) {
        tags$div(class = "sim-list", child_rows)
      }
    )
  })
  tags$div(class = "version-tree", rows)
}

.shiny_warn_version_artifacts <- function(project, version_id, ws_root) {
  runs <- nm_workspace_list_runs(project, version_id, root = ws_root)
  sims <- nm_workspace_list_sims(project, version_id, root = ws_root)
  if (nrow(runs) == 0L && nrow(sims) == 0L) {
    return(invisible(NULL))
  }
  parts <- character()
  if (nrow(runs) > 0L) {
    parts <- c(parts, paste0(nrow(runs), " estimation run(s)"))
  }
  if (nrow(sims) > 0L) {
    parts <- c(parts, paste0(nrow(sims), " simulation(s)"))
  }
  showModal(modalDialog(
    title = "Existing fits or simulations",
    paste0(
      "This model version already has ",
      paste(parts, collapse = " and "),
      ". Editing the control file may invalidate those results."
    ),
    footer = modalButton("OK"),
    easyClose = FALSE
  ))
  invisible(NULL)
}

.shiny_param_table_display <- function(pt, fit = NULL, include_gradient = TRUE,
                                        se_col_name = NULL) {
  if (is.null(se_col_name)) {
    se_col_name <- if (!is.null(fit) && identical(fit$method, "BAYES")) {
      "Posterior SD (95% CI)"
    } else {
      "SE (CI: 2.5% - 97.5%)"
    }
  }
  if (is.null(pt) || nrow(pt) == 0L) {
    return(data.frame(
      type = character(),
      name = character(),
      initial = character(),
      estimate = character(),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  sh <- if (!is.null(fit)) {
    tryCatch(nm_shrinkage(fit), error = function(e) NULL)
  } else {
    NULL
  }
  sh_pct <- if (!is.null(sh) && nrow(sh) > 0L) {
    stats::setNames(round(100 * sh$shrinkage, 1), sh$ETA)
  } else {
    NULL
  }
  has_initial <- "initial" %in% names(pt)
  out <- data.frame(
    type = pt$type,
    name = pt$name,
    initial = if (has_initial) {
      vapply(pt$initial, .shiny_fmt_param_num, character(1L))
    } else {
      rep("", nrow(pt))
    },
    estimate = character(nrow(pt)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out[[se_col_name]] <- character(nrow(pt))
  if (isTRUE(include_gradient)) {
    out$gradient <- character(nrow(pt))
  }
  se_col_idx <- which(names(out) == se_col_name)
  for (i in seq_len(nrow(pt))) {
    est <- .shiny_fmt_param_num(pt$estimate[[i]])
    if (identical(pt$type[[i]], "OMEGA") && !is.null(sh_pct)) {
      eta_nm <- sub("^OMEGA", "ETA", pt$name[[i]])
      if (eta_nm %in% names(sh_pct)) {
        est <- paste0(est, " [Shr. %: ", sh_pct[[eta_nm]], "]")
      }
    }
    out$estimate[[i]] <- est
    if ("se" %in% names(pt) && is.finite(pt$se[[i]])) {
      ci_txt <- if (is.finite(pt$ci_low[[i]]) && is.finite(pt$ci_high[[i]])) {
        paste0(
          " (CI: ",
          .shiny_fmt_sig3(pt$ci_low[[i]]),
          " - ",
          .shiny_fmt_sig3(pt$ci_high[[i]]),
          ")"
        )
      } else {
        ""
      }
      se_txt <- if (pt$se[[i]] == 0 &&
        "gradient" %in% names(pt) &&
        is.finite(pt$gradient[[i]]) &&
        abs(pt$gradient[[i]]) > 0.01) {
        paste0("0 (singular Hessian)", ci_txt)
      } else {
        paste0(.shiny_fmt_sig3(pt$se[[i]]), ci_txt)
      }
      out[[se_col_idx]][[i]] <- se_txt
    }
    if (isTRUE(include_gradient) && "gradient" %in% names(pt)) {
      g <- pt$gradient[[i]]
      if (is.finite(g) && abs(g) < 1e-8) {
        out$gradient[[i]] <- sprintf(
          '<span style="color:#c0392b;font-weight:600;">%s</span>',
          htmltools::htmlEscape(.shiny_fmt_param_num(g))
        )
      } else {
        out$gradient[[i]] <- .shiny_fmt_param_num(g)
      }
    }
  }
  if (!has_initial) {
    out$initial <- NULL
  }
  if (!isTRUE(include_gradient) && "gradient" %in% names(out)) {
    out$gradient <- NULL
  }
  out
}

.shiny_fmt_param_num <- function(x, digits = 5L) {
  if (length(x) == 0L || is.na(x) || !is.finite(x)) {
    return("")
  }
  x <- as.numeric(x)
  if (x != 0 && abs(x) < 0.0001) {
    return(format(x, digits = 3, scientific = TRUE, trim = TRUE))
  }
  format(round(x, digits), nsmall = min(2L, digits), trim = TRUE)
}

.shiny_fmt_sig3 <- function(x) {
  if (length(x) == 0L || is.na(x) || !is.finite(x)) {
    return("")
  }
  x <- as.numeric(x)
  if (x == 0) {
    return("0")
  }
  format(signif(x, 3), scientific = abs(x) < 1e-3 | abs(x) >= 1e4, trim = TRUE)
}

.shiny_fit_metric_label <- function(fit) {
  if (!is.null(fit) && identical(fit$method, "BAYES")) {
    "Log posterior"
  } else {
    "Objective (OFV)"
  }
}

.shiny_fit_metric_value <- function(fit) {
  if (is.null(fit)) {
    return("\u2014")
  }
  if (identical(fit$method, "BAYES")) {
    lp <- fit$log_posterior
    if (!is.null(lp) && is.finite(lp)) {
      return(round(lp, 4))
    }
    return("\u2014")
  }
  obj <- fit$objective
  if (!is.null(obj) && is.finite(obj)) {
    return(round(obj, 4))
  }
  "\u2014"
}

.shiny_param_table_for_fit <- function(model, fit) {
  compute_se <- is.null(fit$par_se)
  if (is.null(fit$model) && !is.null(model)) {
    fit$model <- model
  }
  LibeRation::nm_fit_param_table(model, fit, compute_se = compute_se)
}

.shiny_append_bootstrap_se <- function(display, fit) {
  if (is.null(display) || nrow(display) == 0L || is.null(fit) || is.null(fit$bootstrap)) {
    return(display)
  }
  bs <- fit$bootstrap$se
  if (is.null(bs)) {
    return(display)
  }
  n_ok <- as.integer(fit$bootstrap$n_ok %||% 0L)
  n_ok_col <- fit$bootstrap$n_ok_col
  display$`Bootstrap SE` <- ""
  for (i in seq_len(nrow(display))) {
    nm <- display$name[[i]]
    if (!nm %in% names(bs)) {
      next
    }
    n_rep <- if (!is.null(n_ok_col) && nm %in% names(n_ok_col)) {
      as.integer(n_ok_col[[nm]])
    } else {
      n_ok
    }
    if (n_rep < 2L) {
      display$`Bootstrap SE`[[i]] <- paste0("(", n_rep, " ok)")
    } else if (is.finite(bs[[nm]])) {
      display$`Bootstrap SE`[[i]] <- .shiny_fmt_sig3(bs[[nm]])
    }
  }
  display
}

.shiny_fit_has_npc <- function(fit) {
  fn <- get0(".nm_fit_has_npc", envir = asNamespace("LibeRation"), inherits = FALSE)
  if (is.function(fn)) {
    return(fn(fit))
  }
  FALSE
}

.shiny_fit_has_npde <- function(fit) {
  fn <- get0(".nm_fit_has_npde", envir = asNamespace("LibeRation"), inherits = FALSE)
  if (is.function(fn)) {
    return(fn(fit))
  }
  FALSE
}

.shiny_job_param_detail_table <- function(pt, fit = NULL, include_gradient = TRUE) {
  display <- .shiny_param_table_display(pt, fit = fit, include_gradient = include_gradient)
  display <- .shiny_append_bootstrap_se(display, fit)
  if (nrow(display) == 0L) {
    return(tags$p(class = "text-muted", style = "font-size: 11px;", "No parameters."))
  }
  cols <- names(display)
  header <- tags$tr(lapply(cols, tags$th))
  body_rows <- lapply(seq_len(nrow(display)), function(i) {
    tags$tr(lapply(cols, function(col) {
      val <- display[[col]][[i]]
      if (identical(col, "gradient") && grepl("<span", val, fixed = TRUE)) {
        tags$td(HTML(val))
      } else {
        tags$td(val)
      }
    }))
  })
  tags$table(
    class = "job-detail-table right-param-table",
    tags$thead(header),
    tags$tbody(body_rows)
  )
}

.shiny_jobs_tree_ui <- function(df, selected_job, expanded_jobs, job_root) {
  if (nrow(df) == 0L) {
    return(tags$p(
      class = "text-muted",
      style = "font-size: 11px; padding: 6px;",
      "No jobs yet."
    ))
  }
  rows <- lapply(seq_len(nrow(df)), function(i) {
    jid <- df$id[[i]]
    jtype <- df$job_type[[i]]
    is_sim <- identical(jtype, "sim")
    status <- df$status[[i]]
    can_expand <- !is_sim && status %in% c("success", "error")
    expanded <- jid %in% expanded_jobs
    arrow <- if (can_expand) {
      if (expanded) "\u25BC" else "\u25B6"
    } else {
      ""
    }
    sel <- identical(jid, selected_job)
    label <- as.character(df$label[[i]])
    if (is.na(label) || !nzchar(label)) {
      label <- jid
    }
    type_lbl <- if (is_sim) "simulation" else "estimation"
    meta <- if (is_sim) {
      sid <- as.character(df$sim_id[[i]])
      ns <- df$n_sim[[i]]
      paste0(
        if (!is.na(sid) && nzchar(sid) && !identical(sid, "NA")) paste0("sim: ", sid) else "",
        if (!is.na(ns)) {
          paste0(if (nzchar(sid) && !identical(sid, "NA")) ", " else "", ns, " rep", if (ns != 1L) "s" else "")
        } else {
          ""
        }
      )
    } else {
      as.character(df$method[[i]])
    }
    started <- .shiny_format_job_time(df$started[[i]])
    finished <- .shiny_format_job_time(df$finished[[i]])
    duration_lbl <- .shiny_job_duration_label(df$started[[i]], df$finished[[i]], status)
    detail <- if (can_expand && expanded) {
      if (identical(status, "error")) {
        st_full <- tryCatch(nm_job_status(jid, job_root), error = function(e) NULL)
        err_txt <- if (!is.null(st_full$error)) as.character(st_full$error) else ""
        if (!nzchar(trimws(err_txt))) {
          err_txt <- "Estimation failed (see worker log)."
        }
        tags$div(
          class = "job-detail-panel job-detail-error",
          tags$p(
            style = "font-size: 11px; margin: 4px 0 6px; color: #c0392b;",
            strong("Error:"),
            htmltools::htmlEscape(err_txt)
          )
        )
      } else {
      fit <- tryCatch(nm_job_result(jid, job_root), error = function(e) NULL)
      if (is.null(fit)) {
        tags$div(class = "job-detail-panel", tags$p(class = "text-muted", "Could not load fit result."))
      } else {
        pt <- tryCatch(
          .shiny_param_table_for_fit(fit$model, fit),
          error = function(e) NULL
        )
        tags$div(
          class = "job-detail-panel",
          tags$p(
            style = "font-size: 11px; margin: 4px 0 6px;",
            strong(.shiny_fit_metric_label(fit), ":"),
            .shiny_fit_metric_value(fit)
          ),
          if (!is.null(pt)) {
            .shiny_job_param_detail_table(pt, fit = fit, include_gradient = FALSE)
          } else {
            tags$p(class = "text-muted", "Parameter details unavailable.")
          }
        )
      }
      }
    } else {
      NULL
    }
    tagList(
      tags$div(
        class = paste(
          "job-row", paste0("job-row-", status),
          if (sel) "selected" else NULL
        ),
        `data-job` = jid,
        if (can_expand) {
          tags$span(class = "job-toggle", `data-job` = jid, arrow)
        } else {
          tags$span(class = "job-toggle-spacer")
        },
        tags$span(class = "job-id", jid),
        tags$span(class = "job-label", label),
        tags$span(class = paste("job-status", paste0("job-status-", status)), status),
        tags$span(class = "job-type", type_lbl),
        tags$span(class = "job-meta", meta),
        tags$span(
          class = "job-time",
          paste0(
            if (nzchar(started)) paste0("started: ", started) else "",
            if (nzchar(duration_lbl)) {
              paste0(if (nzchar(started)) " \u00b7 " else "", duration_lbl)
            } else if (nzchar(finished)) {
              paste0(if (nzchar(started)) " \u00b7 finished: " else "finished: ", finished)
            } else {
              ""
            }
          )
        )
      ),
      if (!is.null(detail)) detail
    )
  })
  tags$div(class = "version-tree job-tree", rows)
}

.shiny_run_workspace_simulation <- function(project, version_id, label, seed, n_sim,
                                            use_fit = FALSE, est_run_id = NULL,
                                            root = nm_workspace_root()) {
  parsed <- nm_workspace_parse_model(project, version_id, root = root)
  if (is.null(parsed$model) || is.null(parsed$data)) {
    stop("Could not parse model or dataset for this version.")
  }
  theta <- omega <- sigma <- NULL
  src_run <- NULL
  if (isTRUE(use_fit)) {
    fit <- if (!is.null(est_run_id) && nzchar(est_run_id)) {
      nm_workspace_load_run_fit(project, version_id, est_run_id, root = root)
    } else {
      nm_workspace_load_fit(project, version_id, root = root)
    }
    if (!is.null(fit)) {
      theta <- fit$theta
      omega <- fit$omega
      sigma <- fit$sigma
      src_run <- est_run_id
    }
  }
  n_sim <- max(1L, as.integer(n_sim))
  seed <- as.integer(seed)
  sim_args <- list(theta = theta, omega = omega, sigma = sigma)
  sim_args <- sim_args[!vapply(sim_args, is.null, logical(1L))]

  sim_id <- nm_workspace_new_sim_id(project, version_id, root = root)
  if (n_sim <= 1L) {
    sim_dat <- do.call(
      nm_task,
      c(list("sim", parsed$model, parsed$data, seed = seed), sim_args)
    )
    sim_out <- structure(list(data = sim_dat), class = "nm_dataset")
  } else {
    sims <- do.call(
      nm_simulate,
      c(list(parsed$model, parsed$data, n_sim = n_sim, seed = seed), sim_args)
    )
    sim_out <- LibeRation:::.nm_sim_pack_output(sims)
  }
  nm_workspace_save_sim(
    project, version_id, sim_id, sim_out,
    root = root,
    label = if (nzchar(label)) label else NULL,
    seed = seed,
    n_sim = n_sim,
    use_fit = isTRUE(use_fit),
    est_run_id = src_run
  )
  sim_id
}

.shiny_sim_workload_threshold <- 500000L

.shiny_clean_job_text <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return("")
  }
  txt <- as.character(x)[1L]
  if (is.na(txt) || !nzchar(txt) || identical(txt, "NA")) {
    ""
  } else {
    txt
  }
}

.shiny_default_n_cores <- function() {
  nc <- parallel::detectCores(logical = TRUE)
  if (is.na(nc) || nc < 2L) {
    1L
  } else {
    max(1L, nc - 1L)
  }
}

.shiny_sim_dataset_key <- function(version_id, sim_id) {
  paste0("__sim__/", version_id, "/", sim_id)
}

.shiny_est_dataset_key <- function(version_id, run_id) {
  paste0("__est__/", version_id, "/", run_id)
}

.shiny_model_dataset_key <- function(version_id) {
  paste0("__model__/", version_id)
}

.shiny_fit_prediction_table <- function(fit) {
  if (is.null(fit) || is.null(fit$model)) {
    return(NULL)
  }
  pred <- tryCatch(predict(fit, type = "ipred"), error = function(e) NULL)
  if (is.null(pred)) {
    return(NULL)
  }
  as.data.frame(pred)
}

.shiny_ensure_fit_eta <- function(fit) {
  if (is.null(fit) || is.null(fit$model)) {
    return(fit)
  }
  n_eta <- LibeRation:::.nm_n_eta(fit$model)
  if (n_eta == 0L) {
    return(fit)
  }
  eta_ok <- is.matrix(fit$eta) &&
    nrow(fit$eta) > 0L &&
    ncol(fit$eta) >= n_eta &&
    any(is.finite(fit$eta)) &&
    max(abs(fit$eta), na.rm = TRUE) > 1e-10
  if (isTRUE(eta_ok)) {
    return(fit)
  }
  fit$eta <- LibeRation:::.nm_fit_eta_matrix(fit, data = fit$data)
  fit
}

.shiny_explore_dataset_choices <- function(project, root) {
  groups <- list()
  vers <- nm_workspace_list_versions(project, root = root)
  if (length(vers) > 0L) {
    mdl <- character()
    for (ver in vers) {
      dp <- nm_workspace_model_data_path(project, ver, root = root)
      lbl <- if (!is.null(dp)) {
        paste0(ver, " — ", basename(dp))
      } else {
        paste0(ver, " (model input)")
      }
      mdl[lbl] <- .shiny_model_dataset_key(ver)
    }
    groups[["Model input datasets"]] <- mdl
  }
  ds <- nm_workspace_list_datasets(project, root = root)
  if (nrow(ds) > 0L) {
    groups[["Project datasets"]] <- stats::setNames(ds$path, ds$file)
  }
  runs <- nm_workspace_list_project_runs(project, root = root)
  if (nrow(runs) > 0L) {
    run_labels <- paste0(
      runs$run_id, " — ", runs$version_id,
      vapply(runs$label, function(l) {
        if (is.null(l) || !nzchar(l)) "" else paste0(" (", l, ")")
      }, character(1L)),
      " [", runs$method, "]"
    )
    run_vals <- mapply(
      .shiny_est_dataset_key,
      runs$version_id,
      runs$run_id,
      USE.NAMES = FALSE
    )
    groups[["Estimation outputs"]] <- stats::setNames(run_vals, run_labels)
  }
  sims <- nm_workspace_list_project_sims(project, root = root)
  if (nrow(sims) > 0L) {
    sim_labels <- paste0(
      sims$sim_id, " — ", sims$version_id,
      vapply(sims$label, function(l) {
        if (is.null(l) || !nzchar(l)) "" else paste0(" (", l, ")")
      }, character(1L)),
      " [", sims$n_sim, " rep]"
    )
    sim_vals <- mapply(
      .shiny_sim_dataset_key,
      sims$version_id,
      sims$sim_id,
      USE.NAMES = FALSE
    )
    groups[["Simulation outputs"]] <- stats::setNames(sim_vals, sim_labels)
  }
  if (length(groups) == 0L) {
    return(c("(model version dataset)" = ""))
  }
  c("(model version dataset)" = "", groups)
}

.shiny_sim_primary_df <- function(sim_obj) {
  if (is.null(sim_obj)) {
    return(NULL)
  }
  ds <- .shiny_sim_dataset(sim_obj)
  if (is.null(ds)) {
    return(NULL)
  }
  as.data.frame(ds$data)
}

.shiny_start_simulation_job <- function(session, state, proj, ver, label, seed, n_sim,
                                        use_fit, est_run, n_cores, design, vpc,
                                        sim_compute_npc = FALSE, sim_compute_npde = FALSE,
                                        diag_n_sim = 50L, diag_refit_eta = TRUE,
                                        ws_root, job_root) {
  parsed <- nm_workspace_parse_model(proj, ver, root = ws_root)
  if (is.null(parsed$model) || is.null(parsed$data)) {
    stop("Could not parse model or dataset for this version.")
  }
  if (is.null(design)) {
    design <- list(
      n_sub = .shiny_n_subjects(parsed$data),
      use_design = FALSE,
      seed = seed
    )
  } else {
    design$seed <- seed
  }
  template_data <- nm_sim_template_data(parsed$model, parsed$data, design)
  theta <- omega <- sigma <- NULL
  src_run <- NULL
  if (isTRUE(use_fit)) {
    fit <- if (!is.null(est_run) && nzchar(est_run)) {
      nm_workspace_load_run_fit(proj, ver, est_run, root = ws_root)
    } else {
      nm_workspace_load_fit(proj, ver, root = ws_root)
    }
    if (!is.null(fit)) {
      theta <- fit$theta
      omega <- fit$omega
      sigma <- fit$sigma
      src_run <- est_run
    }
  }
  n_sim <- max(1L, as.integer(n_sim))
  seed <- as.integer(seed)
  n_cores <- max(1L, as.integer(n_cores))
  run_diag <- (isTRUE(sim_compute_npc) || isTRUE(sim_compute_npde)) && isTRUE(use_fit)
  diag_only <- run_diag && !isTRUE(vpc) && n_sim <= 1L
  sim_id <- if (diag_only) {
    paste0("_diag_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  } else {
    nm_workspace_new_sim_id(proj, ver, root = ws_root)
  }
  if (isTRUE(use_fit) && is.null(src_run)) {
    runs <- nm_workspace_list_runs(proj, ver, root = ws_root)
    if (nrow(runs) > 0L) {
      src_run <- runs$run_id[[1L]]
    }
  }
  job <- nm_job_submit_sim(
    model = parsed$model,
    data = template_data,
    project = proj,
    version_id = ver,
    sim_id = sim_id,
    n_sim = n_sim,
    seed = seed,
    n_cores = n_cores,
    pk_engine = "cpp",
    theta = theta,
    omega = omega,
    sigma = sigma,
    label = if (nzchar(label)) label else NULL,
    use_fit = isTRUE(use_fit),
    est_run_id = src_run,
    design = design,
    workspace_root = ws_root,
    job_root = job_root,
    vpc = isTRUE(vpc),
    sim_compute_npc = isTRUE(sim_compute_npc) && isTRUE(use_fit),
    sim_compute_npde = isTRUE(sim_compute_npde) && isTRUE(use_fit),
    diag_n_sim = max(1L, as.integer(diag_n_sim)),
    diag_refit_eta = isTRUE(diag_refit_eta),
    diag_only = diag_only
  )
  state$handles[[job$id]] <- job$process
  state$selected_job <- job$id
  state$active_job_type <- "sim"
  state$active_job_project <- proj
  state$active_job_version <- ver
  state$active_job_sim_id <- sim_id
  state$active_job_est_run <- if (isTRUE(use_fit) && !is.null(src_run) && nzchar(src_run)) {
    src_run
  } else {
    NULL
  }
  state$active_job_diag_only <- diag_only
  state$jobs_rev <- state$jobs_rev + 1L
  message(
    "[LibeRation] Simulation started: job=", job$id,
    "  sim_id=", sim_id,
    "  replicates=", n_sim
  )
  showNotification(
    paste("Simulation job submitted:", job$id),
    type = "message",
    duration = 5L
  )
  if (run_diag) {
    parts <- c(
      if (isTRUE(sim_compute_npc)) "NPC",
      if (isTRUE(sim_compute_npde)) "NPDE"
    )
    showNotification(
      paste(
        paste(parts, collapse = "/"),
        "will run on the linked estimation run",
        if (diag_only) "(no simulation entry)." else "after simulation.",
        sep = " "
      ),
      type = "message",
      duration = 8L
    )
  }
  invisible(list(job_id = job$id, sim_id = sim_id))
}

.shiny_std_colnames <- function(nms) {
  if (is.null(nms) || length(nms) == 0L) {
    return(nms)
  }
  out <- toupper(nms)
  ce <- which(out == "CEVID")
  if (length(ce)) {
    out[ce] <- "cEVID"
  }
  out
}

.shiny_as_explore_df <- function(dat) {
  if (is.null(dat)) {
    return(NULL)
  }
  df <- as.data.frame(dat)
  names(df) <- .shiny_std_colnames(names(df))
  df
}

.shiny_format_job_time <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return("")
  }
  txt <- as.character(x)[1L]
  if (is.na(txt) || !nzchar(txt)) {
    return("")
  }
  pt <- suppressWarnings(as.POSIXct(txt))
  if (is.na(pt)) {
    return(txt)
  }
  format(pt, "%Y-%m-%d %H:%M:%S")
}

.shiny_fmt_duration <- function(secs) {
  secs <- as.numeric(secs)
  if (!is.finite(secs) || secs < 0) {
    return("")
  }
  if (secs < 60) {
    return(paste0(round(secs, 1), " s"))
  }
  if (secs < 3600) {
    return(paste0(round(secs / 60, 1), " min"))
  }
  paste0(round(secs / 3600, 2), " h")
}

.shiny_job_duration_label <- function(started, finished, status = NULL) {
  if (is.null(started) || !nzchar(as.character(started)[1L])) {
    return("")
  }
  t0 <- suppressWarnings(as.POSIXct(started))
  if (is.na(t0)) {
    return("")
  }
  running <- identical(status, "running") || identical(status, "queued")
  t1 <- if (!is.null(finished) && nzchar(as.character(finished)[1L])) {
    suppressWarnings(as.POSIXct(finished))
  } else if (running) {
    Sys.time()
  } else {
    NA
  }
  if (is.na(t1)) {
    return("")
  }
  secs <- as.numeric(difftime(t1, t0, units = "secs"))
  lbl <- .shiny_fmt_duration(secs)
  if (!nzchar(lbl)) {
    return("")
  }
  if (running && (is.null(finished) || !nzchar(as.character(finished)[1L]))) {
    paste0("elapsed: ", lbl)
  } else {
    paste0("duration: ", lbl)
  }
}

.shiny_run_compare_key <- function(version_id, run_id) {
  paste(version_id, run_id, sep = "|||")
}

.shiny_parse_run_compare_key <- function(key) {
  parts <- strsplit(key, "\\|\\|\\|", fixed = TRUE)[[1L]]
  list(version_id = parts[[1L]], run_id = parts[[2L]])
}

.shiny_project_run_choices <- function(project, ws_root) {
  if (!.valid_id(project)) {
    return(stats::setNames(character(), character()))
  }
  vers <- nm_workspace_list_versions(project, root = ws_root)
  keys <- character()
  labels <- character()
  for (ver in vers) {
    runs <- nm_workspace_list_runs(project, ver, root = ws_root)
    if (nrow(runs) == 0L) {
      next
    }
    for (i in seq_len(nrow(runs))) {
      rid <- runs$run_id[[i]]
      key <- .shiny_run_compare_key(ver, rid)
      method <- runs$method[[i]]
      obj <- runs$objective[[i]]
      obj_txt <- if (!is.na(obj) && is.finite(obj)) {
        paste0(", OFV=", format(round(obj, 2), nsmall = 2, trim = TRUE))
      } else {
        ""
      }
      keys <- c(keys, key)
      labels <- c(labels, paste0(ver, " / ", rid, " (", method, obj_txt, ")"))
    }
  }
  stats::setNames(keys, labels)
}

.shiny_compare_picker_modal <- function(choices) {
  empty <- length(choices) == 0L
  modalDialog(
    title = "Compare model runs",
    size = if (empty) "s" else "m",
    if (empty) {
      tags$p(
        class = "text-muted",
        "No estimation runs in this project yet. Run estimation on one or more model versions first."
      )
    } else {
      tagList(
        tags$p(class = "text-muted", style = "font-size: 12px;",
               "Choose two estimation runs (from any model versions) to compare side by side."),
        selectInput("compare_pick_a", "Run 1", choices = c("Select a run…" = "", choices),
                    width = "100%"),
        selectInput("compare_pick_b", "Run 2", choices = c("Select a run…" = "", choices),
                    width = "100%")
      )
    },
    footer = if (empty) {
      modalButton("Close")
    } else {
      tagList(
        modalButton("Cancel"),
        actionButton("compare_pick_submit", "Compare", class = "btn-primary")
      )
    },
    easyClose = TRUE
  )
}

.shiny_load_compare_entries <- function(project, keys, ws_root) {
  if (length(keys) < 2L) {
    return(list())
  }
  entries <- lapply(keys, function(key) {
    parts <- .shiny_parse_run_compare_key(key)
    ver <- parts$version_id
    rid <- parts$run_id
    fit <- tryCatch(
      nm_workspace_load_run_fit(project, ver, rid, root = ws_root),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(NULL)
    }
    list(
      key = key,
      label = paste0(ver, " / ", rid),
      version_id = ver,
      run_id = rid,
      fit = fit
    )
  })
  Filter(Negate(is.null), entries)
}

.shiny_synthetic_tpl_choices <- function() {
  cat <- nm_synthetic_catalog()
  stats::setNames(vapply(cat, `[[`, character(1L), "id"), vapply(cat, `[[`, character(1L), "label"))
}

.shiny_synthetic_tpl_description <- function(id) {
  cat <- nm_synthetic_catalog()
  if (!id %in% names(cat)) {
    return("")
  }
  cat[[id]]$description
}

.shiny_compare_runs_modal_body <- function(entries) {
  if (length(entries) < 2L) {
    return(tags$p(class = "text-muted", "Select at least two runs to compare."))
  }
  param_blocks <- lapply(entries, function(ent) {
    pt <- tryCatch(
      .shiny_param_table_for_fit(ent$fit$model, ent$fit),
      error = function(e) NULL
    )
    disp <- if (!is.null(pt)) .shiny_param_table_display(pt, fit = ent$fit, include_gradient = FALSE) else NULL
    tags$div(
      class = "compare-run-col",
      tags$h5(ent$label),
      tags$p(class = "text-muted", style = "font-size: 11px;",
             paste0("Method: ", ent$fit$method,
                    " | ", .shiny_fit_metric_label(ent$fit), ": ",
                    .shiny_fit_metric_value(ent$fit))),
      if (!is.null(disp) && nrow(disp) > 0L) {
        tags$table(
          class = "table table-condensed compare-param-table",
          tags$thead(tags$tr(lapply(names(disp), tags$th))),
          tags$tbody(lapply(seq_len(nrow(disp)), function(i) {
            tags$tr(lapply(names(disp), function(col) tags$td(disp[[col]][[i]])))
          }))
        )
      } else {
        tags$p(class = "text-muted", "Parameters unavailable.")
      }
    )
  })
  tagList(
    tags$div(class = "compare-runs-grid", param_blocks),
    tags$hr(),
    tags$h5("Goodness of fit"),
    plotOutput("compare_gof_grid", height = "420px")
  )
}

.shiny_meta_label <- function(label, fallback = "") {
  if (is.null(label) || length(label) == 0L) {
    return(fallback)
  }
  txt <- trimws(as.character(label)[1L])
  if (is.na(txt) || !nzchar(txt)) fallback else txt
}

.shiny_dataset_table <- function(obj) {
  if (is.null(obj)) {
    return(NULL)
  }
  if (inherits(obj, "nm_dataset")) {
    dat <- obj$data
    if (is.data.frame(dat) || inherits(dat, "data.table")) {
      return(dat)
    }
    return(NULL)
  }
  if (is.data.frame(obj) || inherits(obj, "data.table")) {
    return(obj)
  }
  if (is.list(obj) && !is.null(obj$data)) {
    if (is.data.frame(obj$data) || inherits(obj$data, "data.table")) {
      return(obj$data)
    }
    if (inherits(obj$data, "nm_dataset")) {
      return(.shiny_dataset_table(obj$data))
    }
  }
  NULL
}

.shiny_advan_choice_labels <- function() {
  c(
    "1" = "ADVAN1 — IV, 1-compartment",
    "2" = "ADVAN2 — oral, 1-compartment",
    "3" = "ADVAN3 — IV, 2-compartment",
    "4" = "ADVAN4 — oral, 2-compartment",
    "6" = "ADVAN6 — general ODE",
    "10" = "ADVAN10 — Michaelis-Menten",
    "11" = "ADVAN11 — IV, 3-compartment",
    "12" = "ADVAN12 — oral, 3-compartment",
    "13" = "ADVAN13 — general ODE (NM7+)"
  )
}

.shiny_ctl_model_help_body <- function(advan, trans = NULL) {
  info <- nm_ctl_model_info(advan, trans)
  diagram <- .shiny_pk_diagram(info)
  param_tbl <- if (nrow(info$parameters) > 0L) {
    tags$table(
      class = "table table-condensed pk-param-table",
      style = "font-size: 12px; margin-top: 8px;",
      tags$thead(tags$tr(tags$th("Symbol"), tags$th("Meaning"))),
      tags$tbody(lapply(seq_len(nrow(info$parameters)), function(i) {
        tags$tr(
          tags$td(tags$code(info$parameters$symbol[[i]])),
          tags$td(info$parameters$meaning[[i]])
        )
      }))
    )
  } else {
    NULL
  }
  tagList(
    tags$p(style = "font-size: 12px; margin-bottom: 6px;", tags$strong("Route: "), info$route),
    tags$p(style = "font-size: 12px;", info$summary),
    diagram,
    if (!is.null(param_tbl)) tagList(tags$h5("Parameters"), param_tbl)
  )
}

.shiny_pk_flow_label <- function(flows, from, to) {
  hit <- Filter(function(f) identical(f$from, from) && identical(f$to, to), flows)
  if (length(hit) == 0L) {
    return("Q")
  }
  hit[[1L]]$label
}

.shiny_pk_diagram <- function(info) {
  ncomp <- length(info$compartments)
  if (ncomp == 0L) {
    return(NULL)
  }
  boxes <- lapply(seq_len(ncomp), function(i) {
    tags$div(class = "pk-box", info$compartments[[i]])
  })
  flows <- info$flows %||% list()
  has_depot <- grepl("[Dd]epot", info$compartments[[1L]])
  central_idx <- if (has_depot) 2L else 1L
  peripheral_idxs <- setdiff(seq_len(ncomp), c(if (has_depot) 1L else integer(), central_idx))
  elim <- Filter(function(f) is.null(f$to), flows)

  if (length(peripheral_idxs) <= 1L) {
    parts <- list()
    for (i in seq_len(ncomp)) {
      parts[[length(parts) + 1L]] <- boxes[[i]]
      if (i < ncomp) {
        fwd <- .shiny_pk_flow_label(flows, i, i + 1L)
        parts[[length(parts) + 1L]] <- tags$div(class = "pk-arrow", paste0("\u2192 ", fwd, " \u2192"))
      }
    }
    if (length(elim) > 0L && central_idx <= ncomp) {
      parts[[length(parts) + 1L]] <- tags$div(class = "pk-elim", paste0("\u2193 ", elim[[1L]]$label))
    }
    return(tags$div(class = "pk-diagram", parts))
  }

  depot_section <- if (has_depot) {
    ka_lab <- .shiny_pk_flow_label(flows, 1L, 2L)
    tagList(
      boxes[[1L]],
      tags$div(class = "pk-arrow", paste0("\u2192 ", ka_lab, " \u2192"))
    )
  } else {
    NULL
  }
  elim_html <- if (length(elim) > 0L) {
    tags$div(class = "pk-elim", paste0("\u2193 ", elim[[1L]]$label))
  } else {
    NULL
  }
  peripheral_rows <- lapply(peripheral_idxs, function(pi) {
    fwd <- .shiny_pk_flow_label(flows, central_idx, pi)
    tags$div(
      class = "pk-branch",
      tags$div(class = "pk-arrow-v", paste0("\u2194 ", fwd)),
      boxes[[pi]]
    )
  })
  tags$div(
    class = "pk-diagram pk-diagram-hub",
    if (!is.null(depot_section)) tags$div(class = "pk-depot-row", depot_section),
    tags$div(
      class = "pk-hub-row",
      tags$div(class = "pk-central-col", boxes[[central_idx]], elim_html),
      tags$div(class = "pk-peripheral-col", peripheral_rows)
    )
  )
}

.shiny_has_successful_fit <- function(fit) {
  if (is.null(fit) || is.null(fit$method)) {
    return(FALSE)
  }
  if (identical(fit$method, "BAYES")) {
    ch <- fit$chains
    return(!is.null(ch) && is.matrix(ch$theta) && nrow(ch$theta) > 0L)
  }
  if (!is.null(fit$par) && any(is.finite(unname(fit$par)))) {
    obj <- fit$objective
    if (!is.null(obj) && is.finite(obj)) {
      return(TRUE)
    }
    if (!is.null(fit$theta) && any(is.finite(fit$theta))) {
      return(TRUE)
    }
  }
  FALSE
}

.shiny_load_run_vpc <- function(project, version_id, est_run_id, ws_root) {
  if (!.valid_id(project) || !.valid_id(version_id) || !.valid_id(est_run_id)) {
    return(NULL)
  }
  fn <- get0("nm_workspace_load_run_vpc", envir = asNamespace("LibeRation"), inherits = FALSE)
  if (!is.null(fn)) {
    return(fn(project, version_id, est_run_id, root = ws_root))
  }
  sims <- tryCatch(
    nm_workspace_list_sims(project, version_id, root = ws_root),
    error = function(e) NULL
  )
  if (is.null(sims) || nrow(sims) == 0L ||
      !all(c("vpc", "est_run_id") %in% names(sims))) {
    return(NULL)
  }
  hits <- sims[sims$vpc %in% TRUE & sims$est_run_id == est_run_id, , drop = FALSE]
  if (nrow(hits) == 0L) {
    return(NULL)
  }
  hits <- hits[order(hits$created, hits$sim_id, decreasing = TRUE), , drop = FALSE]
  sim_id <- hits$sim_id[[1L]]
  sim_obj <- nm_workspace_load_sim(project, version_id, sim_id, root = ws_root)
  if (is.null(sim_obj) || !is.list(sim_obj) || !isTRUE(sim_obj$vpc_mode)) {
    return(NULL)
  }
  list(
    sim_id = sim_id,
    vpc = sim_obj$vpc,
    vpc_obs = sim_obj$vpc_obs
  )
}

.shiny_load_run_vpc_state <- function(state, project, version_id, est_run_id, ws_root) {
  vpc_info <- .shiny_load_run_vpc(project, version_id, est_run_id, ws_root)
  if (is.null(vpc_info)) {
    state$sim_vpc_data <- NULL
    state$sim_vpc_obs <- NULL
    state$vpc_sim_id <- NULL
  } else {
    state$sim_vpc_data <- vpc_info$vpc
    state$sim_vpc_obs <- vpc_info$vpc_obs
    state$vpc_sim_id <- vpc_info$sim_id
  }
  invisible(state)
}

.shiny_sync_center_tabs <- function(session, fit, vpc_data) {
  if (.shiny_has_successful_fit(fit)) {
    shiny::showTab("center_tabs", "GOF", session = session)
  } else {
    shiny::hideTab("center_tabs", "GOF", session = session)
  }
  has_vpc <- !is.null(vpc_data) && is.data.frame(vpc_data) && nrow(vpc_data) > 0L
  if (has_vpc) {
    shiny::showTab("center_tabs", "VPC", session = session)
  } else {
    shiny::hideTab("center_tabs", "VPC", session = session)
  }
  has_npc <- .shiny_fit_has_npc(fit)
  has_npde <- .shiny_fit_has_npde(fit)
  if (has_npc) {
    shiny::showTab("center_tabs", "NPC", session = session)
  } else {
    shiny::hideTab("center_tabs", "NPC", session = session)
  }
  if (has_npde) {
    shiny::showTab("center_tabs", "NPDE", session = session)
  } else {
    shiny::hideTab("center_tabs", "NPDE", session = session)
  }
}

.shiny_ctl_model_help_modal <- function(advan, trans = NULL, focus = c("advan", "trans")) {
  focus <- match.arg(focus)
  info <- nm_ctl_model_info(advan, trans)
  title <- if (focus == "trans" && nm_ctl_show_trans(advan)) {
    paste0("TRANS ", trans, " — ", info$title)
  } else {
    info$title
  }
  modalDialog(
    title = title,
    easyClose = TRUE,
    size = "l",
    .shiny_ctl_model_help_body(advan, trans),
    footer = modalButton("Close")
  )
}

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body { padding-top: 0; }
      .app-header {
        background: #2b579a;
        color: #fff;
        padding: 8px 16px;
        display: flex;
        align-items: center;
        justify-content: space-between;
      }
      .app-header h3 { margin: 0; font-size: 18px; font-weight: 600; }
      .app-header .ws-path { font-size: 11px; opacity: 0.85; max-width: 50%; text-align: right; }
      .ribbon-wrap {
        background: #f3f2f1;
        border-bottom: 1px solid #d1d1d1;
        padding: 0 8px;
      }
      .ribbon-wrap .nav-tabs { border-bottom: none; margin-bottom: 0; }
      .ribbon-wrap .nav-tabs > li > a {
        border-radius: 0;
        padding: 10px 20px;
        color: #333;
        border: none;
        border-bottom: 3px solid transparent;
        margin-right: 2px;
      }
      .ribbon-wrap .nav-tabs > li.active > a,
      .ribbon-wrap .nav-tabs > li.active > a:focus,
      .ribbon-wrap .nav-tabs > li.active > a:hover {
        background: #fff;
        border: none;
        border-bottom: 3px solid #2b579a;
        color: #2b579a;
        font-weight: 600;
      }
      .ribbon-page {
        padding: 12px 16px;
      }
      .ribbon-page .data-tab-controls .selectize-control,
      .ribbon-page .data-tab-controls .selectize-input,
      .ribbon-page .data-tab-controls .selectize-dropdown,
      .ribbon-page .data-tab-controls select,
      .ribbon-page .data-tab-controls label,
      .ribbon-page .data-tab-controls .checkbox,
      .ribbon-page .data-tab-controls .radio {
        font-size: 11px;
      }
      .ribbon-page .data-tab-controls .form-group {
        margin-bottom: 8px;
      }
      .ribbon-page .data-tab-controls .shiny-input-container {
        margin-bottom: 8px;
      }
      .ribbon-page .data-tab-table-wrap {
        max-height: 320px;
        overflow: auto;
        margin-top: 6px;
      }
      .ribbon-page .data-tab-sidebar {
        max-height: calc(100vh - 180px);
        overflow-y: auto;
        overflow-x: hidden;
        padding-right: 4px;
      }
      .data-dataset-section {
        margin: 4px 0 12px;
        padding-bottom: 10px;
        border-bottom: 1px solid #e8e8e8;
      }
      .data-show-dataset-row {
        display: block;
      }
      .data-show-dataset-row .checkbox {
        margin-top: 8px;
        margin-bottom: 4px;
      }
      .data-show-dataset-row .checkbox label {
        font-weight: 600;
        color: #333;
      }
      .data-plot-section-title {
        font-size: 12px;
        font-weight: 600;
        color: #444;
        margin: 4px 0 8px;
      }
      .ribbon-panel {
        background: #fff;
        border-bottom: 1px solid #e1e1e1;
        padding: 12px 16px;
        min-height: 52px;
      }
      .main-workspace {
        padding: 12px 8px;
        height: calc(100vh - 168px);
        min-height: 320px;
        overflow: hidden;
      }
      .main-workspace > .row {
        height: 100%;
        margin-left: -8px;
        margin-right: -8px;
      }
      .main-workspace > .row > [class*='col-'] {
        height: 100%;
        min-height: 0;
        padding-left: 8px;
        padding-right: 8px;
      }
      .center-panel-wrap { height: 100%; min-height: 0; }
      .app-log-wrap {
        border-bottom: 1px solid #d1d1d1;
        background: #fafafa;
        font-size: 12px;
      }
      .app-log-banner {
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: 12px;
        padding: 6px 12px;
        min-height: 36px;
      }
      .app-log-banner.ok { background: #f5f5f5; color: #666; border-left: 4px solid #bdbdbd; }
      .app-log-banner.error { background: #fdecea; border-left: 4px solid #c62828; }
      .app-log-banner.warning { background: #fff8e1; border-left: 4px solid #f9a825; }
      .app-log-banner.message { background: #e8f4fd; border-left: 4px solid #1565c0; }
      .app-log-banner.info { background: #f5f5f5; border-left: 4px solid #757575; }
      .app-log-banner-main {
        flex: 1;
        min-width: 0;
        display: flex;
        align-items: baseline;
        gap: 8px;
        overflow: hidden;
      }
      .app-log-type {
        font-weight: 700;
        font-size: 11px;
        letter-spacing: 0.04em;
        flex-shrink: 0;
      }
      .app-log-time { color: #666; font-size: 11px; flex-shrink: 0; }
      .app-log-source { color: #888; font-size: 11px; flex-shrink: 0; }
      .app-log-msg {
        white-space: nowrap;
        overflow: hidden;
        text-overflow: ellipsis;
      }
      .app-log-actions {
        flex-shrink: 0;
        display: flex;
        align-items: center;
        gap: 10px;
        font-size: 11px;
      }
      .app-log-actions .checkbox { margin: 0; }
      .app-log-history {
        max-height: 160px;
        overflow-y: auto;
        border-top: 1px solid #e8e8e8;
        background: #fff;
      }
      .app-log-table {
        width: 100%;
        border-collapse: collapse;
        font-size: 11px;
        font-family: Consolas, Monaco, 'Courier New', monospace;
      }
      .app-log-table th,
      .app-log-table td {
        padding: 3px 8px;
        border-bottom: 1px solid #f0f0f0;
        text-align: left;
        vertical-align: top;
      }
      .app-log-table th { background: #fafafa; font-weight: 600; position: sticky; top: 0; }
      .log-row-error td { color: #c62828; }
      .log-row-warning td { color: #e65100; }
      .log-row-message td { color: #1565c0; }
      .panel-box {
        border: 1px solid #ddd;
        border-radius: 4px;
        background: #fff;
        height: 100%;
        min-height: 0;
        max-height: 100%;
        display: flex;
        flex-direction: column;
        overflow: hidden;
      }
      .panel-box-header {
        background: #fafafa;
        border-bottom: 1px solid #e8e8e8;
        padding: 8px 12px;
        font-weight: 600;
        font-size: 13px;
        flex-shrink: 0;
      }
      .panel-box-body {
        padding: 8px 10px;
        flex: 1;
        overflow: auto;
      }
      .panel-box-body-fill {
        flex: 1;
        overflow-y: auto;
        overflow-x: hidden;
        display: flex;
        flex-direction: column;
        min-height: 0;
      }
      .center-panel-tabs {
        flex: 1;
        display: flex;
        flex-direction: column;
        min-height: 0;
        overflow: hidden;
      }
      .center-panel-tabs > .tab-content {
        flex: 1;
        display: flex;
        flex-direction: column;
        min-height: 0;
        overflow: hidden;
      }
      .center-panel-tabs > .tab-content > .tab-pane.active {
        flex: 1;
        display: flex;
        flex-direction: column;
        min-height: 0;
        overflow: hidden;
      }
      .model-editor textarea {
        font-family: 'Consolas', 'Monaco', 'Courier New', monospace;
        font-size: 12px;
        line-height: 1.35;
        width: 100%;
        min-height: 100px;
        resize: vertical;
      }
      .model-editor-flex textarea { flex: 1; min-height: 80px; }
      .list-table .dataTables_wrapper { font-size: 12px; }
      .toolbar-row { margin-bottom: 6px; }
      .status-saved { color: #2e7d32; font-weight: 600; font-size: 12px; }
      .status-dirty { color: #b85c00; font-weight: 600; font-size: 12px; }
      .ctl-editor-scroll {
        flex: 1;
        overflow-y: auto;
        min-height: 0;
        padding-right: 4px;
        display: flex;
        flex-direction: column;
      }
      .ctl-section { margin-bottom: 10px; flex-shrink: 0; }
      .ctl-section-flex { flex: 1; min-height: 120px; display: flex; flex-direction: column; }
      .ctl-params-row {
        display: flex;
        gap: 6px;
        align-items: flex-start;
      }
      .ctl-params-row > .ctl-param-col {
        flex: 1;
        min-width: 0;
      }
      .ctl-params-row .dataTables_wrapper { font-size: 11px; }
      .ctl-code-row {
        display: flex;
        gap: 8px;
        flex: 1;
        min-height: 140px;
      }
      .ctl-code-row > .ctl-code-col {
        flex: 1;
        min-width: 0;
        display: flex;
        flex-direction: column;
      }
      .ctl-columns-table { max-height: 160px; overflow-y: auto; }
      .ctl-col-picker-table {
        width: 100%; font-size: 12px; border-collapse: collapse;
      }
      .ctl-col-picker-table th,
      .ctl-col-picker-table td {
        text-align: left; padding: 3px 8px; border-bottom: 1px solid #eee;
      }
      .ctl-col-picker-table th:nth-child(2),
      .ctl-col-picker-table th:nth-child(3),
      .ctl-col-picker-table td:nth-child(2),
      .ctl-col-picker-table td:nth-child(3) {
        text-align: center; width: 72px;
      }
      .ctl-col-picker-table .req-tag {
        color: #888; font-size: 10px; margin-left: 4px;
      }
      .ctl-col-picker-table input[type=checkbox] {
        width: 16px; height: 16px; cursor: pointer;
      }
      .ctl-dt-left .dt-left { text-align: left !important; }
      .ctl-section h5 {
        font-size: 12px; font-weight: 600; color: #444;
        margin: 0 0 4px; text-transform: uppercase; letter-spacing: 0.03em;
      }
      .right-panel-tabs {
        flex: 1;
        display: flex;
        flex-direction: column;
        min-height: 0;
        overflow: hidden;
        margin-top: 6px;
      }
      .right-panel-tabs > .tab-content {
        flex: 1;
        overflow: auto;
        min-height: 0;
      }
      .model-version-title { font-weight: 600; font-size: 13px; }
      .model-version-id { color: #666; font-size: 11px; margin-left: 6px; font-weight: 400; }
      .sidebar-actions { margin: 4px 0 8px; }
      .sidebar-actions .btn { margin-bottom: 4px; }
      .sidebar-actions-divider {
        border: none;
        border-top: 1px solid #d1d1d1;
        margin: 8px 0;
      }
      .ctl-inline-label {
        font-size: 12px;
        font-weight: 600;
        color: #444;
        white-space: nowrap;
        flex-shrink: 0;
        line-height: 28px;
      }
      .ctl-inline-row {
        display: flex;
        align-items: center;
        gap: 6px;
        flex-wrap: wrap;
      }
      .ctl-advan-trans-row {
        display: flex;
        align-items: center;
        flex-wrap: nowrap;
        gap: 6px;
      }
      .ctl-advan-trans-row .ctl-advan-select .shiny-input-container,
      .ctl-advan-trans-row .ctl-trans-select .shiny-input-container {
        margin-bottom: 0;
      }
      .ctl-advan-trans-row .shiny-conditional-panel {
        display: inline-flex;
        align-items: center;
        gap: 6px;
        margin: 0;
        padding: 0;
        vertical-align: middle;
      }
      .ctl-advan-trans-row .ctl-advan-select,
      .ctl-advan-trans-row .ctl-trans-select {
        width: 72px;
        flex-shrink: 0;
      }
      .ctl-advan-trans-row .ctl-advan-select .selectize-control,
      .ctl-advan-trans-row .ctl-trans-select .selectize-control {
        width: 72px !important;
      }
      .ctl-advan-trans-row .selectize-input {
        min-height: 28px;
        padding: 4px 8px;
      }
      .ctl-advan-trans-row .ctl-help-btn {
        margin: 0;
        padding: 0 4px;
        line-height: 1;
        align-self: center;
      }
      .ctl-advan-trans-row .ctl-ode-wrap {
        display: none;
        align-items: center;
        gap: 6px;
      }
      .ctl-advan-trans-row .ctl-ode-wrap .shiny-input-container {
        margin-bottom: 0;
        width: 80px;
      }
      .ctl-problem-row {
        display: flex;
        align-items: center;
        gap: 8px;
        margin-bottom: 8px;
        min-height: 28px;
      }
      .ctl-problem-text {
        flex: 1;
        font-size: 13px;
        color: #333;
        line-height: 1.35;
        min-width: 0;
      }
      .ctl-dataset-row .ctl-dataset-select { flex: 1; min-width: 120px; }
      .ctl-dataset-row .btn { flex-shrink: 0; }
      .right-param-table { font-size: 10px; }
      .right-param-table th,
      .right-param-table td { padding: 2px 5px; }
      .right-panel-tabs .dataTables_wrapper { font-size: 10px; }
      .right-panel-tabs table.dataTable { font-size: 10px; }
      .ctl-select-wrap { flex: 1; min-width: 0; }
      .ctl-help-btn {
        flex-shrink: 0;
        padding: 2px 6px;
        margin-bottom: 10px;
        color: #2b579a;
        font-size: 16px;
        line-height: 1;
      }
      .ctl-help-btn:hover { color: #1a3d6d; }
      .pk-diagram {
        display: flex;
        flex-wrap: wrap;
        align-items: center;
        gap: 8px;
        margin: 12px 0;
      }
      .pk-diagram-hub {
        flex-direction: column;
        align-items: stretch;
        gap: 10px;
      }
      .pk-depot-row, .pk-hub-row {
        display: flex;
        flex-wrap: wrap;
        align-items: center;
        gap: 8px;
      }
      .pk-hub-row { align-items: flex-start; }
      .pk-central-col {
        display: flex;
        flex-direction: column;
        align-items: center;
        gap: 4px;
      }
      .pk-peripheral-col {
        display: flex;
        flex-wrap: wrap;
        gap: 12px;
        align-items: flex-start;
      }
      .pk-branch {
        display: flex;
        flex-direction: column;
        align-items: center;
        gap: 4px;
      }
      .pk-box {
        border: 1px solid #2b579a;
        border-radius: 6px;
        padding: 8px 10px;
        background: #f7f9fc;
        font-size: 12px;
        text-align: center;
        min-width: 90px;
      }
      .pk-arrow { font-size: 13px; color: #555; white-space: nowrap; }
      .pk-arrow-v { font-size: 12px; color: #555; white-space: nowrap; }
      .pk-elim { font-size: 12px; color: #888; white-space: nowrap; }
      .pk-param-table td:first-child { width: 90px; }
      .run-vpc-badge, .run-flag-badge {
        margin-left: 6px;
        padding: 1px 5px;
        border-radius: 3px;
        font-size: 10px;
        font-weight: 600;
      }
      .run-vpc-badge {
        background: #e8f0fa;
        color: #2b579a;
      }
      .run-flag-bootstrap {
        background: #eef6ee;
        color: #2d6a2d;
      }
      .run-flag-npc {
        background: #f3eef8;
        color: #5a3d7a;
      }
      .run-flag-npde {
        background: #fef3e8;
        color: #9a5b2b;
      }
      .btn-action-secondary {
        background-color: #eef2f7;
        border-color: #c5cdd8;
        color: #2b3a4a;
      }
      .btn-action-secondary:hover,
      .btn-action-secondary:focus {
        background-color: #dfe6ef;
        border-color: #aeb9c8;
        color: #1a2633;
      }
      .btn-run-estimation {
        background-color: #2b579a;
        border-color: #234a82;
        color: #fff;
      }
      .btn-run-estimation:hover,
      .btn-run-estimation:focus {
        background-color: #234a82;
        border-color: #1c3d6b;
        color: #fff;
      }
      .btn-copy-version {
        background-color: #eef2f7;
        border-color: #c5cdd8;
        color: #2b3a4a;
      }
      .btn-copy-version:hover,
      .btn-copy-version:focus {
        background-color: #dfe6ef;
        border-color: #aeb9c8;
        color: #1a2633;
      }
      .version-tree {
        font-size: 12px;
        line-height: 1.35;
      }
      .version-row, .sim-row, .project-row, .job-row, .run-row {
        display: flex;
        align-items: baseline;
        gap: 6px;
        padding: 5px 6px;
        border-bottom: 1px solid #eee;
        cursor: pointer;
      }
      .version-row:hover, .sim-row:hover, .project-row:hover, .job-row:hover, .run-row:hover {
        background: #f5f8fb;
      }
      .version-row.selected, .sim-row.selected, .project-row.selected, .job-row.selected, .run-row.selected {
        background: #e8f0fa;
        font-weight: 600;
      }
      .version-toggle {
        width: 14px;
        flex: 0 0 14px;
        color: #666;
        cursor: pointer;
        user-select: none;
      }
      .version-id { font-weight: 600; color: #333; }
      .version-label { color: #444; flex: 1; min-width: 0; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
      .version-meta { color: #888; font-size: 10px; white-space: nowrap; }
      .sim-list { margin-left: 18px; border-left: 2px solid #e6e6e6; }
      .sim-row { padding-left: 10px; font-size: 11px; }
      .run-row { padding-left: 10px; font-size: 11px; }
      .run-row, .version-row { display: flex; align-items: center; }
      .compare-runs-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
        gap: 12px;
      }
      .compare-param-table { font-size: 11px; width: 100%; }
      .run-row.run-empty { color: #999; font-style: italic; cursor: default; }
      .run-id { font-weight: 600; color: #555; min-width: 48px; }
      .run-label { color: #666; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
      .sim-row.sim-empty { color: #999; font-style: italic; cursor: default; }
      .sim-id { font-weight: 600; color: #555; min-width: 48px; }
      .sim-label { color: #666; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
      .project-tree { max-height: 120px; overflow: auto; margin-bottom: 6px; }
      .project-id { font-weight: 600; color: #333; }
      .project-label { color: #666; flex: 1; min-width: 0; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
      .job-tree { max-height: 420px; overflow: auto; margin-top: 6px; }
      .job-toggle, .job-toggle-spacer { width: 14px; flex: 0 0 14px; color: #666; user-select: none; }
      .job-toggle { cursor: pointer; }
      .job-id { font-weight: 600; color: #333; min-width: 72px; }
      .job-label { color: #444; flex: 1; min-width: 0; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
      .job-row-error { background-color: #fdecea !important; }
      .job-row-error:hover { background-color: #f8d7da !important; }
      .job-status-error { background: #f8d7da; color: #842029; }
      .job-row-running { background-color: #eef6ff; }
      .job-row-success { background-color: #f6fff8; }
      .job-row-cancelled { background-color: #fff8e6; }
      .job-status { font-size: 10px; text-transform: uppercase; padding: 1px 5px; border-radius: 3px; background: #eee; color: #555; }
      .job-status-success { background: #d4edda; color: #155724; }
      .job-status-error { background: #f8d7da; color: #721c24; }
      .job-status-running { background: #cce5ff; color: #004085; }
      .job-status-queued { background: #e2e3e5; color: #383d41; }
      .job-status-cancelled { background: #fff3cd; color: #856404; }
      .job-type { color: #666; font-size: 10px; white-space: nowrap; }
      .job-meta { color: #888; font-size: 10px; white-space: nowrap; max-width: 140px; overflow: hidden; text-overflow: ellipsis; }
      .job-time { color: #999; font-size: 10px; white-space: nowrap; margin-left: auto; }
      .job-detail-panel { margin: 0 0 6px 20px; padding: 6px 8px; border-left: 2px solid #dce6f2; background: #fafcff; }
      .job-detail-table { width: 100%; font-size: 11px; border-collapse: collapse; }
      .job-detail-table th, .job-detail-table td { padding: 3px 6px; border-bottom: 1px solid #eee; text-align: left; }
      .job-detail-table th { color: #666; font-weight: 600; }
      .btn-sim-green {
        background-color: #2b579a;
        border-color: #234a82;
        color: #fff;
      }
      .btn-sim-green:hover, .btn-sim-green:focus {
        background-color: #234a82;
        border-color: #1c3d6b;
        color: #fff;
      }
    ")),
    tags$script(HTML("
      $(document).on('change', 'input.ctl-col-pick', function() {
        Shiny.setInputValue('ctl_col_picker_event', {
          col: $(this).attr('data-col'),
          kind: $(this).attr('data-kind'),
          value: this.checked,
          ts: new Date().getTime()
        }, {priority: 'event'});
      });
      $(document).on('click', '.version-toggle', function(e) {
        e.stopPropagation();
        Shiny.setInputValue('version_tree_event', {
          action: 'toggle',
          version: $(this).data('version'),
          ts: new Date().getTime()
        }, {priority: 'event'});
      });
      $(document).on('click', '.version-row', function(e) {
        if ($(e.target).hasClass('version-toggle')) { return; }
        Shiny.setInputValue('version_tree_event', {
          action: 'select_version',
          version: $(this).data('version'),
          ts: new Date().getTime()
        }, {priority: 'event'});
      });
      $(document).on('click', '.sim-row:not(.sim-empty)', function(e) {
        Shiny.setInputValue('version_tree_event', {
          action: 'select_sim',
          version: $(this).data('version'),
          sim: $(this).data('sim'),
          ts: new Date().getTime()
        }, {priority: 'event'});
      });
      $(document).on('click', '.run-row:not(.run-empty)', function(e) {
        Shiny.setInputValue('version_tree_event', {
          action: 'select_run',
          version: $(this).data('version'),
          run: $(this).data('run'),
          ts: new Date().getTime()
        }, {priority: 'event'});
      });
      $(document).on('click', '.project-row', function(e) {
        Shiny.setInputValue('project_tree_event', {
          project: $(this).data('project'),
          ts: new Date().getTime()
        }, {priority: 'event'});
      });
      $(document).on('click', '.job-toggle', function(e) {
        e.stopPropagation();
        Shiny.setInputValue('job_tree_event', {
          action: 'toggle',
          job: $(this).data('job'),
          ts: new Date().getTime()
        }, {priority: 'event'});
      });
      $(document).on('click', '.job-row', function(e) {
        if ($(e.target).hasClass('job-toggle')) { return; }
        Shiny.setInputValue('job_tree_event', {
          action: 'select',
          job: $(this).data('job'),
          ts: new Date().getTime()
        }, {priority: 'event'});
      });
      function syncAdvanTransRow() {
        var adv = $('#ctl_advan').val();
        var isOde = adv === '6';
        $('.ctl-advan-trans-row .ctl-trans-label, .ctl-advan-trans-row .ctl-trans-select').toggle(!isOde);
        $('.ctl-advan-trans-row .ctl-ode-wrap').css('display', isOde ? 'inline-flex' : 'none');
      }
      $(document).on('change', '#ctl_advan', syncAdvanTransRow);
      $(document).on('shiny:connected', function() {
        setTimeout(syncAdvanTransRow, 100);
      });
      $(document).on('shiny:value', function(event) {
        if (event.name === 'ctl_advan') {
          setTimeout(syncAdvanTransRow, 50);
        }
      });
    "))
  ),
  div(
    class = "app-header",
    tags$h3("LibeRation"),
    div(class = "ws-path", textOutput("workspace_path", inline = TRUE))
  ),
  div(
    class = "ribbon-wrap",
    tabsetPanel(
      id = "ribbon_tab",
      type = "tabs",
      tabPanel("Home", value = "home"),
      tabPanel("Jobs", value = "jobs"),
      tabPanel("Data", value = "data")
    )
  ),
  div(
    class = "app-log-wrap",
    uiOutput("app_log_banner"),
    conditionalPanel(
      condition = "input.log_history_open == true",
      uiOutput("app_log_history")
    ),
    div(
      class = "app-log-actions",
      style = "padding: 0 12px 6px; justify-content: flex-end;",
      checkboxInput("log_history_open", "Show history", value = FALSE),
      actionButton("clear_log", "Clear log", class = "btn btn-xs btn-default")
    )
  ),
  conditionalPanel(
    condition = "input.ribbon_tab == 'home' || !input.ribbon_tab",
    div(
      class = "main-workspace",
      fluidRow(
        column(
          2L,
          div(
            class = "panel-box",
            div(
              class = "panel-box-body panel-box-body-fill",
              uiOutput("projects_tree"),
              div(class = "sidebar-actions",
                  actionButton("create_project", "+ New project", class = "btn-primary btn-sm btn-block"),
                  actionButton("delete_project", "Delete project", class = "btn-danger btn-sm btn-block")
              ),
              tags$hr(style = "margin: 8px 0;"),
              div(class = "panel-box-header", style = "padding: 4px 0 6px; border: none; background: transparent;",
                  "Model versions"),
              div(class = "list-table", style = "flex: 1; min-height: 100px; overflow: auto;",
                  uiOutput("versions_tree")),
              uiOutput("compare_runs_panel"),
              tags$hr(class = "sidebar-actions-divider"),
              div(class = "sidebar-actions",
                  actionButton("new_version_copy", "Copy to new",
                               class = "btn-primary btn-sm btn-block"),
                  actionButton("new_version_template", "Create from template",
                               class = "btn-primary btn-sm btn-block"),
                  actionButton("run_estimation", "Run estimation",
                               class = "btn-primary btn-sm btn-block"),
                  actionButton("delete_version", "Delete version", class = "btn-danger btn-sm btn-block")
              ),
              tags$hr(class = "sidebar-actions-divider"),
              div(class = "sidebar-actions",
                  actionButton("create_simulation", "Create simulation",
                               class = "btn-sim-green btn-sm btn-block"),
                  actionButton("delete_simulation", "Delete simulation",
                               class = "btn-danger btn-sm btn-block")
              )
            )
          )
        ),
        column(
          7L,
          div(
            class = "center-panel-wrap",
            div(
            class = "panel-box",
            div(
              class = "panel-box-header",
              fluidRow(
                column(5L, uiOutput("model_version_header")),
                column(
                  7L,
                  style = "text-align: right;",
                  actionButton("reload_model", "Reload", class = "btn-xs btn-default"),
                  actionButton("save_model", "Overwrite existing", class = "btn-xs btn-primary"),
                  actionButton("save_model_as_new", "Save as new", class = "btn-xs btn-default")
                )
              )
            ),
            div(
              class = "panel-box-body panel-box-body-fill",
              uiOutput("model_dirty_banner"),
              div(
                class = "center-panel-tabs",
                tabsetPanel(
                id = "center_tabs",
                type = "tabs",
                tabPanel(
                  "Code",
                  div(
                    class = "ctl-editor-scroll",
                    tags$div(style = "display: none;",
                             textInput("ctl_problem", NULL, value = "", width = "100%")),
                    div(
                      class = "ctl-section",
                      uiOutput("ctl_problem_display")
                    ),
                    div(
                      class = "ctl-section ctl-advan-trans-row",
                      tags$span(class = "ctl-inline-label", "ADVAN"),
                      tags$div(
                        class = "ctl-advan-select",
                        selectInput(
                          "ctl_advan", NULL,
                          choices = .shiny_numeric_advan_choices(),
                          selected = "4",
                          width = "100%"
                        )
                      ),
                      tags$span(class = "ctl-inline-label ctl-trans-label", "TRANS"),
                      tags$div(
                        class = "ctl-trans-select",
                        selectInput(
                          "ctl_trans", NULL,
                          choices = .shiny_numeric_trans_choices(4L),
                          selected = nm_ctl_default_trans(4L),
                          width = "100%"
                        )
                      ),
                      actionButton(
                        "ctl_trans_help", NULL, icon = icon("question-circle"),
                        class = "btn-link ctl-help-btn ctl-trans-label",
                        title = "TRANS parameterization information"
                      ),
                      tags$div(
                        class = "ctl-ode-wrap",
                        tags$span(class = "ctl-inline-label", "Compartments"),
                        numericInput(
                          "ctl_ode_ncomp", NULL,
                          value = 2L, min = 1L, max = 10L, step = 1L,
                          width = "100%"
                        )
                      )
                    ),
                    div(
                      class = "ctl-section ctl-inline-row ctl-dataset-row",
                      tags$span(class = "ctl-inline-label", "Dataset ($DATA)"),
                      tags$div(class = "ctl-dataset-select", uiOutput("ctl_dataset_select")),
                      actionButton(
                        "open_column_picker", "Columns…",
                        class = "btn-xs btn-default",
                        icon = icon("columns"),
                        title = "Select $INPUT / $OUTPUT columns"
                      )
                    ),
                    div(
                      class = "ctl-section ctl-code-row",
                      div(
                        class = "ctl-code-col model-editor model-editor-flex",
                        tags$h5("$PK / $PRED"),
                        textAreaInput("ctl_pk", NULL, value = "", rows = 6L, width = "100%")
                      ),
                      conditionalPanel(
                        condition = "input.ctl_advan == '6' || input.ctl_advan == '13'",
                        div(
                          class = "ctl-code-col model-editor model-editor-flex",
                          tags$h5("$DES"),
                          textAreaInput("ctl_des", NULL, value = "", rows = 5L, width = "100%")
                        )
                      ),
                      div(
                        class = "ctl-code-col model-editor model-editor-flex",
                        tags$h5("$ERROR"),
                        textAreaInput("ctl_error", NULL, value = "Y = F", rows = 6L, width = "100%")
                      )
                    ),
                    div(
                      class = "ctl-section ctl-params-row",
                      div(
                        class = "ctl-param-col",
                        tags$h5("THETA"),
                        DT::dataTableOutput("ctl_theta_table")
                      ),
                      div(
                        class = "ctl-param-col",
                        tags$h5("OMEGA"),
                        DT::dataTableOutput("ctl_omega_table")
                      ),
                      div(
                        class = "ctl-param-col",
                        tags$h5("SIGMA"),
                        DT::dataTableOutput("ctl_sigma_table")
                      )
                    )
                  )
                ),
                tabPanel(
                  "GOF",
                  br(),
                  uiOutput("fit_summary_compact"),
                  uiOutput("gof_diag_status"),
                  fluidRow(
                    column(6L, plotOutput("gof_time_combined", height = "220px")),
                    column(6L, plotOutput("gof_dv_ipred", height = "220px"))
                  ),
                  fluidRow(
                    column(4L, plotOutput("gof_wres_time", height = "200px")),
                    column(4L, plotOutput("gof_cwres_time", height = "200px")),
                    column(4L, plotOutput("gof_qq_cwres", height = "200px"))
                  )
                ),
                tabPanel(
                  "NPC",
                  br(),
                  uiOutput("npc_tab_summary"),
                  fluidRow(
                    column(6L, plotOutput("gof_npc_hist", height = "240px")),
                    column(6L, plotOutput("gof_npc_time", height = "240px"))
                  )
                ),
                tabPanel(
                  "NPDE",
                  br(),
                  uiOutput("npde_tab_summary"),
                  fluidRow(
                    column(6L, plotOutput("gof_npde_qq", height = "240px")),
                    column(6L, plotOutput("gof_npde_time", height = "240px"))
                  )
                ),
                tabPanel(
                  "VPC",
                  br(),
                  uiOutput("vpc_tab_summary"),
                  checkboxInput(
                    "vpc_pc_correct",
                    "Prediction-corrected VPC (DV × PRED / IPRED)",
                    value = FALSE
                  ),
                  plotOutput("center_vpc_plot", height = "520px")
                )
              )
              )
            )
          )
          )
        ),
        column(
          3L,
          div(
            class = "panel-box",
            div(class = "panel-box-header", "Results"),
            div(
              class = "panel-box-body panel-box-body-fill",
              div(
                class = "right-panel-tabs",
                tabsetPanel(
                id = "right_tabs",
                type = "tabs",
                tabPanel(
                  "Parameters",
                  br(),
                  uiOutput("param_summary"),
                  DT::dataTableOutput("param_table")
                ),
                tabPanel(
                  "Report",
                  br(),
                  p(style = "font-size: 12px;",
                    "PDF report with JSON manifest for AI-assisted interpretation."),
                  checkboxGroupInput(
                    "report_sections",
                    "Sections",
                    choices = c(
                      "Fit summary" = "summary",
                      "Parameters" = "parameters",
                      "DV vs time (+ IPRED)" = "gof_time",
                      "IPRED/PRED vs time" = "gof_ipred_time",
                      "Scatter (DV vs IPRED/PRED)" = "gof_scatter",
                      "Residual diagnostics" = "gof_residuals",
                      "Shrinkage table" = "diag_shrinkage",
                      "ETA distributions" = "diag_eta",
                      "Interpretation stub" = "narrative_stub"
                    ),
                    selected = c("summary", "parameters", "gof_time", "gof_ipred_time",
                                 "gof_scatter", "gof_residuals", "diag_shrinkage",
                                 "diag_eta", "narrative_stub")
                  ),
                  textInput("report_name", "Filename (no extension)",
                            value = paste0("report_", format(Sys.Date(), "%Y%m%d")),
                            width = "100%"),
                  actionButton("generate_report", "Generate PDF", class = "btn-primary btn-sm"),
                  tags$hr(),
                  uiOutput("report_status")
                )
              )
              )
            )
          )
        )
      )
    )
  ),
  conditionalPanel(
    condition = "input.ribbon_tab == 'jobs'",
    div(
      class = "ribbon-page",
      fluidRow(
        column(
          8L,
          actionButton("refresh_jobs", "Refresh", class = "btn-sm"),
          actionButton("cancel_job", "Cancel selected", class = "btn-warning btn-sm"),
          actionButton("cleanup_jobs", "Clear finished", class = "btn-sm")
        ),
        column(4L, textOutput("jobs_refresh_clock", inline = TRUE))
      ),
      uiOutput("job_status_banner"),
      br(),
      uiOutput("jobs_tree"),
      tags$hr(),
      h4("Worker log"),
      verbatimTextOutput("job_log"),
      tags$p(class = "text-muted", style = "font-size: 11px;", uiOutput("job_root_display"))
    )
  ),
  conditionalPanel(
    condition = "input.ribbon_tab == 'data'",
    div(
      class = "ribbon-page",
      fluidRow(
        column(
          4L,
          div(
            class = "data-tab-sidebar data-tab-controls",
            uiOutput("data_summary"),
            div(
              class = "data-dataset-section",
              uiOutput("data_dataset_picker"),
              div(
                class = "data-show-dataset-row",
                checkboxInput("data_show_dataset", "Show dataset table", value = FALSE)
              ),
              conditionalPanel(
                condition = "input.data_show_dataset == true",
                checkboxInput("data_explore_all_rows", "Show all rows (incl. doses)", value = FALSE),
                div(
                  class = "data-tab-table-wrap list-table",
                  DT::dataTableOutput("data_table")
                )
              )
            ),
            tags$div(class = "data-plot-section-title", "Plot"),
            selectInput("data_explore_x", "X axis", choices = character(), width = "100%"),
            selectInput("data_explore_y", "Y axis", choices = character(), width = "100%"),
            selectInput(
              "data_explore_strat", "Stratify / colour by (optional)",
              choices = c("(none)" = ""), width = "100%"
            ),
            selectInput(
              "data_explore_split", "Split by (facets, optional)",
              choices = c("(none)" = ""), width = "100%"
            ),
            checkboxInput("data_explore_bin_x", "Bin X (continuous)", value = FALSE),
            conditionalPanel(
              condition = "input.data_explore_bin_x",
              radioButtons(
                "data_explore_bin_x_pos", "X bin position",
                choices = c("Equidistant bins" = "equal", "At bin midpoints" = "midpoint"),
                selected = "equal", inline = TRUE
              ),
              checkboxInput("data_explore_bin_x_manual", "Manual X breaks", value = FALSE),
              conditionalPanel(
                condition = "input.data_explore_bin_x && input.data_explore_bin_x_manual",
                textInput(
                  "data_explore_bin_x_breaks", "X break values",
                  value = "", placeholder = "e.g. 0, 2, 5, 10",
                  width = "100%"
                )
              )
            ),
            checkboxInput("data_explore_bin_y", "Bin Y (continuous)", value = FALSE),
            conditionalPanel(
              condition = "input.data_explore_bin_y",
              radioButtons(
                "data_explore_bin_y_pos", "Y bin position",
                choices = c("Equidistant bins" = "equal", "At bin midpoints" = "midpoint"),
                selected = "equal", inline = TRUE
              ),
              checkboxInput("data_explore_bin_y_manual", "Manual Y breaks", value = FALSE),
              conditionalPanel(
                condition = "input.data_explore_bin_y && input.data_explore_bin_y_manual",
                textInput(
                  "data_explore_bin_y_breaks", "Y break values",
                  value = "", placeholder = "e.g. 0, 1, 5, 20",
                  width = "100%"
                )
              )
            ),
            conditionalPanel(
              condition = paste0(
                "(input.data_explore_bin_x && !input.data_explore_bin_x_manual) || ",
                "(input.data_explore_bin_y && !input.data_explore_bin_y_manual) || ",
                "input.data_explore_plot_type == 'mean_se' || ",
                "input.data_explore_plot_type == 'median_q' || ",
                "input.data_explore_plot_type == 'boxplot' || ",
                "input.data_explore_plot_type == 'violin'"
              ),
              numericInput(
                "data_explore_bins", "Number of bins",
                value = 10L, min = 3L, max = 50L, step = 1L, width = "100%"
              )
            ),
            selectInput(
              "data_explore_plot_type", "Plot type",
              choices = c(
                "Points" = "points",
                "Jittered points" = "jitter",
                "Lines" = "lines",
                "Points + lines" = "both",
                "Smooth (loess)" = "smooth",
                "Linear regression" = "regression",
                "Box plot" = "boxplot",
                "Violin" = "violin",
                "Mean ± SE" = "mean_se",
                "Median + quantiles" = "median_q"
              ),
              selected = "points",
              width = "100%"
            ),
            conditionalPanel(
              condition = "input.data_explore_plot_type == 'mean_se' || input.data_explore_plot_type == 'median_q'",
              checkboxInput("data_explore_show_points", "Show individual points", value = FALSE)
            ),
            conditionalPanel(
              condition = "input.data_explore_plot_type == 'both' || input.data_explore_plot_type == 'lines'",
              radioButtons(
                "data_explore_line_mode", "Line mode",
                choices = c(
                  "Each individual" = "individual",
                  "Mean across subjects" = "mean",
                  "Median across subjects" = "median",
                  "No line" = "none"
                ),
                selected = "individual"
              )
            ),
            checkboxInput("data_explore_adjust_plot", "Adjust plot options", value = FALSE),
            conditionalPanel(
              condition = "input.data_explore_adjust_plot",
              textInput(
                "data_explore_title", "Plot title",
                value = "", placeholder = "Auto (Y vs X)", width = "100%"
              ),
              textInput(
                "data_explore_xlab", "X axis label",
                value = "", placeholder = "Auto", width = "100%"
              ),
              textInput(
                "data_explore_ylab", "Y axis label",
                value = "", placeholder = "Auto", width = "100%"
              ),
              conditionalPanel(
                condition = paste0(
                  "input.data_explore_plot_type == 'points' || ",
                  "input.data_explore_plot_type == 'jitter' || ",
                  "input.data_explore_plot_type == 'both' || ",
                  "(input.data_explore_plot_type == 'mean_se' && input.data_explore_show_points) || ",
                  "(input.data_explore_plot_type == 'median_q' && input.data_explore_show_points)"
                ),
                sliderInput(
                  "data_explore_point_size", "Point size",
                  min = 0.4, max = 1.8, value = 0.85, step = 0.05, width = "100%"
                ),
                selectInput(
                  "data_explore_point_shape", "Point shape",
                  choices = c(
                    "Circle" = 16,
                    "Open circle" = 1,
                    "Square" = 15,
                    "Triangle" = 17,
                    "Diamond" = 18
                  ),
                  selected = 16,
                  width = "100%"
                )
              ),
              conditionalPanel(
                condition = "input.data_explore_plot_type == 'median_q'",
                numericInput(
                  "data_explore_q_interval", "Quantile interval (%)",
                  value = 95, min = 50, max = 99.9, step = 0.5, width = "100%"
                ),
                sliderInput(
                  "data_explore_shade_alpha", "Quantile shade intensity",
                  min = 5, max = 70, value = 25, step = 5, width = "100%"
                ),
                tags$p(class = "text-muted", style = "font-size: 11px;",
                       "Solid line = median; dashed lines and shaded band = outer quantiles (default 2.5% and 97.5%).")
              ),
              conditionalPanel(
                condition = paste0(
                  "(input.data_explore_plot_type == 'mean_se' || ",
                  "input.data_explore_plot_type == 'median_q') && input.data_explore_show_points"
                ),
                sliderInput(
                  "data_explore_point_scatter", "Point scatter",
                  min = 0, max = 100, value = 25, step = 5, width = "100%"
                )
              )
            ),
            tags$p(
              class = "text-muted",
              style = "font-size: 11px;",
              "Observations only (MDV=0, EVID=0). Select a project dataset or load a model version."
            )
          )
        ),
        column(
          8L,
          plotOutput("data_explore_plot", height = "480px")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  ws_root <- nm_workspace_root()
  job_root <- nm_job_root()

  state <- reactiveValues(
    selected_project = NULL,
    selected_version_id = NULL,
    selected_est_run_id = NULL,
    est_run_loaded_version_id = NULL,
    selected_sim_id = NULL,
    expanded_versions = character(),
    expanded_jobs = character(),
    saved_text = "",
    dirty = FALSE,
    ctl_data_file = NULL,
    dataset_columns = character(),
    ctl_input_cols = character(),
    ctl_output_cols = character(),
    columns_rev = 0L,
    ctl_thetas = NULL,
    ctl_omegas = NULL,
    ctl_sigmas = NULL,
    model = NULL,
    data = NULL,
    fit = NULL,
    handles = list(),
    selected_job = NULL,
    job_status_cache = list(),
    job_watch_last = "",
    jobs_rev = 0L,
    initialized = FALSE,
    last_report = NULL,
    projects_rev = 0L,
    versions_rev = 0L,
    est_runs_rev = 0L,
    simulations_rev = 0L,
    active_job_version = NULL,
    active_job_est_run = NULL,
    active_job_project = NULL,
    active_job_type = NULL,
    active_job_sim_id = NULL,
    active_job_diag_only = FALSE,
    log_entries = list(),
    log_rev = 0L,
    explore_df = NULL,
    explore_data_rev = 0L,
    sim_vpc_data = NULL,
    sim_vpc_obs = NULL,
    vpc_sim_id = NULL,
    gof_diag_rev = 0L,
    auto_data_sim_key = NULL,
    auto_data_est_key = NULL,
    explore_all_rows = FALSE,
    compare_entries = NULL
  )

  flags <- new.env(parent = emptyenv())
  flags$loading_ctl <- FALSE
  flags$switching <- FALSE
  flags$pending_version_id <- NULL
  flags$pending_est_run_id <- NULL
  flags$pending_project <- NULL
  flags$switch_kind <- NULL
  flags$pending_sim <- NULL
  flags$artifact_warn_shown <- FALSE
  flags$ctl_load_gen <- 0L

  `%||%` <- function(x, y) if (is.null(x)) y else x

  push_log <- function(type, message, source = NULL, notify = NULL) {
    type <- match.arg(tolower(type), c("error", "warning", "message", "info"))
    msg <- trimws(paste(as.character(message), collapse = "\n"))
    if (!nzchar(msg)) {
      return(invisible(NULL))
    }
    entry <- list(
      time = format(Sys.time(), "%H:%M:%S"),
      type = type,
      message = msg,
      source = source %||% ""
    )
    isolate({
      entries <- state$log_entries
      entries <- c(list(entry), entries)
      if (length(entries) > .shiny_log_max_entries) {
        entries <- entries[seq_len(.shiny_log_max_entries)]
      }
      state$log_entries <- entries
      state$log_rev <- state$log_rev + 1L
    })
    if (is.null(notify)) {
      notify <- type %in% c("error", "warning", "message")
    }
    if (isTRUE(notify)) {
      toast <- if (nzchar(entry$source)) {
        paste0("[", entry$source, "] ", msg)
      } else {
        msg
      }
      shiny::showNotification(
        toast,
        type = if (type == "info") "message" else type,
        duration = if (type == "error") NULL else 15
      )
    }
    invisible(entry)
  }

  showNotification <- function(text, type = "default", action = NULL, duration = 15,
                               closeButton = TRUE, id = NULL, ...) {
    log_type <- tolower(type)
    if (log_type %in% c("default", "notification")) {
      log_type <- "message"
    }
    push_log(log_type, text, notify = FALSE)
    shiny::showNotification(
      text,
      type = type,
      action = action,
      duration = duration,
      closeButton = closeButton,
      id = id,
      ...
    )
  }

  safe_run <- function(label, expr) {
    withCallingHandlers(
      tryCatch(
        force(expr),
        error = function(e) {
          if (inherits(e, "shiny.silent.error") || inherits(e, "validation")) {
            stop(e)
          }
          push_log("error", conditionMessage(e), source = label)
          invisible(NULL)
        }
      ),
      warning = function(w) {
        push_log("warning", conditionMessage(w), source = label, notify = FALSE)
        invokeRestart("muffleWarning")
      }
    )
  }

  .dt_left_opts <- function(extra = list()) {
    modifyList(
      list(columnDefs = list(list(className = "dt-left", targets = "_all"))),
      extra
    )
  }

  .current_advan <- function() {
    as.integer(input$ctl_advan %||% "3")
  }

  .current_trans <- function() {
    adv <- .current_advan()
    if (!nm_ctl_show_trans(adv)) {
      return(nm_ctl_effective_trans(adv))
    }
    as.integer(input$ctl_trans %||% nm_ctl_default_trans(adv))
  }

  .sync_advan_input_cols <- function() {
    ds <- state$dataset_columns
    if (length(ds) == 0L) {
      return(invisible(NULL))
    }
    advan <- .current_advan()
    trans <- .current_trans()
    essential <- intersect(nm_ctl_essential_input_cols(advan, trans), ds)
    keep <- intersect(state$ctl_input_cols, ds)
    in_sel <- unique(c(essential, keep))
    in_sel <- ds[ds %in% in_sel]
    state$ctl_input_cols <- in_sel
    state$columns_rev <- state$columns_rev + 1L
    invisible(in_sel)
  }

  .dataset_full_path <- function(project, rel_path) {
    if (!.valid_id(project) || is.null(rel_path) || !nzchar(rel_path)) {
      return(NULL)
    }
    proj_dir <- file.path(ws_root, "projects", project)
    cand <- file.path(proj_dir, rel_path)
    if (file.exists(cand)) {
      return(normalizePath(cand, winslash = "/", mustWork = FALSE))
    }
    cand2 <- file.path(proj_dir, "data", basename(rel_path))
    if (file.exists(cand2)) {
      return(normalizePath(cand2, winslash = "/", mustWork = FALSE))
    }
    NULL
  }

  .refresh_input_cols <- function(project, rel_path) {
    fp <- .dataset_full_path(project, rel_path)
    if (is.null(fp)) {
      return(character())
    }
    nm_ctl_read_columns(fp)
  }

  .sync_input_selection <- function(dataset_cols, current_input = NULL, current_output = NULL,
                                    advan = 3L, trans = 4L) {
    essential <- intersect(nm_ctl_essential_input_cols(advan, trans), dataset_cols)
    if (length(current_input) > 0L) {
      input_sel <- unique(c(essential, intersect(current_input, dataset_cols)))
    } else {
      input_sel <- unique(c(essential, dataset_cols))
    }
    input_sel <- dataset_cols[dataset_cols %in% input_sel]
    output_sel <- if (length(current_output) > 0L) {
      current_output[current_output %in% nm_ctl_picker_columns(dataset_cols)]
    } else {
      character()
    }
    list(
      dataset_columns = dataset_cols,
      input_cols = input_sel,
      output_cols = output_sel
    )
  }

  .apply_column_selection <- function(dataset_cols, input_cols, output_cols) {
    state$dataset_columns <- dataset_cols
    state$ctl_input_cols <- input_cols
    state$ctl_output_cols <- output_cols
    state$columns_rev <- state$columns_rev + 1L
  }

  columns_picker_df <- reactive({
    state$columns_rev
    state$dataset_columns
    state$ctl_input_cols
    state$ctl_output_cols
    input$ctl_advan
    input$ctl_trans
    cols <- state$dataset_columns
    if (length(cols) == 0L) {
      return(data.frame(
        Column = character(), Required = character(),
        IN = logical(), OUT = logical(), in_dataset = logical(),
        stringsAsFactors = FALSE
      ))
    }
    all <- nm_ctl_picker_columns(cols)
    essential <- intersect(nm_ctl_essential_input_cols(.current_advan(), .current_trans()), all)
    data.frame(
      Column = all,
      Required = ifelse(all %in% essential, "req", ""),
      IN = all %in% state$ctl_input_cols | (all %in% essential & all %in% cols),
      OUT = all %in% state$ctl_output_cols,
      in_dataset = all %in% cols,
      stringsAsFactors = FALSE
    )
  })

  .update_cols_from_picker <- function(df) {
    dataset_cols <- state$dataset_columns
    essential <- intersect(
      nm_ctl_essential_input_cols(.current_advan(), .current_trans()),
      dataset_cols
    )
    in_sel <- df$Column[df$IN & df$in_dataset]
    in_sel <- unique(c(essential, in_sel))
    in_sel <- dataset_cols[dataset_cols %in% in_sel]
    out_sel <- df$Column[df$OUT]
    state$ctl_input_cols <- in_sel
    state$ctl_output_cols <- out_sel
    state$columns_rev <- state$columns_rev + 1L
  }

  collect_ctl_parts <- function() {
    advan <- .current_advan()
    trans <- nm_ctl_effective_trans(advan, .current_trans())
    if (!nm_ctl_is_valid_pair(advan, trans)) {
      trans <- as.integer(nm_ctl_effective_trans(advan, nm_ctl_default_trans(advan)))
    }
    list(
      problem = input$ctl_problem,
      advan = advan,
      trans = trans,
      use_ode = nm_ctl_use_ode(advan),
      subroutine = nm_ctl_subroutine_text(advan, trans),
      data_file = state$ctl_data_file,
      input_cols = state$ctl_input_cols,
      output_cols = state$ctl_output_cols,
      thetas = state$ctl_thetas,
      omegas = state$ctl_omegas,
      sigmas = state$ctl_sigmas,
      pk = input$ctl_pk,
      des = if (nm_ctl_use_ode(advan)) input$ctl_des %||% "" else "",
      error = input$ctl_error
    )
  }

  .shiny_refresh_ctl_des <- function(advan, ncomp = NULL) {
    if (!nm_ctl_use_ode(as.integer(advan))) {
      return(invisible(NULL))
    }
    ncomp <- as.integer(ncomp %||% input$ctl_ode_ncomp %||% 2L)
    ncomp <- max(1L, min(10L, ncomp))
    new_des <- .shiny_default_des(as.integer(advan), 1L, ncomp = ncomp, oral = advan != "13")
    cur <- input$ctl_des %||% ""
    if (!nzchar(trimws(cur)) || identical(trimws(cur), trimws(state$ctl_des_last_default %||% ""))) {
      updateTextAreaInput(session, "ctl_des", value = new_des)
      state$ctl_des_last_default <- new_des
    }
    invisible(new_des)
  }

  compose_ctl_baseline <- function(parts = NULL) {
    advan <- if (!is.null(parts)) as.integer(parts$advan) else .current_advan()
    trans <- if (!is.null(parts)) {
      nm_ctl_effective_trans(advan, parts$trans)
    } else {
      nm_ctl_effective_trans(advan, .current_trans())
    }
    if (!nm_ctl_is_valid_pair(advan, trans)) {
      trans <- as.integer(nm_ctl_effective_trans(advan, nm_ctl_default_trans(advan)))
    }
    des_val <- if (!is.null(parts)) parts$des %||% "" else input$ctl_des %||% ""
    if (nm_ctl_use_ode(advan) && !nzchar(trimws(des_val))) {
      ncomp <- 2L
      des_val <- .shiny_default_des(advan, 1L, ncomp = ncomp, oral = advan != 13L)
    }
    list(
      problem = if (!is.null(parts)) parts$problem %||% "" else input$ctl_problem,
      advan = advan,
      trans = trans,
      use_ode = nm_ctl_use_ode(advan),
      subroutine = if (!is.null(parts)) {
        parts$subroutine %||% nm_ctl_subroutine_text(advan, trans)
      } else {
        nm_ctl_subroutine_text(advan, trans)
      },
      data_file = state$ctl_data_file,
      input_cols = state$ctl_input_cols,
      output_cols = state$ctl_output_cols,
      thetas = state$ctl_thetas,
      omegas = state$ctl_omegas,
      sigmas = state$ctl_sigmas,
      pk = if (!is.null(parts)) parts$pk %||% "" else input$ctl_pk,
      des = if (nm_ctl_use_ode(advan)) des_val else "",
      error = if (!is.null(parts)) parts$error %||% "Y = F" else input$ctl_error
    )
  }

  sync_ctl_baseline <- function(parts = NULL) {
    state$saved_text <- nm_ctl_canonical(nm_ctl_compose(compose_ctl_baseline(parts)))
    state$dirty <- FALSE
  }

  sync_dirty <- function() {
    if (isTRUE(flags$loading_ctl) || !.valid_id(state$selected_version_id)) {
      return()
    }
    cur <- nm_ctl_canonical(nm_ctl_compose(collect_ctl_parts()))
    saved <- nm_ctl_canonical(state$saved_text %||% "")
    was_dirty <- isTRUE(state$dirty)
    state$dirty <- !identical(cur, saved)
    if (!was_dirty && isTRUE(state$dirty) && !isTRUE(flags$artifact_warn_shown)) {
      flags$artifact_warn_shown <- TRUE
      .shiny_warn_version_artifacts(
        state$selected_project, state$selected_version_id, ws_root
      )
    }
  }

  push_ctl_to_ui <- function(parts, project, keep_loading = FALSE) {
    flags$loading_ctl <- TRUE
    updateTextInput(session, "ctl_problem", value = parts$problem %||% "")
    advan <- as.character(parts$advan)
    trans <- as.character(nm_ctl_effective_trans(parts$advan, parts$trans))
    trans_named <- .shiny_numeric_trans_choices(parts$advan)
    updateSelectInput(
      session, "ctl_advan",
      choices = .shiny_numeric_advan_choices(),
      selected = advan
    )
    if (!is.null(trans_named)) {
      updateSelectInput(
        session, "ctl_trans",
        choices = trans_named,
        selected = trans
      )
    }
    state$ctl_data_file <- parts$data_file
    if (!is.null(parts$data_file) && nzchar(parts$data_file)) {
      ds_cols <- .refresh_input_cols(project, parts$data_file)
      if (length(ds_cols) > 0L) {
        synced <- .sync_input_selection(
          ds_cols,
          current_input = parts$input_cols,
          current_output = parts$output_cols,
          advan = as.integer(advan),
          trans = as.integer(trans)
        )
        .apply_column_selection(
          synced$dataset_columns,
          synced$input_cols,
          synced$output_cols
        )
      } else {
        state$ctl_input_cols <- parts$input_cols
        state$ctl_output_cols <- parts$output_cols %||% character()
        state$dataset_columns <- parts$input_cols
      }
    } else {
      state$dataset_columns <- parts$input_cols
      state$ctl_input_cols <- parts$input_cols
      state$ctl_output_cols <- parts$output_cols %||% character()
    }
    th <- parts$thetas
    if (!is.null(th) && nrow(th) > 0L) {
      th <- LibeRation:::.nm_ctl_normalize_thetas(th)
      if (length(th$Label) != nrow(th)) {
        th$Label <- paste0("THETA", th$THETA)
      }
    }
    state$ctl_thetas <- th
    state$ctl_omegas <- .shiny_ctl_param_labels(parts$omegas, "OMEGA")
    state$ctl_sigmas <- .shiny_ctl_param_labels(parts$sigmas, "SIGMA")
    updateTextAreaInput(session, "ctl_pk", value = parts$pk %||% "")
    des_val <- parts$des %||% ""
    if (nm_ctl_use_ode(as.integer(advan)) && !nzchar(des_val)) {
      ncomp <- 2L
      des_val <- .shiny_default_des(as.integer(advan), 1L, ncomp = ncomp, oral = advan != "13")
    }
    state$ctl_des_last_default <- des_val
    updateTextAreaInput(session, "ctl_des", value = des_val)
    updateTextAreaInput(session, "ctl_error", value = parts$error %||% "Y = F")
    if (!isTRUE(keep_loading)) {
      flags$loading_ctl <- FALSE
    }
  }

  clear_ctl_ui <- function() {
    flags$loading_ctl <- TRUE
    updateTextInput(session, "ctl_problem", value = "")
    updateTextAreaInput(session, "ctl_pk", value = "")
    updateTextAreaInput(session, "ctl_des", value = "")
    updateTextAreaInput(session, "ctl_error", value = "Y = F")
    state$ctl_data_file <- NULL
    state$dataset_columns <- character()
    state$ctl_input_cols <- character()
    state$ctl_output_cols <- character()
    state$columns_rev <- state$columns_rev + 1L
    state$ctl_thetas <- NULL
    state$ctl_omegas <- NULL
    state$ctl_sigmas <- NULL
    flags$loading_ctl <- FALSE
  }

  datasets_df <- reactive({
    state$projects_rev
    input$confirm_create_project
    input$delete_project
    proj <- state$selected_project
    if (!.valid_id(proj)) {
      return(data.frame(file = character(), path = character(), stringsAsFactors = FALSE))
    }
    nm_workspace_list_datasets(proj, root = ws_root)
  })

  message("[LibeRation] Shiny server — workspace: ", ws_root)

  output$workspace_path <- renderText({
    ws_root
  })

  output$app_log_banner <- renderUI({
    state$log_rev
    entries <- state$log_entries
    n <- length(entries)
    latest <- if (n > 0L) entries[[1L]] else NULL
    banner_cls <- if (is.null(latest)) {
      "app-log-banner ok"
    } else {
      paste("app-log-banner", latest$type)
    }
    tags$div(
      class = banner_cls,
      tags$div(
        class = "app-log-banner-main",
        if (is.null(latest)) {
          tagList(
            tags$span(class = "app-log-type", "OK"),
            tags$span(class = "app-log-msg", "No messages")
          )
        } else {
          tagList(
            tags$span(class = "app-log-type", toupper(latest$type)),
            tags$span(class = "app-log-time", latest$time),
            if (nzchar(latest$source)) {
              tags$span(class = "app-log-source", paste0("[", latest$source, "]"))
            },
            tags$span(class = "app-log-msg", title = latest$message, latest$message)
          )
        }
      ),
      if (n > 0L) {
        tags$span(
          style = "color:#888;font-size:11px;white-space:nowrap;",
          sprintf("%d entr%s", n, if (n == 1L) "y" else "ies")
        )
      }
    )
  })

  output$app_log_history <- renderUI({
    state$log_rev
    entries <- state$log_entries
    if (length(entries) == 0L) {
      return(NULL)
    }
    tags$table(
      class = "app-log-table",
      tags$thead(
        tags$tr(
          tags$th("Time"),
          tags$th("Type"),
          tags$th("Source"),
          tags$th("Message")
        )
      ),
      tags$tbody(
        lapply(entries, function(e) {
          tags$tr(
            class = paste0("log-row-", e$type),
            tags$td(e$time),
            tags$td(e$type),
            tags$td(if (nzchar(e$source)) e$source else "—"),
            tags$td(e$message)
          )
        })
      )
    )
  })

  observeEvent(input$clear_log, {
    isolate({
      state$log_entries <- list()
      state$log_rev <- state$log_rev + 1L
    })
    push_log("info", "Log cleared.", source = "Log", notify = TRUE)
  })

  projects_df <- reactive({
    state$projects_rev
    input$confirm_create_project
    input$delete_project
    df <- data.frame(project = nm_workspace_list_projects(ws_root), stringsAsFactors = FALSE)
    if (nrow(df) == 0L) {
      df <- data.frame(project = character(), stringsAsFactors = FALSE)
    }
    df
  })

  output$job_root_display <- renderUI({
    tags$span(paste("Job directory:", job_root))
  })

  .version_meta <- function(project, version_id) {
    proj_dir <- file.path(ws_root, "projects", project)
    for (sub in c("versions", "models")) {
      meta_path <- file.path(proj_dir, sub, version_id, "meta.json")
      if (file.exists(meta_path)) {
        return(tryCatch(
          if (requireNamespace("jsonlite", quietly = TRUE)) {
            jsonlite::fromJSON(meta_path, simplifyVector = TRUE)
          } else {
            list()
          },
          error = function(e) list()
        ))
      }
    }
    list()
  }

  versions_df <- reactive({
    state$versions_rev
    state$simulations_rev
    input$new_version_copy
    input$confirm_template_version
    input$save_model
    proj <- state$selected_project
    if (!.valid_id(proj)) {
      return(data.frame(
        version = character(), label = character(), runs = integer(),
        stringsAsFactors = FALSE
      ))
    }
    vers <- nm_workspace_list_versions(proj, root = ws_root)
    if (length(vers) == 0L) {
      return(data.frame(
        version = character(), label = character(), runs = integer(),
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      version = vers,
      label = vapply(vers, function(v) {
        meta <- .version_meta(proj, v)
        .shiny_meta_label(meta$label, fallback = v)
      }, character(1L)),
      runs = vapply(vers, function(v) {
        nrow(nm_workspace_list_runs(proj, v, root = ws_root))
      }, integer(1L)),
      stringsAsFactors = FALSE
    )
  })

  est_runs_df <- reactive({
    state$est_runs_rev
    state$fit
    input$submit_job
    proj <- state$selected_project
    ver <- state$selected_version_id
    if (!.valid_id(proj) || !.valid_id(ver)) {
      return(data.frame(
        run_id = character(), method = character(), objective = character(),
        label = character(), stringsAsFactors = FALSE
      ))
    }
    df <- nm_workspace_list_runs(proj, ver, root = ws_root)
    if (nrow(df) == 0L) {
      return(data.frame(
        run_id = character(), method = character(), objective = character(),
        label = character(), stringsAsFactors = FALSE
      ))
    }
    df$objective <- ifelse(
      is.na(df$objective),
      "",
      format(round(df$objective, 4))
    )
    df[, c("run_id", "method", "objective", "label")]
  })

  output$projects_tree <- renderUI({
    state$projects_rev
    state$selected_project
    .shiny_projects_tree_ui(projects_df(), state$selected_project)
  })

  output$versions_tree <- renderUI({
    state$versions_rev
    state$simulations_rev
    state$est_runs_rev
    state$expanded_versions
    state$selected_version_id
    state$selected_sim_id
    state$selected_est_run_id
    proj <- state$selected_project
    if (!.valid_id(proj)) {
      return(tags$p(
        class = "text-muted",
        style = "font-size: 11px; padding: 6px;",
        "Select a project."
      ))
    }
    .shiny_versions_tree_ui(
      versions_df(),
      proj,
      ws_root,
      state$selected_version_id,
      state$selected_sim_id,
      state$selected_est_run_id,
      state$expanded_versions
    )
  })

  output$compare_runs_panel <- renderUI({
    tags$div(
      class = "sidebar-actions",
      style = "margin-top: 0; padding-top: 0;",
      actionButton(
        "open_model_comparison",
        "Compare runs",
        class = "btn-primary btn-sm btn-block"
      ),
      actionButton(
        "delete_est_run",
        "Delete run",
        class = "btn-danger btn-sm btn-block"
      )
    )
  })

  observeEvent(input$open_model_comparison, {
    proj <- state$selected_project
    if (!.valid_id(proj)) {
      showNotification("Select a project first.", type = "warning")
      return()
    }
    choices <- .shiny_project_run_choices(proj, ws_root)
    showModal(.shiny_compare_picker_modal(choices))
  })

  observeEvent(input$compare_pick_submit, {
    proj <- state$selected_project
    key_a <- input$compare_pick_a %||% ""
    key_b <- input$compare_pick_b %||% ""
    if (!nzchar(key_a) || !nzchar(key_b)) {
      showNotification("Select two runs to compare.", type = "warning")
      return()
    }
    if (identical(key_a, key_b)) {
      showNotification("Choose two different runs.", type = "warning")
      return()
    }
    entries <- .shiny_load_compare_entries(proj, c(key_a, key_b), ws_root)
    if (length(entries) < 2L) {
      showNotification("Could not load one or both runs.", type = "error")
      return()
    }
    state$compare_entries <- entries
    showModal(modalDialog(
      title = "Compare model runs",
      size = "l",
      .shiny_compare_runs_modal_body(entries),
      footer = modalButton("Close"),
      easyClose = TRUE
    ))
  })

  output$compare_gof_grid <- renderPlot({
    entries <- state$compare_entries
    req(length(entries) >= 2L)
    safe_run("Compare GOF plots", {
      if (!.shiny_has_ggplot()) {
        plot.new()
        title("Install ggplot2 for comparison plots")
        return(invisible(NULL))
      }
      plots <- lapply(entries, function(ent) {
        pred <- predict(ent$fit, type = "ipred")
        obs <- .shiny_gof_obs(pred)
        if (nrow(obs) == 0L) {
          return(NULL)
        }
        lims <- range(c(obs$IPRED, obs$DV), finite = TRUE)
        ggplot2::ggplot(obs, ggplot2::aes(x = IPRED, y = DV)) +
          ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray60") +
          ggplot2::geom_point(alpha = 0.55, size = 1.2, color = "steelblue") +
          ggplot2::coord_equal(xlim = lims, ylim = lims) +
          ggplot2::labs(title = ent$label, x = "IPRED", y = "DV") +
          ggplot2::theme_bw(base_size = 10)
      })
      plots <- Filter(Negate(is.null), plots)
      if (length(plots) == 0L) {
        plot.new()
        title("No observations for selected runs")
        return(invisible(NULL))
      }
      if (requireNamespace("patchwork", quietly = TRUE)) {
        patchwork::wrap_plots(plots, ncol = min(2L, length(plots)))
      } else if (requireNamespace("gridExtra", quietly = TRUE)) {
        gridExtra::grid.arrange(grobs = lapply(plots, ggplot2::ggplotGrob), ncol = min(2L, length(plots)))
      } else {
        print(plots[[1L]])
      }
    }) %||% .shiny_empty_plot()
  })

  load_est_run <- function(project, version_id, est_run_id) {
    if (!.valid_id(project) || !.valid_id(version_id) || !.valid_id(est_run_id)) {
      return(invisible(FALSE))
    }
    state$selected_est_run_id <- est_run_id
    state$est_run_loaded_version_id <- version_id
    state$selected_sim_id <- NULL
    state$fit <- nm_workspace_load_run_fit(project, version_id, est_run_id, root = ws_root)
    if (!is.null(state$fit)) {
      if (is.null(state$fit$model) && !is.null(state$model)) {
        state$fit$model <- state$model
      }
      if (is.null(state$fit$data) && !is.null(state$data)) {
        state$fit$data <- state$data
      }
      state$fit <- .shiny_ensure_fit_eta(state$fit)
    }
    state$gof_diag_rev <- state$gof_diag_rev + 1L
    state$auto_data_sim_key <- NULL
    state$auto_data_est_key <- .shiny_est_dataset_key(version_id, est_run_id)
    tbl <- if (is.null(state$fit)) NULL else .shiny_fit_prediction_table(state$fit)
    state$explore_df <- if (is.null(tbl)) NULL else .shiny_as_explore_df(tbl)
    state$explore_data_rev <- (state$explore_data_rev %||% 0L) + 1L
    .shiny_load_run_vpc_state(state, project, version_id, est_run_id, ws_root)
    if (.valid_id(project) && .valid_id(version_id)) {
      parsed <- tryCatch(
        nm_workspace_parse_model(project, version_id, root = ws_root),
        error = function(e) NULL
      )
      if (!is.null(parsed) && !is.null(parsed$data)) {
        state$data <- parsed$data
        state$explore_data_rev <- (state$explore_data_rev %||% 0L) + 1L
      }
    }
    invisible(TRUE)
  }

  load_model_version <- function(project, version_id, keep_est_run = FALSE, keep_sim = FALSE) {
    if (!.valid_id(project) || !.valid_id(version_id)) {
      return(invisible(FALSE))
    }
    txt <- nm_workspace_read_model(project, version_id, root = ws_root)
    state$selected_version_id <- version_id
    if (!isTRUE(keep_sim)) {
      state$selected_sim_id <- NULL
    }
    parts <- nm_ctl_parse(txt)
    push_ctl_to_ui(parts, project, keep_loading = TRUE)
    flags$ctl_load_gen <- flags$ctl_load_gen + 1L
    load_gen <- flags$ctl_load_gen
    session$onFlushed(function() {
      if (!identical(load_gen, flags$ctl_load_gen)) {
        return()
      }
      session$onFlushed(function() {
        if (!identical(load_gen, flags$ctl_load_gen)) {
          return()
        }
        isolate({
          sync_ctl_baseline()
          flags$artifact_warn_shown <- FALSE
        })
        flags$loading_ctl <- FALSE
      }, once = TRUE)
    }, once = TRUE)
    parsed <- tryCatch(
      nm_workspace_parse_model(project, version_id, root = ws_root),
      error = function(e) {
        tryCatch(
          {
            ctl <- nm_ctl_compose(collect_ctl_parts())
            tmp <- tempfile(fileext = ".ctl")
            writeLines(ctl, tmp)
            imp <- nm_read_nonmem(tmp, data_path = state$ctl_data_file)
            fp <- .dataset_full_path(project, state$ctl_data_file)
            if (!is.null(fp)) {
              imp$data <- nm_dataset(fp)
              imp$data_path <- fp
            }
            imp
          },
          error = function(e2) {
            push_log("error", conditionMessage(e), source = "Load model")
            NULL
          }
        )
      }
    )
    if (!is.null(parsed)) {
      state$model <- parsed$model
      state$data <- parsed$data
    } else {
      state$model <- NULL
      state$data <- NULL
    }
    state$explore_df <- NULL
    state$explore_data_rev <- (state$explore_data_rev %||% 0L) + 1L
    if (!is.null(input$data_explore_dataset)) {
      updateSelectInput(session, "data_explore_dataset", selected = "")
    }
    if (!isTRUE(keep_est_run)) {
      runs <- nm_workspace_list_runs(project, version_id, root = ws_root)
      if (nrow(runs) > 0L) {
        load_est_run(project, version_id, runs$run_id[[1L]])
      } else {
        state$selected_est_run_id <- NULL
        state$est_run_loaded_version_id <- NULL
        state$fit <- NULL
      }
    }
    invisible(TRUE)
  }

  show_switch_modal <- function(kind) {
    msg <- switch(
      kind,
      project = "Save changes before switching project?",
      version = "Save changes before switching model version?",
      "Save unsaved model changes?"
    )
    showModal(modalDialog(
      title = "Unsaved model changes",
      msg,
      footer = tagList(
        actionButton("cancel_switch", "Cancel"),
        actionButton("discard_switch", "Discard", class = "btn-warning"),
        actionButton("save_and_switch", "Save", class = "btn-primary")
      ),
      easyClose = FALSE
    ))
  }

  select_project <- function(project) {
    if (!.valid_id(project)) {
      return(invisible(FALSE))
    }
    if (identical(project, state$selected_project)) {
      return(invisible(FALSE))
    }
    if (isTRUE(state$dirty)) {
      flags$pending_project <- project
      flags$pending_version_id <- NULL
      flags$pending_est_run_id <- NULL
      flags$switch_kind <- "project"
      flags$switching <- TRUE
      show_switch_modal("project")
      return(invisible(FALSE))
    }
    state$selected_project <- project
    state$expanded_versions <- character()
    state$selected_sim_id <- NULL
    state$sim_vpc_data <- NULL
    state$sim_vpc_obs <- NULL
    vers <- nm_workspace_list_versions(project, root = ws_root)
    if (length(vers) > 0L) {
      load_model_version(project, vers[[1L]])
    } else {
      state$selected_version_id <- NULL
      state$selected_est_run_id <- NULL
      state$est_run_loaded_version_id <- NULL
      state$selected_sim_id <- NULL
      state$expanded_versions <- character()
      state$model <- NULL
      state$data <- NULL
      state$fit <- NULL
      state$explore_df <- NULL
      state$explore_data_rev <- (state$explore_data_rev %||% 0L) + 1L
      clear_ctl_ui()
    }
    state$versions_rev <- state$versions_rev + 1L
    state$est_runs_rev <- state$est_runs_rev + 1L
    state$simulations_rev <- state$simulations_rev + 1L
    invisible(TRUE)
  }

  select_model_version <- function(version_id) {
    if (!.valid_id(version_id) || !.valid_id(state$selected_project)) {
      return(invisible(FALSE))
    }
    if (identical(version_id, state$selected_version_id)) {
      return(invisible(FALSE))
    }
    if (isTRUE(state$dirty)) {
      flags$pending_version_id <- version_id
      flags$pending_est_run_id <- NULL
      flags$pending_project <- NULL
      flags$switch_kind <- "version"
      flags$switching <- TRUE
      show_switch_modal("version")
      return(invisible(FALSE))
    }
    load_model_version(state$selected_project, version_id)
    state$selected_sim_id <- NULL
    state$est_runs_rev <- state$est_runs_rev + 1L
    invisible(TRUE)
  }

  select_simulation <- function(version_id, sim_id, force_reload = FALSE) {
    proj <- state$selected_project
    if (!.valid_id(proj) || !.valid_id(version_id) || !.valid_id(sim_id)) {
      return(invisible(FALSE))
    }
    if (!identical(version_id, state$selected_version_id)) {
      if (isTRUE(state$dirty)) {
        flags$pending_version_id <- version_id
        flags$pending_est_run_id <- NULL
        flags$pending_project <- NULL
        flags$switch_kind <- "version"
        flags$switching <- TRUE
        show_switch_modal("version")
        return(invisible(FALSE))
      }
      load_model_version(proj, version_id, keep_est_run = TRUE)
    }
    sim_obj <- nm_workspace_load_sim(proj, version_id, sim_id, root = ws_root)
    sim_ds <- .shiny_sim_dataset(sim_obj)
    if (is.null(sim_ds)) {
      push_log("error", paste("Could not load simulation:", sim_id), source = "Simulation")
      return(invisible(FALSE))
    }
    if (is.list(sim_obj) && isTRUE(sim_obj$vpc_mode)) {
      meta_sims <- nm_workspace_list_sims(proj, version_id, root = ws_root)
      row <- meta_sims[meta_sims$sim_id == sim_id, , drop = FALSE]
      est_link <- if (nrow(row) > 0L) row$est_run_id[[1L]] else ""
      if (nzchar(est_link %||% "")) {
        if (!is.null(sim_obj$vpc) && is.data.frame(sim_obj$vpc) && nrow(sim_obj$vpc) > 0L) {
          state$sim_vpc_data <- sim_obj$vpc
          state$sim_vpc_obs <- sim_obj$vpc_obs
          state$vpc_sim_id <- sim_id
        }
        select_est_run(est_link, force_reload = force_reload)
        return(invisible(TRUE))
      }
    }
    state$sim_vpc_data <- NULL
    state$sim_vpc_obs <- NULL
    state$vpc_sim_id <- NULL
    state$selected_sim_id <- sim_id
    state$selected_est_run_id <- NULL
    state$est_run_loaded_version_id <- NULL
    state$fit <- NULL
    state$data <- sim_ds
    tbl <- .shiny_sim_primary_df(sim_obj)
    state$explore_df <- if (is.null(tbl)) NULL else .shiny_as_explore_df(tbl)
    state$auto_data_est_key <- NULL
    state$auto_data_sim_key <- .shiny_sim_dataset_key(version_id, sim_id)
    state$explore_data_rev <- (state$explore_data_rev %||% 0L) + 1L
    showNotification(paste("Loaded simulation:", sim_id), type = "message", duration = 5L)
    invisible(TRUE)
  }

  select_est_run <- function(est_run_id, force_reload = FALSE) {
    if (!.valid_id(est_run_id) || !.valid_id(state$selected_project) ||
        !.valid_id(state$selected_version_id)) {
      return(invisible(FALSE))
    }
    if (!isTRUE(force_reload) &&
        identical(est_run_id, state$selected_est_run_id) &&
        identical(state$selected_version_id, state$est_run_loaded_version_id)) {
      .shiny_load_run_vpc_state(
        state, state$selected_project, state$selected_version_id, est_run_id, ws_root
      )
      return(invisible(TRUE))
    }
    load_est_run(state$selected_project, state$selected_version_id, est_run_id)
    invisible(TRUE)
  }

  apply_est_run_fit <- function(project, version_id, est_run_id, fit) {
    if (!.valid_id(project) || !.valid_id(version_id) || !.valid_id(est_run_id) ||
        is.null(fit)) {
      return(invisible(FALSE))
    }
    if (is.null(fit$model) && !is.null(state$model)) {
      fit$model <- state$model
    }
    if (is.null(fit$data) && !is.null(state$data)) {
      fit$data <- state$data
    }
    fit <- .shiny_ensure_fit_eta(fit)
    state$selected_est_run_id <- est_run_id
    state$est_run_loaded_version_id <- version_id
    state$selected_sim_id <- NULL
    state$fit <- fit
    state$auto_data_sim_key <- NULL
    state$auto_data_est_key <- .shiny_est_dataset_key(version_id, est_run_id)
    tbl <- .shiny_fit_prediction_table(fit)
    state$explore_df <- if (is.null(tbl)) NULL else .shiny_as_explore_df(tbl)
    state$gof_diag_rev <- state$gof_diag_rev + 1L
    state$explore_data_rev <- (state$explore_data_rev %||% 0L) + 1L
    .shiny_load_run_vpc_state(state, project, version_id, est_run_id, ws_root)
    invisible(TRUE)
  }

  observeEvent(input$project_tree_event, {
    if (isTRUE(flags$switching)) {
      return()
    }
    ev <- input$project_tree_event
    if (is.null(ev) || is.null(ev$project)) {
      return()
    }
    project <- ev$project
    if (!.valid_id(project)) {
      return()
    }
    if (identical(project, state$selected_project)) {
      return()
    }
    if (isTRUE(state$dirty)) {
      flags$pending_project <- project
      flags$pending_version_id <- NULL
      flags$pending_est_run_id <- NULL
      flags$switch_kind <- "project"
      flags$switching <- TRUE
      show_switch_modal("project")
      return()
    }
    select_project(project)
  }, ignoreInit = TRUE)

  observeEvent(input$version_tree_event, {
    if (isTRUE(flags$switching)) {
      return()
    }
    ev <- input$version_tree_event
    if (is.null(ev) || is.null(ev$action)) {
      return()
    }
    action <- ev$action
    if (action == "toggle") {
      ver <- ev$version
      if (!.valid_id(ver)) {
        return()
      }
      expanded <- state$expanded_versions
      if (ver %in% expanded) {
        state$expanded_versions <- setdiff(expanded, ver)
      } else {
        state$expanded_versions <- c(expanded, ver)
      }
      return()
    }
    if (action == "select_version") {
      version_id <- ev$version
      if (!.valid_id(version_id)) {
        return()
      }
      if (identical(version_id, state$selected_version_id) &&
          !nzchar(state$selected_sim_id %||% "") &&
          !nzchar(state$selected_est_run_id %||% "")) {
        return()
      }
      if (isTRUE(state$dirty)) {
        flags$pending_version_id <- version_id
        flags$pending_est_run_id <- NULL
        flags$pending_project <- NULL
        flags$switch_kind <- "version"
        flags$switching <- TRUE
        show_switch_modal("version")
        return()
      }
      select_model_version(version_id)
      return()
    }
    if (action == "select_run") {
      version_id <- ev$version
      est_run_id <- ev$run
      if (!.valid_id(version_id) || !.valid_id(est_run_id)) {
        return()
      }
      if (!identical(version_id, state$selected_version_id)) {
        if (isTRUE(state$dirty)) {
          flags$pending_version_id <- version_id
          flags$pending_est_run_id <- est_run_id
          flags$pending_project <- NULL
          flags$switch_kind <- "run"
          flags$switching <- TRUE
          show_switch_modal("version")
          return()
        }
        load_model_version(state$selected_project, version_id, keep_est_run = FALSE)
      }
      select_est_run(est_run_id)
      return()
    }
    if (action == "select_sim") {
      select_simulation(ev$version, ev$sim)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$job_tree_event, {
    ev <- input$job_tree_event
    if (is.null(ev) || is.null(ev$action)) {
      return()
    }
    jid <- ev$job
    if (!.valid_id(jid)) {
      return()
    }
    if (ev$action == "toggle") {
      expanded <- state$expanded_jobs
      if (jid %in% expanded) {
        state$expanded_jobs <- setdiff(expanded, jid)
      } else {
        state$expanded_jobs <- c(expanded, jid)
      }
      return()
    }
    if (ev$action == "select") {
      state$selected_job <- jid
    }
  }, ignoreInit = TRUE)

  observeEvent(input$run_estimation, {
    if (!.valid_id(state$selected_project) || !.valid_id(state$selected_version_id)) {
      showNotification("Select a model version first.", type = "warning")
      return()
    }
    showModal(.shiny_run_estimation_modal())
  })

  observeEvent(input$ctl_trans_help, {
    adv <- as.integer(input$ctl_advan %||% 4L)
    req(nm_ctl_show_trans(adv))
    tr <- as.integer(input$ctl_trans %||% nm_ctl_default_trans(adv))
    showModal(.shiny_ctl_model_help_modal(adv, tr, focus = "trans"))
  })

  output$ctl_problem_display <- renderUI({
    input$ctl_problem
    txt <- trimws(input$ctl_problem %||% "")
    body <- if (!nzchar(txt)) {
      tags$span(class = "text-muted", "No problem statement")
    } else {
      tags$span(txt)
    }
    tags$div(
      class = "ctl-problem-row",
      tags$div(class = "ctl-problem-text", body),
      actionButton(
        "edit_ctl_problem", label = NULL, icon = icon("pencil"),
        class = "btn-xs btn-default", title = "Edit $PROBLEM"
      )
    )
  })

  observeEvent(input$edit_ctl_problem, {
    showModal(modalDialog(
      title = "Edit $PROBLEM",
      textInput("ctl_problem_edit", NULL, value = input$ctl_problem %||% "", width = "100%"),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_ctl_problem", "OK", class = "btn-primary")
      )
    ))
  })

  observeEvent(input$confirm_ctl_problem, {
    updateTextInput(session, "ctl_problem", value = trimws(input$ctl_problem_edit %||% ""))
    removeModal()
    sync_dirty()
  })

  output$ctl_dataset_select <- renderUI({
    state$projects_rev
    state$ctl_data_file
    state$selected_project
    df <- datasets_df()
    if (nrow(df) == 0L) {
      return(tags$span(
        class = "text-muted",
        style = "font-size: 12px;",
        "No datasets in project data/ folder."
      ))
    }
    choices <- stats::setNames(df$path, df$file)
    sel <- state$ctl_data_file
    if (is.null(sel) || !sel %in% choices) {
      sel <- choices[[1L]]
    }
    selectInput("ctl_dataset", NULL, choices = choices, selected = sel, width = "100%")
  })

  observeEvent(input$ctl_dataset, {
    if (isTRUE(flags$loading_ctl)) {
      return()
    }
    rel <- input$ctl_dataset
    if (is.null(rel) || !nzchar(rel)) {
      return()
    }
    if (identical(rel, state$ctl_data_file)) {
      return()
    }
    df <- datasets_df()
    if (nrow(df) == 0L || !rel %in% df$path) {
      return()
    }
    state$ctl_data_file <- rel
    if (.valid_id(state$selected_project)) {
      ds_cols <- .refresh_input_cols(state$selected_project, rel)
      if (length(ds_cols) > 0L) {
        synced <- .sync_input_selection(
          ds_cols,
          current_input = state$ctl_input_cols,
          current_output = state$ctl_output_cols,
          advan = .current_advan(),
          trans = .current_trans()
        )
        .apply_column_selection(
          synced$dataset_columns,
          synced$input_cols,
          synced$output_cols
        )
      }
    }
    state$columns_rev <- state$columns_rev + 1L
    sync_dirty()
  }, ignoreInit = TRUE)

  observeEvent(input$open_column_picker, {
    df <- columns_picker_df()
    showModal(modalDialog(
      title = "$INPUT / $OUTPUT columns",
      tags$p(
        class = "text-muted",
        style = "font-size: 11px; margin-bottom: 8px;",
        "Tick columns to include. Required rows cannot be removed from $INPUT."
      ),
      .shiny_columns_picker_ui(df),
      size = "m",
      easyClose = TRUE,
      footer = modalButton("Done")
    ))
  })

  observeEvent(input$create_simulation, {
    if (!.valid_id(state$selected_project) || !.valid_id(state$selected_version_id)) {
      showNotification("Select a model version first.", type = "warning")
      return()
    }
    ver <- state$selected_version_id
    proj <- state$selected_project
    runs <- nm_workspace_list_runs(proj, ver, root = ws_root)
    has_fit <- nrow(runs) > 0L
    run_choices <- if (has_fit) {
      stats::setNames(
        runs$run_id,
        paste0(runs$run_id, " — ", runs$method)
      )
    } else {
      c("(no estimation runs)" = "")
    }
    max_cores <- parallel::detectCores(logical = TRUE)
    if (is.na(max_cores) || max_cores < 1L) {
      max_cores <- 1L
    }
    default_n_sub <- 10L
    parsed <- tryCatch(
      nm_workspace_parse_model(proj, ver, root = ws_root),
      error = function(e) NULL
    )
    if (!is.null(parsed) && !is.null(parsed$data)) {
      default_n_sub <- max(1L, .shiny_n_subjects(parsed$data))
    }
    showModal(modalDialog(
      title = paste("Create simulation —", ver),
      size = "l",
      fluidRow(
        column(6L, textInput("sim_label", "Label", value = "", width = "100%")),
        column(3L, numericInput("sim_seed", "Random seed", value = sample.int(99999L, 1L), min = 1L, step = 1L, width = "100%")),
        column(3L, numericInput("sim_n_cores", "Parallel cores", value = .shiny_default_n_cores(), min = 1L, max = max_cores, step = 1L, width = "100%"))
      ),
      fluidRow(
        column(4L, numericInput("sim_n_subjects", "Individuals", value = default_n_sub, min = 1L, max = 10000L, step = 1L, width = "100%")),
        column(4L, numericInput("sim_n_replicates", "Replications", value = 1L, min = 1L, max = 1000L, step = 1L, width = "100%")),
        column(4L, numericInput("sim_n_days", "Days (TIME horizon)", value = 1L, min = 1L, max = 365L, step = 1L, width = "100%"))
      ),
      checkboxInput(
        "sim_vpc",
        "Visual predictive check (VPC) — requires multiple replications",
        value = FALSE
      ),
      checkboxInput("sim_use_design", "Custom dosing / sampling design", value = FALSE),
      conditionalPanel(
        condition = "input.sim_use_design",
        fluidRow(
          column(4L, selectInput("sim_dose_mode", "Dosing", choices = c("Single dose" = "single", "Repeat doses" = "repeat", "Steady state" = "steady_state"), width = "100%")),
          column(4L, numericInput("sim_dose_amt", "Default dose amount", value = 320, min = 0, step = 1, width = "100%")),
          column(4L, numericInput("sim_dose_cmt", "Dose CMT", value = 1L, min = 1L, step = 1L, width = "100%"))
        ),
        conditionalPanel(
          condition = "input.sim_dose_mode == 'repeat'",
          fluidRow(
            column(4L, numericInput("sim_dose_n", "Number of doses", value = 3L, min = 1L, max = 100L, step = 1L, width = "100%")),
            column(4L, numericInput("sim_dose_ii", "Inter-dose interval (h)", value = 12, min = 0.1, step = 0.5, width = "100%"))
          )
        ),
        conditionalPanel(
          condition = "input.sim_dose_mode == 'steady_state'",
          numericInput("sim_dose_ii", "Dosing interval II (h)", value = 12, min = 0.1, step = 0.5, width = "100%")
        ),
        textAreaInput(
          "sim_dose_table", "Dose amounts (TIME AMT per line, or AMT only)",
          value = "0 320", rows = 3, width = "100%",
          placeholder = "0 320\n12 320"
        ),
        numericInput("sim_obs_per_day", "Default obs points / day", value = 8L, min = 3L, max = 48L, step = 1L, width = "100%")
      ),
      tags$p(class = "text-muted", style = "font-size: 11px;",
             "Without custom design, the linked dataset structure is kept and resampled to the requested number of individuals."),
      if (has_fit) {
        checkboxInput(
          "sim_use_fit", "Use fitted THETA / OMEGA / SIGMA",
          value = TRUE
        )
      } else {
        tagList(
          tags$p(
            class = "text-muted", style = "font-size: 11px; margin-bottom: 4px;",
            "No estimation run yet — simulation uses THETA / OMEGA / SIGMA from the control file."
          )
        )
      },
      if (has_fit) {
        selectInput(
          "sim_est_run", "Estimation run (optional)",
          choices = c("(most recent fit)" = "", run_choices),
          width = "100%"
        )
      },
      if (has_fit) {
        tagList(
          fluidRow(
            column(
              4L,
              checkboxInput("sim_compute_npc", "Compute NPC on fit", value = FALSE)
            ),
            column(
              4L,
              checkboxInput("sim_compute_npde", "Compute NPDE on fit", value = FALSE)
            ),
            column(
              4L,
              conditionalPanel(
                "input.sim_compute_npc || input.sim_compute_npde",
                checkboxInput(
                  "sim_diag_refit_eta", "Re-fit ETAs on each replicate",
                  value = TRUE
                )
              )
            )
          ),
          conditionalPanel(
            "input.sim_compute_npc || input.sim_compute_npde",
            numericInput(
              "sim_diag_n_sim", "Diagnostic simulations",
              value = 50L, min = 10L, max = 500L, step = 10L, width = "100%"
            ),
            tags$p(
              class = "text-muted", style = "font-size: 11px; margin-top: 0;",
              "NPC/NPDE run in the worker and are linked to the estimation run (like VPC). ",
              "No separate simulation entry is created when only diagnostics are requested."
            )
          )
        )
      },
      tags$p(
        class = "text-muted", style = "font-size: 11px;",
        "Simulates observations from the model and dataset linked to this version. ",
        "Runs as a background job. Monitor progress in the log banner or on the Jobs tab."
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_create_simulation", "Run simulation", class = "btn-primary")
      ),
      easyClose = TRUE
    ))
  })

  .launch_simulation_job <- function() {
    pending <- flags$pending_sim
    if (is.null(pending)) {
      return()
    }
    flags$pending_sim <- NULL
    safe_run("Create simulation", {
      .shiny_start_simulation_job(
        session = session,
        state = state,
        proj = pending$proj,
        ver = pending$ver,
        label = pending$label,
        seed = pending$seed,
        n_sim = pending$n_sim,
        use_fit = pending$use_fit,
        est_run = pending$est_run,
        n_cores = pending$n_cores,
        design = pending$design,
        vpc = pending$vpc,
        sim_compute_npc = pending$sim_compute_npc,
        sim_compute_npde = pending$sim_compute_npde,
        diag_n_sim = pending$diag_n_sim,
        diag_refit_eta = pending$diag_refit_eta,
        ws_root = ws_root,
        job_root = job_root
      )
    })
  }

  observeEvent(input$confirm_create_simulation, {
    removeModal()
    ver <- state$selected_version_id
    proj <- state$selected_project
    if (!.valid_id(proj) || !.valid_id(ver)) {
      return()
    }
    est_run <- input$sim_est_run %||% ""
    if (!nzchar(est_run)) {
      est_run <- NULL
    }
    has_fit <- nrow(nm_workspace_list_runs(proj, ver, root = ws_root)) > 0L
    sim_compute_npc <- has_fit && isTRUE(input$sim_compute_npc)
    sim_compute_npde <- has_fit && isTRUE(input$sim_compute_npde)
    run_diag <- sim_compute_npc || sim_compute_npde
    use_fit <- has_fit && (isTRUE(input$sim_use_fit) || run_diag)
    if (run_diag && !isTRUE(input$sim_use_fit)) {
      showNotification(
        "NPC/NPDE requires fitted THETA / OMEGA / SIGMA; using the selected estimation run.",
        type = "message",
        duration = 6L
      )
    }
    n_sim <- max(1L, as.integer(input$sim_n_replicates))
    if (!isTRUE(input$sim_vpc) && !run_diag && n_sim <= 1L) {
      showNotification(
        "Enable VPC, NPC, NPDE, or set replications > 1.",
        type = "warning",
        duration = 8L
      )
      return()
    }
    if (isTRUE(input$sim_vpc) && n_sim < 10L) {
      showNotification("VPC recommends at least 10 replications; using 10.", type = "warning", duration = 8L)
      n_sim <- max(10L, n_sim)
    }
    seed <- as.integer(input$sim_seed)
    design <- .shiny_sim_build_design(input, seed = seed)
    pending <- list(
      proj = proj,
      ver = ver,
      label = trimws(input$sim_label %||% ""),
      seed = seed,
      n_sim = n_sim,
      use_fit = use_fit,
      est_run = est_run,
      n_cores = max(1L, as.integer(input$sim_n_cores %||% 1L)),
      design = design,
      vpc = isTRUE(input$sim_vpc),
      sim_compute_npc = sim_compute_npc,
      sim_compute_npde = sim_compute_npde,
      diag_n_sim = as.integer(input$sim_diag_n_sim %||% 50L),
      diag_refit_eta = isTRUE(input$sim_diag_refit_eta)
    )
    wl <- tryCatch({
      parsed <- nm_workspace_parse_model(proj, ver, root = ws_root)
      if (is.null(parsed$model) || is.null(parsed$data)) {
        NULL
      } else {
        tmpl <- nm_sim_template_data(parsed$model, parsed$data, design)
        nm_sim_workload(parsed$model, tmpl, n_sim = n_sim)
      }
    }, error = function(e) NULL)
    if (!is.null(wl) && wl$total_points >= .shiny_sim_workload_threshold) {
      flags$pending_sim <- pending
      showModal(modalDialog(
        title = "Large simulation",
        size = "m",
        tags$p(
          "This simulation will evaluate approximately ",
          strong(format(wl$total_points, big.mark = ",", scientific = FALSE)),
          " simulated points ",
          "(subjects × observation rows × replicates)."
        ),
        tags$ul(
          tags$li(paste0(wl$n_subjects, " subject(s)")),
          tags$li(paste0(wl$n_obs_rows, " observation row(s) per subject")),
          tags$li(paste0(wl$n_replicates, " replicate(s)"))
        ),
        tags$p(
          class = "text-warning",
          "This may take a long time and use substantial memory. ",
          "Click Cancel to adjust settings, or Continue to submit the job."
        ),
        footer = tagList(
          actionButton("cancel_sim_large", "Cancel", class = "btn-default"),
          actionButton("confirm_sim_large", "Continue", class = "btn-warning")
        ),
        easyClose = FALSE
      ))
      return()
    }
    flags$pending_sim <- pending
    .launch_simulation_job()
  })

  observeEvent(input$confirm_sim_large, {
    removeModal()
    .launch_simulation_job()
  })

  observeEvent(input$cancel_sim_large, {
    removeModal()
    flags$pending_sim <- NULL
  })

  lapply(
    c("ctl_problem", "ctl_advan", "ctl_trans", "ctl_pk", "ctl_des", "ctl_error", "ctl_ode_ncomp"),
    function(id) {
      observeEvent(input[[id]], {
        if (isTRUE(flags$loading_ctl)) {
          return()
        }
        sync_dirty()
      }, ignoreInit = TRUE)
    }
  )

  observeEvent(input$ctl_advan, {
    if (isTRUE(flags$loading_ctl)) {
      return()
    }
    advan <- input$ctl_advan
    .sync_advan_input_cols()
    if (nm_ctl_use_ode(as.integer(advan))) {
      .shiny_refresh_ctl_des(advan)
    } else {
      cur_des <- input$ctl_des %||% ""
      if (!nzchar(trimws(cur_des))) {
        updateTextAreaInput(session, "ctl_des", value = "")
      }
    }
    if (nm_ctl_show_trans(as.integer(advan))) {
      trans_named <- .shiny_numeric_trans_choices(advan)
      if (!is.null(trans_named)) {
        sel <- input$ctl_trans
        if (is.null(sel) || !as.character(sel) %in% names(trans_named)) {
          sel <- nm_ctl_default_trans(advan)
        }
        updateSelectInput(
          session, "ctl_trans",
          choices = trans_named,
          selected = sel
        )
      }
    }
    sync_dirty()
  }, ignoreInit = TRUE)

  observeEvent(input$ctl_ode_ncomp, {
    if (isTRUE(flags$loading_ctl)) {
      return()
    }
    if (.current_advan() == 6L) {
      .shiny_refresh_ctl_des("6")
      sync_dirty()
    }
  }, ignoreInit = TRUE)

  observeEvent(input$ctl_trans, {
    if (isTRUE(flags$loading_ctl)) {
      return()
    }
    .sync_advan_input_cols()
    sync_dirty()
  }, ignoreInit = TRUE)

  save_current_model <- function() {
    req(state$selected_project, state$selected_version_id)
    parts <- collect_ctl_parts()
    ctl <- nm_ctl_compose(parts)
    nm_workspace_write_model(
      state$selected_project,
      state$selected_version_id,
      ctl,
      root = ws_root,
      data_file = parts$data_file
    )
    state$saved_text <- nm_ctl_canonical(ctl)
    state$dirty <- FALSE
    flags$artifact_warn_shown <- FALSE
    parsed <- tryCatch(
      {
        tmp <- tempfile(fileext = ".ctl")
        writeLines(ctl, tmp)
        imp <- nm_read_nonmem(tmp, data_path = parts$data_file)
        fp <- .dataset_full_path(state$selected_project, parts$data_file)
        if (!is.null(fp)) {
          imp$data <- nm_dataset(fp)
          imp$data_path <- fp
        }
        imp
      },
      error = function(e) NULL
    )
    if (is.null(parsed)) {
      parsed <- tryCatch(
        nm_workspace_parse_model(state$selected_project, state$selected_version_id, root = ws_root),
        error = function(e) NULL
      )
    }
    if (!is.null(parsed)) {
      state$model <- parsed$model
      state$data <- parsed$data
    }
    state$versions_rev <- state$versions_rev + 1L
    TRUE
  }

  observeEvent(input$save_model, {
    if (save_current_model()) {
      showNotification("Model overwritten.", type = "message")
    }
  })

  observeEvent(input$save_model_as_new, {
    req(state$selected_project, state$selected_version_id)
    showModal(modalDialog(
      title = "Save as new version",
      size = "m",
      textInput(
        "save_as_new_label", "Label",
        value = paste0(state$selected_version_id, " (edited)"),
        width = "100%"
      ),
      textInput(
        "save_as_new_id", "Version id (optional, auto if blank)",
        value = "", width = "100%"
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_save_as_new", "Save", class = "btn-primary")
      )
    ))
  })

  observeEvent(input$confirm_save_as_new, {
    removeModal()
    req(state$selected_project)
    parts <- collect_ctl_parts()
    ctl <- nm_ctl_compose(parts)
    new_id <- trimws(input$save_as_new_id %||% "")
    if (!nzchar(new_id)) {
      new_id <- NULL
    }
    label <- trimws(input$save_as_new_label %||% "")
    new_ver <- nm_workspace_new_version(
      state$selected_project,
      version_id = new_id,
      root = ws_root,
      template_ctl = ctl,
      label = if (nzchar(label)) label else NULL,
      data_file = parts$data_file
    )
    state$saved_text <- ctl
    state$dirty <- FALSE
    state$versions_rev <- state$versions_rev + 1L
    load_model_version(state$selected_project, new_ver, keep_est_run = FALSE, keep_sim = FALSE)
    showNotification(paste("Saved new version:", new_ver), type = "message")
  })

  observeEvent(input$delete_simulation, {
    proj <- state$selected_project
    ver <- state$selected_version_id
    sim_id <- state$selected_sim_id
    if (!.valid_id(proj) || !.valid_id(ver) || !.valid_id(sim_id)) {
      showNotification("Select a simulation to delete.", type = "warning")
      return()
    }
    showModal(modalDialog(
      title = "Delete simulation",
      paste0("Delete simulation '", sim_id, "' from version '", ver, "'?"),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_delete_simulation", "Delete", class = "btn-danger")
      )
    ))
  })

  observeEvent(input$confirm_delete_simulation, {
    removeModal()
    proj <- state$selected_project
    ver <- state$selected_version_id
    sim_id <- state$selected_sim_id
    if (!.valid_id(proj) || !.valid_id(ver) || !.valid_id(sim_id)) {
      return()
    }
    nm_workspace_delete_sim(proj, ver, sim_id, root = ws_root)
    if (identical(sim_id, state$vpc_sim_id %||% "")) {
      if (.valid_id(state$selected_est_run_id)) {
        .shiny_load_run_vpc_state(state, proj, ver, state$selected_est_run_id, ws_root)
      } else {
        state$sim_vpc_data <- NULL
        state$sim_vpc_obs <- NULL
        state$vpc_sim_id <- NULL
      }
    }
    if (identical(state$auto_data_sim_key, .shiny_sim_dataset_key(ver, sim_id))) {
      state$auto_data_sim_key <- NULL
      state$explore_df <- NULL
    }
    state$selected_sim_id <- NULL
    if (!.valid_id(state$selected_est_run_id)) {
      state$sim_vpc_data <- NULL
      state$sim_vpc_obs <- NULL
      state$vpc_sim_id <- NULL
    }
    state$simulations_rev <- state$simulations_rev + 1L
    state$explore_data_rev <- (state$explore_data_rev %||% 0L) + 1L
    showNotification(paste("Deleted simulation:", sim_id), type = "message")
  })

  observeEvent(input$ctl_col_picker_event, {
    if (isTRUE(flags$loading_ctl)) {
      return()
    }
    ev <- input$ctl_col_picker_event
    req(ev, nzchar(ev$col))
    df <- columns_picker_df()
    if (nrow(df) == 0L) {
      return()
    }
    idx <- match(ev$col, df$Column)
    if (is.na(idx)) {
      return()
    }
    essential <- intersect(
      nm_ctl_essential_input_cols(.current_advan(), .current_trans()),
      state$dataset_columns
    )
    if (identical(ev$kind, "in")) {
      if (ev$col %in% essential && !isTRUE(ev$value)) {
        showNotification(
          paste(ev$col, "is required and must stay in $INPUT."),
          type = "warning", duration = 3L
        )
        state$columns_rev <- state$columns_rev + 1L
        return()
      }
      if (!ev$col %in% state$dataset_columns && isTRUE(ev$value)) {
        showNotification(
          paste(ev$col, "is not in the dataset and cannot be added to $INPUT."),
          type = "warning", duration = 3L
        )
        state$columns_rev <- state$columns_rev + 1L
        return()
      }
      df$IN[idx] <- isTRUE(ev$value)
    } else if (identical(ev$kind, "out")) {
      df$OUT[idx] <- isTRUE(ev$value)
    }
    .update_cols_from_picker(df)
    sync_dirty()
  }, ignoreInit = TRUE)

  output$ctl_theta_table <- DT::renderDataTable({
    if (is.null(state$ctl_thetas) || nrow(state$ctl_thetas) == 0L) {
      return(DT::datatable(
        data.frame(note = "No THETA records."),
        options = list(dom = "t"), rownames = FALSE, selection = "none"
      ))
    }
    df <- state$ctl_thetas
    display <- data.frame(
      Name = df$Label,
      THETA = df$THETA,
      Lower = df$Lower,
      Initial = df$Value,
      Upper = df$Upper,
      FIX = ifelse(df$FIX, 1L, 0L),
      stringsAsFactors = FALSE
    )
    DT::datatable(
      display,
      editable = list(target = "cell"),
      selection = "none",
      options = .dt_left_opts(list(dom = "t", ordering = FALSE)),
      rownames = FALSE,
      class = "compact ctl-dt-left"
    )
  })

  output$ctl_omega_table <- DT::renderDataTable({
    if (is.null(state$ctl_omegas) || nrow(state$ctl_omegas) == 0L) {
      return(DT::datatable(
        data.frame(note = "No OMEGA records."),
        options = list(dom = "t"), rownames = FALSE, selection = "none"
      ))
    }
    df <- .shiny_ctl_param_labels(state$ctl_omegas, "OMEGA")
    display <- data.frame(
      Name = df$Label,
      OMEGA = df$OMEGA,
      Initial = df$Value,
      stringsAsFactors = FALSE
    )
    DT::datatable(
      display,
      editable = list(target = "cell"),
      selection = "none",
      options = .dt_left_opts(list(dom = "t", ordering = FALSE, scrollY = "100px", scrollCollapse = TRUE)),
      rownames = FALSE,
      class = "compact ctl-dt-left"
    )
  })

  output$ctl_sigma_table <- DT::renderDataTable({
    if (is.null(state$ctl_sigmas) || nrow(state$ctl_sigmas) == 0L) {
      return(DT::datatable(
        data.frame(note = "No SIGMA records."),
        options = list(dom = "t"), rownames = FALSE, selection = "none"
      ))
    }
    df <- .shiny_ctl_param_labels(state$ctl_sigmas, "SIGMA")
    display <- data.frame(
      Name = df$Label,
      SIGMA = df$SIGMA,
      Initial = df$Value,
      stringsAsFactors = FALSE
    )
    DT::datatable(
      display,
      editable = list(target = "cell"),
      selection = "none",
      options = .dt_left_opts(list(dom = "t", ordering = FALSE, scrollY = "100px", scrollCollapse = TRUE)),
      rownames = FALSE,
      class = "compact ctl-dt-left"
    )
  })

  observeEvent(input$ctl_theta_table_cell_edit, {
    info <- input$ctl_theta_table_cell_edit
    req(state$ctl_thetas, info)
    # DT with rownames = FALSE: info$row is 1-based; info$col is 0-based (DT PR #480)
    row_idx <- as.integer(info$row)
    col_idx <- as.integer(info$col) + 1L
    if (row_idx < 1L || row_idx > nrow(state$ctl_thetas)) {
      return()
    }
    if (col_idx == 2L) {
      return()
    }
    df <- state$ctl_thetas
    val <- info$value
    if (col_idx == 1L) {
      df$Label[row_idx] <- trimws(as.character(val))
    } else if (col_idx == 3L) {
      df$Lower[row_idx] <- as.numeric(val)
    } else if (col_idx == 4L) {
      df$Value[row_idx] <- as.numeric(val)
    } else if (col_idx == 5L) {
      df$Upper[row_idx] <- as.numeric(val)
    } else if (col_idx == 6L) {
      df$FIX[row_idx] <- as.integer(val) != 0L
    }
    state$ctl_thetas <- df
    sync_dirty()
  })

  observeEvent(input$ctl_omega_table_cell_edit, {
    info <- input$ctl_omega_table_cell_edit
    req(state$ctl_omegas, info)
    row_idx <- as.integer(info$row)
    col_idx <- as.integer(info$col) + 1L
    if (row_idx < 1L || row_idx > nrow(state$ctl_omegas)) {
      return()
    }
    if (col_idx == 2L) {
      return()
    }
    df <- .shiny_ctl_param_labels(state$ctl_omegas, "OMEGA")
    if (col_idx == 1L) {
      df$Label[row_idx] <- trimws(as.character(info$value))
    } else if (col_idx == 3L) {
      df$Value[row_idx] <- as.numeric(info$value)
    }
    state$ctl_omegas <- df
    sync_dirty()
  })

  observeEvent(input$ctl_sigma_table_cell_edit, {
    info <- input$ctl_sigma_table_cell_edit
    req(state$ctl_sigmas, info)
    row_idx <- as.integer(info$row)
    col_idx <- as.integer(info$col) + 1L
    if (row_idx < 1L || row_idx > nrow(state$ctl_sigmas)) {
      return()
    }
    if (col_idx == 2L) {
      return()
    }
    df <- .shiny_ctl_param_labels(state$ctl_sigmas, "SIGMA")
    if (col_idx == 1L) {
      df$Label[row_idx] <- trimws(as.character(info$value))
    } else if (col_idx == 3L) {
      df$Value[row_idx] <- as.numeric(info$value)
    }
    state$ctl_sigmas <- df
    sync_dirty()
  })

  observeEvent(input$delete_project, {
    if (!.valid_id(state$selected_project)) {
      showNotification("Select a project to delete.", type = "warning")
      return()
    }
    showModal(modalDialog(
      title = "Delete project",
      tagList(
        tags$p(paste0(
          "Permanently delete project '", state$selected_project,
          "' and all model versions, runs, and data?"
        )),
        tags$p(class = "text-muted", style = "font-size: 12px;",
               "Close any open files in this project folder before deleting."),
        textInput(
          "delete_project_confirm", "Type YES to confirm",
          value = "", placeholder = "YES", width = "100%"
        )
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton(
          "confirm_delete_project", "Delete",
          class = "btn-danger", disabled = TRUE
        )
      )
    ))
  })

  observeEvent(input$delete_project_confirm, {
    disable <- !identical(trimws(input$delete_project_confirm %||% ""), "YES")
    updateActionButton(session, "confirm_delete_project", disabled = disable)
  }, ignoreNULL = FALSE)

  observeEvent(input$confirm_delete_project, {
    if (!identical(trimws(input$delete_project_confirm %||% ""), "YES")) {
      showNotification("Type YES in the confirmation box to delete the project.", type = "warning")
      return()
    }
    removeModal()
    proj <- state$selected_project
    req(.valid_id(proj))
    tryCatch(
      {
        nm_workspace_delete_project(proj, root = ws_root)
        state$selected_project <- NULL
        state$selected_version_id <- NULL
        state$selected_est_run_id <- NULL
        state$est_run_loaded_version_id <- NULL
        state$selected_sim_id <- NULL
        state$expanded_versions <- character()
        state$model <- NULL
        state$data <- NULL
        state$fit <- NULL
        state$explore_df <- NULL
        state$saved_text <- ""
        state$dirty <- FALSE
        clear_ctl_ui()
        state$projects_rev <- state$projects_rev + 1L
        state$versions_rev <- state$versions_rev + 1L
        state$est_runs_rev <- state$est_runs_rev + 1L
        state$simulations_rev <- state$simulations_rev + 1L
        state$explore_data_rev <- (state$explore_data_rev %||% 0L) + 1L
        remaining <- nm_workspace_list_projects(ws_root)
        if (length(remaining) > 0L) {
          select_project(remaining[[1L]])
        }
        showNotification(paste("Deleted project:", proj), type = "message")
      },
      error = function(e) {
        msg <- conditionMessage(e)
        if (grepl("Could not delete project folder", msg, fixed = TRUE)) {
          showNotification(msg, type = "error", duration = NULL)
        } else {
          showNotification(msg, type = "error")
        }
      }
    )
  })

  observeEvent(input$delete_version, {
    if (!.valid_id(state$selected_project) || !.valid_id(state$selected_version_id)) {
      showNotification("Select a model version to delete.", type = "warning")
      return()
    }
    showModal(modalDialog(
      title = "Delete model version",
      paste0("Delete version '", state$selected_version_id,
             "' and all its estimation runs?"),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_delete_version", "Delete", class = "btn-danger")
      )
    ))
  })

  observeEvent(input$confirm_delete_version, {
    removeModal()
    proj <- state$selected_project
    ver <- state$selected_version_id
    req(.valid_id(proj), .valid_id(ver))
    tryCatch(
      {
        nm_workspace_delete_version(proj, ver, root = ws_root)
        vers <- nm_workspace_list_versions(proj, root = ws_root)
        if (length(vers) > 0L) {
          load_model_version(proj, vers[[1L]])
        } else {
        state$selected_version_id <- NULL
        state$selected_est_run_id <- NULL
        state$est_run_loaded_version_id <- NULL
        state$model <- NULL
          state$data <- NULL
          state$fit <- NULL
          clear_ctl_ui()
        }
        state$versions_rev <- state$versions_rev + 1L
        state$est_runs_rev <- state$est_runs_rev + 1L
        showNotification(paste("Deleted version:", ver), type = "message")
      },
      error = function(e) showNotification(conditionMessage(e), type = "error")
    )
  })

  observeEvent(input$delete_est_run, {
    if (!.valid_id(state$selected_project) || !.valid_id(state$selected_version_id) ||
        !.valid_id(state$selected_est_run_id)) {
      showNotification("Select an estimation run to delete.", type = "warning")
      return()
    }
    showModal(modalDialog(
      title = "Delete estimation run",
      paste0("Delete run '", state$selected_est_run_id, "'?"),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_delete_est_run", "Delete", class = "btn-danger")
      )
    ))
  })

  observeEvent(input$confirm_delete_est_run, {
    removeModal()
    proj <- state$selected_project
    ver <- state$selected_version_id
    run <- state$selected_est_run_id
    req(.valid_id(proj), .valid_id(ver), .valid_id(run))
    tryCatch(
      {
        nm_workspace_delete_run(proj, ver, run, root = ws_root)
        runs <- nm_workspace_list_runs(proj, ver, root = ws_root)
        if (nrow(runs) > 0L) {
          load_est_run(proj, ver, runs$run_id[[1L]])
        } else {
          state$selected_est_run_id <- NULL
          state$est_run_loaded_version_id <- NULL
          state$fit <- NULL
        }
        state$est_runs_rev <- state$est_runs_rev + 1L
        state$versions_rev <- state$versions_rev + 1L
        showNotification(paste("Deleted run:", run), type = "message")
      },
      error = function(e) showNotification(conditionMessage(e), type = "error")
    )
  })

  observeEvent(input$cancel_switch, {
    removeModal()
    flags$pending_version_id <- NULL
    flags$pending_est_run_id <- NULL
    flags$pending_project <- NULL
    flags$switch_kind <- NULL
    flags$switching <- FALSE
  })

  observeEvent(input$save_and_switch, {
    removeModal()
    save_current_model()
    kind <- flags$switch_kind
    target_version <- flags$pending_version_id
    target_proj <- flags$pending_project
    target_run <- flags$pending_est_run_id
    flags$pending_version_id <- NULL
    flags$pending_est_run_id <- NULL
    flags$pending_project <- NULL
    flags$switch_kind <- NULL
    flags$switching <- FALSE
    flags$artifact_warn_shown <- FALSE
    if (identical(kind, "project") && !is.null(target_proj)) {
      state$selected_project <- target_proj
      state$projects_rev <- state$projects_rev + 1L
      vers <- nm_workspace_list_versions(target_proj, root = ws_root)
      if (length(vers) > 0L) {
        load_model_version(target_proj, vers[[1L]])
      }
    } else if (!is.null(target_version)) {
      load_model_version(
        state$selected_project, target_version,
        keep_est_run = !is.null(target_run)
      )
      if (!is.null(target_run)) {
        select_est_run(target_run)
      }
    }
    state$versions_rev <- state$versions_rev + 1L
    state$est_runs_rev <- state$est_runs_rev + 1L
  })

  observeEvent(input$discard_switch, {
    removeModal()
    kind <- flags$switch_kind
    target_version <- flags$pending_version_id
    target_proj <- flags$pending_project
    target_run <- flags$pending_est_run_id
    flags$pending_version_id <- NULL
    flags$pending_est_run_id <- NULL
    flags$pending_project <- NULL
    flags$switch_kind <- NULL
    flags$switching <- FALSE
    flags$artifact_warn_shown <- FALSE
    if (identical(kind, "project") && !is.null(target_proj)) {
      state$selected_project <- target_proj
      vers <- nm_workspace_list_versions(target_proj, root = ws_root)
      if (length(vers) > 0L) {
        load_model_version(target_proj, vers[[1L]])
      } else {
        state$selected_version_id <- NULL
        state$selected_est_run_id <- NULL
        state$est_run_loaded_version_id <- NULL
        clear_ctl_ui()
      }
    } else if (!is.null(target_version)) {
      load_model_version(
        state$selected_project, target_version,
        keep_est_run = !is.null(target_run)
      )
      if (!is.null(target_run)) {
        select_est_run(target_run)
      }
    }
    state$versions_rev <- state$versions_rev + 1L
    state$est_runs_rev <- state$est_runs_rev + 1L
  })

  observeEvent(input$reload_model, {
    req(state$selected_project, state$selected_version_id)
    if (isTRUE(state$dirty)) {
      showModal(modalDialog(
        title = "Reload model",
        "Discard unsaved changes and reload from disk?",
        footer = tagList(
          modalButton("Cancel"),
          actionButton("confirm_reload", "Reload", class = "btn-warning")
        )
      ))
    } else {
      load_model_version(state$selected_project, state$selected_version_id, keep_est_run = TRUE)
      if (.valid_id(state$selected_est_run_id)) {
        load_est_run(
          state$selected_project,
          state$selected_version_id,
          state$selected_est_run_id
        )
      }
      showNotification("Model reloaded.", type = "message")
    }
  })

  .advan_select_choices <- function() {
    setNames(nm_ctl_advan_choices(), .shiny_advan_choice_labels())
  }

  .sync_advan_trans_inputs <- function(advan_id, trans_id, adv = NULL) {
    adv <- adv %||% input[[advan_id]]
    if (is.null(adv) || !nzchar(as.character(adv))) {
      return(invisible(NULL))
    }
    if (!nm_ctl_show_trans(as.integer(adv))) {
      return(invisible(NULL))
    }
    trans_named <- .shiny_named_trans_choices(adv)
    if (is.null(trans_named)) {
      return(invisible(NULL))
    }
    sel <- input[[trans_id]]
    if (is.null(sel) || !as.character(sel) %in% names(trans_named)) {
      sel <- nm_ctl_default_trans(adv)
    }
    updateSelectInput(
      session, trans_id,
      choices = trans_named,
      selected = sel
    )
  }

  .ctl_template_model_inputs <- function(advan_id, trans_id, label_id, problem_id, advan = "4",
                                         ode_ncomp_id = NULL) {
    trans <- nm_ctl_default_trans(advan)
    ode_ui <- if (!is.null(ode_ncomp_id)) {
      conditionalPanel(
        condition = sprintf("input['%s'] == '6'", advan_id),
        numericInput(
          ode_ncomp_id, "Compartments (#)",
          value = 2L, min = 1L, max = 10L, step = 1L, width = "100%"
        )
      )
    } else {
      NULL
    }
    tagList(
      textInput(
        label_id, "Version label (optional)",
        value = "", width = "100%"
      ),
      textInput(
        problem_id, "Problem statement",
        value = "Template model", width = "100%"
      ),
      selectInput(
        advan_id, "ADVAN",
        choices = .advan_select_choices(), selected = advan, width = "100%"
      ),
      conditionalPanel(
        condition = sprintf("input['%s'] != '6'", advan_id),
        selectInput(
          trans_id, "TRANS",
          choices = stats::setNames(
            nm_ctl_trans_choices(advan),
            paste0("TRANS", nm_ctl_trans_choices(advan))
          ),
          selected = trans,
          width = "100%"
        )
      ),
      ode_ui,
      tags$p(
        class = "text-muted", style = "font-size: 11px; margin-top: 8px;",
        "Builds standard $THETA, $OMEGA, $SIGMA, $PK, and $ERROR blocks for the selected ADVAN/TRANS pair."
      )
    )
  }

  .show_new_project_modal <- function() {
    showModal(modalDialog(
      title = "New project",
      size = "m",
      textInput("new_project_name", "Project name", value = "", width = "100%"),
      textAreaInput(
        "new_project_description", "Description (optional)",
        value = "", width = "100%", rows = 2, resize = "vertical"
      ),
      radioButtons(
        "new_project_mode", "Project type",
        choices = c(
          "Empty project" = "empty",
          "Create from template" = "template"
        ),
        selected = "empty",
        inline = TRUE
      ),
      conditionalPanel(
        condition = "input.new_project_mode == 'template'",
        tags$hr(style = "margin: 12px 0;"),
        tags$strong("Initial model version"),
        tags$p(
          class = "text-muted", style = "font-size: 11px; margin-bottom: 8px;",
          "Imports or generates a dataset and creates the first model version in the new project."
        ),
        radioButtons(
          "project_tpl_data_source", "Dataset",
          choices = c(
            "Built-in synthetic example" = "synthetic",
            "Upload CSV file" = "upload"
          ),
          selected = "synthetic"
        ),
        conditionalPanel(
          condition = "input.project_tpl_data_source == 'synthetic'",
          selectInput(
            "project_tpl_synthetic_id", "Example",
            choices = .shiny_synthetic_tpl_choices(),
            selected = "theo",
            width = "100%"
          ),
          numericInput(
            "project_tpl_theo_nsub", "Number of subjects",
            value = 10L, min = 1L, max = 500L, step = 1L, width = "100%"
          ),
          uiOutput("project_tpl_synthetic_desc")
        ),
        conditionalPanel(
          condition = "input.project_tpl_data_source == 'upload'",
          fileInput(
            "project_tpl_dataset_file", "Dataset file",
            accept = c(".csv", ".txt", ".dat", ".tsv"),
            width = "100%"
          ),
          tags$p(
            class = "text-muted", style = "font-size: 11px; margin-bottom: 8px;",
            "Expected NONMEM-style columns (ID, TIME, DV, AMT, EVID, CMT, MDV, …)."
          ),
          .ctl_template_model_inputs(
            advan_id = "project_tpl_advan",
            trans_id = "project_tpl_trans",
            label_id = "project_tpl_label",
            problem_id = "project_tpl_problem",
            advan = "4",
            ode_ncomp_id = "project_tpl_ode_ncomp"
          )
        ),
        conditionalPanel(
          condition = "input.project_tpl_data_source == 'synthetic'",
          textInput(
            "project_tpl_label", "Version label (optional)",
            value = "", width = "100%"
          ),
          textInput(
            "project_tpl_problem", "Problem statement",
            value = "Synthetic demo", width = "100%"
          )
        )
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_create_project", "Create project", class = "btn-primary")
      ),
      easyClose = TRUE
    ))
  }

  .require_clean_for_new_version <- function() {
    if (!.valid_id(state$selected_project)) {
      showNotification("Select a project first.", type = "warning")
      return(FALSE)
    }
    if (isTRUE(state$dirty)) {
      showNotification("Save or discard changes before creating a new version.", type = "warning")
      return(FALSE)
    }
    TRUE
  }

  .finish_new_version <- function(version_id) {
    state$versions_rev <- state$versions_rev + 1L
    load_model_version(state$selected_project, version_id)
    state$est_runs_rev <- state$est_runs_rev + 1L
    showNotification(paste("Created model version:", version_id), type = "message")
  }

  .show_template_version_modal <- function() {
    proj <- state$selected_project
    ds <- nm_workspace_list_datasets(proj, root = ws_root)
    ds_choices <- if (nrow(ds) > 0L) {
      stats::setNames(ds$path, ds$file)
    } else {
      c("(no datasets in project data/ folder)" = "")
    }
    default_ds <- if (nrow(ds) > 0L) ds$path[[1L]] else ""
    if (.valid_id(state$selected_version_id)) {
      cur_meta <- .version_meta(proj, state$selected_version_id)
      if (!is.null(cur_meta$data_file) && nzchar(as.character(cur_meta$data_file))) {
        default_ds <- as.character(cur_meta$data_file)
      }
    }
    advan <- "4"
    showModal(modalDialog(
      title = "New version from template",
      size = "m",
      selectInput(
        "template_dataset", "Dataset",
        choices = ds_choices, selected = default_ds, width = "100%"
      ),
      .ctl_template_model_inputs(
        advan_id = "template_advan",
        trans_id = "template_trans",
        label_id = "template_label",
        problem_id = "template_problem",
        advan = advan,
        ode_ncomp_id = "template_ode_ncomp"
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_template_version", "Create version", class = "btn-primary")
      ),
      easyClose = TRUE
    ))
  }

  observeEvent(input$confirm_reload, {
    removeModal()
    load_model_version(state$selected_project, state$selected_version_id, keep_est_run = TRUE)
    if (.valid_id(state$selected_est_run_id)) {
      load_est_run(
        state$selected_project,
        state$selected_version_id,
        state$selected_est_run_id
      )
    }
    showNotification("Model reloaded.", type = "message")
  })

  observeEvent(input$new_version_copy, {
    if (!.require_clean_for_new_version()) {
      return()
    }
    if (!.valid_id(state$selected_version_id)) {
      showNotification("Select a model version to copy.", type = "warning")
      return()
    }
    has_fit <- !is.null(current_fit())
    showModal(.shiny_copy_version_modal(has_fit = has_fit))
  })

  observeEvent(input$confirm_copy_version, {
    removeModal()
    if (!.valid_id(state$selected_version_id)) {
      return()
    }
    fit_inits <- NULL
    if (isTRUE(input$copy_update_inits)) {
      fit_inits <- current_fit()
      if (is.null(fit_inits)) {
        showNotification(
          "No fit loaded — copying without updated initials.",
          type = "warning"
        )
      }
    }
    safe_run("Copy to new", {
      version_id <- nm_workspace_copy_version(
        state$selected_project,
        state$selected_version_id,
        root = ws_root,
        fit_inits = fit_inits
      )
      .finish_new_version(version_id)
      if (!is.null(fit_inits)) {
        showNotification("Copied with updated parameter initials.", type = "message")
      }
    })
  })

  observeEvent(input$new_version_template, {
    if (!.require_clean_for_new_version()) {
      return()
    }
    .show_template_version_modal()
  })

  observeEvent(input$template_advan, {
    .sync_advan_trans_inputs("template_advan", "template_trans")
  }, ignoreInit = TRUE)

  observeEvent(input$project_tpl_advan, {
    .sync_advan_trans_inputs("project_tpl_advan", "project_tpl_trans")
  }, ignoreInit = TRUE)

  observeEvent(input$create_project, {
    .show_new_project_modal()
  })

  observeEvent(input$confirm_create_project, {
    raw_name <- trimws(input$new_project_name)
    if (!nzchar(raw_name)) {
      showNotification("Enter a project name.", type = "warning")
      return()
    }
    mode <- input$new_project_mode %||% "empty"
    removeModal()
    safe_run("Create project", {
      desc <- trimws(input$new_project_description)
      meta <- nm_workspace_create_project(
        raw_name,
        path = ws_root,
        template = "empty",
        description = if (nzchar(desc)) desc else NULL
      )
      proj <- meta$name
      if (identical(mode, "template")) {
        tryCatch({
          data_source <- input$project_tpl_data_source %||% "synthetic"
          problem <- trimws(input$project_tpl_problem)
          if (!nzchar(problem)) {
            problem <- if (identical(data_source, "synthetic")) "Synthetic demo" else "Template model"
          }
          lbl <- trimws(input$project_tpl_label)
          if (identical(data_source, "synthetic")) {
            syn_id <- input$project_tpl_synthetic_id %||% "theo"
            n_sub <- as.integer(input$project_tpl_theo_nsub %||% 10L)
            if (!is.finite(n_sub) || n_sub < 1L) {
              n_sub <- 10L
            }
            n_sub <- min(500L, n_sub)
            sim <- nm_synthetic_dataset(id = syn_id, n_sub = n_sub, seed = 1L)
            cat <- nm_synthetic_catalog()[[syn_id]]
            csv_name <- cat$csv %||% paste0(syn_id, ".csv")
            tmp_csv <- tempfile(fileext = ".csv")
            on.exit(unlink(tmp_csv), add = TRUE)
            if (requireNamespace("data.table", quietly = TRUE)) {
              data.table::fwrite(sim$data$data, tmp_csv)
            } else {
              write.csv(sim$data$data, tmp_csv, row.names = FALSE)
            }
            data_file <- nm_workspace_import_dataset(proj, tmp_csv, csv_name, root = ws_root)
            ctl <- LibeRation:::.nm_model_to_ctl(
              sim$model,
              data_file = data_file,
              prob = problem,
              method = NULL
            )
          } else {
            up <- input$project_tpl_dataset_file
            if (is.null(up) || is.null(up$datapath) || length(up$datapath) == 0L) {
              stop("Select a dataset file to upload.")
            }
            data_file <- nm_workspace_import_dataset(
              proj, up$datapath[[1L]], up$name[[1L]], root = ws_root
            )
            advan <- as.integer(input$project_tpl_advan)
            trans <- nm_ctl_effective_trans(advan, input$project_tpl_trans)
            if (!nm_ctl_is_valid_pair(advan, trans)) {
              trans <- nm_ctl_effective_trans(advan, nm_ctl_default_trans(advan))
            }
            ode_n <- as.integer(input$project_tpl_ode_ncomp %||% 2L)
            parts <- nm_ctl_template(
              advan = advan,
              trans = trans,
              data_file = data_file,
              problem = problem,
              ode_ncomp = ode_n
            )
            ctl <- nm_ctl_compose(parts)
          }
          nm_workspace_new_version(
            proj,
            root = ws_root,
            template_ctl = ctl,
            data_file = data_file,
            label = if (nzchar(lbl)) lbl else NULL
          )
        }, error = function(e) {
          nm_workspace_delete_project(proj, root = ws_root)
          stop(conditionMessage(e))
        })
      }
      state$projects_rev <- state$projects_rev + 1L
      state$versions_rev <- state$versions_rev + 1L
      select_project(proj)
      msg <- paste("Created project:", proj)
      if (!identical(raw_name, proj)) {
        msg <- paste0(msg, " (name adjusted from \"", raw_name, "\")")
      }
      showNotification(msg, type = "message")
    })
  })

  observeEvent(input$confirm_template_version, {
    if (!.require_clean_for_new_version()) {
      return()
    }
    data_file <- input$template_dataset
    if (is.null(data_file) || !nzchar(data_file)) {
      showNotification("Select a dataset for the new version.", type = "warning")
      return()
    }
    advan <- as.integer(input$template_advan)
    trans <- nm_ctl_effective_trans(advan, input$template_trans)
    if (!nm_ctl_is_valid_pair(advan, trans)) {
      trans <- nm_ctl_effective_trans(advan, nm_ctl_default_trans(advan))
    }
    removeModal()
    safe_run("Template version", {
      parts <- nm_ctl_template(
        advan = advan,
        trans = trans,
        data_file = data_file,
        problem = trimws(input$template_problem),
        ode_ncomp = as.integer(input$template_ode_ncomp %||% 2L)
      )
      lbl <- trimws(input$template_label)
      version_id <- nm_workspace_new_version(
        state$selected_project,
        root = ws_root,
        template_ctl = nm_ctl_compose(parts),
        data_file = data_file,
        label = if (nzchar(lbl)) lbl else NULL
      )
      .finish_new_version(version_id)
    })
  })

  output$model_version_header <- renderUI({
    state$selected_version_id
    state$versions_rev
    state$selected_project
    ver <- state$selected_version_id
    if (!.valid_id(ver)) {
      return(tags$span(class = "model-version-title", "Model version"))
    }
    proj <- state$selected_project
    meta <- if (.valid_id(proj)) .version_meta(proj, ver) else list()
    title <- .shiny_meta_label(meta$label, fallback = ver)
    tags$span(
      class = "model-version-title",
      title,
      if (!identical(title, ver)) {
        tags$span(class = "model-version-id", ver)
      }
    )
  })

  output$model_dirty_banner <- renderUI({
    if (!.valid_id(state$selected_version_id)) {
      return(NULL)
    }
    tags$div(
      style = "padding: 4px 0 6px;",
      if (isTRUE(state$dirty)) {
        tags$span(class = "status-dirty", "● Unsaved changes")
      } else {
        tags$span(class = "status-saved", "● Saved")
      },
      tags$span(
        style = "color: #666; font-size: 11px; margin-left: 8px;",
        paste0(
          state$selected_project, " / ", state$selected_version_id,
          if (.valid_id(state$selected_est_run_id)) {
            paste0(" → ", state$selected_est_run_id)
          } else {
            ""
          }
        )
      )
    )
  })

  current_fit <- reactive({
    state$fit
  })

  output$param_summary <- renderUI({
    req(state$model)
    fit <- current_fit()
    if (is.null(fit)) {
      tags$p(class = "text-muted", style = "font-size: 12px;", "Initial values only.")
    } else {
      tagList(
        tags$p(
          style = "font-size: 12px;",
          strong("Method:"), fit$method,
          tags$br(),
          strong(.shiny_fit_metric_label(fit), ":"), .shiny_fit_metric_value(fit)
        ),
        if (!is.null(fit$bootstrap)) {
          tags$p(
            class = "text-muted",
            style = "font-size: 11px; margin-top: 4px;",
            tags$span(class = "run-flag-badge run-flag-bootstrap", "Bootstrap"),
            " ",
            fit$bootstrap$n_ok, " / ", fit$bootstrap$n_boot %||% "?", " successful replicates."
          )
        }
      )
    }
  })

  output$param_table <- DT::renderDataTable({
    req(state$model)
    fit <- current_fit()
    pt <- if (!is.null(fit)) {
      .shiny_param_table_for_fit(state$model, fit)
    } else {
      nm_workspace_param_table(state$model, fit, compute_se = FALSE)
    }
    display <- .shiny_param_table_display(pt, fit = fit, include_gradient = FALSE)
    display <- .shiny_append_bootstrap_se(display, fit)
    DT::datatable(
      display,
      options = list(dom = "t", pageLength = 30L, scrollY = "320px"),
      rownames = FALSE,
      class = "compact right-param-table",
      escape = FALSE
    )
  })

  output$data_summary <- renderUI({
    df <- .shiny_explore_df()
    if (is.null(df) || nrow(df) == 0L) {
      return(tags$p(class = "text-muted", "Load a project dataset or model version to explore data."))
    }
    tags$p(
      strong(nrow(df)), " observation rows",
      if (!is.null(state$selected_project)) {
        tags$span(" — project: ", strong(state$selected_project))
      }
    )
  })

  .shiny_explore_df <- reactive({
    input$data_explore_dataset
    input$data_explore_all_rows
    state$data
    state$selected_project
    state$explore_data_rev
    all_rows <- isTRUE(input$data_explore_all_rows)
    if (!is.null(state$explore_df) && is.data.frame(state$explore_df) && nrow(state$explore_df) > 0L) {
      df <- .shiny_as_explore_df(state$explore_df)
    } else {
      dat <- .shiny_dataset_table(state$data)
      df <- if (is.null(dat)) NULL else .shiny_as_explore_df(dat)
    }
    if (is.null(df)) {
      return(NULL)
    }
    if (!all_rows && all(c("MDV", "EVID") %in% names(df))) {
      obs <- df[df$MDV == 0L & df$EVID == 0L, , drop = FALSE]
      if (nrow(obs) > 0L) {
        return(obs)
      }
    }
    df
  })

  output$data_dataset_picker <- renderUI({
    state$projects_rev
    state$selected_project
    state$simulations_rev
    state$est_runs_rev
    proj <- state$selected_project
    if (!.valid_id(proj)) {
      return(tags$p(class = "text-muted", style = "font-size: 11px;", "Select a project to browse datasets."))
    }
    choices <- .shiny_explore_dataset_choices(proj, root = ws_root)
    if (length(choices) <= 1L) {
      return(tags$p(
        class = "text-muted", style = "font-size: 11px;",
        "No datasets, estimation outputs, or simulations in this project yet."
      ))
    }
    sel <- state$auto_data_est_key
    if (is.null(sel) || !nzchar(sel) || !(sel %in% choices)) {
      sel <- state$auto_data_sim_key
    }
    if (is.null(sel) || !nzchar(sel) || !(sel %in% choices)) {
      sel <- ""
    }
    selectInput("data_explore_dataset", "Dataset", choices = choices, selected = sel, width = "100%")
  })

  observeEvent(input$data_explore_dataset, {
    path <- input$data_explore_dataset
    if (is.null(path) || !nzchar(path)) {
      state$explore_df <- NULL
      state$explore_data_rev <- (state$explore_data_rev %||% 0L) + 1L
      return()
    }
    proj <- state$selected_project
    if (!.valid_id(proj)) {
      return()
    }
    if (grepl("^__model__/", path)) {
      ver <- sub("^__model__/", "", path)
      safe_run("Load model dataset", {
        parsed <- nm_workspace_parse_model(proj, ver, root = ws_root)
        tbl <- .shiny_dataset_table(parsed$data)
        if (is.null(tbl)) {
          stop("Could not load model input dataset for ", ver, ".")
        }
        state$explore_df <- .shiny_as_explore_df(tbl)
        state$explore_data_rev <- (state$explore_data_rev %||% 0L) + 1L
      })
      return()
    }
    if (grepl("^__est__/", path)) {
      rest <- sub("^__est__/", "", path)
      parts <- strsplit(rest, "/", fixed = TRUE)[[1L]]
      if (length(parts) < 2L) {
        showNotification("Invalid estimation dataset key.", type = "warning")
        return()
      }
      safe_run("Load estimation output", {
        fit <- nm_workspace_load_run_fit(proj, parts[[1L]], parts[[2L]], root = ws_root)
        if (is.null(fit)) {
          stop("Could not load estimation run ", parts[[2L]], ".")
        }
        if (is.null(fit$model)) {
          parsed <- nm_workspace_parse_model(proj, parts[[1L]], root = ws_root)
          if (!is.null(parsed$model)) {
            fit$model <- parsed$model
          }
          if (is.null(fit$data) && !is.null(parsed$data)) {
            fit$data <- parsed$data
          }
        }
        fit <- .shiny_ensure_fit_eta(fit)
        tbl <- .shiny_fit_prediction_table(fit)
        if (is.null(tbl)) {
          stop("Could not build prediction table for ", parts[[2L]], ".")
        }
        state$explore_df <- .shiny_as_explore_df(tbl)
        state$explore_data_rev <- (state$explore_data_rev %||% 0L) + 1L
      })
      return()
    }
    if (grepl("^__sim__/", path)) {
      rest <- sub("^__sim__/", "", path)
      parts <- strsplit(rest, "/", fixed = TRUE)[[1L]]
      if (length(parts) < 2L) {
        showNotification("Invalid simulation dataset key.", type = "warning")
        return()
      }
      safe_run("Load simulation dataset", {
        sim_obj <- nm_workspace_load_sim(proj, parts[[1L]], parts[[2L]], root = ws_root)
        tbl <- .shiny_sim_primary_df(sim_obj)
        if (is.null(tbl)) {
          stop("Could not load simulation ", parts[[2L]], ".")
        }
        state$explore_df <- .shiny_as_explore_df(tbl)
        state$explore_data_rev <- (state$explore_data_rev %||% 0L) + 1L
      })
      return()
    }
    fp <- .dataset_full_path(proj, path)
    if (is.null(fp)) {
      showNotification("Dataset file not found.", type = "warning")
      return()
    }
    safe_run("Load dataset", {
      tbl <- if (requireNamespace("data.table", quietly = TRUE)) {
        as.data.frame(data.table::fread(fp, showProgress = FALSE))
      } else {
        utils::read.csv(fp, check.names = FALSE)
      }
      state$explore_df <- .shiny_as_explore_df(tbl)
      state$explore_data_rev <- (state$explore_data_rev %||% 0L) + 1L
    })
  }, ignoreInit = TRUE)

  observe({
    df <- tryCatch(.shiny_explore_df(), error = function(e) {
      push_log("error", conditionMessage(e), source = "Data explore columns")
      NULL
    })
    cols <- if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) {
      character()
    } else {
      names(df)
    }
    prev_x <- input$data_explore_x
    prev_y <- input$data_explore_y
    prev_s <- input$data_explore_strat
    prev_split <- input$data_explore_split
    x_sel <- if (!is.null(prev_x) && prev_x %in% cols) prev_x else {
      if ("TIME" %in% cols) "TIME" else if (length(cols)) cols[[1L]] else ""
    }
    y_sel <- if (!is.null(prev_y) && prev_y %in% cols) prev_y else {
      if ("DV" %in% cols) "DV" else if (length(cols) > 1L) cols[[2L]] else x_sel
    }
    strat_choices <- c("(none)" = "", stats::setNames(cols, cols))
    s_sel <- if (!is.null(prev_s) && prev_s %in% strat_choices) {
      prev_s
    } else if (!is.null(prev_s) && prev_s %in% cols) {
      prev_s
    } else {
      ""
    }
    split_sel <- if (!is.null(prev_split) && prev_split %in% strat_choices) {
      prev_split
    } else if (!is.null(prev_split) && prev_split %in% cols) {
      prev_split
    } else {
      ""
    }
    updateSelectInput(session, "data_explore_x", choices = cols, selected = x_sel)
    updateSelectInput(session, "data_explore_y", choices = cols, selected = y_sel)
    updateSelectInput(session, "data_explore_strat", choices = strat_choices, selected = s_sel)
    updateSelectInput(session, "data_explore_split", choices = strat_choices, selected = split_sel)
  })

  hideTab("center_tabs", "GOF", session = session)
  hideTab("center_tabs", "VPC", session = session)
  hideTab("center_tabs", "NPC", session = session)
  hideTab("center_tabs", "NPDE", session = session)

  observe({
    state$gof_diag_rev
    state$est_runs_rev
    state$simulations_rev
    .shiny_sync_center_tabs(session, state$fit, state$sim_vpc_data)
  })

  output$vpc_tab_summary <- renderUI({
    if (is.null(state$sim_vpc_data)) {
      return(tags$p(class = "text-muted", style = "font-size: 12px;",
                    "Run a VPC simulation from a fitted run to view the predictive check."))
    }
    run_id <- state$selected_est_run_id %||% ""
    tags$p(
      style = "font-size: 12px;",
      if (nzchar(run_id)) tagList(strong("Estimation run:"), run_id, tags$br()) else NULL,
      "Time is auto-binned (deciles). ",
      "Blue line: observed median; red dashed: observed 10th/90th percentiles. ",
      "Blue shaded bands: 90% simulation CI of predicted 10th/90th; ",
      "red band: 90% simulation CI of predicted median."
    )
  })

  output$center_vpc_plot <- renderPlot({
    state$sim_vpc_data
    state$sim_vpc_obs
    state$selected_est_run_id
    input$vpc_pc_correct
    fit <- current_fit()
    pc_ok <- isTRUE(input$vpc_pc_correct)
    if (pc_ok && is.null(fit)) {
      showNotification(
        "Prediction-corrected VPC requires the linked estimation fit (IPRED).",
        type = "warning", duration = 6L
      )
      pc_ok <- FALSE
    }
    .shiny_vpc_plot(
      state$sim_vpc_data,
      obs = state$sim_vpc_obs,
      fit = fit,
      pc_correct = pc_ok,
      title = if (pc_ok) {
        "Prediction-corrected visual predictive check"
      } else {
        "Visual predictive check"
      }
    )
  })

  output$data_explore_plot <- renderPlot({
    tryCatch({
      input$data_explore_x
      input$data_explore_y
      input$data_explore_strat
      input$data_explore_split
      input$data_explore_all_rows
      input$data_explore_bin_x_pos
      input$data_explore_bin_y_pos
      input$data_explore_bin_x
      input$data_explore_bin_y
      input$data_explore_bin_x_manual
      input$data_explore_bin_y_manual
      input$data_explore_bin_x_breaks
      input$data_explore_bin_y_breaks
      input$data_explore_bins
      input$data_explore_plot_type
      input$data_explore_show_points
      input$data_explore_point_scatter
      input$data_explore_point_size
      input$data_explore_point_shape
      input$data_explore_line_mode
      input$data_explore_q_interval
      input$data_explore_shade_alpha
      input$data_explore_adjust_plot
      input$data_explore_title
      input$data_explore_xlab
      input$data_explore_ylab
      df <- .shiny_explore_df()
      if (is.null(df)) {
        .shiny_empty_plot("Select variables to plot")
        return(invisible(NULL))
      }
      x_breaks <- if (isTRUE(input$data_explore_bin_x) &&
                      isTRUE(input$data_explore_bin_x_manual)) {
        .shiny_parse_breaks(input$data_explore_bin_x_breaks)
      } else {
        NULL
      }
      y_breaks <- if (isTRUE(input$data_explore_bin_y) &&
                      isTRUE(input$data_explore_bin_y_manual)) {
        .shiny_parse_breaks(input$data_explore_bin_y_breaks)
      } else {
        NULL
      }
      p <- .shiny_data_explore_plot(
        df,
        xvar = input$data_explore_x,
        yvar = input$data_explore_y,
        strat = input$data_explore_strat,
        split_by = input$data_explore_split,
        bin_x = isTRUE(input$data_explore_bin_x),
        bin_y = isTRUE(input$data_explore_bin_y),
        bin_x_pos = input$data_explore_bin_x_pos %||% "equal",
        bin_y_pos = input$data_explore_bin_y_pos %||% "equal",
        nbins = as.integer(input$data_explore_bins %||% 10L),
        bin_x_breaks = x_breaks,
        bin_y_breaks = y_breaks,
        plot_type = input$data_explore_plot_type %||% "points",
        show_points = isTRUE(input$data_explore_show_points),
        q_interval = as.numeric(input$data_explore_q_interval %||% 95),
        shade_alpha = as.numeric(input$data_explore_shade_alpha %||% 25) / 100,
        plot_title = if (isTRUE(input$data_explore_adjust_plot)) input$data_explore_title else NULL,
        xlab_custom = if (isTRUE(input$data_explore_adjust_plot)) input$data_explore_xlab else NULL,
        ylab_custom = if (isTRUE(input$data_explore_adjust_plot)) input$data_explore_ylab else NULL,
        point_scatter = as.numeric(input$data_explore_point_scatter %||% 25) / 100,
        point_size = as.numeric(input$data_explore_point_size %||% 0.85),
        point_shape = as.integer(input$data_explore_point_shape %||% 16L),
        line_mode = input$data_explore_line_mode %||% "individual"
      )
      if (is.null(p)) {
        .shiny_empty_plot("Nothing to plot for selected variables")
        return(invisible(NULL))
      }
      p
    }, error = function(e) {
      push_log("error", conditionMessage(e), source = "Data exploration plot")
      .shiny_empty_plot("Plot failed — see log")
      invisible(NULL)
    })
  })

  output$data_table <- DT::renderDataTable({
    input$data_show_dataset
    input$data_explore_all_rows
    if (!isTRUE(input$data_show_dataset)) {
      return(DT::datatable(
        data.frame(note = "Enable \"Show dataset\" to view rows."),
        options = list(dom = "t"), rownames = FALSE, selection = "none"
      ))
    }
    df <- .shiny_explore_df()
    if (is.null(df) || nrow(df) == 0L) {
      return(DT::datatable(
        data.frame(note = "No data loaded."),
        options = list(dom = "t"), rownames = FALSE, selection = "none"
      ))
    }
    df <- as.data.frame(df)
    DT::datatable(
      df,
      options = list(
        pageLength = 25L,
        lengthMenu = list(c(10, 25, 50, 100, 500, 1000, -1), c(10, 25, 50, 100, 500, 1000, "All")),
        scrollX = TRUE,
        dom = "lftip"
      ),
      rownames = FALSE,
      class = "compact"
    )
  })

  build_est_args <- function() {
    infer_on <- isTRUE(input$compute_inference) && (input$method %||% "") != "BAYES"
    infer_meth <- input$inference_method %||% "cov_step"
    infer_hess <- "numeric"
    cov_method <- "auto"
    if (infer_on) {
      if (infer_meth == "cov_step") {
        cov_method <- "auto"
        infer_hess <- "numeric"
      } else if (infer_meth == "hessian_auto") {
        cov_method <- "hessian"
        infer_hess <- "auto"
      } else if (infer_meth == "hessian_ad") {
        cov_method <- "hessian"
        infer_hess <- "ad"
      } else {
        cov_method <- "hessian"
        infer_hess <- "numeric"
      }
    }
    ctl <- list(
      maxit = as.integer(input$maxit),
      n_cores = 1L,
      compute_inference = infer_on,
      cov_method = cov_method,
      infer_hessian = infer_hess,
      min_retries = as.integer(input$min_retries %||% 0L),
      tweak_inits = isTRUE(input$tweak_inits)
    )
    args <- list(
      grad = input$grad,
      pk_engine = input$pk_engine,
      control = ctl
    )
    meth <- input$method
    if (meth %in% c("FOCE", "FOCEI", "LAPLACE") && !is.null(state$model)) {
      args$pk_engine <- "cpp"
      if (LibeRation:::.nm_cpp_capable(state$model)) {
        args$grad <- "cpp"
      }
    }
    if (meth %in% c("FOCE", "FOCEI")) {
      args$max_outer <- as.integer(input$max_outer)
    }
    if (meth == "SAEM") {
      args$n_iter <- as.integer(input$n_iter)
      args$n_burn <- as.integer(input$n_burn_saem %||% 10L)
      args$n_mcmc <- as.integer(input$n_mcmc %||% 1L)
      args$engine <- "cpp"
    }
    if (meth == "LAPLACE") {
      args$engine <- "cpp"
      args$grad <- "cpp"
      args$n_quad <- as.integer(input$n_quad %||% 5L)
    }
    if (meth == "BAYES") {
      args$engine <- "cpp"
      args$sampler <- input$sampler
      args$n_burn <- as.integer(input$n_burn)
      args$n_sample <- as.integer(input$n_sample)
      args$n_thin <- 1L
    }
    if (meth == "IMP") {
      args$n_imp <- as.integer(input$n_imp %||% 50L)
      args$n_quad <- as.integer(input$imp_n_quad %||% 5L)
    }
    if (isTRUE(input$est_bootstrap) && meth != "BAYES") {
      args$bootstrap_n <- as.integer(input$bootstrap_n %||% 30L)
      args$bootstrap_seed <- as.integer(input$bootstrap_seed %||% 1L)
      args$bootstrap_control <- list(
        maxit = as.integer(input$maxit),
        n_cores = 1L,
        compute_inference = FALSE
      )
    } else {
      args$bootstrap_n <- 0L
    }
    args
  }

  observeEvent(input$submit_job, {
    if (isTRUE(state$dirty)) {
      showNotification("Save the model before submitting a job.", type = "warning")
      return()
    }
    req(state$model, state$data, state$selected_project, state$selected_version_id)
    if (!requireNamespace("callr", quietly = TRUE)) {
      showNotification("Install package 'callr' to submit jobs.", type = "error")
      return()
    }
    v_ok <- tryCatch(
      nm_validate_model(state$model, state$data, stop_on_error = FALSE),
      error = function(e) list(ok = FALSE, issues = conditionMessage(e))
    )
    if (!isTRUE(v_ok$ok)) {
      showNotification(
        paste(v_ok$issues, collapse = "\n"),
        type = "error",
        duration = NULL
      )
      return()
    }
    label <- input$job_label
    if (!nzchar(trimws(label))) {
      label <- paste(state$selected_version_id, input$method, format(Sys.time(), "%H:%M:%S"))
    }
    est_run_id <- nm_workspace_new_run_id(
      state$selected_project,
      state$selected_version_id,
      root = ws_root
    )
    extra <- build_est_args()
    job <- do.call(
      nm_job_submit,
      c(
        list(
          model = state$model,
          data = state$data,
          method = input$method,
          label = label,
          job_root = job_root
        ),
        extra
      )
    )
    state$handles[[job$id]] <- job$process
    state$selected_job <- job$id
    state$active_job_type <- "est"
    state$active_job_sim_id <- NULL
    state$active_job_version <- state$selected_version_id
    state$active_job_est_run <- est_run_id
    state$active_job_project <- state$selected_project
    state$jobs_rev <- state$jobs_rev + 1L
    message(
      "[LibeRation] Estimation started: job=", job$id,
      "  method=", input$method,
      "  run_id=", est_run_id
    )
    showNotification(paste("Job submitted:", job$id), type = "message", duration = 5L)
    if (isTRUE(input$est_bootstrap) && input$method != "BAYES") {
      showNotification(
        paste("Bootstrap (", input$bootstrap_n, " replicates) will run after estimation in the worker.", sep = ""),
        type = "message",
        duration = 8L
      )
    }
    removeModal()
  })

  observe({
    input$submit_job
    input$refresh_jobs
    input$ribbon_tab
    state$jobs_rev
    df <- nm_job_list(job_root)
    sig <- nm_job_watch_signature(job_root)
    if (!identical(sig, isolate(state$job_watch_last))) {
      state$job_watch_last <- sig
      state$jobs_rev <- isolate(state$jobs_rev) + 1L
    }
    n_active <- nrow(df) > 0L && any(df$status %in% c("queued", "running"))
    on_jobs_tab <- identical(input$ribbon_tab, "jobs")
    if (n_active) {
      invalidateLater(150L, session)
    } else if (on_jobs_tab && nrow(df) > 0L) {
      invalidateLater(500L, session)
    }
  })

  job_elapsed_timer <- reactive({
    df <- nm_job_list(job_root)
    if (nrow(df) > 0L && any(df$status %in% c("queued", "running"))) {
      invalidateLater(1000L, session)
    }
    Sys.time()
  })

  jobs_df <- reactive({
    state$jobs_rev
    input$submit_job
    input$refresh_jobs
    input$cancel_job
    input$cleanup_jobs
    job_elapsed_timer()
    nm_job_list(job_root)
  })

  observeEvent(jobs_df(), {
    df <- jobs_df()
    cache <- state$job_status_cache
    if (nrow(df) > 0L) {
      for (i in seq_len(nrow(df))) {
        id <- df$id[[i]]
        st <- df$status[[i]]
        prev <- cache[[id]]
        if (is.null(prev) || !identical(prev, st)) {
          cache[[id]] <- st
          st_full <- nm_job_status(id, job_root)
          jtype <- if (is.null(st_full$job_type)) "est" else st_full$job_type
          if (!is.null(prev) && !identical(prev, st)) {
            if (identical(st, "success")) {
              if (identical(jtype, "sim")) {
                message(
                  "[LibeRation] Simulation finished: job=", id,
                  "  sim_id=", st_full$sim_id %||% ""
                )
              } else {
                obj <- st_full$objective
                if (!is.null(obj) && is.finite(obj)) {
                  message(
                    "[LibeRation] Estimation finished: job=", id,
                    "  objective=", round(obj, 4)
                  )
                } else {
                  message("[LibeRation] Estimation finished: job=", id)
                }
              }
            } else if (identical(st, "error")) {
              message("[LibeRation] Job failed: id=", id, " (", jtype, ")")
            } else if (identical(st, "running") && !identical(prev, "running")) {
              if (identical(jtype, "sim")) {
                message("[LibeRation] Simulation running: job=", id)
              } else {
                message("[LibeRation] Estimation running: job=", id)
              }
            }
          }
          if (identical(st, "success") &&
              identical(jtype, "est") &&
              identical(id, state$selected_job) &&
              !is.null(state$active_job_project) &&
              !is.null(state$active_job_version) &&
              !is.null(state$active_job_est_run)) {
            fit <- tryCatch(nm_job_result(id, job_root), error = function(e) NULL)
            if (!is.null(fit)) {
              nm_workspace_save_run(
                state$active_job_project,
                state$active_job_version,
                state$active_job_est_run,
                fit,
                root = ws_root,
                label = fit$method,
                job_id = id
              )
              if (identical(state$selected_project, state$active_job_project) &&
                  identical(state$selected_version_id, state$active_job_version)) {
                apply_est_run_fit(
                  state$active_job_project,
                  state$active_job_version,
                  state$active_job_est_run,
                  fit
                )
              }
              state$est_runs_rev <- state$est_runs_rev + 1L
              state$versions_rev <- state$versions_rev + 1L
              if (!is.null(fit$bootstrap)) {
                showNotification(
                  paste(
                    "Bootstrap complete:",
                    fit$bootstrap$n_ok, "/", fit$bootstrap$n_boot,
                    "replicates. See Parameters tab for Bootstrap SE."
                  ),
                  type = "message",
                  duration = 8L
                )
              }
            }
          }
          if (identical(st, "success") &&
              identical(jtype, "sim") &&
              identical(id, state$selected_job) &&
              !is.null(state$active_job_project)) {
            state$simulations_rev <- state$simulations_rev + 1L
            ver <- st_full$version_id
            sim_id <- st_full$sim_id
            diag_run <- st_full$diag_est_run %||% ""
            fit_diag <- NULL
            if (nzchar(diag_run)) {
              fit_diag <- tryCatch(
                nm_workspace_load_run_fit(
                  state$active_job_project,
                  ver %||% state$active_job_version,
                  diag_run,
                  root = ws_root
                ),
                error = function(e) NULL
              )
              if (!is.null(fit_diag)) {
                parts <- c(
                  if (.shiny_fit_has_npc(fit_diag)) "NPC",
                  if (.shiny_fit_has_npde(fit_diag)) "NPDE"
                )
                if (length(parts) > 0L) {
                  showNotification(
                    paste(
                      paste(parts, collapse = "/"), "complete — see",
                      paste(parts, collapse = " / "), "tab(s)."
                    ),
                    type = "message",
                    duration = 8L
                  )
                }
                state$versions_rev <- state$versions_rev + 1L
                state$est_runs_rev <- state$est_runs_rev + 1L
              }
            }
            if (!is.null(ver) && nzchar(ver) && !(ver %in% state$expanded_versions)) {
              state$expanded_versions <- c(state$expanded_versions, ver)
            }
            if (identical(state$selected_project, state$active_job_project) &&
                !is.null(ver) && nzchar(ver)) {
              est_run <- state$active_job_est_run
              diag_only <- isTRUE(state$active_job_diag_only)
              if (.valid_id(est_run) && (diag_only || nzchar(diag_run))) {
                if (!identical(ver, state$selected_version_id)) {
                  load_model_version(state$active_job_project, ver, keep_est_run = TRUE)
                }
                if (!is.null(fit_diag) && identical(est_run, diag_run)) {
                  apply_est_run_fit(
                    state$active_job_project, ver, est_run, fit_diag
                  )
                } else {
                  load_est_run(state$active_job_project, ver, est_run)
                }
              } else if (!diag_only && !is.null(sim_id) && nzchar(sim_id) &&
                         !grepl("^_", sim_id)) {
                select_simulation(ver, sim_id, force_reload = TRUE)
              } else if (.valid_id(est_run)) {
                if (!identical(ver, state$selected_version_id)) {
                  load_model_version(state$active_job_project, ver, keep_est_run = TRUE)
                }
                load_est_run(state$active_job_project, ver, est_run)
              }
            }
            state$active_job_diag_only <- FALSE
          }
        }
      }
    }
    state$job_status_cache <- cache
  }, ignoreInit = TRUE)

  observe({
    if (length(state$handles) == 0L) {
      return()
    }
    invalidateLater(2000L, session)
    for (id in names(state$handles)) {
      proc <- state$handles[[id]]
      if (is.null(proc)) next
      if (!proc$is_alive()) {
        state$handles[[id]] <- NULL
      }
    }
  })

  observeEvent(input$cancel_job, {
    req(state$selected_job)
    st <- nm_job_status(state$selected_job, job_root)
    if (is.null(st)) {
      showNotification("No job selected.", type = "warning")
      return()
    }
    if (st$status %in% c("success", "error", "cancelled")) {
      showNotification(paste("Job already finished:", st$status), type = "message")
      return()
    }
    proc <- state$handles[[state$selected_job]]
    if (!is.null(proc)) {
      tryCatch(proc$kill(), error = function(e) NULL)
      state$handles[[state$selected_job]] <- NULL
    }
    nm_job_cancel(state$selected_job, job_root)
    showNotification(paste("Cancelled job:", state$selected_job), type = "warning", duration = 4L)
  })

  output$jobs_tree <- renderUI({
    jobs_df()
    state$selected_job
    state$expanded_jobs
    .shiny_jobs_tree_ui(jobs_df(), state$selected_job, state$expanded_jobs, job_root)
  })

  output$project_tpl_synthetic_desc <- renderUI({
    id <- input$project_tpl_synthetic_id %||% "theo"
    tags$p(
      class = "text-muted", style = "font-size: 11px;",
      .shiny_synthetic_tpl_description(id)
    )
  })

  output$jobs_refresh_clock <- renderText({
    df <- jobs_df()
    n_active <- sum(df$status %in% c("queued", "running"))
    elapsed_txt <- if (n_active > 0L && nrow(df) > 0L) {
      active <- df[df$status %in% c("queued", "running"), , drop = FALSE]
      el <- vapply(seq_len(nrow(active)), function(i) {
        .shiny_job_duration_label(active$started[[i]], active$finished[[i]], active$status[[i]])
      }, character(1L))
      el <- el[nzchar(el)]
      if (length(el) > 0L) paste0(" | longest ", el[[1L]]) else ""
    } else {
      ""
    }
    paste(
      "Updated:", format(Sys.time(), "%H:%M:%S"),
      "|", nrow(df), "job(s)",
      if (n_active > 0L) paste0("(", n_active, " active", elapsed_txt, ")") else ""
    )
  })

  observeEvent(input$refresh_jobs, {
    state$job_watch_last <- ""
    state$jobs_rev <- isolate(state$jobs_rev) + 1L
  })

  observeEvent(input$cleanup_jobs, {
    n <- nm_job_cleanup(job_root)
    showNotification(paste("Removed", n, "finished job(s)."), type = "message")
  })

  output$job_status_banner <- renderUI({
    jobs_df()
    req(state$selected_job)
    st <- nm_job_status(state$selected_job, job_root)
    req(st)
    cls <- switch(
      st$status,
      success = "alert-success",
      error = "alert-danger",
      cancelled = "alert-warning",
      running = "alert-info",
      queued = "alert-secondary",
      "alert-light"
    )
    err <- .shiny_clean_job_text(st$error)
    if (identical(st$status, "error") && !nzchar(err)) {
      err_path <- file.path(job_root, state$selected_job, "error.txt")
      if (file.exists(err_path)) {
        err <- paste(readLines(err_path, warn = FALSE), collapse = "\n")
      }
    }
    jtype <- if (is.null(st$job_type)) "est" else st$job_type
    jtype_label <- if (identical(jtype, "sim")) "simulation" else "estimation"
    finished_lbl <- .shiny_format_job_time(st$finished)
    duration_lbl <- .shiny_job_duration_label(st$started, st$finished, st$status)
    tags$div(
      class = paste("alert", cls),
      role = "alert",
      tags$p(
        strong("Job:"), state$selected_job,
        " | ", strong("Type:"), jtype_label,
        " | ", strong("Status:"), st$status,
        if (nzchar(duration_lbl)) {
          tags$span(" | ", strong(if (identical(st$status, "running")) "Elapsed:" else "Duration:"), duration_lbl)
        },
        if (identical(jtype, "sim")) {
          sim_lbl <- .shiny_clean_job_text(st$sim_id)
          nrep <- if (is.null(st$n_sim)) "" else as.character(st$n_sim)
          tagList(
            if (nzchar(sim_lbl)) {
              tags$span(" | ", strong("Sim:"), sim_lbl)
            },
            if (nzchar(nrep)) {
              tags$span(" | ", strong("Replicates:"), nrep)
            }
          )
        } else if (!is.null(st$objective) && is.finite(st$objective)) {
          tags$span(" | ", strong("Objective:"), round(st$objective, 4))
        },
        if (nzchar(finished_lbl)) {
          tags$span(" | ", strong("Finished:"), finished_lbl)
        }
      ),
      if (nzchar(err)) tags$pre(style = "white-space: pre-wrap;", err)
    )
  })

  output$job_log <- renderText({
    jobs_df()
    req(state$selected_job)
    nm_job_log(state$selected_job, tail = 60L, job_root = job_root)
  })

  gof_data <- reactive({
    state$gof_diag_rev
    fit <- current_fit()
    req(fit)
    safe_run("GOF data", {
      if (is.null(fit$model) && !is.null(state$model)) {
        fit$model <- state$model
      }
      if (is.null(fit$data) && !is.null(state$data)) {
        fit$data <- state$data
      }
      fit <- .shiny_ensure_fit_eta(fit)
      ind <- predict(fit, type = "ipred")
      pop <- predict(fit, type = "ppred")
      obs <- .shiny_gof_obs(ind)
      if (!is.null(fit$gof)) {
        gof <- as.data.frame(fit$gof)
        gof <- gof[gof$MDV == 0L & gof$EVID == 0L, , drop = FALSE]
        extra_cols <- setdiff(
          names(gof),
          c(names(obs), "IPRED", "PRED", "RES", "WRES", "IWRES", "CPRED", "CRES")
        )
        if (length(extra_cols) > 0L && all(c("ID", "TIME") %in% names(gof))) {
          obs <- merge(
            obs,
            gof[, c("ID", "TIME", extra_cols), drop = FALSE],
            by = c("ID", "TIME"),
            all.x = TRUE,
            sort = FALSE
          )
        }
      }
      list(fit = fit, ind = ind, pop = pop, obs = obs)
    })
  })

  output$gof_diag_status <- renderUI({
    state$gof_diag_rev
    fit <- current_fit()
    if (is.null(fit)) {
      return(tags$p(class = "text-muted", style = "font-size: 11px;",
                    "Load an estimation run to view GOF plots."))
    }
    msg <- if (!is.null(fit$gof) && "CWRES" %in% names(fit$gof)) {
      "CWRES and basic residual plots."
    } else {
      "Run simulation with NPC or NPDE to compute CWRES, or re-estimate with inference."
    }
    tags$p(class = "text-muted", style = "font-size: 11px; margin-bottom: 8px;", msg)
  })

  output$npc_tab_summary <- renderUI({
    state$gof_diag_rev
    fit <- current_fit()
    if (is.null(fit) || !.shiny_fit_has_npc(fit)) {
      return(tags$p(
        class = "text-muted", style = "font-size: 11px;",
        "No NPC for this run. Use Run simulation with Compute NPC on fit."
      ))
    }
    info <- fit$npc %||% fit$npc_npde
    tags$p(
      style = "font-size: 11px;",
      strong("NPC:"),
      info$n_sim, " simulations (",
      info$n_ok, " observation rows with valid replicates)."
    )
  })

  output$npde_tab_summary <- renderUI({
    state$gof_diag_rev
    fit <- current_fit()
    if (is.null(fit) || !.shiny_fit_has_npde(fit)) {
      return(tags$p(
        class = "text-muted", style = "font-size: 11px;",
        "No NPDE for this run. Use Run simulation with Compute NPDE on fit."
      ))
    }
    info <- fit$npde %||% fit$npc_npde
    tags$p(
      style = "font-size: 11px;",
      strong("NPDE:"),
      info$n_sim, " simulations (",
      info$n_ok, " observation rows with valid replicates)."
    )
  })

  output$fit_summary_compact <- renderUI({
    fit <- current_fit()
    if (is.null(fit)) {
      return(tags$p(class = "text-muted", style = "font-size: 12px;",
                  "No fit — run estimation from the Run tab."))
    }
    tags$p(
      style = "font-size: 12px; margin-bottom: 6px;",
      strong(fit$method), " | objective = ", round(fit$objective, 4)
    )
  })

  output$gof_time_combined <- renderPlot({
    safe_run("GOF time profiles", {
      gd <- gof_data()
      req(gd)
      ind_lines <- .shiny_profile_lines(gd$ind, "IPRED")
      pop_lines <- .shiny_profile_lines(gd$pop, "IPRED")
      obs <- gd$obs
      if (nrow(obs) == 0L) {
        return(NULL)
      }
      if (.shiny_has_ggplot()) {
        ggplot2::ggplot() +
          ggplot2::geom_line(
            data = pop_lines,
            ggplot2::aes(x = TIME, y = IPRED, group = factor(ID)),
            color = .shiny_col_pred, linewidth = 0.7, alpha = 0.9, linetype = "22"
          ) +
          ggplot2::geom_line(
            data = ind_lines,
            ggplot2::aes(x = TIME, y = IPRED, group = factor(ID)),
            color = .shiny_col_ipred, linewidth = 0.7, alpha = 0.9
          ) +
          ggplot2::geom_point(
            data = obs, ggplot2::aes(x = TIME, y = .data[["DV"]]),
            color = .shiny_col_dv, size = 1.8, alpha = 0.85
          ) +
          ggplot2::labs(title = "Time profiles", x = "Time", y = "DV") +
          ggplot2::theme_bw(base_size = 10)
      } else {
        plot(obs$TIME, obs$DV, pch = 16, col = .shiny_col_dv, main = "Time profiles")
      }
    }) %||% .shiny_empty_plot()
  })

  output$gof_dv_ipred <- renderPlot({
    safe_run("GOF DV vs IPRED", {
      gd <- gof_data()
      req(gd)
      .shiny_scatter_gof(gd$obs, "IPRED", "DV vs IPRED", "IPRED")
    }) %||% .shiny_empty_plot()
  })

  output$gof_wres_time <- renderPlot({
    safe_run("GOF WRES vs time", {
      gd <- gof_data()
      req(gd)
      .shiny_wres_vs(gd$obs, "TIME", "WRES vs time", "Time")
    }) %||% .shiny_empty_plot()
  })

  output$gof_cwres_time <- renderPlot({
    safe_run("GOF CWRES vs time", {
      gd <- gof_data()
      req(gd)
      if (!"CWRES" %in% names(gd$obs)) {
        return(.shiny_empty_plot("Run simulation with NPC or NPDE to compute CWRES"))
      }
      .shiny_residual_vs(gd$obs, "CWRES", "CWRES vs time", "CWRES")
    }) %||% .shiny_empty_plot()
  })

  output$gof_qq_cwres <- renderPlot({
    safe_run("GOF CWRES Q-Q", {
      gd <- gof_data()
      req(gd)
      if (!"CWRES" %in% names(gd$obs)) {
        return(.shiny_empty_plot("Run simulation with NPC or NPDE to compute CWRES"))
      }
      .shiny_qq_plot(gd$obs$CWRES, title = "CWRES normal Q-Q plot")
    }) %||% .shiny_empty_plot()
  })

  output$gof_npde_qq <- renderPlot({
    safe_run("GOF NPDE Q-Q", {
      gd <- gof_data()
      req(gd)
      if (!"NPDE" %in% names(gd$obs)) {
        return(.shiny_empty_plot("Enable NPDE in Run simulation"))
      }
      .shiny_qq_plot(gd$obs$NPDE, title = "NPDE normal Q-Q plot")
    }) %||% .shiny_empty_plot()
  })

  output$gof_npde_time <- renderPlot({
    safe_run("GOF NPDE vs time", {
      gd <- gof_data()
      req(gd)
      if (!"NPDE" %in% names(gd$obs)) {
        return(.shiny_empty_plot("Enable NPDE in Run simulation"))
      }
      .shiny_residual_vs(gd$obs, "NPDE", "NPDE vs time", "NPDE")
    }) %||% .shiny_empty_plot()
  })

  output$gof_npc_hist <- renderPlot({
    safe_run("GOF NPC histogram", {
      gd <- gof_data()
      req(gd)
      if (!"NPC" %in% names(gd$obs)) {
        return(.shiny_empty_plot("Enable NPC in Run simulation"))
      }
      x <- gd$obs$NPC[is.finite(gd$obs$NPC)]
      if (length(x) == 0L) {
        return(.shiny_empty_plot("No NPC values"))
      }
      if (.shiny_has_ggplot()) {
        ggplot2::ggplot(data.frame(NPC = x), ggplot2::aes(x = NPC)) +
          ggplot2::geom_histogram(bins = 20, fill = "steelblue", color = "white", alpha = 0.85) +
          ggplot2::geom_vline(xintercept = 0.05, linetype = 2, color = "firebrick") +
          ggplot2::labs(title = "NPC p-values", x = "NPC", y = "Count") +
          ggplot2::theme_bw(base_size = 11)
      } else {
        hist(x, main = "NPC p-values", xlab = "NPC", col = "steelblue", border = "white")
        abline(v = 0.05, lty = 2, col = "firebrick")
      }
    }) %||% .shiny_empty_plot()
  })

  output$gof_npc_time <- renderPlot({
    safe_run("GOF NPC vs time", {
      gd <- gof_data()
      req(gd)
      if (!"NPC" %in% names(gd$obs)) {
        return(.shiny_empty_plot("Enable NPC in Run simulation"))
      }
      obs <- gd$obs
      ok <- is.finite(obs$NPC)
      obs <- obs[ok, , drop = FALSE]
      if (nrow(obs) == 0L) {
        return(.shiny_empty_plot("No NPC values"))
      }
      if (.shiny_has_ggplot()) {
        ggplot2::ggplot(obs, ggplot2::aes(x = TIME, y = NPC)) +
          ggplot2::geom_hline(yintercept = 0.05, linetype = 2, color = "firebrick") +
          ggplot2::geom_point(alpha = 0.65, color = "steelblue", size = 1.5) +
          ggplot2::labs(title = "NPC vs time", x = "Time", y = "NPC") +
          ggplot2::theme_bw(base_size = 11)
      } else {
        plot(obs$TIME, obs$NPC, pch = 16, col = "steelblue", main = "NPC vs time",
             xlab = "Time", ylab = "NPC")
        abline(h = 0.05, lty = 2, col = "firebrick")
      }
    }) %||% .shiny_empty_plot()
  })

  observeEvent(input$generate_report, {
    fit <- current_fit()
    if (is.null(fit)) {
      showNotification("No fit available for report.", type = "warning")
      return()
    }
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      showNotification(
        "ggplot2 not installed — report will use base graphics. Install ggplot2 for nicer plots.",
        type = "message", duration = 6L
      )
    }
    req(state$selected_project, state$selected_version_id)
    name <- gsub("[^A-Za-z0-9._-]", "_", input$report_name)
    if (!nzchar(name)) {
      name <- paste0("report_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    }
    report_dir <- nm_workspace_reports_dir(state$selected_project, root = ws_root)
    dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
    pdf_path <- file.path(report_dir, paste0(name, ".pdf"))
    secs <- input$report_sections
    sections <- list(
      summary = "summary" %in% secs,
      parameters = "parameters" %in% secs,
      gof_time = "gof_time" %in% secs,
      gof_ipred_time = "gof_ipred_time" %in% secs,
      gof_scatter = "gof_scatter" %in% secs,
      gof_residuals = "gof_residuals" %in% secs,
      diag_shrinkage = "diag_shrinkage" %in% secs,
      diag_eta = "diag_eta" %in% secs,
      narrative_stub = "narrative_stub" %in% secs
    )
    tryCatch(
      {
        res <- nm_report_pdf(
          fit,
          pdf_path,
          sections = sections,
          project_meta = list(
            project = state$selected_project,
            version_id = state$selected_version_id,
            run_id = state$selected_est_run_id,
            workspace = ws_root
          )
        )
        state$last_report <- res
        showNotification("Report generated.", type = "message")
      },
      error = function(e) {
        showNotification(conditionMessage(e), type = "error")
      }
    )
  })

  output$report_status <- renderUI({
    lr <- state$last_report
    if (is.null(lr)) {
      return(tags$p(class = "text-muted", style = "font-size: 12px;", "No report yet."))
    }
    tags$div(
      style = "font-size: 12px;",
      tags$p(strong("PDF:"), lr$pdf),
      if (file.exists(lr$manifest)) {
        tags$p(strong("Manifest:"), lr$manifest)
      }
    )
  })

  observe({
    if (isTRUE(state$initialized)) {
      return()
    }
    df <- projects_df()
    if (nrow(df) == 0L) {
      return()
    }
    state$initialized <- TRUE
    project <- df$project[[1L]]
    if (!.valid_id(project)) {
      return()
    }
    state$selected_project <- project
    vers <- nm_workspace_list_versions(project, root = ws_root)
    if (length(vers) > 0L) {
      load_model_version(project, vers[[1L]])
    }
  })

  outputOptions(output, "model_version_header", suspendWhenHidden = FALSE)
  outputOptions(output, "projects_tree", suspendWhenHidden = FALSE)
  outputOptions(output, "versions_tree", suspendWhenHidden = FALSE)
  outputOptions(output, "ctl_problem_display", suspendWhenHidden = FALSE)
  outputOptions(output, "ctl_dataset_select", suspendWhenHidden = FALSE)
  outputOptions(output, "ctl_theta_table", suspendWhenHidden = FALSE)
  outputOptions(output, "ctl_omega_table", suspendWhenHidden = FALSE)
  outputOptions(output, "ctl_sigma_table", suspendWhenHidden = FALSE)
  outputOptions(output, "jobs_tree", suspendWhenHidden = FALSE)
  outputOptions(output, "jobs_refresh_clock", suspendWhenHidden = FALSE)
  outputOptions(output, "compare_gof_grid", suspendWhenHidden = FALSE)
  outputOptions(output, "job_status_banner", suspendWhenHidden = FALSE)
  outputOptions(output, "job_log", suspendWhenHidden = FALSE)
  outputOptions(output, "app_log_banner", suspendWhenHidden = FALSE)
  outputOptions(output, "app_log_history", suspendWhenHidden = FALSE)
}

shinyApp(ui, server)
