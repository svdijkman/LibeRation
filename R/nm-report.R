#' Generate a model analysis report (PDF)
#'
#' Designed for pharmacometric workflows and future AI narrative extensions.
#' Writes a PDF plus a JSON manifest alongside for agent consumption.
#'
#' @param fit An \code{nm_fit} object.
#' @param output_path Output \code{.pdf} path (must be inside workspace when called from GUI).
#' @param sections List of logical flags: \code{summary}, \code{parameters},
#'   \code{gof_time}, \code{gof_ipred_time}, \code{gof_scatter}, \code{gof_residuals},
#'   \code{diag_shrinkage}, \code{diag_eta}, \code{narrative_stub}.
#' @param project_meta Optional list with \code{project}, \code{run_id}, \code{workspace}.
#' @return Invisibly, list with \code{pdf} and \code{manifest} paths.
#' @examples
#' \dontrun{
#' sim <- nm_synthetic_theo(n_sub = 2L, seed = 1L)
#' fit <- nm_est(sim$model, sim$data, method = "FO",
#'               control = list(maxit = 5L, compute_inference = FALSE))
#' out <- tempfile(fileext = ".pdf")
#' nm_report_pdf(fit, out, sections = list(summary = TRUE, parameters = TRUE))
#' }
#' @export
nm_report_pdf <- function(fit,
                          output_path,
                          sections = list(
                            summary = TRUE,
                            parameters = TRUE,
                            gof_time = TRUE,
                            gof_ipred_time = TRUE,
                            gof_scatter = TRUE,
                            gof_residuals = TRUE,
                            diag_shrinkage = TRUE,
                            diag_eta = TRUE,
                            narrative_stub = TRUE
                          ),
                          project_meta = list()) {
  has_ggplot <- requireNamespace("ggplot2", quietly = TRUE)
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  pred <- predict(fit, type = "ipred")
  pop <- predict(fit, type = "ppred")
  obs <- pred[pred$MDV == 0L & pred$EVID == 0L, ]

  grDevices::pdf(output_path, width = 8.5, height = 11, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  grid::grid.newpage()
  grid::grid.text("LibeRation Analysis Report", y = 0.85, gp = grid::gpar(fontsize = 18, fontface = "bold"))
  sub <- c(
    if (!is.null(project_meta$project)) paste("Project:", project_meta$project),
    if (!is.null(project_meta$run_id)) paste("Model run:", project_meta$run_id),
    paste("Method:", fit$method),
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  )
  for (i in seq_along(sub)) {
    grid::grid.text(sub[i], y = 0.75 - 0.04 * i, gp = grid::gpar(fontsize = 11))
  }

  if (isTRUE(sections$summary)) {
    grid::grid.newpage()
    grid::grid.text("Fit summary", y = 0.95, gp = grid::gpar(fontsize = 14, fontface = "bold"))
    lines <- c(
      paste("Objective (-2LL):", round(fit$objective, 4)),
      paste("Convergence:", fit$convergence %||% "n/a"),
      paste("Subjects:", length(unique(obs$ID))),
      paste("Observations:", nrow(obs))
    )
    if (!is.null(fit$method) && fit$method %in% c("FOCE", "FOCEI", "FO", "LAPLACE", "SAEM", "IMP")) {
      if (!is.null(fit$grad)) {
        lines <- c(lines, paste("Gradient:", fit$grad))
      }
      if (!is.null(fit$pk_engine)) {
        lines <- c(lines, paste("PK engine:", fit$pk_engine))
      }
    }
    for (i in seq_along(lines)) {
      grid::grid.text(lines[i], x = 0.1, y = 0.85 - 0.05 * i, just = "left",
                      gp = grid::gpar(fontsize = 10))
    }
  }

  if (isTRUE(sections$parameters)) {
    grid::grid.newpage()
    grid::grid.text("Parameter estimates", y = 0.95, gp = grid::gpar(fontsize = 14, fontface = "bold"))
    pt <- nm_par_table(fit)
    if (nrow(pt) > 0L) {
      y <- 0.88
      grid::grid.text("Parameter", x = 0.1, y = y, just = "left", gp = grid::gpar(fontface = "bold"))
      grid::grid.text("Estimate", x = 0.5, y = y, just = "left", gp = grid::gpar(fontface = "bold"))
      y <- y - 0.04
      for (i in seq_len(nrow(pt))) {
        grid::grid.text(pt$parameter[i], x = 0.1, y = y, just = "left", gp = grid::gpar(fontsize = 9))
        grid::grid.text(format(round(pt$estimate[i], 5)), x = 0.5, y = y, just = "left", gp = grid::gpar(fontsize = 9))
        y <- y - 0.035
        if (y < 0.1) break
      }
    }
  }

  if (isTRUE(sections$diag_shrinkage)) {
    sh <- tryCatch(nm_shrinkage(fit), error = function(e) NULL)
    if (!is.null(sh) && nrow(sh) > 0L) {
      grid::grid.newpage()
      grid::grid.text("Random-effect shrinkage", y = 0.95,
                      gp = grid::gpar(fontsize = 14, fontface = "bold"))
      y <- 0.86
      hdr <- c("ETA", "SD(eta)", "SD(omega)", "Shr. %")
      xs <- c(0.08, 0.28, 0.48, 0.68)
      for (j in seq_along(hdr)) {
        grid::grid.text(hdr[j], x = xs[j], y = y, just = "left", gp = grid::gpar(fontface = "bold", fontsize = 9))
      }
      y <- y - 0.04
      for (i in seq_len(nrow(sh))) {
        vals <- c(
          sh$ETA[i],
          format(round(sh$sd_eta[i], 4), nsmall = 4),
          format(round(sh$sd_omega[i], 4), nsmall = 4),
          format(round(100 * sh$shrinkage[i], 1), nsmall = 1)
        )
        for (j in seq_along(vals)) {
          grid::grid.text(vals[j], x = xs[j], y = y, just = "left", gp = grid::gpar(fontsize = 9))
        }
        y <- y - 0.035
      }
    }
  }

  if (isTRUE(sections$diag_eta) && !is.null(fit$eta) && is.matrix(fit$eta) && ncol(fit$eta) > 0L) {
    eta <- fit$eta
    cn <- paste0("ETA", seq_len(ncol(eta)))
    if (has_ggplot) {
      for (j in seq_len(ncol(eta))) {
        p <- ggplot2::ggplot(data.frame(eta = eta[, j]), ggplot2::aes(x = eta)) +
          ggplot2::geom_histogram(bins = 20, fill = "steelblue", color = "white", alpha = 0.85) +
          ggplot2::labs(title = paste(cn[j], "distribution (post-hoc)"), x = cn[j], y = "Count") +
          ggplot2::theme_bw()
        print(p)
      }
    } else {
      grid::grid.newpage()
      grid::grid.text("Post-hoc ETA summaries", y = 0.95,
                      gp = grid::gpar(fontsize = 14, fontface = "bold"))
      y <- 0.85
      for (j in seq_len(ncol(eta))) {
        line <- paste0(
          cn[j], ": mean=", round(mean(eta[, j]), 4),
          ", SD=", round(stats::sd(eta[, j]), 4)
        )
        grid::grid.text(line, x = 0.1, y = y, just = "left", gp = grid::gpar(fontsize = 10))
        y <- y - 0.05
      }
    }
  }

  if (isTRUE(sections$narrative_stub) && nrow(obs) > 0L) {
    grid::grid.newpage()
    grid::grid.text("Scientific interpretation (stub)", y = 0.95,
                    gp = grid::gpar(fontsize = 14, fontface = "bold"))
    stub <- c(
      "This section is reserved for AI-assisted or manual clinical/pharmacological",
      "interpretation. Suggested prompts for an agent:",
      "- Summarize population PK parameters in clinical context.",
      "- Comment on residual patterns and model adequacy.",
      "- Relate findings to dose selection and label-relevant exposure metrics.",
      "- Flag potential covariate effects or model limitations."
    )
    for (i in seq_along(stub)) {
      grid::grid.text(stub[i], x = 0.08, y = 0.82 - 0.045 * i, just = "left",
                      gp = grid::gpar(fontsize = 9))
    }
  }

  if (nrow(obs) > 0L) {
    ind_lines <- pred[pred$MDV == 0L & pred$EVID == 0L, ]
    ind_lines <- ind_lines[order(ind_lines$ID, ind_lines$TIME), ]
    pop_lines <- pop[pop$MDV == 0L & pop$EVID == 0L, ]
    pop_lines <- pop_lines[order(pop_lines$ID, pop_lines$TIME), ]

    if (isTRUE(sections$gof_time)) {
      if (has_ggplot) {
        p <- ggplot2::ggplot() +
          ggplot2::geom_line(
            data = ind_lines,
            ggplot2::aes(x = TIME, y = IPRED, group = factor(ID)),
            color = "#D95F02", linewidth = 0.5
          ) +
          ggplot2::geom_point(data = obs, ggplot2::aes(x = TIME, y = DV), size = 1.2) +
          ggplot2::labs(title = "DV vs time with IPRED", x = "Time", y = "DV") +
          ggplot2::theme_bw()
        print(p)
      } else {
        .nm_report_plot_time(obs, ind_lines)
      }
    }

    if (isTRUE(sections$gof_ipred_time)) {
      if (has_ggplot) {
        p <- ggplot2::ggplot(ind_lines, ggplot2::aes(x = TIME, y = IPRED, group = factor(ID))) +
          ggplot2::geom_line(color = "#D95F02", linewidth = 0.5) +
          ggplot2::geom_line(
            data = pop_lines,
            ggplot2::aes(x = TIME, y = PRED, group = factor(ID)),
            color = "#1B9E77", linewidth = 0.4, linetype = "dashed"
          ) +
          ggplot2::labs(title = "IPRED vs time (solid) with PRED (dashed)", x = "Time", y = "Prediction") +
          ggplot2::theme_bw()
        print(p)
      }
    }

    if (isTRUE(sections$gof_scatter)) {
      if (has_ggplot) {
        lims <- range(c(obs$IPRED, obs$DV), finite = TRUE)
        p <- ggplot2::ggplot(obs, ggplot2::aes(x = IPRED, y = DV)) +
          ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
          ggplot2::geom_point(alpha = 0.7) +
          ggplot2::coord_equal(xlim = lims, ylim = lims) +
          ggplot2::labs(title = "DV vs IPRED") +
          ggplot2::theme_bw()
        print(p)
        lims2 <- range(c(obs$PRED, obs$DV), finite = TRUE)
        p2 <- ggplot2::ggplot(obs, ggplot2::aes(x = PRED, y = DV)) +
          ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
          ggplot2::geom_point(alpha = 0.7, color = "#7570B3") +
          ggplot2::coord_equal(xlim = lims2, ylim = lims2) +
          ggplot2::labs(title = "DV vs PRED") +
          ggplot2::theme_bw()
        print(p2)
      } else {
        .nm_report_plot_scatter(obs, "IPRED", "DV vs IPRED")
        .nm_report_plot_scatter(obs, "PRED", "DV vs PRED")
      }
    }

    if (isTRUE(sections$gof_residuals)) {
      .nm_report_residual_pages(obs, has_ggplot)
    }
  }

  manifest_path <- sub("\\.pdf$", "_manifest.json", output_path, ignore.case = TRUE)
  manifest <- list(
    report_version = 2L,
    generated = as.character(Sys.time()),
    project = project_meta$project %||% "",
    run_id = project_meta$run_id %||% "",
    method = fit$method,
    objective = fit$objective,
    theta = as.list(stats::setNames(fit$theta, paste0("THETA", fit$model$THETAS$THETA))),
    omega = if (length(fit$omega) > 0L) as.list(fit$omega) else list(),
    sigma = if (length(fit$sigma) > 0L) as.list(fit$sigma) else list(),
    n_subjects = length(unique(obs$ID)),
    n_obs = nrow(obs),
    sections = sections,
    ai_narrative_placeholder = TRUE
  )
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE, pretty = TRUE)
  }

  invisible(list(pdf = output_path, manifest = manifest_path))
}

#' @keywords internal
.nm_report_residual_pages <- function(obs, has_ggplot = TRUE) {
  plots <- list(
    list(var = "WRES", x = "TIME", title = "WRES vs time", xlab = "Time"),
    list(var = "IWRES", x = "TIME", title = "IWRES vs time", xlab = "Time"),
    list(var = "WRES", x = "IPRED", title = "WRES vs IPRED", xlab = "IPRED"),
    list(var = "IWRES", x = "IPRED", title = "IWRES vs IPRED", xlab = "IPRED"),
    list(var = "RES", x = "PRED", title = "RES vs PRED", xlab = "PRED")
  )
  for (pl in plots) {
    v <- pl$var
    if (!v %in% names(obs)) {
      next
    }
    ok <- !is.na(obs[[v]]) & is.finite(obs[[v]])
    if (!any(ok)) {
      next
    }
    sub <- obs[ok, , drop = FALSE]
    if (has_ggplot) {
      p <- ggplot2::ggplot(sub, ggplot2::aes(x = .data[[pl$x]], y = .data[[v]])) +
        ggplot2::geom_hline(yintercept = 0, linetype = 2) +
        ggplot2::geom_point(alpha = 0.7) +
        ggplot2::labs(title = pl$title, x = pl$xlab, y = v) +
        ggplot2::theme_bw()
      print(p)
    } else if (identical(v, "WRES") && identical(pl$x, "TIME")) {
      .nm_report_plot_wres(sub)
    }
  }
  ok_res <- !is.na(obs$RES) & is.finite(obs$RES)
  if (any(ok_res) && has_ggplot) {
    p <- ggplot2::ggplot(obs[ok_res, ], ggplot2::aes(x = RES)) +
      ggplot2::geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.85) +
      ggplot2::labs(title = "RES histogram", x = "RES", y = "Count") +
      ggplot2::theme_bw()
    print(p)
  }
}

#' @keywords internal
.nm_report_plot_time <- function(obs, ind_lines) {
  graphics::plot(
    obs$TIME, obs$DV, pch = 16, col = "#333333", cex = 0.7,
    xlab = "Time", ylab = "DV", main = "DV vs time with IPRED"
  )
  ids <- unique(ind_lines$ID)
  for (id in ids) {
    sub <- ind_lines[ind_lines$ID == id, ]
    graphics::lines(sub$TIME, sub$IPRED, col = "#D95F02", lwd = 1)
  }
}

#' @keywords internal
.nm_report_plot_scatter <- function(obs, xvar, title) {
  x <- obs[[xvar]]
  y <- obs$DV
  lims <- range(c(x, y), finite = TRUE)
  graphics::plot(x, y, pch = 16, col = "steelblue", cex = 0.75,
       xlab = xvar, ylab = "DV", main = title, xlim = lims, ylim = lims)
  graphics::abline(0, 1, lty = 2, col = "gray50")
}

#' @keywords internal
.nm_report_plot_wres <- function(obs) {
  graphics::plot(obs$TIME, obs$WRES, pch = 16, col = "steelblue", cex = 0.75,
       xlab = "Time", ylab = "WRES", main = "WRES vs time")
  graphics::abline(h = 0, lty = 2, col = "gray50")
}
