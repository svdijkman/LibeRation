#' Generate a modelling report and machine-readable manifest
#'
#' Produces a compact PDF report using base graphics and a JSON manifest with
#' fit provenance and numerical summaries. The report is deterministic and
#' does not require a browser or an external document renderer.
#'
#' @param fit An `nm_fit`.
#' @param file Destination PDF path.
#' @param sections Report sections. Supported values are `summary`,
#'   `parameters`, `gof`, `eta`, `vpc`, and `narrative_stub`.
#' @param vpc Optional `nm_vpc` object for the VPC section.
#' @param title Report title.
#' @param manifest Whether to write a sibling JSON manifest.
#' @return An `nm_report` object containing generated paths and metadata.
#' @export
nm_report <- function(fit, file,
                      sections = c("summary", "parameters", "gof", "eta", "narrative_stub"),
                      vpc = NULL, title = "LibeRation model report",
                      manifest = TRUE) {
  if (!inherits(fit, "nm_fit")) .nm_stop("`fit` must be an nm_fit.")
  file <- path.expand(as.character(file))
  if (length(file) != 1L || is.na(file) || !nzchar(file)) .nm_stop("`file` must be one path.")
  if (!grepl("[.]pdf$", file, ignore.case = TRUE)) file <- paste0(file, ".pdf")
  directory <- dirname(file)
  if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE, showWarnings = FALSE)) {
    .nm_stop("Unable to create report directory: ", directory)
  }
  allowed <- c("summary", "parameters", "gof", "eta", "vpc", "narrative_stub")
  sections <- unique(as.character(sections))
  unknown <- setdiff(sections, allowed)
  if (length(unknown)) .nm_stop("Unknown report section(s): ", paste(unknown, collapse = ", "))
  if ("vpc" %in% sections && !inherits(vpc, "nm_vpc")) sections <- setdiff(sections, "vpc")

  gof <- nm_gof(fit)
  etab <- nm_etab(fit)
  parameters <- c(fit$theta, fit$sigma, fit$omega)
  names(parameters) <- .nm_parameter_names(fit$theta, fit$sigma, fit$omega)
  observed <- gof$EVID == 0L & gof$MDV == 0L & is.finite(gof$DV)
  metadata <- list(
    schema = "liber.report/1", generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC"),
    title = as.character(title), method = fit$method, objective = fit$objective,
    convergence = fit$convergence, records = nrow(gof), subjects = length(unique(gof$ID)),
    observations = sum(observed), parameters = as.list(parameters),
    shrinkage = as.list(etab$shrinkage), sections = sections,
    runtime = list(R = R.version.string, platform = R.version$platform,
                   LibeRation = as.character(utils::packageVersion("LibeRation")))
  )

  report_text_page <- function(heading, lines) {
    graphics::plot.new()
    graphics::par(mar = c(2, 2, 3, 2))
    graphics::title(main = heading, adj = 0, cex.main = 1.35)
    y <- 0.9
    for (line in lines) {
      graphics::text(0.02, y, labels = line, adj = c(0, 1), family = "mono", cex = 0.78)
      y <- y - 0.045
      if (y < 0.05) break
    }
  }

  grDevices::pdf(file, width = 8.27, height = 11.69, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  if ("summary" %in% sections) {
    report_text_page(title, c(
      paste("Generated:", metadata$generated),
      paste("Method:", fit$method),
      paste("Objective:", format(fit$objective, digits = 10)),
      paste("Convergence code:", fit$convergence),
      paste("Subjects:", metadata$subjects),
      paste("Records:", metadata$records),
      paste("Observed DV records:", metadata$observations),
      "",
      paste("ADVAN:", fit$model$ADVAN, "TRANS:", fit$model$TRANS),
      paste("Solver:", fit$model$SOLVER),
      paste("Residual model:", fit$model$LIK_CONFIG$error)
    ))
  }
  if ("parameters" %in% sections) {
    lines <- sprintf("%-16s % .9g", names(parameters), as.numeric(parameters))
    report_text_page("Parameter estimates", c("Parameter        Estimate", strrep("-", 42), lines))
  }
  if ("gof" %in% sections) {
    old_par <- graphics::par(no.readonly = TRUE)
    graphics::par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    graphics::plot(gof$PRED[observed], gof$DV[observed], pch = 16, cex = 0.65,
                   xlab = "PRED", ylab = "DV", main = "DV vs PRED")
    graphics::abline(0, 1, lty = 2, col = "#2b579a")
    graphics::plot(gof$IPRED[observed], gof$DV[observed], pch = 16, cex = 0.65,
                   xlab = "IPRED", ylab = "DV", main = "DV vs IPRED")
    graphics::abline(0, 1, lty = 2, col = "#2b579a")
    graphics::plot(gof$TIME[observed], gof$CWRES[observed], pch = 16, cex = 0.65,
                   col = "#b5484d", ylim = c(-5, 5),
                   xlab = "Time", ylab = "CWRES", main = "CWRES vs time")
    graphics::abline(h = 0, lty = 2, col = "grey50")
    graphics::abline(h = c(-2, 2), lty = 3, col = "#c45b61")
    graphics::plot(gof$PRED[observed], gof$CWRES[observed], pch = 16, cex = 0.65,
                   col = "#b5484d", ylim = c(-5, 5),
                   xlab = "PRED", ylab = "CWRES", main = "CWRES vs PRED")
    graphics::abline(h = 0, lty = 2, col = "grey50")
    graphics::abline(h = c(-2, 2), lty = 3, col = "#c45b61")
    graphics::par(old_par)
  }
  if ("eta" %in% sections) {
    old_par <- graphics::par(no.readonly = TRUE)
    graphics::par(mfrow = c(1, 2), mar = c(6, 4, 2, 1))
    if (ncol(fit$eta)) {
      graphics::boxplot(as.data.frame(fit$eta), las = 2, ylab = "ETA", main = "ETA distributions")
      graphics::barplot(etab$shrinkage, las = 2, ylab = "Shrinkage", main = "ETA shrinkage")
      graphics::abline(h = 0, lty = 2)
    } else {
      graphics::plot.new(); graphics::title("No ETA components")
      graphics::plot.new()
    }
    graphics::par(old_par)
  }
  if ("vpc" %in% sections) {
    report_text_page("Visual predictive check", c(
      paste("Simulations:", vpc$nsim), paste("Interval:", vpc$level),
      "Observed and simulated quantiles are included in the JSON manifest."
    ))
    metadata$vpc <- list(observed = vpc$observed, simulated = vpc$simulated,
                         nsim = vpc$nsim, level = vpc$level)
  }
  if ("narrative_stub" %in% sections) {
    report_text_page("Interpretation notes", c(
      "Model adequacy:", "", "Parameter plausibility:", "",
      "Residual diagnostics:", "", "Predictive performance:", "",
      "Decisions and follow-up analyses:"
    ))
  }
  grDevices::dev.off()
  on.exit(NULL, add = FALSE)
  normalized_pdf <- normalizePath(file, winslash = "/", mustWork = TRUE)
  json_path <- NULL
  if (isTRUE(manifest)) {
    json_path <- sub("[.]pdf$", ".json", file, ignore.case = TRUE)
    jsonlite::write_json(metadata, json_path, auto_unbox = TRUE, pretty = TRUE,
                         digits = NA, null = "null", na = "null")
    if (.Platform$OS.type != "windows") Sys.chmod(json_path, mode = "0600")
    json_path <- normalizePath(json_path, winslash = "/", mustWork = TRUE)
  }
  if (.Platform$OS.type != "windows") Sys.chmod(normalized_pdf, mode = "0600")
  structure(list(pdf = normalized_pdf, json = json_path, metadata = metadata), class = "nm_report")
}
