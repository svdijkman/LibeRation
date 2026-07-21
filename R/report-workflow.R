.nm_report_block_types <- c(
  "title", "introduction", "methods", "run", "comparison", "discussion",
  "conclusion", "appendix", "text", "page_break"
)

.nm_report_block_source <- function(type, source) {
  if (type %in% c("run", "comparison")) return("run")
  if (type %in% c("title", "page_break")) return("user")
  match.arg(source, c("user", "ai"))
}

#' Define one block in a visual report workflow
#'
#' @param type Report-block type.
#' @param title Display heading.
#' @param source Narrative source: user-provided or browser-local AI. Run and
#'   comparison blocks always use saved run evidence.
#' @param text Current editable text for narrative blocks.
#' @param file Optional TXT, Markdown, or DOCX source for user text.
#' @param run_ids Immutable model-run ids used by run/comparison blocks.
#' @param elements Run evidence to include, such as summary, parameters, code,
#'   GOF, continuous/categorical/count/event VPCs, NPDE, NPC, covariance, and
#'   run information.
#' @param options Named block options, including AI instructions, templates,
#'   and locally extracted source documents.
#' @param id Optional stable block id.
#' @return A serializable `nm_report_block`.
#' @export
nm_report_block <- function(type, title = NULL, source = c("user", "ai"),
                            text = "", file = NULL, run_ids = NULL,
                            elements = NULL, options = list(), id = NULL) {
  type <- match.arg(tolower(as.character(type)), .nm_report_block_types)
  source <- .nm_report_block_source(type, source)
  if (is.null(id)) {
    id <- paste0("block-", format(Sys.time(), "%Y%m%d%H%M%OS3", tz = "UTC"), "-",
                 sprintf("%06x", sample.int(0xffffff, 1L)))
  }
  id <- .nm_workspace_component(gsub("[^A-Za-z0-9_.-]", "-", id), "report block id")
  default_title <- switch(
    type, title = "Title", introduction = "Introduction", methods = "Methods",
    run = "Model run", comparison = "Model comparison", discussion = "Discussion",
    conclusion = "Conclusion", appendix = "Appendix", text = "Text",
    page_break = "Page break"
  )
  title <- trimws(as.character(title %||% default_title))
  text <- paste(as.character(text %||% ""), collapse = "\n")
  file <- if (is.null(file)) NULL else path.expand(as.character(file))
  run_ids <- unique(as.character(run_ids %||% character()))
  elements <- unique(tolower(as.character(elements %||%
    if (type %in% c("run", "comparison")) c("summary", "parameters") else character())))
  allowed_elements <- c(
    "summary", "parameters", "code", "gof", "vpc", "vpc_categorical",
    "vpc_count", "vpc_tte", "vpc_competing", "vpc_recurrent",
    "npde", "npc", "covariance", "run_info"
  )
  if (length(setdiff(elements, allowed_elements))) {
    .nm_stop("Unknown report evidence element(s): ",
             paste(setdiff(elements, allowed_elements), collapse = ", "), ".")
  }
  if (!is.list(options) || (length(options) && is.null(names(options)))) {
    .nm_stop("Report-block `options` must be a named list.")
  }
  structure(list(
    schema = "liber.report-block/1", version = 1L, id = id, type = type,
    title = title, source = source, text = text, file = file,
    run_ids = run_ids, elements = elements, options = options
  ), class = "nm_report_block")
}

.nm_report_normalize_block <- function(block) {
  if (inherits(block, "nm_report_block")) return(block)
  if (!is.list(block) || is.null(block$type)) {
    .nm_stop("Every report workflow item must be an nm_report_block or contain `type`.")
  }
  arguments <- block[intersect(names(block), names(formals(nm_report_block)))]
  do.call(nm_report_block, arguments)
}

#' Define a visual report workflow
#'
#' @param blocks Ordered report blocks.
#' @param title Report title.
#' @param formats Output formats; DOCX and PDF are supported.
#' @param style Named rendering options.
#' @param id Optional stable design id.
#' @return A serializable `nm_report_design`.
#' @export
nm_report_design <- function(blocks, title = "LibeRation modelling report",
                             formats = c("docx", "pdf"), style = list(), id = NULL) {
  if (!is.list(blocks) || !length(blocks)) .nm_stop("A report design requires at least one block.")
  blocks <- lapply(blocks, .nm_report_normalize_block)
  block_ids <- vapply(blocks, `[[`, character(1), "id")
  if (anyDuplicated(block_ids)) .nm_stop("Report block ids must be unique.")
  formats <- unique(tolower(as.character(formats)))
  if (!length(formats) || length(setdiff(formats, c("docx", "pdf")))) {
    .nm_stop("Report formats must be DOCX and/or PDF.")
  }
  if (!is.list(style) || (length(style) && is.null(names(style)))) {
    .nm_stop("Report `style` must be a named list.")
  }
  if (is.null(id)) {
    id <- paste0("report-", format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"), "-",
                 sprintf("%08x", sample.int(.Machine$integer.max, 1L)))
  }
  id <- .nm_workspace_component(id, "report design id")
  structure(list(
    schema = "liber.report-design/1", version = 1L, id = id,
    title = trimws(as.character(title)), formats = formats, style = style,
    blocks = blocks, updated = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")
  ), class = "nm_report_design")
}

.nm_report_design_directory <- function(workspace, project) {
  directory <- file.path(.nm_project_path(workspace, project), "reports")
  if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE, showWarnings = FALSE)) {
    .nm_stop("Unable to create the project report-design directory.")
  }
  directory
}

#' Save a report workflow with a project
#'
#' @param workspace An [nm_workspace()] or path.
#' @param project Project id.
#' @param design An [nm_report_design()].
#' @return The design id, invisibly.
#' @export
nm_report_design_save <- function(workspace, project, design) {
  workspace <- if (inherits(workspace, "nm_workspace")) workspace else nm_workspace(workspace)
  if (!inherits(design, "nm_report_design")) .nm_stop("`design` must be an nm_report_design.")
  design$updated <- format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")
  path <- file.path(.nm_report_design_directory(workspace, project), paste0(design$id, ".rds"))
  .nm_workspace_atomic_save(design, path)
  invisible(design$id)
}

#' Load or list project report workflows
#'
#' @param workspace An [nm_workspace()] or path.
#' @param project Project id.
#' @param id Optional report-design id. When omitted, metadata are listed.
#' @return An `nm_report_design` or metadata data frame.
#' @export
nm_report_design_load <- function(workspace, project, id = NULL) {
  workspace <- if (inherits(workspace, "nm_workspace")) workspace else nm_workspace(workspace)
  directory <- .nm_report_design_directory(workspace, project)
  if (!is.null(id)) {
    id <- .nm_workspace_component(id, "report design id")
    value <- .nm_workspace_read(file.path(directory, paste0(id, ".rds")))
    if (!inherits(value, "nm_report_design")) .nm_stop("Stored report design is invalid.")
    return(value)
  }
  files <- list.files(directory, pattern = "[.]rds$", full.names = TRUE)
  rows <- lapply(files, function(path) tryCatch({
    value <- .nm_workspace_read(path)
    data.frame(
      id = value$id, title = value$title, updated = value$updated,
      blocks = length(value$blocks), formats = paste(value$formats, collapse = ", "),
      stringsAsFactors = FALSE
    )
  }, error = function(error) NULL))
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(data.frame(
    id = character(), title = character(), updated = character(),
    blocks = integer(), formats = character(), stringsAsFactors = FALSE
  ))
  output <- do.call(rbind, rows)
  output[order(output$updated, decreasing = TRUE), , drop = FALSE]
}

.nm_report_decode_entities <- function(text) {
  replacements <- c("&amp;" = "&", "&lt;" = "<", "&gt;" = ">",
                    "&quot;" = "\"", "&apos;" = "'")
  for (name in names(replacements)) text <- gsub(name, replacements[[name]], text, fixed = TRUE)
  text
}

.nm_report_document_text <- function(path) {
  if (is.null(path) || !length(path) || !nzchar(path)) return("")
  path <- normalizePath(path.expand(path), winslash = "/", mustWork = TRUE)
  extension <- tolower(tools::file_ext(path))
  if (extension %in% c("txt", "md", "markdown", "rmd")) {
    return(paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
  }
  if (extension != "docx") .nm_stop("User report text files must be TXT, Markdown, or DOCX.")
  directory <- tempfile("liber-docx-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE, force = TRUE), add = TRUE)
  extracted <- tryCatch(
    utils::unzip(path, files = "word/document.xml", exdir = directory),
    error = function(error) character()
  )
  xml <- file.path(directory, "word", "document.xml")
  if (!length(extracted) || !file.exists(xml)) .nm_stop("Unable to read DOCX document text.")
  text <- paste(readLines(xml, warn = FALSE, encoding = "UTF-8"), collapse = "")
  text <- gsub("<w:tab[^>]*/>", "\t", text, perl = TRUE)
  text <- gsub("</w:p>", "\n", text, fixed = TRUE)
  text <- gsub("<[^>]+>", "", text, perl = TRUE)
  trimws(.nm_report_decode_entities(text))
}

.nm_report_block_text <- function(block) {
  file_text <- if (!is.null(block$file)) .nm_report_document_text(block$file) else ""
  text <- if (nzchar(file_text)) file_text else block$text
  if (!nzchar(trimws(text)) && block$source == "ai") {
    text <- "[Local AI text has not yet been generated for this block.]"
  }
  text
}

.nm_report_run_evidence <- function(workspace, project, run_id) {
  opened <- nm_project_load(workspace, project, run_id)
  if (!inherits(opened$result, "nm_fit") &&
      !(is.data.frame(opened$result) && "IPRED" %in% names(opened$result))) {
    .nm_stop("Report run `", run_id, "` is not an estimation or simulation run.")
  }
  list(id = run_id, model = opened$model, data = opened$data,
       result = opened$result,
       diagnostics = nm_project_load_diagnostics(workspace, project, run_id))
}

.nm_report_parameter_table <- function(fit) {
  values <- c(fit$theta, fit$sigma, fit$omega)
  names(values) <- .nm_parameter_names(fit$theta, fit$sigma, fit$omega)
  data.frame(Parameter = names(values), Estimate = as.numeric(values),
             stringsAsFactors = FALSE)
}

.nm_report_markdown_table <- function(table) {
  table <- as.data.frame(table, stringsAsFactors = FALSE)
  if (!nrow(table)) return("_No rows available._")
  escape <- function(value) gsub("[|]", "\\|", as.character(value))
  header <- paste0("| ", paste(escape(names(table)), collapse = " | "), " |")
  rule <- paste0("| ", paste(rep("---", ncol(table)), collapse = " | "), " |")
  rows <- apply(table, 1L, function(row) paste0("| ", paste(escape(row), collapse = " | "), " |"))
  paste(c(header, rule, rows), collapse = "\n")
}

.nm_report_draw_gof <- function(fit) {
  gof <- nm_gof(fit)
  observed <- gof$EVID == 0L & gof$MDV == 0L & is.finite(gof$DV)
  graphics::par(mfrow = c(2, 2), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))
  graphics::plot(gof$PRED[observed], gof$DV[observed], pch = 16, cex = .55,
                 xlab = "PRED", ylab = "DV", main = "DV vs PRED")
  graphics::abline(0, 1, lty = 2)
  graphics::plot(gof$IPRED[observed], gof$DV[observed], pch = 16, cex = .55,
                 xlab = "IPRED", ylab = "DV", main = "DV vs IPRED")
  graphics::abline(0, 1, lty = 2)
  graphics::plot(gof$TIME[observed], gof$CWRES[observed], pch = 16, cex = .55,
                 col = "#b5484d", ylim = c(-5, 5),
                 xlab = "TIME", ylab = "CWRES", main = "CWRES vs TIME")
  graphics::abline(h = c(-2, 0, 2), lty = c(3, 2, 3))
  graphics::plot(gof$PRED[observed], gof$CWRES[observed], pch = 16, cex = .55,
                 col = "#b5484d", ylim = c(-5, 5),
                 xlab = "PRED", ylab = "CWRES", main = "CWRES vs PRED")
  graphics::abline(h = c(-2, 0, 2), lty = c(3, 2, 3))
  invisible(gof)
}

.nm_report_plot_gof <- function(fit, file) {
  grDevices::png(file, width = 1600, height = 1200, res = 180)
  on.exit(grDevices::dev.off(), add = TRUE)
  .nm_report_draw_gof(fit)
  grDevices::dev.off()
  on.exit(NULL, add = FALSE)
  file
}

.nm_report_draw_diagnostic <- function(result, kind) {
  kind <- as.character(kind)
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  if (kind == "vpc" && inherits(result, "nm_vpc")) {
    observed <- result$observed; simulated <- result$simulated
    qnames <- grep("_median$", names(simulated), value = TRUE)
    values <- unlist(c(observed[grep("^Q", names(observed))],
                       simulated[grep("_(lo|median|hi)$", names(simulated))]))
    graphics::plot(range(c(observed$TIME, simulated$TIME), finite = TRUE),
                   range(values, finite = TRUE), type = "n",
                   xlab = "Time", ylab = "DV", main = "Visual predictive check")
    colors <- c("#77aee0", "#e7898d", "#77aee0")
    for (index in seq_along(qnames)) {
      base <- sub("_median$", "", qnames[[index]])
      lo <- simulated[[paste0(base, "_lo")]]; hi <- simulated[[paste0(base, "_hi")]]
      graphics::polygon(c(simulated$TIME, rev(simulated$TIME)), c(lo, rev(hi)),
                        col = grDevices::adjustcolor(colors[[min(index, 3L)]], .28), border = NA)
      graphics::lines(simulated$TIME, simulated[[qnames[[index]]]], lty = 2)
      if (base %in% names(observed)) graphics::lines(observed$TIME, observed[[base]], lwd = 1.5)
    }
    if (!is.null(result$points)) graphics::points(result$points$TIME, result$points$DV, pch = 16, cex = .35)
  } else if (kind == "npde" && inherits(result, "nm_npde")) {
    graphics::par(mfrow = c(1, 2))
    graphics::hist(result$table$NPDE, main = "NPDE distribution", xlab = "NPDE", col = "#77aee0")
    graphics::plot(result$table$TIME, result$table$NPDE, pch = 16, cex = .5,
                   xlab = "Time", ylab = "NPDE", main = "NPDE vs time")
    graphics::abline(h = 0, lty = 2)
  } else if (kind == "npc" && inherits(result, "nm_npc")) {
    graphics::par(mfrow = c(1, 2))
    graphics::hist(result$table$PERCENTILE, breaks = seq(0, 1, .1),
                   main = "Predictive percentiles", xlab = "Percentile", col = "#77aee0")
    graphics::plot(result$table$TIME, result$table$PERCENTILE, pch = 16, cex = .5,
                   xlab = "Time", ylab = "Percentile", main = "NPC vs time")
  } else if (kind == "vpc_categorical" && inherits(result, "nm_vpc_categorical")) {
    simulated <- result$simulated; observed <- result$observed
    graphics::plot(range(simulated$TIME), c(0, 1), type = "n", xlab = "Time",
                    ylab = "Proportion", main = "Categorical VPC")
    categories <- unique(as.character(simulated$CATEGORY %||% "outcome"))
    colors <- grDevices::hcl.colors(max(3L, length(categories)), "Dark 3")
    for (index in seq_along(categories)) {
      category <- categories[[index]]
      sim <- simulated[as.character(simulated$CATEGORY %||% "outcome") == category, , drop = FALSE]
      obs <- observed[as.character(observed$CATEGORY %||% "outcome") == category, , drop = FALSE]
      graphics::polygon(c(sim$TIME, rev(sim$TIME)), c(sim$lower, rev(sim$upper)),
                        col = grDevices::adjustcolor(colors[[index]], .2), border = NA)
      graphics::lines(sim$TIME, sim$median, lty = 2, col = colors[[index]])
      graphics::lines(obs$TIME, obs$PROPORTION, type = "b", pch = 16,
                      col = colors[[index]])
    }
    graphics::legend("topright", legend = categories,
                     col = colors[seq_along(categories)], lty = 1, pch = 16,
                     bty = "n", cex = .8)
  } else if (kind == "vpc_count" && inherits(result, "nm_vpc_count")) {
    simulated <- result$simulated; observed <- result$observed
    limits <- range(c(simulated$MEAN_lower, simulated$MEAN_upper, observed$MEAN),
                    finite = TRUE)
    graphics::plot(range(simulated$TIME), limits, type = "n", xlab = "Time",
                   ylab = "Mean count", main = "Count VPC")
    graphics::polygon(c(simulated$TIME, rev(simulated$TIME)),
                      c(simulated$MEAN_lower, rev(simulated$MEAN_upper)),
                      col = grDevices::adjustcolor("#77aee0", .35), border = NA)
    graphics::lines(simulated$TIME, simulated$MEAN_median, lty = 2)
    graphics::lines(observed$TIME, observed$MEAN, type = "b", pch = 16)
  } else if (kind == "vpc_tte" && inherits(result, "nm_vpc_tte")) {
    simulated <- result$simulated; observed <- result$observed
    graphics::plot(range(simulated$TIME), c(0, 1), type = "n", xlab = "Time",
                   ylab = "Event-free survival", main = "Time-to-event VPC")
    graphics::polygon(c(simulated$TIME, rev(simulated$TIME)),
                      c(simulated$lower, rev(simulated$upper)),
                      col = grDevices::adjustcolor("#77aee0", .35), border = NA)
    graphics::lines(simulated$TIME, simulated$median, lty = 2)
    graphics::lines(observed$TIME, observed$SURVIVAL, lwd = 1.5)
  } else if (kind == "vpc_competing" && inherits(result, "nm_vpc_competing")) {
    simulated <- result$simulated; observed <- result$observed
    causes <- unique(as.character(simulated$CAUSE))
    colors <- grDevices::hcl.colors(max(3L, length(causes)), "Dark 3")
    graphics::plot(range(simulated$TIME), c(0, 1), type = "n", xlab = "Time",
                   ylab = "Cumulative incidence", main = "Competing-risk VPC")
    for (index in seq_along(causes)) {
      cause <- causes[[index]]
      sim <- simulated[as.character(simulated$CAUSE) == cause, , drop = FALSE]
      obs <- observed[as.character(observed$CAUSE) == cause, , drop = FALSE]
      graphics::polygon(c(sim$TIME, rev(sim$TIME)), c(sim$lower, rev(sim$upper)),
                        col = grDevices::adjustcolor(colors[[index]], .2), border = NA)
      graphics::lines(sim$TIME, sim$median, lty = 2, col = colors[[index]])
      graphics::lines(obs$TIME, obs$CIF, type = "b", pch = 16,
                      col = colors[[index]])
    }
    graphics::legend("topleft", legend = causes,
                     col = colors[seq_along(causes)], lty = 1, pch = 16,
                     bty = "n", cex = .8)
  } else if (kind == "vpc_recurrent" && inherits(result, "nm_vpc_recurrent")) {
    simulated <- result$simulated; observed <- result$observed
    limits <- range(c(simulated$lower, simulated$upper, observed$MEAN_CUMULATIVE),
                    finite = TRUE)
    graphics::plot(range(simulated$TIME), limits, type = "n", xlab = "Time",
                   ylab = "Mean cumulative events", main = "Recurrent-event VPC")
    graphics::polygon(c(simulated$TIME, rev(simulated$TIME)),
                      c(simulated$lower, rev(simulated$upper)),
                      col = grDevices::adjustcolor("#77aee0", .35), border = NA)
    graphics::lines(simulated$TIME, simulated$median, lty = 2)
    graphics::lines(observed$TIME, observed$MEAN_CUMULATIVE,
                    type = "b", pch = 16)
  } else {
    graphics::plot.new(); graphics::title(main = toupper(gsub("_", " ", kind)))
    graphics::text(.5, .5, "Saved diagnostic is available in the report manifest.")
  }
  invisible(result)
}

.nm_report_plot_diagnostic <- function(result, kind, file) {
  grDevices::png(file, width = 1600, height = 1000, res = 180)
  on.exit(grDevices::dev.off(), add = TRUE)
  .nm_report_draw_diagnostic(result, kind)
  grDevices::dev.off(); on.exit(NULL, add = FALSE)
  file
}

.nm_report_materialize <- function(design, workspace, project, assets) {
  content <- list()
  add <- function(value) content[[length(content) + 1L]] <<- value
  for (block in design$blocks) {
    if (block$type == "page_break") {
      add(list(type = "page_break")); next
    }
    if (!block$type %in% c("run", "comparison")) {
      add(list(type = "narrative", heading = block$title,
               text = .nm_report_block_text(block), level = if (block$type == "title") 1L else 2L))
      next
    }
    if (!length(block$run_ids)) {
      add(list(type = "narrative", heading = block$title,
               text = "[No model run selected.]", level = 2L)); next
    }
    if (block$type == "comparison" && length(block$run_ids) > 1L &&
        "parameters" %in% block$elements) {
      comparison <- lapply(block$run_ids, function(run_id) {
        .nm_report_run_evidence(workspace, project, run_id)$result
      })
      if (all(vapply(comparison, inherits, logical(1), "nm_fit"))) {
        vectors <- lapply(comparison, function(fit) {
          values <- c(fit$theta, fit$sigma, fit$omega)
          names(values) <- .nm_parameter_names(fit$theta, fit$sigma, fit$omega)
          values
        })
        parameter_names <- unique(unlist(lapply(vectors, names), use.names = FALSE))
        table <- data.frame(Parameter = parameter_names, stringsAsFactors = FALSE)
        for (index in seq_along(vectors)) {
          table[[block$run_ids[[index]]]] <- unname(vectors[[index]][parameter_names])
        }
        add(list(type = "table", heading = paste(block$title, "parameter comparison"),
                 table = table))
      }
    }
    for (run_id in block$run_ids) {
      evidence <- .nm_report_run_evidence(workspace, project, run_id)
      fit <- evidence$result
      heading <- paste0(block$title, if (length(block$run_ids) > 1L) paste0(" - ", run_id) else "")
      add(list(type = "heading", heading = heading, level = 2L))
      if (inherits(fit, "nm_fit") && "summary" %in% block$elements) {
        add(list(type = "narrative", heading = NULL, text = paste(
          "Method:", as.character(fit$method), "\nObjective:", format(fit$objective, digits = 10),
          "\nConvergence code:", fit$convergence
        ), level = 3L))
      }
      if (inherits(fit, "nm_fit") && "parameters" %in% block$elements) {
        add(list(type = "table", heading = "Parameter estimates",
                 table = .nm_report_parameter_table(fit)))
      }
      if ("code" %in% block$elements) {
        add(list(type = "code", heading = "Model code",
                 text = paste("$PK / $PRED", evidence$model$PRED,
                              "$DES", evidence$model$DES,
                              "$ERROR", evidence$model$ERROR, sep = "\n")))
      }
      if (inherits(fit, "nm_fit") && "gof" %in% block$elements) {
        figure <- file.path(assets, paste0("gof-", gsub("[^A-Za-z0-9_.-]", "-", run_id), ".png"))
        .nm_report_plot_gof(fit, figure)
        add(list(type = "figure", heading = "Goodness of fit", file = figure, fit = fit))
      }
      if (inherits(fit, "nm_fit") && "covariance" %in% block$elements &&
          !is.null(fit$covariance$covariance)) {
        covariance <- as.data.frame(fit$covariance$covariance, check.names = FALSE)
        covariance <- data.frame(Parameter = rownames(fit$covariance$covariance), covariance,
                                 check.names = FALSE)
        add(list(type = "table", heading = "Covariance matrix", table = covariance))
      }
      if (inherits(fit, "nm_fit") && "run_info" %in% block$elements) {
        timing <- fit$timing %||% list()
        add(list(type = "table", heading = "Run information", table = data.frame(
          Item = c("Method", "Objective", "Convergence", "Iterations",
                   "Model fit seconds", "Covariance seconds", "Total seconds"),
          Value = c(fit$method, fit$objective, fit$convergence,
                    fit$iterations %||% NA, timing$model_fit_seconds %||% NA,
                    timing$covariance_seconds %||% NA, timing$total_seconds %||% NA),
          stringsAsFactors = FALSE
        )))
      }
      diagnostic_kinds <- intersect(
        block$elements,
        c("vpc", "vpc_categorical", "vpc_count", "vpc_tte",
          "vpc_competing", "vpc_recurrent", "npde", "npc")
      )
      for (kind in diagnostic_kinds) {
        diagnostic <- evidence$diagnostics[[kind]]
        if (is.null(diagnostic)) next
        figure <- file.path(
          assets, paste0(kind, "-", gsub("[^A-Za-z0-9_.-]", "-", run_id), ".png")
        )
        .nm_report_plot_diagnostic(diagnostic, kind, figure)
        add(list(
          type = "figure", heading = toupper(gsub("_", " ", kind)),
          file = figure, diagnostic = diagnostic, diagnostic_kind = kind
        ))
      }
    }
  }
  content
}

.nm_report_write_markdown <- function(design, content, file) {
  lines <- c(paste0("# ", design$title), "")
  for (item in content) {
    if (item$type == "page_break") {
      lines <- c(lines, "\\newpage", ""); next
    }
    if (!is.null(item$heading) && nzchar(item$heading)) {
      level <- item$level %||% 2L
      lines <- c(lines, paste0(strrep("#", level), " ", item$heading), "")
    }
    if (item$type %in% c("narrative", "heading")) {
      if (!is.null(item$text)) lines <- c(lines, item$text, "")
    } else if (item$type == "table") {
      lines <- c(lines, .nm_report_markdown_table(item$table), "")
    } else if (item$type == "code") {
      lines <- c(lines, "```r", item$text, "```", "")
    } else if (item$type == "figure") {
      lines <- c(lines, paste0("![](", normalizePath(item$file, winslash = "/"), ")"), "")
    }
  }
  writeLines(lines, file, useBytes = TRUE)
  invisible(file)
}

.nm_report_pandoc <- function() {
  executable <- Sys.which("pandoc")
  if (nzchar(executable)) return(unname(executable))
  if (.Platform$OS.type == "windows") {
    candidates <- c(
      file.path(Sys.getenv("ProgramFiles", "C:/Program Files"), "Pandoc", "pandoc.exe"),
      file.path(Sys.getenv("LOCALAPPDATA", ""), "Pandoc", "pandoc.exe")
    )
    available <- candidates[file.exists(candidates)]
    if (length(available)) return(available[[1L]])
  }
  if (requireNamespace("rmarkdown", quietly = TRUE)) {
    location <- tryCatch(rmarkdown::find_pandoc(), error = function(error) NULL)
    if (!is.null(location) && nzchar(location$dir)) {
      candidate <- file.path(location$dir, if (.Platform$OS.type == "windows") "pandoc.exe" else "pandoc")
      if (file.exists(candidate)) return(candidate)
    }
  }
  ""
}

.nm_report_render_docx <- function(design, content, file, assets) {
  pandoc <- .nm_report_pandoc()
  if (!nzchar(pandoc)) .nm_stop("DOCX report generation requires Pandoc.")
  markdown <- file.path(assets, "report.md")
  .nm_report_write_markdown(design, content, markdown)
  arguments <- c(markdown, "--from=gfm", "--standalone", "--output", file)
  reference <- design$style$reference_docx %||% NULL
  if (!is.null(reference) && file.exists(reference)) {
    arguments <- c(arguments, "--reference-doc", normalizePath(reference, winslash = "/"))
  }
  output <- system2(pandoc, arguments, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status") %||% 0L
  if (status != 0L || !file.exists(file)) {
    .nm_stop("Pandoc could not generate DOCX: ", paste(output, collapse = "\n"))
  }
  file
}

.nm_report_pdf_text <- function(heading, text) {
  lines <- unlist(lapply(strsplit(text %||% "", "\n", fixed = TRUE)[[1L]], function(line) {
    if (!nzchar(line)) "" else strwrap(line, width = 94)
  }), use.names = FALSE)
  if (!length(lines)) lines <- ""
  chunks <- split(lines, ceiling(seq_along(lines) / 44))
  for (chunk in chunks) {
    graphics::par(mfrow = c(1, 1), mar = c(2, 2, 3, 2), oma = c(0, 0, 0, 0))
    graphics::plot.new()
    if (!is.null(heading) && nzchar(heading)) graphics::title(main = heading, adj = 0, cex.main = 1.25)
    y <- .92
    for (line in chunk) {
      graphics::text(.02, y, line, adj = c(0, 1), cex = .78)
      y <- y - .02
    }
    heading <- NULL
  }
}

.nm_report_render_pdf <- function(design, content, file) {
  grDevices::pdf(file, width = 8.27, height = 11.69, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  .nm_report_pdf_text(design$title, "")
  for (item in content) {
    if (item$type == "page_break") next
    if (item$type == "heading") {
      .nm_report_pdf_text(item$heading, "")
      next
    }
    if (item$type %in% c("narrative", "code")) {
      .nm_report_pdf_text(item$heading, item$text)
    } else if (item$type == "table") {
      table <- as.data.frame(item$table, stringsAsFactors = FALSE)
      lines <- utils::capture.output(print(table, row.names = FALSE, right = FALSE))
      .nm_report_pdf_text(item$heading, paste(lines, collapse = "\n"))
    } else if (item$type == "figure") {
      if (!is.null(item$fit) && inherits(item$fit, "nm_fit")) {
        .nm_report_draw_gof(item$fit)
        graphics::mtext(item$heading, side = 3, outer = TRUE, line = -1.2,
                        adj = 0, cex = 1.05, font = 2)
      } else if (!is.null(item$diagnostic)) {
        graphics::par(mfrow = c(1, 1), oma = c(0, 0, 2, 0))
        .nm_report_draw_diagnostic(item$diagnostic, item$diagnostic_kind)
        graphics::mtext(item$heading, side = 3, outer = TRUE, line = -1.2,
                        adj = 0, cex = 1.05, font = 2)
      } else {
        .nm_report_pdf_text(item$heading, paste("Figure:", item$file))
      }
    }
  }
  grDevices::dev.off(); on.exit(NULL, add = FALSE)
  file
}

#' Render a visual report workflow
#'
#' Materializes selected saved-run evidence and writes a DOCX and/or PDF
#' report. Browser-local AI text must already have been drafted into each
#' block; rendering never contacts an AI provider or the internet.
#'
#' @param design An [nm_report_design()].
#' @param workspace An [nm_workspace()] or path.
#' @param project Project id containing referenced runs.
#' @param directory Output directory.
#' @param name Output base name.
#' @param formats DOCX and/or PDF; defaults to the design formats.
#' @param manifest Write a JSON provenance manifest.
#' @return An `nm_report_bundle` containing generated file paths.
#' @export
nm_report_design_render <- function(design, workspace, project,
                                    directory = getwd(),
                                    name = "liberation-report",
                                    formats = design$formats,
                                    manifest = TRUE) {
  if (!inherits(design, "nm_report_design")) {
    .nm_stop("`design` must be an nm_report_design.")
  }
  workspace <- if (inherits(workspace, "nm_workspace")) workspace else nm_workspace(workspace)
  project <- .nm_workspace_component(project, "project id")
  if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE, showWarnings = FALSE)) {
    .nm_stop("Unable to create the report output directory.")
  }
  directory <- normalizePath(directory, winslash = "/", mustWork = TRUE)
  name <- gsub("[^A-Za-z0-9_.-]", "-", trimws(as.character(name)))
  if (!nzchar(name)) .nm_stop("Report `name` must contain a safe filename character.")
  formats <- unique(tolower(as.character(formats)))
  if (!length(formats) || length(setdiff(formats, c("docx", "pdf")))) {
    .nm_stop("Report formats must be DOCX and/or PDF.")
  }

  assets <- tempfile("liber-report-assets-")
  dir.create(assets)
  on.exit(unlink(assets, recursive = TRUE, force = TRUE), add = TRUE)
  content <- .nm_report_materialize(design, workspace, project, assets)
  files <- list(docx = NULL, pdf = NULL, json = NULL)
  if ("docx" %in% formats) {
    files$docx <- .nm_report_render_docx(
      design, content, file.path(directory, paste0(name, ".docx")), assets
    )
  }
  if ("pdf" %in% formats) {
    files$pdf <- .nm_report_render_pdf(
      design, content, file.path(directory, paste0(name, ".pdf"))
    )
  }
  if (isTRUE(manifest)) {
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      .nm_stop("Report manifests require the `jsonlite` package.")
    }
    files$json <- file.path(directory, paste0(name, ".json"))
    plain <- function(value) {
      if (is.data.frame(value) || is.atomic(value) || is.null(value)) return(value)
      if (is.list(value)) {
        output <- lapply(unclass(value), plain)
        names(output) <- names(value)
        return(output)
      }
      unclass(value)
    }
    manifest_value <- list(
      schema = "liber.report-bundle/1", version = 1L,
      generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC"),
      package_version = as.character(utils::packageVersion("LibeRation")),
      project = project, design = plain(design), formats = formats
    )
    jsonlite::write_json(manifest_value, files$json, auto_unbox = TRUE,
                         pretty = TRUE, null = "null", digits = NA)
  }
  output <- structure(c(files, list(
    design = design, project = project,
    generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")
  )), class = "nm_report_bundle")
  if (.Platform$OS.type != "windows") {
    existing <- unlist(files, use.names = FALSE)
    existing <- existing[!is.na(existing) & nzchar(existing) & file.exists(existing)]
    if (length(existing)) Sys.chmod(existing, mode = "0640")
  }
  output
}
