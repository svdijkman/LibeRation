.liber_gui_background_diagnostics <- function(fit, options) {
  selected <- unique(tolower(as.character(unlist(
    options$types %||% character()
  ))))
  selected <- intersect(selected, c(
    "vpc", "npc", "npde", "vpc_categorical", "vpc_count",
    "vpc_tte", "vpc_competing", "vpc_recurrent"
  ))
  if (!length(selected)) .nm_stop("Select at least one diagnostic.")
  nsim <- max(20L, as.integer(options$nsim %||% 200L))
  seed <- as.integer(options$seed %||% 20260713L)
  created <- list()
  if ("vpc" %in% selected) {
    created$vpc <- nm_vpc(
      fit, nsim = nsim, seed = seed,
      pc_correct = isTRUE(options$pcCorrect),
      stratify = options$stratify %||% NULL
    )
  }
  if ("npc" %in% selected) {
    created$npc <- nm_npc(fit, nsim = nsim, seed = seed)
  }
  if ("npde" %in% selected) {
    created$npde <- nm_npde(fit, nsim = nsim, seed = seed)
  }
  if ("vpc_categorical" %in% selected) {
    created$vpc_categorical <- nm_vpc_categorical(
      fit, outcome = as.character(options$categoricalOutcome %||% "DV"),
      nsim = nsim, seed = seed
    )
  }
  if ("vpc_count" %in% selected) {
    created$vpc_count <- nm_vpc_count(
      fit, outcome = as.character(options$countOutcome %||% "DV"),
      dvid = suppressWarnings(as.numeric(options$countDvid %||% NA_real_)),
      nsim = nsim, seed = seed
    )
  }
  if ("vpc_tte" %in% selected) {
    created$vpc_tte <- nm_vpc_tte(
      fit, event = as.character(options$tteEvent %||% "DV"),
      nsim = nsim, seed = seed
    )
  }
  if ("vpc_competing" %in% selected) {
    created$vpc_competing <- nm_vpc_competing(
      fit, event = as.character(options$tteEvent %||% "DV"),
      dvid = suppressWarnings(as.numeric(
        options$competingDvid %||% NA_real_
      )),
      nsim = nsim, seed = seed
    )
  }
  if ("vpc_recurrent" %in% selected) {
    created$vpc_recurrent <- nm_vpc_recurrent(
      fit, event = as.character(options$tteEvent %||% "DV"),
      dvid = suppressWarnings(as.numeric(
        options$recurrentDvid %||% NA_real_
      )),
      nsim = nsim, seed = seed
    )
  }
  created
}

.liber_gui_background_uncertainty <- function(fit, options) {
  selected <- unique(tolower(as.character(unlist(
    options$types %||% character()
  ))))
  selected <- intersect(selected, c("bootstrap", "profile"))
  if (!length(selected)) .nm_stop("Select bootstrap or profile likelihood.")
  created <- list()
  if ("bootstrap" %in% selected) {
    created$bootstrap <- nm_bootstrap(
      fit, n = as.integer(options$replicates %||% 100L),
      seed = as.integer(options$seed %||% 20260713L),
      level = as.numeric(options$level %||% 0.95),
      maxit = as.integer(options$maxit %||% 100L)
    )
  }
  if ("profile" %in% selected) {
    parameters <- trimws(unlist(strsplit(
      as.character(options$parameters %||% ""), "[,;[:space:]]+"
    )))
    parameters <- parameters[nzchar(parameters)]
    created$profile <- nm_profile(
      fit, parameters = if (length(parameters)) parameters else NULL,
      points = as.integer(options$points %||% 9L),
      span = as.numeric(options$span %||% 3),
      level = as.numeric(options$level %||% 0.95),
      maxit = as.integer(options$maxit %||% 100L)
    )
  }
  created
}

.liber_gui_background_scm <- function(fit, options) {
  lines <- strsplit(
    as.character(options$candidates %||% ""), "\r?\n", perl = TRUE
  )[[1L]]
  lines <- trimws(lines[nzchar(trimws(lines))])
  if (!length(lines)) {
    .nm_stop("Enter at least one parameter,covariate candidate.")
  }
  fields <- lapply(
    lines, function(line) trimws(strsplit(line, ",", fixed = TRUE)[[1L]])
  )
  item <- function(field, index, default = "") {
    if (length(field) >= index && nzchar(field[[index]])) {
      field[[index]]
    } else default
  }
  candidates <- do.call(rbind, lapply(fields, function(field) data.frame(
    parameter = item(field, 1L),
    covariate = item(field, 2L),
    form = item(field, 3L, "continuous"),
    reference = item(field, 4L, NA_character_),
    category = item(field, 5L, NA_character_),
    stringsAsFactors = FALSE
  )))
  nm_scm(
    fit, candidates,
    direction = as.character(options$direction %||% "both"),
    p_forward = as.numeric(options$pForward %||% 0.05),
    p_backward = as.numeric(options$pBackward %||% 0.01),
    max_steps = as.integer(options$maxSteps %||% 20L),
    maxit = as.integer(options$maxit %||% 100L)
  )
}

.liber_gui_compare_runs <- function(workspace_path, project, ids) {
  workspace <- nm_workspace(workspace_path)
  ids <- unique(as.character(ids))
  if (length(ids) != 2L) {
    .nm_stop("Select exactly two estimation runs to compare.")
  }
  entries <- lapply(ids, function(id) nm_project_load(workspace, project, id))
  fits <- lapply(entries, `[[`, "result")
  if (any(!vapply(fits, inherits, logical(1), "nm_fit"))) {
    .nm_stop("Both selected runs must contain estimation results.")
  }
  labels <- make.unique(vapply(
    entries, function(entry) entry$label, character(1)
  ))
  gof_frames <- lapply(fits, nm_gof)
  vectors <- lapply(fits, .liber_gui_parameter_values)
  parameter_names <- unique(unlist(lapply(vectors, names), use.names = FALSE))
  parameters <- data.frame(
    Parameter = parameter_names, stringsAsFactors = FALSE
  )
  for (index in seq_along(vectors)) {
    parameters[[labels[[index]]]] <- unname(
      vectors[[index]][parameter_names]
    )
    covariance <- fits[[index]]$covariance
    if (!is.null(covariance$se)) {
      parameters[[paste(labels[[index]], "SE")]] <- unname(
        covariance$se[parameter_names]
      )
      parameters[[paste(labels[[index]], "RSE")]] <- unname(
        covariance$rse[parameter_names]
      )
    }
  }
  gof_summary <- function(fit, gof) {
    observed <- gof$EVID == 0L & gof$MDV == 0L & is.finite(gof$DV)
    gof <- gof[observed, , drop = FALSE]
    n_parameters <- sum(!fit$model$THETAS$FIX) +
      sum(!fit$model$SIGMAS$FIX) + sum(!fit$model$OMEGAS$FIX)
    n_observations <- nrow(gof)
    c(
      OFV = fit$objective,
      AIC = fit$objective + 2 * n_parameters,
      BIC = fit$objective + log(max(1, n_observations)) * n_parameters,
      `Free parameters` = n_parameters,
      Observations = n_observations,
      `Population RMSE` = sqrt(mean(
        (gof$DV - gof$PRED)^2, na.rm = TRUE
      )),
      `Individual RMSE` = sqrt(mean(
        (gof$DV - gof$IPRED)^2, na.rm = TRUE
      )),
      `Mean WRES` = mean(gof$WRES, na.rm = TRUE),
      `SD WRES` = stats::sd(gof$WRES, na.rm = TRUE),
      `Mean IWRES` = mean(gof$IWRES, na.rm = TRUE),
      `SD IWRES` = stats::sd(gof$IWRES, na.rm = TRUE),
      `Mean CWRES` = mean(gof$CWRES, na.rm = TRUE),
      `SD CWRES` = stats::sd(gof$CWRES, na.rm = TRUE)
    )
  }
  summaries <- Map(gof_summary, fits, gof_frames)
  metric_names <- unique(unlist(lapply(summaries, names), use.names = FALSE))
  gof <- data.frame(Metric = metric_names, stringsAsFactors = FALSE)
  for (index in seq_along(summaries)) {
    gof[[labels[[index]]]] <- unname(summaries[[index]][metric_names])
  }
  runs <- data.frame(
    Run = labels,
    Method = vapply(fits, function(fit) fit$method, character(1)),
    Objective = vapply(fits, function(fit) fit$objective, numeric(1)),
    Convergence = vapply(fits, function(fit) fit$convergence, integer(1)),
    stringsAsFactors = FALSE
  )
  plots <- list(gof = unname(Map(function(label, frame) {
    list(
      label = label,
      fit = list(
        available = TRUE, gof_loaded = TRUE,
        gof = .liber_gui_rows(frame, 5000L)
      )
    )
  }, labels, gof_frames)))
  run_diagnostics <- lapply(ids, function(id) {
    nm_project_load_diagnostics(workspace, project, id)
  })
  for (kind in c(
    "vpc", "vpc_categorical", "vpc_count", "vpc_tte",
    "vpc_competing", "vpc_recurrent", "npde", "npc"
  )) {
    if (all(vapply(
      run_diagnostics, function(item) !is.null(item[[kind]]), logical(1)
    ))) {
      plots[[kind]] <- unname(Map(function(label, item) {
        list(label = label, result = .liber_gui_result(item[[kind]]))
      }, labels, run_diagnostics))
    }
  }
  structure(
    list(
      id = paste0(
        format(Sys.time(), "%Y%m%d%H%M%OS3"), "-",
        sample.int(999999L, 1L)
      ),
      parameters = parameters, gof = gof, runs = runs, plots = plots
    ),
    class = "liber_gui_comparison"
  )
}

.liber_gui_background_task <- function(operation, arguments) {
  operation <- match.arg(operation, c(
    "simulate", "estimate", "estimate_sequence", "diagnostics",
    "uncertainty", "scm", "comparison", "gof", "hmm", "kalman",
    "report", "report_design", "ai_context", "report_ai_context"
  ))
  if (!is.list(arguments)) {
    .nm_stop("Background GUI task arguments must be a list.")
  }
  switch(
    operation,
    simulate = do.call(
      nm_simulate,
      c(list(model = arguments$model, data = arguments$data), arguments$args)
    ),
    estimate = do.call(
      nm_est,
      c(list(model = arguments$model, data = arguments$data), arguments$args)
    ),
    estimate_sequence = do.call(
      nm_est_sequence,
      c(list(model = arguments$model, data = arguments$data), arguments$args)
    ),
    diagnostics = .liber_gui_background_diagnostics(
      arguments$fit, arguments$options
    ),
    uncertainty = .liber_gui_background_uncertainty(
      arguments$fit, arguments$options
    ),
    scm = .liber_gui_background_scm(arguments$fit, arguments$options),
    comparison = .liber_gui_compare_runs(
      arguments$workspace_path, arguments$project, arguments$ids
    ),
    gof = nm_gof(arguments$fit),
    hmm = nm_hmm_decode(arguments$fit, method = "all"),
    kalman = nm_kalman_decode(arguments$fit, type = "individual"),
    report = nm_report(
      arguments$fit, arguments$file,
      sections = arguments$sections, vpc = arguments$vpc
    ),
    report_design = nm_report_design_render(
      arguments$design,
      nm_workspace(arguments$workspace_path),
      arguments$project,
      directory = arguments$directory,
      name = arguments$name,
      formats = arguments$design$formats
    ),
    ai_context = .liber_gui_ai_context(
      nm_workspace(arguments$workspace_path),
      arguments$project,
      selected_run = arguments$selected_run,
      max_runs = arguments$max_runs,
      detail = arguments$detail,
      question = arguments$question
    ),
    report_ai_context = .liber_gui_report_ai_context(
      nm_workspace(arguments$workspace_path),
      arguments$project,
      arguments$run_ids
    )
  )
}
