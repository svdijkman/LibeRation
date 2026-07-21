.liber_gui_rows <- function(data, limit = Inf) {
  data <- as.data.frame(data %||% data.frame(), stringsAsFactors = FALSE)
  if (!nrow(data)) return(list())
  if (is.finite(limit)) data <- utils::head(data, as.integer(limit))
  unname(lapply(seq_len(nrow(data)), function(i) {
    row <- as.list(data[i, , drop = FALSE])
    lapply(row, function(value) {
      if (!length(value) || (length(value) == 1L && is.na(value))) NULL else unname(value)
    })
  }))
}

.liber_gui_matrix_rows <- function(value) {
  if (!is.matrix(value) || !nrow(value) || !ncol(value)) return(list())
  frame <- data.frame(
    Parameter = rownames(value) %||% as.character(seq_len(nrow(value))),
    as.data.frame(value, check.names = FALSE),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  .liber_gui_rows(frame)
}

.liber_gui_covariance <- function(covariance) {
  if (is.null(covariance)) return(list(requested = FALSE, status = "not_requested"))
  status <- as.character(covariance$status %||%
    if (!is.null(covariance$covariance)) "completed" else "failed")
  list(
    requested = TRUE, status = status,
    type = as.character(covariance$type %||% ""),
    error = as.character(covariance$error %||% ""),
    condition = covariance$condition %||% NULL,
    regularization = covariance$regularization %||% NULL,
    min_eigenvalue = if (length(covariance$eigenvalues)) min(covariance$eigenvalues) else NULL,
    max_eigenvalue = if (length(covariance$eigenvalues)) max(covariance$eigenvalues) else NULL,
    samples = covariance$samples %||% NULL,
    actual_samples = covariance$actual_samples %||% NULL,
    sampling = covariance$sampling %||% NULL,
    quadrature_order = covariance$quadrature_order %||% NULL,
    quadrature_level = covariance$quadrature_level %||% NULL,
    quadrature_grid = covariance$quadrature_grid %||% NULL,
    objective_backend = covariance$objective_backend %||% NULL,
    seed = covariance$seed %||% NULL,
    covariance = .liber_gui_matrix_rows(covariance$covariance),
    correlation = .liber_gui_matrix_rows(covariance$correlation)
  )
}

.liber_gui_posterior <- function(fit) {
  if (!inherits(fit, "nm_fit") ||
      !fit$method %in% c("BAYES", "HMC", "NUTS") || is.null(fit$posterior)) {
    return(list(available = FALSE, parameters = list()))
  }
  names <- .nm_parameter_names(fit$theta, fit$sigma, fit$omega)
  population <- fit$posterior$population
  if (is.null(population)) {
    population <- list(
      mean = fit$posterior$mean[names],
      sd = fit$posterior$sd[names],
      quantile = fit$posterior$quantile[, names, drop = FALSE]
    )
    if (!is.null(fit$chain) && nrow(fit$chain) > 1L) {
      population$covariance <- stats::cov(fit$chain[, names, drop = FALSE])
      population$correlation <- stats::cor(fit$chain[, names, drop = FALSE])
    }
  }
  quantile <- population$quantile
  rows <- unname(lapply(seq_along(names), function(index) {
    name <- names[[index]]
    mean <- unname(population$mean[[name]])
    sd <- unname(population$sd[[name]])
    list(
      name = name, mean = mean, posterior_sd = sd,
      posterior_cv = 100 * sd / max(abs(mean), 1e-12),
      median = unname(quantile[2L, name]),
      lower_95 = unname(quantile[1L, name]),
      upper_95 = unname(quantile[3L, name]),
      rhat = unname(population$rhat[[name]] %||% NA_real_),
      ess = unname(population$ess[[name]] %||% NA_real_)
    )
  }))
  list(
    available = TRUE, parameters = rows,
    samples = fit$diagnostics$n_sample %||% if (!is.null(fit$chain)) nrow(fit$chain) else NULL,
    burn = fit$diagnostics$n_burn %||% NULL,
    thin = fit$diagnostics$n_thin %||% NULL,
    outer_acceptance = fit$diagnostics$outer_acceptance %||% NULL,
    eta_acceptance = fit$diagnostics$eta_acceptance %||% NULL,
    divergences = fit$diagnostics$divergences %||% NULL,
    mean_acceptance = fit$diagnostics$mean_acceptance %||% NULL,
    max_depth_hits = fit$diagnostics$max_depth_hits %||% NULL,
    covariance = .liber_gui_matrix_rows(population$covariance),
    correlation = .liber_gui_matrix_rows(population$correlation)
  )
}

.liber_gui_nonparametric <- function(fit) {
  distribution <- fit$nonparametric
  if (!inherits(fit, "nm_fit") || is.null(distribution) ||
      !fit$method %in% c("NPML", "NPAG")) {
    return(list(available = FALSE, supports = list()))
  }
  supports <- as.data.frame(distribution$supports, stringsAsFactors = FALSE)
  supports <- data.frame(
    support = seq_len(nrow(supports)), weight = distribution$weights,
    supports, check.names = FALSE
  )
  list(
    available = TRUE, supports = .liber_gui_rows(supports),
    support_count = nrow(supports),
    log_likelihood = distribution$log_likelihood,
    interpretation = distribution$interpretation %||% ""
  )
}

.liber_gui_model <- function(model, output_catalog = NULL) {
  if (is.null(model)) {
    return(list(
      loaded = FALSE, name = "Untitled model", advan = "--", trans = "--",
      solver = "auto", language = "R", pred = "# Define a model to begin",
      error = "Y = F", des = "", alg = "", compartments = list(), flows = list(),
      theta = list(), omega = list(), sigma = list(), n_theta = 0L,
      n_eta = 0L, n_sigma = 0L, n_state = 0L, input = character(),
      output = character(), outputs = list(),
      diagram = list(available = FALSE, graph = NULL, code_changed = FALSE)
    ))
  }
  if (inherits(model, "NMEngine")) model <- model$model
  if (!inherits(model, "nm_model")) .nm_stop("`model` must be an nm_model or NMEngine.")
  compartments <- model$GRAPH$compartments %||% data.frame()
  output_catalog <- output_catalog %||% nm_model_outputs(model)
  diagram <- nm_model_diagram_get(model)
  diagram_payload <- if (is.null(diagram)) {
    list(available = FALSE, graph = NULL, code_changed = FALSE)
  } else {
    generated <- diagram$generated %||% list()
    list(
      available = TRUE,
      graph = list(
        schema = diagram$schema, version = diagram$version,
        title = diagram$title, advan = diagram$advan,
        residual = diagram$residual, covariates = diagram$covariates,
        compartments = .liber_gui_rows(diagram$compartments),
        flows = .liber_gui_rows(diagram$flows),
        parameters = .liber_gui_rows(diagram$parameters),
        generated = generated
      ),
      code_changed = nzchar(as.character(generated$DES %||% "")) &&
        !identical(trimws(model$DES), trimws(as.character(generated$DES)))
    )
  }
  list(
    loaded = TRUE,
    name = attr(model, "name", exact = TRUE) %||% paste0("ADVAN", model$ADVAN, " model"),
    advan = model$ADVAN,
    trans = model$TRANS,
    solver = model$SOLVER,
    language = model$LANGUAGE,
    input = model$INPUT,
    output = model$OUTPUT %||% character(),
    outputs = .liber_gui_rows(output_catalog),
    n_state = model$n_state,
    dose_cmp = model$DOSECMP,
    obs_cmp = model$OBSCMP,
    pred = model$PRED,
    error = model$ERROR,
    des = model$DES,
    alg = model$ALG %||% "",
    experimental = if (isTRUE(model$EXPERIMENTAL$enabled)) list(
      enabled = TRUE, strict = isTRUE(model$EXPERIMENTAL$strict),
      label = model$EXPERIMENTAL$label %||% "",
      features = unname(model$EXPERIMENTAL$features %||% character())
    ) else NULL,
    dde = if (is.null(model$DDE_CONFIG)) NULL else list(
      step = model$DDE_CONFIG$step, interpolation = model$DDE_CONFIG$interpolation,
      lag_count = length(model$DDE_CONFIG$lags),
      delays = unname(vapply(model$DDE_CONFIG$lags, `[[`, character(1), "delay"))
    ),
    dae = if (is.null(model$DAE_CONFIG)) NULL else list(
      variables = unname(model$DAE_CONFIG$variables),
      tolerance = model$DAE_CONFIG$tolerance, maxit = model$DAE_CONFIG$maxit,
      sparse = !is.null(model$DAE_CONFIG$sparsity)
    ),
    components = unname(lapply(model$COMPONENTS %||% list(), function(value) list(
      name = value$name, type = value$type, scope = value$scope %||% "pred",
      outputs = unname(value$outputs), hash = value$hash
    ))),
    qsp = if (is.null(model$QSP_SYSTEM)) NULL else list(
      species = unname(model$QSP_SYSTEM$species),
      reactions = ncol(model$QSP_SYSTEM$stoichiometry)
    ),
    compartments = if (nrow(compartments)) {
      lapply(seq_len(nrow(compartments)), function(i) list(
        id = compartments$id[[i]], name = compartments$name[[i]]
      ))
    } else list(),
    flows = if (inherits(model$GRAPH, "nm_matrix_model") && nrow(model$GRAPH$flows)) {
      lapply(seq_len(nrow(model$GRAPH$flows)), function(i) as.list(model$GRAPH$flows[i, , drop = FALSE]))
    } else list(),
    theta = lapply(seq_len(nrow(model$THETAS)), function(i) as.list(model$THETAS[i, , drop = FALSE])),
    omega = lapply(seq_len(nrow(model$OMEGAS)), function(i) as.list(model$OMEGAS[i, , drop = FALSE])),
    sigma = lapply(seq_len(nrow(model$SIGMAS)), function(i) as.list(model$SIGMAS[i, , drop = FALSE])),
    omega_structure = model$LIK_CONFIG$omega %||% "diagonal",
    likelihood_type = model$LIK_CONFIG$error %||% model$ERROR_TYPE %||% "none",
    outcomes = if (is.null(model$OUTCOMES)) list() else lapply(model$OUTCOMES, function(value) list(
      name = value$name, family = value$family, dvid = value$dvid,
      prediction = value$prediction,
      categories = unname(value$categories %||% numeric()),
      generated_error = isTRUE(model$outcome_error_generated)
    )),
    hmm = if (is.null(model$HMM_CONFIG)) NULL else list(
      states = unname((attr(model$HMM_CONFIG, "semi_markov", exact = TRUE) %||%
        list(states = model$HMM_CONFIG$states))$states),
      model_type = if (inherits(model$HMM_CONFIG, "nm_factorial_hmm_config")) {
        "factorial hidden Markov"
      } else if (inherits(model$HMM_CONFIG, "nm_hsmm_config")) {
        "hidden semi-Markov"
      } else if (identical(model$HMM_CONFIG$transition_type %||% "discrete", "continuous")) {
        "continuous-time hidden Markov"
      } else "hidden Markov",
      transition_type = model$HMM_CONFIG$transition_type %||% "discrete",
      initial_scale = model$HMM_CONFIG$initial_scale,
      transition_scale = model$HMM_CONFIG$transition_scale,
      rate_scale = model$HMM_CONFIG$rate_scale,
      emission_scale = model$HMM_CONFIG$emission_scale,
      by_dvid = isTRUE(model$HMM_CONFIG$by_dvid)
    ),
    kalman = if (is.null(model$KALMAN_CONFIG)) NULL else list(
      states = unname(model$KALMAN_CONFIG$states),
      baseline = model$KALMAN_CONFIG$baseline,
      by_dvid = isTRUE(model$KALMAN_CONFIG$by_dvid),
      filter = model$KALMAN_CONFIG$filter %||% "linear",
      dynamics = model$KALMAN_CONFIG$dynamics %||% "discrete",
      model_type = if (inherits(model$KALMAN_CONFIG, "nm_switching_state_space_config")) {
        "switching state-space"
      } else "state-space",
      regimes = unname(model$KALMAN_CONFIG$switching$regimes %||% character())
    ),
    random_effects = if (is.null(model$RE_CONFIG)) NULL else list(
      cluster = model$RE_CONFIG$cluster %||% "auto",
      blocks = unname(lapply(model$RE_CONFIG$blocks, function(block) list(
        name = block$name, column = block$column, etas = unname(block$etas)
      )))
    ),
    priors = .liber_gui_rows(model$LIK_CONFIG$priors %||% data.frame()),
    n_theta = nrow(model$THETAS),
    n_eta = model$n_eta,
    n_sigma = nrow(model$SIGMAS),
    diagram = diagram_payload
  )
}

.liber_gui_data <- function(data, include_rows = TRUE, run_output = NULL) {
  if (is.null(data)) {
    return(list(loaded = FALSE, records = 0L, subjects = 0L, observations = 0L,
                columns = character(), numeric_columns = character(),
                categorical_columns = character(), preview = list(),
                preview_all = list(), plot_rows = list(), issues = list(),
                payload_loaded = FALSE))
  }
  data <- if (inherits(data, "nm_dataset")) data else nm_dataset(data)
  frame <- as.data.frame(data)
  if (isTRUE(include_rows) && !is.null(run_output)) {
    run_output <- as.data.frame(run_output, stringsAsFactors = FALSE)
    output_names <- setdiff(names(run_output), c(".ROW", names(frame)))
    if (length(output_names)) {
      rows <- if (".ROW" %in% names(run_output)) {
        suppressWarnings(as.integer(run_output$.ROW))
      } else if (nrow(run_output) == nrow(frame)) {
        seq_len(nrow(frame))
      } else integer()
      valid <- is.finite(rows) & rows >= 1L & rows <= nrow(frame)
      for (name in output_names) {
        value <- rep(NA, nrow(frame))
        value[rows[valid]] <- run_output[[name]][valid]
        frame[[name]] <- value
      }
    }
  }
  visible <- setdiff(names(frame), grep("^\\.", names(frame), value = TRUE))
  plot_frame <- frame[, visible, drop = FALSE]
  numeric_columns <- visible[vapply(plot_frame, is.numeric, logical(1))]
  list(
    loaded = TRUE,
    name = attr(data, "name", exact = TRUE) %||% "Current dataset",
    records = nrow(data),
    subjects = length(unique(data$.ID_INDEX)),
    observations = sum(data$EVID == 0L & data$MDV == 0L),
    columns = visible,
    numeric_columns = numeric_columns,
    categorical_columns = setdiff(visible, numeric_columns),
    preview = if (isTRUE(include_rows)) .liber_gui_rows(plot_frame, 12L) else list(),
    preview_all = if (isTRUE(include_rows)) .liber_gui_rows(plot_frame, 100L) else list(),
    plot_rows = if (isTRUE(include_rows)) .liber_gui_rows(plot_frame, 5000L) else list(),
    payload_loaded = isTRUE(include_rows),
    issues = list()
  )
}

.liber_client_settings_path <- function(workspace) {
  file.path(.nm_workspace_path(workspace), ".liberation", "client-settings.rds")
}

.liber_ai_default_help_model <- function() {
  "Qwen2.5-Coder-3B-Instruct-q4f16_1-MLC"
}

.liber_ai_default_report_model <- function() {
  "Qwen2.5-7B-Instruct-q4f16_1-MLC"
}

.liber_ai_default_model <- function() {
  .liber_ai_default_help_model()
}

.liber_ai_context_setting <- function(value, default = "auto") {
  value <- as.character(value %||% default)[[1L]]
  if (identical(tolower(value), "auto")) return("auto")
  numeric_value <- suppressWarnings(as.integer(value))
  if (is.na(numeric_value) || numeric_value < 1024L || numeric_value > 16384L) {
    return(default)
  }
  as.character(numeric_value)
}

.liber_ai_models <- function() {
  list(
    list(
      id = "SmolLM2-360M-Instruct-q4f16_1-MLC",
      label = "SmolLM2 360M - quickest (~0.4 GB)", tier = "minimal",
      vram_mb = 376L,
      description = "Very fast, but only suitable for simple workflow questions."
    ),
    list(
      id = "Qwen2.5-0.5B-Instruct-q4f16_1-MLC",
      label = "Qwen 2.5 0.5B - quick (~0.9 GB)", tier = "minimal",
      vram_mb = 945L,
      description = "Low-memory option; factual reliability is limited."
    ),
    list(
      id = "Llama-3.2-1B-Instruct-q4f16_1-MLC",
      label = "Llama 3.2 1B - light (~0.9 GB)", tier = "light",
      vram_mb = 879L,
      description = "A lightweight general assistant for short questions."
    ),
    list(
      id = "Qwen2.5-1.5B-Instruct-q4f16_1-MLC",
      label = "Qwen 2.5 1.5B - light (~1.6 GB)", tier = "light",
      vram_mb = 1630L,
      description = "Better instruction following than the minimal models."
    ),
    list(
      id = "SmolLM2-1.7B-Instruct-q4f16_1-MLC",
      label = "SmolLM2 1.7B - legacy (~1.8 GB)", tier = "light",
      vram_mb = 1774L,
      description = "Original LibeRation default; fast, but more prone to unsupported claims."
    ),
    list(
      id = "gemma-2-2b-it-q4f16_1-MLC",
      label = "Gemma 2 2B - balanced (~1.9 GB)", tier = "balanced",
      vram_mb = 1895L,
      description = "Compact general-purpose alternative."
    ),
    list(
      id = "Qwen2.5-3B-Instruct-q4f16_1-MLC",
      label = "Qwen 2.5 3B - balanced (~2.5 GB)", tier = "balanced",
      vram_mb = 2505L,
      description = "Recommended balance of reliability, speed, and memory use."
    ),
    list(
      id = "Qwen2.5-Coder-3B-Instruct-q4f16_1-MLC",
      label = "Qwen 2.5 Coder 3B - recommended for Help (~2.5 GB)", tier = "recommended",
      vram_mb = 2505L,
      description = "Prefer for model-code and syntax questions."
    ),
    list(
      id = "Llama-3.2-3B-Instruct-q4f16_1-MLC",
      label = "Llama 3.2 3B - balanced (~2.3 GB)", tier = "balanced",
      vram_mb = 2264L,
      description = "Strong compact general assistant and alternative to Qwen 3B."
    ),
    list(
      id = "Qwen2.5-7B-Instruct-q4f16_1-MLC",
      label = "Qwen 2.5 7B - recommended for Reports (~5.1 GB)", tier = "quality",
      vram_mb = 5107L,
      description = "Best practical quality option when the GPU has sufficient free memory."
    ),
    list(
      id = "Llama-3.1-8B-Instruct-q4f16_1-MLC-1k",
      label = "Llama 3.1 8B - higher quality, 1K context (~4.6 GB)", tier = "quality",
      vram_mb = 4598L,
      description = "Higher-capacity model with a shorter context window."
    ),
    list(
      id = "gemma-2-9b-it-q4f16_1-MLC",
      label = "Gemma 2 9B - largest (~6.4 GB)", tier = "large",
      vram_mb = 6422L,
      description = "Largest option; requires substantial free GPU memory and loads more slowly."
    )
  )
}

.liber_client_settings_read <- function(workspace) {
  path <- .liber_client_settings_path(workspace)
  defaults <- list(version = 5L, selected_queue = "local", remotes = list(),
                   pending_jobs = list(), ai = list(
                     activated = FALSE, consented = FALSE,
                     help_model = .liber_ai_default_help_model(),
                     report_model = .liber_ai_default_report_model(),
                     help_context = "auto", report_context = "auto",
                     model = .liber_ai_default_help_model()
                   ))
  if (!file.exists(path)) return(defaults)
  value <- tryCatch(readRDS(path), error = function(error) NULL)
  if (!is.list(value)) return(defaults)
  remotes <- value$remotes
  if (!is.list(remotes)) remotes <- list()
  remotes <- Filter(function(item) {
    is.list(item) && length(item$url) == 1L && !is.na(item$url) && nzchar(item$url)
  }, remotes)
  pending_jobs <- value$pending_jobs
  if (!is.list(pending_jobs)) pending_jobs <- list()
  pending_jobs <- Filter(function(item) {
    is.list(item) && length(item$job_id) == 1L && !is.na(item$job_id) &&
      nzchar(item$job_id) && length(item$queue_id) == 1L &&
      !is.na(item$queue_id) && nzchar(item$queue_id)
  }, pending_jobs)
  ai <- value$ai
  if (!is.list(ai)) ai <- defaults$ai
  legacy_model <- as.character(ai$model %||% "")[[1L]]
  help_model <- as.character(ai$help_model %||% legacy_model %||%
    defaults$ai$help_model)[[1L]]
  if (!nzchar(help_model)) help_model <- defaults$ai$help_model
  report_model <- as.character(ai$report_model %||%
    defaults$ai$report_model)[[1L]]
  if (!nzchar(report_model)) report_model <- defaults$ai$report_model
  help_context <- .liber_ai_context_setting(
    ai$help_context, defaults$ai$help_context
  )
  report_context <- .liber_ai_context_setting(
    ai$report_context, defaults$ai$report_context
  )
  list(
    version = 5L,
    selected_queue = as.character(value$selected_queue %||% "local")[[1L]],
    remotes = remotes, pending_jobs = pending_jobs,
    ai = list(
      activated = isTRUE(ai$activated), consented = isTRUE(ai$consented),
      help_model = help_model, report_model = report_model,
      help_context = help_context, report_context = report_context,
      model = help_model
    )
  )
}

.liber_client_settings_write <- function(workspace, selected_queue = "local",
                                         remotes = list(), pending_jobs = list(),
                                         ai = list()) {
  path <- .liber_client_settings_path(workspace)
  .nm_workspace_atomic_save(
    list(version = 5L, selected_queue = as.character(selected_queue)[[1L]],
         remotes = remotes, pending_jobs = pending_jobs,
         ai = list(
           activated = isTRUE(ai$activated), consented = isTRUE(ai$consented),
           help_model = as.character(ai$help_model %||% ai$model %||%
             .liber_ai_default_help_model())[[1L]],
           report_model = as.character(ai$report_model %||%
             .liber_ai_default_report_model())[[1L]],
           help_context = .liber_ai_context_setting(ai$help_context),
           report_context = .liber_ai_context_setting(ai$report_context),
           model = as.character(ai$help_model %||% ai$model %||%
             .liber_ai_default_help_model())[[1L]]
         )),
    path
  )
  try(Sys.chmod(path, mode = "0600"), silent = TRUE)
  invisible(path)
}

.liber_default_workspace <- function() {
  home <- path.expand("~")
  if (.Platform$OS.type == "windows") {
    profile <- Sys.getenv("USERPROFILE", unset = home)
    return(file.path(profile, "Documents", "LibeR", "workspace"))
  }
  file.path(home, "LibeR", "workspace")
}

.liber_gui_library <- function() {
  if (!requireNamespace("LibeRary", quietly = TRUE)) {
    return(list(available = FALSE, entries = list(),
                message = "Install LibeRary to browse and import pharmacometric models."))
  }
  tryCatch({
    root <- LibeRary::library_catalog_root()
    entries <- LibeRary::library_list(root = root)
    list(available = TRUE, root = root, entries = .liber_gui_rows(entries),
         ingest_available = TRUE, message = if (nrow(entries)) "" else "The catalogue is empty.")
  }, error = function(error) list(available = FALSE, entries = list(),
                                  message = conditionMessage(error)))
}

.liber_gui_duration <- function(seconds) {
  seconds <- suppressWarnings(as.numeric(seconds))
  if (length(seconds) != 1L || is.na(seconds) || !is.finite(seconds) || seconds < 0) {
    return("")
  }
  if (seconds < 1) return(sprintf("%.3f seconds", seconds))
  if (seconds < 60) return(sprintf("%.2f seconds", seconds))
  hours <- floor(seconds / 3600)
  minutes <- floor((seconds %% 3600) / 60)
  remainder <- seconds %% 60
  if (hours > 0) {
    return(sprintf("%d h %02d min %04.1f s", hours, minutes, remainder))
  }
  sprintf("%d min %04.1f s", minutes, remainder)
}

.liber_gui_workspace <- function(workspace = NULL, current = NULL, version = NULL,
                                 run = NULL, pending_jobs = list(), jobs = NULL,
                                 queue_id = "local") {
  if (is.null(workspace)) {
    return(list(enabled = FALSE, current = NULL, current_version = NULL,
                current_snapshot = NULL, current_run = NULL,
                current_result_type = NULL, projects = list(), versions = list()))
  }
  workspace <- if (inherits(workspace, "nm_workspace")) workspace else nm_workspace(workspace)
  projects <- nm_project_list(workspace)
  records <- if (is.null(current)) data.frame() else nm_project_list(workspace, current)
  versions <- if (nrow(records)) records[records$entry_type == "version", , drop = FALSE] else records
  if (nrow(versions)) versions$version <- rev(seq_len(nrow(versions)))
  jobs <- as.data.frame(jobs %||% data.frame(), stringsAsFactors = FALSE)
  version_rows <- if (!nrow(versions)) list() else unname(lapply(seq_len(nrow(versions)), function(index) {
    item <- as.list(versions[index, , drop = FALSE])
    children <- records[
      records$entry_type == "run" & records$parent_id == versions$id[[index]], , drop = FALSE
    ]
    if (nrow(children)) children <- children[order(children$run_number), , drop = FALSE]
    item$runs <- .liber_gui_rows(children)
    pending <- Filter(function(context) {
      is.list(context) && identical(as.character(context$project %||% ""), as.character(current)) &&
        identical(as.character(context$version %||% ""), as.character(versions$id[[index]])) &&
        !nzchar(as.character(context$run_id %||% ""))
    }, pending_jobs %||% list())
    if (length(pending)) {
      virtual <- unname(lapply(pending, function(context) {
        status <- as.character(context$status %||% "queued")
        if (identical(as.character(context$queue_id), as.character(queue_id)) &&
            nrow(jobs) && all(c("id", "status") %in% names(jobs))) {
          match <- which(as.character(jobs$id) == as.character(context$job_id))
          if (length(match)) status <- as.character(jobs$status[[match[[1L]]]])
        }
        list(
          id = paste0("queue:", context$queue_id, ":", context$job_id),
          label = as.character(context$label %||% "Queued model run"),
          result_type = if (identical(context$type, "simulate")) "simulation" else "estimation",
          method = as.character(context$method %||% context$type %||% "job"),
          queued_job = TRUE, job_id = as.character(context$job_id),
          queue_id = as.character(context$queue_id), job_status = status,
          run_number = NULL
        )
      }))
      item$runs <- c(item$runs, virtual)
    }
    item
  }))
  selected <- if (!is.null(run) && nrow(records)) {
    records[records$id == run & records$entry_type == "run", , drop = FALSE]
  } else records[0L, , drop = FALSE]
  list(
    enabled = TRUE, path = workspace$path, current = current,
    current_version = version, current_snapshot = version, current_run = run,
    current_result_type = if (nrow(selected)) selected$result_type[[1L]] else NULL,
    projects = .liber_gui_rows(projects), versions = version_rows
  )
}

.liber_gui_parameter_values <- function(fit) {
  values <- c(fit$theta, fit$omega, fit$sigma)
  names(values) <- c(
    .nm_numbered_names("THETA", length(fit$theta)),
    .nm_numbered_names("OMEGA", length(fit$omega)),
    .nm_numbered_names("SIGMA", length(fit$sigma))
  )
  values
}

.liber_gui_fit <- function(fit, include_gof = TRUE, gof = NULL) {
  if (!inherits(fit, "nm_fit")) {
    return(list(available = FALSE, parameters = list(), gof = list(),
                gof_loaded = FALSE, run_info = list()))
  }
  etab <- nm_etab(fit)
  parameters <- .liber_gui_parameter_values(fit)
  covariance <- fit$covariance
  posterior <- .liber_gui_posterior(fit)
  nonparametric <- .liber_gui_nonparametric(fit)
  parameter_rows <- unname(lapply(seq_along(parameters), function(i) {
    name <- names(parameters)[[i]]
    row <- list(name = name, value = unname(parameters[[i]]))
    if (!is.null(covariance$se) && name %in% names(covariance$se)) {
      row$se <- unname(covariance$se[[name]])
      row$rse <- unname(covariance$rse[[name]])
    }
    if (isTRUE(posterior$available)) {
      posterior_names <- vapply(posterior$parameters, `[[`, character(1), "name")
      posterior_index <- match(name, posterior_names)
      if (!is.na(posterior_index)) {
        posterior_row <- posterior$parameters[[posterior_index]]
        row$posterior_sd <- posterior_row$posterior_sd
        row$posterior_cv <- posterior_row$posterior_cv
        row$median <- posterior_row$median
        row$lower_95 <- posterior_row$lower_95
        row$upper_95 <- posterior_row$upper_95
      }
    }
    row
  }))
  if (isTRUE(include_gof) && is.null(gof)) gof <- nm_gof(fit)
  if (!isTRUE(include_gof)) gof <- NULL
  timing <- fit$timing %||% list()
  iterations <- suppressWarnings(as.integer(fit$iterations %||% NA_integer_))
  if (!length(iterations) || is.na(iterations)) {
    iterations <- suppressWarnings(as.integer(fit$evaluations[["gradient"]] %||% NA_integer_))
  }
  if (!length(iterations) || is.na(iterations)) {
    iterations <- suppressWarnings(as.integer(fit$evaluations[["function"]] %||% NA_integer_))
  }
  objective_evaluations <- suppressWarnings(as.integer(
    fit$objective_evaluations %||% fit$evaluations[["function"]] %||% NA_integer_
  ))
  list(
    available = TRUE, method = .nm_fit_method_label(fit), final_method = fit$method,
    method_sequence = as.character(fit$method_sequence %||% fit$method),
    objective = fit$objective,
    convergence = fit$convergence,
    parameters = parameter_rows,
    covariance = .liber_gui_covariance(covariance),
    posterior = posterior,
    nonparametric = nonparametric,
    shrinkage = unname(lapply(seq_along(etab$shrinkage), function(i) {
      list(name = names(etab$shrinkage)[[i]], value = etab$shrinkage[[i]])
    })),
    gof = if (is.null(gof)) list() else .liber_gui_rows(gof, 5000L),
    gof_loaded = !is.null(gof),
    run_info = list(
      Method = fit$method, Objective = fit$objective,
      `Convergence code` = fit$convergence,
      `Iterations to convergence` = if (length(iterations) && !is.na(iterations)) iterations else "",
      `Objective evaluations` = if (length(objective_evaluations) &&
        !is.na(objective_evaluations)) objective_evaluations else "",
      `GQ grid` = if (fit$method == "GQ")
        as.character(fit$diagnostics$quadrature_grid %||% "tensor") else "",
      `GQ order / level` = if (fit$method == "GQ") {
        if (identical(fit$diagnostics$quadrature_grid, "smolyak")) {
          paste("level", fit$diagnostics$quadrature_level %||% "")
        } else paste("order", fit$diagnostics$quadrature_order %||% "")
      } else "",
      `GQ points per subject` = if (fit$method == "GQ")
        fit$diagnostics$quadrature_points %||% "" else "",
      `GQ minimum cancellation ratio` = if (fit$method == "GQ" &&
        length(fit$diagnostics$quadrature_cancellation_ratio)) {
        ratios <- fit$diagnostics$quadrature_cancellation_ratio
        if (any(is.finite(ratios))) min(ratios[is.finite(ratios)]) else ""
      } else "",
      `Model fit time` = .liber_gui_duration(timing$model_fit_seconds),
      `Covariance step` = if (is.null(covariance)) "Not requested" else
        if (identical(covariance$status, "failed")) "Failed" else "Completed",
      `Covariance method` = if (is.null(covariance)) "" else
        toupper(as.character(covariance$type %||% "")),
      `Covariance integration` = if (is.null(covariance)) "" else
        as.character(covariance$sampling %||% ""),
      `Covariance integration points` = if (is.null(covariance)) "" else
        covariance$actual_samples %||% "",
      `Covariance step time` = if (is.null(covariance)) "" else
        .liber_gui_duration(timing$covariance_seconds),
      `Total estimation time` = .liber_gui_duration(timing$total_seconds),
      `Posterior samples` = if (isTRUE(posterior$available)) posterior$samples else "",
      `Nonparametric support points` = if (isTRUE(nonparametric$available))
        nonparametric$support_count else "",
      Subjects = length(unique(fit$data$.ID_INDEX)), Records = nrow(fit$data),
      ADVAN = fit$model$ADVAN, TRANS = fit$model$TRANS,
      Solver = fit$model$SOLVER, Language = fit$model$LANGUAGE
    )
  )
}

.liber_gui_ai_context <- function(workspace, project, selected_run = NULL,
                                  max_runs = 20L, max_parameters = 24L,
                                  detail = c("results", "index"), run_ids = NULL) {
  detail <- match.arg(detail)
  unavailable <- function(message = "No saved project results are available.") {
    list(
      available = FALSE, project = as.character(project %||% ""),
      project_name = "", request_id = "", message = message,
      scope = detail, run_count = 0L, included_runs = 0L,
      omitted_runs = 0L, runs = list()
    )
  }
  if (is.null(workspace) || !nzchar(as.character(project %||% ""))) {
    return(unavailable("No project is selected."))
  }
  workspace <- if (inherits(workspace, "nm_workspace")) workspace else nm_workspace(workspace)
  project <- as.character(project)[[1L]]
  records <- nm_project_list(workspace, project)
  if (!nrow(records)) return(unavailable("The selected project has no saved model versions or runs."))
  projects <- nm_project_list(workspace)
  project_index <- match(project, projects$id)
  project_name <- if (!is.na(project_index)) {
    as.character(projects$name[[project_index]] %||% project)
  } else project
  versions <- records[records$entry_type == "version", , drop = FALSE]
  version_labels <- stats::setNames(as.character(versions$label), as.character(versions$id))
  runs <- records[
    records$entry_type == "run" & records$has_result %in% TRUE &
      records$result_type %in% c("estimation", "simulation"), , drop = FALSE
  ]
  requested_runs <- unique(as.character(run_ids %||% character()))
  if (length(requested_runs)) {
    missing_runs <- setdiff(requested_runs, as.character(runs$id))
    if (length(missing_runs)) {
      .nm_stop("Selected report run(s) are unavailable or incomplete: ",
               paste(missing_runs, collapse = ", "), ".")
    }
    runs <- runs[match(requested_runs, as.character(runs$id)), , drop = FALSE]
  }
  if (!nrow(runs)) {
    result <- unavailable("The selected project has no completed estimation or simulation runs.")
    result$project_name <- project_name
    return(result)
  }
  order_index <- if (length(requested_runs)) seq_len(nrow(runs)) else
    order(as.character(runs$created), decreasing = TRUE, na.last = TRUE)
  selected_run <- as.character(selected_run %||% "")
  if (nzchar(selected_run) && selected_run %in% runs$id) {
    selected_index <- match(selected_run, runs$id)
    order_index <- c(selected_index, order_index[order_index != selected_index])
  }
  max_runs <- max(1L, as.integer(max_runs)[[1L]])
  omitted_runs <- max(0L, nrow(runs) - max_runs)
  runs <- runs[utils::head(order_index, max_runs), , drop = FALSE]
  diagnostic_names <- c(
    gof = "has_gof", vpc = "has_vpc", npc = "has_npc", npde = "has_npde",
    vpc_categorical = "has_vpc_categorical", vpc_count = "has_vpc_count",
    vpc_tte = "has_vpc_tte", vpc_competing = "has_vpc_competing",
    vpc_recurrent = "has_vpc_recurrent",
    bootstrap = "has_bootstrap", profile = "has_profile",
    scm = "has_scm", covariance = "has_covariance"
  )
  summaries <- unname(lapply(seq_len(nrow(runs)), function(index) {
    metadata <- as.list(runs[index, , drop = FALSE])
    parent_id <- as.character(metadata$parent_id)
    version_label <- unname(version_labels[parent_id])
    if (!length(version_label) || is.na(version_label) || !nzchar(version_label)) {
      version_label <- parent_id
    }
    saved_diagnostics <- tryCatch(
      nm_project_load_diagnostics(workspace, project, as.character(metadata$id)),
      error = function(error) list()
    )
    base <- list(
      id = as.character(metadata$id), label = as.character(metadata$label),
      model_version = as.character(version_label),
      run_number = suppressWarnings(as.integer(metadata$run_number)),
      result_type = as.character(metadata$result_type),
      method = as.character(metadata$method %||% ""),
      created = as.character(metadata$created),
      diagnostics = stats::setNames(unname(lapply(diagnostic_names, function(column) {
        name <- names(diagnostic_names)[[match(column, diagnostic_names)]]
        isTRUE(metadata[[column]]) || !is.null(saved_diagnostics[[name]])
      })), names(diagnostic_names))
    )
    if (identical(detail, "index")) {
      base$result_available <- TRUE
      return(base)
    }
    opened <- tryCatch(nm_project_load(workspace, project, metadata$id), error = identity)
    if (inherits(opened, "error")) {
      base$result_available <- FALSE
      base$message <- conditionMessage(opened)
      return(base)
    }
    result <- opened$result
    if (inherits(result, "nm_fit")) {
      estimates <- .liber_gui_parameter_values(result)
      keep <- utils::head(seq_along(estimates), max(1L, as.integer(max_parameters)[[1L]]))
      covariance <- result$covariance %||% NULL
      se <- covariance$se %||% numeric()
      rse <- covariance$rse %||% numeric()
      parameters <- unname(lapply(keep, function(parameter) {
        name <- names(estimates)[[parameter]]
        list(
          name = name, estimate = unname(estimates[[parameter]]),
          se = if (name %in% names(se)) unname(se[[name]]) else NULL,
          rse = if (name %in% names(rse)) unname(rse[[name]]) else NULL
        )
      }))
      timing <- result$timing %||% list()
      base$result_available <- TRUE
      base$method <- .nm_fit_method_label(result)
      base$method_sequence <- as.character(result$method_sequence %||% result$method)
      base$objective <- unname(result$objective)
      base$convergence <- unname(result$convergence)
      base$iterations <- suppressWarnings(as.integer(result$iterations %||%
        result$evaluations[["gradient"]] %||% result$evaluations[["function"]] %||% NA_integer_))
      base$parameters <- parameters
      base$parameters_omitted <- max(0L, length(estimates) - length(keep))
      base$covariance <- list(
        status = if (is.null(covariance)) "not requested" else
          as.character(covariance$status %||% "completed"),
        method = if (is.null(covariance)) "" else
          toupper(as.character(covariance$type %||% ""))
      )
      base$timing_seconds <- list(
        model_fit = unname(timing$model_fit_seconds %||% NULL),
        covariance = unname(timing$covariance_seconds %||% NULL),
        total = unname(timing$total_seconds %||% NULL)
      )
      base$output_columns <- names(result$output %||% data.frame())
      return(base)
    }
    if (is.data.frame(result)) {
      id_column <- if (".ID_INDEX" %in% names(result)) ".ID_INDEX" else
        if ("ID" %in% names(result)) "ID" else NULL
      base$result_available <- TRUE
      base$records <- nrow(result)
      base$subjects <- if (is.null(id_column)) NULL else length(unique(result[[id_column]]))
      base$output_columns <- names(result)
      return(base)
    }
    base$result_available <- FALSE
    base$message <- "The saved run has an unsupported result type."
    base
  }))
  list(
    available = TRUE, project = project, project_name = project_name,
    request_id = "", message = "", scope = detail,
    run_count = nrow(records[records$entry_type == "run" & records$has_result %in% TRUE, , drop = FALSE]),
    included_runs = length(summaries), omitted_runs = omitted_runs,
    runs = summaries
  )
}

.liber_gui_report_ai_context <- function(workspace, project, run_ids,
                                         max_parameters = 48L) {
  run_ids <- unique(as.character(run_ids %||% character()))
  if (!length(run_ids)) {
    return(list(
      available = FALSE, project = as.character(project %||% ""),
      project_name = "", request_id = "", run_ids = character(), runs = list(),
      message = "No model runs are selected in the report workflow."
    ))
  }
  context <- .liber_gui_ai_context(
    workspace, project, max_runs = length(run_ids),
    max_parameters = max_parameters, detail = "results", run_ids = run_ids
  )
  if (!isTRUE(context$available)) return(context)
  context$run_ids <- run_ids
  context$runs <- unname(lapply(context$runs, function(summary) {
    opened <- nm_project_load(workspace, project, summary$id)
    model <- opened$model
    result <- opened$result
    data <- opened$data
    summary$model <- list(
      name = attr(model, "name", exact = TRUE) %||% paste0("ADVAN", model$ADVAN, " model"),
      advan = model$ADVAN, trans = model$TRANS, solver = model$SOLVER,
      language = model$LANGUAGE, omega_structure = model$LIK_CONFIG$omega %||% "diagonal",
      theta_definitions = .liber_gui_rows(model$THETAS),
      omega_definitions = .liber_gui_rows(model$OMEGAS),
      sigma_definitions = .liber_gui_rows(model$SIGMAS),
      pred = as.character(model$PRED %||% ""),
      des = as.character(model$DES %||% ""),
      error = as.character(model$ERROR %||% "")
    )
    summary$data <- list(
      records = nrow(data), subjects = length(unique(data$.ID_INDEX)),
      observations = sum(data$EVID == 0L & data$MDV == 0L),
      columns = setdiff(names(data), grep("^[.]", names(data), value = TRUE))
    )
    if (inherits(result, "nm_fit")) {
      gof <- tryCatch(nm_gof(result), error = function(error) NULL)
      if (!is.null(gof)) {
        observed <- gof$EVID == 0L & gof$MDV == 0L & is.finite(gof$DV)
        gof <- gof[observed, , drop = FALSE]
        n_parameters <- sum(!result$model$THETAS$FIX) +
          sum(!result$model$SIGMAS$FIX) + sum(!result$model$OMEGAS$FIX)
        summary$gof_summary <- list(
          observations = nrow(gof), free_parameters = n_parameters,
          aic = result$objective + 2 * n_parameters,
          bic = result$objective + log(max(1, nrow(gof))) * n_parameters,
          population_rmse = sqrt(mean((gof$DV - gof$PRED)^2, na.rm = TRUE)),
          individual_rmse = sqrt(mean((gof$DV - gof$IPRED)^2, na.rm = TRUE)),
          mean_cwres = mean(gof$CWRES, na.rm = TRUE),
          sd_cwres = stats::sd(gof$CWRES, na.rm = TRUE)
        )
      }
    }
    diagnostics <- tryCatch(
      nm_project_load_diagnostics(workspace, project, summary$id),
      error = function(error) list()
    )
    summary$diagnostic_details <- unname(lapply(names(diagnostics), function(name) {
      item <- diagnostics[[name]]
      fields <- c("nsim", "seed", "pc_correct", "stratify", "level", "n", "status")
      details <- item[intersect(fields, names(item))]
      c(list(type = name, class = class(item)[[1L]] %||% "list"), details)
    }))
    summary
  }))
  context
}

.liber_gui_report <- function(report) {
  if (inherits(report, "nm_report_bundle")) {
    return(list(docx = report$docx, pdf = report$pdf, json = report$json,
                design_id = report$design$id %||% NULL))
  }
  if (!inherits(report, "nm_report")) return(NULL)
  list(docx = NULL, pdf = report$pdf, json = report$json)
}

.liber_gui_result <- function(result) {
  if (is.null(result)) return(list(status = "idle", message = "Ready"))
  if (inherits(result, "liber_gui_validation")) {
    return(utils::modifyList(
      list(status = "validated", message = "Model validation completed"),
      unclass(result)
    ))
  }
  if (inherits(result, "liber_gui_diagram_preview")) {
    return(list(
      status = "validated", kind = "diagram_preview",
      message = "Diagram code preview is ready",
      nonce = result$nonce, preview = result$preview,
      code_changed = isTRUE(result$code_changed)
    ))
  }
  if (inherits(result, "liber_gui_report_document")) {
    return(list(
      status = "completed", kind = "report_document",
      message = paste("Loaded report source", result$name),
      nonce = result$nonce, block_id = result$block_id,
      name = result$name, text = result$text
    ))
  }
  if (inherits(result, "liber_gui_control_export")) {
    return(list(status = "completed", kind = "control_export",
                message = paste("NONMEM control stream written to", result$path),
                path = result$path))
  }
  if (inherits(result, "liber_gui_queued")) {
    return(list(
      status = "queued", message = paste("Queued job", result$id), job_id = result$id
    ))
  }
  if (inherits(result, "error")) {
    return(list(status = "error", message = conditionMessage(result)))
  }
  if (inherits(result, "nm_npc")) {
    return(list(
      status = "completed", kind = "npc", message = "Numerical predictive check completed",
      nsim = result$nsim, outside_90 = result$outside_90,
      table = .liber_gui_rows(result$table, 5000L)
    ))
  }
  if (inherits(result, "nm_npde")) {
    return(list(
      status = "completed", kind = "npde", message = "NPDE completed",
      nsim = result$nsim, summary = as.list(result$summary),
      table = .liber_gui_rows(result$table, 5000L)
    ))
  }
  if (inherits(result, "nm_vpc_categorical")) {
    return(list(
      status = "completed", kind = "vpc_categorical",
      message = "Categorical VPC completed", nsim = result$nsim,
      outcome = result$outcome, categories = as.character(result$categories),
      observed = .liber_gui_rows(result$observed),
      simulated = .liber_gui_rows(result$simulated)
    ))
  }
  if (inherits(result, "nm_vpc_count")) {
    return(list(
      status = "completed", kind = "vpc_count", message = "Count VPC completed",
      nsim = result$nsim, outcome = result$outcome, family = result$family,
      dvid = result$dvid,
      observed = .liber_gui_rows(result$observed),
      simulated = .liber_gui_rows(result$simulated)
    ))
  }
  if (inherits(result, "nm_vpc_tte")) {
    return(list(
      status = "completed", kind = "vpc_tte", message = "TTE VPC completed",
      nsim = result$nsim, event = result$event,
      observed = .liber_gui_rows(result$observed),
      simulated = .liber_gui_rows(result$simulated)
    ))
  }
  if (inherits(result, "nm_vpc_competing")) {
    return(list(
      status = "completed", kind = "vpc_competing",
      message = "Competing-risk VPC completed", nsim = result$nsim,
      event = result$event, dvid = result$dvid,
      causes = as.character(result$causes),
      observed = .liber_gui_rows(result$observed),
      simulated = .liber_gui_rows(result$simulated)
    ))
  }
  if (inherits(result, "nm_vpc_recurrent")) {
    return(list(
      status = "completed", kind = "vpc_recurrent",
      message = "Recurrent-event VPC completed", nsim = result$nsim,
      event = result$event, dvid = result$dvid,
      observed = .liber_gui_rows(result$observed),
      simulated = .liber_gui_rows(result$simulated)
    ))
  }
  if (inherits(result, "nm_bootstrap")) {
    return(list(
      status = "completed", kind = "bootstrap", message = "Bootstrap completed",
      n = result$n, successful = result$successful, level = result$level,
      summary = .liber_gui_rows(result$summary), errors = as.character(result$errors)
    ))
  }
  if (inherits(result, "nm_profile")) {
    return(list(
      status = "completed", kind = "profile", message = "Profile likelihood completed",
      level = result$level, intervals = .liber_gui_rows(result$intervals),
      grid = .liber_gui_rows(result$grid, 5000L)
    ))
  }
  if (inherits(result, "nm_scm")) {
    return(list(
      status = "completed", kind = "scm", message = "SCM completed",
      base_objective = result$base_objective, final_objective = result$final_objective,
      selected = .liber_gui_rows(result$selected), steps = .liber_gui_rows(result$steps)
    ))
  }
  if (inherits(result, "liber_gui_comparison")) {
    return(list(
      status = "completed", kind = "comparison", message = "Estimation runs compared",
      comparison_id = result$id %||% "comparison",
      parameters = .liber_gui_rows(result$parameters),
      gof = .liber_gui_rows(result$gof),
      runs = .liber_gui_rows(result$runs),
      plots = result$plots %||% list()
    ))
  }
  if (inherits(result, "nm_report") || inherits(result, "nm_report_bundle")) {
    return(list(status = "completed", kind = "report", message = "Report generated"))
  }
  if (inherits(result, "nm_fit")) {
    etab <- nm_etab(result)
    parameters <- .liber_gui_parameter_values(result)
    return(list(
      status = if (result$convergence == 0L) "completed" else "warning",
      kind = "estimation", message = paste(result$method, "estimation completed"),
      method = result$method, objective = result$objective,
      convergence = result$convergence,
      parameters = unname(lapply(seq_along(parameters), function(i) {
        list(name = names(parameters)[[i]], value = unname(parameters[[i]]))
      })),
      shrinkage = unname(lapply(seq_along(etab$shrinkage), function(i) {
        list(name = names(etab$shrinkage)[[i]], value = etab$shrinkage[[i]])
      }))
    ))
  }
  if (inherits(result, "nm_vpc")) {
    return(list(
      status = "completed", kind = "vpc", message = "Visual predictive check completed",
      nsim = result$nsim, observed = .liber_gui_rows(result$observed),
      simulated = .liber_gui_rows(result$simulated),
      points = .liber_gui_rows(result$points %||% data.frame(), 10000L),
      pc_correct = isTRUE(result$pc_correct),
      stratify = result$stratify %||% NULL,
      stratified = unname(lapply(result$stratified %||% list(), function(item) {
        list(
          level = as.character(item$level),
          observed = .liber_gui_rows(item$observed),
          simulated = .liber_gui_rows(item$simulated),
          points = .liber_gui_rows(item$points %||% data.frame(), 10000L),
          nsim = result$nsim,
          pc_correct = isTRUE(result$pc_correct)
        )
      }))
    ))
  }
  list(
    status = "completed", kind = "simulation", message = "Simulation completed",
    rows = nrow(result), finite_predictions = sum(is.finite(result$IPRED)),
    solver = attr(result, "solver") %||% "C++"
  )
}

.liber_gui_diagnostics <- function(diagnostics = list(), payload = names(diagnostics)) {
  diagnostics <- diagnostics %||% list()
  types <- c("vpc", "npc", "npde", "vpc_categorical", "vpc_count", "vpc_tte",
             "vpc_competing",
             "vpc_recurrent",
             "bootstrap", "profile", "scm")
  payload <- intersect(as.character(payload %||% character()), types)
  list(
    available = stats::setNames(lapply(types, function(name) !is.null(diagnostics[[name]])), types),
    vpc = if (is.null(diagnostics$vpc) || !"vpc" %in% payload) NULL else .liber_gui_result(diagnostics$vpc),
    npc = if (is.null(diagnostics$npc) || !"npc" %in% payload) NULL else .liber_gui_result(diagnostics$npc),
    npde = if (is.null(diagnostics$npde) || !"npde" %in% payload) NULL else .liber_gui_result(diagnostics$npde),
    vpc_categorical = if (is.null(diagnostics$vpc_categorical) || !"vpc_categorical" %in% payload) NULL else .liber_gui_result(diagnostics$vpc_categorical),
    vpc_count = if (is.null(diagnostics$vpc_count) || !"vpc_count" %in% payload) NULL else .liber_gui_result(diagnostics$vpc_count),
    vpc_tte = if (is.null(diagnostics$vpc_tte) || !"vpc_tte" %in% payload) NULL else .liber_gui_result(diagnostics$vpc_tte),
    vpc_competing = if (is.null(diagnostics$vpc_competing) || !"vpc_competing" %in% payload) NULL else .liber_gui_result(diagnostics$vpc_competing),
    vpc_recurrent = if (is.null(diagnostics$vpc_recurrent) || !"vpc_recurrent" %in% payload) NULL else .liber_gui_result(diagnostics$vpc_recurrent),
    bootstrap = if (is.null(diagnostics$bootstrap) || !"bootstrap" %in% payload) NULL else .liber_gui_result(diagnostics$bootstrap),
    profile = if (is.null(diagnostics$profile) || !"profile" %in% payload) NULL else .liber_gui_result(diagnostics$profile),
    scm = if (is.null(diagnostics$scm) || !"scm" %in% payload) NULL else .liber_gui_result(diagnostics$scm)
  )
}

.liber_gui_hmm <- function(decoded = NULL, available = FALSE, limit = 50000L) {
  unloaded <- function() list(
    available = isTRUE(available), loaded = FALSE, states = list(), rows = list(),
    sequence_summary = list(), observations = 0L, sequences = 0L,
    truncated = FALSE, log_likelihood = NULL, eta_type = ""
  )
  if (is.null(decoded)) return(unloaded())
  if (!inherits(decoded, "nm_hmm_decode")) {
    .nm_stop("`decoded` must be returned by `nm_hmm_decode()`.")
  }
  method <- as.character(attr(decoded, "method", exact = TRUE) %||% "")
  if (!identical(method, "all")) {
    .nm_stop("The GUI HMM payload requires `nm_hmm_decode(..., method = \"all\")`.")
  }
  state_names <- as.character(attr(decoded, "states", exact = TRUE) %||% character())
  state_keys <- make.names(state_names, unique = TRUE)
  frame <- as.data.frame(decoded, stringsAsFactors = FALSE)
  observed <- if ("HMM_ROW_NLL" %in% names(frame)) {
    is.finite(suppressWarnings(as.numeric(frame$HMM_ROW_NLL)))
  } else rep(TRUE, nrow(frame))
  frame <- frame[observed, , drop = FALSE]
  frame$SUBJECT <- if ("ID" %in% names(frame)) frame$ID else frame$.ID_INDEX
  frame$SEQUENCE <- if ("DVID" %in% names(frame)) frame$DVID else 1L
  keep <- unique(c(
    "SUBJECT", "SEQUENCE", "TIME", "DVID", "HMM_ROW_NLL",
    grep("^HMM_(FILTER|SMOOTH|VITERBI)_", names(frame), value = TRUE)
  ))
  keep <- intersect(keep, names(frame))
  total_observations <- nrow(frame)
  summary <- as.data.frame(
    attr(decoded, "sequence_summary", exact = TRUE) %||% data.frame(),
    stringsAsFactors = FALSE
  )
  if (nrow(summary)) {
    summary$SUBJECT <- if ("ID" %in% names(summary)) {
      summary$ID
    } else summary$.ID_INDEX
    summary$SEQUENCE <- if ("DVID" %in% names(summary)) {
      summary$DVID
    } else if ("HMM_SEQUENCE" %in% names(summary)) {
      summary$HMM_SEQUENCE
    } else 1L
    summary_keep <- intersect(unique(c(
      "SUBJECT", "SEQUENCE", "LOG_LIKELIHOOD", "VITERBI_LOG_JOINT",
      "VITERBI_LOG_POSTERIOR"
    )), names(summary))
    summary <- summary[, summary_keep, drop = FALSE]
  }
  list(
    available = TRUE, loaded = TRUE,
    states = unname(lapply(seq_along(state_names), function(index) list(
      label = state_names[[index]], key = state_keys[[index]], index = index
    ))),
    rows = .liber_gui_rows(frame[, keep, drop = FALSE], limit),
    sequence_summary = .liber_gui_rows(summary),
    observations = total_observations,
    sequences = nrow(summary),
    truncated = total_observations > as.integer(limit),
    log_likelihood = as.numeric(attr(decoded, "log_likelihood", exact = TRUE)),
    eta_type = as.character(attr(decoded, "eta_type", exact = TRUE) %||% "")
  )
}

.liber_gui_kalman <- function(decoded = NULL, available = FALSE, limit = 50000L) {
  unloaded <- function() list(
    available = isTRUE(available), loaded = FALSE, states = list(), rows = list(),
    observations = 0L, sequences = 0L, truncated = FALSE,
    log_likelihood = NULL, eta_type = "", filter = "", smoother = ""
  )
  if (is.null(decoded)) return(unloaded())
  if (!inherits(decoded, "nm_kalman_decode")) {
    .nm_stop("`decoded` must be returned by `nm_kalman_decode()`.")
  }
  state_names <- as.character(attr(decoded, "states", exact = TRUE) %||% character())
  state_keys <- make.names(state_names, unique = TRUE)
  frame <- as.data.frame(decoded, stringsAsFactors = FALSE)
  observed <- if ("KF_ROW_NLL" %in% names(frame)) {
    is.finite(suppressWarnings(as.numeric(frame$KF_ROW_NLL)))
  } else rep(TRUE, nrow(frame))
  frame <- frame[observed, , drop = FALSE]
  frame$SUBJECT <- if ("ID" %in% names(frame)) frame$ID else frame$.ID_INDEX
  frame$SEQUENCE <- if ("DVID" %in% names(frame)) frame$DVID else 1L
  if ("KF_INNOVATION_VARIANCE" %in% names(frame)) {
    denominator <- sqrt(suppressWarnings(as.numeric(frame$KF_INNOVATION_VARIANCE)))
    frame$KF_STANDARDIZED_INNOVATION <- suppressWarnings(
      as.numeric(frame$KF_INNOVATION) / denominator
    )
    frame$KF_STANDARDIZED_INNOVATION[!is.finite(frame$KF_STANDARDIZED_INNOVATION)] <- NA_real_
  }
  keep <- unique(c(
    "SUBJECT", "SEQUENCE", "TIME", "DVID", "DV", "KF_INNOVATION",
    "KF_INNOVATION_VARIANCE", "KF_STANDARDIZED_INNOVATION", "KF_ROW_NLL",
    grep("^KF_(PRED|FILTER|SMOOTH|FILTER_SD|SMOOTH_SD)_", names(frame), value = TRUE)
  ))
  keep <- intersect(keep, names(frame))
  total_observations <- nrow(frame)
  sequence_keys <- if (nrow(frame)) {
    unique(paste(frame$SUBJECT, frame$SEQUENCE, sep = "\r"))
  } else character()
  list(
    available = TRUE, loaded = TRUE,
    states = unname(lapply(seq_along(state_names), function(index) list(
      label = state_names[[index]], key = state_keys[[index]], index = index
    ))),
    rows = .liber_gui_rows(frame[, keep, drop = FALSE], limit),
    observations = total_observations,
    sequences = length(sequence_keys),
    truncated = total_observations > as.integer(limit),
    log_likelihood = as.numeric(attr(decoded, "log_likelihood", exact = TRUE)),
    eta_type = as.character(attr(decoded, "eta_type", exact = TRUE) %||% ""),
    filter = as.character(attr(decoded, "filter", exact = TRUE) %||% "linear"),
    smoother = as.character(attr(decoded, "smoother", exact = TRUE) %||% "RTS")
  )
}

.liber_gui_report_design <- function(design = NULL) {
  if (!inherits(design, "nm_report_design")) return(NULL)
  list(
    id = design$id,
    title = design$title,
    name = design$style$filename %||% "liberation-report",
    directory = design$style$output_directory %||% "",
    formats = unname(as.character(design$formats)),
    updated = design$updated,
    blocks = unname(lapply(design$blocks, function(block) list(
      id = block$id, type = block$type, title = block$title,
      source = block$source, text = block$text,
      run_ids = unname(as.character(block$run_ids)),
      elements = unname(as.character(block$elements)),
      options = block$options %||% list()
    )))
  )
}

#' React workbench for LibeRation
#'
#' @param model Optional [nm_model()] or [NMEngine].
#' @param data Optional NONMEM-style dataset.
#' @param jobs Optional job table, such as `LibeRQueue$list()`.
#' @param result Optional latest simulation/result summary.
#' @param fit Optional fitted model used for parameter and diagnostic panels.
#' @param report Optional generated [nm_report()].
#' @param report_design Optional saved [nm_report_design()] restored by the
#'   visual report designer.
#' @param diagnostics Optional saved VPC, NPC, and NPDE results for the selected
#'   estimation run.
#' @param hmm Optional decoded HMM GUI payload. When omitted, availability is
#'   inferred from the selected model and fit and rows remain lazy.
#' @param kalman Optional decoded linear state-space GUI payload. When omitted,
#'   availability is inferred from the selected model and fit and rows remain
#'   lazy.
#' @param log Optional application-log payload.
#' @param job_log Optional stdout/stderr lines for the selected queued job.
#' @param server Optional runtime/server status passed to the workbench.
#' @param workspace Optional [nm_workspace()] and current project metadata.
#' @param library Optional LibeRary catalogue payload.
#' @param ai Optional browser-local WebGPU AI settings, independent Help and
#'   Report context-window choices, and model catalogue.
#' @param ai_context Optional compact, on-demand summaries of saved project
#'   runs supplied to the browser-local Help model. Full result datasets are
#'   never included in this payload.
#' @param report_ai_context Optional detailed, on-demand summaries of the runs
#'   selected in the visual report workflow.
#' @param report_directory Default report output directory displayed by the
#'   visual report builder.
#' @param output_catalog Optional draft output catalogue used by the GUI after
#'   validation but before editor changes are applied.
#' @param run_output Optional selected model-run output columns, aligned by a
#'   `.ROW` column or by row position. These are loaded into Data explorer only
#'   when its row payload is requested.
#' @param data_payload Whether dataset rows are included in the browser payload.
#'   When false, only dataset metadata and column names are included.
#' @param gof_payload Whether GOF rows are generated when `fit` is an [nm_fit].
#' @param diagnostic_payload Names of saved diagnostics whose plot data should
#'   be included. Availability metadata is always included.
#' @param input_id Shiny event input prefix.
#' @param width,height Widget dimensions.
#' @param elementId Optional HTML widget element id.
#' @export
liber_workbench <- function(model = NULL, data = NULL, jobs = NULL, result = NULL,
                             fit = NULL, report = NULL, report_design = NULL,
                             diagnostics = NULL, hmm = NULL, kalman = NULL,
                             log = NULL, job_log = NULL,
                             server = NULL,
                             workspace = NULL,
                             library = NULL, ai = NULL, ai_context = NULL,
                             report_ai_context = NULL, report_directory = "",
                             output_catalog = NULL,
                             run_output = if (inherits(fit, "nm_fit")) fit$output else NULL,
                             data_payload = TRUE, gof_payload = TRUE,
                             diagnostic_payload = names(diagnostics),
                             input_id = "liber_workbench", width = NULL,
                            height = "780px", elementId = NULL) {
  jobs <- as.data.frame(jobs %||% data.frame(), stringsAsFactors = FALSE)
  model_payload <- .liber_gui_model(model, output_catalog = output_catalog)
  fit_payload <- if (is.list(fit) && !inherits(fit, "nm_fit") &&
                     !is.null(fit$available)) {
    fit
  } else {
    .liber_gui_fit(fit, include_gof = gof_payload)
  }
  hmm_payload <- if (is.list(hmm) && !inherits(hmm, "nm_hmm_decode") &&
                     !is.null(hmm$available)) {
    hmm
  } else {
    .liber_gui_hmm(
      hmm,
      available = !is.null(model_payload$hmm) && isTRUE(fit_payload$available)
    )
  }
  kalman_payload <- if (is.list(kalman) && !inherits(kalman, "nm_kalman_decode") &&
                        !is.null(kalman$available)) {
    kalman
  } else {
    .liber_gui_kalman(
      kalman,
      available = !is.null(model_payload$kalman) && isTRUE(fit_payload$available)
    )
  }
  content <- reactR::component("LibeRWorkbench", list(
    model = model_payload,
    dataset = .liber_gui_data(
      data, include_rows = data_payload,
      run_output = if (isTRUE(data_payload)) run_output else NULL
    ),
    jobs = unname(lapply(seq_len(nrow(jobs)), function(i) as.list(jobs[i, , drop = FALSE]))),
    result = .liber_gui_result(result),
    diagnostics = .liber_gui_diagnostics(diagnostics, diagnostic_payload),
    fit = fit_payload,
    hmm = hmm_payload,
    kalman = kalman_payload,
    report = .liber_gui_report(report),
    report_design = .liber_gui_report_design(report_design),
    log = log %||% list(level = "info", current = "Workbench ready",
                       history = "Workbench ready"),
    job_log = as.character(job_log %||% character()),
    workspace = workspace %||% .liber_gui_workspace(),
    library = library %||% .liber_gui_library(),
    ai_context = ai_context %||% list(
      available = FALSE, project = "", project_name = "", request_id = "",
      scope = "index", message = "Project result summaries have not been requested.",
      runs = list()
    ),
    report_ai_context = report_ai_context %||% list(
      available = FALSE, project = "", project_name = "", request_id = "",
      run_ids = character(), message = "Report evidence has not been requested.",
      runs = list()
    ),
    report_directory = as.character(report_directory %||% "")[[1L]],
    ai = ai %||% list(
      activated = FALSE, consented = FALSE,
      help_model = .liber_ai_default_help_model(),
      report_model = .liber_ai_default_report_model(),
      help_context = "auto", report_context = "auto",
      model = .liber_ai_default_help_model(),
      worker_url = "", models = list(), secure_context = TRUE
    ),
    server = server %||% list(
      mode = "local", connected = TRUE, platform = R.version$platform,
      worker = "in-process", isolation = "current R session"
    ),
    support = unname(lapply(seq_len(nrow(nm_support_matrix())), function(i) {
      as.list(nm_support_matrix()[i, , drop = FALSE])
    })),
    inputId = input_id
  ))
  htmlwidgets::createWidget(
    name = "liberWorkbench", reactR::reactMarkup(content),
    width = width, height = height, package = "LibeRation", elementId = elementId
  )
}

#' @noRd
widget_html.liberWorkbench <- function(id, style, class, ...) {
  htmltools::attachDependencies(
    htmltools::tags$div(id = id, class = class, style = style),
    list(
      reactR::html_dependency_corejs(),
      reactR::html_dependency_react(),
      reactR::html_dependency_reacttools()
    )
  )
}

#' Shiny output for the LibeR workbench
#' @param outputId Shiny output id.
#' @param width,height CSS dimensions.
#' @export
liberWorkbenchOutput <- function(outputId, width = "100%", height = "100vh") {
  htmlwidgets::shinyWidgetOutput(
    outputId, "liberWorkbench", width, height, package = "LibeRation"
  )
}

.liber_full_page_ui <- function(head, output) {
  head <- htmltools::tagAppendChild(
    head,
    htmltools::tags$style(htmltools::HTML(paste(
      "html, body { width: 100%; height: 100%; margin: 0; overflow: hidden; }",
      "body > .container-fluid { width: 100%; height: 100vh;",
      "min-width: 0; margin: 0; padding: 0; overflow: hidden; }"
    )))
  )
  shiny::fluidPage(
    head,
    htmltools::tags$div(
      class = "liberation-app-root",
      style = paste(
        "width: 100%; height: 100vh; min-width: 0;",
        "margin: 0; padding: 0; overflow: hidden;"
      ),
      output
    )
  )
}

#' Render a LibeR workbench in Shiny
#' @param expr Expression returning [liber_workbench()].
#' @param env Evaluation environment.
#' @param quoted Whether `expr` is quoted.
#' @export
renderLiberWorkbench <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) expr <- substitute(expr)
  htmlwidgets::shinyRenderWidget(expr, liberWorkbenchOutput, env, quoted = TRUE)
}

.liber_parameter_table_update <- function(table, values, label) {
  if (is.null(values)) return(table)
  if (is.data.frame(values)) {
    values <- lapply(seq_len(nrow(values)), function(i) as.list(values[i, , drop = FALSE]))
  }
  number_or <- function(value, default) {
    output <- suppressWarnings(as.numeric(value))
    if (length(output) != 1L || is.na(output) || !is.finite(output)) default else output
  }
  if (identical(label, "OMEGA")) {
    if (!length(values)) return(table[0L, , drop = FALSE])
    rows <- lapply(seq_along(values), function(i) {
      value <- values[[i]]
      data.frame(
        OMEGA = i,
        Value = number_or(value$Value, if (i <= nrow(table)) table$Value[[i]] else 0.1),
        FIX = isTRUE(value$FIX),
        ROW = as.integer(value$ROW %||% value$OMEGA %||% i),
        COL = as.integer(value$COL %||% value$OMEGA %||% i),
        stringsAsFactors = FALSE
      )
    })
    return(do.call(rbind, rows))
  }
  if (!length(values)) {
    return(if (identical(label, "THETA")) {
      data.frame(THETA = integer(), Value = numeric(), FIX = logical(),
                 LOWER = numeric(), UPPER = numeric())
    } else {
      data.frame(SIGMA = integer(), Value = numeric(), FIX = logical())
    })
  }
  rows <- lapply(seq_along(values), function(i) {
    value <- values[[i]]
    initial <- if (i <= nrow(table)) table$Value[[i]] else
      if (identical(label, "THETA")) 1 else 0.1
    output <- data.frame(
      INDEX = i, Value = number_or(value$Value, initial),
      FIX = isTRUE(value$FIX), stringsAsFactors = FALSE
    )
    names(output)[[1L]] <- label
    if (identical(label, "THETA")) {
      output$LOWER <- number_or(value$LOWER, NA_real_)
      output$UPPER <- number_or(value$UPPER, NA_real_)
    }
    output
  })
  do.call(rbind, rows)
}

.liber_code_reference_max <- function(code, functions) {
  code <- paste(as.character(code %||% ""), collapse = "\n")
  code <- gsub("(?s)/\\*.*?\\*/", " ", code, perl = TRUE)
  code <- gsub("(?m)//.*$|#.*$", " ", code, perl = TRUE)
  pattern <- paste0(
    "(?i)\\b(?:", paste(functions, collapse = "|"),
    ")\\s*\\(\\s*([1-9][0-9]*)\\s*\\)"
  )
  matches <- regmatches(code, gregexpr(pattern, code, perl = TRUE))[[1L]]
  if (!length(matches)) return(0L)
  indices <- suppressWarnings(as.integer(sub(pattern, "\\1", matches, perl = TRUE)))
  if (!length(indices) || all(is.na(indices))) 0L else max(indices, na.rm = TRUE)
}

.liber_model_parameter_requirements <- function(arguments) {
  code <- unname(unlist(arguments[c("PRED", "DES", "ALG", "ERROR")], use.names = FALSE))
  list(
    theta = .liber_code_reference_max(code, "THETA"),
    eta = .liber_code_reference_max(code, "ETA"),
    sigma = .liber_code_reference_max(code, c("ERR", "EPS", "SIGMA"))
  )
}

.liber_parameter_table_ensure <- function(table, count, label) {
  count <- max(nrow(table), as.integer(count %||% 0L))
  if (count <= 0L) return(.liber_parameter_table_update(table, list(), label))
  values <- rep(if (identical(label, "THETA")) 1 else 0.1, count)
  fixed <- rep(FALSE, count)
  if (nrow(table)) {
    values[seq_len(nrow(table))] <- table$Value
    fixed[seq_len(nrow(table))] <- table$FIX
  }
  output <- data.frame(INDEX = seq_len(count), Value = values, FIX = fixed,
                       stringsAsFactors = FALSE)
  names(output)[[1L]] <- label
  if (identical(label, "THETA")) {
    output$LOWER <- output$UPPER <- NA_real_
    if (nrow(table)) {
      if ("LOWER" %in% names(table)) output$LOWER[seq_len(nrow(table))] <- table$LOWER
      if ("UPPER" %in% names(table)) output$UPPER[seq_len(nrow(table))] <- table$UPPER
    }
  }
  output
}

.liber_omega_table_ensure <- function(table, count, structure = "diagonal") {
  count <- max(.nm_n_eta(table), as.integer(count %||% 0L))
  if (count <= 0L) {
    return(data.frame(OMEGA = integer(), Value = numeric(), FIX = logical(),
                      ROW = integer(), COL = integer()))
  }
  positions <- if (identical(structure, "full")) {
    do.call(rbind, lapply(seq_len(count), function(row) {
      data.frame(ROW = row, COL = seq_len(row))
    }))
  } else data.frame(ROW = seq_len(count), COL = seq_len(count))
  output <- data.frame(
    OMEGA = seq_len(nrow(positions)),
    Value = ifelse(positions$ROW == positions$COL, 0.1, 0),
    FIX = FALSE, ROW = positions$ROW, COL = positions$COL,
    stringsAsFactors = FALSE
  )
  if (nrow(table)) {
    existing <- match(
      paste(output$ROW, output$COL, sep = ":"),
      paste(table$ROW, table$COL, sep = ":")
    )
    keep <- which(!is.na(existing))
    output$Value[keep] <- table$Value[existing[keep]]
    output$FIX[keep] <- table$FIX[existing[keep]]
  }
  output
}

.liber_prior_table_update <- function(values) {
  if (is.null(values) || !length(values)) return(NULL)
  if (is.data.frame(values)) {
    values <- lapply(seq_len(nrow(values)), function(i) as.list(values[i, , drop = FALSE]))
  }
  number_or <- function(value, default) {
    output <- suppressWarnings(as.numeric(value))
    if (length(output) != 1L || is.na(output)) default else output
  }
  rows <- lapply(values, function(value) {
    nm_prior(
      parameter = as.character(value$parameter %||% ""),
      distribution = as.character(value$distribution %||% "normal"),
      mean = number_or(value$mean, 0), sd = number_or(value$sd, 1),
      shape = number_or(value$shape, NA_real_),
      rate = number_or(value$rate, NA_real_)
    )
  })
  do.call(rbind, rows)
}

.liber_model_from_event <- function(model, event) {
  if (inherits(model, "NMEngine")) model <- model$model
  if (!inherits(model, "nm_model")) .nm_stop("Load a model before editing it.")
  arguments <- model[intersect(names(model), names(formals(nm_model)))]
  old_advan <- model$ADVAN
  arguments$ADVAN <- as.integer(event$advan %||% old_advan)
  arguments$TRANS <- as.integer(event$trans %||% model$TRANS)
  arguments$PRED <- as.character(event$pred %||% model$PRED)
  arguments$ERROR <- as.character(event$error %||% model$ERROR)
  arguments$DES <- as.character(event$des %||% model$DES)
  arguments$ALG <- as.character(event$alg %||% model$ALG %||% "")
  arguments$INPUT <- unique(as.character(unlist(event$input %||% model$INPUT)))
  arguments$OUTPUT <- unique(as.character(unlist(event$output %||% model$OUTPUT %||% character())))
  likelihood <- as.list(model$LIK_CONFIG)
  likelihood$version <- NULL
  inferred_error <- .nm_error_type(arguments$ERROR, "auto")
  if (identical(inferred_error, "likelihood")) {
    likelihood$error <- "likelihood"
  } else if (identical(likelihood$error, "likelihood") &&
             !identical(inferred_error, "none")) {
    likelihood$error <- inferred_error
  }
  omega_structure <- as.character(
    event$omega_structure %||% model$LIK_CONFIG$omega %||% "diagonal"
  )[[1L]]
  requirements <- .liber_model_parameter_requirements(arguments)
  arguments$THETAS <- .liber_parameter_table_ensure(
    .liber_parameter_table_update(model$THETAS, event[["theta", exact = TRUE]], "THETA"),
    requirements$theta, "THETA"
  )
  arguments$OMEGAS <- .liber_omega_table_ensure(
    .liber_parameter_table_update(model$OMEGAS, event[["omega", exact = TRUE]], "OMEGA"),
    requirements$eta, omega_structure
  )
  arguments$SIGMAS <- .liber_parameter_table_ensure(
    .liber_parameter_table_update(model$SIGMAS, event[["sigma", exact = TRUE]], "SIGMA"),
    requirements$sigma, "SIGMA"
  )
  likelihood$omega <- omega_structure
  likelihood$priors <- if (is.null(event$priors)) {
    model$LIK_CONFIG$priors
  } else {
    .liber_prior_table_update(event$priors)
  }
  arguments$LIK_CONFIG <- .nm_lik_config(
    likelihood, likelihood$error, as.integer(likelihood$iov %||% model$IOV)
  )
  if (!identical(arguments$ADVAN, old_advan)) arguments$GRAPH <- NULL
  arguments$ERROR_TYPE <- "auto"
  edited <- do.call(nm_model, arguments)
  requested_states <- as.integer(event$n_state %||% edited$n_state)
  if (edited$ADVAN %in% c(6L, 13L) && requested_states != edited$n_state) {
    .nm_stop(
      "The compartment count is derived from $DES. Define DADT(1) through DADT(",
      requested_states, ") to use that count."
    )
  }
  attr(edited, "name") <- trimws(as.character(
    event$problem %||% attr(model, "name", exact = TRUE) %||% "Untitled model"
  ))
  edited
}

.liber_estimation_arguments <- function(event) {
  method <- toupper(as.character(event$method %||% "FOCEI"))
  arguments <- list(
    method = method,
    maxit = as.integer(event$maxit %||% 200L),
    eta_maxit = as.integer(event$etaMaxit %||% 100L),
    tolerance = as.numeric(event$tolerance %||% 1e-6),
    n_cores = max(1L, as.integer(event$nCores %||% 1L)),
    print_every = max(0L, as.integer(event$printEvery %||% 0L)),
    covariance = isTRUE(event$covariance),
    covariance_type = as.character(event$covarianceType %||% "hessian"),
    covariance_tolerance = as.numeric(event$covarianceTolerance %||% 1e-8),
    covariance_samples = as.integer(event$covarianceSamples %||% 200L),
    covariance_seed = as.integer(event$covarianceSeed %||%
      event$methodSeed %||% 20260713L)
  )
  if (identical(method, "IMP")) {
    arguments$n_imp <- as.integer(event$nImp %||% 200L)
    arguments$seed <- as.integer(event$methodSeed %||% 20260713L)
  } else if (identical(method, "GQ")) {
    arguments$gq_order <- as.integer(event$gqOrder %||% 5L)
    arguments$gq_grid <- as.character(event$gqGrid %||% "auto")
    arguments$gq_level <- as.integer(event$gqLevel %||% 3L)
    arguments$gq_adaptive <- isTRUE(event$gqAdaptive %||% TRUE)
    arguments$gq_max_points <- as.integer(event$gqMaxPoints %||% 100000L)
  } else if (identical(method, "SAEM")) {
    arguments$n_iter <- as.integer(event$nIter %||% event$maxit %||% 200L)
    arguments$burn <- as.integer(event$burn %||% floor(arguments$n_iter / 3))
    arguments$mcmc_steps <- as.integer(event$mcmcSteps %||% 2L)
    arguments$seed <- as.integer(event$methodSeed %||% 20260713L)
  } else if (identical(method, "BAYES")) {
    arguments$n_burn <- as.integer(event$nBurn %||% 500L)
    arguments$n_sample <- as.integer(event$nSample %||% 1000L)
    arguments$n_thin <- as.integer(event$nThin %||% 1L)
    arguments$seed <- as.integer(event$methodSeed %||% 20260713L)
  } else if (method %in% c("HMC", "NUTS")) {
    arguments$n_warmup <- as.integer(event$nBurn %||% 500L)
    arguments$n_sample <- as.integer(event$nSample %||% 1000L)
    arguments$n_thin <- as.integer(event$nThin %||% 1L)
    arguments$n_chains <- as.integer(event$nChains %||% 4L)
    arguments$target_acceptance <- as.numeric(event$targetAcceptance %||% 0.8)
    arguments$max_depth <- as.integer(event$maxTreeDepth %||% 10L)
    arguments$n_leapfrog <- as.integer(event$nLeapfrog %||% 10L)
    arguments$seed <- as.integer(event$methodSeed %||% 20260719L)
  } else if (method %in% c("NPML", "NPAG")) {
    arguments$np_points <- as.integer(event$npPoints %||% 25L)
    arguments$np_cycles <- as.integer(event$npCycles %||% 3L)
    arguments$np_max_support <- as.integer(event$npMaxSupport %||% 100L)
    arguments$np_grid_step <- as.numeric(event$npGridStep %||% 1)
    arguments$seed <- as.integer(event$methodSeed %||% 20260719L)
  }
  arguments
}

.liber_estimation_stages <- function(event) {
  raw <- event$stages %||% list(event)
  if (is.data.frame(raw)) raw <- lapply(seq_len(nrow(raw)), function(i) as.list(raw[i, , drop = FALSE]))
  if (!is.list(raw) || !length(raw)) raw <- list(event)
  lapply(raw, function(stage) {
    controls <- .liber_estimation_arguments(stage)
    method <- controls$method
    controls$method <- NULL
    do.call(nm_est_stage, c(list(method = method), controls))
  })
}

.liber_model_template <- function(advan = 1L, trans = NULL, n_state = NULL,
                                  problem = NULL) {
  advan <- as.integer(advan)
  input <- c("ID", "TIME", "EVID", "AMT", "RATE", "II", "SS", "CMT", "DV", "MDV")
  theta <- function(values) data.frame(
    THETA = seq_along(values), Value = values, FIX = FALSE
  )
  omega <- function(n) data.frame(
    OMEGA = seq_len(n), Value = rep(0.1, n), FIX = FALSE
  )
  sigma <- data.frame(SIGMA = 1L, Value = 0.05, FIX = FALSE)
  specification <- switch(
    as.character(advan),
    `1` = list(
      trans = 2L, values = c(5, 50),
      pred = "CL = THETA(1) * exp(ETA(1))\nV = THETA(2) * exp(ETA(2))\nS1 = V",
      des = "", eta = 2L
    ),
    `2` = list(
      trans = 2L, values = c(1, 5, 50),
      pred = "KA = THETA(1) * exp(ETA(1))\nCL = THETA(2) * exp(ETA(2))\nV = THETA(3) * exp(ETA(3))\nS2 = V",
      des = "", eta = 3L
    ),
    `3` = list(
      trans = 4L, values = c(5, 30, 8, 70),
      pred = "CL = THETA(1) * exp(ETA(1))\nV1 = THETA(2) * exp(ETA(2))\nQ = THETA(3)\nV2 = THETA(4)\nS1 = V1",
      des = "", eta = 2L
    ),
    `4` = list(
      trans = 4L, values = c(1, 5, 30, 8, 70),
      pred = "KA = THETA(1) * exp(ETA(1))\nCL = THETA(2) * exp(ETA(2))\nV1 = THETA(3) * exp(ETA(3))\nQ = THETA(4)\nV2 = THETA(5)\nS2 = V1",
      des = "", eta = 3L
    ),
    `6` = list(
      trans = 1L, values = c(5, 50),
      pred = "CL = THETA(1) * exp(ETA(1))\nV = THETA(2) * exp(ETA(2))\nK = CL / V\nS1 = V",
      des = "DADT(1) = -K * A(1)", eta = 2L
    ),
    `11` = list(
      trans = 4L, values = c(5, 20, 8, 40, 4, 80),
      pred = "CL = THETA(1) * exp(ETA(1))\nV1 = THETA(2) * exp(ETA(2))\nQ2 = THETA(3)\nV2 = THETA(4)\nQ3 = THETA(5)\nV3 = THETA(6)\nS1 = V1",
      des = "", eta = 2L
    ),
    `12` = list(
      trans = 4L, values = c(1, 5, 20, 8, 40, 4, 80),
      pred = "KA = THETA(1) * exp(ETA(1))\nCL = THETA(2) * exp(ETA(2))\nV1 = THETA(3) * exp(ETA(3))\nQ2 = THETA(4)\nV2 = THETA(5)\nQ3 = THETA(6)\nV3 = THETA(7)\nS2 = V1",
      des = "", eta = 3L
    ),
    `13` = list(
      trans = 1L, values = c(5, 30, 8, 70),
      pred = paste(
        "CL = THETA(1) * exp(ETA(1))", "V1 = THETA(2) * exp(ETA(2))",
        "Q = THETA(3)", "V2 = THETA(4)", "K10 = CL / V1",
        "K12 = Q / V1", "K21 = Q / V2", "S1 = V1", sep = "\n"
      ),
      des = "DADT(1) = -(K10 + K12) * A(1) + K21 * A(2)\nDADT(2) = K12 * A(1) - K21 * A(2)",
      eta = 2L
    ),
    .nm_stop("No GUI template is available for ADVAN", advan, ".")
  )
  if (!is.null(trans) && !advan %in% c(6L, 13L)) {
    trans <- as.integer(trans)
    if (length(trans) != 1L || is.na(trans) || !trans %in% 1:6) {
      .nm_stop("`trans` must be an integer from 1 through 6.")
    }
    specification$trans <- trans
  }
  if (advan %in% c(6L, 13L) && !is.null(n_state)) {
    n_state <- as.integer(n_state)
    if (length(n_state) != 1L || is.na(n_state) || n_state < 1L || n_state > 20L) {
      .nm_stop("ODE templates support between 1 and 20 compartments.")
    }
    existing <- length(gregexpr("DADT\\(", specification$des, fixed = FALSE)[[1L]])
    if (n_state != existing) {
      specification$des <- paste(c(
        "DADT(1) = -K * A(1)",
        if (n_state > 1L) paste0("DADT(", 2:n_state, ") = 0") else character()
      ), collapse = "\n")
      specification$values <- c(5, 50)
      specification$pred <- paste(
        "CL = THETA(1) * exp(ETA(1))",
        "V = THETA(2) * exp(ETA(2))", "K = CL / V", "S1 = V", sep = "\n"
      )
      specification$eta <- 2L
    }
  }
  observation_compartment <- if (advan %in% c(2L, 4L, 12L)) 2L else 1L
  model <- nm_model(
    INPUT = input, ADVAN = advan, TRANS = specification$trans,
    DOSECMP = 1L, OBSCMP = observation_compartment,
    PRED = specification$pred, DES = specification$des,
    ERROR = "Y = F + ERR(1)", THETAS = theta(specification$values),
    OMEGAS = omega(specification$eta), SIGMAS = sigma
  )
  problem <- trimws(as.character(problem %||% ""))
  attr(model, "name") <- if (nzchar(problem)) problem else paste0("ADVAN", advan, " model")
  model
}

.liber_builtin_dataset <- function(model, example = "theophylline", n_subjects = 10L,
                                   seed = 1L) {
  n_subjects <- max(1L, min(500L, as.integer(n_subjects)))
  times <- switch(
    as.character(example),
    sparse = c(1, 4, 12, 24),
    rich = c(0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 18, 24),
    c(0.5, 1, 2, 3, 4, 6, 8, 12, 24)
  )
  blocks <- lapply(seq_len(n_subjects), function(id) {
    rbind(
      data.frame(
        ID = id, TIME = 0, EVID = 1L, AMT = 320, RATE = 0, II = 0,
        SS = 0L, CMT = model$DOSECMP, DV = NA_real_, MDV = 1L
      ),
      data.frame(
        ID = id, TIME = times, EVID = 0L, AMT = 0, RATE = 0, II = 0,
        SS = 0L, CMT = model$OBSCMP, DV = 0, MDV = 0L
      )
    )
  })
  design <- nm_dataset(do.call(rbind, blocks))
  simulated <- nm_simulate(
    model, design, random_effects = TRUE, residual = TRUE,
    seed = as.integer(seed)
  )
  design$DV <- simulated$DV
  eta_columns <- grep("^ETA[0-9]+$", names(simulated), value = TRUE)
  if (length(eta_columns)) {
    truth <- unique(simulated[c("ID", ".ID_INDEX", eta_columns)])
    attr(design, "simulation_eta") <- as.data.frame(truth)
  }
  attr(design, "name") <- switch(
    as.character(example), sparse = "Sparse oral PK example",
    rich = "Rich oral PK example", "Theophylline-style example"
  )
  design
}

.liber_simulation_dataset <- function(data, model, event) {
  default_subjects <- if (is.null(data)) 1L else {
    source_data <- if (inherits(data, "nm_dataset")) data else nm_dataset(data)
    length(unique(source_data$.ID_INDEX))
  }
  n_subjects <- max(1L, min(10000L, as.integer(event$nSubjects %||% default_subjects)))
  if (isTRUE(event$useDesign)) {
    days <- max(1, min(365, as.numeric(event$days %||% 1)))
    observations_per_day <- max(3L, min(48L, as.integer(event$obsPerDay %||% 8L)))
    amount <- max(0, as.numeric(event$doseAmt %||% 320))
    dose_compartment <- max(1L, as.integer(event$doseCmt %||% model$DOSECMP))
    interval <- max(0.01, as.numeric(event$doseII %||% 12))
    mode <- as.character(event$doseMode %||% "single")
    dose_times <- switch(
      mode,
      `repeat` = seq(0, by = interval, length.out = max(1L, as.integer(event$doseN %||% 3L))),
      steady_state = 0,
      0
    )
    dose_amounts <- rep(amount, length(dose_times))
    dose_table <- trimws(as.character(event$doseTable %||% ""))
    if (nzchar(dose_table)) {
      rows <- strsplit(dose_table, "[\r\n;]+")[[1L]]
      values <- lapply(rows[nzchar(trimws(rows))], function(row) {
        suppressWarnings(as.numeric(strsplit(trimws(row), "[ ,\t]+")[[1L]]))
      })
      if (length(values) && all(vapply(values, function(x) length(x) %in% c(1L, 2L) &&
                                        all(is.finite(x)), logical(1)))) {
        if (all(lengths(values) == 1L)) {
          dose_amounts <- vapply(values, `[[`, numeric(1), 1L)
          dose_times <- if (identical(mode, "single")) 0 else
            seq(0, by = interval, length.out = length(dose_amounts))
        } else if (all(lengths(values) == 2L)) {
          dose_times <- vapply(values, `[[`, numeric(1), 1L)
          dose_amounts <- vapply(values, `[[`, numeric(1), 2L)
        } else {
          .nm_stop("Dose-table rows must consistently contain AMT or TIME AMT.")
        }
      } else if (length(values)) {
        .nm_stop("Dose table must contain finite AMT or TIME AMT values.")
      }
    }
    observation_times <- seq(0, days * 24, length.out = days * observations_per_day + 1L)
    blocks <- lapply(seq_len(n_subjects), function(id) {
      doses <- data.frame(
        ID = id, TIME = dose_times, EVID = 1L, AMT = dose_amounts, RATE = 0,
        II = if (identical(mode, "steady_state")) interval else 0,
        SS = if (identical(mode, "steady_state")) 1L else 0L,
        CMT = dose_compartment, DV = NA_real_, MDV = 1L
      )
      observations <- data.frame(
        ID = id, TIME = observation_times, EVID = 0L, AMT = 0, RATE = 0,
        II = 0, SS = 0L, CMT = model$OBSCMP, DV = 0, MDV = 0L
      )
      rbind(doses, observations)
    })
    frame <- do.call(rbind, blocks)
    frame <- frame[order(frame$ID, frame$TIME, -frame$EVID), , drop = FALSE]
    rownames(frame) <- NULL
    dataset <- nm_dataset(frame)
    attr(dataset, "name") <- "Custom simulation design"
    return(dataset)
  }
  if (is.null(data)) .nm_stop("Load a dataset or enable a custom simulation design.")
  dataset <- if (inherits(data, "nm_dataset")) data else nm_dataset(data)
  source <- as.data.frame(dataset)
  source[grep("^\\.", names(source), value = TRUE)] <- NULL
  source_ids <- unique(source$ID)
  selected <- sample(source_ids, n_subjects, replace = n_subjects > length(source_ids))
  blocks <- lapply(seq_along(selected), function(index) {
    block <- source[source$ID == selected[[index]], , drop = FALSE]
    block$ID <- index
    block
  })
  result <- nm_dataset(do.call(rbind, blocks))
  attr(result, "name") <- attr(dataset, "name", exact = TRUE) %||% "Resampled simulation design"
  result
}

#' Launch the LibeR React modelling application
#'
#' @param model Optional model loaded into the workbench.
#' @param data Optional dataset loaded into the workbench.
#' @param queue Optional `LibeRQueue`; simulation jobs are submitted to it.
#' @param workspace Optional workspace object or directory for versioned project
#'   snapshots.
#' @param project Optional project id to open when the application starts.
#' @param launch.browser Passed to [shiny::runApp()]. Use `NULL` to return the
#'   Shiny app object without launching it.
#' @param ... Additional arguments passed to [shiny::runApp()].
#' @noRd
.liber_gui_legacy <- function(model = NULL, data = NULL, queue = NULL,
                      workspace = NULL, project = NULL,
                      launch.browser = getOption("shiny.launch.browser", interactive()), ...) {
  workspace <- if (is.null(workspace)) NULL else {
    if (inherits(workspace, "nm_workspace")) workspace else nm_workspace(workspace)
  }
  initial_result <- NULL
  if (!is.null(project)) {
    if (is.null(workspace)) .nm_stop("`project` requires a workspace.")
    snapshot <- nm_project_load(workspace, project)
    model <- model %||% snapshot$model
    data <- data %||% snapshot$data
    initial_result <- snapshot$result
  }
  state <- shiny::reactiveValues(
    model = if (inherits(model, "NMEngine")) model$model else model,
    data = data, result = initial_result,
    fit = if (inherits(initial_result, "nm_fit")) initial_result else NULL,
    jobs = data.frame(), project = project
  )
  ui <- .liber_full_page_ui(
    htmltools::tags$head(
      htmltools::tags$title("LibeRation"),
      htmltools::tags$meta(name = "viewport", content = "width=device-width, initial-scale=1")
    ),
    liberWorkbenchOutput("workbench", height = "100vh")
  )
  server <- function(input, output, session) {
    if (!is.null(queue)) {
      shiny::observe({
        shiny::invalidateLater(750, session)
        queue$poll(start = TRUE)
        state$jobs <- queue$list()
      })
    }
    output$workbench <- renderLiberWorkbench({
      server_info <- if (is.null(queue)) {
        list(mode = "local", connected = TRUE, platform = R.version$platform,
             worker = "in-process", isolation = "current R session")
      } else {
        list(mode = "queue", connected = TRUE, platform = R.version$platform,
             worker = paste(queue$max_workers, "callr worker(s)"),
             isolation = paste("tenant", queue$user))
      }
      liber_workbench(state$model, state$data, jobs = state$jobs, result = state$result,
                      server = server_info,
                      workspace = .liber_gui_workspace(workspace, state$project),
                      input_id = "liber_workbench", height = "100vh")
    })
    shiny::observeEvent(input$liber_workbench_event, {
      event <- input$liber_workbench_event
      if (identical(event$action, "validate")) {
        state$result <- tryCatch({
          if (is.null(state$model)) .nm_stop("Load a model before validation.")
          nm_compile(state$model)
          structure(list(), class = "liber_gui_validation")
        }, error = identity)
      }
      if (identical(event$action, "simulate")) {
        if (is.null(state$model) || is.null(state$data)) {
          state$result <- simpleError("Load both a model and dataset before simulation.")
        } else if (!is.null(queue)) {
          if (!requireNamespace("LibeRties", quietly = TRUE)) {
            state$result <- simpleError("LibeRties is required for queued execution.")
          } else {
            id <- queue$submit(LibeRties::ls_job("simulate", state$model, state$data))
            state$result <- structure(list(id = id), class = "liber_gui_queued")
          }
        } else {
          state$result <- tryCatch(nm_simulate(state$model, state$data), error = identity)
        }
      }
      if (identical(event$action, "estimate")) {
        if (is.null(state$model) || is.null(state$data)) {
          state$result <- simpleError("Load both a model and dataset before estimation.")
        } else {
          arguments <- list(
            method = as.character(event$method %||% "FOCEI"),
            maxit = as.integer(event$maxit %||% 200L),
            eta_maxit = as.integer(event$etaMaxit %||% 100L),
            tolerance = as.numeric(event$tolerance %||% 1e-6),
            n_cores = max(1L, as.integer(event$nCores %||% 1L)),
            print_every = max(0L, as.integer(event$printEvery %||% 0L)),
            covariance = isTRUE(event$covariance),
            covariance_type = as.character(event$covarianceType %||% "hessian"),
            covariance_tolerance = as.numeric(event$covarianceTolerance %||% 1e-8),
            covariance_samples = as.integer(event$covarianceSamples %||% 200L),
            covariance_seed = as.integer(event$covarianceSeed %||%
              event$methodSeed %||% 20260713L)
          )
          if (identical(arguments$method, "IMP")) {
            arguments$n_imp <- as.integer(event$nImp %||% 200L)
            arguments$seed <- as.integer(event$methodSeed %||% 20260713L)
          } else if (identical(arguments$method, "GQ")) {
            arguments$gq_order <- as.integer(event$gqOrder %||% 5L)
            arguments$gq_grid <- as.character(event$gqGrid %||% "auto")
            arguments$gq_level <- as.integer(event$gqLevel %||% 3L)
            arguments$gq_adaptive <- isTRUE(event$gqAdaptive)
            arguments$gq_max_points <- as.integer(event$gqMaxPoints %||% 100000L)
          } else if (identical(arguments$method, "SAEM")) {
            arguments$n_iter <- as.integer(event$nIter %||% 200L)
            arguments$burn <- as.integer(event$burn %||% floor(arguments$n_iter / 3))
            arguments$mcmc_steps <- as.integer(event$mcmcSteps %||% 2L)
            arguments$seed <- as.integer(event$methodSeed %||% 20260713L)
          } else if (identical(arguments$method, "BAYES")) {
            arguments$n_burn <- as.integer(event$nBurn %||% 500L)
            arguments$n_sample <- as.integer(event$nSample %||% 1000L)
            arguments$n_thin <- as.integer(event$nThin %||% 1L)
            arguments$seed <- as.integer(event$methodSeed %||% 20260713L)
          } else if (arguments$method %in% c("HMC", "NUTS")) {
            arguments$n_warmup <- as.integer(event$nBurn %||% 500L)
            arguments$n_sample <- as.integer(event$nSample %||% 1000L)
            arguments$n_thin <- as.integer(event$nThin %||% 1L)
            arguments$n_chains <- as.integer(event$nChains %||% 4L)
            arguments$target_acceptance <- as.numeric(event$targetAcceptance %||% 0.8)
            arguments$max_depth <- as.integer(event$maxTreeDepth %||% 10L)
            arguments$n_leapfrog <- as.integer(event$nLeapfrog %||% 10L)
            arguments$seed <- as.integer(event$methodSeed %||% 20260719L)
          } else if (arguments$method %in% c("NPML", "NPAG")) {
            arguments$np_points <- as.integer(event$npPoints %||% 25L)
            arguments$np_cycles <- as.integer(event$npCycles %||% 3L)
            arguments$np_max_support <- as.integer(event$npMaxSupport %||% 100L)
            arguments$np_grid_step <- as.numeric(event$npGridStep %||% 1)
            arguments$seed <- as.integer(event$methodSeed %||% 20260719L)
          }
          if (!is.null(queue)) {
            if (!requireNamespace("LibeRties", quietly = TRUE)) {
              state$result <- simpleError("LibeRties is required for queued execution.")
            } else {
              id <- queue$submit(LibeRties::ls_job(
                "estimate", state$model, state$data, arguments = arguments,
                label = paste(arguments$method, "estimation")
              ))
              state$result <- structure(list(id = id), class = "liber_gui_queued")
            }
          } else {
            state$result <- tryCatch(do.call(
              nm_est, c(list(model = state$model, data = state$data), arguments)
            ), error = identity)
            if (inherits(state$result, "nm_fit")) state$fit <- state$result
          }
        }
      }
      if (identical(event$action, "load_csv")) {
        state$result <- tryCatch({
          connection <- textConnection(as.character(event$text %||% ""))
          on.exit(close(connection), add = TRUE)
          state$data <- nm_dataset(utils::read.csv(connection, check.names = FALSE))
          structure(list(), class = "liber_gui_validation")
        }, error = identity)
      }
      if (identical(event$action, "update_model")) {
        state$result <- tryCatch({
          if (is.null(state$model)) .nm_stop("Load a model before editing source.")
          arguments <- state$model[intersect(names(state$model), names(formals(nm_model)))]
          arguments$PRED <- as.character(event$pred %||% state$model$PRED)
          arguments$ERROR <- as.character(event$error %||% state$model$ERROR)
          arguments$DES <- as.character(event$des %||% state$model$DES)
          state$model <- do.call(nm_model, arguments)
          nm_compile(state$model)
          structure(list(), class = "liber_gui_validation")
        }, error = identity)
      }
      if (identical(event$action, "update_parameters")) {
        state$result <- tryCatch({
          if (is.null(state$model)) .nm_stop("Load a model before editing parameters.")
          update_table <- function(table, values, label) {
            if (is.null(values)) return(table)
            if (is.data.frame(values)) values <- lapply(seq_len(nrow(values)), function(i) as.list(values[i, ]))
            if (length(values) != nrow(table)) .nm_stop(label, " parameter count changed unexpectedly.")
            for (i in seq_len(nrow(table))) {
              table$Value[[i]] <- as.numeric(values[[i]]$Value)
              table$FIX[[i]] <- isTRUE(values[[i]]$FIX)
            }
            table
          }
          arguments <- state$model[intersect(names(state$model), names(formals(nm_model)))]
          arguments$THETAS <- update_table(state$model$THETAS, event$theta, "THETA")
          arguments$OMEGAS <- update_table(state$model$OMEGAS, event$omega, "OMEGA")
          arguments$SIGMAS <- update_table(state$model$SIGMAS, event$sigma, "SIGMA")
          state$model <- do.call(nm_model, arguments)
          nm_compile(state$model)
          structure(list(), class = "liber_gui_validation")
        }, error = identity)
      }
      if (identical(event$action, "job_cancel") && !is.null(queue)) {
        state$result <- tryCatch({
          queue$cancel(as.character(event$id)); state$jobs <- queue$poll();
          structure(list(), class = "liber_gui_validation")
        }, error = identity)
      }
      if (identical(event$action, "job_result") && !is.null(queue)) {
        state$result <- tryCatch(queue$result(as.character(event$id)), error = identity)
        if (inherits(state$result, "nm_fit")) state$fit <- state$result
      }
      if (identical(event$action, "vpc")) {
        state$result <- tryCatch({
          if (!inherits(state$fit, "nm_fit")) .nm_stop("Run or retrieve an estimation before VPC.")
          nm_vpc(state$fit, nsim = as.integer(event$nsim %||% 100L))
        }, error = identity)
      }
      if (identical(event$action, "project_create")) {
        state$result <- tryCatch({
          if (is.null(workspace)) .nm_stop("Launch with `workspace` to manage projects.")
          created <- nm_project_create(workspace, as.character(event$name %||% "New project"))
          state$project <- created$id
          structure(list(), class = "liber_gui_validation")
        }, error = identity)
      }
      if (identical(event$action, "project_save")) {
        state$result <- tryCatch({
          if (is.null(workspace) || is.null(state$project)) {
            .nm_stop("Create or open a project before saving a snapshot.")
          }
          nm_project_save(
            workspace, state$project, state$model, state$data,
            state$fit %||% state$result,
            label = as.character(event$label %||% "Workbench snapshot")
          )
          structure(list(), class = "liber_gui_validation")
        }, error = identity)
      }
      if (identical(event$action, "project_open")) {
        state$result <- tryCatch({
          if (is.null(workspace)) .nm_stop("Launch with `workspace` to open projects.")
          opened <- nm_project_load(workspace, as.character(event$id))
          state$model <- opened$model
          state$data <- opened$data
          state$fit <- if (inherits(opened$result, "nm_fit")) opened$result else NULL
          state$project <- opened$project
          opened$result %||% structure(list(), class = "liber_gui_validation")
        }, error = identity)
      }
    }, ignoreInit = TRUE)
  }
  app <- shiny::shinyApp(ui, server)
  if (is.null(launch.browser)) return(app)
  shiny::runApp(app, launch.browser = launch.browser, ...)
}
