.liber_scheduler_scalar <- function(value, name, default = "", empty = TRUE,
                                    maximum = 4096L) {
  value <- as.character(value %||% default)
  if (length(value) != 1L || is.na(value) || grepl("[\r\n]", value) ||
      nchar(value, type = "bytes") > maximum ||
      (!isTRUE(empty) && !nzchar(trimws(value)))) {
    .nm_stop("`", name, "` must be a valid single-line value.")
  }
  trimws(value)
}

.liber_scheduler_integer <- function(value, name, minimum = 1L, maximum = 100000L) {
  numeric <- suppressWarnings(as.numeric(value))
  integer <- suppressWarnings(as.integer(numeric))
  if (length(numeric) != 1L || !is.finite(numeric) || numeric != integer ||
      integer < minimum || integer > maximum) {
    .nm_stop("`", name, "` must be an integer between ", minimum, " and ", maximum, ".")
  }
  integer
}

.liber_scheduler_component <- function(value, name) {
  value <- .liber_scheduler_scalar(value, name, empty = FALSE, maximum = 128L)
  if (!grepl("^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$", value)) {
    .nm_stop("`", name, "` may contain only ASCII letters, digits, '.', '_', or '-'.")
  }
  value
}

.liber_direct_scheduler_normalize <- function(config, check_identity = TRUE) {
  if (!is.list(config)) .nm_stop("Direct scheduler settings must be a named list.")
  backend <- tolower(.liber_scheduler_scalar(
    config$backend, "scheduler backend", default = "slurm", empty = FALSE
  ))
  if (!backend %in% c("slurm", "grid_engine")) {
    .nm_stop("Scheduler backend must be `slurm` or `grid_engine`.")
  }
  ssh <- .liber_ssh_tunnel_normalize(
    utils::modifyList(list(
      remote_host = "127.0.0.1", remote_port = 8000L, local_port = 0L,
      auto_start = FALSE
    ), config$ssh %||% list()),
    check_identity = check_identity
  )
  optional <- function(name) .liber_scheduler_scalar(config[[name]], name)
  saved_limits <- config$limits %||% list()
  if (!is.list(saved_limits)) .nm_stop("Direct scheduler limits must be a named list.")
  storage_key <- tolower(.liber_scheduler_scalar(
    config$storage_key, "storage key", empty = FALSE, maximum = 64L
  ))
  if (!grepl("^[a-f0-9]{64}$", storage_key)) {
    .nm_stop("The direct scheduler storage key must contain 64 hexadecimal characters.")
  }
  list(
    backend = backend,
    queue_name = .liber_scheduler_component(config$queue_name %||% "default", "queue name"),
    root = .liber_scheduler_scalar(config$root, "remote queue root"),
    remote_rscript = .liber_scheduler_scalar(
      config$remote_rscript, "remote Rscript", default = "Rscript", empty = FALSE
    ),
    max_workers = .liber_scheduler_integer(
      config$max_workers %||% 8L, "maximum concurrent scheduler jobs", 1L, 10000L
    ),
    max_cores_per_job = .liber_scheduler_integer(
      config$max_cores_per_job %||% 64L, "maximum cores per job", 1L, 100000L
    ),
    partition = optional("partition"), account = optional("account"),
    qos = optional("qos"), constraint = optional("constraint"),
    queue = optional("queue"), project = optional("project"),
    parallel_environment = .liber_scheduler_scalar(
      config$parallel_environment, "parallel environment", default = "smp", empty = FALSE
    ),
    memory_resource = .liber_scheduler_scalar(
      config$memory_resource, "memory resource", default = "mem", empty = FALSE
    ),
    runtime_resource = .liber_scheduler_scalar(
      config$runtime_resource, "runtime resource", default = "h_rt", empty = FALSE
    ),
    tmpfs_resource = .liber_scheduler_scalar(
      config$tmpfs_resource, "tmpfs resource", default = "tmpfs", empty = FALSE
    ),
    memory_per_core = !identical(config$memory_per_core, FALSE),
    tmpfs_mb = if (is.null(config$tmpfs_mb) || !length(config$tmpfs_mb) ||
                     !nzchar(as.character(config$tmpfs_mb[[1L]]))) NULL else
      .liber_scheduler_integer(config$tmpfs_mb, "tmpfs MB", 1L, 100000000L),
    limits = list(
      max_concurrent_jobs = .liber_scheduler_integer(
        config$max_workers %||% 8L, "maximum concurrent scheduler jobs", 1L, 10000L
      ),
      max_queued_jobs = .liber_scheduler_integer(
        config$max_queued_jobs %||% saved_limits$max_queued_jobs %||% 100L,
        "maximum queued jobs", 1L, 1000000L
      ),
      max_runtime_seconds = .liber_scheduler_integer(
        config$max_runtime_seconds %||% saved_limits$max_runtime_seconds %||% 86400L,
        "wall-time limit", 1L, 31536000L
      ),
      max_cpu_seconds = .liber_scheduler_integer(
        config$max_cpu_seconds %||% saved_limits$max_cpu_seconds %||% 604800L,
        "CPU-time limit", 1L, 315360000L
      ),
      max_memory_mb = .liber_scheduler_integer(
        config$max_memory_mb %||% saved_limits$max_memory_mb %||% 4096L,
        "memory limit", 1L, 100000000L
      )
    ),
    user = .liber_scheduler_component(config$user %||% "local", "queue user"),
    storage_key = storage_key,
    ssh = ssh
  )
}

.liber_remote_shell_quote <- function(value) {
  paste0("'", gsub("'", "'\\''", as.character(value), fixed = TRUE), "'")
}

.liber_ssh_scheduler_command <- function(config, ssh = .liber_ssh_tool("ssh")) {
  config <- .liber_direct_scheduler_normalize(config)
  ssh <- .liber_scheduler_scalar(ssh, "OpenSSH executable", empty = FALSE)
  ssh_config_file <- .liber_ssh_jump_config(config$ssh)
  connection <- config$ssh
  args <- c(
    if (nzchar(ssh_config_file)) c("-F", ssh_config_file) else character(),
    "-T", "-o", "BatchMode=yes", "-o", "ConnectTimeout=15",
    "-o", "ServerAliveInterval=30", "-o", "ServerAliveCountMax=3",
    "-o", paste0(
      "StrictHostKeyChecking=", if (connection$accept_new_host_key) "accept-new" else "yes"
    ),
    "-p", as.character(connection$port)
  )
  if (nzchar(connection$identity_file) && !nzchar(ssh_config_file)) {
    args <- c(args, "-i", connection$identity_file, "-o", "IdentitiesOnly=yes")
  }
  if (nzchar(connection$proxy_host) && !nzchar(ssh_config_file)) {
    proxy_user <- connection$proxy_user
    if (!nzchar(proxy_user)) proxy_user <- connection$user
    jump <- paste0(proxy_user, "@", connection$proxy_host)
    if (connection$proxy_port != 22L) jump <- paste0(jump, ":", connection$proxy_port)
    args <- c(args, "-J", jump)
  }
  target <- if (nzchar(ssh_config_file)) "liber-managed-target" else
    paste0(connection$user, "@", connection$host)
  remote_command <- paste(
    "exec", .liber_remote_shell_quote(config$remote_rscript), "--vanilla -e",
    .liber_remote_shell_quote("LibeRties::ls_direct_scheduler_cli()")
  )
  list(command = ssh, args = c(args, target, remote_command),
       ssh_config_file = ssh_config_file, config = config)
}

.liber_ssh_scheduler_envelope <- function(payload) {
  payload_json <- unname(as.character(jsonlite::toJSON(
    payload, auto_unbox = TRUE, null = "null", digits = 17, force = TRUE
  )))
  list(
    schema = "liberties.ssh.request", version = 1L,
    payload_json = payload_json,
    sha256 = digest::digest(charToRaw(enc2utf8(payload_json)), algo = "sha256", serialize = FALSE)
  )
}

.liber_ssh_scheduler_response <- function(stdout) {
  lines <- strsplit(enc2utf8(as.character(stdout %||% "")), "\n", fixed = TRUE)[[1L]]
  begin <- which(trimws(lines) == "LIBERTIES_SSH_RESPONSE_BEGIN")
  end <- which(trimws(lines) == "LIBERTIES_SSH_RESPONSE_END")
  if (!length(begin) || !length(end) || end[[length(end)]] <= begin[[1L]] + 1L) {
    .nm_stop("The remote scheduler returned no valid LibeRties response.")
  }
  encoded <- paste(lines[seq.int(begin[[1L]] + 1L, end[[length(end)]] - 1L)], collapse = "\n")
  envelope <- tryCatch(
    jsonlite::fromJSON(encoded, simplifyVector = FALSE),
    error = function(error) .nm_stop("Invalid remote scheduler response: ", conditionMessage(error))
  )
  if (!identical(as.character(envelope$schema %||% ""), "liberties.ssh.response") ||
      !identical(as.integer(envelope$version %||% 0L), 1L)) {
    .nm_stop("Unsupported remote scheduler response envelope.")
  }
  payload_json <- as.character(envelope$payload_json %||% "")
  actual <- digest::digest(charToRaw(enc2utf8(payload_json)), algo = "sha256", serialize = FALSE)
  if (!identical(tolower(as.character(envelope$sha256 %||% "")), actual)) {
    .nm_stop("Remote scheduler response checksum mismatch.")
  }
  response <- jsonlite::fromJSON(payload_json, simplifyVector = FALSE)
  if (!isTRUE(response$ok)) .nm_stop(as.character(response$error %||% "Remote scheduler request failed."))
  response$payload
}

.liber_ssh_scheduler_call <- function(client, operation, fields = list()) {
  if (!requireNamespace("processx", quietly = TRUE)) {
    .nm_stop("Direct SSH scheduler execution requires the `processx` package.")
  }
  command <- .liber_ssh_scheduler_command(client$config)
  if (!nzchar(command$command) || !file.exists(command$command)) {
    .nm_stop("OpenSSH was not found. Complete SSH readiness setup first.")
  }
  on.exit({
    if (nzchar(command$ssh_config_file) && file.exists(command$ssh_config_file)) {
      unlink(command$ssh_config_file, force = TRUE)
    }
  }, add = TRUE)
  payload <- c(list(
    schema = "liberties.direct-scheduler", version = 1L,
    operation = operation,
    config = command$config[setdiff(names(command$config), c("ssh", "storage_key"))],
    storage_key = command$config$storage_key
  ), fields)
  input <- tempfile("liber-scheduler-request-", fileext = ".json")
  on.exit(unlink(input, force = TRUE), add = TRUE)
  writeLines(jsonlite::toJSON(
    .liber_ssh_scheduler_envelope(payload), auto_unbox = TRUE, null = "null",
    digits = 17, force = TRUE
  ), input, useBytes = TRUE)
  result <- processx::run(
    command$command, command$args, stdin = input, error_on_status = FALSE,
    timeout = as.numeric(client$timeout) * 1000,
    windows_hide_window = TRUE, encoding = "UTF-8"
  )
  if (!identical(as.integer(result$status), 0L)) {
    detail <- trimws(paste(result$stderr %||% "", result$stdout %||% ""))
    .nm_stop(if (nzchar(detail)) detail else "The SSH scheduler request failed.")
  }
  .liber_ssh_scheduler_response(result$stdout)
}

LibeRDirectScheduler <- R6::R6Class(
  "LibeRDirectScheduler",
  public = list(
    config = NULL, timeout = NULL, url = NULL,
    initialize = function(config, timeout = 60) {
      self$config <- .liber_direct_scheduler_normalize(config, check_identity = FALSE)
      self$timeout <- as.numeric(timeout)
      if (length(self$timeout) != 1L || !is.finite(self$timeout) || self$timeout <= 0) {
        .nm_stop("Direct scheduler timeout must be a positive number of seconds.")
      }
      self$url <- paste0(
        "ssh://", self$config$ssh$user, "@", self$config$ssh$host,
        if (self$config$ssh$port == 22L) "" else paste0(":", self$config$ssh$port),
        "/", self$config$backend
      )
    },
    authenticate = function() .liber_ssh_scheduler_call(self, "authenticate"),
    capabilities = function() .liber_ssh_scheduler_call(self, "capabilities"),
    submit = function(job, idempotency_key = NULL) {
      payload <- .liber_ssh_scheduler_call(self, "submit", list(
        job = LibeRties::ls_job_to_wire(job), idempotency_key = idempotency_key
      ))
      as.character(payload$id)
    },
    list = function() {
      jobs <- .liber_ssh_scheduler_call(self, "list")$jobs
      if (!length(jobs)) return(data.frame(
        id = character(), user = character(), type = character(), label = character(),
        status = character(), submitted = character(), started = character(),
        finished = character(), stringsAsFactors = FALSE
      ))
      result <- do.call(rbind, lapply(jobs, function(item) {
        as.data.frame(lapply(item, unlist, use.names = FALSE), stringsAsFactors = FALSE)
      }))
      rownames(result) <- NULL
      result
    },
    status = function(id) .liber_ssh_scheduler_call(self, "status", list(id = id)),
    result = function(id) {
      LibeRties::ls_result_from_wire(
        .liber_ssh_scheduler_call(self, "result", list(id = id))
      )
    },
    logs = function(id, stream = c("stdout", "stderr")) {
      stream <- match.arg(stream)
      as.character(unlist(
        .liber_ssh_scheduler_call(self, "logs", list(id = id, stream = stream))$lines,
        use.names = FALSE
      ))
    },
    cancel = function(id) isTRUE(
      .liber_ssh_scheduler_call(self, "cancel", list(id = id))$cancelled
    )
  )
)

.liber_direct_scheduler <- function(config, timeout = 60) {
  LibeRDirectScheduler$new(config, timeout = timeout)
}

.liber_direct_scheduler_public_config <- function(config) {
  if (!is.list(config)) return(NULL)
  config[setdiff(names(config), "storage_key")]
}
