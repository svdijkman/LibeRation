.liber_ssh_scalar <- function(value, name, default = NULL, empty = FALSE) {
  value <- as.character(value %||% default %||% "")
  if (length(value) != 1L || is.na(value) || grepl("[\r\n]", value) ||
      (!isTRUE(empty) && !nzchar(trimws(value)))) {
    .nm_stop("`", name, "` must be one ", if (empty) "single-line" else "non-empty single-line", " value.")
  }
  trimws(value)
}

.liber_ssh_port <- function(value, name, allow_auto = FALSE) {
  value <- suppressWarnings(as.integer(value))
  minimum <- if (isTRUE(allow_auto)) 0L else 1L
  if (length(value) != 1L || is.na(value) || value < minimum || value > 65535L) {
    .nm_stop(
      "`", name, "` must be an integer between ", minimum,
      " and 65535", if (allow_auto) " (0 selects an available port automatically)" else "", "."
    )
  }
  value
}

.liber_ssh_host <- function(value, name, empty = FALSE) {
  value <- .liber_ssh_scalar(value, name, empty = empty)
  if (!nzchar(value) && isTRUE(empty)) return("")
  if (!grepl("^[A-Za-z0-9][A-Za-z0-9._-]{0,252}$", value)) {
    .nm_stop("`", name, "` must be a DNS name or IPv4-style host name without spaces or options.")
  }
  value
}

.liber_ssh_user <- function(value, name = "SSH user", empty = FALSE) {
  value <- .liber_ssh_scalar(value, name, empty = empty)
  if (!nzchar(value) && isTRUE(empty)) return("")
  if (!grepl("^[A-Za-z_][A-Za-z0-9_.-]{0,127}$", value)) {
    .nm_stop("`", name, "` is not a portable SSH account name.")
  }
  value
}

.liber_ssh_tunnel_normalize <- function(config, check_identity = TRUE) {
  if (!is.list(config)) .nm_stop("SSH tunnel settings must be a named list.")
  identity <- .liber_ssh_scalar(
    config$identity_file, "SSH identity file", default = "", empty = TRUE
  )
  if (nzchar(identity)) {
    identity <- normalizePath(path.expand(identity), winslash = "/", mustWork = FALSE)
    if (isTRUE(check_identity) && !file.exists(identity)) {
      .nm_stop("The selected SSH identity file does not exist: ", identity)
    }
  }
  proxy_host <- .liber_ssh_host(
    config$proxy_host, "SSH jump host", empty = TRUE
  )
  proxy_user <- .liber_ssh_user(
    config$proxy_user, "SSH jump user", empty = TRUE
  )
  if (!nzchar(proxy_host)) proxy_user <- ""
  list(
    host = .liber_ssh_host(config$host, "SSH host"),
    user = .liber_ssh_user(config$user),
    port = .liber_ssh_port(config$port %||% 22L, "SSH port"),
    remote_host = .liber_ssh_host(
      config$remote_host %||% "127.0.0.1", "remote LibeRties host"
    ),
    remote_port = .liber_ssh_port(
      config$remote_port %||% 8000L, "remote LibeRties port"
    ),
    local_port = .liber_ssh_port(
      config$local_port %||% 0L, "local forwarded port", allow_auto = TRUE
    ),
    proxy_host = proxy_host,
    proxy_user = proxy_user,
    proxy_port = .liber_ssh_port(
      config$proxy_port %||% 22L, "SSH jump-host port"
    ),
    identity_file = identity,
    accept_new_host_key = !identical(config$accept_new_host_key, FALSE),
    auto_start = !identical(config$auto_start, FALSE)
  )
}

.liber_ssh_platform <- function() {
  if (.Platform$OS.type == "windows") return("windows")
  if (identical(Sys.info()[["sysname"]], "Darwin")) return("macos")
  "linux"
}

.liber_ssh_tool <- function(name) {
  name <- as.character(name)[[1L]]
  if (!name %in% c("ssh", "ssh-add", "ssh-keygen", "ssh-copy-id", "powershell")) {
    .nm_stop("Unsupported SSH helper executable: ", name)
  }
  executable <- if (.Platform$OS.type == "windows" &&
                    !grepl("[.]exe$", name, ignore.case = TRUE)) {
    paste0(name, ".exe")
  } else name
  candidates <- unname(Sys.which(c(name, executable)))
  if (.Platform$OS.type == "windows") {
    windows <- Sys.getenv("WINDIR", unset = "C:/Windows")
    program_files <- Sys.getenv("ProgramFiles", unset = "C:/Program Files")
    candidates <- c(
      candidates,
      file.path(windows, "System32", "OpenSSH", executable),
      file.path(windows, "System32", "WindowsPowerShell", "v1.0", executable),
      file.path(program_files, "OpenSSH", executable)
    )
  }
  candidates <- unique(candidates[nzchar(candidates)])
  existing <- candidates[file.exists(candidates)]
  if (!length(existing)) return("")
  normalizePath(existing[[1L]], winslash = "/", mustWork = TRUE)
}

.liber_ssh_version <- function(ssh = .liber_ssh_tool("ssh")) {
  if (!nzchar(ssh) || !requireNamespace("processx", quietly = TRUE)) return("")
  result <- tryCatch(
    processx::run(ssh, "-V", error_on_status = FALSE, timeout = 5,
                  windows_hide_window = TRUE),
    error = function(error) NULL
  )
  if (is.null(result)) return("")
  trimws(paste(c(result$stderr, result$stdout), collapse = " "))
}

.liber_ssh_public_key <- function(identity_file) {
  identity_file <- as.character(identity_file %||% "")[[1L]]
  if (!nzchar(identity_file)) return("")
  public <- paste0(path.expand(identity_file), ".pub")
  if (!file.exists(public)) return("")
  value <- trimws(paste(readLines(public, warn = FALSE, encoding = "UTF-8"),
                        collapse = " "))
  if (length(value) != 1L || !grepl(
    "^(ssh-(ed25519|rsa)|ecdsa-sha2-nistp[0-9]+) [A-Za-z0-9+/=]+(?: |$)", value,
    perl = TRUE
  )) return("")
  value
}

.liber_ssh_key_fingerprint <- function(identity_file) {
  keygen <- .liber_ssh_tool("ssh-keygen")
  public <- paste0(identity_file, ".pub")
  if (!nzchar(keygen) || !file.exists(public) ||
      !requireNamespace("processx", quietly = TRUE)) return("")
  result <- tryCatch(
    processx::run(keygen, c("-lf", public), error_on_status = FALSE,
                  timeout = 5, windows_hide_window = TRUE),
    error = function(error) NULL
  )
  if (is.null(result) || !identical(result$status, 0L)) return("")
  trimws(result$stdout)
}

.liber_ssh_key_candidates <- function(home = path.expand("~"), include = "") {
  directory <- file.path(home, ".ssh")
  paths <- if (dir.exists(directory)) {
    list.files(directory, full.names = TRUE, all.files = FALSE)
  } else character()
  paths <- paths[
    file.exists(paths) & !dir.exists(paths) &
      !grepl("[.](pub|cert)$", paths, ignore.case = TRUE) &
      (grepl("^id_", basename(paths)) | grepl("[.](pem|key)$", paths, ignore.case = TRUE))
  ]
  include <- as.character(include %||% "")[[1L]]
  if (nzchar(include) && file.exists(path.expand(include))) {
    paths <- c(path.expand(include), paths)
  }
  paths <- unique(normalizePath(paths, winslash = "/", mustWork = TRUE))
  if (!length(paths)) return(list())
  preferred <- c("id_ed25519_liber", "id_ed25519", "id_ecdsa", "id_rsa")
  order_key <- match(basename(paths), preferred, nomatch = length(preferred) + 1L)
  paths <- paths[order(order_key, basename(paths), method = "radix")]
  unname(lapply(paths, function(path) list(
    path = path,
    name = basename(path),
    fingerprint = .liber_ssh_key_fingerprint(path),
    public_key = .liber_ssh_public_key(path),
    has_public_key = nzchar(.liber_ssh_public_key(path))
  )))
}

.liber_ssh_agent_status <- function() {
  add <- .liber_ssh_tool("ssh-add")
  if (!nzchar(add) || !requireNamespace("processx", quietly = TRUE)) {
    return(list(status = "unavailable", message = "ssh-add is unavailable.", identities = list()))
  }
  result <- tryCatch(
    processx::run(add, "-l", error_on_status = FALSE, timeout = 5,
                  windows_hide_window = TRUE),
    error = function(error) NULL
  )
  if (is.null(result)) {
    return(list(status = "unavailable", message = "Unable to query ssh-agent.", identities = list()))
  }
  detail <- trimws(paste(c(result$stdout, result$stderr), collapse = " "))
  if (identical(result$status, 0L)) {
    identities <- Filter(nzchar, trimws(strsplit(result$stdout, "\n", fixed = TRUE)[[1L]]))
    return(list(status = "ready", message = paste(length(identities), "key(s) loaded."),
                identities = as.list(identities)))
  }
  if (identical(result$status, 1L)) {
    return(list(status = "empty", message = "ssh-agent is running but has no keys loaded.",
                identities = list()))
  }
  list(status = "unavailable", message = if (nzchar(detail)) detail else "ssh-agent is not running.",
       identities = list())
}

.liber_ssh_readiness <- function(identity_file = "", status = "idle", message = "") {
  ssh <- .liber_ssh_tool("ssh")
  keys <- .liber_ssh_key_candidates(include = identity_file)
  recommended <- if (nzchar(identity_file) && file.exists(path.expand(identity_file))) {
    normalizePath(path.expand(identity_file), winslash = "/", mustWork = TRUE)
  } else if (length(keys)) keys[[1L]]$path else ""
  list(
    platform = .liber_ssh_platform(), available = nzchar(ssh),
    ssh_path = ssh, version = .liber_ssh_version(ssh),
    install_supported = identical(.liber_ssh_platform(), "windows"),
    agent = .liber_ssh_agent_status(), keys = keys,
    recommended_key = recommended,
    status = status, message = as.character(message %||% "")[[1L]],
    nonce = as.numeric(Sys.time())
  )
}

.liber_ssh_powershell_encode <- function(script) {
  raw <- iconv(enc2utf8(script), from = "UTF-8", to = "UTF-16LE", toRaw = TRUE)[[1L]]
  jsonlite::base64_enc(raw)
}

.liber_ssh_windows_launch <- function(script, elevated = FALSE) {
  if (.Platform$OS.type != "windows") .nm_stop("This setup action requires Windows.")
  powershell <- .liber_ssh_tool("powershell")
  if (!nzchar(powershell) || !requireNamespace("processx", quietly = TRUE)) {
    .nm_stop("Windows PowerShell is required for this setup action.")
  }
  encoded <- .liber_ssh_powershell_encode(script)
  verb <- if (isTRUE(elevated)) " -Verb RunAs" else ""
  launcher <- paste0(
    "$ErrorActionPreference='Stop'; Start-Process -FilePath 'powershell.exe'",
    verb,
    " -ArgumentList @('-NoProfile','-ExecutionPolicy','Bypass','-EncodedCommand','",
    encoded, "') | Out-Null"
  )
  result <- processx::run(
    powershell,
    c("-NoProfile", "-NonInteractive", "-EncodedCommand",
      .liber_ssh_powershell_encode(launcher)),
    error_on_status = FALSE, timeout = 20, windows_hide_window = TRUE
  )
  if (!identical(result$status, 0L)) {
    detail <- trimws(paste(c(result$stderr, result$stdout), collapse = " "))
    .nm_stop(if (nzchar(detail)) detail else "Windows declined the setup action.")
  }
  invisible(TRUE)
}

.liber_ssh_install_client <- function() {
  if (.Platform$OS.type != "windows") {
    .nm_stop(
      "Automatic OpenSSH installation is currently available on Windows. On Debian/Ubuntu run ",
      "`sudo apt-get install openssh-client`; on Fedora/RHEL run `sudo dnf install openssh-clients`."
    )
  }
  .liber_ssh_windows_launch(paste(
    "$ErrorActionPreference='Stop'",
    "$cap=Get-WindowsCapability -Online -Name 'OpenSSH.Client~~~~0.0.1.0'",
    "if ($cap.State -ne 'Installed') { Add-WindowsCapability -Online -Name 'OpenSSH.Client~~~~0.0.1.0' | Out-Null }",
    "$agent=Get-Service -Name 'ssh-agent' -ErrorAction SilentlyContinue",
    "if ($null -ne $agent) { Set-Service -Name 'ssh-agent' -StartupType Automatic; Start-Service -Name 'ssh-agent' }",
    sep = "; "
  ), elevated = TRUE)
  invisible(TRUE)
}

.liber_ssh_enable_agent <- function() {
  if (.Platform$OS.type != "windows") {
    .nm_stop("Start the desktop ssh-agent for this login session, then click Refresh readiness.")
  }
  .liber_ssh_windows_launch(paste(
    "$ErrorActionPreference='Stop'",
    "Set-Service -Name 'ssh-agent' -StartupType Automatic",
    "Start-Service -Name 'ssh-agent'",
    sep = "; "
  ), elevated = TRUE)
  invisible(TRUE)
}

.liber_ssh_default_identity <- function() {
  normalizePath(file.path(path.expand("~"), ".ssh", "id_ed25519_liber"),
                winslash = "/", mustWork = FALSE)
}

.liber_ssh_generate_key <- function(path = .liber_ssh_default_identity()) {
  keygen <- .liber_ssh_tool("ssh-keygen")
  if (!nzchar(keygen)) .nm_stop("Install OpenSSH before generating a key.")
  path <- normalizePath(path.expand(path), winslash = "/", mustWork = FALSE)
  if (file.exists(path) || file.exists(paste0(path, ".pub"))) {
    .nm_stop("A key already exists at ", path, ". Select it instead of overwriting it.")
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE, mode = "0700")
  if (.Platform$OS.type != "windows") {
    .nm_stop(
      "Run `ssh-keygen -t ed25519 -a 64 -f ", path,
      " -C \"LibeR remote access\"` in a terminal, then refresh readiness."
    )
  }
  quote_ps <- function(value) paste0("'", gsub("'", "''", value, fixed = TRUE), "'")
  script <- paste0(
    "$ErrorActionPreference='Stop'; & ", quote_ps(keygen),
    " -t ed25519 -a 64 -f ", quote_ps(path),
    " -C 'LibeR remote access'; if ($LASTEXITCODE -ne 0) { throw 'ssh-keygen failed.' }; ",
    "Write-Host ''; Write-Host 'Key created. Return to LibeRation and refresh SSH readiness.'; ",
    "Read-Host 'Press Enter to close'"
  )
  .liber_ssh_windows_launch(script, elevated = FALSE)
  invisible(path)
}

.liber_ssh_load_key <- function(identity_file) {
  identity_file <- normalizePath(path.expand(identity_file), winslash = "/", mustWork = TRUE)
  add <- .liber_ssh_tool("ssh-add")
  if (!nzchar(add)) .nm_stop("ssh-add is unavailable.")
  if (.Platform$OS.type != "windows") {
    .nm_stop("Run `ssh-add ", shQuote(identity_file), "` in a terminal, then refresh readiness.")
  }
  quote_ps <- function(value) paste0("'", gsub("'", "''", value, fixed = TRUE), "'")
  script <- paste0(
    "$ErrorActionPreference='Stop'; & ", quote_ps(add), " ", quote_ps(identity_file),
    "; if ($LASTEXITCODE -ne 0) { throw 'ssh-add failed.' }; ",
    "Write-Host 'Key loaded into ssh-agent.'; Read-Host 'Press Enter to close'"
  )
  .liber_ssh_windows_launch(script, elevated = FALSE)
  invisible(TRUE)
}

.liber_ssh_bootstrap_config <- function(config) {
  if (!nzchar(config$proxy_host) || !nzchar(config$identity_file)) return("")
  proxy_user <- if (nzchar(config$proxy_user)) config$proxy_user else config$user
  path <- tempfile("liber-ssh-bootstrap-", fileext = ".conf")
  strict <- if (config$accept_new_host_key) "accept-new" else "yes"
  identity <- encodeString(config$identity_file, quote = '"')
  writeLines(c(
    "Host liber-bootstrap-target",
    paste("  HostName", config$host), paste("  User", config$user),
    paste("  Port", config$port), paste("  IdentityFile", identity),
    "  IdentitiesOnly yes", "  BatchMode no",
    paste("  StrictHostKeyChecking", strict),
    "  ProxyJump liber-bootstrap-gateway",
    "Host liber-bootstrap-gateway",
    paste("  HostName", config$proxy_host), paste("  User", proxy_user),
    paste("  Port", config$proxy_port), paste("  IdentityFile", identity),
    "  IdentitiesOnly yes", "  BatchMode no",
    paste("  StrictHostKeyChecking", strict)
  ), path, useBytes = TRUE)
  Sys.chmod(path, mode = "0600")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.liber_ssh_test_host <- function(config, hop = c("destination", "gateway"),
                                 timeout = 15) {
  hop <- match.arg(hop)
  config <- .liber_ssh_tunnel_normalize(config)
  if (identical(hop, "gateway") && !nzchar(config$proxy_host)) {
    .nm_stop("Enable and configure a jump host before testing it.")
  }
  ssh <- .liber_ssh_tool("ssh")
  if (!nzchar(ssh)) .nm_stop("OpenSSH is not installed or available.")
  temporary <- ""
  on.exit(if (nzchar(temporary)) unlink(temporary, force = TRUE), add = TRUE)
  strict <- if (config$accept_new_host_key) "accept-new" else "yes"
  args <- c(
    "-T", "-o", "BatchMode=yes", "-o", "ConnectTimeout=8",
    "-o", "ConnectionAttempts=1", "-o", paste0("StrictHostKeyChecking=", strict)
  )
  if (identical(hop, "gateway")) {
    args <- c(args, "-p", as.character(config$proxy_port))
    if (nzchar(config$identity_file)) {
      args <- c(args, "-i", config$identity_file, "-o", "IdentitiesOnly=yes")
    }
    user <- if (nzchar(config$proxy_user)) config$proxy_user else config$user
    destination <- paste0(user, "@", config$proxy_host)
  } else if (nzchar(config$proxy_host) && nzchar(config$identity_file)) {
    temporary <- .liber_ssh_jump_config(config)
    args <- c("-F", temporary, args)
    destination <- "liber-managed-target"
  } else {
    args <- c(args, "-p", as.character(config$port))
    if (nzchar(config$identity_file)) {
      args <- c(args, "-i", config$identity_file, "-o", "IdentitiesOnly=yes")
    }
    if (nzchar(config$proxy_host)) {
      user <- if (nzchar(config$proxy_user)) config$proxy_user else config$user
      jump <- paste0(user, "@", config$proxy_host)
      if (config$proxy_port != 22L) jump <- paste0(jump, ":", config$proxy_port)
      args <- c(args, "-J", jump)
    }
    destination <- paste0(config$user, "@", config$host)
  }
  result <- processx::run(
    ssh, c(args, destination, "true"), error_on_status = FALSE,
    timeout = as.numeric(timeout), windows_hide_window = TRUE
  )
  if (!identical(result$status, 0L)) {
    detail <- trimws(paste(c(result$stderr, result$stdout), collapse = " "))
    .nm_stop(if (nzchar(detail)) detail else paste("Unable to authenticate to", hop, "host."))
  }
  paste(if (identical(hop, "gateway")) "Gateway" else "Destination",
        "SSH authentication succeeded.")
}

.liber_ssh_install_public_key <- function(config, hop = c("destination", "gateway")) {
  hop <- match.arg(hop)
  config <- .liber_ssh_tunnel_normalize(config)
  if (identical(hop, "gateway") && !nzchar(config$proxy_host)) {
    .nm_stop("Enable and configure a jump host before installing a gateway key.")
  }
  public <- paste0(config$identity_file, ".pub")
  if (!nzchar(config$identity_file) || !file.exists(public) ||
      !nzchar(.liber_ssh_public_key(config$identity_file))) {
    .nm_stop("Select a private key with a valid matching `.pub` file first.")
  }
  ssh <- .liber_ssh_tool("ssh")
  if (!nzchar(ssh)) .nm_stop("OpenSSH is not installed or available.")
  if (.Platform$OS.type != "windows") {
    destination <- if (identical(hop, "gateway")) {
      user <- if (nzchar(config$proxy_user)) config$proxy_user else config$user
      paste0(user, "@", config$proxy_host)
    } else paste0(config$user, "@", config$host)
    .nm_stop(
      "Install the public key from a terminal using `ssh-copy-id -i ", public,
      " ", destination, "`, adding `-o ProxyJump=...` for the destination when required."
    )
  }
  temporary <- ""
  if (identical(hop, "destination") && nzchar(config$proxy_host)) {
    temporary <- .liber_ssh_bootstrap_config(config)
    args <- c("-F", temporary)
    destination <- "liber-bootstrap-target"
  } else {
    target_host <- if (identical(hop, "gateway")) config$proxy_host else config$host
    target_user <- if (identical(hop, "gateway") && nzchar(config$proxy_user)) {
      config$proxy_user
    } else config$user
    target_port <- if (identical(hop, "gateway")) config$proxy_port else config$port
    args <- c(
      "-p", as.character(target_port), "-i", config$identity_file,
      "-o", "IdentitiesOnly=yes", "-o", "BatchMode=no", "-o",
      paste0("StrictHostKeyChecking=", if (config$accept_new_host_key) "accept-new" else "yes")
    )
    destination <- paste0(target_user, "@", target_host)
  }
  quote_ps <- function(value) paste0("'", gsub("'", "''", value, fixed = TRUE), "'")
  args_ps <- paste(vapply(args, quote_ps, character(1)), collapse = ",")
  remote_command <- paste0(
    "umask 077; mkdir -p \"$HOME/.ssh\"; chmod 700 \"$HOME/.ssh\"; ",
    "touch \"$HOME/.ssh/authorized_keys\"; key=$(cat); ",
    "grep -qxF \"$key\" \"$HOME/.ssh/authorized_keys\" || ",
    "printf \"%s\\n\" \"$key\" >> \"$HOME/.ssh/authorized_keys\"; ",
    "chmod 600 \"$HOME/.ssh/authorized_keys\""
  )
  cleanup <- if (nzchar(temporary)) {
    paste0("Remove-Item -LiteralPath ", quote_ps(temporary), " -Force -ErrorAction SilentlyContinue")
  } else ""
  script <- paste0(
    "$ErrorActionPreference='Stop'; $ssh=", quote_ps(ssh), "; $args=@(", args_ps,
    "); $destination=", quote_ps(destination), "; $remote=", quote_ps(remote_command),
    "; try { Get-Content -Raw -LiteralPath ", quote_ps(public),
    " | & $ssh @args $destination $remote; if ($LASTEXITCODE -ne 0) { throw 'Public-key installation failed.' }; ",
    "Write-Host 'Public key installed successfully.' } finally { ", cleanup,
    " }; Read-Host 'Press Enter to close'"
  )
  launched <- tryCatch(
    .liber_ssh_windows_launch(script, elevated = FALSE),
    error = identity
  )
  if (inherits(launched, "error")) {
    if (nzchar(temporary)) unlink(temporary, force = TRUE)
    stop(launched)
  }
  invisible(TRUE)
}

.liber_tcp_ready <- function(port, host = "127.0.0.1", timeout = 0.15) {
  connection <- suppressWarnings(tryCatch(
    socketConnection(
      host = host, port = as.integer(port), open = "r+b",
      blocking = TRUE, timeout = timeout
    ),
    error = function(error) NULL
  ))
  if (is.null(connection)) return(FALSE)
  close(connection)
  TRUE
}

.liber_tcp_bindable <- function(port) {
  listener <- suppressWarnings(tryCatch(
    serverSocket(as.integer(port)),
    error = function(error) NULL
  ))
  if (is.null(listener)) return(FALSE)
  close(listener)
  TRUE
}

.liber_ssh_random_port <- function(attempts = 200L) {
  for (port in sample.int(65535L - 49152L + 1L, attempts, replace = FALSE) + 49151L) {
    # A connect probe alone cannot detect Windows Hyper-V/WSL excluded port
    # ranges: they look unused but fail when OpenSSH tries to bind them.
    if (.liber_tcp_bindable(port)) return(as.integer(port))
  }
  .nm_stop("Unable to find an available local TCP port for the SSH tunnel.")
}

.liber_ssh_jump_config <- function(config) {
  if (!nzchar(config$proxy_host) || !nzchar(config$identity_file)) return("")
  proxy_user <- config$proxy_user
  if (!nzchar(proxy_user)) proxy_user <- config$user
  path <- tempfile("liber-ssh-proxyjump-", fileext = ".conf")
  strict <- if (config$accept_new_host_key) "accept-new" else "yes"
  identity <- encodeString(config$identity_file, quote = '"')
  writeLines(c(
    "Host liber-managed-target",
    paste("  HostName", config$host),
    paste("  User", config$user),
    paste("  Port", config$port),
    paste("  IdentityFile", identity),
    "  IdentitiesOnly yes",
    "  BatchMode yes",
    paste("  StrictHostKeyChecking", strict),
    "  ProxyJump liber-managed-gateway",
    "Host liber-managed-gateway",
    paste("  HostName", config$proxy_host),
    paste("  User", proxy_user),
    paste("  Port", config$proxy_port),
    paste("  IdentityFile", identity),
    "  IdentitiesOnly yes",
    "  BatchMode yes",
    paste("  StrictHostKeyChecking", strict)
  ), path, useBytes = TRUE)
  Sys.chmod(path, mode = "0600")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.liber_ssh_tunnel_command <- function(config, ssh = unname(Sys.which("ssh"))) {
  config <- .liber_ssh_tunnel_normalize(config)
  ssh <- .liber_ssh_scalar(ssh, "OpenSSH executable")
  ssh_config_file <- .liber_ssh_jump_config(config)
  local_port <- config$local_port
  if (!local_port) local_port <- .liber_ssh_random_port()
  forward <- paste0(
    "127.0.0.1:", local_port, ":", config$remote_host, ":", config$remote_port
  )
  args <- c(
    if (nzchar(ssh_config_file)) c("-F", ssh_config_file) else character(),
    "-N", "-T", "-o", "BatchMode=yes", "-o", "ExitOnForwardFailure=yes",
    "-o", "ServerAliveInterval=30", "-o", "ServerAliveCountMax=3",
    "-o", paste0(
      "StrictHostKeyChecking=",
      if (config$accept_new_host_key) "accept-new" else "yes"
    ),
    "-p", as.character(config$port), "-L", forward
  )
  if (nzchar(config$identity_file) && !nzchar(ssh_config_file)) {
    args <- c(args, "-i", config$identity_file, "-o", "IdentitiesOnly=yes")
  }
  if (nzchar(config$proxy_host) && !nzchar(ssh_config_file)) {
    proxy_user <- config$proxy_user
    if (!nzchar(proxy_user)) proxy_user <- config$user
    jump <- paste0(proxy_user, "@", config$proxy_host)
    if (config$proxy_port != 22L) jump <- paste0(jump, ":", config$proxy_port)
    args <- c(args, "-J", jump)
  }
  args <- c(
    args,
    if (nzchar(ssh_config_file)) "liber-managed-target"
    else paste0(config$user, "@", config$host)
  )
  list(command = ssh, args = args, config = config, local_port = local_port,
       url = paste0("http://127.0.0.1:", local_port),
       ssh_config_file = ssh_config_file)
}

.liber_ssh_tunnel_start <- function(config, wait_seconds = 8) {
  if (!requireNamespace("processx", quietly = TRUE)) {
    .nm_stop("Managed SSH tunnels require the `processx` package.")
  }
  ssh <- .liber_ssh_tool("ssh")
  if (!nzchar(ssh) || !file.exists(ssh)) {
    .nm_stop(
      "OpenSSH was not found in PATH. Install the operating-system OpenSSH client ",
      "and verify that `ssh` works in a terminal."
    )
  }
  command <- .liber_ssh_tunnel_command(config, ssh = ssh)
  remove_ssh_config <- function() {
    if (nzchar(command$ssh_config_file) && file.exists(command$ssh_config_file)) {
      unlink(command$ssh_config_file, force = TRUE)
    }
  }
  if (.liber_tcp_ready(command$local_port)) {
    remove_ssh_config()
    .nm_stop("Local port ", command$local_port, " is already in use.")
  }
  process <- NULL
  keep <- FALSE
  on.exit({
    if (!keep && !is.null(process) &&
        isTRUE(tryCatch(process$is_alive(), error = function(error) FALSE))) {
      try(process$kill_tree(), silent = TRUE)
    }
    if (!keep) remove_ssh_config()
  }, add = TRUE)
  process <- processx::process$new(
    command$command, command$args, stdin = NULL, stdout = "|", stderr = "|",
    cleanup = TRUE, cleanup_tree = TRUE, windows_hide_window = TRUE
  )
  deadline <- Sys.time() + max(1, as.numeric(wait_seconds))
  repeat {
    if (.liber_tcp_ready(command$local_port)) break
    if (!isTRUE(process$is_alive())) {
      detail <- trimws(paste(
        tryCatch(process$read_all_error(), error = function(error) ""),
        tryCatch(process$read_all_output(), error = function(error) "")
      ))
      .nm_stop(
        "The SSH tunnel stopped before port forwarding became available",
        if (nzchar(detail)) paste0(": ", detail) else "."
      )
    }
    if (Sys.time() >= deadline) {
      .nm_stop(
        "Timed out waiting for the SSH tunnel. Confirm that key or ssh-agent ",
        "authentication works non-interactively and that the remote port is reachable."
      )
    }
    Sys.sleep(0.05)
  }
  keep <- TRUE
  structure(
    c(command, list(process = process, started = .nm_workspace_now())),
    class = "liber_ssh_tunnel"
  )
}

.liber_ssh_tunnel_alive <- function(tunnel) {
  inherits(tunnel, "liber_ssh_tunnel") &&
    isTRUE(tryCatch(tunnel$process$is_alive(), error = function(error) FALSE)) &&
    .liber_tcp_ready(tunnel$local_port)
}

.liber_ssh_tunnel_stop <- function(tunnel) {
  if (!inherits(tunnel, "liber_ssh_tunnel")) return(invisible(FALSE))
  alive <- isTRUE(tryCatch(tunnel$process$is_alive(), error = function(error) FALSE))
  if (alive) try(tunnel$process$kill_tree(), silent = TRUE)
  if (nzchar(tunnel$ssh_config_file %||% "") && file.exists(tunnel$ssh_config_file)) {
    unlink(tunnel$ssh_config_file, force = TRUE)
  }
  invisible(alive)
}

.liber_ssh_runtime_allowed <- function(allow_ssh_tunnel = NULL,
                                       session_workspace = FALSE,
                                       host = NULL) {
  .liber_ollama_runtime_allowed(
    allow_ollama = allow_ssh_tunnel,
    session_workspace = session_workspace, host = host
  )
}

.liber_ssh_session_allowed <- function(session, runtime_allowed = FALSE) {
  .liber_ollama_session_allowed(session, runtime_allowed = runtime_allowed)
}
