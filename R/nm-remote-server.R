#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x) || (length(x) == 1L && is.na(x))) y else x
}

#' Config file for remote scheduler servers
#' @keywords internal
.nm_remote_config_path <- function() {
  if (getRversion() >= "4.0") {
    file.path(tools::R_user_dir("LibeRation", which = "config"), "remote_servers.json")
  } else {
    file.path(path.expand("~"), ".LibeRation", "remote_servers.json")
  }
}

#' @keywords internal
.nm_remote_read_config <- function() {
  path <- .nm_remote_config_path()
  if (!file.exists(path)) {
    return(list(servers = list(), default_server = NULL))
  }
  tryCatch(
    jsonlite::read_json(path, simplifyVector = FALSE),
    error = function(e) list(servers = list(), default_server = NULL)
  )
}

#' @keywords internal
.nm_remote_write_config <- function(cfg) {
  path <- .nm_remote_config_path()
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(cfg, path, auto_unbox = TRUE, pretty = TRUE)
  # Best-effort: restrict read access to the file holding remote API tokens.
  if (.Platform$OS.type == "windows") {
    tryCatch(
      system2("icacls", c(path, "/inheritance:r", "/grant:r",
                          paste0(Sys.info()[["user"]], ":(R,W)")),
              stdout = FALSE, stderr = FALSE),
      error = function(e) NULL
    )
  } else {
    tryCatch(Sys.chmod(path, mode = "0600"), error = function(e) NULL)
  }
  invisible(cfg)
}

#' @keywords internal
.nm_remote_require_json <- function() {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop(
      "Package 'jsonlite' is required for remote jobs. ",
      "Install with install.packages('jsonlite').",
      call. = FALSE
    )
  }
}

#' @keywords internal
.nm_remote_require_curl <- function() {
  if (!requireNamespace("curl", quietly = TRUE)) {
    stop(
      "Package 'curl' is required for remote jobs. ",
      "Install with install.packages('curl').",
      call. = FALSE
    )
  }
}

#' List configured remote scheduler servers
#'
#' @return Data frame with id, name, base_url, username (tokens omitted).
#' @examples
#' nm_remote_server_list()
#' @export
nm_remote_server_list <- function() {
  .nm_remote_require_json()
  cfg <- .nm_remote_read_config()
  srv <- cfg$servers
  if (length(srv) == 0L) {
    return(data.frame(
      id = character(),
      name = character(),
      base_url = character(),
      username = character(),
      stringsAsFactors = FALSE
    ))
  }
  ids <- vapply(srv, function(x) as.character(x$id), character(1L))
  data.frame(
    id = ids,
    name = vapply(srv, function(x) as.character(x$name %||% x$id), character(1L)),
    base_url = vapply(srv, function(x) as.character(x$base_url), character(1L)),
    username = vapply(srv, function(x) as.character(x$username %||% ""), character(1L)),
    stringsAsFactors = FALSE
  )
}

#' Add or update a remote scheduler server
#'
#' @param name Display name.
#' @param base_url API base URL (e.g. \code{"http://cluster:8080"}).
#' @param username Scheduler username.
#' @param token API token.
#' @param id Optional server id; auto-generated if omitted.
#' @param set_default If TRUE, make this the default remote target.
#' @param insecure If TRUE, disable TLS certificate verification (development
#'   only, e.g. for a self-signed certificate). Leave FALSE in production.
#' @return Server id (invisibly).
#' @examples
#' \dontrun{
#' nm_remote_server_add("HPC", "https://hpc", "alice", token = "lr_...")
#' }
#' @export
nm_remote_server_add <- function(name,
                                 base_url,
                                 username,
                                 token,
                                 id = NULL,
                                 set_default = FALSE,
                                 insecure = FALSE) {
  .nm_remote_require_json()
  cfg <- .nm_remote_read_config()
  base_url <- sub("/+$", "", as.character(base_url))
  if (is.null(id) || !nzchar(id)) {
    id <- paste0("srv_", substr(digest::digest(paste(name, base_url, Sys.time())), 1L, 8L))
  }
  if (!is.null(cfg$servers[[id]]) &&
      (is.null(token) || length(token) == 0L || !nzchar(as.character(token)))) {
    token <- cfg$servers[[id]]$token
  }
  if (is.null(token) || length(token) == 0L || !nzchar(as.character(token))) {
    stop("API token is required.", call. = FALSE)
  }
  token <- trimws(as.character(token))
  entry <- list(
    id = id,
    name = name,
    base_url = base_url,
    username = username,
    token = token,
    insecure = isTRUE(insecure)
  )
  cfg$servers[[id]] <- entry
  if (isTRUE(set_default) || is.null(cfg$default_server)) {
    cfg$default_server <- id
  }
  .nm_remote_write_config(cfg)
  invisible(id)
}

#' Remove a remote server
#'
#' @param id Server id from \code{\link{nm_remote_server_list}}.
#' @export
nm_remote_server_remove <- function(id) {
  .nm_remote_require_json()
  cfg <- .nm_remote_read_config()
  cfg$servers[[id]] <- NULL
  if (identical(cfg$default_server, id)) {
    ids <- names(cfg$servers)
    cfg$default_server <- if (length(ids) > 0L) ids[[1L]] else NULL
  }
  .nm_remote_write_config(cfg)
  invisible(TRUE)
}

#' Set default remote server
#'
#' @param id Server id.
#' @export
nm_remote_server_set_default <- function(id) {
  cfg <- .nm_remote_read_config()
  if (is.null(cfg$servers[[id]])) {
    stop("Unknown server id: ", id, call. = FALSE)
  }
  cfg$default_server <- id
  .nm_remote_write_config(cfg)
  invisible(id)
}

#' @keywords internal
.nm_remote_get_server <- function(server = NULL) {
  cfg <- .nm_remote_read_config()
  sid <- server %||% cfg$default_server
  if (is.null(sid) || !nzchar(sid)) {
    stop("No remote server configured.", call. = FALSE)
  }
  entry <- cfg$servers[[sid]]
  if (is.null(entry)) {
    stop("Unknown remote server: ", sid, call. = FALSE)
  }
  entry
}

#' Test connection to a remote scheduler
#'
#' Checks API health and validates the stored API token.
#'
#' @param server Server id or NULL for default.
#' @return List with health response, authenticated username, and optional version warnings.
#' @export
nm_remote_server_test <- function(server = NULL) {
  health <- .nm_remote_http("GET", "/v1/health", server = server)
  auth <- tryCatch(
    .nm_remote_http("GET", "/v1/auth", server = server),
    error = function(e) {
      stop(
        "Health check OK but API token was rejected. ",
        "Re-issue a token in LibeRties admin and update the client server entry. ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  c(health, list(username = auth$username))
}

#' @keywords internal
.nm_remote_http <- function(method,
                            path,
                            server = NULL,
                            body = NULL,
                            timeout = 120L) {
  .nm_remote_require_json()
  .nm_remote_require_curl()
  entry <- .nm_remote_get_server(server)
  url <- paste0(entry$base_url, path)
  base_url <- as.character(entry$base_url %||% "")
  is_https <- grepl("^https://", base_url, ignore.case = TRUE)
  is_local <- grepl("^https?://(localhost|127\\.0\\.0\\.1|\\[::1\\])", base_url, ignore.case = TRUE)
  # Loopback plaintext is exempt by default (never hits the network). Set
  # options(LibeRation.warn_plain_http_local = TRUE) to also warn for localhost,
  # e.g. when rehearsing a production setup on a single machine.
  warn_local <- isTRUE(getOption("LibeRation.warn_plain_http_local", FALSE))
  if (!is_https && (!is_local || warn_local) &&
      !isTRUE(getOption("LibeRation.warned_plain_http", FALSE))) {
    warning(
      "Remote server '", entry$name %||% entry$id %||% base_url,
      "' uses plain http://. GDPR-sensitive data should travel over https:// ",
      "(configure a TLS reverse proxy on the server). ",
      call. = FALSE
    )
    options(LibeRation.warned_plain_http = TRUE)
  }
  headers <- c(
    "X-API-Token" = entry$token,
    "Authorization" = paste("Bearer", entry$token),
    "Content-Type" = "application/json",
    "Accept" = "application/json"
  )
  if (requireNamespace("LibeRation", quietly = TRUE)) {
    headers["X-LibeRation-Version"] <- as.character(utils::packageVersion("LibeRation"))
  }
  if (requireNamespace("LibeRtAD", quietly = TRUE)) {
    headers["X-LibeRtAD-Version"] <- as.character(utils::packageVersion("LibeRtAD"))
  }
  h <- curl::new_handle()
  curl::handle_setheaders(h, .list = headers)
  # Keep TLS verification ON. An explicit per-server `insecure = TRUE` allows a
  # self-signed cert for dev only, and is loudly discouraged.
  insecure <- isTRUE(entry$insecure)
  if (insecure && !isTRUE(getOption("LibeRation.warned_insecure_tls", FALSE))) {
    warning(
      "Remote server '", entry$name %||% entry$id,
      "' has insecure = TRUE: TLS certificate verification is DISABLED. ",
      "Use only for development with self-signed certificates.",
      call. = FALSE
    )
    options(LibeRation.warned_insecure_tls = TRUE)
  }
  curl::handle_setopt(
    h,
    timeout = timeout,
    post = FALSE,
    postfields = NULL,
    customrequest = NULL,
    ssl_verifypeer = if (insecure) 0L else 1L,
    ssl_verifyhost = if (insecure) 0L else 2L
  )
  if (!is.null(body)) {
    body_raw <- jsonlite::toJSON(body, auto_unbox = TRUE, null = "null")
    curl::handle_setopt(h, postfields = body_raw)
  }
  m <- toupper(method)
  txt <- tryCatch({
    if (m == "GET") {
      curl::handle_setopt(h, httpget = TRUE, post = FALSE, customrequest = "GET")
      curl::curl_fetch_memory(url, handle = h)
    } else if (m == "POST") {
      curl::handle_setopt(h, customrequest = "POST")
      curl::curl_fetch_memory(url, handle = h)
    } else if (m == "DELETE") {
      curl::handle_setopt(h, customrequest = "DELETE")
      curl::curl_fetch_memory(url, handle = h)
    } else {
      stop("Unsupported HTTP method: ", method, call. = FALSE)
    }
  }, error = function(e) {
    stop("Remote server unreachable (", entry$base_url, "): ",
         conditionMessage(e), call. = FALSE)
  })
  status <- txt$status_code
  raw_txt <- rawToChar(txt$content)
  parsed <- tryCatch(
    jsonlite::fromJSON(raw_txt, simplifyVector = TRUE),
    error = function(e) list(raw = raw_txt)
  )
  if (status >= 400L) {
    msg <- if (!is.null(parsed$error)) {
      as.character(parsed$error)
    } else if (is.character(parsed) && length(parsed) == 1L) {
      parsed
    } else {
      raw_txt
    }
    msg <- sub("^[[:space:]]+", "", sub("[[:space:]]+$", "", msg))
    stop("Remote API error (", status, "): ", msg, call. = FALSE)
  }
  parsed
}

#' @keywords internal
.nm_data_md5 <- function(data) {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(data, tmp)
  as.character(tools::md5sum(tmp))
}

#' @keywords internal
.nm_rds_b64 <- function(x) {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp)
  jsonlite::base64_enc(readBin(tmp, raw(), file.info(tmp)$size))
}

#' @keywords internal
.nm_rds_read_raw <- function(raw) {
  is_gzip <- length(raw) >= 2L &&
    identical(raw[1:2], as.raw(c(0x1f, 0x8b)))
  con <- if (is_gzip) {
    gzcon(rawConnection(raw, "r"))
  } else {
    rawConnection(raw, "r")
  }
  on.exit(close(con), add = TRUE)
  readRDS(con)
}

#' @keywords internal
.nm_b64_rds <- function(b64) {
  raw <- jsonlite::base64_dec(b64)
  .nm_rds_read_raw(raw)
}
