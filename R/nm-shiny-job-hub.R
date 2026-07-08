#' Shared Shiny job push hub (one background worker per app process).
#'
#' Remote HTTP and filesystem watches run in a single \pkg{callr} worker.
#' Connected Shiny sessions receive \code{liberation_job_push} custom messages
#' when job state changes — no per-session remote polling.
#'
#' @keywords internal
.nm_shiny_job_hub <- new.env(parent = emptyenv())

#' @keywords internal
.nm_shiny_job_hub_init <- function() {
  hub <- .nm_shiny_job_hub
  if (isTRUE(hub$initialized)) {
    return(invisible(hub))
  }
  hub$initialized <- TRUE
  hub$sessions <- list()
  hub$rev <- 0L
  hub$event_offset <- 0L
  hub$dir <- file.path(
    tempdir(),
    "LibeRation_shiny_job_hub",
    paste0("pid", Sys.getpid())
  )
  dir.create(hub$dir, recursive = TRUE, showWarnings = FALSE)
  hub$subs_file <- file.path(hub$dir, "subscriptions.json")
  hub$events_file <- file.path(hub$dir, "events.jsonl")
  hub$state_file <- file.path(hub$dir, "state.json")
  hub$poll_file <- file.path(hub$dir, "poll_sec")
  writeLines("{}", hub$subs_file)
  writeLines("", hub$events_file)
  .nm_shiny_job_hub_write_poll(.nm_shiny_job_hub_poll_sec())
  invisible(hub)
}

#' @keywords internal
.nm_shiny_job_hub_poll_sec <- function() {
  val <- getOption("LibeRation.shiny_job_poll_sec", 5)
  sec <- suppressWarnings(as.numeric(val)[1L])
  if (!is.finite(sec) || sec < 1) sec <- 5
  sec
}

#' Update the hub poll interval so the running background worker picks it up.
#'
#' The worker re-reads this each iteration, so changing the interval in the GUI
#' takes effect on the next poll without restarting the worker.
#' @keywords internal
.nm_shiny_job_hub_write_poll <- function(sec) {
  hub <- .nm_shiny_job_hub
  sec <- suppressWarnings(as.numeric(sec)[1L])
  if (!is.finite(sec) || sec < 1) sec <- 5
  options(LibeRation.shiny_job_poll_sec = sec)
  pf <- hub$poll_file
  if (!is.null(pf)) {
    tryCatch(writeLines(as.character(sec), pf), error = function(e) NULL)
  }
  invisible(sec)
}

#' @keywords internal
.nm_shiny_job_hub_write_subscriptions <- function() {
  hub <- .nm_shiny_job_hub
  if (length(hub$sessions) == 0L) {
    writeLines("{}", hub$subs_file)
    return(invisible(0L))
  }
  subs <- lapply(hub$sessions, function(x) {
    list(
      job_root = x$job_root,
      remote_servers = x$remote_servers,
      sync_remote = isTRUE(x$sync_remote)
    )
  })
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    return(invisible(0L))
  }
  jsonlite::write_json(subs, hub$subs_file, auto_unbox = TRUE, pretty = FALSE)
  invisible(length(subs))
}

#' @keywords internal
.nm_shiny_job_hub_subscribe <- function(session, job_root, remote_servers = character(),
                                        sync_remote = TRUE) {
  .nm_shiny_job_hub_init()
  hub <- .nm_shiny_job_hub
  token <- session$token
  hub$sessions[[token]] <- list(
    session = session,
    job_root = job_root,
    remote_servers = unique(as.character(remote_servers)),
    sync_remote = isTRUE(sync_remote)
  )
  .nm_shiny_job_hub_write_subscriptions()
  .nm_shiny_job_hub_start_worker()
  .nm_shiny_job_hub_start_tick()
  invisible(token)
}

#' Update hub subscription for a Shiny session.
#'
#' @param session Shiny session.
#' @param job_root Job directory to watch.
#' @param remote_servers Remote server ids to sync.
#' @param sync_remote If TRUE, refresh remote stubs in the hub worker.
#' @export
nm_shiny_job_hub_update <- function(session, job_root = NULL,
                                    remote_servers = NULL,
                                    sync_remote = NULL) {
  hub <- .nm_shiny_job_hub
  token <- session$token
  if (is.null(hub$sessions[[token]])) {
    .nm_shiny_job_hub_subscribe(
      session,
      job_root %||% nm_job_root(),
      remote_servers = remote_servers %||% character(),
      sync_remote = sync_remote %||% TRUE
    )
    return(invisible(token))
  }
  entry <- hub$sessions[[token]]
  if (!is.null(job_root)) {
    entry$job_root <- job_root
  }
  if (!is.null(remote_servers)) {
    entry$remote_servers <- unique(as.character(remote_servers))
  }
  if (!is.null(sync_remote)) {
    entry$sync_remote <- isTRUE(sync_remote)
  }
  hub$sessions[[token]] <- entry
  .nm_shiny_job_hub_write_subscriptions()
  invisible(token)
}

#' @keywords internal
.nm_shiny_job_hub_unsubscribe <- function(token) {
  hub <- .nm_shiny_job_hub
  hub$sessions[[token]] <- NULL
  .nm_shiny_job_hub_write_subscriptions()
  if (length(hub$sessions) == 0L) {
    .nm_shiny_job_hub_stop_worker()
  }
  invisible(TRUE)
}

#' @keywords internal
.nm_shiny_job_hub_start_worker <- function() {
  hub <- .nm_shiny_job_hub
  if (!is.null(hub$worker) && isTRUE(hub$worker$is_alive())) {
    return(invisible(hub$worker))
  }
  if (!requireNamespace("callr", quietly = TRUE)) {
    return(invisible(NULL))
  }
  # Only treat LIBERATION_PKG_ROOT (set by liberation_shiny in dev mode) as a
  # candidate for pkgload::load_all. We must NOT fall back to the *installed*
  # package directory here: an installed package has no R/*.R sources (only the
  # compiled lazy-load DB), so load_all() errors out, the worker process dies on
  # startup, and no job push events are ever delivered — which is exactly why the
  # Jobs tab stopped auto-updating once the packages were installed rather than
  # loaded from source.
  pkg_root <- Sys.getenv("LIBERATION_PKG_ROOT", "")
  err_log <- file.path(hub$dir %||% tempdir(), "worker.err")
  hub$worker <- callr::r_bg(
    func = function(subs_file, events_file, poll_sec, pkg_root) {
      .is_source_tree <- function(root) {
        nzchar(root) && dir.exists(root) &&
          file.exists(file.path(root, "DESCRIPTION")) &&
          length(list.files(file.path(root, "R"), pattern = "\\.[Rr]$")) > 0L
      }
      loaded <- FALSE
      if (.is_source_tree(pkg_root) && requireNamespace("pkgload", quietly = TRUE)) {
        loaded <- tryCatch({
          pkgload::load_all(pkg_root, quiet = TRUE, compile = FALSE)
          TRUE
        }, error = function(e) FALSE)
      }
      if (!loaded) {
        if (!requireNamespace("LibeRation", quietly = TRUE)) {
          stop("LibeRation package not available in job hub worker.", call. = FALSE)
        }
        suppressPackageStartupMessages(library(LibeRation))
      }
      LibeRation:::.nm_shiny_job_hub_worker_loop(subs_file, events_file, poll_sec)
    },
    args = list(
      subs_file = hub$subs_file,
      events_file = hub$events_file,
      poll_sec = .nm_shiny_job_hub_poll_sec(),
      pkg_root = pkg_root
    ),
    stdout = if (.Platform$OS.type == "windows") "NUL" else "/dev/null",
    stderr = err_log,
    supervise = TRUE
  )
  invisible(hub$worker)
}

#' @keywords internal
.nm_shiny_job_hub_stop_worker <- function() {
  hub <- .nm_shiny_job_hub
  if (!is.null(hub$worker)) {
    tryCatch(hub$worker$kill(), error = function(e) NULL)
    hub$worker <- NULL
  }
  hub$tick_active <- FALSE
  invisible(TRUE)
}

#' @keywords internal
.nm_shiny_job_hub_worker_loop <- function(subs_file, events_file, poll_sec) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite required for job hub worker.", call. = FALSE)
  }
  last <- new.env(parent = emptyenv())
  poll_file <- file.path(dirname(subs_file), "poll_sec")
  read_poll <- function() {
    val <- poll_sec
    if (file.exists(poll_file)) {
      v <- suppressWarnings(as.numeric(readLines(poll_file, warn = FALSE)[1L]))
      if (length(v) == 1L && is.finite(v) && v >= 1) val <- v
    }
    val
  }
  repeat {
    poll_sec <- read_poll()
    if (!file.exists(subs_file)) {
      Sys.sleep(poll_sec)
      next
    }
    subs <- tryCatch(
      jsonlite::read_json(subs_file, simplifyVector = FALSE),
      error = function(e) list()
    )
    if (length(subs) == 0L) {
      Sys.sleep(poll_sec)
      next
    }
    groups <- list()
    for (token in names(subs)) {
      sub <- subs[[token]]
      job_root <- sub$job_root %||% ""
      if (!nzchar(job_root)) {
        next
      }
      # remote_servers arrives from JSON: a single id parses as an atomic scalar,
      # but 0 or 2+ ids parse as a list. Coerce to a plain character vector so
      # downstream sort()/unique() never hit a non-atomic list (which would crash
      # the worker for local-only or multi-server sessions).
      srv <- as.character(sub$remote_servers %||% character())
      sync_remote <- isTRUE(sub$sync_remote)
      gkey <- paste(job_root, paste(sort(srv), collapse = ","), sync_remote, sep = "|")
      if (is.null(groups[[gkey]])) {
        groups[[gkey]] <- list(
          job_root = job_root,
          remote_servers = srv,
          sync_remote = sync_remote,
          tokens = character()
        )
      }
      groups[[gkey]]$tokens <- c(groups[[gkey]]$tokens, token)
    }
    # Keep background HTTP short so a slow/overloaded server can never stall the
    # hub (the burst stress test pegged the server; the old 30s list + 15s/job
    # log timeouts made each poll take minutes, so remote updates appeared to
    # stop). A whole-body tryCatch also stops one bad server from wedging or
    # crashing the loop for every session.
    http_timeout <- as.integer(getOption("LibeRation.shiny_job_hub_timeout", 8L))
    for (gkey in names(groups)) {
      g <- groups[[gkey]]
      sig <- tryCatch({
        if (isTRUE(g$sync_remote) && length(g$remote_servers) > 0L) {
          for (sid in unique(g$remote_servers)) {
            df <- tryCatch(
              nm_remote_job_list(sid, timeout = http_timeout),
              error = function(e) NULL
            )
            if (!is.null(df) && is.null(attr(df, "error"))) {
              .nm_job_sync_remote_stubs(sid, df, g$job_root)
            }
          }
        }
        sig_local <- .nm_job_watch_signature(g$job_root, light = FALSE)
        sig_remote <- .nm_job_remote_log_signature(
          g$job_root, g$remote_servers, timeout = http_timeout
        )
        paste(sig_local, sig_remote, sep = "||")
      }, error = function(e) NULL)
      if (is.null(sig) || identical(last[[gkey]], sig)) {
        next
      }
      last[[gkey]] <- sig
      hub_rev <- as.integer(Sys.time())
      for (token in g$tokens) {
        evt <- list(
          token = token,
          rev = hub_rev,
          sig = sig,
          time = as.character(Sys.time())
        )
        cat(jsonlite::toJSON(evt, auto_unbox = TRUE), "\n",
            file = events_file, append = TRUE)
      }
    }
    Sys.sleep(poll_sec)
  }
}

#' @keywords internal
.nm_shiny_job_hub_read_events <- function() {
  hub <- .nm_shiny_job_hub
  if (!file.exists(hub$events_file)) {
    return(list())
  }
  lines <- readLines(hub$events_file, warn = FALSE)
  if (length(lines) == 0L) {
    return(list())
  }
  offset <- hub$event_offset %||% 0L
  if (offset >= length(lines)) {
    return(list())
  }
  new_lines <- lines[(offset + 1L):length(lines)]
  hub$event_offset <- length(lines)
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    return(list())
  }
  lapply(new_lines, function(ln) {
    if (!nzchar(trimws(ln))) {
      return(NULL)
    }
    tryCatch(jsonlite::fromJSON(ln, simplifyVector = TRUE), error = function(e) NULL)
  })
}

#' @keywords internal
.nm_shiny_job_hub_flush <- function() {
  hub <- .nm_shiny_job_hub
  events <- .nm_shiny_job_hub_read_events()
  events <- events[!vapply(events, is.null, logical(1L))]
  if (length(events) == 0L) {
    return(invisible(0L))
  }
  n <- 0L
  for (evt in events) {
    token <- evt$token
    entry <- hub$sessions[[token]]
    if (is.null(entry) || is.null(entry$session)) {
      next
    }
    hub$rev <- hub$rev + 1L
    n <- n + 1L
    tryCatch(
      entry$session$sendCustomMessage(
        "liberation_job_push",
        list(
          rev = hub$rev,
          sig = evt$sig %||% "",
          time = evt$time %||% ""
        )
      ),
      error = function(e) NULL
    )
  }
  invisible(n)
}

#' @keywords internal
.nm_shiny_job_hub_start_tick <- function() {
  hub <- .nm_shiny_job_hub
  if (isTRUE(hub$tick_active)) {
    return(invisible(FALSE))
  }
  hub$tick_active <- TRUE
  tick <- function() {
    if (length(hub$sessions) == 0L) {
      hub$tick_active <- FALSE
      .nm_shiny_job_hub_stop_worker()
      return(invisible(NULL))
    }
    # Never let a transient error kill the tick chain; if the reschedule below
    # were skipped, all push updates (local and remote) would stop until the
    # GUI restarts. Guard the body and always reschedule.
    tryCatch({
      .nm_shiny_job_hub_write_subscriptions()
      worker_alive <- tryCatch(
        !is.null(hub$worker) && isTRUE(hub$worker$is_alive()),
        error = function(e) FALSE
      )
      if (!is.null(hub$worker) && !worker_alive) {
        hub$worker <- NULL
        .nm_shiny_job_hub_start_worker()
      }
      .nm_shiny_job_hub_flush()
    }, error = function(e) NULL)
    if (requireNamespace("later", quietly = TRUE)) {
      later::later(tick, delay = 0.4)
    } else {
      hub$tick_active <- FALSE
    }
    invisible(NULL)
  }
  if (requireNamespace("later", quietly = TRUE)) {
    later::later(tick, delay = 0.4)
  } else {
    hub$tick_active <- FALSE
  }
  invisible(TRUE)
}

#' Register a Shiny session with the shared job push hub.
#'
#' @param session Shiny session.
#' @param job_root Job directory to watch.
#' @param remote_servers Character vector of remote server ids to sync.
#' @param sync_remote If TRUE, refresh remote stubs in the hub worker.
#' @export
nm_shiny_job_hub_register <- function(session,
                                      job_root = nm_job_root(),
                                      remote_servers = character(),
                                      sync_remote = TRUE) {
  .nm_shiny_job_hub_subscribe(session, job_root, remote_servers, sync_remote)
  session$onSessionEnded(function() {
    .nm_shiny_job_hub_unsubscribe(session$token)
  })
  invisible(session$token)
}

#' Push an immediate job refresh to all connected Shiny sessions (local change).
#'
#' @param job_root Job directory that changed.
#' @export
nm_shiny_job_hub_notify_local <- function(job_root = nm_job_root()) {
  hub <- .nm_shiny_job_hub
  if (length(hub$sessions) == 0L) {
    return(invisible(0L))
  }
  sig <- .nm_job_watch_signature(job_root, light = FALSE)
  hub$rev <- hub$rev + 1L
  n <- 0L
  for (entry in hub$sessions) {
    if (!identical(entry$job_root, job_root)) {
      next
    }
    n <- n + 1L
    tryCatch(
      entry$session$sendCustomMessage(
        "liberation_job_push",
        list(rev = hub$rev, sig = sig, time = as.character(Sys.time()))
      ),
      error = function(e) NULL
    )
  }
  invisible(n)
}
