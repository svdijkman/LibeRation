#' @keywords internal
.nm_shiny_dev_roots <- function() {
  nm_root <- Sys.getenv("LibeRation_ROOT", "")
  ad_root <- Sys.getenv("LibeRation_LIBERTAD_ROOT", "")
  if (!nzchar(nm_root)) {
    shiny_dir <- system.file("shiny", package = "LibeRation")
    if (nzchar(shiny_dir)) {
      cand <- normalizePath(file.path(shiny_dir, "..", ".."), winslash = "/", mustWork = FALSE)
      desc_path <- file.path(cand, "DESCRIPTION")
      if (file.exists(desc_path)) {
        desc <- tryCatch(read.dcf(desc_path), error = function(e) NULL)
        if (!is.null(desc) && desc[1L, "Package"] == "LibeRation") {
          nm_root <- cand
        }
      }
    }
  }
  if (nzchar(nm_root)) {
    nm_root <- normalizePath(nm_root, winslash = "/", mustWork = FALSE)
  }
  if (!nzchar(ad_root) && nzchar(nm_root)) {
    for (cand in c(
      file.path(nm_root, "..", "LibeRtAD"),
      file.path(dirname(nm_root), "LibeRtAD")
    )) {
      if (dir.exists(cand)) {
        ad_root <- normalizePath(cand, winslash = "/", mustWork = FALSE)
        break
      }
    }
  }
  list(nm_root = nm_root, ad_root = ad_root)
}

#' @keywords internal
.nm_shiny_load_dev <- function() {
  roots <- .nm_shiny_dev_roots()
  if (!nzchar(roots$nm_root)) {
    return(invisible(FALSE))
  }
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    message("[LibeRation] pkgload not available; using installed LibeRation.")
    return(invisible(FALSE))
  }
  if (nzchar(roots$ad_root)) {
    pkgload::load_all(roots$ad_root, quiet = TRUE, compile = FALSE)
  }
  pkgload::load_all(roots$nm_root, quiet = TRUE, recompile = FALSE)
  env <- list(
    mode = "dev",
    nm_root = roots$nm_root,
    ad_root = roots$ad_root
  )
  options(LibeRation.job_dev_env = env)
  message("[LibeRation] Loaded dev sources from ", roots$nm_root)
  invisible(TRUE)
}

#' Launch the LibeRation Shiny GUI (async job queue)
#'
#' Opens a browser-based interface for dataset preview, model setup,
#' background estimation via \code{\link{nm_job_submit}}, and GOF plots.
#'
#' Requires Suggested packages \pkg{shiny}, \pkg{callr}, and \pkg{ggplot2}.
#'
#' @param host Passed to \code{\link[shiny]{runApp}}.
#' @param port Passed to \code{\link[shiny]{runApp}}; \code{NULL} picks a random port.
#' @param launch.browser Open the app in a browser when \code{TRUE}.
#' @param poll_ms Job table auto-refresh interval in milliseconds (default \code{2000}).
#' @param workspace Workspace root directory (Pirana-style project tree). Initialized
#'   via \code{\link{nm_workspace_init}} when \code{NULL} uses the configured or default root.
#' @param create_demo_project When initializing workspace, create a THEO demo project if empty.
#' @return Invisibly, the Shiny app object.
#' @examples
#' \dontrun{
#' liberation_shiny(port = 8765L)
#' }
#' @export
liberation_shiny <- function(host = "127.0.0.1",
                     port = NULL,
                     launch.browser = TRUE,
                     poll_ms = 2000L,
                     workspace = NULL,
                     create_demo_project = TRUE) {
  if (!requireNamespace("DT", quietly = TRUE)) {
    stop("Package 'DT' is required. Install with install.packages('DT').")
  }
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required. Install with install.packages('shiny').")
  }
  dev <- .nm_shiny_load_dev()
  options(LibeRation.shiny_poll_ms = as.integer(poll_ms))
  if (is.null(getOption("LibeRation.job_dev_env"))) {
    options(LibeRation.job_dev_env = .nm_job_dev_env())
  }
  ws <- if (is.null(workspace)) nm_workspace_root() else workspace
  nm_workspace_init(ws, create_demo_project = create_demo_project)
  message("[LibeRation] Shiny GUI — workspace: ", nm_workspace_root())
  message("[LibeRation] Shiny GUI — job directory: ", nm_job_root())
  app_dir <- system.file("shiny", package = "LibeRation")
  if (isTRUE(dev)) {
    dev_app <- file.path(.nm_shiny_dev_roots()$nm_root, "inst", "shiny")
    if (dir.exists(dev_app)) {
      app_dir <- dev_app
    }
    Sys.setenv(LIBERATION_PKG_ROOT = .nm_shiny_dev_roots()$nm_root)
  }
  if (!nzchar(app_dir) || !dir.exists(app_dir)) {
    stop("Shiny app not found in LibeRation (inst/shiny). Reinstall or use devtools::load_all().")
  }
  message("[LibeRation] Shiny GUI — job updates via push hub (poll interval configurable on Jobs tab)")
  shiny::runApp(
    app_dir,
    host = host,
    port = port,
    launch.browser = launch.browser
  )
}
