# Generated from tools/shared/liber-paths.R.
# Run `Rscript tools/sync-shared-runtime.R`; do not edit this copy directly.

.liber_shared_home <- function() {
  home <- path.expand("~")
  if (.Platform$OS.type == "windows") {
    profile <- Sys.getenv("USERPROFILE", unset = home)
    if (nzchar(profile)) home <- profile
  }
  home
}

.liber_shared_user_root <- function(
    ..., envvar = NULL, option = NULL, create = FALSE,
    documents_on_windows = TRUE, normalize = FALSE) {
  configured <- ""
  if (!is.null(envvar) && nzchar(envvar)) {
    configured <- Sys.getenv(envvar, unset = "")
  }
  if (!nzchar(configured) && !is.null(option) && nzchar(option)) {
    candidate <- getOption(option, "")
    if (length(candidate) == 1L && !is.na(candidate)) {
      configured <- as.character(candidate)
    }
  }
  if (nzchar(configured)) {
    path <- path.expand(configured)
  } else {
    home <- .liber_shared_home()
    base <- if (.Platform$OS.type == "windows" &&
                isTRUE(documents_on_windows)) {
      if (tolower(basename(home)) == "documents") home else
        file.path(home, "Documents")
    } else {
      home
    }
    path <- file.path(base, "LibeR", ...)
  }
  if (isTRUE(create) && !dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  if (isTRUE(normalize)) {
    normalizePath(path, winslash = "/", mustWork = FALSE)
  } else {
    path
  }
}
