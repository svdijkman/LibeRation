.onLoad <- function(libname, pkgname) {
  if (!requireNamespace("LibeRtAD", quietly = TRUE)) {
    stop("Package 'LibeRtAD' is required.", call. = FALSE)
  }
  if (requireNamespace("shiny", quietly = TRUE)) {
    shiny::registerInputHandler(
      "liberation.codeEditor",
      function(data, shinysession, name) {
        if (is.null(data)) {
          return("")
        }
        as.character(data)
      },
      force = TRUE
    )
  }
}

.onUnload <- function(libpath) {
  library.dynam.unload("LibeRation", libpath)
}
