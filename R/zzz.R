.onLoad <- function(libname, pkgname) {
  if (!requireNamespace("LibeRtAD", quietly = TRUE)) {
    stop("Package 'LibeRtAD' is required.", call. = FALSE)
  }
}

.onUnload <- function(libpath) {
  library.dynam.unload("LibeRation", libpath)
}
