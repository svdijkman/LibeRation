#!/usr/bin/env Rscript
# Remove legacy NONMEM benchmark artifacts from inst/ (one-time cleanup).
# Safe to run: keeps only inst/nonmem_bench/README.md.

pkg_root <- if (length(grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))) {
  normalizePath(
    file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]])), ".."),
    winslash = "/",
    mustWork = FALSE
  )
} else {
  normalizePath(".", winslash = "/", mustWork = FALSE)
}

target <- file.path(pkg_root, "inst", "nonmem_bench")
if (!dir.exists(target)) {
  message("Nothing to clean: ", target)
  quit(save = "no", status = 0)
}

keep <- file.path(target, "README.md")
removed <- 0L
for (nm in list.dirs(target, full.names = FALSE, recursive = FALSE)) {
  if (!nzchar(nm)) next
  unlink(file.path(target, nm), recursive = TRUE)
  removed <- removed + 1L
}
for (f in list.files(target, full.names = TRUE)) {
  if (!identical(normalizePath(f, winslash = "/"), normalizePath(keep, winslash = "/"))) {
    unlink(f)
    removed <- removed + 1L
  }
}
if (!file.exists(keep)) {
  writeLines(
    c(
      "# Moved",
      "",
      "See bench/nonmem_bench/ at the package root.",
      "Delete this folder if it reappears with benchmark run outputs."
    ),
    keep
  )
}
message("Cleaned ", removed, " item(s) under ", target)
