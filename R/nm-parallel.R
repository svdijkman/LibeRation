#' Parallel-ready subject worker (sequential for now)
#'
#' Maps over subject jobs. When \code{n_threads > 1}, a future version can
#' dispatch to a C++ thread pool or \code{parallel} without changing call sites.
#'
#' @param jobs List of per-subject job specs.
#' @param FUN Function applied to each job.
#' @param n_threads Number of workers (\code{1} = sequential).
#' @param label Optional label for profiling.
#' @keywords internal
.nm_map_subjects <- function(jobs, FUN, n_threads = 1L, label = "subject") {
  n_threads <- as.integer(max(1L, n_threads[1L]))
  n <- length(jobs)
  if (n == 0L) {
    return(list())
  }
  prof_on <- isTRUE(getOption("LibeRation.profile", FALSE))
  t0 <- if (prof_on) proc.time() else NULL
  out <- if (n_threads <= 1L || n <= 1L) {
    lapply(jobs, FUN)
  } else {
    # Reserved for future parallel backends; fall back to sequential.
    lapply(jobs, FUN)
  }
  if (prof_on) {
    .nm_profile_add(label, (proc.time() - t0)[3], n)
  }
  out
}

#' @keywords internal
.nm_resolve_n_threads <- function(control = list()) {
  n <- control$n_threads
  if (is.null(n)) {
    return(1L)
  }
  as.integer(max(1L, n[1L]))
}

#' Chunk subject indices for future parallel Laplace/NLL passes
#' @keywords internal
.nm_subject_chunks <- function(n_sub, n_chunks) {
  n_chunks <- as.integer(max(1L, min(n_chunks, n_sub)))
  split(seq_len(n_sub), cut(seq_len(n_sub), n_chunks, labels = FALSE))
}
