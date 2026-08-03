.nm_psd_projection <- function(matrix, floor = 0) {
  decomposition <- eigen((matrix + t(matrix)) / 2, symmetric = TRUE)
  values <- pmax(decomposition$values, floor)
  projected <- decomposition$vectors %*% (values * t(decomposition$vectors))
  (projected + t(projected)) / 2
}

.nm_modified_cholesky <- function(matrix, tolerance) {
  n <- nrow(matrix)
  if (!n) return(list(matrix = matrix, permutation = integer(), d = numeric()))
  scale <- max(abs(matrix), 1)
  delta <- tolerance * scale
  off_diagonal <- if (n > 1L) matrix[row(matrix) != col(matrix)] else 0
  beta_squared <- max(
    max(abs(diag(matrix))),
    if (n > 1L) max(abs(off_diagonal)) / sqrt(max(1, n^2 - 1)) else 0,
    delta
  )
  work <- matrix
  lower <- diag(n)
  diagonal <- numeric(n)
  permutation <- seq_len(n)
  for (column in seq_len(n)) {
    candidates <- column:n
    previous <- if (column > 1L) seq_len(column - 1L) else integer()
    reduced_diagonal <- vapply(candidates, function(index) {
      work[index, index] - if (length(previous)) {
        sum(lower[index, previous]^2 * diagonal[previous])
      } else 0
    }, numeric(1))
    pivot <- candidates[[which.max(abs(reduced_diagonal))]]
    if (pivot != column) {
      work[c(column, pivot), ] <- work[c(pivot, column), ]
      work[, c(column, pivot)] <- work[, c(pivot, column)]
      permutation[c(column, pivot)] <- permutation[c(pivot, column)]
      if (length(previous)) {
        lower[c(column, pivot), previous] <- lower[c(pivot, column), previous, drop = FALSE]
      }
    }
    c_column <- numeric(max(0L, n - column))
    if (column < n) {
      rows <- (column + 1L):n
      c_column <- work[rows, column]
      if (length(previous)) {
        c_column <- c_column - as.vector(
          lower[rows, previous, drop = FALSE] %*%
            (diagonal[previous] * lower[column, previous])
        )
      }
    }
    c_diagonal <- work[column, column] - if (length(previous)) {
      sum(lower[column, previous]^2 * diagonal[previous])
    } else 0
    theta <- if (length(c_column)) max(abs(c_column)) else 0
    diagonal[[column]] <- max(abs(c_diagonal), theta^2 / beta_squared, delta)
    if (column < n) lower[(column + 1L):n, column] <- c_column / diagonal[[column]]
  }
  permuted <- lower %*% diag(diagonal, nrow = n) %*% t(lower)
  repaired <- base::matrix(0, n, n)
  repaired[permutation, permutation] <- permuted
  repaired <- (repaired + t(repaired)) / 2
  list(matrix = repaired, permutation = permutation, d = diagonal)
}

#' Explicit covariance-matrix repair
#'
#' Repairs are never selected automatically by this function. `"clip"` is
#' appropriate only for negative eigenvalues consistent with floating-point
#' round-off. `"jitter"` is an explicit diagonal-loading repair.
#' `"higham"` uses alternating projections to obtain a nearest positive
#' semidefinite matrix and can preserve a strictly positive original diagonal.
#' `"modified_cholesky"` uses a symmetric pivoted modified-Cholesky
#' factorisation to obtain a well-conditioned positive-definite matrix. It is
#' useful when a downstream solve requires positive definiteness, but is not a
#' nearest-matrix method.
#' Material repairs change the statistical model and should be sensitivity
#' analysed rather than treated as an invisible numerical fix.
#'
#' @param matrix Finite square covariance-like matrix.
#' @param method One of `"none"`, `"clip"`, `"jitter"`, `"higham"`, or
#'   `"modified_cholesky"`.
#' @param tolerance Relative eigenvalue/convergence tolerance.
#' @param preserve_diagonal Preserve the original diagonal for Higham repair.
#' @param max_iterations Maximum Higham alternating-projection iterations.
#' @return A list containing `matrix` and an auditable `diagnostics` record.
#' @export
nm_covariance_repair <- function(matrix,
                                 method = c("none", "clip", "jitter", "higham",
                                            "modified_cholesky"),
                                 tolerance = 1e-10,
                                 preserve_diagonal = FALSE,
                                 max_iterations = 100L) {
  method <- match.arg(method)
  matrix <- as.matrix(matrix)
  if (!is.numeric(matrix) || nrow(matrix) != ncol(matrix) || any(!is.finite(matrix))) {
    .nm_stop("`matrix` must be a finite numeric square matrix.")
  }
  tolerance <- as.numeric(tolerance)
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0) {
    .nm_stop("`tolerance` must be one finite positive number.")
  }
  max_iterations <- as.integer(max_iterations)
  if (length(max_iterations) != 1L || is.na(max_iterations) || max_iterations < 1L) {
    .nm_stop("`max_iterations` must be a positive integer.")
  }
  asymmetry <- max(abs(matrix - t(matrix)), 0)
  original <- (matrix + t(matrix)) / 2
  if (!length(original)) return(list(
    matrix = original,
    diagnostics = list(
      method = method, original_eigenvalues = numeric(),
      repaired_eigenvalues = numeric(), adjustment_norm = 0,
      relative_adjustment = 0, rank = 0L, threshold = tolerance,
      diagonal_shift = 0, iterations = 0L, converged = TRUE,
      preserve_diagonal = isTRUE(preserve_diagonal), asymmetry = asymmetry,
      material_indefiniteness = FALSE, permutation = integer(),
      correction_diagonal = numeric(), condition_number = 1
    )
  ))
  original_eigenvalues <- eigen(original, symmetric = TRUE, only.values = TRUE)$values
  scale <- max(abs(original_eigenvalues), max(abs(original)), 1)
  threshold <- tolerance * scale
  material <- min(original_eigenvalues) < -threshold
  iterations <- 0L
  converged <- TRUE
  diagonal_shift <- 0
  permutation <- seq_len(nrow(original))
  repaired <- original
  if (method == "none") {
    if (min(original_eigenvalues) < -threshold) {
      .nm_stop(
        "Matrix is materially indefinite (minimum eigenvalue ",
        format(min(original_eigenvalues), digits = 7),
        "). Select and audit an explicit repair method."
      )
    }
  } else if (method == "clip") {
    if (material) {
      warning(
        "Eigenvalue clipping is repairing material indefiniteness; treat the result as a changed covariance model and run sensitivity analyses.",
        call. = FALSE
      )
    }
    repaired <- .nm_psd_projection(original)
  } else if (method == "jitter") {
    diagonal_shift <- max(0, threshold - min(original_eigenvalues))
    repaired <- original + diag(diagonal_shift, nrow(original))
  } else if (method == "higham") {
    target_diagonal <- diag(original)
    if (isTRUE(preserve_diagonal) && any(target_diagonal <= 0)) {
      .nm_stop("Higham repair can preserve only a strictly positive covariance diagonal.")
    }
    y <- original
    correction <- base::matrix(0, nrow(original), ncol(original))
    converged <- FALSE
    for (iteration in seq_len(max_iterations)) {
      residual <- y - correction
      projected <- .nm_psd_projection(residual)
      correction <- projected - residual
      next_y <- projected
      if (isTRUE(preserve_diagonal)) diag(next_y) <- target_diagonal
      change <- sqrt(sum((next_y - y)^2)) / max(1, sqrt(sum(y^2)))
      y <- (next_y + t(next_y)) / 2
      iterations <- iteration
      if (change <= tolerance) {
        converged <- TRUE
        break
      }
    }
    if (!converged) .nm_stop("Higham covariance repair did not converge.")
    repaired <- .nm_psd_projection(y)
    if (isTRUE(preserve_diagonal)) {
      current <- diag(repaired)
      if (any(current <= 0)) .nm_stop("Higham covariance repair produced a zero diagonal.")
      rescale <- sqrt(target_diagonal / current)
      repaired <- repaired * outer(rescale, rescale)
      repaired <- (repaired + t(repaired)) / 2
    }
  } else {
    if (isTRUE(preserve_diagonal)) {
      .nm_stop("Modified-Cholesky repair cannot preserve the original diagonal.")
    }
    modified <- .nm_modified_cholesky(original, tolerance = tolerance)
    repaired <- modified$matrix
    permutation <- modified$permutation
  }
  repaired_eigenvalues <- eigen(repaired, symmetric = TRUE, only.values = TRUE)$values
  norm_original <- sqrt(sum(original^2))
  adjustment_norm <- sqrt(sum((repaired - original)^2))
  list(
    matrix = repaired,
    diagnostics = list(
      method = method, original_eigenvalues = original_eigenvalues,
      repaired_eigenvalues = repaired_eigenvalues,
      adjustment_norm = adjustment_norm,
      relative_adjustment = adjustment_norm / max(norm_original, .Machine$double.eps),
      rank = sum(repaired_eigenvalues > threshold), threshold = threshold,
      diagonal_shift = diagonal_shift, iterations = iterations,
      converged = converged, preserve_diagonal = isTRUE(preserve_diagonal),
      asymmetry = asymmetry, material_indefiniteness = material,
      permutation = permutation,
      correction_diagonal = diag(repaired - original),
      condition_number = if (min(repaired_eigenvalues) > 0) {
        max(repaired_eigenvalues) / min(repaired_eigenvalues)
      } else Inf
    )
  )
}

#' Compare downstream results across explicit covariance repairs
#'
#' This helper makes the statistical consequence of a repair visible. The
#' caller supplies the downstream calculation (for example a design criterion,
#' target-attainment probability, or uncertainty summary), and the same
#' calculation is evaluated for each requested repair. Failed alternatives are
#' retained in the result rather than silently discarded.
#'
#' @param matrix Finite square covariance-like matrix.
#' @param statistic Function accepting one covariance matrix and returning a
#'   finite numeric scalar or named vector.
#' @param methods Repair methods to compare.
#' @param reference_method Method against which differences are reported. The
#'   first successful method is used when `NULL`.
#' @param ... Additional arguments passed to [nm_covariance_repair()].
#' @return A `nm_covariance_sensitivity` list containing repaired matrices,
#'   repair diagnostics, and a long-form comparison table.
#' @export
nm_covariance_repair_sensitivity <- function(
    matrix, statistic,
    methods = c("jitter", "higham", "modified_cholesky"),
    reference_method = NULL, ...) {
  if (!is.function(statistic)) .nm_stop("`statistic` must be a function.")
  methods <- unique(as.character(methods))
  allowed <- c("none", "clip", "jitter", "higham", "modified_cholesky")
  if (!length(methods) || any(!methods %in% allowed)) {
    .nm_stop("`methods` must contain supported covariance-repair methods.")
  }
  evaluations <- lapply(methods, function(method) {
    warnings <- character()
    tryCatch(
      withCallingHandlers({
        repair <- nm_covariance_repair(matrix, method = method, ...)
        value <- statistic(repair$matrix)
        if (!is.numeric(value) || !length(value) || any(!is.finite(value))) {
          .nm_stop("`statistic` must return a finite numeric scalar or vector.")
        }
        names(value) <- names(value) %||% paste0("statistic_", seq_along(value))
        list(ok = TRUE, repair = repair, value = value, warnings = warnings)
      }, warning = function(condition) {
        warnings <<- c(warnings, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }),
      error = function(condition) list(
        ok = FALSE, repair = NULL, value = numeric(), warnings = warnings,
        error = conditionMessage(condition)
      )
    )
  })
  names(evaluations) <- methods
  successful <- methods[vapply(evaluations, `[[`, logical(1), "ok")]
  if (!length(successful)) .nm_stop("No covariance-repair sensitivity scenario succeeded.")
  reference_method <- reference_method %||% successful[[1L]]
  if (!reference_method %in% successful) {
    .nm_stop("`reference_method` must identify a successful sensitivity scenario.")
  }
  reference <- evaluations[[reference_method]]$value
  rows <- lapply(methods, function(method) {
    evaluation <- evaluations[[method]]
    if (!evaluation$ok) return(data.frame(
      method = method, statistic = NA_character_, estimate = NA_real_,
      difference = NA_real_, relative_difference = NA_real_,
      relative_matrix_adjustment = NA_real_, status = "failed",
      message = evaluation$error %||% "failed", stringsAsFactors = FALSE
    ))
    value <- evaluation$value
    common <- intersect(names(value), names(reference))
    data.frame(
      method = method, statistic = names(value), estimate = as.numeric(value),
      difference = vapply(names(value), function(name) {
        if (name %in% common) value[[name]] - reference[[name]] else NA_real_
      }, numeric(1)),
      relative_difference = vapply(names(value), function(name) {
        if (!name %in% common) return(NA_real_)
        (value[[name]] - reference[[name]]) /
          max(abs(reference[[name]]), .Machine$double.eps)
      }, numeric(1)),
      relative_matrix_adjustment = evaluation$repair$diagnostics$relative_adjustment,
      status = "ok", message = paste(evaluation$warnings, collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
  structure(list(
    reference_method = reference_method, evaluations = evaluations,
    comparison = do.call(rbind, rows)
  ), class = "nm_covariance_sensitivity")
}
