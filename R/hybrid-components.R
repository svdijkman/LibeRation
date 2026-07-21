.nm_component_hash <- function(value) {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(value, path, version = 3)
  unname(tools::md5sum(path))
}

#' Define an offline differentiable model component
#'
#' Components are expanded to restricted scalar expression IR and therefore
#' remain inside C++/CppAD during simulation and estimation. Weight/data payloads
#' are immutable, hashed, serializable, and never access a network.
#'
#' @param name Stable component label.
#' @param type `"dense_nn"`, `"linear_spline"`, or `"gaussian_process"`.
#' @param scope Compile into `$PK/$PRED`, or directly into `$DES` for a
#'   state/time-dependent learned dynamics or correction term.
#' @param inputs Model symbols used as inputs.
#' @param outputs Generated model assignment names.
#' @param weights,biases Dense-network layer matrices and bias vectors.
#' @param activation Hidden-layer activation.
#' @param knots,values Linear-spline knots and values.
#' @param training,alpha,lengthscale,variance,mean Gaussian-process prediction
#'   payload.
#' @export
nm_component <- function(name,
                         type = c("dense_nn", "linear_spline", "gaussian_process"),
                         scope = c("pred", "des"),
                         inputs, outputs, weights = NULL, biases = NULL,
                         activation = c("tanh", "softplus", "relu"),
                         knots = NULL, values = NULL, training = NULL,
                         alpha = NULL, lengthscale = 1, variance = 1, mean = 0) {
  name <- trimws(as.character(name)); type <- match.arg(type); scope <- match.arg(scope)
  inputs <- trimws(as.character(inputs)); outputs <- trimws(as.character(outputs))
  identifiers <- c(inputs, outputs)
  if (length(name) != 1L || is.na(name) || !nzchar(name) || !length(inputs) ||
      !length(outputs) || anyNA(identifiers) || any(!grepl("^[A-Za-z][A-Za-z0-9_]*$", identifiers))) {
    .nm_stop("A component requires a name and valid input/output symbols.")
  }
  payload <- list()
  if (type == "dense_nn") {
    weights <- lapply(weights, as.matrix); biases <- lapply(biases, as.numeric)
    if (!length(weights) || length(weights) != length(biases)) {
      .nm_stop("A dense component requires matching weight matrices and bias vectors.")
    }
    previous <- length(inputs)
    for (layer in seq_along(weights)) {
      if (ncol(weights[[layer]]) != previous || nrow(weights[[layer]]) != length(biases[[layer]]) ||
          any(!is.finite(weights[[layer]])) || any(!is.finite(biases[[layer]]))) {
        .nm_stop("Dense component layer dimensions or values are invalid.")
      }
      previous <- nrow(weights[[layer]])
    }
    if (previous != length(outputs)) .nm_stop("Final dense layer must have one unit per output.")
    payload <- list(weights = weights, biases = biases, activation = match.arg(activation))
  } else if (type == "linear_spline") {
    if (length(inputs) != 1L || length(outputs) != 1L) {
      .nm_stop("A linear spline has exactly one input and output.")
    }
    knots <- as.numeric(knots); values <- as.numeric(values)
    if (length(knots) < 2L || length(values) != length(knots) ||
        any(!is.finite(c(knots, values))) || is.unsorted(knots, strictly = TRUE)) {
      .nm_stop("Spline knots must be strictly increasing with one finite value per knot.")
    }
    payload <- list(knots = knots, values = values)
  } else {
    training <- as.matrix(training); alpha <- as.numeric(alpha)
    lengthscale <- as.numeric(lengthscale); variance <- as.numeric(variance); mean <- as.numeric(mean)
    if (ncol(training) != length(inputs) || nrow(training) != length(alpha) ||
        any(!is.finite(c(training, alpha))) || !length(alpha)) {
      .nm_stop("Gaussian-process training rows, alpha, and inputs are inconsistent.")
    }
    if (length(lengthscale) == 1L) lengthscale <- rep(lengthscale, length(inputs))
    if (length(lengthscale) != length(inputs) || any(!is.finite(lengthscale)) ||
        any(lengthscale <= 0) || length(outputs) != 1L ||
        length(variance) != 1L || !is.finite(variance) || variance <= 0 ||
        length(mean) != 1L || !is.finite(mean)) {
      .nm_stop("Gaussian-process kernel controls are invalid.")
    }
    payload <- list(training = training, alpha = alpha, lengthscale = lengthscale,
                    variance = variance, mean = mean)
  }
  value <- list(version = 1L, name = name, type = type, scope = scope, inputs = inputs,
                outputs = outputs, payload = payload, offline = TRUE,
                differentiable = TRUE)
  value$hash <- .nm_component_hash(value)
  structure(value, class = "nm_component")
}

.nm_number_code <- function(value) format(as.numeric(value), digits = 17L, scientific = TRUE)

#' Generate restricted expression code for a model component
#' @param component An [nm_component()].
#' @export
nm_component_code <- function(component) {
  if (!inherits(component, "nm_component")) .nm_stop("`component` must be an nm_component.")
  prefix <- paste0("CMP_", gsub("[^A-Za-z0-9_]", "_", toupper(component$name)))
  if (component$type == "dense_nn") {
    current <- component$inputs; lines <- character()
    activation <- function(expression, final) {
      if (final) return(expression)
      switch(component$payload$activation,
             tanh = paste0("tanh(", expression, ")"),
             softplus = paste0("log(1 + exp(", expression, "))"),
             relu = paste0("pmax(", expression, ", 0)"))
    }
    for (layer in seq_along(component$payload$weights)) {
      matrix <- component$payload$weights[[layer]]; bias <- component$payload$biases[[layer]]
      next_names <- if (layer == length(component$payload$weights)) component$outputs else
        paste0(prefix, "_L", layer, "_", seq_len(nrow(matrix)))
      for (unit in seq_len(nrow(matrix))) {
        terms <- c(.nm_number_code(bias[[unit]]), vapply(seq_along(current), function(index) {
          paste0(.nm_number_code(matrix[unit, index]), " * ", current[[index]])
        }, character(1)))
        lines <- c(lines, paste0(next_names[[unit]], " = ",
                                 activation(paste(terms, collapse = " + "),
                                            layer == length(component$payload$weights))))
      }
      current <- next_names
    }
    return(paste(lines, collapse = "\n"))
  }
  if (component$type == "linear_spline") {
    x <- component$inputs[[1L]]; knots <- component$payload$knots
    values <- component$payload$values
    expression <- .nm_number_code(values[[length(values)]])
    for (index in rev(seq_len(length(knots) - 1L))) {
      segment <- paste0(.nm_number_code(values[[index]]), " + (",
                        .nm_number_code(values[[index + 1L]] - values[[index]]), ") * (",
                        x, " - ", .nm_number_code(knots[[index]]), ") / ",
                        .nm_number_code(knots[[index + 1L]] - knots[[index]]))
      expression <- paste0("ifelse(", x, " <= ", .nm_number_code(knots[[index + 1L]]),
                           ", ", segment, ", ", expression, ")")
    }
    expression <- paste0("ifelse(", x, " <= ", .nm_number_code(knots[[1L]]), ", ",
                         .nm_number_code(values[[1L]]), ", ", expression, ")")
    return(paste0(component$outputs[[1L]], " = ", expression))
  }
  payload <- component$payload
  kernels <- vapply(seq_len(nrow(payload$training)), function(row) {
    distance <- vapply(seq_along(component$inputs), function(column) {
      paste0("((", component$inputs[[column]], " - ",
             .nm_number_code(payload$training[row, column]), ") / ",
             .nm_number_code(payload$lengthscale[[column]]), ")^2")
    }, character(1))
    paste0(.nm_number_code(payload$alpha[[row]]), " * exp(-0.5 * (",
           paste(distance, collapse = " + "), "))")
  }, character(1))
  paste0(component$outputs[[1L]], " = ", .nm_number_code(payload$mean), " + ",
         .nm_number_code(payload$variance), " * (", paste(kernels, collapse = " + "), ")")
}

.nm_components <- function(components) {
  if (is.null(components)) return(list())
  if (inherits(components, "nm_component")) components <- list(components)
  if (!is.list(components) || any(!vapply(components, inherits, logical(1), "nm_component"))) {
    .nm_stop("COMPONENTS must contain nm_component objects.")
  }
  names <- vapply(components, `[[`, character(1), "name")
  outputs <- unlist(lapply(components, `[[`, "outputs"), use.names = FALSE)
  if (anyDuplicated(names) || anyDuplicated(outputs)) {
    .nm_stop("Component names and generated outputs must be unique.")
  }
  unname(components)
}

#' @export
print.nm_component <- function(x, ...) {
  cat("LibeRation offline component\n")
  cat("  ", x$name, " [", x$type, ", ", x$scope %||% "pred", "] -> ",
      paste(x$outputs, collapse = ", "), "\n", sep = "")
  cat("  hash:", x$hash, "\n")
  invisible(x)
}
