#' Define a finite mixture model
#'
#' The model's PRED block may refer to `MIXNUM` in tape-safe `ifelse()`
#' expressions. Subject likelihoods are combined with the supplied prior
#' component probabilities using an exact log-sum-exp calculation.
#'
#' @param probability Positive component probabilities.
#' @param label Optional component labels.
#' @return A serializable mixture definition for `nm_lik_config()`.
#' @export
nm_mixture <- function(probability, label = NULL) {
  probability <- as.numeric(probability)
  if (length(probability) < 2L || any(!is.finite(probability)) || any(probability <= 0)) {
    .nm_stop("A mixture requires at least two positive finite probabilities.")
  }
  probability <- probability / sum(probability)
  if (is.null(label)) label <- paste0("MIX", seq_along(probability))
  label <- as.character(label)
  if (length(label) != length(probability) || any(!nzchar(label)) || anyDuplicated(label)) {
    .nm_stop("Mixture labels must be unique, non-empty, and match the probabilities.")
  }
  structure(list(version = 1L, probability = probability, label = label),
            class = "nm_mixture")
}

#' Posterior subject mixture probabilities
#'
#' @param model An `nm_model`, `NMEngine`, or fitted `nm_fit`.
#' @param data Event data; defaults to fitted data for an `nm_fit`.
#' @param theta,eta,sigma Parameter values; fit values are used when available.
#' @return One row per subject with component probabilities and most likely class.
#' @export
nm_mixture_posterior <- function(model, data = NULL, theta = NULL,
                                 eta = NULL, sigma = NULL) {
  fit <- if (inherits(model, "nm_fit")) model else NULL
  if (!is.null(fit)) {
    data <- data %||% fit$data
    theta <- theta %||% fit$theta
    eta <- eta %||% fit$eta
    sigma <- sigma %||% fit$sigma
    model <- fit$model
  }
  engine <- if (inherits(model, "NMEngine")) model else nm_compile(model)
  mixture <- engine$model$LIK_CONFIG$mixtures
  if (is.null(mixture)) .nm_stop("The model does not define a finite mixture.")
  if (is.null(data)) .nm_stop("`data` is required.")
  data <- .nm_engine_data(engine$model, data)
  n_subjects <- length(unique(data$.ID_INDEX))
  n_eta <- .nm_eta_columns(engine$model, data)
  if (is.null(eta)) eta <- matrix(0, n_subjects, n_eta)
  eta <- as.matrix(eta)
  if (!identical(dim(eta), c(n_subjects, n_eta))) {
    .nm_stop("`eta` has the wrong subject/effect dimensions.")
  }
  theta <- theta %||% engine$model$THETAS$Value
  sigma <- sigma %||% engine$model$SIGMAS$Value
  nll <- .liberation_mixture_component_nll(
    engine$pointer, data, as.numeric(theta), eta, as.numeric(sigma)
  )
  log_weight <- sweep(-0.5 * nll, 2, log(mixture$probability), "+")
  maximum <- apply(log_weight, 1, max)
  probability <- exp(log_weight - maximum)
  probability <- probability / rowSums(probability)
  colnames(probability) <- paste0("P_", make.names(mixture$label, unique = TRUE))
  ids <- attr(data, "id_levels") %||% unique(as.character(data$ID))
  result <- data.frame(ID = ids, probability, check.names = FALSE)
  result$MIXNUM <- max.col(probability, ties.method = "first")
  result$MIXTURE <- mixture$label[result$MIXNUM]
  result
}
