#' Update and revalidate a pharmacometric model
#'
#' Rebuilds an existing [nm_model()] after replacing selected constructor
#' fields. This is the shared editing boundary used by other LibeR packages:
#' source compilation, parameter bounds, covariance structure, outcome
#' configuration, and route compatibility are revalidated by [nm_model()]
#' rather than duplicated in each user interface.
#'
#' @param model A validated [nm_model()].
#' @param ... Named [nm_model()] constructor fields to replace, such as
#'   `PK_SOURCE`, `PRED_SOURCE`, `DES`, `ERROR`, `THETAS`, `OMEGAS`, or
#'   `SIGMAS`.
#' @param name Optional display name. Omit to retain the current model name.
#' @return A rebuilt, validated `nm_model`.
#' @export
nm_model_update <- function(model, ..., name = NULL) {
  if (!inherits(model, "nm_model")) {
    .nm_stop("`model` must be an nm_model.")
  }
  updates <- list(...)
  if (length(updates) && (is.null(names(updates)) ||
      any(!nzchar(names(updates))) || anyDuplicated(names(updates)))) {
    .nm_stop("Model updates must be uniquely named.")
  }
  constructor_fields <- names(formals(nm_model))
  unknown <- setdiff(names(updates), constructor_fields)
  if (length(unknown)) {
    .nm_stop(
      "Unknown nm_model field(s): ", paste(unknown, collapse = ", "), "."
    )
  }
  arguments <- model[intersect(names(model), constructor_fields)]
  has_error_update <- "ERROR" %in% names(updates)
  outcomes <- model$OUTCOMES %||% list()
  general_ctmc <- length(outcomes) && any(vapply(
    outcomes,
    function(outcome) {
      identical(outcome$family, "continuous_time_markov") &&
        length(outcome$initial %||% numeric()) > 2L
    },
    logical(1)
  ))
  if (general_ctmc) {
    arguments$HMM_CONFIG <- NULL
    if (has_error_update) {
      .nm_stop(
        "The likelihood code for a general continuous-time Markov OUTCOME ",
        "is generated from its declaration and cannot be edited directly."
      )
    }
  }
  generated_outcome_error <- isTRUE(model$outcome_error_generated)
  generated_hmm_error <- !is.null(attr(
    model$HMM_CONFIG, "generated_error", exact = TRUE
  ))
  generated_kalman_error <- !is.null(attr(
    model$KALMAN_CONFIG, "generated_error", exact = TRUE
  ))
  if (!has_error_update && (
      generated_outcome_error || generated_hmm_error ||
        generated_kalman_error
    )) {
    arguments$ERROR <- NULL
  }
  if (has_error_update && generated_hmm_error) {
    attr(arguments$HMM_CONFIG, "generated_error") <- NULL
  }
  if (has_error_update && generated_kalman_error) {
    attr(arguments$KALMAN_CONFIG, "generated_error") <- NULL
  }
  arguments[names(updates)] <- updates
  updated <- do.call(nm_model, arguments)
  code <- c(
    updated$PK_SOURCE %||% "", updated$PRED_SOURCE %||% "",
    updated$PRED %||% "", updated$DES %||% "", updated$ALG %||% "",
    updated$ERROR %||% ""
  )
  requirements <- c(
    THETA = .liber_code_reference_max(code, "THETA"),
    ETA = .liber_code_reference_max(code, "ETA"),
    SIGMA = .liber_code_reference_max(code, c("ERR", "EPS", "SIGMA"))
  )
  available <- c(
    THETA = nrow(updated$THETAS),
    ETA = if (nrow(updated$OMEGAS)) {
      max(c(updated$OMEGAS$ROW, updated$OMEGAS$COL), na.rm = TRUE)
    } else {
      0L
    },
    SIGMA = nrow(updated$SIGMAS)
  )
  missing <- names(requirements)[requirements > available]
  if (length(missing)) {
    detail <- vapply(missing, function(parameter) {
      paste0(
        parameter, "(", requirements[[parameter]],
        ") is referenced but only ", available[[parameter]],
        " ", parameter, " parameter row(s) are defined"
      )
    }, character(1))
    .nm_stop(
      paste(detail, collapse = "; "),
      ". Add the required parameter rows and apply the model again."
    )
  }
  preserved <- setdiff(
    names(attributes(model)),
    c("names", "class", "row.names")
  )
  for (attribute in preserved) {
    attr(updated, attribute) <- attr(model, attribute, exact = TRUE)
  }
  current_name <- attr(model, "name", exact = TRUE)
  if (!is.null(name)) {
    current_name <- trimws(as.character(name))
    if (length(current_name) != 1L || is.na(current_name) ||
        !nzchar(current_name)) {
      .nm_stop("`name` must be one non-empty string.")
    }
  }
  if (!is.null(current_name)) attr(updated, "name") <- current_name
  updated
}
