.nm_model_contract_fields_v1 <- c(
  "INPUT", "OUTPUT", "ADVAN", "TRANS", "SS", "DOSECMP", "OBSCMP", "PRED",
  "ERROR", "DES", "THETAS", "OMEGAS", "SIGMAS", "COVARIATES", "USE_ODE",
  "ODE_CONTROL", "IOV", "LIK_CONFIG", "SOLVER", "ERROR_TYPE", "GRAPH",
  "LAYOUT", "LANGUAGE"
)

.nm_model_contract_fields_v2 <- c(
  "INPUT", "OUTPUT", "ADVAN", "TRANS", "SS", "DOSECMP", "OBSCMP", "PRED",
  "ERROR", "DES", "ALG", "THETAS", "OMEGAS", "SIGMAS", "COVARIATES",
  "USE_ODE", "ODE_CONTROL", "IOV", "LIK_CONFIG", "HMM_CONFIG",
  "KALMAN_CONFIG", "DDE_CONFIG", "DAE_CONFIG", "RE_CONFIG", "COMPONENTS",
  "EXPERIMENTAL", "OUTCOMES", "SOLVER", "ERROR_TYPE", "GRAPH", "LAYOUT",
  "LANGUAGE"
)

.nm_contract_classes <- c(
  "nm_matrix_model", "nm_lik_config", "nm_mixture", "nm_residual_group",
  "nm_hmm_config", "nm_cthmm_config", "nm_hsmm_config",
  "nm_factorial_hmm_config", "nm_factorial_chain", "nm_kalman_config",
  "nm_switching_state_space_config", "nm_sde_config", "nm_arma_config",
  "nm_dde_config", "nm_dae_config", "nm_re_config", "nm_re_block",
  "nm_component", "nm_experimental_config", "nm_outcomes", "nm_outcome",
  "nm_model_diagram"
)

.nm_contract_encode <- function(value, depth = 0L) {
  if (depth > 50L) .nm_stop("Model contract nesting exceeds 50 levels.")
  if (is.function(value) || is.environment(value) || is.language(value) ||
      typeof(value) %in% c("externalptr", "weakref")) {
    .nm_stop("Model contracts cannot contain executable or pointer-backed objects.")
  }
  if (!is.list(value) || is.data.frame(value)) return(value)
  classes <- setdiff(class(value), c("list"))
  unexpected <- setdiff(classes, .nm_contract_classes)
  if (length(unexpected)) {
    .nm_stop("Unsupported model-contract class: ", paste(unexpected, collapse = ", "), ".")
  }
  values <- lapply(unclass(value), .nm_contract_encode, depth = depth + 1L)
  if (!length(classes)) return(values)
  list(
    `..liber_contract_type` = "classed-list",
    classes = classes,
    named = !is.null(names(value)),
    names = names(value) %||% character(),
    values = unname(values)
  )
}

.nm_contract_decode <- function(value, depth = 0L) {
  if (depth > 50L) .nm_stop("Model contract nesting exceeds 50 levels.")
  if (is.function(value) || is.environment(value) || is.language(value) ||
      typeof(value) %in% c("externalptr", "weakref")) {
    .nm_stop("Model contracts cannot contain executable or pointer-backed objects.")
  }
  if (!is.list(value) || is.data.frame(value)) return(value)
  if (!identical(value$`..liber_contract_type` %||% NULL, "classed-list")) {
    unexpected <- setdiff(class(value), "list")
    if (length(unexpected)) {
      .nm_stop("Model contract contains an unsupported class.")
    }
    return(lapply(value, .nm_contract_decode, depth = depth + 1L))
  }
  required <- c("..liber_contract_type", "classes", "named", "names", "values")
  if (!identical(sort(names(value)), sort(required))) {
    .nm_stop("Malformed classed object in model contract.")
  }
  classes <- as.character(value$classes)
  if (!length(classes) || length(setdiff(classes, .nm_contract_classes))) {
    .nm_stop("Model contract contains an unsupported class.")
  }
  values <- lapply(value$values, .nm_contract_decode, depth = depth + 1L)
  if (isTRUE(value$named)) {
    object_names <- as.character(value$names)
    if (length(object_names) != length(values) || anyDuplicated(object_names)) {
      .nm_stop("Model contract contains invalid object names.")
    }
    names(values) <- object_names
  } else if (length(value$names)) {
    .nm_stop("Unnamed model-contract object unexpectedly contains names.")
  }
  structure(values, class = classes)
}

#' Versioned semantic model contracts
#'
#' These helpers are the single source of truth for transporting a model
#' between LibeR packages. Compiled IR and external pointers are deliberately
#' excluded; the receiving process recompiles the semantic model definition.
#'
#' @param version Contract version. Version 1 is retained for compatibility;
#'   version 2 contains every first-class model configuration in LibeRation.
#' @return Character vector of semantic model fields.
#' @export
nm_model_contract_fields <- function(version = 2L) {
  version <- as.integer(version)
  if (length(version) != 1L || is.na(version)) .nm_stop("Invalid model contract version.")
  switch(
    as.character(version),
    `1` = .nm_model_contract_fields_v1,
    `2` = .nm_model_contract_fields_v2,
    .nm_stop("Unsupported model contract version: ", version, ".")
  )
}

#' @rdname nm_model_contract_fields
#' @param model An [nm_model()] object.
#' @return `nm_model_to_contract()` returns a JSON-compatible semantic contract.
#' @export
nm_model_to_contract <- function(model, version = 2L) {
  if (!inherits(model, "nm_model")) .nm_stop("`model` must be an nm_model.")
  fields <- nm_model_contract_fields(version)
  required <- c("INPUT", "ADVAN", "PRED", "THETAS")
  missing <- setdiff(required, names(model))
  if (length(missing)) {
    .nm_stop("Model is missing required contract field(s): ", paste(missing, collapse = ", "), ".")
  }
  list(
    schema = "liberation.model",
    version = as.integer(version),
    fields = .nm_contract_encode(stats::setNames(
      lapply(fields, function(field) model[[field]]), fields
    ))
  )
}

#' @rdname nm_model_contract_fields
#' @param contract A contract returned by [nm_model_to_contract()].
#' @return `nm_model_from_contract()` returns a newly validated and compiled
#'   `nm_model`.
#' @export
nm_model_from_contract <- function(contract) {
  if (!is.list(contract) || !identical(as.character(contract$schema), "liberation.model")) {
    .nm_stop("Invalid LibeRation model contract.")
  }
  version <- as.integer(contract$version)
  allowed <- nm_model_contract_fields(version)
  fields <- .nm_contract_decode(contract$fields)
  if (!is.list(fields) || is.null(names(fields)) || anyDuplicated(names(fields)) ||
      length(setdiff(names(fields), allowed))) {
    .nm_stop("Model contract contains invalid semantic fields.")
  }
  required <- c("INPUT", "ADVAN", "PRED", "THETAS")
  if (length(setdiff(required, names(fields)))) {
    .nm_stop("Model contract is missing required semantic fields.")
  }
  tryCatch(
    do.call(nm_model, fields),
    error = function(error) .nm_stop("Model contract validation failed: ", conditionMessage(error))
  )
}
