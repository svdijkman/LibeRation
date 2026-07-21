#' Define a sparse mechanistic QSP reaction system
#'
#' The stoichiometric network is compiled into the same C++/CppAD `$DES`
#' representation as an ordinary pharmacometric ODE. Optional algebraic
#' constraints use the block-sparse index-1 DAE solver.
#'
#' @param species Unique species/state names.
#' @param stoichiometry Species-by-reaction numeric matrix.
#' @param rates One restricted expression per reaction. Species names in these
#'   expressions are translated to their state amounts.
#' @param dose_species,observation_species Species name or index.
#' @param algebraic Optional `ALG` residual code.
#' @param dae_config Optional [nm_dae_config()].
#' @export
nm_qsp_system <- function(species, stoichiometry, rates,
                          dose_species = 1L, observation_species = 1L,
                          algebraic = "", dae_config = NULL) {
  species <- trimws(as.character(species))
  if (!length(species) || anyNA(species) || any(!nzchar(species)) ||
      anyDuplicated(species) || any(!grepl("^[A-Za-z][A-Za-z0-9_]*$", species))) {
    .nm_stop("QSP species must be unique identifier names.")
  }
  stoichiometry <- as.matrix(stoichiometry)
  storage.mode(stoichiometry) <- "double"
  rates <- trimws(as.character(rates))
  if (nrow(stoichiometry) != length(species) || ncol(stoichiometry) != length(rates) ||
      !length(rates) || any(!is.finite(stoichiometry)) || anyNA(rates) || any(!nzchar(rates))) {
    .nm_stop("QSP stoichiometry and rate dimensions/values are inconsistent.")
  }
  resolve_species <- function(value, label) {
    if (is.character(value)) value <- match(value, species)
    value <- as.integer(value)
    if (length(value) != 1L || is.na(value) || value < 1L || value > length(species)) {
      .nm_stop("`", label, "` must identify one QSP species.")
    }
    value
  }
  translated_rate <- vapply(rates, function(rate) {
    for (index in order(nchar(species), decreasing = TRUE)) {
      rate <- gsub(paste0("\\b", species[[index]], "\\b"),
                   paste0("A(", index, ")"), rate, perl = TRUE)
    }
    rate
  }, character(1))
  number <- function(value) format(value, digits = 17L, scientific = TRUE)
  des <- vapply(seq_along(species), function(state) {
    nonzero <- which(stoichiometry[state, ] != 0)
    expression <- if (!length(nonzero)) "0" else paste(vapply(nonzero, function(reaction) {
      paste0("(", number(stoichiometry[state, reaction]), ") * (",
             translated_rate[[reaction]], ")")
    }, character(1)), collapse = " + ")
    paste0("DADT(", state, ") = ", expression)
  }, character(1))
  graph <- list(
    compartments = data.frame(
      id = seq_along(species), name = species, state = paste0("A", seq_along(species)),
      stringsAsFactors = FALSE),
    source = "QSP", stoichiometry = stoichiometry,
    reactions = colnames(stoichiometry) %||% paste0("reaction", seq_along(rates))
  )
  structure(list(
    version = 1L, species = species, stoichiometry = stoichiometry,
    rates = rates, translated_rates = translated_rate, DES = paste(des, collapse = "\n"),
    ALG = paste(algebraic, collapse = "\n"), DAE_CONFIG = .nm_dae_config(dae_config),
    DOSECMP = resolve_species(dose_species, "dose_species"),
    OBSCMP = resolve_species(observation_species, "observation_species"),
    GRAPH = graph
  ), class = "nm_qsp_system")
}

#' Build a LibeRation model from a QSP reaction system
#'
#' @param system An [nm_qsp_system()].
#' @param ... Arguments passed to [nm_model()], including `INPUT`, `PRED`,
#'   parameters, likelihood, and solver controls.
#' @param EXPERIMENTAL Explicit experimental-engine acknowledgement.
#' @export
nm_qsp_model <- function(system, ..., EXPERIMENTAL = NULL) {
  if (!inherits(system, "nm_qsp_system")) {
    .nm_stop("`system` must be created by `nm_qsp_system()`.")
  }
  experimental <- .nm_experimental_config(EXPERIMENTAL, "QSP reaction networks")
  supplied <- list(...)
  reserved <- intersect(names(supplied), c(
    "ADVAN", "DES", "ALG", "DAE_CONFIG", "DOSECMP", "OBSCMP", "GRAPH", "EXPERIMENTAL"))
  if (length(reserved)) {
    .nm_stop("QSP model arguments are supplied by `system`: ", paste(reserved, collapse = ", "), ".")
  }
  model <- do.call(nm_model, c(supplied, list(
    ADVAN = 6L, DES = system$DES, ALG = system$ALG,
    DAE_CONFIG = system$DAE_CONFIG, DOSECMP = system$DOSECMP,
    OBSCMP = system$OBSCMP, GRAPH = system$GRAPH,
    EXPERIMENTAL = experimental
  )))
  model$EXPERIMENTAL$features <- sort(unique(c(
    model$EXPERIMENTAL$features, "QSP reaction networks")))
  model$QSP_SYSTEM <- system
  model
}

#' @export
print.nm_qsp_system <- function(x, ...) {
  cat("LibeRation experimental QSP reaction system\n")
  cat("  species:", length(x$species), " reactions:", ncol(x$stoichiometry),
      " algebraic variables:", length(x$DAE_CONFIG$variables %||% character()), "\n")
  invisible(x)
}
