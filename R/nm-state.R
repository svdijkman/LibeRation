#' @keywords internal
.nm_state <- new.env(parent = emptyenv())
.nm_state$use_cpp_pk <- FALSE
.nm_state$cpp_pk_ad_mode <- NA_character_
.nm_state$optim_cache <- NULL
.nm_state$profile <- NULL
