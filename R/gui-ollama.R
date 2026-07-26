.liber_ollama_default_url <- function() {
  "http://127.0.0.1:11434"
}

.liber_ollama_url <- function(value = NULL) {
  value <- trimws(as.character(value %||% .liber_ollama_default_url())[[1L]])
  value <- sub("/+$", "", value)
  if (!grepl(
    "^https?://(?:localhost|127\\.0\\.0\\.1|\\[::1\\])(?::[0-9]{1,5})?$",
    value, ignore.case = TRUE, perl = TRUE
  )) {
    .nm_stop(
      "The Ollama endpoint must use localhost, 127.0.0.1, or [::1]. ",
      "Remote Ollama endpoints are intentionally not supported by LibeRation."
    )
  }
  port <- suppressWarnings(as.integer(sub(
    "^.*:([0-9]+)$", "\\1", value, perl = TRUE
  )))
  if (!is.na(port) && (port < 1L || port > 65535L)) {
    .nm_stop("The Ollama endpoint port must be between 1 and 65535.")
  }
  value
}

.liber_ai_backend_setting <- function(value, ollama_allowed = TRUE) {
  value <- tolower(trimws(as.character(value %||% "webllm")[[1L]]))
  if (!value %in% c("webllm", "ollama")) value <- "webllm"
  if (identical(value, "ollama") && !isTRUE(ollama_allowed)) value <- "webllm"
  value
}

.liber_ollama_model_setting <- function(value, default = "") {
  value <- trimws(as.character(value %||% default)[[1L]])
  if (!nzchar(value)) default else substr(value, 1L, 160L)
}

.liber_loopback_host <- function(value) {
  value <- tolower(trimws(as.character(value %||% "")[[1L]]))
  value <- sub("^https?://", "", value)
  value <- sub("/.*$", "", value)
  if (value %in% c("", "localhost", "127.0.0.1", "::1", "[::1]")) return(TRUE)
  value <- sub(":\\d+$", "", value)
  value %in% c("localhost", "127.0.0.1")
}

.liber_ollama_runtime_allowed <- function(allow_ollama = NULL,
                                           session_workspace = FALSE,
                                           host = NULL) {
  if (identical(allow_ollama, FALSE) || isTRUE(session_workspace)) return(FALSE)
  hosted_environment <- any(nzchar(Sys.getenv(c(
    "RSCONNECT_USER", "RSCONNECT_SERVER", "CONNECT_SERVER",
    "SHINY_SERVER_VERSION", "SHINYAPPS_INSTANCE", "POSIT_CONNECT_URL"
  ), unset = "")))
  if (hosted_environment) return(FALSE)
  local_binding <- .liber_loopback_host(
    host %||% getOption("shiny.host", "127.0.0.1")
  )
  if (!local_binding) return(FALSE)
  if (isTRUE(allow_ollama)) return(TRUE)
  local_binding
}

.liber_ollama_session_allowed <- function(session, runtime_allowed = FALSE) {
  if (!isTRUE(runtime_allowed)) return(FALSE)
  if (inherits(session, "MockShinySession")) return(TRUE)
  if (is.null(session$request)) return(FALSE)
  request <- session$request
  forwarded <- c(
    request$HTTP_X_FORWARDED_FOR %||% "",
    request$HTTP_FORWARDED %||% "",
    request$HTTP_X_REAL_IP %||% ""
  )
  if (any(nzchar(trimws(as.character(forwarded))))) return(FALSE)
  .liber_loopback_host(request$REMOTE_ADDR %||% "")
}

.liber_ollama_models <- function(base_url, models = ellmer::models_ollama) {
  base_url <- .liber_ollama_url(base_url)
  result <- models(base_url = base_url)
  if (!is.data.frame(result) || !"id" %in% names(result)) {
    .nm_stop("Ollama returned an invalid model catalogue.")
  }
  ids <- unique(trimws(as.character(result$id)))
  ids <- ids[nzchar(ids)]
  lapply(ids, function(id) {
    row <- result[match(id, result$id), , drop = FALSE]
    size <- if ("size" %in% names(row)) suppressWarnings(as.numeric(row$size[[1L]])) else NA_real_
    list(
      id = id,
      label = if (is.finite(size) && size > 0) {
        sprintf("%s (%.1f GB)", id, size / 1024^3)
      } else {
        id
      },
      description = paste0(
        "Installed Ollama model. Inference runs through the local R session ",
        "and the loopback-only Ollama service."
      ),
      size_bytes = if (is.finite(size)) size else NULL
    )
  })
}

.liber_ollama_default_model <- function(models) {
  if (!length(models)) return("")
  ids <- vapply(models, `[[`, character(1), "id")
  lower <- tolower(ids)
  has_size <- grepl("[0-9]+(?:\\.[0-9]+)?b(?:[^a-z]|$)", lower, perl = TRUE)
  billions <- rep(NA_real_, length(ids))
  billions[has_size] <- suppressWarnings(as.numeric(sub(
    "^.*?([0-9]+(?:\\.[0-9]+)?)b(?:[^a-z]|$).*$", "\\1",
    lower[has_size], perl = TRUE
  )))
  practical <- which(is.finite(billions) & billions >= 3 & billions <= 14)
  if (length(practical)) {
    return(ids[practical[[which.min(abs(billions[practical] - 8))]]])
  }
  ids[[1L]]
}

.liber_ollama_message_spec <- function(messages) {
  if (!is.list(messages) || !length(messages)) {
    .nm_stop("An Ollama request must contain at least one message.")
  }
  if (length(messages) > 48L) {
    .nm_stop("The Ollama request contains too many conversation messages.")
  }
  normalized <- lapply(messages, function(message) {
    if (!is.list(message)) .nm_stop("Each Ollama message must be an object.")
    role <- tolower(as.character(message$role %||% "user")[[1L]])
    if (!role %in% c("system", "user", "assistant")) {
      .nm_stop("Unsupported Ollama message role: ", role)
    }
    content <- enc2utf8(as.character(message$content %||% "")[[1L]])
    if (nchar(content, type = "bytes") > 160000L) {
      .nm_stop("An Ollama message exceeds the 160 kB safety limit.")
    }
    list(role = role, content = content)
  })
  total <- sum(vapply(
    normalized, function(message) nchar(message$content, type = "bytes"),
    integer(1)
  ))
  if (total > 400000L) .nm_stop("The Ollama prompt exceeds the 400 kB safety limit.")

  system <- vapply(
    Filter(function(message) identical(message$role, "system"), normalized),
    `[[`, character(1), "content"
  )
  conversation <- Filter(function(message) !identical(message$role, "system"), normalized)
  if (!length(conversation)) .nm_stop("The Ollama request has no user message.")
  user_positions <- which(vapply(
    conversation, function(message) identical(message$role, "user"), logical(1)
  ))
  if (!length(user_positions)) .nm_stop("The Ollama request has no user message.")
  last_user <- max(user_positions)
  prompt <- conversation[[last_user]]$content
  history <- if (last_user > 1L) conversation[seq_len(last_user - 1L)] else list()
  list(
    system_prompt = if (length(system)) paste(system, collapse = "\n\n") else NULL,
    history = history,
    prompt = prompt
  )
}

.liber_ollama_chat <- function(messages, model, base_url,
                               max_tokens = 1000L, temperature = 0.1,
                               top_p = 0.8,
                               chat_factory = ellmer::chat_ollama) {
  spec <- .liber_ollama_message_spec(messages)
  model <- .liber_ollama_model_setting(model)
  if (!nzchar(model)) .nm_stop("Select an installed Ollama model first.")
  max_tokens <- suppressWarnings(as.integer(max_tokens %||% 1000L))
  if (!length(max_tokens) || is.na(max_tokens)) max_tokens <- 1000L
  max_tokens <- max(64L, min(4096L, max_tokens))
  temperature <- suppressWarnings(as.numeric(temperature %||% 0.1))
  if (!length(temperature) || !is.finite(temperature)) temperature <- 0.1
  temperature <- max(0, min(2, temperature))
  top_p <- suppressWarnings(as.numeric(top_p %||% 0.8))
  if (!length(top_p) || !is.finite(top_p)) top_p <- 0.8
  top_p <- max(0.01, min(1, top_p))
  chat <- chat_factory(
    system_prompt = spec$system_prompt,
    base_url = .liber_ollama_url(base_url),
    model = model,
    params = ellmer::params(
      max_tokens = max_tokens, temperature = temperature, top_p = top_p
    ),
    echo = "none"
  )
  if (length(spec$history)) {
    turns <- lapply(spec$history, function(message) {
      contents <- list(ellmer::ContentText(message$content))
      if (identical(message$role, "assistant")) {
        ellmer::AssistantTurn(contents)
      } else {
        ellmer::UserTurn(contents)
      }
    })
    chat$set_turns(turns)
  }
  list(chat = chat, prompt = spec$prompt)
}
