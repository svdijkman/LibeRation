# ---------------------------------------------------------------------------
# Worker-side crypto for LibeRties-scheduled jobs. Mirrors the envelope format
# used by LibeRties/R/crypto.R so encrypted payloads interoperate. The DEK is
# supplied in memory via the LIBERTIES_JOB_DEK environment variable (hex).
#
# These helpers are only used when running under the LibeRties scheduler with
# at-rest encryption enabled. Local (single-tenant) jobs keep using plaintext
# args.rds / result.rds.
# ---------------------------------------------------------------------------

#' @keywords internal
.nm_job_dek <- function() {
  hex <- Sys.getenv("LIBERTIES_JOB_DEK", "")
  if (!nzchar(hex) || !requireNamespace("sodium", quietly = TRUE)) {
    return(NULL)
  }
  tryCatch(sodium::hex2bin(hex), error = function(e) NULL)
}

#' @keywords internal
.nm_encrypt_to_file <- function(obj, key, path) {
  msg <- serialize(obj, connection = NULL)
  nonce <- sodium::random(24L)
  ct <- sodium::data_encrypt(msg, key, nonce)
  attr(ct, "nonce") <- NULL
  env <- list(v = 1L, alg = "secretbox", nonce = nonce, ct = as.raw(ct))
  .nm_save_rds_safe(env, path)
}

#' @keywords internal
.nm_decrypt_from_file <- function(key, path) {
  env <- .nm_read_rds_safe(path)
  if (is.null(env)) {
    stop("Encrypted file missing or unreadable: ", path, call. = FALSE)
  }
  unserialize(sodium::data_decrypt(env$ct, key, env$nonce))
}

#' Read job args, decrypting args.enc when a DEK is present.
#' @keywords internal
.nm_job_read_args <- function(job_path) {
  enc <- file.path(job_path, "args.enc")
  dek <- .nm_job_dek()
  if (file.exists(enc) && !is.null(dek)) {
    return(.nm_decrypt_from_file(dek, enc))
  }
  .nm_read_rds_safe(file.path(job_path, "args.rds"))
}

#' Write a job result, encrypting to result.enc when a DEK is present.
#' @keywords internal
.nm_job_write_result <- function(obj, job_path) {
  dek <- .nm_job_dek()
  if (!is.null(dek)) {
    .nm_encrypt_to_file(obj, dek, file.path(job_path, "result.enc"))
    return(invisible(TRUE))
  }
  .nm_save_rds_safe(obj, file.path(job_path, "result.rds"))
}
