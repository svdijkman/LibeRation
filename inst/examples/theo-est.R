# THEO-like estimation example (ADVAN 3 TRANS 4, C++ PKADVAN router)
# Load LibeRtAD then LibeRation from source when possible.
.nm_pkg_name_from_desc <- function(desc_path) {
  as.character(read.dcf(desc_path, fields = "Package")[1L, 1L])
}

theo_load_packages <- function() {
  load_pkg <- function(root) {
    if (requireNamespace("devtools", quietly = TRUE) &&
        file.exists(file.path(root, "DESCRIPTION"))) {
      devtools::load_all(root, quiet = TRUE)
      return(TRUE)
    }
    FALSE
  }
  if (file.exists("DESCRIPTION") &&
      identical(.nm_pkg_name_from_desc("DESCRIPTION"), "LibeRation")) {
    for (root in c("../LibeRtAD", "../../LibeRtAD")) {
      load_pkg(root)
    }
    load_pkg(".")
    return(invisible(TRUE))
  }
  for (root in c("..", "../..")) {
    desc <- file.path(root, "DESCRIPTION")
    if (file.exists(desc) &&
        identical(.nm_pkg_name_from_desc(desc), "LibeRation")) {
      for (ad in c(file.path(root, "LibeRtAD"), file.path(dirname(root), "LibeRtAD"))) {
        load_pkg(ad)
      }
      load_pkg(root)
      return(invisible(TRUE))
    }
  }
  if (!requireNamespace("LibeRtAD", quietly = TRUE)) {
    stop("Install LibeRtAD first, or run from a source checkout with devtools.")
  }
  if (!requireNamespace("LibeRation", quietly = TRUE)) {
    stop("Install LibeRation first, or run from a source checkout with devtools.")
  }
  library(LibeRtAD)
  library(LibeRation)
  invisible(FALSE)
}
theo_load_packages()

theo_bench_row <- function(label, fit, elapsed_sec) {
  if (identical(fit$method, "BAYES")) {
    fn_eval <- as.integer(fit$eval$pk)
    gr_eval <- NA_integer_
  } else {
    counts <- LibeRation:::.nm_sum_optim_counts(fit$optim)
    fn_eval <- as.integer(counts[["function"]])
    gr_eval <- as.integer(counts[["gradient"]])
  }
  data.frame(
    label = label,
    method = fit$method,
    grad = fit$grad,
    grad_backend = if (is.null(fit$grad_backend)) "" else fit$grad_backend,
    pk_engine = fit$pk_engine,
    engine = {
      e <- fit[["engine"]]
      if (is.null(e) || !is.character(e)) "" else e
    },
    objective = fit$objective,
    convergence = fit$convergence,
    fn_eval = fn_eval,
    gr_eval = gr_eval,
    time_sec = round(elapsed_sec, 2),
    stringsAsFactors = FALSE
  )
}

theo_bench <- function(label, ...) {
  cat("Running:", label, "...\n")
  t0 <- proc.time()
  fit <- nm_est(...)
  elapsed <- (proc.time() - t0)[3]
  theo_bench_row(label, fit, elapsed)
}

sim <- nm_synthetic_theo(n_sub = 10L, seed = 1L)
th0 <- sim$model$THETAS$Value
om0 <- sim$model$OMEGAS$Value
sg0 <- sim$model$SIGMAS$Value

cat("Initial -2LL (C++ PK):",
    nm_nll(sim$model, sim$data, th0, om0, sg0,
           include_omega_prior = FALSE, pk_engine = "cpp"), "\n\n")

est_control <- list(maxit = 10)
bench_rows <- list()

### Numeric gradients + C++ PK

bench_rows[[length(bench_rows) + 1L]] <- theo_bench(
  "FO numeric / cpp PK",
  sim$model, sim$data, method = "FO",
  grad = "numeric", pk_engine = "cpp", control = est_control
)

bench_rows[[length(bench_rows) + 1L]] <- theo_bench(
  "FOCE numeric / cpp PK",
  sim$model, sim$data, method = "FOCE",
  grad = "numeric", pk_engine = "cpp", max_outer = 5L, control = est_control
)

bench_rows[[length(bench_rows) + 1L]] <- theo_bench(
  "SAEM numeric / cpp PK",
  sim$model, sim$data, method = "SAEM",
  engine = "cpp", grad = "numeric", pk_engine = "cpp",
  n_iter = 10L, seed = 1L
)

bench_rows[[length(bench_rows) + 1L]] <- theo_bench(
  "LAPLACE numeric",
  sim$model, sim$data, method = "LAPLACE",
  engine = "cpp", grad = "numeric", pk_engine = "cpp",
  n_quad = 5L, control = est_control
)

bench_rows[[length(bench_rows) + 1L]] <- theo_bench(
  "BAYES MCMC (cpp)",
  sim$model, sim$data, method = "BAYES",
  engine = "cpp", pk_engine = "cpp",
  n_burn = 50L, n_sample = 100L, n_thin = 2L, seed = 1L
)

### AD gradients (grad = "cpp" -> grad = "ad", backend = "cpp")
# PK uses the R PKADVAN step on the AD tape; pk_engine = "cpp" is overridden to R.

bench_rows[[length(bench_rows) + 1L]] <- theo_bench(
  "FO AD (grad=cpp)",
  sim$model, sim$data, method = "FO",
  grad = "cpp", pk_engine = "cpp", control = est_control
)

bench_rows[[length(bench_rows) + 1L]] <- theo_bench(
  "FOCE AD (grad=cpp)",
  sim$model, sim$data, method = "FOCE",
  grad = "cpp", pk_engine = "cpp", max_outer = 3L, control = est_control
)

bench_rows[[length(bench_rows) + 1L]] <- theo_bench(
  "SAEM AD (grad=cpp)",
  sim$model, sim$data, method = "SAEM",
  engine = "cpp", grad = "cpp", pk_engine = "cpp",
  n_iter = 10L, seed = 1L
)

bench_rows[[length(bench_rows) + 1L]] <- theo_bench(
  "LAPLACE (grad=cpp)",
  sim$model, sim$data, method = "LAPLACE",
  grad = "cpp", pk_engine = "cpp",
  n_quad = 5L, control = est_control
)

bench <- do.call(rbind, bench_rows)
rownames(bench) <- NULL

cat("=== Performance comparison ===\n")
print(bench[, c("label", "grad", "grad_backend", "pk_engine", "engine",
                "objective", "fn_eval", "gr_eval", "time_sec")])

cat("\nNote: fn_eval and gr_eval are optim() counts (L-BFGS-B calls fn and gr\n",
    "together each iteration, so they are usually equal). BAYES fn_eval is PK\n",
    "likelihood evaluations during MCMC (gr_eval is NA).\n")
cat("\nNumeric vs AD (same method):\n")
for (meth in c("FO", "FOCE", "SAEM", "LAPLACE")) {
  sub <- bench[bench$method == meth, ]
  num <- sub[sub$grad == "numeric", , drop = FALSE]
  ad <- sub[sub$grad == "ad", , drop = FALSE]
  if (nrow(num) >= 1L && nrow(ad) == 1L) {
    num <- num[1L, , drop = FALSE]
    t_num <- as.numeric(num[["time_sec"]])
    t_ad <- as.numeric(ad[["time_sec"]])
    cat(sprintf(
      "  %s: numeric %.2f s, AD %.2f s (numeric/AD time ratio=%.2f)\n",
      meth, t_num, t_ad, t_num / t_ad
    ))
  }
}

invisible(bench)

