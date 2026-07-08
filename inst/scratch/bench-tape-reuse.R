ad_root <- Sys.getenv("LIBERTAD_ROOT", unset = normalizePath("../LibeRtAD", winslash = "/"))
pkgload::load_all(ad_root, quiet = TRUE)
pkgload::load_all(normalizePath("."), quiet = TRUE)

bench <- function(label, n, fn) {
  fn()
  t <- system.time(for (i in seq_len(n)) fn())[[3]]
  cat(sprintf("%-22s %5d evals  %7.2fs  %6.1f ms/eval\n", label, n, t, 1000 * t / n))
  invisible(t)
}

sim <- nm_synthetic_theo(n_sub = 10L, seed = 1L)
m <- sim$model
d <- sim$data
par <- .nm_init_par(m)
eta <- matrix(0, 10L, .nm_n_eta(m))
pop <- .nm_build_pop_objective(m, d, eta, pk_engine = "cpp")
labels <- .nm_par_labels(m)
tape_key <- .nm_ad_tape_key(m, d, eta, "pop")

at0 <- stats::setNames(as.list(par), labels)
.nm_state$optim_cache <- NULL
LibeRtAD::reset_tape()
.nm_ad_eval_cached(pop$fn, at0, labels, "cpp", need_grad = TRUE, tape_key = tape_key)

pars <- lapply(seq_len(20), function(i) par * (1 + 0.01 * i))

run_reuse <- function() {
  .nm_state$optim_cache <- NULL
  i <- sample.int(length(pars), 1L)
  at <- stats::setNames(as.list(pars[[i]]), labels)
  .nm_ad_eval_cached(pop$fn, at, labels, "cpp", need_grad = TRUE, tape_key = tape_key)
}

run_fresh <- function() {
  .nm_state$optim_cache <- NULL
  i <- sample.int(length(pars), 1L)
  at <- stats::setNames(as.list(pars[[i]]), labels)
  .nm_do_call_autodiff(pop$fn, at, "cpp", tape_key = NULL)
}

options(LibeRtAD.tape_reuse = TRUE)
cat("THEO pop objective, n_sub=10, grad=AD/cpp, 20 random par points\n\n")
t_reuse <- bench("reuse (cached tape)", 20L, run_reuse)
t_fresh <- bench("fresh tape each time", 20L, run_fresh)

at <- stats::setNames(as.list(pars[[5]]), labels)
.nm_state$optim_cache <- NULL
g_reuse <- .nm_ad_eval_cached(pop$fn, at, labels, "cpp", need_grad = TRUE, tape_key = tape_key)
.nm_state$optim_cache <- NULL
g_fresh <- .nm_do_call_autodiff(pop$fn, at, "cpp", tape_key = NULL)
cat(sprintf("\nmax |grad diff| at one point: %.2e\n", max(abs(unname(g_reuse) - unname(g_fresh[labels])))))
cat(sprintf("speedup (reuse vs fresh): %.2fx\n", t_fresh / t_reuse))
