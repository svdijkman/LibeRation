#' @keywords internal
.nm_stop <- function(..., call. = FALSE) {
  stop(..., call. = call.)
}

#' @keywords internal
.nm_require_cols <- function(dat, cols, label = "data") {
  missing <- setdiff(cols, names(dat))
  if (length(missing) > 0L) {
    .nm_stop("Missing columns in ", label, ": ", paste(missing, collapse = ", "))
  }
}

#' @keywords internal
.nm_split_lines <- function(code) {
  if (is.null(code) || !nzchar(trimws(code))) {
    return(character())
  }
  lines <- unlist(strsplit(code, "\n", fixed = TRUE))
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  lines
}

#' @keywords internal
.nm_logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

#' @keywords internal
.nm_split_est_control <- function(control) {
  print_grad_every <- control$print_grad_every
  if (is.null(print_grad_every)) {
    print_grad_every <- 0L
  }
  print_grad_every <- as.integer(print_grad_every[1L])
  if (print_grad_every < 0L) {
    .nm_stop("control$print_grad_every must be >= 0.")
  }
  n_cores <- control$n_cores
  if (is.null(n_cores)) {
    n_cores <- .nm_default_n_cores()
  }
  n_cores <- as.integer(n_cores[1L])
  if (n_cores < 1L) {
    .nm_stop("control$n_cores must be >= 1.")
  }
  optim_control <- control[!names(control) %in% c(
    "print_grad_every", "n_cores", "compute_inference", "compute_covariance",
    "cov_method", "cov_refit_eta", "infer_hessian",
    "min_retries", "tweak_inits", "par_scale", "mceta", "maxit_eta", "maxit_eta_warm",
    "focei_grad", "focei_eta"
  )]
  list(
    print_grad_every = print_grad_every,
    n_cores = n_cores,
    optim_control = optim_control
  )
}

#' @keywords internal
.nm_default_n_cores <- function() {
  opt <- getOption("LibeRtAD.n_cores", NULL)
  if (!is.null(opt)) {
    return(as.integer(max(1L, opt[1L])))
  }
  if (!isTRUE(getOption("LibeRtAD.n_cores_auto", TRUE))) {
    return(1L)
  }
  if (requireNamespace("parallel", quietly = TRUE)) {
    det <- suppressWarnings(parallel::detectCores(logical = FALSE))
    if (!is.na(det) && det >= 2L) {
      return(as.integer(max(1L, det - 1L)))
    }
  }
  1L
}

#' @keywords internal
.nm_resolve_n_cores <- function(control = list()) {
  n <- control$n_cores
  if (is.null(n)) {
    return(.nm_default_n_cores())
  }
  as.integer(max(1L, n[1L]))
}

#' @keywords internal
.nm_pkg_root <- function(pkg) {
  if (!pkg %in% loadedNamespaces()) {
    return(NULL)
  }
  ip <- system.file("", package = pkg)
  if (!nzchar(ip)) {
    return(NULL)
  }
  normalizePath(ip, winslash = "/", mustWork = FALSE)
}

#' @keywords internal
.nm_pkg_is_installed <- function(pkg) {
  root <- .nm_pkg_root(pkg)
  !is.null(root) && file.exists(file.path(root, "Meta", "package.rds"))
}

#' @keywords internal
.nm_optim_counts_na <- function() {
  stats::setNames(c(NA_integer_, NA_integer_), c("function", "gradient"))
}

#' @keywords internal
.nm_optim_counts_one <- function(step) {
  if (is.null(step) || is.null(step$counts)) {
    return(.nm_optim_counts_na())
  }
  step$counts
}

#' @keywords internal
.nm_sum_optim_counts <- function(optim) {
  if (is.null(optim)) {
    return(.nm_optim_counts_na())
  }
  if (!is.null(optim$counts)) {
    return(.nm_optim_counts_one(optim))
  }
  if ("step1" %in% names(optim)) {
    c1 <- .nm_optim_counts_one(optim$step1)
    c2 <- .nm_optim_counts_one(optim$step2)
    return(stats::setNames(
      c(
        as.integer(sum(c1[["function"]], c2[["function"]], na.rm = TRUE)),
        as.integer(sum(c1[["gradient"]], c2[["gradient"]], na.rm = TRUE))
      ),
      c("function", "gradient")
    ))
  }
  if (is.list(optim) && length(optim) > 0L) {
    first <- optim[[1L]]
    if (is.list(first) && !is.null(first$counts)) {
      fn <- 0
      gr <- 0
      for (o in optim) {
        c0 <- .nm_optim_counts_one(o)
        fn <- fn + c0[["function"]]
        gr <- gr + c0[["gradient"]]
      }
      return(stats::setNames(
        c(as.integer(fn), as.integer(gr)),
        c("function", "gradient")
      ))
    }
  }
  .nm_optim_counts_one(optim)
}

#' @keywords internal
.nm_parallel_start_cluster <- function(n_cores) {
  ad_root <- .nm_pkg_root("LibeRtAD")
  nm_root <- .nm_pkg_root("LibeRation")
  if (is.null(ad_root) || is.null(nm_root)) {
    return(NULL)
  }
  use_load_all <- !.nm_pkg_is_installed("LibeRtAD") || !.nm_pkg_is_installed("LibeRation")
  cl <- parallel::makePSOCKcluster(n_cores)
  if (use_load_all) {
    if (!requireNamespace("pkgload", quietly = TRUE)) {
      parallel::stopCluster(cl)
      return(NULL)
    }
    parallel::clusterExport(cl, c("ad_root", "nm_root"), envir = environment())
    errs <- parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        library(LibeRtAD)
        pkgload::load_all(ad_root, quiet = TRUE, export_all = TRUE, compile = FALSE)
        pkgload::load_all(nm_root, quiet = TRUE, export_all = TRUE, compile = FALSE)
      })
      TRUE
    })
    if (length(errs) != n_cores || !all(vapply(errs, isTRUE, logical(1L)))) {
      parallel::stopCluster(cl)
      return(NULL)
    }
    return(cl)
  }
  pkg_libs <- unique(c(dirname(ad_root), dirname(nm_root)))
  parallel::clusterExport(cl, "pkg_libs", envir = environment())
  errs <- parallel::clusterEvalQ(cl, {
    .libPaths(c(pkg_libs, .libPaths()))
    suppressPackageStartupMessages(library(LibeRtAD))
    suppressPackageStartupMessages(library(LibeRation))
    TRUE
  })
  if (length(errs) != n_cores || !all(vapply(errs, isTRUE, logical(1L)))) {
    parallel::stopCluster(cl)
    return(NULL)
  }
  cl
}

#' @keywords internal
.nm_saem_mstep_control <- function(control, k, n_burn) {
  ctl <- control
  maxit <- control$maxit
  if (is.null(maxit)) {
    maxit <- 100L
  }
  if (k <= n_burn) {
    burn <- control$maxit_burn
    if (is.null(burn)) {
      burn <- max(5L, as.integer(floor(maxit / 3)))
    }
    ctl$maxit <- as.integer(burn[1L])
  }
  if (is.null(ctl$factr)) {
    ctl$factr <- 1e5
  }
  ctl
}

#' @keywords internal
.nm_saem_mh_one_r <- function(model, subj, eta_cur, theta, omega, sigma,
                              n_mcmc, pk_engine) {
  n_eta <- length(eta_cur)
  n_old <- .nm_conditional_nll(
    model, subj, theta, omega, sigma, eta_cur, pk_engine
  )
  for (m in seq_len(n_mcmc)) {
    prop <- eta_cur + rnorm(n_eta, sd = sqrt(pmax(omega, 1e-8)))
    n_new <- .nm_conditional_nll(
      model, subj, theta, omega, sigma, prop, pk_engine
    )
    if (log(runif(1)) < n_old - n_new) {
      eta_cur <- prop
      n_old <- n_new
    }
  }
  eta_cur
}

#' @keywords internal
.nm_saem_mh_r <- function(model, dat, ids, eta_mat, theta, omega, sigma,
                          n_mcmc, pk_engine, control = list()) {
  n_sub <- length(ids)
  n_eta <- ncol(eta_mat)
  if (n_eta == 0L) {
    return(eta_mat)
  }
  n_cores <- .nm_resolve_n_cores(control)
  fit_one <- function(job) {
    .nm_saem_mh_one_r(
      job$model, job$subj, job$eta0, job$theta, job$omega, job$sigma,
      job$n_mcmc, job$pk_engine
    )
  }
  jobs <- lapply(seq_len(n_sub), function(j) {
    list(
      model = model,
      subj = .nm_subject_slice(dat, ids[j]),
      eta0 = eta_mat[j, ],
      theta = theta,
      omega = omega,
      sigma = sigma,
      n_mcmc = n_mcmc,
      pk_engine = pk_engine
    )
  })
  if (n_cores > 1L && n_sub > 1L) {
    res <- .nm_parallel_lapply(jobs, fit_one, n_cores = n_cores)
    for (j in seq_len(n_sub)) {
      eta_mat[j, ] <- res[[j]]
    }
    return(eta_mat)
  }
  for (j in seq_len(n_sub)) {
    eta_mat[j, ] <- fit_one(jobs[[j]])
  }
  eta_mat
}

#' @keywords internal
.nm_parallel_lapply <- function(X, FUN, n_cores = 1L) {
  n_cores <- min(as.integer(n_cores), length(X))
  if (n_cores <= 1L || length(X) <= 1L) {
    return(lapply(X, FUN))
  }
  if (.Platform$OS.type == "windows") {
    cl <- .nm_parallel_start_cluster(n_cores)
    if (is.null(cl)) {
      warning(
        "Parallel subject fits unavailable (install LibeRtAD and LibeRation, ",
        "or use devtools::load_all() with pkgload). ",
        "Running sequentially; set control$n_cores = 1 to silence.",
        call. = FALSE
      )
      return(lapply(X, FUN))
    }
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::parLapplyLB(cl, X, FUN)
  } else {
    parallel::mclapply(X, FUN, mc.cores = n_cores)
  }
}

#' @keywords internal
.nm_print_grad <- function(step, grad, prefix = "") {
  if (!nzchar(prefix)) {
    cat("Gradient evaluation ", step, ":\n", sep = "")
  } else {
    cat(prefix, " gradient evaluation ", step, ":\n", sep = "")
  }
  nm <- names(grad)
  for (i in seq_along(grad)) {
    cat(sprintf("  %-8s % .6g\n", nm[i], grad[i]))
  }
  invisible(grad)
}

#' @keywords internal
.nm_wrap_grad_trace <- function(gr_fn, print_every, grad_names, prefix = "") {
  if (print_every <= 0L) {
    return(gr_fn)
  }
  state <- new.env(parent = emptyenv())
  state$count <- 0L
  state$gr_fn <- gr_fn
  state$print_every <- print_every
  state$grad_names <- grad_names
  state$prefix <- prefix
  function(x) {
    g <- state$gr_fn(x)
    state$count <- state$count + 1L
    if (state$count %% state$print_every == 0L) {
      .nm_print_grad(
        state$count,
        stats::setNames(as.numeric(g), state$grad_names),
        prefix = state$prefix
      )
    }
    g
  }
}

#' @keywords internal
.nm_num_grad <- function(f, x, eps = 1e-5) {
  g <- numeric(length(x))
  fx <- f(x)
  if (!is.finite(fx)) {
    fx <- .Machine$double.xmax
  }
  for (i in seq_along(x)) {
    xp <- x
    xm <- x
    xp[i] <- xp[i] + eps
    xm[i] <- xm[i] - eps
    fp <- f(xp)
    fm <- f(xm)
    if (!is.finite(fp)) {
      fp <- fx
    }
    if (!is.finite(fm)) {
      fm <- fx
    }
    g[i] <- (fp - fm) / (2 * eps)
    if (!is.finite(g[i])) {
      g[i] <- 0
    }
  }
  g
}

#' @keywords internal
.nm_stable_vcov <- function(H) {
  if (is.null(H) || length(H) == 0L) {
    return(NULL)
  }
  H <- (H + t(H)) / 2
  n <- nrow(H)
  d <- diag(H)
  pos <- d[is.finite(d) & d > 0]
  base <- if (length(pos) > 0L) stats::median(pos) else 1
  if (!is.finite(base) || base <= 0) {
    base <- 1
  }
  ridge <- max(1e-8, 1e-3 * base)
  vcov <- tryCatch(
    solve(H + diag(ridge, n)),
    error = function(e) NULL
  )
  if (!is.null(vcov)) {
    return(vcov)
  }
  ev <- tryCatch(eigen(H, symmetric = TRUE), error = function(e) NULL)
  if (is.null(ev)) {
    return(NULL)
  }
  vals <- pmax(ev$values, ridge)
  Hreg <- ev$vectors %*% diag(vals, n) %*% t(ev$vectors)
  vcov <- tryCatch(
    solve(Hreg + diag(ridge, n)),
    error = function(e) NULL
  )
  if (!is.null(vcov)) {
    return(vcov)
  }
  if (requireNamespace("MASS", quietly = TRUE)) {
    return(tryCatch(MASS::ginv(Hreg), error = function(e) NULL))
  }
  NULL
}

#' @keywords internal
.nm_num_hessian <- function(f, x, eps = 1e-5) {
  n <- length(x)
  g0 <- .nm_num_grad(f, x, eps = eps)
  H <- matrix(0, n, n)
  for (i in seq_len(n)) {
    xp <- x
    xp[i] <- xp[i] + eps
    H[i, ] <- (.nm_num_grad(f, xp, eps = eps) - g0) / eps
  }
  (H + t(H)) / 2
}

#' @keywords internal
.nm_effective_n_quad <- function(model, n_quad) {
  n_quad <- as.integer(n_quad)
  n_eta <- .nm_n_eta(model)
  if (n_eta >= 3L && isTRUE(getOption("LibeRtAD.laplace_auto_n_quad", TRUE))) {
    cap <- as.integer(getOption("LibeRtAD.laplace_n_quad_cap", 3L))
    min(n_quad, cap)
  } else {
    n_quad
  }
}

#' @keywords internal
.nm_gh_nodes <- function(n) {
  n <- as.integer(n)
  if (n < 1L || n > 11L) {
    .nm_stop("Gaussian quadrature order must be between 1 and 11.")
  }
  if (requireNamespace("statmod", quietly = TRUE)) {
    gh <- statmod::gauss.quad.prob(n, dist = "normal")
    return(list(nodes = gh$nodes, weights = gh$weights))
  }
  # Fallback: precomputed Gauss-Hermite for standard normal (n = 1, 3, 5, 7, 9, 11)
  tables <- list(
    `1` = list(nodes = 0, weights = 1),
    `3` = list(
      nodes = c(-1.224744871391589, 0, 1.224744871391589),
      weights = c(0.295408973912104, 1.18163590060369, 0.295408973912104)
    ),
    `5` = list(
      nodes = c(-2.020182870456086, -0.958572464698798, 0,
                0.958572464698798, 2.020182870456086),
      weights = c(0.019953242059045, 0.393619323152241, 0.945308720482942,
                  0.393619323152241, 0.019953242059045)
    ),
    `7` = list(
      nodes = c(-2.571356011682409, -1.355626019367357, -0.342901327223704, 0,
                0.342901327223704, 1.355626019367357, 2.571356011682409),
      weights = c(0.002985045148820, 0.099090322734308, 0.357781605736023,
                  0.641880547310781, 0.357781605736023, 0.099090322734308,
                  0.002985045148820)
    ),
    `9` = list(
      nodes = c(-3.029899444060, -1.899171040454, -0.911266815343, 0,
                0.911266815343, 1.899171040454, 3.029899444060),
      weights = c(0.000971781245, 0.033050847, 0.240138611, 0.436077448,
                  0.240138611, 0.033050847, 0.000971781245)
    ),
    `11` = list(
      nodes = c(-3.442257297, -2.161761757, -1.229847829, -0.373333990, 0,
                0.373333990, 1.229847829, 2.161761757, 3.442257297),
      weights = c(0.000266431, 0.008530556, 0.075942449, 0.283830357,
                  0.462910047, 0.283830357, 0.075942449, 0.008530556,
                  0.000266431)
    )
  )
  if (!as.character(n) %in% names(tables)) {
    .nm_stop("Install package 'statmod' for quadrature order ", n, ".")
  }
  tables[[as.character(n)]]
}

#' @keywords internal
.nm_product_grid <- function(nodes, weights, n_eta) {
  if (n_eta == 1L) {
    return(list(
      nodes = matrix(nodes, ncol = 1L),
      weights = weights
    ))
  }
  grids <- expand.grid(replicate(n_eta, nodes, simplify = FALSE))
  w <- rep(1, nrow(grids))
  for (j in seq_len(n_eta)) {
    w <- w * weights[match(grids[[j]], nodes)]
  }
  list(nodes = as.matrix(grids), weights = w)
}
