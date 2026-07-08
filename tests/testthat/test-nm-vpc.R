test_that("quantile column name helper is clean", {
  expect_identical(.nm_vpc_qnames(c(0.05, 0.5, 0.95)), c("q5", "q50", "q95"))
  expect_identical(.nm_vpc_qnames(c(0.025, 0.975)), c("q2p5", "q97p5"))
  expect_false(any(grepl("\\s", .nm_vpc_qnames(c(0.1, 0.9)))))
})

test_that("strata labels combine columns and validate", {
  obs <- data.frame(ID = 1:4, SEX = c(0, 1, 0, 1), GRP = c("a", "a", "b", "b"),
                    stringsAsFactors = FALSE)
  expect_equal(.nm_vpc_strata_labels(obs, NULL), rep("all", 4L))
  lab <- .nm_vpc_strata_labels(obs, c("SEX", "GRP"))
  expect_equal(lab[1], "SEX=0, GRP=a")
  expect_equal(length(unique(lab)), 4L)
  expect_error(.nm_vpc_strata_labels(obs, "NOPE"), "not found")
})

test_that("B12: prediction correction warns and sets NA (never silent fallback)", {
  pred <- c(1, 2, NA, 4, NA, NA)
  strb <- c("a", "a", "a", "a", "b", "b")
  expect_warning(
    res <- .nm_vpc_pc_factor(pred, strb),
    "could not be computed"
  )
  # group b (all NA) -> factor NA, listed as failed
  expect_true(all(is.na(res$factor[strb == "b"])))
  expect_true("b" %in% res$failed)
  # finite preds get a finite factor = median(finite)/pred
  med_a <- stats::median(c(1, 2, 4))
  expect_equal(res$factor[1], med_a / 1)
  expect_equal(res$factor[2], med_a / 2)
  # the NA pred inside group a is left NA (not fabricated)
  expect_true(is.na(res$factor[3]))
})

test_that("nm_vpc computes stratified quantiles and simulation bands", {
  skip_on_cran()
  set.seed(101)
  sim <- nm_synthetic_iv1(n_sub = 16L, seed = 5L)
  m <- sim$model
  d <- data.table::as.data.table(sim$data$data)
  d$SEX <- ifelse(as.integer(d$ID) %% 2L == 0L, 1L, 0L)
  dd <- nm_dataset_from_table(d)
  fit <- nm_est(m, dd, method = "FO",
                control = list(compute_inference = FALSE, n_cores = 1L))

  v <- nm_vpc(fit, n_sim = 40L, n_bins = 5L, seed = 7L, n_cores = 1L)
  expect_s3_class(v, "nm_vpc")
  expect_false(v$pc)
  need <- c("obs_q5", "obs_q50", "obs_q95",
            "sim_q50_med", "sim_q5_lo", "sim_q95_hi", "xmed", "n_obs")
  expect_true(all(need %in% names(v$stats)))
  expect_gt(nrow(v$stats), 0L)
  # observed quantiles ordered
  expect_true(all(v$stats$obs_q5 <= v$stats$obs_q50 + 1e-9, na.rm = TRUE))
  expect_true(all(v$stats$obs_q50 <= v$stats$obs_q95 + 1e-9, na.rm = TRUE))
  # simulation band brackets its own median
  expect_true(all(v$stats$sim_q50_lo <= v$stats$sim_q50_med + 1e-9, na.rm = TRUE))
  expect_true(all(v$stats$sim_q50_med <= v$stats$sim_q50_hi + 1e-9, na.rm = TRUE))

  pcv <- nm_pcvpc(fit, n_sim = 40L, n_bins = 4L, strata = "SEX",
                  seed = 7L, n_cores = 1L)
  expect_true(pcv$pc)
  expect_equal(length(unique(pcv$stats$strat)), 2L)
  expect_true(all(c("SEX=0", "SEX=1") %in% pcv$stats$strat))
})
