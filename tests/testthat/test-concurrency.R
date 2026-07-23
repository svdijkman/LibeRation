test_that("independent R workers can fit concurrently without sharing tapes", {
  skip_if_not_installed("callr")
  worker <- function() {
    model <- LibeRation::nm_model(
      INPUT = c("ID", "TIME", "EVID", "AMT", "DV", "MDV"),
      ADVAN = 1, DOSECMP = 1, OBSCMP = 1,
      PRED = "CL=THETA(1); V=THETA(2); S1=V",
      ERROR = "Y=F+ERR(1)",
      THETAS = data.frame(
        THETA = 1:2, Value = c(2, 20), FIX = c(FALSE, TRUE)
      ),
      SIGMAS = data.frame(SIGMA = 1, Value = 0.1, FIX = TRUE)
    )
    data <- data.frame(
      ID = 1, TIME = c(0, 1, 4),
      EVID = c(1, 0, 0), AMT = c(100, 0, 0),
      DV = c(NA, 4.5, 3.3), MDV = c(1, 0, 0)
    )
    fit <- LibeRation::nm_est(
      model, data, method = "FO", maxit = 2L, n_cores = 1L
    )
    c(objective = fit$objective, theta1 = fit$theta[[1L]])
  }

  first <- callr::r_bg(worker, libpath = .libPaths())
  second <- callr::r_bg(worker, libpath = .libPaths())
  on.exit({
    if (first$is_alive()) first$kill()
    if (second$is_alive()) second$kill()
  }, add = TRUE)
  first$wait(timeout = 60000)
  second$wait(timeout = 60000)

  first_result <- first$get_result()
  second_result <- second$get_result()
  expect_true(all(is.finite(first_result)))
  expect_equal(first_result, second_result, tolerance = 1e-10)
})
