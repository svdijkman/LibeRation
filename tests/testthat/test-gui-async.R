test_that("simulation background entry point returns the engine result", {
  fixture <- estimation_fixture()
  result <- .liber_gui_background_task(
    "simulate",
    list(
      model = fixture$model,
      data = fixture$data,
      args = list(
        nsim = 1L, random_effects = FALSE, residual = FALSE,
        seed = 17L, n_cores = 1L
      )
    )
  )
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), nrow(fixture$data))
  expect_true(all(c("ID", "TIME", "DV") %in% names(result)))
})
