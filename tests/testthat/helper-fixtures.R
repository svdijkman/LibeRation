estimation_fixture <- function(fix = TRUE) {
  model <- nm_model(
    INPUT = c("ID", "TIME", "EVID", "AMT", "DV", "MDV"),
    ADVAN = 1, DOSECMP = 1, OBSCMP = 1,
    PRED = "CL=THETA(1)*exp(ETA(1)); V=THETA(2); S1=V",
    ERROR = "Y=F+ERR(1)",
    THETAS = data.frame(THETA = 1:2, Value = c(2, 20), FIX = fix),
    OMEGAS = data.frame(OMEGA = 1, Value = 0.09, FIX = fix),
    SIGMAS = data.frame(SIGMA = 1, Value = 0.2, FIX = fix)
  )
  times <- c(0, 1, 4, 12)
  eta <- c(-0.15, 0.05, 0.2)
  data <- do.call(rbind, lapply(seq_along(eta), function(id) {
    prediction <- 5 * exp(-2 * exp(eta[[id]]) / 20 * times)
    data.frame(
      ID = id, TIME = times, EVID = c(1, 0, 0, 0),
      AMT = c(100, 0, 0, 0), MDV = c(1, 0, 0, 0),
      DV = c(NA, prediction[-1] + c(0.05, -0.03, 0.02))
    )
  }))
  list(model = model, data = data)
}
