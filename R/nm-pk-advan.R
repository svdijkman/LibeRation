#' PKADVAN-style analytical step functions (R port of src/nm_pk_pkadvan.h)
#' @keywords internal

.ad_sqrt <- function(x) {
  if (.nm_is_ad(x)) {
    .ad_dispatch("sqrt", x)
  } else {
    sqrt(pmax(as.numeric(x), 0))
  }
}

#' @keywords internal
.nm_hybrid_lambdas_2 <- function(k10, k12, k21) {
  E1 <- .ad_add(k10, k12)
  E2 <- k21
  s <- .ad_add(E1, E2)
  disc <- .ad_add(
    .ad_sub(.ad_mul(s, s), .ad_mul(4, .ad_mul(E1, E2))),
    .ad_mul(4, .ad_mul(k12, k21))
  )
  root <- .ad_sqrt(disc)
  half <- 0.5
  list(
    lambda1 = .ad_mul(half, .ad_add(s, root)),
    lambda2 = .ad_mul(half, .ad_sub(s, root)),
    E1 = E1,
    E2 = E2
  )
}

#' @keywords internal
.nm_oral_gut_sum3 <- function(dt, ka, l1, l2, numer_fn) {
  t1 <- .ad_div(
    .ad_mul(.ad_exp(.ad_mul(-dt, ka)), numer_fn(ka)),
    .ad_mul(.ad_sub(l1, ka), .ad_sub(l2, ka))
  )
  t2 <- .ad_div(
    .ad_mul(.ad_exp(.ad_mul(-dt, l1)), numer_fn(l1)),
    .ad_mul(.ad_sub(l2, l1), .ad_sub(ka, l1))
  )
  t3 <- .ad_div(
    .ad_mul(.ad_exp(.ad_mul(-dt, l2)), numer_fn(l2)),
    .ad_mul(.ad_sub(l1, l2), .ad_sub(ka, l2))
  )
  .ad_add(.ad_add(t1, t2), t3)
}

#' @keywords internal
.nm_step_2_oral <- function(dt, ka, k20, k23, k32, prev_gut, prev_a2, prev_a3) {
  h <- .nm_hybrid_lambdas_2(k20, k23, k32)
  l1 <- h$lambda1
  l2 <- h$lambda2
  E2 <- h$E1
  E3 <- h$E2
  d <- .ad_sub(l2, l1)
  d_num <- .ad_scalar_value(d)
  if (abs(d_num) < 1e-12) {
    e <- .ad_exp(.ad_mul(-l1, dt))
    a2 <- .ad_mul(prev_a2, e)
    a3 <- .ad_mul(prev_a3, e)
  } else {
    B2 <- .ad_add(.ad_mul(prev_a2, E3), .ad_mul(prev_a3, k32))
    a2 <- .ad_div(
      .ad_sub(
        .ad_mul(.ad_sub(B2, .ad_mul(prev_a2, l1)), .ad_exp(.ad_mul(-l1, dt))),
        .ad_mul(.ad_sub(B2, .ad_mul(prev_a2, l2)), .ad_exp(.ad_mul(-l2, dt)))
      ),
      d
    )
    B3 <- .ad_add(.ad_mul(prev_a3, E2), .ad_mul(prev_a2, k23))
    a3 <- .ad_div(
      .ad_sub(
        .ad_mul(.ad_sub(B3, .ad_mul(prev_a3, l1)), .ad_exp(.ad_mul(-l1, dt))),
        .ad_mul(.ad_sub(B3, .ad_mul(prev_a3, l2)), .ad_exp(.ad_mul(-l2, dt)))
      ),
      d
    )
  }
  gut_cen <- .nm_oral_gut_sum3(dt, ka, l1, l2, function(x) .ad_sub(E3, x))
  gut_per <- .nm_oral_gut_sum3(dt, ka, l1, l2, function(x) 1)
  a2 <- .ad_add(a2, .ad_mul(.ad_mul(prev_gut, ka), gut_cen))
  a3 <- .ad_add(a3, .ad_mul(.ad_mul(.ad_mul(prev_gut, ka), k23), gut_per))
  gut <- .ad_mul(prev_gut, .ad_exp(.ad_mul(-dt, ka)))
  list(gut = gut, a2 = a2, a3 = a3)
}

#' @keywords internal
.nm_cmp2oral_trans4_state <- function(dat, ka, k20, k23, k32, ad_mode) {
  n <- nrow(dat)
  if (!ad_mode) {
    gut <- rep(0, n)
    a2 <- rep(0, n)
    a3 <- rep(0, n)
    for (i in seq_len(n)) {
      if (i > 1L) {
        dt <- dat$TIME[i] - dat$TIME[i - 1L]
        st <- .nm_step_2_oral(dt, ka, k20, k23, k32, gut[i - 1L], a2[i - 1L], a3[i - 1L])
        gut[i] <- st$gut
        a2[i] <- st$a2
        a3[i] <- st$a3
      }
      if (dat$EVID[i] == 1L) {
        gut[i] <- gut[i] + dat$AMT[i] * dat$F1[i]
      }
    }
    return(list(gut = gut, a2 = a2, a3 = a3))
  }
  gut <- vector("list", n)
  a2 <- vector("list", n)
  a3 <- vector("list", n)
  z <- newConstant(name = "pk_zero", value = 0)
  for (i in seq_len(n)) {
    dose <- if (dat$EVID[i] == 1L) {
      newConstant(name = paste0("dose_", i), value = dat$AMT[i] * dat$F1[i])
    } else {
      z
    }
    if (i == 1L) {
      gut[[i]] <- dose
      a2[[i]] <- z
      a3[[i]] <- z
    } else {
      dt <- dat$TIME[i] - dat$TIME[i - 1L]
      st <- .nm_step_2_oral(dt, ka, k20, k23, k32, gut[[i - 1L]], a2[[i - 1L]], a3[[i - 1L]])
      gut[[i]] <- .ad_add(st$gut, dose)
      a2[[i]] <- st$a2
      a3[[i]] <- st$a3
    }
  }
  list(gut = gut, a2 = a2, a3 = a3)
}
