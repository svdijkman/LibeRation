
### Helper functions

mat_exp <- function(M) {
  if (length(M) == 1L) return(matrix(exp(M), 1, 1))
  e <- eigen(M)
  Re(e$vectors %*% diag(exp(e$values)) %*% solve(e$vectors))
}

propagate <- function(K, A, dt) {
  as.numeric(mat_exp(K * dt) %*% A)
}

propagate_infusion <- function(K, A, dt, rate) {
  n <- length(A)

  M <- rbind(
    cbind(K, rate),
    c(rep(0, n), 0)
  )

  as.numeric(mat_exp(M * dt) %*% c(A, 1))[1:n]
}

ss_bolus_or_oral <- function(K, D, tau) {
  Phi <- mat_exp(K * tau)
  solve(diag(nrow(K)) - Phi, Phi %*% D)
}

ss_infusion <- function(K, rate, dur, tau) {
  n <- nrow(K)

  A0 <- rep(0, n)

  A_after_inf <- propagate_infusion(K, A0, dur, rate)
  A_end_cycle <- propagate(K, A_after_inf, tau - dur)

  Phi <- mat_exp(K * tau)

  solve(diag(n) - Phi, A_end_cycle)
}


### One compartment bolus model SS initial conditions solver

K1 <- function(CL, Vc) {
  matrix(-CL / Vc, 1, 1)
}

ss_1_bolus <- function(dose, tau, times, CL, Vc) {
  K <- K1(CL, Vc)

  D <- c(dose)

  Apre <- ss_bolus_or_oral(K, D, tau)
  A <- Apre + D

  data.frame(
    time = times,
    conc = sapply(times, function(t) propagate(K, A, t)[1] / Vc)
  )
}

### One compartment IV model SS initial conditions solver (uses K1 function same as bolus)

ss_1_infusion <- function(dose, tau, dur, times, CL, Vc) {
  K <- K1(CL, Vc)

  rate <- c(dose / dur)

  Apre <- ss_infusion(K, rate, dur, tau)

  data.frame(
    time = times,
    conc = sapply(times, function(t) {
      if (t <= dur) {
        A <- propagate_infusion(K, Apre, t, rate)
      } else {
        A <- propagate_infusion(K, Apre, dur, rate)
        A <- propagate(K, A, t - dur)
      }

      A[1] / Vc
    })
  )
}

### One compartment oral model SS initial conditions solver

K1_oral <- function(CL, Vc, KA) {
  matrix(c(
    -KA,       0,
     KA, -CL / Vc
  ), 2, 2, byrow = TRUE)
}

ss_1_oral <- function(dose, tau, times, CL, Vc, KA, F = 1) {
  K <- K1_oral(CL, Vc, KA)

  D <- c(F * dose, 0)  # depot, central

  Apre <- ss_bolus_or_oral(K, D, tau)
  A <- Apre + D

  data.frame(
    time = times,
    conc = sapply(times, function(t) propagate(K, A, t)[2] / Vc)
  )
}


### Two compartment bolus model SS initial conditions solver

K2 <- function(CL, Vc, Q, Vp) {
  k10 <- CL / Vc
  k12 <- Q / Vc
  k21 <- Q / Vp

  matrix(c(
    -(k10 + k12),  k21,
      k12,        -k21
  ), 2, 2, byrow = TRUE)
}

ss_2_bolus <- function(dose, tau, times, CL, Vc, Q, Vp) {
  K <- K2(CL, Vc, Q, Vp)

  D <- c(dose, 0)

  Apre <- ss_bolus_or_oral(K, D, tau)
  A <- Apre + D

  data.frame(
    time = times,
    conc = sapply(times, function(t) propagate(K, A, t)[1] / Vc)
  )
}

### Two compartment IV model SS initial conditions solver (uses K2 function same as bolus)

ss_2_infusion <- function(dose, tau, dur, times, CL, Vc, Q, Vp) {
  K <- K2(CL, Vc, Q, Vp)

  rate <- c(dose / dur, 0)

  Apre <- ss_infusion(K, rate, dur, tau)

  data.frame(
    time = times,
    conc = sapply(times, function(t) {
      if (t <= dur) {
        A <- propagate_infusion(K, Apre, t, rate)
      } else {
        A <- propagate_infusion(K, Apre, dur, rate)
        A <- propagate(K, A, t - dur)
      }

      A[1] / Vc
    })
  )
}

### Two compartment oral model SS initial conditions solver

K2_oral <- function(CL, Vc, Q, Vp, KA) {
  k10 <- CL / Vc
  k12 <- Q / Vc
  k21 <- Q / Vp

  matrix(c(
    -KA,       0,       0,
     KA, -(k10+k12),  k21,
      0,      k12,   -k21
  ), 3, 3, byrow = TRUE)
}

ss_2_oral <- function(dose, tau, times, CL, Vc, Q, Vp, KA, F = 1) {
  K <- K2_oral(CL, Vc, Q, Vp, KA)

  D <- c(F * dose, 0, 0)  # depot, central, peripheral

  Apre <- ss_bolus_or_oral(K, D, tau)
  A <- Apre + D

  data.frame(
    time = times,
    conc = sapply(times, function(t) propagate(K, A, t)[2] / Vc)
  )
}

### Three compartment bolus model SS initial conditions solver

K3 <- function(CL, Vc, Q2, Vp2, Q3, Vp3) {
  k10 <- CL / Vc
  k12 <- Q2 / Vc
  k21 <- Q2 / Vp2
  k13 <- Q3 / Vc
  k31 <- Q3 / Vp3

  matrix(c(
    -(k10+k12+k13),  k21,   k31,
       k12,         -k21,     0,
       k13,            0,  -k31
  ), 3, 3, byrow = TRUE)
}

ss_3_bolus <- function(dose, tau, times, CL, Vc, Q2, Vp2, Q3, Vp3) {
  K <- K3(CL, Vc, Q2, Vp2, Q3, Vp3)

  D <- c(dose, 0, 0)

  Apre <- ss_bolus_or_oral(K, D, tau)
  A <- Apre + D

  data.frame(
    time = times,
    conc = sapply(times, function(t) propagate(K, A, t)[1] / Vc)
  )
}


### Three compartment IV model SS initial conditions solver (uses K3 function same as bolus)

ss_3_infusion <- function(dose, tau, dur, times, CL, Vc, Q2, Vp2, Q3, Vp3) {
  K <- K3(CL, Vc, Q2, Vp2, Q3, Vp3)

  rate <- c(dose / dur, 0, 0)

  Apre <- ss_infusion(K, rate, dur, tau)

  data.frame(
    time = times,
    conc = sapply(times, function(t) {
      if (t <= dur) {
        A <- propagate_infusion(K, Apre, t, rate)
      } else {
        A <- propagate_infusion(K, Apre, dur, rate)
        A <- propagate(K, A, t - dur)
      }

      A[1] / Vc
    })
  )
}

### Three compartment oral model SS initial conditions solver


K3_oral <- function(CL, Vc, Q2, Vp2, Q3, Vp3, KA) {
  k10 <- CL / Vc
  k12 <- Q2 / Vc
  k21 <- Q2 / Vp2
  k13 <- Q3 / Vc
  k31 <- Q3 / Vp3

  matrix(c(
    -KA,              0,     0,     0,
     KA, -(k10+k12+k13),  k21,   k31,
      0,            k12,  -k21,    0,
      0,            k13,     0, -k31
  ), 4, 4, byrow = TRUE)
}

ss_3_oral <- function(dose, tau, times, CL, Vc, Q2, Vp2, Q3, Vp3, KA, F = 1) {
  K <- K3_oral(CL, Vc, Q2, Vp2, Q3, Vp3, KA)

  D <- c(F * dose, 0, 0, 0)  # depot, central, peripheral 1, peripheral 2

  Apre <- ss_bolus_or_oral(K, D, tau)
  A <- Apre + D

  data.frame(
    time = times,
    conc = sapply(times, function(t) propagate(K, A, t)[2] / Vc)
  )
}
