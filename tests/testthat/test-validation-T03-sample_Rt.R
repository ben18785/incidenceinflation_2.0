skip_validation("fast")

# T3. sample_Rt under the Poisson renewal model: conjugate gamma posterior
#   R | I ~ Gamma(shape + sum_{t in piece} I_t, rate + sum_{t in piece} Lambda_t),
#   Lambda_t = sum_{s=1}^{40} omega_s I(t - s), omega_s on (s - 1/2, s + 1/2]  (spec A.2).
# A growing history makes the discretisation of omega matter: with this design the
# (s-1, s] rule shifts the posterior mean by z ~ 48 and an unnormalised omega by
# z ~ 31.
#
# T3b. Under the negative binomial model sample_Rt takes Metropolis-Hastings steps
# from Rt_current; the target is pi(R) proportional to
#   prod_t NegBin(I_t | R Lambda_t, kappa) Gamma(R | shape, rate),
# computed on a fine grid (helper-validation-stats.R ref_grid_dist).

t3_serial <- list(mean = 15.3, sd = 9.3)
t3_omega <- ref_serial_weights(15.3, 9.3, 40)
t3_hist <- round(20 * exp(0.06 * (1:60)))           # cases_before, earliest first
t3_cases <- local({                                   # data: roughly R = 1.9 throughout
  I <- numeric(10)
  for (t in 1:10) I[t] <- round(1.9 * ref_lambda(I, t3_hist, t3_omega)[t])
  I
})
t3_index <- rep(1:2, each = 5)
t3_prior <- list(shape = 2, rate = 0.5)

t3_piece <- function(out) col_of(out, "Rt_index")

test_that("T3 sample_Rt (Poisson) draws from the conjugate gamma posterior", {
  set.seed(3003)
  df <- tibble::tibble(time_onset = 1:10, cases_true = t3_cases, Rt_index = t3_index)
  t0 <- proc.time()[["elapsed"]]
  out <- sample_Rt(df, t3_prior, t3_serial, serial_max = 40, ndraws = 4000,
                   is_negative_binomial = FALSE, cases_before = t3_hist)
  note("T3: 4000 draws took %.1f s", proc.time()[["elapsed"]] - t0)
  L <- ref_lambda(t3_cases, t3_hist, t3_omega)
  # alternative discretisation, for diagnosis only
  sh <- ref_gamma_shape(15.3, 9.3); rt <- ref_gamma_rate(15.3, 9.3)
  wB <- stats::pgamma(1:40, sh, rt) - stats::pgamma(0:39, sh, rt); wB <- wB / sum(wB)
  LB <- ref_lambda(t3_cases, t3_hist, wB)
  for (k in 1:2) {
    x <- col_of(out, "Rt")[t3_piece(out) == k]
    expect_equal(length(x), 4000)
    a <- t3_prior$shape + sum(t3_cases[t3_index == k])
    b <- t3_prior$rate + sum(L[t3_index == k])
    bB <- t3_prior$rate + sum(LB[t3_index == k])
    expect_mean_z(x, a / b, sprintf("T3 piece %d posterior mean", k),
                  sprintf("((s-1,s] omega would give %.5f)", a / bB))
    expect_z(stats::sd(x), sqrt(a) / b, stats::sd(x) / sqrt(2 * (length(x) - 1)),
             sprintf("T3 piece %d posterior sd", k))
    expect_ks(x, function(q) stats::pgamma(q, a, b), sprintf("T3 piece %d distribution", k))
  }
})

test_that("T3 sample_Rt (Poisson, maximise = TRUE) returns the posterior mode", {
  df <- tibble::tibble(time_onset = 1:10, cases_true = t3_cases, Rt_index = t3_index)
  out <- sample_Rt(df, t3_prior, t3_serial, serial_max = 40, maximise = TRUE,
                   is_negative_binomial = FALSE, cases_before = t3_hist)
  L <- ref_lambda(t3_cases, t3_hist, t3_omega)
  for (k in 1:2) {
    x <- col_of(out, "Rt")[t3_piece(out) == k]
    a <- t3_prior$shape + sum(t3_cases[t3_index == k])
    b <- t3_prior$rate + sum(L[t3_index == k])
    expect_equal(unique(x), (a - 1) / b, tolerance = 1e-4,
                 label = sprintf("T3 piece %d maximised Rt", k))
  }
})

# ---- T3b: negative binomial, MH kernel ------------------------------------------
t3b_hist <- rep(30, 60)
t3b_cases <- c(41, 50, 47, 39, 52, 44, 49, 46, 43, 51)
t3b_kappa <- 5
t3b_prior <- list(shape = 2, rate = 1)
t3b_target <- local({
  L <- ref_lambda(t3b_cases, t3b_hist, t3_omega)
  x <- seq(1e-4, 8, length.out = 1e5)
  ld <- vapply(x, function(r) sum(ref_dnb(t3b_cases, r * L, t3b_kappa, log = TRUE)), 0) +
    stats::dgamma(x, t3b_prior$shape, t3b_prior$rate, log = TRUE)
  ref_grid_dist(x, ld)
})
t3b_call <- function(r0, nsteps = NULL) {
  df <- tibble::tibble(time_onset = 1:10, cases_true = t3b_cases, Rt_index = 1L)
  args <- list(df, t3b_prior, t3_serial, kappa = t3b_kappa, serial_max = 40, ndraws = 1,
               is_negative_binomial = TRUE,
               Rt_current = tibble::tibble(Rt_index = 1L, Rt = r0),
               cases_before = t3b_hist)
  if (!is.null(nsteps)) args$nsteps <- nsteps
  col_of(do.call(sample_Rt, args), "Rt")[1]
}

test_that("T3b sample_Rt (NB): invariance of the exact posterior, and the kernel moves", {
  set.seed(3103)
  N <- 3000
  r0 <- t3b_target$rand(N)
  t0 <- proc.time()[["elapsed"]]
  r1 <- vapply(r0, t3b_call, numeric(1))
  note("T3b(i): %d calls took %.1f s; posterior mean %.4f", N, proc.time()[["elapsed"]] - t0, t3b_target$mean)
  expect_ks(r1, t3b_target$cdf, "T3b(i) outputs from exact-posterior starts")
  expect_mean_z(r1, t3b_target$mean, "T3b(i) mean of outputs")
  moved <- mean(abs(r1 - r0) > 1e-12)
  note("T3b(ii): share of calls whose output differs from Rt_current: %.3f", moved)
  expect_gt(moved, 0)
})

test_that("T3b sample_Rt (NB): convergence from a poor start (weaker evidence alone)", {
  set.seed(3104)
  N <- 3000
  t0 <- proc.time()[["elapsed"]]
  r1 <- vapply(seq_len(N), function(i) t3b_call(0.2, nsteps = 200), numeric(1))
  note("T3b(iii): %d calls with nsteps = 200 took %.1f s", N, proc.time()[["elapsed"]] - t0)
  expect_ks(r1, t3b_target$cdf, "T3b(iii) outputs from R = 0.2 after 200 steps")
  expect_mean_z(r1, t3b_target$mean, "T3b(iii) mean of outputs")
})
