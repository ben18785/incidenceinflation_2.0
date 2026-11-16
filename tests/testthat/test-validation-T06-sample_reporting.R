skip_validation("fast")

# T6. sample_reporting: Metropolis-Hastings for the mean and sd of the gamma
# reporting delay. Its stationary distribution should be the exact posterior
#   prod_t prod_j Binomial(a_j | N - A_{j-1}, 1 - S_j / S_{j-1})     (spec eq. 17, S_0 = 1)
#   x Gamma(mean | mean = mean_mu, sd = mean_sigma) x Gamma(sd | mean = sd_mu, sd = sd_sigma)
# with S_j = 1 - G(t_j - t + 1/2) (A.4). The prior arguments are read as the mean and
# sd of each gamma prior. The exact posterior is
# integrated on a 400 x 400 grid (no mass at its edges).
# Design: the shape/rate reading of the prior moves the posterior mean of
# the delay mean from 8.02 to 7.55 (posterior sd 0.34).

t6_prior <- list(mean_mu = 10, mean_sigma = 5, sd_mu = 5, sd_sigma = 3)

t6_data <- local({
  set.seed(606)
  onsets <- 1:5; Ntrue <- 40; rdays <- c(2, 4, 6, 9, 13, 18)
  sim_chain <- function(N, onset, rd, mean, sd) {
    surv <- function(d) ref_surv_case(d, mean, sd)
    u <- stats::runif(N); grid <- 0:100; Sg <- surv(grid)
    dl <- vapply(u, function(x) { w <- which(Sg < x); if (length(w)) grid[w[1]] else Inf }, 0)
    vapply(rd, function(r) sum(dl <= r - onset), 0)
  }
  do.call(rbind, lapply(onsets, function(t) {
    rd <- rdays[rdays >= t]
    data.frame(time_onset = t, time_reported = rd,
               cases_reported = sim_chain(Ntrue, t, rd, 8, 4), cases_true = Ntrue)
  }))
})

t6_loglik <- function(m, s) {
  sum(vapply(unique(t6_data$time_onset), function(t) {
    d <- t6_data[t6_data$time_onset == t, ]
    ref_chain_loglik(d$cases_reported, d$time_reported, t, d$cases_true[1],
                     function(x) ref_surv_case(x, m, s))
  }, 0))
}

t6_exact <- function(prior_fun) {
  mg <- seq(2, 25, length.out = 400); sg <- seq(0.5, 20, length.out = 400)
  lp <- outer(mg, sg, Vectorize(t6_loglik)) + prior_fun(mg, sg)
  w <- exp(lp - max(lp)); w <- w / sum(w)
  M <- mg[row(w)]; S <- sg[col(w)]
  c(mean_m = sum(w * M), var_m = sum(w * M^2) - sum(w * M)^2,
    mean_s = sum(w * S), var_s = sum(w * S^2) - sum(w * S)^2)
}

test_that("T6 sample_reporting targets the exact posterior of the delay (mean, sd)", {
  set.seed(6006)
  ex <- t6_exact(function(mg, sg) outer(ref_dgamma_ms(mg, 10, 5, log = TRUE),
                                         ref_dgamma_ms(sg, 5, 3, log = TRUE), "+"))
  alt <- t6_exact(function(mg, sg) outer(stats::dgamma(mg, 10, 5, log = TRUE),
                                          stats::dgamma(sg, 5, 3, log = TRUE), "+"))
  start <- tibble::tibble(reporting_piece_index = 1L, mean = 8, sd = 4)
  logp0 <- t6_loglik(8, 4) + ref_dgamma_ms(8, 10, 5, log = TRUE) + ref_dgamma_ms(4, 5, 3, log = TRUE)
  # the four documented columns: reporting_piece_index may be left out with one piece
  snap <- tibble::as_tibble(t6_data)
  t0 <- proc.time()[["elapsed"]]
  out <- sample_reporting(snap, start, logp0, t6_prior,
                          metropolis_parameters = list(mean_step = 0.25, sd_step = 0.1),
                          ndraws = 20000)
  note("T6: 20000 MH iterations took %.1f s", proc.time()[["elapsed"]] - t0)
  rp <- out$reporting_parameters
  m <- col_of(rp, "mean"); s <- col_of(rp, "sd")
  expect_equal(length(m), 20000)
  keep <- seq(10001, 20000)
  m <- m[keep]; s <- s[keep]
  note("T6: acceptance (share of changed iterations) %.3f", mean(diff(m) != 0))
  diag <- sprintf("(shape/rate reading would give mean %.4f, sd %.4f)", alt[["mean_m"]], alt[["mean_s"]])
  expect_mcmc_mean(m, ex[["mean_m"]], "T6 posterior mean of the delay mean", diag)
  expect_mcmc_mean(s, ex[["mean_s"]], "T6 posterior mean of the delay sd", diag)
  expect_mcmc_mean((m - ex[["mean_m"]])^2, ex[["var_m"]], "T6 posterior variance of the delay mean")
  expect_mcmc_mean((s - ex[["mean_s"]])^2, ex[["var_s"]], "T6 posterior variance of the delay sd")
})

# ---- T6b: onset days with a single report -----------------------------------------
# T6 cannot tell whether sample_reporting uses each onset day's first report. Here every
# onset day has exactly ONE report, so under the spec (A.1: totals at t0 are zero,
# S_0 = 1) the likelihood is prod_t Binomial(C_t | N, 1 - S_C(t | 10)). If the first
# report were ignored the likelihood would be flat and the posterior would be the prior.
# Design: spec posterior mean/sd of the delay mean 10.10/1.10
# and of the delay sd 7.70/1.71; prior-only 10.02/4.84 and 4.99/2.96. Edge mass of the
# spec posterior on the grid used here: 2e-6.
t6b_onsets <- 1:8; t6b_N <- 40; t6b_rday <- 10
t6b_cnt <- local({
  set.seed(616)
  vapply(t6b_onsets, function(t) stats::rbinom(1, t6b_N, ref_cdf_case(t6b_rday - t, 8, 4)), 0)
})
t6b_loglik <- function(m, s)
  sum(stats::dbinom(t6b_cnt, t6b_N, ref_cdf_case(t6b_rday - t6b_onsets, m, s), log = TRUE))

test_that("T6b sample_reporting uses an onset day's first (here only) report", {
  set.seed(6016)
  mg <- seq(2, 30, length.out = 400); sg <- seq(0.3, 20, length.out = 400)
  lp <- outer(mg, sg, Vectorize(t6b_loglik)) +
    outer(ref_dgamma_ms(mg, 10, 5, log = TRUE), ref_dgamma_ms(sg, 5, 3, log = TRUE), "+")
  w <- exp(lp - max(lp)); w <- w / sum(w)
  M <- mg[row(w)]; S <- sg[col(w)]
  ex <- c(mean_m = sum(w * M), var_m = sum(w * M^2) - sum(w * M)^2,
          mean_s = sum(w * S), var_s = sum(w * S^2) - sum(w * S)^2)
  snap <- tibble::tibble(time_onset = as.integer(t6b_onsets), time_reported = as.integer(t6b_rday),
                         cases_reported = as.integer(t6b_cnt), cases_true = as.integer(t6b_N),
                         reporting_piece_index = 1L)
  start <- tibble::tibble(reporting_piece_index = 1L, mean = 10, sd = 7)
  logp0 <- t6b_loglik(10, 7) + ref_dgamma_ms(10, 10, 5, log = TRUE) + ref_dgamma_ms(7, 5, 3, log = TRUE)
  t0 <- proc.time()[["elapsed"]]
  out <- sample_reporting(snap, start, logp0, t6_prior,
                          metropolis_parameters = list(mean_step = 0.5, sd_step = 0.3),
                          ndraws = 60000)
  note("T6b: 60000 MH iterations took %.1f s", proc.time()[["elapsed"]] - t0)
  rp <- out$reporting_parameters
  m <- col_of(rp, "mean")[10001:60000]; s <- col_of(rp, "sd")[10001:60000]
  note("T6b: acceptance %.3f; sample mean/sd of delay mean %.3f/%.3f, of delay sd %.3f/%.3f",
       mean(diff(m) != 0), mean(m), stats::sd(m), mean(s), stats::sd(s))
  diag <- "(prior only, i.e. first report ignored: delay mean 10.02 sd 4.84; delay sd mean 4.99)"
  expect_mcmc_mean(m, ex[["mean_m"]], "T6b posterior mean of the delay mean", diag)
  expect_mcmc_mean(s, ex[["mean_s"]], "T6b posterior mean of the delay sd", diag)
  expect_mcmc_mean((m - ex[["mean_m"]])^2, ex[["var_m"]], "T6b posterior variance of the delay mean", diag)
  expect_mcmc_mean((s - ex[["mean_s"]])^2, ex[["var_s"]], "T6b posterior variance of the delay sd", diag)
})
