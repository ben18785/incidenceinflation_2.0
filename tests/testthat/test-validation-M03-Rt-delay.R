skip_validation("medium")

# M3 (medium tier): mcmc_single with true cases held fixed
# (fixed_cases = TRUE), Poisson renewal, R_t changepoints sampled (priors$Rt_changepoints
# = list(), so the defaults apply: change prior Beta(1, T) and min_length 7, per the
# mcmc_single help page), cases only, no revisions.
# Expected values:
# - R_t: pi(P | I) propto B(1 + K, T + 0 - 1 - K) prod_seg M(u, v), with the closed-form
#   Poisson-gamma segment marginal M (A.10.5, eqs. 11 and 28), by dynamic programming;
#   E[R_t] = sum over segments of P(segment) (shape + I_uv) / (rate + Lambda_uv).
# - Delay (mean, sd): given I the delay is independent of R_t; its posterior is the
#   binomial chain (eq. 17, S_0 = 1) times the gamma priors read as (mean, sd), on a
#   120 x 120 grid (edge mass 4e-95).
# Data from the independent simulator: R = 1.3 on days 1-20 and 0.8 on days 21-40, history of 20
# a day, delay gamma(6, 3), reports every 2 days to day 40.
# Design: P(K = 1) = 0.948, P(K = 2) = 0.051; E[R_1] = 1.272,
# E[R_22] = 0.733; delay posterior mean 6.010 (sd 0.106), delay sd 2.987 (sd 0.094).

m_omega <- ref_serial_weights(15.3, 9.3, 40)

test_that("M3 mcmc_single with fixed cases: R_t changepoints and delay posteriors", {
  skip_validation("medium")
  set.seed(7003)
  hist <- rep(20, 40)
  I <- sim_renewal(c(rep(1.3, 20), rep(0.8, 20)), hist, m_omega)
  d <- sim_reports(I, seq(2, 40, 2), list(mean = 6, sd = 3))
  L <- ref_lambda(I, hist, m_omega)
  dp <- ref_partition_dp(40, ref_rt_log_seg_pois(I, L, 1, 1),
                         function(u, v) (1 + sum(I[u:v])) / (1 + sum(L[u:v])), a = 1, b = 40, min_length = 7)
  ll <- function(m, s) sum(vapply(1:40, function(t) {
    x <- d$cases[d$cases$time_onset == t, ]
    if (!nrow(x)) return(0)
    ref_chain_loglik(x$cases_reported, x$time_reported, t, I[t], function(z) ref_surv_case(z, m, s))
  }, 0))
  mg <- seq(3, 10, length.out = 120); sg <- seq(0.5, 7, length.out = 120)
  lp <- outer(mg, sg, Vectorize(ll)) + outer(ref_dgamma_ms(mg, 8, 4, log = TRUE), ref_dgamma_ms(sg, 4, 2, log = TRUE), "+")
  w <- exp(lp - max(lp)); w <- w / sum(w); M <- mg[row(w)]; S <- sg[col(w)]
  ex <- c(mean_m = sum(w * M), var_m = sum(w * M^2) - sum(w * M)^2,
          mean_s = sum(w * S), var_s = sum(w * S^2) - sum(w * S)^2)

  snap <- d$cases; snap$Rt_index <- 1L
  niter <- 4000
  set.seed(7103)
  t0 <- proc.time()[["elapsed"]]
  out <- mcmc_single(niter, snap,
                     priors = list(Rt = list(shape = 1, rate = 1),
                                   reporting = list(mean_mu = 8, mean_sigma = 4, sd_mu = 4, sd_sigma = 2),
                                   max_cases = 300, Rt_changepoints = list()),
                     serial_parameters = m_serial,
                     initial_cases_true = tibble::tibble(time_onset = 1:40, cases_true = as.integer(I)),
                     initial_reporting_parameters = list(mean = 6, sd = 3),
                     initial_Rt = tibble::tibble(Rt_index = 1L, Rt = 1),
                     print_to_screen = FALSE, cases_before = hist, fixed_cases = TRUE)
  note("M3: %d iterations took %.1f s", niter, proc.time()[["elapsed"]] - t0)

  K <- m_kept(out$Rt_changepoints, "n_changes", niter)
  note("M3: exact P(K = 0..3) %s; sample %s", paste(sprintf("%.4f", dp$pK[1:4]), collapse = " "),
       paste(sprintf("%.4f", vapply(0:3, function(k) mean(K == k), 0)), collapse = " "))
  for (k in 1:2) expect_mcmc_mean(as.numeric(K == k), dp$pK[k + 1], sprintf("M3 P(K = %d)", k))
  for (t in c(1, 15, 20, 21, 22, 30, 40))
    expect_mcmc_mean(m_series(out$Rt, t, "Rt", niter), dp$mean_by_day[t], sprintf("M3 E[R_%d]", t))

  rp <- out$reporting
  m <- m_kept(rp, "mean", niter); s <- m_kept(rp, "sd", niter)
  note("M3: delay mean %.4f (exact %.4f), delay sd %.4f (exact %.4f); acceptance %.3f",
       mean(m), ex[["mean_m"]], mean(s), ex[["mean_s"]], mean(diff(m) != 0))
  expect_mcmc_mean(m, ex[["mean_m"]], "M3 posterior mean of the delay mean")
  expect_mcmc_mean(s, ex[["mean_s"]], "M3 posterior mean of the delay sd")
  expect_mcmc_mean((m - ex[["mean_m"]])^2, ex[["var_m"]], "M3 posterior variance of the delay mean")
  expect_mcmc_mean((s - ex[["mean_s"]])^2, ex[["var_s"]], "M3 posterior variance of the delay sd")
})
