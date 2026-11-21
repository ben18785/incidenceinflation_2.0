skip_validation("medium")

# M6 (medium tier): C5 of the sensitivity analyses for S7, a rate of wrong records which
# varies over the intervals between reports (sensitivity_extensions.tex, eq. s7-entered),
# by Geweke's successive-conditional simulator (m_wrong_geweke,
# helper-validation-medium.R).
#
# Model: 15 onset days with 10 true cases each, none fatal; reports on days 3, 5, 8, 12,
# 16, 20 and 26, so seven intervals; case delay gamma(3, 2). On day k, Poisson(eta_J(k) I
# [S(k-1-t) - S(k-t)]) wrong records are entered for onset day t, where J(k) is the first
# report on or after day k (helper-validation-sim.R, sim_wrong_records_day); each is
# withdrawn with probability rho on each later day and shown dead with probability
# delta. Priors: log eta_1 ~ N(log 0.04, 0.3^2), steps of log eta ~ N(0, tau^2 (days) / 7),
# tau ~ half-normal(0.25), rho ~ Beta(8, 12), delta ~ Beta(2, 3). At most 6 wrong records
# present: averaged over the prior, more are present at some report of some onset day
# with probability 4e-5 (wrong_cap_probability), so the truncated likelihood the sampler
# uses is the model's for this purpose.
# Checks: the chain's moments of log eta_1, of the change of log eta over the walk, of tau,
# rho and delta against the prior's, batch-means z <= 4.

m6_rd <- c(3L, 5L, 8L, 12L, 16L, 20L, 26L)

test_that("M6 (C5) Geweke: the updates of a walk of the rate of wrong records leave its prior unchanged", {
  skip_validation("medium")
  n_days <- 15
  I <- rep(10L, n_days)
  delay <- list(mean = 3, sd = 2)
  mu <- log(0.04); sigma <- 0.3; scale <- 0.25
  priors <- list(rho_alpha = 8, rho_beta = 12)
  settings <- deaths_settings(n_days, list(wrong_dead_alpha = 2, wrong_dead_beta = 3,
                                           n_parameter_steps = 2))
  f_simulate <- function(revision, delta) {
    rates <- exp(revision$wrong_walk$log_rate)
    sim_reports_deaths(I, m6_rd, delay, rep(0, n_days), ref_pdie(), rep(5, n_days),
                       rep(3, n_days), tau = 1, rho = revision$rho, delta = delta,
                       eta_day = function(k) rates[sim_interval_of_day(k, m6_rd)])
  }
  f_record <- function(revision, state) {
    walk <- revision$wrong_walk
    tibble::tibble(first = walk$log_rate[1], change = walk$log_rate[7] - walk$log_rate[1],
                   tau = walk$tau, rho = revision$rho, delta = state$delta)
  }
  niter <- 3000
  set.seed(9601)
  tau <- abs(stats::rnorm(1, 0, scale))
  log_rate <- mu + cumsum(c(stats::rnorm(1, 0, sigma),
                            stats::rnorm(6, 0, tau * sqrt(diff(m6_rd) / 7))))
  revision <- list(rho = stats::rbeta(1, 8, 12),
                   wrong_walk = list(log_rate = log_rate, report_days = m6_rd, tau = tau,
                                     mu = mu, sigma = sigma, tau_scale = scale,
                                     tau_step = 0.5))
  t0 <- proc.time()[["elapsed"]]
  draws <- m_wrong_geweke(niter, I, delay, revision, stats::rbeta(1, 2, 3), priors, settings,
                          f_simulate, f_record, max_erroneous = 6)
  note("M6: %d iterations took %.1f s", niter, proc.time()[["elapsed"]] - t0)
  expect_mcmc_mean(draws$first, mu, "M6 mean of log eta_1")
  expect_mcmc_mean((draws$first - mu)^2, sigma^2, "M6 variance of log eta_1")
  expect_mcmc_mean(draws$change^2, scale^2 * (26 - 3) / 7,
                   "M6 mean square of the change of log eta over the walk")
  expect_mcmc_mean(draws$tau, scale * sqrt(2 / pi), "M6 mean of tau")
  expect_mcmc_mean(draws$tau^2, scale^2, "M6 mean square of tau")
  f_beta_moments <- function(a, b) c(a / (a + b), a * (a + 1) / ((a + b) * (a + b + 1)))
  expect_mcmc_mean(draws$rho, f_beta_moments(8, 12)[1], "M6 mean of rho")
  expect_mcmc_mean(draws$rho^2, f_beta_moments(8, 12)[2], "M6 mean square of rho")
  expect_mcmc_mean(draws$delta, f_beta_moments(2, 3)[1], "M6 mean of delta")
})
