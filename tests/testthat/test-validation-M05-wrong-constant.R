skip_validation("medium")

# M5 (medium tier): C5 of the sensitivity analyses for S6, wrong records at a constant
# rate (sensitivity_extensions.tex, eq. s6-entered), by Geweke's successive-conditional
# simulator (m_wrong_geweke, helper-validation-medium.R).
#
# Model: 15 onset days with 8 true cases each, none fatal; reports on days 3, 5, 8, 12,
# 16, 20 and 26; case delay gamma(3, 2). On each day after an onset day, Poisson(eta_bg)
# wrong records are entered for it, each withdrawn with probability rho on each later
# day and shown dead with probability delta (helper-validation-sim.R,
# sim_wrong_records_constant). Priors: eta_bg ~ Gamma(shape 6.25, rate 62.5) (mean 0.1,
# sd 0.04), rho ~ Beta(8, 12), delta ~ Beta(2, 3). At most 6 wrong records present:
# averaged over the prior, more are present at some report of some onset day with
# probability 2e-5 (wrong_cap_probability), so the truncated likelihood the sampler
# uses is the model's for this purpose.
# Checks: the chain's mean and mean square of eta_bg, rho and delta against the prior's,
# batch-means z <= 4.

m5_rd <- c(3L, 5L, 8L, 12L, 16L, 20L, 26L)

test_that("M5 (C5) Geweke: the updates of a constant rate of wrong records leave its prior unchanged", {
  skip_validation("medium")
  n_days <- 15
  I <- rep(8L, n_days)
  delay <- list(mean = 3, sd = 2)
  priors <- list(eta_bg_mu = 0.1, eta_bg_sigma = 0.04, rho_alpha = 8, rho_beta = 12)
  settings <- deaths_settings(n_days, list(wrong_dead_alpha = 2, wrong_dead_beta = 3,
                                           n_parameter_steps = 2))
  f_simulate <- function(revision, delta)
    sim_reports_deaths(I, m5_rd, delay, rep(0, n_days), ref_pdie(), rep(5, n_days),
                       rep(3, n_days), tau = 1, rho = revision$rho, delta = delta,
                       eta_bg = revision$background_rate)
  f_record <- function(revision, state)
    tibble::tibble(eta_bg = revision$background_rate, rho = revision$rho, delta = state$delta)
  niter <- 3000
  set.seed(9501)
  revision <- list(rho = stats::rbeta(1, 8, 12), background_rate = stats::rgamma(1, 6.25, 62.5))
  t0 <- proc.time()[["elapsed"]]
  draws <- m_wrong_geweke(niter, I, delay, revision, stats::rbeta(1, 2, 3), priors, settings,
                          f_simulate, f_record, max_erroneous = 6)
  note("M5: %d iterations took %.1f s", niter, proc.time()[["elapsed"]] - t0)
  shape <- 6.25; rate <- 62.5
  expect_mcmc_mean(draws$eta_bg, shape / rate, "M5 mean of eta_bg")
  expect_mcmc_mean(draws$eta_bg^2, shape * (shape + 1) / rate^2, "M5 mean square of eta_bg")
  f_beta_moments <- function(a, b) c(a / (a + b), a * (a + 1) / ((a + b) * (a + b + 1)))
  expect_mcmc_mean(draws$rho, f_beta_moments(8, 12)[1], "M5 mean of rho")
  expect_mcmc_mean(draws$rho^2, f_beta_moments(8, 12)[2], "M5 mean square of rho")
  expect_mcmc_mean(draws$delta, f_beta_moments(2, 3)[1], "M5 mean of delta")
  expect_mcmc_mean(draws$delta^2, f_beta_moments(2, 3)[2], "M5 mean square of delta")
})
