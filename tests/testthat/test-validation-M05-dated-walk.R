skip_validation("medium")

# M5 (medium tier): Geweke's test of the update of a fraction dated by death which
# varies by onset week (S9 of the sensitivity analyses, C5 of its test plan).
#
# Model (sensitivity_extensions.tex, S9): onset weeks of 7 days from day 1; logit f_(1)
# ~ N(0, 1.6^2), logit f_(w+1) = logit f_(w) + N(0, tau_f^2), tau_f ~ half-normal(0.3).
# Given the fractions, the data of the update are the labels: of the F(t) fatal cases of
# onset day t, F_mis(t) ~ Bin(F(t), f_t) are dated by death, f_t the fraction of the week
# of t; and the fatal cases from before day 1, G_0 ~ Poisson(r f_(1)).
# Successive-conditional chain: the package's update of the walk and tau_f
# (deaths_dated_fraction_step), then the labels drawn again given the fractions, here,
# from the model. The parameters' marginal is then their prior, whose moments are exact:
# E logit f_(w) = 0, Var logit f_(w) = 1.6^2 + (w - 1) 0.3^2, E tau_f = 0.3 sqrt(2 / pi),
# E tau_f^2 = 0.09; and the chain is also compared, for each, with independent draws
# from the prior (marginal-conditional), by a two-sample z.
# 15 onset days (three weeks, the last of one day), 4 fatal cases a day, r = 3.
# Checks: batch-means z <= 4 on the chain after a burn-in of 1,000.

m5_state <- function(n_fatal) {
  n_days <- length(n_fatal)
  settings <- dated_by_death_settings(n_days, list(fraction = list(walk = TRUE),
                                                   mean = 5, sd = 3))
  list(onset_times = seq_len(n_days), deaths_true = n_fatal, cfr = rep(0.5, n_days),
       dated = c(settings, list(out = rep(0, n_days), inflow = rep(0, n_days),
                                inflow_weight = c(6, rep(0, n_days - 1)))))
}

test_that("M5 the walk of the fraction dated by death passes Geweke's test", {
  skip_validation("medium")
  set.seed(7501)
  n_fatal <- rep(4, 15)
  week <- (seq_along(n_fatal) - 1) %/% 7 + 1
  state <- m5_state(n_fatal)
  r <- 6 * 0.5
  niter <- 400000; burn <- 1000
  out <- matrix(NA, niter, 4)
  t0 <- proc.time()[["elapsed"]]
  for (i in seq_len(niter)) {
    state <- deaths_dated_fraction_step(state)
    x <- state$dated$walk$logit_fraction
    f <- stats::plogis(x[week])
    # the labels given the fractions
    state$dated$out <- stats::rbinom(length(n_fatal), n_fatal, f)
    state$dated$inflow[1] <- stats::rpois(1, r * f[1])
    out[i, ] <- c(x, state$dated$walk$step_sd)
  }
  note("M5: %d iterations took %.1f s", niter, proc.time()[["elapsed"]] - t0)
  out <- out[-seq_len(burn), ]

  # independent draws from the prior
  n_prior <- 100000
  tau <- abs(stats::rnorm(n_prior, 0, 0.3))
  x1 <- stats::rnorm(n_prior, 0, 1.6)
  x2 <- x1 + tau * stats::rnorm(n_prior)
  prior <- cbind(x1, x2, x2 + tau * stats::rnorm(n_prior), tau)

  labels <- c("logit f_(1)", "logit f_(2)", "logit f_(3)", "tau_f")
  exact_mean <- c(0, 0, 0, 0.3 * sqrt(2 / pi))
  exact_square <- c(1.6^2, 1.6^2 + 0.09, 1.6^2 + 2 * 0.09, 0.09)
  for (j in 1:4) {
    expect_mcmc_mean(out[, j], exact_mean[j], paste("M5 mean of", labels[j]))
    expect_mcmc_mean(out[, j]^2, exact_square[j], paste("M5 mean square of", labels[j]))
    se <- sqrt(stat_batch_se(out[, j])^2 + stats::var(prior[, j]) / n_prior)
    expect_z(mean(out[, j]), mean(prior[, j]), se,
             paste("M5", labels[j], "against independent draws from the prior"))
  }
  # the steps standardised by tau_f are standard normal
  steps <- (out[, 2:3] - out[, 1:2]) / out[, 4]
  expect_mcmc_mean(steps[, 1]^2, 1, "M5 standardised first step")
  expect_mcmc_mean(steps[, 2]^2, 1, "M5 standardised second step")
  expect_mcmc_mean(as.numeric(out[, 4] < 0.3 * stats::qnorm(0.75)), 0.5,
                   "M5 tau_f below its prior median")
})
