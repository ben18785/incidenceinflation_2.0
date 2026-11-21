skip_validation("medium")

# M6 (medium tier), C5 of S8 in sensitivity_extensions.tex: Geweke's test of the update
# of the share of the tail of the reporting delay for each onset week (the walk, the sd
# of its steps and the shared mean of the tail; deaths_tail_walk_step).
#
# Model: 15 onset days (onset weeks 1-7, 8-14, 15), 30 cases each, a known CFR of 0.4,
# fatal and other cases with the same reporting delay, no wrong records and no fatal
# cases reported when their death is recorded. The case delay is gamma(5, 3) with
# probability 1 - q_t and exponential with mean mu otherwise (eq. S(t | t') of S8);
# a case of onset day t is in the report of day r if its delay is at most r - t + 1/2.
# Prior: logit q_(1) ~ N(logit 0.036, 1.2^2), logit q_(w+1) - logit q_(w) ~ N(0, tau^2),
# tau ~ half-normal(0.15), mu ~ gamma with mean 20 and sd 8.
# Successive-conditional chain: simulate every case's reporting delay given the current
# (q, mu) with the independent simulator (sim_rdelay_tail), set the counts newly reported
# in each report, then one update of (q, tau, mu). The chain's marginal must be the prior.
# Checks (batch-means z, Bonferroni over the block): the means of logit q_(1),
# logit q_(3), tau and mu; E[(logit q_(1) - m)^2] = 1.2^2; E[(logit q_(3) - m)^2] =
# 1.2^2 + 2 E[tau^2]; E[tau^2] = 0.15^2; the variance of mu.

m05_state <- function(report_days, n_days, cases, delays) {
  snapshot <- expand.grid(time_onset = seq_len(n_days), time_reported = report_days)
  snapshot <- snapshot[snapshot$time_reported >= snapshot$time_onset, ]
  snapshot$cases_reported <- 0L
  snapshot$deaths_reported <- 0L
  settings <- deaths_settings(n_days, list(initial_cfr = 0.4,
                                           initial_delay = list(mean = 5, sd = 3)))
  state <- initial_deaths_state(tibble::as_tibble(snapshot),
                                tibble::tibble(time_onset = seq_len(n_days),
                                               cases_true = cases), settings)
  state$deaths_true <- round(0.4 * cases)
  state
}

# The split of the reports of each onset day, from delays drawn by the simulator: fatal
# and other cases newly reported in each report, and the deaths shown (death recorded
# after a gamma(5, 3) delay from onset, independent of the report; not informative here).
m05_paths <- function(state, cases, delay, share, tail_mean) {
  for (t in seq_along(state$onset_times)) {
    rd <- state$reports[[t]]$time_reported
    n <- cases[t]
    n_fatal <- state$deaths_true[t]
    is_f <- seq_len(n) <= n_fatal
    tc <- sim_rdelay_tail(n, delay, share[t], tail_mean)
    td <- sim_rgamma_ms(n, 5, 3)
    rep_f <- vapply(rd, function(r) sum(is_f & tc <= r - t + 0.5), 0)
    rep_s <- vapply(rd, function(r) sum(!is_f & tc <= r - t + 0.5), 0)
    dead <- vapply(rd, function(r) sum(is_f & tc <= r - t + 0.5 & td <= r - t + 0.5), 0)
    before <- function(x) c(0, x[-length(x)])
    state$paths[[t]] <- list(fatal_new = rep_f - before(rep_f), survivors_new = rep_s - before(rep_s),
                             fatal_before = before(rep_f), survivors_before = before(rep_s),
                             alive_before = before(rep_f) - before(dead),
                             deaths_new = dead - before(dead))
    state$reports[[t]]$cases_reported <- rep_f + rep_s
    state$reports[[t]]$deaths_reported <- dead
  }
  state
}

test_that("M6 (C5) Geweke: the update of the share of the tail for each week keeps the prior", {
  skip_validation("medium")
  set.seed(9501)
  n_days <- 15
  report_days <- c(3, 6, 10, 15, 22, 30, 45)
  cases <- rep(30, n_days)
  delay <- list(mean = 5, sd = 3)
  walk <- tail_walk_settings(TRUE)
  prior <- list(tail_mean_mu = 20, tail_mean_sigma = 8, tail_weight_walk = walk)
  settings <- list(tail_prior = prior, tail_steps = list(tail_mean_step = 0.4),
                   n_parameter_steps = 1)
  state <- m05_state(report_days, n_days, cases, delay)
  reporting_by_onset <- tibble::tibble(mean = rep(delay$mean, n_days), sd = rep(delay$sd, n_days))
  # a start drawn from the prior
  tau <- abs(stats::rnorm(1, 0, walk$step_sd_scale))
  x <- cumsum(c(stats::rnorm(1, walk$first_mean, walk$first_sd), stats::rnorm(2, 0, tau)))
  revision <- list(tail = list(weight = stats::plogis(x),
                               mean = stats::rgamma(1, shape = (20 / 8)^2, rate = 20 / 8^2),
                               day_first = 1, rw_sd = tau))
  state <- deaths_refresh_survivals(state, reporting_by_onset, revision)
  week <- (seq_len(n_days) - 1) %/% 7 + 1
  n_iter <- 20000
  draws <- matrix(NA, n_iter, 4, dimnames = list(NULL, c("x1", "x3", "tau", "mu")))
  t0 <- proc.time()[["elapsed"]]
  for (i in seq_len(n_iter)) {
    tail <- revision$tail
    state <- m05_paths(state, cases, delay, tail$weight[week], tail$mean)
    revision <- deaths_tail_walk_step(state, cases, reporting_by_onset, revision, settings)
    xi <- stats::qlogis(revision$tail$weight)
    draws[i, ] <- c(xi[1], xi[3], revision$tail$rw_sd, revision$tail$mean)
  }
  note("M6: %d Geweke iterations took %.1f min", n_iter, (proc.time()[["elapsed"]] - t0) / 60)
  m <- walk$first_mean
  s2 <- walk$first_sd^2
  t2 <- walk$step_sd_scale^2
  zmax <- z_bonferroni(9)
  f_check <- function(x, exact, label)
    expect_z(mean(x), exact, stat_batch_se(x), paste("M6", label), zmax = zmax)
  f_check(draws[, "x1"], m, "mean of logit q_(1)")
  f_check(draws[, "x3"], m, "mean of logit q_(3)")
  f_check(draws[, "tau"], walk$step_sd_scale * sqrt(2 / pi), "mean of tau")
  f_check(draws[, "mu"], 20, "mean of mu")
  f_check((draws[, "x1"] - m)^2, s2, "E[(logit q_(1) - m)^2]")
  f_check((draws[, "x3"] - m)^2, s2 + 2 * t2, "E[(logit q_(3) - m)^2]")
  f_check(draws[, "tau"]^2, t2, "E[tau^2]")
  f_check((draws[, "mu"] - 20)^2, 64, "variance of mu")
  f_check(draws[, "x3"] - draws[, "x1"], 0, "mean of logit q_(3) - logit q_(1)")
})
