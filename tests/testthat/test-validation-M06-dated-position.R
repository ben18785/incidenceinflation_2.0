skip_validation("medium")

# M6 (medium tier): Geweke's test of the changes to the fatal cases dated by death
# (relabel, redate, unseen, inflow: deaths_dated_step) where they are recorded for a day
# part of the way to death (S5 of the sensitivity analyses, C5 of its test plan).
#
# Model (sensitivity_extensions.tex, S5; appendix, 'Cases recorded for their day of
# death'), in the exact form: 8 onset days; I(t) Poisson renewal with R = 1.1 from a
# history of 4 a day (serial gamma(5, 3), 10 days); F(t) ~ Bin(I(t), 0.5); each fatal
# case misdated with probability f = 0.4 and then recorded for day t + z1, z1 ~ p1, or
# unseen if that is after day 8; those with onset before day 1 recorded for day d,
# G_0(d) ~ Poisson(0.5 f sum_s history_s p1(d - s)). Time from onset to death gamma(4,
# 2): p_die on (z - 1, z] for the correctly dated, p1 = gamma(2, 1) on (z1 - 1, z1] and
# p2 = gamma(2, 1) on (z2, z2 + 1] for the misdated (ref_position_delays, a1 = a / 2).
# Reports of day d on days 2, 3, 5, 7, 8, 10, 12 (those on or after d): a non-fatal or
# correctly dated fatal case is reported after a gamma(3, 2) delay from d, a correctly
# dated one dies z ~ p_die days after d and is shown as dead once reported and its death
# recorded after a gamma(2, 1.5) delay; a misdated case dies z2 ~ p2 days after d and
# enters the records, shown as dead, when its death is recorded after the same delay.
# Many short successive-conditional chains: each starts from a draw of the parameters
# (the true cases and the labels) and the reports from the model, then alternates the
# package's update of the labels (deaths_dated_step, with the reports' probability from
# its exact forward filter with the survival of eq. (s5-mixed)) with drawing the reports
# again given the parameters, here, from the model. Each pair stays a draw from the
# joint distribution, so after the last step the parameters are draws from their prior.
# Checks: for each of several functions of the parameters, the mean over 2,000 chains
# after 10 steps against the mean over independent draws from the prior (two-sample z
# <= 4); and the update moves the labels of most chains.

m6_n_days <- 8
m6_report_days <- c(2, 3, 5, 7, 8, 10, 12)
m6_history <- rep(4, 10)
m6_omega <- ref_serial_weights(5, 3, 10)
m6_p_die <- ref_pdie(4, 2, max(1, ceiling(stats::qgamma(0.999, 4, 1))))
m6_position <- ref_position_delays(4, 2)

# the true cases and the labels, from the model
m6_draw_parameters <- function() {
  n <- m6_n_days
  p1 <- m6_position$p_recorded
  I <- sim_renewal(rep(1.1, n), m6_history, m6_omega)
  fatal <- stats::rbinom(n, I, 0.5)
  lags <- matrix(0, n, length(p1))
  unseen <- rep(0, n)
  for (t in seq_len(n)) {
    z1 <- sample.int(length(p1), stats::rbinom(1, fatal[t], 0.4), replace = TRUE, prob = p1)
    for (z in z1) {
      if (t + z <= n) lags[t, z] <- lags[t, z] + 1 else unseen[t] <- unseen[t] + 1
    }
  }
  s <- seq_along(m6_history) - length(m6_history)          # ..., -1, 0
  weight <- vapply(seq_len(n), function(d) {
    lag <- d - s
    ok <- lag <= length(p1)
    sum(m6_history[ok] * p1[lag[ok]])
  }, 0)
  list(I = I, fatal = fatal, lags = lags, unseen = unseen,
       inflow = stats::rpois(n, weight * 0.5 * 0.4), inflow_weight = weight)
}

# the misdated cases recorded for each day
m6_into <- function(lags, inflow) {
  n <- m6_n_days
  vapply(seq_len(n), function(d) {
    z <- seq_len(d - 1)
    sum(lags[cbind(d - z, z)]) + inflow[d]
  }, 0)
}

# the reports of each recorded day given the cases of each kind recorded for it
m6_draw_reports <- function(survivors, correct, misdated) {
  lapply(seq_len(m6_n_days), function(d) {
    rd <- m6_report_days[m6_report_days >= d]
    x <- rd - d
    t_survivor <- sim_rgamma_ms(survivors[d], 3, 2)
    t_correct <- sim_rgamma_ms(correct[d], 3, 2)
    z <- sample.int(length(m6_p_die), correct[d], replace = TRUE, prob = m6_p_die)
    rec_correct <- z + sim_rgamma_ms(correct[d], 2, 1.5)
    p2 <- m6_position$p_after
    rec_misdated <- sample.int(length(p2), misdated[d], replace = TRUE, prob = p2) - 1 +
      sim_rgamma_ms(misdated[d], 2, 1.5)
    cases <- vapply(x, function(v) sum(t_survivor <= v + 0.5) + sum(t_correct <= v + 0.5) +
                      sum(rec_misdated <= v + 0.5), 0)
    deaths <- vapply(x, function(v) sum(t_correct <= v + 0.5 & rec_correct <= v + 0.5) +
                       sum(rec_misdated <= v + 0.5), 0)
    tibble::tibble(time_reported = rd, cases_reported = cases, deaths_reported = deaths)
  })
}

# the state of the package's model of deaths holding the parameters and reports
m6_state <- function(par, reports, settings) {
  n <- m6_n_days
  out <- rowSums(par$lags) + par$unseen
  into <- m6_into(par$lags, par$inflow)
  p1 <- m6_position$p_recorded
  survivals <- lapply(seq_len(n), function(d) joint_survivals(
    reports[[d]]$time_reported, d, list(mean = 3, sd = 2), list(mean = 2, sd = 1.5), NULL,
    "recording", list(before = list(lag_probabilities = m6_p_die)),
    after_probabilities = settings$position$after_probabilities))
  list(onset_times = seq_len(n), deaths_true = par$fatal, cfr = rep(0.5, n),
       reports = reports, survivals = survivals, paths = vector("list", n),
       wrong = vector("list", n),
       dated = c(settings, list(
         lags = par$lags, unseen = par$unseen, inflow = par$inflow,
         inflow_weight = par$inflow_weight,
         tail = vapply(seq_len(n), function(t) sum(p1[seq_along(p1) > n - t]), 0),
         out = out, into = into)))
}

# the parameters of a state, as a list in the form drawn
m6_parameters <- function(state, cases_true) {
  dated <- state$dated
  list(I = cases_true, fatal = state$deaths_true, lags = dated$lags, unseen = dated$unseen,
       inflow = dated$inflow, inflow_weight = dated$inflow_weight)
}

# the reports drawn again given the parameters
m6_redraw <- function(state, cases_true) {
  dated <- state$dated
  correct <- state$deaths_true - dated$out
  state$reports <- m6_draw_reports(cases_true - state$deaths_true, correct, dated$into)
  state
}

# what is compared: functions of the parameters
m6_statistics <- function(par) {
  z1 <- col(par$lags)
  c("misdated recorded for a day analysed" = sum(par$lags),
    "sum of z1 of those" = sum(z1 * par$lags),
    "unseen" = sum(par$unseen),
    "from before day 1" = sum(par$inflow),
    "misdated recorded for days 6-8" = sum(m6_into(par$lags, par$inflow)[6:8]),
    "I(1)" = par$I[1], "I(4)" = par$I[4], "I(8)" = par$I[8],
    "F(2)" = par$fatal[2], "F(7)" = par$fatal[7],
    "dated by onset" = sum(par$fatal) - sum(par$lags) - sum(par$unseen))
}

test_that("M6 the changes to the fatal cases dated by death pass Geweke's test with the day recorded for part of the way to death", {
  skip_validation("medium")
  set.seed(7601)
  settings <- dated_by_death_settings(
    m6_n_days, list(fraction = 0.4, exact = TRUE, position = TRUE, n_moves = 4 * m6_n_days),
    onset_to_death = list(mean = 4, sd = 2))
  # the package's delays are those of the model
  expect_equal(settings$lag_probabilities, m6_position$p_recorded, tolerance = 1e-12)
  expect_equal(settings$position$after_probabilities, m6_position$p_after, tolerance = 1e-12)
  renewal <- dated_renewal(rep(1.1, m6_n_days), 1, FALSE, list(mean = 5, sd = 3), 10,
                           m6_history)
  n_chains <- 2000; n_steps <- 10
  t0 <- proc.time()[["elapsed"]]
  moved <- 0
  after <- t(vapply(seq_len(n_chains), function(i) {
    par <- m6_draw_parameters()
    correct <- par$fatal - rowSums(par$lags) - par$unseen
    state <- m6_state(par, m6_draw_reports(par$I - par$fatal, correct,
                                           m6_into(par$lags, par$inflow)), settings)
    cases_true <- par$I
    start <- state$dated$lags
    for (step in seq_len(n_steps)) {
      stepped <- deaths_dated_step(state, cases_true, renewal)
      state <- m6_redraw(stepped$state, stepped$cases_true)
      cases_true <- stepped$cases_true
    }
    moved <<- moved + !identical(state$dated$lags, start)
    m6_statistics(m6_parameters(state, cases_true))
  }, numeric(11)))
  note("M6: %d chains of %d steps took %.1f s", n_chains, n_steps, proc.time()[["elapsed"]] - t0)
  note("M6: the labels moved in %d of %d chains", moved, n_chains)
  expect(moved > n_chains / 2, sprintf("M6: the labels moved in only %d of %d chains", moved, n_chains))

  # independent draws from the prior
  n_prior <- 40000
  prior <- t(vapply(seq_len(n_prior), function(i) m6_statistics(m6_draw_parameters()),
                    numeric(11)))
  for (j in seq_len(ncol(after))) {
    se <- sqrt(stats::var(after[, j]) / n_chains + stats::var(prior[, j]) / n_prior)
    expect_z(mean(after[, j]), mean(prior[, j]), se,
             paste("M6", colnames(after)[j], "after the update against the prior"))
  }
})
