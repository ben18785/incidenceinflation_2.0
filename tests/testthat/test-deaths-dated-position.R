# Fatal cases dated by death recorded for a day part of the way from onset to
# death (S5 of the sensitivity analyses): a misdated case with onset on day t
# is recorded for day t + z1 and dies z2 days after that day, z1 and z2
# independent, and is reported when its death is recorded. The labels hold
# z1; z2 is summed out of the survival of each misdated case.

# the shape and rate of the time from onset to death of the fits
position_shape <- 9.3^2 / 5.5^2
position_rate <- 9.3 / 5.5^2

# the two delays written from the note: z1 >= 1 gamma(a1, c) on (z1 - 1, z1]
# and z2 >= 0 gamma(a - a1, c) on (z2, z2 + 1], each truncated at its 99.9th
# percentile and rescaled
position_delays_reference <- function(a1, a=position_shape, c=position_rate) {
  z1_max <- max(1, ceiling(stats::qgamma(0.999, a1, c)))
  p1 <- stats::pgamma(seq_len(z1_max), a1, c) - stats::pgamma(seq_len(z1_max) - 1, a1, c)
  p2 <- 1
  if(a1 < a) {
    z2_max <- ceiling(stats::qgamma(0.999, a - a1, c))
    p2 <- stats::pgamma(seq_len(z2_max), a - a1, c) -
      stats::pgamma(seq_len(z2_max) - 1, a - a1, c)
  }
  list(p1=p1 / sum(p1), p2=p2 / sum(p2))
}

# the survival of the delay from a death on day d to its being recorded, by
# the end of day t, written from the appendix: the gamma delay, with the rate
# of recording on each day after the death multiplied by that day's
# multiplier (that of the first report on or after it: the effect of its day
# of the week times its own effect; one outside the reports), mixed with a
# share b of deaths in a backlog cleared with a constant hazard from the
# middle of the day hold days after the death
reference_survival_recorded <- function(t, d, mean, sd, multiplier=function(i) 1,
                                        share=0, rate=0, hold=0) {
  f_log <- function(x) stats::pgamma(x, shape=mean^2 / sd^2, rate=mean / sd^2,
                                     lower.tail=FALSE, log.p=TRUE)
  log_g <- f_log(0.5)
  for(i in seq_len(t - d) + d)
    log_g <- log_g + multiplier(i) * (f_log(i - d + 0.5) - f_log(i - d - 0.5))
  (1 - share) * exp(log_g) + share * exp(-rate * max(0, t - (d + hold) + 0.5))
}

# eq. (s5-mixed): a misdated case recorded for day d has not been reported by
# the end of day t, summed over the day of death d + z2, with the share of
# the backlog given (that of the day recorded for)
reference_survival_mixed <- function(t, d, p2, ...) {
  z2 <- seq_along(p2) - 1
  sum(p2[z2 > t - d]) +
    sum(vapply(z2[z2 <= t - d], function(z)
      p2[z + 1] * reference_survival_recorded(t, d + z, ...), numeric(1)))
}

# the effects of the day of the week and of each report on the delay to a
# death being recorded, and a backlog whose share differs by day, with a
# holding period
position_effects <- function(report_days, n_days=40) {
  effects <- update_report_multipliers(list(
    report_days=report_days, weekday_first=1,
    weekday_log_effects=stats::setNames(c(0.3, -0.2, 0, 0.5, -0.4, 0.1, 0.2), 1:7),
    report_log_effects=seq(-0.3, 0.4, length.out=length(report_days))))
  share <- 0.1 + 0.6 * (seq_len(n_days) %% 3 == 0) + 0.1 * (seq_len(n_days) %% 2 == 0)
  effects$backlog <- initial_backlog(n_days, report_days, backlog_settings(list(
    share_week=1, initial_share=share, initial_rate=0.15, hold_days=2)))
  effects$before <- list(lag_probabilities=lag_probabilities_gamma(9.3, 5.5))
  effects
}

# the multiplier of each day, from the effects as specified
position_multiplier <- function(effects) {
  report_days <- effects$report_days
  function(i) {
    if(i < min(report_days) | i > max(report_days))
      return(1)
    j <- which(report_days >= i)[1]
    exp(effects$weekday_log_effects[[as.character((report_days[j] - 1) %% 7 + 1)]] +
          effects$report_log_effects[j])
  }
}

test_that("the settings of the position of the day recorded for are checked", {
  onset_to_death <- list(mean=9.3, sd=5.5)
  f_settings <- function(...)
    deaths_settings(30, list(onset_to_death=onset_to_death,
                             dated_by_death=list(fraction=list(), ...)))$dated_by_death
  main <- f_settings()
  # without the setting, or with it FALSE, the settings are as they were
  expect_false("position" %in% names(main))
  expect_identical(f_settings(position=FALSE), main)
  expect_identical(f_settings(position=NULL), main)
  # by default the time to death is divided at a share with a Beta(a / 2, a / 2)
  # distribution
  dated <- f_settings(position=TRUE)
  reference <- position_delays_reference(position_shape / 2)
  expect_equal(dated$position$share_shape, position_shape / 2)
  expect_equal(dated$lag_probabilities, reference$p1, tolerance=1e-12)
  expect_equal(dated$position$after_probabilities, reference$p2, tolerance=1e-12)
  expect_equal(dated$max_lag, length(reference$p1))
  expect_false("lag_probabilities" %in% names(dated$position))
  expect_equal(f_settings(position=list(share_shape=1))$lag_probabilities,
               position_delays_reference(1)$p1, tolerance=1e-12)
  # a longest time to the day recorded for can be given, but that of the
  # time to death of other fatal cases is not taken for it
  expect_equal(deaths_settings(30, list(
    onset_to_death=list(mean=9.3, sd=5.5, max_lag=40),
    dated_by_death=list(fraction=list(), position=TRUE)))$dated_by_death$max_lag,
    length(reference$p1))
  expect_equal(f_settings(position=TRUE, max_lag=8)$max_lag, 8)
  # or the two distributions given
  given <- f_settings(position=list(lag_probabilities=c(0.5, 0.5),
                                    after_probabilities=c(0.2, 0.3, 0.5)))
  expect_equal(given$lag_probabilities, c(0.5, 0.5))
  expect_equal(given$max_lag, 2)
  expect_equal(given$position$after_probabilities, c(0.2, 0.3, 0.5))
  expect_null(given$position$share_shape)

  expect_error(f_settings(position="yes"), "TRUE or a named list")
  expect_error(f_settings(position=list(shape=1)), "Unknown settings for the position")
  expect_error(f_settings(position=list(share_shape=0)), "must lie in")
  expect_error(f_settings(position=list(share_shape=-1)), "must lie in")
  expect_error(f_settings(position=list(share_shape=position_shape * 1.01)), "must lie in")
  expect_error(f_settings(position=list(share_shape=c(1, 2))), "must lie in")
  expect_error(f_settings(position=list(share_shape=1, after_probabilities=1)),
               "not both")
  expect_error(f_settings(position=list(lag_probabilities=c(0.5, 0.5))),
               "both 'lag_probabilities' and 'after_probabilities'")
  expect_error(f_settings(position=list(lag_probabilities=c(0.5, 0.6),
                                        after_probabilities=1)), "summing to one")
  expect_error(f_settings(position=list(lag_probabilities=c(1.5, -0.5),
                                        after_probabilities=1)), "non-negative")
  # it needs the delay from a death to its being recorded
  expect_error(f_settings(position=TRUE, reporting=list(mean=2, sd=1)),
               "needs reporting \"recording\"")
  expect_error(dated_by_death_settings(30, list(fraction=0.3, mean=9.3, sd=5.5,
                                                position=TRUE)),
               "needs reporting \"recording\"")
  expect_error(relaxed_stan_data(
    dplyr::tibble(time_onset=1, time_reported=1, cases_reported=1), list(mean=5, sd=3),
    deaths_df=dplyr::tibble(time_onset=1, time_reported=1, deaths_reported=0),
    dated_by_death=list(fraction=0.3, mean=9.3, sd=5.5, position=TRUE)),
    "part of the way to death")
  expect_error(simulate_cases_deaths(rep(10, 5), 1:10, list(mean=5, sd=3), 0.5,
                                     list(mean=3, sd=2),
                                     dated_by_death=list(fraction=0.3, mean=9.3, sd=5.5,
                                                         position=TRUE)),
               "need onset_to_death")
})

test_that("the two delays sum to one, start at one and zero days and have the means of the discretised gamma distributions", {
  dated <- deaths_settings(30, list(onset_to_death=list(mean=9.3, sd=5.5),
                                    dated_by_death=list(fraction=0.3,
                                                        position=TRUE)))$dated_by_death
  p1 <- dated$lag_probabilities
  p2 <- dated$position$after_probabilities
  expect_equal(sum(p1), 1)
  expect_equal(sum(p2), 1)
  expect_true(all(p1 > 0) & all(p2 > 0))
  # z1 from one day, z2 from zero: the first of each is the probability of
  # the gamma's first day, (0, 1]
  a1 <- position_shape / 2
  expect_equal(p1[1] * stats::pgamma(length(p1), a1, position_rate),
               stats::pgamma(1, a1, position_rate))
  expect_equal(p2[1] * stats::pgamma(length(p2), position_shape - a1, position_rate),
               stats::pgamma(1, position_shape - a1, position_rate))
  # their means: the gamma rounded up and rounded down, by integration over
  # each day, within the 99.9th percentile
  f_mean <- function(shape, n, shift) {
    days <- seq_len(n)
    parts <- vapply(days, function(z) stats::integrate(
      function(x) stats::dgamma(x, shape, position_rate), z - 1, z)$value, numeric(1))
    sum((days - shift) * parts) / sum(parts)
  }
  mean_1 <- sum(seq_along(p1) * p1)
  mean_2 <- sum((seq_along(p2) - 1) * p2)
  expect_equal(mean_1, f_mean(a1, length(p1), 0), tolerance=1e-6)
  expect_equal(mean_2, f_mean(position_shape - a1, length(p2), 1), tolerance=1e-6)
  # each about the gamma's mean, 4.65 days, plus or less half a day
  expect_equal(mean_1, 9.3 / 2 + 0.5, tolerance=0.02)
  expect_equal(mean_2, 9.3 / 2 - 0.5, tolerance=0.02)
  # the time from onset to death of a misdated case, z1 + z2, has close to
  # the mean of the time from onset to death, 9.3 days, which the discretised
  # p_die exceeds by half a day
  convolution <- stats::convolve(p1, rev(p2), type="open")
  z <- seq_along(convolution)
  expect_equal(sum(convolution), 1)
  expect_lt(abs(sum(z * convolution) - 9.3), 0.1)
  expect_equal(sum(z * convolution), mean_1 + mean_2)
  p_die <- lag_probabilities_gamma(9.3, 5.5)
  expect_lt(abs(sum(seq_along(p_die) * p_die) - 9.3 - 0.5), 0.05)
  # with a1 = a there is no second part, and z1 is the time to death
  whole <- deaths_settings(30, list(
    onset_to_death=list(mean=9.3, sd=5.5),
    dated_by_death=list(fraction=0.3, position=list(share_shape=position_shape))))$dated_by_death
  expect_identical(whole$lag_probabilities, p_die)
  expect_identical(whole$position$after_probabilities, 1)
  expect_identical(whole$max_lag, length(p_die))
})

test_that("the survival of a misdated case is the mixture over its day of death, with the effects, weekdays and backlog", {
  report_days <- c(4, 5, 7, 8, 11, 12, 15, 19, 22, 26, 33)
  effects <- position_effects(report_days)
  multiplier <- position_multiplier(effects)
  # the multipliers are those of the effects as specified
  expect_equal(vapply(1:40, multiplier, numeric(1)),
               c(rep(1, 3), effects$report_multipliers$values, rep(1, 40 - 33)))
  delay <- list(mean=3, sd=2.5)
  p2 <- deaths_settings(30, list(onset_to_death=list(mean=9.3, sd=5.5),
                                 dated_by_death=list(fraction=0.3, position=TRUE))
                        )$dated_by_death$position$after_probabilities
  for(day in c(4, 6, 9)) {
    days <- report_days[report_days >= day]
    share <- effects$backlog$share[day]
    log_mixed <- log_dated_recorded(days, day, delay, effects, p2)
    expected <- vapply(days, function(t) reference_survival_mixed(
      t, day, p2, mean=3, sd=2.5, multiplier=multiplier, share=share, rate=0.15,
      hold=2), numeric(1))
    expect_equal(exp(log_mixed), expected, tolerance=1e-10)
    # it is what the forward filter is given for the misdated cases
    survivals <- joint_survivals(days, day, list(mean=6, sd=4), delay, NULL,
                                 "recording", effects, after_probabilities=p2)
    expect_equal(survivals$log_dated, log_mixed)
    expect_equal(joint_path_exact_forward(
      seq_along(days), seq_along(days), 0, length(days), 0, survivals)$arguments$p_dated,
      survival_steps(log_mixed))
    # it lies in [0, 1] and does not increase
    expect_true(all(log_mixed <= 0 & is.finite(log_mixed)))
    expect_true(all(diff(log_mixed) <= 0))
    # with the death on the day recorded for, the survival of the delay from
    # that day
    expect_equal(log_dated_recorded(days, day, delay, effects, 1),
                 log_undetected_days(days, day, delay, effects_from_death(effects)))
    expect_equal(log_dated_recorded(days, day, delay, effects, c(1, 0, 0)),
                 log_undetected_days(days, day, delay, effects_from_death(effects)),
                 tolerance=1e-10)
    expect_equal(exp(log_dated_recorded(days, day, delay, effects, 1)),
                 vapply(days, function(t) reference_survival_recorded(
                   t, day, 3, 2.5, multiplier, share, 0.15, 2), numeric(1)),
                 tolerance=1e-10)
    expect_equal(joint_survivals(days, day, list(mean=6, sd=4), delay, NULL,
                                 "recording", effects)$log_dated,
                 log_dated_recorded(days, day, delay, effects, 1))
  }
  # later deaths are reported later
  days <- report_days[report_days >= 4]
  expect_true(all(log_dated_recorded(days, 4, delay, effects, p2) >=
                    log_dated_recorded(days, 4, delay, effects, 1)))
  # a cache of the delay given gives the same
  cached <- effects
  cached$lag_cache <- new_lag_cache()
  expect_equal(log_dated_recorded(days, 4, delay, cached, p2),
               log_dated_recorded(days, 4, delay, effects, p2))
})

test_that("the share of the backlog of a misdated case is that of the day it is recorded for", {
  report_days <- c(6, 7, 9, 12, 14, 20, 27)
  effects <- position_effects(report_days)
  multiplier <- position_multiplier(effects)
  # recorded for day 6 and dying two days later, on day 8, whose share of
  # the backlog differs
  share <- effects$backlog$share
  expect_gt(abs(share[6] - share[8]), 0.5)
  p2 <- c(0, 0, 1)
  days <- report_days
  f_reference <- function(share)
    vapply(days, function(t) reference_survival_mixed(
      t, 6, p2, mean=3, sd=2.5, multiplier=multiplier, share=share, rate=0.15,
      hold=2), numeric(1))
  log_mixed <- log_dated_recorded(days, 6, list(mean=3, sd=2.5), effects, p2)
  expect_equal(exp(log_mixed), f_reference(share[6]), tolerance=1e-10)
  expect_gt(max(abs(exp(log_mixed) - f_reference(share[8]))), 0.05)
})

# the probability of every sequence of reports of the cases recorded for one
# day, by enumerating what happens to each case: the outcome of a case is the
# report at which it is first reported and that at which it is first shown
# as dead (n_reports + 1 if neither happens by the last report), and its
# probabilities are given as a matrix with a row for each of the first and a
# column for each of the second. The cases are independent, so the
# probability of each sequence is the sum over the outcomes of every case
# which give it of the product of their probabilities.
enumerate_reports <- function(outcomes, n_reports) {
  sequences <- matrix(0, 1, 2 * n_reports)
  p <- 1
  for(outcome in outcomes) {
    cells <- which(outcome > 0, arr.ind=TRUE)
    rows <- rep(seq_len(nrow(sequences)), each=nrow(cells))
    added <- t(apply(cells, 1, function(cell)
      c(seq_len(n_reports) >= cell[1], seq_len(n_reports) >= cell[2])))
    sequences <- sequences[rows, , drop=FALSE] + added[rep(seq_len(nrow(cells)),
                                                          length(p)), , drop=FALSE]
    p <- p[rows] * outcome[cells][rep(seq_len(nrow(cells)), length(p))]
    keys <- apply(sequences, 1, paste, collapse=",")
    p <- as.vector(rowsum(p, keys, reorder=FALSE))
    sequences <- sequences[!duplicated(keys), , drop=FALSE]
  }
  list(sequences=sequences, p=p)
}

# the probability that an event which happens after a delay has happened
# first by each report (the last element: not by the last report), from the
# probability that it has not happened by each report
first_by_report <- function(survival) -diff(c(1, survival, 0))

# the outcomes of the cases recorded for day d with reports on days, written
# from the model: the reporting delay of every case and the delay from a death
# to its being recorded are gamma, as is the time from onset to death of a
# correctly dated fatal case (h, from one day); a misdated case dies z2 days
# after day d (p2, from zero) and is reported, shown as dead, when its death
# is recorded
case_outcomes <- function(d, days, reporting, recording, h, p2_cases, n_fatal,
                          n_survivors) {
  n <- length(days)
  f_survival <- function(x, delay)
    stats::pgamma(x, shape=delay$mean^2 / delay$sd^2, rate=delay$mean / delay$sd^2,
                  lower.tail=FALSE)
  reported <- first_by_report(f_survival(days - d + 0.5, reporting))
  # the death of a case dying on day d + z recorded first by each report
  f_recorded <- function(z)
    first_by_report(ifelse(days - d >= z, f_survival(days - d - z + 0.5, recording), 1))
  survivor <- matrix(0, n + 1, n + 1)
  survivor[, n + 1] <- reported
  fatal <- matrix(0, n + 1, n + 1)
  for(z in seq_along(h)) {
    recorded <- f_recorded(z)
    for(a in seq_len(n + 1))
      for(b in seq_len(n + 1))
        fatal[a, max(a, b)] <- fatal[a, max(a, b)] + h[z] * reported[a] * recorded[b]
  }
  misdated <- lapply(p2_cases, function(p2) {
    outcome <- matrix(0, n + 1, n + 1)
    for(z in seq_along(p2) - 1)
      diag(outcome) <- diag(outcome) + p2[z + 1] * f_recorded(z)
    outcome
  })
  c(misdated, rep(list(fatal), n_fatal), rep(list(survivor), n_survivors))
}

test_that("the forward filter gives the probability of every sequence of reports which enumerating each case's day of death and delays gives", {
  d <- 3
  days <- c(3, 4, 6, 9)
  reporting <- list(mean=2, sd=1.5)
  recording <- list(mean=1.5, sd=1)
  h <- c(0.5, 0.3, 0.2)
  p2 <- c(0.3, 0.4, 0.3)
  survivals <- joint_survivals(days, d, reporting, recording, NULL, "recording",
                               list(before=list(lag_probabilities=h)),
                               after_probabilities=p2)
  f_filter <- function(sequences, n_misdated, n_fatal, n_survivors)
    apply(sequences, 1, function(x) exp(joint_path_exact_forward(
      x[1:4], x[5:8], n_fatal, n_misdated, n_survivors, survivals)$logp))
  for(counts in list(c(3, 1, 2), c(2, 2, 1))) {
    outcomes <- case_outcomes(d, days, reporting, recording, h,
                              rep(list(p2), counts[1]), counts[2], counts[3])
    enumerated <- enumerate_reports(outcomes, 4)
    expect_equal(sum(enumerated$p), 1)
    expect_gt(nrow(enumerated$sequences), 200)
    filtered <- f_filter(enumerated$sequences, counts[1], counts[2], counts[3])
    expect_equal(filtered, enumerated$p, tolerance=1e-10)
    # so every other sequence is impossible under the filter too
    expect_equal(sum(filtered), 1)
  }
  # where z2 depends on z1, the misdated cases of a day do not share a
  # delay, and the filter with their mixed delay is not exact
  p2_by_case <- list(c(0.8, 0.2, 0), c(0.2, 0.6, 0.2), c(0, 0.2, 0.8))
  p2_mixed <- Reduce(`+`, p2_by_case) / 3
  survivals <- joint_survivals(days, d, reporting, recording, NULL, "recording",
                               list(before=list(lag_probabilities=h)),
                               after_probabilities=p2_mixed)
  enumerated <- enumerate_reports(case_outcomes(d, days, reporting, recording, h,
                                                p2_by_case, 1, 2), 4)
  filtered <- f_filter(enumerated$sequences, 3, 1, 2)
  expect_gt(max(abs(filtered - enumerated$p)), 1e-3)
  # and with the same delay for each, they agree again
  enumerated <- enumerate_reports(case_outcomes(d, days, reporting, recording, h,
                                                rep(list(p2_mixed), 3), 1, 2), 4)
  expect_equal(f_filter(enumerated$sequences, 3, 1, 2), enumerated$p, tolerance=1e-10)
})

# reports simulated with misdated cases recorded part of the way to death,
# with the backlog and the time from onset to death of the fits
position_example <- function(n_days=30, seed=41, position=TRUE) {
  set.seed(seed)
  serial_parameters <- list(mean=5, sd=3)
  cases <- true_cases(n_days + 40, function(t) 1, kappa=1000, serial_parameters,
                      initial_parameters=list(mean=30, length=10))
  sim <- simulate_cases_deaths(cases[41:(n_days + 40)], seq_len(n_days + 10),
                               list(mean=6, sd=4), 0.5, list(mean=3, sd=2.5),
                               triggered=0.5,
                               onset_to_death=list(mean=9.3, sd=5.5),
                               dated_by_death=list(fraction=0.35, position=position),
                               recording_backlog=list(share=0.2, rate=0.05,
                                                      hold_days=3))
  list(snapshots=sim$snapshots, cases_before=cases[1:40],
       serial_parameters=serial_parameters,
       initial_cases=sim$snapshots %>%
         dplyr::group_by(time_onset) %>%
         dplyr::summarise(cases_true=max(cases_reported) + 2))
}

position_fit <- function(example, niterations, position=TRUE, exact=TRUE,
                         ridge=NULL, fraction=list(), ...) {
  snapshots <- example$snapshots
  mcmc_single(
    niterations=niterations,
    snapshot_with_Rt_index_df=snapshots %>%
      dplyr::select(time_onset, time_reported, cases_reported) %>%
      dplyr::mutate(Rt_index=1),
    priors=list(Rt=list(shape=1, rate=1),
                reporting=list(mean_mu=6, mean_sigma=5, sd_mu=4, sd_sigma=4),
                max_cases=300,
                overdispersion=list(mean=20, sd=20),
                Rt_changepoints=list(),
                reporting_changepoints=list(min_length=14)),
    serial_parameters=example$serial_parameters,
    initial_cases_true=example$initial_cases,
    initial_reporting_parameters=list(mean=6, sd=4),
    initial_Rt=dplyr::tibble(Rt_index=1, Rt=1),
    print_to_screen=FALSE,
    is_negative_binomial=TRUE,
    initial_overdispersion=10,
    cases_before=example$cases_before,
    deaths=list(
      snapshot=snapshots %>%
        dplyr::select(time_onset, time_reported, deaths_reported),
      settings=list(onset_to_death=list(mean=9.3, sd=5.5), fatal_delay=FALSE,
                    triggered=list(initial=0.5, n_steps=0),
                    recording_backlog=list(hold_days=3),
                    dated_by_death=list(fraction=fraction, exact=exact,
                                        ridge=ridge, position=position))),
    ...)
}

# the state of the model of deaths after a few iterations of a fit with
# misdated cases recorded part of the way to death
position_test_state <- function(position=TRUE, exact=TRUE, ridge=NULL, seed=3) {
  set.seed(seed)
  fit <- position_fit(position_example(), 3, position=position, exact=exact,
                      ridge=ridge, return_state=TRUE)
  v <- fit$state$variables
  list(state=v$deaths_state, cases_true=v$deaths_cases_true,
       reporting_by_onset=v$deaths_running_df, revision=v$revision_current,
       settings=v$deaths_options,
       recorded=v$deaths_cases_true - v$deaths_state$dated$out +
         v$deaths_state$dated$into,
       renewal=dated_renewal(rep(1, length(v$deaths_cases_true)), 10, TRUE,
                             list(mean=5, sd=3), 40, v$cases_before))
}

test_that("in a state, a misdated case has the delay and share of the backlog of the day it is recorded for", {
  # recorded for day 6 and dying on that day or later, with the delay to a
  # death being recorded and the share of the backlog changing on day 8
  test <- position_test_state()
  state <- test$state
  p2 <- state$dated$position$after_probabilities
  state$recording$starts <- c(1, 8)
  state$recording$parameters <- dplyr::tibble(mean=c(2, 6), sd=c(1.5, 3))
  state$recording_effects$backlog$share[] <- 0.1
  state$recording_effects$backlog$share[8] <- 0.8
  state <- deaths_refresh_survivals(state, test$reporting_by_onset, test$revision)
  days <- state$reports[[6]]$time_reported
  f_mixed <- function(delay, share) {
    effects <- state$recording_effects
    effects$backlog$share[] <- share
    log_dated_recorded(days, 6, delay, effects, p2)
  }
  expect_equal(state$survivals[[6]]$log_dated, f_mixed(list(mean=2, sd=1.5), 0.1))
  expect_gt(max(abs(state$survivals[[6]]$log_dated - f_mixed(list(mean=6, sd=3), 0.1))),
            0.05)
  expect_gt(max(abs(state$survivals[[6]]$log_dated - f_mixed(list(mean=2, sd=1.5), 0.8))),
            0.05)
})

# the log-probability of the fatal cases dated by death, written from the
# model: for each onset day the multinomial over dated by onset and dated by
# death recorded after each number of days (or after the last day), with
# the probabilities of the time to the day recorded for given, and the
# binomial of its fatal cases; for each recorded day the Poisson probability
# of those with onset before the first day, and its reports
position_dated_logp <- function(state, cases_true, p1) {
  dated <- state$dated
  n_days <- length(state$onset_times)
  f <- dated$fraction
  logp <- 0
  for(k in seq_len(n_days)) {
    v <- dated$lags[k, ]
    u <- dated$unseen[k]
    fatal <- state$deaths_true[k]
    correct <- fatal - sum(v) - u
    if(correct < 0 | fatal > cases_true[k])
      return(-Inf)
    tail <- sum(p1[seq_along(p1) > n_days - k])
    logp <- logp + lgamma(fatal + 1) - lgamma(correct + 1) - sum(lgamma(v + 1)) -
      lgamma(u + 1) + correct * log(1 - f[k]) + sum(v[v > 0] * log(f[k] * p1[v > 0])) +
      (if(u > 0) u * log(f[k] * tail) else 0) +
      stats::dbinom(fatal, cases_true[k], state$cfr[k], log=TRUE)
  }
  for(t in seq_len(n_days)) {
    logp <- logp + stats::dpois(dated$inflow[t],
                                dated$inflow_weight[t] * state$cfr[1] * f[1], log=TRUE)
    logp <- logp + if(dated$exact) dated_exact_forward(state, cases_true, t)$logp
    else fatal_reports_logp(state$paths[[t]], state$survivals[[t]],
                            state$deaths_true[t] - dated$out[t] + dated$into[t],
                            dated$into[t])
  }
  logp
}

test_that("the labels hold the time to the day recorded for, and the cases unseen and from before the first day use its distribution", {
  for(exact in c(FALSE, TRUE)) {
    test <- position_test_state(exact=exact, seed=4 + exact)
    state <- test$state
    cases_true <- test$cases_true
    dated <- state$dated
    p1 <- position_delays_reference(position_shape / 2)$p1
    expect_equal(dated$lag_probabilities, p1, tolerance=1e-12)
    expect_equal(dated$max_lag, length(p1))
    expect_equal(ncol(dated$lags), length(p1))
    n_days <- length(state$onset_times)
    # a case is unseen if it is recorded for a day after the last
    expect_equal(dated$tail, vapply(seq_len(n_days), function(t)
      sum(p1[seq_along(p1) > n_days - t]), numeric(1)), tolerance=1e-12)
    # those from before the first day recorded for day t: the earlier cases
    # times the probability of being recorded t - s days after onset s
    before <- position_example()$cases_before
    expect_equal(dated$inflow_weight, vapply(seq_len(n_days), function(t) {
      s <- seq_along(before) - length(before)
      lag <- t - s
      sum(before[lag <= length(p1)] * p1[lag[lag <= length(p1)]])
    }, numeric(1)), tolerance=1e-12)
    expect_gt(sum(dated$inflow_weight), 0)
    expect_gt(sum(dated$lags), 0)
    # and the multinomial of the labels uses it
    expect_equal(dated_logp(state, cases_true), position_dated_logp(state, cases_true, p1),
                 tolerance=1e-10)
    # the survivals of each recorded day are those of the misdated cases
    k <- which.max(dated$into)
    days <- state$reports[[k]]$time_reported
    delays <- deaths_delay_by_onset(state)
    expect_equal(state$survivals[[k]]$log_dated,
                 log_dated_recorded(days, k, delay_of_onset(delays, k),
                                    state$recording_effects,
                                    dated$position$after_probabilities))
  }
})

test_that("changes to the fatal cases dated by death give the right change in log-probability with the day recorded for part of the way to death", {
  for(exact in c(FALSE, TRUE)) {
    test <- position_test_state(exact=exact, seed=6 + exact)
    state <- test$state
    cases_true <- test$cases_true
    p1 <- state$dated$lag_probabilities
    n_days <- length(cases_true)
    f_full <- function(state, cases_true)
      position_dated_logp(state, cases_true, p1) +
      dated_renewal_logp(cases_true, test$renewal)
    set.seed(34)
    n_applied <- c(relabel=0, redate=0, unseen=0, inflow=0)
    for(i in 1:(if(exact) 200 else 400)) {
      type <- sample(names(n_applied), 1)
      # the times to the day recorded for are drawn as the sampler draws them
      move <- list(type=type, day=sample.int(n_days, 1),
                   lag=sample.int(length(p1), 1, prob=p1),
                   lag_to=sample.int(length(p1), 1, prob=p1), sign=sample(c(-1, 1), 1))
      proposed <- dated_apply_move(state, cases_true, move)
      if(is.null(proposed))
        next
      logp_proposed <- f_full(proposed$state, proposed$cases_true)
      if(!is.finite(logp_proposed))
        next
      # the log of the acceptance ratio, as deaths_dated_step computes it
      log_ratio <- dated_logp(proposed$state, proposed$cases_true, proposed$cohorts,
                              proposed$days) +
        dated_renewal_logp(proposed$cases_true, test$renewal) -
        dated_logp(state, cases_true, proposed$cohorts, proposed$days) -
        dated_renewal_logp(cases_true, test$renewal)
      expect_equal(log_ratio, logp_proposed - f_full(state, cases_true), tolerance=1e-10)
      reverse <- move
      reverse$sign <- -move$sign
      if(type == "redate") {
        reverse$lag <- move$lag_to
        reverse$lag_to <- move$lag
      }
      back <- dated_apply_move(proposed$state, proposed$cases_true, reverse)
      expect_equal(back$state$dated[c("lags", "unseen", "inflow", "out", "into")],
                   state$dated[c("lags", "unseen", "inflow", "out", "into")])
      expect_equal(back$cases_true, cases_true)
      state <- proposed$state
      cases_true <- proposed$cases_true
      n_applied[type] <- n_applied[type] + 1
    }
    expect_true(all(n_applied > 3))
    expect_equal(cases_true - state$dated$out + state$dated$into, test$recorded)
  }
  # the step itself runs and keeps the cases recorded for each day
  test <- position_test_state()
  set.seed(35)
  stepped <- deaths_dated_step(test$state, test$cases_true, test$renewal)
  expect_false(identical(stepped$cases_true, test$cases_true))
  expect_equal(stepped$cases_true - stepped$state$dated$out + stepped$state$dated$into,
               test$recorded)
  expect_true(max(which(colSums(stepped$state$dated$lags) > 0)) <=
                length(test$state$dated$lag_probabilities))
})

test_that("the update of the delay to a death being recorded uses the survival of the misdated cases", {
  test <- position_test_state()
  # the survivals held, at the current delays
  state <- deaths_refresh_survivals(test$state, test$reporting_by_onset, test$revision)
  delays <- deaths_delay_by_onset(state)
  ks <- which(state$dated$into > 0)
  expect_gt(length(ks), 3)
  for(k in ks[1:3]) {
    # at the current delay, the log-probability from the survivals held
    logp <- deaths_recording_day_logp(state, k, delay_of_onset(delays, k))
    dated <- purrr::map(dated_counts(state), k)
    expect_equal(logp, fatal_reports_logp(state$paths[[k]], state$survivals[[k]],
                                          dated$fatal_recorded, dated$into))
    # which are worked out again, with the misdated cases' own survival
    other <- state
    other$survivals[[k]]$log_dated <- other$survivals[[k]]$log_dated / 2
    expect_equal(deaths_recording_day_logp(other, k, delay_of_onset(delays, k)), logp)
    main <- state
    main$dated$position <- NULL
    expect_false(isTRUE(all.equal(
      deaths_recording_day_logp(main, k, delay_of_onset(delays, k)), logp)))
  }
})

test_that("the joint move along the ridge uses the time to the day recorded for", {
  x <- position_test_state(ridge=TRUE, seed=8)
  state <- x$state
  dated <- state$dated
  p1 <- dated$lag_probabilities
  expect_false(isTRUE(all.equal(p1, lag_probabilities_gamma(9.3, 5.5))))
  # the probability of choosing a change to dated by death
  by_onset <- state$deaths_true - dated$out
  by_onset[1] <- 0
  day <- which.max(by_onset)
  lags <- seq_len(min(dated$max_lag, day - 1))
  for(lag in lags[1:3])
    expect_equal(exp(ridge_relabel_logq(state, list(day=day, lag=lag, sign=1))),
                 by_onset[day] / sum(by_onset) * p1[lag] / sum(p1[lags]))
  total_up <- 0
  for(day in which(by_onset > 0))
    for(lag in seq_len(min(dated$max_lag, day - 1)))
      total_up <- total_up + exp(ridge_relabel_logq(state, list(day=day, lag=lag, sign=1)))
  expect_equal(total_up, 1)
  # the reverse of a joint move returns to the state, with the reverse ratio
  set.seed(11)
  for(sign in c(1, -1)) {
    forward <- ridge_propose(state, x$cases_true, x$reporting_by_onset, x$revision,
                             sign, 5, 0.04)
    expect_false(is.null(forward))
    expect_equal(forward$cases_true - forward$state$dated$out + forward$state$dated$into,
                 x$recorded)
    undo <- purrr::map(rev(forward$moves), ~list(type="relabel", day=.$day, lag=.$lag,
                                                lag_to=NA, sign=-sign))
    back <- ridge_propose(forward$state, forward$cases_true, x$reporting_by_onset,
                          x$revision, -sign, 5, -0.04, undo)
    for(name in c("lags", "out", "into", "unseen", "inflow", "fraction"))
      expect_equal(back$state$dated[[name]], state$dated[[name]], info=name)
    expect_equal(back$cases_true, x$cases_true)
    expect_equal(back$state$survivals, state$survivals)
    expect_equal(back$log_q, -forward$log_q)
    logp <- ridge_logp(state, x$cases_true, x$renewal, x$settings)$logp
    logp_forward <- ridge_logp(forward$state, forward$cases_true, x$renewal,
                               x$settings)$logp
    expect_true(is.finite(logp) & is.finite(logp_forward))
    expect_equal((logp_forward - logp + forward$log_jacobian + forward$log_q) +
                   (logp - logp_forward + back$log_jacobian + back$log_q), 0)
    # the survivals of the proposed state are those of the misdated cases
    k <- which.max(forward$state$dated$into)
    expect_equal(forward$state$survivals[[k]]$log_dated, log_dated_recorded(
      forward$state$reports[[k]]$time_reported, k,
      delay_of_onset(deaths_delay_by_onset(forward$state), k),
      forward$state$recording_effects, dated$position$after_probabilities))
  }
  # and its log-probability is that of the model with the split summed out
  terms <- ridge_logp(state, x$cases_true, x$renewal, x$settings)
  recording <- state$recording
  prior <- reporting_centre_logp_prior(recording$centre, recording$spread,
                                       x$settings$delay_prior,
                                       x$settings$delay_changepoints$spread_sd) +
    sum(reporting_segments_logp_prior(recording$parameters, recording$centre,
                                      recording$spread)) +
    stats::dbeta(dated$fraction[1], dated$fraction_prior$alpha,
                 dated$fraction_prior$beta, log=TRUE)
  expect_equal(terms$logp, prior + position_dated_logp(state, x$cases_true, p1) +
                 dated_renewal_logp(x$cases_true, x$renewal), tolerance=1e-10)
})

test_that("the simulator records a misdated case for a day part of the way to death and reports it when its death is recorded", {
  # every case of onset day 1 is fatal and misdated, recorded two days later
  # and dying two days after that
  set.seed(12)
  sim <- simulate_cases_deaths(
    c(300, rep(0, 9)), 1:20, list(mean=1, sd=0.5), c(1, rep(0, 9)),
    list(mean=0.5, sd=0.2), onset_to_death=list(mean=9.3, sd=5.5),
    dated_by_death=list(fraction=1, position=list(lag_probabilities=c(0, 1),
                                                  after_probabilities=c(0, 0, 1))))
  expect_equal(sim$truth$n_dated_by_death, c(0, 0, 300, rep(0, 7)))
  shown <- dplyr::filter(sim$snapshots, time_onset == 3)
  # none is reported before its death, on day 5, and then each shown as dead
  expect_equal(shown$cases_reported[shown$time_reported < 5], c(0, 0))
  expect_gt(shown$cases_reported[shown$time_reported == 5], 100)
  expect_equal(shown$deaths_reported, shown$cases_reported)
  # without the setting the simulator draws what it drew before
  f_sim <- function(...) {
    set.seed(13)
    simulate_cases_deaths(rep(50, 10), 1:15, list(mean=4, sd=3), 0.5,
                          list(mean=3, sd=2), onset_to_death=list(mean=9.3, sd=5.5),
                          recording_backlog=list(share=0.3, rate=0.1),
                          dated_by_death=c(list(fraction=0.4), list(...)))
  }
  expect_identical(f_sim(position=NULL), f_sim())
  expect_identical(f_sim(position=FALSE), f_sim())
})

test_that("reports simulated with the day recorded for part of the way to death have the frequencies the forward filter gives", {
  # onset day 1: two cases, both fatal and misdated, recorded for day 2 and
  # dying 0 to 2 days later; onset day 2: two cases, each fatal with
  # probability one half, none misdated. The share of the backlog of day 2
  # differs from those of the days of death after it.
  report_days <- c(2, 3, 5, 8)
  share <- c(0.1, 0.6, 0.05, 0.05, 0.05, 0.05)
  p2 <- c(0.3, 0.4, 0.3)
  reporting <- list(mean=2, sd=1.5)
  recording <- list(mean=1.5, sd=1)
  f_simulate <- function()
    simulate_cases_deaths(
      c(2, 2, 0, 0, 0, 0), report_days, reporting, c(1, 0.5, 0, 0, 0, 0), recording,
      onset_to_death=list(mean=3, sd=1.5),
      dated_by_death=list(fraction=c(1, 0, 0, 0, 0, 0),
                          position=list(lag_probabilities=1, after_probabilities=p2)),
      recording_backlog=list(share=share, rate=0.3, hold_days=1))
  set.seed(14)
  n <- 2000
  replicates <- lapply(seq_len(n), function(i) {
    sim <- f_simulate()
    shown <- dplyr::filter(sim$snapshots, time_onset == 2)
    list(fatal=sim$truth$deaths_onset[2], misdated=sim$truth$n_dated_by_death[2],
         cases=shown$cases_reported, deaths=shown$deaths_reported)
  })
  expect_true(all(purrr::map_dbl(replicates, "misdated") == 2))
  effects <- list(
    before=list(lag_probabilities=lag_probabilities_gamma(3, 1.5)),
    backlog=initial_backlog(6, report_days, backlog_settings(list(
      share_week=1, initial_share=share, initial_rate=0.3, hold_days=1))))
  survivals <- joint_survivals(report_days, 2, reporting, recording, NULL, "recording",
                               effects, after_probabilities=p2)
  fatal <- purrr::map_dbl(replicates, "fatal")
  keys <- purrr::map_chr(replicates, ~paste(c(.$cases, .$deaths), collapse=","))
  for(n_fatal in 0:2) {
    # every sequence possible, as enumerated above without the backlog, and
    # its probability from the forward filter
    possible <- enumerate_reports(case_outcomes(
      2, report_days, reporting, recording, c(0.5, 0.5), list(p2, p2), n_fatal,
      2 - n_fatal), 4)$sequences
    probabilities <- apply(possible, 1, function(x) exp(joint_path_exact_forward(
      x[1:4], x[5:8], n_fatal, 2, 2 - n_fatal, survivals)$logp))
    names(probabilities) <- apply(possible, 1, paste, collapse=",")
    expect_equal(sum(probabilities), 1)
    expect_gof(table(keys[fatal == n_fatal]), probabilities,
               paste("sequences of reports with", n_fatal, "fatal cases of onset day 2"))
  }
  # and the deaths shown in each report have their expected number
  dead_mixed <- 1 - exp(survivals$log_dated)
  dead_fatal <- (1 - exp(survivals$log_reported_fatal)) * (1 - exp(survivals$log_recorded))
  for(j in seq_along(report_days)) {
    deaths <- purrr::map_dbl(replicates, ~.$deaths[j])
    expect_mean_z(deaths - 2 * dead_mixed[j] - fatal * dead_fatal[j], 0,
                  paste("deaths shown in report", j))
  }
})

test_that("with the day of death as the day recorded for the model is the main model", {
  test <- position_test_state(position=NULL, seed=9)
  main <- deaths_refresh_survivals(test$state, test$reporting_by_onset, test$revision)
  cases_true <- test$cases_true
  for(after in list(1, c(1, 0, 0))) {
    nested <- main
    nested$dated$position <- list(share_shape=position_shape, after_probabilities=after)
    nested <- deaths_refresh_survivals(nested, test$reporting_by_onset, test$revision)
    expect_equal(nested$survivals, main$survivals, tolerance=1e-10)
    expect_equal(dated_logp(nested, cases_true), dated_logp(main, cases_true),
                 tolerance=1e-10)
    expect_equal(deaths_log_likelihood(nested, cases_true),
                 deaths_log_likelihood(main, cases_true), tolerance=1e-10)
    for(k in c(3, 12, 24)) {
      expect_equal(dated_exact_forward(nested, cases_true, k)$logp,
                   dated_exact_forward(main, cases_true, k)$logp, tolerance=1e-10)
      expect_equal(deaths_recording_day_logp(nested, k, list(mean=3, sd=2)),
                   deaths_recording_day_logp(main, k, list(mean=3, sd=2)),
                   tolerance=1e-10)
    }
  }
  # and on the data sets of the reference fits (helper-reference.R), with
  # the settings of the fits
  for(seed in c(101, 202)) {
    fit <- reference_fit(reference_data(seed), FALSE, 1, seed=seed + 1)
    v <- fit$state$variables
    f_refresh <- function(state)
      deaths_refresh_survivals(state, v$deaths_running_df, v$revision_current)
    main <- f_refresh(v$deaths_state)
    nested <- main
    nested$dated$position <- list(share_shape=position_shape, after_probabilities=c(1, 0))
    nested <- f_refresh(nested)
    cases_true <- v$deaths_cases_true
    expect_equal(nested$survivals, main$survivals, tolerance=1e-10)
    expect_equal(dated_logp(nested, cases_true), dated_logp(main, cases_true),
                 tolerance=1e-10)
    expect_equal(deaths_log_likelihood(nested, cases_true),
                 deaths_log_likelihood(main, cases_true), tolerance=1e-10)
    expect_true(all(is.finite(deaths_log_likelihood(main, cases_true))))
  }
  # a fit with a1 = a draws what the main model draws
  example <- position_example(seed=42, position=NULL)
  set.seed(10)
  fit_main <- position_fit(example, 3, position=NULL)
  set.seed(10)
  fit_nested <- position_fit(example, 3, position=list(share_shape=position_shape))
  for(name in setdiff(names(fit_main), c("state", "iterations_done")))
    expect_equal(fit_nested[[name]], fit_main[[name]], info=name)
})

test_that("a fit with the day recorded for part of the way to death runs and can be taken up again", {
  example <- position_example()
  set.seed(15)
  whole <- position_fit(example, 6)
  set.seed(15)
  first <- position_fit(example, 3, return_state=TRUE)
  set.seed(99)
  stats::runif(5)
  second <- position_fit(example, 3, resume_from=first$state)
  for(name in setdiff(names(whole), c("state", "iterations_done")))
    expect_equal(dplyr::bind_rows(first[[name]], second[[name]]), whole[[name]],
                 info=name)
  expect_equal(first$state$variables$deaths_state$dated$position,
               deaths_settings(30, list(onset_to_death=list(mean=9.3, sd=5.5),
                                        dated_by_death=list(fraction=list(), position=TRUE))
                               )$dated_by_death$position)
  fraction <- whole$deaths_parameters$dated_fraction
  expect_gt(length(unique(fraction)), 1)
  both <- dplyr::left_join(whole$deaths, whole$cases, by=c("iteration", "time_onset"))
  expect_true(all(both$deaths_true <= both$cases_true))
})

# a short fit for predictions, with wrong records (as in test-predict.R)
position_predict_fit <- function(example, position) {
  snapshots <- example$snapshots
  mcmc_single(
    niterations=3,
    snapshot_with_Rt_index_df=snapshots %>%
      dplyr::select(time_onset, time_reported, cases_reported) %>%
      dplyr::mutate(Rt_index=1),
    priors=list(Rt=list(shape=1, rate=1),
                reporting=list(mean_mu=6, mean_sigma=5, sd_mu=4, sd_sigma=4),
                max_cases=300, overdispersion=list(mean=20, sd=20),
                revisions=list(eta_mu=0.1, eta_sigma=0.1, rho_alpha=1, rho_beta=1),
                max_erroneous=3, Rt_changepoints=list(),
                reporting_changepoints=list(min_length=14)),
    serial_parameters=example$serial_parameters,
    initial_cases_true=example$initial_cases,
    initial_reporting_parameters=list(mean=6, sd=4),
    initial_Rt=dplyr::tibble(Rt_index=1, Rt=1),
    initial_revision_parameters=list(eta=0.02, rho=0.1),
    print_to_screen=FALSE, is_negative_binomial=TRUE, initial_overdispersion=10,
    cases_before=example$cases_before, return_state=TRUE,
    deaths=list(
      snapshot=snapshots %>%
        dplyr::select(time_onset, time_reported, deaths_reported),
      settings=list(onset_to_death=list(mean=9.3, sd=5.5), fatal_delay=FALSE,
                    triggered=list(initial=0.5, n_steps=0L),
                    recording_backlog=list(hold_days=3),
                    dated_by_death=list(fraction=list(), exact=TRUE,
                                        position=position))))
}

predict_draw <- function(fit, iteration) {
  lapply(fit[c("cases", "reporting_days", "revisions", "reporting_effects",
               "deaths", "deaths_parameters")],
         function(df) df[df$iteration == iteration, ])
}

test_that("predictions use the survival of the misdated cases, and are those of the main model at the nested case", {
  example <- position_example(position=NULL)
  set.seed(16)
  fit <- position_predict_fit(example, NULL)
  v <- fit$state$variables
  draw <- predict_draw(fit, 3)
  nested <- v
  nested$deaths_state$dated$position <- list(share_shape=position_shape,
                                             after_probabilities=c(1, 0))
  set.seed(17)
  main <- predict_reports(v, draw, 41:45)
  set.seed(17)
  expect_equal(predict_reports(nested, draw, 41:45), main)

  # a fit with the day recorded for part of the way to death
  example <- position_example()
  set.seed(18)
  fit <- position_predict_fit(example, TRUE)
  v <- fit$state$variables
  draw <- predict_draw(fit, 3)
  prediction <- predict_state_from_draw(v, draw)
  k <- which.max(prediction$deaths_state$dated$into)
  days <- c(v$deaths_state$reports[[k]]$time_reported, 41:45)
  terms <- predict_onset_terms(prediction, k, days)
  expect_equal(terms$survivals$log_dated, log_dated_recorded(
    days, k, delay_of_onset(deaths_delay_by_onset(prediction$deaths_state), k),
    prediction$deaths_state$recording_effects,
    v$deaths_state$dated$position$after_probabilities))
  out <- predict_reports(v, draw, 41:45)
  expect_equal(sort(unique(out$time_reported)), 41:45)
  expect_true(all(out$deaths_reported <= out$cases_reported))
})
