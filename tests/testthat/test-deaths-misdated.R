# A supplied time from onset to death, and an estimated fraction of fatal
# cases which are recorded for their day of death

test_that("the survival of a delay which runs from a death is a sum over days of death", {
  h <- lag_probabilities_gamma(9.3, 5.5)
  expect_equal(sum(h), 1)
  f_log_survival <- function(delay)
    stats::pgamma(delay, shape=1.44, rate=0.48, lower.tail=FALSE, log.p=TRUE)
  lags <- c(0, 1, 4, 9, 20, 60)
  direct <- purrr::map_dbl(lags, function(k) {
    z <- seq_along(h)
    sum(h[z > k]) + sum(h[z <= k] * exp(f_log_survival(k - z[z <= k] + 0.5)))
  })
  expect_equal(exp(log_survival_after_event(lags, h, f_log_survival)), direct)
  # nothing can be recorded on the day of onset, and the order does not matter
  expect_equal(log_survival_after_event(0, h, f_log_survival), 0)
  expect_equal(log_survival_after_event(rev(lags), h, f_log_survival),
               rev(log_survival_after_event(lags, h, f_log_survival)))
  # through log_undetected_days, with and without multipliers of one
  effects <- list(before=list(lag_probabilities=h))
  parameters <- list(mean=3, sd=2.5)
  log_s <- log_undetected_days(5 + lags, 5, parameters, effects)
  expect_equal(log_s, log_survival_after_event(
    lags, h, function(delay) delay_log_survival(delay, 3, 2.5)))
  effects$report_multipliers <- list(day_first=1, values=rep(1, 100))
  expect_equal(log_undetected_days(5 + lags, 5, parameters, effects), log_s)
  expect_true(all(diff(log_s) <= 0))
})

test_that("the settings for the time from onset to death are checked", {
  expect_null(deaths_settings(30)$onset_to_death)
  settings <- deaths_settings(30, list(onset_to_death=list(mean=9.3, sd=5.5)))
  expect_equal(settings$onset_to_death$lag_probabilities,
               lag_probabilities_gamma(9.3, 5.5))
  # the delay estimated is then that from death
  expect_equal(settings$initial_delay, list(mean=3, sd=2.5))
  expect_lt(settings$delay_prior$mean_mu, deaths_settings(30)$delay_prior$mean_mu)
  expect_equal(deaths_settings(30, list(onset_to_death=list(mean=9.3, sd=5.5),
                                        initial_delay=list(mean=5, sd=4)))$initial_delay,
               list(mean=5, sd=4))
  expect_error(deaths_settings(30, list(onset_to_death=list(mean=9.3))),
               "positive mean and sd")
  expect_error(deaths_settings(30, list(onset_to_death=list(mean=9, sd=5, shape=2))),
               "Unknown settings")
})

test_that("the fraction dated by death can be given a prior", {
  onset_to_death <- list(mean=9.3, sd=5.5)
  settings <- deaths_settings(30, list(
    onset_to_death=onset_to_death,
    dated_by_death=list(fraction=list(alpha=2, beta=3, initial=0.3))))
  dated <- settings$dated_by_death
  expect_equal(dated$fraction, rep(0.3, 30))
  expect_equal(dated$fraction_prior, list(alpha=2, beta=3))
  # the time from onset to death and the delay are those of other fatal cases
  expect_equal(dated$lag_probabilities, settings$onset_to_death$lag_probabilities)
  expect_equal(dated$reporting, "recording")
  # a known fraction has no prior, and a fraction of zero means none
  expect_null(deaths_settings(30, list(
    onset_to_death=onset_to_death,
    dated_by_death=list(fraction=0.2)))$dated_by_death$fraction_prior)
  expect_null(deaths_settings(30, list(
    onset_to_death=onset_to_death, dated_by_death=list(fraction=0)))$dated_by_death)
  expect_error(deaths_settings(30, list(
    onset_to_death=onset_to_death,
    dated_by_death=list(fraction=list(alpha=0)))), "must be positive")
  expect_error(deaths_settings(30, list(
    onset_to_death=onset_to_death,
    dated_by_death=list(fraction=list(initial=1)))), "strictly between")
  expect_error(deaths_settings(30, list(
    onset_to_death=onset_to_death,
    dated_by_death=list(fraction=list(gamma=1)))), "Unknown settings")
  expect_error(deaths_settings(30, list(
    dated_by_death=list(fraction=0.2, mean=7, sd=3, reporting="recording"))),
    "onset_to_death")
  # with fatal cases reported when their death is recorded, but not exactly
  expect_type(deaths_settings(30, list(
    onset_to_death=onset_to_death, triggered=list(),
    dated_by_death=list(fraction=list())))$triggered, "list")
  expect_error(deaths_settings(30, list(
    onset_to_death=onset_to_death, triggered=list(),
    dated_by_death=list(fraction=0.2, exact=TRUE))), "exact form")
})

test_that("fatal cases dated by death combine with those reported at death", {
  h <- lag_probabilities_gamma(9.3, 5.5)
  days <- c(1, 2, 3, 5, 8, 12, 17, 25, 40)
  f_survivals <- function(triggered)
    joint_survivals(days, 1, list(mean=7, sd=5), list(mean=3, sd=2.5), NULL,
                    "recording", list(before=list(lag_probabilities=h)), 0, triggered)
  survivals <- f_survivals(0.5)
  expect_true(survivals$is_dated_recording)
  expect_false(survivals$is_dated_same)
  # a case dated by death is reported after the delay from death alone
  expect_equal(survivals$log_dated,
               log_undetected_days(days, 1, list(mean=3, sd=2.5)))
  # with none dated by death, as without them
  none <- mixed_probabilities(survivals, 0)
  with <- mixed_probabilities(survivals, c(0, 0.3))
  expect_equal(with$p_fatal[1, ], as.vector(none$p_fatal))
  expect_equal(with$g[1, ], as.vector(none$g))
  expect_true(all(with$g >= 0 & with$g <= 1))
  expect_true(all(with$p_fatal >= 0 & with$p_fatal <= 1))
  # those dated by death are reported sooner, and shown as dead at once
  expect_gt(with$p_fatal[2, 1], with$p_fatal[1, 1])
  expect_gt(with$g[2, 1], with$g[1, 1])
  expect_equal(as.vector(mixed_probabilities(survivals, 1)$g), rep(1, length(days)))
  # as the proportion reported at death goes to zero, the model without it
  expect_equal(mixed_probabilities(f_survivals(1e-9), 0.3),
               mixed_probabilities(f_survivals(0), 0.3), tolerance=1e-6)
  # the share shown as dead in each report follows from the probabilities
  fraction <- 0.3
  p <- mixed_probabilities(survivals, fraction)
  new <- diff(c(0, 1 - cumprod(1 - as.vector(p$p_fatal))))
  dead <- 0
  alive <- 0
  shown <- numeric(length(days))
  for(j in seq_along(days)) {
    dead <- dead + new[j] * p$g[j] + alive * p$q[j]
    alive <- alive * (1 - p$q[j]) + new[j] * (1 - p$g[j])
    shown[j] <- dead
  }
  own <- 1 - exp(survivals$log_reported_fatal)
  recorded <- 1 - exp(survivals$log_recorded)
  expect_equal(shown, (1 - fraction) * recorded * (0.5 + 0.5 * own) +
                 fraction * (1 - exp(survivals$log_dated)))
})

test_that("cases are simulated with a time from onset to death", {
  set.seed(2)
  cases <- rep(400, 40)
  sim <- simulate_cases_deaths(cases, 1:40, list(mean=6, sd=4), 1,
                               list(mean=3, sd=2.5),
                               onset_to_death=list(mean=9.3, sd=5.5))
  snapshots <- sim$snapshots
  # no death can be shown on the day of onset
  expect_equal(sum(snapshots$deaths_reported[
    snapshots$time_reported == snapshots$time_onset]), 0)
  expect_equal(sim$truth$cases_true, cases)
  set.seed(2)
  sim <- simulate_cases_deaths(cases, 1:40, list(mean=6, sd=4), 1,
                               list(mean=3, sd=2.5),
                               onset_to_death=list(mean=9.3, sd=5.5),
                               dated_by_death=list(fraction=0.4))
  truth <- sim$truth
  expect_gt(sum(truth$n_dated_by_death), 0)
  expect_equal(truth$n_dated_by_death[1], 0)
  # cases dated by death are reported sooner after the day they are counted on
  expect_gt(sum(sim$snapshots$deaths_reported[
    sim$snapshots$time_reported == sim$snapshots$time_onset + 1]), 0)
  expect_error(simulate_cases_deaths(cases, 1:40, list(mean=6, sd=4), 1,
                                     list(mean=3, sd=2.5), death_from="report",
                                     onset_to_death=list(mean=9.3, sd=5.5)),
               "death_from")
  expect_error(simulate_cases_deaths(cases, 1:40, list(mean=6, sd=4), 1,
                                     list(mean=3, sd=2.5),
                                     dated_by_death=list(fraction=0.4)),
               "time from onset to death")
})

misdated_example <- function() {
  set.seed(41)
  n_days <- 30
  serial_parameters <- list(mean=5, sd=3)
  cases <- true_cases(n_days + 40, function(t) 1, kappa=1000, serial_parameters,
                      initial_parameters=list(mean=30, length=10))
  sim <- simulate_cases_deaths(cases[41:(n_days + 40)], seq_len(n_days + 10),
                               list(mean=6, sd=4), 0.5, list(mean=3, sd=2.5),
                               triggered=0.5,
                               onset_to_death=list(mean=9.3, sd=5.5),
                               dated_by_death=list(fraction=0.3))
  list(snapshots=sim$snapshots, cases_before=cases[1:40],
       serial_parameters=serial_parameters,
       initial_cases=sim$snapshots %>%
         dplyr::group_by(time_onset) %>%
         dplyr::summarise(cases_true=max(cases_reported) + 2))
}

misdated_fit <- function(example, niterations, fraction=list(), triggered=list(), ...) {
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
                    triggered=triggered,
                    dated_by_death=list(fraction=fraction))),
    ...)
}

test_that("the fraction dated by death is sampled", {
  example <- misdated_example()
  set.seed(5)
  fit <- misdated_fit(example, 8)
  fraction <- fit$deaths_parameters$dated_fraction
  expect_length(fraction, 8)
  expect_true(all(fraction > 0 & fraction < 1))
  expect_gt(length(unique(fraction)), 1)
  expect_true(all(fit$deaths_parameters$triggered > 0))
  # a known fraction stays where it is
  set.seed(5)
  fit <- misdated_fit(example, 3, fraction=0.25, triggered=FALSE)
  expect_equal(fit$deaths_parameters$dated_fraction, rep(0.25, 3))
})

test_that("a fit with an estimated fraction dated by death can be taken up again", {
  example <- misdated_example()
  set.seed(7)
  whole <- misdated_fit(example, 8)
  set.seed(7)
  first <- misdated_fit(example, 3, return_state=TRUE)
  second <- misdated_fit(example, 5, resume_from=first$state)
  for(name in setdiff(names(whole), c("state", "iterations_done")))
    expect_equal(dplyr::bind_rows(first[[name]], second[[name]]), whole[[name]],
                 info=name)
})

test_that("the delay to a death being recorded can have a heavier tail", {
  expect_false(deaths_settings(30)$recording_tail)
  tail <- deaths_settings(30, list(recording_tail=list(mean_mu=30)))$recording_tail
  expect_equal(tail$mean_mu, 30)
  expect_equal(tail$initial, list(weight=0.1, mean=40))
  expect_error(deaths_settings(30, list(recording_tail=list(weight=0.1))),
               "Unknown settings")
  expect_error(deaths_settings(30, list(recording_tail=list(mean_sigma=0))),
               "must be positive")
  expect_error(deaths_settings(30, list(recording_tail=list(
    initial=list(weight=1, mean=40)))), "strictly between")
  # the tail is that of the delay from death, whichever kind of case it is
  h <- lag_probabilities_gamma(9.3, 5.5)
  effects <- list(before=list(lag_probabilities=h), tail=list(weight=0.2, mean=30))
  days <- 1 + c(0, 5, 20, 80)
  survivals <- joint_survivals(days, 1, list(mean=7, sd=5), list(mean=3, sd=2.5),
                               NULL, "recording", effects)
  without <- joint_survivals(days, 1, list(mean=7, sd=5), list(mean=3, sd=2.5),
                             NULL, "recording", effects["before"])
  expect_gt(survivals$log_recorded[4], without$log_recorded[4])
  expect_equal(exp(survivals$log_dated),
               0.8 * exp(without$log_dated) + 0.2 * exp(-(days - 0.5) / 30))
  # simulated deaths are recorded later with it
  set.seed(3)
  f_dead <- function(tail) {
    sim <- simulate_cases_deaths(rep(2000, 5), 1:60, list(mean=6, sd=4), 1,
                                 list(mean=3, sd=2.5),
                                 onset_to_death=list(mean=9.3, sd=5.5),
                                 recording_tail=tail)
    sum(sim$snapshots$deaths_reported[sim$snapshots$time_reported == 35])
  }
  expect_lt(f_dead(list(weight=0.3, mean=40)), f_dead(NULL))
})

test_that("the tail of the delay to a death being recorded is sampled", {
  example <- misdated_example()
  snapshots <- example$snapshots
  set.seed(9)
  fit <- mcmc_single(
    niterations=6,
    snapshot_with_Rt_index_df=snapshots %>%
      dplyr::select(time_onset, time_reported, cases_reported) %>%
      dplyr::mutate(Rt_index=1),
    priors=list(Rt=list(shape=1, rate=1),
                reporting=list(mean_mu=6, mean_sigma=5, sd_mu=4, sd_sigma=4),
                max_cases=300, overdispersion=list(mean=20, sd=20),
                Rt_changepoints=list(),
                reporting_changepoints=list(min_length=14)),
    serial_parameters=example$serial_parameters,
    initial_cases_true=example$initial_cases,
    initial_reporting_parameters=list(mean=6, sd=4),
    initial_Rt=dplyr::tibble(Rt_index=1, Rt=1),
    print_to_screen=FALSE, is_negative_binomial=TRUE, initial_overdispersion=10,
    cases_before=example$cases_before,
    deaths=list(
      snapshot=snapshots %>%
        dplyr::select(time_onset, time_reported, deaths_reported),
      settings=list(onset_to_death=list(mean=9.3, sd=5.5), fatal_delay=FALSE,
                    recording_tail=list(),
                    dated_by_death=list(fraction=list()))))
  parameters <- fit$deaths_parameters
  expect_true(all(parameters$recording_tail_weight > 0 &
                    parameters$recording_tail_weight < 1))
  expect_true(all(parameters$recording_tail_mean > 0))
  expect_gt(length(unique(c(parameters$recording_tail_weight,
                            parameters$recording_tail_mean))), 2)
})

test_that("each report can have its own effect on the delay to a death being recorded", {
  expect_false(deaths_settings(30)$recording_report)
  report <- deaths_settings(30, list(recording_report=TRUE))$recording_report
  expect_equal(report$sd_shape, 3)
  expect_equal(report$initial_sd, sqrt(0.2 / 2))
  expect_equal(deaths_settings(30, list(recording_report=list(step=0.3)))$
                 recording_report$step, 0.3)
  expect_error(deaths_settings(30, list(recording_report=list(sd=1))),
               "Unknown settings")
  expect_error(deaths_settings(30, list(recording_report=list(step=0))),
               "must be positive")
  # a report which records more deaths than usual shows more as dead from then
  effects <- update_report_multipliers(list(
    report_days=c(5, 10, 15, 20), report_log_effects=c(0, 0, 0, 0)))
  more <- effects
  more$report_log_effects[3] <- log(4)
  more <- update_report_multipliers(more)
  f <- function(effects) exp(joint_survivals(
    c(5, 10, 15, 20), 1, list(mean=3, sd=2), list(mean=6, sd=4),
    recording_effects=effects)$log_recorded)
  expect_equal(f(more)[1:2], f(effects)[1:2])
  expect_lt(f(more)[3], f(effects)[3])
  # or on the tail of the delay alone, which then mixes with the rest
  expect_error(deaths_settings(30, list(recording_report=list(tail_only=TRUE))),
               "only if it has a tail")
  h <- lag_probabilities_gamma(9.3, 5.5)
  tail_only <- more
  tail_only$report_tail_only <- TRUE
  tail_only <- update_report_multipliers(tail_only)
  expect_null(tail_only$report_multipliers)
  tail_only$tail <- list(weight=0.3, mean=20)
  tail_only$before <- list(lag_probabilities=h)
  days <- c(5, 10, 15, 20, 40)
  g <- function(effects) exp(log_undetected_days(days, 1, list(mean=3, sd=2), effects))
  gamma_part <- g(list(before=tail_only$before))
  tail_part <- g(list(before=tail_only$before, tail=list(weight=1, mean=20),
                      report_multipliers=more$report_multipliers))
  expect_equal(g(tail_only), 0.7 * gamma_part + 0.3 * tail_part)
  expect_lt(tail_part[3], g(list(before=tail_only$before, tail=list(weight=1, mean=20)))[3])
  # and the effects are sampled
  example <- misdated_example()
  snapshots <- example$snapshots
  set.seed(9)
  fit <- mcmc_single(
    niterations=6,
    snapshot_with_Rt_index_df=snapshots %>%
      dplyr::select(time_onset, time_reported, cases_reported) %>%
      dplyr::mutate(Rt_index=1),
    priors=list(Rt=list(shape=1, rate=1),
                reporting=list(mean_mu=6, mean_sigma=5, sd_mu=4, sd_sigma=4),
                max_cases=300, overdispersion=list(mean=20, sd=20),
                Rt_changepoints=list(),
                reporting_changepoints=list(min_length=14)),
    serial_parameters=example$serial_parameters,
    initial_cases_true=example$initial_cases,
    initial_reporting_parameters=list(mean=6, sd=4),
    initial_Rt=dplyr::tibble(Rt_index=1, Rt=1),
    print_to_screen=FALSE, is_negative_binomial=TRUE, initial_overdispersion=10,
    cases_before=example$cases_before,
    deaths=list(
      snapshot=snapshots %>%
        dplyr::select(time_onset, time_reported, deaths_reported),
      settings=list(onset_to_death=list(mean=9.3, sd=5.5), fatal_delay=FALSE,
                    recording_report=TRUE,
                    dated_by_death=list(fraction=list()))))
  parameters <- fit$deaths_parameters
  columns <- grep("^recording_report_[0-9]", names(parameters), value=TRUE)
  expect_equal(length(columns), length(unique(snapshots$time_reported)))
  expect_true(all(parameters$recording_report_sd > 0))
  expect_gt(length(unique(unlist(parameters[columns]))), 1)
  # and the effects on the tail alone are sampled
  set.seed(9)
  fit <- mcmc_single(
    niterations=3,
    snapshot_with_Rt_index_df=snapshots %>%
      dplyr::select(time_onset, time_reported, cases_reported) %>%
      dplyr::mutate(Rt_index=1),
    priors=list(Rt=list(shape=1, rate=1),
                reporting=list(mean_mu=6, mean_sigma=5, sd_mu=4, sd_sigma=4),
                max_cases=300, overdispersion=list(mean=20, sd=20),
                Rt_changepoints=list(),
                reporting_changepoints=list(min_length=14)),
    serial_parameters=example$serial_parameters,
    initial_cases_true=example$initial_cases,
    initial_reporting_parameters=list(mean=6, sd=4),
    initial_Rt=dplyr::tibble(Rt_index=1, Rt=1),
    print_to_screen=FALSE, is_negative_binomial=TRUE, initial_overdispersion=10,
    cases_before=example$cases_before,
    deaths=list(
      snapshot=snapshots %>%
        dplyr::select(time_onset, time_reported, deaths_reported),
      settings=list(onset_to_death=list(mean=9.3, sd=5.5), fatal_delay=FALSE,
                    recording_tail=list(),
                    recording_report=list(tail_only=TRUE),
                    dated_by_death=list(fraction=list()))))
  parameters <- fit$deaths_parameters
  expect_gt(length(unique(unlist(parameters[columns]))), 1)
  expect_true(all(parameters$recording_tail_weight > 0))
})

test_that("cases dated by death which cannot be in the reports are kept", {
  state <- list(
    onset_times=1:2,
    dated=list(into=c(1, 2), out=c(1, 0), exact=FALSE),
    paths=list(
      # the one fatal case reported was first shown as alive
      list(fatal_new=c(1, 0), deaths_new=c(0, 1)),
      # two fatal cases reported, one of them shown as dead at once
      list(fatal_new=c(2, 0), deaths_new=c(1, 0))))
  expect_equal(deaths_cases_offset(state)$cases_offset, c(0 - 1, 1 - 0))
  state$paths <- NULL
  expect_equal(deaths_cases_offset(state)$cases_offset, c(0, 2))
})

test_that("simulated proportions can differ between onset times", {
  set.seed(6)
  cases <- rep(600, 30)
  sim <- simulate_cases_deaths(cases, 1:60, list(mean=6, sd=4), 1,
                               list(mean=3, sd=2.5),
                               onset_to_death=list(mean=9.3, sd=5.5),
                               dated_by_death=list(fraction=rep(c(0, 0.6), each=15)),
                               triggered=rep(c(0.8, 0), each=15))
  # none of the first fifteen onset times are dated by death
  expect_equal(sim$truth$cases_true[1], 600)
  expect_equal(sum(sim$truth$n_dated_by_death[1:16]), 0)
  expect_gt(sum(sim$truth$n_dated_by_death), 0)
  expect_error(simulate_cases_deaths(cases, 1:60, list(mean=6, sd=4), 1,
                                     list(mean=3, sd=2.5), triggered=c(0.1, 0.2)),
               "single value or one for each onset time")
})
