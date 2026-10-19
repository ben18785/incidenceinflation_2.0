backlog_constant <- function(share, rate, n_days=200, report_days=seq_len(n_days)) {
  settings <- backlog_settings(list(initial_share=share, initial_rate=rate))
  initial_backlog(n_days, report_days, settings)
}

test_that("the settings of the backlog are checked", {
  expect_false(deaths_settings(30)$recording_backlog)
  backlog <- deaths_settings(30, list(recording_backlog=TRUE))$recording_backlog
  expect_equal(backlog$share_week, 7)
  backlog <- deaths_settings(30, list(recording_backlog=list(rate_sd=2)))$recording_backlog
  expect_equal(backlog$rate_sd, 2)
  expect_error(deaths_settings(30, list(recording_backlog=list(rate=1))), "Unknown settings")
  expect_error(deaths_settings(30, list(recording_backlog=list(share_step=0))),
               "must be positive")
  expect_error(deaths_settings(30, list(recording_backlog=list(initial_share=1))),
               "strictly between")
  expect_error(deaths_settings(30, list(recording_backlog=TRUE, recording_tail=list())),
               "takes the place")
  expect_error(deaths_settings(30, list(recording_backlog=TRUE, recording_report=TRUE)),
               "takes the place")
})

test_that("the hazard of the backlog is that of the report interval of each day", {
  backlog <- backlog_constant(0.2, 0.1, n_days=20, report_days=c(5, 8, 12))
  backlog$log_rate <- log(c(0.1, 0.2, 0.4))
  backlog <- refresh_backlog(backlog)
  expect_equal(backlog_daily_rates(backlog), c(rep(0.1, 5), rep(0.2, 3), rep(0.4, 4)))
  expect_equal(backlog_cumulative(backlog, c(0, 5, 12, 14)),
               c(0, 0.5, 0.5 + 0.6 + 1.6, 2.7 + 0.8))
  expect_equal(backlog$share, rep(0.2, 20))
})

test_that("a backlog with a constant hazard is the exponential tail", {
  h <- lag_probabilities_gamma(9.3, 5.5)
  days <- 3 + c(0, 4, 15, 40, 100)
  backlog <- backlog_constant(0.3, 1 / 40)
  with_backlog <- log_undetected_days(days, 3, list(mean=3, sd=2.5),
                                      list(before=list(lag_probabilities=h),
                                           backlog=backlog))
  with_tail <- log_undetected_days(days, 3, list(mean=3, sd=2.5),
                                   list(before=list(lag_probabilities=h),
                                        tail=list(weight=0.3, mean=40)))
  expect_equal(with_backlog, with_tail)
  # and so for a case dated by death, whose delay runs from its day
  without_before <- log_undetected_days(days, 3, list(mean=3, sd=2.5),
                                        list(backlog=backlog))
  expect_equal(exp(without_before),
               0.7 * exp(log_undetected_days(days, 3, list(mean=3, sd=2.5))) +
                 0.3 * exp(-(days - 3 + 0.5) / 40))
  # a share of zero leaves the delay as it was
  backlog$share[] <- 0
  expect_equal(log_undetected_days(days, 3, list(mean=3, sd=2.5),
                                   list(before=list(lag_probabilities=h),
                                        backlog=backlog)),
               log_undetected_days(days, 3, list(mean=3, sd=2.5),
                                   list(before=list(lag_probabilities=h))))
})

test_that("the hazard of the backlog acts on calendar days", {
  # deaths of two onset days wait for the same days of clearing
  backlog <- backlog_constant(0.5, 1e-4, n_days=60, report_days=c(20, 30, 60))
  backlog$log_rate <- log(c(1e-4, 0.3, 1e-4))
  backlog <- refresh_backlog(backlog)
  f <- function(onset) exp(log_backlog_days(c(20, 30, 60), onset, backlog))
  expect_equal(f(5), f(10), tolerance=1e-3)
  expect_equal(f(5)[2], exp(-3), tolerance=1e-2)
})

test_that("simulated deaths in the backlog are recorded with the hazard given", {
  set.seed(11)
  rate <- c(rep(0.002, 30), rep(0.1, 20), rep(0.02, 50))
  sim <- simulate_cases_deaths(c(20000, rep(0, 9)), 1:100, list(mean=0.2, sd=0.1),
                               1, list(mean=1, sd=0.5),
                               onset_to_death=list(mean=9.3, sd=5.5),
                               recording_backlog=list(share=0.4, rate=rate))
  shown <- sim$snapshots %>% dplyr::filter(time_onset == 1)
  backlog <- backlog_constant(0.4, 0.002, n_days=100, report_days=1:100)
  backlog$log_rate <- log(rate)
  backlog <- refresh_backlog(backlog)
  expected <- exp(log_undetected_days(shown$time_reported, 1, list(mean=1, sd=0.5),
                                      list(before=list(lag_probabilities=
                                                         lag_probabilities_gamma(9.3, 5.5)),
                                           backlog=backlog)))
  observed <- 1 - shown$deaths_reported / sim$truth$deaths_true[1]
  is_checked <- shown$time_reported %in% c(20, 35, 50, 80)
  expect_equal(observed[is_checked], expected[is_checked], tolerance=0.03)
  expect_error(simulate_cases_deaths(rep(10, 5), 1:10, list(mean=6, sd=4), 1,
                                     list(mean=3, sd=2.5),
                                     recording_backlog=list(share=0.2, rate=0.1)),
               "time from onset to death")
})

backlog_example <- function() {
  set.seed(41)
  n_days <- 30
  serial_parameters <- list(mean=5, sd=3)
  cases <- true_cases(n_days + 40, function(t) 1, kappa=1000, serial_parameters,
                      initial_parameters=list(mean=30, length=10))
  sim <- simulate_cases_deaths(cases[41:(n_days + 40)], seq_len(n_days + 10),
                               list(mean=6, sd=4), 0.5, list(mean=3, sd=2.5),
                               triggered=0.5,
                               onset_to_death=list(mean=9.3, sd=5.5),
                               dated_by_death=list(fraction=0.3),
                               recording_backlog=list(share=0.3, rate=0.05))
  list(snapshots=sim$snapshots, cases_before=cases[1:40],
       serial_parameters=serial_parameters,
       initial_cases=sim$snapshots %>%
         dplyr::group_by(time_onset) %>%
         dplyr::summarise(cases_true=max(cases_reported) + 2))
}

for(sampler in c("slice", "metropolis")) test_that(paste("the backlog is sampled by", sampler), {
  example <- backlog_example()
  snapshots <- example$snapshots
  set.seed(13)
  fit <- mcmc_single(
    niterations=5,
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
                    triggered=list(), recording_backlog=list(sampler=sampler),
                    dated_by_death=list(fraction=list()))))
  parameters <- fit$deaths_parameters
  rates <- parameters[, grepl("^backlog_rate_[0-9]+$", names(parameters))]
  expect_equal(ncol(rates), length(unique(snapshots$time_reported)))
  expect_true(all(as.matrix(rates) > 0))
  expect_gt(length(unique(rates[[1]])), 1)
  share <- fit$deaths$backlog_share
  expect_true(all(share > 0 & share < 1))
  expect_gt(length(unique(share)), 1)
  expect_true(all(parameters$backlog_share_rw_sd > 0))
})

test_that("the backlog is sampled with the exact form of fatal cases dated by death", {
  example <- backlog_example()
  snapshots <- example$snapshots
  set.seed(14)
  fit <- mcmc_single(
    niterations=5,
    snapshot_with_Rt_index_df=snapshots %>%
      dplyr::select(time_onset, time_reported, cases_reported) %>%
      dplyr::mutate(Rt_index=1),
    priors=list(Rt=list(shape=1, rate=1),
                reporting=list(mean_mu=6, mean_sigma=5, sd_mu=4, sd_sigma=4),
                max_cases=300, Rt_changepoints=list(),
                reporting_changepoints=list(min_length=14)),
    serial_parameters=example$serial_parameters,
    initial_cases_true=example$initial_cases,
    initial_reporting_parameters=list(mean=6, sd=4),
    initial_Rt=dplyr::tibble(Rt_index=1, Rt=1),
    print_to_screen=FALSE,
    cases_before=example$cases_before,
    deaths=list(
      snapshot=snapshots %>%
        dplyr::select(time_onset, time_reported, deaths_reported),
      settings=list(onset_to_death=list(mean=9.3, sd=5.5), fatal_delay=FALSE,
                    triggered=list(), recording_backlog=list(hold_days=5),
                    dated_by_death=list(fraction=list(), exact=TRUE))))
  expect_equal(nrow(fit$deaths), 5 * 30)
  both <- dplyr::left_join(fit$deaths, fit$cases, by=c("iteration", "time_onset"))
  expect_true(all(both$deaths_true <= both$cases_true))
  expect_true(all(fit$deaths$dated_by_death <= fit$deaths$deaths_recorded))
  expect_true(all(is.finite(fit$deaths_parameters$triggered)))
  # without negative binomial renewal, that term is missing
  expect_equal(fit$log_likelihood$iteration, 1:5)
  expect_true(all(is.na(fit$log_likelihood$renewal)))
  expect_true(all(is.finite(fit$log_likelihood$reports)))
})

test_that("the terms of the log-likelihood are those of the state at the end of each iteration", {
  example <- backlog_example()
  snapshots <- example$snapshots
  set.seed(15)
  fit <- mcmc_single(
    niterations=4,
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
    cases_before=example$cases_before, return_state=TRUE,
    deaths=list(
      snapshot=snapshots %>%
        dplyr::select(time_onset, time_reported, deaths_reported),
      settings=list(onset_to_death=list(mean=9.3, sd=5.5), fatal_delay=FALSE,
                    triggered=list(), recording_backlog=list(hold_days=5),
                    dated_by_death=list(fraction=list(), exact=TRUE))))
  logl <- fit$log_likelihood
  expect_equal(logl$iteration, 1:4)
  expect_true(all(is.finite(as.matrix(logl[, c("reports", "recording", "renewal", "fatal")]))))
  v <- fit$state$variables
  last <- logl[4, ]
  terms <- deaths_log_likelihood(v$deaths_state, v$deaths_cases_true)
  expect_equal(last$reports, sum(v$reporting_state$logl))
  expect_equal(last$recording, terms[["recording"]])
  expect_equal(last$fatal, terms[["fatal"]])
  expect_equal(last$renewal, v$logp_overdispersion - v$logp_prior)
  # carried on, the iterations are numbered on
  more <- mcmc_single(
    niterations=2,
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
    cases_before=example$cases_before, resume_from=fit$state,
    deaths=list(
      snapshot=snapshots %>%
        dplyr::select(time_onset, time_reported, deaths_reported),
      settings=list(onset_to_death=list(mean=9.3, sd=5.5), fatal_delay=FALSE,
                    triggered=list(), recording_backlog=list(hold_days=5),
                    dated_by_death=list(fraction=list(), exact=TRUE))))
  expect_equal(more$log_likelihood$iteration, 5:6)
})

test_that("the probability a death in the backlog is not recorded is at most one", {
  backlog <- backlog_constant(0.3, 0.05, n_days=60, report_days=1:60)
  log_s <- log_undetected_days(1:40, 1, list(mean=3, sd=2.5),
                               list(before=list(lag_probabilities=
                                                  lag_probabilities_gamma(9.3, 5.5)),
                                    backlog=backlog))
  expect_true(all(log_s <= 0))
  expect_equal(log_s[1], 0)
})

test_that("elliptical slice sampling with a normal likelihood draws from the posterior", {
  set.seed(5)
  # prior N(1, 2^2), likelihood N(x; 3, 1): posterior N(2.6, 0.8)
  x <- 0
  logl <- stats::dnorm(3, x, 1, log=TRUE)
  draws <- numeric(4000)
  for(i in seq_along(draws)) {
    out <- elliptical_slice(x, 1, function() stats::rnorm(1, 0, 2),
                            function(x) stats::dnorm(3, x, 1, log=TRUE), logl)
    x <- out$x
    logl <- out$logl
    draws[i] <- x
  }
  expect_equal(mean(draws), 2.6, tolerance=0.05)
  expect_equal(stats::var(draws), 0.8, tolerance=0.1)
})

test_that("the backlog can start from a share for each onset day and a rate for each day", {
  settings <- backlog_settings(list(initial_share=rep(c(0.4, 0.1), c(14, 6)),
                                    initial_rate=c(rep(0.01, 10), rep(0.1, 10))))
  backlog <- initial_backlog(20, c(10, 15, 20), settings)
  expect_equal(stats::plogis(backlog$logit_share), c(0.4, 0.4, 0.1))
  expect_equal(exp(backlog$log_rate), c(0.01, 0.1, 0.1))
  expect_error(backlog_settings(list(sampler="gibbs")), "must be 'slice'")
})

test_that("with a holding period, deaths in the backlog are recorded only after it", {
  set.seed(12)
  rate <- c(rep(0.01, 30), rep(0.1, 20), rep(0.02, 50))
  backlog_truth <- list(share=0.4, rate=rate, hold_days=14)
  sim <- simulate_cases_deaths(c(20000, rep(0, 9)), 1:100, list(mean=0.2, sd=0.1),
                               1, list(mean=1, sd=0.5),
                               onset_to_death=list(mean=9.3, sd=5.5),
                               recording_backlog=backlog_truth)
  shown <- sim$snapshots %>% dplyr::filter(time_onset == 1)
  settings <- backlog_settings(list(initial_share=0.4, hold_days=14))
  backlog <- initial_backlog(100, 1:100, settings)
  backlog$log_rate <- log(rate)
  backlog <- refresh_backlog(backlog)
  expected <- exp(log_undetected_days(shown$time_reported, 1, list(mean=1, sd=0.5),
                                      list(before=list(lag_probabilities=
                                                         lag_probabilities_gamma(9.3, 5.5)),
                                           backlog=backlog)))
  observed <- 1 - shown$deaths_reported / sim$truth$deaths_true[1]
  is_checked <- shown$time_reported %in% c(10, 20, 30, 45, 60, 90)
  expect_equal(observed[is_checked], expected[is_checked], tolerance=0.03)
  # nothing of the backlog is recorded within the holding period: by day 14,
  # the deaths not recorded are the backlog and those of the gamma delay
  # still to come
  expect_gt(expected[shown$time_reported == 14], 0.4)
  expect_error(backlog_settings(list(hold_days=-1)), "whole number")
})
