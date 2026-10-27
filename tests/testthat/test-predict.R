predict_example <- function(n_reports=40) {
  set.seed(51)
  n_days <- 30
  serial_parameters <- list(mean=5, sd=3)
  cases <- true_cases(n_days + 40, function(t) 1, kappa=1000, serial_parameters,
                      initial_parameters=list(mean=30, length=10))
  sim <- simulate_cases_deaths(cases[41:(n_days + 40)], seq_len(n_days + 30),
                               list(mean=6, sd=4), 0.5, list(mean=3, sd=2.5),
                               triggered=1,
                               onset_to_death=list(mean=9.3, sd=5.5),
                               dated_by_death=list(fraction=0.3),
                               recording_backlog=list(share=0.3, rate=0.05),
                               wrong=list(eta=0.02, rho=0.1, delta=0.3))
  snapshots <- dplyr::filter(sim$snapshots, time_reported <= n_reports)
  list(snapshots=snapshots, later=sim$snapshots, cases_before=cases[1:40],
       serial_parameters=serial_parameters,
       initial_cases=snapshots %>%
         dplyr::group_by(time_onset) %>%
         dplyr::summarise(cases_true=max(cases_reported) + 2))
}

predict_example_fit <- function(example, niterations=3) {
  snapshots <- example$snapshots
  mcmc_single(
    niterations=niterations,
    snapshot_with_Rt_index_df=snapshots %>%
      dplyr::select(time_onset, time_reported, cases_reported) %>%
      dplyr::mutate(Rt_index=1),
    priors=list(Rt=list(shape=1, rate=1),
                reporting=list(mean_mu=6, mean_sigma=5, sd_mu=4, sd_sigma=4),
                max_cases=300, overdispersion=list(mean=20, sd=20),
                revisions=list(eta_mu=0.1, eta_sigma=0.1, rho_alpha=1, rho_beta=1),
                max_erroneous=4,
                Rt_changepoints=list(),
                reporting_changepoints=list(min_length=14)),
    serial_parameters=example$serial_parameters,
    initial_cases_true=example$initial_cases,
    initial_reporting_parameters=list(mean=6, sd=4),
    initial_Rt=dplyr::tibble(Rt_index=1, Rt=1),
    initial_revision_parameters=list(eta=0.02, rho=0.1, weekday=TRUE),
    print_to_screen=FALSE, is_negative_binomial=TRUE, initial_overdispersion=10,
    cases_before=example$cases_before, return_state=TRUE,
    deaths=list(
      snapshot=snapshots %>%
        dplyr::select(time_onset, time_reported, deaths_reported),
      settings=list(onset_to_death=list(mean=9.3, sd=5.5), fatal_delay=FALSE,
                    triggered=list(initial=1, n_steps=0L),
                    recording_backlog=list(hold_days=5),
                    dated_by_death=list(fraction=list(), exact=TRUE))))
}

predict_draw <- function(fit, iteration) {
  lapply(fit[c("cases", "reporting_days", "revisions", "reporting_effects",
               "deaths", "deaths_parameters")],
         function(df) df[df$iteration == iteration, ])
}

test_that("the state at a kept draw is that of the sampler at that draw", {
  example <- predict_example()
  set.seed(52)
  fit <- predict_example_fit(example)
  v <- fit$state$variables
  state <- v$deaths_state
  prediction <- predict_state_from_draw(v, predict_draw(fit, 3))
  rebuilt <- prediction$deaths_state
  expect_equal(prediction$cases_true, v$deaths_cases_true)
  expect_equal(rebuilt$deaths_true, state$deaths_true)
  expect_equal(rebuilt$dated$out, state$dated$out)
  expect_equal(rebuilt$dated$into, state$dated$into)
  expect_equal(rebuilt$delays$mean, deaths_delay_by_onset(state)$mean)
  expect_equal(rebuilt$recording_effects$backlog$log_rate,
               state$recording_effects$backlog$log_rate)
  expect_equal(rebuilt$recording_effects$backlog$cumulative,
               state$recording_effects$backlog$cumulative)
  expect_equal(prediction$revision_parameters$report_multipliers,
               v$revision_current$report_multipliers)
  expect_equal(rebuilt$recording_effects$report_multipliers,
               state$recording_effects$report_multipliers)
  # the delay to a death being recorded acts on the reports as in the sampler
  for(k in seq_along(state$onset_times)) {
    terms <- predict_onset_terms(prediction, k, state$reports[[k]]$time_reported)
    expect_equal(terms$survivals$log_recorded, state$survivals[[k]]$log_recorded)
    expect_equal(terms$survivals$log_dated, state$survivals[[k]]$log_dated)
  }
  # a draw from another chain's sampler cannot be used without a state
  expect_error(predict_state_from_draw(list(is_deaths=FALSE), predict_draw(fit, 3)),
               "model of deaths")
})

test_that("later reports are simulated from the last report", {
  example <- predict_example()
  set.seed(53)
  fit <- predict_example_fit(example)
  v <- fit$state$variables
  draw <- predict_draw(fit, 3)
  out <- predict_reports(v, draw, 38:45)
  # only days after the last report are simulated
  expect_equal(sort(unique(out$time_reported)), 41:45)
  expect_equal(sort(unique(out$time_onset)), 1:30)
  expect_true(all(out$deaths_reported <= out$cases_reported))
  expect_true(all(out$deaths_reported >= 0))
  # a later report differs from the last by the cases and deaths still to
  # come and the wrong records, of which there are at most max_erroneous
  last <- example$snapshots %>% dplyr::filter(time_reported == 40)
  first <- dplyr::left_join(dplyr::filter(out, time_reported == 41), last,
                            by="time_onset", suffix=c("", "_last"))
  expect_true(all(first$cases_reported >= first$cases_reported_last - 4))
  # far ahead, every case recorded for a day has been reported and every
  # death recorded: what is left are the wrong records
  prediction <- predict_state_from_draw(v, draw)
  dated <- dated_counts(prediction$deaths_state)
  far <- predict_reports(v, draw, 2000, onset_times=c(5, 20))
  recorded <- prediction$cases_true - dated$out + dated$into
  extra <- far$cases_reported - recorded[c(5, 20)]
  expect_true(all(extra >= 0 & extra <= 4))
  extra_dead <- far$deaths_reported - dated$fatal_recorded[c(5, 20)]
  expect_true(all(extra_dead >= 0 & extra_dead <= extra))
  expect_error(predict_reports(v, draw, 30:40), "no later reports")
})

test_that("the effects of later reports are those of their days of the week", {
  example <- predict_example()
  set.seed(54)
  fit <- predict_example_fit(example)
  v <- fit$state$variables
  prediction <- predict_state_from_draw(v, predict_draw(fit, 3))
  extended <- predict_extend_reports(prediction, c(43, 47))
  revision <- extended$revision_parameters
  multipliers <- revision$report_multipliers
  # days up to the last report are as they were
  n_before <- length(prediction$revision_parameters$report_multipliers$values)
  expect_equal(multipliers$values[seq_len(n_before)],
               prediction$revision_parameters$report_multipliers$values)
  # days 41 to 43 take the effect of the report of day 43, and 44 to 47 that
  # of day 47
  effect <- function(day) exp(revision$weekday_log_effects[[
    as.character(weekday_of_days(day, revision$weekday_first))]])
  expect_equal(multipliers$values[n_before + 1:3], rep(effect(43), 3))
  expect_equal(multipliers$values[n_before + 4:7], rep(effect(47), 4))
  # the backlog's hazard is held between the new reports and goes on from the
  # last
  backlog <- extended$deaths_state$recording_effects$backlog
  before <- prediction$deaths_state$recording_effects$backlog
  expect_equal(backlog$report_days, c(before$report_days, 43, 47))
  expect_equal(backlog$log_rate[seq_along(before$log_rate)], before$log_rate)
  expect_equal(backlog$cumulative[seq_along(before$cumulative)], before$cumulative)
  rates <- backlog_daily_rates(backlog)
  expect_equal(length(unique(rates[41:43])), 1)
  expect_equal(length(unique(rates[44:47])), 1)
})

test_that("the simulation of later reports follows the model's transitions", {
  # one report so far: 2 survivors, 1 fatal case dated by onset (shown as
  # alive), 1 dated by death, and one wrong record shown as alive
  path <- list(survivors_before=0, survivors_new=2, fatal_reported=1,
               dated_reported=1, deaths_new=0, wrong_alive=1, wrong_dead=0)
  totals <- list(survivors=5, fatal=3, dated=2)
  # everything happens in the next report, and nothing is entered
  probabilities <- list(p=c(0.5, 1), p_fatal=c(0.5, 1), p_dated=c(0.5, 1),
                        g=c(0, 1), q=c(0, 1))
  wrong <- list(lambda=c(0, 0), w=c(0, 1), delta=0.5, max_erroneous=2)
  out <- predict_onset_simulate(path, probabilities, wrong, totals, 1)
  expect_equal(unname(out[1, ]), c(5 + 3 + 2, 3 + 2))
  # nothing happens
  probabilities <- list(p=c(0.5, 0), p_fatal=c(0.5, 0), p_dated=c(0.5, 0),
                        g=c(0, 0), q=c(0, 0))
  wrong$w <- c(0, 0)
  out <- predict_onset_simulate(path, probabilities, wrong, totals, 1)
  expect_equal(unname(out[1, ]), c(2 + 1 + 1 + 1, 1))
  # with no room for more wrong records, none can be entered
  wrong <- list(lambda=c(0, 50), w=c(0, 0), delta=0.5, max_erroneous=1)
  expect_error(predict_onset_simulate(path, probabilities, wrong, totals, 1,
                                      max_tries=5), "max_erroneous")
  # survivors newly reported are binomial
  set.seed(55)
  probabilities <- list(p=c(0.5, 0.3), p_fatal=c(0.5, 0), p_dated=c(0.5, 0),
                        g=c(0, 0), q=c(0, 0))
  wrong <- list(lambda=c(0, 0), w=c(0, 0), delta=0.5, max_erroneous=2)
  draws <- replicate(4000, predict_onset_simulate(path, probabilities, wrong,
                                                  totals, 1)[1, 1])
  expect_equal(mean(draws - 5), 3 * 0.3, tolerance=0.05)
})
