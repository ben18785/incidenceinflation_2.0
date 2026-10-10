simulate_corrected <- function(seed=1, n_days=40, first_report=20, time_corrected=32,
                               corrected=1:15) {
  set.seed(seed)
  cases_true <- stats::rpois(n_days, 25)
  simulated <- simulate_cases_deaths(
    cases_true, first_report:n_days, list(mean=5, sd=3), 0.4, list(mean=9, sd=5),
    correction=list(time_reported=time_corrected, onset_times=corrected,
                    wrong=0.1, missing=0.03, misdated=0.1))
  c(simulated, list(cases_true=cases_true, time_corrected=time_corrected,
                    corrected=corrected))
}

test_that("a correcting report changes the counts of the onset times it corrects only", {
  set.seed(3)
  cases_true <- stats::rpois(40, 25)
  wrong <- simulate_cases_deaths(
    cases_true, 20:40, list(mean=5, sd=3), 0.4, list(mean=9, sd=5),
    correction=list(time_reported=32, onset_times=1:15, wrong=0.2))
  expect_equal(wrong$truth$cases_true, cases_true)
  # records which should not be there are present until the correction: long
  # after onset the counts exceed the true cases before it and equal them after
  f_counts <- function(day) {
    rows <- wrong$snapshots[wrong$snapshots$time_reported == day, ]
    rows$cases_reported[order(rows$time_onset)]
  }
  expect_true(all(f_counts(31)[1:15] >= cases_true[1:15]))
  expect_gt(sum(f_counts(31)[1:15]), sum(cases_true[1:15]))
  expect_equal(f_counts(40)[1:15], cases_true[1:15])
  expect_true(all(f_counts(31)[16:25] <= cases_true[16:25]))

  simulated <- simulate_corrected()
  snapshots <- simulated$snapshots
  expect_true(all(snapshots$deaths_reported <= snapshots$cases_reported))
  # long after onset, the corrected reports hold the true cases
  latest <- snapshots[snapshots$time_reported == 40 & snapshots$time_onset <= 15, ]
  expect_equal(latest$cases_reported, simulated$cases_true[1:15])
  before <- snapshots[snapshots$time_reported == 31 & snapshots$time_onset <= 15, ]
  expect_false(all(before$cases_reported == simulated$cases_true[1:15]))
})

test_that("split_corrected_reports splits the reports of the onset times corrected", {
  simulated <- simulate_corrected()
  parts <- split_corrected_reports(simulated$snapshots, 32, 1:15)
  expect_equal(nrow(parts$corrected) + nrow(parts$uncorrected),
               nrow(simulated$snapshots))
  expect_true(all(parts$uncorrected$time_onset <= 15 &
                    parts$uncorrected$time_reported < 32))
  expect_equal(min(parts$corrected$time_reported[parts$corrected$time_onset <= 15]), 32)
  expect_equal(min(parts$corrected$time_reported[parts$corrected$time_onset == 16]), 20)
})

fit_corrected <- function(simulated, niterations, is_uncorrected=TRUE) {
  parts <- split_corrected_reports(simulated$snapshots, simulated$time_corrected,
                                   simulated$corrected)
  snapshots <- parts$corrected
  largest <- snapshots %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(cases_true=max(cases_reported) + 2)
  mcmc_single(
    niterations=niterations,
    snapshot_with_Rt_index_df=snapshots %>%
      dplyr::select(time_onset, time_reported, cases_reported) %>%
      dplyr::mutate(Rt_index=1),
    priors=list(Rt=list(shape=1, rate=1),
                reporting=list(mean_mu=5, mean_sigma=10, sd_mu=3, sd_sigma=5),
                max_cases=150, Rt_changepoints=list(),
                reporting_changepoints=list(min_length=14, n_partition_moves=5),
                revisions=list(eta_mu=0.01, eta_sigma=0.02, rho_alpha=2, rho_beta=4),
                max_erroneous=4),
    serial_parameters=list(mean=5, sd=3),
    initial_cases_true=largest,
    initial_reporting_parameters=list(mean=5, sd=3),
    initial_Rt=dplyr::tibble(Rt_index=1, Rt=1),
    reporting_metropolis_parameters=list(mean_step=0.3, sd_step=0.25),
    print_to_screen=FALSE,
    initial_revision_parameters=list(eta=0.01, rho=0.3, weekday=TRUE, weekday_first=1),
    deaths=c(list(snapshot=snapshots %>%
                    dplyr::select(time_onset, time_reported, deaths_reported),
                  settings=list(delay_changepoints=list(min_length=14,
                                                        n_partition_moves=5))),
             if(is_uncorrected) list(uncorrected=parts$uncorrected)))
}

test_that("the model runs with uncorrected series and returns them", {
  simulated <- simulate_corrected(seed=2)
  set.seed(5)
  fit <- fit_corrected(simulated, 6)
  uncorrected <- fit$deaths_uncorrected
  expect_equal(nrow(uncorrected), 6 * 15)
  expect_true(all(uncorrected$deaths_uncorrected <= uncorrected$cases_uncorrected))
  # the uncorrected series hold at least the cases reported before the
  # correction, less the wrong records allowed for
  largest <- simulated$snapshots %>%
    dplyr::filter(time_onset <= 15, time_reported < 32) %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(largest=max(cases_reported))
  joined <- dplyr::left_join(uncorrected, largest, by="time_onset")
  expect_true(all(joined$cases_uncorrected >= joined$largest - 4))
  expect_true(all(is.finite(fit$deaths$cfr)))
  expect_equal(nrow(fit$deaths), 6 * 40)
})

test_that("the uncorrected series enter what the delays are updated with", {
  simulated <- simulate_corrected(seed=3)
  parts <- split_corrected_reports(simulated$snapshots, 32, 1:15)
  settings <- deaths_settings(40)
  largest <- parts$corrected %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(cases_true=max(cases_reported) + 2)
  state <- initial_deaths_state(parts$corrected, largest, settings)
  state$uncorrected <- initial_uncorrected(parts$uncorrected, state, settings)
  reporting <- dplyr::tibble(mean=rep(5, 40), sd=3)
  set.seed(6)
  state <- deaths_paths_step(state, largest$cases_true, reporting)
  state <- deaths_uncorrected_step(state, reporting, max_cases=150)
  without <- state
  without$uncorrected <- NULL
  sub <- uncorrected_state(state)
  expect_equal(sub$delays, deaths_delay_by_onset(state)[1:15, ])

  total <- deaths_reporting_logp(state, largest$cases_true, reporting)
  expect_equal(total,
               deaths_reporting_logp(without, largest$cases_true, reporting) +
                 deaths_reporting_logp(sub, sub$cases_true, reporting[1:15, ]))
  f_with <- deaths_reporting_logp_function(state, largest$cases_true)
  f_without <- deaths_reporting_logp_function(without, largest$cases_true)
  expect_equal(f_with(1:40, list(mean=5, sd=3)), total)
  # onset times which were not corrected are not affected
  expect_equal(f_with(16:40, list(mean=6, sd=3)), f_without(16:40, list(mean=6, sd=3)))
  expect_false(f_with(1:15, list(mean=6, sd=3)) == f_without(1:15, list(mean=6, sd=3)))
  g_with <- deaths_recording_logp_function(state)
  g_without <- deaths_recording_logp_function(without)
  expect_equal(g_with(16:40, list(mean=9, sd=5)), g_without(16:40, list(mean=9, sd=5)))
  expect_equal(g_with(3, list(mean=9, sd=5)),
               g_without(3, list(mean=9, sd=5)) +
                 deaths_recording_day_logp(sub, 3, list(mean=9, sd=5)))
})

test_that("the cases of an uncorrected series are not tied to the case fatality ratio", {
  # with every case long since reported, the number of fatal cases of an
  # uncorrected series is fixed by its reports whatever the ratio elsewhere
  simulated <- simulate_corrected(seed=4)
  parts <- split_corrected_reports(simulated$snapshots, 32, 1:15)
  settings <- deaths_settings(40)
  largest <- parts$corrected %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(cases_true=max(cases_reported) + 2)
  state <- initial_deaths_state(parts$corrected, largest, settings)
  state$uncorrected <- initial_uncorrected(parts$uncorrected, state, settings)
  reporting <- dplyr::tibble(mean=rep(5, 40), sd=3)
  set.seed(7)
  for(i in 1:5)
    state <- deaths_uncorrected_step(state, reporting, max_cases=150)
  expect_true(all(state$uncorrected$cfr > 0 & state$uncorrected$cfr < 1))
  expect_false(all(state$uncorrected$cfr == settings$initial_cfr))
  reported <- purrr::map_dbl(state$uncorrected$reports, ~max(.$cases_reported))
  expect_true(all(state$uncorrected$cases_true >= reported))
  expect_true(all(state$uncorrected$cases_true <= 150))
})
