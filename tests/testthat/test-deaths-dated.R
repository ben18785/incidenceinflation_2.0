# a state of the model of deaths in which some fatal cases are dated by death
dated_test_state <- function(reporting=NULL, cases_before=NULL, n_days=25) {
  cases <- rep(40, n_days)
  sim <- simulate_cases_deaths(
    cases, seq_len(n_days + 5), list(mean=5, sd=3), 0.4, list(mean=9, sd=4),
    dated_by_death=list(fraction=0.4, mean=5, sd=3, reporting=reporting))
  settings <- deaths_settings(n_days, list(dated_by_death=list(
    fraction=0.4, mean=5, sd=3, reporting=reporting)))
  initial_cases <- sim$snapshots %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(cases_true=max(cases_reported) + 5)
  state <- initial_deaths_state(sim$snapshots, initial_cases, settings, cases_before)
  reporting_by_onset <- dplyr::tibble(mean=rep(5, n_days), sd=rep(3, n_days))
  cases_true <- state$cases_onset
  state <- deaths_paths_step(state, cases_true, reporting_by_onset)
  list(state=state, cases_true=cases_true, sim=sim, settings=settings,
       recorded=initial_cases$cases_true,
       renewal=dated_renewal(rep(1, n_days), 5, TRUE, list(mean=5, sd=3), 20,
                             cases_before))
}

# the counts which follow from the fatal cases dated by death
expect_dated_consistent <- function(state, cases_true, recorded) {
  dated <- state$dated
  n_days <- length(state$onset_times)
  expect_equal(dated$out, rowSums(dated$lags) + dated$unseen)
  into <- sapply(seq_len(n_days), function(t) {
    z <- seq_len(min(dated$max_lag, t - 1))
    sum(dated$lags[cbind(t - z, z)])
  }) + dated$inflow
  expect_equal(dated$into, into)
  expect_true(all(dated$lags >= 0) & all(dated$unseen >= 0) & all(dated$inflow >= 0))
  # fatal cases dated by onset, and cases which are not fatal
  expect_true(all(state$deaths_true - dated$out >= 0))
  expect_true(all(cases_true - state$deaths_true >= 0))
  # the cases recorded for each day are not altered
  expect_equal(cases_true - dated$out + dated$into, recorded)
  # no case can die on a day for which there are no reports
  days <- row(dated$lags) + col(dated$lags)
  expect_true(all(dated$lags[days > n_days] == 0))
}

test_that("initial_dated_state gives counts by day of onset which are consistent", {
  set.seed(11)
  test <- dated_test_state()
  dated <- test$state$dated
  expect_dated_consistent(test$state, test$cases_true, test$recorded)
  expect_true(sum(dated$lags) > 0)
  # no earlier cases are known, so none can arrive from before the first day
  expect_equal(dated$inflow_weight, rep(0, 25))
  expect_equal(sum(dated$inflow), 0)
  expect_equal(dated$into[1], 0)
  # deaths after the last day are possible only for the latest onsets
  expect_equal(dated$tail[1], 0)
  expect_equal(dated$tail[25], 1)
  expect_equal(dated$tail[24], 1 - dated$lag_probabilities[1])
  expect_true(is.finite(dated_logp(test$state, test$cases_true)))
  # about the fraction supplied of the fatal cases recorded are dated by death
  counts <- dated_counts(test$state)
  expect_true(abs(sum(counts$into[5:25]) / sum(counts$fatal_recorded[5:25]) - 0.4) < 0.05)
  # before the reports have been gone through, every case dated by death is
  # taken to be in them; afterwards, no more than could have been shown as
  # dead in the first report they were in
  without_paths <- test$state
  without_paths["paths"] <- list(NULL)
  expect_equal(deaths_cases_offset(without_paths)$cases_offset,
               dated$into - dated$out)
  offset <- deaths_cases_offset(test$state)
  expect_true(all(offset$cases_offset <= dated$into - dated$out))
  expect_true(all(offset$cases_offset >= -dated$out))

  # with earlier cases, some recorded for the first days had onset before them
  set.seed(12)
  with_before <- dated_test_state(cases_before=rep(40, 10))
  dated <- with_before$state$dated
  expect_true(all(dated$inflow_weight[1:3] > 0))
  # ten days of earlier cases are known
  expect_equal(dated$inflow_weight[1], 40 * sum(dated$lag_probabilities[1:10]))
  expect_equal(dated$inflow_weight[2], 40 * sum(dated$lag_probabilities[2:11]))
  expect_equal(dated$inflow_weight[dated$max_lag + 1], 0)
  expect_dated_consistent(with_before$state, with_before$cases_true,
                          with_before$recorded)

  # without such cases nothing is held
  state <- initial_deaths_state(
    test$sim$snapshots,
    dplyr::tibble(time_onset=1:25, cases_true=test$recorded),
    deaths_settings(25))
  expect_null(state$dated)
  expect_null(deaths_cases_offset(state))
  expect_equal(dated_counts(state)$into, rep(0, 25))
  expect_equal(state$cases_onset, test$recorded)
})

test_that("dated_renewal_logp is that of the renewal model", {
  cases <- c(5, 8, 3, 0, 6, 9)
  serial <- list(mean=3, sd=2)
  for(cases_before in list(NULL, c(4, 7, 2))) {
    renewal <- dated_renewal(rep(1.3, 6), 4, TRUE, serial, 5, cases_before)
    statistics <- renewal_statistics(cases, serial, 5, cases_before)
    expect_equal(as.vector(renewal$weights %*% c(renewal$cases_earlier, cases)),
                 statistics$infectiousness)
    expect_equal(
      dated_renewal_logp(cases, renewal),
      state_process_nb_logp_all_onsets(
        4, dplyr::tibble(time_onset=1:6, cases_true=cases, Rt=1.3), serial, 5,
        cases_before))
    renewal$is_negative_binomial <- FALSE
    is_informative <- statistics$is_informative
    expect_equal(
      dated_renewal_logp(cases, renewal),
      sum(stats::dpois(cases[is_informative],
                       1.3 * statistics$infectiousness[is_informative], log=TRUE)))
  }
})

test_that("changes to the fatal cases dated by death keep the recorded counts and give the right change in log-probability", {
  set.seed(13)
  test <- dated_test_state(reporting=list(mean=2, sd=1), cases_before=rep(40, 10))
  state <- test$state
  cases_true <- test$cases_true
  f_full <- function(state, cases_true) dated_logp(state, cases_true)
  n_applied <- c(relabel=0, redate=0, unseen=0, inflow=0)
  for(i in 1:400) {
    type <- sample(names(n_applied), 1)
    move <- list(type=type, day=sample.int(25, 1),
                 lag=sample.int(state$dated$max_lag, 1),
                 lag_to=sample.int(state$dated$max_lag, 1), sign=sample(c(-1, 1), 1))
    proposed <- dated_apply_move(state, cases_true, move)
    if(is.null(proposed))
      next
    logp_proposed <- f_full(proposed$state, proposed$cases_true)
    if(!is.finite(logp_proposed)) {
      # impossible states are also found from the parts which change
      expect_equal(dated_logp(proposed$state, proposed$cases_true,
                              proposed$cohorts, proposed$days), -Inf)
      next
    }
    expect_equal(
      dated_logp(proposed$state, proposed$cases_true, proposed$cohorts, proposed$days) -
        dated_logp(state, cases_true, proposed$cohorts, proposed$days),
      logp_proposed - f_full(state, cases_true))
    state <- proposed$state
    cases_true <- proposed$cases_true
    n_applied[type] <- n_applied[type] + 1
  }
  expect_true(all(n_applied > 5))
  expect_dated_consistent(state, cases_true, test$recorded)
  # the fatal cases recorded for each day are not altered either
  expect_equal(dated_counts(state)$fatal_recorded,
               dated_counts(test$state)$fatal_recorded)
  expect_false(identical(cases_true, test$cases_true))
  expect_error(dated_apply_move(state, cases_true, list(type="other", day=1, sign=1)))
  # a case cannot have had its onset before the first day unless earlier
  # cases are known, nor die after the last day if there is time for all to die
  expect_null(dated_apply_move(state, cases_true,
                               list(type="unseen", day=1, sign=1)))
  expect_null(dated_apply_move(state, cases_true,
                               list(type="inflow", day=25, sign=1)))
  expect_null(dated_apply_move(state, cases_true,
                               list(type="relabel", day=2, lag=3, sign=1)))
})

test_that("deaths_dated_step moves cases between days of onset and keeps the recorded counts", {
  set.seed(14)
  test <- dated_test_state()
  stepped <- deaths_dated_step(test$state, test$cases_true, test$renewal)
  expect_dated_consistent(stepped$state, stepped$cases_true, test$recorded)
  expect_false(identical(stepped$cases_true, test$cases_true))
  expect_equal(sum(stepped$state$dated$inflow), 0)
  expect_true(is.finite(dated_logp(stepped$state, stepped$cases_true)))
  # the other steps still run, and leave the recorded counts possible
  state <- deaths_paths_step(stepped$state, stepped$cases_true,
                             dplyr::tibble(mean=rep(5, 25), sd=rep(3, 25)))
  f_logp <- deaths_cases_logp_function(state, 20)
  min_cases <- attr(f_logp, "min_cases")
  expect_true(is.finite(f_logp(stepped$cases_true[20])))
  expect_equal(f_logp(min_cases - 1), -Inf)
  expect_true(is.finite(f_logp(min_cases)))
  state <- deaths_true_step(state, stepped$cases_true)
  expect_dated_consistent(state, stepped$cases_true, test$recorded)
  state <- deaths_parameters_step(state, stepped$cases_true, test$settings)
  expect_true(all(state$cfr > 0 & state$cfr < 1))
  draws <- deaths_draws(state, 1)
  expect_equal(draws$deaths$dated_by_death, state$dated$into)
  expect_equal(draws$deaths$deaths_recorded, dated_counts(state)$fatal_recorded)
  # without such cases the step does nothing
  state$dated <- NULL
  expect_equal(deaths_dated_step(state, stepped$cases_true, test$renewal)$cases_true,
               stepped$cases_true)
})

test_that("mcmc_single runs with fatal cases dated by death", {
  set.seed(15)
  n_days <- 40
  serial_parameters <- list(mean=5, sd=3)
  cases <- true_cases(n_days + 40, function(t) 1, kappa=1000, serial_parameters,
                      initial_parameters=list(mean=30, length=10))
  cases_before <- cases[1:40]
  cases <- cases[41:(n_days + 40)]
  dated_by_death <- list(fraction=0.4, mean=6, sd=3, reporting=list(mean=2, sd=1))
  sim <- simulate_cases_deaths(cases, seq_len(n_days + 30), list(mean=4, sd=3), 0.5,
                               list(mean=8, sd=4), dated_by_death=dated_by_death)
  snapshots <- sim$snapshots
  initial_cases <- snapshots %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(cases_true=max(cases_reported))
  priors <- list(Rt=list(shape=1, rate=1),
                 reporting=list(mean_mu=5, mean_sigma=5, sd_mu=3, sd_sigma=3),
                 max_cases=300,
                 Rt_changepoints=list(),
                 reporting_changepoints=list(min_length=14))
  fit <- mcmc_single(
    niterations=20,
    snapshot_with_Rt_index_df=snapshots %>%
      dplyr::select(time_onset, time_reported, cases_reported) %>%
      dplyr::mutate(Rt_index=1),
    priors=priors,
    serial_parameters=serial_parameters,
    initial_cases_true=initial_cases,
    initial_reporting_parameters=list(mean=4, sd=3),
    initial_Rt=dplyr::tibble(Rt_index=1, Rt=1),
    print_to_screen=FALSE,
    cases_before=cases_before,
    deaths=list(snapshot=snapshots %>%
                  dplyr::select(time_onset, time_reported, deaths_reported),
                settings=list(dated_by_death=dated_by_death,
                              initial_delay=list(mean=8, sd=4))))
  expect_equal(nrow(fit$deaths), 20 * n_days)
  expect_true(all(c("deaths_recorded", "dated_by_death") %in% colnames(fit$deaths)))
  both <- dplyr::left_join(fit$deaths, fit$cases, by=c("iteration", "time_onset"))
  expect_true(all(both$deaths_true <= both$cases_true))
  expect_true(all(both$dated_by_death <= both$deaths_recorded))
  # the fatal cases recorded for early days are all shown as dead by the end
  latest <- snapshots %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(deaths=max(deaths_reported), cases=max(cases_reported))
  last <- dplyr::filter(fit$deaths, iteration == 20)
  expect_equal(last$deaths_recorded[1:15], latest$deaths[1:15])
  # true cases by day of onset differ from the cases recorded for each day
  last_cases <- dplyr::filter(fit$cases, iteration == 20)
  expect_false(identical(last_cases$cases_true[1:15], latest$cases[1:15]))
  # but cases are only moved between days, so that the totals are close
  expect_true(abs(sum(last_cases$cases_true) / sum(latest$cases) - 1) < 0.1)
})

# a state in which the fraction is estimated and some fatal cases recorded
# for the first days had their onset before the first day
dated_collapsed_state <- function() {
  set.seed(21)
  test <- dated_test_state(cases_before=rep(30, 20))
  state <- test$state
  state$dated$fraction_prior <- list(alpha=2, beta=3)
  state$dated$collapse_fraction <- TRUE
  c(test, list(collapsed=state))
}

test_that("the factors left out with the fraction integrated out are those of dated_fraction_terms", {
  test <- dated_collapsed_state()
  state <- test$collapsed
  expect_true(sum(state$dated$inflow_weight) > 0)
  terms <- dated_fraction_terms(state)
  prior <- state$dated$fraction_prior
  collapsed <- dated_logp(state, test$cases_true, collapsed=TRUE)
  for(f in c(0.2, 0.55)) {
    state$dated$fraction[] <- f
    full <- dated_logp(state, test$cases_true)
    expect_equal(full - collapsed,
                 (terms$A - prior$alpha) * log(f) + (terms$B - prior$beta) * log(1 - f) -
                   terms$r * f)
  }
})

test_that("dated_fraction_marginal is the integral over the fraction, and dated_fraction_draw draws from it", {
  test <- dated_collapsed_state()
  state <- test$collapsed
  prior <- state$dated$fraction_prior
  for(r_scale in c(0, 1, 40)) {
    state$dated$inflow_weight <- test$collapsed$dated$inflow_weight * r_scale
    if(r_scale == 0)
      state$dated$inflow[] <- 0
    terms <- dated_fraction_terms(state)
    f_log_density <- function(f) (terms$A - 1) * log(f) + (terms$B - 1) * log(1 - f) -
      terms$r * f - lbeta(prior$alpha, prior$beta)
    # integrated about the mode, scaled to avoid underflow
    mode <- stats::optimize(f_log_density, c(0, 1), maximum=TRUE)$maximum
    peak <- f_log_density(mode)
    integral <- stats::integrate(function(f) exp(f_log_density(f) - peak), 0, 1,
                                 rel.tol=1e-10)$value
    marginal <- dated_fraction_marginal(state)
    expect_equal(marginal$log_marginal, peak + log(integral), tolerance=1e-8)
    mean_f <- stats::integrate(function(f) f * exp(f_log_density(f) - peak), 0, 1,
                               rel.tol=1e-10)$value / integral
    sd_f <- sqrt(stats::integrate(function(f) (f - mean_f)^2 * exp(f_log_density(f) - peak),
                                  0, 1, rel.tol=1e-10)$value / integral)
    set.seed(22)
    draws <- replicate(4000, dated_fraction_draw(state))
    expect_lt(abs(mean(draws) - mean_f), 4 * sd_f / sqrt(4000))
    expect_equal(stats::sd(draws), sd_f, tolerance=0.05)
  }
})

test_that("with the fraction integrated out the changes give the right change in log-probability", {
  test <- dated_collapsed_state()
  state <- test$collapsed
  cases_true <- test$cases_true
  f_full <- function(s, c) dated_logp(s, c, collapsed=TRUE) +
    dated_fraction_marginal(s)$log_marginal
  set.seed(23)
  n_checked <- 0
  for(i in 1:200) {
    type <- sample(c("relabel", "redate", "unseen", "inflow"), 1)
    move <- list(type=type, day=sample.int(25, 1), lag=sample.int(5, 1),
                 lag_to=sample.int(5, 1), sign=sample(c(-1, 1), 1))
    proposed <- dated_apply_move(state, cases_true, move)
    if(is.null(proposed) || !is.finite(f_full(proposed$state, proposed$cases_true)))
      next
    n_checked <- n_checked + 1
    expect_equal(
      dated_logp(proposed$state, proposed$cases_true, proposed$cohorts, proposed$days,
                 TRUE) + dated_fraction_marginal(proposed$state)$log_marginal -
        dated_logp(state, cases_true, proposed$cohorts, proposed$days, TRUE) -
        dated_fraction_marginal(state)$log_marginal,
      f_full(proposed$state, proposed$cases_true) - f_full(state, cases_true))
  }
  expect_gt(n_checked, 20)
  # the step runs, keeps the recorded counts, and the fraction is then drawn
  stepped <- deaths_dated_step(state, cases_true, test$renewal)
  expect_dated_consistent(stepped$state, stepped$cases_true, test$recorded)
  drawn <- deaths_dated_fraction_step(stepped$state)
  expect_false(identical(drawn$dated$fraction, stepped$state$dated$fraction))
  expect_true(all(drawn$dated$fraction == drawn$dated$fraction[1]))
})

test_that("the fraction can be integrated out only if it is estimated", {
  expect_error(dated_by_death_settings(10, list(fraction=0.3, mean=5, sd=3,
                                                collapse_fraction=TRUE)),
               "only if it is estimated")
  settings <- dated_by_death_settings(10, list(fraction=list(alpha=1, beta=1), mean=5,
                                               sd=3, collapse_fraction=TRUE))
  expect_true(settings$collapse_fraction)
  expect_false(dated_by_death_settings(10, list(fraction=list(alpha=1, beta=1), mean=5,
                                                sd=3))$collapse_fraction)
})
