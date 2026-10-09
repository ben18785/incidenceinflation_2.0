f_skip_without_stan <- function() {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_version(), error=function(e) NULL)),
          "CmdStan is not installed")
}

# simulated reports of cases and deaths whose days start at one and end on the
# last onset day, as the relaxed model needs
relaxed_deaths_test_data <- function(cfr, n_days=40, wrong=NULL) {
  serial <- list(mean=5, sd=3)
  cases_all <- true_cases(n_days + 40, function(t) 1, kappa=1000, serial,
                          initial_parameters=list(mean=40, length=10))
  sim <- simulate_cases_deaths(cases_all[41:(n_days + 40)], seq_len(n_days),
                               list(mean=5, sd=3), cfr, list(mean=8, sd=4),
                               wrong=wrong)
  c(sim, list(serial=serial, cases_before=cases_all[1:40]))
}

test_that("relaxed_stan_data lays out the deaths for the model", {
  set.seed(31)
  sim <- relaxed_deaths_test_data(0.4, n_days=20)
  cases_df <- dplyr::select(sim$snapshots, time_onset, time_reported, cases_reported)
  deaths_df <- dplyr::select(sim$snapshots, time_onset, time_reported, deaths_reported)
  without <- relaxed_stan_data(cases_df, sim$serial)
  expect_equal(without$is_deaths, 0)
  expect_equal(length(without$obs_deaths), 0)
  expect_equal(without$n_cfr_segments, 1)
  expect_equal(without$max_dead, rep(0, 20))
  data <- relaxed_stan_data(cases_df, sim$serial, weekday_first=2, deaths_df=deaths_df,
                            cfr_max_changes=2, recording_max_changes=1)
  expect_equal(data$is_deaths, 1)
  expect_equal(as.numeric(data$obs_deaths),
               dplyr::arrange(sim$snapshots, time_onset, time_reported)$deaths_reported)
  expect_equal(data$n_cfr_segments, 3)
  expect_equal(data$n_recording_segments, 2)
  expect_equal(data$is_recording_weekday, 1)
  expect_equal(relaxed_stan_data(cases_df, sim$serial, weekday_first=2,
                                 deaths_df=deaths_df,
                                 recording_weekday=FALSE)$is_recording_weekday, 0)
  # no count of deaths falls, so no wrong records shown as dead are needed
  expect_equal(data$max_dead, rep(0, 20))
  fallen <- deaths_df
  row <- which(fallen$time_onset == 3 & fallen$time_reported == 20)
  fallen$deaths_reported[row] <- fallen$deaths_reported[row] - 1
  expect_equal(relaxed_stan_data(cases_df, sim$serial, deaths_df=fallen,
                                 max_erroneous=2)$max_dead[1:4], c(0, 0, 2, 0))
  inits <- relaxed_inits(data)()
  expect_equal(length(inits$cfr), 3)
  expect_equal(length(inits$recording_mean), 2)
  expect_equal(length(inits$cfr_p_change), 1)
  expect_equal(length(relaxed_inits(without)()$cfr), 0)
  expect_error(relaxed_stan_data(cases_df, sim$serial, deaths_df=deaths_df[-1, ]))
  expect_error(relaxed_stan_data(cases_df, sim$serial, deaths_df=cases_df))
  expect_error(relaxed_stan_data(cases_df, sim$serial, deaths_df=deaths_df,
                                 cfr_max_changes=-1))
  too_many <- deaths_df
  too_many$deaths_reported[5] <- 1000
  expect_error(relaxed_stan_data(cases_df, sim$serial, deaths_df=too_many))
})

test_that("the likelihood of deaths in the relaxed model is exact without wrong records", {
  f_skip_without_stan()
  model <- relaxed_model(compile_standalone=TRUE, force_recompile=TRUE)
  model$expose_functions(global=FALSE)
  f <- model$functions
  days <- c(3L, 4L, 6L, 7L, 10L)
  onset <- 2L
  cases <- c(1L, 3L, 4L, 6L, 6L)
  deaths <- c(0L, 1L, 1L, 2L, 3L)
  cases_true <- 8
  log_unrecorded <- log_undetected_days(onset:10, onset, list(mean=6, sd=4))
  probabilities <- joint_probabilities(days, onset, list(mean=4, sd=3),
                                       list(mean=6, sd=4))
  logp_cases <- sum(stats::dbinom(diff(c(0, cases)), cases_true - c(0, cases[-5]),
                                  probabilities$p, log=TRUE))
  for(cfr in c(0.2, 0.6)) {
    # the fatal cases and the split of each report summed out
    logp_joint <- log(sum(sapply(0:cases_true, function(fatal)
      stats::dbinom(fatal, cases_true, cfr) *
        exp(joint_path_forward(cases, deaths, fatal, cases_true - fatal,
                               probabilities)$logp))))
    expect_equal(
      f$deaths_logp(days, cases, deaths, onset, cfr, log_unrecorded, rep(0, 5), 0.3, 0L),
      logp_joint - logp_cases)
    # allowing for wrong records which cannot be entered changes nothing
    expect_equal(
      f$deaths_logp(days, cases, deaths, onset, cfr, log_unrecorded, rep(0, 5), 0.3, 2L),
      logp_joint - logp_cases)
  }
  # a fall in the deaths is possible only with a wrong record shown as dead
  fallen <- c(0L, 2L, 1L, 2L, 3L)
  expect_equal(f$deaths_logp(days, cases, fallen, onset, 0.4, log_unrecorded,
                             rep(0.2, 5), 0.3, 0L), -Inf)
  expect_true(is.finite(f$deaths_logp(days, cases, fallen, onset, 0.4, log_unrecorded,
                                      rep(0.2, 5), 0.3, 2L)))

  # two quantities which change independently: every pair of paths summed
  set.seed(32)
  ll <- matrix(stats::rnorm(4 * 4), 4, 4)
  paths <- expand.grid(a2=1:2, a3=1:2, a4=1:2, b2=1:2, b3=1:2, b4=1:2)
  f_weight <- function(path, p) {
    path <- c(1, path)
    if(any(diff(path) < 0))
      return(0)
    p^sum(diff(path) == 1) * (1 - p)^sum(diff(path) == 0)
  }
  total <- sum(apply(paths, 1, function(row) {
    a <- c(1, row[1:3])
    b <- c(1, row[4:6])
    f_weight(row[1:3], 0.2) * f_weight(row[4:6], 0.3) *
      exp(sum(ll[cbind(1:4, (a - 1) * 2 + b)]))
  }))
  expect_equal(f$changepoint2_marginal(ll, 2L, 2L, 0.2, 0.3), log(total))
  # with a single segment of each there is nothing to sum over
  expect_equal(f$changepoint2_marginal(ll[, 1, drop=FALSE], 1L, 1L, 0.5, 0.5),
               sum(ll[, 1]))
  # the probabilities that each day belongs to each segment sum to one
  smoothed <- f$changepoint_smoothed(ll[, 1:2], 0.2)
  expect_equal(rowSums(smoothed), rep(1, 4))
  expect_equal(smoothed[1, ], c(1, 0))
})

test_that("the relaxed model with deaths recovers a case fatality ratio which changes", {
  f_skip_without_stan()
  set.seed(33)
  n_days <- 40
  cfr <- c(rep(0.6, 20), rep(0.2, 20))
  sim <- relaxed_deaths_test_data(cfr, n_days)
  data <- relaxed_stan_data(
    dplyr::select(sim$snapshots, time_onset, time_reported, cases_reported),
    sim$serial, cases_before=sim$cases_before, days_sampled=20,
    deaths_df=dplyr::select(sim$snapshots, time_onset, time_reported, deaths_reported),
    cfr_max_changes=1,
    priors=utils::modifyList(relaxed_priors(),
                             list(recording_mean_mu=8, recording_sd_mu=4)))
  fit <- suppressMessages(fit_relaxed(data, "laplace", relaxed_model(), seed=5,
                                      refresh=0, draws=200, show_messages=FALSE))
  draws <- relaxed_draws(fit)
  expect_equal(colnames(draws$deaths),
               c("iteration", "time_onset", "deaths_true", "cfr", "delay_mean", "delay_sd"))
  expect_equal(nrow(draws$deaths), 200 * n_days)
  estimates <- draws$deaths %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(cfr=median(cfr), deaths=median(deaths_true), delay=median(delay_mean))
  expect_true(abs(estimates$cfr[5] - 0.6) < 0.08)
  expect_true(abs(estimates$cfr[35] - 0.2) < 0.08)
  expect_true(mean(draws$parameters$cfr_n_changes == 1) > 0.9)
  expect_true(abs(estimates$delay[20] - 8) < 2)
  # the deaths of early onset times have all been shown
  expect_equal(estimates$deaths[1:10], sim$truth$deaths_true[1:10])
  both <- dplyr::left_join(draws$deaths, draws$cases, by=c("iteration", "time_onset"))
  expect_true(all(both$deaths_true <= round(both$cases_true)))
  # without deaths there are no such draws
  data_cases <- relaxed_stan_data(
    dplyr::select(sim$snapshots, time_onset, time_reported, cases_reported),
    sim$serial, cases_before=sim$cases_before, days_sampled=20)
  fit_cases <- suppressMessages(fit_relaxed(data_cases, "laplace", relaxed_model(),
                                            seed=5, refresh=0, draws=50,
                                            show_messages=FALSE))
  expect_null(relaxed_draws(fit_cases)$deaths)
})
