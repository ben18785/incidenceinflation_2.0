f_relaxed_data <- function(seed=1, days_total=14) {
  set.seed(seed)
  simulated <- simulate_report_multipliers(
    1:days_total, c("1"=0.2, "2"=-0.3, "3"=0.1, "4"=0, "5"=0.15, "6"=-0.1, "7"=-0.05), 0)
  revision_parameters <- list(eta=0.08, rho=0.25,
                              report_multipliers=simulated$report_multipliers)
  df <- generate_snapshots(days_total, function(t) rep(1.2, length(t)),
                           list(mean=3, sd=2), list(mean=4, sd=2.5),
                           initial_parameters=list(mean=40, length=5),
                           revision_parameters=revision_parameters)
  list(df=df, revision_parameters=revision_parameters)
}

test_that("relaxed_stan_data lays out the data for the model", {
  l <- f_relaxed_data()
  df <- l$df %>% dplyr::filter(!(time_reported %in% c(3, 6)))
  data <- relaxed_stan_data(df, list(mean=3, sd=2), cases_before=c(5, 7, 9),
                            weekday_first=3, days_sampled=6, serial_max=10,
                            Rt_period_width=7)
  expect_equal(data$n_days, 14)
  expect_equal(data$n_obs, nrow(df))
  expect_equal(sum(data$obs_n), nrow(df))
  # rows of each onset day are contiguous and in order of report
  for(t in c(1, 5, 14)) {
    rows <- data$obs_start[t] + seq_len(data$obs_n[t]) - 1
    expected <- df %>% dplyr::filter(time_onset == t) %>% dplyr::arrange(time_reported)
    expect_equal(data$obs_day[rows], as.integer(expected$time_reported))
    expect_equal(data$obs_count[rows], as.integer(expected$cases_reported))
  }
  # the latest onset times are those whose true cases are unknown
  expect_equal(data$n_sampled, 6)
  expect_equal(as.numeric(data$sampled_day), 9:14)
  expect_true(all(data$cases_lower >= data$cases_fixed))
  # earlier cases are padded with zeros to the length of the serial interval
  expect_equal(data$cases_before, c(rep(0, 7), 5, 7, 9))
  expect_equal(sum(data$serial_weights), 1)
  # periods end on the last onset time
  expect_equal(data$period, rep(1:2, each=7))
  expect_equal(relaxed_stan_data(df, list(mean=3, sd=2), Rt_period_width=5)$period,
               c(rep(1, 4), rep(2, 5), rep(3, 5)))
  # each day takes the day of the week of the first report on or after it
  expect_equal(data$weekday_report[2], as.integer(weekday_of_days(2, 3)))
  expect_equal(data$weekday_report[3], as.integer(weekday_of_days(4, 3)))
  expect_equal(data$weekday_report[6], as.integer(weekday_of_days(7, 3)))
  expect_equal(data$weekday_report[14], as.integer(weekday_of_days(14, 3)))
  expect_equal(data$is_weekday, 1)
  no_weekday <- relaxed_stan_data(df, list(mean=3, sd=2))
  expect_equal(no_weekday$is_weekday, 0)
  expect_true(all(no_weekday$weekday_report == 0))
  expect_true(all(names(relaxed_priors()) %in% names(data)))

  expect_error(relaxed_stan_data(df %>% dplyr::filter(time_onset != 4),
                                 list(mean=3, sd=2)))
  expect_error(relaxed_stan_data(df %>% dplyr::mutate(time_reported=time_reported + 5),
                                 list(mean=3, sd=2)))
})

test_that("relaxed_inits gives valid starting values", {
  l <- f_relaxed_data()
  data <- relaxed_stan_data(l$df, list(mean=3, sd=2), days_sampled=5)
  init <- relaxed_inits(data)()
  expect_equal(length(init$cases_extra), 5)
  expect_true(all(init$cases_extra > 0))
  expect_equal(sum(init$log_weekday), 0)
  expect_equal(length(init$Rt_step_raw), data$n_periods - 1)
})

f_skip_without_stan <- function() {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_version(), error=function(e) NULL)),
          "CmdStan is not installed")
}

test_that("the likelihood of the relaxed model matches that of the package", {
  f_skip_without_stan()
  l <- f_relaxed_data(2)
  df <- l$df %>% dplyr::filter(!(time_reported %in% c(3, 6, 7, 11)))
  model <- relaxed_model(compile_standalone=TRUE, force_recompile=TRUE)
  model$expose_functions(global=FALSE)
  f <- model$functions
  multipliers <- l$revision_parameters$report_multipliers
  daily <- rep(1, 14)
  daily[multipliers$day_first - 1 + seq_along(multipliers$values)] <- multipliers$values

  # the survival function, found by numerical integration
  for(parameters in list(c(7, 4), c(3, 5), c(15, 12))) {
    shape <- (parameters[1] / parameters[2])^2
    rate <- parameters[1] / parameters[2]^2
    expect_equal(f$gamma_log_survival(60L, shape, rate),
                 pgamma(1:60 - 0.5, shape, rate, lower.tail=FALSE, log.p=TRUE),
                 tolerance=1e-3)
  }
  log_undetected <- f$log_undetected_matrix(14L, 4, 2.5, daily)
  expect_equal(log_undetected[5, 1:10],
               log_undetected_days(5:14, 5, list(mean=4, sd=2.5),
                                   list(report_multipliers=multipliers)),
               tolerance=1e-4)

  # reported counts of each onset day, for whole numbers of true cases
  revision_parameters <- list(eta=0.05, rho=0.3, report_multipliers=multipliers)
  for(t in 1:14) {
    observation_df <- df %>%
      dplyr::filter(time_onset == t) %>%
      dplyr::arrange(time_reported)
    if(nrow(observation_df) <= 1)
      next
    for(extra in c(0, 4)) {
      cases <- observation_df$cases_true[1] + extra
      expect_equal(
        f$observation_logp(cases, as.integer(observation_df$time_reported),
                           as.integer(observation_df$cases_reported), t,
                           log_undetected[t, ], 0.05, 0.3, 3L),
        observation_process_revisions_logp(observation_df, cases, t,
                                           list(mean=4, sd=2.5),
                                           revision_parameters, 3),
        tolerance=1e-3)
    }
  }
  expect_equal(f$neg_binomial_real_logp(7, 5.5, 3), dnbinom(7, mu=5.5, size=3, log=TRUE))
  expect_equal(f$binomial_real_logp(3, 10, 0.3), dbinom(3, 10, 0.3, log=TRUE))
  expect_equal(f$binomial_real_logp(11, 10, 0.3), -Inf)
})

test_that("the relaxed model can be fitted and its draws extracted", {
  f_skip_without_stan()
  set.seed(3)
  df <- generate_snapshots(30, function(t) rep(1.1, length(t)), list(mean=4, sd=2),
                           list(mean=4, sd=2.5),
                           initial_parameters=list(mean=40, length=10),
                           revision_parameters=list(eta=0.03, rho=0.3))
  data <- relaxed_stan_data(dplyr::select(df, -cases_true), list(mean=4, sd=2),
                            days_sampled=12, weekday_first=1)
  model <- relaxed_model()
  fit <- suppressMessages(fit_relaxed(data, "laplace", model, seed=4, refresh=0,
                                      draws=200, show_messages=FALSE))
  draws <- relaxed_draws(fit)
  expect_equal(colnames(draws$cases), c("iteration", "time_onset", "cases_true"))
  expect_equal(nrow(draws$cases), 200 * 30)
  expect_equal(nrow(draws$Rt), 200 * 30)
  expect_true(all(c("delay_mean", "delay_sd", "eta", "rho", "kappa") %in%
                    colnames(draws$parameters)))
  # true cases of old onset times are fixed, those of recent ones exceed reports
  truth <- unique(df[, c("time_onset", "cases_true")])
  medians <- draws$cases %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(cases=median(cases_true), spread=sd(cases_true))
  expect_true(all(medians$spread[1:18] == 0))
  expect_true(all(medians$spread[19:30] > 0))
  expect_true(all(medians$cases[19:30] >= data$cases_lower[19:30]))
  # estimates of recent cases are of the right size
  expect_equal(sum(medians$cases[19:30]), sum(truth$cases_true[19:30]), tolerance=0.3)
  expect_equal(median(draws$parameters$delay_mean), 4, tolerance=0.3)
  expect_error(fit_relaxed(data, "nonsense", model))
})
