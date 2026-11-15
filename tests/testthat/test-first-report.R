# The first report of an onset day is part of the probability of its reports:
# every count is zero at a notional report before it, so the cases reported
# by the first report are binomial given the true cases.

delay <- list(mean = 5, sd = 2)
f_survival <- function(day, day_onset = 1, reporting_parameters = delay)
  exp(log_undetected_prob(day + 0.5, day_onset, reporting_parameters))

test_that("a single report of no cases thins the renewal prior of the true cases", {
  observation_df <- dplyr::tibble(time_reported = 7, cases_reported = 0)
  cases <- 0:200
  f_mean <- function(df) {
    logp <- conditional_cases_logp(cases, df, rep(20, 40), Rt = 1, day_onset = 1,
                                   serial_parameters = list(mean = 6, sd = 3),
                                   reporting_parameters = delay)
    p <- exp(logp - matrixStats::logSumExp(logp))
    sum(cases * p)
  }
  prior_mean <- f_mean(observation_df[0, ])
  expect_equal(prior_mean, 20, tolerance = 1e-6)
  # Poisson(mu) thinned by a report of zero is Poisson(mu * S)
  expect_equal(f_mean(observation_df), prior_mean * f_survival(7), tolerance = 1e-6)
})

test_that("without wrong records the reports of an onset day telescope", {
  observation_df <- dplyr::tibble(time_reported = c(3, 6, 10), cases_reported = c(2, 7, 9))
  cases <- 9:40
  logp <- observation_process_logp(observation_df, cases, day_onset = 1,
                                   reporting_parameters = delay)
  # as a function of the true cases, only the last report matters
  logp_last <- stats::dbinom(9, cases, 1 - f_survival(10), log = TRUE)
  expect_equal(logp - logp_last, rep(logp[1] - logp_last[1], length(cases)), tolerance = 1e-8)
  # and the whole probability is that of the chain of binomials from zero
  days <- c(1, observation_df$time_reported)
  counts <- c(0, observation_df$cases_reported)
  survival <- c(1, f_survival(observation_df$time_reported))
  exact <- sapply(cases, function(n) sum(stats::dbinom(
    diff(counts), n - counts[-4], 1 - survival[-1] / survival[-4], log = TRUE)))
  expect_equal(logp, exact, tolerance = 1e-8)
  # a single report gives the binomial of its count
  expect_equal(observation_process_logp(observation_df[1, ], cases, 1, delay),
               stats::dbinom(2, cases, 1 - f_survival(3), log = TRUE), tolerance = 1e-8)
})

test_that("the forward filter of wrong records includes the first report", {
  revision_parameters <- list(eta = 0.3, rho = 0.2)
  observation_df <- dplyr::tibble(time_reported = c(4, 8), cases_reported = c(3, 4))
  cases <- 2:12
  # by brute force, from the definitions: wrong records present at the first
  # report are those entered on each day after onset and not yet withdrawn
  s <- f_survival(1:8)
  mean_first <- sum((s[1:3] - s[2:4]) * (1 - 0.2)^(4 - 2:4))
  p_withdrawn <- 1 - (1 - 0.2)^4
  f_exact <- function(n) {
    total <- 0
    for(e1 in 0:3) {
      d1 <- 3 - e1
      p1 <- stats::dpois(e1, 0.3 * n * mean_first) * stats::dbinom(d1, n, 1 - s[4])
      if(p1 == 0) next
      for(m in 0:e1) for(entered in 0:20) {
        d <- 4 - 3 - entered + m
        if(d < 0) next
        total <- total + p1 * stats::dbinom(m, e1, p_withdrawn) *
          stats::dpois(entered, 0.3 * n * (s[4] - s[8])) *
          stats::dbinom(d, n - d1, 1 - s[8] / s[4])
      }
    }
    log(total)
  }
  exact <- sapply(cases, f_exact)
  for(use_compiled in c(TRUE, FALSE))
    expect_equal(observation_process_revisions_logp(
      observation_df, cases, 1, delay, revision_parameters, max_erroneous = 15,
      use_compiled = use_compiled), exact, tolerance = 1e-6)
  # a single report: genuine cases and wrong records present at it
  single <- sapply(cases, function(n) log(sum(
    stats::dpois(0:3, 0.3 * n * mean_first) * stats::dbinom(3 - 0:3, n, 1 - s[4]))))
  for(use_compiled in c(TRUE, FALSE))
    expect_equal(observation_process_revisions_logp(
      observation_df[1, ], cases, 1, delay, revision_parameters, max_erroneous = 15,
      use_compiled = use_compiled), single, tolerance = 1e-8)
  # a first report on the onset day has no wrong records yet
  expect_equal(observation_process_revisions_logp(
    dplyr::tibble(time_reported = 1, cases_reported = 1), cases, 1, delay,
    revision_parameters),
    stats::dbinom(1, cases, 1 - f_survival(1), log = TRUE), tolerance = 1e-8)
})

test_that("an onset day with a single report informs the reporting delay", {
  df <- dplyr::tibble(time_onset = c(1, 2), time_reported = c(8, 8),
                      cases_reported = c(6, 1), cases_true = c(10, 10),
                      reporting_piece_index = 1)
  f_logp <- function(mean) observation_process_all_times_logp(
    df, dplyr::tibble(reporting_piece_index = 1, mean = mean, sd = 2))
  expected <- function(mean) {
    rp <- list(mean = mean, sd = 2)
    stats::dbinom(6, 10, 1 - f_survival(8, 1, rp), log = TRUE) +
      stats::dbinom(1, 10, 1 - f_survival(8, 2, rp), log = TRUE)
  }
  expect_equal(f_logp(4), expected(4), tolerance = 1e-8)
  expect_equal(f_logp(9), expected(9), tolerance = 1e-8)
  expect_false(isTRUE(all.equal(f_logp(4), f_logp(9))))
  # the same through the function used with changes in the delay
  f_changepoints <- reporting_logp_function(df)
  expect_equal(f_changepoints(1:2, list(mean = 4, sd = 2)), expected(4), tolerance = 1e-8)
})

test_that("the joint model of cases and deaths includes the first report", {
  # one report: the survivors and fatal cases newly reported are binomial
  # from zero
  probabilities <- list(p = 0.3, g = 0.5, q = 0.2)
  set.seed(1)
  path <- sample_fatal_path(4, 1, 3, 3, probabilities)
  survivals <- list(log_reported = log(1 - probabilities$p),
                    log_recorded = log(1 - probabilities$g), is_dated_same = TRUE)
  survivals$log_dated <- survivals$log_reported
  logp <- joint_cases_deaths_logp(path, survivals, 0.4, 4:9)
  fatal_values <- attr(logp, "deaths_true")
  f_exact <- function(cases_true, deaths_true)
    stats::dbinom(deaths_true, cases_true, 0.4, log = TRUE) +
    stats::dbinom(path$fatal_new, deaths_true, 0.3, log = TRUE) +
    stats::dbinom(path$survivors_new, cases_true - deaths_true, 0.3, log = TRUE)
  for(i in 1:3)
    expect_equal(logp[i + 1, 1], f_exact(4 + i, fatal_values[1]), tolerance = 1e-8)
  # the probability of being reported by the first report is from zero
  expect_equal(survival_steps(log(c(0.6, 0.2))), c(0.4, 1 - 0.2 / 0.6))
})

test_that("the simulators report from a notional report of zero before the first", {
  set.seed(2)
  n <- 2000
  # a report on the onset day holds the cases reported within half a day
  delay_fast <- list(mean = 2, sd = 2)
  counts <- replicate(400, observed_cases_trajectory(n, c(1, 2), 1, delay_fast))
  expect_equal(mean(counts[1, ]), n * (1 - f_survival(1, 1, delay_fast)), tolerance = 0.02)
  expect_equal(mean(counts[2, ]), n * (1 - f_survival(2, 1, delay_fast)), tolerance = 0.02)
  # a first report days after onset, with wrong records entered before it
  revision_parameters <- list(eta = 0.2, rho = 0.3)
  counts <- replicate(400, observed_cases_trajectory(n, 4, 1, delay,
                                                     revision_parameters))
  s <- f_survival(1:4)
  expected <- n * (1 - s[4]) +
    0.2 * n * sum((s[1:3] - s[2:4]) * (1 - 0.3)^(4 - 2:4))
  expect_equal(mean(counts), expected, tolerance = 0.02)
  # the simulator of whole outbreaks, on the onset day and the day after
  sim <- purrr::map_dfr(1:200, function(i) observed_cases_revisions(
    rep(n, 3), delay_fast, list(eta = 0.2, rho = 0.3)))
  by_lag <- sim %>%
    dplyr::filter(time_onset == 1) %>%
    dplyr::group_by(time_reported) %>%
    dplyr::summarise(mean = mean(cases_reported))
  s <- f_survival(1:3, 1, delay_fast)
  expect_equal(by_lag$mean[1], n * (1 - s[1]), tolerance = 0.02)
  expect_equal(by_lag$mean[2], n * (1 - s[2]) + 0.2 * n * (s[1] - s[2]), tolerance = 0.02)
})
