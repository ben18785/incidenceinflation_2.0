f_log_base <- function(days, day_onset, reporting_parameters) {
  purrr::map_dbl(days, ~log(undetected_prob(. + 0.5, day_onset, reporting_parameters)))
}

test_that("log_undetected_days matches undetected_prob without extensions", {
  reporting_parameters <- list(mean=5, sd=3)
  days <- c(4, 5, 8, 20)
  expected <- f_log_base(days, 4, reporting_parameters)
  expect_equal(log_undetected_days(days, 4, reporting_parameters), expected)
  expect_equal(log_undetected_days(days, 4, reporting_parameters, list(eta=0.1, rho=0.2)),
               expected)
  # a tail of zero weight makes no difference
  expect_equal(log_undetected_days(days, 4, reporting_parameters,
                                   list(tail=list(weight=0, mean=10))),
               expected)
  # remains defined long after onset
  expect_true(is.finite(log_undetected_days(500, 4, reporting_parameters)))
})

test_that("log_undetected_days mixes in an exponential tail", {
  reporting_parameters <- list(mean=5, sd=3)
  revision_parameters <- list(tail=list(weight=0.2, mean=12))
  days <- c(4, 6, 10, 30)
  expected <- log(0.8 * exp(f_log_base(days, 4, reporting_parameters)) +
                    0.2 * exp(-(days + 0.5 - 4) / 12))
  expect_equal(log_undetected_days(days, 4, reporting_parameters, revision_parameters),
               expected)
  # the tail makes late detection more likely
  expect_true(log_undetected_days(30, 4, reporting_parameters, revision_parameters) >
                log_undetected_days(30, 4, reporting_parameters))
  # with all weight on the tail the delay is exponential
  expect_equal(log_undetected_days(days, 4, reporting_parameters,
                                   list(tail=list(weight=1, mean=12))),
               -(days + 0.5 - 4) / 12)
})

test_that("log_undetected_days applies multipliers to the rate of detection", {
  reporting_parameters <- list(mean=5, sd=3)
  log_base <- f_log_base(3:8, 3, reporting_parameters)
  rates <- -diff(log_base) # for days 4 to 8

  # multipliers for days 2 to 7
  revision_parameters <- list(
    report_multipliers=list(day_first=2, values=c(9, 9, 0.5, 2, 0, 1.5)))
  # days 4 to 7 have multipliers 0.5, 2, 0, 1.5; day 8 is beyond those given
  expected <- log_base[1] - cumsum(c(0, c(0.5, 2, 0, 1.5, 1) * rates))
  expect_equal(log_undetected_days(3:8, 3, reporting_parameters, revision_parameters),
               expected)
  expect_equal(log_undetected_days(c(8, 5), 3, reporting_parameters, revision_parameters),
               expected[c(6, 3)])
  # no detection on a day with a multiplier of zero
  expect_equal(expected[4], expected[3])
  # on the onset day itself multipliers have no effect
  expect_equal(log_undetected_days(3, 3, reporting_parameters, revision_parameters),
               log_base[1])

  # multipliers of one make no difference
  revision_parameters <- list(report_multipliers=list(day_first=1, values=rep(1, 20)))
  expect_equal(log_undetected_days(3:8, 3, reporting_parameters, revision_parameters),
               log_base)

  # larger multipliers mean earlier detection
  f <- function(multiplier)
    log_undetected_days(8, 3, reporting_parameters,
                        list(report_multipliers=list(day_first=1, values=rep(multiplier, 20))))
  expect_true(f(2) < f(1))
  expect_true(f(0.5) > f(1))
  expect_equal(f(2) - log_base[1], 2 * (f(1) - log_base[1]))
})

test_that("log_undetected_days combines a tail with multipliers", {
  reporting_parameters <- list(mean=5, sd=3)
  tail <- list(weight=0.3, mean=10)
  log_tail <- log_undetected_days(3:6, 3, reporting_parameters, list(tail=tail))
  revision_parameters <- list(
    tail=tail, report_multipliers=list(day_first=4, values=c(2, 0.5, 3)))
  expected <- log_tail[1] - cumsum(c(0, c(2, 0.5, 3) * -diff(log_tail)))
  expect_equal(log_undetected_days(3:6, 3, reporting_parameters, revision_parameters),
               expected)
})

test_that("weekday_of_days cycles through the days of the week", {
  expect_equal(weekday_of_days(1:9), c(1:7, 1, 2))
  expect_equal(weekday_of_days(1:9, weekday_first=6), c(6, 7, 1:7))
  expect_equal(weekday_of_days(c(15, 1), weekday_first=3), c(3, 3))
})

test_that("report_multipliers gives each day the multiplier of the next report", {
  report_days <- c(3, 4, 8, 10)
  report_log_effects <- log(c(2, 0.5, 3, 0.25))
  multipliers <- report_multipliers(report_days, report_log_effects=report_log_effects)
  expect_equal(multipliers$day_first, 3)
  # days 3 to 10: reports on days 3, 4, 8 (covering days 5 to 8) and 10 (9 and 10)
  expect_equal(multipliers$values, c(2, 0.5, 3, 3, 3, 3, 0.25, 0.25))

  # effects of the day of the week: days 3, 4, 8, 10 are days 3, 4, 1, 3 of the week
  weekday_log_effects <- c("1"=log(4), "3"=log(0.5), "4"=log(0.5))
  multipliers <- report_multipliers(report_days, weekday_log_effects=weekday_log_effects)
  expect_equal(multipliers$values, c(0.5, 0.5, 4, 4, 4, 4, 0.5, 0.5))
  # starting the week on a different day
  weekday_log_effects <- c("2"=log(4), "4"=log(0.5), "5"=log(0.5))
  multipliers <- report_multipliers(report_days, weekday_log_effects=weekday_log_effects,
                                    weekday_first=2)
  expect_equal(multipliers$values, c(0.5, 0.5, 4, 4, 4, 4, 0.5, 0.5))

  # both effects multiply
  multipliers <- report_multipliers(report_days,
                                    weekday_log_effects=c("1"=log(4), "3"=log(0.5), "4"=log(0.5)),
                                    report_log_effects=report_log_effects)
  expect_equal(multipliers$values, c(1, 0.25, 12, 12, 12, 12, 0.125, 0.125))

  # neither
  multipliers <- report_multipliers(report_days)
  expect_equal(multipliers$values, rep(1, 8))
})

test_that("update_report_multipliers refreshes multipliers within revision parameters", {
  revision_parameters <- list(eta=0.1, rho=0.2, report_days=c(2, 3, 5))
  updated <- update_report_multipliers(revision_parameters)
  expect_null(updated$report_multipliers)

  revision_parameters$report_log_effects <- log(c(2, 1, 0.5))
  updated <- update_report_multipliers(revision_parameters)
  expect_equal(updated$report_multipliers,
               list(day_first=2, values=c(2, 1, 0.5, 0.5)))
  expect_equal(updated$eta, 0.1)

  # weekday effects with the week starting on day 5
  revision_parameters$weekday_log_effects <- c("1"=log(3), "6"=log(1), "7"=log(1 / 3))
  revision_parameters$weekday_first <- 5
  updated <- update_report_multipliers(revision_parameters)
  # days 2, 3, 5 are days 6, 7, 2... of the week starting at 5: 6, 7, 2
  expect_equal(weekday_of_days(c(2, 3, 5), 5), c(6, 7, 2))
  revision_parameters$weekday_log_effects <- c("2"=log(3), "6"=log(1), "7"=log(1 / 3))
  updated <- update_report_multipliers(revision_parameters)
  expect_equal(updated$report_multipliers$values, c(2, 1 / 3, 1.5, 1.5))

  # removing all effects removes the multipliers
  updated$report_log_effects <- NULL
  updated$weekday_log_effects <- NULL
  expect_null(update_report_multipliers(updated)$report_multipliers)
})

test_that("delay_quantile inverts the probability of remaining undetected", {
  reporting_parameters <- list(mean=5, sd=3)
  expect_equal(delay_quantile(0.9, reporting_parameters),
               qgamma_mean_sd(0.9, 5, 3))
  expect_equal(delay_quantile(0.9, reporting_parameters, list(eta=0.1)),
               qgamma_mean_sd(0.9, 5, 3))

  revision_parameters <- list(tail=list(weight=0.2, mean=15))
  for(p in c(0.5, 0.9, 0.99)) {
    delay <- delay_quantile(p, reporting_parameters, revision_parameters)
    p_undetected <- 0.8 * (1 - pgamma_mean_sd(delay, 5, 3)) + 0.2 * exp(-delay / 15)
    expect_equal(p_undetected, 1 - p, tolerance=1e-4)
  }
  # a tail lengthens the delay by which nearly all cases are detected
  expect_true(delay_quantile(0.99, reporting_parameters, revision_parameters) >
                qgamma_mean_sd(0.99, 5, 3))
  expect_equal(delay_quantile(0, reporting_parameters, revision_parameters), 0)
})

test_that("observation probability with a tail and multipliers matches product
of binomial probabilities when there are no erroneous records", {
  observation_df <- dplyr::tibble(time_reported=1:6,
                                  cases_reported=c(0, 2, 5, 5, 8, 9))
  reporting_parameters <- list(mean=4, sd=2.5)
  cases_true <- 9:15
  l_revision_parameters <- list(
    list(eta=0, rho=0.3, tail=list(weight=0.25, mean=9)),
    list(eta=0, rho=0.3,
         report_multipliers=list(day_first=2, values=c(0.5, 2, 1, 0.2, 1.5))),
    list(eta=0, rho=0.3, tail=list(weight=0.25, mean=9),
         report_multipliers=list(day_first=2, values=c(0.5, 2, 1, 0.2, 1.5))))
  for(revision_parameters in l_revision_parameters) {
    log_s <- log_undetected_days(1:6, 1, reporting_parameters, revision_parameters)
    p_detect <- 1 - exp(diff(log_s))
    expected <- purrr::map_dbl(cases_true, function(case_true) {
      sum(dbinom(diff(observation_df$cases_reported),
                 case_true - observation_df$cases_reported[1:5],
                 p_detect, log=TRUE))
    })
    logp <- observation_process_logp(observation_df, cases_true, 1,
                                     reporting_parameters,
                                     revision_parameters=revision_parameters)
    expect_equal(logp, expected)
  }

  # multipliers of one and a tail of zero weight change nothing
  logp_plain <- observation_process_logp(
    observation_df, cases_true, 1, reporting_parameters,
    revision_parameters=list(eta=0.1, rho=0.3))
  logp_neutral <- observation_process_logp(
    observation_df, cases_true, 1, reporting_parameters,
    revision_parameters=list(eta=0.1, rho=0.3, tail=list(weight=0, mean=9),
                             report_multipliers=list(day_first=1, values=rep(1, 6))))
  expect_equal(logp_neutral, logp_plain)

  # a report in which nothing can be detected cannot add genuine cases
  logp_zero <- observation_process_logp(
    observation_df, cases_true, 1, reporting_parameters,
    revision_parameters=list(eta=0, rho=0.3,
                             report_multipliers=list(day_first=2, values=c(1, 0, 1, 1, 1))))
  expect_equal(logp_zero, rep(-Inf, length(cases_true)))

  # the slower version does not allow for extensions
  expect_error(observation_process_revisions_logp(
    observation_df, cases_true, 1, reporting_parameters,
    list(eta=0.1, rho=0.3, tail=list(weight=0.2, mean=9)), use_compiled=FALSE))
})

test_that("extensions to reporting carry through to the first reporting day
and to misdating", {
  reporting_parameters <- list(mean=4, sd=2.5)
  revision_parameters <- list(eta=0.2, rho=0.3, tail=list(weight=0.25, mean=9),
                              report_multipliers=list(day_first=1, values=c(1, 1, 0.3, 2, 1, 1)))
  # onset on day 2, first reported on day 5
  log_alpha <- erroneous_initial_logp(5, 6, 8, 2, reporting_parameters,
                                      revision_parameters, max_erroneous=3)
  p_undetected <- exp(log_undetected_days(2:5, 2, reporting_parameters, revision_parameters))
  mean_erroneous <- 0.2 * 8 * sum(-diff(p_undetected) * 0.7^(2:0))
  p_detect <- 1 - p_undetected[4] / p_undetected[1]
  probs <- dpois(0:3, mean_erroneous) * dbinom(6 - 0:3, 8, p_detect)
  expect_equal(exp(log_alpha[1, ]), probs / sum(probs))

  # probability of not yet being recorded under the onset day with misdating
  revision_parameters$pi <- 0.3
  revision_parameters$shifts <- c(-1, 1)
  revision_parameters$shift_probs <- c(0.5, 0.5)
  p_newly <- -diff(p_undetected)
  # detected on day 3: only the negative shift is possible; days 4 and 5: both
  p_awaiting <- p_newly[1] * 0.5 * 0.7^2 + p_newly[2] * 0.7 + p_newly[3]
  expect_equal(unrecorded_prob(5, 2, reporting_parameters, revision_parameters),
               p_undetected[4] + 0.3 * p_awaiting)

  # misdated records entered from a neighbouring onset day
  unit <- misdated_entry_unit(c(4, 5, 6), 3, 1, reporting_parameters, revision_parameters)
  p_neighbour <- exp(log_undetected_days(4:6, 2, reporting_parameters, revision_parameters))
  expect_equal(unit$periods, 0.3 * 0.5 * -diff(p_neighbour))
})

test_that("simulated reported cases with a tail and multipliers match the
observation probability", {
  set.seed(50)
  reporting_parameters <- list(mean=3, sd=2)
  revision_parameters <- list(
    eta=0.25, rho=0.5, tail=list(weight=0.3, mean=6),
    report_multipliers=list(day_first=1, values=c(1, 1, 1, 1, 1, 0.3, 2, 1, 0.5)))
  cases_true <- c(0, 0, 0, 0, 4, 0, 0, 0, 0)
  nreplicates <- 20000
  trajectories <- replicate(nreplicates, {
    df <- observed_cases_revisions(cases_true, reporting_parameters,
                                   revision_parameters)
    paste0(df$cases_reported[df$time_onset == 5], collapse="_")
  })
  frequencies <- sort(table(trajectories), decreasing=TRUE) / nreplicates
  l_trajectories <- purrr::map(strsplit(names(frequencies), "_"), as.numeric)
  for(i in 1:12) {
    prob <- exp(observation_process_revisions_logp(
      dplyr::tibble(time_reported=5:9, cases_reported=l_trajectories[[i]]),
      4, 5, reporting_parameters, revision_parameters, max_erroneous=7))
    std_error <- sqrt(prob * (1 - prob) / nreplicates)
    expect_true(abs(frequencies[[i]] - prob) < 4 * std_error)
  }
})

test_that("simulated detection follows the multipliers and the tail", {
  set.seed(51)
  reporting_parameters <- list(mean=3, sd=2)
  cases_true <- c(2000, rep(0, 11))
  f_reported <- function(revision_parameters) {
    df <- observed_cases_revisions(cases_true, reporting_parameters,
                                   revision_parameters, days_max_follow_up=20)
    df$cases_reported[df$time_onset == 1]
  }
  # no detection on days 4 and 5, and faster detection on day 6
  values <- rep(1, 12)
  values[c(4, 5)] <- 0
  values[6] <- 3
  revision_parameters <- list(report_multipliers=list(day_first=1, values=values))
  reported <- f_reported(revision_parameters)
  expect_equal(reported[4], reported[3])
  expect_equal(reported[5], reported[3])
  expect_true(reported[6] > reported[5])
  # proportions detected match the probabilities
  p_detected <- 1 - exp(log_undetected_days(1:12, 1, reporting_parameters, revision_parameters) -
                          log_undetected_days(1, 1, reporting_parameters, revision_parameters))
  expect_equal(reported / 2000, p_detected, tolerance=0.05)

  # with a tail, fewer are detected early and detection continues for longer
  reported_plain <- f_reported(NULL)
  revision_parameters <- list(tail=list(weight=0.5, mean=20))
  reported_tail <- f_reported(revision_parameters)
  expect_true(reported_tail[12] < reported_plain[12])
  p_detected <- 1 - exp(log_undetected_days(1:12, 1, reporting_parameters, revision_parameters) -
                          log_undetected_days(1, 1, reporting_parameters, revision_parameters))
  expect_equal(reported_tail / 2000, p_detected, tolerance=0.05)

  # without extensions the simulation is as before
  set.seed(52)
  reported_1 <- f_reported(list(eta=0.1, rho=0.2))
  set.seed(52)
  reported_2 <- f_reported(list(eta=0.1, rho=0.2, tail=list(weight=0, mean=5),
                                report_multipliers=list(day_first=1, values=rep(1, 12))))
  expect_equal(reported_1, reported_2)
})

test_that("generate_snapshots allows extensions to reporting", {
  set.seed(53)
  days_total <- 28
  Rt_function <- stats::approxfun(1:days_total, rep(1.1, days_total))
  # no reporting at weekends (days 6 and 7 of each week)
  weekday_log_effects <- c("1"=0, "2"=0, "3"=0, "4"=0, "5"=0, "6"=-Inf, "7"=-Inf)
  simulated <- simulate_report_multipliers(1:days_total, weekday_log_effects)
  expect_null(simulated$report_log_effects)
  expect_equal(simulated$report_multipliers$values[1:8], c(1, 1, 1, 1, 1, 0, 0, 1))
  df <- generate_snapshots(
    days_total, Rt_function, list(mean=5, sd=3), list(mean=4, sd=2),
    revision_parameters=list(report_multipliers=simulated$report_multipliers))
  changes <- df %>%
    dplyr::group_by(time_onset) %>%
    dplyr::mutate(change=c(0, diff(cases_reported))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(weekday=weekday_of_days(time_reported))
  expect_equal(sum(changes$change[changes$weekday %in% c(6, 7)]), 0)
  expect_true(sum(changes$change[changes$weekday == 1]) > 0)
})

test_that("simulate_report_multipliers draws random effects of reports", {
  set.seed(54)
  report_days <- seq(1, 400, 2)
  simulated <- simulate_report_multipliers(report_days, report_sd=0.5)
  expect_equal(length(simulated$report_log_effects), 200)
  expect_equal(sd(simulated$report_log_effects), 0.5, tolerance=0.15)
  expect_equal(mean(simulated$report_log_effects), 0, tolerance=0.1)
  expect_equal(simulated$report_multipliers$day_first, 1)
  expect_equal(length(simulated$report_multipliers$values), 399)
  # each day has the multiplier of the next report
  expect_equal(simulated$report_multipliers$values[c(1, 2, 3)],
               exp(simulated$report_log_effects[c(1, 2, 2)]))
  # combined with effects of the day of the week
  weekday_log_effects <- stats::setNames(log(c(1, 2, 3, 4, 5, 6, 7)), 1:7)
  both <- simulate_report_multipliers(1:7, weekday_log_effects, report_sd=0)
  expect_equal(both$report_multipliers$values, 1:7)
})

test_that("sample_report_sd draws from inverse-gamma conditional posterior", {
  set.seed(55)
  report_log_effects <- c(0.3, -0.2, 0.5, -0.4, 0.1, 0, -0.6, 0.2)
  prior <- list(report_sd_shape=3, report_sd_scale=0.2)
  draws <- purrr::map_dbl(1:5000, ~sample_report_sd(report_log_effects, prior))
  shape <- 3 + 4
  scale <- 0.2 + sum(report_log_effects^2) / 2
  # mean of an inverse-gamma distribution for the variance
  expect_equal(mean(draws^2), scale / (shape - 1), tolerance=0.05)
  expect_true(all(draws > 0))
  # more variable effects give a larger standard deviation
  draws_large <- purrr::map_dbl(1:2000, ~sample_report_sd(3 * report_log_effects, prior))
  expect_true(mean(draws_large) > mean(draws))
})

f_reporting_data <- function(seed, days_total=10, weekday_log_effects=NULL,
                             report_sd=0.6, tail=NULL, eta=0.02) {
  set.seed(seed)
  simulated <- simulate_report_multipliers(1:days_total, weekday_log_effects, report_sd)
  revision_parameters <- list(eta=eta, rho=0.3, tail=tail,
                              report_multipliers=simulated$report_multipliers)
  Rt_function <- function(t) rep(1.2, length(t))
  df <- generate_snapshots(days_total, Rt_function, list(mean=3, sd=2),
                           list(mean=3, sd=2),
                           initial_parameters=list(mean=20, length=5),
                           revision_parameters=revision_parameters) %>%
    dplyr::mutate(reporting_piece_index=1)
  list(df=df, simulated=simulated, revision_parameters=revision_parameters,
       r_params=dplyr::tibble(reporting_piece_index=1, mean=3, sd=2))
}

test_that("sample_report_effects targets the conditional posterior of the
effects of reports", {
  # three onset days, each reported on days 1 to 3, so that the reports on
  # days 2 and 3 have effects which matter
  l <- f_reporting_data(56, days_total=3)
  df <- l$df
  grid <- report_log_effect_grid()
  report_days <- 1:3
  report_sd <- 0.6
  max_erroneous <- max_downward_revision(df) + 2
  revision_parameters <- list(eta=0.02, rho=0.3, report_days=report_days,
                              report_log_effects=rep(0, 3), report_sd=report_sd)
  revision_parameters <- update_report_multipliers(revision_parameters)

  # exact posterior of the effects of the second and third reports (the first
  # report has no effect, since no case is detected on its onset day)
  f_logl <- function(log_effects) {
    revision_tmp <- revision_parameters
    revision_tmp$report_log_effects <- log_effects
    revision_tmp <- update_report_multipliers(revision_tmp)
    observation_process_all_times_logp(df, l$r_params, revision_tmp, max_erroneous)
  }
  log_prior <- dnorm(grid, 0, report_sd, log=TRUE)
  index <- seq(1, length(grid), 3) # a coarser grid suffices for the means
  m_logp <- matrix(NA, length(index), length(index))
  for(a in seq_along(index))
    for(b in seq_along(index))
      m_logp[a, b] <- f_logl(c(0, grid[index[a]], grid[index[b]])) +
        log_prior[index[a]] + log_prior[index[b]]
  m_p <- exp(m_logp - matrixStats::logSumExp(m_logp))
  mean_exact <- c(sum(rowSums(m_p) * grid[index]), sum(colSums(m_p) * grid[index]))
  sd_exact <- sqrt(c(sum(rowSums(m_p) * grid[index]^2), sum(colSums(m_p) * grid[index]^2)) -
                     mean_exact^2)

  set.seed(57)
  n_sweeps <- 600
  logl_current <- f_logl(revision_parameters$report_log_effects)
  m_draws <- matrix(NA, n_sweeps, 3)
  n_accepted <- 0
  for(i in 1:n_sweeps) {
    res <- sample_report_effects(df, l$r_params, revision_parameters,
                                 logl_current, max_erroneous, n_blocks=1)
    revision_parameters <- res$revision_parameters
    logl_current <- res$logl
    n_accepted <- n_accepted + res$n_accepted
    m_draws[i, ] <- revision_parameters$report_log_effects
  }
  # values stay on the grid and the log-probability carried is correct
  expect_true(all(round(m_draws, 8) %in% round(grid, 8)))
  expect_equal(logl_current, f_logl(revision_parameters$report_log_effects))
  expect_equal(revision_parameters$report_multipliers$values,
               exp(revision_parameters$report_log_effects))
  # proposals from the approximation are usually accepted
  expect_true(n_accepted / n_sweeps > 0.5)
  z_scores <- (colMeans(m_draws[, 2:3]) - mean_exact) / (sd_exact / sqrt(n_sweeps / 3))
  expect_true(all(abs(z_scores) < 3.5))
})

test_that("metropolis steps for the tail and the day of the week return
values with their log-probability", {
  l <- f_reporting_data(58, days_total=16, tail=list(weight=0.2, mean=8),
                        weekday_log_effects=stats::setNames(log(c(1, 1, 1, 1, 1, 0.2, 0.2) * 0.2^(-2 / 7)), 1:7))
  df <- l$df
  max_erroneous <- max_downward_revision(df) + 2
  revision_parameters <- list(eta=0.02, rho=0.3, tail=list(weight=0.1, mean=10),
                              report_days=1:16, weekday_first=1,
                              weekday_log_effects=stats::setNames(rep(0, 7), 1:7))
  revision_parameters <- update_report_multipliers(revision_parameters)
  f_logl <- function(revision_tmp)
    observation_process_all_times_logp(df, l$r_params, revision_tmp, max_erroneous)
  prior <- list(tail_weight_alpha=1, tail_weight_beta=3, tail_mean_mu=15,
                tail_mean_sigma=10, weekday_sd=1)
  steps <- list(tail_weight_step=0.3, tail_mean_step=0.2, weekday_step=0.2)
  logl <- f_logl(revision_parameters)
  expect_true(is.finite(logl))

  set.seed(59)
  # tail: always accepted if the current log-probability is -Inf
  res <- metropolis_step_tail(df, l$r_params, revision_parameters, -Inf, prior, steps,
                              max_erroneous)
  expect_true(res$revision_parameters$tail$weight != 0.1)
  expect_true(res$revision_parameters$tail$weight > 0 & res$revision_parameters$tail$weight < 1)
  expect_true(res$revision_parameters$tail$mean > 0)
  expect_equal(res$logl, f_logl(res$revision_parameters))
  # never accepted if it is Inf
  res <- metropolis_step_tail(df, l$r_params, revision_parameters, Inf, prior, steps,
                              max_erroneous)
  expect_equal(res$revision_parameters$tail, revision_parameters$tail)
  expect_equal(res$logl, Inf)

  # day of the week
  res <- metropolis_step_weekday(df, l$r_params, revision_parameters, -Inf, prior, steps,
                                 max_erroneous)
  log_effects <- res$revision_parameters$weekday_log_effects
  expect_equal(names(log_effects), as.character(1:7))
  expect_equal(sum(log_effects), 0)
  expect_true(any(log_effects != 0))
  expect_equal(res$logl, f_logl(res$revision_parameters))
  expect_equal(res$revision_parameters$report_multipliers,
               report_multipliers(1:16, log_effects))
  res <- metropolis_step_weekday(df, l$r_params, revision_parameters, Inf, prior, steps,
                                 max_erroneous)
  expect_equal(res$revision_parameters$weekday_log_effects,
               revision_parameters$weekday_log_effects)

  # a single day of the week has no free effect
  revision_single <- revision_parameters
  revision_single$weekday_log_effects <- c("1"=0)
  res <- metropolis_step_weekday(df, l$r_params, revision_single, logl, prior, steps,
                                 max_erroneous)
  expect_equal(res$logl, logl)

  # repeated steps move the effects of the weekend below those of weekdays
  set.seed(60)
  current <- revision_parameters
  logl_current <- logl
  for(i in 1:150) {
    res <- metropolis_step_weekday(df, l$r_params, current, logl_current, prior, steps,
                                   max_erroneous)
    current <- res$revision_parameters
    logl_current <- res$logl
  }
  expect_true(logl_current > logl)
  expect_true(mean(current$weekday_log_effects[c("6", "7")]) <
                mean(current$weekday_log_effects[as.character(1:5)]) - 0.5)
})

test_that("accept_reject rejects proposals whose log-probability is undefined", {
  res <- accept_reject(-Inf, -Inf, "current", "proposed")
  expect_equal(res$parameter, "current")
  res <- accept_reject(-5, NaN, "current", "proposed")
  expect_equal(res$parameter, "current")
  expect_equal(res$logp, -5)
  res <- accept_reject(-Inf, -5, "current", "proposed")
  expect_equal(res$parameter, "proposed")
})

test_that("sample_cases_history extends the period of incomplete reporting
if the delay has a tail", {
  set.seed(61)
  days_total <- 40
  r_list <- list(mean=3, sd=2)
  Rt_function <- stats::approxfun(1:days_total, rep(1.1, days_total))
  s_params <- list(mean=5, sd=3)
  df <- generate_snapshots(days_total, Rt_function, s_params, r_list) %>%
    dplyr::mutate(reporting_piece_index=1)
  cases_current <- df %>%
    dplyr::select(time_onset, cases_true) %>%
    unique()
  cases_current$cases_true <- cases_current$cases_true + 7
  df <- df %>% dplyr::select(-cases_true)
  r_params <- dplyr::tibble(reporting_piece_index=1, mean=3, sd=2)
  latest <- purrr::map_dbl(1:days_total, ~dplyr::last(df$cases_reported[df$time_onset == .]))

  f_n_sampled <- function(revision_parameters) {
    df_est <- sample_cases_history(df, 2000, Rt_function, s_params, r_params,
                                   cases_true_current=cases_current,
                                   revision_parameters=revision_parameters)
    estimated <- df_est %>%
      dplyr::select(time_onset, cases_estimated) %>%
      unique()
    # onset days whose counts are fixed equal the latest reported
    sum(estimated$cases_estimated != latest | (1:days_total) > days_total - 3)
  }
  n_plain <- f_n_sampled(list(eta=0.01, rho=0.3))
  n_tail <- f_n_sampled(list(eta=0.01, rho=0.3, tail=list(weight=0.3, mean=10)))
  expect_true(n_tail > n_plain)
})

f_mcmc_reporting <- function(l, initial_revision_parameters, niter=3, ...) {
  df <- l$df
  initial_cases_true <- df %>%
    dplyr::select(time_onset, cases_true) %>%
    unique()
  days_total <- max(df$time_onset)
  snapshot_with_Rt_index_df <- df %>%
    dplyr::select(-cases_true, -reporting_piece_index) %>%
    dplyr::mutate(Rt_index=ifelse(time_onset <= days_total / 2, 1, 2))
  priors <- list(Rt=list(shape=1, rate=1),
                 reporting=list(mean_mu=5, mean_sigma=10, sd_mu=3, sd_sigma=5),
                 revisions=list(eta_mu=0.05, eta_sigma=0.05, rho_alpha=1, rho_beta=1,
                                pi_alpha=1, pi_beta=4),
                 max_erroneous=max_downward_revision(df) + 2,
                 max_cases=3000)
  mcmc(niterations=niter,
       snapshot_with_Rt_index_df,
       priors,
       list(mean=3, sd=2),
       initial_cases_true,
       list(mean=3, sd=2),
       dplyr::tibble(Rt_index=c(1, 2), Rt=c(1.2, 1.2)),
       print_to_screen=FALSE,
       initial_revision_parameters=initial_revision_parameters,
       ...)
}

test_that("mcmc works with each extension to reporting switched on", {
  l <- f_reporting_data(62, days_total=24, tail=list(weight=0.2, mean=8))
  niter <- 3
  base <- list(eta=0.02, rho=0.3)

  # none: no additional output
  set.seed(63)
  res_none <- f_mcmc_reporting(l, base, niter)
  expect_false("reporting_effects" %in% names(res_none))
  expect_false("report_effects" %in% names(res_none))

  # explicitly switched off is the same as not mentioned
  set.seed(63)
  res_off <- f_mcmc_reporting(l, c(base, list(weekday=FALSE, report_effects=FALSE)), niter)
  expect_equal(res_off$cases, res_none$cases)
  expect_equal(res_off$Rt, res_none$Rt)

  # tail only
  res <- f_mcmc_reporting(l, c(base, list(tail=list(weight=0.1, mean=10))), niter)
  expect_equal(nrow(res$reporting_effects), niter)
  expect_true(all(res$reporting_effects$tail_weight > 0 & res$reporting_effects$tail_weight < 1))
  expect_true(all(res$reporting_effects$tail_mean > 0))
  expect_true(all(is.na(res$reporting_effects$report_sd)))
  expect_false("report_effects" %in% names(res))
  expect_false(any(grepl("weekday", colnames(res$reporting_effects))))

  # day of the week only: the week starts on day 3
  res <- f_mcmc_reporting(l, c(base, list(weekday=TRUE, weekday_first=3)), niter)
  weekday_columns <- paste0("weekday_", 1:7)
  expect_true(all(weekday_columns %in% colnames(res$reporting_effects)))
  products <- apply(as.matrix(res$reporting_effects[, weekday_columns]), 1, prod)
  expect_equal(products, rep(1, niter))
  expect_true(all(is.na(res$reporting_effects$tail_weight)))

  # random effects of reports only
  res <- f_mcmc_reporting(l, c(base, list(report_effects=TRUE)), niter)
  expect_equal(nrow(res$report_effects), niter * 24)
  expect_true(all(c("iteration", "time_reported", "effect", "chain") %in%
                    colnames(res$report_effects)))
  expect_true(all(res$report_effects$effect > 0))
  expect_true(all(res$reporting_effects$report_sd > 0))
  expect_true(length(unique(res$reporting_effects$report_sd)) > 1)

  # all three, with misdating, over two chains
  initial <- c(base, list(pi=0.05, shifts=c(-1, 1), shift_probs=c(0.5, 0.5),
                          tail=list(weight=0.1, mean=10), weekday=TRUE,
                          report_effects=TRUE))
  res <- f_mcmc_reporting(
    l, initial, niter, nchains=2,
    revision_metropolis_parameters=list(eta_step=0.3, rho_step=0.3, pi_step=0.3,
                                        tail_weight_step=0.2, tail_mean_step=0.2,
                                        weekday_step=0.05, report_blocks=3))
  expect_equal(nrow(res$reporting_effects), 2 * niter)
  expect_equal(sort(unique(res$report_effects$chain)), c(1, 2))
  expect_equal(sum(is.na(res$cases$cases_true)), 0)
  expect_true(all(is.finite(res$reporting$mean)))

  # holding the tail fixed
  res <- f_mcmc_reporting(
    l, c(base, list(tail=list(weight=0.1, mean=10))), niter,
    revision_metropolis_parameters=list(eta_step=0.3, rho_step=0.3,
                                        tail_weight_step=0, tail_mean_step=0))
  expect_equal(res$reporting_effects$tail_weight, rep(0.1, niter))
  expect_equal(res$reporting_effects$tail_mean, rep(10, niter))
})

test_that("mcmc recovers an effect of the day of the week", {
  # reports at the weekend bring in cases at a tenth of the usual rate
  weekday_log_effects <- stats::setNames(log(c(1, 1, 1, 1, 1, 0.1, 0.1) * 0.1^(-2 / 7)), 1:7)
  l <- f_reporting_data(64, days_total=28, weekday_log_effects=weekday_log_effects,
                        report_sd=0, eta=0.005)
  set.seed(65)
  res <- f_mcmc_reporting(
    l, list(eta=0.005, rho=0.3, weekday=TRUE), niter=120,
    revision_metropolis_parameters=list(eta_step=0.3, rho_step=0.3, weekday_step=0.1))
  effects <- res$reporting_effects %>%
    dplyr::filter(iteration > 60)
  weekend <- mean(c(effects$weekday_6, effects$weekday_7))
  weekdays <- mean(c(effects$weekday_1, effects$weekday_2, effects$weekday_3,
                     effects$weekday_4, effects$weekday_5))
  expect_true(weekend < 0.5 * weekdays)
})

test_that("with a cache the probabilities of remaining undetected are the same", {
  before <- list(lag_probabilities=diff(stats::pgamma(0:30, 4, 0.4)))
  multipliers <- report_multipliers(c(12, 15, 19, 26, 33, 40),
                                    weekday_log_effects=stats::setNames(c(0.3, -0.1, -0.2, 0.1, 0, 0.2, -0.3), 1:7),
                                    report_log_effects=c(0.2, -0.3, 0.1, 0, 0.4, -0.2))
  backlog <- initial_backlog(40, c(12, 19, 26, 40), backlog_settings(TRUE))
  settings <- list(
    list(),
    list(before=before),
    list(report_multipliers=multipliers),
    list(before=before, report_multipliers=multipliers),
    list(tail=list(weight=0.2, mean=15)),
    list(before=before, tail=list(weight=0.2, mean=15),
         report_multipliers=multipliers),
    list(before=before, report_multipliers=multipliers, backlog=backlog))
  cache <- new_lag_cache()
  for(revision_parameters in settings)
    for(onset in c(3, 10, 25))
      for(parameters in list(list(mean=8, sd=6), list(mean=12, sd=5))) {
        days <- onset + c(0, 2, 5, 9, 20)
        cached <- revision_parameters
        cached$lag_cache <- cache
        expect_identical(log_undetected_days(days, onset, parameters, cached),
                         log_undetected_days(days, onset, parameters, revision_parameters))
      }
  # one set of values for each delay and what alters it, without the
  # multipliers or the backlog
  expect_equal(length(ls(cache)), 2 * 4)
})
