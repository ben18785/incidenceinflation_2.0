# S7 of the sensitivity analyses: the rate at which wrong records are entered,
# in proportion to the genuine cases being reported, varies over the intervals
# between reports, following a random walk on the log scale
# (deaths-wrong-rate.R). The expected values are worked out here day by day
# from the model, not from the package's code.

# the interval of each day: the first report, of all, made on or after it,
# or the last for days after the last report
walk_interval_daily <- function(days, report_days)
  vapply(days, function(d) {
    j <- which(report_days >= d)[1]
    if(is.na(j)) length(report_days) else j
  }, numeric(1))

# the part of the unit of entry of each report from the days of each
# interval, day by day: the genuine cases newly reported on day k, each
# still present at the report with probability (1 - rho)^(t_j - k)
walk_split_daily <- function(days, day_onset, reporting, revision, report_days) {
  rho <- revision$rho
  undetected <- function(day) exp(log_undetected_days(day, day_onset, reporting, revision))
  before <- c(day_onset, days[-length(days)])
  split <- matrix(0, length(days), length(report_days))
  for(j in seq_along(days)) {
    if(days[j] <= before[j])
      next
    for(k in (before[j] + 1):days[j]) {
      J <- walk_interval_daily(k, report_days)
      split[j, J] <- split[j, J] +
        (undetected(k - 1) - undetected(k)) * (1 - rho)^(days[j] - k)
    }
  }
  split
}

walk_revision <- function(log_rate, report_days, rho=0.3, tau=0.4) {
  list(rho=rho, wrong_walk=list(log_rate=log_rate, report_days=report_days, tau=tau,
                                mu=log(0.07), sigma=1.2, tau_scale=0.5, tau_step=0.3))
}

# reports with wrong records in proportion to the genuine cases and a
# correcting report, with the state of the model of deaths after the paths
# are drawn with the revision parameters given
walk_example <- function(f_revision, seed=1, n_days=12, max_erroneous=6) {
  set.seed(seed)
  cases_true <- stats::rpois(n_days, 8)
  sim <- simulate_cases_deaths(
    cases_true, 3:(n_days + 6), list(mean=3, sd=2), 0.3, list(mean=5, sd=3),
    wrong=list(eta=0.15, rho=0.3, delta=0.4),
    correction=list(time_reported=n_days, onset_times=1:4, wrong=0.1, missing=0.03))
  parts <- split_corrected_reports(sim$snapshots, n_days, 1:4)
  settings <- deaths_settings(n_days)
  settings$n_parameter_steps <- 0
  largest <- parts$corrected %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(cases_true=max(cases_reported) + 2)
  state <- initial_deaths_state(parts$corrected, largest, settings)
  state$uncorrected <- initial_uncorrected(parts$uncorrected, state, settings)
  reporting <- dplyr::tibble(mean=rep(3, n_days), sd=2)
  report_days <- sort(unique(c(parts$corrected$time_reported,
                               parts$uncorrected$time_reported)))
  revision <- f_revision(report_days)
  state <- deaths_paths_step(state, largest$cases_true, reporting, revision,
                             max_erroneous)
  state <- deaths_uncorrected_step(state, reporting, revision, max_erroneous,
                                   max_cases=60)
  list(state=state, cases_true=largest$cases_true, reporting=reporting,
       revision=revision, settings=settings, max_erroneous=max_erroneous,
       report_days=report_days)
}

test_that("C1: the setting of a walk of the rate of wrong records is checked", {
  settings <- wrong_rate_settings(list(wrong_rate="walk"))
  expect_equal(settings, list(form="walk", eta_walk_mu=log(0.07), eta_walk_sigma=1.2,
                              tau_eta_scale=0.5))
  settings <- wrong_rate_settings(list(wrong_rate="walk", eta_walk_sigma=0.5,
                                       tau_eta_scale=0.2))
  expect_equal(settings$eta_walk_sigma, 0.5)
  expect_equal(settings$tau_eta_scale, 0.2)
  expect_error(wrong_rate_settings(list(wrong_rate=c("constant", "walk"))),
               "cannot be combined")
  expect_error(wrong_rate_settings(list(wrong_rate="walk", eta_walk_sigma=0)), "positive")
  expect_error(wrong_rate_settings(list(wrong_rate="walk", tau_eta_scale=-1)), "positive")
  expect_error(wrong_rate_settings(list(wrong_rate="walk", eta_walk_mu=Inf)), "finite")
  expect_error(wrong_rate_settings(list(wrong_rate="walk", eta_walk_mu=c(1, 2))), "single")
  expect_error(wrong_rate_settings(list(wrong_rate="walk"), is_deaths=FALSE),
               "model of deaths")
  # the starting values
  report_days <- c(9, 5, 7, 12, 7)
  revision <- initial_wrong_rate(list(eta=0.05, rho=0.3), settings, report_days)
  expect_null(revision$eta)
  expect_equal(wrong_rate_form(revision), "walk")
  expect_true(is_wrong_records(revision))
  walk <- revision$wrong_walk
  expect_equal(walk$report_days, c(5, 7, 9, 12))
  expect_equal(walk$log_rate, rep(log(0.05), 4))
  expect_equal(walk[c("tau", "mu", "sigma", "tau_scale", "tau_step")],
               list(tau=0.2, mu=log(0.07), sigma=0.5, tau_scale=0.2, tau_step=0.3))
  revision <- initial_wrong_rate(list(eta=0.05, rho=0.3, eta_walk=c(0.1, 0.2, 0.3, 0.4),
                                      tau_eta=0.6), settings, report_days, tau_step=0.1)
  expect_equal(names(revision), c("rho", "wrong_walk"))
  expect_equal(revision$wrong_walk$log_rate, log(c(0.1, 0.2, 0.3, 0.4)))
  expect_equal(revision$wrong_walk$tau, 0.6)
  expect_equal(revision$wrong_walk$tau_step, 0.1)
  expect_error(initial_wrong_rate(list(eta=0.05, rho=0.3, eta_walk=c(0.1, 0.2)),
                                  settings, report_days), "one for each of the 4")
  expect_error(initial_wrong_rate(list(eta=-0.05, rho=0.3), settings, report_days),
               "positive")
  expect_error(initial_wrong_rate(list(eta=0.05, rho=0.3, tau_eta=0), settings,
                                  report_days), "positive")
})

test_that("S7.1: the interval of each day is that of the first report of all on or after it", {
  report_days <- c(5, 7, 9, 12)
  days <- 1:16
  expect_equal(wrong_walk_interval(days, report_days),
               walk_interval_daily(days, report_days))
  expect_equal(wrong_walk_interval(c(1, 5, 6, 12, 13, 40), report_days),
               c(1, 1, 2, 4, 4, 4))
  # the intervals of a fit are those of the reports of every onset time,
  # including those made before a correction
  set.seed(13)
  cases_true <- stats::rpois(20, 8)
  sim <- simulate_cases_deaths(
    cases_true, c(seq(4, 26, 3), 27), list(mean=3, sd=2), 0.3, list(mean=5, sd=3),
    wrong=list(eta=0.1, rho=0.3, delta=0.4),
    correction=list(time_reported=19, onset_times=1:4, wrong=0.1))
  parts <- split_corrected_reports(sim$snapshots, 19, 1:4)
  # an uncorrected series reported on a day no other onset time is
  parts$uncorrected <- dplyr::bind_rows(
    parts$uncorrected, dplyr::mutate(
      dplyr::filter(parts$uncorrected, time_onset == 2, time_reported == 16),
      time_reported=17))
  parts$uncorrected <- dplyr::arrange(parts$uncorrected, time_onset, time_reported)
  fit <- mcmc_single(
    niterations=1,
    snapshot_with_Rt_index_df=parts$corrected %>%
      dplyr::select(time_onset, time_reported, cases_reported) %>%
      dplyr::mutate(Rt_index=1),
    priors=list(Rt=list(shape=1, rate=1),
                reporting=list(mean_mu=3, mean_sigma=3, sd_mu=2, sd_sigma=2),
                max_cases=150, Rt_changepoints=list(),
                reporting_changepoints=list(min_length=14),
                revisions=list(wrong_rate="walk", rho_alpha=1, rho_beta=1),
                max_erroneous=6),
    serial_parameters=list(mean=5, sd=3),
    initial_cases_true=parts$corrected %>%
      dplyr::group_by(time_onset) %>%
      dplyr::summarise(cases_true=max(cases_reported) + 2),
    initial_reporting_parameters=list(mean=3, sd=2),
    initial_Rt=dplyr::tibble(Rt_index=1, Rt=1),
    print_to_screen=FALSE,
    initial_revision_parameters=list(eta=0.1, rho=0.3),
    return_state=TRUE,
    deaths=list(snapshot=parts$corrected %>%
                  dplyr::select(time_onset, time_reported, deaths_reported),
                uncorrected=parts$uncorrected))
  walk <- fit$state$variables$revision_current$wrong_walk
  expect_equal(walk$report_days, sort(unique(c(parts$corrected$time_reported,
                                               parts$uncorrected$time_reported))))
  expect_true(17 %in% walk$report_days)
  expect_equal(colnames(fit$revisions),
               c("eta", "rho", "pi", "iteration", "tau_eta",
                 paste0("eta_", walk$report_days)))
  expect_true(is.na(fit$revisions$eta))
})

test_that("S7.2: the split units sum to those of the main model and are those of each interval, day by day", {
  reporting <- list(mean=4, sd=3)
  report_days <- c(5, 7, 9, 12)
  revision <- walk_revision(rep(log(0.1), 4), report_days)
  revision$tail <- list(weight=0.1, mean=20)
  cases <- list(
    every=c(5, 7, 9, 12),
    # left out of the reports of days 7 and 9
    left_out=c(5, 12),
    # first reported at a later report, as after a correction
    corrected=c(9, 12),
    # and on the onset day itself
    onset_day=c(3, 5, 12))
  for(name in names(cases)) {
    days <- cases[[name]]
    units <- wrong_walk_units(days, 3, reporting, revision)
    expect_equal(units$unit,
                 wrong_entry_units(days, 3, reporting, list(eta=0.1, rho=0.3,
                                                            tail=revision$tail))$unit,
                 info=name)
    expect_equal(rowSums(units$split), units$unit, info=name)
    expect_equal(units$split, walk_split_daily(days, 3, reporting, revision, report_days),
                 info=name)
    expect_equal(units$w, 1 - 0.7^c(0, diff(days)), info=name)
    expect_equal(deaths_wrong_units(days, 3, reporting, revision), units, info=name)
  }
  units <- wrong_walk_units(c(5, 12), 3, reporting, revision)
  # the report after a gap spreads over the intervals of the reports missed
  expect_true(all(units$split[2, 2:4] > 0))
  expect_equal(units$split[2, 1], 0)
  units <- wrong_walk_units(c(9, 12), 3, reporting, revision)
  expect_true(all(units$split[1, 1:3] > 0))
  expect_equal(units$split[1, 4], 0)
  # a report after the last of the walk takes the last interval
  units <- wrong_walk_units(c(5, 15), 3, reporting, revision)
  expect_equal(units$split, walk_split_daily(c(5, 15), 3, reporting, revision, report_days))
  expect_true(units$split[2, 4] > 0)
})

test_that("S7.3: with every rate of the walk equal to eta it is the main model", {
  example <- walk_example(function(report_days)
    walk_revision(rep(log(0.15), length(report_days)), report_days))
  state <- example$state
  walk <- example$revision
  main <- list(eta=0.15, rho=0.3)
  set.seed(14)
  state_main <- deaths_paths_step(state, example$cases_true, example$reporting, main,
                                  example$max_erroneous)
  set.seed(14)
  state_walk <- deaths_paths_step(state, example$cases_true, example$reporting, walk,
                                  example$max_erroneous)
  expect_equal(purrr::map(state_walk$wrong, "lambda"),
               purrr::map(state_main$wrong, "lambda"), tolerance=1e-10)
  # with the same lambda the same paths are drawn
  expect_equal(state_walk$paths, state_main$paths)
  # the log-probability of the reports of every onset time, as the forward
  # filter gives it
  for(k in seq_along(state$onset_times)) {
    reports <- state$reports[[k]]
    f_logp <- function(s) joint_path_forward(
      reports$cases_reported, reports$deaths_reported, s$deaths_true[k],
      example$cases_true[k] - s$deaths_true[k],
      path_probabilities(s$survivals[[k]], 0), s$wrong[[k]])$logp
    expect_equal(f_logp(state_walk), f_logp(state_main), tolerance=1e-10)
  }
  # the log-probability given the paths, over every series, and what the
  # updates of the true cases and of rho use
  expect_equal(deaths_reporting_logp(state_main, example$cases_true, example$reporting, walk),
               deaths_reporting_logp(state_main, example$cases_true, example$reporting, main),
               tolerance=1e-10)
  for(k in c(2, 7))
    expect_equal(deaths_wrong_mean(state_walk, k, walk),
                 deaths_wrong_mean(state_main, k, main), tolerance=1e-10)
  priors <- list(rho_alpha=2, rho_beta=3)
  for(rho in c(0.2, 0.5))
    expect_equal(
      deaths_wrong_rho_logp(rho, deaths_wrong_series(state_walk, example$cases_true,
                                                     example$reporting), walk, priors),
      deaths_wrong_rho_logp(rho, deaths_wrong_series(state_main, example$cases_true,
                                                     example$reporting), main, priors),
      tolerance=1e-10)
})

test_that("S7.4: the log-likelihood of the walk is the change in the log-probability of the reports", {
  example <- walk_example(function(report_days) walk_revision(
    log(seq(0.05, 0.3, length.out=length(report_days))), report_days), seed=15)
  state <- example$state
  revision <- example$revision
  series <- deaths_wrong_series(state, example$cases_true, example$reporting)
  data <- wrong_walk_data(series)
  log_rate <- revision$wrong_walk$log_rate
  # the means of the entries, day by day
  means <- unlist(lapply(seq_along(series$paths), function(k) {
    days <- series$reports[[k]]$time_reported
    split <- walk_split_daily(days, series$onset_times[k],
                              delay_of_onset(series$reporting_by_onset, k), revision,
                              revision$wrong_walk$report_days)
    series$cases_recorded[k] * as.vector(split %*% exp(log_rate))
  }))
  entered <- unlist(lapply(series$paths, function(p)
    p$wrong_entered_alive + p$wrong_entered_dead))
  is_used <- means > 0 | entered > 0
  expect_equal(data$entered, entered[is_used])
  expect_equal(as.vector(data$exposure %*% exp(log_rate)), means[is_used], tolerance=1e-10)
  expect_equal(wrong_walk_logl(log_rate, data),
               sum(stats::dpois(entered, means, log=TRUE)), tolerance=1e-10)
  # C3: changing one rate changes the full log-probability of the reports,
  # given the paths, by as much as the log-likelihood of the walk
  f_full <- function(x) {
    proposed <- revision
    proposed$wrong_walk$log_rate <- x
    deaths_reporting_logp(state, example$cases_true, example$reporting, proposed)
  }
  for(J in c(1, 3, length(log_rate))) {
    x <- log_rate
    x[J] <- x[J] + 0.4
    expect_equal(wrong_walk_logl(x, data) - wrong_walk_logl(log_rate, data),
                 f_full(x) - f_full(log_rate), tolerance=1e-10)
  }
})

test_that("S7.4: the update of the true and fatal cases together keeps the entries of wrong records", {
  example <- walk_example(function(report_days) walk_revision(
    log(seq(0.05, 0.3, length.out=length(report_days))), report_days), seed=16)
  state <- example$state
  revision <- example$revision
  k <- 9
  path <- state$paths[[k]]
  probabilities <- path_probabilities(state$survivals[[k]], 0)
  days <- state$reports[[k]]$time_reported
  split <- walk_split_daily(days, state$onset_times[k], list(mean=3, sd=2), revision,
                            revision$wrong_walk$report_days)
  w <- 1 - (1 - 0.3)^c(0, diff(days))
  entered <- path$wrong_entered_alive + path$wrong_entered_dead
  f_full <- function(cases, fatal) {
    sum(dbinom_log_recycled(path$fatal_new, fatal - path$fatal_before,
                            probabilities$p_fatal)) +
      sum(dbinom_log_recycled(path$survivors_new, cases - fatal - path$survivors_before,
                              probabilities$p)) +
      stats::dbinom(fatal, cases, state$cfr[k], log=TRUE) +
      sum(stats::dpois(entered, cases * as.vector(split %*% exp(revision$wrong_walk$log_rate)),
                       log=TRUE)) +
      sum(dbinom_log_recycled(path$wrong_removed, path$wrong_before, w))
  }
  f_scratch <- function(cases) {
    logp <- vapply(0:cases, function(fatal) f_full(cases, fatal), numeric(1))
    max(logp) + log(sum(exp(logp - max(logp))))
  }
  f_update <- deaths_cases_logp_function(state, state$onset_times[k], revision)
  values <- attr(f_update, "min_cases") + 0:8
  update <- f_update(values)
  scratch <- vapply(values, f_scratch, numeric(1))
  expect_true(all(is.finite(update)))
  expect_equal(update - update[1], scratch - scratch[1], tolerance=1e-10)
})

test_that("S7.4 (C4): the elliptical slice update of the walk draws from its posterior, found on a grid", {
  walk <- list(log_rate=log(c(0.1, 0.1, 0.1)), report_days=c(7, 14, 21), tau=0.5,
               mu=log(0.1), sigma=0.8, tau_scale=0.5, tau_step=0.3)
  data <- list(entered=c(3, 1, 2),
               exposure=rbind(c(20, 0, 0), c(0, 15, 5), c(0, 0, 25)))
  # the posterior on a grid of the three log rates
  grid <- seq(-6, 1, length.out=71)
  x <- as.matrix(expand.grid(grid, grid, grid))
  rates <- exp(x)
  means <- rates %*% t(data$exposure)
  logp <- rowSums(stats::dpois(matrix(data$entered, nrow(x), 3, byrow=TRUE), means,
                               log=TRUE)) +
    stats::dnorm(x[, 1], walk$mu, walk$sigma, log=TRUE) +
    stats::dnorm(x[, 2] - x[, 1], 0, walk$tau, log=TRUE) +
    stats::dnorm(x[, 3] - x[, 2], 0, walk$tau, log=TRUE)
  weights <- exp(logp - max(logp))
  weights <- weights / sum(weights)
  set.seed(17)
  n_draws <- 6000
  f_logl <- function(log_rate) wrong_walk_logl(log_rate, data)
  logl <- f_logl(walk$log_rate)
  draws <- matrix(NA, n_draws, 3)
  for(i in seq_len(n_draws)) {
    out <- wrong_walk_slice(walk, f_logl, logl)
    walk <- out$walk
    logl <- out$logl
    draws[i, ] <- walk$log_rate
  }
  expect_equal(logl, f_logl(walk$log_rate))
  thinned <- seq(1, n_draws, 5)
  for(J in 1:3) {
    marginal <- tapply(weights, x[, J], sum)
    mean_exact <- sum(grid * marginal)
    var_exact <- sum(grid^2 * marginal) - mean_exact^2
    expect_mcmc_mean(draws[, J], mean_exact, paste("mean of log rate", J))
    expect_mcmc_mean((draws[, J] - mean_exact)^2, var_exact, paste("variance of log rate", J))
    cdf <- stats::approxfun(grid + diff(grid)[1] / 2, cumsum(marginal), yleft=0, yright=1)
    expect_ks(draws[thinned, J], cdf, paste("log rate", J))
  }
})

test_that("S7.4 (C4): with no wrong records the walk and its sd are drawn from their prior", {
  walk <- list(log_rate=log(c(0.1, 0.1, 0.1)), report_days=c(7, 14, 28), tau=0.3,
               mu=log(0.07), sigma=1.2, tau_scale=0.5, tau_step=0.5)
  data <- list(entered=numeric(0), exposure=matrix(0, 0, 3))
  expect_equal(wrong_walk_logl(walk$log_rate, data), 0)
  set.seed(18)
  n_draws <- 8000
  draws <- matrix(NA, n_draws, 4)
  for(i in seq_len(n_draws)) {
    walk <- wrong_walk_step(walk, data)
    draws[i, ] <- c(walk$log_rate, walk$tau)
  }
  expect_mcmc_mean(draws[, 1], log(0.07), "mean of the first log rate")
  expect_mcmc_mean((draws[, 1] - log(0.07))^2, 1.2^2, "variance of the first log rate")
  # the steps have variance tau^2 (days) / 7, and E tau^2 = 0.5^2
  expect_mcmc_mean((draws[, 3] - draws[, 1])^2, 0.5^2 * 21 / 7,
                   "mean square of the change over the walk")
  expect_mcmc_mean(draws[, 4], 0.5 * sqrt(2 / pi), "mean of tau")
  expect_mcmc_mean(draws[, 4]^2, 0.5^2, "mean square of tau")
  thinned <- seq(1, n_draws, 10)
  expect_ks(draws[thinned, 1], function(x) stats::pnorm(x, log(0.07), 1.2), "first log rate")
  expect_ks(draws[thinned, 4], function(x) 2 * stats::pnorm(x, 0, 0.5) - 1, "tau")
})

test_that("C6: simulate_cases_deaths enters wrong records at the rate of each interval", {
  set.seed(19)
  days <- c(2, 3, 6, 10)
  rates <- c(0.05, 0.4, 0.1, 0.2)
  rho <- 0.25
  delta <- 0.3
  n_sims <- 2000
  counts <- replicate(n_sims, {
    sim <- simulate_cases_deaths(c(30, 20), days, list(mean=3, sd=2), 0,
                                 list(mean=5, sd=3),
                                 wrong=list(eta=rates, rho=rho, delta=delta))
    rbind(sim$snapshots$cases_reported, sim$snapshots$deaths_reported)
  })
  f_expected <- function(t, cases) {
    d <- days[days >= t]
    p_reported <- function(x) stats::pgamma(x - t + 0.5, shape=(3 / 2)^2, rate=3 / 2^2)
    genuine <- cases * p_reported(d)
    wrong <- vapply(d, function(r) {
      if(r <= t)
        return(0)
      k <- (t + 1):r
      sum(rates[walk_interval_daily(k, days)] * cases * (p_reported(k) - p_reported(k - 1)) *
            (1 - rho)^(r - k))
    }, numeric(1))
    list(cases=genuine + wrong, deaths=delta * wrong)
  }
  expected <- list(cases=c(f_expected(1, 30)$cases, f_expected(2, 20)$cases),
                   deaths=c(f_expected(1, 30)$deaths, f_expected(2, 20)$deaths))
  zmax <- z_bonferroni(2 * length(expected$cases))
  for(i in seq_along(expected$cases)) {
    expect_z(mean(counts[1, i, ]), expected$cases[i], stats::sd(counts[1, i, ]) / sqrt(n_sims),
             paste("mean count of cases, row", i), zmax=zmax)
    if(expected$deaths[i] > 0)
      expect_z(mean(counts[2, i, ]), expected$deaths[i],
               stats::sd(counts[2, i, ]) / sqrt(n_sims),
               paste("mean count of deaths, row", i), zmax=zmax)
  }
  # the simulator draws a walk with steps of variance tau^2 (days) / 7
  steps <- t(replicate(3000, {
    sim <- simulate_cases_deaths(c(3, 3), days, list(mean=3, sd=2), 0, list(mean=5, sd=3),
                                 wrong=list(eta=0.1, tau=0.6, rho=rho, delta=delta))
    expect_equal(sim$wrong_rates$time_reported, days)
    diff(c(log(0.1), log(sim$wrong_rates$eta)))
  }))
  expect_true(all(steps[, 1] == 0))
  for(J in 2:4)
    expect_z(mean(steps[, J]^2), 0.6^2 * diff(days)[J - 1] / 7,
             0.6^2 * diff(days)[J - 1] / 7 * sqrt(2 / 3000),
             paste("mean square of step", J))
  expect_error(simulate_cases_deaths(c(3, 3), days, list(mean=3, sd=2), 0,
                                     list(mean=5, sd=3),
                                     wrong=list(eta=c(0.1, 0.2), rho=rho, delta=delta)),
               "rate of wrong records")
})

test_that("C8 (S7.5): predictions take the rate of the last interval after the last report", {
  example <- walk_example(function(report_days) walk_revision(
    log(seq(0.05, 0.3, length.out=length(report_days))), report_days), seed=20)
  state <- example$state
  state$uncorrected <- NULL
  revision <- example$revision
  prediction <- list(deaths_state=state, revision_parameters=revision,
                     reporting_by_onset=example$reporting, cases_true=example$cases_true,
                     max_erroneous=example$max_erroneous)
  k <- 6
  days_before <- state$reports[[k]]$time_reported
  days <- c(days_before, max(example$report_days) + c(2, 5))
  terms <- predict_onset_terms(prediction, k, days)
  n <- length(days)
  units_main <- wrong_entry_units(days, state$onset_times[k], list(mean=3, sd=2),
                                  list(eta=1, rho=0.3))
  rate_last <- exp(revision$wrong_walk$log_rate[length(example$report_days)])
  expect_equal(terms$wrong$lambda[n - 0:1],
               rate_last * example$cases_true[k] * units_main$unit[n - 0:1])
  # and with every rate equal, as the main model
  revision$wrong_walk$log_rate[] <- log(0.2)
  prediction$revision_parameters <- revision
  main <- prediction
  main$revision_parameters <- list(eta=0.2, rho=0.3)
  expect_equal(predict_onset_terms(prediction, k, days),
               predict_onset_terms(main, k, days), tolerance=1e-10)
})

# a short fit with a walk of the rate of wrong records
walk_fit <- function(example, niterations, ...) {
  snapshots <- example$snapshots
  mcmc_single(
    niterations=niterations,
    snapshot_with_Rt_index_df=snapshots %>%
      dplyr::select(time_onset, time_reported, cases_reported) %>%
      dplyr::mutate(Rt_index=1),
    priors=list(Rt=list(shape=1, rate=1),
                reporting=list(mean_mu=6, mean_sigma=5, sd_mu=4, sd_sigma=4),
                max_cases=300, overdispersion=list(mean=20, sd=20),
                revisions=list(wrong_rate="walk", rho_alpha=1, rho_beta=1),
                max_erroneous=8,
                Rt_changepoints=list(),
                reporting_changepoints=list(min_length=14)),
    serial_parameters=example$serial_parameters,
    initial_cases_true=example$initial_cases,
    initial_reporting_parameters=list(mean=6, sd=4),
    initial_Rt=dplyr::tibble(Rt_index=1, Rt=1),
    initial_revision_parameters=list(eta=0.05, rho=0.3, weekday=TRUE),
    revision_metropolis_parameters=list(eta_step=0.5, rho_step=0.5, tau_eta_step=0.4),
    print_to_screen=FALSE, is_negative_binomial=TRUE, initial_overdispersion=10,
    cases_before=example$cases_before,
    deaths=list(
      snapshot=snapshots %>%
        dplyr::select(time_onset, time_reported, deaths_reported),
      settings=list(onset_to_death=list(mean=9.3, sd=5.5), fatal_delay=FALSE,
                    triggered=list(initial=1, n_steps=0L),
                    dated_by_death=list(fraction=list(), exact=TRUE))),
    ...)
}

walk_fit_example <- function() {
  set.seed(21)
  n_days <- 20
  serial_parameters <- list(mean=5, sd=3)
  cases <- true_cases(n_days + 40, function(t) 1, kappa=1000, serial_parameters,
                      initial_parameters=list(mean=20, length=10))
  sim <- simulate_cases_deaths(cases[41:(n_days + 40)], seq(2, n_days + 20, 2),
                               list(mean=6, sd=4), 0.5, list(mean=3, sd=2.5),
                               triggered=1, onset_to_death=list(mean=9.3, sd=5.5),
                               dated_by_death=list(fraction=0.3),
                               wrong=list(eta=0.05, tau=0.5, rho=0.3, delta=0.3))
  snapshots <- dplyr::filter(sim$snapshots, time_reported <= n_days + 10)
  list(snapshots=snapshots, cases_before=cases[1:40],
       serial_parameters=serial_parameters,
       initial_cases=snapshots %>%
         dplyr::group_by(time_onset) %>%
         dplyr::summarise(cases_true=max(cases_reported) + 2))
}

test_that("C7 and C8: a fit with a walk runs, is taken up again and predicts later reports", {
  example <- walk_fit_example()
  set.seed(22)
  whole <- walk_fit(example, 6)
  revisions <- whole$revisions
  report_days <- sort(unique(example$snapshots$time_reported))
  expect_equal(colnames(revisions), c("eta", "rho", "pi", "iteration", "tau_eta",
                                      paste0("eta_", report_days)))
  expect_true(all(is.na(revisions$eta)))
  rates <- as.matrix(revisions[, paste0("eta_", report_days)])
  expect_true(all(rates > 0))
  expect_equal(nrow(unique(rates)), 6)
  expect_true(all(revisions$tau_eta > 0))
  expect_true(is.list(whole$wrong_cap))
  # C7: taken up again, the draws are those of the fit run straight through
  set.seed(22)
  first <- walk_fit(example, 2, return_state=TRUE)
  set.seed(98)
  second <- walk_fit(example, 4, resume_from=first$state)
  names_draws <- setdiff(names(whole), c("state", "iterations_done", "wrong_cap"))
  for(name in names_draws)
    expect_equal(dplyr::bind_rows(first[[name]], second[[name]]), whole[[name]],
                 info=name)
  # C8: a kept draw gives the walk, whose last rate the predictions use
  draw <- lapply(second[c("cases", "reporting_days", "revisions", "reporting_effects",
                          "deaths", "deaths_parameters")],
                 function(df) df[df$iteration == 6, ])
  variables <- second$state$variables
  prediction <- predict_state_from_draw(variables, draw)
  walk <- prediction$revision_parameters$wrong_walk
  expect_equal(walk$report_days, report_days)
  expect_equal(exp(walk$log_rate), as.numeric(rates[6, ]))
  expect_equal(walk$tau, revisions$tau_eta[6])
  expect_equal(walk, variables$revision_current$wrong_walk)
  set.seed(23)
  later <- predict_reports(variables, draw, c(34, 37), onset_times=c(5, 18))
  expect_equal(nrow(later), 4)
  expect_true(all(later$deaths_reported <= later$cases_reported))
})
