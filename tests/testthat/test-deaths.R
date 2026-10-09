test_that("deaths_settings gives defaults, takes changes and checks them", {
  settings <- deaths_settings(50)
  expect_equal(settings$beta, 50)
  expect_equal(settings$min_length, 14)
  expect_equal(settings$delay_changepoints$min_length, 14)
  expect_equal(deaths_settings(50, list(min_length=7))$min_length, 7)
  expect_equal(deaths_settings(50, list(delay_changepoints=FALSE))$
                 delay_changepoints$n_partition_moves, 0)
  expect_error(deaths_settings(50, list(unknown=1)))
  expect_error(deaths_settings(50, list(cfr_alpha=0)))
  expect_error(deaths_settings(50, list(initial_cfr=1)))
  expect_null(settings$dated_by_death)
  expect_true("dated_by_death" %in% names(settings))
})

test_that("dated_by_death_settings gives the times from onset to death", {
  expect_null(dated_by_death_settings(10, NULL))
  expect_null(dated_by_death_settings(10, list(fraction=0, mean=8, sd=4)))
  dated <- deaths_settings(3, list(dated_by_death=list(
    fraction=c(0.1, 0.2, 0.3), mean=8, sd=4)))$dated_by_death
  expect_equal(dated$fraction, c(0.1, 0.2, 0.3))
  expect_equal(dated$max_lag, ceiling(stats::qgamma(0.999, shape=4, rate=0.5)))
  expect_equal(sum(dated$lag_probabilities), 1)
  # a time from onset to death is rounded up, so is at least one day
  expect_equal(dated$lag_probabilities[2] / dated$lag_probabilities[1],
               (stats::pgamma(2, 4, 0.5) - stats::pgamma(1, 4, 0.5)) /
                 stats::pgamma(1, 4, 0.5))
  expect_null(dated$reporting)
  expect_equal(dated$n_moves, 12)
  with_delay <- dated_by_death_settings(5, list(
    fraction=0.2, mean=8, sd=4, max_lag=6, reporting=list(mean=2, sd=1), n_moves=3))
  expect_equal(length(with_delay$lag_probabilities), 6)
  expect_equal(with_delay$reporting, list(mean=2, sd=1))
  expect_equal(with_delay$n_moves, 3)
  expect_error(dated_by_death_settings(5, list(fraction=c(0.1, 0.2), mean=8, sd=4)))
  expect_error(dated_by_death_settings(5, list(fraction=1.2, mean=8, sd=4)))
  expect_error(dated_by_death_settings(5, list(fraction=0.2)))
  expect_error(dated_by_death_settings(5, list(fraction=0.2, mean=8, sd=4, other=1)))
  expect_error(dated_by_death_settings(5, list(fraction=0.2, mean=8, sd=4,
                                               reporting=list(mean=2))))
})

test_that("segment_log_marginals_cfr matches numerical integration", {
  deaths <- c(3, 5, 2, 8)
  cases <- c(10, 12, 9, 11)
  log_marginals <- segment_log_marginals_cfr(deaths, cases, 2, 3)
  f_integrand <- function(cfr, days)
    sapply(cfr, function(x)
      prod(x^deaths[days] * (1 - x)^(cases[days] - deaths[days])) *
        stats::dbeta(x, 2, 3))
  for(s in 1:4)
    for(t in s:4)
      expect_equal(log_marginals[s, t],
                   log(stats::integrate(f_integrand, 0, 1, days=s:t,
                                        rel.tol=1e-10)$value),
                   tolerance=1e-6)
  expect_true(all(log_marginals[lower.tri(log_marginals)] == -Inf))
  expect_error(segment_log_marginals_cfr(c(5, 1), c(4, 2)))
  expect_error(segment_log_marginals_cfr(1, c(4, 2)))
})

test_that("sample_cfr_changepoints finds a large change and recovers the ratios", {
  set.seed(1)
  cases <- rep(60, 40)
  deaths <- c(stats::rbinom(20, 60, 0.6), stats::rbinom(20, 60, 0.2))
  sampled <- sample_cfr_changepoints(deaths, cases, 0.05, min_length=5)
  expect_equal(sampled$starts, c(1, 21))
  expect_equal(length(sampled$cfr), 40)
  expect_true(abs(sampled$cfr[1] - 0.6) < 0.06)
  expect_true(abs(sampled$cfr[40] - 0.2) < 0.06)
  constant <- sample_cfr_changepoints(deaths, cases, 0.05, is_changepoints=FALSE)
  expect_equal(constant$starts, 1)
  expect_equal(length(unique(constant$cfr)), 1)
  expect_error(sample_cfr_changepoints(deaths, cases, 0))
})

test_that("joint_probabilities are those of the two delays", {
  days <- c(3, 5, 9)
  probabilities <- joint_probabilities(days, 2, list(mean=5, sd=3),
                                       list(mean=10, sd=6))
  s_reported <- exp(log_undetected_days(days, 2, list(mean=5, sd=3)))
  s_recorded <- exp(log_undetected_days(days, 2, list(mean=10, sd=6)))
  expect_equal(probabilities$p, 1 - s_reported / c(1, s_reported[-3]))
  expect_equal(probabilities$g, 1 - s_recorded)
  expect_equal(probabilities$q, 1 - s_recorded / c(1, s_recorded[-3]))
  # fatal cases dated by their death are shown as dead once reported
  dated <- joint_probabilities(days, 2, list(mean=5, sd=3), list(mean=10, sd=6),
                               fraction_dated_by_death=0.3)
  expect_equal(dated$g, 0.3 + 0.7 * (1 - s_recorded))
  expect_equal(dated$p, probabilities$p)
  expect_equal(dated$p_fatal, probabilities$p)
  expect_equal(dated$q, probabilities$q)
})

test_that("mixed_probabilities are those of a mixture of the two kinds of fatal case", {
  days <- c(3, 5, 9, 14)
  survivals <- joint_survivals(days, 2, list(mean=6, sd=3), list(mean=10, sd=6),
                               dated_reporting_parameters=list(mean=2, sd=1))
  s_onset <- exp(survivals$log_reported)
  s_dated <- exp(survivals$log_dated)
  s_recorded <- exp(survivals$log_recorded)
  expect_false(survivals$is_dated_same)
  expect_equal(s_dated, exp(log_undetected_days(days, 2, list(mean=2, sd=1))))
  fractions <- c(0, 0.3, 1)
  probabilities <- mixed_probabilities(survivals, fractions)
  expect_equal(dim(probabilities$p_fatal), c(3, 4))
  for(i in 1:3) {
    fraction <- fractions[i]
    mixture <- (1 - fraction) * s_onset + fraction * s_dated
    mixture_before <- c(1, mixture[-4])
    expect_equal(probabilities$p_fatal[i, ], 1 - mixture / mixture_before)
    # the probability a fatal case newly reported is dated by death
    share <- fraction * (c(1, s_dated[-4]) - s_dated) / (mixture_before - mixture)
    expect_equal(probabilities$g[i, ], share + (1 - share) * (1 - s_recorded))
  }
  expect_equal(probabilities$p, 1 - s_onset / c(1, s_onset[-4]))
  expect_equal(probabilities$p_fatal[1, ], probabilities$p)
  # cases dated by death are reported sooner, so are more of those reported early
  share <- (probabilities$g[2, ] - (1 - s_recorded)) / s_recorded
  expect_true(share[1] > 0.3 & share[4] < 0.3)
  single <- path_probabilities(survivals, 0.3)
  expect_equal(single$p_fatal, probabilities$p_fatal[2, ])
  expect_equal(single$g, probabilities$g[2, ])
})

test_that("deaths_new_logp is the convolution of two binomials", {
  f_exact <- function(deaths_new, fatal_new, fatal_alive, g, q)
    log(sum(sapply(0:deaths_new, function(x)
      stats::dbinom(x, fatal_new, g) * stats::dbinom(deaths_new - x, fatal_alive, q))))
  expect_equal(deaths_new_logp(2, 3, 4, 0.3, 0.2), f_exact(2, 3, 4, 0.3, 0.2))
  expect_equal(deaths_new_logp(0, 0, 0, 0.3, 0.2), 0)
  expect_equal(deaths_new_logp(1, 0, 0, 0.3, 0.2), -Inf)
  expect_equal(deaths_new_logp(2, c(0, 3), c(4, 1), 0.3, 0.2),
               c(f_exact(2, 0, 4, 0.3, 0.2), f_exact(2, 3, 1, 0.3, 0.2)))
})

# the probability of every possible split of the reported cases, by enumeration
enumerate_fatal_paths <- function(cases, deaths, deaths_true, survivors_true,
                                  probabilities) {
  n_reports <- length(cases)
  p_fatal <- probabilities$p_fatal
  if(is.null(p_fatal))
    p_fatal <- probabilities$p
  grid <- expand.grid(lapply(cases, function(x) 0:x))
  weights <- apply(grid, 1, function(fatal) {
    fatal <- c(0, fatal)
    cases_all <- c(0, cases)
    deaths_all <- c(0, deaths)
    weight <- 1
    for(j in seq_len(n_reports)) {
      a <- fatal[j + 1] - fatal[j]
      b <- cases_all[j + 1] - cases_all[j] - a
      if(a < 0 | b < 0 | fatal[j + 1] < deaths_all[j + 1] | fatal[j + 1] > deaths_true |
         cases_all[j + 1] - fatal[j + 1] > survivors_true | fatal[j] < deaths_all[j])
        return(0)
      deaths_new <- deaths_all[j + 1] - deaths_all[j]
      weight <- weight *
        stats::dbinom(a, deaths_true - fatal[j], p_fatal[j]) *
        stats::dbinom(b, survivors_true - (cases_all[j] - fatal[j]), probabilities$p[j]) *
        sum(stats::dbinom(0:deaths_new, a, probabilities$g[j]) *
              stats::dbinom(deaths_new:0, fatal[j] - deaths_all[j], probabilities$q[j]))
    }
    weight
  })
  names(weights) <- apply(grid, 1, paste, collapse="-")
  weights
}

test_that("fatal_path_forward gives the probability of the reports", {
  cases <- c(2, 3, 5)
  deaths <- c(0, 1, 2)
  probabilities <- list(p=c(0.3, 0.4, 0.5), g=c(0.2, 0.4, 0.6), q=c(0.2, 0.25, 0.33))
  for(deaths_true in 2:4) {
    weights <- enumerate_fatal_paths(cases, deaths, deaths_true, 3, probabilities)
    forward <- fatal_path_forward(cases, deaths, deaths_true, 3, probabilities)
    expect_equal(forward$logp, log(sum(weights)))
  }
  # fatal cases reported at a different speed from others
  probabilities$p_fatal <- c(0.5, 0.2, 0.7)
  weights <- enumerate_fatal_paths(cases, deaths, 3, 3, probabilities)
  forward <- fatal_path_forward(cases, deaths, 3, 3, probabilities)
  expect_equal(forward$logp, log(sum(weights)))
  probabilities$p_fatal <- NULL
  # impossible: more deaths shown than fatal cases
  expect_equal(fatal_path_forward(cases, deaths, 1, 5, probabilities)$logp, -Inf)
  expect_error(sample_fatal_path(cases, deaths, 1, 5, probabilities))
  expect_error(fatal_path_forward(c(3, 2), c(0, 0), 2, 2, probabilities))
  expect_error(fatal_path_forward(c(1, 2), c(2, 2), 2, 2, probabilities))
})

test_that("sample_fatal_path draws splits with the right probabilities", {
  set.seed(2)
  cases <- c(2, 3, 5)
  deaths <- c(0, 1, 2)
  probabilities <- list(p=c(0.3, 0.4, 0.5), g=c(0.2, 0.4, 0.6), q=c(0.2, 0.25, 0.33))
  weights <- enumerate_fatal_paths(cases, deaths, 3, 3, probabilities)
  weights <- weights[weights > 0] / sum(weights)
  n_draws <- 5000
  draws <- replicate(n_draws, {
    path <- sample_fatal_path(cases, deaths, 3, 3, probabilities)
    expect_equal(path$fatal_new + path$survivors_new, diff(c(0, cases)))
    expect_true(all(path$alive_before >= 0))
    paste(path$fatal_before + path$fatal_new, collapse="-")
  })
  observed <- table(factor(draws, levels=names(weights)))
  expect_equal(sum(observed), n_draws)
  expect_true(stats::chisq.test(observed, p=weights)$p.value > 0.001)
})

test_that("joint_cases_deaths_logp sums to the probability of the split of the reported cases", {
  cases <- c(2, 3, 5)
  deaths <- c(0, 1, 2)
  probabilities <- list(p=c(0.3, 0.4, 0.5), g=c(0.2, 0.4, 0.6), q=c(0.2, 0.25, 0.33))
  set.seed(3)
  path <- sample_fatal_path(cases, deaths, 3, 3, probabilities)
  survivals <- list(log_reported=cumsum(log(1 - probabilities$p)),
                    log_recorded=log(1 - probabilities$g), is_dated_same=TRUE)
  survivals$log_dated <- survivals$log_reported
  logp <- joint_cases_deaths_logp(path, survivals, 0.4, 5:9)
  fatal_values <- attr(logp, "deaths_true")
  expect_equal(nrow(logp), 5)
  expect_equal(ncol(logp), length(fatal_values))
  fatal_reported <- sum(path$fatal_new)
  expect_equal(fatal_values[1], fatal_reported)
  # an element against the definition
  f_exact <- function(cases_true, deaths_true)
    stats::dbinom(deaths_true, cases_true, 0.4, log=TRUE) +
    sum(stats::dbinom(path$fatal_new, deaths_true - path$fatal_before,
                      probabilities$p, log=TRUE)) +
    sum(stats::dbinom(path$survivors_new,
                      cases_true - deaths_true - path$survivors_before,
                      probabilities$p, log=TRUE))
  expect_equal(logp[3, 2], f_exact(7, fatal_values[2]))
  # fewer survivors than those reported is impossible
  expect_equal(logp[1, ncol(logp)], -Inf)

  # with fatal cases dated by death: two of the fatal cases of the onset time
  # are recorded for later days, and one from an earlier onset time for this
  survivals$log_dated <- cumsum(log(1 - c(0.6, 0.5, 0.5)))
  survivals$is_dated_same <- FALSE
  logp <- joint_cases_deaths_logp(path, survivals, 0.4, 5:9, dated_out=2,
                                  dated_in=1, fraction=0.3)
  fatal_values <- attr(logp, "deaths_true")
  expect_equal(fatal_values[1], fatal_reported - 1 + 2)
  f_exact <- function(cases_true, deaths_true) {
    fatal_recorded <- deaths_true - 2 + 1
    probabilities_mixed <- path_probabilities(survivals, 1 / fatal_recorded)
    stats::dbinom(deaths_true, cases_true, 0.4, log=TRUE) +
      stats::dbinom(2, deaths_true, 0.3, log=TRUE) +
      sum(stats::dbinom(path$fatal_new, fatal_recorded - path$fatal_before,
                        probabilities_mixed$p_fatal, log=TRUE)) +
      sum(deaths_new_logp(path$deaths_new, path$fatal_new, path$alive_before,
                          probabilities_mixed$g, probabilities_mixed$q)) +
      sum(stats::dbinom(path$survivors_new,
                        cases_true - deaths_true - path$survivors_before,
                        probabilities$p, log=TRUE))
  }
  expect_equal(logp[4, 2], f_exact(8, fatal_values[2]))
  expect_equal(logp[5, 1], f_exact(9, fatal_values[1]))
  expect_equal(
    fatal_reports_logp(path, survivals, c(4, 0), dated_in=1)[2], -Inf)
})

test_that("simulate_cases_deaths gives reports consistent with its truth", {
  set.seed(4)
  cases <- rep(200, 30)
  sim <- simulate_cases_deaths(cases, 1:80, list(mean=5, sd=3), 0.3,
                               list(mean=10, sd=5))
  expect_true(all(sim$snapshots$deaths_reported <= sim$snapshots$cases_reported))
  latest <- sim$snapshots %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(cases=max(cases_reported), deaths=max(deaths_reported))
  expect_equal(latest$cases, sim$truth$cases_true)
  expect_equal(latest$deaths, sim$truth$deaths_true)
  expect_true(abs(sum(sim$truth$deaths_true) / sum(cases) - 0.3) < 0.03)
  falls <- sim$snapshots %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(fall=any(diff(deaths_reported) < 0) | any(diff(cases_reported) < 0))
  expect_false(any(falls$fall))
  # the ratio can differ between onset times
  sim <- simulate_cases_deaths(cases, 1:80, list(mean=5, sd=3),
                               c(rep(0.6, 15), rep(0.1, 15)), list(mean=10, sd=5))
  expect_true(mean(sim$truth$deaths_true[1:15]) > 3 * mean(sim$truth$deaths_true[16:30]))
  # a death recorded only after the case is reported comes later
  sim_report <- simulate_cases_deaths(cases, 1:80, list(mean=5, sd=3), 0.3,
                                      list(mean=10, sd=5), death_from="report")
  expect_true(mean(sim_report$truth$delay_dead_mean) > mean(sim$truth$delay_dead_mean))
  expect_error(simulate_cases_deaths(cases, 1:80, list(mean=5, sd=3), c(0.1, 0.2),
                                     list(mean=10, sd=5)))
  # fatal cases dated by their death appear later and are dead once reported
  sim_dated <- simulate_cases_deaths(
    cases, 1:120, list(mean=5, sd=3), 0.3, list(mean=10, sd=5),
    dated_by_death=list(fraction=1, mean=8, sd=4))
  expect_equal(sum(sim_dated$truth$n_dated_by_death), sum(sim_dated$truth$deaths_true))
  expect_true(sum(sim_dated$truth$cases_true[1:3]) < sum(cases[1:3]))
  expect_true(all(sim_dated$truth$cases_onset == cases))
  half <- simulate_cases_deaths(
    cases, 1:120, list(mean=5, sd=3), 0.3, list(mean=10, sd=5),
    dated_by_death=list(fraction=0.5, mean=8, sd=4))
  expect_true(abs(sum(half$truth$n_dated_by_death) /
                    sum(half$truth$deaths_true) - 0.5) < 0.1)
  expect_true(all(half$snapshots$deaths_reported <= half$snapshots$cases_reported))
  expect_true(all(half$truth$deaths_onset <= half$truth$cases_onset))
  # a case dated by death is recorded for a day after its onset
  expect_equal(sim_dated$truth$deaths_true[1], 0)
  # with their own short delay these cases are reported sooner
  f_reported_early <- function(reporting) {
    sim <- simulate_cases_deaths(
      cases, 1:120, list(mean=20, sd=5), 1, list(mean=10, sd=5),
      dated_by_death=list(fraction=1, mean=8, sd=4, reporting=reporting))
    same_day <- dplyr::filter(sim$snapshots, time_reported == time_onset + 2)
    sum(same_day$deaths_reported) / sum(sim$truth$deaths_true)
  }
  expect_true(f_reported_early(list(mean=1, sd=0.5)) > 0.9)
  expect_true(f_reported_early(NULL) < 0.01)
})

test_that("mcmc_single with deaths recovers a case fatality ratio which changes", {
  set.seed(5)
  n_days <- 40
  serial_parameters <- list(mean=5, sd=3)
  cases <- true_cases(n_days + 40, function(t) 1, kappa=1000, serial_parameters,
                      initial_parameters=list(mean=30, length=10))
  cases_before <- cases[1:40]
  cases <- cases[41:(n_days + 40)]
  cfr <- c(rep(0.6, 20), rep(0.2, 20))
  sim <- simulate_cases_deaths(cases, seq_len(n_days + 30), list(mean=4, sd=3), cfr,
                               list(mean=6, sd=4))
  snapshots <- sim$snapshots
  initial_cases <- snapshots %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(cases_true=max(cases_reported))
  priors <- list(Rt=list(shape=1, rate=1),
                 reporting=list(mean_mu=5, mean_sigma=5, sd_mu=3, sd_sigma=3),
                 max_cases=300,
                 Rt_changepoints=list(),
                 reporting_changepoints=list(min_length=14))
  f_fit <- function(...) mcmc_single(
    niterations=30,
    snapshot_with_Rt_index_df=snapshots %>%
      dplyr::select(time_onset, time_reported, cases_reported) %>%
      dplyr::mutate(Rt_index=1),
    priors=priors,
    serial_parameters=serial_parameters,
    initial_cases_true=initial_cases,
    initial_reporting_parameters=list(mean=4, sd=3),
    initial_Rt=dplyr::tibble(Rt_index=1, Rt=1),
    print_to_screen=FALSE,
    cases_before=cases_before, ...)
  deaths <- list(snapshot=snapshots %>%
                   dplyr::select(time_onset, time_reported, deaths_reported),
                 settings=list(min_length=7,
                               delay_prior=list(mean_mu=6, mean_sigma=5,
                                                sd_mu=4, sd_sigma=4),
                               initial_delay=list(mean=8, sd=5)))
  progress_file <- tempfile()
  fit <- f_fit(deaths=deaths, progress_file=progress_file, progress_every=10)
  expect_true(grepl("^30 of 30 iterations", readLines(progress_file)))
  expect_error(f_fit(deaths=deaths, progress_file=progress_file, progress_every=0))
  expect_true(all(c("deaths", "deaths_parameters") %in% names(fit)))
  expect_equal(nrow(fit$deaths), 30 * n_days)
  expect_equal(nrow(fit$deaths_parameters), 30)
  kept <- fit$deaths %>% dplyr::filter(iteration > 15)
  estimates <- kept %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(cfr=stats::median(cfr), deaths=stats::median(deaths_true))
  expect_true(abs(estimates$cfr[5] - 0.6) < 0.1)
  expect_true(abs(estimates$cfr[35] - 0.2) < 0.1)
  # all deaths of early onset times have been shown
  expect_equal(estimates$deaths[1:15], sim$truth$deaths_true[1:15])
  # fatal cases never exceed cases
  both <- dplyr::left_join(fit$deaths, fit$cases, by=c("iteration", "time_onset"))
  expect_true(all(both$deaths_true <= both$cases_true))
  expect_true(all(fit$deaths_parameters$n_cfr_changes[16:30] >= 1))

  # the requirements of the model of deaths
  expect_error(f_fit(deaths=list()))
  priors_no_changes <- priors
  priors_no_changes$reporting_changepoints <- NULL
  expect_error(mcmc_single(
    niterations=2,
    snapshot_with_Rt_index_df=snapshots %>%
      dplyr::select(time_onset, time_reported, cases_reported) %>%
      dplyr::mutate(Rt_index=1),
    priors=priors_no_changes, serial_parameters=serial_parameters,
    initial_cases_true=initial_cases,
    initial_reporting_parameters=list(mean=4, sd=3),
    initial_Rt=dplyr::tibble(Rt_index=1, Rt=1), print_to_screen=FALSE,
    deaths=deaths))
})
