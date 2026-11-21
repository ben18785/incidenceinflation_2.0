# S8 of the sensitivity analyses: the share of cases in the heavier tail of
# the reporting delay differs by onset week, following a random walk on the
# logit scale (R/deaths-tail-walk.R). The tests are labelled as in the test
# plan of sensitivity_extensions.tex.

# the survival of the reporting delay with a tail, written out directly
tail_direct_survival <- function(elapsed, mean, sd, weight, tail_mean)
  (1 - weight) * stats::pgamma(elapsed, shape=mean^2 / sd^2, rate=mean / sd^2,
                               lower.tail=FALSE) +
  weight * exp(-elapsed / tail_mean)

# a fit of the reference data (helper-reference.R) with the settings of the
# Bundibugyo fits and, if walk is not FALSE, a share of the tail for each
# onset week
tail_walk_fit <- function(data, niterations, seed=1, walk=TRUE, weight=0.05,
                          resume_from=NULL, is_uncorrected=TRUE, deaths=TRUE) {
  set.seed(seed)
  mcmc_single(
    niterations=niterations,
    snapshot_with_Rt_index_df=data$snapshots %>%
      dplyr::select(time_onset, time_reported, cases_reported) %>%
      dplyr::mutate(Rt_index=1),
    priors=list(Rt=list(shape=1, rate=1),
                reporting=list(mean_mu=8, mean_sigma=5, sd_mu=6, sd_sigma=4),
                max_cases=400, overdispersion=list(mean=20, sd=20),
                revisions=list(eta_mu=0.1, eta_sigma=0.1, rho_alpha=1, rho_beta=1),
                max_erroneous=8,
                reporting_effects=list(tail_weight_alpha=1, tail_weight_beta=19,
                                       tail_mean_mu=40, tail_mean_sigma=30,
                                       tail_weight_walk=walk),
                Rt_changepoints=list(),
                reporting_changepoints=list(min_length=14, n_partition_moves=5)),
    serial_parameters=data$serial_parameters,
    initial_cases_true=data$initial_cases,
    initial_reporting_parameters=list(mean=6, sd=4),
    initial_Rt=dplyr::tibble(Rt_index=1, Rt=1),
    reporting_metropolis_parameters=list(mean_step=0.3, sd_step=0.25),
    print_to_screen=FALSE, is_negative_binomial=TRUE, initial_overdispersion=20,
    overdispersion_metropolis_sd=8,
    initial_revision_parameters=list(eta=0.02, rho=0.2, weekday=TRUE,
                                     weekday_first=1,
                                     tail=list(weight=weight, mean=40)),
    revision_metropolis_parameters=list(eta_step=0.3, rho_step=0.3, pi_step=0.3),
    cases_before=data$cases_before, n_parameter_steps=2, return_state=TRUE,
    resume_from=resume_from,
    deaths=if(deaths) list(
      snapshot=data$snapshots %>%
        dplyr::select(time_onset, time_reported, deaths_reported),
      settings=list(
        onset_to_death=list(mean=9.3, sd=5.5),
        dated_by_death=list(fraction=list(alpha=1, beta=1, initial=0.25),
                            exact=TRUE),
        initial_cfr=0.4, initial_delay=list(mean=15, sd=13),
        delay_prior=list(mean_mu=4, mean_sigma=4, sd_mu=4, sd_sigma=4),
        delay_changepoints=list(min_length=14, n_partition_moves=5,
                                first_length=20),
        triggered=list(initial=1, n_steps=0L), fatal_delay=FALSE,
        recording_backlog=list(hold_days=5L)),
      uncorrected=if(is_uncorrected) data$uncorrected))
}

# the revision parameters with the tail of another onset week than that of
# day_onset
tail_other_week <- function(single, tail, day_onset) {
  week <- tail_week_of_days(tail, day_onset)
  others <- setdiff(seq_along(tail$weight), week)
  differing <- others[tail$weight[others] != tail$weight[week]]
  single$tail$weight <- tail$weight[differing[1]]
  single
}

# the shares of the tail of each onset week, made to differ
tail_with_weeks <- function(tail, weight) {
  tail$weight <- rep(weight, length.out=length(tail$weight))
  tail
}

# the log-probability of everything given the split of each onset time's
# reports which involves the reporting delay, computed afresh: every report's
# cases and deaths newly shown, whether or not they depend on the delay, and
# the wrong records entered and removed, for the onset times and the
# uncorrected series
tail_full_logl <- function(v, revision) {
  f_series <- function(state, cases_true, reporting_by_onset) {
    delays <- deaths_delay_by_onset(state)
    dated <- dated_counts(state)
    sum(vapply(seq_along(state$onset_times), function(k) {
      path <- state$paths[[k]]
      days <- state$reports[[k]]$time_reported
      reporting <- delay_of_onset(reporting_by_onset, k)
      survivals <- joint_survivals(days, state$onset_times[k], reporting,
                                   delay_of_onset(delays, k), revision,
                                   state$dated$reporting, state$recording_effects,
                                   state$nu, state$triggered)
      units <- wrong_entry_units(days, state$onset_times[k], reporting, revision)
      fatal_reports_logp(path, survivals, dated$fatal_recorded[k], dated$into[k],
                         is_deaths_included=TRUE) +
        sum(dbinom_log_recycled(path$survivors_new,
                                cases_true[k] - state$deaths_true[k] -
                                  path$survivors_before,
                                survival_steps(survivals$log_reported))) +
        sum(stats::dpois(path$wrong_entered_alive + path$wrong_entered_dead,
                         revision$eta * (cases_true[k] - dated$out[k] + dated$into[k]) *
                           units$unit, log=TRUE)) +
        sum(dbinom_log_recycled(path$wrong_removed, path$wrong_before, units$w))
    }, numeric(1)))
  }
  state <- v$deaths_state
  out <- f_series(state, v$deaths_cases_true, v$deaths_running_df)
  sub <- uncorrected_state(state)
  if(!is.null(sub))
    out <- out + f_series(sub, sub$cases_true,
                          uncorrected_reporting(state, v$deaths_running_df))
  out
}

# the log-prior of the walk, written out
tail_walk_prior_direct <- function(tail, walk) {
  x <- stats::qlogis(tail$weight)
  stats::dnorm(x[1], walk$first_mean, walk$first_sd, log=TRUE) +
    sum(stats::dnorm(diff(x), 0, tail$rw_sd, log=TRUE))
}

test_that("C1: the settings of the walk of the share of the tail are checked", {
  walk <- tail_walk_settings(TRUE)
  expect_equal(walk$first_mean, stats::qlogis(0.036))
  expect_equal(walk$first_sd, 1.2)
  expect_equal(walk$step_sd_scale, 0.15)
  expect_equal(tail_walk_settings(list(first_sd=2))$first_sd, 2)
  expect_error(tail_walk_settings(list(scale=1)), "Unknown settings for tail_weight_walk: scale")
  expect_error(tail_walk_settings(list(step_sd_scale=-0.1)), "must be positive")
  expect_error(tail_walk_settings(list(first_sd=0)), "must be positive")
  expect_error(tail_walk_settings(list(rw_sd_step=-1)), "zero or more")
  expect_error(tail_walk_settings(list(first_mean=c(1, 2))), "single number")
  expect_error(tail_walk_settings(0.5), "TRUE, FALSE or a named list")
  # the starting shares: one, or one for each onset week
  tail <- initial_tail_walk(list(weight=0.05, mean=40), 1:20, walk)
  expect_equal(tail$weight, rep(0.05, 3))
  expect_equal(tail$day_first, 1)
  expect_equal(tail$rw_sd, walk$initial_rw_sd)
  expect_equal(initial_tail_walk(list(weight=c(0.1, 0.2, 0.3), mean=40), 5:25,
                                 walk)$weight, c(0.1, 0.2, 0.3))
  expect_error(initial_tail_walk(list(weight=c(0.1, 0.2), mean=40), 1:20, walk),
               "one for each of the 3 onset weeks")
  expect_error(initial_tail_walk(list(weight=c(0.1, 1, 0.2), mean=40), 1:20, walk),
               "strictly between 0 and 1")
  expect_error(initial_tail_walk(list(weight=0.1, mean=-1), 1:20, walk),
               "positive number")

  # the walk needs the tail and a model of deaths
  data <- reference_data(101)
  initial <- list(eta=0.02, rho=0.2, tail=list(weight=0.05, mean=40))
  f_fit <- function(revision, deaths=TRUE)
    mcmc_single(
      niterations=1,
      snapshot_with_Rt_index_df=data$snapshots %>%
        dplyr::select(time_onset, time_reported, cases_reported) %>%
        dplyr::mutate(Rt_index=1),
      priors=list(Rt=list(shape=1, rate=1),
                  reporting=list(mean_mu=8, mean_sigma=5, sd_mu=6, sd_sigma=4),
                  max_cases=400, max_erroneous=8,
                  revisions=list(eta_mu=0.1, eta_sigma=0.1, rho_alpha=1, rho_beta=1),
                  reporting_effects=list(tail_weight_walk=TRUE),
                  reporting_changepoints=list()),
      serial_parameters=data$serial_parameters,
      initial_cases_true=data$initial_cases,
      initial_reporting_parameters=list(mean=6, sd=4),
      initial_Rt=dplyr::tibble(Rt_index=1, Rt=1), print_to_screen=FALSE,
      initial_revision_parameters=revision,
      deaths=if(deaths) list(snapshot=data$snapshots %>%
                               dplyr::select(time_onset, time_reported, deaths_reported)))
  expect_error(f_fit(initial[c("eta", "rho")]), "needs the tail")
  expect_error(f_fit(initial, deaths=FALSE), "only with a model of deaths")
  expect_error(f_fit(c(initial[c("eta", "rho")],
                       list(tail=list(weight=c(0.1, 0.2), mean=40)))),
               "one for each of the 5 onset weeks")
})

test_that("S8.1: each onset day reads the share of its week", {
  reporting <- list(mean=6, sd=4)
  tail <- list(weight=c(0.02, 0.2, 0.05), mean=30, day_first=3, rw_sd=0.1)
  # onset days 3 to 9 are week 1, 10 to 16 week 2, 17 on week 3; days before
  # the first and after the last take the shares of those weeks
  expect_equal(tail_week_of_days(tail, c(1, 3, 9, 10, 16, 17, 40)),
               c(1, 1, 1, 2, 2, 3, 3))
  expect_equal(tail_of_onset(tail, 12), list(weight=0.2, mean=30))
  expect_identical(tail_of_onset(list(weight=0.1, mean=30), 12),
                   list(weight=0.1, mean=30))
  for(day_onset in c(5, 12, 40)) {
    q <- tail$weight[tail_week_of_days(tail, day_onset)]
    days <- day_onset + 0:30
    weekly <- log_undetected_days(days, day_onset, reporting, list(tail=tail))
    expect_equal(weekly, log_undetected_days(days, day_onset, reporting,
                                             list(tail=list(weight=q, mean=30))),
                 tolerance=1e-12)
    expect_equal(exp(weekly), tail_direct_survival(days - day_onset + 0.5, 6, 4, q, 30),
                 tolerance=1e-10)
  }
  # two onset days of different weeks differ, as their shares do
  expect_false(isTRUE(all.equal(log_undetected_days(30, 5, reporting, list(tail=tail)),
                                log_undetected_days(37, 12, reporting, list(tail=tail)))))

  # with multipliers for each day, multipliers of the tail alone and a cache
  revision <- update_report_multipliers(list(
    report_days=c(6, 9, 14, 20, 30), weekday_first=2,
    weekday_log_effects=c("1"=0.3, "3"=-0.2, "5"=0.1, "6"=0, "7"=-0.1),
    report_log_effects=c(0.2, -0.3, 0.1, 0, 0.4), report_tail_only=TRUE))
  for(cache in list(NULL, new_lag_cache()))
    for(day_onset in c(5, 12)) {
      q <- tail$weight[tail_week_of_days(tail, day_onset)]
      days <- day_onset + 0:20
      weekly <- revision
      weekly$tail <- tail
      weekly$lag_cache <- cache
      single <- revision
      single$tail <- list(weight=q, mean=30)
      expect_equal(log_undetected_days(days, day_onset, reporting, weekly),
                   log_undetected_days(days, day_onset, reporting, single),
                   tolerance=1e-12)
    }

  # with every share equal it is unchanged
  equal <- tail_with_weeks(tail, 0.05)
  for(day_onset in c(5, 12, 40))
    expect_identical(log_undetected_days(day_onset + 0:30, day_onset, reporting,
                                         list(tail=equal)),
                     log_undetected_days(day_onset + 0:30, day_onset, reporting,
                                         list(tail=list(weight=0.05, mean=30))))
})

test_that("S8.2: every use of the tail takes the share of the onset week", {
  reporting <- list(mean=6, sd=4)
  tail <- list(weight=c(0.02, 0.2, 0.05), mean=30, day_first=1, rw_sd=0.1)
  revision <- list(eta=0.1, rho=0.3, tail=tail)
  days <- c(6, 9, 14, 20, 26, 35)
  for(day_onset in c(3, 9, 16)) {
    q <- tail$weight[tail_week_of_days(tail, day_onset)]
    days_onset <- days[days >= day_onset]
    f_survival <- function(d) tail_direct_survival(d - day_onset + 0.5, 6, 4, q, 30)
    # the entry of wrong records, day by day
    units <- wrong_entry_units(days_onset, day_onset, reporting, revision)
    days_before <- c(day_onset, days_onset[-length(days_onset)])
    unit_direct <- vapply(seq_along(days_onset), function(j) {
      k <- seq_len(days_onset[j] - days_before[j]) + days_before[j]
      if(length(k) == 0)
        return(0)
      sum((f_survival(k - 1) - f_survival(k)) * 0.7^(days_onset[j] - k))
    }, numeric(1))
    expect_equal(units$unit, unit_direct, tolerance=1e-10)
    # the survivals of the cases, fatal or not (with a slower delay for the
    # fatal), from which the binomials of the cases newly reported and the
    # probabilities of being shown as dead are made
    survivals <- joint_survivals(days_onset, day_onset, reporting,
                                 list(mean=12, sd=6), revision, nu=log(1.5))
    expect_equal(exp(survivals$log_reported), f_survival(days_onset), tolerance=1e-10)
    expect_equal(exp(survivals$log_reported_fatal),
                 tail_direct_survival(days_onset - day_onset + 0.5, 9, 6, q, 30),
                 tolerance=1e-10)
    probabilities <- path_probabilities(survivals)
    s <- f_survival(days_onset)
    expect_equal(probabilities$p, 1 - s / c(1, s[-length(s)]), tolerance=1e-10)
  }

  # the window of onset times whose reports are taken as complete uses the
  # share of the latest onset week: the onset times drawn are marked
  tail_window <- list(weight=c(0.01, 0.01, 0.01, 0.01, 0.1), mean=10, day_first=1,
                      rw_sd=0.1)
  observation_df <- dplyr::tibble(time_onset=1:30, time_reported=31,
                                  cases_reported=10, reporting_piece_index=1)
  drawn <- testthat::with_mocked_bindings(
    sample_cases_history(observation_df, 100, function(t) rep(1, length(t)),
                         list(mean=5, sd=3), dplyr::tibble(reporting_piece_index=1,
                                                           mean=5, sd=3),
                         revision_parameters=list(tail=tail_window),
                         f_observation_logp_onset=function(time) NULL),
    sample_true_cases_single_onset=function(...) -1)
  quantile_direct <- stats::uniroot(
    function(d) tail_direct_survival(d, 5, 3, 0.1, 10) - 0.01, c(0, 200),
    tol=1e-10)$root
  expect_equal(sort(unique(drawn$time_onset[drawn$cases_estimated == -1])),
               which(1:30 >= 30 - quantile_direct))
  # which differs from the window with the share of the other weeks
  quantile_other <- stats::uniroot(
    function(d) tail_direct_survival(d, 5, 3, 0.01, 10) - 0.01, c(0, 200),
    tol=1e-10)$root
  expect_gt(floor(quantile_direct), floor(quantile_other))

  # the package's simulator takes a share for each onset time
  set.seed(31)
  weight <- c(0.02, 0.3)
  sim <- simulate_cases_deaths(c(20000, 20000), 1:60, list(mean=4, sd=2), 0.5,
                               list(mean=3, sd=2),
                               reporting_tail=list(weight=weight, mean=20))
  checks <- expand.grid(t=1:2, lag=c(5, 15, 30, 50))
  for(i in seq_len(nrow(checks))) {
    t <- checks$t[i]
    day <- t + checks$lag[i]
    p <- 1 - tail_direct_survival(day - t + 0.5, 4, 2, weight[t], 20)
    x <- sim$snapshots$cases_reported[sim$snapshots$time_onset == t &
                                        sim$snapshots$time_reported == day]
    expect_z(x / 20000, p, sqrt(p * (1 - p) / 20000),
             sprintf("cases of onset day %d reported by day %d", t, day),
             zmax=z_bonferroni(nrow(checks)))
  }
  # and wrong records are entered as cases are first reported, tail and all
  sim <- simulate_cases_deaths(c(20000, 20000), 1:60, list(mean=4, sd=2), 0,
                               list(mean=3, sd=2),
                               reporting_tail=list(weight=weight, mean=20),
                               wrong=list(eta=0.05, rho=0.1, delta=0))
  for(t in 1:2) {
    day <- t + 40
    f_survival <- function(d) tail_direct_survival(d - t + 0.5, 4, 2, weight[t], 20)
    k <- (t + 1):day
    wrong_mean <- 0.05 * 20000 * sum((f_survival(k - 1) - f_survival(k)) * 0.9^(day - k))
    p <- 1 - f_survival(day)
    x <- sim$snapshots$cases_reported[sim$snapshots$time_onset == t &
                                        sim$snapshots$time_reported == day]
    expect_z(x, 20000 * p + wrong_mean, sqrt(20000 * p * (1 - p) + wrong_mean),
             sprintf("records of onset day %d in the report of day %d", t, day),
             zmax=z_bonferroni(2))
  }
  expect_error(simulate_cases_deaths(c(10, 10, 10), 1:10, list(mean=4, sd=2), 0.5,
                                     list(mean=3, sd=2),
                                     reporting_tail=list(weight=c(0.1, 0.2), mean=20)),
               "one for each onset time")
})

test_that("S8.2: every place the tail is read is known", {
  # each function of the package which reads a tail, or its weight: those of
  # the cases' reporting delay take the share of the onset week (checked in
  # the tests above and below) or are not reached with a share for each week
  # (cases alone, and the relaxed model); the rest read the tail of the delay
  # to a death being recorded, or the probability of a death unseen
  known <- c(
    # the cases' reporting delay, a share for each week allowed for
    "reporting.R: log_undetected_days", "reporting.R: tail_of_onset",
    "reporting.R: tail_week_of_days", "reporting.R: lag_cache_values",
    "reporting.R: delay_quantile", "sampling.R: sample_cases_history",
    "sampling.R: mcmc_single", "deaths-parameters.R: deaths_tail_step",
    "deaths-tail-walk.R: initial_tail_walk", "deaths-tail-walk.R: tail_walk_steps",
    "deaths-tail-walk.R: deaths_tail_walk_step", "deaths-tail-walk.R: tail_walk_draws",
    "predict.R: predict_state_from_draw", "deaths.R: simulate_cases_deaths",
    # the cases' reporting delay without a model of deaths, or its presence only
    "reporting.R: metropolis_step_tail",
    "probabilities.R: observation_process_revisions_logp",
    "simulated-data.R: generate_snapshots", "interface.R: fit_incidence_mcmc",
    "interface.R: fit_incidence_relaxed", "relaxed.R: relaxed_inits",
    "relaxed.R: relaxed_draws",
    # the delay to a death being recorded, and the deaths unseen
    "deaths-parameters.R: deaths_recording_tail_step", "deaths.R: deaths_draws",
    "deaths-dated.R: dated_logp", "deaths-dated.R: dated_apply_move",
    "deaths-dated.R: deaths_dated_step")
  directory <- test_path("..", "..", "R")
  skip_if_not(dir.exists(directory), "the source of the package is not here")
  found <- character(0)
  for(path in list.files(directory, pattern="[.]R$", full.names=TRUE)) {
    lines <- readLines(path)
    is_code <- !grepl("^\\s*#", lines)
    is_read <- is_code & grepl("tail\\$weight|\\$tail\\b|tail_weight\\b", lines)
    starts <- grep("^[A-Za-z_.][A-Za-z0-9_.]* <- function", lines)
    for(line in which(is_read)) {
      start <- max(starts[starts <= line])
      found <- c(found, paste0(basename(path), ": ",
                               sub(" <- function.*", "", lines[start])))
    }
  }
  expect_setequal(unique(found), known)
})

test_that("S8.3: the walk's likelihood is the change in the log-posterior", {
  data <- reference_data(101)
  fit <- tail_walk_fit(data, 2, seed=3)
  v <- fit$state$variables
  state <- v$deaths_state
  revision <- v$revision_current
  # the likelihood the update of the walk uses
  f_logl <- NULL
  testthat::with_mocked_bindings(
    deaths_tail_walk_step(state, v$deaths_cases_true, v$deaths_running_df,
                          revision, v$deaths_options),
    tail_walk_steps=function(tail, walk, prior, mean_step, f_logl_walk, n_steps) {
      f_logl <<- f_logl_walk
      tail
    })
  expect_true(is.function(f_logl))
  walk <- v$priors_reporting$tail_weight_walk
  tail <- tail_with_weeks(revision$tail, c(0.03, 0.08, 0.02, 0.12, 0.05))
  f_full <- function(tail) {
    proposed <- revision
    proposed$tail <- tail
    tail_full_logl(v, proposed) + tail_walk_prior_direct(tail, walk)
  }
  full <- f_full(tail)
  logl <- f_logl(tail)
  for(w in seq_along(tail$weight)) {
    changed <- tail
    changed$weight[w] <- stats::plogis(stats::qlogis(tail$weight[w]) + 0.7)
    expect_equal(f_logl(changed) - logl,
                 f_full(changed) - full -
                   (tail_walk_prior_direct(changed, walk) -
                      tail_walk_prior_direct(tail, walk)),
                 tolerance=1e-10, info=paste("week", w))
  }
  # and the mean of the tail, as its update uses it
  changed <- tail
  changed$mean <- 25
  expect_equal(f_logl(changed) - logl, f_full(changed) - full, tolerance=1e-10)
})

test_that("S8.4 (C4): the walk drawn with no likelihood follows its prior", {
  set.seed(41)
  walk <- tail_walk_settings(list(first_mean=-2, first_sd=0.8, step_sd_scale=0.3))
  prior <- list(tail_mean_mu=20, tail_mean_sigma=8)
  tail <- list(weight=stats::plogis(c(-2, -1.5, -2.5)), mean=20, day_first=1,
               rw_sd=0.2)
  n_thin <- 10
  n_draws <- 5000
  draws <- matrix(NA, n_draws, 5)
  for(i in seq_len(n_draws * n_thin)) {
    tail <- tail_walk_steps(tail, walk, prior, 0.6, function(tail) 0)
    if(i %% n_thin == 0) {
      x <- stats::qlogis(tail$weight)
      draws[i / n_thin, ] <- c(x[1], diff(x) / tail$rw_sd, tail$rw_sd, tail$mean)
    }
  }
  # the first share, the standardised steps, their sd and the mean of the tail
  expect_mcmc_mean(draws[, 1], -2, "logit of the first share")
  expect_mcmc_mean((draws[, 1] + 2)^2, 0.64, "variance of the logit of the first share")
  expect_ks(draws[, 1], function(x) stats::pnorm(x, -2, 0.8), "logit of the first share")
  for(j in 2:3)
    expect_ks(draws[, j], stats::pnorm, paste("standardised step", j - 1))
  half_normal_mean <- 0.3 * sqrt(2 / pi)
  expect_mcmc_mean(draws[, 4], half_normal_mean, "sd of the steps")
  expect_ks(draws[, 4], function(x) 2 * stats::pnorm(x, 0, 0.3) - 1, "sd of the steps")
  shape <- 20^2 / 8^2
  rate <- 20 / 8^2
  expect_mcmc_mean(draws[, 5], 20, "mean of the tail")
  expect_mcmc_mean((draws[, 5] - 20)^2, 64, "variance of the mean of the tail")
  expect_ks(draws[, 5], function(x) stats::pgamma(x, shape=shape, rate=rate),
            "mean of the tail")
})

test_that("S8.4 (C4): the walk drawn with a normal likelihood follows its posterior", {
  set.seed(42)
  # the sd of the steps held, so that the posterior is normal
  walk <- tail_walk_settings(list(first_mean=-2, first_sd=1, rw_sd_step=0))
  rw_sd <- 0.5
  y <- c(-1, -2.5, -1.8, -3)
  sigma <- 0.6
  n <- length(y)
  prior_cov <- 1 + rw_sd^2 * (outer(seq_len(n), seq_len(n), pmin) - 1)
  post_cov <- solve(solve(prior_cov) + diag(n) / sigma^2)
  post_mean <- as.vector(post_cov %*% (solve(prior_cov, rep(-2, n)) + y / sigma^2))
  tail <- list(weight=stats::plogis(rep(-2, n)), mean=20, day_first=1, rw_sd=rw_sd)
  f_logl <- function(tail) sum(stats::dnorm(y, stats::qlogis(tail$weight), sigma, log=TRUE))
  n_thin <- 5
  n_draws <- 5000
  draws <- matrix(NA, n_draws, n)
  for(i in seq_len(n_draws * n_thin)) {
    tail <- tail_walk_steps(tail, walk, list(tail_mean_mu=20, tail_mean_sigma=8),
                            0.5, f_logl)
    if(i %% n_thin == 0)
      draws[i / n_thin, ] <- stats::qlogis(tail$weight)
  }
  expect_equal(tail$rw_sd, rw_sd)
  zmax <- z_bonferroni(2 * n)
  for(j in seq_len(n)) {
    label <- paste("logit of the share of week", j)
    expect_z(mean(draws[, j]), post_mean[j], stat_batch_se(draws[, j]), label, zmax=zmax)
    squares <- (draws[, j] - post_mean[j])^2
    expect_z(mean(squares), post_cov[j, j], stat_batch_se(squares),
             paste("variance of the", label), zmax=zmax)
    expect_ks(draws[, j], function(x) stats::pnorm(x, post_mean[j], sqrt(post_cov[j, j])),
              label)
  }
})

test_that("S8.4 (C3): the update of the mean of the tail has the right acceptance ratio", {
  data <- reference_data(202)
  fit <- tail_walk_fit(data, 1, seed=5)
  v <- fit$state$variables
  revision <- v$revision_current
  prior <- v$priors_reporting
  f_logl <- function(tail) {
    proposed <- revision
    proposed$tail <- tail
    deaths_reporting_logp(v$deaths_state, v$deaths_cases_true, v$deaths_running_df,
                          proposed)
  }
  # the full log-posterior on the log scale of the mean
  f_full <- function(mean) {
    proposed <- revision
    proposed$tail$mean <- mean
    tail_full_logl(v, proposed) +
      stats::dgamma(mean, shape=prior$tail_mean_mu^2 / prior$tail_mean_sigma^2,
                    rate=prior$tail_mean_mu / prior$tail_mean_sigma^2, log=TRUE) +
      log(mean)
  }
  set.seed(6)
  current <- revision$tail$mean
  tail <- revision$tail
  logl <- f_logl(tail)
  full <- f_full(current)
  for(i in 1:25) {
    proposed <- current * exp(stats::rnorm(1, 0, 0.5))
    tail_proposed <- tail
    tail_proposed$mean <- proposed
    logl_proposed <- f_logl(tail_proposed)
    log_ratio <- tail_mean_log_ratio(proposed, current, logl_proposed, logl, prior)
    expect_equal(log_ratio, f_full(proposed) - full, tolerance=1e-8)
    # the reverse move
    expect_equal(tail_mean_log_ratio(current, proposed, logl, logl_proposed, prior),
                 -log_ratio)
  }
})

test_that("S8.5 (C2): with every share equal the model is the main model", {
  for(seed in c(101, 202)) {
    data <- reference_data(seed)
    fit <- reference_fit(data, FALSE, 1, seed=seed + 1)
    v <- fit$state$variables
    state <- v$deaths_state
    single <- v$revision_current
    weekly <- single
    weekly$tail <- initial_tail_walk(single$tail, state$onset_times,
                                     tail_walk_settings(TRUE))
    expect_equal(length(weekly$tail$weight), 5)
    delays <- deaths_delay_by_onset(state)
    dated <- dated_counts(state)
    for(k in seq_along(state$onset_times)) {
      days <- state$reports[[k]]$time_reported
      reporting <- delay_of_onset(v$deaths_running_df, k)
      f_terms <- function(revision) {
        survivals <- joint_survivals(days, state$onset_times[k], reporting,
                                     delay_of_onset(delays, k), revision,
                                     state$dated$reporting, state$recording_effects,
                                     state$nu, state$triggered)
        units <- wrong_entry_units(days, state$onset_times[k], reporting, revision)
        wrong <- list(lambda=revision$eta * v$deaths_cases_true[k] * units$unit,
                      w=units$w, delta=state$delta, max_erroneous=v$max_erroneous)
        reports <- state$reports[[k]]
        forward <- joint_path_exact_forward(
          reports$cases_reported, reports$deaths_reported,
          state$deaths_true[k] - dated$out[k], dated$into[k],
          v$deaths_cases_true[k] - state$deaths_true[k], survivals, wrong)
        list(survivals=survivals, units=units, logp=forward$logp)
      }
      expect_identical(f_terms(weekly), f_terms(single))
    }
    expect_identical(
      deaths_reporting_logp(state, v$deaths_cases_true, v$deaths_running_df, weekly),
      deaths_reporting_logp(state, v$deaths_cases_true, v$deaths_running_df, single))
    f_cases <- function(revision)
      deaths_cases_logp_function(deaths_refresh_survivals(state, v$deaths_running_df,
                                                          revision),
                                 state$onset_times[30], revision)(50:80)
    expect_identical(f_cases(weekly), f_cases(single))
  }
})

test_that("S8: a fit with the walk updates the shares and stores them", {
  data <- reference_data(101)
  fit <- tail_walk_fit(data, 3, seed=7)
  effects <- fit$reporting_effects
  v <- fit$state$variables
  tail <- v$revision_current$tail
  names_weeks <- paste0("tail_weight_week_", 1:5)
  expect_true(all(c(names_weeks, "tail_weight_rw_sd") %in% names(effects)))
  expect_true(all(is.na(effects$tail_weight)))
  # the draws of the last iteration are the state at its end
  last <- effects[effects$iteration == 3, ]
  expect_equal(unlist(last[names_weeks], use.names=FALSE), tail$weight)
  expect_equal(last$tail_mean, tail$mean)
  expect_equal(last$tail_weight_rw_sd, tail$rw_sd)
  expect_gt(length(unique(effects$tail_weight_week_2)), 1)
  expect_true(all(unlist(effects[names_weeks]) > 0 & unlist(effects[names_weeks]) < 1))
  # the survivals held are those of each onset time's own week
  state <- v$deaths_state
  for(k in c(3, 10, 20, 34)) {
    q <- tail$weight[tail_week_of_days(tail, state$onset_times[k])]
    single <- v$revision_current
    single$tail <- list(weight=q, mean=tail$mean)
    expect_equal(state$survivals[[k]]$log_reported,
                 log_undetected_days(state$reports[[k]]$time_reported,
                                     state$onset_times[k],
                                     delay_of_onset(v$deaths_running_df, k), single))
  }
})

test_that("C7: a fit with the walk taken up again makes the draws it would have made", {
  data <- reference_data(202)
  whole <- tail_walk_fit(data, 4, seed=8)
  first <- tail_walk_fit(data, 2, seed=8)
  set.seed(99)
  stats::runif(3)
  second <- tail_walk_fit(data, 2, seed=1, resume_from=first$state)
  names_draws <- setdiff(names(whole), c("state", "iterations_done"))
  for(name in names_draws)
    expect_equal(dplyr::bind_rows(first[[name]], second[[name]]), whole[[name]],
                 info=name)
  expect_equal(second$state$variables$revision_current$tail,
               whole$state$variables$revision_current$tail)
})

test_that("C8: predictions read the share of each onset week", {
  set.seed(61)
  n_days <- 30
  serial_parameters <- list(mean=5, sd=3)
  cases <- true_cases(n_days + 40, function(t) 1, kappa=1000, serial_parameters,
                      initial_parameters=list(mean=30, length=10))
  sim <- simulate_cases_deaths(cases[41:(n_days + 40)], seq_len(n_days + 30),
                               list(mean=6, sd=4), 0.5, list(mean=3, sd=2.5),
                               triggered=1, onset_to_death=list(mean=9.3, sd=5.5),
                               dated_by_death=list(fraction=0.3),
                               recording_backlog=list(share=0.3, rate=0.05),
                               wrong=list(eta=0.02, rho=0.1, delta=0.3),
                               reporting_tail=list(weight=rep(c(0.02, 0.08), c(15, 15)),
                                                   mean=30))
  snapshots <- dplyr::filter(sim$snapshots, time_reported <= 40)
  data <- list(snapshots=snapshots, cases_before=cases[1:40],
               serial_parameters=serial_parameters,
               initial_cases=snapshots %>%
                 dplyr::group_by(time_onset) %>%
                 dplyr::summarise(cases_true=max(cases_reported) + 2))
  f_draw <- function(fit, iteration)
    lapply(fit[c("cases", "reporting_days", "revisions", "reporting_effects",
                 "deaths", "deaths_parameters")],
           function(df) df[df$iteration == iteration, ])
  fit <- tail_walk_fit(data, 3, seed=62, is_uncorrected=FALSE)
  v <- fit$state$variables
  draw <- f_draw(fit, 3)
  prediction <- predict_state_from_draw(v, draw)
  expect_equal(prediction$revision_parameters$tail, v$revision_current$tail)
  # at an earlier draw, the shares of that draw
  earlier <- predict_state_from_draw(v, f_draw(fit, 2))$revision_parameters$tail
  expect_equal(earlier$weight,
               unlist(fit$reporting_effects[2, paste0("tail_weight_week_", 1:5)],
                      use.names=FALSE))
  # each onset time's survivals take the share of its week (those held in
  # the state are from before the last update of the reporting delay)
  state <- v$deaths_state
  tail <- prediction$revision_parameters$tail
  for(k in c(2, 12, 29)) {
    days <- state$reports[[k]]$time_reported
    terms <- predict_onset_terms(prediction, k, days)
    single <- prediction$revision_parameters
    single$tail <- list(weight=tail$weight[tail_week_of_days(tail, state$onset_times[k])],
                        mean=tail$mean)
    expect_equal(terms$survivals$log_reported,
                 log_undetected_days(days, state$onset_times[k],
                                     delay_of_onset(prediction$reporting_by_onset, k),
                                     single))
    expect_false(isTRUE(all.equal(
      terms$survivals$log_reported,
      log_undetected_days(days, state$onset_times[k],
                          delay_of_onset(prediction$reporting_by_onset, k),
                          tail_other_week(single, tail, state$onset_times[k])))))
  }
  # later reports, with the share of the last onset week for the latest
  # onset time
  out <- predict_reports(v, draw, 41:44)
  expect_equal(sort(unique(out$time_reported)), 41:44)
  expect_true(all(out$deaths_reported <= out$cases_reported))
  extended <- predict_extend_reports(prediction, 41:44)
  terms <- predict_onset_terms(extended, 30, c(state$reports[[30]]$time_reported, 41:44))
  single <- extended$revision_parameters
  single$tail <- list(weight=dplyr::last(single$tail$weight), mean=single$tail$mean)
  expect_equal(terms$survivals$log_reported,
               log_undetected_days(c(state$reports[[30]]$time_reported, 41:44), 30,
                                   delay_of_onset(prediction$reporting_by_onset, 30),
                                   single))
  # at the nested case, the terms are those of the main model
  nested <- prediction
  nested$revision_parameters$tail <- tail_with_weeks(prediction$revision_parameters$tail,
                                                     0.05)
  main <- prediction
  main$revision_parameters$tail <- list(weight=0.05, mean=prediction$revision_parameters$tail$mean)
  for(k in c(2, 29))
    expect_identical(predict_onset_terms(nested, k, state$reports[[k]]$time_reported),
                     predict_onset_terms(main, k, state$reports[[k]]$time_reported))
})
