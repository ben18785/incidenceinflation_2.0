# A fraction of fatal cases dated by death which varies by onset week (S9 of
# the sensitivity analyses): a random walk on the logit scale over weeks of
# onset, drawn whole by elliptical slice sampling given which fatal cases are
# dated by death.

# a state of the model of deaths in which the fraction dated by death
# differs between onset weeks, with earlier cases known so that some fatal
# cases recorded for the first days had their onset before the first day
walk_test_state <- function(logit_fraction=stats::qlogis(c(0.2, 0.5, 0.3, 0.6)),
                            exact=FALSE, n_days=25, seed=31) {
  set.seed(seed)
  cases <- rep(40, n_days)
  week_of_day <- (seq_len(n_days) - 1) %/% 7 + 1
  fraction <- stats::plogis(logit_fraction[week_of_day])
  sim <- simulate_cases_deaths(
    cases, seq_len(n_days + 5), list(mean=5, sd=3), 0.4, list(mean=9, sd=4),
    dated_by_death=list(fraction=fraction, mean=5, sd=3,
                        reporting=list(mean=2, sd=1)))
  settings <- deaths_settings(n_days, list(dated_by_death=list(
    fraction=list(walk=list(initial=fraction)), mean=5, sd=3,
    reporting=list(mean=2, sd=1), exact=exact)))
  initial_cases <- sim$snapshots %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(cases_true=max(cases_reported) + 5)
  cases_before <- rep(30, 20)
  state <- initial_deaths_state(sim$snapshots, initial_cases, settings, cases_before)
  reporting_by_onset <- dplyr::tibble(mean=rep(5, n_days), sd=rep(3, n_days))
  cases_true <- state$cases_onset
  state <- deaths_paths_step(state, cases_true, reporting_by_onset)
  list(state=state, cases_true=cases_true, settings=settings,
       recorded=initial_cases$cases_true, reporting_by_onset=reporting_by_onset,
       renewal=dated_renewal(rep(1, n_days), 5, TRUE, list(mean=5, sd=3), 20,
                             cases_before))
}

# the same state with the walk and the fraction of each onset time replaced
with_walk <- function(state, logit_fraction) {
  state$dated$walk$logit_fraction <- logit_fraction
  state$dated$fraction <- stats::plogis(logit_fraction[state$dated$walk$week_of_day])
  state
}

# the log-probability of the fatal cases dated by death, written from the
# model (appendix, 'The joint probability of the labels'): for each onset day
# the multinomial over dated by onset and dated by death dying after each
# number of days (or after the last day), with the fraction of the week of
# the onset day, and the binomial of its fatal cases; for each recorded day
# the Poisson probability of those with onset before the first day, with the
# fraction of the first week, and its reports
reference_dated_logp <- function(state, cases_true, logit_fraction) {
  dated <- state$dated
  n_days <- length(state$onset_times)
  f <- stats::plogis(logit_fraction[(seq_len(n_days) - 1) %/% 7 + 1])
  h <- dated$lag_probabilities
  logp <- 0
  for(k in seq_len(n_days)) {
    v <- dated$lags[k, ]
    u <- dated$unseen[k]
    fatal <- state$deaths_true[k]
    correct <- fatal - sum(v) - u
    if(correct < 0 | fatal > cases_true[k])
      return(-Inf)
    logp <- logp + lgamma(fatal + 1) - lgamma(correct + 1) - sum(lgamma(v + 1)) -
      lgamma(u + 1) + correct * log(1 - f[k]) + sum(v[v > 0] * log(f[k] * h[v > 0])) +
      (if(u > 0) u * log(f[k] * dated$tail[k]) else 0) +
      stats::dbinom(fatal, cases_true[k], state$cfr[k], log=TRUE)
  }
  for(t in seq_len(n_days)) {
    logp <- logp + stats::dpois(dated$inflow[t],
                                dated$inflow_weight[t] * state$cfr[1] * f[1], log=TRUE) +
      fatal_reports_logp(state$paths[[t]], state$survivals[[t]],
                         state$deaths_true[t] - dated$out[t] + dated$into[t],
                         dated$into[t])
  }
  logp
}

test_that("the settings of a fraction which varies by onset week are checked", {
  base <- list(mean=5, sd=3)
  settings <- dated_by_death_settings(20, c(list(fraction=list(walk=TRUE)), base))
  walk <- settings$walk
  expect_equal(walk[c("first_mean", "first_sd", "step_sd_scale", "week")],
               list(first_mean=0, first_sd=1.6, step_sd_scale=0.3, week=7))
  expect_equal(walk$week_of_day, rep(1:3, c(7, 7, 6)))
  expect_equal(walk$logit_fraction, rep(stats::qlogis(0.2), 3))
  expect_equal(settings$fraction, rep(0.2, 20))
  expect_null(settings$fraction_prior)
  # a starting value for each onset day: the mean of each week is taken
  initial <- rep(c(0.1, 0.3, 0.5), c(7, 7, 6))
  initial[1] <- 0.2
  settings <- dated_by_death_settings(20, c(list(fraction=list(walk=list(
    initial=initial, first_sd=1, step_sd_scale=0.1, initial_step_sd=0.05))), base))
  expect_equal(settings$walk$logit_fraction,
               stats::qlogis(c(mean(initial[1:7]), 0.3, 0.5)))
  expect_equal(settings$fraction, stats::plogis(settings$walk$logit_fraction)[
    rep(1:3, c(7, 7, 6))])
  expect_equal(settings$walk$step_sd, 0.05)
  expect_equal(settings$walk$first_sd, 1)
  # without a walk the settings are as they were, with nothing added
  expect_false("walk" %in% names(dated_by_death_settings(
    20, c(list(fraction=list(alpha=1, beta=1)), base))))
  expect_false("walk" %in% names(dated_by_death_settings(20, c(list(fraction=0.3), base))))

  f_error <- function(fraction, ...)
    dated_by_death_settings(20, c(list(fraction=fraction), base, list(...)))
  expect_error(f_error(list(walk=TRUE, alpha=2)), "set in 'walk'")
  expect_error(f_error(list(walk=list(other=1))), "Unknown settings for the walk")
  expect_error(f_error(list(walk="yes")), "TRUE or a named list")
  expect_error(f_error(list(walk=list(first_sd=-1))), "must be positive")
  expect_error(f_error(list(walk=list(step_sd_scale=0))), "must be positive")
  expect_error(f_error(list(walk=list(first_sd=c(1, 2)))), "single number")
  expect_error(f_error(list(walk=list(week=3.5))), "whole numbers")
  expect_error(f_error(list(walk=list(initial=1.2))), "strictly between 0 and 1")
  expect_error(f_error(list(walk=list(initial=0))), "strictly between 0 and 1")
  expect_error(f_error(list(walk=list(initial=rep(0.2, 7)))), "one for each onset time")
  # what is not allowed for with a fraction which varies
  expect_error(f_error(list(walk=TRUE), collapse_fraction=TRUE), "cannot be integrated out")
  expect_error(f_error(list(walk=TRUE), exact=TRUE, ridge=TRUE),
               "not allowed for with a fraction which varies")
  expect_error(relaxed_stan_data(
    dplyr::tibble(time_onset=1, time_reported=1, cases_reported=1), list(mean=5, sd=3),
    deaths_df=dplyr::tibble(time_onset=1, time_reported=1, deaths_reported=0),
    dated_by_death=c(list(fraction=list(walk=TRUE)), base)),
    "varies by onset week")
})

test_that("each onset day has the fraction of its week, and the cases from before the first day that of the first week", {
  test <- walk_test_state()
  state <- test$state
  cases_true <- test$cases_true
  logit_fraction <- state$dated$walk$logit_fraction
  expect_equal(state$dated$walk$week_of_day, rep(1:4, c(7, 7, 7, 4)))
  expect_equal(state$dated$fraction,
               stats::plogis(logit_fraction)[rep(1:4, c(7, 7, 7, 4))])
  expect_true(sum(state$dated$inflow) > 0)
  # the probability of the fatal cases recorded for the first days with onset
  # before the first day changes with the fraction of the first week only
  days <- which(state$dated$inflow_weight > 0)
  f_inflow <- function(s) dated_logp(s, cases_true, integer(0), days)
  other <- with_walk(state, logit_fraction + c(0, 1, -1, 2))
  expect_equal(f_inflow(other), f_inflow(state))
  first <- with_walk(state, logit_fraction + c(0.5, 0, 0, 0))
  mean_inflow <- state$dated$inflow_weight[days] * state$cfr[1]
  expect_equal(f_inflow(first) - f_inflow(state),
               sum(stats::dpois(state$dated$inflow[days],
                                mean_inflow * first$dated$fraction[1], log=TRUE) -
                     stats::dpois(state$dated$inflow[days],
                                  mean_inflow * state$dated$fraction[1], log=TRUE)))
  # and dated_logp is that of the model written out, for the fraction of
  # each onset day's week
  expect_equal(dated_logp(state, cases_true),
               reference_dated_logp(state, cases_true, logit_fraction))
})

test_that("with the same fraction in every week the model is the main model", {
  test <- walk_test_state(logit_fraction=rep(stats::qlogis(0.35), 4))
  state <- test$state
  cases_true <- test$cases_true
  single <- state
  single$dated$walk <- NULL
  single$dated$fraction <- rep(0.35, 25)
  single$dated$fraction_prior <- list(alpha=1, beta=1)
  expect_equal(state$dated$fraction, single$dated$fraction)
  expect_equal(dated_logp(state, cases_true), dated_logp(single, cases_true),
               tolerance=1e-10)
  expect_equal(deaths_log_likelihood(state, cases_true),
               deaths_log_likelihood(single, cases_true), tolerance=1e-10)
  for(k in c(3, 12, 24))
    expect_equal(deaths_cases_logp_function(state, k)(cases_true[k] + 0:3),
                 deaths_cases_logp_function(single, k)(cases_true[k] + 0:3),
                 tolerance=1e-10)
  # the walk's log-likelihood is then the log-density of the single
  # fraction given the labels (dated_fraction_terms, with the uniform prior),
  # up to a constant
  terms <- dated_fraction_terms(single)
  counts <- dated_walk_counts(state)
  f_single <- function(f) (terms$A - 1) * log(f) + (terms$B - 1) * log(1 - f) - terms$r * f
  differences <- sapply(c(0.1, 0.35, 0.7), function(f)
    dated_walk_logl(rep(stats::qlogis(f), 4), counts) - f_single(f))
  expect_equal(differences, rep(differences[1], 3), tolerance=1e-10)
  expect_equal(differences[1], 0, tolerance=1e-10)
})

test_that("the counts by week are those of the labels, and sum to those of dated_fraction_terms", {
  test <- walk_test_state()
  state <- test$state
  dated <- state$dated
  counts <- dated_walk_counts(state)
  # counted from the labels: the fatal cases of each onset day dated by death
  # (dying on a day with reports, or after the last) and the rest
  week <- (seq_len(25) - 1) %/% 7 + 1
  n_dated <- rowSums(dated$lags) + dated$unseen
  expect_equal(counts$dated, sapply(1:4, function(w) sum(n_dated[week == w])))
  expect_equal(counts$by_onset,
               sapply(1:4, function(w) sum(state$deaths_true[week == w] - n_dated[week == w])))
  expect_equal(counts$inflow, sum(dated$inflow))
  expect_true(all(counts$dated > 0) & all(counts$by_onset > 0))
  state$dated$fraction_prior <- list(alpha=1, beta=1)
  terms <- dated_fraction_terms(state)
  expect_equal(sum(counts$dated) + counts$inflow, terms$A - 1)
  expect_equal(sum(counts$by_onset), terms$B - 1)
  expect_equal(counts$r, terms$r)
  # a week with no onset days counts nothing
  state$dated$walk$logit_fraction <- c(state$dated$walk$logit_fraction, 0)
  expect_equal(dated_walk_counts(state)$dated[5], 0)
})

test_that("changes to the fatal cases dated by death give the right change in log-probability with fractions which differ by week", {
  for(exact in c(FALSE, TRUE)) {
    test <- walk_test_state(exact=exact, seed=32 + exact)
    state <- test$state
    cases_true <- test$cases_true
    logit_fraction <- state$dated$walk$logit_fraction
    week <- state$dated$walk$week_of_day
    # the full log-posterior, from scratch: the renewal model and the parts of
    # the model which the changes alter (the walk's prior does not change)
    f_full <- function(state, cases_true)
      (if(exact) dated_logp(state, cases_true)
       else reference_dated_logp(state, cases_true, logit_fraction)) +
      dated_renewal_logp(cases_true, test$renewal)
    set.seed(34)
    n_applied <- c(relabel=0, redate=0, unseen=0, inflow=0)
    n_across <- 0
    for(i in 1:(if(exact) 150 else 400)) {
      type <- sample(names(n_applied), 1)
      move <- list(type=type, day=sample.int(25, 1),
                   lag=sample.int(state$dated$max_lag, 1),
                   lag_to=sample.int(state$dated$max_lag, 1), sign=sample(c(-1, 1), 1))
      proposed <- dated_apply_move(state, cases_true, move)
      if(is.null(proposed))
        next
      logp_proposed <- f_full(proposed$state, proposed$cases_true)
      if(!is.finite(logp_proposed))
        next
      # the log of the acceptance ratio, as deaths_dated_step computes it
      log_ratio <- dated_logp(proposed$state, proposed$cases_true, proposed$cohorts,
                              proposed$days) +
        dated_renewal_logp(proposed$cases_true, test$renewal) -
        dated_logp(state, cases_true, proposed$cohorts, proposed$days) -
        dated_renewal_logp(cases_true, test$renewal)
      expect_equal(log_ratio, logp_proposed - f_full(state, cases_true), tolerance=1e-10)
      # the reverse change returns to the state
      reverse <- move
      reverse$sign <- -move$sign
      if(type == "redate") {
        reverse$lag <- move$lag_to
        reverse$lag_to <- move$lag
      }
      back <- dated_apply_move(proposed$state, proposed$cases_true, reverse)
      expect_equal(back$state$dated[c("lags", "unseen", "inflow", "out", "into")],
                   state$dated[c("lags", "unseen", "inflow", "out", "into")])
      expect_equal(back$cases_true, cases_true)
      n_across <- n_across + (length(unique(week[proposed$cohorts])) > 1)
      state <- proposed$state
      cases_true <- proposed$cases_true
      n_applied[type] <- n_applied[type] + 1
    }
    expect_true(all(n_applied > 3))
    # some of the changes move a case between onset days in different weeks
    expect_gt(n_across, 5)
    expect_equal(cases_true - state$dated$out + state$dated$into, test$recorded)
  }
  # the step itself runs and leaves the walk alone
  test <- walk_test_state()
  set.seed(35)
  stepped <- deaths_dated_step(test$state, test$cases_true, test$renewal)
  expect_false(identical(stepped$cases_true, test$cases_true))
  expect_equal(stepped$state$dated$fraction, test$state$dated$fraction)
  expect_equal(stepped$cases_true - stepped$state$dated$out + stepped$state$dated$into,
               test$recorded)
})

test_that("the updates of the fatal cases use the fraction of the onset day's week", {
  test <- walk_test_state()
  state <- test$state
  cases_true <- test$cases_true
  for(k in c(4, 17)) {
    # a single fraction, at the value of the week of onset day k
    single <- state
    single$dated$walk <- NULL
    single$dated$fraction <- rep(state$dated$fraction[k], 25)
    values <- cases_true[k] + 0:4
    expect_equal(deaths_cases_logp_function(state, k)(values),
                 deaths_cases_logp_function(single, k)(values))
    expect_equal(joint_cases_deaths_logp(
      state$paths[[k]], state$survivals[[k]], state$cfr[k], cases_true[k],
      state$dated$out[k], state$dated$into[k], dated_counts(state)$fraction[k]),
      joint_cases_deaths_logp(
        single$paths[[k]], single$survivals[[k]], single$cfr[k], cases_true[k],
        single$dated$out[k], single$dated$into[k], dated_counts(single)$fraction[k]))
  }
  # and the fraction of another week makes no difference
  other <- with_walk(state, state$dated$walk$logit_fraction + c(0, 0, 1, 0))
  expect_equal(deaths_cases_logp_function(other, 4)(cases_true[4] + 0:4),
               deaths_cases_logp_function(state, 4)(cases_true[4] + 0:4))
})

# a state of the walk alone, with the labels given as counts: n_days onset
# days of fatal cases, some dated by death, and those from before the first
# day recorded for day 1
walk_counts_state <- function(n_dated, n_fatal, inflow=0, inflow_weight=0,
                              cfr=0.4, walk=list()) {
  n_days <- length(n_dated)
  settings <- dated_by_death_settings(n_days, list(fraction=list(walk=walk),
                                                   mean=5, sd=3))
  list(onset_times=seq_len(n_days), deaths_true=n_fatal, cfr=rep(cfr, n_days),
       dated=c(settings, list(out=n_dated,
                              inflow=c(inflow, rep(0, n_days - 1)),
                              inflow_weight=c(inflow_weight, rep(0, n_days - 1)))))
}

# the walk and the standard deviation of its steps after n_steps updates
# from each of a set of starting values (a row of the walk followed by the
# standard deviation)
walk_draws <- function(state, starts, n_steps) {
  n_weeks <- ncol(starts) - 1
  t(apply(starts, 1, function(start) {
    state$dated$walk$logit_fraction <- start[seq_len(n_weeks)]
    state$dated$walk$step_sd <- start[n_weeks + 1]
    for(step in seq_len(n_steps))
      state <- deaths_dated_fraction_step(state)
    c(state$dated$walk$logit_fraction, state$dated$walk$step_sd)
  }))
}

# draws which are from the distribution they started from, and have moved
expect_walk_draws <- function(draws, starts, exact_cdfs, label) {
  names <- c("logit f_1", "logit f_2", "tau_f")
  for(i in 1:3) {
    what <- paste(label, names[i])
    expect_ks(draws[, i], exact_cdfs[[i]], what)
    expect_ks(starts[, i], exact_cdfs[[i]], paste(what, "(starting values)"))
    expect_lt(stats::cor(draws[, i], starts[, i]), 0.9, label=what)
  }
}

test_that("the update of the walk draws from its distribution given the labels, against a grid for two weeks", {
  # two weeks of seven onset days, with fixed labels: in the first week 6 of
  # 35 fatal cases are dated by death and in the second 16 of 30, and 2 from
  # before the first day, with 5 expected for each unit of cfr times fraction
  n_fatal <- rep(c(5, 4), c(7, 7)) + c(rep(0, 13), 2)
  n_dated <- c(1, 1, 1, 1, 1, 1, 0, 2, 2, 2, 2, 2, 2, 4)
  state <- walk_counts_state(n_dated, n_fatal, inflow=2, inflow_weight=5)
  expect_equal(sapply(1:2, function(w) sum(n_dated[state$dated$walk$week_of_day == w])),
               c(6, 16))
  # the posterior on a grid, from the model: x1 = logit f_1, x2 = x1 + tau u,
  # with u standard normal, so that the grid is smooth at small tau
  f_logl <- function(x, n_dated, n_correct)
    n_dated * stats::plogis(x, log.p=TRUE) +
    n_correct * stats::plogis(x, lower.tail=FALSE, log.p=TRUE)
  spacing <- c(x1=0.025, u=0.1, tau=0.005)
  grid <- expand.grid(x1=seq(-4.5, 1.5, spacing[["x1"]]), u=seq(-6, 6, spacing[["u"]]),
                      tau=seq(spacing[["tau"]] / 2, 1.6, spacing[["tau"]]))
  f_log_posterior <- function(x1, u, tau) {
    x2 <- x1 + tau * u
    f_logl(x1, 6 + 2, 29) - 5 * 0.4 * stats::plogis(x1) + f_logl(x2, 16, 14) +
      stats::dnorm(x1, 0, 1.6, log=TRUE) + stats::dnorm(u, log=TRUE) +
      stats::dnorm(tau, 0, 0.3, log=TRUE)
  }
  log_w <- f_log_posterior(grid$x1, grid$u, grid$tau)
  w <- exp(log_w - max(log_w))
  w <- w / sum(w)
  f_cdf <- function(x) {
    order <- order(x)
    stats::approxfun(x[order], cumsum(w[order]), yleft=0, yright=1, ties="ordered")
  }
  values <- list(grid$x1, grid$x1 + grid$tau * grid$u, grid$tau)
  # starting values drawn from the grid, each within its cell
  set.seed(41)
  n <- 4000
  cells <- grid[sample.int(nrow(grid), n, replace=TRUE, prob=w), ]
  for(name in names(spacing))
    cells[[name]] <- cells[[name]] + stats::runif(n, -0.5, 0.5) * spacing[[name]]
  starts <- cbind(cells$x1, cells$x1 + cells$tau * cells$u, cells$tau)
  draws <- walk_draws(state, starts, 10)
  for(i in 1:3) {
    label <- c("logit f_1", "logit f_2", "tau_f")[i]
    exact_mean <- sum(w * values[[i]])
    expect_mean_z(draws[, i], exact_mean, paste(label, "mean"))
    expect_mean_z((draws[, i] - exact_mean)^2, sum(w * (values[[i]] - exact_mean)^2),
                  paste(label, "variance"))
  }
  expect_walk_draws(draws, starts, lapply(values, f_cdf), "given the labels:")

  # with no labels, nothing is known of the fractions: the draws are from the
  # prior
  state <- walk_counts_state(rep(0, 14), rep(0, 14))
  set.seed(42)
  n <- 3000
  tau <- abs(stats::rnorm(n, 0, 0.3))
  x1 <- stats::rnorm(n, 0, 1.6)
  starts <- cbind(x1, x1 + tau * stats::rnorm(n), tau)
  draws <- walk_draws(state, starts, 10)
  expect_mean_z(draws[, 1], 0, "logit f_1 under the prior")
  expect_mean_z(draws[, 1]^2, 1.6^2, "variance of logit f_1 under the prior")
  expect_mean_z(draws[, 3], 0.3 * sqrt(2 / pi), "tau_f under the prior")
  expect_mean_z(draws[, 2]^2, 1.6^2 + 0.3^2, "variance of logit f_2 under the prior")
  # logit f_2 is logit f_1 plus a step whose sd is half-normal
  f_cdf_2 <- function(x) sapply(x, function(y) stats::integrate(function(s)
    2 * stats::dnorm(s, 0, 0.3) * stats::pnorm(y, 0, sqrt(1.6^2 + s^2)), 0, Inf)$value)
  expect_walk_draws(draws, starts,
                    list(function(x) stats::pnorm(x, 0, 1.6), f_cdf_2,
                         function(x) 2 * stats::pnorm(x, 0, 0.3) - 1),
                    "under the prior:")
  expect_ks((draws[, 2] - draws[, 1]) / draws[, 3], stats::pnorm,
            "standardised step under the prior")
})

# the state with a single known fraction in place of the walk
dated_by_death_state_single <- function(state) {
  state$dated$walk <- NULL
  state
}

test_that("the update of the walk keeps the state consistent and leaves the labels alone", {
  test <- walk_test_state()
  set.seed(43)
  stepped <- deaths_dated_fraction_step(test$state)
  walk <- stepped$dated$walk
  expect_false(identical(walk$logit_fraction, test$state$dated$walk$logit_fraction))
  expect_equal(stepped$dated$fraction, stats::plogis(walk$logit_fraction[walk$week_of_day]))
  expect_equal(stepped$dated[c("lags", "unseen", "inflow", "out", "into")],
               test$state$dated[c("lags", "unseen", "inflow", "out", "into")])
  # the draws hold the fraction of each onset day and the sd of the steps
  draws <- deaths_draws(stepped, 1)
  expect_equal(draws$deaths$dated_fraction, stepped$dated$fraction)
  expect_equal(draws$deaths_parameters$dated_fraction, stepped$dated$fraction[1])
  expect_equal(draws$deaths_parameters$dated_fraction_step_sd, walk$step_sd)
  # without the walk there are no such columns
  single <- deaths_draws(dated_by_death_state_single(stepped), 1)
  expect_false("dated_fraction" %in% names(single$deaths))
  expect_false("dated_fraction_step_sd" %in% names(single$deaths_parameters))
})

# a small fit with the fraction varying by onset week
walk_fit_example <- function() {
  set.seed(44)
  n_days <- 25
  serial_parameters <- list(mean=5, sd=3)
  cases <- true_cases(n_days + 40, function(t) 1, kappa=1000, serial_parameters,
                      initial_parameters=list(mean=30, length=10))
  fraction <- rep(c(0.2, 0.5), c(14, 11))
  dated_by_death <- list(fraction=fraction, mean=6, sd=3, reporting=list(mean=2, sd=1))
  sim <- simulate_cases_deaths(cases[41:(n_days + 40)], seq_len(n_days + 20),
                               list(mean=4, sd=3), 0.5, list(mean=8, sd=4),
                               dated_by_death=dated_by_death)
  dated_by_death$fraction <- list(walk=list(initial=0.3))
  list(snapshots=sim$snapshots, cases_before=cases[1:40],
       serial_parameters=serial_parameters, dated_by_death=dated_by_death,
       initial_cases=sim$snapshots %>%
         dplyr::group_by(time_onset) %>%
         dplyr::summarise(cases_true=max(cases_reported)))
}

walk_fit <- function(example, niterations, ...) {
  snapshots <- example$snapshots
  mcmc_single(
    niterations=niterations,
    snapshot_with_Rt_index_df=snapshots %>%
      dplyr::select(time_onset, time_reported, cases_reported) %>%
      dplyr::mutate(Rt_index=1),
    priors=list(Rt=list(shape=1, rate=1),
                reporting=list(mean_mu=5, mean_sigma=5, sd_mu=3, sd_sigma=3),
                max_cases=300,
                Rt_changepoints=list(),
                reporting_changepoints=list(min_length=14)),
    serial_parameters=example$serial_parameters,
    initial_cases_true=example$initial_cases,
    initial_reporting_parameters=list(mean=4, sd=3),
    initial_Rt=dplyr::tibble(Rt_index=1, Rt=1),
    print_to_screen=FALSE,
    cases_before=example$cases_before,
    deaths=list(snapshot=snapshots %>%
                  dplyr::select(time_onset, time_reported, deaths_reported),
                settings=list(dated_by_death=example$dated_by_death,
                              initial_delay=list(mean=8, sd=4))),
    ...)
}

test_that("a fit with a fraction which varies by onset week runs and can be taken up again", {
  example <- walk_fit_example()
  set.seed(45)
  whole <- walk_fit(example, 6)
  set.seed(45)
  first <- walk_fit(example, 3, return_state=TRUE)
  set.seed(99)
  stats::runif(5)
  second <- walk_fit(example, 3, resume_from=first$state)
  for(name in setdiff(names(whole), c("state", "iterations_done")))
    expect_equal(dplyr::bind_rows(first[[name]], second[[name]]), whole[[name]],
                 info=name)
  # each onset day has its week's fraction, which moves
  expect_true("dated_fraction" %in% names(whole$deaths))
  expect_true("dated_fraction_step_sd" %in% names(whole$deaths_parameters))
  last <- dplyr::filter(whole$deaths, iteration == 6)
  week <- (last$time_onset - 1) %/% 7 + 1
  expect_equal(as.vector(tapply(last$dated_fraction, week, function(x) diff(range(x)))),
               rep(0, 4))
  expect_gt(length(unique(last$dated_fraction)), 1)
  expect_equal(whole$deaths_parameters$dated_fraction,
               dplyr::filter(whole$deaths, time_onset == 1)$dated_fraction)
  expect_equal(length(unique(whole$deaths_parameters$dated_fraction_step_sd)) > 1, TRUE)
})

# reports for predictions, of the kind fitted (as in test-predict.R), and a
# short fit to them with the exact form of fatal cases dated by death
walk_predict_example <- function() {
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
  snapshots <- dplyr::filter(sim$snapshots, time_reported <= 40)
  list(snapshots=snapshots, cases_before=cases[1:40],
       serial_parameters=serial_parameters,
       initial_cases=snapshots %>%
         dplyr::group_by(time_onset) %>%
         dplyr::summarise(cases_true=max(cases_reported) + 2))
}

walk_predict_fit <- function(example, fraction) {
  snapshots <- example$snapshots
  mcmc_single(
    niterations=3,
    snapshot_with_Rt_index_df=snapshots %>%
      dplyr::select(time_onset, time_reported, cases_reported) %>%
      dplyr::mutate(Rt_index=1),
    priors=list(Rt=list(shape=1, rate=1),
                reporting=list(mean_mu=6, mean_sigma=5, sd_mu=4, sd_sigma=4),
                max_cases=300, overdispersion=list(mean=20, sd=20),
                revisions=list(eta_mu=0.1, eta_sigma=0.1, rho_alpha=1, rho_beta=1),
                max_erroneous=4, Rt_changepoints=list(),
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
                    dated_by_death=list(fraction=fraction, exact=TRUE))))
}

walk_predict_draw <- function(fit, iteration) {
  lapply(fit[c("cases", "reporting_days", "revisions", "reporting_effects",
               "deaths", "deaths_parameters")],
         function(df) df[df$iteration == iteration, ])
}

test_that("predictions read the fraction of each onset week, and are those of the main model when it is the same in every week", {
  example <- walk_predict_example()
  set.seed(46)
  fit <- walk_predict_fit(example, list())
  v <- fit$state$variables
  draw <- walk_predict_draw(fit, 3)
  # the main model, with a walk added at the same fraction in every week
  walk_variables <- v
  dated <- v$deaths_state$dated
  walk <- dated_walk_settings(length(dated$fraction), TRUE)
  walk$logit_fraction[] <- stats::qlogis(dated$fraction[1])
  walk_variables$deaths_state$dated$walk <- walk
  walk_draw <- draw
  walk_draw$deaths$dated_fraction <- dated$fraction[1]
  walk_draw$deaths_parameters$dated_fraction_step_sd <- 0.15
  prediction <- predict_state_from_draw(walk_variables, walk_draw)
  expect_equal(prediction$deaths_state$dated$fraction, dated$fraction)
  expect_equal(prediction$deaths_state$dated$walk$step_sd, 0.15)
  set.seed(47)
  main <- predict_reports(v, draw, 41:45)
  set.seed(47)
  expect_identical(predict_reports(walk_variables, walk_draw, 41:45), main)

  # a fit with the walk: the state at a kept draw has that draw's fractions
  set.seed(48)
  walk_fit <- walk_predict_fit(example, list(walk=TRUE))
  v <- walk_fit$state$variables
  for(iteration in 2:3) {
    draw <- walk_predict_draw(walk_fit, iteration)
    prediction <- predict_state_from_draw(v, draw)
    expect_equal(prediction$deaths_state$dated$fraction, draw$deaths$dated_fraction)
    expect_equal(prediction$deaths_state$dated$walk$step_sd,
                 draw$deaths_parameters$dated_fraction_step_sd)
  }
  expect_equal(prediction$deaths_state$dated$walk$logit_fraction,
               v$deaths_state$dated$walk$logit_fraction)
  out <- predict_reports(v, draw, 41:45)
  expect_equal(sort(unique(out$time_reported)), 41:45)
  expect_true(all(out$deaths_reported <= out$cases_reported))
})
