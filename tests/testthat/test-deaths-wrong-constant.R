# S6 of the sensitivity analyses: wrong records arrive at a constant rate
# eta_bg for each day recorded for and each calendar day from the day after
# it, whatever its cases (deaths-wrong-rate.R). The expected values are
# worked out here day by day from the model, not from the package's code.

# chi_j for each report of a day recorded for, day by day: one wrong record
# entered on each day k = t_{j-1} + 1, ..., t_j, still present at t_j with
# probability (1 - rho)^(t_j - k), with t_0 the day recorded for
constant_units_daily <- function(days, day_onset, rho) {
  before <- c(day_onset, days[-length(days)])
  vapply(seq_along(days), function(j) {
    if(days[j] <= before[j])
      return(0)
    k <- (before[j] + 1):days[j]
    sum((1 - rho)^(days[j] - k))
  }, numeric(1))
}

# reports with wrong records at a constant rate and a correcting report,
# with the state of the model of deaths after the paths are drawn
constant_example <- function(seed=1, n_days=12, eta_bg=0.1, rho=0.3,
                             max_erroneous=6) {
  set.seed(seed)
  cases_true <- stats::rpois(n_days, 6)
  sim <- simulate_cases_deaths(
    cases_true, 3:(n_days + 6), list(mean=3, sd=2), 0.3, list(mean=5, sd=3),
    wrong=list(eta_bg=eta_bg, rho=rho, delta=0.4),
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
  revision <- list(background_rate=eta_bg, rho=rho)
  state <- deaths_paths_step(state, largest$cases_true, reporting, revision,
                             max_erroneous)
  state <- deaths_uncorrected_step(state, reporting, revision, max_erroneous,
                                   max_cases=60)
  list(state=state, cases_true=largest$cases_true, reporting=reporting,
       revision=revision, settings=settings, max_erroneous=max_erroneous)
}

test_that("C1: the setting of the rate of wrong records is checked", {
  expect_equal(wrong_rate_settings(NULL)$form, "cases")
  expect_equal(wrong_rate_settings(list(eta_mu=0.1, eta_sigma=0.1))$form, "cases")
  expect_equal(wrong_rate_settings(list(wrong_rate="cases"))$form, "cases")
  settings <- wrong_rate_settings(list(wrong_rate="constant", eta_bg_mu=0.02,
                                       eta_bg_sigma=0.02))
  expect_equal(settings, list(form="constant", eta_bg_mu=0.02, eta_bg_sigma=0.02))
  expect_error(wrong_rate_settings(list(wrong_rate="background")), "must be one of")
  expect_error(wrong_rate_settings(list(wrong_rate=c("constant", "cases"))),
               "must be one of")
  expect_error(wrong_rate_settings(list(wrong_rate="constant")), "eta_bg_mu")
  expect_error(wrong_rate_settings(list(wrong_rate="constant", eta_bg_mu=0.02,
                                        eta_bg_sigma=-1)), "positive")
  expect_error(wrong_rate_settings(list(wrong_rate="constant", eta_bg_mu=0.02,
                                        eta_bg_sigma=0.02), is_deaths=FALSE),
               "model of deaths")
  # the starting values
  revision <- initial_wrong_rate(list(eta=0.05, rho=0.3), settings)
  expect_equal(revision, list(rho=0.3, background_rate=0.02))
  expect_equal(wrong_rate_form(revision), "constant")
  expect_equal(wrong_rate_form(list(eta=0.05, rho=0.3)), "cases")
  expect_true(is_wrong_records(revision))
  revision <- initial_wrong_rate(list(eta_bg=0.04, rho=0.3), settings)
  expect_equal(revision, list(rho=0.3, background_rate=0.04))
  # nothing is taken for eta by partial matching
  expect_null(revision$eta)
  expect_error(initial_wrong_rate(list(eta_bg=-0.04, rho=0.3), settings), "positive")
  expect_error(initial_wrong_rate(list(eta_bg=c(0.1, 0.2), rho=0.3), settings),
               "single positive")
  expect_error(initial_wrong_rate(list(eta_bg=0.04, rho=1.2), settings), "rho")
})

test_that("S6.1: the units of entry are chi_j, day by day and in closed form", {
  rho <- 0.3
  f_closed <- function(days, day_onset)
    (1 - (1 - rho)^diff(c(day_onset, days))) / rho
  cases <- list(
    consecutive=list(days=5:8, onset=4),
    gaps=list(days=c(6, 9, 13), onset=4),
    # left out of the reports of days 8 and 9: one interval spans them
    left_out=list(days=c(6, 7, 10, 11), onset=4),
    # the first report is made on the day recorded for
    onset_day=list(days=c(4, 5, 9), onset=4))
  for(name in names(cases)) {
    days <- cases[[name]]$days
    onset <- cases[[name]]$onset
    units <- wrong_constant_units(days, onset, rho)
    expect_equal(units$unit, constant_units_daily(days, onset, rho), info=name)
    expect_equal(units$unit, f_closed(days, onset), info=name)
    # the probability of withdrawal is that of the main model
    expect_equal(units$w, wrong_entry_units(days, onset, list(mean=5, sd=3),
                                            list(eta=0.1, rho=rho))$w, info=name)
    expect_equal(deaths_wrong_units(days, onset, list(mean=5, sd=3),
                                    list(background_rate=0.1, rho=rho)),
                 units, info=name)
  }
  expect_equal(wrong_constant_units(c(4, 5), 4, rho)$unit[1], 0)
  # the expected number present at each report followed from report to
  # report is that of entering and withdrawing them day by day
  days <- c(6, 9, 13, 14)
  units <- wrong_constant_units(days, 4, rho)
  present <- Reduce(function(x, j) x * (1 - units$w[j]) + units$unit[j],
                    seq_along(days), accumulate=TRUE, 0)[-1]
  expect_equal(present, vapply(days, function(d) sum((1 - rho)^(d - (5:d))), numeric(1)))
})

test_that("S6.2: with a constant rate the wrong records entered depend neither on the cases nor on eta", {
  example <- constant_example()
  state <- example$state
  revision <- example$revision
  for(k in seq_along(state$onset_times)) {
    days <- state$reports[[k]]$time_reported
    expect_equal(state$wrong[[k]]$lambda,
                 0.1 * constant_units_daily(days, state$onset_times[k], 0.3))
  }
  # the uncorrected series too
  sub <- state$uncorrected
  for(k in seq_along(sub$onset_times))
    expect_equal(sub$wrong[[k]]$lambda,
                 0.1 * constant_units_daily(sub$reports[[k]]$time_reported,
                                            sub$onset_times[k], 0.3))
  set.seed(2)
  more <- deaths_paths_step(state, example$cases_true + 7, example$reporting,
                            revision, example$max_erroneous)
  other <- deaths_paths_step(state, example$cases_true,
                             dplyr::tibble(mean=rep(6, 12), sd=4),
                             c(revision, list(eta=0.5)), example$max_erroneous)
  expect_equal(purrr::map(more$wrong, "lambda"), purrr::map(state$wrong, "lambda"))
  expect_equal(purrr::map(other$wrong, "lambda"), purrr::map(state$wrong, "lambda"))
  expect_equal(wrong_rate_form(c(revision, list(eta=0.5))), "constant")
  # the term of the wrong records entered is left out of the update of the
  # true cases, since it does not depend on them
  expect_null(deaths_wrong_mean(state, 3, revision))
  expect_false(is.null(deaths_wrong_mean(state, 3, list(eta=0.1, rho=0.3))))
  # but is in the probability of the reports used for the reporting delay
  path <- state$paths[[3]]
  entered <- path$wrong_entered_alive + path$wrong_entered_dead
  f_logp <- function(eta_bg)
    deaths_reporting_day_logp(state, 3, example$cases_true, list(mean=3, sd=2),
                              list(background_rate=eta_bg, rho=0.3))
  lambda <- function(eta_bg) eta_bg * constant_units_daily(
    state$reports[[3]]$time_reported, state$onset_times[3], 0.3)
  expect_equal(f_logp(0.2) - f_logp(0.1),
               sum(stats::dpois(entered, lambda(0.2), log=TRUE) -
                     stats::dpois(entered, lambda(0.1), log=TRUE)))
})

test_that("S6.3: the update of the fatal cases uses the wrong records of a constant rate", {
  # C3: the log-probability the update of the true fatal cases uses, from
  # the forward filter with the inputs the update holds, differs between
  # numbers of fatal cases by as much as the probability of the reports
  # found by enumerating every path, with the entries of a constant rate
  snapshot <- dplyr::tibble(time_onset=1, time_reported=c(2, 3, 5),
                            cases_reported=c(2, 4, 3), deaths_reported=c(0, 1, 1))
  settings <- deaths_settings(1)
  state <- initial_deaths_state(snapshot, dplyr::tibble(time_onset=1, cases_true=5),
                                settings)
  revision <- list(background_rate=0.3, rho=0.4)
  set.seed(3)
  state <- deaths_paths_step(state, 5, dplyr::tibble(mean=3, sd=2), revision,
                             max_erroneous=2)
  wrong <- state$wrong[[1]]
  lambda <- 0.3 * constant_units_daily(c(2, 3, 5), 1, 0.4)
  expect_equal(wrong$lambda, lambda)
  probabilities <- path_probabilities(state$survivals[[1]], 0)
  f_update <- function(fatal)
    joint_path_forward(snapshot$cases_reported, snapshot$deaths_reported, fatal,
                       5 - fatal, probabilities, wrong)$logp +
    stats::dbinom(fatal, 5, state$cfr[1], log=TRUE)
  f_scratch <- function(fatal)
    wrong_paths_reference(snapshot$cases_reported, snapshot$deaths_reported, fatal,
                          5 - fatal, probabilities,
                          list(lambda=lambda, w=1 - (1 - 0.4)^c(0, 1, 2),
                               delta=state$delta, max_erroneous=2))$logp +
    stats::dbinom(fatal, 5, state$cfr[1], log=TRUE)
  values <- 1:4
  update <- vapply(values, f_update, numeric(1))
  scratch <- vapply(values, f_scratch, numeric(1))
  expect_true(all(is.finite(update)))
  expect_equal(update - update[1], scratch - scratch[1], tolerance=1e-8)
})

test_that("S6.3: the update of the true and fatal cases together drops only the entries of wrong records", {
  # C3: the log-probability of the reports for possible true cases, with the
  # fatal cases summed out, differs between true cases by as much as the
  # full log-probability given the paths, written out term by term with the
  # entries of a constant rate
  example <- constant_example(seed=4)
  state <- example$state
  k <- 9
  path <- state$paths[[k]]
  survivals <- state$survivals[[k]]
  probabilities <- path_probabilities(survivals, 0)
  days <- state$reports[[k]]$time_reported
  lambda <- 0.1 * constant_units_daily(days, state$onset_times[k], 0.3)
  w <- 1 - (1 - 0.3)^c(0, diff(days))
  entered <- path$wrong_entered_alive + path$wrong_entered_dead
  cfr <- state$cfr[k]
  f_full <- function(cases, fatal) {
    survivors <- cases - fatal
    sum(dbinom_log_recycled(path$fatal_new, fatal - path$fatal_before,
                            probabilities$p_fatal)) +
      sum(dbinom_log_recycled(path$survivors_new, survivors - path$survivors_before,
                              probabilities$p)) +
      stats::dbinom(fatal, cases, cfr, log=TRUE) +
      sum(stats::dpois(entered, lambda, log=TRUE)) +
      sum(dbinom_log_recycled(path$wrong_removed, path$wrong_before, w))
  }
  f_scratch <- function(cases) {
    logp <- vapply(0:cases, function(fatal) f_full(cases, fatal), numeric(1))
    max(logp) + log(sum(exp(logp - max(logp))))
  }
  f_update <- deaths_cases_logp_function(state, state$onset_times[k], example$revision)
  values <- attr(f_update, "min_cases") + 0:8
  update <- f_update(values)
  scratch <- vapply(values, f_scratch, numeric(1))
  expect_true(all(is.finite(update)))
  expect_equal(update - update[1], scratch - scratch[1], tolerance=1e-10)
})

test_that("S6.4: the constant rate is drawn from its gamma conditional distribution, and rho's target is the direct sum", {
  example <- constant_example(seed=5)
  state <- example$state
  priors <- list(eta_bg_mu=0.05, eta_bg_sigma=0.05, rho_alpha=2, rho_beta=3)
  # the wrong records entered and the sum of chi_j, counted from the paths
  # of every series, including the uncorrected ones
  all_paths <- c(state$paths, state$uncorrected$paths)
  expect_true(length(state$uncorrected$paths) > 0)
  n <- sum(vapply(all_paths, function(p) sum(p$wrong_entered_alive) +
                    sum(p$wrong_entered_dead), numeric(1)))
  reports <- c(state$reports, state$uncorrected$reports)
  onsets <- c(state$onset_times, state$uncorrected$onset_times)
  exposure <- sum(vapply(seq_along(reports), function(k)
    sum(constant_units_daily(reports[[k]]$time_reported, onsets[k], 0.3)), numeric(1)))
  expect_true(n > 0)
  shape <- priors$eta_bg_mu^2 / priors$eta_bg_sigma^2 + n
  rate <- priors$eta_bg_mu / priors$eta_bg_sigma^2 + exposure
  set.seed(6)
  n_draws <- 5000
  draws <- vapply(seq_len(n_draws), function(i) {
    updated <- deaths_wrong_step(state, example$cases_true, example$reporting,
                                 example$revision, example$settings, priors)
    updated$revision_parameters$background_rate
  }, numeric(1))
  expect_z(mean(draws), shape / rate, sqrt(shape) / rate / sqrt(n_draws),
           "mean of the constant rate")
  expect_z(stats::var(draws), shape / rate^2,
           shape / rate^2 * sqrt((2 + 6 / shape) / n_draws), "variance of the constant rate")
  expect_ks(draws, function(x) stats::pgamma(x, shape, rate), "constant rate")

  # the target of the update of rho is the direct sum over every series
  series <- deaths_wrong_series(state, example$cases_true, example$reporting)
  for(rho in c(0.1, 0.3, 0.7)) {
    direct <- sum(vapply(seq_along(all_paths), function(k) {
      p <- all_paths[[k]]
      days <- reports[[k]]$time_reported
      sum(dbinom_log_recycled(p$wrong_removed, p$wrong_before,
                              1 - (1 - rho)^c(0, diff(days)))) +
        sum(stats::dpois(p$wrong_entered_alive + p$wrong_entered_dead,
                         0.1 * constant_units_daily(days, onsets[k], rho), log=TRUE))
    }, numeric(1))) +
      stats::dbeta(rho, 2, 3, log=TRUE) + log(rho) + log(1 - rho)
    expect_equal(deaths_wrong_rho_logp(rho, series, example$revision, priors), direct,
                 tolerance=1e-10)
  }
})

# every configuration of wrong records which are present at one or more of
# the reports: those entered on day k and present at the reports from the
# first on or after k to report j1 are Poisson, independently, with mean
# eta_bg [(1 - rho)^(t_j1 - k) - (1 - rho)^(t_(j1 + 1) - k)], split between
# alive and dead; those present at no report sum out
constant_enumeration <- function(days, day_onset, eta_bg, rho, delta, max_erroneous) {
  n_reports <- length(days)
  cells <- list()
  for(k in (day_onset + 1):max(days)) {
    j0 <- which(days >= k)[1]
    for(j1 in j0:n_reports)
      for(is_dead in c(FALSE, TRUE)) {
        present <- (1 - rho)^(days[j1] - k) -
          if(j1 < n_reports) (1 - rho)^(days[j1 + 1] - k) else 0
        cells[[length(cells) + 1]] <- list(
          reports=j0:j1, is_dead=is_dead,
          mean=eta_bg * present * ifelse(is_dead, delta, 1 - delta))
      }
  }
  out <- list()
  f_recurse <- function(i, alive, dead, logp) {
    if(i > length(cells)) {
      out[[length(out) + 1]] <<- list(alive=alive, dead=dead, logp=logp)
      return(invisible())
    }
    cell <- cells[[i]]
    n <- 0
    repeat {
      add <- rep(0, n_reports)
      add[cell$reports] <- n
      a <- alive + if(cell$is_dead) 0 else add
      d <- dead + if(cell$is_dead) add else 0
      if(any(a + d > max_erroneous))
        break
      f_recurse(i + 1, a, d, logp + stats::dpois(n, cell$mean, log=TRUE))
      n <- n + 1
    }
  }
  f_recurse(1, rep(0, n_reports), rep(0, n_reports), 0)
  out
}

test_that("S6.5: the forward filter gives the probability of every sequence of reports found by entering and withdrawing wrong records day by day", {
  days <- c(1, 3, 4)
  day_onset <- 0
  eta_bg <- 0.3
  rho <- 0.35
  delta <- 0.4
  survivors <- 2
  p <- c(0.3, 0.5, 0.4)
  configurations <- constant_enumeration(days, day_onset, eta_bg, rho, delta, 2)
  # with the genuine cases reported in a chain of binomials
  genuine <- expand.grid(g1=0:2, g2=0:2, g3=0:2)
  genuine <- genuine[genuine$g1 <= genuine$g2 & genuine$g2 <= genuine$g3, ]
  probability <- list()
  for(r in seq_len(nrow(genuine))) {
    g <- as.numeric(genuine[r, ])
    p_genuine <- prod(stats::dbinom(diff(c(0, g)), survivors - c(0, g[-3]), p))
    for(x in configurations) {
      key <- paste(c(g + x$alive + x$dead, x$dead), collapse=" ")
      before <- if(is.null(probability[[key]])) 0 else probability[[key]]
      probability[[key]] <- before + p_genuine * exp(x$logp)
    }
  }
  revision <- list(background_rate=eta_bg, rho=rho)
  units <- deaths_wrong_units(days, day_onset, list(mean=3, sd=2), revision)
  wrong <- list(lambda=wrong_entry_mean(units, survivors, revision), w=units$w,
                delta=delta, max_erroneous=2)
  probabilities <- list(p=p, p_fatal=p, g=c(0.5, 0.5, 0.5), q=c(0.5, 0.5, 0.5))
  expect_true(length(probability) > 50)
  for(key in names(probability)) {
    counts <- as.numeric(strsplit(key, " ")[[1]])
    filtered <- joint_path_forward(counts[1:3], counts[4:6], 0, survivors,
                                   probabilities, wrong)$logp
    expect_equal(filtered, log(probability[[key]]), tolerance=1e-8, info=key)
  }
  # what is left out is the probability of more wrong records than are allowed
  # for at some report (S6.6)
  expect_equal(1 - sum(unlist(probability)),
               wrong_cap_probability(wrong$lambda, wrong$w, 2), tolerance=1e-10)
})

test_that("S6.6: the probability of more wrong records than are allowed for agrees with simulation", {
  set.seed(7)
  day_onset <- 0
  days <- seq(2, 20, 2)
  eta_bg <- 0.5
  rho <- 0.3
  max_erroneous <- 3
  units <- wrong_constant_units(days, day_onset, rho)
  exact <- wrong_cap_probability(eta_bg * units$unit, units$w, max_erroneous)
  expect_true(exact > 0.05 & exact < 0.95)
  n_sims <- 20000
  # day by day: entered on each day after the day recorded for, each
  # surviving each later day with probability 1 - rho
  is_over <- vapply(seq_len(n_sims), function(i) {
    k <- (day_onset + 1):max(days)
    entered <- rep(k, stats::rpois(length(k), eta_bg))
    life <- stats::rgeom(length(entered), rho)
    present <- vapply(days, function(d) sum(entered <= d & entered + life >= d), numeric(1))
    any(present > max_erroneous)
  }, logical(1))
  expect_prop_z(sum(is_over), n_sims, exact, "probability of exceeding the limit")
  # none can be exceeded without wrong records, and a larger limit is
  # exceeded less often
  expect_equal(wrong_cap_probability(rep(0, 3), c(0, 0.5, 0.5), 2), 0)
  expect_lt(wrong_cap_probability(eta_bg * units$unit, units$w, 6), exact)

  # for every series of a state, and at the posterior mean with a warning
  example <- constant_example(seed=8)
  by_day <- deaths_wrong_cap(example$state, example$cases_true, example$reporting,
                             example$revision, 2)
  expect_equal(nrow(by_day), length(example$state$onset_times) +
                 length(example$state$uncorrected$onset_times))
  expect_equal(sum(by_day$is_uncorrected), length(example$state$uncorrected$onset_times))
  k <- 5
  units <- wrong_constant_units(example$state$reports[[k]]$time_reported,
                                example$state$onset_times[k], 0.3)
  expect_equal(by_day$probability[k], wrong_cap_probability(0.1 * units$unit, units$w, 2))
  variables <- list(revision_current=example$revision, is_deaths=TRUE,
                    deaths_state=example$state, deaths_cases_true=example$cases_true,
                    deaths_running_df=example$reporting, max_erroneous=2)
  revisions <- dplyr::tibble(eta=NA, rho=c(0.2, 0.4), pi=0, iteration=1:2,
                             eta_bg=c(0.05, 0.15))
  expect_warning(cap <- wrong_records_cap(variables, revisions), "max_erroneous")
  expect_equal(cap$by_day, by_day)
  expect_equal(cap$largest, max(by_day$probability))
  variables$max_erroneous <- 12
  expect_silent(cap <- wrong_records_cap(variables, revisions))
  expect_lt(cap$largest, 0.01)
  variables$revision_current <- list(eta=0.1, rho=0.3)
  expect_error(wrong_records_cap(variables, revisions), "other than 'cases'")
})

test_that("C6: simulate_cases_deaths enters wrong records at a constant rate", {
  set.seed(9)
  days <- c(2, 3, 6, 10)
  eta_bg <- 0.4
  rho <- 0.25
  delta <- 0.3
  n_sims <- 2000
  counts <- replicate(n_sims, {
    sim <- simulate_cases_deaths(c(3, 30), days, list(mean=3, sd=2), 0,
                                 list(mean=5, sd=3),
                                 wrong=list(eta_bg=eta_bg, rho=rho, delta=delta))
    rbind(sim$snapshots$cases_reported, sim$snapshots$deaths_reported)
  })
  # the genuine cases reported, and the wrong records present, entered
  # day by day and withdrawn day by day
  f_expected <- function(t, cases) {
    d <- days[days >= t]
    genuine <- cases * stats::pgamma(d - t + 0.5, shape=(3 / 2)^2, rate=3 / 2^2)
    wrong <- vapply(d, function(r) if(r <= t) 0 else
      sum(eta_bg * (1 - rho)^(r - ((t + 1):r))), numeric(1))
    list(cases=genuine + wrong, deaths=delta * wrong)
  }
  expected <- list(cases=c(f_expected(1, 3)$cases, f_expected(2, 30)$cases),
                   deaths=c(f_expected(1, 3)$deaths, f_expected(2, 30)$deaths))
  zmax <- z_bonferroni(2 * length(expected$cases))
  for(i in seq_along(expected$cases)) {
    expect_z(mean(counts[1, i, ]), expected$cases[i], stats::sd(counts[1, i, ]) / sqrt(n_sims),
             paste("mean count of cases, row", i), zmax=zmax)
    if(expected$deaths[i] > 0)
      expect_z(mean(counts[2, i, ]), expected$deaths[i],
               stats::sd(counts[2, i, ]) / sqrt(n_sims),
               paste("mean count of deaths, row", i), zmax=zmax)
  }
  # the count of the onset day with few cases has as many wrong records as
  # that with many
  expect_true(all(counts[2, , ] <= counts[1, , ]))
})

# a short fit with wrong records at a constant rate
constant_fit <- function(example, niterations, ...) {
  snapshots <- example$snapshots
  mcmc_single(
    niterations=niterations,
    snapshot_with_Rt_index_df=snapshots %>%
      dplyr::select(time_onset, time_reported, cases_reported) %>%
      dplyr::mutate(Rt_index=1),
    priors=list(Rt=list(shape=1, rate=1),
                reporting=list(mean_mu=6, mean_sigma=5, sd_mu=4, sd_sigma=4),
                max_cases=300, overdispersion=list(mean=20, sd=20),
                revisions=list(wrong_rate="constant", eta_bg_mu=0.05,
                               eta_bg_sigma=0.05, rho_alpha=1, rho_beta=1),
                max_erroneous=8,
                Rt_changepoints=list(),
                reporting_changepoints=list(min_length=14)),
    serial_parameters=example$serial_parameters,
    initial_cases_true=example$initial_cases,
    initial_reporting_parameters=list(mean=6, sd=4),
    initial_Rt=dplyr::tibble(Rt_index=1, Rt=1),
    initial_revision_parameters=list(eta=0.02, rho=0.3, weekday=TRUE),
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

constant_fit_example <- function() {
  set.seed(10)
  n_days <- 20
  serial_parameters <- list(mean=5, sd=3)
  cases <- true_cases(n_days + 40, function(t) 1, kappa=1000, serial_parameters,
                      initial_parameters=list(mean=20, length=10))
  sim <- simulate_cases_deaths(cases[41:(n_days + 40)], seq(2, n_days + 20, 2),
                               list(mean=6, sd=4), 0.5, list(mean=3, sd=2.5),
                               triggered=1, onset_to_death=list(mean=9.3, sd=5.5),
                               dated_by_death=list(fraction=0.3),
                               wrong=list(eta_bg=0.05, rho=0.3, delta=0.3))
  snapshots <- dplyr::filter(sim$snapshots, time_reported <= n_days + 10)
  list(snapshots=snapshots, cases_before=cases[1:40],
       serial_parameters=serial_parameters,
       initial_cases=snapshots %>%
         dplyr::group_by(time_onset) %>%
         dplyr::summarise(cases_true=max(cases_reported) + 2))
}

test_that("C7 and C8: a fit with a constant rate runs, is taken up again and predicts later reports", {
  example <- constant_fit_example()
  set.seed(11)
  whole <- constant_fit(example, 6)
  revisions <- whole$revisions
  expect_true(all(is.na(revisions$eta)))
  expect_true(all(revisions$eta_bg > 0))
  expect_equal(length(unique(revisions$eta_bg)), 6)
  expect_true(all(revisions$rho > 0 & revisions$rho < 1))
  expect_true(all(c("largest", "expected", "by_day") %in% names(whole$wrong_cap)))
  expect_equal(whole$wrong_cap$means[["eta_bg"]], mean(revisions$eta_bg))
  # C7: taken up again, the draws are those of the fit run straight through
  set.seed(11)
  first <- constant_fit(example, 3, return_state=TRUE)
  set.seed(99)
  second <- constant_fit(example, 3, resume_from=first$state)
  names_draws <- setdiff(names(whole), c("state", "iterations_done", "wrong_cap"))
  for(name in names_draws)
    expect_equal(dplyr::bind_rows(first[[name]], second[[name]]), whole[[name]],
                 info=name)
  # C8: a kept draw gives the rate, which the predictions use after the
  # last report
  draw <- lapply(second[c("cases", "reporting_days", "revisions", "reporting_effects",
                          "deaths", "deaths_parameters")],
                 function(df) df[df$iteration == 6, ])
  variables <- second$state$variables
  prediction <- predict_state_from_draw(variables, draw)
  expect_equal(prediction$revision_parameters$background_rate, draw$revisions$eta_bg)
  expect_null(prediction$revision_parameters$eta)
  days <- c(example$snapshots$time_reported[example$snapshots$time_onset == 5],
            34, 37)
  terms <- predict_onset_terms(prediction, 5, days)
  expect_equal(terms$wrong$lambda, draw$revisions$eta_bg *
                 constant_units_daily(days, 5, draw$revisions$rho))
  set.seed(12)
  later <- predict_reports(variables, draw, c(34, 37), onset_times=c(5, 18))
  expect_equal(nrow(later), 4)
  expect_true(all(later$deaths_reported <= later$cases_reported))
  # the cap is reported for the fit's draws, here combined
  cap <- wrong_records_cap(variables, dplyr::bind_rows(first$revisions,
                                                       second$revisions))
  expect_equal(cap$means[["rho"]], mean(whole$revisions$rho))
  # a constant rate needs a model of deaths
  expect_error(mcmc_single(
    niterations=1,
    snapshot_with_Rt_index_df=example$snapshots %>%
      dplyr::select(time_onset, time_reported, cases_reported) %>%
      dplyr::mutate(Rt_index=1),
    priors=list(Rt=list(shape=1, rate=1),
                reporting=list(mean_mu=6, mean_sigma=5, sd_mu=4, sd_sigma=4),
                max_cases=300,
                revisions=list(wrong_rate="constant", eta_bg_mu=0.05, eta_bg_sigma=0.05,
                               eta_mu=0.05, eta_sigma=0.05, rho_alpha=1, rho_beta=1)),
    serial_parameters=example$serial_parameters,
    initial_cases_true=example$initial_cases,
    initial_reporting_parameters=list(mean=6, sd=4),
    initial_Rt=dplyr::tibble(Rt_index=1, Rt=1),
    initial_revision_parameters=list(eta=0.02, rho=0.3),
    print_to_screen=FALSE), "model of deaths")
})
