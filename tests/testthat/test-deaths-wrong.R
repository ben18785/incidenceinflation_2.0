test_that("wrong_entry_units are those of the model of wrong records", {
  reporting <- list(mean=5, sd=3)
  revision <- list(eta=0.1, rho=0.3)
  days <- c(6, 7, 9)
  units <- wrong_entry_units(days, 4, reporting, revision)
  undetected <- function(day) exp(log_undetected_days(day, 4, reporting))
  expect_equal(units$w, c(0, 0.3, 1 - 0.7^2))
  # entered on each day since the report before and still present at the report
  expect_equal(units$unit[2:3],
               c(undetected(6) - undetected(7),
                 (undetected(7) - undetected(8)) * 0.7 + (undetected(8) - undetected(9))))
  # entered on each day before the first report and still present at it
  expect_equal(units$unit[1],
               (undetected(4) - undetected(5)) * 0.7 +
                 (undetected(5) - undetected(6)))
  # the first report is made on the onset day itself
  expect_equal(wrong_entry_units(c(4, 5), 4, reporting, revision)$unit[1], 0)
})

test_that("joint_path_forward agrees with the sampler written in R when there are no wrong records", {
  cases <- c(2, 3, 5, 5)
  deaths <- c(0, 1, 2, 3)
  probabilities <- list(p=c(0.3, 0.4, 0.5, 0.2), p_fatal=c(0.5, 0.2, 0.7, 0.3),
                        g=c(0.2, 0.4, 0.6, 0.7), q=c(0.2, 0.25, 0.33, 0.4))
  for(fatal in 3:5)
    for(survivors in 2:4)
      expect_equal(
        joint_path_forward(cases, deaths, fatal, survivors, probabilities)$logp,
        fatal_path_forward(cases, deaths, fatal, survivors, probabilities)$logp)
  expect_equal(joint_path_forward(cases, deaths, 2, 5, probabilities)$logp, -Inf)
  expect_equal(joint_path_forward(cases, deaths, -1, 5, probabilities)$logp, -Inf)
  # without wrong records a count cannot fall
  expect_equal(joint_path_forward(c(3, 2), c(0, 0), 3, 3, probabilities)$logp, -Inf)
  expect_error(joint_path_forward(c(1, 2), c(2, 2), 3, 3, probabilities))
  set.seed(21)
  forward <- joint_path_forward(cases, deaths, 4, 3, probabilities)
  path <- joint_path_backward(forward)
  expect_equal(path$fatal_new + path$survivors_new, diff(c(0, cases)))
  expect_equal(path$deaths_new, diff(c(0, deaths)))
  expect_equal(path$wrong_before, rep(0, 4))
  expect_true(all(path$alive_before >= 0))
})

test_that("joint_path_forward gives the probability of reports in which counts fall", {
  # the case count falls in the third report, and the death count in the fourth
  cases <- c(3, 5, 4, 5)
  deaths <- c(1, 2, 2, 1)
  probabilities <- list(p=c(0.3, 0.4, 0.5, 0.2), p_fatal=c(0.5, 0.2, 0.7, 0.3),
                        g=c(0.2, 0.4, 0.6, 0.7), q=c(0.2, 0.25, 0.33, 0.4))
  wrong <- list(lambda=c(0.3, 0.2, 0.1, 0.4), w=c(0, 0.3, 0.5, 0.3), delta=0.4,
                max_erroneous=2)
  for(fatal in c(2, 4))
    for(survivors in c(3, 5)) {
      reference <- wrong_paths_reference(cases, deaths, fatal, survivors,
                                         probabilities, wrong)
      expect_equal(
        joint_path_forward(cases, deaths, fatal, survivors, probabilities, wrong)$logp,
        reference$logp, tolerance=1e-8)
    }
  # impossible unless enough wrong records are allowed
  wrong$max_erroneous <- 0
  expect_equal(joint_path_forward(cases, deaths, 4, 5, probabilities, wrong)$logp, -Inf)
})

test_that("joint_path_backward draws paths with the right probabilities", {
  set.seed(22)
  cases <- c(3, 2, 3)
  deaths <- c(1, 1, 1)
  probabilities <- list(p=c(0.4, 0.4, 0.5), p_fatal=c(0.5, 0.3, 0.6),
                        g=c(0.3, 0.4, 0.6), q=c(0.2, 0.25, 0.33))
  wrong <- list(lambda=c(0.4, 0.3, 0.3), w=c(0, 0.4, 0.5), delta=0.4,
                max_erroneous=2)
  reference <- wrong_paths_reference(cases, deaths, 3, 3, probabilities, wrong)
  states <- reference$states
  n_states <- nrow(states)
  # the probability of every path of states
  paths <- expand.grid(seq_len(n_states), seq_len(n_states), seq_len(n_states))
  weights <- reference$transitions[[1]][reference$start, paths[, 1]] *
    reference$transitions[[2]][cbind(paths[, 1], paths[, 2])] *
    reference$transitions[[3]][cbind(paths[, 2], paths[, 3])]
  f_label <- function(r, alive, dead) paste(r, alive, dead, sep="/", collapse=" ")
  labels <- apply(paths, 1, function(index)
    f_label(states$r[index], states$alive[index], states$dead[index]))
  is_possible <- weights > 0
  expected <- stats::setNames(weights[is_possible] / sum(weights), labels[is_possible])

  forward <- joint_path_forward(cases, deaths, 3, 3, probabilities, wrong)
  n_draws <- 6000
  draws <- replicate(n_draws, {
    path <- joint_path_backward(forward)
    fatal <- path$fatal_before + path$fatal_new
    wrong_dead <- deaths - (cumsum(path$deaths_new))
    wrong_total <- cases - fatal - (path$survivors_before + path$survivors_new)
    # those entered and removed account for the wrong records present
    expect_equal(wrong_total,
                 path$wrong_before - path$wrong_removed +
                   path$wrong_entered_alive + path$wrong_entered_dead)
    expect_true(all(path$wrong_removed <= path$wrong_before))
    f_label(fatal, wrong_total - wrong_dead, wrong_dead)
  })
  expect_true(all(draws %in% names(expected)))
  observed <- table(factor(draws, levels=names(expected)))
  # paths too rare to be expected in the draws are pooled
  is_common <- expected * n_draws >= 5
  counts <- c(observed[is_common], sum(observed[!is_common]))
  probabilities_pooled <- c(expected[is_common], sum(expected[!is_common]))
  expect_true(length(counts) > 5)
  expect_true(stats::chisq.test(counts, p=probabilities_pooled)$p.value > 0.001)
})

test_that("sample_spike_slab draws from the prior when nothing depends on the value", {
  set.seed(23)
  value <- 0
  draws <- numeric(6000)
  for(i in seq_along(draws)) {
    value <- sample_spike_slab(value, function(x) 0, probability_zero=0.3, sd=0.5,
                               step=0.3, n_steps=3)
    draws[i] <- value
  }
  expect_true(abs(mean(draws == 0) - 0.3) < 0.04)
  expect_true(abs(stats::sd(draws[draws != 0]) - 0.5) < 0.05)
  expect_true(abs(mean(draws[draws != 0])) < 0.05)
  # and moves to where the data put it
  value <- 0
  for(i in 1:200)
    value <- sample_spike_slab(value, function(x) stats::dnorm(x, -0.8, 0.05, log=TRUE))
  expect_true(abs(value + 0.8) < 0.2)
  for(i in 1:50)
    value <- sample_spike_slab(value, function(x) stats::dnorm(x, 0, 0.01, log=TRUE))
  expect_equal(value, 0)
})

test_that("metropolis_weekday_effects keeps effects which sum to zero", {
  set.seed(24)
  effects <- stats::setNames(rep(0, 5), 1:5)
  target <- c(0.4, -0.2, 0.1, -0.4, 0.1)
  f_logl <- function(x) sum(stats::dnorm(x, target, 0.05, log=TRUE))
  current <- list(log_effects=effects, logl=f_logl(effects))
  for(i in 1:600)
    current <- metropolis_weekday_effects(current$log_effects, f_logl,
                                          current$logl, 1, 0.05)
  expect_equal(sum(current$log_effects), 0)
  expect_equal(names(current$log_effects), as.character(1:5))
  expect_true(all(abs(current$log_effects - target) < 0.15))
  expect_equal(current$logl, f_logl(current$log_effects))
  single <- metropolis_weekday_effects(c("1"=0), f_logl, 0, 1, 0.05)
  expect_equal(single$log_effects, c("1"=0))
})

test_that("simulate_cases_deaths adds wrong records which are later removed", {
  set.seed(25)
  cases <- rep(100, 30)
  sim <- simulate_cases_deaths(cases, 1:90, list(mean=5, sd=3), 0.3,
                               list(mean=10, sd=5),
                               wrong=list(eta=0.05, rho=0.2, delta=0.5))
  falls <- sim$snapshots %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(cases=sum(diff(cases_reported) < 0),
                     deaths=sum(diff(deaths_reported) < 0),
                     cases_last=dplyr::last(cases_reported),
                     deaths_last=dplyr::last(deaths_reported))
  # a fall is seen only when more wrong records are removed than cases added
  expect_true(sum(falls$cases) > 5)
  expect_true(sum(falls$deaths) > 5)
  # all wrong records have been removed by the end
  expect_equal(falls$cases_last, sim$truth$cases_true)
  expect_equal(falls$deaths_last, sim$truth$deaths_true)
  expect_true(all(sim$snapshots$deaths_reported <= sim$snapshots$cases_reported))
})

test_that("mcmc_single runs with wrong records, effects of the day of the week and a delay for fatal cases", {
  set.seed(26)
  n_days <- 35
  serial_parameters <- list(mean=5, sd=3)
  cases <- true_cases(n_days + 40, function(t) 1, kappa=1000, serial_parameters,
                      initial_parameters=list(mean=30, length=10))
  cases_before <- cases[1:40]
  cases <- cases[41:(n_days + 40)]
  sim <- simulate_cases_deaths(cases, seq_len(n_days + 30), list(mean=4, sd=3), 0.5,
                               list(mean=8, sd=4),
                               fatal_reporting_parameters=list(mean=2, sd=1.5),
                               wrong=list(eta=0.05, rho=0.3, delta=0.5))
  snapshots <- sim$snapshots
  initial_cases <- snapshots %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(cases_true=max(cases_reported))
  largest_fall <- max_downward_revision(snapshots)
  expect_true(largest_fall >= 1)
  priors <- list(Rt=list(shape=1, rate=1),
                 reporting=list(mean_mu=5, mean_sigma=5, sd_mu=3, sd_sigma=3),
                 max_cases=300,
                 revisions=list(eta_mu=0.05, eta_sigma=0.05, rho_alpha=2, rho_beta=4),
                 max_erroneous=largest_fall + 1,
                 Rt_changepoints=list(),
                 reporting_changepoints=list(min_length=14))
  f_fit <- function(initial_revisions) mcmc_single(
    niterations=25,
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
    initial_revision_parameters=initial_revisions,
    deaths=list(snapshot=snapshots %>%
                  dplyr::select(time_onset, time_reported, deaths_reported),
                settings=list(initial_delay=list(mean=8, sd=4))))
  fit <- f_fit(list(eta=0.05, rho=0.3, weekday=TRUE, weekday_first=3))
  expect_equal(nrow(fit$revisions), 25)
  expect_true(all(fit$revisions$eta > 0 & fit$revisions$rho > 0 & fit$revisions$rho < 1))
  expect_true(length(unique(fit$revisions$eta)) > 20)
  expect_true(all(paste0("weekday_", 1:7) %in% colnames(fit$reporting_effects)))
  parameters <- fit$deaths_parameters
  expect_true(all(c("fatal_delay_ratio", "wrong_dead",
                    paste0("recording_weekday_", 1:7)) %in% colnames(parameters)))
  expect_true(all(parameters$wrong_dead > 0 & parameters$wrong_dead < 1))
  # effects of the day of the week multiply to one
  expect_equal(apply(as.matrix(parameters[, paste0("recording_weekday_", 1:7)]), 1, prod),
               rep(1, 25))
  # fatal cases are reported sooner
  expect_true(all(parameters$fatal_delay_ratio[15:25] < 0.9))
  both <- dplyr::left_join(fit$deaths, fit$cases, by=c("iteration", "time_onset"))
  expect_true(all(both$deaths_true <= both$cases_true))
  # what is not allowed for alongside deaths
  expect_error(f_fit(list(eta=0.05, rho=0.3, pi=0.1)))
})

test_that("wrong records are withdrawn daily, between reports as before the first", {
  # the expected number present at each report, followed from report to
  # report, equals that of entering and withdrawing them day by day
  reporting <- list(mean=8, sd=6)
  revision <- list(eta=0.2, rho=0.1)
  day_onset <- 3
  days <- day_onset + c(5, 13, 27)
  units <- wrong_entry_units(days, day_onset, reporting, revision)
  present <- numeric(length(days))
  for(j in seq_along(days))
    present[j] <- (if(j > 1) present[j - 1] * (1 - units$w[j]) else 0) + units$unit[j]
  undetected <- exp(log_undetected_days(day_onset:max(days), day_onset, reporting))
  newly <- -diff(undetected)
  daily <- sapply(days, function(d) {
    k <- seq_len(d - day_onset)
    sum(newly[k] * (1 - revision$rho)^(d - day_onset - k))
  })
  expect_equal(present, daily)
  # and the simulator of the reports of one onset day agrees: its mean count
  # is the genuine cases reported plus the wrong records present
  set.seed(4)
  n <- 1000
  counts <- replicate(1000, observed_cases_trajectory(n, days, day_onset, reporting,
                                                      revision))
  expected <- n * (1 - undetected[days - day_onset + 1]) + revision$eta * n * daily
  expect_equal(rowMeans(counts), expected, tolerance=0.01)
})
