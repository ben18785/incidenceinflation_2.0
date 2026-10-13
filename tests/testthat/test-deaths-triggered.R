# fatal cases which are reported when their death is recorded, if not before

test_that("no fatal cases reported at death gives the model without them", {
  days <- c(3, 4, 6, 9, 14, 20)
  reporting <- list(mean=7, sd=5)
  death <- list(mean=11, sd=9)
  expect_identical(joint_probabilities(days, 1, reporting, death, triggered=0),
                   joint_probabilities(days, 1, reporting, death))
  # only fatal cases are affected, and not how deaths are recorded later
  with <- joint_probabilities(days, 1, reporting, death, triggered=0.5)
  without <- joint_probabilities(days, 1, reporting, death)
  expect_equal(with$p, without$p)
  expect_equal(with$q, without$q)
  # fatal cases are reported sooner, and more of those in the first report
  # are dead. Later the fatal cases still to be reported are mostly those
  # whose deaths have not been recorded, so fewer can be dead on being reported.
  expect_true(all(with$p_fatal[1:3] > without$p_fatal[1:3]))
  expect_gt(with$g[1], without$g[1])
  expect_true(all(with$g >= 0 & with$g <= 1))
})

test_that("the probabilities with fatal cases reported at death match individual cases", {
  set.seed(4)
  days <- c(3, 4, 6, 9, 10, 14, 20, 30)
  reporting <- list(mean=7, sd=5)
  death <- list(mean=11, sd=9)
  triggered <- 0.6
  n <- 400000
  f_delay <- function(p) stats::rgamma(n, p$mean^2 / p$sd^2, p$mean / p$sd^2)
  delay_reported <- f_delay(reporting)
  delay_recorded <- f_delay(death)
  is_triggered <- stats::runif(n) < triggered
  delay_reported[is_triggered] <- pmin(delay_reported, delay_recorded)[is_triggered]
  # a delay d is in a report made at least d - 0.5 days after onset
  elapsed <- days - 1 + 0.5
  is_reported <- cbind(FALSE, sapply(elapsed, function(e) delay_reported <= e))
  is_dead <- cbind(FALSE, sapply(elapsed, function(e)
    pmax(delay_reported, delay_recorded) <= e))
  simulated <- purrr::map_dfr(seq_along(days), function(j) {
    is_new <- is_reported[, j + 1] & !is_reported[, j]
    is_alive <- is_reported[, j] & !is_dead[, j]
    dplyr::tibble(p_fatal=sum(is_new) / sum(!is_reported[, j]),
                  g=sum(is_new & is_dead[, j + 1]) / sum(is_new),
                  q=sum(is_alive & is_dead[, j + 1]) / sum(is_alive))
  })
  probabilities <- joint_probabilities(days, 1, reporting, death, triggered=triggered)
  expect_equal(simulated$p_fatal, probabilities$p_fatal, tolerance=0.02)
  expect_equal(simulated$g, probabilities$g, tolerance=0.02)
  expect_equal(simulated$q[-1], probabilities$q[-1], tolerance=0.02)
})

test_that("a case reported at its death is shown as dead at once in simulated reports", {
  set.seed(8)
  cases <- rep(40, 20)
  # the reporting delay is so long that only deaths bring cases into the reports
  sim <- simulate_cases_deaths(cases, 1:30, list(mean=500, sd=10), 0.5,
                               list(mean=6, sd=3), triggered=1)
  expect_gt(sum(sim$snapshots$cases_reported), 0)
  expect_equal(sim$snapshots$cases_reported, sim$snapshots$deaths_reported)
  expect_error(simulate_cases_deaths(cases, 1:30, list(mean=5, sd=3), 0.5,
                                     list(mean=6, sd=3), triggered=1.5),
               "between 0 and 1")
  expect_error(simulate_cases_deaths(cases, 1:30, list(mean=5, sd=3), 0.5,
                                     list(mean=6, sd=3), triggered=0.5,
                                     death_from="report"),
               "death_from")
})

test_that("the settings for fatal cases reported at death are checked", {
  expect_false(deaths_settings(30)$triggered)
  expect_equal(deaths_settings(30, list(triggered=list(alpha=2)))$triggered,
               list(alpha=2, beta=1, initial=0.5, step=0.3, n_steps=3))
  expect_error(deaths_settings(30, list(triggered=list(gamma=1))), "Unknown settings")
  expect_error(deaths_settings(30, list(triggered=list(alpha=0))), "must be positive")
  expect_error(deaths_settings(30, list(triggered=list(initial=1))), "strictly between")
})

triggered_example <- function() {
  set.seed(31)
  n_days <- 30
  serial_parameters <- list(mean=5, sd=3)
  cases <- true_cases(n_days + 40, function(t) 1, kappa=1000, serial_parameters,
                      initial_parameters=list(mean=30, length=10))
  sim <- simulate_cases_deaths(cases[41:(n_days + 40)], seq_len(n_days + 10),
                               list(mean=6, sd=4), 0.5, list(mean=8, sd=5),
                               triggered=0.6)
  list(snapshots=sim$snapshots, cases_before=cases[1:40],
       serial_parameters=serial_parameters,
       initial_cases=sim$snapshots %>%
         dplyr::group_by(time_onset) %>%
         dplyr::summarise(cases_true=max(cases_reported)))
}

triggered_fit <- function(example, niterations, triggered=list(), ...) {
  snapshots <- example$snapshots
  mcmc_single(
    niterations=niterations,
    snapshot_with_Rt_index_df=snapshots %>%
      dplyr::select(time_onset, time_reported, cases_reported) %>%
      dplyr::mutate(Rt_index=1),
    priors=list(Rt=list(shape=1, rate=1),
                reporting=list(mean_mu=6, mean_sigma=5, sd_mu=4, sd_sigma=4),
                max_cases=300,
                overdispersion=list(mean=20, sd=20),
                Rt_changepoints=list(),
                reporting_changepoints=list(min_length=14)),
    serial_parameters=example$serial_parameters,
    initial_cases_true=example$initial_cases,
    initial_reporting_parameters=list(mean=6, sd=4),
    initial_Rt=dplyr::tibble(Rt_index=1, Rt=1),
    print_to_screen=FALSE,
    is_negative_binomial=TRUE,
    initial_overdispersion=10,
    cases_before=example$cases_before,
    deaths=list(
      snapshot=snapshots %>%
        dplyr::select(time_onset, time_reported, deaths_reported),
      settings=list(initial_delay=list(mean=8, sd=5), triggered=triggered)),
    ...)
}

test_that("the proportion of fatal cases reported at death is sampled", {
  example <- triggered_example()
  set.seed(5)
  fit <- triggered_fit(example, 8)
  triggered <- fit$deaths_parameters$triggered
  expect_length(triggered, 8)
  expect_true(all(triggered > 0 & triggered < 1))
  expect_gt(length(unique(triggered)), 1)
  # without it the proportion stays at zero
  set.seed(5)
  fit <- triggered_fit(example, 3, triggered=FALSE)
  expect_equal(fit$deaths_parameters$triggered, rep(0, 3))
})

test_that("a fit with fatal cases reported at death can be taken up again", {
  example <- triggered_example()
  set.seed(7)
  whole <- triggered_fit(example, 8)
  set.seed(7)
  first <- triggered_fit(example, 3, return_state=TRUE)
  second <- triggered_fit(example, 5, resume_from=first$state)
  for(name in setdiff(names(whole), c("state", "iterations_done")))
    expect_equal(dplyr::bind_rows(first[[name]], second[[name]]), whole[[name]],
                 info=name)
})
