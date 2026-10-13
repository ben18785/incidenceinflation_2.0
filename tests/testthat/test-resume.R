# taking sampling up where an earlier call stopped

resume_example <- function() {
  set.seed(26)
  n_days <- 35
  serial_parameters <- list(mean=5, sd=3)
  cases <- true_cases(n_days + 40, function(t) 1, kappa=1000, serial_parameters,
                      initial_parameters=list(mean=30, length=10))
  cases_before <- cases[1:40]
  cases <- cases[41:(n_days + 40)]
  sim <- simulate_cases_deaths(cases, seq_len(n_days + 30), list(mean=4, sd=3), 0.5,
                               list(mean=8, sd=4),
                               wrong=list(eta=0.05, rho=0.3, delta=0.5))
  snapshots <- sim$snapshots
  list(snapshots=snapshots, cases_before=cases_before,
       serial_parameters=serial_parameters,
       initial_cases=snapshots %>%
         dplyr::group_by(time_onset) %>%
         dplyr::summarise(cases_true=max(cases_reported)),
       largest_fall=max_downward_revision(snapshots))
}

resume_fit <- function(example, niterations, deaths=TRUE, delay_changepoints=list(),
                       ...) {
  snapshots <- example$snapshots
  mcmc_single(
    niterations=niterations,
    snapshot_with_Rt_index_df=snapshots %>%
      dplyr::select(time_onset, time_reported, cases_reported) %>%
      dplyr::mutate(Rt_index=1),
    priors=list(Rt=list(shape=1, rate=1),
                reporting=list(mean_mu=5, mean_sigma=5, sd_mu=3, sd_sigma=3),
                max_cases=300,
                overdispersion=list(mean=20, sd=20),
                revisions=list(eta_mu=0.05, eta_sigma=0.05, rho_alpha=2, rho_beta=4),
                max_erroneous=example$largest_fall + 1,
                Rt_changepoints=list(),
                reporting_changepoints=list(min_length=14)),
    serial_parameters=example$serial_parameters,
    initial_cases_true=example$initial_cases,
    initial_reporting_parameters=list(mean=4, sd=3),
    initial_Rt=dplyr::tibble(Rt_index=1, Rt=1),
    print_to_screen=FALSE,
    is_negative_binomial=TRUE,
    initial_overdispersion=10,
    cases_before=example$cases_before,
    initial_revision_parameters=list(eta=0.05, rho=0.3, weekday=TRUE, weekday_first=3,
                                     tail=list(weight=0.05, mean=20)),
    deaths=if(deaths) list(
      snapshot=snapshots %>%
        dplyr::select(time_onset, time_reported, deaths_reported),
      settings=list(initial_delay=list(mean=8, sd=4),
                    delay_changepoints=delay_changepoints)),
    ...)
}

expect_same_draws <- function(whole, first, second) {
  names_draws <- setdiff(names(whole), c("state", "iterations_done"))
  expect_setequal(setdiff(names(second), c("state", "iterations_done")), names_draws)
  for(name in names_draws)
    expect_equal(dplyr::bind_rows(first[[name]], second[[name]]), whole[[name]],
                 info=name)
}

test_that("a fit of cases and deaths taken up again makes the draws it would have made", {
  example <- resume_example()
  set.seed(3)
  whole <- resume_fit(example, 12)
  set.seed(3)
  first <- resume_fit(example, 6, return_state=TRUE)
  expect_equal(first$iterations_done, 6)
  expect_equal(first$state$iterations, 6)
  # the random numbers drawn in between make no difference
  set.seed(99)
  stats::runif(5)
  second <- resume_fit(example, 6, resume_from=first$state)
  expect_equal(range(second$cases$iteration), c(7, 12))
  expect_equal(second$state$iterations, 12)
  expect_same_draws(whole, first, second)
})

test_that("a fit with a drifting delay to a death being recorded can be taken up again", {
  example <- resume_example()
  drift <- list(drift=7, spread_sd=0.25, mean_step=0.12, sd_step=0.12)
  set.seed(4)
  whole <- resume_fit(example, 8, delay_changepoints=drift)
  set.seed(4)
  first <- resume_fit(example, 3, delay_changepoints=drift, return_state=TRUE)
  second <- resume_fit(example, 5, delay_changepoints=drift, resume_from=first$state)
  expect_same_draws(whole, first, second)
})

test_that("a fit of cases alone can be taken up again, more than once", {
  example <- resume_example()
  set.seed(5)
  whole <- resume_fit(example, 12, deaths=FALSE)
  set.seed(5)
  first <- resume_fit(example, 4, deaths=FALSE, return_state=TRUE)
  second <- resume_fit(example, 4, deaths=FALSE, resume_from=first$state)
  third <- resume_fit(example, 4, deaths=FALSE, resume_from=second$state)
  expect_equal(range(third$cases$iteration), c(9, 12))
  for(name in setdiff(names(whole), c("state", "iterations_done")))
    expect_equal(dplyr::bind_rows(first[[name]], second[[name]], third[[name]]),
                 whole[[name]], info=name)
})

test_that("sampling stops at a time limit and can be taken up from there", {
  example <- resume_example()
  set.seed(6)
  whole <- resume_fit(example, 6, deaths=FALSE)
  set.seed(6)
  first <- resume_fit(example, 6, deaths=FALSE, time_limit=1e-6)
  expect_equal(first$iterations_done, 1)
  expect_equal(unique(first$cases$iteration), 1)
  expect_equal(unique(first$Rt$iteration), 1)
  second <- resume_fit(example, 5, deaths=FALSE, resume_from=first$state)
  expect_same_draws(whole, first, second)
  # a fit with no time limit and no state asked for returns what it did before
  expect_false("state" %in% names(whole))
})

test_that("a state is not taken up with different reports", {
  example <- resume_example()
  set.seed(7)
  first <- resume_fit(example, 2, deaths=FALSE, return_state=TRUE)
  other <- example
  other$snapshots$cases_reported[1] <- other$snapshots$cases_reported[1] + 1
  other$initial_cases <- other$snapshots %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(cases_true=max(cases_reported))
  expect_error(resume_fit(other, 2, deaths=FALSE, resume_from=first$state),
               "different reports")
  expect_error(resume_fit(example, 2, deaths=FALSE, resume_from=list(a=1)))
  expect_error(resume_fit(example, 2, deaths=FALSE, time_limit=-1))
})
