# reports of ten onset days of 20 cases each, all reported by day 5, of which
# 5 are shown as dead until the last report, which shows deaths_last
start_reports <- function(day_last, deaths_last, n=10) {
  days <- c(5L, 10L, 20L, as.integer(day_last))
  cases <- expand.grid(time_onset=1:n, time_reported=days)
  cases <- cases[cases$time_reported >= cases$time_onset, ]
  cases$cases_reported <- 20L
  cases$Rt_index <- 1L
  deaths <- cases[, c("time_onset", "time_reported")]
  deaths$deaths_reported <- ifelse(deaths$time_reported == day_last, deaths_last, 5L)
  list(cases=tibble::as_tibble(cases), deaths=tibble::as_tibble(deaths), n=n)
}

start_fit <- function(reports, niterations=30, settings=list(),
                      initial_revision_parameters=NULL, priors_extra=list()) {
  settings <- utils::modifyList(
    list(cfr_changepoints=FALSE, delay_changepoints=FALSE,
         initial_delay=list(mean=5, sd=2), triggered=FALSE), settings)
  priors <- c(list(Rt=list(shape=1, rate=1),
                   reporting=list(mean_mu=5, mean_sigma=3, sd_mu=3, sd_sigma=2),
                   max_cases=100, reporting_changepoints=list()),
              priors_extra)
  mcmc_single(
    niterations, reports$cases, priors=priors,
    serial_parameters=list(mean=6, sd=3),
    initial_cases_true=tibble::tibble(time_onset=seq_len(reports$n), cases_true=20L),
    initial_reporting_parameters=list(mean=2, sd=1),
    initial_Rt=tibble::tibble(Rt_index=1L, Rt=1),
    print_to_screen=FALSE, cases_before=rep(20, 20), fixed_cases=TRUE,
    initial_revision_parameters=initial_revision_parameters,
    deaths=list(snapshot=reports$deaths, settings=settings))
}

test_that("the sampler starts with more fatal cases than initial_cfr gives where the reports need them", {
  set.seed(1)
  fit <- start_fit(start_reports(30, 15))
  expect_true(all(fit$deaths$deaths_true >= 15))
  later <- fit$deaths[fit$deaths$iteration > 15, ]
  expect_true(abs(mean(later$cfr) - 0.75) < 0.1)
})

test_that("the sampler starts from fatal cases the reports allow when deaths can no longer be recorded", {
  # 300 days after onset, a death not yet recorded has probability zero under
  # the starting delay, so the start of initial_cfr = 0.5 (10 fatal cases of
  # 20) is impossible; the fatal cases must start at the 5 shown as dead
  for(day_last in c(60, 300)) {
    set.seed(2)
    fit <- start_fit(start_reports(day_last, 5))
    first <- fit$deaths[fit$deaths$iteration == 1, ]
    expect_true(all(first$deaths_true == 5))
    later <- fit$deaths[fit$deaths$iteration > 15, ]
    expect_true(abs(mean(later$cfr) - 0.25) < 0.1)
  }
})

test_that("the start is drawn from fatal cases the reports allow with wrong records", {
  set.seed(3)
  fit <- start_fit(start_reports(300, 5), niterations=10,
                   initial_revision_parameters=list(eta=0.05, rho=0.3),
                   priors_extra=list(revisions=list(eta_mu=0.05, eta_sigma=0.05,
                                                    rho_alpha=2, rho_beta=4),
                                     max_erroneous=2))
  first <- fit$deaths[fit$deaths$iteration == 1, ]
  # every fatal case is shown as dead by the last report, and wrong records
  # shown as dead can stand in for at most max_erroneous of the deaths
  expect_true(all(first$deaths_true >= 3 & first$deaths_true <= 5))
})

test_that("the start is drawn from fatal cases the reports allow with cases dated by death", {
  dated <- list(fraction=0.4, mean=6, sd=3, reporting=list(mean=2, sd=1))
  for(exact in c(FALSE, TRUE)) {
    set.seed(4)
    settings <- list(dated_by_death=if(exact) c(dated, list(exact=TRUE)) else dated)
    fit <- start_fit(start_reports(300, 5), niterations=10, settings=settings)
    first <- fit$deaths[fit$deaths$iteration == 1, ]
    # by day 300 every fatal case recorded for a day is shown as dead
    expect_true(all(first$deaths_recorded == 5))
    expect_true(all(first$dated_by_death <= first$deaths_recorded))
  }
})

test_that("possible_fatal_start draws only numbers the reports allow, and NULL if there are none", {
  f_forward <- function(x) list(logp=if(x %in% c(3, 4)) log(x) else -Inf)
  f_logp <- function(x, forward) forward$logp
  set.seed(5)
  draws <- replicate(400, possible_fatal_start(f_forward, f_logp, 0, 6)$deaths_true)
  expect_true(all(draws %in% c(3, 4)))
  expect_true(abs(mean(draws == 4) - 4 / 7) < 0.08)
  expect_null(possible_fatal_start(function(x) list(logp=-Inf), f_logp, 0, 6))
  expect_null(possible_fatal_start(f_forward, f_logp, 5, 4))
})

test_that("an onset time whose reports no number of fatal cases allows gives an error saying what to change", {
  reports <- start_reports(30, 15)
  # more deaths shown than cases reported is caught before sampling; here the
  # deaths fall, which no fatal cases can explain without wrong records
  reports$deaths$deaths_reported[reports$deaths$time_onset == 1 &
                                   reports$deaths$time_reported == 20] <- 18L
  expect_error(start_fit(reports, niterations=2), "initial_delay")
})
