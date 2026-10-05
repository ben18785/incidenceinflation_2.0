test_that("compiled recursion matches R recursion without batches", {
  set.seed(30)
  reporting_parameters <- list(mean=4, sd=3)
  l_observations <- list(
    dplyr::tibble(time_reported=c(1, 3, 4, 6, 9), cases_reported=c(0, 2, 3, 3, 5)),
    dplyr::tibble(time_reported=c(4, 5, 7, 8, 30), cases_reported=c(3, 6, 5, 7, 6)),
    dplyr::tibble(time_reported=c(60, 70, 80), cases_reported=c(40, 40, 37)))
  l_revision_parameters <- list(
    list(eta=0.1, rho=0.3),
    list(eta=0, rho=0.3),
    list(eta=0.05, rho=0.02, pi=0.2, shifts=c(-1, 1, 2), shift_probs=c(0.5, 0.3, 0.2)))
  for(observation_df in l_observations) {
    cases_true <- (max(observation_df$cases_reported) - 6):(max(observation_df$cases_reported) + 15)
    cases_true <- cases_true[cases_true >= 0]
    num_periods <- nrow(observation_df) - 1
    l_entries <- list(
      NULL,
      list(initial=0.4, periods=seq(0.05, 0.3, length.out=num_periods)),
      list(initial=runif(length(cases_true)),
           periods=matrix(runif(length(cases_true) * num_periods),
                          nrow=length(cases_true))))
    for(revision_parameters in l_revision_parameters) {
      for(misdated_entries in l_entries) {
        for(max_erroneous in c(0, 3, 6)) {
          f_logp <- function(use_compiled)
            observation_process_revisions_logp(
              observation_df, cases_true, 1, reporting_parameters,
              revision_parameters, max_erroneous, misdated_entries,
              use_compiled=use_compiled)
          logp_r <- f_logp(FALSE)
          logp_compiled <- f_logp(TRUE)
          expect_equal(is.finite(logp_compiled), is.finite(logp_r))
          is_finite <- is.finite(logp_r)
          expect_equal(logp_compiled[is_finite], logp_r[is_finite])
        }
      }
    }
  }
})

test_that("compiled recursion is unaffected by a batch which has no effect", {
  observation_df <- dplyr::tibble(time_reported=c(4, 5, 7, 8, 30),
                                  cases_reported=c(3, 6, 5, 7, 6))
  reporting_parameters <- list(mean=4, sd=3)
  revision_parameters <- list(eta=0.1, rho=0.3)
  cases_true <- 3:20
  logp <- observation_process_revisions_logp(
    observation_df, cases_true, 1, reporting_parameters, revision_parameters, 4)

  # detected before first reporting day
  logp_batch <- observation_process_revisions_logp(
    observation_df, cases_true, 1, reporting_parameters, revision_parameters, 4,
    batch=list(phi=0.5, time_detected=4))
  expect_equal(logp_batch, logp)

  # zero intensity
  for(time_detected in c(5, 8, 100)) {
    logp_batch <- observation_process_revisions_logp(
      observation_df, cases_true, 1, reporting_parameters, revision_parameters, 4,
      batch=list(phi=0, time_detected=time_detected), max_batch=10)
    expect_equal(logp_batch, logp)
  }

  # batches need compiled code
  expect_error(observation_process_revisions_logp(
    observation_df, cases_true, 1, reporting_parameters, revision_parameters, 4,
    batch=list(phi=0.5, time_detected=8), use_compiled=FALSE))
})

test_that("compiled recursion with a batch matches enumeration over all
paths", {
  reporting_parameters <- list(mean=4, sd=3)
  revision_parameters <- list(eta=0.2, rho=0.4)
  day_onset <- 1
  max_erroneous <- 2
  max_batch <- 6
  phi <- 0.4

  # probability of reported cases by enumerating genuine reports, erroneous
  # entries and withdrawals, and batch entries, in each interval
  f_prob <- function(days, cases, case_true, time_detected) {
    f_undetected <- function(day) undetected_prob(day + 0.5, day_onset, reporting_parameters)
    f_path <- function(period, genuine, erroneous, batch) {
      if(period > length(days) - 1)
        return(1)
      p_detect <- detected_after_unobserved_prob(days[period + 1], days[period],
                                                 day_onset, reporting_parameters)
      p_new <- f_undetected(days[period]) - f_undetected(days[period + 1])
      p_withdrawn <- 1 - (1 - revision_parameters$rho)^(days[period + 1] - days[period])
      change <- cases[period + 1] - cases[period]
      is_detected <- days[period + 1] >= time_detected
      was_detected <- days[period] >= time_detected
      prob <- 0
      for(withdrawn in 0:erroneous) {
        for(entered in 0:(max_erroneous - erroneous + withdrawn)) {
          batch_changes <- if(was_detected) 0
          else if(is_detected) -batch
          else 0:(max_batch - batch)
          for(batch_change in batch_changes) {
            newly_reported <- change - entered + withdrawn - batch_change
            if(newly_reported < 0 | newly_reported > case_true - genuine)
              next
            prob_batch <- if(was_detected | is_detected) 1
            else dpois(batch_change, phi * case_true * p_new)
            prob_step <- dbinom(withdrawn, erroneous, p_withdrawn) *
              dpois(entered, revision_parameters$eta * case_true * p_new) *
              prob_batch *
              dbinom(newly_reported, case_true - genuine, p_detect)
            prob <- prob + prob_step * f_path(period + 1,
                                              genuine + newly_reported,
                                              erroneous + entered - withdrawn,
                                              batch + batch_change)
          }
        }
      }
      prob
    }
    # first reported on onset day: no erroneous or batch records
    f_path(1, cases[1], 0, 0)
  }

  days <- c(1, 2, 3, 5, 6, 8)
  l_cases <- list(c(0, 2, 4, 6, 2, 3), c(0, 1, 3, 3, 0, 1), c(0, 3, 5, 6, 6, 7))
  cases_true <- 1:9
  for(cases in l_cases) {
    for(time_detected in c(3, 6, 7, 20)) {
      prob_enumerated <- purrr::map_dbl(cases_true,
                                        ~f_prob(days, cases, ., time_detected))
      logp <- observation_process_revisions_logp(
        dplyr::tibble(time_reported=days, cases_reported=cases),
        cases_true, day_onset, reporting_parameters, revision_parameters,
        max_erroneous,
        batch=list(phi=phi, time_detected=time_detected), max_batch=max_batch)
      expect_equal(exp(logp), prob_enumerated)
    }
  }
})

test_that("compiled recursion with a batch matches frequencies of simulated
reported cases", {
  set.seed(31)
  reporting_parameters <- list(mean=3, sd=2)
  revision_parameters <- list(eta=0.1, rho=0.4)
  batches <- dplyr::tibble(start=2, width=2, phi=0.5, time_detected=6)
  # onset day 3 lies in footprint; days 2 to 8 are simulated
  cases_true <- c(0, 0, 5, 0, 0, 0, 0, 0)
  nreplicates <- 20000
  trajectories <- replicate(nreplicates, {
    df <- observed_cases_revisions(cases_true, reporting_parameters,
                                   revision_parameters, batches=batches)
    paste0(df$cases_reported[df$time_onset == 3], collapse="_")
  })
  frequencies <- sort(table(trajectories), decreasing=TRUE) / nreplicates
  l_trajectories <- purrr::map(strsplit(names(frequencies), "_"), as.numeric)

  # all batch records are withdrawn on detection, so counts fall then
  expect_true(mean(purrr::map_lgl(l_trajectories, ~.[4] < .[3])) > 0.5)

  for(i in 1:12) {
    prob <- exp(observation_process_revisions_logp(
      dplyr::tibble(time_reported=3:8, cases_reported=l_trajectories[[i]]),
      5, 3, reporting_parameters, revision_parameters, max_erroneous=7,
      batch=list(phi=0.5, time_detected=6), max_batch=20))
    std_error <- sqrt(prob * (1 - prob) / nreplicates)
    expect_true(abs(frequencies[[i]] - prob) < 4 * std_error)
  }
})

test_that("compiled recursion with a batch is unaffected by whether early
reporting days are observed", {
  reporting_parameters <- list(mean=3, sd=2)
  revision_parameters <- list(eta=0.2, rho=0.4)
  case_true <- 3
  max_erroneous <- 7
  max_batch <- 12
  batch <- list(phi=0.6, time_detected=5)
  max_reported <- case_true + max_erroneous + max_batch
  f_prob <- function(days, cases) {
    exp(observation_process_revisions_logp(
      dplyr::tibble(time_reported=days, cases_reported=cases),
      case_true, 1, reporting_parameters, revision_parameters, max_erroneous,
      batch=batch, max_batch=max_batch))
  }
  f_prob_long_way <- function(cases) {
    numerator <- 0
    denominator <- 0
    for(cases_hidden in 0:max_reported) {
      numerator <- numerator + f_prob(1:6, c(0, cases_hidden, cases))
      denominator <- denominator + f_prob(1:3, c(0, cases_hidden, cases[1]))
    }
    numerator / denominator
  }
  l_cases <- list(c(3, 4, 3, 3), c(5, 6, 2, 3), c(2, 3, 1, 2), c(4, 4, 1, 1))
  for(cases in l_cases) {
    prob <- f_prob(3:6, cases)
    expect_true(prob > 0)
    expect_equal(prob, f_prob_long_way(cases), tolerance=1e-4)
  }
})
