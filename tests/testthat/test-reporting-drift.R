test_that("reporting_drift_starts gives fixed segments, the first long enough", {
  settings <- reporting_changepoint_settings(
    118, list(min_length=14, first_length=69, drift=7))
  expect_equal(reporting_drift_starts(118, settings),
               c(1, 70, 77, 84, 91, 98, 105, 112))
  # a short remainder is joined to the segment before
  expect_equal(reporting_drift_starts(114, settings),
               c(1, 70, 77, 84, 91, 98, 105))
  settings <- reporting_changepoint_settings(20, list(min_length=7, drift=7))
  expect_equal(reporting_drift_starts(20, settings), c(1, 8, 15))
  expect_equal(reporting_drift_starts(5, settings), 1)
  expect_error(reporting_changepoint_settings(20, list(drift=0)))
})

test_that("sample_reporting_drift returns its prior when there are no data", {
  set.seed(4)
  n_days <- 40
  settings <- reporting_changepoint_settings(
    n_days, list(min_length=7, drift=7, spread_sd=0.25, mean_step=0.3,
                 sd_step=0.3))
  prior <- list(mean_mu=14, mean_sigma=10, sd_mu=8, sd_sigma=8)
  state <- initial_reporting_changepoints(list(mean=14, sd=8))
  f_logp <- function(indices, reporting_parameters) 0
  draws <- matrix(NA, 4000, 2)
  for(i in seq_len(nrow(draws))) {
    state <- sample_reporting_changepoints(state, f_logp, n_days, prior, settings)
    draws[i, ] <- c(state$spread, state$parameters$mean[1])
  }
  draws <- draws[-(1:500), ]
  expect_equal(state$starts, c(1, 8, 15, 22, 29, 36))
  expect_equal(mean(draws[, 1]), 0.25 * sqrt(2 / pi), tolerance=0.15)
  # the first mean has prior sd 10 and moves slowly, so only a loose check
  expect_equal(mean(draws[, 2]), 14, tolerance=0.35)
})

test_that("sample_reporting_drift follows a delay which shortens", {
  set.seed(5)
  n_days <- 56
  truth <- seq(20, 10, length.out=n_days)
  settings <- reporting_changepoint_settings(
    n_days, list(min_length=7, drift=7, spread_sd=0.25, mean_step=0.12,
                 sd_step=0.12, n_segment_steps=2))
  prior <- list(mean_mu=14, mean_sigma=10, sd_mu=8, sd_sigma=8)
  f_logp <- function(indices, reporting_parameters)
    sum(stats::dnorm(log(reporting_parameters$mean), log(truth[indices]), 0.2, log=TRUE) +
          stats::dnorm(log(reporting_parameters$sd), log(truth[indices]), 0.2, log=TRUE))
  state <- initial_reporting_changepoints(list(mean=30, sd=5))
  means <- matrix(NA, 1500, 8)
  for(i in seq_len(nrow(means))) {
    state <- sample_reporting_changepoints(state, f_logp, n_days, prior, settings)
    means[i, ] <- state$parameters$mean
  }
  estimate <- colMeans(means[-(1:500), ])
  expect_equal(estimate[1], mean(truth[1:7]), tolerance=0.1)
  expect_equal(estimate[8], mean(truth[50:56]), tolerance=0.1)
  expect_true(all(diff(estimate) < 0))
})
