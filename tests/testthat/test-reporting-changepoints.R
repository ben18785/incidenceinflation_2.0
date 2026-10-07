f_delay_data <- function(seed, days_total=40, starts=c(1, 21), means=c(3, 8),
                         sds=c(2, 4), revision_parameters=NULL,
                         days_max_follow_up=30) {
  set.seed(seed)
  cases_true <- true_cases(days_total, function(t) rep(1.1, length(t)), 1000,
                           list(mean=4, sd=2),
                           initial_parameters=list(mean=60, length=5))
  reporting_df <- piecewise_reporting_parameters(starts, means, sds, days_total)
  df <- observed_cases(cases_true, reporting_df, days_max_follow_up,
                       revision_parameters)
  list(df=df, reporting_df=reporting_df)
}

f_state <- function(l_or_n, starts=1, means=5, sds=3) {
  initial_reporting_changepoints(dplyr::tibble(mean=means, sd=sds), starts)
}

prior_parameters <- list(mean_mu=5, mean_sigma=10, sd_mu=3, sd_sigma=5)

test_that("reporting_logp_function agrees with the log-probability over all onset times", {
  l <- f_delay_data(1, days_total=12, starts=c(1, 7))
  f_logp <- reporting_logp_function(l$df)
  df <- l$df %>% dplyr::mutate(reporting_piece_index=ifelse(time_onset < 7, 1, 2))
  r_params <- dplyr::tibble(reporting_piece_index=1:2, mean=c(3, 8), sd=c(2, 4))
  expect_equal(
    f_logp(1:6, list(mean=3, sd=2)) + f_logp(7:12, list(mean=8, sd=4)),
    observation_process_all_times_logp(df, r_params))
  expect_equal(f_logp(integer(0), list(mean=3, sd=2)), 0)
  expect_equal(f_logp(1:3, list(mean=-1, sd=2)), -Inf)
  expect_equal(f_logp(1:3, list(mean=1, sd=0)), -Inf)
  # the parameters which generated a segment are more likely than the others'
  expect_true(f_logp(1:6, list(mean=3, sd=2)) > f_logp(1:6, list(mean=8, sd=4)))
  expect_true(f_logp(7:12, list(mean=8, sd=4)) > f_logp(7:12, list(mean=3, sd=2)))
})

test_that("reporting_logp_function allows for revisions and misdating", {
  revision_parameters <- list(eta=0.05, rho=0.3)
  l <- f_delay_data(2, days_total=10, starts=c(1, 6),
                    revision_parameters=revision_parameters)
  df <- l$df %>% dplyr::mutate(reporting_piece_index=ifelse(time_onset < 6, 1, 2))
  r_params <- dplyr::tibble(reporting_piece_index=1:2, mean=c(3, 8), sd=c(2, 4))
  max_erroneous <- max_downward_revision(l$df) + 2
  f_logp <- reporting_logp_function(l$df, revision_parameters, max_erroneous)
  expect_equal(
    f_logp(1:5, list(mean=3, sd=2)) + f_logp(6:10, list(mean=8, sd=4)),
    observation_process_all_times_logp(df, r_params, revision_parameters,
                                       max_erroneous))
  revision_parameters <- list(eta=0.05, rho=0.3, pi=0.1, shifts=c(-1, 1),
                              shift_probs=c(0.5, 0.5))
  f_logp <- reporting_logp_function(l$df, revision_parameters, max_erroneous + 6)
  expect_equal(
    f_logp(1:5, list(mean=3, sd=2)) + f_logp(6:10, list(mean=8, sd=4)),
    observation_process_all_times_logp(df, r_params, revision_parameters,
                                       max_erroneous + 6))
  # more than one true case count for an onset time is an error
  df_bad <- l$df
  df_bad$cases_true[1] <- df_bad$cases_true[1] + 1
  expect_error(reporting_logp_function(df_bad))
})

test_that("priors on segment parameters and their centre are as described", {
  parameters <- dplyr::tibble(mean=c(4, 6), sd=c(2, 5))
  centre <- list(mean=5, sd=3)
  expect_equal(
    reporting_segments_logp_prior(parameters, centre, 0.4),
    dlnorm(c(4, 6), log(5), 0.4, log=TRUE) + dlnorm(c(2, 5), log(3), 0.4, log=TRUE))
  expect_equal(
    reporting_centre_logp_prior(centre, 0.4, prior_parameters, 0.5),
    dgamma_mean_sd(5, 5, 10, log=TRUE) + dgamma_mean_sd(3, 3, 5, log=TRUE) +
      dnorm(0.4, 0, 0.5, log=TRUE))
  expect_equal(reporting_centre_logp_prior(list(mean=-1, sd=3), 0.4,
                                           prior_parameters, 0.5), -Inf)
  expect_equal(reporting_centre_logp_prior(centre, 0, prior_parameters, 0.5), -Inf)
})

test_that("n_splits counts the ways of splitting a segment", {
  expect_equal(n_splits(10, 1), 9)
  expect_equal(n_splits(10, 5), 1)
  expect_equal(n_splits(10, 6), 0)
  expect_equal(n_splits(14, 7), 1)
  expect_equal(n_splits(20, 7), 7)
})

test_that("initial_reporting_changepoints and settings are built correctly", {
  state <- initial_reporting_changepoints(list(mean=5, sd=3))
  expect_equal(state$starts, 1)
  expect_equal(state$parameters, dplyr::tibble(mean=5, sd=3))
  expect_equal(state$centre, list(mean=5, sd=3))
  state <- f_state(starts=c(1, 9), means=c(2, 8), sds=c(1, 4))
  expect_equal(state$centre$mean, 4)
  expect_equal(state$centre$sd, 2)
  expect_error(initial_reporting_changepoints(list(mean=5, sd=3), starts=c(1, 4)))

  settings <- reporting_changepoint_settings(50)
  expect_equal(settings$beta, 50)
  expect_equal(settings$min_length, 7)
  settings <- reporting_changepoint_settings(50, list(min_length=3, mean_step=1))
  expect_equal(settings$min_length, 3)
  expect_equal(settings$mean_step, 1)
  expect_error(reporting_changepoint_settings(50, list(nonsense=1)))
  expect_error(reporting_changepoint_settings(50, list(alpha=-1)))
  expect_error(reporting_changepoint_settings(50, list(min_length=0)))
})

test_that("moves leave the state consistent", {
  l <- f_delay_data(3, days_total=30, starts=c(1, 16))
  f_logp <- reporting_logp_function(l$df)
  settings <- reporting_changepoint_settings(30, list(min_length=3))
  state <- f_state(starts=c(1, 16), means=c(3, 8), sds=c(2, 4))
  set.seed(4)
  for(i in 1:40) {
    state <- sample_reporting_changepoints(state, f_logp, 30, prior_parameters,
                                           settings)
    expect_true(is_valid_partition(state$starts, 30, 3))
    expect_equal(nrow(state$parameters), length(state$starts))
    expect_true(all(state$parameters$mean > 0 & state$parameters$sd > 0))
    # log-probabilities held for each segment are up to date
    ends <- c(state$starts[-1] - 1, 30)
    logl <- purrr::map_dbl(seq_along(state$starts), function(j) {
      f_logp(state$starts[j]:ends[j], as.list(state$parameters[j, ]))
    })
    expect_equal(state$logl, logl)
  }
  expect_error(sample_reporting_changepoints(
    f_state(starts=c(1, 3), means=c(3, 8), sds=c(2, 4)), f_logp, 30,
    prior_parameters, reporting_changepoint_settings(30)))
})

test_that("with no data the moves sample the prior", {
  set.seed(5)
  n_days <- 12
  min_length <- 2
  settings <- reporting_changepoint_settings(
    n_days, list(alpha=2, beta=6, min_length=min_length, mean_step=2, sd_step=1.5))
  f_logp <- function(indices, reporting_parameters) 0
  state <- f_state()
  n_iterations <- 30000
  n_changes <- vector(length=n_iterations)
  first_mean <- vector(length=n_iterations)
  centre_mean <- vector(length=n_iterations)
  log_ratio <- vector(length=n_iterations)
  for(i in seq_len(n_iterations)) {
    state <- sample_reporting_changepoints(state, f_logp, n_days,
                                           prior_parameters, settings)
    n_changes[i] <- length(state$starts) - 1
    centre_mean[i] <- state$centre$mean
    log_ratio[i] <- log(state$parameters$mean[1] / state$centre$mean) / state$spread
  }
  # the number of changes follows the prior over valid partitions
  k <- 0:5
  probs <- purrr::map_dbl(k, ~n_valid_partitions(., n_days, min_length)) *
    exp(changepoints_logp_prior(k, n_days, 2, 6))
  probs <- probs / sum(probs)
  freqs <- as.numeric(table(factor(n_changes, levels=k))) / n_iterations
  expect_true(max(abs(freqs - probs)) < 0.03)
  # the centre follows its gamma prior
  expect_equal(mean(centre_mean), 5, tolerance=0.15)
  # parameters of segments are log-normal about the centre
  thinned <- log_ratio[seq(1, n_iterations, 10)]
  expect_equal(mean(thinned), 0, tolerance=0.1)
  expect_equal(sd(thinned), 1, tolerance=0.1)
})

test_that("a change in the delay is found, and none where there is none", {
  skip_on_cran()
  l <- f_delay_data(6, days_total=40, starts=c(1, 21), means=c(3, 8), sds=c(2, 4))
  f_logp <- reporting_logp_function(l$df)
  settings <- reporting_changepoint_settings(40)
  set.seed(7)
  state <- f_state()
  n_changes <- c()
  starts_second <- c()
  means_first <- c()
  means_last <- c()
  for(i in 1:500) {
    state <- sample_reporting_changepoints(state, f_logp, 40, prior_parameters,
                                           settings)
    if(i > 150) {
      n_changes <- c(n_changes, length(state$starts) - 1)
      starts_second <- c(starts_second, state$starts[2])
      means_first <- c(means_first, state$parameters$mean[1])
      means_last <- c(means_last, utils::tail(state$parameters$mean, 1))
    }
  }
  expect_true(mean(n_changes >= 1) > 0.95)
  expect_true(abs(median(starts_second, na.rm=TRUE) - 21) <= 3)
  expect_equal(median(means_first), 3, tolerance=0.2)
  expect_equal(median(means_last), 8, tolerance=0.3)

  l <- f_delay_data(8, days_total=40, starts=1, means=5, sds=3)
  f_logp <- reporting_logp_function(l$df)
  state <- f_state()
  n_changes <- c()
  for(i in 1:400) {
    state <- sample_reporting_changepoints(state, f_logp, 40, prior_parameters,
                                           settings)
    n_changes <- c(n_changes, length(state$starts) - 1)
  }
  expect_true(mean(n_changes == 0) > 0.8)
})

test_that("piecewise_reporting_parameters expands segments to onset times", {
  df <- piecewise_reporting_parameters(c(1, 4), c(3, 8), c(2, 4), 6)
  expect_equal(df$time_onset, 1:6)
  expect_equal(df$mean, c(3, 3, 3, 8, 8, 8))
  expect_equal(df$sd, c(2, 2, 2, 4, 4, 4))
  expect_error(piecewise_reporting_parameters(c(1, 4), 3, c(2, 4), 6))
})

test_that("simulated reports follow the delay of their segment", {
  l <- f_delay_data(9, days_total=60, starts=c(1, 31), means=c(2, 12), sds=c(1, 5),
                    days_max_follow_up=5)
  # proportion of cases reported within five days of onset
  proportions <- l$df %>%
    dplyr::filter(time_reported == time_onset + 5) %>%
    dplyr::mutate(segment=ifelse(time_onset < 31, 1, 2)) %>%
    dplyr::group_by(segment) %>%
    dplyr::summarise(proportion=sum(cases_reported) / sum(cases_true))
  expect_true(proportions$proportion[1] > 0.9)
  expect_true(proportions$proportion[2] < 0.25)
})

f_mcmc_delay <- function(l, niter, reporting_changepoints=list(),
                         Rt_changepoints=NULL, initial_revision_parameters=NULL,
                         ...) {
  df <- l$df
  initial_cases_true <- df %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(cases_true=max(cases_reported))
  snapshot_with_Rt_index_df <- df %>%
    dplyr::select(time_onset, time_reported, cases_reported) %>%
    dplyr::mutate(Rt_index=ceiling(time_onset / 10))
  priors <- list(Rt=list(shape=1, rate=1),
                 reporting=prior_parameters,
                 max_cases=5000)
  priors$reporting_changepoints <- reporting_changepoints
  priors$Rt_changepoints <- Rt_changepoints
  if(!is.null(initial_revision_parameters)) {
    priors$revisions <- list(eta_mu=0.05, eta_sigma=0.05, rho_alpha=1, rho_beta=1)
    priors$max_erroneous <- max_downward_revision(df) + 2
  }
  mcmc_single(niterations=niter,
              snapshot_with_Rt_index_df,
              priors,
              list(mean=4, sd=2),
              initial_cases_true,
              list(mean=5, sd=3),
              dplyr::tibble(Rt_index=1:max(snapshot_with_Rt_index_df$Rt_index), Rt=1),
              print_to_screen=FALSE,
              initial_revision_parameters=initial_revision_parameters,
              ...)
}

test_that("mcmc_single samples the times of changes in the reporting delay", {
  l <- f_delay_data(10, days_total=30, starts=c(1, 16))
  niter <- 4
  set.seed(11)
  res <- f_mcmc_delay(l, niter, list(min_length=3))
  expect_equal(colnames(res$reporting_days),
               c("iteration", "time_onset", "reporting_piece_index", "mean", "sd"))
  expect_equal(nrow(res$reporting_days), niter * 30)
  expect_equal(colnames(res$reporting_changepoints),
               c("iteration", "n_changes", "centre_mean", "centre_sd", "spread"))
  expect_equal(nrow(res$reporting_changepoints), niter)
  n_pieces <- res$reporting %>% dplyr::count(iteration)
  expect_equal(n_pieces$n, res$reporting_changepoints$n_changes + 1)
  expect_equal(nrow(res$cases), niter * 30)
  expect_equal(nrow(res$Rt), niter * 3)

  # without the prior on changes there is a single delay throughout
  res_fixed <- f_mcmc_delay(l, niter, NULL)
  expect_false("reporting_days" %in% names(res_fixed))
  expect_equal(nrow(res_fixed$reporting), niter)
  expect_error(f_mcmc_delay(l, 2, list(nonsense=1)))
  expect_error(f_mcmc_delay(l, 2, maximise=TRUE))
})

test_that("mcmc_single runs with changes in both Rt and the delay, and with revisions", {
  revision_parameters <- list(eta=0.03, rho=0.3)
  l <- f_delay_data(12, days_total=24, starts=c(1, 13),
                    revision_parameters=revision_parameters)
  set.seed(13)
  res <- f_mcmc_delay(l, 3, list(min_length=3), Rt_changepoints=list(min_length=3),
                      initial_revision_parameters=revision_parameters)
  expect_equal(nrow(res$reporting_days), 3 * 24)
  expect_equal(nrow(res$Rt), 3 * 24)
  expect_equal(nrow(res$revisions), 3)
  expect_true(all(c("Rt_changepoints", "reporting_changepoints") %in% names(res)))
})

test_that("mcmc recovers a change in the reporting delay", {
  skip_on_cran()
  l <- f_delay_data(14, days_total=50, starts=c(1, 21), means=c(3, 9), sds=c(2, 4))
  set.seed(15)
  res <- f_mcmc_delay(l, 300)
  summary_df <- res$reporting_days %>%
    dplyr::filter(iteration > 100) %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(mean=median(mean))
  expect_true(all(abs(summary_df$mean[1:17] - 3) < 0.8))
  expect_true(all(abs(summary_df$mean[25:40] - 9) < 2.5))
  expect_true(mean(res$reporting_changepoints$n_changes[-(1:100)] >= 1) > 0.9)
  combined <- combine_chains(list(res, res))
  expect_equal(nrow(combined$reporting_changepoints), 600)
})
