f_revision_parameters <- function(pi=0.2, eta=0.05, rho=0.3) {
  list(eta=eta, rho=rho, pi=pi,
       shifts=c(-2, -1, 1, 2), shift_probs=c(0.1, 0.4, 0.4, 0.1))
}

test_that("misdating_prob returns zero unless specified", {
  expect_equal(misdating_prob(list(eta=0.1, rho=0.2)), 0)
  expect_equal(misdating_prob(f_revision_parameters(0.3)), 0.3)
})

test_that("unrecorded_prob equals probability undetected without misdating", {
  reporting_parameters <- list(mean=4, sd=3)
  days <- c(5, 6, 9, 20)
  p_undetected <- purrr::map_dbl(days, ~undetected_prob(. + 0.5, 5, reporting_parameters))
  expect_equal(unrecorded_prob(days, 5, reporting_parameters, list(eta=0.1, rho=0.2)),
               p_undetected)
  expect_equal(unrecorded_prob(days, 5, reporting_parameters, f_revision_parameters(0)),
               p_undetected)
  # misdated cases are recorded later
  p_unrecorded <- unrecorded_prob(days, 5, reporting_parameters, f_revision_parameters(0.3))
  expect_equal(p_unrecorded[1], p_undetected[1])
  expect_true(all(p_unrecorded[-1] > p_undetected[-1]))
  expect_true(all(diff(p_unrecorded) < 0))
  expect_true(all(p_unrecorded > 0 & p_unrecorded < 1))
})

test_that("unrecorded_prob matches direct calculation", {
  reporting_parameters <- list(mean=4, sd=3)
  revision_parameters <- f_revision_parameters(0.3)
  f_s <- function(day) undetected_prob(day + 0.5, 1, reporting_parameters)
  # detected on day 2: only negative shifts possible (prob 0.5); day 3: shift of
  # one also possible (0.9); day 4: all possible
  p_awaiting <- (f_s(1) - f_s(2)) * 0.5 * 0.7^2 +
    (f_s(2) - f_s(3)) * 0.9 * 0.7 +
    (f_s(3) - f_s(4)) * 1
  expect_equal(unrecorded_prob(4, 1, reporting_parameters, revision_parameters),
               f_s(4) + 0.3 * p_awaiting)
})

test_that("misdated_entry_unit matches direct calculation", {
  reporting_parameters <- list(mean=4, sd=3)
  revision_parameters <- f_revision_parameters(0.3)
  observation_days <- c(6, 7, 9)
  # cases arising on day 3 recorded under day 4
  unit <- misdated_entry_unit(observation_days, 4, 1, reporting_parameters,
                              revision_parameters)
  f_s <- function(day) undetected_prob(day + 0.5, 3, reporting_parameters)
  expect_equal(unit$periods, 0.3 * 0.4 * c(f_s(6) - f_s(7), f_s(7) - f_s(9)))
  expect_equal(unit$initial,
               0.3 * 0.4 * ((f_s(4) - f_s(5)) * 0.7 + (f_s(5) - f_s(6))))

  # cases arising on a later day cannot be detected before they arise
  unit <- misdated_entry_unit(c(4, 5, 6, 7), 4, -2, reporting_parameters,
                              revision_parameters)
  f_s <- function(day) undetected_prob(day + 0.5, 6, reporting_parameters)
  expect_equal(unit$initial, 0)
  expect_equal(unit$periods[1], 0)
  expect_equal(unit$periods[2], 0.3 * 0.1 * (1 - f_s(6)))
  expect_equal(unit$periods[3], 0.3 * 0.1 * (f_s(6) - f_s(7)))
})

test_that("misdated_entries sums over neighbouring onset days", {
  reporting_parameters <- list(mean=4, sd=3)
  revision_parameters <- f_revision_parameters(0.3)
  observation_days <- c(6, 7, 9)
  cases_by_onset <- c("3"=10, "4"=7, "5"=20, "6"=4)
  entries <- misdated_entries(observation_days, 4, cases_by_onset,
                              reporting_parameters, revision_parameters)
  f_unit <- function(shift)
    misdated_entry_unit(observation_days, 4, shift, reporting_parameters,
                        revision_parameters)
  # neighbours of day 4 present are days 3 (shift 1), 5 (shift -1), 6 (shift -2)
  expect_equal(entries$periods,
               10 * f_unit(1)$periods + 20 * f_unit(-1)$periods + 4 * f_unit(-2)$periods)
  expect_equal(entries$initial,
               10 * f_unit(1)$initial + 20 * f_unit(-1)$initial + 4 * f_unit(-2)$initial)

  entries_excluding <- misdated_entries(observation_days, 4, cases_by_onset,
                                        reporting_parameters, revision_parameters,
                                        exclude=5)
  expect_equal(entries_excluding$periods,
               10 * f_unit(1)$periods + 4 * f_unit(-2)$periods)
})

test_that("simulated reported cases match unrecorded_prob and misdated entries", {
  set.seed(10)
  reporting_parameters <- list(mean=4, sd=3)
  revision_parameters <- f_revision_parameters(0.3, eta=0)
  # cases arise only on day 5 of 12
  cases_true <- c(0, 0, 0, 0, 400, 0, 0, 0, 0, 0, 0, 0)
  nreplicates <- 60
  l_reported <- purrr::map(1:nreplicates, ~observed_cases_revisions(
    cases_true, reporting_parameters, revision_parameters, days_max_follow_up=20))
  f_mean_reported <- function(time_onset, time_reported)
    mean(purrr::map_dbl(l_reported, function(df)
      df$cases_reported[df$time_onset == time_onset &
                          df$time_reported == time_reported]))

  # cases recorded under their own onset day
  p_unrecorded <- unrecorded_prob(c(5, 8, 12), 5, reporting_parameters,
                                  revision_parameters)
  for(i in 2:3) {
    expected <- 400 * (1 - p_unrecorded[i] / p_unrecorded[1])
    expect_equal(f_mean_reported(5, c(8, 12)[i - 1]), expected, tolerance=0.02)
  }

  # misdated records present under day 6 (shift of one) and day 3 (shift of minus two)
  for(shift in c(1, -2)) {
    unit <- misdated_entry_unit(c(10, 11), 5 + shift, shift, reporting_parameters,
                                revision_parameters)
    expect_equal(f_mean_reported(5 + shift, 10), 400 * unit$initial, tolerance=0.1)
  }

  # cases are conserved: every detected case is reported under some onset day
  totals <- purrr::map_dbl(l_reported, function(df)
    sum(df$cases_reported[df$time_reported == 12]))
  expect_true(all(totals <= 400))
  p_detected <- 1 - undetected_prob(12.5, 5, reporting_parameters) /
    undetected_prob(5.5, 5, reporting_parameters)
  expect_equal(mean(totals), 400 * p_detected, tolerance=0.02)
})

test_that("observation_process_revisions_logp handles misdated entries", {
  observation_matrix <- dplyr::tibble(time_reported=c(4, 5, 7, 8),
                                      cases_reported=c(3, 6, 5, 7))
  reporting_parameters <- list(mean=4, sd=3)
  cases_true <- 4:15

  # no misdating and no entries is as before
  logp <- observation_process_revisions_logp(
    observation_matrix, cases_true, 2, reporting_parameters, list(eta=0.1, rho=0.3))
  logp_zero <- observation_process_revisions_logp(
    observation_matrix, cases_true, 2, reporting_parameters, list(eta=0.1, rho=0.3),
    misdated_entries=list(initial=0, periods=rep(0, 3)))
  expect_equal(logp_zero, logp)

  # entries which vary by true cases match calculating each separately
  revision_parameters <- f_revision_parameters(0.2)
  initial <- seq(0.1, 1.2, length.out=12)
  periods <- matrix(seq(0.01, 0.9, length.out=36), nrow=12)
  logp_matrix <- observation_process_revisions_logp(
    observation_matrix, cases_true, 2, reporting_parameters, revision_parameters,
    misdated_entries=list(initial=initial, periods=periods))
  logp_long_way <- purrr::map_dbl(seq_along(cases_true), function(i)
    observation_process_revisions_logp(
      observation_matrix, cases_true[i], 2, reporting_parameters, revision_parameters,
      misdated_entries=list(initial=initial[i], periods=periods[i, ])))
  expect_equal(logp_matrix, logp_long_way)
  expect_true(all(is.finite(logp_matrix)))

  # misdated entries make falls in reported cases more likely
  logp_none <- observation_process_revisions_logp(
    observation_matrix, 8, 2, reporting_parameters, f_revision_parameters(0.2, eta=1e-6))
  logp_some <- observation_process_revisions_logp(
    observation_matrix, 8, 2, reporting_parameters, f_revision_parameters(0.2, eta=1e-6),
    misdated_entries=list(initial=0.5, periods=rep(0.3, 3)))
  expect_true(logp_some > logp_none)
})

f_misdating_data <- function(days_total=9, seed=20) {
  set.seed(seed)
  revision_parameters <- f_revision_parameters(0.25, eta=0.05, rho=0.3)
  Rt_function <- function(t) rep(1.3, length(t))
  df <- generate_snapshots(days_total, Rt_function, list(mean=2, sd=1),
                           list(mean=3, sd=2),
                           initial_parameters=list(mean=6, length=5),
                           revision_parameters=revision_parameters) %>%
    dplyr::mutate(reporting_piece_index=1)
  list(df=df, revision_parameters=revision_parameters, Rt_function=Rt_function,
       r_params=dplyr::tibble(reporting_piece_index=1, mean=3, sd=2))
}

test_that("misdating_observation_logp_function matches changes in
observation_process_all_times_logp", {
  l <- f_misdating_data()
  df <- l$df
  max_erroneous <- 8
  onset_times <- sort(unique(df$time_onset))
  l_observations <- split(as.data.frame(df[, c("time_reported", "cases_reported")]),
                          df$time_onset)
  cases_df <- df %>% dplyr::select(time_onset, cases_true) %>% unique()
  cases_by_onset <- stats::setNames(cases_df$cases_true + 2, cases_df$time_onset)
  l_reporting <- stats::setNames(purrr::map(onset_times, ~list(mean=3, sd=2)),
                                 onset_times)

  for(day_onset in c(1, 5, 9)) {
    f_logp <- misdating_observation_logp_function(
      day_onset, l_observations, cases_by_onset, l_reporting,
      l$revision_parameters, max_erroneous)
    cases <- cases_by_onset[[as.character(day_onset)]] + 0:4
    f_all_times <- function(case) {
      cases_tmp <- cases_by_onset
      cases_tmp[as.character(day_onset)] <- case
      df_tmp <- df %>%
        dplyr::mutate(cases_true=cases_tmp[as.character(time_onset)])
      observation_process_all_times_logp(df_tmp, l$r_params,
                                         l$revision_parameters, max_erroneous)
    }
    logp_all <- purrr::map_dbl(cases, f_all_times)
    logp <- f_logp(cases)
    expect_true(all(is.finite(logp)))
    expect_equal(diff(logp), diff(logp_all))
  }
})

test_that("sample_cases_history with misdating targets joint posterior of
case counts", {
  l <- f_misdating_data()
  df <- l$df %>% dplyr::select(-cases_true)
  s_params <- list(mean=2, sd=1)
  r_list <- list(mean=3, sd=2)
  max_erroneous <- max_downward_revision(df) + 2
  max_cases <- 26
  p_gamma_cutoff <- 0.3
  days_total <- max(df$time_onset)

  onset_times <- sort(unique(df$time_onset))
  latest_reported <- purrr::map_dbl(
    onset_times, ~dplyr::last(df$cases_reported[df$time_onset==.]))
  max_reported <- purrr::map_dbl(
    onset_times, ~max(df$cases_reported[df$time_onset==.]))
  days_uncertain <- onset_times[
    onset_times >= days_total - max_uncertain_days(p_gamma_cutoff, r_list)]
  expect_equal(days_uncertain, c(8, 9))

  # exact joint posterior by enumeration
  l_possible <- purrr::map(
    days_uncertain,
    ~max(0, max_reported[onset_times==.] - max_erroneous):max_cases)
  grid <- as.matrix(expand.grid(l_possible))
  cases_fixed <- latest_reported[onset_times < min(days_uncertain)]
  logp_joint <- rep(NA, nrow(grid))
  for(r in seq_len(nrow(grid))) {
    cases_all <- c(cases_fixed, grid[r, ])
    df_tmp <- df %>% dplyr::mutate(cases_true=cases_all[time_onset])
    logp <- observation_process_all_times_logp(
      df_tmp, l$r_params, l$revision_parameters, max_erroneous)
    for(k in seq_along(days_uncertain)) {
      day <- days_uncertain[k]
      logp <- logp + state_process_logp(
        cases_all[day], rev(cases_all[1:(day - 1)]), l$Rt_function(day), s_params)
    }
    logp_joint[r] <- logp
  }
  p_joint <- exp(logp_joint - matrixStats::logSumExp(logp_joint))
  mean_exact <- colSums(grid * p_joint)
  sd_exact <- sqrt(colSums(grid^2 * p_joint) - mean_exact^2)

  set.seed(21)
  nsweeps <- 800
  cases_current <- dplyr::tibble(time_onset=onset_times,
                                 cases_true=max_reported + 2)
  m_draws <- matrix(nrow=nsweeps, ncol=length(days_uncertain))
  for(i in 1:nsweeps) {
    df_est <- sample_cases_history(
      df, max_cases, l$Rt_function, s_params, l$r_params,
      p_gamma_cutoff=p_gamma_cutoff,
      cases_true_current=cases_current,
      revision_parameters=l$revision_parameters,
      max_erroneous=max_erroneous)
    cases_current <- df_est %>%
      dplyr::select(time_onset, cases_true=cases_estimated) %>%
      unique()
    m_draws[i, ] <- cases_current$cases_true[
      match(days_uncertain, cases_current$time_onset)]
  }
  z_scores <- (colMeans(m_draws) - mean_exact) / (sd_exact / sqrt(nsweeps / 3))
  expect_true(all(abs(z_scores) < 3))

  # current cases are needed
  expect_error(sample_cases_history(
    df, max_cases, l$Rt_function, s_params, l$r_params,
    revision_parameters=l$revision_parameters, max_erroneous=max_erroneous))
})

test_that("revision parameter proposals, prior and Jacobian include misdating
probability", {
  set.seed(22)
  current <- f_revision_parameters(0.2)
  steps <- list(eta_step=0.3, rho_step=0.3, pi_step=0.3)
  proposed <- propose_revision_parameters(current, steps)
  expect_true(proposed$pi != current$pi)
  expect_true(proposed$pi > 0 & proposed$pi < 1)
  expect_equal(proposed$shifts, current$shifts)
  expect_equal(proposed$shift_probs, current$shift_probs)
  # not proposed without step size
  proposed <- propose_revision_parameters(current, list(eta_step=0.3, rho_step=0.3))
  expect_equal(proposed$pi, current$pi)
  # nor without misdating
  proposed <- propose_revision_parameters(list(eta=0.1, rho=0.2), steps)
  expect_equal(sort(names(proposed)), c("eta", "rho"))

  prior <- list(eta_mu=0.01, eta_sigma=0.01, rho_alpha=2, rho_beta=5,
                pi_alpha=1, pi_beta=9)
  expect_equal(prior_revision_parameters(current, prior),
               dexp(0.05, 100, log=TRUE) + dbeta(0.3, 2, 5, log=TRUE) +
                 dbeta(0.2, 1, 9, log=TRUE))
  expect_equal(log_jacobian_revision_parameters(current),
               log(0.05 * 0.3 * 0.7 * 0.2 * 0.8))
})

test_that("mcmc works when cases can be misdated", {
  l <- f_misdating_data(days_total=30, seed=23)
  niter <- 3
  df <- l$df
  initial_cases_true <- df %>%
    dplyr::select(time_onset, cases_true) %>%
    unique()
  snapshot_with_Rt_index_df <- df %>%
    dplyr::select(-cases_true, -reporting_piece_index) %>%
    dplyr::mutate(Rt_index=ifelse(time_onset <= 15, 1, 2))
  priors <- list(Rt=list(shape=1, rate=1),
                 reporting=list(mean_mu=5, mean_sigma=10, sd_mu=3, sd_sigma=5),
                 revisions=list(eta_mu=0.1, eta_sigma=0.1,
                                rho_alpha=1, rho_beta=1,
                                pi_alpha=1, pi_beta=4),
                 max_erroneous=max_downward_revision(df) + 2,
                 max_cases=2000)
  f_mcmc <- function(initial_revision_parameters, ...) {
    mcmc(niterations=niter,
         snapshot_with_Rt_index_df,
         priors,
         list(mean=2, sd=1),
         initial_cases_true,
         list(mean=3, sd=2),
         dplyr::tibble(Rt_index=c(1, 2), Rt=c(1.3, 1.3)),
         print_to_screen=FALSE,
         initial_revision_parameters=initial_revision_parameters,
         revision_metropolis_parameters=list(eta_step=0.3, rho_step=0.3, pi_step=0.3),
         ...)
  }
  res <- f_mcmc(l$revision_parameters)
  expect_equal(nrow(res$cases), 30 * niter)
  expect_equal(sum(is.na(res$cases$cases_true)), 0)
  expect_equal(nrow(res$revisions), niter)
  expect_true(all(res$revisions$pi > 0 & res$revisions$pi < 1))
  expect_true(length(unique(c(l$revision_parameters$pi, res$revisions$pi))) > 1)

  # errors
  expect_error(f_mcmc(l$revision_parameters, condition_on_future_cases=FALSE))
  wrong <- l$revision_parameters
  wrong$shift_probs <- c(0.5, 0.5)
  expect_error(f_mcmc(wrong))
  wrong <- l$revision_parameters
  wrong$shifts <- c(-1, 0, 1, 2)
  expect_error(f_mcmc(wrong))
})
