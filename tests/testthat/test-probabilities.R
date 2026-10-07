test_that("observation_process_single_logp returns reasonable log-probs", {
  reporting_parameters <- list(mean=5, sd=4)
  day_2 <- 4
  day_1 <- 2
  day_onset <- 0
  I_true <- 10
  I_day_1 <- 3

  # I_day_2 > I_true should be impossible
  I_day_2 <- 11
  logp <- observation_process_single_logp(I_true, I_day_2, I_day_1,
                                   day_2, day_1, day_onset, reporting_parameters)
  expect_equal(-Inf, logp)

  # I_day_2 < I_day_1 should be impossible
  I_day_2 <- 2
  logp <- observation_process_single_logp(I_true, I_day_2, I_day_1,
                                   day_2, day_1, day_onset, reporting_parameters)
  expect_equal(-Inf, logp)

  # I_day_2 near I_true should be more likely with shorter delays
  I_day_2 <- 10
  logp <- observation_process_single_logp(I_true, I_day_2, I_day_1,
                                   day_2, day_1, day_onset, reporting_parameters)
  reporting_parameters1 <- list(mean=1, sd=1)
  logp1 <- observation_process_single_logp(I_true, I_day_2, I_day_1,
                                    day_2, day_1, day_onset, reporting_parameters1)
  expect_true(logp1 > logp)

  # prob 1 if I_day_1=I_true
  I_day_2 <- 10
  I_day_1 <- 10
  logp <- observation_process_single_logp(I_true, I_day_2, I_day_1,
                                   day_2, day_1, day_onset, reporting_parameters)
  expect_equal(0, logp)
})

test_that("observation_process_single_logp can return vectors", {
  reporting_parameters <- list(mean=5, sd=4)
  day_2 <- 4
  day_1 <- 2
  day_onset <- 0
  I_true <- 10:20
  I_day_1 <- 3
  I_day_2 <- 5
  logps <- observation_process_single_logp(I_true, I_day_2, I_day_1,
                                          day_2, day_1, day_onset, reporting_parameters)
  n <- length(I_true)
  expect_equal(length(logps), n)
  log_p_long_way <- vector(length=n)
  for(i in 1:n)
    log_p_long_way[i] <- observation_process_single_logp(
      I_true[i], I_day_2, I_day_1, day_2, day_1, day_onset, reporting_parameters)
  expect_true(all.equal(logps, log_p_long_way))
})

test_that("observation_process_logp returns reasonable log prob values", {
  observation_matrix <- dplyr::tibble(time_reported=c(1, 20, 21),
                                      cases_reported=c(3, 5, 7))
  reporting_parameters <- list(mean=1, sd=1)
  logp <- observation_process_logp(observation_df=observation_matrix,
                                   cases_true = 10,
                                   day_onset = 1,
                                   reporting_parameters=reporting_parameters)
  expect_true(logp < 0)

  # test that prob declines as delay increases
  observation_matrix <- dplyr::tibble(time_reported=c(1, 20, 21),
                                      cases_reported=c(1, 1, 1))
  logp1 <- observation_process_logp(observation_df=observation_matrix,
                                   cases_true = 10,
                                   day_onset = 0,
                                   reporting_parameters=reporting_parameters)
  expect_true(logp1 < logp)

  observation_matrix <- dplyr::tibble(time_reported=c(1, 1, 5),
                                      cases_reported=c(3, 5, 7))
  expect_error(observation_process_logp(observation_df=observation_matrix,
                                   cases_true = 10,
                                   day_onset = 1,
                                   reporting_parameters=reporting_parameters))
  observation_matrix <- dplyr::tibble(time_reported=c(1, 2, 5),
                                      cases_reported=c(3, 2, 7))
  expect_error(observation_process_logp(observation_df=observation_matrix,
                                        cases_true = 10,
                                        day_onset = 1,
                                        reporting_parameters=reporting_parameters))
})

test_that("state_process_logp produces reasonable log prob values", {
  s_params <- list(mean=5, sd=2)
  Rt <- 1
  t_max <- 20
  cases_history <- rep(1, t_max)
  w <- weights_series(t_max, s_params)
  a_sum <- Rt * sum(w * cases_history)
  cases_true <- 1
  expect_equal(state_process_logp(cases_true, cases_history, Rt, s_params), -1)

  Rt <- 2
  a_sum <- Rt * sum(w * cases_history)
  cases_true <- 2
  expect_equal(state_process_logp(cases_true, cases_history, Rt, s_params), log(2 * exp(-2)))
})

test_that("state_process_logp produces works fine with vectors of true cases", {
  s_params <- list(mean=5, sd=2)
  Rt <- 1
  t_max <- 20
  cases_history <- rep(1, t_max)
  w <- weights_series(t_max, s_params)
  a_sum <- Rt * sum(w * cases_history)
  cases_true <- 1:5
  logps <- state_process_logp(cases_true, cases_history, Rt, s_params)
  n <- length(cases_true)
  log_p_long_way <- vector(length=n)
  for(i in 1:n)
    log_p_long_way[i] <- state_process_logp(cases_true[i], cases_history, Rt, s_params)
  expect_true(all.equal(logps, log_p_long_way))
})

test_that("state_process_logp works on with a negative binomial model", {
  s_params <- list(mean=5, sd=2)
  Rt <- 1
  t_max <- 20
  cases_history <- rep(1, t_max)
  w <- weights_series(t_max, s_params)
  a_sum <- Rt * sum(w * cases_history)
  cases_true <- 1
  logp <- state_process_logp(cases_true, cases_history, Rt, s_params,
                             is_negative_binomial=TRUE, kappa=2)
  expect_equal(logp, -1.216395, tolerance = 10^-3)
})

test_that("state_process_logp produces works fine with vectors of true cases with NB model", {
  s_params <- list(mean=5, sd=2)
  Rt <- 1
  t_max <- 20
  cases_history <- rep(1, t_max)
  w <- weights_series(t_max, s_params)
  a_sum <- Rt * sum(w * cases_history)
  cases_true <- 1:5
  logps <- state_process_logp(cases_true, cases_history, Rt, s_params,
                              is_negative_binomial = TRUE, kappa=2)
  n <- length(cases_true)
  log_p_long_way <- vector(length=n)
  for(i in 1:n)
    log_p_long_way[i] <- state_process_logp(cases_true[i], cases_history, Rt, s_params,
                                            is_negative_binomial = TRUE, kappa=2)
  expect_true(all.equal(logps, log_p_long_way))
})

test_that("conditional_cases_logp returns reasonable values", {
  observation_matrix <- dplyr::tibble(time_reported=c(1, 3, 5),
                                      cases_reported=c(1, 1, 1))
  reporting_parameters <- list(mean=1, sd=1)
  cases_true <- 80
  s_params <- list(mean=5, sd=2)
  Rt <- 1
  t_max <- 20
  cases_history <- rep(1, t_max)
  day_onset <- 0
  logp <- conditional_cases_logp(observation_df=observation_matrix,
                         cases_true=cases_true,
                         cases_history=cases_history,
                         Rt=Rt,
                         day_onset=day_onset,
                         serial_parameters=s_params,
                         reporting_parameters=reporting_parameters)
  expect_true(logp < 0)
  reporting_parameters <- list(mean=20, sd=10)
  logp1 <- conditional_cases_logp(observation_df=observation_matrix,
                                 cases_true=cases_true,
                                 cases_history=cases_history,
                                 Rt=Rt,
                                 day_onset=day_onset,
                                 serial_parameters=s_params,
                                 reporting_parameters=reporting_parameters)
  expect_true(logp1 > logp)
})

test_that("conditional_cases_logp works fine with vectorised true cases", {
  observation_matrix <- dplyr::tibble(time_reported=c(1, 3, 5),
                                      cases_reported=c(1, 1, 1))
  reporting_parameters <- list(mean=1, sd=1)
  cases_true <- 10:20
  s_params <- list(mean=5, sd=2)
  Rt <- 1
  t_max <- 20
  cases_history <- rep(1, t_max)
  day_onset <- 0
  logps <- conditional_cases_logp(observation_df=observation_matrix,
                                 cases_true=cases_true,
                                 cases_history=cases_history,
                                 Rt=Rt,
                                 day_onset=day_onset,
                                 serial_parameters=s_params,
                                 reporting_parameters=reporting_parameters)
  n <- length(cases_true)
  log_p_long_way <- vector(length=n)
  for(i in 1:n)
    log_p_long_way[i] <- conditional_cases_logp(observation_df=observation_matrix,
                                                cases_true=cases_true[i],
                                                cases_history=cases_history,
                                                Rt=Rt,
                                                day_onset=day_onset,
                                                serial_parameters=s_params,
                                                reporting_parameters=reporting_parameters)
  expect_true(all.equal(logps, log_p_long_way))
})

days_total <- 100
Rt_1 <- 1.5
Rt_2 <- 1.0
Rt_3 <- 1.3
v_Rt <- c(rep(Rt_1, 40), rep(Rt_2, 20), rep(Rt_3, 40))
Rt_function <- stats::approxfun(1:days_total, v_Rt)
s_params <- list(mean=5, sd=3)
r_params <- list(mean=2, sd=4)
kappa <- 1000
days_total <- 3
df <- generate_snapshots(days_total, Rt_function, s_params, r_params,
                         kappa=kappa)

test_that("observation_process_all_times_logp works ok for single reporting parameter index", {

  df <- df %>%
    dplyr::mutate(reporting_piece_index=1)

  # reporting parameters must be in tibble form when it reaches observation_process_all_times_logp
  r_params <- list(mean=2, sd=4)
  expect_error(observation_process_all_times_logp(df, r_params))

  # reporting parameters must contain reporting_piece_index column
  r_params <- dplyr::tibble(mean=r_params$mean, sd=r_params$sd)
  expect_error(observation_process_all_times_logp(df, r_params))

  # check that log probabilities of subperiods add up to overall
  r_params <- r_params %>%
    dplyr::mutate(reporting_piece_index=1)
  logp <- observation_process_all_times_logp(df, r_params)
  df1 <- df %>% dplyr::filter(time_onset==1)
  logp1 <- observation_process_logp(df1, df1$cases_true[1],
                                    1, r_params)
  df1 <- df %>% dplyr::filter(time_onset==2)
  logp2 <- observation_process_logp(df1, df1$cases_true[1],
                                    2,r_params)
  expect_equal(logp, logp1 + logp2)

  # check -inf if negative mu or sd
  r_params <- r_params %>%
    dplyr::mutate(mean=-1, sd=3)
  logp <- observation_process_all_times_logp(df, r_params)
  expect_equal(logp, -Inf)
  r_params <- r_params %>%
    dplyr::mutate(mean=1, sd=-3)
  logp <- observation_process_all_times_logp(df, r_params)
  expect_equal(logp, -Inf)

  # check passing more than one true cases per onset time
  df <- dplyr::tribble(
    ~time_onset, ~time_reported, ~cases_reported, ~cases_true,
    1, 1, 0, 49,
    1, 2, 0, 49,
    1, 3, 0, 50,
    2, 2, 0, 55,
    2, 3, 0, 55,
    3, 3, 0, 54
  ) %>%
    dplyr::mutate(reporting_piece_index=1)
  r_params <- r_params %>%
    dplyr::mutate(mean=1, sd=3)
  expect_error(
    observation_process_all_times_logp(df, r_params)
  )
})

test_that("observation_process_all_times_logp works ok for multiple reporting parameter index", {

  # check that log-p is additive over different pieces
  r_params <- dplyr::tibble(mean=c(5, 10), sd=c(3, 2)) %>%
    dplyr::mutate(reporting_piece_index=c(1, 2))
  df <- df %>%
    dplyr::mutate(reporting_piece_index=c(rep(1, 3), rep(2, 3)))
  logp <- observation_process_all_times_logp(df, r_params)

  df_1 <- df %>%
    dplyr::filter(reporting_piece_index==1)
  r_1 <- r_params %>%
    dplyr::filter(reporting_piece_index==1)
  logp1 <- observation_process_all_times_logp(df_1, r_1)
  df_2 <- df %>%
    dplyr::filter(reporting_piece_index==2)
  r_2 <- r_params %>%
    dplyr::filter(reporting_piece_index==2)
  logp2 <- observation_process_all_times_logp(df_2, r_2)
  expect_equal(logp, logp1 + logp2)

  # check that logp -Inf if any of means or sds < 0
  r_params_e <- r_params %>%
    dplyr::mutate(mean=c(1, -2))
  logp <- observation_process_all_times_logp(df, r_params_e)
  expect_equal(logp, -Inf)

  r_params_e <- r_params %>%
    dplyr::mutate(sd=c(1, -2))
  logp <- observation_process_all_times_logp(df, r_params_e)
  expect_equal(logp, -Inf)
})

test_that("state_process_nb_logp_all_onsets works ok", {
  s_params <- list(mean=5, sd=2)
  cases_history_rt_df <- dplyr::tribble(
    ~time_onset, ~cases_true, ~Rt,
    1, 2, 1.5,
    2, 4, 1.3,
    3, 5, 1.1
  )
  kappa <- 3
  logp_1 <- state_process_logp(cases_history_rt_df$cases_true[2],
                               cases_history_rt_df$cases_true[1],
                               kappa=kappa,
                               serial_parameters = s_params,
                               Rt=cases_history_rt_df$Rt[2],
                               is_negative_binomial=TRUE)
  logp_2 <- state_process_logp(cases_history_rt_df$cases_true[3],
                               rev(cases_history_rt_df$cases_true[1:2]),
                               kappa=kappa,
                               serial_parameters = s_params,
                               Rt=cases_history_rt_df$Rt[3],
                               is_negative_binomial=TRUE)
  logp_overall <- state_process_nb_logp_all_onsets(kappa=kappa,
                                   cases_history_rt_df=cases_history_rt_df,
                                   serial_parameters = s_params)
  expect_equal(logp_1 + logp_2, logp_overall)

  # returns -Inf if kappa 0 or negative
  val <- state_process_nb_logp_all_onsets(kappa=0,
                                   cases_history_rt_df=cases_history_rt_df,
                                   serial_parameters = s_params)
  expect_equal(val, -Inf)

  val <- state_process_nb_logp_all_onsets(kappa=-1,
                                          cases_history_rt_df=cases_history_rt_df,
                                          serial_parameters = s_params)
  expect_equal(val, -Inf)
})

test_that("future_cases_logp returns zeros when there are no later onset days", {
  s_params <- list(mean=5, sd=2)
  cases_history <- c(4, 3, 5, 2)
  logps <- future_cases_logp(cases_true=0:10,
                             cases_history=cases_history,
                             cases_future=numeric(0),
                             Rt_future=numeric(0),
                             serial_parameters=s_params)
  expect_equal(logps, rep(0, 11))
})

test_that("future_cases_logp throws error if Rt and cases different lengths", {
  s_params <- list(mean=5, sd=2)
  expect_error(future_cases_logp(cases_true=0:10,
                                 cases_history=c(4, 3, 5, 2),
                                 cases_future=c(3, 2),
                                 Rt_future=1.5,
                                 serial_parameters=s_params))
})

test_that("future_cases_logp matches sum of state_process_logp over later days", {
  s_params <- list(mean=3, sd=2)
  cases_history <- c(4, 3, 5, 2, 6)
  cases_future <- c(5, 0, 7, 3)
  Rt_future <- c(1.2, 1.2, 0.7, 0.7)
  cases_true <- 0:12

  f_long_way <- function(kappa, is_negative_binomial) {
    logps <- rep(0, length(cases_true))
    for(i in seq_along(cases_true)) {
      for(k in seq_along(cases_future)) {
        history <- c(rev(cases_future[seq_len(k - 1)]),
                     cases_true[i],
                     cases_history)
        logps[i] <- logps[i] + state_process_logp(
          cases_true=cases_future[k],
          cases_history=history,
          Rt=Rt_future[k],
          serial_parameters=s_params,
          kappa=kappa,
          is_negative_binomial=is_negative_binomial)
      }
    }
    logps
  }

  # Poisson
  logps <- future_cases_logp(cases_true, cases_history, cases_future,
                             Rt_future, s_params)
  expect_equal(logps, f_long_way(NULL, FALSE))

  # negative binomial
  kappa <- 3
  logps_nb <- future_cases_logp(cases_true, cases_history, cases_future,
                                Rt_future, s_params,
                                kappa=kappa, is_negative_binomial=TRUE)
  expect_equal(logps_nb, f_long_way(kappa, TRUE))
  expect_true(max(abs(logps - logps_nb)) > 0.01)
})

test_that("conditional_cases_logp with later cases is proportional to joint
probability of observations and all case counts", {
  observation_matrix <- dplyr::tibble(time_reported=c(6, 7, 9),
                                      cases_reported=c(1, 2, 3))
  reporting_parameters <- list(mean=3, sd=2)
  s_params <- list(mean=3, sd=2)
  day_onset <- 6
  cases_history <- c(4, 3, 5, 2, 6) # days 5, 4, ..., 1
  cases_future <- c(5, 0, 7) # days 7, 8, 9
  Rt <- 1.1
  Rt_future <- c(1.2, 0.7, 0.7)
  cases_true <- 3:20

  # joint log-probability of observations for day_onset and of the case counts
  # on all days after the first
  f_joint <- function(case_true) {
    cases_all <- c(rev(cases_history), case_true, cases_future)
    Rt_all <- c(rep(NA, length(cases_history)), Rt, Rt_future)
    logp <- observation_process_logp(observation_matrix, case_true,
                                     day_onset, reporting_parameters)
    for(day in day_onset:length(cases_all))
      logp <- logp + state_process_logp(
        cases_true=cases_all[day],
        cases_history=rev(cases_all[1:(day - 1)]),
        Rt=Rt_all[day],
        serial_parameters=s_params)
    logp
  }
  logp_joint <- purrr::map_dbl(cases_true, f_joint)

  logps <- conditional_cases_logp(cases_true, observation_matrix, cases_history,
                                  Rt, day_onset, s_params, reporting_parameters,
                                  cases_future=cases_future,
                                  Rt_future=Rt_future)
  expect_equal(logps, logp_joint)

  # without later cases it is not proportional to joint
  logps_past_only <- conditional_cases_logp(
    cases_true, observation_matrix, cases_history,
    Rt, day_onset, s_params, reporting_parameters)
  differences <- logps_past_only - logp_joint
  expect_true(diff(range(differences)) > 0.1)
})

test_that("erroneous_transitions enumerates valid changes in erroneous records", {
  df <- erroneous_transitions(0)
  expect_equal(nrow(df), 1)
  expect_equal(unlist(df[1, ], use.names=FALSE), rep(0, 4))

  max_erroneous <- 3
  df <- erroneous_transitions(max_erroneous)
  expect_equal(nrow(df), 30)
  expect_equal(nrow(unique(df)), 30)
  expect_true(all(df$withdrawn <= df$erroneous_previous))
  expect_true(all(df$entered >= 0))
  expect_true(all(df$erroneous_current >= 0 & df$erroneous_current <= max_erroneous))
  expect_equal(df$erroneous_current,
               df$erroneous_previous - df$withdrawn + df$entered)

  expect_error(erroneous_transitions(-1))
})

test_that("dbinom_log_safe returns -Inf outside of support", {
  x <- c(2, -1, 3, 0, 1)
  size <- c(4, 4, 2, -1, 1)
  logp <- dbinom_log_safe(x, size, 0.3)
  expect_equal(logp[1], dbinom(2, 4, 0.3, log=TRUE))
  expect_equal(logp[5], dbinom(1, 1, 0.3, log=TRUE))
  expect_equal(logp[2:4], rep(-Inf, 3))
})

test_that("observation_process_revisions_logp matches observation_process_logp
when no erroneous records are entered", {
  observation_matrix <- dplyr::tibble(time_reported=c(1, 3, 4, 6, 9),
                                      cases_reported=c(0, 2, 3, 3, 5))
  reporting_parameters <- list(mean=4, sd=2)
  cases_true <- 5:30
  logp <- observation_process_logp(observation_matrix, cases_true, 1,
                                   reporting_parameters)
  for(rho in c(0.01, 0.3, 0.9)) {
    for(max_erroneous in c(0, 1, 3)) {
      logp_revisions <- observation_process_logp(
        observation_matrix, cases_true, 1, reporting_parameters,
        revision_parameters=list(eta=0, rho=rho),
        max_erroneous=max_erroneous)
      expect_equal(logp_revisions, logp)
    }
  }

  # also when called directly
  logp_revisions <- observation_process_revisions_logp(
    observation_matrix, cases_true, 1, reporting_parameters,
    revision_parameters=list(eta=0, rho=0.5))
  expect_equal(logp_revisions, logp)

  # true cases below those reported are impossible without erroneous records
  logp_revisions <- observation_process_logp(
    observation_matrix, 4, 1, reporting_parameters,
    revision_parameters=list(eta=0, rho=0.5))
  expect_equal(logp_revisions, -Inf)
})

test_that("observation_process_revisions_logp matches enumeration over all
paths of genuine reports, erroneous entries and withdrawals", {
  observation_matrix <- dplyr::tibble(time_reported=c(2, 3, 5, 6),
                                      cases_reported=c(1, 3, 2, 4))
  reporting_parameters <- list(mean=4, sd=3)
  revision_parameters <- list(eta=0.3, rho=0.4)
  day_onset <- 1
  max_erroneous <- 2

  days <- observation_matrix$time_reported
  cases <- observation_matrix$cases_reported
  f_path_prob <- function(case_true, period, genuine, erroneous) {
    if(period > length(days) - 1)
      return(1)
    p_detect <- detected_after_unobserved_prob(days[period + 1], days[period],
                                               day_onset, reporting_parameters)
    p_newly_detected <-
      undetected_prob(days[period] + 0.5, day_onset, reporting_parameters) -
      undetected_prob(days[period + 1] + 0.5, day_onset, reporting_parameters)
    p_withdrawn <- 1 - (1 - revision_parameters$rho)^(days[period + 1] - days[period])
    change <- cases[period + 1] - cases[period]
    prob <- 0
    for(withdrawn in 0:erroneous) {
      for(entered in 0:(max_erroneous - erroneous + withdrawn)) {
        newly_reported <- change - entered + withdrawn
        if(newly_reported < 0 | newly_reported > case_true - genuine)
          next
        prob_step <- dbinom(withdrawn, erroneous, p_withdrawn) *
          dpois(entered, revision_parameters$eta * case_true * p_newly_detected) *
          dbinom(newly_reported, case_true - genuine, p_detect)
        prob <- prob + prob_step * f_path_prob(
          case_true, period + 1,
          genuine + newly_reported,
          erroneous + entered - withdrawn)
      }
    }
    prob
  }

  # distribution of erroneous records on first reporting day, given that one
  # day has passed since onset
  f_initial_prob <- function(case_true) {
    p_newly_detected <-
      undetected_prob(day_onset + 0.5, day_onset, reporting_parameters) -
      undetected_prob(days[1] + 0.5, day_onset, reporting_parameters)
    p_detect <- detected_after_unobserved_prob(days[1], day_onset, day_onset,
                                               reporting_parameters)
    erroneous <- 0:max_erroneous
    probs <- dpois(erroneous, revision_parameters$eta * case_true * p_newly_detected) *
      ifelse(cases[1] - erroneous >= 0,
             dbinom(cases[1] - erroneous, case_true, p_detect), 0)
    probs / sum(probs)
  }
  f_prob <- function(case_true) {
    probs_initial <- f_initial_prob(case_true)
    prob <- 0
    for(erroneous in 0:max_erroneous)
      if(probs_initial[erroneous + 1] > 0)
        prob <- prob + probs_initial[erroneous + 1] *
          f_path_prob(case_true, 1, cases[1] - erroneous, erroneous)
    prob
  }

  cases_true <- 1:12
  prob_enumerated <- purrr::map_dbl(cases_true, f_prob)
  logp <- observation_process_revisions_logp(
    observation_matrix, cases_true, day_onset, reporting_parameters,
    revision_parameters, max_erroneous)
  expect_equal(exp(logp), prob_enumerated)

  # true cases can be fewer than max reported but not by more than max erroneous
  expect_true(all(prob_enumerated[cases_true >= 2] > 0))
  expect_equal(prob_enumerated[1], 0)
})

test_that("observation_process_revisions_logp gives probabilities which sum to
one across all possible reported cases", {
  reporting_parameters <- list(mean=3, sd=2)
  revision_parameters <- list(eta=0.3, rho=0.4)
  case_true <- 3
  max_erroneous <- 7
  max_reported <- case_true + max_erroneous
  total <- 0
  for(cases_1 in 0:max_reported) {
    for(cases_2 in 0:max_reported) {
      observation_matrix <- dplyr::tibble(time_reported=c(1, 2, 4),
                                          cases_reported=c(0, cases_1, cases_2))
      total <- total + exp(observation_process_revisions_logp(
        observation_matrix, case_true, 1, reporting_parameters,
        revision_parameters, max_erroneous))
    }
  }
  expect_true(total <= 1 + 1e-10)
  expect_true(abs(total - 1) < 1e-5)
})

test_that("observation_process_revisions_logp matches frequencies of
simulated reported cases", {
  set.seed(42)
  reporting_parameters <- list(mean=3, sd=2)
  revision_parameters <- list(eta=0.25, rho=0.5)
  case_true <- 4
  day_onset <- 5
  days_reporting <- c(5, 6, 7, 9)
  nreplicates <- 20000
  trajectories <- replicate(
    nreplicates,
    paste0(observed_cases_trajectory(case_true, days_reporting, day_onset,
                                     reporting_parameters, revision_parameters),
           collapse="_"))
  frequencies <- sort(table(trajectories), decreasing=TRUE) / nreplicates

  # some reported cases were revised downwards and some exceeded true cases
  l_trajectories <- purrr::map(strsplit(names(frequencies), "_"), as.numeric)
  expect_true(any(purrr::map_lgl(l_trajectories, ~any(diff(.) < 0))))
  expect_true(any(purrr::map_lgl(l_trajectories, ~max(.) > case_true)))

  for(i in 1:10) {
    observation_matrix <- dplyr::tibble(time_reported=days_reporting,
                                        cases_reported=l_trajectories[[i]])
    prob <- exp(observation_process_revisions_logp(
      observation_matrix, case_true, day_onset, reporting_parameters,
      revision_parameters, max_erroneous=7))
    std_error <- sqrt(prob * (1 - prob) / nreplicates)
    expect_true(abs(frequencies[[i]] - prob) < 4 * std_error)
  }
})

test_that("observation_process_logp only allows downward revisions if
revision parameters supplied", {
  observation_matrix <- dplyr::tibble(time_reported=c(1, 2, 5, 7),
                                      cases_reported=c(3, 4, 3, 7))
  reporting_parameters <- list(mean=3, sd=2)
  expect_error(observation_process_logp(observation_matrix, 10, 1,
                                        reporting_parameters))
  logp <- observation_process_logp(observation_matrix, 10, 1,
                                   reporting_parameters,
                                   revision_parameters=list(eta=0.1, rho=0.2))
  expect_true(is.finite(logp))
  expect_true(logp < 0)

  # days must still be increasing
  observation_matrix$time_reported <- c(1, 2, 2, 7)
  expect_error(observation_process_logp(observation_matrix, 10, 1,
                                        reporting_parameters,
                                        revision_parameters=list(eta=0.1, rho=0.2)))
})

test_that("observation_process_revisions_logp returns -Inf for invalid
revision parameters", {
  observation_matrix <- dplyr::tibble(time_reported=c(1, 2, 5),
                                      cases_reported=c(3, 4, 3))
  reporting_parameters <- list(mean=3, sd=2)
  f <- function(eta, rho)
    observation_process_revisions_logp(observation_matrix, 5:7, 1,
                                       reporting_parameters,
                                       list(eta=eta, rho=rho))
  expect_equal(f(-0.1, 0.5), rep(-Inf, 3))
  expect_equal(f(0.1, 0), rep(-Inf, 3))
  expect_equal(f(0.1, 1), rep(-Inf, 3))
  expect_true(all(is.finite(f(0.1, 0.5))))
})

test_that("conditional_cases_revisions_logp matches conditional_cases_logp
for all non-negligible case counts", {
  observation_matrix <- dplyr::tibble(time_reported=c(6, 7, 9, 10),
                                      cases_reported=c(10, 22, 20, 25))
  reporting_parameters <- list(mean=3, sd=2)
  revision_parameters <- list(eta=0.05, rho=0.3)
  s_params <- list(mean=3, sd=2)
  cases_history <- c(40, 30, 50, 20, 60)
  cases_future <- c(50, 30, 70)
  cases_true <- 22:2000

  logps_full <- conditional_cases_logp(
    cases_true, observation_matrix, cases_history,
    Rt=1.1, day_onset=6, s_params, reporting_parameters,
    cases_future=cases_future, Rt_future=c(1.2, 0.7, 0.7),
    revision_parameters=revision_parameters)
  logps <- conditional_cases_revisions_logp(
    cases_true, observation_matrix, cases_history,
    Rt=1.1, day_onset=6, s_params, reporting_parameters,
    revision_parameters,
    cases_future=cases_future, Rt_future=c(1.2, 0.7, 0.7))

  is_evaluated <- is.finite(logps)
  expect_true(sum(is_evaluated) < 0.5 * length(cases_true))
  expect_equal(logps[is_evaluated], logps_full[is_evaluated])
  probs_full <- exp(logps_full - matrixStats::logSumExp(logps_full))
  expect_true(sum(probs_full[!is_evaluated]) < 1e-15)

  # single observation contributes nothing
  logps_single <- conditional_cases_revisions_logp(
    cases_true, observation_matrix[1, ], cases_history,
    Rt=1.1, day_onset=6, s_params, reporting_parameters,
    revision_parameters)
  expect_equal(logps_single,
               state_process_logp(cases_true, cases_history, 1.1, s_params))
})

test_that("observation_process_all_times_logp sums over onset times when
revisions are allowed", {
  reporting_parameters <- dplyr::tibble(reporting_piece_index=1, mean=3, sd=2)
  revision_parameters <- list(eta=0.1, rho=0.3)
  df <- dplyr::tibble(
    time_onset=c(1, 1, 1, 2, 2, 2, 3),
    time_reported=c(1, 2, 3, 2, 3, 4, 3),
    cases_reported=c(0, 3, 2, 1, 4, 4, 0),
    cases_true=c(4, 4, 4, 6, 6, 6, 2),
    reporting_piece_index=1)
  logp <- observation_process_all_times_logp(df, reporting_parameters,
                                             revision_parameters)
  r_list <- list(mean=3, sd=2)
  logp_long_way <- observation_process_revisions_logp(
    df[1:3, ], 4, 1, r_list, revision_parameters) +
    observation_process_revisions_logp(
      df[4:6, ], 6, 2, r_list, revision_parameters)
  expect_equal(logp, logp_long_way)

  # without revision parameters downwards revisions throw error
  expect_error(observation_process_all_times_logp(df, reporting_parameters))
})

test_that("erroneous_initial_logp gives distribution of erroneous records on
first reporting day", {
  reporting_parameters <- list(mean=4, sd=3)
  revision_parameters <- list(eta=0.3, rho=0.2)
  cases_true <- c(2, 4, 5, 8, 20)
  max_erroneous <- 3

  # no erroneous records if first reported on onset day or none are entered
  log_alpha <- erroneous_initial_logp(1, 2, cases_true, 1, reporting_parameters,
                                      revision_parameters, max_erroneous)
  expect_equal(dim(log_alpha), c(5, 4))
  expect_equal(exp(log_alpha[, 1]), rep(1, 5))
  expect_equal(exp(log_alpha[, -1]), matrix(0, nrow=5, ncol=3))
  log_alpha <- erroneous_initial_logp(4, 5, cases_true, 1, reporting_parameters,
                                      list(eta=0, rho=0.2), max_erroneous)
  expect_equal(exp(log_alpha[, 1]), rep(1, 5))

  # otherwise a distribution conditional on cases first reported
  log_alpha <- erroneous_initial_logp(4, 5, cases_true, 1, reporting_parameters,
                                      revision_parameters, max_erroneous)
  expect_equal(rowSums(exp(log_alpha)), rep(1, 5))
  # five reported cases with two true cases requires three erroneous records
  expect_equal(exp(log_alpha[1, ]), c(0, 0, 0, 1))
  # more true cases make erroneous records less likely to explain those reported
  expect_true(all(diff(exp(log_alpha[3:5, 1])) > 0))

  # matches direct calculation
  f_undetected <- function(day) undetected_prob(day + 0.5, 1, reporting_parameters)
  mean_erroneous <- 0.3 * 8 * (
    (f_undetected(1) - f_undetected(2)) * 0.8^2 +
      (f_undetected(2) - f_undetected(3)) * 0.8 +
      (f_undetected(3) - f_undetected(4)))
  p_detect <- 1 - f_undetected(4) / f_undetected(1)
  probs <- dpois(0:3, mean_erroneous) * dbinom(5 - 0:3, 8, p_detect)
  expect_equal(exp(log_alpha[4, ]), probs / sum(probs))

  # impossible if reported cases exceed true cases by more than max erroneous
  log_alpha <- erroneous_initial_logp(4, 6, cases_true, 1, reporting_parameters,
                                      revision_parameters, max_erroneous)
  expect_equal(log_alpha[1, ], rep(-Inf, 4))
  expect_true(all(is.finite(rowSums(exp(log_alpha[-1, ])))))
})

test_that("observation_process_revisions_logp is unaffected by whether early
reporting days are observed", {
  reporting_parameters <- list(mean=3, sd=2)
  revision_parameters <- list(eta=0.3, rho=0.4)
  case_true <- 3
  day_onset <- 1
  max_erroneous <- 8
  max_reported <- case_true + max_erroneous
  f_prob <- function(days, cases) {
    exp(observation_process_revisions_logp(
      dplyr::tibble(time_reported=days, cases_reported=cases),
      case_true, day_onset, reporting_parameters,
      revision_parameters, max_erroneous))
  }

  # probability of cases reported on days 4 to 6 given those on day 3, found
  # by summing over the cases which would have been reported on day 2
  f_prob_long_way <- function(cases) {
    numerator <- 0
    denominator <- 0
    for(cases_hidden in 0:max_reported) {
      numerator <- numerator + f_prob(1:6, c(0, cases_hidden, cases))
      denominator <- denominator + f_prob(1:3, c(0, cases_hidden, cases[1]))
    }
    numerator / denominator
  }

  l_cases <- list(c(2, 3, 3, 3), c(3, 2, 3, 3), c(4, 3, 3, 3), c(1, 0, 2, 1),
                  c(0, 1, 3, 4), c(5, 5, 4, 3))
  for(cases in l_cases) {
    prob <- f_prob(3:6, cases)
    expect_true(prob > 0)
    expect_equal(prob, f_prob_long_way(cases), tolerance=1e-5)
  }
})

test_that("observation_process_revisions_logp allows reported cases to fall
below those first reported long after onset", {
  observation_matrix <- dplyr::tibble(time_reported=c(60, 70, 80),
                                      cases_reported=c(40, 40, 37))
  reporting_parameters <- list(mean=4, sd=3)
  cases_true <- 34:42
  logp <- observation_process_revisions_logp(
    observation_matrix, cases_true, 1, reporting_parameters,
    list(eta=0.05, rho=0.02), max_erroneous=5)
  # true cases must be within max_erroneous of those reported
  expect_equal(logp[1], -Inf)
  expect_true(all(is.finite(logp[-1])))
  # true cases unlikely to exceed those latest reported so long after onset
  expect_true(all(logp[cases_true > 37] < logp[cases_true == 37] - 10))

  # but not if first reported on onset day
  observation_matrix$time_reported <- c(1, 70, 80)
  logp <- observation_process_revisions_logp(
    observation_matrix, cases_true, 1, reporting_parameters,
    list(eta=0.05, rho=0.02), max_erroneous=5)
  expect_equal(logp, rep(-Inf, length(cases_true)))
})

test_that("probability of detection remains defined long after onset", {
  reporting_parameters <- list(mean=5, sd=3)
  # matches direct calculation where that is possible
  p_direct <- (undetected_prob(3.5, 1, reporting_parameters) -
                 undetected_prob(6.5, 1, reporting_parameters)) /
    undetected_prob(3.5, 1, reporting_parameters)
  expect_equal(detected_after_unobserved_prob(6, 3, 1, reporting_parameters),
               p_direct)
  expect_equal(exp(log_undetected_prob(c(3.5, 6.5), 1, reporting_parameters)),
               c(undetected_prob(3.5, 1, reporting_parameters),
                 undetected_prob(6.5, 1, reporting_parameters)))

  # the probability of remaining undetected underflows to zero
  expect_equal(undetected_prob(400.5, 1, reporting_parameters), 0)
  p_late <- detected_after_unobserved_prob(410, 400, 1, reporting_parameters)
  expect_true(is.finite(p_late))
  expect_true(p_late > 0 & p_late <= 1)

  # so that observations long after onset have a defined probability
  observation_matrix <- dplyr::tibble(time_reported=c(400, 410, 420),
                                      cases_reported=c(30, 30, 30))
  logp <- observation_process_logp(observation_matrix, 30:32, 1,
                                   reporting_parameters)
  expect_equal(logp[1], 0)
  expect_true(all(is.finite(logp)))
  expect_true(all(diff(logp) < 0))
  logp_revisions <- observation_process_logp(
    observation_matrix, 30:32, 1, reporting_parameters,
    revision_parameters=list(eta=0, rho=0.2))
  expect_equal(logp_revisions, logp)
})

test_that("state_process_logp is zero for cases with no earlier cases", {
  logp <- state_process_logp(cases_true=0:5, cases_history=numeric(0), Rt=1.5,
                             serial_parameters=list(mean=5, sd=3))
  expect_equal(logp, rep(0, 6))
  logp <- state_process_logp(cases_true=3, cases_history=numeric(0), Rt=1.5,
                             serial_parameters=list(mean=5, sd=3),
                             kappa=2, is_negative_binomial=TRUE)
  expect_equal(logp, 0)
})

test_that("recent_history keeps the latest days and pads with zeros", {
  expect_equal(recent_history(c(5, 4, 3), 5), c(5, 4, 3, 0, 0))
  expect_equal(recent_history(c(5, 4, 3), 2), c(5, 4))
  expect_equal(recent_history(c(5, 4, 3), 3), c(5, 4, 3))
  expect_equal(recent_history(numeric(0), 2), c(0, 0))
})

test_that("state_process_logp uses the weights of the serial interval over serial_max days", {
  s_params <- list(mean=5, sd=3)
  w <- weights_series(10, s_params)
  # a history shorter than serial_max: earlier days contribute nothing and
  # the weights are not rescaled to the days available
  history <- c(7, 3, 9)
  expect_equal(
    state_process_logp(4, history, 1.5, s_params, serial_max=10),
    dpois(4, 1.5 * sum(w[1:3] * history), log=TRUE))
  # a history longer than serial_max: older days are ignored
  history <- c(7, 3, 9, 2, 8, 6, 4, 5, 1, 3, 20, 30)
  expect_equal(
    state_process_logp(4, history, 1.5, s_params, serial_max=10),
    dpois(4, 1.5 * sum(w * history[1:10]), log=TRUE))
  expect_equal(
    state_process_logp(4, history, 1.5, s_params, kappa=3,
                       is_negative_binomial=TRUE, serial_max=10),
    dnbinom(4, mu=1.5 * sum(w * history[1:10]), size=3, log=TRUE))
  # no history at all says nothing
  expect_equal(state_process_logp(c(4, 5), numeric(0), 1.5, s_params), c(0, 0))
})

test_that("the renewal mean is the same wherever it is calculated", {
  set.seed(1)
  s_params <- list(mean=6, sd=4)
  serial_max <- 12
  cases <- rpois(30, 20)
  Rt <- runif(30, 0.5, 2)
  statistics <- renewal_statistics(cases, s_params, serial_max)
  f_logp <- function(t) {
    state_process_logp(cases[t], rev(cases[seq_len(t - 1)]), Rt[t], s_params,
                       serial_max=serial_max)
  }
  # the case and overdispersion updates agree with the Rt update, both early
  # in the series and after more than serial_max days
  for(t in c(2, 3, 8, 13, 20, 30))
    expect_equal(f_logp(t),
                 dpois(cases[t], Rt[t] * statistics$infectiousness[t], log=TRUE))

  # the contribution of later days to the update of a case count
  for(t in c(1, 4, 15)) {
    days_future <- (t + 1):min(30, t + serial_max)
    logp <- future_cases_logp(cases[t], rev(cases[seq_len(t - 1)]),
                              cases[days_future], Rt[days_future], s_params,
                              serial_max=serial_max)
    expect_equal(logp, sum(purrr::map_dbl(days_future, f_logp)))
    # and as a function of the count: differences match recalculation
    cases_alt <- cases
    cases_alt[t] <- cases[t] + 5
    statistics_alt <- renewal_statistics(cases_alt, s_params, serial_max)
    logp_alt <- future_cases_logp(cases_alt[t], rev(cases[seq_len(t - 1)]),
                                  cases[days_future], Rt[days_future], s_params,
                                  serial_max=serial_max)
    expect_equal(
      logp_alt,
      sum(dpois(cases[days_future],
                Rt[days_future] * statistics_alt$infectiousness[days_future],
                log=TRUE)))
  }
  # days more than serial_max later do not depend on the count
  far <- future_cases_logp(c(cases[1], cases[1] + 50), numeric(0),
                           cases[2:30], Rt[2:30], s_params, serial_max=serial_max)
  near <- future_cases_logp(c(cases[1], cases[1] + 50), numeric(0),
                            cases[2:13], Rt[2:13], s_params, serial_max=serial_max)
  expect_equal(diff(far), diff(near))

  # the negative binomial log-probability over all onset times
  df <- dplyr::tibble(time_onset=1:30, cases_true=cases, Rt=Rt)
  informative <- statistics$infectiousness > 0 & seq_len(30) > 1
  expect_equal(
    state_process_nb_logp_all_onsets(4, df, s_params, serial_max=serial_max),
    sum(dnbinom(cases[informative],
                mu=Rt[informative] * statistics$infectiousness[informative],
                size=4, log=TRUE)))

  # and the likelihood of a piece used when Rt is updated
  w <- weights_series(serial_max, s_params)
  expect_equal(
    nb_log_likelihood_Rt_piece(1.3, 4, w, 21:30, df),
    sum(dnbinom(cases[21:30], mu=1.3 * statistics$infectiousness[21:30],
                size=4, log=TRUE)))
})
