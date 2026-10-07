test_that("sample_true_cases_single_onset produces reasonable case draws", {
  observation_matrix <- dplyr::tibble(time_reported=c(1, 3, 5),
                                      cases_reported=c(1, 1, 1))
  reporting_parameters <- list(mean=5, sd=3)
  max_cases <- 100
  s_params <- list(mean=10, sd=1)
  Rt <- 2
  t_max <- 30
  cases_history <- rep(4, t_max)
  day_onset <- 0
  w <- weights_series(t_max, s_params)
  mean_cases <- expected_cases(Rt, w, cases_history)
  case <- sample_true_cases_single_onset(observation_df=observation_matrix,
                                 cases_history=cases_history,
                                 max_cases=max_cases,
                                 Rt=Rt,
                                 day_onset=day_onset,
                                 serial_parameters=s_params,
                                 reporting_parameters=reporting_parameters)
  expect_true(case >= min(observation_matrix$cases_reported))
  expect_true(case <= max_cases)

  # when reporting mean is low, expect to have pretty much seen all cases
  reporting_parameters <- list(mean=1, sd=1)
  cases <- sample_true_cases_single_onset(observation_df=observation_matrix,
                            cases_history=cases_history,
                            max_cases=max_cases,
                            Rt=Rt,
                            day_onset=day_onset,
                            serial_parameters=s_params,
                            reporting_parameters=reporting_parameters,
                            ndraws=20)
  expect_true(any(cases <= 5))

  # check that cases jump if serial interval looks further back
  cases_history <- c(rep(4, t_max / 2), rep(1000, t_max / 2))
  s_params <- list(mean=1, sd=1)
  cases <- sample_true_cases_single_onset(observation_df=observation_matrix,
                             cases_history=cases_history,
                             max_cases=max_cases,
                             Rt=Rt,
                             day_onset=day_onset,
                             serial_parameters=s_params,
                             reporting_parameters=reporting_parameters,
                             ndraws=20)
  s_params <- list(mean=25, sd=1)
  cases1 <- sample_true_cases_single_onset(observation_df=observation_matrix,
                             cases_history=cases_history,
                             max_cases=max_cases,
                             Rt=Rt,
                             day_onset=day_onset,
                             serial_parameters=s_params,
                             reporting_parameters=reporting_parameters,
                             ndraws=20)
  expect_true(mean(cases1) > mean(cases))

  max_cases <- 2
  expect_warning(sample_true_cases_single_onset(observation_df=observation_matrix,
                                                cases_history=cases_history,
                                                max_cases=max_cases,
                                                Rt=Rt,
                                                day_onset=day_onset,
                                                serial_parameters=s_params,
                                                reporting_parameters=reporting_parameters,
                                                ndraws=20))
  max_cases <- 0
  expect_error(sample_true_cases_single_onset(
    observation_df=observation_matrix,
    cases_history=cases_history,
    max_cases=max_cases,
    Rt=Rt,
    day_onset=day_onset,
    serial_parameters=s_params,
    reporting_parameters=reporting_parameters,
    ndraws=20))
})

test_that("sample_true_cases_single_onset produces consistent maximum", {
  observation_matrix <- dplyr::tibble(time_reported=c(1, 3, 5),
                                      cases_reported=c(1, 1, 1))
  reporting_parameters <- list(mean=5, sd=3)
  max_cases <- 100
  s_params <- list(mean=10, sd=1)
  Rt <- 2
  t_max <- 30
  cases_history <- rep(4, t_max)
  day_onset <- 0
  max_observed_cases <- max(observation_matrix$cases_reported)
  possible_cases <- max_observed_cases:max_cases
  logps <- conditional_cases_logp(possible_cases, observation_matrix, cases_history,
                                  Rt, day_onset, s_params, reporting_parameters)

  map_val <- possible_cases[which.max(logps)]
  case <- sample_true_cases_single_onset(
    observation_df=observation_matrix,
    cases_history=cases_history,
    max_cases=max_cases,
    Rt=Rt,
    day_onset=day_onset,
    serial_parameters=s_params,
    reporting_parameters=reporting_parameters,
    maximise=T)
  expect_equal(case, map_val)
})

test_that("max_uncertain_days selects reasonable gamma distribution quanties", {
  r_params <- list(mean=7, sd=0.01)
  expect_equal(round(max_uncertain_days(0.5, r_params)), 7)
  expect_true(max_uncertain_days(0.95, r_params) > 7)

  r_params <- list(mean=7, sd=5)
  expect_true(max_uncertain_days(0.95, r_params) > 10)
  expect_equal(max_uncertain_days(0, r_params), 0)
})

test_that("sample_cases_history adds cases_estimated that look reasonable", {
  days_total <- 100
  kappa <- 1000
  r_params <- list(mean=10, sd=3)
  v_Rt <- c(rep(1.5, 40), rep(0.4, 20), rep(1.5, 40))
  Rt_function <- stats::approxfun(1:days_total, v_Rt)
  s_params <- list(mean=5, sd=3)
  df <- generate_snapshots(days_total, Rt_function,
                           s_params, r_params,
                           kappa=kappa, thinned=T) %>%
    dplyr::mutate(reporting_piece_index=1)
  max_cases <- 5000
  r_params <- dplyr::tibble(mean=10, sd=3, reporting_piece_index=1)
  df_est <- sample_cases_history(df, max_cases, Rt_function, s_params, r_params)
  expect_true(all.equal(df_est %>% dplyr::select(-cases_estimated),
                        df))
  expect_true(max(df_est$cases_estimated) < max_cases)
  expect_true(min(df_est$cases_estimated) > 0)
  expect_equal(sum(is.na(df_est$cases_estimated)), 0)

  df_group <- df_est %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(
      cases_reported=max(cases_reported),
      cases_estimated=min(cases_estimated),
      cases_true=mean(cases_true)) %>%
    dplyr::mutate(diff=cases_estimated-cases_reported) %>%
    dplyr::mutate(diff_true=cases_estimated-cases_true)
  expect_equal(sum(df_group$diff < 0), 0)
  expect_true(max(abs(df_group$diff_true)) < 300)

  # throws error when reporting_piece_index not in observation_onset_df
  df_tmp <- df %>%
    dplyr::select(-reporting_piece_index)
  expect_error(sample_cases_history(df_tmp, max_cases, Rt_function, s_params, r_params))
})

test_that("sample_cases_history yields a single case history when
maximising", {
  days_total <- 30
  kappa <- 1000
  r_params <- list(mean=10, sd=3)
  v_Rt <- c(rep(1.5, 10), rep(0.4, 10), rep(1.5, 10))
  Rt_function <- stats::approxfun(1:days_total, v_Rt)
  s_params <- list(mean=5, sd=3)
  df <- generate_snapshots(days_total, Rt_function,
                           s_params, r_params,
                           kappa=kappa, thinned=T) %>%
    dplyr::mutate(reporting_piece_index=1)
  max_cases <- 5000
  r_params <- dplyr::tibble(
    reporting_piece_index=1,
    mean=r_params$mean,
    sd=r_params$sd
  )
  f_est <- function(i) {
    df_est <- sample_cases_history(df, max_cases, Rt_function,
                         s_params, r_params,
                         maximise = T)
    df_est$cases_estimated
  }
  cases <- purrr::map(seq(1, 4, 1), f_est)
  expect_true(all.equal(cases[[1]], cases[[2]]))
  expect_true(all.equal(cases[[2]], cases[[3]]))
  expect_true(all.equal(cases[[3]], cases[[4]]))
})

test_that("sample_or_maximise_from_gamma returns either draws or maximum of
          gamma", {
  shape <- 5
  rate <- 5
  ndraws <- 1
  val <- sample_or_maximise_gamma(shape, rate, ndraws)
  expect_equal(length(val), ndraws)
  ndraws <- 20
  vals <- sample_or_maximise_gamma(shape, rate, ndraws)
  expect_equal(length(vals), ndraws)

  # maximise
  val <- sample_or_maximise_gamma(shape, rate, ndraws,
                                  maximise=T)
  expect_equal(val, (shape - 1) / rate)
})

# tests for sampling Rt
days_total <- 100
Rt_1 <- 1.5
Rt_2 <- 1.0
Rt_3 <- 1.3
v_Rt <- c(rep(Rt_1, 40), rep(Rt_2, 20), rep(Rt_3, 40))
Rt_function <- stats::approxfun(1:days_total, v_Rt)
s_params <- list(mean=5, sd=3)
r_params <- list(mean=10, sd=3)
kappa <- 1000
df <- generate_snapshots(days_total, Rt_function, s_params, r_params,
                         kappa=kappa) %>%
  dplyr::select(time_onset, cases_true)
Rt_indices <- unlist(purrr::map(seq(1, 5, 1), ~rep(., 20)))
Rt_index_lookup <- dplyr::tibble(
  time_onset=seq_along(Rt_indices),
  Rt_index=Rt_indices)
df <- df %>%
  dplyr::left_join(Rt_index_lookup, by = "time_onset") %>%
  dplyr::select(time_onset, cases_true, Rt_index) %>%
  unique()
Rt_prior <- list(shape=1, rate=1)

test_that("sample_Rt_single_piece returns reasonable values: these are
          basically functional tests", {

  ndraws <- 42
  Rt_vals <- sample_Rt_single_piece(
    2, df,
    Rt_prior, s_params,
    ndraws = ndraws,
    serial_max = 20)
  expect_equal(ndraws, length(Rt_vals))

  f_Rt <- function(piece) {
    sample_Rt_single_piece(
      piece, df,
      Rt_prior, s_params,
      ndraws = 1000,
      serial_max = 20)
  }

  Rt_vals <- sample_Rt_single_piece(
    2, df,
    Rt_prior, s_params,
    ndraws = 1000,
    serial_max = 40)

  # test for maximisation
  f_Rt <- function(i) sample_Rt_single_piece(
    2, df,
    Rt_prior, s_params,
    ndraws = 1000,
    serial_max = 20,
    maximise = T)
  expect_equal(length(f_Rt(1)), 1)
  Rt_vals1 <- purrr::map_dbl(seq(1, 3, 1), f_Rt)
  # no variation in maximum so should be no sd
  expect_equal(sd(Rt_vals1), 0)

  # mode of sampling distribution should be close to max
  mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  Rt_vals <- sample_Rt_single_piece(
    2, df,
    Rt_prior, s_params,
    ndraws = 1000,
    serial_max = 20)
  expect_true(abs(mode(Rt_vals) - Rt_vals[1]) < 0.2)
})

test_that("sample_Rt_single_piece works ok with NB renewal model", {
  ndraws <- 420
  Rt_vals <- sample_Rt_single_piece(
    2, df,
    Rt_prior, s_params,
    ndraws = ndraws,
    serial_max = 20)

  # for large kappa draws should be similar
  kappa <- 10000
  Rt_vals_nb <- sample_Rt_single_piece(
    2, df,
    Rt_prior, s_params,
    ndraws = ndraws,
    serial_max = 20,
    kappa=kappa,
    is_negative_binomial=TRUE)
  expect_true(abs(mean(Rt_vals_nb) - mean(Rt_vals)) < 0.5)
  expect_true(abs(sd(Rt_vals_nb) - sd(Rt_vals)) < 0.1)

  # for small kappa the NB sd should exceed Poisson
  kappa <- 0.5
  Rt_vals_nb <- sample_Rt_single_piece(
    2, df,
    Rt_prior, s_params,
    ndraws = ndraws,
    serial_max = 20,
    kappa=kappa,
    is_negative_binomial=TRUE)
  expect_true(sd(Rt_vals_nb) > sd(Rt_vals))
})

test_that("sample_Rt returns sensible values", {
  ndraws <- 300
  Rt_df <- sample_Rt(df,
                     Rt_prior, s_params,
                     ndraws = ndraws,
                     serial_max = 20)

  # check shapes of outputs
  cnames <- colnames(Rt_df)
  expect_true(all.equal(cnames,
                        c("Rt_index", "draw_index", "Rt")))
  npieces <- max(df$Rt_index)
  expect_equal(npieces * ndraws, nrow(Rt_df))

  # check substance of outputs
  f_Rt <- function(piece) {
    val <- sample_Rt_single_piece(
      piece, df,
      Rt_prior, s_params,
      ndraws = 1000,
      serial_max = 20)
    mean(val)
  }
  indices <- seq(1, 5, 1)
  Rt_vals <- purrr::map_dbl(indices, f_Rt)
  single_df <- dplyr::tibble(Rt_index=indices,
                             Rt_single=Rt_vals)
  Rt_df <- Rt_df %>%
    dplyr::group_by(.data$Rt_index) %>%
    dplyr::summarise(Rt=mean(Rt)) %>%
    dplyr::left_join(single_df, by = "Rt_index") %>%
    dplyr::mutate(diff=Rt-Rt_single)
  expect_true(abs(sum(Rt_df$diff)) < 0.2)

  # check maximisation
  f_Rt <- function(piece) {
    val <- sample_Rt_single_piece(
      piece, df,
      Rt_prior, s_params,
      serial_max = 20,
      maximise = T)
    mean(val)
  }
  indices <- seq(1, 5, 1)
  Rt_vals <- purrr::map_dbl(indices, f_Rt)
  single_df <- dplyr::tibble(Rt_index=indices,
                             Rt_single=Rt_vals)
  Rt_df <- sample_Rt(df,
                     Rt_prior, s_params,
                     ndraws = ndraws,
                     serial_max = 20,
                     maximise = T)
  expect_equal(nrow(Rt_df), max(df$Rt_index))
  Rt_df <- Rt_df %>%
    dplyr::left_join(single_df, by = "Rt_index") %>%
    dplyr::mutate(diff=Rt-Rt_single)
  expect_true(abs(sum(Rt_df$diff)) < 0.0001)
})

test_that("sample_Rt works for negative binomial renewal model", {

  ndraws <- 300
  Rt_df <- sample_Rt(df,
                     Rt_prior, s_params,
                     ndraws = ndraws,
                     serial_max = 20)

  Rt_df_nb <- sample_Rt(df,
                     Rt_prior, s_params,
                     ndraws = ndraws,
                     serial_max = 20,
                     kappa=1,
                     is_negative_binomial=TRUE)
  expect_equal(nrow(Rt_df), nrow(Rt_df_nb))
  expect_equal(ncol(Rt_df), ncol(Rt_df_nb))

  # errors thrown if not providing kappa
  expect_error(sample_Rt(df,
                         Rt_prior, s_params,
                         ndraws = ndraws,
                         serial_max = 20,
                         is_negative_binomial=TRUE))

  # errors thrown if providing invalid kappa
  expect_error(sample_Rt(df,
                         Rt_prior, s_params,
                         ndraws = ndraws,
                         serial_max = 20,
                         kappa=0,
                         is_negative_binomial=TRUE))
})

test_that("accept_reject works ok", {

  logp_current <- 0
  logp_proposed <- -10000
  current_parameters <- 0
  proposed_parameters <- 1
  list_param_logp <- accept_reject(logp_current, logp_proposed,
    current_parameters,
    proposed_parameters)
  expect_equal(current_parameters, list_param_logp$parameter)
  expect_equal(logp_current, list_param_logp$logp)
})

test_that("prior_reporting_parameters works ok", {
  current_reporting_parameters <- list(mean=3, sd=2, reporting_piece_index=1)
  prior_params <- list(mean_mu=3, mean_sigma=2,
                       sd_mu=5, sd_sigma=1)
  val <- prior_reporting_parameters(current_reporting_parameters,
                                    prior_params)
  val1 <- dgamma_mean_sd(current_reporting_parameters$mean,
                         prior_params$mean_mu,
                         prior_params$mean_sigma,
                         log=TRUE)
  val2 <- dgamma_mean_sd(current_reporting_parameters$sd,
                         prior_params$sd_mu,
                         prior_params$sd_sigma,
                         log=TRUE)
  expect_equal(val, val1 + val2)
})

test_that("propose_reporting_parameters works ok", {
  r_params <- dplyr::tibble(mean=4, sd=4, reporting_piece_index=1)
  met_params <- list(mean_step=0.1, sd_step=0.2)
  r_params1 <- propose_reporting_parameters(
    r_params, met_params)
  expect_true(abs(r_params$mean - r_params1$mean) < 3)
  expect_true(abs(r_params$sd - r_params1$sd) < 3)

  met_params <- list(mean_step=0.000001, sd_step=0.2)
  r_params1 <- propose_reporting_parameters(
    r_params, met_params)
  expect_true(abs(r_params$mean - r_params1$mean) < 0.1)
})

days_total <- 30
df <- generate_snapshots(days_total, Rt_function,
                         s_params, r_params,
                         kappa=kappa, thinned=T) %>%
  dplyr::mutate(reporting_piece_index=1)

test_that("metropolis_step works as expected", {
  met_params <- list(mean_step=0.01, sd_step=0.01)
  prior_params <- list(mean_mu=2, mean_sigma=100,
                       sd_mu=2, sd_sigma=100)
  r_params <- dplyr::tibble(mean=r_params$mean,
                     sd=r_params$sd,
                     reporting_piece_index=1)
  logp_current <- -1000
  list_param_logp <- metropolis_step(
    df, r_params, logp_current,
    prior_params, met_params)
  r_params1 <- list_param_logp$parameter
  expect_true(abs(r_params1$mean - r_params$mean) < 0.2)
  expect_true(abs(r_params1$sd - r_params$sd) < 0.2)
})

test_that("metropolis_steps returns multiple steps", {

  met_params <- list(mean_step=0.01, sd_step=0.01)
  ndraws <- 10
  rep_prior_params <- list(mean_mu=5, sd_mu=3,
                           mean_sigma=5, sd_sigma=3)
  r_params <- dplyr::tibble(mean=r_params$mean,
                     sd=r_params$sd,
                     reporting_piece_index=1)
  logp_current <- -1000
  list_draws_logp <- metropolis_steps(
    snapshot_with_true_cases_df=df,
    current_reporting_parameters=r_params,
    logp_current=logp_current,
    prior_parameters=rep_prior_params,
    metropolis_parameters=met_params,
    ndraws=ndraws)
  expect_equal(names(list_draws_logp)[2], "logp")
  output <- list_draws_logp$reporting_parameters
  expect_equal(nrow(output), ndraws)
  expect_true(all.equal(colnames(output),
                        c("reporting_piece_index", "draw_index", "mean", "sd")))
})

test_that("propose_overdispersion_parameter works ok", {
  current <- 0
  sd <- 2
  new <- propose_overdispersion_parameter(current, sd)
  expect_true(current != new)
})

test_that("metropolis_step_overdispersion works ok", {

  overdispersion_current <- 100
  logp_current <- -100000
  cases_history_rt_df <- df %>%
    dplyr::select("time_onset", "cases_true") %>%
    unique() %>%
    dplyr::mutate(Rt=2)

  list_param_logp <- metropolis_step_overdispersion(
    overdispersion_current,
    logp_current,
    cases_history_rt_df,
    serial_parameters=list(mean=2, sd=3),
    prior_overdispersion_parameter=list(mean=10, sd=100),
    overdispersion_metropolis_sd=1)
  expect_true(list_param_logp$overdispersion != overdispersion_current)
  expect_true(list_param_logp$logp > logp_current)
})

test_that("metropolis_step_overdispersion uses log of prior density", {

  cases_history_rt_df <- df %>%
    dplyr::select("time_onset", "cases_true") %>%
    unique() %>%
    dplyr::mutate(Rt=2)
  s_params <- list(mean=2, sd=3)
  prior <- list(mean=10, sd=100)

  # always accepts
  list_param_logp <- metropolis_step_overdispersion(
    overdispersion_current=100,
    logp_current=-Inf,
    cases_history_rt_df,
    serial_parameters=s_params,
    prior_overdispersion_parameter=prior,
    overdispersion_metropolis_sd=1)
  overdispersion <- list_param_logp$overdispersion
  logp_expected <- state_process_nb_logp_all_onsets(
    overdispersion, cases_history_rt_df, s_params) +
    dgamma_mean_sd(overdispersion, prior$mean, prior$sd, log=TRUE)
  expect_equal(list_param_logp$logp, logp_expected)
})

test_that("maximise_reporting_logp maximises prob", {
  prior_params <- list(mean_mu=2, mean_sigma=100,
                       sd_mu=2, sd_sigma=100)
  r_params <- dplyr::tibble(reporting_piece_index=1,
                            mean=r_params$mean,
                            sd=r_params$sd)
  df <- df %>%
    dplyr::mutate(reporting_piece_index=1)
  output <- maximise_reporting_logp(df, r_params, prior_params)
  expect_true(abs(output$mean - r_params$mean) < 0.4)
  expect_true(abs(output$sd - r_params$sd) < 0.4)
  expect_equal(nrow(output), 1)
})

test_that("maximise_reporting_logp throws errors with piece info missing", {

  df <- df %>%
    dplyr::select(-"reporting_piece_index")
  prior_params <- list(mean_mu=2, mean_sigma=100,
                       sd_mu=2, sd_sigma=100)

  # r_params doesn't contain reporting_piece_index
  df_temp <- df %>%
    dplyr::mutate(reporting_piece_index=1)
  expect_error(maximise_reporting_logp(df_temp, r_params, prior_params))

  # df doesn't contain reporting_piece_index
  r_params <- dplyr::tibble(reporting_piece_index=1,
                            mean=r_params$mean,
                            sd=r_params$sd)
  expect_error(maximise_reporting_logp(df, r_params, prior_params))

})

test_that("maximise_reporting_logp works ok with multiple pieces", {

  days_total <- 30
  r_params <- dplyr::tibble(mean=c(rep(10, 15), rep(3, 15)),
                     sd=c(rep(3, 15), rep(1, 15)),
                     time_onset=seq_along(mean)) %>%
    dplyr::mutate(reporting_piece_index=c(rep(1, 15), rep(2, 15)))
  df <- generate_snapshots(days_total, Rt_function,
                           s_params, r_params,
                           kappa=kappa, thinned=T)
  r_params_short <- r_params %>%
    dplyr::select(time_onset, reporting_piece_index) %>%
    unique()

  df <- df %>%
    dplyr::left_join(r_params_short, by="time_onset")
  r_params <- dplyr::tibble(reporting_piece_index=c(1, 2),
                     mean=c(8, 4),
                     sd=c(4, 2))
  prior_params <- list(mean_mu=2, mean_sigma=100,
                       sd_mu=2, sd_sigma=100)
  output <- maximise_reporting_logp(df, r_params, prior_params)

  # try optimisation on each subperiod: should yield same
  r_params_1 <- r_params %>%
    dplyr::filter(reporting_piece_index==1)
  df_1 <- df %>%
    dplyr::filter(reporting_piece_index==1)
  output_1 <- maximise_reporting_logp(df_1, r_params_1, prior_params)
  r_params_2 <- r_params %>%
    dplyr::filter(reporting_piece_index==2)
  df_2 <- df %>%
    dplyr::filter(reporting_piece_index==2)
  output_2 <- maximise_reporting_logp(df_2, r_params_2, prior_params)

  expect_equal(output$mean[1], output_1$mean[1])
  expect_equal(output$mean[2], output_2$mean[1])
  expect_equal(output$sd[1], output_1$sd[1])
  expect_equal(output$sd[2], output_2$sd[1])
  expect_equal(output$reporting_piece_index[1], output_1$reporting_piece_index[1])
  expect_equal(output$reporting_piece_index[2], output_2$reporting_piece_index[1])
})

test_that("sample_reporting produces output of correct shape", {

  prior_params <- list(mean_mu=2, mean_sigma=100,
                       sd_mu=2, sd_sigma=100)
  met_params <- list(mean_step=0.01, sd_step=0.01)
  r_params <- dplyr::tibble(mean=r_params$mean,
                     sd=r_params$sd,
                     reporting_piece_index=1)
  logp_current <- -1000
  result <- sample_reporting(df, r_params, logp_current, prior_params, met_params)
  expect_equal(names(result)[2], "logp")
  output <- result$reporting_parameters
  expect_equal(nrow(output), 1)
  expect_equal(max(output$draw_index), 1)

  ndraws <- 2
  result <- sample_reporting(df, r_params, logp_current, prior_params, met_params,
                             ndraws=ndraws)
  output <- result$reporting_parameters
  expect_equal(nrow(output), ndraws)
  expect_equal(max(output$draw_index), ndraws)

  result <- sample_reporting(df, r_params, logp_current, prior_params, met_params,
                             maximise=T)
  output <- result$reporting_parameters
  expect_equal(nrow(output), 1)
  expect_equal(max(output$draw_index), 1)
})


test_that("mcmc produces outputs of correct shape", {

  niter <- 2
  days_total <- 100
  r_params <- list(mean=10, sd=3)
  s_params <- list(mean=5, sd=3)
  v_Rt <- c(rep(1.5, 40), rep(0.4, 20), rep(1.5, 40))
  Rt_function <- stats::approxfun(1:days_total, v_Rt)
  Rt_prior <- list(shape=1, rate=1)
  kappa <- 10
  df <- generate_snapshots(days_total, Rt_function, s_params, r_params,
                           kappa=kappa)

  snapshot_with_Rt_index_df <- df
  initial_cases_true <- df %>%
    dplyr::select(time_onset, cases_true) %>%
    unique()
  snapshot_with_Rt_index_df <- snapshot_with_Rt_index_df %>%
    dplyr::select(-cases_true)
  initial_Rt <- tidyr::tribble(~Rt_index, ~Rt,
                        1, 1.5,
                        2, 1.5,
                        3, 0.4,
                        4, 1.5,
                        5, 1.5)
  Rt_indices <- unlist(purrr::map(seq(1, 5, 1), ~rep(., 20)))

  Rt_index_lookup <- tidyr::tibble(
    time_onset=seq_along(Rt_indices),
    Rt_index=Rt_indices)
  snapshot_with_Rt_index_df <- snapshot_with_Rt_index_df %>%
    dplyr::left_join(Rt_index_lookup, by="time_onset") %>%
    dplyr::mutate(reporting_piece_index=1)

  initial_reporting_parameters <- dplyr::tibble(mean=5, sd=3, reporting_piece_index=1)
  serial_parameters <- list(mean=5, sd=3)
  priors <- list(Rt=Rt_prior,
                 reporting=list(mean_mu=5,
                                mean_sigma=10,
                                sd_mu=3,
                                sd_sigma=5),
                 max_cases=5000)

  # test throws errors
  ## wrongly named cols
  wrong_df <- snapshot_with_Rt_index_df %>%
    dplyr::rename(time_onset_wrong=time_onset)
  expect_error(mcmc(niterations=niter,
                    wrong_df,
                    priors,
                    serial_parameters,
                    initial_cases_true,
                    initial_reporting_parameters,
                    initial_Rt,
                    reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
                    serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE))

  ## too many cols
  wrong_df <- snapshot_with_Rt_index_df %>%
    dplyr::mutate(unnecessary_col="hi")
  expect_error(mcmc(niterations=niter,
                    wrong_df,
                    priors,
                    serial_parameters,
                    initial_cases_true,
                    initial_reporting_parameters,
                    initial_Rt,
                    reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
                    serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE))

  # five columns but not reporting_piece_index
  wrong_df <- snapshot_with_Rt_index_df %>%
    dplyr::rename(reporting_piece_index_wrong=reporting_piece_index)
  expect_error(mcmc(niterations=niter,
                    wrong_df,
                    priors,
                    serial_parameters,
                    initial_cases_true,
                    initial_reporting_parameters,
                    initial_Rt,
                    reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
                    serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE))

  # overdispersion parameter <= 0
  expect_error(mcmc(niterations=niter,
                   snapshot_with_Rt_index_df,
                   priors,
                   serial_parameters,
                   initial_cases_true,
                   initial_reporting_parameters,
                   initial_Rt,
                   reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
                   serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE,
                   initial_overdispersion = -1))

  # maximisation with a NB model doesn't work
  priors_tmp <- priors
  priors_tmp$overdispersion <- list(mean=5, sd=10)
  expect_error(mcmc(niterations=niter,
                    snapshot_with_Rt_index_df,
                    priors_tmp,
                    serial_parameters,
                    initial_cases_true,
                    initial_reporting_parameters,
                    initial_Rt,
                    reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
                    serial_max=40, p_gamma_cutoff=0.99, maximise=TRUE,
                    is_negative_binomial=TRUE))


  # test MCMC sampling
  niter <- 5
  res <- mcmc(niterations=niter,
              snapshot_with_Rt_index_df,
              priors,
              serial_parameters,
              initial_cases_true,
              initial_reporting_parameters,
              initial_Rt,
              reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
              serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE)

  ## check outputs
  ### overall
  expect_equal(length(res), 3)

  ### cases
  cases_df <- res$cases
  expect_true(all.equal(c("time_onset", "cases_true", "iteration", "chain"),
                        colnames(cases_df)))
  expect_equal(min(cases_df$iteration), 1)
  expect_equal(max(cases_df$iteration), niter)
  expect_equal(min(cases_df$time_onset), min(df$time_onset))
  expect_equal(max(cases_df$time_onset), max(df$time_onset))
  expect_equal(min(cases_df$chain), 1)
  expect_equal(max(cases_df$chain), 1)

  ## Rt
  rt_df <- res$Rt
  expect_true(all.equal(c("iteration", "Rt_index", "Rt", "chain"),
                        colnames(rt_df)))
  expect_equal(min(rt_df$iteration), 1)
  expect_equal(max(rt_df$iteration), niter)
  expect_equal(min(rt_df$Rt_index), min(initial_Rt$Rt_index))
  expect_equal(max(rt_df$Rt_index), max(initial_Rt$Rt_index))
  expect_equal(min(rt_df$chain), 1)
  expect_equal(max(rt_df$chain), 1)

  # reporting delays
  reporting_df <- res$reporting
  expect_true(all.equal(c("reporting_piece_index", "mean", "sd", "iteration", "chain"),
                        colnames(reporting_df)))
  expect_equal(min(reporting_df$iteration), 1)
  expect_equal(max(reporting_df$iteration), niter)
  expect_equal(min(reporting_df$chain), 1)
  expect_equal(max(reporting_df$chain), 1)

  # test optimisation
  niter <- 2 # needed since maximisation is iterative
  res <- mcmc(niterations=niter,
              snapshot_with_Rt_index_df,
              priors,
              serial_parameters,
              initial_cases_true,
              initial_reporting_parameters,
              initial_Rt,
              reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
              serial_max=40, p_gamma_cutoff=0.99, maximise=TRUE)
  ### overall
  expect_equal(length(res), 3)

  ### cases
  cases_df <- res$cases
  expect_true(all.equal(c("time_onset", "cases_true", "iteration", "chain"),
                        colnames(cases_df)))
  expect_equal(min(cases_df$iteration), 1)
  expect_equal(max(cases_df$iteration), niter)
  expect_equal(min(cases_df$time_onset), min(df$time_onset))
  expect_equal(max(cases_df$time_onset), max(df$time_onset))
  expect_equal(min(cases_df$chain), 1)
  expect_equal(max(cases_df$chain), 1)

  ## Rt
  rt_df <- res$Rt
  expect_true(all.equal(c("iteration", "Rt_index", "Rt", "chain"),
                        colnames(rt_df)))
  expect_equal(min(rt_df$iteration), 1)
  expect_equal(max(rt_df$iteration), niter)
  expect_equal(min(rt_df$Rt_index), min(initial_Rt$Rt_index))
  expect_equal(max(rt_df$Rt_index), max(initial_Rt$Rt_index))
  expect_equal(min(rt_df$chain), 1)
  expect_equal(max(rt_df$chain), 1)

  # reporting delays
  reporting_df <- res$reporting
  expect_true(all.equal(c("reporting_piece_index", "mean", "sd", "iteration", "chain"),
                        colnames(reporting_df)))
  expect_equal(min(reporting_df$iteration), 1)
  expect_equal(max(reporting_df$iteration), niter)
  expect_equal(min(reporting_df$chain), 1)
  expect_equal(max(reporting_df$chain), 1)


  # test MCMC sampling with NB model
  niter <- 5

  ## prior not provided on overdispersion parameter
  expect_error(mcmc(niterations=niter,
              snapshot_with_Rt_index_df,
              priors,
              serial_parameters,
              initial_cases_true,
              initial_reporting_parameters,
              initial_Rt,
              reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
              serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE,
              is_negative_binomial=TRUE))

  # sampling works if given overdispersion prior
  priors$overdispersion <- list(mean=5, sd=2)
  res <- mcmc(niterations=niter,
              snapshot_with_Rt_index_df,
              priors,
              serial_parameters,
              initial_cases_true,
              initial_reporting_parameters,
              initial_Rt,
              reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
              serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE,
              is_negative_binomial=TRUE)
  expect_equal(length(res), 4)
  expect_true("overdispersion" %in% names(res))
  overdispersion <- res$overdispersion
  expect_equal(nrow(overdispersion), 5)
  expect_true("overdispersion" %in% names(overdispersion))
  expect_true("iteration" %in% names(overdispersion))
})

test_that("multiple chains works", {

  days_total <- 100
  r_params <- list(mean=10, sd=3)
  s_params <- list(mean=5, sd=3)
  v_Rt <- c(rep(1.5, 40), rep(0.4, 20), rep(1.5, 40))
  Rt_function <- stats::approxfun(1:days_total, v_Rt)
  Rt_prior <- list(shape=1, rate=1)
  kappa <- 10
  df <- generate_snapshots(days_total, Rt_function, s_params, r_params,
                           kappa=kappa)

  snapshot_with_Rt_index_df <- df
  initial_cases_true <- df %>%
    dplyr::select(time_onset, cases_true) %>%
    unique()
  snapshot_with_Rt_index_df <- snapshot_with_Rt_index_df %>%
    dplyr::select(-cases_true)
  initial_Rt <- tidyr::tribble(~Rt_index, ~Rt,
                               1, 1.5,
                               2, 1.5,
                               3, 0.4,
                               4, 1.5,
                               5, 1.5)
  Rt_indices <- unlist(purrr::map(seq(1, 5, 1), ~rep(., 20)))

  Rt_index_lookup <- tidyr::tibble(
    time_onset=seq_along(Rt_indices),
    Rt_index=Rt_indices)
  snapshot_with_Rt_index_df <- snapshot_with_Rt_index_df %>%
    dplyr::left_join(Rt_index_lookup)

  initial_reporting_parameters <- list(mean=5, sd=3)
  serial_parameters <- list(mean=5, sd=3)
  priors <- list(Rt=Rt_prior,
                 reporting=list(mean_mu=5,
                                mean_sigma=10,
                                sd_mu=3,
                                sd_sigma=5),
                 max_cases=5000)

  # multiple chains in serial
  niter <- 2
  nchains <- 2
  res <- mcmc(niterations=niter,
              snapshot_with_Rt_index_df,
              priors,
              serial_parameters,
              initial_cases_true,
              initial_reporting_parameters,
              initial_Rt,
              reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
              serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE,
              nchains=nchains)
  cases_df <- res$cases
  Rt_df <- res$Rt
  rep_df <- res$reporting

  expect_equal(max(cases_df$chain), nchains)
  expect_equal(max(Rt_df$chain), nchains)
  expect_equal(max(rep_df$chain), nchains)

  nchains <- 3
  res <- mcmc(niterations=niter,
              snapshot_with_Rt_index_df,
              priors,
              serial_parameters,
              initial_cases_true,
              initial_reporting_parameters,
              initial_Rt,
              reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
              serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE,
              nchains=nchains)
  cases_df <- res$cases
  Rt_df <- res$Rt
  rep_df <- res$reporting

  expect_equal(max(cases_df$chain), nchains)
  expect_equal(max(Rt_df$chain), nchains)
  expect_equal(max(rep_df$chain), nchains)

  # multiple chains in parallel
  library(doParallel)
  cl <- makeCluster(2)
  registerDoParallel(cl)
  res <- mcmc(niterations=niter,
              snapshot_with_Rt_index_df,
              priors,
              serial_parameters,
              initial_cases_true,
              initial_reporting_parameters,
              initial_Rt,
              reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
              serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE,
              nchains=nchains, is_parallel=TRUE, is_negative_binomial=FALSE)
  stopCluster(cl)

  expect_equal(max(cases_df$chain), nchains)
  expect_equal(max(Rt_df$chain), nchains)
  expect_equal(max(rep_df$chain), nchains)

  # multiple chains in parallel using NB model
  cl <- makeCluster(2)
  registerDoParallel(cl)
  priors$overdispersion <- list(mean=5, sd=10)
  res <- mcmc(niterations=niter,
              snapshot_with_Rt_index_df,
              priors,
              serial_parameters,
              initial_cases_true,
              initial_reporting_parameters,
              initial_Rt,
              reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
              serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE,
              nchains=nchains, is_parallel=TRUE, is_negative_binomial=TRUE)
  stopCluster(cl)

  cases_df <- res$cases
  Rt_df <- res$Rt
  rep_df <- res$reporting
  overdispersion_df <- res$overdispersion

  expect_equal(max(cases_df$chain), nchains)
  expect_equal(max(Rt_df$chain), nchains)
  expect_equal(max(rep_df$chain), nchains)
  expect_equal(max(overdispersion_df$chain), nchains)
})

test_that("construct_w_matrix works ok", {
  w <- c(0.5, 0.25, 0.25)
  piece_width <- 5
  m_w <- construct_w_matrix(w, piece_width)

  wmax <- length(w)
  expect_equal(nrow(m_w), piece_width)
  expect_equal(ncol(m_w), wmax + piece_width - 1)

  row_sums <- rowSums(m_w)
  expect_true(sum(row_sums==1) == piece_width)
})

test_that("nb_log_likelihood_Rt_piece works as expected", {

  w <- c(0.5, 0.25, 0.25)
  Rt <- 2.2
  kappa <- 3.5
  onset_times <- c(4, 5)
  cases_df <- dplyr::tribble(
    ~time_onset, ~cases_true,
    1, 1,
    2, 3,
    3, 6,
    4, 5,
    5, 3
  )

  # check throws error if cases_df doesn't have onset times at bottom
  cases_df_wrong <- cases_df %>%
    dplyr::arrange(desc(time_onset))
  expect_error(nb_log_likelihood_Rt_piece(
    Rt, kappa, w, onset_times, cases_df_wrong)
  )

  # check works out logp pointwise
  logp <- nb_log_likelihood_Rt_piece(
    Rt, kappa, w, onset_times, cases_df
  )
  mu <- sum(Rt * w * rev(cases_df$cases_true[1:3]))
  logp_1 <- dnbinom(cases_df$cases_true[4], mu=mu, size=kappa,
                    log = TRUE)
  mu <- sum(Rt * w * rev(cases_df$cases_true[2:4]))
  logp_2 <- dnbinom(cases_df$cases_true[5], mu=mu, size=kappa,
                    log = TRUE)
  expect_equal(logp_1 + logp_2, logp)

  # check onset_time 1 contributes no logp
  onset_times <- c(1, 2)
  cases_df <- cases_df %>%
    dplyr::filter(time_onset %in% onset_times)
  logp <- nb_log_likelihood_Rt_piece(
    Rt, kappa, w, onset_times, cases_df
  )
  mu <- sum(Rt * w * c(cases_df$cases_true[1], 0, 0))
  logp_2 <- dnbinom(cases_df$cases_true[2], mu=mu, size=kappa,
                    log = TRUE)
  expect_equal(logp, logp_2)

})

test_that("sample_nb_Rt_piece works as expected", {

  prior_shape <- 1
  prior_rate <- 0.2

  w <- c(0.5, 0.25, 0.25)
  onset_times <- c(4, 5)
  cases_df <- dplyr::tribble(
    ~time_onset, ~cases_true,
    1, 10,
    2, 30,
    3, 60,
    4, 50,
    5, 30
  )

  posterior_shape <- prior_shape + sum(cases_df$cases_true[4:5])
  posterior_rate <- prior_rate + sum(w * rev(cases_df$cases_true[1:3])) +
    sum(w * rev(cases_df$cases_true[2:4]))

  ndraws <- 3000
  set.seed(1)

  # for large kappa the NB and Poisson distributions should be similar
  kappa <- 10000
  Rt_nb <- sample_nb_Rt_piece(
    prior_shape, prior_rate,
    posterior_shape, posterior_rate,
    kappa,
    w,
    onset_times,
    cases_df,
    ndraws)
  expect_equal(length(Rt_nb), ndraws)
  expect_equal(mean(Rt_nb), posterior_shape / posterior_rate, tolerance=0.02)
  expect_equal(sd(Rt_nb), sqrt(posterior_shape) / posterior_rate, tolerance=0.1)

  # for small kappa draws follow the conditional distribution found numerically
  kappa <- 1
  f_density <- function(Rt, power=0) {
    purrr::map_dbl(Rt, function(r) {
      r^power * exp(nb_log_likelihood_Rt_piece(r, kappa, w, onset_times, cases_df)) *
        stats::dgamma(r, prior_shape, prior_rate)
    })
  }
  constant <- integrate(f_density, 0, Inf)$value
  mean_exact <- integrate(f_density, 0, Inf, power=1)$value / constant
  sd_exact <- sqrt(integrate(f_density, 0, Inf, power=2)$value / constant -
                     mean_exact^2)
  Rt_nb <- sample_nb_Rt_piece(
    prior_shape, prior_rate,
    posterior_shape, posterior_rate,
    kappa,
    w,
    onset_times,
    cases_df,
    ndraws)
  expect_equal(mean(Rt_nb), mean_exact, tolerance=0.05)
  expect_equal(sd(Rt_nb), sd_exact, tolerance=0.1)
  expect_true(sd(Rt_nb) > 2 * sqrt(posterior_shape) / posterior_rate)
  probs <- c(0.1, 0.5, 0.9)
  f_cdf <- function(q) integrate(f_density, 0, q)$value / constant
  expect_equal(purrr::map_dbl(stats::quantile(Rt_nb, probs), f_cdf),
               probs, tolerance=0.1, ignore_attr=TRUE)

  # a single step from the current value either stays there or moves
  Rt_one <- sample_nb_Rt_piece(
    prior_shape, prior_rate, posterior_shape, posterior_rate, kappa, w,
    onset_times, cases_df, ndraws=1, nsteps=0, Rt_current=0.7)
  expect_equal(Rt_one, 0.7)
  # the starting value is forgotten, even if it is very unlikely or invalid
  for(Rt_current in list(500, 1e-6, -1, NA, NULL)) {
    Rt_far <- sample_nb_Rt_piece(
      prior_shape, prior_rate, posterior_shape, posterior_rate, kappa, w,
      onset_times, cases_df, ndraws=200, nsteps=10, Rt_current=Rt_current)
    expect_true(all(is.finite(Rt_far) & Rt_far > 0))
    expect_true(abs(median(Rt_far[-(1:20)]) - mean_exact) < 2 * sd_exact)
  }

  # very strong overdispersion does not cause failure
  Rt_small <- sample_nb_Rt_piece(
    prior_shape, prior_rate, posterior_shape, posterior_rate, 0.001, w,
    onset_times, cases_df, ndraws=50)
  expect_true(all(is.finite(Rt_small) & Rt_small > 0))
})

test_that("sample_Rt starts negative binomial updates from current values", {
  set.seed(2)
  cases_df <- dplyr::tibble(time_onset=1:20,
                            cases_true=c(10, rpois(19, 30)),
                            Rt_index=rep(1:2, each=10))
  prior <- list(shape=1, rate=1)
  s_params <- list(mean=4, sd=2)
  # with no steps taken the current values are returned
  Rt_df <- sample_Rt(cases_df, prior, s_params, kappa=5, serial_max=10,
                     is_negative_binomial=TRUE, nsteps=0,
                     Rt_current=dplyr::tibble(Rt_index=c(2, 1), Rt=c(0.6, 1.7)))
  expect_equal(Rt_df$Rt[Rt_df$Rt_index == 1], 1.7)
  expect_equal(Rt_df$Rt[Rt_df$Rt_index == 2], 0.6)
  Rt_df <- sample_Rt(cases_df, prior, s_params, kappa=5, serial_max=10,
                     is_negative_binomial=TRUE,
                     Rt_current=dplyr::tibble(Rt_index=c(2, 1), Rt=c(0.6, 1.7)))
  expect_true(all(Rt_df$Rt > 0))
})

test_that("sample_true_cases_single_onset uses later cases if supplied", {
  observation_matrix <- dplyr::tibble(time_reported=c(6, 7),
                                      cases_reported=c(1, 2))
  reporting_parameters <- list(mean=10, sd=3)
  s_params <- list(mean=2, sd=1)
  cases_history <- rep(5, 5)
  max_cases <- 200
  f_map <- function(cases_future) {
    sample_true_cases_single_onset(
      observation_df=observation_matrix,
      cases_history=cases_history,
      max_cases=max_cases,
      Rt=1,
      day_onset=6,
      serial_parameters=s_params,
      reporting_parameters=reporting_parameters,
      maximise=TRUE,
      cases_future=cases_future,
      Rt_future=rep(1, length(cases_future)))
  }
  # many later cases imply more cases on onset day than few later cases
  expect_true(f_map(c(60, 60)) > f_map(c(5, 5)))
  expect_true(f_map(c(5, 5)) > f_map(c(0, 0)))
})

test_that("sample_cases_history conditioning on later cases targets joint
posterior of case counts", {
  set.seed(11)
  days_total <- 8
  s_params <- list(mean=2, sd=1)
  r_list <- list(mean=3, sd=2)
  Rt_function <- function(t) rep(1.4, length(t))
  df <- generate_snapshots(days_total, Rt_function, s_params, r_list,
                           initial_parameters=list(mean=6, length=5)) %>%
    dplyr::filter(time_reported <= days_total) %>%
    dplyr::select(-cases_true) %>%
    dplyr::mutate(reporting_piece_index=1)
  r_params <- dplyr::tibble(reporting_piece_index=1, mean=3, sd=2)
  max_cases <- 18
  p_gamma_cutoff <- 0.55

  onset_times <- sort(unique(df$time_onset))
  max_reported <- purrr::map_dbl(
    onset_times, ~max(df$cases_reported[df$time_onset==.]))
  days_uncertain <- onset_times[
    onset_times >= days_total - max_uncertain_days(p_gamma_cutoff, r_list)]
  expect_equal(days_uncertain, c(6, 7, 8))

  # exact joint posterior of uncertain case counts by enumeration
  observations <- purrr::map(
    days_uncertain,
    ~df[df$time_onset==., c("time_reported", "cases_reported")])
  cases_fixed <- max_reported[onset_times < min(days_uncertain)]
  grid <- as.matrix(expand.grid(
    purrr::map(days_uncertain, ~max_reported[onset_times==.]:max_cases)))
  logp_joint <- rep(0, nrow(grid))
  for(r in seq_len(nrow(grid))) {
    for(k in seq_along(days_uncertain)) {
      day <- days_uncertain[k]
      cases_earlier <- c(cases_fixed, grid[r, seq_len(k - 1)])
      if(nrow(observations[[k]]) > 1)
        logp_joint[r] <- logp_joint[r] + observation_process_logp(
          observations[[k]], grid[r, k], day, r_list)
      logp_joint[r] <- logp_joint[r] + state_process_logp(
        grid[r, k], rev(cases_earlier), Rt_function(day), s_params)
    }
  }
  p_joint <- exp(logp_joint - matrixStats::logSumExp(logp_joint))
  mean_exact <- colSums(grid * p_joint)
  sd_exact <- sqrt(colSums(grid^2 * p_joint) - mean_exact^2)

  # Gibbs sweeps
  nsweeps <- 1000
  cases_current <- dplyr::tibble(time_onset=onset_times,
                                 cases_true=max_reported + 3)
  m_draws <- matrix(nrow=nsweeps, ncol=length(days_uncertain))
  for(i in 1:nsweeps) {
    # the upper limit on cases is kept small so that the exact posterior can
    # be enumerated, which the sampler warns about
    df_est <- suppressWarnings(sample_cases_history(
      df, max_cases, Rt_function, s_params, r_params,
      p_gamma_cutoff=p_gamma_cutoff,
      cases_true_current=cases_current))
    expect_equal(sum(is.na(df_est$cases_estimated)), 0)
    cases_current <- df_est %>%
      dplyr::select(time_onset, cases_true=cases_estimated) %>%
      unique()
    m_draws[i, ] <- cases_current$cases_true[
      match(days_uncertain, cases_current$time_onset)]
  }

  # certain days are unaffected by values supplied
  cases_certain <- cases_current$cases_true[onset_times < min(days_uncertain)]
  expect_equal(cases_certain, cases_fixed)

  # sampled means within Monte Carlo error of exact (allowing for autocorrelation)
  z_scores <- (colMeans(m_draws) - mean_exact) / (sd_exact / sqrt(nsweeps / 3))
  expect_true(all(abs(z_scores) < 3))
})

test_that("sample_cases_history throws errors for malformed current cases", {
  days_total <- 30
  r_params <- list(mean=10, sd=3)
  Rt_function <- stats::approxfun(1:days_total, rep(1.5, days_total))
  s_params <- list(mean=5, sd=3)
  df <- generate_snapshots(days_total, Rt_function,
                           s_params, r_params, thinned=T) %>%
    dplyr::mutate(reporting_piece_index=1)
  cases_current <- df %>%
    dplyr::select(time_onset, cases_true) %>%
    unique()
  df <- df %>% dplyr::select(-cases_true)
  r_params <- dplyr::tibble(mean=10, sd=3, reporting_piece_index=1)
  max_cases <- 5000

  # works fine
  df_est <- sample_cases_history(df, max_cases, Rt_function, s_params, r_params,
                                 cases_true_current=cases_current)
  expect_equal(sum(is.na(df_est$cases_estimated)), 0)

  # wrongly named columns
  expect_error(sample_cases_history(
    df, max_cases, Rt_function, s_params, r_params,
    cases_true_current=dplyr::rename(cases_current, cases=cases_true)))

  # missing onset times
  expect_error(sample_cases_history(
    df, max_cases, Rt_function, s_params, r_params,
    cases_true_current=cases_current[1:(days_total - 2), ]))
})

test_that("mcmc_single runs both with and without conditioning on later cases", {
  niter <- 2
  days_total <- 40
  r_params <- list(mean=5, sd=3)
  s_params <- list(mean=5, sd=3)
  Rt_function <- stats::approxfun(1:days_total, rep(1.3, days_total))
  df <- generate_snapshots(days_total, Rt_function, s_params, r_params)
  initial_cases_true <- df %>%
    dplyr::select(time_onset, cases_true) %>%
    unique()
  snapshot_with_Rt_index_df <- df %>%
    dplyr::select(-cases_true) %>%
    dplyr::mutate(Rt_index=ifelse(time_onset <= 20, 1, 2))
  initial_Rt <- dplyr::tibble(Rt_index=c(1, 2), Rt=c(1.3, 1.3))
  priors <- list(Rt=list(shape=1, rate=1),
                 reporting=list(mean_mu=5,
                                mean_sigma=10,
                                sd_mu=3,
                                sd_sigma=5),
                 max_cases=5000)

  for(condition_on_future_cases in c(TRUE, FALSE)) {
    res <- mcmc_single(niterations=niter,
                       snapshot_with_Rt_index_df,
                       priors,
                       s_params,
                       initial_cases_true,
                       r_params,
                       initial_Rt,
                       print_to_screen=FALSE,
                       condition_on_future_cases=condition_on_future_cases)
    expect_equal(nrow(res$cases), days_total * niter)
    expect_equal(sum(is.na(res$cases$cases_true)), 0)
    expect_equal(nrow(res$Rt), 2 * niter)
  }
})

test_that("sample_true_cases_single_onset allows true cases below max reported
if revision parameters supplied", {
  # reported cases long after onset with a short delay settle at three
  observation_matrix <- dplyr::tibble(time_reported=c(1, 3, 5, 20, 21, 22),
                                      cases_reported=c(0, 5, 5, 3, 3, 3))
  reporting_parameters <- list(mean=2, sd=1)
  s_params <- list(mean=2, sd=1)
  cases_history <- rep(4, 5)
  f_map <- function(revision_parameters) {
    sample_true_cases_single_onset(
      observation_df=observation_matrix,
      cases_history=cases_history,
      max_cases=100,
      Rt=1,
      day_onset=1,
      serial_parameters=s_params,
      reporting_parameters=reporting_parameters,
      maximise=TRUE,
      revision_parameters=revision_parameters)
  }
  expect_error(f_map(NULL))
  expect_equal(f_map(list(eta=0.1, rho=0.2)), 3)

  # draws never fall more than max_erroneous below max reported
  draws <- sample_true_cases_single_onset(
    observation_df=observation_matrix[1:3, ],
    cases_history=cases_history,
    max_cases=100,
    Rt=1,
    day_onset=1,
    serial_parameters=s_params,
    reporting_parameters=reporting_parameters,
    ndraws=500,
    revision_parameters=list(eta=0.5, rho=0.2),
    max_erroneous=2)
  expect_true(min(draws) >= 3)
  expect_true(any(draws < 5))
})

test_that("sample_cases_history uses latest reported cases for onset times
outside of uncertain period", {
  set.seed(3)
  days_total <- 40
  r_list <- list(mean=3, sd=2)
  revision_parameters <- list(eta=0.2, rho=0.3)
  Rt_function <- stats::approxfun(1:days_total, rep(1.2, days_total))
  s_params <- list(mean=5, sd=3)
  df <- generate_snapshots(days_total, Rt_function, s_params, r_list,
                           revision_parameters=revision_parameters) %>%
    dplyr::select(-cases_true) %>%
    dplyr::mutate(reporting_piece_index=1)
  r_params <- dplyr::tibble(reporting_piece_index=1, mean=3, sd=2)

  max_erroneous <- max_downward_revision(df) + 1
  df_est <- sample_cases_history(df, 2000, Rt_function, s_params, r_params,
                                 revision_parameters=revision_parameters,
                                 max_erroneous=max_erroneous)
  expect_equal(sum(is.na(df_est$cases_estimated)), 0)

  df_summary <- df_est %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(cases_latest=cases_reported[which.max(time_reported)],
                     cases_max=max(cases_reported),
                     cases_estimated=cases_estimated[1])
  days_uncertain <- max_uncertain_days(0.99, r_list)
  df_certain <- df_summary %>%
    dplyr::filter(time_onset < days_total - days_uncertain)
  expect_true(any(df_certain$cases_latest < df_certain$cases_max))
  expect_equal(df_certain$cases_estimated, df_certain$cases_latest)

  # without revision parameters these data cannot be handled
  expect_error(sample_cases_history(df, 2000, Rt_function, s_params, r_params))
})

test_that("sample_cases_history with downward revisions targets joint
posterior of case counts", {
  set.seed(15)
  days_total <- 8
  s_params <- list(mean=2, sd=1)
  r_list <- list(mean=3, sd=2)
  revision_parameters <- list(eta=0.3, rho=0.3)
  max_erroneous <- 2
  Rt_function <- function(t) rep(1.4, length(t))
  # reported cases for the three latest onset times include downward revisions
  df <- dplyr::tibble(
    time_onset=c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5,
                 6, 6, 6, 7, 7, 8),
    time_reported=c(1, 8, 2, 8, 3, 8, 4, 8, 5, 6, 8,
                    6, 7, 8, 7, 8, 8),
    cases_reported=c(0, 4, 0, 3, 0, 5, 0, 4, 0, 4, 3,
                     0, 3, 2, 0, 2, 0),
    reporting_piece_index=1)
  r_params <- dplyr::tibble(reporting_piece_index=1, mean=3, sd=2)
  max_cases <- 24
  p_gamma_cutoff <- 0.55

  onset_times <- sort(unique(df$time_onset))
  max_reported <- purrr::map_dbl(
    onset_times, ~max(df$cases_reported[df$time_onset==.]))
  latest_reported <- purrr::map_dbl(
    onset_times, ~dplyr::last(df$cases_reported[df$time_onset==.]))
  days_uncertain <- onset_times[
    onset_times >= days_total - max_uncertain_days(p_gamma_cutoff, r_list)]
  expect_equal(days_uncertain, c(6, 7, 8))

  # exact joint posterior of uncertain case counts by enumeration
  cases_fixed <- latest_reported[onset_times < min(days_uncertain)]
  l_possible <- purrr::map(
    days_uncertain,
    ~max(0, max_reported[onset_times==.] - max_erroneous):max_cases)
  l_logp_observation <- purrr::map(seq_along(days_uncertain), function(k) {
    observations <- df[df$time_onset==days_uncertain[k],
                       c("time_reported", "cases_reported")]
    if(nrow(observations) > 1)
      observation_process_revisions_logp(
        observations, l_possible[[k]], days_uncertain[k], r_list,
        revision_parameters, max_erroneous)
    else
      rep(0, length(l_possible[[k]]))
  })
  grid <- as.matrix(expand.grid(l_possible))
  grid_index <- as.matrix(expand.grid(purrr::map(l_possible, seq_along)))
  logp_joint <- rep(0, nrow(grid))
  for(r in seq_len(nrow(grid))) {
    for(k in seq_along(days_uncertain)) {
      cases_earlier <- c(cases_fixed, grid[r, seq_len(k - 1)])
      logp_joint[r] <- logp_joint[r] +
        l_logp_observation[[k]][grid_index[r, k]] +
        state_process_logp(grid[r, k], rev(cases_earlier),
                           Rt_function(days_uncertain[k]), s_params)
    }
  }
  p_joint <- exp(logp_joint - matrixStats::logSumExp(logp_joint))
  mean_exact <- colSums(grid * p_joint)
  sd_exact <- sqrt(colSums(grid^2 * p_joint) - mean_exact^2)

  # Gibbs sweeps
  nsweeps <- 1000
  cases_current <- dplyr::tibble(time_onset=onset_times,
                                 cases_true=max_reported + 3)
  m_draws <- matrix(nrow=nsweeps, ncol=length(days_uncertain))
  for(i in 1:nsweeps) {
    df_est <- sample_cases_history(
      df, max_cases, Rt_function, s_params, r_params,
      p_gamma_cutoff=p_gamma_cutoff,
      cases_true_current=cases_current,
      revision_parameters=revision_parameters,
      max_erroneous=max_erroneous)
    cases_current <- df_est %>%
      dplyr::select(time_onset, cases_true=cases_estimated) %>%
      unique()
    m_draws[i, ] <- cases_current$cases_true[
      match(days_uncertain, cases_current$time_onset)]
  }

  z_scores <- (colMeans(m_draws) - mean_exact) / (sd_exact / sqrt(nsweeps / 3))
  expect_true(all(abs(z_scores) < 3))
})

test_that("propose_revision_parameters proposes values within bounds", {
  current <- list(eta=0.01, rho=0.2)
  steps <- list(eta_step=2, rho_step=2)
  l_proposed <- purrr::map(1:200, ~propose_revision_parameters(current, steps))
  etas <- purrr::map_dbl(l_proposed, "eta")
  rhos <- purrr::map_dbl(l_proposed, "rho")
  expect_true(all(etas > 0))
  expect_true(all(rhos > 0 & rhos < 1))
  expect_true(length(unique(etas)) > 1)
  expect_true(length(unique(rhos)) > 1)

  # zero step size holds parameter fixed
  proposed <- propose_revision_parameters(current, list(eta_step=1, rho_step=0))
  expect_true(proposed$eta != current$eta)
  expect_equal(proposed$rho, current$rho)
})

test_that("prior_revision_parameters and log_jacobian_revision_parameters
return correct values", {
  revision_parameters <- list(eta=0.02, rho=0.3)
  prior <- list(eta_mu=0.01, eta_sigma=0.01, rho_alpha=2, rho_beta=5)
  expect_equal(prior_revision_parameters(revision_parameters, prior),
               dexp(0.02, 100, log=TRUE) + dbeta(0.3, 2, 5, log=TRUE))
  expect_equal(log_jacobian_revision_parameters(revision_parameters),
               log(0.02 * 0.3 * 0.7))
})

test_that("metropolis_step_revisions returns parameters and their log-probability", {
  set.seed(4)
  days_total <- 30
  r_list <- list(mean=3, sd=2)
  revision_parameters <- list(eta=0.2, rho=0.3)
  Rt_function <- stats::approxfun(1:days_total, rep(1.2, days_total))
  df <- generate_snapshots(days_total, Rt_function, list(mean=5, sd=3), r_list,
                           thinned=TRUE,
                           revision_parameters=revision_parameters) %>%
    dplyr::mutate(reporting_piece_index=1)
  r_params <- dplyr::tibble(reporting_piece_index=1, mean=3, sd=2)
  prior <- list(eta_mu=0.1, eta_sigma=0.1, rho_alpha=1, rho_beta=1)
  steps <- list(eta_step=0.2, rho_step=0.2)
  max_erroneous <- max_downward_revision(df) + 1
  f_logp <- function(parameters)
    observation_process_all_times_logp(df, r_params, parameters, max_erroneous) +
    prior_revision_parameters(parameters, prior) +
    log_jacobian_revision_parameters(parameters)

  # always accepts
  current <- list(eta=0.1, rho=0.5)
  res <- metropolis_step_revisions(df, r_params, current, -Inf, prior, steps,
                                   max_erroneous)
  expect_true(res$revision_parameters$eta != current$eta)
  expect_true(is.finite(res$logp))
  expect_equal(res$logp, f_logp(res$revision_parameters))

  # never accepts
  res <- metropolis_step_revisions(df, r_params, current, Inf, prior, steps,
                                   max_erroneous)
  expect_equal(res$revision_parameters, current)
  expect_equal(res$logp, Inf)

  # moves towards values with higher probability more often than not
  current <- list(eta=0.001, rho=0.5)
  logp_current <- f_logp(current)
  steps <- list(eta_step=0.2, rho_step=0)
  etas <- purrr::map_dbl(1:40, ~metropolis_step_revisions(
    df, r_params, current, logp_current, prior, steps,
    max_erroneous)$revision_parameters$eta)
  expect_true(sum(etas > current$eta) > sum(etas < current$eta))
})

test_that("mcmc works when downward revisions are allowed", {
  set.seed(5)
  niter <- 3
  days_total <- 40
  r_params <- list(mean=4, sd=2)
  s_params <- list(mean=5, sd=3)
  revision_parameters <- list(eta=0.1, rho=0.3)
  Rt_function <- stats::approxfun(1:days_total, rep(1.3, days_total))
  df <- generate_snapshots(days_total, Rt_function, s_params, r_params,
                           thinned=TRUE,
                           revision_parameters=revision_parameters)
  initial_cases_true <- df %>%
    dplyr::select(time_onset, cases_true) %>%
    unique()
  snapshot_with_Rt_index_df <- df %>%
    dplyr::select(-cases_true) %>%
    dplyr::mutate(Rt_index=ifelse(time_onset <= 20, 1, 2))
  initial_Rt <- dplyr::tibble(Rt_index=c(1, 2), Rt=c(1.3, 1.3))
  priors <- list(Rt=list(shape=1, rate=1),
                 reporting=list(mean_mu=5,
                                mean_sigma=10,
                                sd_mu=3,
                                sd_sigma=5),
                 revisions=list(eta_mu=0.1, eta_sigma=0.1,
                                rho_alpha=1, rho_beta=1),
                 max_erroneous=max_downward_revision(df) + 1,
                 max_cases=5000)
  f_mcmc <- function(priors, initial_revision_parameters, ...) {
    mcmc(niterations=niter,
         snapshot_with_Rt_index_df,
         priors,
         s_params,
         initial_cases_true,
         r_params,
         initial_Rt,
         print_to_screen=FALSE,
         initial_revision_parameters=initial_revision_parameters,
         ...)
  }

  # data with downward revisions need revision parameters
  expect_error(f_mcmc(priors, NULL))

  res <- f_mcmc(priors, revision_parameters)
  expect_equal(nrow(res$cases), days_total * niter)
  expect_equal(sum(is.na(res$cases$cases_true)), 0)
  expect_equal(nrow(res$revisions), niter)
  expect_true(all(c("eta", "rho", "iteration", "chain") %in% colnames(res$revisions)))
  expect_true(all(res$revisions$eta > 0))
  expect_true(all(res$revisions$rho > 0 & res$revisions$rho < 1))

  # fixing rho
  res <- f_mcmc(priors, revision_parameters,
                revision_metropolis_parameters=list(eta_step=0.5, rho_step=0))
  expect_equal(res$revisions$rho, rep(revision_parameters$rho, niter))

  # multiple chains
  res <- f_mcmc(priors, revision_parameters, nchains=2)
  expect_equal(nrow(res$revisions), niter * 2)
  expect_equal(sort(unique(res$revisions$chain)), c(1, 2))

  # errors
  priors_wrong <- priors
  priors_wrong$revisions <- NULL
  expect_error(f_mcmc(priors_wrong, revision_parameters))
  expect_error(f_mcmc(priors, revision_parameters, maximise=TRUE))
  expect_error(f_mcmc(priors, list(eta=0, rho=0.3)))
  expect_error(f_mcmc(priors, list(eta=0.1, rho=1)))
  priors_wrong <- priors
  priors_wrong$max_erroneous <- max_downward_revision(df) - 1
  expect_error(f_mcmc(priors_wrong, revision_parameters), "max_erroneous")
})

test_that("mcmc works when cases fall below those first reported long after
onset", {
  set.seed(6)
  niter <- 3
  days_total <- 40
  day_first_report <- 25
  r_params <- list(mean=4, sd=2)
  s_params <- list(mean=5, sd=3)
  revision_parameters <- list(eta=0.1, rho=0.05)
  Rt_function <- stats::approxfun(1:days_total, rep(1.3, days_total))
  # reports only available from some time after earliest onset days
  df <- generate_snapshots(days_total, Rt_function, s_params, r_params,
                           revision_parameters=revision_parameters) %>%
    dplyr::filter(time_reported >= day_first_report)
  df_first <- df %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(is_below_first=any(cases_reported < cases_reported[1]),
                     day_first=min(time_reported))
  expect_true(any(df_first$is_below_first &
                    df_first$day_first > df_first$time_onset))

  initial_cases_true <- df %>%
    dplyr::select(time_onset, cases_true) %>%
    unique()
  snapshot_with_Rt_index_df <- df %>%
    dplyr::select(-cases_true) %>%
    dplyr::mutate(Rt_index=ifelse(time_onset <= 20, 1, 2))
  priors <- list(Rt=list(shape=1, rate=1),
                 reporting=list(mean_mu=5,
                                mean_sigma=10,
                                sd_mu=3,
                                sd_sigma=5),
                 revisions=list(eta_mu=0.1, eta_sigma=0.1,
                                rho_alpha=1, rho_beta=1),
                 max_erroneous=max_downward_revision(df) + 1,
                 max_cases=5000)
  res <- mcmc(niterations=niter,
              snapshot_with_Rt_index_df,
              priors,
              s_params,
              initial_cases_true,
              r_params,
              dplyr::tibble(Rt_index=c(1, 2), Rt=c(1.3, 1.3)),
              print_to_screen=FALSE,
              initial_revision_parameters=revision_parameters)
  expect_equal(nrow(res$cases), days_total * niter)
  expect_equal(sum(is.na(res$cases$cases_true)), 0)
  expect_true(all(is.finite(res$reporting$mean)))
  expect_equal(nrow(res$revisions), niter)
})

test_that("sample_cases_history uses known cases from before the first onset time", {
  # the first onset time has a renewal mean only if earlier cases are known,
  # which then changes the conditional distribution of its case count
  s_params <- list(mean=3, sd=2)
  r_params <- dplyr::tibble(reporting_piece_index=1, mean=6, sd=3)
  df <- dplyr::tibble(time_onset=c(1, 1, 2, 2), time_reported=c(1, 2, 2, 3),
                      cases_reported=c(0, 2, 0, 1), reporting_piece_index=1)
  Rt_function <- function(t) rep(1, length(t))
  f_mean <- function(cases_before, seed) {
    set.seed(seed)
    mean(replicate(300, {
      suppressWarnings(sample_cases_history(
        df, 400, Rt_function, s_params, r_params,
        cases_before=cases_before)$cases_estimated[1])
    }))
  }
  low <- f_mean(rep(5, 10), 1)
  high <- f_mean(rep(60, 10), 2)
  expect_true(high > low + 10)
  # no history: no error and no renewal term for the first day
  expect_true(is.finite(f_mean(NULL, 3)))
})
