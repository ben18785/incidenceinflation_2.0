#' Samples a case count arising on a given onset day
#'
#' The distribution being drawn from is given by:
#' \deqn{p(cases_true_t|data, Rt, reporting_params, serial_params) \propto p(data|cases_true, reporting_params)
#'  p(cases_true_t|cases_true_t_1, cases_true_t_2, ..., Rt, serial_params)}
#'
#' If cases_future is supplied, the distribution additionally conditions on the
#' case counts arising on later onset days (see \code{conditional_cases_logp}).
#'
#' If revision_parameters is supplied, reported case counts can be revised
#' downwards, which means the true case count can be fewer than the maximum
#' reported (by up to max_erroneous).
#'
#' @param max_cases maximum possible cases thought to arise on a given day
#' @inheritParams conditional_cases_logp
#' @inheritParams conditional_cases_revisions_logp
#' @param ndraws number of draws of cases
#' @param maximise rather than sample a case count give the case count with the
#' maximum probability (by default is FALSE)
#'
#' @return a sampled case count arising on a given onset day
sample_true_cases_single_onset <- function(
  observation_df, cases_history, max_cases,
  Rt, day_onset, serial_parameters, reporting_parameters,
  ndraws=1, maximise=FALSE, kappa=NULL, is_negative_binomial=FALSE,
  cases_future=NULL, Rt_future=NULL,
  revision_parameters=NULL, max_erroneous=3,
  f_observation_logp=NULL, serial_max=40) {

  max_observed_cases <- max(observation_df$cases_reported)
  if(max_observed_cases > max_cases)
    stop("Max possible cases should be (much) greater than max observed cases.")
  if(is.null(revision_parameters) & is.null(f_observation_logp)) {
    possible_cases <- max_observed_cases:max_cases
    logps <- conditional_cases_logp(possible_cases, observation_df, cases_history,
                                    Rt, day_onset, serial_parameters, reporting_parameters,
                                    kappa, is_negative_binomial,
                                    cases_future, Rt_future,
                                    serial_max=serial_max)
  } else {
    # reported cases can exceed true cases by the erroneous records, including
    # those of any batch whose footprint contains the onset day
    max_excess <- ifelse(is.null(revision_parameters), 0, max_erroneous)
    if(!is.null(batch_for_onset(day_onset, revision_parameters$batches)))
      max_excess <- max_excess +
        ifelse(is.null(revision_parameters$max_batch), 30,
               revision_parameters$max_batch)
    min_cases <- max(0, max_observed_cases - max_excess)
    # the function supplied can say that fewer (or no fewer) are possible
    if(!is.null(attr(f_observation_logp, "min_cases")))
      min_cases <- attr(f_observation_logp, "min_cases")
    possible_cases <- min_cases:max_cases
    logps <- conditional_cases_revisions_logp(
      possible_cases, observation_df, cases_history,
      Rt, day_onset, serial_parameters, reporting_parameters,
      revision_parameters,
      kappa, is_negative_binomial,
      cases_future, Rt_future,
      max_erroneous,
      f_observation_logp=f_observation_logp,
      serial_max=serial_max)
  }
  probs <- exp(logps - matrixStats::logSumExp(logps))
  if(dplyr::last(probs) > 0.01)
    warning(paste0("Cases too few for onset day: ", day_onset,
                   ". Increase max_cases."))
  if(maximise)
    possible_cases[which.max(probs)]
  else
    sample(possible_cases, ndraws, prob=probs, replace=TRUE)
}

#' Calculates max number of days we are uncertain about reporting
#'
#' @param p_gamma_cutoff a p value (0 <= p <= 1) indicating the threshold above which
#' we deem certainty
#' @inheritParams conditional_cases_logp
#'
#' @return a number of days
max_uncertain_days <- function(p_gamma_cutoff, reporting_parameters) {
  r_mean <- reporting_parameters$mean
  r_sd <- reporting_parameters$sd
  days_from_end <- delay_quantile_mean_sd(p_gamma_cutoff, r_mean, r_sd)
  days_from_end
}

#' Draws a possible history (or histories) of cases
#'
#' The distribution being drawn from at each onset time t is given by:
#' \deqn{p(cases_true_t|data, Rt, reporting_params, serial_params) \propto p(data|cases_true, reporting_params)
#'  p(cases_true_t|cases_true_t_1, cases_true_t_2, ..., Rt, serial_params)}
#' If cases_true_current is supplied, the distribution at each onset time t
#' additionally conditions on the case counts arising on later onset days,
#' which depend on cases_true_t through the renewal equation:
#' \deqn{\prod_tau>t p(cases_true_tau|cases_true_tau-1, cases_true_tau-2, ..., Rt_tau, serial_params)}
#' Onset times are updated from the earliest to the latest, so that the case
#' counts for later onset times are those given in cases_true_current.
#'
#' @param observation_onset_df a tibble with four columns: time_onset, time_reported, cases_reported, reporting_piece_index
#' @param reporting_parameters a tibble with three columns: reporting_piece_index,
#' mean, sd.
#' @inheritParams observed_cases
#' @inheritParams sample_true_cases_single_onset
#' @inheritParams true_cases
#' @inheritParams max_uncertain_days
#' @param cases_true_current a tibble with two columns: "time_onset" and
#' "cases_true", which represents current estimates of the true number of cases
#' with each onset time. If NULL (the default) each case count is drawn
#' conditional only on case counts with earlier onset times.
#' @param cases_before a vector of known true case counts for the days
#' immediately before the first onset time, from the earliest to the latest
#' (optional). These enter the renewal equation as the history of cases for
#' the first onset times but are not themselves estimated, and nothing is
#' assumed about how they were reported.
#' @param future_max maximum number of later onset days conditioned on when
#' cases_true_current is supplied. Defaults to 40.
#' @param f_observation_logp_onset an optional function of an onset time
#' which returns a function giving the log-probability of the reports of that
#' onset time for a vector of possible true case counts. This replaces the
#' model of reporting otherwise used (see \code{deaths_cases_logp_function}).
#' If the function returned has an attribute 'min_cases', that is the smallest
#' true case count considered.
#' @param cases_offset an optional tibble with columns time_onset and
#' cases_offset: the amount by which the cases recorded for each onset time
#' exceed the true cases whose onset was on it (see
#' \code{deaths_cases_offset}). For onset times taken to be fully reported,
#' the true cases are the latest reported count less this amount.
#'
#' @return a tibble with an extra cases_estimated column
#' @export
#' @importFrom rlang .data
sample_cases_history <- function(
  observation_onset_df, max_cases,
  Rt_function, serial_parameters, reporting_parameters,
  p_gamma_cutoff=0.99,
  maximise=FALSE,
  kappa=NULL,
  is_negative_binomial=FALSE,
  cases_true_current=NULL,
  future_max=40,
  revision_parameters=NULL,
  max_erroneous=3,
  cases_before=NULL,
  f_observation_logp_onset=NULL,
  cases_offset=NULL) {

  if(!"reporting_piece_index" %in% colnames(observation_onset_df))
    stop("observation_onset_df must contain a column: 'reporting_piece_index'.")

  # TODO make this uncertain period determination better; put into separate function
  latest_onset_time <- max(observation_onset_df$time_onset)
  df_latest_onset_time <- observation_onset_df %>%
    dplyr::filter(time_onset==latest_onset_time)
  reporting_index_latest_onset_time <- df_latest_onset_time$reporting_piece_index[1]
  reporting_latest_onset_time <- reporting_parameters %>%
    dplyr::filter(reporting_piece_index == reporting_index_latest_onset_time)
  reporting_latest_onset_time <- list(mean=reporting_latest_onset_time$mean[1],
                                      sd=reporting_latest_onset_time$sd[1])
  uncertain_period <- max_uncertain_days(p_gamma_cutoff, reporting_latest_onset_time)
  # a heavier-tailed delay extends the period for which reporting is incomplete
  if(!is.null(revision_parameters$tail))
    uncertain_period <- delay_quantile(p_gamma_cutoff, reporting_latest_onset_time,
                                       revision_parameters)

  start_uncertain_period <- latest_onset_time - uncertain_period

  observation_history_df <- observation_onset_df %>%
    dplyr::group_by(.data$time_onset) %>%
    dplyr::mutate(cases_estimated=ifelse(
      .data$time_onset < start_uncertain_period,
      # latest count, since counts can be revised downwards
      .data$cases_reported[which.max(.data$time_reported)], NA)) %>%
    dplyr::ungroup()
  # the cases whose onset was on a day can differ from those recorded for it
  if(!is.null(cases_offset))
    observation_history_df$cases_estimated <-
      observation_history_df$cases_estimated -
      cases_offset$cases_offset[match(observation_history_df$time_onset,
                                      cases_offset$time_onset)]
  onset_times <- unique(observation_history_df$time_onset)
  onset_times_uncertain_period <- onset_times[onset_times >= start_uncertain_period]

  condition_on_future <- !is.null(cases_true_current)
  if(condition_on_future) {
    if(!all(c("time_onset", "cases_true") %in% colnames(cases_true_current)))
      stop("cases_true_current must contain columns: 'time_onset' and 'cases_true'.")
    if(!all(onset_times_uncertain_period %in% cases_true_current$time_onset))
      stop("cases_true_current must contain a case count for each onset time.")
    index_current <- match(observation_history_df$time_onset,
                           cases_true_current$time_onset)
    observation_history_df$cases_estimated <- ifelse(
      is.na(observation_history_df$cases_estimated),
      cases_true_current$cases_true[index_current],
      observation_history_df$cases_estimated)
  }

  is_misdating <- !is.null(revision_parameters) &&
    misdating_prob(revision_parameters) > 0
  if(is_misdating) {
    if(!condition_on_future)
      stop("cases_true_current must be supplied if cases can be misdated.")
    observation_sorted_df <- observation_history_df %>%
      dplyr::arrange(.data$time_onset, .data$time_reported)
    l_observations <- split(
      as.data.frame(observation_sorted_df[, c("time_reported", "cases_reported")]),
      observation_sorted_df$time_onset)
    reporting_by_onset_df <- observation_sorted_df %>%
      dplyr::select("time_onset", "reporting_piece_index") %>%
      unique() %>%
      dplyr::left_join(reporting_parameters, by="reporting_piece_index")
    l_reporting <- stats::setNames(
      purrr::map(seq_len(nrow(reporting_by_onset_df)),
                 ~list(mean=reporting_by_onset_df$mean[.],
                       sd=reporting_by_onset_df$sd[.])),
      reporting_by_onset_df$time_onset)
    cases_estimated_df <- observation_sorted_df %>%
      dplyr::select("time_onset", "cases_estimated") %>%
      unique()
    cases_by_onset <- stats::setNames(cases_estimated_df$cases_estimated,
                                      cases_estimated_df$time_onset)
  }

  for(i in seq_along(onset_times_uncertain_period)) {

    onset_time <- onset_times_uncertain_period[i]
    snapshots_at_onset_time_df <- observation_history_df %>%
      dplyr::filter(.data$time_onset==onset_time)
    reporting_index <- snapshots_at_onset_time_df$reporting_piece_index[1]
    reporting_parameters_tmp <- reporting_parameters %>%
      dplyr::filter(reporting_piece_index == reporting_index)
    reporting_parameters_tmp <- list(
      mean=reporting_parameters_tmp$mean,
      sd=reporting_parameters_tmp$sd
      )
    snapshots_at_onset_time_df <- snapshots_at_onset_time_df %>%
      dplyr::select("time_reported", "cases_reported")
    pre_observation_df <- observation_history_df %>%
      dplyr::filter(.data$time_onset < onset_time) %>%
      dplyr::select("time_onset", "cases_estimated") %>%
      unique() %>%
      dplyr::arrange(dplyr::desc(.data$time_onset))
    # followed by any known cases from before the first onset time
    cases_history <- c(pre_observation_df$cases_estimated, rev(cases_before))
    Rt <- Rt_function(onset_time)

    cases_future <- NULL
    Rt_future <- NULL
    if(condition_on_future) {
      post_observation_df <- observation_history_df %>%
        dplyr::filter(.data$time_onset > onset_time,
                      .data$time_onset <= onset_time + future_max) %>%
        dplyr::select("time_onset", "cases_estimated") %>%
        unique() %>%
        dplyr::arrange(.data$time_onset)
      cases_future <- post_observation_df$cases_estimated
      Rt_future <- Rt_function(post_observation_df$time_onset)
    }

    f_observation_logp <- NULL
    if(is_misdating)
      f_observation_logp <- misdating_observation_logp_function(
        onset_time, l_observations, cases_by_onset, l_reporting,
        revision_parameters, max_erroneous)
    if(!is.null(f_observation_logp_onset))
      f_observation_logp <- f_observation_logp_onset(onset_time)

    case <- sample_true_cases_single_onset(
      observation_df=snapshots_at_onset_time_df,
      cases_history=cases_history,
      max_cases=max_cases,
      Rt=Rt,
      day_onset=onset_time,
      serial_parameters=serial_parameters,
      reporting_parameters=reporting_parameters_tmp,
      ndraws=1,
      maximise=maximise,
      kappa=kappa,
      is_negative_binomial=is_negative_binomial,
      cases_future=cases_future,
      Rt_future=Rt_future,
      revision_parameters=revision_parameters,
      max_erroneous=max_erroneous,
      f_observation_logp=f_observation_logp,
      serial_max=future_max)
    index_onset_time <- which(observation_history_df$time_onset==onset_time)
    observation_history_df$cases_estimated[index_onset_time] <- case
    if(is_misdating)
      cases_by_onset[as.character(onset_time)] <- case
  }
  observation_history_df
}

#' Constructs a matrix of w vectors
#'
#' @param w a vector of weights
#' @param piece_width the width of an Rt piece
#'
#' @return a matrix of ws
construct_w_matrix <- function(w, piece_width) {
  wmax <- length(w)
  m_w <- matrix(nrow = piece_width,
                ncol = (wmax + piece_width - 1))
  for(i in 1:piece_width) {
    n_before_zeros <- i - 1
    n_trailing_zeros <- piece_width - n_before_zeros - 1
    m_w[i, ] <- c(rep(0, i - 1), w, rep(0, n_trailing_zeros))
  }
  m_w
}

#' Calculates negative-binomial log-likelhood within a single Rt piece
#'
#' @param Rt an Rt value for the piece
#' @param kappa overdispersion parameter
#' @param w weights corresponding to the generation times
#' @param onset_times the onset times corresponding to the piece
#' @param cases_df a tibble with 'cases_true' as a column which has been ordered
#' so that the latest onset times are at the bottom
#'
#' @return a log-likelihood of the piece
nb_log_likelihood_Rt_piece <- function(Rt, kappa, w, onset_times, cases_df) {

  matches <- match(onset_times, cases_df$time_onset)
  if(max(matches) != nrow(cases_df))
    stop("Onset times must be last entries in cases_df.")

  piece_width <- length(onset_times)
  m_w <- construct_w_matrix(w, piece_width)
  wmax <- length(w)
  v_i <- rev(cases_df$cases_true)[1:(wmax + piece_width)]
  v_i <- ifelse(!is.na(v_i), v_i, 0) # necessary in case we reach back before start of data
  v_i_dep <- v_i[1:piece_width]
  v_i_ind <- v_i[2:(wmax + piece_width)]
  log_prob <- stats::dnbinom(v_i_dep, mu=(Rt * m_w %*% v_i_ind), size=kappa,
          log=TRUE)
  if(onset_times[1] == min(cases_df$time_onset)) { # first case can't be generated from nothing
    log_prob <- log_prob[1:(length(log_prob) - 1)]
  }

  sum(log_prob)
}

#' Samples Rt for a single piece under a negative binomial renewal model using
#' Metropolis-Hastings steps
#'
#' Under a negative binomial renewal model the conditional distribution of Rt
#' has no closed form. Each step proposes a new value in one of two ways,
#' chosen at random: independently of the current value, from a gamma
#' distribution centred on the conditional distribution under a Poisson model
#' and widened to allow for overdispersion; or by a random walk on the log of
#' the current value. The proposal is accepted or rejected in the usual way, so
#' that the draws target the conditional distribution exactly.
#'
#' @param prior_shape Rt prior shape parameter
#' @param prior_rate Rt prior rate parameter
#' @param posterior_shape shape parameter of the conditional distribution of Rt
#' under a Poisson renewal model
#' @param posterior_rate rate parameter of the conditional distribution of Rt
#' under a Poisson renewal model
#' @inheritParams nb_log_likelihood_Rt_piece
#' @param ndraws number of draws of Rt to return
#' @param nsteps number of Metropolis-Hastings steps taken for each draw
#' @param Rt_current the current value of Rt, from which the steps start; if
#' NULL they start from the mean under a Poisson model
#'
#' @return a vector of Rt draws
sample_nb_Rt_piece <- function(prior_shape, prior_rate,
          posterior_shape, posterior_rate,
          kappa,
          w,
          onset_times,
          cases_df,
          ndraws,
          nsteps=10,
          Rt_current=NULL) {

  # proposal: the Poisson conditional with its sd inflated in line with the
  # variance of the negative binomial relative to the Poisson
  mu <- posterior_shape / posterior_rate
  sd <- sqrt(posterior_shape) / posterior_rate
  mean_daily_cases <- (posterior_shape - prior_shape) / length(onset_times)
  sd <- 1.5 * sd * sqrt(1 + mean_daily_cases / kappa)
  shape_q <- max(mu^2 / sd^2, 1)
  rate_q <- shape_q / mu
  step_log <- 1 / sqrt(shape_q)

  f_target <- function(Rt) {
    if(!is.finite(Rt) || Rt <= 0)
      return(-Inf)
    nb_log_likelihood_Rt_piece(Rt, kappa, w, onset_times, cases_df) +
      stats::dgamma(Rt, prior_shape, prior_rate, log=TRUE)
  }

  Rt <- Rt_current
  if(is.null(Rt) || !is.finite(Rt) || Rt <= 0)
    Rt <- mu
  logp <- f_target(Rt)
  draws <- vector(length=ndraws)
  for(i in seq_len(ndraws)) {
    for(j in seq_len(nsteps)) {
      if(stats::runif(1) < 0.5) {
        Rt_proposed <- stats::rgamma(1, shape_q, rate_q)
        log_q <- stats::dgamma(Rt, shape_q, rate_q, log=TRUE) -
          stats::dgamma(Rt_proposed, shape_q, rate_q, log=TRUE)
      } else {
        Rt_proposed <- Rt * exp(stats::rnorm(1, 0, step_log))
        # Jacobian of proposing on the log scale
        log_q <- log(Rt_proposed) - log(Rt)
      }
      logp_proposed <- f_target(Rt_proposed)
      log_r <- logp_proposed - logp + log_q
      # a current value of zero probability is always left
      if(!is.finite(logp) && is.finite(logp_proposed))
        log_r <- Inf
      if(!is.na(log_r) && log(stats::runif(1)) < log_r) {
        Rt <- Rt_proposed
        logp <- logp_proposed
      }
    }
    draws[i] <- Rt
  }
  draws
}

#' Draws from the gamma distribution or returns the value which maximises
#' it
#'
#' @param shape the shape parameter of a gamma distribution
#' @param rate the rate parameter of a gamma distribution
#' @param ndraws number of draws if maximise=FALSE
#' @param maximise whether to return the mode of the gamma distribution
#'
#' @return a value or (if ndraws > 1) a vector of values
sample_or_maximise_gamma <- function(shape, rate, ndraws, maximise=FALSE) {
  if(maximise)
    (shape - 1) / rate
  else
    stats::rgamma(ndraws, shape, rate)
}


#' Sample a single Rt value corresponding to a single piecewise-
#' constant element of an Rt vector
#'
#' @param Rt_piece_index the index of the Rt piece being sampled
#' @param cases_history_df a tibble with three columns: time_onset, cases_true
#' and Rt_index
#' @param Rt_prior_parameters a list with elements 'shape' and 'rate' describing
#' the gamma prior for Rt
#' @inheritParams sample_cases_history
#' @param ndraws number of draws of Rt
#' @inheritParams true_cases
#' @param nsteps number of Metropolis-Hastings steps taken for each draw of Rt
#' if assuming a negative binomial model (see \code{sample_nb_Rt_piece})
#' @param Rt_current the current value of Rt for the piece, from which these
#' steps start (optional)
#'
#' @return a draw (or draws) for Rt
#' @importFrom rlang .data
sample_Rt_single_piece <- function(
    Rt_piece_index, cases_history_df,
    Rt_prior_parameters, serial_parameters,
    kappa=NULL,
    serial_max=40, ndraws=1,
    maximise=FALSE, is_negative_binomial=FALSE,
    nsteps=10, Rt_current=NULL, cases_before=NULL) {

  short_df <- cases_history_df %>%
    dplyr::filter(.data$Rt_index <= Rt_piece_index)
  # known cases from before the first onset time form part of the history
  if(length(cases_before) > 0) {
    time_first <- min(cases_history_df$time_onset)
    short_df <- dplyr::bind_rows(
      dplyr::tibble(
        time_onset=time_first - rev(seq_along(cases_before)),
        cases_true=cases_before,
        Rt_index=NA),
      short_df)
  }
  time_max_post_initial_period <- max(short_df$time_onset) - serial_max

  posterior_shape <- Rt_prior_parameters$shape
  posterior_rate <- Rt_prior_parameters$rate

  short_df <- short_df %>%
    dplyr::mutate(time_after_start = .data$time_onset - serial_max)
  onset_times <- short_df %>%
    dplyr::filter(.data$Rt_index == Rt_piece_index) %>%
    dplyr::pull(.data$time_onset)

  w <- weights_series(serial_max, serial_parameters)
  for(i in seq_along(onset_times)) {

    onset_time <- onset_times[i]
    true_cases <- short_df %>%
      dplyr::filter(.data$time_onset == onset_time) %>%
      dplyr::pull(.data$cases_true)
    posterior_shape <- posterior_shape + true_cases
    cases_history <- short_df %>%
      dplyr::filter(.data$time_onset < onset_time) %>%
      dplyr::arrange(dplyr::desc(.data$time_onset)) %>%
      dplyr::pull(.data$cases_true)
    diff_time <- serial_max - length(cases_history)
    if(diff_time <= 0)
      cases_history <- cases_history[1:serial_max]
    else
      cases_history <- c(cases_history, rep(0, diff_time))
    posterior_rate <- posterior_rate + sum(w * cases_history)
  }

  if(!is_negative_binomial) {
    sample_or_maximise_gamma(
      posterior_shape, posterior_rate, ndraws, maximise)
  } else {
    sample_nb_Rt_piece(Rt_prior_parameters$shape, Rt_prior_parameters$rate,
              posterior_shape, posterior_rate,
              kappa,
              w,
              onset_times,
              short_df,
              ndraws,
              nsteps,
              Rt_current)
  }
}


#' Sample piecewise-constant Rt values
#'
#' If the renewal model is specified by a Poisson (the default):
#' \deqn{cases_true_t ~ Poisson(Rt * \sum_tau=1^t_max w_t cases_true_t-tau)}
#' If an Rt value is given a gamma prior, this results in a posterior
#' distribution:
#' \deqn{Rt ~ gamma(alpha + cases_true_t, beta + \sum_tau=1^t_max w_t cases_true_t-tau))}
#' where alpha and beta are the shape and rate parameters of the gamma
#' prior distribution. Here, we assume that Rt is constant over a set of onset
#' times 'onset_time_set'. This means that the posterior for a single Rt value is given
#' by:
#' \deqn{Rt ~ gamma(alpha + \sum_{t in onset_time_set} cases_true_t,
#'       beta + \sum_{t in onset_time_set}\sum_tau=1^t_max w_t cases_true_t-tau))}
#' This function either returns a draw (or draws if ndraws>1) from
#' this posterior, or it returns the Rt set that maximises it
#' (if maximise=TRUE).
#' Alternatively the renewal equation may be specified by a negative binomial distribution:
#' \deqn{cases_true_t ~ NB(Rt * \sum_tau=1^t_max w_t cases_true_t-tau, kappa)}
#' where kappa is the overdispersion parameter. In this case Rt is updated by
#' Metropolis-Hastings steps (see \code{sample_nb_Rt_piece}).
#'
#' @inheritParams sample_Rt_single_piece
#' @inheritParams sample_cases_history
#' @param Rt_current a tibble with columns 'Rt_index' and 'Rt' giving the
#' current value of Rt for each piece, from which Metropolis-Hastings steps
#' start under a negative binomial model (optional)
#' @return a tibble with three columns: "Rt_piece_index", "draw_index", "Rt"
#' @export
sample_Rt <- function(cases_history_df,
                      Rt_prior_parameters,
                      serial_parameters,
                      kappa=NULL,
                      serial_max=40,
                      ndraws=1,
                      maximise=FALSE,
                      is_negative_binomial=FALSE,
                      nsteps=10,
                      Rt_current=NULL,
                      cases_before=NULL) {

  if(is_negative_binomial) {
    if(is.null(kappa)) {
      stop("Overdispersion parameter must not be null if using a negative binomial model.")
    } else if(kappa <= 0){
      stop("Overdispersion parameter must be positive.")
    }
  }
  Rt_piece_indices <- unique(cases_history_df$Rt_index)
  num_Rt_pieces <- length(Rt_piece_indices)
  if(maximise)
    ndraws <- 1
  draw_indices <- seq(1, ndraws, 1)
  m_draws <- matrix(nrow = num_Rt_pieces * ndraws,
                    ncol = 3)
  k <- 1
  for(i in seq_along(Rt_piece_indices)) {
    Rt_piece_index <- Rt_piece_indices[i]
    Rt_piece_current <- NULL
    if(!is.null(Rt_current))
      Rt_piece_current <- Rt_current$Rt[match(Rt_piece_index, Rt_current$Rt_index)]
    Rt_vals <- sample_Rt_single_piece(
      Rt_piece_index, cases_history_df,
      Rt_prior_parameters, serial_parameters,
      kappa,
      serial_max, ndraws, maximise=maximise,
      is_negative_binomial=is_negative_binomial,
      nsteps=nsteps,
      Rt_current=Rt_piece_current,
      cases_before=cases_before)
    for(j in 1:ndraws) {
      m_draws[k, ] <- c(Rt_piece_index, j, Rt_vals[j])
      k <- k + 1
    }
  }
  colnames(m_draws) <- c("Rt_index", "draw_index", "Rt")
  m_draws <- m_draws %>%
    dplyr::as_tibble()
  m_draws
}

#' Propose new reporting parameters using normal kernel
#' centered at current values
#'
#' @param current_reporting_parameters a tibble with column names: "reporting_piece_index", "mean", "sd"
#' @param metropolis_parameters named list of 'mean_step', 'sd_step' containing
#' step sizes for Metropolis step
#'
#' @return a tibble with column names: "reporting_piece_index", "mean", "sd"
propose_reporting_parameters <- function(
  current_reporting_parameters,
  metropolis_parameters) {

  mean_now <- current_reporting_parameters$mean
  sd_now <- current_reporting_parameters$sd
  mean_stepsize <- metropolis_parameters$mean_step
  sd_stepsize <- metropolis_parameters$sd_step
  mean_proposed <- purrr::map_dbl(mean_now, ~stats::rnorm(1, ., mean_stepsize))
  sd_proposed <- purrr::map_dbl(sd_now, ~stats::rnorm(1, ., sd_stepsize))

  current_reporting_parameters %>%
    dplyr::mutate(
      mean=mean_proposed,
      sd=sd_proposed
    )
}

#' Gamma prior for reporting parameters
#'
#' @inheritParams propose_reporting_parameters
#' @param prior_parameters named list with elements 'mean_mu', 'mean_sigma', 'sd_mu',
#' 'sd_sigma' representing the gamma prior parameters for the mean and sd
#' parameters of the reporting parameters (itself described by a gamma
#' distribution)
#'
#' @return a log-probability density
prior_reporting_parameters <- function(
  current_reporting_parameters,
  prior_parameters) {

  num_reporting_parameters <- max(current_reporting_parameters$reporting_piece_index)

  # assumes that all reporting parameters have same prior
  mean <- current_reporting_parameters$mean
  sd <- current_reporting_parameters$sd
  logp_mean <- sum(dgamma_mean_sd(mean,
                                  prior_parameters$mean_mu,
                                  prior_parameters$mean_sigma,
                                  log=TRUE))
  logp_sigma <- sum(dgamma_mean_sd(sd,
                                   prior_parameters$sd_mu,
                                   prior_parameters$sd_sigma,
                                   log=TRUE))
  logp_mean + logp_sigma
}

#' Performs Metropolis accept-reject step
#'
#' @param logp_current current log-probability value
#' @param logp_proposed proposed log-probability value
#' @param current_parameters current parameter(s) value
#' @param proposed_parameters proposed parameter(s) value
#'
#' @return a named list with two elements: 'parameter', a parameter value (which
#' may be a non-scalar); and 'logp', the log-probability corresponding to the
#' parameter value returned
accept_reject <- function(
    logp_current, logp_proposed,
    current_parameters,
    proposed_parameters) {

  log_r <- logp_proposed - logp_current
  log_u <- log(stats::runif(1))

  # a proposal with zero probability is rejected, including where the current
  # value also has zero probability
  if(is.na(log_r))
    log_r <- -Inf

  # nocov start
  if(log_r > log_u) {
    new_parameters <- proposed_parameters
    logp <- logp_proposed
  } else {
    new_parameters <- current_parameters
    logp <- logp_current
  }
  # nocov end
  list(
    parameter=new_parameters,
    logp=logp
  )
}


#' Sample reporting parameters using a single Metropolis step
#'
#' @inheritParams observation_process_all_times_logp
#' @inheritParams propose_reporting_parameters
#' @inheritParams prior_reporting_parameters
#' @param logp_current the value fo the log-probability at current parameter values
#'
#' @return a named list with two elements: 'parameter', a parameter value (which
#' may be a non-scalar); and 'logp', the log-probability corresponding to the
#' parameter value returned
metropolis_step <- function(snapshot_with_true_cases_df,
                            current_reporting_parameters,
                            logp_current,
                            prior_parameters,
                            metropolis_parameters,
                            revision_parameters=NULL,
                            max_erroneous=3) {

  proposed_reporting_parameters <- propose_reporting_parameters(
    current_reporting_parameters,
    metropolis_parameters)
  logp_proposed <- observation_process_all_times_logp(
    snapshot_with_true_cases_df=snapshot_with_true_cases_df,
    reporting_parameters=proposed_reporting_parameters,
    revision_parameters=revision_parameters,
    max_erroneous=max_erroneous
  ) + prior_reporting_parameters(proposed_reporting_parameters,
                                 prior_parameters)

  list_parameters_logp <- accept_reject(
    logp_current, logp_proposed,
    current_reporting_parameters,
    proposed_reporting_parameters)

  list_parameters_logp
}

#' Sample reporting parameters using Metropolis MCMC
#'
#' @inheritParams metropolis_step
#' @param ndraws number of iterates of the Markov chain to simulate
#'
#' @return a named list with two elements: 'reporting_parameters', a tibble with
#' four columns: "reporting_piece_index", "draw_index", "mean, "sd"; and
#' 'logp' a value of the log-probability for the current parameter values (which
#' is NULL when maximising since this is not used)
metropolis_steps <- function(
  snapshot_with_true_cases_df,
  current_reporting_parameters,
  logp_current,
  prior_parameters,
  metropolis_parameters,
  ndraws,
  revision_parameters=NULL,
  max_erroneous=3) {

  reporting_parameters <- current_reporting_parameters
  num_reporting_parameters <- max(reporting_parameters$reporting_piece_index)
  m_reporting <- matrix(ncol = 4, nrow = ndraws * num_reporting_parameters)
  k <- 1
  for(i in 1:ndraws) {
    list_parameters_logp <- metropolis_step(
      snapshot_with_true_cases_df,
      reporting_parameters,
      logp_current,
      prior_parameters,
      metropolis_parameters,
      revision_parameters,
      max_erroneous
    )
    reporting_parameters <- list_parameters_logp$parameter
    logp_current <- list_parameters_logp$logp
    for(j in 1:num_reporting_parameters) {
      m_reporting[k, ] <- c(j,
                            i,
                            reporting_parameters$mean[j],
                            reporting_parameters$sd[j])
      k <- k + 1
    }
  }
  colnames(m_reporting) <- c("reporting_piece_index", "draw_index", "mean", "sd")
  m_reporting <- m_reporting %>%
    dplyr::as_tibble()

  list(
    reporting_parameters=m_reporting,
    logp=logp_current
  )
}

#' Select reporting parameters by maximising log-probability
#'
#' @inheritParams metropolis_step
#'
#' @return a tibble with four columns: "reporting_piece_index", "draw_index", "mean, "sd"
maximise_reporting_logp <- function(
  snapshot_with_true_cases_df,
  current_reporting_parameters,
  prior_parameters) {

  if(!"reporting_piece_index" %in% colnames(snapshot_with_true_cases_df))
    stop("snapshot_with_true_cases_df must contain a column: 'reporting_piece_index'.")

  if(!"reporting_piece_index" %in% colnames(current_reporting_parameters))
    stop("current_reporting_parameters must contain a column: 'reporting_piece_index'.")

  num_reporting_indices <- dplyr::n_distinct(snapshot_with_true_cases_df$reporting_piece_index)

  if(num_reporting_indices == 1) {

    objective_function <- function(theta) {
      -observation_process_all_times_logp(
        snapshot_with_true_cases_df,
        current_reporting_parameters %>%
          dplyr::mutate(mean=theta[1],
                        sd=theta[2])) +
        prior_reporting_parameters(
          list(mean=theta[1],
               sd=theta[2],
               reporting_piece_index=current_reporting_parameters$reporting_piece_index),
          prior_parameters)
    }

    start_point <- c(current_reporting_parameters$mean,
                     current_reporting_parameters$sd)
    theta <- stats::optim(start_point, objective_function)$par
    overall_reporting_parameters <- current_reporting_parameters %>%
      dplyr::mutate(mean=theta[1],
                    sd=theta[2],
                    draw_index=1)

  } else {

    # since reporting delays determined only by corresponding onset times
    # we can optimise each piece separately
    for(i in 1:num_reporting_indices) {

      df_temp <- snapshot_with_true_cases_df %>%
        dplyr::filter(reporting_piece_index==i)
      reporting_piece <- current_reporting_parameters %>%
        dplyr::filter(reporting_piece_index==i)
      present_reporting <- maximise_reporting_logp(df_temp, reporting_piece, prior_parameters)

      if(i == 1)
        overall_reporting_parameters <- present_reporting
      else
        overall_reporting_parameters <- overall_reporting_parameters %>%
        dplyr::bind_rows(present_reporting)
    }

  }

  overall_reporting_parameters
}

#' Propose new overdispersion parameter by sampling from a normal centered
#' on the current value
#'
#' @param overdispersion_current an overdispersion parameter value (should exceed 0)
#' @param overdispersion_metropolis_sd the standard deviation of the proposal kernel
#'
#' @return a proposed overdispersion parameter value
propose_overdispersion_parameter <- function(
    overdispersion_current,
    overdispersion_metropolis_sd) {

  stats::rnorm(1, overdispersion_current,
               overdispersion_metropolis_sd)
}

#' Performs a single Metropolis step to update overdispersion parameter
#'
#' @param overdispersion_current current overdispersion parameter value (must exceed 0)
#' @param logp_current current log-probability value from previous update
#' @inheritParams state_process_nb_logp_all_onsets
#' @param prior_overdispersion_parameter a named list with elements: 'mean' and
#' 'sd' denoting the mean and sd of a gamma distribution
#' @inheritParams propose_overdispersion_parameter
#'
#' @return a named list with two elements: 'overdispersion', an overdispersion parameter
#' value and 'logp', the log-probability corresponding to the returned parameter
#' set
metropolis_step_overdispersion <- function(
    overdispersion_current,
    logp_current,
    cases_history_rt_df,
    serial_parameters,
    prior_overdispersion_parameter,
    overdispersion_metropolis_sd,
    serial_max=40,
    cases_before=NULL
    ) {

  overdispersion_proposed <- propose_overdispersion_parameter(
    overdispersion_current,
    overdispersion_metropolis_sd
  )

  logp_prior <- dgamma_mean_sd(overdispersion_proposed,
                               prior_overdispersion_parameter$mean,
                               prior_overdispersion_parameter$sd,
                               log=TRUE)
  logp_proposed <- state_process_nb_logp_all_onsets(
    overdispersion_proposed, cases_history_rt_df, serial_parameters,
    serial_max=serial_max, cases_before=cases_before
  ) + logp_prior

  list_parameters_logp <- accept_reject(
    logp_current, logp_proposed,
    overdispersion_current, overdispersion_proposed)

  list(
    overdispersion=list_parameters_logp$parameter,
    logp=list_parameters_logp$logp
  )
}

#' Draw reporting parameter values either by sampling or by
#' maximising
#'
#' @inheritParams metropolis_steps
#' @param maximise if true choose reporting parameters by maximising
#' log-probability; else (default) use Metropolis MCMC
#' to draw parameters
#' @param logp_current the value of the log-probability at the current parameter values
#'
#' @return a named list with two elements: 'reporting_parameters', a tibble with
#' four columns: "reporting_piece_index", "draw_index", "mean, "sd"; and
#' 'logp' a value of the log-probability for the current parameter values (which
#' is NULL when maximising since this is not used)
#' @export
sample_reporting <- function(
  snapshot_with_true_cases_df,
  current_reporting_parameters,
  logp_current,
  prior_parameters,
  metropolis_parameters,
  maximise=FALSE,
  ndraws=1,
  revision_parameters=NULL,
  max_erroneous=3) {

  if(maximise) {
    reporting_parameters <- maximise_reporting_logp(
      snapshot_with_true_cases_df,
      current_reporting_parameters,
      prior_parameters)
    list_parameter_logp <- list(
      reporting_parameters=reporting_parameters,
      logp=NULL)
  } else {
    list_parameter_logp <- metropolis_steps(
      snapshot_with_true_cases_df,
      current_reporting_parameters,
      logp_current,
      prior_parameters,
      metropolis_parameters,
      ndraws=ndraws,
      revision_parameters=revision_parameters,
      max_erroneous=max_erroneous)
  }

  list_parameter_logp
}

#' Propose new revision parameters
#'
#' Proposals are made using normal kernels centered at the current values of
#' log(eta) and logit(rho), so that proposed values always lie within bounds.
#'
#' @param current_revision_parameters named list of 'eta' and 'rho' and, if
#' cases can be misdated, 'pi', 'shifts' and 'shift_probs' (see
#' \code{unrecorded_prob})
#' @param metropolis_parameters named list of 'eta_step' and 'rho_step' (and
#' optionally 'pi_step') containing step sizes (on the log, logit and logit
#' scales respectively). A step size of zero holds a parameter fixed at its
#' current value.
#'
#' @return a named list like current_revision_parameters with proposed values
propose_revision_parameters <- function(
    current_revision_parameters,
    metropolis_parameters) {

  log_eta <- stats::rnorm(1, log(current_revision_parameters$eta),
                          metropolis_parameters$eta_step)
  logit_rho <- stats::rnorm(1, stats::qlogis(current_revision_parameters$rho),
                            metropolis_parameters$rho_step)
  proposed_revision_parameters <- current_revision_parameters
  proposed_revision_parameters$eta <- exp(log_eta)
  proposed_revision_parameters$rho <- stats::plogis(logit_rho)

  # misdating probability (if cases can be misdated)
  if(misdating_prob(current_revision_parameters) > 0 &
     !is.null(metropolis_parameters$pi_step)) {
    logit_pi <- stats::rnorm(1, stats::qlogis(current_revision_parameters$pi),
                             metropolis_parameters$pi_step)
    proposed_revision_parameters$pi <- stats::plogis(logit_pi)
  }
  proposed_revision_parameters
}

#' Prior for revision parameters
#'
#' A gamma prior is assumed for eta and a beta prior for rho (and for pi, if
#' cases can be misdated).
#'
#' @inheritParams propose_revision_parameters
#' @param prior_parameters named list with elements 'eta_mu', 'eta_sigma',
#' representing the mean and sd of the gamma prior for eta; 'rho_alpha',
#' 'rho_beta' representing the shape parameters of the beta prior for rho; and,
#' if cases can be misdated, 'pi_alpha', 'pi_beta' representing the shape
#' parameters of the beta prior for pi
#'
#' @return a log-probability density
prior_revision_parameters <- function(
    current_revision_parameters,
    prior_parameters) {

  logp <- dgamma_mean_sd(current_revision_parameters$eta,
                         prior_parameters$eta_mu,
                         prior_parameters$eta_sigma,
                         log=TRUE) +
    stats::dbeta(current_revision_parameters$rho,
                 prior_parameters$rho_alpha,
                 prior_parameters$rho_beta,
                 log=TRUE)
  if(misdating_prob(current_revision_parameters) > 0)
    logp <- logp + stats::dbeta(current_revision_parameters$pi,
                                prior_parameters$pi_alpha,
                                prior_parameters$pi_beta,
                                log=TRUE)
  logp
}

#' Log-determinant of Jacobian of transformation from (log(eta), logit(rho))
#' to (eta, rho)
#'
#' Revision parameters are proposed on the transformed scale, so this is added
#' to the log-probability density of (eta, rho) in a Metropolis step.
#'
#' @inheritParams propose_revision_parameters
#'
#' @return a log-determinant
log_jacobian_revision_parameters <- function(current_revision_parameters) {
  eta <- current_revision_parameters$eta
  rho <- current_revision_parameters$rho
  log_jacobian <- log(eta) + log(rho) + log(1 - rho)
  p_misdated <- misdating_prob(current_revision_parameters)
  if(p_misdated > 0)
    log_jacobian <- log_jacobian + log(p_misdated) + log(1 - p_misdated)
  log_jacobian
}

#' Sample revision parameters using a single Metropolis step
#'
#' @inheritParams observation_process_all_times_logp
#' @inheritParams propose_revision_parameters
#' @inheritParams prior_revision_parameters
#' @param logp_current the value of the log-probability at current parameter
#' values: this is the sum of \code{observation_process_all_times_logp},
#' \code{prior_revision_parameters} and \code{log_jacobian_revision_parameters}
#'
#' @return a named list with two elements: 'revision_parameters', a named list
#' of 'eta' and 'rho'; and 'logp', the log-probability corresponding to the
#' parameter values returned
metropolis_step_revisions <- function(
    snapshot_with_true_cases_df,
    reporting_parameters,
    current_revision_parameters,
    logp_current,
    prior_parameters,
    metropolis_parameters,
    max_erroneous=3) {

  proposed_revision_parameters <- propose_revision_parameters(
    current_revision_parameters,
    metropolis_parameters)

  logp_proposed <- observation_process_all_times_logp(
    snapshot_with_true_cases_df=snapshot_with_true_cases_df,
    reporting_parameters=reporting_parameters,
    revision_parameters=proposed_revision_parameters,
    max_erroneous=max_erroneous
  ) + prior_revision_parameters(proposed_revision_parameters,
                                prior_parameters) +
    log_jacobian_revision_parameters(proposed_revision_parameters)

  list_parameters_logp <- accept_reject(
    logp_current, logp_proposed,
    current_revision_parameters,
    proposed_revision_parameters)

  list(
    revision_parameters=list_parameters_logp$parameter,
    logp=list_parameters_logp$logp
  )
}

#' Runs MCMC or optimisation to estimate Rt, cases and reporting parameters
#'
#' @param niterations number of MCMC iterations to run or number of iterative maximisations to run
#' @param snapshot_with_Rt_index_df a tibble with four columns: time_onset, time_reported,
#' cases_reported, Rt_index; an optional fifth column can be provided named reporting_piece_index,
#' which specifies which reporting distribubtion each onset time corresponds to
#' @param priors a named list with: 'Rt', 'reporting', 'max_cases' (and optionally
#' 'overdispersion' if using a negative binomial model). These take
#' the form: 'Rt' is a named list with elements 'shape' and 'rate' describing
#' the gamma prior for each Rt; 'reporting' is a named list with elements
#' 'mean_mu', 'mean_sigma', 'sd_mu', 'sd_sigma' representing the gamma
#' prior parameters for the mean and sd parameters of the reporting parameters
#' (itself described by a gamma distribution); max cases controls the upper
#' limit of the discrete uniform distribution representing the prior on true
#' cases; 'overdispersion' is a named list specifying the mean and sd of a gamma
#' prior on this parameter. If the list contains an element 'Rt_changepoints'
#' (which can be an empty list), the times at which Rt changes are treated as
#' unknown and are sampled (see \code{sample_Rt_changepoints}): the 'Rt_index'
#' column of the data and initial_Rt are then used only to start the sampler.
#' Its optional elements are: 'alpha' and 'beta', the parameters of the beta
#' prior on the probability that Rt changes on each onset day (which default
#' to 1 and the number of onset days, see \code{changepoints_logp_prior});
#' 'min_length', the minimum number of days between changes (default 7); and
#' 'grid', the Rt values over which Rt is integrated numerically under a
#' negative binomial model (see \code{Rt_grid}). The 'Rt' element of the
#' result then has a row for each onset time in each iteration, and an element
#' 'Rt_changepoints' gives the number of changes in each iteration.
#' Likewise, if the list contains an element 'reporting_changepoints' (which
#' can be an empty list), the onset times at which the reporting delay changes
#' are unknown and are sampled (see \code{sample_reporting_changepoints}): any
#' 'reporting_piece_index' column in the data is ignored and the sampler
#' starts from a single segment. Its optional elements are described in
#' \code{reporting_changepoint_settings}. The result then has elements
#' 'reporting_days' (the delay parameters for each onset time in each
#' iteration) and 'reporting_changepoints' (the number of changes and the
#' centre and spread of the delay parameters of segments in each iteration).
#' @inheritParams sample_Rt_single_piece
#' @inheritParams max_uncertain_days
#' @param reporting_metropolis_parameters named list of 'mean_step', 'sd_step' containing
#' step sizes for Metropolis step
#' @param maximise whether to estimate MAP values of parameters (if true) or
#' sample parameter values using MCMC (if false). By default this is false.
#' @param initial_cases_true a tibble with two columns: "time_onset" and "cases_true", which represents initial
#' estimates of the true number of cases with each onset time.
#' @param initial_reporting_parameters provides initial guesses of the mean and
#' sd of the reporting delay distribution(s). These can be either a named with two named elements
#' ('mean', 'sd') for a time-invariant reporting delay or a tibble with three columns:
#' 'reporting_piece_index', 'mean', 'sd' (where the number of indices corresponds to the number
#' provided in the data frame).
#' @param initial_Rt initial guess of the Rt values in each of the piecewise segments.
#' Provided in the form of a tibble with columns: 'Rt_index' and 'Rt'
#' @param print_to_screen prints progress of MCMC sampling to screen. Defaults to true.
#' @param initial_overdispersion the initial value of the overdispersion parameter if
#' assuming a negative binomial sampling model (default to 5).
#' @inheritParams metropolis_step_overdispersion
#' @param condition_on_future_cases if TRUE (the default) each true case count
#' is drawn conditional on the current case counts at all other onset times,
#' including later ones which depend on it through the renewal equation. If
#' FALSE each case count is drawn conditional only on case counts with earlier
#' onset times, which was the behaviour in versions before 2.0.
#' @param initial_revision_parameters named list of 'eta' and 'rho' providing
#' initial guesses of the parameters governing downward revisions of reported
#' case counts (see \code{observation_process_revisions_logp}). If NULL (the
#' default), reported case counts are assumed never to be revised downwards. If
#' supplied, priors must contain an element 'revisions' (see
#' \code{prior_revision_parameters}) and can optionally contain 'max_erroneous'
#' (which defaults to 3).
#' If cases can be misdated the list also contains 'pi', 'shifts' and
#' 'shift_probs' (see \code{unrecorded_prob}) and priors$revisions contains
#' 'pi_alpha' and 'pi_beta'. Batches of erroneous records are allowed for if
#' priors contains an element 'batches': a named list of 'beta', 'zeta',
#' 'phi_mu', 'phi_sigma' and 'rho_batch' (see \code{sample_batches}), which
#' can also contain priors for these (see \code{sample_batch_parameters}),
#' 'max_batch' (default 30) and 'min_total_fall' (see \code{batch_candidates}).
#' Initial batches can be given as an element 'batches' of this list; by
#' default they are chosen from falls in reported cases which coincide.
#' Three extensions to the model of reporting are switched on through further
#' elements of this list: 'tail' (a named list of initial 'weight' and 'mean'
#' of an exponential tail to the delay distribution); 'weekday' (TRUE for an
#' effect of the day of the week of each report, with 'weekday_first' giving
#' the day of the week, from 1 to 7, of time 1); and 'report_effects' (TRUE for
#' a random effect of each report). See \code{log_undetected_days}. Their
#' priors can be set through priors$reporting_effects (see
#' \code{metropolis_step_tail}, \code{metropolis_step_weekday} and
#' \code{sample_report_sd}).
#' @param revision_metropolis_parameters named list of 'eta_step' and 'rho_step'
#' containing step sizes for Metropolis step for revision parameters (see
#' \code{propose_revision_parameters})
#' @param cases_before a vector of known true case counts for the days
#' immediately before the first onset time in the data, from the earliest to
#' the latest (optional). These enter the renewal equation as the history of
#' cases for the first onset times, so that Rt there is not overstated, but
#' are not themselves estimated and nothing is assumed about how they were
#' reported. This is intended for series which begin after an outbreak has
#' started, where counts for earlier onset times are known well enough but
#' the data needed to model their reporting are not available.
#' @param n_parameter_steps the number of Metropolis steps taken in each
#' iteration for each group of continuous parameters updated in that way: the
#' reporting delay parameters of each segment (if their changes are sampled),
#' the revision parameters, the extensions to the model of reporting and the
#' overdispersion parameter. These steps are cheap beside the update of the
#' true cases, so taking several improves mixing at little cost.
#' @param deaths an optional named list which adds a model of deaths among
#' cases: 'snapshot' is a tibble with columns time_onset, time_reported and
#' deaths_reported (the number of the cases with each onset time shown as
#' having died in each report), and 'settings' is an optional named list of
#' settings (see \code{deaths_settings}). The true fatal cases of each onset
#' time are binomial given the true cases, with a case fatality ratio which
#' can change at unknown onset times. A fatal case is shown as dead once it
#' has been reported and its death recorded, the latter after its own delay
#' from onset, which can also change. The results then include "deaths" (the
#' true fatal cases, case fatality ratio and delay for each onset time) and
#' "deaths_parameters". The list can also contain 'uncorrected': a tibble
#' with columns time_onset, time_reported, cases_reported and deaths_reported
#' holding reports made before a report which corrected the records of those
#' onset times (see \code{split_corrected_reports}); these inform the delays
#' and the parameters of wrong records only, fatal cases dated by death being
#' ignored in them, and the results then include "deaths_uncorrected". Wrong records (initial_revision_parameters with 'eta'
#' and 'rho') and effects of the day of the week ('weekday') are allowed for.
#' This requires priors$reporting_changepoints.
#' @param progress_file an optional path of a file to which progress is
#' written (the number of iterations completed and the time taken), so that
#' a chain run in the background or in parallel can be followed. The file is
#' overwritten each time.
#' @param progress_every the number of iterations between each time progress
#' is written. Defaults to 50.
#' @param fixed_cases if TRUE, the true cases are held at initial_cases_true
#' throughout and everything else is sampled given them. Defaults to FALSE.
#' @param return_state if TRUE, the result has an element 'state' from which
#' sampling can be taken up again, and 'iterations_done'. Defaults to FALSE.
#' @param resume_from the 'state' returned by an earlier call. Sampling
#' carries on from where that call stopped, as if it had not stopped: the
#' draws are those which the earlier call would have gone on to make. All
#' other arguments must be as in the earlier call, apart from niterations
#' (the number of further iterations), time_limit and those for progress;
#' the starting values are not used. The iterations returned are numbered
#' on from those of the earlier call.
#' @param time_limit an optional number of minutes after which sampling
#' stops at the end of the iteration under way, returning the draws made so
#' far (and the state, as with return_state).
#' @return a named list of three tibbles: "cases", "Rt" and "reporting" which contain estimates of the model parameters;
#' and "revisions" if initial_revision_parameters is supplied
#' @importFrom rlang .data
#' @importFrom methods is
#' @export
mcmc_single <- function(
  niterations,
  snapshot_with_Rt_index_df,
  priors,
  serial_parameters,
  initial_cases_true,
  initial_reporting_parameters,
  initial_Rt,
  reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
  serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE, print_to_screen=TRUE,
  initial_overdispersion=5,
  is_negative_binomial=FALSE,
  overdispersion_metropolis_sd=0.5,
  condition_on_future_cases=TRUE,
  initial_revision_parameters=NULL,
  revision_metropolis_parameters=list(eta_step=0.5, rho_step=0.5),
  cases_before=NULL,
  n_parameter_steps=1,
  deaths=NULL,
  progress_file=NULL,
  progress_every=50,
  fixed_cases=FALSE,
  return_state=FALSE,
  resume_from=NULL,
  time_limit=NULL) {

  if(!is.null(time_limit) && (!is.numeric(time_limit) || time_limit <= 0))
    stop("time_limit must be a positive number of minutes.")
  if(!is.null(resume_from) &&
     !all(c("variables", "seed", "iterations", "check") %in% names(resume_from)))
    stop("resume_from must be the 'state' returned by an earlier call.")
  if(fixed_cases & (maximise | !condition_on_future_cases))
    stop("fixed_cases needs maximise FALSE and condition_on_future_cases TRUE.")
  if(!is.null(progress_file) && progress_every < 1)
    stop("progress_every must be at least one.")
  time_started <- Sys.time()

  if(n_parameter_steps < 1)
    stop("n_parameter_steps must be at least one.")

  is_deaths <- !is.null(deaths)
  if(is_deaths) {
    if(maximise)
      stop("Maximisation with a model of deaths is not yet allowed.")
    if(!"snapshot" %in% names(deaths))
      stop("deaths must contain 'snapshot': a tibble of reported deaths.")
    if(!"reporting_changepoints" %in% names(priors))
      stop("With a model of deaths, priors$reporting_changepoints must be supplied (changes in the reporting delay can be switched off by setting its n_partition_moves to zero).")
    if(any(c("report_effects", "pi", "batches") %in%
           names(initial_revision_parameters)) | "batches" %in% names(priors))
      stop("With a model of deaths, wrong records, a heavier tail of the reporting delay and effects of the day of the week are allowed for, but not random effects of reports, misdating or batches.")
    if(!condition_on_future_cases)
      stop("condition_on_future_cases must be TRUE with a model of deaths.")
    deaths_n_days <- length(unique(snapshot_with_Rt_index_df$time_onset))
    deaths_options <- deaths$settings
    if(is.null(deaths_options))
      deaths_options <- list()
    if(is.null(deaths_options$n_parameter_steps))
      deaths_options$n_parameter_steps <- n_parameter_steps
    deaths_options <- deaths_settings(deaths_n_days, deaths_options)
    deaths_state <- initial_deaths_state(
      dplyr::left_join(
        snapshot_with_Rt_index_df[, c("time_onset", "time_reported", "cases_reported")],
        deaths$snapshot[, c("time_onset", "time_reported", "deaths_reported")],
        by=c("time_onset", "time_reported")),
      initial_cases_true, deaths_options, cases_before)
    # with fatal cases dated by death, the true cases by day of onset differ
    # from those recorded for each day
    initial_cases_true <- dplyr::tibble(time_onset=deaths_state$onset_times,
                                        cases_true=deaths_state$cases_onset)
    is_uncorrected <- !is.null(deaths$uncorrected)
    if(is_uncorrected) {
      deaths_state$uncorrected <- initial_uncorrected(deaths$uncorrected,
                                                      deaths_state, deaths_options)
      uncorrected_samples <- vector("list", niterations)
    }
    deaths_samples <- vector("list", niterations)
    deaths_parameters_samples <- vector("list", niterations)
  }

  if(!is.null(cases_before)) {
    if(any(is.na(cases_before)) || any(cases_before < 0))
      stop("cases_before must be non-negative case counts.")
    if(maximise)
      stop("Maximisation with known cases before the first onset time is not yet allowed.")
  }

  cnames <- colnames(snapshot_with_Rt_index_df)
  expected_names <- c("time_onset", "time_reported",
                      "cases_reported", "Rt_index")

  if(sum(cnames %in% expected_names) != 4)
    stop("Incorrect column names in snapshot_with_Rt_index_df")

  if(length(cnames) != 4 & length(cnames) != 5)
    stop("There should be either four or five columns in snapshot_with_Rt_index_df")

  if(length(cnames) == 5 & !("reporting_piece_index" %in% cnames))
    stop("The reporting delay indices must be provided through a column named 'reporting_piece_index'")

  # assume if column not supplied then reporting delay same throughout
  if(!"reporting_piece_index" %in% cnames)
    snapshot_with_Rt_index_df$reporting_piece_index <- 1

  if(initial_overdispersion <= 0)
    stop("Initial overdispersion value (for a negative binomial renewal model) must be positive.")
  overdispersion_current <- initial_overdispersion

  if(is_negative_binomial)
    if(!"overdispersion" %in% names(priors))
      stop("If negative binomial model is used a prior must be provided for the overdispersion parameter.")

  if(is_negative_binomial)
    if(maximise)
      stop("Maximisation with a negative binomial renewal model is not yet allowed.")

  is_revisions <- !is.null(initial_revision_parameters)
  revision_current <- initial_revision_parameters
  max_erroneous <- 3
  if(is_revisions) {
    if(!"revisions" %in% names(priors))
      stop("If revisions to reported cases are allowed a prior must be provided for the revision parameters.")
    if(maximise)
      stop("Maximisation when revisions to reported cases are allowed is not yet allowed.")
    if(revision_current$eta <= 0 | revision_current$rho <= 0 | revision_current$rho >= 1)
      stop("Initial revision parameters must satisfy: eta > 0 and 0 < rho < 1.")
    if("max_erroneous" %in% names(priors))
      max_erroneous <- priors$max_erroneous
    # larger falls in reported cases are possible if there are batches
    if(!"batches" %in% names(priors))
      check_revisions_possible(snapshot_with_Rt_index_df, max_erroneous)
    if(misdating_prob(revision_current) > 0) {
      if(!condition_on_future_cases)
        stop("condition_on_future_cases must be TRUE if cases can be misdated.")
      if(revision_current$pi >= 1)
        stop("Initial misdating probability must be less than one.")
      if(is.null(revision_current$shifts) |
         length(revision_current$shifts) != length(revision_current$shift_probs))
        stop("If cases can be misdated, 'shifts' and 'shift_probs' of equal length must be given.")
      if(any(revision_current$shifts == 0) |
         abs(sum(revision_current$shift_probs) - 1) > 1e-8)
        stop("Shifts must be non-zero and their probabilities must sum to one.")
    }
  }

  # extensions to the model of reporting
  is_tail <- is_revisions && !is.null(revision_current$tail)
  is_weekday <- is_revisions && isTRUE(revision_current$weekday)
  is_report_effects <- is_revisions && isTRUE(revision_current$report_effects)
  is_reporting_extensions <- is_tail | is_weekday | is_report_effects
  if(is_reporting_extensions) {
    priors_reporting <- list(tail_weight_alpha=1, tail_weight_beta=3,
                             tail_mean_mu=15, tail_mean_sigma=10,
                             weekday_sd=1,
                             report_sd_shape=3, report_sd_scale=0.2)
    if("reporting_effects" %in% names(priors))
      priors_reporting[names(priors$reporting_effects)] <- priors$reporting_effects
    steps_reporting <- list(tail_weight_step=0.3, tail_mean_step=0.2,
                            weekday_step=0.05, report_blocks=2)
    names_steps <- intersect(names(steps_reporting),
                             names(revision_metropolis_parameters))
    steps_reporting[names_steps] <- revision_metropolis_parameters[names_steps]
  }
  if(is_weekday | is_report_effects) {
    revision_current$report_days <- sort(unique(
      snapshot_with_Rt_index_df$time_reported))
    if(is.null(revision_current$weekday_first))
      revision_current$weekday_first <- 1
    if(is_weekday) {
      weekdays_present <- sort(unique(weekday_of_days(
        revision_current$report_days, revision_current$weekday_first)))
      revision_current$weekday_log_effects <- stats::setNames(
        rep(0, length(weekdays_present)), weekdays_present)
    }
    if(is_report_effects) {
      revision_current$report_log_effects <- rep(
        0, length(revision_current$report_days))
      revision_current$report_sd <- sqrt(
        priors_reporting$report_sd_scale / (priors_reporting$report_sd_shape - 1))
    }
    revision_current <- update_report_multipliers(revision_current)
  }
  # the delay to a death being recorded has its own effects of the day of
  # the week
  if(is_deaths) {
    if(is.null(deaths_options$recording_weekday))
      deaths_options$recording_weekday <- is_weekday
    if(isTRUE(deaths_options$recording_weekday))
      deaths_state$recording_effects <- initial_recording_effects(
        snapshot_with_Rt_index_df$time_reported,
        ifelse(is.null(initial_revision_parameters$weekday_first), 1,
               initial_revision_parameters$weekday_first))
    if(!identical(deaths_options$recording_tail, FALSE))
      deaths_state$recording_effects$tail <- deaths_options$recording_tail$initial
    # a death is recorded after the time from onset to death and then a delay
    if(!is.null(deaths_options$onset_to_death))
      deaths_state$recording_effects$before <- deaths_options$onset_to_death[
        "lag_probabilities"]
    if(is_reporting_extensions) {
      if(is.null(deaths$settings$weekday_sd))
        deaths_options$weekday_sd <- priors_reporting$weekday_sd
      if(is.null(deaths$settings$weekday_step))
        deaths_options$weekday_step <- steps_reporting$weekday_step
      deaths_options$tail_prior <- priors_reporting
      deaths_options$tail_steps <- steps_reporting
    }
  }

  # batches of erroneous records
  is_batches <- "batches" %in% names(priors)
  if(is_batches) {
    if(!is_revisions)
      stop("Batches are only allowed for if initial_revision_parameters is supplied.")
    names_batch <- c("beta", "zeta", "phi_mu", "phi_sigma", "rho_batch")
    if(!all(names_batch %in% names(priors$batches)))
      stop("priors$batches must contain: 'beta', 'zeta', 'phi_mu', 'phi_sigma' and 'rho_batch'.")
    batch_parameters <- priors$batches[names_batch]
    revision_current$max_batch <- 30
    if("max_batch" %in% names(priors$batches))
      revision_current$max_batch <- priors$batches$max_batch
    min_total_fall <- max_erroneous + 1
    if("min_total_fall" %in% names(priors$batches))
      min_total_fall <- priors$batches$min_total_fall
    candidate_batches <- batch_candidates(snapshot_with_Rt_index_df, min_total_fall)
    # start only with those batches without which reported cases are
    # impossible: candidates containing a fall exceeding max_erroneous
    if(is.null(revision_current$batches)) {
      large_falls_df <- snapshot_with_Rt_index_df %>%
        dplyr::group_by(.data$time_onset) %>%
        dplyr::arrange(.data$time_reported, .by_group=TRUE) %>%
        dplyr::mutate(fall=c(0, -diff(.data$cases_reported))) %>%
        dplyr::ungroup() %>%
        dplyr::filter(.data$fall > max_erroneous)
      is_needed <- purrr::map_lgl(seq_len(nrow(candidate_batches)), function(k) {
        any(large_falls_df$time_reported == candidate_batches$time_detected[k] &
              large_falls_df$time_onset >= candidate_batches$start[k] &
              large_falls_df$time_onset < candidate_batches$start[k] +
              candidate_batches$width[k])
      })
      revision_current$batches <- initial_batches(candidate_batches[is_needed, ],
                                                  batch_parameters)
    }
  }

  # changes in Rt at unknown times
  is_Rt_changepoints <- "Rt_changepoints" %in% names(priors)
  if(is_Rt_changepoints) {
    if(maximise)
      stop("Maximisation when the times of changes in Rt are unknown is not yet allowed.")
    Rt_changepoints <- list(
      alpha=1,
      beta=length(unique(snapshot_with_Rt_index_df$time_onset)),
      min_length=7,
      grid=Rt_grid())
    Rt_changepoints[names(priors$Rt_changepoints)] <- priors$Rt_changepoints
    if(Rt_changepoints$alpha <= 0 | Rt_changepoints$beta <= 0)
      stop("The parameters of the prior on the probability of a change in Rt must be positive.")
    if(Rt_changepoints$min_length < 1)
      stop("The minimum length of a segment of constant Rt must be at least one day.")
    Rt_p_change <- Rt_changepoints$alpha /
      (Rt_changepoints$alpha + Rt_changepoints$beta)
    Rt_days_samples <- vector("list", niterations)
    Rt_changepoints_samples <- vector("list", niterations)
  }

  reporting_current <- initial_reporting_parameters
  if(methods::is(reporting_current, "list")) # if only a single list provided
    reporting_current <- dplyr::tibble(
      reporting_piece_index=1,
      mean=reporting_current$mean,
      sd=reporting_current$sd
    )
  num_reporting_parameters <- max(reporting_current$reporting_piece_index)

  # changes in the reporting delay at unknown onset times
  is_reporting_changepoints <- "reporting_changepoints" %in% names(priors)
  if(is_reporting_changepoints) {
    if(maximise)
      stop("Maximisation when the times of changes in the reporting delay are unknown is not yet allowed.")
    onset_times_all <- sort(unique(snapshot_with_Rt_index_df$time_onset))
    reporting_changepoints <- reporting_changepoint_settings(
      length(onset_times_all),
      c(priors$reporting_changepoints,
        reporting_metropolis_parameters[setdiff(
          names(reporting_metropolis_parameters),
          names(priors$reporting_changepoints))]))
    # the sampler starts from a single segment
    reporting_current <- reporting_current[1, ]
    reporting_current$reporting_piece_index <- 1
    snapshot_with_Rt_index_df$reporting_piece_index <- 1
    if(!"n_segment_steps" %in% names(priors$reporting_changepoints))
      reporting_changepoints$n_segment_steps <- n_parameter_steps
    reporting_state <- initial_reporting_changepoints(reporting_current)
    reporting_days_samples <- vector("list", niterations)
    reporting_changepoints_samples <- vector("list", niterations)
  }

  df_running <- snapshot_with_Rt_index_df %>%
    dplyr::left_join(initial_cases_true, by = "time_onset") %>%
    dplyr::left_join(initial_Rt, by = "Rt_index")

  num_Rts <- nrow(initial_Rt)
  Rt_samples <- matrix(nrow = num_Rts * niterations,
                       ncol = 3)
  num_cases_true <- unique(df_running$time_onset)

  if(print_to_screen) {
    progress_bar <- utils::txtProgressBar(
      min = 0,
      max = niterations,
      style = 3,
      width = niterations, # Needed to avoid multiple printings
      char = "=")

    init <- numeric(niterations)
    end <- numeric(niterations)
  }

  max_cases <- priors$max_cases

  # taking sampling up where an earlier call stopped: everything which is
  # carried from one iteration to the next is put back as it was, with the
  # state of the random number generator
  state_check <- c(nrow(snapshot_with_Rt_index_df),
                   sum(snapshot_with_Rt_index_df$cases_reported),
                   if(is_deaths) sum(deaths$snapshot$deaths_reported))
  iteration_offset <- 0
  if(!is.null(resume_from)) {
    if(!isTRUE(all.equal(resume_from$check, state_check)))
      stop("resume_from is the state of a fit to different reports.")
    list2env(resume_from$variables, envir=environment())
    assign(".Random.seed", resume_from$seed, envir=globalenv())
    iteration_offset <- resume_from$iterations
  }
  iterations_done <- niterations

  k <- 1
  for(i in 1:niterations) {

    if(print_to_screen)
      init[i] <- Sys.time()

    # sample which of the cases in each report are fatal
    if(is_deaths) {
      deaths_running_df <- df_running %>%
        dplyr::select("time_onset", "cases_true", "reporting_piece_index") %>%
        unique() %>%
        dplyr::left_join(reporting_current[, c("reporting_piece_index", "mean", "sd")],
                         by="reporting_piece_index")
      deaths_running_df <- deaths_running_df[
        match(deaths_state$onset_times, deaths_running_df$time_onset), ]
      deaths_state <- deaths_paths_step(deaths_state, deaths_running_df$cases_true,
                                        deaths_running_df, revision_current,
                                        max_erroneous)
      # the reports made before a correction, for the onset times corrected
      deaths_state <- deaths_uncorrected_step(deaths_state, deaths_running_df,
                                              revision_current, max_erroneous,
                                              max_cases)
      # which fatal cases are recorded for their day of death, which moves
      # true cases between days of onset
      if(!is.null(deaths_state$dated)) {
        deaths_Rt_df <- unique(df_running[, c("time_onset", "Rt")])
        deaths_dated <- deaths_dated_step(
          deaths_state, deaths_running_df$cases_true,
          dated_renewal(
            deaths_Rt_df$Rt[match(deaths_state$onset_times, deaths_Rt_df$time_onset)],
            overdispersion_current, is_negative_binomial, serial_parameters,
            serial_max, cases_before))
        deaths_state <- deaths_dated_fraction_step(deaths_dated$state)
        df_running$cases_true <- deaths_dated$cases_true[
          match(df_running$time_onset, deaths_state$onset_times)]
      }
    }

    # sample incidence
    Rt_current <- df_running %>%
      dplyr::select(time_onset, Rt) %>%
      unique()
    Rt_function <- stats::approxfun(Rt_current$time_onset,
                                    Rt_current$Rt)
    df_current <- df_running %>%
      dplyr::select(time_onset, time_reported, cases_reported, reporting_piece_index)

    cases_true_current <- NULL
    if(condition_on_future_cases)
      cases_true_current <- df_running %>%
        dplyr::select(time_onset, cases_true) %>%
        unique()

    # with fixed_cases the true cases stay at their initial values
    df_temp <- if(fixed_cases) dplyr::left_join(
      df_current,
      dplyr::rename(cases_true_current, cases_estimated="cases_true"),
      by="time_onset")
    else sample_cases_history(
      df_current, max_cases,
      Rt_function, serial_parameters, reporting_current,
      p_gamma_cutoff=p_gamma_cutoff,
      maximise=maximise,
      kappa=overdispersion_current,
      is_negative_binomial=is_negative_binomial,
      cases_true_current=cases_true_current,
      future_max=serial_max,
      revision_parameters=revision_current,
      max_erroneous=max_erroneous,
      cases_before=cases_before,
      f_observation_logp_onset=if(is_deaths)
        function(time) deaths_cases_logp_function(deaths_state, time,
                                                  revision_current),
      cases_offset=if(is_deaths) deaths_cases_offset(deaths_state))

    # sample the true fatal cases, the delay from onset to a death being
    # recorded and the case fatality ratio
    if(is_deaths) {
      cases_by_onset_df <- unique(df_temp[, c("time_onset", "cases_estimated")])
      deaths_cases_true <- cases_by_onset_df$cases_estimated[
        match(deaths_state$onset_times, cases_by_onset_df$time_onset)]
      deaths_state <- deaths_true_step(deaths_state, deaths_cases_true)
      # what acts on the reporting delay or on wrong records, then what acts
      # on the delay to a death being recorded
      deaths_effects <- deaths_reporting_effects_step(
        deaths_state, deaths_cases_true, deaths_running_df, revision_current,
        deaths_options, priors$revisions,
        ifelse(is.null(revision_metropolis_parameters$rho_step), 0.5,
               revision_metropolis_parameters$rho_step))
      deaths_state <- deaths_effects$state
      revision_current <- deaths_effects$revision_parameters
      deaths_state <- deaths_parameters_step(deaths_state, deaths_cases_true,
                                             deaths_options)
      deaths_state <- deaths_refresh_survivals(deaths_state, deaths_running_df,
                                               revision_current)
      deaths_tmp <- deaths_draws(deaths_state, i)
      deaths_samples[[i]] <- deaths_tmp$deaths
      deaths_parameters_samples[[i]] <- deaths_tmp$deaths_parameters
      if(is_uncorrected)
        uncorrected_samples[[i]] <- uncorrected_draws(deaths_state, i)
      if(is_revisions) {
        revision_tmp <- dplyr::tibble(eta=revision_current$eta,
                                      rho=revision_current$rho, pi=0, iteration=i)
        revision_samples <- if(i == 1) revision_tmp
        else dplyr::bind_rows(revision_samples, revision_tmp)
      }
      if(is_reporting_extensions) {
        reporting_effects_tmp <- dplyr::tibble(
          iteration=i,
          tail_weight=ifelse(is_tail, revision_current$tail$weight, NA),
          tail_mean=ifelse(is_tail, revision_current$tail$mean, NA), report_sd=NA)
        for(weekday in names(revision_current$weekday_log_effects))
          reporting_effects_tmp[[paste0("weekday_", weekday)]] <-
            exp(revision_current$weekday_log_effects[[weekday]])
        reporting_effects_samples <- if(i == 1) reporting_effects_tmp
        else dplyr::bind_rows(reporting_effects_samples, reporting_effects_tmp)
      }
    }

    # sample Rt
    cases_history_df <- df_running %>%
      dplyr::select(time_onset, Rt_index) %>%
      dplyr::left_join(df_temp, by = "time_onset", relationship = "many-to-many") %>%
      dplyr::select(-c("cases_reported", "time_reported")) %>%
      dplyr::rename(cases_true=cases_estimated) %>%
      unique()

    if(is_Rt_changepoints) {
      # the segments within which Rt is constant are sampled along with Rt
      df_Rt_days <- sample_Rt_changepoints(
        cases_history_df %>%
          dplyr::select("time_onset", "cases_true") %>%
          unique(),
        priors$Rt,
        serial_parameters,
        p_change=Rt_p_change,
        min_length=Rt_changepoints$min_length,
        serial_max=serial_max,
        kappa=overdispersion_current,
        is_negative_binomial=is_negative_binomial,
        grid=Rt_changepoints$grid,
        cases_before=cases_before)
      Rt_n_changes <- max(df_Rt_days$Rt_index) - 1
      Rt_p_change <- sample_change_probability(
        Rt_n_changes, nrow(df_Rt_days),
        Rt_changepoints$alpha, Rt_changepoints$beta)
      cases_history_df <- cases_history_df %>%
        dplyr::select(-"Rt_index") %>%
        dplyr::left_join(df_Rt_days[, c("time_onset", "Rt_index")],
                         by="time_onset")
      df_Rt <- df_Rt_days %>%
        dplyr::select("Rt_index", "Rt") %>%
        unique() %>%
        dplyr::mutate(draw_index=1)
      Rt_days_samples[[i]] <- df_Rt_days %>%
        dplyr::mutate(iteration=i)
      Rt_changepoints_samples[[i]] <- dplyr::tibble(
        iteration=i,
        n_changes=Rt_n_changes,
        p_change=Rt_p_change)
    } else {
      df_Rt <- sample_Rt(cases_history_df,
                         priors$Rt,
                         serial_parameters,
                         kappa=overdispersion_current,
                         serial_max=serial_max,
                         ndraws=1,
                         maximise=maximise,
                         is_negative_binomial=is_negative_binomial,
                         Rt_current=unique(df_running[, c("Rt_index", "Rt")]),
                         cases_before=cases_before)
      # nocov start
      if(nrow(df_Rt) != num_Rts)
        stop("Number of Rts outputted not equal to initial Rt dims.")
      # nocov end
      # store Rts
      for(j in 1:num_Rts) {
        Rt_index <- df_Rt$Rt_index[j]
        Rt_temp <- df_Rt$Rt[j]
        Rt_samples[k, ] <- c(i, Rt_index, Rt_temp)
        k <- k + 1
      }
    }

    # sample reporting parameters
    df_temp <- df_temp %>%
      dplyr::rename(cases_true=cases_estimated)
    # recalculated each iteration since it depends on the true cases, which
    # have just been updated
    logp_current_reporting <- NULL
    # not needed with a model of deaths, where the reports of fatal and other
    # cases are handled separately
    if(!maximise & !is_deaths) {
      logp_current_reporting <- observation_process_all_times_logp(
        snapshot_with_true_cases_df=df_temp,
        reporting_parameters=reporting_current,
        revision_parameters=revision_current,
        max_erroneous=max_erroneous
      ) + prior_reporting_parameters(reporting_current,
                                     priors$reporting)
      if(is_revisions & !is.finite(logp_current_reporting))
        stop("Reported cases have zero probability at current parameter values: max_erroneous (or max_batch) may need to be increased.")
    }
    if(is_reporting_changepoints) {
      # the segments within which the delay is constant are sampled too
      reporting_state <- sample_reporting_changepoints(
        reporting_state,
        if(is_deaths)
          deaths_reporting_logp_function(deaths_state, deaths_cases_true,
                                         revision_current)
        else
          reporting_logp_function(df_temp, revision_current, max_erroneous),
        length(onset_times_all),
        priors$reporting,
        reporting_changepoints)
      reporting_current <- reporting_state$parameters %>%
        dplyr::mutate(reporting_piece_index=dplyr::row_number())
      piece_by_onset <- dplyr::tibble(
        time_onset=onset_times_all,
        reporting_piece_index=segment_indices(reporting_state$starts,
                                              length(onset_times_all)))
      df_temp <- df_temp %>%
        dplyr::select(-"reporting_piece_index") %>%
        dplyr::left_join(piece_by_onset, by="time_onset")
      cases_history_df <- cases_history_df %>%
        dplyr::select(-"reporting_piece_index") %>%
        dplyr::left_join(piece_by_onset, by="time_onset")
      # in the form returned by sample_reporting
      list_reporting_logp <- list(
        reporting_parameters=reporting_current,
        logp=sum(reporting_state$logl) +
          prior_reporting_parameters(reporting_current, priors$reporting))
      reporting_days_samples[[i]] <- piece_by_onset %>%
        dplyr::left_join(reporting_current, by="reporting_piece_index") %>%
        dplyr::mutate(iteration=i)
      reporting_changepoints_samples[[i]] <- dplyr::tibble(
        iteration=i,
        n_changes=nrow(reporting_current) - 1,
        centre_mean=reporting_state$centre$mean,
        centre_sd=reporting_state$centre$sd,
        spread=reporting_state$spread)
    } else {
      list_reporting_logp <- sample_reporting(
        snapshot_with_true_cases_df=df_temp,
        current_reporting_parameters=reporting_current,
        logp_current=logp_current_reporting,
        prior_parameters=priors$reporting,
        metropolis_parameters=reporting_metropolis_parameters,
        maximise=maximise,
        ndraws=1,
        revision_parameters=revision_current,
        max_erroneous=max_erroneous)
      reporting_temp <- list_reporting_logp$reporting_parameters
      reporting_temp$draw_index <- NULL
      reporting_current <- reporting_temp
    }

    # sample revision parameters (with a model of deaths, done above)
    if(is_revisions & !is_deaths) {
      # log-probability of observations at reporting parameters just sampled
      logp_observation <- list_reporting_logp$logp -
        prior_reporting_parameters(reporting_current, priors$reporting)
      logp_current_revisions <- logp_observation +
        prior_revision_parameters(revision_current, priors$revisions) +
        log_jacobian_revision_parameters(revision_current)
      for(step in seq_len(n_parameter_steps)) {
        list_revision_logp <- metropolis_step_revisions(
          snapshot_with_true_cases_df=df_temp,
          reporting_parameters=reporting_current,
          current_revision_parameters=revision_current,
          logp_current=logp_current_revisions,
          prior_parameters=priors$revisions,
          metropolis_parameters=revision_metropolis_parameters,
          max_erroneous=max_erroneous)
        revision_current <- list_revision_logp$revision_parameters
        logp_current_revisions <- list_revision_logp$logp
      }

      revision_tmp <- dplyr::tibble(
        eta=revision_current$eta,
        rho=revision_current$rho,
        pi=misdating_prob(revision_current),
        iteration=i)
      if(i == 1)
        revision_samples <- revision_tmp
      else
        revision_samples <- revision_samples %>%
          dplyr::bind_rows(revision_tmp)
    }

    # sample the extensions to the model of reporting
    if(is_reporting_extensions & !is_deaths) {
      # log-probability of observations at the values just sampled
      logl_reporting <- list_revision_logp$logp -
        prior_revision_parameters(revision_current, priors$revisions) -
        log_jacobian_revision_parameters(revision_current)
      for(step in seq_len(n_parameter_steps)) {
        if(is_tail) {
          list_tail <- metropolis_step_tail(
            df_temp, reporting_current, revision_current, logl_reporting,
            priors_reporting, steps_reporting, max_erroneous)
          revision_current <- list_tail$revision_parameters
          logl_reporting <- list_tail$logl
        }
        if(is_weekday) {
          list_weekday <- metropolis_step_weekday(
            df_temp, reporting_current, revision_current, logl_reporting,
            priors_reporting, steps_reporting, max_erroneous)
          revision_current <- list_weekday$revision_parameters
          logl_reporting <- list_weekday$logl
        }
      }
      if(is_report_effects) {
        list_effects <- sample_report_effects(
          df_temp, reporting_current, revision_current, logl_reporting,
          max_erroneous, steps_reporting$report_blocks)
        revision_current <- list_effects$revision_parameters
        revision_current$report_sd <- sample_report_sd(
          revision_current$report_log_effects, priors_reporting)
      }

      reporting_effects_tmp <- dplyr::tibble(
        iteration=i,
        tail_weight=ifelse(is_tail, revision_current$tail$weight, NA),
        tail_mean=ifelse(is_tail, revision_current$tail$mean, NA),
        report_sd=ifelse(is_report_effects, revision_current$report_sd, NA))
      if(is_weekday)
        for(weekday in names(revision_current$weekday_log_effects))
          reporting_effects_tmp[[paste0("weekday_", weekday)]] <-
            exp(revision_current$weekday_log_effects[[weekday]])
      report_effects_tmp <- NULL
      if(is_report_effects)
        report_effects_tmp <- dplyr::tibble(
          iteration=i,
          time_reported=revision_current$report_days,
          effect=exp(revision_current$report_log_effects))
      if(i == 1) {
        reporting_effects_samples <- reporting_effects_tmp
        report_effects_samples <- report_effects_tmp
      } else {
        reporting_effects_samples <- reporting_effects_samples %>%
          dplyr::bind_rows(reporting_effects_tmp)
        report_effects_samples <- report_effects_samples %>%
          dplyr::bind_rows(report_effects_tmp)
      }
    }

    # sample batches and the parameters governing them
    if(is_batches) {
      lists <- observation_lists(df_temp, reporting_current)
      revision_current$batches <- sample_batches_step(
        batches=revision_current$batches,
        onset_times=lists$onset_times,
        l_observations=lists$l_observations,
        cases_by_onset=lists$cases_by_onset,
        l_reporting=lists$l_reporting,
        revision_parameters=revision_current,
        max_erroneous=max_erroneous,
        batch_parameters=batch_parameters,
        candidates=candidate_batches)
      batch_parameters <- sample_batch_parameters(
        revision_current$batches, batch_parameters, priors$batches,
        lists$onset_times)

      batches_tmp <- revision_current$batches %>%
        dplyr::mutate(iteration=i)
      batch_parameters_tmp <- dplyr::tibble(
        beta=batch_parameters$beta,
        zeta=batch_parameters$zeta,
        rho_batch=batch_parameters$rho_batch,
        n_batches=nrow(revision_current$batches),
        iteration=i)
      if(i == 1) {
        batches_samples <- batches_tmp
        batch_parameters_samples <- batch_parameters_tmp
      } else {
        batches_samples <- batches_samples %>%
          dplyr::bind_rows(batches_tmp)
        batch_parameters_samples <- batch_parameters_samples %>%
          dplyr::bind_rows(batch_parameters_tmp)
      }
    }
    reporting_current$iteration <- i
    if(i == 1)
      reporting_samples <- reporting_current
    else
      reporting_samples <- reporting_samples %>% dplyr::bind_rows(reporting_current)

    # update main df used in sampling
    df_running <- df_temp %>%
      dplyr::select(-cases_true) %>%
      dplyr::left_join(cases_history_df, by = c("time_onset", "reporting_piece_index")) %>%
      dplyr::left_join(df_Rt, by = "Rt_index") %>%
      dplyr::select(-draw_index)

    if(is_negative_binomial) {
      df_nb_tmp <- df_running %>%
        dplyr::select(time_onset, cases_true, Rt) %>%
        unique()

      # recalculated each iteration since it depends on the true cases and Rt,
      # which have just been updated
      logp_prior <- dgamma_mean_sd(overdispersion_current,
                                   priors$overdispersion$mean,
                                   priors$overdispersion$sd,
                                   log=TRUE)
      logp_overdispersion <- state_process_nb_logp_all_onsets(
        overdispersion_current, df_nb_tmp, serial_parameters,
        serial_max=serial_max, cases_before=cases_before
      ) + logp_prior

      for(step in seq_len(n_parameter_steps)) {
        list_parameter_logp <- metropolis_step_overdispersion(
          overdispersion_current,
          logp_overdispersion,
          df_nb_tmp,
          serial_parameters,
          priors$overdispersion,
          overdispersion_metropolis_sd,
          serial_max=serial_max,
          cases_before=cases_before
        )
        overdispersion_current <- list_parameter_logp$overdispersion
        logp_overdispersion <- list_parameter_logp$logp
      }

      overdipersion_tmp <- dplyr::tibble(
        overdispersion=overdispersion_current,
        iteration=i)

      if(i == 1)
        overdispersion_samples <- overdipersion_tmp
      else
        overdispersion_samples <- overdispersion_samples %>%
          dplyr::bind_rows(overdipersion_tmp)
    }

    # store cases
    cases_history_df <- cases_history_df %>%
      dplyr::select(-Rt_index) %>%
      dplyr::mutate(iteration=i)
    if(i == 1) {
      cases_history_samples <- cases_history_df
    } else {
      cases_history_samples <- cases_history_samples %>%
        dplyr::bind_rows(cases_history_df)
    }

    if(!is.null(progress_file) &&
       (i %% progress_every == 0 | i == niterations))
      writeLines(sprintf(
        "%d of %d iterations, %.1f minutes", i, niterations,
        as.numeric(Sys.time() - time_started, units="mins")), progress_file)

    if(print_to_screen) {
      end[i] <- Sys.time()
      utils::setTxtProgressBar(progress_bar, i)
      time <- round(lubridate::seconds_to_period(sum(end - init)), 0)

      # Estimated remaining time based on the
      # mean time that took to run the previous iterations
      est <- niterations * (mean(end[end != 0] - init[init != 0])) - time
      remaining <- round(lubridate::seconds_to_period(est), 0)

      cat(paste(" // Execution time:", tolower(as.character(time)),
                " // Estimated time remaining:", tolower(as.character(remaining))), "")
    }

    if(!is.null(time_limit) &&
       as.numeric(Sys.time() - time_started, units="mins") >= time_limit) {
      iterations_done <- i
      break
    }
  }

  Rt_samples <- Rt_samples[seq_len(k - 1), , drop=FALSE] %>%
    as.data.frame()
  colnames(Rt_samples) <- c("iteration", "Rt_index", "Rt")
  if(is_Rt_changepoints)
    Rt_samples <- dplyr::bind_rows(Rt_days_samples) %>%
      dplyr::relocate("iteration", "time_onset", "Rt_index", "Rt")

  cases_samples <- cases_history_samples %>%
    dplyr::select(-reporting_piece_index)

  reporting_samples <- reporting_samples %>%
    dplyr::relocate("reporting_piece_index", "mean", "sd", "iteration")

  list_results <-list(
    cases=cases_samples,
    Rt=Rt_samples,
    reporting=reporting_samples)

  if(is_negative_binomial)
    list_results$overdispersion <- overdispersion_samples

  if(is_revisions)
    list_results$revisions <- revision_samples

  if(is_batches) {
    list_results$batches <- batches_samples
    list_results$batch_parameters <- batch_parameters_samples
  }

  if(is_Rt_changepoints)
    list_results$Rt_changepoints <- dplyr::bind_rows(Rt_changepoints_samples)

  if(is_reporting_changepoints) {
    list_results$reporting_days <- dplyr::bind_rows(reporting_days_samples) %>%
      dplyr::relocate("iteration", "time_onset", "reporting_piece_index",
                      "mean", "sd")
    list_results$reporting_changepoints <- dplyr::bind_rows(
      reporting_changepoints_samples)
  }

  if(is_reporting_extensions) {
    list_results$reporting_effects <- reporting_effects_samples
    if(is_report_effects)
      list_results$report_effects <- report_effects_samples
  }

  if(is_deaths) {
    list_results$deaths <- dplyr::bind_rows(deaths_samples)
    list_results$deaths_parameters <- dplyr::bind_rows(deaths_parameters_samples)
    if(is_uncorrected)
      list_results$deaths_uncorrected <- dplyr::bind_rows(uncorrected_samples)
  }

  # iterations are numbered on from those of the call being taken up
  if(iteration_offset > 0)
    list_results <- lapply(list_results, function(df) {
      df$iteration <- df$iteration + iteration_offset
      df
    })

  if(return_state | !is.null(resume_from) | !is.null(time_limit)) {
    # everything carried from one iteration to the next: all that the
    # function holds apart from its arguments, the draws and its functions
    variables <- as.list(environment(), all.names=TRUE)
    names_left_out <- c(
      names(formals(mcmc_single)), "list_results", "variables", "progress_bar",
      "init", "end", "time", "est", "remaining", "time_started", "i", "k",
      "state_check", "iteration_offset", "iterations_done",
      grep("_samples$", names(variables), value=TRUE))
    variables <- variables[setdiff(names(variables), names_left_out)]
    variables <- variables[!vapply(variables, is.function, logical(1))]
    list_results$state <- list(
      variables=variables,
      seed=get(".Random.seed", envir=globalenv()),
      iterations=iteration_offset + iterations_done,
      check=state_check)
    list_results$iterations_done <- iterations_done
  }

  list_results
}

#' Combines Markov chains across multiple runs of mcmc_single
#'
#' @param list_of_results a list of results, where each element is a result of running mcmc_single
#' @inheritParams mcmc_single
#' @return a named list of three tibbles: "cases", "Rt" and "reporting" (and if running
#' a negative binomial model an additional 'overdispersion' element) which
#' contain estimates of the model parameters with chain index  included
combine_chains <- function(list_of_results, is_negative_binomial=FALSE) {

  for(i in seq_along(list_of_results)) {
    res <- list_of_results[[i]]
    cases_df <- res$cases
    Rt_df <- res$Rt
    reporting_df <- res$reporting
    cases_df$chain <- i
    Rt_df$chain <- i
    reporting_df$chain <- i

    if(is_negative_binomial) {
      overdispersion_df <- res$overdispersion
      overdispersion_df$chain <- i
    }
    is_revisions <- "revisions" %in% names(res)
    if(is_revisions) {
      revisions_df <- res$revisions
      revisions_df$chain <- i
      if(i == 1)
        revisions_overall <- revisions_df
      else
        revisions_overall <- revisions_overall %>% dplyr::bind_rows(revisions_df)
    }
    if(i == 1) {
      cases_overall <- cases_df
      Rt_overall <- Rt_df
      reporting_overall <- reporting_df
      if(is_negative_binomial) {
        overdispersion_overall <- overdispersion_df
      }
    } else {
      cases_overall <- cases_overall %>% dplyr::bind_rows(cases_df)
      Rt_overall <- Rt_overall %>% dplyr::bind_rows(Rt_df)
      reporting_overall <- reporting_overall %>% dplyr::bind_rows(reporting_df)
      if(is_negative_binomial) {
        overdispersion_overall <- overdispersion_overall %>% dplyr::bind_rows(overdispersion_df)
      }
    }
  }
  list_combined <- list(
    cases=cases_overall,
    Rt=Rt_overall,
    reporting=reporting_overall
  )

  if(is_negative_binomial)
    list_combined$overdispersion <- overdispersion_overall

  if(is_revisions)
    list_combined$revisions <- revisions_overall

  for(name in c("batches", "batch_parameters", "reporting_effects",
               "report_effects", "Rt_changepoints", "reporting_days",
               "reporting_changepoints", "deaths", "deaths_parameters")) {
    if(name %in% names(list_of_results[[1]]))
      list_combined[[name]] <- purrr::map_dfr(
        seq_along(list_of_results),
        ~dplyr::mutate(list_of_results[[.]][[name]], chain=.))
  }

  list_combined
}

#' Runs MCMC or optimisation to estimate Rt, cases and reporting parameters
#'
#' @param niterations number of MCMC iterations to run or number of iterative maximisations to run
#' @param data a tibble with five columns: time_onset, time_reported, cases_reported, Rt_index, reporting_index
#' @param priors a named list with: 'Rt', 'reporting', 'max_cases'. These take
#' the form: 'Rt' is a named list with elements 'shape' and 'rate' describing
#' the gamma prior for each Rt; 'reporting' is a named list with elements
#' 'mean_mu', 'mean_sigma', 'sd_mu', 'sd_sigma' representing the gamma
#' prior parameters for the mean and sd parameters of the reporting parameters
#' (itself described by a gamma distribution); max cases controls the upper
#' limit of the discrete uniform distribution representing the prior on true
#' cases
#' @inheritParams sample_Rt_single_piece
#' @inheritParams max_uncertain_days
#' @inheritParams mcmc_single
#' @param print_to_screen prints progress of MCMC sampling to screen. Defaults to true. Disabled when is_parallel is TRUE.
#' @param nchains number of Markov chains to run. Defaults to 1
#' @param is_parallel Boolean to indicate whether or not to run chains in parallel. Defaults to FALSE.
#' @return a named list of three tibbles: "cases", "Rt" and "reporting" which contain estimates of the model parameters
#' @export
#' @importFrom rlang .data
mcmc <- function(
    niterations,
    data,
    priors,
    serial_parameters,
    initial_cases_true,
    initial_reporting_parameters,
    initial_Rt,
    reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
    serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE, print_to_screen=TRUE,
    nchains=1, is_parallel=FALSE,
    initial_overdispersion=5,
    is_negative_binomial=FALSE,
    overdispersion_metropolis_sd=0.25,
    condition_on_future_cases=TRUE,
    initial_revision_parameters=NULL,
    revision_metropolis_parameters=list(eta_step=0.5, rho_step=0.5),
    cases_before=NULL,
    n_parameter_steps=1) {

  if(nchains==1) {
    res <- mcmc_single(niterations,
                       data,
                       priors,
                       serial_parameters,
                       initial_cases_true,
                       initial_reporting_parameters,
                       initial_Rt,
                       reporting_metropolis_parameters,
                       serial_max,
                       p_gamma_cutoff,
                       maximise,
                       print_to_screen,
                       initial_overdispersion,
                       is_negative_binomial,
                       overdispersion_metropolis_sd,
                       condition_on_future_cases,
                       initial_revision_parameters,
                       revision_metropolis_parameters,
                       cases_before=cases_before,
                       n_parameter_steps=n_parameter_steps)
    res$cases$chain <- 1
    res$Rt$chain <- 1
    res$reporting$chain <- 1
    if(is_negative_binomial)
      res$overdispersion$chain <- 1
    if(!is.null(initial_revision_parameters))
      res$revisions$chain <- 1
    for(name in c("batches", "batch_parameters", "reporting_effects",
                  "report_effects", "Rt_changepoints", "reporting_days",
                  "reporting_changepoints"))
      if(name %in% names(res))
        res[[name]] <- dplyr::mutate(res[[name]], chain=1)
  } else {

    list_of_results <- vector(mode = "list", length = nchains)
    if(is_parallel) {

      if (requireNamespace("foreach", quietly = TRUE)) {

        f_run_single <- function() {
          mcmc_single(niterations,
                      data,
                      priors,
                      serial_parameters,
                      initial_cases_true,
                      initial_reporting_parameters,
                      initial_Rt,
                      reporting_metropolis_parameters,
                      serial_max,
                      p_gamma_cutoff,
                      maximise,
                      print_to_screen,
                      initial_overdispersion,
                      is_negative_binomial,
                      overdispersion_metropolis_sd,
                      condition_on_future_cases,
                      initial_revision_parameters,
                      revision_metropolis_parameters,
                      cases_before=cases_before,
                      n_parameter_steps=n_parameter_steps)
        }

        list_of_results <- foreach::foreach(i=1:nchains, .export = "mcmc_single") %dopar% {
          res <- f_run_single()
        }

      } else {
        # nocov start
        warning("The foreach package must be installed to use this functionality")
        #Either exit or do something without rgl
        return(NULL)
        # nocov end
      }

    } else {
      for(i in seq_along(1:nchains)) {
        res <- mcmc_single(niterations,
                           data,
                           priors,
                           serial_parameters,
                           initial_cases_true,
                           initial_reporting_parameters,
                           initial_Rt,
                           reporting_metropolis_parameters,
                           serial_max,
                           p_gamma_cutoff,
                           maximise,
                           print_to_screen,
                           initial_overdispersion,
                           is_negative_binomial,
                           overdispersion_metropolis_sd,
                           condition_on_future_cases,
                           initial_revision_parameters,
                           revision_metropolis_parameters,
                           cases_before=cases_before,
                       n_parameter_steps=n_parameter_steps)
        list_of_results[[i]] <- res
      }
    }

    res <- combine_chains(list_of_results, is_negative_binomial)
  }

  res
}
