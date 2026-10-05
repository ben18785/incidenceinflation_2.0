#' Probability case remains undetected over time
#'
#' This is effectively a survival function, which shows the
#' probability that an actual case originating on day t_onset
#' remains unreported at t_1.
#'
#' @param t_1 observation time subsequent to time when case originates: t1>=t0
#' @param t_onset time when case originates
#' @param reporting_parameters named list of 'mean' and 'sd' of gamma distribution
#' characterising the reporting delay distribution
#'
#' @return a probability
#' @examples
#' incidenceinflation:::undetected_prob(2, 1, list(mean=10, sd=5))
undetected_prob <- function(t_1, t_onset, reporting_parameters){
  if(t_1 < t_onset)
    stop("t_1 must equal or exceed t_onset")
  delay_mean <- reporting_parameters$mean
  delay_sd <- reporting_parameters$sd
  1 - pgamma_mean_sd(t_1 - t_onset, delay_mean, delay_sd)
}

#' Probability a case is detected by a later time given it was unobserved
#' earlier
#'
#' This is the probability that a case originating on day day_onset is detected
#' between days day_1 and day_2 given that it was unobserved on day day_1. This
#' probability is calculated using the survival functions (which indicate the
#' probability a case remains undetected by a certain time):
#' \deqn{prob(detect) = \frac{S(day_1|day_onset) - S(day_2|day_onset)}{S(day_1|day_onset)}}
#'
#' @param day_2 day of last observation: day_2 > day_1 >= day_onset
#' @param day_1 day of first observation: day_1 >= day_onset
#' @param day_onset day when case originates
#' @inheritParams undetected_prob
#'
#' @return a probability
detected_after_unobserved_prob <- function(day_2, day_1, day_onset,
                                           reporting_parameters){
  if(day_2 < day_1)
    stop("second observation date must be at same time or after first")
  (undetected_prob(day_1 + 0.5, day_onset, reporting_parameters) -
      undetected_prob(day_2 + 0.5, day_onset, reporting_parameters)) /
    undetected_prob(day_1 + 0.5, day_onset, reporting_parameters)
}

#' Cases arising on a given day which are reported between two subsequent days
#'
#' The number of cases originating on day_onset which are subsequently reported
#' between day_1 and day_2 is modelled as a binomial distribution:
#' \deqn{cases_reported\sim binomial(cases_true-cases_observed, p_detect)}
#' where cases_observed is the number of cases originating on day_onset reported on
#' day_1; p_detect is the probability a case was undetected on day_1 but
#' is detected by day_2.
#'
#' @inheritParams detected_after_unobserved_prob
#' @param cases_true true case count originating on day_onset
#' @param cases_observed observed case count for cases originating on day_onset on
#' day_1
#'
#' @return a count representing number of cases reported between day_1 and day_2
observed_cases_single <- function(cases_observed, cases_true, day_2, day_1, day_onset,
                                  reporting_parameters){
  if(cases_true < cases_observed)
    stop("true case count must exceed reported")
  I_remaining <- cases_true - cases_observed
  if(I_remaining == 0)
    cases <- 0
  else
    cases <- stats::rbinom(1, I_remaining,
                           detected_after_unobserved_prob(day_2, day_1, day_onset,
                                                          reporting_parameters))
  cases
}


#' Generate trajectory of reported cases for a given count originating on a
#' particular day
#'
#' If revision_parameters is supplied, the reported case counts also include
#' erroneous records, which are entered and later withdrawn, so that counts can
#' be revised downwards (see \code{observation_process_revisions_logp}). No
#' erroneous records are entered before the first day on which cases are counted.
#'
#' @inheritParams observed_cases_single
#' @inheritParams observation_process_revisions_logp
#' @param days_reporting a vector of days on which cases were counted
#'
#' @return a vector of reported case counts of same length as days_reporting
observed_cases_trajectory <- function(cases_true, days_reporting, day_onset,
                                      reporting_parameters,
                                      revision_parameters=NULL){

  I_observed <- vector(length = length(days_reporting))
  # genuine cases reported and erroneous records not yet withdrawn
  I_genuine <- 0
  I_erroneous <- 0
  for(i in 1:length(days_reporting)){
    if(i == 1) {
      day_previous_report <- day_onset
    } else {
      day_previous_report <- days_reporting[i - 1]
    }
    new_cases <- observed_cases_single(I_genuine, cases_true,
                                       days_reporting[i],
                                       day_previous_report,
                                       day_onset,
                                       reporting_parameters)
    if(!is.null(revision_parameters) & i > 1) {
      p_newly_detected <-
        undetected_prob(day_previous_report + 0.5, day_onset, reporting_parameters) -
        undetected_prob(days_reporting[i] + 0.5, day_onset, reporting_parameters)
      p_withdrawn <- 1 - (1 - revision_parameters$rho)^(days_reporting[i] - day_previous_report)
      withdrawn <- stats::rbinom(1, I_erroneous, p_withdrawn)
      entered <- stats::rpois(1, revision_parameters$eta * cases_true * p_newly_detected)
      I_erroneous <- I_erroneous + entered - withdrawn
    }
    I_genuine <- I_genuine + new_cases
    I_observed[i] <- I_genuine + I_erroneous
  }
  I_observed
}

#' Generates true cases for a single day
#'
#' Cases are assumed to be generated by a negative-binomial renewal
#' process:
#' \deqn{cases_t \sim neg-binomial(Rt \sum_s=1^\infty w(s) cases_t-s, kappa)}
#' where kappa is the over-dispersion parameter. Here, we truncate this
#' sum at a user specified lag (determined by the length of the weight vector).
#'
#' @param Rt effective reproduction number
#' @param kappa over-dispersion parameter of negative binomial
#' @param cases_history a vector containing history of cases arranged from
#' recent to past
#' @param weights a vector containing weights representing the serial interval
#' distribution
#'
#' @return a count of cases
true_cases_single <- function(Rt, kappa, cases_history, weights){
  if(length(cases_history) != length(weights))
    stop("case history vector must be same length as weights")
  mean_cases <- expected_cases(Rt=Rt, weights=weights, cases_history=cases_history)
  stats::rnbinom(1, mu=mean_cases, size=kappa)
}

#' Generate true case series
#'
#' Cases are assumed to be generated by a negative-binomial renewal
#' process:
#' \deqn{cases_t \sim neg-binomial(Rt \sum_s=1^\infty w(s) cases_t-s, kappa)}
#' where kappa is the over-dispersion parameter. Here, we truncate this
#' sum at a user specified lag (determined by the length of the weight vector).
#'
#' @param days_total number of days to simulate process for (post-initial
#' period)
#' @param Rt_function takes day as an input and outputs an Rt value
#' @inheritParams true_cases_single
#' @inheritParams gamma_discrete_pmf
#' @param initial_parameters a named list of 'mean' and 'length'
#' which is used to generate seed cases by sampling from a negative binomial
#' distribution
#' @param serial_max maximum point at which to truncate sum in renewal process
#'
#' @return vector of true cases of length days_total
#' @export
#'
#' @examples
#' library(incidenceinflation)
#' # generate case series for 100 days with time-varying Rt
#' days_total <- 100
#'
#' # make Rt interpolation function
#' v_Rt <- c(rep(1.3, 25), rep(1, 25), rep(2, 50))
#' Rt_function <- stats::approxfun(1:days_total, v_Rt)
#'
#' # serial interval parameters
#' s_params <- list(mean=5, sd=1)
#'
#' # renewal parameter
#' kappa <- 2
#' cases <- true_cases(days_total, Rt_function, kappa, s_params)
true_cases <- function(days_total, Rt_function, kappa, serial_parameters,
                       initial_parameters=list(mean=3, length=5),
                       serial_max=40){

  w <- weights_series(serial_max, serial_parameters)
  initial_length <- initial_parameters$length
  I_true <- vector(length=(days_total + initial_length))

  # seed epidemic with initial values
  I_true[1:initial_length] <- stats::rnbinom(n=initial_length,
                                             mu=initial_parameters$mean,
                                             size=kappa)

  # generate further incidences using state equation
  for(d in 1:days_total){
    d1 <- d + initial_length
    if(d1 > serial_max) {
      I_history <- I_true[(d1 - serial_max):(d1 - 1)]
    } else {
      I_history <- I_true[1:d1]
      d_remaining <- serial_max - d1
      I_history <- c(rep(0, d_remaining), I_history)
    }
    I_history_rev <- rev(I_history)
    Rt <- Rt_function(d)
    I_true[d1] <- true_cases_single(Rt, kappa, I_history_rev, w)
  }
  I_true[-(1:initial_length)]
}

#' Creates a reporting parameter tibble with a single set of reporting
#' parameters
#'
#' @param time_onsets a vector of symptom onset times
#' @inheritParams undetected_prob
#'
#' @return a tibble with three columns: 'time_onset', 'mean', 'sd'
create_reporting_from_single_parameters_df <- function(time_onsets,
    reporting_parameters) {

  dplyr::tibble(time_onset=time_onsets,
                mean=reporting_parameters$mean,
                sd=reporting_parameters$sd
                )
}

#' Generate reported case trajectories for each day when cases appear
#'
#' @param cases_true a vector of true cases originating each day
#' @param reporting_parameters either a named list of 'mean' and 'sd' of gamma
#' distribution indicating a single reporting delay distribution for the whole period
#' or a tibble with columns: 'time_onset', 'mean', 'sd' indicating the reporting
#' delay distribution for each time period
#' @param days_max_follow_up max days at which to simulate reporting case
#' trajectory
#' @inheritParams observed_cases_trajectory
#'
#' @return a tibble with observed case trajectories for each time of onset
#' @export
#' @importFrom magrittr "%>%"
#' @examples
#' library(incidenceinflation)
#' observed_cases(stats::rpois(5, 5), list(mean=5, sd=1))
observed_cases <- function(cases_true, reporting_parameters,
                           days_max_follow_up=30,
                           revision_parameters=NULL){

  if(methods::is(reporting_parameters, "list"))
    reporting_parameters <- create_reporting_from_single_parameters_df(
      seq_along(cases_true),
      reporting_parameters)

  if(nrow(reporting_parameters) != length(cases_true))
    stop("Number of rows in reporting parameters tibble must match number of cases.")

  d_max <- length(cases_true)
  for(t in 1:d_max){

    a_max <- min(c(d_max, t + days_max_follow_up))
    obs_time <- seq(t, a_max, 1)
    reporting_parameters_current <- list(mean=reporting_parameters$mean[t],
                                         sd=reporting_parameters$sd[t])
    cases_obs_trajec <- observed_cases_trajectory(
      cases_true=cases_true[t],
      days_reporting=obs_time,
      day_onset=t,
      reporting_parameters_current,
      revision_parameters=revision_parameters)

    # stack into data frame containing identifying info
    I_obs_single_onset <- dplyr::tibble(
      time_onset=rep(t, length(cases_obs_trajec)),
      time_reported=obs_time,
      cases_reported=cases_obs_trajec,
      cases_true=cases_true[t])
    if(t == 1)
      cases_obs <- I_obs_single_onset
    else
      cases_obs <- cases_obs %>%
      dplyr::bind_rows(I_obs_single_onset)
  }
  cases_obs
}

#' Thins cases for a single onset time
#'
#' Only those data points around where the case count changes contain
#' information. This method searches and keeps only those data points.
#'
#' @param case_obs_single a tibble produced by observed_cases but subsetted to
#' a single onset time
#'
#' @return a thinned tibble
thin_single_series <- function(case_obs_single) {
  if(dplyr::n_distinct(case_obs_single$time_onset) > 1)
    stop("Only single onset time allowed.")
  changepoints <- which(diff(case_obs_single$cases_reported) != 0)
  postchangepoints <- changepoints + 1
  num_points <- nrow(case_obs_single)
  interesting_points <- sort(unique(c(1, postchangepoints, changepoints,
                                      num_points)))
  case_obs_single[interesting_points, ]
}

#' Thins cases
#'
#' Only those data points around where the case count changes contain
#' information. This method searches and keeps only those data points.
#'
#' @param case_obs a tibble produced by observed_cases
#'
#' @return a thinned tibble
#' @importFrom magrittr "%>%"
#' @export
thin_series <- function(case_obs) {
  onset_times <- unique(case_obs$time_onset)
  for(i in seq_along(onset_times)) {
    single_df <- case_obs %>%
      dplyr::filter(case_obs$time_onset == onset_times[i])
    thinned_single_df <- thin_single_series(single_df)
    if(i == 1)
      thinned_df <- thinned_single_df
    else
      thinned_df <- thinned_df %>%
      dplyr::bind_rows(thinned_single_df)
  }
  thinned_df
}

#' Generates case data with snapshots of reported cases at each onset day
#'
#' @inheritParams true_cases
#' @inheritParams observed_cases
#' @inheritParams observed_cases_revisions
#' @param thinned a Boolean indicating whether to thin cases to only informative time points
#'
#' @return a tibble with observed case trajectories for each time of onset
#' @export
#' @examples
#' library(incidenceinflation)
#' # generate case series for 100 days with time-varying Rt
#' days_total <- 100
#'
#' # make Rt interpolation function
#' v_Rt <- c(rep(1.3, 25), rep(1, 25), rep(2, 50))
#' Rt_function <- stats::approxfun(1:days_total, v_Rt)
#'
#' # serial interval parameters
#' s_params <- list(mean=5, sd=1)
#'
#' # renewal over-dispersion parameter
#' kappa <- 2
#'
#' # reporting delays
#' r_params <- list(mean=10, sd=5)
#' reported_cases <- generate_snapshots(
#'     days_total, Rt_function,
#'     s_params, r_params)
generate_snapshots <- function(days_total, Rt_function, serial_parameters,
                               reporting_parameters,
                               thinned=FALSE,
                               kappa=1000,
                               days_max_follow_up=30,
                               initial_parameters=list(mean=30, length=20),
                               serial_max=40,
                               revision_parameters=NULL,
                               batches=NULL) {

  check_parameter_names(reporting_parameters,
                        serial_parameters)

  real_cases <- true_cases(
    days_total=days_total, Rt_function=Rt_function,
    kappa=kappa, serial_parameters=serial_parameters,
    initial_parameters=initial_parameters,
    serial_max=serial_max
  )

  is_misdating <- !is.null(revision_parameters) &&
    misdating_prob(revision_parameters) > 0
  if(is_misdating | !is.null(batches))
    reported_cases <- observed_cases_revisions(
      cases_true=real_cases,
      reporting_parameters=reporting_parameters,
      revision_parameters=revision_parameters,
      days_max_follow_up=days_max_follow_up,
      batches=batches)
  else
    reported_cases <- observed_cases(
      cases_true=real_cases,
      reporting_parameters=reporting_parameters,
      days_max_follow_up=days_max_follow_up,
      revision_parameters=revision_parameters)

  if(thinned)
    reported_cases <- thin_series(reported_cases)

  reported_cases
}

#' Draws batches of erroneous records from their prior
#'
#' Onset days are considered in order. An onset day lying within the footprint
#' of an earlier batch starts no batch; any other onset day starts a batch with
#' probability beta. A batch has a footprint of width 1 + Poisson(zeta) onset
#' days, a gamma-distributed intensity and is detected a geometrically
#' distributed number of days after its footprint starts.
#'
#' @param days_total number of onset days
#' @param batch_parameters named list of 'beta', the probability a batch starts
#' on an onset day; 'zeta', the mean number of onset days beyond the first in a
#' footprint; 'phi_mu' and 'phi_sigma', the mean and sd of the gamma
#' distribution of intensities; and 'rho_batch', the daily probability a batch
#' is detected
#'
#' @return a tibble with a row for each batch and columns: 'start', the first
#' onset day of its footprint; 'width', the number of onset days in its
#' footprint; 'phi', its intensity; and 'time_detected', the day on which it
#' is detected
#' @export
sample_batches <- function(days_total, batch_parameters) {
  starts <- c()
  widths <- c()
  day <- 1
  while(day <= days_total) {
    if(stats::runif(1) < batch_parameters$beta) {
      width <- 1 + stats::rpois(1, batch_parameters$zeta)
      starts <- c(starts, day)
      widths <- c(widths, width)
      day <- day + width
    } else {
      day <- day + 1
    }
  }
  n_batches <- length(starts)
  dplyr::tibble(
    start=as.numeric(starts),
    width=as.numeric(widths),
    phi=stats::rgamma(n_batches,
                      shape=batch_parameters$phi_mu^2 / batch_parameters$phi_sigma^2,
                      rate=batch_parameters$phi_mu / batch_parameters$phi_sigma^2),
    time_detected=as.numeric(starts) + 1 +
      stats::rgeom(n_batches, batch_parameters$rho_batch))
}

#' Finds the batch (if any) whose footprint contains an onset day
#'
#' @param day_onset an onset day
#' @param batches a tibble as returned by \code{sample_batches} or NULL
#'
#' @return a named list of 'phi' and 'time_detected' or NULL if the onset day
#' is in no footprint
batch_for_onset <- function(day_onset, batches) {
  if(is.null(batches) || nrow(batches) == 0)
    return(NULL)
  index <- which(batches$start <= day_onset &
                   day_onset < batches$start + batches$width)
  if(length(index) == 0)
    return(NULL)
  if(length(index) > 1)
    stop("Footprints of batches must not overlap.")
  list(phi=batches$phi[index], time_detected=batches$time_detected[index])
}

#' Generate reported case trajectories allowing for all types of revision
#'
#' Reported cases for each onset day are simulated day by day. As well as
#' genuine cases being detected, reported counts can include: spurious records,
#' which are entered and later deleted; misdated records, which are genuine
#' cases recorded under a nearby onset day until corrected, when they move to
#' their true onset day; and records belonging to a batch, which accumulate for
#' the onset days in its footprint until the batch is detected, when all are
#' withdrawn. No erroneous record is entered on its onset day itself, and a
#' record cannot be given an onset day on or after the day it is entered.
#'
#' @inheritParams observed_cases
#' @param revision_parameters named list of 'eta', 'rho' and optionally 'pi',
#' 'shifts' and 'shift_probs' (see \code{unrecorded_prob}); or NULL
#' @param batches a tibble of batches as returned by \code{sample_batches}; or NULL
#'
#' @return a tibble with observed case trajectories for each time of onset
#' @export
observed_cases_revisions <- function(cases_true, reporting_parameters,
                                     revision_parameters=NULL,
                                     days_max_follow_up=30,
                                     batches=NULL) {

  if(methods::is(reporting_parameters, "list"))
    reporting_parameters <- create_reporting_from_single_parameters_df(
      seq_along(cases_true),
      reporting_parameters)
  if(nrow(reporting_parameters) != length(cases_true))
    stop("Number of rows in reporting parameters tibble must match number of cases.")

  eta <- 0
  rho <- 0
  p_misdated <- 0
  shifts <- c()
  shift_probs <- c()
  if(!is.null(revision_parameters)) {
    eta <- revision_parameters$eta
    rho <- revision_parameters$rho
    p_misdated <- misdating_prob(revision_parameters)
    if(p_misdated > 0) {
      shifts <- revision_parameters$shifts
      shift_probs <- revision_parameters$shift_probs
    }
  }

  d_max <- length(cases_true)
  n_shifts <- length(shifts)
  I_detected <- rep(0, d_max) # genuine cases detected, wherever recorded
  I_genuine <- rep(0, d_max) # genuine cases recorded under their onset day
  I_spurious <- rep(0, d_max)
  I_batch <- rep(0, d_max)
  # misdated records by true onset day (rows) and shift (columns)
  m_misdated <- matrix(0, nrow=d_max, ncol=max(n_shifts, 1))
  l_batch <- purrr::map(seq_len(d_max), ~batch_for_onset(., batches))
  m_reported <- matrix(NA, nrow=d_max, ncol=d_max)
  diag(m_reported) <- 0

  for(day in seq_len(d_max)[-1]) {
    for(t in 1:(day - 1)) {
      r_params <- list(mean=reporting_parameters$mean[t],
                       sd=reporting_parameters$sd[t])
      p_newly_detected <- undetected_prob(day - 1 + 0.5, t, r_params) -
        undetected_prob(day + 0.5, t, r_params)

      # corrections of misdated records present on previous day
      if(n_shifts > 0) {
        corrected <- stats::rbinom(n_shifts, m_misdated[t, ], rho)
        m_misdated[t, ] <- m_misdated[t, ] - corrected
        I_genuine[t] <- I_genuine[t] + sum(corrected)
      }

      # newly detected genuine cases, some of which are misdated
      new_cases <- observed_cases_single(I_detected[t], cases_true[t],
                                         day, day - 1, t, r_params)
      I_detected[t] <- I_detected[t] + new_cases
      new_misdated <- 0
      if(n_shifts > 0 & new_cases > 0) {
        is_possible <- shifts <= day - t - 1
        if(any(is_possible)) {
          new_misdated <- stats::rbinom(1, new_cases,
                                        p_misdated * sum(shift_probs[is_possible]))
          if(new_misdated > 0)
            m_misdated[t, is_possible] <- m_misdated[t, is_possible] +
              as.vector(stats::rmultinom(1, new_misdated, shift_probs[is_possible]))
        }
      }
      I_genuine[t] <- I_genuine[t] + new_cases - new_misdated

      # spurious records
      withdrawn <- stats::rbinom(1, I_spurious[t], rho)
      entered <- stats::rpois(1, eta * cases_true[t] * p_newly_detected)
      I_spurious[t] <- I_spurious[t] + entered - withdrawn

      # records belonging to a batch
      batch <- l_batch[[t]]
      if(!is.null(batch)) {
        if(day < batch$time_detected)
          I_batch[t] <- I_batch[t] +
            stats::rpois(1, batch$phi * cases_true[t] * p_newly_detected)
        else
          I_batch[t] <- 0
      }
    }

    # misdated records are reported under the onset day they were given
    I_misdated <- rep(0, d_max)
    for(i in seq_len(n_shifts)) {
      days_given <- seq_len(d_max) + shifts[i]
      is_within <- days_given >= 1 & days_given <= d_max
      I_misdated[days_given[is_within]] <- I_misdated[days_given[is_within]] +
        m_misdated[is_within, i]
    }
    onset_days <- 1:(day - 1)
    m_reported[onset_days, day] <- (I_genuine + I_spurious + I_batch +
                                      I_misdated)[onset_days]
  }

  purrr::map_dfr(seq_len(d_max), function(t) {
    obs_time <- seq(t, min(c(d_max, t + days_max_follow_up)), 1)
    dplyr::tibble(
      time_onset=rep(t, length(obs_time)),
      time_reported=obs_time,
      cases_reported=m_reported[t, obs_time],
      cases_true=cases_true[t])
  })
}
