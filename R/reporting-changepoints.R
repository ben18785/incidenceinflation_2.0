#' Creates a function giving the observation process log-probability of sets
#' of onset times which share a reporting delay distribution
#'
#' @param snapshot_with_true_cases_df a tibble with columns: time_onset,
#' time_reported, cases_reported, cases_true
#' @inheritParams observation_process_all_times_logp
#'
#' @return a function of 'indices' (which of the onset times, in increasing
#' order, to include) and 'reporting_parameters' (a named list of 'mean' and
#' 'sd') which returns a log-probability
#' @importFrom rlang .data
reporting_logp_function <- function(snapshot_with_true_cases_df,
                                    revision_parameters=NULL,
                                    max_erroneous=3) {
  df <- snapshot_with_true_cases_df %>%
    dplyr::arrange(.data$time_onset, .data$time_reported)
  by_onset_df <- df %>%
    dplyr::select("time_onset", "cases_true") %>%
    unique()
  onset_times <- by_onset_df$time_onset
  if(length(onset_times) != dplyr::n_distinct(onset_times))
    stop("There must be only one true case measurement per each onset time.")
  l_observations <- split(
    as.data.frame(df[, c("time_onset", "time_reported", "cases_reported")]),
    df$time_onset)
  cases_by_onset <- stats::setNames(by_onset_df$cases_true, onset_times)
  is_misdating <- !is.null(revision_parameters) &&
    misdating_prob(revision_parameters) > 0

  function(indices, reporting_parameters) {
    if(reporting_parameters$mean <= 0 | reporting_parameters$sd <= 0)
      return(-Inf)
    sum(purrr::map_dbl(indices, function(k) {
      onset_time <- onset_times[k]
      observation_df <- l_observations[[as.character(onset_time)]]
      if(nrow(observation_df) <= 1)
        return(0)
      if(is_misdating)
        observation_process_revisions_logp(
          observation_df, cases_by_onset[[k]], onset_time,
          reporting_parameters, revision_parameters, max_erroneous,
          misdated_entries(observation_df$time_reported, onset_time,
                           cases_by_onset, reporting_parameters,
                           revision_parameters))
      else
        observation_process_logp(
          observation_df, cases_by_onset[[k]], onset_time,
          reporting_parameters, revision_parameters, max_erroneous)
    }))
  }
}

#' Log-prior density of the reporting delay parameters of segments
#'
#' The mean and the sd of the delay in each segment are log-normally
#' distributed about a centre shared by all segments, with a shared spread
#' (the sd of their logs).
#'
#' @param parameters a tibble with a row for each segment and columns 'mean'
#' and 'sd'
#' @param centre a named list of 'mean' and 'sd': the values about which the
#' parameters of segments are distributed
#' @param spread the sd of the logs of the parameters of segments about the
#' logs of the centre
#'
#' @return a vector with a log-density for each segment
reporting_segments_logp_prior <- function(parameters, centre, spread) {
  stats::dlnorm(parameters$mean, log(centre$mean), spread, log=TRUE) +
    stats::dlnorm(parameters$sd, log(centre$sd), spread, log=TRUE)
}

#' Log-prior density of the centre and spread of the reporting delay
#' parameters of segments
#'
#' The centre has the gamma priors otherwise used for the reporting
#' parameters themselves, and the spread has a half-normal prior.
#'
#' @inheritParams reporting_segments_logp_prior
#' @param prior_parameters named list with elements 'mean_mu', 'mean_sigma',
#' 'sd_mu', 'sd_sigma'
#' @param spread_sd the sd of the half-normal prior on the spread
#'
#' @return a log-density
reporting_centre_logp_prior <- function(centre, spread, prior_parameters,
                                        spread_sd) {
  if(centre$mean <= 0 | centre$sd <= 0 | spread <= 0)
    return(-Inf)
  dgamma_mean_sd(centre$mean, prior_parameters$mean_mu,
                 prior_parameters$mean_sigma, log=TRUE) +
    dgamma_mean_sd(centre$sd, prior_parameters$sd_mu,
                   prior_parameters$sd_sigma, log=TRUE) +
    stats::dnorm(spread, 0, spread_sd, log=TRUE)
}

#' Number of ways of splitting a segment in two
#'
#' @param length number of days in the segment
#' @inheritParams is_valid_partition
#'
#' @return a count
n_splits <- function(length, min_length=1) {
  max(length - 2 * min_length + 1, 0)
}

#' An initial state for sampling reporting delays whose changes are unknown
#'
#' @param reporting_parameters a named list of 'mean' and 'sd' for a single
#' segment covering all onset times; or a tibble with columns 'mean' and 'sd'
#' and a row for each segment
#' @param starts a vector of the onset times (counted from one, in increasing
#' order) on which segments start
#' @param spread initial value of the spread of the parameters of segments
#'
#' @return a named list of 'starts', 'parameters' (a tibble with columns
#' 'mean' and 'sd'), 'centre' (a named list of 'mean' and 'sd') and 'spread'
#' @export
initial_reporting_changepoints <- function(reporting_parameters, starts=1,
                                           spread=0.3) {
  parameters <- dplyr::tibble(mean=reporting_parameters$mean,
                              sd=reporting_parameters$sd)
  if(nrow(parameters) != length(starts))
    stop("There must be reporting parameters for each segment.")
  list(
    starts=starts,
    parameters=parameters,
    centre=list(mean=exp(mean(log(parameters$mean))),
                sd=exp(mean(log(parameters$sd)))),
    spread=spread)
}

#' Default settings for sampling reporting delays whose changes are unknown
#'
#' @param n_days number of onset times
#' @param settings a named list of any settings to change
#'
#' @return a named list of: 'alpha' and 'beta' (parameters of the beta prior
#' on the probability of a change on each day, see
#' \code{changepoints_logp_prior}); 'min_length' (the minimum number of onset
#' times in a segment); 'spread_sd' (see \code{reporting_centre_logp_prior});
#' 'mean_step' and 'sd_step' (step sizes for the parameters of a segment);
#' 'birth_sd' (the sd of the log of the ratio of the parameters proposed for a
#' new segment to those of the segment it is split from); 'shift_max' (the
#' furthest a change is proposed to move); 'hyper_step' and 'n_hyper_steps'
#' (the step size and number of updates of the centre and spread);
#' 'n_partition_moves' (the number of times in each sweep that a split or
#' merge, and then a shift, is proposed); 'n_segment_steps' (the number of
#' times in each sweep that the parameters of each segment are updated)
reporting_changepoint_settings <- function(n_days, settings=list()) {
  defaults <- list(alpha=1, beta=n_days, min_length=7, spread_sd=0.5,
                   mean_step=0.25, sd_step=0.1, birth_sd=0.3, shift_max=3,
                   hyper_step=0.3, n_hyper_steps=5, n_partition_moves=5,
                   n_segment_steps=1)
  unknown <- setdiff(names(settings), names(defaults))
  if(length(unknown) > 0)
    stop(paste0("Unknown settings: ", paste(unknown, collapse=", "), "."))
  defaults[names(settings)] <- settings
  if(defaults$alpha <= 0 | defaults$beta <= 0)
    stop("The parameters of the prior on the probability of a change must be positive.")
  if(defaults$min_length < 1)
    stop("The minimum length of a segment must be at least one day.")
  defaults
}

#' Updates the reporting delay parameters of each segment in turn
#'
#' @param state a named list as returned by
#' \code{initial_reporting_changepoints}, with an additional element 'logl': a
#' vector of the observation process log-probability of each segment
#' @param f_logp a function as returned by \code{reporting_logp_function}
#' @param n_days number of onset times
#' @param settings a named list as returned by
#' \code{reporting_changepoint_settings}
#'
#' @return the updated state
reporting_segments_step <- function(state, f_logp, n_days, settings) {
  ends <- c(state$starts[-1] - 1, n_days)
  for(j in seq_along(state$starts)) {
    current <- state$parameters[j, ]
    proposed <- dplyr::tibble(
      mean=stats::rnorm(1, current$mean, settings$mean_step),
      sd=stats::rnorm(1, current$sd, settings$sd_step))
    if(proposed$mean <= 0 | proposed$sd <= 0)
      next
    logl_proposed <- f_logp(state$starts[j]:ends[j], as.list(proposed))
    log_r <- logl_proposed - state$logl[j] +
      reporting_segments_logp_prior(proposed, state$centre, state$spread) -
      reporting_segments_logp_prior(current, state$centre, state$spread)
    if(!is.na(log_r) && log(stats::runif(1)) < log_r) {
      state$parameters[j, ] <- proposed
      state$logl[j] <- logl_proposed
    }
  }
  state
}

#' Proposes splitting a segment of constant reporting delay in two
#'
#' A segment and a day within it are chosen uniformly. The days on one side
#' (chosen with equal probability) keep the parameters of the segment and
#' those on the other are given new parameters, proposed log-normally about
#' them.
#'
#' @inheritParams reporting_segments_step
#'
#' @return the updated state
reporting_birth_step <- function(state, f_logp, n_days, settings) {
  starts <- state$starts
  ends <- c(starts[-1] - 1, n_days)
  n_segments <- length(starts)
  j <- sample.int(n_segments, 1)
  n_ways <- n_splits(ends[j] - starts[j] + 1, settings$min_length)
  if(n_ways < 1)
    return(state)
  start_new <- starts[j] + settings$min_length - 1 + sample.int(n_ways, 1)
  current <- state$parameters[j, ]
  proposed <- dplyr::tibble(
    mean=stats::rlnorm(1, log(current$mean), settings$birth_sd),
    sd=stats::rlnorm(1, log(current$sd), settings$birth_sd))
  log_q <- stats::dlnorm(proposed$mean, log(current$mean), settings$birth_sd, log=TRUE) +
    stats::dlnorm(proposed$sd, log(current$sd), settings$birth_sd, log=TRUE)
  is_new_left <- stats::runif(1) < 0.5
  if(is_new_left)
    days_new <- starts[j]:(start_new - 1)
  else
    days_new <- start_new:ends[j]
  logl_new <- f_logp(days_new, as.list(proposed))
  logl_old <- f_logp(days_new, as.list(current))
  n_changes <- n_segments - 1
  log_r <- logl_new - logl_old +
    changepoints_logp_prior(n_changes + 1, n_days, settings$alpha, settings$beta) -
    changepoints_logp_prior(n_changes, n_days, settings$alpha, settings$beta) +
    reporting_segments_logp_prior(proposed, state$centre, state$spread) -
    log_q + log(n_ways)
  if(is.na(log_r) || log(stats::runif(1)) >= log_r)
    return(state)

  logl_kept <- state$logl[j] - logl_old
  if(is_new_left) {
    parameters_pair <- dplyr::bind_rows(proposed, current)
    logl_pair <- c(logl_new, logl_kept)
  } else {
    parameters_pair <- dplyr::bind_rows(current, proposed)
    logl_pair <- c(logl_kept, logl_new)
  }
  state$starts <- append(starts, start_new, after=j)
  state$parameters <- dplyr::bind_rows(
    state$parameters[seq_len(j - 1), ],
    parameters_pair,
    state$parameters[seq_len(n_segments - j) + j, ])
  state$logl <- c(state$logl[seq_len(j - 1)], logl_pair,
                  state$logl[seq_len(n_segments - j) + j])
  state
}

#' Proposes merging two neighbouring segments of constant reporting delay
#'
#' A change is chosen uniformly, and the segments either side of it are
#' merged, taking the parameters of one of them (chosen with equal
#' probability). This is the reverse of \code{reporting_birth_step}.
#'
#' @inheritParams reporting_segments_step
#'
#' @return the updated state
reporting_death_step <- function(state, f_logp, n_days, settings) {
  starts <- state$starts
  ends <- c(starts[-1] - 1, n_days)
  n_segments <- length(starts)
  n_changes <- n_segments - 1
  if(n_changes < 1)
    return(state)
  i <- sample.int(n_changes, 1)
  is_left_kept <- stats::runif(1) < 0.5
  j_kept <- ifelse(is_left_kept, i, i + 1)
  j_dropped <- ifelse(is_left_kept, i + 1, i)
  kept <- state$parameters[j_kept, ]
  dropped <- state$parameters[j_dropped, ]
  logl_merged <- f_logp(starts[j_dropped]:ends[j_dropped], as.list(kept))
  log_q <- stats::dlnorm(dropped$mean, log(kept$mean), settings$birth_sd, log=TRUE) +
    stats::dlnorm(dropped$sd, log(kept$sd), settings$birth_sd, log=TRUE)
  n_ways <- n_splits(ends[i + 1] - starts[i] + 1, settings$min_length)
  log_r <- logl_merged - state$logl[j_dropped] +
    changepoints_logp_prior(n_changes - 1, n_days, settings$alpha, settings$beta) -
    changepoints_logp_prior(n_changes, n_days, settings$alpha, settings$beta) -
    reporting_segments_logp_prior(dropped, state$centre, state$spread) +
    log_q - log(n_ways)
  if(is.na(log_r) || log(stats::runif(1)) >= log_r)
    return(state)

  state$logl[j_kept] <- state$logl[j_kept] + logl_merged
  state$starts <- starts[-(i + 1)]
  state$parameters <- state$parameters[-j_dropped, ]
  state$logl <- state$logl[-j_dropped]
  state
}

#' Proposes moving a change in the reporting delay by a few days
#'
#' @inheritParams reporting_segments_step
#'
#' @return the updated state
reporting_shift_step <- function(state, f_logp, n_days, settings) {
  starts <- state$starts
  n_changes <- length(starts) - 1
  if(n_changes < 1)
    return(state)
  i <- sample.int(n_changes, 1)
  shift <- sample(c(-1, 1), 1) * sample.int(settings$shift_max, 1)
  starts_proposed <- starts
  starts_proposed[i + 1] <- starts[i + 1] + shift
  if(!is_valid_partition(starts_proposed, n_days, settings$min_length))
    return(state)
  if(shift > 0) {
    # days pass from the later segment to the earlier one
    days <- starts[i + 1]:(starts_proposed[i + 1] - 1)
    j_from <- i + 1
    j_to <- i
  } else {
    days <- starts_proposed[i + 1]:(starts[i + 1] - 1)
    j_from <- i
    j_to <- i + 1
  }
  logl_to <- f_logp(days, as.list(state$parameters[j_to, ]))
  logl_from <- f_logp(days, as.list(state$parameters[j_from, ]))
  log_r <- logl_to - logl_from
  if(is.na(log_r) || log(stats::runif(1)) >= log_r)
    return(state)
  state$starts <- starts_proposed
  state$logl[j_to] <- state$logl[j_to] + logl_to
  state$logl[j_from] <- state$logl[j_from] - logl_from
  state
}

#' Updates the centre and spread of the reporting delay parameters of segments
#'
#' Random walk Metropolis steps are taken on the logs of the centre of the
#' means, the centre of the sds and the spread.
#'
#' @inheritParams reporting_segments_step
#' @inheritParams reporting_centre_logp_prior
#'
#' @return the updated state
reporting_centre_step <- function(state, prior_parameters, settings) {
  f_target <- function(centre, spread) {
    logp <- reporting_centre_logp_prior(centre, spread, prior_parameters,
                                        settings$spread_sd)
    if(!is.finite(logp))
      return(-Inf)
    # the last terms are the Jacobian of taking logs
    logp + sum(reporting_segments_logp_prior(state$parameters, centre, spread)) +
      log(centre$mean) + log(centre$sd) + log(spread)
  }
  logp_current <- f_target(state$centre, state$spread)
  for(i in seq_len(settings$n_hyper_steps)) {
    centre <- list(
      mean=state$centre$mean * exp(stats::rnorm(1, 0, settings$hyper_step)),
      sd=state$centre$sd * exp(stats::rnorm(1, 0, settings$hyper_step)))
    spread <- state$spread * exp(stats::rnorm(1, 0, settings$hyper_step))
    logp_proposed <- f_target(centre, spread)
    log_r <- logp_proposed - logp_current
    if(!is.na(log_r) && log(stats::runif(1)) < log_r) {
      state$centre <- centre
      state$spread <- spread
      logp_current <- logp_proposed
    }
  }
  state
}

#' Samples reporting delays which change at unknown onset times
#'
#' Onset times are partitioned into segments, within each of which the mean
#' and sd of the reporting delay are constant. One sweep of
#' Metropolis-Hastings moves is made: the parameters of each segment are
#' updated; several times over, a segment is split in two or two are merged
#' (with equal probability) and a change is moved; and the centre and spread
#' of the distribution of the parameters of segments are updated.
#'
#' The prior probability of a partition is described in
#' \code{changepoints_logp_prior}, with the probability of a change on each
#' day integrated out.
#'
#' @param state a named list as returned by
#' \code{initial_reporting_changepoints}
#' @inheritParams reporting_segments_step
#' @inheritParams reporting_centre_logp_prior
#'
#' @return the updated state, with an element 'logl' giving the observation
#' process log-probability of each segment
#' @export
sample_reporting_changepoints <- function(state, f_logp, n_days,
                                          prior_parameters, settings) {
  if(!is_valid_partition(state$starts, n_days, settings$min_length))
    stop("The segments are not a valid partition of the onset times.")
  ends <- c(state$starts[-1] - 1, n_days)
  # recalculated on each call since the true cases will have been updated
  state$logl <- purrr::map_dbl(seq_along(state$starts), function(j) {
    f_logp(state$starts[j]:ends[j], as.list(state$parameters[j, ]))
  })

  for(i in seq_len(settings$n_segment_steps))
    state <- reporting_segments_step(state, f_logp, n_days, settings)
  for(i in seq_len(settings$n_partition_moves)) {
    if(stats::runif(1) < 0.5)
      state <- reporting_birth_step(state, f_logp, n_days, settings)
    else
      state <- reporting_death_step(state, f_logp, n_days, settings)
    state <- reporting_shift_step(state, f_logp, n_days, settings)
  }
  reporting_centre_step(state, prior_parameters, settings)
}

#' Reporting delay parameters for each onset time given those of segments
#'
#' @inheritParams piecewise_values
#' @param means a vector of the mean reporting delay in each segment
#' @param sds a vector of the sd of the reporting delay in each segment
#'
#' @return a tibble with a row for each onset time and columns: 'time_onset',
#' 'mean' and 'sd'
#' @export
piecewise_reporting_parameters <- function(starts, means, sds, n_days) {
  dplyr::tibble(
    time_onset=seq_len(n_days),
    mean=piecewise_values(starts, means, n_days),
    sd=piecewise_values(starts, sds, n_days))
}
