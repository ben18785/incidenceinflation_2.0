#' Onset times lying in the footprint of a batch
#'
#' @param start first onset day of the footprint
#' @param width number of onset days in the footprint
#' @param onset_times vector of onset times in the data
#'
#' @return those onset times which lie in the footprint
footprint_onset_times <- function(start, width, onset_times) {
  onset_times[onset_times >= start & onset_times < start + width]
}

#' Determines whether the footprints of a set of batches overlap
#'
#' @param batches a tibble of batches as returned by \code{sample_batches}
#'
#' @return TRUE if any two footprints share an onset day
is_overlapping_batches <- function(batches) {
  if(nrow(batches) < 2)
    return(FALSE)
  batches <- batches[order(batches$start), ]
  any(batches$start[-1] < (batches$start + batches$width)[-nrow(batches)])
}

#' Log prior probability of a set of batches
#'
#' Considering onset times in order, an onset time in the footprint of an
#' earlier batch starts no batch and any other starts one with probability
#' beta, so that footprints cannot overlap. The prior probability of the
#' footprints is then:
#' \deqn{beta^K (1 - beta)^F \prod_b Poisson(width_b - 1 | zeta)}
#' where K is the number of batches and F is the number of onset times in no
#' footprint. Intensities have independent gamma priors and each batch is
#' detected with probability rho_batch on each day after its footprint starts.
#'
#' @inheritParams is_overlapping_batches
#' @inheritParams sample_batches
#' @inheritParams footprint_onset_times
#'
#' @return a log-probability, which is -Inf if the batches are not valid
batches_logp_prior <- function(batches, batch_parameters, onset_times) {
  n_onset <- length(onset_times)
  if(nrow(batches) == 0)
    return(n_onset * log(1 - batch_parameters$beta))

  is_valid <- all(batches$width >= 1) &
    all(batches$start %in% onset_times) &
    all(batches$time_detected >= batches$start + 1) &
    all(batches$phi > 0) &
    !is_overlapping_batches(batches)
  if(!is_valid)
    return(-Inf)

  n_covered <- sum(purrr::map2_dbl(
    batches$start, batches$width,
    ~length(footprint_onset_times(.x, .y, onset_times))))
  nrow(batches) * log(batch_parameters$beta) +
    (n_onset - n_covered) * log(1 - batch_parameters$beta) +
    sum(stats::dpois(batches$width - 1, batch_parameters$zeta, log=TRUE)) +
    sum(dgamma_mean_sd(batches$phi, batch_parameters$phi_mu,
                       batch_parameters$phi_sigma, log=TRUE)) +
    sum(stats::dgeom(batches$time_detected - batches$start - 1,
                     batch_parameters$rho_batch, log=TRUE))
}

#' Finds candidate batches from falls in reported cases which coincide
#'
#' A batch which is detected produces falls in the reported cases for adjacent
#' onset times at the same reporting time. For each reporting time, this finds
#' runs of onset times whose reported cases fell (allowing gaps of up to
#' max_gap onset times within a run) and retains those runs whose total fall
#' is at least min_total_fall.
#'
#' @inheritParams max_downward_revision
#' @param min_total_fall smallest total fall in reported cases across a run of
#' onset times for it to be a candidate
#' @param max_gap largest number of consecutive onset times without a fall
#' allowed within a run
#'
#' @return a tibble with a row for each candidate and columns 'start', 'width'
#' and 'time_detected' (as in \code{sample_batches}) and 'total_fall'
#' @export
#' @importFrom rlang .data
batch_candidates <- function(snapshot_df, min_total_fall, max_gap=1) {
  falls_df <- snapshot_df %>%
    dplyr::group_by(.data$time_onset) %>%
    dplyr::arrange(.data$time_reported, .by_group=TRUE) %>%
    dplyr::mutate(fall=c(0, -diff(.data$cases_reported))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data$fall > 0) %>%
    dplyr::arrange(.data$time_reported, .data$time_onset)

  candidates <- dplyr::tibble(start=numeric(0), width=numeric(0),
                              time_detected=numeric(0), total_fall=numeric(0))
  if(nrow(falls_df) == 0)
    return(candidates)

  falls_df <- falls_df %>%
    dplyr::group_by(.data$time_reported) %>%
    dplyr::mutate(run=cumsum(c(1, diff(.data$time_onset) > max_gap + 1))) %>%
    dplyr::group_by(.data$time_reported, .data$run) %>%
    dplyr::summarise(start=min(.data$time_onset),
                     width=max(.data$time_onset) - min(.data$time_onset) + 1,
                     total_fall=sum(.data$fall),
                     .groups="drop") %>%
    dplyr::filter(.data$total_fall >= min_total_fall)

  dplyr::tibble(start=as.numeric(falls_df$start),
                width=as.numeric(falls_df$width),
                time_detected=as.numeric(falls_df$time_reported),
                total_fall=as.numeric(falls_df$total_fall))
}

#' Chooses initial batches from a set of candidates
#'
#' Candidates are taken in decreasing order of their total fall in reported
#' cases, omitting any whose footprint overlaps one already chosen.
#'
#' @param candidates a tibble as returned by \code{batch_candidates}
#' @inheritParams sample_batches
#'
#' @return a tibble of batches as returned by \code{sample_batches}, each with
#' intensity phi_mu
initial_batches <- function(candidates, batch_parameters) {
  batches <- dplyr::tibble(start=numeric(0), width=numeric(0),
                           phi=numeric(0), time_detected=numeric(0))
  candidates <- candidates[order(-candidates$total_fall), ]
  for(i in seq_len(nrow(candidates))) {
    proposed <- dplyr::bind_rows(
      batches,
      dplyr::tibble(start=candidates$start[i], width=candidates$width[i],
                    phi=batch_parameters$phi_mu,
                    time_detected=candidates$time_detected[i]))
    if(!is_overlapping_batches(proposed))
      batches <- proposed
  }
  batches[order(batches$start), ]
}

#' Splits snapshots into the lists used when sampling batches
#'
#' @param snapshot_with_true_cases_df a tibble with columns: time_onset,
#' time_reported, cases_reported, cases_true, reporting_piece_index
#' @param reporting_parameters a tibble with columns: 'reporting_piece_index', 'mean', 'sd'
#'
#' @return a named list of 'onset_times'; 'l_observations', 'l_reporting' and
#' 'cases_by_onset' (as in \code{misdating_observation_logp_function})
#' @importFrom rlang .data
observation_lists <- function(snapshot_with_true_cases_df, reporting_parameters) {
  df <- snapshot_with_true_cases_df %>%
    dplyr::arrange(.data$time_onset, .data$time_reported)
  by_onset_df <- df %>%
    dplyr::select("time_onset", "cases_true", "reporting_piece_index") %>%
    unique() %>%
    dplyr::left_join(
      reporting_parameters[, c("reporting_piece_index", "mean", "sd")],
      by="reporting_piece_index")
  list(
    onset_times=by_onset_df$time_onset,
    l_observations=split(
      as.data.frame(df[, c("time_reported", "cases_reported")]), df$time_onset),
    l_reporting=stats::setNames(
      purrr::map(seq_len(nrow(by_onset_df)),
                 ~list(mean=by_onset_df$mean[.], sd=by_onset_df$sd[.])),
      by_onset_df$time_onset),
    cases_by_onset=stats::setNames(by_onset_df$cases_true, by_onset_df$time_onset))
}

#' Observation process log-probability for each of a set of onset times
#'
#' @param onset_times the onset times for which to calculate log-probabilities
#' @inheritParams misdating_observation_logp_function
#' @param revision_parameters named list of revision parameters, which can
#' contain an element 'batches' (see \code{observation_process_revisions_logp})
#'
#' @return a vector of log-probabilities
observation_logp_by_onset <- function(onset_times, l_observations, cases_by_onset,
                                      l_reporting, revision_parameters,
                                      max_erroneous=3) {
  is_misdating <- misdating_prob(revision_parameters) > 0
  purrr::map_dbl(onset_times, function(onset_time) {
    name_onset <- as.character(onset_time)
    observation_df <- l_observations[[name_onset]]
    if(nrow(observation_df) <= 1)
      return(0)
    entries <- NULL
    if(is_misdating)
      entries <- misdated_entries(observation_df$time_reported, onset_time,
                                  cases_by_onset, l_reporting[[name_onset]],
                                  revision_parameters)
    observation_process_revisions_logp(
      observation_df, cases_by_onset[[name_onset]], onset_time,
      l_reporting[[name_onset]], revision_parameters, max_erroneous, entries)
  })
}

#' Updates batches of erroneous records using Metropolis-Hastings moves
#'
#' One move is made which either adds a batch or removes one, followed by
#' moves which change the intensity, the detection time and the footprint of
#' each batch in turn. A batch to be added is drawn either from the prior or
#' (with equal probability, if there are any) from a set of candidates
#' suggested by the data.
#'
#' @param batches a tibble of current batches as returned by \code{sample_batches}
#' @inheritParams observation_logp_by_onset
#' @param onset_times vector of all onset times in the data
#' @inheritParams sample_batches
#' @param candidates a tibble as returned by \code{batch_candidates}
#' @param phi_step step size for proposals of the log of each intensity
#'
#' @return a tibble of batches
#' @export
sample_batches_step <- function(batches, onset_times, l_observations,
                                cases_by_onset, l_reporting,
                                revision_parameters, max_erroneous,
                                batch_parameters, candidates, phi_step=0.3) {

  f_logp <- function(days, batches_tmp) {
    if(length(days) == 0)
      return(0)
    revision_parameters$batches <- batches_tmp
    sum(observation_logp_by_onset(days, l_observations, cases_by_onset,
                                  l_reporting, revision_parameters, max_erroneous))
  }
  f_footprint <- function(batches_tmp, index)
    footprint_onset_times(batches_tmp$start[index], batches_tmp$width[index],
                          onset_times)
  # a ratio of Inf arises if reported cases are impossible under current batches
  f_is_accepted <- function(log_ratio)
    !is.na(log_ratio) && log(stats::runif(1)) < log_ratio

  n_onset <- length(onset_times)
  n_candidates <- nrow(candidates)
  p_candidate <- if(n_candidates > 0) 0.5 else 0
  f_log_proposal <- function(start, width, time_detected) {
    is_candidate <- any(candidates$start == start & candidates$width == width &
                          candidates$time_detected == time_detected)
    p_prior <- stats::dpois(width - 1, batch_parameters$zeta) *
      stats::dgeom(time_detected - start - 1, batch_parameters$rho_batch) / n_onset
    log((1 - p_candidate) * p_prior +
          ifelse(is_candidate, p_candidate / n_candidates, 0))
  }
  # prior probability of adding a batch, relative to not, excluding intensity
  f_log_prior_added <- function(start, width, time_detected, n_covered) {
    log(batch_parameters$beta) - n_covered * log(1 - batch_parameters$beta) +
      stats::dpois(width - 1, batch_parameters$zeta, log=TRUE) +
      stats::dgeom(time_detected - start - 1, batch_parameters$rho_batch, log=TRUE)
  }

  # add or remove a batch
  n_batches <- nrow(batches)
  if(stats::runif(1) < 0.5) {
    if(stats::runif(1) < p_candidate) {
      index <- sample(n_candidates, 1)
      start <- candidates$start[index]
      width <- candidates$width[index]
      time_detected <- candidates$time_detected[index]
    } else {
      start <- onset_times[sample(n_onset, 1)]
      width <- 1 + stats::rpois(1, batch_parameters$zeta)
      time_detected <- start + 1 + stats::rgeom(1, batch_parameters$rho_batch)
    }
    phi <- stats::rgamma(
      1, shape=batch_parameters$phi_mu^2 / batch_parameters$phi_sigma^2,
      rate=batch_parameters$phi_mu / batch_parameters$phi_sigma^2)
    proposed <- dplyr::bind_rows(
      batches, dplyr::tibble(start=start, width=width, phi=phi,
                             time_detected=time_detected))
    if(!is_overlapping_batches(proposed)) {
      days <- footprint_onset_times(start, width, onset_times)
      log_ratio <- f_logp(days, proposed) - f_logp(days, batches) +
        f_log_prior_added(start, width, time_detected, length(days)) -
        log(n_batches + 1) - f_log_proposal(start, width, time_detected)
      if(f_is_accepted(log_ratio))
        batches <- proposed
    }
  } else if(n_batches > 0) {
    index <- sample(n_batches, 1)
    proposed <- batches[-index, ]
    days <- f_footprint(batches, index)
    log_ratio <- f_logp(days, proposed) - f_logp(days, batches) -
      f_log_prior_added(batches$start[index], batches$width[index],
                        batches$time_detected[index], length(days)) +
      log(n_batches) +
      f_log_proposal(batches$start[index], batches$width[index],
                     batches$time_detected[index])
    if(f_is_accepted(log_ratio))
      batches <- proposed
  }

  # change each batch in turn
  for(index in seq_len(nrow(batches))) {
    # intensity
    proposed <- batches
    proposed$phi[index] <- batches$phi[index] * exp(stats::rnorm(1, 0, phi_step))
    days <- f_footprint(batches, index)
    log_ratio <- f_logp(days, proposed) - f_logp(days, batches) +
      batches_logp_prior(proposed, batch_parameters, onset_times) -
      batches_logp_prior(batches, batch_parameters, onset_times) +
      log(proposed$phi[index] / batches$phi[index])
    if(f_is_accepted(log_ratio))
      batches <- proposed

    # detection time
    proposed <- batches
    proposed$time_detected[index] <- batches$time_detected[index] +
      sample(c(-3:-1, 1:3), 1)
    log_ratio <- batches_logp_prior(proposed, batch_parameters, onset_times) -
      batches_logp_prior(batches, batch_parameters, onset_times)
    if(is.finite(log_ratio))
      log_ratio <- log_ratio + f_logp(days, proposed) - f_logp(days, batches)
    if(f_is_accepted(log_ratio))
      batches <- proposed

    # footprint: extend or shrink at either end
    proposed <- batches
    move <- sample(4, 1)
    if(move == 1) {
      proposed$start[index] <- batches$start[index] - 1
      proposed$width[index] <- batches$width[index] + 1
    } else if(move == 2) {
      proposed$start[index] <- batches$start[index] + 1
      proposed$width[index] <- batches$width[index] - 1
    } else if(move == 3) {
      proposed$width[index] <- batches$width[index] + 1
    } else {
      proposed$width[index] <- batches$width[index] - 1
    }
    log_ratio <- batches_logp_prior(proposed, batch_parameters, onset_times) -
      batches_logp_prior(batches, batch_parameters, onset_times)
    if(is.finite(log_ratio)) {
      days <- union(f_footprint(batches, index), f_footprint(proposed, index))
      log_ratio <- log_ratio + f_logp(days, proposed) - f_logp(days, batches)
    }
    if(f_is_accepted(log_ratio))
      batches <- proposed
  }
  batches
}

#' Samples the parameters governing batches conditional on the batches
#'
#' With a beta prior for beta, a gamma prior for zeta and a beta prior for
#' rho_batch, the conditional posteriors given K batches are:
#' \deqn{beta \sim beta(a + K, b + F)}
#' where F is the number of onset times in no footprint;
#' \deqn{zeta \sim gamma(shape + \sum_b (width_b - 1), rate + K)}
#' and
#' \deqn{rho_batch \sim beta(a + K, b + \sum_b (time_detected_b - start_b - 1))}
#' Any of these whose prior is not given is left unchanged.
#'
#' @inheritParams batches_logp_prior
#' @param batch_priors named list which can contain: 'beta_alpha' and
#' 'beta_beta'; 'zeta_shape' and 'zeta_rate'; 'rho_batch_alpha' and 'rho_batch_beta'
#'
#' @return batch_parameters with sampled values
sample_batch_parameters <- function(batches, batch_parameters, batch_priors,
                                    onset_times) {
  n_batches <- nrow(batches)
  n_covered <- sum(purrr::map2_dbl(
    batches$start, batches$width,
    ~length(footprint_onset_times(.x, .y, onset_times))))
  if(!is.null(batch_priors$beta_alpha))
    batch_parameters$beta <- stats::rbeta(
      1, batch_priors$beta_alpha + n_batches,
      batch_priors$beta_beta + length(onset_times) - n_covered)
  if(!is.null(batch_priors$zeta_shape))
    batch_parameters$zeta <- stats::rgamma(
      1, shape=batch_priors$zeta_shape + sum(batches$width - 1),
      rate=batch_priors$zeta_rate + n_batches)
  if(!is.null(batch_priors$rho_batch_alpha))
    batch_parameters$rho_batch <- stats::rbeta(
      1, batch_priors$rho_batch_alpha + n_batches,
      batch_priors$rho_batch_beta +
        sum(batches$time_detected - batches$start - 1))
  batch_parameters
}
