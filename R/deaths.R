# Deaths among cases, and a case fatality ratio which changes at unknown
# onset times.
#
# Of the true cases arising on an onset time, a number eventually die: the
# true deaths. These are binomial given the true cases, with a probability
# (the case fatality ratio) which is constant within segments of onset times.
# Deaths come to be reported in the same way as cases, after their own delay
# from onset, so that the model of reporting for cases is used for them with
# the true deaths in place of the true cases.

#' Default settings for the model of deaths
#'
#' @param n_days number of onset times
#' @param settings a named list of any settings to change
#'
#' @return a named list of: 'cfr_alpha' and 'cfr_beta' (parameters of the beta
#' prior for the case fatality ratio of each segment); 'alpha', 'beta' and
#' 'min_length' (the prior on changes in the case fatality ratio, see
#' \code{changepoints_logp_prior}); 'cfr_changepoints' (FALSE holds the case
#' fatality ratio constant); 'initial_cfr'; 'initial_delay' (a named list of
#' the 'mean' and 'sd' of the delay from onset to a death being reported);
#' 'delay_prior' (as priors$reporting in \code{mcmc_single});
#' 'delay_changepoints' (a named list of settings passed to
#' \code{reporting_changepoint_settings}, or FALSE for a delay which does not
#' change); 'revisions' (NULL if reported deaths never fall, or a named list
#' of initial 'eta' and 'rho', optionally with 'weekday' and 'weekday_first',
#' as initial_revision_parameters in \code{mcmc_single});
#' 'revisions_prior' (as priors$revisions in \code{mcmc_single});
#' 'revision_steps' (a named list of 'eta_step', 'rho_step' and
#' 'weekday_step'); 'weekday_sd' (the sd of the prior of the log effect of
#' each day of the week); 'max_erroneous'; 'p_gamma_cutoff' (true deaths are
#' sampled for onset times recent enough that less than this proportion of
#' deaths would have been reported); and 'n_parameter_steps'
#' @export
deaths_settings <- function(n_days, settings=list()) {
  defaults <- list(
    cfr_alpha=1, cfr_beta=1,
    alpha=1, beta=n_days, min_length=14,
    cfr_changepoints=TRUE,
    initial_cfr=0.5,
    initial_delay=list(mean=14, sd=8),
    delay_prior=list(mean_mu=14, mean_sigma=10, sd_mu=8, sd_sigma=8),
    delay_changepoints=list(min_length=14),
    revisions=NULL,
    revisions_prior=list(eta_mu=0.005, eta_sigma=0.005, rho_alpha=1, rho_beta=1),
    revision_steps=list(eta_step=0.5, rho_step=0.5, weekday_step=0.05),
    weekday_sd=1,
    max_erroneous=3,
    p_gamma_cutoff=0.99,
    n_parameter_steps=1)
  unknown <- setdiff(names(settings), names(defaults))
  if(length(unknown) > 0)
    stop(paste0("Unknown settings for deaths: ", paste(unknown, collapse=", "), "."))
  defaults[names(settings)] <- settings
  if(defaults$cfr_alpha <= 0 | defaults$cfr_beta <= 0)
    stop("The parameters of the prior for the case fatality ratio must be positive.")
  if(defaults$alpha <= 0 | defaults$beta <= 0)
    stop("The parameters of the prior on the probability of a change must be positive.")
  if(defaults$min_length < 1)
    stop("The minimum length of a segment must be at least one day.")
  if(defaults$initial_cfr <= 0 | defaults$initial_cfr >= 1)
    stop("The initial case fatality ratio must lie strictly between 0 and 1.")
  is_delay_changes <- !identical(defaults$delay_changepoints, FALSE)
  delay_settings <- if(is_delay_changes) defaults$delay_changepoints else list()
  if(is.null(delay_settings$min_length))
    delay_settings$min_length <- 14
  if(!is_delay_changes)
    delay_settings$n_partition_moves <- 0
  if(is.null(delay_settings$n_segment_steps))
    delay_settings$n_segment_steps <- defaults$n_parameter_steps
  defaults$delay_changepoints <- reporting_changepoint_settings(n_days, delay_settings)
  defaults
}

#' Log marginal likelihoods of true deaths for every possible segment of
#' constant case fatality ratio
#'
#' With a beta prior for the case fatality ratio of a segment and binomial
#' deaths, the ratio integrates out to give a beta-binomial. Terms which are
#' the same for every partition (the binomial coefficients) are left out.
#'
#' @param deaths_true a vector of true deaths on consecutive onset times
#' @param cases_true a vector of true cases on the same onset times
#' @param cfr_alpha first parameter of the beta prior
#' @param cfr_beta second parameter of the beta prior
#'
#' @return a matrix whose element (s, t) is the log marginal likelihood of a
#' segment which starts on day s and ends on day t (-Inf if s > t)
segment_log_marginals_cfr <- function(deaths_true, cases_true, cfr_alpha=1,
                                      cfr_beta=1) {
  if(length(deaths_true) != length(cases_true))
    stop("There must be true deaths and true cases for each onset time.")
  if(any(deaths_true > cases_true) | any(deaths_true < 0))
    stop("True deaths must lie between zero and true cases.")
  n_days <- length(deaths_true)
  c_deaths <- c(0, cumsum(deaths_true))
  c_survivors <- c(0, cumsum(cases_true - deaths_true))
  f_sum <- function(cumulative)
    outer(seq_len(n_days), seq_len(n_days),
          function(s, t) pmax(cumulative[t + 1] - cumulative[s], 0))
  log_marginals <- lbeta(cfr_alpha + f_sum(c_deaths), cfr_beta + f_sum(c_survivors)) -
    lbeta(cfr_alpha, cfr_beta)
  log_marginals[lower.tri(log_marginals)] <- -Inf
  log_marginals
}

#' Samples the segments of constant case fatality ratio and the ratio in each
#'
#' The partition is drawn exactly from its distribution given the true deaths
#' and true cases, with the ratio of each segment integrated out (see
#' \code{sample_partition}). The ratio of each segment is then drawn from its
#' beta conditional distribution.
#'
#' @inheritParams segment_log_marginals_cfr
#' @inheritParams partition_forward
#' @param is_changepoints if FALSE a single segment is used
#'
#' @return a named list of 'starts' (the onset times, counted from one, on
#' which segments start) and 'cfr' (the ratio on each onset time)
#' @export
sample_cfr_changepoints <- function(deaths_true, cases_true, p_change,
                                    min_length=1, cfr_alpha=1, cfr_beta=1,
                                    is_changepoints=TRUE) {
  n_days <- length(deaths_true)
  starts <- 1
  if(is_changepoints) {
    if(p_change <= 0 | p_change >= 1)
      stop("The probability of a change must lie strictly between 0 and 1.")
    starts <- sample_partition(
      segment_log_marginals_cfr(deaths_true, cases_true, cfr_alpha, cfr_beta),
      p_change, min_length)
  }
  ends <- c(starts[-1] - 1, n_days)
  values <- purrr::map_dbl(seq_along(starts), function(j) {
    days <- starts[j]:ends[j]
    deaths <- sum(deaths_true[days])
    stats::rbeta(1, cfr_alpha + deaths, cfr_beta + sum(cases_true[days]) - deaths)
  })
  list(starts=starts, cfr=values[segment_indices(starts, n_days)])
}

#' Samples the true deaths of one onset time
#'
#' The distribution drawn from is
#' \deqn{p(deaths_true|...) \propto binomial(deaths_true|cases_true, cfr)
#' p(deaths_reported|deaths_true, delay)}
#'
#' @param observation_df a data frame with columns 'time_reported' and
#' 'cases_reported' holding the reported deaths in order of time_reported
#' @param cases_true the true cases of the onset time
#' @param cfr the case fatality ratio of the onset time
#' @param deaths_current the current true deaths
#' @inheritParams observation_process_logp
#' @param log_negligible values whose binomial log-probability is this far
#' below the largest are given zero probability without evaluating the
#' probability of the reports
#'
#' @return a count
sample_true_deaths_single_onset <- function(observation_df, cases_true, cfr,
                                            deaths_current, day_onset,
                                            reporting_parameters,
                                            revision_parameters=NULL,
                                            max_erroneous=3, log_negligible=50) {
  max_reported <- max(observation_df$cases_reported)
  lower <- max_reported
  if(!is.null(revision_parameters))
    lower <- max(0, max_reported - max_erroneous)
  if(lower > cases_true)
    return(deaths_current)
  possible <- lower:cases_true
  logps_upper <- stats::dbinom(possible, cases_true, cfr, log=TRUE)
  if(nrow(observation_df) > 1) {
    f_observation <- function(deaths)
      observation_process_logp(observation_df, deaths, day_onset,
                               reporting_parameters, revision_parameters,
                               max_erroneous)
    index_pilot <- unique(c(which.max(logps_upper),
                            which.min(abs(possible - deaths_current))))
    logps_pilot <- logps_upper[index_pilot] + f_observation(possible[index_pilot])
    is_evaluated <- logps_upper >= max(logps_pilot) - log_negligible
    logps <- rep(-Inf, length(possible))
    logps[is_evaluated] <- logps_upper[is_evaluated] +
      f_observation(possible[is_evaluated])
  } else {
    logps <- logps_upper
  }
  if(!any(is.finite(logps)))
    return(deaths_current)
  probs <- exp(logps - matrixStats::logSumExp(logps))
  possible[sample.int(length(possible), 1, prob=probs)]
}

#' Starting state of the model of deaths
#'
#' @param deaths_snapshot_df a tibble with columns: time_onset, time_reported
#' and deaths_reported, with the same onset times as the snapshots of cases
#' @param cases_true_df a tibble with columns time_onset and cases_true
#' @param settings as returned by \code{deaths_settings}
#'
#' @return a named list holding the data in the form used by the model of
#' reporting ('snapshot'), the true deaths and case fatality ratio by onset
#' time, the segments of the ratio and of the delay, and the revision
#' parameters
#' @importFrom rlang .data
initial_deaths_state <- function(deaths_snapshot_df, cases_true_df, settings) {
  if(!all(c("time_onset", "time_reported", "deaths_reported") %in%
          colnames(deaths_snapshot_df)))
    stop("The snapshots of deaths must contain columns: 'time_onset', 'time_reported' and 'deaths_reported'.")
  snapshot <- deaths_snapshot_df %>%
    dplyr::select("time_onset", "time_reported", cases_reported="deaths_reported") %>%
    dplyr::arrange(.data$time_onset, .data$time_reported)
  onset_times <- sort(unique(snapshot$time_onset))
  cases_true_df <- dplyr::arrange(unique(cases_true_df[, c("time_onset", "cases_true")]),
                                  .data$time_onset)
  if(!identical(as.numeric(onset_times), as.numeric(cases_true_df$time_onset)))
    stop("The snapshots of deaths must have the same onset times as those of cases.")
  if(any(diff(onset_times) != 1))
    stop("Onset times must be consecutive.")
  n_days <- length(onset_times)

  revisions <- settings$revisions
  if(!is.null(revisions)) {
    if(revisions$eta <= 0 | revisions$rho <= 0 | revisions$rho >= 1)
      stop("Initial revision parameters for deaths must satisfy: eta > 0 and 0 < rho < 1.")
    check_revisions_possible(snapshot, settings$max_erroneous)
    if(isTRUE(revisions$weekday)) {
      revisions$report_days <- sort(unique(snapshot$time_reported))
      if(is.null(revisions$weekday_first))
        revisions$weekday_first <- 1
      weekdays_present <- sort(unique(weekday_of_days(
        revisions$report_days, revisions$weekday_first)))
      revisions$weekday_log_effects <- stats::setNames(
        rep(0, length(weekdays_present)), weekdays_present)
      revisions <- update_report_multipliers(revisions)
    }
  } else if(max_downward_revision(snapshot) > 0) {
    stop("Reported deaths fall, so revisions must be allowed for deaths (see deaths_settings).")
  }

  # start from the larger of the deaths reported and those expected, without
  # exceeding the true cases
  reported <- snapshot %>%
    dplyr::group_by(.data$time_onset) %>%
    dplyr::summarise(largest=max(.data$cases_reported))
  deaths_true <- pmin(cases_true_df$cases_true,
                      pmax(reported$largest,
                           round(settings$initial_cfr * cases_true_df$cases_true)))
  list(
    snapshot=snapshot,
    onset_times=onset_times,
    deaths_true=deaths_true,
    cfr=rep(settings$initial_cfr, n_days),
    cfr_starts=1,
    p_change=settings$alpha / (settings$alpha + settings$beta),
    reporting=initial_reporting_changepoints(settings$initial_delay),
    revisions=revisions)
}

#' Reporting delay parameters of deaths for each onset time
#'
#' @param state as returned by \code{initial_deaths_state}
#'
#' @return a tibble with columns 'mean' and 'sd' and a row for each onset time
deaths_delay_by_onset <- function(state) {
  n_days <- length(state$onset_times)
  state$reporting$parameters[segment_indices(state$reporting$starts, n_days), ]
}

#' Log-probability of the true deaths given possible true cases of an onset time
#'
#' This is the term which the model of deaths adds to the conditional
#' distribution of the true cases.
#'
#' @inheritParams deaths_delay_by_onset
#' @param time_onset an onset time
#' @param cases_true a vector of possible true cases
#'
#' @return a vector of log-probabilities
deaths_cases_logp <- function(state, time_onset, cases_true) {
  k <- match(time_onset, state$onset_times)
  stats::dbinom(state$deaths_true[k], cases_true, state$cfr[k], log=TRUE)
}

#' Updates everything in the model of deaths given the true cases
#'
#' In turn: the true deaths of each recent onset time; the delay from onset to
#' a death being reported, and the onset times at which it changes; the
#' revision parameters of deaths and the effects of the day of the week (if
#' reported deaths can fall); and the case fatality ratio, with the onset
#' times at which it changes.
#'
#' @inheritParams deaths_delay_by_onset
#' @param cases_true a vector of the true cases on each onset time
#' @inheritParams initial_deaths_state
#'
#' @return the updated state
deaths_step <- function(state, cases_true, settings) {
  n_days <- length(state$onset_times)
  if(length(cases_true) != n_days)
    stop("There must be true cases for each onset time.")
  max_erroneous <- settings$max_erroneous
  l_observations <- split(
    as.data.frame(state$snapshot[, c("time_reported", "cases_reported")]),
    state$snapshot$time_onset)

  # true deaths
  delays <- deaths_delay_by_onset(state)
  uncertain_period <- delay_quantile(
    settings$p_gamma_cutoff, list(mean=delays$mean[n_days], sd=delays$sd[n_days]))
  time_latest <- max(state$onset_times)
  for(k in seq_len(n_days)) {
    observation_df <- l_observations[[as.character(state$onset_times[k])]]
    if(state$onset_times[k] < time_latest - uncertain_period) {
      # reporting is taken to be complete
      state$deaths_true[k] <- min(dplyr::last(observation_df$cases_reported),
                                  cases_true[k])
    } else {
      state$deaths_true[k] <- sample_true_deaths_single_onset(
        observation_df, cases_true[k], state$cfr[k], state$deaths_true[k],
        state$onset_times[k], list(mean=delays$mean[k], sd=delays$sd[k]),
        state$revisions, max_erroneous)
    }
  }

  # delay from onset to a death being reported
  by_onset <- dplyr::tibble(time_onset=state$onset_times,
                            cases_true=state$deaths_true)
  df <- dplyr::left_join(state$snapshot, by_onset, by="time_onset")
  state$reporting <- sample_reporting_changepoints(
    state$reporting,
    reporting_logp_function(df, state$revisions, max_erroneous),
    n_days, settings$delay_prior, settings$delay_changepoints)

  # revision parameters and effects of the day of the week
  if(!is.null(state$revisions)) {
    reporting_parameters <- state$reporting$parameters %>%
      dplyr::mutate(reporting_piece_index=dplyr::row_number())
    df$reporting_piece_index <- segment_indices(
      state$reporting$starts, n_days)[match(df$time_onset, state$onset_times)]
    logl <- sum(state$reporting$logl)
    logp <- logl + prior_revision_parameters(state$revisions, settings$revisions_prior) +
      log_jacobian_revision_parameters(state$revisions)
    for(step in seq_len(settings$n_parameter_steps)) {
      result <- metropolis_step_revisions(
        df, reporting_parameters, state$revisions, logp,
        settings$revisions_prior, settings$revision_steps, max_erroneous)
      state$revisions <- result$revision_parameters
      logp <- result$logp
    }
    if(!is.null(state$revisions$weekday_log_effects)) {
      logl <- logp - prior_revision_parameters(state$revisions, settings$revisions_prior) -
        log_jacobian_revision_parameters(state$revisions)
      for(step in seq_len(settings$n_parameter_steps)) {
        result <- metropolis_step_weekday(
          df, reporting_parameters, state$revisions, logl,
          list(weekday_sd=settings$weekday_sd), settings$revision_steps,
          max_erroneous)
        state$revisions <- result$revision_parameters
        logl <- result$logl
      }
    }
  }

  # case fatality ratio
  sampled <- sample_cfr_changepoints(
    state$deaths_true, cases_true, state$p_change, settings$min_length,
    settings$cfr_alpha, settings$cfr_beta, settings$cfr_changepoints)
  state$cfr_starts <- sampled$starts
  state$cfr <- sampled$cfr
  if(settings$cfr_changepoints)
    state$p_change <- sample_change_probability(
      length(sampled$starts) - 1, n_days, settings$alpha, settings$beta)
  state
}

#' The values held in the state of the model of deaths, for storage
#'
#' @inheritParams deaths_delay_by_onset
#' @param iteration the iteration
#'
#' @return a named list of two tibbles: 'deaths' (a row for each onset time
#' with the true deaths, the case fatality ratio and the delay) and
#' 'deaths_parameters' (a single row)
deaths_draws <- function(state, iteration) {
  delays <- deaths_delay_by_onset(state)
  parameters <- dplyr::tibble(
    iteration=iteration,
    n_cfr_changes=length(state$cfr_starts) - 1,
    p_change=state$p_change,
    n_delay_changes=length(state$reporting$starts) - 1)
  if(!is.null(state$revisions)) {
    parameters$eta <- state$revisions$eta
    parameters$rho <- state$revisions$rho
    for(weekday in names(state$revisions$weekday_log_effects))
      parameters[[paste0("weekday_", weekday)]] <-
        exp(state$revisions$weekday_log_effects[[weekday]])
  }
  list(
    deaths=dplyr::tibble(
      iteration=iteration,
      time_onset=state$onset_times,
      deaths_true=state$deaths_true,
      cfr=state$cfr,
      delay_mean=delays$mean,
      delay_sd=delays$sd),
    deaths_parameters=parameters)
}

#' Simulates reported cases and deaths from individual cases
#'
#' Each case is reported after a delay from onset. A case dies with a
#' probability (the case fatality ratio of its onset time) and, if it does,
#' the death is recorded after a further delay from onset. A case is shown as
#' dead in a report once it has been reported and its death recorded, so that
#' reported deaths never exceed reported cases. As elsewhere in the package,
#' an event with delay d from onset is taken to be in a report made at least
#' d - 0.5 days after onset.
#'
#' @param cases_true a vector of true cases on onset times 1, 2, ...
#' @param report_days the times at which reports are made
#' @param reporting_parameters the delay from onset to a case being reported:
#' a named list of 'mean' and 'sd' or a tibble with a row for each onset time
#' @param cfr the case fatality ratio: a single value or one for each onset time
#' @param death_parameters the delay from onset to a death being recorded, in
#' the same form as reporting_parameters
#'
#' @return a named list of 'snapshots' (a tibble with columns time_onset,
#' time_reported, cases_reported and deaths_reported) and 'truth' (a tibble
#' with columns time_onset, cases_true, deaths_true, cfr, and the mean and sd
#' of the delay from onset to being shown as dead among the cases simulated)
#' @export
#' @importFrom rlang .data
simulate_cases_deaths <- function(cases_true, report_days, reporting_parameters,
                                  cfr, death_parameters) {
  n_days <- length(cases_true)
  f_expand <- function(parameters) {
    if(length(parameters$mean) == 1)
      return(dplyr::tibble(mean=rep(parameters$mean, n_days),
                           sd=rep(parameters$sd, n_days)))
    if(length(parameters$mean) != n_days)
      stop("There must be delay parameters for each onset time.")
    dplyr::tibble(mean=parameters$mean, sd=parameters$sd)
  }
  reporting <- f_expand(reporting_parameters)
  death <- f_expand(death_parameters)
  if(length(cfr) == 1)
    cfr <- rep(cfr, n_days)
  if(length(cfr) != n_days)
    stop("There must be a case fatality ratio for each onset time.")
  f_delay <- function(n, mean, sd)
    stats::rgamma(n, shape=mean^2 / sd^2, rate=mean / sd^2)

  results <- purrr::map(seq_len(n_days), function(t) {
    n <- cases_true[t]
    delay_reported <- f_delay(n, reporting$mean[t], reporting$sd[t])
    is_fatal <- stats::runif(n) < cfr[t]
    delay_dead <- pmax(delay_reported, f_delay(n, death$mean[t], death$sd[t]))
    days <- report_days[report_days >= t]
    elapsed <- days - t + 0.5
    list(
      snapshots=dplyr::tibble(
        time_onset=t,
        time_reported=days,
        cases_reported=purrr::map_dbl(elapsed, ~sum(delay_reported <= .)),
        deaths_reported=purrr::map_dbl(elapsed, ~sum(delay_dead[is_fatal] <= .))),
      truth=dplyr::tibble(
        time_onset=t, cases_true=n, deaths_true=sum(is_fatal), cfr=cfr[t],
        delay_dead_mean=mean(delay_dead[is_fatal]),
        delay_dead_sd=stats::sd(delay_dead[is_fatal])))
  })
  list(snapshots=purrr::map_dfr(results, "snapshots"),
       truth=purrr::map_dfr(results, "truth"))
}
