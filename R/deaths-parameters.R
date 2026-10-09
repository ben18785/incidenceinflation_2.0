# Updates of the parameters of the joint model of cases and deaths which act
# on the delays or on wrong records: the effects of the day of the week on
# each delay, the ratio of the reporting delay of fatal cases to that of
# others, and the parameters governing wrong records. All are made given how
# many of the cases in each report are fatal and how many are wrong records.

#' Proposes new effects of the day of the week and accepts or rejects them
#'
#' The log effects sum to zero. All but the first have independent normal
#' priors with mean zero, and are proposed together using a normal kernel
#' (as in \code{metropolis_step_weekday}).
#'
#' @param log_effects the current log effects, named by day of the week
#' @param f_logl a function of the log effects giving the log-probability of
#' what depends on them
#' @param logl_current its value at the current log effects
#' @param weekday_sd the standard deviation of the prior for each log effect
#' @param weekday_step the standard deviation of the changes proposed
#'
#' @return a named list of 'log_effects' and 'logl'
metropolis_weekday_effects <- function(log_effects, f_logl, logl_current,
                                       weekday_sd, weekday_step) {
  n_free <- length(log_effects) - 1
  if(n_free < 1)
    return(list(log_effects=log_effects, logl=logl_current))
  free <- log_effects[-1] + stats::rnorm(n_free, 0, weekday_step)
  proposed <- stats::setNames(c(-sum(free), free), names(log_effects))
  logl_proposed <- f_logl(proposed)
  f_log_prior <- function(x) sum(stats::dnorm(x[-1], 0, weekday_sd, log=TRUE))
  log_ratio <- logl_proposed + f_log_prior(proposed) -
    logl_current - f_log_prior(log_effects)
  if(!is.na(log_ratio) && log(stats::runif(1)) < log_ratio)
    return(list(log_effects=proposed, logl=logl_proposed))
  list(log_effects=log_effects, logl=logl_current)
}

#' Starting effects of the day of the week on the delay to a death being
#' recorded
#'
#' @param report_days the days on which reports are made
#' @inheritParams weekday_of_days
#'
#' @return a named list of 'report_days', 'weekday_first',
#' 'weekday_log_effects' (zero for each day of the week on which reports are
#' made) and 'report_multipliers' (see \code{report_multipliers})
initial_recording_effects <- function(report_days, weekday_first=1) {
  report_days <- sort(unique(report_days))
  weekdays_present <- sort(unique(weekday_of_days(report_days, weekday_first)))
  update_report_multipliers(list(
    report_days=report_days,
    weekday_first=weekday_first,
    weekday_log_effects=stats::setNames(rep(0, length(weekdays_present)),
                                        weekdays_present)))
}

#' Updates the effects of the day of the week on the delay to a death being
#' recorded
#'
#' @inheritParams deaths_delay_by_onset
#' @param settings as returned by \code{deaths_settings}, with elements
#' 'weekday_sd' and 'weekday_step' added
#'
#' @return the state with element 'recording_effects' updated
deaths_recording_weekday_step <- function(state, settings) {
  effects <- state$recording_effects
  if(is.null(effects))
    return(state)
  delays <- deaths_delay_by_onset(state)
  f_effects <- function(log_effects) {
    effects$weekday_log_effects <- log_effects
    update_report_multipliers(effects)
  }
  f_logl <- function(log_effects) {
    proposed <- f_effects(log_effects)
    sum(purrr::map_dbl(seq_along(state$onset_times), function(k)
      deaths_recording_day_logp(state, k, delay_of_onset(delays, k), proposed)))
  }
  current <- list(log_effects=effects$weekday_log_effects)
  current$logl <- f_logl(current$log_effects)
  for(step in seq_len(settings$n_parameter_steps))
    current <- metropolis_weekday_effects(
      current$log_effects, f_logl, current$logl, settings$weekday_sd,
      settings$weekday_step)
  state$recording_effects <- f_effects(current$log_effects)
  state
}

#' Updates a parameter which is zero with a fixed prior probability and
#' otherwise has a normal prior with mean zero
#'
#' From zero, a value drawn from that normal prior is proposed. From any other
#' value, either zero is proposed or a small change is, with equal
#' probability. Each proposal is accepted with the Metropolis-Hastings
#' probability.
#'
#' @param value the current value
#' @param f_logl a function of the value giving the log-probability of what
#' depends on it
#' @param probability_zero the prior probability that the value is zero
#' @param sd the standard deviation of the normal prior otherwise
#' @param step the standard deviation of the small changes proposed
#' @param n_steps the number of proposals
#'
#' @return the updated value
sample_spike_slab <- function(value, f_logl, probability_zero=0.5, sd=0.5,
                              step=0.1, n_steps=1) {
  log_odds_zero <- stats::qlogis(probability_zero)
  logl <- f_logl(value)
  for(i in seq_len(n_steps)) {
    if(value == 0) {
      proposed <- stats::rnorm(1, 0, sd)
      # the proposal is the prior, which cancels; zero is proposed half the
      # time from other values
      log_ratio_extra <- -log_odds_zero + log(0.5)
    } else if(stats::runif(1) < 0.5) {
      proposed <- 0
      log_ratio_extra <- log_odds_zero - log(0.5)
    } else {
      proposed <- stats::rnorm(1, value, step)
      log_ratio_extra <- stats::dnorm(proposed, 0, sd, log=TRUE) -
        stats::dnorm(value, 0, sd, log=TRUE)
    }
    logl_proposed <- f_logl(proposed)
    log_ratio <- logl_proposed - logl + log_ratio_extra
    if(!is.na(log_ratio) && log(stats::runif(1)) < log_ratio) {
      value <- proposed
      logl <- logl_proposed
    }
  }
  value
}

#' Updates the ratio of the reporting delay of fatal cases to that of others
#'
#' The log of the ratio, nu, is zero with a fixed prior probability and
#' otherwise has a normal prior with mean zero (see
#' \code{sample_spike_slab}).
#'
#' @inheritParams deaths_paths_step
#' @param fatal_delay the element of that name of the result of
#' \code{deaths_settings}; nothing is done if it is FALSE
#'
#' @return the state with element 'nu' updated
deaths_fatal_delay_step <- function(state, cases_true, reporting_by_onset,
                                    revision_parameters=NULL, fatal_delay) {
  if(identical(fatal_delay, FALSE))
    return(state)
  state$nu <- sample_spike_slab(
    state$nu,
    function(nu) deaths_reporting_logp(state, cases_true, reporting_by_onset,
                                       revision_parameters, nu, is_fatal_only=TRUE),
    fatal_delay$probability_same, fatal_delay$sd, fatal_delay$step,
    fatal_delay$n_steps)
  state
}

#' Updates the parameters governing wrong records
#'
#' Given the wrong records entered and removed between reports: the
#' probability that a wrong record is shown as dead is drawn from its beta
#' conditional distribution; eta, the expected number of wrong records
#' entered for each genuine case reported, is drawn from its gamma
#' conditional distribution; and rho, the daily probability a wrong record is
#' removed, is updated by Metropolis steps on the logit scale.
#'
#' @inheritParams deaths_paths_step
#' @param settings as returned by \code{deaths_settings}
#' @param prior_parameters named list with elements 'eta_mu', 'eta_sigma',
#' 'rho_alpha' and 'rho_beta' (see \code{prior_revision_parameters})
#' @param rho_step the standard deviation of the changes proposed in the
#' logit of rho
#'
#' @return a named list of the updated 'state' and 'revision_parameters'
deaths_wrong_step <- function(state, cases_true, reporting_by_onset,
                              revision_parameters, settings, prior_parameters,
                              rho_step=0.5) {
  if(!is_wrong_records(revision_parameters))
    return(list(state=state, revision_parameters=revision_parameters))
  n_days <- length(state$onset_times)
  dated <- dated_counts(state)
  cases_recorded <- cases_true - dated$out + dated$into
  f_total <- function(name) sum(purrr::map_dbl(state$paths, ~sum(.[[name]])))
  entered_alive <- f_total("wrong_entered_alive")
  entered_dead <- f_total("wrong_entered_dead")
  state$delta <- stats::rbeta(1, settings$wrong_dead_alpha + entered_dead,
                              settings$wrong_dead_beta + entered_alive)

  exposure <- sum(cases_recorded *
                    purrr::map_dbl(state$wrong_units, ~sum(.$unit)))
  revision_parameters$eta <- stats::rgamma(
    1,
    shape=prior_parameters$eta_mu^2 / prior_parameters$eta_sigma^2 +
      entered_alive + entered_dead,
    rate=prior_parameters$eta_mu / prior_parameters$eta_sigma^2 + exposure)

  # rho acts on the removal of wrong records and on how many of those entered
  # before the first report are still present at it
  f_logp <- function(rho) {
    proposed <- revision_parameters
    proposed$rho <- rho
    logl <- sum(purrr::map_dbl(seq_len(n_days), function(k) {
      path <- state$paths[[k]]
      units <- wrong_entry_units(state$reports[[k]]$time_reported,
                                 state$onset_times[k],
                                 delay_of_onset(reporting_by_onset, k), proposed)
      sum(dbinom_log_recycled(path$wrong_removed, path$wrong_before, units$w)) +
        stats::dpois(path$wrong_entered_alive[1] + path$wrong_entered_dead[1],
                     proposed$eta * cases_recorded[k] * units$unit[1], log=TRUE)
    }))
    logl + stats::dbeta(rho, prior_parameters$rho_alpha,
                        prior_parameters$rho_beta, log=TRUE) +
      log(rho) + log1p(-rho)
  }
  rho <- revision_parameters$rho
  logp <- f_logp(rho)
  for(step in seq_len(settings$n_parameter_steps)) {
    proposed <- stats::plogis(stats::rnorm(1, stats::qlogis(rho), rho_step))
    if(proposed <= 0 | proposed >= 1)
      next
    logp_proposed <- f_logp(proposed)
    if(is.finite(logp_proposed) && log(stats::runif(1)) < logp_proposed - logp) {
      rho <- proposed
      logp <- logp_proposed
    }
  }
  revision_parameters$rho <- rho
  list(state=state, revision_parameters=revision_parameters)
}

#' Updates the effects of the day of the week on the reporting delay in the
#' joint model of cases and deaths
#'
#' @inheritParams deaths_paths_step
#' @param settings as returned by \code{deaths_settings}, with elements
#' 'weekday_sd' and 'weekday_step' added
#'
#' @return the revision parameters with elements 'weekday_log_effects' and
#' 'report_multipliers' updated
deaths_reporting_weekday_step <- function(state, cases_true, reporting_by_onset,
                                          revision_parameters, settings) {
  if(is.null(revision_parameters$weekday_log_effects))
    return(revision_parameters)
  f_parameters <- function(log_effects) {
    revision_parameters$weekday_log_effects <- log_effects
    update_report_multipliers(revision_parameters)
  }
  f_logl <- function(log_effects)
    deaths_reporting_logp(state, cases_true, reporting_by_onset,
                          f_parameters(log_effects))
  current <- list(log_effects=revision_parameters$weekday_log_effects)
  current$logl <- f_logl(current$log_effects)
  for(step in seq_len(settings$n_parameter_steps))
    current <- metropolis_weekday_effects(
      current$log_effects, f_logl, current$logl, settings$weekday_sd,
      settings$weekday_step)
  f_parameters(current$log_effects)
}

#' Updates everything which acts on the reporting delay or on wrong records
#' in the joint model of cases and deaths, other than the reporting delay
#' itself
#'
#' These are the parameters governing wrong records, the effects of the day
#' of the week on the reporting delay and the ratio of the reporting delay of
#' fatal cases to that of others. What is held about the delays of each onset
#' time is then recalculated.
#'
#' @inheritParams deaths_wrong_step
#' @param settings as returned by \code{deaths_settings}, with elements
#' 'weekday_sd' and 'weekday_step' added
#'
#' @return a named list of the updated 'state' and 'revision_parameters'
deaths_reporting_effects_step <- function(state, cases_true, reporting_by_onset,
                                          revision_parameters, settings,
                                          prior_parameters=NULL, rho_step=0.5) {
  updated <- deaths_wrong_step(state, cases_true, reporting_by_onset,
                               revision_parameters, settings, prior_parameters,
                               rho_step)
  state <- updated$state
  revision_parameters <- deaths_reporting_weekday_step(
    state, cases_true, reporting_by_onset, updated$revision_parameters, settings)
  state <- deaths_fatal_delay_step(state, cases_true, reporting_by_onset,
                                   revision_parameters, settings$fatal_delay)
  list(state=deaths_refresh_survivals(state, reporting_by_onset,
                                      revision_parameters),
       revision_parameters=revision_parameters)
}
