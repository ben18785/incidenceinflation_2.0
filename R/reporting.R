# Extensions to the model of reporting delays: a heavier-tailed delay
# distribution, an effect of the day of the week of each report, and a random
# effect of each report. Each is optional, and is carried within the revision
# parameters:
#
# - 'tail': a named list of 'weight' and 'mean'. A proportion 'weight' of cases
#   have a delay which is exponential with mean 'mean' rather than gamma.
# - 'report_multipliers': a named list of 'day_first' and 'values', giving for
#   each day from day_first onwards a multiplier on the rate at which cases not
#   yet detected are detected on that day. The multiplier for a day is that of
#   the first report made on or after it: the product of the effect of the
#   report's day of the week and the report's own random effect.

#' Log-probability a case remains undetected at the end of each of a set of days
#'
#' The probability that a case arising on day_onset remains undetected at the
#' end of a day is S(day + 0.5 - day_onset), where S is the survival function
#' of the reporting delay. The delay is gamma, optionally mixed with an
#' exponential tail. If multipliers are supplied, the rate of detection on each
#' day is multiplied by that day's multiplier, so that the log-probability of
#' remaining undetected becomes:
#' \deqn{log S(0.5) - \sum_k m_k (log S(k - 1) - log S(k))}
#' where the sum is over days after the onset day.
#'
#' @param days vector of days, each a whole number of days on or after day_onset
#' @param day_onset the onset day
#' @param reporting_parameters named list of 'mean' and 'sd' of the gamma delay
#' @param revision_parameters named list which can contain 'tail' and
#' 'report_multipliers' (see above); or NULL
#'
#' @return a vector of log-probabilities
log_undetected_days <- function(days, day_onset, reporting_parameters,
                                revision_parameters=NULL) {
  shape <- reporting_parameters$mean^2 / reporting_parameters$sd^2
  rate <- reporting_parameters$mean / reporting_parameters$sd^2
  tail <- revision_parameters$tail
  f_log_base <- function(days_tmp) {
    delay <- days_tmp + 0.5 - day_onset
    log_gamma <- stats::pgamma(delay, shape=shape, rate=rate,
                               lower.tail=FALSE, log.p=TRUE)
    if(is.null(tail) || tail$weight == 0)
      return(log_gamma)
    a <- log1p(-tail$weight) + log_gamma
    b <- log(tail$weight) - delay / tail$mean
    pmax(a, b) + log1p(exp(-abs(a - b)))
  }

  multipliers <- revision_parameters$report_multipliers
  n_days <- round(max(days) - day_onset)
  if(is.null(multipliers) | n_days < 1)
    return(f_log_base(days))

  days_after <- day_onset + seq_len(n_days)
  log_base <- f_log_base(c(day_onset, days_after))
  index <- days_after - multipliers$day_first + 1
  values <- rep(1, n_days)
  is_within <- index >= 1 & index <= length(multipliers$values)
  values[is_within] <- multipliers$values[index[is_within]]
  log_effective <- log_base[1] - c(0, cumsum(values * -diff(log_base)))
  log_effective[round(days - day_onset) + 1]
}

#' Day of the week of each of a set of days
#'
#' @param days vector of days
#' @param weekday_first day of the week (a number from 1 to 7) of day 1
#'
#' @return a vector of numbers from 1 to 7
weekday_of_days <- function(days, weekday_first=1) {
  (weekday_first - 1 + days - 1) %% 7 + 1
}

#' Builds the multiplier on the rate of detection for each day
#'
#' The multiplier for a day is that of the first report made on or after it;
#' days after the last report have a multiplier of one.
#'
#' @param report_days sorted vector of the days on which reports are made
#' @param weekday_log_effects log of the effect of each day of the week on
#' which reports are made, named by the number of the day of the week; or NULL
#' @param report_log_effects log of the random effect of each report, in the
#' order of report_days; or NULL
#' @inheritParams weekday_of_days
#'
#' @return a named list of 'day_first' and 'values'
report_multipliers <- function(report_days, weekday_log_effects=NULL,
                               report_log_effects=NULL, weekday_first=1) {
  log_effects <- rep(0, length(report_days))
  if(!is.null(weekday_log_effects))
    log_effects <- log_effects + weekday_log_effects[
      as.character(weekday_of_days(report_days, weekday_first))]
  if(!is.null(report_log_effects))
    log_effects <- log_effects + report_log_effects
  day_first <- min(report_days)
  days <- seq(day_first, max(report_days), 1)
  # index of the first report made on or after each day
  index_report <- findInterval(days - 1e-8, report_days) + 1
  list(day_first=day_first,
       values=unname(exp(log_effects[index_report])))
}

#' Refreshes the daily multipliers held within the revision parameters
#'
#' @param revision_parameters named list containing 'report_days' and
#' optionally 'weekday_log_effects', 'report_log_effects' and 'weekday_first'
#'
#' @return revision_parameters with element 'report_multipliers' updated
update_report_multipliers <- function(revision_parameters) {
  has_weekday <- !is.null(revision_parameters$weekday_log_effects)
  has_report <- !is.null(revision_parameters$report_log_effects)
  if(!has_weekday & !has_report) {
    revision_parameters$report_multipliers <- NULL
    return(revision_parameters)
  }
  weekday_first <- revision_parameters$weekday_first
  if(is.null(weekday_first))
    weekday_first <- 1
  revision_parameters$report_multipliers <- report_multipliers(
    revision_parameters$report_days,
    revision_parameters$weekday_log_effects,
    revision_parameters$report_log_effects,
    weekday_first)
  revision_parameters
}

#' Delay by which a given proportion of cases have been detected
#'
#' Multipliers on the rate of detection are not allowed for.
#'
#' @param p a proportion
#' @inheritParams log_undetected_days
#'
#' @return a delay in days
delay_quantile <- function(p, reporting_parameters, revision_parameters=NULL) {
  if(is.null(revision_parameters$tail) || revision_parameters$tail$weight == 0)
    return(qgamma_mean_sd(p, reporting_parameters$mean, reporting_parameters$sd))
  if(p <= 0)
    return(0)
  rev_tmp <- list(tail=revision_parameters$tail)
  f <- function(delay)
    exp(log_undetected_days(delay - 0.5, 0, reporting_parameters, rev_tmp)) - (1 - p)
  upper <- 10
  while(f(upper) > 0 & upper < 1e4)
    upper <- upper * 2
  stats::uniroot(f, c(0, upper))$root
}

#' Grid of values of the log of a report's random effect
#'
#' Random effects of reports are restricted to this grid, which allows them to
#' be proposed from an approximation to their conditional distribution.
#'
#' @return a vector of values
report_log_effect_grid <- function() {
  seq(-3, 3, length.out=61)
}

#' Samples the random effects of reports
#'
#' The log of each report's effect has a normal prior with mean zero and
#' standard deviation report_sd, restricted to \code{report_log_effect_grid}.
#' New values are proposed for a block of reports at a time, from an
#' approximation to their conditional distribution in which reported cases are
#' taken never to be revised downwards, so that the cases newly reported in
#' each interval are binomial. The proposal is accepted or rejected using the
#' exact probability of the observations, so that the approximation affects
#' only how often proposals are accepted.
#'
#' @param snapshot_with_true_cases_df a tibble with columns: time_onset,
#' time_reported, cases_reported, cases_true, reporting_piece_index
#' @param reporting_parameters a tibble with columns: 'reporting_piece_index', 'mean', 'sd'
#' @param revision_parameters named list of revision parameters containing
#' 'report_days', 'report_log_effects' and 'report_sd'
#' @param logl_current log-probability of the observations at current values
#' @inheritParams erroneous_transitions
#' @param n_blocks number of blocks into which reports are divided
#'
#' @return a named list of 'revision_parameters', 'logl' and 'n_accepted'
#' @importFrom rlang .data
sample_report_effects <- function(snapshot_with_true_cases_df, reporting_parameters,
                                  revision_parameters, logl_current,
                                  max_erroneous=3, n_blocks=2) {
  grid <- report_log_effect_grid()
  report_days <- revision_parameters$report_days

  # intervals between reports, with counts made non-decreasing
  intervals_df <- snapshot_with_true_cases_df %>%
    dplyr::left_join(
      reporting_parameters[, c("reporting_piece_index", "mean", "sd")],
      by="reporting_piece_index") %>%
    dplyr::group_by(.data$time_onset) %>%
    dplyr::arrange(.data$time_reported, .by_group=TRUE) %>%
    dplyr::mutate(
      cases=rev(cummin(rev(.data$cases_reported))),
      time_previous=dplyr::lag(.data$time_reported),
      newly=.data$cases - dplyr::lag(.data$cases),
      remaining=pmax(.data$cases_true, max(.data$cases)) - dplyr::lag(.data$cases)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(.data$time_previous), .data$remaining > 0)

  # rate of detection accumulated over each interval, apart from the random
  # effect of the report which ends it
  revision_tmp <- revision_parameters
  revision_tmp$report_log_effects <- NULL
  revision_tmp <- update_report_multipliers(revision_tmp)
  hazard <- purrr::pmap_dbl(
    list(intervals_df$time_previous, intervals_df$time_reported,
         intervals_df$time_onset, intervals_df$mean, intervals_df$sd),
    function(day_1, day_2, day_onset, mean, sd) {
      log_s <- log_undetected_days(c(day_1, day_2), day_onset,
                                   list(mean=mean, sd=sd), revision_tmp)
      max(log_s[1] - log_s[2], 1e-12)
    })

  # approximate conditional distribution of each report's log effect
  index_report <- match(intervals_df$time_reported, report_days)
  m_logq <- matrix(0, nrow=length(report_days), ncol=length(grid))
  for(g in seq_along(grid)) {
    log_not <- -exp(grid[g]) * hazard
    cell <- ifelse(intervals_df$newly == 0,
                   intervals_df$remaining * log_not,
                   intervals_df$newly * log(-expm1(log_not)) +
                     (intervals_df$remaining - intervals_df$newly) * log_not)
    sums <- tapply(cell, index_report, sum)
    m_logq[as.integer(names(sums)), g] <- sums
  }
  log_prior <- stats::dnorm(grid, 0, revision_parameters$report_sd, log=TRUE)
  m_logq <- sweep(m_logq, 2, log_prior, "+")
  m_logq <- m_logq - matrixStats::rowLogSumExps(m_logq)

  f_logl <- function(revision_tmp) {
    observation_process_all_times_logp(
      snapshot_with_true_cases_df=snapshot_with_true_cases_df,
      reporting_parameters=reporting_parameters,
      revision_parameters=revision_tmp,
      max_erroneous=max_erroneous)
  }

  n_reports <- length(report_days)
  blocks <- split(seq_len(n_reports),
                  sample(rep(seq_len(n_blocks), length.out=n_reports)))
  n_accepted <- 0
  for(block in blocks) {
    index_current <- match(round(revision_parameters$report_log_effects[block], 8),
                           round(grid, 8))
    index_proposed <- purrr::map_int(
      block, ~sample(length(grid), 1, prob=exp(m_logq[., ])))
    proposed <- revision_parameters
    proposed$report_log_effects[block] <- grid[index_proposed]
    proposed <- update_report_multipliers(proposed)
    logl_proposed <- f_logl(proposed)
    log_ratio <- logl_proposed - logl_current +
      sum(log_prior[index_proposed] - log_prior[index_current]) +
      sum(m_logq[cbind(block, index_current)] - m_logq[cbind(block, index_proposed)])
    if(!is.na(log_ratio) && log(stats::runif(1)) < log_ratio) {
      revision_parameters <- proposed
      logl_current <- logl_proposed
      n_accepted <- n_accepted + 1
    }
  }
  list(revision_parameters=revision_parameters, logl=logl_current,
       n_accepted=n_accepted)
}

#' Samples the standard deviation of the log of the random effects of reports
#'
#' With an inverse-gamma prior for the variance, the conditional posterior is:
#' \deqn{sigma^2 \sim inverse-gamma(shape + n / 2, scale + \sum_r (log z_r)^2 / 2)}
#' (ignoring the restriction of effects to a grid).
#'
#' @param report_log_effects vector of the log of each report's effect
#' @param prior_parameters named list of 'report_sd_shape' and
#' 'report_sd_scale', the parameters of the inverse-gamma prior for the variance
#'
#' @return a standard deviation
sample_report_sd <- function(report_log_effects, prior_parameters) {
  shape <- prior_parameters$report_sd_shape + length(report_log_effects) / 2
  scale <- prior_parameters$report_sd_scale + sum(report_log_effects^2) / 2
  sqrt(1 / stats::rgamma(1, shape=shape, rate=scale))
}

#' Samples the tail of the delay distribution using a Metropolis step
#'
#' The weight has a beta prior and the mean a gamma prior. Proposals are made
#' on the logit and log scales respectively.
#'
#' @inheritParams sample_report_effects
#' @param prior_parameters named list of 'tail_weight_alpha',
#' 'tail_weight_beta', 'tail_mean_mu' and 'tail_mean_sigma'
#' @param metropolis_parameters named list of 'tail_weight_step' and 'tail_mean_step'
#'
#' @return a named list of 'revision_parameters' and 'logl'
metropolis_step_tail <- function(snapshot_with_true_cases_df, reporting_parameters,
                                 revision_parameters, logl_current,
                                 prior_parameters, metropolis_parameters,
                                 max_erroneous=3) {
  f_log_prior <- function(tail) {
    stats::dbeta(tail$weight, prior_parameters$tail_weight_alpha,
                 prior_parameters$tail_weight_beta, log=TRUE) +
      dgamma_mean_sd(tail$mean, prior_parameters$tail_mean_mu,
                     prior_parameters$tail_mean_sigma, log=TRUE) +
      # Jacobian of transformation to logit and log scales
      log(tail$weight) + log(1 - tail$weight) + log(tail$mean)
  }
  current <- revision_parameters$tail
  proposed <- revision_parameters
  proposed$tail <- list(
    weight=stats::plogis(stats::rnorm(1, stats::qlogis(current$weight),
                                      metropolis_parameters$tail_weight_step)),
    mean=exp(stats::rnorm(1, log(current$mean),
                          metropolis_parameters$tail_mean_step)))
  logl_proposed <- observation_process_all_times_logp(
    snapshot_with_true_cases_df=snapshot_with_true_cases_df,
    reporting_parameters=reporting_parameters,
    revision_parameters=proposed,
    max_erroneous=max_erroneous)
  log_ratio <- logl_proposed + f_log_prior(proposed$tail) -
    logl_current - f_log_prior(current)
  if(!is.na(log_ratio) && log(stats::runif(1)) < log_ratio)
    return(list(revision_parameters=proposed, logl=logl_proposed))
  list(revision_parameters=revision_parameters, logl=logl_current)
}

#' Samples the effects of the day of the week using a Metropolis step
#'
#' The log effects of the days of the week on which reports are made sum to
#' zero. All but the first have independent normal priors with mean zero, and
#' are proposed together using a normal kernel.
#'
#' @inheritParams sample_report_effects
#' @param prior_parameters named list containing 'weekday_sd', the standard
#' deviation of the prior for each log effect
#' @param metropolis_parameters named list containing 'weekday_step'
#'
#' @return a named list of 'revision_parameters' and 'logl'
metropolis_step_weekday <- function(snapshot_with_true_cases_df, reporting_parameters,
                                    revision_parameters, logl_current,
                                    prior_parameters, metropolis_parameters,
                                    max_erroneous=3) {
  current <- revision_parameters$weekday_log_effects
  n_free <- length(current) - 1
  if(n_free < 1)
    return(list(revision_parameters=revision_parameters, logl=logl_current))
  free <- current[-1] + stats::rnorm(n_free, 0, metropolis_parameters$weekday_step)
  proposed <- revision_parameters
  proposed$weekday_log_effects <- stats::setNames(c(-sum(free), free), names(current))
  proposed <- update_report_multipliers(proposed)
  logl_proposed <- observation_process_all_times_logp(
    snapshot_with_true_cases_df=snapshot_with_true_cases_df,
    reporting_parameters=reporting_parameters,
    revision_parameters=proposed,
    max_erroneous=max_erroneous)
  f_log_prior <- function(log_effects)
    sum(stats::dnorm(log_effects[-1], 0, prior_parameters$weekday_sd, log=TRUE))
  log_ratio <- logl_proposed + f_log_prior(proposed$weekday_log_effects) -
    logl_current - f_log_prior(current)
  if(!is.na(log_ratio) && log(stats::runif(1)) < log_ratio)
    return(list(revision_parameters=proposed, logl=logl_proposed))
  list(revision_parameters=revision_parameters, logl=logl_current)
}
