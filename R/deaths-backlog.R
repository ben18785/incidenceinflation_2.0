# A backlog of deaths cleared in calendar time (docs/bvd_deaths.tex of the
# paper, 'A backlog of deaths cleared in calendar time'). A share b_t of the
# deaths of the cases of onset day t join a backlog instead of being recorded
# after the gamma delay from death. A death in the backlog is recorded on day
# s with hazard lambda_s, the same for every onset day:
#
#   S^U(d | t') = (1 - b_t) S_gamma(t' - d) + b_t exp(-Lambda(d, t')),
#
# where Lambda(d, t') is the hazard summed from the death (the middle of day
# d) to the end of day t'. The share is a random walk on the logit scale over
# onset weeks. The hazard is constant between consecutive reports, and its log
# is a random walk over reports with variance proportional to the time
# between them; the days up to the first report share one value.
#
# The backlog is carried in the effects on recording as the element
# 'backlog', a named list of:
# - 'logit_share': the logit of the share for each onset week;
# - 'week_of_day': the onset week of each onset day 1, 2, ...;
# - 'share': the share for each onset day 1, 2, ...;
# - 'log_rate': the log of the hazard for each report interval (the first is
#   that of every day up to the first report);
# - 'report_days': the days of the reports;
# - 'cumulative': the hazard summed over days 1, 2, ..., to the last report;
# - 'share_rw_sd' and 'rate_rw_sd': the standard deviations of the steps of
#   the two random walks (that of the hazard for seven days).

#' Starting values of the backlog
#'
#' @param n_days the number of onset days
#' @param report_days the days of the reports
#' @param backlog the element 'recording_backlog' of the result of
#' \code{deaths_settings}
#'
#' @return the element 'backlog' of the effects on recording
initial_backlog <- function(n_days, report_days, backlog) {
  report_days <- sort(unique(report_days))
  week_of_day <- (seq_len(n_days) - 1) %/% backlog$share_week + 1
  refresh_backlog(list(
    logit_share=rep(stats::qlogis(backlog$initial_share), max(week_of_day)),
    week_of_day=week_of_day,
    log_rate=rep(log(backlog$initial_rate), length(report_days)),
    report_days=report_days,
    share_rw_sd=backlog$initial_share_rw_sd,
    rate_rw_sd=backlog$initial_rate_rw_sd))
}

#' Recomputes the share of each onset day and the summed hazard of each day
#' from the values for each week and each report interval
#'
#' @param backlog the element 'backlog' of the effects on recording
#'
#' @return backlog with elements 'share' and 'cumulative' updated
refresh_backlog <- function(backlog) {
  backlog$share <- stats::plogis(backlog$logit_share[backlog$week_of_day])
  backlog$cumulative <- cumsum(backlog_daily_rates(backlog))
  backlog
}

#' The hazard of a death in the backlog being recorded on each day
#'
#' @inheritParams refresh_backlog
#' @param n_days_out the number of days for which it is wanted (by default to
#' the last report); days after the last report have the hazard of its
#' interval
#'
#' @return a vector of hazards for days 1, 2, ...
backlog_daily_rates <- function(backlog, n_days_out=max(backlog$report_days)) {
  days <- seq_len(n_days_out)
  # the first report made on or after each day
  index <- pmin(findInterval(days - 1e-8, backlog$report_days) + 1,
                length(backlog$report_days))
  exp(backlog$log_rate[index])
}

#' The hazard summed over days 1 to each of a set of days
#'
#' @inheritParams refresh_backlog
#' @param days whole numbers of days (one or more)
#'
#' @return a vector of the same length as days
backlog_cumulative <- function(backlog, days) {
  cumulative <- backlog$cumulative
  n <- length(cumulative)
  last <- exp(backlog$log_rate[length(backlog$log_rate)])
  out <- rep(0, length(days))
  is_within <- days >= 1 & days <= n
  out[is_within] <- cumulative[days[is_within]]
  is_after <- days > n
  out[is_after] <- cumulative[n] + (days[is_after] - n) * last
  out
}

#' Log-probability that a death in the backlog has not been recorded by the
#' end of each of a set of days
#'
#' The death happens z days after onset with probability h_z (or on the onset
#' day itself, where there is no 'before'). Its record waits from the middle
#' of the day of death to the end of each day:
#' \deqn{\sum_{z > k} h_z + \sum_{z \le k} h_z \exp(-[L(t') - L(d) + \lambda_d / 2])}
#' with k = t' - t, d = t + z and L the summed hazard.
#'
#' @inheritParams log_undetected_days
#' @param backlog the element 'backlog' of the effects on recording
#' @param before see \code{log_undetected_days}; NULL if the death is on the
#' onset day
#'
#' @return a vector of log-probabilities
log_backlog_days <- function(days, day_onset, backlog, before=NULL) {
  lags <- round(days - day_onset)
  h <- if(is.null(before)) 1 else c(0, before$lag_probabilities)
  # z = 0, 1, 2, ...: the probability of death z days after onset
  z_max <- length(h) - 1
  days_death <- day_onset + 0:z_max
  rates_death <- backlog_cumulative(backlog, days_death) -
    backlog_cumulative(backlog, days_death - 1)
  # the summed hazard from the middle of the day of death
  start <- backlog_cumulative(backlog, days_death) - rates_death / 2
  end <- backlog_cumulative(backlog, days)
  vapply(seq_along(days), function(i) {
    k <- lags[i]
    z <- 0:min(k, z_max)
    alive <- if(k < z_max) sum(h[(k + 2):(z_max + 1)]) else 0
    p <- alive + sum(h[z + 1] * exp(-(end[i] - start[z + 1])))
    log(max(p, 1e-300))
  }, numeric(1))
}

#' Log-probability that a death has not been recorded, where a share of
#' deaths join the backlog
#'
#' @inheritParams log_undetected_days
#'
#' @return a vector of log-probabilities
log_undetected_days_backlog <- function(days, day_onset, reporting_parameters,
                                        revision_parameters) {
  backlog <- revision_parameters$backlog
  share <- backlog$share[min(max(round(day_onset), 1), length(backlog$share))]
  revision_parameters$backlog <- NULL
  a <- log1p(-share) + log_undetected_days(days, day_onset, reporting_parameters,
                                           revision_parameters)
  if(share <= 0)
    return(a)
  b <- log(share) + log_backlog_days(days, day_onset, backlog,
                                     revision_parameters$before)
  pmax(a, b) + log1p(exp(-abs(a - b)))
}

#' Updates the backlog: the share for each onset week, the hazard for each
#' report interval and the standard deviations of the two random walks
#'
#' The logit of the share of each week is changed in turn, with only the
#' onset days of that week in the likelihood. The log hazard is changed in
#' blocks of consecutive report intervals, by a common shift (which keeps the
#' steps within the block) and by independent changes; only onset days up to
#' the end of the block are affected. The standard deviations are updated
#' given the steps of their random walks, which have half-normal priors.
#'
#' @inheritParams deaths_delay_by_onset
#' @param settings the element 'recording_backlog' of the result of
#' \code{deaths_settings}; nothing is done if it is FALSE
#'
#' @return the state with element 'backlog' of 'recording_effects' updated
deaths_recording_backlog_step <- function(state, settings) {
  effects <- state$recording_effects
  if(identical(settings, FALSE) | is.null(settings) | is.null(effects$backlog))
    return(state)
  delays <- deaths_delay_by_onset(state)
  sub <- uncorrected_state(state)
  f_logl <- function(effects, ks) {
    out <- sum(purrr::map_dbl(ks, function(k)
      deaths_recording_day_logp(state, k, delay_of_onset(delays, k), effects)))
    ks_sub <- which(sub$day %in% ks)
    if(length(ks_sub) > 0)
      out <- out + sum(purrr::map_dbl(ks_sub, function(k)
        deaths_recording_day_logp(sub, k, delay_of_onset(sub$delays, k), effects)))
    out
  }
  backlog <- effects$backlog
  onset_times <- state$onset_times
  week_of_onset <- backlog$week_of_day[pmin(pmax(round(onset_times), 1),
                                            length(backlog$week_of_day))]
  f_with <- function(backlog) {
    effects$backlog <- refresh_backlog(backlog)
    effects
  }

  # log-prior of a random walk: its first value and its steps
  f_log_prior_share <- function(x, sd_step)
    stats::dnorm(x[1], settings$share_mu, settings$share_sd, log=TRUE) +
    sum(stats::dnorm(diff(x), 0, sd_step, log=TRUE))
  gaps <- c(NA, diff(backlog$report_days))
  f_log_prior_rate <- function(x, sd_step)
    stats::dnorm(x[1], settings$rate_mu, settings$rate_sd, log=TRUE) +
    sum(stats::dnorm(diff(x), 0, sd_step * sqrt(gaps[-1] / 7), log=TRUE))
  f_accept <- function(log_ratio)
    !is.na(log_ratio) && log(stats::runif(1)) < log_ratio

  for(step in seq_len(settings$n_steps)) {
    # the share of each onset week
    for(w in sample.int(length(backlog$logit_share))) {
      ks <- which(week_of_onset == w)
      if(length(ks) == 0)
        next
      proposed <- backlog
      proposed$logit_share[w] <- backlog$logit_share[w] +
        stats::rnorm(1, 0, settings$share_step)
      log_ratio <- f_logl(f_with(proposed), ks) - f_logl(f_with(backlog), ks) +
        f_log_prior_share(proposed$logit_share, backlog$share_rw_sd) -
        f_log_prior_share(backlog$logit_share, backlog$share_rw_sd)
      if(f_accept(log_ratio))
        backlog <- proposed
    }
    # the hazard of each report interval, in blocks of consecutive intervals
    n_reports <- length(backlog$log_rate)
    n_blocks <- min(settings$n_blocks, n_reports)
    for(kind in c("shift", "each")) {
      for(block in split(seq_len(n_reports),
                         cut(seq_len(n_reports), n_blocks, labels=FALSE))) {
        ks <- which(onset_times <= backlog$report_days[max(block)])
        if(length(ks) == 0)
          next
        proposed <- backlog
        change <- if(kind == "shift") stats::rnorm(1, 0, settings$rate_step)
          else stats::rnorm(length(block), 0, settings$rate_step / 2)
        proposed$log_rate[block] <- backlog$log_rate[block] + change
        log_ratio <- f_logl(f_with(proposed), ks) - f_logl(f_with(backlog), ks) +
          f_log_prior_rate(proposed$log_rate, backlog$rate_rw_sd) -
          f_log_prior_rate(backlog$log_rate, backlog$rate_rw_sd)
        if(f_accept(log_ratio))
          backlog <- proposed
      }
    }
    # the standard deviations of the steps, with half-normal priors
    f_sd <- function(current, f_log_prior_walk, x, scale) {
      proposed <- current * exp(stats::rnorm(1, 0, settings$rw_sd_step))
      log_ratio <- f_log_prior_walk(x, proposed) - f_log_prior_walk(x, current) +
        stats::dnorm(proposed, 0, scale, log=TRUE) -
        stats::dnorm(current, 0, scale, log=TRUE) + log(proposed) - log(current)
      if(f_accept(log_ratio)) proposed else current
    }
    backlog$share_rw_sd <- f_sd(backlog$share_rw_sd, f_log_prior_share,
                                backlog$logit_share, settings$share_rw_sd_scale)
    backlog$rate_rw_sd <- f_sd(backlog$rate_rw_sd, f_log_prior_rate,
                               backlog$log_rate, settings$rate_rw_sd_scale)
  }
  state$recording_effects$backlog <- refresh_backlog(backlog)
  state
}

#' Default settings of the backlog
#'
#' @param supplied TRUE or a named list of settings to change
#'
#' @return a named list of: 'share_mu' and 'share_sd' (normal prior of the
#' logit of the share of the first onset week), 'share_rw_sd_scale' (scale of
#' the half-normal prior for the sd of its weekly steps), 'share_week' (days in
#' an onset week), 'rate_mu' and 'rate_sd' (normal prior of the log hazard
#' up to the first report), 'rate_rw_sd_scale' (scale of the half-normal prior
#' for the sd of its steps over seven days), 'initial_share', 'initial_rate',
#' 'initial_share_rw_sd', 'initial_rate_rw_sd', the steps proposed
#' ('share_step', 'rate_step', 'rw_sd_step'), 'n_blocks' and 'n_steps'
backlog_settings <- function(supplied) {
  backlog <- list(share_mu=stats::qlogis(0.2), share_sd=1, share_rw_sd_scale=0.5,
                  share_week=7, rate_mu=log(1 / 60), rate_sd=1,
                  rate_rw_sd_scale=0.5, initial_share=0.2, initial_rate=1 / 60,
                  initial_share_rw_sd=0.2, initial_rate_rw_sd=0.2,
                  share_step=0.4, rate_step=0.3, rw_sd_step=0.3, n_blocks=6,
                  n_steps=1)
  if(isTRUE(supplied))
    supplied <- list()
  unknown <- setdiff(names(supplied), names(backlog))
  if(length(unknown) > 0)
    stop(paste0("Unknown settings for recording_backlog: ",
                paste(unknown, collapse=", "), "."))
  backlog[names(supplied)] <- supplied
  positive <- c("share_sd", "share_rw_sd_scale", "share_week", "rate_sd",
                "rate_rw_sd_scale", "initial_rate", "initial_share_rw_sd",
                "initial_rate_rw_sd", "share_step", "rate_step", "rw_sd_step",
                "n_blocks", "n_steps")
  if(any(unlist(backlog[positive]) <= 0))
    stop("The scales, starting values, steps and blocks of the backlog must be positive.")
  if(backlog$initial_share <= 0 | backlog$initial_share >= 1)
    stop("The share of deaths in the backlog must start strictly between 0 and 1.")
  backlog
}

#' The share and hazard of the backlog, for storage
#'
#' @param backlog the element 'backlog' of the effects on recording
#'
#' @return a named list of values: 'backlog_share_rw_sd', 'backlog_rate_rw_sd'
#' and 'backlog_rate_<day>' for each report day
backlog_draws <- function(backlog) {
  out <- list(backlog_share_rw_sd=backlog$share_rw_sd,
              backlog_rate_rw_sd=backlog$rate_rw_sd)
  for(r in seq_along(backlog$report_days))
    out[[paste0("backlog_rate_", backlog$report_days[r])]] <- exp(backlog$log_rate[r])
  out
}
