# Predictions of later reports from a fit of the joint model of cases and
# deaths.
#
# Each draw kept by mcmc_single, with the state the sampler holds at its end,
# gives the values of everything the reports depend on but the split of the
# reports of each day: how many of the cases in each report are fatal (of
# each kind), how many are wrong records and how many of each are shown as
# dead. That split is drawn from its distribution given the rest (as the
# sampler draws it), and the reports still to come are then simulated from
# the last one, with the effects of reports not yet made drawn from their
# priors. The result is a draw from the posterior predictive distribution of
# the later reports.

#' The state of the model of deaths at a kept draw
#'
#' @param variables the variables of the sampler's state at the end of a
#' chain (element 'variables' of the 'state' returned by \code{mcmc_single}
#' with return_state TRUE), which supply what does not change between draws
#' @param draw a named list of the draws of mcmc_single for one iteration:
#' 'cases', 'reporting_days', 'revisions', 'reporting_effects', 'deaths' and
#' 'deaths_parameters', each with the rows of that iteration
#'
#' @return a named list of 'deaths_state', 'revision_parameters',
#' 'reporting_by_onset' (a tibble of 'mean' and 'sd', a row for each onset
#' time), 'cases_true' (the true cases by day of onset), 'max_erroneous' and
#' 'weekday_sd' (named 'cases' and 'recording')
#' @export
predict_state_from_draw <- function(variables, draw) {
  if(!isTRUE(variables$is_deaths))
    stop("Predictions of later reports need a fit with a model of deaths.")
  state <- variables$deaths_state
  onsets <- state$onset_times
  f_onset <- function(df, name) {
    values <- df[[name]][match(onsets, df$time_onset)]
    if(any(is.na(values)))
      stop(paste0("The draw has no '", name, "' for some onset times."))
    values
  }
  parameters <- draw$deaths_parameters
  if(nrow(parameters) != 1 || nrow(draw$revisions) != 1)
    stop("The draw must be of a single iteration.")

  # the reporting delay of each onset time and the true cases
  reporting_by_onset <- dplyr::tibble(mean=f_onset(draw$reporting_days, "mean"),
                                      sd=f_onset(draw$reporting_days, "sd"))
  cases_true <- f_onset(draw$cases, "cases_true")

  # wrong records and the effects of the day of the week on reporting
  revision <- variables$revision_current
  revision$lag_cache <- NULL
  revision$eta <- draw$revisions$eta
  revision$rho <- draw$revisions$rho
  values_drawn <- c(as.list(draw$reporting_effects), as.list(parameters))
  f_weekday <- function(prefix, current) {
    values <- vapply(names(current), function(w)
      log(values_drawn[[paste0(prefix, w)]]), numeric(1))
    stats::setNames(values, names(current))
  }
  if(!is.null(revision$weekday_log_effects)) {
    revision$weekday_log_effects <- f_weekday("weekday_", revision$weekday_log_effects)
    revision <- update_report_multipliers(revision)
  }
  if(!is.null(revision$tail)) {
    revision$tail$weight <- draw$reporting_effects$tail_weight
    revision$tail$mean <- draw$reporting_effects$tail_mean
  }

  # the true fatal cases, those dated by death and the delay to a death
  # being recorded
  deaths <- draw$deaths
  state$deaths_true <- f_onset(deaths, "deaths_true")
  state$cfr <- f_onset(deaths, "cfr")
  state$delays <- dplyr::tibble(mean=f_onset(deaths, "delay_mean"),
                                sd=f_onset(deaths, "delay_sd"))
  if(!is.null(state$dated)) {
    into <- f_onset(deaths, "dated_by_death")
    state$dated$into <- into
    state$dated$out <- state$deaths_true - f_onset(deaths, "deaths_recorded") + into
  }
  state$nu <- log(parameters$fatal_delay_ratio)
  state$triggered <- parameters$triggered
  state$delta <- parameters$wrong_dead

  # what acts on the delay to a death being recorded
  effects <- state$recording_effects
  if(!is.null(effects$weekday_log_effects)) {
    effects$weekday_log_effects <- f_weekday("recording_weekday_",
                                             effects$weekday_log_effects)
    effects <- update_report_multipliers(effects)
  }
  if(!is.null(effects$tail)) {
    effects$tail$weight <- parameters$recording_tail_weight
    effects$tail$mean <- parameters$recording_tail_mean
  }
  if(!is.null(effects$report_log_effects))
    stop("Predictions with effects of each report on the delay to a death being recorded are not yet allowed for.")
  backlog <- effects$backlog
  if(!is.null(backlog)) {
    backlog$log_rate <- log(vapply(backlog$report_days, function(day)
      parameters[[paste0("backlog_rate_", day)]], numeric(1)))
    # the share is held for each onset week
    share <- f_onset(deaths, "backlog_share")
    days <- pmin(pmax(round(onsets), 1), length(backlog$week_of_day))
    weeks <- backlog$week_of_day[days]
    backlog$logit_share[weeks] <- stats::qlogis(share)
    backlog$share_rw_sd <- parameters$backlog_share_rw_sd
    backlog$rate_rw_sd <- parameters$backlog_rate_rw_sd
    effects$backlog <- refresh_backlog(backlog)
  }
  state$recording_effects <- effects
  # what was held about the reports of each onset time is made again
  state$paths <- NULL
  state$survivals <- NULL
  state$wrong_units <- NULL
  state$wrong <- NULL

  weekday_sd <- c(cases=variables$priors_reporting$weekday_sd,
                  recording=variables$deaths_options$weekday_sd)
  list(deaths_state=state, revision_parameters=revision,
       reporting_by_onset=reporting_by_onset, cases_true=cases_true,
       max_erroneous=variables$max_erroneous, weekday_sd=weekday_sd)
}

#' The survivals and wrong records of one onset time at a set of report days
#'
#' @param prediction as returned by \code{predict_state_from_draw}
#' @param k the onset time, counted from one
#' @param days the report days
#'
#' @return a named list of 'survivals' (see \code{joint_survivals}) and
#' 'wrong' (as in \code{joint_path_forward}, or NULL without wrong records)
predict_onset_terms <- function(prediction, k, days) {
  state <- prediction$deaths_state
  revision <- prediction$revision_parameters
  reporting_parameters <- delay_of_onset(prediction$reporting_by_onset, k)
  survivals <- joint_survivals(
    days, state$onset_times[k], reporting_parameters,
    delay_of_onset(deaths_delay_by_onset(state), k), revision,
    state$dated$reporting, state$recording_effects, state$nu, state$triggered)
  wrong <- NULL
  if(is_wrong_records(revision)) {
    dated <- dated_counts(state)
    units <- wrong_entry_units(days, state$onset_times[k], reporting_parameters,
                               revision)
    wrong <- list(
      lambda=revision$eta *
        (prediction$cases_true[k] - dated$out[k] + dated$into[k]) * units$unit,
      w=units$w, delta=state$delta, max_erroneous=prediction$max_erroneous)
  }
  list(survivals=survivals, wrong=wrong)
}

#' Draws the split of the reports of one onset time given the rest of the
#' state
#'
#' @inheritParams predict_onset_terms
#'
#' @return as \code{joint_path_exact_backward}
predict_onset_path <- function(prediction, k) {
  state <- prediction$deaths_state
  reports <- state$reports[[k]]
  terms <- predict_onset_terms(prediction, k, reports$time_reported)
  dated <- dated_counts(state)
  deaths_true <- state$deaths_true[k]
  survivors <- prediction$cases_true[k] - deaths_true
  if(isTRUE(state$dated$exact)) {
    # the wrong records are not held fixed: the split is drawn whole, from
    # its distribution given the rest
    forward <- joint_path_exact_forward(
      reports$cases_reported, reports$deaths_reported, deaths_true - dated$out[k],
      dated$into[k], survivors, terms$survivals, terms$wrong)
    if(!is.finite(forward$logp))
      stop(paste0("The reports of onset time ", state$onset_times[k],
                  " are impossible at the draw."))
    return(joint_path_exact_backward(forward))
  }
  stop("Predictions need the exact form of fatal cases dated by death.")
}

#' Adds reports not yet made to what acts on the delays, drawing the effects
#' they have from their priors
#'
#' Days after the last report take the effect of the day of the week of the
#' first report made on or after them, as days before it do; a day of the week
#' on which no report had been made has its effect drawn from its prior. The
#' log of the backlog's hazard carries on its random walk over the new
#' reports.
#'
#' @inheritParams predict_onset_terms
#' @param days_future the days of the later reports
#'
#' @return prediction with its revision parameters and the effects on the
#' recording of deaths extended
predict_extend_reports <- function(prediction, days_future) {
  f_extend <- function(effects, sd) {
    if(is.null(effects$report_days))
      return(effects)
    days <- sort(unique(c(effects$report_days, days_future)))
    effects$report_days <- days
    if(!is.null(effects$weekday_log_effects)) {
      weekday_first <- ifelse(is.null(effects$weekday_first), 1, effects$weekday_first)
      missing <- setdiff(as.character(weekday_of_days(days, weekday_first)),
                         names(effects$weekday_log_effects))
      for(w in missing)
        effects$weekday_log_effects[[w]] <- stats::rnorm(1, 0, sd)
      effects <- update_report_multipliers(effects)
    }
    effects
  }
  prediction$revision_parameters <- f_extend(prediction$revision_parameters,
                                             prediction$weekday_sd[["cases"]])
  state <- prediction$deaths_state
  effects <- state$recording_effects
  backlog <- effects$backlog
  effects$backlog <- NULL
  effects <- f_extend(effects, prediction$weekday_sd[["recording"]])
  if(!is.null(backlog)) {
    if(!is.null(backlog$rate_days))
      stop("Predictions with a backlog hazard held over blocks of days are not yet allowed for.")
    new_days <- sort(setdiff(days_future, backlog$report_days))
    new_days <- new_days[new_days > max(backlog$report_days)]
    if(length(new_days) > 0) {
      gaps <- diff(c(max(backlog$report_days), new_days))
      steps <- stats::rnorm(length(new_days), 0, backlog$rate_rw_sd * sqrt(gaps / 7))
      backlog$log_rate <- c(backlog$log_rate,
                            backlog$log_rate[length(backlog$log_rate)] + cumsum(steps))
      backlog$report_days <- c(backlog$report_days, new_days)
      backlog <- refresh_backlog(backlog)
    }
    effects$backlog <- backlog
  }
  state$recording_effects <- effects
  prediction$deaths_state <- state
  prediction
}

#' Simulates the later reports of one onset time from the last one
#'
#' The reports go on as the model has them: of the cases not yet reported,
#' each kind is reported in each later report with its probability; a fatal
#' case dated by onset which is newly reported is shown as dead with
#' probability g, and one shown as alive is later shown as dead with
#' probability q; a fatal case dated by death is shown as dead at once; and
#' wrong records are removed and entered as between earlier reports. Since
#' the model allows at most max_erroneous wrong records at once, later
#' reports are simulated again until none has more.
#'
#' @param path the split of the reports so far (as returned by
#' \code{predict_onset_path})
#' @param probabilities as returned by \code{path_probabilities}, with an
#' element 'p_dated', for the reports so far followed by the later ones
#' @param wrong as in \code{joint_path_forward} for the same reports, or NULL
#' @param totals a named list of 'survivors', 'fatal' (dated by onset) and
#' 'dated' (dated by death): the cases recorded for the day of each kind
#' @param n_before the number of reports so far
#' @param max_tries the most times the later reports are simulated
#'
#' @return a matrix with a row for each later report and columns
#' 'cases_reported' and 'deaths_reported'
predict_onset_simulate <- function(path, probabilities, wrong, totals, n_before,
                                   max_tries=1000) {
  later <- seq(n_before + 1, length(probabilities$p))
  f_last <- function(x) x[n_before]
  # what the last report holds: the cases of each kind reported, the fatal
  # cases dated by onset shown as dead and the wrong records of each kind
  survivors <- f_last(path$survivors_before + path$survivors_new)
  fatal <- f_last(path$fatal_reported)
  dated <- f_last(path$dated_reported)
  alive <- fatal - sum(path$deaths_new)
  is_wrong <- !is.null(wrong)
  max_erroneous <- if(is_wrong) wrong$max_erroneous else 0
  out <- matrix(NA_real_, length(later), 2,
                dimnames=list(NULL, c("cases_reported", "deaths_reported")))
  for(try in seq_len(max_tries)) {
    s <- survivors; f <- fatal; d <- dated; a <- alive
    wa <- f_last(path$wrong_alive); wd <- f_last(path$wrong_dead)
    is_over <- FALSE
    for(i in seq_along(later)) {
      j <- later[i]
      s <- s + stats::rbinom(1, totals$survivors - s, probabilities$p[j])
      d <- d + stats::rbinom(1, totals$dated - d, probabilities$p_dated[j])
      f_new <- stats::rbinom(1, totals$fatal - f, probabilities$p_fatal[j])
      a_dead <- stats::rbinom(1, a, probabilities$q[j])
      f_dead <- stats::rbinom(1, f_new, probabilities$g[j])
      f <- f + f_new
      a <- a + f_new - f_dead - a_dead
      if(is_wrong) {
        wa <- wa - stats::rbinom(1, wa, wrong$w[j]) +
          stats::rpois(1, wrong$lambda[j] * (1 - wrong$delta))
        wd <- wd - stats::rbinom(1, wd, wrong$w[j]) +
          stats::rpois(1, wrong$lambda[j] * wrong$delta)
        if(wa + wd > max_erroneous) {
          is_over <- TRUE
          break
        }
      }
      out[i, ] <- c(s + f + d + wa + wd, f - a + d + wd)
    }
    if(!is_over)
      return(out)
  }
  stop("The later reports kept having more wrong records than max_erroneous.")
}

#' Simulates later reports from a kept draw of a fit of the joint model of
#' cases and deaths
#'
#' For each onset time, the split of its reports (how many of the cases are
#' fatal, of each kind, and how many are wrong records) is drawn given the
#' draw, as the sampler draws it, and the reports on days_future are then
#' simulated from the last report. The effects of reports not yet made are
#' drawn from their priors, and the backlog's hazard carries on its random
#' walk. Given a draw from the posterior, the result is a draw from the
#' posterior predictive distribution of the later reports.
#'
#' Only the exact form of fatal cases dated by death is allowed for, and only
#' onset times already in the reports are simulated.
#'
#' @inheritParams predict_state_from_draw
#' @param days_future the days of the later reports (those not after the last
#' report are left out)
#' @param onset_times the onset times to simulate (by default all of them)
#'
#' @return a tibble with columns time_onset, time_reported, cases_reported
#' and deaths_reported
#' @export
predict_reports <- function(variables, draw, days_future, onset_times=NULL) {
  prediction <- predict_state_from_draw(variables, draw)
  state <- prediction$deaths_state
  if(!isTRUE(state$dated$exact))
    stop("Predictions need the exact form of fatal cases dated by death.")
  if(!is.null(state$uncorrected))
    stop("Predictions with reports made before a correction are not yet allowed for.")
  day_last <- max(unlist(lapply(state$reports, function(r) r$time_reported)))
  days_future <- sort(unique(days_future[days_future > day_last]))
  if(length(days_future) == 0)
    stop("There are no later reports to simulate.")
  # the split of the reports of each onset time depends on the rest alone,
  # so those of the onset times wanted can be drawn by themselves
  wanted <- if(is.null(onset_times)) seq_along(state$onset_times)
    else which(state$onset_times %in% onset_times)
  extended <- predict_extend_reports(prediction, days_future)
  dated <- dated_counts(state)
  out <- lapply(wanted, function(k) {
    path <- predict_onset_path(prediction, k)
    days_before <- state$reports[[k]]$time_reported
    days <- c(days_before, days_future)
    terms <- predict_onset_terms(extended, k, days)
    probabilities <- path_probabilities(terms$survivals, 0)
    probabilities$p_dated <- survival_steps(terms$survivals$log_dated)
    totals <- list(survivors=prediction$cases_true[k] - state$deaths_true[k],
                   fatal=state$deaths_true[k] - dated$out[k],
                   dated=dated$into[k])
    simulated <- predict_onset_simulate(path, probabilities, terms$wrong,
                                        totals, length(days_before))
    dplyr::tibble(time_onset=state$onset_times[k], time_reported=days_future,
                  cases_reported=simulated[, "cases_reported"],
                  deaths_reported=simulated[, "deaths_reported"])
  })
  dplyr::bind_rows(out)
}
