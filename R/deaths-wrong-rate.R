# The rate at which wrong records are entered, in the joint model of cases
# and deaths. In the main model, wrong records of the cases recorded for a day
# are entered in proportion to the genuine cases being reported: on day k the
# number entered is Poisson with mean eta times the genuine cases of the day
# expected to be reported on day k (see wrong_entry_units). Each is withdrawn
# with probability rho on each later day, so those entered since the report
# before which are still present at report j are Poisson with mean
# eta Itilde x_j, where Itilde is the genuine cases recorded for the day.
#
# The setting priors$revisions$wrong_rate of mcmc_single chooses another form:
#
# - "cases" (the default): the main model.
# - "constant": wrong records arrive at a constant rate eta_bg for each day
#   recorded for and each calendar day from the day after it, whatever its
#   cases. Those entered since report j - 1 still present at report j are
#   Poisson with mean eta_bg chi_j, where
#     chi_j = sum_{k = t_{j-1} + 1}^{t_j} (1 - rho)^(t_j - k)
#           = [1 - (1 - rho)^(t_j - t_{j-1})] / rho,
#   with t_0 the day recorded for. The revision parameters then hold
#   'background_rate', eta_bg, in place of 'eta' (a name which 'eta' does not
#   partially match).
# - "walk": as in the main model, but on day k the rate is eta_{J(k)}, where
#   J(k) is the first report, of those of every day recorded for, made on or
#   after day k (every day up to the first report has the first rate, and
#   days after the last report the last). The mean of those entered since
#   report j - 1 still present at report j is then
#     Itilde sum_J eta_J x_{j, J},
#   where x_{j, J} is the part of x_j from the days of interval J. The log of
#   the rates is a random walk over the intervals, with steps of variance
#   tau_eta^2 (days in the interval) / 7. The revision parameters then hold
#   'wrong_walk' in place of 'eta', a named list of:
#   - 'log_rate': the log of the rate of each interval;
#   - 'report_days': the day of the report which ends each interval;
#   - 'tau': tau_eta, the standard deviation of the steps over seven days;
#   - 'mu' and 'sigma': the normal prior of the log of the first rate;
#   - 'tau_scale': the scale of the half-normal prior of tau_eta;
#   - 'tau_step': the standard deviation of the changes proposed in its log.
#
# In every form the withdrawal of wrong records, the share shown as dead and
# the most present at once are as in the main model.

#' Settings of the rate at which wrong records are entered
#'
#' @param prior_parameters the element 'revisions' of the priors of
#' \code{mcmc_single}, which can contain 'wrong_rate': "cases" (the default),
#' "constant" or "walk" (see the start of deaths-wrong-rate.R). With
#' "constant" it must contain 'eta_bg_mu' and 'eta_bg_sigma', the mean and sd
#' of the gamma prior of the rate. With "walk" it can contain 'eta_walk_mu'
#' and 'eta_walk_sigma' (the normal prior of the log of the first rate;
#' log(0.07) and 1.2 by default) and 'tau_eta_scale' (the scale of the
#' half-normal prior of the sd of the steps of the walk over seven days; 0.5
#' by default).
#' @param is_deaths whether there is a model of deaths
#'
#' @return a named list with element 'form' and the prior of the rate
wrong_rate_settings <- function(prior_parameters, is_deaths=TRUE) {
  form <- prior_parameters$wrong_rate
  if(is.null(form))
    return(list(form="cases"))
  if(!is.character(form) || length(form) != 1 ||
     !form %in% c("cases", "constant", "walk"))
    stop("priors$revisions$wrong_rate must be one of 'cases', 'constant' or 'walk' (a constant rate and a walk cannot be combined).")
  if(form == "cases")
    return(list(form="cases"))
  if(!is_deaths)
    stop("A rate of wrong records other than 'cases' needs a model of deaths.")
  if(form == "constant") {
    prior <- prior_parameters[c("eta_bg_mu", "eta_bg_sigma")]
    if(any(vapply(prior, is.null, logical(1))))
      stop("With a constant rate of wrong records, priors$revisions must contain 'eta_bg_mu' and 'eta_bg_sigma'.")
    if(any(unlist(prior) <= 0))
      stop("The mean and sd of the prior of the constant rate of wrong records must be positive.")
    return(c(list(form="constant"), prior))
  }
  prior <- list(eta_walk_mu=log(0.07), eta_walk_sigma=1.2, tau_eta_scale=0.5)
  names_given <- intersect(names(prior), names(prior_parameters))
  prior[names_given] <- prior_parameters[names_given]
  if(any(lengths(prior) != 1))
    stop("The prior of the walk of the rate of wrong records must have single values.")
  if(!is.finite(prior$eta_walk_mu))
    stop("The mean of the prior of the log of the first rate of wrong records must be finite.")
  if(prior$eta_walk_sigma <= 0 | prior$tau_eta_scale <= 0)
    stop("The sd of the prior of the log of the first rate of wrong records, and the scale of the prior of the sd of its steps, must be positive.")
  c(list(form="walk"), prior)
}

#' Starting values of the rate at which wrong records are entered, where it
#' is not that of the main model
#'
#' @param revision_parameters the starting revision parameters of
#' \code{mcmc_single}. With a constant rate these can contain 'eta_bg' (by
#' default the mean of its prior); with a walk, 'eta_walk' (one rate, or one
#' for each interval between reports; by default 'eta') and 'tau_eta' (by
#' default 0.2).
#' @param settings as returned by \code{wrong_rate_settings}
#' @param report_days the days of the reports of every day recorded for,
#' including those made before a correction (for a walk)
#' @param tau_step the standard deviation of the changes proposed in the log
#' of tau_eta (for a walk)
#'
#' @return the revision parameters, with 'background_rate' or 'wrong_walk' in
#' place of 'eta'
initial_wrong_rate <- function(revision_parameters, settings, report_days=NULL,
                               tau_step=0.3) {
  rho <- revision_parameters$rho
  if(is.null(rho) || rho <= 0 | rho >= 1)
    stop("Initial revision parameters must satisfy: 0 < rho < 1.")
  if(settings$form == "constant") {
    eta_bg <- revision_parameters$eta_bg
    if(is.null(eta_bg))
      eta_bg <- settings$eta_bg_mu
    if(length(eta_bg) != 1 || eta_bg <= 0)
      stop("The initial constant rate of wrong records must be a single positive value.")
    revision_parameters[c("eta", "eta_bg")] <- NULL
    revision_parameters$background_rate <- eta_bg
    return(revision_parameters)
  }
  report_days <- sort(unique(report_days))
  initial <- revision_parameters[["eta_walk"]]
  if(is.null(initial))
    initial <- revision_parameters[["eta"]]
  if(!length(initial) %in% c(1, length(report_days)))
    stop(paste0("The initial rates of wrong records must be one value or one for each of the ",
                length(report_days), " intervals between reports."))
  if(any(initial <= 0))
    stop("The initial rates of wrong records must be positive.")
  tau <- revision_parameters[["tau_eta"]]
  if(is.null(tau))
    tau <- 0.2
  if(length(tau) != 1 || tau <= 0 | tau_step <= 0)
    stop("The initial sd of the steps of the walk of the rate of wrong records, and the step proposed in its log, must be positive.")
  revision_parameters[c("eta", "eta_walk", "tau_eta")] <- NULL
  revision_parameters$wrong_walk <- list(
    log_rate=log(rep_len(initial, length(report_days))),
    report_days=report_days, tau=tau, mu=settings$eta_walk_mu,
    sigma=settings$eta_walk_sigma, tau_scale=settings$tau_eta_scale,
    tau_step=tau_step)
  revision_parameters
}

#' The form of the rate at which wrong records are entered
#'
#' @inheritParams log_undetected_days
#'
#' @return "cases", "constant" or "walk" (see the start of deaths-wrong-rate.R)
wrong_rate_form <- function(revision_parameters) {
  if(!is.null(revision_parameters$background_rate))
    return("constant")
  if(!is.null(revision_parameters$wrong_walk))
    return("walk")
  "cases"
}

#' The interval between reports into which each of a set of days falls
#'
#' @param days whole numbers of days
#' @param report_days the days of the reports, of every day recorded for, in
#' increasing order
#'
#' @return for each day, the index of the first report made on or after it;
#' days after the last report have the index of the last
wrong_walk_interval <- function(days, report_days) {
  pmin(findInterval(days - 1e-8, report_days) + 1, length(report_days))
}

#' What governs the wrong records of the cases recorded for one day, where
#' they arrive at a constant rate
#'
#' @inheritParams joint_survivals
#' @param rho the daily probability a wrong record is removed
#'
#' @return as \code{wrong_entry_units}, where 'unit' is for each unit of the
#' rate (not of the rate times the genuine cases): chi_j, the expected number
#' entered since the report before (or, for the first report, since the day
#' recorded for) which are still present at the report, when one is entered
#' on each day
wrong_constant_units <- function(days, day_onset, rho) {
  gaps <- round(diff(c(day_onset, days)))
  list(unit=-expm1(gaps * log1p(-rho)) / rho, w=c(0, 1 - (1 - rho)^diff(days)))
}

#' What governs the wrong records of the cases recorded for one day, where
#' the rate varies over the intervals between reports
#'
#' @inheritParams wrong_entry_units
#'
#' @return as \code{wrong_entry_units}, with an element 'split': a matrix with
#' a row for each report and a column for each interval between reports,
#' holding the part of 'unit' from the days of each interval. 'unit' is the
#' sum of each row.
wrong_walk_units <- function(days, day_onset, reporting_parameters,
                             revision_parameters) {
  rho <- revision_parameters$rho
  report_days <- revision_parameters$wrong_walk$report_days
  days_all <- seq(day_onset, max(day_onset, days[length(days)]))
  undetected <- exp(log_undetected_days(days_all, day_onset, reporting_parameters,
                                        revision_parameters))
  newly <- -diff(undetected)
  # the interval of each day after the onset day
  interval <- wrong_walk_interval(day_onset + seq_along(newly), report_days)
  days_before <- c(day_onset, days[-length(days)])
  split <- matrix(0, length(days), length(report_days))
  for(j in seq_along(days)) {
    k <- seq_len(round(days[j] - days_before[j])) + round(days_before[j] - day_onset)
    if(length(k) == 0)
      next
    sums <- rowsum(newly[k] * (1 - rho)^(round(days[j] - day_onset) - k), interval[k])
    split[j, as.integer(rownames(sums))] <- sums[, 1]
  }
  list(unit=rowSums(split), split=split, w=c(0, 1 - (1 - rho)^diff(days)))
}

#' What governs the wrong records of the cases recorded for one day, in the
#' form of the rate in use
#'
#' @inheritParams wrong_entry_units
#'
#' @return as \code{wrong_entry_units} for the main model,
#' \code{wrong_constant_units} for a constant rate and \code{wrong_walk_units}
#' for a walk
deaths_wrong_units <- function(days, day_onset, reporting_parameters,
                               revision_parameters) {
  form <- wrong_rate_form(revision_parameters)
  if(form == "constant")
    return(wrong_constant_units(days, day_onset, revision_parameters$rho))
  if(form == "walk")
    return(wrong_walk_units(days, day_onset, reporting_parameters,
                            revision_parameters))
  wrong_entry_units(days, day_onset, reporting_parameters, revision_parameters)
}

#' Expected wrong records entered since the report before which are still
#' present at each report
#'
#' @param units as returned by \code{deaths_wrong_units}
#' @param cases_recorded the genuine cases recorded for the day
#' @inheritParams log_undetected_days
#'
#' @return a vector with an element for each report (lambda of
#' \code{joint_path_forward})
wrong_entry_mean <- function(units, cases_recorded, revision_parameters) {
  form <- wrong_rate_form(revision_parameters)
  if(form == "constant")
    return(revision_parameters$background_rate * units$unit)
  if(form == "walk")
    return(cases_recorded *
             as.vector(units$split %*% exp(revision_parameters$wrong_walk$log_rate)))
  revision_parameters$eta * cases_recorded * units$unit
}

#' The wrong records of every series of reports, for the updates of what
#' governs them
#'
#' The uncorrected series are further series of reports with wrong records.
#'
#' @inheritParams deaths_paths_step
#'
#' @return a named list of 'paths', 'reports', 'onset_times', 'wrong_units',
#' 'cases_recorded' (the genuine cases recorded for the day) and
#' 'reporting_by_onset', with an element (or row) for each onset time
#' followed by each uncorrected series
deaths_wrong_series <- function(state, cases_true, reporting_by_onset) {
  dated <- dated_counts(state)
  sub <- uncorrected_state(state)
  list(
    paths=c(state$paths, sub$paths),
    reports=c(state$reports, sub$reports),
    onset_times=c(state$onset_times, sub$onset_times),
    wrong_units=c(state$wrong_units, sub$wrong_units),
    cases_recorded=c(cases_true - dated$out + dated$into, sub$cases_true),
    reporting_by_onset=rbind(
      as.data.frame(reporting_by_onset[, c("mean", "sd")]),
      as.data.frame(reporting_by_onset[sub$day, c("mean", "sd")])))
}

#' The log of the conditional density of rho, the daily probability a wrong
#' record is removed, up to a constant, on the logit scale
#'
#' rho acts on the removal of wrong records present at a report and on how
#' many of those entered since the report before are still present.
#'
#' @param rho a value of rho
#' @param series as returned by \code{deaths_wrong_series}
#' @inheritParams deaths_wrong_step
#'
#' @return a log-density
deaths_wrong_rho_logp <- function(rho, series, revision_parameters,
                                  prior_parameters) {
  proposed <- revision_parameters
  proposed$rho <- rho
  logl <- sum(purrr::map_dbl(seq_along(series$onset_times), function(k) {
    path <- series$paths[[k]]
    units <- deaths_wrong_units(series$reports[[k]]$time_reported,
                                series$onset_times[k],
                                delay_of_onset(series$reporting_by_onset, k), proposed)
    sum(dbinom_log_recycled(path$wrong_removed, path$wrong_before, units$w)) +
      sum(stats::dpois(path$wrong_entered_alive + path$wrong_entered_dead,
                       wrong_entry_mean(units, series$cases_recorded[k], proposed),
                       log=TRUE))
  }))
  logl + stats::dbeta(rho, prior_parameters$rho_alpha,
                      prior_parameters$rho_beta, log=TRUE) +
    log(rho) + log1p(-rho)
}

#' Draws the constant rate at which wrong records arrive from its gamma
#' conditional distribution
#'
#' Given the paths, the wrong records entered since each report before and
#' still present are Poisson with means eta_bg chi_j, so with a gamma prior of
#' shape k and rate theta the conditional distribution is gamma with shape
#' k + n and rate theta + X, where n is the number entered and X the sum of
#' the chi_j.
#'
#' @param series as returned by \code{deaths_wrong_series}
#' @param prior_parameters named list of 'eta_bg_mu' and 'eta_bg_sigma', the
#' mean and sd of the gamma prior
#'
#' @return a draw of eta_bg
wrong_constant_rate_draw <- function(series, prior_parameters) {
  entered <- sum(purrr::map_dbl(series$paths, ~sum(.$wrong_entered_alive +
                                                     .$wrong_entered_dead)))
  exposure <- sum(purrr::map_dbl(series$wrong_units, ~sum(.$unit)))
  stats::rgamma(
    1,
    shape=prior_parameters$eta_bg_mu^2 / prior_parameters$eta_bg_sigma^2 + entered,
    rate=prior_parameters$eta_bg_mu / prior_parameters$eta_bg_sigma^2 + exposure)
}

#' The wrong records entered, and what each interval's rate contributes to
#' their means, for the update of the walk
#'
#' @param series as returned by \code{deaths_wrong_series}, with units from
#' \code{wrong_walk_units}
#'
#' @return a named list of 'entered' (the wrong records entered since the
#' report before and still present, for each report of each series which can
#' have any) and 'exposure' (a matrix with a row for each of these and a
#' column for each interval: the genuine cases recorded times the split
#' units), so that the means are exposure times the rates
wrong_walk_data <- function(series) {
  entered <- unlist(lapply(series$paths, function(path)
    path$wrong_entered_alive + path$wrong_entered_dead))
  exposure <- do.call(rbind, lapply(seq_along(series$paths), function(k)
    series$cases_recorded[k] * series$wrong_units[[k]]$split))
  is_used <- rowSums(exposure) > 0 | entered > 0
  list(entered=entered[is_used], exposure=exposure[is_used, , drop=FALSE])
}

#' Log-likelihood of the log rates of the walk given the wrong records
#' entered
#'
#' @param log_rate the log of the rate of each interval
#' @param data as returned by \code{wrong_walk_data}
#'
#' @return a log-likelihood
wrong_walk_logl <- function(log_rate, data) {
  sum(stats::dpois(data$entered, as.vector(data$exposure %*% exp(log_rate)),
                   log=TRUE))
}

#' Log-prior of the walk of the log rates of wrong records
#'
#' @param log_rate the log of the rate of each interval
#' @param tau the standard deviation of the steps over seven days
#' @param walk the element 'wrong_walk' of the revision parameters
#'
#' @return a log-density
wrong_walk_log_prior <- function(log_rate, tau, walk) {
  gaps <- diff(walk$report_days)
  stats::dnorm(log_rate[1], walk$mu, walk$sigma, log=TRUE) +
    sum(stats::dnorm(diff(log_rate), 0, tau * sqrt(gaps / 7), log=TRUE))
}

#' Draws the walk of the log rates of wrong records by elliptical slice
#' sampling, given the standard deviation of its steps
#'
#' @param walk the element 'wrong_walk' of the revision parameters
#' @param f_logl the log-likelihood of the log rates
#' @param logl the log-likelihood of the current log rates
#'
#' @return a list of the walk updated ('walk') and its log-likelihood ('logl')
wrong_walk_slice <- function(walk, f_logl, logl) {
  n <- length(walk$log_rate)
  sd_steps <- walk$tau * sqrt(diff(walk$report_days) / 7)
  out <- elliptical_slice(
    walk$log_rate, rep(walk$mu, n),
    function() cumsum(c(stats::rnorm(1, 0, walk$sigma),
                        stats::rnorm(n - 1, 0, sd_steps))),
    f_logl, logl)
  walk$log_rate <- out$x
  list(walk=walk, logl=out$logl)
}

#' Updates the walk of the rate of wrong records and the standard deviation
#' of its steps
#'
#' Given the paths, the wrong records entered since each report before and
#' still present are Poisson with means Itilde sum_J eta_J x_{j, J}. The walk
#' of log eta_J is drawn whole by elliptical slice sampling given its
#' standard deviation, which is then changed by Metropolis steps on the log
#' scale, both with the walk held (which needs no likelihood) and with its
#' steps scaled with it (non-centred), as for the backlog
#' (\code{backlog_slice_steps}). The standard deviation has a half-normal
#' prior.
#'
#' @param walk the element 'wrong_walk' of the revision parameters
#' @param data as returned by \code{wrong_walk_data}
#' @param n_steps the number of times each update is made
#'
#' @return the walk updated
wrong_walk_step <- function(walk, data, n_steps=1) {
  f_logl <- function(log_rate) wrong_walk_logl(log_rate, data)
  f_accept <- function(log_ratio)
    !is.na(log_ratio) && log(stats::runif(1)) < log_ratio
  f_half_normal <- function(x) stats::dnorm(x, 0, walk$tau_scale, log=TRUE)
  # a walk with its steps scaled by r, from the same first value
  f_scale <- function(x, r) c(x[1], x[1] + cumsum(diff(x) * r))
  logl <- f_logl(walk$log_rate)
  for(step in seq_len(n_steps)) {
    out <- wrong_walk_slice(walk, f_logl, logl)
    walk <- out$walk
    logl <- out$logl
    if(length(walk$log_rate) == 1)
      next
    # its standard deviation with the walk held: no likelihood is needed
    current <- walk$tau
    proposed <- current * exp(stats::rnorm(1, 0, walk$tau_step))
    log_ratio <- wrong_walk_log_prior(walk$log_rate, proposed, walk) -
      wrong_walk_log_prior(walk$log_rate, current, walk) +
      f_half_normal(proposed) - f_half_normal(current) +
      log(proposed) - log(current)
    if(f_accept(log_ratio))
      walk$tau <- proposed
    # and with its steps scaled with it, the standardised steps held
    current <- walk$tau
    proposed <- current * exp(stats::rnorm(1, 0, walk$tau_step))
    log_rate <- f_scale(walk$log_rate, proposed / current)
    logl_proposed <- f_logl(log_rate)
    log_ratio <- logl_proposed - logl +
      f_half_normal(proposed) - f_half_normal(current) +
      log(proposed) - log(current)
    if(f_accept(log_ratio)) {
      walk$tau <- proposed
      walk$log_rate <- log_rate
      logl <- logl_proposed
    }
  }
  walk
}

#' Updates the rate at which wrong records are entered, where it is not that
#' of the main model
#'
#' @param series as returned by \code{deaths_wrong_series}
#' @inheritParams deaths_wrong_step
#' @param n_steps the number of updates of a walk
#'
#' @return the revision parameters updated
deaths_wrong_rate_step <- function(series, revision_parameters, prior_parameters,
                                   n_steps=1) {
  if(wrong_rate_form(revision_parameters) == "constant") {
    revision_parameters$background_rate <- wrong_constant_rate_draw(series,
                                                                    prior_parameters)
    return(revision_parameters)
  }
  revision_parameters$wrong_walk <- wrong_walk_step(
    revision_parameters$wrong_walk, wrong_walk_data(series), n_steps)
  revision_parameters
}

#' The rate at which wrong records are entered, for storage
#'
#' @inheritParams log_undetected_days
#' @param iteration the iteration
#'
#' @return a tibble with one row and columns 'eta' (NA, since the rate is not
#' that of the main model), 'rho', 'pi' (0) and 'iteration', with 'eta_bg'
#' for a constant rate, or 'tau_eta' and 'eta_<day>' (the rate of the
#' interval ending with the report on that day) for a walk
wrong_rate_draws <- function(revision_parameters, iteration) {
  out <- dplyr::tibble(eta=NA_real_, rho=revision_parameters$rho, pi=0,
                       iteration=iteration)
  if(wrong_rate_form(revision_parameters) == "constant") {
    out$eta_bg <- revision_parameters$background_rate
    return(out)
  }
  walk <- revision_parameters$wrong_walk
  out$tau_eta <- walk$tau
  for(r in seq_along(walk$report_days))
    out[[paste0("eta_", walk$report_days[r])]] <- exp(walk$log_rate[r])
  out
}

#' The rate at which wrong records are entered at a kept draw
#'
#' @inheritParams log_undetected_days
#' @param revisions the draws of the revision parameters (as returned by
#' \code{mcmc_single}) for one iteration
#'
#' @return the revision parameters with the rate of the draw (and, for a walk,
#' the sd of its steps, if the draw has it)
wrong_rate_from_draw <- function(revision_parameters, revisions) {
  form <- wrong_rate_form(revision_parameters)
  if(form == "cases")
    revision_parameters$eta <- revisions$eta
  else if(form == "constant")
    revision_parameters$background_rate <- revisions$eta_bg
  else {
    walk <- revision_parameters$wrong_walk
    walk$log_rate <- log(vapply(walk$report_days, function(day)
      revisions[[paste0("eta_", day)]], numeric(1)))
    if(!is.null(revisions$tau_eta))
      walk$tau <- revisions$tau_eta
    revision_parameters$wrong_walk <- walk
  }
  revision_parameters
}

#' Probability that more wrong records are present than are allowed for at
#' some report of one series of reports
#'
#' The wrong records present at report j are those present at the report
#' before which have not been removed, binomial, together with those entered
#' since and still present, Poisson and independent of them. This follows
#' their distribution from report to report over the numbers allowed for.
#'
#' @param lambda the expected number of wrong records entered since the
#' report before which are still present, for each report
#' @param w the probability that one present at the report before has been
#' removed, for each report
#' @param max_erroneous the most allowed for at once
#'
#' @return the probability that more than max_erroneous are present at one or
#' more of the reports
wrong_cap_probability <- function(lambda, w, max_erroneous) {
  states <- 0:max_erroneous
  p <- c(1, rep(0, max_erroneous))
  for(j in seq_along(lambda)) {
    transition <- matrix(0, length(states), length(states))
    for(from in states)
      for(kept in 0:from)
        transition[from + 1, (kept:max_erroneous) + 1] <-
          transition[from + 1, (kept:max_erroneous) + 1] +
          stats::dbinom(kept, from, 1 - w[j]) *
          stats::dpois(0:(max_erroneous - kept), lambda[j])
    p <- as.vector(p %*% transition)
  }
  max(0, 1 - sum(p))
}

#' Probability, for each series of reports, that more wrong records are
#' present than are allowed for at some report
#'
#' @inheritParams deaths_paths_step
#'
#' @return a tibble with a row for each onset time, followed by each
#' uncorrected series, and columns 'time_onset', 'is_uncorrected' and
#' 'probability' (see \code{wrong_cap_probability})
deaths_wrong_cap <- function(state, cases_true, reporting_by_onset,
                             revision_parameters, max_erroneous=3) {
  series <- deaths_wrong_series(state, cases_true, reporting_by_onset)
  n_corrected <- length(state$onset_times)
  probability <- vapply(seq_along(series$onset_times), function(k) {
    units <- deaths_wrong_units(series$reports[[k]]$time_reported,
                                series$onset_times[k],
                                delay_of_onset(series$reporting_by_onset, k),
                                revision_parameters)
    wrong_cap_probability(
      wrong_entry_mean(units, series$cases_recorded[k], revision_parameters),
      units$w, max_erroneous)
  }, numeric(1))
  dplyr::tibble(time_onset=series$onset_times,
                is_uncorrected=seq_along(series$onset_times) > n_corrected,
                probability=probability)
}

#' Probability that more wrong records are present than are allowed for, at
#' the posterior mean of what governs them
#'
#' With a rate other than that of the main model, the model allows at most
#' max_erroneous wrong records to be present at once in the reports of a day
#' recorded for. This gives, at the posterior mean of the rate (or of each
#' rate of the walk) and of rho, the probability under the model that more
#' are present at some report, for each day recorded for (and each
#' uncorrected series), with the true cases, delays and reports of the state
#' given. A warning is given if it exceeds one per cent for any of them,
#' since the limit is then not negligible.
#'
#' @param variables the variables of the sampler's state (element 'variables'
#' of the 'state' returned by \code{mcmc_single} with return_state TRUE)
#' @param revisions the draws of the revision parameters (element 'revisions'
#' of the result of \code{mcmc_single}, possibly of several chains or calls)
#'
#' @return a named list of 'largest' (the largest of the probabilities),
#' 'expected' (their sum, the expected number of days recorded for at which
#' the limit is exceeded), the posterior means used ('eta_bg' or 'eta_<day>',
#' and 'rho') and 'by_day' (as returned by \code{deaths_wrong_cap})
#' @export
wrong_records_cap <- function(variables, revisions) {
  revision <- variables$revision_current
  if(!isTRUE(variables$is_deaths) || wrong_rate_form(revision) == "cases")
    stop("The probability of exceeding the limit on wrong records is given for a model of deaths with a rate of wrong records other than 'cases'.")
  means <- colMeans(as.matrix(revisions[, setdiff(
    names(revisions), c("eta", "pi", "iteration", "chain", "tau_eta"))]))
  revision <- wrong_rate_from_draw(revision, as.list(means))
  revision$rho <- means[["rho"]]
  by_day <- deaths_wrong_cap(variables$deaths_state, variables$deaths_cases_true,
                             variables$deaths_running_df, revision,
                             variables$max_erroneous)
  largest <- max(by_day$probability)
  if(largest > 0.01)
    warning(sprintf("At the posterior mean, more than %d wrong records are present at some report with probability %.3f for one day recorded for: max_erroneous may need to be increased.",
                    variables$max_erroneous, largest))
  list(largest=largest, expected=sum(by_day$probability), means=means,
       by_day=by_day)
}
