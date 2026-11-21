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
#
# In every form the withdrawal of wrong records, the share shown as dead and
# the most present at once are as in the main model.

#' Settings of the rate at which wrong records are entered
#'
#' @param prior_parameters the element 'revisions' of the priors of
#' \code{mcmc_single}, which can contain 'wrong_rate': "cases" (the default),
#' or "constant" (see the start of deaths-wrong-rate.R). With "constant" it
#' must contain 'eta_bg_mu' and 'eta_bg_sigma', the mean and sd of the gamma
#' prior of the rate.
#' @param is_deaths whether there is a model of deaths
#'
#' @return a named list with element 'form' and the prior of the rate
wrong_rate_settings <- function(prior_parameters, is_deaths=TRUE) {
  form <- prior_parameters$wrong_rate
  if(is.null(form))
    return(list(form="cases"))
  if(!is.character(form) || length(form) != 1 ||
     !form %in% c("cases", "constant"))
    stop("priors$revisions$wrong_rate must be one of 'cases' or 'constant'.")
  if(form == "cases")
    return(list(form="cases"))
  if(!is_deaths)
    stop("A rate of wrong records other than 'cases' needs a model of deaths.")
  prior <- prior_parameters[c("eta_bg_mu", "eta_bg_sigma")]
  if(any(vapply(prior, is.null, logical(1))))
    stop("With a constant rate of wrong records, priors$revisions must contain 'eta_bg_mu' and 'eta_bg_sigma'.")
  if(any(unlist(prior) <= 0))
    stop("The mean and sd of the prior of the constant rate of wrong records must be positive.")
  c(list(form="constant"), prior)
}

#' Starting values of the rate at which wrong records are entered, where it
#' is not that of the main model
#'
#' @param revision_parameters the starting revision parameters of
#' \code{mcmc_single}, which can contain 'eta_bg' (by default the mean of its
#' prior)
#' @param settings as returned by \code{wrong_rate_settings}
#'
#' @return the revision parameters, with 'background_rate' in place of 'eta'
initial_wrong_rate <- function(revision_parameters, settings) {
  rho <- revision_parameters$rho
  if(is.null(rho) || rho <= 0 | rho >= 1)
    stop("Initial revision parameters must satisfy: 0 < rho < 1.")
  eta_bg <- revision_parameters$eta_bg
  if(is.null(eta_bg))
    eta_bg <- settings$eta_bg_mu
  if(length(eta_bg) != 1 || eta_bg <= 0)
    stop("The initial constant rate of wrong records must be a single positive value.")
  revision_parameters[c("eta", "eta_bg")] <- NULL
  revision_parameters$background_rate <- eta_bg
  revision_parameters
}

#' The form of the rate at which wrong records are entered
#'
#' @inheritParams log_undetected_days
#'
#' @return "cases" or "constant" (see the start of deaths-wrong-rate.R)
wrong_rate_form <- function(revision_parameters) {
  if(!is.null(revision_parameters$background_rate))
    return("constant")
  "cases"
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

#' What governs the wrong records of the cases recorded for one day, in the
#' form of the rate in use
#'
#' @inheritParams wrong_entry_units
#'
#' @return as \code{wrong_entry_units} for the main model and
#' \code{wrong_constant_units} for a constant rate
deaths_wrong_units <- function(days, day_onset, reporting_parameters,
                               revision_parameters) {
  form <- wrong_rate_form(revision_parameters)
  if(form == "constant")
    return(wrong_constant_units(days, day_onset, revision_parameters$rho))
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

#' The rate at which wrong records are entered, for storage
#'
#' @inheritParams log_undetected_days
#' @param iteration the iteration
#'
#' @return a tibble with one row and columns 'eta' (NA, since the rate is not
#' that of the main model), 'rho', 'pi' (0), 'iteration' and 'eta_bg'
wrong_rate_draws <- function(revision_parameters, iteration) {
  out <- dplyr::tibble(eta=NA_real_, rho=revision_parameters$rho, pi=0,
                       iteration=iteration)
  out$eta_bg <- revision_parameters$background_rate
  out
}

#' The rate at which wrong records are entered at a kept draw
#'
#' @inheritParams log_undetected_days
#' @param revisions the draws of the revision parameters (as returned by
#' \code{mcmc_single}) for one iteration
#'
#' @return the revision parameters with the rate of the draw
wrong_rate_from_draw <- function(revision_parameters, revisions) {
  if(wrong_rate_form(revision_parameters) == "constant")
    revision_parameters$background_rate <- revisions$eta_bg
  else
    revision_parameters$eta <- revisions$eta
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
#' recorded for. This gives, at the posterior mean of the rate and of rho,
#' the probability under the model that more are present at some report, for
#' each day recorded for (and each uncorrected series), with the true cases,
#' delays and reports of the state given. A warning is given if it exceeds one per cent for any of them, since
#' the limit is then not negligible.
#'
#' @param variables the variables of the sampler's state (element 'variables'
#' of the 'state' returned by \code{mcmc_single} with return_state TRUE)
#' @param revisions the draws of the revision parameters (element 'revisions'
#' of the result of \code{mcmc_single}, possibly of several chains or calls)
#'
#' @return a named list of 'largest' (the largest of the probabilities),
#' 'expected' (their sum, the expected number of days recorded for at which
#' the limit is exceeded), the posterior means used ('eta_bg' and 'rho') and
#' 'by_day' (as returned by \code{deaths_wrong_cap})
#' @export
wrong_records_cap <- function(variables, revisions) {
  revision <- variables$revision_current
  if(!isTRUE(variables$is_deaths) || wrong_rate_form(revision) == "cases")
    stop("The probability of exceeding the limit on wrong records is given for a model of deaths with a rate of wrong records other than 'cases'.")
  means <- colMeans(as.matrix(revisions[, setdiff(
    names(revisions), c("eta", "pi", "iteration", "chain"))]))
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
