# Deaths among cases, and a case fatality ratio which changes at unknown
# onset times.
#
# Of the true cases arising on an onset time, a number eventually die: the
# true fatal cases. These are binomial given the true cases, with a
# probability (the case fatality ratio) which is constant within segments of
# onset times. Every case is reported after a delay from onset. The death of
# a fatal case is recorded after a second delay from onset, and the case is
# shown as dead in a report once it has been reported and its death recorded.
# A report gives, for each onset time, the cases reported and how many of
# them are shown as dead; how many of the others will die is not observed and
# is sampled.
#
# Optionally some fatal cases are recorded for their day of death, not their
# day of onset (see deaths-dated.R). The reports then count cases by the day
# they are recorded for, while transmission and the case fatality ratio are
# by true day of onset.

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
#' change); 'dated_by_death' (NULL, or what is supplied about fatal cases
#' recorded for their day of death: see \code{dated_by_death_settings});
#' 'fatal_delay' (FALSE if fatal cases are reported as others are; otherwise
#' a named list of 'probability_same', the prior probability that they are,
#' 'sd', the standard deviation of the normal prior for the log of the ratio
#' of the two delays if they are not, 'step', the standard deviation of the
#' changes proposed in it, and 'n_steps'); 'recording_weekday' (whether the
#' delay to a death being recorded has its own effects of the day of the
#' week: by default it does if the reporting delay has them);
#' 'wrong_dead_alpha' and 'wrong_dead_beta' (parameters of the beta prior for
#' the probability a wrong record is shown as dead); 'weekday_sd' and
#' 'weekday_step' (the standard deviation of the prior for the log of each
#' effect of the day of the week, and of the changes proposed in it); and
#' 'n_parameter_steps'.
#' The delay here is that from onset to a death being recorded.
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
    dated_by_death=NULL,
    fatal_delay=list(),
    recording_weekday=NULL,
    wrong_dead_alpha=1, wrong_dead_beta=1,
    weekday_sd=1, weekday_step=0.05,
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
  defaults["dated_by_death"] <- list(
    dated_by_death_settings(n_days, defaults$dated_by_death))
  if(!identical(defaults$fatal_delay, FALSE)) {
    fatal_delay <- list(probability_same=0.5, sd=0.5, step=0.1, n_steps=3)
    unknown <- setdiff(names(defaults$fatal_delay), names(fatal_delay))
    if(length(unknown) > 0)
      stop(paste0("Unknown settings for fatal_delay: ",
                  paste(unknown, collapse=", "), "."))
    fatal_delay[names(defaults$fatal_delay)] <- defaults$fatal_delay
    if(fatal_delay$probability_same <= 0 | fatal_delay$probability_same >= 1)
      stop("The prior probability that fatal cases are reported as others are must lie strictly between 0 and 1.")
    if(fatal_delay$sd <= 0 | fatal_delay$step <= 0)
      stop("The standard deviations for the delay of fatal cases must be positive.")
    defaults$fatal_delay <- fatal_delay
  }
  if(defaults$wrong_dead_alpha <= 0 | defaults$wrong_dead_beta <= 0)
    stop("The parameters of the prior for the probability a wrong record is shown as dead must be positive.")
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


#' Logs of the probabilities that the events of the cases recorded for one
#' day have not happened by each report
#'
#' @param days the times of the reports, in increasing order
#' @param day_onset the day for which the cases are recorded
#' @param reporting_parameters named list of the 'mean' and 'sd' of the delay
#' from onset to a case being reported
#' @param death_parameters named list of the 'mean' and 'sd' of the delay from
#' onset to a death being recorded
#' @inheritParams log_undetected_days
#' @param dated_reporting_parameters named list of the 'mean' and 'sd' of the
#' delay from death to being reported, for fatal cases recorded for their day
#' of death; if NULL (the default) this is the delay of other fatal cases
#' counted from onset
#' @param recording_effects a named list which can contain
#' 'report_multipliers' for the delay to a death being recorded, as
#' revision_parameters can for the reporting delay; or NULL
#' @param nu the log of the ratio of the reporting delay of fatal cases to
#' that of others (both its mean and its sd)
#'
#' @return a named list of vectors with an element for each report:
#' 'log_reported' (a case which is not fatal has not been reported),
#' 'log_reported_fatal' (a fatal case dated by onset has not been reported),
#' 'log_recorded' (the death of a fatal case has not been recorded) and
#' 'log_dated' (a fatal case recorded for its day of death has not been
#' reported); and 'is_dated_same' (whether the last is the same as the second)
joint_survivals <- function(days, day_onset, reporting_parameters,
                            death_parameters, revision_parameters=NULL,
                            dated_reporting_parameters=NULL,
                            recording_effects=NULL, nu=0) {
  f <- function(parameters, effects=revision_parameters)
    log_undetected_days(days, day_onset, parameters, effects)
  log_reported <- f(reporting_parameters)
  log_reported_fatal <- log_reported
  if(nu != 0)
    log_reported_fatal <- f(list(mean=reporting_parameters$mean * exp(nu),
                                 sd=reporting_parameters$sd * exp(nu)))
  is_dated_same <- is.null(dated_reporting_parameters)
  list(log_reported=log_reported,
       log_reported_fatal=log_reported_fatal,
       log_recorded=f(death_parameters, recording_effects),
       log_dated=if(is_dated_same) log_reported_fatal
       else f(dated_reporting_parameters),
       is_dated_same=is_dated_same)
}

#' Probability that an event which had not happened by the report before has
#' happened by each report
#'
#' @param log_s the log of the probability the event has not happened by each
#' report
#'
#' @return a vector of probabilities
survival_steps <- function(log_s) {
  -expm1(pmin(log_s - c(0, log_s[-length(log_s)]), 0))
}

#' Probabilities governing the reports of the cases recorded for one day,
#' where a fraction of the fatal cases are recorded for their day of death
#'
#' Each fatal case recorded for the day is taken to be dated by its death
#' with probability 'fraction', independently of the others. The probability
#' that such a case has not been reported is then a mixture over the two
#' kinds, and a fatal case newly reported is one of those dated by death with
#' a probability which depends on the report. A case dated by death is shown
#' as dead in the first report it appears in; one dated by onset only if its
#' death has been recorded by then.
#'
#' @param survivals as returned by \code{joint_survivals}
#' @param fraction the fraction of the fatal cases recorded for the day which
#' are dated by death: a single value or a vector of possible values
#'
#' @return a named list of: 'p' and 'q' (vectors with an element for each
#' report: the probability a non-fatal case not reported at the report before
#' is reported by this one, and the probability a death not recorded at the
#' report before is recorded by this one); and 'p_fatal' and 'g' (matrices
#' with a row for each of 'fraction' and a column for each report: the
#' probability a fatal case not reported at the report before is reported by
#' this one, and the probability that a fatal case newly reported in this
#' report is shown as dead in it)
mixed_probabilities <- function(survivals, fraction=0) {
  n <- length(fraction)
  n_reports <- length(survivals$log_reported)
  f_rows <- function(x) matrix(x, n, n_reports, byrow=TRUE)
  p <- survival_steps(survivals$log_reported)
  # fatal cases dated by onset can be reported at a different speed
  log_fatal <- survivals$log_reported_fatal
  if(is.null(log_fatal))
    log_fatal <- survivals$log_reported
  p_fatal <- f_rows(survival_steps(log_fatal))
  share <- matrix(fraction, n, n_reports)
  if(!survivals$is_dated_same && any(fraction > 0)) {
    log_mixture <- matrix(
      log_add_exp(outer(log1p(-fraction), log_fatal, "+"),
                  outer(log(fraction), survivals$log_dated, "+")),
      n, n_reports)
    log_mixture_before <- cbind(rep(0, n), log_mixture[, -n_reports, drop=FALSE])
    p_fatal <- -expm1(pmin(log_mixture - log_mixture_before, 0))
    log_dated_new <- c(0, survivals$log_dated[-n_reports]) +
      log(survival_steps(survivals$log_dated))
    share_new <- exp(outer(log(fraction), log_dated_new, "+") -
                       (log_mixture_before + log(p_fatal)))
    # no fatal case of either kind can be newly reported
    is_undefined <- !is.finite(share_new)
    share_new[is_undefined] <- share[is_undefined]
    share <- pmin(share_new, 1)
  }
  list(p=p,
       q=survival_steps(survivals$log_recorded),
       p_fatal=p_fatal,
       g=share + (1 - share) * f_rows(-expm1(survivals$log_recorded)))
}

#' Probabilities governing the reports of the cases recorded for one day in
#' the joint model of cases and deaths
#'
#' @inheritParams joint_survivals
#' @param fraction_dated_by_death the fraction of the fatal cases recorded for
#' the day which are dated by death
#'
#' @return a named list of vectors with an element for each report: 'p' (the
#' probability a non-fatal case not yet reported at the report before is
#' reported by this one), 'p_fatal' (the same for a fatal case), 'g' (the
#' probability that a fatal case newly reported in this report is shown as
#' dead in it: either it is dated by its death, or its death has been
#' recorded by this report) and 'q' (the probability a death not recorded at
#' the report before is recorded by this one). Before the first report nothing
#' has been reported or recorded.
joint_probabilities <- function(days, day_onset, reporting_parameters,
                                death_parameters, revision_parameters=NULL,
                                fraction_dated_by_death=0,
                                dated_reporting_parameters=NULL,
                                recording_effects=NULL, nu=0) {
  path_probabilities(
    joint_survivals(days, day_onset, reporting_parameters, death_parameters,
                    revision_parameters, dated_reporting_parameters,
                    recording_effects, nu),
    fraction_dated_by_death)
}

#' As \code{mixed_probabilities} for a single fraction, with vectors only
#'
#' @inheritParams mixed_probabilities
#'
#' @return a named list of vectors 'p', 'p_fatal', 'g' and 'q'
path_probabilities <- function(survivals, fraction=0) {
  probabilities <- mixed_probabilities(survivals, fraction)
  list(p=probabilities$p,
       p_fatal=as.vector(probabilities$p_fatal),
       g=as.vector(probabilities$g),
       q=probabilities$q)
}

#' Binomial log-probability which is -Inf wherever the count is impossible
#'
#' @param x a count or vector of counts
#' @param size a number of trials or vector of them
#' @param prob a probability or vector of them
#'
#' @return a vector of log-probabilities, as long as the longest argument
dbinom_log_recycled <- function(x, size, prob) {
  n <- max(length(x), length(size), length(prob))
  x <- rep_len(x, n)
  size <- rep_len(size, n)
  prob <- rep_len(prob, n)
  logp <- rep(-Inf, n)
  is_valid <- size >= 0 & x >= 0 & x <= size
  logp[is_valid] <- stats::dbinom(x[is_valid], size[is_valid], prob[is_valid],
                                  log=TRUE)
  logp
}

#' Log of the sum of exponentials of two vectors, element by element
#'
#' @param a a vector
#' @param b a vector
#'
#' @return a vector
log_add_exp <- function(a, b) {
  larger <- pmax(a, b)
  out <- larger + log1p(exp(-abs(a - b)))
  out[!is.finite(larger)] <- larger[!is.finite(larger)]
  out
}

#' Log-probability of the deaths newly shown in a report
#'
#' Deaths newly shown are those of fatal cases newly reported which are shown
#' as dead at once (each with probability g) and those of fatal cases
#' reported earlier and shown as alive whose death has now been recorded (each
#' with probability q). All arguments can be vectors.
#'
#' @param deaths_new the number of deaths newly shown
#' @param fatal_new the number of fatal cases newly reported
#' @param fatal_alive the number of fatal cases reported earlier whose deaths
#' had not been recorded
#' @param g a probability
#' @param q a probability
#'
#' @return a vector of log-probabilities
deaths_new_logp <- function(deaths_new, fatal_new, fatal_alive, g, q) {
  n <- max(length(deaths_new), length(fatal_new), length(fatal_alive),
           length(g), length(q))
  if(n == 0)
    return(numeric(0))
  logp <- rep(-Inf, n)
  for(x in 0:min(max(deaths_new), max(fatal_new)))
    logp <- log_add_exp(
      logp,
      dbinom_log_recycled(x, fatal_new, g) +
        dbinom_log_recycled(deaths_new - x, fatal_alive, q))
  logp
}

#' Samples, for one onset time, how many of the cases in each report are fatal
#'
#' The cases in a report are those which will survive and those which will
#' die, of which only those whose deaths have been recorded are shown as dead.
#' The number of fatal cases among those reported is drawn for every report
#' together, exactly, from its distribution given the true numbers of fatal
#' and non-fatal cases and the reports, by filtering forwards over the reports
#' and sampling backwards.
#'
#' @param cases_reported the cases in each report
#' @param deaths_reported the deaths in each report
#' @param deaths_true the true number of fatal cases
#' @param survivors_true the true number of non-fatal cases
#' @param probabilities as returned by \code{joint_probabilities}; if it has
#' no element 'p_fatal', fatal cases are reported as others are
#'
#' @return a named list of vectors with an element for each report:
#' 'fatal_new' and 'survivors_new' (the fatal and non-fatal cases newly
#' reported), 'fatal_before' and 'survivors_before' (those reported by the
#' report before), 'alive_before' (fatal cases in the report before whose
#' deaths had not been recorded) and 'deaths_new' (the deaths newly shown)
sample_fatal_path <- function(cases_reported, deaths_reported, deaths_true,
                              survivors_true, probabilities) {
  forward <- fatal_path_forward(cases_reported, deaths_reported, deaths_true,
                                survivors_true, probabilities)
  if(!is.finite(forward$logp))
    stop("The reports of cases and deaths are impossible given the true cases and deaths.")
  fatal_path_backward(forward)
}

#' Filters forwards over the reports of one onset time
#'
#' @inheritParams sample_fatal_path
#'
#' @return a named list holding the log-weights of each number of fatal cases
#' reported by each report, the log-weights of the changes between reports,
#' and 'logp': the log-probability of the reports given the true numbers of
#' fatal and non-fatal cases, with the number of fatal cases among those
#' reported summed out
fatal_path_forward <- function(cases_reported, deaths_reported, deaths_true,
                               survivors_true, probabilities) {
  n_reports <- length(cases_reported)
  cases_before <- c(0, cases_reported[-n_reports])
  deaths_before <- c(0, deaths_reported[-n_reports])
  cases_new <- cases_reported - cases_before
  deaths_new <- deaths_reported - deaths_before
  if(any(cases_new < 0) | any(deaths_new < 0))
    stop("Reported cases and deaths must not fall in the joint model of cases and deaths.")
  if(any(deaths_reported > cases_reported))
    stop("Reported deaths must not exceed reported cases.")
  p_fatal <- probabilities$p_fatal
  if(is.null(p_fatal))
    p_fatal <- probabilities$p
  r_max <- min(max(cases_reported), deaths_true)
  r_values <- 0:r_max

  # forwards: log-weights of each number of fatal cases reported
  log_alpha <- matrix(-Inf, n_reports + 1, r_max + 1)
  log_alpha[1, 1] <- 0
  l_log_w <- vector("list", n_reports)
  for(j in seq_len(n_reports)) {
    p <- probabilities$p[j]
    is_valid <- is.finite(log_alpha[j, ]) &
      r_values >= deaths_before[j] &
      r_values <= cases_before[j] &
      cases_before[j] - r_values <= survivors_true
    log_w <- matrix(-Inf, cases_new[j] + 1, r_max + 1)
    if(any(is_valid)) {
      r <- r_values[is_valid]
      for(a in 0:cases_new[j]) {
        log_w[a + 1, is_valid] <-
          dbinom_log_recycled(a, deaths_true - r, p_fatal[j]) +
          dbinom_log_recycled(cases_new[j] - a,
                          survivors_true - (cases_before[j] - r), p) +
          deaths_new_logp(deaths_new[j], a, r - deaths_before[j],
                          probabilities$g[j], probabilities$q[j])
      }
    }
    l_log_w[[j]] <- log_w
    for(a in 0:cases_new[j]) {
      index_to <- which(r_values - a >= 0)
      index_from <- index_to - a
      log_alpha[j + 1, index_to] <- log_add_exp(
        log_alpha[j + 1, index_to],
        log_alpha[j, index_from] + log_w[a + 1, index_from])
    }
    # the numbers reported must be possible given the reports
    is_possible <- r_values >= deaths_reported[j] & r_values <= cases_reported[j] &
      cases_reported[j] - r_values <= survivors_true
    log_alpha[j + 1, !is_possible] <- -Inf
  }
  list(log_alpha=log_alpha, l_log_w=l_log_w, r_values=r_values,
       cases_new=cases_new, cases_before=cases_before,
       deaths_new=deaths_new, deaths_before=deaths_before,
       logp=matrixStats::logSumExp(log_alpha[n_reports + 1, ]))
}

#' Samples backwards over the reports of one onset time
#'
#' @param forward as returned by \code{fatal_path_forward}
#'
#' @return as \code{sample_fatal_path}
fatal_path_backward <- function(forward) {
  log_alpha <- forward$log_alpha
  l_log_w <- forward$l_log_w
  r_values <- forward$r_values
  cases_new <- forward$cases_new
  cases_before <- forward$cases_before
  deaths_new <- forward$deaths_new
  deaths_before <- forward$deaths_before
  n_reports <- length(cases_new)
  f_sample <- function(logps) {
    probs <- exp(logps - max(logps))
    sample.int(length(logps), 1, prob=probs)
  }
  fatal <- rep(0, n_reports + 1)
  fatal[n_reports + 1] <- r_values[f_sample(log_alpha[n_reports + 1, ])]
  for(j in rev(seq_len(n_reports))) {
    a_values <- 0:min(cases_new[j], fatal[j + 1])
    index_from <- fatal[j + 1] - a_values + 1
    logps <- log_alpha[j, index_from] +
      l_log_w[[j]][cbind(a_values + 1, index_from)]
    fatal[j] <- fatal[j + 1] - a_values[f_sample(logps)]
  }
  fatal_before <- fatal[-(n_reports + 1)]
  fatal_new <- diff(fatal)
  list(
    fatal_new=fatal_new,
    survivors_new=cases_new - fatal_new,
    fatal_before=fatal_before,
    survivors_before=cases_before - fatal_before,
    alive_before=fatal_before - deaths_before,
    deaths_new=deaths_new)
}

#' Log-probability of the cases newly reported, of one kind, in each report
#'
#' @param counts_new the cases of the kind newly reported in each report
#' @param counts_before the cases of the kind reported by the report before
#' @param totals a vector of possible true numbers of cases of the kind
#' @param p the probability, for each report, that a case not reported by
#' the report before is reported by it
#'
#' @return a vector with a log-probability for each of totals
reports_stream_logp <- function(counts_new, counts_before, totals, p) {
  logp <- rep(0, length(totals))
  for(j in seq_along(counts_new))
    logp <- logp + dbinom_log_recycled(counts_new[j], totals - counts_before[j], p[j])
  logp
}

#' What governs the wrong records of the cases recorded for one day
#'
#' Wrong records are entered in proportion to the genuine cases being
#' reported and each is later removed with a fixed probability each day, as in
#' \code{observation_process_revisions_logp}.
#'
#' @inheritParams joint_survivals
#' @param revision_parameters named list containing 'rho', the daily
#' probability a wrong record is removed, and anything which alters the
#' reporting delay (see \code{log_undetected_days})
#'
#' @return a named list of vectors with an element for each report: 'unit'
#' (the expected number of wrong records entered since the report before, for
#' each genuine case recorded for the day and each unit of eta; for the first
#' report, of those entered since the onset day which are still present) and
#' 'w' (the probability a wrong record present at the report before has been
#' removed)
wrong_entry_units <- function(days, day_onset, reporting_parameters,
                              revision_parameters) {
  rho <- revision_parameters$rho
  undetected <- exp(log_undetected_days(days, day_onset, reporting_parameters,
                                        revision_parameters))
  unit_first <- 0
  n_days <- floor(days[1] - day_onset)
  if(n_days >= 1) {
    days_before <- day_onset + seq_len(n_days)
    undetected_before <- exp(log_undetected_days(
      c(day_onset, days_before), day_onset, reporting_parameters,
      revision_parameters))
    unit_first <- sum(-diff(undetected_before) * (1 - rho)^(days[1] - days_before))
  }
  list(unit=c(unit_first, -diff(undetected)),
       w=c(0, 1 - (1 - rho)^diff(days)))
}

#' Filters forwards over the reports of the cases recorded for one day,
#' allowing for wrong records
#'
#' @inheritParams sample_fatal_path
#' @param deaths_true the number of fatal cases recorded for the day
#' @param survivors_true the number of other cases recorded for the day
#' @param wrong NULL if there are no wrong records; otherwise a named list of
#' 'lambda' and 'w' (vectors with an element for each report: the expected
#' number of wrong records entered since the report before and the
#' probability one present at the report before has been removed), 'delta'
#' (the probability a wrong record is shown as dead) and 'max_erroneous' (the
#' largest number present at once)
#'
#' @return a named list of 'logp' (the log-probability of the reports, with
#' the fatal cases among those reported and the wrong records summed out),
#' 'alpha' (see \code{joint_path_forward_cpp}) and 'arguments'
joint_path_forward <- function(cases_reported, deaths_reported, deaths_true,
                               survivors_true, probabilities, wrong=NULL) {
  n_reports <- length(cases_reported)
  if(any(deaths_reported > cases_reported))
    stop("Reported deaths must not exceed reported cases.")
  if(is.null(wrong))
    wrong <- list(lambda=rep(0, n_reports), w=rep(0, n_reports), delta=0,
                  max_erroneous=0)
  p_fatal <- probabilities$p_fatal
  if(is.null(p_fatal))
    p_fatal <- probabilities$p
  arguments <- list(
    cases=as.numeric(cases_reported), deaths=as.numeric(deaths_reported),
    fatal=as.integer(deaths_true), survivors=as.integer(survivors_true),
    p=probabilities$p, p_fatal=p_fatal, g=probabilities$g, q=probabilities$q,
    lambda=as.numeric(wrong$lambda), w=as.numeric(wrong$w), delta=wrong$delta,
    max_erroneous=as.integer(wrong$max_erroneous))
  if(deaths_true < 0 | survivors_true < 0)
    return(list(logp=-Inf, alpha=NULL, arguments=arguments))
  forward <- do.call(joint_path_forward_cpp, arguments)
  # states of negligible probability are dropped for speed; if that leaves
  # nothing, none are dropped
  if(!is.finite(forward$logp))
    forward <- do.call(joint_path_forward_cpp, c(arguments, list(negligible=0)))
  forward$arguments <- arguments
  forward
}

#' Samples backwards over the reports of the cases recorded for one day
#'
#' @param forward as returned by \code{joint_path_forward}
#'
#' @return as \code{sample_fatal_path}, where the cases and deaths are the
#' genuine ones, with vectors added for the wrong records: 'wrong_before'
#' (those present at the report before), 'wrong_removed' (those of them
#' removed), and 'wrong_entered_alive' and 'wrong_entered_dead' (those entered
#' which are shown as alive and as dead)
joint_path_backward <- function(forward) {
  arguments <- forward$arguments
  drawn <- do.call(joint_path_backward_cpp,
                   c(list(alpha=forward$alpha), arguments))
  f_before <- function(x) c(0, x[-length(x)])
  fatal <- drawn[, 1]
  wrong_alive <- drawn[, 2]
  wrong_dead <- drawn[, 3]
  wrong <- wrong_alive + wrong_dead
  survivors <- arguments$cases - fatal - wrong
  deaths <- arguments$deaths - wrong_dead
  list(
    fatal_new=fatal - f_before(fatal),
    survivors_new=survivors - f_before(survivors),
    fatal_before=f_before(fatal),
    survivors_before=f_before(survivors),
    alive_before=f_before(fatal) - f_before(deaths),
    deaths_new=deaths - f_before(deaths),
    wrong_before=f_before(wrong),
    wrong_removed=drawn[, 4] + drawn[, 5],
    wrong_entered_alive=wrong_alive - f_before(wrong_alive) + drawn[, 4],
    wrong_entered_dead=wrong_dead - f_before(wrong_dead) + drawn[, 5])
}

#' Filters forwards over the reports of the cases recorded for one day,
#' following the fatal cases dated by onset and by death separately
#'
#' This is the exact form of the model for fatal cases recorded for their day
#' of death: the number of them among the cases recorded for the day is a
#' fixed number, not a probability for each fatal case.
#'
#' @inheritParams joint_path_forward
#' @param fatal_onset the number of fatal cases recorded for the day which
#' are dated by onset
#' @param fatal_dated the number which are dated by death
#' @param survivals as returned by \code{joint_survivals}
#'
#' @return as \code{joint_path_forward}
joint_path_exact_forward <- function(cases_reported, deaths_reported, fatal_onset,
                                     fatal_dated, survivors_true, survivals,
                                     wrong=NULL) {
  n_reports <- length(cases_reported)
  if(any(deaths_reported > cases_reported))
    stop("Reported deaths must not exceed reported cases.")
  if(is.null(wrong))
    wrong <- list(lambda=rep(0, n_reports), w=rep(0, n_reports), delta=0,
                  max_erroneous=0)
  probabilities <- path_probabilities(survivals, 0)
  arguments <- list(
    cases=as.numeric(cases_reported), deaths=as.numeric(deaths_reported),
    fatal=as.integer(fatal_onset), dated=as.integer(fatal_dated),
    survivors=as.integer(survivors_true),
    p=probabilities$p, p_fatal=probabilities$p_fatal,
    p_dated=survival_steps(survivals$log_dated),
    g=probabilities$g, q=probabilities$q,
    lambda=as.numeric(wrong$lambda), w=as.numeric(wrong$w), delta=wrong$delta,
    max_erroneous=as.integer(wrong$max_erroneous))
  if(fatal_onset < 0 | fatal_dated < 0 | survivors_true < 0)
    return(list(logp=-Inf, alpha=NULL, arguments=arguments))
  forward <- do.call(joint_path_exact_forward_cpp, arguments)
  # states of negligible probability are dropped for speed; if that leaves
  # nothing, none are dropped
  if(!is.finite(forward$logp))
    forward <- do.call(joint_path_exact_forward_cpp, c(arguments, list(negligible=0)))
  forward$arguments <- arguments
  forward
}

#' Samples backwards over the reports of the cases recorded for one day,
#' following the fatal cases dated by onset and by death separately
#'
#' @param forward as returned by \code{joint_path_exact_forward}
#'
#' @return as \code{joint_path_backward}, where the fatal cases and the
#' deaths are those dated by onset, with vectors 'dated_new' and
#' 'dated_before' added for the fatal cases dated by death (those newly
#' reported in each report and those reported by the report before), and
#' 'is_exact' (TRUE)
joint_path_exact_backward <- function(forward) {
  arguments <- forward$arguments
  drawn <- do.call(joint_path_exact_backward_cpp,
                   c(list(alpha=forward$alpha), arguments))
  f_before <- function(x) c(0, x[-length(x)])
  fatal <- drawn[, 1]
  dated <- drawn[, 2]
  wrong_alive <- drawn[, 3]
  wrong_dead <- drawn[, 4]
  wrong <- wrong_alive + wrong_dead
  survivors <- arguments$cases - fatal - dated - wrong
  deaths <- arguments$deaths - wrong_dead - dated
  list(
    fatal_new=fatal - f_before(fatal),
    survivors_new=survivors - f_before(survivors),
    fatal_before=f_before(fatal),
    survivors_before=f_before(survivors),
    alive_before=f_before(fatal) - f_before(deaths),
    deaths_new=deaths - f_before(deaths),
    wrong_before=f_before(wrong),
    wrong_removed=drawn[, 5] + drawn[, 6],
    wrong_entered_alive=wrong_alive - f_before(wrong_alive) + drawn[, 5],
    wrong_entered_dead=wrong_dead - f_before(wrong_dead) + drawn[, 6],
    dated_new=dated - f_before(dated),
    dated_before=f_before(dated),
    is_exact=TRUE)
}

#' Log-probability, given how many of the cases in each report are fatal, of
#' the fatal cases newly reported and the deaths newly shown in each report
#'
#' @param path as returned by \code{sample_fatal_path}
#' @param survivals as returned by \code{joint_survivals}
#' @param fatal_recorded a vector of possible numbers of fatal cases recorded
#' for the day
#' @param dated_in the number of them which are dated by death
#' @param is_deaths_included if FALSE the deaths newly shown are left out
#'
#' @return a vector with a log-probability for each of fatal_recorded
fatal_reports_logp <- function(path, survivals, fatal_recorded, dated_in=0,
                               is_deaths_included=TRUE) {
  n <- length(fatal_recorded)
  n_reports <- length(path$fatal_new)
  if(isTRUE(path$is_exact)) {
    # the two kinds of fatal case are followed separately
    probabilities <- path_probabilities(survivals, 0)
    logp <- reports_stream_logp(path$fatal_new, path$fatal_before,
                                fatal_recorded - dated_in, probabilities$p_fatal) +
      sum(dbinom_log_recycled(path$dated_new, dated_in - path$dated_before,
                              survival_steps(survivals$log_dated)))
    is_informative <- path$fatal_new > 0 | path$alive_before > 0
    if(is_deaths_included)
      logp <- logp + sum(deaths_new_logp(
        path$deaths_new[is_informative], path$fatal_new[is_informative],
        path$alive_before[is_informative], probabilities$g[is_informative],
        probabilities$q[is_informative]))
    return(logp)
  }
  fraction <- ifelse(fatal_recorded > 0,
                     pmin(dated_in / pmax(fatal_recorded, 1), 1), 0)
  probabilities <- mixed_probabilities(survivals, fraction)
  logp <- rowSums(matrix(
    dbinom_log_recycled(rep(path$fatal_new, each=n),
                        as.vector(outer(fatal_recorded, path$fatal_before, "-")),
                        as.vector(probabilities$p_fatal)),
    n, n_reports))
  is_informative <- path$fatal_new > 0 | path$alive_before > 0
  if(is_deaths_included & any(is_informative))
    logp <- logp + rowSums(matrix(
      deaths_new_logp(rep(path$deaths_new[is_informative], each=n),
                      rep(path$fatal_new[is_informative], each=n),
                      rep(path$alive_before[is_informative], each=n),
                      as.vector(probabilities$g[, is_informative]),
                      rep(probabilities$q[is_informative], each=n)),
      n, sum(is_informative)))
  logp[fatal_recorded < dated_in] <- -Inf
  logp
}

#' Log-probability of the reports of one onset time for possible true cases
#' and true fatal cases
#'
#' The true cases and true fatal cases are those whose onset was on the day.
#' Of the fatal cases, dated_out are recorded for a later day (their day of
#' death) and so are not among the cases recorded for the day, while dated_in
#' fatal cases with earlier onsets are.
#'
#' @inheritParams fatal_reports_logp
#' @param cfr the case fatality ratio of the onset time
#' @param cases_true a vector of possible true cases
#' @param dated_out the number of the fatal cases of the onset time which are
#' dated by death
#' @param fraction the probability that a fatal case of the onset time is
#' dated by death
#' @param wrong_mean NULL if there are no wrong records; otherwise, for each
#' report, the expected number of wrong records entered since the report
#' before for each genuine case recorded for the day
#'
#' @return a matrix with a row for each of cases_true and a column for each
#' possible number of fatal cases, holding the log of the joint probability
#' of the fatal cases, of dated_out of them being dated by death and of the
#' reports (-Inf where impossible). Its attribute 'deaths_true' holds the
#' numbers of fatal cases to which the columns correspond.
joint_cases_deaths_logp <- function(path, survivals, cfr, cases_true,
                                    dated_out=0, dated_in=0, fraction=0,
                                    wrong_mean=NULL) {
  n_reports <- length(path$fatal_new)
  fatal_reported <- path$fatal_before[n_reports] + path$fatal_new[n_reports]
  survivors_reported <- path$survivors_before[n_reports] +
    path$survivors_new[n_reports]
  # fatal cases of the onset time which are dated by onset: in the exact
  # form the path counts these alone
  onset_min <- max(fatal_reported - ifelse(isTRUE(path$is_exact), 0, dated_in), 0)
  onset_max <- max(max(cases_true) - dated_out - survivors_reported, onset_min)
  fatal_values <- (onset_min:onset_max) + dated_out
  survivors_max <- max(max(cases_true) - fatal_values[1], survivors_reported)
  survivors_values <- survivors_reported:survivors_max
  logp_fatal <- fatal_reports_logp(path, survivals,
                                   fatal_values - dated_out + dated_in, dated_in,
                                   is_deaths_included=dated_in > 0) +
    stats::dbinom(dated_out, fatal_values, fraction, log=TRUE)
  logp_survivors <- reports_stream_logp(path$survivors_new, path$survivors_before,
                                        survivors_values,
                                        survival_steps(survivals$log_reported))
  cases <- matrix(cases_true, length(cases_true), length(fatal_values))
  fatal <- matrix(fatal_values, length(cases_true), length(fatal_values), byrow=TRUE)
  survivors <- cases - fatal
  is_possible <- survivors >= survivors_reported
  logp <- matrix(-Inf, length(cases_true), length(fatal_values))
  logp[is_possible] <- stats::dbinom(fatal[is_possible], cases[is_possible], cfr,
                                     log=TRUE) +
    logp_fatal[fatal[is_possible] - fatal_values[1] + 1] +
    logp_survivors[survivors[is_possible] - survivors_reported + 1]
  # wrong records are entered in proportion to the cases recorded for the day
  if(!is.null(wrong_mean)) {
    entered <- path$wrong_entered_alive + path$wrong_entered_dead
    logp <- logp + purrr::map_dbl(
      cases_true - dated_out + dated_in,
      ~sum(stats::dpois(entered, wrong_mean * max(., 0), log=TRUE)))
  }
  attr(logp, "deaths_true") <- fatal_values
  logp
}

#' Starting state of the model of deaths
#'
#' @param snapshot_df a tibble with columns: time_onset, time_reported,
#' cases_reported and deaths_reported
#' @param cases_true_df a tibble with columns time_onset and cases_true: the
#' cases taken at the start to be recorded for each day
#' @param settings as returned by \code{deaths_settings}
#' @inheritParams sample_cases_history
#'
#' @return a named list holding the reports of each onset time, the true
#' fatal cases and case fatality ratio by onset time, and the segments of the
#' ratio and of the delay from onset to a death being recorded. Its element
#' 'cases_onset' holds the true cases by day of onset with which to start,
#' which differ from those supplied if some fatal cases are dated by death
#' (see \code{initial_dated_state}).
#' @importFrom rlang .data
initial_deaths_state <- function(snapshot_df, cases_true_df, settings,
                                 cases_before=NULL) {
  if(!all(c("time_onset", "time_reported", "cases_reported", "deaths_reported") %in%
          colnames(snapshot_df)))
    stop("The snapshots must contain columns: 'time_onset', 'time_reported', 'cases_reported' and 'deaths_reported'.")
  snapshot <- snapshot_df %>%
    dplyr::select("time_onset", "time_reported", "cases_reported", "deaths_reported") %>%
    dplyr::arrange(.data$time_onset, .data$time_reported)
  onset_times <- sort(unique(snapshot$time_onset))
  cases_true_df <- dplyr::arrange(unique(cases_true_df[, c("time_onset", "cases_true")]),
                                  .data$time_onset)
  if(!identical(as.numeric(onset_times), as.numeric(cases_true_df$time_onset)))
    stop("There must be initial true cases for each onset time.")
  if(any(diff(onset_times) != 1))
    stop("Onset times must be consecutive.")
  if(any(snapshot$deaths_reported > snapshot$cases_reported))
    stop("Reported deaths must not exceed reported cases.")
  n_days <- length(onset_times)
  reports <- split(as.data.frame(snapshot[, c("time_reported", "cases_reported",
                                             "deaths_reported")]),
                   snapshot$time_onset)[as.character(onset_times)]
  largest <- purrr::map_dbl(reports, ~max(.$deaths_reported))
  cases_true <- cases_true_df$cases_true
  if(any(cases_true < purrr::map_dbl(reports, ~max(.$cases_reported))))
    stop("Initial true cases must be at least the cases reported.")
  deaths_true <- pmin(cases_true,
                      pmax(largest, round(settings$initial_cfr * cases_true)))
  state <- list(
    reports=reports,
    onset_times=onset_times,
    deaths_true=unname(deaths_true),
    cases_onset=cases_true,
    cfr=rep(settings$initial_cfr, n_days),
    cfr_starts=1,
    p_change=settings$alpha / (settings$alpha + settings$beta),
    recording=initial_reporting_changepoints(settings$initial_delay),
    recording_effects=NULL,
    nu=0,
    delta=settings$wrong_dead_alpha /
      (settings$wrong_dead_alpha + settings$wrong_dead_beta),
    dated=NULL,
    paths=NULL,
    survivals=NULL,
    wrong_units=NULL)
  if(!is.null(settings$dated_by_death))
    state <- initial_dated_state(state, settings$dated_by_death, cases_before)
  state
}

#' Delay from onset to a death being recorded for each onset time
#'
#' @param state as returned by \code{initial_deaths_state}
#'
#' @return a tibble with columns 'mean' and 'sd' and a row for each onset time
deaths_delay_by_onset <- function(state) {
  n_days <- length(state$onset_times)
  state$recording$parameters[segment_indices(state$recording$starts, n_days), ]
}

#' Fatal cases dated by death, by the day of their onset and by the day for
#' which they are recorded
#'
#' @inheritParams deaths_delay_by_onset
#'
#' @return a named list of vectors with an element for each day: 'out' (the
#' fatal cases with onset on the day which are dated by death, and so are
#' recorded for a later day), 'into' (the fatal cases dated by death which
#' are recorded for the day), 'fraction' (the probability that a fatal case
#' with onset on the day is dated by death) and 'fatal_recorded' (the fatal
#' cases recorded for the day). Without such cases the first three are zero.
dated_counts <- function(state) {
  zeros <- rep(0, length(state$onset_times))
  counts <- list(out=zeros, into=zeros, fraction=zeros)
  if(!is.null(state$dated))
    counts <- state$dated[c("out", "into", "fraction")]
  counts$fatal_recorded <- state$deaths_true - counts$out + counts$into
  counts
}

#' The delay parameters of one onset time, as a named list
#'
#' @param by_onset a tibble with columns 'mean' and 'sd' and a row for each
#' onset time
#' @param k the onset time, counted from one
#'
#' @return a named list of 'mean' and 'sd'
delay_of_onset <- function(by_onset, k) {
  list(mean=by_onset$mean[k], sd=by_onset$sd[k])
}

#' Whether wrong records are allowed for
#'
#' @inheritParams log_undetected_days
#'
#' @return TRUE or FALSE
is_wrong_records <- function(revision_parameters) {
  !is.null(revision_parameters$eta)
}

#' Samples how many of the cases in each report are fatal, and how many are
#' wrong records, for every onset time
#'
#' @inheritParams deaths_delay_by_onset
#' @param cases_true a vector of the true cases of each onset time
#' @param reporting_by_onset a tibble with columns 'mean' and 'sd' and a row
#' for each onset time: the delay from onset to a case being reported
#' @param revision_parameters NULL, or a named list which can contain 'eta'
#' and 'rho' (if wrong records are allowed for: see
#' \code{observation_process_revisions_logp}) and anything which alters the
#' reporting delay (see \code{log_undetected_days})
#' @inheritParams erroneous_transitions
#'
#' @param n_proposals the number of times a change in the true fatal cases
#' of each onset time is proposed
#' @param step_max the largest change proposed
#'
#' @details Once all the cases of an onset time have been reported, the
#' number of fatal cases among those reported equals the true number of fatal
#' cases, so that neither can change while the other is held fixed. The true
#' number of fatal cases is therefore first updated with the split of the
#' reported cases summed out: a change is proposed and accepted with the
#' Metropolis probability, using the probability of the reports from
#' \code{joint_path_forward}. The split is then drawn given the result. The
#' fatal cases dated by death, into and out of the day, are held fixed.
#'
#' @return the state with elements 'deaths_true', 'paths' (see
#' \code{joint_path_backward}), 'survivals' (see \code{joint_survivals}) and
#' 'wrong_units' (see \code{wrong_entry_units}) updated
deaths_paths_step <- function(state, cases_true, reporting_by_onset,
                              revision_parameters=NULL, max_erroneous=3,
                              n_proposals=2, step_max=3) {
  delays <- deaths_delay_by_onset(state)
  dated <- dated_counts(state)
  is_wrong <- is_wrong_records(revision_parameters)
  is_exact <- isTRUE(state$dated$exact)
  n_days <- length(state$onset_times)
  paths <- vector("list", n_days)
  l_survivals <- vector("list", n_days)
  l_units <- vector("list", n_days)
  l_wrong <- vector("list", n_days)
  for(k in seq_len(n_days)) {
    reports <- state$reports[[k]]
    reporting_parameters <- delay_of_onset(reporting_by_onset, k)
    survivals <- joint_survivals(
      reports$time_reported, state$onset_times[k], reporting_parameters,
      delay_of_onset(delays, k), revision_parameters, state$dated$reporting,
      state$recording_effects, state$nu)
    wrong <- NULL
    if(is_wrong) {
      l_units[[k]] <- wrong_entry_units(reports$time_reported, state$onset_times[k],
                                        reporting_parameters, revision_parameters)
      wrong <- list(
        lambda=revision_parameters$eta *
          (cases_true[k] - dated$out[k] + dated$into[k]) * l_units[[k]]$unit,
        w=l_units[[k]]$w, delta=state$delta, max_erroneous=max_erroneous)
    }
    f_forward <- function(deaths_true) {
      fatal_recorded <- deaths_true - dated$out[k] + dated$into[k]
      if(is_exact)
        return(joint_path_exact_forward(
          reports$cases_reported, reports$deaths_reported,
          deaths_true - dated$out[k], dated$into[k], cases_true[k] - deaths_true,
          survivals, wrong))
      joint_path_forward(
        reports$cases_reported, reports$deaths_reported, fatal_recorded,
        cases_true[k] - deaths_true,
        path_probabilities(survivals, ifelse(fatal_recorded > 0,
                                             dated$into[k] / fatal_recorded, 0)),
        wrong)
    }
    f_logp <- function(deaths_true, forward)
      forward$logp +
      stats::dbinom(deaths_true, cases_true[k], state$cfr[k], log=TRUE) +
      stats::dbinom(dated$out[k], deaths_true, dated$fraction[k], log=TRUE)
    deaths_true <- state$deaths_true[k]
    forward <- f_forward(deaths_true)
    if(!is.finite(forward$logp))
      stop("The reports of cases and deaths are impossible given the true cases and deaths.")
    for(proposal in seq_len(n_proposals)) {
      proposed <- deaths_true + sample(c(-1, 1), 1) * sample.int(step_max, 1)
      if(proposed < dated$out[k] | proposed > cases_true[k])
        next
      forward_proposed <- f_forward(proposed)
      log_ratio <- f_logp(proposed, forward_proposed) - f_logp(deaths_true, forward)
      if(is.finite(log_ratio) && log(stats::runif(1)) < log_ratio) {
        deaths_true <- proposed
        forward <- forward_proposed
      }
    }
    state$deaths_true[k] <- deaths_true
    paths[[k]] <- if(is_exact) joint_path_exact_backward(forward)
    else joint_path_backward(forward)
    l_survivals[[k]] <- survivals
    l_wrong[k] <- list(wrong)
  }
  state$paths <- paths
  state$survivals <- l_survivals
  state$wrong_units <- if(is_wrong) l_units
  # what the wrong records of each onset time were given, for the exact form
  state$wrong <- l_wrong
  state
}

#' Recalculates what is held about the delays of each onset time
#'
#' @inheritParams deaths_paths_step
#'
#' @return the state with elements 'survivals' and 'wrong_units' updated
deaths_refresh_survivals <- function(state, reporting_by_onset,
                                     revision_parameters=NULL) {
  delays <- deaths_delay_by_onset(state)
  is_wrong <- is_wrong_records(revision_parameters)
  for(k in seq_along(state$onset_times)) {
    days <- state$reports[[k]]$time_reported
    reporting_parameters <- delay_of_onset(reporting_by_onset, k)
    state$survivals[[k]] <- joint_survivals(
      days, state$onset_times[k], reporting_parameters, delay_of_onset(delays, k),
      revision_parameters, state$dated$reporting, state$recording_effects,
      state$nu)
    if(is_wrong)
      state$wrong_units[[k]] <- wrong_entry_units(
        days, state$onset_times[k], reporting_parameters, revision_parameters)
  }
  state
}

#' Creates the function giving the log-probability of the reports of an onset
#' time for possible true cases, with the true fatal cases summed out
#'
#' @inheritParams deaths_delay_by_onset
#' @param time_onset an onset time
#' @inheritParams deaths_paths_step
#'
#' @return a function of a vector of possible true cases. Its attribute
#' 'min_cases' is the smallest number of true cases which is possible.
deaths_cases_logp_function <- function(state, time_onset,
                                       revision_parameters=NULL) {
  k <- match(time_onset, state$onset_times)
  path <- state$paths[[k]]
  survivals <- state$survivals[[k]]
  cfr <- state$cfr[k]
  dated <- purrr::map(dated_counts(state), k)
  wrong_mean <- deaths_wrong_mean(state, k, revision_parameters)
  f <- function(cases_true) {
    logp <- joint_cases_deaths_logp(path, survivals, cfr, cases_true,
                                    dated$out, dated$into, dated$fraction,
                                    wrong_mean)
    matrixStats::rowLogSumExps(logp)
  }
  n_reports <- length(path$fatal_new)
  attr(f, "min_cases") <- dated$out +
    max(path$fatal_before[n_reports] + path$fatal_new[n_reports] -
          ifelse(isTRUE(path$is_exact), 0, dated$into), 0) +
    path$survivors_before[n_reports] + path$survivors_new[n_reports]
  f
}

#' Expected wrong records entered between reports for each genuine case
#' recorded for an onset time
#'
#' @inheritParams deaths_paths_step
#' @param k the onset time, counted from one
#'
#' @return NULL if wrong records are not allowed for; otherwise a vector with
#' an element for each report
deaths_wrong_mean <- function(state, k, revision_parameters=NULL) {
  if(!is_wrong_records(revision_parameters))
    return(NULL)
  revision_parameters$eta * state$wrong_units[[k]]$unit
}

#' Samples the true fatal cases of every onset time given the true cases
#'
#' @inheritParams deaths_paths_step
#'
#' @return the state with element 'deaths_true' updated
deaths_true_step <- function(state, cases_true) {
  dated <- dated_counts(state)
  for(k in seq_along(state$onset_times)) {
    # wrong records do not depend on how many of the true cases are fatal
    logp <- joint_cases_deaths_logp(state$paths[[k]], state$survivals[[k]],
                                    state$cfr[k], cases_true[k], dated$out[k],
                                    dated$into[k], dated$fraction[k])
    values <- attr(logp, "deaths_true")
    logp <- as.vector(logp)
    if(!any(is.finite(logp)))
      stop("The true cases of an onset time are fewer than the cases reported.")
    state$deaths_true[k] <- values[sample.int(length(values), 1,
                                              prob=exp(logp - max(logp)))]
  }
  state
}

#' Log-probability of what depends on the reporting delay in the reports of
#' one onset time, given how many of the cases in each report are fatal and
#' how many are wrong records
#'
#' This is the probability of the fatal and other cases newly reported in
#' each report, of the wrong records entered and removed and, where fatal
#' cases dated by death have their own delay, of the deaths newly shown
#' (since that delay and the reporting delay together decide how many of the
#' fatal cases newly reported are dated by death).
#'
#' @inheritParams deaths_paths_step
#' @param k the onset time, counted from one
#' @param reporting_parameters named list of the 'mean' and 'sd' of the
#' reporting delay
#' @param nu the log of the ratio of the reporting delay of fatal cases to
#' that of others
#' @param is_fatal_only if TRUE only the parts for fatal cases are included
#'
#' @return a log-probability
deaths_reporting_day_logp <- function(state, k, cases_true, reporting_parameters,
                                      revision_parameters=NULL, nu=state$nu,
                                      is_fatal_only=FALSE) {
  path <- state$paths[[k]]
  days <- state$reports[[k]]$time_reported
  dated <- purrr::map(dated_counts(state), k)
  f <- function(parameters)
    log_undetected_days(days, state$onset_times[k], parameters, revision_parameters)
  survivals <- state$survivals[[k]]
  survivals$log_reported <- f(reporting_parameters)
  survivals$log_reported_fatal <- survivals$log_reported
  if(nu != 0)
    survivals$log_reported_fatal <- f(list(mean=reporting_parameters$mean * exp(nu),
                                           sd=reporting_parameters$sd * exp(nu)))
  if(survivals$is_dated_same)
    survivals$log_dated <- survivals$log_reported_fatal
  else
    survivals$log_dated <- f(state$dated$reporting)
  # with a single delay, which fatal cases are shown as dead at once does
  # not depend on it
  logp <- fatal_reports_logp(
    path, survivals, dated$fatal_recorded, dated$into,
    is_deaths_included=dated$into > 0 & !survivals$is_dated_same &
      !isTRUE(path$is_exact))
  if(is_fatal_only)
    return(logp)
  logp <- logp + reports_stream_logp(
    path$survivors_new, path$survivors_before,
    cases_true[k] - state$deaths_true[k], survival_steps(survivals$log_reported))
  if(is_wrong_records(revision_parameters)) {
    units <- wrong_entry_units(days, state$onset_times[k], reporting_parameters,
                               revision_parameters)
    logp <- logp +
      sum(stats::dpois(path$wrong_entered_alive + path$wrong_entered_dead,
                       revision_parameters$eta *
                         (cases_true[k] - dated$out + dated$into) * units$unit,
                       log=TRUE)) +
      sum(dbinom_log_recycled(path$wrong_removed, path$wrong_before, units$w))
  }
  logp
}

#' As \code{deaths_reporting_day_logp}, summed over the onset times
#'
#' @inheritParams deaths_reporting_day_logp
#' @inheritParams deaths_paths_step
#'
#' @return a log-probability
deaths_reporting_logp <- function(state, cases_true, reporting_by_onset,
                                  revision_parameters=NULL, nu=state$nu,
                                  is_fatal_only=FALSE) {
  sum(purrr::map_dbl(seq_along(state$onset_times), function(k)
    deaths_reporting_day_logp(state, k, cases_true,
                              delay_of_onset(reporting_by_onset, k),
                              revision_parameters, nu, is_fatal_only)))
}

#' Creates a function giving the log-probability of the cases newly reported
#' for sets of onset times which share a reporting delay distribution
#'
#' This takes the place of \code{reporting_logp_function} in the joint model
#' of cases and deaths, where fatal and non-fatal cases are reported
#' separately.
#'
#' @inheritParams deaths_paths_step
#'
#' @return a function of 'indices' and 'reporting_parameters' (see
#' \code{reporting_logp_function})
deaths_reporting_logp_function <- function(state, cases_true,
                                           revision_parameters=NULL) {
  function(indices, reporting_parameters) {
    if(reporting_parameters$mean <= 0 | reporting_parameters$sd <= 0)
      return(-Inf)
    sum(purrr::map_dbl(indices, function(k)
      deaths_reporting_day_logp(state, k, cases_true, reporting_parameters,
                                revision_parameters)))
  }
}

#' Log-probability of the deaths newly shown in the reports of one onset
#' time, given how many of the cases in each report are fatal
#'
#' @inheritParams deaths_reporting_day_logp
#' @param death_parameters named list of the 'mean' and 'sd' of the delay
#' from onset to a death being recorded
#' @param recording_effects see \code{joint_survivals}
#'
#' @return a log-probability
deaths_recording_day_logp <- function(state, k, death_parameters,
                                      recording_effects=state$recording_effects) {
  path <- state$paths[[k]]
  dated <- purrr::map(dated_counts(state), k)
  survivals <- state$survivals[[k]]
  survivals$log_recorded <- log_undetected_days(
    state$reports[[k]]$time_reported, state$onset_times[k], death_parameters,
    recording_effects)
  probabilities <- path_probabilities(
    survivals, ifelse(dated$fatal_recorded > 0 & !isTRUE(path$is_exact),
                      dated$into / dated$fatal_recorded, 0))
  is_informative <- path$fatal_new > 0 | path$alive_before > 0
  sum(deaths_new_logp(path$deaths_new[is_informative],
                      path$fatal_new[is_informative],
                      path$alive_before[is_informative],
                      probabilities$g[is_informative],
                      probabilities$q[is_informative]))
}

#' Creates a function giving the log-probability of the deaths newly shown
#' for sets of onset times which share a delay from onset to a death being
#' recorded
#'
#' @inheritParams deaths_paths_step
#'
#' @return a function of 'indices' and 'reporting_parameters' (here the
#' parameters of the delay from onset to a death being recorded)
deaths_recording_logp_function <- function(state) {
  function(indices, reporting_parameters) {
    if(reporting_parameters$mean <= 0 | reporting_parameters$sd <= 0)
      return(-Inf)
    sum(purrr::map_dbl(indices, function(k)
      deaths_recording_day_logp(state, k, reporting_parameters)))
  }
}

#' Updates the delay from onset to a death being recorded and the case
#' fatality ratio, with the onset times at which each changes
#'
#' The effects of the day of the week on the delay to a death being recorded
#' are updated too, if there are any.
#'
#' @inheritParams deaths_paths_step
#' @inheritParams initial_deaths_state
#'
#' @return the updated state
deaths_parameters_step <- function(state, cases_true, settings) {
  n_days <- length(state$onset_times)
  state$recording <- sample_reporting_changepoints(
    state$recording,
    deaths_recording_logp_function(state),
    n_days, settings$delay_prior, settings$delay_changepoints)
  state <- deaths_recording_weekday_step(state, settings)

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
#' with the true fatal cases whose onset was on the day, the case fatality
#' ratio, the delay from onset to a death being recorded, the fatal cases
#' recorded for the day and how many of those are dated by death) and
#' 'deaths_parameters' (a single row, which includes 'fatal_delay_ratio', the
#' ratio of the reporting delay of fatal cases to that of others,
#' 'wrong_dead', the probability a wrong record is shown as dead, and the
#' effect of each day of the week on the delay to a death being recorded)
deaths_draws <- function(state, iteration) {
  delays <- deaths_delay_by_onset(state)
  dated <- dated_counts(state)
  parameters <- dplyr::tibble(
    iteration=iteration,
    n_cfr_changes=length(state$cfr_starts) - 1,
    p_change=state$p_change,
    n_delay_changes=length(state$recording$starts) - 1,
    delay_centre_mean=state$recording$centre$mean,
    delay_centre_sd=state$recording$centre$sd,
    delay_spread=state$recording$spread,
    fatal_delay_ratio=exp(state$nu),
    wrong_dead=state$delta)
  for(weekday in names(state$recording_effects$weekday_log_effects))
    parameters[[paste0("recording_weekday_", weekday)]] <-
      exp(state$recording_effects$weekday_log_effects[[weekday]])
  list(
    deaths=dplyr::tibble(
      iteration=iteration,
      time_onset=state$onset_times,
      deaths_true=state$deaths_true,
      cfr=state$cfr,
      delay_mean=delays$mean,
      delay_sd=delays$sd,
      deaths_recorded=dated$fatal_recorded,
      dated_by_death=dated$into),
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
#' Optionally a fraction of fatal cases are recorded for their day of death,
#' not their day of onset. Such a case is counted on the day it died and is
#' shown as dead as soon as it is reported. It is reported after a delay
#' counted from that day: one with its own supplied mean and sd or, if none
#' is supplied, the reporting delay of other fatal cases.
#'
#' @param cases_true a vector of true cases on onset times 1, 2, ...
#' @param report_days the times at which reports are made
#' @param reporting_parameters the delay from onset to a case being reported:
#' a named list of 'mean' and 'sd' or a tibble with a row for each onset time
#' @param cfr the case fatality ratio: a single value or one for each onset time
#' @param death_parameters the delay from onset to a death being recorded, in
#' the same form as reporting_parameters
#' @param fatal_reporting_parameters the delay from onset to a fatal case
#' being reported, if it differs from that of other cases (optional), in the
#' same form as reporting_parameters
#' @param death_from "onset" (the default) if the delay to a death being
#' recorded runs from onset independently of when the case is reported, or
#' "report" if it runs from when the case is reported
#' @param dated_by_death an optional named list of 'fraction' (the fraction of
#' fatal cases which are recorded for their day of death), 'mean' and 'sd' (of
#' the gamma-distributed time from onset to death of those cases, which is
#' rounded up to a whole number of days) and, optionally, 'reporting' (a named
#' list of the 'mean' and 'sd' of the delay from the day of death to those
#' cases being reported)
#' @param wrong an optional named list of 'eta', 'rho' and 'delta', which adds
#' wrong records to the reports. On each day after the day for which cases
#' are recorded, a Poisson number of wrong records is entered, with mean eta
#' times the cases recorded for the day times the probability a case is first
#' reported on that day. Each is shown as dead with probability delta and is
#' removed with probability rho on each later day.
#'
#' @return a named list of 'snapshots' (a tibble with columns time_onset,
#' time_reported, cases_reported and deaths_reported, where time_onset is the
#' day for which cases are recorded) and 'truth' (a tibble with a row for
#' each such day and columns: time_onset, cases_true and deaths_true (the
#' cases and fatal cases recorded for the day), cases_onset and deaths_onset
#' (the cases and fatal cases whose onset was on the day), cfr,
#' n_dated_by_death (the fatal cases recorded for the day which are dated by
#' death), and the mean and sd of the delay from the day to being shown as
#' dead)
#' @export
#' @importFrom rlang .data
simulate_cases_deaths <- function(cases_true, report_days, reporting_parameters,
                                  cfr, death_parameters,
                                  fatal_reporting_parameters=NULL,
                                  death_from=c("onset", "report"),
                                  dated_by_death=NULL, wrong=NULL) {
  death_from <- match.arg(death_from)
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
  fatal_reporting <- reporting
  if(!is.null(fatal_reporting_parameters))
    fatal_reporting <- f_expand(fatal_reporting_parameters)
  if(length(cfr) == 1)
    cfr <- rep(cfr, n_days)
  if(length(cfr) != n_days)
    stop("There must be a case fatality ratio for each onset time.")
  f_delay <- function(mean, sd)
    stats::rgamma(length(mean), shape=mean^2 / sd^2, rate=mean / sd^2)

  # one row for each case
  cases_by_onset <- cases_true
  onset <- rep(seq_len(n_days), cases_true)
  n <- length(onset)
  is_fatal <- stats::runif(n) < cfr[onset]
  is_dated_by_death <- rep(FALSE, n)
  day <- onset
  if(!is.null(dated_by_death)) {
    is_dated_by_death <- is_fatal & stats::runif(n) < dated_by_death$fraction
    time_to_death <- f_delay(rep(dated_by_death$mean, n), rep(dated_by_death$sd, n))
    # rounded up, so that a case dated by its death is at least a day late
    day[is_dated_by_death] <- onset[is_dated_by_death] +
      ceiling(time_to_death[is_dated_by_death])
  }
  # delays are counted from the day for which the case is recorded
  delay_reported <- ifelse(
    is_fatal,
    f_delay(fatal_reporting$mean[onset], fatal_reporting$sd[onset]),
    f_delay(reporting$mean[onset], reporting$sd[onset]))
  if(!is.null(dated_by_death$reporting)) {
    n_dated <- sum(is_dated_by_death)
    delay_reported[is_dated_by_death] <- f_delay(
      rep(dated_by_death$reporting$mean, n_dated),
      rep(dated_by_death$reporting$sd, n_dated))
  }
  delay_recorded <- f_delay(death$mean[onset], death$sd[onset])
  if(death_from == "report")
    delay_recorded <- delay_reported + delay_recorded
  delay_dead <- ifelse(is_dated_by_death, delay_reported,
                       pmax(delay_reported, delay_recorded))

  results <- purrr::map(seq_len(n_days), function(t) {
    is_day <- day == t
    days <- report_days[report_days >= t]
    elapsed <- days - t + 0.5
    wrong_alive <- rep(0, length(days))
    wrong_dead <- rep(0, length(days))
    if(!is.null(wrong) && max(days) > t) {
      # the day each wrong record is entered and the day it is removed
      days_after <- seq(t + 1, max(days), 1)
      shape <- reporting$mean[t]^2 / reporting$sd[t]^2
      rate <- reporting$mean[t] / reporting$sd[t]^2
      p_first_reported <- diff(stats::pgamma(c(t, days_after) + 0.5 - t,
                                             shape=shape, rate=rate))
      entered <- rep(days_after, stats::rpois(
        length(days_after), wrong$eta * sum(is_day) * p_first_reported))
      removed <- entered + stats::rgeom(length(entered), wrong$rho) + 1
      is_dead <- stats::runif(length(entered)) < wrong$delta
      f_present <- function(is_kind)
        purrr::map_dbl(days, ~sum(entered[is_kind] <= . & removed[is_kind] > .))
      wrong_alive <- f_present(!is_dead)
      wrong_dead <- f_present(is_dead)
    }
    list(
      snapshots=dplyr::tibble(
        time_onset=t,
        time_reported=days,
        cases_reported=purrr::map_dbl(elapsed, ~sum(delay_reported[is_day] <= .)) +
          wrong_alive + wrong_dead,
        deaths_reported=purrr::map_dbl(
          elapsed, ~sum(delay_dead[is_day & is_fatal] <= .)) + wrong_dead),
      truth=dplyr::tibble(
        time_onset=t, cases_true=sum(is_day), deaths_true=sum(is_day & is_fatal),
        cases_onset=cases_by_onset[t], deaths_onset=sum(onset == t & is_fatal),
        cfr=cfr[t],
        n_dated_by_death=sum(is_day & is_dated_by_death),
        delay_dead_mean=mean(delay_dead[is_day & is_fatal]),
        delay_dead_sd=stats::sd(delay_dead[is_day & is_fatal])))
  })
  list(snapshots=purrr::map_dfr(results, "snapshots"),
       truth=purrr::map_dfr(results, "truth"))
}
