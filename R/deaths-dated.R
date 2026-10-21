# Fatal cases which are recorded for their day of death.
#
# For some of those who die, the day the reports give for the case is the day
# of death and not the day of onset. The reports then count cases by the day
# they are recorded for, while transmission and the case fatality ratio are
# by true day of onset. The model holds, for each day of onset, how many of
# its fatal cases are dated by death and when they died, and moves these
# cases between the two sets of counts.
#
# The distribution of the time from onset to death of these cases is supplied.
# The fraction of fatal cases which are dated by death is supplied or, with a
# prior, estimated. Their delay from death to being reported is supplied, or
# is that of other fatal cases counted from the day of death, or is the delay
# from a death to its being recorded.

#' Settings for fatal cases which are recorded for their day of death
#'
#' @param n_days number of onset times
#' @param dated NULL, or a named list of: 'fraction' (the fraction of the
#' fatal cases with onset on a day which are recorded for their day of death:
#' a single value or one for each onset time if it is known; or, if it is to
#' be estimated, a named list of 'alpha' and 'beta', the parameters of its
#' beta prior, and optionally 'initial', the value to start from. An
#' estimated fraction is the same for every onset time); 'mean' and 'sd' (of
#' the time from onset to death of those cases, which has a gamma distribution
#' rounded up to a whole number of days, so is at least one day; by default
#' those of onset_to_death); and optionally
#' 'max_lag' (the longest such time allowed for, by default the 99.9th
#' percentile), 'reporting' (a named list of the 'mean' and 'sd' of the delay
#' from the day of death to those cases being reported, or "recording" if
#' they are reported after the delay from a death to its being recorded,
#' which needs onset_to_death and is the default with it; if absent
#' otherwise this is the
#' reporting delay of other cases, counted from the day of death) and
#' 'n_moves' (the number of changes proposed in each iteration, by default
#' four for each onset time, or one in the exact form) and 'exact' (if TRUE,
#' the fatal cases dated by onset and by death among those reported are
#' followed separately, which is exact and slower; by default each fatal
#' case recorded for a day is taken to be dated by death with a probability)
#' and 'collapse_fraction' (if TRUE, and the fraction is estimated, the
#' changes to which fatal cases are dated by death are made with the fraction
#' integrated out, and the fraction is then drawn exactly given them: see
#' \code{dated_fraction_marginal}. By default FALSE, the changes are made at
#' the current fraction)
#' @param onset_to_death NULL, or the element of that name of the result of
#' \code{deaths_settings}: the time from onset to death of every fatal case
#'
#' @return NULL if no fatal cases are dated by death; otherwise the list with
#' 'fraction' given for each onset time, defaults filled in,
#' 'lag_probabilities' added (the probability of each time from onset to
#' death, from one day to max_lag) and 'fraction_prior' (NULL if the fraction
#' is known; otherwise a named list of 'alpha' and 'beta')
dated_by_death_settings <- function(n_days, dated, onset_to_death=NULL) {
  if(is.null(dated))
    return(NULL)
  if(!is.null(onset_to_death)) {
    for(name in c("mean", "sd", "max_lag"))
      if(is.null(dated[[name]]))
        dated[[name]] <- onset_to_death[[name]]
    if(is.null(dated$reporting))
      dated$reporting <- "recording"
  }
  fraction_prior <- NULL
  if(is.list(dated$fraction)) {
    fraction_prior <- list(alpha=1, beta=1, initial=0.2)
    unknown <- setdiff(names(dated$fraction), names(fraction_prior))
    if(length(unknown) > 0)
      stop(paste0("Unknown settings for the fraction dated by death: ",
                  paste(unknown, collapse=", "), "."))
    fraction_prior[names(dated$fraction)] <- dated$fraction
    if(fraction_prior$alpha <= 0 | fraction_prior$beta <= 0)
      stop("The parameters of the prior for the fraction dated by death must be positive.")
    if(fraction_prior$initial <= 0 | fraction_prior$initial >= 1)
      stop("The initial fraction dated by death must lie strictly between 0 and 1.")
    dated$fraction <- fraction_prior$initial
    fraction_prior$initial <- NULL
  }
  unknown <- setdiff(names(dated),
                     c("fraction", "mean", "sd", "max_lag", "reporting", "n_moves",
                       "exact", "collapse_fraction"))
  if(length(unknown) > 0)
    stop(paste0("Unknown settings for dated_by_death: ",
                paste(unknown, collapse=", "), "."))
  fraction <- dated$fraction
  if(is.null(fraction) || !length(fraction) %in% c(1, n_days))
    stop("The fraction dated by death must be a single value or one for each onset time.")
  if(any(is.na(fraction)) || any(fraction < 0 | fraction > 1))
    stop("The fraction dated by death must lie between 0 and 1.")
  if(all(fraction == 0) & is.null(fraction_prior))
    return(NULL)
  is_collapsed <- isTRUE(dated$collapse_fraction)
  if(is_collapsed & is.null(fraction_prior))
    stop("The fraction dated by death can be integrated out only if it is estimated.")
  if(is.null(dated$mean) || is.null(dated$sd) || dated$mean <= 0 || dated$sd <= 0)
    stop("A positive mean and sd of the time from onset to death must be given.")
  shape <- dated$mean^2 / dated$sd^2
  rate <- dated$mean / dated$sd^2
  max_lag <- dated$max_lag
  if(is.null(max_lag))
    max_lag <- max(1, ceiling(stats::qgamma(0.999, shape=shape, rate=rate)))
  if(max_lag < 1)
    stop("max_lag must be at least one day.")
  reporting <- dated$reporting
  if(identical(reporting, "recording")) {
    if(is.null(onset_to_death))
      stop("Fatal cases dated by death can be reported after the delay from a death to its being recorded only if onset_to_death is supplied.")
  } else if(!is.null(reporting) &&
     (!is.list(reporting) || is.null(reporting$mean) || is.null(reporting$sd) ||
      reporting$mean <= 0 || reporting$sd <= 0))
    stop("The delay from death to being reported must have a positive mean and sd.")
  lag_probabilities <- diff(stats::pgamma(0:max_lag, shape=shape, rate=rate))
  is_exact <- isTRUE(dated$exact)
  n_moves <- dated$n_moves
  if(is.null(n_moves))
    n_moves <- ifelse(is_exact, 1, 4) * n_days
  list(fraction=rep_len(fraction, n_days),
       fraction_prior=fraction_prior,
       mean=dated$mean, sd=dated$sd,
       max_lag=max_lag,
       lag_probabilities=lag_probabilities / sum(lag_probabilities),
       reporting=reporting,
       n_moves=n_moves,
       exact=is_exact,
       collapse_fraction=is_collapsed)
}

#' Adds fatal cases dated by death to the starting state of the model of
#' deaths
#'
#' The state as built holds the cases and fatal cases recorded for each day.
#' A share of the fatal cases recorded for each day, equal to the fraction
#' supplied, is taken to be dated by death, with onsets on earlier days drawn
#' from the distribution of the time from onset to death. The true cases and
#' fatal cases by day of onset follow.
#'
#' @param state the state as built by \code{initial_deaths_state}
#' @param dated as returned by \code{dated_by_death_settings}
#' @inheritParams sample_cases_history
#'
#' @return the state with elements 'deaths_true' and 'cases_onset' now by day
#' of onset, and an element 'dated': a named list holding what is supplied
#' and
#' \itemize{
#' \item 'lags': a matrix whose element (k, z) is the number of fatal cases
#' with onset on day k which are dated by death and died z days later, on a
#' day for which there are reports;
#' \item 'unseen': for each day, the number of its fatal cases which are dated
#' by death and die after the last day, so that they are in no report;
#' \item 'inflow': for each day, the number of fatal cases recorded for it
#' whose onset was before the first day;
#' \item 'inflow_weight': for each day, the number of such cases expected for
#' each unit of the product of the case fatality ratio and the fraction dated
#' by death (zero if no earlier cases are known);
#' \item 'tail': for each day, the probability that a fatal case with onset
#' on it which is dated by death dies after the last day;
#' \item 'out' and 'into': see \code{dated_counts}.
#' }
initial_dated_state <- function(state, dated, cases_before=NULL) {
  n_days <- length(state$onset_times)
  max_lag <- dated$max_lag
  h <- dated$lag_probabilities
  fatal_recorded <- state$deaths_true
  survivors <- state$cases_onset - fatal_recorded

  # deaths after the last day
  days_left <- n_days - seq_len(n_days)
  tail <- purrr::map_dbl(days_left, ~sum(h[seq_len(max_lag) > .]))
  # onsets before the first day: cases_before ends on the day before it
  n_before <- length(cases_before)
  inflow_weight <- purrr::map_dbl(seq_len(n_days), function(t) {
    lags <- seq_len(max_lag)
    index <- n_before + t - lags
    is_earlier <- lags >= t & index >= 1
    sum(cases_before[index[is_earlier]] * h[is_earlier])
  })

  lags <- matrix(0, n_days, max_lag)
  inflow <- rep(0, n_days)
  for(t in seq_len(n_days)) {
    # a case dated by death is shown as dead once reported, so to start with
    # there are no more of them than the deaths in the latest report
    n_dated <- min(round(dated$fraction[t] * fatal_recorded[t]),
                   utils::tail(state$reports[[t]]$deaths_reported, 1))
    reports <- state$reports[[t]]
    # in the exact form each enters as a case and a death in the same report,
    # so there are no more than the reports' increases in both allow
    if(isTRUE(dated$exact))
      n_dated <- min(n_dated, sum(pmin(
        pmax(diff(c(0, reports$cases_reported)), 0),
        pmax(diff(c(0, reports$deaths_reported)), 0))))
    # where a death is recorded only after the time from onset to death, a
    # death shown on the day itself must be that of a case dated by death (in
    # the exact form, every such death)
    shown_same_day <- max(c(0, reports$deaths_reported[
      reports$time_reported == state$onset_times[t]]))
    if(identical(dated$reporting, "recording") && shown_same_day > 0)
      n_dated <- max(n_dated, ifelse(isTRUE(dated$exact), shown_same_day, 1))
    weights <- h * (seq_len(max_lag) < t)
    weight_inflow <- (inflow_weight[t] > 0) * sum(h[seq_len(max_lag) >= t])
    if(n_dated == 0 | sum(weights) + weight_inflow == 0)
      next
    counts <- as.vector(stats::rmultinom(1, n_dated, c(weights, weight_inflow)))
    inflow[t] <- counts[max_lag + 1]
    for(z in which(counts[seq_len(max_lag)] > 0))
      lags[t - z, z] <- lags[t - z, z] + counts[z]
  }
  into <- purrr::map_dbl(seq_len(n_days), function(t) {
    z <- seq_len(min(max_lag, t - 1))
    sum(lags[cbind(t - z, z)])
  }) + inflow
  by_onset <- fatal_recorded - into
  seen <- by_onset + rowSums(lags)
  p_unseen <- dated$fraction * tail
  unseen <- ifelse(p_unseen < 1, round(seen * p_unseen / (1 - p_unseen)), 0)
  out <- rowSums(lags) + unseen

  state$deaths_true <- by_onset + out
  state$cases_onset <- survivors + state$deaths_true
  state$dated <- c(
    dated,
    list(lags=lags, unseen=unseen, inflow=inflow, inflow_weight=inflow_weight,
         tail=tail, out=out, into=into))
  state
}

#' What the renewal model needs to give the log-probability of a series of
#' true cases
#'
#' @param Rt a vector of the reproduction number on each onset time
#' @param kappa the overdispersion parameter of a negative binomial renewal
#' model
#' @param is_negative_binomial whether the renewal model is negative binomial
#' (if FALSE it is Poisson)
#' @inheritParams renewal_statistics
#'
#' @return a named list of Rt, kappa, is_negative_binomial, 'cases_earlier'
#' (the serial_max days of cases before the first, of which any not known are
#' zero) and 'weights' (a matrix which, multiplied by the earlier cases
#' followed by the true cases, gives the infectiousness on each onset time)
dated_renewal <- function(Rt, kappa, is_negative_binomial, serial_parameters,
                          serial_max=40, cases_before=NULL) {
  n_days <- length(Rt)
  w <- weights_series(serial_max, serial_parameters)
  weights <- matrix(0, n_days, serial_max + n_days)
  for(t in seq_len(n_days))
    weights[t, t:(t + serial_max - 1)] <- rev(w)
  list(Rt=Rt, kappa=kappa, is_negative_binomial=is_negative_binomial,
       cases_earlier=rev(recent_history(rev(cases_before), serial_max)),
       weights=weights)
}

#' Log-probability of a series of true cases under the renewal model
#'
#' Days with no infectiousness before them are left out, as in
#' \code{renewal_statistics}.
#'
#' @param cases_true a vector of the true cases of each onset time
#' @param renewal as returned by \code{dated_renewal}
#'
#' @return a log-probability
dated_renewal_logp <- function(cases_true, renewal) {
  infectiousness <- as.vector(
    renewal$weights %*% c(renewal$cases_earlier, cases_true))
  is_informative <- infectiousness > 0
  mean_cases <- (renewal$Rt * infectiousness)[is_informative]
  if(renewal$is_negative_binomial)
    sum(stats::dnbinom(cases_true[is_informative], mu=mean_cases,
                       size=renewal$kappa, log=TRUE))
  else
    sum(stats::dpois(cases_true[is_informative], lambda=mean_cases, log=TRUE))
}

#' x log(y), which is zero where x is zero
#'
#' @param x a vector
#' @param y a vector
#'
#' @return a vector
x_log_y <- function(x, y) {
  ifelse(x == 0, 0, x * log(y))
}

#' Log-probability of the parts of the model which a change to the fatal
#' cases dated by death can alter, other than the renewal model
#'
#' For each day of onset given, these are the probability of its true fatal
#' cases given its true cases, and the probability of how its fatal cases
#' divide into those dated by onset, those dated by death which died after
#' each number of days, and those dated by death which die after the last
#' day. For each recorded day given, they are the probability of the fatal
#' cases recorded for it whose onset was before the first day, and the
#' probability of its reports given how many of the cases in each are fatal
#' (the cases which are not fatal are left out, since no such change alters
#' them).
#'
#' @inheritParams deaths_delay_by_onset
#' @param cases_true a vector of the true cases of each onset time
#' @param cohorts the days of onset to include, counted from one
#' @param days the recorded days to include, counted from one
#' @param collapsed if TRUE, the factors in which the fraction dated by death
#'   appears are left out (they are those of \code{dated_fraction_terms},
#'   which \code{dated_fraction_marginal} integrates over the fraction)
#'
#' @return a log-probability (-Inf if the state is impossible)
dated_logp <- function(state, cases_true, cohorts=seq_along(state$onset_times),
                       days=seq_along(state$onset_times), collapsed=FALSE) {
  dated <- state$dated
  if(any(cases_true < 0))
    return(-Inf)
  logp <- 0
  for(k in cohorts) {
    fatal <- state$deaths_true[k]
    lags <- dated$lags[k, ]
    unseen <- dated$unseen[k]
    by_onset <- fatal - sum(lags) - unseen
    if(by_onset < 0 | fatal > cases_true[k] | any(lags < 0) | unseen < 0)
      return(-Inf)
    # with the fraction left out, the factors f^(dated) (1 - f)^(by onset)
    fraction <- if(collapsed) 1 else dated$fraction[k]
    logp <- logp +
      lgamma(fatal + 1) - lgamma(by_onset + 1) - sum(lgamma(lags + 1)) -
      lgamma(unseen + 1) +
      (if(collapsed) 0 else x_log_y(by_onset, 1 - fraction)) +
      sum(x_log_y(lags, fraction * dated$lag_probabilities)) +
      x_log_y(unseen, fraction * dated$tail[k]) +
      stats::dbinom(fatal, cases_true[k], state$cfr[k], log=TRUE)
  }
  for(t in days) {
    if(dated$inflow[t] < 0)
      return(-Inf)
    if(dated$inflow[t] > 0 | dated$inflow_weight[t] > 0) {
      mean_inflow <- dated$inflow_weight[t] * state$cfr[1]
      # with the fraction left out, the factors f^(inflow) exp(-mean f)
      logp <- logp + if(collapsed)
        x_log_y(dated$inflow[t], mean_inflow) - lgamma(dated$inflow[t] + 1)
      else
        stats::dpois(dated$inflow[t], mean_inflow * dated$fraction[1], log=TRUE)
    }
    # in the exact form a change of kind is impossible given how many of each
    # kind are in each report, so that is summed out
    if(dated$exact)
      logp <- logp + dated_exact_forward(state, cases_true, t)$logp
    else
      logp <- logp + fatal_reports_logp(
        state$paths[[t]], state$survivals[[t]],
        state$deaths_true[t] - dated$out[t] + dated$into[t], dated$into[t])
  }
  logp
}

#' Filters forwards over the reports of one onset time in the exact form, at
#' the current numbers of fatal cases of each kind
#'
#' @inheritParams dated_logp
#' @param t the onset time, counted from one
#'
#' @return as \code{joint_path_exact_forward}
dated_exact_forward <- function(state, cases_true, t) {
  reports <- state$reports[[t]]
  wrong <- state$wrong[[t]]
  # the wrong records are held at those of the current path (see
  # deaths_paths_step)
  fixed <- if(!is.null(wrong) && wrong$max_erroneous > 0)
    list(wrong=exact_fixed_wrong(state$paths[[t]]))
  joint_path_exact_forward(
    reports$cases_reported, reports$deaths_reported,
    state$deaths_true[t] - state$dated$out[t], state$dated$into[t],
    cases_true[t] - state$deaths_true[t], state$survivals[[t]], wrong, fixed)
}

#' Applies a change to the fatal cases dated by death
#'
#' There are four kinds of change, none of which alters the cases or the
#' fatal cases recorded for any day:
#' \itemize{
#' \item 'relabel': a fatal case recorded for a day changes from one with
#' onset on that day, dated by onset, to one with onset 'lag' days earlier
#' which died on that day, dated by death (sign 1), or the reverse (sign -1).
#' This moves a true case between the two days of onset.
#' \item 'redate': a fatal case dated by death and recorded for a day changes
#' from one which died 'lag' days after onset to one which died 'lag_to' days
#' after onset, which moves a true case between two days of onset.
#' \item 'unseen': a fatal case with onset on a day which is dated by death
#' and dies after the last day is added (sign 1) or removed (sign -1).
#' \item 'inflow': a fatal case recorded for a day changes from one with
#' onset on that day to one with onset before the first day (sign 1), or the
#' reverse (sign -1).
#' }
#'
#' @inheritParams dated_logp
#' @param move a named list of 'type', 'day', 'lag', 'lag_to' and 'sign'
#'
#' @return NULL if the change is not possible; otherwise a named list of the
#' changed 'state' and 'cases_true', and the 'cohorts' and 'days' which
#' \code{dated_logp} needs to give the change in log-probability
dated_apply_move <- function(state, cases_true, move) {
  n_days <- length(state$onset_times)
  day <- move$day
  sign <- move$sign
  f_onset <- function(k, change) {
    state$deaths_true[k] <<- state$deaths_true[k] + change
    cases_true[k] <<- cases_true[k] + change
  }
  f_lag <- function(k, lag, change) {
    state$dated$lags[k, lag] <<- state$dated$lags[k, lag] + change
    state$dated$out[k] <<- state$dated$out[k] + change
    f_onset(k, change)
  }
  if(move$type == "relabel") {
    k <- day - move$lag
    if(k < 1)
      return(NULL)
    f_lag(k, move$lag, sign)
    f_onset(day, -sign)
    state$dated$into[day] <- state$dated$into[day] + sign
    cohorts <- c(k, day)
    days <- day
  } else if(move$type == "redate") {
    k <- day - move$lag
    k_to <- day - move$lag_to
    if(k < 1 | k_to < 1 | k == k_to)
      return(NULL)
    f_lag(k, move$lag, -1)
    f_lag(k_to, move$lag_to, 1)
    cohorts <- c(k, k_to)
    days <- integer(0)
  } else if(move$type == "unseen") {
    if(state$dated$tail[day] <= 0)
      return(NULL)
    state$dated$unseen[day] <- state$dated$unseen[day] + sign
    state$dated$out[day] <- state$dated$out[day] + sign
    f_onset(day, sign)
    cohorts <- day
    days <- integer(0)
  } else if(move$type == "inflow") {
    if(state$dated$inflow_weight[day] <= 0)
      return(NULL)
    state$dated$inflow[day] <- state$dated$inflow[day] + sign
    state$dated$into[day] <- state$dated$into[day] + sign
    f_onset(day, -sign)
    cohorts <- day
    days <- day
  } else {
    stop("Unknown kind of change to the fatal cases dated by death.")
  }
  list(state=state, cases_true=cases_true, cohorts=cohorts, days=days)
}

#' Updates which fatal cases are dated by death
#'
#' Changes of the kinds described in \code{dated_apply_move} are proposed and
#' accepted with the Metropolis probability, which involves the renewal model
#' and \code{dated_logp}. Since the changes move true cases between days of
#' onset, the true cases are updated too.
#'
#' @inheritParams dated_logp
#' @inheritParams dated_renewal_logp
#'
#' @return a named list of the updated 'state' and 'cases_true'
deaths_dated_step <- function(state, cases_true, renewal) {
  dated <- state$dated
  if(is.null(dated))
    return(list(state=state, cases_true=cases_true))
  n_days <- length(state$onset_times)
  types <- c("relabel", "redate", "unseen", "inflow")
  days_unseen <- which(dated$tail > 0)
  days_inflow <- which(dated$inflow_weight > 0)
  type_weights <- c(0.5, 0.2, 0.2 * (length(days_unseen) > 0),
                    0.1 * (length(days_inflow) > 0))
  f_one <- function(values) values[sample.int(length(values), 1)]
  f_lag <- function() sample.int(dated$max_lag, 1, prob=dated$lag_probabilities)
  logp_renewal <- dated_renewal_logp(cases_true, renewal)
  # with the fraction integrated out, each change alters the total numbers of
  # each kind, and so the marginal of the fraction's factors
  collapsed <- isTRUE(dated$collapse_fraction) & !is.null(dated$fraction_prior)
  logp_marginal <- if(collapsed) dated_fraction_marginal(state)$log_marginal else 0
  for(i in seq_len(dated$n_moves)) {
    type <- types[sample.int(4, 1, prob=type_weights)]
    move <- list(
      type=type,
      day=switch(type, unseen=f_one(days_unseen), inflow=f_one(days_inflow),
                 sample.int(n_days, 1)),
      lag=f_lag(), lag_to=f_lag(), sign=f_one(c(-1, 1)))
    proposed <- dated_apply_move(state, cases_true, move)
    if(is.null(proposed))
      next
    logp_proposed <- dated_logp(proposed$state, proposed$cases_true,
                                proposed$cohorts, proposed$days, collapsed)
    if(!is.finite(logp_proposed))
      next
    logp_renewal_proposed <- dated_renewal_logp(proposed$cases_true, renewal)
    logp_marginal_proposed <- if(collapsed)
      dated_fraction_marginal(proposed$state)$log_marginal else 0
    log_ratio <- logp_proposed + logp_renewal_proposed + logp_marginal_proposed -
      dated_logp(state, cases_true, proposed$cohorts, proposed$days, collapsed) -
      logp_renewal - logp_marginal
    if(is.finite(log_ratio) && log(stats::runif(1)) < log_ratio) {
      state <- proposed$state
      cases_true <- proposed$cases_true
      logp_renewal <- logp_renewal_proposed
      logp_marginal <- logp_marginal_proposed
      # in the exact form, how many of each kind are in each report is drawn
      # again for the onset time whose numbers changed
      if(dated$exact)
        for(t in proposed$days)
          state$paths[[t]] <- joint_path_exact_redraw(
            dated_exact_forward(state, cases_true, t))
    }
  }
  list(state=state, cases_true=cases_true)
}

#' The factors of the model in which the fraction dated by death appears
#'
#' Given how many fatal cases are dated by death and by onset, and how many
#' recorded for the first days had their onset before the first day, these
#' factors, with the fraction's beta prior, are
#'   f^(A - 1) (1 - f)^(B - 1) exp(-r f),
#' where A is the prior's alpha plus the fatal cases dated by death and those
#' from before the first day, B its beta plus the fatal cases dated by onset,
#' and r the expected number from before the first day for each unit of the
#' fraction.
#'
#' @inheritParams deaths_delay_by_onset
#'
#' @return a named list of 'A', 'B' and 'r'
dated_fraction_terms <- function(state) {
  dated <- state$dated
  prior <- dated$fraction_prior
  n_dated <- sum(dated$out)
  list(A=prior$alpha + n_dated + sum(dated$inflow),
       B=prior$beta + sum(state$deaths_true) - n_dated,
       r=sum(dated$inflow_weight) * state$cfr[1])
}

#' The fraction dated by death integrated out
#'
#' The integral over the fraction of the factors of
#' \code{dated_fraction_terms} (with the beta prior's normalising constant)
#' is exact: since exp(-r f) = exp(-r) sum_j r^j (1 - f)^j / j!,
#'   integral = exp(-r) sum_j r^j / j! B(A, B + j) / B(alpha, beta),
#' a beta function times Kummer's function 1F1(B; A + B; r) exp(-r), whose
#' terms are all positive and summed until the rest is negligible. The terms
#' also give the distribution of the fraction as a mixture of beta
#' distributions, Beta(A, B + j) with weights proportional to them.
#'
#' @inheritParams deaths_delay_by_onset
#'
#' @return a named list of 'log_marginal' (the log of the integral), 'A', 'B'
#' and 'log_weights' (the log of the unnormalised weight of each j from 0)
dated_fraction_marginal <- function(state) {
  terms <- dated_fraction_terms(state)
  A <- terms$A
  B <- terms$B
  r <- terms$r
  prior <- state$dated$fraction_prior
  log_weights <- 0
  if(r > 0) {
    # the terms rise while (B + j) r / ((A + B + j) (j + 1)) > 1, so at most
    # until j = r, then fall faster than r^j / j!
    n <- ceiling(r + 12 * sqrt(r) + 40)
    repeat {
      j <- 0:n
      log_weights <- lgamma(B + j) - lgamma(B) - lgamma(A + B + j) +
        lgamma(A + B) + j * log(r) - lgamma(j + 1)
      if(log_weights[n + 1] < max(log_weights) - 40)
        break
      n <- 2 * n
    }
  }
  log_sum <- max(log_weights) + log(sum(exp(log_weights - max(log_weights))))
  list(log_marginal=lbeta(A, B) - lbeta(prior$alpha, prior$beta) - r + log_sum,
       A=A, B=B, log_weights=log_weights)
}

#' Draws the fraction dated by death exactly given which fatal cases are
#'
#' The distribution is the mixture of beta distributions of
#' \code{dated_fraction_marginal}: j is drawn with probability proportional
#' to its weight, then the fraction from Beta(A, B + j).
#'
#' @inheritParams deaths_delay_by_onset
#'
#' @return a value of the fraction
dated_fraction_draw <- function(state) {
  marginal <- dated_fraction_marginal(state)
  weights <- exp(marginal$log_weights - max(marginal$log_weights))
  j <- sample.int(length(weights), 1, prob=weights) - 1
  value <- stats::rbeta(1, marginal$A, marginal$B + j)
  # a draw at either end (with a very large count) is kept inside
  min(max(value, .Machine$double.eps), 1 - .Machine$double.eps)
}

#' Updates the fraction of fatal cases which are dated by death
#'
#' Given how many of the fatal cases of each day of onset are dated by death,
#' the fraction (one value for all onset times, with a beta prior) has a beta
#' distribution but for the fatal cases recorded for the first days whose
#' onset was before the first day, whose expected number is proportional to
#' it. A value is drawn from the beta distribution with the same mode and
#' curvature there as the distribution of the fraction, and is accepted with
#' the Metropolis-Hastings probability.
#'
#' @inheritParams deaths_delay_by_onset
#' @param n_steps the number of values proposed
#'
#' @return the state with element 'fraction' of its element 'dated' updated
deaths_dated_fraction_step <- function(state, n_steps=2) {
  dated <- state$dated
  prior <- dated$fraction_prior
  if(is.null(dated) || is.null(prior))
    return(state)
  # after changes made with the fraction integrated out, the fraction must be
  # drawn exactly from its distribution given them
  if(isTRUE(dated$collapse_fraction)) {
    state$dated$fraction <- rep(dated_fraction_draw(state), length(dated$fraction))
    return(state)
  }
  n_dated <- sum(dated$out)
  # log density: a log(f) + b log(1 - f) - rate f
  a <- prior$alpha - 1 + n_dated + sum(dated$inflow)
  b <- prior$beta - 1 + sum(state$deaths_true) - n_dated
  rate <- sum(dated$inflow_weight) * state$cfr[1]
  f_logp <- function(x) x_log_y(a, x) + x_log_y(b, 1 - x) - rate * x
  shape_1 <- a + 1
  shape_2 <- b + 1
  if(rate > 0 & a > 0 & b > 0) {
    mode <- (a + b + rate - sqrt((a + b + rate)^2 - 4 * rate * a)) / (2 * rate)
    if(is.finite(mode) && mode > 0 && mode < 1) {
      n <- (a / mode^2 + b / (1 - mode)^2) * mode * (1 - mode)
      shape_1 <- mode * n + 1
      shape_2 <- (1 - mode) * n + 1
    }
  }
  f_weight <- function(x) f_logp(x) - stats::dbeta(x, shape_1, shape_2, log=TRUE)
  current <- dated$fraction[1]
  for(step in seq_len(n_steps)) {
    proposed <- stats::rbeta(1, shape_1, shape_2)
    if(proposed <= 0 | proposed >= 1)
      next
    log_ratio <- f_weight(proposed) - f_weight(current)
    if(is.finite(log_ratio) && log(stats::runif(1)) < log_ratio)
      current <- proposed
  }
  state$dated$fraction <- rep(current, length(dated$fraction))
  state
}

#' The amount by which the cases recorded for each day exceed the true cases
#' whose onset was on it
#'
#' A fatal case dated by death which is recorded for a day but is not yet in
#' the reports is left out, so that an onset time taken to be fully reported
#' keeps any such case: in the exact form these are known, and otherwise no
#' more of the cases dated by death are taken to be in the reports than the
#' fatal cases which could have been shown as dead in the first report they
#' were in.
#'
#' @inheritParams deaths_delay_by_onset
#'
#' @return NULL if no fatal cases are dated by death; otherwise a tibble with
#' columns time_onset and cases_offset
deaths_cases_offset <- function(state) {
  if(is.null(state$dated))
    return(NULL)
  into <- state$dated$into
  if(state$dated$exact & !is.null(state$paths))
    into <- purrr::map_dbl(state$paths, ~sum(.$dated_new))
  else if(!is.null(state$paths))
    # no more of them can be in the reports than the fatal cases which could
    # have been shown as dead in the first report they were in
    into <- pmin(into, purrr::map_dbl(
      state$paths, ~sum(pmin(.$fatal_new, pmax(.$deaths_new, 0)))))
  dplyr::tibble(time_onset=state$onset_times,
                cases_offset=into - state$dated$out)
}
