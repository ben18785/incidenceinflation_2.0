# A joint move along the ridge between the fraction of fatal cases dated by
# death and the delay from onset to a death being recorded.
#
# A fatal case dated by death is recorded for its day of death and reported
# after the delay to its death being recorded, while one dated by onset is
# reported after the time from onset to death and then that delay. The
# reports therefore say much more about the mean delay from a listed onset to
# a death being shown than about how it divides between the two: more cases
# dated by death go with a longer delay to a death being recorded. The other
# updates change one of these at a time, given how many of each kind are in
# each report, which holds the others nearly fixed, so the chains move along
# the ridge slowly. This move changes them together: several fatal cases
# change kind at once, the fraction moves by as much, and the means of the
# delay to a death being recorded are scaled, with how many of each kind are
# in each report summed out; these are then drawn again given the new state.

#' Settings for the joint move of the fraction dated by death and the delay
#' to a death being recorded
#'
#' @param supplied NULL or FALSE (no such move), TRUE (the defaults) or a
#' named list of any of: 'n_moves' (the number proposed in each iteration,
#' default 2); 'max_relabel' (the most fatal cases which change kind in one
#' move, default 40); 'p_shift_only' (the probability that none change kind,
#' so that only the delay is scaled, default 0.25); 'slope' (the change in
#' the log of the mean delay for each unit of change in the fraction, default
#' 3) and 'sd' (the sd of a further normal change in that log, default 0.05)
#'
#' @return NULL or the named list of settings
ridge_settings <- function(supplied) {
  if(is.null(supplied) || identical(supplied, FALSE))
    return(NULL)
  ridge <- list(n_moves=2, max_relabel=40, p_shift_only=0.25, slope=3, sd=0.05)
  if(!isTRUE(supplied)) {
    unknown <- setdiff(names(supplied), names(ridge))
    if(length(unknown) > 0)
      stop(paste0("Unknown settings for the joint move: ",
                  paste(unknown, collapse=", "), "."))
    ridge[names(supplied)] <- supplied
  }
  if(ridge$n_moves < 1 | ridge$max_relabel < 1 | ridge$sd <= 0 |
     ridge$p_shift_only < 0 | ridge$p_shift_only >= 1)
    stop("The joint move needs at least one move and one case, a positive sd and a probability of a shift alone in [0, 1).")
  ridge
}

#' Log-probability of choosing a change of kind of one fatal case, as the
#' joint move chooses it
#'
#' A change to dated by death (sign 1) takes a fatal case dated by onset,
#' each equally likely, and a time from onset to death from its distribution
#' over the times which keep the onset within the days analysed. A change to
#' dated by onset (sign -1) takes a fatal case dated by death which died
#' within the days analysed, each equally likely.
#'
#' @inheritParams deaths_delay_by_onset
#' @param move a named list of 'day' (the day the case is recorded for), 'lag'
#' and 'sign', as for a change of type 'relabel' in \code{dated_apply_move}
#'
#' @return a log-probability (-Inf if the change cannot be chosen)
ridge_relabel_logq <- function(state, move) {
  dated <- state$dated
  if(move$sign == 1) {
    by_onset <- state$deaths_true - dated$out
    by_onset[1] <- 0  # no earlier day to move its onset to
    lags <- seq_len(min(dated$max_lag, move$day - 1))
    if(by_onset[move$day] <= 0 | !move$lag %in% lags)
      return(-Inf)
    log(by_onset[move$day]) - log(sum(by_onset)) +
      log(dated$lag_probabilities[move$lag]) -
      log(sum(dated$lag_probabilities[lags]))
  } else {
    k <- move$day - move$lag
    if(k < 1 || dated$lags[k, move$lag] <= 0)
      return(-Inf)
    log(dated$lags[k, move$lag]) - log(sum(dated$lags))
  }
}

#' Chooses a change of kind of one fatal case, as described in
#' \code{ridge_relabel_logq}
#'
#' @inheritParams deaths_delay_by_onset
#' @param sign 1 for a change to dated by death, -1 for the reverse
#'
#' @return NULL if there is no such case; otherwise a move as for
#' \code{dated_apply_move}
ridge_relabel_choose <- function(state, sign) {
  dated <- state$dated
  if(sign == 1) {
    by_onset <- state$deaths_true - dated$out
    by_onset[1] <- 0  # no earlier day to move its onset to
    if(sum(by_onset) <= 0)
      return(NULL)
    day <- sample.int(length(by_onset), 1, prob=by_onset)
    lags <- seq_len(min(dated$max_lag, day - 1))
    lag <- lags[sample.int(length(lags), 1, prob=dated$lag_probabilities[lags])]
  } else {
    if(sum(dated$lags) <= 0)
      return(NULL)
    cell <- sample.int(length(dated$lags), 1, prob=as.vector(dated$lags))
    k <- (cell - 1) %% nrow(dated$lags) + 1
    lag <- (cell - 1) %/% nrow(dated$lags) + 1
    day <- k + lag
  }
  list(type="relabel", day=day, lag=lag, lag_to=NA, sign=sign)
}

#' Log-probability of the parts of the model which the joint move alters,
#' with how many of each kind of fatal case are in each report summed out
#'
#' @inheritParams deaths_delay_by_onset
#' @inheritParams dated_logp
#' @inheritParams dated_renewal_logp
#' @param settings as returned by \code{deaths_settings}
#'
#' @return a named list of 'logp', 'forwards' (the forward filters of each
#' onset time) and 'forwards_uncorrected' (those of the uncorrected series)
ridge_logp <- function(state, cases_true, renewal, settings) {
  dated <- state$dated
  recording <- state$recording
  logp <- reporting_centre_logp_prior(
    recording$centre, recording$spread, settings$delay_prior,
    settings$delay_changepoints$spread_sd) +
    sum(reporting_segments_logp_prior(recording$parameters, recording$centre,
                                      recording$spread))
  if(!is.null(dated$fraction_prior)) {
    f <- dated$fraction[1]
    if(f <= 0 | f >= 1)
      return(list(logp=-Inf))
    logp <- logp + stats::dbeta(f, dated$fraction_prior$alpha,
                                dated$fraction_prior$beta, log=TRUE)
  }
  if(!is.finite(logp))
    return(list(logp=-Inf))
  # the fatal cases of each kind given the true cases, the cases recorded for
  # each day whose onset was before the first, and the renewal model
  logp <- logp + dated_logp(state, cases_true, days=integer(0))
  if(!is.finite(logp))
    return(list(logp=-Inf))
  for(t in which(dated$inflow > 0 | dated$inflow_weight > 0))
    logp <- logp + stats::dpois(dated$inflow[t],
                                dated$inflow_weight[t] * state$cfr[1] * dated$fraction[1],
                                log=TRUE)
  logp <- logp + dated_renewal_logp(cases_true, renewal)
  # the reports, with the split of each among the kinds summed out
  forwards <- lapply(seq_along(state$onset_times),
                     function(t) dated_exact_forward(state, cases_true, t))
  logp <- logp + sum(purrr::map_dbl(forwards, "logp"))
  forwards_uncorrected <- NULL
  sub <- uncorrected_state(state)
  if(!is.null(sub)) {
    forwards_uncorrected <- lapply(seq_along(sub$onset_times), function(k) {
      reports <- sub$reports[[k]]
      joint_path_forward(reports$cases_reported, reports$deaths_reported,
                         sub$deaths_true[k], sub$cases_true[k] - sub$deaths_true[k],
                         path_probabilities(sub$survivals[[k]], 0), sub$wrong[[k]])
    })
    logp <- logp + sum(purrr::map_dbl(forwards_uncorrected, "logp"))
  }
  list(logp=logp, forwards=forwards, forwards_uncorrected=forwards_uncorrected)
}

#' Proposes a joint move, as described in \code{deaths_ridge_step}
#'
#' @inheritParams deaths_delay_by_onset
#' @inheritParams dated_logp
#' @inheritParams deaths_paths_step
#' @param sign the direction: 1 (towards more fatal cases dated by death) or -1
#' @param m the number of fatal cases which change kind
#' @param noise the further change in the log of the mean delays
#' @param moves NULL (the changes of kind are chosen) or a list of the m
#' changes of kind to make, as for \code{dated_apply_move}
#'
#' @return NULL if the move is impossible; otherwise a named list of the
#' proposed 'state' (with the survivals given the new delays) and
#' 'cases_true', 'moves' (the changes of kind made), 'log_q' (the log of the
#' probability of choosing the reverse changes less that of choosing these)
#' and 'log_jacobian'
ridge_propose <- function(state, cases_true, reporting_by_onset,
                          revision_parameters, sign, m, noise, moves=NULL) {
  ridge <- state$dated$ridge
  n_fatal <- sum(state$deaths_true)
  if(n_fatal <= 0)
    return(NULL)
  log_q <- 0
  made <- vector("list", m)
  for(j in seq_len(m)) {
    move <- if(is.null(moves)) ridge_relabel_choose(state, sign) else moves[[j]]
    if(is.null(move) || move$sign != sign)
      return(NULL)
    log_q <- log_q - ridge_relabel_logq(state, move)
    applied <- dated_apply_move(state, cases_true, move)
    if(is.null(applied) || !is.finite(log_q))
      return(NULL)
    state <- applied$state
    cases_true <- applied$cases_true
    # the reverse move undoes this change from the state it leads to
    log_q <- log_q + ridge_relabel_logq(
      state, list(day=move$day, lag=move$lag, sign=-sign))
    made[[j]] <- move
  }
  if(!is.finite(log_q))
    return(NULL)
  if(!is.null(state$dated$fraction_prior))
    state$dated$fraction <- state$dated$fraction + sign * m / n_fatal
  log_scale <- ridge$slope * sign * m / n_fatal + noise
  state$recording$parameters$mean <- state$recording$parameters$mean * exp(log_scale)
  state$recording$centre$mean <- state$recording$centre$mean * exp(log_scale)
  list(state=deaths_refresh_survivals(state, reporting_by_onset, revision_parameters),
       cases_true=cases_true, moves=made, log_q=log_q,
       log_jacobian=(nrow(state$recording$parameters) + 1) * log_scale)
}

#' Moves the fraction dated by death and the delay to a death being recorded
#' together
#'
#' Each move chooses a direction s (1 or -1, equally likely) and a number m
#' of fatal cases (none with probability p_shift_only, otherwise from one to
#' max_relabel, equally likely). Then m fatal cases change kind one after
#' another, to dated by death if s is 1 and the reverse otherwise, each chosen
#' as in \code{ridge_relabel_logq}; if the fraction is estimated it changes by
#' s m / N, where N is the number of fatal cases, which is what it would be
#' were it the share of them dated by death; and the means of the delay to a
#' death being recorded of every segment, and their centre, are multiplied
#' by exp(slope s m / N + e), with e normal with mean zero and sd 'sd'. The
#' reverse move has direction -s, the same m and -e, and undoes the changes
#' of kind in reverse order, so the ratio of the probabilities of proposing
#' the two is the product of those of choosing each change of kind; the
#' Jacobian of the scaling is exp((n + 1) times the change in log), for n
#' segments. The move is accepted with the Metropolis-Hastings probability,
#' for the model with how many of each kind of fatal case are in each report
#' summed out (\code{ridge_logp}); if it is accepted these are then drawn
#' given the new state.
#'
#' @inheritParams deaths_delay_by_onset
#' @inheritParams dated_logp
#' @inheritParams dated_renewal_logp
#' @inheritParams deaths_paths_step
#' @inheritParams ridge_logp
#'
#' @return a named list of the updated 'state' and 'cases_true' and
#' 'accepted', the number of moves accepted
deaths_ridge_step <- function(state, cases_true, renewal, reporting_by_onset,
                              revision_parameters, settings) {
  ridge <- state$dated$ridge
  if(is.null(ridge) || !is.null(settings$delay_changepoints$drift))
    return(list(state=state, cases_true=cases_true, accepted=0))
  if(isTRUE(state$dated$collapse_fraction))
    stop("The joint move cannot be combined with the fraction integrated out.")
  current <- ridge_logp(state, cases_true, renewal, settings)
  accepted <- 0
  for(i in seq_len(ridge$n_moves)) {
    sign <- sample(c(-1, 1), 1)
    m <- if(stats::runif(1) < ridge$p_shift_only) 0
      else sample.int(ridge$max_relabel, 1)
    proposed <- ridge_propose(state, cases_true, reporting_by_onset,
                              revision_parameters, sign, m,
                              stats::rnorm(1, 0, ridge$sd))
    if(is.null(proposed))
      next
    candidate <- ridge_logp(proposed$state, proposed$cases_true, renewal, settings)
    log_ratio <- candidate$logp - current$logp + proposed$log_jacobian +
      proposed$log_q
    if(is.finite(log_ratio) && log(stats::runif(1)) < log_ratio) {
      state <- proposed$state
      cases_true <- proposed$cases_true
      current <- candidate
      accepted <- accepted + 1
    }
  }
  # how many of each kind are in each report, given the new state
  if(accepted > 0) {
    state$paths <- lapply(current$forwards, joint_path_exact_redraw)
    if(!is.null(current$forwards_uncorrected))
      state$uncorrected$paths <- lapply(current$forwards_uncorrected,
                                        joint_path_backward)
  }
  list(state=state, cases_true=cases_true, accepted=accepted)
}
