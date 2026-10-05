#' Calculates observation process log probability density for a pair of days
#'
#' The probability of observations is a function of the number of cases remaining
#' to be reported:
#' \deqn{I_remaining = I_true - I_day_1}
#' and the number of cases observed between day_1 and day_2:
#' \deqn{I_obs = I_day_2 - I_day_1}
#' It is given by:
#' \deqn{binomial(I_obs|I_remaining, p_detect)}
#' where p_detect is the probability a thus undetected case is detected between
#' day_1 and day_2.
#'
#' @param cases_day_2 number of cases arising on day_onset observed by day_2
#' @param cases_day_1 number of cases arising on day_onset observed by day_1 < day_2
#' @param cases_true true case count(s) originating on day_onset: can be a single value or vector
#' @inheritParams detected_after_unobserved_prob
#' @inheritParams observed_cases_single
#'
#' @return a log-probability or vector of log-probabilities
observation_process_single_logp <- function(cases_true, cases_day_2, cases_day_1, day_2, day_1,
                                            day_onset, reporting_parameters){
  cases_observed <- cases_day_2 - cases_day_1
  cases_remaining <- cases_true - cases_day_1
  p_detect <- detected_after_unobserved_prob(day_2, day_1, day_onset,
                                             reporting_parameters)
  stats::dbinom(cases_observed, cases_remaining, p_detect, log = T)
}

#' Enumerates possible changes in the count of erroneous records between two
#' reporting days
#'
#' Between two reporting days, erroneous records present on the first day can
#' be withdrawn and new erroneous records can be entered. This enumerates all
#' combinations of these which keep the count of erroneous records between zero
#' and max_erroneous.
#'
#' @param max_erroneous maximum number of erroneous records that can be held
#' for a single onset day at any one time
#'
#' @return a data frame with columns: 'erroneous_previous', the count on the
#' first day; 'withdrawn', the number of these withdrawn; 'entered', the number
#' of new erroneous records; and 'erroneous_current', the count on the second day
erroneous_transitions <- function(max_erroneous) {
  if(max_erroneous < 0)
    stop("max_erroneous must be non-negative.")
  df <- expand.grid(erroneous_current=0:max_erroneous,
                    withdrawn=0:max_erroneous,
                    erroneous_previous=0:max_erroneous)
  df$entered <- df$erroneous_current - df$erroneous_previous + df$withdrawn
  df <- df[df$withdrawn <= df$erroneous_previous & df$entered >= 0, ]
  rownames(df) <- NULL
  df[, c("erroneous_previous", "withdrawn", "entered", "erroneous_current")]
}

#' Binomial log-probability which is -Inf (rather than NaN) outside of support
#'
#' @param x vector of counts
#' @param size vector of numbers of trials, which can be negative
#' @param prob probability of success
#'
#' @return a vector of log-probabilities
dbinom_log_safe <- function(x, size, prob) {
  logp <- rep(-Inf, length(x))
  is_valid <- size >= 0 & x >= 0 & x <= size
  logp[is_valid] <- stats::dbinom(x[is_valid], size[is_valid], prob, log=TRUE)
  logp
}

#' Calculates log-probability of count of erroneous records on the first
#' reporting day
#'
#' If cases are first reported on the onset day, no erroneous records have yet
#' been entered. Otherwise, erroneous records are assumed to have been entered
#' and withdrawn on each day between the onset day and the first reporting day
#' in the same way as between later reporting days. The count of erroneous
#' records on the first reporting day then has prior distribution:
#' \deqn{E \sim Poisson(eta * I_true * \sum_k (S(k - 1) - S(k)) * (1 - rho)^(day_first - k))}
#' where the sum is over days after the onset day up to the first reporting day.
#' Independently, the genuine cases reported by then are:
#' \deqn{D \sim binomial(I_true, p_detect)}
#' This calculates the distribution of E conditional on the cases first
#' reported, C = D + E, with E truncated at max_erroneous.
#'
#' @param day_first first day on which cases were reported
#' @param cases_first count of cases reported on day_first
#' @inheritParams observation_process_revisions_logp
#'
#' @return a matrix of log-probabilities with a row for each element of
#' cases_true and a column for each count of erroneous records from zero to
#' max_erroneous. If cases_first is impossible for an element of cases_true,
#' its row is -Inf.
erroneous_initial_logp <- function(day_first, cases_first, cases_true, day_onset,
                                   reporting_parameters, revision_parameters,
                                   max_erroneous=3) {
  n_cases <- length(cases_true)
  log_alpha <- matrix(-Inf, nrow=n_cases, ncol=max_erroneous + 1)

  n_days <- floor(day_first - day_onset)
  if(n_days < 1 | revision_parameters$eta == 0) {
    log_alpha[, 1] <- 0
    return(log_alpha)
  }

  # expected erroneous records entered each day which remain on day_first
  days <- day_onset + seq_len(n_days)
  p_undetected <- 1 - pgamma_mean_sd(c(day_onset, days) + 0.5 - day_onset,
                                     reporting_parameters$mean,
                                     reporting_parameters$sd)
  p_newly_detected <- -diff(p_undetected)
  p_remaining <- (1 - revision_parameters$rho)^(day_first - days)
  mean_erroneous <- revision_parameters$eta * cases_true *
    sum(p_newly_detected * p_remaining)

  p_detect <- detected_after_unobserved_prob(day_first, day_onset, day_onset,
                                             reporting_parameters)
  for(e in 0:max_erroneous)
    log_alpha[, e + 1] <- stats::dpois(e, mean_erroneous, log=TRUE) +
      dbinom_log_safe(rep(cases_first - e, n_cases), cases_true, p_detect)

  # condition on cases first reported
  log_normaliser <- matrixStats::rowLogSumExps(log_alpha)
  log_alpha <- log_alpha - log_normaliser
  log_alpha[!is.finite(log_normaliser), ] <- -Inf
  log_alpha
}

#' Calculates observation process log probability allowing reported case counts
#' to be revised downwards
#'
#' The reported count for an onset day is assumed to be the sum of genuine
#' cases reported so far and erroneous records not yet withdrawn:
#' \deqn{C = D + E}
#' Between a pair of reporting days, day_1 < day_2, there are three independent changes:
#' \deqn{d \sim binomial(I_true - D_day_1, p_detect)}
#' genuine cases are newly reported, exactly as in \code{observation_process_single_logp};
#' \deqn{n \sim Poisson(eta * I_true * (S(day_1) - S(day_2)))}
#' erroneous records are entered, where S is the probability a case remains undetected; and
#' \deqn{m \sim binomial(E_day_1, 1 - (1 - rho)^(day_2 - day_1))}
#' erroneous records are withdrawn. The change in reported cases is then:
#' \deqn{C_day_2 - C_day_1 = d + n - m}
#' E is unobserved so is summed over using a forward recursion, with E
#' truncated at max_erroneous. Its distribution on the first reporting day is
#' given by \code{erroneous_initial_logp}.
#'
#' @param observation_df a two-column data frame with columns 'time_reported' and 'cases_reported'
#' @param cases_true true case count(s) originating on day_onset: can be a single value or vector
#' @inheritParams detected_after_unobserved_prob
#' @param revision_parameters named list of 'eta' (>= 0), the expected number of
#' erroneous records entered per genuine case reported; and 'rho' (between 0 and 1),
#' the daily probability an erroneous record is withdrawn. If NULL (the default
#' where this is optional), reported cases cannot be revised downwards.
#' @inheritParams erroneous_transitions
#'
#' @return a log-probability or vector of log-probabilities
observation_process_revisions_logp <- function(observation_df, cases_true,
                                               day_onset, reporting_parameters,
                                               revision_parameters,
                                               max_erroneous=3) {
  observation_days <- observation_df$time_reported
  observation_cases <- observation_df$cases_reported
  if(any(diff(observation_days) <= 0))
    stop("days in observation matrix must be increasing.")

  eta <- revision_parameters$eta
  rho <- revision_parameters$rho
  n_cases <- length(cases_true)
  if(eta < 0 | rho <= 0 | rho >= 1)
    return(rep(-Inf, n_cases))

  transitions <- erroneous_transitions(max_erroneous)
  n_transitions <- nrow(transitions)
  erroneous_previous <- rep(transitions$erroneous_previous, each=n_cases)
  entered <- rep(transitions$entered, each=n_cases)
  change_erroneous <- rep(transitions$erroneous_current - transitions$erroneous_previous,
                          each=n_cases)
  cases_true_long <- rep(cases_true, n_transitions)

  # log-probability of observations so far and each count of erroneous records
  log_alpha <- erroneous_initial_logp(observation_days[1], observation_cases[1],
                                      cases_true, day_onset,
                                      reporting_parameters, revision_parameters,
                                      max_erroneous)

  num_periods <- nrow(observation_df) - 1
  for(i in seq_len(num_periods)) {
    day_1 <- observation_days[i]
    day_2 <- observation_days[i + 1]
    p_detect <- detected_after_unobserved_prob(day_2, day_1, day_onset,
                                               reporting_parameters)
    p_newly_detected <- undetected_prob(day_1 + 0.5, day_onset, reporting_parameters) -
      undetected_prob(day_2 + 0.5, day_onset, reporting_parameters)
    p_withdrawn <- 1 - (1 - rho)^(day_2 - day_1)
    change_cases <- observation_cases[i + 1] - observation_cases[i]

    logp_withdrawn <- stats::dbinom(transitions$withdrawn,
                                    transitions$erroneous_previous,
                                    p_withdrawn, log=TRUE)
    logp_entered <- stats::dpois(entered,
                                 eta * cases_true_long * p_newly_detected,
                                 log=TRUE)
    # genuine cases remaining are those not among the non-erroneous records
    logp_genuine <- dbinom_log_safe(
      change_cases - change_erroneous,
      cases_true_long - observation_cases[i] + erroneous_previous,
      p_detect)

    log_terms <- log_alpha[, transitions$erroneous_previous + 1, drop=FALSE] +
      matrix(logp_entered + logp_genuine, nrow=n_cases) +
      matrix(logp_withdrawn, nrow=n_cases, ncol=n_transitions, byrow=TRUE)

    for(e in 0:max_erroneous)
      log_alpha[, e + 1] <- matrixStats::rowLogSumExps(
        log_terms[, transitions$erroneous_current == e, drop=FALSE])
  }
  matrixStats::rowLogSumExps(log_alpha)
}

#' Calculates observation process log probability density
#'
#' The probability of observations is a function of the number of cases remaining
#' to be reported:
#' \deqn{I_remaining = I_true - I_day_1}
#' and the number of cases observed between day_1 and day_2:
#' \deqn{I_obs = I_day_2 - I_day_1}
#' It is given by:
#' \deqn{binomial(I_obs|I_remaining, p_detect)}
#' where p_detect is the probability a thus undetected case is detected between
#' day_1 and day_2.
#' This function calculates the combined probability of all observed cases pertaining
#' to a particular onset date.
#'
#' If revision_parameters is supplied, reported case counts are allowed to be
#' revised downwards and the probability is instead calculated using
#' \code{observation_process_revisions_logp}.
#'
#' @param observation_df a two-column data frame with columns 'time_reported' and 'cases_reported'
#' @inheritParams observation_process_single_logp
#' @inheritParams observation_process_revisions_logp
#'
#' @return a log-probability or vector of log-probabilities
observation_process_logp <- function(observation_df, cases_true,
                                     day_onset, reporting_parameters,
                                     revision_parameters=NULL,
                                     max_erroneous=3){
  if(!is.null(revision_parameters))
    return(observation_process_revisions_logp(
      observation_df=observation_df,
      cases_true=cases_true,
      day_onset=day_onset,
      reporting_parameters=reporting_parameters,
      revision_parameters=revision_parameters,
      max_erroneous=max_erroneous))

  observation_days <- observation_df$time_reported
  observation_cases <- observation_df$cases_reported
  change_days <- diff(observation_days)
  if(any(change_days <= 0))
    stop("days in observation matrix must be increasing.")
  change_cases <- diff(observation_cases)
  if(any(change_cases < 0))
    stop("reported cases in observation matrix must be increasing.")

  num_periods <- nrow(observation_df) - 1
  logp <- 0
  for(i in 1:num_periods) {
    logp_day <- observation_process_single_logp(cases_true=cases_true,
                                         cases_day_2=observation_cases[i + 1],
                                         cases_day_1=observation_cases[i],
                                         day_2=observation_days[i + 1],
                                         day_1=observation_days[i],
                                         day_onset=day_onset,
                                         reporting_parameters=reporting_parameters)
    logp <- logp + logp_day
  }
  logp
}

#' Determines the probability of a given number of cases given a history of them
#'
#' If a Poisson renwal is being used, this calculates:
#' \deqn{Pois(cases_true_t|Rt * \sum_tau=1^t_max w_t cases_true_t-tau)}
#' or the following if a negative-binomial model is used:
#' \deqn{NB(cases_true_t|Rt * \sum_tau=1^t_max w_t cases_true_t-tau, kappa)}
#'
#' @inheritParams observed_cases_single
#' @inheritParams true_cases_single
#' @inheritParams gamma_discrete_pmf
#' @param kappa overdispersion parameter
#' @param is_negative_binomial if negative-binomial renewal model specified (defaults to FALSE)
#'
#' @return a log-probability or vector of log-probabilities
state_process_logp <- function(cases_true, cases_history, Rt, serial_parameters,
                               kappa=NULL, is_negative_binomial=FALSE){
  t_max <- length(cases_history)
  w <- weights_series(t_max, serial_parameters)
  mean_cases <- expected_cases(Rt=Rt, weights=w, cases_history=cases_history)

  if(!is_negative_binomial)
    stats::dpois(cases_true, lambda=mean_cases, log=TRUE)
  else
    stats::dnbinom(cases_true, mu=mean_cases, size=kappa, log=TRUE)
}

#' Calculates the log-probability of case counts arising after a given onset day
#' as a function of the true case count on that day
#'
#' Cases arising on a given onset day form part of the history of cases for
#' all later onset days. This calculates the contribution of these later days:
#' \deqn{\prod_tau=t+1^t+F p(cases_true_tau|cases_true_tau-1, cases_true_tau-2, ..., Rt_tau, serial_params)}
#' where each term is Poisson or negative binomial, as in \code{state_process_logp}.
#'
#' @param cases_true true case count(s) originating on day_onset: can be a single value or vector
#' @param cases_history a vector containing history of cases arranged from recent to past
#' (all arising before day_onset)
#' @param cases_future a vector of true case counts arising on the onset days
#' immediately following day_onset, arranged from the earliest to the latest
#' @param Rt_future a vector of Rt values for the onset days corresponding to cases_future
#' @inheritParams state_process_logp
#'
#' @return a log-probability or vector of log-probabilities
future_cases_logp <- function(cases_true, cases_history, cases_future, Rt_future,
                              serial_parameters,
                              kappa=NULL, is_negative_binomial=FALSE) {
  n_future <- length(cases_future)
  if(length(Rt_future) != n_future)
    stop("cases_future and Rt_future must be same length.")

  logp <- rep(0, length(cases_true))
  if(n_future == 0)
    return(logp)

  n_history <- length(cases_history)
  w_longest <- weights_series(n_history + n_future, serial_parameters)

  for(k in 1:n_future) {
    # same (truncated then normalised) weights as used in state_process_logp
    w <- w_longest[1:(n_history + k)]
    w <- w / sum(w)

    # history for this later day from recent to past excluding day_onset,
    # which is k days before it
    cases_other <- c(rev(cases_future[seq_len(k - 1)]), cases_history)
    mean_cases <- Rt_future[k] * (w[k] * cases_true + sum(w[-k] * cases_other))

    if(!is_negative_binomial)
      logp <- logp + stats::dpois(cases_future[k], lambda=mean_cases, log=TRUE)
    else
      logp <- logp + stats::dnbinom(cases_future[k], mu=mean_cases, size=kappa, log=TRUE)
  }
  logp
}

#' Calculates (unnormalised) log-probability of a true case count pertaining to cases arising
#' on a particular onset day
#'
#' There are two components to this:
#' \deqn{p(cases_true_t|data, Rt, reporting_params, serial_params) \propto p(data|cases_true, reporting_params)
#'  p(cases_true_t|cases_true_t_1, cases_true_t_2, ..., Rt, serial_params)}
#' If cases_future is supplied, a third component is included which accounts
#' for the case counts arising on later onset days, since these depend on cases_true_t:
#' \deqn{\prod_tau>t p(cases_true_tau|cases_true_tau-1, cases_true_tau-2, ..., Rt_tau, serial_params)}
#'
#' @inheritParams state_process_logp
#' @inheritParams observation_process_logp
#' @inheritParams future_cases_logp
#'
#' @return an (unnormalised) log-probability or vector of such log-probabilities
conditional_cases_logp <- function(cases_true, observation_df, cases_history,
                                   Rt, day_onset, serial_parameters, reporting_parameters,
                                   kappa=NULL, is_negative_binomial=FALSE,
                                   cases_future=NULL, Rt_future=NULL,
                                   revision_parameters=NULL, max_erroneous=3) {

  logp_observation <- 0
  if(nrow(observation_df) > 1) # handling cases arising today which were reported today
    logp_observation <- observation_process_logp(observation_df=observation_df,
                                                 cases_true=cases_true,
                                                 day_onset=day_onset,
                                                 reporting_parameters=reporting_parameters,
                                                 revision_parameters=revision_parameters,
                                                 max_erroneous=max_erroneous)
  logp_state <- state_process_logp(cases_true=cases_true,
                                   cases_history=cases_history,
                                   Rt=Rt,
                                   serial_parameters=serial_parameters,
                                   kappa=kappa,
                                   is_negative_binomial=is_negative_binomial)
  logp_future <- 0
  if(!is.null(cases_future))
    logp_future <- future_cases_logp(cases_true=cases_true,
                                     cases_history=cases_history,
                                     cases_future=cases_future,
                                     Rt_future=Rt_future,
                                     serial_parameters=serial_parameters,
                                     kappa=kappa,
                                     is_negative_binomial=is_negative_binomial)
  logp_observation + logp_state + logp_future
}

#' Calculates (unnormalised) log-probability of true case counts pertaining to
#' cases arising on a particular onset day when reported case counts can be
#' revised downwards
#'
#' This is the same quantity as calculated by \code{conditional_cases_logp} but,
#' since the observation process is more costly to evaluate when revisions are
#' allowed, this is only evaluated for case counts whose probability is
#' non-negligible. Since the observation process log-probability cannot exceed
#' zero, the remaining components bound the log-probability from above; any case
#' count whose bound is more than log_negligible below the log-probability of
#' a pilot case count is assigned a log-probability of -Inf.
#'
#' @inheritParams conditional_cases_logp
#' @param log_negligible difference in log-probability from a pilot case count
#' beyond which case counts are assigned zero probability
#'
#' @return a vector of (unnormalised) log-probabilities
conditional_cases_revisions_logp <- function(cases_true, observation_df, cases_history,
                                             Rt, day_onset, serial_parameters, reporting_parameters,
                                             revision_parameters,
                                             kappa=NULL, is_negative_binomial=FALSE,
                                             cases_future=NULL, Rt_future=NULL,
                                             max_erroneous=3,
                                             log_negligible=50) {

  # a single observation contributes nothing to the log-probability, so this
  # gives the components due to the renewal process only
  logps_upper <- conditional_cases_logp(cases_true, observation_df[1, ], cases_history,
                                        Rt, day_onset, serial_parameters, reporting_parameters,
                                        kappa, is_negative_binomial,
                                        cases_future, Rt_future)
  if(nrow(observation_df) <= 1)
    return(logps_upper)

  f_observation <- function(cases) {
    observation_process_logp(observation_df=observation_df,
                             cases_true=cases,
                             day_onset=day_onset,
                             reporting_parameters=reporting_parameters,
                             revision_parameters=revision_parameters,
                             max_erroneous=max_erroneous)
  }

  cases_last_reported <- dplyr::last(observation_df$cases_reported)
  index_pilot <- unique(c(which.max(logps_upper),
                          which.min(abs(cases_true - cases_last_reported))))
  logps_pilot <- logps_upper[index_pilot] + f_observation(cases_true[index_pilot])

  is_evaluated <- logps_upper >= max(logps_pilot) - log_negligible
  logps <- rep(-Inf, length(cases_true))
  logps[is_evaluated] <- logps_upper[is_evaluated] +
    f_observation(cases_true[is_evaluated])
  logps
}

#' Observation probability across all onset times
#'
#' @param snapshot_with_true_cases_df a tibble with four columns:
#' time_onset, time_reported, cases_reported, cases_true
#' @param reporting_parameters a tibble with columns: 'reporting_piece_index', 'mean', 'sd'
#' @inheritParams observation_process_revisions_logp
#'
#' @return a log probability
#' @importFrom rlang .data
observation_process_all_times_logp <- function(
  snapshot_with_true_cases_df,
  reporting_parameters,
  revision_parameters=NULL,
  max_erroneous=3){

  if(methods::is(reporting_parameters, "list"))
    stop("Reporting parameters should be a tibble not a list.")

  if(!"reporting_piece_index" %in% colnames(reporting_parameters))
    stop("reporting_parameters must contain a column: 'reporting_piece_index'.")

  reporting_parameters <- snapshot_with_true_cases_df %>%
    dplyr::left_join(reporting_parameters, by="reporting_piece_index") %>%
    dplyr::select("time_onset", "mean", "sd") %>%
    unique()

  logp <- 0
  onset_times <- unique(snapshot_with_true_cases_df$time_onset)

  # check that only one true case per onset time
  n_onset <- length(onset_times)
  test_df <- snapshot_with_true_cases_df %>%
    dplyr::select("time_onset", "cases_true") %>%
    unique()
  if(n_onset != nrow(test_df))
    stop("There must be only one true case measurement per each onset time.")

  mu <- reporting_parameters$mean
  sd <- reporting_parameters$sd
  if(sum(mu < 0) > 0 | sum(sd < 0) > 0)
    return(-Inf)

  for(i in seq_along(onset_times)) {
    onset_time <- onset_times[i]
    short_df <- snapshot_with_true_cases_df %>%
      dplyr::filter(.data$time_onset==onset_time)
    cases_true <- short_df$cases_true[1]
    reporting_parameters_current <- list(mean=reporting_parameters$mean[i],
                                         sd=reporting_parameters$sd[i])
    if(nrow(short_df) > 1) {# handling cases arising today which were reported today
      logp_new <- observation_process_logp(
        short_df,
        cases_true=cases_true,
        day_onset=onset_time,
        reporting_parameters=reporting_parameters_current,
        revision_parameters=revision_parameters,
        max_erroneous=max_erroneous)
      logp <- logp + logp_new
    }
  }
  logp
}


#' Calculates the overall log-probability of cases conditional on their history
#' and Rt values, assuming a negative-binomial renewal model
#'
#' Specifically, this assumes the probability is of the form:
#' \deqn{\prod_tau=1^T NB(cases_true_tau|Rt * \sum_s=1^tau-1 w_s cases_true_tau-s, kappa)}
#'
#' @param kappa an overdispersion parameter value (if <= 0, returns -Inf)
#' @param cases_history_rt_df a tibble with three columns: 'time_onset',
#' 'cases_true' and 'Rt'
#' @inheritParams state_process_logp
#'
#' @return a log-probability value
state_process_nb_logp_all_onsets <- function(
    kappa, cases_history_rt_df, serial_parameters) {

  if(kappa <= 0)
    return(-Inf)

  onset_times <- unique(cases_history_rt_df$time_onset)

  log_p <- 0
  for(i in seq_along(onset_times)) {

    onset_time <- onset_times[i]
    if(onset_time > 1) {

      cases_true <- cases_history_rt_df$cases_true[i]
      pre_observation_df <- cases_history_rt_df %>%
        dplyr::filter(.data$time_onset < onset_time) %>%
        dplyr::arrange(dplyr::desc(.data$time_onset))
      cases_history <- pre_observation_df$cases_true
      Rt <- cases_history_rt_df$Rt[i]

      logp_single_onset_time <- state_process_logp(
        cases_true=cases_true,
        cases_history=cases_history,
        Rt=Rt,
        serial_parameters=serial_parameters,
        kappa=kappa,
        is_negative_binomial=TRUE)

      log_p <- log_p + logp_single_onset_time
    }
  }
  log_p
}

