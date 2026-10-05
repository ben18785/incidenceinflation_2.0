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

#' Probability that a genuine case is misdated
#'
#' @param revision_parameters named list which may contain an element 'pi'
#'
#' @return the misdating probability, which is zero if not specified
misdating_prob <- function(revision_parameters) {
  if(is.null(revision_parameters$pi)) 0 else revision_parameters$pi
}

#' Probability that a case has not yet been recorded under its onset day
#'
#' Without misdating this is the probability that a case remains undetected,
#' S. With misdating, a detected case may be recorded under a nearby onset day
#' until it is corrected, which happens with probability rho each day after it
#' is entered, so that the probability is:
#' \deqn{S(day) + pi * \sum_k g(k) * (S(k - 1) - S(k)) * (1 - rho)^(day - k)}
#' where the sum is over days after the onset day up to day and g(k) is the
#' proportion of shifts possible for a case detected on day k: a record cannot
#' be given an onset day on or after the day on which it is entered.
#'
#' @param days vector of days (each a whole number of days after day_onset)
#' @inheritParams detected_after_unobserved_prob
#' @param revision_parameters named list of 'eta', 'rho' and optionally 'pi',
#' the probability that a genuine case is misdated, along with 'shifts' (a
#' vector of non-zero integers) and 'shift_probs' (their probabilities)
#'
#' @return a vector of probabilities, evaluated at the end of each day
unrecorded_prob <- function(days, day_onset, reporting_parameters,
                            revision_parameters) {
  f_undetected <- function(day)
    1 - pgamma_mean_sd(day + 0.5 - day_onset,
                       reporting_parameters$mean, reporting_parameters$sd)
  p_misdated <- misdating_prob(revision_parameters)
  n_days <- round(max(days) - day_onset)
  if(p_misdated == 0 | n_days < 1)
    return(f_undetected(days))

  days_after <- seq_len(n_days)
  p_undetected <- f_undetected(day_onset + c(0, days_after))
  p_newly_detected <- -diff(p_undetected)
  p_shift_possible <- purrr::map_dbl(
    days_after,
    ~sum(revision_parameters$shift_probs[revision_parameters$shifts <= . - 1]))
  p_awaiting <- numeric(n_days + 1)
  for(k in days_after)
    p_awaiting[k + 1] <- p_awaiting[k] * (1 - revision_parameters$rho) +
    p_shift_possible[k] * p_newly_detected[k]
  p_unrecorded <- p_undetected + p_misdated * p_awaiting
  p_unrecorded[round(days - day_onset) + 1]
}

#' Expected misdated records entered for an onset day per case arising on a
#' neighbouring onset day
#'
#' A case with onset day (day_onset - shift) which is detected and misdated by
#' shift is recorded under day_onset. For each interval between reporting days,
#' the expected number of such records entered per case is:
#' \deqn{pi * h(shift) * (S(day_1) - S(day_2))}
#' where S is the probability a case arising on the neighbouring day remains
#' undetected. Also calculated is the expected number present on the first
#' reporting day: those entered on days after day_onset which have not yet
#' been corrected.
#'
#' @param observation_days vector of days on which cases for day_onset were reported
#' @param day_onset the onset day under which misdated records are entered
#' @param shift difference between day_onset and the true onset day of the cases
#' @inheritParams unrecorded_prob
#'
#' @return a named list of 'initial', the expected number present on the first
#' reporting day; and 'periods', a vector of the expected numbers entered
#' between successive reporting days
misdated_entry_unit <- function(observation_days, day_onset, shift,
                                reporting_parameters, revision_parameters) {
  p_shift <- revision_parameters$shift_probs[revision_parameters$shifts == shift]
  day_true <- day_onset - shift
  f_undetected <- function(day)
    1 - pgamma_mean_sd(pmax(day + 0.5 - day_true, 0),
                       reporting_parameters$mean, reporting_parameters$sd)
  rate <- misdating_prob(revision_parameters) * p_shift

  initial <- 0
  n_days <- floor(observation_days[1] - day_onset)
  if(n_days >= 1) {
    days <- day_onset + seq_len(n_days)
    p_newly_detected <- f_undetected(days - 1) - f_undetected(days)
    initial <- rate * sum(
      p_newly_detected * (1 - revision_parameters$rho)^(observation_days[1] - days))
  }
  list(initial=initial,
       periods=rate * -diff(f_undetected(observation_days)))
}

#' Expected misdated records entered for an onset day from all neighbouring
#' onset days
#'
#' @inheritParams misdated_entry_unit
#' @param cases_by_onset a vector of true case counts named by onset time
#' @param exclude an onset time whose cases are not to be included
#'
#' @return a named list of 'initial' and 'periods' as returned by
#' \code{misdated_entry_unit}, summed across neighbouring onset days
misdated_entries <- function(observation_days, day_onset, cases_by_onset,
                             reporting_parameters, revision_parameters,
                             exclude=NULL) {
  entries <- list(initial=0, periods=rep(0, length(observation_days) - 1))
  for(shift in revision_parameters$shifts) {
    day_true <- day_onset - shift
    cases <- cases_by_onset[as.character(day_true)]
    if(is.na(cases) | identical(as.numeric(day_true), as.numeric(exclude)))
      next
    unit <- misdated_entry_unit(observation_days, day_onset, shift,
                                reporting_parameters, revision_parameters)
    entries$initial <- entries$initial + cases[[1]] * unit$initial
    entries$periods <- entries$periods + cases[[1]] * unit$periods
  }
  entries
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
#' @param mean_erroneous_extra expected number of misdated records present on
#' the first reporting day (a single value or one for each element of cases_true)
#' @param phi_batch intensity of a batch whose footprint contains the onset
#' day and which has not been detected by the first reporting day (zero if none).
#' The count of its records, X, on the first reporting day has prior distribution:
#' \deqn{X \sim Poisson(phi * I_true * (S(day_onset) - S(day_first)))}
#' and the distribution calculated is then that of (E, X) conditional on the
#' cases first reported, C = D + E + X.
#' @param max_batch maximum count of batch records (zero if there is no batch)
#'
#' @return a matrix of log-probabilities with a row for each element of
#' cases_true and a column for each count of erroneous records from zero to
#' max_erroneous. If max_batch exceeds zero, there is a column for each
#' combination of counts, with (e, x) in column 1 + e + (max_erroneous + 1) * x.
#' If cases_first is impossible for an element of cases_true, its row is -Inf.
erroneous_initial_logp <- function(day_first, cases_first, cases_true, day_onset,
                                   reporting_parameters, revision_parameters,
                                   max_erroneous=3, mean_erroneous_extra=0,
                                   phi_batch=0, max_batch=0) {
  n_cases <- length(cases_true)
  n_erroneous <- max_erroneous + 1
  log_alpha <- matrix(-Inf, nrow=n_cases, ncol=n_erroneous * (max_batch + 1))

  n_days <- floor(day_first - day_onset)
  if(n_days < 1 | (revision_parameters$eta == 0 & all(mean_erroneous_extra == 0) &
                   phi_batch == 0)) {
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
    sum(p_newly_detected * p_remaining) + mean_erroneous_extra

  p_unrecorded <- unrecorded_prob(c(day_onset, day_first), day_onset,
                                  reporting_parameters, revision_parameters)
  p_detect <- 1 - p_unrecorded[2] / p_unrecorded[1]
  mean_batch <- phi_batch * cases_true * (p_undetected[1] - p_undetected[n_days + 1])
  for(x in 0:max_batch) {
    logp_batch <- if(max_batch > 0) stats::dpois(x, mean_batch, log=TRUE) else 0
    for(e in 0:max_erroneous)
      log_alpha[, 1 + e + n_erroneous * x] <-
        stats::dpois(e, mean_erroneous, log=TRUE) + logp_batch +
        dbinom_log_safe(rep(cases_first - e - x, n_cases), cases_true, p_detect)
  }

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
#' where this is optional), reported cases cannot be revised downwards. The
#' list can also contain 'pi', 'shifts' and 'shift_probs' if cases can be
#' misdated (see \code{unrecorded_prob}); and 'batches', a tibble of batches of
#' erroneous records as returned by \code{sample_batches}, along with
#' 'max_batch', which are used if the batch argument is not supplied.
#' @inheritParams erroneous_transitions
#' @param misdated_entries expected numbers of misdated records entered for
#' this onset day from neighbouring onset days: a named list of 'initial' and
#' 'periods' as returned by \code{misdated_entries}. Each can alternatively
#' vary across elements of cases_true: 'initial' as a vector and 'periods' as a
#' matrix with a row for each element. If NULL (the default) there are none.
#' @param batch a named list of 'phi' and 'time_detected' describing a batch of
#' erroneous records whose footprint contains the onset day (see
#' \code{sample_batches}); or NULL (the default) if there is none. Records
#' belonging to the batch are entered at rate phi * I_true * (S(day_1) - S(day_2))
#' until the batch is detected, and all are withdrawn on the first reporting
#' day on or after time_detected. Their count is summed over along with E.
#' @param max_batch maximum count of batch records for the onset day
#' @param use_compiled if TRUE (the default) the recursion is carried out by
#' compiled code; the alternative is slower and does not allow for batches
#'
#' @return a log-probability or vector of log-probabilities
observation_process_revisions_logp <- function(observation_df, cases_true,
                                               day_onset, reporting_parameters,
                                               revision_parameters,
                                               max_erroneous=3,
                                               misdated_entries=NULL,
                                               batch=NULL,
                                               max_batch=30,
                                               use_compiled=TRUE) {
  observation_days <- observation_df$time_reported
  observation_cases <- observation_df$cases_reported
  if(any(diff(observation_days) <= 0))
    stop("days in observation matrix must be increasing.")

  eta <- revision_parameters$eta
  rho <- revision_parameters$rho
  n_cases <- length(cases_true)
  if(eta < 0 | rho <= 0 | rho >= 1)
    return(rep(-Inf, n_cases))

  # batches can alternatively be carried within the revision parameters
  if(is.null(batch) & !is.null(revision_parameters$batches)) {
    batch <- batch_for_onset(day_onset, revision_parameters$batches)
    if(!is.null(revision_parameters$max_batch))
      max_batch <- revision_parameters$max_batch
  }

  if(use_compiled)
    return(observation_process_revisions_compiled_logp(
      observation_days, observation_cases, cases_true, day_onset,
      reporting_parameters, revision_parameters, max_erroneous,
      misdated_entries, batch, max_batch))
  if(!is.null(batch))
    stop("Batches are only allowed for if use_compiled is TRUE.")

  transitions <- erroneous_transitions(max_erroneous)
  n_transitions <- nrow(transitions)
  erroneous_previous <- rep(transitions$erroneous_previous, each=n_cases)
  entered <- rep(transitions$entered, each=n_cases)
  change_erroneous <- rep(transitions$erroneous_current - transitions$erroneous_previous,
                          each=n_cases)
  cases_true_long <- rep(cases_true, n_transitions)

  # log-probability of observations so far and each count of erroneous records
  num_periods <- nrow(observation_df) - 1
  mean_initial_extra <- 0
  mean_entered_extra <- matrix(0, nrow=n_cases, ncol=num_periods)
  if(!is.null(misdated_entries)) {
    mean_initial_extra <- misdated_entries$initial
    mean_entered_extra <- misdated_entries$periods
    if(!is.matrix(mean_entered_extra))
      mean_entered_extra <- matrix(mean_entered_extra, nrow=n_cases,
                                   ncol=num_periods, byrow=TRUE)
  }

  log_alpha <- erroneous_initial_logp(observation_days[1], observation_cases[1],
                                      cases_true, day_onset,
                                      reporting_parameters, revision_parameters,
                                      max_erroneous, mean_initial_extra)

  is_misdating <- misdating_prob(revision_parameters) > 0
  if(is_misdating)
    p_unrecorded <- unrecorded_prob(observation_days, day_onset,
                                    reporting_parameters, revision_parameters)

  for(i in seq_len(num_periods)) {
    day_1 <- observation_days[i]
    day_2 <- observation_days[i + 1]
    if(is_misdating)
      p_detect <- (p_unrecorded[i] - p_unrecorded[i + 1]) / p_unrecorded[i]
    else
      p_detect <- detected_after_unobserved_prob(day_2, day_1, day_onset,
                                                 reporting_parameters)
    p_newly_detected <- undetected_prob(day_1 + 0.5, day_onset, reporting_parameters) -
      undetected_prob(day_2 + 0.5, day_onset, reporting_parameters)
    p_withdrawn <- 1 - (1 - rho)^(day_2 - day_1)
    change_cases <- observation_cases[i + 1] - observation_cases[i]

    logp_withdrawn <- stats::dbinom(transitions$withdrawn,
                                    transitions$erroneous_previous,
                                    p_withdrawn, log=TRUE)
    logp_entered <- stats::dpois(
      entered,
      eta * cases_true_long * p_newly_detected +
        rep(mean_entered_extra[, i], n_transitions),
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

#' Prepares inputs for and calls the compiled forward recursion
#'
#' @param observation_days vector of days on which cases were reported
#' @param observation_cases vector of cases reported on those days
#' @inheritParams observation_process_revisions_logp
#'
#' @return a vector of log-probabilities
observation_process_revisions_compiled_logp <- function(
    observation_days, observation_cases, cases_true, day_onset,
    reporting_parameters, revision_parameters, max_erroneous,
    misdated_entries, batch, max_batch) {

  n_cases <- length(cases_true)
  num_periods <- length(observation_days) - 1
  mean_initial_extra <- 0
  mean_entered_extra <- matrix(0, nrow=n_cases, ncol=num_periods)
  if(!is.null(misdated_entries)) {
    mean_initial_extra <- misdated_entries$initial
    mean_entered_extra <- misdated_entries$periods
    if(!is.matrix(mean_entered_extra))
      mean_entered_extra <- matrix(mean_entered_extra, nrow=n_cases,
                                   ncol=num_periods, byrow=TRUE)
  }

  # a batch detected by the first reporting day has no effect
  phi <- 0
  period_detected <- 0
  if(!is.null(batch) && observation_days[1] < batch$time_detected) {
    phi <- batch$phi
    is_after <- observation_days[-1] >= batch$time_detected
    period_detected <- if(any(is_after)) which(is_after)[1] else num_periods + 1
  } else {
    max_batch <- 0
  }

  log_alpha <- erroneous_initial_logp(observation_days[1], observation_cases[1],
                                      cases_true, day_onset,
                                      reporting_parameters, revision_parameters,
                                      max_erroneous, mean_initial_extra,
                                      phi, max_batch)

  p_undetected <- 1 - pgamma_mean_sd(observation_days + 0.5 - day_onset,
                                     reporting_parameters$mean,
                                     reporting_parameters$sd)
  p_unrecorded <- p_undetected
  if(misdating_prob(revision_parameters) > 0)
    p_unrecorded <- unrecorded_prob(observation_days, day_onset,
                                    reporting_parameters, revision_parameters)

  revisions_logp_cpp(
    cases_true=as.numeric(cases_true),
    observation_cases=as.numeric(observation_cases),
    p_detect=-diff(p_unrecorded) / p_unrecorded[-(num_periods + 1)],
    p_newly_detected=-diff(p_undetected),
    p_withdrawn=1 - (1 - revision_parameters$rho)^diff(observation_days),
    eta=revision_parameters$eta,
    mean_entered_extra=mean_entered_extra,
    log_alpha_initial=log_alpha,
    max_erroneous=max_erroneous,
    phi=phi,
    period_detected=period_detected,
    max_batch=max_batch)
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
#' @param f_observation_logp a function taking a vector of true case counts and
#' returning the log-probability of all observations which depend on them. If
#' NULL (the default) this is \code{observation_process_logp} for this onset
#' day. With misdating, observations for neighbouring onset days also depend on
#' the true case count (see \code{misdating_observation_logp_function}).
#'
#' @return a vector of (unnormalised) log-probabilities
conditional_cases_revisions_logp <- function(cases_true, observation_df, cases_history,
                                             Rt, day_onset, serial_parameters, reporting_parameters,
                                             revision_parameters,
                                             kappa=NULL, is_negative_binomial=FALSE,
                                             cases_future=NULL, Rt_future=NULL,
                                             max_erroneous=3,
                                             log_negligible=50,
                                             f_observation_logp=NULL) {

  # a single observation contributes nothing to the log-probability, so this
  # gives the components due to the renewal process only
  logps_upper <- conditional_cases_logp(cases_true, observation_df[1, ], cases_history,
                                        Rt, day_onset, serial_parameters, reporting_parameters,
                                        kappa, is_negative_binomial,
                                        cases_future, Rt_future)
  if(nrow(observation_df) <= 1 & is.null(f_observation_logp))
    return(logps_upper)

  f_observation <- f_observation_logp
  if(is.null(f_observation))
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

#' Builds a function giving the log-probability of all observations which
#' depend on the true cases arising on an onset day when cases can be misdated
#'
#' With misdating, cases arising on an onset day can be recorded under
#' neighbouring onset days, so that the observations for those days also depend
#' on its true case count. The function returned sums the observation process
#' log-probability for the onset day itself and for each neighbouring onset day.
#'
#' @param day_onset the onset day whose true case count is to be varied
#' @param l_observations a list, named by onset time, of data frames with
#' columns 'time_reported' and 'cases_reported' ordered by time_reported
#' @param cases_by_onset a vector of current true case counts named by onset time
#' @param l_reporting a list, named by onset time, of reporting parameters
#' (each a named list of 'mean' and 'sd')
#' @inheritParams observation_process_revisions_logp
#'
#' @return a function taking a vector of true case counts for day_onset and
#' returning a vector of log-probabilities
misdating_observation_logp_function <- function(day_onset, l_observations,
                                                cases_by_onset, l_reporting,
                                                revision_parameters,
                                                max_erroneous=3) {
  name_onset <- as.character(day_onset)
  observation_df <- l_observations[[name_onset]]
  entries_own <- misdated_entries(observation_df$time_reported, day_onset,
                                  cases_by_onset, l_reporting[[name_onset]],
                                  revision_parameters)

  # neighbouring onset days under which cases arising on day_onset can be recorded
  l_neighbours <- list()
  for(shift in revision_parameters$shifts) {
    name_neighbour <- as.character(day_onset + shift)
    neighbour_df <- l_observations[[name_neighbour]]
    if(is.null(neighbour_df) || nrow(neighbour_df) <= 1)
      next
    l_neighbours[[name_neighbour]] <- list(
      day=day_onset + shift,
      observation_df=neighbour_df,
      cases_true=cases_by_onset[[name_neighbour]],
      reporting_parameters=l_reporting[[name_neighbour]],
      entries_other=misdated_entries(
        neighbour_df$time_reported, day_onset + shift, cases_by_onset,
        l_reporting[[name_neighbour]], revision_parameters, exclude=day_onset),
      entries_unit=misdated_entry_unit(
        neighbour_df$time_reported, day_onset + shift, shift,
        l_reporting[[name_neighbour]], revision_parameters))
  }

  function(cases) {
    logp <- rep(0, length(cases))
    if(nrow(observation_df) > 1)
      logp <- observation_process_revisions_logp(
        observation_df, cases, day_onset, l_reporting[[name_onset]],
        revision_parameters, max_erroneous, entries_own)
    for(neighbour in l_neighbours) {
      entries <- list(
        initial=neighbour$entries_other$initial +
          cases * neighbour$entries_unit$initial,
        periods=outer(cases, neighbour$entries_unit$periods) +
          matrix(neighbour$entries_other$periods, nrow=length(cases),
                 ncol=length(neighbour$entries_unit$periods), byrow=TRUE))
      logp <- logp + observation_process_revisions_logp(
        neighbour$observation_df, rep(neighbour$cases_true, length(cases)),
        neighbour$day, neighbour$reporting_parameters,
        revision_parameters, max_erroneous, entries)
    }
    logp
  }
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

  is_misdating <- !is.null(revision_parameters) &&
    misdating_prob(revision_parameters) > 0
  if(is_misdating)
    cases_by_onset <- stats::setNames(test_df$cases_true, test_df$time_onset)

  for(i in seq_along(onset_times)) {
    onset_time <- onset_times[i]
    short_df <- snapshot_with_true_cases_df %>%
      dplyr::filter(.data$time_onset==onset_time)
    cases_true <- short_df$cases_true[1]
    reporting_parameters_current <- list(mean=reporting_parameters$mean[i],
                                         sd=reporting_parameters$sd[i])
    if(nrow(short_df) > 1) {# handling cases arising today which were reported today
      if(is_misdating)
        logp_new <- observation_process_revisions_logp(
          short_df,
          cases_true=cases_true,
          day_onset=onset_time,
          reporting_parameters=reporting_parameters_current,
          revision_parameters=revision_parameters,
          max_erroneous=max_erroneous,
          misdated_entries=misdated_entries(
            short_df$time_reported, onset_time, cases_by_onset,
            reporting_parameters_current, revision_parameters))
      else
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

