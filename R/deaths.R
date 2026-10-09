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
#' change); and 'n_parameter_steps'. The delay here is that from onset to a
#' death being recorded.
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

#' Probabilities governing the reports of one onset time in the joint model
#' of cases and deaths
#'
#' @param days the times of the reports, in increasing order
#' @param day_onset the onset time
#' @param reporting_parameters named list of the 'mean' and 'sd' of the delay
#' from onset to a case being reported
#' @param death_parameters named list of the 'mean' and 'sd' of the delay from
#' onset to a death being recorded
#' @inheritParams log_undetected_days
#'
#' @return a named list of vectors with an element for each report: 'p' (the
#' probability a case not yet reported at the report before is reported by
#' this one), 'g' (the probability a death has been recorded by this report)
#' and 'q' (the probability a death not recorded at the report before is
#' recorded by this one). Before the first report nothing has been reported
#' or recorded.
joint_probabilities <- function(days, day_onset, reporting_parameters,
                                death_parameters, revision_parameters=NULL) {
  f <- function(parameters) {
    log_s <- log_undetected_days(days, day_onset, parameters, revision_parameters)
    list(log_s=log_s, step=-expm1(pmin(log_s - c(0, log_s[-length(log_s)]), 0)))
  }
  reported <- f(reporting_parameters)
  recorded <- f(death_parameters)
  list(p=reported$step, g=-expm1(recorded$log_s), q=recorded$step)
}

#' Binomial log-probability which is -Inf wherever the count is impossible
#'
#' @param x a count or vector of counts
#' @param size a number of trials or vector of them (recycled against x)
#' @param prob a probability
#'
#' @return a vector of log-probabilities
dbinom_log_recycled <- function(x, size, prob) {
  n <- max(length(x), length(size))
  x <- rep_len(x, n)
  size <- rep_len(size, n)
  logp <- rep(-Inf, n)
  is_valid <- size >= 0 & x >= 0 & x <= size
  logp[is_valid] <- stats::dbinom(x[is_valid], size[is_valid], prob, log=TRUE)
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
#' Deaths newly shown are those of fatal cases newly reported whose death had
#' already been recorded (each with probability g) and those of fatal cases
#' reported earlier and shown as alive whose death has now been recorded (each
#' with probability q).
#'
#' @param deaths_new the number of deaths newly shown
#' @param fatal_new the number of fatal cases newly reported (a vector)
#' @param fatal_alive the number of fatal cases reported earlier whose deaths
#' had not been recorded (a vector)
#' @param g a probability
#' @param q a probability
#'
#' @return a vector of log-probabilities
deaths_new_logp <- function(deaths_new, fatal_new, fatal_alive, g, q) {
  n <- max(length(fatal_new), length(fatal_alive))
  fatal_new <- rep_len(fatal_new, n)
  fatal_alive <- rep_len(fatal_alive, n)
  logp <- rep(-Inf, n)
  for(x in 0:min(deaths_new, max(fatal_new)))
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
#' @param probabilities as returned by \code{joint_probabilities}
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
          dbinom_log_recycled(a, deaths_true - r, p) +
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

#' Log-probability of the reports of one onset time for possible true cases
#' and true fatal cases
#'
#' @param path as returned by \code{sample_fatal_path}
#' @param p as in \code{reports_stream_logp}
#' @param cfr the case fatality ratio of the onset time
#' @param cases_true a vector of possible true cases
#'
#' @return a matrix with a row for each of cases_true and a column for each
#' number of fatal cases from that reported to the largest possible, holding
#' the log of the joint probability of the fatal cases and the reports
#' (-Inf where impossible). Its attribute 'deaths_true' holds the numbers of
#' fatal cases to which the columns correspond.
joint_cases_deaths_logp <- function(path, p, cfr, cases_true) {
  n_reports <- length(path$fatal_new)
  fatal_reported <- path$fatal_before[n_reports] + path$fatal_new[n_reports]
  survivors_reported <- path$survivors_before[n_reports] +
    path$survivors_new[n_reports]
  fatal_max <- max(max(cases_true) - survivors_reported, fatal_reported)
  fatal_values <- fatal_reported:fatal_max
  survivors_max <- max(max(cases_true) - fatal_reported, survivors_reported)
  survivors_values <- survivors_reported:survivors_max
  logp_fatal <- reports_stream_logp(path$fatal_new, path$fatal_before,
                                    fatal_values, p)
  logp_survivors <- reports_stream_logp(path$survivors_new, path$survivors_before,
                                        survivors_values, p)
  cases <- matrix(cases_true, length(cases_true), length(fatal_values))
  fatal <- matrix(fatal_values, length(cases_true), length(fatal_values), byrow=TRUE)
  survivors <- cases - fatal
  is_possible <- survivors >= survivors_reported
  logp <- matrix(-Inf, length(cases_true), length(fatal_values))
  logp[is_possible] <- stats::dbinom(fatal[is_possible], cases[is_possible], cfr,
                                     log=TRUE) +
    logp_fatal[fatal[is_possible] - fatal_reported + 1] +
    logp_survivors[survivors[is_possible] - survivors_reported + 1]
  attr(logp, "deaths_true") <- fatal_values
  logp
}

#' Starting state of the model of deaths
#'
#' @param snapshot_df a tibble with columns: time_onset, time_reported,
#' cases_reported and deaths_reported
#' @param cases_true_df a tibble with columns time_onset and cases_true
#' @param settings as returned by \code{deaths_settings}
#'
#' @return a named list holding the reports of each onset time, the true
#' fatal cases and case fatality ratio by onset time, and the segments of the
#' ratio and of the delay from onset to a death being recorded
#' @importFrom rlang .data
initial_deaths_state <- function(snapshot_df, cases_true_df, settings) {
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
  list(
    reports=reports,
    onset_times=onset_times,
    deaths_true=deaths_true,
    cfr=rep(settings$initial_cfr, n_days),
    cfr_starts=1,
    p_change=settings$alpha / (settings$alpha + settings$beta),
    recording=initial_reporting_changepoints(settings$initial_delay),
    paths=NULL)
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

#' Samples how many of the cases in each report are fatal, for every onset time
#'
#' @inheritParams deaths_delay_by_onset
#' @param cases_true a vector of the true cases of each onset time
#' @param reporting_by_onset a tibble with columns 'mean' and 'sd' and a row
#' for each onset time: the delay from onset to a case being reported
#' @inheritParams log_undetected_days
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
#' \code{fatal_path_forward}. The split is then drawn given the result.
#'
#' @return the state with elements 'deaths_true' and 'paths' updated (see
#' \code{sample_fatal_path}; each path also holds the probability 'p' used)
deaths_paths_step <- function(state, cases_true, reporting_by_onset,
                              revision_parameters=NULL, n_proposals=2,
                              step_max=3) {
  delays <- deaths_delay_by_onset(state)
  paths <- vector("list", length(state$onset_times))
  for(k in seq_along(state$onset_times)) {
    reports <- state$reports[[k]]
    probabilities <- joint_probabilities(
      reports$time_reported, state$onset_times[k],
      list(mean=reporting_by_onset$mean[k], sd=reporting_by_onset$sd[k]),
      list(mean=delays$mean[k], sd=delays$sd[k]), revision_parameters)
    f_forward <- function(deaths_true)
      fatal_path_forward(reports$cases_reported, reports$deaths_reported,
                         deaths_true, cases_true[k] - deaths_true, probabilities)
    f_logp <- function(deaths_true, forward)
      forward$logp + stats::dbinom(deaths_true, cases_true[k], state$cfr[k], log=TRUE)
    deaths_true <- state$deaths_true[k]
    forward <- f_forward(deaths_true)
    if(!is.finite(forward$logp))
      stop("The reports of cases and deaths are impossible given the true cases and deaths.")
    largest <- max(reports$deaths_reported)
    for(proposal in seq_len(n_proposals)) {
      proposed <- deaths_true + sample(c(-1, 1), 1) * sample.int(step_max, 1)
      if(proposed < largest | proposed > cases_true[k])
        next
      forward_proposed <- f_forward(proposed)
      log_ratio <- f_logp(proposed, forward_proposed) - f_logp(deaths_true, forward)
      if(is.finite(log_ratio) && log(stats::runif(1)) < log_ratio) {
        deaths_true <- proposed
        forward <- forward_proposed
      }
    }
    state$deaths_true[k] <- deaths_true
    path <- fatal_path_backward(forward)
    path$p <- probabilities$p
    paths[[k]] <- path
  }
  state$paths <- paths
  state
}

#' Creates the function giving the log-probability of the reports of an onset
#' time for possible true cases, with the true fatal cases summed out
#'
#' @inheritParams deaths_delay_by_onset
#' @param time_onset an onset time
#'
#' @return a function of a vector of possible true cases
deaths_cases_logp_function <- function(state, time_onset) {
  k <- match(time_onset, state$onset_times)
  path <- state$paths[[k]]
  cfr <- state$cfr[k]
  function(cases_true) {
    logp <- joint_cases_deaths_logp(path, path$p, cfr, cases_true)
    matrixStats::rowLogSumExps(logp)
  }
}

#' Samples the true fatal cases of every onset time given the true cases
#'
#' @inheritParams deaths_paths_step
#'
#' @return the state with element 'deaths_true' updated
deaths_true_step <- function(state, cases_true) {
  for(k in seq_along(state$onset_times)) {
    path <- state$paths[[k]]
    logp <- joint_cases_deaths_logp(path, path$p, state$cfr[k], cases_true[k])
    values <- attr(logp, "deaths_true")
    logp <- as.vector(logp)
    if(!any(is.finite(logp)))
      stop("The true cases of an onset time are fewer than the cases reported.")
    state$deaths_true[k] <- values[sample.int(length(values), 1,
                                              prob=exp(logp - max(logp)))]
  }
  state
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
  survivors_true <- cases_true - state$deaths_true
  function(indices, reporting_parameters) {
    if(reporting_parameters$mean <= 0 | reporting_parameters$sd <= 0)
      return(-Inf)
    sum(purrr::map_dbl(indices, function(k) {
      path <- state$paths[[k]]
      log_s <- log_undetected_days(state$reports[[k]]$time_reported,
                                   state$onset_times[k], reporting_parameters,
                                   revision_parameters)
      p <- -expm1(pmin(log_s - c(0, log_s[-length(log_s)]), 0))
      reports_stream_logp(path$fatal_new, path$fatal_before,
                          state$deaths_true[k], p) +
        reports_stream_logp(path$survivors_new, path$survivors_before,
                            survivors_true[k], p)
    }))
  }
}

#' Creates a function giving the log-probability of the deaths newly shown
#' for sets of onset times which share a delay from onset to a death being
#' recorded
#'
#' @inheritParams deaths_paths_step
#'
#' @return a function of 'indices' and 'reporting_parameters' (here the
#' parameters of the delay from onset to a death being recorded)
deaths_recording_logp_function <- function(state, revision_parameters=NULL) {
  function(indices, reporting_parameters) {
    if(reporting_parameters$mean <= 0 | reporting_parameters$sd <= 0)
      return(-Inf)
    sum(purrr::map_dbl(indices, function(k) {
      path <- state$paths[[k]]
      log_s <- log_undetected_days(state$reports[[k]]$time_reported,
                                   state$onset_times[k], reporting_parameters,
                                   revision_parameters)
      g <- -expm1(log_s)
      q <- -expm1(pmin(log_s - c(0, log_s[-length(log_s)]), 0))
      is_informative <- path$fatal_new > 0 | path$alive_before > 0
      sum(purrr::map_dbl(which(is_informative), function(j) {
        deaths_new_logp(path$deaths_new[j], path$fatal_new[j],
                        path$alive_before[j], g[j], q[j])
      }))
    }))
  }
}

#' Updates the delay from onset to a death being recorded and the case
#' fatality ratio, with the onset times at which each changes
#'
#' @inheritParams deaths_paths_step
#' @inheritParams initial_deaths_state
#'
#' @return the updated state
deaths_parameters_step <- function(state, cases_true, settings,
                                   revision_parameters=NULL) {
  n_days <- length(state$onset_times)
  state$recording <- sample_reporting_changepoints(
    state$recording,
    deaths_recording_logp_function(state, revision_parameters),
    n_days, settings$delay_prior, settings$delay_changepoints)

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
#' with the true fatal cases, the case fatality ratio and the delay from
#' onset to a death being recorded) and 'deaths_parameters' (a single row)
deaths_draws <- function(state, iteration) {
  delays <- deaths_delay_by_onset(state)
  list(
    deaths=dplyr::tibble(
      iteration=iteration,
      time_onset=state$onset_times,
      deaths_true=state$deaths_true,
      cfr=state$cfr,
      delay_mean=delays$mean,
      delay_sd=delays$sd),
    deaths_parameters=dplyr::tibble(
      iteration=iteration,
      n_cfr_changes=length(state$cfr_starts) - 1,
      p_change=state$p_change,
      n_delay_changes=length(state$recording$starts) - 1,
      delay_centre_mean=state$recording$centre$mean,
      delay_centre_sd=state$recording$centre$sd,
      delay_spread=state$recording$spread))
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
#'
#' @return a named list of 'snapshots' (a tibble with columns time_onset,
#' time_reported, cases_reported and deaths_reported) and 'truth' (a tibble
#' with columns time_onset, cases_true, deaths_true, cfr, and the mean and sd
#' of the delay from onset to being shown as dead among the cases simulated)
#' @export
#' @importFrom rlang .data
simulate_cases_deaths <- function(cases_true, report_days, reporting_parameters,
                                  cfr, death_parameters,
                                  fatal_reporting_parameters=NULL,
                                  death_from=c("onset", "report")) {
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
  f_delay <- function(n, mean, sd)
    stats::rgamma(n, shape=mean^2 / sd^2, rate=mean / sd^2)

  results <- purrr::map(seq_len(n_days), function(t) {
    n <- cases_true[t]
    is_fatal <- stats::runif(n) < cfr[t]
    delay_reported <- ifelse(
      is_fatal,
      f_delay(n, fatal_reporting$mean[t], fatal_reporting$sd[t]),
      f_delay(n, reporting$mean[t], reporting$sd[t]))
    delay_recorded <- f_delay(n, death$mean[t], death$sd[t])
    if(death_from == "report")
      delay_recorded <- delay_reported + delay_recorded
    delay_dead <- pmax(delay_reported, delay_recorded)
    days <- report_days[report_days >= t]
    elapsed <- days - t + 0.5
    list(
      snapshots=dplyr::tibble(
        time_onset=t,
        time_reported=days,
        cases_reported=purrr::map_dbl(elapsed, ~sum(delay_reported <= .)),
        deaths_reported=purrr::map_dbl(elapsed, ~sum(delay_dead[is_fatal] <= .))),
      truth=dplyr::tibble(
        time_onset=t, cases_true=n, deaths_true=sum(is_fatal), cfr=cfr[t],
        delay_dead_mean=mean(delay_dead[is_fatal]),
        delay_dead_sd=stats::sd(delay_dead[is_fatal])))
  })
  list(snapshots=purrr::map_dfr(results, "snapshots"),
       truth=purrr::map_dfr(results, "truth"))
}
