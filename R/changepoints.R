#' Checks that a set of segments is a valid partition of a run of days
#'
#' Segments are contiguous runs of days, described by the day on which each
#' starts. The first segment starts on day 1. Every segment must contain at
#' least min_length days, unless there is only a single segment.
#'
#' @param starts a vector of the days on which segments start
#' @param n_days the number of days partitioned
#' @param min_length the minimum number of days in a segment
#'
#' @return TRUE if the segments are a valid partition and FALSE otherwise
is_valid_partition <- function(starts, n_days, min_length=1) {
  if(length(starts) == 0 || starts[1] != 1)
    return(FALSE)
  if(any(diff(starts) <= 0) || max(starts) > n_days)
    return(FALSE)
  if(length(starts) == 1)
    return(TRUE)
  all(diff(c(starts, n_days + 1)) >= min_length)
}

#' The segment to which each day belongs
#'
#' @inheritParams is_valid_partition
#'
#' @return a vector of length n_days of segment indices
segment_indices <- function(starts, n_days) {
  findInterval(seq_len(n_days), starts)
}

#' Log-prior probability (up to a constant) of a partition with a given number
#' of changes
#'
#' Each day after the first starts a new segment with probability p, and p is
#' given a beta(alpha, beta) prior. With p integrated out, the prior
#' probability of a partition depends only on its number of changes K:
#' \deqn{p(partition) \propto B(alpha + K, beta + n - K) / B(alpha, beta)}
#' where n is one fewer than the number of days. Partitions with a segment
#' shorter than the minimum length have zero probability, which is not checked
#' here (see \code{is_valid_partition}).
#'
#' @param n_changes number of changes (one fewer than the number of segments)
#' @inheritParams is_valid_partition
#' @param alpha first parameter of the beta prior on the probability of a
#' change on each day
#' @param beta second parameter of the beta prior on the probability of a
#' change on each day
#'
#' @return a log-probability
changepoints_logp_prior <- function(n_changes, n_days, alpha, beta) {
  n <- n_days - 1
  lbeta(alpha + n_changes, beta + n - n_changes) - lbeta(alpha, beta)
}

#' Number of partitions of a run of days with a given number of changes in
#' which every segment has at least a minimum length
#'
#' @inheritParams changepoints_logp_prior
#' @inheritParams is_valid_partition
#'
#' @return a count
n_valid_partitions <- function(n_changes, n_days, min_length=1) {
  if(n_changes == 0)
    return(1)
  n_spare <- n_days - (n_changes + 1) * min_length
  if(n_spare < 0)
    return(0)
  choose(n_spare + n_changes, n_changes)
}

#' Samples the probability of a change on each day given the number of changes
#'
#' @inheritParams changepoints_logp_prior
#'
#' @return a probability
sample_change_probability <- function(n_changes, n_days, alpha, beta) {
  stats::rbeta(1, alpha + n_changes, beta + n_days - 1 - n_changes)
}

#' Calculates, for each day, the log of the summed weight of all partitions of
#' the days up to and including it
#'
#' The weight of a partition is the product over its segments of the marginal
#' likelihood of the segment, and over days after the first of p (if the day
#' starts a segment) or 1 - p (if not).
#'
#' @param log_marginals a matrix whose element (s, t) is the log marginal
#' likelihood of a segment which starts on day s and ends on day t (for s <= t)
#' @param p_change probability that each day after the first starts a segment
#' @inheritParams is_valid_partition
#'
#' @return a vector of length one more than the number of days, whose first
#' element (zero) corresponds to no days
partition_forward <- function(log_marginals, p_change, min_length=1) {
  n_days <- nrow(log_marginals)
  log_forward <- c(0, rep(-Inf, n_days))
  for(t in seq_len(n_days))
    log_forward[t + 1] <- matrixStats::logSumExp(
      partition_last_segment_logp(t, log_forward, log_marginals, p_change,
                                  min_length))
  log_forward
}

#' Log-weights of the possible starts of a segment which ends on a given day
#'
#' @param day_end the day on which the segment ends
#' @param log_forward the result of \code{partition_forward} (of which only
#' the elements for days before day_end are used)
#' @inheritParams partition_forward
#'
#' @return a vector of length day_end of log-weights of the segment starting
#' on each of days 1 to day_end
partition_last_segment_logp <- function(day_end, log_forward, log_marginals,
                                        p_change, min_length=1) {
  n_days <- nrow(log_marginals)
  starts <- seq_len(day_end)
  log_w <- log_forward[starts] + log_marginals[starts, day_end] +
    ifelse(starts > 1, log(p_change), 0) + (day_end - starts) * log1p(-p_change)
  is_long_enough <- (day_end - starts + 1) >= min_length |
    (starts == 1 & day_end == n_days)
  log_w[!is_long_enough] <- -Inf
  log_w
}

#' Samples a partition of a run of days into segments
#'
#' Draws exactly from the distribution over partitions whose weights are
#' described in \code{partition_forward}, by sampling the start of the last
#' segment, then the start of the segment before it, and so on.
#'
#' @inheritParams partition_forward
#'
#' @return a vector of the days on which segments start
sample_partition <- function(log_marginals, p_change, min_length=1) {
  n_days <- nrow(log_marginals)
  log_forward <- partition_forward(log_marginals, p_change, min_length)
  starts <- c()
  day_end <- n_days
  while(day_end >= 1) {
    log_w <- partition_last_segment_logp(day_end, log_forward, log_marginals,
                                         p_change, min_length)
    probs <- exp(log_w - matrixStats::logSumExp(log_w))
    start <- sample.int(day_end, size=1, prob=probs)
    starts <- c(start, starts)
    day_end <- start - 1
  }
  starts
}

#' Samples a partition of a run of days from its prior given the probability
#' of a change on each day
#'
#' @inheritParams is_valid_partition
#' @inheritParams partition_forward
#'
#' @return a vector of the days on which segments start
#' @export
sample_changepoints <- function(n_days, p_change, min_length=1) {
  sample_partition(matrix(0, n_days, n_days), p_change, min_length)
}

#' Expands values held by segments into a value for each day
#'
#' @inheritParams is_valid_partition
#' @param values a vector with one value for each segment
#'
#' @return a vector of length n_days
#' @export
piecewise_values <- function(starts, values, n_days) {
  if(length(starts) != length(values))
    stop("There must be one value for each segment.")
  if(!is_valid_partition(starts, n_days))
    stop("Segments must start on day 1 and on increasing days no later than n_days.")
  values[segment_indices(starts, n_days)]
}

#' The quantities through which each onset day enters the renewal likelihood
#'
#' These are the true case count on the day and the total infectiousness of
#' earlier cases, \deqn{\Lambda_t = \sum_s w_s cases_true_t-s}. The renewal
#' model says nothing about days with no infectiousness before them (such as
#' the first), which are marked as uninformative.
#'
#' @param cases_true a vector of true case counts on consecutive onset days,
#' from the earliest to the latest
#' @inheritParams sample_Rt_single_piece
#'
#' @return a tibble with a row for each day and columns: 'cases_true',
#' 'infectiousness' and 'is_informative'
renewal_statistics <- function(cases_true, serial_parameters, serial_max=40) {
  w <- weights_series(serial_max, serial_parameters)
  n_days <- length(cases_true)
  cases_padded <- c(rep(0, serial_max), cases_true)
  infectiousness <- purrr::map_dbl(seq_len(n_days), function(t) {
    sum(w * rev(cases_padded[t:(t + serial_max - 1)]))
  })
  dplyr::tibble(
    cases_true=cases_true,
    infectiousness=infectiousness,
    is_informative=infectiousness > 0)
}

#' Log marginal likelihood of every possible segment of constant Rt under a
#' Poisson renewal model
#'
#' With a gamma(shape, rate) prior on the Rt of a segment, its marginal
#' likelihood has a closed form. Terms which are the same for every partition
#' are left out:
#' \deqn{shape log(rate) - log Gamma(shape) + log Gamma(shape + Y) - (shape + Y) log(rate + L)}
#' where Y and L are the sums of the case counts and of the infectiousness
#' over the informative days of the segment.
#'
#' @param statistics the result of \code{renewal_statistics}
#' @param Rt_prior_parameters a list with elements 'shape' and 'rate'
#' describing the gamma prior for Rt
#'
#' @return a matrix whose element (s, t) is the log marginal likelihood of a
#' segment which starts on day s and ends on day t (-Inf if s > t)
segment_log_marginals_poisson <- function(statistics, Rt_prior_parameters) {
  shape <- Rt_prior_parameters$shape
  rate <- Rt_prior_parameters$rate
  sums <- segment_sums(statistics)
  # elements below the diagonal are not segments
  sums$cases[lower.tri(sums$cases)] <- 0
  sums$infectiousness[lower.tri(sums$infectiousness)] <- 0
  log_marginals <- shape * log(rate) - lgamma(shape) +
    lgamma(shape + sums$cases) - (shape + sums$cases) * log(rate + sums$infectiousness)
  log_marginals[lower.tri(log_marginals)] <- -Inf
  log_marginals
}

#' Sums of case counts and infectiousness over the informative days of every
#' possible segment
#'
#' @inheritParams segment_log_marginals_poisson
#'
#' @return a named list of two matrices, 'cases' and 'infectiousness', whose
#' elements (s, t) are sums over days s to t
segment_sums <- function(statistics) {
  y <- ifelse(statistics$is_informative, statistics$cases_true, 0)
  lambda <- ifelse(statistics$is_informative, statistics$infectiousness, 0)
  c_y <- c(0, cumsum(y))
  c_lambda <- c(0, cumsum(lambda))
  n_days <- length(y)
  list(
    cases=outer(seq_len(n_days), seq_len(n_days),
                function(s, t) c_y[t + 1] - c_y[s]),
    infectiousness=outer(seq_len(n_days), seq_len(n_days),
                         function(s, t) c_lambda[t + 1] - c_lambda[s]))
}

#' A grid of Rt values used to integrate over Rt numerically
#'
#' Values are equally spaced on the log scale. Each comes with the width of
#' the interval of Rt which it represents.
#'
#' @param n number of values
#' @param lower smallest value
#' @param upper largest value
#'
#' @return a tibble with columns 'Rt', 'lower' and 'upper', the last two being
#' the ends of the interval which each value represents
Rt_grid <- function(n=400, lower=0.02, upper=20) {
  log_step <- (log(upper) - log(lower)) / (n - 1)
  log_Rt <- seq(log(lower), log(upper), length.out=n)
  dplyr::tibble(
    Rt=exp(log_Rt),
    lower=exp(log_Rt - log_step / 2),
    upper=exp(log_Rt + log_step / 2))
}

#' Log-likelihood of each onset day at each of a grid of Rt values
#'
#' @inheritParams segment_log_marginals_poisson
#' @param grid the result of \code{Rt_grid}
#' @param kappa overdispersion parameter of a negative binomial renewal model
#' @param is_negative_binomial whether the renewal model is negative binomial
#' (if FALSE it is Poisson)
#'
#' @return a matrix with a row for each day and a column for each Rt value;
#' rows of uninformative days are zero
renewal_logl_grid <- function(statistics, grid, kappa=NULL,
                              is_negative_binomial=TRUE) {
  mean_cases <- outer(statistics$infectiousness, grid$Rt)
  # placeholder for uninformative days, whose rows are set to zero below
  mean_cases[!statistics$is_informative, ] <- 1
  cases <- matrix(statistics$cases_true, nrow(statistics), nrow(grid))
  if(is_negative_binomial)
    logl <- stats::dnbinom(cases, mu=mean_cases, size=kappa, log=TRUE)
  else
    logl <- stats::dpois(cases, lambda=mean_cases, log=TRUE)
  logl <- matrix(logl, nrow(statistics), nrow(grid))
  logl[!statistics$is_informative, ] <- 0
  logl
}

#' Log of the prior probability of the interval of Rt represented by each
#' value of a grid
#'
#' @inheritParams renewal_logl_grid
#' @inheritParams segment_log_marginals_poisson
#'
#' @return a vector of log-probabilities
Rt_grid_log_prior <- function(grid, Rt_prior_parameters) {
  log(stats::pgamma(grid$upper, Rt_prior_parameters$shape, Rt_prior_parameters$rate) -
        stats::pgamma(grid$lower, Rt_prior_parameters$shape, Rt_prior_parameters$rate))
}

#' Log marginal likelihood of every possible segment, integrating over the
#' value held by the segment numerically
#'
#' @param logl a matrix of the log-likelihood of each day (rows) at each of a
#' grid of values (columns)
#' @param log_prior a vector of the log of the prior probability of each value
#' of the grid
#'
#' @return a matrix whose element (s, t) is the log marginal likelihood of a
#' segment which starts on day s and ends on day t (-Inf if s > t)
segment_log_marginals_grid <- function(logl, log_prior) {
  n_days <- nrow(logl)
  n_grid <- ncol(logl)
  c_logl <- rbind(rep(0, n_grid), apply(logl, 2, cumsum))
  if(n_days == 1)
    c_logl <- rbind(rep(0, n_grid), logl)
  log_marginals <- matrix(-Inf, n_days, n_days)
  for(t in seq_len(n_days)) {
    # row s holds the summed log-likelihood of days s to t at each value
    logl_segments <- matrix(c_logl[t + 1, ], t, n_grid, byrow=TRUE) -
      c_logl[seq_len(t), , drop=FALSE]
    log_marginals[seq_len(t), t] <- matrixStats::rowLogSumExps(
      sweep(logl_segments, 2, log_prior, "+"))
  }
  log_marginals
}

#' Samples Rt for a segment from a grid of values
#'
#' A value of the grid is drawn with probability proportional to its prior
#' probability and the likelihood of the days of the segment. Rt is then drawn
#' uniformly from the interval which that value represents.
#'
#' @param logl_segment a vector of the summed log-likelihood of the days of
#' the segment at each value of the grid
#' @inheritParams segment_log_marginals_grid
#' @inheritParams renewal_logl_grid
#'
#' @return an Rt value
sample_Rt_grid <- function(logl_segment, log_prior, grid) {
  log_w <- logl_segment + log_prior
  probs <- exp(log_w - matrixStats::logSumExp(log_w))
  k <- sample.int(nrow(grid), size=1, prob=probs)
  stats::runif(1, grid$lower[k], grid$upper[k])
}

#' Samples a piecewise-constant Rt whose changes are unknown
#'
#' Onset days are partitioned into segments, within each of which Rt is
#' constant. The partition and the Rt of each segment are drawn jointly from
#' their distribution given the true case counts and the probability of a
#' change on each day: the Rt values are integrated out to give the marginal
#' likelihood of every possible segment; a partition is drawn exactly (see
#' \code{sample_partition}); and Rt is then drawn for each of its segments.
#'
#' Under a Poisson renewal model the integral over Rt has a closed form. Under
#' a negative binomial renewal model it is calculated numerically over a grid
#' of Rt values (see \code{Rt_grid}).
#'
#' @param cases_history_df a tibble with columns 'time_onset' and
#' 'cases_true', with a row for each of a run of consecutive onset times
#' @inheritParams segment_log_marginals_poisson
#' @inheritParams sample_Rt_single_piece
#' @inheritParams partition_forward
#' @inheritParams renewal_logl_grid
#' @param use_grid whether to integrate over Rt numerically under a Poisson
#' renewal model (always the case under a negative binomial model)
#'
#' @return a tibble with a row for each onset time and columns: 'time_onset',
#' 'Rt_index' (the segment to which the onset time belongs) and 'Rt'
#' @export
sample_Rt_changepoints <- function(cases_history_df,
                                   Rt_prior_parameters,
                                   serial_parameters,
                                   p_change,
                                   min_length=1,
                                   serial_max=40,
                                   kappa=NULL,
                                   is_negative_binomial=FALSE,
                                   grid=Rt_grid(),
                                   use_grid=FALSE) {

  cases_history_df <- cases_history_df %>%
    dplyr::arrange(.data$time_onset)
  if(any(diff(cases_history_df$time_onset) != 1))
    stop("Onset times must be consecutive.")
  if(is_negative_binomial) {
    if(is.null(kappa))
      stop("Overdispersion parameter must not be null if using a negative binomial model.")
    else if(kappa <= 0)
      stop("Overdispersion parameter must be positive.")
  }
  if(p_change <= 0 | p_change >= 1)
    stop("The probability of a change must lie strictly between 0 and 1.")

  n_days <- nrow(cases_history_df)
  statistics <- renewal_statistics(cases_history_df$cases_true,
                                   serial_parameters, serial_max)
  is_grid <- is_negative_binomial | use_grid
  if(is_grid) {
    logl <- renewal_logl_grid(statistics, grid, kappa, is_negative_binomial)
    log_prior <- Rt_grid_log_prior(grid, Rt_prior_parameters)
    log_marginals <- segment_log_marginals_grid(logl, log_prior)
  } else {
    log_marginals <- segment_log_marginals_poisson(statistics, Rt_prior_parameters)
  }

  starts <- sample_partition(log_marginals, p_change, min_length)
  ends <- c(starts[-1] - 1, n_days)
  Rt <- purrr::map_dbl(seq_along(starts), function(j) {
    days <- starts[j]:ends[j]
    if(is_grid) {
      sample_Rt_grid(colSums(logl[days, , drop=FALSE]), log_prior, grid)
    } else {
      is_used <- statistics$is_informative[days]
      stats::rgamma(
        1,
        Rt_prior_parameters$shape + sum(statistics$cases_true[days][is_used]),
        Rt_prior_parameters$rate + sum(statistics$infectiousness[days][is_used]))
    }
  })

  indices <- segment_indices(starts, n_days)
  dplyr::tibble(
    time_onset=cases_history_df$time_onset,
    Rt_index=indices,
    Rt=Rt[indices])
}
