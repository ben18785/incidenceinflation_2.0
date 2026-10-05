#' Generate weights series
#'
#' @param t_max maximum time to generate weights until
#' @inheritParams gamma_discrete_pmf
#'
#' @return a vector of normalised weights of length t_max
weights_series <- function(t_max, serial_parameters) {
  day_series <- seq(1, t_max, 1)
  w <- purrr::map_dbl(day_series, ~gamma_discrete_pmf(., serial_parameters))
  # due to truncation, normalise series
  w <- w / sum(w)
  w
}

#' Gives cases expected given history of cases and Rt
#'
#' @inheritParams state_process_logp
#' @param weights a vector of normalised weights
#'
#' @return an expected number of cases
expected_cases <- function(Rt, weights, cases_history) {
  if(length(weights) != length(cases_history))
    stop("weights and history of cases must be same length.")
  Rt * sum(weights * cases_history)
}

#' Reparameterisation of gamma cdf using mean and sd
#'
#' @param x a value to be evaluated at
#' @param mu a mean
#' @param sigma a standard deviation
#'
#' @return a cumulative density
pgamma_mean_sd <- function(x, mu, sigma) {
  shape <- mu^2 / sigma^2
  rate <- mu / sigma^2
  stats::pgamma(x, shape, rate)
}

#' Reparameterisation of gamma pdf using mean and sd
#'
#' @param x a value to be evaluated at
#' @param mu a mean
#' @param sigma a standard deviation
#' @param ... other arguments passed to method
#' @return a density
dgamma_mean_sd <- function(x, mu, sigma, ...) {
  shape <- mu^2 / sigma^2
  rate <- mu / sigma^2
  stats::dgamma(x, shape, rate, ...)
}

#' Reparameterisation of gamma inverse-cdf using mean and sd
#'
#' @param x a value to be evaluated at
#' @param mu a mean
#' @param sigma a standard deviation
#'
#' @return a positive value
qgamma_mean_sd <- function(x, mu, sigma) {
  shape <- mu^2 / sigma^2
  rate <- mu / sigma^2
  stats::qgamma(x, shape, rate)
}

#' Discrete gamma probability mass function
#'
#' @param day day to evaluate pmf
#' @param serial_parameters named list of 'mean' and 'sd' of gamma distribution
#' characterising the serial interval distribution
#'
#' @return a probability
gamma_discrete_pmf <- function(day, serial_parameters){
  delay_mean <- serial_parameters$mean
  delay_sd <- serial_parameters$sd
  pgamma_mean_sd(day + 0.5, delay_mean, delay_sd) -
    pgamma_mean_sd(day - 0.5, delay_mean, delay_sd)
}

check_parameter_names <- function(reporting_parameters, serial_parameters) {
  correct_names <- c("mean", "sd")
  if(!all(correct_names %in% names(reporting_parameters)))
    stop("reporting_parameters must contain 'mean' and 'sd'.")
  if(!all(correct_names %in% names(serial_parameters)))
    stop("serial_parameters must contain 'mean' and 'sd'.")
}

#' Calculates the largest downward revision in reported cases
#'
#' For each onset time, this is the largest fall in reported cases from any
#' earlier count. When downward revisions are allowed, the maximum number of
#' erroneous records (max_erroneous) must be at least as large as this, and
#' should exceed it.
#'
#' @param snapshot_df a tibble with columns: 'time_onset', 'time_reported',
#' 'cases_reported'
#'
#' @return the largest downward revision across all onset times
#' @export
#' @importFrom rlang .data
max_downward_revision <- function(snapshot_df) {
  df <- snapshot_df %>%
    dplyr::group_by(.data$time_onset) %>%
    dplyr::arrange(.data$time_reported, .by_group=TRUE) %>%
    dplyr::summarise(
      revision=max(cummax(.data$cases_reported) - .data$cases_reported))
  max(df$revision)
}

#' Checks reported cases are possible when downward revisions are allowed
#'
#' Reported cases for an onset time cannot fall by more than max_erroneous from
#' any earlier count. Also, if cases are first reported on the onset day itself
#' (when all records are assumed genuine) they cannot subsequently fall below
#' the count first reported.
#'
#' @inheritParams max_downward_revision
#' @inheritParams erroneous_transitions
#'
#' @return TRUE (invisibly) if reported cases are possible, otherwise throws an error
#' @importFrom rlang .data
check_revisions_possible <- function(snapshot_df, max_erroneous) {
  if(max_downward_revision(snapshot_df) > max_erroneous)
    stop("Reported cases are revised downwards by more than max_erroneous: increase max_erroneous.")

  df <- snapshot_df %>%
    dplyr::group_by(.data$time_onset) %>%
    dplyr::arrange(.data$time_reported, .by_group=TRUE) %>%
    dplyr::summarise(
      is_below_first=.data$time_reported[1] - .data$time_onset[1] < 1 &
        any(.data$cases_reported < .data$cases_reported[1]))
  if(any(df$is_below_first))
    stop("Reported cases fall below the count first reported on the onset day, which is not allowed.")
  invisible(TRUE)
}
