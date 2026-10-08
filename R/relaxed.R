# A relaxed version of the model, written in Stan, in which nothing is
# discrete: see inst/stan/relaxed.stan. It can be fitted by Hamiltonian Monte
# Carlo or by variational approximations, using the cmdstanr package.

#' Default priors of the relaxed model
#'
#' @return a named list of: the means and sds of the gamma priors on the mean
#' and sd of the reporting delay, on eta and on the overdispersion parameter
#' kappa; the parameters of the beta prior on rho; the sd of the normal prior
#' on the log of each effect of the day of the week; the sd of the normal
#' prior on the log of Rt in the first period and the scale of the
#' half-normal prior on the sd of the steps in the log of Rt between periods
#' (if Rt follows a random walk); the shape and rate of the gamma prior on Rt
#' (if it is independent between periods); the parameters of the beta prior on
#' the weight of the tail of the delay and the mean and sd of the gamma prior
#' on its mean; the shape and scale of the inverse-gamma prior on the variance
#' of the logs of the effects of reports; and the parameters of the beta prior
#' on the probability that a case is misdated
#' @export
relaxed_priors <- function() {
  list(delay_mean_mu=7, delay_mean_sigma=10, delay_sd_mu=5, delay_sd_sigma=5,
       eta_mu=0.005, eta_sigma=0.005, rho_alpha=1, rho_beta=4,
       kappa_mu=20, kappa_sigma=20, weekday_sd=1,
       log_Rt_first_sd=1, Rt_step_sd=0.1, Rt_shape=1, Rt_rate=1,
       tail_weight_alpha=1, tail_weight_beta=3, tail_mean_mu=15, tail_mean_sigma=10,
       report_sd_shape=3, report_sd_scale=0.2,
       pi_alpha=1, pi_beta=9)
}

#' Prepares the data for the relaxed model
#'
#' @param snapshot_df a tibble with columns: time_onset, time_reported,
#' cases_reported. Onset times must be consecutive and start at 1, and no
#' report can be later than the last onset time.
#' @inheritParams state_process_logp
#' @inheritParams sample_cases_history
#' @inheritParams observation_process_revisions_logp
#' @param weekday_first the day of the week (from 1 to 7) of time 1, if there
#' is to be an effect of the day of the week of each report; NULL for none
#' @param days_sampled the number of the latest onset times whose true cases
#' are unknown. The true cases of earlier onset times are fixed at their
#' latest reported counts.
#' @param Rt_period_width number of onset times in each period within which
#' Rt is constant; periods are aligned to end on the last onset time
#' @param Rt_index the period of each onset time (a vector of whole numbers
#' from 1, non-decreasing), in place of periods of equal width
#' @param Rt_random_walk if TRUE (the default) the log of Rt follows a random
#' walk between periods; if FALSE Rt is independent between periods, with a
#' gamma prior, as when periods are given to \code{mcmc}
#' @param reporting_piece_index the period of each onset time for the
#' reporting delay, which has its own mean and sd in each (a vector of whole
#' numbers from 1); by default there is a single delay
#' @param is_negative_binomial whether the renewal model is negative binomial
#' (the default) or Poisson
#' @param tail whether a proportion of cases have a delay which is
#' exponential rather than gamma (see \code{log_undetected_days})
#' @param report_effects whether each report has its own random effect on the
#' rate of reporting
#' @param misdating NULL if cases are never misdated; or a named list of
#' 'shifts' and 'shift_probs' (see \code{unrecorded_prob}), in which case the
#' probability that a case is misdated is estimated
#' @param priors a named list as returned by \code{relaxed_priors}
#'
#' @return a named list of data for the Stan model
#' @export
#' @importFrom rlang .data
relaxed_stan_data <- function(snapshot_df, serial_parameters, cases_before=NULL,
                              max_erroneous=3, weekday_first=NULL,
                              days_sampled=35, serial_max=40, Rt_period_width=7,
                              Rt_index=NULL, Rt_random_walk=TRUE,
                              reporting_piece_index=NULL,
                              is_negative_binomial=TRUE, tail=FALSE,
                              report_effects=FALSE, misdating=NULL,
                              priors=relaxed_priors()) {
  df <- snapshot_df %>%
    dplyr::arrange(.data$time_onset, .data$time_reported)
  onset_times <- sort(unique(df$time_onset))
  n_days <- length(onset_times)
  if(!all(onset_times == seq_len(n_days)))
    stop("Onset times must be consecutive and start at 1.")
  if(max(df$time_reported) > n_days)
    stop("No report can be later than the last onset time.")

  obs_n <- as.integer(table(factor(df$time_onset, levels=onset_times)))
  obs_start <- as.integer(match(onset_times, df$time_onset))
  summary_df <- df %>%
    dplyr::group_by(.data$time_onset) %>%
    dplyr::summarise(largest=max(.data$cases_reported),
                     latest=dplyr::last(.data$cases_reported))
  is_sampled <- onset_times > n_days - days_sampled

  # the first report made on or after each day
  report_days <- sort(unique(df$time_reported))
  index_report <- findInterval(seq_len(n_days) - 1e-8, report_days) + 1
  has_report <- index_report <= length(report_days)
  report_index <- ifelse(has_report, index_report, 0L)
  weekday_report <- rep(0L, n_days)
  if(!is.null(weekday_first))
    weekday_report[has_report] <- as.integer(weekday_of_days(
      report_days[index_report[has_report]], weekday_first))

  f_check_index <- function(index, name) {
    if(length(index) != n_days || index[1] != 1 || any(!diff(index) %in% c(0, 1)))
      stop(paste0(name, " must give a period for each onset time, starting at 1 ",
                  "and rising by at most 1 at a time."))
    as.integer(index)
  }
  if(is.null(Rt_index)) {
    offset <- (Rt_period_width - n_days %% Rt_period_width) %% Rt_period_width
    period <- as.integer(ceiling((seq_len(n_days) + offset) / Rt_period_width))
  } else {
    period <- f_check_index(Rt_index, "Rt_index")
  }
  delay_period <- rep(1L, n_days)
  if(!is.null(reporting_piece_index))
    delay_period <- f_check_index(reporting_piece_index, "reporting_piece_index")

  shifts <- integer(0)
  shift_probs <- numeric(0)
  if(!is.null(misdating)) {
    shifts <- as.integer(misdating$shifts)
    shift_probs <- misdating$shift_probs
    if(length(shifts) == 0 || length(shifts) != length(shift_probs) ||
       any(shifts == 0) || abs(sum(shift_probs) - 1) > 1e-8)
      stop("Shifts must be non-zero and their probabilities must sum to one.")
  }
  unknown <- setdiff(names(priors), names(relaxed_priors()))
  if(length(unknown) > 0)
    stop(paste0("Unknown priors: ", paste(unknown, collapse=", "), "."))
  priors_all <- relaxed_priors()
  priors_all[names(priors)] <- priors

  c(list(
    n_days=n_days,
    n_obs=nrow(df),
    obs_day=as.integer(df$time_reported),
    obs_count=as.integer(df$cases_reported),
    obs_start=obs_start,
    obs_n=obs_n,
    max_erroneous=max_erroneous,
    n_sampled=sum(is_sampled),
    sampled_day=as.array(which(is_sampled)),
    cases_lower=as.numeric(summary_df$largest),
    cases_fixed=as.numeric(summary_df$latest),
    serial_max=serial_max,
    serial_weights=weights_series(serial_max, serial_parameters),
    cases_before=rev(recent_history(rev(cases_before), serial_max)),
    is_negative_binomial=as.integer(is_negative_binomial),
    n_periods=max(period),
    period=period,
    is_Rt_random_walk=as.integer(Rt_random_walk),
    n_delay_periods=max(delay_period),
    delay_period=delay_period,
    is_tail=as.integer(tail),
    is_weekday=as.integer(!is.null(weekday_first)),
    weekday_report=weekday_report,
    is_report_effects=as.integer(report_effects),
    n_reports=length(report_days),
    report_index=as.integer(report_index),
    is_misdating=as.integer(!is.null(misdating)),
    n_shifts=length(shifts),
    shifts=as.array(shifts),
    shift_probs=as.array(shift_probs)),
    priors_all)
}

#' Compiles the relaxed model
#'
#' @param ... passed to \code{cmdstanr::cmdstan_model}
#' @return a compiled model
#' @export
relaxed_model <- function(...) {
  if(!requireNamespace("cmdstanr", quietly=TRUE))
    stop("The cmdstanr package must be installed to use the relaxed model.") # nocov
  path <- system.file("stan", "relaxed.stan", package="incidenceinflation")
  cmdstanr::cmdstan_model(path, ...)
}

#' Initial values for the relaxed model
#'
#' @param data as returned by \code{relaxed_stan_data}
#' @return a function returning a named list of initial values
relaxed_inits <- function(data) {
  function() {
    days_before_end <- data$n_days - data$sampled_day
    list(
      cases_extra=as.array((1 + data$cases_lower[data$sampled_day]) *
        2 * exp(-days_before_end / 4) * stats::runif(data$n_sampled, 0.5, 1.5) + 1),
      delay_mean=as.array(stats::runif(data$n_delay_periods, 5, 9)),
      delay_sd=as.array(stats::runif(data$n_delay_periods, 3, 6)),
      log_weekday=rep(0, 7),
      eta=0.01,
      rho=0.2,
      kappa=as.array(stats::runif(data$is_negative_binomial, 5, 30)),
      log_Rt_raw=as.array(rep(0, data$n_periods)),
      Rt_step=as.array(rep(0.05, data$is_Rt_random_walk)),
      tail_weight=as.array(rep(0.05, data$is_tail)),
      tail_mean=as.array(rep(15, data$is_tail)),
      report_raw=as.array(rep(0, data$is_report_effects * data$n_reports)),
      report_sd=as.array(rep(0.2, data$is_report_effects)),
      pi_misdated=as.array(rep(0.05, data$is_misdating)))
  }
}

#' Fits the relaxed model
#'
#' @param data as returned by \code{relaxed_stan_data}
#' @param method one of: "nuts" (Hamiltonian Monte Carlo, which samples the
#' posterior of the relaxed model); "pathfinder", "advi" (mean-field) or
#' "advi_fullrank" (variational approximations); or "laplace" (a normal
#' approximation at the posterior mode)
#' @param model a compiled model as returned by \code{relaxed_model}
#' @param ... passed to the corresponding method of the model
#'
#' @return a cmdstanr fit
#' @export
fit_relaxed <- function(data, method=c("nuts", "pathfinder", "advi",
                                       "advi_fullrank", "laplace"),
                        model=relaxed_model(), ...) {
  method <- match.arg(method)
  init <- relaxed_inits(data)
  switch(
    method,
    nuts=model$sample(data=data, init=init, ...),
    pathfinder=model$pathfinder(data=data, init=init, ...),
    advi=model$variational(data=data, init=init, algorithm="meanfield", ...),
    advi_fullrank=model$variational(data=data, init=init, algorithm="fullrank", ...),
    laplace=model$laplace(data=data, init=init, ...))
}

#' Draws from a fit of the relaxed model, in the form of those of \code{mcmc}
#'
#' @param fit as returned by \code{fit_relaxed}
#'
#' @return a named list of tibbles: 'cases' (iteration, time_onset,
#' cases_true), 'Rt' (iteration, time_onset, Rt) and 'parameters' (a column
#' for each of the other parameters; where the delay differs between periods
#' its mean and sd have a column for each)
#' @export
#' @importFrom rlang .data
relaxed_draws <- function(fit) {
  f_long <- function(variable, value) {
    m <- fit$draws(variable, format="draws_matrix")
    dplyr::tibble(
      iteration=rep(seq_len(nrow(m)), ncol(m)),
      time_onset=rep(seq_len(ncol(m)), each=nrow(m)),
      value=as.numeric(m)) %>%
      dplyr::rename(!!value := "value")
  }
  # parameters which a model option leaves out are skipped
  names_parameters <- c("delay_mean", "delay_sd", "eta", "rho", "kappa", "Rt_step",
                        "tail_weight", "tail_mean", "report_sd", "pi_misdated",
                        "weekday_effect", "report_effect")
  names_fit <- unique(sub("\\[.*", "", posterior::variables(fit$draws())))
  parameters <- fit$draws(intersect(names_parameters, names_fit),
                          format="draws_df") %>%
    dplyr::as_tibble() %>%
    dplyr::select(-dplyr::any_of(c(".chain", ".iteration"))) %>%
    dplyr::rename(iteration=".draw")
  # those of which there is only one are named without an index
  for(name in c("delay_mean", "delay_sd", "kappa", "Rt_step", "tail_weight",
                "tail_mean", "report_sd", "pi_misdated")) {
    indexed <- grep(paste0("^", name, "\\["), colnames(parameters), value=TRUE)
    if(length(indexed) == 1)
      colnames(parameters)[colnames(parameters) == indexed] <- name
  }
  list(cases=f_long("cases_true", "cases_true"),
       Rt=f_long("Rt", "Rt"),
       parameters=parameters)
}
