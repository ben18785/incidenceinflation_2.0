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
#' of the logs of the effects of reports; the parameters of the beta prior
#' on the probability that a case is misdated; and the parameters of the beta
#' prior on the probability that each onset time starts a new segment, where
#' changes are at unknown times (the second defaults to the number of onset
#' times, as in \code{changepoints_logp_prior}); and, where deaths are
#' modelled, the parameters of the beta prior on the case fatality ratio, the
#' means and sds of the gamma priors on the mean and sd of the delay from
#' onset to a death being recorded, and the parameters of the beta prior on
#' the probability that a wrong record is shown as dead
#' @export
relaxed_priors <- function() {
  list(delay_mean_mu=7, delay_mean_sigma=10, delay_sd_mu=5, delay_sd_sigma=5,
       eta_mu=0.005, eta_sigma=0.005, rho_alpha=1, rho_beta=4,
       kappa_mu=20, kappa_sigma=20, weekday_sd=1,
       log_Rt_first_sd=1, Rt_step_sd=0.1, Rt_shape=1, Rt_rate=1,
       tail_weight_alpha=1, tail_weight_beta=3, tail_mean_mu=15, tail_mean_sigma=10,
       report_sd_shape=3, report_sd_scale=0.2,
       pi_alpha=1, pi_beta=9,
       change_alpha=1, change_beta=NA,
       cfr_alpha=1, cfr_beta=1,
       recording_mean_mu=14, recording_mean_sigma=10,
       recording_sd_mu=8, recording_sd_sigma=8,
       wrong_dead_alpha=1, wrong_dead_beta=1)
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
#' @param Rt_max_changes if greater than zero, Rt changes at unknown onset
#' times, at most this number of times, and is constant between changes with
#' the gamma prior of \code{relaxed_priors}. Rt_index, Rt_period_width and
#' Rt_random_walk are then ignored. The prior on the changes is that of
#' \code{changepoints_logp_prior}, limited to this number of changes and
#' without a minimum length of segment. The default, zero, switches this off.
#' @param delay_max_changes as Rt_max_changes, for the mean and sd of the
#' reporting delay; reporting_piece_index is then ignored. Each change allowed
#' adds about as much computation as the whole model of reporting without
#' changes. This cannot be combined with misdating.
#' @param deaths_df NULL, or a tibble with columns time_onset, time_reported
#' and deaths_reported (the number of the cases in each report which are
#' shown as dead), with a row for each row of snapshot_df. The deaths are then
#' modelled too: a case dies with a probability (the case fatality ratio) and
#' its death is recorded after a delay from onset which is independent of the
#' reporting delay. Fatal cases are taken to be reported as others are, and
#' all to be recorded for their day of onset.
#' @param cfr_max_changes the largest number of changes in the case fatality
#' ratio, which are at unknown onset times (as Rt_max_changes)
#' @param recording_max_changes the same for the mean and sd of the delay
#' from onset to a death being recorded. The changes in the two are summed
#' over together, so the computation for deaths grows with the product of one
#' more than each.
#' @param recording_weekday whether the delay to a death being recorded has
#' its own effects of the day of the week of each report (by default, if the
#' reporting delay has them)
#' @param fatal_delay FALSE (the default) if fatal cases are reported as
#' others are; TRUE, or a named list of 'probability_same' and 'sd', if their
#' reporting delay can differ: it is then that of other cases with its mean
#' and sd multiplied by a ratio whose log is zero with prior probability
#' probability_same (by default 0.5) and otherwise normal with mean zero and
#' standard deviation sd (by default 0.5)
#' @param dated_by_death NULL, or what is supplied about fatal cases which
#' are recorded for their day of death and not their day of onset (see
#' \code{dated_by_death_settings}). The true cases estimated are then those
#' recorded for each day, from which the cases by day of onset follow.
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
                              Rt_max_changes=0, delay_max_changes=0,
                              deaths_df=NULL, cfr_max_changes=0,
                              recording_max_changes=0,
                              recording_weekday=!is.null(weekday_first),
                              fatal_delay=FALSE, dated_by_death=NULL,
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
  # changes at unknown onset times: the periods become segments in order
  if(Rt_max_changes < 0 | delay_max_changes < 0)
    stop("The maximum number of changes cannot be negative.")
  if(Rt_max_changes > 0) {
    period <- as.integer(pmin(seq_len(n_days), Rt_max_changes + 1))
    Rt_random_walk <- FALSE
  }
  if(delay_max_changes > 0) {
    if(!is.null(misdating))
      stop("Changes in the delay at unknown times cannot be combined with misdating.")
    delay_period <- as.integer(pmin(seq_len(n_days), delay_max_changes + 1))
  }
  # deaths among the reported cases
  is_deaths <- !is.null(deaths_df)
  obs_deaths <- integer(0)
  max_dead <- rep(0L, n_days)
  if(is_deaths) {
    if(!all(c("time_onset", "time_reported", "deaths_reported") %in% colnames(deaths_df)))
      stop("deaths_df must contain columns: time_onset, time_reported and deaths_reported.")
    if(cfr_max_changes < 0 | recording_max_changes < 0)
      stop("The maximum number of changes cannot be negative.")
    joined <- dplyr::left_join(
      df[, c("time_onset", "time_reported", "cases_reported")],
      deaths_df[, c("time_onset", "time_reported", "deaths_reported")],
      by=c("time_onset", "time_reported"))
    if(any(is.na(joined$deaths_reported)))
      stop("There must be reported deaths for each report of each onset time.")
    if(any(joined$deaths_reported > joined$cases_reported))
      stop("Reported deaths must not exceed reported cases.")
    obs_deaths <- as.integer(joined$deaths_reported)
    # wrong records shown as dead are allowed for where a count of deaths falls
    falls <- joined %>%
      dplyr::group_by(.data$time_onset) %>%
      dplyr::summarise(is_fall=any(diff(.data$deaths_reported) < 0))
    max_dead <- as.integer(ifelse(falls$is_fall, max_erroneous, 0))
  }
  # fatal cases reported at a different speed, or recorded for their day of
  # death
  is_fatal_delay <- !identical(fatal_delay, FALSE)
  fatal_delay_settings <- list(probability_same=0.5, sd=0.5)
  if(is.list(fatal_delay))
    fatal_delay_settings[names(fatal_delay)] <- fatal_delay
  dated <- dated_by_death_settings(n_days, dated_by_death)
  is_dated <- !is.null(dated)
  if((is_fatal_delay | is_dated) & !is_deaths)
    stop("deaths_df must be supplied if fatal cases are reported differently or dated by death.")
  if((is_fatal_delay | is_dated) & !is.null(misdating))
    stop("Misdating cannot be combined with fatal cases which are reported differently or dated by death.")
  inflow_cases <- rep(0, n_days)
  if(is_dated) {
    # cases with onset before the first day which, if fatal and dated by
    # death, die on each day: cases_before ends on the day before the first
    n_before <- length(cases_before)
    lags <- seq_len(dated$max_lag)
    inflow_cases <- purrr::map_dbl(seq_len(n_days), function(t) {
      index <- n_before + t - lags
      is_earlier <- lags >= t & index >= 1
      sum(cases_before[index[is_earlier]] * dated$lag_probabilities[is_earlier])
    })
  }
  unknown <- setdiff(names(priors), names(relaxed_priors()))
  if(length(unknown) > 0)
    stop(paste0("Unknown priors: ", paste(unknown, collapse=", "), "."))
  priors_all <- relaxed_priors()
  priors_all[names(priors)] <- priors
  if(is.na(priors_all$change_beta))
    priors_all$change_beta <- n_days

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
    is_Rt_changepoints=as.integer(Rt_max_changes > 0),
    n_delay_periods=max(delay_period),
    delay_period=delay_period,
    is_delay_changepoints=as.integer(delay_max_changes > 0),
    is_tail=as.integer(tail),
    is_weekday=as.integer(!is.null(weekday_first)),
    weekday_report=weekday_report,
    is_report_effects=as.integer(report_effects),
    n_reports=length(report_days),
    report_index=as.integer(report_index),
    is_misdating=as.integer(!is.null(misdating)),
    n_shifts=length(shifts),
    shifts=as.array(shifts),
    shift_probs=as.array(shift_probs),
    is_deaths=as.integer(is_deaths),
    obs_deaths=as.array(obs_deaths),
    n_cfr_segments=as.integer(ifelse(is_deaths, cfr_max_changes + 1, 1)),
    n_recording_segments=as.integer(ifelse(is_deaths, recording_max_changes + 1, 1)),
    max_dead=max_dead,
    is_recording_weekday=as.integer(is_deaths & recording_weekday),
    is_fatal_delay=as.integer(is_fatal_delay),
    fatal_delay_same_prob=fatal_delay_settings$probability_same,
    fatal_delay_sd=fatal_delay_settings$sd,
    is_dated=as.integer(is_dated),
    dated_fraction=if(is_dated) dated$fraction else rep(0, n_days),
    n_lags=if(is_dated) dated$max_lag else 0L,
    lag_probs=as.array(if(is_dated) dated$lag_probabilities else numeric(0)),
    inflow_cases=inflow_cases,
    is_dated_delay=as.integer(is_dated && !is.null(dated$reporting)),
    dated_delay_mean=if(is_dated && !is.null(dated$reporting)) dated$reporting$mean else 1,
    dated_delay_sd=if(is_dated && !is.null(dated$reporting)) dated$reporting$sd else 1),
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
      log_Rt_raw=as.array(stats::rnorm(data$n_periods, 0, 0.1)),
      Rt_p_change=as.array(rep(1 / data$n_days, data$is_Rt_changepoints)),
      delay_p_change=as.array(rep(1 / data$n_days, data$is_delay_changepoints)),
      Rt_step=as.array(rep(0.05, data$is_Rt_random_walk)),
      tail_weight=as.array(rep(0.05, data$is_tail)),
      tail_mean=as.array(rep(15, data$is_tail)),
      report_raw=as.array(rep(0, data$is_report_effects * data$n_reports)),
      report_sd=as.array(rep(0.2, data$is_report_effects)),
      pi_misdated=as.array(rep(0.05, data$is_misdating)),
      cfr=as.array(stats::runif(data$is_deaths * data$n_cfr_segments, 0.3, 0.5)),
      recording_mean=as.array(stats::runif(
        data$is_deaths * data$n_recording_segments, 10, 14)),
      recording_sd=as.array(stats::runif(
        data$is_deaths * data$n_recording_segments, 5, 8)),
      log_weekday_recording=rep(0, 7),
      wrong_dead=as.array(rep(0.5, data$is_deaths)),
      cfr_p_change=as.array(rep(
        1 / data$n_days, data$is_deaths * (data$n_cfr_segments > 1))),
      recording_p_change=as.array(rep(
        1 / data$n_days, data$is_deaths * (data$n_recording_segments > 1))),
      nu=as.array(rep(-0.1, data$is_deaths * data$is_fatal_delay)))
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
#' cases_true), 'Rt' (iteration, time_onset, Rt), 'delay' (iteration,
#' time_onset, mean, sd) and 'parameters' (a column
#' for each of the other parameters; where the delay differs between periods
#' its mean and sd have a column for each). If deaths are modelled there is
#' also 'deaths' (iteration, time_onset, deaths_true, cfr, delay_mean and
#' delay_sd, the last two being of the delay from onset to a death being
#' recorded). If some fatal cases are dated by death, 'cases' holds the true
#' cases by day of onset and 'cases_recorded' those recorded for each day.
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
                        "Rt_p_change", "delay_p_change", "Rt_n_changes",
                        "delay_n_changes", "weekday_effect", "report_effect",
                        "wrong_dead", "cfr_p_change", "recording_p_change",
                        "cfr_n_changes", "recording_n_changes",
                        "weekday_recording_effect", "fatal_delay_ratio")
  names_fit <- unique(sub("\\[.*", "", posterior::variables(fit$draws())))
  parameters <- fit$draws(intersect(names_parameters, names_fit),
                          format="draws_df") %>%
    dplyr::as_tibble() %>%
    dplyr::select(-dplyr::any_of(c(".chain", ".iteration"))) %>%
    dplyr::rename(iteration=".draw")
  # those of which there is only one are named without an index
  for(name in c("delay_mean", "delay_sd", "kappa", "Rt_step", "tail_weight",
                "tail_mean", "report_sd", "pi_misdated", "Rt_p_change",
                "delay_p_change", "wrong_dead", "cfr_p_change",
                "recording_p_change")) {
    indexed <- grep(paste0("^", name, "\\["), colnames(parameters), value=TRUE)
    if(length(indexed) == 1)
      colnames(parameters)[colnames(parameters) == indexed] <- name
  }
  out <- list(cases=f_long("cases_true", "cases_true"),
              Rt=f_long("Rt", "Rt"),
              delay=dplyr::left_join(f_long("delay_mean_day", "mean"),
                                     f_long("delay_sd_day", "sd"),
                                     by=c("iteration", "time_onset")),
              parameters=parameters)
  # deaths, if they are modelled
  if("cfr_day" %in% names_fit && ncol(fit$draws("cfr_day", format="draws_matrix")) > 0) {
    by <- c("iteration", "time_onset")
    out$deaths <- f_long("deaths_true", "deaths_true") %>%
      dplyr::left_join(f_long("cfr_day", "cfr"), by=by) %>%
      dplyr::left_join(f_long("recording_mean_day", "delay_mean"), by=by) %>%
      dplyr::left_join(f_long("recording_sd_day", "delay_sd"), by=by)
  }
  # where some fatal cases are dated by death, the cases are those by day of
  # onset and the cases recorded for each day are given too
  if("cases_onset" %in% names_fit &&
     ncol(fit$draws("cases_onset", format="draws_matrix")) > 0) {
    out$cases_recorded <- out$cases
    out$cases <- f_long("cases_onset", "cases_true")
  }
  out
}
