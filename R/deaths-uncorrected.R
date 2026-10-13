# Reports made before a report which corrected the records of some onset
# times. For each such onset time the reports are split into two series: the
# reports from the correction on, which are handled like those of any other
# onset time, and the reports before it (the "uncorrected series"), which are
# handled here.
#
# The uncorrected series of an onset time follows the same model as any other
# series of reports, with the delays of that onset time, but has its own
# numbers of cases and of fatal cases. These have flat priors and enter
# neither the renewal model nor the case fatality ratio, so the series informs
# only the delays and the parameters of wrong records.

#' Splits snapshots at a report which corrected the records of some onset
#' times
#'
#' @param snapshot_df a tibble with columns 'time_onset' and 'time_reported'
#' and any others
#' @param time_corrected the time of the report which made the correction
#' @param onset_times the onset times whose records it corrected
#'
#' @return a named list of two tibbles: 'corrected' (all reports of the other
#' onset times and, for the onset times corrected, the reports from the
#' correction on) and 'uncorrected' (the reports of the onset times corrected
#' which were made before the correction)
#' @export
#' @importFrom rlang .data
split_corrected_reports <- function(snapshot_df, time_corrected, onset_times) {
  is_uncorrected <- snapshot_df$time_onset %in% onset_times &
    snapshot_df$time_reported < time_corrected
  list(corrected=snapshot_df[!is_uncorrected, ],
       uncorrected=snapshot_df[is_uncorrected, ])
}

#' Initial state of the uncorrected series
#'
#' @param snapshot_df a tibble with columns 'time_onset', 'time_reported',
#' 'cases_reported' and 'deaths_reported': the reports made before the
#' correction, for the onset times it corrected
#' @inheritParams deaths_delay_by_onset
#' @param settings as returned by \code{deaths_settings}
#'
#' @return a named list with an element for each series in 'day' (which of
#' the onset times of the state it belongs to), 'onset_times', 'reports',
#' 'cases_true' and 'deaths_true' (what the uncorrected records would have
#' come to) and 'cfr' (the share of those cases which are fatal, which is
#' free for each series)
#' @importFrom rlang .data
initial_uncorrected <- function(snapshot_df, state, settings) {
  if(!all(c("time_onset", "time_reported", "cases_reported", "deaths_reported") %in%
          colnames(snapshot_df)))
    stop("The uncorrected reports must contain columns: 'time_onset', 'time_reported', 'cases_reported' and 'deaths_reported'.")
  snapshot <- snapshot_df %>%
    dplyr::select("time_onset", "time_reported", "cases_reported", "deaths_reported") %>%
    dplyr::arrange(.data$time_onset, .data$time_reported)
  onset_times <- sort(unique(snapshot$time_onset))
  if(!all(onset_times %in% state$onset_times))
    stop("The uncorrected reports must be for onset times which are in the snapshots.")
  if(any(snapshot$deaths_reported > snapshot$cases_reported))
    stop("Reported deaths must not exceed reported cases.")
  reports <- split(as.data.frame(snapshot[, c("time_reported", "cases_reported",
                                             "deaths_reported")]),
                   snapshot$time_onset)[as.character(onset_times)]
  cases_true <- unname(purrr::map_dbl(reports, ~max(.$cases_reported)))
  largest <- unname(purrr::map_dbl(reports, ~max(.$deaths_reported)))
  list(day=match(onset_times, state$onset_times),
       onset_times=onset_times,
       reports=reports,
       cases_true=cases_true,
       deaths_true=pmin(cases_true,
                        pmax(largest, round(settings$initial_cfr * cases_true))),
       cfr=rep(settings$initial_cfr, length(onset_times)),
       paths=NULL, survivals=NULL, wrong_units=NULL, wrong=NULL)
}

#' The uncorrected series as a state of the model of deaths
#'
#' The result can be passed to the functions which act on the reports of each
#' onset time, with the parameters shared with the other reports taken from
#' the state.
#'
#' @inheritParams deaths_delay_by_onset
#' @param nu the log of the ratio of the reporting delay of fatal cases to
#' that of others
#'
#' @return NULL if there are no uncorrected series; otherwise a state with an
#' element for each uncorrected series where the state has one for each onset
#' time
uncorrected_state <- function(state, nu=state$nu) {
  uncorrected <- state$uncorrected
  if(is.null(uncorrected))
    return(NULL)
  c(uncorrected[c("day", "onset_times", "reports", "cases_true", "deaths_true",
                  "cfr", "paths", "survivals", "wrong_units", "wrong")],
    list(delays=deaths_delay_by_onset(state)[uncorrected$day, ],
         recording_effects=state$recording_effects,
         nu=nu, triggered=state$triggered, delta=state$delta, dated=NULL))
}

#' Rows of the reporting delay by onset time for the uncorrected series
#'
#' @inheritParams deaths_paths_step
#'
#' @return a tibble with a row for each uncorrected series
uncorrected_reporting <- function(state, reporting_by_onset) {
  reporting_by_onset[state$uncorrected$day, c("mean", "sd")]
}

#' Updates the uncorrected series
#'
#' For each series: the number of fatal cases is updated with the split of
#' the reported cases summed out, and the split is drawn, as in
#' \code{deaths_paths_step}; the number of cases is drawn from its conditional
#' distribution with the number of fatal cases summed out, and the number of
#' fatal cases is drawn given it; and the share of cases which are fatal is
#' drawn from its beta conditional distribution. That share has a uniform
#' prior for each series, so that the number of fatal cases is uniform given
#' the number of cases; the number of cases is uniform up to max_cases.
#'
#' @inheritParams deaths_paths_step
#' @param max_cases the largest possible number of cases
#'
#' @return the state with element 'uncorrected' updated
deaths_uncorrected_step <- function(state, reporting_by_onset,
                                    revision_parameters=NULL, max_erroneous=3,
                                    max_cases=1000) {
  if(is.null(state$uncorrected))
    return(state)
  sub <- uncorrected_state(state)
  sub <- deaths_paths_step(sub, sub$cases_true,
                           uncorrected_reporting(state, reporting_by_onset),
                           revision_parameters, max_erroneous)
  for(k in seq_along(sub$onset_times)) {
    f_logp <- deaths_cases_logp_function(sub, sub$onset_times[k],
                                         revision_parameters)
    cases_min <- attr(f_logp, "min_cases")
    values <- cases_min:max(max_cases, cases_min)
    logp <- f_logp(values)
    sub$cases_true[k] <- values[sample.int(length(values), 1,
                                           prob=exp(logp - max(logp)))]
  }
  sub <- deaths_true_step(sub, sub$cases_true)
  sub$cfr <- stats::rbeta(length(sub$cfr), 1 + sub$deaths_true,
                          1 + sub$cases_true - sub$deaths_true)
  names_updated <- c("cases_true", "deaths_true", "cfr", "paths", "survivals",
                     "wrong_units", "wrong")
  state$uncorrected[names_updated] <- sub[names_updated]
  state
}

#' The values held for the uncorrected series, for storage
#'
#' @inheritParams deaths_delay_by_onset
#' @param iteration the iteration
#'
#' @return a tibble with a row for each uncorrected series
uncorrected_draws <- function(state, iteration) {
  dplyr::tibble(iteration=iteration,
                time_onset=state$uncorrected$onset_times,
                cases_uncorrected=state$uncorrected$cases_true,
                deaths_uncorrected=state$uncorrected$deaths_true)
}
