# Shared settings and output readers for the medium tier (M*).

m_serial <- list(mean = 15.3, sd = 9.3)
# (m_omega, the serial weights, is set in each M test file: helper-validation-spec.R loads after this file)

# Values of column 'col' at onset day t, one per iteration, ordered by iteration,
# from the iterations after the first half.
m_series <- function(df, t, col, niter) {
  x <- df[col_of(df, "time_onset") == t & col_of(df, "iteration") > niter / 2, , drop = FALSE]
  col_of(x, col)[order(col_of(x, "iteration"))]
}
m_kept <- function(df, col, niter) {
  x <- df[col_of(df, "iteration") > niter / 2, , drop = FALSE]
  col_of(x, col)[order(col_of(x, "iteration"))]
}

# C5 of the sensitivity analyses (sensitivity_extensions.tex): Geweke's successive-
# conditional simulator for the wrong records of the deaths model. The true cases I of
# each onset day are held and none is fatal (CFR 0), so the counts of deaths are the
# wrong records shown dead. Each iteration simulates the reports given the current
# parameters (f_simulate, a function of the revision parameters and delta, from the
# independent simulator), then updates the split of
# the reports (deaths_paths_step, exact given the parameters) and the parameters of
# the wrong records given it (deaths_wrong_step). The parameters then have their prior
# as stationary distribution. f_record gives the values kept from the revision
# parameters and the state.
m_wrong_geweke <- function(niter, I, delay, revision, delta, priors, settings, f_simulate,
                           f_record, max_erroneous) {
  n_days <- length(I)
  f_reports <- function(sim)
    dplyr::left_join(sim$cases, sim$deaths, by = c("time_onset", "time_reported"))
  state <- initial_deaths_state(f_reports(f_simulate(revision, delta)),
                                # only to start: the counts include wrong records
                                tibble::tibble(time_onset = seq_len(n_days),
                                               cases_true = I + max_erroneous),
                                settings)
  state$cfr <- rep(0, n_days)
  state$deaths_true <- rep(0, n_days)
  state$delta <- delta
  reporting <- tibble::tibble(mean = rep(delay$mean, n_days), sd = delay$sd)
  out <- vector("list", niter)
  for (m in seq_len(niter)) {
    d <- f_reports(f_simulate(revision, state$delta))
    state$reports <- split(as.data.frame(d[, c("time_reported", "cases_reported", "deaths_reported")]),
                           d$time_onset)[as.character(seq_len(n_days))]
    state <- deaths_paths_step(state, I, reporting, revision, max_erroneous)
    updated <- deaths_wrong_step(state, I, reporting, revision, settings, priors)
    state <- updated$state
    revision <- updated$revision_parameters
    out[[m]] <- f_record(revision, state)
  }
  dplyr::bind_rows(out)
}
