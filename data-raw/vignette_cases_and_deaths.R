# Precomputes the fit shown in vignettes/cases_and_deaths.Rmd, which takes
# too long to run when the vignette is built. Run by the workflow
# vignette-precompute.yaml:
#   Rscript data-raw/vignette_cases_and_deaths.R chain <chain> <iterations>
# fits one chain and saves its draws to vignette_chain_<chain>.rds;
#   Rscript data-raw/vignette_cases_and_deaths.R combine
# summarises the chains into vignettes/cases_and_deaths_results.rds.
#
# The simulation is the same as in the vignette, which checks that the
# reports it simulates are those fitted here.

library(incidenceinflation)
library(dplyr)

simulate_vignette_outbreak <- function() {
  set.seed(1)
  n_days <- 60
  serial <- list(mean = 5, sd = 3)
  Rt <- function(t) ifelse(t < 30, 1.2, 0.9)
  # 40 days before the first onset day give the renewal its history
  cases_all <- true_cases(n_days + 40, function(t) Rt(t - 40), kappa = 1000, serial,
                          initial_parameters = list(mean = 4, length = 10))
  cfr <- ifelse(seq_len(n_days) < 30, 0.5, 0.3)
  sim <- simulate_cases_deaths(
    cases_all[41:(n_days + 40)],
    report_days = seq(3, n_days, by = 3),
    reporting_parameters = list(mean = 7, sd = 4),
    cfr = cfr,
    death_parameters = list(mean = 4, sd = 3),
    onset_to_death = list(mean = 9, sd = 5),
    dated_by_death = list(fraction = 0.25))
  c(sim, list(n_days = n_days, serial = serial, cases_before = cases_all[1:40],
              Rt = Rt(seq_len(n_days)), cfr = cfr))
}

fit_vignette_chain <- function(outbreak, chain, niterations) {
  set.seed(100 + chain)
  snapshots <- outbreak$snapshots
  # chains start from different values: the largest count reported for each
  # onset day, inflated for recent days by different amounts, and different
  # case fatality ratios and delays
  inflation <- c(2, 1, 3.5, 1.5)[chain]
  largest <- snapshots %>%
    group_by(time_onset) %>%
    summarise(cases_true = max(cases_reported))
  days_before_end <- outbreak$n_days - largest$time_onset
  initial_cases <- largest %>%
    mutate(cases_true = ceiling(cases_true * (1 + inflation * exp(-days_before_end / 4)) +
                                  5 * (days_before_end < 5)))
  mcmc_single(
    niterations = niterations,
    snapshot_with_Rt_index_df = snapshots %>%
      select(time_onset, time_reported, cases_reported) %>%
      mutate(Rt_index = 1),
    priors = list(
      Rt = list(shape = 1, rate = 1),
      reporting = list(mean_mu = 7, mean_sigma = 10, sd_mu = 5, sd_sigma = 5),
      max_cases = ceiling(4 * max(initial_cases$cases_true) + 50),
      Rt_changepoints = list(),
      reporting_changepoints = list(min_length = 14, n_partition_moves = 20)),
    serial_parameters = outbreak$serial,
    initial_cases_true = initial_cases,
    initial_reporting_parameters = list(list(mean = 7, sd = 4), list(mean = 11, sd = 7),
                                        list(mean = 4, sd = 2.5), list(mean = 9, sd = 5))[[chain]],
    initial_Rt = tibble(Rt_index = 1, Rt = 1),
    print_to_screen = FALSE,
    cases_before = outbreak$cases_before,
    deaths = list(
      snapshot = snapshots %>% select(time_onset, time_reported, deaths_reported),
      settings = list(
        initial_cfr = c(0.4, 0.2, 0.7, 0.3)[chain],
        onset_to_death = list(mean = 9, sd = 5),
        dated_by_death = list(fraction = list(alpha = 1, beta = 1, initial = 0.25),
                              exact = TRUE),
        initial_delay = list(list(mean = 4, sd = 3), list(mean = 8, sd = 5),
                             list(mean = 2, sd = 2), list(mean = 6, sd = 4))[[chain]])))
}

summarise_vignette_chains <- function(outbreak, files) {
  chains <- lapply(files, readRDS)
  f_bind <- function(name) bind_rows(lapply(chains, function(x) x$draws[[name]]))
  f_days <- function(df, column) df %>%
    group_by(time_onset) %>%
    summarise(median = median(.data[[column]]), lower = quantile(.data[[column]], 0.05),
              upper = quantile(.data[[column]], 0.95), .groups = "drop")
  deaths <- f_bind("deaths")
  f_rhat <- function(df, column) {
    m <- tidyr::pivot_wider(select(df, "chain", "iteration", all_of(column)),
                            names_from = "chain", values_from = all_of(column))
    tryCatch(posterior::rhat(as.matrix(select(m, -"iteration"))), error = function(e) NA)
  }
  parameters <- f_bind("deaths_parameters")
  list(
    snapshots = outbreak$snapshots,
    iterations = vapply(chains, function(x) x$niterations, numeric(1)),
    cfr = f_days(deaths, "cfr"),
    Rt = f_days(f_bind("Rt"), "Rt"),
    cases = f_days(f_bind("cases"), "cases_true"),
    deaths = f_days(deaths, "deaths_true"),
    dated_fraction = parameters$dated_fraction,
    n_cfr_changes = parameters$n_cfr_changes,
    rhat = c(dated_fraction = f_rhat(parameters, "dated_fraction"),
             delay_centre_mean = f_rhat(parameters, "delay_centre_mean"),
             Rt_last = f_rhat(filter(f_bind("Rt"), time_onset == outbreak$n_days), "Rt"),
             cfr_first = f_rhat(filter(deaths, time_onset == 1), "cfr"),
             cfr_last = f_rhat(filter(deaths, time_onset == outbreak$n_days), "cfr")))
}

args <- commandArgs(TRUE)
outbreak <- simulate_vignette_outbreak()
if(args[1] == "chain") {
  chain <- as.integer(args[2])
  niterations <- as.integer(args[3])
  fit <- fit_vignette_chain(outbreak, chain, niterations)
  # the second half of the iterations, thinned to 500 draws
  kept <- unique(round(seq(niterations / 2 + 1, niterations, length.out = 500)))
  draws <- lapply(fit[c("cases", "Rt", "deaths", "deaths_parameters")], function(df)
    df %>% filter(iteration %in% kept) %>% mutate(chain = chain))
  saveRDS(list(niterations = niterations, draws = draws),
          sprintf("vignette_chain_%d.rds", chain))
} else if(args[1] == "combine") {
  files <- sort(list.files(".", "^vignette_chain_[0-9]+[.]rds$", recursive = TRUE, full.names = TRUE))
  out <- summarise_vignette_chains(outbreak, files)
  print(out$rhat)
  saveRDS(out, "vignettes/cases_and_deaths_results.rds", compress = "xz")
}
