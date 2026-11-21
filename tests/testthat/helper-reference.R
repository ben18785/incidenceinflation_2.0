# Reference fits of the main model (test-reference.R), recorded before the
# changes for the sensitivity analyses and checked against after each of
# them. The settings are those of the fits to the reports of the Bundibugyo
# outbreak, on small simulated data sets: the exact form of fatal cases dated
# by death with an estimated fraction (with and without the joint move along
# the ridge), the backlog of deaths with a holding period, a time from onset
# to death, fatal cases reported when their death is recorded, wrong records
# with effects of the day of the week, the heavier tail of the reporting
# delay and a report correcting some onset times.

reference_data <- function(seed) {
  set.seed(seed)
  n_days <- 35
  serial_parameters <- list(mean=8, sd=5)
  cases <- true_cases(n_days + 40, function(t) 1.05 - 0.2 * (t > 60), kappa=50,
                      serial_parameters,
                      initial_parameters=list(mean=4, length=10))
  report_days <- sort(unique(c(10, cumsum(sample(1:3, 30, replace=TRUE)) + 10)))
  report_days <- report_days[report_days <= n_days + 5]
  time_corrected <- report_days[length(report_days) - 3]
  corrected <- 5:10
  sim <- simulate_cases_deaths(
    cases[41:(n_days + 40)], report_days, list(mean=6, sd=4), 0.45,
    list(mean=3, sd=2.5), triggered=1,
    onset_to_death=list(mean=9.3, sd=5.5),
    dated_by_death=list(fraction=0.3),
    recording_backlog=list(share=0.2, rate=0.05, hold_days=5),
    wrong=list(eta=0.05, rho=0.3, delta=0.3),
    correction=list(time_reported=time_corrected, onset_times=corrected,
                    wrong=0.1, missing=0.03, misdated=0.1))
  parts <- split_corrected_reports(sim$snapshots, time_corrected, corrected)
  list(snapshots=parts$corrected, uncorrected=parts$uncorrected,
       cases_before=cases[1:40], serial_parameters=serial_parameters,
       initial_cases=parts$corrected %>%
         dplyr::group_by(time_onset) %>%
         dplyr::summarise(cases_true=ceiling(1.3 * max(cases_reported)) + 3))
}

reference_fit <- function(data, ridge, niterations, seed, resume_from=NULL) {
  set.seed(seed)
  mcmc_single(
    niterations=niterations,
    snapshot_with_Rt_index_df=data$snapshots %>%
      dplyr::select(time_onset, time_reported, cases_reported) %>%
      dplyr::mutate(Rt_index=1),
    priors=list(Rt=list(shape=1, rate=1),
                reporting=list(mean_mu=8, mean_sigma=5, sd_mu=6, sd_sigma=4),
                max_cases=400, overdispersion=list(mean=20, sd=20),
                revisions=list(eta_mu=0.1, eta_sigma=0.1, rho_alpha=1, rho_beta=1),
                max_erroneous=8,
                reporting_effects=list(tail_weight_alpha=1, tail_weight_beta=19,
                                       tail_mean_mu=40, tail_mean_sigma=30),
                Rt_changepoints=list(),
                reporting_changepoints=list(min_length=14, n_partition_moves=5)),
    serial_parameters=data$serial_parameters,
    initial_cases_true=data$initial_cases,
    initial_reporting_parameters=list(mean=6, sd=4),
    initial_Rt=dplyr::tibble(Rt_index=1, Rt=1),
    reporting_metropolis_parameters=list(mean_step=0.3, sd_step=0.25),
    print_to_screen=FALSE, is_negative_binomial=TRUE, initial_overdispersion=20,
    overdispersion_metropolis_sd=8,
    initial_revision_parameters=list(eta=0.02, rho=0.2, weekday=TRUE,
                                     weekday_first=1,
                                     tail=list(weight=0.05, mean=40)),
    revision_metropolis_parameters=list(eta_step=0.3, rho_step=0.3, pi_step=0.3),
    cases_before=data$cases_before, n_parameter_steps=2, return_state=TRUE,
    resume_from=resume_from,
    deaths=list(
      snapshot=data$snapshots %>%
        dplyr::select(time_onset, time_reported, deaths_reported),
      settings=list(
        onset_to_death=list(mean=9.3, sd=5.5),
        dated_by_death=c(list(fraction=list(alpha=1, beta=1, initial=0.25),
                              exact=TRUE),
                         if(ridge) list(ridge=TRUE)),
        initial_cfr=0.4, initial_delay=list(mean=15, sd=13),
        delay_prior=list(mean_mu=4, mean_sigma=4, sd_mu=4, sd_sigma=4),
        delay_changepoints=list(min_length=14, n_partition_moves=5,
                                first_length=20),
        triggered=list(initial=1, n_steps=0L), fatal_delay=FALSE,
        recording_backlog=list(hold_days=5L)),
      uncorrected=data$uncorrected))
}

# the cases of the reference fits: two data sets, each without and with the
# joint move along the ridge
reference_cases <- function()
  expand.grid(seed=c(101, 202), ridge=c(FALSE, TRUE))

# what is compared: the draws (less the times taken) and the state at the end
reference_summary <- function(fit) {
  draws <- fit[setdiff(names(fit), c("state", "iterations_done"))]
  state <- fit$state$variables
  list(draws=draws,
       deaths_state=state$deaths_state[setdiff(names(state$deaths_state),
                                               c("timings"))],
       cases_true=state$deaths_cases_true,
       reporting_logl=state$reporting_state$logl,
       log_likelihood=deaths_log_likelihood(state$deaths_state,
                                            state$deaths_cases_true))
}

reference_iterations <- 8

# draws from a fixed seed are compared exactly only on the platform they were
# recorded on, since another can round differently
reference_platform <- function()
  paste(R.version$platform, R.version$version.string, sessionInfo()$BLAS)
