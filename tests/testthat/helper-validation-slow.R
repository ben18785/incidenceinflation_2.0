# Shared fits and checks for the recovery tests of the sensitivity analyses in the slow
# tier (C9 of sensitivity_extensions.tex: S6, S7).

# a fit of the deaths model to the S6 or S7 data with the form of the rate given
s_wrong_fit <- function(dat, niter, seed, revisions, initial_revisions) {
  snap <- dat$cases; snap$Rt_index <- 1L
  last <- tapply(snap$cases_reported, snap$time_onset, function(x) x[length(x)])
  init_cases <- tibble::tibble(time_onset = 1:dat$Tn, cases_true = as.integer(ceiling(1.3 * last) + 5))
  set.seed(seed)
  mcmc_single(
    niter, snap,
    priors = list(Rt = list(shape = 1, rate = 1),
                  reporting = list(mean_mu = 5, mean_sigma = 3, sd_mu = 3, sd_sigma = 2),
                  max_cases = 300, Rt_changepoints = list(), reporting_changepoints = list(),
                  revisions = revisions, max_erroneous = 8),
    serial_parameters = list(mean = 15.3, sd = 9.3),
    initial_cases_true = init_cases,
    initial_reporting_parameters = list(mean = 5, sd = 2.5),
    initial_Rt = tibble::tibble(Rt_index = 1L, Rt = 1),
    print_to_screen = FALSE, cases_before = dat$hist,
    initial_revision_parameters = initial_revisions,
    deaths = list(snapshot = dat$deaths,
                  settings = list(cfr_changepoints = FALSE, initial_cfr = 0.45,
                                  onset_to_death = list(mean = 9.3, sd = 5.5, max_lag = 36),
                                  initial_delay = list(mean = 4, sd = 2),
                                  delay_prior = list(mean_mu = 5, mean_sigma = 4, sd_mu = 3, sd_sigma = 2),
                                  delay_changepoints = list(min_length = 20, first_length = 20),
                                  triggered = list(alpha = 1, beta = 1, initial = 1, step = 0.1, n_steps = 0),
                                  wrong_dead_alpha = 1, wrong_dead_beta = 1)))
}

# four chains in parallel where the platform allows it
s_wrong_chains <- function(f_chain, n_chains = 4) {
  if (.Platform$OS.type == "unix")
    parallel::mclapply(seq_len(n_chains), f_chain,
                       mc.cores = min(n_chains, parallel::detectCores()))
  else lapply(seq_len(n_chains), f_chain)
}

# the kept draws of a quantity in each chain (the second half), as a list of vectors
s_wrong_kept <- function(fits, f_values, niter)
  lapply(fits, function(fit) {
    x <- f_values(fit)
    x[seq_along(x) > niter / 2]
  })

s_wrong_check <- function(chains, truth, label) {
  rhat <- split_rhat(unlist(chains), rep(seq_along(chains), lengths(chains)))
  x <- unlist(chains)
  ci <- stats::quantile(x, c(0.005, 0.995), names = FALSE)
  note("%s: truth %.4g, posterior mean %.4g, 99%% interval [%.4g, %.4g], split R-hat %.3f",
       label, truth, mean(x), ci[1], ci[2], rhat)
  expect(is.finite(rhat) && rhat < 1.2, sprintf("%s: split R-hat %.3f is not below 1.2", label, rhat))
  expect(truth >= ci[1] && truth <= ci[2],
         sprintf("%s: truth %.4g outside the 99%% interval [%.4g, %.4g] (posterior mean %.4g)",
                 label, truth, ci[1], ci[2], mean(x)))
}

# the total of the true cases of the last 14 onset days, by iteration
s_wrong_last14 <- function(fit, Tn) {
  cs <- fit$cases[col_of(fit$cases, "time_onset") >= Tn - 13, ]
  as.numeric(tapply(col_of(cs, "cases_true"), col_of(cs, "iteration"), sum))
}
