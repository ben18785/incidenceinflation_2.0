skip_validation("slow")

# S9 (slow tier): recovery with a fraction dated by death which varies by onset week
# (S9 of the sensitivity analyses, test S9.7 and C9 of its test plan), on one dataset from
# the independent simulator (helper-validation-sim.R, sim_reports_deaths_dated).
#
# Data: T = 80 onset days; Poisson renewal with R = 1.2 on days 1-30 and 0.9 on days
# 31-80, history of 30 a day; case delay gamma(5, 2.5); constant CFR 0.45; time from onset
# to death p_die (gamma 9.3, 5.5 on (z-1, z], z = 1..36); recording delay from death
# gamma(3, 2); triggered reports with tau = 1 (held there); wrong records eta = 0.05,
# rho = 0.3, shown dead with probability delta = 0.3; no backlog. A fatal case is
# misdated with probability 0.2 for onset weeks 1-6 (days 1-42) and 0.5 for weeks 7-12
# (days 43-80), and those of the history with 0.2. Reports with gaps drawn from 1-7 days
# from day 2 to day 80 (about 25).
# Fit: four chains of mcmc_single with misdating in the exact form and the fraction a
# random walk over onset weeks (dated_by_death fraction = list(walk = TRUE)), the
# recording delay constant; at most 4 wrong records present (the largest fall in the data
# is 2; the exact form is slow with many more). About 0.5-1 s an iteration on a Mac.
# Checks: max split R-hat < 1.2 over the CFR, the fraction of onset days 28 and 56, the
# sd of the walk's steps and the true cases of the last 14 onset days in total; the
# central 99% interval (second half of each chain) contains the truth for the CFR, the
# fraction before the step (onset day 28, week 4) and after it (onset day 56, week 8),
# and the true cases of the last 14 onset days in total. A single dataset, so this is
# not a calibration test; 99% intervals keep false alarms low.

s9_data <- function() {
  set.seed(8901)
  Tn <- 80; hist <- rep(30, 40); omega <- ref_serial_weights(15.3, 9.3, 40)
  R <- c(rep(1.2, 30), rep(0.9, 50))
  I <- sim_renewal(R, hist, omega)
  rd <- sim_report_days(2, Tn, 1:7)
  fraction <- c(rep(0.2, 42), rep(0.5, 38))
  sim <- sim_reports_deaths_dated(I, rd, list(mean = 5, sd = 2.5), rep(0.45, Tn), ref_pdie(),
                                  rep(3, Tn), rep(2, Tn), fraction, history = hist, tau = 1,
                                  eta = 0.05, rho = 0.3, delta = 0.3)
  list(Tn = Tn, hist = hist, R = R, I = I, rd = rd, fraction = fraction, sim = sim)
}

s9_chain <- function(dat, niter, seed, time_limit) {
  snap <- dat$sim$cases; snap$Rt_index <- 1L
  last <- tapply(snap$cases_reported, snap$time_onset, function(x) x[length(x)])
  init_cases <- tibble::tibble(time_onset = 1:dat$Tn, cases_true = as.integer(ceiling(1.3 * last) + 5))
  set.seed(seed)
  mcmc_single(
    niter, snap,
    priors = list(Rt = list(shape = 1, rate = 1),
                  reporting = list(mean_mu = 5, mean_sigma = 3, sd_mu = 3, sd_sigma = 2),
                  max_cases = 400, Rt_changepoints = list(), reporting_changepoints = list(),
                  revisions = list(eta_mu = 0.05, eta_sigma = 0.05, rho_alpha = 1, rho_beta = 4),
                  max_erroneous = 4),
    serial_parameters = list(mean = 15.3, sd = 9.3),
    initial_cases_true = init_cases,
    initial_reporting_parameters = list(mean = 5, sd = 2.5),
    initial_Rt = tibble::tibble(Rt_index = 1L, Rt = 1),
    print_to_screen = FALSE, cases_before = dat$hist,
    initial_revision_parameters = list(eta = 0.05, rho = 0.3),
    time_limit = time_limit,
    deaths = list(snapshot = dat$sim$deaths,
                  settings = list(cfr_changepoints = FALSE, initial_cfr = 0.45,
                                  onset_to_death = list(mean = 9.3, sd = 5.5, max_lag = 36),
                                  initial_delay = list(mean = 3, sd = 2),
                                  delay_prior = list(mean_mu = 3, mean_sigma = 3, sd_mu = 2, sd_sigma = 2),
                                  delay_changepoints = FALSE,
                                  dated_by_death = list(fraction = list(walk = TRUE), exact = TRUE),
                                  triggered = list(alpha = 1, beta = 1, initial = 1, step = 0.1, n_steps = 0),
                                  wrong_dead_alpha = 1, wrong_dead_beta = 1)))
}

test_that("S9 a fraction dated by death which varies by onset week is recovered (99% intervals)", {
  skip_validation("slow")
  dat <- s9_data()
  note("S9 data: %d reports (days %s); %d cases, %d fatal, %d misdated (%d of %d fatal in days 1-42, %d of %d in 43-80)",
       length(dat$rd), paste(dat$rd, collapse = " "), sum(dat$I), sum(dat$sim$fatal),
       sum(dat$sim$misdated), sum(dat$sim$misdated[1:42]), sum(dat$sim$fatal[1:42]),
       sum(dat$sim$misdated[43:80]), sum(dat$sim$fatal[43:80]))
  niter <- 5000; n_chains <- 4
  t0 <- proc.time()[["elapsed"]]
  # in parallel where there are the cores; each chain stops at the time limit (minutes)
  fits <- parallel::mclapply(seq_len(n_chains), function(chain)
    s9_chain(dat, niter, 8902 + chain, time_limit = 270),
    mc.cores = min(n_chains, parallel::detectCores()))
  note("S9: %d chains took %.1f min", n_chains, (proc.time()[["elapsed"]] - t0) / 60)
  for (fit in fits) if (inherits(fit, "try-error")) stop(fit)
  done <- min(vapply(fits, function(fit) max(fit$deaths$iteration), 0))
  note("S9: %d iterations in each chain", done)
  expect(done >= 500, sprintf("S9: only %d iterations done in each chain", done))

  # the second half of each chain, to the iterations all have done
  f_kept <- function(df) df[col_of(df, "iteration") > done / 2 & col_of(df, "iteration") <= done, ]
  f_series <- function(element, t, col)
    lapply(fits, function(fit) {
      x <- f_kept(fit[[element]])
      if (!is.null(t)) x <- x[col_of(x, "time_onset") == t, ]
      col_of(x, col)[order(col_of(x, "iteration"))]
    })
  f_last_cases <- function(fit) {
    x <- f_kept(fit$cases)
    x <- x[col_of(x, "time_onset") >= dat$Tn - 13, ]
    as.numeric(tapply(col_of(x, "cases_true"), col_of(x, "iteration"), sum))
  }
  series <- list(
    "CFR" = f_series("deaths", 1, "cfr"),
    "fraction, onset day 28" = f_series("deaths", 28, "dated_fraction"),
    "fraction, onset day 56" = f_series("deaths", 56, "dated_fraction"),
    "sd of the steps of the walk" = f_series("deaths_parameters", NULL, "dated_fraction_step_sd"),
    "true cases of the last 14 onset days" = lapply(fits, f_last_cases))
  rhat <- vapply(series, function(x) split_rhat(unlist(x), rep(seq_along(x), lengths(x))), 0)
  for (name in names(rhat)) note("S9 split R-hat, %s: %.3f", name, rhat[[name]])
  rhat[is.na(rhat)] <- Inf   # too few draws to tell
  expect(max(rhat) < 1.2, sprintf("S9: max split R-hat %.3f (%s) is not below 1.2",
                                  max(rhat), names(rhat)[which.max(rhat)]))

  ci_check <- function(x, truth, label) {
    x <- unlist(x)
    ci <- stats::quantile(x, c(0.005, 0.995), names = FALSE)
    note("S9 %s: truth %.4g, posterior mean %.4g, 99%% interval [%.4g, %.4g]", label, truth, mean(x), ci[1], ci[2])
    expect(truth >= ci[1] && truth <= ci[2],
           sprintf("S9 %s: truth %.4g outside the 99%% interval [%.4g, %.4g] (posterior mean %.4g)",
                   label, truth, ci[1], ci[2], mean(x)))
  }
  ci_check(series[["CFR"]], 0.45, "CFR")
  ci_check(series[["fraction, onset day 28"]], 0.2, "fraction before the step, onset day 28")
  ci_check(series[["fraction, onset day 56"]], 0.5, "fraction after the step, onset day 56")
  ci_check(series[["true cases of the last 14 onset days"]], sum(dat$I[(dat$Tn - 13):dat$Tn]),
           "true cases of the last 14 onset days")
})
