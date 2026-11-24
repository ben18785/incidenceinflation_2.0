# S3 (slow tier): recovery of the main model's misdating (fatal cases dated by death in
# the exact form with a single estimated fraction, as the fits to the Bundibugyo reports)
# on three datasets from the independent simulator (helper-validation-sim.R,
# sim_reports_deaths_dated), each its own job (test-validation-S03-1, -2 and -3). The
# recovery tests of S5 and S9 (sensitivity analyses since dropped) found a second mode in
# the exact form, with the gamma delay from death to recording collapsed (shape near 0)
# and too few fatal cases misdated, which chains from cold starts
# reach and do not leave; S3 asks whether the main model, from starts as the real fits
# use, finds the truth.
#
# Data: T = 80 onset days; Poisson renewal with R = 1.2 on days 1-30 and 0.9 on days
# 31-80, history of 30 a day; case delay gamma(5, 2.5); constant CFR 0.45; time from onset
# to death p_die (gamma 9.3, 5.5 on (z-1, z], z = 1..36); recording delay from death
# gamma(3, 2); triggered reports with tau = 1 (held there); wrong records eta = 0.05,
# rho = 0.3, shown dead with probability delta = 0.3; no backlog. A fatal case is misdated
# with probability 0.35 (those of the history too; the real fits estimate about 0.38).
# Reports with gaps drawn from 1-7 days from day 2 to day 80 (about 25).
# Fit: four chains of mcmc_single in parallel, each stopped at 5 h (GitHub's limit for a
# job is 6 h), with the settings of the real fits' misdating: the exact form, the
# fraction with a uniform prior, the joint move along the ridge, the recording delay
# constant with the real fits' prior (mean_mu = 4, mean_sigma = 4, sd_mu = 4,
# sd_sigma = 4), starting from a mean of 3 and sd of 2.5; cold starts, the true cases
# 1.1-1.5 times the last reported and the fraction 0.1, 0.25 (the real fits' start), 0.5
# and 0.7; at most 4 wrong records present (as the S5 and S9 tests had).
# Checks: the central 99% interval (second half of each chain) contains the truth for the
# CFR, the fraction misdated, the mean delay from death to recording and the true cases
# of the last 14 onset days in total; the split R-hat of each is noted. The chain means
# and the shape of the recording delay (mean^2 / sd^2, 2.25 in truth) are noted, to tell
# the collapsed mode. Each dataset alone is not a calibration test; 99% intervals keep
# false alarms low.

s3_seeds <- function() c(8301, 8311, 8321)

s3_data <- function(replicate) {
  set.seed(s3_seeds()[replicate])
  Tn <- 80; hist <- rep(30, 40); omega <- ref_serial_weights(15.3, 9.3, 40)
  R <- c(rep(1.2, 30), rep(0.9, 50))
  I <- sim_renewal(R, hist, omega)
  rd <- sim_report_days(2, Tn, 1:7)
  sim <- sim_reports_deaths_dated(I, rd, list(mean = 5, sd = 2.5), rep(0.45, Tn), ref_pdie(),
                                  rep(3, Tn), rep(2, Tn), rep(0.35, Tn), history = hist,
                                  tau = 1, eta = 0.05, rho = 0.3, delta = 0.3)
  list(Tn = Tn, hist = hist, R = R, I = I, rd = rd, sim = sim)
}

s3_chain <- function(dat, niter, chain, seed, time_limit) {
  snap <- dat$sim$cases; snap$Rt_index <- 1L
  last <- tapply(snap$cases_reported, snap$time_onset, function(x) x[length(x)])
  init_cases <- tibble::tibble(time_onset = 1:dat$Tn,
                               cases_true = as.integer(ceiling(c(1.3, 1.1, 1.5, 1.2)[chain] * last) + 5))
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
                                  initial_delay = list(mean = 3, sd = 2.5),
                                  delay_prior = list(mean_mu = 4, mean_sigma = 4, sd_mu = 4, sd_sigma = 4),
                                  delay_changepoints = FALSE,
                                  dated_by_death = list(fraction = list(alpha = 1, beta = 1,
                                                                        initial = c(0.1, 0.25, 0.5, 0.7)[chain]),
                                                        exact = TRUE, ridge = TRUE),
                                  triggered = list(alpha = 1, beta = 1, initial = 1, step = 0.1, n_steps = 0),
                                  wrong_dead_alpha = 1, wrong_dead_beta = 1)))
}

# the recovery test of one dataset
s3_test <- function(replicate) {
  label <- sprintf("S3-%d", replicate)
  dat <- s3_data(replicate)
  note("%s data: %d reports (days %s); %d cases, %d fatal, %d misdated", label,
       length(dat$rd), paste(dat$rd, collapse = " "), sum(dat$I), sum(dat$sim$fatal),
       sum(dat$sim$misdated))
  niter <- 5000; n_chains <- 4
  t0 <- proc.time()[["elapsed"]]
  # in parallel where there are the cores; each chain stops at the time limit (minutes)
  fits <- parallel::mclapply(seq_len(n_chains), function(chain)
    s3_chain(dat, niter, chain, s3_seeds()[replicate] + chain, time_limit = 300),
    mc.cores = min(n_chains, parallel::detectCores()))
  note("%s: %d chains took %.1f min", label, n_chains, (proc.time()[["elapsed"]] - t0) / 60)
  for (fit in fits) if (inherits(fit, "try-error")) stop(fit)
  done <- min(vapply(fits, function(fit) max(fit$deaths$iteration), 0))
  note("%s: %d iterations in each chain", label, done)
  expect(done >= 500, sprintf("%s: only %d iterations done in each chain", label, done))

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
    "fraction misdated" = f_series("deaths_parameters", NULL, "dated_fraction"),
    "mean delay from death to recording" = f_series("deaths", 1, "delay_mean"),
    "true cases of the last 14 onset days" = lapply(fits, f_last_cases))
  shape <- Map(function(m, s) (m / s)^2, series[["mean delay from death to recording"]],
               f_series("deaths", 1, "delay_sd"))
  for (name in names(series))
    note("%s chain means, %s: %s", label, name,
         paste(sprintf("%.4g", vapply(series[[name]], mean, 0)), collapse = " "))
  note("%s chain means, shape of the recording delay (truth 2.25): %s", label,
       paste(sprintf("%.3g", vapply(shape, mean, 0)), collapse = " "))
  rhat <- vapply(series, function(x) split_rhat(unlist(x), rep(seq_along(x), lengths(x))), 0)
  for (name in names(rhat)) note("%s split R-hat, %s: %.3f", label, name, rhat[[name]])

  ci_check <- function(x, truth, what) {
    x <- unlist(x)
    ci <- stats::quantile(x, c(0.005, 0.995), names = FALSE)
    note("%s %s: truth %.4g, posterior mean %.4g, 99%% interval [%.4g, %.4g]", label, what, truth,
         mean(x), ci[1], ci[2])
    expect(truth >= ci[1] && truth <= ci[2],
           sprintf("%s %s: truth %.4g outside the 99%% interval [%.4g, %.4g] (posterior mean %.4g)",
                   label, what, truth, ci[1], ci[2], mean(x)))
  }
  ci_check(series[["CFR"]], 0.45, "CFR")
  ci_check(series[["fraction misdated"]], 0.35, "fraction misdated")
  ci_check(series[["mean delay from death to recording"]], 3, "mean delay from death to recording")
  ci_check(series[["true cases of the last 14 onset days"]], sum(dat$I[(dat$Tn - 13):dat$Tn]),
           "true cases of the last 14 onset days")
}
