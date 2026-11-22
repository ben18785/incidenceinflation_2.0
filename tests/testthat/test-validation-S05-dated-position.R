skip_validation("slow")

# S5 (slow tier): recovery with misdated cases recorded for a day part of the way from
# onset to death (S5 of the sensitivity analyses, test S5.8 and C9 of its test plan), on
# one dataset from the independent simulator (helper-validation-sim.R,
# sim_reports_deaths_dated with p_recorded and p_after).
#
# Data: T = 80 onset days; Poisson renewal with R = 1.2 on days 1-30 and 0.9 on days
# 31-80, history of 30 a day; case delay gamma(5, 2.5); constant CFR 0.45; time from onset
# to death p_die (gamma 9.3, 5.5 on (z-1, z], z = 1..36) for the correctly dated fatal
# cases; recording delay from death gamma(3, 2); triggered reports with tau = 1 (held
# there); wrong records eta = 0.05, rho = 0.3, shown dead with probability delta = 0.3; no
# backlog. A fatal case is misdated with probability 0.35 (those of the history too), and
# is then recorded for the day z1 days after onset and dies z2 days after that day:
# gamma(a / 2, c) on (z1 - 1, z1] and on (z2, z2 + 1], a and c the shape and rate of the
# time to death (ref_position_delays). Reports with gaps drawn from 1-7 days from day 2 to
# day 80 (about 25).
# Fit: four chains of mcmc_single with misdating in the exact form, an estimated fraction
# and dated_by_death position = TRUE (the default a1 = a / 2), the recording delay
# constant; at most 4 wrong records present (as in S9: the exact form is slow with many
# more); three changes of the labels proposed for each onset time in each iteration (the
# labels mix slowly). Each chain stops at 4.5 hours. The chains start near the truth (the
# true cases, floored at the most reported, and a fraction of 0.35), as fits of S5 start
# from the main fit: from overdispersed starting values every chain fell into a region
# where the gamma delay from death to recording collapses (shape near 0, about half the
# fatal cases misdated; joint log density 1,000-3,000 below the truth) and did not leave it
# in 4.5 hours, though a chain started at the truth stays there.
# Checks: max split R-hat < 1.2 over the CFR, the fraction misdated, the mean delay from
# death to recording and the true cases of the last 14 onset days in total; the central
# 99% interval (second half of each chain) contains the truth for each of them. A single
# dataset, so this is not a calibration test; 99% intervals keep false alarms low.

s5_data <- function() {
  set.seed(8501)
  Tn <- 80; hist <- rep(30, 40); omega <- ref_serial_weights(15.3, 9.3, 40)
  R <- c(rep(1.2, 30), rep(0.9, 50))
  I <- sim_renewal(R, hist, omega)
  rd <- sim_report_days(2, Tn, 1:7)
  position <- ref_position_delays(9.3, 5.5)
  sim <- sim_reports_deaths_dated(I, rd, list(mean = 5, sd = 2.5), rep(0.45, Tn), ref_pdie(),
                                  rep(3, Tn), rep(2, Tn), rep(0.35, Tn), history = hist,
                                  tau = 1, eta = 0.05, rho = 0.3, delta = 0.3,
                                  p_recorded = position$p_recorded,
                                  p_after = position$p_after)
  list(Tn = Tn, hist = hist, R = R, I = I, rd = rd, sim = sim)
}

s5_chain <- function(dat, niter, seed, time_limit) {
  snap <- dat$sim$cases; snap$Rt_index <- 1L
  most <- tapply(snap$cases_reported, snap$time_onset, max)
  init_cases <- tibble::tibble(time_onset = 1:dat$Tn,
                               cases_true = as.integer(pmax(dat$I, most[as.character(1:dat$Tn)],
                                                            na.rm = TRUE)))
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
                                  dated_by_death = list(fraction = list(initial = 0.35), exact = TRUE,
                                                        position = TRUE, n_moves = 3 * dat$Tn),
                                  triggered = list(alpha = 1, beta = 1, initial = 1, step = 0.1, n_steps = 0),
                                  wrong_dead_alpha = 1, wrong_dead_beta = 1)))
}

test_that("S5 misdated cases recorded part of the way to death: CFR, fraction and delay recovered (99% intervals)", {
  skip_validation("slow")
  dat <- s5_data()
  note("S5 data: %d reports (days %s); %d cases, %d fatal, %d misdated, %d recorded for a day analysed",
       length(dat$rd), paste(dat$rd, collapse = " "), sum(dat$I), sum(dat$sim$fatal),
       sum(dat$sim$misdated), sum(dat$sim$recorded_misdated))
  niter <- 5000; n_chains <- 4
  t0 <- proc.time()[["elapsed"]]
  # in parallel where there are the cores; each chain stops at the time limit (minutes)
  fits <- parallel::mclapply(seq_len(n_chains), function(chain)
    s5_chain(dat, niter, 8502 + chain, time_limit = 270),
    mc.cores = min(n_chains, parallel::detectCores()))
  note("S5: %d chains took %.1f min", n_chains, (proc.time()[["elapsed"]] - t0) / 60)
  for (fit in fits) if (inherits(fit, "try-error")) stop(fit)
  done <- min(vapply(fits, function(fit) max(fit$deaths$iteration), 0))
  note("S5: %d iterations in each chain", done)
  expect(done >= 500, sprintf("S5: only %d iterations done in each chain", done))

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
  rhat <- vapply(series, function(x) split_rhat(unlist(x), rep(seq_along(x), lengths(x))), 0)
  for (name in names(rhat)) note("S5 split R-hat, %s: %.3f", name, rhat[[name]])
  rhat[is.na(rhat)] <- Inf   # too few draws to tell
  expect(max(rhat) < 1.2, sprintf("S5: max split R-hat %.3f (%s) is not below 1.2",
                                  max(rhat), names(rhat)[which.max(rhat)]))

  ci_check <- function(x, truth, label) {
    x <- unlist(x)
    ci <- stats::quantile(x, c(0.005, 0.995), names = FALSE)
    note("S5 %s: truth %.4g, posterior mean %.4g, 99%% interval [%.4g, %.4g]", label, truth, mean(x), ci[1], ci[2])
    expect(truth >= ci[1] && truth <= ci[2],
           sprintf("S5 %s: truth %.4g outside the 99%% interval [%.4g, %.4g] (posterior mean %.4g)",
                   label, truth, ci[1], ci[2], mean(x)))
  }
  ci_check(series[["CFR"]], 0.45, "CFR")
  ci_check(series[["fraction misdated"]], 0.35, "fraction misdated")
  ci_check(series[["mean delay from death to recording"]], 3, "mean delay from death to recording")
  ci_check(series[["true cases of the last 14 onset days"]], sum(dat$I[(dat$Tn - 13):dat$Tn]),
           "true cases of the last 14 onset days")
})
