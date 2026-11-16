skip_validation("slow")

# S2 (slow tier): the deaths model of mcmc_single on one realistic dataset from the
# independent simulator (helper-validation-sim.R).
#
# S2 data: T = 80 onset days; Poisson renewal with R = 1.2 on days 1-30 and 0.9 on days
# 31-80, history of 30 a day; case delay gamma(5, 2.5); constant CFR 0.45; time from
# onset to death p_die (gamma 9.3, 5.5 on (z-1, z], z = 1..36); recording delay from
# death gamma(3, 2) for onset days 1-40 and gamma(6, 3) for 41-80; triggered reports with
# tau = 1 (held there); wrong records eta = 0.05, rho = 0.3, shown dead with probability
# delta = 0.3, entered and withdrawn day by day (eq. 23); no misdating, no
# backlog. Reports with gaps drawn from 1-7 days from day 2 to day 80 (about 25), plus a
# report on day 94 that is simulated but not given to the fit.
# S2 check: the central 99% interval (second half of the chain) contains the truth for
# the CFR, the recording delay mean on each side of its change (onset days 20 and 60),
# eta, R_t on days 10, 40 and 70, and the true cases of the last 14 onset days in total.
# A single dataset, so this is not a calibration test; 99% intervals keep false alarms
# low.
# Self-check of the simulator (plain R): it matches the closed-form E[C] and E[D] at every report (|z| < 1.4, 300 replicates).

s2_data <- function() {
  set.seed(8201)
  Tn <- 80; hist <- rep(30, 40); omega <- ref_serial_weights(15.3, 9.3, 40)
  R <- c(rep(1.2, 30), rep(0.9, 50))
  I <- sim_renewal(R, hist, omega)
  rd <- sim_report_days(2, Tn, 1:7)
  rec_mean <- c(rep(3, 40), rep(6, 40)); rec_sd <- c(rep(2, 40), rep(3, 40))
  sim <- sim_reports_deaths(I, c(rd, 94L), list(mean = 5, sd = 2.5), rep(0.45, Tn), ref_pdie(),
                            rec_mean, rec_sd, tau = 1, eta = 0.05, rho = 0.3, delta = 0.3)
  list(Tn = Tn, hist = hist, R = R, I = I, rd = rd, fatal = sim$fatal,
       cases = sim$cases[sim$cases$time_reported <= Tn, ], deaths = sim$deaths[sim$deaths$time_reported <= Tn, ])
}

test_that("S2 mcmc_single deaths model recovers the truth on one dataset (99% intervals)", {
  skip_validation("slow")
  dat <- s2_data()
  note("S2 data: %d reports (days %s); %d cases, %d fatal", length(dat$rd), paste(dat$rd, collapse = " "),
       sum(dat$I), sum(dat$fatal))
  snap <- dat$cases; snap$Rt_index <- 1L
  last <- tapply(snap$cases_reported, snap$time_onset, function(x) x[length(x)])
  init_cases <- tibble::tibble(time_onset = 1:dat$Tn, cases_true = as.integer(ceiling(1.3 * last) + 5))
  niter <- 3000
  set.seed(8202)
  t0 <- proc.time()[["elapsed"]]
  out <- mcmc_single(
    niter, snap,
    priors = list(Rt = list(shape = 1, rate = 1),
                  reporting = list(mean_mu = 5, mean_sigma = 3, sd_mu = 3, sd_sigma = 2),
                  max_cases = 400, Rt_changepoints = list(), reporting_changepoints = list(),
                  revisions = list(eta_mu = 0.05, eta_sigma = 0.05, rho_alpha = 1, rho_beta = 4),
                  max_erroneous = 8),
    serial_parameters = list(mean = 15.3, sd = 9.3),
    initial_cases_true = init_cases,
    initial_reporting_parameters = list(mean = 5, sd = 2.5),
    initial_Rt = tibble::tibble(Rt_index = 1L, Rt = 1),
    print_to_screen = FALSE, cases_before = dat$hist,
    initial_revision_parameters = list(eta = 0.05, rho = 0.3),
    return_state = TRUE,
    deaths = list(snapshot = dat$deaths,
                  settings = list(cfr_changepoints = FALSE, initial_cfr = 0.45,
                                  onset_to_death = list(mean = 9.3, sd = 5.5, max_lag = 36),
                                  initial_delay = list(mean = 4, sd = 2),
                                  delay_prior = list(mean_mu = 5, mean_sigma = 4, sd_mu = 3, sd_sigma = 2),
                                  delay_changepoints = list(min_length = 20, first_length = 20),
                                  triggered = list(alpha = 1, beta = 1, initial = 1, step = 0.1, n_steps = 0),
                                  wrong_dead_alpha = 1, wrong_dead_beta = 1)))
  note("S2: %d iterations took %.1f min", niter, (proc.time()[["elapsed"]] - t0) / 60)

  ci_check <- function(x, truth, label) {
    ci <- stats::quantile(x, c(0.005, 0.995), names = FALSE)
    note("S2 %s: truth %.4g, posterior mean %.4g, 99%% interval [%.4g, %.4g]", label, truth, mean(x), ci[1], ci[2])
    expect(truth >= ci[1] && truth <= ci[2],
           sprintf("S2 %s: truth %.4g outside the 99%% interval [%.4g, %.4g] (posterior mean %.4g)",
                   label, truth, ci[1], ci[2], mean(x)))
  }
  ci_check(m_series(out$deaths, 1, "cfr", niter), 0.45, "CFR")
  ci_check(m_series(out$deaths, 20, "delay_mean", niter), 3, "recording delay mean, onset day 20")
  ci_check(m_series(out$deaths, 60, "delay_mean", niter), 6, "recording delay mean, onset day 60")
  ci_check(m_kept(out$revisions, "eta", niter), 0.05, "eta")
  for (t in c(10, 40, 70)) ci_check(m_series(out$Rt, t, "Rt", niter), dat$R[t], sprintf("R_%d", t))
  cs <- out$cases[col_of(out$cases, "iteration") > niter / 2 & col_of(out$cases, "time_onset") >= dat$Tn - 13, ]
  tot <- tapply(col_of(cs, "cases_true"), col_of(cs, "iteration"), sum)
  ci_check(as.numeric(tot), sum(dat$I[(dat$Tn - 13):dat$Tn]), "true cases of the last 14 onset days")
})
