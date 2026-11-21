skip_validation("slow")

# S6 (slow tier): C9 of the sensitivity analyses for S6, wrong records at a constant
# rate (sensitivity_extensions.tex, eq. s6-entered): recovery by the deaths model of
# mcmc_single, with priors$revisions$wrong_rate = "constant", of data from the
# independent simulator (helper-validation-sim.R).
#
# Data: T = 80 onset days; Poisson renewal with R = 1.15 on days 1-30 and 0.9 on days
# 31-80, history of 12 a day; case delay gamma(5, 2.5); constant CFR 0.45; time from onset
# to death p_die (gamma 9.3, 5.5 on (z-1, z], z = 1..36); recording delay from death
# gamma(3, 2) for onset days 1-40 and gamma(6, 3) for 41-80; triggered reports with tau = 1
# (held there); wrong records entered at the constant rate eta_bg = 0.02 for each onset
# day on each later day, withdrawn with probability rho = 0.3 on each later day and shown
# dead with probability delta = 0.3. Reports with gaps drawn from 1-7 days from day 2 to
# day 80 (about 25).
# Fit: four chains of 1500 iterations, run in parallel; the second half of each kept.
# Checks: max split R-hat < 1.2 over the quantities checked; the central 99% interval of
# the pooled draws contains the truth for the CFR, eta_bg, rho and the true cases of the
# last 14 onset days in total.

s6_data <- function() {
  set.seed(8601)
  Tn <- 80; hist <- rep(12, 40); omega <- ref_serial_weights(15.3, 9.3, 40)
  R <- c(rep(1.15, 30), rep(0.9, 50))
  I <- sim_renewal(R, hist, omega)
  rd <- sim_report_days(2, Tn, 1:7)
  rec_mean <- c(rep(3, 40), rep(6, 40)); rec_sd <- c(rep(2, 40), rep(3, 40))
  sim <- sim_reports_deaths(I, rd, list(mean = 5, sd = 2.5), rep(0.45, Tn), ref_pdie(),
                            rec_mean, rec_sd, tau = 1, rho = 0.3, delta = 0.3, eta_bg = 0.02)
  list(Tn = Tn, hist = hist, I = I, rd = rd, fatal = sim$fatal, cases = sim$cases,
       deaths = sim$deaths)
}

test_that("S6 (C9) a constant rate of wrong records is recovered with the CFR (99% intervals, four chains)", {
  skip_validation("slow")
  dat <- s6_data()
  note("S6 data: %d reports (days %s); %d cases, %d fatal; largest fall %d", length(dat$rd),
       paste(dat$rd, collapse = " "), sum(dat$I), sum(dat$fatal), max_downward_revision(dat$cases))
  niter <- 1500
  t0 <- proc.time()[["elapsed"]]
  fits <- s_wrong_chains(function(chain) s_wrong_fit(
    dat, niter, 8610 + chain,
    revisions = list(wrong_rate = "constant", eta_bg_mu = 0.02, eta_bg_sigma = 0.02,
                     rho_alpha = 1, rho_beta = 1),
    initial_revisions = list(eta_bg = 0.02, rho = 0.3)))
  note("S6: four chains of %d iterations took %.1f min", niter, (proc.time()[["elapsed"]] - t0) / 60)
  expect_true(all(vapply(fits, function(fit) is.list(fit) && "revisions" %in% names(fit), logical(1))))
  for (fit in fits) note("S6: probability of more wrong records than allowed for, largest %.2g",
                         fit$wrong_cap$largest)
  f_deaths <- function(name) function(fit) {
    d <- fit$deaths[col_of(fit$deaths, "time_onset") == 1, ]
    col_of(d, name)[order(col_of(d, "iteration"))]
  }
  s_wrong_check(s_wrong_kept(fits, f_deaths("cfr"), niter), 0.45, "S6 CFR")
  s_wrong_check(s_wrong_kept(fits, function(fit) fit$revisions$eta_bg, niter), 0.02, "S6 eta_bg")
  s_wrong_check(s_wrong_kept(fits, function(fit) fit$revisions$rho, niter), 0.3, "S6 rho")
  s_wrong_check(s_wrong_kept(fits, function(fit) s_wrong_last14(fit, dat$Tn), niter),
                sum(dat$I[(dat$Tn - 13):dat$Tn]), "S6 true cases of the last 14 onset days")
})
