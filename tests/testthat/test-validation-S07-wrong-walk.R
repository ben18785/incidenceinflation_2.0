skip_validation("slow")

# S7 (slow tier): C9 of the sensitivity analyses for S7, a rate of wrong records which
# varies over the intervals between reports (sensitivity_extensions.tex, eq.
# s7-entered): recovery by the deaths model of mcmc_single, with
# priors$revisions$wrong_rate = "walk", of data from the independent simulator
# (helper-validation-sim.R).
#
# Data: as for S6 (test-validation-S06-wrong-constant.R), but with wrong records entered
# in proportion to the genuine cases being reported, at the rate eta_J(k) of the interval
# between reports containing day k: 0.06 for the intervals ending on or before day 40 and
# 0.12, double, after; withdrawn with probability rho = 0.3 on each later day and shown
# dead with probability delta = 0.3.
# Fit: four chains of 1500 iterations, run in parallel, with the default priors of the
# walk; the second half of each kept.
# Checks: max split R-hat < 1.2 over the quantities checked; the central 99% interval of
# the pooled draws contains the truth for the CFR, rho, the mean rate of the intervals
# ending on days 12 to 38 (before the change) and of those ending on days 50 to 78
# (after it), and the true cases of the last 14 onset days in total.

s7_data <- function() {
  set.seed(8701)
  Tn <- 80; hist <- rep(12, 40); omega <- ref_serial_weights(15.3, 9.3, 40)
  R <- c(rep(1.15, 30), rep(0.9, 50))
  I <- sim_renewal(R, hist, omega)
  rd <- sim_report_days(2, Tn, 1:7)
  rates <- ifelse(rd <= 40, 0.06, 0.12)
  rec_mean <- c(rep(3, 40), rep(6, 40)); rec_sd <- c(rep(2, 40), rep(3, 40))
  sim <- sim_reports_deaths(I, rd, list(mean = 5, sd = 2.5), rep(0.45, Tn), ref_pdie(),
                            rec_mean, rec_sd, tau = 1, rho = 0.3, delta = 0.3,
                            eta_day = function(k) rates[sim_interval_of_day(k, rd)])
  list(Tn = Tn, hist = hist, I = I, rd = rd, rates = rates, fatal = sim$fatal,
       cases = sim$cases, deaths = sim$deaths)
}

test_that("S7 (C9) a rate of wrong records which doubles is recovered with the CFR (99% intervals, four chains)", {
  skip_validation("slow")
  dat <- s7_data()
  note("S7 data: %d reports (days %s); %d cases, %d fatal; largest fall %d", length(dat$rd),
       paste(dat$rd, collapse = " "), sum(dat$I), sum(dat$fatal), max_downward_revision(dat$cases))
  niter <- 1500
  t0 <- proc.time()[["elapsed"]]
  fits <- s_wrong_chains(function(chain) s_wrong_fit(
    dat, niter, 8710 + chain,
    revisions = list(wrong_rate = "walk", rho_alpha = 1, rho_beta = 1),
    initial_revisions = list(eta = 0.07, rho = 0.3)))
  note("S7: four chains of %d iterations took %.1f min", niter, (proc.time()[["elapsed"]] - t0) / 60)
  expect_true(all(vapply(fits, function(fit) is.list(fit) && "revisions" %in% names(fit), logical(1))))
  f_deaths <- function(name) function(fit) {
    d <- fit$deaths[col_of(fit$deaths, "time_onset") == 1, ]
    col_of(d, name)[order(col_of(d, "iteration"))]
  }
  f_rate <- function(days) function(fit) {
    days_walk <- dat$rd[dat$rd >= min(days) & dat$rd <= max(days)]
    rowMeans(as.matrix(fit$revisions[, paste0("eta_", days_walk), drop = FALSE]))
  }
  s_wrong_check(s_wrong_kept(fits, f_deaths("cfr"), niter), 0.45, "S7 CFR")
  s_wrong_check(s_wrong_kept(fits, function(fit) fit$revisions$rho, niter), 0.3, "S7 rho")
  s_wrong_check(s_wrong_kept(fits, f_rate(c(12, 38)), niter), 0.06,
                "S7 mean rate of the intervals ending on days 12 to 38")
  s_wrong_check(s_wrong_kept(fits, f_rate(c(50, 78)), niter), 0.12,
                "S7 mean rate of the intervals ending on days 50 to 78")
  s_wrong_check(s_wrong_kept(fits, function(fit) s_wrong_last14(fit, dat$Tn), niter),
                sum(dat$I[(dat$Tn - 13):dat$Tn]), "S7 true cases of the last 14 onset days")
})
