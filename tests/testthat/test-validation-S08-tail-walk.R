skip_validation("slow")

# S8 (slow tier), C9 and S8.6 of sensitivity_extensions.tex: the deaths model of
# mcmc_single with a share of the tail of the reporting delay for each onset week, on one
# dataset from the independent simulator (helper-validation-sim.R).
#
# Data: as S2 (test-validation-S02-deaths.R), T = 80 onset days, Poisson renewal with
# R = 1.2 on days 1-30 and 0.9 on days 31-80, history of 30 a day; CFR 0.45; p_die;
# recording delay gamma(3, 2) for onset days 1-40 and gamma(6, 3) for 41-80; triggered
# reports with tau = 1; wrong records eta = 0.05, rho = 0.3, delta = 0.3. The case delay
# is gamma(5, 2.5) with probability 1 - q_t and exponential with mean 30 otherwise, where
# q_t is the share of the onset week of day t (weeks of 7 days from day 1), rising
# linearly from 0.02 in week 1 to 0.10 in week 12. Reports with gaps drawn from 1-7 days
# from day 2 to day 80.
# Fit: four chains of 2000 iterations, run in parallel, each stopped at 5 h (GitHub's limit
# for a job is 6 h), the second halves taken to the iterations all have done; the tail starting from a single share of 0.036 and
# a mean of 30; prior of the walk as in the note (logit q_(1) ~ N(logit 0.036, 1.2^2),
# tau_q ~ half-normal(0.15)); gamma prior of the mean of the tail with mean 30 and sd 20.
# Checks, on the second half of each chain: the central 99% interval contains the truth
# for the CFR, the shares of weeks 2 and 6, eta and the true cases of the last 14 onset
# days in total; the mean of the tail (informed by the walk and the other parameters
# together) is only noted. The split R-hat of each is noted, not checked.

s08_shares <- function() seq(0.02, 0.10, length.out = 12)

s08_data <- function() {
  set.seed(8801)
  Tn <- 80; hist <- rep(30, 40); omega <- ref_serial_weights(15.3, 9.3, 40)
  R <- c(rep(1.2, 30), rep(0.9, 50))
  I <- sim_renewal(R, hist, omega)
  rd <- sim_report_days(2, Tn, 1:7)
  rec_mean <- c(rep(3, 40), rep(6, 40)); rec_sd <- c(rep(2, 40), rep(3, 40))
  week <- (seq_len(Tn) - 1) %/% 7 + 1
  sim <- sim_reports_deaths(I, rd, list(mean = 5, sd = 2.5), rep(0.45, Tn), ref_pdie(),
                            rec_mean, rec_sd, tau = 1, eta = 0.05, rho = 0.3, delta = 0.3,
                            tail_share = s08_shares()[week], tail_mean = 30)
  list(Tn = Tn, hist = hist, R = R, I = I, rd = rd, fatal = sim$fatal,
       cases = sim$cases, deaths = sim$deaths)
}

s08_fit <- function(dat, niter, chain, time_limit = NULL) {
  snap <- dat$cases; snap$Rt_index <- 1L
  last <- tapply(snap$cases_reported, snap$time_onset, function(x) x[length(x)])
  init_cases <- tibble::tibble(time_onset = 1:dat$Tn,
                               cases_true = as.integer(ceiling(c(1.3, 1.1, 1.5, 1.2)[chain] * last) + 5))
  set.seed(8802 + chain)
  mcmc_single(
    niter, snap, time_limit = time_limit,
    priors = list(Rt = list(shape = 1, rate = 1),
                  reporting = list(mean_mu = 5, mean_sigma = 3, sd_mu = 3, sd_sigma = 2),
                  max_cases = 400, Rt_changepoints = list(), reporting_changepoints = list(),
                  revisions = list(eta_mu = 0.05, eta_sigma = 0.05, rho_alpha = 1, rho_beta = 4),
                  reporting_effects = list(tail_mean_mu = 30, tail_mean_sigma = 20,
                                           tail_weight_walk = TRUE),
                  max_erroneous = 8),
    serial_parameters = list(mean = 15.3, sd = 9.3),
    initial_cases_true = init_cases,
    initial_reporting_parameters = list(mean = 5, sd = 2.5),
    initial_Rt = tibble::tibble(Rt_index = 1L, Rt = 1),
    print_to_screen = FALSE, cases_before = dat$hist,
    initial_revision_parameters = list(eta = 0.05, rho = 0.3,
                                       tail = list(weight = 0.036, mean = 30)),
    deaths = list(snapshot = dat$deaths,
                  settings = list(cfr_changepoints = FALSE, initial_cfr = 0.45,
                                  onset_to_death = list(mean = 9.3, sd = 5.5, max_lag = 36),
                                  initial_delay = list(mean = 4, sd = 2),
                                  delay_prior = list(mean_mu = 5, mean_sigma = 4, sd_mu = 3, sd_sigma = 2),
                                  delay_changepoints = list(min_length = 20, first_length = 20),
                                  triggered = list(alpha = 1, beta = 1, initial = 1, step = 0.1, n_steps = 0),
                                  wrong_dead_alpha = 1, wrong_dead_beta = 1)))
}

test_that("S8 (C9) the deaths model with a share of the tail for each onset week recovers the shares", {
  skip_validation("slow")
  dat <- s08_data()
  note("S8 data: %d reports (days %s); %d cases, %d fatal", length(dat$rd), paste(dat$rd, collapse = " "),
       sum(dat$I), sum(dat$fatal))
  niter <- 2000
  n_chains <- 4
  t0 <- proc.time()[["elapsed"]]
  # each chain stops at the time limit (minutes), inside the job's limit of six hours
  fits <- parallel::mclapply(seq_len(n_chains), function(chain) s08_fit(dat, niter, chain, 300),
                             mc.cores = min(n_chains, parallel::detectCores()))
  failed <- vapply(fits, function(f) inherits(f, "try-error"), logical(1))
  expect(!any(failed), paste("S8: a chain failed:", paste(unlist(fits[failed]), collapse = "; ")))
  if (any(failed)) return(invisible())
  # the iterations all the chains have done
  niter <- min(vapply(fits, function(fit) max(col_of(fit$deaths, "iteration")), 0))
  note("S8: %d chains of %d iterations took %.1f min", n_chains, niter,
       (proc.time()[["elapsed"]] - t0) / 60)
  expect(niter >= 500, sprintf("S8: only %d iterations done in each chain", niter))
  fits <- lapply(fits, function(fit) lapply(fit, function(x)
    if (is.data.frame(x) && "iteration" %in% names(x)) x[col_of(x, "iteration") <= niter, ] else x))

  # the second half of each chain of a quantity, and its chain
  f_kept <- function(f_draws) {
    values <- lapply(fits, f_draws)
    list(x = unlist(values), chain = rep(seq_along(values), lengths(values)))
  }
  shares <- s08_shares()
  last14 <- function(fit) {
    cs <- fit$cases[col_of(fit$cases, "iteration") > niter / 2 &
                      col_of(fit$cases, "time_onset") >= dat$Tn - 13, ]
    as.numeric(tapply(col_of(cs, "cases_true"), col_of(cs, "iteration"), sum))
  }
  quantities <- list(
    list(label = "CFR", truth = 0.45, f = function(fit) m_series(fit$deaths, 1, "cfr", niter)),
    list(label = "share of the tail, onset week 2", truth = shares[2],
         f = function(fit) m_kept(fit$reporting_effects, "tail_weight_week_2", niter)),
    list(label = "share of the tail, onset week 6", truth = shares[6],
         f = function(fit) m_kept(fit$reporting_effects, "tail_weight_week_6", niter)),
    list(label = "eta", truth = 0.05, f = function(fit) m_kept(fit$revisions, "eta", niter)),
    list(label = "true cases of the last 14 onset days", truth = sum(dat$I[(dat$Tn - 13):dat$Tn]),
         f = last14),
    list(label = "mean of the tail", truth = 30, interval = FALSE,
         f = function(fit) m_kept(fit$reporting_effects, "tail_mean", niter)))
  for (q in quantities) {
    kept <- f_kept(q$f)
    rhat <- split_rhat(kept$x, kept$chain)
    ci <- stats::quantile(kept$x, c(0.005, 0.995), names = FALSE)
    note("S8 %s: truth %.4g, posterior mean %.4g, 99%% interval [%.4g, %.4g], split R-hat %.3f",
         q$label, q$truth, mean(kept$x), ci[1], ci[2], rhat)
    if (!isFALSE(q$interval))
      expect(q$truth >= ci[1] && q$truth <= ci[2],
             sprintf("S8 %s: truth %.4g outside the 99%% interval [%.4g, %.4g] (posterior mean %.4g)",
                     q$label, q$truth, ci[1], ci[2], mean(kept$x)))
  }
})
