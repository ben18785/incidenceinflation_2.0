skip_validation("sbc")

# S1 (VALIDATION_TIER "all"): simulation-based calibration of the
# cases-only discrete model, fit_incidence(model = "discrete", method = "gibbs").
#
# 30 replicates. In each, every parameter is drawn from the priors passed to the fit
# (deliberately narrow), and the data from the independent simulator (helper-validation-sim.R):
# - partition of R_t: product-partition prior, p ~ Beta(1, T) integrated out, segments
#   >= 7 days (model_priors change_alpha = 1, change_beta = NA -> T; model_options
#   Rt_min_length = 7); R in each segment ~ Gamma(shape 40, rate 40);
# - Poisson renewal, history of 30 a day (cases_before), serial interval (15.3, 9.3);
# - delay mean ~ Gamma(mean 6, sd 0.5), delay sd ~ Gamma(mean 3, sd 0.4) (A3: mean/sd);
# - wrong records: eta ~ Gamma(mean 0.03, sd 0.008), rho ~ Beta(15, 35), entered and
#   withdrawn day by day (eqs. 23 and 24, A.7 first report);
#   max_erroneous = 8, so the cap is practically never reached;
# - reports with gaps drawn from 1-7 days, the last on day T = 40.
# Statistics, on the replicates that converge (split R-hat <
# 1.05 for every quantity; the others are excluded and reported):
# - SBC ranks of R_10, R_35, the delay mean, eta and the true cases of the last onset
#   day (I_40, the nowcast; added when the test was written, before any run), among
#   120 posterior draws (2 chains x 60): uniformity by chi-square on 5 bins, p >= 1e-3;
# - coverage of central 90% intervals: fails if the binomial lower tail p < 0.01.
# Package errors in a replicate are counted (and must be zero) with the message only.
# Self-checks of the simulator (plain R): it matches the forward filter of
# helper-validation-hmm.R for wrong records, and the prior sampler matches P(K) under
# the partition prior.

s1_settings <- list(T = 40, hist = rep(30, 40), R_shape = 40, R_rate = 40,
                    dm = c(6, 0.5), ds = c(3, 0.4), eta = c(0.03, 0.008), rho = c(15, 35))

s1_draw_truth <- function() {
  s <- s1_settings
  st <- sim_partition_prior(s$T, 1, s$T, 7)
  sg <- ref_segments(st, s$T)
  R <- numeric(s$T)
  for (r in seq_len(nrow(sg))) R[sg$u[r]:sg$v[r]] <- stats::rgamma(1, s$R_shape, s$R_rate)
  list(R = R,
       delay_mean = sim_rgamma_ms(1, s$dm[1], s$dm[2]), delay_sd = sim_rgamma_ms(1, s$ds[1], s$ds[2]),
       eta = sim_rgamma_ms(1, s$eta[1], s$eta[2]), rho = stats::rbeta(1, s$rho[1], s$rho[2]))
}

s1_rank <- function(draws, truth) {
  lo <- sum(draws < truth); ties <- sum(draws == truth)
  lo + if (ties > 0) sample(0:ties, 1) else 0
}

test_that("S1 simulation-based calibration of fit_incidence (discrete, cases only)", {
  skip_validation("sbc")
  s <- s1_settings; n_rep <- 30; chains <- 2; draws <- 60; L <- chains * draws
  omega <- ref_serial_weights(15.3, 9.3, 40)
  qs <- c("R_10", "R_35", "delay_mean", "eta", "I_40")
  ranks <- matrix(NA_real_, n_rep, length(qs), dimnames = list(NULL, qs))
  covered <- ranks; rhat <- ranks
  n_err <- 0
  pri <- model_priors(delay_mean_mu = s$dm[1], delay_mean_sigma = s$dm[2],
                      delay_sd_mu = s$ds[1], delay_sd_sigma = s$ds[2],
                      eta_mu = s$eta[1], eta_sigma = s$eta[2],
                      rho_alpha = s$rho[1], rho_beta = s$rho[2],
                      Rt_shape = s$R_shape, Rt_rate = s$R_rate, change_alpha = 1)
  opt <- model_options(is_negative_binomial = FALSE, Rt_model = "changepoints", Rt_min_length = 7,
                       max_erroneous = 8)
  t_all <- proc.time()[["elapsed"]]
  for (i in seq_len(n_rep)) {
    set.seed(9100 + i)
    tr <- s1_draw_truth()
    I <- sim_renewal(tr$R, s$hist, omega)
    rd <- sim_report_days(1, s$T, 1:7)
    snap <- sim_reports_wrong(I, rd, list(mean = tr$delay_mean, sd = tr$delay_sd), tr$eta, tr$rho)
    t0 <- proc.time()[["elapsed"]]
    fit <- tryCatch(fit_incidence(snap, list(mean = 15.3, sd = 9.3), model = "discrete",
                                  options = opt, priors = pri, cases_before = s$hist,
                                  control = list(niterations = 600, chains = chains, draws = draws,
                                                 cores = 2, seed = 9200 + i)),
                    error = function(e) { note("S1 replicate %d: the package raised an error: %s", i, conditionMessage(e)); NULL })
    if (is.null(fit)) { n_err <- n_err + 1; next }
    get <- function(el, col, t = NULL) {
      x <- fit[[el]]
      if (!is.null(t)) x <- x[col_of(x, "time_onset") == t, , drop = FALSE]
      x <- x[order(col_of(x, "chain"), col_of(x, "draw")), , drop = FALSE]
      list(v = col_of(x, col), ch = col_of(x, "chain"))
    }
    vals <- list(R_10 = get("Rt", "Rt", 10), R_35 = get("Rt", "Rt", 35),
                 delay_mean = get("delay", "mean", 1), eta = get("parameters", "eta"),
                 I_40 = get("cases", "cases_true", s$T))
    truth <- c(R_10 = tr$R[10], R_35 = tr$R[35], delay_mean = tr$delay_mean, eta = tr$eta, I_40 = I[s$T])
    for (q in qs) {
      v <- vals[[q]]$v
      ranks[i, q] <- s1_rank(v, truth[[q]])
      ci <- stats::quantile(v, c(0.05, 0.95), names = FALSE)
      covered[i, q] <- truth[[q]] >= ci[1] && truth[[q]] <= ci[2]
      rhat[i, q] <- split_rhat(v, vals[[q]]$ch)
    }
    note("S1 replicate %d: %.0f s; K = %d; ranks %s; R-hat %s", i, proc.time()[["elapsed"]] - t0,
         sum(diff(tr$R) != 0), paste(ranks[i, ], collapse = " "), paste(sprintf("%.3f", rhat[i, ]), collapse = " "))
  }
  note("S1: %d replicates took %.1f min; %d raised errors", n_rep, (proc.time()[["elapsed"]] - t_all) / 60, n_err)
  expect(n_err == 0, sprintf("S1: the package raised an error in %d of %d replicates", n_err, n_rep))
  ok <- stats::complete.cases(ranks) & apply(rhat, 1, function(r) all(is.finite(r) & r < 1.05))
  note("S1: %d of %d replicates converged (split R-hat < 1.05 for all quantities); excluded: %s",
       sum(ok), n_rep, paste(which(!ok & stats::complete.cases(ranks)), collapse = ", "))
  for (q in qs) {
    r <- ranks[ok, q]
    bins <- tabulate(floor(r / (L + 1) * 5) + 1, 5)
    p <- if (length(r) >= 10) stats::pchisq(sum((bins - length(r) / 5)^2 / (length(r) / 5)), 4, lower.tail = FALSE) else NA
    note("S1 %s: rank bins %s; chi-square p = %.3g; 90%% coverage %d of %d", q, paste(bins, collapse = " "), p,
         sum(covered[ok, q]), sum(ok))
    expect(is.finite(p) && p >= 1e-3, sprintf("S1 %s: SBC ranks not uniform (bins %s, chi-square p = %.3g < 0.001, n = %d)",
                                              q, paste(bins, collapse = " "), p, length(r)))
    expect_coverage(sum(covered[ok, q]), sum(ok), 0.9, sprintf("S1 %s 90%% intervals", q))
  }
})
