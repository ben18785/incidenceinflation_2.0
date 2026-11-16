skip_validation("medium")

# M1, M2 (medium tier): whole-sampler checks of the CFR with the deaths model
# of mcmc_single, in a regime where the CFR posterior is known exactly.
#
# Data from the independent simulator (helper-validation-sim.R): 40 onset days, reports every 2 days to day
# 40 and then on days 50, 70, 100, 140 and 190. Case delay gamma(4, 2); a fatal case's
# death is recorded after gamma(8, 3) from onset. At the last report (>= 150 days after
# every onset) every case and death is recorded: survival at 150 days is 0 to double
# precision for the true delays and 3e-6 even for a delay of mean 30, sd 15. So the true
# fatal cases are pinned, F(t) = D(t | t_J), and with fixed_cases = TRUE at the true I:
#   M1 (constant CFR):  CFR | data ~ Beta(alpha + sum F, beta + sum (I - F))   (A.10.4)
#   M2 (CFR changepoints): pi(P | F, I) propto B(a + K, b + T - 1 - K)
#        prod_seg B(alpha + F_uv, beta + I_uv - F_uv) on segments >= min_length
#        (eqs. 11 and 28), by dynamic programming (checked against full enumeration of
#        the 819 partitions: agreement to 1e-14).
# Checks: batch-means z <= 4 on the second half of the chain.
# Design: M1 exact Beta(340, 804), mean 0.2972, sd
# 0.0135. Power: swapping the prior's alpha and beta moves the mean by z ~ 6 at 500 kept
# near-independent draws. M2: P(K = 1) = 0.872, P(K = 2) = 0.123.

m_rdays <- c(seq(2, 40, 2), 50, 70, 100, 140, 190)
m_omega <- ref_serial_weights(15.3, 9.3, 40)
m_cfr_data <- function(seed, cfr) {
  set.seed(seed)
  I <- sim_renewal(rep(1, 40), rep(30, 40), m_omega)
  d <- sim_reports(I, m_rdays, list(mean = 4, sd = 2), cfr = cfr, death_delay = list(mean = 8, sd = 3))
  FJ <- d$deaths$deaths_reported[d$deaths$time_reported == max(m_rdays)]
  stopifnot(all(FJ == d$fatal))   # pinned in the data themselves
  list(I = I, d = d, F = FJ)
}

m_cfr_run <- function(dat, niter, settings) {
  snap <- dat$d$cases; snap$Rt_index <- 1L
  mcmc_single(niter, snap,
              priors = list(Rt = list(shape = 1, rate = 1),
                            reporting = list(mean_mu = 5, mean_sigma = 3, sd_mu = 3, sd_sigma = 2),
                            max_cases = 300, reporting_changepoints = list()),
              serial_parameters = m_serial,
              initial_cases_true = tibble::tibble(time_onset = seq_along(dat$I), cases_true = as.integer(dat$I)),
              initial_reporting_parameters = list(mean = 4, sd = 2),
              initial_Rt = tibble::tibble(Rt_index = 1L, Rt = 1),
              print_to_screen = FALSE, cases_before = rep(30, 40), fixed_cases = TRUE,
              deaths = list(snapshot = dat$d$deaths, settings = settings))
}

m_check_pinned <- function(out, dat, niter, label) {
  dd <- out$deaths
  dd <- dd[col_of(dd, "iteration") > niter / 2, , drop = FALSE]
  bad <- sum(col_of(dd, "deaths_true") != dat$F[col_of(dd, "time_onset")])
  expect(bad == 0, sprintf("%s: deaths_true differs from the final reported deaths in %d of %d rows", label, bad, nrow(dd)))
}

test_that("M1 mcmc_single with deaths: constant CFR posterior is the exact beta", {
  skip_validation("medium")
  dat <- m_cfr_data(7001, rep(0.3, 40))
  a <- 3 + sum(dat$F); b <- 7 + sum(dat$I - dat$F)
  niter <- 1000
  set.seed(7101)
  t0 <- proc.time()[["elapsed"]]
  out <- m_cfr_run(dat, niter, list(cfr_changepoints = FALSE, cfr_alpha = 3, cfr_beta = 7,
                                     delay_changepoints = FALSE, initial_delay = list(mean = 8, sd = 3),
                                     triggered = FALSE))
  note("M1: %d iterations took %.1f s", niter, proc.time()[["elapsed"]] - t0)
  m_check_pinned(out, dat, niter, "M1")
  x <- m_series(out$deaths, 1, "cfr", niter)
  expect_equal(length(x), niter / 2)
  ex_mean <- a / (a + b); ex_var <- a * b / ((a + b)^2 * (a + b + 1))
  note("M1: exact Beta(%g, %g) mean %.5f sd %.5f; sample mean %.5f sd %.5f", a, b, ex_mean, sqrt(ex_var), mean(x), stats::sd(x))
  diag <- sprintf("(swapped prior would give mean %.5f)", (7 + sum(dat$F)) / (10 + sum(dat$I)))
  expect_mcmc_mean(x, ex_mean, "M1 posterior mean of the CFR", diag)
  expect_mcmc_mean((x - ex_mean)^2, ex_var, "M1 posterior variance of the CFR")
  for (q in c(0.05, 0.5, 0.95))
    expect_mcmc_mean(as.numeric(x <= stats::qbeta(q, a, b)), q, sprintf("M1 P(CFR <= exact %.2f quantile)", q))
})

test_that("M2 mcmc_single with deaths: CFR changepoints, exact partition posterior", {
  skip_validation("medium")
  dat <- m_cfr_data(7002, c(rep(0.45, 20), rep(0.25, 20)))
  al <- 2; be <- 5; ml <- 7
  dp <- ref_partition_dp(40, ref_cfr_log_seg(dat$F, dat$I, al, be),
                         function(u, v) (al + sum(dat$F[u:v])) / (al + be + sum(dat$I[u:v])),
                         a = 1, b = 40, min_length = ml)
  niter <- 1500
  set.seed(7102)
  t0 <- proc.time()[["elapsed"]]
  # change prior alpha and beta left at their defaults (1 and the number of onset days)
  out <- m_cfr_run(dat, niter, list(cfr_changepoints = TRUE, cfr_alpha = al, cfr_beta = be,
                                     min_length = ml, delay_changepoints = FALSE,
                                     initial_delay = list(mean = 8, sd = 3), triggered = FALSE))
  note("M2: %d iterations took %.1f s", niter, proc.time()[["elapsed"]] - t0)
  m_check_pinned(out, dat, niter, "M2")
  K <- m_kept(out$deaths_parameters, "n_cfr_changes", niter)
  note("M2: exact P(K = 0..3) %s; sample %s", paste(sprintf("%.4f", dp$pK[1:4]), collapse = " "),
       paste(sprintf("%.4f", vapply(0:3, function(k) mean(K == k), 0)), collapse = " "))
  for (k in 1:2) expect_mcmc_mean(as.numeric(K == k), dp$pK[k + 1], sprintf("M2 P(K = %d)", k))
  for (t in c(1, 18, 20, 22, 40))
    expect_mcmc_mean(m_series(out$deaths, t, "cfr", niter), dp$mean_by_day[t], sprintf("M2 E[CFR_%d]", t))
})
