skip_validation("fast")

# T5. sample_cases_history: posterior of the true cases of an onset day.
# Without wrong records the chain of binomials depends on N only through
# Binomial(A_J | N, 1 - S_J), with S_J = 1 - G(t_J - t + 1/2) (spec A.10.2, Steps 1-3;
# A.4). With the renewal prior this gives closed forms:
#   Poisson:  I = C_J + U,  U ~ Poisson(R Lambda (1 - p)),           p = G(t_J - t + 1/2)
#   NB(kappa): U ~ NegBin(size kappa + C_J, mean (kappa + C_J)(1 - p)/(kappa/mu + p))
# Power: using G(t_J - t) or G(t_J - t + 1) instead moves the mean of U
# by z ~ 25; omitting the later onset day's renewal term in (c) by z ~ 22.

t5_serial <- list(mean = 15.3, sd = 9.3)
t5_omega <- ref_serial_weights(15.3, 9.3, 40)
t5_rep <- tibble::tibble(reporting_piece_index = 1L, mean = 8, sd = 4)
t5_obs <- tibble::tibble(time_onset = 1L, time_reported = c(2L, 4L, 6L),
                         cases_reported = c(1L, 4L, 9L), reporting_piece_index = 1L)

t5_estimates <- function(out, onset) {
  x <- out[col_of(out, "time_onset") == onset, , drop = FALSE]
  v <- unique(col_of(x, "cases_estimated"))
  if (length(v) != 1) stop(sprintf("expected one cases_estimated for onset %d, got %d", onset, length(v)))
  v
}

t5_draws <- function(N, obs, onset, ...) {
  vapply(seq_len(N), function(i) t5_estimates(sample_cases_history(obs, ...), onset), numeric(1))
}

test_that("T5(a) Poisson renewal: I = C_J + Poisson(R Lambda (1 - p))", {
  set.seed(5001)
  R <- 1.2; hist <- rep(30, 40)
  mu <- R * ref_lambda(0, hist, t5_omega)[1]
  p <- ref_cdf_case(5, 8, 4)
  t0 <- proc.time()[["elapsed"]]
  I <- t5_draws(4000, t5_obs, 1, max_cases = 300, Rt_function = function(t) rep(R, length(t)),
                serial_parameters = t5_serial, reporting_parameters = t5_rep,
                is_negative_binomial = FALSE, cases_before = hist)
  note("T5(a): 4000 calls took %.1f s", proc.time()[["elapsed"]] - t0)
  U <- I - 9
  sup <- 0:291
  expect_pmf(U, sup, stats::dpois(sup, mu * (1 - p)), "T5(a) pmf of I - C_J",
             sprintf("(G(delta) would give mean %.2f, G(delta+1) %.2f)",
                     mu * (1 - ref_pgamma_ms(5, 8, 4)), mu * (1 - ref_pgamma_ms(6, 8, 4))))
  expect_mean_z(U, mu * (1 - p), "T5(a) mean of I - C_J")
})

test_that("T5(b) negative binomial renewal (kappa = 3): gamma-Poisson closed form", {
  set.seed(5002)
  R <- 1.2; hist <- rep(30, 40); kap <- 3; C <- 9
  mu <- R * ref_lambda(0, hist, t5_omega)[1]
  p <- ref_cdf_case(5, 8, 4)
  I <- t5_draws(4000, t5_obs, 1, max_cases = 300, Rt_function = function(t) rep(R, length(t)),
                serial_parameters = t5_serial, reporting_parameters = t5_rep,
                kappa = kap, is_negative_binomial = TRUE, cases_before = hist)
  U <- I - C
  m <- (kap + C) * (1 - p) / (kap / mu + p)
  sup <- 0:291
  expect_pmf(U, sup, stats::dnbinom(sup, size = kap + C, mu = m), "T5(b) pmf of I - C_J")
  expect_mean_z(U, m, "T5(b) mean of I - C_J")
})

test_that("T5(c) conditioning on the later onset day (cases_true_current)", {
  set.seed(5003)
  serial <- list(mean = 3, sd = 1.5)
  om <- ref_serial_weights(3, 1.5, 40)
  R <- 1.2; hist <- rep(30, 40); kap <- 10
  obs <- tibble::tibble(time_onset = c(1L, 2L), time_reported = c(3L, 3L),
                        cases_reported = c(3L, 2L), reporting_piece_index = 1L)
  I2 <- 70
  current <- tibble::tibble(time_onset = 1:2, cases_true = c(40, I2))
  # exact conditional of I_1 given I_2 = 70 and the reports of day 1
  Iv <- 3:400
  L1 <- ref_lambda(c(0, 0), hist, om)[1]
  L2 <- vapply(Iv, function(i1) ref_lambda(c(i1, 0), hist, om)[2], 0)
  lw <- ref_dnb(Iv, R * L1, kap, log = TRUE) +
    stats::dbinom(3, Iv, ref_cdf_case(2, 8, 4), log = TRUE) +
    ref_dnb(I2, R * L2, kap, log = TRUE)
  w <- exp(lw - max(lw)); w <- w / sum(w)
  lw0 <- lw - ref_dnb(I2, R * L2, kap, log = TRUE); w0 <- exp(lw0 - max(lw0)); w0 <- w0 / sum(w0)
  I1 <- t5_draws(4000, obs, 1, max_cases = 400, Rt_function = function(t) rep(R, length(t)),
                 serial_parameters = serial, reporting_parameters = t5_rep,
                 kappa = kap, is_negative_binomial = TRUE, cases_true_current = current,
                 cases_before = hist)
  expect_pmf(I1, Iv, w, "T5(c) pmf of I_1 given I_2",
             sprintf("(without the later day's renewal term the mean would be %.2f)", sum(Iv * w0)))
  expect_mean_z(I1, sum(Iv * w), "T5(c) mean of I_1 given I_2")
})

test_that("T5(d) onset days past p_gamma_cutoff return their latest count", {
  set.seed(5004)
  obs <- tibble::tibble(time_onset = 1L, time_reported = c(10L, 40L),
                        cases_reported = c(18L, 25L), reporting_piece_index = 1L)
  expect_gt(ref_cdf_case(39, 8, 4), 0.99)
  I <- t5_draws(20, obs, 1, max_cases = 300, Rt_function = function(t) rep(1.2, length(t)),
                serial_parameters = t5_serial, reporting_parameters = t5_rep,
                p_gamma_cutoff = 0.99, is_negative_binomial = FALSE, cases_before = rep(30, 40))
  expect_true(all(I == 25), label = "T5(d) every draw equals the latest count 25")
})

# T5(e) tests eq. (23): wrong records entered between reports are counted day by day and
# thinned to those still present at t_j,
#   n_j ~ Poisson(eta Itilde sum_{k=t_{j-1}+1}^{t_j} [S(k-1) - S(k)] (1-rho)^(t_j-k)),
# matching A.7's first report. Self-check (plain R): this filter matches simulated
# report probabilities of the daily process (|z| <= 2.3 over 8 patterns, 2e5 draws).
# Design (exact posterior mean 23.683): counting every entry since t_{j-1} as present
# at t_j gives 23.457 (z ~ 5); withdrawal probability rho per interval instead of
# 1-(1-rho)^gap moves the mean by z ~ 5.6; no wrong records at the first report z ~ 4.3;
# doubling eta z ~ 13; no half day z ~ 27; cap 6 vs 40 changes the mean by 2e-4.
test_that("T5(e) with wrong records: forward filter (eqs. 23, 24, A.7)", {
  set.seed(5005)
  R <- 1.0; hist <- rep(30, 40); eta <- 0.15; rho <- 0.3
  rd <- c(2L, 4L, 6L, 9L); cnt <- c(4L, 9L, 8L, 12L)    # one fall, at day 6
  obs <- tibble::tibble(time_onset = 1L, time_reported = rd, cases_reported = cnt,
                        reporting_piece_index = 1L)
  mu <- R * ref_lambda(0, hist, t5_omega)[1]
  surv <- function(d) ref_surv_case(d, 8, 4)
  Nv <- 0:200
  ll <- vapply(Nv, function(n) ref_hmm_loglik_cases(cnt, rd, 1, n, surv, eta, rho, e_max = 6,
                                                    entry_daily = TRUE), 0)
  lw <- stats::dpois(Nv, mu, log = TRUE) + ll
  w <- exp(lw - max(lw)); w <- w / sum(w)
  t0 <- proc.time()[["elapsed"]]
  I <- t5_draws(4000, obs, 1, max_cases = 200, Rt_function = function(t) rep(R, length(t)),
                serial_parameters = t5_serial, reporting_parameters = t5_rep,
                is_negative_binomial = FALSE,
                revision_parameters = list(eta = eta, rho = rho), max_erroneous = 6,
                cases_before = hist)
  note("T5(e): 4000 calls took %.1f s; exact posterior mean %.3f ",
       proc.time()[["elapsed"]] - t0, sum(Nv * w))
  expect_pmf(I, Nv, w, "T5(e) pmf of I with wrong records")
  expect_mean_z(I, sum(Nv * w), "T5(e) mean of I with wrong records")
})
