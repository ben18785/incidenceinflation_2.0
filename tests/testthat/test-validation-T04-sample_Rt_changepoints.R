skip_validation("fast")

# T4. sample_Rt_changepoints: the partition of R_t and each segment's R are drawn
# exactly from their conditional distribution (spec A.10.5 with eq. 28), given p.
#   (a) Poisson, closed-form segment marginal likelihood;
#   (b) Poisson on the grid (spec: 400 values on [0.02, 20]);
#   (c) negative binomial, kappa = 5, segment marginal by quadrature.
# Design: 34 partitions, largest posterior probability 0.39 (Poisson),
# 0.36 (NB). A 400-point grid on [0.02, 20], linear or log spaced, changes every
# partition probability above 1e-3 by less than 0.02% relative, so grid error is
# negligible beside Monte Carlo error.

t4_serial <- list(mean = 15.3, sd = 9.3)
t4_omega <- ref_serial_weights(15.3, 9.3, 40)
t4_hist <- rep(20, 60)
t4_cases <- c(19, 24, 18, 22, 21, 31, 27, 35, 29, 33)
t4_L <- ref_lambda(t4_cases, t4_hist, t4_omega)
t4_prior <- list(shape = 1, rate = 1)
t4_p <- 0.3; t4_l <- 2; t4_n <- 10

t4_run <- function(N, ...) {
  df <- tibble::tibble(time_onset = 1:t4_n, cases_true = t4_cases)
  lapply(seq_len(N), function(i) sample_Rt_changepoints(
    df, t4_prior, t4_serial, p_change = t4_p, min_length = t4_l, serial_max = 40,
    cases_before = t4_hist, ...))
}

t4_check <- function(draws, log_seg, seg_mean, label) {
  keys <- vapply(draws, function(d) {
    idx <- col_of(d, "Rt_index"); tt <- col_of(d, "time_onset")
    o <- order(tt); idx <- idx[o]
    ref_partition_key(c(1, which(diff(idx) != 0) + 1))
  }, "")
  R <- do.call(rbind, lapply(draws, function(d) col_of(d, "Rt")[order(col_of(d, "time_onset"))]))
  parts <- ref_enumerate_partitions(t4_n, t4_l)
  pr <- ref_partition_posterior(parts, t4_n, log_seg,
                                function(st) ref_log_prior_partition_p(st, t4_n, t4_p))
  names(pr) <- vapply(parts, ref_partition_key, "")
  expect_gof(table(keys), pr, sprintf("%s partition frequencies", label))
  ex <- ref_expected_by_day(parts, pr, t4_n, seg_mean)
  for (t in seq_len(t4_n)) expect_mean_z(R[, t], ex[t], sprintf("%s E[R_%d]", label, t))
}

test_that("T4(a) sample_Rt_changepoints, Poisson closed form", {
  set.seed(4001)
  t0 <- proc.time()[["elapsed"]]
  draws <- t4_run(10000, is_negative_binomial = FALSE, use_grid = FALSE)
  note("T4(a): 10000 draws took %.1f s", proc.time()[["elapsed"]] - t0)
  t4_check(draws, ref_rt_log_seg_pois(t4_cases, t4_L, t4_prior$shape, t4_prior$rate),
           function(u, v) (t4_prior$shape + sum(t4_cases[u:v])) / (t4_prior$rate + sum(t4_L[u:v])),
           "T4(a)")
})

test_that("T4(b) sample_Rt_changepoints, Poisson on the grid", {
  set.seed(4002)
  t0 <- proc.time()[["elapsed"]]
  draws <- t4_run(10000, is_negative_binomial = FALSE, use_grid = TRUE)
  note("T4(b): 10000 draws took %.1f s", proc.time()[["elapsed"]] - t0)
  t4_check(draws, ref_rt_log_seg_pois(t4_cases, t4_L, t4_prior$shape, t4_prior$rate),
           function(u, v) (t4_prior$shape + sum(t4_cases[u:v])) / (t4_prior$rate + sum(t4_L[u:v])),
           "T4(b)")
})

test_that("T4(c) sample_Rt_changepoints, negative binomial (kappa = 5)", {
  set.seed(4003)
  t0 <- proc.time()[["elapsed"]]
  draws <- t4_run(10000, is_negative_binomial = TRUE, kappa = 5)
  note("T4(c): 10000 draws took %.1f s", proc.time()[["elapsed"]] - t0)
  t4_check(draws, ref_rt_log_seg_nb(t4_cases, t4_L, t4_prior$shape, t4_prior$rate, 5),
           function(u, v) ref_rt_seg_mean_nb(t4_cases, t4_L, t4_prior$shape, t4_prior$rate, 5, u, v),
           "T4(c)")
})
