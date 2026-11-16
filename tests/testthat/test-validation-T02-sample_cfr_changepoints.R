skip_validation("fast")

# T2. sample_cfr_changepoints draws the partition of the CFR exactly from
#   pi(P | F, I, p) proportional to p^K (1-p)^(T-1-K) prod_seg B(a+F_uv, b+I_uv-F_uv)/B(a, b)
# and then each segment's CFR from Beta(a + F_uv, b + I_uv - F_uv)
# (spec A.10.4, Steps 1-3, eq. 28). Expected values by enumeration.
#
# Design: posterior spread over 13 partitions (largest 0.32);
# swapping alpha and beta moves partition probabilities by up to 0.37.

t2_I <- c(20, 25, 0, 30, 22, 28, 26, 24)      # day 3 has no cases (zero-count edge case)
t2_F <- c(5, 7, 0, 9, 11, 12, 12, 11)
t2_alpha <- 2; t2_beta <- 5; t2_p <- 0.3; t2_l <- 2

test_that("T2 sample_cfr_changepoints: exact partition and CFR draw (eq. 28)", {
  set.seed(2002)
  n_days <- length(t2_I); N <- 20000
  t0 <- proc.time()[["elapsed"]]
  draws <- lapply(seq_len(N), function(i) sample_cfr_changepoints(
    deaths_true = t2_F, cases_true = t2_I, p_change = t2_p, min_length = t2_l,
    cfr_alpha = t2_alpha, cfr_beta = t2_beta, is_changepoints = TRUE))
  note("T2: %d draws took %.1f s", N, proc.time()[["elapsed"]] - t0)

  starts <- lapply(draws, function(d) sort(unique(d$starts)))
  expect_true(all(vapply(starts, function(st) st[1] == 1, TRUE)),
              label = "T2 every partition has a segment starting on day 1")
  keys <- vapply(starts, ref_partition_key, "")
  cfr <- do.call(rbind, lapply(draws, function(d) as.numeric(d$cfr)))
  expect_equal(ncol(cfr), n_days)

  parts <- ref_enumerate_partitions(n_days, t2_l)
  pr <- ref_partition_posterior(parts, n_days, ref_cfr_log_seg(t2_F, t2_I, t2_alpha, t2_beta),
                                function(st) ref_log_prior_partition_p(st, n_days, t2_p))
  names(pr) <- vapply(parts, ref_partition_key, "")
  expect_gof(table(keys), pr, "T2 partition frequencies")

  # CFR constant within each drawn segment (deterministic)
  const_ok <- vapply(seq_len(N), function(i) {
    sg <- ref_segments(starts[[i]], n_days)
    all(vapply(seq_len(nrow(sg)), function(r) diff(range(cfr[i, sg$u[r]:sg$v[r]])) < 1e-12, TRUE))
  }, TRUE)
  expect(all(const_ok), sprintf("CFR varies within a drawn segment in %d draws", sum(!const_ok)))

  # E[CFR_t] for every day
  segmean <- function(u, v) (t2_alpha + sum(t2_F[u:v])) / (t2_alpha + t2_beta + sum(t2_I[u:v]))
  ex <- ref_expected_by_day(parts, pr, n_days, segmean)
  swapped <- ref_expected_by_day(parts, pr, n_days, function(u, v)
    (t2_beta + sum(t2_F[u:v])) / (t2_alpha + t2_beta + sum(t2_I[u:v])))
  for (t in seq_len(n_days)) expect_mean_z(cfr[, t], ex[t], sprintf("T2 E[CFR_%d]", t),
                                          sprintf("(if alpha/beta swapped within segments: %.4f)", swapped[t]))

  # within the modal partition, each segment's CFR is Beta(a + F, b + I - F)
  modal <- names(pr)[which.max(pr)]
  sel <- which(keys == modal)
  sg <- ref_segments(as.numeric(strsplit(modal, "-")[[1]]), n_days)
  for (r in seq_len(nrow(sg))) {
    Fs <- sum(t2_F[sg$u[r]:sg$v[r]]); Is <- sum(t2_I[sg$u[r]:sg$v[r]])
    a <- t2_alpha + Fs; b <- t2_beta + Is - Fs
    expect_ks(cfr[sel, sg$u[r]], function(q) stats::pbeta(q, a, b),
              sprintf("T2 CFR of segment %d-%d in modal partition", sg$u[r], sg$v[r]))
  }
})

test_that("T2 sample_cfr_changepoints with is_changepoints = FALSE is the conjugate beta", {
  set.seed(2003)
  N <- 4000
  draws <- lapply(seq_len(N), function(i) sample_cfr_changepoints(
    deaths_true = t2_F, cases_true = t2_I, p_change = t2_p, min_length = t2_l,
    cfr_alpha = t2_alpha, cfr_beta = t2_beta, is_changepoints = FALSE))
  cfr <- do.call(rbind, lapply(draws, function(d) as.numeric(d$cfr)))
  expect(all(apply(cfr, 1, function(x) diff(range(x)) < 1e-12)), "CFR differs between days with a single segment")
  a <- t2_alpha + sum(t2_F); b <- t2_beta + sum(t2_I) - sum(t2_F)
  expect_ks(cfr[, 1], function(q) stats::pbeta(q, a, b), "T2 single-segment CFR ~ Beta(a+sumF, b+sumI-sumF)")
  expect_mean_z(cfr[, 1], a / (a + b), "T2 single-segment E[CFR]")
})
