skip_validation("fast")

# T1. sample_changepoints draws a partition from its prior given p (spec A.3, eq. 9):
#   pi(P | p) proportional to p^K (1 - p)^(T - 1 - K) on A_l
# (partitions whose segments all have at least l days). Expected values come from
# enumerating every partition (helper-validation-spec.R), checked in plain R against the
# spec's closed form N_l(k).

# the days on which segments start, so day 1 is always one of them
t1_draw_keys <- function(n, n_days, p, l) {
  vapply(seq_len(n), function(i) {
    st <- sort(unique(sample_changepoints(n_days, p, min_length = l)))
    if (length(st) == 0 || st[1] != 1) return("no start on day 1")
    ref_partition_key(st)
  }, character(1))
}

for (cfg in list(c(l = 1, p = 0.2), c(l = 1, p = 0.5), c(l = 3, p = 0.2), c(l = 3, p = 0.5))) {
  test_that(sprintf("T1 sample_changepoints follows eq. 9 (T = 9, l = %d, p = %.1f)", cfg[["l"]], cfg[["p"]]), {
    set.seed(1000 + 10 * cfg[["l"]] + 100 * cfg[["p"]])
    n_days <- 9; l <- cfg[["l"]]; p <- cfg[["p"]]; N <- 20000
    t0 <- proc.time()[["elapsed"]]
    keys <- t1_draw_keys(N, n_days, p, l)
    note("T1 (l=%d, p=%.1f): %d draws took %.1f s", l, p, N, proc.time()[["elapsed"]] - t0)

    parts <- ref_enumerate_partitions(n_days, l)
    lw <- vapply(parts, ref_log_prior_partition_p, numeric(1), n_days = n_days, p = p)
    probs <- stats::setNames(exp(lw) / sum(exp(lw)), vapply(parts, ref_partition_key, ""))
    # every draw must be an allowed partition
    bad <- setdiff(unique(keys), names(probs))
    expect(length(bad) == 0, sprintf("partitions outside A_l drawn: %s", paste(head(bad), collapse = "; ")))
    expect_gof(table(keys), probs, sprintf("T1 partitions (l=%d, p=%.1f)", l, p))

    # marginal of K against N_l(k) p^k (1-p)^(T-1-k)
    K <- vapply(strsplit(keys, "-"), length, 1) - 1
    ks <- 0:(n_days - 1)
    pk <- ref_n_partitions(n_days, l, ks) * p^ks * (1 - p)^(n_days - 1 - ks)
    pk <- pk / sum(pk)
    for (k in ks[pk > 1e-3]) expect_prop_z(sum(K == k), N, pk[k + 1], sprintf("T1 P(K=%d) (l=%d, p=%.1f)", k, l, p))
  })
}
