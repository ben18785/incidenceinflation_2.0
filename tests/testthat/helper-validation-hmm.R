# Independent forward filter for the reports of ONE onset day, cases only, with
# wrong records (spec A.7 and A.8). Written from the spec; nothing here touches the
# package under test.
#
# Hidden state at report j: E_j, the wrong records present. The genuine cases
# reported are then G_j = C_j - E_j (the observed count less the wrong records).
#
# First report (made on day t1 > t):
#   G_1 ~ Binomial(N, 1 - S(t | t1))                            (eq. 17, S_0 = 1)
#   E_1 ~ Poisson(eta N sum_{k=t+1}^{t1} [S(k-1) - S(k)] (1-rho)^(t1-k))
#                                               (A.7, "The first report")
# Report j - 1 -> j:
#   G_j - G_{j-1} ~ Binomial(N - G_{j-1}, 1 - S_j / S_{j-1})       (eq. 17)
#   n_j ~ Poisson(eta N sum_{k=t_{j-1}+1}^{t_j} [S(k-1) - S(k)] (1-rho)^(t_j-k))
#     with entry_daily = TRUE (eq. 23: wrong records entered since the report
#     before and still present), or, with entry_daily = FALSE, an earlier form
#     of eq. 23 in which every entry since t_{j-1} is present at t_j:
#   n_j ~ Poisson(eta N [S_{j-1} - S_j])
#   e_j ~ Binomial(E_{j-1}, 1 - (1 - rho)^(t_j - t_{j-1}))           (eq. 24)
#   E_j = E_{j-1} + n_j - e_j
# States with E_j > e_max are given probability zero (one reading of the cap;
# tests using this filter are designed so the cap is practically never reached).
#
# Returns log L_t(N) = log Pr(C_1, ..., C_J | N).

ref_hmm_loglik_cases <- function(counts, report_days, onset, N, surv, eta, rho,
                                 e_max = 50, entry_daily = FALSE) {
  J <- length(counts)
  S <- surv(report_days - onset)
  delta1 <- report_days[1] - onset
  # --- first report
  states <- 0:e_max
  alpha <- rep(0, length(states))
  mu1 <- ref_wrong_mean_daily(delta1, eta, N, rho, surv)
  for (i in seq_along(states)) {
    E <- states[i]; G <- counts[1] - E
    if (G < 0 || G > N) next
    alpha[i] <- stats::dbinom(G, N, 1 - S[1]) * stats::dpois(E, mu1)
  }
  if (J > 1) for (j in 2:J) {
    p_rep <- 1 - S[j] / S[j - 1]
    mu_n <- if (entry_daily) {
      k <- (report_days[j - 1] + 1):report_days[j]
      eta * N * sum((surv(k - 1 - onset) - surv(k - onset)) * (1 - rho)^(report_days[j] - k))
    } else eta * N * (S[j - 1] - S[j])
    r_w <- 1 - (1 - rho)^(report_days[j] - report_days[j - 1])
    new <- rep(0, length(states))
    for (i_new in seq_along(states)) {
      E <- states[i_new]; G <- counts[j] - E
      if (G < 0 || G > N) next
      tot <- 0
      for (i_old in seq_along(states)) {
        if (alpha[i_old] == 0) next
        Eo <- states[i_old]; Go <- counts[j - 1] - Eo
        a <- G - Go
        if (a < 0) next
        pa <- stats::dbinom(a, N - Go, p_rep)
        if (pa == 0) next
        # sum over withdrawals e: E = Eo + n - e, n >= 0
        e <- 0:Eo
        n <- E - Eo + e
        ok <- n >= 0
        pe <- sum(stats::dbinom(e[ok], Eo, r_w) * stats::dpois(n[ok], mu_n))
        tot <- tot + alpha[i_old] * pa * pe
      }
      new[i_new] <- tot
    }
    alpha <- new
  }
  log(sum(alpha))
}

# Simulate the reports of one onset day from the same process (used only to
# self-test the filter). mode "pdf" follows eqs. 23-24 between reports and the
# day-by-day first-report formula; mode "daily" withdraws records day by day
# throughout.
ref_sim_one_onset_cases <- function(N, report_days, onset, surv, eta, rho,
                                    mode = c("pdf", "daily")) {
  mode <- match.arg(mode)
  last <- max(report_days)
  days <- (onset + 1):last
  # genuine report day of each case: first day k (>= onset) with T_C <= k - onset + 1/2
  u <- stats::runif(N)
  # S is non-increasing in delta; find the first delta with S(delta) < u
  grid <- 0:(last - onset)
  Sg <- surv(grid)
  rep_delta <- vapply(u, function(x) { w <- which(Sg < x); if (length(w)) grid[w[1]] else Inf }, numeric(1))
  genuine <- vapply(report_days, function(r) sum(rep_delta <= r - onset), numeric(1))
  # wrong records: entries on each day k > onset, Poisson(eta N [S(k-1) - S(k)])
  ent_day <- c();
  for (k in days) {
    m <- eta * N * (surv(k - 1 - onset) - surv(k - onset))
    n <- stats::rpois(1, m)
    if (n > 0) ent_day <- c(ent_day, rep(k, n))
  }
  wrong <- numeric(length(report_days))
  if (mode == "daily") {
    # each record survives each later day with prob 1 - rho
    life <- if (length(ent_day)) stats::rgeom(length(ent_day), rho) else numeric(0)
    # record entered on day k, withdrawn on day k + life + 1 (present through k + life)
    for (j in seq_along(report_days)) {
      r <- report_days[j]
      wrong[j] <- sum(ent_day <= r & ent_day + life >= r)
    }
  } else {
    # pdf: daily thinning up to the first report; afterwards only records present at
    # the previous report can be withdrawn, each with prob 1 - (1-rho)^gap.
    r1 <- report_days[1]
    pre <- ent_day[ent_day <= r1]
    present <- sum(stats::runif(length(pre)) < (1 - rho)^(r1 - pre))
    wrong[1] <- present
    if (length(report_days) > 1) for (j in 2:length(report_days)) {
      gap <- report_days[j] - report_days[j - 1]
      removed <- stats::rbinom(1, present, 1 - (1 - rho)^gap)
      added <- sum(ent_day > report_days[j - 1] & ent_day <= report_days[j])
      present <- present - removed + added
      wrong[j] <- present
    }
  }
  genuine + wrong
}
