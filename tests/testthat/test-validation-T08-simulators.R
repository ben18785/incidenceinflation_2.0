skip_validation("fast")

# T8. The package's own simulators, tested against the generative model of the spec
# (permitted by the brief). Expected values are derived from the spec, not from the
# package.

# ---- T8a true_cases: NB renewal (eq. 8) -----------------------------------------
# Pearson residuals (I_t - R Lambda_t) / sqrt(R Lambda_t (1 + R Lambda_t / kappa)) for
# t >= 41 (so Lambda_t involves only returned days) have mean 0 and mean square 1.
# Power: generating with the (s-1, s] serial weights gives z ~ 29.
test_that("T8a true_cases follows the NB renewal equation with the spec's omega", {
  set.seed(8001)
  sp <- list(mean = 6, sd = 3); R <- 1.5; kap <- 50; Tn <- 70
  om <- ref_serial_weights(6, 3, 40)
  t0 <- proc.time()[["elapsed"]]
  series <- lapply(1:500, function(i) as.numeric(true_cases(
    Tn, function(t) rep(R, length(t)), kap, sp,
    initial_parameters = list(mean = 10, length = 10), serial_max = 40)))
  note("T8a: 500 series took %.1f s", proc.time()[["elapsed"]] - t0)
  expect_equal(unique(vapply(series, length, 0)), Tn)
  r <- unlist(lapply(series, function(I) {
    L <- ref_lambda(I, numeric(0), om); t <- 41:Tn
    (I[t] - R * L[t]) / sqrt(R * L[t] * (1 + R * L[t] / kap))
  }))
  expect_mean_z(r, 0, "T8a mean Pearson residual")
  expect_mean_z(r^2, 1, "T8a mean squared Pearson residual")
})

# ---- T8b/c observed_cases ----------------------------------------------------------
t8_rows <- function(df) {
  data.frame(t = col_of(df, "time_onset"), r = col_of(df, "time_reported"),
             c = col_of(df, "cases_reported"))
}

test_that("T8b observed_cases: E[C(t|r)] = I G(r - t + 1/2), counts non-decreasing", {
  set.seed(8002)
  I <- 20000
  d <- t8_rows(observed_cases(rep(I, 5), list(mean = 6, sd = 3), days_max_follow_up = 30))
  note("T8b: observed_cases with 5 onset days returned %d rows, report times %d to %d",
       nrow(d), min(d$r), max(d$r))
  d <- d[d$r >= d$t, ]
  expect_gt(nrow(d), 0)
  p <- ref_cdf_case(d$r - d$t, 6, 3)
  rows <- which(p > 1e-6 & p < 1 - 1e-6)
  zk <- z_bonferroni(length(rows))
  for (i in rows)
    expect_z(d$c[i], I * p[i], sqrt(I * p[i] * (1 - p[i])),
             sprintf("T8b onset %d report %d", d$t[i], d$r[i]), zmax = zk)
  mono <- tapply(seq_len(nrow(d)), d$t, function(ix) { o <- ix[order(d$r[ix])]; all(diff(d$c[o]) >= 0) })
  expect_true(all(mono), label = "T8b counts non-decreasing in report time")
})

test_that("T8b observed_cases: counts are binomial (variance)", {
  set.seed(8003)
  I <- 200; lag <- 5
  # no report is made after the last onset day (help page of observed_cases), so 10
  # onset days are passed, the later ones with zero cases, for lag 5 of onset day 1
  x <- vapply(1:1000, function(i) {
    d <- t8_rows(observed_cases(c(I, rep(0, 9)), list(mean = 6, sd = 3), days_max_follow_up = 30))
    v <- d$c[d$t == 1 & d$r == 1 + lag]
    if (length(v) != 1) stop(sprintf("expected one row for onset 1, report %d; got %d (report times %s)",
                                      1 + lag, length(v), paste(range(d$r), collapse = "-")))
    v
  }, 0)
  p <- ref_cdf_case(lag, 6, 3)
  expect_mean_z(x, I * p, "T8b mean at lag 5")
  expect_mean_z((x - I * p)^2, I * p * (1 - p), "T8b variance at lag 5")
})

t8c_check <- function(fun, label) {
  set.seed(8004)
  I <- 200; eta <- 0.3; rho <- 0.1
  # no report is made after the last onset day, so 30 onset days are passed, the later
  # ones with zero cases
  sims <- lapply(1:2000, function(i) t8_rows(fun(c(I, rep(0, 29)))))
  all_r <- sort(unique(unlist(lapply(sims, function(d) d$r[d$t == 1]))))
  note("%s: onset day 1 reported at times %d to %d", label, min(all_r), max(all_r))
  expect(any(all_r >= 2 & all_r <= 25), sprintf("%s: no reports of onset 1 between days 2 and 25", label))
  surv <- function(dl) ref_surv_case(dl, 6, 3)
  for (r in all_r[all_r >= 2 & all_r <= 25]) {
    x <- vapply(sims, function(d) { v <- d$c[d$t == 1 & d$r == r]; if (length(v)) v else NA }, 0)
    ex <- I * (1 - surv(r - 1)) + ref_wrong_mean_daily(r - 1, eta, I, rho, surv)
    expect_mean_z(x[!is.na(x)], ex, sprintf("%s mean count at report %d", label, r),
                  sprintf("(no wrong records: %.2f)", I * (1 - surv(r - 1))))
  }
}

test_that("T8c observed_cases with eta, rho: wrong records entered and withdrawn daily (A.7)", {
  t8c_check(function(I) observed_cases(I, list(mean = 6, sd = 3), days_max_follow_up = 30,
                                       revision_parameters = list(eta = 0.3, rho = 0.1)),
            "T8c observed_cases")
})

test_that("T8c observed_cases_revisions with eta, rho: same expectation", {
  t8c_check(function(I) observed_cases_revisions(I, list(mean = 6, sd = 3),
                                                 revision_parameters = list(eta = 0.3, rho = 0.1),
                                                 days_max_follow_up = 30),
            "T8c observed_cases_revisions")
})

# ---- T8d simulate_cases_deaths ---------------------------------------------------------
t8d_rows <- function(sim) {
  s <- sim$snapshots
  data.frame(t = col_of(s, "time_onset"), r = col_of(s, "time_reported"),
             c = col_of(s, "cases_reported"), d = col_of(s, "deaths_reported"))
}

# z checks of every row against exact probabilities of being counted as a case and as
# a death; var_extra_* add the variance of Poisson wrong records when present.
# Rows with no wrong records (extra = 0) are checked with the exact binomial tail at the
# same Bonferroni level (two-sided 2 pnorm(-zk)), because the normal approximation fails
# where few cases remain unreported (T8d.1). Rows with wrong records (T8d.6) keep the z
# check.
t8d_check <- function(d, I, pC, pD, label, extra_C = 0, extra_D = 0) {
  zk <- z_bonferroni(sum(pC > 1e-6 & pC < 1 - 1e-6) + sum(pD > 1e-6 & pD < 1 - 1e-6))
  one_row <- function(x, p, extra, what, i) {
    lab <- sprintf("%s %s, onset %d report %d", label, what, d$t[i], d$r[i])
    if (extra == 0) expect_binom_exact(x, I, p, lab, p_min = 2 * stats::pnorm(-zk))
    else expect_z(x, I * p + extra, sqrt(I * p * (1 - p) + extra), lab, zmax = zk)
  }
  for (i in seq_len(nrow(d))) {
    if (pC[i] > 1e-6 && pC[i] < 1 - 1e-6)
      one_row(d$c[i], pC[i], extra_C[min(i, length(extra_C))], "cases", i)
    if (pD[i] > 1e-6 && pD[i] < 1 - 1e-6)
      one_row(d$d[i], pD[i], extra_D[min(i, length(extra_D))], "deaths", i)
  }
  expect_true(all(d$d <= d$c), label = sprintf("%s deaths never exceed cases", label))
}

test_that("T8d.1 simulate_cases_deaths: independent delays from onset", {
  set.seed(8011)
  I <- 20000; cfr <- 0.4
  sim <- simulate_cases_deaths(rep(I, 3), report_days = c(1:45, 400),
                               reporting_parameters = list(mean = 6, sd = 3), cfr = cfr,
                               death_parameters = list(mean = 12, sd = 5), death_from = "onset",
                               triggered = 0)
  d <- t8d_rows(sim); d <- d[d$r >= d$t & d$r <= 45, ]
  dl <- d$r - d$t
  pC <- ref_cdf_case(dl, 6, 3)
  pD <- cfr * pC * ref_cdf_case(dl, 12, 5)
  t8d_check(d, I, pC, pD, "T8d.1")
  # truth: fatal cases binomial in the cases (chi-square over days), far report = truth
  tr <- sim$truth
  co <- col_of(tr, "cases_onset"); fo <- col_of(tr, "deaths_onset")
  ok <- co > 0
  x2 <- sum((fo[ok] - cfr * co[ok])^2 / (co[ok] * cfr * (1 - cfr)))
  expect_gte(stats::pchisq(x2, sum(ok), lower.tail = FALSE), P_MIN)
  far <- t8d_rows(sim); far <- far[far$r == 400, ]
  far <- far[order(far$t), ]
  expect_equal(far$c, col_of(tr, "cases_true")[match(far$t, col_of(tr, "time_onset"))])
  expect_equal(far$d, col_of(tr, "deaths_true")[match(far$t, col_of(tr, "time_onset"))])
})

test_that("T8d.2 simulate_cases_deaths: triggered = 1 (A.5, Step 1)", {
  set.seed(8012)
  I <- 20000; cfr <- 0.4
  sim <- simulate_cases_deaths(rep(I, 3), report_days = 1:45,
                               reporting_parameters = list(mean = 10, sd = 5), cfr = cfr,
                               death_parameters = list(mean = 6, sd = 3), death_from = "onset",
                               triggered = 1)
  d <- t8d_rows(sim); d <- d[d$r >= d$t, ]
  dl <- d$r - d$t
  SC <- ref_surv_case(dl, 10, 5); SD <- ref_surv_case(dl, 6, 3)
  pC <- (1 - cfr) * (1 - SC) + cfr * (1 - SC * SD)
  pD <- cfr * (1 - SD)
  t8d_check(d, I, pC, pD, "T8d.2")
})

test_that("T8d.3 simulate_cases_deaths: onset_to_death and recording delay (eq. 15)", {
  set.seed(8013)
  I <- 20000; cfr <- 0.4
  pd <- ref_pdie(9.3, 5.5, 36)
  sim <- simulate_cases_deaths(rep(I, 3), report_days = 1:60,
                               reporting_parameters = list(mean = 6, sd = 3), cfr = cfr,
                               death_parameters = list(mean = 3, sd = 2), death_from = "onset",
                               triggered = 0, onset_to_death = list(mean = 9.3, sd = 5.5))
  d <- t8d_rows(sim); d <- d[d$r >= d$t, ]
  dl <- d$r - d$t
  pC <- ref_cdf_case(dl, 6, 3)
  pD <- cfr * pC * (1 - ref_surv_death(dl, pd, 3, 2))
  t8d_check(d, I, pC, pD, "T8d.3")
})

# Eq. (16) in its half-day form: backlog survival exp(-lambda max(0, x - m + 1/2)), with
# the hazard counted from the middle of day d + m.
test_that("T8d.4 simulate_cases_deaths: backlog of deaths (eqs. 15-16)", {
  set.seed(8014)
  I <- 1e5; cfr <- 0.5; b <- 0.5; lam <- 0.5; m <- 5
  pd <- ref_pdie(9.3, 5.5, 36)
  sim <- simulate_cases_deaths(I, report_days = 1:70,
                               reporting_parameters = list(mean = 4, sd = 2), cfr = cfr,
                               death_parameters = list(mean = 3, sd = 2), death_from = "onset",
                               triggered = 0, onset_to_death = list(mean = 9.3, sd = 5.5),
                               recording_backlog = list(share = b, rate = lam, hold_days = m))
  d <- t8d_rows(sim); d <- d[d$r >= d$t, ]
  dl <- d$r - d$t
  pC <- ref_cdf_case(dl, 4, 2)
  pD <- cfr * pC * (1 - ref_surv_death(dl, pd, 3, 2, backlog_share = b, backlog_rate = lam,
                                       hold = m, backlog_half_day = TRUE))
  t8d_check(d, I, pC, pD, "T8d.4")
})

test_that("T8d.5 simulate_cases_deaths: cases dated by death (A.6)", {
  set.seed(8015)
  f <- 0.3; cfr <- 0.4
  pd <- ref_pdie(9.3, 5.5, 36)
  sim <- simulate_cases_deaths(rep(3000, 30), report_days = c(1:30, 400),
                               reporting_parameters = list(mean = 6, sd = 3), cfr = cfr,
                               death_parameters = list(mean = 3, sd = 2), death_from = "onset",
                               onset_to_death = list(mean = 9.3, sd = 5.5),
                               dated_by_death = list(fraction = f))
  tr <- sim$truth
  day <- col_of(tr, "time_onset"); Fo <- col_of(tr, "deaths_onset"); nd <- col_of(tr, "n_dated_by_death")
  Fo_by_day <- function(x) { v <- Fo[match(x, day)]; ifelse(is.na(v), 0, v) }
  zk <- z_bonferroni(26)
  for (dd in 5:30) {
    z <- seq_along(pd)
    q <- f * pd
    ex <- sum(Fo_by_day(dd - z) * q)
    vr <- sum(Fo_by_day(dd - z) * q * (1 - q))
    expect_z(nd[match(dd, day)], ex, sqrt(vr), sprintf("T8d.5 cases dated by death recorded for day %d", dd), zmax = zk)
  }
})

test_that("T8d.6 simulate_cases_deaths: wrong records (A.7)", {
  set.seed(8016)
  I <- 20000; cfr <- 0.4; eta <- 0.3; rho <- 0.1; delta <- 0.4
  sim <- simulate_cases_deaths(I, report_days = 1:40,
                               reporting_parameters = list(mean = 6, sd = 3), cfr = cfr,
                               death_parameters = list(mean = 12, sd = 5), death_from = "onset",
                               triggered = 0, wrong = list(eta = eta, rho = rho, delta = delta))
  d <- t8d_rows(sim); d <- d[d$r >= d$t, ]
  dl <- d$r - d$t
  surv <- function(x) ref_surv_case(x, 6, 3)
  pC <- ref_cdf_case(dl, 6, 3)
  pD <- cfr * pC * ref_cdf_case(dl, 12, 5)
  w <- vapply(dl, function(x) ref_wrong_mean_daily(x, eta, I, rho, surv), 0)
  t8d_check(d, I, pC, pD, "T8d.6", extra_C = w, extra_D = delta * w)
})

# ---- T8e sample_batches [API-only] ---------------------------------------------------------
test_that("T8e sample_batches matches its help page [API-only]", {
  set.seed(8020)
  bp <- list(beta = 0.05, zeta = 2, phi_mu = 5, phi_sigma = 2, rho_batch = 0.2)
  days <- 2000
  bs <- lapply(1:30, function(i) sample_batches(days, bp))
  all_b <- do.call(rbind, lapply(bs, as.data.frame))
  # footprints do not overlap and are in order
  ok <- vapply(bs, function(b) {
    b <- as.data.frame(b); if (nrow(b) < 2) return(TRUE)
    o <- order(b$start); s <- b$start[o]; w <- b$width[o]
    all(s[-1] > s[-length(s)] + w[-length(w)] - 1)
  }, TRUE)
  expect_true(all(ok), label = "T8e footprints do not overlap")
  inner <- all_b[all_b$start <= days - 30, ]
  sup <- 0:20
  expect_pmf(inner$width - 1, sup, stats::dpois(sup, bp$zeta), "T8e width - 1 ~ Poisson(zeta)")
  expect_ks(all_b$phi, function(q) ref_pgamma_ms(q, bp$phi_mu, bp$phi_sigma), "T8e intensity ~ gamma(mean, sd)")
  # number of batches: an independent simulation of the help page's description
  my_count <- function() {
    n <- 0; d <- 1
    while (d <= days) {
      if (stats::runif(1) < bp$beta) { n <- n + 1; d <- d + 1 + stats::rpois(1, bp$zeta) } else d <- d + 1
    }
    n
  }
  ref <- replicate(4000, my_count())
  cnt <- vapply(bs, function(b) nrow(as.data.frame(b)), 0)
  expect_z(mean(cnt), mean(ref), stats::sd(ref) / sqrt(length(cnt)), "T8e mean number of batches")
})
