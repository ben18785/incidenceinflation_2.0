# Independent individual-level simulator, written from spec Appendix A. Used as the
# only source of synthetic data in the medium (M*) and slow (S*) tiers; the
# package's simulators are only ever tested, never used to make data.
#
# Conventions (helper-validation-spec.R):
# - a case with onset t and reporting delay T_C (gamma, mean/sd) is in the report made
#   on day r >= t if T_C <= r - t + 1/2 (A.4, half-day convention);
# - without onset_to_death, a fatal case's death is recorded after its own gamma delay
#   T_D from onset (same convention), and it is shown as dead in report r once it is
#   both reported and its death recorded (API: mcmc_single 'deaths'; A.5).

# Draw continuous gamma delays by mean and sd.
sim_rgamma_ms <- function(n, mean, sd) stats::rgamma(n, shape = (mean / sd)^2, rate = mean / sd^2)

# Reporting delays of n cases of one onset day with a heavier tail (sensitivity_extensions.tex,
# S8): with probability share the delay is exponential with mean tail_mean, otherwise
# gamma, so that P(T > d + 1/2) = (1 - q_t){1 - G(d + 1/2)} + q_t exp{-(d + 1/2)/mu}.
sim_rdelay_tail <- function(n, delay, share = 0, tail_mean = 1) {
  tc <- sim_rgamma_ms(n, delay$mean, delay$sd)
  is_tail <- stats::runif(n) < share
  tc[is_tail] <- stats::rexp(sum(is_tail), 1 / tail_mean)
  tc
}

# Poisson or NB renewal (A.2): I_t ~ Pois(R_t Lambda_t) or NB(R_t Lambda_t, kappa).
sim_renewal <- function(R, history, omega, kappa = Inf) {
  I <- numeric(length(R))
  for (t in seq_along(R)) {
    mu <- R[t] * ref_lambda(c(I[seq_len(t - 1)], 0), history, omega)[t]   # Lambda_t uses days < t only
    I[t] <- if (is.finite(kappa)) stats::rnbinom(1, size = kappa, mu = mu) else stats::rpois(1, mu)
  }
  I
}

# Reports of cases (and optionally deaths) for onset days 1..length(I).
# report_days: the days on which reports are made; onset t appears in reports r >= t.
# cfr: NULL for cases only, or a CFR per onset day; death_delay: list(mean, sd).
# Returns list(cases = tibble(time_onset, time_reported, cases_reported),
#              deaths = tibble(time_onset, time_reported, deaths_reported) or NULL,
#              fatal = true fatal cases per onset day).
sim_reports <- function(I, report_days, delay, cfr = NULL, death_delay = NULL) {
  rows_c <- list(); rows_d <- list(); fatal <- integer(length(I))
  for (t in seq_along(I)) {
    rd <- report_days[report_days >= t]
    if (!length(rd)) next
    n <- I[t]
    tc <- sim_rgamma_ms(n, delay$mean, delay$sd)
    rows_c[[t]] <- data.frame(time_onset = t, time_reported = rd,
                              cases_reported = vapply(rd, function(r) sum(tc <= r - t + 0.5), 0))
    if (!is.null(cfr)) {
      dead <- stats::runif(n) < cfr[t]
      fatal[t] <- sum(dead)
      td <- sim_rgamma_ms(n, death_delay$mean, death_delay$sd)
      shown <- vapply(rd, function(r) sum(dead & tc <= r - t + 0.5 & td <= r - t + 0.5), 0)
      rows_d[[t]] <- data.frame(time_onset = t, time_reported = rd, deaths_reported = shown)
    }
  }
  out <- list(cases = tibble::as_tibble(do.call(rbind, rows_c)), deaths = NULL, fatal = fatal)
  out$cases[] <- lapply(out$cases, as.integer)
  if (!is.null(cfr)) {
    out$deaths <- tibble::as_tibble(do.call(rbind, rows_d))
    out$deaths[] <- lapply(out$deaths, as.integer)
  }
  out
}

# ---- Wrong records (A.7, eq. 23) --------------------------------------------------
# For an onset day t with I genuine cases and case-delay survival surv(delta)
# (= S_C(t | t + delta)): on each day k = t+1, ..., last, Poisson(eta I [S(k-1-t) -
# S(k-t)]) wrong records are entered (none on the onset day itself). Each is withdrawn
# with probability rho on each later day, independently, so a record entered on day k
# is present at report r >= k with probability (1 - rho)^(r - k). Each is shown as dead
# with probability delta (deaths model only). Returns, for each report day in rd, the
# wrong records present and those of them shown as dead.
sim_wrong_records <- function(t, I, rd, surv, eta, rho, delta = 0) {
  out <- list(present = numeric(length(rd)), dead = numeric(length(rd)))
  last <- max(rd)
  if (eta <= 0 || last < t + 1) return(out)
  k <- (t + 1):last
  n <- stats::rpois(length(k), eta * I * (surv(k - 1 - t) - surv(k - t)))
  if (!sum(n)) return(out)
  day <- rep(k, n)
  life <- stats::rgeom(length(day), rho)          # days survived after entry
  dead <- stats::runif(length(day)) < delta
  for (j in seq_along(rd)) {
    here <- day <= rd[j] & day + life >= rd[j]
    out$present[j] <- sum(here); out$dead[j] <- sum(here & dead)
  }
  out
}

# Cases-only reports with wrong records (the S1 generative model).
sim_reports_wrong <- function(I, report_days, delay, eta, rho) {
  surv <- function(d) ref_surv_case(d, delay$mean, delay$sd)
  rows <- lapply(seq_along(I), function(t) {
    rd <- report_days[report_days >= t]
    if (!length(rd)) return(NULL)
    tc <- sim_rgamma_ms(I[t], delay$mean, delay$sd)
    g <- vapply(rd, function(r) sum(tc <= r - t + 0.5), 0)
    w <- sim_wrong_records(t, I[t], rd, surv, eta, rho)
    data.frame(time_onset = t, time_reported = rd, cases_reported = g + w$present)
  })
  d <- tibble::as_tibble(do.call(rbind, rows)); d[] <- lapply(d, as.integer); d
}

# ---- Full deaths process (A.4-A.5, A.7), no misdating, no backlog --------------------
# Each case of onset day t is fatal with probability cfr[t]. A fatal case dies z ~ pdie
# days after onset; its death is recorded after a gamma delay Trec (rec_mean[t],
# rec_sd[t]) from the day of death, by report r if z <= r - t and
# Trec <= r - t - z + 1/2 (helper-validation-spec.R ref_surv_rec). A case is reported by r if its
# own delay T_C <= r - t + 1/2 or, if fatal and of the triggered kind (probability tau),
# once its death is recorded (A.5). A fatal case is shown as dead once reported and its
# death recorded. Wrong records as above, entered in proportion to S_C, each shown dead
# with probability delta. S8: tail_share (one per onset day, or NULL) gives each case of
# day t its case delay from sim_rdelay_tail, and S_C of the wrong records has the tail.
sim_reports_deaths <- function(I, report_days, delay, cfr, pdie, rec_mean, rec_sd,
                               tau = 1, eta = 0, rho = 0.5, delta = 0,
                               tail_share = NULL, tail_mean = 1) {
  q <- if (is.null(tail_share)) rep(0, length(I)) else tail_share
  rows <- list(); fatal <- integer(length(I))
  for (t in seq_along(I)) {
    rd <- report_days[report_days >= t]
    if (!length(rd)) next
    n <- I[t]
    surv <- function(d) ref_surv_case(d, delay$mean, delay$sd, tail_weight = q[t],
                                      tail_mean = tail_mean)
    tc <- if (is.null(tail_share)) sim_rgamma_ms(n, delay$mean, delay$sd)
    else sim_rdelay_tail(n, delay, q[t], tail_mean)
    is_f <- stats::runif(n) < cfr[t]; fatal[t] <- sum(is_f)
    z <- sample.int(length(pdie), n, replace = TRUE, prob = pdie)
    trec <- sim_rgamma_ms(n, rec_mean[t], rec_sd[t])
    trig <- stats::runif(n) < tau
    cr <- dr <- numeric(length(rd))
    for (j in seq_along(rd)) {
      x <- rd[j] - t
      own <- tc <= x + 0.5
      rec <- is_f & z <= x & trec <= x - z + 0.5
      rep <- own | (rec & trig)
      cr[j] <- sum(rep); dr[j] <- sum(rep & rec)
    }
    w <- sim_wrong_records(t, n, rd, surv, eta, rho, delta)
    rows[[t]] <- data.frame(time_onset = t, time_reported = rd,
                            cases_reported = cr + w$present, deaths_reported = dr + w$dead)
  }
  d <- tibble::as_tibble(do.call(rbind, rows)); d[] <- lapply(d, as.integer)
  list(cases = d[, c("time_onset", "time_reported", "cases_reported")],
       deaths = d[, c("time_onset", "time_reported", "deaths_reported")], fatal = fatal)
}

# ---- Full deaths process with misdated cases (appendix A, 'Cases recorded for their
# day of death'; sensitivity_extensions.tex, S9) ----------------------------------------
# As sim_reports_deaths, except that each fatal case of onset day t is misdated with
# probability fraction[t] (one value per onset day, so a fraction varying by onset week
# is given week by week). A misdated case dying z ~ pdie days after onset is recorded for
# its day of death d = t + z: it is in no report if d is after the last onset day
# (n_days); otherwise it enters the records when its death does, already shown as dead,
# in report r >= d if its recording delay Trec (gamma rec_mean[d], rec_sd[d], those of
# the day it is recorded for) has Trec <= r - d + 1/2. Fatal cases with onset before
# day 1 and dying on day d >= 1, misdated: G_0(d) ~ Poisson(sum_s history_s cfr[1]
# fraction[1] pdie(d - s)), history being the cases of the days before day 1 (the last
# is day 0), recorded as other misdated cases. Wrong records arise in proportion to the
# genuine cases recorded for the day, Itilde(d) = I(d) - (its misdated cases) +
# (misdated cases recorded for d), eq. (A.wrong-entered). Returns also, for each onset
# day, the true fatal cases and those misdated.
sim_reports_deaths_dated <- function(I, report_days, delay, cfr, pdie, rec_mean, rec_sd,
                                     fraction, history = NULL, tau = 1, eta = 0,
                                     rho = 0.5, delta = 0) {
  n_days <- length(I)
  surv <- function(d) ref_surv_case(d, delay$mean, delay$sd)
  zmax <- length(pdie)
  fatal <- misdated <- integer(n_days)
  # the misdated cases recorded for each day: their recording delays
  arrivals <- vector("list", n_days)
  f_arrive <- function(d, n) {
    if (n > 0) arrivals[[d]] <<- c(arrivals[[d]], sim_rgamma_ms(n, rec_mean[d], rec_sd[d]))
  }
  correct <- vector("list", n_days)
  for (t in seq_len(n_days)) {
    n <- I[t]
    tc <- sim_rgamma_ms(n, delay$mean, delay$sd)
    is_f <- stats::runif(n) < cfr[t]; fatal[t] <- sum(is_f)
    is_mis <- is_f & stats::runif(n) < fraction[t]; misdated[t] <- sum(is_mis)
    z <- sample.int(zmax, n, replace = TRUE, prob = pdie)
    for (d in (t + z[is_mis])[t + z[is_mis] <= n_days]) f_arrive(d, 1)
    keep <- !is_mis
    correct[[t]] <- list(tc = tc[keep], is_f = is_f[keep], z = z[keep],
                         trec = sim_rgamma_ms(sum(keep), rec_mean[t], rec_sd[t]),
                         trig = stats::runif(sum(keep)) < tau)
  }
  if (length(history)) {
    days_before <- seq_along(history) - length(history)       # ..., -1, 0
    for (d in seq_len(min(n_days, zmax))) {
      lag <- d - days_before
      ok <- lag >= 1 & lag <= zmax
      f_arrive(d, stats::rpois(1, sum(history[ok] * cfr[1] * fraction[1] * pdie[lag[ok]])))
    }
  }
  rows <- list()
  for (t in seq_len(n_days)) {
    rd <- report_days[report_days >= t]
    if (!length(rd)) next
    k <- correct[[t]]; a <- arrivals[[t]]
    cr <- dr <- numeric(length(rd))
    for (j in seq_along(rd)) {
      x <- rd[j] - t
      own <- k$tc <= x + 0.5
      rec <- k$is_f & k$z <= x & k$trec <= x - k$z + 0.5
      rep <- own | (rec & k$trig)
      mis <- sum(a <= x + 0.5)
      cr[j] <- sum(rep) + mis; dr[j] <- sum(rep & rec) + mis
    }
    w <- sim_wrong_records(t, length(k$tc) + length(a), rd, surv, eta, rho, delta)
    rows[[t]] <- data.frame(time_onset = t, time_reported = rd,
                            cases_reported = cr + w$present, deaths_reported = dr + w$dead)
  }
  d <- tibble::as_tibble(do.call(rbind, rows)); d[] <- lapply(d, as.integer)
  list(cases = d[, c("time_onset", "time_reported", "cases_reported")],
       deaths = d[, c("time_onset", "time_reported", "deaths_reported")],
       fatal = fatal, misdated = misdated,
       recorded_misdated = vapply(arrivals, length, 0))
}

# Report days with gaps drawn from gaps (with replacement), from first up to and
# including last.
sim_report_days <- function(first, last, gaps) {
  r <- first
  while (r[length(r)] < last) r <- c(r, r[length(r)] + sample(gaps, 1))
  r <- r[r < last]
  as.integer(c(r, last))
}

# Draw a partition of days 1..n_days from the product-partition prior with p ~
# Beta(a, b) integrated out and segments of at least min_length days (A.3, eq. 11):
# P(K) propto N_l(K) B(a + K, b + T - 1 - K), then uniform among partitions with K
# changes. Returns the segment starts.
sim_partition_prior <- function(n_days, a, b, min_length) {
  parts <- ref_enumerate_partitions(n_days, min_length)
  K <- vapply(parts, length, 0) - 1
  w <- exp(lbeta(a + K, b + n_days - 1 - K) - max(lbeta(a + K, b + n_days - 1 - K)))
  parts[[sample.int(length(parts), 1, prob = w)]]
}
