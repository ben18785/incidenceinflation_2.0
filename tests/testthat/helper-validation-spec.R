# Independent reference implementations of the model in Appendix A of the epi
# paper (bvd_cfr_paper/paper_epi; "spec" in these files). They were written
# from the spec alone, by a checker who never saw the package's code, and they
# never call the package: expected values in the validation tests come from
# here, from helper-validation-sim.R and from closed forms, never from package
# output. Section and equation numbers refer to the spec as it was then.

# ---- Gamma parameterised by mean and sd ---------------------------------------
# The spec writes every gamma as (mean, sd), e.g. Table 3; and the API describes
# delays by 'mean' and 'sd'. shape = m^2/s^2, rate = m/s^2.
ref_gamma_shape <- function(mean, sd) (mean / sd)^2
ref_gamma_rate  <- function(mean, sd) mean / sd^2
ref_pgamma_ms <- function(x, mean, sd) {
  stats::pgamma(x, shape = ref_gamma_shape(mean, sd), rate = ref_gamma_rate(mean, sd))
}
ref_dgamma_ms <- function(x, mean, sd, log = FALSE) {
  stats::dgamma(x, shape = ref_gamma_shape(mean, sd), rate = ref_gamma_rate(mean, sd), log = log)
}

# ---- Serial interval (A.2) ----------------------------------------------------
# omega_s = mass of the gamma serial interval on (s - 1/2, s + 1/2),
# s = 1..serial_max, renormalised to sum to one.
ref_serial_weights <- function(mean, sd, serial_max = 40) {
  s <- seq_len(serial_max)
  w <- ref_pgamma_ms(s + 0.5, mean, sd) - ref_pgamma_ms(s - 0.5, mean, sd)
  w / sum(w)
}

# Lambda_t = sum_{s=1}^{serial_max} omega_s I(t - s)   (eq. 8)
# 'history' is the cases before onset day 1, earliest first (as cases_before);
# 'cases' are the cases on onset days 1..T. Missing history counts as zero.
ref_lambda <- function(cases, history, omega) {
  x <- c(history, cases)
  off <- length(history)
  vapply(seq_along(cases), function(t) {
    idx <- off + t - seq_along(omega)
    ok <- idx >= 1
    sum(omega[ok] * x[idx[ok]])
  }, numeric(1))
}

# ---- Negative binomial with mean mu and overdispersion kappa (A.2) -------------
# Variance mu (1 + mu / kappa), i.e. size = kappa.
ref_dnb <- function(x, mu, kappa, log = FALSE) {
  stats::dnbinom(x, size = kappa, mu = mu, log = log)
}

# ---- Time from onset to death (A.4, Table 3) -----------------------------------
# p_die(z) = mass of the gamma on (z - 1, z], z = 1..zmax, rescaled to sum to one.
ref_pdie <- function(mean = 9.3, sd = 5.5, zmax = 36) {
  z <- seq_len(zmax)
  p <- ref_pgamma_ms(z, mean, sd) - ref_pgamma_ms(z - 1, mean, sd)
  p / sum(p)
}

# S5 of the sensitivity analyses (sensitivity_extensions.tex): a misdated case is
# recorded for the day z1 days after onset and dies z2 days after that day, z1 and z2
# independent. With T_die ~ Gamma(a, c) (mean, sd) and the share of the way to death
# at which it is recorded Beta(a1, a - a1), z1 is Gamma(a1, c) on (z1 - 1, z1], z1 >= 1,
# and z2 is Gamma(a - a1, c) on (z2, z2 + 1], z2 >= 0, each truncated at its 99.9th
# percentile and rescaled. a1 = share a (by default a / 2).
ref_position_delays <- function(mean = 9.3, sd = 5.5, share = 1 / 2) {
  a <- ref_gamma_shape(mean, sd); c <- ref_gamma_rate(mean, sd); a1 <- share * a
  f_days <- function(shape) {
    n <- max(1, ceiling(stats::qgamma(0.999, shape, c)))
    p <- stats::pgamma(seq_len(n), shape, c) - stats::pgamma(seq_len(n) - 1, shape, c)
    p / sum(p)
  }
  list(p_recorded = f_days(a1), p_after = if (a1 < a) f_days(a - a1) else 1)
}

# ---- Reporting of cases (A.4) ---------------------------------------------------
# Probability a case with onset on day t has NOT been reported by the end of day
# t' = t + delta. Onset is taken at the middle of day t, so
#   S(t | t') = 1 - G(t' - t + 1/2),
# with the optional exponential tail of eq. (13).
ref_surv_case <- function(delta, mean, sd, tail_weight = 0, tail_mean = 1) {
  x <- delta + 0.5
  (1 - tail_weight) * (1 - ref_pgamma_ms(x, mean, sd)) +
    tail_weight * exp(-x / tail_mean)
}
ref_cdf_case <- function(delta, mean, sd, ...) 1 - ref_surv_case(delta, mean, sd, ...)

# ---- Death and its recording (A.4, eqs. 15 and 16) -------------------------------
# Survival of the recording delay of a death that occurred on day d, at the end of
# day t' (delta_rec = t' - d). Prompt group: gamma delay counted from the death,
# with the half-day convention of eq. (14) ("applied as in eq. (14)"); backlog
# group: hazard lambda per day from day d + m + 1 on (constant lambda here):
#   S_rec,B = exp(-sum_{s=d+m+1}^{t'} lambda) = exp(-lambda * max(0, t' - d - m)).
# Mixture with share b (eq. 16).
ref_surv_rec <- function(delta_rec, mean, sd, backlog_share = 0, backlog_rate = 0,
                         hold = 0, backlog_half_day = FALSE) {
  prompt <- 1 - ref_pgamma_ms(delta_rec + 0.5, mean, sd)
  expo <- pmax(0, delta_rec - hold + if (backlog_half_day) 0.5 else 0)
  backlog <- exp(-backlog_rate * expo)
  (1 - backlog_share) * prompt + backlog_share * backlog
}

# S_D(t | t') for an onset-to-recording delay made of the time to death (p_die)
# and the recording delay from the day of death (eq. 15).
ref_surv_death <- function(delta, pdie, rec_mean, rec_sd, ...) {
  vapply(delta, function(dl) {
    z <- seq_along(pdie)
    alive <- sum(pdie[z > dl])
    zz <- z[z <= dl]
    dead_unrec <- if (length(zz)) sum(pdie[zz] * ref_surv_rec(dl - zz, rec_mean, rec_sd, ...)) else 0
    alive + dead_unrec
  }, numeric(1))
}

# ---- Wrong records present at the first report (A.7, "The first report") --------
# E(t | t1) ~ Poisson(eta * Itilde * sum_{k=t+1}^{t1} [S(t|k-1) - S(t|k)] (1-rho)^(t1-k))
# surv_fun(delta) gives S(t | t + delta).
ref_wrong_mean_daily <- function(delta1, eta, n_genuine, rho, surv_fun) {
  if (delta1 < 1) return(0)
  k <- seq_len(delta1)                  # k - t = 1..delta1
  inc <- surv_fun(k - 1) - surv_fun(k)
  eta * n_genuine * sum(inc * (1 - rho)^(delta1 - k))
}

# ---- Partitions (A.3) -------------------------------------------------------------
# All partitions of days 1..n_days into consecutive segments each of at least
# min_length days (first segment at least first_min days). Each partition is
# returned as the vector of its segment starts (always including day 1).
ref_enumerate_partitions <- function(n_days, min_length = 1, first_min = min_length) {
  out <- list()
  rec <- function(starts, next_start) {
    # current segment begins at the last element of starts
    cur <- starts[length(starts)]
    need <- if (length(starts) == 1) first_min else min_length
    # option 1: current segment runs to the end
    if (n_days - cur + 1 >= need) out[[length(out) + 1]] <<- starts
    # option 2: start a new segment at s
    lo <- cur + need
    hi <- n_days - min_length + 1
    if (lo <= hi) for (s in lo:hi) rec(c(starts, s))
  }
  rec(1)
  out
}

ref_partition_key <- function(starts) paste(starts, collapse = "-")

# N_l(k) = choose(T - l(k+1) + k, k)  (A.3, Step 3)
ref_n_partitions <- function(n_days, min_length, k) {
  E <- n_days - min_length * (k + 1)
  ifelse(E < 0, 0, choose(E + k, k))
}

# Segment bounds of a partition given by its starts.
ref_segments <- function(starts, n_days) {
  ends <- c(starts[-1] - 1, n_days)
  data.frame(u = starts, v = ends)
}

# Log prior of a partition given p (eq. 9, conditional on p).
ref_log_prior_partition_p <- function(starts, n_days, p) {
  K <- length(starts) - 1
  K * log(p) + (n_days - 1 - K) * log1p(-p)
}

# Log marginal prior with p ~ Beta(a, b) integrated out (eq. 11).
ref_log_prior_partition_marg <- function(K, n_days, a, b) {
  lbeta(a + K, b + n_days - 1 - K)
}

# Exact posterior over an enumerated list of partitions, given a function
# log_seg(u, v) giving the log marginal likelihood of one segment, and a log prior
# function of the starts. Returns normalised probabilities.
ref_partition_posterior <- function(parts, n_days, log_seg, log_prior) {
  # cache segment terms
  cache <- new.env()
  seg_val <- function(u, v) {
    key <- paste(u, v)
    if (!is.null(cache[[key]])) return(cache[[key]])
    val <- log_seg(u, v)
    assign(key, val, envir = cache)
    val
  }
  lw <- vapply(parts, function(st) {
    sg <- ref_segments(st, n_days)
    log_prior(st) + sum(mapply(seg_val, sg$u, sg$v))
  }, numeric(1))
  w <- exp(lw - max(lw))
  w / sum(w)
}

# ---- CFR segment marginal likelihood (A.10.4, Step 1) ------------------------------
# Up to the binomial coefficients, which are common to every partition:
#   log B(alpha + F_uv, beta + I_uv - F_uv) - log B(alpha, beta)
ref_cfr_log_seg <- function(deaths, cases, alpha, beta) {
  function(u, v) {
    Fs <- sum(deaths[u:v]); Is <- sum(cases[u:v])
    lbeta(alpha + Fs, beta + Is - Fs) - lbeta(alpha, beta)
  }
}

# ---- R segment marginal likelihood, Poisson renewal (A.10.5, closed form) -----------
# int prod_t Pois(I_t | R Lambda_t) Gamma(R | shape, rate) dR
ref_rt_log_seg_pois <- function(cases, lambda, shape, rate) {
  function(u, v) {
    I <- cases[u:v]; L <- lambda[u:v]
    sum(I * log(L) - lgamma(I + 1)) + shape * log(rate) - lgamma(shape) +
      lgamma(shape + sum(I)) - (shape + sum(I)) * log(rate + sum(L))
  }
}

# Same under the negative binomial renewal, by numerical quadrature.
ref_rt_log_seg_nb <- function(cases, lambda, shape, rate, kappa) {
  function(u, v) {
    I <- cases[u:v]; L <- lambda[u:v]
    f_log <- function(R) vapply(R, function(r)
      sum(ref_dnb(I, r * L, kappa, log = TRUE)) +
        stats::dgamma(r, shape, rate, log = TRUE), numeric(1))
    # locate the mode for scaling
    opt <- stats::optimize(function(r) -f_log(r), c(1e-6, 50))
    m <- -opt$objective
    val <- stats::integrate(function(r) exp(f_log(r) - m), 0, Inf,
                            rel.tol = 1e-10, subdivisions = 2000L)$value
    log(val) + m
  }
}

# Posterior mean of R in a segment, NB renewal, by quadrature.
ref_rt_seg_mean_nb <- function(cases, lambda, shape, rate, kappa, u, v) {
  I <- cases[u:v]; L <- lambda[u:v]
  f_log <- function(R) vapply(R, function(r)
    sum(ref_dnb(I, r * L, kappa, log = TRUE)) +
      stats::dgamma(r, shape, rate, log = TRUE), numeric(1))
  opt <- stats::optimize(function(r) -f_log(r), c(1e-6, 50))
  m <- -opt$objective
  num <- stats::integrate(function(r) r * exp(f_log(r) - m), 0, Inf, rel.tol = 1e-10,
                          subdivisions = 2000L)$value
  den <- stats::integrate(function(r) exp(f_log(r) - m), 0, Inf, rel.tol = 1e-10,
                          subdivisions = 2000L)$value
  num / den
}

# Expected value of a per-day quantity under an enumerated partition posterior,
# where seg_mean(u, v) gives the conditional mean within a segment.
ref_expected_by_day <- function(parts, probs, n_days, seg_mean) {
  out <- numeric(n_days)
  for (i in seq_along(parts)) {
    sg <- ref_segments(parts[[i]], n_days)
    for (r in seq_len(nrow(sg))) {
      out[sg$u[r]:sg$v[r]] <- out[sg$u[r]:sg$v[r]] + probs[i] * seg_mean(sg$u[r], sg$v[r])
    }
  }
  out
}

# ---- Partition posterior with p integrated out, by dynamic programming -------------
# pi(P | data) proportional to B(a + K, b + T - 1 - K) prod_seg exp(log_seg(u, v)),
# segments at least min_length long (first at least first_min). Returns the
# posterior of K and the posterior mean of a per-day quantity whose conditional
# mean within segment (u, v) is seg_mean(u, v). (Eqs. 11 and 28.)
ref_partition_dp <- function(n_days, log_seg, seg_mean, a, b, min_length = 1,
                             first_min = min_length, k_max = n_days) {
  Tn <- n_days
  LS <- matrix(-Inf, Tn, Tn); MS <- matrix(NA_real_, Tn, Tn)
  for (u in 1:Tn) for (v in u:Tn) {
    need <- if (u == 1) first_min else min_length
    if (v - u + 1 >= need) { LS[u, v] <- log_seg(u, v); MS[u, v] <- seg_mean(u, v) }
  }
  lse <- function(x) { m <- max(x); if (!is.finite(m)) return(-Inf); m + log(sum(exp(x - m))) }
  # forward: f[v, k+1] = log sum over partitions of 1..v with k changes
  f <- matrix(-Inf, Tn, k_max + 1)
  for (v in 1:Tn) {
    f[v, 1] <- LS[1, v]
    if (v >= 2) for (k in 1:min(k_max, v - 1)) {
      terms <- vapply(2:v, function(u) f[u - 1, k] + LS[u, v], numeric(1))
      f[v, k + 1] <- lse(terms)
    }
  }
  # backward: g[u, k+1] = log sum over partitions of u..T (u starting a segment,
  # u > 1 so min_length applies) with k further changes after the segment at u
  g <- matrix(-Inf, Tn + 1, k_max + 1)
  for (u in Tn:2) {
    g[u, 1] <- LS[u, Tn]
    for (k in 1:k_max) {
      if (u + 1 > Tn) break
      terms <- vapply((u + 1):Tn, function(w) LS[u, w - 1] + g[w, k], numeric(1))
      g[u, k + 1] <- lse(terms)
    }
  }
  lprior <- lbeta(a + (0:k_max), b + Tn - 1 - (0:k_max))
  logZk <- f[Tn, ] + lprior
  logZ <- lse(logZk)
  pK <- exp(logZk - logZ)
  # E[quantity on day t] = sum over segments (u, v) containing t of
  #   P(segment u..v is in the partition) * seg_mean(u, v)
  # P(u..v) = sum_{k1, k2} exp(f_before(u-1, k1) + LS[u,v] + g_after(v+1, k2)
  #           + lprior[k1 + k2 + (u>1) + (v<T)]) / Z
  Ed <- numeric(Tn)
  for (u in 1:Tn) for (v in u:Tn) {
    if (!is.finite(LS[u, v])) next
    lw <- c()
    k1_range <- if (u == 1) 0 else 0:(k_max)
    for (k1 in k1_range) {
      before <- if (u == 1) 0 else f[u - 1, k1 + 1]
      if (!is.finite(before)) next
      k2_range <- if (v == Tn) 0 else 0:k_max
      for (k2 in k2_range) {
        after <- if (v == Tn) 0 else g[v + 1, k2 + 1]
        if (!is.finite(after)) next
        K <- k1 + k2 + (u > 1) + (v < Tn)
        if (K > k_max) next
        lw <- c(lw, before + LS[u, v] + after + lprior[K + 1])
      }
    }
    if (!length(lw)) next
    pr <- exp(lse(lw) - logZ)
    Ed[u:v] <- Ed[u:v] + pr * MS[u, v]
  }
  list(pK = pK, mean_by_day = Ed)
}

# ---- Binomial reporting chain (A.5, eq. 17) ----------------------------------------
# Log-likelihood of the cumulative reported counts 'counts' at report days
# 'report_days' for one onset day 'onset', given N true cases and a survival
# function surv(delta). S_0 = 1 before the first report (A.10.2).
ref_chain_loglik <- function(counts, report_days, onset, N, surv) {
  S <- surv(report_days - onset)
  Sprev <- c(1, S[-length(S)])
  A <- counts; Aprev <- c(0, A[-length(A)])
  a <- A - Aprev
  if (any(a < 0) || A[length(A)] > N) return(-Inf)
  # once S has reached 0 every remaining case is certainly reported
  p <- ifelse(Sprev > 0, 1 - S / Sprev, 1)
  ll <- sum(stats::dbinom(a, N - Aprev, p, log = TRUE))
  if (is.nan(ll)) -Inf else ll
}
