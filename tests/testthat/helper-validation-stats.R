# Statistical checks, with thresholds fixed before the tests were first run.
#   z checks: |z| <= 4.
#   chi-square / KS: p >= 1e-4.
#   coverage: lower-tail binomial p >= 0.01.

Z_MAX <- 4
P_MIN <- 1e-4

# z statistic of an estimate against its exact value.
stat_z <- function(estimate, exact, se) (estimate - exact) / se

# Expect |z| <= Z_MAX; label describes the quantity, extra is printed on failure.
expect_z <- function(estimate, exact, se, label, extra = "", zmax = Z_MAX) {
  z <- stat_z(estimate, exact, se)
  testthat::expect(
    is.finite(z) && abs(z) <= zmax,
    sprintf("%s: estimate %.6g vs exact %.6g (SE %.3g), z = %.2f exceeds %.2f. %s",
            label, estimate, exact, se, z, zmax, extra)
  )
  invisible(z)
}

# Bonferroni threshold for a block of k z checks made row by row over a table, so
# that the block as a whole has false-alarm rate 1e-4.
z_bonferroni <- function(k) max(Z_MAX, stats::qnorm(1 - 1e-4 / (2 * max(k, 1))))

# Exact two-sided binomial check of one count x ~ Binomial(n, p): the doubled smaller
# tail, p_two = min(1, 2 min(P(X <= x), P(X >= x))), must be at least p_min. Use p_min =
# 2 pnorm(-zmax) to keep the level of the z check, which it replaces where the normal
# approximation fails (n min(p, 1 - p) small).
expect_binom_exact <- function(x, n, p, label, p_min = 2 * stats::pnorm(-Z_MAX)) {
  lo <- stats::pbinom(x, n, p)
  hi <- stats::pbinom(x - 1, n, p, lower.tail = FALSE)
  p_two <- min(1, 2 * min(lo, hi))
  testthat::expect(
    is.finite(p_two) && p_two >= p_min,
    sprintf("%s: count %.0f vs expected %.6g (n = %.0f, p = %.6g); exact two-sided p = %.3g < %.3g",
            label, x, n * p, n, p, p_two, p_min)
  )
  invisible(p_two)
}

# z check of a sample mean against an exact mean, SE from the sample.
expect_mean_z <- function(x, exact, label, extra = "") {
  expect_z(mean(x), exact, stats::sd(x) / sqrt(length(x)), label, extra)
}

# z check of a sample proportion against an exact probability.
expect_prop_z <- function(successes, n, p_exact, label, extra = "") {
  se <- sqrt(p_exact * (1 - p_exact) / n)
  if (se == 0) se <- 1 / n
  expect_z(successes / n, p_exact, se, label, extra)
}

# Chi-square goodness of fit of observed counts against exact probabilities.
# Cells with expected count below 5 are pooled into one cell; if that pooled cell
# is itself below 5 it is merged with the smallest remaining cell. Observations
# falling outside the support given (names not in 'probs') make the test fail.
stat_chisq_gof <- function(observed, probs) {
  # observed: named counts; probs: named probabilities summing to ~1
  keys <- union(names(probs), names(observed))
  o <- as.numeric(observed[keys]); o[is.na(o)] <- 0
  p <- as.numeric(probs[keys]); p[is.na(p)] <- 0
  n <- sum(o)
  if (any(o > 0 & p == 0)) return(list(p.value = 0, stat = Inf, df = NA,
                                       note = "observations outside the support"))
  e <- n * p
  small <- e < 5
  oc <- o[!small]; ec <- e[!small]
  if (any(small)) {
    po <- sum(o[small]); pe <- sum(e[small])
    if (pe < 5 && length(ec)) {
      i <- which.min(ec); po <- po + oc[i]; pe <- pe + ec[i]; oc <- oc[-i]; ec <- ec[-i]
    }
    oc <- c(oc, po); ec <- c(ec, pe)
  }
  stat <- sum((oc - ec)^2 / ec)
  df <- length(oc) - 1
  list(p.value = stats::pchisq(stat, df, lower.tail = FALSE), stat = stat, df = df,
       note = "")
}

expect_gof <- function(observed, probs, label, extra = "") {
  r <- stat_chisq_gof(observed, probs)
  testthat::expect(
    r$p.value >= P_MIN,
    sprintf("%s: chi-square = %.2f on %s df, p = %.3g < %g. %s %s",
            label, r$stat, r$df, r$p.value, P_MIN, r$note, extra)
  )
  invisible(r)
}

# Discrete draws against an exact pmf given on a support vector.
expect_pmf <- function(draws, support, pmf, label, extra = "") {
  probs <- stats::setNames(pmf / sum(pmf), as.character(support))
  obs <- table(factor(as.character(draws), levels = unique(c(names(probs), as.character(draws)))))
  expect_gof(obs, probs, label, extra)
}

# KS test against an exact continuous CDF.
expect_ks <- function(x, cdf, label, extra = "") {
  r <- suppressWarnings(stats::ks.test(x, cdf))
  testthat::expect(
    r$p.value >= P_MIN,
    sprintf("%s: KS D = %.4f, p = %.3g < %g (n = %d). %s",
            label, r$statistic, r$p.value, P_MIN, length(x), extra)
  )
  invisible(r)
}

# Batch-means standard error of the mean of an MCMC series.
stat_batch_se <- function(x, n_batches = 20) {
  n <- length(x)
  b <- floor(n / n_batches)
  x <- x[seq_len(b * n_batches)]
  m <- colMeans(matrix(x, nrow = b))
  stats::sd(m) / sqrt(n_batches)
}

expect_mcmc_mean <- function(x, exact, label, extra = "") {
  expect_z(mean(x), exact, stat_batch_se(x), label, extra)
}

# Lower-tail binomial check for coverage counts.
expect_coverage <- function(n_covered, n, nominal, label) {
  p <- stats::pbinom(n_covered, n, nominal)
  testthat::expect(
    p >= 0.01,
    sprintf("%s: %d of %d intervals contain the truth (nominal %.2f); lower-tail p = %.3g < 0.01",
            label, n_covered, n, nominal, p)
  )
  invisible(p)
}

# CDF and sampler from a density given on a fine grid (trapezoid rule).
ref_grid_dist <- function(x, log_dens) {
  d <- exp(log_dens - max(log_dens))
  cw <- c(0, cumsum((d[-1] + d[-length(d)]) / 2 * diff(x)))
  cw <- cw / cw[length(cw)]
  mean_x <- sum(((x[-1] * d[-1] + x[-length(x)] * d[-length(d)]) / 2) * diff(x)) /
    sum((d[-1] + d[-length(d)]) / 2 * diff(x))
  list(
    cdf = stats::approxfun(x, cw, yleft = 0, yright = 1, ties = "ordered"),
    rand = function(n) stats::approx(cw, x, xout = stats::runif(n), ties = "ordered")$y,
    mean = mean_x
  )
}

# Pull a column out of a data frame, failing informatively if it is absent.
col_of <- function(df, name) {
  if (!name %in% names(df)) stop(sprintf("column '%s' not found; columns are: %s",
                                         name, paste(names(df), collapse = ", ")))
  df[[name]]
}

# Print a line into the test output (informational probes and diagnostics).
note <- function(...) cat("\n[NOTE] ", sprintf(...), "\n", sep = "")
