skip_validation("fast")

# T9. Deterministic utilities with hand-computed expected values.

test_that("T9 piecewise_values expands segment values", {
  expect_equal(as.numeric(piecewise_values(c(1, 4, 8), c(10, 20, 30), 10)),
               c(10, 10, 10, 20, 20, 20, 20, 30, 30, 30))
  expect_equal(as.numeric(piecewise_values(1, 5, 4)), rep(5, 4))
})

test_that("T9 piecewise_reporting_parameters gives each onset time its segment's delay", {
  out <- piecewise_reporting_parameters(c(1, 3), c(5, 9), c(2, 4), 5)
  expect_equal(as.numeric(col_of(out, "time_onset")), 1:5)
  expect_equal(as.numeric(col_of(out, "mean")), c(5, 5, 9, 9, 9))
  expect_equal(as.numeric(col_of(out, "sd")), c(2, 2, 4, 4, 4))
})

test_that("T9 max_downward_revision is the largest fall from ANY earlier count", {
  df <- tibble::tibble(
    time_onset = c(1, 1, 1, 1, 2, 2, 2),
    time_reported = c(2, 3, 4, 5, 3, 4, 5),
    cases_reported = c(10, 9, 8, 7, 5, 6, 4))
  # onset 1: 10 -> 7 is a fall of 3 (consecutive falls are only 1); onset 2: 6 -> 4 = 2
  expect_equal(as.numeric(max_downward_revision(df)), 3)
  up <- tibble::tibble(time_onset = c(1, 1, 2, 2), time_reported = c(2, 3, 3, 4),
                       cases_reported = c(1, 2, 3, 3))
  expect_equal(as.numeric(max_downward_revision(up)), 0)
})

test_that("T9 split_corrected_reports (spec 2.1 and A.7)", {
  df <- tibble::tibble(time_onset = rep(1:3, each = 3), time_reported = rep(c(5, 10, 15), 3),
                       cases_reported = 1:9)
  out <- split_corrected_reports(df, time_corrected = 10, onset_times = c(1, 2))
  key <- function(d) sort(paste(col_of(d, "time_onset"), col_of(d, "time_reported")))
  expect_equal(key(out$uncorrected), c("1 5", "2 5"))
  expect_equal(key(out$corrected), sort(c("1 10", "1 15", "2 10", "2 15", "3 5", "3 10", "3 15")))
  # the other columns are carried along
  expect_equal(sort(c(col_of(out$corrected, "cases_reported"), col_of(out$uncorrected, "cases_reported"))), 1:9)
})

test_that("T9 thin_series keeps the likelihood of the binomial chain (A.10.2 telescoping)", {
  df <- tibble::tibble(time_onset = 1, time_reported = 1:12,
                       cases_reported = c(0, 0, 2, 2, 2, 5, 5, 5, 5, 6, 6, 6))
  th <- thin_series(df)
  tr <- col_of(th, "time_reported"); cr <- col_of(th, "cases_reported")
  # rows kept are rows of the input
  expect_true(all(paste(tr, cr) %in% paste(df$time_reported, df$cases_reported)))
  for (N in c(6, 10)) for (pars in list(c(4, 2), c(8, 5))) {
    surv <- function(d) ref_surv_case(d, pars[1], pars[2])
    full <- ref_chain_loglik(df$cases_reported, df$time_reported, 1, N, surv)
    o <- order(tr)
    thin <- ref_chain_loglik(cr[o], tr[o], 1, N, surv)
    expect_equal(thin, full, tolerance = 1e-10,
                 label = sprintf("T9 thinned log-likelihood (N=%d, mean=%g, sd=%g)", N, pars[1], pars[2]))
  }
})

# BDA3 split R-hat: each chain split into halves; with m half-chains of n draws,
# B = n var(half-chain means), W = mean within variances, var+ = (n-1)/n W + B/n,
# R-hat = sqrt(var+ / W). Even chain lengths so the split is unambiguous.
ref_split_rhat <- function(values, chain) {
  halves <- unlist(lapply(split(values, chain), function(x) {
    n <- length(x) / 2; list(x[1:n], x[(n + 1):(2 * n)])
  }), recursive = FALSE)
  n <- length(halves[[1]])
  means <- vapply(halves, mean, 0); vars <- vapply(halves, stats::var, 0)
  B <- n * stats::var(means); W <- mean(vars)
  sqrt(((n - 1) / n * W + B / n) / W)
}

test_that("T9 split_rhat matches the BDA3 split R-hat", {
  set.seed(9001)
  v <- c(stats::rnorm(200), stats::rnorm(200, 0.3), cumsum(stats::rnorm(200)) / 10)
  ch <- rep(1:3, each = 200)
  expect_equal(as.numeric(split_rhat(v, ch)), ref_split_rhat(v, ch), tolerance = 1e-8)
  v2 <- stats::rnorm(400); ch2 <- rep(1:2, each = 200)
  expect_equal(as.numeric(split_rhat(v2, ch2)), ref_split_rhat(v2, ch2), tolerance = 1e-8)
  expect_true(is.na(split_rhat(c(1, 2, 3), c(1, 2, 3))), label = "T9 split_rhat NA with one draw per chain")
})

test_that("T9 batch_candidates finds coincident falls [API-only]", {
  df <- tibble::tibble(
    time_onset = rep(1:10, each = 2), time_reported = rep(c(5, 10), 10),
    cases_reported = c(5, 5,  5, 5,  6, 4,  6, 5,  5, 5,  7, 5,  5, 5,  5, 5,  4, 3,  5, 5))
  # report 10: falls at onsets 3 (2), 4 (1), 6 (2), 9 (1); onset 5 has no fall (gap 1)
  out <- batch_candidates(df, min_total_fall = 3, max_gap = 1)
  expect_equal(nrow(out), 1)
  expect_equal(as.numeric(col_of(out, "start")), 3)
  expect_equal(as.numeric(col_of(out, "width")), 4)
  expect_equal(as.numeric(col_of(out, "time_detected")), 10)
  expect_equal(as.numeric(col_of(out, "total_fall")), 5)
})

test_that("T9 fit_methods lists each model's default method first", {
  fm <- fit_methods()
  expect_true(all(c("discrete", "relaxed") %in% names(fm)))
  expect_equal(fm$discrete[[1]], "gibbs")
  expect_equal(fm$relaxed[[1]], "nuts")
})
