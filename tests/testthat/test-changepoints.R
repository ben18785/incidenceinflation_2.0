# all partitions of n_days into segments of at least min_length days
f_all_partitions <- function(n_days, min_length=1) {
  l_partitions <- list()
  for(k in 0:(2^(n_days - 1) - 1)) {
    is_start <- c(TRUE, as.logical(bitwAnd(k, 2^(0:(n_days - 2)))))
    if(n_days == 1)
      is_start <- TRUE
    starts <- which(is_start)
    if(is_valid_partition(starts, n_days, min_length))
      l_partitions[[length(l_partitions) + 1]] <- starts
  }
  l_partitions
}

f_partition_logw <- function(starts, log_marginals, p_change) {
  n_days <- nrow(log_marginals)
  ends <- c(starts[-1] - 1, n_days)
  n_changes <- length(starts) - 1
  sum(log_marginals[cbind(starts, ends)]) + n_changes * log(p_change) +
    (n_days - 1 - n_changes) * log(1 - p_change)
}

f_random_marginals <- function(n_days) {
  m <- matrix(rnorm(n_days^2), n_days, n_days)
  m[lower.tri(m)] <- -Inf
  m
}

test_that("is_valid_partition checks starts and minimum lengths", {
  expect_true(is_valid_partition(1, 10))
  expect_true(is_valid_partition(1, 3, min_length=7))
  expect_true(is_valid_partition(c(1, 4, 8), 10, min_length=3))
  expect_false(is_valid_partition(c(1, 4, 9), 10, min_length=3))
  expect_false(is_valid_partition(c(1, 3, 8), 10, min_length=3))
  expect_false(is_valid_partition(c(2, 5), 10))
  expect_false(is_valid_partition(c(1, 5, 5), 10))
  expect_false(is_valid_partition(c(1, 11), 10))
  expect_false(is_valid_partition(c(), 10))
})

test_that("segment_indices and piecewise_values expand segments to days", {
  expect_equal(segment_indices(c(1, 4, 5), 7), c(1, 1, 1, 2, 3, 3, 3))
  expect_equal(segment_indices(1, 3), c(1, 1, 1))
  expect_equal(piecewise_values(c(1, 3), c(2, 0.5), 5), c(2, 2, 0.5, 0.5, 0.5))
  expect_error(piecewise_values(c(1, 3), 2, 5))
  expect_error(piecewise_values(c(2, 3), c(1, 2), 5))
})

test_that("n_valid_partitions counts partitions correctly", {
  for(min_length in 1:3) {
    l_partitions <- f_all_partitions(9, min_length)
    n_changes <- purrr::map_dbl(l_partitions, length) - 1
    for(k in 0:8)
      expect_equal(n_valid_partitions(k, 9, min_length), sum(n_changes == k))
  }
})

test_that("changepoints_logp_prior is a beta-binomial over partitions", {
  n_days <- 12
  alpha <- 1.5
  beta <- 7
  # summed over all partitions, probabilities are one
  log_p <- purrr::map_dbl(0:(n_days - 1), function(k) {
    log(choose(n_days - 1, k)) + changepoints_logp_prior(k, n_days, alpha, beta)
  })
  expect_equal(sum(exp(log_p)), 1)
  # matches integrating the probability of a change out numerically
  f <- function(p) p^3 * (1 - p)^(n_days - 1 - 3) * dbeta(p, alpha, beta)
  expect_equal(exp(changepoints_logp_prior(3, n_days, alpha, beta)),
               integrate(f, 0, 1)$value, tolerance=1e-6)
  # with beta equal to the number of days and alpha of one, no changes has
  # prior probability of about a half
  expect_equal(exp(changepoints_logp_prior(0, 100, 1, 100)), 0.5, tolerance=0.01)
})

test_that("sample_change_probability draws from its beta conditional", {
  set.seed(1)
  draws <- replicate(4000, sample_change_probability(3, 50, 1, 20))
  expect_true(all(draws > 0 & draws < 1))
  expect_equal(mean(draws), (1 + 3) / (1 + 20 + 49), tolerance=0.05)
})

test_that("partition_forward sums the weights of all partitions", {
  set.seed(2)
  n_days <- 9
  log_marginals <- f_random_marginals(n_days)
  for(min_length in c(1, 2, 4)) {
    l_partitions <- f_all_partitions(n_days, min_length)
    log_w <- purrr::map_dbl(l_partitions, f_partition_logw, log_marginals, 0.3)
    log_forward <- partition_forward(log_marginals, 0.3, min_length)
    expect_equal(log_forward[n_days + 1], matrixStats::logSumExp(log_w))
  }
  expect_equal(partition_forward(log_marginals, 0.3)[1], 0)
  # a single segment is always possible, however short
  log_forward <- partition_forward(log_marginals[1:3, 1:3], 0.3, min_length=7)
  expect_equal(log_forward[4], log_marginals[1, 3] + 2 * log(0.7))
})

test_that("sample_partition draws partitions with the right probabilities", {
  set.seed(3)
  n_days <- 7
  log_marginals <- f_random_marginals(n_days)
  for(min_length in c(1, 2)) {
    l_partitions <- f_all_partitions(n_days, min_length)
    log_w <- purrr::map_dbl(l_partitions, f_partition_logw, log_marginals, 0.4)
    probs <- exp(log_w - matrixStats::logSumExp(log_w))
    names_partitions <- purrr::map_chr(l_partitions, paste, collapse="-")
    n_draws <- 20000
    draws <- replicate(n_draws, paste(
      sample_partition(log_marginals, 0.4, min_length), collapse="-"))
    expect_true(all(draws %in% names_partitions))
    freqs <- as.numeric(table(factor(draws, levels=names_partitions))) / n_draws
    expect_true(max(abs(freqs - probs)) < 0.015)
  }
})

test_that("sample_partition respects the minimum length of segments", {
  set.seed(4)
  log_marginals <- f_random_marginals(30)
  for(i in 1:50) {
    starts <- sample_partition(log_marginals, 0.5, min_length=5)
    expect_true(is_valid_partition(starts, 30, 5))
  }
})

test_that("sample_changepoints draws from the prior", {
  set.seed(5)
  n_days <- 20
  # with no minimum length the number of changes is binomial
  n_changes <- replicate(4000, length(sample_changepoints(n_days, 0.15)) - 1)
  expect_equal(mean(n_changes), 19 * 0.15, tolerance=0.05)
  expect_equal(var(n_changes), 19 * 0.15 * 0.85, tolerance=0.1)
  # with a minimum length its distribution follows from counting partitions
  n_changes <- replicate(6000, length(sample_changepoints(n_days, 0.15, 4)) - 1)
  k <- 0:5
  probs <- purrr::map_dbl(k, ~n_valid_partitions(., n_days, 4)) *
    0.15^k * 0.85^(19 - k)
  probs <- probs / sum(probs)
  freqs <- as.numeric(table(factor(n_changes, levels=k))) / 6000
  expect_true(max(abs(freqs - probs)) < 0.02)
  expect_equal(sample_changepoints(1, 0.5), 1)
})

test_that("alternating draws of partition and change probability give a beta-binomial", {
  set.seed(6)
  n_days <- 15
  alpha <- 2
  beta <- 6
  p_change <- 0.2
  n_changes <- vector(length=6000)
  for(i in seq_along(n_changes)) {
    n_changes[i] <- length(sample_changepoints(n_days, p_change)) - 1
    p_change <- sample_change_probability(n_changes[i], n_days, alpha, beta)
  }
  k <- 0:(n_days - 1)
  probs <- choose(n_days - 1, k) *
    exp(changepoints_logp_prior(k, n_days, alpha, beta))
  freqs <- as.numeric(table(factor(n_changes, levels=k))) / length(n_changes)
  expect_true(max(abs(freqs - probs)) < 0.03)
})

test_that("renewal_statistics gives infectiousness from earlier cases", {
  serial_parameters <- list(mean=4, sd=2)
  cases_true <- c(3, 0, 5, 2, 7, 1)
  statistics <- renewal_statistics(cases_true, serial_parameters, serial_max=4)
  w <- weights_series(4, serial_parameters)
  expect_equal(statistics$cases_true, cases_true)
  expect_equal(statistics$infectiousness[1], 0)
  expect_equal(statistics$infectiousness[2], w[1] * 3)
  expect_equal(statistics$infectiousness[4], w[1] * 5 + w[2] * 0 + w[3] * 3)
  expect_equal(statistics$infectiousness[6],
               w[1] * 7 + w[2] * 2 + w[3] * 5 + w[4] * 0)
  expect_equal(statistics$is_informative, c(FALSE, rep(TRUE, 5)))
  # agrees with the mean used in the state process
  expect_equal(
    2 * statistics$infectiousness[5],
    expected_cases(2, w, rev(cases_true[1:4])))
  # days following only zeros are uninformative
  statistics <- renewal_statistics(c(0, 0, 2, 1), serial_parameters, serial_max=4)
  expect_equal(statistics$is_informative, c(FALSE, FALSE, FALSE, TRUE))
})

test_that("segment_sums sums over informative days of each segment", {
  statistics <- dplyr::tibble(
    cases_true=c(4, 2, 3, 6),
    infectiousness=c(0, 1.5, 2.5, 3),
    is_informative=c(FALSE, TRUE, TRUE, TRUE))
  sums <- segment_sums(statistics)
  expect_equal(sums$cases[1, 4], 11)
  expect_equal(sums$cases[2, 3], 5)
  expect_equal(sums$cases[1, 1], 0)
  expect_equal(sums$infectiousness[1, 4], 7)
  expect_equal(sums$infectiousness[3, 4], 5.5)
})

test_that("segment_log_marginals_poisson matches numerical integration", {
  serial_parameters <- list(mean=4, sd=2)
  cases_true <- c(5, 6, 4, 9, 12, 10, 15)
  prior <- list(shape=2, rate=1.5)
  statistics <- renewal_statistics(cases_true, serial_parameters, serial_max=10)
  log_marginals <- segment_log_marginals_poisson(statistics, prior)
  expect_true(all(log_marginals[lower.tri(log_marginals)] == -Inf))
  f_check <- function(s, t) {
    days <- s:t
    days <- days[statistics$is_informative[days]]
    f <- function(Rt) {
      purrr::map_dbl(Rt, function(r) {
        exp(sum(dpois(cases_true[days], r * statistics$infectiousness[days],
                      log=TRUE))) * dgamma(r, prior$shape, prior$rate)
      })
    }
    # terms left out of the closed form as they do not depend on the partition
    log_const <- sum(cases_true[days] * log(statistics$infectiousness[days]) -
                       lgamma(cases_true[days] + 1))
    expect_equal(log_marginals[s, t] + log_const,
                 log(integrate(f, 0, Inf, rel.tol=1e-10)$value), tolerance=1e-5)
  }
  f_check(2, 5)
  f_check(1, 7)
  f_check(4, 4)
  f_check(6, 7)
  # a segment of only uninformative days has marginal likelihood of one
  expect_equal(log_marginals[1, 1], 0)
})

test_that("Rt_grid covers its range with adjoining intervals", {
  grid <- Rt_grid(n=50, lower=0.1, upper=10)
  expect_equal(nrow(grid), 50)
  expect_equal(range(grid$Rt), c(0.1, 10))
  expect_equal(grid$upper[-50], grid$lower[-1])
  expect_true(all(grid$lower < grid$Rt & grid$Rt < grid$upper))
  # the prior probabilities of the intervals sum to (almost) one
  log_prior <- Rt_grid_log_prior(Rt_grid(), list(shape=1, rate=0.5))
  expect_equal(sum(exp(log_prior)), 1, tolerance=0.02)
})

test_that("renewal_logl_grid gives the likelihood of each day at each Rt", {
  serial_parameters <- list(mean=4, sd=2)
  cases_true <- c(5, 6, 4, 9)
  statistics <- renewal_statistics(cases_true, serial_parameters, serial_max=10)
  grid <- Rt_grid(n=20)
  logl <- renewal_logl_grid(statistics, grid, kappa=3)
  expect_equal(dim(logl), c(4, 20))
  expect_true(all(logl[1, ] == 0))
  expect_equal(logl[3, 7], dnbinom(4, mu=grid$Rt[7] * statistics$infectiousness[3],
                                   size=3, log=TRUE))
  logl_poisson <- renewal_logl_grid(statistics, grid, is_negative_binomial=FALSE)
  expect_equal(logl_poisson[4, 12],
               dpois(9, grid$Rt[12] * statistics$infectiousness[4], log=TRUE))
  # the negative binomial tends to the Poisson as overdispersion vanishes
  expect_equal(renewal_logl_grid(statistics, grid, kappa=1e8), logl_poisson,
               tolerance=1e-5)
})

test_that("segment_log_marginals_grid matches the closed form for a Poisson model", {
  serial_parameters <- list(mean=4, sd=2)
  cases_true <- c(5, 6, 4, 9, 12, 10, 15, 22, 18)
  prior <- list(shape=2, rate=1.5)
  statistics <- renewal_statistics(cases_true, serial_parameters, serial_max=10)
  grid <- Rt_grid(n=2000, lower=0.001, upper=30)
  logl <- renewal_logl_grid(statistics, grid, is_negative_binomial=FALSE)
  log_marginals_grid <- segment_log_marginals_grid(
    logl, Rt_grid_log_prior(grid, prior))
  log_marginals <- segment_log_marginals_poisson(statistics, prior)
  log_const <- ifelse(
    statistics$is_informative,
    cases_true * log(statistics$infectiousness) - lgamma(cases_true + 1), 0)
  for(s in 1:9)
    for(t in s:9)
      expect_equal(log_marginals_grid[s, t],
                   log_marginals[s, t] + sum(log_const[s:t]), tolerance=1e-3)
  expect_true(all(log_marginals_grid[lower.tri(log_marginals_grid)] == -Inf))
  # a single day
  expect_equal(dim(segment_log_marginals_grid(logl[1, , drop=FALSE],
                                              Rt_grid_log_prior(grid, prior))),
               c(1, 1))
})

test_that("segment_log_marginals_grid matches numerical integration for a negative binomial model", {
  serial_parameters <- list(mean=4, sd=2)
  cases_true <- c(5, 6, 4, 9, 12, 10, 15)
  prior <- list(shape=2, rate=1.5)
  kappa <- 4
  statistics <- renewal_statistics(cases_true, serial_parameters, serial_max=10)
  grid <- Rt_grid(n=2000, lower=0.001, upper=30)
  log_marginals <- segment_log_marginals_grid(
    renewal_logl_grid(statistics, grid, kappa), Rt_grid_log_prior(grid, prior))
  f <- function(Rt) {
    purrr::map_dbl(Rt, function(r) {
      exp(sum(dnbinom(cases_true[3:6], mu=r * statistics$infectiousness[3:6],
                      size=kappa, log=TRUE))) * dgamma(r, prior$shape, prior$rate)
    })
  }
  expect_equal(log_marginals[3, 6], log(integrate(f, 0, Inf, rel.tol=1e-10)$value),
               tolerance=1e-3)
})

test_that("sample_Rt_grid draws from the posterior for a segment", {
  set.seed(7)
  serial_parameters <- list(mean=4, sd=2)
  cases_true <- c(20, 25, 31, 40, 52, 60, 75)
  prior <- list(shape=1, rate=0.5)
  statistics <- renewal_statistics(cases_true, serial_parameters, serial_max=10)
  grid <- Rt_grid()
  logl <- renewal_logl_grid(statistics, grid, is_negative_binomial=FALSE)
  draws <- replicate(4000, sample_Rt_grid(colSums(logl), Rt_grid_log_prior(grid, prior),
                                          grid))
  shape <- prior$shape + sum(cases_true[-1])
  rate <- prior$rate + sum(statistics$infectiousness[-1])
  expect_equal(mean(draws), shape / rate, tolerance=0.01)
  expect_equal(sd(draws), sqrt(shape) / rate, tolerance=0.08)
})

f_simulate_cases <- function(Rt, serial_parameters, initial=30, kappa=NULL) {
  n_days <- length(Rt)
  w <- weights_series(40, serial_parameters)
  cases <- c(rep(0, 40), rep(NA, n_days))
  cases[41] <- initial
  for(t in 2:n_days) {
    mean_cases <- Rt[t] * sum(w * rev(cases[t:(t + 39)]))
    if(is.null(kappa))
      cases[t + 40] <- rpois(1, mean_cases)
    else
      cases[t + 40] <- rnbinom(1, mu=mean_cases, size=kappa)
  }
  dplyr::tibble(time_onset=seq_len(n_days), cases_true=cases[-(1:40)])
}

test_that("sample_Rt_changepoints returns a piecewise-constant Rt", {
  set.seed(8)
  serial_parameters <- list(mean=5, sd=3)
  df <- f_simulate_cases(rep(1.2, 40), serial_parameters)
  df_Rt <- sample_Rt_changepoints(df, list(shape=1, rate=0.5), serial_parameters,
                                  p_change=0.2, min_length=3)
  expect_equal(colnames(df_Rt), c("time_onset", "Rt_index", "Rt"))
  expect_equal(df_Rt$time_onset, 1:40)
  expect_equal(df_Rt$Rt_index[1], 1)
  expect_true(all(diff(df_Rt$Rt_index) %in% c(0, 1)))
  expect_true(all(table(df_Rt$Rt_index) >= 3) | max(df_Rt$Rt_index) == 1)
  # one value of Rt for each segment
  expect_equal(nrow(unique(df_Rt[, c("Rt_index", "Rt")])), max(df_Rt$Rt_index))
  expect_true(all(df_Rt$Rt > 0))
  # rows can be supplied in any order
  df_Rt <- sample_Rt_changepoints(df[40:1, ], list(shape=1, rate=0.5),
                                  serial_parameters, p_change=0.2)
  expect_equal(df_Rt$time_onset, 1:40)
})

test_that("sample_Rt_changepoints checks its inputs", {
  serial_parameters <- list(mean=5, sd=3)
  df <- dplyr::tibble(time_onset=1:10, cases_true=10:19)
  prior <- list(shape=1, rate=0.5)
  expect_error(sample_Rt_changepoints(df[-4, ], prior, serial_parameters, 0.1))
  expect_error(sample_Rt_changepoints(df, prior, serial_parameters, 0))
  expect_error(sample_Rt_changepoints(df, prior, serial_parameters, 1))
  expect_error(sample_Rt_changepoints(df, prior, serial_parameters, 0.1,
                                      is_negative_binomial=TRUE))
  expect_error(sample_Rt_changepoints(df, prior, serial_parameters, 0.1,
                                      is_negative_binomial=TRUE, kappa=-1))
})

test_that("sample_Rt_changepoints finds a change in Rt and where it happens", {
  set.seed(9)
  serial_parameters <- list(mean=5, sd=3)
  prior <- list(shape=1, rate=0.5)
  Rt <- piecewise_values(c(1, 31), c(1.5, 0.6), 60)
  df <- f_simulate_cases(Rt, serial_parameters, initial=20)
  draws <- purrr::map_dfr(1:200, function(i) {
    sample_Rt_changepoints(df, prior, serial_parameters, p_change=1 / 60,
                           min_length=7) %>%
      dplyr::mutate(iteration=i)
  })
  summary_df <- draws %>%
    dplyr::group_by(.data$time_onset) %>%
    dplyr::summarise(Rt=median(.data$Rt))
  expect_true(all(abs(summary_df$Rt[10:25] - 1.5) < 0.25))
  expect_true(all(abs(summary_df$Rt[36:60] - 0.6) < 0.15))
  n_changes <- draws %>%
    dplyr::group_by(.data$iteration) %>%
    dplyr::summarise(n=max(.data$Rt_index) - 1)
  expect_true(mean(n_changes$n >= 1) > 0.95)
  # changes are placed near the true one
  starts <- draws %>%
    dplyr::group_by(.data$iteration) %>%
    dplyr::filter(.data$Rt_index == 2) %>%
    dplyr::summarise(start=min(.data$time_onset))
  expect_true(abs(median(starts$start) - 31) <= 3)
})

test_that("sample_Rt_changepoints rarely finds changes when there are none", {
  set.seed(10)
  serial_parameters <- list(mean=5, sd=3)
  prior <- list(shape=1, rate=0.5)
  df <- f_simulate_cases(rep(1.1, 60), serial_parameters, initial=20)
  n_changes <- replicate(200, {
    max(sample_Rt_changepoints(df, prior, serial_parameters, p_change=1 / 60,
                               min_length=7)$Rt_index) - 1
  })
  expect_true(mean(n_changes == 0) > 0.8)
})

test_that("sample_Rt_changepoints agrees between closed form and grid", {
  set.seed(11)
  serial_parameters <- list(mean=5, sd=3)
  prior <- list(shape=1, rate=0.5)
  Rt <- piecewise_values(c(1, 21), c(1.4, 0.8), 40)
  df <- f_simulate_cases(Rt, serial_parameters, initial=15)
  f_draws <- function(use_grid, ...) {
    purrr::map_dfr(1:300, function(i) {
      sample_Rt_changepoints(df, prior, serial_parameters, p_change=0.05,
                             min_length=5, use_grid=use_grid, ...) %>%
        dplyr::mutate(iteration=i)
    }) %>%
      dplyr::group_by(.data$time_onset) %>%
      dplyr::summarise(Rt=mean(.data$Rt), n=mean(.data$Rt_index))
  }
  exact <- f_draws(FALSE)
  grid <- f_draws(TRUE)
  expect_true(max(abs(exact$Rt - grid$Rt)) < 0.08)
  expect_true(max(abs(exact$n - grid$n)) < 0.2)
  # a negative binomial model with little overdispersion gives the same
  nb <- f_draws(FALSE, is_negative_binomial=TRUE, kappa=1e6)
  expect_true(max(abs(exact$Rt - nb$Rt)) < 0.08)
})

test_that("sample_Rt_changepoints finds fewer changes when overdispersion is allowed for", {
  set.seed(12)
  serial_parameters <- list(mean=5, sd=3)
  prior <- list(shape=1, rate=0.5)
  df <- f_simulate_cases(rep(1.05, 80), serial_parameters, initial=40, kappa=3)
  f_n_changes <- function(...) {
    mean(purrr::map_dbl(1:100, function(i) {
      max(sample_Rt_changepoints(df, prior, serial_parameters, p_change=1 / 80,
                                 min_length=3, ...)$Rt_index) - 1
    }))
  }
  expect_true(f_n_changes(is_negative_binomial=TRUE, kappa=3) < f_n_changes())
})

f_changepoint_data <- function(seed, days_total=50, starts=c(1, 26),
                               values=c(1.5, 0.7), kappa=1000) {
  set.seed(seed)
  Rt <- piecewise_values(starts, values, days_total)
  Rt_function <- function(t) Rt[pmin(pmax(round(t), 1), days_total)]
  df <- generate_snapshots(days_total, Rt_function, list(mean=4, sd=2),
                           list(mean=3, sd=2), kappa=kappa,
                           initial_parameters=list(mean=30, length=5))
  list(df=df, Rt=Rt)
}

f_mcmc_changepoints <- function(l, niter, Rt_changepoints=list(), ...) {
  df <- l$df
  initial_cases_true <- df %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(cases_true=max(cases_reported))
  snapshot_with_Rt_index_df <- df %>%
    dplyr::select(time_onset, time_reported, cases_reported) %>%
    dplyr::mutate(Rt_index=1)
  priors <- list(Rt=list(shape=1, rate=1),
                 reporting=list(mean_mu=5, mean_sigma=10, sd_mu=3, sd_sigma=5),
                 overdispersion=list(mean=10, sd=10),
                 max_cases=3000)
  if(!is.null(Rt_changepoints))
    priors$Rt_changepoints <- Rt_changepoints
  mcmc_single(niterations=niter,
              snapshot_with_Rt_index_df,
              priors,
              list(mean=4, sd=2),
              initial_cases_true,
              list(mean=3, sd=2),
              dplyr::tibble(Rt_index=1, Rt=1),
              print_to_screen=FALSE,
              ...)
}

test_that("mcmc_single samples the times of changes in Rt", {
  l <- f_changepoint_data(20)
  niter <- 4
  set.seed(21)
  res <- f_mcmc_changepoints(l, niter)
  n_days <- max(l$df$time_onset)
  expect_equal(colnames(res$Rt), c("iteration", "time_onset", "Rt_index", "Rt"))
  expect_equal(nrow(res$Rt), niter * n_days)
  expect_equal(colnames(res$Rt_changepoints), c("iteration", "n_changes", "p_change"))
  expect_equal(nrow(res$Rt_changepoints), niter)
  expect_true(all(res$Rt_changepoints$p_change > 0 & res$Rt_changepoints$p_change < 1))
  # the number of changes recorded matches the segments
  n_changes <- res$Rt %>%
    dplyr::group_by(iteration) %>%
    dplyr::summarise(n=max(Rt_index) - 1)
  expect_equal(n_changes$n, res$Rt_changepoints$n_changes)
  # segments are at least a week long by default
  lengths <- res$Rt %>%
    dplyr::count(iteration, Rt_index)
  expect_true(all(lengths$n >= 7))
  expect_equal(nrow(res$cases), niter * n_days)

  # without the prior on changes the pieces given in the data are used
  res_fixed <- f_mcmc_changepoints(l, niter, Rt_changepoints=NULL)
  expect_equal(colnames(res_fixed$Rt), c("iteration", "Rt_index", "Rt"))
  expect_false("Rt_changepoints" %in% names(res_fixed))
})

test_that("mcmc_single with changes in Rt honours its settings", {
  l <- f_changepoint_data(22, days_total=30, starts=c(1, 16))
  # a minimum length exceeding the series allows no changes
  res <- f_mcmc_changepoints(l, 3, list(min_length=40))
  expect_true(all(res$Rt_changepoints$n_changes == 0))
  # a prior which insists on changes produces them
  res <- f_mcmc_changepoints(l, 3, list(alpha=500, beta=1, min_length=1))
  expect_true(all(res$Rt_changepoints$n_changes > 15))
  expect_error(f_mcmc_changepoints(l, 2, list(alpha=0)))
  expect_error(f_mcmc_changepoints(l, 2, list(min_length=0)))
  expect_error(f_mcmc_changepoints(l, 2, maximise=TRUE))
})

test_that("mcmc_single with changes in Rt runs under a negative binomial model", {
  l <- f_changepoint_data(23, days_total=30, starts=c(1, 16), kappa=5)
  res <- f_mcmc_changepoints(l, 3, list(grid=Rt_grid(n=100)),
                             is_negative_binomial=TRUE)
  expect_equal(nrow(res$Rt), 3 * 30)
  expect_equal(nrow(res$overdispersion), 3)
  expect_true(all(res$Rt$Rt > 0))
})

test_that("mcmc recovers a change in Rt and combines chains", {
  skip_on_cran()
  l <- f_changepoint_data(24, days_total=60, starts=c(1, 31), values=c(1.4, 0.7))
  set.seed(25)
  res <- f_mcmc_changepoints(l, 120)
  summary_df <- res$Rt %>%
    dplyr::filter(iteration > 40) %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(Rt=median(Rt))
  expect_true(all(abs(summary_df$Rt[12:25] - 1.4) < 0.3))
  expect_true(all(abs(summary_df$Rt[37:50] - 0.7) < 0.25))
  expect_true(mean(res$Rt_changepoints$n_changes[-(1:40)] >= 1) > 0.9)

  combined <- combine_chains(list(res, res))
  expect_equal(nrow(combined$Rt_changepoints), 240)
  expect_equal(sort(unique(combined$Rt_changepoints$chain)), c(1, 2))
  expect_true("time_onset" %in% colnames(combined$Rt))
})

test_that("the negative binomial Rt update uses the overdispersion parameter", {
  # with strong overdispersion Rt is less certain than under a Poisson model
  set.seed(26)
  serial_parameters <- list(mean=5, sd=3)
  df <- f_simulate_cases(rep(1.2, 30), serial_parameters, initial=40) %>%
    dplyr::mutate(Rt_index=1)
  prior <- list(shape=1, rate=1)
  draws_small <- sample_Rt(df, prior, serial_parameters, kappa=0.5, ndraws=400,
                           is_negative_binomial=TRUE, nresamples=2000)
  draws_large <- sample_Rt(df, prior, serial_parameters, kappa=1000, ndraws=400,
                           is_negative_binomial=TRUE, nresamples=2000)
  expect_true(sd(draws_small$Rt) > 2 * sd(draws_large$Rt))
})
