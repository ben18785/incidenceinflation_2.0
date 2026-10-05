f_batch_parameters <- function(beta=0.05) {
  list(beta=beta, zeta=3, phi_mu=0.4, phi_sigma=0.2, rho_batch=0.1)
}

test_that("footprints and overlaps are identified", {
  onset_times <- 1:20
  expect_equal(footprint_onset_times(5, 3, onset_times), 5:7)
  expect_equal(footprint_onset_times(19, 5, onset_times), 19:20)
  expect_equal(length(footprint_onset_times(25, 2, onset_times)), 0)

  batches <- dplyr::tibble(start=c(10, 3), width=c(4, 5), phi=0.3,
                           time_detected=c(30, 20))
  expect_false(is_overlapping_batches(batches))
  expect_false(is_overlapping_batches(batches[1, ]))
  expect_false(is_overlapping_batches(batches[0, ]))
  batches$width[2] <- 8
  expect_true(is_overlapping_batches(batches))

  expect_null(batch_for_onset(9, batches[1, ]))
  expect_equal(batch_for_onset(12, batches[1, ]),
               list(phi=0.3, time_detected=30))
  expect_null(batch_for_onset(12, NULL))
  expect_error(batch_for_onset(10, batches))
})

test_that("batches_logp_prior matches direct calculation", {
  batch_parameters <- f_batch_parameters()
  onset_times <- 1:20
  empty <- dplyr::tibble(start=numeric(0), width=numeric(0), phi=numeric(0),
                         time_detected=numeric(0))
  expect_equal(batches_logp_prior(empty, batch_parameters, onset_times),
               20 * log(0.95))

  batches <- dplyr::tibble(start=c(3, 18), width=c(4, 5), phi=c(0.3, 0.5),
                           time_detected=c(12, 19))
  # second footprint covers three onset times within the data
  logp <- 2 * log(0.05) + (20 - 4 - 3) * log(0.95) +
    sum(dpois(c(3, 4), 3, log=TRUE)) +
    sum(dgamma(c(0.3, 0.5), shape=4, rate=10, log=TRUE)) +
    sum(dgeom(c(8, 0), 0.1, log=TRUE))
  expect_equal(batches_logp_prior(batches, batch_parameters, onset_times), logp)

  # invalid batches
  f_invalid <- function(column, value) {
    batches_tmp <- batches
    batches_tmp[[column]][1] <- value
    batches_logp_prior(batches_tmp, batch_parameters, onset_times)
  }
  expect_equal(f_invalid("width", 0), -Inf)
  expect_equal(f_invalid("width", 16), -Inf) # overlaps
  expect_equal(f_invalid("start", 0), -Inf)
  expect_equal(f_invalid("time_detected", 3), -Inf)
  expect_equal(f_invalid("phi", 0), -Inf)
})

test_that("sample_batches draws batches which do not overlap", {
  set.seed(40)
  batch_parameters <- f_batch_parameters(0.2)
  l_batches <- purrr::map(1:200, ~sample_batches(40, batch_parameters))
  expect_false(any(purrr::map_lgl(l_batches, is_overlapping_batches)))
  expect_true(all(purrr::map_lgl(
    l_batches, ~all(.$time_detected >= .$start + 1) & all(.$width >= 1))))
  expect_true(all(is.finite(purrr::map_dbl(
    l_batches, ~batches_logp_prior(., batch_parameters, 1:40)))))
  widths <- unlist(purrr::map(l_batches, "width"))
  expect_equal(mean(widths), 4, tolerance=0.1)
  phis <- unlist(purrr::map(l_batches, "phi"))
  expect_equal(mean(phis), 0.4, tolerance=0.1)
})

test_that("batch_candidates finds falls in reported cases which coincide", {
  df <- tidyr::expand_grid(time_onset=1:8, time_reported=1:12) %>%
    dplyr::filter(time_reported >= time_onset) %>%
    dplyr::mutate(cases_reported=pmin(10, 10 * (time_reported - time_onset)))
  expect_equal(nrow(batch_candidates(df, 2)), 0)

  # onset times 2, 3 and 5 fall together at time 9; onset time 7 falls at time 11
  f_fall <- function(df, onsets, time, size)
    dplyr::mutate(df, cases_reported=ifelse(
      time_onset %in% onsets & time_reported >= time,
      cases_reported - size, cases_reported))
  df <- df %>% f_fall(c(2, 3, 5), 9, 2) %>% f_fall(7, 11, 1)
  candidates <- batch_candidates(df, min_total_fall=4)
  expect_equal(nrow(candidates), 1)
  expect_equal(candidates$start, 2)
  expect_equal(candidates$width, 4)
  expect_equal(candidates$time_detected, 9)
  expect_equal(candidates$total_fall, 6)

  expect_equal(nrow(batch_candidates(df, min_total_fall=1)), 2)
  # without allowing a gap the run is split
  candidates <- batch_candidates(df, min_total_fall=2, max_gap=0)
  expect_equal(candidates$start, c(2, 5))
  expect_equal(candidates$width, c(2, 1))
})

test_that("initial_batches chooses candidates which do not overlap", {
  candidates <- dplyr::tibble(start=c(2, 4, 10), width=c(4, 3, 2),
                              time_detected=c(9, 12, 15),
                              total_fall=c(6, 9, 3))
  batches <- initial_batches(candidates, f_batch_parameters())
  expect_equal(batches$start, c(4, 10))
  expect_equal(batches$time_detected, c(12, 15))
  expect_equal(batches$phi, c(0.4, 0.4))
  expect_equal(nrow(initial_batches(candidates[0, ], f_batch_parameters())), 0)
})

test_that("sample_batch_parameters draws from conditional posteriors", {
  set.seed(41)
  batches <- dplyr::tibble(start=c(3, 18), width=c(4, 5), phi=c(0.3, 0.5),
                           time_detected=c(12, 19))
  batch_parameters <- f_batch_parameters()
  expect_equal(sample_batch_parameters(batches, batch_parameters, list(), 1:20),
               batch_parameters)

  batch_priors <- list(beta_alpha=1, beta_beta=9, zeta_shape=2, zeta_rate=1,
                       rho_batch_alpha=1, rho_batch_beta=1)
  draws <- purrr::map(1:2000, ~sample_batch_parameters(
    batches, batch_parameters, batch_priors, 1:20))
  # beta(1 + 2, 9 + 13); gamma(2 + 7, 1 + 2); beta(1 + 2, 1 + 8)
  expect_equal(mean(purrr::map_dbl(draws, "beta")), 3 / 25, tolerance=0.05)
  expect_equal(mean(purrr::map_dbl(draws, "zeta")), 9 / 3, tolerance=0.05)
  expect_equal(mean(purrr::map_dbl(draws, "rho_batch")), 3 / 12, tolerance=0.05)
  expect_equal(draws[[1]]$phi_mu, 0.4)
})

test_that("sample_batches_step targets the prior when there are no data", {
  set.seed(42)
  batch_parameters <- f_batch_parameters(0.08)
  onset_times <- 1:25
  # a single observation for each onset time carries no information
  l_observations <- stats::setNames(
    purrr::map(onset_times, ~data.frame(time_reported=., cases_reported=0)),
    onset_times)
  cases_by_onset <- stats::setNames(rep(5, 25), onset_times)
  l_reporting <- stats::setNames(purrr::map(onset_times, ~list(mean=3, sd=2)),
                                 onset_times)
  revision_parameters <- list(eta=0.1, rho=0.3)

  n_prior <- 4000
  l_prior <- purrr::map(1:n_prior, ~sample_batches(25, batch_parameters))
  f_summaries <- function(l_batches) {
    c(n_batches=mean(purrr::map_dbl(l_batches, nrow)),
      n_covered=mean(purrr::map_dbl(l_batches, ~sum(purrr::map2_dbl(
        .$start, .$width, function(s, w) length(footprint_onset_times(s, w, onset_times)))))),
      mean_phi=mean(unlist(purrr::map(l_batches, "phi"))),
      mean_wait=mean(unlist(purrr::map(l_batches, ~.$time_detected - .$start))))
  }
  summaries_prior <- f_summaries(l_prior)

  # with and without candidates suggested by data
  l_candidates <- list(
    dplyr::tibble(start=numeric(0), width=numeric(0), time_detected=numeric(0),
                  total_fall=numeric(0)),
    dplyr::tibble(start=c(4, 15), width=c(3, 5), time_detected=c(9, 30),
                  total_fall=c(5, 5)))
  for(candidates in l_candidates) {
    batches <- sample_batches(25, batch_parameters)
    n_sweeps <- 12000
    l_draws <- vector("list", n_sweeps)
    for(i in 1:n_sweeps) {
      batches <- sample_batches_step(batches, onset_times, l_observations,
                                     cases_by_onset, l_reporting,
                                     revision_parameters, 3, batch_parameters,
                                     candidates)
      l_draws[[i]] <- batches
    }
    expect_false(any(purrr::map_lgl(l_draws, is_overlapping_batches)))
    summaries <- f_summaries(l_draws[-(1:500)])
    expect_equal(summaries[["n_batches"]], summaries_prior[["n_batches"]], tolerance=0.1)
    expect_equal(summaries[["n_covered"]], summaries_prior[["n_covered"]], tolerance=0.12)
    expect_equal(summaries[["mean_phi"]], summaries_prior[["mean_phi"]], tolerance=0.1)
    expect_equal(summaries[["mean_wait"]], summaries_prior[["mean_wait"]], tolerance=0.12)
  }
})

f_batch_data <- function(seed=43, days_total=30) {
  set.seed(seed)
  revision_parameters <- list(eta=0.01, rho=0.2)
  batches <- dplyr::tibble(start=8, width=5, phi=0.4, time_detected=22)
  Rt_function <- stats::approxfun(1:days_total, rep(1.1, days_total))
  df <- generate_snapshots(days_total, Rt_function, list(mean=5, sd=3),
                           list(mean=4, sd=2),
                           revision_parameters=revision_parameters,
                           batches=batches) %>%
    dplyr::mutate(reporting_piece_index=1)
  list(df=df, revision_parameters=revision_parameters, batches=batches,
       r_params=dplyr::tibble(reporting_piece_index=1, mean=4, sd=2))
}

test_that("sample_batches_step finds a batch which has been detected", {
  l <- f_batch_data()
  max_erroneous <- 3
  batch_parameters <- f_batch_parameters()
  candidates <- batch_candidates(l$df, max_erroneous + 1)
  expect_true(any(candidates$time_detected == 22 & candidates$start <= 9 &
                    candidates$start + candidates$width >= 12))

  lists <- observation_lists(l$df, l$r_params)
  revision_parameters <- l$revision_parameters
  revision_parameters$max_batch <- 40

  # without a batch the reported cases are impossible for some onset times
  logp_none <- observation_logp_by_onset(
    8:12, lists$l_observations, lists$cases_by_onset, lists$l_reporting,
    revision_parameters, max_erroneous)
  expect_true(any(!is.finite(logp_none)))
  revision_parameters$batches <- l$batches
  logp_batch <- observation_logp_by_onset(
    8:12, lists$l_observations, lists$cases_by_onset, lists$l_reporting,
    revision_parameters, max_erroneous)
  expect_true(all(is.finite(logp_batch)))

  # starting from no batches
  batches <- l$batches[0, ]
  n_sweeps <- 150
  l_draws <- vector("list", n_sweeps)
  for(i in 1:n_sweeps) {
    batches <- sample_batches_step(batches, lists$onset_times,
                                   lists$l_observations, lists$cases_by_onset,
                                   lists$l_reporting, revision_parameters,
                                   max_erroneous, batch_parameters, candidates)
    l_draws[[i]] <- batches
  }
  l_draws <- l_draws[-(1:50)]
  # a batch detected at time 22 whose footprint includes the middle of the true one
  is_found <- purrr::map_lgl(l_draws, ~any(
    .$time_detected == 22 & .$start <= 10 & .$start + .$width > 10))
  expect_true(mean(is_found) > 0.95)
  starts <- purrr::map_dbl(l_draws, ~.$start[which(.$time_detected == 22)[1]])
  expect_true(abs(mean(starts, na.rm=TRUE) - 8) <= 1.5)
  phis <- purrr::map_dbl(l_draws, ~.$phi[which(.$time_detected == 22)[1]])
  expect_true(mean(phis, na.rm=TRUE) > 0.15 & mean(phis, na.rm=TRUE) < 0.9)
})

test_that("mcmc works when there are batches of erroneous records", {
  l <- f_batch_data(seed=44)
  niter <- 4
  df <- l$df
  initial_cases_true <- df %>%
    dplyr::select(time_onset, cases_true) %>%
    unique()
  snapshot_with_Rt_index_df <- df %>%
    dplyr::select(-cases_true, -reporting_piece_index) %>%
    dplyr::mutate(Rt_index=ifelse(time_onset <= 15, 1, 2))
  priors <- list(Rt=list(shape=1, rate=1),
                 reporting=list(mean_mu=5, mean_sigma=10, sd_mu=3, sd_sigma=5),
                 revisions=list(eta_mu=0.05, eta_sigma=0.05,
                                rho_alpha=1, rho_beta=1),
                 batches=c(f_batch_parameters(),
                           list(beta_alpha=1, beta_beta=19, max_batch=40)),
                 max_erroneous=3,
                 max_cases=2000)
  f_mcmc <- function(priors, initial_revision_parameters, ...) {
    mcmc(niterations=niter,
         snapshot_with_Rt_index_df,
         priors,
         list(mean=5, sd=3),
         initial_cases_true,
         list(mean=4, sd=2),
         dplyr::tibble(Rt_index=c(1, 2), Rt=c(1.1, 1.1)),
         print_to_screen=FALSE,
         initial_revision_parameters=initial_revision_parameters,
         ...)
  }

  # falls this large are impossible without batches
  priors_none <- priors
  priors_none$batches <- NULL
  expect_error(f_mcmc(priors_none, l$revision_parameters), "max_erroneous")

  res <- f_mcmc(priors, l$revision_parameters)
  expect_equal(nrow(res$cases), 30 * niter)
  expect_equal(sum(is.na(res$cases$cases_true)), 0)
  expect_equal(nrow(res$batch_parameters), niter)
  expect_true(all(res$batch_parameters$n_batches >= 1))
  expect_true(all(c("start", "width", "phi", "time_detected", "iteration", "chain")
                  %in% colnames(res$batches)))
  expect_true(all(res$batches$time_detected[res$batches$start <= 10] == 22))
  expect_true(length(unique(res$batch_parameters$beta)) > 1)
  expect_equal(res$batch_parameters$zeta, rep(3, niter))

  # multiple chains
  res <- f_mcmc(priors, l$revision_parameters, nchains=2)
  expect_equal(sort(unique(res$batches$chain)), c(1, 2))
  expect_equal(nrow(res$batch_parameters), 2 * niter)

  # errors
  priors_wrong <- priors
  priors_wrong$batches$phi_mu <- NULL
  expect_error(f_mcmc(priors_wrong, l$revision_parameters))
  expect_error(f_mcmc(priors, NULL))
})
