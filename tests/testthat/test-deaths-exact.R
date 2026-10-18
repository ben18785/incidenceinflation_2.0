# the probability of the reports of one onset time when the fatal cases dated
# by onset and by death are followed separately, worked out in R
exact_paths_reference <- function(cases, deaths, fatal, dated, survivors,
                                  arguments, wrong) {
  n_reports <- length(cases)
  states <- expand.grid(r=0:max(cases), d=0:max(deaths),
                        alive=0:wrong$max_erroneous, dead=0:wrong$max_erroneous)
  states <- states[states$alive + states$dead <= wrong$max_erroneous, ]
  cases_all <- c(0, cases)
  deaths_all <- c(0, deaths)
  f_transition <- function(j, from, to) {
    a <- to$r - from$r
    a_dated <- to$d - from$d
    survivors_from <- cases_all[j] - from$r - from$d - from$alive - from$dead
    survivors_to <- cases_all[j + 1] - to$r - to$d - to$alive - to$dead
    b <- survivors_to - survivors_from
    # fatal cases dated by onset which are shown as dead
    deaths_from <- deaths_all[j] - from$dead - from$d
    deaths_to <- deaths_all[j + 1] - to$dead - to$d
    deaths_new <- deaths_to - deaths_from
    if(a < 0 | a_dated < 0 | b < 0 | survivors_from < 0 | to$r > fatal | to$d > dated |
       survivors_to > survivors | deaths_from < 0 | deaths_new < 0 |
       deaths_to > to$r | deaths_from > from$r)
      return(0)
    p_wrong <- 0
    for(removed_alive in 0:from$alive)
      for(removed_dead in 0:from$dead) {
        entered_alive <- to$alive - from$alive + removed_alive
        entered_dead <- to$dead - from$dead + removed_dead
        if(entered_alive < 0 | entered_dead < 0)
          next
        p_wrong <- p_wrong +
          stats::dbinom(removed_alive, from$alive, wrong$w[j]) *
          stats::dbinom(removed_dead, from$dead, wrong$w[j]) *
          stats::dpois(entered_alive, wrong$lambda[j] * (1 - wrong$delta)) *
          stats::dpois(entered_dead, wrong$lambda[j] * wrong$delta)
      }
    x <- 0:deaths_new
    p_wrong *
      stats::dbinom(a, fatal - from$r, arguments$p_fatal[j]) *
      stats::dbinom(a_dated, dated - from$d, arguments$p_dated[j]) *
      stats::dbinom(b, survivors - survivors_from, arguments$p[j]) *
      sum(stats::dbinom(x, a, arguments$g[j]) *
            stats::dbinom(deaths_new - x, from$r - deaths_from, arguments$q[j]))
  }
  n_states <- nrow(states)
  alpha <- rep(0, n_states)
  alpha[states$r == 0 & states$d == 0 & states$alive == 0 & states$dead == 0] <- 1
  for(j in seq_len(n_reports)) {
    transition <- matrix(0, n_states, n_states)
    for(i in which(alpha > 0))
      for(k in seq_len(n_states))
        transition[i, k] <- f_transition(j, states[i, ], states[k, ])
    alpha <- as.vector(alpha %*% transition)
  }
  log(sum(alpha))
}

test_that("joint_path_exact_forward gives the probability of the reports", {
  days <- c(3, 4, 6, 8)
  cases <- c(2, 4, 5, 6)
  deaths <- c(1, 2, 2, 3)
  survivals <- joint_survivals(days, 2, list(mean=4, sd=3), list(mean=6, sd=4),
                               dated_reporting_parameters=list(mean=1.5, sd=1),
                               nu=-0.3)
  none <- list(lambda=rep(0, 4), w=rep(0, 4), delta=0, max_erroneous=0)
  for(fatal in 2:3)
    for(dated in 0:2) {
      forward <- joint_path_exact_forward(cases, deaths, fatal, dated, 3, survivals)
      expect_equal(forward$logp,
                   exact_paths_reference(cases, deaths, fatal, dated, 3,
                                         forward$arguments, none),
                   tolerance=1e-8)
    }
  # with none dated by death it is the model in which they are not followed
  expect_equal(
    joint_path_exact_forward(cases, deaths, 3, 0, 4, survivals)$logp,
    joint_path_forward(cases, deaths, 3, 4, path_probabilities(survivals, 0))$logp)
  # with wrong records, and a death count which falls
  deaths_fall <- c(1, 2, 1, 3)
  wrong <- list(lambda=c(0.3, 0.2, 0.1, 0.3), w=c(0, 0.3, 0.5, 0.4), delta=0.4,
                max_erroneous=1)
  forward <- joint_path_exact_forward(cases, deaths_fall, 3, 1, 3, survivals, wrong)
  expect_equal(forward$logp,
               exact_paths_reference(cases, deaths_fall, 3, 1, 3, forward$arguments,
                                     wrong),
               tolerance=1e-8)
  # with fatal cases reported when their death is recorded too, and none
  # dated by death, it is again the model in which they are not followed
  h <- lag_probabilities_gamma(4, 3)
  survivals_triggered <- joint_survivals(
    days, 2, list(mean=4, sd=3), list(mean=3, sd=2.5), NULL, "recording",
    list(before=list(lag_probabilities=h)), 0, 0.6)
  expect_equal(
    joint_path_exact_forward(cases, deaths, 3, 0, 4, survivals_triggered)$logp,
    joint_path_forward(cases, deaths, 3, 4, path_probabilities(survivals_triggered, 0))$logp)
  for(dated in 1:2) {
    forward <- joint_path_exact_forward(cases, deaths, 3, dated, 3, survivals_triggered)
    expect_equal(forward$logp,
                 exact_paths_reference(cases, deaths, 3, dated, 3, forward$arguments, none),
                 tolerance=1e-8)
  }
  expect_equal(joint_path_exact_forward(cases, deaths, -1, 1, 3, survivals)$logp, -Inf)
  expect_error(joint_path_exact_forward(c(1, 2), c(2, 2), 2, 1, 3, survivals))

  # a path drawn accounts for every count in the reports
  set.seed(41)
  forward <- joint_path_exact_forward(cases, deaths, 3, 2, 3, survivals)
  for(i in 1:20) {
    path <- joint_path_exact_backward(forward)
    expect_true(path$is_exact)
    expect_equal(cumsum(path$fatal_new + path$dated_new + path$survivors_new), cases)
    expect_equal(cumsum(path$deaths_new + path$dated_new), deaths)
    expect_true(all(path$alive_before >= 0) & all(path$dated_new >= 0))
    expect_true(sum(path$dated_new) <= 2 & sum(path$fatal_new) <= 3)
    # the probability of the path given the numbers of each kind is finite
    expect_true(is.finite(fatal_reports_logp(path, survivals, 5, 2)))
    expect_equal(fatal_reports_logp(path, survivals, sum(path$fatal_new) - 1 +
                                      2, 2), -Inf)
  }
})

test_that("mcmc_single runs with the exact form for fatal cases dated by death", {
  set.seed(42)
  n_days <- 30
  serial_parameters <- list(mean=5, sd=3)
  cases <- true_cases(n_days + 40, function(t) 1, kappa=1000, serial_parameters,
                      initial_parameters=list(mean=25, length=10))
  cases_before <- cases[1:40]
  cases <- cases[41:(n_days + 40)]
  dated_by_death <- list(fraction=0.4, mean=6, sd=3, reporting=list(mean=2, sd=1))
  sim <- simulate_cases_deaths(cases, seq_len(n_days + 30), list(mean=4, sd=3), 0.5,
                               list(mean=8, sd=4), dated_by_death=dated_by_death)
  snapshots <- sim$snapshots
  initial_cases <- snapshots %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(cases_true=max(cases_reported))
  settings <- deaths_settings(n_days, list(dated_by_death=c(dated_by_death,
                                                            list(exact=TRUE))))
  expect_true(settings$dated_by_death$exact)
  expect_equal(settings$dated_by_death$n_moves, n_days)
  fit <- mcmc_single(
    niterations=15,
    snapshot_with_Rt_index_df=snapshots %>%
      dplyr::select(time_onset, time_reported, cases_reported) %>%
      dplyr::mutate(Rt_index=1),
    priors=list(Rt=list(shape=1, rate=1),
                reporting=list(mean_mu=5, mean_sigma=5, sd_mu=3, sd_sigma=3),
                max_cases=300, Rt_changepoints=list(),
                reporting_changepoints=list(min_length=14)),
    serial_parameters=serial_parameters,
    initial_cases_true=initial_cases,
    initial_reporting_parameters=list(mean=4, sd=3),
    initial_Rt=dplyr::tibble(Rt_index=1, Rt=1),
    print_to_screen=FALSE,
    cases_before=cases_before,
    deaths=list(snapshot=snapshots %>%
                  dplyr::select(time_onset, time_reported, deaths_reported),
                settings=list(dated_by_death=c(dated_by_death, list(exact=TRUE)),
                              initial_delay=list(mean=8, sd=4))))
  expect_equal(nrow(fit$deaths), 15 * n_days)
  both <- dplyr::left_join(fit$deaths, fit$cases, by=c("iteration", "time_onset"))
  expect_true(all(both$deaths_true <= both$cases_true))
  expect_true(all(both$dated_by_death <= both$deaths_recorded))
  latest <- snapshots %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(deaths=max(deaths_reported))
  last <- dplyr::filter(fit$deaths, iteration == 15)
  expect_equal(last$deaths_recorded[1:10], latest$deaths[1:10])
  # the number dated by death changes between iterations
  totals <- fit$deaths %>%
    dplyr::group_by(iteration) %>%
    dplyr::summarise(dated=sum(dated_by_death))
  expect_true(length(unique(totals$dated)) > 3)
})

test_that("the exact form with the wrong records or the counts held sums to the free form", {
  days <- c(3, 4, 6, 8)
  cases <- c(2, 4, 5, 6)
  deaths <- c(1, 2, 1, 3)
  survivals <- joint_survivals(days, 2, list(mean=4, sd=3), list(mean=6, sd=4),
                               dated_reporting_parameters=list(mean=1.5, sd=1),
                               nu=-0.3)
  wrong <- list(lambda=c(0.3, 0.2, 0.1, 0.3), w=c(0, 0.3, 0.5, 0.4), delta=0.4,
                max_erroneous=1)
  f_sum <- function(logps) {
    logps <- logps[is.finite(logps)]
    max(logps) + log(sum(exp(logps - max(logps))))
  }
  free <- joint_path_exact_forward(cases, deaths, 3, 1, 3, survivals, wrong)$logp
  # over every path of the wrong records (three pairs each report)
  paths_wrong <- as.matrix(expand.grid(rep(list(0:2), 4)))
  logps <- apply(paths_wrong, 1, function(x)
    joint_path_exact_forward(cases, deaths, 3, 1, 3, survivals, wrong,
                             list(wrong=x))$logp)
  expect_equal(f_sum(logps), free, tolerance=1e-8)
  # over every path of the fatal cases of each kind reported
  paths_r <- as.matrix(expand.grid(rep(list(0:3), 4)))
  paths_r <- paths_r[apply(paths_r, 1, function(x) all(diff(x) >= 0)), ]
  paths_d <- as.matrix(expand.grid(rep(list(0:1), 4)))
  paths_d <- paths_d[apply(paths_d, 1, function(x) all(diff(x) >= 0)), ]
  logps <- unlist(lapply(seq_len(nrow(paths_r)), function(i)
    sapply(seq_len(nrow(paths_d)), function(j)
      joint_path_exact_forward(cases, deaths, 3, 1, 3, survivals, wrong,
                               list(r=paths_r[i, ], d=paths_d[j, ]))$logp)))
  expect_equal(f_sum(logps), free, tolerance=1e-8)

  # a path drawn in two parts keeps what it holds fixed consistent
  set.seed(3)
  forward <- joint_path_exact_forward(cases, deaths, 3, 1, 3, survivals, wrong)
  path <- joint_path_exact_backward(forward)
  for(i in 1:20) {
    forward <- joint_path_exact_forward(cases, deaths, 3, 1, 3, survivals, wrong,
                                        list(wrong=exact_fixed_wrong(path)))
    expect_true(is.finite(forward$logp))
    path <- joint_path_exact_redraw(forward)
    expect_equal(cumsum(path$fatal_new + path$dated_new + path$survivors_new) +
                   path$wrong_alive + path$wrong_dead, cases)
    expect_equal(cumsum(path$deaths_new + path$dated_new) + path$wrong_dead, deaths)
  }
})
