# The joint move of the fraction dated by death and the delay to a death
# being recorded

ridge_example <- function() {
  set.seed(41)
  n_days <- 30
  serial_parameters <- list(mean=5, sd=3)
  cases <- true_cases(n_days + 40, function(t) 1, kappa=1000, serial_parameters,
                      initial_parameters=list(mean=30, length=10))
  sim <- simulate_cases_deaths(cases[41:(n_days + 40)], seq_len(n_days + 10),
                               list(mean=6, sd=4), 0.5, list(mean=3, sd=2.5),
                               triggered=0.5,
                               onset_to_death=list(mean=9.3, sd=5.5),
                               dated_by_death=list(fraction=0.3))
  list(snapshots=sim$snapshots, cases_before=cases[1:40],
       serial_parameters=serial_parameters,
       initial_cases=sim$snapshots %>%
         dplyr::group_by(time_onset) %>%
         dplyr::summarise(cases_true=max(cases_reported) + 2))
}

ridge_fit <- function(example, niterations, ridge=TRUE, fraction=list(), ...) {
  snapshots <- example$snapshots
  mcmc_single(
    niterations=niterations,
    snapshot_with_Rt_index_df=snapshots %>%
      dplyr::select(time_onset, time_reported, cases_reported) %>%
      dplyr::mutate(Rt_index=1),
    priors=list(Rt=list(shape=1, rate=1),
                reporting=list(mean_mu=6, mean_sigma=5, sd_mu=4, sd_sigma=4),
                max_cases=300,
                overdispersion=list(mean=20, sd=20),
                Rt_changepoints=list(),
                reporting_changepoints=list(min_length=14)),
    serial_parameters=example$serial_parameters,
    initial_cases_true=example$initial_cases,
    initial_reporting_parameters=list(mean=6, sd=4),
    initial_Rt=dplyr::tibble(Rt_index=1, Rt=1),
    print_to_screen=FALSE,
    is_negative_binomial=TRUE,
    initial_overdispersion=10,
    cases_before=example$cases_before,
    deaths=list(
      snapshot=snapshots %>%
        dplyr::select(time_onset, time_reported, deaths_reported),
      settings=list(onset_to_death=list(mean=9.3, sd=5.5), fatal_delay=FALSE,
                    triggered=list(initial=0.5, n_steps=0),
                    dated_by_death=list(fraction=fraction, exact=TRUE,
                                        ridge=ridge))),
    ...)
}

# the state of the model of deaths after a few iterations, with what the
# joint move needs
ridge_state <- function() {
  set.seed(3)
  fit <- ridge_fit(ridge_example(), 4, return_state=TRUE)
  v <- fit$state$variables
  list(state=v$deaths_state, cases_true=v$deaths_cases_true,
       reporting_by_onset=v$deaths_running_df, revision=v$revision_current,
       settings=v$deaths_options,
       renewal=dated_renewal(rep(1, length(v$deaths_cases_true)), 10, TRUE,
                             list(mean=5, sd=3), 40, v$cases_before))
}

test_that("the settings of the joint move are checked", {
  expect_null(ridge_settings(NULL))
  expect_null(ridge_settings(FALSE))
  expect_equal(ridge_settings(TRUE)$max_relabel, 40)
  expect_equal(ridge_settings(list(slope=2))$slope, 2)
  expect_error(ridge_settings(list(step=2)), "Unknown settings")
  expect_error(ridge_settings(list(sd=0)), "positive sd")
  onset_to_death <- list(mean=9.3, sd=5.5)
  expect_equal(deaths_settings(30, list(
    onset_to_death=onset_to_death,
    dated_by_death=list(fraction=list(), exact=TRUE, ridge=TRUE)))$dated_by_death$ridge,
    ridge_settings(TRUE))
  expect_null(deaths_settings(30, list(
    onset_to_death=onset_to_death,
    dated_by_death=list(fraction=list(), exact=TRUE)))$dated_by_death$ridge)
  expect_error(deaths_settings(30, list(
    onset_to_death=onset_to_death,
    dated_by_death=list(fraction=list(), ridge=TRUE))), "exact form")
})

test_that("the changes of kind are chosen with the probabilities given", {
  x <- ridge_state()
  state <- x$state
  dated <- state$dated
  # over every change which can be chosen, the probabilities sum to one
  by_onset <- state$deaths_true - dated$out
  total_up <- 0
  for(day in which(by_onset > 0))
    for(lag in seq_len(min(dated$max_lag, day - 1)))
      total_up <- total_up + exp(ridge_relabel_logq(
        state, list(day=day, lag=lag, sign=1)))
  expect_equal(total_up, 1)
  cells <- which(dated$lags > 0, arr.ind=TRUE)
  expect_gt(nrow(cells), 0)
  total_down <- sum(purrr::map_dbl(seq_len(nrow(cells)), function(i)
    exp(ridge_relabel_logq(state, list(day=cells[i, 1] + cells[i, 2],
                                       lag=cells[i, 2], sign=-1)))))
  expect_equal(total_down, 1)
  # changes which cannot be made cannot be chosen
  expect_equal(ridge_relabel_logq(state, list(day=1, lag=1, sign=1)), -Inf)
  empty <- which(dated$lags == 0, arr.ind=TRUE)[1, ]
  expect_equal(ridge_relabel_logq(
    state, list(day=empty[1] + empty[2], lag=empty[2], sign=-1)), -Inf)
  # the changes chosen are those with positive probability, as often as that
  set.seed(1)
  chosen <- purrr::map(1:4000, ~ridge_relabel_choose(state, -1))
  keys <- purrr::map_chr(chosen, ~paste(.$day, .$lag))
  top <- names(sort(table(keys), decreasing=TRUE))[1]
  move <- chosen[[match(top, keys)]]
  expect_equal(mean(keys == top), exp(ridge_relabel_logq(state, move)),
               tolerance=0.1)
  expect_true(all(purrr::map_dbl(chosen[1:200], ~ridge_relabel_logq(state, .)) > -Inf))
  up <- purrr::map(1:200, ~ridge_relabel_choose(state, 1))
  expect_true(all(purrr::map_dbl(up, ~ridge_relabel_logq(state, .)) > -Inf))
})

test_that("the reverse of a joint move returns to the state, with the reverse ratio", {
  x <- ridge_state()
  set.seed(11)
  for(sign in c(1, -1)) {
    forward <- ridge_propose(x$state, x$cases_true, x$reporting_by_onset,
                             x$revision, sign, 6, 0.04)
    expect_false(is.null(forward))
    expect_equal(forward$state$dated$fraction - x$state$dated$fraction,
                 rep(sign * 6 / sum(x$state$deaths_true), length(x$cases_true)))
    expect_equal(sum(forward$state$dated$out) - sum(x$state$dated$out), sign * 6)
    expect_equal(log(forward$state$recording$parameters$mean /
                       x$state$recording$parameters$mean),
                 rep(x$state$dated$ridge$slope * sign * 6 / sum(x$state$deaths_true) +
                       0.04, nrow(x$state$recording$parameters)))
    # the cases recorded for each day are not altered
    expect_equal(forward$cases_true - forward$state$dated$out + forward$state$dated$into,
                 x$cases_true - x$state$dated$out + x$state$dated$into)
    undo <- purrr::map(rev(forward$moves), ~list(type="relabel", day=.$day, lag=.$lag,
                                                lag_to=NA, sign=-sign))
    back <- ridge_propose(forward$state, forward$cases_true, x$reporting_by_onset,
                          x$revision, -sign, 6, -0.04, undo)
    for(name in c("lags", "out", "into", "unseen", "inflow", "fraction"))
      expect_equal(back$state$dated[[name]], x$state$dated[[name]], info=name)
    expect_equal(back$state$deaths_true, x$state$deaths_true)
    expect_equal(back$cases_true, x$cases_true)
    expect_equal(back$state$recording$parameters, x$state$recording$parameters)
    expect_equal(back$state$recording$centre, x$state$recording$centre)
    expect_equal(back$state$survivals, x$state$survivals)
    expect_equal(back$log_q, -forward$log_q)
    expect_equal(back$log_jacobian, -forward$log_jacobian)
    # so the ratios of the move and of its reverse cancel
    logp <- ridge_logp(x$state, x$cases_true, x$renewal, x$settings)$logp
    logp_forward <- ridge_logp(forward$state, forward$cases_true, x$renewal,
                               x$settings)$logp
    expect_true(is.finite(logp) & is.finite(logp_forward))
    expect_equal((logp_forward - logp + forward$log_jacobian + forward$log_q) +
                   (logp - logp_forward + back$log_jacobian + back$log_q), 0)
  }
  # with none changing kind, only the delay is scaled
  shift <- ridge_propose(x$state, x$cases_true, x$reporting_by_onset, x$revision,
                         1, 0, 0.1)
  expect_equal(shift$state$dated$fraction, x$state$dated$fraction)
  expect_equal(shift$log_q, 0)
  expect_equal(shift$log_jacobian, 0.1 * (nrow(x$state$recording$parameters) + 1))
})

test_that("the log-probability of the joint move is that of the model with the split of the reports summed out", {
  x <- ridge_state()
  state <- x$state
  # the reports and the fatal cases recorded for each day whose onset was
  # before the first, as in dated_logp
  terms <- ridge_logp(state, x$cases_true, x$renewal, x$settings)
  reports <- dated_logp(state, x$cases_true, cohorts=integer(0))
  rest <- dated_logp(state, x$cases_true, days=integer(0))
  recording <- state$recording
  prior <- reporting_centre_logp_prior(recording$centre, recording$spread,
                                       x$settings$delay_prior,
                                       x$settings$delay_changepoints$spread_sd) +
    sum(reporting_segments_logp_prior(recording$parameters, recording$centre,
                                      recording$spread)) +
    stats::dbeta(state$dated$fraction[1], state$dated$fraction_prior$alpha,
                 state$dated$fraction_prior$beta, log=TRUE)
  expect_equal(terms$logp, prior + reports + rest +
                 dated_renewal_logp(x$cases_true, x$renewal))
  # a fraction outside (0, 1) is impossible
  state$dated$fraction[] <- 1.2
  expect_equal(ridge_logp(state, x$cases_true, x$renewal, x$settings)$logp, -Inf)
})

test_that("a fit with the joint move runs and can be taken up again", {
  example <- ridge_example()
  set.seed(5)
  fit <- ridge_fit(example, 15)
  accepted <- fit$deaths_parameters$ridge_accepted
  expect_length(accepted, 15)
  expect_gt(sum(accepted), 0)
  both <- dplyr::left_join(fit$deaths, fit$cases, by=c("iteration", "time_onset"))
  expect_true(all(both$deaths_true <= both$cases_true))
  expect_true(all(both$dated_by_death <= both$deaths_recorded))
  # without it there is no such column
  set.seed(5)
  expect_false("ridge_accepted" %in%
                 names(ridge_fit(example, 2, ridge=NULL)$deaths_parameters))
  # with a known fraction, the fraction stays where it is
  set.seed(5)
  known <- ridge_fit(example, 4, fraction=0.3)
  expect_equal(known$deaths_parameters$dated_fraction, rep(0.3, 4))
  # taken up again
  set.seed(7)
  whole <- ridge_fit(example, 6)
  set.seed(7)
  first <- ridge_fit(example, 2, return_state=TRUE)
  second <- ridge_fit(example, 4, resume_from=first$state)
  for(name in setdiff(names(whole), c("state", "iterations_done")))
    expect_equal(dplyr::bind_rows(first[[name]], second[[name]]), whole[[name]],
                 info=name)
})
