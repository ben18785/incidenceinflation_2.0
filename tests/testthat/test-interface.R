f_interface_data <- function(seed=1, days_total=24) {
  set.seed(seed)
  df <- generate_snapshots(days_total, function(t) rep(1.1, length(t)),
                           list(mean=4, sd=2), list(mean=4, sd=2.5),
                           initial_parameters=list(mean=40, length=10),
                           revision_parameters=list(eta=0.03, rho=0.3))
  dplyr::select(df, -cases_true)
}

f_check_form <- function(fit, n_days, n_draws) {
  expect_s3_class(fit, "incidenceinflation_fit")
  expect_equal(colnames(fit$cases), c("chain", "draw", "time_onset", "cases_true"))
  expect_equal(colnames(fit$Rt), c("chain", "draw", "time_onset", "Rt"))
  expect_equal(colnames(fit$delay), c("chain", "draw", "time_onset", "mean", "sd"))
  expect_equal(nrow(fit$cases), n_days * n_draws)
  expect_equal(nrow(fit$Rt), n_days * n_draws)
  expect_equal(nrow(fit$delay), n_days * n_draws)
  expect_equal(nrow(fit$parameters), n_draws)
  expect_true(all(c("chain", "draw", "eta", "rho") %in% colnames(fit$parameters)))
  expect_true(all(fit$Rt$Rt > 0))
  expect_true(fit$seconds > 0)
}

test_that("model_options and model_priors check what they are given", {
  options <- model_options()
  expect_equal(options$Rt_model, "random_walk")
  expect_equal(options$delay_model, "single")
  expect_null(options$weekday_first)
  expect_error(model_options(Rt_model="nonsense"))
  expect_error(model_options(delay_model="periods"))
  expect_error(model_options(misdating=list(shifts=c(-1, 1))))
  expect_equal(model_options(delay_model="periods",
                             reporting_piece_index=rep(1:2, 5))$delay_model, "periods")
  priors <- model_priors(eta_mu=0.02, Rt_shape=2)
  expect_equal(priors$eta_mu, 0.02)
  expect_equal(priors$Rt_shape, 2)
  expect_equal(priors$rho_beta, relaxed_priors()$rho_beta)
  expect_error(model_priors(nonsense=1))
})

test_that("options which a version of the model does not have give an error", {
  expect_error(check_options_model(model_options(), "discrete"), "random walk")
  expect_null(check_options_model(model_options(Rt_model="periods"), "discrete"))
  expect_null(check_options_model(model_options(), "relaxed"))
  expect_error(check_options_model(model_options(batches=list(beta=0.01)), "relaxed"),
               "discrete model")
  misdating <- list(shifts=c(-1, 1), shift_probs=c(0.5, 0.5))
  expect_error(check_options_model(
    model_options(delay_model="changepoints", misdating=misdating), "relaxed"))
  expect_null(check_options_model(
    model_options(Rt_model="periods", delay_model="changepoints", misdating=misdating),
    "discrete"))
  df <- f_interface_data()
  expect_error(fit_incidence(df, list(mean=4, sd=2), "discrete"))
  expect_error(fit_incidence(df, list(mean=4, sd=2), "nonsense"))
  # each version of the model has its own methods
  expect_equal(fit_methods()$discrete, "gibbs")
  expect_equal(fit_methods()$relaxed[1], "nuts")
  periods <- model_options(Rt_model="periods")
  expect_error(fit_incidence(df, list(mean=4, sd=2), "discrete", "nuts", periods),
               "'gibbs'")
  expect_error(fit_incidence(df, list(mean=4, sd=2), "relaxed", "gibbs"), "'nuts'")
  expect_error(fit_incidence(df, list(mean=4, sd=2), "relaxed", c("nuts", "laplace")))
  expect_error(fit_incidence(df[, 1:2], list(mean=4, sd=2), "discrete",
                             model_options(Rt_model="periods")))
  expect_error(fit_incidence(dplyr::filter(df, time_onset != 5), list(mean=4, sd=2),
                             "discrete", options=model_options(Rt_model="periods")))
})

test_that("period_index aligns periods to the last onset time", {
  expect_equal(period_index(14, 7), rep(1:2, each=7))
  expect_equal(period_index(10, 4), c(1, 1, 2, 2, 2, 2, 3, 3, 3, 3))
  expect_equal(period_index(4, 7, index=c(1, 1, 2, 2)), c(1, 1, 2, 2))
  expect_error(period_index(4, 7, index=1:3))
})

test_that("split_rhat detects disagreement between and within chains", {
  set.seed(2)
  chain <- rep(1:3, each=200)
  expect_equal(split_rhat(rnorm(600), chain), 1, tolerance=0.05)
  expect_true(split_rhat(rnorm(600) + 3 * (chain == 2), chain) > 1.5)
  # a drift within a single chain
  expect_true(split_rhat(rnorm(200) + seq(0, 5, length.out=200), rep(1, 200)) > 1.5)
  expect_equal(split_rhat(rep(1, 100), rep(1:2, each=50)), 1)
  expect_true(is.na(split_rhat(1:2, c(1, 1))))
})

test_that("the discrete model is fitted through the common interface", {
  df <- f_interface_data()
  control <- list(niterations=20, chains=2, draws=5, n_parameter_steps=1)
  set.seed(3)
  fit <- fit_incidence(df, list(mean=4, sd=2), "discrete",
                       options=model_options(Rt_model="periods", Rt_period_width=8),
                       control=control)
  f_check_form(fit, 24, 10)
  expect_equal(c(fit$model, fit$method), c("discrete", "gibbs"))
  expect_equal(sort(unique(fit$cases$chain)), 1:2)
  expect_true("kappa" %in% colnames(fit$parameters))
  # Rt is constant within each of three periods
  Rt <- fit$Rt %>% dplyr::filter(chain == 1, draw == 1) %>% dplyr::arrange(time_onset)
  expect_equal(Rt$Rt[1], Rt$Rt[8])
  expect_equal(Rt$Rt[9], Rt$Rt[16])
  expect_equal(length(unique(Rt$Rt)), 3)
  # a single delay throughout
  delay <- fit$delay %>% dplyr::filter(chain == 1, draw == 1)
  expect_equal(length(unique(delay$mean)), 1)
  expect_equal(nrow(fit$diagnostics$rhat), 4 + ncol(fit$parameters) - 2)
  expect_null(fit$raw)

  s <- summary(fit)
  expect_equal(nrow(s$cases), 24)
  expect_equal(colnames(s$Rt), c("time_onset", "q05", "q50", "q95"))
  expect_true(all(s$cases$q05 <= s$cases$q95))
  expect_true("eta" %in% s$parameters$parameter)
  expect_output(print(fit), "discrete model by method 'gibbs'")
  expect_output(print(fit), "Largest split R-hat")

  # other options: changes at unknown times, effect of the day of the week,
  # Poisson renewal, misdating, known earlier cases
  fit <- fit_incidence(
    df, list(mean=4, sd=2), "discrete", "gibbs",
    model_options(Rt_model="changepoints", delay_model="changepoints",
                  weekday_first=2, is_negative_binomial=FALSE, delay_min_length=5,
                  misdating=list(shifts=c(-1, 1), shift_probs=c(0.5, 0.5)),
                  max_erroneous=8),
    cases_before=rep(40, 10), control=control)
  f_check_form(fit, 24, 10)
  expect_true(all(c("Rt_n_changes", "delay_n_changes", "weekday_1", "pi_misdated") %in%
                    colnames(fit$parameters)))
  expect_false("kappa" %in% colnames(fit$parameters))
  expect_error(fit_incidence(df, list(mean=4, sd=2), "discrete",
                             options=model_options(Rt_model="periods"),
                             control=list(nonsense=1)))
  # delay periods given in advance
  fit <- fit_incidence(
    df, list(mean=4, sd=2), "discrete",
    options=model_options(Rt_model="periods", delay_model="periods",
                  reporting_piece_index=rep(1:2, each=12), tail=TRUE),
    control=control)
  delay <- fit$delay %>% dplyr::filter(chain == 1, draw == 1) %>% dplyr::arrange(time_onset)
  expect_equal(delay$mean[1], delay$mean[12])
  expect_equal(delay$mean[13], delay$mean[24])
  expect_true(all(c("tail_weight", "tail_mean") %in% colnames(fit$parameters)))
})

test_that("the relaxed model is fitted through the common interface, with results in the same form", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_version(), error=function(e) NULL)),
          "CmdStan is not installed")
  df <- f_interface_data()
  model <- relaxed_model()
  options <- model_options(Rt_model="periods", Rt_period_width=8, weekday_first=2,
                           days_sampled=10)
  fit <- suppressMessages(fit_incidence(
    df, list(mean=4, sd=2), "relaxed", "laplace", options,
    control=list(model=model, draws=50, seed=4, show_messages=FALSE)))
  f_check_form(fit, 24, 50)
  expect_equal(c(fit$model, fit$method), c("relaxed", "laplace"))
  expect_false(is.null(fit$raw))
  expect_null(fit$diagnostics$rhat)
  expect_output(print(fit), "relaxed model by method 'laplace'")

  # the same options in the discrete model give the same parameters
  fit_mcmc <- fit_incidence(df, list(mean=4, sd=2), "discrete", options=options,
                            control=list(niterations=20, chains=1, draws=5,
                                         n_parameter_steps=1))
  expect_setequal(colnames(fit$parameters), colnames(fit_mcmc$parameters))
  expect_equal(names(summary(fit)), names(summary(fit_mcmc)))

  # Hamiltonian Monte Carlo, with a random walk for Rt
  fit <- suppressWarnings(suppressMessages(fit_incidence(
    df, list(mean=4, sd=2), "relaxed", options=model_options(days_sampled=10),
    control=list(model=model, chains=2, parallel_chains=2, iter_warmup=30,
                 iter_sampling=20, seed=5, show_messages=FALSE,
                 show_exceptions=FALSE))))
  f_check_form(fit, 24, 40)
  expect_equal(fit$method, "nuts")
  expect_equal(sort(unique(fit$Rt$chain)), 1:2)
  expect_equal(max(fit$Rt$draw), 20)
  expect_true("Rt_step" %in% colnames(fit$parameters))
  expect_false(is.null(fit$diagnostics$rhat))
  expect_true(is.numeric(fit$diagnostics$divergences))
  expect_error(fit_incidence(df, list(mean=4, sd=2), "relaxed", "laplace",
                             model_options(batches=list(beta=0.01))))
})
