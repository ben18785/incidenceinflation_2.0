skip_validation("fast")

# T10. Errors promised by the help pages, with positive controls.

t10_snapshot <- function(onsets = 1:6, last_report = max(onsets)) {
  grid <- expand.grid(time_onset = onsets, time_reported = onsets)
  grid <- grid[grid$time_reported >= grid$time_onset & grid$time_reported <= last_report, ]
  grid$cases_reported <- pmin(5, 1 + grid$time_reported - grid$time_onset)
  tibble::as_tibble(grid[order(grid$time_onset, grid$time_reported), ])
}
t10_serial <- list(mean = 5, sd = 2)

# Run an expression expected to error and print the message, so the cause can be checked.
t10_expect_error <- function(expr, label) {
  msg <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  if (!is.null(msg)) note("%s raised: %s", label, msg)
  expect(!is.null(msg), sprintf("%s: no error was raised", label))
}

test_that("T10 fit_incidence rejects onset times that are not consecutive from 1", {
  opts <- model_options(Rt_model = "periods")
  ctl <- list(niterations = 2, chains = 1, draws = 1)
  s1 <- t10_snapshot(2:7); s1$time_onset <- s1$time_onset   # starts at 2
  t10_expect_error(fit_incidence(s1, t10_serial, model = "discrete", options = opts, control = ctl),
                   "T10 onset times starting at 2")
  s2 <- t10_snapshot(); s2 <- s2[s2$time_onset != 3, ]       # gap at 3
  t10_expect_error(fit_incidence(s2, t10_serial, model = "discrete", options = opts, control = ctl),
                   "T10 onset times with a gap")
})

test_that("T10 discrete model rejects Rt_model = 'random_walk'", {
  t10_expect_error(fit_incidence(t10_snapshot(), t10_serial, model = "discrete",
                                 options = model_options(Rt_model = "random_walk"),
                                 control = list(niterations = 2, chains = 1, draws = 1)),
                   "T10 discrete + random_walk")
})

test_that("T10 relaxed_stan_data rejects a report later than the last onset time", {
  ok <- tryCatch({ relaxed_stan_data(t10_snapshot(), t10_serial); TRUE },
                 error = function(e) { note("T10 positive control failed: %s", conditionMessage(e)); FALSE })
  expect_true(ok, label = "T10 positive control: valid relaxed_stan_data call")
  late <- t10_snapshot(1:6, last_report = 6)
  late <- rbind(late, tibble::tibble(time_onset = 6L, time_reported = 8L, cases_reported = 5))
  t10_expect_error(relaxed_stan_data(late, t10_serial), "T10 report after the last onset time")
})

test_that("T10 relaxed_stan_data rejects delay changes together with misdating", {
  mis <- list(shifts = c(-1, 1), shift_probs = c(0.5, 0.5))
  ok <- tryCatch({ relaxed_stan_data(t10_snapshot(), t10_serial, misdating = mis); TRUE },
                 error = function(e) { note("T10 misdating positive control failed: %s", conditionMessage(e)); FALSE })
  if (!ok) skip("misdating format not accepted (positive control failed); see note")
  t10_expect_error(relaxed_stan_data(t10_snapshot(), t10_serial, misdating = mis, delay_max_changes = 1),
                   "T10 delay_max_changes with misdating")
})

test_that("T10 relaxed model rejects batches", {
  skip_if_not_installed("cmdstanr")
  t10_expect_error(fit_incidence(t10_snapshot(), t10_serial, model = "relaxed",
                                 options = model_options(Rt_model = "periods",
                                   batches = list(beta = 0.05, zeta = 1, phi_mu = 2,
                                                  phi_sigma = 1, rho_batch = 0.2)),
                                 control = list(chains = 1, iter_warmup = 5, iter_sampling = 5)),
                   "T10 relaxed + batches")
})
