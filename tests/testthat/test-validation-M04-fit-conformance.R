skip_validation("medium")

# M4 (medium tier): output conformance of fit_incidence (discrete, gibbs).
# A short run on data from the independent simulator (helper-validation-sim.R); the checks are those of the
# fit_incidence help page (Value) and its control settings, plus basic invariants. No
# statistical expectation is made: 200 iterations are far too few for that.

test_that("M4 fit_incidence (discrete) returns the documented object", {
  skip_validation("medium")
  set.seed(7004)
  omega <- ref_serial_weights(15.3, 9.3, 40)
  hist <- rep(20, 40); n_days <- 30
  I <- sim_renewal(c(rep(1.2, 15), rep(0.9, 15)), hist, omega)
  d <- sim_reports(I, seq(2, n_days, 2), list(mean = 5, sd = 2.5))
  snap <- d$cases[d$cases$time_reported <= n_days, ]
  chains <- 2; draws <- 40
  t0 <- proc.time()[["elapsed"]]
  fit <- fit_incidence(snap, list(mean = 15.3, sd = 9.3), model = "discrete",
                       options = model_options(is_negative_binomial = FALSE, Rt_model = "changepoints"),
                       cases_before = hist,
                       control = list(niterations = 200, chains = chains, draws = draws, seed = 7104))
  note("M4: fit_incidence (2 chains x 200 iterations) took %.1f s", proc.time()[["elapsed"]] - t0)

  expect_s3_class(fit, "incidenceinflation_fit")
  for (nm in c("cases", "Rt", "delay", "parameters", "diagnostics", "model", "method", "options",
               "priors", "seconds"))
    expect(nm %in% names(fit), sprintf("M4: element '%s' missing", nm))
  expect_identical(fit$model, "discrete")
  expect_identical(fit$method, "gibbs")
  expect(all(c("chain", "draw", "time_onset", "cases_true") %in% names(fit$cases)), "M4: cases columns")
  expect(all(c("chain", "draw", "time_onset", "Rt") %in% names(fit$Rt)), "M4: Rt columns")
  expect(all(c("chain", "draw", "time_onset", "mean", "sd") %in% names(fit$delay)), "M4: delay columns")
  expect(all(c("chain", "draw") %in% names(fit$parameters)), "M4: parameters columns")

  # draws per chain and one row per onset time in each draw
  for (el in c("cases", "Rt", "delay")) {
    x <- fit[[el]]
    per_chain <- tapply(col_of(x, "draw"), col_of(x, "chain"), function(v) length(unique(v)))
    expect(length(per_chain) == chains && all(per_chain == draws),
           sprintf("M4: %s has draws per chain %s (expected %d in each of %d chains)", el,
                   paste(per_chain, collapse = ", "), draws, chains))
    cnt <- table(paste(col_of(x, "chain"), col_of(x, "draw")))
    expect(all(cnt == n_days), sprintf("M4: %s does not have one row per onset time in every draw", el))
    days <- tapply(col_of(x, "time_onset"), paste(col_of(x, "chain"), col_of(x, "draw")),
                   function(v) identical(sort(as.integer(v)), seq_len(n_days)))
    expect(all(days), sprintf("M4: %s onset times are not 1..%d in every draw", el, n_days))
  }
  expect_equal(nrow(fit$parameters), chains * draws)

  # invariants
  ct <- col_of(fit$cases, "cases_true")
  expect(all(ct >= 0 & ct == round(ct)), "M4: cases_true not whole non-negative numbers")
  expect(all(col_of(fit$Rt, "Rt") > 0), "M4: Rt not positive")
  expect(all(col_of(fit$delay, "mean") > 0 & col_of(fit$delay, "sd") > 0), "M4: delay mean/sd not positive")
  expect("Rt_n_changes" %in% names(fit$parameters),
         sprintf("M4: Rt_model = 'changepoints' but no Rt_n_changes; parameters has: %s",
                 paste(names(fit$parameters), collapse = ", ")))
  if ("Rt_n_changes" %in% names(fit$parameters)) {
    k <- fit$parameters$Rt_n_changes
    expect(all(k >= 0 & k == round(k)), "M4: Rt_n_changes not a whole non-negative number")
    # R_t is piecewise constant with Rt_n_changes changes in each draw
    rt <- fit$Rt[order(fit$Rt$chain, fit$Rt$draw, fit$Rt$time_onset), ]
    nch <- tapply(rt$Rt, paste(rt$chain, rt$draw), function(v) sum(abs(diff(v)) > 0))
    key <- paste(fit$parameters$chain, fit$parameters$draw)
    expect(all(nch[key] == k), sprintf("M4: changes in R_t differ from Rt_n_changes in %d draws",
                                       sum(nch[key] != k)))
  }
  note("M4: parameters columns: %s; diagnostics: %s", paste(names(fit$parameters), collapse = ", "),
       paste(names(fit$diagnostics), collapse = ", "))
})
