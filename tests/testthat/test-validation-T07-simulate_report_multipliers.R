skip_validation("fast")

# T7. simulate_report_multipliers. spec eq. (14): the factor of day s is that of the
# first day on or after s on which a report is published, m_s = psi_{w(r(s))}.
# API: each report's multiplier is its weekday effect times its own random effect,
# whose log is N(0, report_sd^2); weekday of day d = ((weekday_first - 1 + d - 1) mod 7) + 1.

t7_weekday <- function(d, first) ((first - 1 + d - 1) %% 7) + 1

test_that("T7 day multipliers come from the first report on or after the day (eq. 14)", {
  rd <- c(3, 5, 6, 10, 17, 18)
  eff <- stats::setNames(c(0.1, -0.2, 0.3, -0.4, 0.5, -0.6, 0.7), as.character(1:7))
  out <- simulate_report_multipliers(rd, weekday_log_effects = eff, report_sd = 0, weekday_first = 3)
  # with no random effects every report's log effect is zero
  expect_equal(as.numeric(out$report_log_effects), rep(0, length(rd)), tolerance = 1e-12)
  rm <- out$report_multipliers
  days <- rm$day_first + seq_along(rm$values) - 1
  keep <- days >= 1 & days <= max(rd)
  r_of <- vapply(days[keep], function(d) rd[rd >= d][1], 0)
  expected <- exp(eff[as.character(t7_weekday(r_of, 3))])
  expect(sum(keep) > 0, "no day multipliers returned for days up to the last report")
  expect_equal(as.numeric(rm$values[keep]), as.numeric(expected), tolerance = 1e-10,
               label = "T7 multipliers of days 1..last report")
})

test_that("T7 random report effects are N(0, sd^2) and set each day's multiplier", {
  set.seed(7007)
  rd <- 1:2000
  out <- simulate_report_multipliers(rd, weekday_log_effects = NULL, report_sd = 0.5, weekday_first = 1)
  le <- as.numeric(out$report_log_effects)
  expect_equal(length(le), length(rd))
  expect_ks(le, function(q) stats::pnorm(q, 0, 0.5), "T7 report log effects ~ N(0, 0.5^2)")
  rm <- out$report_multipliers
  days <- rm$day_first + seq_along(rm$values) - 1
  keep <- days >= 1 & days <= max(rd)
  expect_equal(as.numeric(rm$values[keep]), exp(le[days[keep]]), tolerance = 1e-10,
               label = "T7 day multipliers with a report every day")
})
