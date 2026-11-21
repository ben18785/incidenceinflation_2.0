# B2 of the test plan: with the settings of the sensitivity analyses at their
# defaults, the main model is unchanged. The draws of short fits from fixed
# seeds, and the state at their end, are those recorded before the changes
# (fixtures/make-reference.R), so the new code also draws no random numbers
# when it is not used.

test_that("the main model gives the draws and state recorded before the changes", {
  path <- test_path("fixtures", "reference.rds")
  skip_if_not(file.exists(path), "no reference fits recorded")
  reference <- readRDS(path)
  skip_if_not(identical(reference$platform, reference_platform()),
              "the reference fits were recorded on another platform")
  skip_on_cran()
  cases <- reference$cases
  for(i in seq_len(nrow(cases))) {
    fit <- reference_fit(reference_data(cases$seed[i]), cases$ridge[i],
                         reference_iterations, seed=cases$seed[i] + 1)
    now <- reference_summary(fit)
    before <- reference$fits[[i]]
    label <- sprintf("data set %d, ridge %s", cases$seed[i], cases$ridge[i])
    for(name in names(before$draws))
      expect_identical(now$draws[[name]], before$draws[[name]],
                       label=paste(label, "draws", name))
    for(name in names(before$deaths_state))
      expect_identical(now$deaths_state[[name]], before$deaths_state[[name]],
                       label=paste(label, "state", name))
    expect_identical(now$cases_true, before$cases_true, label=label)
    expect_identical(now$reporting_logl, before$reporting_logl, label=label)
    expect_identical(now$log_likelihood, before$log_likelihood, label=label)
  }
})
