# Records the reference fits of the main model (helper-reference.R) at the
# current commit, for test-reference.R. Run from the package's directory, on
# the commit before a change to the model:
#   Rscript tests/testthat/fixtures/make-reference.R
devtools::load_all(".", quiet=TRUE)
library(dplyr)
source("tests/testthat/helper-reference.R")
cases <- reference_cases()
fits <- lapply(seq_len(nrow(cases)), function(i)
  reference_summary(reference_fit(reference_data(cases$seed[i]), cases$ridge[i],
                                  reference_iterations, seed=cases$seed[i] + 1)))
reference <- list(
  commit=system("git rev-parse HEAD", intern=TRUE),
  platform=reference_platform(), cases=cases, fits=fits)
saveRDS(reference, "tests/testthat/fixtures/reference.rds", compress="xz")
