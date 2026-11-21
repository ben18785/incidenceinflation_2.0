# The validation tests (test-validation-*.R) check the package against the
# model as specified, with expected values from independent code (see
# helper-validation-spec.R). They are slow, so they run only when
# VALIDATION_TIER asks for them:
#   fast    T1-T10, the samplers and simulators one at a time (~15 min)
#   medium  also M1-M5, whole-sampler runs with exact posteriors (~15 min more)
#   slow    also S2 and S9, recovery of deaths-model fits (~2.5 h more each)
#   all     also S1, simulation-based calibration (~1.5 h more)
# Thresholds are fixed in advance (helper-validation-stats.R): a failure is a
# finding to investigate, not a tolerance to retune.

validation_tiers <- c("fast", "medium", "slow", "all")

skip_validation <- function(tier) {
  level <- match(Sys.getenv("VALIDATION_TIER"), validation_tiers, nomatch = 0)
  needed <- match(tier, c("fast", "medium", "slow", "sbc"))
  testthat::skip_if(level < needed,
                    sprintf("validation tier '%s' (set VALIDATION_TIER)", validation_tiers[needed]))
}
