skip_validation("slow")

# S3, dataset 1 of 3: recovery of the main model's misdating (the exact form with a single
# estimated fraction); data, fit and checks in helper-validation-S03.R.

test_that("S3 dataset 1: the main model's misdating recovers the CFR, fraction and delay (99% intervals)", {
  skip_validation("slow")
  s3_test(1)
})
