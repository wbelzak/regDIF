context("test-regDIF")

test_that("Stop for negative tuning values.", {
  expect_error(regDIF(ida[,1:6], ida[,7:9], itemtypes = "categorical", penalty = -1, anchor = 1), "Penalty values must be non-negative.")
})
