context("test-regDIF")

test_that("Warning for not dichotomous responses works", {
  expect_error(regDIF(mtcars[,1:4], mtcars[,6], 0), "not dichotomously scored")
})
