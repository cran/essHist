context("checkHistogram")

test_that("checkHistogram with weird inputs", {
  x = rnorm(10)
  expect_error(checkHistogram(hist(x), c(x,1)))
})