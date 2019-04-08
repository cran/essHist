context("genIntv")

test_that("genIntv with weird inputs", {
  expect_error(genIntv(0))
  expect_error(genIntv(numeric(0)))
  expect_error(genIntv(Inf))
  expect_error(genIntv(NA))
  expect_error(genIntv(NaN))
  expect_error(genIntv(1))
  expect_error(genIntv(2, "Hello"))
})