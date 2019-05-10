context("msQuantile")

test_that("msQuantile with weird inputs", {
  expect_error(msQuantile())
  expect_error(msQuantile(NA))
  expect_error(msQuantile(1))
  expect_error(msQuantile(numeric(0)))
  expect_error(msQuantile(c(3,5)))
  expect_error(msQuantile(100,alpha=10))
  expect_error(msQuantile(100,alpha=c(0.1,1.1)))
  expect_error(msQuantile(10,alpha=NA))
  expect_error(msQuantile(10,alpha=numeric(0)))
  expect_error(msQuantile(10,intv=numeric(0)))
  intv = list('left'=c(0,1,2,4),'right'=c(1,2,3))
  expect_error(msQuantile(10,intv=intv))
  intv = list('left'=c(0,1,2,4),'right'=c(0,1,2,3))
  expect_error(msQuantile(10,intv=intv))
})

test_that("the validity of msQuantile", {
  n = 1e3
  
  expect_equal(msQuantile(n), msQuantile(n))
  
  minThd <- function (n, intv) {
    theta = unique(intv$right - intv$left + (intv$left == 1))/n
    max(-sqrt(2 - 2*log (theta*(1-theta))))
  } 
  expect_less_than(minThd(n, genIntv(n)), msQuantile(n))
  
  expect_less_than(msQuantile(n, mode="Con"), msQuantile(n, mode="Gen"))
})
