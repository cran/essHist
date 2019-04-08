context("essHistogram")

test_that("essHistogram with weird inputs", {
  expect_message(essHistogram(0))
  expect_error(essHistogram(Inf))
  expect_error(essHistogram(NaN))
  expect_error(essHistogram(NA))
  expect_error(essHistogram(rnorm(10),alpha=1.1))
  expect_error(essHistogram(rnorm(10),q=Inf))
  expect_warning(essHistogram(rnorm(10),alpha=c(0.1,0.2)))
  expect_warning(essHistogram(rnorm(10),q=0.1,alpha=0.1))
})


test_that("the validity of essHistogram", {
  # check whether the histogram is valid or not
  is.valid.hist <- function (h, data) {
    is.valid = TRUE
    if (length(data) != sum(h$counts)) {
      message('Invalid: Number of data points does NOT match!')
      is.valid = FALSE
    }
    if (length(h$counts) != length(h$density) && length(h$counts) != length(h$mids) 
        && length(h$counts) != length(h$breaks) -1 ) {
      message("Invalid: Components of 'h' have incompatable lengths!")
      is.valid = FALSE
    }
    if (sum(data <= h$breaks[2] & data >= h$breaks[1]) != h$counts[1]) {
      message('Invalid: Counts do NOT match break points!')
      is.valid = FALSE
    }
    if (length(h$counts) > 1) {
      for (k in 2:length(h$counts)) {
        if (sum(data <= h$breaks[k+1] & data > h$breaks[k]) != h$counts[k]) {
          message('Invalid: Counts do NOT match break points!')
          is.valid = FALSE
          break
        }
      }
    }
    if (is.valid) {
      message('The histogram is valid!\n')
    }
    is.valid
  }
  
  # home made examples
  data = c(rep(1,5),rep(2,3), rep(3,3))
  expect_true(is.valid.hist(essHistogram(data),data))
  
  data = c(rep(-5,100), pmax(rnorm(100),-4))
  expect_true(is.valid.hist(essHistogram(data),data))

  data = c(rep(-3,300),rep(3,600))
  expect_true(is.valid.hist(essHistogram(data),data))
  
  data = c(rep(-3,300),rep(-6,500),pmax(rnorm(1e3),-2.5))
  expect_true(is.valid.hist(essHistogram(data),data))
  
  # Thomas Staut's example
  data = c(5,7,3,8,1,1)
  expect_true(is.valid.hist(essHistogram(data),data))
  
  data = c(5,8,8)
  expect_true(is.valid.hist(essHistogram(data),data))
})