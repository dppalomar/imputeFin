context("Miscelaneous functions")
#library(testthat)


test_that("is_inner_NA() works",{
  y <- c(NA, NA, 1, 2, NA, NA, 3, 4, NA, NA, 5, 6, NA, NA)
  y_inner_NA <- c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
  expect_identical(is_inner_NA(y), y_inner_NA)
  expect_identical(is_inner_NA(as.matrix(y)), as.matrix(y_inner_NA))
  
  y_zoo <- zoo::zoo(y, seq(as.Date("2016-01-01"), length = 14, by = "days"))
  expect_identical(is_inner_NA(y_zoo), y_inner_NA)
  
  y_xts <- xts::as.xts(y_zoo)
  expect_identical(is_inner_NA(y_xts), as.matrix(y_inner_NA))  
}) 
