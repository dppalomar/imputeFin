context("Function \"imputeAR1t()\"")

load("t_data.RData")

test_that("error control works",{
  expect_error(imputeAR1t(y = median), "\"y\" must be coercible to a vector or matrix.")
  expect_error(imputeAR1t(y = "Hongkong"), "\"y\" only allows numerical or NA values.")
  expect_error(imputeAR1t(y = c(1, 2, NA, 3, 4, NA)), "Each time series in \"y\" must have at least five observations.")
  expect_error(imputeAR1t(y = y_missing_numeric, n_samples = 1.5), "\"n_samples\" must be a positive integer.")
  expect_error(imputeAR1t(y = y_missing_numeric, n_burn = -1), "\"n_burn\" must be a positive integer.")
  expect_error(imputeAR1t(y = y_missing_numeric, n_thin = -1), "\"n_thin\" must be a positive integer.")
}) 


test_that("imputation works", {
  
  # time series without NA's at the head and tail
  imputation_t <- imputeAR1t(y_missing_numeric[,1])
  expect_true(!anyNA(imputation_t))
  
  # time series with NA's at the head and tail
  imputation_t <- imputeAR1t(y_missing_numeric[,2], impute_head_NAs = TRUE, impute_tail_NAs = TRUE)
  expect_true(!anyNA(imputation_t))
  
  
})





