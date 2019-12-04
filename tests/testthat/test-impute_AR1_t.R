context("Function \"impute_AR1_t()\"")
#library(testthat)

load("t_data.RData")

test_that("error control works",{
  expect_error(impute_AR1_t(y = median), "\"y\" must be coercible to a vector or matrix.")
  expect_error(impute_AR1_t(y = "Hongkong"), "\"y\" only allows numerical or NA values.")
  expect_error(impute_AR1_t(y = c(1, 2, NA, 3, 4, NA)), "Each time series in \"y\" must have at least five observations.")
  expect_error(impute_AR1_t(y = y_missing_numeric, n_samples = 1.5), "\"n_samples\" must be a positive integer.")
  expect_error(impute_AR1_t(y = y_missing_numeric, n_burn = -1), "\"n_burn\" must be a positive integer.")
  expect_error(impute_AR1_t(y = y_missing_numeric, n_thin = -1), "\"n_thin\" must be a positive integer.")
}) 


test_that("imputation works", {
  # time series without NA's at the head and tail
  imputation_t <- impute_AR1_t(y_missing_numeric[, 1])
  expect_true(!anyNA(imputation_t))
  
  # time series with NA's at the head and tail
  imputation_t <- impute_AR1_t(y_missing_numeric[, 2], impute_leading_NAs = TRUE, impute_trailing_NAs = TRUE)
  expect_true(!anyNA(imputation_t))
})





