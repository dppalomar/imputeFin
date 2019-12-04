context("Function \"impute_AR1_Gaussian()\"")
#library(testthat)

load("Gaussian_data.RData")

test_that("error control works",{
  expect_error(impute_AR1_Gaussian(y = median), "\"y\" must be coercible to a vector or matrix.")
  expect_error(impute_AR1_Gaussian(y = "Hongkong"), "\"y\" only allows numerical or NA values.")
  expect_error(impute_AR1_Gaussian(y = c(1, 2, NA, 3, 4, NA)), "Each time series in \"y\" must have at least five observations.")
  expect_error(impute_AR1_Gaussian(y = y_missing_numeric, n_samples = 1.5), "\"n_samples\" must be a positive integer.")
}) 


test_that("imputation works", {
  # time series without NA's at the head and tail
  imputation_Gaussian <- impute_AR1_Gaussian(y_missing_numeric[, 1])
  expect_true(!anyNA(imputation_Gaussian))
  
  # time series with NA's at the head and tail
  imputation_Gaussian <- impute_AR1_Gaussian(y_missing_numeric[, 2], impute_leading_NAs = TRUE, impute_trailing_NAs = TRUE)
  expect_true(!anyNA(imputation_Gaussian))
})





