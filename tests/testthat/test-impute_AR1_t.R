context("Function \"impute_AR1_t()\"")
#library(testthat)

data(ts_AR1_t)
y_missing         <- ts_AR1_t$y_missing
y_missing_numeric <- as.matrix(y_missing)


test_that("error control works",{
  expect_error(impute_AR1_t(y = median), "\"y\" must be coercible to a vector or matrix.")
  expect_error(impute_AR1_t(y = "Hongkong"), "\"y\" only allows numerical or NA values.")
  expect_error(impute_AR1_t(y = c(1, 2, NA, 3, 4, NA)), "Each time series in \"y\" must have at least 5 observations.")
  expect_error(impute_AR1_t(y = y_missing_numeric, n_samples = 1.5), "\"n_samples\" must be a positive integer.")
  expect_error(impute_AR1_t(y = y_missing_numeric, n_burn = -1), "\"n_burn\" must be a positive integer.")
  expect_error(impute_AR1_t(y = y_missing_numeric, n_thin = -1), "\"n_thin\" must be a positive integer.")
}) 


test_that("imputation works", {
  # time series without NA's at the head and tail
  imputation_t <- impute_AR1_t(y_missing_numeric[, 1])
  expect_false(anyNA(imputation_t))

  # imputation does not screw up the numerical values
  y_imputed <- impute_AR1_t(y_missing[, 3])
  expect_equivalent(y_missing[!is.na(y_missing[, 3]), 3], y_imputed[!is.na(y_missing[, 3])])
  #plot_imputed(y_imputed)
})


test_that("imputation plus outliers work", {
  y_outlier <- y_missing[, 3, drop = FALSE]
  idx_outliers <- c(100, 222)
  val_outliers <- c(100,  50)
  y_outlier[idx_outliers] <- val_outliers
  
  idx_not_NA <- which(!is.na(y_missing[, 3]))
  idx_not_NA_or_outlier <- setdiff(idx_not_NA, idx_outliers)
  
  expect_equal(y_missing[idx_not_NA_or_outlier, 3], 
               y_outlier[idx_not_NA_or_outlier])
      
  y_clean <- impute_AR1_t(y_outlier, remove_outliers = TRUE, outlier_prob_th = 0.002)
  expect_equivalent(y_missing[idx_not_NA_or_outlier, 3], 
                    y_clean[idx_not_NA_or_outlier])  

  expect_equivalent(y_missing[idx_not_NA, 3], 
                    y_clean[idx_not_NA], tolerance = 0.2)

  #plot_imputed(y_missing[, 3])
  #plot_imputed(y_outlier)
  #plot_imputed(y_clean)
})

