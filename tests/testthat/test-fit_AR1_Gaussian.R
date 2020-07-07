context("Function \"fit_AR1_Gaussian()\"")
#library(testthat)

data(ts_AR1_Gaussian)
y_missing         <- ts_AR1_Gaussian$y_missing
y_missing_numeric <- as.matrix(y_missing)


test_that("error control works",{
  expect_error(fit_AR1_Gaussian(y = median), "\"y\" must be coercible to a vector or matrix.")
  expect_error(fit_AR1_Gaussian(y = "Hongkong"), "\"y\" only allows numerical or NA values.")
  expect_error(fit_AR1_Gaussian(y = c(1, 2, NA, 3, 4, NA)), "Each time series in \"y\" must have at least 5 observations.")
  expect_error(fit_AR1_Gaussian(y = y_missing_numeric, tol = -0.1), "\"tol\" must be greater than 0.")
  expect_error(fit_AR1_Gaussian(y = y_missing_numeric, maxiter = -2), "\"maxiter\" must be greater than 1.")
}) 


test_that("default mode works", {
  # estimation_Gaussian_uni_check <- fit_AR1_Gaussian(y_missing_numeric[, 1])
  # estimation_Gaussian_multi_check <- fit_AR1_Gaussian(y_missing_numeric)
  # save(estimation_Gaussian_uni_check, estimation_Gaussian_multi_check, file = "estimation_Gaussian_check.RData", version = 2, compress = "xz")
  load("estimation_Gaussian_check.RData")
  
  # check the uni time series & numeric 
  estimation_Gaussian_uni <- fit_AR1_Gaussian(y_missing_numeric[, 1])
  expect_equal(estimation_Gaussian_uni, estimation_Gaussian_uni_check)
  
  # check the multiple time series & numeric  
  estimation_Gaussian_multi <- fit_AR1_Gaussian(y_missing_numeric)
  expect_equal(estimation_Gaussian_multi, estimation_Gaussian_multi_check)
  
  # check the uni time series & zoo 
  estimation_Gaussian_uni <- fit_AR1_Gaussian(y_missing[, 1])
  expect_equal(estimation_Gaussian_uni, estimation_Gaussian_uni_check)
  
  # echeck the multiple time series & zoo
  estimation_Gaussian_zoo <- fit_AR1_Gaussian(y_missing)
  expect_equal(estimation_Gaussian_zoo, estimation_Gaussian_multi_check)
})


test_that("time series with zero mean works", {
  # estimation_Gaussian_zero_mean_check <- fit_AR1_Gaussian(y_missing_numeric[, 1], zero_mean = TRUE)
  # save(estimation_Gaussian_zero_mean_check, file = "estimation_Gaussian_zero_mean_check.RData", version = 2, compress = "xz")
  load("estimation_Gaussian_zero_mean_check.RData")
  
  estimation_Gaussian_zero_mean <- fit_AR1_Gaussian(y_missing_numeric[, 1], zero_mean = TRUE)
  expect_equal(estimation_Gaussian_zero_mean, estimation_Gaussian_zero_mean_check)
})


test_that("time series following random walk works", {
  # estimation_Gaussian_random_walk_check <- fit_AR1_Gaussian(y_missing_numeric[, 1], random_walk = TRUE)
  # save(estimation_Gaussian_random_walk_check, file = "estimation_Gaussian_random_walk_check.RData", version = 2, compress = "xz")
  load("estimation_Gaussian_random_walk_check.RData")

  estimation_Gaussian_random_walk <- fit_AR1_Gaussian(y_missing_numeric[, 1], random_walk = TRUE)
  expect_equal(estimation_Gaussian_random_walk, estimation_Gaussian_random_walk_check)
})


test_that("remove_outliers = TRUE does not screw up", {
  fitted1 <- fit_AR1_Gaussian(y_missing[, 3, drop = FALSE], remove_outliers = FALSE)
  fitted2 <- fit_AR1_Gaussian(y_missing[, 3, drop = FALSE], remove_outliers = TRUE)
  expect_equal(fitted1, fitted2[1:3])
  
  fitted1 <- fit_AR1_Gaussian(y_missing, remove_outliers = FALSE)
  fitted2 <- fit_AR1_Gaussian(y_missing, remove_outliers = TRUE)
  expect_equal(fitted1[-c(1, 2, 3)], fitted2[-c(1, 2, 3)])
  expect_equal(fitted1[[1]], fitted1[[1]][1:3])
  expect_equal(fitted1[[2]], fitted1[[2]][1:3])
  expect_equal(fitted1[[3]], fitted1[[3]][1:3])
})


test_that("remove_outliers = TRUE detects outliers", {
  y_outlier <- y_missing[, 3, drop = FALSE]
  idx_outliers <- c(100, 222)
  val_outliers <- c(100,  50)
  y_outlier[idx_outliers] <- val_outliers
  
  fitted1 <- fit_AR1_t(y_outlier, remove_outliers = TRUE)  
  expect_equal(fitted1$index_outliers, idx_outliers)
  
  fitted2 <- fit_AR1_t(y_missing[, 3, drop = FALSE], remove_outliers = FALSE)
  expect_equal(fitted1$phi0,   fitted2$phi0,   tolerance = 0.001)
  expect_equal(fitted1$phi1,   fitted2$phi1,   tolerance = 0.001)
  expect_equal(fitted1$sigma2, fitted2$sigma2, tolerance = 0.003)
  expect_equal(fitted1$sigma2, fitted2$sigma2, tolerance = 0.05)
  
  # plot_imputed(impute_AR1_t(y_outlier, remove_outliers = FALSE))
  # plot_imputed(impute_AR1_t(y_outlier, remove_outliers = TRUE))
  
  # imputation_result1 <- impute_AR1_t(y_outlier, n_samples = 3, remove_outliers = FALSE)
  # imputation_result2 <- impute_AR1_t(y_outlier, n_samples = 3, remove_outliers = TRUE)
  # plot_imputed(imputation_result1$y_imputed.1)
  # plot_imputed(imputation_result2$y_imputed.1)
})


test_that("attributes work in multivariate case", {
  # create outlier
  y_outlier <- y_missing
  y_outlier[100, 1] <- 2*y_outlier[100, 1]
  # true indices of NAs and outliers
  idx_NAs <- apply(y_outlier, MARGIN = 2, FUN = function(x) which(is.na(as.numeric(x))))
  idx_outliers <- list(100L, NULL, NULL)
  names(idx_outliers) <- names(idx_NAs) <- colnames(y_outlier)
  idx_NAs[[1]] <- sort(union(idx_NAs[[1]], idx_outliers[[1]]))  # because the outlier will become an NA

  
  y_imputed <- impute_AR1_Gaussian(y_outlier, n_samples = 1, return_estimates = FALSE, remove_outliers = TRUE)
  expect_identical(attr(y_imputed, "index_miss"),     idx_NAs)
  expect_identical(attr(y_imputed, "index_outliers"), idx_outliers)
  
  res <- impute_AR1_Gaussian(y_outlier, n_samples = 1, return_estimates = TRUE, remove_outliers = TRUE)
  expect_identical(attr(res$y_imputed, "index_miss"),     idx_NAs)
  expect_identical(attr(res$y_imputed, "index_outliers"), idx_outliers)
  
  res <- impute_AR1_Gaussian(y_outlier, n_samples = 3, return_estimates = TRUE, remove_outliers = TRUE)
  expect_identical(attr(res$y_imputed.1, "index_miss"),     idx_NAs)
  expect_identical(attr(res$y_imputed.1, "index_outliers"), idx_outliers)
})

