context("Function \"fit_AR1_Gaussian()\"")
#library(testthat)

# # generate the data
# set.seed(123)
# phi0 <- 0
# phi1 <- 1
# sigma2 <- 0.01
# n <- 300
# n_miss <- 0.1*n
# n_drop <- 10
# n_total <- n + n_drop
# data <- vector(length = n_total)
# epsilon <- vector(length = n_total)  # innovations
# data[1] <- 0
# for (i in 2:n_total) {
#   epsilon[i] <- rnorm(1, 0, sqrt(sigma2))
#   data[i] <- phi0 + phi1 * data[i - 1] + epsilon[i]
# }
# data <- data[(n_drop + 1):n_total]  # drop the first n_drop to reduce the influence of initial point
# 
# m <- 3
# y_orig<- matrix(rep(data, m), nrow = n, ncol = m)
# colnames(y_orig) <- c("a", "b", "c")
# 
# y_missing_numeric <- y_orig
# index_miss1 <- sort(sample(2:(n - 1), n_miss))
# index_miss2 <- c(1, sort(sample(2:(n - 1), n_miss - 2)), n)
# y_missing_numeric[index_miss1, 1] <- NA
# y_missing_numeric[index_miss2, 2] <- NA
# y_missing <- zoo::zoo(y_missing_numeric, seq(as.Date("2016-01-01"), length = n, by = "days"))
# save(y_missing_numeric, y_missing, phi0, phi1, sigma2, file = "Gaussian_data.RData", version = 2, compress = "xz")

load("Gaussian_data.RData")


test_that("error control works",{
  expect_error(fit_AR1_Gaussian(y = median), "\"y\" must be coercible to a vector or matrix.")
  expect_error(fit_AR1_Gaussian(y = "Hongkong"), "\"y\" only allows numerical or NA values.")
  expect_error(fit_AR1_Gaussian(y = c(1, 2, NA, 3, 4, NA)), "Each time series in \"y\" must have at least five observations.")
  expect_error(fit_AR1_Gaussian(y = y_missing_numeric, tol = -0.1), "\"tol\" must be greater than 0.")
  expect_error(fit_AR1_Gaussian(y = y_missing_numeric, maxiter = -2), "\"maxiter\" must be greater than 1.")
}) 


test_that("default mode works", {
  # estimation_Gaussian_uni_check <- fit_AR1_Gaussian(y_missing_numeric[,1])
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
  # estimation_Gaussian_zero_mean_check <- fit_AR1_Gaussian(y_missing_numeric[,1], zero_mean = TRUE)
  # save(estimation_Gaussian_zero_mean_check, file = "estimation_Gaussian_zero_mean_check.RData", version = 2, compress = "xz")
  load("estimation_Gaussian_zero_mean_check.RData")
  
  estimation_Gaussian_zero_mean <- fit_AR1_Gaussian(y_missing_numeric[, 1], zero_mean = TRUE)
  expect_equal(estimation_Gaussian_zero_mean, estimation_Gaussian_zero_mean_check)
})


test_that("time series following random walk works", {
  # estimation_Gaussian_random_walk_check <- fit_AR1_Gaussian(y_missing_numeric[,1], random_walk = TRUE)
  # save(estimation_Gaussian_random_walk_check, file = "estimation_Gaussian_random_walk_check.RData", version = 2, compress = "xz")
  load("estimation_Gaussian_random_walk_check.RData")

  estimation_Gaussian_random_walk <- fit_AR1_Gaussian(y_missing_numeric[, 1], random_walk = TRUE)
  expect_equal(estimation_Gaussian_random_walk, estimation_Gaussian_random_walk_check)
})


