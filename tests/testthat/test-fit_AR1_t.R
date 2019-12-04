context("Function \"fit_AR1_t()\"")
#library(testthat)

# # generate the data
# phi0 <- 0
# phi1 <- 1
# sigma2 <- 0.01
# nu <- 2
# n <- 300
# n_miss <- 0.1 * n
# n_drop <- 100
# n_total <- n + n_drop
# data <- vector(length = n_total)
# epsilon <- vector(length = n_total - 1)# innovations
# data[1] <- 0
# for (i in 2:n_total) {
#   epsilon[i] <- rt(1, nu) * sqrt(sigma2)
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
# index_miss2 <- c(1,sort(sample(2:(n - 1), n_miss - 2)), n)
# y_missing_numeric[index_miss1, 1] <- NA
# y_missing_numeric[index_miss2, 2] <- NA
# y_missing <- zoo::zoo(y_missing_numeric, seq(as.Date("2016-01-01"), length = n, by = "days"))
# save(y_missing_numeric, y_missing, phi0, phi1, sigma2, nu, file = "t_data.RData", version = 2, compress = "xz")

load("t_data.RData")


test_that("error control works",{
  expect_error(fit_AR1_t(y = median), "\"y\" must be coercible to a vector or matrix.")
  expect_error(fit_AR1_t(y = "Hongkong"), "\"y\" only allows numerical or NA values.")
  expect_error(fit_AR1_t(y = c(1, 2, NA, 3, 4, NA)), "Each time series in \"y\" must have at least five observations.")
  expect_error(fit_AR1_t(y = y_missing_numeric, tol = -0.1), "\"tol\" must be greater than 0.")
  expect_error(fit_AR1_t(y = y_missing_numeric, maxiter = -2), "\"maxiter\" must be greater than 1.")
  expect_error(fit_AR1_t(y = y_missing_numeric, n_chain = -1), "\"n_chain\" must be a positive integer.")
  expect_error(fit_AR1_t(y = y_missing_numeric, n_thin = -1), "\"n_thin\" must be a positive integer.")
  expect_error(fit_AR1_t(y = y_missing_numeric, K = 0), "\"K\" must be a positive integer.")
}) 


test_that("default mode works", {
  # estimation_t_uni_check <- fit_AR1_t(y_missing_numeric[,1])
  # estimation_t_multi_check <- fit_AR1_t(y_missing_numeric)
  # save(estimation_t_uni_check, estimation_t_multi_check, file = "estimation_t_check.RData", version = 2, compress = "xz")
  load("estimation_t_check.RData")
  
  # check the uni time series & numeric 
  estimation_t_uni <- fit_AR1_t(y_missing_numeric[, 1])
  expect_equal(estimation_t_uni, estimation_t_uni_check)
  
  # check the multiple time series & numeric  
  estimation_t_multi <- fit_AR1_t(y_missing_numeric)
  expect_equal(estimation_t_multi, estimation_t_multi_check)
  
  # check the uni time series & zoo 
  estimation_t_uni <- fit_AR1_t(y_missing[, 1])
  expect_equal(estimation_t_uni, estimation_t_uni_check)
  
  # echeck the multiple time series & zoo
  estimation_t_zoo <- fit_AR1_t(y_missing)
  expect_equal(estimation_t_zoo, estimation_t_multi_check)
})


test_that("time series with zero mean works", {
  # estimation_t_zero_mean_check <- fit_AR1_t(y_missing_numeric[,1], zero_mean = TRUE)
  # save(estimation_t_zero_mean_check, file = "estimation_t_zero_mean_check.RData", version = 2, compress = "xz")
  load("estimation_t_zero_mean_check.RData")
  
  estimation_t_zero_mean <- fit_AR1_t(y_missing_numeric[, 1], zero_mean = TRUE)
  expect_equal(estimation_t_zero_mean, estimation_t_zero_mean_check)
})


test_that("time series following random walk works", {
  # estimation_t_random_walk_check <- fit_AR1_t(y_missing_numeric[,1], random_walk = TRUE)
  # save(estimation_t_random_walk_check, file = "estimation_t_random_walk_check.RData", version = 2, compress = "xz")
  load("estimation_t_random_walk_check.RData")
  
  estimation_t_random_walk <- fit_AR1_t(y_missing_numeric[, 1], random_walk = TRUE)
  expect_equal(estimation_t_random_walk, estimation_t_random_walk_check)

})

test_that("SAEM-MCMC algorithm works", {
  estimation_t_SAEM <- fit_AR1_t(y_missing_numeric[,1], fast_and_heuristic = FALSE)
  expect_equal(estimation_t_SAEM$phi0, phi0, tolerance = 0.2)
  expect_equal(estimation_t_SAEM$phi1, phi1, tolerance = 0.1)
  expect_equal(estimation_t_SAEM$sigma2, sigma2, tolerance = 0.003)
  expect_equal(estimation_t_SAEM$sigma2, sigma2, tolerance = 0.4)
})
