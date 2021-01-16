context("Function \"fit_VAR_t()\"")

data(ts_VAR_t)
Y <- ts_VAR_t$Y



test_that("error control works",{
  expect_error(fit_VAR_t(install.packages), "\"Y\" must be coercible to a matrix.")
  expect_error(fit_VAR_t(Y, -1), "\"p\" must be a positive integer.")
  expect_error(fit_VAR_t(Y, 2, omit_missing = 3), "\"omit_missing\" must be a logical value.")
  expect_error(fit_VAR_t(Y, 2, parallel_max_cores = -1), "\"parallel_max_cores\" must be a positive integer.")
  expect_error(fit_VAR_t(Y, 2, verbose = 3), "\"verbose\" must be a logical value.")
  expect_error(fit_VAR_t(Y, 2, return_iterates = 3), "\"return_iterates\" must be a logical value.")
  expect_error(fit_VAR_t(Y, 2, L = -1), "\"L\" must be a positive integer.")
  expect_error(fit_VAR_t(Y, 2, maxiter = -1), "\"maxiter\" must be greater than 1.")
  expect_error(fit_VAR_t(Y, 2, ptol = -1), "\"ptol\" must be greater than 0.")
  expect_error(fit_VAR_t(Y, 2, partition_groups = 3), "\"partition_groups\" must be a logical value.")
  expect_error(fit_VAR_t(Y, 2, K = -1), "\"K\" must be a positive integer.")
}) 


test_that("SAEM-MCMC works", {
  # fitted_VAR_check <- fit_VAR_t(Y, 2)[1:4]
  # save(fitted_VAR_check, file = "fitted_VAR_check.RData", version = 2, compress = "xz")
  load("fitted_VAR_check.RData")
  
  fitted_VAR <- fit_VAR_t(Y, 2)[1:4]
  expect_equal(fitted_VAR, fitted_VAR_check, tolerance = 1e-1)
})


test_that("omit-variable method works", {
  # fitted_VAR_omit_check <- fit_VAR_t(Y, 2, omit_missing = TRUE)[1:4]
  # save(fitted_VAR_omit_check, file = "fitted_VAR_omit_check.RData", version = 2, compress = "xz")
  load("fitted_VAR_omit_check.RData")
  
  fitted_VAR_omit <- fit_VAR_t(Y, 2, omit_missing = TRUE)[1:4]
  expect_equal(fitted_VAR_omit, fitted_VAR_omit_check, tolerance = 1e-3)
})
