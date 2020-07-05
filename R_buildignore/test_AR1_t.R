library(imputeFin)
# source("estimate_impute_AR1_Gaussian.R")
# source("estimate_impute_AR1_t.R")
# source("outlier_and_plot.R")
data(ts_AR1_t)

# test the estimation function
#y_missing <- ts_AR1_t$y_missing
y_missing <- ts_AR1_t$y_missing[, 3, drop = FALSE]
#y_missing <- c(1,2,3.4,4.1,NA,NA,7,8,15,10) # outlier test
estimation_result1 <- fit_AR1_t(y_missing, random_walk = FALSE, zero_mean = FALSE, fast_and_heuristic = TRUE, remove_outlier = TRUE,
                                   return_iterates = FALSE, return_condMean_Gaussian = FALSE,
                                  tol = 1e-10,  maxiter = 100, n_chain = 10, n_thin = 1, K = 30)
estimation_result2 <- fit_AR1_t(y_missing, random_walk = FALSE, zero_mean = FALSE, fast_and_heuristic = FALSE,
                                   return_iterates = FALSE, return_condMean_Gaussian = FALSE,
                                   tol = 1e-10,  maxiter = 100, n_chain = 10, n_thin = 1, K = 30)

#test outlier detection
y_missing[100] <- 100
estimation_result3 <- fit_AR1_t(y_missing, random_walk = FALSE, zero_mean = FALSE, fast_and_heuristic = TRUE, remove_outlier = TRUE,
                                return_iterates = FALSE, return_condMean_Gaussian = FALSE,
                                tol = 1e-10,  maxiter = 100, n_chain = 10, n_thin = 1, K = 30)
estimation_result4 <- fit_AR1_t(y_missing, random_walk = FALSE, zero_mean = FALSE, fast_and_heuristic = TRUE, remove_outlier = FALSE,
                                return_iterates = FALSE, return_condMean_Gaussian = FALSE,
                                tol = 1e-10,  maxiter = 100, n_chain = 10, n_thin = 1, K = 30)

# test the imputation function and plot function
imputation_result <- impute_AR1_t(y_missing, n_samples = 3, random_walk = FALSE, zero_mean = FALSE, 
                                fast_and_heuristic = FALSE, return_estimates = TRUE, tol = 1e-11,  maxiter = 90,K = 20,
                                n_burn = 100, n_thin = 50)
y_imputed = imputation_result$y_imputed.1
plotImputed(y_imputed, column = 1, type = "simple")



