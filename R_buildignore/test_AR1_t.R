library(imputeFin)
data(ts_AR1_t)

# test the estimation function
#y_missing <- ts_AR1_t$y_missing
y_missing <- ts_AR1_t$y_missing[, 3, drop = FALSE]
#y_missing <- c(1,2,3.4,4.1,NA,NA,7,8,15,10)  # outlier test
estimation_result1 <- fit_AR1_t(y_missing, random_walk = FALSE, zero_mean = FALSE, fast_and_heuristic = TRUE, remove_outliers = TRUE,
                                return_iterates = FALSE, return_condMean_Gaussian = FALSE,
                                tol = 1e-10,  maxiter = 100, n_chain = 10, n_thin = 1, K = 30)
estimation_result2 <- fit_AR1_t(y_missing, random_walk = FALSE, zero_mean = FALSE, fast_and_heuristic = TRUE,
                                return_iterates = FALSE, return_condMean_Gaussian = FALSE,
                                tol = 1e-10,  maxiter = 100, n_chain = 10, n_thin = 1, K = 30)
all.equal(estimation_result1, estimation_result2)


#test outlier detection
y_missing[100] <- 100
estimation_result3 <- fit_AR1_t(y_missing, random_walk = FALSE, zero_mean = FALSE, fast_and_heuristic = TRUE, remove_outliers = TRUE,
                                return_iterates = FALSE, return_condMean_Gaussian = FALSE,
                                tol = 1e-10,  maxiter = 100, n_chain = 10, n_thin = 1, K = 30)
estimation_result4 <- fit_AR1_t(y_missing, random_walk = FALSE, zero_mean = FALSE, fast_and_heuristic = TRUE, remove_outliers = FALSE,
                                return_iterates = FALSE, return_condMean_Gaussian = FALSE,
                                tol = 1e-10,  maxiter = 100, n_chain = 10, n_thin = 1, K = 30)

# test the imputation function and plot function
imputation_result1 <- impute_AR1_t(y_missing, n_samples = 3, random_walk = FALSE, zero_mean = FALSE, 
                                   fast_and_heuristic = TRUE, remove_outliers = TRUE,
                                   return_estimates = TRUE, tol = 1e-11,  maxiter = 90, K = 20,
                                   n_burn = 100, n_thin = 50)
imputation_result2 <- impute_AR1_t(y_missing, n_samples = 3, random_walk = FALSE, zero_mean = FALSE, 
                                   fast_and_heuristic = TRUE, 
                                   return_estimates = TRUE, tol = 1e-11,  maxiter = 90, K = 20,
                                   n_burn = 100, n_thin = 50)


plot_imputed(imputation_result1$y_imputed.1)
plot_imputed(imputation_result2$y_imputed.1)
