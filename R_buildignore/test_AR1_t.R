library(imputeFin)
data(AR1_t)

# test the estimation function
y_missing <- AR1_t$y_missing
estimation_result1 <- estimateAR1t(y_missing, random_walk = FALSE, zero_mean = FALSE, method = "heuristic",
                                  iterates = fALSE, condMean_Gaussian = FALSE,
                                  tol = 1e-10,  maxiter = 100, n_chain = 10, n_thin = 1, K = 30)
estimation_result2 <- estimateAR1t(y_missing, random_walk = FALSE, zero_mean = FALSE, method = "stEM",
                                   iterates = FALSE, condMean_Gaussian = FALSE,
                                   tol = 1e-10,  maxiter = 100, n_chain = 10, n_thin = 1, K = 30)



# test the imputation function and plot function
imputation_result <- imputeAR1t(y_missing, n_samples = 3, random_walk = FALSE, zero_mean = FALSE, 
                                method = "heuristic", estimates = FALSE,
                                n_burn = 100, n_thin = 50)
y_imputed = imputation_result$y_imputed.1
plotImputed(y_imputed, column = 1, type = "simple")



