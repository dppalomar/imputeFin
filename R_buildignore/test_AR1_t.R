library(imputeFin)
data(AR1_t)

# test the estimation function
# y_missing <- AR1_Gaussian$y_missing_numeric
y_missing <- AR1_t$y_missing_xts
estimation_result <- estimateAR1t(y_missing, random_walk = FALSE, zero_mean = FALSE,
                                  iterates = TRUE, condMean_Gaussian = TRUE,
                                  n_chain = 10, n_thin = 1, n_iter = 100, K = 30)



# test the imputation function and plot function
imputation_result <- imputeAR1t(y_missing, n_samples = 3, random_walk = FALSE, zero_mean = FALSE,
                                n_burn = 100, n_thin = 50,
                                estimates = TRUE)
y_imputed = imputation_result$y_imputed.1
plotImputed(y_imputed, column = 1, type = "simple")




