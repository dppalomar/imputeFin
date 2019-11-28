library(imputeFin)
#source("estimate_impute_AR1_Gaussian.R")
data(ts_AR1_Gaussian)

# test the estimation function
# y_missing <- ts_AR1_Gaussian$y_missing[, 4, drop = FALSE]
y_missing <- ts_AR1_Gaussian$y_missing
estimation_result <- estimateAR1Gaussian(y_missing, random_walk = FALSE, zero_mean = FALSE,
                                         return_iterates = TRUE, return_condMeanCov = TRUE,
                                         tol = 1e-10,  maxiter = 1000)

# y_missing <- ts_AR1_Gaussian$y_missing[,1, drop = FALSE]
# estimation_result <- estimateAR1Gaussian(y_missing)

# test the imputation function
imputation_result <- imputeAR1Gaussian(y_missing, n_samples = 3, random_walk = FALSE, zero_mean = FALSE,
                                       return_estimates = TRUE)
y_imputed <- imputation_result$y_imputed.1
plotImputed(y_imputed, column = 1, type = "simple")




#########################################################
# comparison with other methods
# library("imputeTS")
# y_imputed_km = na.kalman(y)
# #plot 3
# { plot(y, main = "Imputated by kalman method")
#   lines(y_imputed_km[index_miss_p], col = "blue", lwd = 2) }
# 
# library("knitr")
# rst_orig <- ar(y_orig, order.max = 1, demean = TRUE)
# est_orig <- list("phi0" = rst_orig$x.mean * (1 - sum(rst_orig$ar)),
#                  "phi1" = rst_orig$ar,
#                  "sigma2" = rst_orig$var.pred)
# rst_km <- ar(y_imputed_km, order.max = 1, demean = TRUE)
# est_km <- list("phi0" = rst_km$x.mean * (1 - sum(rst_km$ar)),
#                  "phi1" = rst_km$ar,
#                  "sigma2" = rst_km$var.pred)
# rst_pr <- ar(y_imputed[, 1], order.max = 1, demean = TRUE)
# est_pr <- list("phi0" = rst_pr$x.mean * (1 - sum(rst_pr$ar)),
#                "phi1" = rst_pr$ar,
#                "sigma2" = rst_pr$var.pred)
# rst = cbind(est_orig, est_km, est_pr)
# print( kable(rst))
