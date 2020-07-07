library(imputeFin)
data(ts_AR1_t)

# test the estimation function
#y_missing <- ts_AR1_t$y_missing
y_missing <- ts_AR1_t$y_missing[, 3, drop = FALSE]
#y_missing <- c(1,2,3.4,4.1,NA,NA,7,8,15,10)  # outlier test
estimation_result1 <- fit_AR1_t(y_missing, remove_outliers = TRUE)
estimation_result2 <- fit_AR1_t(y_missing)
all.equal(estimation_result1, estimation_result2)


#test outlier detection
y_missing <- ts_AR1_t$y_missing[, 3, drop = FALSE]
y_missing[100] <- 100
estimation_result3 <- fit_AR1_t(y_missing, remove_outliers = TRUE)
estimation_result4 <- fit_AR1_t(y_missing, remove_outliers = FALSE)

# test the imputation function and plot function
imputation_result0 <- impute_AR1_t(y_missing, remove_outliers = TRUE)
imputation_result1 <- impute_AR1_t(y_missing, n_samples = 3, remove_outliers = TRUE)
imputation_result2 <- impute_AR1_t(y_missing, n_samples = 3)


plot_imputed(y_missing, title = "Original signal")
plot_imputed(imputation_result0)
plot_imputed(imputation_result1$y_imputed.1)
plot_imputed(imputation_result2$y_imputed.1)

#res <- fit_AR1_t(ts_AR1_t$y_missing[, 3, drop = FALSE], remove_outliers = TRUE)
#res <- fit_AR1_t(ts_AR1_Gaussian$y_missing, remove_outliers = TRUE)
#res <- impute_AR1_t(ts_AR1_Gaussian$y_missing[, 3, drop = FALSE], remove_outliers = TRUE)
#res <- impute_AR1_t(ts_AR1_Gaussian$y_missing, remove_outliers = TRUE)


# price jump instead of outlier (due to split of dividend)
y_missing <- ts_AR1_t$y_missing[, 3, drop = FALSE]
y_missing[100:length(y_missing)] <- 10 + y_missing[100:length(y_missing)]
plot_imputed(y_missing, title = "Original signal")
y_imputed <- impute_AR1_t(y_missing, remove_outliers = TRUE)
plot_imputed(y_imputed)

