library(imputeFin)
data(ts_AR1_Gaussian)
y_missing <- ts_AR1_Gaussian$y_missing
y_missing[100, 1] <- 2*y_missing[100, 1]

# indices of NAs
idx_NAs     <- apply(y_missing, MARGIN = 2, FUN = function(x) which(is.na(as.numeric(x))))
idx_outlier <- 100
idx_NAs <- union(idx_NAs, idx_outlier)


y_imputed <- impute_AR1_t(y_missing, n_samples = 1, return_estimates = FALSE, remove_outliers = TRUE)
attr(y_imputed, "index_miss")
attr(y_imputed, "index_outliers")

res <- impute_AR1_t(y_missing, n_samples = 1, return_estimates = TRUE, remove_outliers = TRUE)
attr(res$y_imputed, "index_miss")
attr(res$y_imputed, "index_outliers")

res <- impute_AR1_t(y_missing, n_samples = 3, return_estimates = TRUE, remove_outliers = TRUE)
attr(res$y_imputed.1, "index_miss")
attr(res$y_imputed.1, "index_outliers")


y_imputed <- impute_AR1_Gaussian(y_missing, n_samples = 1, return_estimates = FALSE, remove_outliers = TRUE)
attr(y_imputed, "index_miss")
attr(y_imputed, "index_outliers")

res <- impute_AR1_Gaussian(y_missing, n_samples = 1, return_estimates = TRUE, remove_outliers = TRUE)
attr(res$y_imputed, "index_miss")
attr(res$y_imputed, "index_outliers")

res <- impute_AR1_Gaussian(y_missing, n_samples = 3, return_estimates = TRUE, remove_outliers = TRUE)
attr(res$y_imputed.1, "index_miss")
attr(res$y_imputed.1, "index_outliers")
