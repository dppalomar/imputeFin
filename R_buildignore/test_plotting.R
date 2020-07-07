library(imputeFin)
data(ts_AR1_Gaussian)

y_missing <- ts_AR1_Gaussian$y_missing[, 3]
plot_imputed(y_missing, title = "Original time series with missing values")
y_imputed <- impute_AR1_Gaussian(y_missing)
plot_imputed(y_imputed)

y_missing[1:10]
y_imputed[1:10]
