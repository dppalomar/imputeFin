## Changes in imputeFin version 0.1.2 (2021-02-19)

* New function `fit_VAR_t()` added.

* New wrapper function `impute_rolling_AR1_Gaussian()` to impute on a rolling window basis for long time series.

* New wrapper function `impute_OHLC()` to impute the OHLC conveniently.


## Changes in imputeFin version 0.1.1 (2020-07-11)

* Included detection of outliers in all functions with the argument: remove_outliers = FALSE

* Removed the option to impute leading and trailing NAs as it doesn't make much sense and makes
  the code simpler.
  
* Vignette and README revised to illustrate the detection and removal of outliers.

* Additional argument to output summary messages: verbose = TRUE.


## Changes in imputeFin version 0.1.0 (2019-12-05)

* Initial release is on CRAN.
