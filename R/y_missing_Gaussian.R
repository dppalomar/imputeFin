#' Synthetic Gaussian AR(1) time series containing missing values for estimation and imputation testing purposes
#'
#' @docType data
#'
#' @usage data(AR1_Gaussian)
#'
#' @format A list consisting of the values of parameters of a Gaussian AR(1) model and a zoo object 
#' containing three time series generted from this model. The first and second time series has 10% 
#' missing values, while the third time series is complete. The locations of missing values in the first 
#' time series are consecutive, while the locations of missing values in the second are randomly distributed.
#'
#' \describe{
#'   \item{phi0}{value of phi0 used to generate the time series}
#'   \item{phi1}{value of phi1 used to generate the time series}
#'   \item{sigma2}{value of sigma2 used to generate the time series}
#'   \item{y_missing}{zoo object consisting of three time series with missing values}
#'}
"AR1_Gaussian"
