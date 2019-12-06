#' Synthetic AR(1) Student's t time series with missing values.
#' 
#' Synthetic AR(1) Student's t time series with missing values for 
#' estimation and imputation testing purposes.
#' 
#' @docType data
#'
#' @usage data(ts_AR1_t)
#'
#' @format List with the following elements:
#' \describe{
#'   \item{y_missing}{300 x 3 \code{zoo} object with three AR(1) Student's t time 
#'                    series along the columns: the first column contains a
#'                    time series with 10\% consecutive missing values; the
#'                    second column contains a time series with 10\% missing 
#'                    values randomly distributed; and the third column contains 
#'                    the union of the previous missing values.}
#'   \item{phi0}{Value of \code{phi0} used to generate the time series.}
#'   \item{phi1}{Value of \code{phi1} used to generate the time series.}
#'   \item{sigma2}{Value of \code{sigma2} used to generate the time series.}
#'   \item{nu}{Value of \code{nu} used to generate the time series.}
#' }
#' 
#' @keywords dataset
#' 
"ts_AR1_t"