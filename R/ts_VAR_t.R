#' Synthetic Student's t VAR data with missing values
#' 
#' Synthetic Student's t VAR data with missing values for 
#' estimation and imputation testing purposes.
#' 
#' @docType data
#'
#' @usage data(ts_VAR_t)
#'
#' @format List with the following elements:
#' \describe{
#'   \item{Y}{200 x 3 \code{zoo} object as a Student's t VAR time series.}
#'   \item{phi0}{True value of the constant vector in the VAR model.}
#'   \item{Phii}{True value of the coefficient matrix in the VAR model.}
#'   \item{scatter}{True value of the scatter matrix (of the noise distribution) in the VAR model.}
#'   \item{nu}{True value of the degrees of freedom (of the noise distribution) in the VAR model.}
#' }
#' 
#' @keywords dataset
#' 
"ts_VAR_t"