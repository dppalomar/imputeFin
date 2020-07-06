#' imputeFin: Imputation of Financial Time Series with Missing Values.
#'
#' Missing values often occur in financial data due to a variety 
#' of reasons (errors in the collection process or in the processing stage, 
#' lack of asset liquidity, lack of reporting of funds, etc.). However, 
#' most data analysis methods expect complete data and cannot be employed 
#' with missing values. One convenient way to deal with this issue without 
#' having to redesign the data analysis method is to impute the missing 
#' values. This package provides an efficient way to impute the missing 
#' values based on modeling the time series with a random walk or an 
#' autoregressive (AR) model, convenient to model log-prices and log-volumes 
#' in financial data. In the current version, the imputation is 
#' univariate-based (so no asset correlation is used). In addition,
#' outliers can be detected and removed.
#' 
#' @section Functions:
#' \code{\link{fit_AR1_Gaussian}}, \code{\link{impute_AR1_Gaussian}},
#' \code{\link{fit_AR1_t}}, \code{\link{impute_AR1_t}},
#' \code{\link{plot_imputed}}
#'
#' @section Data:
#' \code{\link{ts_AR1_Gaussian}}, \code{\link{ts_AR1_t}}
#'
#' @section Help:
#' For a quick help see the README file:
#' \href{https://github.com/dppalomar/imputeFin/blob/master/README.md}{GitHub-README}.
#'
#' For more details see the vignette:
#' \href{https://CRAN.R-project.org/package=imputeFin/vignettes/ImputeFinancialTimeSeries.html}{CRAN-vignette}.
#'
#' @author Junyan LIU and Daniel P. Palomar
#' 
#' @references
#' J. Liu, S. Kumar, and D. P. Palomar, "Parameter estimation of heavy-tailed AR model with missing 
#' data via stochastic EM," IEEE Trans. on Signal Processing, vol. 67, no. 8, pp. 2159-2172, 15 April, 2019.
#' <https://doi.org/10.1109/TSP.2019.2899816>
#'
#' @docType package
#' @name imputeFin-package
NULL
