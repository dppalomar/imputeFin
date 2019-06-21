#' imputeFin: Imputation of Financial Time Series with Missing Values.
#'
#' Missing values often occur in financial data due to a variety 
#' of reasons (errors in the collection process or in the processing stage, 
#' lack of asset liquidity, lack of reporting of funds, etc.). However, 
#' most data analysis methods expect complete data and cannot be employed 
#' with missing values. One convenient way to deal with this issue without 
#' having to redesign the data analysis method is to impute the missing 
#' values. This package provides an efficient way to impute the missing 
#' values based on modeling the time series with an autoregressive (AR) model,
#' convenient to model log-prices and log-volume in financial data. In the
#' current version, the imputation is univariate-based (so no asset correlation 
#' is used).
#' 
#' @section Functions:
#' \code{\link{estimateAR1Gaussian}}, \code{\link{imputeAR1Gaussian}},
#' \code{\link{estimateAR1t}}, \code{\link{imputeAR1t}},
#' \code{\link{plotImputed}}
#'
#' @section Data:
#' \code{\link{AR1_Gaussian}}, \code{\link{AR1_t}}
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
#' J. Liu, S. Kumar, and D. P. Palomar (2019). Parameter Estimation of
#' Heavy-Tailed AR Model With Missing Data Via Stochastic EM. \emph{IEEE Trans. on Signal Processing},
#' vol. 67, no. 8, pp. 2159-2172. <https://doi.org/10.1109/TSP.2019.2899816>
#'
#' @docType package
#' @name imputeFin-package
NULL
