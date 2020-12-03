#' @title Impute missing values of an OHLC time series on a rolling window basis based on a Gaussian AR(1) model
#'
#' @description Impute inner missing values (excluding leading and trailing ones) 
#'              of an OHLC time series on a rolling window basis. This is a wrapper 
#'              of the functions \code{\link{impute_AR1_Gaussian}} and 
#'              \code{\link{impute_rolling_AR1_Gaussian}}.
#' 
#' @param y Time series object coercible to a numeric matrix (e.g., \code{zoo} or \code{xts}) 
#'          with four columns denoting the prices Op, Hi, Lo, Cl.
#' @inheritParams impute_rolling_AR1_Gaussian
#' 
#' @return Imputed OHLC prices.
#' 
#' @author Daniel P. Palomar
#' 
#' @seealso \code{\link{impute_AR1_Gaussian}}, \code{\link{impute_rolling_AR1_Gaussian}}
#' 
#' @export
impute_OHLC <- function(y_OHLC, rolling_window = 252, 
                        remove_outliers = FALSE, outlier_prob_th = 1e-3,
                        tol = 1e-10, maxiter = 100) {
  logy_OHLC_imputed <- logy_OHLC <- log(y_OHLC)
  
  # impute Cl
  logy_OHLC_imputed[, 4] <- impute_rolling_AR1_Gaussian(logy_OHLC[, 4], rolling_window, 
                                                        random_walk = TRUE, zero_mean = FALSE, 
                                                        remove_outliers, outlier_prob_th, 
                                                        tol, maxiter)
  
  # impute Op-Cl and then form Op
  Op_Cl <- impute_rolling_AR1_Gaussian(logy_OHLC[, 1] - logy_OHLC_imputed[, 4], rolling_window, 
                                       random_walk = FALSE, zero_mean = TRUE, 
                                       remove_outliers, outlier_prob_th, tol, maxiter)
  logy_OHLC_imputed[, 1] <- logy_OHLC_imputed[, 4] + Op_Cl
  
  # impute Hi-Cl and then form Hi
  Hi_Cl <- impute_rolling_AR1_Gaussian(logy_OHLC[, 2] - logy_OHLC_imputed[, 4], rolling_window, 
                                       random_walk = FALSE, zero_mean = TRUE, 
                                       remove_outliers, outlier_prob_th, tol, maxiter)
  Hi_Cl <- pmax(Hi_Cl, 0)
  logy_OHLC_imputed[, 2] <- logy_OHLC_imputed[, 4] + Hi_Cl
  
  # impute Lo-Cl and then form Lo
  Lo_Cl <- impute_rolling_AR1_Gaussian(logy_OHLC[, 3] - logy_OHLC_imputed[, 4], rolling_window, 
                                       random_walk = FALSE, zero_mean = TRUE, 
                                       remove_outliers, outlier_prob_th, tol, maxiter)
  Lo_Cl <- pmin(Lo_Cl, 0)
  logy_OHLC_imputed[, 3] <- logy_OHLC_imputed[, 4] + Lo_Cl  
  
  return(exp(logy_OHLC_imputed))
}