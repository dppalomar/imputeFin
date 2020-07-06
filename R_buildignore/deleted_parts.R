#
# Deleted parts from function impute_AR1_t:
#
#' @param impute_leading_NAs Logical value indicating if the leading missing values of time 
#'                           series are to be imputed (default is \code{FALSE}).
#' @param impute_trailing_NAs Logical value indicating if the trailing missing values of time 
#'                            series are to be imputed (default is \code{FALSE}).                      

impute_AR1_t <- function(y, n_samples = 1, impute_leading_NAs = FALSE, impute_trailing_NAs = FALSE,
                         random_walk = FALSE, zero_mean = FALSE, 
                         fast_and_heuristic = TRUE, remove_outliers = FALSE,
                         return_estimates = FALSE,
                         tol = 1e-10,  maxiter = 100, K = 30,
                         n_burn = 100, n_thin = 50) {
  
 
  
  # if there are missing values at the head of the time series and impute_leading_NAs == TRUE, impute them.
  if (index_obs_min > 1 & impute_leading_NAs) { 
    for (j in (index_obs_min - 1):1 )
      y_imputed[j, ] <- ( y_imputed[j + 1, ] - rt(n_samples, nu) * sqrt(sigma2) - phi0 )/phi1
  }
  
  # if there are missing values at the tail of the time series and impute_trailing_NAs == TRUE, impute them.
  if (index_obs_max < length(y) & impute_trailing_NAs){
    for (i in (index_obs_max + 1):length(y))
      y_imputed[i, ] <- phi0 + phi1 * y_imputed[i - 1, ] +  rt(n_samples, nu) * sqrt(sigma2)
  }
  
   
}