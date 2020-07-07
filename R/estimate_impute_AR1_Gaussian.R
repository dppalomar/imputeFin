#' @title Fit Gaussian AR(1) model to time series with missing values and/or outliers
#'
#' @description Estimate the parameters of a univariate Gaussian AR(1) model 
#'              to fit the given time series with missing values and/or outliers. 
#'              For multivariate time series, the function will perform a 
#'              number of individual univariate fittings without attempting 
#'              to model the correlations among the time series.
#'              If the time series does not contain missing values, the 
#'              maximum likelihood (ML) estimation is done in one shot.
#'              With missing values, the iterative EM algorithm is employed 
#'              for the estimation until converge is achieved.
#'
#' @param y Time series object coercible to either a numeric vector or numeric matrix 
#'          (e.g., \code{zoo} or \code{xts}) with missing values denoted by \code{NA}. 
#' @param random_walk Logical value indicating if the time series is assumed to be a random walk so that \code{phi1 = 1} 
#'                    (default is \code{FALSE}).
#' @param zero_mean Logical value indicating if the time series is assumed zero-mean so that \code{phi0 = 0} 
#'                  (default is \code{FALSE}).
#' @param remove_outliers Logical value indicating whether to detect and remove outliers.
#' @param return_iterates Logical value indicating if the iterates are to be returned (default is \code{FALSE}).
#' @param return_condMeanCov Logical value indicating if the conditional mean and covariance matrix of the 
#'                           time series (excluding the leading and trailing missing values) given the observed data are to be returned (default is \code{FALSE}).
#' @param tol Positive number denoting the relative tolerance used as stopping criterion (default is \code{1e-8}).
#' @param maxiter Positive integer indicating the maximum number of iterations allowed (default is \code{100}).
#' 
#' @return If the argument \code{y} is a univariate time series (i.e., coercible to a numeric vector), then this 
#'         function will return a list with the following elements:
#' \item{\code{phi0}}{The estimate for \code{phi0} (real number).}
#' \item{\code{phi1}}{The estimate for \code{phi1} (real number).}
#' \item{\code{sigma2}}{The estimate for \code{sigma^2} (positive number).}
#' \item{\code{phi0_iterates}}{Numeric vector with the estimates for \code{phi0} at each iteration
#'                            (returned only when \code{return_iterates = TRUE}).}
#' \item{\code{phi1_iterates}}{Numeric vector with the estimates for \code{phi1} at each iteration
#'                            (returned only when \code{return_iterates = TRUE}).}
#' \item{\code{sigma2_iterates}}{Numeric vector with the estimates for \code{sigma^2} at each iteration
#'                              (returned only when \code{return_iterates = TRUE}).}
#' \item{\code{f_iterates}}{Numeric vector with the objective values at each iteration
#'                          (returned only when \code{return_iterates = TRUE}).}
#' \item{\code{cond_mean_y}}{Numeric vector (of same length as argument \code{y}) with the conditional mean of the time series 
#'                           (excluding the leading and trailing missing values)
#'                           given the observed data (returned only when \code{return_condMeanCov = TRUE}).}
#' \item{\code{cond_cov_y}}{Numeric matrix (with number of columns/rows equal to the length of the argument \code{y})
#'                          with the conditional covariance matrix of the time series (excluding the leading and trailing missing values) 
#'                          given the observed data (returned only when \code{return_condMeanCov = TRUE}).}
#' \item{\code{index_outliers}}{Indices of outliers detected (returned only when \code{remove_outliers = TRUE} and outliers are detected).}
#' 
#' If the argument \code{y} is a multivariate time series (i.e., with multiple columns and coercible to a numeric matrix), 
#' then this function will return a list with each element as in the case of univariate \code{y} corresponding to each
#' of the columns (i.e., one list element per column of \code{y}), with the following additional elements that combine the 
#' estimated values in a convenient vector form:
#' \item{\code{phi0_vct}}{Numeric vector (with length equal to the number of columns of \code{y})
#'                        with the estimates for \code{phi0} for each of the univariate time series.}
#' \item{\code{phi1_vct}}{Numeric vector (with length equal to the number of columns of \code{y})
#'                        with the estimates for \code{phi1} for each of the univariate time series.}
#' \item{\code{sigma2_vct}}{Numeric vector (with length equal to the number of columns of \code{y})
#'                        with the estimates for \code{sigma2} for each of the univariate time series.}
#' 
#' @author Junyan Liu and Daniel P. Palomar
#' 
#' @seealso \code{\link{impute_AR1_Gaussian}}, \code{\link{fit_AR1_t}}
#' 
#' @references
#' R. J. Little and D. B. Rubin, Statistical Analysis with Missing Data, 2nd ed. Hoboken, N.J.: John Wiley & Sons, 2002.
#' 
#' J. Liu, S. Kumar, and D. P. Palomar, "Parameter estimation of heavy-tailed AR model with missing 
#' data via stochastic EM," IEEE Trans. on Signal Processing, vol. 67, no. 8, pp. 2159-2172, 15 April, 2019. 
#' 
#' @examples 
#' library(imputeFin)
#' data(ts_AR1_Gaussian)
#' y_missing <- ts_AR1_Gaussian$y_missing
#' fitted <- fit_AR1_Gaussian(y_missing)
#' 
#' @import zoo
#' @export
fit_AR1_Gaussian <- function(y, random_walk = FALSE, zero_mean = FALSE, remove_outliers = FALSE,
                             return_iterates = FALSE, return_condMeanCov = FALSE,
                             tol = 1e-8, maxiter = 100) {
  # error control
  if (!is.matrix(try(as.matrix(y), silent = TRUE))) stop("\"y\" must be coercible to a vector or matrix.")
  if (tol <= 0) stop("\"tol\" must be greater than 0.")
  if (maxiter < 1) stop("\"maxiter\" must be greater than 1.")
  
  # manage multiple columns
  if (NCOL(y) > 1) {
    estimation_list <- apply(y, MARGIN = 2, FUN = fit_AR1_Gaussian, random_walk, zero_mean, remove_outliers, 
                             return_iterates, return_condMeanCov, tol, maxiter)
    phi0_vct   <- unlist(lapply(estimation_list, function(x) x$phi0))
    phi1_vct   <- unlist(lapply(estimation_list, function(x) x$phi1))
    sigma2_vct <- unlist(lapply(estimation_list, function(x) x$sigma2))
    return(c(estimation_list, list("phi0_vct"   = phi0_vct,
                                   "phi1_vct"   = phi1_vct,
                                   "sigma2_vct" = sigma2_vct)))
  }
  
  #
  #   code for y single-column
  #
  # error control
  if (!is.numeric(y)) stop("\"y\" only allows numerical or NA values.")
  if (sum(!is.na(y)) < 5L) stop("Each time series in \"y\" must have at least 5 observations.")
  y <- as.numeric(y)
  
  # remove the missing values at the head and tail of the time series since they do not affect the estimation result
  index_obs <- which(!is.na(y))
  y <- y[min(index_obs):max(index_obs)]
  idx_offset <- min(index_obs) - 1L
  index_obs <- which(!is.na(y))

  # outlier detection
  if(remove_outliers) {
    fitted_with_outliers <- if (!anyNA(y)) fit_AR1_Gaussian_complete(y, random_walk, zero_mean)
                            else fit_AR1_Gaussian(y, random_walk, zero_mean, remove_outliers = FALSE, return_iterates, return_condMeanCov, tol, maxiter)
    idx_outliers <- find_outliers_AR1_Gaussian(y, index_obs, fitted_with_outliers, p_confidence = 0.95)
    if (!is.null(idx_outliers))
      y[idx_outliers] <- NA  # substitute outliers with NAs
  }
  
  # estimation (after possibly setting outliers to NA)
  if (!anyNA(y))   # trivial case with no NAs
    results <- fit_AR1_Gaussian_complete(y, random_walk, zero_mean)
  else {
    # if there are NAs find the missing blocks
    n <- index_obs <- index_miss <- n_obs <- y_obs <- delta_index_obs <- n_block <- n_in_block <- 
      first_index_in_block <- last_index_in_block <- previous_obs_before_block <- next_obs_after_block <- NULL
    list2env(findMissingBlock(y), envir = environment())

    # objective function, the observed data log-likelihood
    if (return_iterates)
      obj <- function(phi0, phi1, sigma2) {
        sum_phi1 <- sum2_phi1 <- rep(NA, n_obs-1)
        for (i in 1:(n_obs-1)) {
          sum_phi1[i] <- sum(phi1^(0:(delta_index_obs[i] - 1)))
          sum2_phi1[i] <- sum(phi1^((0:(delta_index_obs[i] - 1))*2))
        }
        return(-sum((y_obs[2:n_obs] - phi1^delta_index_obs * y_obs[1:(n_obs - 1)] 
                     - phi0 * sum_phi1)^2 / (2 * sum2_phi1 * sigma2))
               -0.5 * (n_obs - 1) * log(sigma2) - 0.5 * sum(log(sum2_phi1)))
      }
    
    # initialize the estimates
    phi0 <- phi1 <- sigma2 <- f <- c()
    estimation_heuristic <- fit_AR1_Gaussian_heuristic(y, index_miss, random_walk, zero_mean)
    phi1[1] <- estimation_heuristic$phi1
    phi0[1] <- estimation_heuristic$phi0
    sigma2[1] <- estimation_heuristic$sigma2
    if (return_iterates) f[1] <- obj(phi0[1], phi1[1], sigma2[1])    
    
    # loop
    for (k in 1:maxiter) {
      # E-step
      # computation of mean and covariance of y conditional on all the observed data
      cond <- condMeanCov(y_obs, index_obs, n, n_block, n_in_block, 
                          first_index_in_block, last_index_in_block, previous_obs_before_block, next_obs_after_block, 
                          phi0[k], phi1[k], sigma2[k], full_cov = FALSE)
      # computation of sufficient statistics
      s_y2   <- sum(cond$mean_y[-1])
      s_y1   <- sum(cond$mean_y[-n])
      s_y2y2 <- sum(cond$cov_y_diag[-1] + cond$mean_y[-1]^2)  #slow version: s_y2y2 <- sum(diag(cond_cov_y[2:n, 2:n]) + cond_mean_y[2:n]^2)
      s_y1y1 <- sum(cond$cov_y_diag[-n]  + cond$mean_y[-n]^2)  #slow version: s_y1y1 <- sum(diag(cond_cov_y[1:(n - 1), 1:(n - 1)])  + cond_mean_y[1:(n - 1)]^2)
      s_y2y1 <- sum(cond$cov_y_diag1 + cond$mean_y[-1] * cond$mean_y[-n])  #slow version: s_y2y1 <- sum(diag(cond_cov_y[1:(n - 1), 2:n]) + cond_mean_y[2:n] * cond_mean_y[1:(n - 1)])
      
      # M-step (update the estimates)
      if (!random_walk && !zero_mean) {
        phi1[k + 1] <- (s_y2y1 - s_y2 * s_y1 / (n - 1)) / (s_y1y1 - s_y1 * s_y1 / (n - 1)) 
        phi0[k + 1] <- (s_y2 - phi1[k +1] * s_y1) / (n - 1)
      } else if (random_walk && !zero_mean){
        phi1[k + 1] <- 1
        phi0[k + 1] <- (s_y2 - s_y1) / (n - 1)
      } else if (!random_walk && zero_mean){
        phi1[k + 1] <- s_y2y1 / s_y1y1 
        phi0[k + 1] <- 0
      } else{
        phi1[k + 1] <- 1
        phi0[k + 1] <- 0
      }
      sigma2[k + 1] <- ((s_y2y2 + (n - 1) * phi0[k + 1]^2 + phi1[k + 1]^2 * s_y1y1 -
                           2 * phi0[k + 1] * s_y2 - 2 * phi1[k + 1] * s_y2y1 + 2 * phi0[k + 1] * phi1[k + 1] * s_y1) / (n - 1))
      
      # computation of the objective function
      if (return_iterates) f[k + 1] <- obj(phi0[k + 1], phi1[k + 1], sigma2[k + 1])
      
      # termination criterion    
      if (abs(phi0[k + 1] - phi0[k]) <= tol * (abs(phi0[k + 1]) + abs(phi0[k]))/2
          && abs(phi1[k + 1] - phi1[k]) <= tol * (abs(phi1[k + 1]) + abs(phi1[k]))/2
          && abs(sigma2[k + 1] - sigma2[k]) <= tol * (abs(sigma2[k + 1]) + abs(sigma2[k]))/2) 
        break
    }
    
    # collect results to return
    results <- list("phi0"   = phi0[k + 1],
                    "phi1"   = phi1[k + 1],
                    "sigma2" = sigma2[k + 1])
    if (return_iterates) 
      results <- c(results, list("phi0_iterates"   = phi0,
                                 "phi1_iterates"   = phi1,
                                 "sigma2_iterates" = sigma2,
                                 "nu_iterates"     = nu))
    if (return_condMeanCov) {
      cond <- condMeanCov(y_obs, index_obs, n, n_block, n_in_block, 
                          first_index_in_block, last_index_in_block, previous_obs_before_block, next_obs_after_block, 
                          phi0[k+1], phi1[k+1], sigma2[k+1], full_cov = TRUE)
      results <- c(results, list("cond_mean_y" = cond$mean_y,
                                 "cond_cov_y"  = cond$cov_y))
    }
  }
  
  if(remove_outliers)
    results <- c(results, list("index_outliers" = if(is.null(idx_outliers)) NULL
                                                  else idx_outliers + idx_offset))
  return(results)
}



# extracts the diagonal on top of the main diagonal
diag1 <- function(X) {
  m <- min(dim(X))
  X[1 + dim(X)[1L] + 0L:(m - 2L) * (dim(X)[1L] + 1)]  # main diag: x[1 + 0L:(m - 1L) * (dim(x)[1L] + 1)]
}



# trivial case with no NAs (no need for EM algorithm)
fit_AR1_Gaussian_complete <- function(y, random_walk = FALSE, zero_mean = FALSE) {
  n <- length(y)
  s_y2 <- sum(y[-1])
  s_y1 <- sum(y[-n])
  s_y2y2 <- sum(y[-1]^2) 
  s_y1y1 <- sum(y[-n]^2) 
  s_y2y1 <- sum(y[-1]*y[-n])
  if (!random_walk && !zero_mean) {
    phi1 <- (s_y2y1 - s_y2 * s_y1 / (n - 1)) / (s_y1y1 - s_y1 * s_y1 / (n - 1)) 
    phi0 <- (s_y2 - phi1 * s_y1) / (n - 1)
  } else if (random_walk && !zero_mean){
    phi1 <- 1
    phi0 <- (s_y2 - s_y1) / (n - 1)
  } else if (!random_walk && zero_mean){
    phi1 <- s_y2y1 / s_y1y1 
    phi0 <- 0
  } else{
    phi1 <- 1
    phi0 <- 0
  }
  sigma2 <- (( s_y2y2 + (n - 1) * phi0^2 + phi1^2 * s_y1y1 
               - 2 * phi0 * s_y2 - 2 * phi1 * s_y2y1 + 2 * phi0 * phi1 * s_y1 ) / (n - 1))
  return(list("phi0"   = phi0,
              "phi1"   = phi1,
              "sigma2" = sigma2))
}



find_outliers_AR1_Gaussian <- function(y, index_obs, fitted, p_confidence = 0.95) {
  index_outliers <- NULL
  for (i in 2:length(index_obs)) {
    delta_i <- index_obs[i] - index_obs[i-1]
    mu_expected <- sum(fitted$phi1^( 0:(delta_i - 1) )) * fitted$phi0 + fitted$phi1^delta_i * y[index_obs[i-1]]
    upper <- max(y[index_obs[i]], 2*mu_expected - y[index_obs[i]])
    lower <- min(y[index_obs[i]], 2*mu_expected - y[index_obs[i]])
    # compute the probability of the observation in [lower, upper]
    p <- mvtnorm::pmvt(lower = lower, upper = upper, delta = mu_expected, df = 10, sigma = (10 - 2)/10 * fitted$sigma2) 
    #browser()
    #p <- mvtnorm::pmvnorm(lower = lower, upper = upper, mean = mu_expected, sigma = fitted$sigma2) 
    
    # if p larger than 95%, then the observation lies outside the 95% confidence interval.
    if (p > 0.95 ) {
      index_outliers <- c(index_outliers, index_obs[i])
      #browser()
    }
  }
  return(index_outliers)
}  







#' @title Impute missing values of time series based on a Gaussian AR(1) model
#'
#' @description Impute missing values of time series by drawing samples from 
#'              the conditional distribution of the missing values given the 
#'              observed data based on a Gaussian AR(1) model as estimated 
#'              with the function \code{\link{fit_AR1_Gaussian}}. Outliers 
#'              can be detected and removed.
#' 
#' @inheritParams fit_AR1_Gaussian
#' @param n_samples Positive integer indicating the number of imputations (default is \code{1}).
#' @param return_estimates Logical value indicating if the estimates of the model parameters 
#'                         are to be returned (default is \code{FALSE}).
#'                         
#' @return By default (i.e., for \code{n_samples = 1} and \code{return_estimates = FALSE}), 
#'         the function will return an imputed time series of the same class and dimensions 
#'         as the argument \code{y} with one new attribute recording the locations of missing 
#'         values (the function \code{\link{plot_imputed}} will make use of such information
#'         to indicate the imputed values), as well as locations of outliers removed.
#'         
#'         If \code{n_samples > 1}, the function will return a list consisting of \code{n_sample} 
#'         imputed time series with names: y_imputed.1, y_imputed.2, etc. 
#'         
#'         If \code{return_estimates = TRUE}, in addition to the imputed time series \code{y_imputed}, 
#'         the function will return the estimated model parameters:
#'         \item{\code{phi0}}{The estimate for \code{phi0} (numeric scalar or vector depending 
#'                            on the number of time series).}
#'         \item{\code{phi1}}{The estimate for \code{phi1} (numeric scalar or vector depending 
#'                            on the number of time series).}
#'         \item{\code{sigma2}}{The estimate for \code{sigma2} (numeric scalar or vector depending 
#'                              on the number of time series).}
#' 
#' @author Junyan Liu and Daniel P. Palomar
#' 
#' @seealso \code{\link{plot_imputed}}, \code{\link{fit_AR1_Gaussian}}, \code{\link{impute_AR1_t}}
#' 
#' @references
#' R. J. Little and D. B. Rubin, Statistical Analysis with Missing Data, 2nd ed. Hoboken, N.J.: John Wiley & Sons, 2002.
#' 
#' J. Liu, S. Kumar, and D. P. Palomar, "Parameter estimation of heavy-tailed AR model with missing 
#' data via stochastic EM," IEEE Trans. on Signal Processing, vol. 67, no. 8, pp. 2159-2172, 15 April, 2019. 
#' 
#' @examples
#' library(imputeFin)
#' data(ts_AR1_Gaussian) 
#' y_missing <- ts_AR1_Gaussian$y_missing
#' y_imputed <- impute_AR1_Gaussian(y_missing)
#' plot_imputed(y_imputed)
#' 
#' @import zoo
#' @import MASS
#' @export
impute_AR1_Gaussian <- function(y, n_samples = 1, random_walk = FALSE, zero_mean = FALSE, 
                                remove_outliers = FALSE,
                                return_estimates = FALSE, tol = 1e-10, maxiter = 1000) { 
  # error control
  if (!is.matrix(try(as.matrix(y), silent = TRUE))) stop("\"y\" must be coercible to a vector or matrix.")
  if (round(n_samples)!=n_samples | n_samples<=0) stop("\"n_samples\" must be a positive integer.")

  # manage multiple columns  
  if (NCOL(y) > 1) {
    results_list <- lapply(c(1:NCOL(y)), FUN = function(i) {
      impute_AR1_Gaussian(y[, i, drop = FALSE], n_samples, random_walk, zero_mean, remove_outliers, return_estimates, tol, maxiter)
      })
    names(results_list) <- colnames(y)
    if (n_samples == 1 && !return_estimates) {  # return directly a matrix like y
      results <- do.call(cbind, results_list)
      index_miss_list     <- lapply(results_list, FUN = function(res) attr(res, "index_miss"))
      index_outliers_list <- lapply(results_list, FUN = function(res) attr(res, "index_outliers"))
      attr(results, "index_miss")     <- index_miss_list
      attr(results, "index_outliers") <- index_outliers_list
    } else if (n_samples == 1 && return_estimates) {
      results <- do.call(mapply, c("FUN" = cbind, results_list, "SIMPLIFY" = FALSE))
      index_miss_list     <- lapply(results_list, FUN = function(res) attr(res$y_imputed, "index_miss"))
      index_outliers_list <- lapply(results_list, FUN = function(res) attr(res$y_imputed, "index_outliers"))
      attr(results$y_imputed, "index_miss")     <- index_miss_list
      attr(results$y_imputed, "index_outliers") <- index_outliers_list
    } else {
      results <- do.call(mapply, c("FUN" = cbind, results_list, "SIMPLIFY" = FALSE))
      index_miss_list     <- lapply(results_list, FUN = function(res) attr(res$y_imputed.1, "index_miss"))
      index_outliers_list <- lapply(results_list, FUN = function(res) attr(res$y_imputed.1, "index_outliers"))
      for (i in 1:n_samples) {
        attr(results[[i]], "index_miss")     <- index_miss_list
        attr(results[[i]], "index_outliers") <- index_outliers_list
      }
      if (return_estimates) {
        results$phi0   <- as.vector(results$phi0)
        results$phi1   <- as.vector(results$phi1)
        results$sigma2 <- as.vector(results$sigma2)
      }
    }
    return(results)
  }
  
  #
  #   code for y single-column
  #  
  # error control
  if (!is.numeric(y)) stop("\"y\" only allows numerical or NA values.")
  if (sum(!is.na(y)) < 5) stop("Each time series in \"y\" must have at least 5 observations.")
  
  y_attrib <- attributes(y)
  y <- as.numeric(y)
  y_imputed <- matrix(rep(y, times = n_samples), ncol = n_samples)
  
  if (remove_outliers) {
    fitted <- fit_AR1_Gaussian(y, random_walk, zero_mean, remove_outliers = TRUE, tol = tol, maxiter = maxiter)
    if (!is.null(index_outliers <- fitted$index_outliers))
      y[index_outliers] <- NA
  }  
  
  # imputation
  if (!any_inner_NA(y)) {  # trivial case with no inner NAs: do nothing
    index_miss <- which(is.na(y))
    if (return_estimates && !remove_outliers) 
      fitted <- fit_AR1_Gaussian(y, random_walk, zero_mean, remove_outliers = FALSE, tol = tol, maxiter = maxiter)
  } else {
    fitted <- fit_AR1_Gaussian(y, random_walk, zero_mean, remove_outliers = FALSE, return_condMeanCov = TRUE, tol = tol, maxiter = maxiter)
    index_miss <- which(is.na(y))  # indexes of missing values
    index_obs <- which(!is.na(y))
    index_obs_min <- min(index_obs)
    index_obs_max <- max(index_obs)
    
    index_miss_middle <- index_miss[index_miss>index_obs_min & index_miss<index_obs_max]
    if (length(index_miss_middle) > 0) {
      index_miss_deleted <- index_miss_middle - (index_obs_min - 1)
      y_imputed[index_miss_middle, ] <- t(MASS::mvrnorm(n = n_samples, fitted$cond_mean_y[index_miss_deleted], fitted$cond_cov_y[index_miss_deleted, index_miss_deleted]))
    }
  }
  
  # prepare results
  if (n_samples == 1) {
    attributes(y_imputed) <- y_attrib
    attr(y_imputed, "index_miss") <- index_miss
    if(remove_outliers) attr(y_imputed, "index_outliers") <- index_outliers
    results <- if (!return_estimates) y_imputed else list("y_imputed" = y_imputed)
  } else {
    y_imputed <-lapply(split(y_imputed, col(y_imputed)), FUN = function(x) { attributes(x) <- y_attrib
                                                                             attr(x, "index_miss") <- index_miss
                                                                             if(remove_outliers) attr(x, "index_outliers") <- index_outliers
                                                                             return(x) })
    results <- c("y_imputed" = y_imputed)
  }
  if (return_estimates)  results <- c(results, list("phi0"   = fitted$phi0,
                                                    "phi1"   = fitted$phi1,
                                                    "sigma2" = fitted$sigma2))
  return(results)
}



#
# Compute the conditional mean and covariance matrix of y given the observed data and current estimates
#
condMeanCov <- function(y_obs, index_obs, n, n_block, n_in_block, 
                        first_index_in_block, last_index_in_block, previous_obs_before_block, next_obs_after_block, 
                        phi0, phi1, sigma2, full_cov = FALSE) {
  
  cond_mean_y <- rep(NA, n)  # mean of y conditional on all the observed data
  cond_mean_y[index_obs] <- y_obs
  if (full_cov) cond_cov_y <- matrix(0, nrow = n, ncol = n)  # covariance of y conditional on all the observed data
  cond_cov_y_diag <- rep(0, n)
  cond_cov_y_diag1 <- rep(0, n-1)
  
  for (d in 1:n_block) {  # for each missing block
    n_d <- n_in_block[d]  # number of missing values in the d-th missing block  
    index_d <- first_index_in_block[d]:last_index_in_block[d]  # indexes of missing values in the d-th missing block
    phi1_exp <- phi1^(0:(n_d+1))  # phi1[k] to the power of 0, 1, ..., (n_d+1)
    sum_phi1_exp <- cumsum(phi1_exp) 
    cond_mean_block_obs <- sum_phi1_exp[1:(n_d+1)] * phi0 + phi1_exp[2:(n_d+2)] * previous_obs_before_block[d]  # mean of the d-th missing block conditional on all the previous observed samples
    # computation of cond_cov_block_obs
    cond_cov_block_obs_diag <- c(1, rep(NA, n_d))
    for (i in 1:n_d)  cond_cov_block_obs_diag[i+1] <- cond_cov_block_obs_diag[i]*phi1^2 + 1
    cond_cov_block_obs_diag1 <- cond_cov_block_obs_diag[1:n_d] * phi1
    cond_cov_block_obs_lastcol <- cond_cov_block_obs_diag * rev(phi1_exp[1:(n_d+1)])
    if (full_cov) {
      cond_cov_block_obs <- matrix(nrow = n_d+1, ncol = n_d+1)  # covariance of the d-th missing block and the next observation conditional on all the previous observed samples  
      diag(cond_cov_block_obs) <- cond_cov_block_obs_diag
      for (i in 1:n_d)
        cond_cov_block_obs[i, (i+1):(n_d+1)] <-
        cond_cov_block_obs[(i+1):(n_d+1), i] <- cond_cov_block_obs_diag[i] * phi1_exp[2:(n_d+1 - i+1)]
      #sanity check
      #if (sum(abs(cond_cov_block_obs_lastcol - cond_cov_block_obs[, n_d+1])) > 1e-9 ||
      #    sum(abs(cond_cov_block_obs_diag1 - diag1(cond_cov_block_obs))) > 1e-9 ||
      #    sum(abs(cond_cov_block_obs_diag - diag(cond_cov_block_obs))) > 1e-9) {
      #  message("Error in computation of cond_cov_block_obs...")
      #  browser()
      #}
    }

    # mean of the d-th missing block conditional on all the observed data
    cond_mean_block <- cond_mean_block_obs[1:n_d] + cond_cov_block_obs_lastcol[1:n_d] / cond_cov_block_obs_lastcol[n_d+1] *
                                                    (next_obs_after_block[d] - cond_mean_block_obs[n_d+1])
    # covariance of the d-th missing block conditional on all the observed data
    if (full_cov) cond_cov_block <- sigma2 * (cond_cov_block_obs[1:n_d, 1:n_d] - tcrossprod(cond_cov_block_obs[1:n_d, n_d+1])/cond_cov_block_obs[n_d+1, n_d+1])  #slower: cond_cov_block <- sigma2 * (cond_cov_block_obs[1:n_d, 1:n_d] - cond_cov_block_obs[1:n_d, n_d+1] %*% t(cond_cov_block_obs[1:n_d, n_d+1])/cond_cov_block_obs[n_d+1, n_d+1])
    cond_cov_block_diag <- sigma2 * (cond_cov_block_obs_diag[-(n_d+1)] - cond_cov_block_obs_lastcol[1:n_d]^2/cond_cov_block_obs_lastcol[n_d+1])
    cond_cov_block_diag1 <- sigma2 * (cond_cov_block_obs_diag1[-n_d] - cond_cov_block_obs_lastcol[1:(n_d-1)]*cond_cov_block_obs_lastcol[2:n_d]/cond_cov_block_obs_lastcol[n_d+1])

    # update global mean and cov matrix
    cond_mean_y[index_d] <- cond_mean_block
    if (full_cov) cond_cov_y[index_d, index_d] <- cond_cov_block
    cond_cov_y_diag[index_d] <- cond_cov_block_diag  # cond_cov_y_diag[index_d] <- diag(cond_cov_block)
    cond_cov_y_diag1[index_d[-n_d]] <- cond_cov_block_diag1  # cond_cov_y_diag1[index_d[-n_d]] <- diag1(cond_cov_block)
  }
  # sanity check
  #if (full_cov)
  #  if (sum(abs(diag(cond_cov_y) - cond_cov_y_diag)) > 1e-9 ||
  #      sum(abs(diag1(cond_cov_y) - cond_cov_y_diag1)) > 1e-9) {
  #    message("Error in computation of cond_cov_y...")
  #    browser()
  #  }
  
  results <- list("mean_y"      = cond_mean_y,
                  "cov_y_diag"  = cond_cov_y_diag,
                  "cov_y_diag1" = cond_cov_y_diag1)
  if (full_cov) 
    results <- c(results, list("cov_y" = cond_cov_y))
  return(results)
}



#
# A heuristic method to compute the parameters of Gaussian AR(1) model from incomplete time series
#
fit_AR1_Gaussian_heuristic <- function(y, index_miss, random_walk = FALSE, zero_mean = TRUE) {
  
  index_miss_p <- c(0, index_miss, length(y) + 1)
  delta_index_miss_p <- diff(index_miss_p)
  index_delta_index_miss_p <- which(delta_index_miss_p > 2)
  n_obs_block <- length(index_delta_index_miss_p)  # number of observation blocks with more than 1 sample
  n_in_obs_block <- delta_index_miss_p[index_delta_index_miss_p] - 1  # number of observed samples in each qualified observation block
  
  m <- 0
  y_obs2 <- y_obs1 <- c()
  for (i in 1:n_obs_block) {
    y_obs1[(m + 1):(m + n_in_obs_block[i] - 1)] <- y[(index_miss_p[index_delta_index_miss_p[i]] + 1):(index_miss_p[index_delta_index_miss_p[i] + 1] - 2)]
    y_obs2[(m + 1):(m + n_in_obs_block[i] - 1)] <- y[(index_miss_p[index_delta_index_miss_p[i]] + 2):(index_miss_p[index_delta_index_miss_p[i] + 1] - 1)]
    m <- m + n_in_obs_block[i] - 1
  }
  n_y_obs1 <- length(y_obs1)
  
  s_y_obs1 <- sum(y_obs1)
  s_y_obs2 <- sum(y_obs2)
  s_y_obs1_square <- sum(y_obs1^2)
  s_y_obs2_square <- sum(y_obs2^2)
  s_y_obs2_y_obs1 <- sum(y_obs1 * y_obs2)
  
  # compute the estimates
  if (!random_walk && !zero_mean) {
    phi1 <- (s_y_obs2_y_obs1 - s_y_obs2 * s_y_obs1 / n_y_obs1) / (s_y_obs1_square - s_y_obs1 * s_y_obs1 / n_y_obs1) 
    phi0 <- (s_y_obs2 - phi1 * s_y_obs1) / n_y_obs1
  } else if (random_walk && !zero_mean){
    phi1 <- 1
    phi0 <- (s_y_obs2 - phi1 * s_y_obs1) / n_y_obs1
  } else if (!random_walk && zero_mean){
    phi1 <- s_y_obs2_y_obs1 / s_y_obs1_square 
    phi0 <- 0
  } else{
    phi1 <- 1
    phi0 <- 0
  }
  sigma2 <- (( s_y_obs2_square + n_y_obs1 * phi0^2 + phi1^2 * s_y_obs1_square 
               - 2 * phi0 * s_y_obs2 - 2 * phi1 * s_y_obs2_y_obs1 + 2 * phi0 * phi1 * s_y_obs1 ) / n_y_obs1)
  
  return(list("phi0"   = phi0,
              "phi1"   = phi1,
              "sigma2" = sigma2))
}



#
# find the missing blocks in a time series with missing values
#
findMissingBlock <- function(y) {
  n <- length(y)  # length of the time series
  index_obs <- which(!is.na(y))  # indexes of observed values
  index_miss <- setdiff(1:n, index_obs)  # indexes of missing values
  n_obs <- length(index_obs)
  y_obs <- y[index_obs]  # observed values
  delta_index_obs <- diff(index_obs)
  index_delta_index_obs <- which(delta_index_obs > 1)
  n_block <- length(index_delta_index_obs)  # number of missing blocks
  n_in_block <- delta_index_obs[index_delta_index_obs] - 1  # number of missing values in each block
  first_index_in_block <- index_obs[index_delta_index_obs] + 1  # index of the first missing value in each block
  last_index_in_block <- index_obs[index_delta_index_obs] + n_in_block  # index of the last missing value in each block
  previous_obs_before_block <- y[first_index_in_block - 1]  # previous observed value before each block
  next_obs_after_block <- y[last_index_in_block + 1]  # next observed value after each block
  
  return(list("n"                         = n,
              "index_obs"                 = index_obs,
              "index_miss"                = index_miss,
              "n_obs"                     = n_obs,
              "y_obs"                     = y_obs,
              "delta_index_obs"           = delta_index_obs,
              "n_block"                   = n_block,
              "n_in_block"                = n_in_block,
              "first_index_in_block"      = first_index_in_block,
              "last_index_in_block"       = last_index_in_block,
              "previous_obs_before_block" = previous_obs_before_block,
              "next_obs_after_block"      = next_obs_after_block))
}

