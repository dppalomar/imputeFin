#' @title Fit Student's t AR(1) Models to Time Series with Missing values
#'
#' @description Estimate the parameters of a Student's t AR(1) model from a time series with missing values
#'
#' @param y a xts object indicating time series with missing values. The first and last one should not be NA.
#' @param random_walk logical. If TRUE, y is a random walk time series, and phi1 = 1. If FALSE, y is a general AR(1) time series, and phi1 is unknown. The default value is FALSE.
#' @param zero_mean logical. If TRUE, y is a zero-mean time series, and phi0 = 1. If FALSE, y is a general AR(1) time series, and phi0 is unknown. The default value is FALSE.
#' @param n_chain a positive integer indicating the number of the parallel Markov chains used (default \code{10}).
#' @param n_thin  a positive integer indicating the sampling period of the Gibbs sampling (default \code{1}). Every n_thin-th samples is used. This is aimed to reduce the dependence of the samples.
#' @param n_iter a positive integer indicating the number of the iterations (default \code{100}).
#' @param K a positive number controlling the values the step sizes (default \code{30}).
#' @param iterates logical. If TRUE, then the iterates are outputted. If FALSE, they are ignored. The default value is FALSE.
#' @return A list containing the following elements:
#' \item{\code{phi0}}{real number, the estimate for phi0}
#' \item{\code{phi1}}{real number, the estimate for phi1}
#' \item{\code{sigma2}}{positive number, the estimate for sigma^2}
#' \item{\code{nu}}{positive number, the estimate for nu}
#' \item{\code{phi0_iterate}}{a vector of real numbers, the estimates for phi0 in each iteration}
#' \item{\code{phi1_iterate}}{a vector of real numbers, the estimates for phi1 in each iteration}
#' \item{\code{sigma2_iterate}}{a vector of positive numbers, the estimates for sigma^2 in each iteration}
#' \item{\code{nu_iterate}}{a vector of positive numbers, the estimates for nu in each iteration}
#' @author Junyan Liu and Daniel P. Palomar
#' @examples 
#' library(imputeFin)
#' library(xts)
#' # generate a complete Student's t AR(1) time series
#' phi0 <- 0
#' phi1 <- 1
#' sigma2 <- 0.01 
#' nu <- 1
#' n <- 200
#' n_miss <- 25 
#' n_drop <- 100
#' n_total <- n + n_drop
#' data <- vector(length = n_total)
#' epsilon <- vector(length = n_total - 1)# innovations
#' data[1] <- 0
#' for (i in 2:n_total) {
#'   epsilon[i-1] <- rt(1, nu) * sqrt(sigma2)
#'   data[i] <- phi0 + phi1 * data[i-1] + epsilon[i-1]
#' }
#' data <- data[(n_drop + 1):n_total] # drop the first n_drop to reduce the influence of initial point
#' dates <- seq(as.Date("2016-01-01"), length = n, by = "days") 
#' y_orig <- xts(data,  dates)
#' 
#' # creat missing values
#' index_miss <- sample( 2:(n - 1), n_miss, FALSE)
#' index_miss <- sort(index_miss)
#' y <- y_orig
#' y[index_miss] <- NA
#' 
#' # estimate the parameters from this incomplete time series
#' estimation_result <- estimateAR1t(y)
#' @import xts
#' @export
estimateAR1t <- function(y, random_walk = FALSE, zero_mean = FALSE, 
                         iterates = FALSE, condMean_Gaussian = TRUE,
                         n_chain = 10, n_thin = 1, n_iter = 100, K = 30) {
  
  if (NCOL(y) > 1){
    estimation_list <- apply(y, MARGIN = 2, FUN = estimateAR1t, random_walk, zero_mean, iterates, condMean_Gaussian, n_chain, n_thin, n_iter, K)
    phi0 <- unlist(lapply(estimation_list, function(x){x$phi0}))
    phi1 <- unlist(lapply(estimation_list, function(x){x$phi1}))
    sigma2 <- unlist(lapply(estimation_list, function(x){x$sigma2}))
    nu <- unlist(lapply(estimation_list, function(x){x$nu}))
    return(c(estimation_list, list("phi0" = phi0,
                                   "phi1" = phi1,
                                   "sigma2" = sigma2,
                                   "nu" = nu)))
  }
  
  y <- as.numeric(y)
# trivial case with no NAs
  if (!anyNA(y)) return(estimateAR1tComplete(y, random_walk, zero_mean, iterates))
  
# find the missing blocks
  list2env(findMissingBlock(y), envir = environment())
  
# initialize the estimates and some parameters
  phi0 <- phi1 <- sigma2 <- nu <- gamma <- c()
  estimation_Gaussian <- estimateAR1Gaussian(y, random_walk, zero_mean, condMeanCov = TRUE)
  phi0[1] <- estimation_Gaussian$phi0
  phi1[1] <- estimation_Gaussian$phi1
  sigma2[1] <- estimation_Gaussian$sigma2
  nu[1] <- 3
  y_samples <- matrix(estimation_Gaussian$cond_mean_y, n, n_chain)
  tau_samples <- matrix(NA, n, n_chain)
  s <- s_approx <- rep(0, 7)  # approximations of the sufficient statistics
  
  for (k in 1:n_iter) {
    # draw realizations of the missing values from their posterior distribution
    for (j in 1:n_chain) {
      sample <- samplingLatentVariables(y_sample_init = y_samples[, j], n_thin, n_block, n_in_block,
                                        first_index_in_block, last_index_in_block, previous_obs_before_block, next_obs_after_block,
                                        phi0[k], phi1[k], sigma2[k], nu[k])
      y_samples[, j] <- sample$y
      tau_samples[, j] <- sample$tau
    }
    # # Daniel: implement loop above in parallel
    # lmd <- t(sapply(1:n_chain, FUN = function(k) {
    #   Ut <- eigen(cov(X[-c(t0[k]:t1[k]), ]), symmetric = TRUE)$vectors
    #   colMeans((X[c(t0[k]:t1[k]), , drop = FALSE] %*% Ut)^2)
    # }))
    
    # approximate the sufficient statistics
    if (k <= K) {
      gamma[k] <- 1
    } else
      gamma[k] <- 1/(k - K)
    
    s[1] <- sum(log(tau_samples[2:n,]) - tau_samples[2:n,]) / n_chain
    s[2] <- sum(tau_samples[2:n,] * y_samples[2:n,]^2) / n_chain
    s[3] <- sum(tau_samples[2:n,] ) / n_chain
    s[4] <- sum(tau_samples[2:n,] * y_samples[1:(n-1),]^2) / n_chain
    s[5] <- sum(tau_samples[2:n,] * y_samples[2:n,]) / n_chain
    s[6] <- sum(tau_samples[2:n,] * y_samples[2:n,] * y_samples[1:(n - 1),]) / n_chain
    s[7] <- sum(tau_samples[2:n,] * y_samples[1:(n-1),]) / n_chain
    s_approx <- s_approx + gamma[k] * (s - s_approx)
    
    # update the estimates
    if (!random_walk && !zero_mean) {
      phi1[k+1] <- ( s_approx[3] * s_approx[6] - s_approx[5] * s_approx[7] ) / ( s_approx[3] * s_approx[4] -  s_approx[7]^2 )
      phi0[k+1] <- (s_approx[5] - phi1[k+1] * s_approx[7] ) / s_approx[3]
    } else if (random_walk && !zero_mean){
      phi1[k+1] <- 1
      phi0[k+1] <- (s_approx[5] -  s_approx[7] ) / s_approx[3]
    } else if (!random_walk && zero_mean){
      phi1[k+1] <- s_approx[6] / s_approx[4] 
      phi0[k+1] <- 0
    } else{
      phi1[k+1] <- 1
      phi0[k+1] <- 0
    }
    sigma2[k+1] <- (s_approx[2] + phi0[k+1]^2 * s_approx[3] + phi1[k+1]^2 * s_approx[4] - 2 * phi0[k+1] * s_approx[5]
                    - 2 * phi1[k+1] * s_approx[6] + 2 * phi0[k+1] * phi1[k+1] * s_approx[7]) / (n - 1)
    
    f_nu <- function(nu, n, s_approx1)
      return(-sum(0.5 * nu * s_approx1 +  (0.5 * nu * log(0.5 * nu) - lgamma(0.5 * nu)) * (n - 1)))
    
    optimation_result <- optimize(f_nu, c(1e-6, 1e6), n, s_approx[1])
    nu[k + 1] <- optimation_result$minimum
    
  }
  
  results <- list("phi0" = phi0[k + 1],
                  "phi1" = phi1[k + 1],
                  "sigma2" = sigma2[k + 1],
                  "nu" = nu[k + 1])
  if (iterates) 
    results <- c(results, list("phi0_iterate" = phi0,
                               "phi1_iterate" = phi1,
                               "sigma2_iterate" = sigma2,
                               "nu_iterate" = nu))
  if(condMean_Gaussian)
    results <- c(results, list("cond_mean_Gaussian" = estimation_Gaussian$cond_mean_y))
  
  return(results) 
  
}



#' @title Imputate Missing Values in  Incomplete Student's t AR(1) Time Series 
#'
#' @description Estimate the parameters of the Student's t AR(1) model from a time series with missing values and impute the missing values based on the estimates
#'
#' @param y a xts object indicating time series with missing values. The first and last one should not be NA.
#' @param n_sample a positive integer indicating the number of imputations (default \code{1}).
#' @param param a list consisting of the paramters of the Student's t AR(1) time series y if known. The default value is FALSE.
#' @param random_walk logical. If TRUE, y is a random walk time series, and phi1 = 1. If FALSE, y is a general AR(1) time series, and phi1 is unknown. The default value is FALSE.
#' @param zero_mean logical. If TRUE, y is a zero-mean time series, and phi0 = 1. If FALSE, y is a general AR(1) time series, and phi0 is unknown.
#' @param n_burn a positive integer controlling the length of the burn-in period of the Gibb sampling (default \code{100}). The first (n_burn * n_thin) samples generated will be ignored.
#' @param n_thin a positive integer indicating the sampling period of Gibbs sampling. After then burn-in perid (default \code{50}), every n_thin-th sample is used. This is aimed to reduce the dependence of the samples.
#' @return 
#' \item{\code{y_imputed}  }{a xts object, each column is a imputed complete time series}
#' @author Junyan Liu and Daniel P. Palomar
#' @examples
#' library(imputeFin)
#' library(xts)
#' phi0 <- 1
#' phi1 <- 0.5 
#' sigma2 <- 0.01 
#' nu <- 1
#' n <- 200
#' n_miss <- 25 
#' n_drop <- 100
#' n_total <- n + n_drop
#' data <- vector(length = n_total)
#' epsilon <- vector(length = n_total - 1)# innovations
#' data[1] <- 0
#' for (i in 2:n_total) {
#'   epsilon[i - 1] <- rt(1, nu) * sqrt(sigma2)
#'   data[i] <- phi0 + phi1 * data[i - 1] + epsilon[i - 1]
#' }
#' data <- data[(n_drop + 1):n_total] # drop the first n_drop to reduce the influence of initial point
#' dates <- seq(as.Date("2016-01-01"), length = n, by = "days") 
#' y_orig <- xts(data,  dates)
#' 
#' # creat missing values
#' index_miss <- sample( 2:(n - 1), n_miss, FALSE)
#' index_miss <- sort(index_miss)
#' y <- y_orig
#' y[index_miss] <- NA
#' 
#' # impute the missing values and generate n_sample complete time series
#' y_imputed <- imputeAR1t( y_miss, n_sample = 3) # if the parameters are unknown
#' param = list("phi0" = phi0,
#'              "phi1" = phi1,
#'              "sigma2" = sigma2,
#'              "nu" = nu)
#' y_imputed <- imputeAR1t(y_miss, n_sample = 3, param) # if the parameters are unknown
#' @import  xts
#' @export
imputeAR1t <- function(y, n_samples = 1, random_walk = FALSE, zero_mean = FALSE,
                       n_burn = 100, n_thin = 50,
                       estimates = FALSE) {
  
  if (NCOL(y) > 1) {
    results_list <- lapply(c(1:NCOL(y)), FUN = function(i){imputeAR1t(y[, i], n_samples, random_walk, zero_mean, n_burn, n_thin, estimates)})
    if (n_samples == 1 && !estimates) {
      index_miss_list <- lapply(results_list, FUN = function(result){attributes(result)$index_miss})
      results <- do.call(cbind, results_list)
      attr(results, "index_miss") = index_miss_list
    } else if (n_samples == 1 && estimates) {
      index_miss_list <- lapply(results_list, FUN = function(result){attributes(result$y_imputed)$index_miss})
      results <- do.call(mapply, c("FUN" = cbind, results_list, "SIMPLIFY" = FALSE))
      attr(results$y_imputed, "index_miss") = index_miss_list
    } else {
      index_miss_list <- lapply(results_list, FUN = function(result){attributes(result$y_imputed.1)$index_miss})
      results <- do.call(mapply, c("FUN" = cbind, results_list, "SIMPLIFY" = FALSE))
      for (i in 1:n_samples) {
        attr(results[[i]], "index_miss") = index_miss_list  
      }
      if (estimates) {
        results$phi0 <- as.vector(results$phi0)
        results$phi1 <- as.vector(results$phi1)
        results$sigma2 <- as.vector(results$sigma2)
        results$nu <- as.vector(results$nu)
      }
    }
    return(results)
  }  
  
##############################################################
  y_attrib <- attributes(y)
  y <- as.numeric(y)
  
  # trivial case with no NAs
  if (!anyNA(y)){
    y_imputed <- matrix(rep(y, times = n_samples), ncol = n_samples)
    if (estimates) estimation_result <- estimateAR1t(y, random_walk, zero_mean)
    index_miss = NULL
  } else {
    estimation_result <- estimateAR1t(y, random_walk, zero_mean, condMean_Gaussian = TRUE)
    list2env(findMissingBlock(y), envir = environment())
    y_tmp <- estimation_result$cond_mean_Gaussian
    phi0 <- estimation_result$phi0
    phi1 <- estimation_result$phi1
    sigma2 <- estimation_result$sigma2
    nu <- estimation_result$nu
    y_imputed <- matrix(nrow = n, ncol = n_samples)

    # browser()
    
    # burn-in period
    for (i in 1:n_burn) {
      sample <- samplingLatentVariables(y_tmp, n_thin = 1, n_block, n_in_block,
                                        first_index_in_block, last_index_in_block, previous_obs_before_block, next_obs_after_block,
                                        phi0, phi1, sigma2, nu) 
      y_tmp <- sample$y
    }
    # sample every n_thin-th sample
    for (j in 1:n_samples) {
      sample <- samplingLatentVariables(y_tmp, n_thin, n_block, n_in_block,
                                        first_index_in_block, last_index_in_block, previous_obs_before_block, next_obs_after_block,
                                        phi0, phi1, sigma2, nu) 
      y_imputed[, j] <- sample$y
    }
  }
  
  if (n_samples == 1) {
    attributes(y_imputed) <- y_attrib
    attr(y_imputed, "index_miss") <- index_miss
    if (!estimates) {
      results <- y_imputed
    } else
      results <- list("y_imputed" = y_imputed)
  } else {
    y_imputed <-lapply(split(y_imputed, col(y_imputed)), FUN = function(x){attributes(x) <- y_attrib
                                                                           attr(x, "index_miss") <- index_miss
                                                                           return(x)})
    results <- c("y_imputed" = y_imputed)
  }
  
  if (estimates)  results <- c(results, list("phi0" = estimation_result$phi0,
                                             "phi1" = estimation_result$phi1,
                                             "sigma2" = estimation_result$sigma2,
                                             "nu" = estimation_result$nu))
  return(results)
}


# Sampling the latent variables y and tau from the conditional distribution via Gibbs sampling.
#   y_sample_init: a vector of real numbers indicating the initialization for y.
#   n_burn:  a positive integer indicating controling the length of the burn-in period of the Gibb sampling. The first (n_burn * n_thin) samples generated will be ignored.
#   n_thin:  a positive integer indicating the sampling period of Gibbs sampling. After then burn-in perid, every n_thin-th sample is used. This is aimed to reduce the dependence of the samples.
#   n_block: a positive integer indicating the number of the missing blocks.
#   n_in_block: a positive integer indicating the numbers of the missing values in each block.
#   first_index_in_block: a vector of positive integers consisting of  the indexes of the first missing value in each missing block.
#   last_index_in_block: a vector of positive integers consisting of  the indexes of the last missing value in each missing block.
#   previous_obs_before_block: a vector of real numbers consisting of the observed values before each missing block.
#   next_obs_after_block: a vector of real numbers consisting of the observed values after each missing block.
#   phi0: a real number, a parameter of the student's t AR(1) model.
#   phi1: a real number, a parameter of the student's t AR(1) model.
#   sigma2: a positive number, a parameter of the student's t AR(1) model.
#   nu: a positive number indicating the degree of freedom, a parameter of the student's t AR(1) model.


samplingLatentVariables <- function( y_sample_init, n_thin, n_block, n_in_block,
                                    first_index_in_block, last_index_in_block, previous_obs_before_block, next_obs_after_block,
                                    phi0, phi1, sigma2, nu) {
  n <- length(y_sample_init)
  tau_tmp <- vector( length = n )
  y_tmp <-  y_sample_init
  
  max_n_in_block <- max( n_in_block ) # maximum number of missing values in the blocks
  phi1_exponential <- phi1^( 0:(max_n_in_block + 1) )
  sum_phi1_exponential <- cumsum( phi1_exponential )
  
  for(j in 1:n_thin){
    # sample tau
    for (i_tau in 2:n) {
      tau_tmp[i_tau] <- rgamma(n = 1, shape = 0.5 * nu + 0.5,
                               rate = 0.5 * ( (y_tmp[i_tau] - phi0 - phi1 * y_tmp[i_tau - 1])^2 / sigma2 + nu) )
    }
    # sample y
    for (d in 1:n_block )  {
      n_in_d_block <- n_in_block[d] # number of missing values in the d-th missing block
      mu_cd <- ( sum_phi1_exponential[1:(n_in_d_block + 1)] * phi0 
                 + phi1_exponential[2:(n_in_d_block + 2)] *  previous_obs_before_block[d] )
      mu1 <- mu_cd[1:n_in_d_block]
      mu2 <- mu_cd[n_in_d_block + 1]
      
      sigma_cd <- matrix( nrow = n_in_d_block + 1, ncol = n_in_d_block + 1)
      
      for(i in 1 : (n_in_d_block + 1) ){
        if(i == 1){ sigma_cd[1, 1] <- sigma2/tau_tmp[ first_index_in_block[d] ]
        } else {
          sigma_cd[i, i] <- sigma_cd[ i - 1, i - 1 ] * phi1^2 + sigma2/tau_tmp[ first_index_in_block[d] + i - 1]
        }
        
        if( i != n_in_d_block + 1){
          sigma_cd[ i, (i + 1) : (n_in_d_block + 1)] <- sigma_cd[ i, i ] * phi1_exponential[ 2:(n_in_d_block + 1 - i + 1) ]
          sigma_cd[ (i + 1) : (n_in_d_block + 1), i ] <- sigma_cd[ i, i ] * phi1_exponential[ 2:(n_in_d_block + 1 - i + 1) ]
        }
      }
      sigma11 <- sigma_cd[ 1 : n_in_d_block, 1 : n_in_d_block]
      sigma12 <- sigma_cd[ 1 : n_in_d_block, n_in_d_block + 1]
      sigma22 <- sigma_cd[ n_in_d_block + 1, n_in_d_block + 1]
      
      mu_d <- mu1 + sigma12 / sigma22  * ( next_obs_after_block[d] - mu2 )
      sigma_d <- sigma11 - sigma12 %*% t( sigma12 )/sigma22  
      
      y_tmp[ first_index_in_block[d] : last_index_in_block[d]] <- MASS::mvrnorm( n = 1, mu = mu_d, Sigma = sigma_d )
    }
  }
  return(list("y" = y_tmp, 
              "tau" = tau_tmp))
}




estimateAR1tComplete <- function(y, random_walk = FALSE, zero_mean = FALSE, 
                                 iterates = FALSE,
                                 tol = 1e-10,  maxiter = 1000) {
 
  phi0 <- phi1 <- sigma2 <- nu <- c()  
  estimation_Gaussian <- estimateAR1Gaussian(y, random_walk, zero_mean, condMeanCov = TRUE)
  phi0[1] <- estimation_Gaussian$phi0
  phi1[1] <- estimation_Gaussian$phi1
  sigma2[1] <- estimation_Gaussian$sigma2
  nu[1] <- 3
  n <- length(y)
  tmp <- (y[-1] - phi0[1] - phi1[1] * y[-n])^2/sigma2[1]
  exp_tau <- vector( length = n )
  
  if (iterates) {
    f = vector()
    f[1] = sum( log(  gamma( 0.5 * (nu[1] + 1) )/gamma( 0.5 * nu[1] )/sqrt( pi * nu[1] * sigma2[1] ) )
                + - 0.5 * (nu[1] + 1) * log( (y[-1] - phi0[1] - phi1[1] * y[-n])^2/sigma2[1]/nu[1] + 1 ) ) 
  }

  
  for ( k in 1:maxiter) {
    exp_tau = (nu[k] + 1)/( nu[k] + tmp )
    s_tau = sum( exp_tau )
    s_tau_y2 = sum( exp_tau * y[-1] )
    s_tau_y1 = sum( exp_tau * y[-n] )
    s_tau_y1y2 = sum( exp_tau * y[-n] * y[-1] )
    s_tau_y1y1 = sum( exp_tau * y[-n] * y[-n] )
  
    if (!random_walk && !zero_mean) {
      phi1[k+1] <- (s_tau * s_tau_y1y2 - s_tau_y2 * s_tau_y1 )/(s_tau * s_tau_y1y1 - s_tau_y1^2)
      phi0[k+1] <- (s_tau_y2 - phi1[k+1] * s_tau_y1)/s_tau
    } else if (random_walk && !zero_mean){
      phi1[k+1] <- 1
      phi0[k+1] <- (s_tau_y2 - s_tau_y1)/s_tau
    } else if (!random_walk && zero_mean){
      phi1[k+1] <- s_tau_y1y2 / s_tau_y1y1 
      phi0[k+1] <- 0
    } else{
      phi1[k+1] <- 1
      phi0[k+1] <- 0
    }
    sigma2[k+1] = sum( exp_tau * (y[-1] - phi0[k+1] - phi1[k+1] * y[-n])^2 )/(n - 1)
    tmp = (y[-1] - phi0[k+1] - phi1[k+1] * y[-n])^2/sigma2[k+1]
    
    # minus log-likelihood about nu,  with mu and sigma fixed as mu[k+1] and sigma[k+1]
    f_nu = function( nu){
      f_nu = - sum ( - 0.5 * (nu + 1) * log( nu + tmp)
                     + lgamma ( 0.5*(nu + 1) ) - lgamma (0.5*nu) + 0.5 * nu * log(nu ) )
      return( f_nu )
    }
    
    opt_rst = optimise ( f_nu, c(1e-6, 1e6) )
    nu[k+1] = opt_rst$minimum
    
#    g[k] = sum( log(  gamma( 0.5 * (nu[k] + 1) )/gamma( 0.5 * nu[k] )/sqrt( pi * nu[k] * sigma[k+1] ) )
#                + - 0.5 * (nu[k] + 1) * log( tmp/nu[k] + 1 ) )   
    if (iterates) f[k+1] = sum( log(  gamma( 0.5 * (nu[k+1] + 1) )/gamma( 0.5 * nu[k+1] )/sqrt( pi * nu[k+1] * sigma2[k+1] ) )
                   - 0.5 * (nu[k+1] + 1) * log( tmp/nu[k+1] + 1 ) ) 
    
    if (abs(phi0[k + 1] - phi0[k]) <= tol * (abs(phi0[k + 1]) + abs(phi0[k]))/2
        && abs(phi1[k + 1] - phi1[k]) <= tol * (abs(phi1[k + 1]) + abs(phi1[k]))/2
        && abs(sigma2[k + 1] - sigma2[k]) <= tol * (abs(sigma2[k + 1]) + abs(sigma2[k]))/2
        && abs(nu[k + 1] - nu[k]) <= tol * (abs(nu[k + 1]) + abs(nu[k]))/2) 
    break
    
  }
  
  results <- list("phi0" = phi0[k+1], 
                  "phi1" = phi1[k+1], 
                  "sigma2" = sigma2[k+1], 
                  "nu" = nu[k+1])
 if (iterates) 
    results <- c(results, list("phi0_iterate" = phi0,
                               "phi1_iterate" = phi1,
                               "sigma2_iterate" = sigma2,
                               "nu_iterate" = nu,
                               "f_iterate" = f))
  return(results)
}


