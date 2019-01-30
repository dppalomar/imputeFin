#' @title Fit Student's t AR(1) Models to Time Series with Missing values
#'
#' @description Estimate the parameters of a Student's t AR(1) model from a time series with missing values
#'
#' @param y a xts object indicating time series with missing values. The first and last one should not be NA.
#' @param n_chain a positive integer indicating the number of the parallel Markov chains used.
#' @param n_thin  a positive integer indicating the sampling period of the Gibbs sampling. Every n_thin-th samples is used. This is aimed to reduce the dependence of the samples.
#' @param n_iteration a positive integer indicating the number of the iterations.
#' @param K a positive number controlling the values the step sizes.
#' @param estimates_init a vector indicating the initialization of the estimates phi0, phi1, sigma^2, and nu.
#' @param y_sample_init a vector indicating the initialization of complete y for Gibbs sampling. It has the same observed values with y while the missing values are imputed.
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
#' library(MASS)
#' library(imputeTS)
#' 
#' #generate a complete Student's t AR(1) time series
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
#'   epsilon[i-1] <- rt(1, nu) * sqrt(sigma2)
#'   data[i] <- phi0 + phi1 * data[i-1] + epsilon[i-1]
#' }
#' data <- data[(n_drop + 1):n_total] # drop the first n_drop to reduce the influence of initial point
#' dates <- seq(as.Date("2016-01-01"), length = n, by = "days") 
#' y_orig <- xts(data,  dates)
#' 
#' #creat missing values
#' index_miss <- sample( 2:(n - 1), n_miss, FALSE)
#' index_miss <- sort(index_miss)
#' y_miss <- y_orig
#' y_miss[index_miss] <- NA
#' 
#' # estimate the parameters from this incomplete time series
#' estimation_result <- estimateAR1(y_miss)
#' 
#' @export
estimateAR1 <- function(y, n_chain = 10, n_thin = 1, n_iteration = 100, K = 30,
                        estimates_init = NULL, y_sample_init = NULL) {
 # find the missing blocks
  n <- length(y)                 # length of the time series
  index_obs <- which(!is.na(y))  # indexes of observed values
  delta_index_obs <- diff(index_obs)
  index_tmp <- which(delta_index_obs > 1)
  n_block <- length(index_tmp)                             # number of missing blocks
  n_in_block <- delta_index_obs[index_tmp] - 1             # number of missing values in each block
  first_index_in_block <- index_obs[index_tmp] + 1         # index of the first missing value in each block
  last_index_in_block <- index_obs[index_tmp] + n_in_block # index of the last missing value in each block
  previous_obs_before_block <- as.numeric( y[first_index_in_block - 1] ) # previous observed value before each block
  next_obs_after_block <- as.numeric(y[last_index_in_block + 1])       # next observed value after each block
  
 # initialize the estimates and some parameters
  phi0 <- phi1 <- sigma2 <- nu <- gamma <- c()
  if(is.null(estimates_init)) {  # if no input for estimates_init, use Gaussian AR(1) estimates as the initialization
    Gaussian_estimates <- stats::arima(y, order = c(1, 0, 0), include.mean = TRUE)
    phi1[1] <- Gaussian_estimates$coef["ar1"]
    phi0[1] <- Gaussian_estimates$coef["intercept"] * (1-phi1[1])
    sigma2[1] <- Gaussian_estimates$sigma2
    nu[1] <- 3
  } else {
    phi0[1] <- estimates_init[["phi0"]]
    phi1[1] <- estimates_init[["phi1"]]
    sigma2[1] <- estimates_init[["sigma2"]]
    nu[1] <- estimates_init[["nu"]]
  }

  s <- s_approx <- rep(0, 7)  # approximations of the sufficient statistics
  if(is.null(y_sample_init))  # if no input for y_sample_init, use imputed result generated from package imputeTS as initialization
    y_sample_init <- imputeTS::na.kalman(y)

  y_samples <- matrix(y_sample_init, n, n_chain)
  tau_samples <- matrix(NA, n, n_chain)
  for (k in 1:n_iteration) {
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
    if (k <= K)
      gamma[k] <- 1
    else
      gamma[k] <- 1/(k - K)
    s[1] <- sum(log(tau_samples[2:n,]) - tau_samples[2:n,] )/n_chain
    s[2] <- sum(tau_samples[2:n,] * y_samples[2:n,]^2)/n_chain
    s[3] <- sum(tau_samples[2:n,] )/n_chain
    s[4] <- sum(tau_samples[2:n,] * y_samples[1:(n-1),]^2)/n_chain
    s[5] <- sum(tau_samples[2:n,] * y_samples[2:n,])/n_chain
    s[6] <- sum(tau_samples[2:n,] * y_samples[2:n,] * y_samples[1:(n - 1),])/n_chain
    s[7] <- sum(tau_samples[2:n,] * y_samples[1:(n-1),])/n_chain
    s_approx <- s_approx + gamma[k] * (s - s_approx)

    # update the estimates
    phi1[k+1] <- ( s_approx[3] * s_approx[6] - s_approx[5] * s_approx[7] ) / ( s_approx[3] * s_approx[4] -  s_approx[7]^2 )
    phi0[k+1] <- (s_approx[5] - phi1[k+1] * s_approx[7] ) / s_approx[3]
    sigma2[k+1] <- (s_approx[2] + phi0[k+1]^2 * s_approx[3] + phi1[k+1]^2 * s_approx[4] - 2 * phi0[k+1] * s_approx[5]
                    - 2 * phi1[k+1] * s_approx[6] + 2 * phi0[k+1] * phi1[k+1] * s_approx[7]) / (n-1)
    
    fn_surrogate_nu <- function(nu, n, s_approx1)
      return(-sum(0.5*nu*s_approx1 +  (0.5*nu*log(0.5*nu) - lgamma(0.5*nu)) * (n - 1)))
  
    optimation_result <- optimize(fn_surrogate_nu, c(1e-6, 1e6), n, s_approx[1])
    nu[k+1] <- optimation_result$minimum
    
    # # check convergence on parameters and objective function
    # werr <- sum(abs(w_next - wk)) / max(1, sum(abs(wk)))
    # ferr <- abs(fun_next - fun_k) / max(1, abs(fun_k))
    # if (k > 1 && (werr < wtol || ferr < ftol))
    #   break    
  }
  
  return(list("phi0" = phi0[k + 1],
              "phi1" = phi1[k + 1],
              "sigma2" = sigma2[k + 1],
              "nu" = nu[k + 1],
              "n_iter" = k + 1,
              "phi0_iterate" = phi0,
              "phi1_iterate" = phi1,
              "sigma2_iterate" = sigma2, 
              "nu_iterate" = nu))
}



#' @title Imputate Missing Values in  Incomplete Student's t AR(1) Time Series 
#'
#' @description Estimate the parameters of the Student's t AR(1) model from a time series with missing values and impute the missing values based on the estimates
#'
#' @param y a xts object indicating time series with missing values. The first and last one should not be NA.
#' @param n_sample a positive integer indicating the number of imputations.
#' @param n_burn a positive integer controling the length of the burn-in period of the Gibb sampling. The first (n_burn * n_thin) samples generated will be ignored.
#' @param n_thin a positive integer indicating the sampling period of Gibbs sampling. After then burn-in perid, every n_thin-th sample is used. This is aimed to reduce the dependence of the samples.
#' @param parameters a vector consisting the paramters of the Student's t AR(1) if known. 
#' @return 
#' \item{\code{y_imputed}  }{a xts object, each column is a imputed complete time series}
#' @author Junyan Liu and Daniel P. Palomar
#' @examples
#' library(imputeFin)
#' library(MASS)
#' library(imputeTS)
#' 
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
#' y_miss <- y_orig
#' y_miss[index_miss] <- NA
#' 
#' # impute the missing values and generate n_sample complete time series
#' y_imputed <- imputeAR1( y_miss, n_sample = 1) # if the parameters are unknown
#' y_imputed <- imputeAR1(y_miss, n_sample = 1, parameters = c(phi0, phi1, sigma2, nu)) # if the parameters are unknown
#' 
#' @export
imputeAR1 <- function(y, n_sample = 1, n_burn = 100, n_thin = 50, parameters = NA) {
  # if the parameters are known, then estimate the parameters.
  if(any(is.na(parameters))){
    estimation_result <- estimateAR1(y)
    phi0 <- estimation_result$phi0
    phi1 <- estimation_result$phi1
    sigma2 <- estimation_result$sigma2
    nu <- estimation_result$nu
  }else{
    phi0 <- parameters[1]
    phi1 <- parameters[2]
    sigma2 <- parameters[3]
    nu <- parameters[4]
  }

  # impute the missing y and generate complete data sets
    y_imputed <-  imputeMissingValues(y, n_sample, n_burn, n_thin, y_sample_init = na.kalman(y),
                                      phi0, phi1, sigma2, nu)
    return(y_imputed)
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
#
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
      
      y_tmp[ first_index_in_block[d] : last_index_in_block[d]] <- MASS::mvrnorm  ( n = 1, mu = mu_d, Sigma = sigma_d )
    }
  }
  return(list('y' = y_tmp, 'tau' = tau_tmp))
}

# impute y_miss and generate complete time series y's  via Gibbs sampling.
#   y: a xts object indicating time series with missing values. The first and last one should not be NA.
#   n_sample: the number of complete time series that need to be generated.
#   n_burn:  a positive integer controling the length of the burn-in period of the Gibb sampling. The first (n_burn * n_thin) samples generated will be ignored.
#   n_thin:  a positive integer indicating the sampling period of Gibbs sampling. After then burn-in perid, every n_thin-th sample is used. This is aimed to reduce the dependence of the samples.
#   y_sample_init: a vector of real numbers indicating the initialization for y.
#   phi0: a real number, a parameter of the student's t AR(1) model.
#   phi1: a real number, a parameter of the student's t AR(1) model.
#   sigma2: a positive number, a parameter of the student's t AR(1) model.
#   nu: a positive number indicating the degree of freedom, a parameter of the student's t AR(1) model.
#
imputeMissingValues <- function(y, n_sample, n_burn, n_thin, y_sample_init,
                                phi0, phi1, sigma2, nu) {
  # find the missing blocks
  n <- length(y)               # length of the time series
  index_obs <- which( !is.na(y))# indexes of observed values
  delta_index_obs <- diff(index_obs)
  index_tmp <- which( delta_index_obs>1)
  n_block <- length(index_tmp)                             # number of missing blocks
  n_in_block <- delta_index_obs[index_tmp] - 1             # number of missing values in each block
  first_index_in_block <- index_obs[index_tmp] + 1         # index of the first missing value in each block
  last_index_in_block <- index_obs[index_tmp] + n_in_block # index of the last missing value in each block
  previous_obs_before_block <- as.numeric( y[first_index_in_block - 1] ) # previous observed value before each block
  next_obs_after_block <- as.numeric( y[last_index_in_block + 1] )       # next observed value after each block
  
  if (any(is.na(y_sample_init)))
    y_sample_init <- imputeTS::na.kalman(y)
  y_tmp <- as.numeric( y_sample_init)
  y_imputed <- matrix(nrow = n, ncol = n_sample)
  
  # burn-in period
  for (i in 1:n_burn) {
    sample <- samplingLatentVariables(y_tmp, n_thin = 1, n_block, n_in_block,
                                     first_index_in_block, last_index_in_block, previous_obs_before_block, next_obs_after_block,
                                     phi0, phi1, sigma2, nu) 
    y_tmp <- sample$y
  }
  # sample every n_thin-th sample
  for (j in 1:n_sample) {
    sample <- samplingLatentVariables(y_tmp, n_thin, n_block, n_in_block,
                                      first_index_in_block, last_index_in_block, previous_obs_before_block, next_obs_after_block,
                                      phi0, phi1, sigma2, nu) 
    y_imputed[, j] <- sample$y
  }
  y_imputed <- xts(y_imputed, index(y))
  return(y_imputed)
}
