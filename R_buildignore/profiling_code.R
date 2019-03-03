library(imputeFin)
library(xts)

# generate the complete time series
phi0 <- 0
phi1 <- 1
sigma2 <- 0.01 
n <- 500000
n_miss <- 0.3*n
n_drop <- 10
n_total <- n + n_drop
data <- vector(length = n_total)
epsilon <- vector(length = n_total)  # innovations
data[1] <- 0
for (i in 2:n_total) {
  epsilon[i] <- rnorm(1, 0, sqrt(sigma2)) 
  data[i] <- phi0 + phi1 * data[i - 1] + epsilon[i]
}
data <- data[(n_drop + 1):n_total]  # drop the first n_drop to reduce the influence of initial point
y_orig <- xts(data,  seq(as.Date("2016-01-01"), length = n, by = "days"))

# creat missing values
index_miss <- sort(sample(2:(n - 1), n_miss))
y <- y_orig
y[index_miss] <- NA

# test the estimation function
estimation_result <- estimateAR1Gaussian(y, random_walk = FALSE, zero_mean = FALSE, ftol = 1e-10,  
                                         maxiter = 1000, output_iterates = TRUE)  

# profiling
library(profvis)
diag1 <- function(X) {
  m <- min(dim(X))
  X[1 + dim(X)[1L] + 0L:(m - 2L) * (dim(X)[1L] + 1)]  # main diag: x[1 + 0L:(m - 1L) * (dim(x)[1L] + 1)]
}

profvis({
  estimateAR1Gaussian <- function(y, random_walk = FALSE, zero_mean = TRUE, ftol = 1e-10,  
                                  maxiter = 1000, output_iterates = FALSE) {
    #y <- as.vector(y)
    y <- as.numeric(y)

    # find the missing blocks
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
    previous_obs_before_block <- as.numeric(y[first_index_in_block - 1] )  # previous observed value before each block
    next_obs_after_block <- y[last_index_in_block + 1]  # next observed value after each block
    
    
    # objective function, the observed data log-likelihood
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
    estimation_heuristic <- estimateAR1GaussianHeuristic(y, index_miss, random_walk, zero_mean)
    phi1[1] <- estimation_heuristic$phi1
    phi0[1] <- estimation_heuristic$phi0
    sigma2[1] <- estimation_heuristic$sigma2
    f[1] <- obj(phi0[1], phi1[1], sigma2[1])
    
    for (k in 1:maxiter) {
      # E-step
      # computation of mean and covariance of y conditional on all the observed data
      cond <- condMeanCov(y_obs, index_obs, n, n_block, n_in_block, 
                          first_index_in_block, last_index_in_block, previous_obs_before_block, next_obs_after_block, 
                          phi0[k], phi1[k], sigma2[k], full_cov = FALSE)
      # computation of sufficient statistics
      s_y2 <- sum(cond$mean_y[-1])
      s_y1 <- sum(cond$mean_y[-n])
      s_y2y2 <- sum(cond$cov_y_diag[-1] + cond$mean_y[-1]^2)  #slow version: s_y2y2 <- sum(diag(cond_cov_y[2:n, 2:n]) + cond_mean_y[2:n]^2)
      s_y1y1 <- sum(cond$cov_y_diag[-n]  + cond$mean_y[-n]^2)  #slow version: s_y1y1 <- sum(diag(cond_cov_y[1:(n - 1), 1:(n - 1)])  + cond_mean_y[1:(n - 1)]^2)
      s_y2y1 <- sum(cond$cov_y_diag1 + cond$mean_y[-1] * cond$mean_y[-n])  #slow version: s_y2y1 <- sum(diag(cond_cov_y[1:(n - 1), 2:n]) + cond_mean_y[2:n] * cond_mean_y[1:(n - 1)])
      
      # M-step (update the estimates
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
      sigma2[k + 1] <- (( s_y2y2 + (n - 1) * phi0[k + 1]^2 + phi1[k + 1]^2 * s_y1y1 
                          - 2 * phi0[k + 1] * s_y2 - 2 * phi1[k + 1] * s_y2y1 + 2 * phi0[k + 1] * phi1[k + 1] * s_y1 ) / (n - 1))
      
      # computation of the objective function
      f[k + 1] <- obj(phi0[k + 1], phi1[k + 1], sigma2[k + 1])
      
      # stop when the iterates do not change much
      if (abs(f[k + 1] - f[k]) <= ftol * (abs(f[k + 1]) + abs(f[k]))/2) 
        break
    }
    
    if (output_iterates)
      return(list("phi0" = phi0[k + 1],
                  "phi1" = phi1[k + 1],
                  "sigma2" = sigma2[k + 1],
                  "phi0_iterate" = phi0,
                  "phi1_iterate" = phi1,
                  "sigma2_iterate" = sigma2,
                  "f_iterate" = f))
    else
      return(list("phi0" = phi0[k + 1],
                  "phi1" = phi1[k + 1],
                  "sigma2" = sigma2[k + 1]))
  }  
  estimation_result <- estimateAR1Gaussian(y, random_walk = FALSE, zero_mean = FALSE, ftol = 1e-10,  
                                           maxiter = 1000, output_iterates = TRUE)  
})
