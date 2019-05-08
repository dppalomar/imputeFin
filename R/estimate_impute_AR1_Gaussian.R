#' @title Fit Gaussian AR(1) Model to Time Series with Missing Values
#'
#' @description Estimate the parameters of a Gaussian AR(1) model from a time series with missing values
#'
#' @param y a xts object indicating time series with missing values. The first and last one should not be NA.
#' @param random_walk logical. If TRUE, y is a random walk time series, and phi1 = 1. If FALSE, y is a general AR(1) time series, and phi1 is unknown. The default value is FALSE. 
#' @param zero_mean logical. If TRUE, y is a zero-mean time series, and phi0 = 1. If FALSE, y is a general AR(1) time series, and phi0 is unknown. The default value is FALSE.
#' @param ftol a positive number controlling the stopping criterion (default \code{1e-8}).
#' @param maxiter a positive integer indicating the maximum number of iterations allowed (default \code{1000}).
#' @param iterates logical. If TRUE, then the iterates are outputted. If FALSE, they are ignored. The default value is FALSE.
#' @return A list containing the following elements:
#' \item{\code{phi0}}{real number, the estimate for phi0}
#' \item{\code{phi1}}{real number, the estimate for phi1}
#' \item{\code{sigma2}}{positive number, the estimate for sigma^2}
#' \item{\code{phi0_iterate}}{a numerical vector, the estimates for phi0 in each iteration}
#' \item{\code{phi1_iterate}}{a numerical vector, the estimates for phi1 in each iteration}
#' \item{\code{sigma2_iterate}}{a vector of positive numbers, the estimates for sigma^2 in each iteration}
#' \item{\code{f_iterate}}{a numerical vector, the objective values in each iteration}
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
#'   epsilon[i-1] <- rnorm(1, nu) * sqrt(sigma2)
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
#' estimation_result <- estimateAR1Gaussian(y)
#' 
#' @references
#' R. J. Little and D. B. Rubin, Statistical Analysis with Missing Data, 2nd ed. Hoboken, N.J.: John Wiley & Sons, 2002.
#' 
#' @export
estimateAR1Gaussian <- function(y, random_walk = FALSE, zero_mean = FALSE,
                                iterates = FALSE, condMeanCov = FALSE,
                                tol = 1e-10,  maxiter = 1000) {

  if (NCOL(y) > 1) {
    estimation_list <- apply(y, MARGIN = 2, FUN = estimateAR1Gaussian, random_walk, zero_mean, iterates, condMeanCov, tol, maxiter)
    phi0 <- unlist(lapply(estimation_list, function(x){x$phi0}))
    phi1 <- unlist(lapply(estimation_list, function(x){x$phi1}))
    sigma2 <- unlist(lapply(estimation_list, function(x){x$sigma2}))
    return(c(estimation_list, list("phi0" = phi0,
                                   "phi1" = phi1,
                                   "sigma2" = sigma2)))
  }
  y <- as.numeric(y)
  
  # trivial case with no NAs
  if (!anyNA(y)) {
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
    return(list("phi0" = phi0,
                "phi1" = phi1,
                "sigma2" = sigma2))
  }
  
  # find the missing blocks
  list2env(findMissingBlock(y), envir = environment())
  
  # objective function, the observed data log-likelihood
  if (iterates)
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
  if (iterates) f[1] <- obj(phi0[1], phi1[1], sigma2[1])

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
    sigma2[k + 1] <- (( s_y2y2 + (n - 1) * phi0[k + 1]^2 + phi1[k + 1]^2 * s_y1y1 
                      - 2 * phi0[k + 1] * s_y2 - 2 * phi1[k + 1] * s_y2y1 + 2 * phi0[k + 1] * phi1[k + 1] * s_y1 ) / (n - 1))
    
    # computation of the objective function
    if (iterates) f[k + 1] <- obj(phi0[k + 1], phi1[k + 1], sigma2[k + 1])

    # termination criterion    
    if (abs(phi0[k + 1] - phi0[k]) <= tol * (abs(phi0[k + 1]) + abs(phi0[k]))/2
        && abs(phi1[k + 1] - phi1[k]) <= tol * (abs(phi1[k + 1]) + abs(phi1[k]))/2
        && abs(sigma2[k + 1] - sigma2[k]) <= tol * (abs(sigma2[k + 1]) + abs(sigma2[k]))/2) 
    break
    
  }
  
  results <- list("phi0" = phi0[k + 1],
                  "phi1" = phi1[k + 1],
                  "sigma2" = sigma2[k + 1])
  
  if (iterates) 
    results <- c(results, list("phi0_iterate" = phi0,
                               "phi1_iterate" = phi1,
                               "sigma2_iterate" = sigma2,
                               "f_iterate" = f))

  if (condMeanCov) {
    cond <- condMeanCov(y_obs, index_obs, n, n_block, n_in_block, 
                        first_index_in_block, last_index_in_block, previous_obs_before_block, next_obs_after_block, 
                        phi0[k+1], phi1[k+1], sigma2[k+1], full_cov = TRUE)
    results <- c(results, list("cond_mean_y" = cond$mean_y,
                               "cond_cov_y" = cond$cov_y))
  }
  return(results)
}


#
# Extracts the diagonal on top of the main diagonal
#
diag1 <- function(X) {
  m <- min(dim(X))
  X[1 + dim(X)[1L] + 0L:(m - 2L) * (dim(X)[1L] + 1)]  # main diag: x[1 + 0L:(m - 1L) * (dim(x)[1L] + 1)]
}


#' @title Imputate Missing Values in  Incomplete Gaussian AR(1) Time Series 
#'
#' @description Estimate the parameters of the Gaussian AR(1) model from a time series with missing values and impute the missing values based on the estimates
#'
#' @param y a xts object indicating time series with missing values. The first and last one should not be NA.
#' @param n_samples a positive integer indicating the number of imputations (default \code{1}).
#' @param param a list consisting of the paramters of the Student's t AR(1) time series y if known. The default value is FALSE.
#' @param random_walk logical. If TRUE, y is a random walk time series, and phi1 = 1. If FALSE, y is a general AR(1) time series, and phi1 is unknown. The default value is FALSE.
#' @param zero_mean logical. If TRUE, y is a zero-mean time series, and phi0 = 1. If FALSE, y is a general AR(1) time series, and phi0 is unknown.
#' @return 
#' \item{\code{y_imputed}  }{a numerical matrix, each column is a imputed complete time series}
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
#'   epsilon[i - 1] <- rnorm(1, nu) * sqrt(sigma2)
#'   data[i] <- phi0 + phi1 * data[i - 1] + epsilon[i - 1]
#' }
#' data <- data[(n_drop + 1):n_total] # drop the first n_drop to reduce the influence of initial point
#' dates <- seq(as.Date("2016-01-01"), length = n, by = "days") 
#' y_orig <- xts(data,  dates)
#' 
#' # creat missing values
#' index_miss <- sample(2:(n - 1), n_miss, FALSE)
#' index_miss <- sort(index_miss)
#' y <- y_orig
#' y[index_miss] <- NA
#' 
#' # impute the missing values and generate n_samples complete time series
#' y_imputed <- imputeAR1Gaussian( y_miss, n_samples = 3) # if the parameters are unknown
#' param = list("phi0" = phi0,
#'              "phi1" = phi1,
#'              "sigma2" = sigma2,
#'              "nu" = nu)
#' y_imputed <- imputeAR1Gaussian(y_miss, n_samples = 3, param) # if the parameters are unknown
#' @export
imputeAR1Gaussian <- function(y, n_samples = 1, random_walk = FALSE, zero_mean = FALSE,
                              estimates = FALSE) {

  if (NCOL(y) > 1) {
    results_list <- lapply(c(1:NCOL(y)), FUN = function(i){imputeAR1Gaussian(y[, i], n_samples, random_walk, zero_mean, estimates)})
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
      }
    }
    return(results)
  }
  
  y_attrib <- attributes(y)
  y <- as.numeric(y)
  
  # trivial case with no NAs
  if (!anyNA(y)){
    y_imputed <- matrix(rep(y, times = n_samples), ncol = n_samples)
    if (estimates) estimation_result <- estimateAR1Gaussian(y, random_walk, zero_mean)
    index_miss <- NULL
  } else {
    estimation_result <- estimateAR1Gaussian(y, random_walk, zero_mean, condMeanCov = TRUE)
    cond_mean_y <- estimation_result$cond_mean_y
    cond_cov_y <- estimation_result$cond_cov_y
    n <- length(y)  # length of the time series
    index_obs <- which(!is.na(y))  # indexes of observed values
    index_miss <- setdiff(1:n, index_obs)  # indexes of missing values
    y_obs <- y[index_obs]  # observed values
    # impute the missing values by drawing samples from its conditional distribution
    y_imputed <- matrix(nrow = n, ncol = n_samples)
    y_imputed[index_miss, ] <- t(MASS::mvrnorm(n = n_samples, cond_mean_y[index_miss], cond_cov_y[index_miss, index_miss]))
    y_imputed[index_obs, ] <- rep(y_obs, times = n_samples)
  }

  if (n_samples == 1) {
    attributes(y_imputed) <- y_attrib
    attr(y_imputed, "index_miss") <- index_miss
    if (!estimates) {
      results <- y_imputed
    } else
      results <- list("y_imputed" = y_imputed)
  } else {
    y_imputed <-lapply(split(y_imputed, col(y_imputed)), FUN = function(x) { attributes(x) <- y_attrib
                                                                             attr(x, "index_miss") <- index_miss
                                                                             return(x) })
    results <- c("y_imputed" = y_imputed)
  }
  
  if (estimates)  results <- c(results, list("phi0" = estimation_result$phi0,
                                             "phi1" = estimation_result$phi1,
                                             "sigma2" = estimation_result$sigma2))
  return(results)
}


# Compute the conditional mean and covariance matrix of y given the observed data and current estimates

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
  
  results <- list("mean_y" = cond_mean_y,
                  "cov_y_diag" = cond_cov_y_diag,
                  "cov_y_diag1" = cond_cov_y_diag1)
  if (full_cov) 
    results <- c(results, list("cov_y" = cond_cov_y))
  return(results)
}


# A heuristic method to compute the parameters of Gaussian AR(1) model from incomplete time series
 
estimateAR1GaussianHeuristic <- function(y, index_miss, random_walk = FALSE, zero_mean = TRUE) {
  index_miss_p <- c(0, index_miss, length(y) + 1)
  delta_index_miss_p <- diff(index_miss_p)
  index_delta_index_miss_p <- which(delta_index_miss_p > 2)
  n_obs_block <- length(index_delta_index_miss_p)  # number of observation blocks with more than 1 sample
  n_in_obs_block <- delta_index_miss_p[index_delta_index_miss_p] - 1  # number of observed samples in each qualified observation block
  
  m <- 0
  y_obs2 <- y_obs1 <- c()
  for (i in 1:n_obs_block) {
    y_obs1[(m + 1):(m + n_in_obs_block[i] - 1)] <- as.numeric( y[(index_miss_p[index_delta_index_miss_p[i]] + 1):(index_miss_p[index_delta_index_miss_p[i] + 1] - 2)])
    y_obs2[(m + 1):(m + n_in_obs_block[i] - 1)] <- as.numeric( y[(index_miss_p[index_delta_index_miss_p[i]] + 2):(index_miss_p[index_delta_index_miss_p[i] + 1] - 1)])
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
  
  return(list("phi0" = phi0,
              "phi1" = phi1,
              "sigma2" = sigma2))
}


# find the missing blocks

findMissingBlock <- function(y){

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
  
  return(list("n" = n,
              "index_obs" = index_obs,
              "index_miss" = index_miss,
              "n_obs" =  n_obs,
              "y_obs" = y_obs,
              "delta_index_obs" =  delta_index_obs,
              "n_block" = n_block,
              "n_in_block" = n_in_block,
              "first_index_in_block" = first_index_in_block,
              "last_index_in_block" = last_index_in_block,
              "previous_obs_before_block" = previous_obs_before_block,
              "next_obs_after_block" = next_obs_after_block))
}

