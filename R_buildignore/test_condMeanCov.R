library(imputeFin)
library(xts)

# generate the complete time series
phi0 <- 0
phi1 <- 1
sigma2 <- 0.01 
n <- 3000
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

# find the missing blocks
y <- as.numeric(y)
index_obs <- setdiff(1:n, index_miss)  # indexes of missing values
y_obs <- y[index_obs]  # observed values
delta_index_obs <- diff(index_obs)
index_delta_index_obs <- which(delta_index_obs > 1)
n_block <- length(index_delta_index_obs)  # number of missing blocks
n_in_block <- delta_index_obs[index_delta_index_obs] - 1  # number of missing values in each block
first_index_in_block <- index_obs[index_delta_index_obs] + 1  # index of the first missing value in each block
last_index_in_block <- index_obs[index_delta_index_obs] + n_in_block  # index of the last missing value in each block
previous_obs_before_block <- y[first_index_in_block - 1]  # previous observed value before each block
next_obs_after_block <- y[last_index_in_block + 1]  # next observed value after each block

# test the function condMeanCov
cond <- condMeanCov(y_obs, index_obs, n, n_block, n_in_block, 
                    first_index_in_block, last_index_in_block, previous_obs_before_block, next_obs_after_block, 
                    phi0, phi1, sigma2, full_cov = TRUE)

# check the diag & diag1
sum(abs(cond$cov_y_diag - diag(cond$cov_y)))
sum(abs(cond$cov_y_diag1 - diag1(cond$cov_y))) 
