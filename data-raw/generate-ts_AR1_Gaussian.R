# generate complete data
phi0 <- 1
phi1 <- 0.5
sigma2 <- 0.01 
n <- 300
n_miss <- 0.1*n
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

# repeat the time series m times
m <- 3
y_orig <- matrix(rep(data, m), nrow = n, ncol = m)
colnames(y_orig) <- c("a", "b", "c")

# create missing values
index_miss1 <- round(n/2) + 1:n_miss
index_miss2 <- sort(sample(2:(n - 1), n_miss))
index_miss3 <- union(index_miss1, index_miss2)
#index_miss3 <- c(1,sort(sample(2:(n - 1), n_miss - 2)), n)
y_missing_numeric <- y_orig
y_missing_numeric[index_miss1, 1] <- NA
y_missing_numeric[index_miss2, 2] <- NA
y_missing_numeric[index_miss3, 3] <- NA
y_missing <- zoo::zoo(y_missing_numeric, seq(as.Date("2016-01-01"), length = n, by = "days"))


# create variable
ts_AR1_Gaussian <- list("y_missing" = y_missing,
                        "phi0"      = phi0,
                        "phi1"      = phi1,
                        "sigma2"    = sigma2)

# save the data to the package
save(ts_AR1_Gaussian, file = "data/ts_AR1_Gaussian.RData", version = 2)
