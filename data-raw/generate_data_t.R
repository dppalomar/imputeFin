phi0 <- 0
phi1 <- 1 
sigma2 <- 0.01 
nu <- 2
n <- 300
n_miss <- 0.1 * n 
n_drop <- 100
n_total <- n + n_drop
data <- vector(length = n_total)
epsilon <- vector(length = n_total - 1)# innovations
data[1] <- 0
for (i in 2:n_total) {
  epsilon[i] <- rt(1, nu) * sqrt(sigma2)
  data[i] <- phi0 + phi1 * data[i - 1] + epsilon[i]
}
data <- data[(n_drop + 1):n_total]  # drop the first n_drop to reduce the influence of initial point

m <- 3
y_orig<- matrix(rep(data, m), nrow = n, ncol = m)
colnames(y_orig) <- c("a", "b", "c")

# creat missing values
y_missing_numeric <- y_orig
index_miss1 <- round(n/2) + 1:n_miss
index_miss2 <- sort(sample(2:(n - 1), n_miss))
index_miss3 <- union(index_miss1, index_miss2)
#index_miss3 <- c(1,sort(sample(2:(n - 1), n_miss - 2)), n)

y_missing_numeric[index_miss1, 1] <- NA
y_missing_numeric[index_miss2, 2] <- NA
# y_missing_numeric[index_miss3, 3] <- NA
# index_miss3 <- c(1,sort(sample(2:(n - 1), n_miss - 2)), n)
y_missing <- zoo(y_missing_numeric, seq(as.Date("2016-01-01"), length = n, by = "days"))
ts_AR1_t = list("y_missing" = y_missing,
                "phi0"      = phi0,
                "phi1"      = phi1,
                "sigma2"    = sigma2,
                "nu"        = nu)
#index_miss1 <- which(is.na(ts_AR1_t$y_missing[, 1]))
#index_miss2 <- which(is.na(ts_AR1_t$y_missing[, 2]))
#index_miss3 <- union(index_miss1, index_miss2)
#ts_AR1_t$y_missing[index_miss3, 3] <- NA



# save the data to the package
save(ts_AR1_t, file = "data/ts_AR1_t.RData", version = 2)
