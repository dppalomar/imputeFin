# generate the data 
library(xts)
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

m <- 3
data_mtr <- matrix(rep(data, m), nrow = n, ncol = m)

y_orig <- xts(data_mtr, seq(as.Date("2016-01-01"), length = n, by = "days"))
#y_orig <- data_mtr
colnames(y_orig) <- c("a", "b", "c")

#numerical vector
#y_orig <- data_mtr

# creat missing values
index_miss <- sort(sample(2:(n - 1), n_miss))
# index_miss <- round(n/2) + 1:n_miss
y_missing <- y_orig
y_missing[index_miss, 1] <- NA
y_missing[c(5,10,12), 2] <- NA

#
# this is to save the data to the package
#
#save(y_missing, file = "data-raw/y_missing.RData")
usethis::use_data_raw()
usethis::use_data(y_missing, overwrite = TRUE)



AR1_Gaussian$y_missing
            $phi0
AR1_Studentt

