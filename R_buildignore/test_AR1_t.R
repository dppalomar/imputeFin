library(imputeFin)
library(xts)

# generate the complete time series
phi0 <- 0
phi1 <- 1 
sigma2 <- 0.01 
nu <- 2
n <- 300
n_miss <- 30 
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
data_mtr <- matrix(rep(data, m), nrow = n, ncol = m)

y_orig <- xts(data_mtr,  seq(as.Date("2016-01-01"), length = n, by = "days"))

# creat missing values
# index_miss <- sort(sample(2:(n - 1), n_miss))
index_miss <- round(n/2) + 1:n_miss
y <- y_orig
y[index_miss,] <- NA

# test the estimation function
estimation_result <- estimateAR1t(y, random_walk = FALSE, zero_mean = FALSE, output_iterates = TRUE, 
                                  n_chain = 10, n_thin = 1, n_iter = 100, K = 30)
library(ggplot2)
library(gridExtra)
n_iter <- length(estimation_result$phi0_iterate)
p1 <- qplot(1:n_iter, estimation_result$phi0_iterate, geom = c("line", "point"), xlab = "", ylab = expression(phi[0]))  
p2 <- qplot(1:n_iter, estimation_result$phi1_iterate, geom = c("line", "point"), xlab = "", ylab = expression(phi[1]))
p3 <- qplot(1:n_iter, estimation_result$sigma2_iterate, geom = c("line", "point"), xlab = "", ylab = expression(sigma^2)) 
p4 <- qplot(1:n_iter, estimation_result$nu_iterate, geom = c("line", "point"), xlab = "", ylab = expression(nu))  
grid.arrange(p1, p2, p3, p4, ncol = 1)

# test the imputation function and compare with the Gaussian imputed result
imputation_result <- imputeAR1t(y, n_sample = 2, param = NULL, random_walk = FALSE, zero_mean = FALSE, n_burn = 1000, n_thin = 50)
y_imputed <- imputation_result$y_imputed
index_miss_p <- (min(index_miss)-1):(max(index_miss) + 1)

p1 <- ggplot() +
  geom_line(data = fortify(y, melt = TRUE), aes(Index, Value), color = I("black")) +
  geom_line(data = fortify(y_orig[index_miss_p], melt = TRUE), aes(Index, Value), color = I("red")) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 month") +
  labs(title = "Original", x = "", y = "") + 
  theme_bw()
p2 <- ggplot() +
  geom_line(data = fortify(y, melt = TRUE), aes(Index, Value), color = I("black")) +
  geom_line(data = fortify(y_imputed[index_miss_p, 1], melt = TRUE), aes(Index, Value), color = I("red")) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 month") +
  labs(title = "Imputed", x = "", y = "") + 
  theme_bw()
grid.arrange(p1, p2, ncol = 1)

