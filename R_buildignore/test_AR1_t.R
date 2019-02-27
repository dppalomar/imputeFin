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
data <- data[(n_drop + 1):n_total] # drop the first n_drop to reduce the influence of initial point
y_orig <- xts(data, seq(as.Date("2016-01-01"), length = n, by = "days"))
 
# creat missing values
index_miss <- sort(sample(2:(n - 1), n_miss, FALSE))
y <- y_orig
y[index_miss] <- NA

# test the estimation function
estimation_result <- estimateAR1t(y, random_walk = FALSE, zero_mean = FALSE, n_chain = 10, n_thin = 1, 
                         n_iter = 100, K = 30, output_iterates = TRUE)
layout(matrix(c(1, 2, 3, 4)), heights = c(1,1,1,1))
par(mar=c(4, 5, 0.1, 0.5))
n_iter = length(estimation_result$phi0_iterate)
plot(1:n_iter, estimation_result$phi0_iterate, pch = 19, cex = 0.8,  cex.lab = 1.5, cex.axis = 1.2, xaxt="n",  xlab =" ", type='o', col = "green",   ylab = expression(phi[0]))
plot(1:n_iter, estimation_result$phi1_iterate,  pch = 19, cex.lab = 1.5, cex.axis = 1.2, xaxt="n",  xlab =" ",type='o', col = "green", lty = 1,  ylab = expression(phi[1]))
plot(1:n_iter, estimation_result$sigma2_iterate, pch = 19, cex.lab = 1.5, cex.axis = 1.2, xaxt="n",  xlab =" ",type='o', col = "green",    ylab = expression(sigma^2))
plot(1:n_iter, estimation_result$nu_iterate, pch = 19, cex.lab = 1.5, cex.axis = 1.2,  type='o',col = "green",  xlab ="Iteration k", ylab = expression(nu))

# test the imputation function and compare with the Gaussian imputed result
imputation_result <- imputeAR1t(y, n_sample = 1, random_walk = FALSE, zero_mean = FALSE, param = estimation_result, n_burn = 100, n_thin = 50)
layout(1, heights = 1)
par(mar=c(4, 5, 4, 5))
plot(1:n, imputation_result, type = 'o', pch = 19,  ylab = 'y')
points(index_miss, imputation_result[index_miss, 1], col = 'red', pch = 19)
legend( x="bottomright", col = c("black","red"), pch = c(19, 19), legend = c("observed values", "imputed values") )
