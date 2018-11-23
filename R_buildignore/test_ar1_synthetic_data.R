library(imputeFin)
library(MASS)
library(imputeTS)

# generate a complete Student's t AR(1) time series
phi0 <- 1
phi1 <- 0.5 
sigma2 <- 0.01 
nu <- 1
n <- 200
n_miss <- 25 
n_drop <- 100
n_total <- n + n_drop
y <- vector(length = n_total)
epsilon <- vector(length = n_total - 1)# innovations
y[1] <- 0
for (i in 2:n_total) {
  epsilon[i-1] <- rt(1, nu) * sqrt(sigma2)
  y[i] <- phi0 + phi1 * y[i-1] + epsilon[i-1]
}
y <- y[(n_drop + 1):n_total] # drop the first n_drop to reduce the influence of initial point
  
# creat missing values

index_miss <- sample( 2:(n - 1), n_miss, FALSE)
index_miss <- sort(index_miss)
y[index_miss] <- NA

# estimate the parameters from this incomplete time series
estimation_result <- estimateAR1(y)

# plot the estimates versus iterations
layout(matrix(c(1, 2, 3, 4)), heights = c(1,1,1,1))
par(mar=c(4, 5, 0.1, 0.5))
k= 101
plot(1:k,estimation_result$phi0_iterate, pch=19, cex = 0.8,  cex.lab = 1.5, cex.axis = 1.2, xaxt="n",  xlab =" ", type='o', col = "green",   ylab = expression(phi[0]))
plot(1:k,estimation_result$phi1_iterate,  pch=19, cex.lab = 1.5, cex.axis = 1.2, xaxt="n",  xlab =" ",type='o', col = "green", lty = 1,  ylab = expression(phi[1]))
plot(1:k,estimation_result$sigma2_iterate, pch=19, cex.lab = 1.5, cex.axis = 1.2, xaxt="n",  xlab =" ",type='o', col = "green",    ylab = expression(sigma^2) )
plot(1:k,estimation_result$nu_iterate, pch=19, cex.lab = 1.5, cex.axis = 1.2,  type='o',col = "green",  xlab ="Iteration k", ylab = expression(nu) )

# impute the missing values and generate n_sample complete time series
#y_imputed <- imputeAR1( y, n_sample = 1) # if the parameters are unknown
y_imputed <- imputeAR1(y, n_sample = 1, parameters = c(phi0, phi1, sigma2, nu)) # if the parameters are unknown
layout(1, heights = 1)
par(mar=c(4, 5, 4, 5))
plot(1:n, y_imputed, type = 'o', pch = 19,  ylab = 'y')
points(index_miss, y_imputed[index_miss, 1], col = 'red', pch = 19)
legend( x="topright", col = c("black","red"), pch = c(19, 19), legend = c("observed values", "imputed values") )

