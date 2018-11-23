library(imputeFin)
library(MASS)
library(imputeTS)
library(xts)
library(quantmod)

# download data from YahooFinance
hsi_price <- Ad(getSymbols("^HSI", from = "2018-01-01", to = "2018-6-30", auto.assign = FALSE))
y <- as.vector(hsi_price)

# creat NA's
miss_percentage <- 0.2
n <- length(y)
n_miss <- floor(miss_percentage * n)
index_miss <- sample(2:(n - 1), n_miss, FALSE)
index_miss <- sort(index_miss)
y[index_miss] <- NA
index_obs <- setdiff(1:n, index_miss)

# estimate the parameters from this incomplete time series
estimation_result <- estimateAR1(y,  n_chain = 10, n_thin = 1, n_iteration = 100, K = 30,
                     estimates_init = NA,  y_sample_init = NA) 
# plot the estimates versus iteration
layout(matrix(c(1, 2, 3, 4)), heights = c(1,1,1,1))
par(mar=c(4, 5, 0.1, 0.5))
k= 101
plot(1:k,estimation_result$phi0_iterate, pch=19, cex = 0.8,  cex.lab = 1.5, cex.axis = 1.2, xaxt="n",  xlab =" ", type='o', col = "green",   ylab = expression(phi[0]))
plot(1:k,estimation_result$phi1_iterate,  pch=19, cex.lab = 1.5, cex.axis = 1.2, xaxt="n",  xlab =" ",type='o', col = "green", lty = 1,  ylab = expression(phi[1]))
plot(1:k,estimation_result$sigma2_iterate, pch=19, cex.lab = 1.5, cex.axis = 1.2, xaxt="n",  xlab =" ",type='o', col = "green",    ylab = expression(sigma^2) )
plot(1:k,estimation_result$nu_iterate, pch=19, cex.lab = 1.5, cex.axis = 1.2,  type='o',col = "green",  xlab ="Iteration k", ylab = expression(nu) )

# impute the missing values and generate n_sample complete time series
y_imputed <- imputeAR1( y, n_sample = 1, n_burn = 200, n_thin = 50) 
hsi_price_imputed = xts(x = y_imputed, order.by = index(hsi_price))

# plot the imputed data set
layout(1, heights = 1)
par(mar=c(4, 5, 4, 5))
plot(hsi_price_imputed, type = 'o', pch = 19, main = 'Heng Seng Index')
points(hsi_price_imputed[index_miss],  col = "red", cex = 0.95,  pch=19)
addLegend(legend.loc = "topright", legend.names = c("observed values", "imputed values"), pch = c(19, 19),  col = c("black","red"))

