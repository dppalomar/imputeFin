library(imputeFin)
library(MASS)
library(imputeTS)
library(xts)
library(quantmod)
library(tictoc)

# download data from YahooFinance
#y_orig <- log(Ad(getSymbols("^HSI", from = "2018-01-01", to = "2018-6-30", auto.assign = FALSE)))
y_orig <- log(Ad(getSymbols("^GSPC", from = "2012-01-01", to = "2015-07-08", auto.assign = FALSE)))
T <- nrow(y_orig)

# create NA's
miss_pct <- 0.2
n_miss <- floor(miss_pct*T)
# block of consecutive NAs
idx_miss <- round(T/2) + 1:n_miss
y_missing <- y_orig
y_missing[idx_miss] <- NA
# # random pattern of NAs
# idx_miss <- sample(2:(T-1), n_miss)
# idx_miss <- sort(idx_miss)
# y_missing[idx_miss] <- NA
idx_obs <- setdiff(1:T, idx_miss)

# estimate the parameters from this incomplete time series
tic("parameter estimation")
estimation_result <- estimateAR1(y_missing)
toc()
estimation_result$phi0
estimation_result$phi1
estimation_result$sigma2
estimation_result$nu

# plot convergence of estimates versus iteration
layout(matrix(c(1, 2, 3, 4)), heights = c(1, 1, 1, 1))
par(mar=c(4, 5, 0.1, 0.5))
k <- estimation_result$n_iter
plot(1:k, estimation_result$phi0_iterate, pch=19, cex = 0.8,  cex.lab = 1.5, cex.axis = 1.2, xaxt="n",  xlab =" ", type='o', col = "blue",   ylab = expression(phi[0]))
plot(1:k, estimation_result$phi1_iterate, pch=19, cex.lab = 1.5, cex.axis = 1.2, xaxt="n",  xlab =" ",type='o', col = "blue", lty = 1,  ylab = expression(phi[1]))
plot(1:k, estimation_result$sigma2_iterate, pch=19, cex.lab = 1.5, cex.axis = 1.2, xaxt="n",  xlab =" ",type='o', col = "blue",    ylab = expression(sigma^2) )
plot(1:k, estimation_result$nu_iterate, pch=19, cex.lab = 1.5, cex.axis = 1.2,  type='o',col = "blue",  xlab ="Iteration k", ylab = expression(nu) )

# impute the missing values and generate n_sample complete time series
y_imputed <- imputeAR1(y_missing, n_sample = 4, n_burn = 200, n_thin = 100)  #use: num_imputations


# # plot the imputed data set
# layout(1, heights = 1)
# par(mar=c(4, 5, 4, 5))
# plot(y_imputed, type = 'o', pch = 19, main = 'Heng Seng Index')
# points(y_imputed[idx_miss],  col = "red", cex = 0.95,  pch=19)
# addLegend(legend.loc = "topright", legend.names = c("observed values", "imputed values"), pch = c(19, 19),  col = c("black","red"))




idx_miss_bis <- (min(idx_miss)-1):(max(idx_miss)+1)
#setEPS()
#postscript("~/Downloads/SCM_eigenvalues_histogram.eps", width=12, height=8)
par(mfrow=c(2,2))
#plot 1
{ plot(y_missing, main="Imputation or true?")
  lines(y_orig[idx_miss_bis, 1], col="blue", lwd=2) }
#plot 2
{ plot(y_missing, main="Imputation or true?")
  lines(y_imputed[idx_miss_bis, 2], col="blue", lwd=2) }
#plot 3
{ plot(y_missing, main="Imputation or true?")
  lines(y_imputed[idx_miss_bis, 3], col="blue", lwd=2) }
#plot 4
{ plot(y_missing, main="Imputation or true?")
  lines(y_imputed[idx_miss_bis, 4], col="blue", lwd=2) }
#dev.off()

