library(imputeFin)
library(xts)
# generate the complete time series
phi0 <- 1
phi1 <- 0.5
sigma2 <- 0.01 
n <- 500
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

m <- 1
data_mtr <- matrix(rep(data, m), nrow = n, ncol = m)

y_orig <- xts(data_mtr,  seq(as.Date("2016-01-01"), length = n, by = "days"))

# creat missing values
# index_miss <- sort(sample(2:(n - 1), n_miss))
index_miss <- round(n/2) + 1:n_miss
y <- y_orig
y[index_miss,] <- NA


# test the estimation function
estimation_result <- estimateAR1Gaussian(y, random_walk = FALSE, zero_mean = FALSE,
                                         output_iterates = TRUE, condMeanCov = FALSE,
                                         ftol = 1e-10,  maxiter = 1000)

# layout(matrix(c(1, 2, 3, 4)), heights = c(1,1,1,1))
# par(mar = c(4, 5, 0.1, 0.5))
# n_iter <- length(estimation_result$phi0_iterate)
# plot(1:n_iter,estimation_result$phi0_iterate, pch=19, cex = 0.8,  cex.lab = 1.5, cex.axis = 1.2, xaxt="n",  xlab =" ", type='o', col = "green",   ylab = expression(phi[0]))
# plot(1:n_iter,estimation_result$phi1_iterate,  pch=19, cex.lab = 1.5, cex.axis = 1.2, xaxt="n",  xlab =" ",type='o', col = "green", lty = 1,  ylab = expression(phi[1]))
# plot(1:n_iter,estimation_result$sigma2_iterate, pch=19, cex.lab = 1.5, cex.axis = 1.2, xaxt="n",  xlab =" ",type='o', col = "green",    ylab = expression(sigma^2) )
# plot(1:n_iter,estimation_result$f_iterate, pch=19, cex.lab = 1.5, cex.axis = 1.2,  type='o',col = "green",  xlab ="Iteration k", ylab = "obj" )

# plot using ggplot2::qplot
library(ggplot2)
library(gridExtra)
n_iter <- length(estimation_result$phi0_iterate)
p1 <- qplot(1:n_iter, estimation_result$phi0_iterate, geom = c("line", "point"), xlab = "", ylab = expression(phi[0])) + 
  scale_x_continuous(breaks = 1:n_iter)
p2 <- qplot(1:n_iter, estimation_result$phi1_iterate, geom = c("line", "point"), xlab = "", ylab = expression(phi[1])) +
  scale_x_continuous(breaks = 1:n_iter)
p3 <- qplot(1:n_iter, estimation_result$sigma2_iterate, geom = c("line", "point"), xlab = "", ylab = expression(sigma^2)) + 
  scale_x_continuous(breaks = 1:n_iter)
p4 <- qplot(1:n_iter, estimation_result$f_iterate, geom = c("line", "point"), col = I("darkblue"), xlab = "iteration", ylab = "obj") + 
  scale_x_continuous(breaks = 1:n_iter)
grid.arrange(p1, p2, p3, p4, ncol = 1)

# # plot using ggplot2::ggplot
# library(ggplot2)
# library(gridExtra)
# library(tidyr)
# library(dplyr)
# library(tibble)
# 
# as_tibble(estimation_result) %>%
#   mutate("iteration" = c(1:n_iter)) %>%
#   select("iteration", "phi0") %>%
#   ggplot(aes(x = iteration, y = phi0)) +
#   geom_line() + geom_point()
# 
# as_tibble(estimation_result) %>%
#   mutate("iteration" = c(1:n_iter)) %>%
#   select("iteration", "phi0", "phi1") %>%
#   gather("phi0", "phi1", key = "iterate_type", value = "value") %>%
#   ggplot(aes(x = iteration, y = value)) +
#   geom_line() + geom_point() +
#   facet_grid(iterate_type ~ .)


# test the imputation function
param <- list("phi0" = phi0,
              "phi1" = phi1,
              "sigma2" = sigma2)
imputation_result <- imputeAR1Gaussian(y, n_sample = 2, param, random_walk = FALSE, zero_mean = TRUE)
y_imputed <- imputation_result$y_imputed
  
index_miss_p <- (min(index_miss)-1):(max(index_miss) + 1)
par(mfrow=c(2,1))
#plot 1
{ plot(y, main = "Original")
  lines(y_orig[index_miss_p], col = "blue", lwd = 2) }
#plot 2
{ plot(y, main = "Imputed")
  lines(y_imputed[index_miss_p, 1], col = "blue", lwd = 2) }
#plot 3
#{ plot(y, main = "Imputed")
#  lines(y_imputed[index_miss_p, 2], col = "blue", lwd = 2) }

# par(mfrow=c(1,1))


# plot using ggplot2 (http://www.sthda.com/english/articles/32-r-graphics-essentials/128-plot-time-series-data-using-ggplot/)
# autoplot(y) + labs(title = "Original", x = "", y = "") + theme_bw()
# qplot(Index, Value, data = fortify(y, melt=TRUE), geom = "line",
#       main = "Original", xlab = "", ylab = "") + theme_bw()
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


# comparison with other methods
library("imputeTS")
y_imputed_km = na.kalman(y)
#plot 3
{ plot(y, main = "Imputated by kalman method")
  lines(y_imputed_km[index_miss_p], col = "blue", lwd = 2) }

library("knitr")
rst_orig <- ar(y_orig, order.max = 1, demean = TRUE)
est_orig <- list("phi0" = rst_orig$x.mean * (1 - sum(rst_orig$ar)),
                 "phi1" = rst_orig$ar,
                 "sigma2" = rst_orig$var.pred)
rst_km <- ar(y_imputed_km, order.max = 1, demean = TRUE)
est_km <- list("phi0" = rst_km$x.mean * (1 - sum(rst_km$ar)),
                 "phi1" = rst_km$ar,
                 "sigma2" = rst_km$var.pred)
rst_pr <- ar(y_imputed[, 1], order.max = 1, demean = TRUE)
est_pr <- list("phi0" = rst_pr$x.mean * (1 - sum(rst_pr$ar)),
               "phi1" = rst_pr$ar,
               "sigma2" = rst_pr$var.pred)
rst = cbind(est_orig, est_km, est_pr)
print( kable(rst))
