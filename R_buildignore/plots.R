

# plot using ggplot2::qplot
# library(ggplot2)
# library(gridExtra)
# n_iter <- length(estimation_result$phi0_iterate)
# p1 <- qplot(1:n_iter, estimation_result$phi0_iterate, geom = c("line", "point"), xlab = "", ylab = expression(phi[0])) + 
#  scale_x_continuous(breaks = 1:n_iter)
# p2 <- qplot(1:n_iter, estimation_result$phi1_iterate, geom = c("line", "point"), xlab = "", ylab = expression(phi[1])) +
#  scale_x_continuous(breaks = 1:n_iter)
# p3 <- qplot(1:n_iter, estimation_result$sigma2_iterate, geom = c("line", "point"), xlab = "", ylab = expression(sigma^2)) + 
#  scale_x_continuous(breaks = 1:n_iter)
# p4 <- qplot(1:n_iter, estimation_result$f_iterate, geom = c("line", "point"), col = I("darkblue"), xlab = "iteration", ylab = "obj") + 
#  scale_x_continuous(breaks = 1:n_iter)
# grid.arrange(p1, p2, p3, p4, ncol = 1)


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



# index_miss_p <- (min(index_miss)-1):(max(index_miss) + 1)
# par(mfrow=c(2,1))
# plot 1
# { plot(y[, 1], main = "Original")
#  lines(y_orig[index_miss_p, 1], col = "blue", lwd = 2) }
# plot 2
# { plot(y[ ,1], main = "Imputed")
#  lines(y_imputed[index_miss_p, 1], col = "blue", lwd = 2) }


# 
# # plot using ggplot2 (http://www.sthda.com/english/articles/32-r-graphics-essentials/128-plot-time-series-data-using-ggplot/)
# # autoplot(y) + labs(title = "Original", x = "", y = "") + theme_bw()
# # qplot(Index, Value, data = fortify(y, melt=TRUE), geom = "line",
# #       main = "Original", xlab = "", ylab = "") + theme_bw()
# p1 <- ggplot() +
#   geom_line(data = fortify(y[, 1], melt = TRUE), aes(Index, Value), color = I("black")) +
#   geom_line(data = fortify(y_orig[index_miss_p, 1], melt = TRUE), aes(Index, Value), color = I("red")) +
#   scale_x_date(date_labels = "%b %Y", date_breaks = "3 month") +
#   labs(title = "Original", x = "", y = "") + 
#   theme_bw()
# p2 <- ggplot() +
#   geom_line(data = fortify(y[, 1], melt = TRUE), aes(Index, Value), color = I("black")) +
#   geom_line(data = fortify(y_imputed[index_miss_p, 1], melt = TRUE), aes(Index, Value), color = I("red")) +
#   scale_x_date(date_labels = "%b %Y", date_breaks = "3 month") +
#   labs(title = "Imputed", x = "", y = "") + 
#   theme_bw()
# grid.arrange(p1, p2, ncol = 1)

