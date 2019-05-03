#' @export
plotImputed <- function(y_imputed, i){
  y_imputed_i <- y_imputed[,i]
  index_miss <- attributes(y_imputed)$index_miss[[i]]
  index_obs <- setdiff(c(1:length(y_imputed_i)), index_miss)
  
# #  using ggplot
  # data_frm  <- data.frame(index = index(y_imputed_i),
  #                       value = as.vector(y_imputed_i))
  # ggplot() +
  #   geom_line(data = data_frm, aes(x = index, y = value), col = "grey") +
  #   geom_point(data = data_frm[index_obs,], aes(x = index, y = value), col = "#0066CC") +
  #   geom_point(data = data_frm[index_miss,], aes(x = index, y = value), col = "#FF3300") +
  #   ggtitle("imputed time series") + xlab("time") + ylab("value") +
  #   theme_bw()
  
  # using just plot
  index <- index(y_imputed_i)
  plot(index, y_imputed_i,  type = "l",  col = "#666666", xlab = "time", ylab = "values", main = "imputed time series")
  grid()
  points(index[index_obs], y_imputed_i[index_obs],  pch = 19, col = "#0066CC")
  points(index[index_miss], y_imputed_i[index_miss], pch = 19, col = "#FF3300")
  # legend(x = "topleft", legend= c("real values", "imputed values"), col = c("#0066CC", "#FF3300"), pch = 19)
}






