#' @export
plotImputed <- function(y_imputed, column = 1, type = c("ggplot2", "simple")){
  if (is.matrix(y_imputed)) {
    y_imputed_i <- y_imputed[, column]
    index_miss <- attributes(y_imputed)$index_miss[[column]]
  } else {
    y_imputed_i <- y_imputed
    index_miss <- attributes(y_imputed)$index_miss
  }
  index <- index(y_imputed_i)
  
  switch(match.arg(type),
         "simple" = {
           plot(index, y_imputed_i,  type = "l",  col = "black", xlab = "time", ylab = "values", main = "imputed time series")
           grid()
           # points(index[index_obs], y_imputed_i[index_obs],  pch = 19, col = "#0066CC")
           points(index[index_miss], y_imputed_i[index_miss], pch = 19, col = "red")
           # legend(x = "topleft", legend= c("real values", "imputed values"), col = c("#0066CC", "#FF3300"), pch = 19)
         },
         "ggplot2" = {
           data_frm  <- data.frame(index, value = as.vector(y_imputed_i))
           ggplot2::ggplot() +
             ggplot2::geom_line(data = data_frm, ggplot2::aes(x = index, y = value), col = "black") +
             # ggplot2::geom_point(data = data_frm[index_obs,], ggplot2::aes(x = index, y = value), col = "#0066CC") +
             ggplot2::geom_point(data = data_frm[index_miss,], ggplot2::aes(x = index, y = value), col = "#FF3300") +
             ggplot2::ggtitle("imputed time series") + ggplot2::xlab("time") + ggplot2::ylab("value") +
             ggplot2::theme_bw()
         },
         stop("Table type unknown"))
}



