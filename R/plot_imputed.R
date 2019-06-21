#' @title Plot Imputed Time Series
#'
#' @description Plot imputed time series.
#'
#' @param y_imputed imputed time seies, which can be a numeric vector, numeric matrix, or zoo object with an attribute recording the locations of the missing values.
#' @param column a positive integer indicating the number of column of y_imputed to be plotted (default \code{1}).
#' @param type character string specifying the method to plot the time series, "ggplot2" or "simple". The default value is "simple".
#' 
#' @author Junyan Liu and Daniel P. Palomar
#' 
#' @examples
#' library(imputeFin)
#' data(ts_AR1_t) 
#' y_missing <- ts_AR1_t$y_missing
#' y_imputed <- imputeAR1t(y_missing)
#' plotImputed(y_imputed)
#' 
#' @import zoo
#' @export
plotImputed <- function(y_imputed, column = 1, type = c("ggplot2", "simple")) {
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
           plot(index, y_imputed_i, type = "l",  col = "black", xlab = "time", ylab = "value", main = "Imputed time series")
           grid()
           points(index[index_miss], y_imputed_i[index_miss], pch = 19, col = "red")
           legend(x = "topleft", legend = "imputed values", col = "red", pch = 19)
         },
         "ggplot2" = {
           if (!requireNamespace("ggplot2", quietly = TRUE)) 
             stop("Please install package \"ggplot2\" or choose another plot type.", call. = FALSE)
           data_frm  <- data.frame(index, value = as.vector(y_imputed_i))
           index <- value <- NULL  # ugly hack to deal with CRAN note
           # ggplot2::ggplot() +
           #   ggplot2::geom_line(data = data_frm, ggplot2::aes(x = index, y = value), col = "black") +
           #   ggplot2::geom_point(data = data_frm[index_miss, ], ggplot2::aes(x = index, y = value), col = "#DD3344") +
           #   ggplot2::labs(title = "Imputed time series", x = "time")  # + ggplot2::theme_bw()
           ggplot2::ggplot() +
             ggplot2::geom_line(data = data_frm, ggplot2::aes(x = index, y = value), col = "black") +
             ggplot2::geom_point(data = data_frm[index_miss, ], ggplot2::aes(x = index, y = value, col = "imputed values")) +
             ggplot2::labs(title = "Imputed time series", x = "time") +
             ggplot2::theme(legend.title = ggplot2::element_blank(), legend.justification = c(0.9, 0), legend.position = c(0.98, 0.02)) +
             ggplot2::scale_colour_manual(name = "", values = c("imputed values" = "#DD3344"))
         },
         stop("Plot type unknown"))
}
