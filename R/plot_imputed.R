#' @title Plot imputed time series.
#'
#' @description Plot single imputed time series (as returned by functions 
#'              \code{\link{impute_AR1_Gaussian}} and \code{\link{impute_AR1_t}}),
#'              indicating the imputed values with red dots.
#'
#' @param y_imputed Imputed time series (can be any object coercible to a numeric vector 
#'                  or a numeric matrix). If it has the attribute \code{"index_miss"} (as
#'                  returned by any of the imputation functions 
#'                  \code{\link{impute_AR1_Gaussian}} and \code{\link{impute_AR1_t}}), then
#'                  it will indicate the imputed values with red dots.
#' @param column Positive integer indicating the column index to be plotted (only valid if 
#'               the argument \code{y_imputed} is coercible to a matrix with more than one 
#'               column). Default is \code{1}.
#' @param type Type of plot. Valid options: \code{"ggplot2"} and \code{"simple"}. Default is 
#'             \code{"ggplot2"} (the package \code{ggplot2} must be installed).
#' @param title Title of the plot (default is \code{"red"}).
#' @param color_imputed Color for the imputed values (default is \code{"Imputed time series"}).
#' 
#' @author Daniel P. Palomar and Junyan Liu
#' 
#' @examples
#' library(imputeFin)
#' data(ts_AR1_t) 
#' y_missing <- ts_AR1_t$y_missing
#' y_imputed <- impute_AR1_t(y_missing)
#' plot_imputed(y_missing)
#' plot_imputed(y_imputed)
#' 
#' @import zoo
#' @export
plot_imputed <- function(y_imputed, column = 1, type = c("ggplot2", "simple"), title = NULL, color_imputed = "red") {
  any_index_miss <- !is.null(attr(y_imputed, "index_miss"))
  if (is.null(title)) title <- "Imputed time series"
  
  # extract the column to be plotted
  if (NCOL(y_imputed) > 1) {
    y_imputed_i <- y_imputed[, column]
    if (any_index_miss)
      index_miss <- attributes(y_imputed)$index_miss[[column]]
  } else {
    y_imputed_i <- y_imputed
    if (any_index_miss)
      index_miss <- attributes(y_imputed)$index_miss
  }
  
  index_y <- index(y_imputed_i)
  switch(match.arg(type),
         "simple" = {
           #p <- plot(index_y, y_imputed_i, type = "l",  col = "black", xlab = "", ylab = "", main = title)
           #grid()
           p <- plot(y_imputed_i)
           if (any_index_miss) {
             #p <- lines(index_y[index_miss], y_imputed_i[index_miss], col = color_imputed)
             # points(index_y[index_miss], y_imputed_i[index_miss], pch = 19, col = "red")
             # legend(x = "topleft", legend = "imputed values", col = "red", pch = 19)
             p <- lines(y_imputed_i[index_miss], col = color_imputed)
           }
           p
         },
         "ggplot2" = {
           if (!requireNamespace("ggplot2", quietly = TRUE)) 
             stop("Please install package \"ggplot2\" or choose another plot type.", call. = FALSE)
           data_frm  <- data.frame(index = index_y, value = as.vector(y_imputed_i))
           # ggplot2::ggplot() +
           #   ggplot2::geom_line(data = data_frm, ggplot2::aes(x = index, y = value), col = "black") +
           #   ggplot2::geom_point(data = data_frm[index_miss, ], ggplot2::aes(x = index, y = value), col = "#DD3344") +
           #   ggplot2::labs(title = "Imputed time series", x = "time")  # + ggplot2::theme_bw()
           # http://zevross.com/blog/2014/08/04/beautiful-plotting-in-r-a-ggplot2-cheatsheet-3/#manually-adding-legend-items-guides-override.aes
           p <- ggplot2::ggplot() +
             ggplot2::geom_line(data = data_frm, ggplot2::aes_string(x = "index", y = "value"), col = "black") +
             ggplot2::labs(title = title, x = NULL, y = NULL)
             #ggplot2::scale_x_date(date_breaks = "6 months", date_labels = "%b %Y")
           if (any_index_miss)
             p <- p + ggplot2::geom_line(data = data_frm[index_miss, ], ggplot2::aes(x = index, y = value), col = color_imputed)
             #ggplot2::geom_point(data = data_frm[index_miss, ], ggplot2::aes(x = index, y = value), col = "red", size = 0.1)
             # p <- p + ggplot2::geom_point(data = data_frm[index_miss, ], ggplot2::aes(x = index, y = value, col = "imputed values")) +
             #   ggplot2::theme(legend.title = ggplot2::element_blank(), legend.justification = c(0.9, 0), legend.position = c(0.98, 0.02)) +
             #   ggplot2::scale_colour_manual(name = "", values = c("imputed values" = "#DD3344"))
           p
         },
         stop("Plot type unknown"))
}
