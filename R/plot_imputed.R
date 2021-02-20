#' @title Plot imputed time series.
#'
#' @description Plot single imputed time series (as returned by functions 
#'              \code{\link{impute_AR1_Gaussian}} and \code{\link{impute_AR1_t}}),
#'              highlighting the imputed values in a different color.
#'
#' @param y_imputed Imputed time series (can be any object coercible to a numeric vector 
#'                  or a numeric matrix). If it has the attribute \code{"index_miss"} (as
#'                  returned by any of the imputation functions 
#'                  \code{\link{impute_AR1_Gaussian}} and \code{\link{impute_AR1_t}}), then
#'                  it will highlight the imputed values in a different color.
#' @param column Positive integer indicating the column index to be plotted (only valid if 
#'               the argument \code{y_imputed} is coercible to a matrix with more than one 
#'               column). Default is \code{1}.
#' @param type Type of plot. Valid options: \code{"ggplot2"} and \code{"simple"}. Default is 
#'             \code{"ggplot2"} (the package \code{ggplot2} must be installed).
#' @param title Title of the plot (default is \code{"Imputed time series"}).
#' @param color_imputed Color for the imputed values (default is \code{"red"}).
#' 
#' @author Daniel P. Palomar
#' 
#' @examples
#' library(imputeFin)
#' data(ts_AR1_t) 
#' y_missing <- ts_AR1_t$y_missing
#' y_imputed <- impute_AR1_t(y_missing)
#' plot_imputed(y_missing, title = "Original time series with missing values")
#' plot_imputed(y_imputed)
#' 
#' @import zoo
#' @import graphics
#' @export
plot_imputed <- function(y_imputed, column = 1, 
                         title = "Imputed time series", color_imputed = "red", type = c("ggplot2", "simple")) {

  # extract the column to be plotted
  index_miss <- NULL
  any_index_miss <- !is.null(attr(y_imputed, "index_miss"))
  if (NCOL(y_imputed) > 1) {
    if (any_index_miss)
      index_miss <- attributes(y_imputed)$index_miss[[column]]
    y_imputed <- y_imputed[, column]
  } else {
    if (any_index_miss)
      index_miss <- attributes(y_imputed)$index_miss
  }

  # obtain some indices convenient for plotting
  # 1) separate missing values into isolated and nonisolated (containing L neighbors on the left or right)
  index_miss_nonisolated <- NULL
  L <- 2
  if (any_index_miss)
    for (i in index_miss)
      if (all(c(i-L, i) %in% index_miss) || all(c(i, i+L) %in% index_miss))  # any neighbor
        index_miss_nonisolated <- c(index_miss_nonisolated, i)
  index_miss_isolated <- setdiff(index_miss, index_miss_nonisolated)
  # 2) obtain expanded indexes for the plot
  index_miss_expanded <- index_miss
  for (i in index_miss) {
    if (i > 1)
      index_miss_expanded <- union(index_miss_expanded, i - 1)
    if (i < length(y_imputed))
      index_miss_expanded <- union(index_miss_expanded, i + 1)
  }
  index_miss_expanded <- sort(index_miss_expanded)
  
  # plot
  switch(match.arg(type),
         "simple" = {
           #p <- plot(index_y, y_imputed, type = "l",  col = "black", xlab = "", ylab = "", main = title)
           #grid()
           p <- plot(y_imputed, main = title, xlab = "", ylab = "")
           if (any_index_miss) {
             p <- points(y_imputed[index_miss], pch = 20, size = 5, col = color_imputed)
             #p <- lines(index_y[index_miss], y_imputed[index_miss], col = color_imputed)
             # points(index_y[index_miss], y_imputed[index_miss], pch = 20, col = "red")
             # legend(x = "topleft", legend = "imputed values", col = "red", pch = 20)
           }
           invisible(p)
         },
         "ggplot2" = {
           if (!requireNamespace("ggplot2", quietly = TRUE)) 
             message("Please install package \"ggplot2\" or choose a basic plot type with the argument type = \"simple\".")
           else {
             df_all  <- data.frame(index = index(y_imputed), value = as.vector(y_imputed))
             df_obs  <- df_all; df_obs$value[index_miss] <- NA
             df_miss <- df_all; df_miss$value[setdiff(1:nrow(df_all), index_miss_expanded)] <- NA
             
             index <- value <- NULL  # dirty hack to avoid CRAN note
             p <- ggplot2::ggplot() +
               ggplot2::geom_line(data = df_miss, ggplot2::aes(x = index, y = value), col = color_imputed) +
               ggplot2::geom_line(data = df_obs, ggplot2::aes_string(x = "index", y = "value"), col = "black") +
               ggplot2::geom_point(data = df_all[index_miss_isolated, ], ggplot2::aes(x = index, y = value), col = color_imputed, size = 0.8) +
               #ggplot2::geom_point(data = df_all[index_miss, ], ggplot2::aes(x = index, y = value), col = color_imputed, size = 0.8) +
               #ggplot2::scale_x_date(date_breaks = "6 months", date_labels = "%b %Y") +
               ggplot2::labs(title = title, x = NULL, y = NULL)
             
             # p <- ggplot2::ggplot() +
             #   ggplot2::geom_line(data = df_all, ggplot2::aes_string(x = "index", y = "value"), col = "black") +
             #   ggplot2::labs(title = title, x = NULL, y = NULL)
             #   #ggplot2::scale_x_date(date_breaks = "6 months", date_labels = "%b %Y")
             # if (any_index_miss) {
             #   p <- p + ggplot2::geom_point(data = df_all[index_miss_isolated, ], ggplot2::aes(x = index, y = value), col = color_imputed, size = 0.8)
             #   if (!is.null(index_miss_nonisolated)) {
             #     df_all$value[-index_miss_nonisolated] <- NA
             #     p <- p + ggplot2::geom_line(data = df_all, ggplot2::aes(x = index, y = value), col = color_imputed)
             #   }
             # }
             
             suppressWarnings(p)
           }
         },
         stop("Plot type unknown"))
}

# http://zevross.com/blog/2014/08/04/beautiful-plotting-in-r-a-ggplot2-cheatsheet-3/#manually-adding-legend-items-guides-override.aes
# p <- p + ggplot2::geom_point(data = data_frm[index_miss, ], ggplot2::aes(x = index, y = value, col = "imputed values")) +
#   ggplot2::theme(legend.title = ggplot2::element_blank(), legend.justification = c(0.9, 0), legend.position = c(0.98, 0.02)) +
#   ggplot2::scale_colour_manual(name = "", values = c("imputed values" = "#DD3344"))

