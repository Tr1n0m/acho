#' Plots convex hull (2D)
#'
#' This function plots the data and the convex hull for it.
#'
#' @param data input data.
#' @param color color of the convex hull lines.
#' @param points_color color of the data points.
#' @param line_width line width of the convex hull.
#' @param line_type line type of the convex hull.
#' @return Plot with data points and their convex hull.
#' @export
#' @examples
#' n <- 100
#' X <- cbind(rnorm(n, 0, 2),rnorm(n,0,1))
#'
#' plotCHull(X)
#' plotCHull(X, color="red", points_color="blue", line_width=2, line_type=2)
plotCHull <- function(data, color="blue", points_color="black", line_width=1, line_type=1){
  if(dim(data)[2] != 2)
    return("Data has an incorrect number of dimensions!")

  chull_pts <- chull(data[,1],data[,2])
  chull_pts <- c(chull_pts,chull_pts[1])

  plot(data, col=points_color, xlab="", ylab="", xlim=c(min(data[,1]),max(data[,1])), ylim=c(min(data[,2]),max(data[,2])), pch=20, bty="n")
  abline(h=0, v=0)
  lines(data[chull_pts,1], data[chull_pts,2], col=color, lwd=line_width, lty=line_type)
}
