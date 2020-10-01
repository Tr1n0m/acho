#' Draws convex hull (2D)
#'
#' This function draws the convex hull for the input data onto an existing plot.
#'
#' @param data input data.
#' @param color color of the convex hull lines.
#' @param line_width line width of the convex hull.
#' @param line_type line type of the convex hull.
#' @return None (convex hull polygon).
#' @export
#' @examples
#' n <- 100
#' X <- cbind(rnorm(n, 0, 2),rnorm(n,0,1))
#'
#' plot(X, xlab="", ylab="", pch=20, bty="n")
#' abline(h=0, v=0)
#' drawCHull(X)
#'
#' plot(X, xlab="", ylab="", pch=20, bty="n")
#' abline(h=0, v=0)
#' drawCHull(X, color="red", line_width=2, line_type=2)
drawCHull <- function(data, color="blue", line_width=1, line_type=1){
  if(dim(data)[2] != 2)
    return("Data has an incorrect number of dimensions!")

  chull_pts <- chull(data[,1],data[,2])
  chull_pts <- c(chull_pts,chull_pts[1])

  lines(data[chull_pts,1], data[chull_pts,2], col=color, lwd=line_width, lty=line_type)
}
