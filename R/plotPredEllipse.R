#' Plot prediction ellipses (2D)
#'
#' This function plots the data and the prediction ellipses for it.
#'
#' @param mu mean vector of the normal distribution.
#' @param sig covariance matrix of the normal distribution.
#' @param alphas quantiles of the prediction ellipses to be drawn.
#' @param color color of the prediction ellipses.
#' @param points_color color of the data points.
#' @param line_width line width of the prediction ellipses.
#' @param line_type line type of the prediction ellipses.
#' @return Plot with data points and their prediction ellipses.
#' @export
#' @examples
#' s1 <- 2
#' s2 <- 4
#' rho<- 3/4
#'
#' m <- c(3,5)
#' s <- matrix(c(s1^2,rho*s1*s2,rho*s1*s2,s2^2), ncol=2)
#'
#' plotPredEllipse(mu=m, sig=s)
#' plotPredEllipse(mu=m, sig=s, alphas=c(0.25,0.5,0.75,0.9), color="red", line_width=1, line_type=1)
plotPredEllipse <- function(mu, sig, alphas=c(0.1,0.50,0.99), color="blue", points_color="black", line_width=2, line_type=2){
  if(length(mu) != 2)
    return("Mean vector is not valid!")

  if(!(is.matrix(sig) && isSymmetric.matrix(sig)))
    return("Variance matrix is not valid!")

  if(!all(alphas<1 & alphas>0))
    return("Alphas are not valid!")

  for(i in sort(alphas, decreasing=TRUE)){
    r <- eigen(sig, symmetric=TRUE)
    angle <- atan(r$vectors[2,1]/r$vectors[1,1])
    a <- sqrt(r$values[1]*qchisq(i, df=2))
    b <- sqrt(r$values[2]*qchisq(i, df=2))

    t <- seq(0, 2*pi, by=pi/512)
    x <- mu[1] + a*cos(t)*cos(angle) - b*sin(t)*sin(angle)
    y <- mu[2] + a*cos(t)*sin(angle) + b*sin(t)*cos(angle)

    if(i == max(alphas)){
      plot(x, y, col=points_color, main="", xlab="", ylab="", xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), type="n")
      abline(h=0, v=0)
    }
    lines(x, y, col=color, lwd=line_width, lty=line_type)
  }
}
