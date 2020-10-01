#' Draws prediction ellipses (2D)
#'
#' This function draws the prediction ellipses for the input data onto an existing plot.
#'
#' @param mu mean vector of the normal distribution.
#' @param sig covariance matrix of the normal distribution.
#' @param alphas quantiles of the prediction ellipses to be drawn.
#' @param color color of the prediction ellipses.
#' @param line_width line width of the prediction ellipses.
#' @param line_type line type of the prediction ellipses.
#' @return None (prediction ellipses polygon).
#' @export
#' @examples
#' library("mvtnorm")
#' s1 <- 2
#' s2 <- 4
#' rho<- 3/4
#'
#' n <- 500
#' m <- c(0,0)
#' s <- matrix(c(s1^2,rho*s1*s2,rho*s1*s2,s2^2), ncol=2)
#' X <- rmvnorm(n, m, s, method="chol")
#'
#' plot(X, main="", xlab="", ylab="", xlim=c(-10,10), ylim=c(-15,15), type="n")
#' abline(h=-3:3*5, v=-3:3*5, col="grey", lty=2)
#' abline(h=0, v=0)
#' points(X, col="green", pch=20)
#'
#' drawPredEllipse(mu=m, sig=s, alphas=c(0.5,0.75), color="blue")
#' drawPredEllipse(mu=m, sig=s, alphas=c(0.25,0.9), color="red", line_width=1, line_type=1)
drawPredEllipse <- function(mu, sig, alphas=c(0.10,0.50,0.99), color="blue", line_width=2, line_type=2){
  if(length(mu) != 2)
    return("Mean vector is not valid!")

  if(!(is.matrix(sig) && isSymmetric.matrix(sig)))
    return("Covariance matrix is not valid!")

  if(!all(alphas<1 & alphas>0))
    return("Alphas are not valid!")

  for(i in alphas){
    r <- eigen(sig, symmetric=TRUE)
    angle <- atan(r$vectors[2,1]/r$vectors[1,1])
    a <- sqrt(r$values[1]*qchisq(i, df=2))
    b <- sqrt(r$values[2]*qchisq(i, df=2))

    t <- seq(0, 2*pi, by=pi/512)
    x <- mu[1] + a*cos(t)*cos(angle) - b*sin(t)*sin(angle)
    y <- mu[2] + a*cos(t)*sin(angle) + b*sin(t)*cos(angle)

    lines(x, y, col=color, lwd=line_width, lty=line_type)
  }
}
