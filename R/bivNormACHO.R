#' The ACHO function for given distributions (2D)
#'
#' This function calculates a value similar to the "Averaged Convex Hull Overlap" based on the prediction ellipses of the given normal distributions.
#'
#' @param mu1 mean vector of the first set.
#' @param sig1 covariance matrix for the first set.
#' @param mu2 mean vector of the second set.
#' @param sig2 covariance matrix for the second set.
#' @return The general "ACHO" value of the normal distributions.
#' @export
#' @examples
#' library("mvtnorm")
#' n <- 300
#' s1_x <- 2
#' s2_x <- 4
#' rho_x<- 3/4
#' m_x <- c(0,0)
#' s_x <- matrix(c(s1_x^2,rho_x*s1_x*s2_x,rho_x*s1_x*s2_x,s2_x^2), ncol=2)
#' X <- rmvnorm(n, m_x, s_x, method="chol")
#'
#' s1_y <- 2
#' s2_y <- 4
#' rho_y<- -3/4
#' m_y <- c(2,3)
#' s_y <- matrix(c(s1_y^2,rho_y*s1_y*s2_y,rho_y*s1_y*s2_y,s2_y^2), ncol=2)
#' Y <- rmvnorm(n, m_y, s_y, method="chol")
#'
#' plot(X, main="", xlab="", ylab="", xlim=c(-8,8), ylim=c(-15,15), type="n", bty="n")
#' abline(h=0, v=0)
#' points(X, col="blue", pch=20)
#' points(Y, col="red", pch=20)
#'
#' bivNormACHO(m_x, s_x, m_y, s_y)
bivNormACHO <- function(mu1, sig1, mu2, sig2){
  if(length(mu1) != 2 || length(mu2) != 2)
    return("Mean vectors are not valid!")

  if( !(is.matrix(sig1) && isSymmetric.matrix(sig1) && is.matrix(sig2) && isSymmetric.matrix(sig2)) )
    return("Covariance matrices are not valid!")

  alpha <- seq(0.01,0.99,by=0.01)
  l <- length(alpha)
  df <- data.frame(a1=rep(0,l))

  r <- eigen(sig1, symmetric=TRUE)
  df$a1 <- sqrt(r$values[1]*qchisq(alpha, df=2))
  df$b1 <- sqrt(r$values[2]*qchisq(alpha, df=2))
  df$area1 <- df$a1 * df$b1 * pi
  df$angle1 <- rep(atan(r$vectors[2,1]/r$vectors[1,1]),l)

  r <- eigen(sig2, symmetric=TRUE)
  df$a2 <- sqrt(r$values[1]*qchisq(alpha, df=2))
  df$b2 <- sqrt(r$values[2]*qchisq(alpha, df=2))
  df$area2 <- df$a2 * df$b2 * pi
  df$angle2 <- rep(atan(r$vectors[2,1]/r$vectors[1,1]),l)

  df$intersection <- rep(0,l)
  df$acho <- rep(0,l)
  i <- l
  while(i > 0){
    t <- seq(0,2*pi,by=pi/512)

    x1 <- mu1[1] + df$a1[i]*cos(t)*cos(df$angle1[1]) - df$b1[i]*sin(t)*sin(df$angle1[1])
    y1 <- mu1[2] + df$a1[i]*cos(t)*sin(df$angle1[1]) + df$b1[i]*sin(t)*cos(df$angle1[1])

    x2 <- mu2[1] + df$a2[i]*cos(t)*cos(df$angle2[1]) - df$b2[i]*sin(t)*sin(df$angle2[1])
    y2 <- mu2[2] + df$a2[i]*cos(t)*sin(df$angle2[1]) + df$b2[i]*sin(t)*cos(df$angle2[1])

    ellipse1 <- cbind(x1,y1)
    ellipse1 <- rbind(ellipse1,ellipse1[1,])
    ellipse1_poly <- sp::Polygon(ellipse1)
    ellipse1_poly <- sp::Polygons(list(ellipse1_poly), ID="Data1")
    ellipse1_spatial_poly <- sp::SpatialPolygons(list(ellipse1_poly))

    ellipse2 <- cbind(x2,y2)
    ellipse2 <- rbind(ellipse2,ellipse2[1,])
    ellipse2_poly <- sp::Polygon(ellipse2)
    ellipse2_poly <- sp::Polygons(list(ellipse2_poly), ID="Data2")
    ellipse2_spatial_poly <- sp::SpatialPolygons(list(ellipse2_poly))

    intersection <- rgeos::gIntersection(ellipse1_spatial_poly, ellipse2_spatial_poly, checkValidity=TRUE)
    if(!is.null(intersection)){
      df$intersection[i] <- unlist(sapply(slot(intersection,"polygons"),function(p) sapply(slot(p,"Polygons"),slot,"area")))
      df$acho[i] <- 2*df$intersection[i]/(df$area1[i] + df$area2[i])
    }
    else{
      return(mean(df$acho, na.rm=TRUE))
    }
    i <- i-1
  }

  return(mean(df$acho, na.rm=TRUE))
}
