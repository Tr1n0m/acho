#' The ACHO function for given distributions (1D)
#'
#' This function calculates a value similar to the "Averaged Convex Hull Overlap" based on the centralized quantiles of the given normal distributions.
#'
#' @param mu1 mean vector of the first set.
#' @param sig1 standard deviation for the first set.
#' @param mu2 mean vector of the second set.
#' @param sig2 standard deviation for the second set.
#' @return The general "ACHO" value of the normal distributions.
#' @export
#' @examples
#' n <- 1000
#' s <- seq(-5, 6, by=0.01)
#'
#' m1 <- 0
#' s1 <- 1
#' x1 <- dnorm(s, m1, s1)
#' m2 <- 1
#' s2 <- 1.5
#' x2 <- dnorm(s, m2, s2)
#'
#' plot(s, x1, xlab="", ylab="", type="n", bty="n")
#' lines(s, x1, col="blue")
#' lines(s, x2, col="red")
#'
#' normACHO(m1, s1, m2, s2)
normACHO <- function(mu1, sig1, mu2, sig2){
  if( !(sig1>0 && sig2>0) )
    return("Standard deviations are not valid!")

  lvl <- 1:99/200
  q1 <- qnorm(lvl, mu1, sig1)
  q2 <- qnorm(lvl, mu2, sig2)

  df <- data.frame(min1=q1, max1=2*mu1-q1, min2=q2, max2=2*mu2-q2)
  df$range1 <- df$max1 - df$min1
  df$range2 <- df$max2 - df$min2

  df$intsection <- pmin(df$max1,df$max2) - pmax(df$min1,df$min2)
  df$intsection <- ifelse(df$intsection < 0, 0, df$intsection)
  df$cho <- 2*df$intsection/(df$range1 + df$range2)

  return(mean(df$cho))
}
