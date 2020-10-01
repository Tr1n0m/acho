#Installation Polygon Intersection Libraries
#install.packages(c("sp",""rgeos"))

#Polygon
#library("sp")

#gIntersections
#gArea
#library("rgeos")

reduceDataBoundery <- function(data, reduction){
  #reduction only allowed between 0 and 1
  if(reduction <= 0){
    return(c())
  }
  if(reduction >= 1){
    return(data)
  }

  #1 dimensional case
  if(is.null(dim(data)) || dim(data)[2]<2){
    reduced_data <- data
    len <- length(reduced_data)
    threshold <- reduction*len

    while(len > threshold){
      #reduced to the point of no data left
      if(len <= 2){
        return(c())
      }

      bounds <- c(which(reduced_data==min(reduced_data)), which(reduced_data==max(reduced_data)))
      reduced_data <- reduced_data[-bounds]
      len <- len-2
    }

    return(reduced_data)
  }

  #2 dimensional case
  reduced_data <- data
  len <- length(reduced_data[,1])
  threshold <- reduction*len
  r <- reduced_data[,1]

  while(length(r) > threshold){
    #reduced to the point of no data left
    if(is.null(dim(reduced_data)) || length(r)==0){
      return(c())
    }

    c <- chull(reduced_data)
    reduced_data <- reduced_data[-c,]
    r <- r[-c]
  }

  return(reduced_data)
}

getAreaOfIntersection <- function(data1, data2, alpha){
  #Sorensen-Dice overlap
  dice <- 0

  #reducing both datasets via alpha
  data1 = reduceDataBoundery(data1, alpha)
  data2 = reduceDataBoundery(data2, alpha)

  #reduction caused one dataset to be to small
  if(length(data1) <= 1 || length(data2) <= 1){
    return(0)
  }

  # 1 dimensional data
  if(is.null(dim(data1)) || dim(data1)[2]<2){

    data1_min <- min(data1)
    data1_max <- max(data1)
    data2_min <- min(data2)
    data2_max <- max(data2)

    if(data1_min==data1_max || data2_min==data2_max){
      return(dice)
    }

    if(data2_min >= data1_max || data1_min >= data2_max){
      return(dice)
    }
    else{
      o_min <- max(data1_min,data2_min)
      o_max <- min(data1_max,data2_max)

      data1_vol <- data1_max - data1_min
      data2_vol <- data2_max - data2_min

      #jaccard <- (o_max-o_min)/(x_vol+y_vol-(o_max-o_min))
      dice <- 2*(o_max-o_min)/(data1_vol+data2_vol)

      return(dice)
    }
  }

  #reduction caused one dataset to be to small for 2D
  if(is.null(dim(data1)[1]) || is.null(dim(data2)[1])){
    return(0)
  }
  if(dim(data1)[1] <= 2 || dim(data2)[1] <= 2){
    return(0)
  }

  #calculate the convex hulls
  data1_chull <- data1[chull(data1[,1],data1[,2]),]
  data1_chull <- rbind(data1_chull,data1_chull[1,])
  data2_chull <- data2[chull(data2[,1],data2[,2]),]
  data2_chull <- rbind(data2_chull,data2_chull[1,])

  #calculate the spatial polygons of the convex hulls
  data1_poly <- sp::Polygon(data1_chull)
  data1_poly <- sp::Polygons(list(data1_poly), ID="Data1")
  data1_spatial_poly <- sp::SpatialPolygons(list(data1_poly))

  data2_poly <- sp::Polygon(data2_chull)
  data2_poly <- sp::Polygons(list(data2_poly), ID="Data2")
  data2_spatial_poly <- sp::SpatialPolygons(list(data2_poly))

  #calculate the intersection of the convex hulls
  intersection_of_chulls <- rgeos::gIntersection(data1_spatial_poly, data2_spatial_poly, checkValidity = TRUE)

  #does the intersection exist?
  if(!is.null(intersection_of_chulls)){
    #extract the area and the border coordinates of the intesection
    intersection_vol <- unlist(sapply(slot(intersection_of_chulls,"polygons"),function(p) sapply(slot(p,"Polygons"),slot,"area")))

    data1_vol <- rgeos::gArea(data1_spatial_poly)
    data2_vol <- rgeos::gArea(data2_spatial_poly)

    dice <- 2*intersection_vol/(data1_vol + data2_vol)

    return(dice)
  }
  return(dice)
}

#' The base ACHO function
#'
#' This function calculates the "Averaged Convex Hull Overlap" for the input and can also show a plot for the behaviour of the reducing overlaps.
#'
#' @param data1 first dataset (1D/2D).
#' @param data2 second dataset (1D/2D).
#' @param acho_plot Determines if the convex hull overlap behaviour plot is shown. Default TRUE.
#' @return The "ACHO" value of the datasets.
#' @export
#' @examples
#' n <- 1000
#' x1 <- rnorm(n, 0, 1)
#' x2 <- rnorm(n, 1, 1.5)
#'
#' plot(density(x1), main="", xlab="", ylab="", type="n", bty="n")
#' lines(density(x1), col="blue")
#' lines(density(x2), col="red")
#'
#' getACHO(x1, x2)
#' getACHO(x1, x2, acho_plot=TRUE)
getACHO <- function(data1, data2, acho_plot=FALSE){
  #minimal length of the data
  len_min <- 0
  if(is.null(dim(data1)) || dim(data1)[2]<2){
    len_min <- min(length(data1),length(data2))
  }
  else{
    len_min <- min(length(data1[,1]),length(data2[,1]))
  }

  #dynamic averaging grid based on minimal length
  d <- 0
  if(len_min <= 100){
    d <- 10
  } else if(len_min <= 400){
    d <- 5
  } else{
    d <- 1
  }

  alpha <- seq(0,100,by=d)/100
  alpha <- alpha[-1]
  intersections <- rep(0,length(alpha))

  for(i in length(intersections):1){
    intersections[i] <- getAreaOfIntersection(data1, data2, alpha[i])

    #no further calculation needed once it drops to zero
    if(intersections[i] == 0) break
  }

  if(acho_plot){
    #ACHO visualization
    plot(alpha,intersections,col="blue",xlab="Alpha",ylab="Overlap",xlim=c(0,1),ylim=c(0,1),type="n",cex.lab=1.5,cex.axis=1.5)
    abline(h=1:10/10,v=1:10/10,col="gray",lty=3)
    abline(h=c(0,1),v=c(0,1),lwd=2)
    abline(h=mean(intersections,na.rm = TRUE),col="green",lwd=3)
    lines(alpha, intersections,col="blue",lwd=3)
  }

  return(mean(intersections,na.rm = TRUE))
}
