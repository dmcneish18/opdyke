#' @title Convert polar angles to a correlation matrix
#'
#' @description This function transforms a lower triangle matrix of polar angles to a symmetric correlation matrix.
#' Input is a lower triangle matrix of polar angles.  
#'
#' @param polar A lower triangle matrix of polar angles in radians where each element is between 0 and pi.  
#'
#' @author Daniel McNeish
#'
#' Maintainer: Daniel McNeish <dmcneish@asu.edu>
#'
#' @rdname polar2r
#'
#' @return Correlation matrix 
#' @export
#'
#' @examples
#' #assign dataset to “dat”
#' dat<-lavaan::HolzingerSwineford1939
#' #assign correlation matrix to “A”
#' R<-cor(dat[,7:12])
#' #create lower triangle matrix of polar angles
#' P<-r2polar(cor=R, type="cor")
#' #use polar2r to transform the polar angles back to a correlation matrix
#' polar2r(polar=P)

polar2r<-function(polar=NULL){
  
  if(!(all(ifelse(polar[lower.tri(polar)] >= 0& polar[lower.tri(polar)]<=pi,TRUE,FALSE)))){
    stop("Error: Input polar angles are outside of [0,pi]")
  }
  
  if(!(all(ifelse(polar[upper.tri(polar)]==0,TRUE, FALSE)))){
    stop("Error: Input polar angle matrix is not lower triangle")
  }
  
  A5<-matrix(0,nrow=nrow(polar),ncol=ncol(polar))
  A5[1,1]<-1
  for (i in 2:ncol(A5)){
    A5[i,1]=cos(polar[i,1])
    p<-cumprod(sin(polar[i,]))
    A5[i,i]=p[i-1]
    for(j in 2:ncol(A5)){
      if(j<i) {  
        A5[i,j]<-cos(polar[i,j])*p[j-1]   
      }
    }
    
  }
  
  A6<-A5%*%t(A5)
  return(A6)
}