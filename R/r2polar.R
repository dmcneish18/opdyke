#' @title Convert correlation matrix to polar angles
#'
#' @description This function transforms a symmetric correlation matrix to a lower triangle of polar angles
#' Input can be a raw correlation matrix or a lavaan object.
#'
#' @param cor An object containing the observed correlation matrix. Either a raw correlation or a lavaan object. Object type is specified in the 'type' option.
#' this is a required argument for the function. 
#' @param type Specifies the type of object in the cor argument. The default is type='lavaan'. The other acceptable option is type='cor' if the object is a 
#' raw correlation matrix. 
#'
#' @import dplyr lavaan 
#'
#' @author Daniel McNeish
#'
#' Maintainer: Daniel McNeish <dmcneish@asu.edu>
#'
#' @rdname r2polar
#'
#' @return Lower triangle of polar angles. 
#' @export
#'
#' @examples
#' #assign dataset to “dat”
#' dat<-lavaan::HolzingerSwineford1939
#' #assign correlation matrix to “A”
#' R<-cor(dat[,7:12])
#' #apply r2polar function with type= “cor”
#' r2polar(cor=R, type="cor")

r2polar<-function(cor=NULL, type="lavaan"){
  
  if (!(type %in% c("lavaan", "cor"))) {
    stop("Error: r2polar only accepts lavaan objects or correlation matrices.
         
         Please use type='lavaan' (the default) or type='cor'")
  }
  
  if(type=="lavaan"){
    S1<-lavaan::lavInspect(cor, what="sampstat")
    S<-cov2cor(matrix(S1$cov, nrow=nrow(S1$cov),ncol=nrow(S1$cov)))
  }
  
  if(type=="cor"){
    S<-cor
    
    if (!(all(diag(S) == rep(1,nrow(S))))){
      stop("Error: input matrix does not have diagonals equal to 1")
    }
    
    if (!(all(abs(S[lower.tri(S, diag=F)])<1))){
      stop("Error: input matrix has off-diagonal terms with absolute value greater than 1")
    }
  }
  
  A3<-t(chol(S))
  A4<-matrix(0, nrow=nrow(A3), ncol=ncol(A3))
  # get polar coordinates from root of observed
  A4[1,1]<-1
  for (i in 2:ncol(A3)){
    A4[i,1]=acos(A3[i,1])
    p<-cumprod(asin(A3[i,]))
    q<-cumprod(sin(A3[i,]))
    #A4[i,i]=p[i-1]
    for(j in 2:ncol(A3)){
      if(j<i) { 
        q<-cumprod(sin(A4[i,]))
        A4[i,j]<-acos(A3[i,j]/q[j-1])   
      }
    }
  }
  return(A4)
}