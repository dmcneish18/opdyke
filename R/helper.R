#' @import matrixcalc lavaan dplyr
#' @importFrom stats cov2cor
#' @importFrom gsl hyperg_2F1
#' @importFrom pracma csc cot

## Polar decomposition;
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

##### polar coordinates to correlation matrix
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

##Opdyke PDF and CDF
opdyke<-function(cor=NULL, type="lavaan", precision="less"){

  if (!(type %in% c("lavaan", "cor"))) {
    stop("Error: opdyke only accepts lavaan objects or correlation matrices.

         Please use type='lavaan' (the default) or type='cor'")
  }

  if (!(precision %in% c("less", "more"))) {
    stop("Error: precision must be either 'less' or 'more'.

          'less' uses two-decimal point accuracy for polar angles
          'more' uses three-decimal point accuracy for polar angles")
  }

  if(type=="lavaan"){
    S1<-lavaan::lavInspect(cor, what="sampstat")
    A<-cov2cor(matrix(S1$cov, nrow=nrow(S1$cov),ncol=nrow(S1$cov)))
  }
  if(type=="cor"){
    A<-cor

    if (!(matrixcalc::is.positive.definite(A))){
      stop("Error: input matrix is not positive definite.")
    }

    if (!(all(diag(A) == rep(1,nrow(A))))){
      stop("Error: input matrix does not have diagonals equal to 1")
    }

    if (!(all(abs(A[lower.tri(A, diag=F)])<1))){
      stop("Error: input matrix has off-diagonal terms with absolute value greater than 1")
    }

  }

  I<-c(1:ncol(A))

  #all permutations slotted into lower right
  I5<-matrix(0,ncol=length(I),nrow=1)
  for (r2 in 1:(length(I)-1)){
    for(r1 in 0:(r2-1)){

      I1<-c((length(I)-r2),length(I)-r1)
      I2<-I[-c((length(I)-r2),length(I)-r1)]

      I4<-c(I2,I1)
      I5<-rbind(I5,I4)
    }
  }
  com3<-I5[-1,]
  colnames(com3)<-colnames(A)

  pdf<-list()
  for(c in 1:nrow(com3)){
    A1<-A[com3[c,],com3[c,]]
    X<-r2polar(A1,type="cor")
    k=ncol(X)-1
    theta<-X[ncol(X),k]
    ck<-gamma((.5*k+1))/(sqrt(pi)*gamma(.5*k+.5))
    d<-data.frame()

    if(precision=="more"){
      dummy<-10
    }

    if(precision=="less"){
      dummy<-1
    }

    for (r in 1:(314*dummy)){
      z=r/(100*dummy)

      d[r,1]<-z
      d[r,2]<-ck*pracma::csc(z)^2*(1+(pracma::cot(z)-pracma::cot(theta))^2)^(-1-k/2)

      X2<-X
      X2[ncol(X),k]<-z
      Q<-polar2r(X2)
      d[r,3]<-Q[ncol(X),k]
      d[r,4]<-1-(0.5+ck*(pracma::cot(theta)-pracma::cot(z))*gsl::hyperg_2F1(a=.5,b=(1+0.5*k),c=1.5,x=-1*(pracma::cot(z)-pracma::cot(theta))^2))
    }
    cdf<-cumsum(d[,1])
    d[,5]<-1-(cdf/cdf[length(cdf)])
    d[,6]<-theta
    d[,7]<-com3[c,ncol(com3)]
    d[,8]<-com3[c,ncol(com3)-1]


    colnames(d)<-c("z","pdf","r","cdf_exact","cdf_approx","theta","row","column")
    pdf[[c]]<-d
  }
  return(pdf)
}
