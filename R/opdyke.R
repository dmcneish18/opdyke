#' @title Calculate the Opdyke distribution PDF and CDF  for each element in a correlation matrix.
#'
#' @description This function calculates the probability density function and cumulative distribution function of the Opdyke distribution for
#' value between 0 and pi. The output is a list where each list element is a different correlation matrix from the data. This function is useful for creating
#' custom plots or custom calculations that involves the Opdyke distribution.
#'
#' @param cor An object containing the observed correlation matrix. Either a raw correlation or a lavaan object. Object type is specified in the 'type' option.
#' this is a required argument for the function.
#' @param type Specifies the type of object in the cor argument. The default is type='lavaan'. The other acceptable option is type='cor' if the object is a
#' raw correlation matrix.
#' @param precision Controls the precision of the probability density function and cumulative distribution function calculations. The default is precision= “less”
#'  which calculates the PDF and CDF for polar angles between (0,pi) in .01 increments. precision= “more” calculates the PDF and CDF for polar angles between (0,pi)
#'  in .001 increments, which takes considerably longer, especially if there are many correlation elements.
#'
#' @import lavaan matrixcalc
#' @importFrom gsl hyperg_2F1
#' @importFrom pracma csc cot
#'
#' @author Daniel McNeish
#'
#' Maintainer: Daniel McNeish <dmcneish@asu.edu>
#'
#' @rdname opdyke
#'
#' @return list of Opdyke distribution PDF and CDF values for each element of a correlation matrix
#' @export
#'
#' @examples
#' #assign dataset to “dat”
#' dat<-lavaan::HolzingerSwineford1939
#' #assign correlation matrix to “A”
#' R<-cor(dat[,7:12])
#' #create lower triangle matrix of polar angles
#' D<-opdyke::opdyke(cor=R, type="cor")
#' #print rows 100 to 125 of the 10th list element
#' D[[10]][100:125,]

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
