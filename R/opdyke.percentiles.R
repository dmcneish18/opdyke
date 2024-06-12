#' @title Calculate the Opdyke percentile of each predicted correlation from a structural equation model
#'
#' @description This function aids is determining the discrepancy between an observed correlation and a predicted correlation from a structural equation model.
#' The function calculates the percentile of the predicted correlation within the observed correlation's Opdyke distribution. The output shows elements with
#' large discrepancies, the percentiles associated with all elements, the bound of each element for which the correlation matrix remains positive definite, and the
#' interval of correlation residuals that are between the range specified by the "upper" and "lower" arguments.
#'
#' @param fit a lavaan object containing results from a fitted structural equation model. This will contain the observed and predicted covariance
#' and correlation matrices needed to calculate Opdyke percentiles.
#' @param precision Controls the precision of the probability density function and cumulative distribution function calculations. The default is precision= “less”
#'  which calculates the PDF and CDF for polar angles between (0,pi) in .01 increments. precision= “more” calculates the PDF and CDF for polar angles between (0,pi)
#'  in .001 increments, which takes considerably longer, especially if there are many correlation elements.
#' @param lower The minimum percentile of the Opdyke distribution that is considered to be acceptably close to the observed correlation (represented by the Opdyke
#'  distribution median). The percentile of the predicted correlation will be compared to this value in the simplified output to determine if the model reasonably
#'  reproduces an observed correlation. The default is 0.40.
#' @param upper the maximum percentile of the Opdyke distribution that is considered to be acceptably close to the observed correlation (represented by the Opdyke
#' distribution median). The percentile of the predicted correlation will be compared to this value in the simplified output to determine if the model reasonably
#' reproduces an observed correlation. The default is 0.60.
#'
#' @import lavaan
#'
#' @author Daniel McNeish
#'
#' Maintainer: Daniel McNeish <dmcneish@asu.edu>
#'
#' @rdname opdyke.percentiles
#'
#' @return Opdyke percentiles of predicted correlations from a structural equation model
#' @export
#'
#' @examples
#' #assign dataset to “dat”
#' dat<-lavaan::HolzingerSwineford1939
#' #define a two-factor model with correlated factors
#' model<-"f1=~ x1+x2+x3
#' f2 =~ x4+x5+x6
#' f1~~f2"
#' #fit the two-factor model in lavaan
#' f<-lavaan::cfa(data=dat, model=model)
#' #calculate Opdyke percentiles for each predicted correlation element
#' OP<-opdyke::opdyke.percentiles(fit=f, lower=.40, upper=.60)

opdyke.percentiles<-function(fit=NULL,lower=0.40,upper=0.60, precision="less"){

  if (!(precision %in% c("less", "more"))) {
    stop("Error: precision must be either 'less' or 'more'.

          'less' uses two-decimal point accuracy for polar angles
          'more' uses three-decimal point accuracy for polar angles")
  }

  if (upper<= lower) {
    stop("Error: uppper percentile bound must be greater than lower percentile bound")
  }

  if (upper>= 1) {
    stop("Error: uppper percentile bound must be less than 1.

         Upper bounds near 1 are not recommended.")
  }

  if (lower<= 0) {
    stop("Error: lower percentile bound must be greater than 0.

         Lower bounds near 0 are not recommended.")
  }

  #predicted
  P1<-lavaan::lavInspect(fit, what="cor.ov")
  P<-matrix(P1,nrow=nrow(P1),ncol=ncol(P1))
  colnames(P)<-colnames(P1)
  #observed
  S1<-lavaan::lavInspect(fit, what="sampstat")
  S<-cov2cor(matrix(S1$cov, nrow=nrow(P1),ncol=nrow(P1)))
  colnames(S)<-colnames(S1$cov)

  #get pdf
  b<-opdyke(fit, type="lavaan", precision=precision)

  #get cdf value for each correlation
  e<-data.frame()
  for(k in 1:length(b)){
    kk<-dplyr::filter(b[[k]][which.min(abs(b[[k]]$r-P[b[[k]]$row,b[[k]]$column])),])
    e<-rbind(e,kk)
  }

  #get lower %ile bound for each correlation (default is 40%)
  l<-data.frame()
  for(k in 1:length(b)){
    ll<-dplyr::filter(b[[k]][which.min(abs(b[[k]]$cdf_exact-lower)),])
    ll$l<-ll$r-S[b[[k]]$row[1],b[[k]]$column[1]]
    l<-rbind(l,ll)
  }

  #get upper %ile bound for each correlation (default is 60%)
  u<-data.frame()
  for(k in 1:length(b)){
    uu<-dplyr::filter(b[[k]][which.min(abs(b[[k]]$cdf_exact-upper)),])
    uu$u<-uu$r-S[b[[k]]$row[1],b[[k]]$column[1]]
    u<-rbind(u,uu)
  }

  #get lowerPD bound for each correlation
  rl<-data.frame()
  for(k in 1:length(b)){
    rll<-dplyr::filter(b[[k]][which.min(b[[k]]$r),])
    rl<-rbind(rl,rll)
  }

  #get upperrPD bound for each correlation
  ru<-data.frame()
  for(k in 1:length(b)){
    ruu<-dplyr::filter(b[[k]][which.max(b[[k]]$r),])
    ru<-rbind(ru,ruu)
  }

  #matrix of percentiles
  W<-matrix(NA,nrow=nrow(S),ncol=ncol(S))
  W[lower.tri(W,diag=F)]<-round(e[order(e$column,e$row),]$cdf_exact,2)
  W1<-W
  is.na(W) <-  W >= .4 & W <= .6
  colnames(W)<-colnames(S1$cov)
  rownames(W)<-colnames(S1$cov)

  #full
  WW<-as.table(W, na="", quote=F)
  colnames(WW)=rownames(WW)=colnames(S)

  #simplified
  WW1<-as.table(W1, na="", quote=F)
  colnames(WW1)=rownames(WW1)=colnames(S)

  #matrix of 40/60 intervals
  I<-matrix(NA,nrow=nrow(S),ncol=ncol(S))
  I[lower.tri(I,diag=F)]<-paste0("[",round(l[order(l$column,l$row),]$l,2),",",round(u[order(u$column,u$row),]$u,2),"]")

  II<-as.table(I, na="", quote=F)
  colnames(II)=rownames(II)=colnames(S)

  #PD interval
  PDI<-matrix(NA,nrow=nrow(S),ncol=ncol(S))
  PDI[lower.tri(PDI,diag=F)]<-paste0("[",round(rl[order(rl$column,rl$row),]$r,3),",",round(ru[order(ru$column,ru$row),]$r,3),"]")

  PDI1<-as.table(PDI, na="", quote=F)
  colnames(PDI1)=rownames(PDI1)=colnames(S)

  P1<-matrix(P, nrow=nrow(P), ncol=ncol(P))
  P1[upper.tri(P1, diag=T)]<-NA
  PP1<-as.table(round(P1,3), na="", quote=F)
  colnames(PP1)=rownames(PP1)=colnames(P)

  S2<-matrix(S, nrow=nrow(S), ncol=ncol(S))
  S2[upper.tri(S2, diag=T)]<-NA
  SS2<-as.table(round(S2,3), na="", quote=F)
  colnames(SS2)=rownames(SS2)=colnames(S)

  out<-list()

  out$full.percentile.table<-WW1
  out$simple.percentile.table<-WW
  out$residual.interval<-II
  out$posdef.interval<-PDI1
  out$predicted<-PP1
  out$observed<-SS2

  base::cat(paste("Opdyke percentiles between", lower, "and", upper, ": \n\n"))
  print(W, na.print = "")

  return(out)

}
