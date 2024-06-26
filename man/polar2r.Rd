% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polar2r.R
\name{polar2r}
\alias{polar2r}
\title{Convert polar angles to a correlation matrix}
\usage{
polar2r(polar = NULL)
}
\arguments{
\item{polar}{A lower triangle matrix of polar angles in radians where each element is between 0 and pi.}
}
\value{
Correlation matrix
}
\description{
This function transforms a lower triangle matrix of polar angles to a symmetric correlation matrix.
Input is a lower triangle matrix of polar angles.
}
\examples{
#assign dataset to “dat”
dat<-lavaan::HolzingerSwineford1939
#assign correlation matrix to “A”
R<-cor(dat[,7:12])
#create lower triangle matrix of polar angles
P<-r2polar(cor=R, type="cor")
#use polar2r to transform the polar angles back to a correlation matrix
polar2r(polar=P)
}
\author{
Daniel McNeish

Maintainer: Daniel McNeish \href{mailto:dmcneish@asu.edu}{dmcneish@asu.edu}
}
