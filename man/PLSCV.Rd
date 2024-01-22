\name{PLSCV}
\alias{PLSCV}
\title{PLSCV}
\usage{
PLSCV(M,V,ncompMax=20,validation="repeatedcv")
}
\description{
peform PLS on hyperspec object (M), with outcome column name (V) and provide the optimum number of PC and its Rsquared
}
\arguments{
\item{M}{A Hyperspec object.}
\item{V}{A string of the outcome column name.}
\item{ncompMax}{maximum number of principle componants used}
\item{validation}{a string of the name validation method used}
}
\value{
Two numerics: 1.the optimum number of PC and 2.its Rsquared.
}
\examples{
PLSCV(M=NIR,V="GAGs",ncompMax=20,validation="repeatedcv")
}
