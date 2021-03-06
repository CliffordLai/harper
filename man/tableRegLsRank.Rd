% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{tableRegLsRank}
\alias{tableRegLsRank}
\title{The result for the RSR method with respect to 
the normal likelihood ratio test to the rank transformed series}
\format{A list object with the following components:

\code{model_matrix}: A matrix with the estimated coefficents 
                     for the Response Surface Regressions (one row per quantile). 

\code{nc}: a vectpr of the lengths of the simulated series.

\code{pc}: A vector containing all probabilities at which  
           the quantiles are estimated.

\code{r}: the rate in the basis function g(n)=(1/n)^r.

\code{q}: the number of parameters (except the intercept) used in 
          the response surface regression, the intercept term is removed if negative.}
\usage{
tableRegLsRank
}
\description{
This data set is used to build up the response surface regressions 
in order to obtain the p-values for the normal likelihood ratio test.
The simulation set-up of the response surface regression method is 
i.i.d. statistics of length 10^6 with 1000 replications and the
interpolated series sizes are c(8:19,seq(20,50,2),seq(55,100,5)).
}
\references{
MacKinnon, James (2001) : 
Computing numerical distribution functions 
in econometrics, Queen's Economics Department Working Paper, No. 1037.
}
\keyword{internal}
