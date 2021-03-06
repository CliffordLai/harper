% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ABToRPhase.R
\name{ABToRPhase}
\alias{ABToRPhase}
\title{Convert harmonic regression parameters to A*cos+B*sin form}
\usage{
ABToRPhase(AB)
}
\arguments{
\item{AB}{Vector containing the A and B coefficients respectively.}
}
\value{
vector c(R, phi)
}
\description{
The harmonic regression is usually fitted as 
   \eqn{z(t) = mu + A*cos(2*pi*f*t) + B*sin(2*pi*f*t)}
but is more meaningfully expressed as
\eqn{z(t) = mu + R*cos(2*pi*f*t + 2*pi*phi)},
where R and \eqn{-0.5<phi=0.5} are the amplitude and phase.
This function converts the parameters to the amplitude-phase 
parameterization.
}
\examples{
#Bloomfield ozone example, p.17
harper:::ABToRPhase(c(-4.205,1.075))
                                                        
}
\author{
Yuanhao Lai and A. I. McLeod
}
\keyword{internal}
