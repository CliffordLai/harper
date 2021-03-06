% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ptestReg.R
\name{ptestReg}
\alias{ptestReg}
\title{Test short time series for periodicity with maximum likelihood ratio tests}
\usage{
ptestReg(z, t = 1:n, lambdalist = (ceiling(1001/length(z)):500)/1001,
  method = c("LS", "L1", "rank", "MRobust"), returnPvalue = TRUE)
}
\arguments{
\item{z}{A vector or a matrix containg series as columns.}

\item{t}{Length of a series.}

\item{lambdalist}{A vector containg all plausible values of frequencies.}

\item{method}{The statistical test to be used. See details for more information.}

\item{returnPvalue}{Whether the p-value of the test is returned.}
}
\value{
Object of class "Htest" produced.

An object of class "Htest" is a list containing the following components:

\item{obsStat}{Vector containing the observed test statistics.}
\item{pvalue}{Vector containing the p-values of the selected tests.}
\item{freq}{Vector containing the estimated frequencies.}
}
\description{
This function is used to test the existence of the periodicity for 
a short time series (length<=100). Likelihood ratio tests under 
the Gaussian or the Laplace assumptions are provided with 
the response surface method implemented for efficiently obtaining 
accurate p-values under the default setting.
}
\details{
The null hypothesis is set as no peridicities, H0: f=0. 
Discriptions of different test statistics (methods) are as follow:

\code{LS}: The -2 loglikelihood ratio test statistic based on 
the likelihood ratio test with normal noises, 
where the p-values are efficiently computed by 
the response surface method.

\code{L1}: The -2 loglikelihood ratio test statistic based on 
the likelihood ratio test with Laplace noises, 
where the p-values are efficiently computed by 
the response surface method.
}
\examples{
# Simulate the harmonic regression model with standard Gaussian error terms
set.seed(193)
# Non-Fourier frequency
z <- shreg(10, f=2.5/10, hpar=list(snr=10), model="Gaussian")
ptestReg(z,method = "LS") #Normal likelihood ratio test
ptestReg(z,method = "L1") #Laplace likelihood ratio test  
fitHRegLS(z)
fitHRegL1(z)   

z <- shreg(11, f=0.42, hpar=list(snr=10), model="Laplace")
ptestReg(z,method = "LS") #Normal likelihood ratio test
ptestReg(z,method = "rank")
                                                                                  
# Performe tests on the alpha factor experiment
data(alpha)
## Eliminate genes with missing observations
alpha.nonNA <- alpha[complete.cases(alpha),]
## Transpose the data set so that each column stands for a gene
alpha.nonNA <- t(alpha.nonNA)
result <- ptestReg(alpha.nonNA[,1:10], method = "LS") 
str(result)       


# The movtivating example: gene ORF06806 in Cc
data(Cc)
x <- Cc[which(rownames(Cc)=="ORF06806"),]
plot(1:length(x),x,type="b", main="ORF06806",
     xlab="time",ylab="Gene expression")
ptestg(x,method="Fisher") #Fail to detect the periodicity
ptestReg(x,method="LS") #The periodicity is significantly not zero
ptestReg(x,method="rank") #The periodicity is significantly not zero
ptestReg(x,method="L1") #The periodicity is significantly not zero

}
\references{
Islam, M.S. (2008). Peridocity, Change Detection and Prediction in Microarrays. 
Ph.D. Thesis, The University of Western Ontario. 

Li, T. H. (2010). A nonlinear method for robust spectral analysis. 
Signal Processing, IEEE Transactions on, 58(5), 2466-2474.

MacKinnon, James (2001) : 
Computing numerical distribution functions 
in econometrics, Queen's Economics Department Working Paper, No. 1037.
}
\seealso{
\code{\link{ptestg}}, \code{\link{hreg}}.
}
\author{
Yuanhao Lai and A.I. McLeod
}
\keyword{ts}
