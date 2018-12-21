#' Simulate harmonic regression models
#' 
#' @description  
#' Simulates a harmonic regression or a semiparametric harmonic regression.  
#' Possible types of models are normal, t(5), Laplace, cubic and AR1.
#' 
#' @param n Length of series.
#' @param mu constant term
#' @param f Frequency.
#' @param hpar structural parameters that specify A and B. See details.
#' @param model The model used for generating the error term. See details.
#' @param sig The standard error of the series.
#' @param phi Only used if AR1 error distribution is selected.
#' @param c Only used if Contaminated Normal error distribution is selected.
#' @param g A function used to transform the harmonic component. See details.
#' 
#' @return Vector of length n, simulated harmonic series.
#' 
#' @details{ Generate a harmonic series y with length n, 
#' where 
#' \deqn{y_t = mu+g( A*cos(2*pi*f*t)+B*sin(2*pi*f*t) )+sig*e_t,\ t=1,...,n,}
#' and \eqn{e_t} comes from one of the following specified distributions 
#' with mean 0 and standard error 1 that are specified 
#' in the argument \code{model}:
#' 
#' \code{Gaussian}: A standard normal distribution (i.i.d.).
#' 
#' \code{ContaminatedNormal}: Generates a random sample from a normal mixture 
#' \eqn{p*N(0,c^2)+(1-p)*N(0,1)} where \eqn{p=1/(1+c^2))}.
#' 
#' \code{t5}: A t distribution with 5 degrees of freedom 
#' (i.i.d., standardized to mean 0 and variance 1).
#' 
#' \code{Laplace}: A Laplace (double exponential) distribution
#' (i.i.d., standardized to mean 0 and variance 1).
#' 
#' \code{cubic}: A standard normal distribution for e, 
#' but the cube transform is used on the time series.
#' 
#' \code{AR1}: An AR(1) series with autocorrelation paramater phi
#' (standardized to mean 0 and variance 1).
#' 
#' The argument \code{hpar} is a list that specifies the sinusoid. If components
#' \code{A} and \code{B} are defined, these are used. Otherwise if components
#' \code{R} and \code{zeta} are defined, these specify the amplitude and phase.
#' Finally, \code{snr} may be used to generate a model with a specified 
#' signal-to-noise ratio and random phase.
#' }
#' 
#' @author A.I. McLeod and Yuanhao Lai
#' 
#' @seealso \code{\link{hreg}}, \code{\link{ptestReg}}
#' 
#' @examples 
#' #Simulate the harmonic regression model with standard Gaussian error terms
#' z <- shreg(10, f=2.5/10, hpar=list(snr=10))
#' plot(1:10,z,type="b")
#' 
#' #Simulate the AR(1) errors
#' z <- shreg(50, f=2.5/10, hpar=list(snr=2),model="AR1", phi=0.5, mu=100)
#' plot(z, type="o")
#' ans <- hreg(z)
#' plot(ans)
#' plot(ans, type="cp")
#' 
#' #Contaminated normal with large outlier
#' set.seed(77233)
#' z <- shreg(100, f=4/10, hpar=list(snr=10), model="ContaminatedNormal", c=25)
#' (out <- hreg(z))
#' plot(out)
#' (out <- hreg(z, method="MM",K=100))
#' plot(out)
#' 
#' @keywords ts

shreg <- function(n, f=0.0, 
                  model = c("Gaussian", "ContaminatedNormal","t5", "Laplace", 
                            "cubic", "AR1"), 
                  hpar=list(A=NA, B=NA, R=NA, zeta=NA, snr=NA), 
                  mu = 0, sig = 1, phi = 0.5, c=3.0,
                  g = function(x){x} ){
  model <- match.arg(model)
  if(phi>=1 | phi<=-1) 
     stop("phi should be in the range (-1,1) for stationarity!")
  if(sig<=0) 
     stop("sig should be positive!")
  t<-1:n
  x <- 0
  ABQ <- with(hpar, exists("A")&&exists("B"))
  if (ABQ) {
   x <- hpar$A*cos(2*pi*f*t)+hpar$B*sin(2*pi*f*t)
  }
  if (identical(x,0)) {
   RZetaQ <- with(hpar,exists("R") && exists("zeta"))
   if (RZetaQ) {
    x <- hpar$R*cos(2*pi*f*t+hpar$zeta)
   }
  }
  if (identical(x,0)) {
    snrQ <- with(hpar, exists("snr"))
    if(snrQ) {
     R <- sqrt(hpar$snr)*sig
     zeta <- runif(1, min=-pi/2, max=pi/2)
     x <- R*cos(2*pi*f*t+zeta)
   }
  }
  if(model=="cubic"){
    e <- rnorm(n)
    e <- sig*e
    z<-(x+e)^3
  } else {
    e<-switch(model,
              Gaussian=rnorm(n), 
              t5=rt(n,df = 5)/sqrt(5/3), 
              Laplace=simLaplace(n),
              AR1=simAR1(n, phi),
              ContaminatedNormal = rcono(n, c)
              )
    e <- sig*e
    z<- mu + g(x) + e
  }
  z
}