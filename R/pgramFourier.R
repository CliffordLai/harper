#' Periodogram computation, Fourier method
#' 
#' @description  
#' The periodogram is computed at either Fourier or non-Fourier frequencies
#' 
#' @param z time series vector of length n, say.
#' @param numFreq use "default" for usual Fourier frequencies, 
#' {1/n, ..., floor(n/2)/n}. Set fr = N, to evaluate the periodogram at the 
#' Fourier frequencies corresponding to a time series of length N. The 
#' frequencies are in cycles per unit time.
#' @param demean whether the sample mean is subtracted from the series.
#' 
#' @return Periodogram
#' 
#' @details Uses FFT.
#' So if the length of z is a highly composite number, the computation
#' is very efficient. Otherwise the usual DFT is used.
#' 
#' @author A.I. McLeod and Yuanhao Lai
#' 
#' @examples 
#' z<-sunspot.year
#' n<-length(z) 
#' I<-pgramFourier(z)
#' f<-I[,1]
#' I <- I[,2]
#' plot(f, I, xlab="f", ylab="f", type="l") 
#' title(main="Periodogram for Annual Sunpots, 1700-1988") 
#' #
#' z<-c(0.42, 0.89, 1.44, 1.98, 2.21, 2.04, 0.82, 0.62, 0.56, 0.8, 1.33)
#' pgramFourier(z)
#' ans <- pgramFourier(z, numFreq=101)
#' plot(ans[,1], ans[,2], type="l", xlab="frequency", ylab=
#'  "periodogram")
#'   
#'                                                         
#' @keywords ts

pgramFourier <-
  function(z, numFreq ="default", demean=TRUE)
  {
    n <- length(z)
    if (identical(numFreq, "default")) {
     ans <- cbind((1:floor(n/2))/n, (Mod(fft(z))^2/n)[2:(n %/% 2 + 1)])
     } else { 
     stopifnot(numFreq > n)
     madj <- numFreq-n
     if (demean) z <- z-mean(z)
     x <- c(z, rep(0, madj))
     nf <- numFreq
     ans <- cbind((1:floor(nf/2))/nf, (Mod(fft(x))^2 /n)[2:(nf%/%2+1)])
     } 
     colnames(ans) <- c("frequency", "periodogram")
     ans
  }

