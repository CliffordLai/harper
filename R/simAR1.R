#' Simulate AR(1) series
#' 
#' @description  
#' An AR(1) series with mean zero and variance 1 and
#' with autocorrelation paramater phi is simulated. 
#' 
#' @param n Length of series.
#' @param phi Autocorrelation parameter.
#'
#' @return Series of length n.
#' 
#' @details The model equation is:
#' z[t] = phi*z[t-1]+a[t],
#' where z[1] is N(0,1) and a[t] are NID(0, siga),
#' \eqn{siga=\sqrt(1/(1-phi^2))}.
#' 
#' @author A.I. McLeod
#' 
#' @references McLeod, A.I., Yu, Hao and Krougly, Z. (2007),  
#' Algorithms for Linear Time 
#' Series Analysis: With R Package, Journal of Statistical Software  23, 5 1-26.
#'                                                         
#' @keywords internal


simAR1 <- function(n, phi = 0.3){
  sde <- sqrt(1-phi^2)
  e<-rnorm(n, sd=sde) #variance of time series will be 1
  u<-numeric(n)
  u[1]<-rnorm(1)
  for (i in 2:n)
    u[i] <- phi*u[i-1] + e[i]
  u
}

