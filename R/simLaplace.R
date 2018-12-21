#' Simulate series with the Laplace distribution
#' 
#' @description  
#' A Laplace (double exponential) series with mean zero and variance 1 is simulated. 
#' 
#' @param n Length of series.
#'
#' @return vector of length n, simulated harmonic series.
#' 
#' @details The model equation is: 
#' \deqn{e=(-1)^(u1\ge0.5) * ln(u2) / sqrt(2)}
#' where u1 and u2 are independent uniform (0,1) random variables.
#' 
#' @author Yuanhao Lai
#'                                                         
#' @keywords internal


simLaplace <- function(n){
  u1 <- runif(n)
  u2 <- runif(n)
  e <- ifelse(u1>=0.5, 1/sqrt(2), -1/sqrt(2)) * log(u2) #variance of time series will be 1
  e
}