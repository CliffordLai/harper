#' Fit the semiparametric harmonic regression with monotone index model
#' 
#' @description  Fit the semiparametric harmonic regression with monotone 
#' generalzied additive index model by maximizing non-parametric likelihood.
#' It utilizes the function \code{\link[scar]{scair}} in the R package \code{scar}.
#' 
#' @param y A vector containg series as columns.
#' @param t A vector of sampled time points.
#' @param lambda0 The frequency of the harmonic compoenent.
#' @param shape The shape constraint for the function g. 
#' See \code{\link[scar]{scair}} in the R package \code{scar}.
#' @param family a description of the error distribution and link function to be used in the model.
#' Default is gaussian(). See \code{\link[scar]{scair}}.
#' @param ... Additional parameters passed to \code{\link[scar]{scair}}.
#' 
#' @return Object of class "semihregScair" produced.
#' 
#' An object of class "semihregScair" is a list 
#' containing the following components:
#' 
#' \item{harScair}{ An object of class \code{scair}. 
#' See \code{\link[scar]{scair}} in the R package \code{scar}.}
#' \item{lambda0}{The frequency of the harmonic compoenent.}
#' \item{phi0}{The estimated phase/2pi of the harmonic compoenent.}
#' \item{esty}{estimated mean for the observations.}
#' 
#' @details 
#' Given frequency \eqn{\lambda},
#' fit a monotonic single index model with the likelihood method, i.e.,
#' \deqn{E[y_t] = \mu+g( A\cos(2\pi\lambda t)+B\sin(2\pi \lambda t) ),\ t=1,...,n,}
#' where \eqn{y_{t}} consists of independet random variables, 
#' \eqn{a_{1}^2+a_{2}^2=1},
#' \eqn{g(\cdot)} is an unknown strictly monotone function in the range [0,1], 
#' \eqn{0<f<0.5}, \eqn{A^2 + B^2 =1} and \eqn{\sigma>0}.
#' 
#' This is done by using the R package \code{scar}.
#' For more details, see \code{\link{scair}} 
#' in the R package \code{scar}.
#' 
#' @author Yuanhao Lai
#' 
#' @seealso \code{\link[scar]{scair}}.
#' 
#' @references 
#' Chen, Y. and Samworth, R. J. (2014). Generalised additive and index models with shape constraints. arXiv:1404.2957.
#'                                                         
#' @examples
#' # Simulate the series
#' set.seed(193)
#' g <- function(x){-2*exp(-x^3)}
#' n <- 50
#' t <- 1:n
#' lambda0 <- 0.123
#' phi0 <- 0.2*2*pi
#' u <- cos(2*pi*lambda0*t+phi0)
#' e <- rnorm(n)
#' y <- g(u)+e
#' 
#' # Plot ther series
#' par(mfrow=c(1,2))
#' plot(t,g(u),type="b")
#' plot(t,y,type="b")
#' par(mfrow=c(1,1))
#' 
#' # Pre-estimate the frequency and the phase
#' lambdaRange <- seq(1/n,0.5-1/n,length.out = 5*n)
#' paraRankLS <- GetFitRankLS(y,t=1:n,lambdaRange)
#' paraRankLS    
#' 
#' # Fit the semi-harmonic regression
#' fit1 <- semihregScair(y,t,lambda0=paraRankLS[2],shape="in",iter=100)
#' names(fit1)
#' fit1$esty
#' 
#' # Prediction for the continuous time.
#' newt <- seq(1,n,0.05) 
#' head( predy <- predict(fit1,newt)) 
#' 
#' plot(newt,g(cos(2*pi*lambda0*newt+phi0)),
#'      type="b",pch=19,main="Semiparametric with scair",
#'      xlab="t", ylab="g(u)",
#'      ylim=c(min(predy[,2])-1,max(predy[,2])+0.5) )
#' lines(newt,predy[,2],col="red",type="b")
#' 
#' 
#' @keywords ts
semihregScair <- function(y,t,lambda0,
                          shape=c("l","in","de",
                                  "cvxin","cvxde",
                                  "ccvin","ccvde"),
                          family=gaussian(),
                          ...){
  shape <- match.arg(shape)
  
  X <- cbind(X1=cos(2*pi*lambda0*t),
             X2=sin(2*pi*lambda0*t))
  harScair <- scair(x=X,y,shape="in", family=family,...)
  phi0 <- atan2(-harScair$index[2],harScair$index[1])/2/pi
  
  esty <- predict(harScair,X)
  
  ans <- list(harScair=harScair,
              lambda0=lambda0,
              phi0=ifelse(phi0>=0,phi0,1+phi0),
              esty=esty)
  class(ans) <- "semihregScair"
  ans
}

#' @title  Prediction of the semiparametric harmonic regression with monotone index model
#' 
#' @description  Return predicted values of the
#' semiparametric harmonic regression with GAM.
#' The prediction procedure utilizes the function \code{predict.scair}
#' in the package \code{scar}.
#' 
#' @param object An object of class "semihregScair".
#' @param newt A vector of sampled time points for predicted values.
#' @param ... Other parameters to be passed through to predicting functions.
#' 
#' @return A matrix of two column, where
#' the first column contains the sample time points and
#' the second column contains the corresponding predicted values.
#' 
#' @author Yuanhao Lai
#' 
#' @seealso \code{\link{semihregScair}}, \code{\link[scar]{predict.scair}}.
#' 
#' @examples
#' # Simulate the series
#' set.seed(193)
#' g <- function(x){-2*exp(-x^3)}
#' n <- 50
#' t <- 1:n
#' lambda0 <- 0.123
#' phi0 <- 0.2*2*pi
#' u <- cos(2*pi*lambda0*t+phi0)
#' e <- rnorm(n)
#' y <- g(u)+e
#' 
#' # Fit the semi-harmonic regression
#' fit1 <- semihregScair(y,t,lambda0=lambda0,shape="in",iter=100)
#' names(fit1)
#' fit1$esty
#' 
#' # Prediction for the continuous time.
#' newt <- seq(1,n,0.05) 
#' head( predy <- predict(fit1,newt)) 
#' 
#' plot(newt,g(cos(2*pi*lambda0*newt+phi0)),
#'      type="b",pch=19,main="Semiparametric with scair",
#'      xlab="t", ylab="g(u)",
#'      ylim=c(min(predy[,2])-1,max(predy[,2])+0.5) )
#' lines(newt,predy[,2],col="red",type="b")
#' 
predict.semihregScair <- function(object,newt,...){
  predX <- cbind(X1=cos(2*pi*object$lambda0*newt),
                 X2=sin(2*pi*object$lambda0*newt))
  
  yhatScair <- predict(object$harScair,predX,type="response")
  newx <- cbind(t=newt, yhat = yhatScair)
  newx
}
