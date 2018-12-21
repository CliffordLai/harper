#' Fit the semiparametric harmonic regression with monotone 
#' generalzied additive model
#' 
#' @description  Fit the semiparametric harmonic regression with monotone 
#' generalzied additive model.
#' It utilizes the R package \code{scam}.
#' The estimation procedure iterates between the coefficient of the index term 
#' and the smooth monotone P-spline function. 
#' 
#' @param y A vector containg series as columns.
#' @param t A vector of sampled time points.
#' @param lambda0 The frequency of the harmonic compoenent.
#' @param k the dimension of the basis used to represent the smooth term.
#'          See \code{\link[scam]{smooth.construct.mpi.smooth.spec}} in the R package \code{scam}.
#' @param m m+1 is the order of the B-spline basis.
#'          See \code{\link[scam]{smooth.construct.mpi.smooth.spec}} in the R package \code{scam}.
#' @param bs specify the usage of monotone increasing/decreasing P-splines ("mpi"/"mpd") 
#'           as the penalized smoothing basis.
#'          See \code{\link[scam]{shape.constrained.smooth.terms}} in the R package \code{scam}.   
#' @param family A family object specifying the distribution of \code{z} and link to use in fitting. 
#' See \code{\link[stats]{glm}} and \code{\link[stats]{family}} for more details.         
#'                 
#'                               
#' @return Object of class "semihregScam" produced.
#' 
#' An object of class "semihregScam" is a list 
#' containing the following components:
#' 
#' \item{harScam}{ An object of class \code{scam}. 
#' See \code{\link[scam]{scam}} in the R package \code{scam}.}
#' \item{lambda0}{The frequency of the harmonic compoenent.}
#' \item{phi0}{The estimated phase/2pi of the harmonic compoenent.}
#' \item{esty}{estimated mean for the observations.}
#' 
#' 
#' @details 
#' Given frequency \eqn{\lambda} and phase \eqn{\phi},
#' fit a monotonic smoothing term of sinuosoid with 
#' a monotone increasing/decreasing P-splines, i.e.,
#' \deqn{E[y_t] = g(a_{1} cos(2\pi\lambda t)+a_{2} sin(2\pi\lambda t) ),\ t=1,...,n,}
#' where \eqn{y_{t}} consists of independet random variables, 
#' \eqn{a_{1}^2+a_{2}^2=1},
#' \eqn{g(\cdot)} is an unknown strictly monotone function in the range [-1,1], 
#' \eqn{0<f<0.5}, \eqn{\phi\in [0,2\pi]} and \eqn{\sigma>0}.
#' This can be done by using the R package \code{scam}
#' in the R package \code{scam}.
#' Since the frequency \eqn{\lambda} is prespecified in the semiparametric harmonic regression model 
#' but the phase \eqn{\phi} is unknown, then the estimations of 
#' the function \eqn{g(\cdot)} and \eqn{\phi} is conducted in a iterative way.
#' The general optimizer \code{\link[stats]{optim}} is used to minimize the
#' UBRE of the fitted with respect to \eqn{\phi}.
#' 
#' @author Yuanhao Lai
#' 
#' @seealso \code{\link[scam]{scam}}.
#' 
#' @references 
#'  Pya, N., & Wood, S. N. (2015). 
#'  Shape constrained additive models. Statistics and Computing, 25(3), 543-559.
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
#' fit1 <- semihregScam(y,t,lambda0=paraRankLS[2],bs="mpi",family=gaussian())
#' names(fit1)
#' fit1$esty
#' 
#' # Prediction for the continuous time.
#' newt <- seq(1,n,0.05) 
#' head( predy <- predict(fit1,newt)) 
#' 
#' plot(newt,g(cos(2*pi*lambda0*newt+phi0)),
#'      type="b",pch=19,main="Semiparametric with scam",
#'      ylim=c(min(predy[,2])-1,max(predy[,2])+0.5) )
#' lines(newt,predy[,2],col="red",type="l")
#' 
#' 
#' @keywords ts
semihregScam <- function(y,t, lambda0,k=-1,m=2,bs=c("mpi","mpd"),family=gaussian()){
  #
  si <- function(a2,y,x,opt=TRUE,k=-1,m=m,bs="mpi",family=gaussian()) {
    alpha <- c(1,a2) ## constrained alpha defined using free theta
    alpha <- alpha/sqrt(sum(alpha^2))  ## Ensure ||alpha||=1 for identifiability
    sx <- x%*%alpha     ## argument of smooth
    dat <- data.frame(sx=sx)
    semiFit <- scam::scam(y~s(sx,k=k,bs=bs,m=m),family=family,data = dat) ## fit model
    
    if (opt) return(semiFit$gcv.ubre) else {
      esty <- scam::predict.scam(semiFit,dat,type="response")
      ans <- list(alpha=alpha,model=semiFit,esty=esty)
      return(ans)
    }
  } 
  
  #Build up the parametric part
  bs <- match.arg(bs)
  X <- cbind(X1=cos(2*pi*lambda0*t),
             X2=sin(2*pi*lambda0*t))
  phiRange <- seq(0,0.99,by=0.01) 
  paraRankLS <- GetFitRankLS(y,t=t,lambda0,phiRange)
  a2 <- -tan(paraRankLS[3,1]*2*pi)
  
  #Start the iterative estimation
  f1 <- optim(a2,si,y=y,x=X,k=k,m=m,bs=bs,family = family,method = "BFGS")
  a2.est <-f1$par 
  
  msim <- si(a2 = a2.est,y = y,x = X,opt=FALSE,
             k=k,m=m,bs=bs,family = family)
  
  phi0 <- atan2(-msim$alpha[2],msim$alpha[1])/2/pi
  
  esty <- msim$esty
  
  ans <- list(harScam=msim$model,
              lambda0=lambda0,
              phi0=ifelse(phi0>=0,phi0,1+phi0),
              esty=esty)
  class(ans) <- "semihregScam"
  ans
}



#' @title  Prediction of the semiparametric harmonic regression with scam
#' 
#' @description  Return predicted values of the
#' semiparametric harmonic regression with scam,
#' The prediction procedure utilizes the function \code{predict.scam}
#' in the package \code{scam}.
#' 
#' @param object An object of class "semihregScam".
#' @param newt A vector of sampled time points for predicted values.
#' @param ... Other parameters to be passed through to predicting functions.
#' 
#' @return A matrix of two column, where
#' the first column contains the sample time points and
#' the second column contains the corresponding predicted values.
#' 
#' @author Yuanhao Lai
#' 
#' @seealso \code{\link{semihregScam}}, \code{\link[scam]{scam}}.
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
#' fit1 <- semihregScam(y,t,lambda0=lambda0,bs="mpi",family=gaussian())
#' names(fit1)
#' fit1$esty
#' 
#' # Prediction for the continuous time.
#' newt <- seq(1,n,0.05) 
#' head( predy <- predict(fit1,newt)) 
#' 
#' plot(newt,g(cos(2*pi*lambda0*newt+phi0)),
#'      type="b",pch=19,main="Semiparametric with scam",
#'      xlab="t", ylab="g(u)",
#'      ylim=c(min(predy[,2])-1,max(predy[,2])+0.5) )
#' lines(newt,predy[,2],col="red",type="b")
#' 
predict.semihregScam <- function(object,newt,...){
  lambda0 <- object$lambda0
  phi0 <- object$phi0
  newx <- cos(2*pi*lambda0*newt+2*pi*phi0)
  newx <- data.frame(sx=newx)
  pv <- scam::predict.scam(object$harScam,newx,type="response")
  newx <- cbind(t=newt, yhat = pv)
  newx
}
