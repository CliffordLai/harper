#' L1 Fit of Four Parameter Harmonic Regression
#' 
#' @description  
#' Estimates mu, A, B and lambda in the harmonic regression,
#' \eqn{y(t)=mu+A*cos(2*pi*lambda*t)+B*sin(2*pi*lambda*t)+e(t)},
#' where e(t) is assumed IID mean zero and constant variance.
#' The absolute sum of errors is computed for each lambda on a grid of values
#' and lambdaHat corresponds to the smallest. An efficient linear programming
#' method is used to fit the model to equi-spaced frequencies in (0,0.5).
#' The fitting algorithm minimizes the sum of the absolute residuals.
#' The frequency is estimated by the model that has the best fit on the grid.
#' The accuracy is controlled K, usually K=500, is sufficient.
#' 
#' @details 
#' The three parameter harmonic regression model is fitted by least absolute
#' deviation (LAD) for each frequency in the grid which is a subset of the 
#' equi-spaced point \eqn{{j/(2K+1), j=1,\dots,K}} possibly excluding low and
#' high frequencies as specified by the lambdaRange argument. Frequencies lower
#' than the first Fourier frequency, 1/n, where n is the length time series,
#' n=length(z) are usually due to trend while high frequencies are often not of
#' interest and the numerical algorithm for fitting the 3-parameter harmonic
#' regression is sometimes inaccurate at high. 
#' Note that LAD is MLE when the errors are IID double exponential.
#' The LAD fitting uses the function \code{rq.fit.br} from the 
#' package \code{quantreg}.
#' 
#' @param z series.
#' @param t Time points.
#' @param K number of subintervals.
#' @param ncpu number of compute nodes for use with parallel.
#' @param lambdaRange range of frequencies inside (0,0.5), see Details.
#' 
#' @return list with the following components
#' \itemize{
#' \item coefficients const, A, B, lambda
#' \item residuals
#' \item periodogram
#' \item z
#' \item t
#' \item K
#' \item model
#' \item lambdaRange
#' }
#' 
#' @author Yuanhao Lai and A.I. McLeod
#' 
#' @examples 
#' z<-c(0.42, 0.89, 1.44, 1.98, 2.21, 2.04, 0.82, 0.62, 0.56, 0.8, 1.33)
#' fitHRegL1(z)
#' #
#' #on multicore pcs, more the package parallel may be used for the
#' #grid computation but unless n is very large this is not recommended.
#' \dontrun{ #adjust ncpu
#' system.time(ans1 <- fitHRegL1(z)) #0.06 sec on my computer
#' system.time(ans2 <- fitHRegL1(z, ncpu=8)) #4.11 sec
#' identical(ans1, ans2)
#' }
#' @keywords ts

fitHRegL1 <- function(z, t=1:n, K=500, ncpu=1, lambdaRange=c(1/n, 0.45)) {
 n <- length(z)
 stopifnot(n >= 8)
 stopifnot(length(lambdaRange)==2)
 stopifnot(ceiling(ncpu)==ncpu && ncpu>0)
 K <- max(K, floor(length(z)/2)) #reset if needed
 stopifnot(length(K)==1 && ceiling(K)==K && K>= floor(length(z)/2))
 stopifnot(is.numeric(t) && length(t)==length(z))
 lambda <- (1:K)/(2*K+1)
 lambda <- lambda[lambda >= lambdaRange[1] & lambda <= lambdaRange[2]]
 nlambda <- length(lambda)
 SNull <- sum(abs(z-median(z))) 
 SErr <- pg <- numeric(nlambda)
 if (ncpu==1) {
  for (i in 1:nlambda) {
   ans <- quantreg::rq.fit.br(
            x=cbind(rep(1,n), cos(2*pi*lambda[i]*t), sin(2*pi*lambda[i]*t)), 
            y=z, tau=0.5)
   pg[i] <- sum(coef(ans)[2:3]^2)
   SErr[i] <- sum(abs(ans$residuals))
  }
 } else { # start *parallel*
  ind <- parallel::splitIndices(nlambda, ncpu)
  lams <- vector("list", ncpu)
  for (i in 1:ncpu) {
   lams[[i]] <- lambda[ind[[i]]]
  }
  getpg <- function(z, t, fr) {#utility function
   SError <- pg <- numeric(length(fr))
   for (i in 1:length(fr)) {
    ans <-  quantreg::rq.fit.br(
     x=cbind(rep(1,length(z)), cos(2*pi*fr[i]*t), sin(2*pi*fr[i]*t)), 
     y=z, tau=0.5)
    pg[i] <- sum(coef(ans)[2:3]^2)
    SError[i] <- sum(abs(ans$residuals))
   }
   list(SError=SError, pg=pg)
  }
  #use clusterMap
  cl <- parallel::makeCluster(ncpu)
  parallel::clusterExport(cl, c("z", "t", "getpg"), 
                          envir = environment())
  out <- parallel::clusterMap(cl, fun=function(L) getpg(z,t,L), lams)
  parallel::stopCluster(cl)
  lpg <- lSErr <- vector("list", ncpu)
  for (i in 1:ncpu) {
   lpg[[i]] <- out[[i]]$pg
   lSErr[[i]] <- out[[i]]$SError
  }
  pg <- unlist(lpg)
  SErr <- unlist(lSErr)  
 }# end *parallel* Next: optimal frequency
 RSq <- 1-SErr/SNull
 iOpt <- which.max(RSq)
 iOpt2 <- which.max(pg)
 RSqOpt <- RSq[iOpt]
 fOpt <- lambda[iOpt]
#ready output
 pg <- cbind(frequency=lambda, periodogram=pg, RSq=RSq) #for output
#refit with fOpt
 ans <- quantreg::rq.fit.br(
          x=cbind(rep(1,n), cos(2*pi*fOpt*t), sin(2*pi*fOpt*t)), y=z, tau=0.5)
 co <- c(coef(ans), fOpt)
 names(co) <- c("mu", "A", "B", "lambda")
 ans<-list(coefficients=co, RSq=RSqOpt, residuals=as.vector(resid(ans)), 
           periodogram=pg, z=z, t=t, K=K,  model="L1", lambdaRange=lambdaRange)
 class(ans)<-"hreg"
 ans
}
