#' Robust Fitting of Four Parameter Harmonic Regression
#' 
#' @description  
#' Estimates mu, A, B and lambda in the harmonic regression,
#' \eqn{y(t)=mu+A*cos(2*pi*lambda*t)+B*sin(2*pi*lambda*t)+e(t)},
#' where e(t) is assumed IID mean zero and constant variance.
#' The \code{robustbase::lmrob} function is used to fit the three parameter 
#' harmonic regression of a frequency grid. The optimal frequency corresponds 
#' to the harmonic regression with the smallest deviance.
#' The accuracy is controlled K, usually K=500, is sufficient.
#' 
#' @details 
#' The three parameter harmonic regression model is fitted 
#' for each frequency in the grid which is a subset of the 
#' equi-spaced point \eqn{{j/(2K+1), j=1,\dots,K}} possibly excluding low and
#' high frequencies as specified by the lambdaRange argument. Frequencies lower
#' than the first Fourier frequency, 1/n, where n is the length time series,
#' n=length(z) are usually due to trend while high frequencies are often not of
#' interest and the numerical algorithm for fitting the 3-parameter harmonic
#' regression is sometimes inaccurate at high. 
#' The function \link[robustbase]{lmrob} uses the MM robust method
#' which is also implemented in \link[MASS]{rlm}. 
#' The algorithm implemented \link[robustbase]{lmrob} provides 
#' efficient computation. The MM robust
#' regression method is robust as well as resistant to very large outliers.
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
#' \item RSq
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
#' @references 
#' Venables, W. N. and Ripley, B. D. (2002) 
#' Modern Applied Statistics with S. Fourth edition. Springer. 
#'  
#' @seealso \link[robustbase]{lmrob}, \link[MASS]{rlm}
#' 
#' @examples 
#' z<-c(0.42, 0.89, 1.44, 1.98, 2.21, 2.04, 0.82, 0.62, 0.56, 0.8, 1.33)
#' fitHRegLmrob(z)
#' #
#' #on multicore pcs, more the package parallel may be used for the
#' #grid computation but unless n is very large this is not recommended.
#' \dontrun{ #adjust ncpu
#' system.time(ans1 <- fitHRegLmrob(z)) #2.89 sec on my computer
#' system.time(ans2 <- fitHRegLmrob(z, ncpu=8)) #2.65 sec 
#' identical(ans1, ans2)
#' }
#' @keywords ts

fitHRegLmrob <- function(z, t=1:n, K=max(500, floor(length(z)/2)), ncpu=1, 
                         lambdaRange=c(1/n, 0.45)) {
 n <- length(z)
 stopifnot(n >= 8)
 stopifnot(length(lambdaRange)==2)
 stopifnot(ceiling(ncpu)==ncpu && ncpu>0)
 K <- ifelse(K>=floor(length(z)/2), K, floor(length(z)/2)) #reset if needed
 stopifnot(length(K)==1 && ceiling(K)==K && K>= floor(length(z)/2))
 stopifnot(is.numeric(t) && length(t)==length(z))
 lambda <- (1:K)/(2*K+1)
 lambda <- lambda[lambda >= lambdaRange[1] & lambda <= lambdaRange[2]]
 nlambda <- length(lambda)
 if (ncpu==1) {
  RSq <- pg <- numeric(nlambda)
  for (i in 1:nlambda) {
    A <- cos(2*pi*lambda[i]*t)
    B <- sin(2*pi*lambda[i]*t)
    ans <- robustbase::lmrob(z ~ A+B, setting = "KS2011")
    pg[i] <- sum(coef(ans)[2:3]^2)
    RSq[i] <- summary(ans)$r.squared
   }
  } else { # start *parallel*
   ind <- parallel::splitIndices(nlambda, ncpu)
   lams <- vector("list", ncpu)
   for (i in 1:ncpu) {
    lams[[i]] <- lambda[ind[[i]]]
   }
   getR <- function(z, t, fr) {#utility function
    pg <- R2 <- numeric(length(fr))
    for (i in 1:length(fr)) {
     A <- cos(2*pi*fr[i]*t)
     B <- sin(2*pi*fr[i]*t)
     ans <- robustbase::lmrob(z ~ A+B, setting = "KS2011")
     R2[i] <- summary(ans)$r.squared
     pg[i] <- sum(coef(ans)[2:3]^2)
    }
    list(R2=R2, pg=pg)
   }#end getR
   #use clusterMap
   cl <- parallel::makeCluster(ncpu)
   parallel::clusterExport(cl, c("z", "t", "getR"), 
                           envir = environment())
   out <- parallel::clusterMap(cl, fun=function(L) getR(z,t,L), lams)
   parallel::stopCluster(cl)
   lpg <- lRSq <- vector("list", ncpu)
   for (i in 1:ncpu) {
    lpg[[i]] <- out[[i]]$pg
    lRSq[[i]] <- out[[i]]$R2
   }
   pg <- unlist(lpg)
   RSq <- unlist(lRSq)  
  }# end *parallel* Next: optimal frequency
 iOpt <- which.max(RSq)
 iOpt2 <- which.max(pg)
 RSqOpt <- RSq[iOpt]
 fOpt <- lambda[iOpt]
 #ready output
 pg <- cbind(frequency=lambda, periodogram=pg, RSq=RSq) #for output
 #refit with fOpt 
 A <- cos(2*pi*fOpt*t)
 B <- sin(2*pi*fOpt*t)
 ans <- robustbase::lmrob(z ~ A+B, setting = "KS2011")
 co <- c(coef(ans), fOpt)
 names(co) <- c("mu", "A", "B", "lambda")
 ans<-list(coefficients=co, RSq=RSqOpt, residuals=as.vector(resid(ans)), 
            periodogram=pg, z=z, t=t, K=K,  model="MM", lambdaRange=lambdaRange)
  class(ans)<-"hreg"
  ans
 }

 