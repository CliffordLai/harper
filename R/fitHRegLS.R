#' LS Fit for Four Parameter Harmonic Regression
#' 
#' @description  
#' Estimates mu, A, B and lambda in the harmonic regression,
#' \eqn{y(t)=mu+A*cos(2*pi*lambda*t)+B*sin(2*pi*lambda*t)+e(t)},
#' where e(t) is assumed NID mean zero and constant variance. 
#' The default algorithm is enumerative. 
#' 
#' @details The QR decomposition is used to efficiently compute the residual
#' sum of squares on the grid. When \code{exactQ=TRUE}, a NLS
#' algorithm is used and NLS is indicated in the output title.
#'
#' @param z series.
#' @param t Time points.
#' @param K number of subintervals.
#' @param ncpu number of compute nodes for use with parallel.
#' @param lambdaRange range of frequencies inside (0,0.5), see Details.
#' @param exactQ, default setting is FALSE.
#' It indicates whether nolinear optimizer is used to 
#' obtain the final estimate of the frequency.
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
#' @seealso 
#' \code{\link[base]{qr}}
#' 
#' @examples 
#' z<-c(0.42, 0.89, 1.44, 1.98, 2.21, 2.04, 0.82, 0.62, 0.56, 0.8, 1.33)
#' fitHRegLS(z)
#' #
#' #on multicore pcs, more the package parallel may be used for the
#' #grid computation but unless n is very large this is not recommended.
#' \dontrun{ #adjust ncpu
#' system.time(ans1 <- fitHRegLS(z)) #0.06 sec on my computer
#' system.time(ans2 <- fitHRegLS(z, ncpu=8)) #1.67 sec 
#' identical(ans1, ans2)
#' }
#' 
#' @keywords ts

fitHRegLS <- function (z, t=1:n, K=max(500, floor(length(z)/2)), ncpu=1,
              lambdaRange=c(1/n, 0.45), exactQ=FALSE) {
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
#efficient computation using QR decomposition directly
    if (ncpu==1) {
     pg <- SumSq <- numeric(nlambda)
     for(i in 1:nlambda) {
      QR <- qr(cbind(1, cos(2*pi*lambda[i]*t), sin(2*pi*lambda[i]*t)))
      SumSq[i] <- sum(qr.resid(QR, z)^2)
      pg[i] <- sum(qr.coef(QR, z)[2:3]^2)
     }
    } else { # start *parallel*
     ind <- parallel::splitIndices(nlambda, ncpu)
     lams <- vector("list", ncpu)
     for (i in 1:ncpu) {
      lams[[i]] <- lambda[ind[[i]]]
     }
     getpg <- function(z, t, fr) {#utility function
      pg <- SSE <- numeric(length(fr))
      for (i in 1:length(fr)) {
       QR <- qr(cbind(1, cos(2*pi*lambda[i]*t), sin(2*pi*lambda[i]*t)))
       SSE[i] <- sum(qr.resid(QR, z)^2)
       pg[i] <- sum(qr.coef(QR, z)[2:3]^2)
      }
      list(SSE=SSE, pg=pg)
     }#end getpg
     #use clusterMap
     cl <- parallel::makeCluster(ncpu)
     parallel::clusterExport(cl, c("z", "t", "getpg"), 
                             envir = environment())
     out <- parallel::clusterMap(cl, fun=function(L) getpg(z,t,L), lams)
     parallel::stopCluster(cl)
     lpg <- lSSE <- vector("list", ncpu)
     for (i in 1:ncpu) {
      lpg[[i]] <- out[[i]]$pg
      lSSE[[i]] <- out[[i]]$SSE
     }
     pg <- unlist(lpg)
     SumSq <- unlist(lSSE)  
    }# end *parallel*
    RSq <-1 - SumSq/sum((z-mean(z))^2)
    #Next: optimal frequency
    iOpt <- which.max(RSq)
    iOpt2 <- which.max(pg)
    RSqOpt <- RSq[iOpt]
    fOpt <- lambda[iOpt]
    #ready output
    pg <- cbind(frequency=lambda, periodogram=pg, RSq=RSq) #for output
#if exact computation, use optim() to obtain final ***fopt***
    if(exactQ) fopt <- optim(par = fOpt, fn=function(fr){sum(qr.resid(
                 qr(cbind(1, cos(2*pi*fr*t), sin(2*pi*fr*t))), z)^2)},
                 method = "L-BFGS-B",lower = 1/n, upper = 0.5)$par
#re-fit final model with frequency fopt
    ansLM <- lm(z~cbind(cos(2*pi*fOpt*t), sin(2*pi*fOpt*t)))
    co <- c(coef(ansLM), lambdaHat=fOpt)
    names(co)<-c("mu", "A", "B", "lambda")
    res <- resid(ansLM)
    model <- ifelse(exactQ, "NLS", "LS")
    ans<-list(coefficients=co, RSq=RSqOpt, residuals=res, periodogram=pg, z=z, 
              t=t, K=K, model=model, lambdaRange=lambdaRange)
    class(ans)<-"hreg"
    ans
  }