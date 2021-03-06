#' Four Parameter Harmonic Regression
#' 
#' @description  
#' Estimates mu, A, B and lambda in the harmonic regression,
#' \eqn{y(t)=mu+A*cos(2*pi*lambda*t)+B*sin(2*pi*lambda*t)+e(t)},
#' where e(t) is assumed IID mean zero and constant variance. 
#' Four estimation methods are available least-squares (LS),
#' least absolute deviation (L1) and a robust regression methods (MM).
#'
#' @param z series.
#' @param K number of subintervals.
#' @param method one of "LS", "L1", "MM"/
#' @param ncpu default 1.
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
#' @details 
#' The recommended setting of ncpu is 1 - see note. 
#' The default setting for K is the 500 or the number of Fourier frequencies -
#' whichever is larger. 
#' LS regression uses the QR algorithm directly to efficiently compute the
#' residual sum of squares. 
#' L1-regression is fitted using the function \code{\link{rq.fit.br}}
#' in the package \code{quantreg}. 
#' Robust regression uses the  \code{robustbase}
#' package function \code{\link{lmrob}}.
#' The parallel package may be used to improve the efficiency for long time 
#' series but when the series length is less than 100 using ncpu=1 with K=500
#' is much faster due to the overhead - despite using \code{\link{splitIndices}}
#' to balance the load on each compute note.
#' 
#' @note 
#' The \code{parallel} package may be used to speed up the computation on the
#' frequency grid but on Windows machines the overhead for computation is so
#' large that this is not usually effective unless the series length is quite
#' long. See \code{nottem} in the example section. It was surprising that the
#' overhead was so large despite using \code{\link{splitIndices}} to balance
#' the load on each cpu.
#' 
#' @author Yuanhao Lai and A. I. McLeod
#' 
#' @seealso 
#' \code{\link{qr}} 
#' 
#' @examples 
#' z<-c(0.42, 0.89, 1.44, 1.98, 2.21, 2.04, 0.82, 0.62, 0.56, 0.8, 1.33)
#' hreg(z)
#' #
#' #on multicore pcs, more the package parallel may be used for the
#' #grid computation but unless n is very large this is not recommended.
#' \dontrun{ #adjust ncpu
#' system.time(ans1 <- hreg(z)) #0.06 sec on my computer
#' system.time(ans2 <- hreg(z, ncpu=8)) #1.67 sec 
#' system.time(ans3 <- hreg(z, method="L1")) #0.06 sec on my computer
#' system.time(ans4 <- hreg(z, method="L1", ncpu=8)) #4.11 sec
#' system.time(ans5 <- hreg(z, method="MM")) #2.89 sec on my computer
#' system.time(ans6 <- hreg(z, method="MM", ncpu=8)) #2.65 sec 
#' }
#' 
#' @keywords ts

hreg <- function(z, method=c("LS", "L1", "MM"), K=500, ncpu=1, exactQ=FALSE) {
 MTHD <- match.arg(method)
 switch(MTHD, 
        LS = fitHRegLS(z, K=K, ncpu=ncpu, exactQ=exactQ),
        L1 = fitHRegL1(z, K=K, ncpu=ncpu),
        MM = fitHRegLmrob(z, K=K, ncpu=ncpu)
 )
}
