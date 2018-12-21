#' Print Method for "hreg" Object
#' 
#' @description Print estimates of the fitted four pararmeter harmonic regression model.
#' 
#' @param x S3 class "hreg" object.
#' @param ... Other parameters to be passed through to printing functions.
#'
#' @return A terse summary is displayed.
#' 
#' @author A.I. McLeod and Yuanhao Lai.
#' 
#' @seealso \code{\link{hreg}}.
#' 
#' @examples
#' ans <- hreg(nottem)
#' ans     
#'                                                  
#' @keywords ts
#' 
print.hreg <- function(x, ...) {
 res <- x$residual
 sdRes <- sqrt(mean(res^2))
 sdL1 <- mean(abs(res))*sqrt(pi/2)
 sdMad <- mad(res)
 snr <- sum(x$coefficients[2:3]^2)
 z <- x$z
 co <- x$coefficients
 RP <- as.vector(ABToRPhase(co[2:3]))
 co2 <- c(co, amplitude=RP[1], phase=RP[2])
 ttl1 <- "Four Parameter Harmonic Regression, "
 ttl1 <- paste0(ttl1, x$model, " Fit, K = ", x$K, ", n = ", length(x$res), 
               ", RSq = ", round(100*x$RSq,2), "%")
 cat(ttl1, fill=TRUE)
 print(co2)
 cat(paste0("s = ", round(sdRes,6), ", snr = ", round(snr/sdRes^2, 4)))
 cat(paste0(" /  mad = ", round(sdMad,6), ", snr = ", round(snr/sdMad^2,4)), 
     fill=TRUE)
 cat("summary(residual):", fill=TRUE)
 print(summary(res))
}
