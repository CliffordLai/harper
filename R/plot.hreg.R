#' Plot Diagnostics for "hreg" Object
#' 
#' @description Four plots of dignostic checks are available to check the adequacy of the
#' fitted four parameter harmonic regression model.
#' 
#' @param x S3 class "hreg" object.
#' @param type Type of diagnostic plot - see details.
#' @param ... Other parameters to be passed through to predicting functions.
#' 
#' @details 
#' The diagnostic plot types are: "fit", "pg", "acf", "nplot", "cp", "rf" corresponding
#' to fitted model and data, periodogram, autocorrelation, normal probability plot, 
#' cumulative periodogram test, residual-fit spread plot.
#'
#' @return No value is returned. Plots are produced as side-effect.
#' 
#' @author A.I. McLeod and Yuanhao Lai.
#' 
#' @seealso \code{\link{hreg}}.
#' 
#' @examples
#' ans <- hreg(nottem)
#' plot(ans)
#'                                                        
#' @keywords ts
#'
plot.hreg <- function(x, type=c("fit", "pg", "acf", "nplot"), ...){
 multipanelQ <- length(type)==4
 if (multipanelQ) {
  layout(matrix(1:4, nrow=2))
 } else {
  if (length(type)==1) {
   readkey <- function() NULL
  }
  if (length(type) > 1) {
   readkey <- function()
   {
    cat ("Press [enter] to continue")
    line <- readline()
   }
  }
 }
 if ("fit" %in% type) {
  if (!multipanelQ) readkey()
  yf <- x$z - x$residual
  ylim <- range(c(x$z, yf))
  ylim <- ylim + sign(ylim)*c(-1,1)*0.045/diff(ylim)*ylim
  fopt <- x$coefficients[4]
  X <- cbind(cos(2*pi*fopt*x$t), sin(2*pi*fopt*x$t))
  dimnames(X)[[2]] <- c("A", "B")
  X <- as.data.frame.matrix(X)
  out <- lm(x$z ~ A+B, data=X)
  tf <- seq(min(x$t), max(x$t), length.out=250)
  Xnew <- cbind(cos(2*pi*fopt*tf), sin(2*pi*fopt*tf))
  dimnames(Xnew)[[2]] <- c("A", "B")
  Xnew <- as.data.frame.matrix(Xnew)
  yf <- predict(out, newdata=Xnew)
  plot(tf, yf, type="l", xlab="t", ylab="z", ylim=ylim, lwd=2, 
       col=rgb(1,0,0,0.5))
  points(x$t, x$z, pch=16)
  title(main = "Model Fit and Data")
 }
 if ("pg" %in% type) {
  if (!multipanelQ) readkey()
  xlab <- "frequency, f"
  ylab <- "I(f)"
  ttl <- "Quasi-Periodogram"
  ylab <- expression({R^{2}})
  plot(x$periodogram, lwd=2, col="blue", type="l", ylab=ylab, xlab=xlab)
  title(main = ttl)
 }
 if ("nplot" %in% type) {
  if (!multipanelQ) readkey()
  ylim <- range(c(x$residual))
  ylim <- sign(ylim)*c(-1,1)*1.06*ylim
  qqnorm(x$residual, pch=16, ylim=ylim, ylab="empirical", xlab="expected")
  qqline(x$residual)
 }
 if ("rf" %in% type) {
  if (!multipanelQ) readkey()
  yf <- x$z - x$residual
  ares <- sqrt(abs(x$residual))
  plot.default(yf, ares, xlab="fit", ylab="spread", pch=16)
  out <- loess(ares ~ yf, span=0.8, degree=1)
  xl <- seq(min(yf), max(yf), length.out=100)
  yl <- predict(out, newdata=xl)
  lines(xl, yl, lwd=2, col=rgb(1,0,0,0.6))
 } 
 if ("cp" %in% type) {
  if (!multipanelQ) readkey()
  cpgram(x$residual, ci.col="red", 
         main = "Cumulative Periodogram Test with 1% Limit")
  title(ylab="C(f)", sub="f")
 } 
 if ("acf" %in% type) {
  if (!multipanelQ) readkey()
  acf(x$residuals, main="Residual Autocorrelations with 95% Limits")
 } 
 layout(1)
}
