#' Compute the Robust g test statistic
#' 
#' @description  
#' The Robust g test statistic is computed for testing for periodicity.
#' 
#' @param y Vector containing the series.
#' @param freqSet The set to search frequencies.
#' 
#' @return Vector of length 2 containing the Robust g test statistic 
#' and the frequency for the maximum periodogram.
#' 
#' @details This is the a modified version of the function robust.g.test() from the GeneCycle package.
#' 
#' @author Yuanhao Lai
#' 
#' @references  Ahdesmaki, M., Lahdesmaki, H., Pearson, R., Huttunen, H., 
#' and Yli-Harja O.(2005). 
#'  \emph{BMC Bioinformatics} \bold{6}:117. 
#'  \url{http://www.biomedcentral.com/1471-2105/6/117}
#'                                                         
#' @keywords internal

# These programs are from the R package 'GeneCycle', but they are modified a
# little bit for simplification. Compute the Robust g test statistic
GetFitRobustG <- function(y) {
  Dpgram <- rankSpectrum(y)
  Dpgram <- Dpgram[-1,]
  maxL <- which.max(Dpgram[,2])
  gr <- Dpgram[maxL,2]/sum(Dpgram[,2])
  ans <- c(gr,Dpgram[maxL,1]) #statistic and frequency
  names(ans) <- c("gr","freq")
  ans
}


####################################################
# This function implements the robust, rank-based spectral
# estimator introduced in Pearson et al. 2003
#####################################################
#' The robust rank-based spectral
#' 
#' @description  This function implements the robust, rank-based spectral.
#' 
#' @param x Series.
#' @param freqSet The set to search frequencies.
#' 
#' @return A numerical value of the robust rank-based spectral.
#'                                                         
#' @keywords internal

rankSpectrum <- function(x) 
{
  ##############################################
  # Some adjustable parameters
  ##############################################
  # Length of the original sequence
  n <- length(x)
  
  # Let us define the maximum lag for the correlation coefficient:
  maxM <- n-2
  
  ##############################################
  # Correlation coefficient
  ##############################################
  # Reserve space
  Rsm <- matrix(NA, nrow = maxM+1, ncol = 1)
  nonMissing <- rep(TRUE,n)
  nmi <- 1:n # indices for nonmissing values
  
  # Mean removal
  x <- x- mean(x)
  
  # Run through all the lags
  for (lags in 0:maxM)
  {
    # Modified Spearman's method
    indexes <- 1:(n-lags)	# Initial indices
    ends <- n
    
    # Values in both the original and shifted vectors must be present:
    temp <- (nonMissing[1:(ends-lags)] + nonMissing[(lags+1):ends]) >= 2
    indexPresent <- which(temp)
    indexes <- indexes[indexPresent]	# The indices that are present in
    # both sequences
    Rsm[lags+1] <- ifelse(
      length(indexes)<=1 , 
      0 ,
      cor(x[indexes], x[indexes+lags], method="spearman" ) * length( x[indexes] )/n)   
    if(is.na(Rsm[lags+1])) Rsm[lags+1]<-0
  }
  
  # Zero-padding
  # Length of the zero-padded one-sided "rho"(=Rsm)
  zp <- 2*length(x)
  Rsm[(length(Rsm)+1):zp] <- 0
  fftemp <- fft(Rsm)[1:floor(zp/2)]
  freq <- (1:floor(zp/2)-1)/zp


  # The following implementation is as in (Ahdesmaki, Lahdesmaki et al., 2005)
  Ssm <- abs( 2*Re(fftemp) - Rsm[1] )

  # Return the spectral content, frequencies [0,pi)
  ans <- cbind(freq,Ssm)
  colnames(ans) <- c("frequency", "robustPeriodogram")
  return(ans)
}
