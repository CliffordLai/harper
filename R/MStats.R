#' Compute the test statistic based on the robust regression
#' 
#' @description  test statistic based on rank
#' 
#' @param z A vector or a matrix containg series as columns.
#' @param lambdaRange A vector containing candidate frequencies.
#'
#' @return matrix with two columns, where the first column contains test statistics
#' and the second column contains estimated frequencies.
#' 
#' @details 
#' TBA
#' 
#' @author Yuanhao Lai
#' 
#' @references  TBA
#' @keywords internal

MStats <- function(z,lambdaRange=seq(1/N, 0.5-1/N,length.out = 2*N)) {
  #Check the validity of the arguments 
  if(class(z)=="matrix"){
    N <- dim(z)[1]
    M <- dim(z)[2]
  }else if(class(z)=="numeric"){
    N <- length(z)
    M <- 1
    z <- as.matrix(z,ncol=1,drop=FALSE)
  }else{
    stop("z must be a numerical matrix or vector")
  }
  
  
  #Initialize
  t <- 1:N
  nlambdaR <- length(lambdaRange)
  typeRLik <- typeScale <- numeric(nlambdaR)
  FMfitAll <- lapply(1:nlambdaR,FUN = function(x){NULL})
  MLik <- numeric(M)
  freq <- numeric(M)
  
  ##Set the design matrix
  RRX <- array(0,
               dim = c(N,3,nlambdaR),
               dimnames=list(NULL, c("Xi1","Xi2","Xi3"),NULL))
  for(i in 1:nlambdaR){
    X1 <- cos(2*pi*lambdaRange[i]*t)
    X2 <- sin(2*pi*lambdaRange[i]*t) 
    RRX[,,i] <- cbind(1,X1,X2)  
  }
  
  for(m in 1:M){
    y <- z[,m]
    
    #---Select the optimal frequency---#
    for(i in 1:nlambdaR){
      Xi <- RRX[,,i]
      FMfit <- robustbase::lmrob(y~Xi2+Xi3, data=data.frame(Xi,y=y), 
                                 control=lmrob.control(max.it=3000,
                                                       maxit.scale =3000,
                                                       k.max=3000,
                                                       psi="lqq") )
      FMfitAll[[i]] <- FMfit
      typeScale[i] <- FMfit$scale
    }
    #-----------------------------------#
    
    FMfit <- FMfitAll[[which.min(typeScale)]]     
    s0 <- FMfit$scale
    fCtrl <- FMfit$control
    RMfit0 <- robustbase::lmrob(y~1, data=data.frame(Xi,y=y),
                                control=robustbase::lmrob.control(max.it=5000,
                                                                  maxit.scale =5000,
                                                                  k.max=5000,psi="lqq") )
    
    RMfit <- robustbase::lmrob..M..fit(x = Xi[,1,drop=FALSE], y = y,
                                       beta.initial = RMfit0$coefficients, scale = s0,
                                       control = fCtrl, method = fCtrl$method)
    
    psi <- function(u, deriv = 0){
      Mpsi(u, cc = fCtrl$tuning.psi,
           psi = fCtrl$psi, deriv)
    }
    FMres <- FMfit$resid
    RMres <- RMfit$resid 
    FM_sRho <- sum(psi(FMres/s0, deriv = -1))
    RM_sRho <- sum(psi(RMres/s0, deriv = -1))
    tauStar <- mean(psi(FMres/s0,	deriv = 1)) / mean(psi(FMres/s0)^2, deriv = 0)
    MLik[m] <- 2*tauStar*(RM_sRho - FM_sRho)
    freq[m] <- lambdaRange[which.min(typeScale)]
  }
  
  cbind(MLik=MLik, freq=freq)
}


