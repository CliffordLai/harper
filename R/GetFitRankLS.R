#' Robust Estimator of the frequency in the semiparametric harmonic regression
#' 
#' @description  Conduct an estimation for the frequency parameter and the
#' phase parmeter in the semiparametric harmonic regression model.
#' 
#' @param z A vector or a matrix containg series as columns.
#' @param t A vector of sampled time points.
#' @param lambdaRange A vector containing candidate frequencies.
#' @param phiRange A bector containing candidate phases/2pi.
#' @param InvPI A boolean value indicate whether 
#' the inverse parabolic interpolation should be used 
#' for the final estimated frequency.
#' @param Exact Default is TRUE. which indicates that the optimization
#' is performed by enumerating all the combinations of 
#' given candidate parameters;
#' FALSE indicates that the optim() funtion is used with the 'Brent' algorithm
#' to find the optimal phi given all enumerated frequencies within lambdaRange.
#' This latter appoarch is efficient for a single series.
#' Notice that directly applying nonlinear optimizer on both lambda and phi 
#' easily gets trapped in  a local optimal.
#'
#' @return Vector of length 3 containing 
#' the maximal standardized Spearman's rank covariance, 
#' the estimated frequency and the estimated phase/2pi.
#' 
#' @details 
#' Let \eqn{z_t,\ (t=1,\cdots,n)} be a series from,
#' \deqn{z_{t}=g(\cos(2\pi \lambda t+\phi))+\epsilon_{t},\ t=1,\ldots,n,} 
#' where \eqn{0<\lambda<0.5}, \eqn{0\le\phi<2\pi}, and \eqn{\epsilon_{t}}'s are from
#' indenpendent identical distributions.
#' 
#' The estimated frequency and phase are obtained by maximizing
#' the standardized Spearman's rank covariance between the observed series and the fitted series
#' with respect to the parameters \eqn{\lambda} and \eqn{\phi}.
#' 
#' @author Yuanhao Lai
#'                                                         
#' @examples 
#' # Simulate the harmonic regression model with standard Gaussian error terms
#' set.seed(193)
#' lambdaR <- seq(0.1,0.45,length.out=50)
#' z <- shreg(15, f=2.5/15, hpar=list(snr=2,zeta=2*pi*0.3), model="Gaussian")
#' ptestReg(z,method = "LS")$freq #Normal likelihood ratio test
#' ptestg(z)$freq
#' GetFitRankLS(z,t=1:15,lambdaRange=lambdaR)
#' ## Finalize the frequency estimate with inverse parabolic interpolation
#' GetFitRankLS(z,t=1:15,lambdaRange=lambdaR,InvPI=TRUE) # Same result
#' GetFitRankLS(z,t=1:15,lambdaRange=lambdaR,InvPI=FALSE,Exact=FALSE)
#' 
#' #Simulate the extreme case with Cauchy errors
#' set.seed(193)
#' n <- 100
#' t <- 1:n
#' z <- cos(2*pi*0.21*t+2*pi*0.3)+rcauchy(n)
#' ptestReg(z,method = "LS")$freq #Normal likelihood ratio test
#' ptestg(z)$freq
#' GetFitRankLS(z,t=1:n,lambdaRange=lambdaR)
#' ## Finalize the frequency estimate with inverse parabolic interpolation
#' GetFitRankLS(z,t=1:n,lambdaRange=lambdaR,InvPI=TRUE) # Same result
#' GetFitRankLS(z,t=1:n,lambdaRange=lambdaR,InvPI=FALSE,Exact=FALSE)
#' 
#' #Apply the estimation on multiple series
#' data(alpha)
#' ## Eliminate genes with missing observations
#' alpha.nonNA <- alpha[complete.cases(alpha),]
#' ## Using the multiple option to do the test for all the genes
#' ## Transpose the data set so that each column stands for a gene
#' alpha.nonNA <- t(alpha.nonNA)
#' n <- nrow(alpha.nonNA)
#' result <- GetFitRankLS(alpha.nonNA[,1:5],t=1:n,lambdaRange=lambdaR,Exact=FALSE)
#' result     
#' 
#' 
#' 
#' @keywords ts
GetFitRankLS <- function(z,t=1:N,
                         lambdaRange,
                         phiRange = seq(0,0.99,by = 0.01), 
                         InvPI=FALSE,
                         Exact=TRUE){
  #Check the validity of the arguments 
  if(class(z)=="matrix"){
    N <- nrow(z)
    M <- ncol(z)
  }else if(class(z)=="numeric"){
    N <- length(z)
    M <- 1
    z <- as.matrix(z,ncol=1,drop=FALSE)
  }else{
    stop("z must be a numerical matrix or vector!")
  }
  
  if(length(t)!=N) stop("Invalid sampled time points!")
  if(length(lambdaRange)<3 & InvPI) 
    stop("The number of candidate lambda should be at least 3 for doing InvPI!")
  if(length(phiRange)<3 & InvPI) 
    stop("The number of candidate phi should be at least 3 for doing InvPI!")
  if(class(InvPI)!="logical")
    stop("InvPI must be either TRUE or FALSE!")
  if((!Exact) & InvPI){
    warning("InvPI is set to FALSE when Exact=FALSE")
    InvPI <- FALSE
  }

  
  #Begin estimation
  if(Exact){
    ans <- CMGetFitRankLoss(z,t,lambdaRange,phiRange,InvPI)
  }else{
    # Use nonlinear optimizer
    ans <- sapply(1:M,FUN = function(i){
      zi <- z[,i,drop=FALSE]
      RLS <- function(para){
        CMGetFitRankLoss(zi,t,lambdaRange,para,InvPI)[1]
      }
      
      fopt <- optim(par=0.5,fn=RLS,
                    method = "Brent",lower = 0,upper = 0.99,
                    control = list(maxit = 100,fnscale=-1))
      c( CMGetFitRankLoss(zi,t,lambdaRange,fopt$par,InvPI) )
    })
  }

  row.names(ans) <- c("Stats","Frequency","Phi")
  return(ans)
}
