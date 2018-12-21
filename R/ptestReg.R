#' Test short time series for periodicity with maximum likelihood ratio tests
#' 
#' @description  
#' This function is used to test the existence of the periodicity for 
#' a short time series (length<=100). Likelihood ratio tests under 
#' the Gaussian or the Laplace assumptions are provided with 
#' the response surface method implemented for efficiently obtaining 
#' accurate p-values under the default setting.
#' 
#' @param z A vector or a matrix containg series as columns.
#' @param t Length of a series.
#' @param method The statistical test to be used. See details for more information.
#' @param lambdalist A vector containg all plausible values of frequencies.
#' @param returnPvalue Whether the p-value of the test is returned.
#' 
#' @return Object of class "Htest" produced.
#' 
#' An object of class "Htest" is a list containing the following components:
#' 
#' \item{obsStat}{Vector containing the observed test statistics.}
#' \item{pvalue}{Vector containing the p-values of the selected tests.}
#' \item{freq}{Vector containing the estimated frequencies.}
#' 
#' @details The null hypothesis is set as no peridicities, H0: f=0. 
#' Discriptions of different test statistics (methods) are as follow:
#' 
#' \code{LS}: The -2 loglikelihood ratio test statistic based on 
#' the likelihood ratio test with normal noises, 
#' where the p-values are efficiently computed by 
#' the response surface method.
#' 
#' \code{L1}: The -2 loglikelihood ratio test statistic based on 
#' the likelihood ratio test with Laplace noises, 
#' where the p-values are efficiently computed by 
#' the response surface method.
#' 
#' @author Yuanhao Lai and A.I. McLeod
#' 
#' @references
#' 
#' Islam, M.S. (2008). Peridocity, Change Detection and Prediction in Microarrays. 
#' Ph.D. Thesis, The University of Western Ontario. 
#' 
#' Li, T. H. (2010). A nonlinear method for robust spectral analysis. 
#' Signal Processing, IEEE Transactions on, 58(5), 2466-2474.
#' 
#' MacKinnon, James (2001) : 
#' Computing numerical distribution functions 
#' in econometrics, Queen's Economics Department Working Paper, No. 1037.
#' 
#' @seealso \code{\link{ptestg}}, \code{\link{hreg}}.
#' 
#' @examples 
#' # Simulate the harmonic regression model with standard Gaussian error terms
#' set.seed(193)
#' # Non-Fourier frequency
#' z <- shreg(10, f=2.5/10, hpar=list(snr=10), model="Gaussian")
#' ptestReg(z,method = "LS") #Normal likelihood ratio test
#' ptestReg(z,method = "L1") #Laplace likelihood ratio test  
#' fitHRegLS(z)
#' fitHRegL1(z)   
#' 
#' z <- shreg(11, f=0.42, hpar=list(snr=10), model="Laplace")
#' ptestReg(z,method = "LS") #Normal likelihood ratio test
#' ptestReg(z,method = "rank")
#'                                                                                   
#' # Performe tests on the alpha factor experiment
#' data(alpha)
#' ## Eliminate genes with missing observations
#' alpha.nonNA <- alpha[complete.cases(alpha),]
#' ## Transpose the data set so that each column stands for a gene
#' alpha.nonNA <- t(alpha.nonNA)
#' result <- ptestReg(alpha.nonNA[,1:10], method = "LS") 
#' str(result)       
#'
#'
#' # The movtivating example: gene ORF06806 in Cc
#' data(Cc)
#' x <- Cc[which(rownames(Cc)=="ORF06806"),]
#' plot(1:length(x),x,type="b", main="ORF06806",
#'      xlab="time",ylab="Gene expression")
#' ptestg(x,method="Fisher") #Fail to detect the periodicity
#' ptestReg(x,method="LS") #The periodicity is significantly not zero
#' ptestReg(x,method="rank") #The periodicity is significantly not zero
#' ptestReg(x,method="L1") #The periodicity is significantly not zero
#'
#' @keywords ts
ptestReg <- function(z, t=1:n,
                     lambdalist = (ceiling(1001/length(z)):500)/1001,
                     method=c("LS", "L1","rank","MRobust"),
                     returnPvalue=TRUE){
  #Check the validity of the arguments 
  method <- match.arg(method)
  if(class(z)=="matrix"){
    nm <- dim(z)
    n <- nm[1]
    m <- nm[2]
  }else if(class(z)=="numeric"){
    n <- length(z)
    m <- 1
    z <- as.matrix(z,ncol=1,drop=FALSE)
  }else{
    stop("z must be a numerical matrix or vector")
  }
  
  
  if(n>1000){
    warning("Memory may run out because length is too large. 
             ptestg() with Fisher's g test is suggested.") 
    info <- readline("Press N and Enter for continue, otherwise stop.") 
    if(info!="N"){
      stop()
    }
  }

  obsStat <- numeric(m)
  freq <- numeric(m)
  pvalue <- NULL
  
  #Define the function for QR decompostion with the adapted format
  QRX <- function(n){
   # t <- 1:n
    #Set potential frequencies
    # K <- 500
    # aF <- 2.0*K+1
    # lambda <- (ceiling(aF/n):K)/aF
    lambda <- lambdalist
    nlambda <- length(lambda)
    #Set QR decomposition of the design matrix, names are important!
    #Store it in a large matrix
    #Append it by columns and delete 3 rows
    tQslong <- matrix(0,nrow = n-3,ncol = n*nlambda) 
    for(nl in 1:nlambda){
      X1 <- cos(2*pi*lambda[nl]*t)
      X2 <- sin(2*pi*lambda[nl]*t) 
      s1 <- 1+n*(nl-1)
      s2 <- s1+n-1
      tQslong[,s1:s2] <- t(qr.Q(qr(cbind(1,X1,X2)),complete = TRUE))[4:n,]
    }
    list(nlambda=nlambda,lambda=lambda,tQslong=tQslong)
  }
  
  if(method=="LS"){ #the likelihood ratio test with normal noises
    tableRegLs <- NULL
    
    #Compute the SSE under the Null
    SSE1 <- apply(z,2,FUN = function(x){sum((x-mean(x))^2)})
    
    #Compute the minimal SSE under the Alternative
    #The two rows store the minimal SSE and its position respectively
    nlambdaQ <- QRX(n)
    SSE2 <- minSSEC(nlambdaQ$tQslong,z) 
    optIndex <- SSE2[2,]
    SSE2 <- SSE2[1,]
    
    #Test statistic
    obsStat <- n*log(SSE1/SSE2)
    freq <- nlambdaQ$lambda[optIndex]    
    
    #Compute pvalues
    if(returnPvalue){
      data("tableRegLs", package = "harper", envir = environment())
      pvalue <- pvalueRSR(n, obsStat, tableRegLs)
    }
    
  }else if(method=="rank"){ #the likelihood ratio test with normal noises
    tableRegLsRank <- NULL
    
    #Transform the original series to ranks
    RZ <- apply(z,2,FUN = rank, ties.method = "first")
    
    #Compute the SSE under the Null
    SSE1 <- apply(RZ,2,FUN = function(x){sum((x-mean(x))^2)})
    
    #Compute the minimal SSE under the Alternative
    #The two rows store the minimal SSE and its position respectively
    nlambdaQ <- QRX(n)
    SSE2 <- minSSEC(nlambdaQ$tQslong,RZ) 
    optIndex <- SSE2[2,]
    SSE2 <- SSE2[1,]
    
    #Test statistic
    obsStat <- n*log(SSE1/SSE2)
    freq <- nlambdaQ$lambda[optIndex] 
    
    if(returnPvalue){
      data("tableRegLsRank", package = "harper", envir = environment())
      pvalue <- pvalueRSR(n, obsStat, tableRegLsRank)
    }
    
  }else if(method=="L1"){ #the likelihood ratio test with Laplace noises
    tableRegL1LaplaceEven <- NULL
    tableRegL1LaplaceOdd <- NULL
    for (i in 1:m) {
      ansHReg <- fitHRegL1(z[,i],K = 500,lambdaRange = c(1/n,0.5))
      SSE1 <- sum(abs(z[,i]-mean(z[,i]))) 
      SSE2 <- sum(abs(ansHReg$residuals)) 
      obsStat[i] <- n*log( SSE1/SSE2 )
      freq[i] <- ansHReg$coefficients[4]
    }
    
    if(returnPvalue){
      if(identical(n%%2,0)){
        data("tableRegL1LaplaceEven", package = "harper", envir = environment())
        pvalue <- pvalueRSR(n,  obsStat, tableRegL1LaplaceEven)
      }else{
        data("tableRegL1LaplaceOdd", package = "harper", envir = environment())
        pvalue <- pvalueRSR(n, obsStat, tableRegL1LaplaceOdd)
      }
    }
    
  }else if(method=="MRobust"){ 
    ans <- MStats(z) 
    obsStat <- ans[,1]
    freq <- ans[,2]
    pvalue <- NULL
  }
  
  ans <- list(obsStat=obsStat, pvalue=pvalue, freq=freq)
  class(ans)<-"Htest"
  return(ans)
}
