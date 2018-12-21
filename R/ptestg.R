#' Test short time series for periodicity based on periodograms
#' 
#' @description  
#' This function is used to test the existence of the periodicity for 
#' a short time series (length<=100). 
#' Several methods based on periodograms are provided with 
#' the response surface method implemented for efficiently obtaining 
#' accurate p-values under the default setting.
#' 
#' @param z A series or a matrix containg series as columns
#' @param method The statistical test to be used. See details for more information.
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
#' \code{Fisher}: The Fisher's g test statistic. The p-value is computed 
#' directly from the exact distribution.
#' 
#' \code{robust}: The robust g test proprosed in Ahdesmaki et al. (2005), 
#' where the p-value is computed by the response surface regression method.
#' 
#' \code{FisherRSR}: Only for experimental purposes, the Fisher;s g test with 
#' p-value computed form the response surface regression method.
#' 
#' @author Yuanhao Lai and A.I. McLeod 
#' 
#' @references
#' Fisher, R.A. (1929). Tests of significance in harmonic analysis. 
#' Proc. Roy. Soc. A, 125, 54-59.
#' 
#' Ahdesmaki, M., Lahdesmaki, H., Pearson, R., Huttunen, H., and Yli-Harja O.(2005). 
#' \emph{BMC Bioinformatics} \bold{6}:117. 
#' \url{http://www.biomedcentral.com/1471-2105/6/117}.
#' 
#' MacKinnon, James (2001) : 
#' Computing numerical distribution functions 
#' in econometrics, Queen's Economics Department Working Paper, No. 1037.
#' 
#' @seealso \code{\link{ptestReg}}, \code{\link{hreg}}.
#' 
#' @examples 
#' # Simulate the harmonic regression model with standard Gaussian error terms
#' set.seed(193)
#' ## Non-Fourier frequency
#' z <- shreg(10, f=2.5/10, hpar=list(snr=10), model="Gaussian")
#' ptestg(z,method="Fisher")
#' ptestg(z,method="robust")
#' ptestg(z,method="FisherRSR")
#' 
#' # Performe tests on the alpha factor experiment
#' data(alpha)
#' ## Eliminate genes with missing observations
#' alpha.nonNA <- alpha[complete.cases(alpha),]
#' ## Using the multiple option to do the test for all the genes
#' ## Transpose the data set so that each column stands for a gene
#' alpha.nonNA <- t(alpha.nonNA)
#' result <- ptestg(alpha.nonNA[,1:10], method = "Fisher") 
#' str(result)              
#' 
#'                                                                         
#' # The movtivating example: gene ORF06806 in Cc
#' data(Cc)
#' x <- Cc[which(rownames(Cc)=="ORF06806"),]
#' plot(1:length(x),x,type="b", main="ORF06806",
#'      xlab="time",ylab="Gene expression")
#' ptestg(x,method="Fisher") #Fail to detect the periodicity
#' ptestg(x,method="robust") 
#' ptestg(x,method="FisherRSR") 
#' 
#' @keywords ts

ptestg <- function(z, 
                   method=c("Fisher", "robust", "FisherRSR")){
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

  if(n>100 & method!="Fisher"){ stop("length is larger than 100! Fisher's g test is suggested.") }
  obsStat <- numeric(m)
  freq <- numeric(m)
  
  if(method=="Fisher"){ #the Fisher's g test
    for (i in 1:m) {
      obsStatFreq <- GetFitFisherG(z[,i])
      obsStat[i] <- obsStatFreq[1]
      freq[i] <- obsStatFreq[2]
    }
    pvalue <- 1-pFisherg(g=obsStat,n=n, method = "exact")
    
  }else if(method=="robust"){ #the robust rank g test
    tableRgAn <- NULL 
    for (i in 1:m) {
      obsStatFreq <- GetFitRobustG(z[,i])
      obsStat[i] <- obsStatFreq[1]
      freq[i] <- obsStatFreq[2]
    }
    data("tableRgAn", package = "harper", envir = environment())
    pvalue <- pvalueRSR(n, obsStat, tableRgAn)
    
  }else if(method=="FisherRSR"){ #the Fisher'g based on RSR
    tablegEven <- NULL 
    tablegOdd <- NULL
    
    for (i in 1:m) {
      obsStatFreq <- GetFitFisherG(z[,i])
      obsStat[i] <- obsStatFreq[1]
      freq[i] <- obsStatFreq[2]
    }
    pvalue <- 1-pFisherg(g=obsStat,n=n, method = "RSR")
  }
  
  ans <- list(obsStat=obsStat, pvalue=pvalue,freq=freq)
  class(ans)<-"Htest"
  return(ans)
}
