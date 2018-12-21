#' P-value computation for the RSR model
#' 
#' @description  Compute the estimated p-value of the observed statistics from a given RSR model
#' 
#' @param n Length of series.
#' @param Statistic A list of observed statistics.
#' @param tableRSR The RSR fitted models for the required statistic.
#' 
#' @return A list of p-values.
#' 
#' @author Yuanhao Lai
#'                                                         
#' @keywords internal
pvalueRSR <- function(n,Statistic,tableRSR){ 
  pc <- tableRSR$pc
  
  qa_est <- predictRSR(n,tableRSR)
  lx <- qa_est
  ly <- qnorm(pc)
  
  #Select an optimal smoothing parameter using the AIC criterion
  K <- length(pc)
  spanList <- (5:15)/K
  R <- length(spanList)
  ICLoess <- numeric(R)
  
  for(i in 1:R){
    CDF <- loess(ly~lx,degree = 2, span = spanList[i],
                 control = loess.control(surface = "direct") )
    nLoess <- CDF$n
    sigma2<- sum(CDF$residuals^2)/nLoess
    ICLoess[i] <- nLoess*log(sigma2)+2*CDF$enp #AIC
    #ICLoess[i] <- nLoess*log(sigma2)+log(nLoess)*CDF$enp #BIC
  }
  
  #Fit loess
  CDF <- loess(ly~lx,degree = 2, span = spanList[which.min(ICLoess)],
               control = loess.control(surface = "direct") )
  pvalue <- 1-pnorm(predict(CDF,data.frame(lx=Statistic)))

  pvalue
}

#' Quantile computation for the RSR model
#' 
#' @description  Compute the estimated quantiles of the required statistics from a given RSR model
#' 
#' @param n Length of series.
#' @param tableRSR The RSR fitted models for the required statistic.
#' 
#' @return a vector of all the estimated quantiles at pc for 
#'         specidifed sample sizes n from a given RSR model.
#' 
#' @author Yuanhao Lai
#'                                                         
#' @keywords internal
predictRSR <- function(n,tableRSR){
  model_matrix <- tableRSR$model_matrix
  r <- tableRSR$r
  q <- tableRSR$q
  
  if(q>0){
    New <- rep(1, length(n))
  }else{
    New <- NULL
  }
  
  I1 <- 1/n^r
  for(i in 1:abs(q)){
    New <- rbind( New, I1^(i*0.5+0.5) )
  }
  
  return(model_matrix%*%New)
}