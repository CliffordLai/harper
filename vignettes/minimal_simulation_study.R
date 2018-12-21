## ----generatorMSE--------------------------------------------------------
library(stabledist)
library(harper)

simMSE <- function(freq,SNR,dist){
  ################################
  ### Simulate required series ###
  ################################
  n <- 12     #series length
  M <- 100    #number of simulations
  t <- 1:n
  g <- function(x){x} #identify g function
  phi <- runif(1, min=0, max=2*pi)
  u <- cos(2*pi*freq*t+phi)   # Underlying true series
  
  # Simulate error term
  if(dist=="NID"){
    e <- matrix(rnorm(n*M,sd = 1),ncol = M)/sqrt(SNR)
  }else if(dist=="Stable"){
    e <- matrix(rstable(n = n*M, alpha = 1.5, beta = 0),ncol=M)/sqrt(SNR)
  }
  y <- g(u)+e   # Observed series

  ######################################
  ### Generate estimated frequencies ###
  ######################################
  # LS estimate for frequency
  ansLS <- ptestReg(y,method = "LS",returnPvalue = FALSE)
  freqLS <- ansLS$freq

  # Rank-based estimate for frequency
  lambdaRange <- seq(1/n,0.5-1/5/n,length.out = 5*n)
  ansRank <- GetFitRankLS(y,
                            t=1:n,
                            lambdaRange,
                            phiRange = seq(0,0.99,by = 0.01),
                            Exact=FALSE)
  freqRank <- ansRank[2,]
  
  #######################
  ### Compute the MSE ###
  #######################
  MSELS <- mean((freqLS-freq)^2)
  MSERank <- mean((freqRank-freq)^2)
  
  sdMSELS <- sd((freqLS-freq)^2)/sqrt(100)
  sdMSERank <- sd((freqRank-freq)^2)/sqrt(100)
    
  return(rbind(c(MSELS,sdMSELS),
               c(MSERank,sdMSERank)))
}

## ----simulationStudy-----------------------------------------------------
set.seed(193)
dist0 <- c("NID","Stable")
SNR0 <- c(1,2,3)
freq0 <- c(0.1,0.2,0.3)
ansTab <- expand.grid(Method=c("LS","Rank"),dist=dist0, SNR=SNR0, freq=freq0,stringsAsFactors = FALSE)
ansTab <- cbind(ansTab,MSE=0,sdMSE=0)

N <- nrow(ansTab)
ptm <- proc.time()
for(i in seq(1,N-1,2)){
  ansTab[i:(i+1),5:6] <- simMSE(ansTab$freq[i],
                                ansTab$SNR[i],
                                ansTab$dist[i])
}
proc.time() - ptm

## ---- fig.cap = "MSE Comparison when series length is 12", fig.width=7,fig.height=3----
library(ggplot2)
ansTab$dist <- ordered(ansTab$dist, levels=c("NID", "Stable"))
ansTab$Method <- ordered(ansTab$Method, levels=c("LS", "Rank"))

ggplot(ansTab, mapping=aes(x=SNR, y=MSE, colour=Method)) +
  geom_point(size=1) +
  geom_line(size=1) +
  facet_grid(~dist*factor(freq)) +
#  geom_errorbar(aes(SNR, ymin=MSE-1.96*sdMSE, ymax=MSE+1.96*sdMSE), width=.2,size=0.4)+
  scale_x_continuous(breaks=c(0,2))



