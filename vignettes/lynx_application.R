## ----lynx, fig.show ='hold',fig.cap = "Lynx series and its log transformation", fig.width=6,fig.height=6----
#----------Read Data----------#
y <- as.vector(lynx)
ly <- log(y)
n <- length(y)
t <- 1821:1934

par(mfrow=c(2,1))
plot(t,y,type = "b", xlab = "year",
     main = "lynx")
plot(t,ly,type = "b", xlab = "year",
     main="Log lynx")
par(mfrow=c(1,1))

## ----LS------------------------------------------------------------------
#----------Fit the harmonic regression----------#
## Estimate the frequency (LS)
library(harper)
lambdaRange <- seq(1/n,0.5-1/n,length.out = 5*n)
system.time(
LSFit <- ptestReg(log(y),t=t,method = "LS",
                  lambdalist=lambdaRange,
                  returnPvalue = FALSE)
)
LSFit$freq

## Modeling given the estimated frequency
X <- data.frame(y=y,
                x1=cos(2*pi*LSFit$freq*t),
                x2=sin(2*pi*LSFit$freq*t))
fitHarG <- glm(y~x1+x2,data=X,family=gaussian())  #the Gaussian case
fitHarP <- glm(y~x1+x2,data=X,family=poisson())   #the Poisson case

## Fitted values at continuous time
newt <- seq(t[1],t[length(t)],0.1)
newX <- data.frame(x1=cos(2*pi*LSFit$freq*newt),
                   x2=sin(2*pi*LSFit$freq*newt))
predHarG <- predict(fitHarG,newX,type="response")
predHarP <- predict(fitHarP,newX,type="response")

## ----Rank----------------------------------------------------------------
#----------Fit the semiparametric harmonic regression----------#
## Estimate the frequency (Rank-based)
phiRange <- seq(0,0.99,by=0.01) 
system.time(paraRankLS <- GetFitRankLS(y,t=t,lambdaRange,phiRange)) 
paraRankLS[2,1]

## Modeling given the estimated frequency
fitSHarG <- semihregScam(y, t, lambda0=paraRankLS[2,1],
                          family = gaussian(), bs = "mpd") #the Gaussian case
fitSHarP <- semihregScam(y, t, lambda0=paraRankLS[2,1],
                          family = poisson(), bs = "mpd") #the Poisson case

## Fitted values at continuous time
newt <- seq(t[1],t[length(t)],0.1)
predSHarG <- predict(fitSHarG,newt)[,2]
predSHarP <- predict(fitSHarP,newt)[,2]


## ----visualFit, fig.cap = "Comparision of the fitted values", fig.width=6,fig.height=8----
library(tidyverse)
nt <- length(newt)
LynxVis <- data.frame(t=rep(c(t,newt),4),
                      y=c(y,predHarG,y,predHarP,
                          y,predSHarG,y,predSHarP),
                      Method=rep(c("Harmonic(G)","Harmonic(P)",
                                   "Semi(G)","Semi(P)"),each=n+nt),
                      Type=rep(c(rep("observed",n),rep("fitted",nt)),4))
LynxVis$Type <- ordered(LynxVis$Type, levels=c("observed", "fitted"))
LynxVis$Method <- ordered(LynxVis$Method, levels=c("Harmonic(G)", "Semi(G)",
                                                   "Harmonic(P)","Semi(P)"))

ggplot(LynxVis, mapping=aes(x=t, y=y,color=Type)) +
  geom_point(aes(shape=Type))+
  geom_line(aes(linetype=Type),lwd=1) +
  scale_colour_manual(values=c("black", rgb(1,0,0,0.6) )) +
  scale_linetype_manual(values=c("blank", "solid")) +
  scale_shape_manual(values=c(20,NA)) +
  facet_grid(Method~.) +
  labs(x="Year",y="Lynx") +
  theme(legend.position="top")

## ----boot, eval=FALSE----------------------------------------------------
#  library(boot)
#  y <- as.vector(lynx)
#  n <- length(y)
#  t <- 1821:1934
#  lambdaRange <- seq(1/n,0.5-1/n,length.out = 5*n)
#  phiRange <- seq(0,0.99,by=0.02)
#  
#  ## Estimation
#  freqfun1 <- function(d,inds,lambdaRange){
#    LSFit <- ptestReg(d$y[inds],t=d$t[inds],method = "LS",
#                      lambdalist=lambdaRange,
#                      returnPvalue = FALSE)
#    LSFit$freq
#  }
#  
#  freqfun2 <- function(d,inds,lambdaRange,phiRange){
#    paraRankLS <- GetFitRankLS(d$y[inds],t=d$t[inds],lambdaRange,phiRange)
#    paraRankLS[2,]
#  }
#  
#  ## Bootstrap using multicores
#  library(parallel)
#  ncore <- 4 #detectCores()
#  cl <- makeCluster(ncore)
#  clusterExport(cl,"ptestReg")
#  clusterExport(cl,"GetFitRankLS")
#  
#  ## LS method
#  tc <- proc.time()
#  set.seed(193)
#  bfreqLS <- boot(data=data.frame(y=log(y),t=t), statistic=freqfun1,
#                  R=500, sim = "ordinary",
#                  parallel="snow", ncpus=ncore,cl=cl,
#                  lambdaRange=lambdaRange)
#  
#  proc.time()-tc #
#  # user  system elapsed
#  # 0.84    0.18  138.17
#  bfreqLS
#  # Bootstrap Statistics :
#  #   original       bias    std. error
#  # t1* 0.1037369 0.0002357167 0.0004685198
#  
#  
#  # Rank-based method
#  tc <- proc.time()
#  set.seed(193)
#  bfreqRank <- boot(data=data.frame(y=y,t=t), statistic=freqfun2,
#                  R=500, sim = "ordinary",
#                  parallel="snow", ncpus=ncore,cl=cl,
#                  lambdaRange=lambdaRange)
#  
#  proc.time()-tc #
#  # user  system elapsed
#  # 6.47    0.32 1499.32
#  bfreqRank
#  # Bootstrap Statistics :
#  #   original       bias     std. error
#  # t1* 0.1037369 0.0001305769 0.0005192553
#  
#  
#  stopCluster(cl = cl) #Close the clusters

