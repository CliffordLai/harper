## ----gene, fig.show ='hold',fig.cap = "The observed gene expressions", fig.width=7,fig.height=3----
#----------Read Data----------#
library(harper)
data(alpha)
str(alpha)

y1 <- alpha[row.names(alpha)=="YCL055W",]
y2 <- alpha[row.names(alpha)=="YMR198W",]
n <- length(y1)
t <- 1:n

par(mfrow=c(1,2))
plot(t,y1,type="p",pch=19,
     xlab = "time", ylab = "gene expression",
     main="YCL055W")
plot(t,y1,type="p",pch=19,
     xlab = "time", ylab = "gene expression",
     main="YMR198W")
par(mfrow=c(1,1))

## ----dist----------------------------------------------------------------
n <- 18
M <- 500
lambdaRange <- (ceiling(1001/n):500)/1001
Z <- matrix(rnorm(n*M),ncol=M)
#Takes 168s
set.seed(193)
system.time({
  distLS  <- ptestReg(Z,t=t,method = "LS",
                  lambdalist=lambdaRange,
                  returnPvalue = FALSE)$obsStat
})

system.time({
  distRankLS  <- GetFitRankLS(Z,t=1:n,
                              lambdaRange,
                              phiRange = seq(0,0.99,by = 0.01),
                              InvPI=FALSE,
                              Exact=TRUE)[1,]
})

## ----LS------------------------------------------------------------------
#----------Fit the harmonic regression----------#
## Estimate the frequency (LS)
lambdaRange <- (ceiling(1001/n):500)/1001
system.time(
LSFit1 <- ptestReg(y1,t=t,method = "LS",
                  lambdalist=lambdaRange,
                  returnPvalue = FALSE)
)
system.time(
LSFit2 <- ptestReg(y2,t=t,method = "LS",
                  lambdalist=lambdaRange,
                  returnPvalue = FALSE)
)
LSFit1$freq #gene YCL055W
LSFit2$freq #gene YMR198W

## Modeling given the estimated frequency
X1 <- data.frame(y1=y1,
                x1=cos(2*pi*LSFit1$freq*t),
                x2=sin(2*pi*LSFit1$freq*t))
X2 <- data.frame(y2=y2,
                x1=cos(2*pi*LSFit2$freq*t),
                x2=sin(2*pi*LSFit2$freq*t))
fitHar1 <- glm(y1~x1+x2,data=X1,family=gaussian())  
fitHar2 <- glm(y2~x1+x2,data=X2,family=gaussian()) 

## Fitted values at continuous time
newt <- seq(t[1],t[length(t)],0.1)
newX1 <- data.frame(x1=cos(2*pi*LSFit1$freq*newt),
                   x2=sin(2*pi*LSFit1$freq*newt))
newX2 <- data.frame(x1=cos(2*pi*LSFit2$freq*newt),
                   x2=sin(2*pi*LSFit2$freq*newt))
predHar1 <- predict(fitHar1,newX1,type="response")
predHar2 <- predict(fitHar2,newX2,type="response")

## P-values of the likelihood ratio test
pvalueLS1 <- mean(distLS>=LSFit1$obsStat) 
pvalueLS2 <- mean(distLS>=LSFit2$obsStat) 

## ----Rank----------------------------------------------------------------
#----------Fit the semiparametric harmonic regression----------#
## Estimate the frequency (Rank-based)
phiRange <- seq(0,0.99,by=0.01) 
system.time(paraRankLS1 <- GetFitRankLS(y1,t=t,lambdaRange,phiRange)) 
system.time(paraRankLS2 <- GetFitRankLS(y2,t=t,lambdaRange,phiRange)) 
paraRankLS1[2,1]
paraRankLS2[2,1]

## Modeling given the estimated frequency
fitSHar1 <- semihregScam(y1, t, lambda0=paraRankLS1[2,1],
                          family = gaussian(), bs = "mpi") #the Gaussian case
fitSHar2 <- semihregScam(y2, t, lambda0=paraRankLS2[2,1],
                          family = gaussian(), bs = "mpi") #the Poisson case

## Fitted values at continuous time
newt <- seq(t[1],t[length(t)],0.1)
predSHar1 <- predict(fitSHar1,newt)[,2]
predSHar2 <- predict(fitSHar2,newt)[,2]

## P-values of the likelihood ratio test
pvalueRank1 <- mean(distRankLS>=paraRankLS1[1,1])
pvalueRank2 <- mean(distRankLS>=paraRankLS2[1,1])

## ----geneFit, fig.show ='hold',fig.cap = "The observed and fitted gene expressions", fig.width=7,fig.height=3----
par(mfrow=c(1,2))

#YCL055W
plot(t,y1,type="p",pch=19,
     xlab = "time", ylab = "gene expression",
     main="YCL055W")
points(t,fitHar1$fitted.values,col="blue",pch=3)
lines(newt,predHar1,col="blue")
points(t,fitSHar1$esty,col="red",pch=6)
lines(newt,predSHar1,col="red",lty=3)

#YMR198W
plot(t,y1,type="p",pch=19,
     xlab = "time", ylab = "gene expression",
     main="YMR198W")
points(t,fitHar2$fitted.values,col="blue",pch=3)
lines(newt,predHar2,col="blue")
points(t,fitSHar2$esty,col="red",pch=6)
lines(newt,predSHar2,col="red",lty=3)

par(mfrow=c(1,1))

