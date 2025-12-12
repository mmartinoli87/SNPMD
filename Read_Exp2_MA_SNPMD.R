rm(list=ls())
library(xtable)

wd <- "~/SSSUP Dropbox/Mario Martinoli/SNPMD/Results_MC_exp/MA_model/"
setwd(wd)

## Parameters
iS <- 10000
iN <- 100#0
iM <- 100#0
vK <- 2:10
iP <- seq(10, 1000, 10)
lA <- list()
for (p in 1:length(iP)) {
  lA[[p]] <- (-iP[p]:iP[p]) / (1 + iP[p])
}
vP <- as.numeric(unlist(lapply(lA,length)))
dTheta <- 0.5

## Read data II
vEst1.II <- read.csv(paste(wd,"Est_MA_II_N100_M100.csv",sep=""),header = T)[,-1]
vEst2.II <- read.csv(paste(wd,"Est_MA_II_N100_M1000.csv",sep=""),header = T)[,-1]
vEst3.II <- read.csv(paste(wd,"Est_MA_II_N1000_M1000.csv",sep=""),header = T)[,-1]

dBias1.II <- mean(vEst1.II)-dTheta
dVar1.II <- var(vEst1.II)
dMSE1.II <- dBias1.II^2+dVar1.II

dBias2.II <- mean(vEst2.II)-dTheta
dVar2.II <- var(vEst2.II)
dMSE2.II <- dBias2.II^2+dVar2.II

dBias3.II <- mean(vEst3.II)-dTheta
dVar3.II <- var(vEst3.II)
dMSE3.II <- dBias3.II^2+dVar3.II

## Read data MLE
vEst1.mle <- read.csv(paste(wd,"Est_MA_MLE_N100.csv",sep=""),header = T)[,-1]
vEst2.mle <- read.csv(paste(wd,"Est_MA_MLE_N1000.csv",sep=""),header = T)[,-1]

dBias1.mle <- mean(vEst1.mle)-dTheta
dVar1.mle <- var(vEst1.mle)
dMSE1.mle <- dBias1.mle^2+dVar1.mle

dBias2.mle <- mean(vEst2.mle)-dTheta
dVar2.mle <- var(vEst2.mle)
dMSE2.mle <- dBias2.mle^2+dVar2.mle

## Read data SMM
vEst1.SMM <- read.csv(paste(wd,"Est_MA_SMM_N100_M100.csv",sep=""),header = T)[,-1]
vEst2.SMM <- read.csv(paste(wd,"Est_MA_SMM_N100_M1000.csv",sep=""),header = T)[,-1]
vEst3.SMM <- read.csv(paste(wd,"Est_MA_SMM_N1000_M1000.csv",sep=""),header = T)[,-1]

dBias1.SMM <- mean(vEst1.SMM)-dTheta
dVar1.SMM <- var(vEst1.SMM)
dMSE1.SMM <- dBias1.SMM^2+dVar1.SMM

dBias2.SMM <- mean(vEst2.SMM)-dTheta
dVar2.SMM <- var(vEst2.SMM)
dMSE2.SMM <- dBias2.SMM^2+dVar2.SMM

dBias3.SMM <- mean(vEst3.SMM)-dTheta
dVar3.SMM <- var(vEst3.SMM)
dMSE3.SMM <- dBias3.SMM^2+dVar3.SMM
 
## Read data and make plot
aEst1 <- readRDS(paste(wd,"Est_MA_SNPMD_N100_M100_all.rds",sep=""))
mEst1 <- matrix(0,length(iP),length(vK))
mVar1 <- matrix(0,length(iP),length(vK))
for(p in seq_along(iP)) {
  mEst1[p,] <- apply(aEst1[,,p],2,mean)
  mVar1[p,] <- apply(aEst1[,,p],2,var)
}
mBias1 <- mEst1-dTheta
mMSE1  <- mVar1+mBias1^2

aEst2 <- readRDS(paste(wd,"Est_MA_SNPMD_N100_M1000_all.rds",sep=""))
mEst2 <- matrix(0,length(iP),length(vK))
mVar2 <- matrix(0,length(iP),length(vK))
for(p in seq_along(iP)) {
  mEst2[p,] <- apply(aEst2[,,p],2,mean)
  mVar2[p,] <- apply(aEst2[,,p],2,var)
}
mBias2 <- mEst2-dTheta
mMSE2  <- mVar2+mBias2^2

aEst3 <- readRDS(paste(wd,"Est_MA_SNPMD_N1000_M1000_all.rds",sep=""))
mEst3 <- matrix(0,length(iP),length(vK))
mVar3 <- matrix(0,length(iP),length(vK))
for(p in seq_along(iP)) {
  mEst3[p,] <- apply(aEst3[,,p],2,mean)
  mVar3[p,] <- apply(aEst3[,,p],2,var)
}
mBias3 <- mEst3-dTheta
mMSE3  <- mVar3+mBias3^2

par(mfrow = c(3, 1), 
    mar = c(0, 0, 0, 0),  # Reduced margins: bottom=5, left=4, top=2, right=1
    oma = c(4, 4, 1, 1),  # Outer margin at bottom for legend
    xpd = NA)

plot(mMSE1[,1],type="l",lwd=2,ylim=c(0,0.5),xlab="",ylab="",xaxt="n")
for (k in 1:5) {
  lines(mMSE1[,k],lwd=2,lty=k)
}
segments(x0=1,y0=dMSE1.II,x1=100,y1=dMSE1.II,lwd=2,col="darkgray")
segments(x0=1,y0=dMSE1.mle,x1=100,y1=dMSE1.mle,lwd=2,lty=2,col="darkgray")
segments(x0=1,y0=dMSE1.SMM,x1=100,y1=dMSE1.SMM,lwd=2,lty=3,col="darkgray")
legend("topright",legend=c(paste0("K=",2:6),"Ind. Inf.","MLE","SMM"),
       lty = c(1:5, 1, 2, 3),col = c(rep("black",5), rep("darkgray",3)),lwd = 2,
       bty = "n", cex = 0.8)

plot(mMSE2[,1],type="l",lwd=2,ylim=c(0,0.5),xlab="",ylab="",xaxt="n")
for (k in 1:5) {
  lines(mMSE2[,k],lwd=2,lty=k)
}
segments(x0=1,y0=dMSE2.II,x1=100,y1=dMSE2.II,lwd=2,col="darkgray")
segments(x0=1,y0=dMSE1.mle,x1=100,y1=dMSE1.mle,lwd=2,lty=2,col="darkgray")
segments(x0=1,y0=dMSE2.SMM,x1=100,y1=dMSE2.SMM,lwd=2,lty=3,col="darkgray")

plot(mMSE3[,1],type="l",lwd=2,ylim=c(0,0.5),xlab="",ylab="",xaxt="n")
axis(side = 1, at = c(1,10,20,30,40,50,60,70,80,90,100),
     labels=c(vP[1],vP[10],vP[20],vP[30],vP[40],vP[50],
              vP[60],vP[70],vP[80],vP[90],vP[100]))
for (k in 1:5) {
  lines(mMSE3[,k],lwd=2,lty=k)
}
segments(x0=1,y0=dMSE3.II,x1=100,y1=dMSE3.II,lwd=2,col="darkgray")
segments(x0=1,y0=dMSE2.mle,x1=100,y1=dMSE2.mle,lwd=2,lty=2,col="darkgray")
segments(x0=1,y0=dMSE3.SMM,x1=100,y1=dMSE3.SMM,lwd=2,lty=3,col="darkgray")

# Reset layout to default
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 0, 0), xpd = FALSE)

