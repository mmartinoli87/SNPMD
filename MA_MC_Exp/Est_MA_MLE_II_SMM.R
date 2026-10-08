rm(list=ls())
library(tseries)
library(xtable)

wd <- "~/SSSUP Dropbox/Mario Martinoli/SNPMD/Codes/"
setwd(wd)

# Set seed for replication
set.seed(123)

# Source functions
source("Mom_functions.R")

# Define variables
iN <- 1000
iM <- 1000
#iR <- 10000
iS <- 10000

dTheta <- 0.5

###########################
# True process is a MA(1) #
###########################

mEps <- matrix(rnorm(iN*iS),iN,iS)

mY <-mat.or.vec(iN,iS)
mY[1,] <- array(0,c(1,iS))

for(s in 1:iS){
  for(n in 2:iN){
    mY[n,s] <- mEps[n,s]-dTheta*mEps[n-1,s]
	}
}
# plot.ts(mY[,1])

#####################################
# MLE estimation with Kalman filter #
#####################################

vTheta.mle <- numeric(iS)
for (s in 1:iS) {
  lEst.mle <- arma(mY[,s],order = c(0,1))
  vTheta.mle[s] <- -lEst.mle$coef[1]
  print(s)
}

dBias.mle <- mean(vTheta.mle)-dTheta
dVar.mle <- var(vTheta.mle)
dMSE.mle <- dVar.mle+(dBias.mle^2)

nFold <- "~/SSSUP Dropbox/Mario Martinoli/SNPMD/Results_MC_exp/MA_model/"
pFile <- paste(nFold,"Est_MA_MLE_N",iN,".csv",sep="")
write.csv(vTheta.mle,file=pFile)

####################################
# Indirect Inference with an AR(3) #
####################################

vY <- mY[,1]
lEst.aux <- arma(vY,order=c(3,0))
vBeta.hat <- -lEst.aux$coef[1:3]

# mBeta.hat <- mat.or.vec(iS,3)
# for (s in 1:iS){
#   lEst.aux <- arma(mY[,s],order=c(3,0))
# 	mBeta.hat[s,] <- -lEst.aux$coef[1:3]
# }

vBind.rw <- vBeta.hat
# vBind.rw <- c(mean(mBeta.hat[,1]),mean(mBeta.hat[,2]),mean(mBeta.hat[,3]))

mEps.sim <- matrix(rnorm(iM*iS),iM,iS)

fnII <- function(theta0,error,M,bind){
  vZ <- numeric(M)
	vZ[1] <- 0
	for (m in 2:M){
		vZ[m] <- error[m]-theta0*error[m-1]
	}
	lFit.sim <- arma(vZ,order=c(3,0))
	dBeta.sim1 <- -lFit.sim$coef[1]
	dBeta.sim2 <- -lFit.sim$coef[2]
	dBeta.sim3 <- -lFit.sim$coef[3]
	ssr <- (bind[1]-dBeta.sim1)^2+(bind[2]-dBeta.sim2)^2+(bind[3]-dBeta.sim3)^2
}

vTheta.hat <- numeric(iS)
for (s in 1:iS){
	lEst.II <- optimize(fnII,c(0,1),error=mEps.sim[,s],M=iM,bind=vBind.rw)
	vTheta.hat[s] <- lEst.II$minimum
	print(s)
}

nFold <- "~/SSSUP Dropbox/Mario Martinoli/SNPMD/Results_MC_exp/MA_model/"
pFile <- paste(nFold,"Est_MA_II_N",iN,"_M",iM,".csv",sep="")
write.csv(vTheta.hat,file=pFile)

#################################################
# SMM estimation (see Michaelides and Ng, 2000) #
#################################################

# Benchmark series
vY <- mY[,1]

# Benchmark statistics
dS1 <- fnMom1(vY)
dS2 <- fnMom2(vY)
dS3 <- fnACOV1(vY)
dS4 <- fnACOV2(vY)
dS5 <- fnACOV3(vY)

mEps.sim <- matrix(rnorm(iM*iS),iM,iS)

fnSMM <- function(theta0,error,M) {
  vZ <- numeric(M)
  vZ[1] <- 0
  for (m in 2:M){
    vZ[m] <- error[m]-theta0*error[m-1]
  }
  dS1.sim <- fnMom1(vZ)
  dS2.sim <- fnMom2(vZ)
  dS3.sim <- fnACOV1(vZ)
  dS4.sim <- fnACOV2(vZ)
  dS5.sim <- fnACOV3(vZ)
  of1 <- (dS1-dS1.sim)^2
  of2 <- (dS2-dS2.sim)^2
  of3 <- (dS3-dS3.sim)^2
  of4 <- (dS4-dS4.sim)^2
  of5 <- (dS5-dS5.sim)^2
  return(sum(of1+of2+of3+of4+of5))
}

vTheta.smm <- numeric(iS)
for (s in 1:iS){
  lEst.SMM <- optimize(fnSMM,c(0,1),error=mEps.sim[,s],M=iM)
  vTheta.smm[s] <- lEst.SMM$minimum
  print(s)
}

nFold <- "~/SSSUP Dropbox/Mario Martinoli/SNPMD/Results_MC_exp/MA_model/"
pFile <- paste(nFold,"Est_MA_SMM_N",iN,"_M",iM,".csv",sep="")
write.csv(vTheta.smm,file=pFile)

#########################################
# Compute statistics and produce tables #
#########################################

vEst <- read.csv(pFile,header=T)[,2]
vTheta.hat <- vEst

dTheta.hat <- mean(vTheta.hat)
dTheta.hat
dVAR <- var(vTheta.hat)
dSD <- sd(vTheta.hat)
dSD

dBias <- mean(vTheta.hat)-dTheta
dBias
dMSE <- mean((vTheta.hat-dTheta)^2)
dMSE
dRMSE <- sqrt(mean((vTheta.hat-dTheta)^2))
dRMSE

dAlpha <- 0.1 # 90% confidence interval

dLow1 <- quantile(vEst,probs=dAlpha/2)
dUp1 <- quantile(vEst,probs=1-dAlpha/2)
dLow.Hall1 <- 2*mean(vEst)-dUp1
dUp.Hall1 <- 2*mean(vEst)-dLow1
dAcc1 <- length(which(vEst>dLow.Hall1&vEst<dUp.Hall1))/iR*100

dfDec <- data.frame(iN,iM,dBias,dSD,dRMSE)
colnames(dfDec) <- c("N","M","Bias","SD","RMSE")

nTex <- paste(nFold,"Dec_II_N",iN,"_M",iM,".tex",sep="")
print(xtable(round(dfDec,digits=4),type="latex",digits=4),file=nTex,include.rownames=T)




