rm(list=ls())

wd <- "/Users/Mario/PhD/Papers/EstSimMod/New/Codes"
setwd(wd)

source("FnBHmod.R")
 
iS <- 10000
iK <- 10 # Max degree of polynomials

# Load data
# mData <- read.csv("SP500.csv")
# vData <- read.csv("SP500.csv",header=T)[,2]
mData <- read.csv("EURO_STOXX_50.csv")
vData <- read.csv("EURO_STOXX_50.csv",header=T)[,2]
vData <- as.numeric(vData)

# Function to replace NA values with moving averages
fnFill <- function(prices,ma_window=NULL) {
  # Loop through each element of the vector
  for (i in seq_along(prices)) {
    if (is.na(prices[i])) { # If the value is missing
      # Subset the previous non-missing values
      previous_values <- prices[1:(i-1)][!is.na(prices[1:(i-1)])]
      
      # Compute the moving average (use all previous or last N values based on ma_window)
      if (!is.null(ma_window) && length(previous_values) >= ma_window) {
        prices[i] <- mean(tail(previous_values, ma_window))
      } else if (length(previous_values) > 0) {
        prices[i] <- mean(previous_values)
      }
    }
  }
  return(prices)
}

# # SP500 from 2011-12-31 to 2012-12-31
# vPr <- fnFill(vData,ma_window=1)[which(mData$DATE=="2012-01-02"):which(mData$DATE=="2012-12-31")]
# vY <- diff(log(vPr))

# EURO STOXX 50 from 2011-12-30 to 2012-12-31
vPr <- fnFill(vData,ma_window=1)
vY <- diff(log(vPr))

# Sample size
iTau <- 5000
vN <- 1:length(vPr)
iN <- length(vPr)
vM <- 1:length(vPr)
iM <- length(vPr)

# Grids for parameters
vA1 <- seq(1,30,length.out=101) # Beta
vA2 <- seq(-1,1,length.out=101) # g
mG <- expand.grid(vA1,vA2)
iR <- nrow(mG)

# Model parameters
iTau <- 5000
dN1_0 <- 0.5
iW <- 1
dR <- 0.1/250
dSigma2 <- 0.1
dBeta <- 2.25
dG <- 2
dAlpha <- 23
dP_star <- 0.7 #0.8
dEps.mu <- 0.1
dEps.sigma <- 0.1

# True parameter values
vTheta.star <- c(dBeta,dG)

vTheta.hat1 <- rep(0,iS)
vTheta.hat2 <- rep(0,iS)

for (s in 1:iS) {
  # set.seed(s)
  
  ## Simulated paths
  mZ <- matrix(0,nrow=iR,ncol=iM)
  for(j in 1:nrow(mG)) {
  	vP.sim <- fnBHmod(iTau,dEps.mu,dEps.sigma,dN1_0,iW,dR,dSigma2,mG[j,1],mG[j,2],dAlpha,dP_star)[[1]][(iTau-iM):iTau]
  	mZ[j,] <- diff(log(vP.sim))
  	# print(j)
  }
  
  # Create grids and polynomials coefficients and compute Gray's distance
  dP1 <- mean(vY>0)
  vQ1 <- apply(1*(mZ>0),1,mean)
  logsum <- function(vV,vSeq) sum(vV[2:(length(vV)-1)]==vSeq[1] & vV[3:length(vV)]==vSeq[2])
  vP <- rep(0,4)
  mQ <- matrix(0,nrow=iR,ncol=4)
  vComb <- expand.grid(0:1,0:1)
  for (i in 1:4) {
    vP[i] <- logsum(1*(vY>0),as.numeric(vComb[i,]))/iN
    mQ[,i] <- apply(1*(mZ>0),1,logsum,vSeq=as.numeric(vComb[i,]))/iM
  }
  vDist <- 2*abs(dP1-vQ1)+apply(abs(mQ-rep(1,iR)%*%t(vP)),1,sum)
  mDist <- matrix(vDist,length(vA1),length(vA2))
  
  # Loop over different polynomial degrees
  mAIC <- data.frame(Degree=integer(),AIC=numeric())
  
  for (k in 1:iK) {  # Test polynomial degrees from 1 to 10
    lAIC <- lm(vDist ~ poly(as.matrix(mG),k,raw = TRUE))
    lStep <- step(lAIC,trace=0) # Use stepwise selection
    mAIC <- rbind(mAIC,data.frame(Degree=k,AIC=AIC(lStep)))
  }
  
  # Find the degree with the lowest AIC
  iP.opt <- mAIC$Degree[which.min(mAIC$AIC)]
  
  # Fit polynomial regression
  polyFit <- step(lm(vDist ~ poly(as.matrix(mG),iP.opt,raw=TRUE)))
  vC <- polyFit$coefficients
  vfp <- polyFit$fitted.values
  mfp <- matrix(vfp,length(vA1),length(vA2))
  
  ## Estimates
  min.WMD  <- which(mDist == min(mDist), arr.ind = TRUE)
  vTheta.hat1[s] <- vA1[min.WMD[1]]
  vTheta.hat2[s] <- vA2[min.WMD[2]]
  
  print(s)
}

## Save files
nPath <- "/Users/Mario/PhD/Papers/EstSimMod/New/Results_MC_exp/BH/"
# pFile <- paste(nPath,"Est_BH_emp_SP500_S10000.csv",sep="")
pFile <- paste(nPath,"Est_BH_emp_EU50_S10000.csv",sep="")
mEst <- cbind(vTheta.hat1,vTheta.hat2) 
write.csv(mEst,file=pFile)

## Compute sample mean
dBeta <- mean(mEst[,1])
dG <- mean(mEst[,2])

## Compute SE
dSE.beta <- sd(mEst[,1])/sqrt(iS)
dSE.g <- sd(mEst[,2])/sqrt(iS)

# Persp
mObj <- mfp
mObj2 <- mDist
nbcol <- 200
color <- topo.colors(nbcol)
mObj <- mObj[-1,-1]+mObj[-1,-ncol(mObj)]+mObj[-nrow(mObj),-1]+mObj[-nrow(mObj),-ncol(mObj)]
mObj2 <- mObj2[-1,-1]+mObj2[-1,-ncol(mObj2)]+mObj2[-nrow(mObj2),-1]+mObj2[-nrow(mObj2),-ncol(mObj2)]
facetcol <- cut(mObj,nbcol)
facetcol2 <- cut(mObj2,nbcol)
#
pdf(file=paste(nPath,"Obj_fc_BH.pdf",sep=""))
persp(vA1[-1],vA2[-1],mObj,theta=60,phi=30,xlab="Beta",
      ylab="g",zlab="",main="",col=color[facetcol],ticktype="detailed")
dev.off()
pdf(file=paste(nPath,"Obj_fc_BH_rough.pdf",sep=""))
persp(vA1[-1],vA2[-1],mObj2,theta=60,phi=30,xlab="Beta",
      ylab="g",zlab="",main="",col=color[facetcol2],ticktype="detailed")
dev.off()
