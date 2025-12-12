rm(list=ls())

iS <- 1000
nPath <- "/Users/Mario/PhD/Papers/EstSimMod/New/Results_MC_exp/BH/"
mData <- read.csv(paste(nPath,"Est_SVM_emp_SP500_S10000.csv",sep=""),header=T)
head(mData)
# mData <- read.csv(paste(nPath,"Est_SVM_emp_EU50_S10000.csv",sep=""))

vBeta <- mData[,2]
vG <- mData[,3]

## Compute mean
dBeta <- mean(vBeta)
dG <- mean(vG)

## Compute SEdSE.beta <- sd(vBeta)/sqrt(iS)dSE.g <- sd(vG)/sqrt(iS)