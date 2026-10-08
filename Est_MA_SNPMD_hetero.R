rm(list = ls())

wd <- "~"
setwd(wd)

iS <- 10000
iN <- 100#0
iM <- 100#0
iK <- 6
iP <- 100

dTheta <- 0.5 # true MA(1) parameter

iLim <- iP
vA <- (-iLim:iLim) / (1 + iLim)
iR <- length(vA)
mDist <- matrix(0, iS, iR)
mfp <- matrix(0, iS, iR)
vEst <- rep(0, iS)
  
# We'll reuse some objects inside s-loop to avoid reallocating repeatedly
mY <- matrix(0, nrow = iR, ncol = iM)

## Simulate benchmark MA(1) for conditional experiment
vX.cond <- arima.sim(model = list(ma = dTheta), n = iN, innov = rnorm(iN), start.innov = 0)
for (s in 1:iS) {
  # ## Simulate benchmark MA(1) for unconditional experiment
  # vX <- arima.sim(model = list(ma = dTheta), n = iN, innov = rnorm(iN), start.innov = 0)
  vX <- vX.cond
    
  ## Simulated paths for each vA value (each row is a path)
  # NOTE: innovations length must match n = iM
  for (j in 1:iR) {
    mY[j, ] <- arima.sim(model = list(ma = vA[j]), n = iM, innov = rnorm(iM), start.innov = 0)
  }
    
  # Create grids and polynomials coefficients and compute Gray's distance
  dP1 <- mean(vX > 0)
  vQ1 <- apply(1 * (mY > 0), 1, mean)
    
  logsum <- function(vV, vSeq) sum(vV[2:(length(vV) - 1)] == vSeq[1] & vV[3:length(vV)] == vSeq[2])
    
  vP <- numeric(4)
  mQ <- matrix(0, nrow = iR, ncol = 4)
  vComb <- expand.grid(0:1, 0:1)
  for (i in 1:4) {
    vP[i] <- logsum(1 * (vX > 0), as.numeric(vComb[i, ])) / iN
    mQ[, i] <- apply(1 * (mY > 0), 1, logsum, vSeq = as.numeric(vComb[i, ])) / iM
  }
    
  mDist[s, ] <- 2 * abs(dP1 - vQ1) + apply(abs(mQ - rep(1, iR) %*% t(vP)), 1, sum)
  
  # Loop over different polynomial degrees
  mAIC <- data.frame(Degree=integer(),AIC=numeric())
  
  for (k in 1:iK) {  # Test polynomial degrees from 1 to 10
    lAIC <- lm(mDist[s, ] ~ poly(vA,k,raw = TRUE))
    lStep <- step(lAIC,trace=0) # Use stepwise selection
    mAIC <- rbind(mAIC,data.frame(Degree=k,AIC=AIC(lStep)))
  }
  
  # Find the degree with the lowest AIC
  vP.opt <- mAIC$Degree[which.min(mAIC$AIC)]
  
  # Fit polynomial regression
  polyFit <- step(lm(mDist[s, ] ~ poly(vA, vP.opt, raw = TRUE)),trace = 0)
  mfp[s, ] <- polyFit$fitted.values
  
  vEst[s] <- vA[which.min(mfp[s, ])]
    
  print(s)
} # end s-loop

## Save files
# save(mDist,mfp,vEst, file = paste(wd,"/output/Est_MA_SNPMD_N",iN,"_M",iM,"_P",iR,"_het.RData",sep=""))
save(mDist,mfp,vEst, file = paste(wd,"/output/Est_MA_SNPMD_N",iN,"_M",iM,"_P",iR,"_het_cond.RData",sep=""))

