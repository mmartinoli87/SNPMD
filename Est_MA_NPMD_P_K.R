rm(list = ls())

wd <- "~/SNPMD/"
setwd(wd)

iS <- 10000
iN <- 100#0
iM <- 1000
vK <- 2:10
# iP <- c(100,200,400,600,800)
iP <- seq(10, 1000, 10)

# Preallocate: array with dimensions (simulations, length(vK), length(iP))
aEstA <- array(0, dim = c(iS, length(vK), length(iP)),
               dimnames = list(NULL, paste0("K", vK), paste0("P", iP)))

dTheta <- 0.5 # true MA(1) parameter

for (p in seq_along(iP)) {
  iLim <- iP[p]
  
  # vA depends on iLim
  vA <- (-iLim:iLim) / (1 + iLim)
  iR <- length(vA)
  
  # We'll reuse some objects inside s-loop to avoid reallocating repeatedly
  mY <- matrix(0, nrow = iR, ncol = iM)
  
  cat("Starting iP =", iP[p], " (index", p, "of", length(iP), ")\n")
  
  for (s in 1:iS) {
    ## Simulate benchmark MA(1)
    vX <- arima.sim(model = list(ma = dTheta), n = iN, innov = rnorm(iN), start.innov = 0)
    
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
    
    vDist <- 2 * abs(dP1 - vQ1) + apply(abs(mQ - rep(1, iR) %*% t(vP)), 1, sum)
    
    for (k in 1:length(vK)) {
      # Fit polynomial of degree vK[k] and take the vA at which fitted value is minimum
      polyFit <- step(lm(vDist ~ poly(vA, vK[k], raw = TRUE)), trace = 0)
      vfp <- polyFit$fitted.values
      aEstA[s, k, p] <- vA[which.min(vfp)]
    }
    
    if (s %% 1000 == 0) cat("  finished s =", s, "\n")
  } # end s-loop
  
  print(p)
} # end p-loop

# Optionally save the whole 3D array as an RDS for later use
saveRDS(aEstA, file = file.path(wd, "Est_MA_SNPMD_N100_M1000_all.rds"))
