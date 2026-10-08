rm(list = ls())

wd <- "~"
setwd(wd)

iS <- 10000
iN <- 1000
iM <- 1000
vK <- 2:6
## Fifty logarithmically spaced actual grid sizes in [21, 2001]
vGridSize <- 2L * round(
  (exp(seq(log(21), log(2001), length.out = 50L)) - 1) / 2
) + 1L

stopifnot(
  length(vGridSize) == 50L,
  length(unique(vGridSize)) == 50L,
  min(vGridSize) == 21L,
  max(vGridSize) == 2001L,
  all(vGridSize %% 2L == 1L)
)

## Half-width needed by the existing symmetric-grid construction
vHalfGrid <- (vGridSize - 1L) %/% 2L

# Preallocate: array with dimensions (simulations, length(vK), length(iP))
aEstA <- array(
  NA_real_,
  dim = c(iS, length(vK), length(vGridSize)),
  dimnames = list(
    replication = NULL,
    degree = paste0("K", vK),
    grid_size = paste0("P", vGridSize)
  )
)

dTheta <- 0.5 # true MA(1) parameter

for (p in seq_along(vGridSize)) {
  iLim <- vHalfGrid[p]
  
  vA <- (-iLim:iLim) / (1 + iLim)
  iR <- length(vA)
  
  stopifnot(iR == vGridSize[p])
  
  mY <- matrix(0, nrow = iR, ncol = iM)
  
  cat(
    "Starting P =", vGridSize[p],
    "(index", p, "of", length(vGridSize), ")\n"
  )
  
  for (s in 1:iS) {
    ## Simulate benchmark MA(1)
    vX <- arima.sim(model = list(ma = dTheta), n = iN, innov = rnorm(iN), start.innov = 0)
    
    ## Simulated paths for each vA value (each row is a path)
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
saveRDS(
  list(
    estimates = aEstA,
    grid_size = vGridSize,
    half_grid = vHalfGrid,
    degree = vK,
    replications = iS,
    N = iN,
    M = iM,
    theta_true = dTheta
  ),
  file = file.path(
    wd,
    paste0("Est_MA_SNPMD_N", iN, "_M", iM, "_logP.rds")
  )
)
