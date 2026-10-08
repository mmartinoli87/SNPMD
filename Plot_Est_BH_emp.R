rm(list = ls())

## Function for surface colours
fnSurfaceColours <- function(mZ, palette = "YlOrRd", nColours = 100L) {
  
  ## Average height of each rectangular facet
  mFacet <- (
    mZ[-1L, -1L] +
      mZ[-1L, -ncol(mZ)] +
      mZ[-nrow(mZ), -1L] +
      mZ[-nrow(mZ), -ncol(mZ)]
  ) / 4
  
  vPalette <- hcl.colors(nColours, palette)
  
  vFinite <- is.finite(mFacet)
  
  if (!any(vFinite)) {
    return(matrix("grey80", nrow(mFacet), ncol(mFacet)))
  }
  
  vRange <- range(mFacet[vFinite])
  
  if (diff(vRange) == 0) {
    mIndex <- matrix(
      ceiling(nColours / 2),
      nrow(mFacet),
      ncol(mFacet)
    )
  } else {
    vIndex <- rep(NA_integer_, length(mFacet))
    
    vIndex[vFinite] <- as.integer(cut(
      mFacet[vFinite],
      breaks = seq(vRange[1L], vRange[2L], length.out = nColours + 1L),
      include.lowest = TRUE
    ))
    
    mIndex <- matrix(
      vIndex,
      nrow = nrow(mFacet),
      ncol = ncol(mFacet)
    )
  }
  
  mColours <- matrix(
    "grey80",
    nrow(mFacet),
    ncol(mFacet)
  )
  
  mColours[is.finite(mIndex)] <-
    vPalette[mIndex[is.finite(mIndex)]]
  
  mColours
}

### SNPMD surface plots

## Select the empirical configuration
nMarket <- "SP500"       # alternatively: "EURO_STOXX_50"
iYears  <- 4L            # 1, 2, or 4
iGrid   <- 50L

nRDS <- paste0(
  "~/BH_SNPMD_empirical_diagnostics.rds"
)

lResults <- readRDS(nRDS)

param_names <- c("Beta", "Sigma", "b_2", "g_2")
mBounds <- lResults$bounds
mSobol <- lResults$sobol_design

## Locate the selected market/sample-size diagnostics
iJob <- which(
  lResults$jobs$market == nMarket &
    lResults$jobs$years == iYears
)

if (length(iJob) != 1L) {
  stop("The selected market/year configuration was not found.")
}

lDiagnostic <- lResults$last_replication_diagnostics[[iJob]]

## Monte Carlo mean estimate: the other two parameters are held here
vThetaHat <- unlist(
  lResults$mean_estimates[
    lResults$mean_estimates$market == nMarket &
      lResults$mean_estimates$years == iYears,
    param_names,
    drop = FALSE
  ],
  use.names = TRUE
)

## Degree selected by AIC in the saved diagnostic replication
iDegree <- lDiagnostic$AIC$Degree[
  which.min(lDiagnostic$AIC$AIC)
]

## Refit the polynomial-smoothed Gray-distance surface
mBasis <- poly(
  as.matrix(mSobol[, param_names]),
  degree = iDegree,
  raw = TRUE
)

lSurfaceFit <- lm.fit(
  x = cbind(`(Intercept)` = 1, mBasis),
  y = lDiagnostic$raw_distance
)

vCoefficient <- lSurfaceFit$coefficients
vCoefficient[is.na(vCoefficient)] <- 0

## All six parameter pairs
mPairs <- combn(param_names, 2L)

fnSNPMDSurface <- function(xName, yName) {
  
  vX <- seq(
    mBounds[xName, "lower"],
    mBounds[xName, "upper"],
    length.out = iGrid
  )
  
  vY <- seq(
    mBounds[yName, "lower"],
    mBounds[yName, "upper"],
    length.out = iGrid
  )
  
  mPlotGrid <- expand.grid(x = vX, y = vY)
  
  ## Hold the other parameters at their MC mean estimates
  mThetaGrid <- matrix(
    rep(vThetaHat, each = nrow(mPlotGrid)),
    nrow = nrow(mPlotGrid),
    dimnames = list(NULL, param_names)
  )
  
  mThetaGrid[, xName] <- mPlotGrid$x
  mThetaGrid[, yName] <- mPlotGrid$y
  
  ## Evaluate the fitted response surface on the regular plotting grid
  mNewBasis <- poly(
    as.matrix(mThetaGrid[, param_names]),
    degree = iDegree,
    raw = TRUE
  )
  
  vObjective <- drop(
    cbind(`(Intercept)` = 1, mNewBasis) %*% vCoefficient
  )
  
  mObjective <- matrix(
    vObjective,
    nrow = length(vX),
    ncol = length(vY)
  )
  
  list(
    x = vX,
    y = vY,
    z = mObjective,
    xName = xName,
    yName = yName
  )
}

lSurfaces <- lapply(
  seq_len(ncol(mPairs)),
  function(i) fnSNPMDSurface(mPairs[1L, i], mPairs[2L, i])
)

nOutput <- file.path(
  dirname(nRDS),
  sprintf("SNPMD_pairwise_objective_persp_%s_%dyr.pdf", nMarket, iYears)
)

pdf(nOutput, width = 11, height = 14)

par(
  mfrow = c(3, 2),
  mar = c(2.8, 2.8, 3.0, 1.0),
  oma = c(0, 0, 2.0, 0)
)

for (lSurface in lSurfaces) {
  
  mZ <- lSurface$z
  
  persp(
    x = lSurface$x,
    y = lSurface$y,
    z = mZ,
    theta = 45,
    phi = 25,
    expand = 0.70,
    col = fnSurfaceColours(
      mZ,
      palette = "YlOrRd"
    ),
    border = "grey35",
    shade = 0.25,
    ticktype = "detailed",
    xlab = lSurface$xName,
    ylab = lSurface$yName,
    zlab = "Gray distance",
    main = paste(
      lSurface$xName,
      "and",
      lSurface$yName
    )
  )
}

mtext(
  sprintf(
    "%s: smoothed empirical SNPMD objective, %d-year sample (degree %d)",
    nMarket,
    iYears,
    iDegree
  ),
  outer = TRUE,
  cex = 1.1
)

dev.off()

cat("SNPMD perspective figure saved to:\n", nOutput, "\n")

### NPSMLE surface plot

## Select the empirical configuration
nMarket <- "SP500"       # alternatively: "EURO_STOXX_50"
iYears  <- 4L            # 1, 2, or 4
iGrid   <- 50L
iCores  <- 5L

nRoot <- paste0(
  "~"
)

## Change I100 below if BH_MAXIT was changed
nRunDirectory <- file.path(
  nRoot,
  "...",
  paste0(
    "..."
  )
)

nRDS <- file.path(
  nRunDirectory,
  "BH_NPSMLE_empirical_results.rds"
)

if (!file.exists(nRDS)) {
  stop("NPSMLE results file not found: ", nRDS)
}

source(file.path(nRoot, "BH_functions.R"))

lResults <- readRDS(nRDS)

param_names <- c("Beta", "Sigma", "b_2", "g_2")
mBounds <- matrix(c(0,0,-1,-2,10,1,1,2),
                  length(param_names),2)
colnames(mBounds) <- c("lower","upper")
rownames(mBounds) <- param_names

## Monte Carlo mean estimate
dEstimate <- lResults$mean_estimates[
  lResults$mean_estimates$market == nMarket &
    lResults$mean_estimates$years == iYears,
  ,
  drop = FALSE
]

if (nrow(dEstimate) != 1L) {
  stop("The selected market/year configuration was not found.")
}

vThetaHat <- unlist(
  dEstimate[1L, param_names, drop = FALSE],
  use.names = TRUE
)

iN <- as.integer(dEstimate$N)
iM <- as.integer(dEstimate$M)
dScale <- as.numeric(dEstimate$data_scale)

## Recover the empirical estimation sample
nDataFile <- file.path(
  nRunDirectory,
  paste0("BH_NPSMLE_", nMarket, "_analysis_data.csv")
)

dData <- read.csv(
  nDataFile,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dData$Date <- as.Date(dData$Date)

dTraining <- dData[
  dData$Date >= as.Date(dEstimate$training_start) &
    dData$Date <= as.Date(dEstimate$training_end),
  ,
  drop = FALSE
]

dTraining <- tail(dTraining, iN)

if (nrow(dTraining) != iN) {
  stop("The reconstructed empirical sample does not contain N observations.")
}

## Use the same scaling employed in estimation
vY <- dTraining$deviation / dScale

dRiskFreeRate <- lResults$settings$risk_free_rate
lBH <- brock_hommes(seed = 3L, r = dRiskFreeRate)

## Fixed common random numbers for all parameter evaluations

fnKernelLogLikelihood <- function(vY, vTheta, vEpsilon, iNumLags = 3L) {
  
  iT <- length(vY)
  iDraws <- length(vEpsilon)
  
  vContribution <- numeric(iT - iNumLags)
  
  for (t in seq.int(iNumLags + 1L, iT)) {
    
    vSimulated <- lBH$step(
      vY[(t - iNumLags):(t - 1L)],
      vTheta,
      vEpsilon
    )
    
    vDifference <- vSimulated - vY[t]
    
    dBandwidth <- population_sd(vDifference) *
      (4 / (3 * iDraws))^(1 / 5)
    
    if (!is.finite(dBandwidth) || dBandwidth <= 0) {
      return(-Inf)
    }
    
    vLogKernel <- dnorm(
      vDifference,
      mean = 0,
      sd = dBandwidth,
      log = TRUE
    )
    
    dMaximum <- max(vLogKernel)
    
    vContribution[t - iNumLags] <-
      dMaximum + log(mean(exp(vLogKernel - dMaximum)))
  }
  
  ## Average log-likelihood, making values comparable across N
  mean(vContribution)
}

## All six parameter pairs
mPairs <- combn(param_names, 2L)

iBaseInnovationSeed <- 20260917L

fnNPSMLESurface <- function(xName, yName, iPair) {
  
  cat("Computing:", xName, "and", yName, "\n")
  
  vX <- seq(
    mBounds[xName, "lower"],
    mBounds[xName, "upper"],
    length.out = iGrid
  )
  
  vYGrid <- seq(
    mBounds[yName, "lower"],
    mBounds[yName, "upper"],
    length.out = iGrid
  )
  
  mPlotGrid <- expand.grid(
    x = vX,
    y = vYGrid
  )
  
  ## Parameters not displayed are held at their MC mean estimates
  mThetaGrid <- matrix(
    rep(vThetaHat, each = nrow(mPlotGrid)),
    nrow = nrow(mPlotGrid),
    dimnames = list(NULL, param_names)
  )
  
  mThetaGrid[, xName] <- mPlotGrid$x
  mThetaGrid[, yName] <- mPlotGrid$y
  
  fnEvaluate <- function(i) {
    
    ## A different innovation vector is generated at every grid point.
    ## The seed depends only on the pair and grid point, so results are
    ## reproducible and independent of parallel scheduling.
    dSeed <- (
      as.double(iBaseInnovationSeed) +
        as.double(iPair) * 1000003 +
        as.double(i) * 10007
    ) %% 2147483646
    
    set.seed(as.integer(dSeed + 1))
    
    vEpsilon_i <- qnorm(runif(iM))
    
    fnKernelLogLikelihood(
      vY = vY,
      vTheta = mThetaGrid[i, param_names],
      vEpsilon = vEpsilon_i
    )
  }
  
  if (
    iCores > 1L &&
    .Platform$OS.type != "windows"
  ) {
    
    vLogLikelihood <- unlist(
      parallel::mclapply(
        seq_len(nrow(mThetaGrid)),
        fnEvaluate,
        mc.cores = iCores,
        mc.preschedule = TRUE,
        mc.set.seed = FALSE
      ),
      use.names = FALSE
    )
    
  } else {
    
    vLogLikelihood <- vapply(
      seq_len(nrow(mThetaGrid)),
      fnEvaluate,
      numeric(1L)
    )
  }
  
  vLogLikelihood[!is.finite(vLogLikelihood)] <- NA_real_
  
  mLogLikelihood <- matrix(
    vLogLikelihood,
    nrow = length(vX),
    ncol = length(vYGrid)
  )
  
  list(
    x = vX,
    y = vYGrid,
    z = mLogLikelihood,
    xName = xName,
    yName = yName
  )
}

## Build surfaces
lSurfaces <- lapply(
  seq_len(ncol(mPairs)),
  function(i) {
    fnNPSMLESurface(
      xName = mPairs[1L, i],
      yName = mPairs[2L, i],
      iPair = i
    )
  }
)

## Plot surfaces
nOutput <- file.path(
  nRunDirectory,
  sprintf(
    "NPSMLE_pairwise_loglikelihood_%s_%dyr.pdf",
    nMarket,
    iYears
  )
)

pdf(nOutput, width = 11, height = 14)

par(
  mfrow = c(3, 2),
  mar = c(2.8, 2.8, 3.0, 1.0),
  oma = c(0, 0, 2.0, 0)
)

for (lSurface in lSurfaces) {
  
  mZ <- lSurface$z
  
  ## persp() requires a finite plotting range. Invalid likelihood
  ## evaluations, for example at Sigma = 0, are placed at the lowest
  ## finite height solely for graphical purposes.
  bFinite <- is.finite(mZ)
  
  if (!any(bFinite)) {
    stop(
      "No finite likelihood values for ",
      lSurface$xName,
      " and ",
      lSurface$yName
    )
  }
  
  dLowerSurface <- min(mZ[bFinite])
  
  mZ[!bFinite] <- dLowerSurface
  
  persp(
    x = lSurface$x,
    y = lSurface$y,
    z = mZ,
    theta = 45,
    phi = 25,
    expand = 0.70,
    col = fnSurfaceColours(
      mZ,
      palette = "YlOrRd"
    ),
    border = "grey35",
    shade = 0.25,
    ticktype = "detailed",
    xlab = lSurface$xName,
    ylab = lSurface$yName,
    zlab = "Mean log-likelihood",
    main = paste(
      lSurface$xName,
      "and",
      lSurface$yName
    )
  )
}

mtext(
  sprintf(
    "%s: empirical NPSMLE log-likelihood, %d-year sample",
    nMarket,
    iYears
  ),
  outer = TRUE,
  cex = 1.1
)

dev.off()

cat("NPSMLE perspective figure saved to:\n", nOutput, "\n")
