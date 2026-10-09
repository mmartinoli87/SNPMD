rm(list = ls())

# Parallel Monte Carlo SNPMD estimation for multiple Brock-Hommes benchmark
# parameter combinations and sample sizes.
#
# Examples:
#   BH_COMBINATIONS=all BH_CORES=8 Rscript Est_BH_SNPMD_MC_exp.R
#   BH_COMBINATIONS=1,4,9 BH_CORES=4 Rscript Est_BH_SNPMD_MC_exp.R
#   BH_COMBINATIONS=5:8 BH_REPLICATIONS=100 BH_CORES=6 Rscript Est_BH_SNPMD_MC_exp.R
#   BH_SAMPLE_SIZES=250,500,1000 BH_CORES=5 Rscript Est_BH_SNPMD_MC_exp.R
#   BH_SAMPLE_SIZES=500 BH_CORES=5 Rscript Est_BH_SNPMD_MC_exp.R

# Resolve the working directory from this script
vCommandArgs <- commandArgs(trailingOnly = FALSE)
vFileArg <- grep("^--file=", vCommandArgs, value = TRUE)
if (length(vFileArg) > 0L) {
  nScriptPath <- sub("^--file=", "", vFileArg[1L])
  nScriptPath <- gsub("~+~", " ", nScriptPath, fixed = TRUE)
  wd <- dirname(normalizePath(nScriptPath, mustWork = TRUE))
} else {
  wd <- normalizePath(getwd(), mustWork = TRUE)
}
setwd(wd)

if (!requireNamespace("randtoolbox", quietly = TRUE)) {
  stop("Package 'randtoolbox' is required. Install it with install.packages('randtoolbox').")
}

fnEnvInteger <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (length(value) != 1L || !is.finite(value) || value < 1L) {
    stop(name, " must be a positive integer")
  }
  value
}

fnEnvIntegerVector <- function(name, default) {
  specification <- Sys.getenv(name, unset = "")
  if (!nzchar(specification)) return(as.integer(default))

  tokens <- trimws(strsplit(specification, ",", fixed = TRUE)[[1L]])
  values <- suppressWarnings(as.integer(tokens))
  if (!length(values) || any(!nzchar(tokens)) || anyNA(values) || any(values < 1L)) {
    stop(name, " must contain comma-separated positive integers")
  }
  unique(values)
}

fnCombinationIndices <- function(specification, numberAvailable) {
  if (!nzchar(specification) || tolower(specification) == "all") {
    return(seq_len(numberAvailable))
  }

  tokens <- trimws(strsplit(specification, ",", fixed = TRUE)[[1L]])
  indices <- unlist(lapply(tokens, function(token) {
    if (grepl("^[0-9]+:[0-9]+$", token)) {
      endpoints <- as.integer(strsplit(token, ":", fixed = TRUE)[[1L]])
      seq.int(endpoints[1L], endpoints[2L])
    } else if (grepl("^[0-9]+$", token)) {
      as.integer(token)
    } else {
      stop("Invalid BH_COMBINATIONS token: '", token,
           "'. Use 'all', comma-separated indices, or ranges such as 1:4.")
    }
  }), use.names = FALSE)

  indices <- unique(indices)
  if (!length(indices) || any(indices < 1L | indices > numberAvailable)) {
    stop("BH_COMBINATIONS contains indices outside 1:", numberAvailable)
  }
  indices
}

bQuick <- identical(Sys.getenv("BH_QUICK", unset = "0"), "1")
iS <- fnEnvInteger("BH_REPLICATIONS", if (bQuick) 1L else 1000L)
iK <- fnEnvInteger("BH_MAX_POLYNOMIAL_DEGREE", if (bQuick) 2L else 10L)
iP <- fnEnvInteger("BH_SOBOL_POINTS", if (bQuick) 64L else 2^12)
iTau <- fnEnvInteger("BH_TIME_HORIZON", if (bQuick) 300L else 5000L)
iDetectedCores <- parallel::detectCores(logical = FALSE)
if (!is.finite(iDetectedCores)) iDetectedCores <- 1L
iRequestedCores <- fnEnvInteger("BH_CORES", max(1L, iDetectedCores - 1L))

nSampleSizeSpecification <- Sys.getenv("BH_SAMPLE_SIZES", unset = "")
nLegacyN <- Sys.getenv("BH_EMPIRICAL_SIZE", unset = "")
nLegacyM <- Sys.getenv("BH_SIMULATED_SIZE", unset = "")

if (nzchar(nSampleSizeSpecification)) {
  if (nzchar(nLegacyN) || nzchar(nLegacyM)) {
    warning("BH_SAMPLE_SIZES overrides BH_EMPIRICAL_SIZE and BH_SIMULATED_SIZE")
  }
  vSampleSize <- fnEnvIntegerVector(
    "BH_SAMPLE_SIZES",
    if (bQuick) 100L else c(250L, 500L, 1000L)
  )
} else if (nzchar(nLegacyN) || nzchar(nLegacyM)) {
  iLegacyN <- fnEnvInteger(
    "BH_EMPIRICAL_SIZE",
    if (nzchar(nLegacyM)) as.integer(nLegacyM) else if (bQuick) 100L else 500L
  )
  iLegacyM <- fnEnvInteger(
    "BH_SIMULATED_SIZE",
    if (nzchar(nLegacyN)) as.integer(nLegacyN) else if (bQuick) 100L else 500L
  )
  if (iLegacyN != iLegacyM) {
    stop("This multi-sample driver requires N = M. Set equal legacy sizes or use BH_SAMPLE_SIZES.")
  }
  vSampleSize <- iLegacyN
} else {
  vSampleSize <- if (bQuick) 100L else c(250L, 500L, 1000L)
}

if (any(vSampleSize > iTau + 3L)) {
  stop("Every BH_SAMPLE_SIZES value must not exceed BH_TIME_HORIZON + 3")
}

# Parameter search space.
mBounds <- rbind(
  Beta  = c(lower = 0.0,  upper = 10.0),
  Sigma = c(lower = 0.0,  upper = 1.0),
  b_2   = c(lower = -0.3, upper = 0.3),
  g_2   = c(lower = -0.4, upper = 0.4)
)
param_names <- rownames(mBounds)
iNumPar <- nrow(mBounds)

# One common, reproducible Sobol search design for every truth/replication.
iSobolSeed <- 19871028L
mSobolUnit <- randtoolbox::sobol(
  n = iP, dim = iNumPar, init = TRUE,
  scrambling = 0, normal = FALSE
)
set.seed(iSobolSeed)
vSobolShift <- runif(iNumPar)
mSobolUnit <- sweep(mSobolUnit, 2L, vSobolShift, "+") %% 1
mG <- sweep(mSobolUnit, 2L, mBounds[, "upper"] - mBounds[, "lower"], "*")
mG <- sweep(mG, 2L, mBounds[, "lower"], "+")
colnames(mG) <- param_names
iR <- nrow(mG)

# Fixed model inputs and benchmark combinations.
dR <- 0.01
dSigma2 <- 1
dBeta <- c(0, 0.5, 3, 10)
dG <- c(-0.4, 0.4)
dB <- c(-0.3, 0.3)

mGbench <- expand.grid(
  Beta = dBeta,
  g_2 = dG,
  b_2 = dB,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
mThetaTruth <- data.frame(
  Beta = mGbench$Beta,
  Sigma = rep(dSigma2, nrow(mGbench)),
  b_2 = mGbench$b_2,
  g_2 = mGbench$g_2,
  check.names = FALSE
)

vCombination <- fnCombinationIndices(
  Sys.getenv("BH_COMBINATIONS", unset = "all"),
  nrow(mThetaTruth)
)

#' Construct a Brock and Hommes (1998) simulator.
#' Parameter order: c(beta, sigma, b_2, g_2).
fnBHmod <- function(seed = 0L, r = 0.01) {
  seed <- as.integer(seed)
  gross_rate <- 1 + as.numeric(r)

  unpack <- function(theta) {
    theta <- as.numeric(theta)
    if (length(theta) < 4L || (length(theta) - 2L) %% 2L != 0L) {
      stop("theta must be c(beta, sigma, b_2..b_H, g_2..g_H)")
    }
    num_strategies <- 1L + (length(theta) - 2L) %/% 2L
    list(
      beta = theta[1L],
      sigma = theta[2L],
      b = c(0, theta[seq.int(3L, num_strategies + 1L)]),
      g = c(0, theta[seq.int(num_strategies + 2L, length(theta))])
    )
  }

  step <- function(y_lag, theta, shock) {
    y_lag <- as.numeric(y_lag)
    if (length(y_lag) != 3L) stop("The Brock-Hommes step requires three lags")
    p <- unpack(theta)

    beliefs <- p$g * y_lag[3L] + p$b
    lag_belief_term <- p$g * y_lag[1L] + p$b - gross_rate * y_lag[2L]
    utility <- p$beta * (y_lag[3L] - gross_rate * y_lag[2L]) * lag_belief_term
    exp_utility <- exp(pmax(-400, pmin(400, utility)))
    fractions <- exp_utility / sum(exp_utility)

    as.numeric((sum(beliefs * fractions) + p$sigma * shock) / gross_rate)
  }

  simulate <- function(T, theta, rep = 0L, shocks = NULL) {
    T <- as.integer(T)
    if (T < 3L) stop("T must be at least 3")
    unpack(theta)

    if (is.null(shocks)) {
      set.seed(seed + as.integer(rep))
      shocks <- qnorm(runif(T))
    }
    shocks <- as.numeric(shocks)
    if (length(shocks) != T) stop("shocks must have length T")

    y <- numeric(T)
    if (T > 3L) {
      for (t in 4L:T) y[t] <- step(y[(t - 3L):(t - 1L)], theta, shocks[t])
    }
    y
  }

  structure(
    list(seed = seed, r = as.numeric(r), R = gross_rate,
         step = step, simulate = simulate),
    class = "brock_hommes"
  )
}

# Five probabilities entering the Gray distance: Pr(Y > 0), followed by the
# four two-state transition probabilities. The denominator matches the source.
fnProcessProbabilities <- function(vData) {
  vSign <- as.integer(vData > 0)
  iT <- length(vSign)
  if (iT < 3L) stop("At least three observations are required")

  vCombinations <- as.matrix(expand.grid(0:1, 0:1))
  vTransitions <- apply(vCombinations, 1L, function(vState) {
    sum(
      vSign[2L:(iT - 1L)] == vState[1L] &
        vSign[3L:iT] == vState[2L]
    ) / iT
  })
  c(positive = mean(vSign), vTransitions)
}

# Benchmark series and target process probabilities are computed once per
# selected truth and sample size.
lTargetProbabilities <- list()
for (iCombination in vCombination) {
  vThetaStar <- as.numeric(mThetaTruth[iCombination, param_names])
  lBHbench <- fnBHmod(seed = 3L, r = dR)
  vBench <- lBHbench$simulate(iTau + 3L, vThetaStar)
  for (iSampleSize in vSampleSize) {
    nTargetKey <- paste(iCombination, iSampleSize, sep = "_")
    lTargetProbabilities[[nTargetKey]] <- fnProcessProbabilities(
      tail(vBench, iSampleSize)
    )
  }
}

# All truth x replication combinations are independent parallel jobs.
mJobs <- expand.grid(
  combination = vCombination,
  replication = seq_len(iS),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

fnRunJob <- function(iJob) {
  iCombination <- mJobs$combination[iJob]
  s <- mJobs$replication[iJob]

  # Avoid storing full iR x M matrices. Each candidate path is simulated only
  # once, and its tails supply the five probabilities needed for every sample size.
  lSimulatedProbabilities <- setNames(
    lapply(vSampleSize, function(iSampleSize) {
      matrix(NA_real_, nrow = iR, ncol = 5L)
    }),
    as.character(vSampleSize)
  )

  for (j in seq_len(iR)) {
    iSimulationSeed <- 3L + (s - 1L) * iR + j
    lBHsim <- fnBHmod(seed = iSimulationSeed, r = dR)
    vSim <- lBHsim$simulate(iTau + 3L, as.numeric(mG[j, param_names]))
    for (iSampleSize in vSampleSize) {
      lSimulatedProbabilities[[as.character(iSampleSize)]][j, ] <-
        fnProcessProbabilities(tail(vSim, iSampleSize))
    }
  }

  vTruth <- as.numeric(mThetaTruth[iCombination, param_names])
  names(vTruth) <- param_names

  lSampleResults <- vector("list", length(vSampleSize))
  lDiagnostics <- setNames(
    vector("list", length(vSampleSize)),
    as.character(vSampleSize)
  )

  for (iSizeIndex in seq_along(vSampleSize)) {
    iSampleSize <- vSampleSize[iSizeIndex]
    vTarget <- lTargetProbabilities[[paste(iCombination, iSampleSize, sep = "_")]]
    mSimulatedProbabilities <-
      lSimulatedProbabilities[[as.character(iSampleSize)]]

    vDist <- 2 * abs(vTarget[1L] - mSimulatedProbabilities[, 1L]) +
      rowSums(abs(sweep(
        mSimulatedProbabilities[, -1L, drop = FALSE],
        2L, vTarget[-1L], "-"
      )))

    mAIC <- data.frame(Degree = integer(), AIC = numeric())
    for (k in seq_len(iK)) {
      lAIC <- lm(vDist ~ poly(as.matrix(mG), k, raw = TRUE))
      lStep <- step(lAIC, trace = 0)
      mAIC <- rbind(mAIC, data.frame(Degree = k, AIC = AIC(lStep)))
    }
    iOptimalDegree <- mAIC$Degree[which.min(mAIC$AIC)]

    polyFit <- step(
      lm(vDist ~ poly(as.matrix(mG), iOptimalDegree, raw = TRUE)),
      trace = 0
    )
    vfp <- polyFit$fitted.values
    iMinWMD <- which.min(vfp)
    vEstimate <- as.numeric(mG[iMinWMD, param_names])
    names(vEstimate) <- param_names

    lSampleResults[[iSizeIndex]] <- data.frame(
      sample_size = iSampleSize,
      N = iSampleSize,
      M = iSampleSize,
      combination = iCombination,
      replication = s,
      true_Beta = vTruth["Beta"],
      true_Sigma = vTruth["Sigma"],
      true_b_2 = vTruth["b_2"],
      true_g_2 = vTruth["g_2"],
      Beta = vEstimate["Beta"],
      Sigma = vEstimate["Sigma"],
      b_2 = vEstimate["b_2"],
      g_2 = vEstimate["g_2"],
      polynomial_degree = iOptimalDegree,
      min_raw_distance = min(vDist),
      min_fitted_distance = min(vfp),
      row.names = NULL,
      check.names = FALSE
    )

    if (s == iS) {
      lDiagnostics[[as.character(iSampleSize)]] <- list(
        vDist = vDist,
        vfp = vfp,
        AIC = mAIC,
        best_index = iMinWMD
      )
    }
  }

  list(
    result = do.call(rbind, lSampleResults),
    diagnostics = lDiagnostics
  )
}

iCores <- min(iRequestedCores, nrow(mJobs))
cat(sprintf(
  paste0(
    "Running %d job(s): %d truth combination(s) x %d replication(s), ",
    "with %d sample size(s), using %d core(s).\n"
  ),
  nrow(mJobs), length(vCombination), iS, length(vSampleSize), iCores
))
cat("Selected combinations:", paste(vCombination, collapse = ", "), "\n")
cat("Sample sizes (N = M):", paste(vSampleSize, collapse = ", "), "\n")

dStart <- proc.time()[[3L]]
vJobIndex <- seq_len(nrow(mJobs))

if (iCores == 1L) {
  lResults <- lapply(vJobIndex, fnRunJob)
} else if (.Platform$OS.type != "windows") {
  lResults <- parallel::mclapply(
    vJobIndex,
    fnRunJob,
    mc.cores = iCores,
    mc.preschedule = TRUE,
    mc.set.seed = FALSE
  )
} else {
  cl <- parallel::makeCluster(iCores)
  lResults <- tryCatch(
    {
      parallel::clusterExport(
        cl,
        varlist = c(
          "mJobs", "mG", "iR", "vSampleSize", "iTau", "iK", "iS", "dR",
          "param_names", "mThetaTruth", "lTargetProbabilities",
          "fnBHmod", "fnProcessProbabilities", "fnRunJob"
        ),
        envir = environment()
      )
      parallel::parLapply(cl, vJobIndex, fnRunJob)
    },
    finally = parallel::stopCluster(cl)
  )
}

if (any(vapply(lResults, inherits, logical(1L), what = "try-error"))) {
  stop("At least one parallel job failed. Rerun with BH_CORES=1 for a direct error trace.")
}

dElapsed <- proc.time()[[3L]] - dStart
mResults <- do.call(rbind, lapply(lResults, `[[`, "result"))
mResults <- mResults[
  order(mResults$sample_size, mResults$combination, mResults$replication),
]
rownames(mResults) <- NULL

nPath <- path.expand(Sys.getenv(
  "BH_OUTPUT_DIR",
  unset = file.path(wd, "output_parallel")
))
dir.create(nPath, recursive = TRUE, showWarnings = FALSE)

nSampleSizeTag <- paste(vSampleSize, collapse = "-")
pAllCombinedFile <- file.path(
  nPath,
  sprintf(
    "Est_BH_SNPMD_sample_sizes_%s_combined.csv",
    nSampleSizeTag
  )
)
write.csv(mResults, pAllCombinedFile, row.names = FALSE)

for (iSampleSize in vSampleSize) {
  dSize <- mResults[
    mResults$sample_size == iSampleSize,
    ,
    drop = FALSE
  ]
  pSizeCombinedFile <- file.path(
    nPath,
    sprintf(
      "Est_BH_SNPMD_N%d_M%d_combined.csv",
      iSampleSize, iSampleSize
    )
  )
  write.csv(dSize, pSizeCombinedFile, row.names = FALSE)

  for (iCombination in vCombination) {
    dCombination <- dSize[
      dSize$combination == iCombination,
      ,
      drop = FALSE
    ]
    pCombinationFile <- file.path(
      nPath,
      sprintf(
        "Est_BH_SNPMD_N%d_M%d_comb%d.csv",
        iSampleSize, iSampleSize, iCombination
      )
    )
    write.csv(dCombination, pCombinationFile, row.names = FALSE)

    iDiagnosticJob <- which(
      mJobs$combination == iCombination & mJobs$replication == iS
    )
    lDiagnostic <- lResults[[iDiagnosticJob]]$diagnostics[[as.character(iSampleSize)]]
    saveRDS(
      list(
        estimates = dCombination,
        truth = mThetaTruth[iCombination, param_names, drop = FALSE],
        bounds = mBounds,
        sobol_design = mG,
        last_distance = lDiagnostic$vDist,
        last_fitted_distance = lDiagnostic$vfp,
        AIC = lDiagnostic$AIC,
        best_index = lDiagnostic$best_index,
        combination = iCombination,
        configuration = list(
          replications = iS, N = iSampleSize, M = iSampleSize,
          sample_sizes = vSampleSize, time_horizon = iTau,
          sobol_points = iP, maximum_polynomial_degree = iK,
          cores = iCores, sobol_seed = iSobolSeed,
          sobol_shift = vSobolShift
        )
      ),
      file.path(
        nPath,
        sprintf(
          "Diagnostics_BH_SNPMD_N%d_M%d_comb%d.rds",
          iSampleSize, iSampleSize, iCombination
        )
      )
    )
  }
}

saveRDS(
  list(
    results = mResults,
    truths = mThetaTruth[vCombination, , drop = FALSE],
    selected_combinations = vCombination,
    jobs = mJobs,
    elapsed_seconds = dElapsed,
    configuration = list(
      replications = iS, sample_sizes = vSampleSize, time_horizon = iTau,
      sobol_points = iP, maximum_polynomial_degree = iK,
      requested_cores = iRequestedCores, used_cores = iCores,
      sobol_seed = iSobolSeed, sobol_shift = vSobolShift
    )
  ),
  file.path(
    nPath,
    sprintf("Run_BH_SNPMD_sample_sizes_%s.rds", nSampleSizeTag)
  )
)

cat(sprintf("Completed in %.2f seconds.\n", dElapsed))
cat("All-size combined estimates:", normalizePath(pAllCombinedFile), "\n")
cat("Output directory:", normalizePath(nPath), "\n")
