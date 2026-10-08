rm(list = ls())

# Empirical four-parameter SNPMD estimation of the Brock-Hommes model.
#
# The script estimates the BH model on nested 1-, 2-, and 4-year windows
# (250 observations per trading year by default), always imposing M = N.
# The complete estimator is repeated 1,000 times with new simulation innovations
# to measure Monte Carlo uncertainty.
#
# Production example:
#   BH_CORES=5 Rscript Est_BH_SNPMD_emp.R
#
# Useful overrides:
#   BH_YEARS=1,2,4 BH_DAYS_PER_YEAR=250 BH_ESTIMATION_END=2012-12-31 \
#   Rscript Est_BH_SNPMD_empirical_4par_parallel.R
#
# The observed index is expressed as a deviation from a 61-observation moving average.

# -----------------------------------------------------------------------------
# Paths and environment options
# -----------------------------------------------------------------------------

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

fnEnvInteger <- function(name, default, allowZero = FALSE) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  lower <- if (allowZero) 0L else 1L
  if (length(value) != 1L || is.na(value) || value < lower) {
    stop(name, " must be ", if (allowZero) "a non-negative" else "a positive", " integer")
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

fnEnvFlag <- function(name, default = TRUE) {
  value <- tolower(trimws(Sys.getenv(name, unset = if (default) "1" else "0")))
  if (!value %in% c("0", "1", "false", "true", "no", "yes")) {
    stop(name, " must be one of 0, 1, false, true, no, or yes")
  }
  value %in% c("1", "true", "yes")
}

fnEnvChoice <- function(name, default, choices) {
  value <- tolower(trimws(Sys.getenv(name, unset = default)))
  if (!value %in% choices) {
    stop(name, " must be one of: ", paste(choices, collapse = ", "))
  }
  value
}

bQuick <- identical(Sys.getenv("BH_QUICK", unset = "0"), "1")
iS <- fnEnvInteger("BH_REPLICATIONS", if (bQuick) 2L else 1000L)
vYears <- fnEnvIntegerVector("BH_YEARS", c(1L, 2L, 4L))
iDaysPerYear <- fnEnvInteger("BH_DAYS_PER_YEAR", if (bQuick) 50L else 250L)
vSampleSize <- vYears * iDaysPerYear
names(vSampleSize) <- as.character(vYears)
iValidationSize <- fnEnvInteger(
  "BH_VALIDATION_SIZE", if (bQuick) 50L else iDaysPerYear, allowZero = TRUE
)
iK <- fnEnvInteger("BH_MAX_POLYNOMIAL_DEGREE", if (bQuick) 2L else 10L)
iP <- fnEnvInteger("BH_SOBOL_POINTS", if (bQuick) 64L else 2^12)
iBurnIn <- fnEnvInteger("BH_BURN_IN", if (bQuick) 100L else 4000L)
iDetectedCores <- parallel::detectCores(logical = FALSE)
if (!is.finite(iDetectedCores)) iDetectedCores <- 1L
iRequestedCores <- fnEnvInteger("BH_CORES", max(1L, iDetectedCores - 1L))
dEstimationEnd <- as.Date(Sys.getenv("BH_ESTIMATION_END", unset = "2012-12-31"))
if (is.na(dEstimationEnd)) stop("BH_ESTIMATION_END must have format YYYY-MM-DD")
iHalfWindow <- fnEnvInteger("BH_HALF_WINDOW", 30L)
nFundamentalMode <- fnEnvChoice(
  "BH_FUNDAMENTAL_MODE", "trailing", c("centered", "trailing")
)
bAutoDownload <- fnEnvFlag("BH_AUTO_DOWNLOAD", TRUE)
bMakePlots <- fnEnvFlag("BH_MAKE_PLOTS", TRUE)
bResume <- fnEnvFlag("BH_RESUME", TRUE)

nOutputPath <- path.expand(Sys.getenv(
  "BH_OUTPUT_DIR", unset = file.path(wd, "output_empirical_parallel")
))
dir.create(nOutputPath, recursive = TRUE, showWarnings = FALSE)
nCheckpointPath <- file.path(nOutputPath, "checkpoints")
dir.create(nCheckpointPath, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Market data
# -----------------------------------------------------------------------------

mMarketSpecification <- data.frame(
  market = c("SP500", "EURO_STOXX_50"),
  label = c("S&P 500", "Euro Stoxx 50"),
  yahoo_symbol = c("^GSPC", "^STOXX50E"),
  default_file = c(
    file.path(wd, "data", "SP500.csv"),
    file.path(wd, "data", "EURO_STOXX_50.csv")
  ),
  stringsAsFactors = FALSE
)

mMarketSpecification$input_file <- c(
  path.expand(Sys.getenv("BH_SP500_FILE", unset = mMarketSpecification$default_file[1L])),
  path.expand(Sys.getenv(
    "BH_EURO_STOXX50_FILE", unset = mMarketSpecification$default_file[2L]
  ))
)

nMarketSelection <- trimws(strsplit(
  Sys.getenv("BH_MARKETS", unset = "SP500,EURO_STOXX_50"), ",", fixed = TRUE
)[[1L]])
nMarketSelection <- toupper(gsub("[- ]", "_", nMarketSelection))
nMarketSelection[nMarketSelection %in% c("EU50", "EUROSTOXX50", "STOXX50E")] <-
  "EURO_STOXX_50"
if (any(!nMarketSelection %in% mMarketSpecification$market)) {
  stop("BH_MARKETS must contain SP500 and/or EURO_STOXX_50")
}
mMarketSpecification <- mMarketSpecification[
  match(unique(nMarketSelection), mMarketSpecification$market), , drop = FALSE
]

fnParseDate <- function(x) {
  if (inherits(x, "Date")) return(x)
  x <- trimws(as.character(x))
  formats <- c("%Y-%m-%d", "%d/%m/%Y", "%m/%d/%Y", "%Y/%m/%d")
  candidates <- lapply(formats, function(fmt) as.Date(x, format = fmt))
  counts <- vapply(candidates, function(z) sum(!is.na(z)), integer(1L))
  candidates[[which.max(counts)]]
}

fnReadPriceCSV <- function(file, market) {
  if (!file.exists(file)) stop("File does not exist: ", file)
  dRaw <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(dRaw)) stop("No observations in ", file)

  nCleanNames <- gsub("[^a-z0-9]", "", tolower(names(dRaw)))
  iDate <- which(nCleanNames %in% c("date", "datetime", "timestamp"))[1L]
  if (is.na(iDate)) {
    stop("The file for ", market, " needs an explicit Date column")
  }

  dDate <- fnParseDate(dRaw[[iDate]])
  if (sum(!is.na(dDate)) < 0.9 * nrow(dRaw)) {
    stop("Could not parse at least 90% of the dates in ", file)
  }

  vPreferred <- if (market == "SP500") {
    c("sp500", "adjusted", "adjclose", "close", "price", "value")
  } else {
    c("eurostoxx50", "stoxx50e", "sx5e", "adjusted", "adjclose", "close", "price", "value")
  }
  iPrice <- match(vPreferred, nCleanNames, nomatch = 0L)
  iPrice <- iPrice[iPrice > 0L & iPrice != iDate]
  if (!length(iPrice)) {
    iNumeric <- which(vapply(dRaw, function(z) {
      converted <- suppressWarnings(as.numeric(as.character(z)))
      mean(is.finite(converted)) > 0.9
    }, logical(1L)))
    iNumeric <- setdiff(iNumeric, iDate)
    if (!length(iNumeric)) stop("Could not identify a price column in ", file)
    iPrice <- iNumeric[1L]
  } else {
    iPrice <- iPrice[1L]
  }

  vPrice <- suppressWarnings(as.numeric(as.character(dRaw[[iPrice]])))
  dData <- data.frame(Date = dDate, Price = vPrice)
  dData <- dData[is.finite(dData$Price) & !is.na(dData$Date), , drop = FALSE]
  dData <- dData[order(dData$Date), , drop = FALSE]
  dData <- dData[!duplicated(dData$Date, fromLast = TRUE), , drop = FALSE]
  rownames(dData) <- NULL
  if (nrow(dData) < 3L) stop("Too few usable observations in ", file)
  dData
}

fnDownloadPrices <- function(symbol, market, file) {
  if (!requireNamespace("quantmod", quietly = TRUE)) {
    stop(
      "The local ", market, " file is missing or too short and package 'quantmod' ",
      "is unavailable. Install quantmod or provide a long dated CSV through the ",
      if (market == "SP500") "BH_SP500_FILE" else "BH_EURO_STOXX50_FILE",
      " environment variable."
    )
  }

  iPastDays <- ceiling((max(vSampleSize) + 2L * iHalfWindow + 50L) * 365.25 / 250)
  iFutureDays <- ceiling((iValidationSize + 2L * iHalfWindow + 50L) * 365.25 / 250)
  dFrom <- dEstimationEnd - iPastDays
  dTo <- dEstimationEnd + iFutureDays
  cat(sprintf("Downloading %s (%s) from %s to %s ...\n",
              market, symbol, dFrom, dTo))

  x <- suppressWarnings(quantmod::getSymbols(
    symbol, src = "yahoo", from = dFrom, to = dTo + 1L,
    auto.assign = FALSE, warnings = FALSE
  ))
  p <- tryCatch(quantmod::Ad(x), error = function(e) NULL)
  if (is.null(p) || all(!is.finite(as.numeric(p)))) p <- quantmod::Cl(x)
  dData <- data.frame(
    Date = as.Date(zoo::index(p)), Price = as.numeric(p),
    stringsAsFactors = FALSE
  )
  dData <- dData[is.finite(dData$Price) & !is.na(dData$Date), , drop = FALSE]
  if (!nrow(dData)) stop("No usable prices were downloaded for ", market)

  nDownloadedFile <- file.path(nOutputPath, paste0(market, "_downloaded_prices.csv"))
  write.csv(dData, nDownloadedFile, row.names = FALSE)
  cat("Saved downloaded prices to ", nDownloadedFile, "\n", sep = "")
  dData
}

fnConstructDeviations <- function(dPrice) {
  y <- dPrice$Price
  n <- length(y)
  w <- iHalfWindow

  if (nFundamentalMode == "centered") {
    if (n <= 2L * w) stop("The price series is too short for the centred MA window")
    centres <- seq.int(w + 1L, n - w)
    fundamental <- vapply(
      centres, function(i) mean(y[(i - w):(i + w)]), numeric(1L)
    )
  } else {
    iWindow <- 2L * w + 1L
    if (n < iWindow) stop("The price series is too short for the trailing MA window")
    centres <- seq.int(iWindow, n)
    fundamental <- vapply(
      centres, function(i) mean(y[(i - iWindow + 1L):i]), numeric(1L)
    )
  }

  data.frame(
    Date = dPrice$Date[centres],
    observed = y[centres],
    fundamental = fundamental,
    deviation = y[centres] - fundamental,
    stringsAsFactors = FALSE
  )
}

fnExtractWindows <- function(dDeviation, market) {
  dPre <- dDeviation[dDeviation$Date <= dEstimationEnd, , drop = FALSE]
  if (nrow(dPre) < max(vSampleSize)) {
    stop(
      market, " has only ", nrow(dPre), " usable deviations on or before ",
      dEstimationEnd, "; at least ", max(vSampleSize), " are required"
    )
  }

  lTraining <- setNames(lapply(vSampleSize, function(iN) tail(dPre, iN)), names(vSampleSize))
  if (iValidationSize > 0L) {
    dPost <- dDeviation[dDeviation$Date > dEstimationEnd, , drop = FALSE]
    if (nrow(dPost) < iValidationSize) {
      stop(
        market, " has only ", nrow(dPost), " usable post-estimation deviations; ",
        iValidationSize, " are required for validation"
      )
    }
    dValidation <- head(dPost, iValidationSize)
  } else {
    dValidation <- NULL
  }
  list(training = lTraining, validation = dValidation)
}

fnLoadMarket <- function(iMarket) {
  specification <- mMarketSpecification[iMarket, ]
  nCachedFile <- file.path(
    nOutputPath, paste0(specification$market, "_downloaded_prices.csv")
  )
  dPrice <- tryCatch(
    fnReadPriceCSV(specification$input_file, specification$market),
    error = function(e) structure(list(message = conditionMessage(e)), class = "data_error")
  )

  fnTryPrepare <- function(data) {
    tryCatch({
      dDeviation <- fnConstructDeviations(data)
      lWindows <- fnExtractWindows(dDeviation, specification$market)
      list(prices = data, deviations = dDeviation, windows = lWindows)
    }, error = function(e) structure(list(message = conditionMessage(e)), class = "data_error"))
  }

  lPrepared <- if (inherits(dPrice, "data_error")) dPrice else fnTryPrepare(dPrice)
  nDataSource <- specification$input_file
  if (inherits(lPrepared, "data_error") && file.exists(nCachedFile)) {
    dCachedPrice <- tryCatch(
      fnReadPriceCSV(nCachedFile, specification$market),
      error = function(e) structure(list(message = conditionMessage(e)), class = "data_error")
    )
    lCachedPrepared <- if (inherits(dCachedPrice, "data_error")) {
      dCachedPrice
    } else {
      fnTryPrepare(dCachedPrice)
    }
    if (!inherits(lCachedPrepared, "data_error")) {
      dPrice <- dCachedPrice
      lPrepared <- lCachedPrepared
      nDataSource <- nCachedFile
      cat("Using cached downloaded prices: ", nCachedFile, "\n", sep = "")
    }
  }
  if (inherits(lPrepared, "data_error")) {
    if (!bAutoDownload) stop(lPrepared$message)
    warning(
      "Local data for ", specification$market, " cannot support this design (",
      lPrepared$message, "). Downloading a dated series from Yahoo Finance."
    )
    dPrice <- fnDownloadPrices(
      specification$yahoo_symbol, specification$market, specification$input_file
    )
    lPrepared <- fnTryPrepare(dPrice)
    if (inherits(lPrepared, "data_error")) stop(lPrepared$message)
    nDataSource <- nCachedFile
  }

  lPrepared$market <- specification$market
  lPrepared$label <- specification$label
  lPrepared$source_file <- nDataSource
  lPrepared
}

lMarkets <- lapply(seq_len(nrow(mMarketSpecification)), fnLoadMarket)
names(lMarkets) <- vapply(lMarkets, `[[`, character(1L), "market")

# -----------------------------------------------------------------------------
# Brock-Hommes model, Gray distance, and Sobol design
# -----------------------------------------------------------------------------

mBounds <- rbind(
  Beta  = c(lower = 0.0,  upper = 10.0),
  Sigma = c(lower = 0.0,  upper = 1.0),
  b_2   = c(lower = -1, upper = 1),
  g_2   = c(lower = -2, upper = 2)
)
param_names <- rownames(mBounds)
iNumPar <- nrow(mBounds)
dR <- as.numeric(Sys.getenv("BH_RISK_FREE_RATE", unset = "0.0001"))
if (length(dR) != 1L || !is.finite(dR) || dR <= -1) {
  stop("BH_RISK_FREE_RATE must be a finite number greater than -1")
}

iSobolSeed <- fnEnvInteger("BH_SOBOL_SEED", 19871028L)
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

fnBHmod <- function(seed = 0L, r = 0.0001) {
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

fnProcessProbabilities <- function(vData) {
  vSign <- as.integer(vData > 0)
  iT <- length(vSign)
  if (iT < 3L) stop("At least three observations are required")
  mCombinations <- as.matrix(expand.grid(previous = 0:1, current = 0:1))
  vTransitions <- apply(mCombinations, 1L, function(vState) {
    sum(
      vSign[2L:(iT - 1L)] == vState[1L] &
        vSign[3L:iT] == vState[2L]
    ) / iT
  })
  setNames(
    c(mean(vSign), vTransitions),
    c("Pr(Y>0)", "Pr(0,0)", "Pr(1,0)", "Pr(0,1)", "Pr(1,1)")
  )
}

fnGrayDistance <- function(vTarget, mCandidate) {
  if (is.null(dim(mCandidate))) mCandidate <- matrix(mCandidate, nrow = 1L)
  2 * abs(vTarget[1L] - mCandidate[, 1L]) +
    rowSums(abs(sweep(mCandidate[, -1L, drop = FALSE], 2L, vTarget[-1L], "-")))
}

# -------------------
# Repeated estimation
# -------------------

mJobs <- expand.grid(
  market = names(lMarkets),
  years = vYears,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
mJobs$N <- unname(vSampleSize[as.character(mJobs$years)])
mJobs$M <- mJobs$N

fnFitResponseSurface <- function(vDist) {
  mAIC <- data.frame(Degree = seq_len(iK), AIC = NA_real_)
  lModels <- vector("list", iK)
  for (k in seq_len(iK)) {
    lBase <- lm(vDist ~ poly(as.matrix(mG), k, raw = TRUE))
    lModels[[k]] <- step(lBase, trace = 0)
    mAIC$AIC[k] <- AIC(lModels[[k]])
  }
  iOptimalDegree <- mAIC$Degree[which.min(mAIC$AIC)]
  lFit <- lModels[[iOptimalDegree]]
  list(
    degree = iOptimalDegree,
    AIC = mAIC,
    fit = lFit,
    fitted = as.numeric(fitted(lFit))
  )
}

# Within every replication, each Sobol point receives its own innovation stream.
# Seeds depend on both replication and theta, so the entire SNPMD estimate is
# independently repeated while remaining invariant to parallel scheduling.
iMaximumSampleSize <- max(vSampleSize)
fnSimulateDesignPoint <- function(j, s) {
  iSimulationSeed <- 3L + (s - 1L) * iR + j
  lBHsim <- fnBHmod(seed = iSimulationSeed, r = dR)
  vSim <- lBHsim$simulate(
    iBurnIn + iMaximumSampleSize + 3L,
    as.numeric(mG[j, param_names])
  )
  setNames(lapply(vSampleSize, function(iM) {
    fnProcessProbabilities(tail(vSim, iM))
  }), names(vSampleSize))
}

fnEstimateEmpiricalJob <- function(iJob, s, lSimulatedProbabilities,
                                   keepDiagnostics = FALSE) {
  nMarket <- mJobs$market[iJob]
  iYears <- mJobs$years[iJob]
  iN <- mJobs$N[iJob]
  lMarket <- lMarkets[[nMarket]]
  dTraining <- lMarket$windows$training[[as.character(iYears)]]
  stopifnot(nrow(dTraining) == iN, iN == mJobs$M[iJob])

  mCandidate <- lSimulatedProbabilities[[as.character(iYears)]]
  vDist <- fnGrayDistance(vTarget, mCandidate)
  lSurface <- fnFitResponseSurface(vDist)
  vfp <- lSurface$fitted
  iRawMinimum <- which.min(vDist)
  iSmoothMinimum <- which.min(vfp)
  vEstimate <- as.numeric(mG[iSmoothMinimum, param_names])
  names(vEstimate) <- param_names

  mResult <- data.frame(
    market = nMarket,
    market_label = lMarket$label,
    years = iYears,
    replication = s,
    training_start = min(dTraining$Date),
    training_end = max(dTraining$Date),
    N = iN,
    M = iN,
    Beta = vEstimate["Beta"],
    Sigma = vEstimate["Sigma"],
    b_2 = vEstimate["b_2"],
    g_2 = vEstimate["g_2"],
    polynomial_degree = lSurface$degree,
    min_raw_distance = min(vDist),
    min_fitted_distance = min(vfp),
    raw_min_Beta = mG[iRawMinimum, "Beta"],
    raw_min_Sigma = mG[iRawMinimum, "Sigma"],
    raw_min_b_2 = mG[iRawMinimum, "b_2"],
    raw_min_g_2 = mG[iRawMinimum, "g_2"],
    row.names = NULL,
    check.names = FALSE
  )

  list(
    result = mResult,
    diagnostics = if (keepDiagnostics) list(
      target_probabilities = vTarget,
      raw_distance = vDist,
      fitted_distance = vfp,
      AIC = lSurface$AIC,
      raw_minimum_index = iRawMinimum,
      smooth_minimum_index = iSmoothMinimum
    ) else NULL
  )
}

vRunSignature <- list(
  markets = names(lMarkets), years = vYears, sample_sizes = vSampleSize,
  estimation_end = as.character(dEstimationEnd), fundamental_mode = nFundamentalMode,
  half_window = iHalfWindow, risk_free_rate = dR, sobol_points = iP,
  sobol_seed = iSobolSeed, sobol_shift = vSobolShift, burn_in = iBurnIn,
  maximum_polynomial_degree = iK
)

fnCheckpointFile <- function(s) {
  file.path(nCheckpointPath, sprintf("replication_%04d.rds", s))
}

fnReadCheckpoint <- function(s) {
  if (!bResume) return(NULL)
  nFile <- fnCheckpointFile(s)
  if (!file.exists(nFile)) return(NULL)
  lSaved <- tryCatch(readRDS(nFile), error = function(e) NULL)
  if (is.null(lSaved) || !identical(lSaved$signature, vRunSignature) ||
      !identical(lSaved$replication, as.integer(s))) {
    return(NULL)
  }
  lSaved$value
}

fnRunReplication <- function(s) {
  lDesignProbabilities <- lapply(
    seq_len(iR), fnSimulateDesignPoint, s = s
  )
  lSimulatedProbabilities <- setNames(
    lapply(names(vSampleSize), function(nYear) {
      do.call(rbind, lapply(lDesignProbabilities, `[[`, nYear))
    }),
    names(vSampleSize)
  )

  lJobResults <- lapply(seq_len(nrow(mJobs)), function(iJob) {
    fnEstimateEmpiricalJob(
      iJob, s, lSimulatedProbabilities,
      keepDiagnostics = s == iS
    )
  })
  lValue <- list(
    estimates = do.call(rbind, lapply(lJobResults, `[[`, "result")),
    diagnostics = if (s == iS) lapply(lJobResults, `[[`, "diagnostics") else NULL
  )

  saveRDS(
    list(signature = vRunSignature, replication = as.integer(s), value = lValue),
    fnCheckpointFile(s)
  )
  if (s == 1L || s %% 10L == 0L || s == iS) {
    cat(sprintf("Completed replication %d/%d.\n", s, iS))
  }
  lValue
}

dStart <- proc.time()[[3L]]
lReplicationResults <- vector("list", iS)
vPending <- integer()
for (s in seq_len(iS)) {
  lSaved <- fnReadCheckpoint(s)
  if (is.null(lSaved)) {
    vPending <- c(vPending, s)
  } else {
    lReplicationResults[[s]] <- lSaved
  }
}

iRestored <- iS - length(vPending)
iCores <- min(iRequestedCores, max(1L, length(vPending)))
cat(sprintf(
  paste0(
    "Running %d SNPMD replication(s) with %d Sobol points and N=M in {%s}; ",
    "%d checkpoint(s) restored and %d pending, using %d core(s).\n"
  ),
  iS, iR, paste(vSampleSize, collapse = ", "),
  iRestored, length(vPending), iCores
))

if (length(vPending)) {
  if (iCores == 1L) {
    lNewResults <- lapply(vPending, fnRunReplication)
  } else if (.Platform$OS.type != "windows") {
    lNewResults <- parallel::mclapply(
      vPending, fnRunReplication,
      mc.cores = iCores, mc.preschedule = TRUE, mc.set.seed = FALSE
    )
  } else {
    cl <- parallel::makeCluster(iCores)
    lNewResults <- tryCatch({
      parallel::clusterExport(
        cl,
        c(
          "iR", "iS", "iBurnIn", "iMaximumSampleSize", "vSampleSize",
          "mG", "param_names", "dR", "iK", "mJobs", "lMarkets",
          "nCheckpointPath", "vRunSignature",
          "fnBHmod", "fnProcessProbabilities", "fnGrayDistance",
          "fnFitResponseSurface", "fnSimulateDesignPoint",
          "fnEstimateEmpiricalJob", "fnCheckpointFile", "fnRunReplication"
        ),
        envir = environment()
      )
      parallel::parLapply(cl, vPending, fnRunReplication)
    }, finally = parallel::stopCluster(cl))
  }
  lReplicationResults[vPending] <- lNewResults
}

if (any(vapply(lReplicationResults, inherits, logical(1L), what = "try-error"))) {
  stop("At least one replication failed. Rerun with BH_CORES=1 for a direct error trace.")
}

mEstimates <- do.call(rbind, lapply(lReplicationResults, `[[`, "estimates"))
mEstimates <- mEstimates[
  order(mEstimates$market, mEstimates$years, mEstimates$replication),
  , drop = FALSE
]
rownames(mEstimates) <- NULL

lEstimateGroups <- split(
  mEstimates,
  interaction(mEstimates$market, mEstimates$years, drop = TRUE)
)
lMCSummary <- lapply(lEstimateGroups, function(dGroup) {
  do.call(rbind, lapply(param_names, function(nParameter) {
    vEstimate <- dGroup[[nParameter]]
    data.frame(
      market = dGroup$market[1L],
      market_label = dGroup$market_label[1L],
      years = dGroup$years[1L],
      N = dGroup$N[1L],
      M = dGroup$M[1L],
      parameter = nParameter,
      replications = length(vEstimate),
      mean = mean(vEstimate),
      sd = sd(vEstimate),
      mcse = sd(vEstimate) / sqrt(length(vEstimate)),
      median = median(vEstimate),
      q025 = unname(quantile(vEstimate, 0.025)),
      q975 = unname(quantile(vEstimate, 0.975)),
      row.names = NULL,
      check.names = FALSE
    )
  }))
})
mMCSummary <- do.call(rbind, lMCSummary)
mMCSummary <- mMCSummary[
  order(mMCSummary$market, mMCSummary$years, match(mMCSummary$parameter, param_names)),
  , drop = FALSE
]
rownames(mMCSummary) <- NULL

lConfigurationSummary <- lapply(lEstimateGroups, function(dGroup) {
  dResult <- data.frame(
    market = dGroup$market[1L], market_label = dGroup$market_label[1L],
    years = dGroup$years[1L], N = dGroup$N[1L], M = dGroup$M[1L],
    replications = nrow(dGroup), row.names = NULL, check.names = FALSE
  )
  for (nParameter in param_names) {
    vEstimate <- dGroup[[nParameter]]
    dResult[[paste0("mean_", nParameter)]] <- mean(vEstimate)
    dResult[[paste0("sd_", nParameter)]] <- sd(vEstimate)
    dResult[[paste0("mcse_", nParameter)]] <- sd(vEstimate) / sqrt(nrow(dGroup))
  }
  dResult$median_polynomial_degree <- median(dGroup$polynomial_degree)
  dResult$mean_min_raw_distance <- mean(dGroup$min_raw_distance)
  dResult$mean_min_fitted_distance <- mean(dGroup$min_fitted_distance)
  dResult
})
mConfigurationSummary <- do.call(rbind, lConfigurationSummary)
mConfigurationSummary <- mConfigurationSummary[
  order(mConfigurationSummary$market, mConfigurationSummary$years), , drop = FALSE
]
rownames(mConfigurationSummary) <- NULL

# The replication mean is the reported empirical point estimate
mPointEstimates <- do.call(rbind, lapply(lEstimateGroups, function(dGroup) {
  dPoint <- dGroup[1L, c(
    "market", "market_label", "years", "training_start", "training_end", "N", "M"
  ), drop = FALSE]
  for (nParameter in param_names) dPoint[[nParameter]] <- mean(dGroup[[nParameter]])
  dPoint
}))
mPointEstimates <- mPointEstimates[
  order(mPointEstimates$market, mPointEstimates$years), , drop = FALSE
]
rownames(mPointEstimates) <- NULL

for (nMarket in names(lMarkets)) {
  write.csv(
    lMarkets[[nMarket]]$deviations,
    file.path(nOutputPath, paste0("BH_SNPMD_", nMarket, "_analysis_data.csv")),
    row.names = FALSE
  )
}
write.csv(mEstimates, file.path(nOutputPath, "BH_SNPMD_empirical_estimates.csv"), row.names = FALSE)
write.csv(
  mPointEstimates,
  file.path(nOutputPath, "BH_SNPMD_empirical_mean_estimates.csv"),
  row.names = FALSE
)
write.csv(
  mMCSummary,
  file.path(nOutputPath, "BH_SNPMD_empirical_mc_standard_errors.csv"),
  row.names = FALSE
)
write.csv(
  mConfigurationSummary,
  file.path(nOutputPath, "BH_SNPMD_empirical_mc_summary_wide.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    estimates = mEstimates,
    mean_estimates = mPointEstimates,
    mc_summary = mMCSummary,
    mc_summary_wide = mConfigurationSummary,
    jobs = mJobs,
    last_replication_diagnostics = lReplicationResults[[iS]]$diagnostics,
    bounds = mBounds,
    sobol_design = mG,
    configuration = list(
      replications = iS,
      years = vYears,
      days_per_year = iDaysPerYear,
      sample_sizes = vSampleSize,
      estimation_end = dEstimationEnd,
      fundamental_mode = nFundamentalMode,
      risk_free_rate = dR,
      sobol_points = iP,
      sobol_seed = iSobolSeed,
      sobol_shift = vSobolShift,
      burn_in = iBurnIn,
      maximum_polynomial_degree = iK,
      cores = iRequestedCores,
      resume = bResume
    )
  ),
  file.path(nOutputPath, "BH_SNPMD_empirical_diagnostics.rds")
)
