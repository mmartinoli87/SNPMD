rm(list=ls())
library(xtable)

## Folder containing the Monte Carlo results
nPath <- paste0(
  "~"
)

## Locate combination-specific CSV files
vFiles <- list.files(
  path = nPath,
  pattern = "^Est_BH_SNPMD_N[0-9]+_M[0-9]+_comb[0-9]+\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)

if (length(vFiles) == 0L) {
  stop("No combination-specific CSV files found in: ", nPath)
}

## Columns to retain from every CSV
vHead <- c(
  "Beta",
  "Sigma",
  "b_2",
  "g_2",
  "polynomial_degree"
)

## Import files and retain the same five columns
lEstimates <- setNames(
  lapply(vFiles, function(nFile) {
    df <- read.csv(
      nFile,
      header = TRUE,
      check.names = FALSE
    )
    
    vMissing <- setdiff(vHead, names(df))
    
    if (length(vMissing) > 0L) {
      stop(
        "Missing columns in ", basename(nFile), ": ",
        paste(vMissing, collapse = ", ")
      )
    }
    
    df[, vHead, drop = FALSE]
  }),
  basename(vFiles)
)

## Verify dimensions
stopifnot(
  all(vapply(lEstimates, ncol, integer(1L)) == length(vHead)),
  all(vapply(
    lEstimates,
    function(df) identical(names(df), vHead),
    logical(1L)
  ))
)

## Combine all files
dfEstimates <- do.call(rbind, lEstimates)
rownames(dfEstimates) <- NULL

str(lEstimates)
dim(dfEstimates)

## Bias, SD, RMSE and median polynomial degree
vParameters <- c("Beta", "Sigma", "b_2", "g_2")

## True parameter configurations
mBenchmark <- expand.grid(
  Beta = c(0, 0.5, 3, 10),
  g_2  = c(-0.4, 0.4),
  b_2  = c(-0.3, 0.3),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

dfTruth <- data.frame(
  combination = seq_len(nrow(mBenchmark)),
  Beta  = mBenchmark$Beta,
  Sigma = 1,
  b_2   = mBenchmark$b_2,
  g_2   = mBenchmark$g_2,
  check.names = FALSE
)

## Extract N, M and combination number from a filename
fnParseConfiguration <- function(nFile) {
  nPattern <- "_N([0-9]+)_M([0-9]+)_comb([0-9]+)\\.csv$"
  lMatch <- regmatches(
    basename(nFile),
    regexec(nPattern, basename(nFile))
  )[[1L]]
  
  if (length(lMatch) != 4L) {
    stop("Cannot identify N, M and combination from: ", nFile)
  }
  
  c(
    N = as.integer(lMatch[2L]),
    M = as.integer(lMatch[3L]),
    combination = as.integer(lMatch[4L])
  )
}

## Compute statistics for one configuration
fnConfigurationSummary <- function(df, nFile) {
  vConfiguration <- fnParseConfiguration(nFile)
  
  iCombination <- unname(vConfiguration["combination"])
  
  if (!(iCombination %in% dfTruth$combination)) {
    stop("Invalid parameter combination in: ", nFile)
  }
  
  vTruth <- as.numeric(
    dfTruth[
      dfTruth$combination == iCombination,
      vParameters,
      drop = TRUE
    ]
  )
  names(vTruth) <- vParameters
  
  dfResult <- data.frame(
    N = unname(vConfiguration["N"]),
    M = unname(vConfiguration["M"]),
    combination = iCombination,
    replications = nrow(df),
    median_polynomial_degree = median(
      df$polynomial_degree,
      na.rm = TRUE
    ),
    source_file = basename(nFile),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  for (nParameter in vParameters) {
    vEstimate <- as.numeric(df[[nParameter]])
    vEstimate <- vEstimate[is.finite(vEstimate)]
    
    dTruth <- vTruth[nParameter]
    
    dfResult[[paste0(nParameter, "_true")]] <- dTruth
    dfResult[[paste0(nParameter, "_bias")]] <-
      mean(vEstimate - dTruth)
    dfResult[[paste0(nParameter, "_sd")]] <-
      sd(vEstimate)
    dfResult[[paste0(nParameter, "_rmse")]] <-
      sqrt(mean((vEstimate - dTruth)^2))
  }
  
  dfResult
}

## Apply to every element of lEstimates
lSummary <- Map(
  fnConfigurationSummary,
  df = lEstimates,
  nFile = names(lEstimates)
)

## One row for each N-M-parameter configuration
dfSummary <- do.call(rbind, lSummary)
rownames(dfSummary) <- NULL

dfSummary <- dfSummary[
  order(
    dfSummary$N,
    dfSummary$M,
    dfSummary$combination
  ),
]

## Display results
print(dfSummary)

## Save results
write.csv(
  dfSummary,
  file = file.path(nPath, "Summary_BH_SNPMD_all_configurations.csv"),
  row.names = FALSE
)
