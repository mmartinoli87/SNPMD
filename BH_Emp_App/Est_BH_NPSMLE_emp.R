# Empirical NPSMLE estimation of the four-parameter Brock-Hommes
#
# Default design:
#   * S&P 500 and Euro Stoxx 50;
#   * nested 1-, 2-, and 4-year samples ending on 2012-12-31;
#   * N = M in {250, 500, 1000};
#   * 1,000 repetitions with new NPSMLE innovations in each repetition.
#
# Production example:
#   BH_REPLICATIONS=1000 BH_CORES=5 \
#   Rscript Est_BH_NPSMLE_emp.R

rm(list = ls())

args <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args || identical(Sys.getenv("BH_QUICK"), "1")
overwrite <- "--overwrite" %in% args || identical(Sys.getenv("BH_OVERWRITE"), "1")

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[1L]) else ""
script_file <- gsub("~+~", " ", script_file, fixed = TRUE)
default_root <- if (nzchar(script_file)) {
  normalizePath(file.path(dirname(script_file), "..."), mustWork = FALSE)
} else {
  path.expand(paste0(
    "~"
  ))
}
root_dir <- path.expand(Sys.getenv("BH_ROOT_DIR", unset = default_root))

BH_file <- file.path(root_dir, "BH_functions.R")
if (!file.exists(BH_file)) {
  BH_file <- file.path(root_dir, "BH_functions.R")
}
if (!file.exists(BH_file)) stop("Cannot find BH_functions.R")
source(BH_file)

env_integer <- function(name, default, minimum = 1L) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (length(value) != 1L || is.na(value) || value < minimum) {
    stop(name, " must be an integer >= ", minimum)
  }
  value
}

env_numeric <- function(name, default, lower = -Inf, strict = FALSE) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(default))))
  valid <- length(value) == 1L && is.finite(value)
  valid <- valid && if (strict) value > lower else value >= lower
  if (!valid) stop(name, if (strict) " must be > " else " must be >= ", lower)
  value
}

env_flag <- function(name, default = TRUE) {
  value <- tolower(trimws(Sys.getenv(name, unset = if (default) "1" else "0")))
  if (!value %in% c("0", "1", "false", "true", "no", "yes")) {
    stop(name, " must be one of 0, 1, false, true, no, or yes")
  }
  value %in% c("1", "true", "yes")
}

env_choice <- function(name, default, choices) {
  value <- tolower(trimws(Sys.getenv(name, unset = default)))
  if (!value %in% choices) {
    stop(name, " must be one of: ", paste(choices, collapse = ", "))
  }
  value
}

parse_integer_vector <- function(name, default) {
  specification <- Sys.getenv(name, unset = "")
  if (!nzchar(specification)) return(as.integer(default))
  tokens <- trimws(strsplit(specification, ",", fixed = TRUE)[[1L]])
  values <- suppressWarnings(as.integer(tokens))
  if (!length(values) || any(!nzchar(tokens)) || anyNA(values) || any(values < 1L)) {
    stop(name, " must contain comma-separated positive integers")
  }
  unique(values)
}

parse_numeric_vector <- function(name, default, expected_length) {
  text <- Sys.getenv(name, unset = default)
  values <- suppressWarnings(as.numeric(strsplit(text, ",", fixed = TRUE)[[1L]]))
  if (length(values) != expected_length || any(!is.finite(values))) {
    stop(name, " must contain ", expected_length, " comma-separated numbers")
  }
  values
}

# -----------------------------------------------------------------------------
# Design and paths
# -----------------------------------------------------------------------------

param_names <- c("Beta", "Sigma", "b_2", "g_2")
mc_replications <- env_integer("BH_REPLICATIONS", if (quick) 2L else 1000L)
years <- parse_integer_vector("BH_YEARS", c(1L, 2L, 4L))
days_per_year <- env_integer("BH_DAYS_PER_YEAR", if (quick) 50L else 250L)
sample_sizes <- years * days_per_year
names(sample_sizes) <- as.character(years)
validation_size <- env_integer(
  "BH_VALIDATION_SIZE", if (quick) 50L else days_per_year, minimum = 0L
)
validation_draws <- env_integer("BH_VALIDATION_DRAWS", if (quick) 20L else 1000L)
estimation_end <- as.Date(Sys.getenv("BH_ESTIMATION_END", unset = "2012-12-31"))
if (is.na(estimation_end)) stop("BH_ESTIMATION_END must have format YYYY-MM-DD")
half_window <- env_integer("BH_HALF_WINDOW", 30L)
fundamental_mode <- env_choice(
  "BH_FUNDAMENTAL_MODE", "trailing", c("centered", "trailing")
)
standardize_data <- env_flag("BH_STANDARDIZE_DATA", TRUE)
auto_download <- env_flag("BH_AUTO_DOWNLOAD", TRUE)
make_plots <- env_flag("BH_MAKE_PLOTS", TRUE)

detected_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
if (!is.finite(detected_cores)) detected_cores <- 2L
cores <- env_integer("BH_CORES", max(1L, detected_cores - 3L))
batch_size <- env_integer("BH_BATCH_SIZE", cores)
base_seed <- env_integer("BH_SEED", 400L, minimum = 0L)
maxit <- env_integer("BH_MAXIT", if (quick) 3L else 100L)
tolerance <- env_numeric("BH_TOLERANCE", if (quick) 1e-2 else 1e-3, 0, strict = TRUE)
risk_free_rate <- env_numeric("BH_RISK_FREE_RATE", 0.0001, -1, strict = TRUE)
burn_in <- env_integer("BH_BURN_IN", if (quick) 100L else 4000L)

parameter_bounds <- rbind(
  Beta  = c(lower = 0.0,  upper = 10.0),
  Sigma = c(lower = 0.0,  upper = 1.0),
  b_2   = c(lower = -1, upper = 1),
  g_2   = c(lower = -2, upper = 2)
)
initial_params <- parse_numeric_vector("BH_INITIAL_PARAMS", "3,1,0,0", 4L)
names(initial_params) <- param_names
if (any(initial_params < parameter_bounds[, "lower"] |
        initial_params > parameter_bounds[, "upper"])) {
  stop("BH_INITIAL_PARAMS must lie inside the aligned parameter bounds")
}

output_root <- path.expand(Sys.getenv(
  "BH_OUTPUT_DIR",
  unset = file.path(root_dir, "output_empirical_parallel_npsmle")
))
run_tag <- paste0(
  "Emp_est_NPSMLE",
  if (standardize_data) "_std" else "_unscaled"
)
experiment_dir <- file.path(output_root, run_tag)
dir.create(experiment_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Market data and empirical windows
# -----------------------------------------------------------------------------

market_specification <- data.frame(
  market = c("SP500", "EURO_STOXX_50"),
  label = c("S&P 500", "Euro Stoxx 50"),
  yahoo_symbol = c("^GSPC", "^STOXX50E"),
  default_file = c(
    file.path(root_dir, "data", "SP500.csv"),
    file.path(root_dir, "data", "EURO_STOXX_50.csv")
  ),
  stringsAsFactors = FALSE
)
market_specification$input_file <- c(
  path.expand(Sys.getenv("BH_SP500_FILE", unset = market_specification$default_file[1L])),
  path.expand(Sys.getenv(
    "BH_EURO_STOXX50_FILE", unset = market_specification$default_file[2L]
  ))
)

market_selection <- trimws(strsplit(
  Sys.getenv("BH_MARKETS", unset = "SP500,EURO_STOXX_50"), ",", fixed = TRUE
)[[1L]])
market_selection <- toupper(gsub("[- ]", "_", market_selection))
market_selection[market_selection %in% c("EU50", "EUROSTOXX50", "STOXX50E")] <-
  "EURO_STOXX_50"
if (any(!market_selection %in% market_specification$market)) {
  stop("BH_MARKETS must contain SP500 and/or EURO_STOXX_50")
}
market_specification <- market_specification[
  match(unique(market_selection), market_specification$market), , drop = FALSE
]

parse_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  x <- trimws(as.character(x))
  formats <- c("%Y-%m-%d", "%d/%m/%Y", "%m/%d/%Y", "%Y/%m/%d")
  candidates <- lapply(formats, function(fmt) as.Date(x, format = fmt))
  counts <- vapply(candidates, function(z) sum(!is.na(z)), integer(1L))
  candidates[[which.max(counts)]]
}

read_price_csv <- function(file, market) {
  if (!file.exists(file)) stop("File does not exist: ", file)
  raw <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(raw)) stop("No observations in ", file)
  clean_names <- gsub("[^a-z0-9]", "", tolower(names(raw)))
  date_column <- which(clean_names %in% c("date", "datetime", "timestamp"))[1L]
  if (is.na(date_column)) stop("The file for ", market, " needs a Date column")
  dates <- parse_date(raw[[date_column]])
  if (sum(!is.na(dates)) < 0.9 * nrow(raw)) stop("Could not parse dates in ", file)

  preferred <- if (market == "SP500") {
    c("sp500", "adjusted", "adjclose", "close", "price", "value")
  } else {
    c("eurostoxx50", "stoxx50e", "sx5e", "adjusted", "adjclose", "close", "price", "value")
  }
  price_column <- match(preferred, clean_names, nomatch = 0L)
  price_column <- price_column[price_column > 0L & price_column != date_column]
  if (!length(price_column)) {
    numeric_columns <- which(vapply(raw, function(z) {
      converted <- suppressWarnings(as.numeric(as.character(z)))
      mean(is.finite(converted)) > 0.9
    }, logical(1L)))
    numeric_columns <- setdiff(numeric_columns, date_column)
    if (!length(numeric_columns)) stop("Could not identify a price column in ", file)
    price_column <- numeric_columns[1L]
  } else {
    price_column <- price_column[1L]
  }

  prices <- suppressWarnings(as.numeric(as.character(raw[[price_column]])))
  data <- data.frame(Date = dates, Price = prices)
  data <- data[is.finite(data$Price) & !is.na(data$Date), , drop = FALSE]
  data <- data[order(data$Date), , drop = FALSE]
  data <- data[!duplicated(data$Date, fromLast = TRUE), , drop = FALSE]
  rownames(data) <- NULL
  if (nrow(data) < 3L) stop("Too few usable observations in ", file)
  data
}

download_prices <- function(symbol, market) {
  if (!requireNamespace("quantmod", quietly = TRUE)) {
    stop("Install 'quantmod' or provide a sufficiently long dated CSV for ", market)
  }
  past_days <- ceiling((max(sample_sizes) + 2L * half_window + 50L) * 365.25 / 250)
  future_days <- ceiling((validation_size + 2L * half_window + 50L) * 365.25 / 250)
  from <- estimation_end - past_days
  to <- estimation_end + future_days
  cat(sprintf("Downloading %s (%s) from %s to %s ...\n", market, symbol, from, to))
  x <- suppressWarnings(quantmod::getSymbols(
    symbol, src = "yahoo", from = from, to = to + 1L,
    auto.assign = FALSE, warnings = FALSE
  ))
  p <- tryCatch(quantmod::Ad(x), error = function(e) NULL)
  if (is.null(p) || all(!is.finite(as.numeric(p)))) p <- quantmod::Cl(x)
  data <- data.frame(Date = as.Date(zoo::index(p)), Price = as.numeric(p))
  data <- data[is.finite(data$Price) & !is.na(data$Date), , drop = FALSE]
  if (!nrow(data)) stop("No usable prices downloaded for ", market)
  cache_file <- file.path(experiment_dir, paste0(market, "_downloaded_prices.csv"))
  write.csv(data, cache_file, row.names = FALSE)
  data
}

construct_deviations <- function(price_data) {
  y <- price_data$Price
  n <- length(y)
  w <- half_window
  if (fundamental_mode == "centered") {
    if (n <= 2L * w) stop("Price series is too short for the centred MA")
    centres <- seq.int(w + 1L, n - w)
    fundamental <- vapply(
      centres, function(i) mean(y[(i - w):(i + w)]), numeric(1L)
    )
  } else {
    window <- 2L * w + 1L
    if (n < window) stop("Price series is too short for the trailing MA")
    centres <- seq.int(window, n)
    fundamental <- vapply(
      centres, function(i) mean(y[(i - window + 1L):i]), numeric(1L)
    )
  }
  data.frame(
    Date = price_data$Date[centres], observed = y[centres],
    fundamental = fundamental, deviation = y[centres] - fundamental
  )
}

extract_windows <- function(deviations, market) {
  pre <- deviations[deviations$Date <= estimation_end, , drop = FALSE]
  if (nrow(pre) < max(sample_sizes)) {
    stop(market, " has ", nrow(pre), " pre-estimation observations; ",
         max(sample_sizes), " are required")
  }
  training <- setNames(lapply(sample_sizes, function(n) tail(pre, n)), names(sample_sizes))
  if (validation_size > 0L) {
    post <- deviations[deviations$Date > estimation_end, , drop = FALSE]
    if (nrow(post) < validation_size) {
      stop(market, " has ", nrow(post), " holdout observations; ",
           validation_size, " are required")
    }
    validation <- head(post, validation_size)
  } else {
    validation <- NULL
  }
  list(training = training, validation = validation)
}

load_market <- function(i) {
  specification <- market_specification[i, ]
  cache_file <- file.path(
    experiment_dir, paste0(specification$market, "_downloaded_prices.csv")
  )
  prepare <- function(prices) {
    deviations <- construct_deviations(prices)
    list(
      prices = prices, deviations = deviations,
      windows = extract_windows(deviations, specification$market)
    )
  }
  attempt <- function(file) tryCatch(
    prepare(read_price_csv(file, specification$market)),
    error = function(e) structure(list(message = conditionMessage(e)), class = "data_error")
  )

  prepared <- attempt(specification$input_file)
  source_file <- specification$input_file
  if (inherits(prepared, "data_error") && file.exists(cache_file)) {
    cached <- attempt(cache_file)
    if (!inherits(cached, "data_error")) {
      prepared <- cached
      source_file <- cache_file
    }
  }
  if (inherits(prepared, "data_error")) {
    if (!auto_download) stop(prepared$message)
    warning("Local data for ", specification$market, " are unsuitable (",
            prepared$message, "); downloading a dated Yahoo Finance series.")
    prices <- download_prices(specification$yahoo_symbol, specification$market)
    prepared <- prepare(prices)
    source_file <- cache_file
  }
  prepared$market <- specification$market
  prepared$label <- specification$label
  prepared$source_file <- source_file
  prepared
}

markets <- lapply(seq_len(nrow(market_specification)), load_market)
names(markets) <- vapply(markets, `[[`, character(1L), "market")

jobs <- expand.grid(
  market = names(markets), years = years,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
jobs$N <- unname(sample_sizes[as.character(jobs$years)])
jobs$M <- jobs$N

# Training-scale normalization preserves every sign statistic used by SNPMD
# while making the [0,1] Sigma bound meaningful for index-point deviations.
job_data <- lapply(seq_len(nrow(jobs)), function(i) {
  market <- markets[[jobs$market[i]]]
  training <- market$windows$training[[as.character(jobs$years[i])]]
  scale <- if (standardize_data) population_sd(training$deviation) else 1
  if (!is.finite(scale) || scale <= 0) stop("Invalid scale for empirical job ", i)
  list(
    training = training,
    y = training$deviation / scale,
    validation = if (is.null(market$windows$validation)) NULL else
      market$windows$validation$deviation / scale,
    scale = scale
  )
})

# -----------------------------------------------------------------------------
# NPSMLE
# -----------------------------------------------------------------------------

make_seed <- function(job, replication, stream) {
  modulus <- 2147483646
  value <- (
    as.double(base_seed) + as.double(job) * 1000003 +
      as.double(replication) * 101 + as.double(stream) * 15485863
  ) %% modulus
  as.integer(value + 1)
}

empty_fit <- function(message, elapsed = NA_real_) {
  data.frame(
    Beta_hat = NA_real_, Sigma_hat = NA_real_, b_2_hat = NA_real_, g_2_hat = NA_real_,
    objective = NA_real_, convergence = NA_integer_,
    objective_evaluations = NA_integer_, gradient_evaluations = NA_integer_,
    elapsed_seconds = elapsed, status = "error",
    message = substr(as.character(message), 1L, 500L),
    stringsAsFactors = FALSE, check.names = FALSE
  )
}

kernel_log_likelihood <- function(y, simulator_step, theta, epsilon, num_lags = 3L) {
  y <- as.numeric(y)
  epsilon <- as.numeric(epsilon)
  num_draws <- length(epsilon)
  if (num_draws < 2L) stop("At least two likelihood draws are required")
  if (length(y) <= num_lags) stop("The empirical series is too short")
  contributions <- numeric(length(y) - num_lags)
  for (i in seq.int(num_lags + 1L, length(y))) {
    simulated <- as.numeric(simulator_step(
      y[(i - num_lags):(i - 1L)], theta, epsilon
    ))
    if (length(simulated) != num_draws || any(!is.finite(simulated))) {
      return(rep(-Inf, length(contributions)))
    }
    difference <- simulated - y[i]
    bandwidth <- population_sd(difference) * (4 / (3 * num_draws))^(1 / 5)
    if (!is.finite(bandwidth) || bandwidth <= 0) {
      return(rep(-Inf, length(contributions)))
    }
    log_kernel <- dnorm(difference, 0, bandwidth, log = TRUE)
    maximum <- max(log_kernel)
    contributions[i - num_lags] <- maximum + log(mean(exp(log_kernel - maximum)))
  }
  contributions
}

fit_npsmle <- function(y, sample_size, job, replication) {
  start <- proc.time()[[3L]]
  tryCatch({
    set.seed(make_seed(job, replication, 2L))
    # M=N draws; fixed within the optimizer run and renewed across repetitions.
    epsilon <- qnorm(runif(sample_size))
    simulator <- brock_hommes(seed = 3L, r = risk_free_rate)
    objective <- function(a) {
      contributions <- kernel_log_likelihood(y, simulator$step, a, epsilon, 3L)
      value <- -mean(contributions)
      if (is.finite(value)) value else 1e100
    }
    # fit <- optim(
    #   par = initial_params, fn = objective, method = "L-BFGS-B",
    #   lower = parameter_bounds[, "lower"], upper = parameter_bounds[, "upper"],
    #   control = list(maxit = maxit, factr = 1e7, pgtol = tolerance)
    # )
    fit <- optim(
      par = initial_params, fn = objective, method = "BFGS",
      control = list(maxit = maxit, factr = 1e7, pgtol = tolerance)
    )
    estimates <- as.numeric(fit$par)
    data.frame(
      Beta_hat = estimates[1L], Sigma_hat = estimates[2L],
      b_2_hat = estimates[3L], g_2_hat = estimates[4L],
      objective = as.numeric(fit$value), convergence = as.integer(fit$convergence),
      objective_evaluations = unname(as.integer(fit$counts["function"])),
      gradient_evaluations = unname(as.integer(fit$counts["gradient"])),
      elapsed_seconds = proc.time()[[3L]] - start,
      status = if (all(is.finite(estimates))) "ok" else "nonfinite",
      message = if (is.null(fit$message)) "" else as.character(fit$message),
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }, error = function(e) empty_fit(conditionMessage(e), proc.time()[[3L]] - start))
}

run_replication <- function(replication, job) {
  specification <- jobs[job, ]
  market <- markets[[specification$market]]
  data <- job_data[[job]]
  fit <- fit_npsmle(data$y, specification$N, job, replication)
  cbind(
    data.frame(
      replication = as.integer(replication), market = specification$market,
      market_label = market$label, years = specification$years,
      training_start = min(data$training$Date), training_end = max(data$training$Date),
      N = specification$N, M = specification$M, estimator = "NPSMLE",
      stringsAsFactors = FALSE, check.names = FALSE
    ),
    fit,
    data.frame(
      data_scale = data$scale, standardized = standardize_data,
      npsmle_draws = specification$M,
      empirical_sample = "fixed_within_market_window",
      likelihood_aggregation = "log_mean_kernel_density",
      fundamental_mode = fundamental_mode,
      estimation_end = as.character(estimation_end),
      risk_free_rate = risk_free_rate, maxit = maxit, tolerance = tolerance,
      stringsAsFactors = FALSE, check.names = FALSE
    )
  )
}

write_csv_atomic <- function(x, path) {
  temporary <- tempfile(pattern = "empirical_", tmpdir = dirname(path), fileext = ".csv")
  write.csv(x, temporary, row.names = FALSE)
  if (!file.rename(temporary, path)) stop("Could not write ", path)
  invisible(path)
}

completed_replications <- function(results) {
  if (is.null(results) || !nrow(results)) return(integer())
  as.integer(unique(results$replication[
    !is.na(results$replication) & results$replication >= 1L
  ]))
}

combine_without_duplicates <- function(old, new) {
  if (is.null(old) || !nrow(old)) return(new)
  both <- rbind(old, new)
  both <- both[!duplicated(both$replication, fromLast = TRUE), , drop = FALSE]
  both[order(both$replication), , drop = FALSE]
}

parallel_apply <- function(values, job, cluster = NULL) {
  if (cores <= 1L) return(lapply(values, run_replication, job = job))
  if (.Platform$OS.type != "windows") {
    return(parallel::mclapply(
      values, run_replication, job = job,
      mc.cores = cores, mc.preschedule = FALSE, mc.set.seed = FALSE
    ))
  }
  parallel::parLapply(cluster, values, run_replication, job = job)
}

cat("\nEmpirical Brock-Hommes NPSMLE experiment\n")
cat("  markets:", paste(names(markets), collapse = ", "), "\n")
cat("  years:", paste(years, collapse = ", "), "\n")
cat("  sample sizes (N=M):", paste(sample_sizes, collapse = ", "), "\n")
cat("  repetitions per cell:", mc_replications, "\n")
cat("  standardized by training SD:", standardize_data, "\n")
cat("  workers:", cores, "\n")
cat("  output:", experiment_dir, "\n\n")

cluster <- NULL
if (cores > 1L && .Platform$OS.type == "windows") {
  cluster <- parallel::makePSOCKcluster(cores)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterExport(cluster, c("BH_file"), envir = .GlobalEnv)
  parallel::clusterEvalQ(cluster, source(BH_file))
  parallel::clusterExport(
    cluster,
    c(
      "jobs", "markets", "job_data", "standardize_data", "fundamental_mode",
      "estimation_end", "risk_free_rate", "maxit", "tolerance", "base_seed",
      "initial_params", "parameter_bounds", "make_seed", "empty_fit",
      "kernel_log_likelihood", "fit_npsmle", "run_replication"
    ),
    envir = .GlobalEnv
  )
}

for (job in seq_len(nrow(jobs))) {
  specification <- jobs[job, ]
  cell_file <- file.path(
    experiment_dir,
    sprintf(
      "cell_%s_Y%d_N%d_M%d.csv", specification$market,
      specification$years, specification$N, specification$M
    )
  )
  existing <- NULL
  if (file.exists(cell_file) && !overwrite) {
    existing <- read.csv(cell_file, stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(existing)) {
      compatible <-
        all(existing$N == specification$N) && all(existing$M == specification$M) &&
        all(existing$market == specification$market) &&
        all(existing$fundamental_mode == fundamental_mode) &&
        all(existing$standardized == standardize_data) &&
        all(abs(existing$risk_free_rate - risk_free_rate) < .Machine$double.eps^0.5)
      if (!compatible) stop("Existing checkpoint is incompatible: ", cell_file)
    }
  }
  done <- completed_replications(existing)
  pending <- setdiff(seq_len(mc_replications), done)
  cat(sprintf(
    "%s, %d year(s), N=M=%d: %d completed, %d pending\n",
    specification$market, specification$years, specification$N,
    length(done), length(pending)
  ))
  if (length(pending)) {
    batches <- split(pending, ceiling(seq_along(pending) / batch_size))
    for (batch in batches) {
      start <- proc.time()[[3L]]
      batch_results <- do.call(rbind, parallel_apply(batch, job, cluster))
      existing <- combine_without_duplicates(existing, batch_results)
      write_csv_atomic(existing, cell_file)
      cat(sprintf(
        "  saved through repetition %d (%d/%d; %.1f seconds)\n",
        max(batch), length(completed_replications(existing)), mc_replications,
        proc.time()[[3L]] - start
      ))
    }
  }
}

# -----------------------------------------------------------------------------
# Monte Carlo summaries
# -----------------------------------------------------------------------------

cell_files <- list.files(
  experiment_dir,
  pattern = "^cell_(SP500|EURO_STOXX_50)_Y[0-9]+_N[0-9]+_M[0-9]+[.]csv$",
  full.names = TRUE
)
all_results <- do.call(rbind, lapply(
  cell_files, read.csv, stringsAsFactors = FALSE, check.names = FALSE
))
all_results <- all_results[
  order(all_results$market, all_results$years, all_results$replication), , drop = FALSE
]
rownames(all_results) <- NULL

estimate_names <- paste0(param_names, "_hat")
groups <- split(all_results, interaction(all_results$market, all_results$years, drop = TRUE))
summary_long <- do.call(rbind, lapply(groups, function(block) {
  complete_fit <- block$status == "ok" &
    rowSums(is.finite(as.matrix(block[, estimate_names, drop = FALSE]))) ==
      length(estimate_names)
  convergence_rate <- mean(complete_fit & block$convergence == 0L, na.rm = TRUE)
  do.call(rbind, lapply(seq_along(param_names), function(j) {
    valid <- is.finite(block[[estimate_names[j]]])
    estimate <- block[[estimate_names[j]]][valid]
    data.frame(
      market = block$market[1L], market_label = block$market_label[1L],
      years = block$years[1L], N = block$N[1L], M = block$M[1L],
      parameter = param_names[j], replications_requested = mc_replications,
      replications_completed = length(unique(block$replication)),
      n_finite = length(estimate), mean = mean(estimate),
      sd = if (length(estimate) > 1L) sd(estimate) else NA_real_,
      mcse = if (length(estimate) > 1L) sd(estimate) / sqrt(length(estimate)) else NA_real_,
      median = median(estimate), q025 = unname(quantile(estimate, 0.025)),
      q975 = unname(quantile(estimate, 0.975)),
      convergence_rate = convergence_rate, failure_rate = mean(!complete_fit),
      median_objective_evaluations = median(block$objective_evaluations, na.rm = TRUE),
      median_elapsed_seconds = median(block$elapsed_seconds, na.rm = TRUE),
      row.names = NULL, check.names = FALSE
    )
  }))
}))
summary_long <- summary_long[
  order(summary_long$market, summary_long$years,
        match(summary_long$parameter, param_names)), , drop = FALSE
]
rownames(summary_long) <- NULL

summary_wide <- do.call(rbind, lapply(groups, function(block) {
  result <- data.frame(
    market = block$market[1L], market_label = block$market_label[1L],
    years = block$years[1L], N = block$N[1L], M = block$M[1L],
    replications = nrow(block), row.names = NULL, check.names = FALSE
  )
  for (parameter in param_names) {
    estimate <- block[[paste0(parameter, "_hat")]]
    estimate <- estimate[is.finite(estimate)]
    result[[paste0("mean_", parameter)]] <- mean(estimate)
    result[[paste0("sd_", parameter)]] <- sd(estimate)
    result[[paste0("mcse_", parameter)]] <- sd(estimate) / sqrt(length(estimate))
  }
  result$median_objective <- median(block$objective, na.rm = TRUE)
  result
}))
summary_wide <- summary_wide[order(summary_wide$market, summary_wide$years), ]
rownames(summary_wide) <- NULL

mean_estimates <- do.call(rbind, lapply(groups, function(block) {
  result <- block[1L, c(
    "market", "market_label", "years", "training_start", "training_end", "N", "M",
    "data_scale"
  ), drop = FALSE]
  for (parameter in param_names) {
    estimate <- block[[paste0(parameter, "_hat")]]
    result[[parameter]] <- mean(estimate[is.finite(estimate)])
  }
  result
}))
mean_estimates <- mean_estimates[order(mean_estimates$market, mean_estimates$years), ]
rownames(mean_estimates) <- NULL

# # -----------------------------------------------------------------------------
# # Held-out distributional and likelihood validation at the mean estimate
# # -----------------------------------------------------------------------------
# 
# process_probabilities <- function(x) {
#   sign <- as.integer(x > 0)
#   n <- length(sign)
#   combinations <- as.matrix(expand.grid(previous = 0:1, current = 0:1))
#   transitions <- apply(combinations, 1L, function(state) {
#     sum(sign[2L:(n - 1L)] == state[1L] & sign[3L:n] == state[2L]) / n
#   })
#   setNames(
#     c(mean(sign), transitions),
#     c("Pr(Y>0)", "Pr(0,0)", "Pr(1,0)", "Pr(0,1)", "Pr(1,1)")
#   )
# }
# 
# gray_distance <- function(target, candidate) {
#   if (is.null(dim(candidate))) candidate <- matrix(candidate, nrow = 1L)
#   2 * abs(target[1L] - candidate[, 1L]) +
#     rowSums(abs(sweep(candidate[, -1L, drop = FALSE], 2L, target[-1L], "-")))
# }
# 
# validate_job <- function(job) {
#   specification <- jobs[job, ]
#   data <- job_data[[job]]
#   if (is.null(data$validation)) return(NULL)
#   point <- mean_estimates[
#     mean_estimates$market == specification$market &
#       mean_estimates$years == specification$years, , drop = FALSE
#   ]
#   theta <- as.numeric(point[1L, param_names])
#   target <- process_probabilities(data$validation)
#   simulated_probabilities <- do.call(rbind, lapply(seq_len(validation_draws), function(b) {
#     simulator <- brock_hommes(
#       seed = make_seed(job, b, 8L), r = risk_free_rate
#     )
#     simulated <- simulator$simulate(burn_in + validation_size + 3L, theta)
#     process_probabilities(tail(simulated, validation_size))
#   }))
#   model_mean <- colMeans(simulated_probabilities)
#   intervals <- apply(
#     simulated_probabilities, 2L, quantile,
#     probs = c(0.025, 0.975), names = FALSE
#   )
#   distances <- gray_distance(target, simulated_probabilities)
# 
#   # One-step-ahead holdout log score: the last three training observations are
#   # supplied as initial lags, and the likelihood uses M=validation_size draws.
#   set.seed(make_seed(job, 1L, 9L))
#   epsilon <- qnorm(runif(validation_size))
#   simulator <- brock_hommes(seed = 3L, r = risk_free_rate)
#   evaluation_series <- c(tail(data$y, 3L), data$validation)
#   holdout_contributions <- kernel_log_likelihood(
#     evaluation_series, simulator$step, theta, epsilon, 3L
#   )
# 
#   list(
#     empirical = target, model_mean = model_mean,
#     lower = intervals[1L, ], upper = intervals[2L, ],
#     distance_to_model_mean = gray_distance(target, model_mean),
#     distance_mean = mean(distances), distance_median = median(distances),
#     distance_q025 = unname(quantile(distances, 0.025)),
#     distance_q975 = unname(quantile(distances, 0.975)),
#     mean_holdout_log_score = mean(holdout_contributions),
#     mean_holdout_negative_log_score = -mean(holdout_contributions)
#   )
# }
# 
# validation_results <- lapply(seq_len(nrow(jobs)), validate_job)
# probability_validation <- do.call(rbind, lapply(seq_len(nrow(jobs)), function(i) {
#   result <- validation_results[[i]]
#   if (is.null(result)) return(NULL)
#   data.frame(
#     market = jobs$market[i], years = jobs$years[i], N = jobs$N[i], M = jobs$M[i],
#     statistic = names(result$empirical), empirical = as.numeric(result$empirical),
#     simulated_mean = as.numeric(result$model_mean),
#     simulated_q025 = as.numeric(result$lower), simulated_q975 = as.numeric(result$upper),
#     empirical_inside_95_interval =
#       result$empirical >= result$lower & result$empirical <= result$upper,
#     row.names = NULL, check.names = FALSE
#   )
# }))
# validation_summary <- do.call(rbind, lapply(seq_len(nrow(jobs)), function(i) {
#   result <- validation_results[[i]]
#   if (is.null(result)) return(NULL)
#   market <- markets[[jobs$market[i]]]
#   data.frame(
#     market = jobs$market[i], years = jobs$years[i], N = jobs$N[i], M = jobs$M[i],
#     validation_start = min(market$windows$validation$Date),
#     validation_end = max(market$windows$validation$Date),
#     validation_N = validation_size, validation_M = validation_size,
#     distance_to_simulated_mean = result$distance_to_model_mean,
#     mean_distance_across_draws = result$distance_mean,
#     median_distance_across_draws = result$distance_median,
#     distance_q025 = result$distance_q025, distance_q975 = result$distance_q975,
#     mean_holdout_log_score = result$mean_holdout_log_score,
#     mean_holdout_negative_log_score = result$mean_holdout_negative_log_score,
#     statistics_inside_95_interval = sum(
#       result$empirical >= result$lower & result$empirical <= result$upper
#     ),
#     total_statistics = length(result$empirical), row.names = NULL,
#     check.names = FALSE
#   )
# }))
# 
# window_summary <- do.call(rbind, lapply(seq_len(nrow(jobs)), function(i) {
#   data <- job_data[[i]]
#   data.frame(
#     market = jobs$market[i], years = jobs$years[i],
#     start = min(data$training$Date), end = max(data$training$Date),
#     observations = nrow(data$training), data_scale = data$scale,
#     sd_deviation = sd(data$training$deviation),
#     share_positive = mean(data$training$deviation > 0),
#     row.names = NULL, check.names = FALSE
#   )
# }))

for (market in names(markets)) {
  write.csv(
    markets[[market]]$deviations,
    file.path(experiment_dir, paste0("BH_NPSMLE_", market, "_analysis_data.csv")),
    row.names = FALSE
  )
}
write_csv_atomic(all_results, file.path(experiment_dir, "BH_NPSMLE_empirical_estimates.csv"))
write_csv_atomic(mean_estimates, file.path(experiment_dir, "BH_NPSMLE_empirical_mean_estimates.csv"))
write_csv_atomic(summary_long, file.path(experiment_dir, "BH_NPSMLE_empirical_mc_standard_errors.csv"))
write_csv_atomic(summary_wide, file.path(experiment_dir, "BH_NPSMLE_empirical_mc_summary_wide.csv"))

saveRDS(
  list(
    estimates = all_results, mean_estimates = mean_estimates,
    mc_summary = summary_long, mc_summary_wide = summary_wide,
    bounds = parameter_bounds,
    settings = list(
      replications = mc_replications, years = years, days_per_year = days_per_year,
      sample_sizes = sample_sizes, estimation_end = estimation_end,
      validation_size = validation_size, validation_draws = validation_draws,
      fundamental_mode = fundamental_mode, standardized = standardize_data,
      risk_free_rate = risk_free_rate, burn_in = burn_in,
      maxit = maxit, tolerance = tolerance, initial_params = initial_params,
      cores = cores
    )
  ),
  file.path(experiment_dir, "BH_NPSMLE_empirical_results.rds")
)

cat("\nOutputs written to:", experiment_dir, "\n")
cat("Mean estimates:", file.path(experiment_dir, "BH_NPSMLE_empirical_mean_estimates.csv"), "\n")
cat("Monte Carlo standard errors:",
    file.path(experiment_dir, "BH_NPSMLE_empirical_mc_standard_errors.csv"), "\n")
