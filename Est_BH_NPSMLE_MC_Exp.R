
# Monte Carlo NPSMLE experiment for the four-parameter Brock-Hommes model.
#
# Alignment choices:
#   * the same 16 truth configurations and N = M in {250, 500, 1000};
#   * the same r = 0.01 and 5,000-period simulation horizon;
#   * one fixed empirical realization (seed 3) per truth configuration;
#   * M = N one-step innovations in the simulated likelihood;
#   * the same bounded parameter space as SNPMD; and
#   * the conventional log of the average Gaussian-kernel density.

rm(list = ls())

args <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args
overwrite <- "--overwrite" %in% args || identical(Sys.getenv("BH_OVERWRITE"), "1")

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[1L]) else ""
default_root <- if (nzchar(script_file)) {
  normalizePath(file.path(dirname(script_file), ".."), mustWork = FALSE)
} else {
  path.expand(paste0(
    "~"
  ))
}
root_dir <- path.expand(Sys.getenv("BH_ROOT_DIR", unset = default_root))

source_dir <- file.path(root_dir, "R")
if (!dir.exists(source_dir)) source_dir <- root_dir
required_files <- file.path(
  source_dir, c("BH_functions.R", "data_helpers.R")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Required source file(s) not found: ", paste(missing_files, collapse = ", "))
}
invisible(lapply(required_files, source))

env_integer <- function(name, default, minimum = 1L) {
  value <- Sys.getenv(name, unset = as.character(default))
  value <- suppressWarnings(as.integer(value))
  if (length(value) != 1L || is.na(value) || value < minimum) {
    stop(name, " must be an integer >= ", minimum)
  }
  value
}

env_numeric <- function(name, default, lower = -Inf, strict = FALSE) {
  value <- Sys.getenv(name, unset = as.character(default))
  value <- suppressWarnings(as.numeric(value))
  valid <- length(value) == 1L && is.finite(value)
  valid <- valid && if (strict) value > lower else value >= lower
  if (!valid) stop(name, if (strict) " must be > " else " must be >= ", lower)
  value
}

parse_integer_set <- function(text, valid = NULL, name = "value") {
  text <- trimws(text)
  if (!nzchar(text) || tolower(text) == "all") {
    if (is.null(valid)) stop(name, " cannot be empty")
    return(as.integer(valid))
  }
  pieces <- strsplit(text, ",", fixed = TRUE)[[1L]]
  out <- integer()
  for (piece in pieces) {
    piece <- trimws(piece)
    if (grepl("^[0-9]+-[0-9]+$", piece)) {
      endpoints <- as.integer(strsplit(piece, "-", fixed = TRUE)[[1L]])
      out <- c(out, seq.int(endpoints[1L], endpoints[2L]))
    } else {
      number <- suppressWarnings(as.integer(piece))
      if (is.na(number)) stop("Invalid ", name, ": ", piece)
      out <- c(out, number)
    }
  }
  out <- unique(out)
  if (!is.null(valid) && any(!out %in% valid)) {
    stop(name, " must be drawn from: ", paste(valid, collapse = ", "))
  }
  out
}

parse_numeric_vector <- function(text, expected_length, name) {
  out <- suppressWarnings(as.numeric(strsplit(text, ",", fixed = TRUE)[[1L]]))
  if (length(out) != expected_length || any(!is.finite(out))) {
    stop(name, " must contain ", expected_length, " comma-separated numbers")
  }
  out
}

# Parameters configurations
param_names <- c("Beta", "Sigma", "b_2", "g_2")

truth_base <- expand.grid(
  Beta = c(0, 0.5, 3, 10),
  g_2 = c(-0.4, 0.4),
  b_2 = c(-0.3, 0.3),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
truth_grid <- data.frame(
  configuration = seq_len(nrow(truth_base)),
  Beta = truth_base$Beta,
  Sigma = 1,
  b_2 = truth_base$b_2,
  g_2 = truth_base$g_2,
  check.names = FALSE
)

default_replications <- if (quick) 2L else 1000L
default_sizes <- if (quick) "250" else "250,500,1000"
default_maxit <- if (quick) 3L else 100L

mc_replications <- env_integer("BH_REPLICATIONS", default_replications)
sample_sizes <- parse_integer_set(
  Sys.getenv("BH_SAMPLE_SIZES", unset = default_sizes),
  name = "BH_SAMPLE_SIZES"
)
combination_ids <- parse_integer_set(
  Sys.getenv("BH_COMBINATIONS", unset = "all"),
  valid = truth_grid$configuration,
  name = "BH_COMBINATIONS"
)

estimators <- "NPSMLE"

detected_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
if (!is.finite(detected_cores)) detected_cores <- 2L
cores <- env_integer("BH_CORES", max(1L, detected_cores - 1L))
batch_size <- env_integer("BH_BATCH_SIZE", cores)
base_seed <- env_integer("BH_SEED", 400L, minimum = 0L)
maxit <- env_integer("BH_MAXIT", default_maxit)
tolerance <- env_numeric("BH_TOLERANCE", if (quick) 1e-2 else 1e-3, 0, strict = TRUE)
risk_free_rate <- env_numeric("BH_RISK_FREE_RATE", 0.01, 0)
time_horizon <- env_integer("BH_TIME_HORIZON", if (quick) 300L else 5000L)
if (any(sample_sizes > time_horizon)) {
  stop("Every BH_SAMPLE_SIZES value must be <= BH_TIME_HORIZON")
}

parameter_bounds <- rbind(
  Beta  = c(lower = 0.0,  upper = 10.0),
  Sigma = c(lower = 0.0,  upper = 1.0),
  b_2   = c(lower = -0.3, upper = 0.3),
  g_2   = c(lower = -0.4, upper = 0.4)
)
initial_params <- parse_numeric_vector(
  Sys.getenv("BH_INITIAL_PARAMS", unset = "2,1,0,0"),
  length(param_names), "BH_INITIAL_PARAMS"
)
names(initial_params) <- param_names
if (any(initial_params < parameter_bounds[, "lower"] |
        initial_params > parameter_bounds[, "upper"])) {
  stop("BH_INITIAL_PARAMS must lie inside the common SNPMD parameter bounds")
}

output_root <- path.expand(Sys.getenv(
  "BH_OUTPUT_DIR", unset = file.path(root_dir, "output_parallel_npsmle")
))
estimator_tag <- paste(sort(estimators), collapse = "-")
run_tag <- paste0(
  "MC_R", mc_replications,
  "_M_equals_N",
  "_I", maxit,
  "_", estimator_tag,
  if (quick) "_quick" else ""
)
experiment_dir <- file.path(output_root, run_tag)
dir.create(experiment_dir, recursive = TRUE, showWarnings = FALSE)

# Deterministic seeds make results invariant to the number of workers and to
# whether a run is resumed or a subset of configurations is requested.
make_seed <- function(sample_size, configuration, replication, stream) {
  modulus <- 2147483646
  value <- (as.double(base_seed) +
    as.double(sample_size) * 1000003 +
    as.double(configuration) * 10007 +
    as.double(replication) * 101 +
    as.double(stream) * 15485863) %% modulus
  as.integer(value + 1)
}

# For each truth configuration, generate one benchmark series with seed 3 and
# reuse its final N observations in all Monte Carlo replications.
benchmark_data <- list()
for (configuration in combination_ids) {
  truth <- as.numeric(truth_grid[
    truth_grid$configuration == configuration, param_names, drop = TRUE
  ])
  benchmark_simulator <- brock_hommes(seed = 3L, r = risk_free_rate)
  benchmark_path <- benchmark_simulator$simulate(time_horizon + 3L, truth)
  for (sample_size in sample_sizes) {
    benchmark_data[[paste(configuration, sample_size, sep = "_")]] <-
      tail(benchmark_path, sample_size)
  }
}

empty_fit <- function(estimator, message, elapsed = NA_real_) {
  data.frame(
    estimator = estimator,
    Beta_hat = NA_real_, Sigma_hat = NA_real_,
    b_2_hat = NA_real_, g_2_hat = NA_real_,
    objective = NA_real_, convergence = NA_integer_,
    objective_evaluations = NA_integer_, gradient_evaluations = NA_integer_,
    elapsed_seconds = elapsed,
    status = "error", message = substr(as.character(message), 1L, 500L),
    stringsAsFactors = FALSE, check.names = FALSE
  )
}

# Gaussian-kernel simulated log likelihood. For observation t, the conditional
# density estimate is the arithmetic average of the M kernel ordinates, and the
# likelihood contribution is the log of that average. A log-sum-exp calculation
# avoids numerical underflow.
kernel_log_likelihood <- function(y, simulator_step, theta, epsilon,
                                  num_lags = 3L) {
  y <- as.numeric(y)
  epsilon <- as.numeric(epsilon)
  num_lags <- as.integer(num_lags)
  num_draws <- length(epsilon)
  if (num_draws < 2L) stop("At least two likelihood draws are required")
  if (length(y) <= num_lags) stop("The empirical series is too short")

  contributions <- numeric(length(y) - num_lags)
  for (i in seq.int(num_lags + 1L, length(y))) {
    y_simulated <- as.numeric(simulator_step(
      y[(i - num_lags):(i - 1L)], theta, epsilon
    ))
    if (length(y_simulated) != num_draws || any(!is.finite(y_simulated))) {
      return(rep(-Inf, length(contributions)))
    }
    difference <- y_simulated - y[i]
    bandwidth <- population_sd(difference) * (4 / (3 * num_draws))^(1 / 5)
    if (!is.finite(bandwidth) || bandwidth <= 0) {
      return(rep(-Inf, length(contributions)))
    }
    log_kernel <- dnorm(difference, mean = 0, sd = bandwidth, log = TRUE)
    maximum <- max(log_kernel)
    contributions[i - num_lags] <-
      maximum + log(mean(exp(log_kernel - maximum)))
  }
  contributions
}

fit_npsmle <- function(y, simulator, sample_size, configuration, replication) {
  start <- proc.time()[[3L]]
  tryCatch({
    set.seed(make_seed(sample_size, configuration, replication, 2L))
    epsilon <- qnorm(runif(sample_size))
    objective <- function(a) {
      contributions <- kernel_log_likelihood(
        y, simulator$step, a, epsilon, num_lags = 3L
      )
      value <- -mean(contributions)
      if (is.finite(value)) value else 1e100
    }
    fit <- optim(
      par = initial_params,
      fn = objective,
      method = "L-BFGS-B",
      lower = parameter_bounds[, "lower"],
      upper = parameter_bounds[, "upper"],
      control = list(
        maxit = maxit,
        factr = 1e7,
        pgtol = tolerance
      )
    )
    elapsed <- proc.time()[[3L]] - start
    estimates <- as.numeric(fit$par)
    data.frame(
      estimator = "NPSMLE",
      Beta_hat = estimates[1L], Sigma_hat = estimates[2L],
      b_2_hat = estimates[3L], g_2_hat = estimates[4L],
      objective = as.numeric(fit$value), convergence = as.integer(fit$convergence),
      objective_evaluations = unname(as.integer(fit$counts["function"])),
      gradient_evaluations = unname(as.integer(fit$counts["gradient"])),
      elapsed_seconds = elapsed,
      status = if (all(is.finite(estimates))) "ok" else "nonfinite",
      message = if (is.null(fit$message)) "" else as.character(fit$message),
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }, error = function(e) {
    empty_fit("NPSMLE", conditionMessage(e), proc.time()[[3L]] - start)
  })
}

run_replication <- function(replication, sample_size, truth_row, y) {
  configuration <- as.integer(truth_row$configuration)
  truth <- as.numeric(truth_row[1L, param_names])
  names(truth) <- param_names

  estimation_simulator <- brock_hommes(
    seed = 3L,
    r = risk_free_rate
  )

  fits <- list()
  if ("NPSMLE" %in% estimators) {
    fits[["NPSMLE"]] <- fit_npsmle(
      y, estimation_simulator, sample_size, configuration, replication
    )
  }
  out <- do.call(rbind, fits)
  rownames(out) <- NULL
  cbind(
    data.frame(
      replication = as.integer(replication),
      configuration = configuration,
      N = as.integer(sample_size), M = as.integer(sample_size),
      Beta_true = unname(truth["Beta"]), Sigma_true = unname(truth["Sigma"]),
      b_2_true = unname(truth["b_2"]), g_2_true = unname(truth["g_2"]),
      stringsAsFactors = FALSE, check.names = FALSE
    ),
    out,
    data.frame(
      npsmle_draws = as.integer(sample_size),
      empirical_sample = "fixed_within_DGP",
      likelihood_aggregation = "log_mean_kernel_density",
      risk_free_rate = risk_free_rate,
      time_horizon = time_horizon,
      maxit = maxit, tolerance = tolerance,
      stringsAsFactors = FALSE
    )
  )
}

write_csv_atomic <- function(x, path) {
  temporary <- tempfile(pattern = "mc_", tmpdir = dirname(path), fileext = ".csv")
  write.csv(x, temporary, row.names = FALSE)
  if (!file.rename(temporary, path)) stop("Could not write ", path)
  invisible(path)
}

completed_replications <- function(results) {
  if (is.null(results) || !nrow(results)) return(integer())
  table_present <- table(results$replication, results$estimator)
  missing_estimators <- setdiff(estimators, colnames(table_present))
  if (length(missing_estimators)) return(integer())
  row_ids <- rownames(table_present)[rowSums(table_present[, estimators, drop = FALSE] > 0L) == length(estimators)]
  as.integer(row_ids)
}

combine_without_duplicates <- function(old, new) {
  if (is.null(old) || !nrow(old)) return(new)
  both <- rbind(old, new)
  key <- paste(both$replication, both$estimator, sep = "::")
  both <- both[!duplicated(key, fromLast = TRUE), , drop = FALSE]
  both[order(both$replication, both$estimator), , drop = FALSE]
}

parallel_apply <- function(values, fun, sample_size, truth_row, y,
                           cluster = NULL) {
  if (cores <= 1L) return(lapply(values, fun, sample_size, truth_row, y))
  if (.Platform$OS.type != "windows") {
    return(parallel::mclapply(
      values, fun, sample_size, truth_row, y,
      mc.cores = cores, mc.preschedule = FALSE, mc.set.seed = FALSE
    ))
  }
  parallel::parLapply(cluster, values, fun, sample_size, truth_row, y)
}

summarise_results <- function(results) {
  if (!nrow(results)) return(data.frame())
  keys <- unique(results[c("configuration", "N", "M", "estimator")])
  keys <- keys[order(keys$configuration, keys$N, keys$estimator), , drop = FALSE]
  estimates <- paste0(param_names, "_hat")
  truths <- paste0(param_names, "_true")
  rows <- vector("list", nrow(keys) * length(param_names))
  index <- 0L
  for (i in seq_len(nrow(keys))) {
    keep <- results$configuration == keys$configuration[i] &
      results$N == keys$N[i] & results$M == keys$M[i] &
      results$estimator == keys$estimator[i]
    block <- results[keep, , drop = FALSE]
    complete_fit <- block$status == "ok" &
      rowSums(is.finite(as.matrix(block[, estimates, drop = FALSE]))) == length(estimates)
    convergence_rate <- mean(complete_fit & block$convergence == 0L, na.rm = TRUE)
    for (j in seq_along(param_names)) {
      index <- index + 1L
      valid <- is.finite(block[[estimates[j]]])
      estimate <- block[[estimates[j]]][valid]
      truth <- unique(block[[truths[j]]])
      truth <- truth[is.finite(truth)][1L]
      error <- estimate - truth
      rows[[index]] <- data.frame(
        configuration = keys$configuration[i], N = keys$N[i], M = keys$M[i],
        estimator = keys$estimator[i], parameter = param_names[j], truth = truth,
        mean_estimate = if (length(estimate)) mean(estimate) else NA_real_,
        bias = if (length(error)) mean(error) else NA_real_,
        sd = if (length(estimate) > 1L) stats::sd(estimate) else NA_real_,
        rmse = if (length(error)) sqrt(mean(error^2)) else NA_real_,
        n_requested = mc_replications,
        n_completed = length(unique(block$replication)),
        n_finite = sum(valid),
        convergence_rate = convergence_rate,
        failure_rate = mean(!complete_fit),
        median_objective_evaluations = stats::median(
          block$objective_evaluations, na.rm = TRUE
        ),
        median_elapsed_seconds = stats::median(block$elapsed_seconds, na.rm = TRUE),
        stringsAsFactors = FALSE, check.names = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

cat("\nBrock-Hommes Monte Carlo experiment\n")
cat("  configurations:", paste(combination_ids, collapse = ", "), "\n")
cat("  sample sizes (N = M):", paste(sample_sizes, collapse = ", "), "\n")
cat("  replications per cell:", mc_replications, "\n")
cat("  estimators:", paste(estimators, collapse = ", "), "\n")
cat("  design: fixed empirical sample, M = N, r =", risk_free_rate, "\n")
cat("  likelihood: log of mean Gaussian-kernel density\n")
cat("  workers:", cores, "\n")
cat("  results:", experiment_dir, "\n\n")

write.csv(
  truth_grid, file.path(experiment_dir, "parameter_configurations.csv"),
  row.names = FALSE
)

cluster <- NULL
if (cores > 1L && .Platform$OS.type == "windows") {
  cluster <- parallel::makePSOCKcluster(cores)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterExport(
    cluster, varlist = "required_files", envir = .GlobalEnv
  )
  parallel::clusterEvalQ(cluster, invisible(lapply(required_files, source)))
  parallel::clusterExport(
    cluster,
    varlist = c(
      "param_names", "base_seed", "risk_free_rate", "time_horizon",
      "estimators", "parameter_bounds",
      "maxit", "tolerance", "initial_params", "make_seed", "empty_fit",
      "kernel_log_likelihood", "fit_npsmle", "run_replication"
    ),
    envir = .GlobalEnv
  )
}

for (configuration in combination_ids) {
  truth_row <- truth_grid[truth_grid$configuration == configuration, , drop = FALSE]
  for (sample_size in sample_sizes) {
    y <- benchmark_data[[paste(configuration, sample_size, sep = "_")]]
    cell_file <- file.path(
      experiment_dir,
      sprintf("cell_config%02d_N%d_M%d.csv", configuration, sample_size, sample_size)
    )
    existing <- NULL
    if (file.exists(cell_file) && !overwrite) {
      existing <- read.csv(cell_file, stringsAsFactors = FALSE, check.names = FALSE)
    }
    done <- completed_replications(existing)
    pending <- setdiff(seq_len(mc_replications), done)

    cat(sprintf(
      "Configuration %02d, N = M = %d: %d completed, %d pending\n",
      configuration, sample_size, length(done), length(pending)
    ))

    if (length(pending)) {
      batches <- split(pending, ceiling(seq_along(pending) / batch_size))
      for (batch in batches) {
        batch_start <- proc.time()[[3L]]
        batch_results <- parallel_apply(
          batch, run_replication, sample_size, truth_row, y, cluster
        )
        batch_results <- do.call(rbind, batch_results)
        existing <- combine_without_duplicates(existing, batch_results)
        write_csv_atomic(existing, cell_file)
        cat(sprintf(
          "  saved through replication %d (%d/%d; %.1f seconds for batch)\n",
          max(batch), length(completed_replications(existing)), mc_replications,
          proc.time()[[3L]] - batch_start
        ))
      }
    }
  }
}

cell_files <- list.files(
  experiment_dir, pattern = "^cell_config[0-9]+_N[0-9]+_M[0-9]+[.]csv$",
  full.names = TRUE
)
all_results <- if (length(cell_files)) {
  do.call(rbind, lapply(cell_files, read.csv, stringsAsFactors = FALSE, check.names = FALSE))
} else {
  data.frame()
}
if (nrow(all_results)) {
  all_results <- all_results[
    order(all_results$configuration, all_results$N,
          all_results$replication, all_results$estimator),
    , drop = FALSE
  ]
  rownames(all_results) <- NULL
  write_csv_atomic(all_results, file.path(experiment_dir, "all_estimates.csv"))
  summary_results <- summarise_results(all_results)
  write_csv_atomic(summary_results, file.path(experiment_dir, "mc_bias_sd_rmse.csv"))
  saveRDS(
    list(
      estimates = all_results,
      summary = summary_results,
      configurations = truth_grid,
      settings = list(
        replications = mc_replications, sample_sizes = sample_sizes,
        combinations = combination_ids, estimators = estimators, cores = cores,
        npsmle_draws = "M = N",
        empirical_sample = "fixed within DGP and reused over replications",
        likelihood_aggregation = "log of mean Gaussian-kernel density",
        maxit = maxit, tolerance = tolerance,
        risk_free_rate = risk_free_rate, time_horizon = time_horizon,
        parameter_bounds = parameter_bounds,
        initial_params = initial_params
      )
    ),
    file.path(experiment_dir, "mc_results.rds")
  )
}

cat("\nExperiment outputs written to:", experiment_dir, "\n")
