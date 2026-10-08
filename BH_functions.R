# Adaptation in R of workshopFunctions.py from:
# https://github.com/Sylvain-Barde/abm-validation-cef-2023
# Original code Copyright (c) Sylvain Barde, distributed under the MIT License.

# Population standard deviation (NumPy's default ddof = 0).
population_sd <- function(x) {
  x <- as.numeric(x)
  sqrt(mean((x - mean(x))^2))
}

#' Construct a Brock and Hommes (1998) simulator.
#'
#' The parameter vector is c(beta, sigma, b_2, ..., b_H, g_2, ..., g_H).
#' Strategy 1 has b_1 = g_1 = 0, as in the Python workshop implementation.
#'
#' @param seed Integer base seed.
#' @param r Risk-free rate; the gross rate is 1 + r.
#' @return A list with step() and simulate() methods.
brock_hommes <- function(seed = 0L, r = 0.01) {
  seed <- as.integer(seed)
  gross_rate <- 1 + as.numeric(r)

  unpack <- function(params) {
    params <- as.numeric(params)
    if (length(params) < 4L || (length(params) - 2L) %% 2L != 0L) {
      stop("params must be c(beta, sigma, b_2..b_H, g_2..g_H)")
    }
    num_strategies <- 1L + (length(params) - 2L) %/% 2L
    list(
      beta = params[1L],
      sigma = params[2L],
      b = c(0, params[seq.int(3L, num_strategies + 1L)]),
      g = c(0, params[seq.int(num_strategies + 2L, length(params))])
    )
  }

  step <- function(y_lag, params, shock) {
    y_lag <- as.numeric(y_lag)
    if (length(y_lag) != 3L) stop("The Brock-Hommes step requires three lags")
    p <- unpack(params)

    beliefs <- p$g * y_lag[3L] + p$b
    lag_belief_term <- p$g * y_lag[1L] + p$b - gross_rate * y_lag[2L]
    utility <- p$beta * (y_lag[3L] - gross_rate * y_lag[2L]) * lag_belief_term
    exp_utility <- exp(pmax(-400, pmin(400, utility)))
    fractions <- exp_utility / sum(exp_utility)

    as.numeric((sum(beliefs * fractions) + p$sigma * shock) / gross_rate)
  }

  simulate <- function(T, params, rep = 0L, shocks = NULL) {
    T <- as.integer(T)
    if (T < 3L) stop("T must be at least 3")
    unpack(params) # Validate once before the loop.

    if (is.null(shocks)) {
      set.seed(seed + as.integer(rep))
      shocks <- qnorm(runif(T))
    }
    shocks <- as.numeric(shocks)
    if (length(shocks) != T) stop("shocks must have length T")

    y <- numeric(T)
    if (T > 3L) {
      for (t in 4L:T) y[t] <- step(y[(t - 3L):(t - 1L)], params, shocks[t])
    }
    y
  }

  structure(
    list(seed = seed, r = as.numeric(r), R = gross_rate, step = step,
         simulate = simulate),
    class = "brock_hommes"
  )
}

# Python-compatible alias used in the workshop material.
brockHommes <- brock_hommes

#' Non-parametric simulated maximum-likelihood helper.
npsmle <- function(emp_data, param_names = character()) {
  state <- new.env(parent = emptyenv())
  state$iter_count <- 0L
  state$timer <- proc.time()[[3L]]
  state$y <- as.numeric(emp_data)
  state$names <- as.character(param_names)

  log_like <- function(fun, num_lags, theta, epsilon) {
    epsilon <- as.numeric(epsilon)
    n_shocks <- length(epsilon)
    if (n_shocks < 2L) stop("epsilon must contain at least two draws")
    num_lags <- as.integer(num_lags)
    log_vec <- numeric(length(state$y))

    if (length(state$y) > num_lags) {
      for (i in seq.int(num_lags + 1L, length(state$y))) {
        y_sim <- as.numeric(fun(state$y[(i - num_lags):(i - 1L)], theta, epsilon))
        if (length(y_sim) != n_shocks) {
          stop("The one-step simulator must return one value per epsilon draw")
        }
        y_diff <- y_sim - state$y[i]
        bandwidth <- population_sd(y_diff) * (4 / (3 * n_shocks))^(1 / 5)
        if (!is.finite(bandwidth) || bandwidth <= 0) {
          stop("The kernel bandwidth is not positive at observation ", i)
        }
        # sklearn KernelDensity fitted on the singleton {0} is exactly this
        # Gaussian log density. Summing it preserves the upstream objective.
        # log_vec[i] <- sum(dnorm(y_diff, mean = 0, sd = bandwidth, log = TRUE))
        kernel_values <- dnorm(
          y_diff,
          mean = 0,
          sd = bandwidth
        )
        
        log_vec[i] <- log(
          pmax(mean(kernel_values), .Machine$double.xmin)
        )
      }
    }
    log_vec
  }

  callback <- function(xk) {
    now <- proc.time()[[3L]]
    if (state$iter_count == 0L) {
      names_out <- if (length(xk) == length(state$names)) {
        state$names
      } else {
        paste0("param-", seq_along(xk))
      }
      cat(sprintf("%11s  %11s  %s\n", "iteration", "time",
                  paste(sprintf("%11s", names_out), collapse = "  ")))
      cat(strrep("-", 13L * (length(xk) + 2L)), "\n", sep = "")
    }
    state$iter_count <- state$iter_count + 1L
    cat(sprintf("%11d  %11.4e  %s\n", state$iter_count,
                now - state$timer,
                paste(sprintf("%11.4e", xk), collapse = "  ")))
    state$timer <- now
    invisible(NULL)
  }

  structure(
    list(
      state = state,
      logLike = log_like,
      log_like = log_like,
      callback = callback
    ),
    class = "npsmle"
  )
}

#' Uncentred autocorrelation used by the Python workshop code.
autocorr <- function(x) {
  x <- as.numeric(x)
  n <- length(x)
  if (n < 1L) return(numeric())
  # convolve() uses an FFT for long vectors. Its centre and right-hand side
  # match numpy.correlate(x, x, mode = "full") at non-negative lags.
  full <- convolve(x, x, type = "open")
  raw <- full[seq.int(n, 2L * n - 1L)]
  raw / raw[1L]
}

formatTableText <- function(array, format_str = "%8.3f") {
  array <- as.matrix(array)
  matrix(sprintf(format_str, as.numeric(array)), nrow = nrow(array),
         ncol = ncol(array))
}

print_estimation_table <- function(values, row_names = NULL, col_names = NULL,
                                   digits = 3L, title = NULL) {
  values <- as.matrix(values)
  if (!is.null(row_names)) rownames(values) <- row_names
  if (!is.null(col_names)) colnames(values) <- col_names
  if (!is.null(title)) cat("\n", title, "\n", strrep("-", nchar(title)), "\n", sep = "")
  print(round(values, digits = digits))
  invisible(values)
}

# Calls fn through optim while printing a lightweight, deterministic progress
# report every `report_every` objective evaluations.
optim_bfgs <- function(fn, init, tolerance = 1e-3, maxit = 100L,
                       reporter = NULL, report_every = 20L) {
  evaluations <- 0L
  wrapped <- function(par) {
    evaluations <<- evaluations + 1L
    value <- fn(par)
    if (!is.null(reporter) && evaluations %% as.integer(report_every) == 0L) {
      reporter(par)
    }
    value
  }
  optim(
    par = as.numeric(init), fn = wrapped, method = "BFGS",
    control = list(reltol = tolerance, maxit = as.integer(maxit))
  )
}
