read_sp500_deviations <- function(csv_path, start_date = "1994-02-23",
                                  end_date = "2013-12-31", half_window = 30L) {
  raw <- read.csv(csv_path, stringsAsFactors = FALSE)
  required <- c("Date", "SP500")
  if (!all(required %in% names(raw))) {
    stop("SP500 CSV must contain Date and SP500 columns")
  }
  raw$Date <- as.Date(raw$Date, format = "%d/%m/%Y")
  keep <- raw$Date >= as.Date(start_date) & raw$Date <= as.Date(end_date)
  selected <- raw[keep, required]
  selected <- selected[order(selected$Date), ]

  y <- as.numeric(selected$SP500)
  dates <- selected$Date
  w <- as.integer(half_window)
  num_raw_obs <- length(y)
  num_obs <- num_raw_obs - 2L * w - 1L
  if (num_obs < 1L) stop("The selected series is too short for the MA window")

  centres <- seq.int(w + 1L, num_raw_obs - w - 1L)
  fundamental <- vapply(
    centres,
    function(i) mean(y[(i - w):(i + w)]),
    numeric(1L)
  )
  observed <- y[centres]

  list(
    dates = dates[centres],
    observed = observed,
    fundamental = fundamental,
    deviation = observed - fundamental,
    num_raw_obs = num_raw_obs,
    num_obs = num_obs,
    half_window = w
  )
}
