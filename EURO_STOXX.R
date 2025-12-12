rm(list=ls())

library(quantmod)

wd <- "/Users/Mario/PhD/Papers/EstSimMod/New/Codes"
setwd(wd)

# Define the ticker symbol and time range
ticker_symbol <- "^STOXX50E"
start_date <- "2011-12-30"
end_date <- "2013-01-02"

# Download historical data
getSymbols(ticker_symbol, src = "yahoo", from = start_date, to = end_date)
write.csv(as.numeric(STOXX50E$STOXX50E.Adjusted),file="EURO_STOXX_50.csv")
