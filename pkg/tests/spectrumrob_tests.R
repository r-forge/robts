library(robts)

set.seed(1654)

## Creating some test data:
depmodel <- list(ar=c(0.6, 0.3))
x <- arima.sim(model = depmodel, n = 100)


## Check the different methods:

spectrumrob(x, method = "pgram", plot = FALSE)

spectrumrob(x, method = "acf", plot = FALSE)


## Check cases with missing data:

x_missing <- c(x[1:2], NA, NA, x[5:100])

# Use longest complete strech of observations:

spectrumrob(x_missing, method = "pgram", plot = FALSE, na.action = na.contiguous)

spectrumrob(x_missing, method = "acf", plot = FALSE, na.action = na.contiguous)

# Set missings to extreme values:

spectrumrob(x_missing, method = "pgram", plot = FALSE, na.action = na.extremify)

spectrumrob(x_missing, method = "acf", plot = FALSE, na.action = na.extremify)
